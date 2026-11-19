# Palier 5 — PR 2 : congés à rayon variable, topologie, Python, OCCT

| | |
|---|---|
| PR | branche `feat/brep-variable-fillet-topology` |
| Document de référence | [brep_palier5_variable_fillets.md](brep_palier5_variable_fillets.md), § 3.5, § 4 à § 6 et § 9 (bilan) |
| Fichiers | `gbs-brep/blend.h`, `gbs-brep/fillet.h`, `python/gbsBindBrepOps.cpp`, `tests/tests_brep_variable_fillet.cpp`, `tests/tests_brep_variable_fillet_solid.cpp` (nouveau), `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/examples/gbs_examples.py`, `python/tests/test_brep_ops.py`, `docs/sources/design/brep_palier5_variable_fillets.md`, `News.md`, cette note |
| Statut | À valider (§ 9) ; revue de code : § 10 |

## 1. Objectif

Poser sur le solide les congés à rayon variable de la PR 1 :
`fillet(m, s, edges, RadiusLaw, options)`. Il faut pour cela :

- une loi le long des chaînes ;
- des morceaux, des fuseaux et des arêtes vives en faces et en arêtes ;
- les extrémités, les boucles et les coins ;
- le pied d'aube ;
- Python, la galerie, la comparaison à OCCT et le bilan du palier.

En route, deux décisions de conception. Le fuseau de la PR 1 n'existe qu'aux
**creux** de la loi. Aux **sommets**, la boule qui roule laisse une **arête
vive** : choix de l'utilisateur, le 19 nov. 2026, de la traiter exactement dans
cette PR. Le même choix vaut aux angles de la ligne des centres dus à l'arête
elle-même.

## 2. Sources

- Rossignac et Requicha (1984), Choi et Ju (1989) : la boule qui roule, son
  enveloppe.
- Vida, Martin et Várady (1994), § 5 : extrémités, chaînes, sommets.
- CATIA et SolidWorks pour l'usage. OCCT (`BRepFilletAPI_MakeFillet::Add(R1,
  R2, E)`) en boîte noire.

## 3. API

```cpp
auto fillet(Model<T> &, SolidId, std::span<const EdgeId> edges, const RadiusLaw<T> &, const FilletOptions<T> & = {}) -> BuildResult<SolidId>;
```

La loi court le long de chaque chaîne d'arêtes tangentes données : ses points
sur les arêtes de la chaîne, rangés par l'abscisse le long de la chaîne. Une
chaîne de même rayon partout donne le congé de rayon constant. `EdgeBlend` gagne
`creases` (`BlendCrease` : paramètre, courbe, pcourbes sur les deux morceaux,
écart).

Python : `brep.RadiusLaw` (avec `RadiusLaw.linear`), `brep.RadiusPoint`,
`brep.RadiusVariation`, `brep.fillet(m, s, edges, law)`,
`brep.fillet_geometry(m, s, edge, law)`, `BlendPiece`, `BlendLune`,
`BlendCrease`.

## 4. Démarche

### 4.1 Loi le long d'une chaîne

`law_chains` (`blend.h`) groupe les arêtes données en chaînes. Deux arêtes se
suivent si elles sont seules à se toucher en un sommet, ouvertes, et si l'une
y repart dans le sens où l'autre arrive (cos < −(1 − 10⁻⁶) entre les
directions qui quittent le sommet ; un rebroussement n'est pas une suite). Une
arête fermée seule forme une boucle. `chain_radius` (`fillet.h`) prend les
points de la chaîne, avec l'abscisse S = décalage de l'arête + s ou L − s
selon son sens, et construit la loi de la chaîne. Il en tire la loi de chaque
arête (`edge_radius_in`) :

- ses nœuds et angles strictement dans l'arête (son sommet aussi pour une
  arête fermée) ;
- au-delà des sommets d'une arête ouverte, la loi de l'arête continuée avec sa
  pente au sommet, prise sur l'arête : un angle de la loi au sommet d'une
  chaîne est tranché une seule fois, par `joint` (§ 4.4) ;
- son plus grand rayon, pris sur l'arête seule : aux sommets et aux nœuds de la
  loi dans l'arête (une loi linéaire est droite entre ses angles, une loi lisse
  monotone entre ses points) ;
- les sommets où la chaîne continue (`joined`) : le congé y est prolongé d'au
  moins 1,5 rayon, quelle que soit l'option `extension`, la place de l'arête
  vive avec le congé suivant.

`edge_radius` (PR 1) est la chaîne d'une seule arête.

### 4.2 Fuseau ou arête vive (§ 3.5 du palier)

En chaque point candidat, la géométrie décide entre rien, fuseau et arête vive.
Les points candidats sont les angles de la loi, le sommet d'une arête fermée et
les nœuds de l'arête moins que C2, tous dans l'arête. Les stations du test
partent de celles de la loi déjà cheminées (§ 4.3), sur la branche du congé.

1. Elle calcule les deux stations du point, chacune prise d'un côté : pente de
   la loi, tangente et dérivée seconde de l'arête juste à l'intérieur.
2. Elle compare leurs deux arcs. S'ils sont à `tol` l'un de l'autre, il n'y a
   pas de coupe.
3. Sinon, elle regarde si le milieu du fuseau entre les deux arcs, sur la
   boule du point, est dans une boule à 10⁻⁴ (b − a) avant ou après. Recouvert,
   c'est une arête vive ; découvert, un fuseau.

La PR 1 tranchait par la loi seule et posait toujours un fuseau. La sonde de la
réunion des boules a montré qu'à un sommet, ce fuseau est enfoncé jusqu'à
3,5·10⁻³ dans les boules voisines.

Le cheminement de l'arête vive (`crease_curve`) a demandé deux réglages. Les
points sont resserrés vers les contacts selon Tchebychev. Un point est accepté
à tol / 100 près, car près des contacts les deux surfaces deviennent tangentes
et Gauss-Newton n'atteint plus 10⁻¹¹. Avec 1024 points répartis
uniformément, l'écart restait à 5·10⁻⁵ sur l'oblong de la galerie. Avec ces
réglages, il tombe à 2·10⁻⁸.

### 4.3 Stations sur la bonne branche

- **Cheminement préalable** (`prime`). Les stations de chaque loi (celle de
  l'arête, celle de chaque morceau prolongé) sont d'abord cheminées d'un bout à
  l'autre en 512 pas. Le pas est réduit là où Gauss-Newton échoue, ou là où un
  contact bouge plus que 20 fois le chemin de l'arête. Au bord de fuite du pied
  d'aube, de rayon de courbure 0,012 pour un congé de 0,02 à 0,04, le plan
  normal à l'arête coupe aussi l'autre côté de l'aube. Une station partie de
  loin y sautait sur cette autre branche, puis ses voisines la suivaient. Les
  dérivées aux nœuds passaient alors de 3 à 4·10⁴.
- **Rayon des morceaux prolongés**. Un morceau prolongé de 1,5 rayon au-delà
  d'un sommet de la loi garde sa pente, et son rayon y dépasse le plus grand
  rayon de la loi : 0,41 pour un sommet de 0,3 sur l'arête de boîte. Les
  supports, leur portée et le contrôle de repli des offsets sont refaits pour
  ce rayon (le congé recalculé avec 5 % de plus). Sans cela, un dièdre de 30°
  à loi linéaire refusait le congé dès un sommet de 0,2.
- **Stations réutilisées**. Chaque morceau reprend les stations de la loi
  cheminées strictement à l'intérieur, et ne chemine que ses prolongements, à
  partir de ses bouts. Le recalage de la boucle sur la couture décale les
  stations d'une période au lieu de les recheminer.
- **Tangente prise du côté du morceau**. À un angle de l'arête, le plan d'une
  station est pris juste à l'intérieur du morceau, de même que la dérivée
  seconde de l'arête dans la direction de la ligne des centres. C'est le cas
  au sommet du pied d'aube : sa tangente y tourne de 1,2·10⁻³ rad, et les
  boules de part et d'autre diffèrent de 3,4·10⁻⁵.

### 4.4 Topologie

`blend()` (`blend.h`) découpe un congé en morceaux, une `Blend` par morceau.
Les coupes tombent dans l'arête : le premier et le dernier morceau d'une arête
ouverte courent au-delà de ses sommets.

- Entre deux morceaux d'une même arête, les contacts se rejoignent en un sommet
  partagé. Le fuseau est une face de sphère sur les deux arcs, chacun une arête
  avec ses pcourbes sur son morceau et sur la sphère (`sphere_pcurve`).
  L'arête vive est une seule arête, avec ses pcourbes sur les deux morceaux.
  Le dernier morceau d'une boucle et le premier partagent leur arc, ou leur
  fuseau, ou leur arête vive.
- Au sommet d'une chaîne, `joint` partage l'arc s'il est commun. Sinon, pour
  des congés à loi seulement (`EdgeBlend::variable` ; des chanfreins ou des
  congés de rayon constant qui ne se raccordent pas restent refusés), il
  applique le test du § 4.2 : un fuseau sur la boule du sommet, ou l'arête vive
  cherchée entre les deux congés, prolongés au moins de 1,5 rayon au-delà du
  sommet (§ 4.1).
- Une face de fuseau est orientée contre la face du morceau le long de leur
  arc commun, une fois les faces des morceaux posées. Une règle tirée de la
  normale de la sphère donnait une coque non orientable sur le bord de
  cylindre.
- Coin : les bouts des congés sont des cercles de la boule du coin. Un rayon
  qui décroît vers le coin décale ces cercles vers le coin, et les boules
  voisines recouvrent la sphère : refus.
- Boucle autour d'une B-spline fermée non périodique : elle part du croisement
  de la couture, affiné par bisection sur les stations. L'arête fermée est lue
  par `LoopCurve`.
- La boucle des extrémités libres saute les bouts en fuseau et en arête vive.

## 5. Refus

| Situation | Code | Message |
|---|---|---|
| point de la loi sur une arête non congée | `InvalidBounds` | `a point of the law on an edge not filleted` |
| chaîne sans point | `InvalidBounds` | `no radius on a chain of edges` |
| rayon qui décroît vers un coin | `NonGeneric` | `fillets at a corner whose radius drops toward it (the balls next to the corner cover its sphere)` |
| congés de deux arêtes d'une chaîne sans arête vive trouvée | `NonGeneric` | `the fillets of two edges of a chain do not meet along a crease at their vertex` |
| arête vive introuvable entre deux morceaux | `NonGeneric` | `the fillets on either side of a peak of the radius law do not meet along a crease` |
| stations impossibles à suivre le long de l'arête | `NonGeneric` | `the balls of the fillet cannot be followed along the edge` |
| arête dégénérée donnée | `InvalidId` | `a degenerate edge (no curve)` |
| arcs différents au sommet de deux congés sans loi (chanfreins, rayon constant) | `NonGeneric` | `blended edges meeting at a vertex but not tangentially, or blends that do not join there (not handled)` |

## 6. Tests

`tests/tests_brep_variable_fillet_solid.cpp`, 10 cas :

| Test | Contenu |
|---|---|
| `box_edge_volume` | arête de boîte, lois lisse et linéaire, sommet (0,1 ; 0,3 ; 0,15 : arête vive, 8 faces) et creux (0,3 ; 0,1 ; 0,25 : fuseau, 9 faces) ; volume enlevé contre l'intégration indépendante de la réunion des boules (dans chaque plan x, le coin au-dessus des disques, Gauss-Legendre composite, maximum sur s échantillonné puis section d'or), à 2·10⁻⁷ près |
| `constant_law_is_constant_fillet` | même rayon partout : le volume du congé de rayon constant à 10⁻¹² |
| `cylinder_rim` | bord de cylindre, loi lisse fermée (4 faces) ; linéaire, trois morceaux, deux fuseaux (creux, dont un au sommet) et une arête vive (sommet) : 8 faces ; volume entre ceux des rayons constants extrêmes |
| `oblong_chain` | dessus d'un oblong, chaîne fermée de quatre arêtes, lois lisse et linéaire ; lisse : deux fuseaux, aux raccords où la pente de la loi n'est pas nulle ; linéaire : quatre fuseaux et une arête vive ; le même solide avec `extension` = 0 ; volume entre ceux des rayons constants extrêmes |
| `corner` | coin de boîte, rayon croissant vers le coin (10 faces, une sphère, volume entre ceux des rayons 0,15 et 0,25) ; décroissant : refusé |
| `refusals` | point sur une arête non congée, arête congée sans point, congé plus large qu'une face ; modèle inchangé |
| `blade_root` | pied d'aube, loi lisse fermée de 0,02 au sommet à 0,04 à mi-chemin ; boucle partie de la couture, deux morceaux et un fuseau mince au sommet G0 ; volume entre ceux des rayons 0,02 et 0,04 |
| `peak_on_an_acute_edge` | dièdre de 30°, loi linéaire 0,02 ; 0,25 ; 0,02 : morceaux prolongés jusqu'à 0,34, deux morceaux et une arête vive, solide valide |
| `largest_radius_per_edge` | oblong à demi-cercles de rayon 0,3, loi à 0,15 sauf 0,31 au milieu d'un côté droit : accepté (le plus grand rayon de chaque arête est le sien) |
| `degenerate_edge` | arête au pôle d'une boule : `InvalidId` |

`tests/tests_brep_variable_fillet.cpp` (géométrie, PR 1) : arête vive au lieu
du fuseau sur l'arête de boîte linéaire (un sommet). Sur le bord de cylindre
linéaire, deux fuseaux et une arête vive. Un cas de plus, `hollow_and_peak`. Le
contrôle commun vérifie les arêtes vives : de contact à contact, sur les deux
morceaux, écart sous `approx_tol`.

`gbs-occt/tests/tests_brep2_occt.cpp`, `variable_fillet_solids`, à côté de
`Add(R1, R2, E)` sur deux arêtes de boîte. Le résultat de gbs relu par OCCT est
valide, a autant de faces et le même volume. Les volumes d'OCCT sont seulement
rapportés : sur l'arête du dessus, de 0,1 à 0,3, OCCT enlève 0,01895, des arcs
de rayon r(x) dans les plans normaux 0,01860 (calcul à la main), la boule qui
roule de gbs 0,01854. Sa loi linéaire n'est ni l'une ni l'autre, et on ne lit
pas ses sources.

Python : `test_variable_fillet` (`test_brep_ops.py`), et trois cas de galerie
(boîte à deux lois linéaires, sommet et creux ; oblong à loi lisse ; pied d'aube
à loi lisse), valides et fermés une fois maillés.

Une sonde, hors des tests, a vérifié que les points des fuseaux (creux) et des
arêtes vives (sommets) sont sur le bord de la réunion des boules à 10⁻⁹ près.

Suites complètes : 845 tests C++, 140 Python.

## 7. Limites

- Coins : seulement trois congés de même rayon au sommet, le rayon non
  décroissant vers le coin.
- Une arête G0 (le pied d'aube à sa fermeture) donne un fuseau mince de la
  largeur du saut des boules (3,4·10⁻⁵), plutôt qu'un raccord à l'échelle
  de la tolérance.
- Le pied d'aube à rayon variable demande environ 3 s, pour le cheminement
  préalable et le raffinement au bord de fuite.
- Les volumes d'OCCT ne sont pas comparés.

## 8. Suite

Fin du palier 5 (bilan : § 9 du document du palier). La suite de la feuille de
route : raccords de sommet quelconques, chanfreins à distance variable,
sections coniques.

## 9. Points à valider

1. **Arête vive aux sommets** de la loi (décidée le 19 nov. 2026), **fuseau aux
   creux**, choisis par un test de recouvrement géométrique et non par la loi.
2. Le même test **aux angles de la ligne des centres dus à l'arête** (raccord
   droite-cercle d'une chaîne, sommet G0 d'une arête fermée), pour une loi
   lisse aussi.
3. **Loi le long des chaînes** par l'abscisse curviligne de la chaîne
   (`law_chains`, `chain_radius`).
4. **Orientation d'un fuseau** tirée de la face de son morceau.
5. **Cheminement préalable** des stations (512 pas, pas réduit sur échec ou
   saut).
6. **Arête vive** : points Tchebychev, acceptés à tol / 100.
7. **Coin** : bouts en petits cercles quand le rayon croît vers le coin, refus
   quand il décroît.
8. **Volumes d'OCCT rapportés, non comparés.**
9. **Volume de référence** par intégration indépendante de la réunion des
   boules (arête de boîte).

## 10. Revue de code

Revue du 19 nov. 2026, 15 points, tous corrigés :

1. **Rayon des morceaux prolongés** au-delà d'un sommet de la loi, plus grand
   que le plus grand rayon de la loi : supports, portée et repli refaits pour
   lui (§ 4.3). Sonde : un dièdre de 30° à loi linéaire était refusé dès un
   sommet de 0,2 ; il passe jusqu'à 0,35. Test `peak_on_an_acute_edge`.
2. **Plus grand rayon d'une arête de chaîne** : il était pris sur toute la
   chaîne. Sonde : des demi-cercles de rayon 0,3, à 0,15 le long d'eux, étaient
   refusés dès que la chaîne atteignait 0,31. Test `largest_radius_per_edge`.
3. **`joint`** ne pose un fuseau ou une arête vive qu'entre congés à loi ; les
   chanfreins et les congés de rayon constant qui ne se raccordent pas restent
   refusés (`EdgeBlend::variable`).
4. **Sommet de chaîne tranché une fois**, par `joint` : plus de coupe au sommet
   d'une arête ouverte, la loi y continue avec sa pente sur l'arête, et le
   congé est prolongé d'au moins 1,5 rayon où la chaîne continue, quelle que
   soit `extension` (§ 4.1). Test : l'oblong avec `extension` = 0.
5. **Arête vive calculée puis jetée** au sommet d'une chaîne : disparue avec le
   point 4.
6. **Stations du test fuseau ou arête vive** prises sur celles déjà cheminées.
7. **Arête dégénérée** : refusée (`InvalidId`) au lieu de « map::at ». Test
   `degenerate_edge`.
8. **Doc de `t1`, `t2`** : pour une arête fermée, les bouts de la boucle (C++ et
   Python).
9. **`piece_view`** vide aussi les arêtes vives.
10. **`law_chains`** : tangence orientée, un rebroussement ne prolonge pas la
    chaîne, une arête fermée non plus.
11. **Cheminement par morceau** : stations de la loi reprises, prolongements
    seuls cheminés ; recalage sur la couture sans nouveau cheminement (§ 4.3).
12. **Nœuds faibles** : une fonction `weak_knots` sur les nœuds plats, pour les
    surfaces et pour l'arête (`edge_weak_knots`), au lieu de trois copies.
13. **`continued_edge`** appelée une fois par passe ; les bornes du congé en
    sont tirées.
14. **`lune_face`** : la variable morte retirée, une seule sphère en argument ;
    `Blend::lune_at` retiré (plus de fuseau calculé au sommet d'une arête
    ouverte).
15. **Tests** : nombres de fuseaux épinglés (`oblong_chain` : 2 et 4 ;
    `blade_root` : 1).

Hors revue : le nouveau test OCCT employait `Standard_True`, déprécié en OCCT 8
(un avertissement du job occt8 depuis la PR 236) ; remplacé par `true`.
