# Palier 5 — PR 2 : congés à rayon variable, topologie, Python, OCCT

| | |
|---|---|
| PR | branche `feat/brep-variable-fillet-topology` |
| Document de référence | [brep_palier5_variable_fillets.md](brep_palier5_variable_fillets.md), § 3.5, § 4 à § 6 et § 9 (bilan) |
| Fichiers | `gbs-brep/blend.h`, `gbs-brep/fillet.h`, `python/gbsBindBrepOps.cpp`, `tests/tests_brep_variable_fillet.cpp`, `tests/tests_brep_variable_fillet_solid.cpp` (nouveau), `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/examples/gbs_examples.py`, `python/tests/test_brep_ops.py`, `docs/sources/design/brep_palier5_variable_fillets.md`, `News.md`, cette note |
| Statut | À valider (§ 9) |

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
suivent si elles sont seules à se toucher en un sommet et tangentes là
(|cos| > 1 − 10⁻⁶). Une arête fermée seule forme une boucle. `chain_radius`
(`fillet.h`) prend les points de la chaîne, avec l'abscisse
S = décalage de l'arête + s ou L − s selon son sens, et construit la loi de la
chaîne. Il en tire la loi de chaque arête (`edge_radius_in`), dont les nœuds et
angles tombent sur l'arête, ses sommets compris. `edge_radius` (PR 1) est la
chaîne d'une seule arête.

### 4.2 Fuseau ou arête vive (§ 3.5 du palier)

En chaque point candidat, la géométrie décide entre rien, fuseau et arête vive.
Les points candidats sont les angles de la loi, le sommet d'une arête fermée et
les nœuds de l'arête moins que C2.

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
- **Tangente prise du côté du morceau**. À un angle de l'arête, le plan d'une
  station est pris juste à l'intérieur du morceau, de même que la dérivée
  seconde de l'arête dans la direction de la ligne des centres. C'est le cas
  au sommet du pied d'aube : sa tangente y tourne de 1,2·10⁻³ rad, et les
  boules de part et d'autre diffèrent de 3,4·10⁻⁵.

### 4.4 Topologie

`blend()` (`blend.h`) découpe un congé en morceaux, une `Blend` par morceau
dans l'intervalle de l'arête. Les morceaux tout entiers au-delà d'un sommet ne
sont que le prolongement, et sont laissés.

- Entre deux morceaux d'une même arête, les contacts se rejoignent en un sommet
  partagé. Le fuseau est une face de sphère sur les deux arcs, chacun une arête
  avec ses pcourbes sur son morceau et sur la sphère (`sphere_pcurve`).
  L'arête vive est une seule arête, avec ses pcourbes sur les deux morceaux.
  Le dernier morceau d'une boucle et le premier partagent leur arc, ou leur
  fuseau, ou leur arête vive.
- Au sommet d'une chaîne, `joint` partage l'arc s'il est commun. Sinon, il
  applique le test du § 4.2 : un fuseau sur la boule du sommet, ou l'arête vive
  cherchée entre les deux congés prolongés.
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

## 6. Tests

`tests/tests_brep_variable_fillet_solid.cpp`, 7 cas :

| Test | Contenu |
|---|---|
| `box_edge_volume` | arête de boîte, lois lisse et linéaire, sommet (0,1 ; 0,3 ; 0,15 : arête vive, 8 faces) et creux (0,3 ; 0,1 ; 0,25 : fuseau, 9 faces) ; volume enlevé contre l'intégration indépendante de la réunion des boules (dans chaque plan x, le coin au-dessus des disques, Gauss-Legendre composite, maximum sur s échantillonné puis section d'or), à 2·10⁻⁷ près |
| `constant_law_is_constant_fillet` | même rayon partout : le volume du congé de rayon constant à 10⁻¹² |
| `cylinder_rim` | bord de cylindre, loi lisse fermée (4 faces) ; linéaire, trois morceaux, deux fuseaux (creux, dont un au sommet) et une arête vive (sommet) : 8 faces ; volume entre ceux des rayons constants extrêmes |
| `oblong_chain` | dessus d'un oblong, chaîne fermée de quatre arêtes, lois lisse et linéaire ; fuseaux ou arêtes vives aux raccords droite-cercle ; volume entre ceux des rayons constants extrêmes |
| `corner` | coin de boîte, rayon croissant vers le coin (10 faces, une sphère, volume entre ceux des rayons 0,15 et 0,25) ; décroissant : refusé |
| `refusals` | point sur une arête non congée, arête congée sans point, congé plus large qu'une face ; modèle inchangé |
| `blade_root` | pied d'aube, loi lisse fermée de 0,02 au sommet à 0,04 à mi-chemin ; boucle partie de la couture, fuseau mince au sommet G0 ; volume entre ceux des rayons 0,02 et 0,04 |

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

Suites complètes : 842 tests C++, 140 Python.

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
