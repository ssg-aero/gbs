# Palier 8 — PR 3 : congé tritangent et bilan du palier

| | |
|---|---|
| PR | branche `feat/brep-full-round` (sur master e8c3114, après #249) |
| Document de référence | [brep_palier8_face_blends.md](brep_palier8_face_blends.md), § 3 (questions 8, 9 et 10 validées le 23 nov. 2026) |
| Fichiers | `gbs-brep/full_round.h` (nouveau), `gbs-brep/brep`, `gbs-brep/fillet.h` (arc par un point), `gbs-brep/blend.h` (orientation de la face du congé, message de fin libre), `tests/tests_brep_full_round.cpp` (nouveau), `gbs-occt/tests/tests_brep2_occt.cpp`, `python/gbsBindBrepOps.cpp`, `python/tests/test_brep_ops.py`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier8_face_blends.md` (bilan, § 9), `News.md`, cette note |
| Statut | À valider (§ 9) |

## 1. Objectif

Le congé tritangent, comme le « Tritangent Fillet » de CATIA et le « Full round
fillet » de SolidWorks : la face du milieu M, entre deux faces F₁ et F₂, est
remplacée par le congé tangent aux trois.

```cpp
fillet_full_round(m, s, faces1, middle, faces2, options);     // Python : brep.fillet_full_round
full_round_geometry(m, s, face1, middle, face2, options);      // la géométrie seule (EdgeBlend)
```

- **Aucun rayon n'est donné** : c'est celui de la boule qui touche les trois
  faces, et il suit la géométrie le long de M (la moitié de la largeur d'une
  nervure à flancs parallèles).
- Les cas visés : le dessus d'une nervure arrondi, avec ou sans dépouille ; le
  fond d'une rainure (concave, le congé ajoute de la matière) ; le **bord de
  fuite épais d'une aube** entre ses faces B-spline.
- Le bilan du palier 8 (§ 9 de [brep_palier8_face_blends.md](brep_palier8_face_blends.md)).

## 2. Sources

- Vida, Martin et Várady, « A survey of blending methods that use parametric
  surfaces », *Computer-Aided Design* 26(5), 1994, § 5 : congés à trois
  surfaces, boule tangente aux trois.
- L'enveloppe d'une famille de sphères (surface canal) : une boule de rayon
  variable touche son enveloppe le long d'un **cercle caractéristique**, qui
  n'est pas un grand cercle (palier 5, § 3).
- Le comportement de CATIA et de SolidWorks : la face du milieu disparaît, les
  faces de bout coupent le congé, aucun rayon à saisir.
- Le code des paliers 4 et 5 (`skinned`, `free_end`, `blend_on_copy`) et la
  station à rayon inconnu de la corde constante (palier 7, § 4.3).

## 3. API

```cpp
template <std::floating_point T>
auto fillet_full_round(Model<T> &m, SolidId s, std::span<const FaceId> faces1, std::span<const FaceId> middle,
                       std::span<const FaceId> faces2, const FilletOptions<T> &o = {}) -> BuildResult<SolidId>;
// (et la surcharge sur des std::vector)

template <std::floating_point T>
auto full_round_geometry(const Model<T> &m, SolidId s, FaceId face1, FaceId middle, FaceId face2,
                         const FilletOptions<T> &o = {}) -> BuildResult<EdgeBlend<T>>;
```

- La signature de `fillet_full_round` est celle du palier (§ 0) : des ensembles
  de faces de chaque côté et au milieu. **Écart** : cette PR n'en prend qu'une
  par ensemble, et refuse les autres (§ 9, point 4).
- `full_round_geometry` est ajoutée sur le modèle de `fillet_geometry` : la
  surface, la ligne des centres (`spine` ; le rayon en t est la distance de la
  ligne des centres à un contact), les contacts sur F₁ et F₂, tout sur le
  paramètre de l'arête M ∩ F₁. Elle sert aux tests et à la galerie.
- Les options sont celles des congés : `tol`, `approx_tol`, `angular_tol`,
  `extension`. La section doit être l'arc de cercle (le palier exclut les
  sections conique et G2 du tritangent) ; les retraits sont refusés.
- Python : `brep.fillet_full_round(m, s, faces1, middle, faces2, options)` et
  `brep.full_round_geometry(m, s, face1, middle, face2, options)`.

## 4. Démarche

### 4.1 Station (§ 3.1 du palier)

Le **guide** est l'arête g = M ∩ F₁. Une station se pose en g(t), dans le plan
normal à g en ce point. Les inconnues sont :

- les paramètres (u, v) des trois contacts ;
- le centre c ;
- le rayon r.

Soit 10 inconnues. Les équations sont :

- c = Sₖ(uₖ, vₖ) ∓ r Nₖ(uₖ, vₖ) pour F₁, M, F₂ (9 équations, la surface décalée
  de r vers la boule) ;
- (c − g(t)) · g′(t) = 0, le centre dans le plan de la station.

Gauss-Newton (`thick_detail::gauss_newton`), depuis la station voisine.

**Départ.** La première station part du **cercle inscrit** dans les tangentes
des trois faces, dans le plan de la station. On prend :

- le point de l'autre arête M ∩ F₂ dans ce plan (Newton sur sa projection) ;
- la largeur w de M entre les deux arêtes ;
- les directions des côtés, qui s'éloignent de M.

Le triangle formé par M et les deux côtés donne un rayon et un centre fermés.
Les contacts sont déplacés d'autant dans le plan tangent de chaque face. Le
même départ sert si une station voisine saute sur une autre racine (des
contacts qui bougent plus que le guide ne peut les porter : le même critère
que `fillet_law`).

**Côté de la boule.** La boule est dans la matière si M est **convexe** le
long de ses deux arêtes avec les côtés : le congé ôte de la matière (une
nervure). Elle est dehors si M est **concave** le long des deux : le fond d'une
rainure, le congé ajoute de la matière. Le côté est lu par `fillet_detail::setup`
sur chacune des deux arêtes. *Écart* : le palier (§ 3.1) ne prévoyait que la
boule dans la matière ; le cas concave ne coûte rien de plus (§ 9, point 3).

**Marche.** Les stations partent du milieu du guide vers ses deux sommets
(32 pas de chaque côté), puis le continuent au-delà de `extension` r, tant
qu'une boule existe. Au-delà d'un bout où M finit en pointe, il n'y en a plus :
le congé s'arrête à la dernière station résolue.

### 4.2 Section et surface (§ 3.2 du palier)

**Section.** La section est l'arc du **cercle caractéristique**, le cercle qui
passe par les trois contacts. Il va du contact sur F₁ au contact sur F₂ en
passant par le contact P_M.

- À rayon constant, c'est un grand cercle de la boule : l'arc de la boule par
  P_M du palier.
- À rayon variable, une boule touche son enveloppe le long d'un cercle qui
  n'est pas un grand cercle. Ses trois contacts sont sur ce cercle, dont le
  centre est hors de c. *Précision au palier* : le palier disait « l'arc de la
  boule », ce qui n'est exact qu'à rayon constant (§ 9, point 2).

L'arc dépasse 180° quand les flancs se rapprochent en s'éloignant de M. Il fait
exactement 180° entre des flancs parallèles, alors que l'arc le plus court
entre deux points opposés n'est pas défini. Une fonction nouvelle de
`fillet.h`, `arc_through(c, p1, p2, via)`, en tient compte. C'est l'arc
rationnel quadratique en deux segments, de p1 à p2 par le côté de `via`, joint
en son milieu, refusé au-delà de 300°. `Station` gagne un champ facultatif
`via` : `skinned` prend alors `arc_through` au lieu de `arc`. Les congés
d'arête ne posent pas `via` et ne changent pas.

**Surface.** Deux cas :

- **cylindre exact** si les trois faces sont planes, le guide droit et le rayon
  constant. La nervure droite, avec ou sans dépouille, est dans ce cas. Le
  rayon vaut alors (w/2) cos α / (1 − sin α), et l'arc couvre π − 2α. Le guide
  est jugé droit sur ses points, quelle que soit sa courbe (un `LineCurve` ou
  une B-spline de degré 1 : une arête d'extrusion est souvent la seconde) ;
- sinon, la **peau de Hermite** des sections (`skinned`), sous `approx_tol` et
  `angular_tol`, avec une station là où un contact croise un nœud d'une face sous C³
  et aux nœuds du guide sous C³.

**Emboîtement** (§ 3.2 du palier). Entre deux stations voisines, si
|Δr| ≥ |Δc|, une boule contient l'autre : il n'y a pas d'enveloppe. C'est
refusé ; le contrôle est fait sur les 64 intervalles du congé.

### 4.3 Topologie (§ 3.3 du palier) : M réduite à son guide

Le palier décrit le résultat : M ôtée, F₁ et F₂ bornées par les contacts, une
fin libre sur chaque face de bout. C'est exactement ce que fait le congé d'une
arête vive entre F₁ et F₂. La PR obtient donc cette topologie par la machinerie
des congés d'arête, sur une copie :

1. **Copie** du solide.
2. **M réduite à son guide** (`collapse`) :
   - M est ôtée de sa coque ;
   - dans F₂, l'arête M ∩ F₂ est remplacée par le guide (la pcourbe : le
     segment entre les bouts de l'ancienne co-arête) ;
   - les arêtes de M avec ses faces de bout sont réduites au sommet du guide de
     leur côté, et ôtées des contours des faces de bout.
3. **Congé du guide** sur cette copie (`blend_detail::blend`), avec la géométrie
   du § 4.1–4.2 calculée sur le solide donné. Les contacts coupent F₁ et F₂ ;
   chaque face de bout coupe le congé (`free_end`) ; le résultat est contrôlé
   (`check`).
4. La copie réduite est effacée ; seul le résultat reste.

La copie réduite a la **topologie** d'un solide à arête vive, mais pas sa
géométrie : F₂ est à une largeur de M du guide. Aucune opération n'en lit la
géométrie entre F₂ et le guide : la découpe de F₂ suit le contact, qui est loin
du guide, et la fin libre coupe le congé par la face de bout prolongée. Le
contrôle final porte sur le résultat, qui est géométriquement juste.

*Écart* : le palier ne disait pas comment construire cette topologie. Écrire
une topologie à part aurait dupliqué la découpe des faces, les fins libres et
leurs contrôles (§ 9, point 1).

### 4.4 Orientation de la face du congé (`blend.h`)

`blend_face` orientait la face du congé par la **somme des normales des deux
faces**, comparée à la normale du congé au milieu de sa section. Entre deux
flancs parallèles, cette somme est nulle ; au-delà de 180°, elle pointe à
l'envers. Le premier essai donnait une coque non orientable.

L'orientation est maintenant lue **à chaque contact** : la normale du congé y
est celle de la face qu'il touche (le congé lui est tangent). Pour tous les
congés d'arête, les deux règles donnent le même résultat (leur arc fait moins
de 180°) : la suite des tests le confirme (§ 6).

### 4.5 Contacts dans leurs faces

- **Sur F₁ et F₂** : le contrôle des congés d'arête (`check_on_face`), le
  contact sur F₂ lu sur l'arête M ∩ F₂.
- **Sur M** : en 17 points du guide, le contact P_M doit être dans M. Il ne peut
  en sortir que par une arête de bout, comme un congé qui dépasse un coin. Au
  milieu, il doit être dedans.

### 4.6 Refus (§ 3.4 du palier)

Avec des messages qui disent quoi faire ; le modèle reste inchangé.

`InvalidId` :
- un côté ou le milieu vide, une face morte, qui n'est pas du solide, ou donnée
  deux fois ;
- une face du milieu sans arête avec un côté (le palier : « M non voisine de
  F₁ ou de F₂ »).

`InvalidBounds` :
- une section conique ou G2 ; des retraits ; des tolérances non positives.

`NonGeneric` :
- **aucune boule** tangente aux trois faces (Gauss-Newton sans solution, ou
  r ≤ tol) : « fillet its two edges with the sides instead » ;
- un **contact hors de sa face** : une nervure plus basse que le rayon (la
  boule touche les flancs sous leur pied), ou une boule qui sort de M ;
- **boules emboîtées** (§ 4.2) ;
- une **face de bout qui ne coupe pas** le congé (trop oblique ou trop
  étroite) : « the face at that end of the middle face does not cut the full
  round » ;
- M **convexe** le long d'un côté et **concave** le long de l'autre ;
- non traités : un côté ou un milieu de plusieurs faces, M qui fait le tour
  (le bord d'un disque), M trouée, M qui finit en pointe, plusieurs faces à un
  bout de M (une nervure en proue), M qui touche un côté par plusieurs arêtes.

## 5. Ce qui ne change pas

- Les congés d'arête et de faces : leur géométrie est inchangée (`via` n'est
  posé que par le tritangent). L'orientation de la face du congé (§ 4.4) donne
  le même résultat pour eux.
- Le message d'une fin libre qui ne rencontre pas la face de bout ne change que
  pour `fillet_full_round`.

## 6. Tests

`tests_brep_full_round` (9 cas, 0,1 s). Les volumes sont comparés aux formes
fermées.

| Cas | Attendu | Mesuré |
|---|---|---|
| nervure droite, w = 0,4 et 0,8 | cylindre de rayon w/2, arc de π ; aire ôtée w²/2 − π w²/8 par unité de longueur ; solide donné intact, rien de laissé dans le modèle | rayon à 10⁻¹², volume à 2,7·10⁻⁹ ; 2 ms |
| nervure à dépouille α = 5°, 15°, −10° | cylindre de rayon (w/2) cos α / (1 − sin α), arc π − 2α ; aire ôtée 2 r² (tan β − β), β = π/4 − α/2 | rayon à 10⁻¹², volume à 1,2·10⁻⁹ ; 2 ms |
| nervure droite, `tol` = 10⁻⁹ | même volume | 9,5·10⁻¹³ |
| fond de rainure (concave) | cylindre de rayon w/2, deux coins ajoutés | 1,0·10⁻⁹ ; 2 ms |
| nervure à bout oblique (30°, 60°) | la face de bout coupe le demi-cylindre selon une ellipse ; volume des deux coins prismatiques | 8,5·10⁻¹⁰ ; 4 à 5 ms |
| nervure qui se rétrécit (0,6 → 0,2, flancs verticaux) | peau approchée ; les centres sur la droite issue du sommet A des trois plans, le rayon proportionnel à la distance à A ; le congé sur le cône de sommet A | écart au cône 1,9·10⁻¹² ; 28 ms |
| bord de fuite d'aube (faces B-spline cubiques, épaisseur 0,08 → 0,044, vrillée, en flèche) | peau sous `approx_tol` ; contacts tangents aux côtés (sinus de l'angle des normales ≤ `angular_tol`) ; rayon lisse, entre 0,02 et 0,06 | rayon 0,044 → 0,024, différences secondes 1,7·10⁻⁵ ; déviation 1,4·10⁻⁷ ; angle 2,8·10⁻⁷ ; 20 ms |
| refus des faces données | `InvalidId` (côté vide, face deux fois, sans arête avec un côté, face d'un autre solide), `InvalidBounds` (section conique), `NonGeneric` (plusieurs faces ; convexe et concave ; plusieurs faces à un bout) ; modèle inchangé | messages vérifiés |
| refus de forme | nervure plus basse que le rayon (contact hors de sa face) ; dessus en pointe ; nervure en proue ; bord d'un disque ; face de bout à 88° | messages vérifiés, modèle inchangé |

**Tolérance des volumes.** Le palier demandait 10⁻⁹ ; les tests prennent 10⁻⁸
à la tolérance par défaut (10⁻⁷). L'écart vient des courbes où les faces de bout
coupent le congé : des B-splines cubiques à `tol`/10 du demi-cercle exact, comme
au bout de tout congé, mais sur un demi-cercle au lieu d'un quart. À
`tol` = 10⁻⁹, l'écart tombe à 10⁻¹² (test `tolerance`).

**OCCT** (`tests_brep2_occt.full_round_ribs`). OCCT n'a pas de congé
tritangent. Le cas équivalent : la nervure abaissée jusqu'aux contacts,
fusionnée (`BRepAlgoAPI_Fuse`) avec le cylindre tangent, ses faces d'un même
domaine fondues. Nervure droite et nervure à dépouille de 10° : même solide,
mêmes faces (10), volume à 1,0·10⁻⁹ et 8,0·10⁻¹⁰. *Écart* : le palier (§ 6)
ne prévoyait que les formes fermées pour le tritangent.

**Python** : `test_fillet_full_round` (nervure droite : volume, cylindre,
rayon de `full_round_geometry` ; dépouille de 10° : rayon ; deux refus ;
docstring). La galerie gagne `full_round_rib` (dépouille de 8°, volume en forme
fermée) et `full_round_blade` (bord de fuite d'aube), vérifiés par
`test_examples.py` (solide valide, maillage fermé, volume).

## 7. Illustration

Note d'illustration (montrée à la validation, hors du dépôt), rendus réels de
la branche : la nervure droite et à dépouille avant / après, le fond de
rainure, le bout oblique, la nervure qui se rétrécit, le bord de fuite d'aube
avant / après et le rayon le long de l'aube.

## 8. Estimation

| Partie | Heures |
|---|---|
| station à trois faces, section par P_M, cylindre exact, peau | 6 |
| topologie (M réduite au guide), orientation, refus | 5 |
| tests (nervures, rainure, bouts, aube), OCCT | 4 |
| Python, galerie, bilan du palier, illustration | 3 |
| **Total** | **18** |

## 9. Points à valider

1. **Topologie par réduction de M à son guide** (§ 4.3), sur une copie, puis
   la machinerie des congés d'arête (découpe des faces, fins libres, contrôle).
   *Recommandé* : le résultat est celui que décrit le palier, sans dupliquer la
   topologie des paliers 4 à 6.
2. **Section par le cercle caractéristique** des trois contacts (§ 4.2) plutôt
   que « l'arc de la boule » : identique à rayon constant, seul juste à rayon
   variable (palier 5). *Recommandé.*
3. **Côté concave** (fond de rainure) accepté en plus de la boule dans la
   matière ; M convexe d'un côté et concave de l'autre refusée. *Recommandé* :
   c'est le même code, et SolidWorks l'accepte.
4. **Un seul côté et un seul milieu** par appel (pas d'ensembles de faces), et
   pas de M qui fait le tour (bord d'un disque, bord d'attaque d'une aube
   fermée sur elle-même), ni de M qui finit en pointe ou en proue : refusés avec
   « not handled ». C'est un écart au palier (§ 3.1, § 3.3). *Recommandé* pour
   cette PR : les cas visés (nervures, bord de fuite épais) y sont tous ; les
   ensembles demandent une chaîne de guides et le bord fermé une boucle sans
   face de bout, à reprendre dans une suite (bilan, § 9 du palier).
5. **Orientation de la face du congé par ses contacts** (§ 4.4), dans
   `blend.h`, pour tous les congés : même résultat sur les congés d'arête,
   seule juste au-delà de 180°. *Recommandé.*
6. **`arc_through`** et `Station::via` dans `fillet.h`, l'arc en deux segments
   par un point, jusqu'à 300°. *Recommandé.*
7. **Plan des stations normal au guide M ∩ F₁** (question 8 validée) : le
   paramètre du congé est celui de cette arête. Le résultat ne dépend pas du
   côté choisi comme F₁ pour une nervure symétrique ; il en dépend un peu
   (l'échantillonnage, pas la surface exacte) pour une aube. *Recommandé* tel
   que validé.
8. **Tolérance des volumes à 10⁻⁸** à la tolérance par défaut (§ 6), 10⁻¹² à
   `tol` = 10⁻⁹. *Recommandé* : l'écart est celui des fins libres de tout congé.
9. **`full_round_geometry`** ajoutée à l'API, comme `fillet_geometry`, et la
   **comparaison à OCCT** par la fusion d'une nervure abaissée et d'un
   cylindre, en plus des formes fermées. *Recommandé.*
