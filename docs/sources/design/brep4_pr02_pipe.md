# Palier 4 — PR 2 : repère mobile et balayage le long d'un chemin

| | |
|---|---|
| PR | branche `feat/brep-pipe` |
| Document de référence | [brep_palier4.md](brep_palier4.md), § 2 |
| Fichiers | `gbs/movingframe.h` (nouveau), `gbs-brep/pipe.h` (nouveau), `gbs-brep/sweep.h`, `gbs-brep/edit.h`, `gbs-brep/brep`, `python/gbsBindBrepOps.cpp`, `tests/tests_movingframe.cpp` (nouveau), `tests/tests_brep_pipe.cpp` (nouveau), `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/examples/gbs_examples.py`, `docs/sources/examples/gallery.py`, `python/tests/test_examples.py`, `python/tests/test_brep_ops.py`, `docs/sources/design/brep_palier4.md`, `News.md`, cette note |
| Statut | Validé (§ 8, 17 nov. 2026) ; revue de code § 9 |

## 1. Objectif

Balayer un profil (une arête, un contour ou une face) le long d'un chemin (une
arête ou un contour G1) : `make_pipe`. C'est la première des trois opérations du
palier, et celle qui ne dépend d'aucune autre. Le repère qui porte le profil
vient de `gbs/` : il servira aussi aux congés (PR 5, déplacement de l'arc le long
de la ligne des centres).

## 2. Sources

- Rotation minimale par double réflexion : W. Wang, B. Jüttler, D. Zheng,
  Y. Liu, *Computation of rotation minimizing frames*, ACM Transactions on
  Graphics 27(1), 2008.
- Repère à rotation minimale pour le balayage : F. Klok, *Two moving
  coordinate frames for sweeping along a 3D trajectory*, CAGD 3, 1986 ;
  R. Bishop, *There is more than one way to frame a curve*, Amer. Math.
  Monthly 82, 1975.
- *Skinning* et surfaces balayées : L. Piegl, W. Tiller, *The NURBS Book*,
  2e éd., § 10.3 et § 10.4.
- Repli d'une surface canal : condition classique de régularité
  1 − κ ⟨q, n⟩ > 0 (vitesse le long du chemin d'un point du profil).

OpenCascade ne sert que de comparaison en boîte noire dans les tests
(`BRepOffsetAPI_MakePipeShell`).

## 3. Repère mobile (`gbs/movingframe.h`)

`MovingFrame<T>(courbe, u1, u2, loi, angle_tol, cassures, fermé, binormale)`
donne en chaque u un repère orthonormé direct : o sur la courbe, z la tangente
unitaire, x et y dans le plan normal. Les lois (`FrameLaw`) :

- **Rotation minimale** (par défaut) : double réflexion de station en station.
  Une réflexion dans le plan médiateur de la corde, puis une deuxième qui ramène
  la tangente réfléchie sur la tangente suivante ; ordre 4. Les stations sont
  doublées (8 par morceau au départ) jusqu'à ce que le repère, aux stations
  communes à deux passes, bouge de moins de `angle_tol`. Entre deux stations,
  une double réflexion depuis la station précédente : le repère est continu et
  défini partout, sans interpolation d'angle. Au-delà de 65 536 stations sans
  atteindre `angle_tol`, `std::runtime_error` : c'est une limite de précision,
  pas un refus géométrique. `angle_tol` doit être positif
  (`std::invalid_argument`).
- **Chemin fermé** : le repère ne revient pas sur lui-même en général (holonomie
  de la courbe). L'angle résiduel φ est réparti le long du chemin, linéairement
  en longueur des cordes : x tourne de −φ w(u) autour de z, avec w = 0 au départ
  et w = 1 à l'arrivée. `closure_angle()` le donne.
- **Frenet** : x la normale principale. Refusé (`std::domain_error`) où la
  courbure s'annule (contrôlé sur 32 points par morceau), sur une courbe sans
  dérivées secondes, et là où la binormale se retourne.
- **Binormale fixe** : y la projection unitaire d'une direction donnée sur le
  plan normal. Refusé où la tangente est le long de cette direction, et là où y
  se retourne.

**Continuité.** Une inflexion, ou un passage de la tangente par la binormale
fixe, peut tomber entre deux des 32 points : le repère s'y retourne sans que la
courbure ou l'angle mesurés aux points le montrent. Le vecteur qui porte la loi
(binormale de Frenet, binormale fixe projetée) est donc suivi d'un point au
suivant : l'intervalle est coupé en deux tant que ce vecteur y tourne de plus de
π/8, et un retournement qui reste après 40 coupes est refusé. Une rotation
rapide mais continue (hélice de dix spires en un morceau) passe.

Les cassures (jonctions des arêtes d'un contour, seulement G1) sont des stations.
À une cassure, le repère est celui du morceau qui suit ; celui de la fin du
morceau qui précède est la limite à gauche, en `std::nextafter(u, u1)`. Frenet
est refusé si la binormale se retourne d'un côté à l'autre d'une cassure (deux
arcs qui tournent en sens opposés).

## 4. Balayage (`gbs-brep/pipe.h`)

### 4.1 API

```cpp
using PipeFrame = gbs::FrameLaw;   // RotationMinimizing, Frenet, FixedBinormal

template <std::floating_point T>
struct PipeOptions
{
    T tol{brep_default_tolerance<T>};   // faces, couture, coins du chemin
    T approx_tol{T(1e-5)};              // surfaces approchées : écart au balayage exact
    PipeFrame frame{PipeFrame::RotationMinimizing};
    point<T, 3> binormal{};             // pour FixedBinormal
};

auto make_pipe(Model<T> &, const ShapeId &profile, const ShapeId &path, PipeOptions<T> = {})
    -> BuildResult<SweepResult<T>>;
```

Le profil est pris tel quel : le repère au début du chemin est la référence, et
le mouvement rigide qui porte le repère de 0 à s porte le profil. Il n'a donc
pas à être dans le plan normal, ni à toucher le chemin.

Le chemin est vu comme une seule courbe : l'arête k du contour sur [k, k + 1],
dans le sens du contour (une arête renversée est parcourue de u2 à u1). Les faces
d'une arête du chemin sont balayées par le repère de cette arête, jusqu'à sa fin
prise en limite à gauche : un coin toléré (§ 4.5) ou un saut du repère à une
jonction ne passe pas d'une arête à la suivante.

Codes d'erreur : `InvalidId` pour un profil qui n'est ni arête, ni contour, ni
face (vérifié avant toute mesure) et pour un chemin qui n'est ni arête ni
contour ; `InvalidBounds` pour `tol` ou `approx_tol` non positifs ;
`NonGeneric` pour les refus géométriques du § 4.5 ; `PCurveApproximation` pour
une précision inatteignable (repère ou surface approchée).

### 4.2 Surfaces

Pour chaque arête du profil et chaque arête du chemin, une face. Le cas exact
est reconnu sur le **mouvement**, pas sur la loi du repère :

| Arête du chemin | Mouvement le long d'elle | Surface |
|---|---|---|
| segment (reconnu sur ses points, à `tol`) | translation (contrôlée en 4 points) | `prism_surface` de l'arête placée au début du segment : plan, cylindre, extrusion NURBS exacte |
| arc de cercle (reconnu sur ses points, à `tol`) | rotation autour de l'axe du cercle (contrôlée en 4 points) | `revolution_surface` de l'arête placée : cylindre, cône, sphère, tore, révolution NURBS exacte |
| autre | — | B-spline approchée (§ 4.3) |

L'arc est reconnu géométriquement, comme le segment, quel que soit le type de
sa courbe : le cercle passe par trois des 49 points de l'arête, tous à `tol`
de lui (distance au centre et au plan), et le chemin le parcourt dans un seul
sens. Un arc donné en B-spline rationnelle (STEP, `to_nurbs`) donne donc le même
tore qu'un `EllipseCurve`. L'angle de rotation aux points de contrôle est
mesuré sur la position du point, pas sur le paramètre, qui n'est pas l'angle
pour une NURBS.

Ainsi la rotation minimale et Frenet donnent des surfaces exactes le long d'un
arc, et la binormale fixe aussi quand elle est l'axe de l'arc. Sur un chemin
fermé non plan, la correction de l'holonomie fait tourner le repère le long des
arcs : le contrôle le voit et passe à l'approximation. Le plan prévoyait les
seuls segments et cercles comme profils le long d'un arc ; tout profil y est
exact, comme dans `make_revolution`.

### 4.3 Surface approchée

L'arête du profil, en NURBS (`with_nurbs` de `sweep.h`, exacte pour les
primitives et les B-splines), est placée par le repère en n + 1 stations de
l'arête du chemin, puis passée par *skinning* (`loft`, cubique en travers, les
v des stations imposés). Le profil placé garde ses nœuds et ses poids (mouvement
rigide), le *skinning* ne fait donc qu'interpoler la trajectoire de chaque pôle.

L'écart est mesuré entre les stations (aux quarts, aux milieux et aux trois
quarts), en 5 points par intervalle de nœuds du profil, à la surface exacte :
le profil placé par le repère au même paramètre. n part de 4 et double tant
que l'écart dépasse `approx_tol`, jusqu'à 4096 ; au-delà, refus. Les surfaces
voisines partagent exactement le profil placé à leur jonction, et les bouchons
sont des copies exactes du profil : la couture se fait à `tol` (à l'écart de
l'approximation du profil près, quand il n'est ni primitive ni B-spline).

La précision du repère à rotation minimale suit : `angle_tol = approx_tol / ρ`,
où ρ est la plus grande distance du profil au début du chemin. Le repère
approché ne déplace donc pas le profil de plus de `approx_tol`.

Les placements du profil aux stations et entre elles sont calculés une fois par
arête du chemin et par nombre de stations, puis partagés par toutes les arêtes
du profil. Le placement de départ est inversé une fois pour toutes, et la nature
de la surface (prisme, révolution ou B-spline) est décidée une fois par arête du
chemin.

### 4.4 Topologie

Celle de `make_prism` et `make_revolution`, dont `make_pipe` réutilise la
construction (`sweep_detail::sweep`, généralisée à plusieurs surfaces par arête
du profil) :

- une arête le long d'une arête : une face ; le long d'un contour : une coque ;
- un contour : une coque ;
- une face : un solide, fermé par le profil et sa copie placée au bout ; les
  trous du profil sont des vides le long d'un chemin fermé ;
- chemin fermé : pas de bouchons.

Aux jonctions (et à la fermeture), le profil placé à la fin d'une arête et au
début de la suivante peut différer d'un écart toléré (§ 4.5) : la couture se
fait au moins à cet écart, mesuré (`sweep_detail::sweep` prend une tolérance
de couture en plus).

En échec, le modèle est rendu tel qu'il était.

`gbs-brep/edit.h` : la mise en place d'une courbe par un placement de
`copy_shape` devient une fonction libre, `moved_curve`, que `make_pipe` utilise.

### 4.5 Refus (`NonGeneric`)

- **Coin** : à une jonction du chemin (fermé : aussi entre la fin et le début),
  l'angle des tangentes fois ρ dépasse la tolérance de la jonction, le plus
  grand de `tol` et de la tolérance de son sommet (un contour lu en STEP a
  souvent des sommets plus larges que `tol`). C'est l'écart que le coin
  laisserait entre deux faces voisines.
- **Saut du repère à une jonction** : le profil placé à la fin d'une arête
  (limite à gauche) et au début de la suivante s'écartent de plus de la
  tolérance de la jonction. C'est le cas de Frenet entre deux arcs de plans
  différents (la normale principale tourne d'un coup) ; un coin toléré passe.
- **Repli** : en 65 points par arête du chemin et pour chaque point q du profil
  placé, la vitesse de q le long du chemin a pour composante tangente
  1 − κ⃗ · (q − c) ; elle doit rester positive (seuil 10⁻⁹). Au-delà, le profil
  dépasse le rayon de courbure et la surface se replie.
  Sur un chemin sans dérivées secondes (une arête d'une face décalée
  d'une B-spline), la courbure vient de différences de la tangente.
- **Auto-intersection** : le balayage se recoupe loin des replis (spires d'une
  hélice qui se chevauchent, chemin qui repasse près de lui-même). Les
  polylignes du profil (un segment par arête droite, 16 sinon) sont placées à
  des stations espacées d'environ ρ/2 le long du chemin, au plus 200 000
  quadrilatères ; deux quadrilatères dont les stations sont à plus de π ρ (plus
  deux pas) l'une de l'autre le long du chemin (autour de lui s'il est fermé) ne doivent pas se
  couper (test segment-triangle, grille de cellules de la taille du plus grand
  quadrilatère). Plus près, c'est le repli qui décide. Un recouvrement plus
  petit que l'erreur des polylignes (environ 2 % de ρ pour un cercle) peut
  échapper.
- **Repère indéfini** : Frenet où la courbure s'annule ou la binormale se
  retourne, sans dérivées secondes ; binormale fixe le long de la tangente (une
  binormale nulle : `InvalidBounds`).
- **Chemin fermé** sur lequel le repère ne ramène pas le profil sur lui-même à
  la tolérance du sommet de fermeture près (Frenet sur une courbe fermée
  seulement G1, par exemple).

## 5. Python

`gbs.brep.make_pipe(model, profile, path, options)`, `PipeOptions(tol,
approx_tol, frame, binormal)`, `PipeFrame` ; `BRepError` sur les refus.

## 6. Galerie et exemples

Section « Sweeps » de la galerie (`gallery.py`), cas dans `gbs_examples.SWEEPS` :

| Cas | Faces | Volume |
|---|---|---|
| disque le long d'un coude (segment, quart de cercle, segment) | 5 : 2 cylindres, 1 tore, 2 disques | π r² L (Pappus) |
| carré le long d'une hélice, rotation minimale | 6 : 4 B-splines, 2 carrés | aire × L |
| le même, binormale fixe selon z (le carré reste à plat) | 6 | aire × L |
| disque autour d'un stade fermé (deux demi-cercles, deux segments) | 4 : 2 tores, 2 cylindres, pas de bouchon | π r² L |

`test_examples.py::test_sweeps` vérifie validité, volume (2·10⁻³) et maillage
fermé.

## 7. Tests

- `tests_movingframe` : rotation minimale sur une hélice (l'angle avec Frenet
  suit −τ s à 10⁻⁷), sur une courbe plane (la normale du plan reste fixe dans
  le repère), sur un trèfle fermé (holonomie non nulle, repère refermé à 10⁻⁸,
  continu) ; précision atteinte contre un repère fin ; cassures sans effet sur
  une courbe lisse, une cassure donnée deux fois comptée une fois ; refus ;
  retournements entre deux points (inflexion d'une B-spline en S, S de deux arcs
  à la cassure, tangente par la binormale fixe), hélice de dix spires en un
  morceau acceptée ; courbe sans dérivées secondes ; `angle_tol` nul, NaN ou
  inatteignable.
- `tests_brep_pipe` :
  - disque le long d'un cercle : un tore exact, 2π² R r² ; Frenet aussi ;
    quart de cercle avec deux bouchons ;
  - carré le long d'un segment : une boîte ; disque le long d'un coude :
    cylindres et tore exacts, π r² L ; profil contour ou arête : coque ;
  - carré et disque le long d'une hélice B-spline : faces B-splines, écart sous
    `approx_tol`, aire × L ; le disque donne un tube dont tous les points sont à
    r du chemin à 1,1·10⁻⁵ près ; binormale fixe ;
  - stade fermé : tores et cylindres exacts, sans bouchon ;
  - courbe gauche fermée : repère refermé, tube à r près, Frenet aussi ;
  - refus : coin, repli, Frenet sur un segment, binormale le long de la
    tangente ou nulle, chemin qui n'est ni arête ni contour, profil solide ou
    sommet (`InvalidId`), tolérances non positives (`InvalidBounds`), précision
    du repère inatteignable (`PCurveApproximation`) ;
  - modèle inchangé après un échec survenu une fois une face construite (un
    contour dont la deuxième arête, une tente C0, ne s'approche pas à 10⁻¹²) ;
  - auto-intersection : carré le long d'une hélice de pas 0,2 refusé, de pas
    0,5 accepté ; anneau le long d'un huit qui se croise refusé ;
  - jonctions : coin de 10⁻⁸ rad, 10 plans exacts ; S de deux arcs, deux tores
    exacts en rotation minimale, Frenet refusé ; arcs tangents dans des plans à
    30°, Frenet refusé, deux tores en rotation minimale ; stade dont la dernière
    arête finit à 5·10⁻⁷ de la première, sommet à 10⁻⁶ : fermé et valide ;
  - chemin quelconque : coude dont l'arc est une B-spline rationnelle, tore
    exact ; arête d'une face décalée (sans dérivées secondes), balayée sans
    exception, Frenet refusé.
- `tests_brep2_occt.pipes` : coude et hélice contre
  `BRepOffsetAPI_MakePipeShell` (Frenet corrigé) : mêmes volumes, mêmes
  nombres de faces et de solides, résultat de gbs lu par OCCT valide.

## 8. Points validés

1. **Cas exact reconnu sur le mouvement** (translation le long d'un segment,
   rotation autour de l'axe le long d'un arc), quelle que soit la loi du repère ;
   tout profil le long d'un arc donne une révolution exacte.
2. **Repère entre les stations par une double réflexion** depuis la station
   précédente, plutôt que l'interpolation de l'angle prévue au § 2.2 du plan :
   même précision, pas d'angle de référence à choisir.
3. **Écart mesuré sur la surface**, entre les stations, au profil placé par le
   même repère ; précision du repère `approx_tol / ρ`.
4. **Coin toléré** quand l'angle fois la taille du profil reste sous `tol` ;
   refus au-delà.
5. **Repli refusé** dès que 1 − κ⃗ · (q − c) ≤ 10⁻⁹ en un point du profil :
   un tore à collier nul (R = r) est donc refusé, comme la surface qui se
   replie.
6. **Auto-intersection** détectée sur des polylignes placées à des stations
   espacées de ρ/2, les paires à moins de π ρ le long du chemin laissées au
   repli (§ 4.5) : un recouvrement plus petit que l'erreur des polylignes peut
   échapper ; la tolérance d'un contact entre deux parties du balayage n'est pas
   celle de `tol`.
7. **Tolérance d'une jonction** : le plus grand de `tol` et de celle de son
   sommet, pour le coin, le saut du repère et la fermeture ; la couture suit
   l'écart mesuré.

## 9. Revue de code

Revue de la PR (14 points), tous traités dans cette PR :

| # | Point | Traitement |
|---|---|---|
| 1 | Auto-intersection globale non détectée : un carré 0,4 le long d'une hélice de pas 0,2 donnait un solide déclaré valide, le volume recouvert compté deux fois | contrôle par polylignes placées (§ 4.5), refus `NonGeneric` ; tests hélice et huit |
| 2 | Frenet discontinu à une jonction G1 (S de deux arcs) : *skinning* jusqu'à 4096 stations, puis `PCurveApproximation` | `MovingFrame` refuse le retournement de la binormale à une cassure ; `make_pipe` refuse un saut du repère au-delà de la tolérance de la jonction (arcs de plans différents) |
| 3 | Inflexion de Frenet, ou tangente par la binormale fixe, entre deux des 32 points : même échec lent | suivi du vecteur de la loi, intervalles coupés en deux (§ 3) |
| 4 | Coin toléré de 10⁻⁸ rad : les plans de l'arête qui précède devenaient des B-splines (repère pris sur l'arête suivante) | fin de chaque arête en limite à gauche ; couture à l'écart mesuré |
| 5 | Chemin sans dérivées secondes (arête d'une face décalée) : `std::runtime_error` sortait de `make_pipe` | courbure par différences de la tangente ; `PathCurve` transmet `has_second_derivatives` ; Frenet refusé |
| 6 | `failure_leaves_the_model` échouait avant toute face : le retour arrière n'était pas testé | profil dont la deuxième arête échoue une fois la première face faite |
| 7 | Précision inatteignable du repère rendue `NonGeneric`, `approx_tol` nul ou NaN non refusé | `std::runtime_error` → `PCurveApproximation` ; tolérances validées (`InvalidBounds`) |
| 8 | Arc reconnu par son type seulement (`EllipseCurve`) | reconnaissance géométrique, angle mesuré sur la position (§ 4.2) |
| 9 | Type du profil vérifié après les mesures : solide → `NonGeneric`, sommet → `EmptyInput` | vérifié d'entrée, `InvalidId` |
| 10 | Fermeture et coins mesurés à `tol`, pas à la tolérance des sommets | tolérance de la jonction (§ 4.5, point 7 du § 8) |
| 11 | `motion` réévaluait le repère au départ et son inverse à chaque appel ; placements refaits pour chaque arête du profil | `Motion` (inverse calculé une fois), `StationCache` par arête du chemin et niveau, nature de la surface décidée une fois par arête du chemin |
| 12 | Conversion OCCT et aides de test réécrites | `occt_utils::BSplineCurve`, `gbs::length`, `extrema_curve_point` |
| 13 | Doc de `SweepResult::shape` : une arête le long de plusieurs arêtes donne une coque | mise à jour |
| 14 | `j` redondant, nœuds échantillonnés deux fois, lambda `unit` dupliqué, cassures en double | simplifiés ; cassures dédoublonnées (`std::ranges::unique`) |
