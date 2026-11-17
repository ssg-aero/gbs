# Palier 4 — PR 2 : repère mobile et balayage le long d'un chemin

| | |
|---|---|
| PR | branche `feat/brep-pipe` |
| Document de référence | [brep_palier4.md](brep_palier4.md), § 2 |
| Fichiers | `gbs/movingframe.h` (nouveau), `gbs-brep/pipe.h` (nouveau), `gbs-brep/sweep.h`, `gbs-brep/edit.h`, `gbs-brep/brep`, `python/gbsBindBrepOps.cpp`, `tests/tests_movingframe.cpp` (nouveau), `tests/tests_brep_pipe.cpp` (nouveau), `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/examples/gbs_examples.py`, `docs/sources/examples/gallery.py`, `python/tests/test_examples.py`, `python/tests/test_brep_ops.py`, `docs/sources/design/brep_palier4.md`, `News.md`, cette note |
| Statut | À valider (§ 8) |

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
  défini partout, sans interpolation d'angle.
- **Chemin fermé** : le repère ne revient pas sur lui-même en général (holonomie
  de la courbe). L'angle résiduel φ est réparti le long du chemin, linéairement
  en longueur des cordes : x tourne de −φ w(u) autour de z, avec w = 0 au départ
  et w = 1 à l'arrivée. `closure_angle()` le donne.
- **Frenet** : x la normale principale. Refusé (`std::domain_error`) où la
  courbure s'annule, contrôlé sur 32 points par morceau.
- **Binormale fixe** : y la projection unitaire d'une direction donnée sur le
  plan normal. Refusé où la tangente est le long de cette direction.

Les cassures (jonctions des arêtes d'un contour, seulement G1) sont des stations.

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
dans le sens du contour (une arête renversée est parcourue de u2 à u1).

### 4.2 Surfaces

Pour chaque arête du profil et chaque arête du chemin, une face. Le cas exact
est reconnu sur le **mouvement**, pas sur la loi du repère :

| Arête du chemin | Mouvement le long d'elle | Surface |
|---|---|---|
| segment | translation (contrôlée en 4 points) | `prism_surface` de l'arête placée au début du segment : plan, cylindre, extrusion NURBS exacte |
| arc de cercle | rotation autour de l'axe du cercle (contrôlée en 4 points) | `revolution_surface` de l'arête placée : cylindre, cône, sphère, tore, révolution NURBS exacte |
| autre | — | B-spline approchée (§ 4.3) |

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

### 4.4 Topologie

Celle de `make_prism` et `make_revolution`, dont `make_pipe` réutilise la
construction (`sweep_detail::sweep`, généralisée à plusieurs surfaces par arête
du profil) :

- une arête le long d'une arête : une face ; le long d'un contour : une coque ;
- un contour : une coque ;
- une face : un solide, fermé par le profil et sa copie placée au bout ; les
  trous du profil sont des vides le long d'un chemin fermé ;
- chemin fermé : pas de bouchons.

En échec, le modèle est rendu tel qu'il était.

`gbs-brep/edit.h` : la mise en place d'une courbe par un placement de
`copy_shape` devient une fonction libre, `moved_curve`, que `make_pipe` utilise.

### 4.5 Refus (`NonGeneric`)

- **Coin** : à une jonction du chemin (fermé : aussi entre la fin et le début),
  l'angle des tangentes fois ρ dépasse `tol`. C'est l'écart que le coin
  laisserait entre deux faces voisines.
- **Repli** : en 65 points par arête du chemin et pour chaque point q du profil
  placé, la vitesse de q le long du chemin a pour composante tangente
  1 − κ⃗ · (q − c) ; elle doit rester positive (seuil 10⁻⁹). Au-delà, le profil
  dépasse le rayon de courbure et la surface se replie.
- **Repère indéfini** : Frenet où la courbure s'annule, binormale fixe le long
  de la tangente (une binormale nulle : `InvalidBounds`).
- **Chemin fermé** sur lequel le repère ne ramène pas le profil sur lui-même à
  `tol` près (Frenet sur une courbe fermée seulement G1, par exemple).

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
  une courbe lisse ; refus.
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
    tangente ou nulle, chemin qui n'est ni arête ni contour ; modèle inchangé
    après un échec.
- `tests_brep2_occt.pipes` : coude et hélice contre
  `BRepOffsetAPI_MakePipeShell` (Frenet corrigé) : mêmes volumes, mêmes
  nombres de faces et de solides, résultat de gbs lu par OCCT valide.

## 8. Points à valider

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
