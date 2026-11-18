# Palier 4 — PR 4 : décalage de solide et coques épaisses

| | |
|---|---|
| PR | branche `feat/brep-offset-solid` |
| Document de référence | [brep_palier4.md](brep_palier4.md), § 3 |
| Fichiers | `gbs-brep/thick.h` (nouveau), `gbs-brep/brep`, `python/gbsBindBrepOps.cpp`, `tests/tests_brep_thick.cpp` (nouveau), `tests/tests_brep_sew.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/examples/gbs_examples.py`, `docs/sources/examples/gallery.py`, `python/tests/test_brep_ops.py`, `python/tests/test_examples.py`, `docs/sources/design/brep_palier4.md`, `docs/sources/design/brep4_pr03_offset_surfaces.md`, `News.md`, cette note |
| Statut | À valider (§ 9) |

## 1. Objectif

Deux opérations, sur les surfaces décalées de la PR 3 :

- `offset_solid` : le solide décalé d'une distance d, à joints vifs ;
- `make_thick` : le solide évidé à une épaisseur t, ouvert par des faces
  retirées.

Les deux suivent le § 3 du plan, questions 8 (joints vifs) et 9 (coque par
booléen).

## 2. Sources

- Décalage de solides à joints vifs, par intersection des surfaces décalées
  voisines : J. R. Rossignac, A. A. G. Requicha, *Offsetting operations in
  solid modelling*, CAGD 3, 1986 ; T. Maekawa, *An overview of offset curves
  and surfaces*, CAD 31(3), 1999, § 5.
- Coques épaisses par différence d'un solide et de son décalé intérieur :
  Rossignac et Requicha, *ibid.*, § 4 (*shelling*).
- Prolongement d'une B-spline par le polynôme de son dernier intervalle :
  `extrapolated_v` (`gbs/bsstools.h`, blossoms de Ramshaw).

OpenCascade ne sert que de comparaison en boîte noire dans les tests
(`BRepOffsetAPI_MakeOffsetShape` et `BRepOffsetAPI_MakeThickSolid`, joints
`GeomAbs_Intersection`).

## 3. API

```cpp
template <std::floating_point T>
struct OffsetOptions
{
    T tol{brep_default_tolerance<T>};  // tolérance du résultat ; arêtes tangentes en deçà
    IntersectOptions<T> intersect{};   // intersections des surfaces décalées
    T max_extension{T(100)};           // plus grand prolongement d'une face au-delà d'une arête vive, en |d|
};

auto offset_solid(Model<T> &, SolidId, T d, const OffsetOptions<T> & = {}) -> BuildResult<SolidId>;
auto make_thick(Model<T> &, SolidId, std::span<const FaceId> removed, T t, const OffsetOptions<T> & = {})
    -> BuildResult<SolidId>;
```

- **Signe de d** : positif, le solide grossit le long de sa normale
  extérieure ; négatif, il rétrécit.
- **Signe de t** : positif, la coque est creusée vers l'intérieur ; négatif,
  elle est ajoutée vers l'extérieur.
- **Sans face retirée** : `make_thick` rend un solide creux fermé, avec un
  vide.

## 4. Décalage (`offset_solid`)

Toutes les faces du solide (coque extérieure et vides) prennent une distance
le long de leur normale extérieure. La topologie reste celle du solide, sur de
nouvelles entités.

### 4.1 Arêtes : tangentes ou vives

L'arête est échantillonnée en 9 points. En chaque point, on prend les normales
extérieures n₁ et n₂ des deux faces, et leurs distances d₁ et d₂ :

- **Arête tangente** : |d₁ n₁ − d₂ n₂| ≤ `tol` en tous les points. Cela couvre
  les faces tangentes à la même distance, une couture, ou deux faces retirées
  par `make_thick`. Les décalés des deux faces se touchent le long d'une même
  courbe. La nouvelle arête est la pcurve de l'arête sur la surface décalée
  (`CurveOnSurface`), prolongée au-delà de ses bouts si l'arête s'allonge. Le
  prolongement suit le polynôme du dernier intervalle de la pcurve. L'écart
  entre les deux décalés est contrôlé le long de la nouvelle arête.
- **Arête vive** : aucun point tangent. La nouvelle arête est l'intersection
  des deux surfaces décalées (`intersect` : formes fermées entre primitives,
  cheminement sinon), sur les boîtes prolongées des faces. On garde la branche
  qui passe à 10 `tol` des deux nouveaux sommets ; s'il y en a plusieurs, celle
  dont le milieu est le plus proche de celui de l'arête d'origine. Pour une
  arête fermée (un cercle), la boucle entière, reparamétrée depuis son sommet.
- **Arête en partie tangente** : refusée.

### 4.2 Prolongement des faces

Au-delà d'une arête vive d'angle a entre les normales, la nouvelle arête
s'écarte de l'ancienne de |d₂ − d₁ cos a| / sin a le long de la première face.
Chaque face est prolongée du double du plus grand de ces écarts, plus 10 `tol`,
car le sommet de trois faces va plus loin que ses arêtes :

- **Primitive** : sur une boîte plus large. Un paramètre périodique est limité
  à une période, une sphère à ses pôles, un cône s'arrête avant son sommet.
- **B-spline** : prolongée par le polynôme de ses intervalles de bord, en v
  (`extrapolated_v`) puis en u (par transposition).
- **Autres surfaces** (révolution, surface déplacée) : laissées telles quelles.

Un écart au-delà de `max_extension` |d| est refusé.

### 4.3 Sommets

Un sommet est le point commun des surfaces décalées autour de lui :

- **Toutes ses arêtes tangentes** : la surface décalée de n'importe laquelle de
  ses faces, au même paramètre. C'est le cas des pôles d'une sphère et du
  sommet unique d'un tore. L'accord entre les faces est contrôlé.
- **Une arête tangente et des vives** : le point de la pcurve prolongée de
  l'arête tangente où les autres faces la rencontrent. Le système est résolu
  par Gauss-Newton sur le paramètre de l'arête et les (u, v) des autres faces.
  Les autres arêtes tangentes du sommet doivent passer par ce point, à
  `tol` / 10 près : par exemple la couture d'un cône et celle d'un disque.
- **Toutes ses arêtes vives** : on cherche x et les (u, v) de chaque face par
  Gauss-Newton. Le départ est la solution aux moindres carrés des plans
  tangents décalés. Trois faces donnent un système carré. Au-delà de trois, les
  faces doivent s'accorder à `tol` / 10 près (le sommet d'une pyramide
  régulière, par symétrie ; cas non testé) ; sinon, refus.

### 4.4 Contrôles

- **Face qui disparaît** : une arête tangente dont les paramètres s'inversent,
  une arête vive retournée par rapport à l'ancienne ou réduite à un point.
  Exemple : une boîte de 1 décalée de −0,6.
- **Résultat** : `check` (validité, volume positif).
- **Auto-intersection** : deux faces qui ne partagent aucun sommet et se
  rencontrent (section de chaque paire, ses débris effacés). Exemple : les deux
  bouts d'un anneau fendu décalé vers l'extérieur se croisent.

En échec, le modèle est rendu tel qu'il était.

## 5. Coque épaisse (`make_thick`)

`make_thick(s, retirées, t)` décale le solide de −t, mais les faces retirées
de 0 (elles restent en place). Le résultat est `cut(s, intérieur)`, ou
`cut(extérieur, s)` pour t < 0 ; le décalé intermédiaire est effacé.

Les faces retirées sont coïncidentes dans les deux opérandes : le booléen des
faces coïncidentes les ôte, et la coque s'ouvre là. Une arête entre deux faces
retirées est tangente au sens du § 4.1 (d₁ = d₂ = 0) : elle reste en place, et
un sommet entre deux faces retirées et une paroi est le point où cette arête
rencontre la paroi décalée.

## 6. Refus (`NonGeneric`)

- **Repli** d'une surface décalée (`offset_surface`). Le sommet d'un cône tombe
  dans ce cas, puisque sa boîte atteint l'apex.
- **Face qui disparaît**, et **auto-intersection** (§ 4.4).
- **Sommet** où les surfaces décalées ne se rencontrent pas en un point :
  - plus de trois faces vives qui ne s'accordent pas ;
  - une arête tangente fermée qui passe par un sommet vif ;
  - des arêtes tangentes qui ne passent pas par le même point.
- **Arête** :
  - en partie tangente ;
  - tangente entre deux faces à des distances différentes (une face retirée
    tangente à une paroi) ;
  - d'une face avec elle-même ;
  - vive, sans branche d'intersection qui joigne ses sommets.
- **Joint trop aigu** pour être prolongé (`max_extension`).
- **Résultat invalide** pour `check`.

`InvalidId` : solide mort, face retirée qui n'est pas du solide.
`InvalidBounds` : tolérance non positive, distance non finie ou sous la
tolérance.

## 7. Python, galerie

- **Python** : `gbs.brep.offset_solid(model, solid, distance, options)`,
  `make_thick(model, solid, removed, thickness, options)` et
  `OffsetOptions(tol, max_extension)`.
- **Galerie** : section « Offset and thick solids » (`gbs_examples.THICK`), où
  les arêtes du solide d'origine sont en rouge :
  - une boîte décalée de 0,25 ;
  - une dalle arrondie décalée de 0,2 (arêtes tangentes) ;
  - un carré balayé le long d'une hélice, décalé de 0,05 (faces B-splines) ;
  - une boîte évidée à 0,1 et ouverte en haut ;
  - un gobelet.

  `test_examples.py::test_thick` vérifie la validité, le maillage fermé et le
  volume exact à 10⁻³ près.

## 8. Tests

Les volumes sont mesurés sur le maillage (déflexion 5·10⁻⁵, à 3·10⁻⁴ près),
car `signed_volume` se trompe de 0,4 % sur une face percée (signalé à part).

- **`tests_brep_thick`** :
  - boîte : ± d, 6 plans, volume exact ; −0,6 refusé, modèle inchangé ;
  - cylindre : coutures et cercles ;
  - sphère : pôles ; −1,2 replié ;
  - tore : un sommet, deux coutures ;
  - tronc de cône : cône exact, son rayon décalé de d / cos a ;
  - dalle arrondie : arêtes tangentes plan/cylindre ;
  - carré le long d'une hélice : 4 faces B-splines, arêtes vives par
    cheminement ;
  - coques épaisses : boîte ouverte vers l'intérieur et vers l'extérieur,
    boîte creuse fermée avec un vide, gobelet ;
  - refus : sommet d'un cône, arguments ;
  - auto-intersection : anneau fendu.
- **`tests_brep2_occt.offsets`** : boîte, cylindre et dalle (± d), carré le
  long d'une hélice (± 0,05), boîte et gobelet évidés, contre
  `BRepOffsetAPI_MakeOffsetShape` et `BRepOffsetAPI_MakeThickSolid` (joints
  `GeomAbs_Intersection`) : mêmes volumes, mêmes nombres de faces, résultat de
  gbs lu par OCCT valide.
- **Python** : `test_brep_ops.py::test_offset_solid_and_thick` et
  `test_examples.py::test_thick`.
- `tests_brep_sew` qualifie un appel `brep::detail::` : `gbs-brep/brep` inclut
  maintenant `gbs/bsstools.h` (prolongement des B-splines), qui définit
  `gbs::detail`, ambigu sous `using namespace gbs; using namespace gbs::brep;`.

## 9. Points à valider

1. **Classement d'une arête** par |d₁ n₁ − d₂ n₂| ≤ `tol` en 9 points (tangente
   sinon vive), une arête en partie tangente refusée.
2. **Arête vive par l'intersection existante** : formes fermées entre
   primitives, cheminement sinon. La branche est choisie par les deux nouveaux
   sommets (à 10 `tol`), puis par la proximité de l'arête d'origine.
3. **Prolongement des faces** : deux fois l'écart du joint, plus 10 `tol`,
   plafonné à `max_extension` |d| (100 par défaut). Seules les primitives et
   les B-splines sont prolongées.
4. **Sommets** : sur l'arête tangente quand il y en a une, sinon par
   Gauss-Newton sur toutes les faces. Plus de trois faces vives doivent
   s'accorder à `tol` / 10 près, sinon refus. Il n'y a pas de face de raccord
   au sommet, comme le prévoient les joints vifs.
5. **Auto-intersection détectée par sections** entre faces sans sommet commun.
   Une face qui se recoupe elle-même (deux spires d'une hélice qui se touchent)
   ou deux faces voisines qui se recoupent loin de leur arête commune ne sont
   pas vues.
6. **Sommet d'un cône refusé** (le décalé se replie à l'apex). CATIA le
   prolongerait jusqu'au sommet du cône décalé : à faire plus tard.
7. **`make_thick` par booléen** (question 9 du plan). Une face retirée tangente
   à une paroi est refusée.
