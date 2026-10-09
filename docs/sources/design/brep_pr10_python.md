# BREP palier 1 — PR 10 : bindings Python, documentation, ancienne topologie

| | |
|---|---|
| PR | branche `feat/brep-python` |
| Document de référence | [brep_core.md](brep_core.md), sections 8.2 et 10 (PR 10), question 9 ; notes précédentes [PR 1](brep_pr01_model.md) à [PR 9](brep_pr09_iges.md) |
| Fichiers | `python/gbsBindBrep.cpp` (nouveau), `python/gbsbind.cpp` (appel du sous-module, `IgesWriter.add_brep`, `IgesExportReport`), `python/CMakeLists.txt`, `python/tests/test_brep.py` (nouveau), `docs/Doxyfile` (`INPUT` += `../gbs-brep`), `docs/index.rst` (section gbs-brep), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Rendre le noyau BREP utilisable depuis Python avec la même logique qu'en C++,
documenter ses en-têtes, et clore le palier 1. Le retrait de l'ancienne
topologie, prévu au plan, **n'est pas fait** : il dépasse ce qui était prévu
(§ 5) et demande une décision.

## 2. API Python : sous-module `gbs.brep`

```python
import pygbs.gbs as gbs
brep = gbs.brep

m = brep.Model()
faces = [brep.make_face(m, srf) for srf in surfaces]          # FaceId
rep = brep.sew(m, faces, tol=1e-6)                             # SewReport
solid = brep.make_solid(m, rep.shells[0])                      # SolidId
assert brep.check(m, solid).ok()
for e in brep.explore_edges(m, solid):
    print(m.edge(e).curve, brep.TopologyIndex(m, solid).faces_of(e))
w = gbs.IgesWriter(); w.add_brep(m, solid, "part"); w.write("part.igs")
```

| Catégorie | Contenu |
|---|---|
| énumérations | `ShapeType`, `Orientation` (+ `reverse`), `Issue` (noms en `snake_case` tirés de `to_string`) |
| identifiants | `VertexId` … `CompoundId` : `index`, `valid()`, `int()`, `==`, `<`, `hash`, `repr` ; `ShapeId` converti automatiquement (variant de pybind11), `shape_type` |
| entités | `Vertex`, `Edge`, `CoEdge`, `Wire`, `Face`, `FaceUse`, `Shell`, `Solid`, `Compound`, attributs lisibles et modifiables |
| modèle | `Model` : `vertex(id)` … `compound(id)` (copies), `set_vertex(id, v)` … `set_compound`, `ids(type)`, `alive`, `count`, `empty`, `erase`, `compact` (→ `IdRemap.map`), `append`, `repr` |
| requêtes | `explore_vertices` … `explore_compounds`, `TopologyIndex` (`faces_of`, `edges_of`, `coedges_of`, `face_of`, `shells_of`), `is_closed`, `is_chained`, `is_manifold`, `is_orientable`, `free_edges`, `bounding_box`, `signed_volume`, `uv_signed_area`, `surface_closure` |
| builders | `make_vertex`, `make_edge` (5 formes), `make_degenerate_edge`, `make_wire`, `make_face` (surface entière, ou surface + wires), `make_face_from_pcurves`, `MakeFaceOptions`, `sew` (→ `SewReport`), `make_solid`, `make_compound` |
| validation | `check(model, shape, geometry=True)` → `CheckReport` (`ok()`, `bool`, `entries`, `count`, `has`) |
| IGES | `IgesWriter.add_brep(model, shape, name, approx_tol, scale)` → `IgesExportReport` |

## 3. Choix d'architecture et justification

### 3.1 Sous-module `gbs.brep`, `T = double`, sans suffixe `_3d`

Conforme à la question 9 du document : un espace de noms à part évite la
collision avec les noms de la géométrie (`gbs.Line`, `gbs.Wire` potentiels) et
le BREP est 3D par construction.

### 3.2 Exceptions Python plutôt que `std::expected`

En C++, les builders renvoient `BuildResult<Id>` (PR 3). En Python, un échec
lève `gbs.brep.BRepError` (sous-classe de `RuntimeError`) dont le message
contient le code d'erreur et le détail (`"wire not closed: …"`). C'est
l'idiome Python ; la garantie forte reste vraie (le modèle est inchangé). Les
accès à un identifiant mort lèvent la même exception.

### 3.3 Accesseurs par copie

`m.edge(id)` renvoie une **copie** de l'entité, et `m.set_edge(id, edge)` l'écrit.
Rendre une référence Python vers une case de l'arène serait dangereux : les
tables sont des `std::vector`, et toute création d'entité peut les réallouer ;
une référence conservée en Python deviendrait pendante et pourrait faire
planter l'interpréteur. La copie est bon marché (agrégats de quelques champs,
la géométrie est partagée par `shared_ptr`). Conséquence à connaître :
`m.shell(s).faces[0].orient = …` ne modifie rien ; il faut modifier la copie
puis `m.set_shell(s, sh)`. Le test `test_ids_and_model` le documente.

### 3.4 Géométrie partagée avec les classes Python existantes

Les courbes et surfaces passent par les classes de base déjà exposées
(`Curve3d`, `Curve2d`, `Surface3d`, avec holder `shared_ptr`) : toute courbe ou
surface Python existante (`BSCurve3d`, `BSSurfaceRational3d`, `CurveOnSurface3d`,
`SurfaceOfRevolution`…) est acceptée, et `edge.curve` revient sous son type
réel.

### 3.5 `make_face_from_pcurves` : un nom distinct

En C++, `make_face(srf, outer_pcurves, holes)` est une surcharge. En Python, une
liste vide ou une liste de courbes rendrait la résolution des surcharges
ambiguë face à `make_face(srf, outer_wire, holes)` : la variante à boucles 2D
porte donc son propre nom.

## 4. Documentation

`docs/Doxyfile` lit désormais `../gbs-brep` ; `docs/index.rst` gagne une
section « Native BREP core » avec les sept en-têtes. Les notes d'architecture
restent en Markdown dans `docs/sources/design/` (Sphinx n'a pas d'extension
Markdown configurée ; les ajouter à la table des matières demanderait
`myst-parser`).

## 5. Ancienne topologie : retrait non fait, décision demandée

Le plan prévoyait de supprimer `inc/topology/{basetopo,vertex,edge,wire}.h` et
`tests/tests_topo.cpp` dans cette PR. Deux constats l'en empêchent :

1. **Ces classes ont d'autres utilisateurs que `tests_topo.cpp`**, que la
   recherche de la PR 3 avait manqués :
   - `tests/tests_halfedgemesh.cpp` construit des `Wire<T,2>` et `Edge<T,2>`
     **2D** pour discrétiser le bord du mailleur Delaunay (`mesh_wire_uniform`,
     `make_boundary2d_1/2/3`) ;
   - `tests/tests_topo_halfEdgeMeshData.cpp` utilise `Vertex<double,3>` et
     `Edge<double,3>` ;
   - `src/topology/edge.cpp` instancie `Edge<…>` explicitement (fichier qui
     n'est compilé par aucune cible CMake).
   Le BREP natif est 3D seulement : remplacer ces usages 2D, c'est soit porter
   les bords du mailleur en wires 3D dans le plan `z = 0`, soit garder une petite
   structure 2D dédiée au maillage.
2. **La suppression de fichiers suivis** a été refusée par le garde-fou de
   permissions de la session de développement dès la PR 3.

Les classes restent marquées `[[deprecated]]` (PR 3). Options :

| Option | Effet |
|---|---|
| A. Garder l'ancienne topologie comme outil 2D de bord de maillage, la renommer (ex. `MeshBoundary2d`) et retirer `[[deprecated]]` | aucun risque pour le mailleur, deux structures assumées pour deux rôles |
| B. Porter `tests_halfedgemesh` et `tests_topo_halfEdgeMeshData` sur `gbs::brep` (wires 3D à `z = 0`), puis supprimer les quatre en-têtes, `tests_topo.cpp` et `src/topology/edge.cpp` | une seule topologie ; PR dédiée d'environ 300 lignes |
| C. Laisser en l'état (déprécié) jusqu'au palier 2 | rien à faire maintenant |

Recommandation : **B** dans une PR dédiée, car l'API 2D de l'ancien `Wire` se
remplace par `make_wire` sur des segments à `z = 0` et `make_points` sur les
courbes des arêtes, comme le fait déjà `discretize_wire` dans `tests_brep_wire`.

## 6. Tests

`python/tests/test_brep.py` (6 tests) : identifiants et modèle (copie et
`set_*`), boîte de six faces → `sew` → `make_solid` → `check` → explorateur,
`TopologyIndex`, boîte englobante, export IGES et comptage des entités 144/142
dans le fichier ; face depuis un wire d'arêtes données dans le désordre et face
trouée depuis des boucles 2D ; faces naturelles et `surface_closure` ;
exceptions (`invalid tolerance`, `wire not closed`, `empty input`, identifiant
mort) ; `check` signalant `solid_not_outward`. Toute la suite Python passe
(31 tests) ; les suites C++ ne changent pas.

## 7. Points à valider

1. **Sous-module `gbs.brep`**, `T = double`, sans suffixe `_3d`.
2. **Exceptions `BRepError`** en Python à la place de `std::expected`.
3. **Accesseurs par copie** avec `set_*`, plutôt que des références vers l'arène.
4. **`make_face_from_pcurves`** comme nom distinct de la variante à boucles 2D.
5. **Ancienne topologie** : option A, B ou C du § 5.
