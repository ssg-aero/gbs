# BREP palier 1 — PR 11 : retrait de l'ancienne topologie (option B)

| | |
|---|---|
| PR | branche `refactor/brep-legacy-topology` |
| Document de référence | [brep_core.md](brep_core.md), section 2.1 et question 1 ; [note de la PR 10](brep_pr10_python.md), section 5 (option B retenue) |
| Fichiers | supprimés : `inc/topology/basetopo.h`, `vertex.h`, `edge.h`, `wire.h`, `tests/tests_topo.cpp`, `src/topology/edge.cpp` ; modifiés : `tests/tests_halfedgemesh.cpp`, `docs/sources/design/brep_core.md` (correction § 2.1), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Ne garder qu'une topologie dans gbs : le noyau BREP natif. Les classes
`BaseTopo`, `Vertex`, `Edge`, `Wire` de `inc/topology/`, dépréciées depuis la
PR 3, sont supprimées, après le portage de leur dernier utilisateur.

## 2. Inventaire des utilisateurs

| Fichier | Usage | Traitement |
|---|---|---|
| `tests/tests_topo.cpp` | tests de l'ancienne topologie | supprimé ; ses trois cas sont portés dans `tests_brep_wire` depuis la PR 3 |
| `tests/tests_halfedgemesh.cpp` | `Wire<T,2>` / `Edge<T,2>` pour les bords du mailleur Delaunay | **porté** sur `gbs::brep` (§ 3) |
| `tests/tests_topo_halfEdgeMeshData.cpp` | des noms `Vertex`, `Edge`… | **aucun changement** : ce sont des alias locaux vers `HalfEdgeVertex` et `HalfEdge` du maillage, pas l'ancienne topologie |
| `src/topology/edge.cpp` | instanciation explicite de `Edge<…>` | supprimé ; aucune cible CMake ne le compilait |

La note de la PR 10 comptait `tests_topo_halfEdgeMeshData.cpp` parmi les
utilisateurs : c'était une erreur de lecture, corrigée ici.

## 3. Portage des bords du mailleur

Le mailleur Delaunay 2D reçoit son bord sous forme de points `std::array<T,2>`.
Les trois générateurs du test (`make_boundary2d_1` carré, `_2` ellipse, `_3`
polygone non convexe) passaient par un `Wire<T,2>` 2D. Ils construisent
désormais un **wire `gbs::brep` dans le plan `z = 0`** :

- `make_polygon_wire(model, points2d)` : un segment `make_edge` par côté, puis
  `make_wire`, qui fusionne les sommets et ferme le contour ;
- l'ellipse est construite directement en 3D par `build_ellipse<T, 3>`, qui la
  place dans le plan `xy` ;
- `mesh_wire_uniform(model, wire, dm)` parcourt les co-arêtes dans le sens du
  wire, échantillonne chaque arête à pas `dm` en longueur d'arc et renvoie les
  coordonnées `(x, y)`, exactement comme l'ancienne fonction.

`mesh_hed_wire_uniform`, jamais appelée, est retirée.

**Équivalence vérifiée** : l'ancienne version du test, compilée depuis
`master`, et la nouvelle donnent le même nombre d'assertions (5727, qui dépend
du nombre de points de bord générés), toutes réussies.

Le BREP reste 3D seulement (question 4 du document) : un bord 2D est un wire
dans le plan `z = 0`. Aucune structure 2D dédiée n'est réintroduite.

## 4. Défaut existant rencontré

`add_dimension` appliqué à une courbe **rationnelle** (`gbs/bsctools.h`, et
probablement `gbs/bsstools.h` pour les surfaces) ajoute la nouvelle coordonnée
à la place du poids, en dernière position des pôles homogènes : une ellipse
rationnelle 2D relevée en 3D aurait un poids nul. Le portage l'évite en
construisant l'ellipse en 3D ; la correction fait l'objet d'une tâche séparée.

## 5. Tests

Construction complète (toutes les cibles, y compris le module Python) puis
`ctest` sur toutes les suites C++ ; `tests_halfedgemesh` : 5727 assertions,
identique à `master`.

## 6. Points à valider

1. **Bords du mailleur en wires `gbs::brep` à `z = 0`**, sans structure 2D dédiée.
2. **Suppression de `src/topology/edge.cpp`**, qui n'était compilé par aucune cible.
3. **Correction séparée** de `add_dimension` pour la géométrie rationnelle.
