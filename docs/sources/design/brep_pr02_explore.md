# BREP palier 1 — PR 2 : explorateur, index d'adjacence, requêtes de base (`gbs-brep/explore.h`)

| | |
|---|---|
| PR | ssg-aero/gbs#114, branche `feat/brep-explore` (empilée sur #113) |
| Document de référence | [brep_core.md](brep_core.md), section 5 ; note précédente [brep_pr01_model.md](brep_pr01_model.md) |
| Fichiers | `gbs-brep/explore.h` (491 l.), `gbs-brep/brep` (agrégateur), `tests/tests_brep_helpers.h` (boîte partagée), `tests/tests_brep_explore.cpp`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Lire le modèle sans le modifier : descendre dans un shape par type
(`explore`), remonter d'une entité vers ce qui l'utilise (`TopologyIndex`),
et répondre aux questions de base sur un wire ou un shell (fermé, variété,
orienté de façon cohérente, arêtes libres) et sur l'encombrement
(`bounding_box`). Tout est `const` sur le `Model`. Ces briques servent au
sewing (PR 7), au solide (PR 8), à `check()` (PR 6) et à l'export IGES (PR 9).

## 2. Vue d'ensemble

```
explore<Sub>(m, shape)        descente DFS, tables « vu » par type, résultat sans doublon
        │
        ▼
TopologyIndex(m, root)        vecteurs denses indexés par handle, remplis depuis explore :
   edges_of(VertexId)            vertex  → edges
   coedges_of(EdgeId)            edge    → CoEdgeRef { wire, index }
   faces_of(EdgeId)              edge    → faces (uniques)
   face_of(WireId)               wire    → face propriétaire
   shells_of(FaceId)             face    → shells

shell_edge_uses(m, shell)     edge → [ EdgeUse { face, coedge, orient effectif } ]   (dégénérées ignorées)
   ├─ is_manifold     tous ≤ 2
   ├─ is_closed       tous == 2
   ├─ is_orientable   ≤ 2 et sens opposés quand 2
   ├─ free_edges      == 1          (triées)
   └─ non_manifold_edges > 2        (triées)

is_closed / is_chained (wire)  coedge_end(i) == coedge_start(i+1)
coedge_point(m, wire, i, t)    point 3D le long de la co-arête (sens appliqué)
bounding_box(m, shape, n)      sommets ± tol, n points par arête ± tol, grille n×n des faces à bornes naturelles ± tol
```

## 3. Choix d'architecture et justification

### 3.1 Explorateur : visiteur DFS avec tables de visite denses

```cpp
detail::Visitor<T> { const Model<T>& m; std::array<std::vector<uint8_t>, 7> seen; visit_handle<K, Sub>(Handle<K>, std::vector<Sub>&); }
```

- Une table « vu » par type, dimensionnée par `capacity<K>()` (indices denses,
  tombstones compris) : le dédoublonnage est une lecture d'octet, sans
  `unordered_set`. Coût O(taille des tables) en mémoire pour l'allocation
  initiale, O(nombre d'entités atteintes) en temps.
- La descente suit strictement la hiérarchie du document (§ 5.1) ; le
  branchement par type est un `if constexpr` sur `K`, le résultat est poussé
  quand `K == Sub::type`. Pas de remontée : `explore<ShellId>(m, face)` est
  vide, par construction.
- Ordre de découverte déterministe (DFS dans l'ordre des vecteurs) : les
  sommets d'une face sortent dans l'ordre de son wire, ce que les tests
  vérifient et dont l'export IGES profitera.
- Les entités mortes et les handles invalides sont ignorés silencieusement,
  pas d'exception : un modèle en cours de mutation reste explorable.
- Le shape racine est inclus s'il est du type demandé (`explore<SolidId>(m, solid)`
  renvoie le solide), comme `TopExp_Explorer` avec le type du shape lui-même.

### 3.2 Index remontant : vecteurs denses, résultats en `std::span`

- Les cinq tables sont des `std::vector<std::vector<Id>>` indexés par l'indice
  du handle, dimensionnés par `capacity<K>()`. Lecture O(1), retour par
  `std::span<const Id>` sans copie ; un indice hors table donne un span vide
  plutôt qu'une exception, pour pouvoir interroger un index partiel avec des
  handles extérieurs à sa racine.
- `CoEdgeRef { WireId wire; uint32_t index; }` localise une co-arête sans
  introduire de handle de co-arête dans le modèle (décision § 2.5 du document) ;
  c'est exactement ce dont le sewing aura besoin pour rediriger `coedges[i].edge`.
- `faces_of` est dédoublonné (`push_unique`, linéaire sur une liste de 1 à 3
  éléments) : une arête seam, utilisée deux fois par la même face, donne une
  face et deux `CoEdgeRef`.
- L'index n'a **aucun lien de vie** avec le modèle : il est invalidé par toute
  mutation, à la charge de l'appelant (comme `TopExp::MapShapesAndAncestors`).
  Il n'est pas stocké dans `Model` (document § 2.3, § 5.2).

### 3.3 Requêtes de shell : un seul calcul d'usages

`shell_edge_uses` est l'unique point où l'on compose les orientations :
`orient_effectif = compose(coedge.orient, faceuse.orient)`. Les cinq requêtes
sont des filtres sur ce résultat, donc elles ne peuvent pas diverger entre
elles, et le sewing réutilisera la même fonction pour trouver les arêtes
libres et vérifier l'orientation après fusion.

- `unordered_map<EdgeId, vector<EdgeUse>>` : la plupart des shells ont
  beaucoup moins d'arêtes que le modèle n'a d'arêtes au total (plusieurs
  solides dans un même `Model`), d'où une map plutôt qu'un vecteur dense ici.
- Les arêtes **dégénérées** sont exclues du comptage : un pôle de sphère ou un
  apex de cône n'est utilisé qu'une fois et ne doit pas rendre le shell ouvert
  (document § 3.4).
- `is_orientable` vérifie la **cohérence actuelle** : un shell dont une face
  est retournée échoue (le sewing la réparera en basculant le `FaceUse`), un
  ruban de Möbius échoue aussi (irréparable). Le nom suit le document ; voir
  point à valider 2.
- `free_edges` et `non_manifold_edges` sont triées par indice : sortie
  déterministe pour les rapports et les tests.

### 3.4 Boîte englobante conservative

`BoundingBox<T>` est une valeur avec `add(point, tol)`, `add(box)`, `inflate`,
`contains`, `intersects`, `diagonal` ; vide par défaut (`min > max`).
`bounding_box` échantillonne : sommets gonflés de `tol`, `n` points par arête
non dégénérée gonflés de `tol`, et une grille `n × n` de la surface pour les
faces à `natural_bounds`. Limite assumée et documentée dans le code : une face
**découpée** qui bombe au-delà de ses arêtes peut dépasser la boîte. Pour le
filtre grossier du sewing (paires d'arêtes candidates) cela suffit ; une
boîte exacte de face viendra avec le maillage ou l'intégration sur le
contour (PR 8).

## 4. Écarts par rapport au document de conception

| Document (§ 5) | Implémentation | Raison |
|---|---|---|
| `wires_of(EdgeId) -> (WireId, index)` | `coedges_of(EdgeId) -> std::span<const CoEdgeRef>` | nom plus juste : on désigne une co-arête ; même information |
| — | `face_of(WireId)` | nécessaire pour remonter d'une co-arête à sa face |
| — | `is_chained(wire)` | distinguer « ouvert mais chaîné » de « cassé » pour `make_wire` (PR 3) |
| — | `shell_edge_uses`, `EdgeUse`, `non_manifold_edges` | briques exposées pour le sewing et `check()` |
| `bounding_box -> pair<point, point>` | type `BoundingBox<T>` | il faut `intersects`/`inflate` pour le filtre du sewing |
| `signed_volume` | non inclus | prévu en PR 8 avec `make_solid`, comme dans le plan |

## 5. Invariants et complexité

- Toutes les fonctions sont `const` sur le modèle, ne lèvent que via les
  accesseurs typés (handle mort dans un `Compound`, par exemple) ou
  `coedge_point` sur un indice hors du wire (`std::vector::at`).
- `explore` : O(entités atteintes) ; `TopologyIndex` : O(entités sous la
  racine) à la construction, O(1) par requête ; `shell_edge_uses` : O(co-arêtes
  du shell) ; `bounding_box` : O(n · arêtes + n² · faces naturelles) évaluations.

## 6. Tests (`tests_brep_explore`, 219 assertions)

| Test | Couvre |
|---|---|
| `explore_by_type` | comptes depuis solide / face / arête / sommet ; pas de remontée ; ordre de découverte ; handle invalide ⇒ vide ; face morte ignorée mais sommets encore atteints par les autres faces |
| `compound_and_duplicates` | compound mélangeant solide, face, arête, sommet : aucun doublon ; compound imbriqué |
| `topology_index` | 2 faces et 2 co-arêtes par arête, 3 arêtes par sommet, 1 shell par face, `face_of` ; spans vides hors table ; index restreint à une face |
| `wire_closure` | wires de la boîte fermés ; wire tronqué chaîné mais ouvert ; wire permuté cassé ; co-arête retournée casse la chaîne ; wire vide ; `coedge_point` : fin de i = début de i+1 |
| `shell_queries` | boîte fermée/variété/orientée ; boîte ouverte : 4 arêtes libres = arêtes de la face retirée ; face retournée ⇒ non orientable ; shell entièrement retourné ⇒ orientable ; 7ᵉ face ⇒ non variété, 4 arêtes à 3 usages ; arête dégénérée ignorée ; shell vide non fermé |
| `bounding_box` | boîte `[−tol, 1+tol]³`, diagonale, `contains` ; arête ; sommet ; `intersects`/`inflate`/`add` ; nappe quadratique bombée (z atteint 1 au centre) ; face découpée sans arêtes ⇒ boîte vide |

La boîte de référence est déplacée dans `tests/tests_brep_helpers.h`
(namespace `brep_tests`) et partagée avec `tests_brep_model`.

## 7. Points à valider

1. **Pas de remontée dans `explore`** (`explore<ShellId>(m, face)` vide) : la
   remontée passe par `TopologyIndex`. Alternative : un `explore_up`.
2. **Sémantique de `is_orientable`** = cohérence d'orientation actuelle (une
   face retournée ⇒ `false`). Le nom suit le document ; `is_consistently_oriented`
   serait plus exact mais plus long. Garder le nom ?
3. **Arêtes dégénérées exclues des usages** : un shell « sphère » (seam + deux
   dégénérées) sera fermé avec une seule face. Confirmer que c'est le
   comportement voulu pour `check()` et IGES.
4. **`BoundingBox` conservative mais non garantie pour les faces découpées
   bombées** : acceptable pour le filtre du sewing ; la boîte exacte viendra
   plus tard.
5. **`TopologyIndex` sans invalidation automatique** : responsabilité de
   l'appelant, comme OCCT. Alternative : un compteur de version dans `Model`
   et une assertion dans l'index.
6. **`CoEdgeRef` comme localisation de co-arête** (pas de `CoEdgeId` dans le
   modèle) : conforme au document § 2.5 ; à confirmer avant que le sewing s'en
   serve.
