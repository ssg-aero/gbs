# BREP palier 1 — PR 3 : builders de sommets, arêtes et wires (`gbs-brep/builders.h`)

| | |
|---|---|
| PR | branche `feat/brep-builders-wire` |
| Document de référence | [brep_core.md](brep_core.md), section 6.0 et question 13 ; notes précédentes [PR 1](brep_pr01_model.md), [PR 2](brep_pr02_explore.md) |
| Fichiers | `gbs-brep/builders.h` (455 l.), `gbs-brep/brep` (agrégateur), `tests/tests_brep_wire.cpp` (337 l.), `inc/topology/{basetopo,vertex,edge,wire}.h` (dépréciation), `tests/tests_topo.cpp` (avertissements de dépréciation neutralisés), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Fournir les premiers constructeurs : un sommet, une arête (depuis une courbe,
deux points ou deux sommets existants), une arête dégénérée, et un wire libre
obtenu en chaînant des arêtes données dans n'importe quel ordre et n'importe
quel sens, avec fusion des sommets à tolérance près. C'est le remplaçant
natif de `Wire::addEdge` de l'ancienne topologie et de
`BRepBuilderAPI_MakeEdge` / `BRepBuilderAPI_MakeWire` côté OCCT.

C'est aussi la première PR écrite pour C++23 (#115) : les builders renvoient
un `std::expected`.

## 2. Vue d'ensemble

```
BuildResult<Id> = std::expected<Id, BuildError>
BuildError { BuildErrc code; std::string message; std::vector<ShapeId> shapes; }
unwrap(result) -> Id            // lève BRepError si erreur

make_vertex(m, p, tol)                                   -> VertexId
make_edge(m, curve, tol)                                 -> EdgeId   bornes = curve->bounds()
make_edge(m, curve, u1, u2, tol)                         -> EdgeId   sommets créés (1 seul si fermée)
make_edge(m, curve, u1, u2, v1, v2, tol)                 -> EdgeId   sommets imposés, tolérances remontées
make_edge(m, p1, p2, tol)                                -> EdgeId   segment, paramétrage par longueur d'arc
make_edge(m, v1, v2, tol)                                -> EdgeId   segment entre sommets existants
make_degenerate_edge(m, v, u1 = 0, u2 = 1)               -> EdgeId   curve = nullptr, degenerate = true
make_wire(m, span<const EdgeId> | vector<EdgeId>, tol)   -> WireId   chaînage + fusion, pcurves nulles
```

## 3. Choix d'architecture et justification

### 3.1 `std::expected` plutôt qu'exceptions

Le document de conception (question 13) recommandait des exceptions pour les
préconditions et des rapports pour les diagnostics. Avec C++23,
`std::expected` couvre les deux usages d'une seule façon :

- un échec **prévisible** (arêtes non chaînables, sommet hors de la courbe)
  est une valeur que l'appelant teste, sans `try` ; c'est le cas courant
  pour un lecteur STEP ou le sewing qui essaient et passent à la suite ;
- `unwrap()` donne le comportement « exception » à qui le préfère, avec un
  `BRepError` dont le message contient le code et le détail ;
- `BuildError::shapes` désigne les entités en cause (le sommet trop loin,
  l'arête dupliquée), utile pour un rapport ou une sélection dans une IHM.

Les exceptions restent réservées aux **erreurs de programmation** : accéder
à un identifiant mort par `m.edge(id)` lève toujours `BRepError` (PR 1).

`BuildErrc` est une énumération fermée ; `to_string()` la rend lisible. Les
builders des PR suivantes l'étendront (faces, sewing).

### 3.2 Garantie forte : un échec ne modifie pas le modèle

Chaque builder valide tout **avant** la première écriture. Pour `make_wire`,
la fusion des sommets et le chaînage sont calculés sur une copie locale
(union-find sur les sommets, graphe sommet → arêtes) ; le modèle n'est
modifié que si le chaînage réussit. Les tests vérifient les comptes avant et
après chaque échec.

### 3.3 Paramètres non déduits (`std::type_identity_t`)

`T` est déduit **uniquement** du `Model<T>&`. Les autres paramètres
(`shared_ptr<Curve<T,3>>`, bornes, tolérance, points) sont en contexte non
déduit. Sans cela, `make_edge(m, circle)` avec un
`shared_ptr<BSCurveRational<double,3>>` ne compile pas (la déduction de `T`
échoue sur la conversion dérivée → base), et `make_edge(m, crv, 0, 1)` avec
des littéraux entiers non plus. L'appelant passe donc n'importe quelle courbe
dérivée de `Curve<T,3>` et n'importe quel scalaire convertible.

### 3.4 Arêtes

- **Bornes** : `u1 < u2`, comprises dans `curve->bounds()` à `knot_eps` près.
  Une courbe non bornée (`Line::bounds()` renvoie les limites du type) est
  détectée par un seuil (`|u| ≥ max/16`) et exige des bornes explicites.
- **Arête fermée** : si `‖C(u1) − C(u2)‖ ≤ tol`, un seul sommet (cercle complet).
- **Arête dégénérée refusée** dans `make_edge` : si la courbe reste à `tol` de
  son départ sur quatre échantillons, `DegenerateEdge`. Les arêtes dégénérées
  se créent explicitement par `make_degenerate_edge` (pôles, apex), que
  `make_face` (PR 4) utilisera.
- **Sommets imposés** : chaque sommet doit être à `max(tol, vertex.tol)` de
  l'extrémité de courbe ; sa tolérance est ensuite remontée pour contenir
  l'extrémité et rester `≥ edge.tol` (règle vertex ≥ edge du § 4.1).
- **Segment** : `build_segment` avec paramétrage par longueur d'arc, comme
  l'ancien `Edge(pt1, pt2)`.
- La géométrie n'est jamais modifiée ni copiée : l'arête partage le
  `shared_ptr` reçu.

### 3.5 Wire : fusion puis chaînage

1. **Validation** : liste non vide, arêtes vivantes, sans doublon, sans arête
   dégénérée (elles n'ont de sens que dans le contour d'une face, que les
   builders de faces construisent eux-mêmes).
2. **Fusion des sommets** (union-find) : deux sommets d'extrémité fusionnent si
   `‖p_a − p_b‖ ≤ max(tol, tol_a + tol_b)`, c'est-à-dire s'ils sont proches au
   sens du paramètre **ou** si leurs boules de tolérance se recouvrent déjà.
   Le survivant est celui de plus petit indice ; son point n'est **pas
   déplacé** ; sa tolérance est remontée pour contenir les extrémités réelles
   des courbes qui y aboutissent. Coût O(s²) sur les s sommets du wire,
   négligeable pour un contour ; le sewing (PR 7) utilisera un filtre par
   boîtes.
3. **Chaînage** sur les représentants : un sommet partagé par plus de deux
   arêtes ⇒ `Branching` ; un nombre d'extrémités libres différent de 0 ou 2,
   ou des arêtes non atteintes ⇒ `Disconnected`. Zéro extrémité ⇒ wire fermé.
4. **Sens** : la première arête fournie est toujours `Forward` ; pour un cycle,
   le wire commence par elle ; pour une chaîne ouverte, la chaîne est
   retournée si nécessaire. Le résultat est donc déterministe et prévisible
   pour l'appelant.
5. **Application** : toutes les arêtes **du modèle** qui pointaient vers un
   sommet absorbé sont redirigées vers le survivant (y compris celles hors du
   wire, sinon elles pointeraient vers une tombstone), puis les absorbés sont
   effacés.

Le wire produit est **libre** : `pcurve == nullptr` sur chaque co-arête.

Différence avec l'ancien `Wire::addEdge` : celui-ci ajoutait une arête à la
fois, seulement aux extrémités, et modifiait l'arête passée
(`setVertex1/2`) ; ici l'ensemble est traité d'un coup, dans n'importe quel
ordre, et seules les références aux sommets changent.

### 3.6 Dépréciation de l'ancienne topologie

`BaseTopo`, `Vertex`, `Edge`, `Wire` de `inc/topology/` portent
`[[deprecated("legacy topology, use the native BREP core gbs::brep")]]`. Leur
seul utilisateur est `tests/tests_topo.cpp`, dont les trois cas sont portés
dans `tests_brep_wire.cpp` (`vertex`, `wire_square_any_order`,
`discretize_wire`) ; il neutralise l'avertissement en attendant sa
suppression. `tests_halfedgemesh.cpp` inclut ces en-têtes sans instancier les
classes : aucun avertissement. La suppression des en-têtes et de
`tests_topo.cpp` reste prévue en PR 10.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| Exceptions pour les préconditions (Q13) | `std::expected<Id, BuildError>` + `unwrap()` | C++23 adopté (#115) ; une seule convention, garantie forte |
| `make_edge(curve, u1, u2, v1, v2)` sans tolérance | paramètre `tol` ajouté (défaut `brep_default_tolerance`) | il fixe `edge.tol` et le seuil d'acceptation des sommets |
| — | `make_edge(m, v1, v2)` | remplace l'ancien `Edge(vtx1, vtx2)` |
| Critère de fusion « à tol près » | `≤ max(tol, tol_a + tol_b)` | des sommets dont les boules se recouvrent sont déjà confondus au sens du modèle |
| Arêtes dégénérées dans `make_wire` | refusées | elles appartiennent aux contours de faces, construits par les builders de faces |

## 5. Invariants garantis

Après un succès :

- toute arête créée a `u1 < u2`, une courbe non nulle (sauf `make_degenerate_edge`)
  et des sommets dont la tolérance contient l'extrémité de courbe et est
  `≥ edge.tol` ;
- tout wire créé vérifie `is_chained` (ouvert) ou `is_closed` (fermé) de la
  PR 2, et aucune arête du modèle ne pointe vers un sommet effacé.

Après un échec : le modèle est inchangé.

## 6. Tests (`tests_brep_wire`, 10 cas, 144 assertions)

| Test | Couvre |
|---|---|
| `vertex` | tolérance explicite et par défaut ; tolérance nulle, point NaN ; `unwrap` lève |
| `edge_from_curve` | segment (bornes `[0, 5]`, sommets) ; cercle fermé ⇒ 1 sommet ; demi-cercle ; `Line` non bornée refusée puis acceptée avec bornes ; courbe nulle, bornes inversées ou hors domaine, tolérance nulle, segment trop court, courbe réduite à un point ; modèle inchangé après les échecs |
| `edge_on_vertices` | segment entre sommets ; extrémité à 0,5 tol acceptée et tolérance remontée ; tolérance de sommet plus large acceptée ; extrémité trop loin ⇒ `VertexOffCurve` avec le sommet en cause ; même sommet ; sommet mort ; arête fermée sur un sommet donné |
| `degenerate_edge` | attributs, point, bornes invalides, identifiant invalide |
| `wire_square_any_order` | port de `tests_topo.edge_wire` : 4 segments indépendants mélangés et de sens mixtes ⇒ wire fermé, 4 sommets fusionnés, première arête en tête et `Forward`, sens des autres |
| `wire_open` | chaîne ouverte donnée par le milieu, réorientée ; ordre et premier sommet ; wire d'une seule arête |
| `wire_closed_single_edge` | cercle ⇒ wire fermé d'une co-arête |
| `wire_vertex_merge` | écart 0,4 tol : fusion, arête hors wire redirigée, survivant non déplacé, toutes les extrémités dans leur boule ; écart 1e-4 refusé à la tolérance par défaut (modèle inchangé) puis accepté à 1e-3 |
| `wire_errors` | branchement en étoile (modèle inchangé), liste vide, doublon, identifiant invalide, tolérance nulle, arête dégénérée, deux cercles disjoints, segment isolé |
| `discretize_wire` | port de `tests_topo.mesh_wire` : 40 points de bord à pas 0,1 sur le carré, en suivant le sens des co-arêtes |

`tests_topo`, `tests_halfedgemesh`, `tests_brep_model` et `tests_brep_explore`
passent toujours.

## 7. Points à valider

1. **`std::expected` pour tous les builders**, exceptions réservées aux
   identifiants morts passés aux accesseurs. Cela remplace la recommandation
   de la question 13 du document de conception.
2. **Garantie forte** (échec ⇒ modèle inchangé) comme règle pour tous les
   builders à venir, y compris le sewing.
3. **Critère de fusion** `max(tol, tol_a + tol_b)` et **survivant non déplacé**
   (plus petit indice) ; alternative OCCT : placer le survivant au barycentre
   et prendre une tolérance englobante.
4. **Redirection de toutes les arêtes du modèle** vers le survivant, y compris
   celles hors du wire (coût O(arêtes du modèle) à chaque fusion).
5. **Arêtes dégénérées interdites** dans un wire libre.
6. **Première arête en tête et `Forward`** comme règle de sens du wire.
7. **Suppression de `tests_topo.cpp`** : elle a été bloquée par le garde-fou de
   permissions de la session ; le fichier est conservé avec l'avertissement
   neutralisé. Le supprimer dans cette PR (son contenu est porté) ou
   attendre la PR 10 ?
