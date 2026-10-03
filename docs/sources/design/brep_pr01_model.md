# BREP palier 1 — PR 1 : modèle de données (`gbs-brep/model.h`)

| | |
|---|---|
| PR | ssg-aero/gbs#113, branche `feat/brep-model` |
| Document de référence | [brep_core.md](brep_core.md), sections 2 (modèle), 3 (géométrie), 4 (tolérances) |
| Fichiers | `gbs-brep/model.h` (591 l.), `gbs-brep/brep` (agrégateur), `gbs-brep/CMakeLists.txt`, `gbs/gbsconstants.h` (+2 constantes), `gbs/exceptions.h` (`BRepError`), `tests/tests_brep_model.cpp`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Poser le socle de données du noyau BREP : les identifiants, les sept entités
topologiques et les deux usages (`CoEdge`, `FaceUse`), l'arène `Model<T>` qui
les possède, et les opérations de cycle de vie (`add`, `erase`, `compact`,
`append`). Aucun algorithme : pas de builder, pas de validation, pas de
parcours. Tout ce qui est construit l'est « à la main » dans les tests, comme
le ferait un lecteur STEP.

## 2. Vue d'ensemble

```
Model<T>  (arène : une table par type, indices denses)
│
├─ Table<Vertex<T>>    pnt, tol
├─ Table<Edge<T>>      curve (shared_ptr<Curve<T,3>>), u1 < u2, v1, v2, tol, degenerate, same_parameter
├─ Table<Wire<T>>      coedges[ ] ─► CoEdge { edge, orient, pcurve (shared_ptr<Curve<T,2>>) }, closed
├─ Table<Face<T>>      surface (shared_ptr<Surface<T,3>>), wires[0]=extérieur, wires[1..]=trous, tol, natural_bounds
├─ Table<Shell>        faces[ ] ─► FaceUse { face, orient }, closed
├─ Table<Solid>        outer, voids[ ]
└─ Table<Compound>     shapes[ ] : ShapeId (variant des 7 identifiants)

Descente seulement : Compound → Solid → Shell → Face → Wire → Edge → Vertex.
Aucun pointeur arrière ; l'adjacence remontante est un index calculé (PR 2).
Géométrie hors arène, partagée par shared_ptr avec le reste de gbs.
```

## 3. Choix d'architecture et justification

### 3.1 Identifiants : un seul template `Id<ShapeType>`

Le nom `Id` remplace `Handle`, utilisé dans une première version de cette PR :
dans OpenCascade, `Handle(Geom_Curve)` est un pointeur intelligent intrusif qui
**possède** l'objet et compte ses références. Nos identifiants font l'inverse :
un entier sans propriété, valide tant que le `Model` vit. Le mot `Id` dit ce
qu'il est, et s'accorde avec les alias `VertexId`, `EdgeId`, `ShapeId` et avec
`IdRemap`.

```cpp
enum class ShapeType : uint8_t { Vertex, Edge, Wire, Face, Shell, Solid, Compound };
template <ShapeType K> struct Id { static constexpr ShapeType type = K; uint32_t index{npos}; bool valid() const; <=> };
using VertexId = Id<ShapeType::Vertex>;  // … EdgeId, WireId, FaceId, ShellId, SolidId, CompoundId
using ShapeId  = std::variant<VertexId, …, CompoundId>;
```

- Le document de conception écrivait sept `struct` distinctes ; l'implémentation
  les dérive d'un seul template. Même sémantique (types forts, non
  convertibles entre eux), mais le type porte son `ShapeType`, ce qui permet
  d'écrire une fois `Model::alive<K>`, `count<K>`, `ids<Id>`, `erase<K>` et le
  `std::hash`, et au `Visitor` de la PR 2 de brancher sur `Sub::type`.
- `uint32_t` + sentinelle `npos` : 4 octets, un identifiant invalide par défaut, 4
  milliards d'entités par type, largement au-delà du besoin.
- Comparaison trois voies par défaut : les identifiants se trient, ce qui donne des
  sorties déterministes (`free_edges` triées en PR 2) et des clés de `std::map`.
- `std::hash<Id<K>>` mélange l'indice et le type : deux identifiants de types
  différents et de même indice ne collisionnent pas systématiquement si on les
  met dans une table de `ShapeId`.
- `ShapeId` est un `std::variant` : conversion implicite depuis tout identifiant
  (`explore<EdgeId>(m, FaceId{3})` compile), `shape_type()`, `shape_index()`,
  `valid()` par `std::visit`. Le coût (9 octets + discriminant) n'est payé que
  dans `Compound::shapes` et les signatures génériques.

### 3.2 Arène : `Table<E>` privée à `Model`

```cpp
template <typename E> struct Table { std::vector<E> items; std::vector<uint8_t> alive; size_t n_alive; add / is_alive / erase / compact };
```

- `items` et `alive` sont deux vecteurs parallèles : l'entité reste un agrégat
  pur (pas de champ `alive` dedans), et le test de vie est une lecture d'octet.
  `std::vector<uint8_t>` plutôt que `vector<bool>` pour avoir des références
  et un stockage prévisible.
- `n_alive` est maintenu incrémentalement : `count()` est O(1).
- `erase` = tombstone, jamais de déplacement : aucun identifiant n'est invalidé par
  une suppression. C'est ce qui rend les mutations du sewing (PR 7) sûres :
  on fusionne, on redirige, on marque mort, et on compacte à la fin si on
  veut.
- `compact` renvoie `old → new` (`npos` pour les morts), déplace les survivants
  (`std::move`), puis `Model::apply_remap` réécrit **tous** les identifiants stockés
  dans toutes les tables. C'est la seule opération qui invalide des identifiants,
  et elle rend l'`IdRemap` pour que l'appelant traduise les siens.
- `table<K>()` est un `if constexpr` sur `K` : dispatch à la compilation,
  pas de tableau de `std::any` ni de hiérarchie virtuelle d'entités.

### 3.3 Entités : agrégats, comportement en fonctions libres

Conforme au document (§ 2.5) et à la convention de la géométrie gbs. Points
d'implémentation :

- Les tolérances ont une valeur par défaut `brep_default_tolerance<T>` dans
  l'agrégat, donc `Vertex<T>{pnt}` est valide sans préciser `tol`.
- `Edge::curve` peut être `nullptr` : c'est le marqueur d'une arête dégénérée
  (avec `degenerate = true`). `edge_point()` renvoie alors le point du sommet
  au lieu de déréférencer.
- `CoEdge::pcurve` peut être `nullptr` (wire libre) ; `coedge_uv()` lève
  `BRepError` dans ce cas plutôt que de renvoyer une valeur fausse.
- `Shell`, `Solid`, `Compound`, `FaceUse` ne sont pas templates : ils ne
  contiennent que des identifiants.

### 3.4 Convention d'orientation appliquée dans les helpers

Les quatre helpers de `model.h` fixent la convention que les PR suivantes
réutilisent :

| Helper | Convention |
|---|---|
| `coedge_start/coedge_end(m, ce)` | `Forward` : `v1 → v2` ; `Reversed` : `v2 → v1` |
| `edge_point(m, e, u)` | évalue `curve(u)`, `u` dans le repère de l'arête |
| `coedge_uv(m, ce, t)` | `t` croît **le long de la co-arête** ; pour `Reversed`, la pcurve (paramétrée comme l'arête, § 3.2 du document) est lue en `u1 + u2 − t` |

Ainsi une pcurve est stockée une fois, dans le sens de l'arête, et le
renversement est purement arithmétique à la lecture.

### 3.5 `append` : copie topologique, géométrie partagée

`append(other)` copie uniquement les entités vivantes d'`other`, réécrit
leurs identifiants avec l'`IdRemap` renvoyé, et **partage** les `shared_ptr` de
courbes et surfaces (pas de copie profonde). C'est le comportement voulu
pour assembler des morceaux construits séparément ; une copie profonde de la
géométrie serait un `clone()` à ajouter si le besoin apparaît (il n'y en a pas
dans le palier 1).

### 3.6 Erreurs

`BRepError : std::runtime_error` dans `gbs/exceptions.h`, préfixe `"brep: "`.
Les accesseurs typés (`vertex(id)`, `edge(id)`…) lèvent sur un identifiant invalide
**ou mort**, avec le type et l'indice dans le message. Les méthodes `alive()`
permettent de tester sans lever.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| Sept `struct` de identifiants | `Id<ShapeType>` + alias | factorisation (§ 3.1) |
| `Model::add(Vertex)` … | identique ; ajout de `capacity<K>()`, `ids<Id>()`, `empty()`, `alive(ShapeId)`, `erase(ShapeId)`, `count(ShapeType)` | nécessaires aux tests et à l'explorateur |
| `coedge_point` prévu dans l'explorateur | `coedge_uv` et `edge_point` dans `model.h`, `coedge_point` (3D, par wire et indice) dans `explore.h` | les helpers sur la convention d'orientation vont avec le modèle |
| Tables de noms et d'identifiants externes (§ 9.2) | non incluses | prévues avec la lecture STEP ; `compact`/`append` devront alors les remapper aussi |
| `Model::append` | ajouté tel que décrit (§ 2.3), copie des vivants seulement | — |

## 5. Invariants : garantis ou non à ce stade

Cette PR **ne vérifie aucun invariant topologique** (c'est le rôle de `check()`
en PR 6 et des builders). Elle garantit seulement :

- un identifiant renvoyé par `add` est valide et vivant jusqu'à `erase` ;
- `erase` n'invalide aucun autre identifiant ;
- après `compact`, tous les identifiants **stockés dans le modèle** sont valides et
  désignent les mêmes entités qu'avant ; les identifiants externes se traduisent par
  l'`IdRemap` ;
- `append` ne modifie pas les entités déjà présentes.

Le modèle accepte donc des états incohérents (arête dont un sommet est mort,
wire non fermé) : c'est voulu pendant une mutation ; `check()` les détectera.

## 6. Tests (`tests_brep_model`, 376 assertions)

| Test | Couvre |
|---|---|
| `ids` | défaut invalide, égalité, ordre, `unordered_map` à clés identifiants, `ShapeId`, `reverse`/`compose` |
| `box_counts_and_access` | boîte unitaire construite à la main (8/12/6/1/1) ; chaque arête a deux usages de sens opposés ; wires chaînés ; `coedge_uv` ∘ surface = `edge_point` à `tol` près sur 4 paramètres par co-arête ; extrémités sur les sommets ; accès à un identifiant invalide lève |
| `erase_and_compact` | tombstones (capacité inchangée, compte décrémenté, accès lève) ; `compact` renumérote, remap correct pour vivants et morts, identifiants internes réécrits, géométrie des arêtes inchangée ; `compact` idempotent |
| `append` | boîte translatée ajoutée : comptes doublés, identifiants d'origine intacts, identifiants copiés réécrits, tombstone non copié, surfaces partagées |

La boîte de test oriente ses six faces normale sortante (choix des axes
`(u, v)` par face), ce qui en fait aussi le cas de référence « shell fermé,
variété, orienté » des PR suivantes.

## 7. Points à valider

1. **`Id<ShapeType>` template** plutôt que sept structures écrites à la
   main : même usage, moins de code ; les messages d'erreur du compilateur
   affichent `Id<ShapeType::Edge>` au lieu de `EdgeId`. Le nom `Id` (et non
   `Handle`, qui désigne un pointeur possédant chez OCCT) est validé.
2. **`erase` sans cascade** : supprimer une face ne supprime ni ses wires ni ses
   arêtes. Les builders feront le ménage ; faut-il offrir un
   `erase_recursive(ShapeId)` de confort dès maintenant ?
3. **`append` partage la géométrie** (pas de copie profonde). OK pour le
   palier 1 ?
4. **Tolérance par défaut dans l'agrégat** (`Vertex<T>{pnt}` ⇒ `tol = 1e-6`)
   plutôt qu'une construction obligatoirement explicite.
5. **`Edge::same_parameter` par défaut `true`** : une arête créée à la main est
   supposée exacte tant qu'un builder n'a pas mesuré d'écart.
6. **Helpers d'orientation dans `model.h`** (`coedge_uv` lit `u1 + u2 − t` pour
   `Reversed`) : c'est la convention que toutes les PR suivantes appliqueront ;
   à valider maintenant car elle est coûteuse à changer après.
