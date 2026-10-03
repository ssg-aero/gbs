# Noyau BREP natif — conception du palier 1

| | |
|---|---|
| Statut | Validé le 3 octobre 2026 : les quinze questions ouvertes de la section 11 sont tranchées selon les recommandations |
| Périmètre | Palier 1 : modèle de données BREP et constructeurs de base |
| Dépend de | `gbs/` (courbes, surfaces, extrema), `gbs-io/iges.h` (libIGES) |
| Prépare | Palier 2 (intersections, découpe, sweep, offsets), lecture STEP AP203/AP214 |
| Hors périmètre définitif | Booléens généraux, healing, congés : restent délégués à OCCT via `gbs-occt` |

Ce document décrit *quoi* construire et *pourquoi*. Il ne contient pas de code
exécutable : les pseudo-déclarations C++ sont illustratives et pourront bouger
à l'implémentation tant que les décisions de fond (section 0) sont respectées.

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Existant `inc/topology/{basetopo,vertex,edge,wire}.h` | **Remplacer** par un nouveau noyau dans `gbs-brep/`, namespace `gbs::brep` ; les anciens en-têtes sont retirés en fin de palier après migration des tests | Étendre l'embryon actuel |
| Propriété / identité des entités | **Arène** `Model<T>` + identifiants typés (`VertexId`, `EdgeId`, …) ; la géométrie reste en `shared_ptr<Curve>` / `shared_ptr<Surface>` | Graphe de `shared_ptr` + `weak_ptr` |
| Orientation | **Par usage** : `CoEdge` (arête + sens + pcurve) dans un `Wire`, `FaceUse` (face + sens) dans un `Shell` | Orientation portée par la valeur du shape (OCCT `TopoDS_Shape`) |
| Dimension | **3D seulement** : `template <std::floating_point T>`, pas de paramètre `dim` | `<T, dim>` comme la géométrie |
| Pcurves | Portées par la `CoEdge`, une par couple (arête, face), deux pour une arête seam | Liste de représentations sur l'arête (OCCT `BRep_TEdge`) |
| Surfaces élémentaires | La topologie référence `Surface<T,3>` polymorphe ; les élémentaires STEP seront converties en NURBS rationnelles exactes (décision de pilotage) ; aucune classe analytique de surface n'est ajoutée au palier 1 | Classes `Plane`, `Cylinder`, … |
| Tolérances | Une tolérance géométrique par vertex / edge / face, règle `vertex ≥ edge ≥ face` ; la « tolérance d'approximation » de `BaseTopo` devient un paramètre des builders | Deux tolérances sur chaque entité |
| Adjacence | Calculée à la demande par un index (`TopologyIndex`), pas de pointeurs arrière stockés | Pointeurs arrière dans les entités |
| Export IGES | Faces en 144/142/102/126/128 via l'API `DLL_IGES` déjà utilisée ; 186/514/510/508/504/502 reportés (API cœur de libIGES seulement) | Écrire 186 au palier 1 |
| Python | Sous-module `gbs.brep`, `T = double` | Classes à plat dans `gbs` |

---

## 1. Objectifs et non-objectifs

### 1.1 Objectifs du palier 1

1. Un **modèle de données BREP** complet en entités : `Vertex`, `Edge`, `Wire`,
   `Face`, `Shell`, `Solid`, `Compound`, avec orientation, pcurves et tolérances.
2. Un **explorateur topologique** : parcours des sous-entités par type, adjacence
   (arêtes d'une face, faces d'une arête), fermeture, variété (manifold),
   bornes, validité minimale (« BRepCheck-lite »).
3. Des **builders** : arête depuis une courbe, wire depuis des arêtes (fusion des
   sommets à tolérance près), face depuis une surface entière, face depuis un
   wire posé sur une surface (pcurves extraites ou projetées), shell par
   **sewing simple** de faces partageant des arêtes à tolérance près, solide
   depuis un shell fermé, compound.
4. Un **export IGES** des faces (surface découpée 144), des shells et solides
   (ensemble de 144 nommées) avec libIGES déjà en dépendance.
5. Une **API Python** cohérente avec les bindings existants.
6. Reproduire en natif ce que `gbs-occt/brepbuilders.h` délègue aujourd'hui à
   OCCT : `to_edge`, `to_wire`, `to_face`, `to_shell` (sewing), `explode`.

### 1.2 Ce que le palier 1 prépare sans l'implémenter

- Palier 2 : intersections courbe-surface et surface-surface, découpe de faces
  (faces à plusieurs loops, trous), sweep et révolution en BREP, offsets de faces.
- Lecture STEP AP203/AP214 : `advanced_face`, `edge_loop`, `oriented_edge`,
  `closed_shell`, `manifold_solid_brep`, `brep_with_voids`, identifiants, noms,
  unités.

### 1.3 Non-objectifs

- Booléens généraux, healing (recollement avec découpe d'arêtes, jonctions en T),
  congés : restent délégués à OCCT via `gbs-occt`.
- Historique de construction, attributs graphiques, assemblages avec placements.
- Maillage des faces BREP (le mailleur half-edge de `inc/topology/halfEdgeMesh*.h`
  reste indépendant ; voir § 2.2).
- Compatibilité binaire avec `TopoDS_Shape` : aucune conversion native ↔ OCCT au
  palier 1 (candidate évidente pour plus tard, côté `gbs-occt`).

---

## 2. Modèle de données

### 2.1 Bilan de l'embryon actuel et décision

L'existant (`inc/topology/basetopo.h`, `vertex.h`, `edge.h`, `wire.h`) a été
écrit pour alimenter le mailleur, pas pour décrire un solide :

| Constat | Conséquence pour un BREP |
|---|---|
| `BaseTopo` impose `virtual void tessellate() = 0` | La topologie dépend du maillage ; une face BREP n'a pas à savoir se tesseler |
| `Vertex` stocke un `shared_ptr<HalfEdgeVertex>` | Couplage direct au half-edge mesh ; un sommet BREP ne porte qu'un point et une tolérance |
| `Edge` stocke `he_vertices` (sa discrétisation) et possède ses deux `Vertex` par `shared_ptr` | Pas de partage contrôlé d'une arête entre deux faces ; la discrétisation n'a rien à faire dans la topologie |
| `Wire::addEdge` fusionne les sommets en **mutant** l'arête passée (`setVertex1/2`) | Effet de bord sur une arête potentiellement partagée par une autre face |
| Pas d'orientation, pas de pcurve, pas de face | Tout l'essentiel manque |
| Deux tolérances (`precision`, `approximation`) sur chaque entité | Sémantique floue : l'approximation est un paramètre d'algorithme, pas une propriété d'un sommet |
| `Edge`, `Wire` ne sont référencés nulle part hors `tests/tests_topo.cpp` | Le remplacement n'a pas d'impact utilisateur |

**Décision : remplacer.** Le nouveau noyau vit dans un module header-only
`gbs-brep/` (au même rang que `gbs-io/`, `gbs-mesh/`, installé par
`INSTALL_HEADERS`), namespace `gbs::brep`. Le namespace évite la collision des
noms `Vertex`/`Edge`/`Wire` pendant la transition et sépare nettement le BREP
des structures `HalfEdge*` qui restent dans `gbs`. Les trois tests de
`tests/tests_topo.cpp` sont réécrits sur la nouvelle API (PR 3) ; les anciens
en-têtes sont supprimés dans la dernière PR du palier. Ce qui est **repris** de
l'embryon : l'idée d'une tolérance par entité, la fusion de sommets à tolérance
près dans le wire, et la mesure de coïncidence point–courbe par
`extrema_curve_point`.

### 2.2 BREP et half-edge mesh : deux structures, deux rôles

| | `gbs::brep` (ce document) | `gbs::HalfEdge*` (`inc/topology/halfEdgeMesh*.h`, `gbs-mesh/`) |
|---|---|---|
| Rôle | Frontière exacte d'un solide : surfaces et courbes paramétriques | Maillage polygonal : triangles/quads, lissage, qualité, Delaunay |
| Primitive | Face = surface + loops de co-arêtes avec pcurves | Face = cycle de demi-arêtes ; sommet = coordonnées |
| Orientation | Par usage (`CoEdge`, `FaceUse`) | Par la demi-arête elle-même |
| Lien entre les deux | Palier ultérieur : un mailleur de face BREP produira un `HalfEdgeMesh` par face à partir des pcurves (le mailleur 2D Delaunay existe déjà) | — |

Le mot *half-edge* est donc réservé au maillage. La `CoEdge` BREP joue un rôle
analogue (une arête vue depuis une face) mais n'a ni `next`/`previous`/`opposite`
ni coordonnées propres : elle n'est qu'un **usage** d'une `Edge` dans un `Wire`.

### 2.3 Propriété et identité : arène et identifiants

Trois modèles existent :

| | OCCT | Parasolid | Graphe `shared_ptr` (embryon actuel) |
|---|---|---|---|
| Identité | `TopoDS_Shape` = `Handle(TopoDS_TShape)` + `TopLoc_Location` + `TopAbs_Orientation` ; identité par le `TShape` partagé | *Tag* entier dans une *partition* ; les entités sont des enregistrements du noyau | Adresse du `shared_ptr` |
| Partage | Un `TShape` référencé par plusieurs parents ; descente seulement (parent → enfants) | Pointeurs dans les deux sens (edge → fins, fin → loop → face) | Descente seulement, cycles impossibles sans `weak_ptr` |
| Adjacence | Index calculé à la demande (`TopExp::MapShapesAndAncestors`) | Directe par les pointeurs | À calculer |
| Sérialisation / STEP | Via maps shape → id | Tag ≈ id | Via maps |

**Décision : arène + identifiants typés.**

*Arène* (de l'anglais *arena allocator*) : un conteneur unique, le `Model`, qui
**possède toutes les entités** et les range dans des tableaux contigus, un par
type. Une entité n'est pas un objet alloué individuellement et pointé par un
`shared_ptr` ; c'est une case d'un `std::vector`, désignée par son indice
enveloppé dans un type fort (l'*identifiant typé*, par exemple `EdgeId{12}`). Les
relations entre entités sont des identifiants, jamais des pointeurs. La durée de vie
des entités est celle du `Model` ; la géométrie, elle, reste hors de l'arène
en `shared_ptr` partagé avec le reste de gbs.

```cpp
namespace gbs::brep {

    // Identifiants : entiers forts, invalides par défaut. Aucune arithmétique.
    struct VertexId   { std::uint32_t index{npos}; };
    struct EdgeId     { std::uint32_t index{npos}; };
    struct WireId     { std::uint32_t index{npos}; };
    struct FaceId     { std::uint32_t index{npos}; };
    struct ShellId    { std::uint32_t index{npos}; };
    struct SolidId    { std::uint32_t index{npos}; };
    struct CompoundId { std::uint32_t index{npos}; };

    // Identifiant générique : ce qu'un Compound contient et ce que l'explorateur renvoie.
    using ShapeId = std::variant<VertexId, EdgeId, WireId, FaceId, ShellId, SolidId, CompoundId>;

    template <std::floating_point T>
    class Model
    {
        std::vector<Vertex<T>>   m_vertices;
        std::vector<Edge<T>>     m_edges;
        std::vector<Wire<T>>     m_wires;
        std::vector<Face<T>>     m_faces;
        std::vector<Shell>       m_shells;
        std::vector<Solid>       m_solids;
        std::vector<Compound>    m_compounds;
        // attributs optionnels (nom, id externe) : tables annexes, voir § 9
    public:
        auto vertex(VertexId) const -> const Vertex<T> &;
        auto vertex(VertexId)       ->       Vertex<T> &;
        // idem edge(), wire(), face(), shell(), solid(), compound()
        auto add(Vertex<T>) -> VertexId; // idem pour chaque type
        auto erase(ShapeId) -> void;     // tombstone, voir ci-dessous
        auto compact() -> IdRemap;       // supprime les tombstones, renvoie l'ancienne → nouvelle numérotation
    };
}
```

Pourquoi l'arène :

- **Mutations topologiques sûres.** Le sewing (palier 1) et la découpe (palier 2)
  remplacent une arête par une autre dans toutes les faces qui l'utilisent. Avec
  des identifiants, c'est une substitution d'entiers dans les `CoEdge` concernées, sans
  risque de cycle ni de pointeur pendant.
- **Identité stable et sérialisable.** Les ids STEP (`#123`) et les pointeurs DE
  d'IGES se mappent sur des indices ; Python manipule des entiers et non des
  adresses.
- **Copie et comparaison triviales.** Un `FaceId` se copie, se hache, se trie ; on
  peut mettre des identifiants dans des `std::unordered_map` sans écrire de hash
  pour des pointeurs.
- **Cohérence avec la géométrie.** Les courbes et surfaces restent des
  `std::shared_ptr<Curve<T,3>>` / `std::shared_ptr<Surface<T,3>>`, partagées
  comme partout dans gbs (`CurveOnSurface`, `CurveTrimmed`, `SurfaceOfRevolution`
  font déjà ainsi). Une même surface peut être référencée par plusieurs faces
  (après découpe au palier 2).

Coûts assumés :

- Une entité n'existe que dans son `Model`. Pas d'« arête libre » hors modèle :
  on crée toujours un `Model` d'abord. Pour l'usage courant (une pièce = un
  modèle) ce n'est pas une gêne ; pour combiner deux modèles on prévoit
  `Model::append(const Model&) -> IdRemap`.
- Suppression par *tombstone* (drapeau `alive = false`) puis `compact()` explicite
  qui renumérote ; les identifiants détenus par l'utilisateur sont invalidés par
  `compact()`, jamais par `erase()`. Le palier 1 n'appelle `erase()` que dans le
  sewing (arêtes et sommets fusionnés).

Comparaison : c'est le modèle Parasolid (partition + tags) plutôt que celui
d'OCCT (`Handle(...)` à comptage de références). On ne reprend pas de Parasolid les pointeurs arrière
stockés : l'adjacence est un index reconstruit à la demande (§ 5), comme
`TopExp::MapShapesAndAncestors`, ce qui garde le modèle acyclique et les
mutations locales.

### 2.4 Orientation : par usage, avec des co-arêtes

| | OCCT | Parasolid | Proposition |
|---|---|---|---|
| Où est l'orientation d'une arête dans une face ? | Dans la valeur `TopoDS_Edge` listée par le wire (`FORWARD`/`REVERSED`, plus `INTERNAL`/`EXTERNAL`) | Dans le *fin* (une entité par couple loop–edge) | Dans la `CoEdge` (un enregistrement par couple wire–edge) |
| Où est la pcurve ? | Sur le `BRep_TEdge`, dans une liste de représentations indexée par (surface, location) ; cas seam : `BRep_CurveOnClosedSurface` avec deux pcurves, choisie selon l'orientation de l'arête | Sur le fin (SP-curve) | Sur la `CoEdge` |
| Orientation d'une face dans un shell | Dans la valeur `TopoDS_Face` listée par le shell | Région / sens du shell | Dans le `FaceUse` du shell |

**Décision : orientation par usage via `CoEdge` et `FaceUse`.** Une `CoEdge` est
exactement un *fin* Parasolid, un *coedge* ACIS, ou un `oriented_edge` STEP
auquel on attache la pcurve. Justification :

- La pcurve **appartient au couple (arête, face)** : la mettre sur la `CoEdge`
  rend ce fait structurel. OCCT le contourne en indexant par surface, ce qui
  oblige `BRep_Tool::CurveOnSurface(edge, face)` à une recherche, et le cas seam
  à un type spécial.
- Une **arête seam** (cylindre, révolution fermée) apparaît deux fois dans le même
  wire, avec deux `CoEdge` de sens opposés et deux pcurves différentes
  (`u = u1` et `u = u2`) : aucun cas particulier.
- Une arête d'intersection (palier 2) portée par deux faces a naturellement une
  pcurve sur chacune.
- STEP se lit et s'écrit sans traduction : `edge_loop` = `Wire`, `oriented_edge` =
  `CoEdge`, `oriented_face` / `same_sense` = `FaceUse`.

On n'introduit pas les orientations `INTERNAL`/`EXTERNAL` d'OCCT (arêtes et
sommets internes à une face) : hors besoin du palier 1 et 2 ; une `enum class
Orientation : uint8_t { Forward, Reversed }` suffit et laisse la place.

### 2.5 Entités et attributs

Toutes les entités sont des `struct` de données simples (agrégats), sans
hiérarchie virtuelle ; les comportements sont des fonctions libres prenant le
`Model`. C'est la convention de la géométrie gbs (classes de données +
fonctions libres `build_*`, `extrema_*`, `loft`…).

```cpp
namespace gbs::brep {

    enum class Orientation : std::uint8_t { Forward, Reversed };

    template <std::floating_point T>
    struct Vertex
    {
        point<T,3> pnt;
        T          tol;          // rayon de la boule de confusion, > 0
    };

    template <std::floating_point T>
    struct Edge
    {
        std::shared_ptr<Curve<T,3>> curve;   // nullptr ssi degenerate
        T        u1, u2;                     // bornes sur curve, u1 < u2 (sens de l'arête = sens de la courbe)
        VertexId v1, v2;                     // v1 = curve(u1), v2 = curve(u2) à tol près ; v1 == v2 autorisé (arête fermée ou dégénérée)
        T        tol;                        // tube de tolérance autour de la courbe
        bool     degenerate{false};          // arête réduite à un point (pôle de sphère, apex de cône)
        bool     same_parameter{true};       // les pcurves sont paramétrées comme curve (voir § 3.4)
    };

    template <std::floating_point T>
    struct CoEdge
    {
        EdgeId      edge;
        Orientation orient;                  // Forward : parcours de v1 vers v2
        std::shared_ptr<Curve<T,2>> pcurve;  // dans (u,v) de la surface de la face ; nullptr pour un wire libre
    };

    template <std::floating_point T>
    struct Wire
    {
        std::vector<CoEdge<T>> coedges;      // ordonnées, bout à bout (fin de i == début de i+1 au sens de l'orientation)
        bool closed{false};
    };

    template <std::floating_point T>
    struct Face
    {
        std::shared_ptr<Surface<T,3>> surface;
        std::vector<WireId> wires;           // wires[0] = contour extérieur, suivants = trous ; tous fermés
        T    tol;
        bool natural_bounds{false};          // wires[0] coïncide avec le rectangle paramétrique (optimisation IGES N1 = 0)
    };

    struct FaceUse
    {
        FaceId      face;
        Orientation orient;                  // Forward : la normale de la surface pointe vers l'extérieur du shell
    };

    struct Shell
    {
        std::vector<FaceUse> faces;
        bool closed{false};                  // chaque arête utilisée exactement deux fois, en sens opposés
    };

    struct Solid
    {
        ShellId              outer;
        std::vector<ShellId> voids;          // cavités ; vide au palier 1, prévu pour STEP brep_with_voids
    };

    struct Compound
    {
        std::vector<ShapeId> shapes;
    };
}
```

Remarques :

- Il n'y a **pas d'entité `Loop` distincte** : `Wire` sert à la fois de chaîne
  libre (entrée des builders, équivalent `TopoDS_Wire`) et de contour de face
  (`edge_loop` STEP, loop Parasolid). Un wire attaché à une face doit être fermé
  et ses `CoEdge` doivent porter une pcurve ; un wire libre peut être ouvert et
  sans pcurves. Cela évite une conversion wire → loop à chaque `make_face`.
- Il n'y a **pas d'entité `Fin`/`CoEdge` adressable** par identifiant : la `CoEdge`
  est une valeur dans le wire. On y accède par `(WireId, index)`. Si le palier 2
  en a besoin (découpe modifiant une co-arête depuis une autre face), on pourra
  introduire un `CoEdgeId = {WireId, uint32_t}` sans changer le stockage.
- Pas de **placement** (`TopLoc_Location`) : une face transformée est une face
  dont la surface est transformée (`gbs/transform.h` sait le faire pour les
  NURBS). Les assemblages STEP avec placements seront traités à la lecture par
  application de la transformation à la géométrie. Question ouverte Q5.
- `Compound` peut contenir n'importe quel identifiant, y compris d'autres compounds
  (comme `TopoDS_Compound`). On n'introduit pas `CompSolid`.

### 2.6 Invariants du modèle

Ces invariants sont ceux que `check()` (§ 5.4) vérifie et que les builders
garantissent :

1. `Edge` : `u1 < u2` ; `curve` non nul sauf `degenerate` ; `dist(curve(u1), v1.pnt) ≤ v1.tol` et idem en `u2` ; `tol ≤ min(v1.tol, v2.tol)` (§ 4).
2. `CoEdge` dans un wire de face : `pcurve` non nulle, bornes de la pcurve = `[u1, u2]` de l'arête (*SameRange*), déviation `‖S(pcurve(t)) − C(t)‖ ≤ edge.tol` sur un échantillon (*SameParameter*).
3. `Wire` fermé : le sommet de fin de la co-arête *i* (selon son orientation) est le sommet de début de la co-arête *i+1*, cycliquement.
4. `Face` : `wires[0]` tourne dans le sens direct du plan (u,v) (aire signée des pcurves > 0), les trous dans le sens indirect et à l'intérieur de `wires[0]`.
5. `Shell` : chaque arête est utilisée par au plus deux `CoEdge` au travers des `FaceUse` (variété) ; `closed` ssi exactement deux, et de sens globalement opposés une fois les `FaceUse` appliqués.
6. `Solid` : `outer` fermé, orienté normales sortantes ; chaque cavité fermée, orientée normales entrantes.

---

## 3. Géométrie associée

### 3.1 Courbe 3D d'arête

- `Edge::curve` est un `shared_ptr<Curve<T,3>>` polymorphe : `BSCurve`,
  `BSCurveRational`, `Line`, `Circle3d`, `CurveOnSurface`, `CurveTrimmed`,
  `CurveReversed`, `CurveComposite`… tout ce qui dérive de `Curve<T,3>` est
  accepté. La topologie n'exige que `value(u, d)` et `bounds()`.
- Les bornes `[u1, u2]` de l'arête sont stockées dans l'arête, pas dans la
  courbe : une même courbe peut porter plusieurs arêtes (après découpe). La
  courbe peut être non bornée (`Line::bounds()` renvoie `[min, max]`), c'est
  l'arête qui borne.
- **Sens** : l'arête est toujours parcourue de `v1` à `v2` dans le sens croissant
  de `u`. On ne crée jamais d'arête « inversée » : l'inversion est le rôle de
  `CoEdge::orient`. Si l'utilisateur fournit une courbe à inverser, le builder
  l'enveloppe dans `CurveReversed` (qui existe déjà) plutôt que de la modifier.
- **Arête fermée** (cercle complet, iso d'une surface périodique) : `v1 == v2`,
  autorisé. Elle n'est pas dégénérée.

Comparaison : c'est `BRep_Tool::Curve(edge, first, last)` d'OCCT, sans la liste
de représentations alternatives (`Polygon3D`, `PolygonOnTriangulation`) qui
relèvent du maillage.

### 3.2 Pcurves

- `CoEdge::pcurve` est un `shared_ptr<Curve<T,2>>` dans l'espace `(u, v)` de
  `Face::surface`. Concrètement une `BSCurve<T,2>` (ou rationnelle) la plupart
  du temps, une `Line<T,2>` pour les iso-lignes des faces à bornes naturelles.
- **Convention de paramétrage** : la pcurve est paramétrée **comme l'arête**
  (mêmes bornes `[u1, u2]`, même sens `v1 → v2`), **quelle que soit
  l'orientation de la co-arête**. Le renversement est appliqué à la lecture
  (`CoEdge::orient == Reversed` ⇒ parcourir de `u2` vers `u1`). Une seule
  convention à retenir, et la même pcurve sert si la face change de sens dans
  un shell.
- Une pcurve vaut pour un couple (arête, face). Deux faces partageant une arête
  ont chacune leur co-arête, donc leur pcurve. Une arête seam a deux co-arêtes
  dans le même wire, donc deux pcurves.
- **Exactitude** : quand la courbe 3D de l'arête est une `CurveOnSurface` posée
  sur la même surface, la pcurve est sa `basisCurve()` et l'arête est exacte
  (`same_parameter` sans déviation). C'est le cas des faces à bornes
  naturelles, et ce sera celui des courbes d'intersection du palier 2 si on
  les représente par leurs pcurves. Sinon la pcurve est approchée (§ 6.2).

### 3.3 Surfaces

- `Face::surface` est un `shared_ptr<Surface<T,3>>` polymorphe : `BSSurface`,
  `BSSurfaceRational`, `SurfaceOfRevolution`, `SurfaceOffset`… La topologie
  n'exige que `value(u, v, du, dv)` et `bounds()`.
- **Surfaces élémentaires** (plan, cylindre, cône, sphère, tore) : conformément à
  la décision de pilotage, elles seront représentées par des **NURBS rationnelles
  exactes** à la lecture STEP. gbs sait déjà construire des cercles et ellipses
  rationnels exacts (`build_circle`, `build_ellipse`) ; les surfaces
  élémentaires s'obtiennent par extrusion/révolution de ces courbes (le loft
  rationnel existe, `SurfaceOfRevolution` existe). Le palier 1 **n'ajoute pas**
  de classes `Plane`/`Cylinder`/… : la polymorphie de `Surface<T,3>` les
  accueillerait si on les voulait plus tard pour des intersections analytiques
  (question ouverte Q6).
- **Périodicité et fermeture** : `Surface<T,3>` n'expose ni `isUClosed()` ni
  `isUPeriodic()`. Le palier 1 a besoin de savoir si `S(u1, v) == S(u2, v)` pour
  fabriquer les seams (§ 6.1). Proposition : une fonction libre
  `surface_closure(srf, tol) -> {closed_u, closed_v, degenerate_u1, degenerate_u2, degenerate_v1, degenerate_v2}`
  qui échantillonne les quatre bords et compare à `tol`. C'est robuste pour les
  NURBS et `SurfaceOfRevolution` et n'oblige pas à toucher la hiérarchie
  géométrique. On peut plus tard ajouter des méthodes virtuelles avec cette
  détection numérique comme implémentation par défaut (question ouverte Q7).

### 3.4 SameParameter, arêtes dégénérées, seams

**SameParameter.** Pour un paramètre `t ∈ [u1, u2]`, `S_f(pcurve_f(t))` et
`C(t)` ne coïncident en général qu'à `edge.tol` près (la pcurve est approchée).
Comme OCCT (`BRep_TEdge::SameParameter`), on garde un drapeau et on fait
porter l'écart par la tolérance de l'arête : les builders mesurent la
déviation sur un échantillon et font `edge.tol = max(edge.tol, dev)`. On ne
reprend pas le drapeau `SameRange` : l'égalité des bornes est un invariant
imposé (une pcurve dont les bornes diffèrent est reparamétrée affinement par
`CurveReparametrized` / `changeBounds` à la construction).

**Arêtes dégénérées.** Un bord de surface réduit à un point (pôle d'une
sphère NURBS, apex d'un cône, bord `u1` d'une révolution dont la génératrice
touche l'axe) donne une arête avec `degenerate = true`, `curve == nullptr`,
`v1 == v2`, et **une pcurve dans chaque face qui l'utilise** (une iso-ligne
dans (u,v)). C'est la convention OCCT (`BRep_Tool::Degenerated`) et STEP
(`edge_curve` dégénérée avec `vertex_point` identiques). L'explorateur et IGES
traitent ce cas : la pcurve est exportée, la courbe 3D est remplacée par le
point.

**Seams.** Une surface fermée en `u` produit une arête iso `u = u1` (≡ `u = u2`)
utilisée deux fois dans le wire de la face à bornes naturelles : une co-arête
`Forward` avec pcurve `u = u1` (parcourue dans le sens de `v` croissant) et une
co-arête `Reversed` avec pcurve `u = u2`. Le wire reste une chaîne fermée
cohérente. Pour une surface fermée en `u` et `v` (tore) il y a deux seams. Pour
une sphère NURBS, un seam et deux arêtes dégénérées. Le sewing ne doit jamais
tenter d'apparier les deux co-arêtes d'un seam entre elles (elles sont déjà
sur la même arête) ; il les reconnaît par l'identité d'`EdgeId`.

### 3.5 Récapitulatif des conventions de paramétrage

| Élément | Paramètre | Sens | Bornes |
|---|---|---|---|
| `Edge` | `u` de `curve` | `v1 → v2` = `u` croissant | `[u1, u2]`, `u1 < u2`, portées par l'arête |
| `CoEdge::pcurve` | même `u` que l'arête | même sens que l'arête | mêmes `[u1, u2]` |
| `CoEdge::orient` | — | `Reversed` ⇒ lecture de `u2` vers `u1` | — |
| `Face::surface` | `(u, v)` | contour extérieur direct dans (u,v) | `bounds()` de la surface ; la face n'en couvre qu'une partie |
| `FaceUse::orient` | — | `Reversed` ⇒ normale `∂S/∂u ∧ ∂S/∂v` rentrante | — |

---

## 4. Tolérances

### 4.1 Modèle

Une tolérance **géométrique** `tol` par `Vertex`, `Edge`, `Face`, avec la
sémantique OCCT :

- `Vertex::tol` : rayon de la boule centrée sur `pnt` contenant les extrémités
  réelles de toutes les arêtes incidentes.
- `Edge::tol` : rayon du tube autour de `curve` contenant les images
  `S_f(pcurve_f(t))` de toutes les faces incidentes.
- `Face::tol` : écart maximal entre la surface et la frontière « vraie » ; en
  pratique `max` des `tol` de ses arêtes à la construction, rarement plus.

**Règle de propagation** : pour toute incidence, `vertex.tol ≥ edge.tol ≥ face.tol`.
Elle est **imposée à la construction** (les builders remontent les tolérances)
et **vérifiée par `check()`**. C'est la règle de `BRepCheck`. Parasolid a une
notion voisine avec les *tolerant edges/vertices* ; on retient l'approche
OCCT, plus simple et suffisante pour le sewing.

Les tolérances ne se **réduisent jamais** automatiquement ; `set_tolerance()`
accepte seulement d'augmenter, sauf appel explicite `force`.

### 4.2 Lien avec `BaseTopo::precision` et `approximation`

| `BaseTopo` actuel | Devient |
|---|---|
| `precision` (défaut `1e-6`) : coïncidence de deux points, de deux sommets | `Vertex::tol` / `Edge::tol` / `Face::tol` par entité, défaut `brep_default_tolerance<T>` |
| `approximation` (défaut `1e-5`) : seuil accepté pour fusionner des sommets, projeter… | **Paramètre des builders** (`sew(…, tol)`, `make_face(…, pcurve_tol)`), plus propriété d'une entité |

Deux constantes rejoignent `gbs/gbsconstants.h` (« single source of truth ») :

```cpp
template <typename T> inline constexpr T brep_default_tolerance<T> = T(1e-6);   // tol géométrique par défaut d'un vertex/edge/face créé par un builder
template <typename T> inline constexpr T brep_pcurve_approx_tol<T> = T(1e-5);   // déviation acceptée pour l'approximation d'une pcurve projetée
```

Le défaut `1e-6` reprend `BaseTopo` (OCCT utilise `Precision::Confusion() = 1e-7`
en mm ; gbs n'impose pas d'unité, question ouverte Q11). Une tolérance
`float` est légitime (`T = float`) mais les bindings Python restent en `double`
comme ailleurs.

### 4.3 Tolérance de sewing

`sew(model, faces, tol)` : deux arêtes libres sont appariées si leurs
extrémités sont à `tol` l'une de l'autre **et** si la distance maximale entre
les deux courbes sur un échantillon est `≤ tol`. L'arête conservée reçoit
`edge.tol = max(tol_e1, tol_e2, d_mesurée)` et les sommets fusionnés
`vertex.tol = max(tol_v1, tol_v2, dist(p1, p2)/2 + …)`, de façon que les
invariants de § 4.1 restent vrais sans déplacer la géométrie. C'est le
principe de `BRepBuilderAPI_Sewing` : on absorbe les écarts dans les
tolérances, on ne modifie pas les surfaces.

---

## 5. Explorateur et requêtes

### 5.1 Parcours par type

```cpp
namespace gbs::brep {
    // Sous-entités d'un shape, par type, sans doublon, dans l'ordre de découverte.
    template <typename SubId, std::floating_point T>
    auto explore(const Model<T> &m, ShapeId shape) -> std::vector<SubId>;

    // ex. : auto edges = explore<EdgeId>(m, FaceId{3});
    //       auto faces = explore<FaceId>(m, solid);
}
```

Équivalent de `TopExp_Explorer` + `TopExp::MapShapes`. La descente suit
`Compound → {Solid, Shell, Face, Wire, Edge, Vertex}`, `Solid → Shell → FaceUse → Face → Wire → CoEdge → Edge → Vertex`.
Les seams n'apparaissent qu'une fois dans `explore<EdgeId>` (dédoublonnage par
`EdgeId`), mais deux fois si on itère les co-arêtes d'un wire, ce qui est le
comportement attendu dans les deux cas. L'implémentation est un parcours avec
un `std::vector<bool>` de visite par type (les ids étant denses), sans
allocation de map.

### 5.2 Adjacence : index remontant

```cpp
namespace gbs::brep {
    template <std::floating_point T>
    class TopologyIndex   // construit à la demande à partir d'un shape racine
    {
    public:
        TopologyIndex(const Model<T> &m, ShapeId root);
        auto faces_of(EdgeId) const -> std::span<const FaceId>;       // 0, 1 (bord libre), 2 (variété), >2 (non-variété)
        auto edges_of(VertexId) const -> std::span<const EdgeId>;
        auto wires_of(EdgeId) const -> std::span<const std::pair<WireId, std::uint32_t>>; // (wire, index de la co-arête)
        auto shells_of(FaceId) const -> std::span<const ShellId>;
    };
}
```

Équivalent de `TopExp::MapShapesAndAncestors`. L'index est **invalidé par toute
mutation** du modèle ; il est construit une fois par builder qui en a besoin
(sewing, `make_solid`, `check`) et jamais stocké dans les entités (§ 2.3).

### 5.3 Requêtes de base

| Requête | Signature illustrative | Définition |
|---|---|---|
| Bornes | `bounding_box(m, shape, n_samples = 10) -> std::pair<point<T,3>, point<T,3>>` | Boîte englobante des sommets et d'un échantillon des arêtes (et des surfaces pour les faces à bornes naturelles), gonflée des tolérances. Réutilise `getCoordsMinMax` de `inc/topology/baseGeom.h`. |
| Fermeture d'un wire | `is_closed(m, WireId) -> bool` | Invariant 3 de § 2.6 |
| Fermeture d'un shell | `is_closed(m, ShellId) -> bool` | Toute arête non dégénérée a exactement deux usages |
| Variété | `is_manifold(m, ShellId) -> bool` | Toute arête a au plus deux usages |
| Orientabilité | `is_orientable(m, ShellId) -> bool` | Les deux usages de chaque arête sont de sens opposés après application des `FaceUse` |
| Arêtes libres | `free_edges(m, ShellId) -> std::vector<EdgeId>` | Usage unique |
| Point d'une co-arête | `coedge_point(m, wire, i, t)` | Évalue la pcurve en tenant compte de l'orientation |
| Volume signé | `signed_volume(m, ShellId, n_u, n_v)` | Théorème de la divergence sur un échantillonnage des faces ; sert à orienter le solide (§ 6.4) |

### 5.4 Validité minimale (« BRepCheck-lite »)

```cpp
namespace gbs::brep {
    enum class Issue : std::uint8_t {
        InvalidTolerance, EdgeWithoutCurve, EdgeBoundsInverted, VertexOffCurve,
        ToleranceOrderViolated,                // vertex < edge ou edge < face
        CoEdgeWithoutPCurve, PCurveRangeMismatch, SameParameterViolated,
        WireNotClosed, WireNotChained,
        OuterWireNotDirect, InnerWireNotIndirect, InnerWireOutsideOuter,
        NonManifoldEdge, ShellNotClosed, ShellNotOrientable, SolidNotOutward,
    };
    template <std::floating_point T>
    struct CheckReport { std::vector<std::pair<ShapeId, Issue>> issues; bool ok() const; };

    template <std::floating_point T>
    auto check(const Model<T> &m, ShapeId shape, CheckOptions opts = {}) -> CheckReport<T>;
}
```

Ce que `check()` **ne fait pas** (volontairement, c'est le « lite ») :
auto-intersection d'un wire dans (u,v), intersection entre faces d'un shell,
validité de la surface elle-même (ce sont `BRepCheck_Wire::SelfIntersect` et
`BOPAlgo_ArgumentAnalyzer` côté OCCT, hors périmètre). Les tests de § 10
construisent des cas invalides pour chaque `Issue`.

---

## 6. Builders

Tous les builders sont des fonctions libres `snake_case` prenant le `Model<T>&`
en premier argument et renvoyant un identifiant (ou un identifiant + un rapport), à
l'image de `build_segment`, `loft`, `interpolate`.

### 6.0 Sommets, arêtes, wires

```cpp
auto make_vertex(Model<T>&, const point<T,3>&, T tol = brep_default_tolerance<T>) -> VertexId;

// Bornes = curve->bounds(), sommets créés aux extrémités (fusionnés si curve fermée à tol près).
auto make_edge(Model<T>&, std::shared_ptr<Curve<T,3>> curve, T tol = …) -> EdgeId;
auto make_edge(Model<T>&, std::shared_ptr<Curve<T,3>> curve, T u1, T u2, T tol = …) -> EdgeId;
// Sommets imposés : vérifie dist(curve(u1), v1) ≤ max(tol, v1.tol), remonte v1.tol si besoin.
auto make_edge(Model<T>&, std::shared_ptr<Curve<T,3>> curve, T u1, T u2, VertexId v1, VertexId v2) -> EdgeId;
auto make_edge(Model<T>&, const point<T,3> &p1, const point<T,3> &p2) -> EdgeId; // segment, comme l'ancien Edge(pt1, pt2)
auto make_degenerate_edge(Model<T>&, VertexId v) -> EdgeId;

// Chaîne les arêtes bout à bout, fusionne les sommets à tol près (remplace Wire::addEdge),
// détermine l'orientation de chaque co-arête, ferme si le dernier sommet rejoint le premier.
// Les arêtes ne sont jamais modifiées ; seuls les sommets sont fusionnés (les arêtes pointent vers le survivant).
auto make_wire(Model<T>&, std::span<const EdgeId>, T tol = …) -> WireId;
```

`make_wire` reprend la logique de `Wire::addEdge` (recherche par les deux bouts,
fusion à tolérance près) mais sans muter les arêtes : la fusion consiste à
remplacer, dans les arêtes concernées, le `VertexId` absorbé par le survivant
(union-find sur les sommets), puis à marquer l'absorbé comme tombstone. Les
arêtes non chaînables sont signalées (exception `WireNotChained` ou rapport,
question ouverte Q13).

### 6.1 Face depuis une surface entière

```cpp
auto make_face(Model<T>&, std::shared_ptr<Surface<T,3>> srf, T tol = …) -> FaceId;
```

Algorithme :

1. `[u1,u2,v1,v2] = srf->bounds()` ; `c = surface_closure(*srf, tol)` (§ 3.3).
2. Quatre pcurves : `Line<T,2>` ou `BSCurve<T,2>` de degré 1 sur les bords
   `v = v1`, `u = u2`, `v = v2`, `u = u1`, orientées pour tourner dans le sens
   direct de (u,v).
3. Quatre courbes 3D : `CurveOnSurface<T,3>(pcurve, srf)` (exactes, `same_parameter`
   sans déviation). Pour une `BSSurface`, on peut préférer `isoU`/`isoV` qui
   renvoient une vraie `BSCurve` (meilleure pour IGES) ; les deux sont
   équivalentes géométriquement.
4. Sommets : les quatre coins, fusionnés deux à deux si `c.closed_u`/`c.closed_v`,
   réduits à un seul si bord dégénéré.
5. Arêtes : un bord dégénéré donne `make_degenerate_edge` ; un bord fermé
   (`closed_u`) donne **une** arête pour `u = u1` réutilisée en `Reversed` avec la
   pcurve `u = u2` (seam, § 3.4) ; sinon quatre arêtes ordinaires.
6. Wire fermé de 4 co-arêtes (toujours 4, certaines dégénérées ou seam),
   `Face{srf, {wire}, tol, natural_bounds = true}`.

Cas couverts par les tests : plan NURBS (4 arêtes), cylindre (seam + 2 arêtes),
sphère NURBS (seam + 2 dégénérées), tore (2 seams), `SurfaceOfRevolution`
complète (seam) et partielle (4 arêtes).

Comparaison : `BRepBuilderAPI_MakeFace(surface)` ; OCCT traite le seam par
`BRep_CurveOnClosedSurface`, ici par deux co-arêtes.

### 6.2 Face depuis un wire posé sur une surface

```cpp
struct MakeFaceOptions { T pcurve_tol = brep_pcurve_approx_tol<T>; size_t pcurve_degree = 3; size_t n_samples_min = 10; bool check_orientation = true; };

// Wire déjà fermé, chaque arête à pcurve_tol de la surface.
auto make_face(Model<T>&, std::shared_ptr<Surface<T,3>> srf, WireId outer, MakeFaceOptions = {}) -> FaceId;
auto make_face(Model<T>&, std::shared_ptr<Surface<T,3>> srf, WireId outer, std::span<const WireId> holes, MakeFaceOptions = {}) -> FaceId;

// Variante « palier 2 » : contours donnés en 2D ; les courbes 3D sont des CurveOnSurface exactes.
auto make_face(Model<T>&, std::shared_ptr<Surface<T,3>> srf, std::span<const std::shared_ptr<Curve<T,2>>> outer_pcurves, MakeFaceOptions = {}) -> FaceId;
```

Calcul de la pcurve de chaque co-arête, par ordre de préférence :

1. **Extraction exacte.** Si `edge.curve` est (éventuellement à travers
   `CurveTrimmed`, `CurveReversed`) une `CurveOnSurface<T,3>` dont
   `p_basisSurface()` est le même `shared_ptr` que `srf` : `pcurve = basisCurve`
   (reparamétrée/inversée pour respecter § 3.2 ; `CurveReversed` existe,
   `CurveReparametrized` aussi). Déviation nulle.
2. **Projection + approximation.** Sinon :
   - échantillonner l'arête sur `[u1,u2]` (paramètres `t_i`, nombre adapté par la
     déviation, `discretize(crv, n, dev)` existe) ;
   - pour chaque `t_i`, inverser `(u,v)_i = extrema_surf_pnt(srf, C(t_i), …)` avec
     **amorce** au `(u,v)` précédent pour rester sur la bonne branche et éviter
     le balayage global ; sur une surface fermée, **dérouler** le paramètre
     périodique (si `|u_i − u_{i−1}| > période/2`, ajouter `±période`) ;
   - rejeter si `‖S(u_i,v_i) − C(t_i)‖ > pcurve_tol` (l'arête n'est pas sur la
     surface) ;
   - interpoler une `BSCurve<T,2>` aux `(t_i, (u,v)_i)` (paramètres imposés =
     `t_i`, ce qui donne SameParameter par construction aux nœuds ; `interpolate`
     avec paramètres existe) ou l'approximer si le nombre de points est grand
     (`approx`) ;
   - mesurer la déviation aux milieux des intervalles ; raffiner (nouveaux
     échantillons) jusqu'à `≤ pcurve_tol` ; `edge.tol = max(edge.tol, dev)`.
3. **Cas d'une arête sur un seam** : si tous les `(u,v)_i` sont à `tol` près de
   `u = u1 ≡ u2`, le côté est ambigu. On choisit le côté qui rend le wire
   **continu** avec la co-arête précédente (dont le `u` final est connu). Si le
   wire entier est sur le seam, échec explicite : c'est le rôle de § 6.1.

Puis :

4. **Orientation du wire.** Aire signée du polygone des pcurves (échantillon) :
   positive pour `outer`, sinon on **inverse l'ordre et le sens des co-arêtes**
   (jamais la face). Négative exigée pour chaque trou, inversé de même. Test
   d'inclusion des trous dans `outer` par parité de croisements sur le polygone
   (point-dans-polygone, `inc/topology/baseGeom.h` a `orient_2d`).
5. `natural_bounds` détecté si le wire extérieur coïncide avec le rectangle
   paramétrique (quatre co-arêtes iso aux bornes).

Comparaison : `BRepBuilderAPI_MakeFace(surface, wire)` suivi de
`BRepLib::BuildPCurveForEdgeOnPlane` / `ShapeFix_Edge::FixAddPCurve` ; OCCT
projette avec `ProjLib`/`GeomProjLib` puis `Approx_CurveOnSurface`. Ici la
projection repose sur `extrema_surf_pnt` et l'interpolation/approximation gbs.

### 6.3 Shell par sewing

```cpp
struct SewOptions { T tol = brep_default_tolerance<T>; size_t n_samples = 10; bool orient = true; bool allow_non_manifold = false; };

struct SewReport {
    std::vector<ShellId> shells;                       // un par composante connexe
    std::vector<EdgeId>  free_edges;                   // bords restés libres
    std::vector<EdgeId>  non_manifold_edges;           // > 2 usages
    std::vector<std::pair<EdgeId, EdgeId>> rejected;   // candidats proches mais non appariés (distance échantillonnée > tol)
    bool orientable{true};
};

auto sew(Model<T>&, std::span<const FaceId>, SewOptions = {}) -> SewReport;
```

Algorithme :

1. **Candidats.** `TopologyIndex` sur les faces ; arêtes **libres** (un seul usage,
   non dégénérées) = candidats. Les deux co-arêtes d'un seam partagent déjà un
   `EdgeId` : elles ne sont pas libres, donc jamais candidates.
2. **Filtre grossier.** Boîte englobante de chaque candidat (sommets + échantillon)
   gonflée de `tol` ; paires dont les boîtes se coupent. Tri des boîtes sur un
   axe (sweep) pour rester en `O(n log n)` ; pas de kd-tree au palier 1.
3. **Appariement par extrémités.** Pour `(e1, e2)` : `v1₁ ≈ v1₂ ∧ v2₁ ≈ v2₂`
   (même sens) ou `v1₁ ≈ v2₂ ∧ v2₁ ≈ v1₂` (sens inverse), à `tol` près (distance
   des points, pas identité des identifiants). Arêtes fermées (`v1 == v2`) : on
   compare les points et on tranche le sens à l'étape 4.
4. **Appariement par échantillonnage.** `n_samples` points `C1(t_i)` ;
   `d_i = extrema_curve_point(C2, C1(t_i))` (existe, amorçable) ; **et
   réciproquement** `C2 → C1` pour rejeter un recouvrement partiel (jonction en
   T : une arête courte sur une longue passe le test dans un sens seulement).
   Accepté si `max d ≤ tol`. Les paramètres des pieds de projection donnent le
   sens relatif (croissant ou décroissant).
5. **Fusion.** `e1` (plus petit id) survit ; `e1.tol = max(e1.tol, e2.tol, max d)`.
   Toutes les co-arêtes référençant `e2` pointent vers `e1`, avec
   `orient` basculé si le sens relatif est inverse ; leur pcurve est reparamétrée
   affinement sur `[u1, u2]` de `e1` (et inversée par `CurveReversed` si besoin)
   pour respecter § 3.2. Les sommets correspondants sont fusionnés (union-find)
   avec remontée des tolérances (§ 4.3). `e2` et les sommets absorbés deviennent
   tombstones. Si `e1` a déjà deux usages et `allow_non_manifold == false`, la
   paire est refusée et signalée dans `non_manifold_edges`.
6. **Orientation cohérente.** Graphe faces–faces par arêtes partagées ; parcours
   en largeur par composante connexe depuis une face de référence
   (`FaceUse::Forward`). Pour un voisin atteint par `e` : cohérent ssi les deux
   co-arêtes de `e` sont parcourues en **sens opposés** une fois les `FaceUse`
   appliqués ; sinon le `FaceUse` du voisin est basculé en `Reversed`. Si une
   face déjà visitée est atteinte avec une contrainte contradictoire, la
   composante est non orientable (ruban de Möbius) : `orientable = false`, le
   shell est tout de même produit.
7. **Sortie.** Un `Shell` par composante connexe, `closed` évalué ; `free_edges`
   et `rejected` remplis. L'appelant décide (lever, continuer, faire un
   compound).

Limites assumées (restent côté OCCT `BRepBuilderAPI_Sewing` avec *cutting*) :
pas de découpe d'arête pour les jonctions en T, pas de sommet sur une arête,
pas de fusion de faces coplanaires, pas de déplacement de géométrie.

### 6.4 Solide depuis un shell fermé, compound

```cpp
auto make_solid(Model<T>&, ShellId outer) -> SolidId;                  // lève si !is_closed ou !is_orientable
auto make_solid(Model<T>&, ShellId outer, std::span<const ShellId> voids) -> SolidId; // cavités : stockées, non vérifiées au palier 1 au-delà de closed
auto make_compound(Model<T>&, std::span<const ShapeId>) -> CompoundId;
```

Orientation du solide : `signed_volume(shell)` par le théorème de la divergence
(`V = ⅓ ∮ S·n dA`) sur un échantillonnage `(n_u × n_v)` de chaque face avec
intégration en (u,v) restreinte au contour par le test point-dans-polygone des
pcurves ; si `V < 0`, **tous** les `FaceUse` du shell sont basculés. L'estimation
n'a pas besoin d'être précise, seulement de bon signe ; pour une face à bornes
naturelles l'intégrale est directe. OCCT fait `BRepLib::OrientClosedSolid` via
un classifieur de point à l'infini ; la divergence est plus simple à écrire
avec ce que gbs offre et robuste sur des faces découpées.

---

## 7. Export IGES

### 7.1 Ce que libIGES permet aujourd'hui

`gbs-io/iges.h` écrit déjà **126** (courbe NURBS), **128** (surface NURBS),
**120 + 110** (surface de révolution + axe) via l'API `DLL_IGES`. Dans la
version installée (headers `iges/api/`), l'API `DLL_IGES` expose en plus :
**100** (arc), **102** (courbe composite), **104**, **110** (ligne), **122**,
**124** (matrice), **142** (courbe sur surface paramétrique), **144** (surface
paramétrique découpée), **308/408** (sous-figure), **314** (couleur),
**406** (propriété, dont le nom en forme 15 déjà utilisé par `SetLabel`).

Les entités BREP **186** (solide variété), **514** (shell), **510** (face),
**508** (loop), **504** (liste d'arêtes), **502** (liste de sommets) existent
dans l'API **cœur** (`iges/core/entity186.h`…) mais **pas** dans `DLL_IGES`.
Leur usage obligerait soit à passer par la classe cœur `IGES` (symboles non
garantis exportés sur Windows), soit à étendre libIGES en amont. L'entité
**402** (groupe) n'est pas disponible.

### 7.2 Correspondance proposée au palier 1

| Entité gbs | IGES | Détail |
|---|---|---|
| `Face` | **144** | `PTS` = surface, `N1 = natural_bounds ? 0 : 1`, `PTO` = 142 du wire extérieur, `PTI[]` = 142 des trous |
| `Wire` de face | **142** | `SPTR` = surface, `BPTR` = **102** composite des pcurves en **126** (coordonnées `(u, v, 0)`), `CPTR` = **102** composite des courbes 3D en **126**, `PREF = 3` (les deux représentations équivalentes) ; `CRTN = 1` si la pcurve a été projetée, `0` sinon |
| `CoEdge` `Reversed` | — | pcurve et courbe 3D exportées **inversées** (`BSCurveGeneral::reversed()` ; sinon `CurveReversed` + approximation) |
| Arête dégénérée | — | pcurve en 126 dans `BPTR` ; dans `CPTR` une 126 réduite à un point (pôles confondus) : convention acceptée par les lecteurs usuels, à vérifier à l'import OCCT dans les tests |
| `Face::surface` | **128** si `BSSurfaceGeneral` ; **120** si `SurfaceOfRevolution` (chemin existant) ; sinon **approximation** en `BSSurface` (`bssapprox.h`) puis 128 | Les courbes 3D non NURBS sont approchées en `BSCurve` (`approx`) avant 126, comme le fait déjà `add_geom(SurfaceOfRevolution)` pour sa génératrice |
| `Shell`, `Solid` | ensemble de **144**, une par face, nommées `"<nom>/face_<i>"` via 406 forme 15 | 186/514/… reportés (question ouverte Q10) |
| `Compound` | récursion sur ses éléments | — |
| `Vertex`, `Edge`, `Wire` libres | 116 (point) non exposé par DLL : ignoré ; arête → 126 ; wire → 102 | — |

Entrée dans `IgesWriter` : `add_geometry(const brep::Model<T>&, ShapeId, const std::string &name = "")`
à côté des surcharges existantes. Unités : la section globale IGES porte
l'unité ; gbs n'en impose pas. On ajoute un paramètre `scale` (défaut 1) comme
`gbs-occt/export.h` (qui exporte en mm avec `scale = 1000`), appliqué aux pôles
et tolérances à l'écriture.

Le test de référence : exporter la boîte, le cylindre, la sphère et un solide
cousu, relire avec OCCT dans `gbs-occt/tests` (`IGESControl_Reader`) quand
`GBS_USE_OCCT_UTILS` est actif, et vérifier nombre de faces et validité
`BRepCheck_Analyzer`. Sans OCCT, test de relecture par libIGES (`DLL_IGES::Read`)
et comptage d'entités.

---

## 8. API C++ et API Python

### 8.1 Conventions C++

- **Dimension** : `template <std::floating_point T>` uniquement, 3D fixé.
  Justification : `Surface<T,3>` est la seule surface qui ait un sens
  (`SurfaceOfRevolution` est déjà `Surface<T,3>`), IGES et STEP sont 3D, une
  face « 2D » est une face sur un plan `z = 0`. Garder `dim` doublerait les
  instanciations Python sans usage, et rendrait les pcurves (`Curve<T,2>`)
  ambiguës en 2D. Les builders de wires acceptent en revanche `Curve<T,3>`
  construites à partir de 2D par `add_dimension` si besoin.
- **Nommage** : classes `CamelCase` (`Model`, `Face`), identifiants `XxxId`, fonctions
  libres `snake_case` (`make_face`, `sew`, `explore`, `check`, `bounding_box`),
  méthodes `camelCase` sur `Model` (`vertex()`, `addFace()`…) comme `knotsFlats()`.
- **Erreurs** : exceptions dérivées de `std::runtime_error` dans `gbs/exceptions.h`
  pour les préconditions (`BRepError`, `WireNotChained`, `ShellNotClosed`) ;
  rapports (`SewReport`, `CheckReport`) pour les diagnostics non bloquants.
- **Header-only** : `gbs-brep/*.h` avec `#include <gbs-brep/...>`, en-tête
  agrégateur `<gbs-brep/brep>` sur le modèle de `<gbs/curves>`. Pas de raison
  forte de compiler : aucune dépendance nouvelle, les algorithmes sont
  templates sur `T`. Compatible `GBS_USE_PCH`.
- **Fichiers** : `model.h` (identifiants, entités, `Model`), `explore.h`
  (`explore`, `TopologyIndex`, requêtes), `check.h`, `builders.h` (`make_*`),
  `sew.h`, `pcurve.h` (extraction/projection), `closure.h` (`surface_closure`,
  `signed_volume`), `gbs-io/iges.h` étendu pour l'export.

Exemple d'usage visé :

```cpp
using namespace gbs;
brep::Model<double> m;
auto cyl  = std::make_shared<BSSurfaceRational<double,3>>(/* cylindre exact */);
auto top  = std::make_shared<BSSurface<double,3>>(/* disque */);
auto bot  = std::make_shared<BSSurface<double,3>>(/* disque */);
auto f1 = brep::make_face(m, cyl);
auto f2 = brep::make_face(m, top, brep::make_wire(m, {brep::make_edge(m, circle_top)}));
auto f3 = brep::make_face(m, bot, brep::make_wire(m, {brep::make_edge(m, circle_bot)}));
auto rep = brep::sew(m, std::array{f1, f2, f3}, {.tol = 1e-6});
auto so  = brep::make_solid(m, rep.shells.front());
assert(brep::check(m, so).ok());
IgesWriter<double> w; w.add_geometry(m, so, "cylinder"); w.write("cylinder.igs");
```

### 8.2 API Python

- Sous-module `gbs.brep` (`m.def_submodule("brep")` dans `gbsbind.cpp`, fichier
  `python/gbsBindBrep.cpp` + `.h` sur le modèle de `gbsBindCurves`), `T = double`.
- `Model` liée en `std::shared_ptr` ; identifiants liés comme petites classes
  (`VertexId`…, `__int__`, `__eq__`, `__hash__`, `__repr__`) ; `ShapeId` converti
  automatiquement par un caster `std::variant` de pybind11.
- Les entités (`Vertex`, `Edge`, `Face`…) sont exposées en lecture
  (`m.vertex(vid).pnt`, `m.edge(eid).curve`, `m.face(fid).wires`) ; la géométrie
  revient comme les classes Python existantes (`BSCurve3d`, `BSSurface3d`,
  `CurveOnSurface3d`…) grâce au polymorphisme déjà en place sur
  `Curve<T,3>`/`Surface<T,3>`.
- Fonctions libres identiques au C++ : `make_vertex`, `make_edge`, `make_wire`,
  `make_face`, `sew`, `make_solid`, `make_compound`, `explore_edges(m, shape)` /
  `explore_faces(m, shape)` / `explore_vertices(m, shape)` (le template
  `explore<Sub>` devient une fonction par type), `faces_of`, `edges_of`,
  `is_closed`, `is_manifold`, `bounding_box`, `check`.
- `IgesWriter.add_geometry(model, shape, name="")` ajoutée.
- Pas de suffixe `_3d` : le sous-module est 3D par construction.
- `__repr__` JSON via `repr.h` pour `Model` (compte d'entités) ; la sérialisation
  JSON complète du modèle (`gbs-io/tojson.h`) est hors palier 1.

```python
from pygbs import gbs
m = gbs.brep.Model()
f = gbs.brep.make_face(m, srf)                       # face à bornes naturelles
rep = gbs.brep.sew(m, [f1, f2, f3], tol=1e-6)
so = gbs.brep.make_solid(m, rep.shells[0])
assert gbs.brep.check(m, so).ok()
for e in gbs.brep.explore_edges(m, so):
    print(m.edge(e).curve.bounds(), gbs.brep.faces_of(m, so, e))
w = gbs.IgesWriter(); w.add_geometry(m, so, "part"); w.write("part.igs")
```

Tests Python dans `python/tests/test_brep.py` (pytest, comme les existants).

---

## 9. Points de contact avec le palier 2 et la lecture STEP

### 9.1 Réservé pour le palier 2

| Besoin du palier 2 | Ce que le palier 1 prévoit |
|---|---|
| Courbes d'intersection portées par deux faces | `CoEdge::pcurve` par face ; `Edge::curve` peut être une `CurveOnSurface` sur l'une des deux surfaces ou une `BSCurve` approchée, `same_parameter` et `tol` absorbent l'écart |
| Découpe d'une face : nouveaux wires, trous | `Face::wires` à plusieurs éléments avec convention extérieur/trous (§ 2.6 inv. 4), `make_face(srf, pcurves2d)` pour reconstruire une face depuis des contours 2D exacts |
| Découpe d'une arête partagée | Identifiants : remplacer `EdgeId` dans les `CoEdge` des deux faces ; `TopologyIndex::wires_of(edge)` donne où ; `Model::erase` + `compact()` |
| Partage d'une surface entre les morceaux d'une face découpée | `Face::surface` en `shared_ptr` partagé, bornes de la face données par ses wires et non par la surface |
| Sweep / révolution en BREP | `make_face(srf)` sur la surface balayée + `sew` avec les faces d'extrémités ; `surface_closure` détecte le seam de révolution |
| Offsets de faces | `SurfaceOffset` existe ; une face offset = `Face{SurfaceOffset(srf), pcurves identiques}` si les pcurves sont réutilisées telles quelles (même paramétrage), avec remontée de `tol` |
| Classification point / solide | `signed_volume` et le point-dans-polygone (u,v) sont les briques ; un classifieur par lancer de rayon viendra avec les intersections courbe-surface |

### 9.2 Réservé pour la lecture STEP

| STEP AP203/AP214 | gbs::brep |
|---|---|
| `cartesian_point`, `vertex_point` | `Vertex` |
| `edge_curve(start, end, curve, same_sense)` | `Edge` ; `same_sense = false` ⇒ `CurveReversed` enveloppe la courbe pour garder `u1 < u2` dans le sens `v1 → v2` |
| `oriented_edge(edge, orientation)` dans `edge_loop` | `CoEdge{edge, orient}` dans `Wire` |
| `pcurve` / `surface_curve` / `seam_curve` | `CoEdge::pcurve` ; `seam_curve` ⇒ deux co-arêtes sur la même `EdgeId` |
| `face_outer_bound`, `face_bound` (+ `orientation`) | `Face::wires[0]`, suivants ; l'orientation du `face_bound` renverse l'ordre du wire |
| `advanced_face(bounds, surface, same_sense)` | `Face` ; `same_sense` → `FaceUse::orient` dans le shell |
| `plane`, `cylindrical_surface`, … | conversion en `BSSurfaceRational<T,3>` exactes (décision de pilotage), bornes paramétriques issues des pcurves ou des bornes naturelles |
| `closed_shell`, `open_shell` | `Shell` |
| `manifold_solid_brep(outer)`, `brep_with_voids(voids)` | `Solid{outer, voids}` (le champ `voids` existe dès le palier 1) |
| `#id`, `name` des entités | Tables annexes `Model::setName(ShapeId, string)` / `name()`, `setExternalId(ShapeId, int64)` / `externalId()` : `std::unordered_map` par type, vides par défaut, sans coût pour les modèles natifs. L'export IGES utilise `name()` pour 406 |
| Unités (`length_unit`, `si_unit`) | `Model::unitScale` (facteur vers l'unité de travail, défaut 1) renseigné par le lecteur, utilisé par l'export IGES |
| Tolérances `uncertainty_measure_with_unit` | Initialise `brep_default_tolerance` du lecteur ; les tolérances par entité sont ensuite recalculées par `check`/builders |

Les identifiants STEP ne sont **pas** les indices de l'arène (ils ne sont pas
denses) : la table `externalId` fait le lien dans les deux sens.

---

## 10. Plan de développement en PR

Chaque PR est header-only, testée dans `tests/tests_brep_*.cpp` (doctest via
`testing/doctest_gtest.hpp`, macro `TEST(suite, nom)`), et ajoute une entrée à
`News.md` sous « Unreleased ». Les lignes indiquées comptent code + tests.
« h » = heures de session de développement, revue comprise.

| # | PR | Contenu | Tests | Lignes | h |
|---|---|---|---|---|---|
| 1 | `brep/model` | `gbs-brep/model.h` : identifiants, `Orientation`, entités de § 2.5, `Model` (accès, `add`, `erase`, `compact`, `append`), constantes dans `gbsconstants.h`, CMake `gbs-brep/` + `INSTALL_HEADERS`, en-tête agrégateur | Construction manuelle d'une boîte (8 sommets, 12 arêtes, 6 faces planes NURBS, 1 shell, 1 solide) ; `compact()` renumérote ; `append` | 450 | 4 |
| 2 | `brep/explore` | `explore.h` : `explore<Sub>`, `TopologyIndex`, `is_closed`, `is_manifold`, `is_orientable`, `free_edges`, `bounding_box`, `coedge_point` | Sur la boîte de PR 1 : comptes (8/12/6), `faces_of` = 2 partout, boîte ouverte (une face retirée) → 4 arêtes libres, bbox | 400 | 3 |
| 3 | `brep/builders-wire` | `builders.h` : `make_vertex`, `make_edge` (4 surcharges), `make_degenerate_edge`, `make_wire` avec fusion de sommets ; réécriture de `tests/tests_topo.cpp` en `tests_brep_wire.cpp` ; anciens en-têtes marqués `[[deprecated]]` | Les 3 tests actuels portés ; wire ouvert, fermé, arête fermée (cercle), arêtes dans le désordre, arête non chaînable | 400 | 3 |
| 4 | `brep/face-natural` | `closure.h` (`surface_closure`), `make_face(srf)` : seams, dégénérées, 4 cas | Plan, cylindre rationnel, sphère NURBS, tore, révolution complète/partielle : nombre d'arêtes distinctes, `check` ok, pcurves aux bornes | 450 | 4 |
| 5 | `brep/face-wire` | `pcurve.h` (extraction exacte, projection + interpolation/approximation, déroulement périodique, côté du seam), `make_face(srf, wire[, holes])`, `make_face(srf, pcurves2d)`, orientation des wires, `natural_bounds` | Disque sur plan, contour quelconque sur surface libre (déviation ≤ `pcurve_tol`), face avec trou (orientation corrigée automatiquement), arête traversant le seam d'un cylindre, arête hors surface → rejet | 550 | 5 |
| 6 | `brep/check` | `check.h` : toutes les `Issue` de § 5.4, `CheckReport` | Un cas construit pour chaque `Issue` ; les sorties de PR 1–5 passent | 400 | 3 |
| 7 | `brep/sew` | `sew.h` : candidats, filtre par boîtes, extrémités, échantillonnage bidirectionnel, fusion avec reparamétrage des pcurves, orientation par parcours, `SewReport` | 6 faces de boîte en désordre et orientations aléatoires → 1 shell fermé orientable ; cylindre + 2 disques ; deux composantes → 2 shells ; jonction en T → `rejected` ; Möbius → `orientable = false` ; faces avec écart `0.5·tol` → tolérances remontées | 600 | 6 |
| 8 | `brep/solid` | `signed_volume`, `make_solid` (orientation sortante, cavités stockées), `make_compound` | Boîte retournée → volume > 0 après ; cylindre cousu ; shell ouvert → exception ; compound mixte exploré | 350 | 3 |
| 9 | `brep/iges` | `gbs-io/iges.h` : 102, 142, 144, approximation des géométries non NURBS, arêtes dégénérées, `add_geometry(model, shape, name)`, `scale` | Export boîte / cylindre / sphère / solide cousu ; relecture libIGES (comptes) ; relecture OCCT + `BRepCheck_Analyzer` sous `GBS_USE_OCCT_UTILS` | 450 | 4 |
| 10 | `brep/python` | `python/gbsBindBrep.{h,cpp}`, sous-module `gbs.brep`, `IgesWriter.add_geometry`, stubs ; `python/tests/test_brep.py` ; suppression des anciens `inc/topology/{basetopo,vertex,edge,wire}.h` ; documentation (`docs/index.rst` : `doxygenfile` des nouveaux en-têtes, `Doxyfile` `INPUT += ../gbs-brep`) ; `News.md` | pytest : boîte cousue, exploration, export IGES, erreurs levées | 500 | 4 |

Total : environ 4 550 lignes, **39 h**. Ordre imposé : 1 → 2 → 3 → 4 → 5 → 6 →
7 → 8 → 9 → 10 ; les PR 6 et 9 peuvent démarrer dès la PR 5 fusionnée.

---

## 11. Questions ouvertes

Pour chaque question : **R** = recommandation, **A** = alternatives.

**Décision du 3 octobre 2026 : toutes les recommandations sont retenues.** Les
alternatives restent listées pour mémoire.

1. **Remplacer ou étendre l'embryon `inc/topology/{basetopo,vertex,edge,wire}.h` ?**
   R : remplacer (§ 2.1), supprimer les anciens en-têtes à la PR 10 après
   migration de `tests_topo.cpp`. A : les garder en façade dépréciée sur
   `gbs::brep` pendant une version ; les étendre en place (déconseillé :
   couplage au half-edge mesh et `tessellate()` pur virtuel).

2. **Arène + identifiants ou graphe de `shared_ptr` ?**
   R : arène `Model<T>` + identifiants typés (§ 2.3). A : `shared_ptr` descendants +
   `weak_ptr` remontants (plus proche de l'existant et des classes
   géométriques, mais mutations topologiques et identité plus fragiles).

3. **Orientation par usage (`CoEdge`/`FaceUse`) ou par valeur de shape (OCCT) ?**
   R : par usage (§ 2.4), pcurve sur la `CoEdge`. A : `Shape = (id, orientation)`
   et pcurves listées sur l'arête, indexées par face.

4. **BREP 3D seulement ou `<T, dim>` ?**
   R : 3D seulement (§ 8.1). A : garder `dim` pour un BREP 2D (régions planes
   pour le mailleur 2D) ; la face sur plan `z = 0` couvre ce besoin.

5. **Placement / transformation par shape (`TopLoc_Location`) ?**
   R : aucun au palier 1 ; les transformations s'appliquent à la géométrie.
   A : un `Transform` optionnel par `FaceUse`/`Compound` pour les instances
   d'assemblage STEP (économise de la mémoire sur les pièces répétées).

6. **Classes de surfaces analytiques (`Plane`, `Cylinder`, `Cone`, `Sphere`, `Torus`) ?**
   R : non au palier 1 ; NURBS rationnelles exactes, polymorphie `Surface<T,3>`
   préservée pour les ajouter plus tard (§ 3.3). A : les ajouter dès maintenant
   pour préparer des intersections analytiques au palier 2.

7. **Détection de fermeture/périodicité des surfaces : numérique ou virtuelle ?**
   R : fonction libre numérique `surface_closure(srf, tol)` (§ 3.3), sans toucher
   `Surface<T,3>`. A : méthodes virtuelles `isUClosed()`/`isVClosed()` sur
   `Surface<T,3>` avec implémentation par défaut numérique, surchargées par
   `BSSurfaceGeneral` (pôles) et `SurfaceOfRevolution` (angle = 2π).

8. **Emplacement et namespace : `gbs-brep/` + `gbs::brep`, ou `inc/topology/` + `gbs` ?**
   R : `gbs-brep/` et `gbs::brep` (§ 2.1), au même rang que `gbs-io`/`gbs-mesh`.
   A : `inc/topology/brep/` dans `gbs` avec renommage des structures half-edge
   pour lever les ambiguïtés.

9. **Python : sous-module `gbs.brep` ou classes à plat ?**
   R : sous-module (§ 8.2), sans suffixe `_3d`. A : à plat avec préfixe `BRep`
   (`BRepModel`, `make_brep_face`…), cohérent avec l'absence actuelle de
   sous-modules.

10. **Export IGES 186 (solide BREP) au palier 1 ?**
    R : non ; faces en 144 nommées, 186 reporté (§ 7). A : utiliser l'API cœur
    de libIGES (`IGES`, `IGES_ENTITY_186/514/510/508/504/502`) au risque de
    symboles non exportés sur Windows ; ou proposer en amont l'exposition
    `DLL_IGES` de ces entités.

11. **Tolérance par défaut : `1e-6` (BaseTopo actuel) ou `1e-7` (OCCT, en mm) ?**
    R : `1e-6` dans `gbsconstants.h`, surchargeable partout (§ 4.2). A : `1e-7` ;
    ou une tolérance relative à la boîte englobante du modèle.

12. **Sewing : jonctions en T et recouvrements partiels au palier 1 ?**
    R : non, signalés dans `SewReport::rejected`, délégués à OCCT (§ 6.3).
    A : découper l'arête longue au pied de projection (petit pas vers le
    healing, mais ouvre la porte à la découpe d'arêtes avant le palier 2).

13. **Échec des builders : exception ou rapport ?**
    R : exception pour une précondition violée (`make_solid` sur shell ouvert,
    `make_wire` non chaînable, arête hors surface), rapport pour les
    diagnostics partiels (`sew`, `check`). A : tout en rapports avec identifiant
    invalide, style `BRepBuilderAPI_MakeShape::IsDone()`.

14. **`Wire` unique ou `Wire` + `Loop` ?**
    R : un seul type `Wire` (libre ou contour de face, § 2.5). A : `Loop`
    distinct avec pcurves obligatoires, `Wire` libre sans pcurves (plus strict,
    une conversion de plus à chaque `make_face`).

15. **Cavités (`Solid::voids`) au palier 1 ?**
    R : le champ existe et est stocké, `make_solid(outer, voids)` ne vérifie que
    la fermeture des cavités ; l'inclusion et l'orientation entrante sont
    vérifiées au palier 2 avec la classification point/solide. A : retirer le
    champ jusqu'à la lecture STEP.

---

## Annexe A — Correspondance des vocabulaires

| gbs::brep | OCCT | Parasolid | STEP |
|---|---|---|---|
| `Model` | — (graphe de `TShape`) | partition | fichier / `shape_representation` |
| `Compound` | `TopoDS_Compound` | — (plusieurs bodies) | `shape_representation` à plusieurs items |
| `Solid` | `TopoDS_Solid` | body (solid) / region | `manifold_solid_brep`, `brep_with_voids` |
| `Shell` | `TopoDS_Shell` | shell | `closed_shell`, `open_shell` |
| `FaceUse` | orientation de la `TopoDS_Face` dans le shell | sens de la face dans le shell | `oriented_face`, `advanced_face.same_sense` |
| `Face` | `TopoDS_Face` / `BRep_TFace` | face | `advanced_face` |
| `Wire` | `TopoDS_Wire` | loop | `edge_loop`, `face_bound` |
| `CoEdge` | orientation de la `TopoDS_Edge` dans le wire + `BRep_CurveOnSurface` | fin (+ SP-curve) | `oriented_edge` + `pcurve` |
| `Edge` | `TopoDS_Edge` / `BRep_TEdge` | edge | `edge_curve` |
| `Vertex` | `TopoDS_Vertex` / `BRep_TVertex` | vertex | `vertex_point` |
| `tol` | `BRep_Tool::Tolerance` | tolerant entities | `uncertainty_measure_with_unit` (global) |
| `TopologyIndex` | `TopExp::MapShapesAndAncestors` | pointeurs arrière natifs | — |
| `check` | `BRepCheck_Analyzer` (partiel) | `PK_BODY_check` (partiel) | — |
| `sew` | `BRepBuilderAPI_Sewing` sans *cutting* | `PK_BODY_sew_bodies` (partiel) | — |
