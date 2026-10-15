# Palier 2 — PR 9 : booléens (configurations génériques)

| | |
|---|---|
| PR | branche `feat/brep-boolean` |
| Document de référence | [brep_palier2.md](brep_palier2.md), § 6 ; questions 1 et 5 |
| Fichiers | `gbs-brep/boolean.h` (nouveau), `gbs-brep/builders.h` (`BuildErrc::NonGeneric`), `gbs-brep/brep` (agrégateur), `tests/tests_brep_boolean.cpp` (nouveau), `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Fusion, différence et intersection de solides, pour les configurations
**génériques** : les faces des deux opérandes se coupent transversalement.

Le document de conception mettait aussi dans cette PR les faces planes
coplanaires et les primitives coaxiales confondues (§ 6.3). **Je les sépare
dans une PR suivante (9 b)**, avant la PR 10. Elles demandent une découpe 2D
commune des faces confondues : projection des bords de l'une sur l'autre,
découpe croisée des arêtes, et fusion des arêtes superposées. Cela fait une
PR à part entière. En attendant, elles rendent `NonGeneric`, comme les autres
cas dégénérés.

La PR 8 (extrusion, révolution, faces décalées), indépendante, viendra
ensuite.

## 2. API

```cpp
enum class BooleanOp { Fuse, Cut, Common };
struct BooleanOptions { T tol = brep_default_tolerance; IntersectOptions<T> intersect; };
enum class FaceOrigin { A, B };
struct BooleanResult { ShapeId shape; vector<SolidId> solids; vector<pair<FaceId, FaceOrigin>> faces; };

auto boolean(m, BooleanOp op, ShapeId a, ShapeId b, BooleanOptions = {}) -> BuildResult<BooleanResult>;
auto fuse(m, a, b, opts), cut(m, a, b, opts), common(m, a, b, opts);
```

- `a` et `b` : un solide, ou un compound de solides.
- `shape` : un solide, ou un compound s'il y en a plusieurs. Un résultat vide
  donne un compound vide.

## 3. Choix d'architecture et justification

### 3.1 Sur des copies

Les opérandes sont copiés (`copy_shape`, placement identité). **Ils ne sont
jamais modifiés.** Le résultat est fait des faces des copies. En cas d'échec,
les copies sont effacées et le modèle retrouve ses formes d'origine.

### 3.2 Étapes

1. **Section** des copies (PR 6) : chaque arête de section est créée une seule
   fois.
2. **Découpe** des faces **des deux copies** par les mêmes arêtes (PR 7), avec
   la pcurve de chaque côté. Après la découpe, chaque arête de section est
   utilisée par quatre faces : deux morceaux de chaque côté.
3. **Classification** de chaque face par rapport à l'autre opérande, par un
   point intérieur (PR 5). Pour un compound, une face est dedans si elle est
   dans l'un de ses solides.
4. **Sélection.**

   | Opération | Faces gardées |
   |---|---|
   | fusion | A dehors de B, B dehors de A |
   | intersection | A dans B, B dans A |
   | différence | A dehors de B, B dans A **retournées** |

   Chaque arête de section garde exactement deux faces, une de chaque côté.
5. **Assemblage, sans sewing.**
   - Les faces gardées partagent déjà leurs arêtes : les shells sont leurs
     composantes connexes par arêtes, et chacun doit être fermé.
   - Les shells de volume positif sont des solides. Ceux de volume négatif
     sont des cavités, données au plus petit solide dont la boîte les
     contient.
   - `make_solid` vérifie et oriente les shells.
6. **Nettoyage.** Les faces non gardées, les wires de section et les shells et
   solides des copies sont effacés. Les arêtes et sommets inutilisés restent
   dans l'arène, comme après toute édition, jusqu'à un `compact()`.

### 3.3 Configurations refusées (`BuildErrc::NonGeneric`)

| Constat | Exemple |
|---|---|
| faces confondues (`coincident`) | deux boîtes posées sur le même plan |
| arête ou sommet d'un opérande sur une face de l'autre (`boundary_edges`, `points`) | un sommet de boîte exactement sur la sphère |
| faces tangentes le long d'une arête de section | cylindre posé sur un plan |
| intersection singulière ou incomplète | Steinmetz, Viviani |
| coupe pendante, face posée sur le bord de l'autre opérande | — |

C'est la recommandation de la question 5 : une erreur explicite plutôt qu'un
résultat faux. En testant, un premier cas sphère/boîte a ainsi été refusé
à raison : un sommet de la boîte était exactement sur la sphère.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| coplanaires et coaxiaux confondus dans cette PR | PR 9 b | découpe 2D commune et fusion d'arêtes : une PR à part entière |
| — | opérandes copiés | entrées intactes, retour arrière simple en cas d'échec |
| — | `BuildErrc::NonGeneric` ajouté à la fin de l'énumération | les codes existants ne changent pas |

## 5. Tests (`tests_brep_boolean`, 3 cas, 45 assertions, 2,4 s)

`check()` passe sur tous les résultats. Les volumes sont mesurés par
`signed_volume`, à 2e-3 près en relatif : c'est la précision de cet
échantillonnage sur des faces découpées.

| Test | Couvre |
|---|---|
| `boxes` | boîtes 0..2 et 1..3 : fusion (15, 12 faces), intersection (1, 6 faces), différence (7, 9 faces) ; opérandes intacts et valides ; boîtes disjointes (deux solides, intersection vide) ; boîte intérieure retranchée : **cavité** |
| `cylinders_and_spheres` | boîte 0..4 en B-spline et sphère R 3 centrée sur un coin, axe incliné (**cheminement**) : intersection (le huitième de la sphère), différence, fusion ; cylindre R 5 moins une boîte intérieure : cavité |
| `non_generic` | boîtes sur le même plan : `NonGeneric`, opérandes intacts ; formes mortes ou sans solide |

Les suites de découpe, de section, de lecture et d'écriture STEP et de
contrôle passent.

## 6. Points à valider

1. **Report des coplanaires et coaxiaux** en PR 9 b.
2. **Opérandes copiés**, jamais modifiés.
3. **Découpe des deux côtés par les mêmes arêtes**, assemblage sans sewing.
4. **Refus explicite** (`NonGeneric`) des configurations dégénérées.
5. **Cavités** données au plus petit solide dont la boîte les contient.
