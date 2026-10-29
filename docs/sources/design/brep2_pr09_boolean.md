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

## 7. Revue de code (a posteriori, lot 5 de la revue du noyau BREP)

Chaque commentaire a été reproduit par une sonde, sauf mention contraire ;
chaque test cité échoue sur la version précédente. Les commentaires 1 et 3
remplacent le § 3.2, point 5 (cavités données par les boîtes, `make_solid`), et
le commentaire 2 le § 3.1 (copies effacées en cas d'échec).

| Commentaire | Suite donnée |
|---|---|
| 1. Une cavité était donnée au plus petit solide dont la boîte contient la sienne. Un cube de 2 creusé d'une cavité de 1, entouré d'un anneau carré mince (3,2 × 3,2 × 1,2, paroi de 0,1, volume 1,488), fusionnés : la cavité allait à l'anneau. Le cube valait 8, l'anneau 0,488, et `check()` passait sur les deux | corrigé : un point d'une face de la cavité est classé par rapport à chaque solide, sans ses cavités, et la cavité va au plus petit qui le contient. Testé (`cavity_inside_a_ring`) |
| 2. Un échec tardif ne remettait pas le modèle en état : un shell ouvert, une cavité hors de tout solide ou un refus de `make_solid` laissaient les faces gardées et les shells du résultat. Un `BRepError` levé en cours de route (une classification dont tous les rayons sont plafonnés, dans la section, `same_side` ou les coïncidences partielles) sautait le nettoyage. Ces cas sont vus à la lecture : aucune des sondes ne les a déclenchés. Un opérande que `copy_shape` ne sait pas copier (une arête sur une courbe qui n'est ni primitive ni B-spline) faisait lever une exception à `boolean`, et, en second opérande, laissait la copie du premier dans le modèle | corrigé : les tables du modèle ne font que croître, et l'opération retient leur taille au départ. En cas d'échec, tout ce qui a été ajouté depuis est effacé (`rollback`), aussi quand une étape lève `BRepError`, rendu comme `NonGeneric`. Les opérandes ne sont jamais touchés, puisque l'opération travaille sur leurs copies. Testé (`failure_leaves_the_model` : l'opérande non copiable en premier et en second) |
| 3. Coût : le volume de chaque shell du résultat était calculé deux fois, pour son signe puis par `make_solid`. Sur la fusion de deux boîtes découpées en 880 faces (1 472 faces au résultat), ce double calcul faisait 95 % des 6,5 s | corrigé : les solides sont faits directement des shells, déjà orientés par ce signe (vers l'extérieur pour un solide, vers la cavité pour une cavité). La fermeture, l'orientabilité et un volume non nul sont vérifiés comme dans `make_solid`. La fusion passe à 4,3 s. `face_volume_integral` reste en O(n³) par face (test du polygone par cellule) : lot 6 |
| 4. Coût : les coïncidences partielles classaient chaque échantillon par `classify` (boîtes recalculées, projection sur toutes les arêtes). Deux plaques superposées de n trous dont les trous se croisent : 1 trou, 0,17 s ; 4 trous, 1,0 s ; 16 trous, 7,5 s ; 36 trous, 31 s | corrigé avec `FaceClassifier` ([brep2_pr05_classify.md](brep2_pr05_classify.md) § 7) : 0,10 ; 0,33 ; 1,8 et 6,3 s |
| Non traité ici : des faces à quelques tolérances l'une de l'autre (1,5 et 2,5 tol) restent refusées (lot 4, [brep2_pr06_section.md](brep2_pr06_section.md) § 8). Les faces sont confondues à `tol` près, les sommets à 10 tol | à reprendre dans une PR de suite : le même rayon de confusion pour les faces, dans la section, `coincident.h` et `partial.h` |
| Vérifié | 1 920 booléens de contacts (boîte, cylindre, sphère et cône posés sur les faces, les arêtes et les coins de la boîte l'un de l'autre, avec ou sans recouvrement de 0,3), en fusion, intersection et différence : tous réussis après la correction de la régression de #198 ([brep2_pr07_split.md](brep2_pr07_split.md) § 7), aucun ne laisse le modèle changé |
| C++23 | deducing this pour la recherche de racine des shells, au lieu de `std::function` ; `views::iota \| ranges::to` pour l'ordre radial et les parents ; `append_range` pour les coupes et les sommets ; `ranges::contains` avec une projection pour les coupes en double |
