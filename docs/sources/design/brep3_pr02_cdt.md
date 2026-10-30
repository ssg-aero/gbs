# Palier 3 — PR 2 : triangulation de Delaunay contrainte

| | |
|---|---|
| PR | branche `feat/brep-mesh-cdt` |
| Document de référence | [brep_palier3_mesh.md](brep_palier3_mesh.md), § 4.3 |
| Fichiers | `gbs-mesh/cdt.h` (nouveau), `inc/topology/robustPredicates.h`, `tests/tests_mesh_cdt.cpp` (nouveau), `docs/sources/design/brep_palier3_mesh.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

La triangulation de Delaunay contrainte du plan, générique (sans lien avec le
BREP), écrite dans le dépôt (question 1 du palier, décidée) : elle triangulera
les faces dans leur espace de paramètres à la PR 3.

## 2. `ConstrainedDelaunay<T>` (`gbs-mesh/cdt.h`)

| Méthode | Rôle |
|---|---|
| constructeur (boîte) | super-triangle 20 fois plus grand que la boîte des points |
| `insert(p)` | Bowyer-Watson : localisation par marche depuis le dernier triangle ; un point hors du super-triangle ou non fini lève `std::invalid_argument` ; la cavité (triangles dont le cercle circonscrit contient p, marqués par une époque) ne traverse pas une arête contrainte, sauf celle sur laquelle p tombe, qui est coupée en deux arêtes contraintes ; l'éventail lié en une passe ; un point à moins de `eps` (1 000 epsilon de T fois la taille de la boîte) d'un sommet est ce sommet |
| `constrain(a, b)` | arête déjà là : marquée ; sinon **marche** depuis a (Sloan) : le triangle de l'éventail de a où entre le segment, puis à travers chaque arête croisée ; un sommet sur le segment le coupe ; les arêtes croisées sont retournées (un quadrilatère non convexe attend son tour), puis les arêtes créées ramenées à Delaunay où elles le peuvent ; une contrainte croisée lève `std::runtime_error` |
| `remove_outside()` | parcours en largeur 0-1 depuis le super-triangle, un de plus à chaque arête contrainte traversée : sont gardés les triangles de parité impaire (trous compris) ; plus d'insertion ensuite |
| `triangles()`, `points()`, `constrained(a, b)`, `valid()` | triangles dans le sens direct ; cohérence des voisins, de l'orientation et des contraintes (pour les tests) |

Prédicats **exacts** de Shewchuk : `orient2d` (déjà dans
`inc/topology/robustPredicates.h`) et `incircle`, ajouté à côté (filtre en
virgule flottante, puis calcul exact par expansions pour les cas presque
dégénérés, comme des points cocycliques). Les décisions sont les mêmes sur
toutes les plateformes.

La Delaunay en demi-arêtes de `inc/topology` n'est pas reprise (voir
[brep_palier3_mesh.md](brep_palier3_mesh.md), § 4.3).

## 3. Défaut trouvé en écrivant les tests

La récupération d'une contrainte bouclait dès vingt points : une arête à
retourner était lue dans l'ordre de la file, alors que le triangle peut la voir
dans l'autre sens, et le quadrilatère testé pour la convexité était faux. Les
extrémités sont lues sur le triangle.

## 4. Tests (`tests_mesh_cdt`)

| Test | Cas |
|---|---|
| `delaunay_of_random_points` | 500 points : structure valide, aucun point dans un cercle circonscrit, un point redonné est le même ; dans le carré aux côtés contraints : 2n − 2 − 4 triangles, aire 1, toujours Delaunay |
| `grid` | 10 × 10 points (cocirculaires partout) : 162 triangles, aire 81 |
| `constraints` | segment à travers 300 points : récupéré, structure valide ; diagonale d'une grille passant par des sommets : coupée en quatre ; point inséré sur une contrainte : la coupe ; contraintes croisées : refusées |
| `domains` | L concave : aire 3, 4 triangles, insertion refusée ensuite ; carré troué aux contours raffinés, points intérieurs : aire 12, aucun triangle dans le trou |
| `inner_constraints_and_nesting` | carré avec un trou, une île dans le trou (parité 3) et une ligne contrainte qui ne borde rien ; 400 points insérés après les contraintes : aire 68, triangles des deux côtés de la ligne |
| `errors_and_float` | point hors du super-triangle, point non fini : refusés ; `float` : un point sur une arête contrainte la coupe |
| `exact_incircle` | points cocycliques : exactement 0, signes justes à un ulp près, aussi décalés de 1e8 |

Performance (Release, ordre aléatoire) : 10 000 points en 0,03 s, 100 000 en
0,43 s (0,83 s avant la remise à zéro par époque), une contrainte à travers un
million de points en 1,7 ms. Un million de points en ordre aléatoire prend 42 s :
la marche depuis le dernier triangle parcourt de l'ordre de √n triangles par
point ; le maillage insère des points voisins, et un tri de Hilbert le règlerait
pour une insertion en lot.

Les 506 tests C++ passent.

### 4.1 Revue de code

| Commentaire | Suite donnée |
|---|---|
| `orient2d` exact déjà dans le dépôt, et une Delaunay en demi-arêtes que le document ne citait pas | `orient2d` utilisé ; le document du palier dit pourquoi la Delaunay en demi-arêtes n'est pas reprise (mesure : 6,8 s contre 0,05 s pour 20 000 points ; ni trous ni contraintes intérieures) |
| `incircle` en `long double`, double précision sur MSVC et Apple arm64 | `incircle` exact ajouté à `robustPredicates.h` ; testé sur des points cocycliques |
| `constrain` en O(n) par contrainte | marche de Sloan depuis le premier bout ; 1,7 ms à travers un million de points |
| vecteur remis à zéro à chaque insertion ; éventail lié en O(k²) | époques ; deux tables pour l'éventail |
| un côté du super-triangle pris pour l'intérieur | refusé (`std::invalid_argument`), comme un point non fini ; testé |
| `eps` fixe, sans effet en `float` | 1 000 epsilon de T fois la taille de la boîte, pour la fusion des points seulement : les décisions d'orientation sont exactes ; testé en `float` |
| cas non couverts | ligne contrainte isolée, île dans un trou, points après les contraintes, point hors de la boîte, `float` |

## 5. Points à valider

1. La structure indexée plutôt que les demi-arêtes de `inc/topology`.
2. Le super-triangle fini : sans contour contraint, il peut manquer une arête
   de l'enveloppe convexe (absorbée par les triangles du super-triangle). Le
   maillage des faces contraint toujours leur contour.
3. L'insertion d'un grand nombre de points en ordre quelconque (√n pas de
   marche par point) : un tri de Hilbert pour une insertion en lot, si un usage
   le demande.
