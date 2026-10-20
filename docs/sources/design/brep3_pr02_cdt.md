# Palier 3 — PR 2 : triangulation de Delaunay contrainte

| | |
|---|---|
| PR | branche `feat/brep-mesh-cdt` |
| Document de référence | [brep_palier3_mesh.md](brep_palier3_mesh.md), § 4.3 |
| Fichiers | `gbs-mesh/cdt.h` (nouveau), `tests/tests_mesh_cdt.cpp` (nouveau), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

La triangulation de Delaunay contrainte du plan, générique (sans lien avec le
BREP), écrite dans le dépôt (question 1 du palier, décidée) : elle triangulera
les faces dans leur espace de paramètres à la PR 3.

## 2. `ConstrainedDelaunay<T>` (`gbs-mesh/cdt.h`)

| Méthode | Rôle |
|---|---|
| constructeur (boîte) | super-triangle 20 fois plus grand que la boîte des points |
| `insert(p)` | Bowyer-Watson : localisation par marche depuis le dernier triangle (balayage en secours) ; la cavité (triangles dont le cercle circonscrit contient p) ne traverse pas une arête contrainte, sauf celle sur laquelle p tombe, qui est coupée en deux arêtes contraintes ; un point à moins de `eps` (1e-12 de la boîte) d'un sommet est ce sommet |
| `constrain(a, b)` | arête déjà là : marquée ; un sommet sur le segment le coupe ; sinon les arêtes qui le croisent sont retournées (Sloan : un quadrilatère non convexe attend son tour), puis les arêtes créées sont ramenées à Delaunay où elles le peuvent ; des contraintes qui se croisent lèvent `std::runtime_error` |
| `remove_outside()` | parcours en largeur 0-1 depuis le super-triangle, un de plus à chaque arête contrainte traversée : sont gardés les triangles de parité impaire (trous compris) ; plus d'insertion ensuite |
| `triangles()`, `points()`, `constrained(a, b)`, `valid()` | triangles dans le sens direct ; cohérence des voisins, de l'orientation et des contraintes (pour les tests) |

Prédicats : orientation en double précision avec une marge (aire d'un triangle
de hauteur `eps`), cercle circonscrit en `long double`. Ce ne sont pas les
prédicats exacts de Shewchuk : les points d'un maillage de face sont bien
séparés à la tolérance près ; à revoir si un cas réel le demande.

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

Performance (Release) : 10 000 points en 0,03 s, 100 000 en 0,8 s, insérés en
ordre aléatoire.

Les 503 tests C++ passent.

## 5. Points à valider

1. Les prédicats à marge plutôt qu'exacts (§ 2).
2. Le super-triangle fini : sans contour contraint, il peut manquer une arête
   de l'enveloppe convexe (absorbée par les triangles du super-triangle). Le
   maillage des faces contraint toujours leur contour.
3. La recherche des arêtes qui croisent une contrainte, et des sommets sur son
   segment, parcourt tous les triangles et points : en pratique, l'arête de
   contour est presque toujours déjà là après l'insertion de ses points. À
   remplacer par une marche si la PR 3 le demande.
