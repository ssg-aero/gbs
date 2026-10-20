# Palier 3 — PR 1 : discrétisation des arêtes et contours des faces

| | |
|---|---|
| PR | branche `feat/brep-mesh-edges` |
| Document de référence | [brep_palier3_mesh.md](brep_palier3_mesh.md), § 2 et § 3 ; part avec ce document et le bilan du palier 2 bis ([brep_palier2bis.md](brep_palier2bis.md), § 12) |
| Fichiers | `gbs-brep/mesh.h` (nouveau, dans `gbs-brep/brep`), `gbs-brep/builders.h`, `tests/tests_brep_mesh.cpp` (nouveau), `docs/sources/design/brep_palier3_mesh.md` (nouveau), `docs/sources/design/brep_palier2bis.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Première brique du maillage : les points du maillage le long des arêtes, partagés
entre les faces, et le contour de chaque face dans son espace de paramètres, prêt
pour la triangulation (PR 2 et 3).

## 2. Contenu

| Élément | Rôle |
|---|---|
| `MeshOptions` | flèche, angle (0,5 rad), longueur maximale ; `resolved` fixe la flèche par défaut à 1e-3 de la diagonale de la forme |
| `curve_parameters(c, a, b, o)` | paramètres d'une courbe : ses bouts et ses nœuds intérieurs (B-spline, rationnelle ou non), chaque intervalle coupé en deux tant que, à ses quarts, la courbe s'écarte de la corde de plus que la flèche, que les tangentes aux bouts tournent de plus que l'angle, ou que la corde dépasse la longueur maximale |
| `MeshBuilder` | sommets du maillage : un par sommet du modèle, un par point intérieur d'arête ; chaque arête discrétisée **une fois** |
| `MeshBuilder::boundary(f)` | contours (u, v) de la face, l'extérieur d'abord : pcurves des co-arêtes aux paramètres de leur arête, dans l'ordre et le sens du contour, chaque point avec son sommet |
| `BuildErrc::MissingPCurve` | une co-arête de face sans pcurve |

**Couture** : elle apparaît deux fois dans le contour (u1 et u2), avec les mêmes
sommets. **Arête dégénérée** (pôle, sommet de cône) : pas de courbe 3D ; ses
points sont pris sur sa pcurve, autant que l'angle en demande pour un tour
complet autour du point, en proportion de la part du côté de la face qu'elle
couvre (son côté entier : un tour, 13 segments pour 0,5 rad) ; tous sont le
sommet de l'arête.

## 3. Tests (`tests_brep_mesh`)

| Test | Cas |
|---|---|
| `curve_parameters` | segment : ses deux bouts, ou des pas de 0,25 sous une longueur maximale de 0,3 ; cercle rationnel de rayon 2 aux flèches 1e-2, 1e-3, 1e-5 : écart à la corde (20 points par intervalle) sous la flèche, tangentes sous l'angle, nœuds gardés, nombre de points du même ordre que le minimum r (1 − cos(φ/2)) |
| `face_boundaries` | quadrilatère : 4 points ; cylindre analytique et NURBS : chaque point sur la surface, les deux bouts de la couture deux fois, de part et d'autre de la couture ; sphère : chaque pôle un seul sommet, 13 segments au moins |
| `shared_edges` | cube à faces B-spline, longueur maximale 0,5 : 16 points par face, chaque segment de contour utilisé par deux faces, 8 + 12 × 3 sommets |

Les 498 tests C++ passent.

## 4. Points à valider

1. Le test de la flèche aux quarts de chaque intervalle (une ondulation plus
   fine qu'un quart d'intervalle pourrait passer ; les nœuds en départ limitent
   le risque pour les B-splines).
2. La discrétisation d'une arête dégénérée par un tour complet proportionnel :
   provisoire, à revoir à la PR 3 avec la triangulation autour des pôles.
