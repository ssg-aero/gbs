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
| `curve_parameters(c, a, b, o)` | paramètres d'une courbe : les bouts de ses morceaux (`CurvePieces` : ses segments de Bézier, coupés à [a, b]), chaque morceau coupé en deux (de Casteljau) tant que ses pôles s'écartent de la corde de plus que la flèche (borne **certifiée** par l'enveloppe convexe), que son polygone de contrôle tourne de plus que l'angle, ou que la corde dépasse la longueur maximale ; une courbe qui n'est pas une B-spline est échantillonnée (8 points, marge sur les milieux). Une flèche ou un angle non positifs lèvent `std::invalid_argument` |
| `MeshBuilder` | sommets du maillage : un par sommet du modèle, un par point intérieur d'arête ; chaque arête discrétisée **une fois** ; construit sur une forme, il résout les options (`resolved`), sinon il exige une flèche et un angle positifs |
| `MeshBuilder::boundary(f)` | contours (u, v) de la face, l'extérieur d'abord : pcurves des co-arêtes aux paramètres de leur arête, dans l'ordre et le sens du contour, chaque point avec son sommet |
| `BuildErrc::MissingPCurve` | une co-arête de face sans pcurve |
| `BuildErrc::NotSameParameter` | une arête dont les pcurves n'ont pas son paramètre (`same_parameter` faux) : refusée, les points (u, v) étant pris aux paramètres de la courbe 3D |

**Couture** : elle apparaît deux fois dans le contour (u1 et u2), avec les mêmes
sommets. **Arête dégénérée** (pôle, sommet de cône) : pas de courbe 3D ; ses
points sont pris sur sa pcurve, autant que l'angle en demande pour un tour
complet autour du point, en proportion de la part du côté de la face qu'elle
couvre (son côté entier : un tour, 13 segments pour 0,5 rad) ; tous sont le
sommet de l'arête.

## 3. Tests (`tests_brep_mesh`)

| Test | Cas |
|---|---|
| `curve_parameters` | segment : ses deux bouts, ou des pas de 0,25 sous une longueur maximale de 0,3 ; cercle rationnel de rayon 2 aux flèches 1e-2, 1e-3, 1e-5 : écart à la corde (20 points par intervalle) sous la flèche, tangentes sous l'angle, nœuds gardés, nombre de points du même ordre que le minimum r (1 − cos(φ/2)) ; longueur bornée sur le cercle ; **onde** quintique qui croise sa corde aux quarts (l'angle hors jeu) : sous la flèche, alors que le test aux quarts laissait 1,05e-2 ; nœud C0 (pli) jamais à l'intérieur d'un intervalle ; flèche ou angle nuls refusés |
| `holes_and_options` | face carrée trouée d'un cercle (une arête fermée) : contour extérieur d'abord et direct, trou indirect, le sommet de l'arête fermée une fois ; flèche par défaut depuis la forme ; options non résolues refusées ; `MissingPCurve`, `NotSameParameter` |
| `face_boundaries` | quadrilatère : 4 points ; cylindre analytique et NURBS : chaque point sur la surface, les deux bouts de la couture deux fois, de part et d'autre de la couture ; sphère : chaque pôle un seul sommet, 13 segments au moins |
| `shared_edges` | cube à faces B-spline, longueur maximale 0,5 : 16 points par face, chaque segment de contour utilisé par deux faces, 8 + 12 × 3 sommets |

Les 499 tests C++ passent.

### 3.1 Revue de code

| Commentaire | Suite donnée |
|---|---|
| la flèche par défaut (0) n'était appliquée que par `resolved`, que rien n'appelait : options par défaut, chaque intervalle coupé jusqu'au pas minimal | `MeshBuilder` construit sur une forme résout les options ; sinon, flèche ou angle non positifs lèvent `std::invalid_argument` (aussi `curve_parameters`) ; testé |
| flèche testée aux quarts seulement | bornée par les pôles des morceaux de Bézier (`CurvePieces`, enveloppe convexe), coupés par de Casteljau ; échantillonnage pour les autres courbes ; testé par l'onde (la mutation aux quarts échoue) |
| `to_segment` et `inner_knots` refont `distance_to_segment` et `knots()` | supprimés : `CurvePieces` découpe aux nœuds, `distance_to_segment` est utilisé |
| `same_parameter` ignoré | refusé (`NotSameParameter`) ; testé |
| tables en `std::map` sur des indices denses ; boîte des pcurves calculée à chaque face | vecteurs indexés ; boîte calculée à la première arête dégénérée |
| cas non couverts | trou, arête fermée, `MissingPCurve`, `resolved`, nœud C0, longueur sur un cercle : `holes_and_options`, `curve_parameters` |

## 4. Points à valider

1. La flèche certifiée par les pôles des B-splines ; les autres courbes
   (`CurveOnSurface`, primitives) échantillonnées sur 8 points avec une marge.
2. La discrétisation d'une arête dégénérée par un tour complet proportionnel :
   provisoire, à revoir à la PR 3 avec la triangulation autour des pôles.
