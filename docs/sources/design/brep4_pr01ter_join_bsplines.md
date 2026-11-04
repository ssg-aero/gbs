# Palier 4 — PR 1 ter : faces sur des B-splines voisines réunies

| | |
|---|---|
| PR | branche `feat/join-bspline-faces` |
| Document de référence | [brep_palier4.md](brep_palier4.md), § 0 ter ; [brep4_pr01_unify.md](brep4_pr01_unify.md), § 4 |
| Fichiers | `gbs/bssjoin.h` (nouveau), `gbs-brep/unify.h`, `python/gbsBindSurfaceTools.cpp`, `python/gbsBindBrepOps.cpp`, `tests/tests_bssjoin.cpp` (nouveau), `tests/tests_brep_unify.cpp`, `python/tests/test_BSSurface.py`, `python/tests/test_brep_ops.py`, `News.md`, cette note |
| Statut | À valider (§ 6) |

## 1. Objectif

La PR 1 réunit les faces d'un même domaine : même objet surface, primitives
égales, plans sous forme B-spline. Deux faces voisines sur **deux B-splines
distinctes** restaient deux faces, même quand leurs surfaces se prolongent
l'une l'autre sans défaut. C'est le cas :

- d'un modèle découpé en carreaux (un STEP écrit carreau par carreau, une
  surface coupée par `split_u` ou `split_v`) ;
- d'une forme construite par morceaux rectangulaires, chacun une B-spline à
  part, C2 avec ses voisins.

CATIA et SolidWorks montrent une face là où la surface est d'un seul tenant.
Cette PR ajoute à gbs la fonction générique qui manquait, **joindre deux
B-splines le long d'un côté commun**, et s'en sert dans `unify_faces`.

## 2. Contenu

### 2.1 `join_surfaces(a, b, tol)` (`gbs/bssjoin.h`)

Elle cherche un côté de a et un côté de b qui sont une même courbe à tol près
(u = u_min, u_max, v = v_min ou v_max de chacun, dans un sens ou dans l'autre),
puis construit une seule surface.

- **Les paramètres de a sont gardés** : le point (u, v) de a est le point (u, v)
  de la surface jointe, à l'écart rapporté près. b est placé au-delà du côté
  commun ; ses paramètres sont envoyés sur ceux de la surface jointe par une
  carte affine (échange de u et v, inversions, échelle, décalage), rendue avec
  le résultat. C'est ce qui permet à `unify_faces` de garder les pcurves de a
  et de transformer exactement celles de b.
- **Le long du côté**, les degrés sont portés au plus grand et chaque surface
  reçoit les nœuds de l'autre (élévation de degré et insertion de nœuds,
  exactes), après que le paramètre de b a été envoyé sur celui de a.
- **En travers du côté**, les degrés sont portés au plus grand, et b est mis à
  l'échelle pour que la longueur de sa dérivée première en travers rejoigne
  celle de a (moindres carrés le long du côté). Sans cela, la jonction ne peut
  pas être plus que C0. La rangée commune de pôles devient la moyenne des deux.
- **Le nœud de la jonction**, de multiplicité p (C0), est retiré autant de fois
  que toutes les rangées le permettent à tol près : l'algorithme A5.8 du *NURBS
  Book*, appliqué le même nombre de fois à chaque rangée, comme le fait déjà
  `reduce_degree_rows`. On le retire une fois pour une jonction C1, et p fois
  quand les deux surfaces sont des morceaux d'un même polynôme.
  `knots_removed` donne ce nombre.
- **`deviation`** borne la distance entre la surface jointe et a ou b aux mêmes
  paramètres. C'est la somme de deux termes :
  - la moitié de l'écart des rangées communes ;
  - le déplacement des pôles par le retrait, mesuré en réinsérant le nœud
    (opération exacte) et en comparant aux pôles d'avant, puisqu'une B-spline
    reste dans l'enveloppe de ses pôles.

La fonction ne traite que des surfaces non rationnelles, serrées sur le côté
commun, et dont les directions extérieures concordent : la carte de b doit
avoir un déterminant positif, pour que la normale de b soit celle de la surface
jointe. Dans tout autre cas, elle ne rend rien.

`join_curves(a, b, tol)` fait de même pour deux courbes bout à bout, b partant
de la fin de a à tol près. b est décalé pour partir de la fin de a, à sa propre
vitesse : les paramètres des deux courbes sont gardés, si bien qu'une arête et
ses pcurves restent paramétrées de la même façon. Le pôle commun devient leur
moyenne, et le nœud de la jonction est retiré tant que tol le permet.

Le `join(srf1, srf2)` de `bsstools.h` existait déjà, mais il ne joint que le
long de u, sans vérifier le côté commun ni retirer le nœud, et il lève « Not
implemented » dès que les nœuds en v diffèrent. Il reste tel quel.
`join_surfaces` est une nouvelle fonction, pas une modification de `join` :
les appelants de `join` (l'extension d'une surface jusqu'à une courbe) en
dépendent.

### 2.2 `unify_faces` (`gbs-brep/unify.h`)

- **`to_joined(a, b)`** (`put_on_joined` dans la boucle) : quand
  `merge_pair` échoue et que les deux faces sont sur deux `BSSurface`
  différentes qui se joignent **au moins C1** (`knots_removed ≥ 1`), les deux
  faces passent sur la surface jointe. a garde ses pcurves, celles de b sont
  transformées exactement par la carte affine. Les tolérances des faces, de
  leurs arêtes et de leurs sommets sont relevées à leur distance à la surface
  jointe, comme pour `to_plane`. `merge_pair` réunit ensuite les deux faces,
  désormais sur la même surface. Deux surfaces jointes seulement C0, c'est-à-dire
  un pli, restent deux faces : réunies, elles donneraient une face qui porte une
  arête vive en son milieu.
- **`join_smooth_edges`** recolle les arêtes qui se suivent sans coin. Il s'agit
  de deux arêtes bout à bout en un sommet qu'aucune autre arête n'utilise,
  consécutives dans chaque contour qui les utilise, et dont les pcurves
  B-splines se joignent (`join_curves`). L'arête recollée est :
  - sur la courbe jointe, si les deux courbes sont des B-splines jointes au moins
    C1 ;
  - sinon, la courbe de sa pcurve jointe sur la surface de la première face
    (`CurveOnSurface`, exacte sur celle-ci), si les pcurves se joignent au moins
    C1 dans chaque contour.

  Les côtés des faces naturelles (`make_face`) sont des `CurveOnSurface`, que
  `join_edges` ne recolle pas, puisqu'il demande la même courbe. Sans cette
  étape, une face réunie garderait ses côtés coupés là où les morceaux se
  touchaient.
- `UnifyReport::surfaces_joined` donne le nombre de paires de faces passées sur
  la jonction de leurs surfaces, et l'historique reçoit deux entrées par paire
  {nouvelle face, ancienne face, invalide}, comme `to_plane`.
- Python : `join_surfaces`, `SurfaceJoin3d`, `join_curves`, `CurveJoin2d`,
  `CurveJoin3d`, `UnifyReport.surfaces_joined`.

## 3. Ce qui réutilise gbs

Rien n'est réimplémenté :

- l'élévation de degré (`increase_degree`) ;
- l'insertion (`insert_knots`) et le retrait de nœuds (`remove_knot`, avec et
  sans tolérance) ;
- la multiplicité (`multiplicity`, `knots_and_mults`) ;
- `BSSurface`, `BSCurve`, `CurveOnSurface` ;
- dans `unify.h`, le reste de la fusion : `merge_pair`, `mapped`,
  `replace_in_shells`, `hold_tolerance`.

La mesure des tolérances de `put_on_plane` devient une aide,
`raise_to_surface`, partagée avec `put_on_joined` et `join_smooth_edges`.

## 4. Tests

`tests/tests_bssjoin.cpp` (5 cas) :

| Test | Contenu |
|---|---|
| `pieces_joined_back` | une surface bicubique coupée entre deux nœuds : nœud retiré 3 fois, la surface d'origine, ses nœuds, la carte identité ; coupée sur un nœud : retiré 2 fois ; le long de v ; le côté de a en u_min |
| `the_parameters_of_b_turned` | b aux paramètres échangés et inversés, ou inversés deux fois : joint, par la carte (déterminant positif) ; b inversé une seule fois (tourné de l'autre côté) : refusé |
| `other_knots_degrees_and_speed` | b d'un autre degré dans les deux sens, avec un nœud à lui en v, parcouru 4 fois plus lentement : degrés portés à 4, échelle 1/4 retrouvée, jonction C2 au moins |
| `a_kink_stays_a_knot` | deux morceaux plans à angle droit : joints C0, nœud gardé, écart nul |
| `sides_apart_or_close` | b écarté de 10⁻³ : refusé ; de 10⁻⁹ : écart borné au demi-écart |

`tests/tests_brep_unify.cpp`, deux cas de plus :

- `bspline_pieces_joined` : les 4 faces d'une grille 2 × 2 de morceaux d'une
  même B-spline, dont l'un a ses paramètres inversés, deviennent une face sur
  une surface, valide (`check`). Elle a 4 côtés (les arêtes recollées), garde
  sa tolérance et reproduit la surface d'origine à 10⁻¹² près ;
- `bspline_faces_at_a_kink` : deux faces B-splines à angle droit restent deux
  faces.

Python : `test_join_surfaces`, `test_join_curves` (`test_BSSurface.py`), et des
morceaux de B-spline réunis par `unify_faces` (`test_brep_ops.py`).

Sur la branche (macOS arm64, clang de conda-forge, Release, sans OCCT) :

- les 717 tests de `ctest` passent, booléens compris : `unify_faces` les suit
  tous, et aucun résultat de la galerie ne change ;
- les tests Python passent aussi : 93, et la limite connue reste en `xfail`.

## 5. Limites

| Limite | Suite possible |
|---|---|
| Surfaces rationnelles non traitées | même construction sur les pôles pondérés, quand un cas le demande |
| Un côté de a qui touche plusieurs faces voisines (jonction en T) : `sew` ne recoud pas une arête entière contre deux morceaux, si bien que les faces ne sont pas voisines au sens de `unify_faces` | couture en T dans `sew` (couper l'arête longue), utile aussi hors de cette PR |
| Une jonction seulement G1, c'est-à-dire des vitesses différentes le long du côté (une échelle qui varie), reste C0 dans la représentation : le nœud n'est pas retiré, et les faces restent deux | reparamétrage de b le long du côté, ou ré-approximation de la paire à tol |
| Les arêtes recollées sont jointes à leur propre vitesse (`join_curves` ne met pas à l'échelle) : deux arêtes de vitesses différentes ne se recollent pas | mettre à l'échelle la courbe et les pcurves de la seconde arête ensemble |

## 6. Points à valider

1. **Nouvelle fonction** `join_surfaces` plutôt qu'une modification de `join`
   (§ 2.1), dans un nouvel en-tête `gbs/bssjoin.h`.
2. **Seuil de réunion** : au moins C1 (un nœud retiré). Un pli joint C0 reste
   deux faces.
3. **Tolérance du retrait** : celle de `unify_faces`, les tolérances relevées
   à l'écart mesuré, comme pour `to_plane`. Pour des morceaux d'un même
   polynôme, l'écart est un arrondi, de l'ordre de 10⁻¹⁶.
4. **Arêtes recollées sur une `CurveOnSurface`** de la surface de la première
   face, quand leurs courbes ne sont pas des B-splines (§ 2.2).
