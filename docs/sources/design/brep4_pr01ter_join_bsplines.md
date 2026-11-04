# Palier 4 — PR 1 ter : faces sur des B-splines voisines réunies

| | |
|---|---|
| PR | branche `feat/join-bspline-faces` |
| Document de référence | [brep_palier4.md](brep_palier4.md), § 0 ter ; [brep4_pr01_unify.md](brep4_pr01_unify.md), § 4 |
| Fichiers | `gbs/bssjoin.h` (nouveau), `gbs/bsctools.h`, `gbs-brep/unify.h`, `python/gbsBindSurfaceTools.cpp`, `python/gbsBindBrepOps.cpp`, `tests/tests_bssjoin.cpp` (nouveau), `tests/tests_brep_unify.cpp`, `python/tests/test_BSSurface.py`, `python/tests/test_brep_ops.py`, `News.md`, cette note |
| Statut | À valider (§ 6) ; revue de code faite (§ 7) |

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

- **Le côté joint.** Un côté réduit à un point (le pôle d'une sphère, l'apex
  d'un cône ou d'un dôme : tous ses pôles à tol du premier) n'est pas un côté à
  joindre. Parmi les couples de côtés qui se joignent, la fonction garde le plus
  lisse : le plus de nœuds retirés, puis le plus petit écart. Deux surfaces qui
  partagent deux côtés (les moitiés d'une bande fermée) se joignent donc le long
  du côté lisse, quel que soit l'ordre de recherche. `join_surfaces(a, b,
  side_a, side_b, tol)` joint le long des côtés donnés (`SurfaceSide` :
  `u_min`, `u_max`, `v_min`, `v_max`) ; le résultat dit lesquels ont été joints
  (`side_a`, `side_b`). Les coins des deux côtés sont comparés sur les surfaces
  d'origine avant toute copie : c'est le refus le plus fréquent.
- **Les paramètres de a sont gardés** : le point (u, v) de a est le point (u, v)
  de la surface jointe, à l'écart rapporté près, et les bornes de la surface
  jointe sont celles de a, au bit près, hors du côté commun. Les inversions sont
  faites par u → −u, exactes, et défaites de même. b est placé au-delà du côté
  commun ; ses paramètres sont envoyés sur ceux de la surface jointe par une
  carte affine (échange de u et v, inversions, échelle, décalage), rendue avec
  le résultat. C'est ce qui permet à `unify_faces` de garder les pcurves de a
  et de transformer exactement celles de b. Les bornes de b tombent sur celles
  de la surface jointe à l'arrondi de la carte près, quelques ulps d'un côté ou
  de l'autre : `unify_faces` ramène sur les bornes un pôle de pcurve que cet
  arrondi en fait sortir (§ 2.2).
- **Le long du côté**, les degrés sont portés au plus grand et chaque surface
  reçoit les nœuds de l'autre (élévation de degré et insertion de nœuds,
  exactes), après que le paramètre de b a été envoyé sur celui de a.
- **En travers du côté**, les degrés sont portés au plus grand, et b est mis à
  l'échelle pour que la longueur de sa dérivée première en travers rejoigne
  celle de a (moindres carrés le long du côté). Sans cela, la jonction ne peut
  pas être plus que C0. La rangée commune de pôles devient la moyenne des deux.
- **Le nœud de la jonction**, de multiplicité p (C0), est retiré autant de fois
  que toutes les rangées le permettent, le même nombre de fois à chaque rangée,
  comme le fait déjà `reduce_degree_rows`. On le retire une fois pour une
  jonction C1, et p fois quand les deux surfaces sont des morceaux d'un même
  polynôme. `knots_removed` donne ce nombre.
- **`deviation`** borne la distance entre la surface jointe et a ou b aux mêmes
  paramètres. C'est la somme de deux termes :
  - la moitié de l'écart des rangées communes ;
  - le déplacement des pôles par le retrait, mesuré en réinsérant le nœud
    (opération exacte) et en comparant aux pôles d'avant, puisqu'une B-spline
    reste dans l'enveloppe de ses pôles.

  **`deviation` ≤ tol.** L'algorithme A5.8 du *NURBS Book* teste chaque passe
  contre tol, mais les passes s'additionnent. Il ne sert donc qu'à dire combien
  de retraits essayer, avec ce que la moyenne des rangées laisse de tol
  (tol − écart / 2). Le retrait est ensuite mesuré comme ci-dessus, et repris
  avec un retrait de moins tant que la mesure dépasse ce reste.

La fonction ne traite que des surfaces non rationnelles, serrées à toutes leurs
bornes, dans les deux sens (l'élévation de degré et les échantillons le
demandent), et dont les directions extérieures concordent : la carte de b doit
avoir un déterminant positif, pour que la normale de b soit celle de la surface
jointe. Dans tout autre cas, elle ne rend rien.

`join_curves(a, b, tol)` fait de même pour deux courbes bout à bout, serrées,
b partant de la fin de a à tol près. C'est `join(a, b)` de `bsctools.h`, suivi
du retrait du nœud de la jonction tant que tol le permet (même mesure, même
reste). b est décalé pour partir de la fin de a, à sa propre vitesse : les
paramètres des deux courbes sont gardés, si bien qu'une arête et ses pcurves
restent paramétrées de la même façon. Le pôle commun devient leur moyenne.
`join` décale désormais les nœuds de b en `fin de a + (k − début de b)` : la
borne de fin est celle que calcule le lecteur STEP pour une courbe à travers sa
fermeture (`k + fin de a − début de b`, arrondi deux fois, pouvait s'en écarter
d'un ulp).

Le `join(srf1, srf2)` de `bsstools.h` existait déjà, mais il ne joint que le
long de u, sans vérifier le côté commun ni retirer le nœud, et il lève « Not
implemented » dès que les nœuds en v diffèrent. Il reste tel quel.
`join_surfaces` est une nouvelle fonction, pas une modification de `join` :
les appelants de `join` (l'extension d'une surface jusqu'à une courbe) en
dépendent.

### 2.2 `unify_faces` (`gbs-brep/unify.h`)

- **Faces réunies sur la jonction** (`merge_joined` dans la boucle) : quand
  `merge_pair` échoue et que les deux faces sont sur deux `BSSurface`
  différentes, de même orientation dans la coque, chaque arête qu'elles
  partagent dit le côté de chaque surface qu'elle longe (ses pcurves sur une
  borne). Pour chacun de ces couples de côtés, si les surfaces s'y joignent
  **au moins C1** (`knots_removed ≥ 1`) :
  - les deux faces passent sur la surface jointe. a garde ses pcurves, celles de
    b sont transformées exactement par la carte affine, puis ramenées sur les
    bornes là où l'arrondi de la carte les en fait sortir (16 ulps au plus). Les
    tolérances des faces, de leurs arêtes et de leurs sommets sont relevées à
    leur distance à la surface jointe, comme pour `to_plane` ;
  - `merge_pair` réunit les deux faces, désormais sur la même surface ;
  - s'il refuse (une arête partagée qu'une autre face utilise, une surface
    jointe qui se referme sur elle-même...), tout est défait : les anciennes
    faces reprennent leur place dans les coques, les nouvelles sont effacées,
    les tolérances reprennent leurs valeurs. `merge_pair` n'écrit rien quand il
    échoue.

  Deux surfaces jointes seulement C0, c'est-à-dire un pli, restent deux faces :
  réunies, elles donneraient une face qui porte une arête vive en son milieu.
- **`to_joined(a, b)`** passe seulement les deux faces sur la surface jointe,
  le long du côté de leur arête commune, ou, sans arête commune, du côté de la
  jointure la plus lisse. `put_on_plane` et lui finissent par la même aide,
  `rehost` : nouvelle face sur la surface, tolérances tenues, nom, place dans
  les coques.
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
  touchaient. Les tables (arêtes de chaque sommet, contours de chaque arête)
  sont faites une fois et tenues à jour à chaque recollement ; seuls les
  sommets aux bouts de l'arête recollée sont essayés de nouveau. Un morceau
  parcouru à rebours est inversé par `BSCurve::reverse`, qui garde ses bouts
  exacts.
- `join_edges`, `join_seams` puis `join_smooth_edges`, dans cet ordre : trois
  instructions, l'ordre d'évaluation d'une somme n'étant pas fixé en C++.
- `UnifyReport::surfaces_joined` donne le nombre de paires de faces réunies sur
  la jonction de leurs surfaces (comptées aussi dans `faces_merged`). Comme
  pour toute fusion, l'historique reçoit une entrée {face réunie, a, b}.
- Python : `join_surfaces` (les deux formes), `SurfaceSide`, `SurfaceJoin3d`
  (avec `side_a`, `side_b`), `join_curves`, `CurveJoin2d`, `CurveJoin3d`,
  `UnifyReport.surfaces_joined`.

## 3. Ce qui réutilise gbs

Rien n'est réimplémenté :

- l'élévation de degré (`increase_degree`) ;
- l'insertion (`insert_knots`) et le retrait de nœuds (`remove_knot`, avec et
  sans tolérance) ;
- la multiplicité (`multiplicity`, `knots_and_mults`) ;
- `BSSurface`, `BSCurve`, `CurveOnSurface` ;
- `join` de `bsctools.h` pour `join_curves`, et `BSCurve::reverse` ;
- dans `unify.h`, le reste de la fusion : `merge_pair`, `mapped`,
  `replace_in_shells`, `hold_tolerance`.

La mesure des tolérances de `put_on_plane` devient une aide,
`raise_to_surface`, partagée avec `put_on_joined`, `merge_joined` et
`join_smooth_edges`, et la fin de `put_on_plane` une autre, `rehost`.

## 4. Tests

`tests/tests_bssjoin.cpp` (11 cas) :

| Test | Contenu |
|---|---|
| `pieces_joined_back` | une surface bicubique coupée entre deux nœuds : nœud retiré 3 fois, la surface d'origine, ses nœuds, la carte identité ; coupée sur un nœud : retiré 2 fois ; le long de v ; le côté de a en u_min |
| `the_parameters_of_b_turned` | b aux paramètres échangés et inversés, ou inversés deux fois : joint, par la carte (déterminant positif) ; b inversé une seule fois (tourné de l'autre côté) : refusé |
| `other_knots_degrees_and_speed` | b d'un autre degré dans les deux sens, avec un nœud à lui en v, parcouru 4 fois plus lentement : degrés portés à 4, échelle 1/4 retrouvée, jonction C2 au moins |
| `a_kink_stays_a_knot` | deux morceaux plans à angle droit : joints C0, nœud gardé, écart nul |
| `sides_apart_or_close` | b écarté de 10⁻³ : refusé ; de 10⁻⁹ : écart borné au demi-écart |
| `large_parameters` | 120 tirages, paramètres de 100 à 1550, côté de a en u_min ou u_max, b tourné de trois façons : rien n'est évalué hors bornes, les bornes de a sont gardées au bit près, celles de b tombent à 4 ulps des autres |
| `not_clamped_refused` | a non serré le long du côté : refusé (l'élévation de degré lisait hors du tableau) ; une courbe non serrée à son bout éloigné : refusée |
| `collapsed_side_skipped` | un dôme coupé le long d'un méridien : joint le long du méridien, nœud retiré 2 fois, le dôme retrouvé ; le long de l'apex, demandé : refusé |
| `smoothest_side_chosen` | les moitiés d'une bande fermée, lisse d'un côté, en coin de l'autre (essayé en premier) : jointes du côté lisse ; du côté en coin quand on le demande (C0) ; des côtés qui ne sont pas une même courbe : refusés |
| `deviation_within_tol` | 400 courbes et 100 surfaces dont les pôles de b au-delà de la jonction sont déplacés de 0,6 tol au plus : `deviation` ≤ tol, et la distance mesurée ≤ `deviation` |
| `curves` | `join_curves` : les morceaux d'une courbe rejoints ; b ailleurs et d'un degré de plus, décalé à sa vitesse ; un coin C0 ; écartées : refusées ; `join` : la borne de fin calculée comme le lecteur STEP |

`tests/tests_brep_unify.cpp`, quatre cas de plus :

- `bspline_pieces_joined` : les 4 faces d'une grille 2 × 2 de morceaux d'une
  même B-spline, dont l'un a ses paramètres inversés, deviennent une face sur
  une surface, valide (`check`). Elle a 4 côtés (les arêtes recollées), garde
  sa tolérance et reproduit la surface d'origine à 10⁻¹² près ;
- `bspline_faces_at_a_kink` : deux faces B-splines à angle droit restent deux
  faces ;
- `bspline_edges_joined` : les arêtes des morceaux sur leurs isos B-splines
  exactes, comme après un import STEP : les morceaux du bas et du haut sont
  recollés sur la courbe jointe (la branche des courbes 3D de
  `join_smooth_edges`), et la face a quatre côtés, chacun une B-spline ;
- `bspline_faces_left_when_they_would_not_merge` : les moitiés de la bande
  fermée ci-dessus, cousues : jointes du côté lisse, elles donneraient une face
  qui fait le tour d'une surface fermée, ce que `merge_pair` refuse. Les faces
  restent telles qu'elles étaient : mêmes faces, mêmes surfaces, même place dans
  la coque, mêmes tolérances.

Python : `test_join_surfaces` (avec les côtés joints et la forme à côtés
donnés), `test_join_curves` (`test_BSSurface.py`), et des morceaux de B-spline
réunis par `unify_faces` (`test_brep_ops.py`).

Sur la branche (macOS arm64, clang 21 de conda-forge, Release, comparaison OCCT
comprise) :

- les 726 tests de `ctest` passent, booléens compris : `unify_faces` les suit
  tous ;
- `pytest python/tests` : 91 passent hors notebooks, avec et sans `CI`, et la
  limite connue reste en `xfail`.

## 5. Limites

| Limite | Suite possible |
|---|---|
| Surfaces rationnelles non traitées | même construction sur les pôles pondérés, quand un cas le demande |
| Un côté de a qui touche plusieurs faces voisines (jonction en T) : `sew` ne recoud pas une arête entière contre deux morceaux, si bien que les faces ne sont pas voisines au sens de `unify_faces` | couture en T dans `sew` (couper l'arête longue), utile aussi hors de cette PR |
| Une jonction seulement G1, c'est-à-dire des vitesses différentes le long du côté (une échelle qui varie), reste C0 dans la représentation : le nœud n'est pas retiré, et les faces restent deux | reparamétrage de b le long du côté, ou ré-approximation de la paire à tol |
| Les arêtes recollées sont jointes à leur propre vitesse (`join_curves` ne met pas à l'échelle) : deux arêtes de vitesses différentes ne se recollent pas | mettre à l'échelle la courbe et les pcurves de la seconde arête ensemble |
| Avec des paramètres de quelques centaines, l'évaluation d'une `CurveOnSurface` (les côtés des faces naturelles), `check` et `sew` lèvent encore : un point de pcurve sur une borne sort d'un ulp, au-delà de `knot_eps`, une tolérance absolue (2,2·10⁻¹⁴, moins d'un ulp au-delà de 128). Ce défaut est antérieur à cette PR et se voit sur master ; `to_joined` passe ces cas quand on le corrige (sonde : 299 jonctions sur 299) | un test de bornes relatif à l'échelle des paramètres, ou des pcurves ramenées sur les bornes, dans une PR à part (celle de `sew`) |

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

## 7. Revue de code

Quatorze remarques, toutes traitées dans cette PR.

| # | Remarque | Suite |
|---|---|---|
| 1 | Avec des paramètres de quelques centaines, `join_surfaces` lève (l'échantillon dépasse v1 d'un ulp), les bornes de a reviennent raccourcies (miroir `c − (c − k)`), les coins de b passent hors des bornes par la carte. | Échantillon borné à v1 ; inversions par u → −u, exactes, si bien que les bornes de a reviennent au bit près ; pcurves de b ramenées sur les bornes là où l'arrondi de la carte les en fait sortir ; `raise_to_surface` fait de même pour ses points (§ 2.1, § 2.2). Test `large_parameters`. Reste le défaut antérieur des `CurveOnSurface` (§ 5). |
| 2 | Un vecteur non serré le long du côté fait lire hors du tableau à l'élévation de degré ; `join_curves` ne teste que la jonction. | Surfaces serrées à toutes leurs bornes, courbes aux deux bouts, sinon rien (§ 2.1). Test `not_clamped_refused`. |
| 3 | Le premier couple de côtés qui concorde est rendu, même une jonction C0 par un apex ; deux côtés communs : peut-être pas celui de l'arête partagée. | Côtés réduits à un point écartés, la jointure la plus lisse gardée, forme à côtés donnés ; `unify_faces` joint le long du côté de l'arête partagée (§ 2.1, § 2.2). Tests `collapsed_side_skipped`, `smoothest_side_chosen`. |
| 4 | `deviation` peut dépasser tol (jusqu'à 1,22 tol mesuré) : A5.8 teste chaque passe contre tol, les passes s'additionnent. | Retrait mesuré par réinsertion, contre tol − écart / 2, repris avec un retrait de moins tant qu'il dépasse (§ 2.1). Test `deviation_within_tol`. |
| 5 | Les faces passent sur la surface jointe avant de savoir si `merge_pair` les réunira. | Essai puis retour : faces, coques et tolérances remises comme avant quand `merge_pair` refuse (§ 2.2). Test `bspline_faces_left_when_they_would_not_merge`. |
| 6 | `join_edges(…) + join_seams(…) + join_smooth_edges(…)` : ordre d'évaluation non fixé. | Trois instructions. |
| 7 | `bspline_piece` inverse à la main (`u1 + u2 − x`). | `BSCurve::reverse`, qui garde les bouts exacts. |
| 8 | Pas de test de la branche des courbes 3D de `join_smooth_edges`, ni de `join_curves` en C++, ni des cas limites. | Tests `bspline_edges_joined`, `curves` et ceux des points 1 à 4. |
| 9 | `join_smooth_edges` refait toutes ses tables et reparcourt tous les contours après chaque recollement. | Tables faites une fois, tenues à jour ; seuls les bouts de l'arête recollée sont essayés de nouveau ; un seul passage par contour (§ 2.2). |
| 10 | 32 copies de grilles avant le test des coins. | Coins comparés sur les surfaces d'origine, a tourné une fois par côté (§ 2.1). |
| 11 | `join_curves` refait `join` de `bsctools.h`. | `join_curves` = `join` puis le retrait ; `join` décale les nœuds de b exactement (§ 2.1). |
| 12 | La fin de `put_on_joined` recopie celle de `put_on_plane`. | Aide commune `rehost` (§ 3). |
| 13 | Commentaire « un millionième », seuil 10⁻⁹. | Commentaire corrigé : les pcurves se rejoignent à l'arrondi près. |
| 14 | Boucles réécrivables en C++23. | `std::ranges::fill_n`, `std::views::zip` (`zip_transform` manque à la libc++ de conda-forge), un seul passage par contour. |
