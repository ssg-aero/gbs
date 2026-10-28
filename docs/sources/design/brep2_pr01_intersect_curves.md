# Palier 2 — PR 1 : boîtes englobantes, intersection courbe-courbe

| | |
|---|---|
| PR | branche `feat/intersect-curves` |
| Document de référence | [brep_palier2.md](brep_palier2.md), § 2.1 à 2.3 (commité dans cette PR) |
| Fichiers | `gbs/intersect/box.h`, `gbs/intersect/curve_pieces.h`, `gbs/intersect/curves.h` (nouveaux), `tests/tests_intersect_curves.cpp` (nouveau), `docs/sources/design/brep_palier2.md` (nouveau), `docs/index.rst`, `News.md` |
| Statut | À valider, **avec le document de conception du palier 2** |

## 1. Objectif de la PR

Poser la première brique du palier 2 : trouver **toutes** les intersections de
deux courbes sur des plages bornées. La découpe des faces (PR 7) s'en servira
pour croiser, dans l'espace (u, v), les pcurves des arêtes de section avec
celles du bord. La classification d'un point dans une face (PR 5) s'en servira
aussi.

Le document de conception du palier 2 entre dans `master` avec cette PR :
il n'y a pas de PR de documentation seule.

## 2. API

```cpp
template <T> struct IntersectOptions { T tol = extrema_tol; T angular_tol = 1e-3; size_t max_pairs = 200000; };
template <T, dim> struct CurveCurveHit   { T t1, t2; point<T, dim> p; bool tangent; };
template <T>      struct CurveOverlap    { T t1a, t1b, t2a, t2b; };   // t2a > t2b : sens opposés
template <T, dim> struct CurveCurveResult { vector<CurveCurveHit> hits; vector<CurveOverlap> overlaps; bool complete; };

auto intersect(c1, a1, b1, c2, a2, b2, IntersectOptions = {}) -> CurveCurveResult;   // plages données
auto intersect(c1, c2, IntersectOptions = {}) -> CurveCurveResult;                   // bornes des courbes

template <T, dim> struct Box { min, max; add, inflate, overlaps(b, gap), contains(b, gap), diagonal };
```

## 3. Choix d'architecture et justification

### 3.1 Morceaux de courbe (`CurvePieces`)

| Courbe | Morceaux | Boîte | Moitiés |
|---|---|---|---|
| `BSCurve`, `BSCurveRational` | **segments de Bézier** (restriction à la plage, puis insertion des nœuds jusqu'à la multiplicité du degré) | boîte des pôles cartésiens : **garantie** (enveloppe convexe, poids positifs) | **de Casteljau** en coordonnées homogènes : exact |
| toute autre courbe (cercle générique, `CurveOnSurface`…) | un morceau échantillonné (8 intervalles) | boîte des échantillons, élargie de deux fois le plus grand écart de la courbe aux cordes : **heuristique** | nouvel échantillonnage |

Chaque morceau connaît sa **planéité** : la distance maximale de ses points
(pôles ou échantillons) à la corde de ses extrémités.

### 3.2 Algorithme

1. **Filtre.** Les paires de morceaux dont les boîtes sont plus éloignées que
   la tolérance sont écartées.
2. **Portion commune.** Si le plus petit morceau tient dans la boîte de
   l'autre, on projette cinq de ses points sur l'autre courbe. S'ils sont tous
   dans la tolérance, le morceau est une **graine de recouvrement** et n'est
   plus découpé. C'est ce qui évite l'explosion du nombre de paires le long
   d'une portion commune.
3. **Candidat.** Deux morceaux plats (écart à la corde ≤ 1e-3 × corde)
   donnent un point de départ, pris aux points les plus proches de leurs
   cordes. Il est affiné par **Gauss-Newton** sur |c1(s) − c2(t)|², sur les
   courbes entières, et gardé si la distance atteinte est sous la tolérance.
4. **Presque parallèles.** Deux morceaux plats presque parallèles (sinus de
   l'angle < 0,1) peuvent se croiser deux fois. Ils sont encore découpés
   jusqu'à 10 × la tolérance avant de donner un candidat.
5. **Découpe.** Sinon, on coupe en deux le plus grand des morceaux non plats.

Après le parcours :
- les graines qui se suivent sur les deux courbes forment les
  **recouvrements** ;
- les candidats sont triés et **dédoublonnés**. Deux solutions sont les
  mêmes si les courbes restent confondues entre elles, ce qui garde les deux
  solutions distinctes d'une courbe fermée ;
- les croisements aux bouts d'un recouvrement **l'étendent** ; ceux qu'il
  contient sont retirés ;
- une courbe fermée qui coupe l'autre en son point de départ ne compte qu'**une**
  solution.

### 3.3 Tangence ou recouvrement

C'est le seul point délicat, et il faut le valider. Une droite tangente à un
cercle de rayon R reste dans la tolérance sur une longueur d'environ
2√(2 R tol), soit 3 mm pour R = 1 m et tol = 1e-6 m. Par le seul critère de
distance, c'est une portion commune.

On les distingue par le **profil de distance** le long de la portion. Sur un
vrai recouvrement, la distance varie très peu. Sur une tangence, elle part de
presque zéro au contact et monte jusqu'à la tolérance aux bouts. Si l'écart
entre la plus grande et la plus petite distance dépasse tol / 4, la portion
est une tangence. Elle est alors remplacée par un croisement au point le plus
proche, marqué `tangent` (sinus de l'angle < `angular_tol`).

Ce critère classe aussi correctement deux courbes parallèles décalées de
tol / 2 : elles forment un recouvrement, et non une tangence.

### 3.4 Dimension quelconque

L'algorithme ne dépend pas de la dimension. Le document de conception parlait
de 2D (pcurves). En 3D, il donne les points où deux courbes se rencontrent dans
la tolérance (droites gauches : aucune solution). Cela servira à la section
pour fusionner les sommets des arêtes.

### 3.5 Dépendances

Aucune nouvelle. Pas de nlopt : Gauss-Newton 2 × 2 et la projection d'un point
par Newton sont écrits sur place. Les extrema existants (`extrema.h`) ne sont
pas touchés (question 11 du document de conception).

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| courbe-courbe en 2D | toute dimension | même code ; utile à la section en 3D |
| formule fermée pour les droites | pas de cas particulier | un segment de degré 1 est un Bézier plat : une seule paire, Newton converge en un pas |
| — | critère tangence / recouvrement (§ 3.3) | un recouvrement de √tol sinon pour toute tangence |
| — | `complete = false` au-delà de `max_pairs` | garde-fou contre une configuration pathologique |

## 5. Tests (`tests_intersect_curves`, 7 cas, 251 assertions, 50 ms)

| Test | Couvre |
|---|---|
| `boxes` | boîte vide, ajout, écart, inclusion |
| `segments` | croisement exact, parallèles, jonction en T à une extrémité, segments colinéaires dans le même sens et en sens opposé (recouvrement exact) |
| `circles_and_lines` | cercle NURBS et cercle générique (chemin échantillonné) : sécante, **tangente** (une solution marquée tangente), diamètre passant par le point de départ (deux solutions, pas trois) ; deux cercles ; arc posé sur le cercle (recouvrement) ; courbe avec elle-même |
| `bsplines` | sinusoïde cubique contre l'axe (5 solutions à kπ) et contre des horizontales ; deux sinusoïdes ; restriction à une demi-plage |
| `curves_3d` | droites concourantes et gauches, cercle incliné et une de ses cordes |
| `robustness` | cercle générique sur cercle NURBS (recouvrement sur toute la plage, à travers le point de départ) ; croisements rasants jusqu'à y = 0,9999 ; carré coupé par un cercle (8 solutions) |
| `limits` | plage restreinte, tolérance plus large, plafond de paires |

## 6. Points à valider

1. **Le document de conception du palier 2** ([brep_palier2.md](brep_palier2.md)),
   en particulier la question 1 (booléens) ; les recommandations sont
   appliquées tant qu'elles ne sont pas tranchées. Cette PR ne dépend que de
   la question 2 (emplacement `gbs/intersect/`).
2. **Bézier et de Casteljau** pour les B-splines, **échantillonnage** pour le
   reste, avec une boîte heuristique.
3. **Critère tangence / recouvrement** par le profil de distance (§ 3.3).
4. **Dimension quelconque** au lieu du 2D seul.
5. **Solutions aux extrémités** incluses, une seule solution au point de
   départ d'une courbe fermée.

## 7. Revue de code (a posteriori, lot 2 de la revue du noyau BREP)

Chaque commentaire a été reproduit par une sonde avant d'être corrigé. Les points
2, 3 et 5 remplacent le § 3.2 (étapes 2 à 4) et le critère du § 3.3 (écart de
distance de tol / 4). Les commentaires 8 et 9 portent sur du code des PR #191 et
#189 : `smooth_spans` et son emploi dans `deviation_based_params` sont entrés avec
#189 (`fix/area-properties-followup`), sous le master de #191.

| Commentaire | Suite donnée |
|---|---|
| 1. Graines fusionnées à travers la couture d'une c2 fermée (les points se suivent, t2 saute de b2 à a2) ; le test de tangence projetait alors sur le complément du recouvrement : le même cercle, coutures à 0 et 1 rad, donnait 9 croisements tangents et aucun recouvrement | Corrigé. Deux graines ne se suivent que si c2 va de l'une à l'autre sans faire le tour (milieu du chemin sur c2). Chaque portion est prolongée jusqu'au bout de la plage de c2 (la couture ferme l'une et ouvre l'autre) et mesurée sur sa propre plage de c2 : deux recouvrements jointifs. Test `review_seam_overlap` |
| 2. Courbes rasantes dans la tolérance : la graine de tangence partait d'un point où le gradient de Gauss-Newton s'annule, `refine` rendait n'importe quel point à moins de tol, le rayon de dédoublonnage (10 tol) était trop petit : y = 1 − 1e-6 contre le cercle unité, 24 croisements | Corrigé. `refine` va jusqu'à un minimum local : pas divisé par deux tant que la distance ne baisse pas, système 2 × 2 singulier par rapport à sa propre diagonale (des vitesses de 900 et 1 faisaient amortir un croisement rasant par Levenberg-Marquardt dans `minimize`), sortie d'un point selle (le sommet de la bosse entre deux croisements) le long de la courbure négative, des deux côtés. Le long d'une portion de tangence, les minima locaux de la distance (1D) sont les candidats. Deux croisements entre lesquels les courbes restent dans la tolérance ne font qu'un contact, tangent, au point le plus proche : sur une même portion de contact, ou dans la tolérance et parallèles en quelques points hors de toute grille, sans recouvrement entre eux ; les croisements d’un contact se suivent le long de c1 (autour de sa couture si elle est fermée). Les minima d'une graine et les croisements trouvés dedans sont donc un seul contact, ce qui remplace le retrait de la graine. y = 1 − 9e-7 … 1 + 5e-7 : un croisement ; y = 1 − 1e-6 (= tol, la limite) : deux. Test `review_grazing` |
| 3. Recouvrement pris pour une tangence dès que la distance varie de plus de tol / 4 : une droite et la même inclinée de 5e-7 (dans la tolérance sur toute sa longueur), un tronçon exact de 0,01 suivi d'un arc tangent | Corrigé. Chaque portion est prolongée jusqu'où d = tol (bissection) ou jusqu'au bout d'une plage. C'est un recouvrement, si elle dépasse dix tolérances, quand l'écart c1 − c2 n'est pas une quadratique du paramètre (résidu des moindres carrés > tol / 8 : partie exacte, approximation qui ondule) ou quand il reste, outre les départs tangents à ses bouts (distance au sommet de la parabole, d' / d''), un cœur de plus d'un vingtième. Sinon c'est une tangence (ou un croisement sous un petit angle). Un bout de plage où les courbes se quittent déjà en parabole, depuis un sommet intérieur, compte aussi comme un départ : un cercle qui touche le bord d'une face, restreint à la boîte de la face avec la tolérance, est une tangence (vu en combinant cette branche avec celle des formes fermées : la réunion d'un cylindre couché sur une boîte, couture sur le dessus, donnait une coque ouverte). Tests `review_overlap_within_tol`, `review_clipped_tangent_arc`. Limite : un contact d'ordre 3 (tangence en un point d'inflexion) est un recouvrement |
| 4. Croisement à la couture d'une courbe fermée rendu deux fois (`key(lo) == a` en égalité exacte ; le dédoublonnage prenait le milieu des paramètres, l'autre côté sur une courbe fermée) : 157 configurations fausses sur 800 | Corrigé : positions comparées, milieu pris par le chemin le plus court, autour de la couture, sur une courbe fermée (0 sur 800). Test `review_seam_crossing` |
| 5. Courbes presque parallèles juste au-delà de tol : `max_pairs` épuisé (paires en √(R / tol)), les croisements ailleurs perdus, `complete` ignoré des appelants | Corrigé. Test de séparation (*fat line*) le long de la normale à la corde de p. Surtout, deux morceaux plats presque parallèles ne sont plus coupés jusqu'à 10 tol : la distance de p à q est échantillonnée (17 points), Newton part de ses minima locaux et des bouts de q proches de p ; ils ne sont coupés que s'ils restent dans la tolérance sur deux échantillons de suite. Budget partagé entre les paires de morceaux initiaux, examinées à tour de rôle (256 paires). `complete` respecté : `imprint_pair` (`coincident.h`) échoue ; `boundary_crossings` (`section.h`, et son appel dans `partial.h`) et les croisements des `branch_spans` mettent `SectionResult::complete` à faux (booléen `NonGeneric`). Cercles concentriques à 1,1 à 3 tol jusqu’à R = 1000 : complets en moins de 0,02 s (incomplets dès R = 10, en 0,2 s). Respecter `complete` a montré que `coincident_faces` et `coincident_closed_bsplines` (`tests_brep_boolean`) passaient sur un résultat plafonné : deux interpolations d'un même cercle (257 et 129 pôles, écart 8e-7, tol 1,5e-8), 470 croisements et `complete = false` sur master ; désormais complet, 128 croisements, 0,07 s ; dans `coincident_closed_bsplines`, une B-spline de 202 pôles contre le cercle exact (écart 6e-8, tol 3,8e-9) : 689 croisements, incomplet, 0,06 s sur master ; 127 croisements et 4 recouvrements, complet, 0,18 s. Ces deux cas allongent `tests_brep_boolean` de 24,0 à 24,8 s Test `review_parallel_budget` |
| 6. `try_overlap` déclare une portion commune sur 5 échantillons : une courbe de degré 5 qui coupe l'axe exactement en ces 5 points donnait 1 croisement au lieu de 5 | Corrigé : les tangentes doivent aussi être parallèles aux échantillons (sinus ≤ 32 tol / taille du morceau). Test `review_false_overlap` |
| 7. Une B-spline perd ses segments de Bézier quand une borne est à moins de `knot_eps` d'un nœud sans y être (`knots()` fusionne, `bezier_poles` non) : `pieces(nextafter(0.5, 0), 1)` échantillonné | Corrigé : une borne à moins de `knot_eps` × max(1, \|a\|, \|b\|) d'un nœud est mise sur ce nœud avant `trim`, et les segments sont comptés comme `bezier_poles` (nœuds différents). `trim_begin` et `trim_end` (`knotsfunctions.ixx`) gardent leur `knot_eps` absolu, hors du lot. Test `review_pieces_ulp_bounds` |
| 8. `smooth_spans` compare à `knot_eps` absolu (2,2e-14) : un intervalle d'un ulp revient dès \|u\| ≥ 128 (`smooth_spans(b, 0, nextafter(500))`, un nœud en double dans le maillage d'une arête) | Corrigé : `knot_eps` × max(1, \|a\|, \|b\|), pour les bornes et entre ruptures voisines. Test `smooth_spans_scale` (`tests_curves`) |
| 9. `deviation_based_params` répète une rupture qui tombe sur l'un des n paramètres uniformes (nœud à 0,5 avec n = 3 ; composite de deux segments égaux : 0 1 1 2, `discretize` rend un point deux fois) | Corrigé : voisins à moins de la même tolérance relative retirés après le tri. Test `deviation_based_params_breaks_once` (`tests_bscanalysis`) |
