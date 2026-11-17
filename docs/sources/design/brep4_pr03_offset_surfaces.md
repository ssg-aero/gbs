# Palier 4 — PR 3 : surfaces décalées

| | |
|---|---|
| PR | branche `feat/brep-offset-surface` |
| Document de référence | [brep_palier4.md](brep_palier4.md), § 3.2 et § 7 |
| Fichiers | `gbs/offsets.h`, `gbs/curveoffset.h`, `gbs-brep/offset.h`, `python/gbsBindSurfaces.cpp`, `tests/tests_offsets.cpp` (nouveau), `tests/tests_first_derivatives.h` (nouveau), `tests/tests_brep_sweep.cpp`, `tests/tests_extrema.cpp`, `tests/tests_intersect_marching.cpp`, `tests/tests_brep_pipe.cpp`, `gbs-occt/tests/tests_surfbuild.cpp`, `docs/sources/examples/gbs_examples.py`, `docs/sources/examples/gallery.py`, `python/tests/test_primitives.py`, `python/tests/test_examples.py`, `docs/sources/design/brep_palier4.md`, `docs/sources/design/brep4_pr02_pipe.md`, `News.md`, cette note |
| Statut | À valider (§ 11) ; revue de code § 12 |

## 1. Objectif

Préparer le décalage des solides (PR 4) et les congés (PR 5), qui reposent sur
les surfaces décalées : une surface décalée à **dérivées exactes** (le
cheminement des intersections en a besoin, `brep2_pr08` notait la limite des
différences finies), le **décalage exact** quand il existe, et la **détection du
repli** (décalage au-delà d'un rayon de courbure).

Aussi : réécrire dans la notation de gbs les blocs de commentaires `//occt` de
`gbs/offsets.h`, qui décrivaient les dérivées des courbes décalées.

## 2. Sources

- Surfaces et courbes décalées, régularité et repli : T. Maekawa, *An overview
  of offset curves and surfaces*, Computer-Aided Design 31(3), 1999.
- Courbes planes décalées : R. T. Farouki, C. A. Neff, *Analytic properties of
  plane offset curves*, CAGD 7, 1990.
- Dérivées d'un quotient par récurrence (Leibniz), comme pour les B-splines
  rationnelles : L. Piegl, W. Tiller, *The NURBS Book*, 2e éd., § 4.5.
- Courbures principales : M. do Carmo, *Differential Geometry of Curves and
  Surfaces*, § 3.2 et 3.3.

OpenCascade ne sert que de comparaison en boîte noire dans les tests
(`Geom_OffsetSurface`).

## 3. Surface décalée (`gbs/offsets.h`)

`OffsetSurface<T>(S, d)` : P = S + d n, avec n = N / w, N = S_u × S_v et w = |N|.
Elle passe de `gbs::brep` à `gbs` (la géométrie n'a pas à dépendre du BREP) ;
`gbs::brep::OffsetSurface` reste un alias.

**Dérivées exactes, à tout ordre.** Avec X^(a,b) la dérivée d'ordre a en u et b
en v, la règle de Leibniz sur N = S_u × S_v, w² = N · N et N = w n donne, par
a + b croissant :

```
N^(a,b) = Σ C(a,i) C(b,j) S^(i+1,j) × S^(a-i,b-j+1)
w^(a,b) = ( (N·N)^(a,b) − Σ' C(a,i) C(b,j) w^(i,j) w^(a-i,b-j) ) / 2w     (Σ' sans (0,0) ni (a,b))
n^(a,b) = ( N^(a,b) − Σ_{(i,j)≠(0,0)} C(a,i) C(b,j) w^(i,j) n^(a-i,b-j) ) / w
P^(a,b) = S^(a,b) + d n^(a,b)
```

P^(a,b) ne demande que les dérivées de S jusqu'à (a + 1, b + 1) : chaque (a, b)
n'utilise que les ordres inférieurs dans les deux variables, et le rectangle
demandé se calcule seul, en une passe sur S (`Surface::derivatives`, évaluateur
en une passe des B-splines). `derivatives(u, v, 2, 2)` calcule ainsi le
rectangle 3 × 3, pas tout l'ordre 4, dans des tampons sur la pile. La même
récurrence, en une variable, sert aux courbes décalées (§ 7). Avant, P
n'avait que ses dérivées premières, par différences centrées, et
`has_second_derivatives()` était faux : le cheminement ne pouvait ni
reconnaître une tangence le long d'une courbe ni certifier un contact.

- **Base sans dérivées secondes** : dérivées premières par différences
  centrées, comme avant, et pas au-delà.
- **Pôle** (rangée de pôles confondus, sommet) : le point reste celui de la
  normale limite (`normal_direction`).
  - Le point est dégénéré quand S_u lu ou S_v lv (lu, lv les longueurs des
    bornes) est plus court que 10⁻¹⁰ de l'autre : à un pôle NURBS, S_u n'est
    que l'arrondi d'une dérivée nulle (6·10⁻¹⁷), que le test relatif à
    |S_u| |S_v| laissait passer.
  - Près d'un pôle, à une distance e des bornes, les dérivées d'ordre k de n
    perdent les chiffres de e^k : la récurrence divise par |N|, qui est en e. À
    10⁻⁶ d'un pôle de sphère NURBS, P_vv était faux de 0,45. Un point est donc
    pris pour proche d'une dégénérescence quand S_u ou S_v s'annulerait à moins
    de 10⁻⁴ des bornes (|S_u| ≤ 10⁻⁴ |S_uv| lv, ou de même pour S_v).
  - En un tel point, les dérivées sont extrapolées (cubique) depuis quatre
    points distants de 5·10⁻⁴ des bornes, sur le chemin du centre des bornes,
    ou d'un coin quand ce chemin reste près de la ligne dégénérée (un cône
    traversant son sommet au milieu de v) ; puis 5·10⁻³ et 5·10⁻². Aux pôles
    d'une sphère NURBS, à 5·10⁻⁸ près de leurs limites, au lieu de 2·10¹⁶.

## 4. Décalage exact (`offset_surface`)

`offset_surface(S, d, tol, boîte)` rend la surface décalée et le décalage de v
(`ExactOffset{surface, dv}` : S'(u, v + dv) = S(u, v) + d n) :

| Surface | Décalée |
|---|---|
| primitive | la primitive décalée (`primitive_offset`, qui passe aussi de `gbs-brep/offset.h` à `gbs/offsets.h`) : plan déplacé ; cylindre, sphère, cône d'un autre rayon (le cône décale son v) ; tore d'un autre petit rayon |
| `OffsetSurface` | le décalé de sa base à la somme des distances, si son propre décalage ne se replie pas (les normales coïncident alors) : exact pour une base primitive, le repli mesuré sur la base et ses nœuds |
| autre | `OffsetSurface` |

Elle rend `nullopt` là où le décalage se replie sur la boîte : une primitive
au-delà de son rayon de courbure, une autre surface là où `offset_margin` (§ 5)
n'est pas au-dessus de 10⁻⁹, le seuil du repli de `make_pipe`.

Le repli d'une primitive est mesuré sur la boîte (la surface rendue garde les
bornes de la primitive) :

- **Tore** : le petit rayon r + d reste positif, et le cercle des centres
  R + r cos v reste loin de l'axe sur la boîte, pour le tore comme pour son
  décalé. Au-delà de R − r vers l'extérieur, l'équateur intérieur traverse
  l'axe : un tore fuseau était rendu, sur un tore entier ; sur sa moitié
  extérieure, le décalage reste accepté.
- **Cône** : la boîte reste sur une nappe, et le rayon décalé reste positif aux
  deux bouts de la boîte décalée de dv. Le rayon de référence (en v = 0) n'est
  plus testé : il pouvait refuser un décalage valide loin du sommet, ou en
  accepter un qui traverse le sommet près de lui. Quand il est négatif, le
  repère du cône décalé glisse le long de l'axe jusqu'au début de la boîte,
  et dv en tient compte.

Les deux concordent avec `offset_margin` sur la NURBS de la même primitive.

## 5. Repli (`offset_margin`)

Avec k₁ et k₂ les courbures principales de S, comptées positives là où S se
courbe vers n (−1/R sur une sphère de normale extérieure),
P_u × P_v = (1 − d k₁)(1 − d k₂) N : le décalage est régulier tant que
1 − d k > 0 pour les deux, et se replie au-delà (Maekawa, § 3).

`offset_margin(S, d, boîte)` rend le plus petit 1 − d k sur la boîte, et où.

- **Courbures** : la deuxième forme fondamentale est écrite sur une base
  orthonormée du plan tangent (e₁ le long de S_u), où elle est symétrique. Ses
  valeurs propres sont (p + r)/2 ± hypot((p − r)/2, q), et la demi-différence
  est donc une somme de carrés. La forme habituelle H ± √(H² − K) perd la
  moitié des chiffres en un ombilic : sur une sphère NURBS, 10⁻⁸ au lieu de
  10⁻¹⁶.
- **Recherche** :
  - échantillonnage : 33 valeurs par direction, plus 4 par intervalle de nœuds
    d'une B-spline, ou de la base d'une surface décalée ;
  - affinage : recherche par motifs autour des 4 échantillons les plus bas,
    jusqu'à 10⁻⁶ d'un pas ;
  - les points à moins de 10⁻⁴ des bornes d'une dégénérescence (§ 3) sont
    sautés : les courbures y perdent leurs chiffres et sont continues.

## 6. Intersections de surfaces décalées

Aucun algorithme nouveau n'était nécessaire :

- **Entre primitives décalées** : le décalé d'une primitive est une
  primitive, et les formes fermées du palier 2 s'appliquent telles quelles.
- **Avec une surface décalée générique** : le cheminement s'applique. Il a
  maintenant les dérivées secondes, avec lesquelles il distingue une tangence
  le long d'une courbe d'un contact isolé et certifie les contacts (PR 1 bis).

Exemple : un plan tangent à un cylindre NURBS décalé le long d'une génératrice.
- Avant, avec les différences finies, le cheminement rendait deux points et un
  résultat déclaré singulier.
- Maintenant, il rend la génératrice, marquée tangente.

Deux sphères NURBS décalées qui se touchent donnent leur point de contact.

## 7. Courbes décalées (`gbs/curveoffset.h`)

Les blocs `//occt` de `gbs/offsets.h` sont supprimés. Les courbes décalées
C + f n suivent les mêmes récurrences, en une variable, avec
M = J C' dans le plan (J la rotation de +90°, n à gauche) et M = D × C' dans
l'espace (D une direction fixe) :

```
P^(k) = C^(k) + Σ_{i≤k} C(k,i) f^(i) n^(k−i)
```

L'en-tête de `curveoffset.h` les écrit dans la notation de gbs, avec les
sources.

- **`CurveOffset2D`** : la même valeur et les mêmes dérivées, à l'arrondi près.
  Les ordres au-delà de 2, qui levaient une exception, sont maintenant calculés.
- **Tangente nulle** (|M| sous 10⁻¹⁰ de |C'| et de |C''| fois la longueur des
  bornes) : à un bout de la courbe dont la tangente s'annule (pôles d'extrémité
  confondus), le point prend la normale de sa limite depuis l'intérieur (le
  long de C'' au début, de −C'' à la fin). Ailleurs, et pour les dérivées,
  `std::domain_error`, comme pour une tangente le long de la direction dans
  l'espace. Avant, un NaN, sans exception.
- **`CurveOffset3D` : correction.** Ses dérivées étaient fausses : la dérivée
  première oubliait celle de la normale, et l'écart à des différences finies
  était de l'ordre de 1 sur une B-spline. Elles sont maintenant exactes, à tout
  ordre.

## 8. BREP (`gbs-brep/offset.h`)

`offset_face` passe par `offset_surface`, sur la boîte des pcurves de la face
(`face_uv_box`) :

- **Face B-spline qui se replierait** : refusée (`InvalidBounds`, comme une
  primitive au-delà de son rayon de courbure ; avant, elle était acceptée).
- **Face sur une surface décalée** : sa décalée est prise sur la même base,
  sans empiler les décalages.
- **Erreur** : le modèle reste inchangé.

## 9. Python

- `gbs.OffsetSurface(basis, distance)`, avec `basis()` et `distance()`.
- `gbs.offset_surface(surface, distance, tol, box)` rend `(surface, dv)`, ou
  `None` si le décalage se replie.
- `gbs.offset_margin(surface, distance, box)` rend `(marge, (u, v))`.

## 10. Galerie et tests

Section « Offset surfaces » de la galerie, cas dans `gbs_examples.OFFSETS` :

| Cas | Contenu |
|---|---|
| `offset_wavy` | une face B-spline bicubique et ses décalées à ± 0,4 |
| `offset_paraboloid` | un paraboloïde (courbure 1 à l'origine) décalé de 0,5 et 0,9 vers son centre de courbure ; 1,0 refusé |
| `offset_tangent_plane` | cylindre NURBS décalé et plan tangent : la génératrice de contact en rouge |

Tests :

- **`tests_offsets`** (nouveau) :
  - dérivées jusqu'à l'ordre 3 contre les primitives décalées exactes, sur 5
    primitives et 2 distances, à 10⁻¹¹ ;
  - B-splines polynomiales et rationnelles : chaque ordre est la dérivée du
    précédent, et l'évaluation en une passe est identique ;
  - cylindre NURBS décalé au rayon r + d ; sphère NURBS décalée : 1,5 S et ses
    dérivées jusqu'à l'ordre 2, aux pôles et à 10⁻⁹, 10⁻⁶ et 5·10⁻⁵ d'eux
    (10⁻⁷), plus loin (10⁻¹¹) ;
  - cône NURBS traversant son sommet au milieu de v : dérivées finies au
    sommet ;
  - repli des primitives sur leur boîte : tore entier ou à moitié, cône près et
    loin de son sommet, comme `offset_margin` sur leurs NURBS ;
  - repli : paraboloïde (marge 0,5 à l'origine pour d = 0,5, d = 1 refusé,
    accepté sur une boîte loin de l'origine), sphère NURBS (marge 1 + d/R à
    10⁻⁹) ;
  - décalages exacts et décalé d'un décalé (cylindre exact, décalé d'un décalé
    replié refusé) ;
  - intersections : cercle d'un cylindre décalé par un plan, courbe de deux
    B-splines décalées, plan tangent (génératrice tangente), contact de deux
    sphères décalées ;
  - courbes décalées, planes et gauches : cercle décalé contre le cercle
    exact à tout ordre, dérivées contre des différences finies ; tangente nulle
    à un bout (limite), dérivée là et tangente le long de la direction
    (exception).
- **`tests_brep_sweep.offset_of_a_bspline_face_folds`** :
  - face d'un paraboloïde : décalages de 1 et 1,2 refusés, modèle inchangé ;
  - décalé d'un décalé pris sur la base ;
  - −5 accepté.
- **Tests existants** : `tests_extrema.no_second_derivatives`,
  `tests_intersect_marching.tangency_limits` et
  `tests_brep_pipe.any_curve_of_the_path` prenaient `OffsetSurface` comme
  surface sans dérivées secondes, ce qu'elle n'est plus. Ils passent maintenant
  par `FirstDerivativesOnly` (`tests_first_derivatives.h`), qui masque les
  dérivées secondes de sa base.
- **`tests_surfacebuild.offset_surface_derivatives`** (OCCT) : le point et les
  dérivées jusqu'à l'ordre 3 de `Geom_OffsetSurface` sur une B-spline
  bicubique, à 10⁻¹⁰ près, pour d = 0,2 et −0,15.

## 11. Points à valider

1. **`OffsetSurface` et `primitive_offset` passent dans `gbs/offsets.h`** (la
   géométrie), avec un alias dans `gbs::brep`.
2. **Dérivées près d'un pôle** extrapolées (cubique, pas de 5·10⁻⁴ des
   bornes) dans les 10⁻⁴ des bornes autour d'une dégénérescence, plutôt
   qu'une limite exacte : à 10⁻⁷ près environ.
3. **Repli détecté par échantillonnage** puis recherche par motifs, au seuil
   10⁻⁹ (celui de `make_pipe`), sur la boîte des pcurves de la face : un repli
   plus étroit qu'un pas d'échantillonnage, loin des 4 échantillons les plus
   bas, peut échapper. Dans `offset_face`, il rend `InvalidBounds`, comme une
   primitive. Le décalage de solide (PR 4) le rendra `NonGeneric`, comme le
   prévoit le § 3.2 du plan.
4. **Décalé d'un décalé** : il est pris sur la base, à la somme des distances.
5. **`CurveOffset3D` corrigée.** Ses résultats changent (ses dérivées étaient
   fausses), et les deux courbes décalées ont maintenant leurs dérivées à tout
   ordre.
6. **Pas d'algorithme d'intersection propre aux surfaces décalées** : les
   formes fermées entre primitives et le cheminement suffisent, maintenant que
   les dérivées sont exactes (§ 6).
7. **Repli d'une primitive mesuré sur la boîte** (de la face pour
   `offset_face`) ; un cône dont la boîte traverse le sommet (deux nappes) est
   refusé, comme à la lecture STEP.

## 12. Revue de code

Revue de la PR (11 points), dix traités dans cette PR :

| # | Point | Traitement |
|---|---|---|
| 1 | Pôle non détecté : le seuil relatif à \|S_u\| \|S_v\| laissait passer un S_u d'arrondi ; aux pôles d'une sphère NURBS, \|P_vv\| = 2·10¹⁶ au lieu de 1,3, P_v faux de 10 à 20 % | `regular` sur S_u lu et S_v lv ; près d'un pôle, extrapolation (§ 3) |
| 2 | Tore décalé au-delà de R − r : un tore fuseau accepté | cercle des centres loin de l'axe sur la boîte (§ 4) |
| 3 | Repli du cône testé en v = 0, pas sur la boîte : accepté près du sommet, refusé loin de lui | rayon décalé aux deux bouts de la boîte, repère glissé le long de l'axe (§ 4) |
| 4 | Pôles de la sphère NURBS testés seulement en finitude | comparés à 1,5 S, aux pôles et près d'eux |
| 5 | Repli bloqué sur une rangée dégénérée passant au centre | points vers les coins (§ 3) ; test du cône traversant son sommet |
| 6 | Décalé d'un décalé : base primitive non exacte, repli mesuré sur P sans les nœuds de la base | `offset_surface` de la base à la somme des distances, après celui de la base à la première |
| 7 | `derivatives(2, 2)` calculait tout l'ordre 4 (1,2 µs) | rectangle demandé, tampons sur la pile (0,83 µs) |
| 8 | Courbe décalée : une tangente nulle levait une exception, non documentée | limite à un bout ; exception documentée ailleurs (§ 7) |
| 9 | Convention de signe de k inversée dans l'en-tête de `offsets.h` | corrigée (k positif là où S se courbe vers n) |
| 10 | Récurrence de la normale unitaire écrite deux fois | une seule, sur un rectangle (courbes : B = 0) |
| 11 | Lecture STEP : OFFSET_SURFACE d'une base non élémentaire refusé, code du décalage exact dupliqué | laissé à une PR de `gbs-io` (nouvelle définition de surface, écriture, fichiers de test) |
