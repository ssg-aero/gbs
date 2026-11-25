# Évaluation des B-splines de degré élevé

| | |
|---|---|
| PR | branche `fix/bspline-high-degree` |
| Origine | revue de la PR 252 (`feat/brep-serialization`, 24 novembre 2026) : un fichier hostile de 2 Ko (degré 40) bloquait chaque évaluation pendant des heures. La PR 252 refuse les degrés au-delà de 24 en lecture et en écriture (`brep_json_max_degree`) ; la cause, dans `gbs/`, est traitée ici |
| Fichiers | `gbs/basisfunctions.ixx`, `gbs/modules/basisfunctions.cppm`, `gbs/knotsfunctions.ixx`, `gbs/bscinterp.h`, `gbs/bssinterp.h`, `gbs/bssapprox.h`, `gbs/bscshaping.h`, `gbs/bssshaping.h`, `gbs/bscurve.h`, `gbs/gbsconstants.h`, `gbs/gbslib.h`, `tests/tests_basisfunctions.cpp`, `tests/tests_bsinterp.cpp`, `tests/tests_curve_derivatives.cpp`, `tests/tests_curves_shapping.cpp`, `tests/tests_surfaces_shapping.cpp`, `python/tests/test_BSCurve.py`, `News.md` |
| Hors périmètre | la borne de `gbs-io/brep_json.h` (point 4) |
| Statut | À valider |

## 1. Constat

Les évaluateurs de `gbs/basisfunctions.ixx` calculent les p+1 fonctions de
base non nulles d'un intervalle de nœuds par l'algorithme A2.2 de *The NURBS
Book* (valeurs) ou A2.3 (dérivées), en O(p²), sur des tableaux de taille fixe
`bspline_stack_max_degree + 1` (24 + 1) sur la pile. Au-delà de 24, ils se
rabattaient sur `basis_function`, la récurrence de Cox-de Boor écrite
récursivement : chaque N(i,p) appelle deux N(·,p−1), soit 2^p feuilles. Le coût
double à chaque degré. `BSCurve::value` sur une courbe de quatre intervalles
(clang 21, -O3, Apple M) :

| Degré | master | branche |
|---|---|---|
| 24 | 0,24 µs | 0,24 µs |
| 25 | 0,41 s | 0,33 µs |
| 26 | 0,80 s | 0,35 µs |
| 30 | 12,8 s | 0,43 µs |
| 40 | ≈ 3,6 h (extrapolé, ×2 par degré) | 0,68 µs |

La valeur et les dérivées 1 à 3 prennent 1,0 µs par point au degré 26, 1,7 µs
au degré 40. La revue de la PR 252 avait mesuré 0,4 s au degré 26, 6,7 s à 30
et environ deux heures à 40 sur une autre courbe.

Le repli récursif était à huit endroits : `eval_value_decasteljau` et
`eval_ders_decasteljau` (courbe et surface), `eval_rational_ders` (courbe,
pour un ordre de dérivation au-delà de 24 : la récursion par ordre de
`eval_rational_value_simple`, 2^n évaluations), `fill_basis_row` et
`fill_basis_band` (`knotsfunctions.ixx`, matrices de collocation des
interpolations et approximations), `assemble_collocation_sparse`
(`bscinterp.h`), `build_poles` (`bssinterp.h`) et `approx` sous contraintes
(`bssapprox.h`). Une interpolation de degré 26 sur 45 points prenait des
minutes ; elle prend 0,2 ms.

## 2. Décision

Un seul algorithme, A2.2 et A2.3, pour tous les degrés. Sa mémoire de travail
est sur la pile jusqu'au degré 24, comme avant, et sur le tas au-delà. Le coût
est O(p²) par évaluation (O(p(p+n)) pour n dérivées), la mémoire O(p²) pour
A2.3 (la table `ndu`, (p+1)²).

L'algorithme est écrit une fois, dans `gbs::basis_detail`, sur un type de
mémoire de travail passé en paramètre de modèle (*policy*) :

```
template <typename T, size_t N> struct stack_scratch   // T v[N] : tableau local
template <typename T, size_t N> struct heap_scratch    // std::vector<T> v(n)

basis_detail::basis_funcs<Scratch>       A2.2
basis_detail::ders_basis_funs<Scratch>   A2.3
basis_detail::basis_ders<Scratch>        dérivée d'ordre d des p+1 fonctions
basis_detail::eval_value<Scratch>        corps de eval_value_decasteljau (courbe, surface)
basis_detail::ders_span<Scratch>         dérivées 0..n sur un intervalle (courbe, surface)
basis_detail::eval_rational_ders<Scratch>  corps de eval_rational_ders (courbe)
basis_detail::*_heap                     instanciations heap_scratch, GBS_NOINLINE
```

Les évaluateurs de dérivées gardent le prologue qu'ils avaient (sortie mise
à zéro, intervalle de nœuds trouvé), puis choisissent ; ceux de valeur
choisissent d'abord (§ 11, point 1) :

```
eval_ders_decasteljau(u, k, poles, p, n, CK)
 ├─ prologue : CK = 0, i = find_span(…)
 ├─ p ≤ 24 : basis_detail::ders_span<stack_scratch>   en ligne, aucune allocation
 └─ p > 24 : basis_detail::ders_span_heap              hors ligne, rend ses points
               └─ basis_detail::ders_span<heap_scratch>   (l'appelant les recopie)

eval_value_decasteljau(u, k, poles, p, d)
 ├─ p ≤ 24 : basis_detail::eval_value<stack_scratch>
 └─ p > 24 : basis_detail::eval_value_heap             hors ligne, rend sa valeur
```

Les entrées `*_heap` rendent leurs points par valeur (un `std::array` ou un
`std::vector`) et ne reçoivent jamais le pointeur de sortie de l'appelant ;
le § 11, point 1 dit pourquoi.

`basis_funcs`, `ders_basis_funs` et `basis_ders` (publiques) choisissent de
même. Les assemblages (`fill_basis_row`, `fill_basis_band`,
`assemble_collocation_sparse`, `build_poles`, `approx` sous contraintes, le
façonnage) ont un seul chemin : un tampon sur la pile jusqu'au degré 24, un
`std::vector` au-delà, puis la fonction publique, qui choisit elle-même.
`band_lu` (`bscinterp.h`) refuse toujours les degrés au-delà de 24 et passe la
main au chemin creux pivoté, qui a maintenant le même algorithme (point 3).

`eval_rational_ders` (surface) avait déjà une table sur le tas au-delà de
24 ordres ; elle n'a pas changé. `eval_rational_ders` (courbe) en a une
maintenant, au lieu de la récursion par ordre.

`GBS_NOINLINE` (`gbslib.h`) : `__declspec(noinline)` sous MSVC,
`[[gnu::noinline]]` ailleurs (clang, clang-cl, GCC).

## 3. Pourquoi cette forme

Deux formes plus simples ont été essayées et mesurées, puis écartées.

1. **Un tampon « pile ou tas » choisi à l'exécution** (un seul chemin de
   code, `basis_buffer<T, N>(n)` : tableau local si n ≤ N, sinon le tas). Les
   évaluateurs de bas degré sont devenus 1,3 à 2,5 fois plus lents (la valeur
   d'une courbe de degré 3 : 7,2 → 18 ns), et 5 % des valeurs ont changé
   dans leurs derniers bits. Sur arm64, le vectoriseur de boucles de LLVM
   traite la somme `pt[c] += N * pole[c]` comme une réduction ordonnée et
   sépare alors la multiplication-addition fusionnée (FMA) que le compilateur
   a le droit de former (`-ffp-contract=on`) ; il le fait ou non selon le code
   autour. Changer la forme de la fonction change ce choix.
2. **Un `std::vector` dans la branche « au-delà de 24 » de chaque
   évaluateur.** Le code d'allocation, en ligne dans la fonction chaude,
   l'empêche d'être mise en ligne chez ses appelants : 2,3 fois plus lent pour
   la valeur d'une courbe de degré 3.

Avec la *policy*, l'instanciation `stack_scratch` est l'ancien code (les mêmes
tableaux locaux, les mêmes boucles), et le compilateur voit les mêmes
fonctions ; les instanciations `heap_scratch` sont hors ligne et ne pèsent pas
sur les décisions de mise en ligne. Deux copies de l'algorithme (une pile, une
tas) auraient le même effet, au prix d'une double maintenance.

L'espace de noms est `basis_detail` et non `detail` : un code qui écrit
`using namespace gbs; using namespace gbs::brep;` puis `detail::…` (c'est le
cas de `tests/tests_brep_boolean.cpp`) devenait ambigu dès qu'un `gbs::detail`
était toujours déclaré. Même convention que `collocation_detail` et
`intersect_detail`.

## 4. Changements

| Fichier | Changement |
|---|---|
| `basisfunctions.ixx` | `basis_detail` : `stack_scratch` (avec `assert(n <= N)`), `heap_scratch`, A2.2, A2.3, `basis_ders` et le travail des évaluateurs sur l'intervalle, sur la *policy* ; entrées `*_heap` hors ligne, résultats par valeur. Les fonctions publiques (`basis_funcs`, `ders_basis_funs`, `basis_ders`, `eval_value_decasteljau`, `eval_ders_decasteljau`, `eval_rational_ders` de courbe et de surface) choisissent pile ou tas après leur prologue. La déclaration anticipée de `eval_rational_value_simple` (pour le repli par ordre) disparaît ; ses arguments par défaut passent sur sa définition |
| `modules/basisfunctions.cppm` | `<cassert>` dans le fragment global |
| `knotsfunctions.ixx` | `fill_basis_row`, `fill_basis_band` : un seul chemin, tampon sur la pile ou le tas selon p |
| `bscinterp.h` | `assemble_collocation_sparse` : idem ; `band_lu` inchangé |
| `bssinterp.h` | `build_poles` : idem |
| `bssapprox.h` | `approx` sous contraintes : idem, tampons choisis une fois avant la boucle |
| `bscshaping.h`, `bssshaping.h` | `move_to_constraints_delta` : la matrice des contraintes par `basis_ders` sur l'intervalle de nœuds (les seules colonnes non nulles), au lieu de `basis_function` récursif par colonne (§ 11, point 2) |
| `bscurve.h` | `BSCurveRational::value(u, d)` : pour d > 0, l'ordre d de `eval_rational_ders` (une passe) au lieu de `eval_rational_value_simple` (§ 11, point 3) |
| `gbsconstants.h` | commentaire de `bspline_stack_max_degree` : au-delà, le même algorithme sur le tas |
| `gbslib.h` | `GBS_NOINLINE` |

## 5. Résultats jusqu'au degré 24

Vérifié par un banc hors dépôt qui compile les en-têtes de `master` et ceux
de la branche dans deux exécutables, écrit tous les résultats en binaire et
les compare octet par octet (degrés 1 à 24, nœuds avec et sans nœuds doubles,
81 points par courbe, 121 par surface, toutes les dérivées) :

| Chemin | clang 21 arm64, -O2, -O3, `-ffp-contract=off`, `-ffp-contract=fast` |
|---|---|
| `BSCurve::value(u, d)`, `derivatives` ; `BSCurveRational::value(u)`, `derivatives` (jusqu'à l'ordre 24) | identique au bit |
| `BSSurface::value(u, v, du, dv)`, `derivatives` ; `BSSurfaceRational` | identique au bit |
| `eval_value_decasteljau` appelée par un pointeur de fonction (comme les liaisons Python) | identique au bit |
| `interpolate` et `approx` de courbes, `interpolate` de surface | pôles identiques au bit |
| `basis_funcs`, `ders_basis_funs` appelées directement | identique au bit |

Garanti : le même code et les mêmes opérations dans le même ordre. Pas
garanti : un appel direct d'un évaluateur que le compilateur met maintenant en
ligne chez un appelant où il ne l'était pas. Le compilateur peut y fusionner
autrement les multiplications-additions, comme il le faisait déjà dans `master`
d'un appelant à l'autre. Dans le banc, la valeur d'une courbe appelée depuis
`main` diffère ainsi au dernier bit pour 81 % des points, de 6,3e-15 en
relatif au plus. Non vérifié : GCC, MSVC, x86. Sur la CI, les compilateurs ne
forment pas de FMA : clang x86-64 sans `-march` (pas de FMA dans le jeu
d'instructions de base), MSVC en `/fp:precise`. L'arithmétique y est donc
celle du source dans tous les contextes.

Changent volontairement, à tout degré, depuis la revue (§ 11) :
`BSCurveRational::value(u, d)` pour d > 0 (l'ordre d de la passe A4.2 au lieu
de la récursion par ordre) et le façonnage (`moved_to_constraints` de courbe et
de surface, A2.3 au lieu de la récursion). Écart à `master` sur les pôles
façonnés et les dérivées rationnelles d'un banc (degrés 1 à 9) : 4,3e-14
relatif au plus.

## 6. Performances jusqu'au degré 24

Mesures finales, après la revue (§ 11), contre `master` a8eb9cb ; clang 21,
-O3, Apple M, machine au repos (charge < 5), moyenne de trois passages du
meilleur de 7.

À travers les classes (appels virtuels), millisecondes pour 10⁶ évaluations :

| Cas | master | branche | rapport |
|---|---|---|---|
| courbe, valeur, p = 3 | 19,8 | 18,2 | 0,92 |
| courbe, dérivée 1, p = 3 | 23,6 | 22,3 | 0,94 |
| courbe, `derivatives(u, 2)`, p = 3 | 39,1 | 37,8 | 0,97 |
| courbe, valeur, p = 9 | 59,5 | 56,7 | 0,95 |
| courbe, `derivatives(u, 2)`, p = 9 | 109,7 | 108,5 | 0,99 |
| courbe rationnelle, valeur, p = 3 | 21,8 | 19,8 | 0,91 |
| courbe rationnelle, dérivée 1, p = 5 | 137,5 | 49,3 | 0,36 (§ 11, point 3) |
| courbe rationnelle, `derivatives(u, 2)`, p = 9 | 116,0 | 112,9 | 0,97 |
| courbe rationnelle, valeur, p = 9 | 58,1 | 66,5 | **1,14** |
| surface, valeur, 3 × 3 | 28,5 | 24,6 | 0,87 |
| surface, `derivatives(u, v, 1, 1)`, 3 × 3 | 86,7 | 89,3 | 1,03 |
| surface rationnelle, `derivatives(u, v, 1, 1)` | 112,7 | 114,9 | 1,02 |
| interpolation de 2000 points, p = 3 (× 20) | 2,90 | 2,65 | 0,91 |

Appels directs, le type vu par le compilateur (bancs du relecteur, § 11,
point 1), nanosecondes par évaluation :

| Cas | master | branche |
|---|---|---|
| `BSCurve::value`, p = 1 / 3 / 5 / 9 | 3,2 / 8,0 / 18,4 / 48,6 | 3,6 / 8,5 / 17,9 / 50,0 |
| `BSCurve::derivatives(u, 1)`, p = 3 / 9 | 18,6 / 63,6 | 18,5 / 63,2 |
| `BSCurve::derivatives(u, 2)`, p = 3 / 9 | 25,5 / 91,4 | 25,2 / 91,1 |
| `eval_ders_decasteljau<double, 4>`, n = 2, p = 3 / 9 | 23,3 / 89,1 | 24,1 / 88,5 |
| `eval_rational_ders<double, 3>`, n = 2, p = 3 / 9 | 25,5 / 90,3 | 25,4 / 93,7 |
| surface, valeur, p = 3 / 9 | 22,9 / 186,3 | 21,0 / 168,4 |
| surface, dérivées (2, 2), p = 3 / 9 | 138,5 / 871,1 | 138,8 / 872,7 |

Un seul cas ralentit nettement, la valeur d'une courbe rationnelle de degré 9
(+14 %, stable sur trois mesures) ; les degrés 3 et 5 du même chemin sont plus
rapides. Son code est celui de la première version de la PR (+10 % alors) :
la structure des appels y est celle de `master`, c'est le code généré qui
diffère. La valeur d'une courbe par appel direct coûte de 0 à 11 % de plus
aux degrés 1 à 4 (moins d'une demi-nanoseconde), autant au-delà.

## 7. Précision au-delà de 24

| Comparaison | Écart relatif maximal |
|---|---|
| degrés 25 et 26 contre la récursion (`basis_function`, l'ancien chemin) : 21 points sur quatre intervalles (nœud double, bouts), ordres 0 à 4, valeurs, `derivatives` et chaque fonction de base | 2,5e-14 (ordre 4), 2e-15 jusqu'à l'ordre 2 |
| surface 26 × 3 contre la récursion, 15 points, ordres (0..2, 0..2) | 3,1e-14 |
| degrés 26, 30 et 40 contre l'algorithme de de Boor (combinaisons affines des pôles, dérivées par hodographes), 101 points, ordres 0 à 3 | 2,5e-15 |
| partition de l'unité au degré 40 | 2,7e-15 |
| dérivées d'une courbe rationnelle jusqu'à l'ordre 30 (table sur le tas) contre la récursion par ordre (ordres 8, 10, 12) | 7e-16 ; ordres 0 à 6 identiques au bit à ceux de la table sur la pile |

## 8. Limites, suite

- Le façonnage et `BSCurveRational::value(u, d)`, laissés ici à la première
  version, sont traités après la revue (§ 11, points 2 et 3).
- `eval_rational_value_simple` reste récursif sur l'ordre (2^d évaluations) :
  il n'a plus d'appelant dans `gbs/` ; la liaison Python
  `gbs.eval_rational_value_simple` l'expose.
- Code mort laissé tel quel : `eval_value_par` (récursif, sans appelant) et
  la surcharge `eval_value_decasteljau(…) -> T` (appelle `eval_index_span`,
  commenté).
- `basis_function` reste l'API récursive exposée (référence des tests,
  liaison Python `gbs.basis_function`).

## 9. Tests

| Test | Cas | Attendu |
|---|---|---|
| `high_degree_matches_recursive_basis` | degré 26, un intervalle de Bézier ; la récursion (0,4 s par évaluation) calculée une fois par ordre et partagée | fonctions de base (A2.2, A2.3), valeur et dérivées 1 et 2 de la courbe, de la surface 26 × 2 (ordres 0..2 × 0..2) : 1e-12 relatif |
| `high_degree_curves_against_de_boor` | degrés 30 et 40, quatre intervalles, 101 points (nœuds et bouts compris), par `BSCurve` | valeur et dérivées 1 à 3 par `value` et `derivatives` : 1e-12 de l'échelle de chaque ordre ; partition de l'unité à 1e-13 (quelques ms en tout ; des heures avant) |
| `high_degree_surfaces_against_de_boor` | degré 30 en u (3 en v) et 30 en v (2 en u), 81 points, `BSSurface` | ordres (0..2, 0..2) par `value` et `derivatives` : 1e-12 relatif |
| `high_degree_rational` | courbe rationnelle de degré 30 ; dérivées d'une cubique rationnelle jusqu'à l'ordre 30 | valeur et dérivée 1 contre A/w de de Boor : 1e-12 ; ordres 0 à 6 identiques au bit à `derivatives(u, 6)` ; ordres 8 et 10 contre la récursion par ordre `eval_rational_value_simple` : 1e-12 ; `value(u, d)` identique au bit à l'ordre d de `derivatives` (d = 1, 2, 10, 30) ; tout fini |
| `interp_approx_high_degree` (`tests_bsinterp.cpp`) | degré 26, données polynomiales : interpolation de courbe dense (45 points) et creuse (120), approximation de courbe (80 points, 35 pôles), surface par une grille, sous contraintes, approximation sous contraintes | passage aux points à 1e-12 ; reproduction du polynôme à 1e-7 (interpolation), 1e-5 (approximation de courbe : équations normales au degré 26, 3e-7 mesuré), 1e-8 et 1e-10 (surfaces) |
| `move_to_constraints_high_degree` (`tests_curves_shapping.cpp`, `tests_surfaces_shapping.cpp`) | façonnage d'une courbe de degré 30 (un point déplacé, une tangente gardée, un second point) et d'une surface de degré 30 × 2 (un point, ses deux dérivées premières) | contraintes tenues à 1e-8 relatif |
| `derivatives_match_per_order_value_rational` (`tests_curve_derivatives.cpp`, modifié) | la comparaison d'une passe contre la règle par ordre | contre `eval_rational_value_simple`, puisque `value(u, d)` est désormais la passe elle-même |
| Python `test_high_degree_evaluation` | `BSCurve3d` de degré 40, pôles aux abscisses de Greville | x(u) = u à 1e-12, x'(u) = 1 à 1e-10 ; `pytest-timeout` arrêterait un retour au coût exponentiel |

Les tests existants (ctest, pytest) passent.

## 10. Points à valider

1. La forme : un algorithme sur une *policy* de mémoire de travail
   (`stack_scratch`, `heap_scratch`), instanciations « tas » hors ligne,
   `GBS_NOINLINE` dans `gbslib.h`. L'alternative à un seul chemin (§ 3.1)
   coûtait jusqu'à 2,5 fois en vitesse et changeait des résultats au dernier
   bit.
2. La portée du « au bit près » (§ 5) : à travers les classes, les liaisons
   Python, les interpolations et approximations ; pas pour un appel direct
   que le compilateur met en ligne autrement (6e-15 relatif au plus).
3. `band_lu` (`bscinterp.h`) garde sa limite de 24 : au-delà, l'interpolation
   passe par le LU creux pivoté. Elle pourrait accepter les degrés élevés
   maintenant ; je la laisse, la matrice de collocation à ces degrés étant mal
   conditionnée et le chemin pivoté le plus sûr.
4. `brep_json_max_degree` (PR 252, `gbs-io/brep_json.h`) n'est plus imposé
   par le coût de l'évaluation : le relever est une décision à part.
5. Les dérivées d'une courbe rationnelle au-delà de l'ordre 24 changent de
   méthode (table sur le tas au lieu de la récursion par ordre). Leurs
   résultats changent donc ; ils étaient hors de portée au-delà de 25 environ.
6. Le façonnage et `BSCurveRational::value(u, d)` (d > 0) changent de
   méthode à tout degré (§ 5, § 11) : 4,3e-14 relatif au plus.

## 11. Revue de code

Revue de la PR 254 (25 novembre 2026). Un constat important, deux défauts
antérieurs, trois détails.

1. **Important, performances.** Quand le compilateur voit le type concret et
   met l'appel en ligne, les dérivées de courbe de tous ordres étaient de 15 à
   35 % plus lentes que sur `master` aux degrés 2 à 9 ; à travers un
   `Curve*` opaque (le cas du § 6), à ±4 %. Deux causes, trouvées par
   expériences successives :
   - le pointeur de sortie `CK` passé à l'entrée « tas » hors ligne
     s'échappe ; le compilateur ne peut plus garder les sommes de la branche
     « pile » dans des registres et les écrit en mémoire à chaque tour ;
   - le choix pile ou tas fait avant le prologue (sortie mise à zéro,
     `find_span`) change la façon dont le compilateur déroule et vectorise la
     boucle en ligne.

   Fait : les entrées `*_heap` rendent leurs points par valeur et l'appelant
   les recopie ; chaque évaluateur de dérivées garde le prologue de `master`
   et ne choisit qu'ensuite, le travail sur l'intervalle restant écrit une
   fois sur la *policy* (`ders_span`). Les évaluateurs de valeur gardent leur
   première forme (choix d'abord, corps entier sur la *policy*) : leur résultat
   revient par valeur, rien ne s'échappe, et leur donner le prologue les
   rendait deux fois plus lents en appel direct (degré 3 : 8,0 → 18,0 ns),
   mesuré avant de les remettre. Bancs du relecteur, avant et après
   (clang 21 -O3, arm64, ns) :

   | Cas | master | avant | après |
   |---|---|---|---|
   | `BSCurve::derivatives(u, 2)`, appel direct, p = 3 | 25,5 | 35,2 | 25,2 |
   | idem, p = 9 | 91,4 | 103,6 | 91,1 |
   | `BSCurve::derivatives(u, 1)`, p = 5 | 31,3 | 40,4 | 31,4 |
   | idem, toutes les sorties lues, p = 3, n = 2 | 30,1 | 40,1 | 29,7 |
   | `eval_ders_decasteljau<double, 4>`, n = 2, p = 3 | 23,3 | 35,2 | 24,1 |
   | `eval_rational_ders<double, 3>`, n = 2, p = 3 | 25,5 | 35,4 | 25,4 |
   | appel virtuel, n = 2, p = 3 | 36,9 | 35,8 | 36,3 |
   | surface, valeur, p = 3 | 22,9 | 21,4 | 21,0 |

   Le § 6 donne les mesures à travers les classes et selon le degré, appel
   direct compris.
2. **Antérieur, façonnage.** `bscshaping.h` et `bssshaping.h` remplissaient
   leur matrice de contraintes par `basis_function`, récursif, à tout degré.
   Fait : `basis_ders` (qui choisit pile ou tas) sur l'intervalle de nœuds de
   chaque contrainte, les seules colonnes non nulles, le reste à zéro.
   Test : `move_to_constraints_high_degree`, courbe de degré 30 et surface
   30 × 2. Les pôles façonnés changent d'au plus 4,3e-14 relatif à tout degré
   (§ 5).
3. **Antérieur, `BSCurveRational::value(u, d)`** passait par
   `eval_rational_value_simple`, 2^d évaluations. Fait : pour d > 0, l'ordre d
   d'une passe `eval_rational_ders` (sur la pile jusqu'à l'ordre 24) ; d = 0
   inchangé, identique au bit. Test : `value(u, d)` identique au bit à
   `derivatives` (d = 1, 2, 10, 30) ; `derivatives_match_per_order_value_rational`
   compare désormais à `eval_rational_value_simple`. Les dérivées changent
   d'au plus 4,3e-14 relatif (§ 5) ; d = 20 coûtait 2^20 évaluations, c'est
   une passe.
4. **Détail.** `stack_scratch` ignorait sa taille. Fait :
   `assert(n <= N)` dans son constructeur (`<cassert>` ajouté à l'en-tête et
   au fragment global du module).
5. **Détail.** `fill_basis_row`, `fill_basis_band`,
   `assemble_collocation_sparse` et `build_poles` avaient deux branches
   identiques. Fait : un seul chemin, un tampon choisi selon p (tableau sur
   la pile ou `std::vector`), la fonction publique qui choisit elle-même ;
   `bssapprox.h` et le façonnage suivent le même modèle, par direction.
6. **Détail.** `high_degree_curves_against_de_boor` mesurait le temps
   (`ASSERT_LT(ms, 2000.)`) et écrivait sur `std::cout`. Fait : les deux
   retirés.

Vérifié en local : ctest 956/956, pytest 192/192 ; résultats identiques au
bit à `master` jusqu'au degré 24 par les mêmes bancs qu'au § 5, hors les deux
changements voulus (points 2 et 3).
