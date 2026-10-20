# Minimisation sans dérivées dans une boîte : Nelder-Mead

| | |
|---|---|
| PR | branche `feat/nelder-mead` |
| Document de référence | [nlopt_removal.md](nlopt_removal.md), § 1.2 et 8 |
| Fichiers | `gbs/minimize.h` (nouveau), `tests/tests_minimize.cpp` (nouveau), `docs/sources/design/nlopt_removal.md`, `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif

La suppression de NLopt a retiré `solve_N_nlop`, l'optimiseur générique sans
dérivées : il n'avait plus d'appelant dans gbs (`nlopt_removal.md`, § 1.2). Il
en avait ailleurs : tbs l'utilise pour ajuster les lignes de cambrure à leurs
angles métal (2 à 5 variables, somme de valeurs absolues), pour le col d'un
canal (distance de deux côtés dans le plan m-rθ, avec le pas) et dans ses
tests. Ces objectifs n'ont pas de dérivées commodes et ne sont pas lisses ;
COBYLA, sans dérivées, leur convenait.

Cette PR rend à gbs un minimiseur sans dérivées sous bornes, sans dépendance,
pour les bibliothèques qui s'appuient sur gbs (tbs, tbf). Les projections
restent le rôle de `gbs/extrema.h` (Newton, plus rapide et exact), les
fonctions d'une variable celui de Brent (`boost::math::tools::brent_find_minima`).

## 2. API

```cpp
template <std::floating_point T, std::size_t n, typename F>
auto nelder_mead(const F &f, std::array<T, n> x0, const std::array<T, n> &lo, const std::array<T, n> &hi, T tol_x,
                 std::size_t max_evaluations = 1000 * n, T initial_step = T(0.1)) -> MinimizeResult<T, n>; // {x, f, evaluations, converged}
```

- `std::array` de taille fixe, comme `gbs/intersect/projection.h` : les
  appelants connaissent n à la compilation, le simplexe tient sur la pile.
  `T` est déduit des tableaux seuls (`std::type_identity_t<T>` pour `tol_x` et
  `initial_step`) : `nelder_mead(f, x0, lo, hi, 1e-6)` compile avec des
  tableaux `float`.
- `f` est pris par référence universelle et appelé non const : une lambda
  `mutable` ou un foncteur qui compte ses appels sont acceptés. Sa valeur
  doit être convertible en T, ce que la clause `requires` vérifie à l'appel.
- Bornes infinies : variable libre ; `lo[i] == hi[i]` : variable tenue.
- `f` rend NaN là où il n'est pas défini : le point est pris comme le pire
  (+∞). `converged` est faux quand le budget est épuisé ou que f n'a été
  défini nulle part où le simplexe est allé.
- `initial_step` : les pas du simplexe initial (et des redémarrages), relatifs
  à la boîte sur chaque variable (à max(1, |x0ᵢ|) sans borne), un dixième par
  défaut. Il choisit parmi les bassins autour de x0 (§ 3.4).
- Arguments invalides (`tol_x <= 0`, `lo > hi`, `x0` non fini, budget nul,
  `initial_step` non positif ou non fini) : `std::invalid_argument`.

## 3. Algorithme

Nelder-Mead avec les coefficients d'origine (réflexion 1, expansion 2,
contraction ½, réduction ½), contraction extérieure ou intérieure selon que
le point réfléchi améliore ou non le pire sommet. Simplexe initial : un pas
d'un dixième de la boîte sur chaque variable, vers l'intérieur (un dixième de
max(1, |x0ᵢ|) sans borne). Arrêt quand chaque sommet est à moins de
tol_x (1 + |xᵢ|) du meilleur sur chaque coordonnée, en x.

### 3.1 Bornes par changement de variables

Le simplexe évolue en z, sans bornes : x = lo + (hi − lo)(1 + sin z)/2 entre
deux bornes, x = lo + z² ou hi − z² avec une seule, x = z sans. Ramener les
points d'essai sur les bornes (`clamp`) a été essayé d'abord : une réflexion
rognée pose des sommets sur une face de la boîte, le centroïde y passe, et le
simplexe, aplati sur la face, n'en sort plus. Exemple (test
`corner_with_kink`) : minimum en (0,02 ; 1), départ (0,3 ; 0,5) ; avec
`clamp` la méthode s'arrête en (0 ; 1), f = 0,031 au lieu de 0,010, même avec
les redémarrages.

### 3.1 bis Branches du changement de variables

Le sinus est périodique et symétrique : z, π − z et z + 2π donnent le même x
(et ±z avec une seule borne). Des sommets écartés d'une période (z = 0,1 et
0,1 + 2π) ont leur centroïde en 0,1 + π, point miroir de la boîte (revue de
la PR). Chaque point d'essai est donc décalé d'un nombre entier de périodes
vers la période du meilleur sommet, x inchangé.

Les symétriques (π − z, −z) sont laissés tels quels. Les replier sur la
branche principale (`asin(sin z)`, `|z|`) a été essayé : un minimum sur une
borne est un minimum lisse en z (sin z stationnaire en π/2), que le simplexe
enjambe ; replié, il s'écrase contre la borne. Sur les campagnes du § 3.2, le
repli coûte 5 à 20 % d'évaluations de plus et laisse 5 problèmes non
convergés en dimension 5 au lieu de 2 (écart maximal 1,6e-5 au lieu de
4e-6) ; choisir parmi z et π − z le plus proche du meilleur sommet revient au
même. Le seul décalage de période ne change rien aux campagnes (une exécution
sur 2 400 l'utilise).

### 3.2 Redémarrages

Nelder-Mead peut converger vers un point qui n'est pas un minimum quand le
simplexe dégénère (aplati sur une droite, contre un pli de f). Une fois
convergé, il repart du meilleur point avec un simplexe neuf tant que f baisse
(au plus 5 fois). Sur des sommes séparables w (x − xs)² + |x − xs| dans
[0, 1]ⁿ, minimum connu (xs ramené dans la boîte), 200 problèmes par
dimension :

| n | échecs sans redémarrage | échecs avec | évaluations sans / avec |
|---:|---:|---:|---:|
| 2 | 0 | 0 | 113 / 226 |
| 3 | 15 (jusqu'à 4,4e-4) | 0 | 251 / 544 |
| 4 | 33 (jusqu'à 2,7e-2) | 0 | 522 / 1 110 |
| 5 | 97 (jusqu'à 8,6e-2) | 2 (budget épuisé, 4e-6) | 1 032 / 2 288 |

Le redémarrage double à peu près le coût ; c'est le prix d'un résultat juste.

### 3.3 Coefficients adaptatifs écartés

Les coefficients de Gao et Han (2012), fonctions de n, ont été essayés : avec
les redémarrages ils n'améliorent ni l'exactitude ni le coût jusqu'à n = 12
(quadratique couplée : 11 170 évaluations contre 6 550 avec les coefficients
d'origine), et pour n ≤ 2 ils leur sont égaux. Gardés simples.

### 3.4 Pas initial

Sur un objectif à plusieurs minima, le pas initial choisit le bassin : un grand
pas regarde plus loin de x0, un petit reste près de lui. Aucune valeur n'est
meilleure en général : sur les campagnes du § 3.2, un cinquième de la boîte
au lieu d'un dixième fait autant d'échecs à quelques problèmes près. Sur les
lignes de cambrure de tbs (ajustement de 2 à 5 paramètres aux angles métal,
somme de valeurs absolues, plusieurs minima), un dixième laisse une ligne sur
quatorze à 17° de l'angle visé, un cinquième les amène toutes à moins de 0,1° ;
COBYLA (NLopt) en manquait six sur quinze de plus de 1°, jusqu'à 17°. Le pas est donc un paramètre, un
dixième par défaut.

## 4. Tests

`tests/tests_minimize.cpp` : Rosenbrock, quadratique couplée en 5D, bornes
actives (coin), objectif non lisse à une variable, départ hors de la boîte,
variable libre et variable tenue, pas initial (double puits incliné), décalage de période (cas de la revue), lambda `mutable`, foncteur à état, tolérance `double` avec des tableaux `float`, refus à l'appel d'une valeur non convertible, NaN sur une partie du domaine et au départ,
f indéfini partout, budget, `float`, arguments invalides, accord avec
`extrema_curve_curve` sur la distance de deux courbes, redémarrages (40
problèmes séparables tirés par `std::mt19937`, sans la distribution de la
bibliothèque standard pour rester identiques partout). Mutations vérifiées :
sans redémarrage, sans le traitement de NaN, sans expansion, chaque fois un
test échoue.

## 5. Hors périmètre

- Liaison Python : SciPy (`scipy.optimize.minimize(method="Nelder-Mead",
  bounds=...)`) couvre le besoin côté Python.
- Contraintes autres que des bornes.
