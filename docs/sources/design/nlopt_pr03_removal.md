# Suppression de NLopt — PR 3 : retrait de la dépendance

| | |
|---|---|
| PR | branche `feat/nlopt-pr3-removal` |
| Document de référence | [nlopt_removal.md](nlopt_removal.md), § 4.2 et 5 |
| Fichiers | `gbs/solvers.h` (supprimé), `gbs/extrema.h`, `gbs/bscanalysis.h`, `python/gbsbindextrema.h`, `python/gbsbind.cpp`, `CMakeLists.txt`, `cmake/GbsConfig.cmake.in`, `python/CMakeLists.txt`, `tests/CMakeLists.txt`, `gbs-occt/tests/CMakeLists.txt`, `gbs-render/examples/CMakeLists.txt`, `tests_3rdlibs/CMakeLists.txt`, `bench/CMakeLists.txt`, `bench/bench_extrema.cpp`, `.github/workflows/ci.yml`, `ReadMe.md`, `tests/tests_solvers.cpp` et `tests_3rdlibs/tests_nlopt_example.cpp` (supprimés), `tests/tests_extrema.cpp`, `python/tests/test_BSCurve.py`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Après les PR 1 et 2, NLopt ne calcule plus rien dans gbs : il ne restait que
des déclarations (le type du paramètre `solver`, ignoré), les solveurs
génériques sans appelant, leurs tests, et le build. Cette PR les retire : gbs
ne dépend plus de NLopt, ni à la compilation ni à l'exécution.

## 2. Changements

### 2.1 API C++ (rupture, annoncée dans `News.md`)

- `gbs/solvers.h` est supprimé (`solve_N_nlop`, `solve_D_nlop`,
  `nlopt_optimize`, `nlopt_bounds`) ; `default_nlopt_algo` aussi.
- Paramètre `solver` retiré de toutes les fonctions `extrema_*` et de
  `max_curvature_pos`.
- Paramètres `n_bracket` retirés des extrema de courbes
  (`extrema_curve_point`, `extrema_curve_curve`, `extrema_curve_curve_lst`),
  sans effet depuis la subdivision de la PR 1. Ceux des extrema de surfaces
  restent : ils règlent la grille des surfaces hors 3D.
- Sans le paramètre `solver` intercalé, `extrema_surf_pnt(srf, p, u0, v0, tol)`
  et la version à grille avaient le même nombre d'arguments : des flottants
  pour la grille (`extrema_surf_pnt(srf, p, 1e-6, 50., 50.)`) partaient
  silencieusement de (1e-6, 50) avec tol = 50. La grille des surfaces est
  donc un tableau : `extrema_surf_pnt(srf, p, tol, std::array<std::size_t, 2>
  n_bracket = {30, 30})`, `extrema_surf_curve(srf, crv, tol,
  std::array<std::size_t, 3> n_bracket = {30, 30, 30})` (u, v de la surface,
  paramètre de la courbe) ; aucun argument flottant ne peut y arriver.

### 2.2 Python

L'énumération `nlopt_algorithm` et l'argument `solver` disparaissent, ainsi
que `n_bracket` de `extrema_curve_point` et `n_bracket_u` / `n_bracket_v` de
`extrema_curve_curve`. Pour que les anciennes formes positionnelles soient
refusées plutôt que lues autrement :

- la grille des surfaces (`n_bracket_u`, `n_bracket_v`, `n_bracket_w`) est
  en arguments nommés seulement (`py::kw_only`) ;
- les anciennes formes positionnelles avec des entiers
  (`extrema_curve_point(crv, p, 1e-6, 30)`, `extrema_curve_curve(c1, c2, tol,
  30, 30)`, `extrema_surf_pnt(srf, p, tol, 30, 30)`, `extrema_surf_curve(…,
  30, 30, 30)`) ont l'arité des versions à départ et partiraient de u0 = tol
  avec une tolérance de 30. `noconvert` sur la tolérance ne suffit pas :
  pybind11 3 accepte un entier pour un flottant même sans conversion (tour
  numérique de la PEP 484). Des surcharges « garde », enregistrées en premier,
  prennent ces formes avec des entiers stricts et lèvent un `TypeError` qui dit
  quoi écrire à la place.

Cinq flottants positionnels pour `extrema_surf_pnt` restent un départ et une
tolérance : c'est la signature de la version à départ, Python ne peut pas
les distinguer d'une grille ; la grille ne s'écrit plus que nommée.

Les fonctions mortes de `python/gbsbind.cpp` (`extrema_curve_point(py::args)`,
`extrema_curve_point_`, `extrema_curve_point_3d`, `conv_arr` et les blocs
commentés `extrema_PC_`), jamais enregistrées, sont supprimées. Le test de dépréciation devient
`test_extrema_without_nlopt` : `solver=`, `n_bracket=` et les anciennes
formes positionnelles lèvent `TypeError`, la grille nommée des surfaces
fonctionne. Une copie du bloc
`-W error` de la PR 1, tombée par erreur dans `test_extrema_true_distance`, est
retirée.

### 2.3 Build et CI

- `find_package(NLopt)` retiré de `CMakeLists.txt` et `find_dependency(NLopt)`
  de `GbsConfig.cmake.in` ; `NLopt::nlopt` retiré des cibles `gbs`,
  `gbs_pch_base`, Python, tests, tests OCCT, exemples de rendu,
  `tests_3rdlibs` et du bench.
- `tests/tests_solvers.cpp` et `tests_3rdlibs/tests_nlopt_example.cpp`
  supprimés.
- CI : `nlopt` retiré des trois blocs d'environnement conda (Windows, Unix
— qui sert aussi au job macOS Intel — et OCCT) et des commentaires.
- `ReadMe.md` : `nlopt` retiré des paquets nécessaires.
- `bench/bench_extrema.cpp` ne passe plus que des arguments communs aux deux
  versions (il compile encore contre une version d'avant le retrait, liée à
  `-lnlopt`) : les 20 requêtes B-spline varient par une translation de la
  seconde courbe au lieu de `n_bracket`.

### 2.4 Hors de ce dépôt

La recette conda-forge (`gbs-feedstock`) déclare probablement `nlopt` dans ses
dépendances : à retirer lors de la prochaine version publiée.

### 2.5 Revue de code

| Commentaire | Suite donnée |
|---|---|
| même arité pour les versions à départ et à grille de `extrema_surf_pnt` / `extrema_surf_curve` : des flottants pour la grille partent d'un mauvais départ | grille en `std::array<std::size_t, N>` en C++ ; en Python, grille nommée seulement, et surcharges « garde » qui refusent les anciennes formes avec des entiers (`noconvert` ne suffit pas avec pybind11 3). Contrairement à ce que suppose la remarque, cinq flottants positionnels restent valides en Python : c'est la version à départ |
| fonctions mortes de `gbsbind.cpp` modifiées au lieu d'être supprimées | supprimées (environ 120 lignes) |
| le test Python ne vérifie que `solver=None` | anciennes formes positionnelles et nommées des courbes et des surfaces, grille nommée |
| commentaire de `ci.yml` qui justifie l'absence de Qt par NLopt | la phrase sur NLopt est retirée |
| deux inexactitudes de la note (résultats du bench, blocs conda) | corrigées (§ 2.3 et § 3) |

## 3. Tests

- Build complet dans un répertoire neuf (tests, rendu, bindings Python) :
  CMake ne cherche plus NLopt, aucune source n'inclut `nlopt.hpp`, et
  `otool -L` ne montre `libnlopt` ni dans les tests ni dans le module Python.
- Les 489 tests ctest passent (494 avant : les 4 de `tests_solvers` et
  `solver_ignored` sont retirés), ainsi que les 54 tests Python, dont
  `test_extrema_without_nlopt` (macOS arm64, clang conda-forge, Release).
- `bench/bench_extrema.cpp` compile et s'exécute sans `-lnlopt`. Ses
  résultats sont ceux de la note de la PR 1, sauf la ligne « B-splines,
  global » : ses 20 requêtes ont changé (seconde courbe translatée de 0,05 i
  au lieu d'un `n_bracket` ignoré), ses temps et distances ne sont donc pas
  comparables à ceux de cette note.
- Le binaire OCCT (configuration séparée) n'a pas été relancé : seule sa
  ligne de lien change.
