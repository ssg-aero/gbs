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
  et `extrema_surf_pnt(srf, p, tol, nu, nv)` ont le même nombre d'arguments :
  `T` est déduit de la surface, donc des entiers ne peuvent aller qu'à la
  seconde, et des flottants préfèrent la première (correspondance exacte). Même
  chose pour `extrema_surf_curve`.

### 2.2 Python

L'énumération `nlopt_algorithm` et l'argument `solver` disparaissent, ainsi
que `n_bracket` de `extrema_curve_point` et `n_bracket_u` / `n_bracket_v` de
`extrema_curve_curve`. Le test de dépréciation devient
`test_extrema_without_nlopt` (`solver=` lève `TypeError`). Une copie du bloc
`-W error` de la PR 1, tombée par erreur dans `test_extrema_true_distance`, est
retirée.

### 2.3 Build et CI

- `find_package(NLopt)` retiré de `CMakeLists.txt` et `find_dependency(NLopt)`
  de `GbsConfig.cmake.in` ; `NLopt::nlopt` retiré des cibles `gbs`,
  `gbs_pch_base`, Python, tests, tests OCCT, exemples de rendu,
  `tests_3rdlibs` et du bench.
- `tests/tests_solvers.cpp` et `tests_3rdlibs/tests_nlopt_example.cpp`
  supprimés.
- CI : `nlopt` retiré des trois environnements conda (matrice principale,
  macOS Intel, OCCT) et des commentaires.
- `ReadMe.md` : `nlopt` retiré des paquets nécessaires.
- `bench/bench_extrema.cpp` ne passe plus que des arguments communs aux deux
  versions (il compile encore contre une version d'avant le retrait, liée à
  `-lnlopt`) : les 20 requêtes B-spline varient par une translation de la
  seconde courbe au lieu de `n_bracket`.

### 2.4 Hors de ce dépôt

La recette conda-forge (`gbs-feedstock`) déclare probablement `nlopt` dans ses
dépendances : à retirer lors de la prochaine version publiée.

## 3. Tests

- Build complet dans un répertoire neuf (tests, rendu, bindings Python) :
  CMake ne cherche plus NLopt, aucune source n'inclut `nlopt.hpp`, et
  `otool -L` ne montre `libnlopt` ni dans les tests ni dans le module Python.
- Les 489 tests ctest passent (494 avant : les 4 de `tests_solvers` et
  `solver_ignored` sont retirés), ainsi que les 54 tests Python, dont
  `test_extrema_without_nlopt` (macOS arm64, clang conda-forge, Release).
- `bench/bench_extrema.cpp` compile et s'exécute sans `-lnlopt`, avec les
  mêmes résultats que dans la note de la PR 1.
- Le binaire OCCT (configuration séparée) n'a pas été relancé : seule sa
  ligne de lien change.
