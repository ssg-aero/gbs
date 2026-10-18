# Suppression de NLopt — PR 1 : extrema par Newton

| | |
|---|---|
| PR | branche `feat/nlopt-pr1-extrema-newton` |
| Document de référence | [nlopt_removal.md](nlopt_removal.md), § 3, 4 et 5 |
| Fichiers | `gbs/intersect/projection.h` (nouveau), `gbs/extrema.h`, `gbs/intersect/curves.h`, `gbs/intersect/curve_surface.h`, `gbs/execution.h`, `gbs/bssbuild.h`, `gbs-mesh/tfi.h`, `gbs/bscurve.h`, `gbs/bssurf.h`, `gbs/curveonsurface.h`, `gbs-brep/offset.h`, `python/gbsbindextrema.h`, `tests/tests_extrema.cpp`, `tests/tests_bscapprox.cpp`, `tests/tests_mesh.cpp`, `tests/tests_brep_geom.h`, `python/tests/test_BSCurve.py`, `bench/bench_extrema.cpp`, `bench/CMakeLists.txt`, `bench/python/compare_extrema.py`, `docs/index.rst`, `docs/sources/design/nlopt_removal.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Les extrema de `gbs/extrema.h` ne passent plus par NLopt : ils sont calculés
par la méthode de Newton de `gbs/intersect/projection.h`, partagée avec les
projections des intersections, et rendent la vraie distance au lieu de son
carré. NLopt reste lié : `gbs/extrema.h` inclut encore `gbs/solvers.h` pour le
type du paramètre `solver`, ignoré, et `bscanalysis.h` l'utilise toujours
(PR 2).

## 2. Changements

### 2.1 `gbs/intersect/projection.h`

- `projection_detail::minimize<T, n>` : Newton amorti sur f = ½‖r‖², n = 1 à 3
  (§ 3.2 du document de référence). Pas résolu par Cholesky ; hessienne
  complète, puis JᵀJ si elle n'est pas définie positive, puis
  Levenberg-Marquardt, puis pas diagonal ; variables tenues sur une borne
  quand le gradient ou le pas sort de la boîte ; pas divisé par deux tant que
  f ne décroît pas ; descente de gradient projetée si x ne peut pas bouger ;
  arrêt sur max(tol_x², 4 ε)(1 + |x|) ou 100 itérations ; reprise de l'autre
  côté d'une couture.
- `Curve::has_second_derivatives()` et `Surface::has_second_derivatives()`
  (vrai par défaut) : `OffsetSurface` (`gbs-brep/offset.h`, dérivées premières
  seulement) rend faux, `CurveOnSurface` rend la réponse de sa surface. Newton
  utilise alors Gauss-Newton.
- Locales : `newton_point_curve`, `newton_point_surface`,
  `newton_curve_curve`, `newton_curve_surface`.
- Globales, depuis les deux meilleurs échantillons : `project_point_curve`,
  `project_point_surface`, `closest_curve_curve`, `closest_curve_surface`.
  Les courbes ne sont évaluées qu'une fois par échantillon (30 + 30 au lieu de
  2 × 900 évaluations en courbe-courbe, 30 + 900 au lieu de 2 × 27 000 en
  courbe-surface) ; les distances des grilles sont calculées en parallèle
  au-delà de `parallel_min_size` (par ligne en courbe-courbe et
  courbe-surface, puis fusion déterministe des deux meilleurs de chaque
  ligne).

### 2.2 `gbs/extrema.h`

Chaque fonction appelle la fonction locale ou globale correspondante ; types
de retour inchangés, `d` est la distance. `solver` est marqué
`[[maybe_unused]]` et documenté comme ignoré. Les surcharges sans départ sur
une surface vérifient que ses bornes sont finies. `extrema_curve_curve_lst`
passe ses `n_bracket_u` / `n_bracket_v`.

### 2.3 Intersections

`intersect_detail::project` (17 échantillons) et
`intersect_detail::project_on_surface` (grille 5 × 5) appellent
`project_point_curve` / `project_point_surface`. Le second passe ainsi de
Gauss-Newton à Newton. `refine` et `refine_cs` ne changent pas.

### 2.4 Python

L'argument `solver` vaut `None` par défaut ; s'il est passé, un
`DeprecationWarning` est émis et il est ignoré. Si les avertissements sont
changés en erreurs (`-W error`), l'exception est propagée. L'enum `nlopt_algorithm` reste
jusqu'à la PR 3.

### 2.5 Appelants

Avec la vraie distance, les contrôles de `build_intersections` et de
`lattice_intersections` comparent enfin une distance à leur tolérance, comme
le faisait déjà la grille de `gordon`. **Rupture assumée** : à `tol` = 1e-6,
des courbes jusqu'à 1e-3 l'une de l'autre passaient, seulement jusqu'à 1e-6
désormais (`News.md`) ; le seuil le plus strict a été retenu plutôt qu'un
paramètre de contact séparé.

Une exception qui sort d'un algorithme parallèle appelle `std::terminate` :
`gbs/execution.h` gagne `transform_rethrow`, qui capture les exceptions de
chaque élément et relance la première après la région parallèle (comme
`build_batch`). `build_intersections` calcule les extrema en parallèle puis
contrôle les distances en séquentiel ; `get_tg_from_spine`,
`extrema_curve_curve_lst` et `project_points_on_surface` passent par
`transform_rethrow`. `closest_curve` / `closest_p_curve` (`bscanalysis.h`)
ont le même défaut, et ne compilent pas une fois instanciés (`midPoint()`) :
ils sont repris en PR 2.

### 2.6 Écarts au document de référence

Deux ajouts en cours de route, tous deux pour des échantillons confondus :

- **reprise de l'autre côté d'une couture** (§ 3.2, point 5) : sur un tore,
  les quatre coins du domaine sont le même point, et les deux meilleurs
  échantillons tombaient du même côté de la couture en v ; Newton restait
  contre la borne (distance 0,300002 au lieu de 0,3) ;
- **second départ distinct du premier** : les 30 échantillons du bord v = π/2
  d'une sphère sont tous le pôle. Pour un point à 0,039 rad du pôle, les deux
  départs étaient deux copies du pôle, où la dérivée en u est nulle et v contre
  sa borne : Newton n'en sortait pas (0,3562 au lieu de 0,3534 ; NLopt faisait
  la même erreur). Un échantillon à la même distance que le meilleur n'est plus
  pris comme second départ (les copies d'une couture sont laissées à la
  reprise).

## 3. Tests

`tests/tests_extrema.cpp`, en `double` et en `float` :

| Test | Cas | Attendu |
|---|---|---|
| `true_distance_PC` | cercle rationnel, points à 0,25 dehors, 0,4 dedans, 0,3 au-dessus, aux angles 0, ±1e-3 (couture), 1, π ; départ local décalé ; arc d'un quart de tour et point avant son début | la distance (pas son carré), le pied attendu ; le début de l'arc |
| `true_distance_PS` | cylindre (point à 0,5 ; point au-delà du bord v = 2), sphère (point à 1), tore (point à 0,3 près des deux coutures) | distance ; sur le bord : v = 2 et pied exact en u (borne active) |
| `true_distance_CC` | segments gauches à 0,3 ; cercle et segment au-dessus de sa couture, global et local ; `extrema_curve_curve_lst` avec `n_bracket` | distances et pieds |
| `true_distance_CS` | sphère et segment à 0,5, autour de la couture, global et local | distance et pieds |
| `solver_ignored` | `LD_MMA` contre défaut | même résultat |
| `build_intersections_threshold` | courbes à 1e-4 avec `tol` = 1e-6 | `std::invalid_argument` (l'ancien code les acceptait) |

`python/tests/test_BSCurve.py` : vraie distance point-courbe et courbe-courbe,
`DeprecationWarning` quand `solver` est passé.

### 2.7 Revue de code

| Commentaire | Suite donnée |
|---|---|
| CI rouge (Ubuntu, Windows) : `build_intersections` lève dans un `std::transform(par)`, donc `std::terminate` ; macOS passe car `GBS_PAR_EXEC` y est vide | extrema en parallèle par `transform_rethrow`, contrôle hors de la boucle ; même traitement des autres boucles parallèles qui appellent les extrema (§ 2.5) |
| retour de `PyErr_WarnEx` ignoré : `SystemError` sous `-W error` | `throw py::error_already_set()` ; test Python sous `simplefilter("error")` |
| pas de Newton rogné sans mouvement : arrêt sur une borne avec un gradient projeté non nul | variable tenue quand son pas sort de la boîte, pas recalculé ; descente de gradient projetée si x ne bouge pas ; test `newton_step_out_of_the_box` (pavé cisaillé, minimum sur le bord u = 0, départs global et locaux) |
| les extrema globaux lèvent sur une surface non bornée, `pcurve.h` attend `OffSurface` | le cas ne se présente pas : une primitive refuse une boîte infinie à la construction (`primitives.h`), et aucune surface de gbs n'est non bornée ; le contrôle reste pour les surfaces d'utilisateur, qui donnaient des NaN |
| `second_order` avale toute `std::exception` ; une exception par appel sur une surface décalée | `has_second_derivatives()` sur `Curve` et `Surface`, plus de `try` ; test `no_second_derivatives` (surface décalée d'un cylindre). Les géométries de test `Analytic` / `AnalyticCurve` (`tests_brep_geom.h`, dérivées premières seulement) faisaient échouer deux tests BREP une fois le `try` retiré : elles redéfinissent la méthode. **Rupture** pour une classe d'utilisateur dans ce cas, qui marchait avec l'ancien Gauss-Newton : annoncée dans `News.md` |
| seuils de contact cent fois plus stricts sans changement des défauts | décision : le seuil le plus strict, `tol` est une distance partout ; rupture dans `News.md` ; test `lattice_contact_threshold` (treillis à 1e-4 refusé à 1e-6, accepté à 1e-3) |
| balayage séquentiel des 27 000 triplets en courbe-surface | deux meilleurs par ligne en parallèle, puis fusion |
| tests manquants : départ au pôle, surface sans dérivées secondes, repli Gauss-Newton, pas diagonal, point sur l'axe d'un cylindre | `pole_and_axis_starts`, `no_second_derivatives`, `gauss_newton_and_diagonal_fallbacks` (quart de cercle de poids 1 / 7,07 / 100 vu de loin, segments parallèles). Les segments parallèles en départ local **échouaient** : le pas diagonal permutait les deux paramètres à distance égale, accepté sans fin ; d'où Levenberg-Marquardt et le refus d'un grand pas à f égale |

### 3.1 Test existant corrigé

`tests_bscapprox.approx_curve` vérifiait que les points d'une approximation
d'un cercle étaient à moins de 1e-4 du cercle. L'écart réel atteint 1,95e-4 :
seul son carré (3,8e-8) passait le seuil. Le test vérifie désormais que la
distance rendue est |‖p‖ − 1| à 1e-9 près, et que l'écart reste sous 3e-4.
Aucun autre test n'a dû être modifié : les seuils de `build_intersections`,
des treillis TFI et de `pcurve.h` passent avec la vraie distance sur les jeux
de données existants.

### 3.2 Performances

`bench/bench_extrema.cpp`, compilé avec les mêmes options (`-O3`, clang
conda-forge, macOS arm64) contre les en-têtes de master (NLopt, COBYLA) et de
cette branche. Le temps est le meilleur de trois passes ; les évaluations
comptent les appels à `value()` de la géométrie (indépendant de la machine) ;
l'erreur est l'écart entre la vraie distance, recalculée depuis les paramètres
rendus (master rend d²), et la distance exacte.

| Cas | NLopt µs | Newton µs | NLopt éval. | Newton éval. | NLopt > 1e-6 | Newton > 1e-6 |
|---|---:|---:|---:|---:|---:|---:|
| point-cercle, global | 8,9 | 3,8 | 112 | 72 | 2/200 (8,7e-3) | 0/200 |
| point-cercle, local | 7,0 | 1,7 | 55 | 21 | 1/200 | 0/200 |
| point-B-spline, global | 8,5 | 2,4 | 113 | 77 | — | — |
| point-sphère, global | 51,9 | 30,6 | 1 910 | 978 | 3/200 | 0/200 |
| point-sphère, local | 73,9 | 2,8 | 547 | 35 | 1/200 | 0/200 |
| point-tore, global | 402,9 | 29,1 | 4 490 | 979 | 3/200 | 0/200 |
| cercle-segment, global | 680,3 | 4,3 | 11 120 | 151 | — | — |
| B-splines, global | 45,8 | 4,2 | 2 568 | 166 | — | — |
| sphère-segment, global | 1 172 | 52,4 | 108 171 | 1 051 | **17/20 exceptions** | 0/20 |
| sphère-segment, local | 659,3 | 3,1 | 6 101 | 66 | 5/200 (3,9e-2) | 0/200 |

- Newton atteint la distance exacte à 1e-15 près partout ; requête par
  requête (`bench/python/compare_extrema.py`), il n'est jamais moins bon que
  NLopt, et meilleur sur 17 des 1 640 requêtes (jusqu'à 0,21 sur la B-spline).
- Gain de temps de 1,7 (point-sphère global, dominé par les 900 échantillons)
  à 210 (sphère-segment local). Les surcharges globales coûtent surtout
  l'échantillonnage, qui n'a pas changé de densité ; en courbe-courbe et
  courbe-surface, les courbes ne sont plus évaluées qu'une fois par
  échantillon.
- Sur master, 17 des 20 appels globaux sphère-segment **lèvent une
  exception** : COBYLA évalue le segment hors de ses bornes
  (`std::invalid_argument` de `BSCurve`).

### 3.3 Résultats des tests

Build complet (tests, rendu, bindings Python ; clang conda-forge, Ninja,
Release, macOS arm64) : les 485 tests ctest passent, après fusion de master et la revue (dont les suites
`intersect`, `brep`, `step`, `mesh`, `bssbuild`), ainsi que les 53 tests
Python. Le binaire OCCT (configuration séparée) n'a pas été relancé ; les
intersections qu'il exerce passent par `project` et `project_on_surface`,
couverts par les suites `intersect` et `brep`.
