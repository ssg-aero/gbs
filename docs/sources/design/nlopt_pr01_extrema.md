# Suppression de NLopt — PR 1 : extrema par Newton

| | |
|---|---|
| PR | branche `feat/nlopt-pr1-extrema-newton` |
| Document de référence | [nlopt_removal.md](nlopt_removal.md), § 3, 4 et 5 |
| Fichiers | `gbs/intersect/projection.h` (nouveau), `gbs/extrema.h`, `gbs/intersect/curves.h`, `gbs/intersect/curve_surface.h`, `python/gbsbindextrema.h`, `tests/tests_extrema.cpp`, `tests/tests_bscapprox.cpp`, `python/tests/test_BSCurve.py`, `docs/index.rst`, `docs/sources/design/nlopt_removal.md`, `News.md` |
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
  complète, puis JᵀJ si elle n'est pas définie positive, puis pas diagonal ;
  variables tenues sur une borne active ; pas divisé par deux tant que f
  remonte ; arrêt sur max(tol_x², 4 ε)(1 + |x|) ou 100 itérations ; reprise
  de l'autre côté d'une couture.
- Les dérivées secondes d'une géométrie qui ne les fournit pas
  (`OffsetSurface` de `gbs-brep/offset.h` lève une exception au-delà de
  l'ordre 1) sont abandonnées à la première exception : Gauss-Newton.
- Locales : `newton_point_curve`, `newton_point_surface`,
  `newton_curve_curve`, `newton_curve_surface`.
- Globales, depuis les deux meilleurs échantillons : `project_point_curve`,
  `project_point_surface`, `closest_curve_curve`, `closest_curve_surface`.
  Les courbes ne sont évaluées qu'une fois par échantillon (30 + 30 au lieu de
  2 × 900 évaluations en courbe-courbe, 30 + 900 au lieu de 2 × 27 000 en
  courbe-surface) ; les distances des grilles de surface sont calculées en
  parallèle au-delà de `parallel_min_size`.

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
`DeprecationWarning` est émis et il est ignoré. L'enum `nlopt_algorithm` reste
jusqu'à la PR 3.

### 2.5 Appelants

Aucun changement de code : avec la vraie distance, les contrôles de
`build_intersections`, `lattice_intersections` et
`detail::project_point_global` comparent enfin une distance à leur tolérance.

### 2.6 Écart au document de référence

La reprise de l'autre côté d'une couture (§ 3.2, point 5) a été ajoutée en
cours de route : sur un tore, les quatre coins du domaine sont le même point,
et les deux meilleurs échantillons tombaient du même côté de la couture en v ;
Newton restait alors contre la borne (distance 0,300002 au lieu de 0,3).

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

### 3.1 Test existant corrigé

`tests_bscapprox.approx_curve` vérifiait que les points d'une approximation
d'un cercle étaient à moins de 1e-4 du cercle. L'écart réel atteint 1,95e-4 :
seul son carré (3,8e-8) passait le seuil. Le test vérifie désormais que la
distance rendue est |‖p‖ − 1| à 1e-9 près, et que l'écart reste sous 3e-4.
Aucun autre test n'a dû être modifié : les seuils de `build_intersections`,
des treillis TFI et de `pcurve.h` passent avec la vraie distance sur les jeux
de données existants.

### 3.2 Résultats

Build complet (tests, rendu, bindings Python ; clang conda-forge, Ninja,
Release, macOS arm64) : les 480 tests ctest passent (après fusion de master) (dont les suites
`intersect`, `brep`, `step`, `mesh`, `bssbuild`), ainsi que les 53 tests
Python. Le binaire OCCT (configuration séparée) n'a pas été relancé ; les
intersections qu'il exerce passent par `project` et `project_on_surface`,
couverts par les suites `intersect` et `brep`.
