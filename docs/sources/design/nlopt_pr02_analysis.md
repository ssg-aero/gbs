# Suppression de NLopt — PR 2 : analyse des courbes

| | |
|---|---|
| PR | branche `feat/nlopt-pr2-analysis` |
| Document de référence | [nlopt_removal.md](nlopt_removal.md), § 3.4 et 5 |
| Fichiers | `gbs/bscanalysis.h`, `tests/tests_bscanalysis.cpp`, `python/tests/test_BSCurve.py`, `docs/sources/design/nlopt_removal.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Les deux derniers calculs confiés à NLopt, dans `gbs/bscanalysis.h`, passent
aux outils de Boost.Math (en-têtes seuls) : après cette PR, NLopt ne calcule
plus rien dans gbs. Il reste lié pour le type du paramètre `solver`, ignoré,
jusqu'à la PR 3.

## 2. Changements

### 2.1 `max_curvature_pos`

L'ancienne version minimisait 1/‖C''‖² avec COBYLA, mise au carré une fois de
plus par `solve_D_nlop`, et rendait ‖C''‖² au point trouvé : ni la position
du maximum de la courbure, ni sa valeur. Sur une ellipse de rayons 1 et 0,2,
elle rendait u = 0,434 et 1 654 ; la courbure y vaut 0, et son maximum est 25
aux sommets (u = 0,5).

- `curvature(crv, u)` (nouvelle) : ‖C' × C''‖ / ‖C'‖³ en toute dimension,
  calculée par √(‖C'‖²‖C''‖² − (C'·C'')²) / ‖C'‖³ ; 0 où C' s'annule.
- `max_curvature_pos` : le meilleur de 64 échantillons, puis Brent
  (`boost::math::tools::brent_find_minima`) entre ses voisins ; le nombre de
  bits vient de `tol_x` rapporté à l'intervalle, au plus la moitié de la
  mantisse (limite de Brent). Le résultat de Brent n'est gardé que s'il ne
  fait pas moins bien que l'échantillon (un maximum à un nœud, où C'' saute).
- Rend `{u, courbure}` : **rupture** du sens de la seconde valeur, annoncée
  dans `News.md`. `solver` est ignoré jusqu'à la PR 3.

### 2.2 `compute_rolling_ball`

Le cercle tangent à `crv1` en u₁ et à `crv2` : son centre est sur les deux
normales, à égales distances l₁ et l₂ des deux points de contact. L'ancienne
version minimisait |l₁² − l₂²| avec COBYLA depuis la projection de crv1(u₁)
sur crv2.

- L'écart des rayons |l₁| − |l₂| change de signe à la solution : ses
  changements de signe sont cherchés depuis la projection, par pas doublant
  des deux côtés, et le plus proche est résolu par **TOMS 748**
  (`boost::math::tools::toms748_solve`, encadrement garanti) ;
- un changement de signe à travers des normales parallèles (rayons infinis,
  pôle de l'écart) n'est pas une racine : la racine est vérifiée (rayons finis,
  écart sous 1e-8 relatif) ;
- sans changement de signe, Brent minimise |écart| sur tout l'intervalle ;
- le cas des normales parallèles à la projection (courbes parallèles) ne
  change pas.

### 2.3 `closest_curve`, `closest_p_curve`

Ils appelaient `midPoint()`, qui n'existe pas : ils ne compilaient pas une
fois instanciés, et n'avaient aucun appelant. Leur comparateur, qui projetait
deux fois chaque courbe, pouvait lever une exception dans un `min_element`
parallèle (`std::terminate`). Ils projettent maintenant le point **globalement**
(subdivision de la PR 1) sur chaque courbe, une fois, en parallèle par
`transform_rethrow`, et rendent l'itérateur de la plus proche.

### 2.4 En-têtes

`bscanalysis.h` n'inclut plus `nlopt.hpp`, mais `boost/math/tools/minima.hpp`
et `toms748_solve.hpp`.

## 3. Tests

| Test | Cas | Attendu |
|---|---|---|
| `max_curvature_pos` | ellipse de rayons 1 et 0,2, sur [0,1 ; 0,9] et sur tout le domaine ; cercle de rayon 2 ; parabole y = x² (B-spline quadratique exacte) | u = 0,5 et 25 ; un sommet et 25 ; 0,5 ; le sommet x = 0 et 2 |
| `compute_rolling_ball` | droites y = 0 et y = x, trois points de contact ; cercles concentriques ; cercle de rayon 1 et droite y = 3 | centre (a, a/(1 + √2)) et point de contact sur y = x à 1e-9 ; normales parallèles, centre au rayon 2 ; centre tangent aux deux, sur la normale du cercle |
| `closest_curve` | trois segments, un point ; en pointeurs, dans un autre ordre ; une droite non bornée dans la liste | la courbe du milieu ; `std::invalid_argument` remonté (pas de `std::terminate`) |
| Python `test_compute_rolling_ball` | droites y = 0 et y = x | centre attendu |

Comparaison à master, même programme compilé contre les deux versions :

| Cas | master (NLopt) | cette PR |
|---|---|---|
| boule entre y = 0 et y = x (3 points) | centre à 1e-9 … 4e-7, 14 à 58 µs | centre à 1e-16, 3 à 5 µs |
| boule entre un cercle et y = 3 (3 points) | tangence à 1e-6 … 7e-6, 15 à 31 µs | tangence à 6e-15, 5 à 7 µs |
| maximum de courbure de l'ellipse sur [0,1 ; 0,9] | u = 0,434, « 1 654 » | u = 0,5, 25 |

Les 491 tests ctest et les 54 tests Python passent (macOS arm64, clang
conda-forge, Release).
