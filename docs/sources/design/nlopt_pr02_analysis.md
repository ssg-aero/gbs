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
- `max_curvature_pos` : le meilleur d'échantillons pris sur chaque tronçon
  de nœuds d'une B-spline (8 par tronçon, au moins 64 ; 64 pour une autre
  courbe, par `CurvePieces`), puis Brent
  (`boost::math::tools::brent_find_minima`) entre ses voisins ; le nombre de
  bits vient de `tol_x` rapporté à l'intervalle, au plus la moitié de la
  mantisse (limite de Brent). Le résultat de Brent n'est gardé que s'il ne
  fait pas moins bien que l'échantillon (un maximum à un nœud, où C'' saute).
- `std::invalid_argument` sauf si u1 < u2 (intervalle fini).
- Rend `{u, courbure}` : **rupture** du sens de la seconde valeur, annoncée
  dans `News.md`. `solver` est ignoré jusqu'à la PR 3.
- `curvature` : norme du produit vectoriel exacte en 2D (déterminant) et en
  3D ; l'identité de Lagrange, qui s'effondre quand C' et C'' sont presque
  parallèles, ne sert qu'en dimension supérieure.

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
  pôle de l'écart) n'est pas une racine : l'encadrement final [lo, hi] n'est
  accepté que si le sinus de l'angle des normales y garde son signe (pas de
  pôle dedans) et que l'écart y est petit relativement aux rayons (ce dernier
  test n'écarte que le saut d'une normale à un nœud de multiplicité pleine) ;
- `crv2` non bornée (`Line`) : `std::invalid_argument` (`check_bracket_bounds`),
  la recherche porte sur son intervalle ;
- sans changement de signe valide : `std::runtime_error`, aucun cercle n'est
  tangent aux deux courbes ;
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
et `toms748_solve.hpp`, et `gbs/intersect/curve_pieces.h` directement :
`extrema.h` inclut `<gbs/curves>`, donc `curveoffset.h`, donc
`bscanalysis.h`, qui ne voit alors pas encore ce que `extrema.h` déclare.

### 2.5 Revue de code

| Commentaire | Suite donnée |
|---|---|
| `crv2` non bornée : Brent sur [−∞, ∞] ne s'arrête jamais | `check_bracket_bounds(crv2.bounds())` en tête ; `News.md` (passer un segment). La fenêtre finie de la revue (rayon ≤ d / (1 − cos θ)) est laissée pour plus tard |
| la tolérance 1e-8 (1 + l₁) grandit avec les rayons qu'elle doit écarter : le saut au pôle passe pour une racine | continuité de l'encadrement : le sinus des normales garde son signe sur [lo, hi] ; l'écart relatif ne sert plus qu'au nœud de multiplicité pleine. Cercles face à face : **avant, jusqu'à 10 points sur 39 donnaient un centre à environ 1e16**, plus aucun |
| 64 échantillons fixes : un pic plus étroit que l'intervalle / 64 est manqué | 8 échantillons par tronçon de nœuds (`CurvePieces`), au moins 64 ; écrit avec `std::views::transform`, `std::ranges::to` et `std::ranges::max_element`. B-spline de 300 tronçons à coin arrondi sur 2 : **avant, courbure rendue 2e-14 au lieu de 2,83** |
| intervalle non vérifié (u1 == u2 : conversion de l'infini en entier) | `std::invalid_argument` sauf si u1 < u2 et l'intervalle fini. La convention citée (`arc_length_distrib_params` lève `std::length_error` pour u1 ≥ u2) n'existe pas : ce fichier lève `std::invalid_argument` pour un argument invalide |
| sans changement de signe, l'approche la plus proche rendue comme une solution | décision : `std::runtime_error` (le plus strict) ; rupture dans `News.md` |
| l'identité de Lagrange s'effondre aux faibles courbures | déterminant en 2D, `norm(cross)` en 3D, identité au-delà ; Bézier presque droite (courbure 4e-10) : **avant 0**, maintenant exact à 5e-16 relatif |
| `closest_curve` et `closest_p_curve` en double | `analysis_detail::closest` avec une projection (identité ou déréférencement), `std::ranges` |
| tests manquants | `compute_rolling_ball_review` (cercles face à face, `Line`, aucune boule), `max_curvature_pos_narrow_peak` (et intervalles vides ou inversés), `curvature_exact`, `closest_curve` sur une plage vide (rend la fin) |

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

Après la revue, les 494 tests ctest et les 54 tests Python passent (macOS
arm64, clang conda-forge, Release).
