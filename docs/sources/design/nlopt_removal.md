# Suppression de la dépendance à NLopt — conception

| | |
|---|---|
| Statut | **À valider** |
| Périmètre | Remplacer les appels à NLopt (`gbs/solvers.h`, `gbs/extrema.h`, `gbs/bscanalysis.h`) par des méthodes de Newton écrites dans gbs, corriger les défauts trouvés en chemin, puis retirer NLopt du build, de la CI et des bindings Python |
| Dépend de | les projections de Newton de `gbs/intersect` (palier 2 du noyau BREP) |
| Hors périmètre | un optimiseur générique sous contraintes ; l'unification des raffinements `intersect_detail::refine` / `refine_cs` des intersections avec les extrema (possible plus tard, voir § 8) |

Ce document décrit *quoi* construire et *pourquoi*. Chaque PR a sa note
d'architecture `nlopt_prNN_*.md`.

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Algorithme des extrema | **Newton amorti** sur ½‖r‖², hessienne complète (termes de courbure compris), repli **Gauss-Newton** quand elle n'est pas définie positive, **recherche linéaire** (la distance ne remonte jamais), **bornes actives** figées | COBYLA (NLopt), sans dérivées, lent, distance rendue au carré |
| Où vit le code | nouvel en-tête `gbs/intersect/projection.h`, en dimension quelconque ; `intersect_detail::project` et `project_on_surface` l'appellent | une copie de plus dans `extrema.h` |
| Amorçage | grille d'échantillons comme aujourd'hui (30, 30 × 30, 30³), Newton depuis les **deux** meilleurs échantillons, le meilleur résultat gagne | un seul départ (une couture de courbe fermée peut piéger le départ contre la mauvaise borne) |
| Distance rendue | la **vraie distance** ‖r‖ | d² rendu aujourd'hui par erreur |
| Paramètre `solver` (C++) | gardé et **ignoré** en PR 1 et 2, retiré en PR 3 | rupture d'API immédiate |
| Python `solver=`, enum `nlopt_algorithm` | gardés en PR 1 et 2 (`DeprecationWarning` si `solver` est passé), retirés en PR 3 | retrait immédiat |
| Dépendances | supprimer les dépendances **à l'exécution** (bibliothèques à installer et lier) ; les dépendances **à la compilation** en en-têtes seuls sont permises (Boost, Eigen) | réécrire ce que les en-têtes de Boost fournissent |
| Extrema 1D de `bscanalysis.h` | **Brent** (section dorée et interpolation parabolique, `boost::math::tools::brent_find_minima`, déjà utilisé dans `extrema.h`) ; racines encadrées par **TOMS 748** (`boost::math::tools::toms748_solve`) | Newton 1D (dérivées d'ordre 3 ou 4 nécessaires pour la courbure) |
| Newton à plusieurs variables | écrit dans gbs (`projection.h`) : Boost.Math n'a que des outils à une variable (racines, minima) et, pour plusieurs variables, des optimiseurs globaux sans dérivées (évolution différentielle, CMA-ES…), sans bornes actives ni convergence quadratique | Boost.Math `optimization`, Eigen `LevenbergMarquardt` (module non supporté, sans bornes) |
| `solve_N_nlop`, `solve_D_nlop` | supprimés en PR 3 avec `gbs/solvers.h` et leurs tests | en-tête optionnel |

---

## 1. Objectifs et non-objectifs

### 1.1 Objectifs

1. Ne plus dépendre de NLopt : une dépendance de moins à installer (conda,
   vcpkg), un souci de résolution conda en moins sur osx-arm64 (commentaire de
   `ci.yml`), un en-tête lourd de moins dans les unités de compilation. Le
   but est de supprimer les dépendances **à l'exécution** ; les bibliothèques
   en en-têtes seuls restent permises à la compilation (Boost, déjà utilisé
   pour `brent_find_minima`, `newton_raphson_iterate` et les quadratures, et
   Eigen).
2. Des extrema **justes** : vraie distance, convergence quadratique au lieu des
   centaines d'évaluations de COBYLA, comportement défini aux bornes et aux
   coutures.
3. Une seule implémentation des projections de Newton, partagée par les
   intersections (`gbs/intersect`) et les extrema.

### 1.2 Non-objectifs

- Garder un optimiseur générique (`solve_N_nlop` sans dérivées) : il n'a plus
  d'appelant après la PR 2.
- Trouver **tous** les extrema locaux : les fonctions `extrema_*` rendent le
  minimum de la distance trouvé à partir d'une grille, comme aujourd'hui. Les
  intersections complètes restent le rôle de `gbs::intersect`.

---

## 2. État des lieux

### 2.1 Utilisations de NLopt

| Lieu | Usage | Dimension |
|---|---|---|
| `gbs/solvers.h` | `solve_N_nlop` (sans dérivées), `solve_D_nlop` (somme des carrés des résidus) ; seulement bornes et `xtol_rel`, pas de `maxeval` | n |
| `gbs/extrema.h` | `extrema_curve_point`, `extrema_surf_pnt`, `extrema_curve_curve`, `extrema_curve_curve_lst`, `extrema_surf_curve`, avec `default_nlopt_algo = LN_COBYLA` | 1, 2, 2, 3 |
| `gbs/bscanalysis.h` | `max_curvature_pos` (`solve_D_nlop`), `compute_rolling_ball` (`solve_N_nlop`) | 1 |
| `python/gbsbindextrema.h`, `python/gbsbind.cpp` | argument `solver` de tous les extrema, enum `nlopt_algorithm` | — |
| build | `find_package(NLopt REQUIRED)`, liens vers `NLopt::nlopt` (cible `gbs`, PCH, Python, tests, OCCT, exemples de rendu, `tests_3rdlibs`), `GbsConfig.cmake.in`, paquet conda en CI, ReadMe | — |

Appelants des extrema : `dev_from_points`, `closest_curve`,
`closest_p_curve`, `compute_rolling_ball` (`bscanalysis.h`) ;
`get_tg_from_spine`, `loft_approx`, `build_intersections` donc `gordon`
(`bssbuild.h`) ; `project_points_on_surface`, `lattice_intersections`
(`gbs-mesh/tfi.h`) ; `detail::project_point_global` (`gbs-brep/pcurve.h`) ;
l'exemple `gbs-render/examples/vtk_render.cpp`.

### 2.2 Défauts constatés

1. **Distance au carré.** L'objectif rend `‖r‖²`, `solve_D_nlop` le met encore
   au carré (‖r‖⁴ minimisé) et les fonctions rendent `sqrt(minf) = ‖r‖²`. Les
   tests ne vérifient que des distances nulles. Conséquences :
   - `build_intersections` (`bssbuild.h:491`) et `lattice_intersections`
     (`tfi.h:743, 756`) comparent d² à `tol` : deux courbes à 1e-4 l'une de
     l'autre passent le contrôle à 1e-6 (1e-8 < 1e-6) ;
   - `detail::project_point_global` (`pcurve.h:107`) compare la distance de
     Newton à d² : sous 1, d² est plus petit, la solution NLopt (moins précise)
     est retenue à tort et son « écart » sous-estimé ;
   - `dev_from_points` rend des écarts au carré.
2. **Gradient surface-courbe faux.** `srf(x[1],1,0)` et `srf(x[2],0,1)`
   (`extrema.h:286-287`) appellent `value(u = x[i], v = 1, du = 0)` : sans
   effet avec COBYLA (sans dérivées), faux avec un algorithme `LD_*`.
3. **`extrema_curve_curve_lst`** ignore `n_bracket_u` / `n_bracket_v`.
4. **`max_curvature_pos`** minimise 1/‖C''‖², pas la courbure ; même mise au
   carré, la valeur rendue est ‖C''‖² (PR 2).
5. **Arrêt** : NLopt s'arrête sur `xtol_rel = tol_x²` sans nombre maximal
   d'évaluations ; COBYLA peut en faire des centaines.

---

## 3. Algorithme commun

### 3.1 Formulation

Pour chaque problème, un résidu r(x) ∈ ℝ^dim de n paramètres x (n = 1 à 3) :

| Problème | x | r(x) | Jacobien J |
|---|---|---|---|
| point-courbe | t | C(t) − P | [C'] |
| point-surface | (u, v) | S(u, v) − P | [Sᵤ, Sᵥ] |
| courbe-courbe | (s, t) | C₁(s) − C₂(t) | [C₁', −C₂'] |
| courbe-surface | (t, u, v) | C(t) − S(u, v) | [C', −Sᵤ, −Sᵥ] |

On minimise f(x) = ½‖r(x)‖² sur la boîte des bornes. Gradient g = Jᵀr ;
hessienne H = JᵀJ + Σₖ rₖ ∇²rₖ, soit pour la surface
H = [[Sᵤ·Sᵤ + r·Sᵤᵤ, Sᵤ·Sᵥ + r·Sᵤᵥ], [·, Sᵥ·Sᵥ + r·Sᵥᵥ]], et de même pour les
autres (signe moins devant les dérivées secondes du terme soustrait).

Les termes de courbure r·∇²r rendent la convergence **quadratique même quand
la distance n'est pas nulle**. Gauss-Newton (JᵀJ seul), utilisé par les
intersections où le résidu est nul à la solution, n'y converge que
linéairement, d'autant plus lentement que la distance approche du rayon de
courbure.

### 3.2 Itération

À partir d'un départ x₀ dans la boîte :

1. **Bornes actives** : une variable sur sa borne dont le gradient pousse vers
   l'extérieur est figée pour ce pas (ligne et colonne remplacées par
   l'identité, composante du gradient nulle). Sans cela, le pas couplé est
   rogné sur la variable bornée et l'autre variable converge mal (point
   au-delà d'un bord de surface). Une variable sur sa borne dont le **pas** de
   Newton sort de la boîte (par le couplage de H, gradient vers l'intérieur)
   est figée elle aussi et le pas recalculé : rogné, il serait annulé.
2. **Pas de Newton** : résolution de H Δx = −g par Cholesky (n ≤ 3). Si H n'est
   pas définie positive (loin de la courbe, paramétrage rationnel à vitesse
   très variable), ou si la géométrie ne fournit pas ses dérivées secondes
   (`has_second_derivatives()` faux : surface décalée), pas de
   **Gauss-Newton** JᵀJ. Si JᵀJ est singulière (tangentes parallèles en
   courbe-courbe), **Levenberg-Marquardt** JᵀJ + μ diag(JᵀJ), μ = 1e-6 ; en
   dernier recours (dérivée nulle), pas diagonal −gᵢ/(JᵀJ)ᵢᵢ.
3. **Recherche linéaire** : x + αΔx ramené dans la boîte, α = 1, ½, ¼… jusqu'à
   ce que f décroisse, ou reste égale sur un petit pas (au plus 30
   réductions) : un grand pas à f égale ferait permuter sans fin entre deux
   minima égaux (segments parallèles). Si x ne peut pas bouger alors que le
   gradient projeté n'est pas nul, une **descente de gradient** mise à
   l'échelle (−gᵢ/(JᵀJ)ᵢᵢ) est essayée ; sinon arrêt, x est le minimum à la
   précision machine.
4. **Arrêt** : chaque |Δxᵢ| ≤ τ (1 + |xᵢ|) avec τ = max(tol_x², 4 ε), ce qui
   reprend le `xtol_rel = tol_x²` de l'appel NLopt actuel ; ou f = 0 ; ou 100
   itérations. En pratique, 3 à 6 itérations.
5. **Coutures** : arrêté sur une borne où la borne opposée donne la même
   distance (couture d'une courbe ou d'une surface fermée), Newton repart de
   l'autre côté et garde le meilleur résultat, au plus n fois. Sur un tore,
   les quatre coins (0 ou 2π, 0 ou 2π) du domaine sont le même point : les deux
   meilleurs échantillons peuvent tomber du même côté de la couture en v, et
   Newton reste alors contre la borne.

La distance rendue est ‖r(x)‖ au point final.

### 3.3 Amorçage

Les surcharges sans départ échantillonnent les bornes comme aujourd'hui
(`n_bracket` 30 par paramètre ; 30³ pour courbe-surface) et lancent Newton
depuis les **deux** échantillons les plus proches ; le meilleur résultat est
rendu. Sur une courbe ou une surface fermée, un point proche de la couture a
ses deux meilleurs échantillons aux deux bouts du domaine ; un seul départ
peut rester contre la mauvaise borne. Les intersections font déjà ainsi
(correctif de `intersect_detail::project`, palier 2 bis). Le second départ
doit être à une autre distance que le premier : les échantillons d'un bord
dégénéré (pôle d'une sphère) sont tous le même point, où Newton reste bloqué
(dérivée nulle) ; les copies d'un point de part et d'autre d'une couture sont
laissées à la reprise du § 3.2.

Les surcharges avec départ (`u0`, `(u0, v0)`…) lancent Newton depuis ce départ
seul : c'est un raffinement local, comme aujourd'hui.

Un domaine non borné (droite `Line`) ne peut pas être échantillonné :
`check_bracket_bounds` lève `std::invalid_argument`, comme aujourd'hui pour
les courbes ; la vérification est étendue aux surfaces.

### 3.4 Problèmes 1D de `bscanalysis.h` (PR 2)

- `max_curvature_pos` : maximum de la **courbure** κ = ‖C' × C''‖ / ‖C'‖³
  (formule en dimension quelconque : √(‖C'‖²‖C''‖² − (C'·C'')²) / ‖C'‖³) par
  Brent sur −κ, depuis le meilleur de quelques échantillons pour éviter un
  maximum local ; la valeur rendue est κ (aujourd'hui ‖C''‖², ce qui n'est la
  courbure que pour un paramétrage par l'abscisse curviligne). Changement de
  sens à annoncer dans `News.md`.
- `compute_rolling_ball` : la différence des rayons
  ‖L₁(u₁') − L₁(0)‖ − ‖L₂(u₂') − L₂(0)‖ change de signe au contact de la
  boule : racine encadrée par TOMS 748 (`toms748_solve`, en-têtes Boost) sur
  l'intervalle de `crv2` autour du départ donné par la projection, plutôt
  qu'un minimum de sa valeur absolue (non dérivable à la racine, moins
  précis). À défaut d'encadrement, Brent sur la valeur absolue.

---

## 4. API

### 4.1 `gbs/intersect/projection.h` (nouveau)

En dimension quelconque, `T` flottant :

```cpp
// raffinement local depuis un départ ; bornes = intervalle ou boîte de paramètres
auto newton_point_curve(const Curve<T, dim> &c, const point<T, dim> &p, T t0,
                        std::array<T, 2> range, T tol_x) -> std::array<T, 2>;       // {t, d}
auto newton_point_surface(const Surface<T, dim> &s, const point<T, dim> &p, T u0, T v0,
                          std::array<T, 4> box, T tol_x) -> std::array<T, 3>;       // {u, v, d}
auto newton_curve_curve(const Curve<T, dim> &c1, const Curve<T, dim> &c2, T s0, T t0,
                        std::array<T, 2> range1, std::array<T, 2> range2, T tol_x)
    -> std::array<T, 3>;                                                             // {s, t, d}
auto newton_curve_surface(const Curve<T, dim> &c, const Surface<T, dim> &s, T t0, T u0, T v0,
                          std::array<T, 2> range, std::array<T, 4> box, T tol_x)
    -> std::array<T, 4>;                                                             // {t, u, v, d}

// global : grille d'échantillons, Newton depuis les deux meilleurs
auto project_point_curve(c, p, range, n, tol_x);
auto project_point_surface(s, p, box, nu, nv, tol_x);
auto closest_curve_curve(c1, c2, range1, range2, n1, n2, tol_x);
auto closest_curve_surface(c, s, range, box, n, nu, nv, tol_x);
```

Le moteur commun est `projection_detail::minimize<T, n>` (§ 3.2), qui reçoit
une fonction rendant (f, g, H, JᵀJ) et une fonction rendant f seul (recherche
linéaire).

`intersect_detail::project` (curves.h) et `intersect_detail::project_on_surface`
(curve_surface.h) deviennent des appels à `project_point_curve` /
`project_point_surface` avec leurs grilles actuelles (17 et 5 × 5 points).

### 4.2 `gbs/extrema.h`

Signatures et types de retour **inchangés** (`std::array` : `{u, d}`,
`{u, v, d}`, `{u1, u2, d}`, `{u, v, uc, d}`) ; seul `d` change de sens : la
vraie distance. Le paramètre `solver` est ignoré (documenté comme tel) jusqu'à
la PR 3. `tol_x` garde son rôle de tolérance sur les paramètres (§ 3.2).

`extrema_curve_curve(Line<T, 2>, Line<T, 2>)` (forme fermée) ne change pas.

---

## 5. Découpage en PR

| PR | Contenu | Note |
|---|---|---|
| 1 | `gbs/intersect/projection.h` ; `extrema.h` réécrit dessus (NLopt n'y sert plus qu'au type du paramètre ignoré) ; `intersect_detail::project` et `project_on_surface` branchés dessus ; correction de la distance au carré et des seuils de ses appelants ; `extrema_curve_curve_lst` passe ses `n_bracket` ; tests ; ce document | `nlopt_pr01_extrema.md` |
| 2 | `bscanalysis.h` : `max_curvature_pos` (vraie courbure, Brent), `compute_rolling_ball` (Brent) ; `closest_curve` / `closest_p_curve` réécrits (ils appellent `midPoint()`, qui n'existe pas, et leur comparateur peut lever une exception dans un `min_element` parallèle) ; tests | `nlopt_pr02_analysis.md` |
| 3 | retrait de NLopt : `gbs/solvers.h`, `default_nlopt_algo`, paramètres `solver`, enum Python `nlopt_algorithm`, `tests_solvers.cpp`, `tests_nlopt_example.cpp`, CMake (racine, PCH, Python, tests, OCCT, rendu, `tests_3rdlibs`, `GbsConfig.cmake.in`), CI (trois jobs), ReadMe, `News.md` | `nlopt_pr03_removal.md` |

Après la PR 1, `gbs/extrema.h` inclut encore `gbs/solvers.h` pour le type
`nlopt::algorithm` ; après la PR 2, plus aucun calcul n'utilise NLopt ; la PR 3
ne touche plus que des déclarations, le build et la documentation.

---

## 6. Tests

`tests/tests_extrema.cpp`, en `double` et en `float` :

1. **vraies distances** non nulles et connues : point à distance h d'un
   cercle (rationnel), d'un cylindre, d'une sphère ; droites gauches ;
   courbe parallèle à un plan ; on vérifie d = h et non h² ;
2. **bornes** : projection hors de l'intervalle (le point rendu est
   l'extrémité, d est la distance à l'extrémité) ; point au-delà d'un bord de
   surface (u sur la borne **et v exact**, cas des bornes actives) ;
3. **coutures** : point près de la couture d'un cercle complet, d'un cylindre
   et d'un tore complets, des deux côtés ;
4. **raffinement local** depuis un départ, convergence en peu d'itérations ;
5. le paramètre `solver` est accepté et sans effet ;
6. `extrema_curve_curve_lst` avec `n_bracket` ;
7. les suites existantes (`intersect`, `bssbuild`, `mesh`, `bscapprox`,
   `brep`, Python) passent ; `build_intersections` lève bien sur des courbes à
   1e-4 l'une de l'autre avec `tol` = 1e-6 (régression du seuil).

---

## 7. Risques

| Risque | Parade |
|---|---|
| Seuils des appelants : avec la vraie distance, `build_intersections` et `lattice_intersections` deviennent plus stricts (d > d² pour d < 1) | décision : on garde le seuil strict, `tol` est une distance partout comme dans la grille de `gordon` ; rupture annoncée dans `News.md`, tests plus sélectifs |
| Changement de résultat des intersections (`project`, `project_on_surface` branchés sur le nouveau moteur) | suites `intersect`, `brep`, `step` et OCCT relancées |
| Newton converge vers un autre minimum local que COBYLA | même grille, deux départs au lieu d'un : au moins aussi bon que le meilleur échantillon |
| `float` : τ ne peut descendre sous la précision | plancher 4 ε dans le critère d'arrêt |

---

## 8. Suites possibles

- `intersect_detail::refine` et `refine_cs` (Gauss-Newton des intersections,
  arrêt sur `tol`) pourraient utiliser `minimize` ; ils sont laissés tels quels
  ici pour ne pas changer les intersections au-delà de la projection.
- Une variante qui rend tous les minima locaux à partir de la grille (au lieu
  du meilleur) si un appelant en a besoin.
