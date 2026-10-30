# Palier 2 — PR 2 : intersection courbe-surface

| | |
|---|---|
| PR | branche `feat/intersect-curve-surface` |
| Document de référence | [brep_palier2.md](brep_palier2.md), § 2.2 et 2.4 |
| Fichiers | `gbs/intersect/curve_surface.h`, `gbs/intersect/surface_pieces.h`, `gbs/intersect/polynomial.h` (nouveaux) ; `gbs/intersect/box.h` (`segment_overlaps`), `gbs/intersect/curves.h` (test segment/boîte, écart de fusion) ; `tests/tests_intersect_curve_surface.cpp` (nouveau), `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Trouver tous les points où une courbe traverse une surface, et les portions de
la courbe posées sur la surface. C'est la brique :
- du lancer de rayon de la classification point/solide (PR 5) ;
- des départs du cheminement surface-surface, qui intersecte les bords de
  chaque domaine avec l'autre surface (PR 4) ;
- de la section, quand une arête perce une face (PR 6).

## 2. API

```cpp
template <T> struct CurveSurfaceHit     { T t, u, v; point<T, 3> p; bool tangent; };
template <T> struct CurveSurfacePortion { T ta, tb, ua, va, ub, vb; };   // courbe posée sur la surface
template <T> struct CurveSurfaceResult  { vector<CurveSurfaceHit> hits; vector<CurveSurfacePortion> portions; bool complete; };

auto intersect(c, a, b, s, box /* u1, u2, v1, v2 */, IntersectOptions = {}) -> CurveSurfaceResult;
auto intersect(c, s, IntersectOptions = {}) -> CurveSurfaceResult;        // bornes de la courbe et de la surface

auto polynomial_roots(coefs, a, b) -> PolynomialRoots { roots, extremes };   // degré quelconque, sur [a, b]
bool segment_overlaps(a, b, box, gap);                                         // segment contre boîte (dalles)
```

Mêmes options et mêmes conventions que la courbe-courbe : tolérance, sinus de
tangence, plafond de paires.

## 3. Choix d'architecture et justification

### 3.1 Formules fermées

| Courbe | Surface | Équation en t |
|---|---|---|
| droite (`LineCurve`, ou B-spline de degré 1 à deux pôles) | plan | degré 1 |
| droite | cylindre, cône (les deux nappes), sphère | degré 2 |
| droite | tore | degré 4 : (\|x\|² + R² − r²)² = 4R²ρ² |
| cercle ou ellipse (`EllipseCurve`) | plan | d0 + d1 cos t + d2 sin t = 0 |

- **Racines.** `polynomial_roots` isole les racines entre celles de la dérivée
  (récursivement), puis les affine par **bichotomie**. C'est sûr pour le petit
  degré, sans formule de Ferrari ni matrice compagnon.
- **Racines doubles.** Elles ne changent pas de signe : ce sont des extrema,
  rendus à part. Chaque candidat est vérifié géométriquement, par l'inverse en
  formule fermée de la primitive (`parameters`), puis poli par Gauss-Newton.
  Les tangences sont ainsi trouvées sans seuil algébrique.
- **Domaine.** Les paramètres (u, v) sont ramenés dans le tour de la boîte
  demandée, et filtrés par cette boîte.
- **Équation nulle** (droite posée sur le plan ou sur une génératrice, cercle
  dans le plan) : on passe à la subdivision, qui sait rendre les portions et
  leur découpe par la boîte.

### 3.2 Subdivision

| Surface | Carreaux | Boîte | Moitiés |
|---|---|---|---|
| `BSSurface`, `BSSurfaceRational` | **carreaux de Bézier** (restriction à la boîte, puis nœuds insérés jusqu'au degré dans les deux directions) | pôles cartésiens : **garantie** | de Casteljau sur les lignes ou les colonnes |
| autres (primitives, révolution, offset…) | un carreau échantillonné (6 × 6 cellules) | échantillons élargis de deux fois l'écart du centre de chaque cellule à la moyenne de ses coins : **heuristique** | nouvel échantillonnage |

- **Planéité** : distance maximale des pôles (ou échantillons) au plan passant
  par le centre, de normale le produit des diagonales.
- **Direction de découpe** : celle où le carreau est le plus long, mesurée le
  long des **lignes de pôles**. Sur un carreau fermé (cylindre complet), les
  coins coïncident, et une mesure bout à bout coupait toujours dans l'autre
  sens.
- **Candidat** : la corde du morceau de courbe contre le plan tangent au centre
  du carreau, puis Gauss-Newton sur (t, u, v) avec le système normal 3 × 3.
  Près d'une tangence, on projette alternativement sur la surface et sur la
  courbe.
- **Portion** : un morceau de courbe dont cinq points se projettent sur le
  carreau dans la tolérance. Les tangences sont distinguées par le profil de
  distance, comme pour la courbe-courbe.

### 3.3 Segment contre boîte

L'axe d'une boîte englobante ne suit pas les objets. La boîte d'un long segment
incliné couvre un grand volume, et recoupe beaucoup de carreaux qui en sont
loin. Le cas réel est une droite qui traverse un cylindre sur un repère
incliné : il dépassait les 200 000 paires.

Un morceau plat est donc testé comme un **tube autour de sa corde** : le
segment contre la boîte de l'autre, élargie du rayon du tube (méthode des
dalles). Le test est prudent : il ne rejette jamais une paire qui se touche.
Il profite aussi à la courbe-courbe (PR 1).

### 3.4 Correction de la PR 1

Les graines de recouvrement laissent des trous au bord des carreaux, et ces
trous atteignent deux fois la taille minimale d'un morceau. L'écart de fusion
passe de 20 tol à **50 tol**, dans les deux algorithmes. Sans cela, une
génératrice posée sur un cylindre donnait 65 portions au lieu d'une.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| boîtes des primitives en formule fermée | carreaux échantillonnés | même code pour toutes les surfaces non B-spline ; les primitives ont leurs formules fermées pour les droites et les coniques |
| cercle contre plan | cercle **et ellipse** contre plan | même équation |
| — | `polynomial_roots` (dérivées et bichotomie) | plus sûr que Ferrari pour le tore |
| — | `segment_overlaps` | explosion du nombre de paires sur les repères inclinés |

## 5. Tests (`tests_intersect_curve_surface`, 6 cas, 407 assertions, 0,3 s)

| Test | Couvre |
|---|---|
| `polynomial_roots` | cubique à trois racines, sous-intervalle, racine double presque nulle (extremum) et exacte, quartique, coefficients de tête nuls |
| `segment_and_box` | traversée, segment intérieur, raté, raté sauvé par l'écart, diagonale passant à côté d'un coin, boîte vide |
| `lines_and_primitives` | sur un repère incliné : droite contre cylindre, cône, sphère, tore (4 points exacts), plan parallèle ; **formule fermée comparée à la subdivision** (même droite en B-spline à nœud intérieur) ; plan oblique, point hors de la boîte ; droites **tangentes** au cylindre et à la sphère |
| `curves_on_surfaces` | génératrice d'un cylindre (portion bornée par la boîte en v), droite dans un plan bornée par la boîte, parallèle d'un cylindre (portion sur le tour, à travers le seam) |
| `circles_and_planes` | cercle coupant un plan en son point de départ (une solution) et en π ; même cercle en NURBS exacte ; ellipse qui le manque ; cercle tangent ; cercle dans le plan (portion) |
| `bsplines` | sinusoïde cubique contre un demi-cylindre, primitive et NURBS (5 points) ; boucle B-spline contre un tore incliné, primitive et NURBS (4 points) ; droite contre une surface B-spline interpolée |

Les suites `tests_intersect_curves` et `tests_intersect_curve_surface` passent.

## 6. Points à valider

1. **Formules fermées** pour droite/quadriques, droite/tore, conique/plan, avec
   les racines par dérivées et bichotomie, et une vérification géométrique.
2. **Équation nulle → subdivision**, qui rend les portions bornées par la
   boîte.
3. **Carreaux de Bézier** pour les B-splines, **échantillonnés** pour le reste.
4. **Segment contre boîte** pour les morceaux plats.
5. **Écart de fusion** porté à 50 tol (correction de la PR 1).

## 7. Revue de code (a posteriori, lot 2 de la revue du noyau BREP)

| Commentaire | Suite donnée |
|---|---|
| 1. Une courbe posée sur une surface courbe ou inclinée fait exploser la subdivision : une graine n'est semée que si la boîte d'un seul carreau contient le morceau, les paires ne sont rejetées que par des boîtes alignées sur les axes. Génératrice d'un cylindre : `complete` faux et portion tronquée [−1,35 ; 5] au lieu de [−5 ; 5] ; cylindre primitive incliné : 2 816 points ; cubique à 1e-4 d'un plan B-spline incliné : `complete` faux. `classify.h` ignore `complete` | Paires rejetées aussi par la dalle du plan du carreau (`lower_bound` de `projection.h`) ; graine sans test d'inclusion, ses points projetés sur toute la boîte depuis le carreau ; paires dont le morceau est déjà dans une graine sautées. Génératrice de 40 : 0,04 s au lieu de 88 s sans aboutir. `classify` ne compte jamais un rayon dont une intersection (courbe-surface ou courbe-courbe) est plafonnée : autre direction, comme pour une ambiguïté ; si tous les rayons le sont, `BRepError` (`nullopt` pour une face contre un solide, déjà traité comme un échec par les booléens). `interior_points` saute une ligne plafonnée. Testé (`lying_curves_on_tilted_surfaces` ; `capped_intersections` dans `tests_brep_classify`, avec une face et une pcourbe rugueuses à toute échelle) |
| 2 et 9. Contacts tangents et rasants par subdivision : chaque fragment passait le test de profil (cercle NURBS tangent à un plan : 1 point et 14 portions ; à 5e-7 du plan : 69 points tangents) ; une courbe à 0,4 tol d'un plan sur 8 unités, qui le touche en x = kπ : 1 point au lieu d'une portion | Les graines se fusionnent quand la courbe reste dans la tolérance entre elles ; chaque tronçon est prolongé jusqu'où la courbe quitte la tolérance (ou jusqu'au bout de son intervalle : bout ouvert). Il est un contact quand la distance ne reste près de son minimum (à moins d'un huitième du chemin vers la tolérance) que sur moins de la moitié de sa longueur : parabole d'une tangence, V d'une traversée rasante, W d'une courbe qui plonge sous la surface de moins que la tolérance ; une seule traversée alors, au minimum poli par `newton_curve_surface`. Sinon, ou s'il est ouvert aux deux bouts, c'est une portion. Une tangence et les traversées voisines, la courbe dans la tolérance entre elles, ne font qu'un point, tangent. Testé (`tangent_and_grazing_contacts`) |
| 3. Formules fermées dépendantes de l'origine du paramètre de la droite (règle d'abandon des coefficients de tête comparant des puissances différentes, annulations) : rayons partant à 1e4 d'un tore, 182 sur 300 faux (parité de `classify` fausse dans un modèle de 10 m en mm) ; droite de nœuds [1e7, 1e7 + 2] à travers un cylindre unité : aucun point | Droite reparamétrée par l'abscisse curviligne depuis le pied de l'origine du repère ; `polynomial_roots` compare chaque terme \|c_k\| h^k au plus grand max \|c_i\| h^i, h le plus grand \|t\| de l'intervalle. Testé (`lines_far_from_the_origin`) |
| 4. Formules fermées sans tolérance (\|coef\| ≤ 1e-12 absolu ; ellipse contre plan : rien quand \|x\| > 1 + 1e-12) : droite à 1e-9 ou 1e-7 d'un plan (tol 1e-6) : rien ; génératrice d'un cylindre de rayon 100, bâtie par `iso_u`, non reconnue 157 fois sur 200 | « Sur la surface » décidé géométriquement : distance exacte à la primitive en cinq points de la droite (≤ tol : subdivision, qui rend la portion) ; cercle dans un plan si \|d0\| + ρ ≤ tol ; extrema φ et φ + π toujours candidats (contact à moins de la tolérance). Testé (`closed_forms_within_tolerance`) |
| 5. Formules fermées : traversées aux bouts de [a, b] perdues dans 40 % des cas (borne stricte de la racine linéaire, signe arrondi aux bouts), et la subdivision ne prend pas le relais (vecteur vide, pas `nullopt`) | Les bouts a et b sont candidats dès que la courbe y est dans la tolérance ; l'affinage et le filtre décident. Testé (`crossings_at_range_ends`) |
| 6. Gauss-Newton s'arrête au pôle d'une sphère NURBS (pas ramené sur v = π/2, Su = 0, déterminant nul) : 4 traversées proches du pôle sur 300 perdues ; une projection s'y arrête aussi | `refine_cs_restart` : arrêté hors de la tolérance sur un point singulier du paramétrage, l'affinage repart de 8 valeurs du paramètre qui s'annule, l'autre écarté de la ligne singulière d'environ la distance restante, des deux côtés. Testé (`crossings_near_poles`) |
| 7. Une courbe fermée qui traverse la surface en son point de départ y est rendue deux fois (test exact `t == a` et `t == b`) : 314 à 496 cas sur 500 | Les points sont comparés : premier et dernier au même point, au départ de la courbe fermée. Testé (`closed_curve_crossing_at_start`) |
| 8. Au pôle ou au sommet, `parameters()` donne un u arbitraire que le filtre de boîte rejette avant l'affinage ; sur l'autre nappe d'un cône u vaut θ + π : axe d'un fuseau u ∈ [1, 2] par ses pôles, aucun point ; droite sur la nappe v < 0, 1 point au lieu de 2 | Paramètres ramenés dans la boîte au lieu d'être filtrés : l'affinage et la tolérance décident ; sur un cône, de u et u + π, le plus proche (la correction de `ConicalSurface::parameters` est faite à part, dans `primitives.h`). Testé (`poles_apexes_and_nappes`) |

Les sondes de la revue donnent toutes le résultat attendu après correction (génératrices : 1 portion [−5 ; 5], `complete` vrai ; contacts : 1 point tangent, aucune portion ; bouts : 0 perdu sur 2 000 ; pôles : 0 sur 300 ; tore lointain : 0 sur 300). Durées : `tests_brep_boolean` 24,4 s avant, 24,9 s après ; `tests_brep_classify` 0,19 s, 0,66 s (l'essentiel pour le nouveau test, dont chaque rayon plafonné parcourt 200 000 paires) ; `tests_intersect_curve_surface` 0,35 s, 0,27 s avec huit cas de plus.
