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
