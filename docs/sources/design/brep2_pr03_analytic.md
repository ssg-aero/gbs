# Palier 2 — PR 3 : surface-surface, paires de primitives

| | |
|---|---|
| PR | branche `feat/intersect-analytic` |
| Document de référence | [brep_palier2.md](brep_palier2.md), § 2.1, 2.5 et 2.7 |
| Fichiers | `gbs/intersect/surface_surface.h`, `gbs/intersect/analytic_surfaces.h` (nouveaux), `tests/tests_intersect_surfaces.cpp` (nouveau), `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Intersecter deux surfaces quand ce sont deux primitives dont l'intersection a
une forme fermée : droites, cercles, ellipses. Ce sont les cas les plus
fréquents des pièces mécaniques (perçages, faces planes, révolutions
coaxiales). Les autres paires viennent avec le cheminement de la PR 4.

## 2. API

```cpp
template <T> struct IntersectionCurve {
    shared_ptr<Curve<T, 3>> curve; T a, b;        // primitive restreinte à [a, b]
    shared_ptr<Curve<T, 2>> pc1, pc2;             // pcurves, même paramètre
    T deviation; bool exact, closed, tangent;
};
template <T> struct SurfaceSurfaceHit    { T u1, v1, u2, v2; point<T, 3> p; };
template <T> struct SurfaceSurfaceResult { vector<IntersectionCurve> curves; vector<SurfaceSurfaceHit> points; bool coincident, complete; };

auto intersect(s1, box1, s2, box2, IntersectOptions = {}) -> SurfaceSurfaceResult;
auto intersect(s1, s2, IntersectOptions = {}) -> SurfaceSurfaceResult;   // bornes des surfaces
```

`complete = false` signale une paire sans forme fermée. C'est le cas de toute
paire qui n'est pas faite de deux primitives, jusqu'à la PR 4.

## 3. Choix d'architecture et justification

### 3.1 Deux étapes

1. **Courbes porteuses** (`analytic_surfaces.h`). La paire donne des droites
   non bornées, des cercles et des ellipses sur un tour complet, et des points
   de contact isolés. Elle peut aussi constater que les surfaces sont
   confondues. On ne regarde pas encore les boîtes.
2. **Restriction** (`surface_surface.h`). On garde de chaque porteuse les
   parties situées **dans les deux boîtes**, coupées là où un paramètre
   périodique traverse le **seam** d'une boîte. Chaque pcurve reste ainsi
   continue dans sa boîte, ce que la découpe des faces (PR 7) attend.

### 3.2 Paires traitées

| Paire | Méthode | Résultat |
|---|---|---|
| plan / plan | produit vectoriel des normales | une droite, ou parallèles, ou confondus |
| surfaces **coaxiales** : cylindre, cône, sphère, tore entre eux, sphère centrée sur l'axe, plan normal à l'axe, deux sphères, plan et sphère | intersection 2D des **méridiennes** dans le demi-plan (ρ ≥ 0, h) : droites (cylindre, cône à deux nappes, plan) et cercles (sphère, tore) | cercles autour de l'axe, points sur l'axe, cercles tangents, surfaces confondues |
| plan / cylindre, axe parallèle au plan | distance à l'axe | deux génératrices, une génératrice tangente, rien |
| plan / cylindre ou cône, oblique | **conique du plan** (forme quadratique en coordonnées du plan) : centre, axes propres | ellipse ; paraboles et hyperboles laissées à la PR 4 |
| plan passant par le sommet d'un cône | directions annulant la forme quadratique | deux génératrices, une tangente, ou le sommet |
| plan contenant l'axe d'un tore | — | les deux cercles méridiens (isos u exactes) |
| deux cylindres d'axes parallèles | cercles de la section droite | deux droites, une droite tangente, confondus |

La méthode des méridiennes réunit neuf paires en une seule intersection
droite/cercle en 2D. Un cercle porteur est paramétré par l'**angle u de la
surface qui porte l'axe** : sa pcurve sur cette surface est donc un segment.

### 3.3 Restriction aux boîtes

- **Droites** : d'abord coupées par l'intersection des boîtes 3D des deux
  surfaces.
- **Échantillonnage.** On parcourt la porteuse en 512 pas. À chaque pas, on
  calcule ses paramètres sur chaque surface par l'inverse en formule fermée
  (`parameters`), ramenés dans le tour de la boîte.
- **Bichotomie.** Les entrées, sorties et sauts de seam sont localisés par 80
  bichotomies. Un saut est un écart de plus d'un demi-tour sur un paramètre
  périodique.
- **Bord de boîte.** Une valeur posée sur un bord de la boîte (le seam même)
  appartient aux deux côtés. On compare donc avec la dernière valeur prise
  hors des bords : un cercle qui court le long du seam du tore (v = 0) n'est
  pas haché.
- **Nettoyage.**
  - les parties réduites à un point sont retirées ;
  - les extrémités à la tolérance du bout de la plage y sont recollées ;
  - les deux bouts d'un tour complet sont réunis s'ils ne sont pas séparés par
    un seam.

**Limite** : une partie de porteuse plus courte qu'un pas (1/512 de la
plage), entrée et sortie entre deux échantillons, peut être manquée. Pour une
droite, la plage est la région commune des deux surfaces. Pour un cercle,
c'est un tour.

### 3.4 Pcurves

- **Segment exact** quand les paramètres de la surface varient linéairement le
  long de la courbe (écart relatif < 1e-10 sur 9 points). C'est le cas :
  - d'un cercle sur une surface coaxiale (u = t + constante, v constant) ;
  - d'une génératrice ;
  - d'une droite dans un plan ;
  - d'un méridien de tore.
- **Interpolée** sinon : B-spline cubique sur le paramètre de la courbe, avec
  16 à 2 048 points, jusqu'à un écart inférieur à tol / 10, mesuré sur 1 000
  points. Exemples : ellipse dans un plan, cercle sur une sphère dont l'axe
  diffère.
- Les paramètres périodiques sont **déroulés depuis le milieu**, et les bouts
  posés sur un seam suivent leurs voisins.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| plan / tore : deux cercles si perpendiculaire | inclus dans les coaxiaux | même méthode |
| paraboles et hyperboles par cheminement (question 6) | `complete = false` jusqu'à la PR 4 | conforme |
| surfaces confondues : découpe en 2D (§ 6.3) | drapeau `coincident` seul | la région commune viendra avec la découpe |
| — | `complete`, `tangent` dans les résultats | signaler les paires non traitées et les contacts tangents |
| — | courbes coupées aux seams des deux boîtes | pcurves continues dans leur boîte |

## 5. Tests (`tests_intersect_surfaces`, 5 cas, 14 617 assertions, 60 ms)

Chaque résultat est vérifié dans les deux ordres (mêmes nombres). On mesure en
21 points la distance entre la courbe 3D et l'image de chaque pcurve par sa
surface (≤ écart annoncé ≤ tol), et la présence des pcurves dans leur boîte.
Tous les repères sont inclinés.

| Test | Couvre |
|---|---|
| `planes` | droite (longueur 10 dans les boîtes), parallèles, confondus |
| `plane_and_cylinder` | cercle fermé et u = t ; demi-cylindre (demi-cercle) ; boîte privée d'un quart autour de l'origine du paramètre (un arc de 3π, réuni à travers l'origine) ; deux génératrices, génératrice tangente ; ellipse oblique coupée au seam (périmètre de Ramanujan) |
| `plane_sphere_cone` | cercle de la sphère, point tangent, plan éloigné ; cône : cercle (exact sur le cône), ellipse, deux génératrices par le sommet (longueur), parabole → `complete = false` |
| `coaxial` | cylindre/sphère et cylindre/tore (repères tournés : arcs coupés aux seams, exacts), cylindre/cône, plan/tore (deux cercles, puis un cercle tangent sur le seam v = 0), sphère/cône, cylindres confondus |
| `spheres_tori_cylinders` | deux sphères (cercle, point tangent, confondues), plan contenant l'axe d'un tore (méridiens fermés), plan oblique sur tore → incomplet ; cylindres parallèles (deux droites bornées par les deux boîtes, une tangente), cylindres gauches → incomplet |

Les suites `tests_intersect_curves` et `tests_intersect_curve_surface` passent
aussi.

## 6. Points à valider

1. **Méridiennes** pour toutes les paires coaxiales.
2. **Conique du plan** pour plan/cylindre et plan/cône obliques ; paraboles et
   hyperboles renvoyées au cheminement.
3. **Restriction par échantillonnage et bichotomie**, avec sa limite (§ 3.3).
4. **Coupure aux seams des deux boîtes** : une courbe peut être livrée en
   plusieurs arcs, chacun avec des pcurves continues.
5. **Pcurves exactes seulement si linéaires**, interpolées sinon avec écart
   mesuré.
6. **`complete = false`** pour les paires non traitées, jusqu'à la PR 4.

## 7. Revue de code (a posteriori, lot 2 de la revue du noyau BREP)

| Commentaire | Suite donnée |
|---|---|
| plan / plan : droite calculée avec `1 − (n1·n2)²`, qui s'arrondit à 0 entre 1e-10 et 1e-8 rad (droite à l'infini, exception « Curve eval inf » qui termine le booléen ; `clip_line` gardait ±DBL_MAX face à un NaN), et avec des décalages absolus (loin de l'origine : droite hors des plans, ou perdue) | point de la droite par `o + ((n2·(o2 − o))/(d·d)) (d × n1)`, `o` pied du milieu de la région commune sur le plan 1 ; `clip_line` rejette les entrées non finies ; parallèles selon la règle du point 5. Test `planes_at_a_small_angle` |
| tore fuseau (R < r, `DEGENERATE_TOROIDAL_SURFACE`) : le cercle miroir de la méridienne en ρ = −R manquait ; `ConicalSurface::parameters` et `ToroidalSurface::parameters` ignoraient le signe du facteur radial (pcurves décalées de π au-delà du sommet, écart 6) | méridienne : le cercle du tube et, pour un fuseau, son miroir ; `parameters` ajoute π à u où le facteur radial est négatif (cône au-delà du sommet ; tore fuseau quand le cercle de tube centré en −R est le plus proche). Autres utilisateurs : `project_point_global` et `curve_surface.h` s'en servent comme départ, ils y gagnent. Tests `spindle_torus_and_double_cone`, `tests_primitives.inverse_beyond_apex_and_spindle` |
| plans par le sommet d'un cône pointu marqués incomplets : `part[0] < sp.tb` comptait un recouvrement de 1e-32 avec l'intervalle de l'autre nappe (82/385, 111/384, 61/388 à 30°, 45°, 60°) | recouvrement réel exigé (`min − max > eps_t`). Test `planes_through_cone_apex` ; il échoue sur master, mais une fois `ConicalSurface::parameters` corrigé, la pcurve au-delà du sommet est juste et la règle stricte n'est plus qu'une garde |
| tangences et sommet acceptés à tol par les formes fermées, puis la restriction exigeait tol / 10 de chaque pcurve, même exacte (écarts dans ]0,1 tol, tol] : `complete = false` dans cinq configurations) | critère sur l'**ajustement** : `PrimitiveFit::error`, distance entre S(pcurve(t)) et le point de la surface aux paramètres de C(t), ≤ tol / 10, et écart à la courbe ≤ tol + tol / 10. Plus général que `exact ? tol : tol / 10` (couvre aussi une pcurve interpolée d'une porteuse proche de la surface) ; porteuse tangente non recentrée à mi-chemin (inutile ici). Test `tangent_within_tolerance` |
| cercle passant près d'un pôle de sphère : u tourne de π sur une longueur de l'ordre de la distance au pôle, l'échantillonnage uniforme (≤ 2 048) ne suit pas (écart jusqu'à 1,6e-4, incomplet ; 10/200 plan/sphère, 8/200 sphère/sphère) | échantillons **adaptatifs** : 16 au départ, un intervalle est coupé en deux quand l'ajustement manque tol / 10 en son milieu ou en l'un des 1 000 points de contrôle (au plus 4 096 échantillons, 64 tours) ; un point singulier prend un u linéaire en t de ses voisins. Test `circle_near_pole` |
| « parallèle / coaxial » décidé à 1e-10 rad fixe : le même cylindre deux fois, axe incliné de 1e-9 à 1e-7, donnait des points singuliers ou rien, incomplet, 300 ms | décidé par l'écart induit sur la **région commune** des boîtes (`Region`) : axes à moins de tol l'un de l'autre sur la région (`axis_gap`), plan normal à l'axe incliné de a avec sin a × portée ≤ tol, plan parallèle à un axe, plans, cylindres parallèles, plan contenant l'axe d'un tore ; grandeurs prises au milieu de la région. 1e-9 et 1e-8 : confondus. À 1e-7 l'écart vaut tol au bout du cylindre et la marge des boîtes échantillonnées (`SurfacePieces`, 0,27 R) le porte à 1,3e-6 : cheminement, incomplet comme avant. Test `near_parallel` |
| seules les droites étaient coupées à la région commune : une très longue ellipse (plan à 1e-6 de l'axe d'un cylindre) avait ses pcurves ajustées sur le tour entier, tol / 10 manqué | `clip_ellipse` : chaque coordonnée vaut c + a cos t + b sin t, les bornes de la région (élargie de tol) la coupent en φ ∓ acos(…) ; arc à travers l'origine rendu d'un seul tenant ; ajustement sur les arcs gardés seulement. Test `long_ellipse` |
| repère de plan loin de la quadrique : annulation dans le terme constant de la conique (5,7e-6 hors du cône à 1e6) | origine du plan ramenée au pied de l'origine de la quadrique ; dans la méridienne, le plan est pris à sa traversée de l'axe. Test `far_plane_origin` |

Changements induits :

- les ajustements adaptatifs atteignent tol / 10 au plus juste, là où le doublement uniforme le dépassait souvent largement ; une traversée rasante d'une borne glissait alors le long de la courbe de l'erreur divisée par l'angle (`short_parts` : arc de 2,67e-4 au lieu de 2,83e-4). Chaque racine de la pcurve est désormais **polie** sur les paramètres de la courbe elle-même (bissection sur le changement de signe trouvé au-delà du pas de Newton, sinon racine gardée : contact) ;
- des boîtes 3D disjointes donnent un résultat vide sans calcul (ni forme fermée ni cheminement) ;
- sondes de la revue après correction : p1 plus d'exception, droite exacte ou plans confondus ; p2 écart 1e-15, complet ; p3 longueur 21,7656 et 6 courbes ; p4 cylindres confondus jusqu'à 1e-8 ; p5, p6, p7, p8 tous complets ; p10 0 incomplet ; p11 0 signalé ; p12 0 incomplet, les écarts de longueur restants (3 plan/cylindre, 2 plan/cône) viennent du cheminement et sont identiques sur master ;
- temps : `tests_brep_boolean` 25,0 s contre 25,0 s (meilleur de trois, machine chargée), `tests_brep_section` 0,13 s contre 0,12 s, les sept cas d'origine de `tests_intersect_surfaces` 0,09 s contre 0,10 s.
