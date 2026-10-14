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
