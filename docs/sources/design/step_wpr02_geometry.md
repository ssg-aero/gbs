# Écriture STEP — PR 2 : géométrie

| | |
|---|---|
| PR | branche `feat/step-write-geometry` |
| Document de référence | [step_writer.md](step_writer.md), § 3 et 6 ; § 12.1 (primitives) |
| Fichiers | `gbs-io/step/geometry_writer.h` (nouveau), `gbs-io/to_nurbs.h` (nouveau, sorti de `iges_brep.h`), `gbs-io/iges_brep.h`, `gbs/primitives.h` (`iso_v`), `gbs/curveonsurface.h` (accesseurs `const`), `tests/tests_step_write_geometry.cpp` (nouveau), `tests/tests_primitives.cpp`, `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Écrire la géométrie de gbs en entités STEP, au-dessus de l'écrivain Part 21
(PR 1). Ce sont les courbes d'arêtes, les pcurves, les surfaces de faces et le
contexte de représentation (unités, incertitude). La topologie viendra en PR 3.

## 2. API

```cpp
GeometryWriter g{p21_writer, approx_tol};
g.context(unit_mm, uncertainty) -> #id     // contexte 3D : unité du modèle, radian, stéradian, incertitude
g.context2d() -> #id                       // espace paramétrique des pcurves
g.cartesian_point(p), g.direction(d), g.axis2(frame)
g.curve(courbe_3d, a, b) -> #id            // courbe d'arête sur [a, b]
g.curve2d(pcurve, a, b) -> #id
g.surface(surface) -> #id
g.stats()                                  // approximations : nombres, plus grand écart
```

## 3. Choix d'architecture et justification

### 3.1 Que devient chaque géométrie

| gbs | STEP | Paramétrage |
|---|---|---|
| `LineCurve` | `line` (point, `vector` de norme quelconque) | identique |
| `EllipseCurve` | `circle` si r1 = r2, sinon `ellipse` | identique (angle) |
| `Plane`, `CylindricalSurface`, `ConicalSurface`, `SphericalSurface`, `ToroidalSurface` | `plane`, `cylindrical_surface`, `conical_surface`, `spherical_surface`, `toroidal_surface` (`degenerate_toroidal_surface` si R ≤ r) | identique |
| `BSCurve`, `BSCurveRational` (3D et 2D) | `b_spline_curve_with_knots` ; instance complexe avec `rational_b_spline_curve` si rationnelle ; pôles **cartésiens** et poids | identique, **restreinte à [a, b]** |
| `BSSurface`, `BSSurfaceRational` | `b_spline_surface_with_knots` ; complexe si rationnelle ; grille transposée de l'ordre gbs (u le plus rapide) vers `[u][v]` | identique |
| pcurve B-spline de degré 1 à deux pôles | `line` 2D | identique |
| `CurveOnSurface` dont la pcurve est un segment isoparamétrique paramétré comme la courbe (côtés des faces naturelles) | l'**iso exacte** : droite ou cercle sur une primitive (`iso_u`, `iso_v`), B-spline sur une B-spline | identique |
| toute autre courbe ou surface | B-spline cubique interpolée à paramètres conservés, dans `approx_tol` | conservé (approché) |

Le paramétrage est toujours celui de gbs : les pcurves du modèle, exprimées
dans le paramétrage de la surface de leur face, restent valables dans le
fichier.

### 3.2 Approximation commune avec l'IGES (`gbs-io/to_nurbs.h`)

`iges_detail::to_nurbs` (palier 1, PR 9) passe dans `gbs-io/to_nurbs.h`
(namespace `gbs::nurbs_export`), sans dépendance à libIGES, avec ses compteurs
`ApproxStats`. Le rapport IGES en hérite : ses champs gardent leur nom, et
l'export IGES utilise l'en-tête commun. C'est la recommandation de la question
11. Toutes les suites IGES, ainsi que l'aller-retour IGES par OCCT, passent
sans changement.

### 3.3 Unités et contexte

Le contexte décrit l'**unité du modèle** (`Model::unitScale`, en mm par unité) ;
les coordonnées sont écrites telles quelles (question 10).

| mm par unité | Unité STEP |
|---|---|
| 1, 10, 1000, 0,001 | `si_unit` millimètre, centimètre, mètre, micromètre |
| 25,4 et 304,8 | `conversion_based_unit` `'INCH'` et `'FOOT'`, définies en millimètres |
| autre | `conversion_based_unit` `'MODEL UNIT'` définie en millimètres |

S'y ajoutent le radian et le stéradian, et l'incertitude
(`uncertainty_measure_with_unit`, `'distance_accuracy_value'`). Le test relit
chaque contexte par `read_units` de la lecture.

### 3.4 Compléments aux primitives et au cœur

- `ElementarySurface::iso_v(v, u1, u2)` : le parallèle d'une primitive, pour
  écrire exactement les côtés v = constante. C'est une droite pour le plan et
  un cercle pour les autres. Il vaut `nullptr` au pôle d'une sphère, où ce
  côté est une arête dégénérée, sans courbe.
- `CurveOnSurface::p_basisCurve()` et `p_basisSurface()` deviennent `const`
  (lecture seule). Aucun appelant n'est affecté.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| contexte : millimètre, mètre, pouce, unité par conversion | plus centimètre, micromètre, pied | unités fréquentes, sans coût |
| — | `iso_v` ajouté aux primitives | côtés v = constante des faces naturelles écrits exactement |
| — | accesseurs `const` de `CurveOnSurface` | nécessaires pour lire une courbe constante |

## 5. Tests (`tests_step_write_geometry`, 4 cas)

Chaque entité est écrite, relue par l'analyseur et le `GeometryReader` de la
lecture, puis comparée à l'original **au même paramètre**.

| Test | Couvre |
|---|---|
| `curves` | droite (vecteur non unitaire), cercle à cheval sur 0, ellipse, B-spline restreinte à une sous-plage, arcs rationnels (instance complexe), courbe quelconque approchée et comptée |
| `iso_lines_and_pcurves` | génératrice et parallèle d'un cylindre écrits en `line` et `circle`, iso d'une surface B-spline, sans approximation ; pcurve segment écrite en `line` 2D, pcurve B-spline |
| `surfaces` | plan, cylindre, cône, sphère, tore, tore dégénéré, B-spline rationnelle et non rationnelle (grille transposée) ; surface de révolution approchée et comptée |
| `units` | millimètre, centimètre, mètre, micromètre, pouce, pied, unité quelconque, relus par `read_units` ; incertitude ; erreurs |

`tests_primitives` vérifie en plus `iso_v`. Toutes les suites BREP, STEP et
IGES, les tests Python et tout le binaire OCCT passent.

## 6. Points à valider

1. **Primitives en entités analytiques, B-splines exactes**, paramétrage
   toujours identique à celui de gbs.
2. **Courbes d'arêtes restreintes** à la plage de l'arête.
3. **Côtés des faces naturelles écrits comme iso exactes.**
4. **Approximation commune** IGES et STEP (`gbs-io/to_nurbs.h`).
5. **Unités** : celle du modèle, coordonnées non converties.
