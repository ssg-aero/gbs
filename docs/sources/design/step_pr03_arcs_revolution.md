# Lecture STEP — PR 3 : arcs, révolution et extrusion NURBS exacts

| | |
|---|---|
| PR | branche `feat/geom-arcs-revolution` |
| Document de référence | [step_reader.md](step_reader.md), § 3.2 et 3.3 (« Constructions à ajouter à gbs »), question 3 (décidée : surfaces analytiques converties en NURBS rationnelles exactes) |
| Fichiers | `gbs/bselementary.h` (nouveau), `tests/tests_bselementary.cpp` (nouveau), `python/gbsBindBuildSurfaces.cpp`, `python/tests/test_elementary.py` (nouveau), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Donner à gbs les constructions **exactes** dont la PR 4 a besoin pour traduire
la géométrie STEP sans approximation : arc de cercle et d'ellipse entre deux
angles, révolution rationnelle d'une NURBS autour d'un axe, extrusion d'une
NURBS, et les quatre surfaces élémentaires qui s'en déduisent (cylindre, cône,
sphère, tore). Jusqu'ici, gbs ne savait construire que le cercle et l'ellipse
complets, centrés dans le plan xy et paramétrés sur [0, 1].

## 2. API

```cpp
#include <gbs/bselementary.h>

// arcs : C(θ) = O + r1 cos θ X + r2 sin θ Y, θ de theta1 à theta2 (radians)
auto c = build_circle_arc<T>(r, theta1, theta2, ax);          // ax = {centre, axe, direction de départ}
auto e = build_ellipse_arc<T>(r1, r2, theta1, theta2, ax);    // r1 selon la direction de référence
auto c2 = build_circle_arc<T>(r, theta1, theta2, point<T,2>{x, y}); // arc 2D
auto g = build_ellipse_arc<T, dim>(r1, r2, theta1, theta2, O, X, Y); // forme générale, toute dimension

// passage exact angle <-> paramètre
T u = ellipse_arc_parameter(theta1, theta2, theta);
T a = ellipse_arc_angle(theta1, theta2, u);

// surfaces
auto s = build_revolution(generatrix, ax1{origine, direction}, theta1, theta2); // BSSurfaceRational
auto x = build_extrusion(curve, V, v1, v2);       // BSSurface ou BSSurfaceRational selon la courbe
auto cyl = build_cylinder<T>(R, ax, v1, v2, theta1, theta2);
auto co  = build_cone<T>(R, semi_angle, ax, v1, v2, theta1, theta2);
auto sph = build_sphere<T>(R, ax, theta1, theta2, v1, v2);
auto tor = build_torus<T>(R, r, ax, theta1, theta2, v1, v2);
```

Python : mêmes noms dans `gbs` (`build_circle_arc`, `build_ellipse_arc`,
`ellipse_arc_parameter`, `ellipse_arc_angle`, `build_revolution`,
`build_extrusion`, `build_cylinder`, `build_cone`, `build_sphere`,
`build_torus`).

## 3. Choix d'architecture et justification

### 3.1 Un en-tête léger dédié, `gbs/bselementary.h`

Les constructions de courbes vivent d'habitude dans `bscbuild.h` et celles de
surfaces dans `bssbuild.h`. Ce dernier tire le loft, le Gordon et Boost : le
lecteur STEP n'en a pas besoin. Le nouvel en-tête ne dépend que de `bscurve.h`
et `bssurf.h`. `build_circle` et `build_ellipse` restent inchangés.

### 3.2 Découpage en quarts et paramètre égal à l'angle aux nœuds

Un arc d'angle Δ est découpé en n = ⌈Δ / 90°⌉ segments rationnels quadratiques
**d'angle égal** (*The NURBS Book*, § 7.5) : poids 1, cos(Δ/2n), 1, nœuds
doubles entre segments, 2n + 1 pôles.

Les nœuds sont les **angles des bornes de segments, en radians** : le domaine
de l'arc est [θ1, θ2] et non [0, 1]. À chaque nœud, paramètre et angle
coïncident. Entre deux nœuds, une NURBS rationnelle ne peut pas suivre l'angle
linéairement, mais la relation est connue en forme close : pour un segment
d'angle d, de milieu m et de bornes a, a + d,

```
tan((θ − m) / 2) = s · tan(d / 4),   u = a + (s + 1) d / 2,   s ∈ [−1, 1]
```

`ellipse_arc_parameter` et `ellipse_arc_angle` l'appliquent : la conversion est
exacte dans les deux sens. La PR 4 en aura besoin pour un `trimmed_curve` coupé
par angles sur un cercle ou une ellipse. Le découpage dépend seulement de
(θ1, θ2) : les deux fonctions retrouvent les mêmes segments que la construction.

Les segments d'angle égal rendent la courbe **C¹** aux nœuds (le test le
vérifie), pas seulement G¹. Pour l'ellipse, θ est l'**angle excentrique**,
comme le paramètre de l'`ellipse` STEP. L'ellipse est l'image affine du cercle,
donc exacte avec les mêmes poids.

### 3.3 Repères à la manière de STEP

`ax2 = {origine, axe, direction de référence}` est déjà la convention de
`Circle` dans gbs, et c'est celle d'`axis2_placement_3d`. La direction de
référence est rendue orthogonale à l'axe (projection, comme le prescrit
ISO 10303-42), Y = Z × X. Un axe nul ou une référence parallèle à l'axe lève
`std::invalid_argument`.

### 3.4 Révolution : u = rotation, v = génératrice

`build_revolution` (*The NURBS Book*, A8.1) fait tourner chaque pôle de la
génératrice autour de l'axe : pôle (i, j) = O_j + c_i X_j + s_i Y_j, avec
X_j = P_j − O_j (vecteur radial, non normé) et Y_j = Z × X_j. Le poids est le
produit du poids de l'arc et de celui de la génératrice.

- **u est la rotation, v le paramètre de la génératrice**, inchangé : c'est la
  convention de `surface_of_revolution` et des surfaces élémentaires STEP, et
  celle du `nurbs_cylinder` des tests du palier 1. Les pôles sont rangés u le
  plus rapide, comme partout dans gbs.
- Un pôle de la génératrice **sur l'axe** donne X_j = 0 : toute la rangée est
  confondue. C'est le côté dégénéré attendu au pôle d'une sphère ou au sommet
  d'un cône, que `surface_closure` et `make_face` du palier 1 savent déjà
  traiter.
- La génératrice peut être polynomiale (`BSCurve`) ou rationnelle ; la surface
  est toujours rationnelle.

### 3.5 Extrusion

S(u, v) = C(u) + v V, de degré 1 en v, v ∈ [v1, v2], avec V non normé comme
le vecteur de `surface_of_linear_extrusion`. En coordonnées homogènes, chaque
pôle est translaté de w · v · V. La surface est rationnelle si et seulement si
la courbe l'est.

### 3.6 Surfaces élémentaires dans gbs, pas dans le lecteur

Le plan les rangeait dans la PR 4, côté lecteur. Ce sont quatre lignes chacune
au-dessus des constructions précédentes, et elles serviront aussi au palier 2
(primitives BREP). Elles sont donc ici, dans gbs, avec les paramétrisations de
ISO 10303-42 :

| Surface | Construction | S(u, v) |
|---|---|---|
| cylindre | extrusion de l'arc de cercle selon Z | O + R (cos u X + sin u Y) + v Z |
| cône | révolution du segment génératrice | O + (R + v tan α)(cos u X + sin u Y) + v Z |
| sphère | révolution du méridien (arc de −π/2 à π/2) | O + R cos v (cos u X + sin u Y) + R sin v Z |
| tore | révolution du cercle méridien de centre O + R X | O + (R + r cos v)(cos u X + sin u Y) + r sin v Z |

Pour le cylindre et le cône, v est exactement le v de STEP (génératrice de
degré 1). Pour u, et pour le v de la sphère et du tore, le paramètre est
l'angle aux nœuds et se convertit avec `ellipse_arc_parameter`. Comme les
pcurves sont recalculées (question 3), le lecteur n'a besoin de cette
conversion que pour les bornes d'un `rectangular_trimmed_surface`.

### 3.7 Erreurs

Ces fonctions relèvent de la couche géométrique : elles lèvent
`std::invalid_argument` (θ2 ≤ θ1, plus d'un tour, rayon nul ou négatif, axe nul,
demi-angle de cône hors de ]−π/2, π/2[, latitude hors de [−π/2, π/2]), pas
`BRepError`. Le lecteur STEP les attrapera par entité, comme les
`P21AccessError`.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| Surfaces élémentaires construites par le lecteur (PR 4) | dans `gbs/bselementary.h` | réutilisables au palier 2 ; testées ici avec la géométrie |
| — | `ellipse_arc_parameter` / `ellipse_arc_angle` | `trimmed_curve` et `rectangular_trimmed_surface` sont coupés par angles |
| Liaisons Python en PR 10 | ajoutées ici pour la géométrie | petites, et utiles hors STEP |

## 5. Tests

`tests_bselementary` (6 cas, 12 133 assertions), repère incliné non orthonormé
en entrée :

| Test | Couvre |
|---|---|
| `circle_arcs` | arcs de 30° à 360°, angles négatifs, nombre de pôles, bornes = angles, extrémités exactes, point à l'angle via `ellipse_arc_parameter`, aller-retour avec `ellipse_arc_angle`, rayon et planéité en tout paramètre, paramètre = angle et tangente continue aux nœuds, arc 2D |
| `ellipse_arcs` | point à l'angle excentrique, équation implicite, planéité |
| `invalid_arguments` | les dix cas d'erreur |
| `revolution` | génératrice rationnelle et polynomiale, tour complet et partiel : S(u, v) = rotation de G(v) ; rangée dégénérée pour un pôle sur l'axe |
| `extrusion` | courbe rationnelle et polynomiale, type de retour |
| `elementary_surfaces` | cylindre, cône jusqu'à son sommet, sphère et ses pôles, tore partiel : formule STEP au point (angle, v) et équation implicite en tout paramètre (écarts ≤ 1e-12) |

`python/tests/test_elementary.py` : deux tests sur les liaisons.

## 6. Points à valider

1. **En-tête dédié** `gbs/bselementary.h`, sans dépendre de `bssbuild.h`.
2. **Domaine de l'arc en radians [θ1, θ2]**, paramètre égal à l'angle aux
   nœuds, conversion exacte par `ellipse_arc_parameter` / `ellipse_arc_angle`.
3. **Révolution : u = rotation, v = génératrice**, surface toujours rationnelle,
   rangée dégénérée sur l'axe.
4. **Surfaces élémentaires dans gbs** avec les paramétrisations STEP, avancées
   de la PR 4.
5. **`std::invalid_argument`** pour les erreurs de construction géométrique.
