# Écriture STEP — PR préliminaire P1 : primitives analytiques

| | |
|---|---|
| PR | branche `feat/geom-primitives` |
| Document de référence | [step_writer.md](step_writer.md), § 12.1 (décidé : primitives analytiques dans le BREP, à côté des NURBS) |
| Fichiers | `gbs/primitives.h` (nouveau), `gbs-brep/edit.h`, `gbs-brep/pcurve.h`, `gbs-io/step/geometry.h`, `gbs-io/step/reader.h`, `python/gbsBindSurfaces.cpp`, `tests/tests_primitives.cpp` (nouveau), `tests/tests_step_read.cpp`, `python/tests/test_primitives.py` (nouveau), `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Donner à gbs des **primitives analytiques**, utilisables dans le BREP à côté
des NURBS, et faire que le lecteur STEP les **garde** au lieu de les convertir.
L'écrivain pourra ainsi écrire un `cylindrical_surface` pour un cylindre lu ou
construit.

## 2. API (`gbs/primitives.h`)

```cpp
// surfaces : ElementarySurface<T> (Surface<T,3>) sur une boîte finie {u1, u2, v1, v2}
Plane<T>(ax2, box)                         S = O + u X + v Y
CylindricalSurface<T>(R, ax2, box)         S = O + R (cos u X + sin u Y) + v Z
ConicalSurface<T>(R, a, ax2, box)          S = O + (R + v tan a)(cos u X + sin u Y) + v Z
SphericalSurface<T>(R, ax2, box)           S = O + R cos v (cos u X + sin u Y) + R sin v Z
ToroidalSurface<T>(R, r, ax2, box)         S = O + (R + r cos v)(cos u X + sin u Y) + r sin v Z
    value(u, v, du, dv)  bounds()  parameters(p)  to_nurbs()  nurbs_parameters(u, v)  iso_u(u, v1, v2)  moved(f)

// courbes : ElementaryCurve<T> (Curve<T,3>) sur une plage finie [t1, t2]
LineCurve<T>(O, V, t1, t2)                 C = O + t V
EllipseCurve<T>(r1, r2, ax2, t1, t2)       C = O + r1 cos t X + r2 sin t Y   (cercle : r1 = r2)
    value(t, d)  bounds()  parameter(p)  to_nurbs()  nurbs_parameter(t)  reversed()  moved(f)
```

Python : mêmes classes dans `gbs` (`gbs.CylindricalSurface(2., ax, box)`…).

## 3. Choix d'architecture et justification

### 3.1 Paramétrages de STEP, domaine fini

- Les formules sont celles d'ISO 10303-42, dans un repère `ax2` = {origine,
  axe, direction de référence} rendu orthonormé (`Frame3`, Y = Z × X). Ce sont
  les mêmes que celles de la PR 4 de la lecture et de `bselementary.h`.
- Le **domaine est donné à la construction** : la plage d'une arête, la boîte
  d'une face. Le BREP exige des bornes finies, et une face sur un cylindre
  complet a pour boîte [u₀, u₀ + 2π] : sa fermeture se détecte comme pour une
  NURBS.
- **Les angles ne sont pas ramenés dans [0, 2π]** : un arc de 5,9 à 6,6 rad
  chevauche 0. La seule contrainte est une plage d'au plus un tour.
- **Dérivées exactes à tout ordre**, par cos(x + nπ/2) et sin(x + nπ/2).
- **Inverse en formule fermée** (`parameters`), qui rend l'angle dans le tour
  de la boîte (ou de la plage).
- **Conversion exacte** (`to_nurbs`) par `bselementary.h` (PR 3 de la
  lecture), avec la correspondance de paramètres exacte (`nurbs_parameters`).
- **Isoparamétrique analytique** (`iso_u`) : une droite pour le plan, le
  cylindre et le cône, un cercle méridien pour la sphère et le tore. Elle sert
  à l'arête de seam insérée par le lecteur.
- **Déplacement rigide** (`moved`) : le repère est déplacé, la primitive le
  reste.

### 3.2 De nouvelles courbes plutôt que `Line` et `Circle3d`

La conception prévoyait de réutiliser `gbs::Line` et `Circle3d`. Elles ne
conviennent pas aux arêtes :

- `Circle3d` a un domaine figé à [0, 2π] : un arc qui chevauche 0 serait
  refusé par `make_edge`.
- les bornes de `Line` vont de `numeric_limits::min()`, le plus petit réel
  **positif**, à `max()`. Une arête aux paramètres négatifs serait refusée.
  C'est un défaut de `Line`, non corrigé ici pour ne pas changer son
  comportement ailleurs.

`LineCurve` et `EllipseCurve` portent donc une plage explicite et le
paramétrage de STEP : la droite par un vecteur de norme quelconque, l'ellipse
par l'angle excentrique.

`reversed()` donne la même courbe parcourue à rebours, de paramètre t → −t :
le vecteur est opposé pour la droite, et Y, Z sont opposés pour l'ellipse, ce
qui garde un repère direct. Le lecteur s'en sert pour une arête dont
`same_sense` est faux.

### 3.3 Le BREP : rien à changer, sauf deux points

Les builders, la couture, `check`, les volumes et les boîtes passent par
l'interface générique : ils fonctionnent sur les primitives sans
modification. Le test construit des faces naturelles sur cylindre, sphère et
tore, puis un cylindre fermé à partir de boucles ordonnées, dont les cercles
sont des `EllipseCurve` et dont l'aire en (u, v) vaut exactement 2πh.

- **`copy_shape`** déplace aussi les primitives (`moved`), qui restent des
  primitives.
- **`project_point_global`** (première projection d'une pcurve) part, pour une
  primitive, de son **inverse en formule fermée**, puis affine par
  Gauss-Newton. La recherche générique s'égarait près du seam d'un cylindre
  analytique : c'est le même défaut qu'en PR 7 de la lecture, sous une autre
  forme. Le repli générique reste en place si la formule fermée ne suffit pas.

### 3.4 Le lecteur STEP garde les primitives

La question 3 de la lecture est révisée (§ 12.1 de la conception de
l'écriture) :

| Entité STEP | Avant | Maintenant |
|---|---|---|
| `line` | segment B-spline de degré 1 | `LineCurve` (même paramètre) |
| `circle`, `ellipse` | arc NURBS rationnel | `EllipseCurve`, plage de l'arête en angles du fichier |
| `plane`, `cylindrical_surface`, `conical_surface`, `spherical_surface`, `toroidal_surface` (et leurs `rectangular_trimmed_surface`, `offset_surface`) | NURBS rationnelle exacte | primitive sur la boîte de la face |
| B-spline, révolution, extrusion, composite | NURBS | inchangé |

- `CurveDef::to_curve` et `SurfaceDef::to_surface` choisissent entre primitive
  et NURBS ; `to_nurbs` reste disponible.
- Les **paramètres des arêtes** sont désormais ceux du fichier : un cercle
  complet va de t à t + 2π.
- L'**arête de seam insérée** (PR 7 de la lecture) est l'isoparamétrique
  analytique de la primitive.
- Les **pcurves restent recalculées**, comme avant. Elles sont maintenant
  dans le paramétrage de STEP : celle d'un cercle sur un cylindre est un
  segment. Utiliser directement les pcurves du fichier devient possible, mais
  n'est pas fait ici.
- **Export IGES** : une face analytique y est écrite par l'approximation
  existante, dans `approx_tol` (§ 12.1 de la conception de l'écriture).

## 4. Écarts par rapport au document de conception

| Document (§ 12.1) | Implémentation | Raison |
|---|---|---|
| `Line` et `Circle3d` réutilisées, `Ellipse3d` ajoutée | `LineCurve` et `EllipseCurve` nouvelles, à plage explicite | domaines de `Line` et `Circle3d` inadaptés aux arêtes (§ 3.2) |
| — | `iso_u`, `reversed`, `moved` | besoins du lecteur (seam, `same_sense`) et de `copy_shape` |
| — | projection globale par l'inverse en formule fermée | robustesse près du seam (§ 3.3) |

## 5. Tests

| Test | Couvre |
|---|---|
| `tests_primitives` (8 cas, 3 805 assertions) | chaque surface égale à sa NURBS exacte sur une boîte qui chevauche 0 ; formules dans un repère incliné ; dérivées jusqu'à l'ordre 4 contre différences finies ; inverse et isoparamétriques ; ellipse à cheval sur 0, droite de vecteur non unitaire, retournement ; déplacements rigides ; arguments invalides ; faces naturelles, cylindre fermé à partir de boucles ordonnées (aire 2πh, volume πR²h), copie déplacée qui reste en primitives |
| `tests_step_read.elementary_geometry_stays_analytic` | `cylinder.stp` lu en un cylindre, deux plans, deux cercles et une droite, plages d'arête de 2π ; une B-spline reste B-spline |
| `python/tests/test_primitives.py` | liaisons Python |

Toutes les suites BREP et STEP passent sans changement de leurs attentes :
volumes, comptes, validité. Le binaire de tests OCCT passe aussi en entier
(25 cas, dont les 20 formes et l'assemblage écrits par OCCT, désormais lus en
primitives).

## 6. Points à valider

1. **Primitives analytiques sur domaine fini**, paramétrages de STEP, angles
   non ramenés dans [0, 2π].
2. **`LineCurve` et `EllipseCurve`** plutôt que `Line` et `Circle3d`.
3. **Lecteur STEP** : géométrie élémentaire gardée analytique, paramètres
   d'arête du fichier.
4. **Projection globale** partant de l'inverse en formule fermée.
5. **Export IGES** des faces analytiques par approximation, en attendant une
   conversion exacte.

## 7. Revue de code (a posteriori, lot 7 de la revue du noyau BREP)

Pas de défaut trouvé dans les primitives : paramétrages, dérivées de tout
ordre, paramètres inverses, sens inversés et passage aux NURBS exactes. Le lot 2
avait déjà corrigé les paramètres des cônes et des tores
([brep2_pr03_analytic.md](brep2_pr03_analytic.md) § 7).

| Commentaire | Suite donnée |
|---|---|
| `copy_shape` ne savait déplacer que les primitives (`moved`), les B-splines et les courbes sur de telles surfaces, et levait une exception sinon. Les balayages d'un profil bordé par une autre courbe et les booléens d'un opérande qui en porte une échouaient ([brep2_pr08_sweep_offset.md](brep2_pr08_sweep_offset.md) § 7) | corrigé : toute autre courbe ou surface est déplacée exactement par `MovedCurve` / `MovedSurface` (`edit.h`) : point déplacé, dérivées tournées. Deux déplacements se composent, et le placement identité partage la géométrie telle quelle. `Placement::rotate` est ajouté. Testé (`tests_brep_edit.copy_shape_of_any_curve`, et `tests_brep_boolean.failure_leaves_the_model`, dont la fusion passe désormais) |
