# Palier 2 — PR 4 : surface-surface par cheminement

| | |
|---|---|
| PR | branche `feat/intersect-marching` |
| Document de référence | [brep_palier2.md](brep_palier2.md), § 2.6 et 2.7 ; questions 4 et 7 |
| Fichiers | `gbs/intersect/marching.h`, `gbs/intersect/intersection_curve.h` (nouveaux, types sortis de `surface_surface.h`) ; `gbs/intersect/surface_surface.h` (aiguillage), `gbs/intersect/curves.h` (options) ; `tests/tests_intersect_marching.cpp` (nouveau), `tests/tests_intersect_surfaces.cpp` ; `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Intersecter **toute** paire de surfaces : B-splines libres, révolutions,
offsets, et les paires de primitives sans forme fermée (paraboles, tores
coupés obliquement, cylindres gauches). Avec la PR 3, `intersect(s1, box1, s2,
box2)` couvre désormais tous les cas.

## 2. API

Pas de nouvelle fonction publique : `intersect` passe au cheminement quand la
paire n'a pas de forme fermée.

```cpp
struct IntersectOptions { …; bool closed_forms = true; T max_turn = 0.05; };   // forcer le cheminement ; rotation de la tangente par pas
struct SurfaceSurfaceResult { curves; points; coincident; bool singular; bool complete; };
```

- `singular` : une branche s'est arrêtée sur un point singulier (tangence,
  point de branchement), ou les surfaces se touchent.
- `complete = false` : le plafond de paires de carreaux a été atteint.

## 3. Choix d'architecture et justification

### 3.1 Points de départ par subdivision

- **Découpe.** Les deux surfaces sont coupées en carreaux (`SurfacePieces`,
  PR 2). Les paires dont les boîtes se touchent sont coupées en deux jusqu'à
  ce que les deux carreaux soient plats (écart au plan ≤ 2 % de la taille).
- **Départ.** Au centre de chaque paire plate, un Gauss-Newton à **norme
  minimale** sur (u1, v1, u2, v2) pose le point sur l'intersection. La
  résolution passe par une **SVD tronquée** (Eigen) : là où les surfaces sont
  tangentes, J est de rang 2 et le pas va vers les points les plus proches. On
  trouve ainsi les contacts ponctuels, même au pôle d'une sphère.
- **Carreaux presque parallèles** (cônes de normales qui se recouvrent). Ils
  peuvent cacher une petite boucle fermée ou une tangence. S'il n'y a pas de
  départ en leur centre, on les découpe encore, jusqu'à 100 × 10 tol.
- **Pourquoi pas les bords seulement.** Ce critère remplace les départs pris
  sur les bords des domaines (§ 2.6 du document de conception), qui ne trouvent
  pas les boucles intérieures. Les bords sont couverts de toute façon :
  chaque branche est poursuivie jusqu'à eux.

### 3.2 Cheminement

- **Prédiction.** Un pas `h` le long de t = n1 × n2, reporté sur chaque
  surface par moindres carrés dans son plan tangent.
- **Correction.** Newton 4 × 4 : les trois équations S1 = S2, plus le plan
  normal au pas.
- **Pas.** Il est accepté si la tangente tourne de moins de `max_turn`
  (0,05 rad). Il est ensuite ajusté à la rotation mesurée (de × 0,5 à × 2, au
  plus 1/10 de la taille des surfaces), et coupé en deux en cas d'échec.
- **Bord d'une boîte.** Le paramètre qui sort le premier est fixé à sa borne,
  et le point est corrigé par Newton 3 × 3. Le dernier pas obéit à la même
  limite de rotation que les autres.
- **Hors de la surface.** Un point en dehors des bornes propres d'une surface
  n'est jamais évalué, car une B-spline lève une exception : la correction
  rend la main et le bord est traité comme ci-dessus.
- **Boucle fermée.** On revient au départ en 3D **et** en paramètres. Une
  courbe qui rejoint son départ à travers un seam n'est pas une boucle : elle
  s'arrête au seam.
- **Point singulier.** Le sinus de l'angle des normales passe sous
  `angular_tol`, ou le pas tombe sous 10 tol. Le sinus est aussi contrôlé **au
  milieu de chaque pas** : un creux au milieu signale un point de branchement
  enjambé, et le pas est raccourci jusqu'à l'atteindre.
- **Dédoublonnage.** Un départ n'est pas suivi s'il est à 10 tol d'une branche
  déjà tracée. La distance est mesurée à la polyligne, avec une marge pour la
  flèche entre ses points (1 % de chaque segment, majoration de la flèche à
  0,05 rad par pas). La boîte de chaque branche est élargie de cette flèche.

### 3.3 Courbes

- **Interpolation.** La B-spline cubique 3D et les deux pcurves sont
  interpolées sur le **même paramètre** (longueur de corde).
- **Densification.** On mesure l'écart aux quarts de chaque intervalle. Un
  point corrigé est inséré là où l'écart dépasse tol / 2, au plus 8 fois.
- **`closed`** veut dire « fermée en 3D ». C'est le même sens que pour les
  formes fermées : un cercle coupé au seam est fermé en 3D, et ses pcurves vont
  d'un bord à l'autre.
- **Seams.** Les paramètres périodiques ne sont pas déroulés. Une courbe
  s'arrête au seam d'une boîte, comme les courbes des formes fermées.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| départs sur les bords, puis boucles intérieures par subdivision et cônes de normales | départs par subdivision seule, carreaux parallèles redécoupés | couvre les deux cas avec un seul mécanisme |
| espace déroulé pour les surfaces périodiques | arrêt au seam | même convention que la PR 3, pcurves dans leur boîte |
| point singulier : `NonGeneric` aux booléens | drapeau `singular` | l'erreur sera posée par les booléens (PR 9) |
| — | `closed_forms`, `max_turn` dans les options | tests et réglage de la précision |
| — | SVD tronquée pour les départs | contacts tangents, pôles |

## 5. Tests

`tests_intersect_marching` : 6 cas, 12 347 assertions, 0,2 s. Chaque courbe
est vérifiée en 51 points sur les deux surfaces (≤ 2 tol), et ses pcurves dans
leurs boîtes.

| Test | Couvre |
|---|---|
| `against_closed_forms` | cinq paires de primitives (plan/cylindre droit et oblique, génératrices, sphère/plan oblique, cylindre/sphère) par **cheminement forcé** : mêmes longueurs que les formes fermées à 1e-5, aucun point singulier |
| `crossed_cylinders` | rayons 1 et 0,5 à axes perpendiculaires : longueur comparée à l'intégrale de référence (2 × 200 000 pas) ; rayons égaux (Steinmetz) : `singular` |
| `plane_cuts` | parabole sur un cône ; plan oblique sur un tore, primitive et NURBS exacte (mêmes longueurs) ; NURBS d'un cylindre et plan (longueur 4π) |
| `free_surfaces` | bosse B-spline et plan : **boucle fermée intérieure** ; bosse et vague B-spline ; sphère tangente à un plan, contact hors pôle et **au pôle** |
| `singular_points_and_limits` | Viviani (point double) : `singular` ; plafond de paires : `complete = false` |

`tests_intersect_surfaces` passe aussi : ses trois paires qui étaient
« incomplètes » (parabole, tore oblique, cylindres gauches) sont désormais
vérifiées par cheminement.

## 6. Points à valider

1. **Départs par subdivision seule**, carreaux parallèles redécoupés s'ils
   n'en donnent pas.
2. **Arrêt au seam** plutôt que l'espace déroulé.
3. **Pas réglé par la rotation de la tangente** (0,05 rad), et densification
   des courbes jusqu'à tol / 2.
4. **Détection des points singuliers** par le sinus des normales, aux bouts et
   au milieu des pas ; drapeau `singular` plutôt qu'une erreur.
5. **SVD tronquée** pour les départs : contacts tangents et pôles.
6. **`closed` = fermée en 3D**, pour les deux méthodes.
