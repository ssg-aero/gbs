# BREP palier 1 — PR 5 : face depuis un wire posé sur une surface (`gbs-brep/pcurve.h`, `make_face(srf, wire)`)

| | |
|---|---|
| PR | branche `feat/brep-face-wire` |
| Document de référence | [brep_core.md](brep_core.md), section 6.2 ; notes précédentes [PR 1](brep_pr01_model.md) à [PR 4](brep_pr04_face_natural.md) |
| Fichiers | `gbs-brep/pcurve.h` (246 l., nouveau), `gbs-brep/builders.h` (+ trois surcharges de `make_face`, sept codes d'erreur), `gbs-brep/closure.h` (`uv_polygon`, `signed_area`, `inside`), `gbs-brep/brep`, `tests/tests_brep_face_wire.cpp` (308 l., nouveau), `tests/tests_brep_geom.h` (géométrie de test partagée, extraite de `tests_brep_face.cpp`), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Construire une face **découpée** : une surface bornée par un contour extérieur
et des trous, donnés soit comme des wires libres de courbes 3D (le cas du
disque sur un plan, d'un contour importé), soit comme des boucles de courbes
2D dans l'espace paramétrique (le cas du palier 2 : courbes de découpe issues
des intersections). Le cœur du travail est le calcul des **pcurves** : exactes
quand c'est possible, projetées et interpolées sinon, continues à travers le
seam d'une surface fermée. C'est l'équivalent de `BRepBuilderAPI_MakeFace(surface, wire)`
suivi de `ShapeFix_Edge::FixAddPCurve` côté OCCT.

## 2. Vue d'ensemble

```
make_face(m, srf, outer [, holes], opts)                 wires libres de courbes 3D
make_face(m, srf, outer_pcurves2d [, holes2d], opts)     boucles de courbes 2D (chemin du palier 2)

MakeFaceOptions<T> { tol = 1e-6, pcurve_tol = 1e-5, pcurve_degree = 3, n_samples_min = 9, n_samples_max = 1025 }

pour chaque co-arête :
  extract_pcurve(edge, srf)          CurveOnSurface sur CE shared_ptr de surface → pcurve exacte, déviation 0
  sinon project_pcurve(edge, srf…)   échantillons t_i → C(t_i) → projection (u_i, v_i)
                                     → continuité (seam, pôles) → interpolation B-spline aux paramètres t_i
                                     → déviation 3D aux milieux ≤ pcurve_tol, sinon n ← 2n − 1
orientation   : extérieur sens direct, trous sens indirect (ordre et sens des co-arêtes inversés si besoin)
trous         : tous leurs points dans le polygone extérieur (test pair-impair)
écriture      : pcurves posées sur les co-arêtes des wires, tolérances d'arêtes ≥ max(déviation, opts.tol), sommets ajustés
```

## 3. Choix d'architecture et justification

### 3.1 Extraction exacte d'abord

Si la courbe 3D d'une arête est une `CurveOnSurface` **sur le même objet
surface** (même `shared_ptr`), sa courbe 2D est la pcurve, sans calcul ni
erreur. C'est le cas des faces à bornes naturelles (PR 4), du chemin « boucles
2D » de cette PR, et ce sera celui des courbes d'intersection du palier 2.
L'identité est testée par pointeur et non par géométrie : deux surfaces égales
mais distinctes passent par la projection (le test
`projection_on_curved_surface` s'en sert justement).

Limite assumée : une `CurveOnSurface` enveloppée (`CurveTrimmed`,
`CurveReversed`) n'est pas reconnue et passe par la projection, ce qui reste
correct mais approché.

### 3.2 Projection : Gauss-Newton amorcé, recherche globale en secours

- Le premier échantillon est projeté globalement : `extrema_surf_pnt` de gbs
  (grille 30 × 30 puis COBYLA borné), raffiné par Gauss-Newton.
- Les suivants partent du `(u, v)` précédent avec **Gauss-Newton** sur
  `‖S(u,v) − P‖²` (système 2 × 2 de la métrique `Su·Su, Su·Sv, Sv·Sv`),
  borné au rectangle. Convergence quadratique pour un point sur la surface,
  quelques évaluations par point au lieu de quelques centaines pour COBYLA.
- Si l'amorce mène à une mauvaise branche (distance > `pcurve_tol`, typiquement
  après le passage du seam où le point reste bloqué au bord), on refait une
  projection globale.
- Un point dont la distance reste > `pcurve_tol` signifie que l'arête n'est
  pas sur la surface : `EdgeOffSurface`.

### 3.3 Continuité sur les surfaces fermées et aux pôles

Pour chaque direction fermée (`surface_closure`, PR 4) :

- un échantillon situé sur le seam (`u = u1` ou `u = u2`) est **ambigu** ; il
  prend la valeur continue avec ses voisins résolus ;
- si toute l'arête est sur le seam, la valeur la plus proche du point où
  finit la co-arête précédente du wire (`hint`) est retenue. La première
  co-arête, qui n'a pas de précédente, est recalculée en fin de wire avec la
  fin de la dernière co-arête comme repère ;
- un saut de plus d'une demi-période entre deux échantillons consécutifs, ou
  entre la fin d'une co-arête et le début de la suivante, signifie que l'arête
  ou le wire **traverse le seam** : `CrossesSeam`. Il faut alors découper
  l'arête au seam, ce qui relève du palier 2.

Aux points singuliers (pôle de sphère, apex de cône), la coordonnée libre
(celle dont la dérivée s'annule) est copiée du voisin.

On ne « déroule » pas le paramètre au-delà des bornes, contrairement à ce
qu'envisageait le document (§ 6.2) : une B-spline ou une surface de
révolution évaluée hors de son domaine lève une exception dans gbs. Les
pcurves restent donc dans le rectangle paramétrique.

### 3.4 Interpolation aux paramètres de l'arête

Les points `(t_i, uv_i)` sont interpolés par `interpolate(Q, t, p)` de gbs avec
les paramètres **imposés** `t_i` : la pcurve est paramétrée comme l'arête
(convention de la PR 1), et l'égalité `S(pcurve(t_i)) = C(t_i)` est exacte aux
nœuds. La déviation est mesurée aux milieux des intervalles ; tant qu'elle
dépasse `pcurve_tol`, l'échantillonnage est doublé (`n ← 2n − 1`), jusqu'à
`n_samples_max` (sinon `PCurveApproximation`). La déviation finale devient la
tolérance de l'arête (`same_parameter` reste vrai : c'est le contrat
« même paramètre à `tol` près »).

### 3.5 Orientation et trous

Les wires libres n'ont pas de sens imposé par rapport à la surface :
l'aire signée de leur polygone `(u, v)` décide. Le contour extérieur est
mis en sens direct et les trous en sens indirect en **inversant l'ordre et le
sens de leurs co-arêtes**, jamais la surface ni la face (on ne connaît pas
encore le shell). Chaque trou doit avoir tous ses points d'échantillonnage
dans le polygone extérieur (`HoleOutsideOuter` sinon). Le chevauchement de
deux trous n'est pas testé ici ; il relève de `check()` (PR 6).

### 3.6 Les wires deviennent ceux de la face

Le wire libre passé en argument devient le wire de la face : ses co-arêtes
reçoivent les pcurves et éventuellement un nouvel ordre. Il ne doit donc pas
déjà porter de pcurves ni appartenir à une autre face (`WireInUse`). Pour
construire deux faces qui partagent des arêtes (boîte, sewing), on construit
deux wires sur les mêmes arêtes (`make_wire` sur des listes d'arêtes
communes), chacun recevant ses propres pcurves.

### 3.7 Boucles 2D et annulation

`make_face(srf, boucles2d)` crée pour chaque courbe 2D une arête dont la
courbe 3D est `CurveOnSurface(pc, srf)`, chaîne les boucles par `make_wire`
(fusion des sommets), puis appelle `make_face(srf, wires)`, qui extrait les
pcurves exactement. En cas d'échec à n'importe quelle étape, **toutes** les
entités créées sont effacées (elles occupent les indices au-delà des
capacités relevées au départ) : la garantie forte de la PR 3 est tenue.

### 3.8 Garantie forte de `make_face(srf, wires)`

Toutes les pcurves, orientations et vérifications sont calculées sur des
copies des co-arêtes ; le modèle n'est écrit qu'à la fin, quand plus rien ne
peut échouer.

## 4. Écarts par rapport au document de conception

| Document (§ 6.2) | Implémentation | Raison |
|---|---|---|
| Déroulement du paramètre périodique hors bornes | pas de déroulement ; traversée du seam refusée (`CrossesSeam`) | évaluation hors domaine impossible dans gbs ; découpe au seam = palier 2 |
| Projection par `extrema_surf_pnt` amorcé | Gauss-Newton amorcé, `extrema_surf_pnt` pour l'amorce globale et le secours | rapidité et précision sur des centaines d'échantillons |
| `natural_bounds` détecté si le contour est le rectangle | toujours `false` | optimisation IGES seulement ; `make_face(srf)` couvre ce cas |
| Signature `make_face(srf, wire, holes, opts)` | idem, plus surcharge sans trous et surcharge `std::vector<WireId>` | ergonomie |
| — | `uv_polygon`, `signed_area`, `inside` dans `closure.h` | briques réutilisées par `check()` (PR 6) et `signed_volume` (PR 8) |

## 5. Invariants garantis

Après un succès : chaque wire de la face est fermé, toutes ses co-arêtes ont
une pcurve paramétrée comme l'arête ; pour tout paramètre, `S(pcurve(t))` est à
`edge.tol` du point 3D de la co-arête ; `edge.tol ≥ face.tol` et
`vertex.tol ≥ edge.tol` ; le contour extérieur a une aire signée positive, les
trous une aire négative et sont dans le contour extérieur. Après un échec :
modèle inchangé, y compris pour le chemin « boucles 2D ».

## 6. Tests (`tests_brep_face_wire`, 8 cas, 1461 assertions)

Une fonction commune vérifie, pour chaque face produite, la fermeture des
wires, la présence des pcurves, `S(pcurve(t))` contre le point 3D de la
co-arête sur 41 paramètres, et l'ordre des tolérances.

| Test | Couvre |
|---|---|
| `disk_on_plane` | cercle rationnel sur un plan dont `(u, v) = (x, y)` : le wire devient celui de la face, aire `π`, tolérance d'arête ≤ `pcurve_tol` |
| `outer_reoriented_counter_clockwise` | carré donné en sens horaire : ordre et sens inversés, aire +4 |
| `face_with_holes` | carré avec deux trous circulaires donnés en sens direct : retournés, aires −0,16π et −0,09π |
| `hole_outside_outer_leaves_model_unchanged` | trou hors du contour : erreur, aucune pcurve posée, tolérances inchangées, wires réutilisables |
| `projection_on_curved_surface` | carré `CurveOnSurface` sur une **copie** d'une nappe bicubique : projection, déviation ≤ `pcurve_tol`, pcurve égale au segment 2D d'origine à 1e-4 |
| `exact_pcurves_from_2d_loops` | contour rectangulaire et trou triangulaire en 2D : pcurves = courbes données (mêmes objets), aucune déviation, 7 sommets ; échec ⇒ tout ce qui a été créé est effacé |
| `seam_patches_on_cylinder` | demi-cylindres `[0, π]` et `[π, 2π]` : le côté posé sur le seam prend `u = 0` dans le premier cas, `u = 2π` dans le second |
| `errors` | arc traversant le seam, cercle hors surface, wire ouvert, wire déjà utilisé ou donné deux fois, identifiant invalide, surface nulle |

Les aires des contours circulaires sont comparées à π avec 2048 points par co-arête : avec les 16 points par défaut, le polygone inscrit d'un cercle à une seule co-arête sous-estime l'aire de 3 %, ce qui suffit pour le signe (orientation, trous) mais pas pour une mesure.

`tests_brep_model`, `tests_brep_explore`, `tests_brep_wire` et `tests_brep_face` passent toujours ; `tests_brep_face` utilise désormais la géométrie partagée de `tests/tests_brep_geom.h`.

## 7. Points à valider

1. **Identité de surface par pointeur** pour l'extraction exacte.
2. **Pas de déroulement hors bornes** : une arête ou un wire qui traverse le
   seam est refusé (`CrossesSeam`) ; la découpe au seam viendra au palier 2.
3. **Gauss-Newton** pour la projection courante, `extrema_surf_pnt` en amorce
   et en secours.
4. **Interpolation aux paramètres imposés** de l'arête, raffinement par
   doublement jusqu'à `pcurve_tol` (défaut 1e-5) et 1025 échantillons au plus.
5. **Le wire libre devient le wire de la face** (pcurves posées en place,
   ordre éventuellement inversé) ; un wire ne borne qu'une face.
6. **Orientation par l'aire signée en (u, v)**, sans jamais retourner la face.
7. **Trous** : inclusion testée sur tous les points d'échantillonnage ;
   chevauchement entre trous laissé à `check()`.
8. **`natural_bounds` toujours faux** pour une face construite depuis un wire.
