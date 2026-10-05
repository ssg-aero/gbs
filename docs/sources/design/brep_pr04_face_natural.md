# BREP palier 1 — PR 4 : faces à bornes naturelles (`make_face(srf)`, `gbs-brep/closure.h`)

| | |
|---|---|
| PR | branche `feat/brep-face-natural` |
| Document de référence | [brep_core.md](brep_core.md), sections 3.3, 3.4 et 6.1 ; notes précédentes [PR 1](brep_pr01_model.md), [PR 2](brep_pr02_explore.md), [PR 3](brep_pr03_builders_wire.md) |
| Fichiers | `gbs-brep/closure.h` (107 l., nouveau), `gbs-brep/builders.h` (+120 l. : `make_face(srf)`, trois codes d'erreur), `gbs-brep/brep`, `tests/tests_brep_face.cpp` (336 l.), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Construire une face qui couvre **tout** le rectangle paramétrique d'une
surface, avec la topologie correcte quand la surface se referme sur
elle-même : arêtes seam (cylindre, révolution complète, tore) et arêtes
dégénérées (pôles de sphère, apex de cône). C'est l'équivalent de
`BRepBuilderAPI_MakeFace(surface)` d'OCCT, et le cas « face à bornes
naturelles » de l'export IGES (entité 144 avec `N1 = 0`).

## 2. Vue d'ensemble

```
surface_closure(srf, tol, n = 9) -> SurfaceClosure { closed_u, closed_v, degenerate_u1/u2/v1/v2 }
uv_signed_area(m, wire, n = 16)  -> T          (> 0 : sens direct, < 0 : trou)

make_face(m, srf, tol) -> BuildResult<FaceId>

  v2  3 ──────── side 2 (v = v2, u ↓, Reversed) ──────── 2
      │                                                  │
   side 3 (u = u1, v ↓, Reversed)            side 1 (u = u2, v ↑, Forward)
      │                                                  │
  v1  0 ──────── side 0 (v = v1, u ↑, Forward) ───────── 1
      u1                                                 u2

  pcurve d'un côté  : BSCurve<T,2> de degré 1, paramétrée comme l'arête (u pour une iso-v, v pour une iso-u)
  courbe 3D         : CurveOnSurface<T,3>(pcurve, srf), exacte
  côté dégénéré     : make_degenerate_edge (curve = nullptr), pcurve conservée
  direction fermée  : seam = une arête, deux co-arêtes de sens opposés, deux pcurves
  coins confondus   : un seul sommet
```

| Surface | Sommets | Arêtes | dont seams | dont dégénérées | Shell d'une face |
|---|---|---|---|---|---|
| plan | 4 | 4 | 0 | 0 | ouvert, 4 arêtes libres |
| cylindre | 2 | 3 | 1 | 0 | ouvert, 2 cercles libres |
| cône | 2 | 3 | 1 | 1 (apex) | ouvert, 1 cercle libre |
| sphère | 2 | 3 | 1 | 2 (pôles) | **fermé** |
| tore | 1 | 2 | 2 | 0 | **fermé** |
| révolution 2π / π | 2 / 4 | 3 / 4 | 1 / 0 | 0 | ouvert |

## 3. Choix d'architecture et justification

### 3.1 Détection numérique de la fermeture (`surface_closure`)

Conformément au document (§ 3.3, question 7), la hiérarchie `Surface<T,3>`
n'est pas modifiée : la fermeture est mesurée en échantillonnant les quatre
isos de bord (9 points par défaut) et en comparant à `tol`.

- `closed_u` : `max_v ‖S(u1,v) − S(u2,v)‖ ≤ tol` ; idem `closed_v`.
- `degenerate_u1` : `max_v ‖S(u1,v) − S(u1,v1)‖ ≤ tol` ; idem pour les trois autres côtés.
- **Un côté dégénéré n'est jamais un seam** : sur une sphère, les isos `v = v1`
  et `v = v2` sont deux pôles distincts ; si les deux étaient dégénérés au
  même point, on aurait une surface réduite à une courbe.

Avantages : fonctionne pour toutes les surfaces (NURBS, `SurfaceOfRevolution`,
`SurfaceOffset`, surfaces analytiques d'un futur lecteur STEP), sans
connaître leur type. Limite : un échantillonnage peut manquer une
fermeture « presque » vraie entre deux échantillons ; neuf points sur des
isos de surfaces lisses suffisent en pratique, et le paramètre est réglable.

### 3.2 Pcurves de degré 1, courbes 3D `CurveOnSurface`

- Chaque côté a sa pcurve : un segment `BSCurve<T,2>` de degré 1 dont le
  vecteur de nœuds `[t1, t1, t2, t2]` reproduit le paramètre de l'arête (`u`
  pour une iso-v, `v` pour une iso-u). La convention « pcurve paramétrée
  comme l'arête » de la PR 1 est ainsi respectée exactement.
- La courbe 3D est `CurveOnSurface<T,3>(pcurve, srf)` : **exacte** par
  construction (déviation nulle, `same_parameter`), et partage le même
  `shared_ptr` de pcurve. C'est aussi ce que la PR 5 reconnaîtra pour
  extraire une pcurve sans projection.
- Pour une surface NURBS, on aurait pu prendre `isoU`/`isoV`, qui renvoient
  une vraie `BSCurve`. Ce n'est pas fait ici : un seul chemin pour toutes
  les surfaces, et la conversion en NURBS est le travail de l'export IGES
  (PR 9), qui sait faire `isoU`/`isoV` quand la surface est NURBS.

### 3.3 Seams et arêtes dégénérées : sans cas spécial dans le modèle

- **Seam** (`closed_u`) : une seule arête pour `u = u1` (paramètre `v`).
  Le côté 3 l'utilise en `Reversed` avec la pcurve `u = u1`, le côté 1 en
  `Forward` avec la pcurve `u = u2`. Deux co-arêtes, deux pcurves, un `EdgeId`,
  comme prévu au § 2.4 du document. Le shell d'une face cylindrique voit donc
  le seam utilisé deux fois en sens opposés : variété et orienté, comme
  vérifié par `is_orientable`.
- **Arête dégénérée** : `curve = nullptr`, `degenerate = true`, mêmes bornes
  que le côté, et **la pcurve est conservée** sur la co-arête : c'est elle qui
  borne la face dans `(u,v)` et que l'export IGES écrira.
- **Sommets** : les quatre coins sont regroupés s'ils sont à `tol` les uns des
  autres ; la tolérance du sommet couvre l'écart entre les coins regroupés.
  Les seams et pôles en découlent naturellement (sphère : 2 sommets, tore : 1).

### 3.4 Garantie forte et erreurs

Tout est calculé (fermeture, coins, regroupement) avant la première écriture.
Trois codes d'erreur s'ajoutent : `NullSurface`, `UnboundedSurface` (bornes
infinies ou rectangle vide) et `DegenerateSurface` (les quatre côtés
dégénérés : surface réduite à un point).

### 3.5 `uv_signed_area`

Formule du lacet sur 16 points par co-arête, sens des co-arêtes appliqué.
Elle valide ici le sens direct du contour (aire = aire du rectangle), et
servira en PR 5 à orienter les contours extérieurs et les trous d'une face
construite depuis un wire.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| `surface_closure` renvoie aussi les côtés dégénérés | idem, plus la règle « dégénéré ⇒ pas seam » | une sphère aurait sinon un faux seam sur ses pôles |
| Courbe 3D : `CurveOnSurface` ou `isoU`/`isoV` | toujours `CurveOnSurface` | un seul chemin, exact ; conversion NURBS reportée à l'export IGES |
| `make_face(srf)` dans `builders.h` | idem ; `surface_closure` et `uv_signed_area` dans `closure.h` | `closure.h` regroupe les requêtes géométriques sur les faces (et `signed_volume` en PR 8) |

## 5. Invariants garantis

Après un succès, la face vérifie : un seul wire fermé de quatre co-arêtes,
toutes avec pcurve ; aire signée positive et égale à l'aire du rectangle ;
pour chaque co-arête et tout paramètre, `S(pcurve(t))` coïncide avec le point
3D de la co-arête à `edge.tol` près ; tolérance de sommet ≥ tolérance d'arête ;
`natural_bounds = true` ; la surface est partagée (même `shared_ptr`).

## 6. Tests (`tests_brep_face`, 9 cas, 391 assertions)

Les surfaces de test sont analytiques (classe `Analytic` du fichier de test :
cylindre, sphère, tore, cône), plus trois surfaces gbs : plan `BSSurface`,
cylindre **NURBS rationnel exact** (cercle rationnel extrudé) et
`SurfaceOfRevolution` complète et partielle. Pour chaque face, une fonction
commune vérifie les invariants du § 5, compte sommets, arêtes, seams,
arêtes dégénérées, et construit un shell d'une seule face pour tester
`free_edges`, `is_closed` et `is_orientable` de la PR 2.

| Test | Couvre |
|---|---|
| `surface_closure` | les six surfaces : seams, pôles, apex, plan sans fermeture, NURBS rationnel |
| `plane` | 4/4, contour qui part de `(u1,v1)` |
| `cylinder_seam` | analytique et NURBS : 2 sommets, 3 arêtes dont 2 cercles fermés et 1 seam |
| `sphere_poles` | 2 pôles, 2 arêtes dégénérées sans courbe, 1 seam, shell fermé |
| `torus_two_seams` | 1 sommet, 2 seams, shell fermé |
| `cone_apex` | apex dégénéré, seam, cercle de base libre |
| `revolution` | 2π : seam ; π : 4 arêtes ordinaires |
| `shared_surface_and_errors` | surface partagée ; surface nulle, tolérance nulle, surface réduite à un point, bornes infinies, rectangle vide ; modèle inchangé |
| `uv_signed_area_sign` | +1 pour le contour, −1 pour le même contour retourné |

Les suites `tests_brep_model`, `tests_brep_explore`, `tests_brep_wire` passent toujours.

## 7. Points à valider

1. **Fermeture détectée numériquement** (9 échantillons par iso), sans méthode
   virtuelle sur `Surface` (question 7 du document, déjà tranchée ainsi).
2. **Règle « un côté dégénéré n'est jamais un seam »**.
3. **Courbes 3D en `CurveOnSurface`** pour toutes les surfaces, y compris NURBS ;
   conversion en `BSCurve` par `isoU`/`isoV` à l'export IGES.
4. **Pcurves de degré 1 paramétrées comme l'arête** (nœuds `[t1, t1, t2, t2]`).
5. **Arête dégénérée avec pcurve conservée** sur la co-arête.
6. **Ordre du contour** : départ au coin `(u1, v1)` le long de `v = v1`, sens direct.
