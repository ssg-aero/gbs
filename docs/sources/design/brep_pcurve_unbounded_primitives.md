# Surfaces élémentaires non bornées dans leurs directions ouvertes (issue #260)

| | |
|---|---|
| PR | branche `fix/brep-pcurve-plane-bounds` (sur master 2b2b0a5) |
| Document de référence | [brep_pr05_face_wire.md](brep_pr05_face_wire.md), § 3.3 (mis à jour) |
| Fichiers | `gbs-brep/builders.h` (`unbounded_box`, message), `gbs-brep/pcurve.h` (`PCurveErrc::BeyondBounds`, pcurve exacte d'une B-spline sur un plan dans la boîte), `gbs-brep/coincident.h`, `python/gbsBindBrep.cpp`, `python/gbsBindSurfaces.cpp` (docstrings), `tests/tests_brep_face_wire.cpp`, `python/tests/test_brep.py`, `News.md`, cette note |
| Statut | À valider (§ 6) |

## 1. Le défaut

`make_face` sur un `Plane` borné acceptait ou refusait une arête qui sort de sa
boîte selon le type de sa courbe : une B-spline (rationnelle ou non) était
acceptée, un `Circle3d` ou une `EllipseCurve` refusé (« edge is not on the
surface within pcurve_tol »). Rejoué sur master, avec un cercle de rayon 0,1 à
moitié hors de `[0, 1]²` et un carré autour :

| Surface | Courbes | Au-delà de la boîte, sur master |
|---|---|---|
| `Plane` | `BSCurve`, `BSCurveRational` | construit |
| `Plane` | `Circle3d`, `EllipseCurve`, `CurveTrimmed`, `CurveComposite`, `LineCurve` | refusé, « edge off surface » |
| `CylindricalSurface` (au-delà de v) | analytiques ou B-splines | refusé |
| plan B-spline | toutes | refusé |

L'exception était donc le seul couple B-spline sur `Plane`.

## 2. La cause

`project_pcurve` (`pcurve.h`) prend, pour une B-spline sur un plan, la pcurve
exacte (`plane_pcurve` : les pôles dans le repère du plan), sans la comparer à
la boîte. Tout autre couple passe par la projection, bornée au rectangle
paramétrique (`project_point`, puis `clamp_poles`) : un point au-delà se
projette sur le bord, à sa distance de l'arête, et l'arête est refusée comme
hors de la surface. La cause donnée dans l'issue est juste ; elle ne disait pas
que c'est la seule exception.

## 3. La règle

Le lecteur STEP lit les droites et les plans comme infinis et donne à chaque
face une primitive dont la boîte vient de son contour (`parameter_box`) ;
`unify` et `thick` font de même (`with_box`). CATIA, SolidWorks et OCCT tiennent
les plans pour infinis, la face bornée par ses seuls contours. La règle retenue
(option A de l'issue) :

- les **directions ouvertes** des surfaces élémentaires ne sont pas bornées : le
  plan en u et v, le cylindre et le cône en v. Leur boîte n'est qu'un rectangle
  de travail ;
- les **directions fermées** (l'angle d'un cylindre, d'un cône, d'une sphère,
  d'un tore) restent bornées par leur boîte ; la sphère et le tore n'ont pas de
  direction ouverte ;
- le **cône** n'est pas prolongé au-delà de son sommet (v = −R / tan a), où il
  passe dans sa seconde nappe : une arête au-delà est refusée ;
- les **autres surfaces** (B-spline, révolution, décalées) restent bornées.

Une arête au-delà d'une borne est refusée pour toutes les sortes de courbes,
avec un message dédié : « edge beyond the bounds of the surface ». Le code reste
`EdgeOffSurface` (le lecteur STEP s'en sert pour réessayer plus large, et
Python liste les codes). L'invariant du § 3.3 de la PR 5 tient : les pcurves
restent dans le rectangle paramétrique de la surface de la face.

## 4. La correction

- **`make_face`** (le chemin commun des faces bâties sur des fils, `make_face`,
  `make_face_use`, `make_face_from_pcurves`) : `detail::unbounded_box` mesure
  les paramètres du contour sur la primitive (64 points par arête, ceux sur la
  surface à `max(tol, pcurve_tol)`) ; s'il sort de la boîte dans une direction
  ouverte, la face est bâtie sur une copie de la primitive (`with_box`) dont la
  boîte l'englobe, élargie d'un dixième de son étendue de chaque côté (comme la
  boîte d'une face lue en STEP), sans passer le sommet d'un cône. Sinon, la
  surface donnée est partagée, comme avant. La surface donnée n'est jamais
  modifiée.
- **`project_pcurve`** : la pcurve exacte d'une B-spline sur un plan est
  vérifiée dans la boîte (257 points) ; au-delà, `PCurveErrc::BeyondBounds`.
  Une projection qui s'arrête sur le bord de la boîte dans une direction qui ne
  se ferme pas, loin de l'arête, donne aussi `BeyondBounds` (sinon
  `OffSurface`, l'arête vraiment hors de la surface).
- **`coincident.h`** traite `BeyondBounds` comme `OffSurface` (une arête qui
  quitte la face n'est pas une coupe).

## 5. Effets et autres opérations

**Partage du pointeur.** Quand la boîte doit grandir, la face ne partage plus le
`shared_ptr` donné ; quand elle suffit, elle le partage. Ce qui en dépend :

- `unify_faces` : deux faces sur des copies d'un même plan restent reconnues
  sur une même surface (`parameter_map` compare d'abord les pointeurs, puis,
  pour deux primitives du même type, calcule l'application des paramètres par
  leur géométrie) ; testé : deux carrés côte à côte, au-delà de la boîte d'un
  même plan, cousus puis fusionnés en une face ;
- les écrivains STEP et JSON partagent les géométries par pointeur : deux faces
  au-delà de la boîte d'un même plan écrivent chacune son `PLANE` (STEP) ou sa
  surface (JSON), le même plan sur deux boîtes ; la relecture est la même
  (testé : aller-retour d'un solide) ;
- `equal` compare les géométries par valeur ; `edit.h` et `full_round.h`
  regroupent par pointeur des surfaces qu'ils copient eux-mêmes.

**Autres constructeurs et opérations.** Les faces naturelles (`make_face(srf)`)
couvrent la boîte ; `sweep` et `sew` en bâtissent. Le lecteur STEP et `unify`
passent à `make_face_use` une primitive déjà bornée à la face (le pointeur est
gardé). Les booléens, la section, la découpe et les congés (dont
`fillet_faces`) bâtissent leurs faces sur la surface des faces données, ou sur
des surfaces qu'ils créent avec leur boîte, sans `make_face` sur un contour
quelconque : pour eux, la règle est celle de `project_pcurve`, la même
maintenant pour toutes les courbes (la B-spline sur un plan n'est plus acceptée
hors de la boîte). Aucun test des opérations ne change.

## 6. Tests

`tests_brep_face_wire` :

| Cas | Attendu | Mesuré |
|---|---|---|
| cercle (rayon 0,1) et carré, dans `[0, 1]²` puis à moitié dehors, sur un `Plane` : `Circle3d`, `EllipseCurve`, cercle NURBS, `CurveTrimmed`, `CurveComposite`, `LineCurve`, segments B-spline | tous construits, valides | construits (master : B-splines seules au-delà) |
| plan : surface partagée dans la boîte, copie au-delà ; boîte de la copie autour du cercle ; plan donné inchangé ; aire en (u, v) celle du même cercle dedans | | vérifiés, à 10⁻¹² |
| cylindre et cône (v dans `[0, 1]`), pièce jusqu'à v = 1,5, courbes analytiques ou B-splines ; cône jusqu'à v = −1 (avant son sommet, v = −1,73) | construits, valides | vérifiés |
| cône au-delà de son sommet ; cylindre borné à l'angle `[0, π]`, pièce aux angles 4 à 5 | refusés | refus « beyond the bounds » (l'angle) |
| plan B-spline `[0, 1]²` (`Circle3d`, cercle NURBS) ; NURBS d'un demi-cylindre au-delà de v (analytique ou B-spline) | refusés, « beyond the bounds of the surface » | vérifiés |
| deux carrés côte à côte au-delà d'un même plan, cousus | fusionnés en une face valide | vérifié |

Python (`test_brep.py`) : le cercle de l'issue, analytique et NURBS, construit
dedans et dehors ; un pavé 3 × 2 × 1 de six faces sur des `gbs.Plane` de boîte
`[0, 1]²`, arêtes en segments B-spline et en `LineCurve` : valide, volume 6 à
10⁻⁹, aller-retour du fichier BREP exact (`equal`) et de STEP (volume à 10⁻⁶).

`tests_brep_full_round.refusals_of_shape` (une extrémité à 88°, la nervure
jusqu'à x = 16,3 sur un plan de boîte ±10), qui casserait si le plan était
borné, reste inchangé.

## 7. Points à valider

1. **La règle** (§ 3), décidée le 26 nov. 2026 (option A de l'issue).
2. **La marge** d'un dixième de l'étendue du contour, comme le lecteur STEP.
   *Recommandé* : la boîte de la face garde de la place pour une pcurve
   interpolée près du bord.
3. **Le code `EdgeOffSurface`** gardé pour le refus au-delà des bornes, avec un
   message dédié, plutôt qu'un code nouveau. *Recommandé* : le lecteur STEP et
   Python en dépendent ; un code à part changerait leur comportement.
