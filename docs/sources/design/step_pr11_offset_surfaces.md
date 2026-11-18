# Lecture et écriture STEP — PR 11 : surfaces décalées (`offset_surface`)

| | |
|---|---|
| PR | [ssg-aero/gbs#223](https://github.com/ssg-aero/gbs/pull/223), branche `feat/step-offset-surface`, partie de la tête de la PR ssg-aero/gbs#222 (4b48029, fusionnée depuis dans `master`) |
| Documents de référence | [step_reader.md](step_reader.md), § 3.3 (`offset_surface`) ; [step_writer.md](step_writer.md), § 3.2 ; [brep4_pr03_offset_surfaces.md](brep4_pr03_offset_surfaces.md) (surfaces décalées de gbs, point 11 de sa revue) |
| Fichiers | `gbs-io/step/geometry.h`, `gbs-io/step/reader.h`, `gbs-io/step/geometry_writer.h`, `gbs-io/step/writer.h`, `gbs/offsets.h` (règles communes), `tests/tests_step_geometry.cpp`, `tests/tests_step_read.cpp`, `tests/tests_step_write_geometry.cpp`, `tests/tests_step_write.cpp`, `gbs-occt/tests/tests_step_occt.cpp`, `python/tests/test_step.py`, `docs/sources/design/step_reader.md` (§ 3.3, § 13), `News.md`, cette note |
| Statut | À valider (§ 8) ; revue de code au § 9 |

## 1. Objectif

Jusqu'ici, le lecteur ne lisait un `offset_surface` que sur une surface
élémentaire, dont il recodait le décalé exact ; toute autre base levait
`StepUnsupported("OFFSET_SURFACE of a non elementary surface")`, et la face
était perdue. L'écrivain approchait toute surface décalée.

Cette PR :

1. lit un `offset_surface` sur une B-spline, une surface de révolution, une
   extrusion, une surface rognée ou un autre décalé, en une
   `gbs::OffsetSurface` (`gbs/offsets.h`, PR 222) ;
2. aligne le décalé exact des bases élémentaires sur les règles de
   `gbs::primitive_offset` (repli mesuré sur la face, et non plus refus ou
   acceptation sur le seul rayon de référence) ;
3. écrit `offset_surface` pour une face sur une `OffsetSurface` dont la base
   s'écrit exactement.

## 2. Définition d'un décalé (`OffsetDef`)

ISO 10303-42 : S'(u, v) = S(u, v) + d N(u, v), N = S_u × S_v / |S_u × S_v|,
dans les paramètres de la base S.

```cpp
struct OffsetDef
{
    std::shared_ptr<const SurfaceDef> basis; Real d;     // l'entité
    std::shared_ptr<const SurfaceDef> root;  Real total; // la base non décalée, la somme des distances
    std::shared_ptr<const SurfaceDef> exact; Real dv;    // racine élémentaire : son décalé exact, v + dv
    UVBox domain; std::array<bool, 2> periodic;          // la partie de la base où le décalé est régulier
};
```

- **Décalé d'une surface rognée** : lu comme la surface rognée du décalé
  (`TrimmedSurfaceDef` du décalé de la base, mêmes bornes). Mêmes points, même
  normale : le lecteur garde à `rectangular_trimmed_surface` le paramétrage de
  sa base (les sens ne choisissent que l'arc d'une direction périodique), comme
  OCCT.
- **Décalé d'un décalé** : celui de la racine à la somme des distances (les
  normales coïncident tant que le premier ne se replie pas). Le repli de chaque
  niveau reste contrôlé à la construction de la surface (§ 2.2).
- **Normale** (`SurfaceDef::normal`) : `gbs::normal_direction` sur
  `detail::DefSurface`, la définition vue comme une surface gbs dans ses
  paramètres STEP, bornée par son domaine, aux dérivées premières analytiques
  (`SurfaceDef::tangent`, pour chaque type). Elle donne la limite aux pôles et
  aux sommets. Celle d'un décalé est celle de sa racine.

### 2.1 Points et paramètres

- Valeur : `exact` en (u, v + dv) pour une racine élémentaire, sinon
  S(u, v) + total · N(u, v).
- `parameters(p)` sert à la boîte d'une face et aux coutures. Pour une racine
  non élémentaire, on prend d'abord le **pied de p sur la racine**. Les
  inversions existantes le donnent : pour une révolution, le plus proche point
  dans le demi-plan méridien ; pour l'extrusion d'une B-spline, la courbe vue
  le long de l'axe. Pour une racine B-spline, c'est son plus proche point
  (`global_point_surface`, morceaux de Bézier). Suit un **Newton sur le
  décalé** (`newton_point_surface` sur `OffsetSurface(DefSurface(racine),
  total)`), qui rend le résultat exact là où l'inversion de la racine ne donne
  pas le pied (extrusion oblique d'une ellipse).
- Si Newton ne rejoint pas p à max(10⁻³ |d|, 100 fois la précision du
  fichier) près, p est plus près d'une autre partie de la base que de son
  pied : le décalé se recoupe. On cherche alors le plus proche point du décalé
  lui-même (racine bornée), par morceaux échantillonnés. Le second terme évite
  cette recherche pour un point d'arête resté dans la tolérance de son fichier,
  hors de la surface : avec un petit décalé (d = 0,001), chaque inversion la
  lançait (49 µs par point au lieu de 9 sur une B-spline, 130 au lieu de 4,6
  sur une révolution ; 12 et 8 µs maintenant).
- La boîte d'une face sur le décalé d'une surface non élémentaire, rognée ou
  non, est donc celle de ses points (plus 10 %), et non le domaine entier de la
  base ou de la rognée. Le repli est mesuré sur la face, pas sur toute la base.

### 2.2 Surface d'une face

`to_surface(boîte, étendue)` construit une `OffsetSurface` sur la surface de
la base (la B-spline elle-même, ou la NURBS exacte d'une révolution ou d'une
extrusion), aux paramètres de cette surface (`surface_parameters` : ceux de la
base), par `gbs::offset_surface`. Le repli est mesuré sur l'**étendue** de la
face : la boîte de ses points sans la marge de 10 % (`parameter_box` la rend
aussi), dans le tour de la boîte pour une direction périodique, la boîte
entière dans une direction où la face fait le tour. Si `offset_margin` y est au
plus 10⁻⁹, la fonction lève `StepError` (« folds on the face »). La face échoue
alors et le rapport le dit ; les sommets de pôle créés pour sa couture sont
effacés. La marge sert au paramétrage de la surface, pas au contrôle : une face
valide dont la marge atteint un pli est lue.

Pour un décalé de décalé, la surface du premier niveau est construite (et son
repli contrôlé) d'abord ; le second est l'`OffsetSurface` de la surface de la
racine à la somme des distances, dont seul le repli à cette somme est mesuré.
Chaque niveau est contrôlé une fois (0,37 ms au lieu de 0,54 pour deux
niveaux sur une B-spline).

Une face sur le décalé d'une révolution sans arête de couture reçoit pour
couture l'iso de la surface. Cette surface n'étant ni primitive ni NURBS,
l'iso est une `CurveOnSurface` (segment de pcurve) sur l'`OffsetSurface`
(`reader.h`).

`to_nurbs` n'a pas de NURBS exacte d'un décalé non élémentaire : il lève
`StepUnsupported`. Seul `to_surface` l'appelle, et pour les autres définitions.

## 3. Bases élémentaires : les règles de `primitive_offset`

| base | avant | maintenant |
|---|---|---|
| plan | plan déplacé | inchangé |
| cylindre, sphère | rayon r + d > 0, sinon `StepError` | inchangé |
| tore | tore (R, r + d) entier, même quand son cercle des centres R + (r + d) cos v traverse l'axe (tore fuseau, replié) ; tore dégénéré : non pris en charge | **partie régulière** : les v où R + r cos v > 0 et R + (r + d) cos v > 0, soit \|v\| ≤ t, cos t = max(−R/r, −R/(r + d)) ; v n'est plus périodique si ce n'est pas le tour entier. Pomme décalée : prise en charge. Citron : non (la normale y va vers l'intérieur du tube ; `primitive_offset` le refuse aussi) |
| cône | rayon de référence R' = R + d cos α ≥ 0, sinon `StepUnsupported`, quelle que soit la face | **partie régulière** : les v où ni R + v tan α ni R' + v tan α n'est négatif. Si R' < 0, le cône a son repère au sommet du décalé, et v est décalé de dv = R' / tan α |

Le **domaine** du décalé est cette partie régulière. `parameter_box` y borne
la marge de 10 % d'une face, sans faux refus près du pli. `to_surface` lève
`StepError` pour une boîte qui en sort (le bord compris dans la partie : le
sommet du cône décalé, le point du tore sur l'axe, comme une face de cône
jusqu'à son sommet).

**Bornes d'une surface rognée.** STEP voit le décalé d'un tore fermé en v, même
quand seule une partie en est régulière. Les bornes d'une
`rectangular_trimmed_surface` au-dessus du décalé (`SurfaceDef::closed`, la
périodicité de la racine), ou celles de la rognée sous le décalé, sont donc
données dans les tours du tore : v de 5,5 à 6,2 pour [−0,78 ; −0,08], de 5,5 à
0,5 (sens `.T.`) à travers la couture. Elles sont lues avec cette périodicité,
puis ramenées dans le tour de la partie régulière (`detail::into_turn`), où sont
les paramètres du décalé. Avant, elles restaient dans [0, 2π) (face refusée),
ou étaient triées (l'arc complémentaire).

Les règles viennent de `gbs/offsets.h`, communes avec `primitive_offset` :
`torus_offset_cos_min(R, r, d)` (cos t) et `cone_offset_apexes(R, α, d)` (les
sommets du cône et du décalé au même v). **Écart assumé avec
`primitive_offset`** : elle exige une boîte strictement dans la partie (pour le
cône, le sommet du cône admis, pas celui du décalé), alors qu'ici le bord, où
le décalé dégénère en un point, est admis.

Le décalé exact lui-même reste propre aux définitions STEP (point 11 de la
revue de la PR 222). Une définition garde le paramétrage STEP pour tout le
fichier, sans boîte, alors que `primitive_offset` reçoit une boîte et décale le
v d'un cône selon elle ; ici le repère se déplace et v reste.

## 4. Écriture (`geometry_writer.h`)

- Une `OffsetSurface` dont la base s'écrit exactement (primitive, B-spline, ou
  décalé de l'une d'elles) **et est C1 par ses nœuds** devient
  `OFFSET_SURFACE('', base, d, .F.)`. Ses pcurves sont celles de la base,
  puisque c'est l'espace paramétrique d'`offset_surface`.
- Sur une base dont un nœud intérieur est répété autant de fois que le degré,
  le décalé est approché comme avant, et compté. C'est le cas des NURBS d'un
  cercle, d'une sphère ou d'une révolution, G1 mais pas C1. La raison : OCCT
  refuse un tel `offset_surface` (`Geom_OffsetSurface`, « Offset with no C1
  Surface » quand la réduction des nœuds échoue ; à la lecture, « Surface has
  not been created » et la face est perdue). On l'a constaté sur la révolution
  décalée, lue par gbs puis réécrite. L'objectif 1 de l'écrivain est d'être lu
  par CATIA, NX, Creo, SolidWorks et OCCT.
- Le décalé d'une base citron (partie intérieure d'un tore qui se recoupe) est
  approché aussi : gbs ne lit pas un tel `offset_surface` (§ 3), et la face se
  perdait à l'aller-retour gbs → STEP → gbs.
- Une base est écrite une fois par pièce, comme les surfaces des faces
  (`GeometryWriter::surface(s, ids)`, le cache de la pièce) : une face sur une
  B-spline et la face sur son décalé partagent l'entité, le décalé d'un décalé
  s'écrit sur le premier décalé.

## 5. Tests

- `tests_step_geometry` :
  - `offset_cone_keeps_the_parametrization` : un cône décalé vers l'intérieur
    au-delà de son cercle de référence est maintenant lu. Le test vérifie son
    domaine, ses points, sa primitive au sommet, ses paramètres, le refus d'une
    boîte à travers le sommet, et le décalé de ce décalé ;
  - `offset_torus_keeps_its_regular_part` : quatre cas (tore et pomme, vers
    l'extérieur et l'intérieur), avec domaine, périodicité, points, primitive
    et paramètres ; surfaces rognées sous et sur le décalé (bornes de 5,5 à
    6,2 et de 5,5 à 0,5 ramenées dans [−t, t]), une rognée au-delà de t
    refusée ;
  - `offsets_of_non_elementary_surfaces` :
    - bases : B-spline (cylindre parabolique), révolution d'une B-spline,
      extrusion oblique d'une ellipse, extrusion d'une B-spline, B-spline
      rognée, décalé d'un décalé ;
    - S + d N par différences finies à 10⁻⁸ près, paramètres à 10⁻⁹,
      `to_surface` ;
    - repli refusé sur le domaine, admis sur une bande loin du sommet ;
    - boîte d'une face par ses points, y compris des points plus près de
      l'autre bout de la B-spline que de leur pied ;
    - repli mesuré sur l'étendue de la face, non sur sa marge (face jusqu'à
      v = 0,21, pli au-delà de 0,2215, boîte jusqu'à 0,231), et sur la face
      d'un décalé de B-spline rognée, non sur toute la rognée.
- `tests_step_read.faces_on_offset_surfaces` : bande sur le décalé d'une
  révolution entre deux cercles, sans arête de couture (couture en
  `CurveOnSurface` sur l'`OffsetSurface`, aire du tronc de cône décalé à 10⁻⁵) ;
  calotte sur un décalé replié d'un dôme refusée, sans sommet laissé à son pôle.
- `tests_step_write_geometry.offset_surfaces` :
  - B-spline, cylindre et décalé de décalé écrits en `OFFSET_SURFACE`, relus
    au même paramétrage à 10⁻¹² ;
  - décalés d'une révolution non exacte, d'une sphère NURBS (nœuds doubles) et
    d'un tore citron approchés.
- `tests_step_write.faces_on_offset_surfaces` : faces naturelles gbs → STEP →
  gbs, avec des `OFFSET_SURFACE` imbriqués. Mêmes points à 10⁻¹², même aire à
  10⁻⁵, `check` valide. Dans une pièce, une face sur la B-spline et les deux
  faces décalées : une seule `B_SPLINE_SURFACE_WITH_KNOTS`, deux
  `OFFSET_SURFACE`.
- `tests_offsets` (inchangé) couvre `primitive_offset` réécrite sur les règles
  communes.
- OCCT, `tests_step_occt.offset_surfaces` :
  - face sur un `Geom_OffsetSurface` de B-spline (d = 0,3 et −0,25), écrite par
    `STEPControl_Writer` ;
  - lue par gbs sur une `OffsetSurface` : écart aux points d'OCCT < 10⁻¹²,
    aire à 10⁻⁵, face valide ;
  - réécrite par gbs en `OFFSET_SURFACE`, relue par OCCT : valide,
    `Geom_OffsetSurface`, même aire.
- OCCT, `tests_step_occt.offset_solids` : solides décalés par
  `BRepOffsetAPI_MakeOffsetShape::PerformByJoin`, par arcs ou par
  intersections ; on compare les comptes, les volumes (à 2·10⁻³) et la
  validité, avec double aller-retour. Les faces portent sur :
  - un loft (décalé de `B_SPLINE_SURFACE_WITH_KNOTS`) ;
  - une révolution (décalé de `SURFACE_OF_REVOLUTION`, réécrit approché,
    § 4) ;
  - un prisme (décalé de `SURFACE_OF_LINEAR_EXTRUSION`).
- pytest `test_offset_surface_round_trip`.

## 6. Constats sur OCCT (boîte noire)

- `PerformBySimple` couvre les jeux des arêtes vives par des tolérances de 0,2
  à 0,7, inexploitables (le lecteur de gbs refuse au-delà de `max_tolerance`).
  `PerformByJoin` donne des solides valides (tolérances de l'ordre de 10⁻⁵).
- `BRepGProp` se trompe sur une face `Geom_OffsetSurface` d'extrusion presque
  repliée (profil de la PR 8, marge 1 − d k = 0,106).
  - Son aire vaut 26,92 au lieu de 29,00 : l'aire polyédrique des points
    d'OCCT lui-même vaut 28,999, comme l'intégrale de gbs. Le volume est faux
    de 3 %.
  - Ses dérivées `D1` sont justes (écart relatif 7·10⁻⁹) : c'est l'intégration
    qui pèche.
  - D'où, pour les solides décalés, un profil plus doux (rayon de courbure
    2,6 au moins).
- `Geom_OffsetSurface` refuse une base C0 par ses nœuds qu'elle ne peut rendre
  C1 en en retirant (§ 4).

## 7. Écarts par rapport aux documents de conception

- `step_reader.md`, § 3.3 : `offset_surface` sur « `SurfaceOffset` existant, ou
  rapport non supporté ». Il est maintenant lu sur toute base lue par ailleurs.
- `step_writer.md`, § 3.2 : les décalés étaient approchés. Ils sont maintenant
  exacts quand la base l'est et est C1 par ses nœuds.

## 8. Points à valider

1. **Décalé d'une surface rognée** lu comme la rognée du décalé ; **décalé
   d'un décalé** lu comme celui de la racine à la somme des distances, le repli
   de chaque niveau contrôlé à la construction.
2. **Partie régulière comme domaine** d'un tore ou d'un cône décalé, au lieu
   d'accepter un tore replié ou de refuser un cône valide loin de son sommet ;
   bord (point dégénéré) admis, contrairement à `primitive_offset`, avec qui
   les règles sont partagées (`gbs/offsets.h`). Les bornes d'une rognée sont
   lues dans les tours de la racine, puis ramenées dans la partie régulière.
3. **Pomme décalée prise en charge, citron non** ; l'écrivain approche le
   décalé d'un citron.
4. **Repli d'une base non élémentaire mesuré sur l'étendue de la face** (ses
   points, sans la marge de 10 % de sa boîte), rognée ou non.
5. **Paramètres d'un point** : pied sur la racine, puis Newton, puis, en cas de
   recoupement, plus proche point du décalé (seuil max(10⁻³ |d|, 100 fois la
   précision du fichier)). Ce dernier recours est plus coûteux, mais il ne sert
   qu'à un décalé qui se recoupe ou à un point loin de la surface.
6. **Écriture** : `OFFSET_SURFACE` seulement sur une base C1 par ses nœuds, ni
   citron, sinon approximation. Alternative : toujours exact (conforme à
   l'ISO), mais illisible par OCCT sur les NURBS de cercles.
7. **Bases lues** : élémentaires, B-spline, révolution, extrusion, rognée,
   décalée. Le lecteur ne lit pas d'autre surface (`surface_replica`,
   `curve_bounded_surface` restent non prises en charge).

## 9. Revue de code

Douze points, tous traités dans cette PR.

| # | Constat | Correction |
|---|---|---|
| 1 | Tore décalé restreint à [−t, t] sous ou sur une rognée : bornes laissées dans [0, 2π) (face refusée, paramètres hors du domaine) ou triées (arc complémentaire) | `SurfaceDef::closed`, `detail::into_turn` (§ 3) |
| 2 | Décalé d'une surface rognée non élémentaire : repli mesuré sur toute la rognée | `parameter_box` prend la boîte des points (§ 2.1) |
| 3 | Repli mesuré sur la boîte avec sa marge de 10 % : face valide refusée près d'un pli | étendue de la face, rendue par `parameter_box` (§ 2.2) |
| 4 | `OFFSET_SURFACE` d'un citron écrit, refusé à la relecture par gbs | approché (§ 4) |
| 5 | Repli levé après la préparation de la couture : sommets de pôle laissés dans le modèle | effacés (`reader.h`) |
| 6 | Recherche globale lancée par la tolérance du fichier pour un petit décalé (jusqu'à 28 fois plus lent) | seuil max(10⁻³ \|d\|, `OffsetDef::tol`) (§ 2.1) |
| 7 | Règles de partie régulière écrites deux fois (lecteur, `primitive_offset`) | `torus_offset_cos_min`, `cone_offset_apexes` (§ 3) |
| 8 | Décalé de décalé : repli du premier niveau mesuré deux fois | chaque niveau une fois (§ 2.2) |
| 9 | Base d'un décalé réécrite à chaque usage | cache de la pièce (§ 4) |
| 10 | Couture en `CurveOnSurface` sur un décalé sans test | `tests_step_read.faces_on_offset_surfaces` |
| 11 | Dérivées d'un décalé non élémentaire par différences d'ordre 1, code non atteint | celles de l'`OffsetSurface` de la racine ; la normale d'une rognée est celle de sa base |
| 12 | Lambda `at` masqué dans la couture | renommé `uv_at` |
