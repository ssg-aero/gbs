# Lecture et écriture STEP — PR 11 : surfaces décalées (`offset_surface`)

| | |
|---|---|
| PR | branche `feat/step-offset-surface`, sur `feat/brep-offset-surface` (PR ssg-aero/gbs#222, non fusionnée au départ) |
| Documents de référence | [step_reader.md](step_reader.md), § 3.3 (`offset_surface`) ; [step_writer.md](step_writer.md), § 3.2 ; [brep4_pr03_offset_surfaces.md](brep4_pr03_offset_surfaces.md) (surfaces décalées de gbs, point 11 de sa revue) |
| Fichiers | `gbs-io/step/geometry.h`, `gbs-io/step/reader.h`, `gbs-io/step/geometry_writer.h`, `tests/tests_step_geometry.cpp`, `tests/tests_step_write_geometry.cpp`, `tests/tests_step_write.cpp`, `gbs-occt/tests/tests_step_occt.cpp`, `python/tests/test_step.py`, `News.md`, cette note |
| Statut | À valider (§ 8) |

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
- Si Newton ne rejoint pas p à 10⁻³ |d| près, p est plus près d'une autre
  partie de la base que de son pied : le décalé se recoupe. On cherche alors le
  plus proche point du décalé lui-même (racine bornée), par morceaux
  échantillonnés.
- La boîte d'une face sur le décalé d'une B-spline est donc celle de ses points
  (plus 10 %), et non le domaine entier de la B-spline. Le repli est mesuré sur
  la face, pas sur toute la base.

### 2.2 Surface d'une face

`to_surface(boîte)` construit `gbs::offset_surface(base->to_surface(boîte), d,
…, boîte)`, soit une `OffsetSurface` sur la surface de la base (la B-spline
elle-même, ou la NURBS exacte d'une révolution ou d'une extrusion), aux
paramètres de cette surface (`surface_parameters` : ceux de la base). Si
`offset_margin` est au plus 10⁻⁹ sur la boîte, la fonction lève `StepError`
(« folds on the face »). La face échoue alors et le rapport le dit.

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

**Écart assumé avec `primitive_offset`** : elle exige des inégalités strictes,
alors qu'ici le bord, où le décalé dégénère en un point, est admis.

Le code du décalé exact reste propre aux définitions STEP (point 11 de la revue
de la PR 222). Une définition garde le paramétrage STEP pour tout le fichier,
sans boîte, alors que `primitive_offset` reçoit une boîte et peut déplacer le
repère d'un cône selon elle. Les règles sont les mêmes ; elles sont écrites ici
en domaines.

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

## 5. Tests

- `tests_step_geometry` :
  - `offset_cone_keeps_the_parametrization` : un cône décalé vers l'intérieur
    au-delà de son cercle de référence est maintenant lu. Le test vérifie son
    domaine, ses points, sa primitive au sommet, ses paramètres, le refus d'une
    boîte à travers le sommet, et le décalé de ce décalé ;
  - `offset_torus_keeps_its_regular_part` : quatre cas (tore et pomme, vers
    l'extérieur et l'intérieur), avec domaine, périodicité, points, primitive
    et paramètres ;
  - `offsets_of_non_elementary_surfaces` :
    - bases : B-spline (cylindre parabolique), révolution d'une B-spline,
      extrusion oblique d'une ellipse, extrusion d'une B-spline, B-spline
      rognée, décalé d'un décalé ;
    - S + d N par différences finies à 10⁻⁸ près, paramètres à 10⁻⁹,
      `to_surface` ;
    - repli refusé sur le domaine, admis sur une bande loin du sommet ;
    - boîte d'une face par ses points, y compris des points plus près de
      l'autre bout de la B-spline que de leur pied.
- `tests_step_write_geometry.offset_surfaces` :
  - B-spline, cylindre et décalé de décalé écrits en `OFFSET_SURFACE`, relus
    au même paramétrage à 10⁻¹² ;
  - décalés d'une révolution non exacte et d'une sphère NURBS (nœuds doubles)
    approchés.
- `tests_step_write.faces_on_offset_surfaces` : faces naturelles gbs → STEP →
  gbs, avec des `OFFSET_SURFACE` imbriqués. Mêmes points à 10⁻¹², même aire à
  10⁻⁵, `check` valide.
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
   bord (point dégénéré) admis, contrairement à `primitive_offset`.
3. **Pomme décalée prise en charge, citron non.**
4. **Repli d'une base non élémentaire mesuré sur la boîte de la face** (ses
   points, plus 10 %). Une face valide dont la marge de 10 % atteint un pli
   est refusée.
5. **Paramètres d'un point** : pied sur la racine, puis Newton, puis, en cas de
   recoupement, plus proche point du décalé (seuil 10⁻³ |d|). Ce dernier
   recours est plus coûteux, mais il ne sert qu'à un décalé qui se recoupe ou
   à un point loin de la surface.
6. **Écriture** : `OFFSET_SURFACE` seulement sur une base C1 par ses nœuds,
   sinon approximation. Alternative : toujours exact (conforme à l'ISO), mais
   illisible par OCCT sur les NURBS de cercles.
7. **Bases lues** : élémentaires, B-spline, révolution, extrusion, rognée,
   décalée. Le lecteur ne lit pas d'autre surface (`surface_replica`,
   `curve_bounded_surface` restent non prises en charge).
