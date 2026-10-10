# Lecture STEP — PR 4 : géométrie STEP vers géométrie gbs

| | |
|---|---|
| PR | branche `feat/step-geometry` |
| Document de référence | [step_reader.md](step_reader.md), § 3.1 à 3.3, 4 et 6.2 ; questions 2, 3 et 4 |
| Fichiers | `gbs-io/step/units.h` (nouveau), `gbs-io/step/geometry.h` (nouveau), `tests/tests_step_geometry.cpp` (nouveau), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Traduire les entités géométriques d'un fichier STEP en géométrie gbs, **sans
encore de topologie** :

- repères ;
- courbes : droite, cercle, ellipse, B-splines sous toutes leurs formes,
  polyligne, courbe coupée, courbe composite, courbe sur surface ;
- surfaces : plan, cylindre, cône, sphère, tore, révolution, extrusion,
  B-spline, surface coupée, surface décalée ;
- unités du contexte de représentation (longueur, angle) et tolérance du fichier ;
- bornes des surfaces infinies.

La PR 6 s'en servira pour construire arêtes et faces.

## 2. Vue d'ensemble

```
read_units(file, #contexte, unité_cible_mm) -> Units { length, angle, length_name, angle_name, uncertainty }

GeometryReader r{file, units};
r.frame(#id), r.axis1(#id), r.point<3>(#id), r.direction<3>(#id), r.vector<3>(#id)
r.curve<3>(#id) -> shared_ptr<const CurveDef<3>>      (et curve<2> pour les pcurves)
r.surface(#id)  -> shared_ptr<const SurfaceDef>

CurveDef   value(t), parameter(P), domain(), periodic(), reversed(), kind(),
           to_nurbs(t1, t2) -> {NURBS, u1, u2}, nurbs_parameter(t1, t2, t)
SurfaceDef value(u, v), parameters(P), domain(), periodic(), kinds(),
           to_nurbs(box) -> NURBS, nurbs_parameters(box, u, v)
parameter_box(surface, boucles de points, marge = 0.1) -> UVBox
```

## 3. Choix d'architecture et justification

### 3.1 Deux étages : définition STEP, puis NURBS à la demande

Une droite STEP est infinie, un cercle complet. La NURBS exacte d'une arête ne
peut être construite qu'une fois connue sa **plage**, que donnent ses sommets
(PR 6). De même pour une face sur un plan ou un cylindre : il faut d'abord
connaître ses bords. Le lecteur procède donc en deux étages :

1. **Définition** (`CurveDef`, `SurfaceDef`) : l'entité STEP dans **son propre
   paramétrage**, déjà convertie en unités cibles (longueurs) et en radians
   (angles). Une définition évalue des points et **inverse un point en
   paramètre** :
   - analytiquement pour les droites, les coniques, les surfaces élémentaires,
     la révolution et l'extrusion ;
   - par échantillonnage puis Newton pour une B-spline.
2. **NURBS exacte** (`to_nurbs`) sur une plage donnée, avec les constructions de
   la PR 3. Pour un arc, la NURBS va exactement de t1 à t2, ce qui donne des
   bornes d'arête **sans projection numérique** : t1 et t2 sont les angles des
   sommets. `nurbs_parameter` convertit exactement un paramètre STEP en
   paramètre NURBS (angle → paramètre de l'arc rationnel). Une B-spline est
   rendue telle quelle : son paramètre est déjà celui de STEP.

Les définitions sont des `std::variant` (droite, conique, B-spline, courbe
coupée ; neuf genres de surfaces), mises en cache par `#id`. Une même courbe
partagée par plusieurs entités n'est lue qu'une fois.

### 3.2 Bornes des surfaces infinies et seam (question 4)

`parameter_box` reçoit des points échantillonnés sur les boucles de la face et
applique la recommandation de la question 4 : un rectangle **par face**,
englobant ses bords, élargi de **10 %**. Les directions angulaires reçoivent un
traitement propre :

- **La boucle s'enroule autour de l'axe** (cumul des écarts angulaires de ±2π,
  par exemple une bande de cylindre) : la face fait le tour complet. La boîte
  est alors [0, 2π], et le seam de la NURBS coïncide avec celui de la surface
  STEP, donc avec l'arête de seam du fichier.
- **Sinon**, la boîte est l'étendue angulaire des boucles, qui peut chevaucher 0
  (par exemple [−20°, 30°]), et la marge ne dépasse jamais le tour. **Une telle
  face ne rencontre aucun seam** : une arête qui traversait le seam STEP (de
  350° à 10°) n'en traverse plus aucun. La découpe au seam (PR 7) ne reste
  nécessaire que pour les faces qui font le tour complet.
- Le tour retenu est celui dont le centre est dans [0, 2π[.
- Plusieurs boucles sont ramenées au même tour que la première.
- Les points sur l'axe (pôle de sphère, sommet de cône) n'ont pas d'angle et
  sont ignorés en u.

Dans les autres directions, la boîte est l'étendue des points plus 10 %,
limitée au domaine naturel : latitude de la sphère dans [−π/2, π/2], cône
arrêté à son sommet. Une surface B-spline, ou une surface coupée, garde son
propre domaine.

### 3.3 Unités : converties à la lecture (question 6)

`read_units` lit le contexte de représentation (instance complexe
`GEOMETRIC_REPRESENTATION_CONTEXT` + `GLOBAL_UNIT_ASSIGNED_CONTEXT` +
`GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT`). Il comprend :

- les unités SI préfixées (`.MILLI.,.METRE.`) ;
- les unités par conversion (`'INCH'`, `'DEGREE'`), récursivement ;
- les formes simple et complexe de chaque entité.

Le `GeometryReader` applique les facteurs à la lecture :

| Grandeur | Conversion |
|---|---|
| coordonnées, rayons, norme d'un `vector`, décalage | × facteur de longueur |
| demi-angle du cône, paramètres angulaires des coupes | × facteur d'angle (→ radian) |
| paramètres d'une coupe de surface | selon la nature du paramètre (`ParamKind` : angle, longueur ou sans dimension) |
| paramètre d'une droite, nœuds B-spline | inchangés (sans dimension) |
| entités 2D (pcurves) | **aucune** : elles vivent dans l'espace paramétrique STEP de leur surface |

L'unité cible vient des options (défaut : millimètre, recommandation de la
question 6). `Units` garde le nom de l'unité du fichier pour le rapport et pour
`Model::setUnitScale` (PR 2).

### 3.4 B-splines : poids, grille, nœuds non bornés

- **Rationnelles** : pôles cartésiens et poids convertis en coordonnées
  homogènes (la leçon de l'IGES).
- **Grille de surface** : STEP range les pôles `[u][v]`, gbs u le plus rapide.
- **Bézier, quasi uniforme, uniforme** : nœuds engendrés comme le prescrit
  ISO 10303-42.
- **Nœuds non bornés** (`UNIFORM_CURVE`, courbes et surfaces périodiques
  écrites avec des multiplicités 1 aux bouts, comme les écrit OCCT) : gbs
  attend des nœuds bornés (*clamped*). Le lecteur **insère des nœuds**
  (algorithme de Boehm) aux bornes du domaine jusqu'à la multiplicité p + 1,
  puis retire les pôles devenus inutiles. L'opération est exacte : le test la
  compare à une évaluation de Cox-de Boor sur les nœuds d'origine.

### 3.5 Courbes coupées et composites

- `trimmed_curve` : la coupe par point ou par paramètre suit
  `master_representation`, avec repli sur ce qui est présent. Sur une conique,
  la plage suit `sense_agreement` et peut chevaucher 0 (de 300° à 60° vers
  l'avant, c'est [300°, 420°]).
  - `reversed()` signale une courbe parcourue à rebours de sa base ; la PR 6
    le combinera avec le `same_sense` de l'arête.
  - Une coupe d'une coupe n'est pas prise en charge.
- `composite_curve` : chaque segment borné est converti en NURBS rationnelle,
  coupé et orienté selon son `same_sense`. Les segments sont joints en **une
  seule NURBS** par `gbs::join`. Une discontinuité entre segments supérieure à
  1e-6 (relatif) est une erreur.

### 3.6 Surfaces décalées

Un `offset_surface` d'une surface élémentaire est une surface élémentaire, ce
qui reste exact :

| Base | Décalée de d |
|---|---|
| plan | plan translaté de d le long de la normale |
| cylindre | rayon R + d |
| sphère | rayon R + d |
| tore | petit rayon r + d |
| cône | rayon R + d / cos α, même demi-angle |

Le décalage d'une B-spline est rapporté non pris en charge au premier jet, comme
prévu au § 3.3 du document de conception.

### 3.7 Erreurs

| Exception | Cas |
|---|---|
| `StepUnsupported` (`id()`, `type()`) | entité STEP valide mais non prise en charge : `hyperbola`, `parabola`, `offset_curve_3d`, décalage d'une B-spline, entité inconnue, mauvais type à la place d'un repère |
| `StepError` (`id()`) | contenu invalide : rayon nul ou négatif, direction nulle, référence parallèle à l'axe, nœuds incohérents, poids non positif, segments de composite disjoints, unité inconnue |
| `P21AccessError` (PR 1) | argument d'un mauvais genre, référence pendante |

La PR 6 attrapera ces exceptions **par face** pour l'import partiel
(question 8, recommandation). Toutes portent le `#id` fautif pour le rapport.

### 3.8 Emplacement (question 2)

`gbs-io/step/`, namespace `gbs::step`, à côté de l'IGES (recommandation de la
question 2). Les en-têtes ne dépendent que de gbs, de `p21.h` et de
`gbs/bselementary.h` (PR 3).

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| bornes d'arête par projection des sommets sur la courbe convertie (§ 4.3) | inversion **analytique** sur la définition STEP, puis NURBS sur la plage exacte | exacte, sans projection numérique, et la plage d'un arc peut chevaucher 0 |
| découpe des arêtes au seam (PR 7) pour toute arête traversant le seam | la boîte angulaire d'une face qui ne fait pas le tour évite tout seam | la PR 7 ne concernera que les faces qui font le tour complet |
| `offset_surface` : `SurfaceOffset` existant | surfaces élémentaires exactes ; B-spline non prise en charge | exact, et `SurfaceOffset` n'est pas une NURBS |
| environ 600 lignes | environ 1 500 (en-têtes commentés) | inversion analytique et nœuds non bornés |

## 5. Tests (`tests_step_geometry`, 8 cas, 2 817 assertions)

Extraits STEP écrits à la main, chaque entité comparée à sa formule analytique :

| Test | Couvre |
|---|---|
| `units` | mm et radian, incertitude ; pouce et degré par conversion ; mètre en instance simple ; unité cible pouce ; unité inconnue |
| `lines_and_conics` | droite, cercle, ellipse à deux échelles d'unités ; inversion ; arc de 5 à 8 rad (à cheval sur 2π) exact ; plage de plus d'un tour refusée ; cercle 2D non converti |
| `bsplines` | même courbe en `with_knots`, quasi uniforme, Bézier, uniforme (nœuds non bornés), contre Cox-de Boor ; quart de cercle rationnel en instance complexe ; cubique périodique à nœuds non bornés, fermée |
| `trimmed_composite_and_polyline` | coupes en degrés vers l'avant (à cheval sur 0) et à rebours, coupe par point, polyligne, composite avec un segment inversé, `surface_curve`, `hyperbola` non prise en charge |
| `elementary_surfaces` | plan, cylindre, cône jusqu'au sommet, sphère, tore, tore dégénéré, cylindre décalé, cylindre coupé en degrés, en centimètres : formule STEP, NURBS exacte via `nurbs_parameters`, inversion ; cône décalé le long de sa normale |
| `swept_and_bspline_surfaces` | révolution d'une droite et d'un cercle, extrusion d'un cercle et d'une droite, surface B-spline simple et rationnelle complexe (ordre `[u][v]`, poids) |
| `parameter_boxes` | carreau de cylindre à cheval sur 0, bande qui fait le tour, deux boucles de part et d'autre de 0, calotte sphérique avec son pôle, cône jusqu'au sommet, plan |
| `invalid_content` | rayon négatif, demi-angle hors bornes, rayon nul, direction nulle, référence selon l'axe, nœuds incohérents, entité inconnue, mauvais type, référence pendante |

## 6. Points à valider

1. **Deux étages** : définition dans le paramétrage STEP avec inversion
   analytique, puis NURBS exacte sur une plage.
2. **Boîte des surfaces infinies** par face, marge 10 % (recommandation de la
   question 4, sans surface partagée), et **règle angulaire** : tour complet
   seulement si une boucle s'enroule, sinon une plage qui évite le seam.
3. **Unités converties à la lecture**, cible millimètre par défaut
   (recommandation de la question 6), pcurves 2D non converties.
4. **Nœuds non bornés bornés par insertion** (exact).
5. **Composites joints en une NURBS**, décalages limités aux surfaces
   élémentaires.
6. **`gbs-io/step/`, namespace `gbs::step`** (recommandation de la question 2).
