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

## 7. Revue de code (a posteriori, lot 8 de la revue du noyau BREP)

Chaque commentaire a été reproduit par une sonde ; les tests cités sont nouveaux et échouent sur la version précédente, sauf ceux des points 10 et 11, qui y font déborder la pile ou allouer sans borne et n'y ont donc pas été lancés. `tests_step_geometry` compte désormais 19 cas et 6 539 assertions.

| Commentaire | Suite donnée |
|---|---|
| 1. Le facteur du degré écrit sur 10 chiffres (0,01745329252, SolidWorks) met 360° à 2π + 2,04e-11 rad. `TRIMMED_CURVE(cercle, 0., 360., .T.)`, ramené dans [t1 + 1e-12, t1 + 2π + 1e-12[, devenait un arc de 2,04e-11 rad ; le cylindre coupé de 0 à 360° avait u2 − 2π = 2,04e-11, et `to_nurbs` levait « arc: angle span larger than 2 pi » | corrigé : `detail::turn_end` traite un écart à un nombre entier de tours inférieur à `angle_precision` × max(2π, écart) (1e-9 relatif, soit 6,3e-9 rad par tour, ce que couvre l'arrondi du dernier chiffre d'un facteur à 10 chiffres) comme ce nombre de tours, et le tour vaut alors exactement 2π, jamais un peu plus ; `trimmed_curve` et `rectangular_trimmed_surface` s'en servent. Testé (`full_turn_with_a_rounded_degree_factor` : cercle de 0 à 360°, de 90 à 450° et à rebours, cylindre et tore coupés ; des arcs de 359,999° et de 0,001° restent des arcs) |
| 2. `RECTANGULAR_TRIMMED_SURFACE` ignorait `usense` et `vsense` et triait les bornes : sur un cylindre, u de 350° à 10° vers l'avant (20°, à cheval sur le seam) devenait [10°, 350°], le complément ; de même v de 300° à 60° sur un tore | corrigé : dans une direction périodique de la base, la plage va de la première borne à la seconde dans le sens donné, à travers le seam s'il le faut, un tour au plus (ISO 10303-42) : [350°, 370°] ; ailleurs les bornes sont triées. `parameters()` d'une surface coupée rend le tour le plus proche de sa plage. Le sens ne retourne pas le paramétrage (comme OCCT) : la normale reste celle de la base. Testé (`rectangular_trimmed_surface_senses` : quatre combinaisons sur le cylindre, deux sur le tore, un plan) |
| 3. `AXIS2_PLACEMENT_3D` d'axe ±(1, 0, 0) sans `ref_direction` prenait x = (0, 0, 1), alors qu'ISO (`first_proj_axis`) prend (0, 1, 0) : un cercle sur ce repère commençait en (0, 0, 10) au lieu de (0, 10, 0) | corrigé : (0, 1, 0) pour un axe (1, 0, 0) ou (−1, 0, 0). Les autres valeurs par défaut (axe (0, 0, 1), `AXIS2_PLACEMENT_2D` (1, 0), `AXIS1_PLACEMENT` (0, 0, 1)) étaient conformes. Testé (`axis2_placement_defaults`) |
| 4. L'inversion d'un point sur une B-spline partait de 200 échantillons uniformes. Sur une épingle de degré 1 longue de 100 mm, aux branches distantes de 0,005 mm, 85 points sur 400 de la branche de retour s'inversaient sur l'autre branche ; l'arête qu'en tirait `read_step` mesurait 96,74 mm au lieu de 103,26, sous une tolérance relevée à 0,005 | corrigé : p + 2 échantillons par intervalle nodal non vide (64 au moins en tout), puis, depuis les 8 meilleurs minima échantillonnés, Newton protégé par une bissection sur le signe de la dérivée jusqu'à l'échantillon voisin, la dérivée étant prise un ulp à l'intérieur (un coin de degré 1 converge sur son nœud). Le plus proche l'emporte ; les ex æquo (1e-12 relatif) vont au plus petit paramètre : le point de fermeture d'une courbe fermée donne le début du domaine. 0,10 ms par inversion sur une cubique de 1 000 pôles (0,004 ms auparavant, avec 311 inversions fausses sur 1 000). L'arête de l'épingle est lue sans relever de tolérance. Testé (`bspline_inversion_hairpin_and_ties`) |
| 5. `SURFACE_OF_LINEAR_EXTRUSION` d'une B-spline par un vecteur non normal au plan de la courbe : `parameters()` projetait p sur le plan passant par C(a) normal à V, puis prenait le point de la courbe le plus proche dans l'espace. S(0,5 ; 2) donnait (0,263 ; 4,367), à 2,43 mm du point, et la boîte d'une face s'écartait de son bord de 0,88 mm | corrigé : u minimise \|P(C(u) − p)\|, P = I − V Vᵀ / \|V\|², le point le plus proche vu le long de V (`CurveDef::parameter_along`, aussi pour une courbe coupée), puis v suit V. Testé (`extrusion_inversion_with_an_oblique_vector` : B-spline, cercle coupé, B-spline coupée) |
| 6. `SURFACE_OF_REVOLUTION` : u était l'azimut de p compté depuis le premier point de la génératrice hors de l'axe, puis p était ramené dans ce demi-plan. C'est faux dès que la génératrice ne tient pas dans un demi-plan méridien : 2,28 mm d'écart pour une droite gauche (hyperboloïde), 4,24 mm pour une droite qui traverse l'axe ; la boîte d'une face u ∈ [0, 1], v ∈ [0, 3] de l'hyperboloïde avait v ∈ [−0,15 ; 1,65] | corrigé : v minimise la distance dans le demi-plan méridien, (a(v) − a_p)² + (ρ(v) − ρ_p)² (cote et distance à l'axe), par les échantillons et le Newton du point 4 (directement pour une droite non bornée) ; u tourne ensuite le point de la génératrice sur p, et vaut 0 sur l'axe. La révolution est donc inversée numériquement, et non plus analytiquement comme l'annonçait le § 3.1. Testé (`revolution_inversion_skew_generatrix` : hyperboloïde, double cône, cubique gauche, cercle qui traverse l'axe, droite coupée ; boîte de la face) |
| 7. `OFFSET_SURFACE` d'un cône : le décalé gardait l'origine avec le rayon R + d / cos α (§ 3.6) ; c'est la même surface, mais pas le paramétrage d'ISO, S'(u, v) = S(u, v) + d N(u, v). S'(0,4 ; 0) en était à 1,1547 mm, et une coupe par paramètres du décalé tombait à côté (bord v = 0 à z = 0 au lieu de −d sin α = −1) | corrigé : origine déplacée de −d sin α le long de l'axe, rayon R + d cos α, d'où S' = S + d N exactement. Un rayon de référence négatif (décalage vers l'intérieur au-delà du cercle de référence) est rapporté non pris en charge (`StepUnsupported`) et non plus comme une erreur. Testé (`offset_cone_keeps_the_parametrization` : trois décalages, demi-angles positif et négatif, 16 points chacun) |
| 8. `DEGENERATE_TOROIDAL_SURFACE` (petit rayon r > grand rayon R) : `select_outer` était ignoré, le domaine était le tore entier, qui se recoupe, et un point du citron (`.F.`) s'inversait sur la pomme, à 4 et 2,01 mm | corrigé : `TorusDef::part` (pomme pour `.T.`, v ∈ [−t, t] ; citron pour `.F.`, v ∈ [t, 2π − t], où R + r cos v ≤ 0 ; cos t = −R / r) restreint le domaine, où v n'est plus périodique, et choisit la branche de l'inversion (u + π sur le citron). Les deux parts sont exactes : `ToroidalSurface` et `build_torus` acceptent r > R et toute plage de v. Le décalé d'une part n'est pas pris en charge (`StepUnsupported`) : la normale s'inverse sur le citron et la plage de v reste celle de la base. Un grand rayon négatif est une erreur. Testé (`degenerate_torus_apple_and_lemon` : domaine, NURBS, primitive et inversion sur chaque part) |
| 9. `read_units` ne retenait une unité que si elle portait `LENGTH_UNIT` ou `PLANE_ANGLE_UNIT`, ce que n'a pas une instance simple : `SI_UNIT(*,$,.METRE.)` était ignorée et le fichier lu en millimètres (facteur 1 au lieu de 1 000), de même que l'incertitude donnée dans cette unité | corrigé : `detail::unit_kind` classe aussi une instance simple, un `SI_UNIT` par son nom (METRE, RADIAN), un `CONVERSION_BASED_UNIT` par l'unité de son facteur ; les autres unités (stéradian, seconde…) restent ignorées. Testé (`units_as_simple_instances`) |
| 10. Une référence cyclique (`#20=TRIMMED_CURVE('',#20,…)`, un `SURFACE_CURVE` sur lui-même, ou par `COMPOSITE_CURVE`, `RECTANGULAR_TRIMMED_SURFACE`, `OFFSET_SURFACE`) faisait déborder la pile, le cache n'étant rempli qu'après la lecture : `read_step` plantait | corrigé : `curve()` et `surface()` tiennent l'ensemble des #id en cours de lecture et lèvent `StepError(id, "cyclic reference")`, que `read_step` rapporte comme l'échec de l'arête, puis de sa face. Testé (`cyclic_references` : cinq courbes en 2D et en 3D, trois surfaces, lecteur encore utilisable ensuite) |
| 11. `QUASI_UNIFORM_CURVE`, `UNIFORM_CURVE`, `BEZIER_CURVE` (et surfaces) avec moins de p + 1 pôles : `n − p` en `size_t` passait sous zéro et les nœuds s'allouaient sans borne (15 Go) ; de même un `B_SPLINE_CURVE_WITH_KNOTS` de degré 3 à 2 pôles, dont les tailles concordent, et une multiplicité négative ou démesurée | corrigé : un degré d'au moins 1 et au moins degré + 1 pôles sont vérifiés avant tout nœud (`detail::checked_degree`), des multiplicités dans [1, p + 1] et de somme n + p + 1 avant d'allouer ; sinon `StepError`. Testé (`bsplines_with_too_few_poles` : huit courbes, trois surfaces) |
