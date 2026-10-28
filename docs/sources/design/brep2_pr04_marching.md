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

## 7. Revue de code (a posteriori, lot 2 de la revue du noyau BREP)

Chaque commentaire a été reproduit par une sonde ; chaque test ci-dessous
(`tests_intersect_marching`) échoue sur la version précédente. Les détails des
§ 3.1 à 3.3 (marge de 1 %, arrêt sur le départ, 8 tours de densification)
sont remplacés par ce qui suit. Temps inchangés (mêmes machines, passes
alternées) : `tests_brep_boolean` 24,6 s contre 24,7 s, `tests_brep_split`
5,3 s, `tests_brep_section` 0,11 s ; les 11 cas d'origine de
`tests_intersect_marching` 0,24 s au lieu de 0,48 s.

| Commentaire | Suite donnée |
|---|---|
| boucle mince fermée trop tôt, à travers sa largeur : le départ à moins de 1,05 h devant et les paramètres à moins de 1,5 fois le pas suffisaient (crête et plan, ellipse 1,06 × 0,07 le long de la diagonale des boîtes : 2 courbes, 6,50 au lieu de 4,28, écart 0,19) | corrigé : `Marcher::passes_start`, commun aux deux cheminements. Le départ doit être sur le pas : 0 < t·(départ − cur) ≤ \|next − cur\|, à moins de max(10 tol, L × rotation) de la corde, tangentes d'accord (à la rotation du pas plus `max_turn` près) ; testé (`thin_loop`) |
| boucle parcourue 2 à 4 fois quand le départ est près d'un extremum d'un paramètre : le pas qui passe le départ le change à peine | corrigé dans la même fonction : chaque paramètre du départ dans l'intervalle de ceux du pas, élargi de deux fois la flèche admise convertie par `Marcher::rates` (plus grande variation du paramètre par unité de longueur sur sa surface, infinie à un pôle). Un saut de couture, une période, est refusé. Plutôt que \|Δq\| ≤ 2 \|départ − cur\| / \|J\| : sur un long cylindre fin, une courbe le long de la couture a des pas plus longs qu'une période en 3D ; testé (`closing_near_an_extremum` : 402 cheminements autour de l'extremum, 3 cas complets) |
| seconde courbe proche d'une droite prise pour déjà tracée : marge de 1 % de chaque segment, alors que le pas atteint 1/10 de la taille sur une droite (vallée et plan, droites x = ±1e-3 : une courbe de longueur 4 au lieu de 8) | corrigé : la marge est la flèche de chaque pas, gardée par point (`SSPoint::sag`) : celle d'une cubique de Hermite par les tangentes aux bouts, L (a0 + a1) / 6 (une inflexion aussi, qu'une flèche L × rotation / 8 ne voit pas), mesurée au milieu et aux quarts, puis estimée sur la polyligne densifiée (`bound_sags`) ; boîte de la branche élargie de la plus grande. Une tangence qui finit où les surfaces s'écartent garde à son bout la longueur sondée par `parting` : les départs de la bande où l'écart reste dans la tolérance sont sur elle (la marge de 1 % les cachait) ; testé (`close_parallel_curves`) |
| petite boucle manquée selon le paramétrage : le Gauss-Newton de départ, à norme minimale en paramètres, déplace la surface aux grands paramètres (feuille de [−20, 20]² paramétrée sur [0, 1], bosse en (10, 10) ou (15, 15) contre le plan z = 0,2 : aucune courbe) ; une paire plate sans départ était abandonnée | corrigé : inconnues en longueurs (\|J_k\| dx_k), pas le plus court en 3D ; un départ trouvé hors de la paire (boîtes élargies de leur taille) la fait redécouper, la plus grande d'abord, au lieu de l'arrêter. Le découpage préalable de la plus grande quand les tailles diffèrent de plus de 4 fois a été essayé puis écarté : il trouvait des départs à tol / 100 sur des quasi-contacts (plan à 1e-9 sous le bord d'un flanc, au coin d'une couture et d'un bord) d'où le cheminement finit `singular`, et la mise à l'échelle suffit aux cas de la revue ; testé (`small_loop_any_parametrization` : 3 paramétrages × 4 positions, une boucle de même longueur) |
| détour étroit enjambé : pas acceptés sur les seules tangentes aux bouts, polyligne brute gardée pour le dédoublonnage (pente z = 0,3 x et pôle relevé de 0,3 en y = −0,81 : écart 8,2e-2, qui devient la tolérance de l'arête ; relevé de 0,01 : deux courbes identiques) | corrigé : flèche de chaque pas, par Hermite et, au milieu et aux quarts, au point des paramètres interpolés, corrigé sur le plan normal à la corde quand l'écart des surfaces divisé par le sinus de leur angle ne le met pas sur la corde ; pas coupé en deux au-delà de max(tol, L max_turn / 4), ou quand aucun point de la courbe n'est trouvé là. Les points densifiés sont gardés dans `branches`. Densification portée à 16 tours : après un détour suivi à petits pas à côté de longs pas, l'interpolation déborde et 8 tours ne suffisaient pas ; testé (`narrow_detour`, 4 cas). Un détour au milieu d'un long pas droit, ni au milieu ni aux quarts, reste vu par la seule densification |
| faces proches et parallèles : les paires plates parallèles coupées jusqu'au plancher absolu tant que leurs boîtes se touchent (carrés inclinés à 1e-3 : `complete` faux ; cylindres NURBS concentriques de rayons 1 et 1,01 : `complete` faux en 0,22 s, booléen refusé en `NonGeneric`) | corrigé : test de séparation sur les tranches des carreaux (`SurfacePiece::centre`, `normal`, `flatness`) : paire écartée quand tous les pôles de l'un sont au-delà de la tranche de l'autre plus tol, d'un côté (plus sa planéité pour un carreau échantillonné, qui borne sa marge) ; plancher de découpe relatif, max(100 × 10 tol, 1e-4 × taille) ; testé (`close_parallel_faces` : carrés de côté 1 et 100, cylindres à 1e-2 et 1e-3) |
| arrêt du correcteur sur un déplacement absolu de 1e-4 tol, sous l'arrondi des coordonnées loin de l'origine (tol 1e-7 à 3e5 : 11 morceaux, `singular`, 64 faux intervalles de branchement) | corrigé : `Marcher::settled` = max(1e-4 tol, 16 ε (\|p1\| + \|p2\|)) dans les trois correcteurs ; pour les corrections contraintes (plan, paramètre fixé), arrêt aussi quand le résidu ≤ tol / 100 ne décroît plus ; contrôle final gardé, à max(tol / 100, arrondi). Pas pour la recherche de départ à norme minimale : son résidu y décroît lentement vers un quasi-contact, qu'il aurait accepté (fin de tangence de `tangency_limits`) ; testé (`far_from_the_origin` : 3e5 et 1e6 à 1e-7, 1e7 à 1e-6) |
| cercles de rayon inférieur à ~300 tol impossibles à suivre : pas plancher 10 tol, rotation limitée, pas coupé en deux (sphère coupée par son centre : R = 3e-4, des moignons, `singular` ; R = 5e-5, rien) | corrigé : au plancher, un pas convergé dont la flèche h × rotation / 8 ≤ tol / 2 est accepté (pas à un point de branchement) et le pas n'est plus réduit sous le plancher ; un pas rejeté pour sa rotation est réduit à h × 0,8 × max_turn / rotation (entre 1/16 et 0,8) au lieu de moitié ; le plancher est essayé une fois ; testé (`small_circles` : R = 3e-4, 1e-4, 5e-5, longueur 2πR à 1e-5) |
| à vérifier : morceau de 0,75 manqué au-delà de la couture d'un cylindre (cas aléatoires plan/cylindre 36 et 102 de la revue des formes fermées) | reproduit, défaut du cheminement : le petit morceau, entre la couture et la borne en v, part du point de couture, qui est en 3D le bout du grand morceau ; le test « déjà tracé », en 3D seulement, l'écartait. Corrigé : un départ est sur une branche quand ses paramètres sont aussi dans ceux du segment (`Marcher::within`) ; de l'autre côté d'une couture il commence une branche, abandonnée si elle est tout entière en 3D sur les branches tracées (une courbe le long de la couture, tracée des deux côtés) ; un départ près d'un point de branchement ou d'une tangence reste écarté en 3D ; testé (`piece_beyond_a_seam` : les deux cas, mêmes nombres de courbes et longueurs que les formes fermées) |
| à vérifier : hyperbole plan/cône tracée avec un écart de 0,17 et pourtant `complete` (cas aléatoire 811 de la revue des formes fermées) | reproduit, défaut du cheminement : le plan passe à 0,055 du sommet d'un cône dont la boîte le contient (v de −2 à 4, sommet en v = −1,5) ; le dernier pas d'une branche, posé sur la borne v = −2, passait le sommet vers l'autre nappe, accepté sur la seule rotation de la tangente (puis les morceaux étaient retracés : 7 courbes, 40,2 au lieu de 13,5). Corrigé avec le commentaire 5 : la pose sur un bord obéit aux mêmes contrôles que les autres pas (creux du sinus, flèche ; pas de point de la courbe au milieu de la corde : refus). Trois courbes, 13,530827 contre 13,530841 par la formule ; testé (`landing_across_an_apex`) |
| CI : sous Windows (MSVC, clang), la sphère dans un cylindre de son rayon de `singular_points_and_limits` sortait `singular` ; les arrondis de la bibliothèque mathématique en décidaient | reproduit sous macOS en déplaçant repères, origines et rayons de 1e-16 à 1e-12 et de quelques ulps : 53 cas mauvais sur 1000 (21 sur 1000 en placements quelconques ; 171 et 62 sur la version d'avant la revue). Trois décisions limites, aucune dans `passes_start` (épinglé par `closing_on_a_step_end` : un départ sur le bout d'un pas, à un arrondi près, est passé par un des pas voisins et par aucun autre) : (1) le Gauss-Newton de départ sur la tangence, au résidu d'arrondi, erre le long d'elle (dernière valeur singulière juste au-dessus de la troncature) sans converger, les paires sont coupées jusqu'au plafond : il s'arrête aussi quand le résidu à l'arrondi ne décroît plus ; (2) la pose de la tangence sur la couture gardait une marge de 1 % du pas, sous sa flèche : la tangence était tracée deux fois ; les pas de `TangentMarcher` prennent la flèche de Hermite, et `bound_sags` reporte aux bouts de la polyligne l'angle de l'autre bout ; (3) un départ au coin des deux coutures, hors des boîtes une fois posé sur la tangence, devenait un contact isolé avant que la tangence soit tracée : les départs qui ne bougent pas sont décidés après toutes les branches, et rien n'en reste s'ils sont sur une branche. Un départ tangent trouvé loin de sa paire n'en fait plus redécouper la paire. 0 sur 1000 dans les deux familles ; testé (`tangent_sleeve_perturbed` : cinq cas qui échouaient, exacts, et 40 autres) |
| hors du lot : plan à 1e-8 sous le sommet d'une bosse (boucle de 61 tol, `p7_tiny_loop`) : quatre contacts et `singular`, comme avant | non traité : les surfaces y sont dans la tolérance sur toute la boucle, une tangence à la tolérance plutôt qu'une boucle |
