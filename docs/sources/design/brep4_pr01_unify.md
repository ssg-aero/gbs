# Palier 4 — PR 1 : faces d'un même domaine unifiées après les booléens

| | |
|---|---|
| PR | branche `feat/brep-seamless-faces` |
| Document de référence | [brep_palier4.md](brep_palier4.md), § 0 ter |
| Fichiers | `gbs-brep/unify.h` (nouveau), `gbs-brep/boolean.h`, `gbs-brep/brep`, `gbs/primitives.h`, `python/gbsBindBrepOps.cpp`, `tests/tests_brep_boolean.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `python/tests/test_brep_ops.py`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier4.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Un booléen laissait souvent plusieurs faces sur une région qui n'en demande qu'une :
- les morceaux d'une face découpée par une couture (tore de `torus_box`, cylindres de `steinmetz`) ;
- les faces coplanaires ou coaxiales des deux opérandes (`coplanar_boxes`, `stacked_boxes`, `coaxial_cylinders`) ;
- la face d'une sphère qui contient un pôle, avec une couture et une arête dégénérée (`box_sphere`).

CATIA et SolidWorks donnent une face par région, et la galerie de la PR 0 le montrait. Cette PR fait de même à la fin de chaque booléen. C'est le niveau 1 du § 0 ter du palier, sans changer le modèle.

Les contacts ponctuels sur des surfaces B-splines (§ 0 quater) passent dans une PR à part (PR 1 bis) : ce sont des intersections, sans lien avec l'unification.

## 2. Contenu

### 2.1 `unify_faces` (`gbs-brep/unify.h`)

Les paires de faces adjacentes de l'ensemble sont fusionnées tant qu'une fusion réussit. Ensuite, les morceaux d'arêtes restés sont recollés.

**`merge_faces(a, b)`** fusionne deux faces d'un même domaine :
1. **Domaine.**
   - Même objet surface : la carte des paramètres de b vers ceux de a est l'identité.
   - Primitives de même type : la carte est affine, ajustée par moindres carrés sur des points du bord de b, puis vérifiée sur ces points et entre eux.
     - Deux plans coplanaires donnent un déplacement du plan.
     - Deux cylindres coaxiaux donnent un décalage de l'angle et de la hauteur.
   - La carte est vérifiée à tol près, pas au rayon de confusion : les arêtes de b gardent leur tolérance sur la surface de a. Deux faces à quelques tolérances l'une de l'autre, que les booléens traitent comme confondues (`brep_confusion_factor`), restent deux faces.
   - Les paramètres périodiques des échantillons sont déroulés autour d'une première carte, tirée des dérivées en un point. Sans cela, une face qui fait le tour de la période a des échantillons aux deux bouts de la période, et l'ajustement est faux (le cas des cylindres coaxiaux).
2. **Sens.** Les normales extérieures doivent concorder, orientation dans la coque comprise.
3. **Alignement.** Les pcurves de b passent dans les paramètres de a, ce qui est exact pour les B-splines (une carte affine commute avec la base). Elles sont décalées de périodes entières pour qu'une arête commune soit au même endroit dans les deux faces.
4. **Contours.** Les arêtes communes au même endroit sont retirées. Les co-arêtes restantes sont enchaînées en boucles par leurs sommets et leurs paramètres.
   - Il faut une seule boucle d'aire positive, la boucle extérieure ; les autres sont des trous.
5. **Fenêtre.**
   - Une face qui traverse la couture reçoit une fenêtre du paramètre périodique centrée sur elle, par `ElementarySurface::with_box`, et n'a plus de couture.
   - Une face qui fait juste le tour de la période ouvre sa fenêtre sur elle-même. C'est le cas du côté du cylindre de Viviani, pincé au point double de la courbe.
   - Une face qui dépasse la période est refusée.
   - Sur un paramètre non périodique, la fenêtre s'élargit pour une primitive. Pour des morceaux d'une même surface, elle reste celle de la surface, qui les tenait déjà.

**`to_plane`** traite les plans sous forme B-spline. D'autres systèmes écrivent ainsi leurs plans (OCCT après `BRepBuilderAPI_NurbsConvert`), et c'est aussi la partie plane d'une B-spline.
- Le test est sur le rectangle des paramètres de la face : la surface y est affine, S(u, v) = S₀ + u X + v Y, sur une grille de 7 × 7 points.
- La face passe alors sur la primitive `Plane`, ses pcurves étant transformées exactement. Elle peut ensuite fusionner avec une face du même plan au-delà du rectangle de sa surface.
- La conversion n'a lieu que si une face voisine est coplanaire : une face plane isolée garde sa B-spline.

**`reframe_pole`** traite la face d'une sphère qui contient un pôle, et donc une arête dégénérée et la couture qui y mène.
- Condition : la face et le pôle doivent tenir dans un hémisphère ouvert autour de c, la direction moyenne du bord, avec une marge de 0,05 sur le cosinus.
- La face est alors refaite sur la même sphère, tournée : axe normal à c, couture à l'opposé de c, fenêtre centrée sur c. Elle n'a plus ni couture ni pôle, et `make_face_use` projette ses pcurves à nouveau.
- La boucle extérieure est celle qui tourne dans le sens direct autour de la normale de l'usage, mesurée dans le plan tangent en c (projection centrale, bijective sur l'hémisphère).

**`join_edges`** est l'inverse de `split_edge`.
- Deux arêtes sur la même courbe, bout à bout en un sommet que plus rien d'autre n'utilise, et consécutives avec la même pcurve dans chaque contour, redeviennent une arête.

Les sommets aux bouts de l'arête recollée sont portés au moins à sa tolérance, la plus grande des deux morceaux (ordre des tolérances de `check`).

**`join_seams`** recolle les morceaux d'une couture.
- Deux arêtes qui sont la couture d'une même face sur une primitive (iso-lignes c et c + période), bout à bout en un sommet libre, deviennent une arête sur l'iso-courbe de la primitive, avec pour pcurves les deux iso-lignes. C'est exact.
- Les morceaux viennent de faces fusionnées dont les arêtes portent des courbes différentes, que `join_edges` laisse séparées : trois morceaux sur les cylindres coaxiaux, deux sur le tube du tore.

**`UnifyReport`** donne :
- le nombre de fusions, de faces mises sur un plan, de pôles re-repérés et d'arêtes recollées ;
- l'historique dans l'ordre, sous la forme {nouvelle face, morceau, morceau}. Le second morceau est invalide pour une face seulement mise sur une autre surface.

### 2.2 Booléens

- `BooleanOptions::unify` (vrai par défaut) appelle `unify_faces` sur les faces gardées, après la séparation des arêtes non manifold et avant le regroupement en coques.
- L'origine (opérande A ou B) d'une face fusionnée est celle de son premier morceau. `BooleanResult::faces` suit l'historique.
- Python : argument `unify` de `boolean`, `fuse`, `cut` et `common`.

### 2.3 Comparaison à OpenCascade

La comparaison du palier 2 attendait le même nombre de faces que `BRepAlgoAPI`. Elle compare maintenant au résultat d'OCCT unifié par `ShapeUpgrade_UnifySameDomain`, l'équivalent chez OCCT, qui est plus exigeant. Chaque cas a l'une de ces attentes :
- `same` : autant de faces, le cas général ;
- `fewer` : au plus autant, quand gbs fusionne des morceaux qu'OCCT garde séparés.

Si les nombres de solides diffèrent volontairement (`partial_coincidence`), les faces ne sont pas comparées.

## 3. Résultats

### 3.1 Galerie : nombres de faces sans et avec unification

| Cas | fuse | cut | common |
|---|---|---|---|
| cylindres croisés égaux | 8 → 8 | 7 → 6 | 5 → 4 |
| boîtes coplanaires | 14 → 10 | 8 | 6 |
| boîtes empilées | 12 → 10 | 7 → 6 | — |
| cylindres coaxiaux | 5 → 3 | 3 | 3 |
| cylindre couché sur une boîte | 10 → 9 | 7 → 6 | — |
| tore et boîte | 8 | 9 → 8 | 7 → 6 |

Les 13 cas sur primitives exactes restent valides et fermés, avec les mêmes volumes.

Plus aucune face de la galerie n'est séparée de sa voisine par une couture. Les coutures qui restent sont celles des faces qui font vraiment le tour :
- un anneau de cylindre, le tube du tore ;
- la sphère de `box_sphere` en fuse, qui contient ses deux pôles ;
- le sommet d'un cône.

La sphère de `box_sphere` en cut et en common n'a plus ni couture ni pôle.

### 3.2 OpenCascade

On compare 27 paires de solides et 3 opérations au résultat d'OCCT unifié.

| | Cas |
|---|---|
| autant de faces | tous les autres, dont les cas NURBS (plans B-spline mis sur des plans) et `viviani_poles`, `partial_coincidence` (cut) et `partial_coincidence_rising`, qui avaient plus de faces qu'OCCT avant l'unification |
| moins de faces | `cylinder_lying`, `box_on_edge`, `seam_on_face` et leurs versions NURBS en fuse : la face du dessous d'une ligne de contact reste entière, alors qu'OCCT la laisse en deux ; `viviani` en cut : le côté du cylindre pincé au point double est une seule face, deux chez OCCT |
| plus de faces | aucun |

OCCT relit tous ces résultats comme valides (`BRepCheck_Analyzer`).

### 3.3 Tests

Après fusion de `master` (v0.7.0), les 667 tests C++ passent, comparaison à OCCT comprise. Les tests Python passent aussi : 87, et la limite connue reste en `xfail`.

Assertions ajoutées ou changées :
- empilements de boîtes à plans et de boîtes B-splines : 6 faces ;
- cylindres coaxiaux : 3 faces et 3 arêtes, la couture d'un seul tenant ;
- boîte posée sur une arête, coupée : 6 faces ;
- huitième de sphère au coin d'une boîte : sans couture ni pôle (`seamed_faces`) ;
- Python : `fuse(..., unify=False)` garde les 10 faces de l'empilement ;
- boîtes à quelques tolérances l'une de l'autre (`faces_a_few_tolerances_apart`, venu de `master`) : 12 faces en fuse au lieu de 14. Les morceaux d'une face coupée par les arêtes imprimées à moins de d d'elle sont réunis ; les dessus, et les dessous, à d l'un de l'autre restent deux faces (10 à d = 0). Le résultat est vérifié valide.

## 4. Limites

| Limite | Suite |
|---|---|
| Deux B-splines distinctes de même géométrie ne fusionnent pas, sauf les plans (cylindres B-splines fermés empilés de `coincident_closed_bsplines`) | reconnaissance des primitives sous forme B-spline (cylindre, cône, sphère, tore), comme `to_plane` |
| Une face qui fait le tour garde une couture, et une sphère qui contient ses deux pôles garde ses arêtes dégénérées | niveau 2 du § 0 ter (topologie périodique, palier ultérieur) |
| Une calotte de sphère plus grande qu'un hémisphère garde son pôle | un autre repère existe parfois (axe hors de la face) ; à chercher si un cas le demande |
| Deux faces à quelques tolérances l'une de l'autre ne fusionnent pas | fusion en relevant la tolérance des arêtes, comme `ShapeUpgrade_UnifySameDomain` : refusée ici, la qualité visée ne dégrade pas les tolérances |
| Le recollement ne s'applique qu'aux faces d'un même résultat de booléen | `unify_faces` est public : la découpe (`split`) et la couture (`sew`) pourront l'appeler |

## 5. Points à valider

1. **Unification par défaut** dans les booléens (`unify = true`), désactivable.
2. **Comparaison OCCT contre `ShapeUpgrade_UnifySameDomain`** plutôt que contre le résultat brut de `BRepAlgoAPI`.
3. **Plans B-splines convertis en `Plane`** seulement quand une face coplanaire est voisine.
4. **Contacts ponctuels sur B-splines** sortis dans une PR 1 bis.
