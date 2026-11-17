# Palier 4 — PR 1 bis : contacts ponctuels sur des surfaces B-splines

| | |
|---|---|
| PR | branche `feat/brep-bspline-point-contact` |
| Document de référence | [brep_palier4.md](brep_palier4.md), § 0 quater |
| Fichiers | `gbs/intersect/marching.h`, `tests/tests_intersect_marching.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `python/tests/test_examples.py`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier4.md`, `News.md`, cette note |
| Statut | À valider |

## 1. Objectif

Une sphère NURBS posée sur une boîte donnait `NonGeneric` sur les trois
booléens, alors que la même sphère exacte passait. Entre primitives, le point de
contact vient des formes fermées ; entre B-splines, il passe par le cheminement.
Un modèle converti en NURBS (un STEP d'un autre système, la comparaison à OCCT)
échouait donc là où le même modèle en primitives passait. C'était la dernière
limite connue de la galerie (`nurbs_sphere_touching`, en `xfail`).

## 2. Diagnostic

Le cheminement **trouvait** le point de contact, à 4·10⁻⁸ près : la recherche des
départs y converge (Gauss-Newton sur S₁ = S₂). Mais au point de tangence les
normales sont parallèles, la marche lancée de ce départ s'arrête aussitôt
(`End::singular`), et le résultat entier était déclaré singulier. Le booléen
refuse une intersection singulière.

Le point d'arrêt peut être de deux natures :
- un **contact isolé** : les surfaces se touchent sans se traverser (une sphère sur
  un plan, dans une cuvette) ; rien n'est singulier ;
- un **point de branchement** : deux courbes d'intersection s'y croisent (une selle
  sur son plan tangent) ; le refus est justifié.

Le commentaire de `tangent_branch` notait déjà la différence (« both are of the
same size at an isolated contact […], which is not traced »), sans en tirer de
conclusion.

## 3. Contenu

### 3.1 Certification par les secondes formes (`TangentMarcher::isolated_contact`)

Près du contact, l'écart entre les surfaces le long de la normale commune vaut
(l₁ s₁² + l₂ s₂²) / 2, où l₁ et l₂ sont les valeurs propres de II₁ − II₂, la
différence des secondes formes fondamentales dans le plan tangent commun.

- **Définie** (l₁ et l₂ de même signe) : l'écart croît dans toutes les directions,
  le contact est isolé. C'est la certification.
- **Indéfinie** : deux courbes se croisent, point de branchement ; reste singulier.
- **Presque dégénérée** (|l₁| ≤ `angular_tol` |l₂|) : c'est le domaine de
  `tangent_branch`, qui trace une courbe de tangence. Les deux tests sont
  complémentaires, avec le même seuil.
- **Nulle** (|l₂| sous le bruit des courbures, comme `direction`) : contact
  dégénéré, reste singulier.

Le point est ensuite mis sur la tangence par Gauss-Newton sur les quatre
équations de `TangentMarcher::correct` (pieds et normales parallèles). Elles sont
de rang 4 en un contact isolé : la convergence est quadratique, et le point
tombe à 10⁻¹⁵ près au lieu de 4·10⁻⁸. La forme doit rester définie au point
corrigé, et celui-ci dans les boîtes.

`direction` et `isolated_contact` partagent maintenant le calcul de II₁ − II₂
(`form_difference`) ; `direction` est inchangée.

`isolated_contact` rend aussi la **portée** du contact, √(2 tol / |l₁|) : la
distance jusqu'à laquelle les surfaces restent à la tolérance l'une de l'autre
autour de lui (l'écart croît au moins comme |l₁| s² / 2). Elle sert au § 3.3.

### 3.2 Au pôle d'une surface

Une sphère posée sur un plan, son axe vertical, touche par son pôle. Là, S_u = 0 :
la seconde forme n'est pas définie dans les paramètres, et la certification
échouait.

- II est une propriété de la surface, continue à travers un pôle régulier. Une
  surface dont les paramètres dégénèrent au point (plus petite valeur singulière
  de (S_u, S_v) sous 10⁻⁶ de la plus grande) a sa forme prise à des paramètres
  rapprochés du milieu de sa boîte (la boîte de l'intersection) de 10⁻⁴ de leur
  distance à lui, soit 10⁻⁴ de la demi-étendue depuis une borne.
- À travers les pôles, chaque forme est prise sur **sa propre normale**, orientée
  comme la normale commune, et non seulement au point déplacé. Près d'un pôle,
  |S_u| est du même ordre que l'écart d'angle entre les deux normales : avec la
  normale de l'autre surface, la forme le long de u sortait fausse. C'est ce qui
  faisait échouer une sphère légèrement enfoncée, posée en premier (§ 8, point 1).
- La base du plan tangent est prise sur S₁ quand S₂ dégénère.
- Les quatre équations de Newton n'y sont pas définies (la normale par les
  paramètres) : partant d'un point voisin, les itérés échouent au pôle. Le
  contact est alors **le pôle lui-même** (`TangentMarcher::at_pole`), quand il est
  dans la portée du départ : une borne de la boîte le long de laquelle S_u (ou
  S_v) s'annule, le pied du pôle sur l'autre surface à la tolérance, sa normale
  parallèle à la normale limite du pôle (S_uv × S_v, puisque S_u(u, v + h) ≈
  h S_uv), le tout dans les boîtes et II₁ − II₂ définie. Le point rendu est exact
  (le pôle et son pied), là où la première version gardait le départ tel quel,
  sur le seul critère d'un écart dans la tolérance.

### 3.3 Contacts dans la tolérance

Un modèle lu d'un fichier pose rarement la sphère exactement sur la boîte : elle
peut être à une fraction de la tolérance au-dessus, ou enfoncée d'autant. Entre
primitives, les formes fermées donnent alors le point de contact tant que l'écart
ne dépasse pas la tolérance, et une boucle au-delà. Le cheminement fait
maintenant de même.

- **Au-dessus, dans la tolérance** : les surfaces ne se coupent pas, la recherche
  des départs (S₁ = S₂ à 10⁻² tol près) n'en trouve aucun, et le contact était
  perdu. Les paires presque parallèles, descendues au seuil de subdivision sans
  départ, essaient maintenant `isolated_contact` depuis leurs centres : un contact
  certifié (écart dans la tolérance au point corrigé, forme définie) est gardé
  (`touches`). Ce n'est fait qu'une fois par voisinage.
- **Enfoncées dans la tolérance** : l'intersection vraie est une petite boucle de
  rayon √(2 profondeur / |l|), plus petite que la portée du contact. Soit les
  départs y sont arrêtés (normales presque parallèles), soit la boucle est tracée.
  Une boucle tracée dont tous les points sont dans la portée d'un contact certifié
  (le Newton de la tangence part de l'un d'eux) est ce contact. Elle reste une
  branche (ses départs sont écartés), mais n'est plus une courbe. Autour d'un pôle,
  elle s'arrête à la couture : elle n'est pas fermée, et c'est le critère de
  portée qui compte.
- **Départs non certifiés** : un départ arrêté qui ne se certifie pas, mais qui est
  dans la portée d'un contact certifié par un autre départ, est ce contact ; il ne
  rend plus l'intersection singulière.
- Un contact dont la portée contient une courbe tracée est celui de la courbe.

### 3.4 Cheminement

Dans la boucle finale sur les contacts, un départ arrêté (`stopped`) et certifié
devient un point de contact non singulier, à sa position corrigée ; s'y ajoutent
les contacts trouvés sans départ et ceux des boucles (§ 3.3). Sans dérivées
secondes (une surface décalée, par exemple), rien ne change : il reste
singulier.

## 4. Résultats

| | Avant | Après |
|---|---|---|
| `nurbs_sphere_touching` (galerie) | `NonGeneric` sur fuse, cut, common | 7, 6 et 0 faces, comme la sphère exacte ; valide, maillage fermé |
| OCCT : `sphere_touching_nurbs` (boîte et sphère NURBS) | — | mêmes volumes, faces et solides qu'OCCT |
| OCCT : `plane_touching_nurbs_sphere` (boîte exacte, sphère NURBS) | — | idem |
| précision du point de contact | 4·10⁻⁸ | 10⁻¹⁵ |
| sphère NURBS 2·10⁻⁷ à 9·10⁻⁷ au-dessus du plan (tol 10⁻⁶) | aucun point | un point, comme les formes fermées |
| sphère NURBS enfoncée de 10⁻¹² à 5·10⁻⁷, sur son pôle | singulier (1 à 8 points) | un point, non singulier, dans les deux ordres |
| sphère NURBS enfoncée de 6·10⁻⁷ à 9,9·10⁻⁷ | une boucle de rayon 10⁻³ environ | un point, comme les formes fermées ; au-delà de la tolérance, une boucle, comme elles |

## 5. Tests

- `tests_intersect_marching.isolated_contacts` (nouveau) :
  - sphère NURBS sur un plan et sur un carreau plat B-spline, sphère NURBS dans une
    cuvette, cuvette dans une cuvette plus large ;
  - axes de la sphère (1, 1, 1), x et z : contact hors pôle et au pôle ;
  - les deux ordres des surfaces ; un point, aucune courbe, ni singulier ni
    incomplet, à 10⁻⁹ près sur les deux surfaces ;
  - placements dans la tolérance mais pas exacts : sphère au-dessus du plan de
    2·10⁻⁷ et 9·10⁻⁷, enfoncée de 10⁻¹² à 9,9·10⁻⁷, sur son pôle et hors pôle, dans
    les deux ordres. Un point non singulier, au milieu de la sphère et du plan à
    10⁻⁹ près, autant de points et de courbes que les formes fermées. Au-delà de la
    tolérance (±1,5·10⁻⁶), une boucle (enfoncée) ou rien (au-dessus), comme elles.
    Ces cas passent par le pôle lui-même (`at_pole`), l'absorption des départs non
    certifiés, les contacts sans départ et les boucles dans la portée d'un contact.
- `free_surfaces` : la sphère sur un plan par cheminement, à 10⁻⁹ près au lieu de
  10⁻³, non singulière.
- `tangency_limits` : commentaire mis à jour (le contact isolé n'est plus singulier).
- Les cas singuliers existants restent singuliers : surface sans dérivées secondes,
  tangence qui finit en se croisant (`tangency_limits`).
- Comparaison OCCT : deux paires ajoutées (§ 4).
- Python : `nurbs_sphere_touching` n'est plus une limite connue ; nouveau
  `test_nurbs_contact_like_primitives` (mêmes nombres de solides et de faces que la
  sphère exacte, pour les trois opérations). Les booléens de chaque cas de la
  galerie sont calculés une fois pour les tests qui les regardent.

Les 737 tests de `ctest` passent, comparaison à OCCT comprise ; `pytest` : 95 tests
passés sans les notebooks, sans `xfail`, avec et sans `CI`. Les booléens et la
comparaison à OCCT prennent le même temps qu'avant ces corrections.

## 6. Limites

| Limite | Suite |
|---|---|
| Un contact dégénéré (II₁ − II₂ nulle : z = x⁴ + y⁴ sur son plan tangent) reste singulier | ordre supérieur de l'écart ; rare en conception |
| Un contact presque dégénéré (\|l₁\| ≤ `angular_tol` \|l₂\|) passe par le tracé d'une courbe de tangence, comme avant | — |
| Le partage entre contact isolé et courbe de tangence se fait sur le rapport \|l₁\| / \|l₂\|, sans la tolérance : une bille de rayon 1 dans un alésage avec un jeu de 0,002 est un point, alors que les surfaces restent à la tolérance l'une de l'autre sur ±0,032 le long de l'alésage ; avec un jeu de 5·10⁻⁴, c'est une courbe. Les formes fermées font de même | à trancher (§ 7, point 1) : comparer plutôt la portée √(2 tol / \|l₁\|) au pas ou à `small`, entre primitives aussi |
| Sans dérivées secondes (surface décalée), le contact reste singulier | dérivées secondes exactes de `OffsetSurface` (PR 3 du palier) |

## 7. Points à valider

1. **Certification par II₁ − II₂ définie**, avec le seuil `angular_tol` de
   `tangent_branch` pour séparer contact isolé et courbe de tangence. Ce seuil ne
   dépend pas de la tolérance (§ 6) ; il est gardé ici pour rester d'accord avec
   les formes fermées.
2. **Au pôle, formes prises à 10⁻⁴ de la distance au milieu de la boîte**, chacune
   sur la normale propre de sa surface.
3. **Au pôle, le contact est le pôle lui-même**, avec son pied sur l'autre surface,
   quand Newton ne peut pas y être appliqué.
4. **Contacts dans la tolérance** comme entre primitives : un point quand l'écart
   (au-dessus ou enfoncé) ne dépasse pas la tolérance, une boucle au-delà.

## 8. Revue de code

Huit remarques, sept traitées dans cette PR ; la sixième est une question de
conception, laissée au § 7.

| # | Remarque | Suite |
|---|---|---|
| 1 | Au pôle, le contact n'est certifié que si le placement est exact : une sphère NURBS enfoncée de 10⁻¹² dans le plan, posée en premier, reste singulière ; de 10⁻⁹ à 5·10⁻⁷, 7 à 8 points singuliers ; plan en premier, échec dès 10⁻⁷. | Formes prises sur la normale propre de chaque surface à travers les pôles ; contact au pôle lui-même (`at_pole`) quand Newton échoue (§ 3.2). Tests de placements dans la tolérance. |
| 2 | Un départ arrêté non certifié rend singulier un contact certifié par un autre départ (plan en premier, enfoncée de 10⁻⁷ : un point certifié, un départ à 4,5·10⁻⁴, singulier). | Départ absorbé quand il est dans la portée d'un contact certifié (§ 3.3). |
| 3 | Un contact dans la tolérance mais pas exact (sphère 2·10⁻⁷ à 9·10⁻⁷ au-dessus du plan) n'a aucun départ et se perd ; les formes fermées le donnent. (La remarque disait aussi « enfoncée de 9·10⁻⁷ : rien » : c'était une boucle, que la sonde n'affichait pas.) | Contacts cherchés depuis les paires presque parallèles sans départ ; boucle dans la portée d'un contact remplacée par lui (§ 3.3). Parité mesurée avec les formes fermées de −1,5·10⁻⁶ à +1,5·10⁻⁶. |
| 4 | Les tests ne couvrent que des placements exacts ; le repli au pôle n'est pas couvert. | Placements dans la tolérance, au pôle et hors pôle, dans les deux ordres, comparés aux formes fermées (§ 5) ; ils passent par `at_pole`. |
| 5 | Le repli au pôle gardait le départ sur le seul critère d'un écart dans la tolérance, sans la précision annoncée. | Remplacé par le pôle et son pied, exacts (§ 3.2). |
| 6 | Le partage point ou courbe se fait sur un rapport de courbures, sans la tolérance. | Laissé : il est celui des formes fermées ; décrit au § 6 et posé au § 7, point 1. |
| 7 | La doc du déplacement au pôle (« 10⁻⁴ de l'étendue ») ne dit pas ce que fait le code. | Doc et note : 10⁻⁴ de la distance au milieu de la boîte de l'intersection. |
| 8 | Le test Python refait les six booléens de deux cas déjà calculés. | Booléens d'un cas calculés une fois (`functools.cache`), partagés par les tests. |
