# Palier 2 bis — PR 15 : un rayon de confusion pour les faces

| | |
|---|---|
| PR | branche `fix/brep-face-confusion-radius` |
| Document de référence | [brep2_pr06_section.md](brep2_pr06_section.md) § 8 et [brep2_pr09_boolean.md](brep2_pr09_boolean.md) § 7 (revue du noyau BREP, lots 4 et 5), [brep_palier2bis.md](brep_palier2bis.md) § 12 |
| Fichiers | `gbs/gbsconstants.h`, `gbs/intersect/surface_surface.h`, `gbs-brep/section.h`, `gbs-brep/coincident.h`, `gbs-brep/partial.h`, `gbs-brep/boolean.h`, `tests/tests_brep_boolean.cpp`, `tests/tests_intersect_surfaces.cpp`, `docs/sources/design/brep2_pr06_section.md`, `docs/sources/design/brep2_pr09_boolean.md`, `docs/sources/design/brep_palier2bis.md`, `News.md` |
| Statut | Validé le 29 octobre 2026 |

## 1. Objectif

Les lots 4 et 5 de la revue ont laissé une limite : deux faces à quelques
tolérances l'une de l'autre, plus loin que `tol` mais dans le rayon de confusion
des sommets (10 tol), font échouer les booléens. Cas de la revue, `tol` = 1e-6 :
`box.stp` ([0, 2] × [0, 3] × [0, 4]) et sa copie déplacée de (−1 ; 0,5 ; −d),
dont les dessus, et les dessous, sont à d l'un de l'autre.

| d | fusion, commun, différence (master) |
|---|---|
| 0 et 0,5 tol | justes (commun 2,5 (4 − d)) |
| 1 à 2 tol | « a section edge does not split a face » |
| 2,5 à 10 tol | « split_face: the cuts and the boundary do not form closed loops » |
| 10,5 tol et plus | justes : copeaux de largeur d entre les faces |

## 2. Cause

Trois distances de confusion différentes :

- les **faces** sont confondues à `tol` près : formes closes de
  `gbs/intersect/surface_surface.h` (plans parallèles à `tol` sur la région
  commune, écart à `tol` au milieu), points projetés de
  `section_detail::coincident_faces` / `on_surface` (B-splines) ;
- les **sommets** le sont à 10 tol (`section_detail::VertexPool`, les bouts des
  courbes d'intersection interpolées) ;
- l'**empreinte** des faces confondues (`coincident.h`) croise les arêtes à
  `tol`.

Entre `tol` et 10 tol, les dessus ne sont pas confondus, mais la section d'un
flanc d'une boîte avec le dessus de l'autre longe une **partie** d'une arête du
dessus, à d. `edges_along` ne reconnaît que les morceaux d'un sommet à un
autre : jusqu'à 2 tol, le morceau est classé sur le bord (`On`) et abandonné
comme un contact, et les dessus ne sont découpés par rien (« a section edge
does not split a face ») ; au-delà, il devient une coupe à d de l'arête, dont
un bout se confond avec un sommet de l'arête et l'autre non (« do not form
closed loops »).

## 3. Mesures sur master

Sondes : fusion, commun et différence, `check()` et volume à 2e-3 près, pour d
de 0 à 1 000 tol.

| Configuration | Juste sur master | Refusée ou fausse sur master |
|---|---|---|
| `box.stp` et (−1 ; 0,5 ; ±d), (1 ; 1 ; d) : dessus et dessous à d | d ≤ 0,5 tol, d ≥ 10,5 tol (copeaux) | 1 à 10 tol |
| empilées, recouvrement d ((0,5 ; 0,5 ; 4 − d)) | d ≤ 0,5 tol, d ≥ 10,5 tol | 1 à 10 tol (coque ouverte, puis coupes non fermées) |
| empilées, écart d ((0,5 ; 0,5 ; 4 + d)) | toujours (deux solides au-delà de `tol`) | — |
| la même boîte, déplacée de (0 ; 0 ; d) | d ≤ 0,5 tol, d ≥ 2,5 tol (copeaux de d) | 1 à 2 tol |
| cubes à faces B-spline bilinéaires 0..2 et (1 ; 1 ; d) | d ≤ 0,5 tol, d ≥ 10,5 tol | 1,5 à 10 tol |
| cubes B-spline empilés, recouvrement d | d ≤ 0,5 tol, d ≥ 10 tol | 1,5 tol (deux solides), 2,5 à 9,5 tol |
| `box.stp` et sa copie tournée de d rad autour de y (dessus à ± d/2 aux bouts du recouvrement) | d = 0, d ≥ 30 tol | 0,5 à 20 tol |
| cylindres R 5 (primitives), axes à d | d ≤ 0,5 tol | toutes les autres valeurs |
| cylindres à flancs B-spline rationnels, axes à d | d ≤ 0,5 tol | toutes les autres valeurs |

## 4. Choix : un seul rayon de confusion

| Option | Effet | Retenue |
|---|---|---|
| faces confondues au rayon des sommets, R = 10 tol | les faces à moins de R sont la même face, comme leurs sommets : empreinte, puis une seule copie des morceaux confondus | **oui** |
| sommets confondus à 2 tol | les bouts des courbes d'intersection interpolées, à quelques tol, ne se rejoindraient plus ; les copeaux de 2 à 10 tol resteraient dans les résultats | non |
| classer toutes les faces contre l'autre opérande à R | un copeau plus large que R mais plus étroit que 2R n'a plus de point intérieur : empilées à 10,5–20 tol de recouvrement, « a face without inside », juste sur master (mesuré) | non |
| rayon en option (le « fuzzy value » d'OCCT) | hors de portée : R suit `tol`, comme pour les sommets | non |

R est une constante nommée, `brep_confusion_factor` (10) dans
`gbs/gbsconstants.h`, qui remplace les `T(10) * tol` littéraux (sommets,
`edges_along`, `same_side`, `partial.h`). Les résultats sont **tolérants** : une
arête imprimée sur une face à d d'elle porte une tolérance d (≤ R), comme la
pcurve qui l'y pose.

## 5. Réalisation

### 5.1 Décision de la coïncidence (`section.h`, `surface_surface.h`)

| Paire | Avant | Après |
|---|---|---|
| deux plans | parallèles à `tol` sur la région commune (borne 2 · ‖h‖ · angle, toute la diagonale de la région), écart à `tol` au milieu | **plus grand écart** sur la région commune ≤ R : aux coins de la région projetés sur le second plan, l'écart y étant affine (`coincident`) |
| deux quadriques | formes closes à `tol` | formes closes à R (coaxiales à R, même méridienne à R), **confirmées** par les points des faces (`coincident_faces` à R) : chaque test admet R, l'écart pourrait atteindre 2R |
| autres paires | points de chaque face projetés sur l'autre surface, à `tol` ; bords posés (`border_on`) à `tol` | à R |
| section par une surface | formes closes à `tol` | primitives : `coincident` à R ; autres : inchangé |

`coincident(s1, box1, s2, box2, tol)` est une fonction nouvelle de
`surface_surface.h` ; la région commune des boîtes est extraite de `intersect`
(`intersect_detail::common_region`). Les courbes d'intersection restent
calculées à `tol`.

### 5.2 Morceaux le long d'une partie d'arête (`section.h`, `edges_along`)

Un morceau de section **intérieur** à une face (classé `In` à 2 tol), dont les
bouts et 5 points intérieurs sont à moins de R (ou de la tolérance de l'arête)
d'**une** arête du bord, est cette arête : rapportée dans `boundary_edges` et
`edges_on_faces`, imprimée par le booléen sur l'autre face (`imprint_edge`),
sauf si elle est confondue avec une face de la sienne (l'empreinte des faces
confondues l'a déjà fait). Le cas d'un sommet à un autre (lot 4) est inchangé.

Un morceau **sur le bord** en son milieu (`On`) reste un contact : le cercle où
le flanc y = 0 de la boîte coupe le cylindre couché le long de son seam
(`edges_on_faces`) touche l'arête du dessus sur 2e-4, à 1e-9 près ; le
rapporter comme cette arête faisait échouer la fusion (coque ouverte).

### 5.3 Empreinte (`coincident.h`)

À R : croisements des arêtes (courbe-courbe), sommets d'arrivée (`end_at`),
fusion des arêtes superposées (milieux à R), acceptation d'une pcurve (l'arête
posée sur la surface à R : plan, primitive, projection). La tolérance de l'arête
est relevée à l'écart de sa pcurve, comme avant (`cover_deviation`).

### 5.4 Coïncidences partielles (`partial.h`)

Un carreau est posé sur l'autre surface quand ses échantillons sont à R ; la
boîte d'un carreau croisant est écartée de la frontière là où l'écart dépasse R
(10 tol littéral). Sans cela, un carreau à quelques tol serait « hors », et le
morceau qui le porte serait intersecté avec une surface confondue avec lui.

### 5.5 Sélection (`boolean.h`)

Un morceau d'une face confondue est posé sur l'autre opérande au-delà de `tol` :
il est cherché sur les **morceaux des faces avec lesquelles sa face est
confondue**, à 2R (la marge des formes closes), la plus proche (`same_side`),
avant la classification. Les autres faces restent classées à `tol`. Sans la
restriction aux partenaires, le copeau de largeur d d'une face confondue (la
même boîte déplacée de (0 ; 0 ; d), d > R) trouvait la face voisine à d/2, de
normale perpendiculaire, et la fusion perdait le copeau (coque ouverte à
10,5–30 tol, mesuré).

### 5.6 Écart des pcurves (`surface_surface.h`, `fit_primitive_pcurve`)

Une arête posée à d d'une face porte une pcurve dont l'écart, au plus fort, vaut
d ; `check()` le compare à la tolérance de l'arête, relevée à l'écart mesuré.
Deux sous-estimations sont apparues :

| Pcurve | Avant | Après |
|---|---|---|
| exacte (segment) | écart mesuré aux milieux de 1 000 intervalles : un écart qui croît le long de la courbe (plan tourné) était sous-estimé de 0,05 % au bout | 1 001 points, **bouts compris** |
| interpolée | écart mesuré aux milieux et aux points de contrôle, pas aux nœuds d'interpolation, où la pcurve passe par le pied : l'écart de la courbe à la surface y était perdu (cercle sur un cylindre d'axe à 1,5 tol : 1,499993e-6 au lieu de 1,5e-6 en t = 0) | aussi aux **nœuds** (bouts compris) |

## 6. Effet

| Configuration | Après |
|---|---|
| `box.stp` et (−1 ; 0,5 ; ±d), (1 ; 1 ; d) | justes de 0 à 1 000 tol ; jusqu'à R, la topologie de d = 0 (fusion 14 faces, commun 6, différence 8) |
| empilées, recouvrement d | justes, sauf à d = R exactement (§ 8, point 4) |
| empilées, écart d | jusqu'à R, **un** solide (les faces se touchent, comme à d = 0) ; au-delà, deux |
| la même boîte, déplacée de (0 ; 0 ; d) | jusqu'à R, la même boîte (fusion et commun de 6 faces, différence vide) ; au-delà, copeaux de d |
| cubes B-spline (1 ; 1 ; d) ; empilés | justes, sauf à d = R exactement pour les premiers |
| boîte tournée de d rad | justes de 0 à 100 tol, sauf à 20 tol exactement (plus grand écart égal à R) ; jusqu'à 19,5 tol, la topologie de d = 0 |
| cylindres primitifs, axes à d | justes jusqu'à R ; au-delà, refusés comme sur master |
| cylindres B-spline, axes à d | justes jusqu'à 9,5 tol ; au-delà, refusés comme sur master |

Coût : les 16 cas existants de `tests_brep_boolean` prennent 25,6 à 25,8 s,
sur master comme sur la branche (trois passes alternées). Le test de
coïncidence par les formes closes refait la région commune des boîtes (paires
de primitives) ; les quadriques confondues sont confirmées par 72 points au
plus.

## 7. Tests

| Test | Cas | master |
|---|---|---|
| `faces_a_few_tolerances_apart` (`tests_brep_boolean`) | `box.stp` et (−1 ; 0,5 ; −d) à 1,5, 2,5 et 5 tol : volumes, `check()`, 14 faces à la fusion et 6 au commun ; empilées, recouvrement d : un solide de 12 faces, commun vide, différence 24 | « a section edge does not split a face » (1,5 tol), « do not form closed loops » |
| `coincident_within_the_confusion_radius` (`tests_brep_boolean`) | à 5 tol : cubes B-spline, boîte tournée, cylindres primitifs d'axes à d, cylindres B-spline d'axes à d ; volumes et `check()` | « do not form closed loops » |
| `coincident_within_a_distance` (`tests_intersect_surfaces`) | `coincident` : plans à 5 tol (oui à 10 tol, non à 4 tol, `intersect` non) ; plan tourné, écart de 9 tol (oui) et 11 tol (non) ; cylindres coaxiaux de rayons à 5 tol (oui) et 15 tol (non), axes à 5 tol (oui) | fonction nouvelle |
| `pcurve_deviation_at_the_ends` (`tests_intersect_surfaces`) | l'écart d'une pcurve couvre celui des bouts : segment sur un plan tourné (exacte), cercle sur un cylindre d'axe à 5 tol (interpolée) | 4,9975e-6 et 4,99998e-6 pour 5e-6 |

Toutes les suites `tests_brep_*`, `tests_step_*` et `tests_intersect_*`
passent.

## 8. Points à valider

1. **R = 10 tol pour les faces**, le rayon des sommets, facteur fixe
   (`brep_confusion_factor`), pas une option.
2. **Décider à R, retrouver à 2R** : la coïncidence est décidée à R ; la face
   qui porte un morceau confondu est cherchée à 2R, parmi les morceaux des faces
   partenaires seulement.
3. **Changement de comportement** : deux solides dont des faces sont à moins de
   R se touchent. Empilées avec un écart de 1 à 10 tol, la fusion rend un solide
   (deux sur master), comme à écart nul.
4. **Limites**, hors de portée :
   - à d = R exactement, les tests au rayon (coïncidence, morceau le long d'une
     arête, croisements, pcurves) calculent la même distance avec des arrondis
     différents et peuvent se contredire (fenêtre de l'ordre de 1e-16) : refus
     pour les boîtes empilées et les cubes B-spline à 10 tol, la boîte tournée
     de 20 tol. Une marge relative en aval (1e-6) n'en corrigeait qu'une
     partie ; elle n'a pas été gardée ;
   - deux cylindres d'axes à plus de R (10,5 à 100 tol), presque tangents le
     long d'une génératrice : refusés comme sur master.
