# Palier 4 — PR 7 : coins de congés, pied d'aube, bilan

| | |
|---|---|
| PR | branche `feat/brep-fillet-corners` |
| Document de référence | [brep_palier4.md](brep_palier4.md), § 4.6, § 6 et § 10 (bilan) |
| Fichiers | `gbs-brep/blend.h`, `gbs-brep/fillet.h`, `python/gbsBindCurveTools.h`, `testing/blade_fixture.h` (nouveau), `tests/tests_brep_fillet.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/examples/gbs_examples.py`, `docs/sources/examples/gallery.py`, `docs/sources/design/brep_palier4.md`, `docs/sources/design/brep4_pr06_fillet_topology.md`, `News.md`, cette note |
| Statut | Validé (§ 9, 18 nov. 2026) ; revue de code : § 10 |

## 1. Objectif

Dernière PR du palier 4 :

- les **coins** : trois arêtes congées à un même sommet (coin de boîte, boîte
  aux douze arêtes congées, coin intérieur d'une poche), la sphère de coin du
  § 4.6 du plan ;
- un exemple de congé plus exigeant, demandé pour la galerie : une **aube**
  (le profil de `tests_foils`, `type2_blade`) fusionnée à son moyeu, congée en
  pied ;
- la comparaison à OCCT de ces cas, et le **bilan du palier** (§ 10 du plan).

Le pied d'aube a révélé trois manques de la géométrie des congés (PR 5),
corrigés ici (§ 4.2 à § 4.4).

## 2. Sources

- Rossignac et Requicha, « Constant-radius blending in solid modelling »
  (*Computers in Mechanical Engineering*, 1984) : le coin de trois congés de
  même rayon, sphère centrée à r des trois faces.
- Vida, Martin et Várady, « A survey of blending methods that use parametric
  surfaces » (*Computer-Aided Design* 26(5), 1994), § 5 : les raccords de
  sommet (*vertex blends*), dont la sphère est le cas le plus simple.
- Le comportement de CATIA et SolidWorks : trois congés de même rayon en un
  coin se raccordent par une sphère, sans option. OCCT
  (`BRepFilletAPI_MakeFillet`) sert de comparaison, en boîte noire.

## 3. API

Pas de nouvelle fonction : `fillet` traite les coins qu'il refusait. Python :
`gbs.c3_connect(crv1, crv2)` (raccord C3 de deux courbes, ce qui sert aux bords
d'attaque et de fuite du profil) ; dans les exemples, `blade_foil()`,
`blade(m)`, `blade_on_hub(m)`.

## 4. Démarche

### 4.1 Coin sphérique

À un sommet de trois arêtes, toutes congées (`fillet`, un seul rayon), les
trois congés se raccordent par une sphère :

1. **Centre.** Les trois lignes des centres passent par le point à r des trois
   faces : celle des faces A et B est l'intersection de leurs décalés, celle de
   B et C aussi, et le point commun aux trois décalés est sur les deux.
   Gauss-Newton sur les trois paramètres (six équations : ligne 1 = ligne 2,
   ligne 2 = ligne 3, moindres carrés), à partir des paramètres aux sommets.
   Refus si le résidu, ou l'écart d'une ligne au centre, dépasse
   10 tol (+ `approx_tol` par congé approché), ou si un paramètre sort du
   congé prolongé.
2. **Pieds.** Chaque face est bordée par deux des trois congés, dont les
   contacts au centre doivent coïncider : c'est le pied du centre sur la face.
   Un sommet par face, partagé par les deux contacts. Refus si une face n'est
   pas bordée par exactement deux congés, ou si les deux contacts s'écartent.
3. **Arcs.** Chaque congé s'arrête sur son arc au paramètre du centre : arc de
   centre c et de rayon r entre ses deux contacts, donc grand cercle de la
   sphère. L'arête porte la pcourbe de l'arc sur le congé (segment en travers)
   et sur la sphère : n + 1 points projetés, interpolés en cubique, n doublé à
   partir de 32 tant que l'interpolation s'écarte de plus de tol / 10 entre les
   points (1024 au plus). La tolérance de l'arête est relevée à l'écart mesuré
   aux points (`settle`) et entre eux.
4. **Sphère** (`corner_frame`). Le pôle z est, parmi 2048 directions
   réparties sur la sphère, la plus éloignée des arcs, ni elle ni son opposée
   n'étant dans le triangle sphérique du coin. La couture (u = ±π) est au
   milieu du plus grand intervalle de longitudes que laissent les arcs autour
   de z, et x lui est opposé. La face est bornée par les trois arcs chaînés,
   sens trigonométrique en (u, v), et orientée selon la normale de la sphère
   au milieu d'un arc : le solide est à l'intérieur de la sphère à un coin
   convexe, à l'extérieur à un coin concave.

Convexe et concave se traitent de la même façon. Au coin d'une boîte, le
solide est à l'intérieur de la sphère ; au coin d'une poche, il est à
l'extérieur. À un coin pointu (sommet d'un tétraèdre dont les faces
s'écartent de 150° autour de lui), le triangle du coin remplit presque une
demi-sphère : tout pôle est à quelques degrés des arcs (4° au plus pointu
essayé), et c'est l'interpolation plus fine des pcourbes qui garde les arcs
dans la tolérance de leurs arêtes. Les faces des
arêtes (§ 4.4 de la note PR 6) n'ont rien de nouveau : les contacts de deux
congés se rejoignent au pied, comme ceux d'une chaîne tangente, et la suite
de co-arêtes qui passe par le coin est remplacée par les deux contacts.

Le plan (§ 4.6) prévoyait des faces deux à deux orthogonales en première
version. Rien dans la construction ne le demande : le prisme triangulaire,
faces à 60°, est traité et comparé à OCCT.

### 4.2 Face fermée non périodique

La surface de l'aube est une B-spline fermée en u (les sections sont des
courbes fermées), sans être périodique. Le congé du pied en fait le tour :
la ligne des centres cheminée la traverse, et les stations sortaient du
domaine de la B-spline, qui lève une exception au-delà de ses bornes.

- **Vue sans couture** (`fillet_detail::SeamlessSurface`). Le support d'une
  face fermée (`surface_closure`) qui ne se laisse pas évaluer au-delà de sa
  couture (`unrepeated` : une B-spline, un décalé de B-spline) est enveloppé
  d'une vue qui ramène u (ou v) dans les bornes par la largeur du domaine.
  Les décalés, la ligne des centres, les stations et les contacts la voient
  continue. Les périodes (`view_periods`, `periods_of`) s'en déduisent comme
  pour les primitives. Une surface qui se répète déjà au-delà de ses bornes
  (surface de révolution, décalé d'une primitive) n'est pas enveloppée : le
  décalé d'un cylindre reste un cylindre exact pour `offset_surface`.
- **Départ sur la couture.** Pour une arête fermée, la boucle du congé part
  là où le contact sur cette face croise la couture (`crossings`). La
  pcourbe du contact court alors d'une borne à l'autre ; décalée d'un nombre
  entier de périodes, ses pôles ramenés sur les bornes à tol près dans
  l'espace (`into_bounds`), elle reste dans le domaine de la face. Le sommet
  du contact est sur l'arête de couture de la face, qui est coupée là, comme
  pour un bord de cylindre (PR 6).
- **Deux faces fermées.** Le cheminement s'arrête sur les coutures : autour
  de deux faces fermées, la ligne des centres ne se ferme que si leurs
  coutures se rencontrent dans une même section du congé (dans un même plan
  méridien, pour un tronc de cône sous un cylindre). Des coutures écartées le
  long de l'arête la coupent en deux arêtes, congées en chaîne.

### 4.3 Raccords C2 des sections

Le profil est formé de quatre morceaux (deux côtés, bords d'attaque et de
fuite) joints en une B-spline de degré 7 dont les nœuds de jonction ont la
multiplicité 7. Elle n'y est que C2, et le congé qui la longe n'y est que C1.
Deux défauts de l'approximation de la PR 5 en venaient :

- La ligne des centres et les contacts étaient interpolés en **cubique
  global C2** sur les stations : près d'un nœud, l'erreur oscillait et ne
  diminuait plus quand on doublait les stations (3·10⁻⁵ à 512 comme à 2048).
  Ils sont maintenant interpolés en **Hermite cubique** (C1, local), comme la
  surface : valeurs et dérivées le long de la ligne des centres, par
  différences centrées, aux stations.
- Une station est placée là où un contact croise un **nœud de sa face moins
  que C3** (`weak_knots` : multiplicité au moins le degré moins 2), trouvé
  sur la pcourbe cheminée puis affiné par bisection sur les stations
  elles-mêmes.

### 4.4 Raffinement local

Les intervalles étaient tous doublés jusqu'à ce que les points aux quarts de
chacun passent, avec 2048 intervalles au plus. Avec l'Hermite, un intervalle ne
dépend que de ses deux stations : seuls les intervalles hors tolérance sont
coupés en deux. Chacun est contrôlé en son milieu, là où l'Hermite s'écarte le
plus, et en 1/2 ∓ 1/(2√3), là où sa pente s'écarte le plus (les normales). Le
seuil est la moitié des tolérances, pour garder une marge aux points entre les
contrôles : deux tests de tangence de la PR 5 passaient de peu au-dessus de
`angular_tol` sans elle. Le contact du pied balaie le petit bord d'attaque très
vite : c'est là que les intervalles se resserrent.

### 4.5 Tolérances des arêtes

Après la fusion, la face de l'aube a une tolérance plus large que `tol`. Les
contacts, créés à `tol`, violaient l'ordre arête ≥ face de `check`. Avant le
contrôle, chaque arête du résultat est relevée à la tolérance des faces
qu'elle borde, et ses sommets à la sienne.

## 5. Refus

En plus de ceux de la PR 6 (`NonGeneric`) :

| Situation | Message |
|---|---|
| chanfreins à un coin, ou trois arêtes congées à un sommet de plus de trois arêtes | `more than two blended edges at a vertex, or three that are not fillets at a vertex of three edges` |
| coin mêlant arêtes convexes et concaves (haut du coin intérieur d'un L) | `a corner of convex and concave blended edges (not handled)` |
| lignes des centres qui ne se rencontrent pas dans les congés prolongés (centre au-delà d'un sommet de plus que `extension`) | `fillets at a corner whose centre lines do not meet within the blends continued` |
| contacts de deux congés qui ne se rejoignent pas sur une face | `fillets at a corner whose contacts do not meet on a face` |
| aucun pôle de la sphère à plus de 3° des arcs | `a corner whose sphere has no pole clear of its arcs` |
| arête fermée entre deux faces B-splines fermées dont les coutures ne se rencontrent pas dans une même section du congé | `a closed edge between two closed B-spline faces whose seams do not meet in the same section of the fillet (not handled)` |
| faces du coin pas bordées deux à deux | `fillets at a corner whose faces are not shared two by two` |
| contact traversant la couture d'une face non périodique ailleurs qu'en boucle | `a contact crossing the seam of a face that is not periodic (a closed B-spline)` |

Les messages qui renvoyaient à la PR 7 disent maintenant « not handled ».

## 6. Python, galerie

`gbs.c3_connect` (courbes non rationnelles de dimension 1 à 3). Galerie,
section « Fillets and chamfers », trois cas de plus :

| Cas | Faces | Temps |
|---|---|---|
| coin de boîte, ses trois arêtes (r = 0,25) | 10 | 0,01 s |
| boîte arrondie, ses douze arêtes (r = 0,2) | 26 | 0,01 s |
| aube fusionnée à son moyeu, congé en pied (r = 0,03) | 6 | 0,16 s |

Le cas d'aube a une vue rapprochée du bord d'attaque (`FilletCase.camera`) :
le congé est petit devant l'aube.

## 7. Tests

`tests/tests_brep_fillet.cpp`, 8 cas de plus (31 en tout) :

| Test | Contenu |
|---|---|
| `fillet_box_corner` | coin de boîte (r = 0,1 et 0,3) : 1 − 3(1 − π/4) r²(1 − r) − r³(1 − π/6) à 10⁻⁹ ; une seule face sphérique, de rayon r, centrée en (1 − r, 1 − r, 1 − r) |
| `fillet_box_all_edges` | boîte 2 × 1,5 × 1, ses douze arêtes : 26 faces, volume de la boîte arrondie (noyau, dalles, quarts de cylindres, huit huitièmes de boule) à 10⁻⁸ |
| `fillet_concave_corner` | coin intérieur d'une poche cubique dans un bloc : 7 + 3(1 − π/4) r²(1 − r) + r³(1 − π/6) à 10⁻⁹ |
| `corner_refusals` | chanfreins à un coin, modèle inchangé ; haut du coin intérieur d'un L (une arête concave, deux convexes) |
| `fillet_offset_of_cylinder` | cylindre dont le côté est le décalé d'un cylindre (bornes d'un tour : il se ferme, mais se répète au-delà de sa couture) : congé du bord supérieur exact, sur un tore ; π − 2π(1 − x)(1 − π/4) r² à 10⁻⁶ (x : centroïde du coin) |
| `fillet_two_closed_bsplines` | tronc de cône sous un cylindre, côtés B-splines rationnelles fermées : arête concave congée autour des deux faces, coutures dans un même plan (arête fermée, volume à 10⁻⁶ de l'estimation par Pappus) ou écartées (deux arêtes en chaîne, même volume à 10⁻⁸) ; couture du tronc oblique : refus |
| `fillet_sharp_corner` | sommet d'un tétraèdre aux faces écartées de 150°, ses trois arêtes congées : un pied à plus de 90° de la somme des pieds ; arcs de la sphère dans les bornes de ses paramètres et à moins de la tolérance de leur arête sur 1025 points (33 points interpolés laissaient 3·10⁻⁵ pour une arête à 10⁻⁶) |
| `fillet_blade_root` | aube sur moyeu (`testing/blade_fixture.h`) : arête de pied fermée, géométrie sous `approx_tol` et tangente aux faces (`expect_blend`), solide valide à 6 faces, volume ajouté à 10 % de (1 − π/4) r² × longueur du pied (le moyeu se dérobe sous le congé : 5 % de moins), solide donné inchangé |

Les suites existantes passent sans changement, sauf deux tests de tangence de
la PR 5 que l'Hermite local faisait passer de peu au-dessus de `angular_tol`
(1,0 et 1,4·10⁻⁵ pour 10⁻⁵) avant la marge du § 4.4.

`gbs-occt/tests/tests_brep2_occt.cpp`, `fillet_corners_and_blade`, contre
`BRepFilletAPI_MakeFillet` sur le même solide (`expect_blended` : résultat de
gbs relu par OCCT valide, même nombre de faces, même volume à 10⁻⁶ près) :

- coin de boîte ; boîte aux douze arêtes ;
- coin intérieur de poche (trois arêtes concaves) ;
- coin du dessus d'un prisme triangulaire (faces à 60°).

Pour le pied d'aube, le résultat d'OCCT n'est **pas valide**
(`BRepCheck_Analyzer`). Il a 7 faces, et un volume de 6,3560 pour un solide
fusionné de 6,3415, où le congé n'ajoute que 3,8·10⁻⁴. Le test ne compare donc
pas : il vérifie que le résultat de gbs, relu par OCCT, est valide, a 6 faces
et le volume de gbs (6,34191), puis rapporte celui d'OCCT dans un message.

Python : `test_fillets` sur les 9 cas de la galerie, volumes exacts du coin et
de la boîte arrondie.

Suites complètes : 821 tests C++, 136 Python.

## 8. Limites

- Un coin ne raccorde que trois **congés** d'un même rayon, à un sommet de
  trois arêtes, toutes convexes ou toutes concaves. Les chanfreins, les rayons
  différents (successivement), plus de trois arêtes et les coins mixtes
  relèvent des raccords de sommet quelconques (feuille de route, § 9 du
  plan).
- L'arc d'un coin est porté par le congé (`CurveOnSurface`), sa pcourbe sur la
  sphère est interpolée sur 33 points : exacte pour un congé exact, à
  `approx_tol` près pour un congé approché, l'écart absorbé dans la tolérance
  de l'arête.
- Sur une face fermée non périodique, seul un contact en **boucle** peut
  traverser la couture (arête fermée autour de la face). Un contact ouvert
  qui la traverse est refusé. Autour de deux telles faces, une arête fermée
  n'est congée que si leurs coutures se rencontrent dans une même section du
  congé ; une couture oblique est refusée.
- La fixture de l'aube existe en C++ (`testing/blade_fixture.h`) et en Python
  (exemples) : deux copies de la même construction, à garder alignées.

## 9. Points validés

1. **Coin seulement de trois congés** d'un même rayon, à un sommet de trois
   arêtes, convexes ou concaves (pas mêlées) ; le reste refusé et renvoyé aux
   raccords de sommet quelconques.
2. **Faces quelconques au coin** dès cette PR (pas seulement orthogonales), la
   rencontre des lignes des centres et des contacts servant de contrôle.
3. **Repère de la sphère** : x vers les pieds, z perpendiculaire, couture et
   pôles loin du coin.
4. **Vue sans couture** des faces B-splines fermées (`SeamlessSurface`), et
   départ de la boucle du congé sur la couture.
5. **Hermite local** pour la ligne des centres et les contacts, au lieu de
   l'interpolation cubique globale de la PR 5.
6. **Raffinement local** des seuls intervalles hors tolérance, contrôlés à la
   moitié des tolérances (milieu et 1/2 ∓ 1/(2√3)).
7. **Stations aux nœuds moins que C3** des faces.
8. **Tolérances des arêtes relevées** à celles des faces qu'elles bordent.
9. **Pied d'aube non comparé à OCCT**, dont le résultat est invalide : contrôle
   de validité et de volume de notre résultat relu par OCCT, et message.
10. **Fixture d'aube** dans `testing/` (en-tête de test partagé par les tests
    C++ et ceux d'OCCT).
11. **Bilan du palier** (§ 10 du plan) : le palier 4 est terminé ; la suite
    est la feuille de route, à commencer par les congés à rayon variable.

## 10. Revue de code

Douze points, tous corrigés dans cette PR.

| # | Point | Suite |
|---|---|---|
| 1 | La vue sans couture enveloppait toute surface fermée non élémentaire : le décalé d'un cylindre (bornes d'un tour) perdait le décalé exact d'`offset_surface`, et son congé, un tore sur master, devenait une B-spline approchée ; une surface de révolution, déjà périodique, passait par `into_bounds` | Seules les surfaces qui ne se répètent pas au-delà de leurs bornes sont enveloppées (`unrepeated`, § 4.2) ; test `fillet_offset_of_cylinder` |
| 2 | Autour de deux faces B-splines fermées, la boucle ne partait que de la couture de l'une | Coutures dans une même section : la ligne des centres se ferme, le congé passe ; couture oblique : refus explicite au lieu de « the centre line of a closed edge is not closed » ; coutures écartées le long de l'arête : deux arêtes, en chaîne (§ 4.2, § 8) ; test `fillet_two_closed_bsplines` |
| 3 | `into_bounds` n'admettait que 10⁻⁹ période de dépassement, la station de départ n'étant pas recalée sur la couture | Pôles ramenés sur les bornes à tol près dans l'espace (pas de paramètre fois la plus grande vitesse de la surface le long du contact) ; l'écart du contact au congé reste mesuré par `on_blend` |
| 4 | Repère de la sphère sans garantie : x vers la somme des pieds, un pied pouvant être à plus de 90° | `corner_frame` : pôle hors du triangle du coin, le plus loin des arcs ; couture dans le plus grand intervalle de longitudes ; orientation par la normale de la sphère au milieu d'un arc (§ 4.1) |
| 5 | L'écart de la pcourbe d'un arc entre les points interpolés n'était mesuré ni par `settle` ni par `check` | Confirmé sur un coin pointu : 1,4·10⁻⁴ pour une arête à 10⁻⁶, `check` passant. Interpolation doublée jusqu'à tol / 10 entre les points, tolérance de l'arête relevée à l'écart mesuré entre eux ; test `fillet_sharp_corner` |
| 6 | `into_face` recalculait `surface_closure` (39 évaluations) à chaque point | Périodes calculées une fois par `check_on_face` |
| 7 | Documentation de `skinned` restée sur l'ancien algorithme | Réécrite : Hermite, raffinement local, demi-tolérances, `breaks` |
| 8 | Documentation de `fillet` (C++ et Python) : « seulement des chaînes » | Chaînes ou coins de trois arêtes |
| 9 | Message « radii not the same? » alors que `fillet` n'a qu'un rayon | « do not meet within the blends continued », « contacts do not meet on a face » (§ 5) |
| 10 | `fillet_box_corner` ne fixait pas le centre de la sphère | Centre vérifié à 10⁻¹² |
| 11 | Bilan du palier : boîte aux douze arêtes « à 10⁻⁹ », le test à 10⁻⁸ | Bilan corrigé |
| 12 | Coin ou refus décidé par la chaîne `op == "fillet"` | Paramètre `fillets` de `blend_detail::blend` |
