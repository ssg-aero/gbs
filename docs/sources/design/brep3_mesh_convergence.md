# Maillage : convergence du raffinement sur les faces étroites et courbes

| | |
|---|---|
| PR | branche `fix/mesh-refinement-convergence` (sur master a8eb9cb) |
| Document de référence | [brep_palier3_mesh.md](brep_palier3_mesh.md), § 3 (arêtes), § 4.2 (mise à l'échelle), § 4.4 (raffinement), § 9.5 et § 10 ; [brep3_pr03_faces.md](brep3_pr03_faces.md) |
| Fichiers | `gbs-brep/mesh.h`, `tests/tests_brep_mesh.cpp`, `News.md`, `brep_palier3_mesh.md`, `.github/workflows/ci.yml` (§ 8), cette note |
| Statut | À valider (§ 7) |

## 1. Le défaut

CadSmith maille les pièces avec les options par défaut (flèche de 1e-3 de la
diagonale, angle de 0,5, angle minimal de 20°). Sur une pale d'hélice (un modèle
privé, hors du dépôt), deux faces B-spline ne convergeaient pas : environ 400 000
triangles chacune, la limite `max_points` atteinte, près de 5 s chacune, la face
rendue incomplète. Une flèche quatre fois plus fine les maillait en moins d'une
milliseconde. NASA Rotor 37 (`rotor37.step`, hors du dépôt lui aussi) montrait la
même chose sur sa face 3 : 399 982 triangles en 3,6 s à la flèche par défaut,
7 901 triangles (toute la pièce) en 25 ms à une flèche quatre fois plus fine.

Le raffinement n'était donc ni convergent ni monotone : une flèche plus grossière
donnait bien plus de triangles qu'une plus fine. Une sonde (§ 5.3 : 41 flèches de
0,1 à 10 fois le défaut, trois angles, avec et sans angle minimal, sur les STEP
de `tests/in/step`, rotor37 et la pale) trouve 316 maillages arrêtés à
`max_points` sur 4 100, 591 incomplets et 87 repliés. Ces défauts ne sont pas
propres aux pales : `occt_loft.stp` a cinq faces incomplètes à toutes les flèches,
et `assembly.stp`, `torus_band.stp` ou `cylinder.stp` en ont avec un angle de 0,35.

## 2. Les faces en cause

Les faces qui explosent sont des **bandes étroites, très courbes en travers et
longues, courbées aussi dans leur longueur** : le bord de fuite arrondi de
rotor37 (face 3 : 0,26 à 0,42 mm de large, 75 mm de long) et les bords d'attaque et
de fuite de la pale. En travers, l'angle demande des triangles très étroits :
0,1 mm sur rotor37 (la normale y tourne de π sur la largeur). Le long de la bande,
les arêtes sont presque droites : leur discrétisation (§ 3 du palier) ne regarde
que la courbe, et la flèche de 0,74 mm laisse des segments de 30 à 45 mm.

Le diagnostic (des triangles de la face classés, à l'arrêt, par le critère qui
les raffine encore) donne sur rotor37, face 3, à 5 000 points : **1 067 triangles
repliés**, tous collés au bord u = u_max (les trois coins à u = 0,500849 à six
chiffres), et 2 058 triangles de qualité. Aucun ne dépasse la flèche ni l'angle.
Les faces de la pale donnent la même chose.

## 3. La cause

### 3.1 Un triangle replié contre un côté du bord ne se déplie jamais

Le bord d'une face n'est jamais raffiné par sa triangulation : ses points sont
ceux des arêtes, partagés avec les autres faces. La corde d'un segment de bord
s'écarte de la surface d'autant que la flèche le permet, ici de 0,30 mm (le long
de la normale) sur un segment de 45 mm de rotor37. Un point de la surface plus
près du bord que cet écart est **derrière la corde** en 3D : le triangle qu'il
fait avec le segment est replié contre la surface.

Le raffinement donnait à un triangle replié son centre de gravité. Pour un
triangle appuyé sur un côté contraint, c'est le pire choix :

1. le centre de gravité est plus près du côté que le sommet opposé ;
2. la Delaunay prend pour sommet d'un côté contraint le point qui le voit sous le
   plus grand angle, donc le plus proche de son milieu : le nouveau point ;
3. le nouveau triangle sur le côté est donc plus replié encore, et reçoit son
   centre de gravité à son tour.

Les points s'accumulent contre le côté jusqu'à `max_points`. Aucun point
intérieur ne déplie ce triangle : seul un point sur le côté le ferait, et le bord
n'est pas raffiné.

### 3.2 Pourquoi une flèche plus grossière explose et une plus fine converge

Avec la flèche, les cordes des arêtes s'écartent davantage de la surface (comme
L² / R), alors que la largeur des triangles en travers ne bouge pas : l'angle la
fixe (0,1 mm). À un quart de la flèche, la corde du bord de rotor37 s'écarte de
0,14 mm, sous la largeur de la bande : il y a de la place pour un sommet qui ne
replie pas, et le raffinement n'a pas à s'approcher du bord. À la flèche par
défaut, l'écart (0,30 mm) dépasse la largeur des triangles.

L'issue dépend alors du hasard du raffinement : si un point, inséré pour la
qualité ou pour l'angle, tombe dans la bande derrière la corde, la cascade part.
Sur la face 3 de rotor37 (sonde à 20 000 points), selon la flèche :

| flèche / défaut | 0,1 | 0,25 | 0,4 | 0,5 | 0,6 | 0,7 | 0,8 | 0,9 | 1 | 1,25 | 1,5 | 2 | 4 | 8 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| triangles | 166 | 110 | 102 | 39 982* | 242* | 39 982* | 230* | 39 982* | 39 982* | 80 | 258* | 406* | 682* | 1 168* |

(* incomplète.)

### 3.3 La mise à l'échelle ne bornait pas l'allongement des triangles

La triangulation se fait dans (a u, b v), a et b donnés par l'espacement que la
courbure moyenne demande le long de u et de v (§ 4.2 du palier). Sur la bande de
rotor37, cela fait des triangles de 0,1 mm sur 33 mm : un rapport de 330. Le long
de la bande, la corde d'un tel triangle quitte la surface de k L² / 8 (jusqu'à la
flèche), bien plus que sa largeur : les triangles intérieurs aussi se dressent
contre la surface, et se replient. La courbure qui compte pour cette flèche est
celle de l'iso-ligne en 3D, y compris dans le plan tangent (courbure
géodésique) : sur un demi-tube de rayon 0,2 autour d'un arc de rayon 100, le grand
cercle du haut se courbe dans le plan tangent, et la seule rotation des normales,
que la mise à l'échelle mesurait, ne le voit pas.

### 3.4 Les arrêts sans explosion

Les faces 4 et 5 de rotor37 s'arrêtaient incomplètes avec quelques centaines de
triangles (354 et 603) : entre les deux bouts d'un segment de bord, la normale de
la surface tournait de 0,511 et 0,534, plus que l'angle. La discrétisation de
l'arête ne regarde que la rotation de la tangente de la courbe, pas celle de la
normale de la face. Le critère d'angle, entre deux coins d'un côté contraint, ne
pouvait être tenu par aucun point intérieur. Même cas pour `occt_loft.stp` (cinq
faces sur six) et, avec un angle de 0,35, pour des cylindres, des tores et des
sphères des tests.

### 3.5 Le sommet d'un cône pesait dans la courbure moyenne

La courbure moyenne est la rotation totale des normales sur la longueur totale ;
un pôle n'y pèse rien, car sa normale est nulle (« 0 where a normal vanishes »).
Au sommet de `cone_noseam.stp`, le paramètre dépasse le sommet d'un arrondi :
|S_u| vaut 2,5e-17 au lieu de 0, et la normale, longue de 3e-17, est **retournée**.
Chaque colonne de la grille ajoutait un demi-tour : le cône paraissait courbé le
long de ses génératrices, et recevait des triangles isotropes au lieu des
triangles allongés du cylindre. Ce défaut, ancien, devenait visible avec les
critères des arêtes du § 4.1 (la génératrice de couture découpée en 28 segments).

## 4. Le correctif

Le principe : **le bord d'une face doit déjà satisfaire ce que sa triangulation
demande de ses côtés**, puisqu'elle ne le raffine pas, et la métrique ne doit pas
demander de triangles qui se dressent contre la surface.

### 4.1 Les arêtes, face par face (`MeshBuilder::on_faces`)

Chaque segment d'une arête (`curve_parameters`, puis la déviation sur ses faces
comme avant) est coupé en deux tant que, sur l'une de ses faces :

1. **il est plus long que 1 dans l'espace mis à l'échelle de la face**
   (`parameter_scales` sur la boîte de ses pcurves, 33 points par co-arête,
   calculée une fois par face : les échelles de sa triangulation, § 8) : le bord
   est à la taille des triangles que la face veut contre lui ;
2. **les normales de la surface à ses bouts tournent de plus que l'angle** (une
   normale négligeable, à un pôle, ne compte pas) : le critère d'angle des
   triangles le long du côté est tenu ;
3. **sa corde s'écarte de la surface de plus que `tilt` fois l'espacement que la
   face demande là** (en plus de la flèche, comme avant) : `local_spacing`, la
   courbure au milieu du segment (la plus grande rotation de la normale selon u
   ou v sur un pas d'un millième du segment), puis √(8 flèche / k) ou angle / k.
   `tilt` vaut tan(angle / 2) : un triangle aussi large que cet espacement est
   incliné d'au plus la moitié de l'angle sur la surface, loin de se replier
   (avec les normales qui tournent de l'angle en travers, au plus 1,5 fois
   l'angle).

Le critère 3 est local, pas moyen : près du bord d'attaque de rotor37, la face de
pale voisine est bien plus courbée qu'en moyenne, et l'espacement moyen (20 mm)
laissait une corde à 0,29 mm de la surface là où les triangles font 0,1 mm (18
maillages explosés sur 41 avec un angle de 0,2). Le pas de l'estimation est petit
devant le segment : sur toute sa longueur, il lisse la courbure et manque la
bande (5 maillages incomplets de la pale avec un angle de 0,2).

### 4.2 La métrique (`mesh_detail::parameter_scales`)

1. **L'allongement est borné** : l'espacement le long de u est réduit pour que la
   flèche de sa corde sur l'iso-ligne, k_u s_u² / 8, reste sous `tilt` fois
   l'espacement le long de v (de même pour v). k est la courbure en 3D de
   l'iso-ligne (la rotation de S_u, ou celle des normales si elle est plus
   grande), mesurée comme la courbure normale : rotation totale sur longueur
   totale, sur la grille de 9 × 9. Le borne ne joue que là où la flèche est grande
   devant le rayon en travers : sur un cylindre ou une sphère de taille
   ordinaire, rien ne change.
2. **Une normale ou une tangente négligeable** devant les longueurs moyennes
   (moins de 1e-6, comme le pôle de la triangulation) est nulle : un pôle ou un
   sommet ne pèse plus rien, comme le voulait le § 4.2 du palier.

### 4.3 La triangulation (`MeshBuilder::triangulate`)

Le critère d'angle ne compare plus les normales aux deux bouts d'un côté
contraint, comme le critère de flèche ne teste déjà pas le milieu d'un côté du
bord : aucun point intérieur ne les change, et § 4.1 (2) les tient (16 divisions
au plus). Sans cela, un segment que la discrétisation laisse hors de l'angle (au
bout des 16 divisions) arrêtait la face sans rien insérer.

Un triangle sans côté sur le bord, replié à un sommet (la normale de la surface y
fait plus de 90° avec la sienne) sans l'être en son centre, est raffiné comme pour
la flèche : `mesh` signale le repli aux sommets (§ 8).

### 4.4 Essayé puis écarté

| Essai | Résultat |
|---|---|
| Relever `max_points` | ne converge pas : les points s'accumulent contre le bord |
| Ne plus forcer le centre de gravité d'un triangle replié appuyé sur un côté contraint | arrête la cascade du § 3.1, mais pas celle qui suit, entre points intérieurs déjà collés au bord ; face incomplète et repliée ; sans effet une fois les arêtes conformes |
| Refuser tout point qui ferait, avec un côté contraint du triangle qui le contient, un triangle replié | ne se déclenche presque jamais (62 côtés contraints examinés pour 20 000 candidats) : la cascade passe par des triangles qui ne touchent pas le bord |
| Pas de point forcé dans un triangle d'angle minimal sous 1° ou 5° (à l'échelle) | sans effet sur les explosions restantes, des faces incomplètes en plus |
| Critère 3 avec l'espacement moyen de la face | voir § 4.1 |
| Allongement borné par la seule rotation des normales | le demi-tube du test (§ 6) reste replié : la courbure de ses iso-lignes est géodésique |

## 5. Mesures

clang 21, `-O2`, Apple M5 Max, un processus à la fois. Les temps sont ceux du
maillage seul (sans la lecture du STEP).

### 5.1 La pale et rotor37, options par défaut

| Pièce | Avant | Après |
|---|---|---|
| pale, les deux faces de bord | 399 986 et 399 985 triangles, 4 744 et 4 838 ms, incomplètes | 114 et 61 triangles, moins de 1 ms |
| pale, pièce entière (22 faces) | 801 050 triangles, 9 585 ms, 2 faces incomplètes | 1 230 triangles, 4 ms |
| rotor37, face 3 (bord de fuite) | 399 982 triangles, 3 586 ms, incomplète | 292 triangles, 1 ms |
| rotor37, faces 4 et 5 | 354 et 603 triangles, incomplètes | 39 et 226 triangles |
| rotor37, pièce entière (6 faces) | 403 480 triangles, 3 600 ms, 3 faces incomplètes | 2 855 triangles, 9 ms |

À un quart de la flèche, la pale donnait avant 104 et 76 triangles pour ses deux
faces de bord : elle en a maintenant 114 et 61 à la flèche par défaut.

### 5.2 Selon la flèche

Pièce entière, triangles et temps, avant → après (aucune face incomplète ni
repliée après ; avant, de 1 à 4 faces incomplètes) :

| flèche / défaut | 0,25 | 0,5 | 1 | 2 | 4 |
|---|---|---|---|---|---|
| rotor37 | 7 901 → 7 126 (26 → 22 ms) | 404 106 → 4 310 (4 062 → 14 ms) | 403 480 → 2 855 (4 247 → 10 ms) | 402 844 → 1 652 (5 498 → 6 ms) | 403 927 → 1 931 (6 069 → 7 ms) |
| pale | 3 560 → 3 482 (8 → 9 ms) | 401 928 → 1 898 (5 409 → 6 ms) | 801 050 → 1 230 (10 444 → 4 ms) | 800 638 → 884 (10 129 → 4 ms) | 800 422 → 762 (7 839 → 3 ms) |

Les 23 STEP des tests, à la flèche par défaut : mêmes nombres de triangles, sauf
les cônes (1 872 → 206 ; 1 864 → 206 sans couture), les sphères (`sphere.stp`
3 662 → 3 492, `sphere_cut_off_axis.stp` 4 738 → 3 344, `sphere_cut_through_poles.stp`
2 674 → 2 182, `half_spheres.stp` 3 646 → 3 440, `sphere_wedge.stp` 2 496 → 2 376,
`sphere_cap.stp` 1 726 → 1 700) et `occt_loft.stp` (2 794 → 828, complet au lieu
de cinq faces incomplètes, 17,3 → 2,6 ms). Les temps des autres sont ceux d'avant
à 15 % près (le plus lent : `torus_band.stp`, 2,8 → 3,2 ms, les critères des
arêtes).

### 5.3 La sonde

41 flèches de 0,1 à 10 fois le défaut (progression géométrique), pour un angle de
0,5, 0,35 et 0,2 avec l'angle minimal de 20°, et 0,5 sans angle minimal ; limite
de 20 000 points par face ; 25 fichiers (les 23 STEP de `tests/in/step`, rotor37,
la pale), 4 100 maillages. « Explosé » : une face à la limite. « Inversion » : le
plus grand rapport entre les triangles d'une flèche et ceux d'une flèche plus
fine, par fichier et par jeu d'options.

| | Incomplets | Repliés | Explosés | Pire inversion | Temps total |
|---|---|---|---|---|---|
| avant | 591 | 87 | 316 | 86,6 (pale, sans angle minimal) | 401 s |
| après | 2 | 6 | 0 | 2,21 (rotor37) | 17 s |

Sur les 23 STEP des tests seuls, après : aucun maillage incomplet, replié ni
explosé (avant : 288 incomplets, 2 repliés, 88 explosés, dont `occt_loft.stp` à
toutes les flèches) ; sur rotor37 et la pale : 2, 6 et 0 (avant : 303, 85 et 228).

Chaque partie du correctif, retirée seule de la version finale :

| Sans | Incomplets | Repliés | Explosés |
|---|---|---|---|
| (rien : la version finale) | 2 | 6 | 0 |
| l'allongement borné (§ 4.2, 1) | 2 | 226 | 0 |
| la longueur des segments (§ 4.1, 1) | 13 | 4 | 9 |
| l'écart des cordes local (§ 4.1, 3) | 44 | 28 | 40 |
| les normales négligeables (§ 4.2, 2) | 11 | 6 | 9 |
| l'angle aux bouts des segments (§ 4.1, 2) | 2 | 6 | 0 |
| l'angle hors des côtés contraints (§ 4.3) | 2 | 6 | 0 |

Les deux derniers ne changent aucun maillage de la sonde : une fois les segments
à la taille de la face et proches de la surface, leurs normales tournent de moins
que l'angle. Ils restent la garantie de l'option `angle` le long du bord et de
l'arrêt (§ 4.3).

### 5.4 Étanchéité et déviation

- Étanchéité : sur chaque fichier et chaque flèche de 0,25 à 4 fois le défaut,
  les côtés ouverts du maillage sont exactement les segments des arêtes libres
  (0 pour les solides ; 415 pour rotor37 à la flèche par défaut, qui est une
  coque ouverte).
- Déviation (échantillonnée sur 28 points par triangle, les pôles comme dans les
  tests) : inchangée sur les STEP des tests (1,04 fois la flèche au plus, les
  4 % du § 9.5 du palier). Sur la grande face B-spline de la pale (face 15), 1,29
  fois la flèche à 2 fois le défaut (1,10 avant, à 0,5 fois) : le critère ne
  teste que quatre points par triangle, et le point le plus écarté de ces
  triangles allongés n'en est aucun (§ 7, 3).
- Nombres de triangles : inchangés sur les plans, les cylindres, les tores, les
  B-splines des tests. Les cônes en ont 9 fois moins (1 872 → 206 pour
  `cone.stp` à la flèche par défaut, déviation 0,36 fois la flèche), comme le
  cylindre, le long de leurs génératrices (§ 3.5) ; les sphères 5 à 30 % de moins
  (déviation inchangée).

## 6. Tests

`tests_brep_mesh.narrow_curved_strips` : les quatre moitiés d'un tube de rayon 0,2
autour d'un arc de rayon 100 et de longueur 100 (une bande étroite courbée dans sa
longueur, comme le bord d'une pale), de 0,5 à 16 fois leur flèche par défaut
(0,096, la moitié du rayon). Chaque maillage est complet, non replié, à la flèche
près (5 %), de caractéristique d'Euler 1, et n'a pas plus de triangles (à 5 %
près) que celui d'une flèche plus fine ; de la plus fine à la plus grossière, à
moins de 1,5 fois. Avant : la moitié intérieure (ses arêtes se courbant vers elle)
s'arrêtait à la limite de 20 000 points dès la flèche par défaut (39 950
triangles) ; les autres se repliaient, avec de plus en plus de triangles (de 440
à 13 718). Après : de 434 à 520 triangles.

Le test échoue sur la version précédente (repli de la moitié intérieure à 0,5
fois la flèche) ; les 16 autres tests de `tests_brep_mesh` passent avant comme
après. `curved_faces` a ses nombres de triangles mesurés mis à jour (sphère : 1 186
et 12 466 au lieu de 1 276 et 13 540 ; cône : 32 et 128, l'éventail exact du
§ 10 du palier, déjà mesuré alors). `deflection_sweep` donne maintenant 206
triangles pour toutes les flèches du cône (1 782 à 2 184 avant), 3 122 à 3 880 pour
la sphère (3 468 à 4 116).

Toute la suite de la CI (options du job macOS : Release, tests, liaisons Python,
rendu), en local : les 994 tests C++ et les 191 tests Python passent.

Revue (§ 8) : `narrow_curved_strips` prend aussi le rayon arrondi de 3e-13 (la
moitié inférieure repliée avant la revue, comme sur Windows) ;
`edges_within_the_angle_on_their_faces` (une face qui se tord autour d'une arête
droite : seul § 4.1 (2) coupe ses segments là ; sans lui, la face est incomplète)
et `surface_throwing_inside_its_box` (une surface qui lève une exception dans la
boîte de sa face : la face voisine est maillée, la sienne échoue avec
`NonGeneric`) sont nouveaux. `curved_faces` : sphère 1 186 et 12 476.

## 7. Limites et points à valider

1. **Les faces fines entre deux bords proches** : à 4 et 6,3 fois la flèche par
   défaut, le moyeu de rotor37, percé par l'empreinte du pied de pale (épaisse de
   0,3 à quelques millimètres), s'arrête incomplet (3 739 et 3 784 triangles, en
   18 ms). La corde d'un côté de l'empreinte s'écarte de la surface (jusqu'à
   3 mm, dans le plan tangent) plus que l'empreinte n'est épaisse : le contour se
   croise en 3D. Aucun critère local de courbure ne le voit ; il faudrait
   comparer l'écart des cordes à la distance aux autres bords de la face, ou que
   la triangulation d'une face puisse demander la coupe d'un segment de son bord
   (et remailler les faces de l'arête). *Recommandé* : une PR à part si le cas
   se présente à des flèches d'usage ; à la flèche par défaut, rotor37 est
   complet.
2. **Six maillages repliés à la limite** sur 4 100 (rotor37 et la pale, avec un
   angle de 0,35, ou à plus de 1,5 fois la flèche) : des triangles presque
   plats en 3D (côtés 15,4, 32,6 et 48), dont la normale fait 90,1 à 91,4° avec
   celle de la surface à un coin. La triangulation teste le repli au centre de
   gravité, `mesh` le signale aux coins. *Revu* (§ 8) : la triangulation teste
   aussi les sommets d'un triangle sans côté sur le bord (la CI Windows l'a
   montré sur le demi-tube du test) ; restent les triangles appuyés sur un côté
   du bord (rotor37, 4 maillages sur 164 dans la sonde de la revue, dont 3 par
   un côté de 43,6 mm sur le moyeu, la limite 1).
3. **Déviation échantillonnée** jusqu'à 1,29 fois la flèche sur une grande face
   B-spline de la pale (§ 5.4) : le critère à quatre points (§ 9.5 du palier) ;
   le corriger (plus de points, ou les quarts des côtés) changerait tous les
   maillages. *Recommandé* : une PR à part.
4. **Résultats changés** : des arêtes ont plus de points sur les faces très
   anisotropes ; les cônes et les sphères moins de triangles. Les comparaisons à
   OCCT du § 9.2 du palier (cylindre, cône, sphère, à 1e-3 et 1e-2) sont à
   refaire.
5. **Une seule mise à l'échelle par face** (§ 9.5 du palier) : inchangé ; le
   critère 3 du § 4.1 est local, mais la triangulation reste dans l'espace moyen
   de la face.

## 8. Revue de code

Revue de la PR (sept constats), puis la CI, cassée sur trois jobs.

| Constat | Correction |
|---|---|
| **Une exception pouvait sortir de `mesh`** : la métrique d'une face (`parameter_scales` sur la grille de sa boîte) et les dérivées de `on_faces` étaient évaluées pour chaque face de l'arête, hors du `try` de `face`, et `uses_` couvre tout le modèle. Une surface qui lève une exception dans la boîte de sa face (une `OffsetSurface` sans normale là) faisait lever `mesh` d'une autre face, voisine | `metric` rend nul quand la surface lève (mémorisé) ; `on_faces` n'applique alors, pour cette face, que la flèche, comme avant (une exception dans les critères de la face aussi). `face` mesure à nouveau la métrique dans son `try` : la face échoue avec `NonGeneric`. Test `surface_throwing_inside_its_box` (lève avant) |
| **Deux métriques par face** : la discrétisation des arêtes prenait la boîte des pcurves (9 points par co-arête), la triangulation celle des points du bord | une seule métrique, celle de `metric` (33 points par co-arête), prise par la triangulation ; mesuré ci-dessous |
| **Le critère d'angle des côtés du bord, non testé** : retirer § 4.1 (2) ne changeait aucun maillage de la sonde ni des tests ; et sa garantie s'arrête à 16 divisions | test `edges_within_the_angle_on_their_faces` (échoue sans le critère : face incomplète) ; le commentaire dit les 16 divisions |
| Débordement : `previous + previous / 20` au premier tour du test | `std::optional` |
| `News.md` : 3,8 s au lieu de 3,6 s ; `MeshBuilder::on_faces`, privé | 3,6 s ; `MeshBuilder::edge` |
| `brep_palier3_mesh.md` § 3, § 4.2, § 4.4 décrivaient les règles d'avant | révisés, renvoyant à cette note |
| Trois boucles identiques pour annuler les vecteurs négligeables | `mesh_detail::or_null` et `negligible`, `std::views::join` ; partagés avec `normal_or_null` et le pôle de la triangulation |

**La CI.** Trois jobs de la PR ont échoué :

- *Windows (clang)* : `narrow_curved_strips`, la moitié inférieure repliée à la
  flèche par défaut. La libm de Windows arrondit sans doute autrement les points : un
  triangle en écharde (côtés 4,21, 2,33 et 6,54, le long d'une corde qui
  s'écarte de 0,053 de la surface, deux fois la largeur de la bande sous lui) se
  dresse à 85,5° en son centre, ce que la triangulation accepte, et à 92° à un
  sommet, ce que `mesh` signale. Reproduit sur macOS en arrondissant le rayon
  (3e-13 près) : 741 maillages repliés sur 5 808 (4 moitiés, 12 rayons, 121
  flèches de 0,5 à 16 fois le défaut). Correctif : un triangle sans côté sur le
  bord, replié à un sommet, est raffiné (§ 4.3) ; 0 sur 5 808. Un seuil plus
  strict laisse plus de maillages incomplets dans la sonde ci-dessous (6 au lieu
  de 3 à 60°, rotor37 ; 21 à 1,5 fois l'angle, rotor37 et des sphères) : 90°,
  comme `mesh`.
- *Ubuntu* : runner perdu pendant la compilation, 49 min après le départ (Ninja
  lance 6 compilations sur 4 cœurs et 16 Go, celles qui incluent le noyau BREP
  prennent 3 à 4 Go chacune). Trois d'entre elles à la fois, quatre compilations
  en tout (`GBS_HEAVY_COMPILE_JOBS=3`, `build-jobs: 4`), comme Windows / MSVC et
  macOS.
- *macOS* : arrêté au plafond de 60 min pendant les tests Python (51 à 59 min sur
  les dernières PR). Le plafond par défaut passe à 90 min (MSVC garde 120).

**Mesures après la revue** (clang 21, `-O2`, Apple M5 Max ; la sonde du § 5.3 sur
les 23 STEP des tests et rotor37, sans la pale, privée ; 3 936 maillages) :

| | Incomplets | Repliés | Explosés | Temps total |
|---|---|---|---|---|
| avant la revue | 3 | 4 | 0 | 15,6 s |
| après | 2 | 4 | 0 | 15,8 s |

Le critère de repli aux sommets seul ne change aucun de ces 3 936 maillages ; la
métrique unique en change 1 130, 94 de plus de 5 % (`plate_inch.stp`, rotor37,
`closed_bspline_across.stp` surtout), de 0,47 à 1,26 fois. Aux options par défaut
: mêmes nombres de triangles pour 18 fichiers sur 24, sinon `plate_inch.stp`
94 → 108, `sphere.stp` 3 492 → 3 500, `half_spheres.stp` 3 440 → 3 436,
`sphere_wedge.stp` 2 376 → 2 320, `band_cylinder.stp` 446 → 448,
`torus_band.stp` 1 240 → 1 242, rotor37 2 855 → 2 835. Une face incomplète de
rotor37 (3,98 fois le défaut) est maintenant complète. Les mesures de la pale
(§ 5.1, § 5.2) sont celles d'avant la revue.

Toute la suite de la CI en local, après la revue (options du job macOS : Release,
tests, liaisons Python, rendu) : les 997 tests C++ et les 191 tests Python passent.
