# Palier 8 — PR 4 : congé tritangent en boucle, surfaces exactes tirées des stations

| | |
|---|---|
| PR | branche `feat/brep-full-round-2` (sur master a8eb9cb) |
| Document de référence | [brep_palier8_face_blends.md](brep_palier8_face_blends.md), § 3.5 et § 3.7 ; questions 13, 14 et 18 validées le 25 nov. 2026 (le bout arrondi, § 3.6, questions 15 à 17, part dans la PR 5) |
| Fichiers | `gbs-brep/full_round.h`, `gbs-brep/fillet.h` (`SectionOf::creases`, revue de code), `tests/tests_brep_full_round.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `python/gbsBindBrepOps.cpp` (docstring), `python/tests/test_brep_ops.py`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier8_face_blends.md` (statut, § 7, bilan), `News.md`, cette note |
| Statut | À valider (§ 9, § 10) |

## 1. Objectif

Trois ajouts au congé tritangent de la PR 3 (`fillet_full_round`), sans
changement d'API :

1. **La face du milieu M qui fait le tour** : ses deux arêtes avec les côtés
   sont fermées, elle n'a pas de face de bout. Le bord d'un disque ou d'une
   bride, la gorge d'un arbre (concave), l'extrémité d'un tube (M un anneau de
   deux contours), le bord d'une ailette (M la surface fermée d'un profil
   d'aube, une B-spline fermée non périodique).
2. **Les surfaces exactes tirées des stations** : cylindre, cône ou tore lus sur
   les centres et les rayons des boules, quelles que soient les faces, et
   contrôlés entre les stations.
3. **Deux défauts de la PR 3 corrigés** : la marche des stations qui enjambait un
   nez plus court qu'un pas (l'épingle de la revue de conception, refusée) ; la
   pcourbe du guide dans F₂ sur la copie, un segment.

S'y ajoutent les deux refus de repli de la question 18.

## 2. Sources

- Le § 3.5 et le § 3.7 du palier, revus après la revue de conception
  (`brep-review_palier8_pr4_conception.md`), et ses sondes (l'épingle, les
  nœuds C⁰ du profil).
- La marche de `fillet_law` (palier 5) et la vue au-delà d'une couture
  (`SeamlessSurface`), reprises.
- Les enveloppes de sphères : des boules de même rayon le long d'une droite
  enveloppent un cylindre, le long d'un cercle un tore ; de rayon linéaire le
  long d'une droite, un cône.

## 3. Démarche

### 3.1 M qui fait le tour

`setup` admet deux arêtes fermées. M est soit un anneau de deux contours
(chacun une arête fermée : l'extrémité d'un tube), soit un contour fermé par ses
coutures (le bord d'un disque, la gorge d'un arbre). Entre les deux arêtes, les
coutures seules : une autre arête y est refusée. Une arête fermée et l'autre non
est refusée.

La géométrie n'est pas prolongée (aucune face de bout à couper) ; la boucle part
du sommet du guide et se ferme : les boules en a et en b sont comparées, refus
au-delà de 10 `tol`. La topologie est celle d'une arête fermée des paliers 4
et 5 sur la copie réduite (`EndKind::Loop`). Le contrôle du contact sur F₂ se
fait sur toute la boucle.

Le départ à la couture d'un côté (question 14, codé après la revue de code,
§ 10, I-3) : là où le contact d'un côté B-spline fermé non périodique franchit
sa couture, la boucle part de là, comme `fillet_law`, le guide lu au-delà de
son sommet (`LoopCurve`). Tout côté fermé, périodique compris, a une arête de
couture, qui aboutit au sommet de son arête avec M : le contact doit la croiser
en un sommet de la boucle. Deux côtés dont les contacts croisent leurs coutures
à des stations différentes sont refusés (« the sides close across seams that
the contacts cross at different stations of the loop: not handled … align the
seams of the sides ») : la boucle n'a qu'un sommet, là où elle part.

*Écart à la question 14* (§ 9, point 4) : ce refus couvre le cas de la revue,
l'extrémité d'un tube aux deux parois fermées dont les coutures sont tournées
l'une par rapport à l'autre, courant en géométrie importée. `fillet_law` refuse
de même deux faces B-splines fermées dont les coutures ne sont pas dans la même
section.

### 3.2 La marche

Celle de `fillet_law` : 512 pas sur l'arête, 256 de chaque côté depuis son
milieu, chacun divisé par deux (jusqu'à 10 fois) là où Gauss-Newton échoue ou
la boule saute. La garde contre le saut ne grandit pas avec le pas : chaque
contact ne doit pas bouger de plus de 4 fois ce que le guide lui fait parcourir
sur le pas, son déplacement plus sa rotation fois la distance du contact au
guide, plus 10 `tol`, mesuré sur le pas divisé (revue de code, I-1 : au coin
intérieur d'un coude serré, le contact d'en face va bien plus vite que le
guide). Les stations que demandent la peau et les
contrôles partent de la station la plus proche de la marche (jusqu'à 8
divisions). Au-delà des sommets d'un guide ouvert, 32 pas divisés au plus 6 fois.

Les supports qui se ferment par une couture sans être périodiques sont vus
au-delà (`SeamlessSurface`), comme dans `fillet_law`.

La première station part du cercle inscrit dans les tangentes des trois faces,
dans le plan normal au guide en son milieu ; le point de l'autre arête dans ce
plan est, parmi ses traversées du plan (encadrées entre 256 points, puis
coupées en deux), la plus proche. (Un Newton depuis le plus proche de 33 points
partait de l'autre côté d'un tube carré, revue de code, I-1.)

Coût : la nervure droite passe de 2 à 6 ou 9 ms ; l'épingle, longue de 2,
demande 0,3 à 0,4 s ; le bout d'un tube carré, 70 à 140 ms.

### 3.3 Les surfaces exactes

Sur 17 stations, les centres sur une droite (sous `tol`) avec un rayon constant
(sous `tol`) : cylindre ; sur une droite, le rayon linéaire le long d'elle :
cône ; sur un cercle (par trois centres, sous `tol`), le rayon constant : tore.
La primitive est retenue seulement si :

- en 129 stations, une sur deux entre celles de la marche testée, les 9 points
  de l'arc exact sont sur elle à `approx_tol` / 2 et sa normale à
  `angular_tol` / 2 de celle de la boule ;
- les contacts, interpolés par les stations sur leurs faces (à 10⁻⁹ près des
  paramètres résolus), sont sur elle à `tol` en 257 points ;
- la ligne des centres s'interpole à `tol`.

Les pcourbes sont interpolées par morceaux, entre les nœuds du guide sous C³ et
là où un guide sur une surface (l'arête d'une révolution) en croise les nœuds :
la vitesse d'un cercle rationnel y casse, et les stations avec elle. Sinon, la
peau de la PR 3.

**Les plis de la peau** (revue de code, I-2). Le long d'un guide seulement G1
(une droite puis un arc : un saut de courbure), la peau ne convergeait pas
(refus « angular_tol »). Elle plie maintenant, ses dérivées prises d'un seul
côté (une peau C⁰ par morceaux, comme la section G2 de la PR #245), aux nœuds
du guide sous C², là où un guide sur une surface croise les nœuds de celle-ci
sous C², et là où un contact croise un nœud de sa face sous C². `fillet.h` le
prend par `SectionOf::creases`, sans autre changement.

### 3.4 Les refus de repli (question 18)

Aux contacts de 129 stations (celles du contrôle des surfaces, § 3.3), sur les
trois faces, la plus grande courbure principale de la face du côté de la boule
(`most_bending`, revue de code, m-2), comparée à 1/r ; sa direction dit lequel
des deux refus :

- **(a) en travers**, plus proche de l'arc de la section que du guide : la
  boule moins courbe que la face, l'arc en sortirait. Précaution sans cas
  atteignable (§ 6) : la boule inscrite disparaît avant de devenir plus large
  que la face ;
- **(b) le long**, plus proche du guide : la ligne des centres se replierait.

Le message donne le rayon de courbure trouvé et la face. Quand la marche
s'arrête (aucune boule), (b) est d'abord vérifié le long des arêtes des trois
faces au rayon de la première boule (`bending` dans la direction de l'arête) :
une face qui se replie n'a pas de boule là, et c'est ce refus que reçoivent
l'ailette trop épaisse (sur M) et le tube aux coins extérieurs plus petits que
la boule (sur un côté).

*Écart au palier* : le contrôle de repli de `fillet_law` (l'offset de toute la
face) a été essayé et retiré. Sur un U d'une seule face lisse, il voit le repli
du fond du U, loin du congé, et refuse à tort ; sur l'ailette, il ne voit pas le
nez du bord de fuite (son échantillonnage). La courbure aux contacts est plus
juste (§ 9, point 2).

### 3.5 La pcourbe du guide dans F₂

`guide_on_other` remplace le segment de la PR 3. Au paramètre t du guide, le
point de l'ancienne pcourbe où le plan normal au guide en t coupe l'autre arête,
interpolé par 2ᵏ + 1 points, doublés de 16 jusqu'à 10⁻⁹ de l'ancienne pcourbe.

Sur un guide ouvert, la corde de la PR 3 n'avait pas d'effet visible : le
secteur de disque (120° et 300°) passait déjà, ses faces coupées le long des
contacts. Le test écrit d'abord (`curved_open_guide`) échoue sur master par son
volume (une peau au lieu du tore, à 10⁻⁹), pas par la corde. Le défaut ne se voit que sur une boucle,
où le segment est de longueur nulle.

## 4. Ce qui ne change pas

- L'API et les autres congés : `blend.h` et `face_blend.h` ne changent pas ;
  `fillet.h` prend seulement les plis donnés par le tritangent
  (`SectionOf::creases`, § 3.3), vides pour les autres congés.
- Les tests de la PR 3 restent verts, sans changement de tolérance. Trois
  changent d'attente, parce que la PR 4 fait ce qu'ils refusaient :
  - le bord d'un disque, refusé (« closes on itself »), est construit (test
    `disc_rim`) ;
  - un dessus qui se rétrécit jusqu'à 0,005, refusé (« too narrow »), est
    construit (un cône exact) ;
  - la nervure courbe, une peau, est un tore exact (volume à 2,7·10⁻¹⁵ au lieu de
    3,9·10⁻⁸).

## 5. Tests

`tests_brep_full_round` : 36 cas (22 de la PR 3, 14 nouveaux, dont 4 de la
revue de code), à `tol` = 10⁻⁷ pour les volumes.

| Cas | Attendu | Mesuré |
|---|---|---|
| bord d'un disque (R = 1, t = 0,3) | demi-tore exact (rayons 0,85 et 0,15) ; Pappus | 5,4·10⁻¹⁰ ; 15 ms |
| extrémité d'un tube (0,7 à 1) | M fermée par sa couture ; tore exact ; Pappus | 2,0·10⁻¹⁰ |
| gorge d'un arbre (0,4 × 0,3) | concave ; tore exact ; matière ajoutée par Pappus | 1,4·10⁻¹⁰ |
| secteur de disque, 120° et 300° | guide ouvert en arc ; tore exact ; Pappus | 7·10⁻¹², 2·10⁻¹² (master : une peau) |
| nervure courbe (PR 3) | tore exact | 2,7·10⁻¹⁵ |
| épingle de la revue, a = 0,2 / 0,02 / 0,004 | acceptée (master refuse a = 0,004) ; volume ôté le long de la ligne médiane (3,5·10⁻⁶ à 4,5·10⁻⁶) | 3·10⁻¹¹, 1,1·10⁻¹⁰, 4,9·10⁻¹⁰ (la peau) ; 0,3 à 0,4 s |
| ailette au profil de `blade_fixture`, épaisseur 0,004 | boucle à travers la couture du bord de fuite, r = 0,002 partout, peau sous `approx_tol` | déviation 4,8·10⁻⁶ ; 49 ms |
| ailette de 0,01 (r = 0,005 > 0,0035) | refus (b), le long, sur M | message vérifié, modèle inchangé |
| dessus à arête vive au milieu (rayon 0,033 au sommet) | la boule ne touche pas le sommet mais les épaules : construit | solide valide |
| `bending` sur un cylindre de rayon 0,2 | 5 en travers vers l'axe, 0 le long, −5 de l'autre côté | à 10⁻⁹ |
| dessus B-spline bosselé de 2·10⁻⁵ entre deux des 17 stations | pas de cylindre : la peau | `exact` faux, solide valide |
| même face des deux côtés avec une boucle (auge tournée) | refus qui dit de couper la face | message vérifié |
| dessus qui se rétrécit jusqu'à 0,005 (PR 3, refusé) | cône exact ; volume ôté en forme fermée (revue, m-4) | 0,0262543, à 10⁻⁹ |
| nervure qui se rétrécit de 0,6 à 0,2 (PR 3) | volume ôté en forme fermée : le cône sur la section en x = 0 (deux coins d'une ellipse), homothétique depuis le sommet | 0,0373832, à 10⁻⁹ |
| bout d'un tube carré, coins 0,4 / 0,1, 0,4 / 0,02, 0,2 / 0,02, guide sur chaque paroi (revue, I-1) | construit dans les deux ordres ; coins concentriques : centres sur la ligne médiane, r = 0,15 | à 2,1·10⁻⁸ ; 70 à 140 ms |
| tube aux coins extérieurs 0,05 et 0,14 (sous r = 0,15), deux ordres (revue, m-1) | refus (b) sur un côté, nommé « first side » ou « second side » selon l'ordre ; modèle inchangé | message vérifié |
| paroi le long d'une droite puis d'un quart de cercle (guide G1, revue, I-2) | la peau plie au saut de courbure ; centres sur la ligne médiane, r = 0,15 ; volume ôté par Pappus | centres à 10⁻⁷ ; volume à 9,9·10⁻⁶ (la peau, la quadrature) ; 50 ms |
| bout de tube, parois B-splines fermées (carrées ou rondes) ou paroi intérieure cylindrique, coutures alignées ou tournées, deux ordres (revue, I-3) | alignées : construit ; à des stations différentes : refus nommé, modèle inchangé | 10 construits, 6 refusés |
| `most_bending` sur le cylindre de rayon 0,2 (revue, m-2) | 5 en travers vers l'axe, 0 le long de l'autre côté ; 2,5 en biais à 45° | à 10⁻⁹ |

OCCT (`tests_brep2_occt.full_round_loops`) : la gorge contre l'arbre coupé par
l'anneau de la gorge jusqu'aux centres et par le tore, mêmes faces (7), volume à
1,2·10⁻⁸ ; le disque contre le cylindre de rayon 0,85 fusionné avec le tore, à
10⁻⁸, OCCT gardant son tore en deux faces (gbs : une). L'écart de 7,8·10⁻⁷ du
disque (revue de code, m-4) ne venait pas de l'intégration d'OCCT : ses quatre
méthodes (`VolumeProperties` sans et avec `Eps` 10⁻⁹ et 10⁻¹², `VolumePropertiesGK`)
s'accordent à 10⁻¹² sur le solide relu. Il venait de l'écriture STEP : les
contacts, des courbes sur le tore exactes à 4,9·10⁻¹⁰, y sont approchés en
B-splines à `approx_tol` = 10⁻⁵ par défaut, et les plans qu'ils bornent relus
7,8·10⁻⁶ trop petits en aire. Écrit à `approx_tol` = 10⁻¹⁰ (paramètre ajouté à
`to_occt` et `expect_same` du test), le solide relu est à 3·10⁻¹⁰ de la forme
fermée.

Python : `test_fillet_full_round_loop` (le disque : tore, volume). Galerie :
`full_round_disc`.

## 6. Ce qui n'est pas testé

- **(a) en travers, de bout en bout.** Aucun cas naturel trouvé, et la revue de
  code le tient pour inatteignable : la boule inscrite disparaît avant de
  devenir plus large que la face (un dôme circulaire donne toujours r ≤ son
  rayon ; une arête vive est touchée à ses épaules). Gardé comme précaution ;
  le calcul (`bending`, `most_bending`) est testé seul ; le refus, non.
- **Le départ décalé à la couture d'un côté** (question 14) : codé, sans cas de
  bout en bout. Il faut un côté dont le contact croise sa couture ailleurs
  qu'à la station du sommet de son arête avec M (une couture en biais), et
  l'autre côté sans couture croisée ailleurs ; les tubes et les disques des
  tests ont leurs coutures dans le plan normal au guide. Les refus de coutures
  à des stations différentes, eux, sont testés.

## 7. Illustration

Note d'illustration (montrée à la validation, hors du dépôt) : le disque, la
gorge, l'extrémité du tube, l'ailette avant / après ; l'épingle ; le secteur de
disque.

## 8. Estimation

| Partie | Heures |
|---|---|
| boucle (setup, géométrie, topologie, pcourbe du guide) | 4 |
| marche, refus de repli | 4 |
| surfaces exactes contrôlées, ajustements | 4 |
| tests, OCCT, Python, galerie, note | 4 |
| **Total** | **16** |

## 9. Points à valider

1. **La marche de `fillet_law`** (§ 3.2), garde à 4 fois ce que le guide fait
   parcourir à chaque contact (déplacement et rotation) sur le pas divisé, pour
   tous les congés tritangents. *Recommandé* : elle corrige l'épingle et
   l'ailette, et passe les coudes serrés (revue de code, I-1) ; le coût passe de
   2 à 6 ou 9 ms sur une nervure droite.
2. **Les refus de repli par la plus grande courbure principale aux contacts**
   (§ 3.4), en travers et le long, au lieu du contrôle de l'offset de toute la
   face. *Recommandé* : l'offset voit des replis loin du congé (le fond d'un U)
   et manque le nez d'une ailette. Le refus (a), précaution, n'a pas de cas de
   bout en bout (§ 6).
3. **Les surfaces exactes contrôlées** (§ 3.3) : 17 stations pour reconnaître,
   129 pour contrôler, les contacts interpolés à 10⁻⁹ en paramètre, la ligne des
   centres à `tol`. *Recommandé.*
4. **Les coutures des côtés** (question 14, § 3.1) : le départ à la couture
   d'un côté B-spline fermé est codé ; deux côtés fermés dont les contacts
   croisent leurs coutures à des stations différentes sont refusés, avec un
   message qui dit d'aligner les coutures. C'est un **écart** au résultat
   attendu par la revue de code (I-3 : l'extrémité d'un tube aux coutures
   tournées d'un quart de tour, cas courant en géométrie importée) : la boucle
   n'a qu'un sommet, et chaque couture croisée en demande un. Le construire
   demanderait une boucle de deux arêtes (deux congés enchaînés, coupés aux deux
   stations), dans la topologie de `blend.h` : environ 6 h. *Recommandé* : le
   refus nommé dans cette PR, la boucle à deux sommets avec les ensembles de
   faces (PR 6), qui enchaîne déjà plusieurs guides. *Écarté* : un
   recoupement des faces importées par l'utilisateur, qu'il ne peut pas faire
   facilement.
5. **Les trois tests de la PR 3 qui changent d'attente** (§ 4) : ce qu'ils
   refusaient est construit. *Recommandé.*
6. **La comparaison du disque à OCCT à 10⁻⁸** (§ 5), l'écriture STEP du test à
   `approx_tol` = 10⁻¹⁰. *Recommandé* : l'écart de 7,8·10⁻⁷ était
   l'approximation des contacts à l'écriture, pas le congé.

## 10. Revue de code

Revue du commit 563b3d3 (25 nov. 2026) : trois constats importants, cinq
mineurs. Tout est corrigé, sauf I-3 en partie (un écart, § 9, point 4) et m-5,
laissé au bilan du palier comme demandé.

**I-1. La garde contre le saut refusait les coudes serrés**, et le résultat
dépendait de la face donnée comme F₁. Deux causes :

- la garde bornait le mouvement de chaque contact par 4 fois celui du guide ;
  au coin intérieur d'un coude plus serré que les deux tiers de la boule, le
  contact d'en face va plus vite, quel que soit le pas. Elle borne maintenant
  chaque contact par ce que le guide lui fait parcourir : son déplacement plus
  sa rotation fois la distance du contact au guide (§ 3.2) ;
- avec le guide sur la paroi extérieure, la première station elle-même
  échouait : le point de l'autre arête dans le plan normal au guide, cherché
  par Newton depuis le plus proche de 33 points, partait de l'autre côté du
  tube (la paroi opposée, à 1,4 de là). Il est cherché parmi les traversées du
  plan, encadrées entre 256 points (§ 3.2).

Test `square_tube_corners` : coins 0,4 / 0,1 (concentriques), 0,4 / 0,02 et
0,2 / 0,02, le guide sur chaque paroi, tous construits ; pour les coins
concentriques, centres sur la ligne médiane et r = 0,15 à 2,1·10⁻⁸ dans les deux
ordres. Le volume ôté n'y est qu'affiché : `signed_volume` à 512 points par
face s'écarte de 4·10⁻³ sur la peau qui fait le tour (5·10⁻⁴ à 2048 points),
la géométrie est la vraie mesure. L'épingle et l'ailette, que la garde
protégeait, passent toujours, mêmes écarts.

**I-2. Guide G1** : la peau plie aux nœuds sous C² du guide, des surfaces qu'il
suit et des faces des contacts (§ 3.3, `SectionOf::creases`). Test
`wall_along_a_g1_path` : la paroi de la revue (une droite puis un quart de
cercle de rayon 0,1, boule 0,15), construite ; centres sur la ligne médiane et
r = 0,15 à 10⁻⁷ ; volume ôté par Pappus à 9,9·10⁻⁶ (la peau et la quadrature).
Les tubes carrés à coins 0,4 et 0,1, refusés « angular_tol » par la revue,
passent (`square_tube_corners`).

**I-3. Départ à la couture d'un côté** (question 14) : codé (§ 3.1). Il ne
suffit pas au cas de la revue : dans un tube, chaque paroi fermée a une arête de
couture, qui aboutit au sommet de son arête avec le dessus, et le contact doit
la croiser en un sommet de la boucle. Coutures tournées l'une par rapport à
l'autre : deux sommets, refusé avec un message qui le dit (la paroi intérieure
cylindrique périodique comprise, refusée avant par la topologie avec le message
de la revue, qui ne nommait pas la cause). Deux corrections en chemin : un
contact posé exactement sur la couture à une station n'était pas vu (le
changement de signe strict), et le côté du guide, qui croise sa couture au
départ, n'était pas compté. Test `tube_seams` : 8 cas dans les deux ordres,
coutures alignées construites (B-splines carrées ou rondes, cylindre et
B-spline), tournées refusées, modèle inchangé. **Écart** au résultat attendu
(§ 9, point 4) : l'utilisateur tranche.

**m-1. Refus de repli.** (a) en travers gardé comme précaution, dit
inatteignable (§ 3.4, § 6). (b) sur un côté : une fois I-1 corrigé, le tube aux
coins extérieurs de 0,05 à 0,14 reçoit le refus (b) sur la paroi extérieure,
« first side » ou « second side » selon l'ordre. Test
`square_tube_folding_on_a_side` (0,05 et 0,14, deux ordres, modèle inchangé).

**m-2. Courbure aux contacts** : la plus grande courbure principale
(`most_bending`), aux 129 stations du contrôle des surfaces, sa direction
choisissant (a) ou (b) (§ 3.4). Test `bending` : sur le cylindre, 5 vers l'axe
en travers, 0 le long de l'autre côté, 2,5 en biais.

**m-3.** Réglé avant la revue (la primitive reconnue à tort).

**m-4. Volumes.** Le dessus qui se rétrécit jusqu'à 0,005 et la nervure de 0,6
à 0,2 ont leur volume ôté en forme fermée, à 10⁻⁹ (§ 5). Le disque contre OCCT
revient à 10⁻⁸ : l'écart venait de l'écriture STEP des contacts, pas d'OCCT
(§ 5). L'épingle : écarts absolus dans le tableau (3·10⁻¹¹ à 4,9·10⁻¹⁰).

**m-5.** Le congé d'arête du dessus de l'ailette refusé à r ≤ 0,001 (sur
master) n'est pas traité ici ; il reste au bilan du palier.

Tests : `tests_brep_full_round` 36 cas, 9 543 assertions (22 de la PR 3, 14 de
la PR 4 dont 4 de la revue) ; `tests_brep2_occt.full_round_loops` à 10⁻⁸ pour
le disque.
