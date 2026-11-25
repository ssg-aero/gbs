# Palier 8 — PR 4 : congé tritangent en boucle, surfaces exactes tirées des stations

| | |
|---|---|
| PR | branche `feat/brep-full-round-2` (sur master a8eb9cb) |
| Document de référence | [brep_palier8_face_blends.md](brep_palier8_face_blends.md), § 3.5 et § 3.7 ; questions 13, 14 et 18 validées le 25 nov. 2026 (le bout arrondi, § 3.6, questions 15 à 17, part dans la PR 5) |
| Fichiers | `gbs-brep/full_round.h`, `tests/tests_brep_full_round.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `python/gbsBindBrepOps.cpp` (docstring), `python/tests/test_brep_ops.py`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier8_face_blends.md` (statut, § 7, bilan), `News.md`, cette note |
| Statut | À valider (§ 9) |

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

*Écart à la question 14* : le départ à la couture d'un côté B-spline fermé non
périodique, quand un contact la franchit, n'est pas codé. Les côtés des cas
testés sont des plans ou des surfaces de révolution, et la boucle part du sommet
du guide (§ 9, point 4).

### 3.2 La marche

Celle de `fillet_law` : 512 pas sur l'arête, 256 de chaque côté depuis son
milieu, chacun divisé par deux (jusqu'à 10 fois) là où Gauss-Newton échoue ou
la boule saute. La garde contre le saut ne grandit pas avec le pas : les contacts
ne doivent pas bouger de plus de 4 fois ce que porte le guide sur le pas, plus
10 `tol`, mesuré sur le pas divisé. Les stations que demandent la peau et les
contrôles partent de la station la plus proche de la marche (jusqu'à 8
divisions). Au-delà des sommets d'un guide ouvert, 32 pas divisés au plus 6 fois.

Les supports qui se ferment par une couture sans être périodiques sont vus
au-delà (`SeamlessSurface`), comme dans `fillet_law`.

Coût : la nervure droite passe de 2 à 6 ou 9 ms ; l'épingle, longue de 2,
demande 0,4 à 0,5 s.

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

### 3.4 Les refus de repli (question 18)

Aux contacts de 65 stations, sur les trois faces, la courbure normale de la face
dans une direction, du côté de la boule (`bending`), comparée à 1/r :

- **(a) en travers**, le long de l'arc de la section : la boule moins courbe
  que la face, l'arc en sortirait ;
- **(b) le long**, le long du contact d'une station à la voisine : la ligne des
  centres se replierait.

Le message donne le rayon de courbure trouvé et la face. Quand la marche
s'arrête (aucune boule), (b) est d'abord vérifié le long des arêtes des trois
faces au rayon de la première boule : une face qui se replie n'a pas de boule
là, et c'est ce refus que reçoit l'ailette trop épaisse.

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

- L'API et les autres congés : `blend.h`, `fillet.h`, `face_blend.h` ne
  changent pas.
- Les tests de la PR 3 restent verts, sans changement de tolérance. Trois
  changent d'attente, parce que la PR 4 fait ce qu'ils refusaient :
  - le bord d'un disque, refusé (« closes on itself »), est construit (test
    `disc_rim`) ;
  - un dessus qui se rétrécit jusqu'à 0,005, refusé (« too narrow »), est
    construit (un cône exact) ;
  - la nervure courbe, une peau, est un tore exact (volume à 2,7·10⁻¹⁵ au lieu de
    3,9·10⁻⁸).

## 5. Tests

`tests_brep_full_round` : 32 cas (22 de la PR 3, 10 nouveaux), à `tol` = 10⁻⁷
pour les volumes.

| Cas | Attendu | Mesuré |
|---|---|---|
| bord d'un disque (R = 1, t = 0,3) | demi-tore exact (rayons 0,85 et 0,15) ; Pappus | 5,4·10⁻¹⁰ ; 15 ms |
| extrémité d'un tube (0,7 à 1) | M fermée par sa couture ; tore exact ; Pappus | 2,0·10⁻¹⁰ |
| gorge d'un arbre (0,4 × 0,3) | concave ; tore exact ; matière ajoutée par Pappus | 1,4·10⁻¹⁰ |
| secteur de disque, 120° et 300° | guide ouvert en arc ; tore exact ; Pappus | 7·10⁻¹², 2·10⁻¹² (master : une peau) |
| nervure courbe (PR 3) | tore exact | 2,7·10⁻¹⁵ |
| épingle de la revue, a = 0,2 / 0,02 / 0,004 | acceptée (master refuse a = 0,004) ; volume ôté le long de la ligne médiane | 7·10⁻⁶, 3·10⁻⁵, 1,4·10⁻⁴ relatifs (la peau) ; 0,4 à 0,5 s |
| ailette au profil de `blade_fixture`, épaisseur 0,004 | boucle à travers la couture du bord de fuite, r = 0,002 partout, peau sous `approx_tol` | déviation 4,8·10⁻⁶ ; 49 ms |
| ailette de 0,01 (r = 0,005 > 0,0035) | refus (b), le long, sur M | message vérifié, modèle inchangé |
| dessus à arête vive au milieu (rayon 0,033 au sommet) | la boule ne touche pas le sommet mais les épaules : construit | solide valide |
| `bending` sur un cylindre de rayon 0,2 | 5 en travers vers l'axe, 0 le long, −5 de l'autre côté | à 10⁻⁹ |
| dessus B-spline bosselé de 2·10⁻⁵ entre deux des 17 stations | pas de cylindre : la peau | `exact` faux, solide valide |
| même face des deux côtés avec une boucle (auge tournée) | refus qui dit de couper la face | message vérifié |
| dessus qui se rétrécit jusqu'à 0,005 (PR 3, refusé) | cône exact | solide valide |

OCCT (`tests_brep2_occt.full_round_loops`) : la gorge contre l'arbre coupé par
l'anneau de la gorge jusqu'aux centres et par le tore, mêmes faces (7), volume à
1,2·10⁻⁸ ; le disque contre le cylindre de rayon 0,85 fusionné avec le tore.
Pour le disque, OCCT garde son tore en deux faces (gbs : une) et intègre le
solide de gbs relu 7,8·10⁻⁷ sous le sien, là où gbs trouve la forme fermée à
5·10⁻¹⁰ : la comparaison est à 10⁻⁶, avec moins de faces admises. La cause de
l'écart dans OCCT n'est pas établie (la relecture d'un tore d'un tour entier).

Python : `test_fillet_full_round_loop` (le disque : tore, volume). Galerie :
`full_round_disc`.

## 6. Ce qui n'est pas testé

- **(a) en travers, de bout en bout.** Aucun cas naturel trouvé : sur un dessus
  en arc ou à arête vive, la boule touche la face là où elle est moins courbe
  qu'elle (un dôme circulaire donne toujours r ≤ son rayon ; une arête vive est
  touchée à ses épaules). Le calcul (`bending`) est testé seul ; le refus,
  non.
- **(b) sur un côté** : le code est celui de M ; le test le prend sur M.
- **Le départ à la couture d'un côté B-spline fermé** (question 14) : non codé.

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

1. **La marche de `fillet_law`** (§ 3.2), garde à 4 fois la vitesse du guide
   sur le pas divisé, pour tous les congés tritangents. *Recommandé* : elle
   corrige l'épingle et l'ailette ; le coût passe de 2 à 6 ou 9 ms sur une
   nervure droite.
2. **Les refus de repli par la courbure aux contacts** (§ 3.4), en travers et
   le long, au lieu du contrôle de l'offset de toute la face. *Recommandé* :
   l'offset voit des replis loin du congé (le fond d'un U) et manque le nez d'une
   ailette. Le refus (a) n'a pas de cas de bout en bout (§ 6).
3. **Les surfaces exactes contrôlées** (§ 3.3) : 17 stations pour reconnaître,
   129 pour contrôler, les contacts interpolés à 10⁻⁹ en paramètre, la ligne des
   centres à `tol`. *Recommandé.*
4. **Le départ à la couture d'un côté fermé** (question 14) non codé : la boucle
   part toujours du sommet du guide. *Recommandé* pour cette PR : aucun cas de
   test n'en a besoin ; à reprendre avec un cas réel.
5. **Les trois tests de la PR 3 qui changent d'attente** (§ 4) : ce qu'ils
   refusaient est construit. *Recommandé.*
6. **La comparaison du disque à OCCT à 10⁻⁶** (§ 5), l'écart de relecture
   d'OCCT non élucidé. *Recommandé* : l'analytique de gbs est à 5·10⁻¹⁰.
