# Noyau BREP natif — conception : congés qui engloutissent une face

| | |
|---|---|
| Statut | Questions du § 8 validées (23 nov. 2026) ; PR 1 : [brep_swallowed_faces_pr01.md](brep_swallowed_faces_pr01.md) |
| Périmètre | une arête congée ou chanfreinée dont le traitement, continué au-delà d'un sommet, passe **sous** une petite face voisine sans la toucher, et rejoint une autre arête entre les mêmes deux faces : la petite face, les arêtes et les sommets entre elles disparaissent, et le traitement court d'une arête à l'autre (**engloutissement**) ; diagnostic du refus d'origine, détection, contrôles, géométrie, topologie, tests |
| Dépend de | palier 4 : géométrie et topologie des congés et chanfreins (`fillet.h`, `blend.h`), chaînes fermées ; palier 5 : lois ; palier 6 : onglets (`mitre`), fins libres (`reaches`) ; palier 2 : intersection de surfaces (ligne des centres) |
| Hors périmètre | le **débordement** (la boule roule par-dessus une arête vive, son contact passe d'une face à une autre, § 5) ; les congés de faces (*face blend*) entre deux faces sans arête commune ; un mode imposé par l'utilisateur (garder l'arête ou la surface) |

Ce document décrit *quoi* construire et *pourquoi*, comme les documents des
paliers 4 à 7. Il part d'un refus constaté le 23 nov. 2026 (§ 1). Chaque PR aura
sa note d'architecture et sa note d'illustration.

**Provenance.** Le comportement visé est celui de CATIA et de SolidWorks, la
cible de qualité : une arête ou une face prise entièrement dans le volume ôté
par un congé disparaît sans erreur. OpenCascade ne sert que de comparaison en
boîte noire.

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Cause du refus | **géométrique**, pas numérique : les sommets en cause sont de vrais angles (34°), et la boule du congé passe sous la petite face (§ 1) | élargir la tolérance de tangence (`tangent_at`) ou rendre tangents les bouts des arêtes du booléen : sans effet ici, et faux en général |
| Trois régimes | à une extrémité, selon que la boule continuée **touche** la face au-delà (onglet ou fin libre, déjà fait), la **coupe** en partie (débordement, refusé avec un message propre, § 5) ou **passe dessous** (engloutissement, ce document) | un seul traitement général par booléens locaux |
| Détection | une **poche** : entre deux extrémités v et v′ d'arêtes traitées bordées par les **mêmes faces** F₁ et F₂, les bords de F₁ et de F₂ vont de v à v′ par des arêtes non traitées (ou données mais englouties) qui enferment des faces ; la ligne des centres, continuée depuis v, rejoint celle de v′ sur la même branche | chercher les faces englouties par un test de volume sur tout le solide |
| Contrôle | engloutissement accepté seulement si, sur le pont, **les contacts restent dans F₁ et F₂** et si **la surface du traitement ne rencontre aucune face de la poche** (intersection de surfaces), tous les sommets de la poche étant du côté ôté | se fier au seul échec de l'onglet |
| Géométrie | **un seul traitement** le long de la chaîne arête, pont, arête : la ligne des centres est l'intersection des surfaces décalées de F₁ et F₂, continuée à travers le pont (même calcul qu'aujourd'hui, prolongé jusqu'à l'autre arête) ; fermé si la chaîne se referme | deux traitements prolongés et recollés au milieu du pont |
| Topologie | les faces, arêtes et sommets de la poche **effacés** ; dans F₁ et F₂, le parcours de co-arêtes de v à v′ par la poche remplacé par le contact, comme une arête traitée | garder une face dégénérée |
| Arête donnée et engloutie | **acceptée** et sans effet (son traitement disparaît avec elle), comme dans CATIA et SolidWorks | refus (« blend the pocket edges separately ») |
| Arête non donnée et engloutie | **acceptée** aussi : c'est la seule forme valide pour le rayon demandé | refus : l'utilisateur doit la donner |
| Lois et chanfreins | la loi le long de la chaîne, pont compris (abscisse curviligne de la ligne des centres) ; chanfreins mesurés depuis la courbe d'intersection de F₁ et F₂ continuée | congés de rayon constant seulement |

---

## 1. Constat

### 1.1 Le cas

Un cylindre de rayon 1 et de hauteur 2, coupé par le plan incliné de 30° autour
de Y qui passe par (0, 0, 1,5) : z = 1,5 − x tan 30°. Ce plan dépasse z = 2 pour
x < −cos 30° ≈ −0,866 (z = 2,077 en x = −1). Il coupe donc **aussi le dessus** :
la face inclinée est bordée par **deux arcs d'ellipse** (joints à la couture du
cylindre) et par un **segment** (avec le dessus). Il reste du dessus un
**croissant** entre ce segment et le cercle du haut, de largeur 0,134.

`fillet(m, s, bord, 0.14)` sur les trois arêtes du bord est refusé : « the
contacts of two blends at a vertex do not cross on the face between them within
the blends continued ». Sur les deux arcs seulement, il est refusé à leur fin
libre.

### 1.2 Mesures

Directions des arêtes en partant de chaque sommet du bord :

| Sommet | Arêtes | cos | Écart à la tangence |
|---|---|---|---|
| couture (x = 1) | arc, arc | −1 (exact) | 0 |
| (−0,866 ; ±0,5 ; 2) | arc, segment | −0,832 | 0,588 rad (34°) |

L'hypothèse d'une tangence manquée de peu (le seuil de 1e-6 sur le cosinus dans
`tangent_at`) est donc **écartée** : les deux arcs sont joints exactement, et les
deux autres sommets sont de vrais angles, où le cercle du dessus est la
troisième arête. Le congé les traite à juste titre en onglet.

Les angles entre faces le long des arêtes :

- **arcs** (cylindre, plan incliné) : de 120° du côté x > 0 à **60°** du côté
  x < 0 ; la largeur du contact sur le plan, r / tan(α/2), y monte à **1,73 r** ;
- **segment** (plan incliné, dessus) : 150° ; contact à **0,27 r**.

Balayage du rayon, les trois arêtes données :

| r | Résultat aujourd'hui | Ce qui se passe |
|---|---|---|
| 0,02 ; 0,03 | construit, `check` valide, 3 faces de congé | onglets aux deux angles |
| 0,04 à 0,105 | refusé : « a contact of a blend at a vertex does not meet the face beyond the sharp edge » | la boule de l'arc **coupe** le croissant, mais son contact sur le cylindre reste sous z = 2 : elle roule par-dessus le segment (**débordement**, § 5) |
| 0,11 ; 0,14 | refusé : les contacts sur le plan incliné ne se croisent pas | la boule de l'arc **passe sous** le croissant (**engloutissement**) |

Le seuil entre les deux derniers régimes se calcule : en x = −1, le centre de la
boule de l'arc est à z = 2,077 − 1,732 r, son sommet à 2,077 − 0,732 r, sous le
dessus dès que **r > 0,1057**. Le refus mesuré tombe entre 0,105 et 0,11. C'est
aussi le seuil où les contacts sur le plan incliné cessent de se croiser : la
boule tangente au dessus sépare l'onglet de l'engloutissement.

### 1.3 Ce qu'il faut construire pour r = 0,14

La boule qui roule entre le cylindre et le plan incliné le long de l'ellipse
**entière** reste, au-dessus du croissant, sous z = 2, et ses contacts restent
sur les faces réelles (sur le plan incliné en x ≥ −0,79, à droite du segment en
x = −0,866 ; sur le cylindre sous z = 1,83). Le croissant, le segment, l'arc du
cercle du dessus et les deux sommets d'angle sont tous dans le volume ôté. Le
résultat attendu est un **congé fermé** autour de l'ellipse, sans croissant :
celui qu'on obtiendrait si le plan ne coupait pas le dessus. Le congé demandé
sur le segment n'a plus d'objet.

---

## 2. Objectifs et non-objectifs

### 2.1 Objectifs

- Construire les traitements qui engloutissent une **poche** (§ 3.1) : une ou
  plusieurs petites faces, entre deux arêtes traitées bordées par les mêmes
  deux faces, que la boule ou la règle du chanfrein enjambe sans les toucher.
- Cas couverts : le bord du cylindre coupé (§ 1), un méplat peu profond au
  milieu d'une arête de boîte (§ 6), une arête unique dont les deux sommets
  bordent la même poche (le traitement devient fermé).
- Convexe (matière ôtée, la poche dans le volume ôté) et concave (matière
  ajoutée, la poche dans le volume ajouté), congés à rayon constant ou en loi,
  chanfreins.
- Arêtes de la poche données ou non : acceptées, sans effet.
- Python inchangé (mêmes fonctions), galerie, comparaison à OpenCascade.

### 2.2 Non-objectifs

- Le **débordement** (§ 5) : refus avec un message qui le nomme.
- Les congés entre deux faces sans arête commune (un congé de top et de côté qui
  engloutirait un chanfrein entre eux) : c'est un congé de faces.
- Un choix de l'utilisateur entre garder l'arête et garder la surface (le
  *overflow type* de SolidWorks, les *edges to keep* de CATIA).

---

## 3. Détection

### 3.1 Poche

À une extrémité v d'un traitement g, de faces F₁ et F₂, aujourd'hui traitée en
fin libre ou en onglet : on suit le bord de F₁ depuis v, du côté qui ne porte
pas l'arête de g, et de même celui de F₂. S'ils arrivent tous deux à un même
sommet v′ extrémité d'une arête traitée g′ bordée par F₁ et F₂ (g′ = g possible :
les deux sommets d'une même arête), sans passer par une arête traitée qui ne
soit pas dans la poche, les deux parcours ferment une **poche** :

- ses **arêtes** : celles des deux parcours, et celles des faces qu'ils
  enferment ;
- ses **faces** : celles bordées seulement par des arêtes de la poche (le
  croissant ; les trois faces d'un méplat) ;
- ses **sommets**, v et v′ compris.

La poche doit être un disque : ses faces se suivent par leurs arêtes communes,
aucune n'a de trou, aucune ne touche le reste du solide hors des deux parcours.
Sinon, pas d'engloutissement : le traitement actuel s'applique (et refuse s'il
le faut).

### 3.2 Pont

La ligne des centres de g (intersection des surfaces décalées de F₁ et F₂) est
continuée au-delà de v par la marche existante, jusqu'au point le plus proche de
v′, sur la **même branche** que celle de g′ (à tol près, et même sens). Le
morceau entre v et v′ est le **pont**. S'il n'existe pas (la ligne s'arrête,
change de branche, revient sur elle-même), pas d'engloutissement.

### 3.3 Contrôles

L'engloutissement n'est accepté que si **tout** ce qui suit tient (sinon le
traitement actuel s'applique) :

1. **Contacts sur les faces.** Sur le pont, chaque contact est dans sa face (F₁
   ou F₂), à l'intérieur de son bord une fois la poche retirée : il ne coupe
   aucune arête des parcours de la poche. Un contact qui en sort, c'est le
   débordement.
2. **Rien de la poche hors du volume ôté.** La surface du traitement, sur le
   pont, ne rencontre aucune face de la poche (intersection de surfaces, bornée
   aux faces) ; chaque sommet de la poche est du côté ôté (dans la section du
   traitement à la station la plus proche : hors de la boule, entre les deux
   normales aux contacts, du côté de l'arête ; pour un chanfrein, du côté de
   l'arête par rapport à la règle). Le sens est inversé pour une arête concave.
3. **Même convexité** de g et de g′, et même nature (congé ou chanfrein, même
   section).
4. **Lois compatibles** aux deux bouts (§ 4.3).

Le contrôle 2 rend les régimes exclusifs : quand la boule touche une face de la
poche, l'onglet ou la fin libre existent (cas r ≤ 0,03 du § 1.2), ou c'est le
débordement.

### 3.4 Ordre

La détection passe **avant** les fins libres et les onglets, à chaque extrémité
non jointe : la chaîne est d'abord prolongée à travers les poches acceptées,
puis les extrémités restantes sont traitées comme aujourd'hui. Les chaînes se
recalculent après (`law_chains` voit le pont comme un maillon).

---

## 4. Construction

### 4.1 Géométrie

Un seul traitement le long de la chaîne g, pont, g′ (et les maillons suivants) :
la ligne des centres, les contacts et la surface calculés comme aujourd'hui
(`fillet_geometry`), sur une **arête virtuelle** : la courbe d'intersection de
F₁ et F₂ continuée de v à v′, qui sert de paramètre et de référence au chanfrein.
Pas de recollement au milieu du pont : la surface est d'un seul tenant, et
exacte quand elle l'est (tore ou cylindre sur une ligne des centres circulaire
ou droite : le méplat du § 6).

Si la chaîne se referme (le bord du cylindre), le traitement est **fermé**, sans
extrémité, comme celui d'une arête fermée.

### 4.2 Topologie

- Faces, arêtes et sommets de la poche **effacés**.
- Dans F₁ et dans F₂, le parcours de co-arêtes de v à v′ par la poche est
  remplacé, avec les arêtes traitées voisines, par le contact (le mécanisme des
  chaînes : un seul parcours remplacé de bout en bout ; une chaîne fermée
  remplace un fil entier).
- Une arête traitée de la poche (le segment du § 1) n'a pas de traitement
  propre : son congé est retiré de la liste avant la topologie.
- Une couture traversée par le contact est coupée là, comme aujourd'hui.
- Le solide est vérifié (`check`) ; l'échec rend le modèle intact.

### 4.3 Lois

La loi du palier 5 est définie sur les arêtes ; le pont n'en a pas. Proposé : la
loi continue à travers le pont en abscisse curviligne de la ligne des centres,
comme le long d'une chaîne (la loi lisse, C1, par défaut, l'est par
construction) ; une loi dont les valeurs aux deux bords du pont diffèrent varie
le long du pont selon son mode (linéaire ou lisse). Un pli de loi au milieu
d'une poche n'a pas de sens : refusé.

### 4.4 Refus

Ajoutés (`NonGeneric`) : un pont dont la marche échoue ; une poche acceptée par
les contrôles, mais dont le solide obtenu n'est pas valide. Le **débordement**
(contrôle 1 en échec) reçoit son propre message : « a blend that rolls over a
sharp edge next to it (overflow, not handled) », à la place des deux messages
actuels sur les contacts qui ne se croisent pas ou n'atteignent pas la face.

---

## 5. Débordement (hors périmètre)

Pour 0,035 < r < 0,106 dans le cas du § 1, la boule de l'arc coupe le croissant
sans l'engloutir : sur une partie du pont, elle touche le cylindre et le
**dessus**, plus le plan incliné. Le congé passe d'une paire de faces à une autre
en roulant par-dessus le segment, avec un point où la boule touche trois faces.
SolidWorks le règle par son *overflow type* (garder l'arête, garder la surface),
CATIA par les *edges to keep*. C'est un chantier à part (changement de paire de
faces le long d'un congé, raccord à trois faces). Il est seulement **détecté**
et refusé avec un message clair (§ 4.4).

---

## 6. Validation

| Test | Cas | Attendu |
|---|---|---|
| bord du cylindre (§ 1) | r = 0,11 ; 0,14 ; 0,2, les trois arêtes, puis les deux arcs seuls | congé fermé, plus de croissant ; `check` valide ; même solide avec ou sans le segment donné ; **même solide que le congé du cylindre de hauteur 3 coupé par le même plan** (bord elliptique d'une seule arête fermée, sans croissant : le croissant est dans le volume ôté), volumes égaux à approx_tol près |
| régime intermédiaire | r = 0,05 ; 0,08 | refus « overflow » |
| régime bas | r = 0,02 ; 0,03 | inchangé (test `fillet_tilted_cut_rim`) |
| méplat | boîte dont une arête porte, au milieu, un méplat à 45° de profondeur p ; congé r de l'arête (deux morceaux) | pour p < r(√2 − 1) (le méplat sous le congé, mesuré le long de la bissectrice) : **même solide que le congé de la boîte sans méplat** (surface cylindrique exacte, volume exact à 1e-12) ; pour p plus grand : onglet ou refus, jamais d'engloutissement |
| arête unique | une arête ouverte dont les deux sommets bordent la même poche | traitement fermé |
| concave | la poche dans un angle rentrant (matière ajoutée) | idem, volume croissant |
| chanfrein | méplat avec un chanfrein | même solide que le chanfrein de la boîte sans méplat |
| loi | méplat, rayon lisse de r₁ à r₂ le long de l'arête | loi continue à travers le pont |
| Python | le bord du cylindre | même résultat qu'en C++ |
| OpenCascade | les mêmes, `BRepFilletAPI_MakeFillet` | comparaison en boîte noire (volume, réussite) |

Chaque cas a son rendu dans la note d'illustration.

---

## 7. Plan des PR

1. **Engloutissement des congés** : poche, pont, contrôles, géométrie et
   topologie ; congés de rayon constant ou en loi, convexes et concaves ; refus
   « overflow » ; tests du § 6 sauf chanfreins ; le test
   `fillet_tilted_cut_rim` (régime bas) part avec elle.
2. **Chanfreins et Python** : chanfreins sur l'arête virtuelle, galerie,
   comparaison à OpenCascade.

Une seule PR si la première reste raisonnable.

---

## 8. Questions à valider

1. **Diagnostic** (§ 1) : la tolérance de tangence reste telle quelle ; le refus
   vient de la géométrie.
2. **Arête non donnée engloutie** (§ 0) : accepter que le congé des deux arcs
   seuls fasse disparaître le segment, qui n'était pas demandé ? (Proposé : oui,
   c'est la seule forme valide ; CATIA et SolidWorks font de même.)
3. **Arête donnée engloutie** : acceptée sans avertissement ? `BuildResult` n'a
   pas de canal d'avertissement ; en ajouter un serait un chantier commun à toutes
   les opérations. (Proposé : acceptée, documentée.)
4. **Ordre** (§ 3.4) : l'engloutissement avant les onglets et les fins libres,
   le contrôle 2 garantissant qu'ils ne se recouvrent pas.
5. **Contrôle 2** par intersection de la surface du traitement avec les faces de
   la poche, plutôt que par échantillonnage de points.
6. **Lois** (§ 4.3) : la loi continue à travers le pont plutôt qu'un refus des
   lois dans un premier temps.
7. **Débordement** (§ 5) : hors périmètre, refus avec son propre message.
8. **Place dans la feuille de route** : après la PR 2 du palier 7 (en attente),
   comme palier 8, ou comme écart au palier 7 ?

### Sources

- Vida, Martin et Várady, « A survey of blending methods that use parametric
  surfaces » (Computer-Aided Design 26(5), 1994) : interférences entre un
  congé et les faces voisines.
- Rossignac et Requicha, « Constant-radius blending in solid modelling »
  (Computers in Mechanical Engineering, 1984) : volume ôté par la boule.
- SolidWorks, aide du congé : *Overflow type* (*Default*, *Keep edge*, *Keep
  surface*). CATIA, *Edge Fillet* : *Edge(s) to keep*, *Trim ribbons*.
