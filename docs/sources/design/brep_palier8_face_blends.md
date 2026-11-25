# Noyau BREP natif — conception du palier 8 : congés de faces et congés tritangents

| | |
|---|---|
| Statut | Questions 1 à 11 du § 8 validées (23 nov. 2026) ; PR 1 : [brep8_pr01_face_blends.md](brep8_pr01_face_blends.md) (fusionnée, #246) ; ajout du 24 nov. 2026 : bande partielle (§ 2.6), question 12 validée le même jour ; PR 2 : [brep8_pr02_partial_strips_laws.md](brep8_pr02_partial_strips_laws.md) (fusionnée, #249) ; PR 3 : [brep8_pr03_full_round.md](brep8_pr03_full_round.md) (congé tritangent, fusionnée, #251) ; bilan du palier au § 9 ; ajout du 25 nov. 2026 pour la PR 4 : face du milieu qui fait le tour (§ 3.5), bout arrondi (§ 3.6), surfaces exactes tirées des stations (§ 3.7), questions 13 à 18 à valider, revues après la revue de conception |
| Périmètre | **Congé entre deux faces** (*face-face fillet*) : deux faces, ou deux ensembles de faces jointes en tangence, sans arête commune ; la boule roule sur les deux, les faces entre elles sont englouties ; rayon constant ou en loi, sections du palier 7. **Congé tritangent** (*full round*) : la boule touche trois faces, la face du milieu disparaît, le rayon suit la géométrie. Pour chacun : géométrie, topologie, contrôles, Python, galerie, comparaison à OpenCascade. |
| Dépend de | Palier 4 : congés (`fillet.h`, `blend.h`), fins libres. Palier 5 : lois, stations à rayon inconnu. Palier 6 : onglets, raccords de sommet. Palier 7 : sections, `BlendLaw`. Poches englouties ([brep_swallowed_faces.md](brep_swallowed_faces.md)) : arête de remplacement (`Bridge`), copie remplie (`fill`), contrôles d'engloutissement. Palier 2 : intersection de surfaces. |
| Hors périmètre | Le **débordement** (la boule passe d'une paire de faces à une autre par-dessus une arête vive). Deux faces que **rien ne relie** (deux nervures sur une plaque : le congé aurait des bouts ouverts). Le congé tritangent à section G2 ou conique. Les raccords de sommet G2. |

Ce document décrit *quoi* construire et *pourquoi*, comme les documents des
paliers 4 à 7. Chaque PR a sa note d'architecture `brep8_prNN_*.md` et sa note
d'illustration. C'est le dernier point ouvert de la feuille de route du
palier 4 (§ 9, point 3) que chaque bilan depuis a reporté : « congés de
faces ».

**Provenance.** Les algorithmes viennent de la littérature et du comportement de
CATIA (« Face-Face Fillet », « Tritangent Fillet ») et de SolidWorks
(« Face fillet », « Full round fillet »), la cible de qualité ; chaque note
cite ses sources. OpenCascade ne sert que de comparaison en boîte noire dans
les tests : il n'a ni congé de faces ni congé tritangent dans son API
topologique, la comparaison passe par des cas équivalents à un congé d'arête
(§ 6).

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Congé de faces | **une arête virtuelle** : sur une copie du solide, les faces entre F₁ et F₂ (la **bande**) sont ôtées et F₁, F₂ prolongées jusqu'à leur intersection ; le congé de cette arête, avec toute la machinerie des paliers 4 à 7, est le congé de faces. C'est la **généralisation de la poche** de [brep_swallowed_faces.md](brep_swallowed_faces.md) : là, la bande s'arrêtait entre deux arêtes données ; ici, aux faces de bout | une boule roulée sans arête (ligne des centres cheminée depuis une graine, topologie écrite à part) ; un booléen entre le solide et le volume du congé (contacts tangents tout le long : le cas le plus difficile des booléens) |
| Bande | les faces **entre** F₁ et F₂ : composante connexe des faces voisines de F₁ et de F₂ que l'on peut joindre de l'une à l'autre sans les traverser, bordée à ses bouts par des **faces de bout** qui touchent F₁ et F₂ ; une seule bande par appel | toutes les faces du côté ôté (non bornées : le reste du solide) |
| Bouts | les faces de bout, prolongées, coupent l'arête virtuelle : ses sommets sont les points où se rencontrent F₁, F₂ et la face de bout prolongées (Gauss-Newton sur trois surfaces) ; le congé y a une **fin libre** du palier 4, ou un onglet si la face de bout est elle-même un congé | prolonger le congé jusqu'aux bords de F₁ et F₂ (bouts ouverts, solide non fermé) |
| Côté de la boule | **celui de la matière**, comme pour une arête : convexe, le congé ôte ; concave (gorge, rainure), il ajoute | un côté choisi par l'utilisateur |
| Contrôle | celui des poches : accepté seulement si **les contacts restent dans F₁ et F₂** et si **la surface du congé ne rencontre aucune face de la bande**, ses sommets du côté ôté (ou ajouté) ; sinon refus, avec le message du débordement si la boule touche la bande | se fier à la validité du solide final |
| Ensembles de faces | chaque côté peut être **un ensemble de faces jointes en tangence** (G1 le long de leurs arêtes communes) : l'arête virtuelle est alors une **chaîne tangente** (une arête par paire de faces), traitée par les chaînes du palier 4 | une face par côté seulement |
| Faces déjà voisines | F₁ et F₂ qui ont des arêtes communes : le congé de ces arêtes (rien de nouveau), et une bande vide | refus |
| Bande partielle | F₁ et F₂ qui se touchent sur une partie de leur longueur et sont séparées ailleurs par une bande contiguë à leur arête commune (une face en pointe, un chanfrein arrêté) : l'arête virtuelle **prolonge l'arête commune** sous la bande, les deux sont fondues en une arête sur la copie, puis le congé et le contrôle du § 2.4 (§ 2.6) | le congé de l'arête commune seule, refusé à son bout (« a blend wider than a face at an end of its edge ») |
| Lois | la `BlendLaw` du palier 7, ses points donnés **par position** : un point de l'espace, projeté sur l'arête virtuelle (la plus proche), et une valeur ; aussi admis pour les congés d'arête | des points sur les arêtes de la copie, que l'utilisateur ne voit pas |
| Sections | toutes celles du palier 7 (cercle, conique, corde constante, G2), lues sur la copie comme pour les poches (`section_on_copies`) | cercle seulement |
| Congé tritangent | la boule tangente à **trois faces** F₁, M, F₂ ; la station résout les trois contacts, le centre et le **rayon** (Gauss-Newton : 10 inconnues, 9 équations de contact, le plan de la station) ; la section est l'arc de P₁ à P₂ **par** le contact P_M ; M disparaît, le congé prend sa place | une arête virtuelle entre F₁ et F₂ (il n'y en a pas quand les flancs sont parallèles : le cas le plus courant, une nervure sans dépouille) |
| Plan de la station tritangente | **normal à l'arête M ∩ F₁** à son abscisse s (la station du palier 5) ; le paramètre du congé est celui de cette arête | normal à une ligne des centres encore inconnue |
| Topologie tritangente | M ôtée ; F₁ et F₂ bornées par les contacts ; aux bouts de M, les faces de bout coupent le congé (deux fins libres qui se rejoignent sur la même face) | un raccord de sommet à chaque bout |
| Surfaces exactes | congé de faces : celles du congé d'arête sur la copie (cylindre, tore, extrusions de degré (2 ou 4, 1)) ; tritangent : **cylindre** entre deux plans et un troisième (nervure droite, avec ou sans dépouille) ; approchées ailleurs (peau de Hermite sous `approx_tol` / `angular_tol`) | toujours approchées |
| API | `fillet_faces(m, s, faces1, faces2, r ou loi, options)` ; `fillet_full_round(m, s, faces1, middle, faces2, options)` ; Python à l'identique | des options cachées dans `fillet` |

---

## 1. Objectifs et non-objectifs

### 1.1 Objectifs

- **Congé de faces**, comme le « Face-Face Fillet » de CATIA et le « Face
  fillet » de SolidWorks : un congé entre deux faces qui ne se touchent pas,
  par-dessus les faces qui les séparent (un chanfrein, une petite marche, le fond
  d'une gorge). Les cas courants :
  - arrondir par-dessus un chanfrein ou une petite face plate (reprise de modèle
    importé) ;
  - combler une gorge ou une rainure à flancs en V : le congé ajoute la matière
    sous la boule, le fond disparaît ;
  - un pied de bossage chanfreiné repris en congé ;
  - une petite face entre les deux sur une partie seulement de leur longueur :
    une face triangulaire qui finit en pointe sur leur arête commune, un
    chanfrein arrêté au milieu (§ 2.6). C'est un manque d'OpenCascade : un
    congé n'y supprime pas de face.
- **Ensembles de faces** de chaque côté : un flanc fait de plusieurs carreaux
  B-spline jointes en tangence (intrados d'aube, côté d'un bossage en stade).
- Rayon constant, **en loi**, et les **sections** du palier 7.
- **Congé tritangent**, comme le « Tritangent Fillet » de CATIA et le « Full
  round fillet » de SolidWorks : la face du milieu remplacée par un congé
  tangent à elle et à ses deux voisines. Les cas courants :
  - le dessus d'une nervure arrondi en demi-cylindre ;
  - le **bord de fuite** (ou d'attaque) d'une aube : la face étroite du bord
    remplacée par l'arrondi tangent à l'intrados, à l'extrados et au bord. C'est
    le cas qui motive le plus ce palier pour les aubages.
- Python, galerie, comparaison à OCCT par des cas équivalents.

### 1.2 Non-objectifs

- **Débordement** ([brep_swallowed_faces.md](brep_swallowed_faces.md), § 5) : la
  boule qui passe d'une paire de faces à une autre. La station tritangente de ce
  palier en est une brique (le point où la boule touche trois faces), le reste
  est un autre chantier. Il reste détecté et refusé.
- **Deux faces que rien ne relie** : deux nervures en vis-à-vis sur une plaque.
  Le congé entre leurs flancs aurait des bouts ouverts, que CATIA ferme par des
  faces planes et SolidWorks refuse. Refusé avec un message clair.
- **Plusieurs bandes** dans un même appel (F₁ et F₂ qui se font face en deux
  endroits séparés) : un appel par bande.
- **Section G2 ou conique** du congé tritangent : l'arc de cercle seulement.
- **Raccords de sommet G2**, toujours hors périmètre (palier 7).

---

## 2. Congé de faces

### 2.1 Bande et faces de bout

Données : deux ensembles de faces F₁ et F₂ d'un même solide, chacun connexe par
des arêtes de tangence (G1 sous `angular_tol`), sans face commune.

- **Bande** : on part des faces voisines de F₁ qui ne sont pas dans F₂, et l'on
  avance de face en face sans traverser F₁ ni F₂. Une composante qui touche
  aussi F₂ est une bande candidate. Ses **faces de bout** sont celles qui
  touchent F₁, F₂ et la bande, et que le chemin d'une face de F₁ à une face de
  F₂ contourne : aux deux extrémités d'une bande ouverte, aucune pour une bande
  qui fait le tour (pied d'un bossage).
- La bande peut être **vide** : F₁ et F₂ ont des arêtes communes ; le congé de
  faces est alors le congé de ces arêtes (palier 4 à 7).
- Plusieurs candidates : la **plus petite** dont le pont passe le contrôle du
  § 2.4 (comme les deux candidats d'une poche, l'autre étant le reste du
  solide) ; aucune : refus « faces that no strip of faces joins ».

### 2.2 Arête virtuelle et copie

C'est le `Bridge` des poches, généralisé :

1. **Courbe** : l'intersection des surfaces de F₁ et F₂ prolongées (palier 2),
   la branche qui passe sous la bande, du côté opposé à la boule. Avec des
   ensembles de faces, une courbe par paire de faces (F₁ⁱ, F₂ʲ) qui se font
   face, jointes en chaîne tangente aux arêtes internes des ensembles.
2. **Sommets** : aux bouts, le point commun de F₁, F₂ et de la face de bout E,
   prolongées (Gauss-Newton sur trois surfaces depuis le point de la courbe le
   plus proche de E). Pour une bande qui fait le tour : une courbe fermée.
3. **Copie** : la bande ôtée, F₁ et F₂ bornées par l'arête virtuelle (leurs
   co-arêtes le long de la bande remplacées par elle, comme les `runs` d'une
   poche), chaque face de bout recousue (ses co-arêtes avec la bande remplacées
   par le prolongement de ses arêtes avec F₁ et F₂ jusqu'au nouveau sommet,
   ou une seule co-arête si E ne touchait la bande que par une arête).

La copie est un solide valide (`check`) : un modèle à arête vive, celui qu'on
aurait obtenu sans chanfrein, sans gorge ou sans marche.

### 2.3 Le congé

Le congé de l'arête virtuelle sur la copie, avec les mêmes options : rayon ou
loi, section, retraits, tolérances. Toute la topologie des paliers 4 à 7 s'y
applique : fins libres coupées par les faces de bout, chaînes, onglets si la
face de bout est un congé déjà fait, raccord de sommet si elle est traitée dans
le même appel (hors première version).

Le **côté** suit la copie : à une arête virtuelle convexe le congé ôte de la
matière (par-dessus un chanfrein), à une arête concave il en ajoute (comble une
gorge en V).

### 2.4 Contrôle d'engloutissement

Le congé de la copie est celui du solide donné si la bande est tout entière
dans le volume que le congé ôte (ou ajoute) entre sa surface et l'arête
virtuelle. Les contrôles des poches (§ 3.3 de
[brep_swallowed_faces.md](brep_swallowed_faces.md)) :

- les **contacts** restent dans F₁ et F₂ (pas dans la partie prolongée) ;
- la **surface du congé ne rencontre aucune face de la bande** (intersection
  de surfaces, bornée au congé) ;
- tous les **sommets de la bande** sont du côté ôté (ou ajouté).

Échec : refus, avec le message du débordement si la boule touche la bande
(« the blend spills over the faces between: a larger radius, or blend their
edges »), ou du rayon si un contact sort d'une face (« the blend is wider than
the faces »).

### 2.5 Lois par position

Les points d'une `BlendLaw` sont aujourd'hui (arête, fraction de sa longueur).
L'arête virtuelle n'existe que sur la copie. Un point de loi peut donc aussi
être donné **par position** : un point de l'espace et une valeur, projeté sur
l'arête virtuelle (le pied le plus proche, sur la chaîne). Les deux formes se
mêlent dans une même loi ; la forme par position sert aussi aux congés d'arête
(un point posé sur un modèle, comme dans CATIA). La variation (linéaire, lisse)
et l'abscisse curviligne sont celles du palier 5.

### 2.6 Bande partielle : des faces qui se touchent déjà sur une partie

**Cas.** F₁ et F₂ ont une arête commune e (ou une chaîne) sur une partie de
leur longueur ; ailleurs, une bande les sépare, **contiguë à e** :

- une **face en pointe** : triangulaire entre F₁ et F₂, large à un bout,
  finissant en pointe sur un sommet de e (la pointe est le sommet commun de e,
  de la face et de ses arêtes avec F₁ et F₂) ;
- un **chanfrein arrêté** : la bande finit à un sommet de e par une face
  d'extrémité (le triangle qui ferme le chanfrein), qui touche F₁, F₂ et la
  bande ;
- une bande **au milieu** de l'arête : deux arêtes communes e₁ et e₂, la bande
  entre elles.

Aujourd'hui, `fillet_faces` y fait le congé de e seule, et `fillet` de e est
refusé à son bout : le congé y est plus large que la face qui le borne. Les
deux cas ont été construits sur master (un U, paroi et fond, face en pointe ou
chanfrein arrêté à mi-longueur) : refusés pour r = 0,2 à 0,5, alors que la
bande tient tout entière sous le congé pour r ≥ 0,3.

**Bande.** La recherche du § 2.1, avec une règle de plus : une face qui touche
F₁ et F₂ et dont les arêtes avec F₁ et F₂ **finissent sur un sommet de e** est
dans la bande, pas une face de bout (elle est sous le prolongement de e, la
face de bout d'une boîte ne l'est pas). Une face en pointe est dans la bande
par la même règle. Les faces de bout restent celles qui bornent la bande de
l'autre côté (le bout du U).

**Arête virtuelle.** L'intersection de F₁ et F₂ prolongées, comme au § 2.2,
partie du sommet de e (la pointe, ou le sommet où finit le chanfrein) jusqu'à
la face de bout. C'est **le prolongement de la courbe de e** : la même paire de
surfaces. Sur la copie, e et l'arête virtuelle sont **fondues en une arête**
(l'intersection recalculée sur leur réunion, le sommet entre elles ôté), et le
congé est celui d'une seule arête : pas de station de chaîne au milieu, la
surface exacte (cylindre, tore) d'un seul tenant. Bande au milieu : e₁, l'arête
virtuelle et e₂ fondues.

**Congé et contrôle.** Le congé de l'arête fondue sur la copie (§ 2.3), puis le
contrôle du § 2.4 sur la bande seule : contacts dans F₁ et F₂, aucune face de
la bande coupée par le congé, ses sommets du côté ôté ou ajouté. Un rayon trop
petit (le coin du triangle dépasse de l'arc) donne le refus du débordement.

**API.** Aucune nouvelle : `fillet_faces(m, s, F₁, F₂, r)` prend la bande
partielle quand F₁ et F₂ se touchent déjà et qu'une bande est contiguë à leur
arête commune ; sans bande, le congé des arêtes communes comme aujourd'hui. Le
refus de `fillet` d'une arête à son bout, quand une bande y est contiguë, nomme
`fillet_faces`.

**Refus.** Une bande contiguë à e mais dont l'intersection prolongée ne
prolonge pas e (F₁ et F₂ qui changent de surface à la pointe : une chaîne
tangente seulement, traitée comme les ensembles de faces) ; une bande qui ne
tient pas sous le congé (débordement).

---

## 3. Congé tritangent

### 3.1 Station

Données : F₁, M, F₂ (chacun un ensemble de faces jointes en tangence), M voisin
de F₁ et de F₂. Pas d'arête virtuelle : entre deux flancs parallèles (nervure
sans dépouille), F₁ et F₂ ne se coupent pas.

La station se pose sur le **guide** g, l'arête M ∩ F₁ (une chaîne avec des
ensembles), au point g(s) : dans le plan normal à g en g(s), la boule de
centre c et de rayon r tangente aux trois faces. Inconnues : (u, v) des trois
contacts, c, r (10) ; équations : c = Sₖ(uₖ, vₖ) + r Nₖ pour k = 1, M, 2 (9),
c dans le plan (1). Gauss-Newton, comme la station à rayon inconnu de la corde
constante (palier 7, § 4.3), depuis la station voisine ; à la première, depuis
le cercle inscrit dans les tangentes des trois faces au point de g.

Le côté : les normales Nₖ vers l'intérieur de la matière (la boule est dans le
solide, le congé ôte M et ce qui dépasse de la boule).

### 3.2 Section et surface

La section est l'**arc de la boule de P₁ à P₂ qui passe par P_M** (plus de 180°
quand les flancs se rapprochent vers le haut : une B-spline rationnelle en
plusieurs arcs, comme les arcs larges du palier 4). La surface : peau de
Hermite des sections exactes sous `approx_tol` et `angular_tol` ; **cylindre**
exact quand g est droite et les trois faces planes (nervure droite, avec ou
sans dépouille) ; la ligne des centres et r constants le vérifient.

Le rayon r(s) est une sortie : le congé est l'enveloppe des boules de rayon
r(s) (cercle caractéristique), comme les congés à loi du palier 5. Refus là où
‖r′‖ ≥ ‖c′‖ (les boules s'emboîtent : la face M s'amincit trop vite).

### 3.3 Topologie

- M est ôtée ; F₁ et F₂ sont bornées par les contacts P₁ et P₂ (la co-arête
  avec M remplacée par la courbe de contact, comme le côté d'un congé
  d'arête).
- Aux **bouts** de M, chaque face de bout E (voisine de M, de F₁ et de F₂)
  coupe le congé : la courbe E ∩ congé remplace, dans le contour de E, la
  co-arête avec M, et ses deux bouts sont sur les arêtes E ∩ F₁ et E ∩ F₂
  raccourcies. C'est la fin libre du palier 4, sur une face de bout qui touche
  les deux côtés à la fois.
- M qui fait le tour (bord d'un disque, d'une bride, d'une ailette ; le bord
  d'attaque d'une aube n'en est pas un, il finit sur le pied et la tête) : une
  boucle, pas de bout (§ 3.5).
- La partie de la boule hors du solide (P_M ne doit pas dépasser M) : contrôlée
  comme au § 2.4, les contacts dans leurs faces.

### 3.4 Refus

- aucune boule tangente aux trois faces (Gauss-Newton sans solution, ou r ≤ 0) ;
- M non voisine de F₁ ou de F₂ ;
- un contact hors de sa face (la boule déborde) ;
- emboîtement des boules (§ 3.2) ;
- une face de bout qui ne coupe pas le congé.

---

### 3.5 Face du milieu qui fait le tour (PR 4)

*Ajout du 25 nov. 2026, avant la PR 4 ; revu le même jour après la revue de
conception.*

**Cas.** M fait le tour : ses deux arêtes avec F₁ et F₂ sont fermées, et M n'a
pas de face de bout. Le bord d'un disque ou d'une bride (M un cylindre entre
deux plans), la gorge d'un arbre (M le fond d'une rainure circulaire, concave),
le bord d'une ailette mince (M la surface fermée d'un profil d'aube extrudé,
entre ses deux faces planes). Le bord d'attaque d'une aube n'en est pas un : il
finit sur le pied et la tête, ce sont deux faces de bout (§ 3.3). Aujourd'hui
refusé (« a middle face that closes on itself »).

**Station.** Celle du § 3.1, sur le guide fermé, sans prolongement (pas de face
de bout à couper). Trois changements :

- **la marche** est celle de `fillet_law` (palier 5) : 512 pas sur le guide,
  divisés là où Gauss-Newton échoue ou là où les contacts sautent, avec une
  garde contre le saut qui ne grandit pas avec le pas (des contacts qui bougent
  plus que le guide ne peut les porter, mesuré sur un pas divisé) ; les stations
  que demande la peau sont résolues depuis une voisine de la marche seulement.
  La marche de la PR 3 (2 × 32 pas depuis le milieu, garde proportionnelle au
  pas) enjambe un nez plus court qu'un pas : le plan de la station y recoupe M
  de l'autre côté du nez, où une autre boule est tangente aux trois faces, et
  Gauss-Newton y tombe. C'est un **défaut de la PR 3, déjà dans master** (sonde
  de la revue : le dessus d'une paroi mince en épingle, refusé « approx_tol »
  quand le tour du nez, 0,019, est plus court que le pas, 0,031 ; accepté à
  2 × 256 pas). La PR 4 le corrige pour les guides ouverts comme fermés ;
- M et les côtés qui se ferment par une couture sans être périodiques (une
  B-spline fermée : le profil d'une ailette) sont vus au-delà de leur couture
  (`SeamlessSurface`), comme les supports des congés d'arête fermée ;
- la boucle part du sommet du guide (question 14), sauf si un contact franchit
  la couture d'un côté B-spline fermé non périodique : elle part alors de là,
  comme `fillet_law`. Elle se ferme : la station en b est celle en a (même
  centre, mêmes contacts, à un tour de paramètre près), contrôlé.

**Surface.** Le tore exact du § 3.7 sur un disque, une bride ou une gorge ; la
peau approchée ailleurs, fermée sur elle-même.

**Topologie.** Celle du § 3.3 avec une boucle : sur la copie, M est réduite à
son guide (`collapse` : ses coutures ôtées, l'autre arête fermée remplacée par
le guide dans F₂), puis le congé d'une arête fermée des paliers 4 et 5 (fin en
boucle, `EndKind::Loop`) : F₁ et F₂ coupées le long de contacts fermés, la face
du congé fermée par sa couture.

*Pcourbe du guide dans F₂.* `collapse` la pose aujourd'hui comme le segment
entre les points de l'ancienne co-arête à ses deux bouts. Pour une arête fermée,
ces deux points sont le même : une pcourbe de longueur nulle (vérifié à la
lecture de master). C'est faux aussi pour un guide ouvert dès que l'ancienne
pcourbe n'est pas droite dans (u, v) (une arête en arc de cercle sur un plan) :
le segment est une corde. La PR 4 reprend **l'ancienne pcourbe de M ∩ F₂,
reparamétrée sur le guide** (les stations donnent la correspondance).

**Refus.** Ceux du § 3.4, et :

- une boucle qui ne se ferme pas (les stations en a et en b différentes au-delà
  de `tol`) ;
- les deux refus de la question 18 : la boule plus courbe en travers qu'une face
  à son contact ; la ligne des centres qui se replie le long du guide ;
- la même face des deux côtés avec une boucle (une face annulaire entre deux
  arêtes fermées : la ligne de coupe de la PR 3 ne croise pas son contour) :
  refusée, avec le message qui dit de la couper.

**Sondes** (25 nov. 2026 : les nôtres sur une copie du code de la PR 3 sans le
refus ; celles de la revue de conception) :

- *bord d'un disque* (R = 1, épaisseur 0,3) : la station converge sur tout le
  tour, la couture du cylindre comprise ; r = 0,15000000 en 9 points, centres
  sur le cercle de rayon 0,85 à z = 0,15 ; la station en 2π est celle en 0 ;
  11 ms. Peau à 1,1·10⁻⁶ (le tore n'est pas reconnu aujourd'hui : § 3.7) ;
- *ailette* (le profil fermé de `blade_fixture`, une B-spline fermée non
  périodique, extrudé) : sans vue au-delà de la couture, aucune boule à t = 0,
  la couture du bord de fuite ; avec elle, les stations convergent, puis la peau
  échoue (« approx_tol »), même à r = 0,001. **Cause** (revue) : la marche saute
  le nez du bord de fuite (longueur ≈ π · 0,0035 = 0,011, pas ≈ 0,033), comme
  sur l'épingle ; les nœuds C⁰ du paramètre du profil sont hors de cause (ses
  vitesses y sont égales des deux côtés, et un cercle fait de deux arcs à
  vitesses 1:3 se raccorde sans peine). Non rejoué sur l'ailette ;
- *topologie de la boucle* : non sondée, supposée marcher comme pour une arête
  fermée, une fois la pcourbe du guide dans F₂ reprise.

**Ce que font CATIA et SolidWorks** (attendu, non vérifié faute d'accès) : le
« Full round fillet » de SolidWorks prend la face latérale d'un disque entre ses
deux faces planes et donne le demi-tore ; le « Tritangent Fillet » de CATIA de
même.

### 3.6 Bout arrondi (PR 4)

*Ajout du 25 nov. 2026, avant la PR 4 ; revu après la revue de conception.*

**Cas.** Une face de bout E de M tangente aux deux côtés : le bout arrondi
d'une nervure (une nervure en stade, E un demi-cylindre tangent aux flancs ;
avec la même dépouille des deux côtés, un demi-cône). C'est le bout de nervure
le plus courant en pièce moulée. Aujourd'hui refusé (« an end face tangent to a
side »).

**Géométrie (sondée).** La station au sommet du guide où commence E : sur la
nervure en stade (w = 0,5, h = 1, bout d'axe x = 2), le centre de la boule est
**sur l'axe de E** (à 4·10⁻³⁴), r = 0,25 = w/2, et ses contacts avec les flancs
sont **sur les arêtes entre E et les flancs**, à la hauteur du centre. Cette
boule touche E le long d'un cercle : le cercle de E à la hauteur du centre, du
même rayon. Le congé s'y arrête : au-delà, c'est la boule elle-même.

*Correction d'une note antérieure* (PR 3, § 9 point 11) : les contacts ne
rejoignent pas ces arêtes tangentiellement. Ils les coupent (le contact est
horizontal, l'arête verticale) ; ce qui est tangent, c'est E aux côtés, et le
contact sur E (le demi-cercle) au contact sur chaque flanc.

**Condition** (question 16), contrôlée sous `tol` :

1. E tangente aux deux côtés le long de leurs arêtes communes (`smooth` du
   § 2.1, sous `angular_tol`, en 9 points) ;
2. E de révolution (cylindre, cône) autour d'un axe qui passe par le centre de
   la dernière boule ;
3. la distance de ce centre à E égale à r (un bout plus large ou plus étroit que
   la boule est refusé) ;
4. l'arc de contact de la boule avec E, de P₁ à P₂, dans E.

Avec la même dépouille α des deux côtés, E est un cône d'axe vertical qui passe
par le centre ; la boule le touche le long d'un petit cercle, plus bas que le
centre de r sin α (calcul de la revue, à sonder dans la PR). Avec des
dépouilles différentes, l'axe du cône n'est en général pas dans le plan de la
station : refus (condition 2).

**Surface.** Un morceau de la **sphère** de la dernière boule, exact, borné par
l'arc de la section du congé à la dernière station (de P₁ à P₂ par P_M, un grand
cercle) et par l'arc de contact avec E (un grand cercle sur un cylindre, un petit
cercle sur un cône) : un fuseau de sphère entre deux cercles de la boule, comme
les fuseaux des lois de rayon (palier 5, `BlendLune`, `lune_face`). Sur la
nervure droite en stade, un quart de sphère.

**Topologie.** Une fin nouvelle, sur le fuseau. Ce n'est pas une fin libre : sur
la copie réduite, le sommet du guide a trois arêtes et `free_end` couperait le
congé par E, tangente aux côtés (l'échec d'aujourd'hui). À ce bout :

- le congé finit sur l'arc de sa section à la dernière station ;
- la face sphérique y est ajoutée, entre cet arc et l'arc de contact avec E ;
- E est coupée le long de l'arc de contact : sa partie au-dessus disparaît avec
  M, sa partie en dessous reste (E ne disparaît pas : son arête avec M est
  au-dessus de l'arc de contact ; un bout plus bas que r est déjà refusé par le
  contrôle des contacts) ;
- les arêtes entre E et les côtés sont coupées en P₁ et P₂ (le haut disparaît),
  les contacts sur F₁ et F₂ s'arrêtent là.

**Refus.** E tangente à un côté seulement (question 17) ; une des conditions 2 à 4
manquée (non de révolution, non coaxiale, de rayon différent de r, arc hors de
E), avec la condition dans le message.

**Ce que font CATIA et SolidWorks** (attendu, non vérifié) : sur une nervure
en stade, le congé tritangent du dessus finit en quart de sphère au bout
arrondi.

### 3.7 Surfaces exactes tirées des stations (PR 4)

*Ajout du 25 nov. 2026, avant la PR 4 ; revu après la revue de conception.*

La PR 3 reconnaît le cylindre (guide droit, rayon constant) et le cône (trois
plans qui se coupent en un point). Le tore d'une nervure courbe reste une peau
(déviation 1,4·10⁻⁶). La PR 4 tire la surface **des stations** (question 13),
quelles que soient les faces, par trois propriétés des enveloppes de boules :

- centres sur une droite, rayon constant : **cylindre** ;
- centres sur une droite, rayon linéaire le long d'elle : **cône** (son
  sommet où le rayon s'annule) ;
- centres sur un cercle, rayon constant : **tore** (rayon majeur celui du
  cercle des centres, mineur r).

**Contrôle.** Une propriété lue aux stations seules ne suffit pas : une face
B-spline qui s'écarte de son plan ou de sa révolution **entre** deux stations
ferait reconnaître une primitive fausse (c'est la moitié de contrôle du défaut
de `fillet_faces` du mauvais côté). La primitive est retenue seulement si :

- aux stations **et aux milieux des intervalles**, les arcs exacts recalculés
  sont sur elle à `approx_tol` / 2 et ses normales à `angular_tol` / 2 de celles
  de l'enveloppe (le contrôle de `skinned`) ;
- les contacts tracés sur elle sont sur leurs faces à `tol` (déviation mesurée
  comme dans `on_blend`).

Sinon, la peau.

**Contacts.** Le congé garde le paramètre t du guide. Sondé sur la nervure
courbe : l'angle du guide (une courbe sur surface, issue d'une révolution)
s'écarte d'une loi affine en t de 1,6·10⁻² rad, et les (u, v) du contact sur le
flanc de 2,4·10⁻² ; les segments de la PR 3 n'y suffisent pas. Les pcourbes des
contacts sur leurs faces sont donc obtenues **par projection** des contacts aux
stations, interpolées puis contrôlées sous `tol` (pas « v = constante », qui ne
vaut que pour un flanc de révolution coaxial) ; leurs pcourbes sur la
primitive sont ajustées comme aujourd'hui (`primitive_pcurve`). La ligne des
centres est la droite ou le cercle, paramétrée par t de la même façon.

Cas attendus exacts : la nervure courbe de la PR 3 (tore, volume de Pappus à
10⁻⁹ au lieu de 3,9·10⁻⁸), le bord d'un disque et la gorge d'un arbre (tores),
le U d'une seule face (cylindre, déjà exact), la nervure qui se rétrécit (cône,
déjà exact). Cas qui doit rester une peau, ou passer les contrôles : un flanc
B-spline plan à 10⁻⁷ près entre les stations.

## 4. Ce qui ne change pas

- Les congés d'arête : aucune modification de leur géométrie ni de leur
  topologie. Les points de loi par position sont un ajout.
- Les poches : le code des ponts et de la copie remplie est généralisé (une
  bande bornée par des faces de bout, ou par deux arêtes traitées), les tests
  des poches le vérifient.
- Le message de débordement est commun aux poches et aux congés de faces.

---

## 5. Python

`brep.fillet_faces(m, s, faces1, faces2, r, options)` (r un nombre ou une
`BlendLaw`), `brep.fillet_full_round(m, s, faces1, middle, faces2, options)`,
`brep.BlendPoint.at_point(point, value)` pour les points de loi par position
(`at` est déjà le champ de la fraction : nom changé par la PR 2). Les
mêmes avec `fillet_faces_geometry` pour la surface seule, comme
`fillet_geometry`.

---

## 6. Validation

| Contrôle | Attendu |
|---|---|
| validité | `check` sans défaut, coques fermées, tous les cas |
| arête de boîte chanfreinée, congé du dessus au côté | le solide du congé de l'arête **non chanfreinée** : mêmes faces, cylindre exact, volume de la boîte sans chanfrein moins (1 − π/4) r² L, à 10⁻⁹ ; OCCT : congé de l'arête de la boîte sans chanfrein, même volume à 10⁻⁹ |
| même cas, rayon trop petit | refus du débordement (la boule touche le chanfrein), modèle inchangé |
| pied de bossage cylindrique chanfreiné | bande qui fait le tour : tore exact, volume en forme fermée (Pappus) ; OCCT : congé du pied non chanfreiné |
| bossage en stade (ensemble de 4 faces) chanfreiné au pied | chaîne fermée sur l'arête virtuelle ; OCCT : congé du pied non chanfreiné, volume à 10⁻⁶ |
| gorge en V à fond plat | congé concave, fond englouti, matière ajoutée en forme fermée, volume à 10⁻⁹ |
| bande partielle : U, face en pointe entre paroi et fond ; chanfrein arrêté à mi-longueur ; chanfrein au milieu de l'arête | le solide du congé de l'arête vive du U sans la bande : mêmes faces, cylindre exact d'un seul tenant, volume à 10⁻⁹ ; rayon trop petit : refus du débordement ; OCCT : congé de l'arête vive |
| bande partielle sur faces courbes (pied de bossage chanfreiné sur un demi-tour) | tore d'un seul tenant, volume du congé du pied non chanfreiné à 10⁻⁸ |
| loi par position | les points projetés sur l'arête virtuelle aux abscisses attendues ; volume entre ceux des rayons extrêmes ; même loi donnée par (arête, fraction) sur la boîte sans chanfrein : même solide |
| sections | conique et G2 sur l'arête chanfreinée : les volumes du palier 7 (3 − 3 A pour G2) |
| tritangent, nervure droite de largeur w | demi-cylindre de rayon w/2 exact ; aire ôtée w²/2 − π w²/8 par unité de longueur, volume à 10⁻⁹ |
| tritangent, nervure à dépouille α | cylindre exact, r = (w/2) cos α / (1 − sin α), volume en forme fermée à 10⁻⁹ |
| tritangent, bord de fuite d'aube (intrados et extrados B-splines) | contacts tangents (normales à `angular_tol`), r(s) lisse, solide valide ; surface sous `approx_tol` |
| tritangent, bord d'un disque (rayon R, épaisseur t) (PR 4) | demi-tore exact de rayons R − t/2 et t/2 ; volume par Pappus : deux coins d'aire r² (1 − π/4) tournés au rayon de leur centroïde, R − r + r / (6 (1 − π/4)), à 10⁻⁹ ; OCCT : un cylindre de rayon R − t/2 fusionné avec le tore |
| tritangent, ailette (profil fermé lisse extrudé) (PR 4) | boucle fermée au bord de fuite (couture du profil), solide valide, surface sous `approx_tol` |
| tritangent, nervure en stade (bout arrondi) (PR 4) | cylindre et quart de sphère exacts ; volume ôté 2 r² (1 − π/4) L + π r³/6 par bout arrondi, à 10⁻⁹ ; OCCT : la nervure abaissée fusionnée avec le cylindre et la boule du bout |
| tritangent, nervure courbe (PR 4) | tore exact, volume de Pappus à 10⁻⁹ |
| tritangent, gorge d'un arbre (boucle concave) (PR 4) | demi-tore exact, matière ajoutée par Pappus à 10⁻⁹ ; OCCT : l'arbre à gorge abaissée fusionné avec le tore |
| tritangent, dessus d'une paroi mince en épingle (PR 4) | accepté quel que soit le tour du nez (master le refuse quand le tour est plus court qu'un pas de la marche), solide valide |
| tritangent, ailette au profil de `blade_fixture` (PR 4) | boucle fermée au bord de fuite (couture du profil) ; r sous le rayon de courbure du bord de fuite (0,0035) : solide valide, surface sous `approx_tol` |
| refus de la question 18 (PR 4) | (b) l'ailette au-delà de r = 0,0035 : la ligne des centres replie ; (a) une nervure à dessus bombé plus courbe que la boule : l'arc sortirait de M ; chacun avec la valeur et sa direction dans le message |
| surface reconnue à tort (PR 4) | un flanc B-spline plan à 10⁻⁷ près entre les stations : une peau, ou une primitive qui passe les contrôles du § 3.7 |
| bout arrondi en cône, même dépouille des deux côtés (PR 4) | fuseau de sphère entre un grand cercle et un petit cercle, exact ; volume en forme fermée |
| refus du bout arrondi (PR 4) | dépouilles différentes (non coaxial) ; rayon de E différent de r ; tangent d'un seul côté ; modèle inchangé |
| même face des deux côtés avec une boucle (PR 4) | refus, avec le message qui dit de couper la face |
| refus | faces que rien ne relie (deux nervures sur une plaque) ; flancs parallèles en congé de faces (pas d'arête virtuelle : le message propose le congé tritangent) ; plusieurs bandes |
| galerie | arête chanfreinée reprise en congé, gorge en V comblée, nervure arrondie, bord de fuite d'aube arrondi |

OCCT n'a pas de congé de faces ni de congé tritangent : la comparaison passe par
les cas équivalents (le congé de l'arête qui existait avant le chanfrein). Pour
le tritangent, les volumes en forme fermée seulement.

---

## 7. Plan des PR

| # | PR | Contenu | Estimation |
|---|---|---|---|
| 1 | congé de faces | § 2.1 à § 2.4 : bande, faces de bout, arête virtuelle (généralisation de `Bridge` et `fill`), congé sur la copie, contrôles, ensembles de faces, rayon constant ; Python, galerie, OCCT par cas équivalents ; ce document | 18 h |
| 2 | bande partielle, lois et sections des congés de faces | § 2.6 : bande contiguë à une arête commune, arête virtuelle fondue avec elle (10 h) ; § 2.5 : points de loi par position (aussi pour les congés d'arête), lois et sections du palier 7 sur l'arête virtuelle, Python, galerie (8 h ; la gorge concave est déjà faite par la PR 1) | 18 h |
| 3 | congé tritangent et bilan | § 3 : station à trois faces, section par P_M, cylindre exact, topologie (M ôtée, faces de bout), nervures et bord de fuite d'aube, Python, galerie, bilan du palier | 18 h |
| 4 | tritangent : boucle, bout arrondi, tore | § 3.5 : M qui fait le tour, marche de `fillet_law` (corrige aussi l'épingle de la PR 3), pcourbe du guide dans F₂, deux refus de repli (9 à 12 h) ; § 3.6 : bout arrondi, fuseau de sphère, fin nouvelle, condition de coaxialité (11 à 13 h) ; § 3.7 : surfaces exactes tirées des stations, contrôle entre les stations, pcourbes par projection (4 à 5 h) ; Python, galerie, OCCT par cas équivalents | 24 à 30 h |
| 5 | tritangent sur des ensembles de faces (annoncée) | une chaîne de guides : côtés et milieu en faces jointes en tangence (question 6 du congé de faces) ; la réduction de M quand les bords des faces des trois ensembles ne tombent pas aux mêmes endroits | 14 à 18 h |
| | | **Total** | **≈ 94 à 103 h** |

---

## 8. Questions à valider

1. **Périmètre** : congés de faces et congés tritangents dans un même palier ;
   débordement, faces que rien ne relie et plusieurs bandes laissés dehors.
   *Recommandé* : les deux opérations de CATIA et SolidWorks qui manquent, le
   tritangent étant attendu pour les bords d'aube.
2. **Arête virtuelle sur une copie** (généralisation des poches), puis la
   machinerie des congés d'arête. *Recommandé* : toute la topologie, les lois et
   les sections viennent avec, sans booléen à contacts tangents.
3. **Bande bornée par des faces de bout**, une seule par appel, la plus petite
   qui passe le contrôle. *Recommandé.*
4. **Côté de la boule donné par la matière** (convexe : ôte, concave : ajoute),
   sans choix de l'utilisateur. *Recommandé* : c'est le comportement des congés
   d'arête, et le seul qui garde un solide.
5. **Contrôles des poches** (contacts dans F₁ et F₂, congé qui ne touche pas la
   bande) et refus du débordement. *Recommandé.*
6. **Ensembles de faces jointes en tangence** de chaque côté, l'arête virtuelle
   en chaîne. *Recommandé* : indispensable pour les aubes et les bossages.
7. **Points de loi par position**, ouverts aussi aux congés d'arête.
   *Recommandé.*
8. **Station tritangente** dans le plan normal au guide M ∩ F₁, rayon inconnu,
   section par le contact sur M ; pas d'arête virtuelle. *Recommandé* : elle
   marche entre flancs parallèles, et sert plus tard au débordement.
9. **Topologie tritangente** : M ôtée, fins libres sur les faces de bout.
   *Recommandé.*
10. **Validation** par cas équivalents (arête qui existait avant le chanfrein)
    pour OCCT, formes fermées pour les nervures, bord de fuite d'aube sur
    B-splines. *Recommandé.*
11. **Trois PR** (congé de faces ; lois et sections ; tritangent et bilan).
    *Recommandé.*
12. **Bande partielle** (§ 2.6, ajoutée le 24 nov. 2026) : une bande contiguë à
    l'arête commune de F₁ et F₂ (face en pointe, chanfrein arrêté, bande au
    milieu) est engloutie par `fillet_faces`, l'arête virtuelle fondue avec
    l'arête commune ; dans la PR 2, avant les lois. `fillet` d'une arête n'y va
    pas de lui-même, son refus nomme `fillet_faces`. *Recommandé* : c'est le
    manque d'OpenCascade (un congé qui supprime une face), le même code que la
    PR 1, et l'utilisateur choisit les deux faces comme dans CATIA ; faire
    passer `fillet` d'une arête par-dessus une bande à son bout changerait le
    sens d'un appel qui réussit ou échoue aujourd'hui.

13. **Surfaces exactes tirées des stations** (§ 3.7) : cylindre, cône ou tore
    reconnus sur les centres et les rayons des stations, quelles que soient les
    faces ; retenus seulement si, aux stations et entre elles, les arcs exacts
    sont sur la primitive (`approx_tol` / 2, normales à `angular_tol` / 2) et
    les contacts tracés sur elle sur leurs faces (`tol`) ; pcourbes des contacts
    par projection. *Recommandé* : c'est une propriété des enveloppes de boules,
    elle couvre les faces B-spline planes ou de révolution, et le contrôle entre
    les stations écarte une primitive fausse. *Écarté* : reconnaître les surfaces
    par le type des faces (plans, cylindres coaxiaux), qui manque les B-splines
    et oblige à un cas par paire de types ; ou un contrôle aux stations seules.
14. **M qui fait le tour : la boucle part du sommet du guide** (§ 3.5), comme
    un congé d'arête fermée (palier 4), sauf si un contact franchit la couture
    d'un côté B-spline fermé non périodique : elle part alors de là, comme
    `fillet_law`. *Recommandé* : la machinerie des arêtes fermées s'applique
    telle quelle. *Écarté* : partir de la couture de M, qui demanderait de
    recoller le guide à cet endroit.
15. **Bout arrondi : un fuseau de la dernière boule** (§ 3.6), entre l'arc de
    la section et l'arc de contact avec E, E coupée le long de ce dernier, par
    une fin nouvelle. *Recommandé* : la sphère est exacte, et le fuseau reprend
    celui des lois de rayon (palier 5, `lune_face`), à adapter : il y borde deux
    congés, ici un congé et une face. *Écarté* : faire rouler la boule autour du
    bout, E prise comme un troisième côté (un second congé tritangent entre F₁,
    E et F₂) : sur un bout de révolution coaxial, c'est la même sphère, avec une
    station dégénérée (la boule touche E le long d'un cercle).
16. **Bout arrondi accepté sous quatre conditions** (§ 3.6) : E tangente aux deux
    côtés ; de révolution autour d'un axe par le dernier centre ; à la distance r
    de ce centre ; l'arc de contact dans E. Un bout en cône à même dépouille est
    accepté, à dépouilles différentes refusé ; une face de bout B-spline
    tangente refusée. *Recommandé* : le cas des pièces moulées, exact.
    *Écarté* : un raccord de sommet du palier 6 (`fill_n_sided`) pour toute
    face de bout tangente : plus général, approché, à reprendre avec les
    raccords G2 si un cas le demande.
17. **Bout arrondi d'un côté seulement** (E tangente à un côté, vive de
    l'autre) refusé. *Recommandé* : la boule ne touche alors E qu'en un point,
    il n'y a ni fuseau ni fin libre propre ; CATIA et SolidWorks ne le
    proposent pas non plus à notre connaissance (non vérifié). *Écarté* : une
    fin mixte (fuseau d'un côté, fin libre de l'autre).
18. **Deux refus de repli**, que `check` ne voit pas, sur les trois faces à
    leurs contacts, avec la valeur trouvée et sa direction dans le message :
    (a) **en travers** : la courbure normale d'une face à son contact, dans la
    direction de l'arc de la section, plus grande que 1/r (la boule plus courbe
    que la face : l'arc en sort, le congé ajoute de la matière) ; (b) **le long**
    : la ligne des centres qui se replie, le rayon au-delà d'un rayon de courbure
    d'une face décalée (le contrôle de repli de `fillet_law`, `offset_surface`
    des supports au rayon le plus grand). L'ailette relève de (b) : son profil
    est droit en travers, et son bord de fuite (0,0035) est une courbure du
    guide. *Recommandé* : deux défauts distincts, chacun avec son test.
    *Écarté* : laisser le congé se construire et compter sur `check`, qui ne
    voit ni un arc sorti de la matière ni une surface qui se recoupe ; ou un
    seul contrôle, sur la courbure de M seule.

### Sources

- Vida, Martin et Várady, « A survey of blending methods that use parametric
  surfaces », *Computer-Aided Design* 26(5), 1994 : boule roulante entre deux
  surfaces, congés à trois surfaces.
- B. K. Choi, S. Y. Ju, « Constant-radius blending in surface modelling »,
  *Computer-Aided Design* 21(4), 1989 : congés de faces par surfaces décalées.
- J. Hoschek, D. Lasser, *Fundamentals of Computer Aided Geometric Design*,
  § 14 (*blending*).
- Le comportement de CATIA (« Face-Face Fillet », « Tritangent Fillet ») et de
  SolidWorks (« Face fillet », « Full round fillet », *Help points*). OCCT
  (`BRepFilletAPI_MakeFillet`) en boîte noire, par cas équivalents.

---

## 9. Bilan du palier

Les trois PR ont été menées. Les deux opérations de CATIA et de SolidWorks qui
manquaient sont là :

| Opération | Ce qu'elle fait | Surface exacte | PR |
|---|---|---|---|
| congé de faces (`fillet_faces`) | deux faces, ou deux ensembles de faces jointes en tangence, sans arête commune ; la bande entre elles engloutie (chanfrein, marche, fond de gorge) ; convexe ou concave | celles du congé d'arête sur la copie : cylindre, tore | 1 |
| bande partielle | des faces qui se touchent déjà sur une partie de leur longueur : la face en pointe, le chanfrein arrêté, le chanfrein au milieu, engloutis ; l'arête virtuelle fondue avec l'arête commune | un seul tenant (cylindre, tore) | 2 |
| lois et sections | `RadiusLaw` le long de l'arête virtuelle, points par position (`BlendPoint::at_point`, aussi pour les congés d'arête) ; sections conique, corde constante, G2 | celles du palier 7 | 2 |
| congé tritangent (`fillet_full_round`) | la face du milieu remplacée par le congé tangent aux trois faces, le rayon donné par la géométrie ; nervures, fond de rainure, bord de fuite épais d'une aube, d'un carreau ou de trois | cylindre (nervure droite, avec ou sans dépouille) ; cône (nervure qui se rétrécit entre trois plans) | 3 |

Chaque cas est comparé au congé de l'arête qui existait avant la bande (congé
de faces) ou à sa forme fermée (tritangent), et à OpenCascade par des cas
équivalents : le congé de l'arête sans chanfrein (`BRepFilletAPI_MakeFillet`),
ou la fusion d'une nervure abaissée et d'un cylindre.

**Écarts au document.** Validés au fil des notes, ou à valider dans la
dernière :

- PR 1 : plusieurs bandes, la plus petite qui passe le contrôle, sans refus ;
  le congé des arêtes communes seulement si les faces se touchent partout
  (revue de code).
- PR 2 : `BlendPoint::at_point` au lieu de `at` ; `fillet_faces_chordal` ;
  points d'une loi changés en positions avant la copie ; une géométrie commune
  de la corde constante (revue de code).
- PR 3 ([brep8_pr03_full_round.md](brep8_pr03_full_round.md), § 9) : topologie
  par réduction de la face du milieu à son guide sur une copie ; section par le
  cercle caractéristique des trois contacts ; côté concave accepté ;
  orientation de la face du congé par ses contacts (`blend.h`) ; une face par
  côté et au milieu ; la même face des deux côtés coupée sur la copie puis
  refondue ; bout arrondi refusé ; tore approché (revue de code, § 10).

**Ce qui reste hors périmètre**, refusé avec un message clair :

- le **débordement** : la boule qui passe d'une paire de faces à une autre
  par-dessus une arête vive. La station à trois faces de la PR 3 en est la
  brique (le point où la boule touche trois faces) ;
- deux faces **que rien ne relie** (deux nervures sur une plaque) ;
- un chanfrein qui finit sur un autre chanfrein (deux faces de bout, chacune ne
  touchant qu'un côté) ;
- pour le tritangent, en dessous de la cible CATIA / SolidWorks, en trois rangs
  (note de la PR 3, § 9, point 4, validé comme limite provisoire) :
  - **à faire dans la suite** :
    - les **ensembles de faces** par côté ou au milieu (un flanc en plusieurs
      carreaux joints en tangence : une chaîne de guides) ;
    - la face du milieu **qui fait le tour** (le bord d'un disque, d'une bride,
      d'une ailette, la gorge d'un arbre : une boucle sans face de bout) ;
    - le **bout arrondi** : une face de bout tangente aux côtés (une nervure en
      stade, le bout de nervure le plus courant en pièce moulée), où CATIA et
      SolidWorks finissent le congé par un quart de sphère (10 à 12 h, note de
      la PR 3, § 9, point 11) ;
  - **refusé par nature** : la face du milieu qui finit **en pointe** ou en
    **proue** (plusieurs faces à un bout) : le rayon y tend vers zéro ;
  - **refusé pour l'instant, contournable en coupant la face** : la face du
    milieu **trouée**, ou qui touche un côté par **plusieurs arêtes** ;
  - la même face des deux côtés quand ses arêtes avec la face du milieu ne sont
    pas isoparamétriques dans une même direction, quand la ligne qui la coupe
    à mi-chemin croise son contour un autre nombre de fois que deux, ou quand elle se ferme
    par une couture (refusée : la couper avant ; ses bouts, eux, peuvent être
    quelconques, une aube d'un carreau fusionnée à son moyeu passe) ;
  - le **tore exact** d'une nervure courbe (une peau à 10⁻⁶ aujourd'hui) ;
  - les sections conique et G2 ;
- les raccords de sommet G2 (palier 7).

**Suite possible de la feuille de route :**

- la **PR 4** du palier (§ 3.5 à § 3.7, lancée le 25 nov. 2026) : la face du
  milieu **qui fait le tour** (une boucle sans face de bout : le bord d'un
  disque, d'une bride, d'une ailette), le **bout arrondi** (un fuseau de
  sphère), et le **tore exact** avec les surfaces tirées des stations ;
- la **PR 5**, annoncée : le tritangent sur des **ensembles de faces** (une
  chaîne de guides, comme les chaînes du palier 4) ;
- le **débordement** des congés (*rolling ball overflow*), à partir de la
  station à trois faces ;
- la marche des stations du tritangent (PR 3, master) enjambe un nez plus court
  qu'un pas (le dessus d'une paroi mince en épingle, refusé) : corrigée par la
  PR 4 (§ 3.5) ;
- hors du palier, constaté par la revue de conception de la PR 4 : le congé
  d'arête de l'arête du dessus de l'ailette au profil de `blade_fixture` est
  refusé à r = 0,0005 et 0,001 (« the centre line of a closed edge is not
  closed »), accepté à 0,002. Non analysé ;
- la comparaison avec CATIA et SolidWorks des congés tritangents sur des aubes
  réelles (rayon le long du bord, raccord aux faces de bout).

