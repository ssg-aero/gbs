# Noyau BREP natif — conception du palier 8 : congés de faces et congés tritangents

| | |
|---|---|
| Statut | Questions du § 8 validées (23 nov. 2026) ; PR 1 : [brep8_pr01_face_blends.md](brep8_pr01_face_blends.md) ([illustration des choix](https://claude.ai/artifact/HjZDQoSQMBz1MQrWqqnaY4)) |
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
  - un pied de bossage chanfreiné repris en congé.
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
- M qui fait le tour (bord d'un disque, bord d'attaque d'une aube fermée sur
  elle-même) : une boucle, pas de bout.
- La partie de la boule hors du solide (P_M ne doit pas dépasser M) : contrôlée
  comme au § 2.4, les contacts dans leurs faces.

### 3.4 Refus

- aucune boule tangente aux trois faces (Gauss-Newton sans solution, ou r ≤ 0) ;
- M non voisine de F₁ ou de F₂ ;
- un contact hors de sa face (la boule déborde) ;
- emboîtement des boules (§ 3.2) ;
- une face de bout qui ne coupe pas le congé.

---

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
`brep.BlendPoint.at(point, value)` pour les points de loi par position. Les
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
| loi par position | les points projetés sur l'arête virtuelle aux abscisses attendues ; volume entre ceux des rayons extrêmes ; même loi donnée par (arête, fraction) sur la boîte sans chanfrein : même solide |
| sections | conique et G2 sur l'arête chanfreinée : les volumes du palier 7 (3 − 3 A pour G2) |
| tritangent, nervure droite de largeur w | demi-cylindre de rayon w/2 exact ; aire ôtée w²/2 − π w²/8 par unité de longueur, volume à 10⁻⁹ |
| tritangent, nervure à dépouille α | cylindre exact, r = (w/2) cos α / (1 − sin α), volume en forme fermée à 10⁻⁹ |
| tritangent, bord de fuite d'aube (intrados et extrados B-splines) | contacts tangents (normales à `angular_tol`), r(s) lisse, solide valide ; surface sous `approx_tol` |
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
| 2 | lois et sections des congés de faces | § 2.5 : points de loi par position (aussi pour les congés d'arête), lois et sections du palier 7 sur l'arête virtuelle, gorge concave, Python, galerie | 10 h |
| 3 | congé tritangent et bilan | § 3 : station à trois faces, section par P_M, cylindre exact, topologie (M ôtée, faces de bout), nervures et bord de fuite d'aube, Python, galerie, bilan du palier | 18 h |
| | | **Total** | **≈ 46 h** |

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
