# Noyau BREP natif — conception du palier 2 bis

| | |
|---|---|
| Statut | À valider |
| Périmètre | Lever les limites du palier 2 (bilan [brep_palier2.md](brep_palier2.md), § 12.4) : contacts tangents, arêtes posées sur une face, points singuliers, pôles, faces confondues quelconques, restriction exacte des courbes analytiques, nettoyage après opération |
| Dépend de | palier 2 (intersections, section, découpe, booléens) |
| Hors périmètre | congés et chanfreins, balayage le long d'un chemin, coques épaisses, maillage (palier 3) |

Ce document décrit *quoi* construire et *pourquoi*, comme
[brep_palier2.md](brep_palier2.md). Chaque PR a sa note d'architecture
`brep2b_prNN_*.md`.

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Contacts en un point | ils ne découpent rien ; les faces sont classées par **plusieurs points intérieurs**, le premier hors du bord de l'autre opérande fait foi | refus `NonGeneric` (palier 2) |
| Contacts le long d'une ligne | l'arête de tangence est une arête de section comme les autres ; elle peut **ne pas découper** une face qu'elle touche sans la traverser | refus |
| Arêtes de plus de deux faces | les faces sont **appariées autour de l'arête**, chacune avec sa voisine du côté de sa matière ; chaque paire sauf la première reçoit **sa copie de l'arête** : les shells restent des variétés | shell non-variété (arête partagée par quatre faces) |
| Arêtes posées sur une face | l'arête de l'un **imprime** la face de l'autre (découpée à ses traversées du bord, morceaux intérieurs en coupes), comme pour les faces confondues | refus |
| Points de branchement | les directions de départ viennent de la **différence des formes fondamentales secondes** au point ; on chemine sur chaque demi-branche ; le point devient un **sommet de degré 4** | refus `singular` |
| Tangence le long d'une courbe hors formules fermées | **refusée** (`NonGeneric`) | suivi d'une courbe où le jacobien est singulier |
| Pôles | la courbe qui traverse un pôle est **coupée au pôle** (deux arêtes, un sommet) ; l'arête dégénérée du pôle est **découpée** dans l'espace (u, v) | une arête à pcurve discontinue |
| Faces confondues quelconques | coïncidence testée par **projection d'une grille** de points d'une face sur la surface de l'autre ; pcurves des arêtes imprimées par **projection** | primitives seulement (palier 2) |
| Restriction aux boîtes | intersection **courbe-courbe** des pcurves avec les côtés de la boîte (toutes les solutions) | échantillonnage à 512 pas |
| Nettoyage | les arêtes et sommets **créés** par l'opération et inutilisés sont effacés, en cas de succès comme d'échec | `compact()` à la charge de l'appelant |

---

## 1. Objectifs et non-objectifs

### 1.1 Objectifs

Rendre les booléens, la section et la découpe utilisables sur les
configurations que le palier 2 refusait :

1. **contacts** en un point ou le long d'une ligne (solides qui se touchent,
   formes posées l'une sur l'autre) ;
2. **arêtes posées sur une face** de l'autre opérande ;
3. **points singuliers** des courbes d'intersection (points de branchement :
   cylindres égaux qui se croisent, courbe de Viviani) ;
4. **pôles** traversés par une courbe obtenue par cheminement ;
5. **faces confondues** entre surfaces quelconques (NURBS, primitive et NURBS
   de même forme) ;
6. **restriction exacte** des courbes analytiques aux boîtes de paramètres ;
7. **nettoyage** des entités inutilisées après une opération.

### 1.2 Non-objectifs

- Tangence le long d'une courbe hors formules fermées (§ 5.3).
- Coïncidence **partielle** de deux surfaces quelconques (confondues sur une
  partie seulement de leur recouvrement) : refusée.
- Congés, chanfreins, balayage le long d'un chemin, coques épaisses, maillage :
  palier 3.

---

## 2. Inventaire des refus du palier 2

| Refus (`boolean.h`) | Cause | Traité |
|---|---|---|
| « a vertex touches a face of the other operand » | `SectionResult::points` non vide : contact ponctuel | PR 1 |
| « faces tangent along a section edge » | arête de section `tangent` | PR 1 |
| « a section edge does not split a face » | arête pendante : une ligne de contact s'arrête à l'intérieur d'une face | PR 1 (arêtes tangentes) |
| « an edge lies on a face of the other operand » | `boundary_edges` hors faces confondues | PR 2 |
| « singular or incomplete intersection » | cheminement arrêté sur un point singulier | PR 3 (branchements), § 5.3 (tangences) |
| « coincident faces that are not primitives » | `imprint_coincident` limité aux primitives | PR 5 |
| deux morceaux non raccordés au pôle | cheminement à travers un pôle | PR 4 |

---

## 3. Contacts (PR 1)

### 3.1 Contact en un point

Un point de contact (sphère posée sur une boîte, sommet d'un cône sur un
plan) ne découpe aucune face : il suffit de ne plus le refuser. Le risque est
dans la **classification** : le point intérieur d'une face, milieu de sa plus
longue corde, peut tomber exactement sur le point de contact (la face du
dessus de la boîte en son centre). La face serait alors classée `On`.

`interior_points` rend plusieurs points intérieurs (un par ligne de v
constant, plus longues cordes d'abord). `classify(face, solide)` prend le
premier qui n'est pas sur le bord du solide ; la face n'est `On` que si tous
le sont (face confondue).

### 3.2 Contact le long d'une ligne

Les formules fermées donnent la génératrice de tangence d'un plan et d'un
cylindre, marquée `tangent`. C'est une arête de section ordinaire : elle
découpe la face qu'elle traverse de bord à bord (le dessus de la boîte) et
reste pendante dans la face qu'elle touche sans la traverser (le flanc du
cylindre). Une arête pendante est acceptée si elle est tangente, et refusée
sinon : une arête de section qui traverse une face sans la découper trahit
une intersection manquée.

Une correction de la section accompagne ce point. Un morceau de courbe proche
du bord d'une face (`On`) n'est une arête de bord que si l'intersection
courbe-courbe a trouvé un **recouvrement** à cet endroit. Le cercle d'un plan
latéral de la boîte, tangent à l'arête du dessus, n'y laisse qu'un point de
contact. Il était auparavant pris pour une arête posée sur la face.

### 3.3 Arêtes de plus de deux faces

Deux solides qui se touchent le long d'une ligne gardent quatre faces autour
de l'arête commune (fusion de la boîte et du cylindre posé dessus ; deux
copeaux de la différence qui se touchent). Un shell où une arête sert quatre
fois n'est pas une variété (`check` : `NonManifoldEdge`).

**Appariement radial.** Pour chaque usage de l'arête :
- on prend un point de la face à une petite distance h du milieu de l'arête,
  dans la direction qui entre dans la face (n × t, la face étant à gauche de
  son bord vu de l'extérieur), ramenée dans l'espace (u, v) par moindres
  carrés ;
- on mesure son angle autour de l'arête.

À distance finie, deux faces tangentes le long de l'arête se séparent par leur
courbure (écart d'angle ≈ h κ / 2). On trie les usages par angle. Chacun
s'apparie avec son voisin du côté de sa matière, derrière la normale
extérieure. L'appariement doit être symétrique, sinon `NonGeneric`.

Chaque paire sauf la première reçoit une **copie de l'arête** (même courbe,
mêmes sommets). Les shells sont ensuite les composantes connexes par arêtes
partagées, comme au palier 2. La fusion de solides qui se touchent donne deux
solides ; la différence pincée donne deux copeaux distincts. OCCT donne les
mêmes nombres de solides et de faces (§ 9).

### 3.4 Nettoyage

L'opération note la taille des tables d'arêtes et de sommets au départ. À la
fin, en cas de succès comme d'échec, elle efface les arêtes et sommets créés
depuis que rien n'utilise. Cela vise surtout les arêtes des faces non
retenues, les arêtes pendantes et les points de contact.
`erase_shape` listait les faces après avoir effacé les shells qui y mènent :
les faces non retenues des copies restaient dans le modèle. Elles sont
maintenant listées avant tout effacement.

---

## 4. Arêtes posées sur une face (PR 2)

Une boîte posée sur l'arête d'une autre, ou la génératrice d'un cylindre
posée sur une face tangente, laissent une arête d'un opérande **dans** une
face de l'autre, sans faces confondues. La section la repère (`On` le long
d'un recouvrement) et la range dans `boundary_edges`.

Traitement proposé, sur le modèle de l'empreinte des faces confondues
(`coincident.h`) :
1. la section rapporte l'arête **avec la face** qu'elle touche
   (`SectionResult::edges_on_faces`) ;
2. l'arête est découpée là où elle traverse le bord de cette face, ou là où une
   arête de ce bord se termine sur elle (intersection courbe-courbe 3D) ;
   les sommets proches sont fusionnés ;
3. les morceaux intérieurs à la face sont des **coupes** de la face, avec leur
   pcurve calculée sur sa surface (forme fermée pour une primitive, projection
   sinon, § 7.2) ;
4. une arête qui longe une arête de l'autre opérande (deux boîtes qui se
   touchent par une arête) est fusionnée avec elle (étape 2 de l'empreinte).

Après la découpe, l'arête a trois ou quatre faces autour d'elle : l'appariement
radial (§ 3.3) sépare les solides.

Les sommets posés sur une face (coin d'une boîte sur une autre) ne
découpent rien. Un sommet posé sur une arête la découpe (`split_at_vertex`)
quand une arête de section y aboutit, comme au palier 2.

---

## 5. Points singuliers (PR 3)

### 5.1 Directions de départ

En un point singulier p, les normales des deux surfaces sont parallèles. Dans
le plan tangent commun de repère (e1, e2), chaque surface est localement le
graphe z = ½ xᵀ Hₖ x, où Hₖ est sa seconde forme fondamentale exprimée dans ce
repère. Les directions de la courbe d'intersection en p sont les solutions de
xᵀ D x = 0, avec D = H₁ − H₂ :

| det D | Point | Traitement |
|---|---|---|
| < 0 | **branchement** : deux directions, quatre demi-branches | cheminement depuis p + ε d sur chaque demi-branche non encore parcourue |
| > 0 | contact isolé | point de contact (déjà rapporté) |
| ≈ 0, D ≠ 0 | tangence le long d'une courbe | § 5.3 |
| D ≈ 0 | contact d'ordre supérieur | `singular` |

Le cheminement s'arrête aujourd'hui en marquant `singular` dès qu'il approche
un tel point. Il s'arrêtera désormais **sur** le point, et chaque branche qui
y arrive y aura son extrémité.

### 5.2 Topologie

Le point de branchement est une extrémité commune de quatre courbes. La
section le partage par la réserve de sommets (`VertexPool`) : un sommet de
degré 4, laissé en arêtes non chaînées si les wires ne peuvent pas
l'absorber. Dans chaque face, la découpe trouve un nœud de degré 4 où
l'arrangement choisit la demi-arête suivante dans le sens horaire, ce qu'il
sait déjà faire.

Tests :
- cylindres égaux d'axes concourants (Steinmetz) : l'intersection a pour
  volume 16 R³ / 3 ;
- sphère et cylindre de rayon moitié passant par son centre (courbe de
  Viviani, un point double) ;
- comparaison à OCCT.

### 5.3 Tangence le long d'une courbe

Hors des formules fermées (primitives coaxiales, plan et quadrique), la
tangence le long d'une courbe (tore dans un cylindre, deux NURBS tangentes)
rend le jacobien du cheminement singulier en tout point. Elle reste
**refusée** : c'est rare, et les primitives couvrent les cas courants.

---

## 6. Pôles (PR 4)

Au pôle d'une sphère (ou au sommet d'un cône), le paramétrage dégénère :
|S_u × S_v| → 0, toute la ligne v = ±π/2 est le même point. Une courbe
d'intersection qui traverse le pôle reste continue en 3D, mais sa pcurve saute
de u à u + π.

1. **Cheminement.** Quand le pas approche un pôle (|S_u × S_v| petit devant
   |S_u| |S_v|), la courbe est arrêtée au pôle, puis le cheminement reprend de
   l'autre côté. On prédit en 3D le long de la tangente de la courbe, et on
   reprojette avec u choisi du côté opposé. Les deux morceaux sont rapportés
   comme deux courbes qui se terminent au pôle.
2. **Section.** Les deux arêtes partagent le sommet du pôle, qui existe déjà
   dans la face (extrémité de son arête dégénérée).
3. **Découpe.** Dans l'espace (u, v), les coupes aboutissent sur l'arête
   dégénérée du pôle en des u qui ne sont pas ses extrémités. L'arête
   dégénérée est découpée en ces u (pcurve en segments, même sommet), pour
   que l'arrangement ait ses nœuds.

Les formules fermées (plan passant par les pôles d'une sphère) coupent déjà
leurs courbes au pôle : elles servent de référence au test du cheminement.

---

## 7. Faces confondues quelconques (PR 5)

### 7.1 Détection

Le cheminement ne voit pas une coïncidence : il rend des départs partout ou
rien de stable. Avant d'intersecter deux faces dont les boîtes se recouvrent,
une grille de points intérieurs de la plus petite est projetée sur la surface
de l'autre. Si tous les points sont à moins de la tolérance, et que les
normales sont parallèles, les faces sont **confondues** sur leur recouvrement
(`SectionResult::coincident`).

Une coïncidence partielle (une partie des points seulement) est refusée
(§ 1.2).

### 7.2 Empreinte

`imprint_coincident` (`coincident.h`) ne dépend des primitives que par deux
points, qui sont généralisés :
- **pcurve d'une arête sur une face** : forme fermée pour une primitive, sinon
  projection d'échantillons de l'arête et interpolation, avec l'écart mesuré
  (comme les pcurves du palier 1, `pcurve.h`) ;
- **côté d'un seam** : le décalage d'une période se lit sur la fermeture de
  la surface (`surface_closure`, `closure.h` : `closed_u`, `closed_v` et
  largeur des bornes) au lieu du type de primitive.

`same_side` (choix du morceau gardé dans un booléen) projette sur la surface
quelconque (`project_on_surface`) au lieu de `primitive_uv`.

---

## 8. Restriction exacte des courbes analytiques (PR 6)

Les courbes des formules fermées (droites, cercles, ellipses) sont restreintes
aux boîtes de paramètres des deux faces. Le palier 2 échantillonne 512 pas et
affine par dichotomie : un morceau plus court qu'un pas peut être manqué.

La restriction devient une **intersection courbe-courbe** de chaque pcurve
avec les quatre côtés de sa boîte (`gbs/intersect/curves.h`, toutes les
solutions). Les paramètres trouvés découpent la courbe ; chaque morceau est
gardé si son milieu est dans les deux boîtes. Les pcurves exactes (segments)
donnent des solutions exactes, les pcurves interpolées des solutions à leur
écart près.

---

## 9. Validation

Chaque PR ajoute ses cas aux tests de gbs et à la comparaison avec OCCT
(`gbs-occt/tests/tests_brep2_occt.cpp`, mêmes critères qu'au palier 2) :

| PR | Cas comparés à OCCT |
|---|---|
| 1 | sphère posée sur une boîte, cylindre couché sur une boîte, deux copeaux tangents (boîte contre cylindre) |
| 2 | boîte posée sur une arête, boîtes qui se touchent par une arête, génératrice de seam sur un plan tangent |
| 3 | cylindres égaux concourants, Viviani |
| 4 | sphère dont le pôle est traversé par une courbe de cheminement |
| 5 | faces confondues NURBS et NURBS, primitive et NURBS de même forme |
| 6 | courbes qui n'effleurent qu'un coin de boîte |

La dernière PR ajoute les liaisons Python éventuelles et le bilan de ce
document.

---

## 10. Plan des PR

| # | PR | Contenu | Estimation |
|---|---|---|---|
| 1 | contacts (avec ce document) | § 3 : points et lignes de contact, appariement radial, nettoyage | 8 h |
| 2 | arêtes posées sur une face | § 4 | 10 h |
| 3 | points de branchement | § 5 | 12 h |
| 4 | pôles | § 6 | 10 h |
| 5 | faces confondues quelconques | § 7 | 10 h |
| 6 | restriction exacte, bilan | § 8, Python, documentation | 6 h |
| | | **Total** | **≈ 56 h** |

---

## 11. Questions à valider

1. **Solides qui se touchent : séparés.** La fusion de deux solides qui se
   touchent en un point ou le long d'une ligne donne deux solides (un
   composé), comme OCCT. Une arête de contact a une copie par solide.
   *Recommandé.* Alternative : un shell non-variété.
2. **Arêtes pendantes acceptées si tangentes seulement.** Une arête de section
   non tangente qui ne découpe pas une face reste une erreur. *Recommandé.*
3. **Tangence le long d'une courbe hors formules fermées : refusée** (§ 5.3).
   *Recommandé.*
4. **Coïncidence quelconque : totale seulement** sur le recouvrement des deux
   faces (§ 7.1). *Recommandé.*
5. **Courbe coupée au pôle** en deux arêtes plutôt qu'une arête à pcurve
   discontinue (§ 6). *Recommandé.*
6. **Nettoyage automatique** des arêtes et sommets créés par une opération et
   inutilisés (§ 3.4), au lieu d'un `compact()` laissé à l'appelant.
   *Recommandé.*
