# Noyau BREP natif — conception du palier 3 : maillage

| | |
|---|---|
| Statut | **Réalisé** (bilan au § 9) |
| Périmètre | Maillage en triangles des faces, coques et solides du BREP natif, étanche le long des arêtes partagées, à une flèche et un angle donnés ; export STL et OBJ ; Python ; rendu VTK ; comparaison à `BRepMesh_IncrementalMesh` d'OpenCascade |
| Dépend de | palier 1 (modèle, pcurves, coutures, arêtes dégénérées), palier 2 et 2 bis (booléens, dont on maille les résultats) |
| Hors périmètre | maillage volumique (tétraèdres), quadrangles, remaillage d'un maillage existant, réparation de modèles, maillage adapté à un calcul (éléments finis) |

Ce document décrit *quoi* construire et *pourquoi*, comme
[brep_palier2.md](brep_palier2.md) et [brep_palier2bis.md](brep_palier2bis.md).
Chaque PR a sa note d'architecture `brep3_prNN_*.md`.

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Triangulation | **Delaunay contrainte écrite dans le dépôt** (`gbs-mesh/cdt.h`), sur une structure indexée : Bowyer-Watson incrémental, récupération des arêtes contraintes par retournements, retrait des triangles hors de la face ; prédicats **exacts** de Shewchuk (`inc/topology/robustPredicates.h`) | bibliothèque CDT en en-têtes (une dépendance de plus) ; la Delaunay en demi-arêtes de `inc/topology` (`tessellations.h`, § 4.3) |
| Précision | **flèche** (écart entre le maillage et la géométrie) et **angle** entre normales de triangles voisins, comme OCCT ; longueur maximale facultative | nombre de subdivisions fixe |
| Étanchéité | chaque arête est discrétisée **une fois**, ses paramètres partagés par toutes les faces qui l'utilisent (même `t`, donc mêmes points 3D, mêmes sommets du maillage) | coudre a posteriori des maillages de faces indépendants |
| Espace de la triangulation | l'espace (u, v) de la face, **mis à l'échelle** dans chaque direction par la finesse que la courbure, la flèche et l'angle y demandent (révisé à la PR 4 : d'abord les seules longueurs de S_u et S_v) | triangulation 3D directe d'une surface |
| Raffinement | insertion de points (centre circonscrit dans l'espace mis à l'échelle) tant qu'un triangle s'écarte de la surface de plus que la flèche ou plie de plus que l'angle | grille régulière du rectangle des paramètres |
| Coutures, pôles | une couture est deux bords de la triangulation (u = u1 et u = u2) dont les points sont **les mêmes sommets** ; le côté dégénéré d'un pôle est un bord de la triangulation dont tous les points sont **un seul sommet** | maillage spécial des surfaces fermées |
| Résultat | `Mesh` : sommets 3D, triangles, et pour chaque face ses triangles et (u, v) ; normales par sommet et par face | structure en demi-arêtes (`gbs-mesh/halfedgemesh.h`, ébauche) |

---

## 1. Objectifs et non-objectifs

### 1.1 Objectifs

1. **Mailler** en triangles une face, une coque, un solide ou un composé du
   BREP natif, y compris les faces découpées (trous, contours quelconques), les
   surfaces fermées (couture) et les pôles (sphère, cône).
2. **Étanchéité** : deux faces qui partagent une arête partagent les sommets du
   maillage le long de celle-ci ; le maillage d'une coque fermée est fermé
   (chaque arête de triangle utilisée deux fois, en sens opposés).
3. **Précision contrôlée** par la flèche et l'angle, et vérifiée.
4. **Orientation** : les triangles suivent l'orientation des faces dans leur
   coque (normale vers l'extérieur d'un solide).
5. **Export** STL (ASCII et binaire) et OBJ ; **Python** ; **rendu VTK**.
6. **Validation** contre OpenCascade (`BRepMesh_IncrementalMesh`) et par des
   invariants (fermeture, caractéristique d'Euler, volume et aire).

### 1.2 Non-objectifs

- Maillage volumique, quadrangles, maillage pour le calcul (qualité d'éléments
  au-delà d'un angle minimal raisonnable).
- Remaillage ou simplification d'un maillage.
- Réparation : un modèle invalide (pcurves absentes, contours ouverts) est
  refusé avec une erreur, comme les autres opérations du noyau.

---

## 2. Paramètres

| Paramètre | Sens | Défaut |
|---|---|---|
| `deflection` | écart maximal entre le maillage et la géométrie (flèche) | 1e-3 × diagonale de la boîte de la forme |
| `angle` | angle maximal (radians) entre les tangentes de deux segments consécutifs d'une arête, et entre les normales de la surface aux sommets d'un triangle | 0,5 |
| `max_edge` | longueur maximale d'un segment ou d'un côté de triangle (0 : aucune) | 0 |
| `relative` | la flèche est relative à la taille de chaque arête ou face | non |

La flèche d'un segment est la distance maximale de la courbe à la corde ; celle
d'un triangle, la distance maximale de la surface (aux points du triangle dans
l'espace (u, v)) au plan du triangle.

---

## 3. Discrétisation des arêtes

Chaque arête du modèle (pas chaque co-arête) est discrétisée **une fois** :
une suite de paramètres t0 = u1 < t1 < … < tn = u2 de sa courbe 3D.

1. Départ : les morceaux de Bézier de la courbe (pour une B-spline), échantillonnés
   sinon.
2. Subdivision de chaque morceau tant que ses pôles peuvent s'écarter de la corde
   de plus que `deflection`, que son polygone de contrôle tourne de plus que
   `angle`, ou que la corde dépasse `max_edge`.
3. Fusion des morceaux consécutifs tant que ces bornes, sur leurs pôles réunis,
   le permettent (PR 4) : une courbe d'intersection approchée, aux centaines de
   nœuds, n'a pas plus de points que sa forme n'en demande.
4. Subdivision de chaque segment tant que, sur une face de l'arête, la surface aux
   quarts ou au milieu de la corde de la pcurve s'écarte de la corde 3D de plus que
   la flèche (PR 4) : le bord d'une face n'étant jamais raffiné par sa triangulation, les
   cordes des pcurves doivent déjà suivre la surface.
5. Les sommets du maillage aux bouts sont ceux des sommets du modèle.

Les pcurves ont le paramètre de l'arête (`same_parameter`) : les points (u, v)
d'une co-arête sont ses pcurves aux mêmes t. Une arête dégénérée (pôle) n'a pas
de courbe 3D : elle est discrétisée sur sa pcurve, assez finement pour suivre le
bord de la triangulation, et tous ses points sont **un seul sommet**.

---

## 4. Triangulation d'une face

### 4.1 Contours dans l'espace des paramètres

Les contours de la face (le premier extérieur, les suivants des trous) sont des
polygones de (u, v), suite des points des co-arêtes dans l'ordre et le sens du
contour. Une couture apparaît deux fois dans le contour, à u1 et à u2 : deux
bords distincts de la triangulation, mais les mêmes sommets du maillage.

### 4.2 Mise à l'échelle

La triangulation se fait dans (a u, b v), avec a la longueur moyenne de S_u
divisée par l'espacement que les options demandent le long de u (de même b pour
v). Cet espacement est celui qu'il faut à la corde d'une courbe de la courbure
normale moyenne le long de u pour rester dans la flèche (√(8 flèche / k)) et
tourner de moins que l'angle (angle / k), au plus la taille de la face. La
courbure le long de u est la rotation totale des normales sur la longueur totale
entre les échantillons d'une grille 9 × 9 : un pôle, où les normales tournent sur
une longueur nulle, n'y pèse rien.

- Une sphère : a et b dans le rapport des longueurs de S_u et S_v, triangles
  isotropes en 3D (comme la première version).
- Un cylindre : pas de courbure le long des génératrices, triangles allongés le
  long de celles-ci, comme ceux d'OCCT.

La première version (longueurs seules) faisait des triangles isotropes de la
taille que la courbure autour du cylindre demande : 1 830 triangles pour un
cylindre R 5 h 10 à 1e-3 de la diagonale, contre 212 pour OCCT ; 424 maintenant
(§ 9.2).

### 4.3 Delaunay contrainte (`gbs-mesh/cdt.h`)

Le dépôt a déjà une Delaunay en demi-arêtes (`inc/topology/tessellations.h`,
avec la récupération de Sloan dans `edgeRecovery.h`). Elle n'est pas reprise :
ses éléments sont liés par des pointeurs partagés et son insertion parcourt
toutes les faces (6,8 s pour 20 000 points aléatoires, contre 0,05 s pour la
structure indexée ci-dessous), et sa version contrainte ne traite ni les trous
ni les contraintes intérieures. Ses prédicats exacts (`robustPredicates.h`) le
sont : `orient2d`, et `incircle` ajouté à côté.

Générique, en 2D, sans dépendance au BREP :

1. **Bowyer-Watson** incrémental dans un super-triangle : chaque point inséré
   retire les triangles dont le cercle circonscrit le contient et retriangule la
   cavité. Prédicats d'orientation et du cercle **exacts** (Shewchuk) ; points
   à moins de 1 000 arrondis de la boîte fusionnés.
2. **Récupération des contraintes** : chaque segment de contour absent de la
   triangulation est rétabli en retournant les arêtes qui le coupent.
3. **Retrait** : les triangles hors du domaine (propagation depuis le
   super-triangle sans traverser les contraintes, parité des contours pour les
   trous) sont retirés.

### 4.4 Raffinement

Tant qu'un triangle viole un critère, un point est inséré à son centre
circonscrit (dans l'espace mis à l'échelle). Quand ce centre est hors de la face
ou empiète un segment du contour, le contour n'est pas modifié : un point de
contour inséré serait un point de l'arête, partagé par l'autre face, ce que § 3
exclut. Le triangle reçoit alors son centre de gravité, vers l'intérieur. Les
points sont insérés un à un, chacun tant que son triangle existe (Ruppert).

| Critère | Mesure |
|---|---|
| flèche | distance de S(u, v) au triangle, au centre et aux milieux des côtés intérieurs |
| angle | angle entre les normales de la surface aux trois sommets (une normale nulle, à un pôle, ne compte pas) |
| longueur | côté 3D intérieur au-delà de `max_edge` |
| repli | normale du triangle opposée à celle de la surface (près d'un sommet de cône) |
| qualité | angle minimal sous ~20° (Ruppert), pour éviter des triangles dégénérés, dans l'espace mis à l'échelle globalement (celui de la Delaunay ; une métrique locale ne converge pas près des pôles, voir la note de la PR 3) |

Le raffinement est borné (nombre de points par face, `max_points`) ; au-delà, la face est
maillée au mieux et le rapport le signale.

### 4.5 Pôles et sommets de cône

Le côté dégénéré est un bord de la triangulation dans (u, v) : ses points sont
un seul sommet 3D. Les triangles qui y touchent par un côté entier (deux
sommets confondus) sont retirés du maillage final : ils sont d'aire nulle. Ils
ne sont pas raffinés : leurs voisins couvrent la surface autour du pôle. Un
triangle qui n'a qu'un coin au pôle est mesuré là où il est en 3D, ce coin
prenant le paramètre libre des deux autres (§ 10).

---

## 5. Assemblage et orientation

Le maillage d'une forme rassemble les faces : les sommets des arêtes partagés,
les sommets intérieurs propres à chaque face. Les triangles d'une face utilisée
`Reversed` dans sa coque sont retournés. Normales par sommet : celle de la
surface en (u, v) (une par face pour un sommet d'arête, la normale de la surface
étant discontinue d'une face à l'autre).

```cpp
template <std::floating_point T>
struct Mesh
{
    std::vector<point<T, 3>> vertices;
    std::vector<std::array<std::uint32_t, 3>> triangles;
    struct FaceRange { FaceId face; std::uint32_t first, count; };
    std::vector<FaceRange> faces;                 // triangles of each face, in order
    std::vector<point<T, 2>> uv;                  // per corner of each triangle: (u, v) on its face
    std::vector<point<T, 3>> normals;             // per corner: surface normal, oriented
    std::vector<FaceId> incomplete;               // faces whose refinement stopped at max_points
};

auto mesh(const Model<T> &, const ShapeId &, const MeshOptions<T> &) -> BuildResult<Mesh<T>>;
```

---

## 6. Validation

| Contrôle | Attendu |
|---|---|
| fermeture | maillage d'une coque fermée : chaque arête de triangle en deux triangles, sens opposés |
| Euler | V − E + F = 2 − 2g pour une coque fermée de genre g (sphère, tore) |
| flèche | distance des points de la géométrie (échantillons des faces) au maillage ≤ `deflection` |
| volume, aire | convergence vers le volume et l'aire exacts quand la flèche diminue (ordre 2) |
| orientation | volume signé du maillage d'un solide positif |
| OCCT | mêmes formes maillées par `BRepMesh_IncrementalMesh` à la même flèche : volumes et aires proches, nos maillages fermés ; nombres de triangles du même ordre |

Formes de test : primitives (boîte, cylindre, sphère, cône, tore), faces
découpées et trouées, résultats des booléens du palier 2 et 2 bis (comparaison
OCCT existante), formes STEP lues.

---

## 7. Plan des PR

| # | PR | Contenu | Estimation |
|---|---|---|---|
| 1 | arêtes (avec ce document et le bilan du palier 2 bis) | § 2, § 3 : options, discrétisation partagée des arêtes, contours (u, v) des faces avec coutures et pôles | 8 h |
| 2 | Delaunay contrainte | § 4.3 : `gbs-mesh/cdt.h`, générique, tests 2D (Delaunay, contraintes, trous) | 12 h |
| 3 | faces | § 4.1, § 4.2, § 4.4, § 4.5 : triangulation et raffinement d'une face, coutures et pôles | 12 h |
| 4 | formes, export, validation | § 5, § 6 : assemblage étanche, orientation, normales, STL et OBJ, Python, VTK, comparaison OCCT, bilan | 10 h |
| | | **Total** | **≈ 42 h** |

---

## 8. Questions à valider

1. **Delaunay écrite dans le dépôt** (décidé en préparant ce document).
   *Recommandé.*
2. **Flèche et angle comme OCCT**, flèche par défaut relative à la taille de la
   forme. *Recommandé.*
3. **Un point de contour n'est jamais inséré par le raffinement d'une face** :
   la finesse d'une arête est fixée par § 3, pour toutes ses faces. Une face
   très courbe près d'une arête peu courbe est raffinée vers l'intérieur
   seulement. *Recommandé.*
4. **Triangles d'aire nulle aux pôles retirés** du maillage final. *Recommandé.*
5. **Pas de demi-arêtes** : le `Mesh` est une liste de triangles indexés ; une
   structure d'adjacence pourra être construite à côté si un usage le demande.
   *Recommandé.*

---

## 9. Bilan de réalisation

### 9.1 PR et notes d'architecture

| PR | Note | Contenu |
|---|---|---|
| [ssg-aero/gbs#178](https://github.com/ssg-aero/gbs/pull/178) | [brep3_pr01_edges.md](brep3_pr01_edges.md) | options, discrétisation partagée des arêtes, contours (u, v) des faces |
| [ssg-aero/gbs#179](https://github.com/ssg-aero/gbs/pull/179) | [brep3_pr02_cdt.md](brep3_pr02_cdt.md) | Delaunay contrainte indexée, prédicats exacts |
| [ssg-aero/gbs#182](https://github.com/ssg-aero/gbs/pull/182) | [brep3_pr03_faces.md](brep3_pr03_faces.md) | triangulation et raffinement d'une face, coutures, pôles |
| [ssg-aero/gbs#183](https://github.com/ssg-aero/gbs/pull/183) | [brep3_pr04_shapes.md](brep3_pr04_shapes.md) | maillage d'une forme, STL et OBJ, Python, VTK, comparaison OCCT |

### 9.2 Validation par OpenCascade

`gbs-occt/tests/tests_brep3_occt.cpp` : formes construites par OCCT et lues par
gbs (STEP), et booléens de gbs sur elles, maillés par gbs et par
`BRepMesh_IncrementalMesh` à la même flèche absolue et au même angle (0,5) :

- les deux maillages sont fermés ;
- l'écart à la forme est mesuré par `BRepExtrema` aux centres de gravité de
  300 triangles.

Flèche de 1e-3 de la diagonale :

| Forme | gbs : triangles | gbs : écart / flèche | OCCT : triangles | OCCT : écart / flèche |
|---|---|---|---|---|
| boîte | 12 | 0 | 12 | 0 |
| cylindre R 5 h 10 | 424 | 0,31 | 212 | 0,43 |
| sphère R 3 | 3 796 | 0,98 | 2 836 | 0,76 |
| cône à sommet | 1 574 | 0,22 | 511 | 0,57 |
| tore 3, 1 | 3 700 | 0,93 | 3 510 | **1,69** |
| boîte moins sphère (booléen gbs) | 846 | 0,68 | 837 | 0,54 |
| cylindres croisés (booléen gbs) | 444 | 0,67 | 406 | 0,43 |

À 1e-2 de la diagonale, OCCT dépasse aussi la flèche sur la sphère (1,41) ;
gbs reste en dessous partout (0,98 au plus). Les nombres de triangles sont du
même ordre ; gbs en fait plus sur le cylindre et le cône à cause de l'angle
minimal de 20° (sans lui : 252 et 1 260), que l'option `min_angle` règle.

### 9.3 Décisions révisées en cours de route

- **Mise à l'échelle par la courbure** (§ 4.2), au lieu des seules longueurs :
  la comparaison à OCCT a montré des triangles isotropes, de 5 à 9 fois plus
  nombreux sur les cylindres.
- **Fusion des morceaux d'arête** et **raffinement des arêtes sur leurs faces**
  (§ 3) : une courbe d'intersection approchée donnait autant de points que de
  nœuds ; la fusion a fait voir que les cordes des pcurves s'écartaient de la
  surface (1,37 flèche sur une sphère découpée), que le raffinement sur les
  faces corrige.
- **Qualité réglable** (`MeshOptions::min_angle`, 20° par défaut) et **limite de
  points** (`max_points`) en options.
- **Qualité jugée dans l'espace mis à l'échelle globalement** : la métrique
  locale de chaque triangle ne converge pas près des pôles (PR 3).
- **Pas de taille minimale** : une règle « côtés plus courts que la flèche »
  annulait le critère d'angle quand la flèche est lâche (PR 3).

### 9.4 Défauts trouvés en chemin

- Un triangle pouvait se replier près d'un sommet de cône ; les normales
  nulles aux pôles relançaient le raffinement sans fin (PR 3).
- Insérer tous les centres d'une passe d'un coup faisait boucler le raffinement
  de qualité (PR 3).
- Une estimation de courbure moyennée par échantillon explosait près d'un
  sommet de cône (normales différentes à distance nulle) : la rotation totale
  sur la longueur totale la remplace (PR 4).
- Une pcurve ne suit la surface qu'à sa tolérance d'approximation près : le
  raffinement d'une arête sur ses faces retranche cet écart, sans quoi il ne
  s'arrêtait pas (PR 4).

### 9.5 Limites restantes

- **Écart échantillonné** : le critère de flèche teste quatre points par
  triangle ; l'écart mesuré dépasse la flèche de 4 % au plus sur les faces de
  test.
- **Une seule mise à l'échelle par face** : une face dont la courbure varie
  beaucoup (un congé qui s'évase, une B-spline très inégale) est raffinée
  isotropiquement là où sa courbure diffère de la moyenne ; une Delaunay
  anisotrope (métrique par point) le règlerait.
- **Pôles** : les triangles le long d'un pôle sont retirés ; la finesse autour du
  pôle vient des points intérieurs.
- **Pas de maillage volumique**, ni de quadrangles, ni de remaillage (hors
  périmètre).

## 10. Revue de code (a posteriori, lot 3 de la revue du noyau BREP)

Trouvé pendant le lot 3 ([ssg-aero/gbs#196](https://github.com/ssg-aero/gbs/pull/196),
[brep_pr05_face_wire.md](brep_pr05_face_wire.md) § 8) : corriger l'écart des
pcurves projetées déplaçait de 2e-7 en relatif la flèche par défaut de
`cone.stp`, et son maillage passait de 1 948 triangles à 362 000. Le défaut a
été reproduit par une sonde qui balaie la flèche ; le test cité
(`tests_brep_mesh.deflection_sweep`) échoue sur la version précédente. Ce qui
suit remplace la dernière règle du § 2.4 de [brep3_pr03_faces.md](brep3_pr03_faces.md)
(« le critère de repli le raffine »).

| Commentaire | Suite donnée |
|---|---|
| 1. Maillage instable au sommet d'un cône : `cone.stp` à sa flèche par défaut d0 = 6,524e-3 donne 1 948 triangles en 3 ms ; à d0 × (1 + 2e-7), 196 414 triangles en 16 s ; à d0 × 0,999, 70 582 ; faces incomplètes et repliées. La sonde (1 162 flèches de 0,1 à 10 d0, dont 102 à ±5e-7 de d0, limite de 20 000 points) en trouve 170 mauvaises : 31 sur 102 à ±5e-7 de d0, 26 sur 100 ailleurs à ±1 % | **cause**, en trois temps. (a) Un triangle qui a un coin sur le sommet est en 3D le triangle (sommet, B, C), quel que soit le u de ce coin : tous les points du côté dégénéré sont le même sommet. Il était mesuré dans (u, v) : la surface au milieu (u, v) du côté (sommet, B) est sur la génératrice à mi-chemin entre le u du coin et celui de B, pas sur celle de B, où est le côté 3D. Au coin u = 2π, au bout de la couture, les triangles minces le long de la couture avaient ainsi un écart r Δu que leur triangle 3D n'a pas. (b) Leur centre circonscrit étant au-delà de la couture, ils recevaient leur centre de gravité, tiré vers ce coin, au-delà des u de B et de C : le triangle qu'il fait au sommet avec le plus proche des deux est replié, car l'orientation 3D d'un triangle au sommet ne dépend que de l'ordre des u de ses deux autres coins. (c) Aucun point inséré dans un triangle replié au sommet ne le déplie : des deux triangles qu'il fait au sommet, l'un garde l'ordre faux. Chaque centre de gravité en repliait un nouveau plus près du coin, jusqu'à `max_points`. **Corrigé** : un triangle au pôle est mesuré là où il est en 3D (`mesh_detail::combination`) : le coin au pôle prend le paramètre libre des deux autres, le milieu du côté (sommet, B) est sur la génératrice de B, le centre de gravité entre celles de B et de C. De même pour le critère de repli, la normale au sommet dans `mesh` et l'écart échantillonné des tests (l'image du triangle (u, v) n'est pas le triangle 3D). `FaceMesh::free` donne le paramètre libre de chaque point (u au pôle d'une sphère et au sommet d'un cône, −1 ailleurs) |
| 2. Le point forcé d'un triangle au pôle (centre circonscrit hors de la face ou empiétant son bord) | **corrigé** : sur sa médiane depuis le pôle, son paramètre libre entre ceux des deux autres coins ; c'est le centre de gravité quand le coin au pôle est entre eux. L'éventail autour du pôle garde l'ordre de ses génératrices. Un triangle déjà replié garde son centre de gravité (commentaire 3). Sur une sphère limitée de 200 à 3 000 points (281 limites), les maillages ouverts ou repliés à l'arrêt passent de 5 à 2, tous deux sous 250 points |
| 3. Essayé puis écarté : ne plus raffiner un triangle replié au pôle, puisque aucun point inséré dans ce triangle ne le déplie | écarté : sur la même sphère, 149 maillages sur 281 restent ouverts ou repliés, contre 2. Les triangles repliés au pôle sont fréquents en cours de raffinement (23 à 200 points), et ce sont les cavités des points insérés autour, plus larges que le triangle, qui les déplient. La limite `max_points` borne seule le raffinement, comme avant |
| Résultats | sonde : aucune flèche mauvaise sur les 1 162 de `cone.stp`, `sphere.stp`, `cylinder.stp` et `occt_partial_torus.stp` (170 sur le cône avant) ; à ±10 % de d0, le cône a de 1 770 à 2 218 triangles. L'écart exact (formules du cône et de la sphère, dans les deux sens) ne change pas : 0,87 flèche au plus sur le cône, 1,046 sur la sphère (les 4 % du § 9.5). Les nombres de triangles sont ceux d'avant ou moins, de 0,64 à 1,01 fois ceux d'avant de 0,1 à 10 d0 ; à 1e-3 et 20°, le cône en a 15 542 au lieu de 12 500, le hasard du raffinement de qualité (11 940 au lieu de 15 722 à 0,158 d0). Sans `min_angle`, la face d'un cône est l'éventail exact de son sommet à son cercle de base : 126 triangles au lieu de 1 810 à d0, 32 et 128 au lieu de 414 et 4 948 pour le cône de `curved_faces`. Testé (`deflection_sweep` : cône, 103 flèches à ±10 % et ±2e-7 de d0, 1 782 à 2 184 triangles ; sphère, 23 flèches, 3 468 à 4 116). `options_of_shapes` atteint `max_points` sur la sphère : le cône n'a plus besoin de point intérieur. Temps inchangés (deux passes alternées, RelWithDebInfo, clang 21, avec la correction des pcurves) : `tests_brep_boolean` 25,1 et 25,3 s contre 24,9 et 25,1 s, `tests_step_write` 10,2 et 10,0 s contre 10,0 et 10,0 s, les autres suites `tests_brep_*`, `tests_step_*` et `tests_intersect_*` à 0,05 s près ; `tests_brep_mesh` passe de 0,81 s à 1,26 s : 0,47 s pour `deflection_sweep`, 0,82 s pour le reste |
| Hors du lot, trouvés par la sonde, identiques sur la version précédente | non traités : avec `angle` = 0,35 et une flèche lâche (au-delà de 1,8 d0), le raffinement s'arrête sans point à insérer et laisse un triangle replié, 23 flèches sur 362 pour `occt_partial_torus.stp` (sans pôle) et 7 pour `sphere.stp` ; avec `min_angle` = 30°, au-delà des 20,7° qui garantissent l'arrêt de Ruppert, 2 à 3 flèches sur 362 pour `cone.stp` |
