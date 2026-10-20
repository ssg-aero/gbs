# Noyau BREP natif — conception du palier 3 : maillage

| | |
|---|---|
| Statut | **À valider** |
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
| Espace de la triangulation | l'espace (u, v) de la face, **mis à l'échelle** par les longueurs moyennes de S_u et S_v pour que la triangulation y soit à peu près isotrope en 3D | triangulation 3D directe d'une surface |
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

1. Départ : les nœuds de la courbe (pour une B-spline) et quelques points
   réguliers.
2. Subdivision de chaque intervalle tant que la flèche au milieu (et aux quarts)
   dépasse `deflection`, que l'angle entre les tangentes aux bouts dépasse
   `angle`, ou que la corde dépasse `max_edge`.
3. Les sommets du maillage aux bouts sont ceux des sommets du modèle.

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

La triangulation se fait dans (a u, b v), avec a et b les longueurs moyennes de
S_u et S_v sur la face (échantillonnées). Un cylindre de rayon 1 et de hauteur
100 donne ainsi des triangles à peu près isotropes en 3D, au lieu de triangles
très allongés le long de la génératrice.

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
ne sont pas raffinés : leurs voisins couvrent la surface autour du pôle.

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
