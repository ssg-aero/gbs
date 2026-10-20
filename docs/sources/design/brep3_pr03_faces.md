# Palier 3 — PR 3 : triangulation des faces

| | |
|---|---|
| PR | branche `feat/brep-mesh-faces` |
| Document de référence | [brep_palier3_mesh.md](brep_palier3_mesh.md), § 4.1, § 4.2, § 4.4, § 4.5 |
| Fichiers | `gbs-brep/mesh.h`, `gbs-mesh/cdt.h`, `tests/tests_brep_mesh.cpp`, `tests/tests_mesh_cdt.cpp`, `docs/sources/design/brep_palier3_mesh.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Mailler une face du BREP natif en triangles, dans son espace de paramètres :
- son contour (PR 1) est contraint dans la Delaunay de la PR 2 ;
- la triangulation est raffinée à l'intérieur jusqu'à la flèche, l'angle et la longueur demandés ;
- les coutures et les pôles referment le maillage.

L'assemblage des faces en un maillage de forme, avec l'orientation des faces dans leur coque, l'export et la comparaison à OCCT, est l'objet de la PR 4.

## 2. `MeshBuilder::face(f)` (`gbs-brep/mesh.h`)

```cpp
template <std::floating_point T>
struct FaceMesh
{
    std::vector<point<T, 2>> uv;                         // its points, in (u, v)
    std::vector<std::uint32_t> vertices;                 // the vertex of the mesh at each point
    std::vector<std::array<std::uint32_t, 3>> triangles; // counter-clockwise in (u, v): along S_u x S_v
    bool complete{true};                                 // false: the refinement stopped at the limit
};
auto MeshBuilder<T>::face(FaceId) -> BuildResult<FaceMesh<T>>;
```

### 2.1 Étapes

1. **Contour** : les boucles de `boundary(f)` (coutures deux fois, pôles en un seul sommet).
2. **Mise à l'échelle** : l'espace (u, v) est multiplié par les longueurs moyennes de S_u et S_v, échantillonnées sur 9 × 9 points de la boîte du contour. Une longueur nulle (une surface dégénérée dans une direction) prend celle de l'autre.
3. **Delaunay contrainte** :
   - les points du contour sont insérés, ses segments contraints ;
   - un contour qui ne se triangule pas (segments croisés) donne `NonGeneric`.
4. **Raffinement**, par passes. Un triangle du domaine qui n'est pas encore jugé bon reçoit :
   - un point à son centre circonscrit, si ce centre est dans la face et loin de son contour ;
   - à son centre de gravité sinon, mais seulement pour un défaut de géométrie. Un défaut de qualité seul ne reçoit rien.

   Un point n'est inséré que si son triangle existe encore : deux centres voisins insérés dans la même passe créaient des triangles aplatis sans fin. C'est l'algorithme séquentiel de Ruppert, qui garantit l'arrêt du raffinement de qualité.

   Les passes s'arrêtent quand plus aucun point n'est inséré, ou à 200 000 points par face. `complete` dit alors si des défauts de géométrie restent.
5. **Résultat** :
   - les triangles du domaine, sauf ceux qui ont deux coins au même sommet (le long d'un pôle : aire nulle, § 4.5) ;
   - les points intérieurs sont de nouveaux sommets du maillage, S(u, v).

### 2.2 Critères d'un triangle

Ils sont testés dans l'ordre suivant.

| Critère | Mesure | Défaut |
|---|---|---|
| pôle | deux coins au même sommet | aucun : retiré du maillage, ses voisins couvrent la surface |
| repli | normale du triangle 3D contre S_u × S_v au centre de gravité | géométrie |
| taille | tous ses côtés 3D plus courts que la flèche | aucun : il est dans la flèche |
| longueur | côté intérieur (non contraint) au-delà de `max_edge` | géométrie |
| flèche | distance de S aux milieux des côtés intérieurs et au centre de gravité, au triangle 3D | géométrie |
| angle | angle entre les normales de la surface aux coins, au-delà de `angle` | géométrie |
| qualité | plus petit angle sous 20° dans l'espace mis à l'échelle | qualité |

- Les milieux des côtés du contour ne sont pas testés. Le contour n'est jamais raffiné (question 3 du palier) : une flèche sur un côté du contour ne se corrigerait pas, et le raffinement ne finirait pas.
- La distance au triangle tient compte d'un triangle dégénéré : la distance aux côtés quand il n'a pas d'aire.
- Un triangle jugé bon le reste, puisque ses sommets ne bougent pas. Il est gardé dans une table et n'est plus évalué aux passes suivantes.

### 2.3 Pôles et sommets de cône

Le long d'un pôle, la normale S_u × S_v s'annule : son calcul n'y donne que du bruit (dérivées par différences, ou arrondis).
- Une normale de longueur inférieure à 1e-6 fois le produit des échelles est **inconnue** et ne compte pas pour l'angle. Sans cela, le critère d'angle raffinait sans fin vers le pôle.
- Les triangles à deux coins sur le pôle ne sont pas raffinés : ils sont retirés. Raffinés, ils insèrent leur centre de gravité toujours plus près du pôle.
- Près d'un sommet de cône, un triangle (sommet, x, y) dont le coin au sommet n'est pas entre x et y en u peut se **replier** en 3D : son orientation 3D dépend seulement de l'ordre des u de x et de y, l'orientation dans (u, v) aussi de leurs v. Le critère de repli le raffine, avant la règle de taille.

## 3. Ajouts à `ConstrainedDelaunay` (`gbs-mesh/cdt.h`)

Le raffinement doit connaître le domaine avant `remove_outside`, qui interdit ensuite toute insertion.

| Méthode | Rôle |
|---|---|
| `domain()` | les triangles du domaine que `remove_outside` laisserait |
| `place(p)` | `outside` (hors du domaine, sur une contrainte, hors de la boîte, non fini), `inside`, ou `encroaching` : dans le domaine, mais dans le cercle de diamètre une arête contrainte du triangle qui le contient (empiètement de Ruppert) |
| `has_triangle(a, b, c)` | le triangle existe encore |

La classification (parcours 0-1 par parité, déjà utilisé par `remove_outside`) est gardée par triangle :
- une insertion la conserve : la cavité ne traverse pas de contrainte, donc ses nouveaux triangles sont dans la région du triangle où tombe le point ;
- une contrainte ajoutée, ou un point qui en coupe une, l'invalide ; elle est alors recalculée à la demande.

## 4. Tests

| Test | Cas |
|---|---|
| `tests_mesh_cdt.domain_before_removal` | carré troué : aire du domaine 12 ; places dans le trou, dehors, hors de la boîte, sur le contour, en empiètement, non finie ; 300 points tirés, ceux à l'intérieur insérés : le domaine garde son aire, et `remove_outside` laisse les mêmes triangles |
| `tests_brep_mesh.planar_faces` | carré : 2 triangles, χ = 1 ; `max_edge` 0,2 : tous les côtés ≤ 0,2, points sur la surface ; carré troué d'un cercle : aire 4 − π/4, χ = 0, rien dans le trou |
| `tests_brep_mesh.curved_faces` | cylindre analytique et NURBS, sphère, tore, cône, demi-révolution, chacun aux flèches 1e-2 et 1e-3 : `complete`, χ attendu (0, 2, 0, 1, 1), fermé pour la sphère et le tore (chaque côté une fois dans chaque sens), écart échantillonné ≤ 1,5 flèche, aire à 2 flèches relatives, chaque triangle le long de S_u × S_v |
| `tests_brep_mesh.faces_of_a_solid` | les faces d'un cube, orientées par leur coque : maillage fermé, χ = 2 |
| `tests_brep_mesh.faces_of_a_boolean` | boîte 0..4 moins la sphère R 3 en son coin (faces découpées, pôle et couture à l'intérieur des faces) : maillage fermé, χ = 2, écart ≤ 1,5 flèche, volume 64 − 9π/2 à l'aire fois la flèche près |

| Face, flèche 1e-3 | Triangles | Écart mesuré |
|---|---|---|
| cylindre R 1, h 2 | 5 456 | 0,99899e-3 |
| sphère R 1 | 13 564 | 1,028e-3 |
| tore 3, 1 | 56 082 | 1,038e-3 |
| cône | 5 200 | 0,9996e-3 |

L'écart mesuré dépasse la flèche de 4 % au plus : le critère n'échantillonne que quatre points par triangle. Toute la suite `tests_brep_mesh` prend 0,6 s.

Les 511 tests C++ passent.

## 5. Points à valider

1. **Contour jamais raffiné** (question 3 du palier) : les milieux des côtés du contour ne sont pas testés pour la flèche. Une face très courbe le long d'une arête peu courbe (une arête droite sur un cylindre, le long de la génératrice) garde la finesse de l'arête au bord.
2. **Écart échantillonné en quatre points** : jusqu'à 4 % au-delà de la flèche. Plus de points coûteraient autant d'évaluations de la surface.
3. **Taille minimale** : un triangle dont tous les côtés sont plus courts que la flèche n'est pas raffiné, pour l'angle non plus (comme la taille minimale d'OCCT). Un congé de rayon inférieur à la flèche n'est pas suivi en angle.
4. **Discrétisation des arêtes dégénérées** (heuristique de la PR 1 : un tour complet découpé selon l'angle) : gardée. Les triangles le long du pôle étant retirés, la finesse autour du pôle vient des points intérieurs, et le nombre de points du pôle compte peu.
5. **Orientation** : les triangles sont dans le sens de S_u × S_v. Le retournement des faces `Reversed` dans leur coque est fait à l'assemblage (PR 4) ; les tests le font déjà.
