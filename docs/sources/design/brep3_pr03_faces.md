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

   Les passes s'arrêtent quand plus aucun point n'est inséré, ou à `MeshOptions::max_points` points par face (200 000 par défaut). `complete` dit alors si des défauts de géométrie restent.

   Le contour et le raffinement sont dans un même `try` : une erreur de la triangulation (contour qui se croise, marche qui ne finit pas) donne `NonGeneric`, et les points intérieurs déjà créés sont retirés des sommets du maillage. `face()` ne lève donc pas d'exception.
5. **Résultat** :
   - les triangles du domaine, sauf ceux qui ont deux coins au même sommet (le long d'un pôle : aire nulle, § 4.5) ;
   - les points intérieurs sont de nouveaux sommets du maillage, S(u, v).

### 2.2 Critères d'un triangle

Ils sont testés dans l'ordre suivant.

| Critère | Mesure | Défaut |
|---|---|---|
| pôle | deux coins au même sommet | aucun : retiré du maillage, ses voisins couvrent la surface |
| repli | normale du triangle 3D contre S_u × S_v au centre de gravité | géométrie |
| longueur | côté intérieur (non contraint) au-delà de `max_edge` | géométrie |
| flèche | distance de S aux milieux des côtés intérieurs et au centre de gravité, au triangle 3D | géométrie |
| angle | angle entre les normales de la surface aux coins, au-delà de `angle` | géométrie |
| qualité | plus petit angle sous 20° dans l'espace mis à l'échelle | qualité |

- Les milieux des côtés du contour ne sont pas testés. Le contour n'est jamais raffiné (question 3 du palier) : une flèche sur un côté du contour ne se corrigerait pas, et le raffinement ne finirait pas.
- La distance au triangle tient compte d'un triangle dégénéré : la distance aux côtés quand il n'a pas d'aire.
- Un triangle jugé bon le reste, puisque ses sommets ne bougent pas. Il est gardé dans une table et n'est plus évalué aux passes suivantes ; la clé commence par son plus petit sommet, pour qu'un même triangle ait une seule clé.
- La surface au milieu d'un côté intérieur est gardée par côté : les deux triangles du côté la partagent.
- La recherche d'un centre circonscrit dans la triangulation part d'un coin de son triangle, et non du dernier point inséré.

### 2.3 Métrique de la qualité et des centres

La qualité (20°) et les centres circonscrits sont jugés dans l'espace mis à l'échelle **globalement** (§ 4.2 du palier), celui où la triangulation est de Delaunay. Près d'un pôle, où |S_u| tend vers 0, un triangle de 20° y est plus étroit en 3D ; la flèche, l'angle et le repli, mesurés en 3D, restent garantis.

La métrique **locale** (|S_u| et |S_v| au centre de gravité du triangle) a été essayée, avec les mêmes critères :

| Face | Métrique globale | Métrique locale |
|---|---|---|
| sphère, flèche 1e-3 | 13 582 triangles, 0,03 s | 399 818 triangles, limite atteinte (`complete` faux) |
| cône, flèche 1e-3 | 5 200 triangles | 399 842 triangles, limite atteinte |
| cylindre, tore | identiques à 3 % près | identiques à 3 % près |

La triangulation reste de Delaunay dans l'espace global : près d'un pôle, elle refait des triangles allongés en u que la métrique locale juge mauvais, et le raffinement ne s'arrête plus. Juger la qualité en local demanderait une Delaunay anisotrope (métrique par point), hors du périmètre de ce palier.

### 2.4 Pôles et sommets de cône

Le long d'un pôle, la normale S_u × S_v s'annule : son calcul n'y donne que du bruit (dérivées par différences, ou arrondis).
- Une normale de longueur inférieure à 1e-6 fois le produit des échelles est **inconnue** et ne compte pas pour l'angle. Sans cela, le critère d'angle raffinait sans fin vers le pôle.
- Les triangles à deux coins sur le pôle ne sont pas raffinés : ils sont retirés. Raffinés, ils insèrent leur centre de gravité toujours plus près du pôle.
- Près d'un sommet de cône, un triangle (sommet, x, y) dont le coin au sommet n'est pas entre x et y en u peut se **replier** en 3D : son orientation 3D dépend seulement de l'ordre des u de x et de y, l'orientation dans (u, v) aussi de leurs v. Le critère de repli le raffine.

Ces deux règles suffisent à arrêter le raffinement près des pôles. Une règle de taille minimale (pas de raffinement pour un triangle aux côtés plus courts que la flèche), essayée d'abord, ne change aucun maillage des tests ; elle est retirée, car elle annulait le critère d'angle quand la flèche est lâche.

## 3. Ajouts à `ConstrainedDelaunay` (`gbs-mesh/cdt.h`)

Le raffinement doit connaître le domaine avant `remove_outside`, qui interdit ensuite toute insertion.

| Méthode | Rôle |
|---|---|
| `domain()` | les triangles du domaine que `remove_outside` laisserait |
| `place(p, near)` | `outside` (hors du domaine, sur une contrainte, hors de la boîte, non fini), `inside`, ou `encroaching` : dans le domaine, mais dans le cercle de diamètre une arête contrainte du triangle qui le contient (empiètement de Ruppert) ; la marche part du point `near` s'il est donné |
| `has_triangle(a, b, c)` | le triangle existe encore |

La classification (parcours 0-1 par parité, déjà utilisé par `remove_outside`) est gardée par triangle :
- une insertion la conserve : la cavité ne traverse pas de contrainte, donc ses nouveaux triangles sont dans la région du triangle où tombe le point ;
- une contrainte ajoutée, ou un point qui en coupe une, l'invalide ; elle est alors recalculée à la demande.

## 4. Tests

| Test | Cas |
|---|---|
| `tests_mesh_cdt.domain_before_removal` | carré troué : aire du domaine 12 ; places dans le trou, dehors, hors de la boîte, sur le contour, en empiètement, non finie ; 300 points tirés, ceux à l'intérieur insérés : le domaine garde son aire, et `remove_outside` laisse les mêmes triangles |
| `tests_brep_mesh.planar_faces` | carré : 2 triangles, χ = 1 ; `max_edge` 0,2 : tous les côtés ≤ 0,2, points sur la surface ; carré troué d'un cercle : aire 4 − π/4, χ = 0, rien dans le trou |
| `tests_brep_mesh.curved_faces` | cylindre analytique et NURBS, sphère, tore, cône, demi-révolution, chacun aux flèches 1e-2 et 1e-3 : `complete`, χ attendu (0, 2, 0, 1, 1), fermé pour la sphère et le tore (chaque côté une fois dans chaque sens), écart échantillonné ≤ 1,05 flèche, nombre de triangles au plus 25 % au-dessus de celui mesuré, aire à 2 flèches relatives, chaque triangle le long de S_u × S_v |
| `tests_brep_mesh.faces_of_a_solid` | les faces d'un cube, orientées par leur coque : maillage fermé, χ = 2 |
| `tests_brep_mesh.faces_of_a_boolean` | boîte 0..4 moins la sphère R 3 en son coin (faces découpées, pôle et couture à l'intérieur des faces) : maillage fermé, χ = 2, écart ≤ 1,05 flèche, volume 64 − 9π/2 à l'aire fois la flèche près |
| `tests_brep_mesh.face_options_and_failures` | l'angle seul (cylindre, flèche 1, angle 0,1) : normales des coins à moins de 0,1 ; `max_edge` 0,2 sur un tore : tous les côtés ; `max_points` 300 sur une sphère : `complete` faux, maillage fermé, χ = 2 ; contour qui se croise dans (u, v) : `NonGeneric`, aucun sommet ajouté |

| Face, flèche 1e-3 | Triangles | Écart mesuré |
|---|---|---|
| cylindre R 1, h 2 | 5 410 | 0,9998e-3 |
| sphère R 1 | 13 582 | 1,028e-3 |
| tore 3, 1 | 56 060 | 1,038e-3 |
| cône | 5 200 | 0,9996e-3 |

L'écart mesuré dépasse la flèche de 4 % au plus : le critère n'échantillonne que quatre points par triangle. Toute la suite `tests_brep_mesh` prend 0,5 s.

| Face, flèche 1e-4 | Triangles | Avant la revue | Après |
|---|---|---|---|
| sphère | 135 466 | 0,44 s | 0,23 s |
| cylindre | 53 610 | 0,11 s | 0,07 s |
| tore | 398 462 (limite de 200 000 points) | 1,78 s | 0,78 s |

Les 512 tests C++ passent.

### 4.1 Revue de code

| Commentaire | Suite donnée |
|---|---|
| les exceptions de la triangulation sortent de `face()`, qui rend un `BuildResult` | contour et raffinement dans un même `try` (`triangulate`, privée) : `NonGeneric`, et les points intérieurs créés retirés ; testé avec un contour qui se croise |
| coût de chaque passe | clé du cache tournée sur le plus petit sommet, milieux des côtés en cache par côté, marche de `place` depuis un coin du triangle : deux fois plus rapide à 1e-4 (§ 4) ; la liste de travail n'a pas été faite, le parcours du domaine étant une petite part du temps une fois les triangles bons en cache |
| qualité jugée dans l'espace mis à l'échelle globalement | métrique locale essayée et mesurée : le raffinement ne s'arrête plus sur la sphère et le cône ; la métrique globale est gardée, la raison est en § 2.3 |
| borne de 1,5 flèche dans les tests, 4 % dans la note | tests à 1,05 flèche (écart maximal mesuré 1,038) |
| chemins non couverts, tailles non bornées | angle seul, `max_edge` sur un tore, limite de points (devenue l'option `max_points`), contour qui se croise ; nombres de triangles bornés ; la règle de taille minimale, qui annulait l'angle quand la flèche est lâche, est retirée |

## 5. Points à valider

1. **Contour jamais raffiné** (question 3 du palier) : les milieux des côtés du contour ne sont pas testés pour la flèche. Une face très courbe le long d'une arête peu courbe (une arête droite sur un cylindre, le long de la génératrice) garde la finesse de l'arête au bord.
2. **Écart échantillonné en quatre points** : jusqu'à 4 % au-delà de la flèche. Plus de points coûteraient autant d'évaluations de la surface.
3. **Métrique globale** pour la qualité et les centres circonscrits (§ 2.3) : près d'un pôle, des triangles de 20° dans (u, v) mis à l'échelle sont étroits en 3D.
4. **Discrétisation des arêtes dégénérées** (heuristique de la PR 1 : un tour complet découpé selon l'angle) : gardée. Les triangles le long du pôle étant retirés, la finesse autour du pôle vient des points intérieurs, et le nombre de points du pôle compte peu.
5. **Orientation** : les triangles sont dans le sens de S_u × S_v. Le retournement des faces `Reversed` dans leur coque est fait à l'assemblage (PR 4) ; les tests le font déjà.
