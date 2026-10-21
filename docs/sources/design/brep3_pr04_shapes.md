# Palier 3 — PR 4 : maillage des formes, export, validation

| | |
|---|---|
| PR | [ssg-aero/gbs#183](https://github.com/ssg-aero/gbs/pull/183), branche `feat/brep-mesh-shapes` |
| Document de référence | [brep_palier3_mesh.md](brep_palier3_mesh.md), § 3, § 4.2, § 5, § 6, bilan au § 9 |
| Fichiers | `gbs-brep/mesh.h`, `gbs-io/mesh_writer.h` (nouveau), `gbs-render/vtkMeshRender.h` (nouveau), `python/gbsBindBrepOps.cpp`, `python/gbsbindrender.cpp`, `python/vistaplot.py`, `python/tests/test_brep_mesh.py` (nouveau), `tests/tests_brep_mesh.cpp`, `gbs-occt/tests/tests_brep3_occt.cpp` (nouveau), `docs/sources/design/brep_palier3_mesh.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Dernière PR du palier 3 :
- le maillage d'une forme entière ;
- l'export STL et OBJ ;
- Python et le rendu VTK ;
- la comparaison à `BRepMesh_IncrementalMesh` d'OpenCascade, et le bilan du palier.

La comparaison a conduit à réviser la mise à l'échelle des faces et la discrétisation des arêtes (§ 3).

## 2. Maillage d'une forme (`gbs-brep/mesh.h`)

```cpp
template <std::floating_point T>
struct Mesh
{
    struct FaceRange { FaceId face; std::uint32_t first, count; };
    std::vector<point<T, 3>> vertices;
    std::vector<std::array<std::uint32_t, 3>> triangles; // counter-clockwise seen from outside
    std::vector<FaceRange> faces;
    std::vector<point<T, 2>> uv;                         // per corner
    std::vector<point<T, 3>> normals;                    // per corner, unit, oriented
    std::vector<FaceId> incomplete;                      // refinement stopped at max_points
    std::vector<FaceId> folded;                          // a triangle against the surface normal at a corner
};
auto mesh(const Model<T> &, const ShapeId &, const MeshOptions<T> & = {}) -> BuildResult<Mesh<T>>;
T area(const Mesh<T> &);
T volume(const Mesh<T> &);
std::size_t open_sides(const Mesh<T> &);
```

- **Faces** :
  - chaque face de la forme est maillée une fois (`MeshBuilder::face`), les arêtes partagées ;
  - une face prend l'orientation de sa première utilisation dans une coque de la forme, `Forward` si elle n'est dans aucune.
  - `Reversed` retourne ses triangles (deux coins échangés) et ses normales.
- **Normales** : celle de la surface à chaque coin, unitaire et orientée.
  - Là où la surface n'en a pas (pôle, sommet : S_u × S_v petit devant max(|S_u|, |S_v|)²), c'est celle de la surface juste à l'intérieur du triangle, (u, v) déplacé de 1e-3 vers son centre de gravité : la limite de ce côté, radiale au pôle d'une sphère. À défaut, celle du triangle.
  - Une normale de la surface opposée à son triangle est gardée, et la face est signalée dans `folded` : c'est un défaut du maillage ou de l'orientation de la face dans sa coque, que les tests vérifient absent.
- **Sommets** : ceux qu'utilisent les triangles, dans leur ordre de création.
- **Erreurs** : la première face qui ne se maille pas rend son erreur.
- **Mesures** :
  - `area` et `volume` (somme des tétraèdres depuis l'origine) ;
  - `open_sides` compte les côtés qui ne sont pas utilisés une fois dans chaque sens : 0 pour un maillage fermé et orienté de façon cohérente.

Les options gagnent `min_angle` : l'angle minimal de 20° de la PR 3, réglable, 0 pour aucun.

## 3. Révisions apportées par la comparaison à OCCT

### 3.1 Mise à l'échelle par la courbure

Première mesure : les cylindres avaient 5 à 9 fois plus de triangles qu'avec OCCT.
- **Cause** : la mise à l'échelle par les seules longueurs moyennes de S_u et S_v fait des triangles isotropes en 3D. Leur taille est celle que demande la courbure autour du cylindre, alors que le long des génératrices aucune finesse n'est nécessaire.
- **Correction** : chaque direction est divisée en plus par l'espacement qu'y demandent les options :
  - min(√(8 flèche / k), angle / k, taille de la face) ;
  - k est la courbure normale moyenne le long de u ou de v.

Sur une sphère, rien ne change. Sur un cylindre, les triangles s'allongent le long des génératrices.

La courbure est estimée sans dérivées secondes : rotation des normales entre les échantillons d'une grille 9 × 9. On prend la **rotation totale sur la longueur totale**, et non la moyenne des rapports. Avec la moyenne des rapports, deux échantillons presque confondus au sommet d'un cône, aux normales différentes, donnaient une courbure de 1e23.

| Face, flèche 1e-3 (tests analytiques) | Triangles avant | Triangles après |
|---|---|---|
| cylindre R 1, h 2 | 5 410 | 256 |
| tore 3, 1 | 56 060 | 42 760 |
| sphère, cône | 13 582, 5 200 | 13 540, 4 948 |

### 3.2 Fusion des morceaux d'arête

Une courbe d'intersection approchée (cylindres croisés) a des centaines de nœuds, et chaque morceau de Bézier donnait un point. Les morceaux consécutifs sont maintenant fusionnés tant que les bornes de la discrétisation, prises sur leurs pôles réunis, le permettent :
- distance des pôles à la corde (enveloppe convexe) ;
- rotation du polygone de contrôle, jonctions comprises ;
- longueur de la corde.

Un morceau échantillonné (courbe qui n'est pas une B-spline) n'est pas fusionné.

La longueur d'une suite de morceaux fusionnés est doublée tant que les bornes tiennent, puis cherchée par dichotomie : O(k log k) tests de pôles pour une suite de k morceaux, au lieu de O(k²). Toute suite gardée satisfait les bornes.

### 3.3 Raffinement des arêtes sur leurs faces

La fusion a fait voir une erreur que la densité des nœuds masquait : sur une sphère découpée, la corde d'une pcurve s'écartait sur la surface de 1,37 fois la flèche. Le bord d'une face n'est jamais raffiné par sa triangulation (question 3 du palier). Chaque segment d'arête est donc coupé en deux tant que, sur une face de l'arête, la surface aux quarts ou au milieu de la corde de sa pcurve s'écarte de la corde 3D de plus que la flèche.

- **Écart du pcurve** : la pcurve ne suit la surface qu'à sa tolérance d'approximation près. Cet écart, mesuré aux bouts, aux quarts et au milieu du segment, est retranché, sans quoi la subdivision ne s'arrêtait pas.
- **Profondeur** : 16 au plus.
- **Faces de l'arête** : elles sont lues une fois sur tout le modèle.

## 4. Export (`gbs-io/mesh_writer.h`)

| Fonction | Contenu |
|---|---|
| `write_stl(mesh, binary = true, name)` | texte STL : binaire (en-tête de 80 octets « binary STL nom », jamais « solid » en tête, signature du format ASCII ; nombre ; 50 octets par triangle, simple précision, petit-boutiste quel que soit l'hôte) ou ASCII (`solid nom`) ; normale de facette, celle du triangle |
| `write_obj(mesh)` | OBJ : sommets une fois (`v`), normale de chaque coin (`vn`), un groupe par face (`g face<index>`), `f v//vn` |
| `write_stl_file`, `write_obj_file` | fichier ; `InvalidFile` s'il ne s'écrit pas |

- Les nombres sont écrits sous leur forme la plus courte qui se relit exactement, avec un point quelle que soit la locale du processus (`gbs-io/real_format.h`). C'est la fonction du writer STEP, sortie dans un en-tête commun : `std::to_chars`, ou avec libc++, où `to_chars` sur les flottants manque selon la version de macOS ciblée, `snprintf_l` dans la locale « C ».
- Les octets sont écrits par décalages, sans `std::byteswap`.

## 5. Python et VTK

- `gbs.brep` :
  - `MeshOptions(deflection, angle, max_edge, max_points, min_angle)` ;
  - `mesh(model, shape, options)` ;
  - `Mesh` :
    - `vertices`, `triangles`, `uv` et `normals` sont des tableaux numpy, copiés en un bloc ;
    - `faces`, `incomplete` ;
    - méthodes `area()`, `volume()`, `open_sides()` ;
  - `write_stl(mesh, path, binary)` et `write_obj(mesh, path)`, qui lèvent `BRepError`.
- `gbs-render/vtkMeshRender.h` :
  - `make_polydata(mesh)` donne trois points VTK par triangle, avec les normales des coins : l'ombrage suit les surfaces et les arêtes entre faces restent vives ;
  - `make_actor(mesh, couleur)` ;
  - en Python : `gbs.make_actor(mesh, col)`.
- `pygbs.vistaplot` : `mesh_to_polydata(mesh)` (PyVista) et `plot_mesh(mesh)`.

## 6. Tests

| Test | Cas |
|---|---|
| `tests_brep_mesh.mesh_of_shapes` | boîte, cylindre, sphère, cône, tore partiel lus en STEP :<br>- maillages fermés, volumes à l'aire fois la flèche de `signed_volume` ;<br>- normales unitaires du côté de leur triangle ;<br>- faces contiguës et complètes, tous les sommets utilisés.<br>Une face seule : 16 côtés ouverts, aire 1.<br>Une coque cousue (tournée vers l'intérieur) puis retournée : le volume suit le signe de la coque, puis change de signe. |
| `tests_brep_mesh.stl_and_obj` | STL binaire relu (taille, nombre, normales, sommets, volume en simple précision) ; STL ASCII relu (facettes, sommets, volume à 1e-9) ; OBJ (sommets, normales, faces, groupes) ; fichiers et chemin impossible (`InvalidFile`) |
| `tests_brep_mesh.vtk` | polydata : trois points par triangle, normales de la sphère en chaque point, pôles compris |
| `tests_brep_mesh.options_of_shapes` | `min_angle = 0` sur un cône : moins de triangles, fermé, aucune face repliée ; `max_points` 100 : face signalée, maillage fermé |
| `tests_brep_mesh.writers_header_and_locale` | STL binaire nommé « solid_blade » : l'en-tête ne commence pas par « solid » ; sous une locale à virgule (fr_FR ou de_DE, quand l'hôte en a une) : aucune virgule dans le STL ASCII et l'OBJ |
| `tests_brep_mesh.curve_parameters` (ajout) | cubique de 199 morceaux sur une droite : ses deux bouts seulement ; avec `max_edge` 0,25 : 13 segments au plus |
| `tests_brep3_occt.meshes` | 7 formes × 2 flèches contre `BRepMesh_IncrementalMesh` :<br>- les deux maillages fermés ;<br>- volume de gbs à l'aire fois la flèche du volume exact (`BRepGProp`) ;<br>- écart aux faces d'OCCT ≤ 1,05 flèche ;<br>- nombres de triangles à un facteur 6 près. |
| `test_brep_mesh.py` | options, tableaux, volume recalculé en numpy, `max_edge`, erreurs, STL relu par numpy, OBJ, acteur VTK, PyVista |

Les résultats contre OCCT sont au § 9.2 du document du palier. Les bornes de nombres de triangles de `curved_faces` sont mises aux nouvelles valeurs, plus basses.

En local passent :
- les 538 tests C++, dont les comparaisons OCCT ;
- les 58 tests Python.

## 7. Points à valider

1. **Mise à l'échelle par la courbure** (§ 3.1) : elle révise le § 4.2 du palier. Une face de courbure très inégale reste raffinée isotropiquement par endroits (limite au § 9.5 du palier).
2. **Arêtes raffinées sur leurs faces** (§ 3.3) : la finesse d'une arête dépend maintenant aussi des surfaces qui la portent, ce qui reste partagé, donc étanche.
3. **Angle minimal de 20° par défaut** : il double à triple les triangles des cylindres et des cônes par rapport à OCCT (424 contre 212, 1 574 contre 511). `min_angle = 0` donne 252 et 1 260.
4. **`Mesh` par coin** pour `uv` et `normals` (3 par triangle), comme le document du palier le prévoyait. Les sommets restent partagés.

### 7.1 Revue de code

| Commentaire | Suite donnée |
|---|---|
| `snprintf` suit la locale du processus : sous une locale à virgule, STL ASCII et OBJ illisibles | forme la plus courte et exacte, avec un point quelle que soit la locale : la fonction du writer STEP sortie dans `gbs-io/real_format.h`, partagée ; testé sous fr_FR |
| la normale du triangle remplaçait aussi celle de la surface quand elle s'y opposait : un triangle replié passait | repli seulement quand S_u × S_v s'annule ; au pôle, la normale de la surface juste à l'intérieur du triangle ; une opposition est signalée (`Mesh::folded`) et vérifiée absente dans les tests et la comparaison OCCT |
| en-tête binaire qui peut commencer par « solid » | « binary STL » + nom ; testé avec « solid_blade » |
| raffinement des arêtes sur leurs faces : un seul échantillon par corde | quarts et milieu de la corde de la pcurve, écart de la pcurve mesuré aux mêmes paramètres |
| fusion des morceaux en O(k²) | longueur doublée puis dichotomie : O(k log k) |
| cas non couverts | morceaux fusionnés (199 morceaux sur une droite : 2 points), `min_angle = 0`, `incomplete` par `mesh()`, en-tête « solid », locale à virgule |

