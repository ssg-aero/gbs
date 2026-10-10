# Lecture STEP — PR 6 : topologie, lecteur complet

| | |
|---|---|
| PR | branche `feat/step-topology` |
| Document de référence | [step_reader.md](step_reader.md), § 3.4, 5, 6, 8 et 9 ; questions 8 et 9 |
| Fichiers | `gbs-io/step/reader.h` (nouveau), `gbs-brep/builders.h` (deux codes d'erreur), `tests/tests_step_read.cpp` (nouveau), `tests/in/step/*.stp` et `make_test_files.py` (nouveaux), `tests/CMakeLists.txt`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Assembler les PR 1 à 5 en un lecteur : un fichier STEP donne un
`gbs::brep::Model<double>` avec ses solides, shells ouverts, faces, arêtes et
sommets, plus un rapport de lecture. La topologie du fichier est **reprise
telle quelle**, sans sewing.

## 2. API

```cpp
#include <gbs-io/step/reader.h>

gbs::brep::Model<double> m;
auto r = gbs::step::read_step_file(m, "part.stp");            // ou read_step(m, contenu_texte)
if (!r) { r.error().code; r.error().message; }                 // BuildError
r->root;                                                       // la forme lue, ou un compound de plusieurs
r->report;                                                     // StepReadReport

struct StepReadOptions {
  std::optional<double> target_unit_mm;  // défaut : l'unité du modèle (mm pour un modèle neuf)
  bool strict = false;
  MakeFaceOptions<double> face;          // tol relevée à l'incertitude du fichier
  double max_tolerance = 1e-2;           // écart toléré sommet / courbe et arête / surface (unités cibles)
  std::size_t edge_samples = 16;         // points par arête pour borner les surfaces
};

struct StepReadReport {
  schemas, length_unit, angle_unit, length_factor, file_tolerance,
  solids, shells, faces, edges, vertices,
  unsupported, failed,                   // StepIssue { #id, type, message }
  max_tolerance, raised_tolerances, sense_mismatches,
  check                                  // check() du résultat
};
```

Deux codes d'erreur s'ajoutent à `BuildErrc` : `InvalidFile` (fichier illisible,
mal formé ou sans forme) et `UnsupportedEntity` (entité non prise en charge, en
mode strict).

## 3. Choix d'architecture et justification

### 3.1 Des représentations aux formes

Le lecteur parcourt les représentations de forme (`advanced_brep_shape_representation`,
`manifold_surface_shape_representation`, `shape_representation`,
`faceted_brep_shape_representation`). Chacune donne :

- son **contexte**, d'où viennent les unités, avec un `GeometryReader` par
  contexte ;
- le **nom du produit** qui la porte : `shape_definition_representation` →
  `product_definition_shape` → … → `product.name`, directement ou à travers
  une `shape_representation_relationship`.

| Élément de la représentation | Traduction |
|---|---|
| `manifold_solid_brep`, `brep_with_voids` | un solide |
| `shell_based_surface_model` | son shell (un compound s'il en a plusieurs) |
| `mapped_item` (assemblages, PR 9), `faceted_brep` (question 13) | rapportés non pris en charge |
| repère, autre élément | ignoré |

Les racines qu'aucune représentation ne cite sont lues dans le premier
contexte du fichier. Une forme unique est rendue telle quelle, plusieurs dans
un compound.

### 3.2 Sommets et arêtes : partagés par `#id`

- **Sommet** : `vertex_point` → `make_vertex`. Sa tolérance est l'incertitude
  du fichier, au moins `brep_default_tolerance`.
- **Arête** : `edge_curve(v1, v2, courbe, same_sense)`.
  - Les paramètres des sommets sont obtenus par **inversion analytique** sur la
    définition STEP de la courbe (PR 4). Pour une conique, la plage peut
    chevaucher 0 ; pour une arête fermée (v1 = v2), elle couvre un tour entier.
  - Si `same_sense = .F.`, ou si la courbe est une `trimmed_curve` parcourue à
    rebours, la NURBS est **retournée** pour respecter la convention du palier
    1 : `u` croît de `v1` vers `v2`.
  - La tolérance de l'arête est relevée pour contenir ses sommets. Un écart
    supérieur à `max_tolerance` fait échouer l'arête. Les relèvements
    au-delà de l'incertitude du fichier sont comptés (`raised_tolerances`).
- Sommets et arêtes sont **mis en cache par `#id`** : deux faces voisines
  partagent la même arête, comme dans le fichier. Une arête en échec n'est
  tentée qu'une fois.

### 3.3 Faces

1. **Boucles.** Chaque `edge_loop` donne ses co-arêtes dans l'ordre et le sens
   de ses `oriented_edge`. Une borne (`face_bound`) d'orientation `.F.`
   retourne la boucle. Une `vertex_loop` ne crée pas de wire, mais son sommet
   compte pour borner la surface.
2. **Surface.** Les arêtes des boucles sont échantillonnées dans le sens de
   parcours, puis `parameter_box` (PR 4) borne la surface et `to_nurbs` la
   convertit. **Chaque face a sa propre NURBS**, même si plusieurs faces citent
   la même surface STEP (recommandation de la question 4).
3. **Construction** par `make_wire_ordered` et `make_face_use` (PR 5).
   - La boucle de la `face_outer_bound` est l'extérieure. Sans
     `face_outer_bound`, chaque boucle est essayée tour à tour comme
     extérieure.
   - Une face sans boucle d'arêtes (sphère ou tore complet) est une face
     naturelle sur toute la surface, orientée par `same_sense`.
4. **Précision.** Si une arête s'écarte de sa surface de plus que
   `face.pcurve_tol` (`EdgeOffSurface`, `PCurveApproximation`), la face est
   reconstruite avec une précision dix fois plus large, jusqu'à
   `max_tolerance`. Les fichiers industriels ont des écarts bien supérieurs à
   l'incertitude annoncée.
5. **Sens.** Le `FaceUse` rendu par `make_face_use` garde les sens du fichier.
   On le compare à celui qu'annonce `same_sense`, et les désaccords sont
   comptés (`sense_mismatches`). Ils sont attendus pour une face bornée par ses
   seuls seams, dont le sens est fixé par `make_solid`.

### 3.4 Shells et solides

- `closed_shell` / `open_shell` : un `Shell` des usages de faces.
  `oriented_closed_shell` et `oriented_face` retournent les usages si leur
  orientation est `.F.`.
- `manifold_solid_brep` : `make_solid`, qui vérifie le volume et tourne le
  shell vers l'extérieur si besoin. `brep_with_voids` : avec ses cavités.
- **Noms et identifiants.** Le `#id` de chaque entité topologique lue va dans
  `externalId`, son nom non vide dans `name`. Le nom du produit va à la forme
  racine quand elle n'en a pas.

### 3.5 Import partiel par défaut, strict sur option (question 8)

Ce qui suit applique la recommandation de la question 8 :

| Situation | Défaut | `strict = true` |
|---|---|---|
| entité non prise en charge (surface inconnue, `mapped_item`…) | rapportée dans `unsupported`, la face qui en dépend dans `failed` | échec `UnsupportedEntity` |
| face non construite (arête, géométrie, builder) | omise et rapportée ; son shell reste ouvert ; le solide n'est pas créé et **son shell est rendu à sa place** | échec avec le code du builder |
| erreur de syntaxe, fichier absent, aucune forme, unité d'un modèle non vide différente de la cible | échec `InvalidFile` | idem |

En cas d'échec, le **modèle est rendu inchangé** : les entités créées sont
effacées (pierres tombales, que `compact()` retire) et l'unité est restaurée.
Le rapport contient toujours le `check()` du résultat.

### 3.6 Unités et modèle

L'unité cible est par défaut **celle du modèle** :

- pour un modèle neuf, le millimètre (question 6) ;
- `target_unit_mm` la fixe sur un modèle vide ;
- lire dans un modèle non vide d'une autre unité est refusé, comme `append`
  (PR 2).

`Model::unitScale` est renseigné.

### 3.7 AP242 (question 9)

Le schéma n'est pas filtré : la partie BREP d'un fichier AP242 est lue comme
celle d'AP203 ou AP214, et le reste est ignoré (recommandation de la question
9).

## 4. Écarts par rapport au document de conception

| Document (§ 9) | Implémentation | Raison |
|---|---|---|
| `read_step(m, path)` et `read_step(m, string_view)` | `read_step_file(m, path)` et `read_step(m, contenu)` | une chaîne littérale se convertit vers les deux types : appel ambigu |
| `target_unit_mm = 1.` | `std::optional`, défaut : unité du modèle | lire dans un modèle existant sans préciser l'unité |
| `assemblies = true` | absent | les assemblages arrivent en PR 9 ; les `mapped_item` sont rapportés |
| `StepReadReport::failed_faces`, `max_tolerance_raised` | `failed` (faces et solides), `raised_tolerances`, `max_tolerance`, `sense_mismatches` | rapport plus complet |
| fichiers de test « écrits à la main » | engendrés par `tests/in/step/make_test_files.py`, versionnés avec lui | lisibles, reproductibles et modifiables ; le script documente leur contenu |

## 5. Tests (`tests_step_read`, 6 cas, 107 assertions)

Les fichiers de `tests/in/step/` sont lus par leur chemin absolu
(`GBS_TESTS_IN_DIR`, défini pour tous les tests), quel que soit le répertoire
courant.

| Fichier | Contenu | Vérifié |
|---|---|---|
| `box.stp` | boîte 2 × 3 × 4 mm, trois plans d'axe rentrant (`same_sense .F.`), arêtes parcourues à rebours, une face nommée | solide, 6/12/8, volume 24, nom du produit, nom de face, `#id`, aucun désaccord de sens, `check()` |
| `cylinder.stp` | cylindre R 5, h 10, `seam_curve` utilisée deux fois | volume 250π, `check()` |
| `sphere.stp` | sphère bornée par son seul seam | deux arêtes dégénérées insérées, volume 36π |
| `cone.stp` | cône, angles en **degrés**, sommet sans arête | une arête dégénérée, volume πR²H/3 |
| `plate_inch.stp` | plaque trouée en **pouces**, borne du trou `.F.`, `shell_based_surface_model` | shell ouvert, aires en mm², lecture en pouces |
| `bspline_patch.stp` | quart de cylindre en B-spline rationnelle (instance complexe), arcs rationnels | rayon 1 sur toute la face |
| `box_unsupported_face.stp` | boîte dont une face repose sur une surface inconnue | défaut : 5 faces, shell ouvert rendu, rapport ; strict : `UnsupportedEntity`, modèle inchangé |
| texte en mémoire | erreur de syntaxe, aucune forme, fichier absent, unité de modèle différente, deux lectures dans le même modèle | |

Toutes les suites BREP et STEP passent.

## 6. Points à valider

1. **API** `read_step` / `read_step_file`, unité cible par défaut celle du
   modèle.
2. **Topologie partagée par `#id`**, arêtes retournées pour respecter
   `u1 → u2 = v1 → v2`, tolérances relevées et comptées.
3. **Une NURBS par face**, bornée par ses bords.
4. **Échelle de précision** des pcurves (×10 jusqu'à `max_tolerance`).
5. **Import partiel** : shell rendu à la place d'un solide incomplet ; **strict**
   : modèle inchangé (recommandation de la question 8).
6. **AP242 lu sans filtrage** (recommandation de la question 9).
7. **Fichiers de test engendrés par script**, script versionné.
