# Écriture STEP — PR 3 : topologie, pièces, assemblages

| | |
|---|---|
| PR | branche `feat/step-write-topology` |
| Document de référence | [step_writer.md](step_writer.md), § 4, 5, 7, 8 ; § 12.2 (instances) |
| Fichiers | `gbs-io/step/writer.h` (nouveau), `tests/tests_step_write.cpp` (nouveau), `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Écrire des formes du modèle dans un fichier STEP complet :
- la topologie, avec ses pcurves ;
- les pièces et leurs produits ;
- les assemblages, avec leurs pièces partagées et leurs placements ;
- l'API `write_step`.

C'est le premier **aller-retour complet** gbs → STEP → gbs.

## 2. API

```cpp
struct StepWriteOptions { schema = "AP214" | "AP242"; approx_tol = 1e-5; pcurves = true; flatten = false;
                          optional<double> uncertainty; product_name, author, organization; };
struct StepWriteReport  { entities; parts, assemblies, occurrences; solids, shells, faces, edges, vertices;
                          approximated_curves, approximated_surfaces, max_deviation; skipped; };

auto write_step(model, shapes | shape, opts) -> BuildResult<pair<string, StepWriteReport>>;
auto write_step_file(model, shapes | shape, path, opts) -> BuildResult<StepWriteReport>;
```

## 3. Choix d'architecture et justification

### 3.1 Topologie

| gbs | STEP |
|---|---|
| sommet | `vertex_point` |
| arête | `edge_curve(v1, v2, géométrie, .T.)` ; la géométrie est la courbe restreinte à l'arête (PR 2) ; avec ses pcurves, une `surface_curve`, ou une `seam_curve` si une face l'utilise deux fois |
| pcurve d'une co-arête | `pcurve(surface, definitional_representation(courbe 2D, contexte 2D))` ; pour un seam, la pcurve de l'usage `Forward` d'abord, comme OCCT |
| co-arête | `oriented_edge` ; les **arêtes dégénérées sont omises** |
| wire | `edge_loop` dans l'ordre des co-arêtes ; `face_outer_bound` pour `wires[0]`, `face_bound` pour les trous |
| face et son usage | `advanced_face(bornes, surface, same_sense = usage Forward)` ; bornes d'orientation `same_sense` (§ 4.1 de la conception) ; une face utilisée dans les deux sens est écrite une fois, puis en `oriented_face` |
| shell | `closed_shell` s'il est fermé, sinon `open_shell` |
| solide | `manifold_solid_brep` ; avec cavités, `brep_with_voids` dont les cavités sont écrites tournées vers l'extérieur, puis utilisées à rebours (`oriented_closed_shell(.F.)`) |

Sommets, arêtes et surfaces sont mis en cache par pièce, et donc partagés comme
dans le modèle. Les surfaces sont écrites en premier : une pcurve référence la
surface de sa face.

### 3.2 Pièces

| Forme | Représentation |
|---|---|
| solides | `advanced_brep_shape_representation` (solides et repère d'origine) |
| shells, faces | `manifold_surface_shape_representation` (`shell_based_surface_model`) ; une face seule dans un `open_shell` |
| les deux | `shape_representation` liée à une représentation de chaque sorte (`shape_representation_relationship`), comme le lit la PR 9 de la lecture |

Chaque pièce est un produit complet :
- `product`, `product_definition_formation` et `product_definition` ;
- `product_definition_shape` et `shape_definition_representation` ;
- les contextes d'application du schéma choisi ;
- une `product_related_product_category` en fin de fichier.

Le nom du produit est celui de la forme. À défaut, c'est
`StepWriteOptions::product_name`, et sinon `part<n>` ou `assembly<n>`.

### 3.3 Assemblages et instances partagées

Un compound devient un assemblage, sous la forme d'une `shape_representation`
de repères. Pour chaque membre :
- **Pièce ou sous-assemblage.** Le membre est écrit **une fois par clé de
  pièce** (`Model::instance`, PR préliminaire P2), et réutilisé si une autre
  instance a la même clé et la même signature (nombres de faces, d'arêtes et
  de sommets). Sans provenance, le membre est une pièce à lui seul.
- **Repère de la pièce.** Une instance déplacée y est ramenée par l'inverse de
  son placement. Le calcul se fait sur une **copie de travail** du modèle (un
  `append` fait une seule fois, puis `copy_shape`), jamais sur le modèle
  lui-même, qui reste `const`.
- **Occurrence.** Elle comprend :
  - une `next_assembly_usage_occurrence` ;
  - une `context_dependent_shape_representation` portant une
    `representation_relationship_with_transformation` (rep_1 = enfant, rep_2 =
    parent, comme OCCT) ;
  - une `item_defined_transformation` du repère d'origine de l'enfant vers son
    placement dans le parent. Ce placement vaut `P_parent⁻¹ ∘ P_membre`, P
    étant les placements de provenance (l'identité sans provenance).
- **Option `flatten`** : une seule pièce contenant toutes les formes, placées
  dans le repère du modèle.
- **Wires, arêtes et sommets libres** : omis et listés dans `skipped`
  (question 8).

### 3.4 Erreurs : tout ou rien

Les erreurs suivantes rendent une `BuildError` sans rien écrire :
- forme morte ;
- schéma inconnu ;
- face sans surface ;
- dossier non inscriptible pour `write_step_file`.

L'écrivain ne modifie jamais le modèle.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| `write_step(m, span, opts)` | plus une surcharge pour une forme seule | usage le plus courant |
| rapport | plus `parts`, `assemblies`, `occurrences` | suivre le partage des pièces |
| — | copie de travail du modèle pour les instances déplacées | le modèle passé est `const` |

## 5. Tests (`tests_step_write`, 5 cas)

| Test | Couvre |
|---|---|
| `round_trip_of_the_test_files` | les **14 fichiers** de `tests/in/step` (boîte, cylindre, sphère, cône, plaque en pouces, carreau B-spline, bande, calotte, cône sans seam, bande trouée, tore partiel et loft OCCT, assemblage, `mapped_item`) : lus, écrits, relus **sans échec**, même nombre de solides, faces, arêtes et sommets, même volume à 1e-6 près, même nom, même unité, **aucune approximation** |
| `analytic_entities_and_pcurves` | cylindre : 1 `cylindrical_surface`, 2 `plane`, 2 `circle`, 1 `seam_curve`, 2 `surface_curve`, 6 `pcurve` ; sans pcurves ; en-tête AP242 ; schéma refusé |
| `assemblies_share_their_parts` | ASM { BOX, BOX, SUB { CYL } } : 4 produits, 4 occurrences, **2 solides écrits** (BOX partagée) ; relu : mêmes noms, mêmes positions, BOX de nouveau partagée ; `flatten` : une pièce, trois solides |
| `shapes_built_in_gbs` | boîte cousue (faces B-spline), sphère et tore en faces naturelles sur primitives (seams), surface donnée par une formule approchée et comptée, sommet omis ; relu : volumes |
| `units_files_and_errors` | modèle en pouces relu en millimètres (facteur 25,4) ; écriture de fichier ; forme morte, dossier inexistant |

Toutes les suites BREP, STEP et IGES passent.

## 6. Points à valider

1. **Pcurves écrites** en `surface_curve` / `seam_curve`, usage `Forward`
   d'abord pour un seam.
2. **Sens des faces** par `same_sense` et l'orientation des bornes.
3. **Cavités** écrites tournées vers l'extérieur et utilisées à rebours.
4. **Pièces partagées par clé et signature**, ramenées dans leur repère sur une
   copie de travail.
5. **Placement d'occurrence** `P_parent⁻¹ ∘ P_membre`, rep_1 = enfant.
6. **Tout ou rien** pour les erreurs ; formes libres omises et rapportées.
