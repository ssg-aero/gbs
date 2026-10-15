# Palier 2 bis — PR 1 : contacts tangents

| | |
|---|---|
| PR | branche `feat/brep2b-tangent-contacts` |
| Document de référence | [brep_palier2bis.md](brep_palier2bis.md), § 3 (avec ce document) |
| Fichiers | `gbs-brep/boolean.h`, `gbs-brep/classify.h`, `gbs-brep/section.h`, `tests/tests_brep_boolean.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/design/brep_palier2bis.md` (nouveau), `docs/sources/design/brep_palier2.md`, `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Les booléens acceptent les solides qui se touchent en un point ou le long
d'une ligne, au lieu de rendre `NonGeneric`. Ils ne laissent plus d'entités
inutilisées dans le modèle.

## 2. Changements

| Fichier | Changement |
|---|---|
| `boolean.h` | les points de contact et les arêtes tangentes ne sont plus refusés ; une arête tangente peut rester pendante dans une face qu'elle touche sans la traverser |
| `boolean.h` | `pair_radially`, `separate_non_manifold_edges` : à une arête gardée par plus de deux faces, les faces sont appariées autour de l'arête et chaque paire sauf la première reçoit une copie de l'arête (§ 3.3 du document) |
| `boolean.h` | `erase_unused` : les arêtes et sommets créés par l'opération et inutilisés sont effacés, en cas de succès comme d'échec ; les fils et le composé de la section aussi en cas d'échec |
| `boolean.h` | `erase_shape` liste les faces, shells, solides et composés **avant** d'en effacer : les faces non retenues des copies restaient dans le modèle |
| `classify.h` | `interior_points` : un point intérieur par ligne de v constant, plus longues cordes d'abord ; `classify(face, solide)` prend le premier qui n'est pas sur le bord du solide |
| `section.h` | un morceau de courbe proche du bord d'une face n'est une arête de bord que si l'intersection courbe-courbe a trouvé un recouvrement à cet endroit ; sinon c'est un contact ponctuel |

### 2.1 Le contact au centre d'une face

La sphère posée sur la boîte la touche au centre du dessus, exactement au
point intérieur que choisit `interior_point`. Ce point était classé `On`, et
la face jetée. Les autres lignes de v constant donnent d'autres points, et le
premier qui n'est pas sur le bord classe la face.

### 2.2 Le cercle tangent à une arête

Le cylindre couché sur la boîte coupe les faces latérales y = 0 et y = 3
selon des cercles tangents à l'arête du dessus. Restreint à la face, chaque
cercle se réduit à un arc minuscule autour du point de tangence. Cet arc est
`On` (à la tolérance de l'arête) et était rapporté comme une arête posée sur
la face, d'où un refus. Le test de recouvrement le ramène à un contact.

## 3. Tests

| Test | Cas | Attendu |
|---|---|---|
| `tangent_contacts` | sphère R 3 posée sur la boîte 2 × 3 × 4 | fusion : 2 solides, 24 + 36 π ; commun vide ; différence : la boîte ; seules les arêtes du résultat restent dans le modèle |
| | cylindre R 5 couché sur la boîte, tangent le long de x = 1, z = 4 | fusion : 2 solides ; commun vide ; différence : 24 |
| | boîte 3..5 × −1.5..1.5 × 3..7 dont la face x = 5 est tangente au cylindre (tourné pour éloigner le seam) | différence : **deux copeaux** qui se touchent le long de la génératrice, chacun avec sa propre arête ; commun et fusion : volumes calculés en forme fermée |
| `non_generic` | boîte tournée de 45° posée sur une arête du dessus de l'autre (PR 2) | `NonGeneric`, **aucune arête ni aucun sommet laissés** par la tentative |
| OCCT `booleans` | `sphere_touching`, `cylinder_lying`, `tangent_slivers` × 3 opérations | valides, volumes à 2e-3, **mêmes nombres de solides et de faces** que `BRepAlgoAPI` |

Les suites `brep`, `intersect`, `step`, `iges` et le binaire OCCT passent.

## 4. Points à valider

1. Questions 1, 2 et 6 du document (solides séparés, arêtes pendantes
   tangentes seulement, nettoyage automatique).
2. **Pas de l'appariement radial** : h = 1 % de la corde de l'arête. Assez
   pour séparer un plan et un cylindre de rayon 1 (≈ 5e-3 rad), assez petit
   pour rester dans des faces étroites.
3. Le plan des PR et les estimations du § 10.
