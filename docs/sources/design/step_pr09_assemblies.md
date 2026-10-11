# Lecture STEP — PR 9 : assemblages

| | |
|---|---|
| PR | branche `feat/step-assemblies` |
| Document de référence | [step_reader.md](step_reader.md), § 7 et question 7 (décidée : assemblages **aplatis**, géométrie dupliquée par instance) |
| Fichiers | `gbs-brep/edit.h` (`Placement`, `copy_shape`), `gbs-io/step/reader.h`, `tests/tests_brep_edit.cpp`, `tests/tests_step_read.cpp`, `tests/in/step/assembly.stp`, `mapped_items.stp` et générateur, `gbs-occt/tests/tests_step_occt.cpp`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Lire les fichiers d'assemblage : produits, occurrences placées, sous-assemblages.
Conformément à la question 7, le modèle n'a pas de placement. **Chaque
occurrence devient une copie placée** de la forme de son composant, et chaque
assemblage un **compound nommé** d'après son produit.

## 2. API

```cpp
#include <gbs-brep/edit.h>

struct Placement<T> {                  // p -> R p + t, rigide
  static Placement frame(o, x, y, z);  // repère local -> global
  point apply(p); Placement operator*(b) /* b d'abord */; Placement inverse(); bool is_identity();
};
auto copy_shape(Model<T> &m, const ShapeId &shape, const Placement<T> &p) -> ShapeId;
```

Le lecteur n'a pas d'API nouvelle : `read_step` rend la racine de l'assemblage
(un compound nommé) ou, sans assemblage, la forme comme avant.

## 3. Choix d'architecture et justification

### 3.1 `copy_shape` dans le BREP

La copie transformée est une opération générique : elle sert aussi au palier 2
(motifs, symétries à venir).

- Chaque sous-entité est copiée **une fois**. Les arêtes, sommets et faces
  partagés le restent dans la copie : une boîte copiée a toujours 12 arêtes et
  8 sommets.
- **Géométrie** :
  - les NURBS sont transformées par leurs pôles cartésiens, avec la fonction
    `gbs::transform` existante ;
  - une `CurveOnSurface` (côté de face naturelle) suit sa surface transformée ;
  - les **pcurves sont partagées** : un déplacement rigide ne change pas le
    paramétrage ;
  - tout autre type de géométrie lève `BRepError`.
- Noms et identifiants externes sont recopiés.

`Placement` est une transformation rigide écrite sans Eigen (matrice 3 × 3 et
translation), avec composition et inverse exactes pour une rotation.

### 3.2 Des représentations aux instances

Le lecteur ne lit plus seulement des formes : il produit des **instances**,
une forme lue **une seule fois** et son placement dans le repère de la
représentation.

| Source | Instances |
|---|---|
| élément BREP de la représentation | la forme, placement identité |
| `mapped_item(source, cible)` | celles de la représentation projetée, composées avec `cible · origine⁻¹` (`representation_map`) |
| représentation liée sans transformation (`shape_representation_relationship`) | les siennes : la représentation BREP d'une pièce derrière sa `shape_representation` |

Les repères sont lus dans le **contexte de leur représentation**, donc dans
ses unités.

### 3.3 Structure produit et aplatissement

1. Les **racines** sont les produits qui ont des composants
   (`next_assembly_usage_occurrence.relating`) sans être eux-mêmes composants.
2. Le placement d'une occurrence vient de sa
   `context_dependent_shape_representation` : une `representation_relationship_with_transformation`
   portant une `item_defined_transformation(item_1, item_2)`.
   - La transformation amène `item_1` (dans `rep_1`) sur `item_2` (dans
     `rep_2`) : placement `F(item_2) · F(item_1)⁻¹`.
   - `rep_1` est d'habitude celle de l'enfant, comme chez OCCT. Si c'est
     `rep_2` qui appartient à l'enfant, le placement est inversé.
3. La forme d'un produit, avec le placement accumulé `P` :
   - ses instances propres, placées par `P` ;
   - et, récursivement, celles de chaque composant, placées par
     `P · placement de l'occurrence`.

   Un produit sans composant et avec une seule forme rend cette forme. Sinon,
   c'est un **compound au nom du produit**, dont l'identifiant externe est le
   `#id` de sa `product_definition`.
4. **Matérialisation** : la première instance non déplacée d'une forme est la
   forme lue elle-même ; toutes les autres sont des `copy_shape`. Une forme lue
   mais seulement utilisée par des copies est effacée à la fin : le modèle ne
   garde que les instances.
5. Les représentations **non atteintes** sont lues comme avant : pièces seules
   et fichiers sans assemblage. Celles qui sont la cible d'un `mapped_item` ou
   d'une relation avec transformation ne sont lues que comme instances.

Une profondeur de plus de 64 niveaux est refusée, comme garde contre les
cycles. Une occurrence dont le placement est illisible est rapportée, et ses
voisines sont lues.

## 4. Écarts par rapport au document de conception

| Document (§ 7) | Implémentation | Raison |
|---|---|---|
| copie transformée par `transform` de gbs sur les NURBS | `copy_shape` générique dans le BREP, pcurves partagées | réutilisable, et la topologie partagée est conservée |
| un compound par sous-assemblage | idem, et une pièce à forme unique sans composant est rendue telle quelle | pas de compound à un seul élément pour une pièce |
| — | `mapped_item` pris en charge | autre écriture des instances, fréquente hors OCCT |
| liaisons Python | en PR 10 avec `read_step` | regroupées avec la lecture |

## 5. Tests

| Test | Couvre |
|---|---|
| `tests_brep_edit.placement` | repère, composition, inverse, identité |
| `tests_brep_edit.copy_shape_moves_a_solid` | boîte cousue déplacée et tournée : volume, partage conservé (12 arêtes, 8 sommets), noms, original intact, `check()` ; face naturelle (courbes sur surface) |
| `tests_step_read.assemblies` | `assembly.stp` : ASM { BOX déplacée, BOX tournée d'un quart de tour, SUB { CYL } }, noms, positions par boîtes englobantes, volume total, aucun original résiduel ; `mapped_items.stp` : deux instances placées par `mapped_item`, dont une tournée |
| `tests_step_occt.assembly` (job OCCT) | pièce partagée instanciée trois fois, rotations autour d'axes inclinés et sous-assemblage, écrite par OCCT en mode assemblage : même nombre de solides, même volume, même boîte englobante qu'OCCT |

Toutes les suites BREP et STEP passent, ainsi que tout le binaire de tests OCCT
(25 cas) et les tests Python du BREP.

## 6. Points à valider

1. **`Placement` et `copy_shape` dans le BREP**, topologie partagée conservée,
   pcurves partagées.
2. **Instances** : forme lue une fois, première instance non déplacée = la forme,
   autres = copies ; originaux inutilisés effacés.
3. **Placement d'une occurrence** `F(item_2) · F(item_1)⁻¹`, inversé si l'enfant
   est `rep_2`.
4. **Compound au nom du produit**, `externalId` = `#id` de la
   `product_definition`.
5. **`mapped_item`** pris en charge.
