# Écriture STEP — PR préliminaire P2 : provenance des instances

| | |
|---|---|
| PR | branche `feat/brep-instances` |
| Document de référence | [step_writer.md](step_writer.md), § 12.2 (décidé : instances partagées et placements réels à l'écriture, géométrie aplatie dans le modèle) |
| Fichiers | `gbs-brep/placement.h` (nouveau, `Placement` sorti de `edit.h`), `gbs-brep/model.h`, `gbs-brep/edit.h`, `gbs-io/step/reader.h`, `python/gbsBindBrep.cpp`, `tests/tests_brep_model.cpp`, `tests/tests_brep_edit.cpp`, `tests/tests_step_read.cpp`, `python/tests/test_step.py`, `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

La lecture d'un assemblage (PR 9 de la lecture) aplatit chaque occurrence en
une copie placée de sa pièce, puis oublie que ces copies sont la même pièce.
Pour écrire l'assemblage avec ses **pièces partagées** et leurs **placements**,
le modèle doit garder cette provenance, sans changer la géométrie aplatie.

## 2. API

```cpp
template <std::floating_point T>
struct InstanceOf { std::uint64_t part; Placement<T> place; };   // placement : pièce -> modèle

std::uint64_t Model::newPart();                                  // clé unique dans le modèle
void Model::setInstance(const ShapeId &, std::optional<InstanceOf<T>>);
std::optional<InstanceOf<T>> Model::instance(const ShapeId &) const;
```

Python : `brep.new_part(m)`, `brep.set_instance(m, shape, brep.InstanceOf(part, place))`,
`brep.instance(m, shape)` (`None` si absent).

## 3. Choix d'architecture et justification

### 3.1 Un attribut, pas un placement dans la topologie

La provenance est une **table creuse** à côté des noms et des identifiants
externes (PR 2 de la lecture). La géométrie reste dupliquée par instance :
l'exploration, `check`, les volumes, les boîtes, la couture et l'export IGES
ne changent pas. C'est le choix du § 12.2 de la conception ; les vrais
placements dans la topologie restent un chantier de palier 2.

### 3.2 Des clés de pièce données par le modèle

Les clés viennent de `Model::newPart()`, un compteur du modèle. `setInstance`
refuse une clé qui n'en vient pas. On aurait pu prendre le `#id` STEP de la
`product_definition`, mais deux fichiers lus dans le même modèle auraient alors
pu partager des clés par hasard.

### 3.3 Cycle de vie

| Opération | Effet |
|---|---|
| `erase` | supprime la provenance de l'entité |
| `compact` | réindexe |
| `append` | recopie, en **décalant les clés** de l'autre modèle après les siennes |
| `copy_shape` | une instance copiée reste une instance de la même pièce, son placement composé avec celui de la copie (`p ∘ place`) ; vaut pour toutes les entités copiées, donc pour les instances imbriquées dans un compound |

### 3.4 Le lecteur remplit la provenance

La forme de chaque produit lu par `product_shape` porte sa provenance :

- une **clé de pièce par `product_definition`**, lecture par lecture ;
- le placement est le placement composé de l'occurrence, du repère de la
  pièce au repère du modèle ;
- une pièce, un sous-assemblage et la racine (placement identité) ont chacun
  leur provenance. L'écrivain peut ainsi partager aussi un sous-assemblage
  utilisé deux fois.

Les instances placées par un `mapped_item` n'ont pas de provenance : elles
seront écrites comme pièces distinctes. Leur partage pourra être ajouté si un
cas réel le demande.

### 3.5 `Placement` dans son propre en-tête

`Placement` passe de `edit.h` à `gbs-brep/placement.h`, que `model.h` inclut.
L'API ne change pas.

## 4. Écarts par rapport au document de conception

| Document (§ 12.2) | Implémentation | Raison |
|---|---|---|
| clé = `#id` de la `product_definition` | clé donnée par `Model::newPart()` | unicité dans le modèle, plusieurs lectures (§ 3.2) |
| provenance des pièces | aussi des sous-assemblages et de la racine | partage d'un sous-assemblage à l'écriture |
| — | `mapped_item` sans provenance | cas plus rare, gardé simple (§ 3.4) |

## 5. Tests

| Test | Couvre |
|---|---|
| `tests_brep_model.instances` | clés, refus d'une clé inconnue, retrait, `erase`, `compact`, `append` (clés décalées, compteur repris), identifiant mort |
| `tests_brep_edit.copy_shape_composes_instances` | copie d'un compound : l'instance intérieure garde sa pièce, placement composé ; le compound n'en devient pas une |
| `tests_step_read.assembly_provenance` | `assembly.stp` : les deux BOX partagent leur pièce, sommets ramenés à la pièce identiques ; placements (translation, quart de tour) ; cylindre placé à travers SUB ; racine identité ; seconde lecture avec ses propres clés |
| `test_step.py::test_assembly_provenance` | liaisons Python |

Toutes les suites BREP et STEP et tout le binaire OCCT passent.

## 6. Points à valider

1. **Provenance en attribut**, géométrie aplatie inchangée.
2. **Clés données par le modèle** (`newPart`), décalées par `append`.
3. **`copy_shape` compose** le placement d'une instance copiée.
4. **Provenance aussi pour les sous-assemblages et la racine.**
5. **`mapped_item` sans provenance** pour l'instant.

## 7. Revue de code (a posteriori, lot 1 de la revue du noyau BREP)

| Commentaire | Suite donnée |
|---|---|
| composition du placement d'une instance copiée (`p * place`, `place` appliqué d'abord) | vérifiée : correcte |
| clés de pièces lors d'un `append` du modèle dans lui-même | décalées comme pour un autre modèle ; testé (`append_to_itself`) |
| `Placement::frame` ne normalise pas et `inverse` suppose R orthonormée | le lecteur STEP construit des repères orthonormés (`Frame`) : à vérifier avec le lecteur (lot 8) |
