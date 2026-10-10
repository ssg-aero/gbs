# Lecture STEP — PR 2 : attributs du modèle (noms, identifiants externes, unité)

| | |
|---|---|
| PR | branche `feat/brep-attributes` |
| Document de référence | [step_reader.md](step_reader.md), § 6.1 et question 10 ; [brep_core.md](brep_core.md), § 9.2 (tables réservées au palier 1) |
| Fichiers | `gbs-brep/model.h` (attributs), `gbs-io/iges_brep.h` (label IGES = nom de la face), `python/gbsBindBrep.cpp`, `tests/tests_brep_model.cpp`, `tests/tests_brep_iges.cpp`, `python/tests/test_brep.py`, `docs/sources/design/step_reader.md` (correction `from_chars`, commit séparé), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Donner au `Model` les trois attributs que le palier 1 avait réservés sans les
écrire, et dont le lecteur STEP a besoin : le **nom** d'une entité (produit,
face nommée), son **identifiant externe** (le `#id` de l'instance STEP d'origine)
et l'**unité** du modèle.

## 2. API

```cpp
m.setName(shape, "bracket");          // "" retire le nom
std::string_view n = m.name(shape);   // "" si aucun ; vue valide jusqu'au prochain setName de cette entité
m.setExternalId(shape, 1234);         // std::nullopt retire l'identifiant
std::optional<std::int64_t> id = m.externalId(shape);
m.setUnitScale(25.4);                 // longueur d'une unité du modèle en millimètres (1 par défaut)
double s = m.unitScale();
```

Python : `m.set_name`, `m.name`, `m.set_external_id` (`None` retire),
`m.external_id`, propriété `m.unit_scale`.

## 3. Choix d'architecture et justification

### 3.1 Tables creuses à côté de l'arène

Une table de hachage par type d'entité et par attribut (`std::array` de sept
`std::unordered_map<index, valeur>`), à côté des tables d'entités. Un modèle
natif, sans nom ni identifiant, ne paie rien : pas de champ dans les agrégats
`Vertex`, `Edge`…, dont la taille et la sémantique restent celles du palier 1.
Les attributs n'apparaissent donc pas dans les copies d'entités (`m.edge(id)`) :
ils appartiennent au modèle, pas à l'entité.

### 3.2 Cycle de vie aligné sur les entités

| Opération | Effet sur les attributs |
|---|---|
| `erase(id)` | supprime les attributs de l'entité (sinon `compact()` les ferait migrer vers un survivant) |
| `compact()` | réindexe les clés avec l'`IdRemap` ; ceux des entités mortes disparaissent |
| `append(other)` | recopie ceux de `other` avec les identifiants remappés |
| accès sur un identifiant mort ou invalide | `BRepError`, comme les accesseurs d'entités |

### 3.3 L'unité décrit, elle ne convertit pas

`unitScale()` dit ce que valent les coordonnées (1 unité = `unitScale` mm) ;
`setUnitScale` ne touche **pas** à la géométrie. La conversion est le travail du
lecteur STEP (§ 6.2 du document : conversion vers l'unité cible au moment de la
construction). Conséquence pour `append` : deux modèles d'unités différentes ne
sont pas fusionnés silencieusement, `append` lève `BRepError` ; c'est à
l'appelant de convertir l'un des deux.

### 3.4 Noms dans l'export IGES

`add_brep` (palier 1, PR 9) donne à une face son **propre nom** comme label IGES
quand elle en a un, et sinon le nom passé en argument suivi de son rang, comme
avant. Un nom lu dans un fichier STEP se retrouvera donc dans l'IGES exporté.
Le label d'une entrée de répertoire IGES fait 8 caractères, aligné à droite.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| `unit_scale` membre public | `unitScale()` / `setUnitScale()` avec validation (fini, > 0) | refuser une unité nulle ou négative |
| — | `append` refuse des unités différentes | éviter un mélange silencieux de millimètres et de pouces |
| Noms utilisés par l'IGES (palier 1, § 9.2) | fait ici | — |

## 5. Tests

| Test | Couvre |
|---|---|
| `tests_brep_model.attributes` (24 assertions) | noms et identifiants sur solide, face, arête, sommet ; retrait ; identifiants invalides refusés ; unité et validation ; `erase` supprime, `compact` réindexe sans faire migrer le nom d'une face effacée ; `append` recopie et remappe ; unités différentes refusées |
| `tests_brep_iges.box_solid` (+2 assertions) | une face nommée « LID » a ce label dans le fichier, les cinq autres « box2 » à « box6 » |
| `test_brep.py::test_attributes` | API Python, `None` retire, propriété `unit_scale`, erreurs |

Toutes les suites BREP passent, ainsi que les tests Python du BREP.

## 6. Points à valider

1. **Attributs en tables creuses** dans le modèle, pas dans les entités.
2. **`erase` supprime les attributs**, `compact` et `append` les remappent.
3. **`unitScale` déclaratif** (pas de conversion) et **`append` refusé** entre unités différentes.
4. **Label IGES = nom propre de la face** quand il existe.
