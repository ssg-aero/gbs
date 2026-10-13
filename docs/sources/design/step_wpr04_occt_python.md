# Écriture STEP — PR 4 : validation par OpenCascade, Python, documentation

| | |
|---|---|
| PR | branche `feat/step-write-occt-python` |
| Document de référence | [step_writer.md](step_writer.md), § 8 (Python) et § 9 (points 3 et 4) |
| Fichiers | `gbs-occt/tests/tests_step_occt.cpp`, `python/gbsBindBrep.cpp`, `python/tests/test_step.py`, `ReadMe.md`, `docs/sources/design/step_writer.md` (statut, bilan), `News.md` ; en commit séparé, les révisions de conception de la PR 138 qui n'avaient pas atteint `master` |
| Statut | À valider |

## 1. Objectif de la PR

Valider l'écrivain par un lecteur de référence, le rendre utilisable en
Python, et clore le plan de l'écriture.

## 2. Validation par OpenCascade (job OCCT de la CI)

### 2.1 Double aller-retour

Les 20 formes de la comparaison de la lecture (PR 8) font désormais le parcours
complet : **OCCT → STEP → gbs → STEP → OCCT**.

- **Formes** : primitives, booléens, congés, révolutions, extrusion, loft, tore
  incliné, sphère percée, tore partiel.
- **Critères** : le fichier écrit par gbs est relu par `STEPControl_Reader`,
  et la forme obtenue doit être **valide** pour `BRepCheck_Analyzer`. Elle
  doit avoir les **mêmes nombres de faces, d'arêtes non dégénérées et de
  sommets** que la forme d'origine, et le **même volume** (2e-3, 5e-3 pour
  les congés et le loft).
- **Aucune approximation** : tout est primitive ou B-spline exacte.

### 2.2 Assemblage

L'assemblage OCCT de la PR 9 de la lecture est réécrit par gbs. Il contient une
pièce partagée instanciée trois fois, des rotations autour d'axes inclinés et
un sous-assemblage. Seuls **deux solides** sont écrits, la boîte et le pion,
une fois chacun. Relu par OCCT, on retrouve le même nombre de solides
(quatre), le même volume et la même boîte englobante.

### 2.3 Formes construites dans gbs

Ce sont des formes qui ne viennent d'aucun fichier :
- une boîte cousue de faces B-spline bilinéaires ;
- une sphère et un tore en faces naturelles sur primitives, avec leurs seams
  et leurs pôles ;
- un tube ouvert.

Elles sont écrites par gbs et lues par OCCT : formes valides, nombres de faces
et volumes exacts.

### 2.4 Limite de la validation

`STEPControl_Reader` applique par défaut son traitement de forme
(ShapeFix), qui peut corriger de petits défauts. La validité est donc mesurée
**après** ce traitement, comme pour tout utilisateur d'OCCT. Les nombres
d'entités et les volumes identiques montrent qu'il n'a rien eu à reconstruire.

## 3. API Python

```python
report = brep.write_step(model, shapes, path, schema="AP214", approx_tol=1e-5, pcurves=True,
                         flatten=False, uncertainty=None, product_name="", author="", organization="")
text, report = brep.write_step_text(model, shapes, ...)
```

- `shapes` est une forme ou une liste de formes.
- Le rapport `StepWriteReport` donne les nombres de pièces, d'assemblages,
  d'occurrences et d'entités, les approximations et les formes omises.
- Une erreur lève `gbs.brep.BRepError` : forme morte, schéma inconnu, fichier
  non inscriptible.
- La convention est symétrique de `read_step` / `read_step_text` : un chemin
  par défaut, le texte sur demande.

## 4. Documentation

- Le `ReadMe.md` montre l'écriture en Python et résume ce que lisent et
  écrivent le lecteur et l'écrivain.
- `step_writer.md` passe au statut **réalisé**, avec un § 13 de bilan : les PR
  et leurs notes, les décisions révisées et les limites connues.
- **Rattrapage** : les révisions de la conception décidées le 12 octobre (§ 12
  sur les primitives et les instances, notes de révision dans
  `step_reader.md` et `brep_core.md`) avaient été poussées après la fusion de
  la PR 138 et n'étaient pas sur `master`. Elles sont reprises ici dans un
  **commit séparé**, sans modification. Les notes P1 et P2 y renvoient.

## 5. Tests

| Test | Couvre |
|---|---|
| `tests_step_occt` (7 cas, 425 assertions) | double aller-retour des 20 formes, assemblage réécrit, formes construites dans gbs |
| `test_step.py::test_write_and_read_back` | écriture d'un fichier, relecture, volume ; écriture en texte AP242 avec nom de produit |
| `test_step.py::test_write_assembly_and_errors` | assemblage (2 assemblages, 2 solides, 4 occurrences), `flatten`, schéma refusé, dossier inexistant |

Toutes les suites BREP, STEP et IGES, tous les tests Python du BREP et de STEP
et tout le binaire OCCT passent.

## 6. Points à valider

1. **Critères de la validation OCCT** : validité après le traitement d'OCCT,
   nombres d'entités et volumes identiques, aucune approximation.
2. **API Python** symétrique de la lecture.
3. **Rattrapage des révisions de conception** dans un commit séparé.
