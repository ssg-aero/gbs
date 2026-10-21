# Palier 4 — PR 0 : exemples visuels et documentation

| | |
|---|---|
| PR | branche `feat/brep-examples` |
| Document de référence | [brep_palier4.md](brep_palier4.md), § 0 bis (document du palier, livré avec cette PR) |
| Fichiers | `docs/sources/examples/` (nouveau : `gbs_examples.py`, `gallery.py`, `booleans.ipynb`, `meshing.ipynb`, `index.rst`), `docs/index.rst`, `python/tests/test_examples.py` (nouveau), `.github/workflows/ci.yml`, `docs/sources/design/brep_palier4.md` (nouveau), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Voir les booléens et le maillage, au lieu de les juger seulement par des
volumes, des nombres de faces et OpenCascade. C'est aussi la base des exemples
de chaque PR du palier 4 et des rendus de la présentation de GBS.

## 2. Contenu

### 2.1 `gbs_examples.py`

- **Primitives** construites avec l'API du BREP elle-même, sans fichier :
  - `box` : prisme d'une face plane ;
  - `cylinder` : prisme d'un disque ;
  - `cone` : révolution d'un triangle ;
  - `sphere` : révolution d'un demi-disque ;
  - `torus` : révolution d'un disque.

  Elles sont placées par une `Placement` (origine et axe). Avec `exact=True`, les cercles sont des `EllipseCurve` et les balayages donnent des cylindres, sphères et tores exacts. Avec `exact=False`, ce sont leurs NURBS, et les surfaces sont NURBS.
- **Cas de booléens** (`CASES`), repris de la comparaison OCCT des paliers 2 et 2 bis :
  - boîtes ; boîte et cylindre ; boîte et sphère ;
  - cylindres croisés, inégaux et égaux (points de branchement) ;
  - boîtes coplanaires et empilées ; cylindres coaxiaux ;
  - sphère qui touche une boîte ; cylindre couché sur une boîte ;
  - cône coupé par son sommet ; tore et boîte ; sphère centrée sur un cylindre ;
  - une limite connue : la sphère NURBS qui touche une boîte.

  `run(case)` calcule `fuse`, `cut` et `common`, et garde les erreurs.
- **Rendu PyVista** :
  - une couleur par face, ce qui montre les découpes ;
  - les arêtes du BREP en noir ;
  - les opérandes en transparence ;
  - `plot_case` met côte à côte les opérandes et les trois résultats, avec pour chacun le nombre de faces, le volume du maillage et sa fermeture.

### 2.2 Galerie et notebooks

- `gallery.py --out dossier` :
  - une image PNG par cas, plus six maillages vus avec leurs triangles ;
  - `gallery.md`, qui les rassemble avec, pour chaque résultat, le nombre de faces, le volume, la validité (`check`) et la fermeture.
  - Durée : 4,5 s en local.
- **CI** : la galerie est rendue sur Linux (xvfb) après les tests Python et publiée en artefact `gallery` de chaque build.
- `booleans.ipynb` et `meshing.ipynb` sont interactifs, en PyVista avec trame, et en images statiques à défaut.
  - Les deux notebooks sont exécutés sans erreur en local, puis livrés sans leurs sorties.
  - La CI ne les exécute pas (question 11 du palier).
- **Sphinx** : une page « Examples » dans la table des matières.

### 2.3 Test

`python/tests/test_examples.py` rejoue chaque cas :
- chaque résultat est valide, son maillage fermé et sans face repliée ;
- les volumes des maillages vérifient fuse + common = a + b et cut + common = a, à 5e-3 près ;
- la limite connue échoue encore sur les trois opérations : sa levée sera vue.

## 3. Ce que les exemples ont montré

| Constat | Suite |
|---|---|
| Les 13 cas sur primitives exactes donnent des résultats valides et fermés, conformes à la vue | — |
| Une **sphère NURBS** qui touche une boîte en un point donne `NonGeneric` sur les trois opérations, alors que la sphère exacte passe : le contact ponctuel n'est reconnu qu'entre primitives | limite connue, gardée comme cas de la galerie et du test ; à traiter avec les contacts tangents (palier 2 bis, PR 1) |
| La révolution d'un segment perpendiculaire à l'axe (disque du fond d'un cône) donne une surface NURBS, pas un plan : les faces d'une révolution sont naturelles (tout le rectangle des paramètres), et sur un plan ce rectangle serait un carré | voulu ; noté |

## 4. Résultats de la galerie

| Cas | fuse | cut | common |
|---|---|---|---|
| boîtes | 12 faces, V = 42 | 9 faces, V = 18 | 6 faces, V = 6 |
| boîte et sphère au coin | 7 faces, V = 36,30 | 7 faces, V = 22,24 | 4 faces, V = 1,764 |
| cylindres croisés égaux | 8 faces, V = 32,18 | 7 faces, V = 13,46 | 5 faces, V = 5,315 |
| tore et boîte | 8 faces, V = 126,0 | 9 faces, V = 67,24 | 7 faces, V = 12,78 |
| sphère NURBS qui touche une boîte | `NonGeneric` | `NonGeneric` | `NonGeneric` |

Tous les résultats sont valides et fermés, sauf la limite connue. En local passent :
- les 538 tests C++ (sans changement) ;
- les 73 tests Python, dont les 15 nouveaux.

## 5. Points à valider

1. **Emplacement** `docs/sources/examples`, à côté des notebooks de théorie.
2. **Galerie en artefact de CI** plutôt que versionnée dans le dépôt : pas d'images dans l'historique.
3. **Notebooks livrés sans sorties**, non exécutés par la CI.
4. **Limite connue gardée** dans les cas, pour qu'on la voie et que sa levée soit remarquée.
