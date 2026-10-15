# Palier 2 bis — PR 3 : points de branchement

| | |
|---|---|
| PR | branche `feat/brep2b-branch-points` |
| Document de référence | [brep_palier2bis.md](brep_palier2bis.md), § 5 (révisé par cette PR) |
| Fichiers | `gbs/intersect/marching.h`, `gbs/intersect/intersection_curve.h`, `gbs-brep/section.h`, `python/gbsBindBrepOps.cpp`, `tests/tests_intersect_marching.cpp`, `tests/tests_brep_boolean.cpp`, `python/tests/test_intersect.py`, `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/design/brep_palier2bis.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Le cheminement ne s'arrête plus sur un point de branchement, où deux branches
de l'intersection se croisent et où les surfaces sont tangentes. La section y
découpe les courbes, ce qui permet les booléens de cylindres égaux qui se
croisent (Steinmetz) ou d'une sphère et d'un cylindre de Viviani.

## 2. Changements

| Fichier | Changement |
|---|---|
| `marching.h` | `step_over` : quand le pas devient trop court ou que les normales deviennent parallèles, on cherche le plus court pas le long de la tangente qui retombe sur la même branche au-delà (le sinus repasse au-dessus du seuil, avec un creux entre les deux points). Si le point est sur le bord d'une boîte, la branche s'y pose |
| `marching.h` | un pas ordinaire qui enjambe le point se reconnaît au **retournement de n1 × n2** par rapport au sens de parcours (signe suivi le long de la branche). Ce pas est raccourci jusqu'à 64 tol, puis gardé |
| `marching.h` | départs traités du mieux au moins bien conditionné ; prédiction factorisée (`predict`) |
| `intersection_curve.h` | `IntersectionCurve::branch_spans` : intervalles de paramètre franchis, non densifiés (aucun point ne s'y calcule) |
| `section.h` | les intervalles franchis sont croisés deux à deux, sur une même courbe ou entre courbes (courbe-courbe à 100 tol, deux interpolations), et les courbes sont découpées aux croisements |
| `gbsBindBrepOps.cpp` | `branch_spans` en Python |

### 2.1 Révision du document

Le document prévoyait de calculer, au point singulier, les directions des
branches par la seconde forme fondamentale, puis de cheminer sur chaque
demi-branche. C'est inutile ici : une branche traverse le point tout droit, et
la subdivision donne des départs sur toutes les branches. Le § 5 du document
est réécrit en conséquence.

## 3. Tests

| Test | Cas | Attendu |
|---|---|---|
| `crossed_cylinders` (marching) | Steinmetz, seam du second cylindre passant par les points de branchement ou tourné | non singulier ; longueur de deux ellipses de demi-axes √2 et 1 à 1e-5 ; 4 intervalles franchis (seam tourné) |
| `singular_points_and_limits` | Viviani | non singulier ; longueur de la courbe à 1e-5 ; 2 intervalles franchis (le point double, deux fois) ; une sphère dans un cylindre de son rayon (tangence le long d'un cercle, par cheminement) reste `singular` |
| `branch_points` (booléens) | cylindres R 5 de longueur 30, axes z et x | commun : 16 R³ / 3 ; différence : **deux moitiés** qui se touchent aux points de branchement ; fusion |
| OCCT `booleans` | `steinmetz`, `viviani` × 3 opérations | valides, volumes à 2e-3, **mêmes nombres de solides et de faces** que `BRepAlgoAPI` |
| Python | Steinmetz : non singulier, 4 intervalles | |

Les suites `brep`, `intersect`, `step`, `iges`, le binaire OCCT (29 cas) et
les tests Python passent.

## 4. Limites

- **Tangence le long d'une courbe** hors formules fermées : toujours refusée
  (question 3 du document).
- **Viviani d'OCCT** : la sphère d'OCCT a ses pôles sur la courbe quand son
  axe est celui du cylindre. Le test tourne l'axe ; le passage par un pôle est
  l'objet de la PR 4.

## 5. Points à valider

1. La révision du § 5 (franchissement tout droit au lieu des directions de
   départ).
2. La tolérance de croisement des intervalles (100 tol).
