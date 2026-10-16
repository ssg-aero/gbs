# Palier 2 bis — PR 7 : robustesse près des tangences

| | |
|---|---|
| PR | branche `fix/brep-tangent-tolerances` |
| Document de référence | [brep2b_pr06_exact_restriction.md](brep2b_pr06_exact_restriction.md), § 2.2 et point à valider 2 |
| Fichiers | `gbs-brep/split.h`, `gbs/intersect/curve_pieces.h`, `tests/tests_brep_split.cpp`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

La PR 6 avait noté une fragilité : sur le cas OCCT `tangent_slivers` (une face
de boîte tangente à un cylindre), un écart de 1,6e-13 sur l'extrémité d'une
courbe de section suffisait à rendre un shell ouvert. Cette PR en trouve la
cause et la corrige.

## 2. Reproduction

Les paramètres de coupe de la restriction exacte ont été perturbés
volontairement de ±1e-13 à ±1e-7, en valeur relative, le temps de l'essai (le
crochet n'est pas gardé).

| Perturbation | Avant | Après |
|---|---|---|
| +1e-13, +3e-13 | `tangent_slivers` en fusion : shell ouvert | passe |
| −1e-13 … ±1e-7 (comparaison OCCT), ±1e-11 … ±1e-8 (`tests_brep_boolean`) | — | passe |

## 3. Cause

Sur le flanc du cylindre, les deux copeaux sont deux lentilles qui se touchent
le long de la génératrice de tangence. La région extérieure a pour trou le
contour de leur union, qui **longe les arêtes des deux lentilles**.

`split_face` donne chaque boucle-trou à la plus petite boucle extérieure qui
contient son premier point. Ce point est un sommet d'une lentille, donc sur son
contour : selon l'arrondi, le test de parité le comptait dedans, et le trou
était donné à la lentille. La région extérieure perdait son trou, la lentille
devenait une face aberrante classée dedans, et le flanc du cylindre manquait au
résultat.

C'est le même défaut que la boucle jumelle corrigée à la PR 4, dans une forme
plus générale : un trou qui partage des arêtes avec une autre boucle, sans en
avoir exactement les mêmes.

## 4. Correction

Le point testé pour un trou est décalé, depuis le milieu d'une de ses
demi-arêtes, **du côté gauche** : là où se trouve la face qui le possède,
hors de toute boucle qui partage ses arêtes. Le décalage vaut 1e-3 de la marge de
nœud de l'arrangement (1e-6 de la diagonale de la boîte de la face). La règle de
la boucle jumelle est gardée.

## 5. Performance : décomposition de Bézier

En mesurant les tests, `coincident_closed_bsplines` (PR 5) prenait 72 s, dont
63 s pour l'intersection d'un cylindre B-spline rationnel avec sa copie. Les
pcurves du disque plan bordé par le cercle rationnel sont projetées et
interpolées, avec environ 1000 nœuds, car le paramétrage rationnel n'est que C¹
à ses nœuds.

`CurvePieces` les décomposait en Bézier en insérant les nœuds un par un
(O(n²) en mémoire et en temps), et cela pour chaque courbe de section. La
décomposition suit désormais l'algorithme A5.6 du NURBS Book : tous les nœuds
intérieurs sont élevés au degré d'un coup, en O(n p²).

| | Avant | Après |
|---|---|---|
| `coincident_closed_bsplines` | 72 s | 8,8 s |
| `ctest` complet (-j8) | — | 49 s |

## 6. Tests

| Test | Cas |
|---|---|
| `hole_against_pinched_regions` (découpe) | deux losanges qui se touchent par un sommet dans une face carrée, tournés de kπ/8 : trois faces, la face extérieure garde son trou de 8 arêtes. **Test de couverture** : il passe aussi avec l'ancien code, dont le premier point ne tombait pas sur un sommet compté dedans. La régression elle-même est montrée par la perturbation du § 2 |

Les suites C++ (471 tests), le binaire OCCT (29 cas) et les 51 tests Python
passent.

## 7. Points à valider

1. Le décalage du point de test (1e-6 de la diagonale de la boîte de la face).
2. **Suite possible** : des pcurves exactes sur un plan (projection affine des
   pôles, poids gardés) éviteraient les pcurves de 1000 nœuds. Elles seraient
   rationnelles, ce que l'empreinte (`replace_edge`) et l'écriture STEP
   doivent alors accepter.
