# Palier 2 bis — PR 7 : robustesse près des tangences

| | |
|---|---|
| PR | branche `fix/brep-tangent-tolerances` |
| Document de référence | [brep2b_pr06_exact_restriction.md](brep2b_pr06_exact_restriction.md), § 2.2 et point à valider 2 |
| Fichiers | `gbs-brep/split.h`, `gbs/bscurve.h`, `gbs/intersect/curve_pieces.h`, `gbs/intersect/surface_pieces.h`, `tests/tests_brep_split.cpp`, `tests/tests_bscurve.cpp`, `News.md` |
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

Le point testé pour un trou est pris **à côté du trou** : depuis le milieu du
premier segment de son polygone, sur la normale gauche (le côté de la face qui
le possède), **à mi-chemin du premier bord rencontré** dans cette direction,
quelle que soit la boucle. Il est ainsi dans la face voisine du trou, quelle que
soit son épaisseur : pas de traversée d'une région plus mince qu'un décalage
fixe. La règle de la boucle jumelle est gardée.

(La première version de cette PR prenait un décalage fixe de 1e-6 de la
diagonale ; la revue a montré qu'un anneau plus mince le traversait.)

## 5. Performance : décomposition de Bézier

En mesurant les tests, `coincident_closed_bsplines` (PR 5) prenait 72 s, dont
63 s pour l'intersection d'un cylindre B-spline rationnel avec sa copie. Les
pcurves du disque plan bordé par le cercle rationnel sont projetées et
interpolées, avec environ 1000 nœuds, car le paramétrage rationnel n'est que C¹
à ses nœuds.

`CurvePieces` les décomposait en Bézier en insérant les nœuds un par un
(O(n²) en mémoire et en temps), et cela pour chaque courbe de section. La
décomposition suit désormais l'algorithme A5.6 du NURBS Book, en O(n p²) :
`gbs::bezier_poles` (`bscurve.h`, à côté des outils de nœuds), qui rend le
réseau de pôles raffiné dans un seul tampon (n_seg p + 1 pôles), ou rien quand
les nœuds ne sont pas serrés aux bouts ou qu'un nœud intérieur dépasse le degré
(une discontinuité : retour à l'échantillonnage). `SurfacePieces` l'applique
ligne par ligne puis colonne par colonne (A5.7) au lieu d'insérer les nœuds un
par un.

| | Avant | Après |
|---|---|---|
| `coincident_closed_bsplines` | 72 s | 7,6 s |
| `ctest` complet (-j8) | — | 50 s |

## 6. Tests

| Test | Cas | Ancien code | Première version de la PR |
|---|---|---|---|
| `hole_against_pinched_regions` | deux losanges qui se touchent par un sommet dans une face carrée, tournés de kπ/8 ; la face trouée doit avoir le carré comme contour extérieur | **échoue** (le trou va à un losange) | passe |
| `hole_from_a_leftmost_vertex` | le trou commence au sommet le plus à gauche d'un losange, que le test de parité compte dedans | **échoue** | passe |
| `hole_of_a_thin_ring` | île dans un anneau de 1e-7, plus mince que le décalage fixe | passe | **échoue** |
| `bezier_poles` (`tests_bscurve`) | degrés 1 à 5, nœuds de multiplicités 1 à p, rationnelles ou non : chaque segment de Bézier égal à la courbe en 21 points à 1e-12, et pôles égaux à ceux de l'insertion nœud par nœud ; multiplicité p + 1 et nœuds non serrés refusés | — | — |

Les suites C++ (474 tests), le binaire OCCT (29 cas) et les 51 tests Python
passent.

### 6.1 Revue de code

| Commentaire | Suite donnée |
|---|---|
| A5.6 : un nœud intérieur de multiplicité p + 1 sous-déborde (`p − mult`) et donne des pôles nuls | `bezier_poles` refuse les multiplicités au-delà du degré et les nœuds non serrés : retour à l'échantillonnage ; testé |
| décalage fixe de 1e-6 de la diagonale, de l'ordre d'une lamelle à la tolérance | point à mi-chemin du premier bord rencontré sur la normale gauche ; testé par l'anneau mince |
| différence finie au lieu de la dérivée | la direction est maintenant celle du segment du polygone, sur lequel porte le test de parité |
| A5.6 à mettre dans gbs, `SurfacePieces` insère toujours un nœud à la fois | `gbs::bezier_poles` dans `bscurve.h`, utilisé par `CurvePieces` et, ligne par ligne puis colonne par colonne, par `SurfacePieces` |
| `Q` en vecteur de vecteurs | tampon plat de n_seg p + 1 pôles et deux segments de travail |
| aucun test n'épingle les corrections | les trois tests de découpe et `bezier_poles` ci-dessus |

## 7. Points à valider

1. Le point de test du trou à mi-chemin du premier bord.
2. **Suite possible** : des pcurves exactes sur un plan (projection affine des
   pôles, poids gardés) éviteraient les pcurves de 1000 nœuds. Elles seraient
   rationnelles, ce que l'empreinte (`replace_edge`) et l'écriture STEP
   doivent alors accepter.
