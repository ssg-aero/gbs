# Palier 2 — PR 9 b : booléens avec faces confondues

| | |
|---|---|
| PR | branche `feat/brep-boolean-coincident` |
| Document de référence | [brep_palier2.md](brep_palier2.md), § 6.3 ; [brep2_pr09_boolean.md](brep2_pr09_boolean.md), § 1 |
| Fichiers | `gbs-brep/coincident.h` (nouveau), `gbs-brep/boolean.h`, `gbs-brep/classify.h`, `gbs-brep/section.h`, `gbs-brep/explore.h` (`bounding_box`), `gbs-brep/brep`, `tests/tests_brep_boolean.cpp`, `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Traiter dans les booléens les **faces confondues** : faces planes
coplanaires, et primitives coaxiales de même rayon. C'est le complément
annoncé par la PR 9. Exemples : des boîtes posées ou alignées sur un même
plan, des cylindres coaxiaux bout à bout ou qui se chevauchent, une forme
fusionnée avec elle-même.

## 2. Empreinte des faces confondues (`coincident.h`)

Pour chaque paire de faces confondues (signalée par la section, PR 6) :

1. **Croisements des bords.**
   - Les arêtes du bord des deux faces sont croisées en 3D (courbe-courbe,
     PR 1).
   - Une arête est coupée là où une arête de l'autre face la croise, y
     commence ou y finit. Le sommet vient du même réservoir que celui de la
     section : un croisement déjà sommet de section le reste.
   - Deux sommets confondus sont **fusionnés**.
   - On répète jusqu'à ce qu'il n'y ait plus rien à couper.
2. **Arêtes superposées.**
   - Après l'étape 1, deux arêtes superposées ont les mêmes sommets et le
     même milieu : elles sont **fusionnées en une seule**.
   - Les co-arêtes de l'arête absorbée prennent l'arête gardée, avec une
     pcurve recalculée sur leur face (exacte si linéaire).
   - Cette pcurve est **décalée d'une période** pour rester du côté du seam
     de l'ancienne. Sans cela, les deux côtés du seam d'un cylindre se
     retrouvaient au même endroit, et le bord de la face était cassé.
3. **Coupes.** Les arêtes du bord de chaque face situées **dans** l'autre face
   en deviennent des coupes, avec leur pcurve sur cette face.

Seules les primitives sont traitées : c'est la forme fermée de la section qui
détecte la coïncidence, et les pcurves s'y calculent en formule fermée.

## 3. Sélection des morceaux confondus (`boolean.h`)

Un morceau classé **sur** l'autre opérande est apparié à la face de l'autre
qui porte son point intérieur. On compare les normales sortantes (normale de
la surface, retournée selon l'usage de la face) :

| Opération | Même côté | Côtés opposés (contact) |
|---|---|---|
| fusion | le morceau de A, une fois | aucun : la face de contact disparaît |
| intersection | le morceau de A, une fois | aucun |
| différence A − B | aucun : B enlève la matière au-dessus | le morceau de A |

Les morceaux de B confondus ne sont jamais gardés : celui de A les représente.

Les arêtes d'un opérande posées sur une face de l'autre (`boundary_edges` de
la section) sont acceptées si elles bordent une face confondue. Ailleurs,
elles rendent toujours `NonGeneric`.

## 4. Corrections trouvées en route

| Défaut | Effet | Correction |
|---|---|---|
| `bounding_box` d'une face bornée par ses seuls wires = boîte de ses arêtes (palier 1) | une sphère bornée par son seam avait la boîte d'un **demi-cercle** : la section ne voyait pas la paire sphère / face de boîte, et un contact tangent passait inaperçu | la surface est échantillonnée sur la boîte des pcurves de **chaque** face |
| boîtes échantillonnées sous-estimées (flèche entre les points) | un point d'un cylindre en x = −5 était hors d'une boîte en x ≥ −4,76, et déclaré dehors | les filtres de `section` et de `classify` élargissent leur boîte de 5 % de sa taille |
| rayon partant d'une surface | la racine en t ≈ 0 peut tomber à −1e-16, et l'impact est manqué | `classify(solide)` teste d'abord si le point est **sur** une face (projection, puis classement dans la face) |

La section du tore partiel (PR 6) a désormais 3 arêtes au lieu de 4 : le test
place son plan au milieu de la boîte englobante du tore, et cette boîte a
changé. Le wire reste fermé.

## 5. Tests (`tests_brep_boolean`, 4 cas, 88 assertions)

| Test | Couvre |
|---|---|
| `coincident_faces` | boîtes décalées en x et y (dessous et dessus coplanaires) : fusion 40, intersection 8, différence 16, opérandes intacts ; **empilées** face contre face : fusion 48 en un solide de 10 faces, intersection vide, différence 24 ; en escalier : 48 ; **la même boîte deux fois** : fusion 24 en 6 faces, intersection 24, différence vide ; **cylindres coaxiaux** R 5 en z 0..10 et 5..15 : fusion 15 × 25π, intersection et différence 5 × 25π |
| `non_generic` | sphère tangente au dessus d'une boîte : `NonGeneric`, opérandes intacts ; erreurs |
| `boxes`, `cylinders_and_spheres` | inchangés (PR 9) |

Les 27 suites BREP, STEP, IGES, d'intersection et des primitives passent.

## 6. Points à valider

1. **Empreinte croisée** des faces confondues, avec fusion des sommets et des
   arêtes superposées.
2. **Décalage d'une période** des pcurves recalculées sur un seam.
3. **Règles de sélection** des morceaux confondus (tableau du § 3).
4. **`bounding_box`** échantillonnant la surface de toutes les faces : la boîte
   est plus grande et prudente.
5. **Test « sur une face »** avant le lancer de rayon.
