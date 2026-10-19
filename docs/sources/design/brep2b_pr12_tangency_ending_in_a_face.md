# Palier 2 bis — PR 12 : refus d'une tangence qui cesse dans une face

| | |
|---|---|
| PR | branche `test/tangency-ending-in-face` |
| Document de référence | [brep2b_pr11_tangent_curves.md](brep2b_pr11_tangent_curves.md), § 8 et point à valider 2 |
| Fichiers | `tests/tests_brep_boolean.cpp`, `docs/sources/design/brep2b_pr11_tangent_curves.md` |
| Statut | À valider |

## 1. Objectif de la PR

La PR 11 suit les tangences le long d'une courbe. Une tangence qui cesse à
l'intérieur des faces termine la branche en `singular`, et le booléen doit alors
rester `NonGeneric`, sans rien laisser dans le modèle. Seul le `singular` de
l'intersection était testé ; la revue l'a relevé, et c'était la suite consignée
au § 8 de la note de la PR 11.

Pas de changement de code : le comportement attendu est déjà celui de master.

## 2. Le cas

Un solide x ∈ [−1, 1], y ∈ [−1, 1], de z = −2 à la face B-spline
z = −max(x, 0)³ :

| Face | Surface |
|---|---|
| dessus | B-spline cubique en x, nœud simple en 0, pôles aux abscisses de Greville (−1, −2/3, 0, 2/3, 1), cotes (0, 0, 0, 0, −1) : exactement −max(x, 0)³ ; linéaire en y |
| côtés y = ±1 | B-splines réglées entre la droite z = −2 et la courbe du dessus (mêmes nœuds, mêmes pôles) |
| dessous, côtés x = ±1 | quadrilatères plans |

Les six faces naturelles sont cousues (`sew`) puis fermées en solide ; volume
8 − 2 ∫₀¹ x³ dx = 7,5.

Le cylindre de rayon 5 de `cylinder.stp`, couché le long de x à z = 5, touche le
dessus le long de y = 0 : tangent pour x ≤ 0, la face s'en écarte au-delà (écart
x³), la tangence cesse vers x ≈ 0,01.

## 3. Test

| Test | Vérifie |
|---|---|
| `tangency_ending_in_a_face` (`tests_brep_boolean`) | la section est `singular` et complète (le refus vient bien de la tangence, pas d'un plafond) ; fusion, commun et différence rendent `NonGeneric` ; arêtes et sommets du modèle inchangés après chaque refus ; les opérandes restent valides |

Les 493 tests C++ passent.

## 4. Points à valider

1. Le cas construit, avec une face exacte plutôt qu'une interpolation.
