# Palier 2 bis — PR 12 : refus d'une tangence qui cesse dans une face

| | |
|---|---|
| PR | branche `test/tangency-ending-in-face` |
| Document de référence | [brep2b_pr11_tangent_curves.md](brep2b_pr11_tangent_curves.md), § 8 et point à valider 2 |
| Fichiers | `tests/tests_brep_boolean.cpp`, `docs/sources/design/brep2b_pr11_tangent_curves.md`, `docs/sources/design/brep_palier2bis.md` |
| Statut | À valider |

## 1. Objectif de la PR

La PR 11 suit les tangences le long d'une courbe. Une tangence qui cesse à
l'intérieur des faces termine la branche en `singular`, et le booléen doit alors
être refusé (`NonGeneric`) et laisser le modèle tel qu'il était. Seul le
`singular` de l'intersection était testé ; la revue l'a relevé, et c'était la
suite consignée au § 8 de la note de la PR 11.

Pas de changement de code : le comportement attendu est déjà celui de master.

## 2. Le cas

Un solide x ∈ [−1, 1], y ∈ [−1, 1], de z = −2 à la face B-spline
z = −b max(x, 0)³ (b = 1 pour le cas, 0 pour le témoin) :

| Face | Surface |
|---|---|
| dessus | B-spline cubique en x, nœud simple en 0, pôles aux abscisses de Greville (−1, −2/3, 0, 2/3, 1), cotes (0, 0, 0, 0, −b) : exactement −b max(x, 0)³ ; linéaire en y |
| côtés y = ±1 | B-splines réglées entre la droite z = −2 et la courbe du dessus (mêmes nœuds, mêmes pôles) |
| dessous, côtés x = ±1 | quadrilatères plans |

Les six faces naturelles sont cousues (`sew`) puis fermées en solide ; volume
8 − 2b ∫₀¹ x³ dx, soit 7,5.

Le cylindre de rayon 5 de `cylinder.stp`, couché le long de x à z = 5, sur
x ∈ [−5, 5], touche le dessus le long de y = 0. Ses seuls contacts avec le bloc
sont cette ligne et son extrémité (−1, 0, 0) sur l'arête du côté x = −1 : il
passe à 5 − √24 ≈ 0,10 au-dessus des côtés y = ±1, et à 1 au-dessus du côté
x = 1. La tangence va donc **du bord de la face** (−1, 0, 0) **à l'intérieur** :
la face s'écarte au-delà de x = 0 (écart x³), et la tangence cesse vers
x ≈ 0,01, où l'écart atteint la tolérance (1e-6) ; ce point dépend donc de la
tolérance.

## 3. Test

| Test | Vérifie |
|---|---|
| `tangency_ending_in_a_face` (`tests_brep_boolean`), témoin | dessus plat (b = 0), tangent tout du long : fusion 8 + 250π (deux solides), commun vide, différence 8. Le refus du cas vient donc de la fin de la tangence, pas de la construction ni du contact en (−1, 0, 0) |
| section | `singular` et complète (pas un plafond) ; une seule arête, tangente, de (−1, 0, 0) à 0 < x < 0,05 |
| booléens | fusion, commun et différence rendent `NonGeneric`, avec le message « singular or incomplete intersection » (sept refus différents rendent `NonGeneric`) ; après chaque refus, le nombre de formes de **chaque type** (sommets, arêtes, wires, faces, coques, solides, composés) est inchangé ; les opérandes restent valides, le volume du bloc inchangé |

La vérification de tous les types (`counts`) remplace aussi celle des seules
arêtes et sommets dans `bspline_tangency` (PR 11).

Les 493 tests C++ passent.

### 3.1 Revue de code

| Commentaire | Suite donnée |
|---|---|
| remise en état vérifiée sur les arêtes et sommets seulement : une fuite de faces, coques ou solides passait (mutations d'`erase_shape`) | `counts` compare les sept types, ici et dans `bspline_tangency` |
| rien ne montre que le refus vient de la fin de la tangence | témoin à dessus plat qui fusionne ; arête de section unique, de (−1, 0, 0) à l'intérieur ; message du refus |
| § 8 de la note de la PR 11 contradictoire ; point 2 et § 12.5 du palier sans la PR 12 | corrigés ; numérotation du § 7 rétablie |
| note : extrémités de la branche, jeu du cylindre, dépendance à la tolérance, portée de la remise en état | § 2 et § 3 |
| lisibilité : rangées de pôles des côtés implicites, enveloppes `shared_ptr`, `CHECK(a == b)`, second modèle sans commentaire, volume du bloc après refus | rangées explicites (`side`), `bspline`, `CHECK_EQ`, commentaire, volume vérifié. `std::views::zip` non retenu : pas encore utilisé dans le dépôt, à vérifier avec la libc++ d'Apple |

## 4. Points à valider

1. Le cas construit, avec une face exacte plutôt qu'une interpolation.
2. **Suite possible** : accepter une telle tangence comme une arête à bout libre
   dans la face. Sans le drapeau `singular`, la fusion donne déjà 7,5 + 250π
   (mesuré par la revue) : le reste de la chaîne semble le supporter.
