# Palier 2 bis — PR 13 : tangence à bout libre

| | |
|---|---|
| PR | branche `feat/tangency-parting` |
| Document de référence | [brep2b_pr12_tangency_ending_in_a_face.md](brep2b_pr12_tangency_ending_in_a_face.md), point à valider 2 |
| Fichiers | `gbs/intersect/marching.h`, `tests/tests_intersect_marching.cpp`, `tests/tests_brep_boolean.cpp`, `docs/sources/design/brep_palier2bis.md`, `docs/sources/design/brep2b_pr12_tangency_ending_in_a_face.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Une tangence qui cesse à l'intérieur d'une face terminait la branche en
`singular`, et le booléen était refusé (`NonGeneric`, PR 12). La revue de la
PR 12 a mesuré que, sans ce drapeau, la fusion du cas de test donne déjà le bon
volume : la section et la découpe acceptent une arête tangente à bout libre, qui
ne découpe pas la face (palier 2 bis, point 2 des décisions : arêtes pendantes
acceptées si tangentes).

Mais une tangence peut cesser de deux façons :

| Au-delà de la fin | Intersection | Traitement |
|---|---|---|
| les surfaces **s'écartent** | la courbe de tangence s'arrête là | bout libre accepté |
| les surfaces **se croisent** | deux courbes transverses partent du bout de la tangence, où les normales sont parallèles | reste `singular` |

## 2. Le critère

À travers la courbe, l'écart g = (S1 − S2) · n2 croît comme b s² / 2, où b est la
grande valeur propre, **signée**, de II1 − II2 (la même que pour la direction,
PR 11 : `Direction::bend`).

Quand la tangence se perd (le pas tombe sous le pas minimal, la correction ne
trouvant plus d'écart dans la tolérance), `TangentMarcher::parts` regarde un peu
plus loin, à 10 tol du dernier point, dans la direction de la courbe : le même
système (pied et normales parallèles) y a encore une solution, l'**extremum de
l'écart à travers**, cette fois à n'importe quelle distance (`any_distance`).

- Même signe que b : l'écart ne s'annule pas à travers, les surfaces
  s'écartent : fin `parted`, pas `singular`.
- Signe opposé : l'écart s'annule deux fois à travers, les surfaces se croisent
  : `singular`.
- Écart sous la moitié de la tolérance, ou pas de solution : `singular` (prudent).

`Marcher::End` gagne `parted`, que seul `TangentMarcher` rend.

## 3. Effet

| Cas | Avant | Après |
|---|---|---|
| bloc sous z = −max(x, 0)³, cylindre couché dessus (la face s'écarte au-delà de x = 0) | `NonGeneric` | fusion 7,5 + 250π (deux solides), commun vide, différence 7,5 ; une arête de section tangente, de (−1, 0, 0) à x ≈ 0,01 |
| bloc sous z = +max(x, 0)³ (la face monte dans le cylindre) | `NonGeneric` | `NonGeneric`, modèle inchangé |

## 4. Tests

| Test | Cas |
|---|---|
| `tangency_limits` (`tests_intersect_marching`) | face z = −b max(x, 0)³ sous un cylindre de rayon 1 : b = 1, une courbe tangente de x = −1 à la fin du contact, pas `singular` ; b = −1, `singular` ; `parts` au bout de la tangence vrai pour b = 1, faux pour b = −1 |
| `tangency_ending_in_a_face` (`tests_brep_boolean`) | témoin à dessus plat ; b = 1 : section ni `singular` ni incomplète, une arête tangente de (−1, 0, 0) à 0 < x < 0,05, fusion, commun et différence acceptés avec leurs volumes ; b = −1 (volume du bloc 8,5) : section `singular` et complète, les trois booléens `NonGeneric`, chaque type de forme inchangé, opérandes valides et volume inchangé |

Mutations vérifiées :

| Mutation | Attrapée par |
|---|---|
| `parts` toujours faux | les deux tests (b = 1 redevient `singular`, booléens refusés) |
| `parts` toujours vrai | le test direct de `parts` seulement : dans le cas b = −1, les courbes transverses rendent déjà l'intersection `singular` (elles partent d'un point où les normales sont parallèles). Le critère est une seconde garde pour ce cas |
| signe inversé | les deux assertions du test d'intersection |

Les 493 tests C++ (binaire OCCT compris) et les 54 tests Python passent.

## 5. Points à valider

1. Le critère du signe de l'extremum de l'écart, un pas de 10 tol au-delà.
2. Le bout libre accepté dépend de la tolérance : la tangence cesse où l'écart
   atteint la tolérance (x ≈ 0,01 pour x³ et 1e-6).
