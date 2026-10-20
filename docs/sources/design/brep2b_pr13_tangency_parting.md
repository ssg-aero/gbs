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

## 2. Où finit la tangence

La tangence finit là où les surfaces cessent de se courber pareil : la petite
valeur propre de II1 − II2 (`Direction::small`, nulle le long d'une courbe de
contact) dépasse `angular_tol` fois la grande, le même certificat qu'au départ
(PR 11). Ce point ne dépend pas de la tolérance : x = angular_tol / (6 R) pour la
face −max(x, 0)³ sous un cylindre de rayon R (1,7e-4 pour R = 1, 3,3e-5 pour
R = 5), à deux pas minimaux près (20 tol). Avant la revue, la branche finissait
où l'écart atteignait la tolérance, x = tol^(1/3) : 0,046 à 1e-4, 0,01 à 1e-6,
0,002 à 1e-8.

## 3. Ce qui suit la fin

À travers la courbe, l'écart g = (S1 − S2) · n2 croît comme b s² / 2, où b est la
grande valeur propre, **signée**, de II1 − II2 (`Direction::bend`).

`TangentMarcher::parting` regarde au-delà de la fin dans la direction de la
courbe, à 10 tol puis en doublant jusqu'à un pas : le même système (pied et
normales parallèles) y donne l'**extremum de l'écart à travers**, à n'importe
quelle distance (`Correction::any_distance`), jusqu'à ce qu'il dépasse la
moitié de la tolérance.

| Verdict | Condition | Fin de la branche |
|---|---|---|
| `parts` | l'écart a le signe de b : il ne s'annule pas à travers | `parted`, pas `singular` |
| `crosses` | l'écart a le signe opposé : il s'annule deux fois à travers, deux courbes transverses partent de la fin | `singular` |
| `boundary` | la sonde sort des boîtes avant que l'écart soit mesurable : les surfaces restent dans la tolérance jusqu'au bord | la branche se pose sur le bord (`boundary`) |
| `unknown` | pas de solution, ou écart trop petit jusqu'à un pas | `singular` (prudent) |

`Marcher::End` gagne `parted`, que seul `TangentMarcher` rend. Les options de
`TangentMarcher::correct` sont nommées (`Correction{.fixed, .value,
.any_distance}`).

## 4. Effet

| Cas | Avant | Après |
|---|---|---|
| bloc sous z = −max(x, 0)³, cylindre couché dessus (la face s'écarte au-delà de x = 0) | `NonGeneric` | fusion 7,5 + 250π (deux solides), commun vide, différence 7,5 ; une arête de section tangente, de (−1, 0, 0) à x ≈ 3,3e-5 |
| bloc sous z = +max(x, 0)³ (la face monte dans le cylindre) | `NonGeneric` | `NonGeneric`, modèle inchangé |

## 5. Tests

| Test | Cas |
|---|---|
| `tangency_limits` (`tests_intersect_marching`) | face z = −b max(x, 0)³ sous un cylindre de rayon 1, dans **les deux ordres** des surfaces et pour les tolérances 1e-4, 1e-6, 1e-8 : une courbe tangente de x = −1 à x = 1e-3 / 6 (à 20 tol près), `singular` seulement pour b = −1 ; `parting` à la fin : `parts` pour b = 1, `crosses` pour b = −1. Face s'écartant à 0,005 de son bord x = 1 : la branche va jusqu'au bord, pas `singular` ; `parting` depuis un point passé la fin : `boundary`, avec l'itéré hors des boîtes |
| `tangency_ending_in_a_face` (`tests_brep_boolean`) | témoin à dessus plat ; b = 1 : section ni `singular` ni incomplète, une arête tangente de (−1, 0, 0) à 0 < x < 0,05, fusion, commun et différence acceptés avec leurs volumes ; b = −1 (volume du bloc 8,5) : section `singular` et complète, les trois booléens `NonGeneric`, chaque type de forme inchangé, opérandes valides et volume inchangé |

Mutations vérifiées :

| Mutation | Attrapée par |
|---|---|
| `parts` toujours faux | les deux tests (b = 1 redevient `singular`, booléens refusés) |
| `parts` toujours vrai | le test direct de `parting` seulement : dans le cas b = −1, les courbes transverses rendent déjà l'intersection `singular` (elles partent d'un point où les normales sont parallèles). Le critère est une seconde garde pour ce cas |
| signe inversé | les deux assertions du test d'intersection |
| pas de pose au bord après `boundary` | **non attrapée** par l'intersection complète : un pas y passe le bord et se pose par le chemin ordinaire. Le verdict `boundary` est épinglé par le test direct |

Ces mutations ont été faites à la main, le temps de l'essai ; aucune n'est un test.

Les 493 tests C++ (binaire OCCT compris) et les 54 tests Python passent.

### 5.1 Revue de code

| Commentaire | Suite donnée |
|---|---|
| la fin dépend de la tolérance (x = tol^(1/3)) ; la fin géométrique est où `small` quitte zéro | un point n'est accepté que si les surfaces se courbent pareil ; la fin est à x = angular_tol / (6 R) quelle que soit la tolérance, testé pour 1e-4, 1e-6, 1e-8 |
| sonde à 10 tol hors des boîtes : une séparation à moins de 10 tol du bord rendue `singular` | sonde de 10 tol jusqu'à un pas ; sortie des boîtes avant un écart mesurable : verdict `boundary`, la branche se pose sur le bord |
| `-1, T(0), true` en paramètres positionnels | `Correction{.fixed, .value, .any_distance}`, initialisateurs désignés aux appels |
| convention de signe jamais testée avec les surfaces échangées ; refus de `parts` non testés ; mutations faites à la main | les deux ordres testés, directement aussi ; verdict `boundary` testé ; `unknown` (pas de solution) non construit. Les mutations restent faites à la main (tableau ci-dessus) |

## 6. Points à valider

1. Le critère du signe de l'extremum de l'écart, un pas de 10 tol au-delà.
2. La fin géométrique, où les surfaces cessent de se courber pareil (seuil
   `angular_tol` sur le rapport des valeurs propres), au lieu de l'écart à la
   tolérance.
