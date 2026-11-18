# Volume signé exact sur les faces à trous

| | |
|---|---|
| PR | branche `fix/brep-volume-faces-with-holes` |
| Document de référence | [brep_pr08_solid.md](brep_pr08_solid.md), § 2 (`face_volume_integral`, `signed_volume`) |
| Fichiers | `gbs-brep/closure.h`, `tests/tests_brep_solid.cpp`, commentaires de `tests/tests_brep_boolean.cpp`, `tests_brep_split.cpp`, `tests_brep_sweep.cpp`, `tests_brep_mesh.cpp`, `News.md`, cette note |
| Statut | À valider |

## 1. Objectif

`signed_volume(m, shell, n)` sert à `make_solid` et à `check` pour orienter les
solides, et à `boolean.h` pour distinguer solides et cavités. Les tests de
booléens, de balayage et de découpe s'en servent aussi pour mesurer les volumes.
Il était faux sur les faces à trous et sur les faces découpées.

Cas de référence : une boîte 2 × 1,5 × 1, évidée à 0,1 et ouverte sur le dessus
(`cut` par la boîte 0,1..1,9 × 0,1..1,4 × 0,1..1). Le résultat a 11 faces planes,
dont l'anneau du dessus, percé d'un trou. Volume exact : 0,894 (maillage : 0,894).

| n | 64 | 160 | 256 |
|---|---|---|---|
| avant | 0,88103 | 0,89775 | 0,89489 |
| après | 0,894 | 0,894 | 0,894 |

## 2. Diagnostic

`face_volume_integral` calculait ⅓ ∬ S · (Su × Sv) du dv par la règle du point
milieu sur une grille n × n de la boîte (u,v) du contour extérieur. Une cellule
comptait quand son centre était dans la région (contour moins trous).

Le trou était bien soustrait. L'erreur venait de l'appartenance par centre : tout
bord qui ne tombe pas sur une ligne de la grille coûte de l'ordre d'une cellule
d'aire sur sa longueur. Sur l'anneau, paramétré sur [0, 1]², le trou va de 0,05 à
0,95 en u et de 0,0667 à 0,9333 en v. Il ne tombe sur la grille pour aucun n
usuel. La face donnait 0,22375 au lieu de 0,22 à n = 160. L'erreur est en O(1/n),
sans signe constant : elle varie de façon erratique avec n (tableau ci-dessus).

## 3. Contenu

### 3.1 Intégrale sur les contours (théorème de Green)

Soit f(u, v) = S · (Su × Sv) et F(u, v) = ∫ de u₀ à u de f(s, v) ds. Alors
∂F/∂u = f et, par Green, sur une région dont le bord est parcouru dans le sens
direct (contour extérieur trigonométrique, trous horaires) :

```
∬_région f du dv = ∮_bord F(u, v) dv
```

Les trous sont pris par leur propre sens de parcours : il n'y a plus de test
d'appartenance. Les fils sont parcourus comme dans la face (le sens de chaque
co-arête appliqué, une co-arête `Reversed` comptée négativement).

`u₀` est le plus petit u du contour extérieur (polygone de 16 points par
co-arête). N'importe quelle constante conviendrait. Ce choix garde les points
évalués dans la boîte (u,v) des pcurves, c'est-à-dire là où l'ancienne grille
évaluait déjà la surface.

### 3.2 Quadratures

Les deux intégrales utilisent une règle de Gauss-Legendre à 8 points (exacte
jusqu'au degré 15), tabulée dans `detail::gauss_legendre_8`. Soit
k = max(2, n / 16) :

- **le long de chaque pcurve** : sur ses intervalles lisses (`smooth_spans` de
  `gbs/curvebreaks.h` : nœuds d'une B-spline, jonctions d'une composite).
  Chaque intervalle est coupé en morceaux égaux, de 2 à k, au prorata de sa
  longueur en (u,v) rapportée à la taille de la boîte. Un point où dv/dt est nul
  (iso-v, arête dégénérée à un pôle) ne contribue pas, et F n'y est pas calculé ;
- **le long de u (calcul de F)** : sur k panneaux égaux de la largeur de la boîte
  en u, coupés aux nœuds en u si la surface est une B-spline (rationnelle ou non),
  puis le panneau partiel jusqu'à u.

Sur une face plane à pcurves droites, f est polynomiale de bas degré en (u, v)
et l'intégrande le long d'un segment l'est aussi : le résultat est exact à
l'arrondi près, quel que soit n. Ailleurs, l'erreur décroît comme une puissance
élevée de 1/n.

### 3.3 Coût

L'ancien calcul faisait n² évaluations par face (n = 64 : 4 096 ; n = 160 :
25 600 ; n = 256 : 65 536), quel que soit le contour. Le nouveau en fait environ
8 × morceaux par co-arête, fois 4 × k en moyenne pour F, et rien sur les iso-v.
Pour une face de quatre côtés : environ 2 000 à n = 64 et 13 000 à n = 160. Le
coût croît avec le nombre de co-arêtes, mais chaque co-arête courte n'a que deux
morceaux. La suite complète (771 tests) passe en 13,8 s.

### 3.4 API

Inchangée : `face_volume_integral(m, face, n = 64)` et
`signed_volume(m, shell, n = 64)`. Le paramètre n garde son sens de finesse, mais
sa valeur ne correspond plus à une grille.

## 4. Résultats

Écart relatif au volume exact, n = 64 :

| Solide | avant | après |
|---|---|---|
| boîte ouverte (anneau plan) | 1,5e-2 | < 1e-15 |
| sphère (une face naturelle) | ~1e-4 | 8e-11 |
| tore (une face naturelle) | ~1e-10 | 1,4e-10 |
| cylindre fermé par deux disques découpés | ~1e-3 | 5,7e-8 |

Sur le cylindre, l'écart ne dépend plus de n (5,7e-8 à 64 comme à 160). Il vient
des pcurves circulaires, approchées à la tolérance de `make_face`, et non de la
quadrature.

## 5. Tests (`tests_brep_solid`)

| Test | Couvre |
|---|---|
| `planar_face_with_holes` | face plane avec un trou rectangulaire et un trou triangulaire hors de toute grille : ⅓ z × aire à 1e-12 pour n = 2, 7, 64, 160 ; trou circulaire : à 1e-7 |
| `volumes_of_boxes_with_holes` | boîte évidée ouverte (11 faces) et boîte percée de part en part (10 faces, valide par `check`) : volume à 1e-12 pour n = 7, 64, 160 ; volume du maillage de la boîte ouverte |
| `face_volume_converges` | remplace `face_volume_by_rows` (qui comparait la recherche des cellules ligne par ligne au test cellule par cellule, désormais sans objet) : faces de `plate_inch`, `band_with_hole`, `cylinder` et `sphere.stp`, n = 64 contre n = 512 à 1e-9 relatif |
| `curved_solids_volumes` | tolérances resserrées : sphère et tore 1e-8, cylindre à disques découpés 1e-6 (au lieu de 1e-3 et 3e-3) |

Toute la suite passe (771 tests), y compris les volumes des booléens, des
balayages, des découpes et de la lecture et de l'écriture STEP.

## 6. Limites

- La précision sur une face découpée reste celle de ses pcurves : une pcurve
  approchée à 1e-7 donne un volume à 1e-7 près, pas mieux.
- Les panneaux en u sont coupés aux nœuds de la surface, mais les morceaux le
  long d'une pcurve ne le sont pas là où elle traverse un nœud en v. L'ordre
  baisse alors, sans perte d'exactitude sur les faces planes (une B-spline plane
  de degré 1 n'a pas de nœud intérieur utile).
- F est intégré de u₀ à u à travers la boîte, y compris hors de la face. La
  surface doit pouvoir y être évaluée : c'est déjà ce que supposait la grille.

## 7. Points à valider

1. Remplacement complet de la grille par l'intégrale de contour, sans option de
   repli.
2. Correspondance n → k = max(2, n / 16) et la règle à 8 points.
3. Suppression de `face_volume_by_rows`, remplacé par un test de convergence.
4. Tolérances resserrées dans `curved_solids_volumes`.
5. Note séparée plutôt qu'une section de revue dans `brep_pr08_solid.md`.
