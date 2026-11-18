# Volume signé exact sur les faces à trous

| | |
|---|---|
| PR | branche `fix/brep-volume-faces-with-holes` |
| Document de référence | [brep_pr08_solid.md](brep_pr08_solid.md), § 3.1 (`face_volume_integral`, `signed_volume`) |
| Fichiers | `gbs-brep/closure.h`, `gbs/offsets.h` (`surface_flat_knots`), `tests/tests_brep_solid.cpp`, commentaires de `tests/tests_brep_boolean.cpp`, `tests_brep_split.cpp`, `tests_brep_sweep.cpp`, `tests_brep_mesh.cpp`, `brep_pr08_solid.md` (renvoi), `News.md`, cette note |
| Statut | À valider (revue de code : § 8) |

## 1. Objectif

`signed_volume(m, shell, n)` sert à `make_solid` et à `check` pour orienter les
solides, et à `boolean.h` pour distinguer solides et cavités. Le lecteur STEP
l'appelle donc deux fois par solide. Les tests de booléens, de balayage et de
découpe s'en servent aussi pour mesurer les volumes. Il était faux sur les faces
à trous et sur les faces découpées.

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

Il n'y a plus de test d'appartenance. Chaque fil compte avec le sens attendu de
lui : direct pour le contour extérieur, indirect pour un trou, quel que soit son
propre sens. L'aire signée de son polygone (16 points par co-arête) donne ce
sens ; un fil parcouru à l'envers voit sa contribution changer de signe. La
région est ainsi la même que pour l'ancienne grille, qui ignorait les sens, même
sur une face que `check` signale (`OuterWireNotDirect`, `InnerWireNotIndirect`).
Dans un fil, chaque co-arête est prise dans son sens (une co-arête `Reversed`
comptée négativement).

`u₀` est le plus petit u des nœuds des règles le long des fils. N'importe quelle
constante conviendrait ; celle-ci garde F dans la boîte des points où il est
calculé.

### 3.2 Quadratures

Les intégrales utilisent la règle de Gauss-Legendre à 10 points de Boost
(`boost::math::quadrature::gauss<T, 10>`, tabulée, exacte jusqu'au degré 19).
Soit k = max(2, n / 16) et les nœuds de la surface donnés par
`surface_flat_knots` (`gbs/offsets.h`) : ceux d'une B-spline, rationnelle ou non,
ou de la base d'une surface décalée.

- **le long de chaque pcurve** : sur ses intervalles lisses (`smooth_spans` de
  `gbs/curvebreaks.h` : nœuds d'une B-spline, jonctions d'une composite), coupés
  là où la pcurve traverse un nœud de la surface, en u ou en v. Les traversées
  sont repérées sur un polygone de 8 cordes par intervalle, puis placées par
  dichotomie. Chaque morceau est coupé en parts égales, de 1 à k, au prorata de
  sa longueur en (u,v) rapportée à la taille de la boîte. Un point où dv/dt est
  nul (iso-v, arête dégénérée à un pôle) ne donne pas de nœud.
- **le long de u (calcul de F)** : sur k panneaux égaux de la largeur des nœuds
  en u, coupés aux nœuds en u de la surface. Les intégrales de u₀ jusqu'à chaque
  coupure sont **tabulées** en v : aux 16 points de Tchebychev de chaque
  intervalle de v (k intervalles égaux, coupés aux nœuds en v), calculées quand
  un nœud tombe dans l'intervalle pour la première fois. Elles sont interpolées
  en v (formule barycentrique, degré 15). F en un nœud vaut l'interpolée jusqu'à
  la dernière coupure avant u, plus une règle sur le reste du panneau.

Exactitude, à l'arrondi près, quel que soit n et quels que soient les trous :

- **face plane à pcurves droites** : f est polynomiale de bas degré ;
- **B-spline polynomiale de degré 3 au plus, pcurves droites** : dans une cellule
  de nœuds, f est de degré (8, 8), F de degré (9, 8), donc de degré 17 le long
  d'un segment (règle exacte jusqu'à 19). Les intégrales tabulées sont de degré 8
  en v (interpolation exacte jusqu'à 15). Les morceaux étant coupés aux
  traversées de nœuds, chaque règle reste dans une cellule.

Ailleurs, l'erreur décroît comme une puissance élevée de 1/n.

Un paramètre hors des bornes de la surface de 1e-5 de leur étendue au plus est
ramené sur elles : une pcurve à sa tolérance d'un bord en sort de quelques 1e-9,
et une B-spline lève une exception dès 2,2e-14 hors de ses bornes. L'erreur ainsi
introduite est du second ordre en l'écart. Au-delà, la surface est évaluée telle
quelle (une surface fermée prise dans un autre tour).

### 3.3 Coût

Un nœud coûte une règle le long de u (10 évaluations de f, soit 30 de la
surface) et une interpolation à 16 points, quel que soit le nombre de panneaux.
Une table coûte 16 × 10 évaluations de f par panneau, pour chaque intervalle de
v où tombe un nœud. Le coût suit donc le nombre de nœuds le long des fils plus
celui des cellules de nœuds rencontrées, et non leur produit.

Mesures (n = 64, une face ou les fichiers STEP des tests) :

| Cas | grille (avant) | Green sans tables (1re version) | Green avec tables |
|---|---|---|---|
| 24 fichiers STEP des tests | 10,8 ms | 18,0 ms | 13,2 ms |
| B-spline de 100 travées en u, pcurves de 100 travées | 0,28 ms | 284 ms (7,2 millions d'évaluations) | 15 ms (350 000) |
| 100 travées, pcurves de 400 travées | 0,26 ms | 1 153 ms | 29 ms |
| 400 travées, pcurves de 400 travées | 0,29 ms | 4 640 ms | 58 ms |

La grille faisait n² évaluations quelle que soit la face, mais sa précision était
d'ordre 1. La première version recalculait F depuis u₀ en chaque nœud, en
sommant tous les panneaux : son coût suivait le produit des travées des pcurves
par les nœuds en u de la surface. Sur une surface à nombreux nœuds dans les deux
sens, le coût des tables suit le nombre de cellules rencontrées, ce qui est le
minimum pour une intégrale exacte par cellule.

### 3.4 API

Inchangée : `face_volume_integral(m, face, n = 64)` et
`signed_volume(m, shell, n = 64)`. Le paramètre n garde son sens de finesse, mais
sa valeur ne correspond plus à une grille ; en dessous de 48, k vaut 2.

`surface_flat_knots(s)` (`gbs/offsets.h`) devient public. C'était
`offset_detail::knots`, déjà utilisé par `offset_margin`. `closure.h` inclut
désormais `gbs/offsets.h`.

## 4. Résultats

Écart relatif au volume exact, n = 64 :

| Solide | avant | 1re version | après |
|---|---|---|---|
| boîte ouverte (anneau plan) | 1,5e-2 | < 1e-15 | < 1e-15 |
| sphère (une face naturelle) | ~1e-4 | 8e-11 | 7,4e-11 |
| tore (une face naturelle) | ~1e-10 | 1,4e-10 | 1,4e-10 |
| cylindre fermé par deux disques découpés | ~1e-3 | 5,7e-8 | 5,7e-8 |
| face sur un patch bicubique de 3 × 3 travées, pcurves droites | 4,9e-3 | 9,1e-6 | < 1e-14 |

Sur la sphère, le tore et le cylindre, l'écart ne dépend plus de n (le même à 64
et à 160). Il vient de leur représentation et, pour le cylindre, des pcurves
circulaires, approchées à la tolérance de `make_face` ; il ne vient pas de la
quadrature. Sur les solides des fichiers STEP, l'écart au volume exact est
inférieur à 1e-13 (sphère, tore partiel, coude torique, boîte), et de 4,5e-9 pour
le cylindre (pcurves de ses disques).

## 5. Tests (`tests_brep_solid`)

| Test | Couvre |
|---|---|
| `planar_face_with_holes` | face plane avec un trou rectangulaire et un trou triangulaire hors de toute grille : ⅓ z × aire à 1e-12 pour n = 2, 7, 64, 160 ; trou circulaire à 1e-7 ; la même face aux fils parcourus à l'envers (signalée `OuterWireNotDirect`) : même valeur |
| `volumes_of_boxes_with_holes` | boîte évidée ouverte (11 faces) et boîte percée de part en part (10 faces, valide par `check`) : volume à 1e-12 pour n = 7, 64, 160 ; volume du maillage de la boîte ouverte |
| `face_volume_converges` | remplace `face_volume_by_rows` (qui comparait la recherche des cellules ligne par ligne au test cellule par cellule, désormais sans objet) : faces de `plate_inch`, `band_with_hole`, `cylinder` et `sphere.stp`, n = 64 contre n = 512 à 1e-9 relatif |
| `volumes_of_step_solids` | solides de `box`, `sphere`, `cylinder`, `torus_elbow` et `occt_partial_torus.stp` contre leur volume exact : 1e-13 relatif, 1e-8 pour le cylindre (pcurves) |
| `face_volume_on_bspline_surfaces` | losange à pcurves droites sur un patch bicubique ondulé de 3 × 3 et 5 × 5 travées, pcurves d'une ou de 7 travées : même valeur à 1e-13 pour n = 2, 16, 64, 160 et 512 ; pcurves hors des bornes d'un patch bilinéaire de 1e-9 (un côté, un sommet, les quatre côtés) : pas d'exception, valeur exacte à 1e-9 ; surface de 100 travées et pcurves de 100 travées : moins de 600 000 évaluations de la surface (350 000 ; 7,2 millions sans les tables) |
| `curved_solids_volumes` | tolérances resserrées : sphère et tore 1e-8, cylindre à disques découpés 1e-6 (au lieu de 1e-3 et 3e-3) |

Toute la suite passe, y compris les volumes des booléens, des balayages, des
découpes et de la lecture et de l'écriture STEP.

## 6. Limites

- La précision sur une face découpée reste celle de ses pcurves : une pcurve
  approchée à 1e-7 donne un volume à 1e-7 près, pas mieux.
- Les coupures aux nœuds ne connaissent que les B-splines et la base d'une
  surface décalée. Sur une surface de révolution ou une extrusion d'une courbe
  B-spline, les nœuds de la courbe ne coupent ni les panneaux ni les morceaux :
  l'ordre baisse là, sans perte d'exactitude sur les faces planes.
- F est intégré de u₀ à u à travers la boîte des nœuds, y compris hors de la
  face. La surface doit pouvoir y être évaluée : c'est déjà ce que supposait la
  grille.
- Une pcurve qui sort des bornes d'une B-spline de plus de 1e-5 de leur étendue
  fait encore lever une exception (une face hors de sa surface).
- Les tables tiennent 16 valeurs par panneau et par intervalle de v rencontré :
  2,6 millions de réels au pire pour une surface de 400 × 400 travées, libérés en
  sortie.

## 7. Points à valider

1. Remplacement complet de la grille par l'intégrale de contour, sans option de
   repli.
2. Correspondance n → k = max(2, n / 16), la règle à 10 points de Boost et
   l'interpolation de degré 15 en v.
3. Chaque fil pris avec le sens attendu de lui, quel que soit le sien.
4. Paramètres ramenés sur les bornes jusqu'à 1e-5 de leur étendue.
5. `surface_flat_knots` public dans `gbs/offsets.h`, inclus par `closure.h`.
6. Suppression de `face_volume_by_rows`, remplacé par des tests de convergence
   et de volumes exacts.
7. Tolérances resserrées dans `curved_solids_volumes`.
8. Note séparée, avec un renvoi depuis le § 3.1 de `brep_pr08_solid.md`.

## 8. Revue de code

| Point relevé | Suite |
|---|---|
| La surface était évaluée au v exact des pcurves et jusqu'à leur u. Une B-spline lève une exception dès 2,2e-14 hors de ses bornes : une pcurve qui en sortait de 1e-9 faisait échouer le calcul (l'ancienne grille rendait 0,1583333335). Dans le lecteur STEP, la pièce était perdue (`make_solid`) ou toute la lecture échouait (`check`) | corrigé : un paramètre hors des bornes de 1e-5 de leur étendue au plus est ramené sur elles (erreur du second ordre en l'écart) ; au-delà, la surface est évaluée telle quelle (surface fermée dans un autre tour). Testé : un côté, un sommet, les quatre côtés à 1e-9 hors du patch, valeur exacte à 1e-9 |
| Coût : chaque nœud recalculait F depuis u₀ en sommant tous les panneaux. Le coût suivait le produit des travées des pcurves par les nœuds en u : 284 ms au lieu de 0,28 ms pour 100 × 100, 4,6 s pour 400 × 400 | corrigé : intégrales tabulées en v aux points de Tchebychev, interpolées (§ 3.3) : 15 ms et 58 ms ; 13,2 ms sur les fichiers STEP au lieu de 18,0. Testé par le nombre d'évaluations (350 000 au lieu de 7,2 millions) |
| Le volume dépendait du sens des fils en (u,v) : un trou parcouru dans le sens direct était ajouté, un contour horaire changeait de signe. `check` et `make_solid` ne vérifient pas ces sens avant le calcul | corrigé : chaque fil compte avec le sens attendu de lui, d'après l'aire signée de son polygone (§ 3.1). Testé : la face à trous aux fils inversés donne la même valeur |
| Les morceaux le long d'une pcurve n'étaient pas coupés aux traversées des nœuds de la surface : 9e-6 d'erreur à n = 64 sur un patch bicubique de 3 × 3 travées, et pire à n = 160 qu'à n = 64 sur 5 × 5 | corrigé : coupures aux traversées en u et en v, et intervalles des tables coupés aux nœuds en v ; règle à 10 points : exact à l'arrondi sur une B-spline polynomiale de degré 3 au plus à pcurves droites (§ 3.2). Testé à 1e-13 pour n = 2 à 512 |
| `detail::surface_knots_u` réécrivait deux fonctions existantes et ignorait la base d'une surface décalée | corrigé : `offset_detail::knots` devient `surface_flat_knots`, public, et sert aux deux directions |
| `face_volume_converges` ne comparait la méthode qu'à elle-même ; rien sur une B-spline à nœuds traversés, une pcurve hors bornes ni le coût | corrigé : `volumes_of_step_solids` (volumes exacts) et `face_volume_on_bspline_surfaces` |
| La note de la PR 8 décrivait toujours la grille ; cette note se trompait sur le coût et sur les points évalués | corrigé : renvoi dans `brep_pr08_solid.md` § 3.1, §§ 3.1 à 3.3 et 6 réécrits |
| La table de Gauss-Legendre écrite à la main doublait celle de Boost | corrigé : `boost::math::quadrature::gauss<T, 10>` (tabulée par Boost) |
| Boucles par indices sur les coupures et les intervalles | corrigé : `std::views::zip(x, x \| std::views::drop(1))`, l'idiome déjà employé dans gbs (`curve_pieces.h`, `bscanalysis.h`) ; `std::views::pairwise` manque à la libc++ d'Apple, avec laquelle la CI compile sur macOS. `std::ranges::upper_bound` pour la dernière coupure avant u |
