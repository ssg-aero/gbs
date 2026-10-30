# Palier 2 — PR 5 : classification

| | |
|---|---|
| PR | branche `feat/brep-classify` |
| Document de référence | [brep_palier2.md](brep_palier2.md), § 3 ; question 8 |
| Fichiers | `gbs-brep/classify.h` (nouveau), `gbs-brep/brep` (agrégateur), `tests/tests_brep_classify.cpp` (nouveau), `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Dire si un point est dans une face, sur son bord ou dehors, et de même pour un
solide. C'est la troisième étape des booléens (PR 9) : chaque face découpée de
A est classée par rapport à B.

## 2. API

```cpp
enum class Location { In, Out, On };

auto classify(m, FaceId f, point<T, 2> uv, T tol = brep_default_tolerance) -> Location;           // dans l'espace (u, v) de la face
auto classify(m, SolidId s, point<T, 3> p, T tol = …) -> Location;                                  // cavités comprises
auto interior_point(m, FaceId f, T tol = …) -> optional<point<T, 2>>;                               // loin du bord
auto classify(m, FaceId f, SolidId s, T tol = …) -> optional<Location>;                             // face entière, par un point intérieur
```

## 3. Choix d'architecture et justification

### 3.1 Point dans une face

1. **Bord.** Si le point de la surface est à moins de `tol`, ou de la tolérance
   de l'arête, d'une courbe d'arête (projection, PR 1) ou d'un sommet
   dégénéré, la réponse est `On`. Ce test se fait en **3D**, pour que la
   tolérance ait le même sens partout, pôles et seams compris.
2. **Parité.** Une demi-droite partant de (u, v) est croisée avec toutes les
   pcurves des wires (courbe-courbe, PR 1). Un nombre impair de croisements
   veut dire `In`. Les deux côtés d'un seam et les pcurves des arêtes
   dégénérées (pôles) font partie du contour en (u, v), et comptent
   normalement.
3. **Rayon ambigu.** Un rayon qui passe par un sommet (compté par deux
   co-arêtes), tangent à une pcurve, ou le long d'une pcurve, est
   **abandonné**. On en prend un autre, dans une suite fixe de 8 angles
   irréguliers. Le résultat est donc **déterministe**.

Le test des polygones échantillonnés du palier 1 (`inside` dans `closure.h`)
reste utilisé par `check`. Ce nouveau test est exact aux pcurves près.

### 3.2 Point dans un solide

1. **Hors de la boîte englobante** : `Out`, sans autre calcul.
2. **Rayon.** Une demi-droite (`LineCurve`) est croisée avec la surface de
   chaque face, sur la boîte de ses pcurves (courbe-surface, PR 2). Pour les
   primitives, c'est une formule fermée.
3. **Impacts.** Chaque impact est classé dans sa face (§ 3.1) :
   - `Out` : il est ignoré ;
   - impact à moins de `tol` du point : le point est `On` ;
   - impact sur une arête, tangent à la face, ou rayon posé le long d'une
     face : le rayon est **ambigu**.
4. **Directions.** Elles sont prises dans une suite fixe de 16, réparties
   sur la sphère (spirale dorée) et jamais parallèles aux axes.
5. **Cavités.** Leurs faces font partie des faces du solide, et la parité les
   traite sans cas particulier.

### 3.3 Point intérieur d'une face

Sur sept lignes v = constante, réparties dans la boîte des pcurves, on prend
les croisements avec les pcurves. Le milieu de la plus longue corde
intérieure, vérifié par `classify`, est le point choisi. Il est loin du bord,
donc sûr pour classer une face entière.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| `classify(m, f, uv, tol)`, `classify(m, s, p, tol)`, `interior_point` | idem, plus `classify(m, f, s)` | la classification d'une face dont les booléens ont besoin |
| — | test du bord en 3D, avant la parité | tolérance indépendante du paramétrage |

## 5. Tests (`tests_brep_classify`, 4 cas, 141 assertions, 0,2 s)

Les solides des tests de lecture STEP sont réutilisés.

| Test | Couvre |
|---|---|
| `points_in_faces` | plaque en pouces percée : dedans, dans le trou, près du bord du trou, sur le trou, sur le bord, au sommet, dehors, à 25 µm du bord ; points alignés sur les sommets et le trou (rayons ambigus) ; point intérieur |
| `faces_with_seams_and_poles` | face latérale du cylindre (seam utilisé deux fois) : dedans de part et d'autre du seam, sur le seam, sur le cercle, hors de la plage ; sphère bornée par son seam : jusqu'aux pôles, sur le pôle ; **point intérieur de chaque face** de sept fichiers (boîte, cylindre, sphère, cône, bande trouée, tore partiel et loft OCCT) |
| `points_in_solids` | boîte 2 × 3 × 4 : dedans, à 1 µm des faces, dehors, sur une face, une arête, un sommet, grille de points alignés ; cylindre (dedans, coin de la boîte englobante, seam, cercle) ; sphère (centre, près du pôle, pôle) ; cône (axe, dehors, sommet) |
| `cavities_and_faces` | boîte cousue avec une cavité : dedans, dans la cavité (`Out`), sur la cavité ; les 12 faces classées par rapport à une boîte qui entoure la cavité : 6 dedans, 6 dehors |

Les suites `tests_brep_sew`, `tests_brep_check` et `tests_step_write` passent
avec l'agrégateur modifié.

## 6. Points à valider

1. **Bord testé en 3D**, avant la parité.
2. **Rayons ambigus abandonnés**, directions dans une suite fixe
   (déterministe).
3. **Classification d'une face par un point intérieur**, milieu de la plus
   longue corde intérieure.
4. **`classify.h` dans l'agrégateur** `gbs-brep/brep` : il entraîne
   `gbs/intersect/`.

## 7. Revue de code (a posteriori, lot 4 de la revue du noyau BREP)

Chaque commentaire a été reproduit par une sonde ; chaque test cité échoue sur
la version précédente. Le commentaire 1 vaut aussi pour la section
([brep2_pr06_section.md](brep2_pr06_section.md) § 8), qui intersecte les
surfaces sur la même boîte.

| Commentaire | Suite donnée |
|---|---|
| 1. La boîte (u, v) d'une face (`face_uv_box`) était celle de 17 points par pcurve. Elle coupe le bombé d'un bord courbe, et elle sert de borne stricte : un point hors d'elle est `Out`, et les rayons et les projections du test dans un solide s'y limitent. Un disque de rayon 1 dont le cercle commence à 10° a une boîte x ≤ cos 10° = 0,985 : le point (0,99 ; 0) était `Out`. Dans le prisme de ce disque, le point (0,988 ; 0 ; 0,99) aussi, car le premier rayon traverse le dessus là | corrigé : la boîte enveloppe les pcurves. C'est celle des morceaux de Bézier de chaque pcurve (`CurvePieces` : les pôles d'un morceau l'enveloppent), coupés en deux tant qu'ils dépassent de plus de 1e-3 de la taille de la pcurve la boîte de ses points. Une polyligne (degré 1) donne directement ses extrémités et ses pôles. Une courbe qui n'est pas une B-spline garde la boîte estimée de `CurvePieces`. Testé (`near_the_extreme_of_a_curved_boundary` : disque et prisme) |
| Coût : la classification d'un point dans un solide recalculait la boîte (u, v) de chaque face à chaque rayon et à chaque impact. Sur une boîte découpée en 880 faces, un point coûtait 2,3 ms, dont 2,1 ms de boîtes | les boîtes (u, v) sont calculées une fois par point : 2,0 ms, avec la boîte exacte. Restent la boîte du solide et celles des faces en 3D (1,6 ms), recalculées à chaque appel. Les garder d'un appel à l'autre demande un objet de classification, à voir avec les booléens (lot 5) : dans une fusion de deux solides de 880 faces, la classification pèse 15 %, les volumes (`signed_volume`, `make_solid`) 85 % |
| C++23 | `interior_points` trie par `ranges::greater` sur la longueur et rend `found \| views::values \| ranges::to<vector>()` |
| Coût, trouvé au lot 5 : classer beaucoup de points dans une même face recalculait sa boîte (u, v) à chaque point, et projetait chaque point sur toutes ses arêtes. C'est ce que font les coïncidences partielles des booléens (36 échantillons par carreau de Bézier) : deux plaques superposées de 36 trous chacune, dont les trous se croisent, se fusionnaient en 31 s, dont 99 % dans `partial_of` | `FaceClassifier` garde la boîte (u, v) de la face, les boîtes de ses arêtes en 3D et celles de ses pcurves : une arête loin du point n'est pas projetée, une pcurve loin d'un rayon n'est pas croisée. Les coïncidences partielles et le test de faces confondues de la section (`on_surface`) s'en servent : la fusion des plaques passe à 6,3 s (16 trous : de 7,5 s à 1,8 s). Testé (`face_classifier` : mêmes réponses que `classify` sur une grille de 777 points de la plaque percée) |
