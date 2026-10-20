# Palier 2 bis — PR 11 : tangence le long d'une courbe

| | |
|---|---|
| PR | branche `feat/tangent-curves` |
| Document de référence | [brep_palier2bis.md](brep_palier2bis.md), § 5.3 et § 12.5 |
| Fichiers | `gbs/intersect/marching.h`, `tests/tests_intersect_marching.cpp`, `tests/tests_brep_boolean.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/design/brep_palier2bis.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Hors des formules fermées, deux surfaces tangentes le long d'une courbe (un
cylindre couché sur une face B-spline, un tore dans un cylindre, une sphère
dans un cylindre de son rayon) donnaient une intersection `singular`, et les
booléens `NonGeneric`. C'était la dernière grande limite du palier 2 bis
(§ 12.5).

## 2. Pourquoi le cheminement s'arrêtait

Le cheminement suit n1 × n2 et corrige par Newton sur S1 − S2 = 0. Sur une
courbe de tangence, n1 × n2 est nul et le jacobien de S1 − S2 est de rang 2 en
tout point : ni la direction ni la correction ne sont définies.

## 3. Le système de la tangence

Inconnues (u1, v1, u2, v2). Quatre équations :

| Équations | Sens |
|---|---|
| (S1 − S2) · S2_u = 0, (S1 − S2) · S2_v = 0 | S2 est le pied de S1 sur la seconde surface |
| n2 · S1_u = 0, n2 · S1_v = 0 | les normales sont parallèles |

Sur une courbe de tangence, ses solutions forment une courbe : le système est de
rang 3, et la distance |S1 − S2| y est dans la tolérance. Il est bien
conditionné, contrairement à S1 − S2 = 0. Il garde aussi le pied, si bien
qu'une tangence à 1e-9 près (surfaces écartées ou enfoncées d'un arrondi) donne
une seule courbe, et non deux courbes transverses à √tol l'une de l'autre.

Le jacobien demande les dérivées secondes des surfaces ; la dérivée de la
normale unitaire est (N_u − n (n · N_u)) / |N|, avec N_u = S_uu × S_v + S_u × S_uv.
Les lignes sont mises en longueurs : les équations du pied divisées par |S2_u|,
|S2_v|, celles des normales (des sinus) multipliées par la longueur d'un pas.

**Direction.** Au premier ordre, toute direction du plan tangent commun garde
S1 = S2. La courbe suit la direction où les surfaces se courbent pareil : le
vecteur propre de la plus petite valeur propre de II1 − II2, différence des
secondes formes fondamentales sur une base commune du plan tangent. Sans
direction (II1 − II2 nulle dans toutes les directions : surfaces confondues,
contact ombilical), rien n'est suivi.

**Correction.** Chaque point est corrigé sur le plan normal au pas, comme dans
le cheminement transverse : cinq équations, quatre inconnues, par moindres
carrés (SVD). Un itéré qui sort d'une surface est rendu à l'appelant, qui se
pose sur le bord de la boîte, l'un des paramètres fixé sur la borne.

## 4. Intégration au cheminement

`TangentMarcher` (`marching.h`) suit une branche dans les deux sens : arrêt sur
le bord d'une boîte, retour au départ (courbe fermée, les paramètres revenus
aussi, comme dans `Marcher` : à travers une couture, la branche se pose sur la
couture), ou perte de la tangence (`singular`). `marching` l'essaie à chaque
départ où les surfaces sont tangentes (sinus des normales sous `angular_tol`),
si les deux surfaces ont des dérivées secondes (`has_second_derivatives()` ;
sinon, le résultat d'avant, `singular`). La plus petite valeur propre de
II1 − II2 certifie une courbe de contact : nulle le long de la courbe, du même
ordre que l'autre en un contact isolé (une sphère sur un plan), qui n'est pas
suivi. Le départ est classé : pas une courbe (chemin habituel), déjà suivi, ou
branche. La
courbe est construite comme les autres, densifiée par le même correcteur, et
marquée `tangent` : la section et les booléens la traitent comme les tangences
des formules fermées (arête qui peut ne pas découper une face qu'elle touche).

Les départs près d'une tangence s'étalent de part et d'autre, à √(tol R) : le
test « déjà suivie » porte sur le départ remis sur la tangence, et la boîte de
la branche est grossie de la flèche de ses cordes, comme pour les autres.

## 5. Pose au contact d'un bord touché

Les copeaux tangents d'une boîte B-spline contre un cylindre ont montré la
limite notée à la PR 10 : sur la face z = 3, le cercle touche le bord x = 5 au
point même où finit la ligne de tangence. La pose tangente s'arrêtait à
quelques 1e-6 du contact, les sommets ne se rejoignaient pas (« a section edge
does not split a face »).

Quand toutes les poses sont tangentes, la pose gardée est maintenant amenée sur
le contact (`Marcher::touch`) : S1 = S2 et la tangente le long de l'iso,
(n1 × n2) · (n × S_o) = 0, quatre équations bien conditionnées où la courbe
s'écarte de l'iso (la dérivée de la dernière équation par différences
centrées, qui ne demandent que les dérivées premières : `touch` sert pour toutes
les surfaces). Le contact est mis sur la borne s'il en est à la tolérance.

## 6. Tests

| Test | Cas | master |
|---|---|---|
| `tangent_along_a_curve` (`tests_intersect_marching`) | flanc B-spline rationnel d'un cylindre couché sur le plan y = 1 : une courbe tangente, la génératrice x = 0, longueur 2 ; tore dans un cylindre coaxial, tangents le long du cercle de rayon 3 : longueur 6π ; les mêmes le long d'une couture (l'équateur extérieur, couture v = 0 du tore ; la génératrice x = 1, couture u = 0 du flanc) | `singular` |
| `singular_points_and_limits` | la sphère dans un cylindre de son rayon, qui y vérifiait le refus, donne l'équateur, tangent, longueur 2π | `singular` |
| `along_a_border` | plan borné par la boîte que le cercle touche : les huit bouts des quatre quarts sont les contacts à 1e-12 | bouts à quelques 1e-6 |
| `bspline_tangency` (`tests_brep_boolean`, ex-`non_generic`) | cylindre couché sur un cube à faces B-spline : fusion 64 + 250π (deux solides), commun vide, différence 64 ; le refus d'une intersection incomplète (`max_pairs`) et la remise en état du modèle restent vérifiés | `NonGeneric` |
| `tangency_limits` | face B-spline z = −max(x, 0)³ (exacte, nœud en 0) sous un cylindre couché : la tangence cesse en x ≈ 0, branche `singular`, de x = −1 à la fin du contact ; vingt départs au-delà de la couture d'une tangence fermée (tore dans un cylindre) : jamais fermée, paramètres monotones, posée sur la couture ; sphère décalée (`OffsetSurface`, sans dérivées secondes) dans un cylindre de son rayon : `singular` sans exception ; sphère sur un plan : un contact isolé, aucune courbe | exception pour l'offset |
| comparaison OCCT | `cylinder_lying_nurbs`, `seam_on_face_nurbs` (boîte convertie en NURBS par OCCT, cylindre primitif) ; `tangent_slivers_nurbs`, comparé à OCCT sur les primitives : OCCT ne découpe pas sa face NURBS le long de la tangence (un solide au lieu de deux copeaux), gbs donne les nombres de solides et de faces d'OCCT sur les primitives | `NonGeneric` |

Le cylindre reste primitif dans les paires OCCT : converti en NURBS par OCCT, il
est relu par gbs avec des pcurves à 5e-6 de ses arêtes, que la vérification de
lecture refuse (sans lien avec cette PR).

Les 492 tests C++, le binaire OCCT (32 paires, 29 cas) et les 54 tests Python passent.

### 6.1 Revue de code

| Commentaire | Suite donnée |
|---|---|
| dérivées secondes évaluées sans le garde `has_second_derivatives()` : `OffsetSurface` lève une exception | suivi de la tangence seulement si les deux surfaces les ont ; testé avec une sphère décalée (la mutation sans garde lève l'exception) |
| fermeture en 3D sans le contrôle des paramètres de `Marcher::march` | contrôle ajouté. Testé sur vingt départs au-delà d'une couture, mais **non épinglé** : sans le contrôle, la pose sur la couture vient toujours d'abord (la fenêtre où la fermeture passerait est de quelques 1e-4 sur le pas) |
| plus petite valeur propre calculée mais inutilisée : un contact isolé est essayé jusqu'au pas minimal | refus si elle dépasse `angular_tol` fois l'autre ; testé (sphère sur un plan). Le coût évité n'est pas mesurable sur ce cas (0,07 s avec ou sans) |
| six appels à `value()` là où `derivatives(u, v, 2, 2)` fait une passe | `derivatives` : une passe sur les B-splines. `direction` au point suivant refait encore l'évaluation du dernier itéré de `correct` : gardé, une évaluation par pas |
| différences finies décentrées dans `touch` | différences centrées (erreur en h²), décentrées vers l'intérieur contre une borne seulement. Les dérivées exactes de g demanderaient celles des deux normales : `touch` sert aussi aux surfaces sans dérivées secondes, une pose par bout de branche |
| résultat à trois états caché dans un `optional` et un vecteur vide | `TangentBranch::Kind` : `not_a_curve`, `traced`, `branch` |
| cas non couverts | `tangency_limits` ci-dessus |

## 7. Points à valider

1. Le système du pied et des normales parallèles, avec la direction par la
   différence des secondes formes fondamentales.
2. Une tangence qui cesse à l'intérieur des faces (les surfaces s'écartent ou
   se croisent au-delà) termine la branche en `singular` : le booléen reste
   `NonGeneric` dans ce cas (testé au niveau de l'intersection, et du booléen
   à la PR 12).
3. La pose au contact d'un bord touché (§ 5), qui lève le point 2 de la PR 10.
4. Pas de suivi de tangence pour une surface sans dérivées secondes
   (`OffsetSurface`) : `singular` comme avant.

## 8. Suite

- **Fait à la PR 12**
  ([brep2b_pr12_tangency_ending_in_a_face.md](brep2b_pr12_tangency_ending_in_a_face.md)) :
  un booléen sur un solide B-spline dont une face porte une tangence qui cesse à
  l'intérieur (la face z = −max(x, 0)³ sous un cylindre couché) est refusé
  (`NonGeneric`), et le modèle est remis en état. Jusque-là, seul le `singular`
  de l'intersection était testé (`tangency_limits`).
