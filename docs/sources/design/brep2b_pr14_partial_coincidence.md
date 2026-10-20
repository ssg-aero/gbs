# Palier 2 bis — PR 14 : coïncidence partielle

| | |
|---|---|
| PR | branche `feat/partial-coincidence` |
| Document de référence | [brep_palier2bis.md](brep_palier2bis.md), § 1.2, § 7.1 et § 12.5 |
| Fichiers | `gbs-brep/partial.h` (nouveau), `gbs-brep/boolean.h`, `gbs-brep/section.h`, `gbs/intersect/marching.h`, `tests/tests_brep_boolean.cpp`, `tests/tests_intersect_marching.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/design/brep_palier2bis.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Deux faces confondues sur une partie seulement de leur recouvrement (une face
B-spline posée sur un plan sur une partie de son étendue, qui s'en écarte ou le
traverse au-delà) étaient refusées : la détection (§ 7.1 du palier) déclarait
les surfaces différentes dès qu'un point s'écartait, et le cheminement trouvait
alors des tangences partout dans la région commune (`NonGeneric`). C'était la
dernière limite fonctionnelle du palier 2 bis (§ 12.5).

## 2. Principe : la frontière suit les lignes de nœuds

Un carreau de Bézier est analytique : s'il est posé sur une autre surface
analytique sur un ouvert, il l'est sur tout le carreau. La région commune est
donc une **réunion de carreaux** de la face B-spline (la porteuse), et sa
frontière suit des **lignes de nœuds** : des isos exactes, indépendantes de la
tolérance, comme la fin géométrique d'une tangence (PR 13).

L'autre choix, la ligne où l'écart atteint la tolérance, serait plus général
mais ferait bouger le sommet avec la tolérance, et serait mal conditionné où
l'écart croît lentement (x³). Il a été écarté.

## 3. Étapes (`partial.h`, appelées par `boolean` avant la section)

1. **Détection** (`partial_coincidence`) : pour chaque paire de faces dont l'une
   est une B-spline à nœuds intérieurs, les carreaux de la porteuse (entre ses
   lignes de nœuds, dans sa boîte) sont échantillonnés (6 × 6) et projetés sur
   l'autre surface. Un carreau est **posé** (tous ses points en vis-à-vis dans
   la tolérance), **hors** (au moins un au-delà), ou **sans vis-à-vis**. Il y a
   coïncidence partielle quand des carreaux sont posés et d'autres hors.
2. **Frontière** : les segments de lignes de nœuds entre un carreau posé et un
   carreau hors, fusionnés bout à bout.
3. **Découpe** (`split_partial`) : la porteuse est découpée le long de ces
   segments (arêtes sur les isos de sa surface, `row_or_iso` ; pcurves en
   segments ; coupées au contour de la face). Le morceau posé est alors une
   coïncidence **totale** avec l'autre face, que le chemin existant traite
   (`imprint_coincident` imprime la frontière sur l'autre face).
4. **Croisements** : sur les carreaux hors de l'autre surface, le signe de
   l'écart est relevé. Les morceaux hors sont intersectés avec l'autre face
   **seulement sur les carreaux où il change de signe**
   (`SectionOptions::restricted`) ; aucun s'il garde un signe. Près de la
   frontière, l'écart reste sous la tolérance sur une bande (x < tol^(1/3) pour
   x³) où le cheminement trouverait des tangences partout ; les surfaces s'y
   touchent seulement le long de la frontière, déjà arête des deux.

## 4. Correction du cheminement trouvée en chemin

Le cas qui traverse a montré un défaut du cheminement, sans lien avec la
coïncidence : une courbe qui finit sur deux bords des boîtes (le flanc vertical
d'une boîte contre la face B-spline, la courbe allant de son bord inférieur à son
bord latéral) n'était **pas trouvée**. Le départ d'une paire de carreaux plats,
corrigé depuis leurs centres, gardait ses paramètres dans les boîtes en les
tronquant ; le pas poussant un paramètre au-delà d'un bord, l'itéré y restait
bloqué.

Correction (`Marcher::correct`, départs) : un paramètre sur une borne que le pas
pousse dehors est **tenu** sur la borne, et les autres sont recalculés
(ensemble actif) ; l'itéré glisse le long du bord. Un essai qui découpait plus
finement les carreaux plats changeait d'autres résultats et coûtait cher ; un
second départ sur la droite commune aux plans des carreaux ne suffisait pas.

Aussi : la réserve de sommets du booléen (`VertexPool`) saute les sommets que
l'empreinte des faces confondues a fusionnés.

## 5. Effet

| Cas | Avant | Après |
|---|---|---|
| bloc sous z = −max(x, 0)³, boîte posée dessus, x ∈ [−1, 0,5] (la face s'écarte) | `NonGeneric` | fusion 10,5 (un solide), commun vide, différence 7,5 |
| même bloc, z = +max(x, 0)³ (la face monte dans la boîte) | `NonGeneric` | fusion 11,47, commun 0,031, différence 8,47 |
| z = −(x₊³ − 16 (x − 0,5)₊³) : posée, s'écarte, puis traverse à x₀ ≈ 0,829 | `NonGeneric` | fusion, commun et différence aux volumes calculés |
| cheminement : flanc vertical et face B-spline, courbe d'un bord à l'autre | aucune courbe | la courbe, bouts à 1e-6 |

Coût : sur les paires communes de la comparaison OCCT, 2,24 s contre 2,19 s
sur master.

## 6. Tests

| Test | Cas | master |
|---|---|---|
| `partial_coincidence` (`tests_brep_boolean`) | détection : porteuse, frontière u = 0 sur v ∈ [0, 1], écart d'un signe ; faces sans coïncidence non détectées ; les deux cas du bloc (s'écarte, monte) avec leurs volumes ; le cas qui traverse : un carreau à intersecter (u ∈ [0,5, 1]), volumes par les primitives de f | `NonGeneric` |
| `start_on_a_border` (`tests_intersect_marching`) | la face B-spline sur u ∈ [0,5, 1] contre le flanc y = −0,8 : une courbe de (x₀, −0,8, 0) à (0,9, −0,8, 0,295) | **aucune courbe** |
| comparaison OCCT | `partial_coincidence`, `partial_coincidence_rising` (bloc extrudé d'un profil B-spline puis converti en NURBS) | — |

Dans la comparaison OCCT, gbs découpe le dessus du bloc à la ligne de nœuds
(une face de plus, `same_faces = false`). Pour le cas qui s'écarte, OCCT ne voit
pas le contact : sa fusion rend les deux opérandes tels quels (deux solides,
6 + 6 faces, volumes 7,5 et 3) ; gbs les réunit le long de la région commune
(un solide, vérifié dans `partial_coincidence`), d'où `same_solids = false` pour
cette paire. Les volumes et la validité sont comparés.

Les 495 tests C++, le binaire OCCT (29 cas) et les 54 tests Python passent.

## 7. Points à valider

1. La frontière sur les lignes de nœuds de la porteuse (§ 2).
2. L'échantillonnage des carreaux (6 × 6) : un croisement plus étroit qu'un pas
   d'échantillon dans un carreau dont l'écart semble d'un signe serait manqué.
3. Hors de portée : deux B-splines dont les lignes de nœuds se croisent dans la
   région commune (un carreau de la porteuse à moitié posé) ; un carreau hors
   adjacent à la frontière où l'écart change de signe (bande de tolérance) : le
   booléen reste refusé.
4. La correction du cheminement par ensemble actif (§ 4).
