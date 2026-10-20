# Palier 2 bis — PR 14 : coïncidence partielle

| | |
|---|---|
| PR | branche `feat/partial-coincidence` |
| Document de référence | [brep_palier2bis.md](brep_palier2bis.md), § 1.2, § 7.1 et § 12.5 |
| Fichiers | `gbs-brep/partial.h` (nouveau, dans `gbs-brep/brep`), `gbs-brep/boolean.h`, `gbs-brep/section.h`, `gbs/intersect/marching.h`, `tests/tests_brep_boolean.cpp`, `tests/tests_intersect_marching.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/design/brep_palier2bis.md`, `News.md` |
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
   lignes de nœuds, dans sa boîte) dont la boîte de Bézier atteint l'autre face
   sont échantillonnés (6 × 6) et projetés sur l'autre surface ; les autres sont
   sautés. Un carreau est **posé** (tous ses points en vis-à-vis dans la
   tolérance), **hors** (au moins un au-delà), ou **sans vis-à-vis**. Il y a
   coïncidence partielle quand des carreaux sont posés et d'autres hors.
2. **Frontière** : les segments de lignes de nœuds entre un carreau posé et un
   carreau hors, fusionnés bout à bout.
3. **Découpe** (`split_partial_coincidences`) : chaque porteuse est analysée
   contre **toutes** les faces de l'autre forme, puis découpée **une fois** le
   long de toutes les frontières (arêtes sur les isos de sa surface,
   `row_or_iso` ; pcurves en segments ; coupées au contour de la face). Un
   morceau est posé sur une face s'il contient un carreau posé pour elle (testé
   au centre des carreaux de la grille) : un morceau peut tenir des carreaux
   posés et des carreaux sans vis-à-vis, la frontière ne passant qu'entre posés
   et hors. Le morceau posé est alors une coïncidence **totale** avec l'autre
   face, que le chemin existant traite (`imprint_coincident` imprime la
   frontière sur l'autre face). Les restrictions qui nomment une face découpée
   plus tard (une porteuse de l'autre forme) passent à ses morceaux.
4. **Croisements** : sur les carreaux hors de l'autre surface, le signe de
   l'écart est relevé. Contre un **plan**, les distances au plan des pôles de
   Bézier du carreau (les coefficients de Bernstein de l'écart) le
   **certifient** ; contre une autre surface, les échantillons décident. Les
   morceaux hors sont intersectés avec l'autre face **seulement sur les
   carreaux où il change de signe** (`SectionOptions::restricted`) ; aucun s'il
   garde un signe. Près de la frontière, l'écart reste sous la tolérance sur une
   bande (x < tol^(1/3) pour x³) où le cheminement trouverait des tangences
   partout ; les surfaces s'y touchent seulement le long de la frontière, déjà
   arête des deux. La boîte d'un carreau croisant qui borde la frontière en est
   donc écartée, du plus petit décalage où l'écart dépasse 10 tol tout le long
   (un pli : 6e-5 ; une tangence : aucun, l'intersection reste `singular`).

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
(ensemble actif, repris tant qu'un paramètre s'ajoute : quatre passes au plus) ;
l'itéré glisse le long du bord. Un essai qui découpait plus
finement les carreaux plats changeait d'autres résultats et coûtait cher ; un
second départ sur la droite commune aux plans des carreaux ne suffisait pas.

Aussi : la réserve de sommets du booléen (`VertexPool`) saute les sommets que
l'empreinte des faces confondues a fusionnés ; une arête partagée par deux faces
confondues de l'autre opérande, imprimée par les deux, n'est plus une coupe en
double (`split_face` la refusait).

## 5. Effet

| Cas | Avant | Après |
|---|---|---|
| bloc sous z = −max(x, 0)³, boîte posée dessus, x ∈ [−1, 0,5] (la face s'écarte) | `NonGeneric` | fusion 10,5 (un solide), commun vide, différence 7,5 |
| même bloc, z = +max(x, 0)³ (la face monte dans la boîte) | `NonGeneric` | fusion 11,47, commun 0,031, différence 8,47 |
| z = −(x₊³ − 16 (x − 0,5)₊³) : posée, s'écarte, puis traverse à x₀ ≈ 0,829 | `NonGeneric` | fusion, commun et différence aux volumes calculés |
| cheminement : flanc vertical et face B-spline, courbe d'un bord à l'autre | aucune courbe | la courbe, bouts à 1e-6 |

Coût : sur les paires communes de la comparaison OCCT, 2,16 s contre 2,19 s
sur master.

## 6. Tests

| Test | Cas | master |
|---|---|---|
| `partial_coincidence` (`tests_brep_boolean`) | détection : porteuse, frontière u = 0 sur v ∈ [0, 1], écart d'un signe ; faces sans coïncidence non détectées ; les deux cas du bloc (s'écarte, monte) avec leurs volumes ; dessus **transposé** (frontière le long de u, `along_u`), **rationnel**, et sous une boîte depuis x = −0,3 (un morceau à carreaux posés et sans vis-à-vis) ; un solide **en L** dont deux faces du dessous sont posées sur la même porteuse ; un **croisement plus étroit que l'échantillonnage** (pli en 0, puis x ((x − ½)² − 1/200) : bande de 0,14) trouvé par les pôles, boîte écartée de la frontière, copeau du commun au volume calculé ; le cas qui traverse : un carreau à intersecter (u ∈ [0,5, 1]), volumes par les primitives de f | `NonGeneric` |
| `start_on_a_border` (`tests_intersect_marching`) | la face B-spline sur u ∈ [0,5, 1] contre le flanc y = −0,8 : une courbe de (x₀, −0,8, 0) à (0,9, −0,8, 0,295) | **aucune courbe** |
| comparaison OCCT | `partial_coincidence`, `partial_coincidence_rising` (bloc extrudé d'un profil B-spline puis converti en NURBS) | — |

Dans la comparaison OCCT, gbs découpe le dessus du bloc à la ligne de nœuds
(une face de plus, `same_faces = false`). Pour le cas qui s'écarte, OCCT ne voit
pas le contact : sa fusion rend les deux opérandes tels quels (deux solides,
6 + 6 faces, volumes 7,5 et 3) ; gbs les réunit le long de la région commune
(un solide, vérifié dans `partial_coincidence`), d'où `same_solids = false` pour
cette paire. Les volumes et la validité sont comparés.

Les 495 tests C++, le binaire OCCT (29 cas) et les 54 tests Python passent.

### 6.1 Revue de code

| Commentaire | Suite donnée |
|---|---|
| un morceau classé par un seul point intérieur : tombé sur un carreau sans vis-à-vis, la coïncidence était perdue sans refus | classé par la grille : posé s'il contient un carreau posé ; testé (boîte depuis x = −0,3, qui échoue sur la version précédente) |
| signe d'un carreau décidé par 36 échantillons ; contre un plan, les coefficients de Bernstein le certifient | fait pour un plan (primitive ou B-spline plane) ; testé par le croisement étroit (la mutation sans les pôles ne le voit pas). Les autres surfaces gardent l'échantillonnage (point à valider 2) |
| paires fixées avant toute découpe : les morceaux ne sont pas comparés aux autres faces | chaque porteuse analysée contre toutes les faces, découpée une fois ; restrictions reportées sur les morceaux des faces découpées ensuite ; testé (solide en L) |
| coût : 36 projections par carreau, pour tous les carreaux | carreaux dont la boîte de Bézier n'atteint pas l'autre face sautés ; temps des paires OCCT inchangé. Le pré-échantillon 2 × 2 n'a pas été ajouté |
| ensemble actif en une passe | repris tant qu'un paramètre s'ajoute (quatre passes au plus) ; le coin des boîtes n'a pas de test dédié |
| cas non couverts | ci-dessus. En chemin : un carreau croisant qui borde la frontière (le pli du croisement étroit) était `singular` à cause de la frontière elle-même ; sa boîte en est écartée (§ 3, étape 4) |

## 7. Points à valider

1. La frontière sur les lignes de nœuds de la porteuse (§ 2).
2. Contre une surface autre qu'un plan, l'échantillonnage des carreaux
   (6 × 6) : un croisement plus étroit qu'un pas d'échantillon serait manqué
   (contre un plan, les pôles le certifient).
3. Hors de portée : deux B-splines dont les lignes de nœuds se croisent dans la
   région commune (un carreau de la porteuse à moitié posé) ; un carreau
   croisant qui borde la frontière par une tangence (bande de tolérance) : le
   booléen reste refusé.
4. La correction du cheminement par ensemble actif (§ 4).
