# Palier 6 — PR 1 : onglets, extrémités à plus de trois arêtes, sphère à n arcs

| | |
|---|---|
| PR | branche `feat/brep-vertex-blends` |
| Document de référence | [brep_palier6_vertex_blends.md](brep_palier6_vertex_blends.md), § 2, § 3 et § 4.1 |
| Fichiers | `gbs-brep/blend.h`, `tests/tests_brep_vertex_blends.cpp` (nouveau), `tests/tests_brep_fillet.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier6_vertex_blends.md`, les notes `brep5_pr01_variable_geometry.md` et `brep5_pr02_variable_topology.md` (statuts), `News.md`, cette note |
| Statut | À valider (§ 9) ; revue de code : § 10 |

## 1. Objectif

Première PR du palier 6 : les sommets où se rencontrent des arêtes traitées
sans former de chaîne, sauf le raccord de sommet lui-même (PR 2 et 3).

- **Onglet de deux traitements** à un sommet de trois arêtes, la troisième
  vive : deux congés, deux chanfreins, convexes ou concaves (§ 4.2).
- **Trois chanfreins** à un sommet de trois arêtes (§ 4.3).
- **Extrémité libre** d'un traitement à un sommet de plus de trois arêtes
  (§ 4.4).
- **Sphère à n arcs** : les n congés de même rayon à un sommet de n arêtes
  dont les lignes des centres concourent (§ 4.5).

Le palier 4 refusait tous ces cas.

## 2. Sources

- Vida, Martin et Várady, « A survey of blending methods that use parametric
  surfaces » (*Computer-Aided Design* 26(5), 1994), § 5 : les coins où les
  congés se coupent sans raccord.
- Rossignac et Requicha, « Constant-radius blending in solid modelling »
  (1984) : le coin sphérique, ici à n faces.
- Steiner (formule du volume d'un convexe grossi par une boule) pour la
  référence de la pyramide arrondie (§ 6).
- OCCT (`BRepFilletAPI_MakeFillet`, `BRepFilletAPI_MakeChamfer`) en boîte
  noire.

## 3. API

Pas de fonction nouvelle : `fillet` et `chamfer` acceptent les configurations
du § 1. Rien de nouveau en Python ; les exemples de la galerie gagnent cinq cas
(§ 6).

## 4. Démarche

### 4.1 Classement

Pour chaque extrémité d'un traitement à un sommet v (n arêtes à v, k traitées
dont celle-ci) :

| Cas | Traitement |
|---|---|
| k = 1, n = 3 | extrémité libre (palier 4) |
| k = 1, n > 3 | extrémité libre sur n − 2 faces (§ 4.4) |
| k = 2, arêtes tangentes à v | chaîne (palier 4 et 5) |
| k = 2, non tangentes, n = 3 | onglet (§ 4.2) |
| k = n, congés | sphère (§ 4.5) ; refusée si les lignes des centres ne concourent pas (rayons différents : PR 3) |
| k = 3 = n, chanfreins | onglet des trois (§ 4.3) |
| autres | refus |

La tangence est celle de `law_chains` : la direction sortante d'une arête à v
opposée à celle de l'autre (à 10⁻⁶ près).

### 4.2 Onglet de deux traitements

Les traitements A et B ont une face commune F ; c est la troisième arête, vive,
entre leurs autres faces Fa et Fb.

1. **Croisement X** des contacts de A et B sur F (Gauss-Newton sur leurs deux
   paramètres, `crossing`).
2. **Rencontre avec c** : le contact de A sur Fa rencontre Fb (prolongée) sur c,
   celui de B sur Fb rencontre Fa sur c.
3. **Le plus étroit.** Si le contact de A sur Fa rencontre la surface de B,
   sur le traitement B lui-même (son s entre 0 et 1 : pas sur le dos d'un
   cylindre) et plus loin de v que c, A est le plus étroit : il s'arrête là, en
   P. Sinon on essaie B. Si aucun ne coupe l'autre, les deux contacts doivent
   rencontrer c au même point (même largeur), P.
4. **Courbe d'onglet** entre A et B, de P à X, à travers le plus étroit
   (`cross_curve`, la machinerie de l'extrémité libre du palier 4 généralisée à
   un intervalle de s et à une autre surface que la face d'extrémité).
5. **Le plus large** continue au-delà de P, coupé par la face de l'autre contact
   du plus étroit (au-delà de c pour lui) : courbe d'extrémité de P jusqu'à son
   autre contact sur c. Cette courbe, sur la face, se chaîne au contact du plus
   étroit (`face_cuts` dans `cut_face`) ; c est coupée là.

**Point double en X.** Les deux congés sont tangents à F en X, donc tangents
entre eux : leur intersection y a un point double (les deux ellipses de deux
cylindres de même rayon s'y croisent). Gauss-Newton part mal de là. La courbe
est donc cheminée depuis P, régulier, et son dernier point X est donné, non
résolu (`cross_curve(…, last)`). Les chanfreins, qui ne sont pas tangents aux
faces, n'ont pas ce point double.

L'extrémité du plus large a deux arêtes (courbe d'onglet, courbe sur la face) :
`End::pieces`, dans l'ordre du contact 0 au contact 1, que `blend_face`
parcourt.

### 4.3 Trois chanfreins

Chaque paire a une face commune et un croisement X. Le point Q commun aux trois
surfaces est résolu par Gauss-Newton sur les (t, s) des trois (six inconnues),
à partir de la moyenne des X ; il doit être dans chacun (s entre 0 et 1). Les
trois lignes d'onglet vont de chaque X à Q (`cross_curve`, son bout vérifié en
Q). Chaque chanfrein a deux pièces : de son contact 0 à Q, de Q à son
contact 1.

### 4.4 Extrémité libre à plus de trois arêtes

Les arêtes et faces autour de v sont mises dans l'ordre (`ring`) à partir de
l'arête traitée et de la face de son contact 0. Le contact 0 rencontre la face
suivante sur l'arête e₁ ; le traitement rencontre chaque arête intérieure eᵢ
(Gauss-Newton courbe-surface, puis son s sur le traitement par `locate`) ; le
contact 1 rencontre la face précédente sur eₙ₋₁. Les s doivent croître de 0 à 1.
Chaque arête est coupée en son point ; entre deux points, la courbe sur la face
est cheminée (`cross_curve`). Chaque face remplace sa suite d'arêtes autour de v
par sa pièce (`replace_run`), comme la troisième face du palier 4.

### 4.5 Sphère à n arcs

`corner` prend n congés au lieu de trois : Gauss-Newton sur n paramètres
(n − 1 égalités des lignes des centres), n faces partagées deux à deux, n pieds.
Le repère de la sphère (`corner_frame`) teste l'intérieur du polygone sphérique
des pieds, mis dans l'ordre des faces autour de v (`ring`), par le côté de son
milieu pour chaque côté (le triangle du palier 4 prenait son troisième sommet).

## 5. Refus

Nouveaux, ou changés :

- deux arêtes traitées non tangentes à un sommet de plus de trois arêtes ;
- des arêtes traitées à un sommet où d'autres restent vives, hors l'onglet de
  deux ; une chaîne tangente traversant un sommet avec d'autres arêtes
  traitées (deux quelconques d'entre elles tangentes) ;
- un onglet d'une arête convexe et d'une arête concave ;
- un onglet dont les contacts ne se croisent pas sur la face commune dans les
  traitements prolongés, dont les traitements ne se coupent pas ou ne
  rencontrent pas c au même point, dont l'autre contact du plus large (des
  deux, s'ils rencontrent c au même point) ne rencontre pas c dans son
  traitement prolongé, ou dont c est coupée à son autre sommet ;
- trois chanfreins dont le point commun n'est pas dans chacun ;
- une extrémité libre qui rencontre les arêtes autour du sommet dans le
  désordre ;
- trois congés de rayons différents au sommet (lignes des centres non
  concourantes) : toujours refusés, jusqu'au raccord de sommet (PR 3).

Retirés : « more than three edges » (extrémité libre), « three that are not
fillets » (chanfreins au coin), « not tangentially » à un sommet de trois
arêtes.

## 6. Tests

`tests_brep_vertex_blends` (10 cas), volumes de référence établis à part :

| Cas | Référence | Écart |
|---|---|---|
| onglet de deux congés au coin d'une boîte, rayons (0,2 ; 0,2), (0,15 ; 0,3), (0,3 ; 0,15), (0,25 ; 0,1) | prismes (1 − π/4) r² L, moins le volume ôté deux fois : ∫ h₁(q) h₂(q) dq, h la profondeur d'un congé à la distance q de la face commune (Gauss-Legendre, q = u²) | 10⁻⁹ |
| onglet de deux chanfreins | 1 − d² + d³/3 | 10⁻¹² |
| trois chanfreins | 1 − 3d²/2 + 3d³/4 (inclusion-exclusion des trois prismes), le sommet commun en (1 − d/2)³ | 10⁻¹² |
| onglet concave (poche) | 7 + 2 (1 − π/4) r² − ∫ h² dq | 10⁻⁹ |
| extrémité libre à quatre arêtes (dessus coupé en diagonale) | 2 − (1 − π/4) r² | 10⁻⁹ |
| pyramide aux huit arêtes congées | Steiner : la pyramide réduite de (ρ − r)/ρ autour du centre de sa sphère inscrite, grossie de r : V + r A + r² Σ ℓ φ / 2 + 4πr³/3 (φ l'angle extérieur de chaque arête) | 10⁻⁸ |
| un côté de la pyramide | valide, 6 faces, volume entre la pyramide et la pyramide moins le prisme du côté | — |
| onglet à un sommet de 120° et de 150° (prisme sur un quadrilatère) | le même solide avec un prolongement de trois rayons, à 10⁻⁹ ; à 150°, l'autre contact rencontre l'arête vive à 1,73 r au-delà du sommet | 10⁻⁹ |
| refus | trois congés de rayons différents au coin ; deux côtés de la pyramide à son sommet (quatre arêtes, non tangents) ; une pyramide sur une base en flèche, dont un côté rencontre les arêtes autour du sommet dans le désordre ; trois chanfreins de 0,05, 0,05 et 0,6 au coin d'une boîte (distances par arête, par `blend_detail::blend`), leur point commun (0,7 ; 0,7 ; 1,25) hors de deux d'entre eux ; une arête convexe et une concave au coin intérieur d'un L ; modèle inchangé | — |
| chaîne traversant un sommet | `chain_through` sur des arêtes nues, dans tous les ordres : à faces planes, aucune arête qu'on puisse traiter n'arrive à un tel sommet (les faces de part et d'autre d'une chaîne droite contiennent sa droite, une arête entre deux d'entre elles est plate) | — |

`tests_brep_fillet` : les trois tests qui attendaient un refus attendent
maintenant le résultat (onglet de deux congés à volume exact, coin
chanfreiné, extrémité à quatre arêtes déplacée dans le nouveau fichier).

OCCT (`vertex_blends`) : onglet de deux congés et de deux chanfreins, onglet
concave, pyramide aux huit arêtes (sphère à quatre arcs) et un de ses côtés :
résultats de gbs relus valides, même volume à 10⁻⁶, mêmes nombres de faces.

Exemples : `fillet_mitre_box`, `fillet_mitre_pocket`, `chamfer_corner`,
`fillet_pyramid`, `fillet_pyramid_side`, volumes exacts sauf le dernier.

## 7. Comparaison à OCCT : le coin de trois chanfreins

OCCT ne fait pas l'onglet des trois chanfreins : il coupe le coin par le plan
des trois croisements des contacts, une face triangulaire de plus. Mesuré :
10 faces contre 9, et un volume plus petit de d³/12 exactement, le tétraèdre
entre ce plan et le point commun des trois chanfreins (à 10⁻⁹, vérifié par le
test). Ce n'est donc pas une erreur de l'un ou de l'autre, mais deux choix. Le
§ 8 du palier (question 3) a retenu l'onglet ; le comportement de CATIA et
SolidWorks n'est pas vérifié sur ce point (question 1 du § 9 ci-dessous).

## 8. Limites

- Onglet seulement de deux traitements à un sommet de trois arêtes ; à plus
  de trois arêtes, deux traitements non tangents sont refusés.
- Trois traitements à un sommet de trois arêtes : chanfreins (onglet) ou
  congés de même rayon (sphère) ; les congés de rayons différents attendent le
  raccord de sommet (PR 3), comme un mélange de congés et de chanfreins (une
  même opération ne les mêle pas aujourd'hui).
- Coins mêlés (convexes et concaves) : refusés, en sphère comme en onglet.
- L'extrémité libre à plus de trois arêtes suppose que le traitement
  rencontre les arêtes autour du sommet dans l'ordre.

## 9. Points à valider

1. **Coin de trois chanfreins en onglet** (point commun, pas de facette),
   alors qu'OCCT pose une facette triangulaire par les trois croisements. À
   confirmer d'après CATIA ou SolidWorks ; passer à la facette ne changerait que
   `mitre3`.
2. **Courbe d'onglet cheminée depuis P**, son point double X donné.
3. **Le plus étroit choisi par son contact** : il rencontre l'autre traitement
   sur le traitement même, au-delà de l'arête vive ; sinon les deux contacts
   rencontrent l'arête vive au même point.
4. **`End::pieces`** pour les extrémités de plusieurs arêtes, plutôt qu'une
   arête composite.
5. **`cross_curve`** comme machinerie commune (extrémité libre, onglet, face
   d'extrémité à plusieurs pièces).
6. **Pyramide arrondie validée par la formule de Steiner.**

## 10. Revue de code

Revue du 20 nov. 2026, 6 points, tous corrigés :

1. **Onglets et extrémités à plus de trois arêtes prolongés** comme les
   extrémités libres : la boucle qui prolonge un traitement 2, 4 puis 8 fois
   teste aussi, avant toute coupe, l'autre contact d'un onglet sur la face au
   delà de l'arête vive, et les contacts d'une extrémité libre sur les faces
   voisines autour du sommet (`reaches`). Le contact du plus étroit n'a plus à
   rencontrer c dans son traitement : il ne sert qu'à partir. Sonde : un
   onglet à un sommet de 150° était refusé avec le prolongement par défaut ;
   il passe, le même solide qu'avec un prolongement de trois rayons (test
   `mitre_obtuse`).
2. **Chaîne traversant un sommet** : la tangence est testée entre deux
   quelconques des traitements au sommet, et non plus seulement avec le
   premier parcouru ; le refus ne dépend plus de l'ordre des arêtes.
3. **Onglet mêlé** convexe et concave refusé comme tel (« a mitre of convex
   and concave blended edges (not handled) ») au lieu d'un motif de largeurs.
   Testé sur le coin intérieur d'un L.
4. **Trois chanfreins** : les sommets X et Q élargis aux bouts des lignes
   d'onglet acceptés.
5. **Tests** : refus à un sommet de quatre arêtes et de l'onglet mêlé, modèle
   inchangé ; onglets à 120° et 150°. Puis les trois refus restants :
   - extrémité dans le désordre : trouvée par une sonde de 150 pyramides au
     hasard, chaque côté congé seul (600 congés : 2 dans le désordre, sur des
     bases en flèche ; 512 acceptés, tous valides), puis fixée à une base en
     flèche aux valeurs rondes ;
   - point commun de trois chanfreins hors de l'un : inatteignable par
     `chamfer`, une seule distance au coin (300 coins trièdres au hasard :
     jamais) ; atteint par `blend_detail::blend` avec une distance par arête
     (0,05, 0,05, 0,6) ;
   - chaîne traversant un sommet : impossible à faces planes avec des arêtes
     qu'on puisse traiter (voir § 6) ; le test de tangence sorti en
     `tangent_at` et `chain_through` (repris par `law_chains`), testé sur des
     arêtes nues dans tous les ordres.
6. **Doc de `blend()`** : la sphère de n arcs, `mitre`, `mitre3`.
