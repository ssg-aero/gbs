# Noyau BREP natif — conception du palier 6 : raccords de sommet quelconques

| | |
|---|---|
| Statut | Questions du § 8 validées (20 nov. 2026) ; chaque PR valide sa note et a sa note d'illustration |
| Périmètre | ce qui se passe à un sommet où se rencontrent des arêtes congées ou chanfreinées sans former une chaîne tangente : arêtes du sommet arrondies en partie (onglet), chanfreins, extrémité libre à un sommet de plus de trois arêtes, coin sphérique de n congés, raccord de sommet (*vertex blend*) entre des rayons différents, des lois variables, plus de trois arêtes, des arêtes convexes et concaves, avec retraits (*setback*) ; Python, galerie, comparaison à OpenCascade |
| Dépend de | palier 4 : congés et chanfreins (`fillet.h`, `blend.h`), coin sphérique (PR 7) ; palier 5 : lois de rayon (`radius_law.h`), morceaux ; palier 2 : intersection de surfaces (`gbs/intersect/surface_surface.h`) |
| Hors périmètre | chanfreins à distance variable, sections coniques ou à corde constante, congés de faces (*face blend*), raccords de sommet G2 |

Ce document décrit *quoi* construire et *pourquoi*, comme
[brep_palier4.md](brep_palier4.md) et
[brep_palier5_variable_fillets.md](brep_palier5_variable_fillets.md). Chaque PR
a sa note d'architecture `brep6_prNN_*.md`. C'est le deuxième point de la
feuille de route du palier 4 (§ 9), et ce qui restait au bilan du palier 5
(§ 9) : rayons différents au coin, rayon décroissant vers un coin, plus de
trois arêtes.

**Provenance.** Les algorithmes viennent de la littérature et du comportement de
CATIA et SolidWorks, la cible de qualité ; chaque note cite ses sources.
OpenCascade ne sert que de comparaison en boîte noire dans les tests.

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Cas au sommet | trois traitements, choisis par la configuration (§ 2) : **onglet** si une arête du sommet n'est pas traitée, ou s'il y a un chanfrein ; **sphère** si elle est exacte (congés de même rayon au sommet, lignes des centres concourantes, pas de retrait) ; **raccord de sommet** sinon | un mode d'angle choisi par l'utilisateur à chaque sommet (plus tard, s'il le faut) |
| Onglet | découpe **locale** : deux traitements voisins autour du sommet coupés le long de l'intersection de leurs surfaces, la face au-delà d'une arête non traitée coupant chacun comme une extrémité libre ; les points où trois surfaces se rencontrent par Gauss-Newton | booléens entre le solide et les volumes ôtés ou ajoutés par chaque congé (contacts tangents tout le long : le cas le plus difficile des booléens) |
| Raccord de sommet | **une seule face**, une B-spline bicubique découpée (pcourbes droites dans un domaine polygonal), d'**énergie de flexion minimale** sous des contraintes **G1 linéaires** le long des bords (la surface transfinie de Kato, prévue comme guide, écartée par la PR 2) | n carreaux à quatre côtés autour d'un point central (n faces par coin, G1 entre eux à refaire) ; carreau de Gregory rationnel (pas de type de surface nouveau dans gbs) |
| Bords du raccord | les **arcs d'extrémité** des congés (sections aux retraits) et, sur chaque face entre deux congés, une **courbe de raccord** tangente aux deux contacts (Hermite dans les paramètres de la face), réduite à un point quand les contacts se croisent au bout des deux congés | bords imposés par l'utilisateur |
| Retraits | **automatiques et minimaux** : chaque congé s'arrête à la plus éloignée des stations où ses contacts croisent ceux de ses voisins ; un **retrait par arête** au choix de l'utilisateur (`FilletOptions::setbacks`) l'allonge, pour un coin plus doux (coin « valise ») | retrait fixe en rayons |
| Sphère | gardée quand elle est exacte, **généralisée à n arcs** (sommet de pyramide) | toujours le raccord (une surface approchée là où la sphère est exacte) |
| Continuité | **G1** le long de tous les bords, sous `angular_tol` ; position sous `approx_tol`, l'écart absorbé dans la tolérance des arêtes, comme les coins de la PR 7 du palier 4 | G2 (courbure continue) : options de CATIA, plus tard |
| Surface à n côtés | dans le **noyau gbs** (`gbs/nsided.h`), utilisable hors BREP et en Python ; le BREP lui fournit les rubans | dans `gbs-brep` seulement |
| Coins mêlés | convexes et concaves au même sommet : **raccord de sommet** si tout est congé, onglet sinon ; refus documentés là où les contacts ne se croisent pas | refus de tout coin mêlé (situation du palier 4) |

---

## 1. Objectifs et non-objectifs

### 1.1 Objectifs

Toutes les configurations courantes d'un sommet, que le palier 4 refuse :

- **deux arêtes congées** à un sommet de trois arêtes, la troisième vive (coin
  de boîte arrondi sur deux arêtes) : les deux congés se coupent en onglet ;
- **chanfreins** à un coin (trois chanfreins d'une boîte, chanfrein et congé
  mêlés) ;
- **extrémité libre** d'un congé à un sommet de plus de trois arêtes (sommet
  d'une pyramide, une seule arête congée) : coupée par plusieurs faces ;
- **sphère de n congés** de même rayon quand leurs lignes des centres
  concourent (sommet d'une pyramide régulière à quatre faces) ;
- **raccord de sommet** de n congés quelconques : rayons différents
  (r₁, r₂, r₃ au coin d'une boîte), lois variables (rayon décroissant vers le
  coin, refusé au palier 5), quatre arêtes ou plus, arêtes convexes et
  concaves mêlées (coin d'une équerre), avec ou sans retraits ;
- Python, galerie, comparaison à `BRepFilletAPI_MakeFillet` et
  `BRepFilletAPI_MakeChamfer`.

C'est ce que CATIA (« Edge Fillet », options *Blend corner* et *Setback
distance*) et SolidWorks (congés à retrait, *Setback fillet*) proposent ; le
coin sans retrait est leur comportement par défaut.

### 1.2 Non-objectifs

- Raccords G2 (courbure continue aux bords du raccord).
- Chanfreins à distance variable, sections coniques : feuille de route.
- Congés de faces (*face blend*).
- Un mode d'angle imposé par l'utilisateur (onglet là où un raccord est
  possible, ou l'inverse).

---

## 2. Classement des sommets

À un sommet v de n arêtes, dont k traitées (congé ou chanfrein) et qui ne
forment pas une chaîne tangente traversant v :

| Configuration | Traitement | Palier 4-5 |
|---|---|---|
| k = 1, n = 3 | extrémité libre, coupée par la troisième face | fait |
| k = 1, n > 3 | extrémité libre, coupée par les n − 2 faces qui ne bordent pas l'arête (§ 3.3) | refusé |
| k = 2, tangentes | chaîne | fait |
| 2 ≤ k < n | **onglet** (§ 3) | refusé |
| k = n, un chanfrein au moins | **onglet** | refusé |
| k = n, congés, même rayon au sommet, lignes des centres concourantes, pas de retrait | **sphère** à n arcs (§ 4.1) | fait pour n = 3 et rayon croissant vers le coin |
| k = n, congés, autres cas | **raccord de sommet** (§ 4) | refusé |

Une chaîne tangente qui traverse v avec d'autres arêtes traitées au sommet (un
té) est un cas d'onglet ou de raccord : la chaîne y compte pour deux arêtes.

---

## 3. Onglet

### 3.1 Principe

Autour de v, les arêtes et les faces alternent : e₁, F₁, e₂, F₂, …, eₙ, Fₙ
(Fᵢ entre eᵢ et eᵢ₊₁). Chaque arête traitée porte la surface de son
traitement, prolongée au-delà de v comme aujourd'hui. Chaque arête non traitée
reste vive. Le coin est découpé paire par paire :

- **deux traitements voisins** (eᵢ et eᵢ₊₁ traitées, de part et d'autre de Fᵢ)
  se coupent le long de l'intersection de leurs surfaces. Elle part du point où
  leurs contacts se croisent sur Fᵢ (Gauss-Newton sur les deux pcourbes) et se
  chemine (`gbs/intersect/surface_surface.h`, même départ et mêmes pas que les
  sections du palier 2) jusqu'au point suivant ;
- **un traitement voisin d'une arête vive** (eᵢ traitée, eᵢ₊₁ non) est coupé par
  la face au-delà, Fᵢ₊₁, comme une extrémité libre : sa courbe d'extrémité est
  cheminée sur lui ;
- les **points de rencontre** de trois surfaces (deux traitements et une face,
  ou un traitement et deux faces sur l'arête vive) sont résolus par
  Gauss-Newton à trois inconnues, à partir des bouts des courbes déjà
  cheminées. L'arête vive est raccourcie jusqu'à son point.

À un coin convexe, le résultat est l'intersection des solides traités chacun
seul (localement) ; à un coin concave, leur réunion. La découpe locale les
calcule sans booléens : les surfaces des traitements se coupent
transversalement, alors que les contacts, tangents aux faces tout le long,
seraient le cas le plus difficile d'un booléen.

### 3.2 Chanfreins

Un chanfrein se coupe comme un congé : plan ou surface réglée contre son
voisin. Trois chanfreins au coin d'une boîte se rencontrent en un point (trois
faces, trois arêtes d'onglet), à toutes distances tant que ce point est dans
chacun. Pas de facette ajoutée au coin : c'est l'intersection des solides
chanfreinés. OCCT, lui, coupe le coin par le plan des trois croisements des
contacts (une face triangulaire de plus, le tétraèdre d³/12 en moins pour trois
chanfreins de distance d) ; le comportement de CATIA et SolidWorks reste à
vérifier sur ce point (PR 1, § 7 de sa note).

### 3.3 Extrémité libre à plus de trois arêtes

Un traitement seul à un sommet de n > 3 arêtes est coupé successivement par les
n − 2 faces qui ne bordent pas son arête, dans l'ordre autour de v : sa courbe
d'extrémité est faite de n − 2 morceaux, un par face, joints aux points où le
traitement rencontre deux faces voisines (sur leur arête commune, raccourcie
jusque-là). Une face que le traitement ne rencontre pas avant la suivante est
sautée.

### 3.4 Refus

- deux surfaces voisines qui ne se coupent pas, ou dont l'intersection ne
  rejoint pas le point suivant dans le traitement prolongé (un traitement plus
  large que la face entre eux) ;
- une intersection tangente (deux congés de même rayon dont les surfaces se
  touchent sans se couper : c'est un cas de sphère ou de chaîne) ;
- un point de rencontre au-delà d'une arête de face voisine.

---

## 4. Raccord de sommet

### 4.1 Sphère à n arcs

La sphère de la PR 7 du palier 4 se généralise : n lignes des centres qui
concourent en un point c à r des n faces, chaque congé arrêté sur son arc par
c (grand cercle), la face de sphère bornée par les n arcs. Même repère
(`corner_frame`, pôles et couture loin du polygone sphérique), mêmes contrôles
(résidu de Gauss-Newton sur les n − 1 rencontres, contacts communs aux pieds).
Avec un rayon croissant vers le coin, les bouts sont des petits cercles comme
au palier 5.

### 4.2 Retraits et bords

*Complété par la PR 3 ([note](brep6_pr03_vertex_blends.md), § 4.3 et § 4.4) :
où deux contacts finiraient à moins de 7·10⁻³ de la taille du trou l'un de
l'autre, ils finissent au croisement (rayons presque égaux, fins légèrement
obliques), ou les deux congés reculent d'un quart de leur rayon (contacts
obliques, rayons différents : courbures incompatibles au point commun) ; la fin
d'un congé entre deux courbes de raccord de moins d'un dixième du trou recule
jusqu'à ce seuil, et un retrait demandé qui en laisserait de trop courtes est
refusé ; la face s'appuie sur une sphère guide quand ses côtés sont près d'une
sphère.*

Pour les autres cas, chaque congé s'arrête à une station : son **retrait**.

- **Croisements.** Sur chaque face Fᵢ entre deux congés i et i + 1, leurs
  contacts se croisent en Xᵢ (Gauss-Newton sur les deux pcourbes, refus s'ils
  ne se croisent pas dans les congés prolongés).
- **Retrait automatique.** Le congé i s'arrête à la plus éloignée de v des deux
  stations dont un contact passe par un croisement (Xᵢ₋₁ ou Xᵢ). Son arc
  d'extrémité (la section à cette station : arc de cercle à rayon constant,
  arc du cercle caractéristique à rayon variable) est alors hors des boules
  des autres congés.
- **Retrait de l'utilisateur.** `FilletOptions::setbacks`, une distance par
  couple (sommet, arête), mesurée le long de l'arête depuis le sommet ; elle
  ne peut que reculer le retrait automatique (refus sinon). CATIA et
  SolidWorks appellent ce coin à grands retraits un coin « valise »
  (*suitcase corner*).
- **Courbes de raccord.** Sur Fᵢ, le contact du congé i finit en Eᵢ, celui du
  congé i + 1 en E′ᵢ. Si Eᵢ = E′ᵢ (= Xᵢ), le bord du raccord a là un sommet.
  Sinon, une courbe de Fᵢ va de Eᵢ à E′ᵢ, tangente aux deux contacts : un
  segment de Hermite cubique dans les paramètres de la face, tangentes celles
  des pcourbes des contacts, longueurs de tangente égales à la corde. Refus si
  elle sort de la face ou recoupe un contact.

Le raccord a donc m côtés, n ≤ m ≤ 2n : les n arcs d'extrémité et les courbes
de raccord non réduites à un point. À la sphère de même rayon sans retrait,
les croisements sont les pieds de c et les arcs ceux de la sphère.

### 4.3 Rubans

Chaque côté porte un **ruban** : sa courbe P(s) et une dérivée transverse T(s),
tangente à la surface voisine et dirigée vers l'intérieur du raccord.

- **Arc d'extrémité** d'un congé : T est la dérivée de la surface du congé le
  long de son arête, vers v, à la station du retrait.
- **Courbe de raccord** sur une face : T est dans le plan tangent de la face,
  perpendiculaire à la courbe, vers le côté ôté de la face.

Aux sommets du polygone, les deux rubans voisins ont le même plan tangent : à
un croisement Xᵢ, les deux arcs sont tangents à Fᵢ ; au bout d'une courbe de
raccord, l'arc et la courbe sont tangents à la même face. Les rubans sont donc
compatibles en G1 aux coins. Leurs torsions ne le sont pas en général, ce que
la surface transfinie absorbe (fonctions de mélange rationnelles, comme le
carreau de Gregory).

Longueur de T : proportionnelle à la largeur du raccord au droit du côté (la
distance au côté opposé), pour que les rubans ne débordent pas les uns sur les
autres.

### 4.4 Surface transfinie

*Révisé par la PR 2 ([note](brep6_pr02_nsided.md), § 4.1) : la surface de Kato,
mesurée, s'enfonce d'un tiers de la profondeur d'un trou pentagonal dans un
paraboloïde ; elle n'est plus utilisée. La forme vient de l'énergie de flexion
(§ 4.5).*

Le domaine est le **polygone régulier** à m côtés, le paramètre de chaque ruban
porté linéairement sur son côté. Les dérivées transverses des rubans sont
rendues **compatibles aux coins** : le long du côté i, a(s) P′(s) + b(s) X(s), a
et b linéaires, fixés aux coins par les tangentes des deux courbes qui s'y
rencontrent.

### 4.5 Approximation B-spline

Une B-spline bicubique sur un rectangle qui contient le polygone du domaine,
nœuds uniformes N × N, d'**énergie de flexion minimale** (intégrée sur le
polygone par une règle de Gauss sur des triangles ; faible hors du polygone),
par moindres carrés pondérés :

- **position** aux points d'échantillonnage des côtés : S = Pᵢ(s) ;
- **tangence** aux mêmes points : nᵢ(s) · ∂S/∂ξ = 0 et nᵢ(s) · ∂S/∂η = 0 ;
- **dérivées transverses** compatibles, tenues au milieu des côtés, relâchées
  aux coins (où les torsions des deux côtés ne s'accordent pas) ;
- N doublé à partir de 4 jusqu'à ce que l'écart de position aux côtés soit sous
  `approx_tol`, l'angle des normales sous `angular_tol`, et que la forme ne
  bouge plus (2·10⁻³ de la taille du trou) ; 64 au plus, refus au-delà.

La face est la B-spline découpée par le polygone : ses pcourbes sont les côtés,
des segments. Les arêtes communes avec les congés et les faces gardent leur
courbe exacte (arc sur le congé, courbe de raccord sur la face), et l'écart à
la B-spline le long d'elles est absorbé dans leur tolérance (`settle`), comme
les arcs des coins de la PR 7 du palier 4.

### 4.6 Contrôles et refus

- côtés de longueur nulle hors des croisements, polygone de bords qui se
  recoupe (refus) ;
- après approximation : position sous `approx_tol` et normales sous
  `angular_tol` le long de tous les côtés, aucune auto-intersection grossière
  (normales de la B-spline à moins de 90° de celles de la surface guide sur
  une grille) ; refus sinon.

---

## 5. Topologie et Python

- **Onglet** : une nouvelle sorte d'extrémité (`EndKind::Mitre`), l'arête
  d'onglet portant ses pcourbes sur les deux traitements ; les arêtes vives
  raccourcies ; les faces voisines coupées comme aux extrémités libres.
- **Raccord** : `EndKind::VertexBlend`. La face du raccord, ses arêtes (arcs
  d'extrémité partagés avec les congés, courbes de raccord partagées avec les
  faces), orientée selon la face du congé voisin le long de leur arc commun,
  comme les fuseaux du palier 5. Dans le fil de chaque face Fᵢ, la suite de
  co-arêtes autour de v est remplacée par les contacts et la courbe de
  raccord.
- **Sphère à n arcs** : le chemin du coin de la PR 7, n arcs au lieu de trois.

Python : `FilletOptions.setbacks` (liste de `Setback(vertex, edge,
distance)`), et dans le noyau `gbs.fill_n_sided(ribbons, approx_tol,
angular_tol)`, qui rend la B-spline et le polygone de son domaine.

---

## 6. Validation

| Contrôle | Attendu |
|---|---|
| validité | `check` sans défaut, coques fermées, sur tous les cas |
| onglet | volumes de référence par intégration indépendante (un point est dans le résultat s'il est dans chacun des solides traités seuls, à un coin convexe) à 10⁻⁸ ; chanfreins : volume exact du polyèdre |
| sphère à n arcs | volume de référence par intégration de la réunion des boules intérieures (comme au palier 5) ; sommet de pyramide régulière |
| raccord | G1 mesuré le long de tous les côtés (angle des normales sous `angular_tol`), position sous `approx_tol` ; courbures principales bornées (pas d'ondulation : au plus quelques fois 1 / plus petit rayon) |
| continuité | rayons (r, r, r(1 + ε)) au coin d'une boîte, sans retrait : le volume tend vers celui de la sphère quand ε tend vers 0 |
| symétrie | coin de boîte symétrique (mêmes rayons, mêmes retraits), permutations des rayons : volumes égaux à 10⁻⁸ |
| OCCT | `BRepFilletAPI_MakeFillet` (rayons différents aux arêtes d'un coin), `BRepFilletAPI_MakeChamfer` ; résultats de gbs relus par OCCT valides, mêmes nombres de faces quand OCCT construit le même découpage, volumes rapportés (sa surface de coin n'est pas la nôtre) ; onglets : volumes à 10⁻⁶ |
| galerie | coin de boîte à trois rayons, coin valise (retraits), sommet de pyramide (sphère à quatre arcs), deux arêtes d'une boîte (onglet), coin chanfreiné, coin d'équerre (mêlé), congé à rayon décroissant vers un coin |

---

## 7. Plan des PR

| # | PR | Contenu | Estimation |
|---|---|---|---|
| 1 | onglets | § 2, § 3, § 4.1 : classement des sommets, onglets (congés, chanfreins), extrémités libres à plus de trois arêtes, sphère à n arcs ; volumes de référence, comparaison OCCT ; ce document, et les statuts des notes du palier 5 ([note](brep6_pr01_mitres.md)) | 16 h |
| 2 | surface à n côtés | § 4.3 à § 4.6 dans le noyau : rubans, dérivées transverses compatibles, B-spline d'énergie de flexion minimale sous contraintes G1 (`gbs/nsided.h`), Python `gbs.fill_n_sided` ; tests géométriques (G1, position, quadriques retrouvées, sphère, rotation) ([note](brep6_pr02_nsided.md)) | 18 h |
| 3 | raccords de sommet | § 4.2, § 5, § 6 : retraits automatiques et de l'utilisateur, courbes de raccord, faces des raccords, lois variables, coins mêlés, Python, galerie, comparaison OCCT, bilan ([note](brep6_pr03_vertex_blends.md)) | 18 h |
| | | **Total** | **≈ 52 h** |

---

## 8. Questions à valider

1. **Trois traitements choisis par la configuration** (§ 2) : onglet si une
   arête du sommet est vive ou s'il y a un chanfrein, sphère quand elle est
   exacte, raccord sinon ; pas de mode imposé par l'utilisateur à la première
   version. *Recommandé* : c'est le comportement par défaut de CATIA et
   SolidWorks.
2. **Onglet par découpe locale** (intersection des surfaces des traitements),
   plutôt que par booléens. *Recommandé* : les surfaces des traitements se
   coupent transversalement ; les booléens auraient à suivre les contacts,
   tangents tout le long.
3. **Chanfreins en onglet**, sans facette au coin. *Recommandé.*
4. **Raccord en une seule face**, B-spline bicubique découpée par un polygone,
   plutôt que n carreaux autour d'un point central. *Recommandé* : une face par
   coin, comme CATIA et SolidWorks, pas de G1 à reconstruire entre carreaux.
5. **Surface transfinie à rubans** (Kato, ou *Generalized Coons patch*) comme
   forme guide, approchée sous contraintes G1 linéaires et énergie de flexion.
   *Recommandé.* Révisé par la PR 2 (à valider, point 1 de sa note) : pas de
   guide, l'énergie de flexion seule sous les contraintes, la guide de Kato
   s'enfonçant trop.
6. **Retraits automatiques minimaux** (aux croisements des contacts), et retrait
   par arête au choix de l'utilisateur, qui ne peut que l'allonger.
   *Recommandé.* Complété par la PR 3 (à valider, points 3 et 4 de sa note) :
   fins aux croisements ou recul d'un quart du rayon là où deux contacts
   finiraient presque au même point.
7. **Courbes de raccord sur les faces** : Hermite cubique dans les paramètres de
   la face, tangente aux deux contacts. *Recommandé* (une courbe exacte sur la
   face, donc une arête exacte côté face).
8. **Sphère gardée quand elle est exacte**, généralisée à n arcs. *Recommandé* :
   pas de surface approchée là où une surface exacte existe.
9. **G1 seulement**, sous `angular_tol`, l'écart des bords absorbé dans la
   tolérance des arêtes. *Recommandé* ; G2 plus tard.
10. **Coins mêlés** (convexes et concaves) traités, avec refus documentés là où
    les contacts ne se croisent pas. *Recommandé.*
11. **Surface à n côtés dans le noyau gbs** (`gbs/nsided.h`, Python), pas
    seulement dans le BREP. *Recommandé* : utile pour combler un trou hors des
    congés.
12. **Validation du raccord sans volume exact** : G1, courbures, continuité vers
    la sphère, symétrie ; volumes d'OCCT rapportés seulement. *Recommandé* : la
    surface d'un raccord n'est définie que par sa construction.
13. **Trois PR** (onglets ; surface à n côtés ; raccords). *Recommandé.*

### Sources

- T. Várady, A. Rockwood, « Geometric construction for setback vertex
  blending », *Computer-Aided Design* 29(6), 1997 : retraits, courbes de
  raccord, raccord à 2n côtés.
- K. Kato, « Generation of N-sided surface patches with holes »,
  *Computer-Aided Design* 23(10), 1991 : surface transfinie à n côtés.
- P. Salvi, T. Várady, A. Rockwood, « Ribbon-based transfinite surfaces »,
  *Computer Aided Geometric Design* 31(9), 2014 : rubans, paramétrisation sur
  le polygone, *Generalized Coons patch*.
- P. Charrot, J. A. Gregory, « A pentagonal surface patch for computer aided
  geometric design », *Computer Aided Geometric Design* 1(1), 1984 : mélange
  rationnel aux coins, torsions incompatibles.
- Vida, Martin et Várady, « A survey of blending methods that use parametric
  surfaces », *Computer-Aided Design* 26(5), 1994, § 5 : raccords de sommet.
- Le comportement de CATIA (« Edge Fillet », *Blend corner*, *Setback
  distance*) et de SolidWorks (*Setback fillet*). OCCT
  (`BRepFilletAPI_MakeFillet`, `BRepFilletAPI_MakeChamfer`) en boîte noire.

---

## 9. Bilan du palier

Les trois PR ont été menées. Les configurations courantes d'un sommet que les
paliers 4 et 5 refusaient sont maintenant construites :

| Configuration | Traitement | PR |
|---|---|---|
| deux traitements à un sommet de trois arêtes, la troisième vive | onglet | 1 |
| trois chanfreins au coin | onglet en un point | 1 |
| extrémité libre à plus de trois arêtes | coupée par n − 2 faces | 1 |
| n congés de même rayon, lignes des centres concourantes | sphère à n arcs | 1 |
| rayons différents, plus de trois arêtes, rayon décroissant vers le coin, convexe et concave | raccord de sommet | 3 |
| retraits de l'utilisateur | coin valise | 3 |

La surface à n côtés (`gbs/nsided.h`, PR 2) sert aussi hors du BREP. La PR 3 l'a
rendue invariante d'échelle et lui a ajouté une sphère guide.

**Écarts au document.** Les écarts ont été validés au fil des notes, ou sont à
valider dans la dernière :

- PR 2 : la guide de Kato a été remplacée par l'énergie de flexion seule.
- PR 3 : fins aux croisements, recul d'un quart du rayon, recul des fins entre
  deux courbes de raccord courtes (§ 4.2) ; sphère guide pour la continuité
  vers la sphère ; jusqu'à 128 carreaux.

**Ce qui reste**, pour la feuille de route :

- Les raccords G2, sans la couche limite de courbure au bord.
- Un mode de coin imposé par l'utilisateur.
- Les chanfreins à distance variable et les sections coniques.
- Un domaine adapté aux coins (Salvi et Várady), qui réduirait le coût des
  trous de 6 côtés et plus.
- La facette d'OCCT au coin de trois chanfreins, à comparer à CATIA et à
  SolidWorks.
