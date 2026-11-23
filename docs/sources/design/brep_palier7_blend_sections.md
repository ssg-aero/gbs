# Noyau BREP natif — conception du palier 7 : chanfreins variables et sections de congé

| | |
|---|---|
| Statut | Questions du § 8 validées (23 nov. 2026), avec la règle C pour la section G2 (§ 4.4) ; chaque PR valide sa note et a sa note d'illustration |
| Périmètre | Chanfreins dont les distances suivent une loi le long d'une arête ou d'une chaîne, et chanfreins distance-angle. Congés dont la section n'est pas un arc de cercle : section **conique** (paramètre ρ), section **à corde constante**, section **à courbure continue** (G2 avec les faces). Pour chacun : géométrie, topologie (extrémités, chaînes, onglets, coins), Python, galerie, comparaison à OpenCascade. |
| Dépend de | Palier 4 : géométrie et topologie des congés et chanfreins (`fillet.h`, `blend.h`). Palier 5 : lois (`radius_law.h`), morceaux, fuseaux et plis. Palier 6 : onglets, raccords de sommet (`nsided.h`). |
| Hors périmètre | Congés de faces (*face blend*). Raccords de sommet G2. Mode de coin imposé par l'utilisateur. Section définie par une courbe quelconque de l'utilisateur. |

Ce document décrit *quoi* construire et *pourquoi*, comme les documents des
paliers 4 à 6. Chaque PR a sa note d'architecture `brep7_prNN_*.md` et sa note
d'illustration. C'est la suite de la feuille de route laissée par les paliers 5
(§ 1.2) et 6 (§ 9) : chanfreins à distance variable, sections coniques.

**Provenance.** Les algorithmes viennent de la littérature et du comportement de
CATIA et SolidWorks, la cible de qualité ; chaque note cite ses sources.
OpenCascade ne sert que de comparaison en boîte noire dans les tests.

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Loi d'un chanfrein | **la loi du palier 5**, généralisée : `BlendLaw` (points de contrôle sur les arêtes, valeur, variation linéaire ou lisse, abscisse curviligne le long de la chaîne) ; `RadiusLaw` en devient un alias, l'API du palier 5 inchangée | une seconde structure `DistanceLaw` copiée de la première |
| Chanfrein variable | **deux lois**, d₁(s) sur la face donnée et d₂(s) sur l'autre (une seule : symétrique) ; contacts à ces distances dans le plan normal à l'arête, comme au palier 4 ; surface **réglée** entre les contacts | surface réglée entre l'arête décalée et les contacts (pas de sens pour un chanfrein) |
| Chanfrein distance-angle | **d₁(s) et un angle θ constant** mesuré depuis la face de d₁ : d₂ en chaque station par Gauss-Newton (la règle part du contact 1 sous l'angle θ et rencontre l'autre face) | d₂ = d₁ tan θ (juste seulement entre deux plans à angle droit) |
| Point anguleux d'une loi linéaire | **pli** du chanfrein (deux surfaces réglées, une arête vive entre elles), la topologie des plis du palier 5 | lisser la loi |
| Section d'un congé | une **forme de section** (`BlendSection`) dans le plan de l'arc du palier 4 ou 5, entre les **mêmes contacts** (ceux de la boule de rayon r ou r(s)) : arc de cercle (par défaut), conique ρ, à courbure continue ; la corde constante fixe les contacts autrement | une opération par forme de section (code dupliqué) |
| Section conique | l'arc de **conique** tangent aux deux faces aux contacts : B-spline rationnelle quadratique de pôles P₁, S, P₂ (S le point où se coupent les tangentes aux faces dans le plan de la section) et de poids w = ρ / (1 − ρ) en S ; ρ dans ]0, 1[, ρ = 0,5 une parabole ; l'arc de cercle pour w = cos(θ/2), θ l'angle balayé (ρ ≈ 0,414 à une arête à angle droit) ; ρ constant ou en loi | une conique donnée par son excentricité (moins parlant que ρ, le paramètre de CATIA) |
| Corde constante | les contacts à la **distance c l'un de l'autre** au lieu d'être sur une boule de rayon donné : la station résout \|P₁ − P₂\| = c (au lieu de \|C − P₁\| = \|C − P₂\| = r) ; la section un arc de cercle tangent aux faces, le rayon qui en découle variable | le rayon équivalent calculé une fois (faux dès que l'angle des faces varie) |
| Courbure continue (G2) | la section **de degré 4** (5 pôles) tangente aux faces et de même **courbure normale** que chacune au contact ; sa forme la **plus lisse** (∫ κ′² minimale) ; ses contacts placés pour qu'elle ait la **flèche de l'arc de rayon r** (à angle droit, contacts à 1,30 r, pic de courbure 1,23/r, sans saut) | les contacts de la boule de rayon r avec la même flèche (pic de 1,38/r collé aux contacts : pire reflet que l'arc), ou la forme la plus lisse entre ces contacts (pic de 1,60/r, congé plus plein) |
| Surface | le chemin des congés approchés du palier 4 : sections exactes aux stations, peau de Hermite, raffinement sous `approx_tol` et `angular_tol` ; **exacte** quand c'est possible : conique ou G2 le long d'une arête droite entre deux plans = extrusion (degré (2 ou 4, 1)) ; conique autour d'une ligne des centres circulaire = révolution | toujours approchée |
| Topologie | **inchangée** : extrémités libres, chaînes, onglets, raccords de sommet du palier 6 ; un **coin sphérique** seulement pour des arcs de cercle ; des congés coniques ou G2 au même coin font un **raccord de sommet** | coin sphérique déformé (pas une surface exacte) |
| Continuité au coin | un raccord de sommet entre congés G2 reste **G1** (palier 6) ; noté comme limite | G2 des raccords de sommet (hors périmètre) |
| API | `chamfer(m, s, edges, face, BlendLaw d1, BlendLaw d2)`, `chamfer_angle(m, s, edges, face, d1, angle)` ; `FilletOptions::section` (`BlendSection` : forme, ρ, profondeur) ; `fillet_chordal(m, s, edges, c)` ; mêmes formes dans `fillet_geometry`, `chamfer_geometry` ; Python à l'identique | des options cachées dans des surcharges de `fillet` |

---

## 1. Objectifs et non-objectifs

### 1.1 Objectifs

- **Chanfreins variables**, comme le « Chamfer » de CATIA (longueurs variables
  d'un bout à l'autre, *Length1/Length2* et *Length/Angle*) et le chanfrein
  de SolidWorks (*Distance-distance* asymétrique, *Angle-distance*) : les
  distances de part et d'autre de l'arête suivent une loi, ou une distance et
  un angle.
- **Congés coniques** : le « Conic » de CATIA (*Conic parameter* ρ) et le
  congé « Conic Rho » de SolidWorks : une section plus pleine (ρ > ρ_cercle) ou
  plus creuse que l'arc, entre les mêmes contacts.
- **Congés à corde constante** : le « Chordal fillet » de CATIA, la largeur du
  congé fixée quel que soit l'angle des faces.
- **Congés à courbure continue** : l'option « Curvature continuous » de
  SolidWorks et la section G2 de CATIA : la courbure ne saute plus au contact,
  les reflets restent continus.
- Toute la topologie des paliers 4 à 6 pour ces nouvelles formes : chaînes,
  boucles, extrémités libres, onglets, raccords de sommet.
- Python, galerie, comparaison à `BRepFilletAPI_MakeChamfer` (`AddDD`,
  `AddDA`) et à `BRepFilletAPI_MakeFillet` (section circulaire seulement :
  OCCT n'a ni conique ni corde constante).

### 1.2 Non-objectifs

- Congés de faces (*face blend*) : un autre palier.
- Raccords de sommet G2 : un raccord entre congés G2 reste G1 sur ses bords.
- Section quelconque dessinée par l'utilisateur.
- Section G2 à rayon variable dans la première version de la PR 3 : proposée
  avec rayon constant ou en loi comme les autres sections, sans garantie
  supplémentaire sur la courbure le long de l'arête.

---

## 2. Lois

### 2.1 `BlendLaw`

La structure du palier 5 (`RadiusPoint` : arête, fraction de sa longueur,
valeur ; variation linéaire ou lisse ; abscisse curviligne le long d'une chaîne
tangente) sert à toute grandeur d'un traitement : rayon, distance d'un
chanfrein, ρ d'une conique. Elle est renommée `BlendLaw<T>` (points
`BlendPoint<T>` avec un champ `value`), `RadiusLaw` et `RadiusPoint` restant des
alias (et `radius` un nom du champ en Python) : rien ne change pour les
appels du palier 5.

### 2.2 Prolongement et contrôles

Comme au palier 5 : la loi est prolongée avec sa pente au-delà de ses points
extrêmes (linéaire) ou constante (lisse) ; refus où elle n'est plus positive
sur le traitement prolongé. Pour ρ : refus hors de ]0, 1[.

---

## 3. Chanfreins

### 3.1 Distances variables

Le chanfrein du palier 4 a ses contacts à la distance dₖ de l'arête dans le
plan normal à l'arête, sur la face k (Gauss-Newton depuis les stations
voisines, pcourbe interpolée doublée jusqu'à tol / 10). Avec une loi, la
distance en chaque station est dₖ(s(t)) ; le reste est inchangé.

La surface est **réglée** entre les deux contacts, au même paramètre :
primitive quand les contacts le permettent (plan, cylindre, cône, comme au
palier 4), B-spline réglée sinon. Cas exact remarquable : entre deux plans, une
loi linéaire le long d'une arête droite donne deux contacts droits non
parallèles, et la surface réglée est un **paraboloïde hyperbolique** (B-spline
bilinéaire, exacte).

### 3.2 Distance-angle

La distance d₁(s) sur la face donnée ; la règle du chanfrein part du contact 1
dans le plan normal à l'arête, sous l'angle θ avec la face 1 (vers l'intérieur
du solide à une arête convexe), et rencontre l'autre face : c'est le contact 2
(Gauss-Newton : une inconnue de plus par station, le paramètre le long de la
règle). Entre deux plans à angle droit, d₂ = d₁ tan θ ; ailleurs, d₂ varie avec
l'angle des faces même si d₁ est constant. Refus si la règle ne rencontre pas
l'autre face dans le chanfrein prolongé (θ trop grand pour l'angle des faces).

### 3.3 Plis et morceaux

Un point anguleux d'une loi linéaire (sa pente saute) donne un **pli** : le
chanfrein en deux morceaux réglés, l'arête entre eux la règle à la station du
point. C'est la topologie des plis du palier 5 (`BlendCrease`), sans fuseau
(pas de boule).

### 3.4 Topologie

Celle des chanfreins des paliers 4 et 6 : extrémités libres, chaînes tangentes
(la loi continue d'une arête à l'autre : les règles se rejoignent), onglets de
deux chanfreins, trois chanfreins au coin (point commun), chanfrein et congé
mêlés. Aucun cas nouveau : les découpes travaillent sur les surfaces, quelles
que soient les distances.

---

## 4. Sections de congé

### 4.1 Cadre commun

Au palier 4 (rayon constant) comme au palier 5 (loi), chaque station du congé
donne deux contacts P₁, P₂ et le plan de la section (normal à la ligne des
centres, ou celui du cercle caractéristique), et la section est l'arc de
cercle entre eux. La forme de section (`BlendSection`) remplace cet arc par
une autre courbe **du même plan, entre les mêmes points, tangente aux mêmes
faces** :

- les contacts, la ligne des centres et toute la topologie restent ceux du
  congé circulaire de même rayon (les extrémités libres, onglets et raccords
  coupent la surface, quelle qu'elle soit) ;
- seule la surface change : sections aux stations, peau de Hermite, raffinement
  sous `approx_tol` / `angular_tol`, comme au palier 4.

Le point S où se coupent, dans le plan de la section, les tangentes aux deux
faces aux contacts sert à toutes les formes (défini dès que les faces ne sont
pas tangentes, ce que le congé refuse déjà).

### 4.2 Section conique

La section est la conique tangente aux faces en P₁ et P₂ : la B-spline
rationnelle quadratique de pôles P₁, S, P₂, de poids 1, w, 1, avec
w = ρ / (1 − ρ). ρ est le paramètre de CATIA : le point d'épaule de la conique
partage le segment du milieu de [P₁P₂] à S dans le rapport ρ. ρ < 0,5 : ellipse
(plus creuse), ρ = 0,5 : parabole, ρ > 0,5 : hyperbole (plus pleine). L'arc de
cercle est le cas w = cos(θ/2), θ l'angle balayé par l'arc (90° à une arête à
angle droit : ρ = cos 45° / (1 + cos 45°) ≈ 0,414) : la section conique contient
le congé circulaire, et le test de non-régression le vérifie.

ρ constant, ou une `BlendLaw` le long de la chaîne (avec un rayon constant ou
en loi).

**Exact** : le long d'une arête droite entre deux plans, l'extrusion de la
conique (B-spline rationnelle de degré (2, 1)) ; ailleurs approchée comme les
congés du palier 4.

### 4.3 Corde constante

Le congé de corde c : à chaque station, les contacts sont les points des deux
faces dans le plan normal à l'arête où un cercle tangent aux deux faces a une
corde c entre ses contacts. La station du palier 5 (Gauss-Newton sur les
paramètres des deux contacts, rayon donné) devient : même système, le rayon
inconnu, l'équation |P₁ − P₂| = c en plus. Le rayon r(t) qui en sort est celui
d'un congé à loi ; la surface est l'enveloppe des boules de ce rayon
(cercle caractéristique), exactement comme au palier 5 — d'où la même
machinerie. Entre deux plans d'angle constant, le rayon est constant et le
congé est un cylindre : cas de contrôle.

### 4.4 Courbure continue (G2)

La section est une courbe de Bézier de degré 4, de pôles Q₀ = P₁, Q₁, Q₂, Q₃,
Q₄ = P₂ :

- **tangence** : Q₁ sur la tangente à la face 1 en P₁ (vers S), Q₃ sur celle
  de la face 2 ;
- **courbure** : la courbure de la section en P₁ égale la courbure normale de
  la face 1 dans le plan de la section (0 pour un plan), idem en P₂. La
  courbure d'une Bézier à son extrémité vaut κ = (n − 1)/n · h/a², où h est la
  hauteur de Q₂ sur la droite Q₀Q₁ et a = |Q₀Q₁|. Elle fixe donc h. Entre deux
  plans, Q₂ est en S.

**Les degrés de liberté restants.** Les longueurs |P₁Q₁| et |P₂Q₃|, prises
égales par symétrie (λ = |P₁Q₁| / |P₁S|), et la place des contacts. La note
d'illustration du § 8 compare trois règles, faces à angle droit :

| Règle | Contacts | Pic de courbure | Défaut |
|---|---|---|---|
| A : contacts de la boule de rayon r, même flèche que l'arc | r | 1,38/r, collé aux contacts | toute la transition dans 3 % de la section : reflets pires qu'avec l'arc |
| B : contacts de la boule, forme la plus lisse | r | 1,60/r au milieu | congé plus plein (flèche 1,33 fois celle de l'arc) |
| **C : forme la plus lisse, même flèche que l'arc (retenue)** | **1,30 r** | **1,23/r au milieu** | contacts plus écartés que ceux du congé de rayon r |

La forme la plus lisse minimise ∫ κ′(s)² ds le long de la section. C'est un
problème à une inconnue, λ ≈ 0,35 à angle droit, résolu par Brent à chaque
station. La règle C place ensuite les contacts pour que la flèche, c'est-à-dire
la distance du milieu de la section au coin S, soit celle de l'arc de rayon r.
Le congé G2 enlève alors autant de matière au milieu que le congé circulaire de
même rayon : r garde son sens pour l'utilisateur. Les congés « à courbure continue » de SolidWorks semblent aussi plus larges que
le congé circulaire de même rayon ; c'est à vérifier dans SolidWorks, je n'ai
pas pu le confirmer. La station résout ainsi les contacts à la distance
k(φ)·r, où k ne dépend que de l'angle des faces dans le plan de la section.
Pour une face courbe, la courbure normale entre dans les équations.

**Exact** entre deux plans le long d'une arête droite : extrusion de degré
(4, 1). Ailleurs approchée. La courbure de la surface en travers, au contact,
est vérifiée égale à celle de la face (à 10⁻³ relatif près, mesurée sur la
B-spline approchée).

### 4.5 Coins et raccords

Une sphère n'est le coin de n congés que pour des sections circulaires de même
rayon. Dès qu'un congé au coin est conique, à corde constante (rayon variable)
ou G2, le coin est un **raccord de sommet** (palier 6) : ses rubans sont les
sections aux retraits (quelle que soit leur forme) et les courbes de raccord.
Le raccord reste G1 sur ses bords (non-objectif G2). Les onglets et les
extrémités libres ne changent pas.

---

## 5. Topologie et Python

Aucune sorte d'extrémité nouvelle. Les plis des chanfreins variables suivent
les plis du palier 5. `fillet`, `chamfer`, `fillet_geometry`,
`chamfer_geometry` prennent les nouvelles formes ; `fillet_chordal` en plus.

Python : `brep.BlendLaw`, `brep.BlendPoint` (`RadiusLaw`, `RadiusPoint`
gardés), `brep.BlendSection` (`circular()`, `conic(rho)`, `conic(law)`,
`curvature_continuous()`), `FilletOptions.section`,
`brep.chamfer(..., d1_law, d2_law)`, `brep.chamfer_angle`,
`brep.fillet_chordal`.

---

## 6. Validation

| Contrôle | Attendu |
|---|---|
| validité | `check` sans défaut, coques fermées, tous les cas |
| chanfrein variable, arête de boîte | volume exact ∫ d₁(s) d₂(s) / 2 ds (lois linéaire et lisse), à 10⁻⁹ ; paraboloïde hyperbolique exact pour la loi linéaire |
| distance-angle | entre deux plans à angle droit, d₂ = d₁ tan θ, volume exact ; entre un plan et un cylindre, contact 2 vérifié sur la règle sous l'angle θ (10⁻⁹) |
| conique | ρ = ρ_cercle : la surface du congé circulaire (même volume à 10⁻¹², même géométrie à `approx_tol`) ; arête de boîte : aire de la section de conique en forme fermée, volume à 10⁻⁹ ; ρ croissant : volume enlevé décroissant |
| corde constante | entre deux plans : le cylindre de rayon c / (2 sin(φ/2)) ; entre un plan et un cylindre : \|P₁ − P₂\| = c à tol en chaque station mesurée |
| G2 | courbure en travers au contact égale à celle de la face (plan : 0) ; arête de boîte : contacts à 1,30 r, flèche de l'arc, aire de la section quartique en forme fermée, volume à 10⁻⁹ ; bandes zébrées continues **et lisses** au contact (illustration) |
| coins | congés coniques, G2, à corde constante au coin d'une boîte : raccord de sommet, G1 sous les tolérances ; trois chanfreins variables au coin : point commun dans chacun |
| chaînes, plis | chaîne de deux arêtes tangentes, loi continue : règles et sections communes ; loi linéaire à un point intermédiaire : un pli, deux morceaux |
| OCCT | `AddDD` (asymétrique constant) et `AddDA` : mêmes volumes à 10⁻⁹, mêmes faces ; congé circulaire inchangé ; conique, corde constante, G2 : volumes de gbs seulement (OCCT ne les a pas) |
| galerie | boîte à chanfrein variable, chanfrein distance-angle sur un cylindre, congés coniques ρ = 0,3 / 0,5 / 0,7 côte à côte, congé à corde constante entre un plan et un cylindre, congé G2 en bandes zébrées face au congé circulaire, coin de trois congés coniques |

---

## 7. Plan des PR

| # | PR | Contenu | Estimation |
|---|---|---|---|
| 1 | chanfreins variables | § 2, § 3 : `BlendLaw` (alias `RadiusLaw`), chanfreins à deux lois, distance-angle, plis, topologie, Python, galerie, OCCT (`AddDD`, `AddDA`) ; ce document, les statuts des notes du palier 6 ([note](brep7_pr01_variable_chamfers.md)) | 14 h |
| 2 | sections conique et à corde constante | § 4.1 à § 4.3, § 4.5 : `BlendSection`, conique (ρ constant ou en loi, exacte sur l'arête droite entre plans), corde constante, coins en raccord de sommet, Python, galerie ; les statuts des notes du palier 6 et de la PR 1 ([note](brep7_pr02_conic_chordal.md)) | 16 h |
| 3 | section G2 et bilan | § 4.4 : section quartique G2, exacte entre plans, vérification de la courbure, coins, Python, galerie, bilan du palier | 14 h |
| | | **Total** | **≈ 44 h** |

---

## 8. Questions à valider

1. **Périmètre** : chanfreins variables + trois sections de congé (conique,
   corde constante, G2) dans un même palier, les congés de faces renvoyés à un
   palier suivant. *Recommandé* : même machinerie (stations, sections, peau),
   la topologie inchangée.
2. **`BlendLaw` commune**, `RadiusLaw` gardée comme alias. *Recommandé.*
3. **Chanfrein à deux lois** (d₁ sur la face donnée, d₂ sur l'autre) et
   **distance-angle** à angle constant. *Recommandé* (les deux modes de CATIA et
   SolidWorks).
4. **Pli** au point anguleux d'une loi linéaire de chanfrein (comme au palier 5).
   *Recommandé.*
5. **Conique en ρ** (le paramètre de CATIA), constant ou en loi, entre les
   contacts du congé circulaire de même rayon. *Recommandé* : la conique
   contient l'arc (non-régression), et le rayon garde son sens (la largeur).
6. **Corde constante par la station du palier 5** (rayon inconnu, corde
   imposée), surface enveloppe des boules. *Recommandé.*
7. **Section G2 quartique**, symétrique, de forme la plus lisse (∫ κ′²
   minimale), ses contacts écartés pour garder la flèche de l'arc de rayon r
   (règle C du § 4.4 : 1,30 r à angle droit, pic de 1,23/r). *Recommandé* : la
   même matière enlevée au milieu que le congé circulaire, la courbure sans
   saut ni pic. Écartées : les contacts de la boule avec la même flèche (pic
   de 1,38/r collé aux contacts) et la forme la plus lisse entre ces contacts
   (congé plus plein). L'illustration du § 8 les compare.
8. **Surfaces exactes** sur l'arête droite entre deux plans (extrusions
   rationnelles de degré 2 ou 4), approchées ailleurs. *Recommandé.*
9. **Coins non circulaires en raccord de sommet G1** (pas de sphère
   déformée, pas de G2 au coin). *Recommandé.*
10. **Validation** : volumes exacts en forme fermée sur les arêtes de boîte,
    non-régression conique = cercle, courbure au contact pour G2, OCCT pour les
    chanfreins seulement. *Recommandé.*
11. **Trois PR** (chanfreins ; conique et corde ; G2 et bilan). *Recommandé.*

### Sources

- Vida, Martin et Várady, « A survey of blending methods that use parametric
  surfaces », *Computer-Aided Design* 26(5), 1994 : sections de congé, coniques,
  corde constante.
- L. Piegl, W. Tiller, *The NURBS Book*, 2ᵉ éd., § 7.5 : coniques rationnelles
  quadratiques, poids et paramètre d'épaule.
- G. Farin, *Curves and Surfaces for CAGD*, 5ᵉ éd., § 12 (coniques) et § 11
  (courbure aux extrémités d'une Bézier) : la section G2.
- J. Hoschek, D. Lasser, *Fundamentals of Computer Aided Geometric Design*,
  § 14 (*blending*) : congés à courbure continue.
- Le comportement de CATIA (« Chamfer » *Length1/Length2*, *Length/Angle* ;
  « Edge Fillet » *Conic parameter*, *Chordal fillet*) et de SolidWorks
  (chanfrein *Distance-distance*, *Angle-distance* ; congé *Conic Rho*,
  *Curvature continuous*). OCCT (`BRepFilletAPI_MakeChamfer`,
  `BRepFilletAPI_MakeFillet`) en boîte noire.
