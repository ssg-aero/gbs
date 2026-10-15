# Noyau BREP natif — conception du palier 2

| | |
|---|---|
| Statut | **Réalisé** (PR 147 à 157, octobre 2026) ; bilan au § 12. Les recommandations des questions sont appliquées, la question 1 (booléens) comprise |
| Périmètre | Intersections géométriques (courbes 2D, courbe-surface, surface-surface), classification, section, découpe de faces, constructions par extrusion et révolution, faces décalées ; **booléens** en fin de palier (question 1) |
| Dépend de | palier 1 ([brep_core.md](brep_core.md)), primitives analytiques ([step_p1_primitives.md](step_p1_primitives.md)), lecture et écriture STEP (pour la validation par OCCT) |
| Hors périmètre | congés et chanfreins, healing général, offsets de shells et coques épaisses, balayage le long d'un chemin (question 9), maillage |

Ce document décrit *quoi* construire et *pourquoi*, comme
[brep_core.md](brep_core.md) pour le palier 1. Les pseudo-déclarations C++
sont illustratives.

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Booléens | **Fusion, différence, intersection** en fin de palier, pour les configurations **génériques** et les **faces coplanaires** ; les autres cas dégénérés rendent une erreur explicite. Revient sur le § 1.3 du palier 1 (question 1) | booléens toujours délégués à OCCT ; booléens complets d'emblée |
| Emplacement | intersections géométriques dans **`gbs/intersect/`**, sans dépendance au BREP et utilisables seules ; opérations topologiques dans `gbs-brep/` | tout dans `gbs-brep/` |
| Solutions | **toutes** les solutions dans un domaine borné : **subdivision** par boîtes englobantes puis **Newton** ; **formules fermées** pour les paires de primitives | extrema locaux existants (`extrema.h`, nlopt, une seule solution) |
| Courbes d'intersection | **exactes** quand c'est possible (droite, cercle, ellipse, pcurves exactes quand le paramétrage le permet) ; sinon **cheminement** (*marching*), puis B-spline 3D et deux pcurves, avec un écart **mesuré** porté par la tolérance de l'arête | toujours approchées ; courbe implicite sans représentation |
| Cas dégénérés | tangences, surfaces confondues, sommet sur une face : **reconnus** à la tolérance d'opération ; traités pour les plans et les primitives coaxiales, sinon **erreur** `BuildErrc::NonGeneric` désignant les entités | résultat approximatif silencieux |
| Classification | point/face en (u, v) sur les pcurves ; point/solide par **lancer de rayon**, avec une autre direction si le rayon touche une arête ou un sommet ou s'il est tangent | nombre d'enroulement 3D par intégrale |
| Découpe | **arrangement 2D** dans l'espace paramétrique de chaque face (seams et périodicité gérés) ; arêtes de bord découpées par `split_edge`, **partagées** avec la face voisine pour que le shell reste fermé | découpe en 3D |
| Tolérance | une **tolérance d'opération** (défaut `brep_default_tolerance`) ; les tolérances des entités créées **remontent** à l'écart mesuré (règle sommet ≥ arête ≥ face du palier 1) | tolérance fixe |
| Dépendances | **aucune nouvelle** ; Eigen pour les petits systèmes (déjà utilisé) ; pas de nlopt dans le nouveau code | |
| Validation | formules fermées ; **comparaison à OCCT** (`BRepAlgoAPI_Section`, `Common`, `Fuse`, `Cut` : volumes, aires, longueurs, nombres de faces) ; résultats écrits en STEP et relus **valides** par OCCT | tests internes seulement |

---

## 1. Objectifs et non-objectifs

### 1.1 Objectifs

1. **Intersections géométriques** donnant toutes les solutions : courbe-courbe
   en 2D, courbe-surface, surface-surface. Elles servent au BREP et sont
   utilisables seules, en C++ comme en Python.
2. **Classification** d'un point par rapport à une face et à un solide.
3. **Section** : les arêtes d'intersection de deux formes, ou d'une forme et
   d'une surface, avec leurs pcurves sur chaque face.
4. **Découpe** d'une face, puis d'une forme, par des arêtes qui la traversent.
5. **Constructions** : extrusion et révolution d'arêtes, de wires et de faces
   en BREP, faces décalées.
6. **Booléens** (question 1) : fusion, différence et intersection de solides.

### 1.2 Non-objectifs

- Congés et chanfreins.
- Healing général : recollement avec reconstruction, jonctions en T,
  réparation de fichiers.
- Offsets de shells et coques épaisses : ils demandent d'intersecter les faces
  décalées entre elles, c'est-à-dire les booléens, plus le traitement des
  faces qui disparaissent.
- Balayage le long d'un chemin quelconque (question 9).
- Booléens sur les configurations dégénérées autres que celles du § 6.3.

### 1.3 Ce que le palier 1 avait décidé

Le § 1.3 de [brep_core.md](brep_core.md) laissait les booléens à OCCT
« définitivement ». Ce palier construit toutes leurs briques : intersections,
section, découpe et classification. Une fois ces briques en place, les
booléens génériques coûtent deux PR. La question 1 propose donc de les faire,
mais en se limitant aux cas génériques et coplanaires. Les configurations
dégénérées restantes (tangences entre surfaces courbes, surfaces courbes
confondues non coaxiales) restent à OCCT, avec une erreur claire.

---

## 2. Intersections géométriques (`gbs/intersect/`)

### 2.1 Résultats

```cpp
namespace gbs {
template <typename T, size_t dim> struct CurveCurveHit   { T t1, t2; point<T, dim> p; bool tangent; };
template <typename T>             struct CurveOverlap    { T t1a, t1b, t2a, t2b; };          // portion commune
template <typename T>             struct CurveSurfaceHit { T t, u, v; point<T, 3> p; bool tangent; };
template <typename T>             struct SurfaceSurfaceHit { T u1, v1, u2, v2; point<T, 3> p; };

template <typename T> struct IntersectionCurve {
    std::shared_ptr<Curve<T, 3>> curve; T a, b;            // courbe 3D sur [a, b]
    std::shared_ptr<Curve<T, 2>> pc1, pc2;                 // pcurves sur chaque surface, même paramètre
    T deviation;                                           // écart mesuré entre les trois représentations
    bool exact, closed;
};
template <typename T> struct SurfaceSurfaceResult {
    std::vector<IntersectionCurve<T>> curves;
    std::vector<SurfaceSurfaceHit<T>> points;              // contacts ponctuels (tangences isolées)
    bool coincident;                                       // surfaces confondues sur le domaine
};
struct IntersectOptions { T tol; T angular_tol; T max_step; /* … */ };

auto intersect(const Curve<T, 2> &c1, T a1, T b1, const Curve<T, 2> &c2, T a2, T b2, IntersectOptions)
    -> CurveCurveResult;                                   // hits + overlaps
auto intersect(const Curve<T, 3> &c, T a, T b, const Surface<T, 3> &s, Box2 domain, IntersectOptions)
    -> CurveSurfaceResult;                                 // hits + portions sur la surface
auto intersect(const Surface<T, 3> &s1, Box2 d1, const Surface<T, 3> &s2, Box2 d2, IntersectOptions)
    -> SurfaceSurfaceResult;
}
```

Chaque intersection travaille sur un **domaine borné** : la plage de l'arête,
la boîte paramétrique de la face. Les surfaces infinies (plan d'un outil de
coupe) sont bornées par la boîte de l'autre opérande.

### 2.2 Boîtes englobantes

Les boîtes sont la base de la subdivision.

| Géométrie | Boîte |
|---|---|
| B-spline et NURBS | enveloppe des pôles d'un segment de Bézier (exacte, rétrécit avec la subdivision) |
| primitives | formule fermée (cercle, morceau de cylindre, de sphère…) |
| autres (`CurveOnSurface`, révolution, offset…) | échantillonnage élargi par une borne de la dérivée sur le pas ; c'est la seule boîte non garantie, la marge est réglable |

### 2.3 Courbe-courbe en 2D

Elle sert à la découpe : croiser les pcurves des arêtes de section avec les
pcurves du bord d'une face.

- **Droites** (côtés des faces naturelles, pcurves de degré 1) : formule
  fermée.
- **Cas général** : subdivision des deux courbes tant que les boîtes se
  chevauchent, jusqu'à des segments presque droits, puis Newton à deux
  inconnues. Les solutions trop proches sont fusionnées.
- **Portions communes** : si un segment reste dans la tolérance de l'autre
  courbe sur toute sa longueur, c'est une portion commune (deux bords
  confondus), rendue comme `CurveOverlap`.

### 2.4 Courbe-surface

Elle sert au lancer de rayon et aux arêtes qui percent une face.

- **Formules fermées.** Une droite contre un plan, un cylindre, un cône ou une
  sphère donne une équation du second degré. Contre un tore, c'est une
  équation du quatrième degré, résolue puis affinée par Newton. Un cercle
  contre un plan se calcule aussi directement.
- **Cas général** : subdivision du segment de courbe et du carreau de surface,
  puis Newton à trois inconnues (t, u, v).
- **Portion de courbe sur la surface** : détectée comme à la courbe-courbe et
  rendue comme portion.

### 2.5 Surface-surface : paires de primitives

| Paire | Résultat exact |
|---|---|
| plan / plan | une droite |
| plan / cylindre | une ou deux droites (plan parallèle à l'axe), un cercle (plan perpendiculaire), une ellipse (plan oblique) |
| plan / sphère | un cercle ou un point de contact |
| plan / cône | un cercle, une ou deux droites (plan passant par le sommet), une ellipse ; parabole et hyperbole par cheminement (question 6) |
| plan / tore | deux cercles si le plan est perpendiculaire à l'axe ; deux paires de cercles s'il contient l'axe ; sinon cheminement |
| sphère, cylindre, cône, tore **coaxiaux** entre eux | des cercles |
| deux cylindres d'axes parallèles | des droites |
| autres paires | cheminement (§ 2.6) |

**Pcurves exactes** quand le paramétrage est affine le long de la courbe :
- une droite sur un plan ;
- une iso d'une primitive (cercle sur un cylindre, une sphère ou un tore
  coaxial, génératrice d'un cylindre ou d'un cône) : un segment en (u, v).

Sinon, par exemple une ellipse sur un cylindre ou un cercle sur un plan
oblique, la pcurve est **interpolée** sur le paramètre de la courbe 3D, et son
écart est mesuré.

### 2.6 Surface-surface : cas général par cheminement

1. **Points de départ.**
   - *Courbes qui touchent le bord* : on intersecte les bords de chaque
     domaine (isos) avec l'autre surface (§ 2.4). Chaque solution est un départ.
   - *Boucles fermées intérieures* : on subdivise les deux surfaces tant que
     leurs boîtes se chevauchent. On arrête dès que les **cônes de normales**
     de deux carreaux ne se recouvrent pas : aucune boucle fermée n'y est
     possible. Sur les petits carreaux restants, Newton sur la distance donne
     un point de départ.
2. **Cheminement.** La direction est `n1 × n2`. Le pas est adapté à la
   courbure et à la flèche tolérée. Chaque pas est une prédiction, puis une
   correction de Newton sur (u1, v1, u2, v2), contrainte dans le plan normal
   au pas.
3. **Arrêts.** Le cheminement s'arrête au bord d'un domaine, au retour au
   point de départ (courbe fermée) ou sur une **singularité**. Une singularité,
   c'est `|n1 × n2|` sous la tolérance angulaire : tangence ou point de
   branchement. Elle est rendue comme cas dégénéré (§ 2.7).
4. **Courbes.** Les points calculés donnent une B-spline 3D et deux pcurves
   interpolées sur le **même paramètre**. On contrôle l'écart au milieu de
   chaque intervalle, et on raffine tant qu'il dépasse la tolérance.
5. **Surfaces périodiques** : le cheminement se fait dans l'espace déroulé. Les
   pcurves restent continues, et la découpe les ramène dans le domaine (§ 5.1).

### 2.7 Cas dégénérés

| Situation | Traitement |
|---|---|
| deux plans confondus | `coincident` ; la découpe les traite en 2D dans le plan (§ 6.3) |
| primitives coaxiales confondues (même rayon) | `coincident` |
| contact ponctuel (sphère tangente à un plan) | `points` |
| tangence le long d'une courbe entre primitives (cylindre posé sur un plan) | courbe exacte marquée tangente |
| tangence ou branchement entre surfaces quelconques | courbes arrêtées au point singulier et marquées ; les booléens rendent `NonGeneric` |

---

## 3. Classification (`gbs-brep/classify.h`)

```cpp
enum class Location { In, Out, On };
auto classify(const Model<T> &m, FaceId f, point<T, 2> uv, T tol) -> Location;
auto classify(const Model<T> &m, SolidId s, point<T, 3> p, T tol) -> Location;
auto interior_point(const Model<T> &m, FaceId f) -> std::optional<point<T, 2>>;
```

- **Point dans une face** : un rayon horizontal dans l'espace (u, v), croisé
  avec les pcurves des wires (§ 2.3). C'est exact, contrairement aux polygones
  échantillonnés du palier 1, qui restent utilisés par `check`.
- **Point dans un solide** : un rayon depuis le point, croisé avec chaque face
  (§ 2.4), puis classification de chaque impact dans sa face. Les impacts
  sont comptés par parité. Si un impact est tangent ou tombe sur une arête ou
  un sommet, on reprend avec une autre direction, tirée dans une suite fixe
  pour rester **déterministe**. Un point sur une face, à la tolérance près,
  est `On`.
- **Face par rapport à un solide** : un point intérieur à la face, loin de son
  bord, puis la classification de ce point.

---

## 4. Section (`gbs-brep/section.h`)

```cpp
auto section(Model<T> &m, ShapeId a, ShapeId b, SectionOptions) -> BuildResult<SectionResult>;
auto section(Model<T> &m, ShapeId a, std::shared_ptr<Surface<T, 3>> tool, SectionOptions) -> BuildResult<SectionResult>;
// SectionResult : un compound de wires, et pour chaque arête créée les deux faces qui la portent
```

1. **Paires de faces** filtrées par leurs boîtes englobantes.
2. **Intersection des surfaces** sur les boîtes paramétriques des deux faces
   (§ 2).
3. **Restriction aux faces.** On croise chaque pcurve de section avec les
   pcurves du bord de chaque face (§ 2.3). On en tire les intervalles dont le
   milieu est dans les deux faces (§ 3). Chaque intervalle devient une
   **arête** qui porte ses deux pcurves : une co-arête par face, comme le
   prévoyait le § 9.1 du palier 1.
4. **Sommets.** Les extrémités sont fusionnées, à la tolérance d'opération, avec
   les sommets existants et entre elles. Une arête de section qui finit sur le
   bord d'une face partage ce sommet avec la face voisine.
5. **Wires** : les arêtes sont chaînées par leurs sommets (`make_wire`).

La section suffit pour un plan de coupe ou un contour d'usinage. C'est aussi la
première étape des booléens.

---

## 5. Découpe (`gbs-brep/split.h`)

### 5.1 Découpe d'une face

```cpp
auto split_face(Model<T> &m, FaceId f, std::span<const EdgeId> cuts, T tol) -> BuildResult<std::vector<FaceId>>;
```

Les arêtes de coupe ont une pcurve sur la face, et viennent de la section.

1. **Points de coupe sur le bord.** Les arêtes du bord sont coupées par
   `split_edge` aux points d'arrivée des coupes. Elles restent **partagées** :
   la face voisine reçoit les mêmes morceaux, et le shell reste fermé.
2. **Seam.** Une coupe qui traverse le seam d'une surface périodique est coupée
   au seam. Ses morceaux sont décalés d'une période pour revenir dans le
   domaine de la face.
3. **Arrangement.** On construit le graphe des demi-arêtes en (u, v). Autour de
   chaque sommet, les arêtes sont ordonnées par l'angle de leur tangente, puis
   par leur courbure en cas d'égalité. On parcourt ensuite les boucles.
4. **Nouvelles faces.** Les boucles extérieures (aire signée positive)
   deviennent des faces. Chaque trou va dans la boucle extérieure qui le
   contient (§ 3). Toutes les faces partagent la **même surface**
   (`shared_ptr`).
5. **Arêtes pendantes**, qui ne séparent rien : elles sont rapportées et
   ignorées.

### 5.2 Découpe d'une forme

```cpp
auto split(Model<T> &m, ShapeId shape, ShapeId|Surface tool, SplitOptions) -> BuildResult<SplitResult>;
```

C'est la section, puis la découpe de chaque face touchée. Les shells et solides
gardent leurs usages : une face découpée est remplacée par ses morceaux, avec le
même sens. Un solide découpé par une surface reste **un** solide aux faces
subdivisées. Le séparer en deux solides demande des faces de fermeture : c'est
un booléen avec un demi-espace (§ 6).

---

## 6. Booléens (`gbs-brep/boolean.h`, question 1)

### 6.1 API

```cpp
enum class BooleanOp { Fuse, Cut, Common };
auto boolean(Model<T> &m, BooleanOp op, ShapeId a, ShapeId b, BooleanOptions) -> BuildResult<BooleanResult>;
// a et b : solides ou compounds de solides ; BooleanResult : le résultat, et l'origine de chaque face (a, b, nouvelle)
```

### 6.2 Algorithme

C'est le schéma de la « fusion générale », simplifié :

1. **Section** de A par B (§ 4) ; les arêtes de section sont partagées par
   les deux formes.
2. **Découpe** de toutes les faces touchées de A et de B (§ 5).
3. **Classification** de chaque face découpée de A par rapport à B, et
   inversement : dedans, dehors, ou confondue (même sens ou sens opposé).
4. **Sélection.**
   | Opération | Faces gardées |
   |---|---|
   | fusion | A dehors, B dehors, confondues de même sens une fois |
   | intersection | A dedans, B dedans, confondues de même sens une fois |
   | différence | A dehors, B dedans retournées, confondues de sens opposé une fois |
5. **Assemblage.** Les faces gardées partagent déjà leurs arêtes : les shells
   s'obtiennent par parcours, **sans sewing**. Puis `make_solid` (cavités
   comprises) et `check()`.

### 6.3 Configurations dégénérées traitées

- **Faces planes coplanaires**, très fréquentes (boîtes posées l'une sur
  l'autre, poches débouchantes) : la découpe se fait en 2D dans le plan
  commun, avec les bords des deux faces.
- **Primitives coaxiales confondues** (deux cylindres de même axe et de même
  rayon) : même traitement, dans leur paramétrage commun.
- **Contact par une arête ou un sommet** : les sommets et arêtes sont
  fusionnés à la tolérance d'opération ; les solides qui se touchent
  seulement deviennent deux solides du résultat.

Les autres cas dégénérés rendent `NonGeneric`, avec les deux faces en cause.

---

## 7. Constructions (`gbs-brep/sweep.h`, `gbs-brep/offset.h`)

```cpp
auto make_prism(Model<T> &m, ShapeId profile, vector<T, 3> d) -> BuildResult<ShapeId>;          // arête → face, wire → shell, face → solide
auto make_revolution(Model<T> &m, ShapeId profile, ax1<T, 3> axis, T angle) -> BuildResult<ShapeId>;
auto offset_face(Model<T> &m, FaceId f, T d) -> BuildResult<FaceId>;
```

| Profil | Extrusion | Révolution |
|---|---|---|
| droite | plan | plan (perpendiculaire à l'axe), cylindre (parallèle), cône (sinon) |
| cercle | cylindre (perpendiculaire à la direction) | sphère (centre sur l'axe, plan contenant l'axe), tore (plan contenant l'axe) |
| B-spline | B-spline exacte (`build_extrusion`) | NURBS exacte (`build_revolution`) |
| autre | B-spline approchée et rapportée | idem |

- Les faces latérales sont construites avec leurs arêtes partagées. Les faces
  d'extrémité sont des copies déplacées du profil (`copy_shape`). Le shell est
  donc fermé **sans sewing**.
- Une révolution complète a son seam ; un profil qui touche l'axe donne des
  pôles et des arêtes dégénérées (comme les faces naturelles du palier 1).
- **Face décalée.** Sur une primitive, le décalage est exact : plan translaté,
  cylindre, sphère et tore de rayon modifié. Ailleurs, c'est une
  `SurfaceOffset`. Le paramétrage étant le même, les pcurves sont **reprises
  telles quelles** (§ 9.1 du palier 1). Les arêtes 3D sont recalculées depuis
  ces pcurves.

---

## 8. Python

- Sous-module `gbs.intersect` pour les intersections géométriques.
- `gbs.brep` reçoit `classify`, `section`, `split_face`, `split`, `boolean`
  (et `fuse`, `cut`, `common`), `make_prism`, `make_revolution` et
  `offset_face`.
- Les erreurs lèvent `gbs.brep.BRepError`, comme au palier 1.

---

## 9. Stratégie de test

1. **Géométrie.**
   - Références en formule fermée : droite-quadrique, plan-primitive,
     coaxiaux.
   - Symétrie : `intersect(a, b)` et `intersect(b, a)` donnent les mêmes
     courbes.
   - Cheminement : cylindres croisés de rayons différents ou égaux (courbe de
     Steinmetz, branchement), sphère et cylindre décalé (fenêtre de Viviani),
     surfaces B-spline libres. On mesure l'écart des points aux deux surfaces.
   - Tangences et cas dégénérés reconnus.
2. **BREP.**
   - `check()` passe sur tous les résultats.
   - Volumes connus : boîte moins boîte, boîte moins cylindre traversant,
     sphère moins cylindre, fusion de deux cylindres croisés (formule de
     Steinmetz), boîtes coplanaires.
   - Aires et longueurs de sections connues.
3. **OCCT** (job OCCT de la CI) : une quinzaine de configurations, comparées à
   `BRepAlgoAPI_*` en volume, aire et nombre de faces. Les résultats gbs sont
   écrits en STEP, relus par OCCT et doivent être **valides**
   (`BRepCheck_Analyzer`).
4. **Temps de calcul** comparés à OCCT, rapportés mais sans seuil.

---

## 10. Plan de développement en PR

Chaque PR est accompagnée de sa note d'architecture.

| # | PR | Contenu | Tests | Lignes | h |
|---|---|---|---|---|---|
| 1 | `geom/intersect-curves` | `gbs/intersect/` : boîtes englobantes (§ 2.2), courbe-courbe 2D (§ 2.3), portions communes | formules fermées, B-splines croisées, tangences, portions communes | 600 | 7 |
| 2 | `geom/intersect-curve-surface` | courbe-surface (§ 2.4) : droite contre les primitives (quartique du tore), cercle-plan, cas général | formules, B-splines, courbe posée sur la surface | 600 | 7 |
| 3 | `geom/intersect-analytic` | surface-surface des paires de primitives (§ 2.5), pcurves exactes ou interpolées, cas confondus (§ 2.7) | chaque ligne du tableau du § 2.5 | 700 | 8 |
| 4 | `geom/intersect-marching` | cheminement (§ 2.6) : départs, cônes de normales, pas adaptatif, singularités, courbes et pcurves | Steinmetz, Viviani, B-splines libres, tore-cylindre, écarts mesurés | 900 | 12 |
| 5 | `brep/classify` | point/face, point/solide, face/solide (§ 3) | boîte, cylindre, sphère percée, points sur le bord, rayons ambigus | 400 | 5 |
| 6 | `brep/section` | section de deux formes ou par une surface (§ 4) | plan à travers boîte, cylindre et tore ; cylindres croisés ; longueurs | 700 | 9 |
| 7 | `brep/split` | découpe de faces et de formes (§ 5), seams | face découpée par un cercle (trou), cylindre coupé au travers du seam, `check()` | 800 | 10 |
| 8 | `brep/sweep-offset` | extrusion, révolution, face décalée (§ 7) | profils droits, circulaires, B-spline ; révolution touchant l'axe ; volumes | 600 | 7 |
| 9 | `brep/boolean` | booléens (§ 6), coplanaires et coaxiaux | volumes connus, boîtes coplanaires, `NonGeneric` | 900 | 12 |
| 10 | `brep/occt-python` | comparaison OCCT, Python, documentation | job OCCT, pytest | 500 | 6 |

Total : environ 6 700 lignes, **83 h**. Ordre : 1 → 2 → 3 → 4, puis 5 (après
2), 6 (après 3 et 4), 7, 9, 10. La PR 8 ne dépend que du palier 1 et peut
avancer à tout moment. Sans les booléens, on retire la PR 9 et une partie de
la PR 10 : 66 h.

---

## 11. Questions ouvertes

Pour chaque question : **R** = recommandation, **A** = alternatives.

1. **Booléens dans le palier 2 ?**
   R : oui, pour les configurations génériques, les faces coplanaires et les
   primitives coaxiales confondues (§ 6.3). Les autres cas dégénérés rendent
   `NonGeneric` et restent à OCCT.
   A : s'arrêter à la section et à la découpe, comme le prévoyait le palier 1 ;
   ou viser des booléens complets, tangences comprises (bien plus coûteux et
   long à fiabiliser).

2. **Intersections dans `gbs/intersect/`, indépendantes du BREP ?**
   R : oui. Elles servent aussi aux courbes et surfaces seules (Python
   compris).
   A : dans `gbs-brep/`.

3. **Courbes d'intersection exactes quand c'est possible ?**
   R : oui : droites, cercles et ellipses en primitives, pcurves exactes quand
   le paramétrage est affine. Elles sont écrites en entités analytiques en
   STEP, et les opérations suivantes en profitent.
   A : toujours une B-spline, plus simple, mais l'exactitude est perdue.

4. **Cheminement pour le cas général ?**
   R : oui, avec des départs trouvés par subdivision et le critère des cônes
   de normales. C'est la méthode des noyaux industriels, précise et rapide.
   A : subdivision seule jusqu'à des carreaux plans, puis chaînage : plus
   simple et plus robuste aux boucles, mais coûteux et moins précis.

5. **Cas dégénérés non traités : erreur `NonGeneric` ?**
   R : oui, plutôt qu'un résultat faux silencieux.
   A : un meilleur effort avec avertissement dans le rapport.

6. **Paraboles et hyperboles (plan / cône) ?**
   R : B-spline par cheminement, avec l'écart mesuré. Ces cas sont rares, et
   gbs n'a pas ces courbes.
   A : ajouter `ParabolaCurve` et `HyperbolaCurve` aux primitives (STEP les a),
   environ 3 h.

7. **Pcurves des courbes de cheminement interpolées, avec écart mesuré ?**
   R : oui : une B-spline 2D sur le même paramètre que la courbe 3D.
   A : une `CurveOnSurface` sans courbe 3D indépendante (exacte sur une face,
   mais pas sur l'autre).

8. **Point dans un solide par lancer de rayon déterministe ?**
   R : oui.
   A : nombre d'enroulement par intégrale sur les faces (robuste près du bord,
   mais coûteux).

9. **Balayage le long d'un chemin (*pipe*) ?**
   R : plus tard. Il demande un repère mobile (Frenet ou rotation minimale) et
   donne des surfaces approchées.
   A : l'ajouter à la PR 8, environ 4 h.

10. **Tolérance d'opération unique, avec remontée des tolérances des entités ?**
    R : oui, comme le sewing du palier 1.
    A : une tolérance par opérande.

11. **Les extrema existants (`extrema.h`, nlopt) ?**
    R : les garder tels quels pour leurs utilisateurs actuels. Le nouveau code
    ne s'en sert pas.
    A : les réécrire sur le nouveau noyau d'intersection (hors périmètre).

12. **Python : sous-module `gbs.intersect` ?**
    R : oui, avec les opérations BREP dans `gbs.brep`.
    A : tout dans `gbs.brep`.

---

## 12. Bilan de réalisation

### 12.1 PR et notes d'architecture

| # | PR | Note |
|---|---|---|
| 1 | courbe-courbe, boîtes (avec ce document) | [brep2_pr01_intersect_curves.md](brep2_pr01_intersect_curves.md) |
| 2 | courbe-surface | [brep2_pr02_curve_surface.md](brep2_pr02_curve_surface.md) |
| 3 | surface-surface, paires de primitives | [brep2_pr03_analytic.md](brep2_pr03_analytic.md) |
| 4 | surface-surface par cheminement | [brep2_pr04_marching.md](brep2_pr04_marching.md) |
| 5 | classification | [brep2_pr05_classify.md](brep2_pr05_classify.md) |
| 6 | section | [brep2_pr06_section.md](brep2_pr06_section.md) |
| 7 | découpe | [brep2_pr07_split.md](brep2_pr07_split.md) |
| 8 | extrusion, révolution, face décalée | [brep2_pr08_sweep_offset.md](brep2_pr08_sweep_offset.md) |
| 9 | booléens génériques | [brep2_pr09_boolean.md](brep2_pr09_boolean.md) |
| 9 b | booléens avec faces confondues | [brep2_pr09b_coincident.md](brep2_pr09b_coincident.md) |
| 10 | comparaison OCCT, Python, documentation | [brep2_pr10_occt_python.md](brep2_pr10_occt_python.md) |

### 12.2 Validation par OpenCascade

- **21 booléens** (fusion, différence, intersection) de formes construites par
  OCCT :
  - boîtes qui se chevauchent ;
  - boîte percée par un cylindre ;
  - boîte et sphère inclinée ;
  - cylindres croisés, dont les seams se coupent sur la courbe d'intersection ;
  - boîtes coplanaires et boîtes empilées ;
  - cylindres coaxiaux.

  Écrits par gbs et relus par OCCT, les résultats sont valides
  (`BRepCheck_Analyzer`). Ils ont les mêmes volumes (2e-3) et les **mêmes
  nombres de solides et de faces** que `BRepAlgoAPI`, chaque opération en
  quelques dizaines de millisecondes.
- **Sections** d'une boîte, d'un cylindre et d'une sphère inclinée par des plans
  obliques : même longueur que `BRepAlgoAPI_Section` à 1e-4 près, wires fermés.
- **Balayages.** Une plaque percée est extrudée en biais, tournée sur un tour
  complet (avec sa cavité torique) puis sur 2 rad. On obtient les mêmes
  volumes et nombres de faces que `BRepPrimAPI_MakePrism` et `MakeRevol`.

### 12.3 Décisions révisées en cours de route

| Sujet | Document | Réalisé |
|---|---|---|
| découpe des PR | booléens et coplanaires en PR 9 | PR 9 (génériques) et 9 b (faces confondues) |
| départs du cheminement | bords, puis cônes de normales | subdivision seule, carreaux parallèles redécoupés |
| surfaces périodiques | espace déroulé | courbes coupées aux seams des boîtes (pcurves continues) |
| constructions | arêtes partagées par construction | faces naturelles cousues |
| faces décalées | `SurfaceOffset` | primitives exactes, `OffsetSurface` avec dérivées |

### 12.4 Limites connues

- **Restriction des courbes analytiques aux boîtes.** Elle se fait par
  échantillonnage (512 pas) et bichotomie : une partie plus courte qu'un pas
  peut être manquée.
- **Pôles.** Le cheminement à travers le pôle d'une sphère livre deux morceaux
  non raccordés (PR 6, § 6).
- **Configurations dégénérées refusées** (`NonGeneric`) : tangences, contacts
  ponctuels, arêtes d'un opérande posées sur une face non confondue de
  l'autre, points de branchement.
- **Faces confondues** : seulement entre primitives (plans, coaxiales de même
  rayon).
- **Pas de balayage le long d'un chemin**, ni de coques épaisses.
- **Arêtes et sommets inutilisés.** Après une opération, ils restent dans
  l'arène jusqu'à un `compact()`.
