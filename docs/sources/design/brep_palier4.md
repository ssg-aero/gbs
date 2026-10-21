# Noyau BREP natif — conception du palier 4 : balayage, congés, coques épaisses

| | |
|---|---|
| Statut | **À valider** |
| Périmètre | exemples visuels et documentation des booléens et du maillage (PR 0) ; balayage d'un profil le long d'un chemin (*pipe*) ; décalage d'un solide et coques épaisses (évidement) ; congés à rayon constant et chanfreins le long d'arêtes ; Python ; comparaison à OpenCascade ; exemples visuels de chaque opération |
| Dépend de | palier 1 (modèle, constructeurs, couture), palier 2 et 2 bis (intersections, section, découpe, booléens, `offset_face`, `make_prism`, `make_revolution`), palier 3 (maillage, pour les contrôles visuels) |
| Hors périmètre | congés à rayon variable, raccords de sommet quelconques (plus de trois arêtes, rayons différents), congés de faces (*face blend*), balayage à section variable ou guidé par deux rails, réparation de géométries auto-intersectées, décalage arrondi (joints en arcs) |

Ce document décrit *quoi* construire et *pourquoi*, comme
[brep_palier2.md](brep_palier2.md) et [brep_palier3_mesh.md](brep_palier3_mesh.md).
Chaque PR a sa note d'architecture `brep4_prNN_*.md`.

Les trois opérations étaient promises « au palier 3 » par les paliers 2 et 2 bis ;
le palier 3 est devenu le maillage, elles forment donc ce palier 4.

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Congés | **écrits dans le noyau**, ce qui révise `brep_core.md` (« congés : restent délégués à OCCT ») | les laisser à OCCT via `gbs-occt` |
| Repère du balayage | **rotation minimale** (double réflexion de Wang et al.) par défaut ; Frenet et binormale fixe en options ; chemin fermé : la torsion résiduelle répartie le long du chemin | Frenet seul (il tourne aux inflexions et s'annule sur les droites) |
| Surfaces balayées | **exactes** quand un cas connu s'y prête (droite ou cercle le long d'une droite ou d'un cercle : plan, cylindre, cône, tore, révolution) ; **B-spline approchée** sinon, par sections placées puis *skinning*, raffinée jusqu'à `approx_tol` mesurée sur la surface exacte | B-spline approchée toujours |
| Chemin | une arête ou un contour **G1** (tangente continue) ; un coin est refusé à la première version | onglets aux coins (*miter*) |
| Ligne des centres d'un congé | intersection des deux surfaces **décalées** de r (décalages exacts des primitives, formes fermées des intersections du palier 2 ; cheminement sinon) | équations du congé résolues par un cheminement dédié |
| Surface du congé | **exacte** pour les paires de primitives qui le permettent (plan-plan : cylindre ; plan-cylindre ou plan-cône d'axe normal au plan, cylindre-cylindre coaxiaux : tore) ; **B-spline rationnelle approchée** sinon (arcs de cercle exacts, rationnels de degré 2 en travers, approchés le long) | approchée toujours |
| Topologie d'un congé | **locale** : faces adjacentes découpées le long des courbes de contact, bande retirée, face du congé insérée, faces d'extrémité découpées par la surface du congé, puis couture | booléen avec un solide outil (contacts tangents tout le long : la configuration la plus fragile des booléens) |
| Raccords de sommet | arêtes enchaînées tangentes (une seule ligne des centres) ; **coin de trois arêtes de même rayon à angles droits** (sphère) ; refusés sinon (`NonGeneric`) | raccords quelconques (*vertex blends* de n arêtes) |
| Décalage d'un solide | chaque face décalée (`offset_face`), les arêtes refaites par **intersection** des surfaces décalées voisines (joints vifs), surfaces prolongées au besoin | joints arrondis (un congé de rayon d le long de chaque arête) |
| Coque épaisse | `cut(solide, solide décalé vers l'intérieur)`, les faces retirées décalées de 0 : coïncidentes, le booléen du palier 2 les ouvre | construction directe des deux peaux et des faces de tranche |
| Auto-intersections | **détectées et refusées** (rayon ou épaisseur au-delà d'un rayon de courbure, face qui disparaît) | réparées |

---

## 0 bis. Exemples visuels et documentation (PR 0)

Les booléens des paliers 2 et 2 bis sont validés par des volumes, des nombres de
faces et OpenCascade, jamais **vus**. Avant le palier 4, une PR de documentation
et d'exemples les montre ; chaque PR du palier 4 ajoute ensuite les siens.

| Livrable | Contenu |
|---|---|
| `docs/sources/examples/gbs_examples.py` | module commun : construction des cas (opérandes lus en STEP ou construits par `make_prism`, `make_revolution`, placés par `copy_shape`), maillage (palier 3), rendu PyVista : une couleur par face (pour voir les découpes), arêtes du BREP en traits, opérandes en transparence |
| `docs/sources/examples/booleans.ipynb`, `meshing.ipynb` | notebooks interactifs : chaque cas, ses opérandes, `fuse`, `cut`, `common` en 3D ; maillage d'une forme, flèche, comparaison des options |
| `docs/sources/examples/gallery.py` | rendu hors écran de tous les cas en PNG (une image par cas et par opération, les opérandes en vignette) et d'une page `gallery.md` qui les rassemble |
| CI | la galerie générée sur Linux (xvfb) et publiée en artefact du build : les régressions visuelles se voient d'un coup d'œil |
| Sphinx | une page « Exemples » qui renvoie aux notebooks et à la galerie |

Les cas reprennent ceux de la comparaison OCCT des paliers 2 et 2 bis (boîtes,
boîte et cylindre, boîte et sphère, cylindres croisés, opérandes coplanaires,
empilés, coaxiaux, contacts tangents, arêtes sur faces, points de branchement,
pôles, coïncidences partielles) et des maillages du palier 3.

Une **présentation des avancées de GBS** (article LinkedIn) suit cette PR et
s'appuie sur ses rendus ; elle n'est pas dans le dépôt.

---

## 1. Objectifs et non-objectifs

### 1.1 Objectifs

1. **Balayage le long d'un chemin** (`make_pipe`) : une arête, un contour ou une
   face balayés le long d'un contour G1 ; face balayée : solide fermé par deux
   copies du profil, ou par aucune si le chemin est fermé.
2. **Décalage d'un solide** (`offset_solid`) à une distance d (positive vers
   l'extérieur), joints vifs, et **coque épaisse** (`make_thick`) : évidement
   d'un solide à une épaisseur t, ouvert par des faces retirées.
3. **Congés** à rayon constant (`fillet`) et **chanfreins** (`chamfer`, deux
   distances) le long d'arêtes choisies d'un solide, convexes ou concaves, avec
   les enchaînements tangents et les coins de trois arêtes.
4. Résultats **valides** (`check`), fermés, et **comparés à OpenCascade**
   (`BRepOffsetAPI_MakePipeShell`, `BRepOffsetAPI_MakeThickSolid`,
   `BRepFilletAPI_MakeFillet`, `BRepFilletAPI_MakeChamfer`) : volumes, nombres
   de faces.
5. **Python** : toutes les opérations dans `gbs.brep`.

### 1.2 Non-objectifs

- Congés à rayon variable, congés de faces (entre deux faces non adjacentes),
  raccords de sommet au-delà des cas du § 4.6.
- Balayage à section variable, guidé (deux rails), à loi de torsion imposée.
- Décalage arrondi, décalage de coques ouvertes (surfaces épaisses).
- Réparation : une opération qui s'auto-intersecte est refusée, pas corrigée.

---

## 2. Balayage le long d'un chemin

### 2.1 API

```cpp
enum class PipeFrame : std::uint8_t { RotationMinimizing, Frenet, FixedBinormal };

template <std::floating_point T>
struct PipeOptions
{
    T tol{brep_default_tolerance<T>};
    T approx_tol{T(1e-5)};                  // surface approchée : écart à la surface exacte
    PipeFrame frame{PipeFrame::RotationMinimizing};
    point<T, 3> binormal{};                 // pour FixedBinormal
};

auto make_pipe(Model<T> &, const ShapeId &profile, const ShapeId &path, const PipeOptions<T> & = {})
    -> BuildResult<SweepResult<T>>;
```

Le profil est placé tel quel : le repère au début du chemin est le repère de
référence, le profil n'a pas à être dans le plan normal (comme
`BRepOffsetAPI_MakePipe`).

### 2.2 Repère mobile

- **Rotation minimale** : double réflexion (Wang, Jüttler, Zheng, Liu 2008) en
  des points du chemin assez serrés pour que l'écart de la rotation reconstruite
  soit sous `approx_tol` ; entre eux, interpolation de l'angle. Chemin fermé :
  l'angle résiduel au retour est réparti linéairement selon l'abscisse.
- **Frenet** : refusé là où la courbure s'annule (droite, inflexion).
- **Binormale fixe** : le repère garde une direction donnée (comme
  `BRepOffsetAPI_MakePipeShell::SetMode(gp_Dir)`).

### 2.3 Surfaces

Pour chaque arête du profil et chaque arête du chemin, une face :

| Profil | Chemin | Surface |
|---|---|---|
| segment | segment | plan (exact) |
| cercle | segment | cylindre ou cône (exacts, comme `make_prism`) |
| segment ou cercle | cercle, repère de rotation minimale | révolution (`revolution_surface` de `make_revolution` : cylindre, cône, sphère, tore) |
| autre | autre | B-spline : le profil placé par le repère en n stations du chemin, puis *skinning* (`loft`) ; n doublé tant que l'écart, mesuré entre les stations à la surface exacte (profil placé par le repère au même paramètre), dépasse `approx_tol` |

`SweepResult::approximated` compte les faces approchées (comme
`make_prism`, dont `approx_tol` est enfin utilisé).

### 2.4 Topologie

Comme `make_prism` : faces naturelles cousues (`sew`), une coque par profil
connexe ; profil face : deux bouchons (copies du profil placées par le repère
aux deux bouts), solide avec les trous du profil en vides. Un chemin fermé ne
donne pas de bouchons ; le profil doit alors retomber sur lui-même (contrôlé).

### 2.5 Refus

- Chemin non G1 (angle entre tangentes aux jonctions au-delà de `tol` angulaire)
  : `NonGeneric` à la première version.
- Profil plus grand que le rayon de courbure du chemin (la surface se replie) :
  détecté en comparant la distance du profil à l'axe au rayon de courbure
  minimal du chemin ; `NonGeneric`.

---

## 3. Décalage d'un solide et coques épaisses

### 3.1 API

```cpp
auto offset_solid(Model<T> &, SolidId, T d, const OffsetOptions<T> & = {}) -> BuildResult<SolidId>;
auto make_thick(Model<T> &, SolidId, std::span<const FaceId> removed, T thickness, const OffsetOptions<T> & = {})
    -> BuildResult<SolidId>;
```

`thickness` positive : vers l'intérieur (évidement) ; négative : vers
l'extérieur.

### 3.2 Décalage

1. Chaque face décalée de d (`offset_face` : exact pour les primitives,
   `OffsetSurface` sinon). `OffsetSurface` reçoit des **dérivées exactes**
   (premières à partir des secondes de S, quand S les a) : le cheminement en a
   besoin (ses dérivées par différences finies sont la limite notée en
   `brep2_pr08`).
2. Chaque arête du solide devient l'**intersection des deux surfaces décalées**
   de ses faces (formes fermées entre primitives, cheminement sinon), restreinte
   au voisinage de l'arête d'origine (boîte des pcurves élargie de d) ; chaque
   sommet, l'intersection des trois surfaces décalées (ou des arêtes refaites).
3. Les faces sont refaites sur les surfaces décalées avec les nouvelles arêtes,
   **prolongées** quand une surface décalée doit dépasser la boîte de son
   original (joint vif d'une arête convexe décalée vers l'extérieur) : les
   primitives se prolongent d'elles-mêmes ; une B-spline est prolongée par
   extrapolation de ses carreaux de bord (C2), et la face refusée si
   l'intersection n'y est pas trouvée.
4. Couture et `check`.

Une face qui **disparaît** (décalage plus grand que sa taille, une face étroite
entre deux arêtes concaves) ou une surface décalée qui se replie (d au-delà du
rayon de courbure, détecté par les courbures principales échantillonnées) :
refusé (`NonGeneric`).

### 3.3 Coque épaisse

`make_thick(s, removed, t)` : le solide décalé de −t, sauf les faces de
`removed`, décalées de 0 ; puis `cut(s, intérieur)`. Les faces retirées sont
coïncidentes dans les deux opérandes : le booléen des faces coïncidentes
(palier 2, `coincident.h`) les retire, et la coque s'ouvre là. Sans face
retirée : un vide fermé (solide creux).

---

## 4. Congés et chanfreins

### 4.1 API

```cpp
auto fillet(Model<T> &, SolidId, std::span<const EdgeId> edges, T radius, const FilletOptions<T> & = {}) -> BuildResult<SolidId>;
auto chamfer(Model<T> &, SolidId, std::span<const EdgeId> edges, T d1, T d2, const FilletOptions<T> & = {}) -> BuildResult<SolidId>;
```

Le solide donné n'est pas modifié (comme les booléens : le résultat est un
nouveau solide, les faces non touchées partagées par copie).

### 4.2 Chaînes

Les arêtes choisies sont groupées en **chaînes** : deux arêtes qui se suivent à
un sommet avec des tangentes continues et les mêmes deux faces (ou des faces
qui se suivent tangentiellement) forment une seule ligne des centres. Une chaîne
est ouverte (deux extrémités) ou fermée (le bord circulaire d'un cylindre).

### 4.3 Ligne des centres

Pour une arête entre les faces F1 et F2 (convexe ou concave, lu sur l'angle des
normales orientées par la coque) : la ligne des centres est l'intersection de
S1 et S2 décalées de r **du côté de la matière** pour une arête convexe, du côté
extérieur pour une concave. Décalages exacts des primitives (`offset_face`) ;
intersection par le palier 2 (formes fermées : plan-plan, plan-cylindre, etc. ;
cheminement sinon), restreinte au tronçon qui borde l'arête.

Les **courbes de contact** sont les pieds de la ligne des centres sur S1 et S2
(c − r n1, c − r n2), avec leurs pcurves : ce sont les coupes des faces.

### 4.4 Surface du congé

| Paire de faces | Surface |
|---|---|
| deux plans | cylindre de rayon r, axe la ligne des centres |
| plan et cylindre (ou cône) d'axe normal au plan | tore |
| deux cylindres coaxiaux, cylindre et plan parallèle à l'axe | cylindre ou tore |
| autre | B-spline rationnelle : en travers, l'arc de cercle exact (rationnel de degré 2) entre les deux points de contact, de centre sur la ligne des centres ; le long, l'approximation des trois courbes (centre, contacts) à `approx_tol`, les poids et poles intermédiaires interpolés |

Le chanfrein est la surface réglée entre les deux courbes de contact, prises
aux distances d1 et d2 de l'arête sur F1 et F2 (intersections de F1 avec S2
décalée de d2, et réciproquement) : plan entre deux plans, cône ou réglée
B-spline sinon.

### 4.5 Topologie

Pour chaque chaîne :

1. F1 et F2 sont **découpées** le long des courbes de contact (`split_face`),
   et la bande entre la courbe de contact et l'arête retirée.
2. À chaque **extrémité** d'une chaîne ouverte, au sommet où arrive une
   troisième face F3 (sommet de degré 3) : la surface du congé est intersectée
   avec F3 (intersection du palier 2), F3 découpée par cette courbe et son coin
   retiré. La face du congé est bornée par les deux contacts et les deux courbes
   d'extrémité.
3. La face du congé est construite sur sa surface avec ces arêtes ; le tout est
   **cousu** en un nouveau solide (`sew`, `make_solid`), contrôlé (`check`).

### 4.6 Raccords de sommet

- **Chaîne tangente** : traitée par § 4.2 (une ligne des centres).
- **Coin de trois arêtes congées de même rayon** (coin de boîte, trois faces
  deux à deux orthogonales à la première version, puis quelconques si la
  ligne des centres des trois converge) : les trois lignes des centres se
  coupent en un point, centre d'une **sphère** de rayon r ; la face de coin est
  la sphère bornée par les trois arcs d'extrémité des congés.
- Autre configuration (deux arêtes congées non tangentes sans troisième, quatre
  arêtes, rayons différents) : `NonGeneric`.

### 4.7 Refus

- Rayon au-delà du rayon de courbure d'une face du côté du congé (la surface
  décalée se replie) ou de la largeur d'une face (le contact sort de la face ou
  traverse une autre arête) : `NonGeneric`.
- Arête non manifold, arête de couture, arête dégénérée : `InvalidId`.

---

## 5. Python

`gbs.brep.make_pipe`, `PipeOptions`, `PipeFrame`, `offset_solid`, `make_thick`,
`fillet`, `chamfer`, `FilletOptions` ; `BRepError` sur les refus.

---

## 6. Validation

| Contrôle | Attendu |
|---|---|
| validité | `check` sans défaut, coques fermées, `signed_volume` positif |
| volumes exacts | balayage d'un disque le long d'un cercle (tore : 2π² R r²), d'un carré le long d'une hélice (Pappus : aire × longueur du chemin, profil normal) ; congé d'une arête de boîte (retire (1 − π/4) r² L) ; coin de boîte à trois congés (formule fermée) ; boîte évidée (différence de deux boîtes) |
| OCCT | mêmes opérations par `BRepOffsetAPI_MakePipeShell`, `BRepOffsetAPI_MakeThickSolid`, `BRepFilletAPI_MakeFillet`, `BRepFilletAPI_MakeChamfer` : volumes à 1e-3 près, mêmes nombres de faces, résultats de gbs lus par OCCT valides (`BRepCheck_Analyzer`) |
| géométrie | congé : tangent aux deux faces le long des contacts (angle des normales sous la tolérance angulaire), à distance r de la ligne des centres ; balayage approché : écart sous `approx_tol` |
| maillage | les résultats maillés par le palier 3 : fermés |

La boîte à 12 arêtes congées lue au palier STEP (`step_pr08_occt_compare.md`)
sert de référence OCCT pour les congés et les coins.

---

## 7. Plan des PR

| # | PR | Contenu | Estimation |
|---|---|---|---|
| 0 | exemples et documentation (avec ce document) | § 0 bis : module d'exemples, notebooks des booléens et du maillage, galerie PNG générée par la CI, page Sphinx | 8 h |
| 1 | repère mobile et balayage | § 2 : `gbs/` repères (rotation minimale, Frenet, binormale fixe), surfaces balayées exactes et approchées ; `make_pipe` ; tests de volumes | 12 h |
| 2 | surfaces décalées | `OffsetSurface` à dérivées exactes ; intersections de surfaces décalées (formes fermées entre primitives décalées, cheminement) ; détection du repli (courbures principales) | 8 h |
| 3 | décalage de solide et coques épaisses | § 3 : `offset_solid` (joints vifs, prolongement des faces), `make_thick` ; comparaison OCCT | 14 h |
| 4 | congés : géométrie | § 4.3, § 4.4 : ligne des centres, contacts, surfaces exactes et approchées, chanfreins | 12 h |
| 5 | congés : topologie | § 4.2, § 4.5 : chaînes, découpe, extrémités, couture ; un ou plusieurs bords indépendants | 14 h |
| 6 | coins, Python, validation, bilan | § 4.6, § 5, § 6 : coin sphérique, Python, comparaison OCCT, bilan du palier | 10 h |
| | | **Total** | **≈ 78 h** |

Chaque PR de 1 à 6 ajoute ses cas aux notebooks et à la galerie.

Ordre : la PR 0 d'abord ; le balayage est indépendant ; les surfaces décalées (PR 2) servent au
décalage (PR 3) et aux congés (PR 4).

---

## 8. Questions à valider

1. **Congés dans le noyau**, révisant `brep_core.md` (« congés : restent
   délégués à OCCT »). *Recommandé* (décidé en choisissant ce palier).
2. **Repère de rotation minimale par défaut**, Frenet et binormale fixe en
   options. *Recommandé.*
3. **Surfaces exactes quand un cas connu s'y prête, B-spline approchée sinon**,
   l'écart mesuré à la surface exacte. *Recommandé.*
4. **Chemin G1 seulement** à la première version ; les coins (onglets) plus tard.
   *Recommandé.*
5. **Topologie locale pour les congés** plutôt qu'un booléen avec un solide
   outil, dont les contacts tangents tout le long sont la configuration la plus
   fragile des booléens. *Recommandé.*
6. **Raccords de sommet limités** aux chaînes tangentes et au coin de trois
   arêtes de même rayon (sphère) ; les autres refusés. *Recommandé.*
7. **Chanfreins avec les congés** (même découpe, surface réglée). *Recommandé.*
8. **Décalage à joints vifs** (intersection des surfaces décalées prolongées),
   pas de joints arrondis. *Recommandé* (comme `BRepOffsetAPI_MakeThickSolid`
   en mode `GeomAbs_Intersection`).
9. **Coque épaisse par booléen** (solide moins son décalage intérieur, faces
   retirées coïncidentes) plutôt qu'une construction directe des deux peaux.
   *Recommandé* : réutilise les booléens et les faces coïncidentes du palier 2.
10. **Auto-intersections refusées** (`NonGeneric`), jamais réparées.
    *Recommandé.*
11. **Exemples dans `docs/sources/examples`** (notebooks et galerie générée par
    la CI), les notebooks non exécutés par la CI (seule la galerie l'est).
    *Recommandé* : les notebooks restent rapides à ouvrir et la CI ne dépend pas
    de Jupyter.
