# Surfaces de subdivision de Catmull–Clark converties en B-splines — conception

| | |
|---|---|
| Statut | **À valider** |
| Périmètre | cage de contrôle polygonale (sommets, faces, arêtes vives, coins) ; surface limite de Catmull–Clark évaluée exactement ; conversion en surfaces B-splines bicubiques, **exacte** sur les régions régulières regroupées en grandes surfaces, **approchée à tolérance** autour des sommets extraordinaires ; bords et arêtes vives ; coque et solide du BREP natif ; interpolation (la surface limite passe par des points ou des sections imposés) ; Python |
| Dépend de | `gbs/bssurf.h` (insertion de nœuds), `gbs/bssapprox.h` (moindres carrés), `gbs-brep` (faces, coques, solides, `check`), `gbs-io/step` (export) |
| Hors périmètre | schémas de Loop (triangles), Doo–Sabin, √3 ; arêtes semi-vives (§ 9) ; surfaces rationnelles en sortie, donc cercles exacts (§ 9) ; cages non manifold ; édition interactive de la cage ; T-splines |

Ce document décrit *quoi* construire et *pourquoi*, comme
[brep_palier4.md](brep_palier4.md). Chaque PR a sa note d'architecture
`subd_prNN_*.md`.

---

## Pourquoi Catmull–Clark

Les formes à embranchements (raccords de tubes, carters, pièces de fonderie) ne
se prêtent pas aux surfaces à quatre côtés : un *loft* ou une surface de Gordon
ne couvre pas un embranchement, et des surfaces raccordées à la main se
rejoignent en G0, au mieux G1 au prix d'un long travail.

Une **cage** de Catmull–Clark (CC) est un maillage de quadrangles de topologie
quelconque ; sa surface limite est lisse partout, C2 sauf aux sommets
extraordinaires où elle est C1. Elle se modèle en déplaçant quelques dizaines de
sommets, dans gbs ou dans un modeleur (Blender, qui suit les règles
d'OpenSubdiv). La CAO, le BREP et STEP veulent des B-splines : d'où la
conversion.

Elle repose sur un fait : **sur un maillage régulier, CC est l'insertion de
nœuds d'une B-spline bicubique uniforme** (Catmull & Clark 1978). La limite sur
une face régulière est donc exactement un carreau bicubique dont les 16 pôles
sont les sommets de son voisinage. Seuls les sommets de valence différente de 4
(*extraordinaires*, EV) et les faces non quadrangulaires donnent une limite qui
n'est pas une B-spline finie : une infinité d'anneaux de carreaux qui
rapetissent vers l'EV (Stam 1998).

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Règles de subdivision | celles du schéma *Catmark* d'**OpenSubdiv** (Pixar) : CC d'origine, bords « bord et coins », arêtes vives et coins marqués | règles propres à gbs (une cage faite dans Blender ou Maya donnerait une autre surface) |
| Implémentation | **écrite dans le dépôt**, en en-têtes, sans dépendance ; **pas d'OpenSubdiv** dans le build ni en CI : des valeurs de référence calculées une fois par OpenSubdiv, avec un petit programme lancé à la main, sont versionnées dans `tests/in/subd/` (§ 10) | OpenSubdiv lié (bibliothèque compilée, à installer partout) ; parties copiées (les règles tiennent en quelques lignes, l'évaluation de `Far` est grande et approchée aux EV) |
| Cage | sommets, faces en listes d'indices (polygones quelconques, sens trigonométrique vu de l'extérieur), arêtes vives, coins ; connectivité en **demi-arêtes à indices** construite et vérifiée, cage non manifold ou mal orientée **refusée** explicitement | `gbs-mesh/halfedgemesh.h` (ébauche à références, `float`, pour CUDA) ; OpenMesh |
| Évaluation exacte | **subdivision locale** jusqu'à un sous-carreau régulier, évalué comme B-spline ; aux sommets irréguliers, les **masques limites** (formes fermées, sinon vecteurs propres de l'anneau) | décomposition propre de Stam (même résultat, tables par valence ; possible plus tard pour la vitesse) ; masques d'OpenSubdiv (approchés aux darts et en travers des arêtes vives irrégulières) |
| Régions régulières | **exactes** (aux arrondis près), **regroupées** : les faces régulières découpées en rectangles le long des séparatrices des EV, une B-spline bicubique à nœuds uniformes par rectangle, C2 à l'intérieur | un carreau par face (des centaines de faces BREP) |
| Faces irrégulières | une B-spline bicubique par face touchant un EV (n par n-gone), **nœuds gradués vers l'EV** (½, ¼, ⅛…) ; rangées de pôles du côté des régions régulières **imposées** (C2 exact avec elles), plan tangent **exact** à l'EV ; profondeur augmentée jusqu'à `approx_tol` en position et `angle_tol` en normale, **mesurées sur la limite exacte** | carreaux de Gregory ou ACC (Loop–Schaefer) : pas des B-splines ; constructions G1 exactes (Peters) : forme imposée qui s'écarte de CC |
| Continuité | C2 dans et entre les régions régulières ; G0 **exact** partout (bords partagés pôle à pôle) ; G1 exact aux EV, **G1 à `angle_tol` près** le long des arêtes qui en partent, mesuré et rapporté | G1 exact partout (contraintes non linéaires ou degré élevé) |
| BREP | une face par surface, **arêtes partagées construites depuis la topologie de la cage** (pas de couture a posteriori), cage fermée → solide | `sew` des surfaces produites |
| Interpolation | positions limites **linéaires** en la cage : moindres carrés creux (Halstead et al. 1993), avec un terme de lissage | aucune (la cage seule, la surface ne passe par rien) |
| Tolérances | `approx_tol` (1e-5) et `angle_tol` (0,1°) sont des **options de chaque appel**, Python compris, rapportées avec les écarts mesurés : réglables selon la géométrie, et leurs défauts révisables d'une PR à l'autre | constantes |

---

## 1. Objectifs et non-objectifs

### 1.1 Objectifs

1. **Évaluer exactement** la surface limite de CC d'une cage en tout point
   (face, u, v), dérivées premières et secondes comprises (à l'EV même, la
   position et les tangentes seules).
2. **Convertir** la surface limite en surfaces B-splines bicubiques
   (`BSSurface<T, 3>`, non rationnelles) : exactes là où c'est possible, à
   tolérance ailleurs, avec l'écart mesuré et rapporté.
3. **Peu de faces** : les faces régulières regroupées en grandes surfaces, au
   niveau de CATIA et SolidWorks plutôt que d'un carreau par face.
4. **BREP valide** : coque ouverte ou solide fermé, `check()` sans défaut,
   exportable en STEP.
5. **Bords et arêtes vives** : bords libres (courbe de bord B-spline cubique
   exacte), arêtes vives et coins marqués.
6. **Interpolation** : faire passer la surface limite par des points ou des
   courbes imposés.
7. **Python** : cage, évaluation, conversion, BREP.

### 1.2 Non-objectifs

- Un modeleur de cage interactif ; la cage se construit par programme ou se
  lit (OBJ) depuis un modeleur.
- Le rendu en temps réel de la surface de subdivision.
- Les autres schémas (Loop, Doo–Sabin, √3).

---

## 2. La cage

```cpp
namespace gbs::subd
{
    template <std::floating_point T>
    struct Cage
    {
        points_vector<T, 3> vertices;
        std::vector<std::vector<std::size_t>> faces;     ///< sens trigonométrique vu de l'extérieur, 3 sommets ou plus
        std::vector<std::array<std::size_t, 2>> creases; ///< arêtes vives
        std::vector<std::size_t> corners;                ///< sommets interpolés par la limite
    };
}
```

La connectivité (`Topology`) se construit une fois : demi-arêtes à indices,
sommets de bord, valences, étiquettes. Elle refuse, avec une `CageError`
(`ValueError` en Python) qui désigne l'élément fautif :

- une arête partagée par trois faces ou plus, un sommet dont les faces forment
  plusieurs éventails (non manifold) ;
- deux faces voisines d'orientations contraires ;
- une face à sommet répété, un sommet isolé, un indice hors de la cage.

**Règles** (Catmark d'OpenSubdiv) :

| Élément | Règle |
|---|---|
| point de face | centre des sommets de la face |
| point d'arête lisse | (v₀ + v₁ + f₀ + f₁) / 4 |
| point d'arête vive ou de bord | milieu |
| sommet lisse de valence n | ((n − 2) V + Σ eⱼ / n + Σ fⱼ / n) / n |
| sommet sur deux arêtes vives ou de bord | (P₋ + 6 V + P₊) / 8 |
| coin (marqué, trois arêtes vives ou plus, ou sommet de bord d'une seule face) | fixe |

Une courbe d'arêtes vives ou de bord a pour limite la **B-spline cubique
uniforme** de ses sommets, indépendante du reste de la cage.

---

## 3. Évaluation exacte de la limite

- **Face régulière** : un quadrangle dont les quatre sommets sont réguliers
  (intérieur de valence 4 ; de bord à deux faces ; sur une arête vive droite
  avec deux faces de chaque côté ; coin d'une seule face). Sa limite est le
  carreau bicubique uniforme de ses 16 voisins ; le long d'un bord ou d'une
  arête vive, les points au-delà (absents, ou de l'autre côté de l'arête) sont
  remplacés par les **points fantômes** 2B − I, réfléchis dans chaque direction
  à un coin. C'est exact : la règle de bord conserve cette relation à chaque
  étape.
- **Sommet irrégulier** : position et tangentes limites par les masques de
  Halstead et al. (1993) à un sommet lisse ; ailleurs (dart, arête vive, coin)
  par les vecteurs propres de la matrice de subdivision de son secteur, là où
  les masques d'OpenSubdiv ne sont pas exacts ; pas de plan tangent donné
  quand son spectre est dégénéré (`subd_pr01_cage.md`, § 5.5).
- **Point quelconque** (face, u, v) d'une face irrégulière : on subdivise le
  voisinage de la face (les faces qui touchent ses sommets, une dizaine) et on
  descend dans le sous-carreau qui contient (u, v), jusqu'à ce qu'il soit
  régulier ; environ log₂(1 / max(u, v)) étapes, comme chez Stam (1998), qui
  fait le même calcul par la base propre. Les dérivées suivent (facteur 2 par
  étape et par ordre de dérivation). Chaque quart est pris dans un repère qui
  lui est propre, translaté et mis à l'échelle par une puissance de 2 (exact) :
  sinon, près d'un sommet extraordinaire, les arrondis des coordonnées de la
  cage l'emportent, puis les nombres sortent des flottants.
- **n-gone** : n carreaux, un par quadrangle de sa première subdivision, chacun
  sur [0, 1]², comme les faces *ptex* d'OpenSubdiv.

```cpp
template <std::floating_point T>
class LimitSurface
{
public:
    explicit LimitSurface(const Cage<T> &cage);                       // Topology, ou BuildError
    auto patch_count() const -> std::size_t;                          // une par quadrangle, n par n-gone
    auto patch(std::size_t face, std::size_t corner = 0) const -> std::size_t;
    auto value(std::size_t patch, T u, T v, std::size_t du = 0, std::size_t dv = 0) const -> point<T, 3>;
    auto normal(std::size_t patch, T u, T v) const -> point<T, 3>;
    auto vertex_limit(std::size_t vertex) const -> point<T, 3>;      // masque limite
};
template <std::floating_point T>
auto refine(const Cage<T> &cage, std::size_t levels = 1) -> Cage<T>; // les étiquettes suivent
```

---

## 4. Régions régulières : surfaces exactes regroupées

### 4.1 Découpe en rectangles

Les faces régulières sont découpées en **rectangles** de quadrangles : des
séparatrices partent de chaque EV, de chaque coin et des coins des faces
irrégulières, et suivent les lignes du maillage jusqu'à un bord, une arête
vive, une face irrégulière ou une séparatrice déjà tracée (*motorcycle graph*,
Eppstein et al. 2008). Les jonctions en T sont admises : le côté d'un rectangle
peut toucher plusieurs voisins, et le BREP le coupe en autant d'arêtes.

Une bande de faces régulières qui fait le tour (un anneau d'un tube sans EV)
donne une surface **fermée** : une B-spline dont les bords coïncident et une
arête de couture, comme un cylindre.

### 4.2 Surface d'un rectangle

Un rectangle de m × n faces : ses pôles sont la grille de (m + 3) × (n + 3)
sommets de la cage (points fantômes au bord et le long des arêtes vives), ses
nœuds les entiers −3 … m + 3 et −3 … n + 3. Insérer les nœuds 0 et m (puis 0 et
n) jusqu'à la multiplicité 3 la rend **serrée** (*clamped*), avec (m + 3) ×
(n + 3) pôles et des nœuds intérieurs simples : elle est la limite exacte, C2,
aux arrondis de l'insertion près. Une face de cage couvre un intervalle unité
de chaque paramètre.

---

## 5. Faces irrégulières : approximation contrôlée

### 5.1 Structure de la limite près d'un EV

Après une subdivision, une face dont le coin (0, 0) est un EV a trois
sous-carreaux réguliers, donc exactement bicubiques, et un seul irrégulier,
[0, ½]² ; et ainsi de suite. À la profondeur k, la limite est **exactement
bicubique par morceaux** sur toute la face sauf le carré [0, 2⁻ᵏ]² qui touche
l'EV :

```
v
1 ┌───────────────┬───────────────┐
  │     exact     │     exact     │
½ ├───────┬───────┼───────────────┤
  │ exact │ exact │               │
¼ ├───┬───┼───────┤     exact     │
  │ex.│ex.│ exact │               │
⅛ ├───┼───┤       │               │
  │ ◆ │ex.│       │               │
0 └───┴───┴───────┴───────────────┘ u
  0   ⅛   ¼       ½               1
```

À la profondeur 3 : ◆, le carré [0, ⅛]² qui touche l'EV, est le seul où la
limite n'est pas un polynôme bicubique.

### 5.2 Surface d'une face irrégulière

Une `BSSurface` bicubique sur [0, 1]², de nœuds **gradués** vers chaque coin
extraordinaire (0, 2⁻ᵏ, …, ¼, ½, 1 dans chaque direction) : son espace contient
la limite exacte hors du coin, l'approximation ne porte que sur un carré qui
rapetisse de moitié à chaque niveau.

Les pôles sont ajustés aux moindres carrés sur des échantillons de la limite
exacte, sous contraintes **linéaires** :

| Contrainte | Effet |
|---|---|
| côté partagé avec une région régulière : ses trois rangées de pôles sont celles du voisin exact (nœuds insérés) | courbe de bord identique, C2 exact |
| côté partagé avec une autre face irrégulière : une courbe de bord commune, ajustée une fois, imposée aux deux | G0 exact |
| coin EV : pôle = position limite, ses deux voisins dans le plan tangent limite | G1 exact à l'EV |
| échantillons de normale : ∂S/∂u · N = 0, ∂S/∂v · N = 0 avec N la normale limite exacte (pondérés) | G1 le long des arêtes de l'EV, à `angle_tol` près |

Le niveau k augmente tant que l'écart en position dépasse `approx_tol` ou
l'angle entre normales dépasse `angle_tol`, mesurés sur des points **distincts**
de ceux de l'ajustement (milieux et quarts des cellules, bords). La limite n'est
que C1 à l'EV et sa courbure s'y comporte mal (non bornée ou dégénérée selon la
valence, Peters & Reif 2008) : l'approximation suit la limite, elle ne la corrige
pas.

Les deux tolérances se règlent à chaque appel, selon la géométrie (plus serré
pour une pièce usinée, plus lâche pour un avant-projet) ; les écarts atteints
sont rapportés par surface. Les valeurs par défaut sont un point de départ, à
revoir sur des pièces réelles au fil des PR.

```cpp
template <std::floating_point T>
struct ToBSplineOptions
{
    T approx_tol{T(1e-5)};             ///< écart en position à la limite exacte
    T angle_tol{T(0.1) * pi / T(180)}; ///< écart des normales (radians)
    std::size_t max_level{20};          ///< profondeur maximale des nœuds gradués
    bool merge_regions{true};           ///< false : un carreau par face (diagnostic)
};

template <std::floating_point T>
struct SubdPatch
{
    BSSurface<T, 3> surface;
    bool exact;                         ///< région régulière : la limite elle-même
    std::vector<std::size_t> patches;   ///< carreaux de la limite couverts (ligne par ligne)
    T deviation{}, angle{};             ///< mesurés, nuls si exact
};
```

---

## 6. BREP

`gbs-brep/subd.h` :

```cpp
template <std::floating_point T>
[[nodiscard]] auto make_subd_shape(Model<T> &m, const subd::Cage<T> &cage, subd::ToBSplineOptions<T> opts = {})
    -> BuildResult<SubdShape>; // shell ou solid, face de chaque SubdPatch, écarts max, nombre de surfaces exactes et approchées
```

- Une face par `SubdPatch`, bornée par son rectangle de paramètres
  (`natural_bounds`), pcurves exactes (isoparamétriques).
- Les **arêtes** se construisent depuis la cage : un côté de surface est coupé
  aux jonctions en T, chaque morceau est une arête **partagée** par les deux
  faces qui le bordent, sa courbe l'isoparamétrique de l'une d'elles —
  identique pôle à pôle à celle de l'autre (§ 4, § 5). Pas de couture a
  posteriori : la topologie est connue.
- Cage sans bord → `make_solid` ; avec bords → coque ouverte.
- Une erreur ne laisse rien dans le modèle (le rollback des booléens).

---

## 7. Interpolation et constructeurs de cage

- La position limite de chaque sommet est une combinaison linéaire fixe de son
  voisinage, comme tout point (face, u, v) : L = M P, M creuse. Faire passer
  la limite par des points Q aux sommets revient à résoudre M P = Q (Halstead
  et al. 1993) ; des points en (face, u, v) quelconques ajoutent des lignes à M ;
  un terme de lissage (énergie de la cage raffinée) évite les ondulations.
  Moindres carrés creux (`Eigen::SparseQR` ou gradient conjugué ; en-têtes Eigen
  admis).
- **Courbes imposées** à un bord : ses sommets ajustés pour que la B-spline
  cubique de bord approche une courbe donnée (un cercle, § 9) à tolérance.
- **Cage depuis des sections** : une suite d'anneaux de points (sections d'un
  tube), ouverts ou fermés, donne une cage de quadrangles, fermée ou non le long
  du tube.
- **OBJ** : lecture et écriture d'une cage (`gbs-io`), pour la modeler dans
  Blender et la convertir dans gbs.

---

## 8. Applications

gbs garde les outils génériques : cage, évaluation, conversion, BREP,
interpolation. Les générateurs de pièces qui s'en servent vivent dans des
projets à part, qui jugent chaque PR sur leurs propres géométries en plus des
cages de test.

---

## 9. Limites et prolongements

| Limite | Conséquence | Prolongement possible |
|---|---|---|
| CC est polynomial | un cercle (une bride, un raccord de tube) n'est qu'approché, à tolérance par l'ajustement des bords (§ 7) | subdivision compatible NURBS (Cashman et al. 2009) |
| Limite seulement C1 aux EV | G1 à l'EV, courbure mal tenue autour | constructions de Peters (PCCM 2000) ou de Karčiauskas–Peters, G1 exact, en option |
| Arêtes semi-vives (DeRose et al. 1998) non traitées | arête vive ou lisse, rien entre | netteté entière : s étapes de règles vives puis lisses, exactes en B-splines plus denses près de l'arête |
| Faces irrégulières gardées seules | une petite face par carreau autour de chaque EV | regrouper une bande irrégulière avec sa région voisine, réapprochée |

---

## 10. Validation

| Contrôle | Attendu |
|---|---|
| Région régulière contre l'évaluation exacte, (face, u, v) au hasard | écart relatif < 1e-13 |
| Sommets de la cage raffinée k fois contre les surfaces | sur la surface, aux arrondis près |
| Masques limites contre la valeur des surfaces aux sommets réguliers | idem |
| **Valeurs de référence d'OpenSubdiv** sur des cages types (cube, tore, cage à triangles et pentagones, prisme hexagonal, bords, arêtes vives, coins) : cage raffinée une et deux fois, positions et tangentes limites aux sommets (`Far::PrimvarRefiner::Limit`), limite en (face, u, v) (`Far::PatchTable`, raffinement adaptatif au niveau 10) | mêmes valeurs aux arrondis près là où OpenSubdiv est exact : cages raffinées indice pour indice, positions sauf aux darts, normales aux sommets lisses et aux bords réguliers, surface sur ses carreaux réguliers (hors des carrés de 2⁻¹⁰ autour des sommets irréguliers, qu'il ferme par des carreaux approchés) ; ailleurs, contre la subdivision répétée et les normales voisines |
| Faces irrégulières sur échantillons indépendants | écart ≤ `approx_tol`, angle ≤ `angle_tol` |
| BREP | `check()` sans défaut ; cage fermée : solide, volume qui converge vers celui de la cage raffinée ; aller-retour STEP ; `BRepCheck_Analyzer` d'OpenCascade sur le STEP |

Les valeurs de référence viennent de `scripts/subd_reference/` : un programme
C++ d'une centaine de lignes, hors du build de gbs, qui lit les cages types et
écrit les valeurs en JSON dans `tests/in/subd/`. On ne le relance que si les
cages ou les règles changent ; les tests lisent le JSON, la CI n'a pas besoin
d'OpenSubdiv. Aucun code d'OpenSubdiv n'entre dans gbs, seulement des nombres.

OpenSubdiv est installé dans l'environnement de développement, pas en CI :

```bash
micromamba install -n <env> -c conda-forge libopensubdiv   # 3.7.0, seul paquet ajouté
```

- Bibliothèque compilée (`libosdCPU`, `libosdGPU`) et ses en-têtes
  (`opensubdiv/far`, `sdc`, `bfr`…) ; le programme ne lie que
  `OpenSubdiv::osdCPU`.
- La configuration CMake du paquet conda-forge lie `TBB::tbb` sans le chercher :
  `find_package(TBB REQUIRED)` avant `find_package(OpenSubdiv REQUIRED)`.
- `Far::PrimvarRefiner::Limit` exige `fullTopologyInLastLevel = true` dans les
  options du raffinement uniforme.
- Vérifié sur un cube : position limite d'un coin (−0,5 ; −0,5 ; 0,5), celle du
  masque de Halstead et al. sur la cage raffinée une fois.

---

## 11. Découpage en PR

| PR | Contenu | Note |
|---|---|---|
| 1 | `gbs/subd/` : cage, topologie et refus, règles (lisses, bords, arêtes vives, coins), `refine`, masques limites, `LimitSurface` exacte ; valeurs de référence d'OpenSubdiv (`scripts/subd_reference/`) ; Python ; ce document | `subd_pr01_cage.md` |
| 2 | régions régulières : test de régularité, points fantômes, découpe en rectangles, une B-spline exacte par rectangle, bandes fermées | `subd_pr02_regular.md` |
| 3 | faces irrégulières : nœuds gradués, contraintes, ajustement, raffinement à `approx_tol` et `angle_tol`, mesure | `subd_pr03_irregular.md` |
| 4 | `gbs-brep/subd.h` : coque ou solide, arêtes partagées, rapport ; STEP ; Python ; exemples de la galerie | `subd_pr04_brep.md` |
| 5 | interpolation, bords imposés, cage depuis des sections, OBJ | `subd_pr05_fit.md` |

Les PR 1 à 3 restent dans `gbs/` (en-têtes, sans BREP) ; la PR 4 seule touche
`gbs-brep`.

---

## Références

- E. Catmull, J. Clark, *Recursively generated B-spline surfaces on arbitrary
  topological meshes*, Computer-Aided Design 10(6), 1978.
- M. Halstead, M. Kass, T. DeRose, *Efficient, fair interpolation using
  Catmull–Clark surfaces*, SIGGRAPH 1993.
- J. Stam, *Exact evaluation of Catmull–Clark subdivision surfaces at arbitrary
  parameter values*, SIGGRAPH 1998.
- T. DeRose, M. Kass, T. Truong, *Subdivision surfaces in character animation*,
  SIGGRAPH 1998.
- J. Peters, *Patching Catmull–Clark meshes*, SIGGRAPH 2000.
- J. Peters, U. Reif, *Subdivision Surfaces*, Springer, 2008.
- C. Loop, S. Schaefer, *Approximating Catmull–Clark subdivision surfaces with
  bicubic patches*, ACM TOG 27(1), 2008.
- D. Eppstein, M. Goodrich, E. Kim, R. Tamstorf, *Motorcycle graphs: canonical
  quad mesh partitioning*, SGP 2008.
- T. Cashman, U. Augsdörfer, N. Dodgson, M. Sabin, *NURBS with extraordinary
  points: high-degree, non-uniform, rational subdivision schemes*, SIGGRAPH 2009.
- M. Nießner, C. Loop, M. Meyer, T. DeRose, *Feature-adaptive GPU rendering of
  Catmull–Clark subdivision surfaces*, ACM TOG 31(1), 2012.
- OpenSubdiv, Pixar : schéma *Catmark*, règles de bord et d'arêtes vives.
