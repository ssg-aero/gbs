# Lecture STEP — conception

| | |
|---|---|
| Statut | Proposition, à valider avant tout développement |
| Périmètre | Lecture de fichiers STEP (ISO 10303-21) AP203, AP214 et la partie BREP d'AP242 vers le noyau BREP natif `gbs::brep` |
| Dépend de | palier 1 du BREP natif ([brep_core.md](brep_core.md), notes PR 1 à 11) |
| Hors périmètre | écriture STEP, PMI, couleurs et calques, maillages tessellés AP242, BREP facettés |

Ce document décrit *quoi* construire et *pourquoi*, comme [brep_core.md](brep_core.md)
pour le palier 1. Les pseudo-déclarations C++ sont illustratives.

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Analyse du fichier | **Analyseur Part 21 maison**, header-only, sans dépendance, en C++23 — **décidé** (question 1) | STEPcode (générateur de code, très lourd) ; OCCT (dépendance que l'on veut justement éviter) |
| Emplacement | `gbs-io/step/`, namespace `gbs::step` | module `gbs-step/` séparé |
| Surfaces élémentaires | converties en **NURBS rationnelles exactes** (décision du palier 1) ; domaine borné par les bords de la face — **décidé** (question 3) | classes analytiques |
| Pcurves | **recalculées** depuis les courbes 3D (projection de la PR 5) ; les pcurves du fichier ne servent qu'à départager le côté d'un seam — **décidé** (question 3) | réutiliser les pcurves du fichier (paramétrage incompatible avec les NURBS converties, absentes de nombreux fichiers) |
| Topologie | **reprise telle quelle** du fichier : sommets, arêtes, ordre et sens des co-arêtes ; aucun sewing | recoudre des faces indépendantes (perte de l'information exacte du fichier) |
| Compléments au palier 1 | wire **ordonné** acceptant un seam ; **arêtes dégénérées** insérées aux pôles ; **découpe** d'une arête au seam — **décidé** (question 5) | rejeter ces faces |
| Unités | converties vers une unité cible (défaut **millimètre**) ; unité d'origine conservée dans le modèle | garder l'unité du fichier |
| Assemblages | **aplatis** : chaque instance devient une copie transformée de la géométrie, dans un compound — **décidé** (question 7) | instances partagées avec placement (non prévu au palier 1, question 5) |
| Échecs | **import partiel** avec rapport détaillé ; option stricte | tout ou rien |
| Validation | fichiers STEP écrits à la main dans `tests/in/step/` + **comparaison avec le lecteur STEP d'OCCT** dans `gbs-occt/tests` | fichiers industriels seuls |

---

## 1. Objectifs et non-objectifs

### 1.1 Objectifs

1. Lire un fichier STEP Part 21 et en extraire les solides (`manifold_solid_brep`,
   `brep_with_voids`), les shells ouverts (`shell_based_surface_model`) et les
   faces, dans un `gbs::brep::Model<double>`.
2. Couvrir la géométrie des schémas AP203 (`config_control_design`), AP214
   (`automotive_design`) et la partie BREP d'AP242, tels qu'écrits par les
   systèmes courants (CATIA, NX, Creo, SolidWorks, OCCT).
3. Produire un modèle **valide** au sens de `check()` (PR 6), avec un rapport
   de lecture : entités ignorées, faces non construites, écarts de tolérance.
4. Conserver noms de produits, identifiants STEP (`#id`) et unités.
5. API C++ et Python.

### 1.2 Non-objectifs

- Écriture STEP (voir question 12).
- PMI, couleurs (`styled_item`), calques, propriétés de validation.
- Maillages tessellés (`tessellated_*` d'AP242) et BREP facettés (`faceted_brep`,
  `poly_loop`) : voir question 13.
- Réparation de géométrie au-delà de ce que font les builders (pas de healing,
  décision du palier 1).

---

## 2. Analyse du format Part 21

### 2.1 Ce qu'il faut savoir lire

```
ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''),'2;1');
FILE_NAME('part.stp','2026-10-09T10:00:00',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));
ENDSEC;
DATA;
#10=CARTESIAN_POINT('',(0.,0.,0.));
#11=DIRECTION('',(0.,0.,1.));
#20=(BOUNDED_CURVE()B_SPLINE_CURVE(2,(#30,#31,#32),.UNSPECIFIED.,.F.,.F.)
     B_SPLINE_CURVE_WITH_KNOTS((3,3),(0.,1.),.UNSPECIFIED.)CURVE()
     GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE((1.,0.707106781,1.))
     REPRESENTATION_ITEM(''));
#40=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#41,'distance_accuracy_value','');
ENDSEC;
END-ISO-10303-21;
```

| Élément lexical | Exemple | Remarque |
|---|---|---|
| instance simple | `#12=LINE('',#10,#13);` | |
| instance complexe | `#20=(A()B(...)C(...));` | indispensable : NURBS rationnelles, unités |
| référence | `#13` | résolue après lecture complète (références en avant autorisées) |
| réel | `1.E-07`, `-2.5`, `0.` | `strtod_l` (locale « C ») sous libc++, `std::from_chars` ailleurs (voir ci-dessous) |
| entier | `3` | |
| chaîne | `'it''s'` | `''` échappe `'` ; `\X2\…\X0\` (UTF-16) et `\S\` décodés |
| énumération | `.T.`, `.UNSPECIFIED.` | |
| omis / dérivé | `$`, `*` | |
| liste | `(…)`, imbriquée | |
| paramètre typé | `LENGTH_MEASURE(1.E-07)` | |
| commentaire | `/* … */` | |
| binaire | `"0123"` | rare, conservé brut |

### 2.2 Architecture de l'analyseur

```
gbs-io/step/p21.h
  Lexer      : flux de jetons sur un tampon (fichier en mémoire), sans allocation par jeton
  Parser     : instances → table dense  id → Record { type (ou types d'une instance complexe), Value args }
  Value      : variant { Ref, Real, Int, String, Enum, Omitted, Derived, List, Typed }
  P21File    : header (schéma, nom), table des instances, index par type
```

- **Deux passes** : lecture de toutes les instances (le fichier est en mémoire),
  puis résolution des références à la demande. Les références en avant sont
  normales en Part 21.
- **Arène** : les valeurs vivent dans des `std::vector` contigus (même choix que
  le `Model` du palier 1), une instance pointe sur un intervalle. Objectif de
  performance : lire un fichier de 100 Mo en quelques secondes.
- Les noms d'entités sont normalisés en majuscules ; le type d'une instance
  complexe est l'ensemble trié de ses types partiels.
- Erreur de syntaxe : `std::expected<P21File, P21Error>`, l'erreur portant la
  ligne, la colonne et le jeton attendu ; c'est la seule erreur fatale (un
  fichier mal formé n'a pas de sens partiel). Même convention que les builders
  du palier 1 (`BuildResult`).

**Moyens C++23 retenus** (décision de la question 1) — limités à ce que les
trois chaînes de la CI fournissent (clang ≥ 19, MSVC 2022, AppleClang ≥ 17) :

| Besoin | Moyen |
|---|---|
| tampon du fichier sans copie, jetons | `std::string_view`, `std::span` |
| erreurs sans exception | `std::expected` |
| lecture des entiers | `std::from_chars` |
| lecture des réels, exacte et sans dépendance à la locale | `std::from_chars` avec libstdc++ et la STL de MSVC ; `strtod_l` avec la locale « C » sous libc++, dont la version conda-forge déclare `from_chars` sur les flottants indisponible sous macOS (constaté en PR 1, voir [step_pr01_p21.md](step_pr01_p21.md)) |
| valeurs d'un paramètre | `std::variant` et `std::visit` |
| parcours et filtres de la table d'instances | `std::ranges`, `std::views`, `std::ranges::to` |
| états impossibles du lexer | `std::unreachable` |
| énumérations vers leur valeur | `std::to_underlying` |

Écartés parce que la bibliothèque standard d'Apple ne les fournit pas encore :
`std::flat_map` (on utilise un `std::vector` trié ou une table dense indexée
par `#id`), `std::print`, `std::generator`.

---

## 3. Sous-ensemble du schéma lu

Les entités sont celles de la partie 42 (géométrie et topologie) communes à
AP203, AP214 et AP242.

### 3.1 Contexte, produits, unités

| Entité | Usage |
|---|---|
| `product`, `product_definition_formation`, `product_definition`, `product_definition_shape`, `shape_definition_representation` | trouver la représentation de chaque produit, et son nom |
| `advanced_brep_shape_representation`, `manifold_surface_shape_representation`, `shape_representation` | racines des formes |
| `shape_representation_relationship` | représentation de forme → représentation BREP |
| `geometric_representation_context`, `global_unit_assigned_context`, `global_uncertainty_assigned_context` | unités et tolérance |
| `si_unit` (`length_unit`, `plane_angle_unit`), `conversion_based_unit` | conversion mètre/millimètre/pouce, radian/degré |
| `uncertainty_measure_with_unit` | tolérance de distance du fichier |
| `next_assembly_usage_occurrence`, `context_dependent_shape_representation`, `representation_relationship_with_transformation`, `item_defined_transformation` | assemblages (§ 7) |

L'unité d'**angle** compte autant que l'unité de longueur : elle s'applique à
`conical_surface.semi_angle` et aux paramètres de `trimmed_curve` sur un cercle.

### 3.2 Courbes

| Entité STEP | Représentation gbs |
|---|---|
| `cartesian_point`, `direction`, `vector`, `axis1_placement`, `axis2_placement_3d`, `axis2_placement_2d` | points et repères |
| `line` | segment NURBS de degré 1, borné par l'arête |
| `circle`, `ellipse` | arc ou cercle NURBS rationnel exact (degré 2) |
| `b_spline_curve_with_knots` (+ `rational_b_spline_curve` en instance complexe), `bezier_curve`, `quasi_uniform_curve`, `uniform_curve` | `BSCurve` / `BSCurveRational` : multiplicités et nœuds, **pôles cartésiens + poids convertis en coordonnées homogènes** (leçon de l'IGES) |
| `trimmed_curve` | courbe de base restreinte (paramètres ou points de coupe) |
| `composite_curve`, `composite_curve_segment` | `CurveComposite` ou concaténation NURBS |
| `polyline` | NURBS de degré 1 |
| `surface_curve`, `seam_curve` | la `curve_3d` ; les `pcurve` associées sont lues pour le côté du seam |
| `pcurve`, `definitional_representation` | courbe 2D dans le paramétrage **STEP** de la surface (§ 4.3) |
| `offset_curve_3d`, `hyperbola`, `parabola` | rapport « non supporté » au premier jet, puis approximation |

### 3.3 Surfaces

| Entité STEP | Représentation gbs |
|---|---|
| `plane` | NURBS bilinéaire sur le rectangle englobant les bords de la face (§ 4.2) |
| `cylindrical_surface` | NURBS rationnelle exacte : cercle × segment |
| `conical_surface` | NURBS rationnelle exacte : cercle × segment de rayon variable ; sommet = côté dégénéré |
| `spherical_surface` | NURBS rationnelle exacte : cercle × demi-cercle ; pôles dégénérés |
| `toroidal_surface`, `degenerate_toroidal_surface` | NURBS rationnelle exacte : cercle × cercle |
| `b_spline_surface_with_knots` (+ `rational_b_spline_surface`), `bezier_surface` | `BSSurface` / `BSSurfaceRational`, pôles transposés de l'ordre STEP (`[u][v]`) vers l'ordre gbs (`u` le plus rapide) |
| `surface_of_revolution` | révolution rationnelle exacte de la génératrice NURBS |
| `surface_of_linear_extrusion` | extrusion exacte de la génératrice NURBS |
| `rectangular_trimmed_surface` | surface de base restreinte |
| `offset_surface` | `SurfaceOffset` existant, ou rapport « non supporté » |

**Constructions à ajouter à gbs** (aujourd'hui seuls le cercle complet,
l'ellipse et le segment existent) : arc de cercle et d'ellipse NURBS exacts entre
deux angles, révolution rationnelle exacte d'une NURBS autour d'un axe,
extrusion d'une NURBS. Ce sont des constructions classiques (*The NURBS Book*,
§ 7.3 à 8.5) qui serviront aussi au palier 2 (sweep, révolution en BREP).

### 3.4 Topologie

| Entité STEP | `gbs::brep` |
|---|---|
| `vertex_point` | `Vertex` (tolérance = incertitude du fichier) |
| `edge_curve(start, end, geometry, same_sense)` | `Edge` ; si `same_sense = .F.`, la courbe est inversée pour garder `u1 < u2` de `v1` vers `v2` |
| `oriented_edge(edge, orientation)` | `CoEdge` (sens) |
| `edge_loop` | `Wire` **ordonné**, dans l'ordre du fichier |
| `vertex_loop` | arête dégénérée (pôle) |
| `face_outer_bound`, `face_bound` (+ `orientation`) | contour extérieur et trous ; une `orientation = .F.` inverse la boucle |
| `advanced_face(bounds, surface, same_sense)`, `face_surface` | `Face` ; `same_sense = .F.` ⇒ `FaceUse` `Reversed` dans le shell (§ 5.3) |
| `closed_shell`, `open_shell`, `oriented_closed_shell` | `Shell` |
| `manifold_solid_brep` | `Solid` |
| `brep_with_voids` | `Solid` avec cavités |
| `shell_based_surface_model` | shells ouverts dans un `Compound` |

---

## 4. Géométrie : paramétrage et pcurves

### 4.1 Pourquoi convertir change le paramétrage

Une `cylindrical_surface` STEP est paramétrée par l'angle `u` et la hauteur `v`.
Son équivalent NURBS rationnel exact a la même géométrie mais un paramètre `u`
qui n'est **pas** l'angle (le cercle rationnel en quatre quarts n'est pas
paramétré par l'angle). Il en va de même pour le cône, la sphère et le tore. Les
pcurves du fichier, exprimées dans le paramétrage STEP, ne sont donc pas
directement valables sur la surface convertie ; beaucoup de fichiers n'en
contiennent d'ailleurs pas.

### 4.2 Stratégie retenue

1. Les **courbes 3D** des arêtes sont lues exactement.
2. Les **surfaces** sont converties en NURBS exactes. Les surfaces infinies
   (plan, cylindre, cône) sont bornées par le rectangle paramétrique qui
   englobe les bords de la face, avec une marge, puis converties.
3. Les **pcurves sont recalculées** par la projection de la PR 5
   (`extract_pcurve` / `project_pcurve`), qui gère déjà la continuité au seam
   et aux pôles.
4. Les pcurves du fichier servent uniquement à **départager le côté d'un seam**
   quand l'ambiguïté subsiste (arête entièrement sur le seam en début de boucle).

Alternative écartée : transformer exactement les pcurves STEP par la
fonction angle → paramètre rationnel (`t = tan(θ/2)` par quart). Elle est exacte
mais non polynomiale : il faudrait de toute façon réinterpoler, avec plus de
code que la projection déjà écrite et testée.

### 4.3 Arêtes : bornes et sens

- La courbe d'une `edge_curve` est souvent infinie (`line`) ou fermée (`circle`) :
  les bornes `[u1, u2]` sont obtenues en **projetant les deux sommets** sur la
  courbe convertie. Pour une arête fermée (même sommet aux deux bouts), la
  plage est la période entière.
- `same_sense = .F.` : la courbe convertie est inversée (`reverse()` d'une
  NURBS), pour respecter la convention « `u` croissant de `v1` vers `v2` » du
  palier 1.

---

## 5. Topologie : compléments nécessaires au palier 1

La topologie STEP est **explicite** : les arêtes sont partagées par identité,
l'ordre et le sens des co-arêtes sont donnés. On la reprend telle quelle, sans
sewing. Trois situations courantes ne passent pas les builders du palier 1.

### 5.1 Wire ordonné avec seam

`make_wire` (PR 3) réordonne les arêtes et refuse une arête répétée. Une face
cylindrique complète STEP a pour boucle `[cercle bas, seam, cercle haut inversé,
seam inversé]` : la même arête deux fois. Il faut un builder
`make_wire_ordered(m, {(EdgeId, Orientation)…})` qui vérifie le chaînage dans
l'ordre donné et accepte une arête utilisée deux fois en sens opposés. La
projection de la PR 5 donne alors les deux pcurves du seam (côtés `u = u1` et
`u = u2`) grâce au repère de continuité déjà en place.

### 5.2 Pôles : arêtes dégénérées insérées

STEP n'a pas d'arête dégénérée. Une sphère complète est souvent bornée par
`[seam, seam inversé]` : en 3D la boucle est fermée (les deux bouts sont les
pôles), mais en `(u, v)` il manque les côtés `v = ±π/2`. Le builder de face
insère une **arête dégénérée** (PR 3, `make_degenerate_edge`) partout où deux
co-arêtes consécutives se touchent en 3D sur un point singulier de la surface
mais pas en `(u, v)`. Même chose au sommet d'un cône, ou pour un `vertex_loop`.

### 5.3 Sens des faces

STEP oriente une boucle extérieure dans le sens direct **vu depuis la normale
de la face**, qui est opposée à celle de la surface quand
`advanced_face.same_sense = .F.`. Le builder de la PR 5 tourne toujours le
contour extérieur dans le sens direct en `(u, v)` en retournant ses co-arêtes ;
pour garder les sens effectifs du fichier, on pose alors `FaceUse = Reversed`
dans le shell. Les deux retournements se compensent : les usages des arêtes
restent opposés entre faces voisines, comme dans le fichier. `make_solid` (PR 8)
vérifie ensuite le sens extérieur par le volume.

### 5.4 Arêtes qui traversent le seam

Une arête peut traverser le seam d'une surface périodique (un arc de cercle
de 350° à 10° sur un cylindre dont le seam est à 0°). La PR 5 la refuse
(`CrossesSeam`). Le lecteur la **découpe** au point où elle atteint le seam :
recherche du paramètre de l'arête où `u(t)` de la projection passe par la
borne, nouveau sommet, deux arêtes. C'est une opération de type palier 2
(découpe d'arête), limitée ici au seam ; elle est isolée dans une PR (question 5).

---

## 6. Unités, tolérances, noms

### 6.1 Attributs du modèle (prévus au palier 1, § 9.2)

Le palier 1 avait réservé ces tables sans les implémenter. Elles arrivent avec
le lecteur :

```cpp
template <std::floating_point T> class Model {
  …
  void setName(const ShapeId &, std::string);      auto name(const ShapeId &) const -> std::string_view;
  void setExternalId(const ShapeId &, std::int64_t); auto externalId(const ShapeId &) const -> std::optional<std::int64_t>;
  T unit_scale{1};   // facteur de l'unité du modèle vers le millimètre ; 1 si inconnu
};
```

Tables creuses (`unordered_map`), vides pour un modèle natif ; `compact()` et
`append()` les remappent.

### 6.2 Unités

Toutes les longueurs sont converties vers l'unité cible des options (défaut :
millimètre) au moment de la conversion géométrique ; les angles vers le radian.
L'unité du fichier est rapportée et `unit_scale` renseigné.

### 6.3 Tolérances

L'incertitude du fichier (`uncertainty_measure_with_unit`, souvent 1e-7 à
1e-3 mm), convertie, devient la tolérance par défaut des sommets et des arêtes ;
les builders la remontent ensuite si la géométrie réelle s'en écarte (écarts
entre courbes et surfaces des fichiers industriels). Ces remontées sont
rapportées : une tolérance très supérieure à l'incertitude annoncée signale un
fichier de mauvaise qualité.

### 6.4 Noms et identifiants

Nom de produit → nom du solide ou du compound ; nom non vide d'une face ou d'une
arête → son nom ; `#id` de chaque entité topologique lue → `externalId`.

---

## 7. Assemblages

Un fichier d'assemblage décrit des produits, leurs représentations, et des
occurrences avec placement (`item_defined_transformation` entre deux
`axis2_placement_3d`). Conformément à la question 5 du palier 1 (pas de
placement dans le modèle), le lecteur **aplatit** l'arbre : chaque occurrence
devient une copie de la forme du composant, géométrie transformée
(`transform` de gbs sur les NURBS), dans un compound par sous-assemblage, nommé
d'après le produit. La géométrie est dupliquée par instance : acceptable pour le
palier ; un partage par instance demanderait des placements dans le modèle
(question 7).

---

## 8. Rapport et politique d'échec

```cpp
struct StepReadReport {
  std::string schema;                               // AUTOMOTIVE_DESIGN, CONFIG_CONTROL_DESIGN, AP242…
  std::string length_unit; double unit_to_target;   // unité du fichier et facteur appliqué
  double file_tolerance;
  std::vector<std::pair<std::int64_t, std::string>> unsupported;   // #id, type
  std::vector<std::pair<std::int64_t, BuildError>> failed_faces;    // #id de l'advanced_face, cause
  std::size_t solids, shells, faces, edges;
  double max_tolerance_raised;
};
```

- **Partiel par défaut** : une face qui ne peut pas être construite est
  omise, le shell reste ouvert, le solide n'est pas créé et le shell est rendu
  seul ; tout est dans le rapport.
- **Strict** (`StepReadOptions::strict`) : la première face en échec fait
  échouer la lecture (`BuildResult` en erreur).
- Le modèle produit passe toujours `check()` ; le rapport de `check()` est joint.

---

## 9. API

```cpp
namespace gbs::step {
  struct StepReadOptions { double target_unit_mm = 1.; bool strict = false; bool assemblies = true; brep::MakeFaceOptions<double> face{}; };
  struct StepReadResult  { brep::ShapeId root; StepReadReport report; };

  auto read_step(brep::Model<double> &m, const std::filesystem::path &file, StepReadOptions = {}) -> brep::BuildResult<StepReadResult>;
  auto read_step(brep::Model<double> &m, std::string_view content, StepReadOptions = {})          -> brep::BuildResult<StepReadResult>;

  // niveau bas, utilisable seul (inspection, outils)
  auto parse_p21(std::string_view content) -> std::expected<P21File, P21Error>;
}
```

Python : `gbs.brep.read_step(model, path, target_unit_mm=1., strict=False)` →
`(root, report)`.

Header-only comme le reste de gbs, **sans nouvelle dépendance** : le lecteur
s'appuie sur `gbs-brep` et la géométrie gbs.

---

## 10. Stratégie de test

1. **Fichiers écrits à la main** dans `tests/in/step/`, petits et lisibles : boîte,
   cylindre (seam), sphère (pôles), plaque trouée, carreau B-spline rationnel,
   cône, tore, assemblage de deux boîtes, unités en pouces et en degrés. Ils
   tournent sur toutes les plateformes.
2. **Comparaison avec OCCT** dans `gbs-occt/tests` (job OCCT de la CI) : formes
   primitives OCCT (`BRepPrimAPI_*`, booléens, congés) écrites par
   `STEPControl_Writer`, relues par gbs et par `STEPControl_Reader` ; on compare
   nombres de faces, arêtes, sommets, volumes et surfaces, et `check()` doit
   être valide. C'est la validation la plus forte : des fichiers réalistes,
   produits par un écrivain de référence, sans données à versionner.
3. **Fichiers industriels publics** (jeux de test CAx-IF / NIST, publiés pour
   ces essais) : test optionnel, désactivé par défaut, sur un répertoire local.
4. **Robustesse de l'analyseur** : jetons limites, chaînes encodées, instances
   complexes, commentaires, fichier tronqué ⇒ erreur propre.

---

## 11. Plan de développement en PR

Chaque PR est accompagnée de sa note d'architecture (règle du palier 1).

| # | PR | Contenu | Tests | Lignes | h |
|---|---|---|---|---|---|
| 1 | `step/p21` | Lexer et analyseur Part 21, `P21File`, valeurs, instances complexes, décodage des chaînes | extraits écrits à la main, erreurs avec position, fichier de 50 Mo généré (temps) | 550 | 5 |
| 2 | `brep/attributes` | Noms, identifiants externes, `unit_scale` dans `Model` (palier 1 § 9.2), remappés par `compact` / `append` | ajout, suppression, compaction, fusion de modèles | 300 | 3 |
| 3 | `geom/arcs-revolution` | Arc de cercle et d'ellipse NURBS exacts, révolution rationnelle et extrusion d'une NURBS | exactitude des points, rayons, poids | 450 | 4 |
| 4 | `step/geometry` | Repères, courbes, surfaces élémentaires et B-spline → NURBS, unités (longueur, angle), bornes des surfaces infinies | chaque entité contre sa formule analytique, rationnelles comprises | 600 | 6 |
| 5 | `brep/ordered-wire` | `make_wire_ordered` (seam répété), insertion d'arêtes dégénérées aux pôles, sens des faces | cylindre, sphère, cône depuis des boucles ordonnées ; `check()` | 450 | 4 |
| 6 | `step/topology` | Sommets, arêtes (bornes par projection, `same_sense`), boucles, faces, shells, solides, cavités, surfaces ouvertes ; rapport ; mode strict | fichiers `tests/in/step/` : boîte, cylindre, sphère, plaque trouée, carreau B-spline ; `check()` valide | 600 | 6 |
| 7 | `step/seam-split` | Découpe des arêtes qui traversent un seam | arc à cheval sur le seam d'un cylindre ; `check()` | 350 | 4 |
| 8 | `step/occt-compare` | Comparaison systématique avec OCCT (primitives, booléens, congés écrits par OCCT) | volumes, comptes, validité | 400 | 4 |
| 9 | `step/assemblies` | Produits, occurrences, transformations, aplatissement en compounds nommés | assemblage écrit à la main et par OCCT | 450 | 5 |
| 10 | `step/python` | `gbs.brep.read_step`, rapport, documentation, `News.md` | pytest | 250 | 3 |

Total : environ 4 400 lignes, **44 h**. Ordre : 1, 2 et 3 sont indépendantes ;
4 dépend de 3 ; 5 du palier 1 seul ; 6 de 1, 2, 4, 5 ; 7 et 8 de 6 ; 9 de 6 ;
10 en dernier.

---

## 12. Questions ouvertes

Pour chaque question : **R** = recommandation, **A** = alternatives.

**État au 9 octobre 2026** : questions 1, 3, 5, 7 et 12 décidées (voir chaque question) ; questions 2, 4, 6, 8, 9, 10, 11, 13 et 14 encore ouvertes.

1. **Analyseur Part 21 maison ?** — **Décidé : oui, écrit en C++23** (§ 2.2).
   R : oui, header-only, environ 550 lignes ; le format est simple et stable.
   A : STEPcode (bibliothèque générée depuis les schémas EXPRESS, lourde à
   construire et à empaqueter) ; OCCT (contraire à l'objectif d'indépendance).

2. **Emplacement : `gbs-io/step/`, namespace `gbs::step` ?**
   R : oui, à côté de l'IGES. A : module `gbs-step/` séparé.

3. **Surfaces élémentaires en NURBS et pcurves recalculées ?** — **Décidé : conversion en NURBS rationnelles exactes, pcurves recalculées.**
   R : oui (§ 4.2), cohérent avec la question 6 du palier 1. A : classes
   analytiques (`Plane`, `Cylinder`…) qui garderaient les pcurves STEP telles
   quelles et faciliteraient les intersections analytiques du palier 2.

4. **Bornes des surfaces infinies** : rectangle englobant les bords de la face,
   avec une marge de 10 % ? R : oui. A : surface commune bornée par toutes les
   faces qui la partagent (une seule surface pour plusieurs faces).

5. **Arêtes traversant un seam : découpe dans le lecteur (PR 7) ?** — **Décidé : oui, découpe dans le lecteur.**
   R : oui, sinon une part notable des fichiers industriels échoue. A : refuser
   la face et la rapporter ; ou déplacer le seam de la surface convertie pour
   que l'arête ne le traverse plus (possible pour une seule arête, pas en
   général).

6. **Unité cible millimètre par défaut, unité d'origine conservée ?**
   R : oui. A : garder l'unité du fichier telle quelle (gbs n'impose pas
   d'unité, question 11 du palier 1).

7. **Assemblages aplatis (géométrie dupliquée par instance) ?** — **Décidé : oui, aplatis pour cette phase.**
   R : oui pour cette phase. A : introduire des placements dans le modèle
   (`FaceUse` / `Compound` avec transformation), plus économe pour les
   assemblages à nombreuses instances.

8. **Import partiel avec rapport par défaut, option stricte ?**
   R : oui. A : tout ou rien.

9. **AP242 : lire la partie BREP et ignorer le reste ?**
   R : oui, la géométrie et la topologie sont celles de la partie 42.
   A : refuser les fichiers AP242.

10. **Attributs du modèle (noms, identifiants, unité) dans `gbs-brep` (PR 2) ?**
    R : oui, comme prévu au palier 1. A : les garder dans une structure à côté
    du modèle, propre au lecteur.

11. **Couleurs et calques ?** R : hors périmètre de cette phase.
    A : lire `styled_item` / `presentation_layer_assignment` dans une table
    d'attributs.

12. **Écriture STEP ?** — **Décidé : phase suivante.**
    R : phase suivante ; l'analyseur Part 21 donne le
    format et les tables d'entités à réutiliser. A : l'ajouter dès maintenant.

13. **BREP facettés (`faceted_brep`, `poly_loop`) ?** R : plus tard ; ils se
    lisent simplement (faces planes polygonales) mais concernent surtout les
    fichiers de maillage. A : les inclure dans la PR 6.

14. **Fichiers de test industriels versionnés ?** R : non ; fichiers écrits à la
    main et fichiers produits par OCCT dans les tests, fichiers CAx-IF en test
    optionnel local. A : versionner quelques fichiers publics de petite taille.
