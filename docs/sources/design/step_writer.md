# Écriture STEP — conception

| | |
|---|---|
| Statut | Proposition, à valider avant tout développement ; questions 1 et 6 décidées le 12 octobre 2026 (§ 11), deux PR préliminaires ajoutées (§ 12) ; **représentation analytique retenue** pour les primitives (§ 12.1), ce qui révise la question 6 du palier 1 et la question 3 de la lecture |
| Périmètre | Écriture de fichiers STEP (ISO 10303-21) AP214, et AP242 en option, depuis le noyau BREP natif `gbs::brep` |
| Dépend de | lecture STEP réalisée ([step_reader.md](step_reader.md), PR 127 à 137) ; export IGES du palier 1 ([brep_pr09_iges.md](brep_pr09_iges.md)) |
| Hors périmètre | PMI, couleurs et calques, maillages tessellés, BREP facettés, AP203 |

Ce document décrit *quoi* construire et *pourquoi*, comme
[step_reader.md](step_reader.md) pour la lecture. Les pseudo-déclarations C++
sont illustratives.

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Sérialisation | **Écrivain Part 21 maison**, pendant de l'analyseur (`p21.h`), sans dépendance | réutiliser une bibliothèque tierce |
| Emplacement | `gbs-io/step/`, namespace `gbs::step`, à côté du lecteur | module séparé |
| Surfaces | **primitives analytiques de gbs** (plan, cylindre, cône, sphère, tore ; droite, cercle, ellipse), présentes dans le BREP à côté des NURBS et écrites en entités analytiques ; le reste en **B-spline** — **décidé** (question 1, PR préliminaire P1) | tout en B-spline ; reconnaître les primitives dans les NURBS ; NURBS qui portent leur définition analytique |
| Courbes des arêtes | B-spline **restreinte à la plage de l'arête** (insertion de nœuds, exact) | la courbe entière, la plage étant déduite des sommets |
| Pcurves | **écrites** (`surface_curve` / `seam_curve` avec leurs `pcurve`), valables telles quelles : le paramétrage des primitives est celui de STEP | omises (STEP les rend facultatives) |
| Géométrie non NURBS | **approximée** dans `approx_tol` et rapportée, comme à l'export IGES ; code d'approximation **mis en commun** | refuser la forme |
| Sens des faces | `advanced_face.same_sense` = sens de l'usage de la face ; boucles écrites telles quelles, borne d'orientation `.F.` pour une face retournée (comme OCCT) | retourner les boucles ; `oriented_face` |
| Arêtes dégénérées | **omises** (STEP n'en a pas ; le lecteur les réinsère) | `vertex_loop` |
| Produits | **une pièce par solide ou shell de premier niveau** ; un **compound devient un assemblage** nommé ; les **instances** d'une même pièce sont écrites **une fois, avec leurs placements** — **décidé** (question 6, PR préliminaire P2) | une seule pièce contenant toutes les formes ; placements identité seulement |
| Schéma | **AP214** par défaut, AP242 en option | AP203 (entités d'approbation obligatoires) |
| Unités | celles du modèle (`unitScale`) : millimètre, mètre, pouce, ou unité par conversion | toujours le millimètre |
| Échecs | écriture **tout ou rien** pour les erreurs ; formes non écrivables (arêtes et sommets libres) **rapportées et omises** | écriture partielle silencieuse |
| Validation | **aller-retour** gbs → STEP → gbs, et **lecture par OCCT** (validité, volumes, comptes) ; double aller-retour OCCT → gbs → STEP → OCCT sur les 20 formes de la PR 8 de la lecture | relecture par gbs seule |

---

## 1. Objectifs et non-objectifs

### 1.1 Objectifs

1. Écrire les solides (avec cavités), les shells ouverts ou fermés et les
   faces d'un `gbs::brep::Model<double>` dans un fichier STEP lisible par les
   systèmes courants (CATIA, NX, Creo, SolidWorks, OCCT).
2. Géométrie exacte pour tout ce qui est NURBS, c'est-à-dire tout ce que
   produisent le lecteur STEP et les builders ; approximation contrôlée pour le
   reste.
3. Conserver les noms (produits, faces), la structure des compounds et l'unité
   du modèle.
4. Relire à l'identique : un aller-retour gbs → STEP → gbs redonne les mêmes
   nombres d'entités, le même volume, et un `check()` valide.
5. API C++ et Python.

### 1.2 Non-objectifs

- PMI, couleurs (`styled_item`), calques.
- Placements dans le modèle : la géométrie reste dupliquée par instance
  (question 7 de la lecture). Le partage est restitué **à l'écriture**, grâce à
  la provenance des instances (§ 12.2).
- AP203, qui impose des entités d'approbation, de sécurité et de personnes.

---

## 2. Écrivain Part 21

### 2.1 Ce qu'il doit produire

```
ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''),'2;1');
FILE_NAME('part.stp','2026-10-12T10:00:00',(''),(''),'gbs','gbs','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=(BOUNDED_CURVE()B_SPLINE_CURVE(2,(#3,#4,#5),.UNSPECIFIED.,.F.,.F.)
   B_SPLINE_CURVE_WITH_KNOTS((3,3),(0.,1.),.UNSPECIFIED.)CURVE()
   GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE((1.,0.707106781186548,1.))
   REPRESENTATION_ITEM(''));
ENDSEC;
END-ISO-10303-21;
```

### 2.2 Architecture

```
gbs-io/step/p21_writer.h
  P21Writer
    add(type, args…) -> #id           instance simple
    add_complex({(type, args…)…}) -> #id  instance complexe, parties dans l'ordre alphabétique (règle Part 21)
    valeurs : Ref, Int, Real, String, Enum, Logical, Omitted ($), Derived (*), List, Typed
    header(description, name, author, organization, schema)
    str() / write(path)
```

- **Numérotation** dans l'ordre d'ajout. Les références vers l'avant ne sont
  pas nécessaires : on écrit les feuilles d'abord.
- **Réels** : écrits au plus court qui se relise à l'identique, avec toujours
  un point décimal (`1.`, `1.5E-07`), indépendamment de la locale.
  `std::to_chars` sur les flottants est, comme `from_chars`, indisponible avec
  la libc++ de conda-forge sous macOS. On reprend donc la solution de la PR 1
  de la lecture : `snprintf_l` avec la locale « C » sous libc++, `to_chars`
  ailleurs.
- **Chaînes** : apostrophe doublée, antislash doublé ; tout caractère hors ASCII
  imprimable encodé `\X2\…\X0\` (UTF-16), `\X4\` au-delà du plan de base.
  C'est l'inverse exact du décodage de la PR 1.
- **Lignes** coupées pour rester lisibles. Part 21 n'impose pas de longueur
  maximale, mais certains lecteurs anciens en ont une.
- **Contrôle** : tout fichier écrit est relu par `parse_p21` dans les tests.

---

## 3. Géométrie

### 3.1 Courbes et surfaces

| gbs | STEP écrit |
|---|---|
| primitive de gbs (§ 12.1) : droite, cercle, ellipse | `line`, `circle`, `ellipse` |
| primitive de gbs : plan, cylindre, cône, sphère, tore | `plane`, `cylindrical_surface`, `conical_surface`, `spherical_surface`, `toroidal_surface` |
| `BSCurve<3>` | `b_spline_curve_with_knots` (multiplicités et nœuds distincts) |
| `BSCurveRational<3>` | instance complexe `bounded_curve`, `b_spline_curve`, `b_spline_curve_with_knots`, `curve`, `geometric_representation_item`, `rational_b_spline_curve`, `representation_item` ; **pôles cartésiens et poids**, convertis depuis la forme homogène de gbs (leçon de l'IGES) |
| `BSSurface<3>` / `BSSurfaceRational<3>` | `b_spline_surface_with_knots` (complexe si rationnelle) ; pôles **transposés** de l'ordre gbs (u le plus rapide) vers l'ordre STEP `[u][v]` |
| pcurve `BSCurve<2>` / `BSCurveRational<2>` | les mêmes entités en 2D, dans une `definitional_representation` de contexte 2D |
| autre courbe ou surface (`CurveOnSurface`, `SurfaceOfRevolution`, décalages…) | approximée (§ 3.3) |

### 3.2 Une courbe par arête, restreinte à l'arête

L'arête gbs porte une courbe et une plage [u1, u2] qui peut n'en couvrir qu'une
partie. C'est le cas d'un cercle coupé par `split_edge`, ou d'une courbe
partagée par deux arêtes. Un lecteur STEP déduit la plage en projetant les
sommets ; sur une courbe fermée, c'est ambigu. L'écrivain écrit donc la courbe
**restreinte à [u1, u2]**, ce qui est exact pour une B-spline (insertion de
nœuds), avec `same_sense = .T.`. Une arête fermée devient une courbe fermée
entière.

### 3.3 Approximation (mise en commun avec l'IGES)

L'export IGES (palier 1, PR 9) convertit déjà en NURBS cubique, par
interpolation à paramètres conservés, toute courbe ou surface non NURBS,
jusqu'à `approx_tol`. Ce code (`iges_detail::to_nurbs`) vit dans un en-tête qui
dépend de libIGES. Il est **déplacé** dans `gbs-io/to_nurbs.h`, sans
dépendance, et partagé par les deux exports (question 11).

Cas particulier : la `CurveOnSurface` d'un côté de face naturelle, dont la
pcurve est un segment isoparamétrique. Sa courbe 3D est **exactement** l'iso de
la NURBS (`isoU` / `isoV`) : aucune approximation.

### 3.4 Pcurves

Le modèle porte une pcurve par co-arête, paramétrée comme l'arête, dans
l'espace de la NURBS de sa face.

- **Surface B-spline** : c'est exactement l'espace paramétrique de la surface
  écrite ; les pcurves sont **valables telles quelles**.
- **Surface primitive** (§ 12.1) : elle a le paramétrage de STEP (l'angle
  pour un cylindre) ; ses pcurves sont, elles aussi, **valables telles
  quelles**. Le pcurve d'un cercle sur un cylindre est un segment, écrit
  exactement en `line` 2D.

Aucune conversion de paramétrage n'est donc nécessaire à l'écriture.

- Une arête utilisée par deux faces : `surface_curve(courbe_3d, (pcurve_1,
  pcurve_2), .CURVE_3D.)`.
- Une arête de seam, deux fois dans la même face : `seam_curve`, avec ses deux
  pcurves sur la même surface.
- Les pcurves sont restreintes à [u1, u2] comme la courbe 3D.

---

## 4. Topologie

| gbs | STEP |
|---|---|
| `Vertex` | `vertex_point(cartesian_point)` |
| `Edge` (non dégénérée) | `edge_curve(v1, v2, courbe, .T.)` (§ 3.2) |
| `Edge` dégénérée | **omise** des boucles (STEP n'a pas d'arête dégénérée ; le lecteur gbs, comme OCCT, la réinsère au pôle) |
| `CoEdge` | `oriented_edge(*, *, edge_curve, orientation)` |
| `Wire` d'une face | `edge_loop` dans l'ordre des co-arêtes |
| `Face` | `advanced_face(bornes, surface, same_sense)` ; `wires[0]` en `face_outer_bound`, les trous en `face_bound` |
| `Shell` | `closed_shell` si fermé, `open_shell` sinon |
| `Solid` | `manifold_solid_brep`, ou `brep_with_voids` s'il a des cavités (`oriented_closed_shell`) |

### 4.1 Sens des faces

Dans le modèle, la boucle extérieure est dans le sens direct **en (u, v)**,
autour de la normale de la surface, et le `FaceUse` du shell dit si la normale
de la face suit celle de la surface. STEP veut la boucle dans le sens direct
**vu depuis la normale de la face**. L'écrivain écrit donc :

- `same_sense = (usage == Forward)` ;
- les boucles telles quelles, avec une borne d'orientation `.T.` si
  `same_sense`, `.F.` sinon : la boucle est alors parcourue à rebours, ce qui
  la remet dans le sens direct vu de la face.

C'est l'écriture d'OCCT, que la PR 8 de la lecture a vue fonctionner dans les
deux sens. Une face utilisée par deux shells en sens contraires (cas non
manifold) est écrite une fois et utilisée par `oriented_face` dans le second.

### 4.2 Formes de premier niveau

| Forme passée à l'écrivain | Écrite comme |
|---|---|
| `Solid` | pièce : `advanced_brep_shape_representation` contenant un `manifold_solid_brep` |
| `Shell` | pièce : `manifold_surface_shape_representation` contenant un `shell_based_surface_model` |
| `Face` | idem, dans un `open_shell` d'une face |
| `Compound` | assemblage (§ 5) |
| `Wire`, `Edge`, `Vertex` libres | **omis et rapportés** (question 8) |

---

## 5. Produits et assemblages

- Chaque pièce est un produit complet : `product`,
  `product_definition_formation`, `product_definition`,
  `product_definition_shape`, `shape_definition_representation`, et le
  contexte d'application du schéma. Son **nom** est celui de la forme, sinon
  un nom des options, sinon `"part"` numéroté.
- Un **compound devient un assemblage** : un produit nommé comme le compound,
  dont chaque membre est un composant (`next_assembly_usage_occurrence`). Les
  compounds imbriqués donnent des sous-assemblages.
- **Instances** (§ 12.2) : les membres qui sont des instances d'une même pièce
  renvoient à **une seule pièce**, écrite une fois dans son propre repère (la
  géométrie de l'instance ramenée par l'inverse de son placement). Chaque
  occurrence porte son placement (`item_defined_transformation`). Un membre
  sans provenance est une pièce à lui seul, placée par l'identité.
- La lecture (PR 9) aplatit un assemblage en compounds nommés et, avec P2,
  garde la provenance de chaque copie : l'écriture **restitue la hiérarchie,
  les noms, le partage des pièces et les placements**.
- Option `flatten` : une seule pièce contenant toutes les formes, pour les
  logiciels qui n'attendent pas d'assemblage.

---

## 6. Unités, tolérance, noms

- **Longueur** : unité du modèle (`Model::unitScale`, en mm par unité) :
  - 1 → `si_unit(.MILLI.,.METRE.)` ;
  - 1000 → `si_unit($,.METRE.)` ;
  - 25,4 → `conversion_based_unit('INCH', 25.4 mm)` ;
  - autre valeur → unité par conversion vers le millimètre.

  Les coordonnées sont écrites telles quelles, dans l'unité du modèle.
- **Angle** : radian (seuls les paramètres angulaires des entités analytiques en
  dépendraient ; il n'y en a pas avec des B-splines seules).
- **Incertitude** (`uncertainty_measure_with_unit`) : `brep_default_tolerance`
  par défaut, modifiable par les options (question 9).
- **Noms** : nom de la forme pour le produit, nom de la face pour
  `advanced_face`, chaîne vide sinon. Les `externalId` ne sont pas réutilisés
  comme `#id` : la numérotation est celle du nouveau fichier.

---

## 7. Rapport et politique d'échec

```cpp
struct StepWriteReport {
  std::size_t entities;                       // instances écrites
  std::size_t solids, shells, faces, edges, vertices;
  std::size_t approximated_curves, approximated_surfaces;
  double max_deviation;                       // des approximations
  std::vector<brep::ShapeId> skipped;         // formes non écrivables (arêtes, sommets libres)
};
```

- **Tout ou rien** pour les erreurs : forme morte, co-arête sans pcurve, face
  sans surface, fichier non inscriptible. La fonction échoue (`BuildError`) et
  rien n'est écrit.
- Les formes qui n'ont pas d'équivalent BREP STEP (wires, arêtes et sommets
  libres) sont omises et listées dans `skipped`.

---

## 8. API

```cpp
namespace gbs::step {
  struct StepWriteOptions {
    std::string schema = "AP214";           // ou "AP242"
    double approx_tol = 1e-5;               // pour la géométrie non NURBS
    bool pcurves = true;
    bool flatten = false;                   // une seule pièce au lieu d'un assemblage
    std::optional<double> uncertainty;      // défaut : brep_default_tolerance
    std::string product_name, author, organization;
  };

  auto write_step(const brep::Model<double> &m, std::span<const brep::ShapeId> shapes, const StepWriteOptions & = {})
      -> brep::BuildResult<std::pair<std::string, StepWriteReport>>;      // contenu du fichier
  auto write_step_file(const brep::Model<double> &m, std::span<const brep::ShapeId> shapes,
                       const std::filesystem::path &file, const StepWriteOptions & = {})
      -> brep::BuildResult<StepWriteReport>;
}
```

Python : `gbs.brep.write_step(model, shapes, path, schema="AP214", …)` →
`report`, `write_step_text(…)` → `(texte, report)`, symétriques de
`read_step` / `read_step_text`.

Header-only, **sans nouvelle dépendance**.

---

## 9. Stratégie de test

1. **Écrivain Part 21** : valeurs de tous les genres, chaînes non ASCII, réels
   limites (`1e-300`, `-0.`, entiers en réels), instances complexes ; chaque
   texte écrit est relu par `parse_p21` et comparé valeur à valeur.
2. **Aller-retour gbs** : les formes des tests BREP et tous les fichiers de
   `tests/in/step` sont lus, écrits, relus. On compare les nombres de faces,
   d'arêtes non dégénérées et de sommets, les volumes et les noms, et le
   `check()` du résultat doit être valide.
3. **Lecture par OCCT** (job OCCT de la CI) : un fichier écrit par gbs est lu
   par `STEPControl_Reader`. `BRepCheck_Analyzer` doit le trouver valide, avec
   les mêmes volumes et les mêmes nombres d'entités.
4. **Double aller-retour** : les 20 formes de la PR 8 de la lecture et
   l'assemblage OCCT suivent le parcours OCCT → STEP → gbs → STEP → OCCT, sans
   écart de volume ni d'entités. C'est le test le plus exigeant : il valide les
   conventions de sens, les seams et les pôles dans le sens de l'écriture.

---

## 10. Plan de développement en PR

Chaque PR est accompagnée de sa note d'architecture.

| # | PR | Contenu | Tests | Lignes | h |
|---|---|---|---|---|---|
| P1 | `geom/primitives` | primitives analytiques de gbs (§ 12.1) : surfaces, ellipse, domaine périodique des coniques, inverse en formule fermée, conversion NURBS exacte ; `copy_shape` ; le lecteur STEP les garde au lieu de convertir | formules, dérivées, inverse, conversion exacte ; faces, solides et `check()` sur primitives ; tests de lecture inchangés en volumes | 700 | 9 |
| P2 | `brep/instances` | provenance des instances dans le modèle (§ 12.2) ; remplie par le lecteur STEP et par `copy_shape` | `erase` / `compact` / `append` ; lecture d'un assemblage | 300 | 4 |
| 1 | `step/p21-writer` | `P21Writer` : valeurs, instances simples et complexes, réels relisibles, chaînes encodées, en-tête et schémas | relecture par `parse_p21` | 300 | 3 |
| 2 | `step/write-geometry` | `gbs-io/to_nurbs.h` mis en commun avec l'IGES ; primitives en entités analytiques ; courbes et surfaces B-spline (rationnelles, transposition), restriction à une plage, pcurves 2D, contexte, unités, incertitude | aller-retour des entités géométriques par la lecture (PR 4) | 500 | 5 |
| 3 | `step/write-topology` | topologie, sens des faces, seams, arêtes dégénérées omises, pièces, assemblages de compounds avec instances partagées et placements (§ 12.2), `write_step`, rapport | aller-retour gbs → STEP → gbs sur tous les fichiers de test, assemblages compris | 600 | 6 |
| 4 | `step/write-occt-python` | lecture par OCCT, double aller-retour des formes OCCT, liaisons Python, documentation | job OCCT, pytest | 300 | 4 |

Total : environ 2 400 lignes, **28 h**, plus ce document. Ordre : P1 et P2
(indépendantes, sans lien avec l'écriture), puis 1, 2, 3, 4. La PR 1 peut
avancer en parallèle de P1 et P2.

---

## 11. Questions ouvertes

Pour chaque question : **R** = recommandation, **A** = alternatives.

1. **Surfaces toutes en B-spline ?** — **Décidé : non. Des primitives
   analytiques sont ajoutées à gbs (PR P1)**, présentes dans le BREP à côté des
   NURBS, et écrites en entités analytiques.
   R initiale : oui au premier jet. C'est valide et lu partout, et c'est exact pour tout
   ce que gbs produit.
   A : reconnaître les surfaces analytiques. Un plan bilinéaire et les NURBS
   exactes de la PR 3 se reconnaissent sans ambiguïté. Certains logiciels de
   CAO, et les opérations ultérieures de l'utilisateur, préfèrent les entités
   analytiques. Ce serait une PR supplémentaire d'environ 4 h.

2. **Emplacement : `gbs-io/step/`, namespace `gbs::step` ?**
   R : oui, à côté du lecteur.
   A : un en-tête `gbs-io/step_writer.h` hors du dossier.

3. **Courbe d'arête restreinte à la plage de l'arête ?**
   R : oui, pour lever l'ambiguïté d'une courbe fermée.
   A : la courbe entière, plus fidèle au partage de courbes du modèle.

4. **Pcurves écrites ?**
   R : oui par défaut (option `pcurves`). Elles sont exactes et sans coût de
   conversion, et certains lecteurs les exigent.
   A : les omettre, pour des fichiers plus petits.

5. **Sens des faces : `same_sense` de l'usage et borne `.F.` (écriture d'OCCT) ?**
   R : oui.
   A : retourner les boucles et écrire `same_sense = .T.` partout ; ou
   `oriented_face`.

6. **Compounds en assemblages de placements identité ?** — **Décidé :
   compounds en assemblages, avec la gestion des instances (PR P2)** : pièces
   partagées et placements réels ; option `flatten` conservée.
   R initiale : oui, avec l'option `flatten`.
   A : toujours une seule pièce ; ou une pièce par forme sans assemblage.

7. **Schéma AP214 par défaut, AP242 en option ?**
   R : oui. AP214 est le plus largement lu ; l'en-tête AP242 ne change pas les
   entités BREP écrites.
   A : AP242 par défaut, le schéma actuel de l'ISO.

8. **Arêtes et sommets libres omis et rapportés ?**
   R : oui.
   A : les écrire dans un `geometric_curve_set` (wireframe), ce qui ajouterait
   un type de représentation.

9. **Incertitude écrite : `brep_default_tolerance` ?**
   R : oui, modifiable par les options.
   A : la plus grande tolérance d'arête du modèle, plus honnête pour un modèle
   lu dans un fichier imprécis.

10. **Unité : celle du modèle, coordonnées non converties ?**
    R : oui.
    A : toujours le millimètre, en convertissant.

11. **Mettre en commun l'approximation avec l'IGES (`gbs-io/to_nurbs.h`) ?**
    R : oui ; l'export IGES l'inclut au lieu de la définir.
    A : la dupliquer dans l'écrivain STEP.

12. **Couleurs et calques ?**
    R : hors périmètre, comme à la lecture.
    A : les écrire depuis une table d'attributs, à créer dans le modèle.

---

## 12. PR préliminaires

### 12.1 P1 — primitives analytiques de gbs (`geom/primitives`)

**Besoin.** Écrire un `cylindrical_surface` plutôt qu'une B-spline suppose que
la face *repose* sur un cylindre, avec son repère et son rayon. Aujourd'hui, le
lecteur le sait (définitions de la PR 4), puis l'oublie en convertissant en
NURBS.

**Décision (12 octobre 2026) : des surfaces analytiques dans le BREP, à côté
des NURBS.** Le noyau BREP ne suppose nulle part que la géométrie est une NURBS.
Builders, pcurves, couture, `check`, volumes et boîtes englobantes passent par
l'interface `Surface` / `Curve`, et les tests du palier 1 construisent déjà des
faces et des solides sur des surfaces analytiques. Ajouter des primitives ne
demande donc pas de traiter deux familles partout. C'était l'intention de la
question 6 du palier 1 : « non au palier 1, polymorphie préservée pour les
ajouter plus tard ».

```cpp
// gbs/primitives.h — paramétrages d'ISO 10303-42, repère ax2 = {origine, axe, direction de référence}
template <std::floating_point T>
class CylindricalSurface : public Surface<T, 3>
{
public:
    CylindricalSurface(T radius, const ax2<T, 3> &frame, std::array<T, 4> box); // box = domaine (u en angle, v en longueur)
    auto value(T u, T v, size_t du = 0, size_t dv = 0) const -> point<T, 3> override; // dérivées exactes
    auto bounds() const -> std::array<T, 4> override;                                 // la boîte : un domaine fini
    auto parameters(const point<T, 3> &p) const -> std::array<T, 2>;                  // inverse en formule fermée
    auto to_nurbs() const -> BSSurfaceRational<T, 3>;                                  // NURBS exacte (bselementary.h)
    auto radius() const; auto frame() const;
};
// de même : Plane, ConicalSurface, SphericalSurface, ToroidalSurface ;
// courbes : Line et Circle3d existent ; ajout d'Ellipse3d
```

- **Paramétrage de STEP**, donc l'angle pour les directions de révolution. Les
  pcurves d'un fichier deviennent valables sur la surface lue. Celui d'un
  cercle sur un cylindre est un segment.
- **Domaine fini** donné à la construction (la boîte de la face, PR 4 de la
  lecture) : le BREP exige des bornes finies, et la fermeture se détecte comme
  pour une NURBS (u sur 2π).
- **Coniques périodiques** : `Circle3d` est aujourd'hui bornée à [0, 2π]. Un
  arc qui chevauche 0 (de 350° à 10°) doit pouvoir être une arête : le
  domaine des coniques devient la droite réelle, périodique de 2π.
- **Inverse en formule fermée** (`parameters`), déjà écrite dans
  `gbs-io/step/geometry.h` ; la projection générique reste disponible.
- **Conversion NURBS exacte** (`to_nurbs`) par `bselementary.h` (PR 3 de la
  lecture), pour les consommateurs qui n'acceptent que des NURBS.

**Ce qui change autour :**

| Composant | Changement |
|---|---|
| lecteur STEP | plan, cylindre, cône, sphère, tore, droite, cercle, ellipse gardés analytiques (révise la question 3 de la lecture) ; les plages d'arêtes sont les paramètres STEP ; surfaces de révolution, d'extrusion et B-spline converties comme avant ; l'arête de seam insérée (PR 7) devient une `CurveOnSurface` isoparamétrique, valable sur toute surface |
| `copy_shape` | transforme aussi le repère d'une primitive |
| export IGES | une face analytique est écrite par l'approximation existante (interpolation à paramètres conservés, pcurves valables), dans `approx_tol`, et comptée dans le rapport. Une conversion exacte (`to_nurbs` avec pcurves réinterpolées) pourra suivre |
| algorithmes du BREP | inchangés |

**Alternative écartée** : des NURBS qui portent leur définition analytique.
Elle ne touchait à aucun algorithme, mais gardait le paramétrage rationnel. Il
aurait fallu réinterpoler les pcurves à l'écriture STEP, et les intersections
analytiques du palier 2 n'en auraient pas profité.

### 12.2 P2 — provenance des instances (`brep/instances`)

**Besoin.** Écrire un assemblage avec ses pièces partagées et ses placements
suppose de savoir que deux solides du modèle sont deux instances de la même
pièce, et où chacune est placée. Le lecteur le sait (PR 9), puis l'oublie en
aplatissant.

**Proposition : un attribut du modèle, sans placement dans la topologie.**

```cpp
struct InstanceOf { std::uint64_t part; Placement<T> place; };  // pièce (clé), placement de l'instance
m.setInstance(shape, InstanceOf{…});  m.instance(shape) -> std::optional<InstanceOf>;
```

- **Table creuse** à côté des noms et des identifiants externes (PR 2 de la
  lecture), tenue par `erase`, `compact` et `append` de la même façon.
- **La géométrie reste aplatie et dupliquée** : aucun algorithme ne change, et
  la décision de la question 7 de la lecture tient.
- **Le lecteur** remplit l'attribut : clé = `#id` de la `product_definition` de
  la pièce, placement = placement composé de l'instance. `copy_shape` compose le
  placement d'une instance qu'il copie.
- **L'écrivain** regroupe les instances par clé. Il écrit la pièce une fois,
  à partir de la première instance ramenée par l'inverse de son placement.
  Avant d'écrire chaque occurrence, il vérifie que les instances d'une même
  clé concordent (même nombre d'entités) ; sinon, il les écrit comme pièces
  distinctes.
- **Alternative écartée** : de vrais placements dans le modèle (`Compound` à
  membres placés), économes en mémoire pour les gros assemblages. Mais
  l'exploration, `check`, les volumes, les boîtes englobantes, la couture et
  l'export IGES devraient tous appliquer les placements. C'est un chantier de
  palier 2.

