# Fichiers exacts d'un modèle BREP (`gbs-io/brep_json.h`)

| | |
|---|---|
| PR | branche `feat/brep-serialization` |
| Fichiers | `gbs-io/brep_json.h` (nouveau), `gbs-brep/seamless.h` (nouveau, venu de `gbs-brep/fillet.h`), `gbs-brep/fillet.h`, `gbs-brep/brep`, `gbs-brep/model.h` (`nextPart`, `reservePartKeys`), `python/gbsBindBrep.cpp`, `tests/tests_io_brep_json.cpp` (nouveau), `python/tests/test_brep_json.py` (nouveau), `docs/index.rst`, `ReadMe.md`, `News.md` |
| Statut | À valider (§ 11) |

## 1. Objectif

Enregistrer un modèle BREP (`gbs::brep::Model<T>`) et le relire **tel qu'il était**, pas une
copie de ses formes. Une application de CAO (CadSmith, besoin G11) garde dans son document des
identifiants d'entités, des racines, des corps consommés : ouvrir le fichier doit rendre les
mêmes entités aux mêmes identifiants, la même géométrie au bit près.

STEP ne le permet pas, ni à l'écriture ni à la lecture : l'écrivain n'écrit que des racines,
renumérote et approche ce qu'AP214 n'a pas ; le lecteur réordonne les formes, coupe des arêtes
aux coutures, recalcule les pcurves. STEP reste le format d'**échange** ; ce format-ci est celui
du **document**.

## 2. API

```cpp
namespace gbs::io
{
    inline constexpr std::string_view brep_json_format = "gbs-brep";
    inline constexpr std::uint32_t brep_json_version = 1;
    template <typename T> concept brep_json_real = std::same_as<T, double> || std::same_as<T, float>;

    enum class BrepReadErrc : std::uint8_t { Syntax, NotBrep, NewerVersion, Invalid };
    class BrepReadError : public TextError   // code(), where() (pointeur JSON), reason() ; copie sans exception
    class BrepWriteError : public BRepError; // une géométrie que le format ne tient pas, un id hors de sa table

    struct BrepWriteOptions { bool pretty{false}; };

    template <brep_json_real T> auto write_brep(const brep::Model<T> &, const BrepWriteOptions & = {}) -> std::string;
    template <brep_json_real T> void write_brep(const brep::Model<T> &, std::ostream &, const BrepWriteOptions & = {});
    template <brep_json_real T> void write_brep_file(const brep::Model<T> &, const std::filesystem::path &, const BrepWriteOptions & = {});
    template <brep_json_real T = double> auto read_brep(std::string_view text, std::string_view source = "BREP JSON") -> brep::Model<T>;
    template <brep_json_real T = double> auto read_brep(std::istream &, std::string_view source = "BREP JSON stream") -> brep::Model<T>;
    template <brep_json_real T = double> auto read_brep_file(const std::filesystem::path &) -> brep::Model<T>;
}
```

- `write_brep_file` écrit par `write_text_atomic` : le fichier est remplacé entier ou laissé tel
  qu'il était, y compris quand l'écriture est refusée (`BrepWriteError` avant toute écriture).
- `read_brep` rend un modèle neuf : les identifiants exacts n'ont de sens que dans un modèle vide.
- Erreurs comme `json.h` et `obj.h` (des exceptions, pas des `BuildResult`) :
  `BrepReadError` pour un texte qui n'est pas un modèle lisible, `std::filesystem::filesystem_error`
  pour un fichier illisible ou impossible à écrire, `std::ios_base::failure` pour un flux.
  `what()` vaut « source: where: reason », par exemple
  `part.brep: /edges/12/v/0: a vertex index below 40 or null expected`.

`Model` gagne deux fonctions, pour que la clé de pièce suivante (`newPart`) survive au fichier :
`nextPart()` (lecture) et `reservePartKeys(next)` (marque données les clés au-dessous de `next`).

Python (`gbs.brep`) : `write_brep(model, path, pretty=False)`, `write_brep_text(model,
pretty=False)`, `read_brep(path)`, `read_brep_text(text)` ; `BrepReadError` (un `ValueError`, avec
`code` de type `BrepReadErrc`, `where`, `reason`), `BrepWriteError` (un `BRepError`), `OSError`
pour un fichier. Lecture et écriture sans le GIL. Un `Model` se **pickle** par ce texte (exact,
comme `copy.deepcopy`).

## 3. Choix du format : JSON

JSON exact sur RapidJSON (`gbs-io/json.h`), plutôt qu'un binaire, sur mesures :

| Modèle | STEP | Lecture STEP | JSON | Écriture | Lecture (dont analyse JSON) | gzip du JSON | Binaire estimé |
|---|---|---|---|---|---|---|---|
| `rotor37.step` (6 faces, 68 géométries) | 301 Ko | 19 ms | 268 Ko | 0,7 ms | 0,7 ms (0,7) | 97 Ko | 109 Ko |
| Rotor de deux pales, 14 faces, B-splines denses | 18,8 Mo | 9,6 s | 10,7 Mo | 23 ms | 24 ms (24) | 3,5 Mo | 4,4 Mo |
| Synthétique : 1 000 blocs percés, 30 000 faces, 54 000 arêtes, 18 096 géométries | — | — | 16,6 Mo | 46 ms | 57 ms (41) | 1,1 Mo | 6,7 Mo |

Mesures en Release (Apple M5 Max, clang 21), le meilleur de 5 essais ; « binaire estimé » compte 8
octets par réel et 4 par entier du fichier, sans en-têtes ni types : la borne basse d'un
binaire. Le programme de mesure n'est pas dans le dépôt (il lit des pièces qui n'y sont pas).

- **Vitesse :** 0,3 à 0,45 Go/s dans les deux sens. Lire le JSON coûte de 100 à 400 fois moins
  que lire le STEP de la même pièce, et moins que le maillage d'affichage qui suit. Un binaire
  diviserait un temps déjà négligeable.
- **Taille :** de 0,6 à 0,9 fois le STEP. Un binaire serait 2,5 fois plus petit, mais le JSON
  compressé l'est plus encore (3,5 Mo contre 4,4 Mo ; 1,1 Mo contre 6,7 Mo) : une application
  qui veut des fichiers petits compresse son document (une archive), ce que gbs n'a pas à faire
  (aucune dépendance de plus, `project-runtime-deps-only`).
- **Ce que le binaire coûterait :** un second format à maintenir et à valider contre une
  entrée hostile (boutisme, tailles, alignement), illisible pour un humain et pour un outil.
- **Ce que JSON donne :** un fichier qu'on lit et compare (`diff`), qu'un script Python ouvre
  avec `json`, une version et des clés explicites ; un analyseur éprouvé, sans récursion.

## 4. Le format (version 1)

```json
{"format": "gbs-brep", "version": 1, "real": "double", "unit_scale": 1.0, "next_part": 1,
 "geometry": [ {"type": "bspline_curve", "dim": 2, "degree": 2, "knots": [...], "poles": [[...], ...]},
               {"type": "plane", "dim": 3, "frame": [[o], [x], [y], [z]], "box": [u1, u2, v1, v2]},
               {"type": "curve_on_surface", "dim": 3, "pcurve": 0, "surface": 1}, ... ],
 "vertices":  [ {"p": [x, y, z], "tol": 1e-07}, null, ... ],
 "edges":     [ {"curve": 2, "u": [u1, u2], "v": [0, 1], "tol": 1e-07, "degenerate": false, "same_parameter": true}, ... ],
 "wires":     [ {"coedges": [[edge, "F", pcurve], [edge, "R", null], ...], "closed": true}, ... ],
 "faces":     [ {"surface": 1, "wires": [0, 3], "tol": 1e-07, "natural_bounds": false, "name": "top"}, ... ],
 "shells":    [ {"faces": [[face, "F"], ...], "closed": true}, ... ],
 "solids":    [ {"outer": 0, "voids": []}, ... ],
 "compounds": [ {"shapes": [["solid", 0], ["face", 7]], "instance": {"part": 1, "placement": {"r": [[...],[...],[...]], "t": [...]}}}, ... ]}
```

- **Tables :** une par type d'entité, dans l'ordre des indices, **les mortes comprises** (`null`).
  L'indice dans le tableau est l'identifiant ; un identifiant stocké vaut son indice, ou `null`
  s'il est invalide (`npos`).
- **Attributs** dans l'objet de l'entité, seulement s'ils existent : `name`, `external_id`,
  `instance` (`part`, `placement`). Une entité morte n'en a pas (`erase` les efface).
- **Géométrie :** un tableau `geometry` de tous les objets distincts (par adresse), chacun
  avant ceux qui le citent ; les entités et les géométries composées les citent par leur indice.
  `type` et `dim` disent la classe ; les clés suivent la classe (§ 5.3).
- **Orientations :** `"F"` ou `"R"`. **Réels :** la forme décimale la plus courte qui se relit au
  bit près (Grisu de RapidJSON, relu avec `kParseFullPrecisionFlag`) ; `-0.0`, `Infinity`,
  `-Infinity` et `NaN` sont écrits (JSON étendu de RapidJSON, que le module `json` de Python lit).
- `pretty` indente le texte, un tableau de nombres sur une ligne ; le modèle relu est le même.

**Versions.** `version` est un entier ; ce code écrit 1 et lit 1. Une version plus grande est
refusée (`NewerVersion`, « a newer gbs is needed ») **avant** tout autre contrôle : un fichier
plus récent peut tenir ce que ce code ne connaît pas. Toute évolution du format (une clé, un type
de géométrie, un sens nouveau) prend une version de plus ; un lecteur plus récent lit les
anciennes.

## 5. Exactitude

### 5.1 Identifiants et entités mortes

Les tables gardent leur taille et leurs trous : une entité morte est relue morte (ajoutée par
défaut, puis `erase`), donc `capacity`, `count`, `alive` et chaque identifiant sont ceux du
modèle écrit. Le contenu d'une entité morte n'est observable par aucune fonction de `Model` et
n'est pas écrit. Les identifiants stockés sont écrits tels quels, y compris vers une entité morte
(`erase` ne propage pas) ou invalides : le modèle relu a les mêmes incohérences que l'original.
Un identifiant **au-delà** de sa table (une faute de programmation) est refusé à l'écriture
(`BrepWriteError`) : la lecture le refuserait.

**Compactage explicite :** qui veut un fichier sans trous appelle `Model::compact()` avant
d'écrire, et traduit ses identifiants par l'`IdRemap` rendu. L'écrivain ne compacte jamais seul.

### 5.2 Partages

La géométrie est dédupliquée par l'adresse de l'objet le plus dérivé : un objet tenu par
plusieurs entités (une surface partagée par les morceaux d'une face coupée, une courbe par deux
arêtes) ou cité par une autre géométrie (la surface d'une courbe sur surface, la base d'un
décalage) est écrit une fois et relu en un seul `shared_ptr`. Deux objets égaux mais distincts
restent distincts. Les partages topologiques sont des identifiants : une arête partagée par deux
co-arêtes le reste.

Les fonctions de loi (`BSCfunction` d'une courbe décalée ou reparamétrée) sont des valeurs,
écrites dans l'objet qui les tient.

### 5.3 Géométrie : les classes tenues, à l'identique

Le type dynamique **exact** est reconnu (`typeid`), pas un `dynamic_cast` : une classe dérivée
d'une B-spline n'est pas écrite comme une B-spline (elle perdrait ce qu'elle ajoute).

| `type` | Classe | Données écrites |
|---|---|---|
| `bspline_curve`, `nurbs_curve` (2D, 3D) | `BSCurve`, `BSCurveRational` | degré, nœuds à plat, pôles tels que stockés (homogènes `(w x, w y, w z, w)` pour une NURBS), `bounds` si elles diffèrent des nœuds extrêmes |
| `bspline_surface`, `nurbs_surface` (2D, 3D) | `BSSurface`, `BSSurfaceRational` | degrés, nœuds U et V, pôles (u d'abord), `bounds` si elles diffèrent |
| `line`, `ellipse` | `LineCurve`, `EllipseCurve` | origine et vecteur ; rayons et repère ; intervalle |
| `plane`, `cylinder`, `cone`, `sphere`, `torus` | primitives | paramètres, repère `(o, x, y, z)`, boîte |
| `trimmed_curve`, `composite_curve`, `reversed_curve`, `reparametrized_curve`, `extended_curve`, `curve_on_surface` (2D, 3D) | courbes composées | leurs géométries (indices) et paramètres |
| `offset_curve` (2D, 3D) | `CurveOffset2D/3D<T, BSCfunction<T>>` | base, loi, direction (3D) |
| `moved_curve`, `moved_surface` | `brep::MovedCurve`, `brep::MovedSurface` | base, placement |
| `revolution_surface`, `offset_surface`, `seamless_surface` | `SurfaceOfRevolution`, `OffsetSurface`, `brep::SeamlessSurface` | courbe et axe, base et distance, base et périodes |

Ce qui oblige à des soins particuliers, parce qu'un constructeur ne rend pas toujours ce qu'on
lui donne :

- **Pôles rationnels** écrits homogènes : diviser puis remultiplier par le poids ne rendrait pas
  les mêmes bits.
- **Bornes d'une B-spline** distinctes de ses nœuds : `trim(u1, u2, false)` d'une courbe les
  change seules, et le relecteur l'appelle. Pour une surface, `trimU` laisse nœuds et pôles quand
  les bornes sont à moins de `knot_eps` des nœuds extrêmes (son retour anticipé) et pose les
  bornes : c'est le seul cas qui arrive, le relecteur fait de même. Des bornes V distinctes ne
  peuvent pas être rendues (`trimV` passe par `invertUV`, qui les recalcule) ; aucune opération
  n'en produit, l'écrivain les refuse.
- **Repères des primitives** écrits avec leurs quatre vecteurs : celui d'une primitive déplacée
  (`Frame3::moved`) n'est pas exactement orthonormé, et le constructeur depuis un `ax2` le
  recalculerait.
- **Boîte d'une primitive** posée par `with_box` après construction sur une boîte quelconque :
  `with_box` (qu'utilisent les opérations qui déplacent une couture) ne vérifie que la finitude et
  le non-vide, alors que le constructeur de la sphère refuse une latitude au-delà des pôles d'un
  ulp.
- **Placements composés :** le constructeur de `MovedCurve` compose un placement avec celui d'une
  base déjà déplacée ; la base écrite ne l'est jamais, et le relecteur refuse un fichier qui le
  dirait.
- **Direction d'une courbe décalée 3D :** le constructeur la normalise ; l'écrivain refuse une
  direction que sa normalisation ne rend pas au bit près (aucune opération de gbs n'en crée).

**Refusées à l'écriture** (`BrepWriteError`, avec l'entité : « edge 12: curve of class … cannot
be written ») : les anciennes `Line` et `Circle` (leurs constructeurs normalisent l'axe et le
rayon : relues, elles ne seraient plus les mêmes), toute classe d'une application, les quatre
classes internes qui ne restent jamais dans un modèle (`FieldCurve`, `LoopCurve`, `PathCurve`,
`DefSurface`).

Inventaire de ce que gbs met dans un modèle (relu dans le code) : B-splines et NURBS 2D et 3D, les
sept primitives, `CurveOnSurface` (dont la surface peut être une `SeamlessSurface`, une
`OffsetSurface`, une `MovedSurface`), `CurveComposite` de `CurveTrimmed` (2D et 3D, `make_thick`),
`MovedCurve` et `MovedSurface` (`copy_shape` avec un placement), `OffsetSurface`. Les autres
classes du tableau ne viennent que de l'API ; elles sont tenues pour qu'un modèle construit à la
main s'enregistre aussi.

### 5.4 `SeamlessSurface` sort de `fillet.h`

La vue sans couture d'une face B-spline fermée (`fillet_detail::SeamlessSurface`) devient la
surface des courbes de contact d'un congé qui passe la couture (le congé du pied d'une pale) :
elle **reste dans le modèle**. Pour que l'écrivain la connaisse sans inclure le noyau des congés
(4 000 lignes, 3 à 4 Go à compiler), elle passe telle quelle dans `gbs-brep/seamless.h`,
`gbs::brep::SeamlessSurface`. `fillet.h` l'inclut ; ses usages non qualifiés dans `fillet_detail`
la trouvent dans l'espace de noms englobant : aucune autre ligne ne change.

### 5.5 Réels

Bit pour bit pour tout double fini, les sous-normaux, `-0` et les infinis compris (vérifié sur
3,2 millions de doubles tirés au hasard, puis par le test). Un NaN revient NaN, sans son signe ni
sa charge. En `float`, le texte écrit le double exact du float, relu et vérifié float ; `real`
dit le type, et un modèle `float` n'est pas lu en `double` ni l'inverse.

## 6. Lire une entrée hostile

Rien n'est cru avant d'être vérifié ; aucune entrée ne mène à un comportement indéfini :

- **Analyse** sans récursion (`kParseIterativeFlag`) : un million de crochets imbriqués est une
  faute de syntaxe, pas un dépassement de pile. Un octet NUL est refusé.
- **Schéma strict :** chaque objet n'a que ses clés, chacune une fois, les obligatoires présentes ;
  chaque valeur a son type JSON. Une clé inconnue est une faute : une faute de frappe ne perd
  pas un nom en silence (une clé nouvelle vient avec une version nouvelle).
- **Indices :** un identifiant est `null` ou un entier sous la taille de sa table (tailles lues
  avant les entités) ; une géométrie ne cite qu'une géométrie **avant** elle, ni elle-même : aucun
  cycle, aucune récursion à la lecture ; une référence est vérifiée de la bonne sorte
  (« a 3D curve expected, geometry 4 is a 3D surface »).
- **B-splines :** nœuds finis et dans l'ordre, en nombre `pôles + degré + 1`, au moins `degré + 1`
  pôles ; pôles en nombre `nu × nv` ; vérifié avant le constructeur (qui écrirait sur `cerr`).
  Le reste des invariants est vérifié par les constructeurs de gbs, dont le refus
  (`std::invalid_argument`, `std::length_error`) devient une `BrepReadError` à la géométrie.
- **Attributs :** nom non vide, clé de pièce entre 1 et `next_part - 1`, unité finie et > 0.

Le modèle lu est un modèle que l'API de `Model` aurait pu construire, pas un solide valide : la
lecture ne lance pas `check`, comme `Model::add` ne le fait pas.

## 7. Validation

`tests_io_brep_json` (15 cas, 181 150 assertions, 1,5 s) compare le modèle relu au modèle écrit
**sans passer par l'écrivain** : table par table, entité par entité, chaque réel au bit près ;
pour chaque géométrie la classe exacte, ses données (pôles, nœuds, repère, rayons, placements),
ses valeurs et dérivées en des points ; et les partages : les paires d'objets rencontrées forment
une bijection (un objet partagé d'un côté l'est de l'autre, et nulle part ailleurs). Puis le
texte du modèle relu est le texte écrit, compact et indenté.

- `every_geometry_class` : chaque classe du tableau du § 5.3, en 2D et en 3D, imbriquées et
  partagées (plus de 10 partages retrouvés), bornes de B-splines distinctes des nœuds, sphère
  au-delà des pôles par `with_box`, repère non orthonormé ;
- `dead_entities_ids_and_attributes` : entités mortes de chaque table, encore citées ; ids
  invalides ; noms à guillemets, retours, octet NUL, UTF-8 ; `external_id` extrêmes ; instance et
  clé de pièce suivante ; unité ; puis le même modèle compacté ;
- `numbers_bit_for_bit` (30 000 réels tirés au hasard, sous-normaux, ±0, ±∞, NaN), `float_model` ;
- `booleans_and_fillets` : bloc moins cylindre, boîte moins sphère, congé d'un bloc, congé du pied
  d'une pale (des `SeamlessSurface`), relus avec le même `check` et le même volume au bit près ;
- `step_files` : chacun des fichiers STEP de `tests/in/step`, puis tous dans un modèle dont un
  tiers des racines est effacé ;
- `files_and_streams` ; refus : `version_and_format`, `truncated_and_malformed` (chaque préfixe
  du texte, un octet NUL, l'imbrication profonde), `corrupt_values` (22 fautes, chacune à sa place
  et avec sa raison), `forward_reference_refused`, `constructors_refusals_reported`,
  `random_damage` (3 000 textes abîmés au hasard : lus, ou refusés par `BrepReadError`, jamais
  autre chose), `refused_writes` (ancienne `Line`, classe d'application, `Circle3d` imbriqué,
  id hors table ; le fichier en place reste intact).

Une mutation du relecteur (oublier les bornes d'une courbe) fait échouer deux cas : la
comparaison voit ce que l'écrivain et le relecteur perdraient ensemble.

Python : `test_brep_json.py` (fichier, texte, pickle, erreurs de lecture avec `code`, `where`,
`reason`, erreur d'écriture).

## 8. Usage attendu par CadSmith

```cpp
// Enregistrer : le modèle dans le document, à côté des métadonnées de l'application
const std::string brep = gbs::io::write_brep(model);          // BrepWriteError : rien d'écrit
// … le document les met ensemble (par exemple {"cadsmith": {...}, "brep": <texte>} ou une archive),
// puis gbs::write_text_atomic(path, document)
// Ouvrir
gbs::brep::Model<double> m = gbs::io::read_brep<double>(brep_text, path_for_messages);
// BrepReadError : e.code() == NewerVersion -> « fichier d'une version plus récente » ;
// Syntax / NotBrep / Invalid -> « fichier abîmé » avec e.where() et e.reason()
```

Les identifiants (`SolidId`, `ShapeId`…) que le document garde désignent les mêmes entités après
lecture : les racines, les corps consommés et la visibilité se stockent comme des indices. Le
texte du modèle peut aussi s'insérer tel quel dans un JSON plus grand (c'est un objet JSON
complet) ; CadSmith garde alors sa propre version de document à côté de celle de gbs.

## 9. Ce qui n'est pas couvert

- Les classes refusées du § 5.3 (anciennes `Line`, `Circle`, classes d'une application). Une
  application qui en met dans un modèle doit les convertir (`LineCurve`, `EllipseCurve`) ou les
  approcher avant d'écrire.
- Le signe et la charge d'un NaN.
- Le contenu d'une entité morte (inobservable).
- L'enregistrement incrémental ou partiel (un sous-ensemble de formes) : on écrit le modèle entier ;
  un sous-ensemble passe par `copy_shape` dans un modèle neuf.
- La compression : à l'application.
- Le `long double` (`brep_json_real` : `double` et `float`).

## 10. Écarts

Aucun document de conception ne prévoyait ce format ; deux changements hors de `gbs-io` :
`SeamlessSurface` déplacée (§ 5.4), et `Model::nextPart` / `reservePartKeys` (la clé de pièce
suivante n'était pas lisible : sans elle, une pièce créée après relecture pourrait reprendre la
clé d'une pièce effacée avant l'enregistrement).

## 11. À valider

1. JSON plutôt que binaire (§ 3), la compression laissée à l'application.
2. Les entités mortes gardées par défaut, le compactage au choix de l'appelant (§ 5.1).
3. Le schéma strict : une clé inconnue est une faute, toute évolution prend une version (§ 4, § 6).
4. Les refus à l'écriture (anciennes `Line` et `Circle`, classes d'application, ids hors table)
   plutôt qu'une écriture approchée (§ 5.3).
5. `SeamlessSurface` publique dans `gbs-brep/seamless.h` (§ 5.4).
6. Les deux fonctions ajoutées à `Model` (§ 10).
7. Les noms de l'API (`write_brep`, `read_brep`, `BrepReadError`…) et l'espace `gbs::io`, comme
   `obj.h`.
