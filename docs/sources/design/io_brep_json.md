# Fichiers exacts d'un modèle BREP (`gbs-io/brep_json.h`)

| | |
|---|---|
| PR | branche `feat/brep-serialization` |
| Fichiers | `gbs-io/brep_json.h` (nouveau), `gbs-brep/equal.h` (nouveau), `gbs-brep/seamless.h` (nouveau, venu de `gbs-brep/fillet.h`), `gbs-brep/fillet.h`, `gbs-brep/brep`, `gbs-brep/model.h` (`nextPart`, `reservePartKeys`, `BuildErrc::UnsupportedVersion`), `gbs-brep/builders.h` (`unwrap` d'un `BuildResult<void>`), `gbs/curveoffset.h` (`unit_direction`), `gbs/curvecomposite.h` (revue, § 13), `python/gbsBindBrep.cpp`, `tests/tests_io_brep_json.cpp` (nouveau), `tests/tests_curves.cpp`, `python/tests/test_brep_json.py` (nouveau), `docs/index.rst`, `docs/Doxyfile`, `ReadMe.md`, `News.md` |
| Exigences | besoin G11 de CadSmith et sa revue (`CadSmith/docs/revues/plan-v1.md`, section « G11 ») : réponses au § 10 |
| Statut | À valider (§ 12, § 13) |

## 1. Objectif

Enregistrer un modèle BREP (`gbs::brep::Model<T>`) et le relire **tel qu'il était**, pas une
copie de ses formes. Une application de CAO (CadSmith, besoin G11) garde dans son document des
identifiants d'entités (racines, corps consommés, formes masquées) : ouvrir le fichier doit rendre
les mêmes entités aux mêmes identifiants, la même géométrie au bit près, les mêmes partages.

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
    inline constexpr unsigned brep_json_parse_flags; // pleine précision, sans récursion, UTF-8 vérifié
    template <typename T> concept brep_json_real = std::same_as<T, double> || std::same_as<T, float>;
    struct BrepWriteOptions { bool pretty{false}; };

    // écriture : fonction pure d'un modèle const
    auto write_brep(const brep::Model<T> &, const BrepWriteOptions & = {}) -> brep::BuildResult<std::string>;
    auto write_brep_compact(const brep::Model<T> &, const BrepWriteOptions & = {}) -> brep::BuildResult<std::pair<std::string, brep::IdRemap>>;
    auto write_brep(const brep::Model<T> &, std::ostream &, const BrepWriteOptions & = {}) -> brep::BuildResult<void>;
    auto write_brep_file(const brep::Model<T> &, const std::filesystem::path &, const BrepWriteOptions & = {}) -> brep::BuildResult<void>;
    auto write_brep_value(const brep::Model<T> &, W &rapidjson_writer) -> brep::BuildResult<void>;

    // lecture : un modèle neuf
    auto read_brep<T = double>(std::string_view text, std::string_view source = "BREP JSON") -> brep::BuildResult<brep::Model<T>>;
    auto read_brep<T = double>(std::istream &, std::string_view source = "BREP JSON stream") -> brep::BuildResult<brep::Model<T>>;
    auto read_brep_file<T = double>(const std::filesystem::path &) -> brep::BuildResult<brep::Model<T>>;
    auto read_brep_value<T = double>(const rapidjson::Value &, std::string_view source = "BREP JSON value") -> brep::BuildResult<brep::Model<T>>;
}
namespace gbs::brep
{
    auto first_difference(const Model<T> &, const Model<T> &) -> std::optional<std::string>; // gbs-brep/equal.h
    bool equal(const Model<T> &, const Model<T> &);
}
```

- **Erreurs : des `BuildResult`**, comme le lecteur et l'écrivain STEP, aucune exception
  (sauf `std::bad_alloc`). L'écriture refuse avant d'écrire quoi que ce soit, en nommant
  l'entité (`BuildError::shapes` et le message, « edge 12: curve of class … cannot be written ») ;
  la lecture dit où, en pointeur JSON (« part.brep: /edges/12/v/0: a vertex index below 40 or
  null expected »). Codes : `UnsupportedEntity` et `InvalidId` à l'écriture, `InvalidFile` et
  `UnsupportedVersion` (nouveau) à la lecture ; un fichier qui ne s'ouvre pas ou ne s'écrit pas :
  `InvalidFile`.
- `write_brep_file` écrit par `write_text_atomic` : le fichier est remplacé entier ou laissé tel
  qu'il était, y compris quand l'écriture est refusée.
- `write_brep_compact` écrit le modèle tel que `Model::compact()` le rendrait (une copie compactée :
  le modèle donné n'est pas modifié) et rend cet `IdRemap`, pour renuméroter dans le même passage
  les identifiants tenus hors du modèle.
- `write_brep_value` et `read_brep_value` placent le modèle dans un document de l'application :
  une valeur JSON (après une clé de l'enveloppe) écrite par le `rapidjson::Writer` de
  l'application, relue d'un document analysé avec `brep_json_parse_flags` (§ 9).
- `read_brep` rend un modèle neuf ; pour l'ajouter à un autre, `Model::append` et son `IdRemap`.
- `brep::equal` et `brep::first_difference` (§ 7) : l'égalité exacte de deux modèles, pour les
  tests de gbs et ceux de l'application.

`Model` gagne `nextPart()` et `reservePartKeys(next)`, pour que la clé de pièce suivante
(`newPart`) survive au fichier ; `BuildErrc` gagne `UnsupportedVersion` ; `unwrap` accepte un
`BuildResult<void>`.

Python (`gbs.brep`) : `write_brep(model, path, pretty=False)`, `write_brep_text(model,
pretty=False)`, `write_brep_compact_text(model, pretty=False)` → `(texte, IdRemap)`,
`read_brep(path)`, `read_brep_text(text)`, `equal(a, b)`, `first_difference(a, b)` ; une faute
lève `BRepError` avec son `code` (`BuildErrc.invalid_file`, `unsupported_version`,
`unsupported_entity`, `invalid_id`) et `shapes`. Sans le GIL. Un `Model` se **pickle** par ce
texte (exact, comme `copy.deepcopy`).

## 3. Choix du format : JSON

JSON exact sur RapidJSON (`gbs-io/json.h`), plutôt qu'un binaire, sur mesures :

| Modèle | STEP | Lecture STEP | JSON | Écriture | Lecture (dont analyse JSON) | gzip du JSON | Binaire estimé |
|---|---|---|---|---|---|---|---|
| `rotor37.step` (6 faces, 68 géométries) | 301 Ko | 17 ms | 268 Ko | 0,5 ms | 0,7 ms (0,6) | 97 Ko | 109 Ko |
| Rotor de deux pales, 14 faces, B-splines denses | 18,8 Mo | 9,6 s | 10,7 Mo | 20 ms | 24 ms (presque tout) | 3,5 Mo | 4,4 Mo |
| Synthétique : 1 000 blocs percés, 30 000 faces, 54 000 arêtes, 18 096 géométries | — | — | 16,6 Mo | 36 ms | 44 ms (32) | 1,1 Mo | 6,7 Mo |

Mesures en Release (Apple M5 Max, clang 21), le meilleur de 5 essais ; « binaire estimé » compte 8
octets par réel et 4 par entier du fichier, sans en-têtes ni types : la borne basse d'un
binaire. Le programme de mesure n'est pas dans le dépôt (il lit des pièces qui n'y sont pas) ;
d'une série à l'autre, les temps varient de 20 %. `brep::equal` du modèle synthétique et de sa
relecture : 69 ms.

- **Vitesse :** 0,35 à 0,55 Go/s dans les deux sens, l'écriture en flux (l'API SAX de RapidJSON,
  aucun DOM en plus du modèle) ; la lecture analyse un DOM. Une pièce de 50 000 faces s'écrit et se
  relit en moins d'une centaine de millisecondes chacune : le budget de la sauvegarde automatique de
  CadSmith (moins d'une seconde) est tenu de loin. Lire le JSON coûte de 100 à 400 fois moins que
  lire le STEP de la même pièce.
- **Taille :** de 0,6 à 0,9 fois le STEP. Un binaire serait 2,5 fois plus petit, mais le JSON
  compressé l'est plus encore (3,5 Mo contre 4,4 Mo ; 1,1 Mo contre 6,7 Mo) : une application
  qui veut des fichiers petits compresse son document (une archive), ce que gbs n'a pas à faire
  (aucune dépendance de plus à l'exécution).
- **Ce que le binaire coûterait :** un second format à maintenir, à spécifier et à valider contre
  une entrée hostile, illisible pour un humain et pour un outil.
- **Ce que JSON donne :** un fichier qu'on lit et compare (`diff`), qu'un script Python ouvre
  avec `json` (le texte est du JSON standard : ni NaN, ni infini, UTF-8 vérifié), une version et
  des clés explicites ; un analyseur éprouvé, sans récursion.

## 4. Spécification du format, version 1

Un objet JSON, ses membres dans un ordre quelconque à la lecture (l'écrivain les écrit dans cet
ordre) ; chaque membre est obligatoire sauf mention ; aucun autre membre n'est admis.

| Membre | Valeur |
|---|---|
| `format` | `"gbs-brep"` |
| `version` | entier, 1 |
| `real` | `"double"` ou `"float"` : le `T` du modèle |
| `unit_scale` | réel > 0 : millimètres par unité du modèle |
| `next_part` | entier de 1 à 2⁶³ : la clé de pièce que `newPart()` donnera ; toute clé d'instance est en dessous (borné pour que `newPart()` ne puisse pas déborder) |
| `geometry` | tableau des géométries (ci-dessous) ; chacune ne cite que des géométries d'indice inférieur, à une profondeur de 256 au plus (`brep_json_max_depth` : une géométrie qui n'en cite aucune est à la profondeur 1, une autre à un de plus que la plus profonde qu'elle cite) |
| `vertices`, `edges`, `wires`, `faces`, `shells`, `solids`, `compounds` | un tableau par table, à l'indice de l'entité : son objet, ou `null` pour une entité morte |

**Entités.** Un identifiant est l'indice d'une entité de sa table (vivante ou morte), ou `null`
pour un identifiant invalide. Une orientation vaut `"F"` (Forward) ou `"R"` (Reversed).

| Table | Membres |
|---|---|
| `vertices` | `p` : `[x, y, z]` ; `tol` |
| `edges` | `curve` : courbe 3D ou `null` ; `u` : `[u1, u2]` ; `v` : `[v1, v2]` (sommets) ; `tol` ; `degenerate`, `same_parameter` : booléens |
| `wires` | `coedges` : `[[arête, orientation, pcurve 2D ou null], …]` dans l'ordre ; `closed` |
| `faces` | `surface` : surface 3D ou `null` ; `wires` : `[contour, …]` (l'extérieur d'abord) ; `tol` ; `natural_bounds` |
| `shells` | `faces` : `[[face, orientation], …]` dans l'ordre ; `closed` |
| `solids` | `outer` : coque ; `voids` : `[coque, …]` |
| `compounds` | `shapes` : `[[type, indice], …]` dans l'ordre, `type` parmi `"vertex"`, `"edge"`, `"wire"`, `"face"`, `"shell"`, `"solid"`, `"compound"` |

Chaque entité vivante porte en plus, seulement s'ils existent : `name` (chaîne UTF-8 non vide),
`external_id` (entier signé 64 bits), `instance` : `{"part": clé, "placement": placement}`.
Un placement : `{"r": [[3 réels] × 3 lignes de R], "t": [3 réels]}`.

**Géométries.** Un objet `{"type": …, "dim": 2 ou 3, …}` ; une référence à une autre géométrie
est son indice dans `geometry`, d'une sorte vérifiée (courbe ou surface, 2D ou 3D) ; un repère
`frame` : `[o, x, y, z]`, quatre points écrits tels que gbs les tient ; une loi `function` : un
objet `{"degree", "knots", "poles": [réels], "bounds"?}` (une B-spline 1D, écrite en place). Le
degré d'une B-spline (courbe, surface, loi) est au plus 24 (`brep_json_max_degree`, la limite
`bspline_stack_max_degree` de gbs : au-delà, l'évaluation retombe sur une récursion dont le coût
double à chaque degré, deux heures par point au degré 40).

| `type` | `dim` | Membres | Classe |
|---|---|---|---|
| `bspline_curve`, `nurbs_curve` | 2, 3 | `degree`, `knots` (à plat, multiplicités comprises), `poles` (`dim` réels ; `dim + 1` pour une NURBS : homogènes `(w x, w y, w z, w)`, comme gbs les tient), `bounds` : `[u1, u2]` si elles diffèrent des nœuds extrêmes | `BSCurve`, `BSCurveRational` |
| `bspline_surface`, `nurbs_surface` | 2, 3 | `degree` : `[p, q]`, `knots` : `[U, V]`, `poles` (u d'abord : le pôle (i, j) à i + j nu), `bounds` : `[u1, u2, v1, v2]` si U diffère des nœuds | `BSSurface`, `BSSurfaceRational` |
| `line` | 3 | `origin`, `vector`, `range` | `LineCurve` |
| `ellipse` | 3 | `radii` : `[r1, r2]`, `frame`, `range` | `EllipseCurve` |
| `plane` | 3 | `frame`, `box` | `Plane` |
| `cylinder`, `sphere` | 3 | `radius`, `frame`, `box` | `CylindricalSurface`, `SphericalSurface` |
| `cone` | 3 | `radius`, `semi_angle`, `frame`, `box` | `ConicalSurface` |
| `torus` | 3 | `major_radius`, `minor_radius`, `frame`, `box` | `ToroidalSurface` |
| `trimmed_curve` | 2, 3 | `basis`, `range` | `CurveTrimmed` |
| `composite_curve` | 2, 3 | `curves` : au moins une, chacune de bornes croissantes (`u1 ≤ u2` : `CurveComposite` cherche son morceau dans des paramètres croissants) | `CurveComposite` |
| `reversed_curve` | 2, 3 | `basis` | `CurveReversed` |
| `reparametrized_curve` | 2, 3 | `basis`, `function` | `CurveReparametrized` |
| `extended_curve` | 2, 3 | `basis`, `clamped` | `CurveExtended` |
| `curve_on_surface` | 2, 3 | `pcurve` (2D), `surface` (de `dim`) | `CurveOnSurface` |
| `offset_curve` | 2, 3 | `basis`, `function`, `direction` (3D, unitaire à 16 ε près, relue telle quelle) | `CurveOffset2D/3D<T, BSCfunction<T>>` |
| `moved_curve`, `moved_surface` | 3 | `basis` (pas elle-même déplacée), `placement` | `brep::MovedCurve`, `brep::MovedSurface` |
| `revolution_surface` | 3 | `curve` (2D), `axis` : `[o, z, x]`, `angles` | `SurfaceOfRevolution` |
| `offset_surface` | 3 | `basis`, `distance` | `OffsetSurface` |
| `seamless_surface` | 3 | `basis`, `periods` | `brep::SeamlessSurface` |

**Réels :** des nombres JSON, la forme décimale la plus courte qui se relit au bit près (Grisu
de RapidJSON), relus en pleine précision ; finis (ni NaN ni infini : ce ne serait plus du JSON).
**Texte :** UTF-8 ; `pretty` indente, un tableau de nombres sur une ligne.

**Versions.** `version` est un entier ; ce code écrit 1 et lit 1. Une version plus grande est
refusée (`UnsupportedVersion` : « version N of the format, newer than this reader (version M) »)
**avant** tout autre contrôle : un fichier plus récent peut tenir ce que ce code ne connaît pas.
Toute évolution du format (un membre, un type de géométrie, un sens nouveau) prend une version de
plus ; un lecteur plus récent lit toutes les anciennes. La version de gbs est celle de ce texte ;
l'enveloppe d'une application a la sienne.

## 5. Exactitude

### 5.1 Identifiants et entités mortes

Par défaut, les tables gardent leur taille et leurs trous : une entité morte est relue morte
(ajoutée par défaut, puis `erase`), donc `capacity`, `count`, `alive` et chaque identifiant sont
ceux du modèle écrit, **sans renumérotation** : les identifiants que l'application tient restent
justes tels quels. Le contenu d'une entité morte n'est observable par aucune fonction de `Model`
et n'est pas écrit. Les identifiants stockés sont écrits tels quels, y compris vers une entité
morte (`erase` ne propage pas) ou invalides : le modèle relu a les mêmes incohérences que
l'original. Un identifiant **au-delà** de sa table (une faute de programmation) est refusé à
l'écriture (`InvalidId`) : la lecture le refuserait.

`write_brep_compact` donne l'autre choix : le modèle compacté et son `IdRemap`, celui de
`compact()` ; les identifiants vers une entité morte y deviennent invalides, comme dans
`compact()`.

### 5.2 Partages

La géométrie est dédupliquée par l'adresse de l'objet le plus dérivé : un objet tenu par
plusieurs entités (une surface partagée par les morceaux d'une face coupée, une courbe par deux
arêtes, la géométrie partagée par `copy_shape` à placement identité) ou cité par une autre
géométrie (la surface d'une courbe sur surface, la base d'un décalage) est écrit une fois et relu
en un seul `shared_ptr`. Deux objets égaux mais distincts restent distincts. Les géométries
forment un graphe acyclique, pas un arbre ; relu, ce graphe est le même (`brep::equal` le
vérifie), donc les `use_count` comptés dans le modèle et les comparaisons de pointeurs aussi.
Les partages topologiques sont des identifiants : une arête partagée par deux co-arêtes le reste.

Les fonctions de loi (`BSCfunction` d'une courbe décalée ou reparamétrée) sont des valeurs,
écrites dans l'objet qui les tient.

### 5.3 Géométrie : les classes tenues, à l'identique

Le type dynamique **exact** est reconnu (`typeid`), pas un `dynamic_cast` : une classe dérivée
d'une B-spline n'est pas écrite comme une B-spline (elle perdrait ce qu'elle ajoute). Ce qui
oblige à des soins particuliers, parce qu'un constructeur ne rend pas toujours ce qu'on lui donne :

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
- **Direction d'une courbe décalée 3D :** le constructeur la normalise, et normaliser deux fois ne
  rend pas les mêmes bits pour un tiers des directions ((0, 1, 1) parmi elles). `CurveOffset3D` a
  donc un constructeur qui la prend telle quelle (étiquette `gbs::unit_direction`, la direction
  unitaire à 16 ε près, `is_unit`), que le relecteur appelle ; l'écrivain refuse seulement une
  direction qui n'est pas unitaire (aucun constructeur n'en fait, sauf d'un vecteur nul ou d'un
  vecteur dont la norme déborde).

**Refusées à l'écriture** (`UnsupportedEntity`, avec l'entité) : les anciennes `Line` et `Circle`
(leurs constructeurs normalisent l'axe et le rayon : relues, elles ne seraient plus les mêmes),
toute classe d'une application, les classes de calcul internes (`FieldCurve`, qui tient une
lambda, `LoopCurve`, `PathCurve`, `DefSurface`).

**Ce que gbs met dans un modèle** (relu dans le code, puis vérifié par le test `every_operation`,
§ 8) : B-splines et NURBS 2D et 3D, les sept primitives, `CurveOnSurface` (dont la surface peut
être une `SeamlessSurface`, une `OffsetSurface`, une `MovedSurface`), `CurveComposite` de
`CurveTrimmed` (2D et 3D, `make_thick`), `MovedCurve` et `MovedSurface` (`copy_shape` avec un
placement), `OffsetSurface`. Les classes de calcul n'y restent jamais, **sauf `SeamlessSurface`**
(§ 5.4). Les autres classes du tableau du § 4 ne viennent que de l'API ; elles sont tenues pour
qu'un modèle construit à la main s'enregistre aussi.

### 5.4 `SeamlessSurface` sort de `fillet.h`

La vue sans couture d'une face B-spline fermée (`fillet_detail::SeamlessSurface`) devient la
surface des courbes de contact d'un congé qui passe la couture (le congé du pied d'une pale) :
elle **reste dans le modèle**, et elle y a un sens (évaluer le contact au-delà de la couture).
Plutôt que de l'interdire dans un modèle (il faudrait changer les congés et leurs résultats),
elle devient une classe publique : elle passe telle quelle dans `gbs-brep/seamless.h`,
`gbs::brep::SeamlessSurface`, que l'écrivain connaît sans inclure le noyau des congés (4 000
lignes, 3 à 4 Go à compiler). `fillet.h` l'inclut ; ses usages non qualifiés dans `fillet_detail`
la trouvent dans l'espace de noms englobant : aucune autre ligne ne change.

### 5.5 Réels

Bit pour bit pour tout double fini, les sous-normaux et `-0` compris (vérifié sur 3,2 millions de
doubles tirés au hasard, puis par le test). Un réel non fini (NaN, infini) est refusé à l'écriture
avec l'entité qui le tient : aucun modèle valide n'en a, et le texte reste du JSON standard. En
`float`, le texte écrit le double exact du float, relu et vérifié float ; `real` dit le type, et
un modèle `float` n'est pas lu en `double` ni l'inverse.

### 5.6 Déterminisme et pureté

Le même modèle donne le même texte, octet pour octet : les entités dans l'ordre des tables, leurs
attributs avec elles (jamais dans l'ordre d'une table de hachage), les géométries dans l'ordre où
les entités les rencontrent (jamais par adresse) ; deux modèles égaux construits séparément
donnent le même texte. L'écriture est une fonction du modèle `const` seul, sans cache ni état
global : plusieurs fils peuvent écrire le même modèle pendant que d'autres le lisent.

## 6. Lire une entrée hostile

Rien n'est cru avant d'être vérifié ; aucune entrée ne mène à un comportement indéfini (le test
tourne aussi sous AddressSanitizer et UndefinedBehaviorSanitizer, § 8) :

- **Analyse** sans récursion (`kParseIterativeFlag`) : un million de crochets imbriqués est une
  faute de syntaxe, pas un dépassement de pile ; UTF-8 vérifié ; un octet NUL est refusé.
- **Schéma strict :** chaque objet n'a que ses membres, chacun une fois, les obligatoires
  présents ; chaque valeur a son type JSON. Un membre inconnu est une faute : une faute de frappe
  ne perd pas un nom en silence (un membre nouveau vient avec une version nouvelle).
- **Indices :** un identifiant est `null` ou un entier sous la taille de sa table (tailles lues
  avant les entités) ; une géométrie ne cite qu'une géométrie **avant** elle, ni elle-même : aucun
  cycle, aucune récursion à la lecture ; une référence est vérifiée de la bonne sorte
  (« a 3D curve expected, geometry 4 is a 3D surface »).
- **Profondeur :** 256 niveaux au plus (§ 4), comptés pendant la lecture (la profondeur de chaque
  géométrie se déduit de celles, déjà lues, qu'elle cite). Ce qui traverse une géométrie récurse
  une fois par niveau : l'écrivain, `brep::equal`, l'évaluation, la destruction de la chaîne de
  `shared_ptr` ; sans borne, un fichier de 93 Ko (2 000 `reversed_curve` chaînées) faisait
  déborder la pile d'un fil à l'écriture du modèle lu. L'écrivain la borne aussi, avant de
  descendre plus bas (une chaîne bâtie en mémoire est refusée sans être parcourue jusqu'au bout).
- **B-splines :** degré 24 au plus (§ 4 : au-delà, un fichier de 2 Ko bloquerait chaque
  évaluation pendant des heures) ; nœuds finis et dans l'ordre, en nombre `pôles + degré + 1`, au
  moins `degré + 1` pôles ; pôles en nombre `nu × nv` ; vérifié avant le constructeur (qui écrirait
  sur `cerr`).
- **Courbes composites :** chaque morceau de bornes croissantes ; sinon `CurveComposite::value`
  lisait hors de sa liste (`crv_lst_[-1]`), ce que la recherche bornée de son morceau empêche
  désormais aussi pour un paramètre arrondi sous son début.
- **Réels d'un modèle `float` :** dans la plage du float avant la conversion, indéfinie au-delà.
  Le reste des invariants est vérifié par les constructeurs de gbs, dont le refus devient un
  `InvalidFile` à la géométrie.
- **Attributs :** nom non vide et UTF-8 (vérifié ici aussi : `read_brep_value` reçoit un document
  analysé avec les options de l'application, peut-être sans `kParseValidateEncodingFlag`), clé de
  pièce entre 1 et `next_part - 1`, `next_part` au plus 2⁶³, unité finie et > 0.

Le modèle lu est un modèle que l'API de `Model` aurait pu construire, pas un solide valide : la
lecture ne lance pas `check`, comme `Model::add` ne le fait pas.

## 7. Égalité exacte (`gbs-brep/equal.h`)

`first_difference(a, b)` rend la première différence, décrite (« edge 12: B-spline poles
differ », « face 2: alive on one side only », « a geometry shared on one side and not the same
way on the other »), ou rien ; `equal(a, b)` dit s'il n'y en a pas. Elle compare :

- les tables : la même taille, les mêmes entités vivantes, chaque champ, chaque réel au bit près
  (`-0` et `0` diffèrent), les identifiants stockés ;
- les attributs (noms, identifiants externes, instances et leurs placements), la clé de pièce
  suivante, l'unité ;
- la géométrie : la même classe exacte, ses données de définition au bit près (pôles, nœuds,
  degrés, bornes, repères, rayons, placements, lois), et ses valeurs et dérivées en des
  paramètres échantillons, au bit près (5 par courbe jusqu'à la dérivée seconde, 3 × 3 par
  surface avec les dérivées premières) ; une classe que gbs ne connaît pas : ses bornes et ses
  valeurs seulement ;
- les partages : les paires d'objets rencontrées forment une bijection.

Un modèle et sa copie (les copies de `Model` partagent la géométrie) sont égaux ; un modèle relu
est égal au modèle écrit. Les tests de gbs s'en servent, CadSmith aussi (critère de L7s).

## 8. Validation

`tests_io_brep_json` (22 cas, 4 s ; passé aussi à la main sous AddressSanitizer et
UndefinedBehaviorSanitizer, compilé par le clang d'Apple, corrections de la revue comprises : aucune
faute — la CI ne le fait pas, et
la bibliothèque ASan du clang de conda se bloque à son démarrage sur ce macOS) compare le modèle relu au modèle écrit par `brep::equal`, qui ne passe pas par
l'écrivain ; puis le texte du modèle relu est le texte écrit, compact et indenté. Une mutation du
relecteur (oublier les bornes d'une courbe) fait échouer deux cas (« edge 1: curve bounds
differ ») : l'égalité voit ce que l'écrivain et le relecteur perdraient ensemble.

- `every_operation` (test d'acceptation 1) : dans un même modèle, les résultats de chaque
  opération : prisme, révolutions (cylindres, cône, sphère, tore), lecture STEP, tuyau, booléens
  (différence, union, intersection, faces unifiées), section, découpe, congé, chanfrein, congé à
  corde, solide décalé, solide évidé (`make_thick`), face décalée (`OffsetSurface`), copies à
  placement identité et placées (`MovedSurface`), congé du pied d'une pale (`SeamlessSurface`),
  et depuis la revue les congés qui font le plus de géométrie neuve : un coin de trois congés de
  trois rayons (congé de sommet, B-spline), d'un seul rayon (sphère), une section conique, un
  `fillet_faces` sur un chanfrein.
  Les classes attendues y sont, l'aller-retour est égal, et `check` et les volumes sont les mêmes
  au bit près : aucune classe de calcul n'y reste (l'écriture la refuserait) ;
- `every_geometry_class` : chaque classe du § 4, en 2D et en 3D, imbriquées et partagées ; la
  surface partagée par deux faces l'est encore (test d'acceptation 4), le plan écrit une fois ;
- `equality_sees_every_difference` : un bit d'un pôle, `-0`, un partage rompu, une entité morte, un
  nom, l'unité, la clé de pièce ;
- `dead_entities_ids_and_attributes`, `empty_and_absent_names` (test 3), `compacted` (test 2 :
  l'`IdRemap` est celui de `compact()`, le modèle relu égal au modèle compacté, le modèle écrit
  inchangé), `numbers_bit_for_bit`, `not_finite_and_not_utf8_refused`, `float_model` ;
- `step_files` : chacun des fichiers STEP de `tests/in/step`, puis tous dans un modèle dont un
  tiers des racines est effacé ;
- `deterministic_and_pure` (test 6) : deux écritures, deux constructions, huit fils à la fois :
  les mêmes octets ;
- `files_streams_and_values` (fichier, flux, valeur dans le document d'une application) ;
- refus (test 5) : `version_and_format`, `truncated_and_malformed` (chaque préfixe du texte, un
  octet NUL, l'imbrication profonde), `corrupt_values` (22 fautes, chacune à sa place),
  `constructors_refusals_reported`, `random_damage` (3 000 textes abîmés au hasard : lus, ou
  refusés par un `BuildError`, jamais une exception), `refused_writes` (ancienne `Line`, classe
  d'application, `Circle3d` imbriqué, composite de bornes décroissantes, degré 25, id hors table,
  compacté aussi ; le fichier en place reste intact) ;
- bornes du format (revue, § 13) : `geometry_depth_bounded` (256 niveaux lus et écrits, 257
  refusés, une chaîne de 20 000 refusée à la lecture, une de 1 000 bâtie en mémoire refusée à
  l'écriture sans être parcourue), `offset_curve_directions` (204 directions relues au bit près,
  une direction non unitaire refusée), `next_part_bounded` ; `tests_curves.composite_out_of_its_bounds`.

Python : `test_brep_json.py` (fichier, texte, compactage, pickle, égalité, erreurs).

## 9. Usage attendu par CadSmith

```cpp
// Enregistrer (dans un worker, sur l'état immuable) : le modèle dans le document de l'application
rapidjson::StringBuffer b;
rapidjson::Writer<rapidjson::StringBuffer> w(b);
w.StartObject();
w.Key("cadsmith"); /* … version de l'enveloppe, racines, corps consommés, visibilité … */
w.Key("model");
if (auto r = gbs::io::write_brep_value(model, w); !r) /* r.error() : l'entité en cause */;
w.EndObject();
gbs::write_text_atomic(path, {b.GetString(), b.GetSize()});

// Ouvrir
rapidjson::Document d;
d.Parse<gbs::io::brep_json_parse_flags>(text.data(), text.size());
auto m = gbs::io::read_brep_value<double>(d["model"], path_and_member);
// m.error().code == UnsupportedVersion -> « fichier d'une version plus récente de gbs » ;
// InvalidFile -> « fichier abîmé », le message dit où
```

Les identifiants que le document garde désignent les mêmes entités après lecture, sans
renumérotation. Pour un fichier sans entités mortes, `write_brep_compact` et l'`IdRemap` à
appliquer aux métadonnées dans le même passage. Une enveloppe en deux parties (une archive, le
texte de `write_brep` à côté des métadonnées) marche aussi.

**Arrêt coopératif (G7).** Les points où un jeton d'arrêt sera vérifié sont prêts : à chaque
entité des tables (`Writer::table`, `Reader::table`) et à chaque géométrie du réservoir
(`Writer::run`, `Reader::run`), soit quelques milliers de vérifications pour 50 000 faces. Le
mécanisme (le jeton, le code d'erreur) sera celui de G7, commun à gbs ; il n'est pas inventé ici.

## 10. Réponses aux exigences de la revue de CadSmith (G11)

| | Exigence | Réponse |
|---|---|---|
| 1 | Tout le modèle, sans exception silencieuse | Fait : chaque champ des sept tables, les attributs, `next_part`, l'unité ; ce que le format ne tient pas est refusé avec l'entité (`BuildError`). **Nom vide et nom absent :** `Model` ne tient pas de nom vide (`setName(s, "")` l'efface, `name()` rend alors `""` comme pour une forme jamais nommée) ; il n'y a rien à distinguer, et un nom vide dans un fichier est refusé (aucun modèle n'a pu l'écrire). |
| 2 | Identités ; `write` compacte et rend l'`IdRemap` | Fait autrement par défaut : les entités mortes sont gardées, **aucun identifiant ne change**, l'application n'a rien à renuméroter. `write_brep_compact` fait ce que demande la revue (le remap de `compact()`), pour qui veut un fichier sans trous. |
| 3 | Géométrie bit pour bit, évaluation comprise ; NaN et infini refusés | Fait (§ 5.3, § 5.5, § 7). Pas de périodicité dans les B-splines de gbs. |
| 4 | Tous les types que gbs-brep met dans un modèle ; jamais de classe de calcul | Fait, avec le test de chaque opération. **`SeamlessSurface` reste** dans les contacts d'un congé qui passe une couture : elle est rendue publique et sérialisée (§ 5.4) plutôt qu'interdite. `FieldCurve`, `PathCurve`, `LoopCurve` n'y restent jamais ; un type inconnu est refusé (`UnsupportedEntity`), jamais un abandon. |
| 5 | Partages restitués | Fait (§ 5.2) ; vérifié par la bijection de `brep::equal`. |
| 6 | `BuildResult`, aucune exception ; indices vérifiés ; fuzzing ; un nouveau `Model` | Fait. **Indices vers une entité vivante :** non, dans les bornes et du bon type seulement : `erase` ne propage pas, un modèle peut citer une entité morte (un composé dont un membre est effacé), et le refuser rendrait ce modèle impossible à relire. |
| 7 | Version ; plus récente refusée « version N, lecteur version M » | Fait (§ 4). |
| 8 | Déterminisme | Fait (§ 5.6), testé. |
| 9 | Fonction pure d'un modèle `const` | Fait (§ 5.6), testé à huit fils. |
| 10 | `brep::equal(a, b)` | Fait, avec `first_difference` (§ 7). |
| R | Taille et temps (50 000 faces en moins d'une seconde) | Tenu de loin (§ 3) ; écriture en flux. |
| R | Binaire ou compression en option | Non : mesures du § 3 ; la compression revient à l'enveloppe. |
| R | API qui s'imbrique | Fait : `write_brep_value` / `read_brep_value`, ou les octets de `write_brep`. |
| R | Arrêt coopératif (G7) | Points prévus (§ 9), mécanisme laissé à G7. |
| R | Spécification écrite | § 4. |

## 11. Ce qui n'est pas couvert

- Les classes refusées du § 5.3 (anciennes `Line`, `Circle`, classes d'une application). Une
  application qui en met dans un modèle doit les convertir (`LineCurve`, `EllipseCurve`) ou les
  approcher avant d'écrire.
- Les réels non finis (refusés), le contenu d'une entité morte (inobservable).
- L'enregistrement incrémental ou partiel (un sous-ensemble de formes) : on écrit le modèle entier ;
  un sous-ensemble passe par `copy_shape` dans un modèle neuf.
- La compression, l'arrêt coopératif (§ 9), le `long double`.

## 12. À valider

1. JSON plutôt que binaire (§ 3), la compression laissée à l'application.
2. Les entités mortes gardées par défaut, `write_brep_compact` au choix (§ 5.1, § 10 exigence 2).
3. Le schéma strict : un membre inconnu est une faute, toute évolution prend une version (§ 4, § 6).
4. Les refus à l'écriture (anciennes `Line` et `Circle`, classes d'application, réels non finis,
   noms non UTF-8, ids hors table) plutôt qu'une écriture approchée (§ 5.3, § 5.5).
5. `SeamlessSurface` publique dans `gbs-brep/seamless.h` plutôt qu'interdite dans un modèle (§ 5.4).
6. Les références vers une entité morte admises à la lecture (§ 10 exigence 6).
7. Les ajouts hors `gbs-io` : `Model::nextPart` et `reservePartKeys`, `BuildErrc::UnsupportedVersion`,
   `unwrap` d'un `BuildResult<void>`, `gbs-brep/equal.h` ; depuis la revue, le constructeur
   `CurveOffset3D(…, gbs::unit_direction)` et la recherche bornée de `CurveComposite::value`.
8. Les noms de l'API (`write_brep`, `read_brep`, `write_brep_value`…) et l'espace `gbs::io`, comme
   `obj.h` ; les erreurs en `BuildResult`, comme STEP.
9. Les bornes du format, ajoutées à la revue (§ 13) : profondeur de 256, degré de 24, `next_part`
   de 2⁶³ au plus, morceaux d'une composite de bornes croissantes ; refusées à l'écriture comme à
   la lecture.

## 13. Revue de code

Revue de la PR 252 (24 nov. 2026) : douze constats, de 1 à 6 et le 8 confirmés par une sonde, tous
corrigés, chacun avec son test.

1. **Profondeur des géométries non bornée** (important). Un fichier de 93 Ko (2 000
   `reversed_curve` chaînées) était lu ; écrire le modèle lu faisait déborder la pile d'un fil, et
   celle du fil principal à 20 000 niveaux. `brep::equal`, l'évaluation et, avec libstdc++ et MSVC
   (dont `vector` se détruit d'avant en arrière), la destruction de `geom_` récursaient autant.
   La profondeur est bornée à 256 (`brep_json_max_depth`, § 4, § 6), à la lecture et à
   l'écriture. Test `geometry_depth_bounded`.
2. **Composite à plages décroissantes** (important). `CurveComposite::value` lisait
   `crv_lst_[-1]` : `brep::equal` du modèle lu finissait en segfault. Le lecteur et l'écrivain
   refusent un morceau de bornes décroissantes, et `CurveComposite::value` borne l'indice de son
   morceau (un paramètre arrondi sous le début prenait aussi l'indice -1). Tests `corrupt_values`,
   `refused_writes`, `tests_curves.composite_out_of_its_bounds`.
3. **`CurveOffset3D` refusée pour un tiers des directions** (important), (0, 1, 1) parmi elles :
   normaliser deux fois change un ulp. Le constructeur `CurveOffset3D(…, gbs::unit_direction)`
   prend la direction telle quelle, unitaire à 16 ε près (`is_unit` ; une direction normalisée
   l'est à 1,5 ε au pire) ; l'écrivain ne refuse plus qu'une direction non unitaire (§ 5.3).
   Test `offset_curve_directions`.
4. **Degré illimité** (antérieur à la PR) : au-delà de 24, l'évaluation est exponentielle (0,4 s
   par point au degré 26, 6,7 s au degré 30). Le format borne le degré à 24
   (`brep_json_max_degree`), à la lecture et à l'écriture. Rendre le repli d'évaluation non
   exponentiel reste à faire, hors de cette PR. Tests `corrupt_values`, `refused_writes`.
5. `write_brep_compact` n'appliquait pas le refus `InvalidId` (l'id hors table devenait `null`) :
   `check_ids` passe avant `compact()`. Test `refused_writes`.
6. `read_brep_value` ne vérifiait pas l'UTF-8 des noms (le document de l'application peut être
   analysé sans `kParseValidateEncodingFlag`) : le lecteur le vérifie. Test
   `files_streams_and_values`.
7. `next_part` = 2⁶⁴ − 1 était accepté, et `newPart()` débordait ensuite : borné à 2⁶³, à la
   lecture et à l'écriture. Tests `corrupt_values`, `next_part_bounded`.
8. Le message d'un id hors table disait « (1 vertexs) » : le pluriel vient de `table_keys`. Test
   `refused_writes`.
9. Un nombre hors de la plage du float était converti avant d'être refusé (comportement indéfini) :
   la plage est vérifiée d'abord. Test `float_model`.
10. `every_operation` ne couvrait ni congé de sommet, ni `fillet_faces`, ni section conique : ils y
    sont (§ 8).
11. Les répertoires temporaires restaient quand un test échouait : garde `TempDir`.
12. Idiomes C++23 : `std::ranges::contains` (`is_surface_type`), `std::ranges::find` (`fields`),
    `std::countl_one` (`valid_utf8`). Les deux comparaisons au bit (`equal_detail::same`,
    `brep_json_detail::same_bits`) restent séparées : `brep_json.h` n'inclut pas `equal.h`.
