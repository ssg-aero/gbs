# Fichiers OBJ de maillages polygonaux

| | |
|---|---|
| PR | branche `feat/io-obj` |
| Fichiers | `gbs-io/obj.h` (nouveau), `gbs-io/real_format.h` (`parse_real`, venu de `gbs-io/step/p21.h`), `gbs-io/step/p21.h`, `gbs-io/files.h` (`read_text_or_throw`, `TextError`), `gbs-io/json.h`, `python/gbsBindBrepOps.cpp`, `tests/tests_io_obj.cpp` (nouveau), `python/tests/test_brep_mesh.py`, `News.md` |
| Statut | À valider ; revue de code faite (§ 7) |

## 1. Objectif

Lire et écrire un **maillage polygonal** au format Wavefront OBJ, celui qu'exportent
les modeleurs (Blender, Maya, Rhino) : des sommets, des faces d'un nombre
quelconque de côtés, et les deux autres éléments du format qui ne renvoient qu'à
des sommets, les lignes (`l`) et les points (`p`).

gbs-io écrivait déjà de l'OBJ (`mesh_writer.h`), mais seulement les triangles
d'un maillage du BREP, avec leurs normales, et ne lisait rien. Un modèle
polygonal fait dans un modeleur (un maillage de quadrangles ou de
polygones) ne pouvait ni entrer dans gbs ni en sortir sans être triangulé.

## 2. API

```cpp
namespace gbs::io
{
    template <typename T> concept obj_real = std::same_as<T, float> || std::same_as<T, double>;

    struct Polygons // à plat : les sommets du k-ième sont indices[offsets[k], offsets[k + 1])
    {
        std::vector<std::uint32_t> indices;
        std::vector<std::size_t> offsets{0};
        // size(), empty(), operator[](k) -> std::span<const std::uint32_t>, push_back(range),
        // begin() / end() (une plage de spans), ==, construction depuis des listes
    };
    template <obj_real T> struct PolygonMesh
    {
        points_vector<T, 3> vertices;
        Polygons faces, lines;             // indices depuis 0
        std::vector<std::uint32_t> points;
    };
    template <obj_real T> struct ObjRead { PolygonMesh<T> mesh; std::vector<std::pair<std::string, std::size_t>> skipped; };
    class ObjError : public TextError { std::size_t line; }; // TextError (files.h) : reason()

    template <obj_real T = double> auto read_obj(std::string_view text, std::string_view source = "OBJ") -> ObjRead<T>;
    template <obj_real T = double> auto read_obj_file(const std::filesystem::path &) -> ObjRead<T>;
    template <obj_real T> auto write_obj(const PolygonMesh<T> &, std::string_view comment = {}) -> std::string;
    template <obj_real T> void write_obj_file(const PolygonMesh<T> &, const std::filesystem::path &, std::string_view comment = {});
}
```

Les faces et les lignes sont rangées **à plat** (`Polygons`) : une allocation
pour toutes, au lieu d'une par face. Un maillage de quadrangles coûte 16 octets
d'indices et 8 de début par face, contre environ 80 avec un vecteur par face.
Les indices sont sur 32 bits, comme ceux des maillages du BREP (`brep::Mesh`) :
un maillage tient au plus 4 294 967 295 sommets. T est `float` ou `double`, les
types que `parse_real` et `append_real` lisent et écrivent exactement.

Python (`gbs.brep`, à côté de `write_stl` et `write_obj` des maillages du BREP) :
`PolygonMesh(vertices, faces, lines, points)`, `read_obj(path)` et
`read_obj_text(text)` → `(PolygonMesh, [(mot-clé, nombre)])`, `write_obj(mesh,
path, comment="")`, `write_obj_text(mesh, comment="")` ; une faute de lecture
est une `brep.ObjError`, sous-classe de `ValueError`, avec `line` et `reason` ;
un fichier illisible est une `OSError`. Les attributs d'un `PolygonMesh` sont
des copies (des listes, un tableau numpy) : on change un maillage en assignant
une nouvelle valeur entière (`mesh.faces = faces`), pas en modifiant un élément
de la copie, et la docstring le dit. La lecture et l'écriture se font sans le
GIL, et le maillage lu passe à Python sans copie.

`write_obj` et `write_obj_file` sont des surcharges de celles de
`mesh_writer.h` (un autre type de maillage) ; elles lèvent, comme `json.h` et
`files.h`, au lieu de rendre un `BuildResult` : rien du BREP ici.

## 3. Lecture

- `v x y z` : un sommet ; jusqu'à quatre nombres après z (un poids, une
  couleur que certains exportateurs ajoutent) sont admis et non gardés ; tout
  autre mot est une faute. Les coordonnées sont lues au type T, en un seul
  arrondi : en `float`, 1e39 est refusé comme non fini.
- `f` (et `fo`, ancienne forme), `l`, `p` : des références `i`, `i/t`, `i//n` ou
  `i/t/n`, dont seul i est gardé ; depuis 1, ou négatives (relatives aux sommets
  déjà lus : −1 le dernier). Une référence à un sommet défini plus loin est
  admise, et vérifiée à la fin.
- `#` commence un commentaire ; une ligne qui finit par `\` continue sur la
  suivante ; une ligne finit par `\n`, `\r\n` ou un `\r` seul (les fichiers de
  l'ancien Mac OS, qui sinon ne formaient qu'une ligne) ; une marque d'ordre
  des octets UTF-8 en tête (certains outils Windows) est sautée.
- Tout autre énoncé (`vt`, `vn`, `g`, `o`, `s`, `usemtl`, `mtllib`, les courbes
  et surfaces de forme libre…) est sauté et **compté par mot-clé**
  (`ObjRead::skipped`) : rien n'est perdu sans le dire.
- Refusés avec `ObjError` et la ligne (la première d'un énoncé continué) : une
  face de moins de 3 sommets, une ligne de moins de 2, un élément `p` vide, une
  référence nulle, hors des sommets, au-delà de 2³² − 1 ou mal formée (les
  extrêmes des entiers compris, sans débordement), une coordonnée absente, mal
  formée ou non finie en T, plus de quatre nombres après z.
- Les réels sont lus par `parse_real`, celle du lecteur STEP, indépendante de la
  locale et identique sur toutes les plateformes (`strtod_l` ou `strtof_l` dans
  la locale « C » avec la libc++ de conda-forge, `from_chars` ailleurs), en
  `float` ou en `double`. Comme `from_chars`, elle refuse un réel hexadécimal et
  les blancs de tête, que `strtod` lit. Elle passe de `step/p21.h` à
  `real_format.h`, à côté de `append_real` ; le lecteur STEP la reprend par une
  déclaration `using`.
- `read_obj_file` lève par `read_text_or_throw` (`files.h`), partagée avec
  `parse_file` de `json.h` ; `ObjError` et `JsonError` dérivent de `TextError`
  (`files.h`), qui porte `reason()`.

## 4. Écriture

Les lignes de `comment` en commentaires, puis les sommets (`v`, chaque
coordonnée sous la plus courte forme qui se relit exactement **en T**,
`append_real` : 0.1f s'écrit `0.1`, pas les 17 chiffres du double qu'il est),
puis `f`, `l` et un seul `p`, depuis 1. Ce que la lecture refuserait (un
élément trop court, une référence hors des sommets, un sommet non fini) lève
`std::invalid_argument` avant d'écrire quoi que ce soit. Le fichier s'écrit
entier ou pas du tout (`write_text_atomic`).

## 5. Validation

`tests/tests_io_obj.cpp`, 8 cas :

| Test | Contenu |
|---|---|
| `a_modeller_file` | un fichier comme ceux de Blender (objet, groupes, normales, coordonnées de texture, matériaux) : triangles, quadrangle, références relatives, poids et couleur ignorés, lignes et points ; les énoncés sautés comptés dans l'ordre |
| `round_trip_exact` | réels au hasard sur six ordres de grandeur, un dénormal, −0, 1/3, un pentagone : relu égal bit pour bit, réécrit à l'identique ; en `float` aussi, sous la forme la plus courte d'un float (2000 tirages) |
| `whatever_the_locale` | dans une locale à virgule décimale (fr_FR ou de_DE quand le système l'a), le texte n'a pas de virgule et se relit ; sans une telle locale, un avertissement le dit ; la locale est remise même après un échec |
| `continuations_comments_and_line_ends` | `\r\n`, commentaires en fin de ligne, lignes blanches, `\` ; une référence en avant ; un texte vide |
| `byte_order_mark_and_line_ends` | une marque d'ordre des octets sautée (4 sommets, la face juste) ; des `\r` seuls en fin de ligne, et les numéros de ligne avec eux |
| `polygons` | `Polygons` : construit depuis des listes, indexé, parcouru en spans, comparé, vidé |
| `faults_named_by_their_line` | chaque refus, son message et sa ligne ; des mots après z, un réel hexadécimal, 1e39 en `float`, les références extrêmes (−2⁶³, 2⁶³ − 1, au-delà de 2⁶⁴) ; ce que la lecture refuse n'est pas écrit |
| `files` | aller-retour par fichier ; fichier absent ou dossier : `filesystem_error` ; le nom du fichier dans le message d'une faute ; le dossier temporaire retiré même après un échec |

`python/tests/test_brep_mesh.py::test_polygon_mesh_obj` : la même chose par
Python, plus `ObjError` (`line`, `reason`) et l'assignation des attributs. Les tests du lecteur et de l'écrivain STEP passent sans changement
(`tests_step_p21`, `tests_step_read`, `tests_step_write`, `tests_step_p21_writer`).

## 6. À valider

1. **Le maillage polygonal en type de gbs-io** (`PolygonMesh`), à part du
   maillage triangulé du BREP (`brep::Mesh`) : il sert à échanger, pas à
   calculer ; une conversion vers `gbs-mesh` viendrait si on en a besoin.
2. **Lignes et points gardés**, les autres éléments sautés et comptés :
   `vt`, `vn` (attributs des coins) ne servent pas à la géométrie d'un maillage
   polygonal ; les courbes et surfaces de forme libre de l'OBJ sont rares, et
   gbs a STEP et IGES pour elles.
3. **Exceptions** plutôt que `BuildResult`, comme `json.h`.
4. **`parse_real` dans `real_format.h`** (§ 3).

## 7. Revue de code

Quatorze remarques, toutes traitées dans cette PR.

| # | Remarque | Suite |
|---|---|---|
| 1 | Une marque d'ordre des octets UTF-8 fait perdre le premier sommet et décale toutes les faces, sans erreur. | Sautée en tête du texte (§ 3). Test `byte_order_mark_and_line_ends`. |
| 2 | Les mots après z ne sont pas lus : un fichier aux fins de ligne `\r` seules ne donne qu'un sommet, `v 1 2 3 oops` passe. | `\r` seul fin de ligne ; après z, quatre nombres au plus, rien d'autre (§ 3). Mêmes tests. |
| 3 | `read_obj<float>` accepte 1e39 (un sommet infini, une conversion indéfinie) et lit en deux arrondis. | Lecture au type T par `parse_real` en `float` (§ 3). Test `faults_named_by_their_line`. |
| 4 | `-LLONG_MIN` déborde (comportement indéfini). | Comparaison sans négation ; références au-delà de 2³² − 1 refusées. Même test. |
| 5 | `parse_real` lit l'hexadécimal avec la libc++ de conda-forge, pas ailleurs. | Refusé comme par `from_chars`, avec les blancs de tête (§ 3). Même test. |
| 6 | En Python, modifier un élément d'un attribut est sans effet. | Attributs documentés comme des copies à réassigner ; test de l'assignation. |
| 7 | `write_obj` écrit tout T comme un double : 17 chiffres pour un float, `long double` inexact. | `append_real` d'un float, la plus courte en float ; T restreint à `float` et `double` (§ 2, § 4). Test `round_trip_exact`. |
| 8 | En Python, le maillage lu est copié, et le GIL gardé pendant les entrées-sorties. | Maillage déplacé ; lecture et écriture sans le GIL. |
| 9 | Une allocation par face. | `Polygons` à plat, indices sur 32 bits (§ 2). Test `polygons`. |
| 10 | `read_obj_file` recopie l'erreur de fichier de `parse_file`. | `read_text_or_throw` dans `files.h`, pour les deux. |
| 11 | `ObjError` recopie `reason()` de `JsonError`. | Base commune `TextError` dans `files.h`. |
| 12 | Sentinelle −1 et raison en paramètre de sortie. | `std::expected<std::uint32_t, std::string>` ; les lignes de commentaire par `std::views::split`. |
| 13 | Tests sans ces cas ; test de locale muet sans locale à virgule ; dossier temporaire laissé après un échec. | Cas ajoutés ; avertissement ; gardes de locale et de dossier (RAII). |
| 14 | En Python, `ValueError` sans la ligne. | `brep.ObjError` (sous-classe de `ValueError`) avec `line` et `reason`. |
