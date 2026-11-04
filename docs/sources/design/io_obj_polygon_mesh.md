# Fichiers OBJ de maillages polygonaux

| | |
|---|---|
| PR | branche `feat/io-obj` |
| Fichiers | `gbs-io/obj.h` (nouveau), `gbs-io/real_format.h` (`parse_real`, venu de `gbs-io/step/p21.h`), `gbs-io/step/p21.h`, `python/gbsBindBrepOps.cpp`, `tests/tests_io_obj.cpp` (nouveau), `python/tests/test_brep_mesh.py`, `News.md` |
| Statut | À valider |

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
    template <std::floating_point T> struct PolygonMesh
    {
        points_vector<T, 3> vertices;
        std::vector<std::vector<std::size_t>> faces, lines; // indices depuis 0
        std::vector<std::size_t> points;
    };
    template <std::floating_point T> struct ObjRead { PolygonMesh<T> mesh; std::vector<std::pair<std::string, std::size_t>> skipped; };
    class ObjError : public std::runtime_error { std::size_t line; std::string_view reason() const; };

    template <std::floating_point T = double> auto read_obj(std::string_view text, std::string_view source = "OBJ") -> ObjRead<T>;
    template <std::floating_point T = double> auto read_obj_file(const std::filesystem::path &) -> ObjRead<T>;
    template <std::floating_point T> auto write_obj(const PolygonMesh<T> &, std::string_view comment = {}) -> std::string;
    template <std::floating_point T> void write_obj_file(const PolygonMesh<T> &, const std::filesystem::path &, std::string_view comment = {});
}
```

Python (`gbs.brep`, à côté de `write_stl` et `write_obj` des maillages du BREP) :
`PolygonMesh(vertices, faces, lines, points)`, `read_obj(path)` et
`read_obj_text(text)` → `(PolygonMesh, [(mot-clé, nombre)])`, `write_obj(mesh,
path, comment="")`, `write_obj_text(mesh, comment="")` ; une faute de lecture
est une `ValueError` qui nomme la ligne, un fichier illisible une `OSError`.

`write_obj` et `write_obj_file` sont des surcharges de celles de
`mesh_writer.h` (un autre type de maillage) ; elles lèvent, comme `json.h` et
`files.h`, au lieu de rendre un `BuildResult` : rien du BREP ici.

## 3. Lecture

- `v x y z` : un sommet ; des nombres après z (un poids, une couleur que
  certains exportateurs ajoutent) sont admis et non gardés.
- `f` (et `fo`, ancienne forme), `l`, `p` : des références `i`, `i/t`, `i//n` ou
  `i/t/n`, dont seul i est gardé ; depuis 1, ou négatives (relatives aux sommets
  déjà lus : −1 le dernier). Une référence à un sommet défini plus loin est
  admise, et vérifiée à la fin.
- `#` commence un commentaire ; une ligne qui finit par `\` continue sur la
  suivante ; les fins de ligne `\r\n` sont admises.
- Tout autre énoncé (`vt`, `vn`, `g`, `o`, `s`, `usemtl`, `mtllib`, les courbes
  et surfaces de forme libre…) est sauté et **compté par mot-clé**
  (`ObjRead::skipped`) : rien n'est perdu sans le dire.
- Refusés avec `ObjError` et la ligne (la première d'un énoncé continué) : une
  face de moins de 3 sommets, une ligne de moins de 2, un élément `p` vide, une
  référence nulle, hors des sommets ou mal formée, une coordonnée absente, mal
  formée ou non finie.
- Les réels sont lus par `parse_real`, celle du lecteur STEP, indépendante de la
  locale et identique sur toutes les plateformes (`strtod_l` dans la locale « C »
  avec la libc++ de conda-forge, `from_chars` ailleurs). Elle passe de
  `step/p21.h` à `real_format.h`, à côté de `append_real` ; le lecteur STEP la
  reprend par une déclaration `using`, sans changement.

## 4. Écriture

Les lignes de `comment` en commentaires, puis les sommets (`v`, chaque
coordonnée sous la plus courte forme qui se relit exactement, `append_real`),
puis `f`, `l` et un seul `p`, depuis 1. Ce que la lecture refuserait (un
élément trop court, une référence hors des sommets, un sommet non fini) lève
`std::invalid_argument` avant d'écrire quoi que ce soit. Le fichier s'écrit
entier ou pas du tout (`write_text_atomic`).

## 5. Validation

`tests/tests_io_obj.cpp`, 6 cas :

| Test | Contenu |
|---|---|
| `a_modeller_file` | un fichier comme ceux de Blender (objet, groupes, normales, coordonnées de texture, matériaux) : triangles, quadrangle, références relatives, poids et couleur ignorés, lignes et points ; les énoncés sautés comptés dans l'ordre |
| `round_trip_exact` | réels au hasard sur six ordres de grandeur, un dénormal, −0, 1/3, un pentagone : relu égal bit pour bit, réécrit à l'identique ; en `float` aussi |
| `whatever_the_locale` | dans une locale à virgule décimale (fr_FR ou de_DE quand le système l'a), le texte n'a pas de virgule et se relit |
| `continuations_comments_and_line_ends` | `\r\n`, commentaires en fin de ligne, lignes blanches, `\` ; une référence en avant ; un texte vide |
| `faults_named_by_their_line` | chaque refus, son message et sa ligne ; ce que la lecture refuse n'est pas écrit |
| `files` | aller-retour par fichier ; fichier absent ou dossier : `filesystem_error` ; le nom du fichier dans le message d'une faute |

`python/tests/test_brep_mesh.py::test_polygon_mesh_obj` : la même chose par
Python. Les tests du lecteur et de l'écrivain STEP passent sans changement
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
