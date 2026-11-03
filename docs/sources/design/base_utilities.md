# Fonctions de base : JSON exact, écriture atomique, normale de surface, `make_range`

| | |
|---|---|
| PR | branche `feat/base-utilities` |
| Fichiers | `gbs-io/files.h` (nouveau), `gbs-io/json.h` (nouveau), `gbs-io/fromjson.h`, `gbs-io/tojson.h`, `gbs-io/step/p21.h`, `gbs/bssanalysis.h`, `gbs/surfaceoffset.h`, `gbs/math.ixx`, `python/gbsbind.cpp`, `python/gbsBindDiscretize.h`, `tests/tests_files.cpp` (nouveau), `tests/tests_bssanalysis.cpp`, `tests/tests_surfaceoffset.cpp`, `tests/tests_maths.cpp`, `python/tests/test_BSSurface.py`, `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif

Quatre fonctions de base qui manquaient à gbs, ou qu'il avait sous une forme
fausse :

1. **Lecture JSON exacte.** Un fichier écrit par `write_js_doc` et relu par
   `parse_file` ne rend pas tous ses réels : le `Writer` de RapidJSON écrit le
   plus court nombre de chiffres qui redonne le double, mais le lecteur par
   défaut en arrondit certains au voisin. Et `parse_file` ignorait les fautes
   de syntaxe : le document restait vide ou partiel, sans erreur.
2. **Écriture atomique.** `write_js_doc` écrivait sur place : un lecteur qui
   ouvre le fichier pendant l'écriture en voit une partie, et une écriture
   interrompue laisse un fichier tronqué. Il ne signalait pas non plus un
   échec.
3. **Normale d'une surface.** `gbs/surfaceoffset.h` appelle
   `normal_direction(surface, u, v)`, qui n'existait pas : le gabarit ne
   compilait pas dès qu'on l'instanciait (et `value(u)`, `bounds()` à deux
   bornes, `changeBounds` absent de `Surface`). Les normales de surface
   étaient recalculées à la main en une douzaine d'endroits de `gbs-brep`.
4. **`make_range` sans dérive.** La version 1D ajoutait son pas n − 2 fois :
   l'erreur croît le long de la suite (0,30000000000000004 au lieu de 0,3 pour
   n = 11), et ses valeurs différaient de celles des versions 2D et 3D, qui
   calculent chaque valeur à part.

## 2. Ce qui existait

Recherche sur `origin/master` à jour (56df16a) et sur les 110 branches
distantes (`git grep` sur chaque référence), dans `gbs/`, `gbs-io/`,
`gbs-mesh/`, `gbs-render/`, `gbs-brep/`, `tests/`, `python/` ; aucune PR
ouverte.

| Fonction | Trouvé | Décision |
|---|---|---|
| 1. Lecture JSON | `gbs::parse_file` (`fromjson.h`) : sans `kParseFullPrecisionFlag`, faute de syntaxe ignorée. Aucune branche n'utilise `kParseFullPrecisionFlag`, `GetParseError` ni `GetErrorOffset`. | **Étendue** : même nom, même signature (le nom de fichier devient un `std::filesystem::path`, qu'un `const char *` construit) ; `parse_json` pour un texte en mémoire. |
| 2. Écriture atomique | `gbs::write_js_doc` (`tojson.h`) : sur place. `gbs::step::detail::read_text` (`step/p21.h`) lit un fichier d'un bloc. Aucune écriture par renommage dans aucune branche (`filesystem::rename`, `std::rename`, `MoveFileEx`, `ReplaceFile` : rien). | **Étendue** : `write_js_doc` passe par `write_text_atomic` (nouvelle). `read_text` **déplacée** dans `gbs-io/files.h` et réutilisée par la lecture JSON ; `p21.h` la reprend par `using`. |
| 3. Normale | `normal_direction` des courbes 2D et 3D (`bscanalysis.h`) ; aucune pour une surface dans aucune branche. Calculs à la main de S_u × S_v dans `gbs-brep` (`offset.h`, `mesh.h`, `boolean.h`, `unify.h`, `partial.h`) et `gbs/intersect`. La règle des points dégénérés de `gbs-brep/mesh.h` (`normal_at`). | **Ajoutée** dans `bssanalysis.h`, à côté de `discretize(Surface)`, avec le seuil de `gbs-brep/mesh.h`. Les calculs de `gbs-brep` ne sont pas remplacés ici (§ 7). |
| — `discretize` avec normales | `discretize(Surface, nu, nv)` ; `discretize_with_params` pour une courbe. | **Ajoutée** : `discretize_with_normals`, sur le modèle de `discretize_with_params` ; la grille de paramètres est partagée avec `discretize`. |
| 4. `make_range` | Les versions 1D, 2D, 3D dans `gbs/math.ixx` (et le module `math.cppm`, qui l'inclut). | **Corrigée** sur place. |

## 3. API

### 3.1 `gbs-io/files.h` (nouveau, sans géométrie ni bibliothèque)

```cpp
auto read_text(const std::filesystem::path &) -> std::optional<std::string>;   // déplacée de step/p21.h
void write_text_atomic(const std::filesystem::path &, std::string_view text);  // filesystem_error
```

`write_text_atomic` écrit dans `.<nom>.<16 chiffres hexadécimaux>.tmp`, dans
le même dossier (le renommage ne copie jamais d'un système de fichiers à un
autre), puis renomme sur le fichier. Le nom vient de `std::random_device` :
deux écrivains simultanés n'ont jamais le même fichier temporaire, et le
dernier renommé l'emporte, entier. En cas d'échec, le fichier temporaire est
supprimé, le fichier n'est pas touché, et l'erreur porte la raison du système
(`errno`, ou le code du renommage).

### 3.2 `gbs-io/json.h` (nouveau ; inclus par `fromjson.h` et `tojson.h`)

```cpp
class JsonError : public std::runtime_error { std::size_t offset, line, column; std::string reason; };
void parse_json(std::string_view text, rapidjson::Document &, std::string_view source = "JSON text");
void parse_file(const std::filesystem::path &, rapidjson::Document &);           // déplacée de fromjson.h
void write_js_doc(const rapidjson::Document &, const std::filesystem::path &);   // déplacée de tojson.h
```

- Lecture avec `kParseFullPrecisionFlag`. Sur 20 000 doubles de toutes
  grandeurs, le lecteur par défaut en rend une partie au voisin ; avec le
  drapeau, tous reviennent bit pour bit (test `json_file_round_trip_is_exact`).
- `JsonError` : décalage en octets depuis 0, ligne et colonne en octets
  depuis 1, raison de RapidJSON ; `what()` = « source: line L, column C:
  raison ». La colonne compte des octets, pas des caractères UTF-8.
- `write_js_doc` écrit le même texte compact qu'avant, par
  `write_text_atomic`.
- `json.h` ne tire ni les courbes, ni les surfaces, ni Boost.UUID, que
  `tojson.h` inclut : un programme qui ne lit que des paramètres l'inclut
  seul.

### 3.3 `gbs/bssanalysis.h`

```cpp
template <typename T> auto normal_direction(const Surface<T, 3> &, T u, T v) -> point<T, 3>;  // std::domain_error
template <typename T> auto discretize_with_normals(const Surface<T, 3> &, size_t nu, size_t nv); // (points, normales)
```

Python : `gbs.normal_direction(srf, u, v)` (surcharge de celle des courbes),
`gbs.discretize_surface_with_normals(srf, nu, nv)` → deux tableaux (nu nv, 3).
`std::domain_error` devient `ValueError`.

### 3.4 `gbs/surfaceoffset.h`

`SurfaceOffset<T, 3, Func>` compile et marche : `value(u, v)` =
S + f(u, v) n, `bounds()` à quatre bornes. Les dérivées lèvent
`std::runtime_error` (elles le faisaient déjà, en commentaire) ;
`has_second_derivatives()` est faux. `changeBounds`, qui appelait une méthode
que `Surface` n'a pas, est retirée. Un `static_assert` limite la classe à la
3D, la seule où une surface a une normale.

### 3.5 `gbs/math.ixx` : `make_range(v1, v2, n)`

Chaque valeur est `v1 + (v2 − v1) i / (n − 1)`, arrondie comme dans les
versions 2D et 3D ; v1 et v2 sont rendus exactement. Un point est traité
coordonnée par coordonnée (avant : le pas d'un point était multiplié par un
réel, puis ajouté). Les versions 2D et 3D donnent maintenant les mêmes
valeurs intérieures que la 1D (test `make_range_without_drift`).

## 4. Choix

### 4.1 La normale aux points dégénérés

Là où S_u × S_v s'annule devant max(|S_u|, |S_v|)² (seuil 1e-6, celui de
`gbs-brep/mesh.h`), `normal_direction` rend **la limite de la normale depuis
l'intérieur du domaine**. Elle prend la normale en trois points, à 1, 2 et 3
pas de 1e-4 du chemin de (u, v) au centre des bornes, et extrapole la
parabole qui passe par eux jusqu'à (u, v). L'erreur est en h³.

Mesures :

- au pôle d'une sphère, primitive ou NURBS à rangée de pôles confondus :
  5e-11 ;
- ailleurs : 1e-15 ;
- au sommet d'un cône : la normale de la génératrice du côté du centre.

Pourquoi pas les autres comportements :

- **Lever une exception, ou rendre NaN ou l'infini** comme la normale 3D des
  courbes : l'offset d'une sphère, d'une surface de révolution fermée sur son
  axe ou d'une surface B-spline à pôles confondus serait faux ou planterait en
  ses pôles, alors que la surface y a une normale. Seule la paramétrisation est
  dégénérée.
- **Un seul point décalé**, comme le maillage de `gbs-brep` : l'erreur serait
  de l'ordre du pas. Il faudrait 1e-4 pour rester au-dessus du seuil, soit une
  normale fausse de 1e-4 rad, trop pour un offset.
- **Les dérivées secondes** (n ∝ S_uv × S_v au pôle d'une ligne v = cste) : ne
  couvre qu'un cas, et suppose `has_second_derivatives()`.

Le côté du centre est arbitraire au sommet d'un cône, où la surface n'a pas de
normale : c'est documenté. Si les trois points sont eux aussi dégénérés
(surface réduite à un point, paramètres non finis), la fonction lève
`std::domain_error`.

### 4.2 Windows et le renommage

Sur Windows, `std::filesystem::rename` (`MoveFileExW`, avec
`MOVEFILE_REPLACE_EXISTING`) échoue si un autre processus tient le fichier
ouvert sans `FILE_SHARE_DELETE`. C'est le cas d'un lecteur qui l'ouvre avec
`fopen` ou `std::ifstream`, comme le font les bibliothèques C et C++ de
Microsoft. L'erreur est « accès refusé » ou « violation de partage ».
`write_text_atomic` lève alors `filesystem_error`, et le fichier garde son
ancien texte entier. Elle ne réessaie pas : c'est à l'appelant de savoir si un
lecteur peut tenir le fichier, et combien de temps attendre (point 2).

### 4.3 Atomique, pas durable

Rien n'est forcé sur le disque (pas de `fsync` du fichier ni du dossier) : une
coupure de courant juste après peut perdre le nouveau texte. Sur un système de
fichiers journalisé, elle ne mélange pas l'ancien et le nouveau (point 3).

### 4.4 Parallélisme et exceptions

`discretize_with_normals` calcule les normales par `transform_rethrow`
au-delà de `parallel_min_size`, et en série en dessous. Une `domain_error` ne
sort jamais d'un algorithme parallèle (sinon `std::terminate`) : le test
`discretize_with_normals` la vérifie sur une grille de 100 × 100.

## 5. Changements de comportement

- `parse_file` lève `JsonError` sur une faute de syntaxe, au lieu de laisser
  un document vide ou partiel. Il lit le fichier en binaire : sous Windows,
  les fins de ligne `\r\n` ne sont plus converties, ce que RapidJSON prend
  comme des espaces.
- `write_js_doc` lève `filesystem_error` si le fichier ne peut pas être écrit
  (avant : rien), et laisse un fichier temporaire le temps de l'écriture.
- `make_range(v1, v2, n)` : les valeurs intérieures changent. Mesuré sur
  [0, 1], de n = 3 à 1000, contre l'ancienne somme :

  | n | Valeurs intérieures qui changent | Écart maximal |
  |---|---|---|
  | 3 et 5 | aucune | — |
  | 11 | 3 sur 9 | 1 ULP |
  | 100 | 85 sur 98 | 18 ULP |
  | 1000 | — | 241 ULP, soit 2,7e-14 |

  Sur 2π, l'écart va jusqu'à 1,5e-13. Ce sont les erreurs de l'ancienne
  somme ; les nouvelles valeurs sont arrondies une seule fois. Une suite
  symétrique, de −1 à 1 par exemple, passe maintenant exactement par 0.

  Tout ce qui échantillonne par `make_range` se déplace d'autant :
  `discretize` de courbes et de surfaces, `abs_curv`, `offsets.h`,
  `bssbuild.h`, `loftBase.h`, les maillages de `gbs-mesh`.

## 6. Tests

| Test | Ce qu'il vérifie |
|---|---|
| `tests_files` (nouveau) | 20 000 doubles relus bit pour bit, alors que le lecteur par défaut en manque ; une courbe B-spline écrite et relue à l'identique ; la position et la raison d'une faute de syntaxe, un fichier absent, un texte vide ; `read_text` sur des octets, un fichier absent, un dossier. |
| `tests_files`, écriture atomique | Remplacement, texte vide, 1 Mio, aucun fichier temporaire laissé. Échec du renommage sur un dossier non vide (le dossier intact, le temporaire supprimé). Échec de création dans un dossier absent. Quatre écrivains simultanés : le fichier final est entier. |
| `tests_maths.make_range_without_drift` | i / (n − 1) arrondi une fois, en `double` et en `float`, pour un point, égal aux versions 2D et 3D ; 0,3 pour n = 11 ; une suite de 4001 valeurs sur [−1000, 1000]. |
| `tests_bssanalysis` | La normale d'une sphère tournée, primitive et NURBS, à ses pôles compris (1e-9) ; un plan ; une surface réduite à un point et un paramètre NaN (`domain_error`) ; le sommet d'un cône ; `discretize_with_normals` (mêmes points que `discretize`, normales de `normal_direction`, exception propagée en série et en parallèle). |
| `tests_surfaceoffset` | L'offset d'une sphère est une sphère, pôles compris. |
| `python/tests/test_BSSurface.py` | `normal_direction` et `discretize_surface_with_normals` depuis Python. |

Les tests qui comparent à des valeurs de référence ont été repassés avec le
nouveau `make_range` : `ctest` et `pytest` passent tous, sans changer une
tolérance ni une référence (§ 6 bis).

### 6 bis. Tests à référence et `make_range`

Suite complète sur macOS (clang 21.1, Release, environnement conda-forge),
avant et après :

| | `master` (56df16a) | Cette branche |
|---|---|---|
| `ctest` | 682 / 682 | 695 / 695 (13 nouveaux) |
| `pytest python/tests` | — | 90 réussis, 1 `xfail` déjà marqué |

Aucune référence ni tolérance n'a été changée. Les tests qui passent par
`make_range` ou `discretize` et comparent par égalité ou à 4 ULP
(`ASSERT_EQ`, `ASSERT_DOUBLE_EQ`) sont dans `tests_basisfunctions`,
`tests_bscanalysis`, `tests_bscapprox`, `tests_bscurve`, `tests_bsinterp`,
`tests_bssanalysis`, `tests_bssbuild`, `tests_bssurf`, `tests_curves`,
`tests_knotsfunctions`, `tests_maths` et `tests_mesh`. Tous passent. Les tests
Python de `discretize` passent aussi.

## 7. Hors périmètre

Remarqué, non traité ici :

- Les calculs à la main de S_u × S_v dans `gbs-brep` et `gbs/intersect`
  pourraient appeler `normal_direction`. Mais chacun a sa propre règle aux
  points dégénérés (un `optional`, un repli sur le triangle, rien), et la
  changer déplacerait des résultats de booléens. C'est une PR à part, si l'on
  veut.
- `discretize(Surface)` évalue par `transform_threshold` : une surface qui
  lève en un point termine le processus au-delà de `parallel_min_size`.
  `transform_rethrow` le corrigerait, au prix d'un passage parallèle même pour
  les petites grilles.
- `python/gbsbind.cpp` lie `tangential_direction` à `normal_direction` pour
  les courbes : en Python, `tangential_direction` rend la normale.
- `mesh_writer.h` et les écrivains STEP et IGES écrivent aussi sur place ; ils
  pourraient passer par `write_text_atomic`.

## 8. Points à valider

1. **Le comportement aux points dégénérés** (§ 4.1) : la limite depuis
   l'intérieur, extrapolée, plutôt qu'une exception. Recommandation : la
   garder, c'est ce dont un offset a besoin.
2. **Pas de nouvel essai sous Windows** (§ 4.2) quand un lecteur tient le
   fichier. Recommandation : laisser l'appelant décider. Un nombre d'essais en
   option peut s'ajouter plus tard sans casser l'API.
3. **Pas de `fsync`** (§ 4.3). Recommandation : en rester là pour des
   fichiers de travail ; une option `durable` le jour où gbs écrit des
   fichiers qu'une coupure ne doit pas perdre.
4. **`parse_file` lève sur une faute de syntaxe** : c'est un changement de
   comportement. Recommandation : le garder, l'ancien comportement cachait
   l'erreur jusqu'au premier accès au document.
5. **Les noms** : `parse_json`, `JsonError`, `read_text`,
   `write_text_atomic`, `discretize_with_normals`. Ils suivent
   `discretize_with_params` et `read_p21` ; `read_text` passe de
   `gbs::step::detail` à `gbs`.
6. **`make_range` change des résultats** (§ 5) au-delà d'un ULP pour les
   grands n. Recommandation : l'accepter, c'est l'erreur de l'ancienne somme
   qui disparaît, et c'est signalé dans `News.md`.
7. **Les extrémités des versions 2D et 3D** ne sont pas forcées à v2 comme en
   1D : (v2 − v1)(n − 1)/(n − 1) peut s'écarter de v2 − v1 d'un ULP.
   Recommandation : les forcer aussi, dans une PR à part, car cela change
   d'autres résultats.
