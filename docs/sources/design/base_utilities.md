# Fonctions de base : JSON exact, écriture atomique, normale de surface, `make_range`

| | |
|---|---|
| PR | branche `feat/base-utilities` |
| Fichiers | `gbs-io/files.h` (nouveau), `gbs-io/json.h` (nouveau), `gbs-io/fromjson.h`, `gbs-io/tojson.h`, `gbs-io/step/p21.h`, `gbs/bssanalysis.h`, `gbs/surfaceoffset.h` (retiré), `gbs/surfaces`, `gbs-brep/offset.h`, `gbs/math.ixx`, `python/gbsbind.cpp`, `python/gbsBindDiscretize.h`, `tests/tests_files.cpp` (nouveau), `tests/tests_bssanalysis.cpp`, `tests/tests_surfaceoffset.cpp` (retiré), `tests/tests_brep_sweep.cpp`, `tests/tests_maths.cpp`, `python/tests/test_BSSurface.py`, `docs/index.rst`, `News.md` |
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
3. **Normale d'une surface.** `gbs/surfaceoffset.h` appelait
   `normal_direction(surface, u, v)`, qui n'existait pas : le gabarit ne
   compilait pas dès qu'on l'instanciait (et `value(u)`, `bounds()` à deux
   bornes, `changeBounds` absent de `Surface`). La surface décalée de
   `gbs-brep` (`OffsetSurface`) divisait par |S_u × S_v| : aux pôles d'une
   B-spline à rangée de pôles confondus, ce produit est de l'ordre de
   l'arrondi, et sa direction est du bruit (des points jusqu'à 0,23 à côté sur
   une sphère de rayon 2 décalée de 0,5 ; NaN s'il s'annule). Les normales de surface
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
| 3. Normale | `normal_direction` des courbes 2D et 3D (`bscanalysis.h`) ; aucune pour une surface dans aucune branche. Calculs à la main de S_u × S_v dans `gbs-brep` (`offset.h`, `mesh.h`, `boolean.h`, `unify.h`, `partial.h`) et `gbs/intersect`. La règle des points dégénérés de `gbs-brep/mesh.h` (`normal_at`). | **Ajoutée** dans `bssanalysis.h`, à côté de `discretize(Surface)`, avec la règle de `gbs-brep/mesh.h` appliquée aux dérivées multipliées par les longueurs des bornes (§ 4.1). `OffsetSurface` l'utilise ; `SurfaceOffset` est retirée (§ 3.4). Les autres calculs de `gbs-brep` ne sont pas remplacés ici (§ 7). |
| — `discretize` avec normales | `discretize(Surface, nu, nv)` ; `discretize_with_params` pour une courbe. | **Ajoutée** : `discretize_with_normals`, sur le modèle de `discretize_with_params` ; la grille de paramètres est partagée avec `discretize`, et une seule boucle donne le point et la normale de chaque nœud. |
| 4. `make_range` | Les versions 1D, 2D, 3D dans `gbs/math.ixx` (et le module `math.cppm`, qui l'inclut). | **Corrigée** sur place. |

## 3. API

### 3.1 `gbs-io/files.h` (nouveau, sans géométrie ni bibliothèque)

```cpp
auto read_text(const std::filesystem::path &) -> std::optional<std::string>;   // déplacée de step/p21.h
void write_text_atomic(const std::filesystem::path &, std::string_view text);  // filesystem_error
auto utf8_string(const std::filesystem::path &) -> std::string;                // un chemin pour un message
```

`write_text_atomic` écrit dans `.<nom>.<16 chiffres hexadécimaux>.tmp`, dans
le même dossier (le renommage ne copie jamais d'un système de fichiers à un
autre), puis renomme sur le fichier. Le nom vient de `std::random_device` :
deux écrivains simultanés n'ont jamais le même fichier temporaire, et le
dernier renommé l'emporte, entier. Le nom du fichier y est coupé à 200 unités
(sans couper un caractère UTF-8 ou UTF-16), pour que le tout tienne dans les
255 d'un nom. En cas d'échec, le fichier temporaire est supprimé, le fichier
n'est pas touché, et l'erreur porte la raison du système (`errno`, ou le code
du renommage).

Comme une écriture sur place :

- un lien symbolique est suivi : le lien reste, sa cible reçoit le texte ;
- le fichier garde ses droits (copiés sur le temporaire avant le renommage) ;
- un fichier qu'on ne peut pas écrire (lecture seule) est refusé : il est
  ouvert en ajout, sans rien écrire, avant tout.

Différences qui restent : un lien physique est rompu (les autres noms du
fichier gardent l'ancien texte) ; le propriétaire et les listes de contrôle
d'accès sont ceux de l'écrivain.

Le nom temporaire est construit sur les chaînes natives des chemins, et les
messages passent par `utf8_string` : sous Windows, `path::string()` convertit
dans la page de code du système et lève pour un nom qu'elle ne contient pas.

### 3.2 `gbs-io/json.h` (nouveau ; inclus par `fromjson.h` et `tojson.h`)

```cpp
class JsonError : public std::runtime_error { std::size_t offset, line, column; auto reason() const noexcept -> std::string_view; };
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
  `reason()` est la fin de `what()` : copier un `JsonError` ne lève jamais.
- Un octet nul est une erreur : RapidJSON arrête le texte là, ce qui suit ne
  serait pas lu. Après une erreur, le document est nul.
- `parse_file` lève `std::filesystem::filesystem_error` (un
  `std::runtime_error`) avec la raison du système quand le fichier ne peut pas
  être lu : absent, dossier, droits.
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

`discretize` et `discretize_with_normals` passent par `build_batch`
(`gbs/execution.h`) avec le seuil `parallel_min_size` : en parallèle au-delà,
en série en deçà, et la première exception, dans l'ordre de la grille, est
relancée après la boucle (§ 4.4).

### 3.4 Surfaces décalées : `OffsetSurface`, et `SurfaceOffset` retirée

`gbs::brep::OffsetSurface` (`gbs-brep/offset.h`, celle de `offset_face`) fait
déjà S + d n, avec les dérivées premières par différences finies. Elle prend
maintenant sa normale de `normal_direction` : aux pôles d'une B-spline à
rangée de pôles confondus, la limite de la normale au lieu d'une direction
d'arrondi.

`gbs/surfaceoffset.h` est retirée, et `gbs/surfaces` ne l'inclut plus.
`SurfaceOffset` ne compilait pas dès
qu'on l'instanciait, n'avait pas de dérivées (ce qui l'excluait des
constructeurs et de `check`, cf. `brep2_pr08_sweep_offset.md`) et doublait
`OffsetSurface`. Son seul apport, une distance variable f(u, v), reviendra sur
`OffsetSurface` si un besoin le demande (point 8).

### 3.5 `gbs/math.ixx` : `make_range(v1, v2, n)`

Chaque valeur est `v1 + (v2 − v1) i / (n − 1)`, arrondie comme dans les
versions 2D et 3D ; v1 et v2 sont rendus exactement. Un point est traité
coordonnée par coordonnée (avant : le pas d'un point était multiplié par un
réel, puis ajouté). Les versions 2D et 3D donnent maintenant les mêmes
valeurs intérieures que la 1D (test `make_range_without_drift`).

## 4. Choix

### 4.1 La normale aux points dégénérés

**La règle.** Les dérivées sont comparées multipliées par les longueurs des
bornes, S_u Δu et S_v Δv (1 pour un côté non borné) : la règle ne dépend ni
des unités ni d'un changement d'échelle des paramètres. Un point est régulier
si |S_u Δu × S_v Δv| > 1e-6 max(|S_u| Δu, |S_v| Δv)², la règle de
`gbs-brep/mesh.h` sur ces dérivées. Sur les dérivées brutes, cette règle
comparait une longueur par radian à une longueur par unité de longueur. Un
cône de 30° haut de 10 mm, en mètres, n'avait pas de normale à son sommet,
alors qu'il en avait une en millimètres ; de même un cône de 1°, un cylindre
de rayon 2e6 en tout point, une bande plane 1e4 × 1e-3 sur [0, 1]².

**La limite.** Ailleurs, `normal_direction` rend **la limite de la normale
depuis l'intérieur du domaine**. Elle prend la normale en trois points, à 1,
2 et 3 pas de 1e-4 du chemin de (u, v) au centre des bornes, et extrapole la
parabole qui passe par eux jusqu'à (u, v). L'erreur est en h³. Si ce chemin
reste sur la ligne dégénérée (du sommet d'un double cône, au milieu du
domaine, vers le centre), elle prend le chemin vers un coin des bornes.

**Les surfaces minces.** Une surface plus de 1e6 fois plus longue que large
(|S_u| Δu < 1e-6 |S_v| Δv autour du point) ne satisfait jamais la règle. Sa
normale est alors prise, au point ou par la même limite, là où S_u et S_v font
un angle de plus de 1e-6 rad. La règle stricte passe d'abord, car elle borne
l'erreur près d'un pôle : la dérivée S_u d'une NURBS y est de l'ordre de
l'erreur d'arrondi des coordonnées, et un angle seul la prendrait pour une
direction.

Mesures :

- au pôle d'une sphère, primitive ou NURBS à rangée de pôles confondus :
  5e-11 ;
- ailleurs : 1e-15 ;
- au sommet d'un cône (de 30°, 1°, 0,5°, 0,1°, et de 10 mm en mètres) : la
  normale de la génératrice du côté du chemin, à 1e-6 ;
- un cylindre de rayon 2e6 : 1e-12 ; la bande 1e4 × 1e-3 : exacte.

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

Le côté du chemin est arbitraire au sommet d'un cône, où la surface n'a pas de
normale : c'est documenté. Si les points de tous les chemins sont eux aussi
dégénérés (surface réduite à un point), ou si (u, v) n'est pas fini, la
fonction lève `std::domain_error`.

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

`discretize` et `discretize_with_normals` évaluent la grille par
`build_batch` avec le seuil `parallel_min_size` : en parallèle au-delà, en
série en dessous. Une exception ne sort jamais d'un algorithme parallèle
(sinon `std::terminate`) : la première, dans l'ordre de la grille, est
relancée après la boucle. Avant, `discretize` terminait le processus sur une
surface qui lève en un nœud, comme une surface décalée où la base n'a pas de
normale. `build_batch` garde deux vecteurs de n entrées (indices et
exceptions) le temps de la boucle. Les tests `discretize_rethrows` et
`discretize_with_normals` le vérifient sur des grilles de 40 × 40 et
100 × 100.

## 5. Changements de comportement

- `parse_file` lève `JsonError` sur une faute de syntaxe, au lieu de laisser
  un document vide ou partiel. Il lit le fichier en binaire : sous Windows,
  les fins de ligne `\r\n` ne sont plus converties, ce que RapidJSON prend
  comme des espaces.
- `parse_file` lève `filesystem_error`, avec la raison du système, si le
  fichier ne peut pas être lu (avant : `std::runtime_error` « File not
  found » ; `filesystem_error` en dérive).
- `write_js_doc` lève `filesystem_error` si le fichier ne peut pas être écrit
  (avant : rien), et laisse un fichier temporaire le temps de l'écriture. Il
  suit les liens symboliques et garde les droits du fichier ; il rompt un lien
  physique.
- `discretize(Surface)` relance l'exception d'une surface qui lève en un nœud,
  au lieu de terminer le processus au-delà de 1000 nœuds.
- `gbs::brep::OffsetSurface` : aux pôles d'une B-spline, le point décalé
  (avant : jusqu'à 0,23 à côté, ou NaN) ; ailleurs d n / |n| au lieu de
  (d / |n|) n, soit un ULP au plus.
- `gbs/surfaceoffset.h` n'existe plus (la classe ne compilait pas).
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
| `tests_files` (nouveau) | 20 000 doubles relus bit pour bit ; une courbe B-spline écrite et relue à l'identique ; la position et la raison d'une faute de syntaxe, un texte vide, un octet nul (document nul après), la copie d'un `JsonError` ; la raison d'un fichier absent ou d'un dossier ; `read_text` sur des octets, un fichier absent, un dossier. |
| `tests_files`, écriture atomique | Remplacement, texte vide, 1 Mio, aucun fichier temporaire laissé. Échec sur un dossier non vide (le dossier intact, le temporaire supprimé). Échec de création dans un dossier absent. Comme sur place : à travers un lien symbolique, les droits 0600 gardés, un fichier en lecture seule refusé. Un nom de 240 octets, un nom en caractères chinois et accentués (écrit puis relu). Quatre écrivains simultanés : le fichier final est entier. |
| `tests_maths.make_range_without_drift` | i / (n − 1) arrondi une fois, en `double` et en `float`, pour un point, égal aux versions 2D et 3D ; 0,3 pour n = 11 ; une suite de 4001 valeurs sur [−1000, 1000]. |
| `tests_bssanalysis` | La normale d'une sphère tournée, primitive et NURBS, à ses pôles compris (1e-9) ; un plan ; une surface réduite à un point et un paramètre NaN (`domain_error`) ; le sommet d'un cône ; sans dépendre des unités : un cône de 10 mm en mètres et en millimètres, des cônes de 1°, 0,5° et 0,1° à leur sommet, un cylindre de rayon 2e6, une bande 1e4 × 1e-3 ; le sommet d'un double cône ; `discretize` et `discretize_with_normals` qui relancent l'exception d'une surface, en série et en parallèle ; `discretize_with_normals` (mêmes points que `discretize`, normales de `normal_direction`). |
| `tests_brep_sweep` | La surface décalée d'une sphère B-spline est une sphère, pôles compris (elle s'en écartait de 0,23). |
| `python/tests/test_BSSurface.py` | `normal_direction` et `discretize_surface_with_normals` depuis Python. |

Les tests qui comparent à des valeurs de référence ont été repassés avec le
nouveau `make_range` : `ctest` et `pytest` passent tous, sans changer une
tolérance ni une référence (§ 6 bis).

### 6 bis. Tests à référence et `make_range`

Suite complète sur macOS (clang 21.1, Release, environnement conda-forge),
avant et après :

| | `master` (56df16a) | Cette branche |
|---|---|---|
| `ctest` | 682 / 682 | 700 / 700 (après la revue, § 9 ; comparaison OCCT comprise) |
| `pytest python/tests` | — | 88 réussis hors notebooks (90 avec), 1 `xfail` déjà marqué |

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
8. **Le retrait de `SurfaceOffset`** (§ 3.4) : elle ne compilait pas et
   doublait `OffsetSurface`. Recommandation : la retirer ; une distance
   variable f(u, v) se mettra sur `OffsetSurface` le jour où elle servira.
9. **La règle des surfaces minces** (§ 4.1) : au-delà d'un rapport de 1e6
   entre les côtés, l'angle seul. Recommandation : la garder, une surface
   longue et fine a une normale.

## 9. Revue de code

Quinze remarques, toutes traitées dans cette PR.

| # | Remarque | Suite |
|---|---|---|
| 1 | Seuil de dégénérescence dépendant des unités et de la paramétrisation : `domain_error` sur un cylindre de rayon 2e6, une bande 1e4 × 1e-3, le sommet d'un cône de 1° ou de 10 mm en mètres. | Dérivées multipliées par les longueurs des bornes ; l'angle seul sur une surface mince (§ 4.1). Tests `normal_direction_whatever_the_units`. |
| 2 | `write_text_atomic` remplace un lien symbolique par un fichier, rompt un lien physique. | Le lien est suivi ; la rupture d'un lien physique est documentée. Test `write_text_atomic_writes_as_in_place`. |
| 3 | Les droits du fichier remplacé sont perdus (0600 → 0644) ; un fichier en lecture seule est remplacé. | Droits copiés sur le temporaire ; un fichier qu'on ne peut écrire est refusé. Même test. |
| 4 | Le repli longe une ligne dégénérée qui passe par le centre (double cône). | Chemin vers un coin quand celui vers le centre échoue. Test `normal_direction_at_the_apex_of_a_double_cone`. |
| 5 | `discretize` termine le processus sur une surface qui lève, ce que fait une surface décalée. | `build_batch` avec le seuil : exception relancée (§ 4.4). Test `discretize_rethrows`. |
| 6 | `path::string()` lève sous Windows pour un nom hors de la page de code. | Nom temporaire en chaînes natives, `utf8_string` pour les messages. Test avec un nom chinois et accentué. |
| 7 | `parse_json` ignore ce qui suit un octet nul. | Un octet nul est une erreur, le document est nul. |
| 8 | Nom temporaire trop long au-delà de 233 octets. | Nom coupé à 200 unités sans couper un caractère. Test d'un nom de 240 octets. |
| 9 | Un test exige que le lecteur par défaut de RapidJSON se trompe. | Assertion retirée. |
| 10 | La copie de `JsonError` peut lever. | `reason()` lit la fin de `what()`. |
| 11 | `parse_file` dit « File not found » pour toute lecture manquée. | `filesystem_error` avec la raison du système. Test `parse_file_says_why_a_file_cannot_be_read`. |
| 12 | Seuil parallèle refait à la main au lieu de `build_batch`. | `build_batch`. |
| 13 | `SurfaceOffset` double `gbs::brep::OffsetSurface`, fausse aux pôles d'une B-spline. | `OffsetSurface` prend `normal_direction` ; `SurfaceOffset` retirée (§ 3.4). Test `offset_surface_at_the_poles_of_a_bspline`. |
| 14 | `discretize_with_normals` construit la grille deux fois et fait deux passes. | Une grille, une boucle. |
| 15 | Boucle hexadécimale du nom temporaire. | `std::format`. |
