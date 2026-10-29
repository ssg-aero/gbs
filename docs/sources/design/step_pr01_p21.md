# Lecture STEP — PR 1 : analyseur Part 21 (`gbs-io/step/p21.h`)

| | |
|---|---|
| PR | branche `feat/step-p21` |
| Document de référence | [step_reader.md](step_reader.md), section 2 et question 1 (décidée : analyseur maison en C++23) |
| Fichiers | `gbs-io/step/p21.h` (nouveau), `tests/tests_step_p21.cpp` (nouveau), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Lire la **syntaxe** d'un fichier STEP (ISO 10303-21) : transformer le texte en
une table d'instances et de valeurs, sans rien interpréter du schéma. Une
`CYLINDRICAL_SURFACE` est ici une instance de ce nom avec trois arguments ; la
traduction en géométrie gbs viendra en PR 4.

## 2. Vue d'ensemble

```
parse_p21(std::string_view) / read_p21(path) -> std::expected<P21File, P21Error{line, column, message}>

P21File
  size(), contains(id), find(id) -> optional<InstanceView>, instance(id) (lève si absente)
  instances()            toutes, dans l'ordre du fichier (vue paresseuse)
  instances_of("TYPE")   index par type ; une instance complexe figure sous chacun de ses types partiels
  header(), schemas()    entités d'en-tête, noms de FILE_SCHEMA

InstanceView   id(), line(), is_complex(), type(), has_type("T"), part_count(), part(i | "T"), args(), args("T")
PartView       type(), args()
ListView       size(), [i], at(i), itération
ValueView      kind(), is_null(), as_ref/as_int/as_real/as_string/as_enum/as_bool/as_logical/as_binary/as_list,
               typed_name()/typed_value(), as_measure() (déballe LENGTH_MEASURE(…))
```

## 3. Choix d'architecture et justification

### 3.1 Trois étages, une seule passe sur le texte

1. **Lexer** : découpe le tampon en jetons (`std::string_view` sur le texte, sans
   copie), suit ligne et colonne, saute espaces et commentaires `/* */`.
2. **Analyseur** à descente récursive : `fichier → sections → instances →
   parties → valeurs`. Il remplit des tableaux contigus. Les valeurs imbriquées
   (listes, paramètres typés) sont lues sur une pile explicite, sans récursion
   (§ 7).
3. **Finalisation** : index des `#id` et des types, schémas de l'en-tête.

Les références en avant (`#5` utilisé avant d'être défini) sont normales en
Part 21 : rien n'est résolu pendant la lecture, une référence est un entier ;
`instance(id)` la résout à la demande et lève `P21AccessError` si elle pend.

### 3.2 Stockage en arènes, accès par vues

Comme le `Model` du palier 1, les données vivent dans des `std::vector`
contigus indexés par entier :

| Tableau | Contenu |
|---|---|
| `nodes_` | une valeur par élément : genre + charge utile (entier, réel, référence, ou intervalle) |
| `items_` | enfants des listes, contigus liste par liste |
| `chars_` | chaînes décodées, énumérations, binaires, bout à bout |
| `parts_`, `instances_`, `header_` | parties (type interné, liste d'arguments) et instances |
| `types_`, `by_type_` | noms de types internés (en majuscules), instances par type |
| `dense_` / `sparse_` | `#id` → instance : table dense si les numéros sont compacts (cas courant), table de hachage sinon |

Les vues (`InstanceView`, `ValueView`, `ListView`, `PartView`) sont deux
entiers et un pointeur : copiables à volonté, valides tant que le `P21File`
existe. Aucune allocation par valeur ; une liste imbriquée est construite sur
une pile de travail réutilisée puis recopiée d'un bloc.

**Mesure** : un fichier synthétique de 9 Mo et 200 000 instances (points,
directions, repères, instances complexes) est lu en **30 ms** (Release, clang
21, Apple M). La lecture n'est pas le goulot d'étranglement : l'interprétation
géométrique le sera.

### 3.3 Erreurs : `std::expected` pour la syntaxe, exception pour le contenu

- Une **erreur de syntaxe** rend `std::unexpected(P21Error{ligne, colonne,
  message})` : jeton attendu, chaîne non terminée, réel mal formé ou hors
  limites, énumération ou binaire mal formé, numéro d'instance trop grand,
  imbrication trop profonde, instance définie deux fois, section inconnue,
  fichier tronqué. Un fichier mal formé n'a pas de lecture partielle utile. Une
  chaîne n'est jamais une erreur : un échappement mal formé reste du texte (§ 7).
- Une **erreur de contenu** (valeur d'un autre genre que celui demandé,
  référence pendante, index hors bornes, par `[i]` comme par `at(i)`) lève
  `P21AccessError` depuis les vues.
  L'interprétation (PR 4 et 6) l'attrapera par face ou par entité et la
  rapportera, conformément à l'import partiel décidé.

### 3.4 Couverture de la syntaxe

| Élément | Traitement |
|---|---|
| instances simples et **complexes** `(A()B(…)…)` | parties ordonnées, accès par type, index par type partiel |
| paramètres **typés** `LENGTH_MEASURE(1.E-07)` | `Typed` (nom + valeur), `as_measure()` les déballe |
| entiers, réels (`1.`, `1.E-07`, `.5`, `+3`), énumérations, `$`, `*`, binaires, listes imbriquées | genres dédiés |
| chaînes : `''`, `\\`, `\S\c`, `\X\hh`, `\X2\…\X0\` (paires de substitution UTF-16 recombinées), `\X4\…\X0\`, `\Pc\` | décodées en **UTF-8** à la lecture ; fins de ligne dans une chaîne ignorées |
| commentaires `/* */` partout | sautés |
| plusieurs sections `DATA`, `DATA('nom',(schéma))` de l'édition 3 | lues |
| sections `ANCHOR`, `REFERENCE`, `SIGNATURE` de l'édition 3 | sautées comme texte brut |
| `FILE_SCHEMA` | noms exposés par `schemas()` (AP203, AP214, AP242) |

### 3.5 Moyens C++23 et une exception de portabilité

Moyens utilisés : `std::expected`, `std::string_view`, `std::span`,
`std::ranges` / `std::views` (itération paresseuse des instances),
`std::unreachable`, `std::from_chars` pour les entiers.

**Exception mesurée** : `std::from_chars` pour les **réels** est déclaré
*indisponible* par la libc++ de conda-forge sous macOS (version minimale de
macOS ciblée trop ancienne), et l'Apple clang système le fournit sans définir la
macro `__cpp_lib_to_chars`. La lecture des réels passe donc par `strtod_l` avec
la locale « C » quand la bibliothèque standard est libc++, et par
`std::from_chars` ailleurs (libstdc++ sous Linux, STL de MSVC sous Windows). Les
deux sont exacts et **indépendants de la locale** (un `strtod` simple lirait mal
`1.5` dans une locale française). Le tableau du § 2.2 du document de conception
est corrigé en conséquence.

## 4. Écarts par rapport au document de conception

| Document (§ 2) | Implémentation | Raison |
|---|---|---|
| `std::from_chars` pour les réels | `strtod_l` (locale « C ») sous libc++, `from_chars` ailleurs | indisponible avec la libc++ conda-forge sous macOS (§ 3.5) |
| `P21File` : table des instances, index par type | idem, plus `find`, `contains`, `header()`, `schemas()` | besoins de la PR 4 et des tests |
| Erreur de contenu | `P21AccessError` levée par les vues | les accès sont nombreux et imbriqués ; l'interprétation attrape par entité |

## 5. Tests (`tests_step_p21`, 8 cas, 111 assertions)

| Test | Couvre |
|---|---|
| `header_and_simple_instances` | schéma, en-tête, coordonnées, références, `$`, recherche insensible à la casse, numéro de ligne, instance absente |
| `values` | entiers signés, réels sous toutes leurs formes, `.T.`/`.F.`/`.U.`, énumérations, `$`, `*`, binaire, listes imbriquées et vides, paramètres typés, accès d'un mauvais genre |
| `strings` | `''`, `\\`, `\X2\`, `\X\`, `\S\`, `\X4\`, paire de substitution UTF-16, fin de ligne dans une chaîne, `\PA\`, chaîne vide |
| `complex_instances_and_type_index` | B-spline rationnelle complexe : parties, arguments par type partiel, index par type, ordre du fichier |
| `comments_sections_and_sparse_ids` | commentaires, section `ANCHOR` avec ses jetons `<a>`, deux sections `DATA`, numéro `#1000000000` (index par hachage) |
| `syntax_errors` | `;` manquant (ligne exacte), chaîne non terminée, virgule en trop (ligne et colonne), double définition, échappement invalide, fichier non STEP, fichier tronqué, instance complexe vide, réel invalide |
| `read_file_and_dangling_reference` | lecture depuis un fichier, référence pendante, fichier absent |
| `large_file_performance` | 200 000 instances, 9 Mo : 30 ms localement, borne large de 10 s pour la CI |

## 6. Points à valider

1. **Arènes et vues** plutôt qu'un arbre d'objets alloués un par un.
2. **`std::expected` pour la syntaxe, exception `P21AccessError` pour le contenu.**
3. **Chaînes décodées en UTF-8 dès la lecture.**
4. **`strtod_l` sous libc++** pour les réels, `std::from_chars` ailleurs.
5. **Sections `ANCHOR` / `REFERENCE` / `SIGNATURE` sautées** sans être lues.

## 7. Revue de code (a posteriori, lot 8 de la revue du noyau BREP)

Chaque commentaire a été reproduit par une sonde ; les tests cités sont nouveaux,
et ils échouent sur la version précédente (le premier y fait déborder la pile),
ou n'y compilent pas pour ce qui appelle `detail::read_text`,
`detail::real_underflows` ou `std::ranges::find_if` sur une `ListView`
temporaire. Les mesures portent sur un fichier synthétique de 95 Mo et
1 000 000 d'instances (Release, clang 21, Apple M).

| Commentaire | Suite donnée |
|---|---|
| 1. `parse_list`, `parse_value` et les paramètres typés se rappelaient sans limite. Un fichier de 4 Ko de `(` imbriquées faisait déborder la pile d'un fil de 512 Kio (macOS) : arrêt du processus, pas d'erreur | corrigé : au-delà de `max_nesting` = 256 niveaux de listes ou de paramètres typés, erreur de syntaxe « nesting too deep » à la position de la parenthèse ou du type fautif ; les fichiers réels imbriquent quelques niveaux. Une limite seule ne suffisait pas : 256 appels récursifs prenaient 640 Kio de pile en Debug (2,5 Ko par niveau) et débordaient encore sous ASan. `parse_value` est donc itératif : les listes et paramètres typés ouverts sont sur une pile explicite (`open_`), la pile d'appel reste de 32 Kio quelle que soit l'imbrication, en Debug comme en Release, pour le même temps de lecture. Testé (`deep_nesting_is_an_error` : 100 000 `(` puis 100 000 `T(` dans un fil, erreur ligne 8 colonne 261 ; 255 niveaux acceptés ; aussi en Debug et sous ASan) |
| 2. `ListView::operator[]` ne vérifiait pas l'index, alors que le § 3.3 promettait `P21AccessError`. `args()[2]` d'une instance à un argument lisait `#7`, argument de l'instance suivante ; `args()[1]` de la dernière liste lisait au-delà de `items_`. Les interprétations de `geometry.h` et `reader.h`, qui attrapent `std::exception`, recevaient une donnée fausse au lieu d'un échec | corrigé : `[i]` lève `P21AccessError` (« list index 2 out of range (size 1) ») ; `at(i)` en est un synonyme. Coût mesuré nul : 1,5 ns par valeur parcourue par index, avant comme après. Testé (`list_index_out_of_range_throws`) |
| 3. Une barre oblique inverse isolée ou un échappement mal formé dans une chaîne faisait échouer tout le fichier : `FILE_NAME('C:\Users\me\part.stp',…)` rendait « malformed escape in string » | corrigé : décodage tolérant, une chaîne n'est jamais une erreur. Une séquence `\` inconnue ou mal formée (`\U`, `\X\G1`, `\X2\00E\X0\`, bloc sans `\X0\`, `\S\` final) reste du texte ; les échappements valides sont décodés comme avant, `\S\''` compris. L'assertion de `syntax_errors` qui exigeait l'échec est retirée. Testé (`malformed_escapes_kept_as_text`) |
| 4. Une marque d'ordre d'octets UTF-8 (EF BB BF) avant `ISO-10303-21;` était refusée (« unexpected character », ligne 1 colonne 1) | corrigé : sautée par le lexer. Testé (`utf8_byte_order_mark_skipped` : lignes et colonnes des erreurs inchangées) |
| 5. Les réels différaient selon la plateforme. La branche `std::from_chars` (libstdc++, MSVC) refusait `+1.5` et `+1.E2`, valides en Part 21, et `1.E-400` ; la branche `strtod_l` (libc++) ne lisait pas `errno` (`1.E400` lu comme l'infini) et refusait les jetons de 128 caractères et plus (tampon fixe) | corrigé : un comportement partout. Le `+` initial est sauté, un dépassement est une erreur (« real out of range »), un sous-dépassement donne 0 (ou le dénormalisé), sans limite de longueur (tampon sur la pile jusqu'à 63 caractères, `std::string` au-delà). Avec `from_chars`, qui laisse la valeur indéfinie hors bornes, le sens du dépassement se lit sur les chiffres (`detail::real_underflows`). La branche `strtod_l` est gardée : `from_chars` flottant reste indisponible avec la libc++ de conda-forge. Testé (`reals_read_alike_on_every_platform`, sur les deux branches : la branche `from_chars` forcée avec l'Apple clang) |
| 6. Le résultat de `std::from_chars` était ignoré pour les numéros d'instance : `#18446744073709551616` était lu `#0`, et deux numéros trop grands donnaient « instance #0 defined twice » | corrigé : « instance number too large », à la définition comme à la référence. Testé (`instance_number_too_large` : 2^64 refusé, 2^64 − 1 accepté) |
| 7. Énumérations et binaires n'étaient pas validés, et leurs fins de ligne pas comptées : `.UNSPEC` + saut de ligne + `IFIED.` était lu comme une énumération contenant un saut de ligne, `"hello world"` comme un binaire, et après un binaire sur trois lignes les erreurs étaient signalées deux lignes trop haut | corrigé selon la grammaire de la norme : ENUMERATION = `.` UPPER { UPPER \| DIGIT } `.`, UPPER = A–Z et `_` (les minuscules restent acceptées et mises en majuscules, comme les mots-clés) ; BINARY = `"` (0 à 3) { 0–9 A–F } `"`, dont les sauts de ligne sont sautés et comptés. Testé (`enumerations_and_binaries_validated`) |
| 8. Le décodeur produisait de l'UTF-8 invalide, alors que l'API promet de l'UTF-8 : `\X2\DC00\X0\` (demi-paire isolée) donnait ED B0 80, une paire coupée en deux blocs `\X2\` adjacents ED B8 80, `\X4\7FFFFFFF\X0\` FF BF BF BF, et l'octet brut E9 d'un fichier ISO 8859-1 restait tel quel | corrigé : la demi-paire haute en attente passe d'un bloc `\X2\` au suivant (même par-dessus une fin de ligne) ; une demi-paire isolée ou un point de code au-delà de U+10FFFF devient U+FFFD ; une suite UTF-8 brute valide est gardée, un octet brut qui n'en forme pas est lu en ISO 8859-1, le codage par défaut de Part 21. Le résultat est toujours de l'UTF-8 valide. Testé (`strings_always_valid_utf8`) |
| 9. Le lexer utilisait `<cctype>` (`isalpha`, `isdigit`, `toupper`…), qui dépend de la locale C globale : dans la locale fr_FR.UTF-8, le mot-clé `CAF` + octet C9 était accepté | corrigé : classes de caractères ASCII (`detail::is_alpha`, `is_digit`, `is_space`, `to_upper`), `<cctype>` n'est plus inclus. Testé (`lexer_ignores_the_global_locale` : C, fr_FR, tr_TR, de_DE quand elles existent) |
| 10. L'itérateur de `ListView` gardait un pointeur sur la `ListView` : celui d'un temporaire (`inst.args()[1].as_list().begin()`) pendait. Même défaut pour `instances()`, vue `transform` dont l'itérateur pointe sur la vue | corrigé : un seul itérateur à accès direct, `detail::IndexIterator`, qui porte le fichier et l'index ; `ListView` et la nouvelle `InstanceRange` (rendue par `instances()`) sont des vues empruntées (`std::ranges::enable_borrowed_range`), à accès direct. Testé (`iterators_outlive_their_view` : itérateurs de temporaires, `std::ranges::find_if` sur un temporaire, `back()`, `instances()[1]`) |
| 11. `read_p21` lisait le fichier par `istreambuf_iterator` : 185 à 212 ms pour 95 Mo | corrigé : `detail::read_text` (taille lue en fin de fichier, une lecture dans une chaîne déjà dimensionnée, lecture par flux si le fichier ne se positionne pas, `nullopt` pour un dossier ou un fichier absent) : 5 à 9 ms. `read_step_file` (`reader.h`) s'en sert aussi. Testé (`read_text_whole_file`) |
| 12. Mémoire : un `Node` de 24 octets, et des tableaux qui croissaient sans réserve. Capacités cumulées 532 Mo, 5,6 fois le fichier (`nodes_` : 16,8 millions de places pour 11 millions de valeurs, 403 Mo) ; empreinte mesurée 802 Mo (8,4 fois) ; lecture en 353 à 363 ms | corrigé : `Node` de 12 octets, la charge utile tenant dans deux mots de 32 bits (`get` / `set` par `std::bit_cast`) ; réserves calculées sur la taille du texte (7 à 9 octets de texte par valeur, 35 à 45 par entité partielle, mesurés sur les fichiers réels du dépôt), index par type compté avant d'être rempli, grand excédent rendu à la fin seulement (la place réservée et non touchée ne coûte pas de mémoire physique). Capacités 279 Mo (2,9 fois), empreinte 274 Mo (2,9 fois), lecture en 249 à 255 ms. Testé (`value_nodes_are_compact`) |
| 13. `type_id` construisait une copie en majuscules du nom à chaque `has_type`, `part` ou `instances_of` : 32 ns par appel, et une allocation pour un nom de plus de 22 caractères | corrigé : recherche hétérogène (`detail::TypeHash` transparent et `std::equal_to<>`, pas de `std::flat_map`, absent de MSVC et de l'Apple clang de la CI), essai exact d'abord, copie en majuscules seulement pour un nom qui contient des minuscules ; `intern` de même pendant la lecture. 5 ns par appel. Testé (`type_lookup_does_not_allocate` : aucune allocation pour un nom de 25 caractères, recherche toujours insensible à la casse) |
