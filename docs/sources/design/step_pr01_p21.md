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
   parties → valeurs`. Il remplit des tableaux contigus.
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
  message})` : jeton attendu, chaîne non terminée, échappement invalide, réel
  mal formé, instance définie deux fois, section inconnue, fichier tronqué. Un
  fichier mal formé n'a pas de lecture partielle utile.
- Une **erreur de contenu** (valeur d'un autre genre que celui demandé,
  référence pendante, index hors bornes) lève `P21AccessError` depuis les vues.
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
