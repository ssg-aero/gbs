# Écriture STEP — PR 1 : écrivain Part 21 (`gbs-io/step/p21_writer.h`)

| | |
|---|---|
| PR | branche `feat/step-p21-writer` |
| Document de référence | [step_writer.md](step_writer.md), § 2 |
| Fichiers | `gbs-io/step/p21_writer.h` (nouveau), `tests/tests_step_p21_writer.cpp` (nouveau), `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Écrire la **syntaxe** d'un fichier STEP : en-tête, puis une table d'instances
simples ou complexes dont les paramètres sont des valeurs Part 21. C'est le
pendant de l'analyseur de la PR 1 de la lecture. Comme lui, il ne connaît
aucun schéma : les entités BREP viendront dans les PR 2 et 3.

## 2. API

```cpp
P21Value::ref(id)  integer(i)  real(x)  string(utf8)  enumeration("NAME")  boolean(b)  logical(opt)
P21Value::omitted()  derived()  list({…})  refs({ids})  reals({xs})  integers({xs})  typed("LENGTH_MEASURE", v)

P21Header { description, implementation_level, name, time_stamp, author, organization,
            preprocessor_version, originating_system, authorization, schemas }
schema::ap214, schema::ap242, schema::ap203

P21Writer w{header};
auto id = w.add("CARTESIAN_POINT", {P21Value::string(""), P21Value::reals({0., 0., 0.})});
auto c  = w.add_complex({{"B_SPLINE_CURVE", {…}}, {"RATIONAL_B_SPLINE_CURVE", {…}}, …});
w.str();  w.write(path);
```

## 3. Choix d'architecture et justification

### 3.1 Valeurs formatées à la création

Une `P21Value` est son texte Part 21, produit dès sa création. Il n'y a ni arbre
de valeurs ni seconde passe, et les erreurs apparaissent au plus tôt, à
l'endroit fautif :

- réel non fini ;
- chaîne UTF-8 invalide ;
- nom d'énumération ou de type invalide ;
- référence `#0`.

Elles lèvent `std::invalid_argument`. Les instances sont numérotées dans
l'ordre d'ajout. L'écrivain d'entités ajoutera les feuilles d'abord, si bien
qu'aucune référence vers l'avant n'est nécessaire.

### 3.2 Réels : forme la plus courte qui se relise à l'identique

- Avec libstdc++ et la STL de MSVC, `std::to_chars` donne la forme la plus
  courte exacte.
- Sous libc++, `to_chars` sur les flottants est indisponible, comme
  `from_chars` (PR 1 de la lecture). On prend alors la plus courte de
  `%.15g`, `%.16g`, `%.17g` qui se relit à l'identique par `strtod_l`, avec
  `snprintf_l` dans la locale « C ».
- Dans les deux cas, la sortie est **indépendante de la locale**.
- La forme est celle de Part 21 : un point décimal toujours présent (`1.`,
  `-0.`), exposant en `E` sans zéros de tête (`1.E-7`, `1.E300`).
- Le test relit **bit à bit** plus de 2 000 réels : valeurs limites,
  dénormalisés, `-0.` et 2 000 tirages sur toute l'étendue des exposants.

### 3.3 Chaînes : l'inverse exact du décodage

| Caractère | Écrit |
|---|---|
| ASCII imprimable | tel quel |
| apostrophe, antislash | `''`, `\\` |
| hors ASCII ou de contrôle, dans le plan de base | suite `\X2\hhhh…\X0\` |
| au-delà du plan de base (émoticônes…) | suite `\X4\hhhhhhhh…\X0\` |

Le fichier produit est donc de l'ASCII pur. Les suites consécutives de même
largeur sont regroupées. L'analyseur les redécode en UTF-8 : le test fait
l'aller-retour sur des accents, des idéogrammes, une émoticône et un retour à
la ligne.

### 3.4 Instances complexes dans l'ordre alphabétique

Part 21 impose l'**ordre alphabétique** des entités partielles d'une instance
complexe. `add_complex` trie : on peut donner les parties dans n'importe quel
ordre. Cet ordre est l'ordre ASCII (`BOUNDED_CURVE` avant `B_SPLINE_CURVE`),
celui qu'écrit OCCT.

### 3.5 Lignes et en-tête

- Les instances longues sont coupées **après une virgule, hors des chaînes**,
  vers 80 colonnes. Part 21 n'impose pas de longueur maximale, mais les
  fichiers restent lisibles et certains lecteurs anciens en ont une.
- L'en-tête suit le standard : `FILE_DESCRIPTION`, `FILE_NAME` et
  `FILE_SCHEMA`.
- L'horodatage vaut l'heure UTC courante par défaut, au format ISO 8601.
  `gmtime_r` (`gmtime_s` sous Windows) est utilisé plutôt que `<chrono>`
  formaté, absent de certaines libc++.

## 4. Écarts par rapport au document de conception

| Document (§ 2.2) | Implémentation | Raison |
|---|---|---|
| valeurs `Ref, Int, Real…` en variante | `P21Value` formatée à la création | plus simple, sans seconde passe ; erreurs au plus tôt |
| — | `refs`, `reals`, `integers` | listes les plus fréquentes des entités géométriques |
| binaires | non écrits | aucune entité BREP n'en utilise |

## 5. Tests (`tests_step_p21_writer`, 5 cas, 2 103 assertions)

Chaque texte écrit est relu par `parse_p21`.

| Test | Couvre |
|---|---|
| `values_read_back` | tous les genres de valeurs, listes imbriquées et vides, paramètre typé, références |
| `reals` | aller-retour bit à bit de plus de 2 000 réels, formes `1.`, `-0.`, `1.E-7`, `1.E300`, refus de l'infini et de NaN |
| `strings` | apostrophe, antislash, accents, idéogrammes, émoticône, retour à la ligne ; fichier en ASCII pur ; UTF-8 invalide refusé |
| `complex_instances_header_and_lines` | ordre alphabétique des parties, en-tête complet (auteur accentué, schéma AP242), lignes coupées hors chaînes, horodatage par défaut |
| `errors_and_file` | noms invalides, instance complexe vide, `#0`, écriture et relecture d'un fichier, dossier inexistant |

## 6. Points à valider

1. **Valeurs formatées à la création**, erreurs par `std::invalid_argument`.
2. **Réels les plus courts exacts**, indépendants de la locale, sous libc++
   comme ailleurs.
3. **Chaînes en ASCII pur** (`\X2\`, `\X4\`).
4. **Tri alphabétique** des parties des instances complexes.
5. **Lignes coupées vers 80 colonnes** hors chaînes.
