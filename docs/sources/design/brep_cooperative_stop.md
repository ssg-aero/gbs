# Arrêt coopératif des opérations longues (`StopRequested`, `BuildErrc::Cancelled`)

| | |
|---|---|
| PR | branche `feat/brep-cooperative-stop` |
| Fichiers | `gbs-brep/model.h` (`StopRequested`, `BuildErrc::Cancelled`, `detail::stop_requested`, `detail::cancelled`, `detail::stop_point`, `Model::releasePartKeys`), `gbs-brep/boolean.h`, `gbs-brep/partial.h` (`split_partial_coincidences`), `gbs-brep/section.h`, `gbs-brep/mesh.h`, `gbs-brep/sweep.h`, `gbs-brep/pipe.h`, `gbs-brep/fillet.h` (champ `stop` de `FilletOptions`), `gbs-brep/blend.h` (quatre points d'arrêt), `gbs-io/step/reader.h`, `python/gbsBindBrep.cpp` (`BuildErrc.cancelled`), `tests/tests_brep_stop.cpp` (nouveau), `docs/index.rst`, `News.md` |
| Exigences | besoin G7 d'une application de CAO (CadSmith) : arrêter un booléen, un maillage, une lecture STEP en moins d'une seconde |
| Statut | Relu ; décisions prises (§ 10) |

## 1. Objectif

Aucune opération longue de gbs ne pouvait être interrompue : une application de CAO (CadSmith) qui la
lance en tâche de fond sur une copie de son modèle ne peut que jeter le résultat à son arrivée ; le
thread reste occupé jusqu'au bout, et quitter attend la fin du calcul. Il faut un arrêt demandé
d'un autre thread, que l'opération voit assez souvent pour s'arrêter en **moins d'une seconde**, en
rendant une erreur distincte, le modèle **inchangé**, une autre opération réussissant ensuite.

## 2. API

```cpp
namespace gbs::brep
{
    enum class BuildErrc : std::uint8_t { /* ... */ UnsupportedVersion, Cancelled };

    /// Appelé par l'opération entre deux étapes de son travail ; vrai : elle s'arrête là et
    /// échoue avec BuildErrc::Cancelled, le modèle tel qu'il était. Vide par défaut.
    using StopRequested = std::function<bool()>;

    template <std::floating_point T> struct BooleanOptions { /* ... */ StopRequested stop{}; };
    template <std::floating_point T> struct SectionOptions { /* ... */ StopRequested stop{}; }; // section, split
    template <std::floating_point T> struct MeshOptions    { /* ... */ StopRequested stop{}; };
    template <std::floating_point T> struct SweepOptions   { /* ... */ StopRequested stop{}; }; // make_prism, make_revolution
    template <std::floating_point T> struct PipeOptions    { /* ... */ StopRequested stop{}; };
    template <std::floating_point T> struct FilletOptions  { /* ... */ StopRequested stop{}; }; // fillet, chamfer, chamfer_angle, fillet_chordal

    template <std::floating_point T> class Model
    {
        // les clés de pièce données depuis `next` rendues (retour arrière de read_step)
        void releasePartKeys(std::uint64_t next);
    };

    namespace detail
    {
        bool stop_requested(const StopRequested &) noexcept; // stop && stop() ; une exception du rappel : vrai
        BuildError cancelled(std::string_view op);         // {Cancelled, "<op>: stopped on request", {}}
        void stop_point(const StopRequested &, std::string_view op); // lève BuildFailure(cancelled(op))
    }
}
namespace gbs::step
{
    struct StepReadOptions { /* ... */ brep::StopRequested stop{}; };
}
```

Usage depuis une application :

```cpp
std::atomic<bool> stop_asked{false};               // ou un std::stop_source
gbs::brep::BooleanOptions<double> o;
o.stop = [&] { return stop_asked.load(std::memory_order_relaxed); };
auto r = gbs::brep::cut(model, a, b, o);           // dans le thread de calcul
// ailleurs : stop_asked = true;  ->  r.error().code == BuildErrc::Cancelled
```

Le champ est **le dernier** de chaque structure d'options et vide par défaut : les appelants
actuels (initialisation par agrégat positionnelle ou désignée, liaisons Python) ne changent pas.

## 3. Choix d'architecture

### 3.1 Un rappel `std::function<bool()>` plutôt que `std::stop_token`

`std::stop_token` aurait été l'idiome C++20, mais libc++ ne le fournit hors du mode expérimental
(`-fexperimental-library`) qu'à partir de sa version 20 : le clang d'Xcode 16, celui du job macOS
de la CI et des postes de nombreux utilisateurs, s'appuie sur une libc++ plus ancienne. Un
rappel couvre tous les cas : un `std::stop_token` s'y emballe en une ligne
(`[t = source.get_token()] { return t.stop_requested(); }`), un drapeau atomique aussi, et un
appelant peut y mettre une échéance. Le coût d'un `std::function` vide testé à chaque point
d'arrêt est négligeable (§ 6).

Contrat du rappel (documenté sur `StopRequested`) : appelé du thread de l'opération, souvent
(jusqu'à quelques milliers de fois par seconde), il doit être bon marché et sûr pendant qu'un autre
thread demande l'arrêt ; une fois vrai, il le reste. Ce dernier point est utile : un congé qui
attrape une erreur pour réessayer autrement (`through_pockets`) retombe aussitôt sur l'arrêt.

Il ne doit pas lever. S'il lève quand même, l'exception vaut une demande d'arrêt :
`detail::stop_requested` l'attrape et rend vrai, l'opération échoue avec `Cancelled`, le modèle tel
qu'il était, et l'exception est perdue. C'est le seul choix qui tienne la garantie dans toutes les
opérations par le même chemin : le maillage, les congés et la lecture STEP attrapaient déjà toute
exception (en échec `NonGeneric` ou `InvalidFile`), le booléen, la section et les balayages ne la
rattrapaient pas et laissaient le modèle modifié. Ces trois-là reviennent aussi maintenant en
arrière avant de laisser passer toute autre exception (`catch (...)`, retour arrière, `throw`).

### 3.2 Un code `Cancelled`, pas d'exception publique

Chaque opération rend `BuildErrc::Cancelled` par son `BuildResult` habituel, sans entités
(`shapes` vide), avec le message `"<opération>: stopped on request"`. À l'intérieur, deux formes :

- là où l'opération rend directement son résultat (boucles de `boolean`, `section`, `mesh`, des
  balayages), un test `detail::stop_requested` suivi d'un `return std::unexpected(detail::cancelled(...))` ;
- là où le point d'arrêt est loin de l'entrée (classification d'une face solide par solide,
  raffinement d'une face, étapes d'un congé), `detail::stop_point` lève `BuildFailure` : chaque
  opération l'attrape déjà pour revenir en arrière (`boolean`, `sweep_detail::sweep`,
  `blend_detail::guarded`, `MeshBuilder::face`) ; la lecture STEP lève son `StrictStop`, qui
  traverse déjà tous ses `catch` et revient en arrière en mode strict comme en mode souple.

### 3.3 Le retour arrière existant

Les opérations travaillent sur des copies et effacent, en cas d'échec, tout ce qu'elles ont ajouté
(`detail::arena_size` / `detail::rollback`) : l'arrêt en profite sans rien ajouter. Deux manques
sont apparus avec les tests (modèle comparé par `brep::equal` après compactage) :

- **`section`** ne revenait pas en arrière : elle n'échouait jamais après avoir créé des arêtes.
  L'arrêt efface ce qu'elle a ajouté ;
- **`read_step`** laissait avancé le compteur des clés de pièce (`Model::newPart`, assemblages lus)
  après un échec en mode strict, et donc après un arrêt. `Model::releasePartKeys(next)` rend les clés
  depuis `next` (refus si une instance vivante en use encore) ; le retour arrière de `read_step`
  l'appelle avec `nextPart()` d'avant la lecture. Un échec strict en profite aussi.

## 4. Points d'arrêt

| Opération | Points d'arrêt |
|---|---|
| `boolean` (`fuse`, `cut`, `common`) | étape 0, les coïncidences partielles (`split_partial_coincidences`) : avant chaque face B-spline porteuse, chaque paire de faces analysée, chaque rangée de carreaux de Bézier d'une paire ; la section (ci-dessous) ; entre deux faces découpées ; entre deux faces classées et, pour chacune, entre deux solides de l'autre opérande ; avant l'unification ; avant chaque coque assemblée |
| `section` (deux formes, ou une forme et une surface), `split` | entre deux paires de faces dont les boîtes se touchent (la découpe qui suit la section, faite en place, ne s'arrête pas) |
| `mesh` | entre deux faces ; dans une face : à chaque passe du raffinement, toutes les 4096 triangles examinés, tous les 1024 points insérés ; toutes les 4096 triangles orientés |
| `MeshBuilder::face` | les points d'une face ci-dessus ; le constructeur, réutilisable après un arrêt |
| `read_step` | avant chaque face construite, chaque instance placée (copie), et avant le contrôle final (`check`) |
| `make_prism`, `make_revolution` | entre deux arêtes du profil balayées ; avant la couture |
| `make_pipe` | entre deux morceaux du chemin balayés par chaque arête du profil ; ceux des balayages (message « make_pipe: stopped on request ») |
| `fillet`, `chamfer`, `chamfer_angle`, `fillet_chordal` | entre deux arêtes (géométrie de chaque congé), avant chaque prolongement, avant les bouts, avant les contacts |

## 5. Garanties

- Arrêt demandé **avant** : l'opération rend `Cancelled` au premier point, le modèle inchangé.
- Arrêt **pendant** : `Cancelled`, le modèle inchangé : mêmes entités vivantes, égales au bit près
  (`brep::equal` après `compact()` des deux côtés : seules les entités ajoutées puis effacées, en
  fin de tables, diffèrent), même unité, même compteur de pièces.
- Arrêt **après** la fin, ou jamais demandé : le résultat est celui d'une opération sans arrêt.
- L'opération suivante réussit (rien n'est laissé à moitié, aucun état global).
- Un rappel qui lève : comme un arrêt demandé (§ 3.1).

Non garanti : un arrêt ne coupe pas un calcul indivisible en cours (une intersection de deux
surfaces, la classification d'une face contre un solide, le volume d'une coque à l'assemblage,
l'approximation d'une surface de congé) ;
le délai est celui du plus long de ces calculs (§ 7).

## 6. Coût sans arrêt

Sans arrêt (champ vide), chaque point d'arrêt teste un `std::function` vide : un booléen de 2 s
du cas lourd (§ 7) en passe 7 034, un maillage d'une seconde 2 187, soit quelques microsecondes
en tout. Avec un rappel qui lit un drapeau atomique, l'écart reste dans le bruit de la mesure :
le `cut` du cas lourd, trois fois de suite chacun, 1,839, 1,848 et 1,848 s sans rappel ; 1,831,
1,856 et 1,851 s avec.

## 7. Délais mesurés

Cas lourds construits par programme (`tests_brep_stop`, cas `delays`) : une plaque de 6 × 6
cellules, une goupille (cylindre) traversant chaque cellule et une sphère à moitié enfoncée entre
elles (62 solides) ; le `cut` de la plaque par leur composé, le maillage du résultat (déflexion
1e-4), la lecture du fichier STEP des opérandes et de cinq copies du résultat (7,7 Mo). L'arrêt est
demandé d'un autre thread 20 ms après le premier point d'arrêt ; le plus long écart entre deux
points d'arrêt d'une exécution complète (début et fin compris) est le délai d'arrêt au pire.
Mac M (Release, clang 21), machine chargée par d'autres compilations :

| Opération | Durée complète | Points d'arrêt | Plus long écart | Arrêt mesuré |
|---|---|---|---|---|
| booléen (`cut`) | 1,9 à 2,0 s | 7 034 | 15 à 17 ms | 0,2 à 0,3 ms |
| maillage | 0,8 à 1,0 s | 2 187 | 21 à 27 ms | 0,3 à 2,6 ms |
| lecture STEP (7,7 Mo) | 0,19 à 0,22 s | 475 | 83 à 97 ms (l'analyse du texte, avant la première face) | 0,2 à 6 ms |
| booléen de faces B-spline (cas `delays_bspline`) | 0,94 s | 129 | 258 ms (le volume de la coque, à l'assemblage) | 9 ms |

Le cas `delays_bspline` : une pale lissée par 24 sections tordues (une face B-spline de nombreux
carreaux), coupée par sa copie tordue et déplacée d'un peu : les faces restent proches sur toute
leur étendue, et l'étape 0 du booléen (les coïncidences partielles) examine chaque carreau de
chacune contre l'autre : 0,53 s à elle seule. Sans ses points d'arrêt, la relecture a mesuré 1,45 s
d'un seul tenant sur une pièce réelle de 22 faces B-spline (un arrêt vu après 1,4 s) ; avec eux, le
plus long écart de ce cas tombe dans l'assemblage.

Le test exige, pour l'écart au pire comme pour l'arrêt mesuré, moins d'une seconde sur un poste,
et moins de 5 s quand la variable `CI` est posée (machines partagées de l'intégration continue,
plus lentes et chargées) ; les délais mesurés sont toujours affichés.

Face isolée très fine (une sphère de rayon 1, déflexion 1e-6) :

| `max_points` | Durée | Plus long écart |
|---|---|---|
| 200 000 (défaut) | 0,76 s | 36 ms |
| 1 000 000 | 5,6 s | 0,29 s |
| 3 000 000 | 19 s | 1,1 à 1,4 s |

Le dernier écart est dans `ConstrainedDelaunay::domain()` et `remove_outside()` (`gbs-mesh/cdt.h`),
d'un seul tenant, sur des millions de triangles : hors de portée sans toucher `gbs-mesh` ; la
borne par défaut de 200 000 points par face le garde loin.

Avant les points d'arrêt dans la classification (solide par solide), une plaque de 12 × 12 cellules
(265 solides) laissait un écart de 6,4 s : la classification d'une seule face (la face supérieure
percée de 265 trous) contre les 265 solides de l'autre opérande. Voir § 9.

## 8. Tests (`tests/tests_brep_stop.cpp`)

- `no_stop_asked` : un rappel jamais vrai donne le même modèle que sans rappel ; `to_string(Cancelled)`.
- `boolean` (fuse, cut, common ; un cylindre, une sphère), `section_and_split` (deux formes, une
  surface, `split`), `mesh` (et `MeshBuilder::face` arrêté puis repris), `read_step` (dans un modèle
  non vide, mode souple et strict, avec un assemblage : compteur de pièces), `sweeps` (prisme,
  demi-tour et tour complet, tuyau le long d'un coude), `blends` (trois congés en coin, chanfreins) :
  pour chacun, `expect_stops` demande l'arrêt avant, à chaque point d'arrêt (tous jusqu'à 100, un
  échantillon de 50 au-delà), et après la fin, et donne un rappel qui lève à mi-chemin ; à chaque
  arrêt `Cancelled` sans entités et le modèle égal à celui d'avant (`brep::equal` après
  compactage) ; puis l'opération sans arrêt réussit.
- `delays`, `delays_bspline` : § 7.

## 9. Constats hors du périmètre

- **Classification lente d'une face à nombreux trous** contre un opérande de nombreux solides :
  `boolean_detail::locate` classe la face contre chaque solide (`classify(m, f, s)`), qui recalcule
  à chaque fois un point intérieur de la face et ses rayons : 25 ms par solide pour une face de 266
  boucles, 6,5 s pour 265 solides. Un point intérieur calculé une fois par face, ou une boîte
  englobante par solide, le réduirait beaucoup. À traiter à part, comme le coût de l'étape 0 sur
  des faces B-spline proches (§ 7).
- **Suite (« G7b ») : opérations non couvertes.** Une application les lance aussi en tâche de
  fond : `check`, les intégrales de volume et d'aire, l'écriture STEP et l'écriture et la lecture
  `brep_json` ; courtes sur une pièce (`check` : environ 0,1 s), elles ne sont pas bornées sur un
  grand assemblage. Puis, si le besoin vient : `offset_solid`, `make_thick` (`OffsetOptions`),
  `sew`, `fillet_faces`, `fillet_full_round`, `unify_faces` seul, les fonctions `*_geometry` des
  congés, la lecture IGES. Le champ `stop` de `FilletOptions` n'est lu que par `fillet`,
  `chamfer`, `chamfer_angle` et `fillet_chordal`.
- L'analyse du texte STEP (`parse_p21`, environ 80 Mo/s) n'a pas de point d'arrêt : moins d'une
  seconde jusqu'à 80 Mo environ.
- Les liaisons Python n'exposent pas `stop` (un rappel Python appelé sans le GIL du thread de calcul
  demanderait un soin à part) ; elles connaissent `BuildErrc.cancelled`.

## 10. Points validés

1. Un rappel `std::function<bool()>` (`StopRequested`) plutôt que `std::stop_token` (§ 3.1) :
   retenu (choix technique, libc++ avant 20).
2. Le nom `stop` du champ, le dernier de chaque structure d'options ; le code `Cancelled` à la fin
   de `BuildErrc` (ses valeurs numériques inchangées) : retenu. Un `switch` exhaustif sans
   `default` d'un appelant devient un avertissement, comme pour tout nouveau code.
3. Les opérations couvertes et celles laissées pour plus tard : `check`, les intégrales,
   l'écriture STEP et `brep_json` attendent une suite, « G7b » (§ 9).
4. `Model::releasePartKeys` et son usage dans le retour arrière de `read_step`, qui corrige aussi
   un échec strict (§ 3.3) : retenu.
5. Le retour arrière ajouté à `section`, sur arrêt et sur toute exception : retenu.
6. Le seuil du test `delays` : une seconde sur un poste, 5 s en intégration continue, les délais
   toujours affichés (§ 7).
7. Après la relecture : des points d'arrêt dans l'étape 0 du booléen (§ 4, § 7) ; un rappel qui
   lève vaut un arrêt (§ 3.1) ; le message d'un tuyau arrêté ; tous les points d'arrêt des petits
   cas essayés (§ 8).
