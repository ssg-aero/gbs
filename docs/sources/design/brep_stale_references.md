# Références au modèle gardées à travers un ajout

| | |
|---|---|
| PR | branche `fix/brep-stale-references` (sur master 2b2b0a5) |
| Fichiers | `gbs-brep/unify.h` (`join_seams`), `gbs-brep/model.h` (option de test `GBS_BREP_REALLOCATE_ON_ADD`), `tests/tests_brep_unify.cpp`, `News.md`, cette note |
| Statut | À valider (§ 5) |

## 1. Le défaut

`Model` range chaque sorte d'entité dans un `std::vector` (`Model::Table`). Un
ajout (`m.add`) peut le réallouer : toute référence, tout pointeur, tout
itérateur ou `span` pris auparavant sur une entité de la même sorte devient
périmé. Le test `start_at_a_slanted_seam` gardait ainsi une référence à une arête
à travers `fillet_full_round` ; il plantait sous Windows (#258, note
[brep8_full_round_loop_order.md](brep8_full_round_loop_order.md), § 6).

La relecture de ce correctif a trouvé le même défaut dans la bibliothèque :
`join_seams` (`unify.h`) prend des références `d1`, `d2` à deux arêtes, ajoute
l'arête qui les joint, puis lit `d1.u1`. Latent : dans les tests, le tableau des
arêtes avait toujours de la place.

## 2. La correction

`d1` et `d2` sont des copies, prises avant l'ajout.

## 3. L'audit de la classe

Par l'outil plutôt que par la lecture seule : une option de compilation réservée
aux tests, `GBS_BREP_REALLOCATE_ON_ADD`, fait réallouer le tableau à chaque
ajout (`Table::add`, par où passent tous les ajouts, `append` compris). Sous
AddressSanitizer, toute référence gardée à travers un ajout lit alors de la
mémoire libérée et est signalée.

- **L'outil vérifié** : avec le `join_seams` de master, le test ci-dessous
  signale un *heap-use-after-free* (lecture `unify.h:1565`, mémoire libérée par
  l'ajout `unify.h:1547`), dès le cas non compacté.
- **Les suites** : les 30 tests `tests_brep_*` et les 9 de `gbs-io`
  (`tests_io_brep_json`, `tests_io_obj`, `tests_tojson`, les 6 `tests_step_*`),
  compilés avec Apple clang 21, `-fsanitize=address`, l'option active :
  **aucune erreur**, tous les cas passent (les tests d'affichage VTK de
  `tests_brep_mesh` et `tests_tojson`, 4 cas, sont exclus : VTK n'est pas
  compilé pour ce compilateur). L'ASan du clang de l'environnement conda ne
  démarre pas sur ce poste ; celui d'Apple clang est utilisé.

L'audit couvre les chemins que les tests parcourent : un chemin non testé peut
encore garder une référence.

## 4. Tests

`tests_brep_unify.seams_joined_when_the_edges_reallocate` : deux tours de
cylindre empilés, chacun avec sa couture, fusionnés en une face dont la couture
court le long de deux arêtes, jointes par `join_seams`. Une fois tel quel, une
fois le modèle compacté après la couture (`Model::compact` refait les tableaux
juste pleins : l'arête ajoutée par `join_seams` réalloue). Une couture, sur
toute la hauteur, face valide. Sous AddressSanitizer, sans l'option de test,
master signale le *heap-use-after-free* (le cas compacté le garantit) ; corrigé,
aucun.

## 5. Point à valider

1. **Le stockage du modèle reste un `std::vector`** ; l'option de test sert à
   trouver les références périmées. Un `std::deque` fermerait la classe de
   défauts (références stables à l'ajout en fin), sans changer l'API ni les
   temps mesurés (quatre suites lourdes, à 2 % près) ; c'est un choix de
   conception, à part. *Recommandé* dans cette PR : la correction et l'option
   de test seulement.
