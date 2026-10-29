# Lecture STEP — PR 7 : seam, découpe et insertion

| | |
|---|---|
| PR | branche `feat/step-seam` |
| Document de référence | [step_reader.md](step_reader.md), § 5.4 et question 5 (décidée : découpe dans le lecteur) ; PR 4 § 3.2 ; PR 5 § 4 (proposition d'insertion de seam) |
| Fichiers | `gbs-brep/edit.h` (nouveau), `gbs-brep/brep`, `gbs-brep/pcurve.h` (projection globale), `gbs-io/step/reader.h`, `python/gbsBindBrep.cpp`, `tests/tests_brep_edit.cpp` (nouveau), `tests/tests_step_read.cpp`, `tests/in/step/` (4 fichiers, générateur), `python/tests/test_brep.py`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Le document prévoyait ici la **découpe des arêtes qui traversent le seam**
d'une surface périodique. Deux choix des PR précédentes ont déplacé le
problème :

- **PR 4** : une face qui ne fait pas le tour de son axe reçoit une boîte
  angulaire qui évite le seam. Ses arêtes ne traversent donc plus aucun seam.
- **PR 5** : un seam présent dans le fichier (même arête deux fois) est
  accepté.

Reste le cas, fréquent hors OCCT, d'une face qui **fait le tour sans arête de
seam** :

- une bande de cylindre bornée par deux cercles ;
- une calotte sphérique bornée par un cercle ;
- un cône borné par son seul cercle de base.

La NURBS convertie a un seam ; le fichier n'en a pas. Cette PR **insère le
seam**, et **découpe** les arêtes qu'il rencontre (proposition de la note de la
PR 5, retenue).

## 2. API

```cpp
#include <gbs-brep/edit.h>   // inclus par <gbs-brep/brep>

struct SplitEdge { EdgeId first, second; VertexId vertex; };
auto split_edge(Model<T> &m, EdgeId e, T u) -> BuildResult<SplitEdge>;
```

Python : `brep.split_edge(model, edge, u)` → `(first, second, vertex)`.

Le lecteur n'a pas d'API nouvelle : l'insertion est automatique.

## 3. Choix d'architecture et justification

### 3.1 `split_edge` : première opération d'édition du BREP

Couper une arête en u est une opération générique du palier 2 ; elle est donc
dans `gbs-brep/edit.h`, pas dans le lecteur.

- L'arête **garde son identifiant** et devient [u1, u]. Une nouvelle arête
  couvre [u, u2]. Les deux **partagent la courbe 3D** et se rejoignent en un
  nouveau sommet, à la tolérance de l'arête.
- **Tous les wires du modèle** qui utilisent l'arête sont mis à jour : une
  co-arête `Forward` devient (première, seconde), une `Reversed` (seconde,
  première). Les deux moitiés **gardent la pcurve de la co-arête** : elle est
  paramétrée comme l'arête, donc toujours exacte sur chaque moitié. Une face
  déjà construite reste valide sans recalcul (test : disque dont le cercle est
  coupé après coup, `check()` valide).
- Nom et identifiant externe sont recopiés sur la seconde moitié.
- Une arête fermée devient deux arêtes entre son sommet et le nouveau.

### 3.2 Insertion du seam dans le lecteur

Avant de construire une face, le lecteur cherche les boucles qui **s'enroulent
en u** (cumul angulaire de ±2π sur les points échantillonnés). L'insertion a
lieu si la surface est fermée en u seulement et si aucune arête n'apparaît deux
fois dans la face.

1. **Angle du seam u_s** : de préférence l'angle d'un sommet d'une boucle
   enroulée, pour éviter une découpe. Sinon, un angle régulier (64
   candidats). Dans les deux cas, u_s doit être **hors de l'étendue angulaire
   des trous** : le seam ne doit pas traverser un trou.
2. **Point de passage sur chaque boucle enroulée** : un sommet à l'angle u_s,
   sinon le point où une arête atteint u_s. Celui-ci est trouvé par bisection
   sur l'angle, et l'arête y est **coupée** (`split_edge`). Une arête
   partagée garde sa cohérence :
   - le cache des arêtes STEP associe désormais à une arête du fichier la
     **suite** de ses morceaux ;
   - une face lue plus tard reçoit les deux moitiés ;
   - une face déjà construite est mise à jour par `split_edge`.
3. **Extrémité du seam** :
   - **Bande** : le point de passage sur la seconde boucle.
   - **Calotte** : le pôle, côté de la face. Une face est à gauche de sa
     boucle vue depuis sa normale ; elle est donc côté +v si la boucle tourne
     vers +u autour de la normale de la surface (avec `same_sense`). Le pôle
     réutilise le sommet d'une `vertex_loop` s'il y en a une, sinon un sommet
     est créé.
4. **Surface et arête de seam** : la boîte de la face est [u_s, u_s + 2π] en u,
   pôle compris en v. L'arête de seam est **l'isoparamétrique exacte** u = u_s
   de la NURBS de la face (`isoU`), entre les paramètres NURBS de ses deux
   extrémités.
5. **Boucle unique** : la boucle A, le seam vers B, la boucle B, le seam
   retourné. Pour une calotte : A, le seam vers le pôle, puis retour. Les trous
   restent des trous. Le builder de la PR 5 place les deux usages du seam de
   part et d'autre et insère l'arête dégénérée au pôle.

**Si la face échoue**, l'arête de seam et le pôle créé sont effacés. Les
découpes restent : géométrie et topologie inchangées, et les autres faces qui
partagent ces arêtes les utilisent déjà.

**Non pris en charge** (la face est rapportée en échec) :

- surface fermée dans les deux directions (tore) ; levé au lot 8 de la revue
  (§ 7) quand les boucles ne s'enroulent que dans une direction ;
- plus de deux boucles enroulées ;
- boucle enroulée unique sur une surface sans pôle.

### 3.3 Correction de la projection globale (palier 1)

Un point situé à 0,09 rad **avant** le seam de la NURBS était projeté sur le
bord opposé. La recherche bornée (`extrema_surf_pnt`) s'arrête dans un minimum
local au bord u1, à 0,2 mm du point. `project_point_global` relance désormais
Gauss-Newton depuis le nœud le plus proche d'une grille 32 × 32, **seulement
si** le premier résultat n'est pas sur la surface. Le coût reste nul dans le
cas courant. Un test de non-régression le couvre.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| PR 7 : découpe des arêtes qui traversent le seam | découpe **et insertion** du seam ; la traversée elle-même est évitée par la boîte angulaire (PR 4) | le cas réel est l'absence de seam, pas sa traversée |
| découpe propre au lecteur (question 5) | `split_edge` dans `gbs-brep/edit.h`, utilisée par le lecteur | opération générique, réutilisable au palier 2 |
| — | projection globale corrigée dans `gbs-brep/pcurve.h` | défaut révélé par une arête proche du seam |

## 5. Tests

| Test | Couvre |
|---|---|
| `tests_brep_edit.split_edge_in_free_wires` | découpe, courbe partagée, nom et identifiant recopiés, deux wires en sens opposés mis à jour, erreurs sans effet |
| `tests_brep_edit.split_edge_of_a_face_keeps_it_valid` | cercle d'un disque coupé après construction : pcurve partagée, aire inchangée, `check()` valide |
| `tests_brep_edit.projection_near_the_seam` | point à 0,09 rad avant le seam projeté exactement |
| `tests_step_read.seam_insertion` | `band_cylinder.stp` (sommets des cercles à 0 et 1 rad : un cercle coupé), `sphere_cap.stp` (calotte, pôle créé), `cone_noseam.stp` (sommet créé) : solides valides de volume exact ; `band_with_hole.stp` (trou de −0,3 à 0,3 rad) : seam placé hors du trou, les deux cercles coupés, 9 arêtes, `check()` valide |
| `test_brep.py::test_split_edge` | liaison Python |

Toutes les suites BREP et STEP passent, ainsi que les tests Python du BREP.

## 6. Points à valider

1. **`split_edge`** dans le BREP (identifiant conservé, courbe et pcurves
   partagées, tous les wires mis à jour).
2. **Insertion automatique du seam** dans le lecteur : angle d'un sommet de
   préférence, hors des trous, découpe des arêtes au besoin.
3. **Pôle d'une calotte** choisi par le sens de la boucle et `same_sense`.
4. **Découpes conservées** si la face échoue ensuite.
5. **Limites** : tore, plus de deux boucles enroulées.
6. **Correction de `project_point_global`** (grille de repli).

## 7. Revue de code (a posteriori, lot 8 de la revue du noyau BREP)

Chaque commentaire a été reproduit par une sonde ; les tests cités sont nouveaux et
échouent sur la version précédente du lecteur.

| Commentaire | Suite donnée |
|---|---|
| 1. Reporté du lot 3 ([step_pr05_ordered_wire.md](step_pr05_ordered_wire.md) § 7) : la boîte (u, v) d'une face déroulait l'angle au passage d'un pôle par le plus petit saut. Un fuseau de 4 rad entre deux méridiens devenait le fuseau complémentaire, le fuseau de π et les deux moitiés d'une sphère coupée par son seam donnaient un shell non orientable | corrigé : au passage d'un pôle, le saut d'angle suit le côté de la face (`detail::angle_step`) : une face à gauche de ses boucles en (u, v) est sous un pôle en haut de v, où le saut est dans ]−2π, 0], au-dessus d'un pôle en bas, où il est dans [0, 2π[. `parameter_box` reçoit ce sens (`same_sense`, la normale STEP étant S_u × S_v), et le lecteur l'applique aux boucles qu'il compte (`loop_turn`). Ailleurs, le plus petit saut reste la règle : la géométrie y décide seule, comme le veut le § 3.3 de [step_pr06_topology.md](step_pr06_topology.md). Testé (`faces_through_the_poles` : `half_spheres.stp`, `sphere_wedge.stp`, volumes à 1e-3 près) |
| 2. Une sphère coupée par un plan hors de son axe (x = 2,4) : la face sphérique contient les deux pôles et sa boucle ne s'enroule pas. La petite calotte était construite à sa place, et le shell n'était pas orientable | corrigé : une face dont toutes les boucles sont des trous (aire signée en (u, v) contraire à son sens) sur une surface fermée par deux pôles reçoit un seam d'un pôle à l'autre, à un angle hors des trous ; ses boucles deviennent des trous. Testé (`sphere_cut_off_axis.stp`) |
| 3. La même sphère coupée par le plan x = 0 : le cercle passe par les deux pôles, où l'angle de ses points saute de π. L'approximation des pcurves échouait | corrigé : une arête qui passe par un pôle de la surface y est coupée (`split_at_poles`), et les boucles sont rééchantillonnées. Le pôle devient un sommet de la boucle, et son passage est réglé comme au commentaire 1. Testé (`sphere_cut_through_poles.stp`) |
| 4. Limite du § 3.2 : un tore sans arête de seam échouait. C'est le cas d'un congé autour d'une arête de révolution (bande qui s'enroule en u), ou d'un coude de tuyau (boucles qui s'enroulent en v). Ces faces sont aussi courantes que les bandes de cylindre sans seam des exporteurs qui n'en écrivent pas | corrigé : le seam est inséré dans la direction où les boucles s'enroulent, u ou v (`SeamPlan::k`) : ligne iso-v, coupe des arêtes et boîte [v_s, v_s + 2π] pour un coude. Reste non pris en charge : des boucles enroulées dans les deux directions. Testé (`torus_without_seam_edge` : solides dont le volume est exact à 1e-4 près) |
