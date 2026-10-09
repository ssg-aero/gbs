# BREP palier 1 — PR 7 : sewing simple (`gbs-brep/sew.h`)

| | |
|---|---|
| PR | branche `feat/brep-sew` |
| Document de référence | [brep_core.md](brep_core.md), sections 4.3 et 6.3, question 12 ; notes précédentes [PR 1](brep_pr01_model.md) à [PR 6](brep_pr06_check.md) |
| Fichiers | `gbs-brep/sew.h` (424 l., nouveau), `gbs-brep/brep`, `tests/tests_brep_sew.cpp` (244 l., nouveau), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Assembler des faces construites indépendamment (chacune avec ses propres
arêtes et sommets) en shells : apparier les arêtes libres qui coïncident à
tolérance près, les fusionner, regrouper les faces par composante connexe et
les orienter de façon cohérente. C'est le remplaçant natif de
`occt_utils::to_shell` (`BRepBuilderAPI_Sewing`) de `gbs-occt`, **sans
découpe d'arêtes** (pas de jonction en T), conformément à la question 12 du
document.

## 2. Vue d'ensemble

```
sew(m, faces, SewOptions{tol = 1e-6, n_samples = 10, n_pcurve_max = 513}) -> BuildResult<SewReport>

SewReport { shells, merged (survivant, absorbé), free_edges, rejected, ambiguous, orientable }

1. candidats     arêtes libres (1 usage dans les faces données, non dégénérées), boîte gonflée de tol, tri sur x
2. paires        boîtes qui se coupent (balayage) → extrémités confondues à tol (même sens ou sens opposé)
                 → distance échantillonnée dans les DEUX sens ≤ tol (sinon « rejected »)
3. fusion        paires triées par distance croissante, chaque arête au plus une fois (sinon « ambiguous »)
                 survivant = plus petit id : garde sa courbe ; la co-arête de l'absorbé est redirigée,
                 son sens composé avec le sens relatif, sa pcurve reconstruite dans le paramètre du survivant
                 sommets appariés fusionnés (union-find), tolérances remontées, absorbés effacés
4. orientation   parcours en largeur des faces adjacentes ; un shell par composante ; drapeau closed calculé
```

## 3. Choix d'architecture et justification

### 3.1 Appariement : extrémités puis distance dans les deux sens

Le test d'extrémités est un filtre rapide et décide du sens relatif ; une arête
fermée (cercle) passe les deux tests de sens et les deux sont essayés. La
distance est ensuite mesurée en projetant `n_samples` points de chaque courbe
sur l'autre (Gauss-Newton amorcé par la correspondance linéaire des
paramètres, recherche grossière sur 33 points en secours). Mesurer **dans les
deux sens** rejette les recouvrements partiels : une arête courte posée sur
une longue passerait dans un seul sens.

Les paires dont les extrémités coïncident mais dont les courbes s'écartent
(deux arcs différents entre les mêmes points) sont listées dans `rejected` ;
elles signalent souvent une erreur de modélisation. Deux arêtes de la **même
face** peuvent être appariées : c'est ainsi qu'une face qui s'enroule sur
elle-même reçoit un seam.

### 3.2 Choix des paires : la plus proche d'abord, résultat variété

Les paires valides sont triées par distance croissante et retenues
gloutonnement, chaque arête au plus une fois : le résultat est toujours
variété. Une paire valide écartée parce qu'une de ses arêtes est déjà cousue
est listée dans `ambiguous` (cas typique : trois faces qui se rencontrent le
long d'une même arête).

### 3.3 Pcurve de l'arête absorbée : composition, pas reparamétrage affine

Le document (§ 6.3) prévoyait de reparamétrer affinement la pcurve de l'arête
absorbée sur l'intervalle du survivant. C'est **faux dès que les deux courbes
n'ont pas le même paramétrage** : un cercle angulaire et un cercle rationnel
de même géométrie ont des paramètres qui ne sont pas liés linéairement. On
calcule donc `q(t) = p₂(φ(t))`, où `φ(t)` est le paramètre de la projection de
`C₁(t)` sur `C₂`, on interpole `q` aux paramètres du survivant, et on mesure la
déviation réelle `|S₂(q(t)) − C₁(t)|` (nœuds et milieux), en doublant
l'échantillonnage jusqu'à atteindre l'écart mesuré entre les courbes ou
`n_pcurve_max`. Comme `p₂` est continue dans l'espace `(u, v)` de sa face, la
nouvelle pcurve l'est aussi, même sur une surface fermée, à condition que `φ` soit elle-même continue et monotone. Les deux extrémités de `φ` sont donc **fixées** (elles ont été appariées à `tol` près) et chaque échantillon intérieur est projeté en partant du précédent. Sans cela, sur une courbe fermée dont début et fin sont le même point 3D, une projection libre peut tomber sur la mauvaise extrémité du paramètre et replier la pcurve ; c'est arrivé sur le runner macOS Intel de la CI, à cause de différences d'arrondi. Le test
`cylinder_closed_by_two_disks` exerce ce cas : disque supérieur bordé par un
cercle rationnel, cylindre analytique.

### 3.4 Fusion des sommets, tolérances

Les sommets aux extrémités appariées sont fusionnés par union-find (le plus
petit id survit, son point n'est pas déplacé), puis toutes les arêtes du modèle
qui pointaient vers un sommet absorbé sont redirigées, comme dans `make_wire`
(PR 3). Les écarts sont absorbés dans les tolérances, sans toucher à la
géométrie (§ 4.3) : `edge.tol = max(tol₁, tol₂, distance mesurée, déviation de
la pcurve)`, puis chaque sommet est agrandi pour contenir les extrémités réelles
des courbes et rester ≥ aux tolérances de ses arêtes.

### 3.5 Orientation

Pour chaque composante connexe, la face de plus petit id est `Forward` ; un
voisin atteint par une arête reçoit le sens qui rend les deux usages de cette
arête opposés (`compose` de la PR 1). Si une face déjà orientée reçoit une
contrainte contraire, la composante n'est pas orientable (ruban de Möbius) : le
shell est quand même produit et `orientable = false`. Le sens **global**
(normales vers l'extérieur) n'est pas décidé ici : c'est le rôle de
`make_solid` (PR 8), qui a besoin du volume signé.

### 3.6 Erreurs et garantie

`sew` n'échoue que sur une entrée invalide (liste vide, tolérance, face morte
ou en double, face sans surface, wire sans pcurve) et laisse alors le modèle
inchangé. Sinon il réussit toujours et décrit le résultat dans `SewReport`,
y compris quand rien n'a pu être cousu.

## 4. Écarts par rapport au document de conception

| Document (§ 6.3) | Implémentation | Raison |
|---|---|---|
| Pcurve reparamétrée affinement | composition `p₂ ∘ φ` puis interpolation, déviation mesurée | exact quel que soit le paramétrage des deux courbes (§ 3.3) |
| `allow_non_manifold` | absent : résultat toujours variété, conflits dans `ambiguous` | un shell non variété n'a pas d'usage au palier 1 |
| `non_manifold_edges` dans le rapport | remplacé par `ambiguous` (paires) | décrit le conflit avec les deux arêtes |
| — | `merged` dans le rapport | traçabilité (survivant, absorbé) |
| Paires internes à une face | autorisées | création de seams |

Limites assumées : pas de découpe d'arête (jonctions en T laissées libres),
pas de fusion d'un sommet sur une arête, arêtes fermées appariées seulement si
leurs points de départ coïncident.

## 5. Tests (`tests_brep_sew`, 8 cas, 57 assertions)

| Test | Couvre |
|---|---|
| `box_from_six_independent_faces` | six faces naturelles, données dans le désordre, trois normales vers l'intérieur : 12 fusions, 24 → 12 arêtes, 24 → 8 sommets, un shell fermé, variété, orienté (3 usages retournés), `check` valide |
| `cylinder_closed_by_two_disks` | cylindre analytique + disque à cercle angulaire + disque à cercle **rationnel** : 2 fusions, shell fermé de 3 arêtes et 2 sommets, `check` valide (SameParameter des pcurves reconstruites compris), aire `(u, v)` de chaque disque égale à celle du cercle (aucun repliement) |
| `components_and_open_shell` | deux boîtes mélangées, l'une ouverte : 2 shells, l'un fermé, l'autre ouvert avec 4 arêtes libres |
| `t_junction_stays_free` | arête longue face à deux arêtes courtes : non cousue, 10 arêtes libres |
| `gap_absorbed_in_tolerances` | écart de 4e-4 : refusé à tol 1e-4, cousu à tol 1e-3 avec tolérances d'arête et de sommets ≥ écart |
| `same_ends_different_curves_rejected` | segments cousus, arcs de mêmes extrémités refusés |
| `moebius_strip_is_not_orientable` | six quadrilatères en ruban de Möbius : 6 fusions, `orientable = false`, `check` signale `ShellNotOrientable` |
| `errors_leave_model_unchanged` | liste vide, tolérance nulle, face morte, face en double |

Les sept suites BREP passent.

## 6. Points à valider

1. **Pas de découpe d'arêtes** : jonctions en T et sommets sur arête laissés libres (question 12).
2. **Distance mesurée dans les deux sens** pour rejeter les recouvrements partiels.
3. **Appariement glouton par distance croissante**, résultat toujours variété, conflits rapportés.
4. **Pcurve reconstruite par composition `p₂ ∘ φ`** au lieu d'un reparamétrage affine.
5. **Paires d'arêtes d'une même face autorisées** (création de seams).
6. **Orientation relative seulement** (face de plus petit id `Forward`) ; le sens extérieur est décidé par `make_solid`.
7. **Le sewing ne lève pas d'erreur** quand rien ne se coud : tout est dans le rapport.
