# Nettoyage des nœuds à l'export STEP et IGES

| | |
|---|---|
| PR | branche `feat/export-knot-cleanup` |
| Issue | [ssg-aero/gbs#213](https://github.com/ssg-aero/gbs/issues/213) |
| Fichiers | `gbs-io/knot_cleanup.h` (nouveau), `gbs-io/step/geometry_writer.h`, `gbs-io/step/writer.h`, `gbs-io/iges_brep.h`, `python/gbsBindBrep.cpp`, `python/gbsbind.cpp`, `tests/tests_knot_cleanup.cpp` (nouveau), `python/tests/test_step.py`, `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif

Une B-spline peut porter des intervalles de nœuds très courts : un nœud à un
ulp d'un autre (issue #213 : un nœud de loft `0.9000000000000001` coupé à
0,9), des nœuds à 1e-13 l'un de l'autre. Pour gbs, la courbe est valide et
exacte. Un logiciel de CAO qui fusionne les nœuds plus proches que sa propre
tolérance y voit une multiplicité supérieure à degré + 1 et écarte l'entité
(SolidWorks écarte la face, sans message).

Chaque logiciel a ses réglages, et aucun ne publie de tolérance de fusion en
paramètre ; un écart en paramètre ne dit d'ailleurs rien de la géométrie. La
tolérance des nœuds de gbs (`knot_eps`) est partagée par tous ses algorithmes
et ne change pas. Le nettoyage se fait **à l'export**, sur demande, avec une
tolérance **en longueur** donnée par l'utilisateur.

## 2. Ce qui existait

- `remove_knot(u, p, num, U, Pw, tol)` (`gbs/knotsfunctions.ixx`) : la
  suppression de nœud bornée du NURBS Book (A5.8), réutilisée telle quelle.
- `insert_knots` : l'insertion (A5.1), réutilisée pour comparer deux courbes.
- `reduce_degree_rows` compare déjà deux courbes en les raffinant à l'union de
  leurs nœuds ; le même principe sert ici à borner l'écart.
- Les deux écrivains BREP passent toutes leurs B-splines par
  `nurbs_export::to_nurbs` (`gbs-io/to_nurbs.h`), puis par
  `GeometryWriter::bspline_curve` / `bspline_surface` (STEP) ou `set_126` /
  `set_128` (IGES).

Aucun nettoyage de nœuds n'existait.

## 3. L'algorithme (`gbs-io/knot_cleanup.h`)

`clean_knots(courbe ou surface, tol, stats)` :

1. **Candidats.** Un intervalle [a, b] est candidat si les points de la
   courbe en a, b et aux quarts sont à moins de `tol` les uns des autres
   (pour une surface : de chaque rangée de pôles, qui encadrent les
   isoparamétriques). Un intervalle court en paramètre mais long dans
   l'espace (la courbe y parcourt une vraie distance) n'est pas touché.
2. **Suppression.** Un des deux nœuds est déplacé sur l'autre, les pôles
   inchangés : b sur a, sinon a sur b ; une extrémité du domaine ne bouge
   jamais. Les paramètres hors de [a, b] gardent donc leurs points, et les
   pcurves restent valides.
   - À une extrémité, la multiplicité devient degré + 1 + e : les e premières
     (ou dernières) fonctions de base sont nulles sur le domaine, leurs pôles
     sont retirés (exact).
   - À l'intérieur, la multiplicité atteinte est ramenée à celle des deux nœuds
     par suppression de nœud (A5.8), sur toutes les rangées à la fois ; si
     elle reste au-dessus du degré alors que la courbe n'avait pas de
     discontinuité, la suppression est refusée.
3. **Vérification.** Le résultat est comparé à l'**original** (pas à l'étape
   précédente) en raffinant les deux à l'union de leurs nœuds : le plus grand
   écart entre pôles borne l'écart entre les deux courbes au même paramètre
   (enveloppe convexe). Au-delà de `tol`, la suppression est annulée et
   l'intervalle gardé (`kept_spans`). Pour une rationnelle, l'écart des pôles
   cartésiens est une estimation, pas une borne.
4. Une surface est nettoyée selon u puis v ; la tolérance est partagée entre
   les deux directions (les écarts s'ajoutent).

Les bornes sont gardées. `KnotCleanupStats` compte les courbes et surfaces
nettoyées, les intervalles retirés et gardés, et l'écart maximal.

Écarté : déplacer le nœud au milieu de l'intervalle (les deux côtés changent),
retirer l'intervalle en reparamétrant la suite (les pcurves ne seraient plus
valides), approximer à nouveau la courbe (perte de l'exactitude là où rien
n'est à corriger).

## 4. Les écrivains

| | STEP (`write_step`) | IGES (`add_brep`) |
|---|---|---|
| Option | `StepWriteOptions::knot_tol` (`std::optional<double>`, absent par défaut) | `IgesExportOptions::knot_tol` (0 par défaut) |
| Nettoyé | courbes 3D des arêtes, surfaces des faces | courbes 3D des arêtes (après inversion), surfaces des faces, y compris les NURBS exactes des primitives |
| Rapport | `StepWriteReport::knot_cleanup` | `IgesExportReport::knot_cleanup` |
| Tolérance du fichier | l'incertitude écrite couvre l'écart du nettoyage | la résolution minimale couvre l'écart du nettoyage |

La tolérance est en unités du modèle (avant le facteur d'échelle de l'IGES).
Les pcurves ne sont pas nettoyées : elles vivent dans l'espace des paramètres
de la surface, où une longueur n'a pas de sens.

Python : argument `knot_tol` de `brep.write_step`, `brep.write_step_text` et
`IgesWriter.add_brep` ; `KnotCleanupStats` lisible dans les rapports.

Désactivé par défaut : sans demande, un B-spline est écrit tel quel, comme
avant.

## 5. Constat en cours de route : l'écrivain STEP fusionne déjà les nœuds

Sans `knot_tol`, la surface de l'issue #213 est écrite avec les nœuds
`(0.9, 1.)` et les multiplicités `(5, 4)` : `GeometryWriter::knots_and_mults`
regroupe les nœuds par `unflat_knots`, qui confond deux nœuds à moins de
`knot_eps`. Le fichier est donc invalide dès l'écriture, et le lecteur STEP de
gbs le rejette lui-même (« knot multiplicity 5 out of [1, degree + 1] »).
SolidWorks n'a pas fusionné les nœuds : il a lu ce que gbs avait écrit.
L'écrivain IGES écrit les nœuds à plat et n'a pas ce défaut.

Point ouvert, à trancher avant fusion :

- regrouper par égalité exacte dans `knots_and_mults` (le fichier décrit
  fidèlement la B-spline ; un lecteur strict peut encore fusionner) ;
- ou activer le nettoyage par défaut, avec l'incertitude du fichier comme
  tolérance.

## 6. Tests

- `tests_knot_cleanup.curve_span_of_an_ulp_at_the_start` : la courbe coupée
  comme dans l'issue, extrémité ramenée à degré + 1, même courbe au même
  paramètre, rien à refaire au second passage.
- `curve_close_interior_knots` : deux nœuds intérieurs à 1e-13, fusionnés,
  écart < tol.
- `long_spans_kept` : intervalles longs, intervalle court en paramètre mais
  long dans l'espace, tolérance nulle : rien ne change.
- `surface_of_issue_213` : la surface de l'issue après la coupe.
- `rational_surface` : surface rationnelle, nœuds proches selon u.
- `step_and_iges_export` : désactivé par défaut ; avec `knot_tol`, la face
  écrite en STEP est relue avec les nœuds v `0.9 ×4, 1 ×4` ; IGES.
- `test_step.test_write_step_knot_tol` : la même chose en Python.
