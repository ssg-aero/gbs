# Nettoyage des nœuds à l'export STEP et IGES

| | |
|---|---|
| PR | branche `feat/export-knot-cleanup` |
| Issue | [ssg-aero/gbs#213](https://github.com/ssg-aero/gbs/issues/213) |
| Fichiers | `gbs-io/knot_cleanup.h` (nouveau), `gbs-io/step/geometry_writer.h` (nœuds regroupés par égalité exacte), `gbs-io/step/writer.h`, `gbs-io/iges_brep.h`, `gbs/knotsfunctions.ixx` (coupe sur le nœud), `python/gbsBindBrep.cpp`, `python/gbsbind.cpp`, `tests/tests_knot_cleanup.cpp` (nouveau), `python/tests/test_step.py`, `python/tests/conftest.py`, `docs/index.rst`, `News.md` |
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
et ne change pas. Deux corrections :

- **à la source**, une coupe à moins de `knot_eps` d'un nœud se fait sur le
  nœud (§ 5) : gbs ne fabrique plus ces intervalles ;
- **à l'export**, le nettoyage retire ceux qui sont déjà dans le modèle (un
  fichier importé, un modèle d'une version antérieure), avec une tolérance
  **en longueur** réglable.

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
   jamais. Les fonctions de base dont le support contient le nœud déplacé
   changent : la courbe bouge sur les p + 1 intervalles de chaque côté, de peu
   pour un intervalle court. Les paramètres gardent leur sens, et les pcurves
   restent valides à l'écart près (point 3).
   - À une extrémité, la multiplicité devient degré + 1 + e : les e premières
     (ou dernières) fonctions de base sont nulles sur le domaine, leurs pôles
     sont retirés (exact).
   - À l'intérieur, la multiplicité atteinte est ramenée à celle des deux nœuds
     par suppression de nœud (A5.8), sur toutes les rangées à la fois ; si
     elle reste au-dessus du degré alors que la courbe n'avait pas de
     discontinuité, la suppression est refusée.
3. **Vérification.** Le résultat est comparé à l'**original** (pas à l'étape
   précédente), sur la fenêtre où il a bougé : de p + 1 nœuds avant
   l'intervalle à p + 1 nœuds après. Hors de cette fenêtre, la courbe est celle
   déjà vérifiée. Les deux morceaux sont raffinés à l'union de leurs nœuds ; le
   plus grand écart entre pôles borne l'écart entre les deux courbes au même
   paramètre (enveloppe convexe). Au-delà de `tol`, la suppression est annulée
   et l'intervalle gardé (`kept_spans`).
   - Pour une rationnelle à poids positifs, la borne est celle des pôles
     homogènes, pris autour du centre c des pôles :
     |C1 − C2| ≤ (max |w1 (P1 − c) − w2 (P2 − c)| + max |P2 − c| max |w1 − w2|)
     / min w1, de C1 − C2 = (A1 − A2) / W1 + (C2 − c)(W2 − W1) / W1 avec
     A = Σ N w (P − c).
   - Le balayage reprend p + 1 intervalles avant le dernier retiré : ceux
     d'avant n'ont pas changé. Le coût devient linéaire en nombre
     d'intervalles courts. Quand chaque retrait reprenait tout le balayage et
     raffinait toute la courbe : 595 ms pour une courbe à 499 intervalles
     courts, 1,25 s pour une surface de 50 rangées à 199 intervalles courts.
     Maintenant : 2,1 ms et 22 ms, avec les mêmes écarts. Sans intervalle
     court : 5 ms pour une surface de 200 × 200 pôles.
4. Une surface est nettoyée selon u puis v ; la tolérance est partagée entre
   les deux directions (les écarts s'ajoutent).
5. Un vecteur de nœuds dont les extrémités ne sont pas de multiplicité
   degré + 1 (non serré, comme une courbe périodique lue d'un autre logiciel)
   est laissé tel quel : la suppression à une extrémité suppose cette
   multiplicité.
6. **Pcurves.** Une pcurve est nettoyée dans les paramètres de sa surface,
   multipliés par √2 fois les plus grandes vitesses |S_u| et |S_v| relevées le
   long d'elle (aux nœuds et aux quarts de ses intervalles). Un écart (du, dv)
   déplace le point de la surface d'au plus |S_u| |du| + |S_v| |dv|, que la
   distance des paramètres ainsi multipliés borne. Les vitesses étant
   relevées, l'écart est une estimation. Une pcurve le long de laquelle la
   surface ne bouge pas dans une direction (sur un pôle) est laissée telle
   quelle.

Les bornes sont gardées. `KnotCleanupStats` compte les courbes 3D, les pcurves
et les surfaces nettoyées, les intervalles retirés et gardés, et l'écart
maximal (une borne, sauf pour les pcurves).

Écarté : déplacer le nœud au milieu de l'intervalle (les deux côtés changent),
retirer l'intervalle en reparamétrant la suite (les pcurves ne seraient plus
valides), approximer à nouveau la courbe (perte de l'exactitude là où rien
n'est à corriger).

## 4. Les écrivains

| | STEP (`write_step`) | IGES (`add_brep`) |
|---|---|---|
| Option | `StepWriteOptions::knot_tol` (`std::optional<double>`) : par défaut l'incertitude donnée, sinon `brep_default_tolerance` ; 0 : aucun nettoyage | `IgesExportOptions::knot_tol` : par défaut `brep_default_tolerance` ; 0 : aucun nettoyage |
| Nettoyé | courbes 3D des arêtes, surfaces des faces, pcurves mesurées sur leur surface | courbes 3D des arêtes (après inversion), surfaces des faces, y compris les NURBS exactes des primitives, pcurves mesurées sur leur surface |
| Rapport | `StepWriteReport::knot_cleanup` | `IgesExportReport::knot_cleanup` |
| Tolérance du fichier | incertitude calculée : tolérance typique + 2 × l'écart du nettoyage ; incertitude donnée : `knot_tol` au plus sa moitié | résolution minimale : tolérance maximale + 2 × l'écart du nettoyage |

La tolérance est en unités du modèle (avant le facteur d'échelle de l'IGES).
Une arête et sa face sont nettoyées chacune de leur côté, jusqu'à `knot_tol`
chacune : l'écart entre la courbe 3D et la surface le long de la pcurve peut
croître de deux fois l'écart du nettoyage. Sans incertitude donnée, celle du
fichier en tient compte. Avec une incertitude donnée, qui est la déclaration de
l'appelant, c'est `knot_tol` qui en tient compte : au plus sa moitié, sa
valeur par défaut. Une valeur plus grande est ramenée à cette moitié.

Les pcurves sont nettoyées comme les courbes 3D, leur longueur mesurée sur la
surface (§ 3, point 6) : une pcurve coupée à un ulp d'un de ses nœuds par une
ancienne version, ou importée ainsi, gardait l'intervalle d'un ulp, que les
lecteurs qui fusionnent les nœuds proches rejettent comme les autres.

Python : argument `knot_tol` de `brep.write_step`, `brep.write_step_text`
(`None` : la valeur par défaut) et `IgesWriter.add_brep` ; `KnotCleanupStats`
lisible dans les rapports.

**Actif par défaut.** Un écart inférieur à l'incertitude déclarée est, par
définition, indiscernable pour le lecteur : retirer un intervalle plus court ne
change rien à ce qu'il peut constater, et rend le fichier sûr quel que soit son
réglage de fusion des nœuds, sans que l'utilisateur ait à connaître le
problème. Sans incertitude donnée, l'écrivain STEP la calcule à la fin (médiane
des tolérances plus deux fois l'écart du nettoyage, relevée des écarts des
approximations) ; le nettoyage, qui en a besoin dès le départ, prend
`brep_default_tolerance`, le plancher de cette incertitude. Avec une
incertitude donnée, il en prend la moitié. Les retraits sont comptés dans le
rapport.

Les arêtes isoparamétriques d'une surface B-spline portent ses nœuds : elles
sont nettoyées comme elle (la face de l'issue : la surface et ses deux arêtes
selon v).

## 5. L'écrivain STEP fusionnait déjà les nœuds

En cours de route : la surface de l'issue #213 était écrite avec les nœuds
`(0.9, 1.)` et les multiplicités `(5, 4)`. `GeometryWriter::knots_and_mults`
regroupait les nœuds par `unflat_knots`, qui confond deux nœuds à moins de
`knot_eps`. Le fichier était invalide dès l'écriture, et le lecteur STEP de gbs
le rejetait lui-même (« knot multiplicity 5 out of [1, degree + 1] ») :
SolidWorks n'a pas fusionné les nœuds, il a lu ce que gbs avait écrit.
L'écrivain IGES écrit les nœuds à plat et n'avait pas ce défaut.

`knots_and_mults` regroupe désormais par égalité exacte
(`std::views::chunk_by`) : l'écrivain transcrit la B-spline telle qu'elle est.
Seul, ce correctif ne suffit pas : un lecteur qui fusionne les nœuds proches
retrouverait la multiplicité 5. Le nettoyage par défaut (§ 4) traite ce cas ;
avec `knot_tol = 0`, le fichier écrit `(4,1,4),(0.9,0.9000000000000001,1.)`
et le lecteur de gbs le relit tel quel.

### 5 bis. La coupe faisait l'intervalle

L'intervalle d'un ulp de l'issue vient de la coupe. `split_v(srf, 0.9)` coupe
par `trim_begin` / `trim_end` (`gbs/knotsfunctions.ixx`), qui inséraient le
paramètre exact 0,9 à côté du nœud 0,9000000000000001. `trim_begin` renonçait
déjà à couper à moins de `knot_eps` du premier nœud, mais pas d'un nœud
intérieur. Le `trim` de `to_nurbs` à l'export faisait de même pour une arête
dont le paramètre tombe à un ulp d'un nœud de sa pcurve.

`snap_to_knot` amène désormais la coupe sur le nœud le plus proche à moins de
`knot_eps`, au début comme à la fin : les deux morceaux de l'issue ont les
nœuds `…, 0.5, 0.9000000000000001 ×4` et `0.9000000000000001 ×4, 1 ×4`. Les
bornes demandées restent celles de l'objet (à un ulp des nœuds). Au passage,
`trim_begin` sur un nœud de multiplicité degré + 1 (un coin) donnait
degré + 2 nœuds au début et un pôle de trop : corrigé, le morceau commence au
coin.

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
- `step_and_iges_export` : par défaut, la face écrite en STEP est relue avec
  les nœuds v `0.9 ×4, 1 ×4` ; avec `knot_tol = 0`, les nœuds écrits un à un et
  relus tels quels ; une tolérance donnée nettoie la surface et ses deux arêtes
  selon v ; IGES de même.
- `cut_made_at_the_knot` : les morceaux de la surface de l'issue sans
  intervalle d'un ulp ; une pcurve coupée à un ulp d'un nœud par `to_nurbs` ;
  une coupe sur un coin.
- `pcurve_on_its_surface` : une pcurve à intervalle d'un ulp sur un cylindre
  de rayon 5, nettoyée, son image à moins de la tolérance ; par l'écrivain STEP
  avec sa surface, pas sans.
- `unclamped_left_as_it_is`, `many_short_spans` (199 intervalles de 1e-12,
  tous retirés, l'écart borné).
- `rational_surface` : l'écart relevé est sous la borne rapportée.
- `step_uncertainty_and_knot_tol` : l'incertitude calculée couvre deux fois
  l'écart ; avec une incertitude donnée, l'écart reste sous sa moitié.
- `test_step.test_write_step_knot_tol` : la même chose en Python, sur la
  surface coupée comme avant (par insertion de nœuds), `split_v` donnant
  désormais des morceaux propres.

Les tests Python n'affichent rien sans `GBS_TEST_PLOT_ON` : `conftest.py`
neutralise `plot_curves`, `plot_curves_2d`, `plot_surfaces`, `render_actors`
et `Plotter.show` pour tous les tests, en local comme en CI.

## 7. Revue de code

Onze remarques, toutes traitées dans cette PR.

| # | Remarque | Suite |
|---|---|---|
| 1 | Les pcurves ne sont pas nettoyées, et l'export leur crée le motif de l'issue (coupe à un ulp d'un nœud). | Nettoyées dans les paramètres de leur surface (§ 3, point 6) ; la coupe se fait sur le nœud (§ 5 bis). Tests `pcurve_on_its_surface`, `cut_made_at_the_knot`. |
| 2 | La cause est la coupe, pas corrigée à la source. | `snap_to_knot` dans `trim_begin` et `trim_end` (§ 5 bis). |
| 3 | Une incertitude donnée plus petite que `knot_tol` ne couvre pas l'écart. | `knot_tol` au plus la moitié de l'incertitude donnée. Test `step_uncertainty_and_knot_tol`. |
| 4 | L'écart arête-face peut croître de 2 × `knot_tol`. | Incertitude calculée et résolution IGES relevées de deux fois l'écart. Même test. |
| 5 | Balayage et écart complets refaits à chaque retrait. | Écart sur la fenêtre, balayage repris p + 1 intervalles avant (§ 3, point 3). Test `many_short_spans`. |
| 6 | « Les paramètres hors de [a, b] gardent leurs points » est faux. | Commentaire et note corrigés (§ 3, point 2). |
| 7 | `e = m − (p + 1)` déborde sur un vecteur non serré. | Vecteur non serré laissé tel quel. Test `unclamped_left_as_it_is`. |
| 8 | Pour une rationnelle, l'écart est une estimation. | Borne des pôles homogènes (§ 3, point 3). Test `rational_surface`. |
| 9 | `mult` et `distinct` refont des fonctions existantes. | `gbs::multiplicity`, `std::ranges::unique` ; le test reprend `distinct`. |
| 10 | Interrupteur d'affichage Python réparti test par test. | Dans `conftest.py` seulement. |
| 11 | Regroupement des nœuds à la main. | `std::views::chunk_by`. |
