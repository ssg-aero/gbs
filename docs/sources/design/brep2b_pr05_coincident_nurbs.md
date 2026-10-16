# Palier 2 bis — PR 5 : faces confondues quelconques

| | |
|---|---|
| PR | branche `feat/brep2b-coincident-nurbs` |
| Document de référence | [brep_palier2bis.md](brep_palier2bis.md), § 7 (précisé par cette PR) |
| Fichiers | `gbs-brep/section.h`, `gbs-brep/coincident.h`, `gbs-brep/boolean.h`, `gbs-brep/classify.h`, `gbs/intersect/marching.h`, `tests/tests_brep_boolean.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/design/brep_palier2bis.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Les booléens traitent les faces confondues entre surfaces quelconques : deux
B-splines, ou une primitive et une B-spline de même forme. Jusqu'ici, seules
les primitives étaient reconnues confondues, par leurs formes fermées. Les
arêtes posées sur une face B-spline (PR 2) en profitent aussi.

## 2. Changements

| Fichier | Changement |
|---|---|
| `section.h` | `on_surface`, `coincident_faces` : points intérieurs d'une face projetés sur la surface de l'autre, dans les deux sens ; décision du § 7.1 du document, y compris les faces qui se touchent par leurs bords (intersection singulière) |
| `coincident.h` | `edge_pcurve` : projection (`project_pcurve`) hors des primitives ; `periods` : périodes lues sur la fermeture de la surface ; une arête sans pcurve parce qu'elle sort de la surface est ignorée si son milieu est hors de la face |
| `boolean.h` | `same_side` projette sur toute surface ; messages d'erreur mis à jour |
| `marching.h` | courbes le long du bord d'une surface B-spline (§ 3) |
| `classify.h` | évaluation ramenée dans les bornes (§ 3) |

## 3. Défauts trouvés en chemin

Les booléens de boîtes B-spline font longer à leurs courbes d'intersection
les bords des surfaces, ce que les primitives (planes, infinies) ne faisaient
jamais :

| Défaut | Correction |
|---|---|
| les pcurves interpolées d'une branche le long d'un bord dépassaient le rectangle des paramètres, et l'évaluation de la B-spline levait une exception | pôles ramenés dans le rectangle (enveloppe convexe), comme les pcurves projetées de `pcurve.h` |
| un paramètre au-delà d'une borne d'un arrondi (le paramètre fixe du bord) faisait prendre le point pour hors de la surface | ramené sur la borne (1e-10 de l'étendue) |
| la pose sur une boîte choisissait un paramètre hors de sa boîte d'un simple arrondi, avec un rapport négatif | ignoré : seuls comptent les paramètres vraiment sortis, avec un rapport positif |
| `classify` (PR 4) admettait des paramètres dans la marge d'arrondi des bornes, puis évaluait la surface au-delà | évaluation ramenée dans les bornes |

Ces branches finissaient `singular`, et les booléens de boîtes B-spline
étaient refusés.

## 4. Tests

| Test | Cas | Attendu |
|---|---|---|
| `coincident_bsplines` | cubes à faces B-spline bilinéaires décalés (fonds et dessus confondus), empilés (une face contre l'autre, flancs coplanaires qui se touchent) ; boîte à plans contre cube B-spline (fonds sur le même plan) | volumes des trois opérations ; 10 faces pour les cubes empilés fusionnés |
| `non_generic` | cylindre couché sur un cube B-spline : tangence le long d'une ligne trouvée par cheminement | toujours `NonGeneric` ; rien laissé dans le modèle. L'ancien cas (arête sur une face B-spline) est maintenant traité |
| OCCT `booleans` | `coplanar_nurbs`, `stacked_nurbs` (boîtes converties en NURBS par OCCT), `plane_and_nurbs` | valides, volumes, **mêmes nombres de solides et de faces** que `BRepAlgoAPI` |

Les suites `brep`, `intersect`, `step`, `iges`, le binaire OCCT (29 cas) et
les tests Python passent. La relecture par gbs d'une forme d'OCCT affiche
maintenant les problèmes du contrôle quand il échoue.

## 5. Limites

- **Cylindres convertis en NURBS par OCCT** : leurs pcurves sont à 5e-6 de
  leurs arêtes, au-delà de la tolérance du contrôle à la relecture par gbs.
  Le cas a été retiré de la comparaison : c'est la qualité de l'entrée qui
  est en cause, pas le booléen.
- **Coïncidence partielle** de deux surfaces quelconques : refusée (question 4
  du document).
- **Coût** : le test de coïncidence projette jusqu'à 72 points par paire de
  faces dont l'une n'est pas une primitive. Il s'arrête au premier point qui
  contredit la coïncidence, ce qui est le cas courant.

## 6. Revue de code

| Commentaire | Suite donnée |
|---|---|
| `section.h` : « pas de point en vis-à-vis + intersection singulière ⇒ confondues » peut absorber des surfaces différentes tangentes, et jette les courbes trouvées | règle restreinte (`touching_along_borders`) : une arête de bord d'une face longe l'autre surface (au moins 3 de 9 échantillons) **et** toutes les courbes trouvées sont sur un bord ; sinon l'intersection est gardée et rapportée `singular` |
| `coincident.h` : `edge_pcurve` perd le code d'erreur de `project_pcurve` ; une arête qui traverse le seam d'une B-spline fermée est ignorée | `edge_pcurve` rend `std::expected<PCurveFit, PCurveErrc>`. `CrossesSeam` : l'arête est **découpée au seam** (`seam_crossing`, `split_at_seam`) ; seul `OffSurface` avec un milieu hors de la face est ignoré (`cuts_on`, commun à l'empreinte des faces et des arêtes) |
| `marching.h` : la pose rejette les rapports négatifs, un point courant hors de sa boîte d'un arrondi finit `singular` | rapport ramené dans [0, 1] |
| `marching.h` : bloc de pose dupliqué entre `step_over` et `march` | `Marcher::land` |
| `marching.h` : clamp des pôles dupliqué avec `pcurve.h` | `clamp_poles` dans `gbs/bscinterp.h`, utilisé par les deux |
| `section.h` : 72 `classify` par paire avant toute projection | projection d'abord ; `classify` seulement pour les points qui ne tombent pas sur le bord de l'autre boîte |
| `section.h` : sentinelles entières | `on_surface` rend `std::optional<int>`, `coincident_faces` une `enum class Coincidence { No, Yes, Undecided }` |
| `boolean.h` : `same_side` réévalue la surface ; `surface_closure` recalculée par arête | distance de la projection réutilisée ; fermeture calculée une fois par face (`cuts_on`, `replace_edge`) |
| tests : ni B-spline fermée, ni seam, ni recouvrement partiel | `coincident_closed_bsplines` : cylindre dont le flanc est une B-spline rationnelle fermée (disque bordé par le cercle rationnel, extrudé), contre lui-même, décalé et tourné (seams différents, cercles qui traversent le seam de l'autre), empilé (contact par le bord), contre un cylindre primitif de même rayon, et une boîte posée dessus |

### 6.1 Défauts trouvés par ces tests

| Défaut | Correction |
|---|---|
| `intersect_detail::project` (point sur courbe) s'arrêtait sur l'échantillon le plus proche quand le dénominateur de Newton était négatif (paramétrage rationnel), et, sur une courbe fermée, partait de la mauvaise extrémité contre sa borne : un point du cercle de bord était à 0,39 de lui, et le cercle commun d'un couvercle et d'un flanc ressortait en 110 arêtes de section | 16 échantillons, pas de Gauss-Newton quand Newton est indéfini, Newton depuis les **deux** meilleurs échantillons |
| `project_on_surface` : même défaut sur une surface fermée (nœud de grille du mauvais côté du seam) | Gauss-Newton depuis les deux meilleurs nœuds |
| la pcurve d'un cercle primitif sur un flanc B-spline rationnel ne descend pas sous 3e-6 (paramétrage rationnel seulement C¹ aux nœuds) | pcurve projetée à la tolérance d'approximation des pcurves (1e-5) ; la tolérance de l'arête et de ses sommets est relevée à l'écart (`cover_deviation`) |
| un point à moins de la tolérance de deux surfaces perpendiculaires est jusqu'à √2 tol de leur bord commun | morceaux de section classés à 2 tol |

## 7. Points à valider

1. La règle du § 7.1 pour les faces qui se touchent par leurs bords
   (intersection singulière sans point en vis-à-vis : confondues).
2. Le retrait du cas des cylindres convertis par OCCT.
