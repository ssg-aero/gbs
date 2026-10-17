# Palier 2 bis — PR 8 : pcurves exactes sur un plan

| | |
|---|---|
| PR | branche `feat/brep-exact-plane-pcurves` |
| Document de référence | [brep2b_pr07_tangent_tolerances.md](brep2b_pr07_tangent_tolerances.md), point à valider 2 |
| Fichiers | `gbs-brep/pcurve.h`, `gbs-brep/sew.h`, `gbs-brep/coincident.h`, `gbs-brep/offset.h`, `tests/tests_brep_boolean.cpp`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

La pcurve d'une arête B-spline sur un plan était projetée puis interpolée. Pour
un cercle rationnel, dont le paramétrage n'est que C¹ à ses nœuds, il fallait
des centaines d'échantillons : 513 pôles pour les couvercles d'un cylindre
B-spline. Sur un plan, cette pcurve peut être **exacte**.

## 2. Changements

| Fichier | Changement |
|---|---|
| `pcurve.h` | `plane_pcurve` : les pôles de la courbe sont envoyés dans les paramètres du plan, S(u, v) = O + u X + v Y. Une courbe rationnelle garde ses coordonnées homogènes (w (P − O)·X, w (P − O)·Y, w). Les nœuds sont les mêmes, donc la pcurve a le paramètre de l'arête. L'écart est la plus grande distance d'un pôle au plan (enveloppe convexe). Rien si un pôle s'écarte de plus de la tolérance |
| `pcurve.h` | `detail::as_bspline` : une `CurveOnSurface` qui suit une iso d'une surface B-spline, comme les arêtes des faces naturelles, redevient une B-spline. On prend la ligne de pôles d'un bord, ou `isoU`/`isoV` à l'intérieur, reparamétrée linéairement (renversée si son paramètre décroît) |
| `pcurve.h` | `project_pcurve` essaie d'abord `plane_pcurve` |
| `sew.h` | la couture, qui recomposait par interpolation (jusqu'à 513 échantillons) la pcurve de l'arête absorbée, prend la pcurve exacte sur un plan |
| `coincident.h` | `edge_pcurve` essaie `plane_pcurve` ; `replace_edge` décale une pcurve B-spline, rationnelle ou non (`translated`) |
| `offset.h` | le décalage d'une pcurve en v accepte une B-spline rationnelle |

Consommateurs des pcurves rationnelles vérifiés : l'écriture STEP et l'export
IGES passent par `to_nurbs`, qui les écrit exactement ; `CurvePieces` les
décompose ; la découpe et la classification ne font que les évaluer.

## 3. Effet

- Couvercles du cylindre B-spline : pcurves rationnelles de **9 pôles** au lieu
  de B-splines interpolées de 513 pôles, exactes à 1e-12.
- `coincident_closed_bsplines` : 7,6 s → 6,1 s. Le reste vient des 20 morceaux
  de cheminement que produit une courbe qui longe le bord d'une face B-spline,
  et de leurs croisements avec les bords : c'est hors du sujet de cette PR.

## 4. Tests

| Test | Cas |
|---|---|
| `exact_plane_pcurves` | couvercles du cylindre B-spline : pcurves rationnelles de 9 pôles, S(pc(t)) = C(t) à 1e-12 en 51 points ; aller-retour STEP sans courbe approchée, volume retrouvé ; face plane bordée de segments B-spline : pcurves de 2 pôles |

Les 475 tests C++, le binaire OCCT (29 cas) et les 51 tests Python passent.

## 5. Points à valider

1. Une `CurveOnSurface` iso est reconnue comme B-spline seulement si sa pcurve
   est un segment exactement à u ou v constant (égalité stricte des
   coordonnées).
2. Les cercles des formules fermées (`EllipseCurve`) sur un plan gardent leur
   pcurve interpolée : leur paramètre est l'angle. Une ellipse 2D exacte
   demanderait une primitive 2D.
