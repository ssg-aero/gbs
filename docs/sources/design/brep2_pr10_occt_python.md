# Palier 2 — PR 10 : validation par OpenCascade, Python, documentation

| | |
|---|---|
| PR | branche `feat/brep2-occt-python` |
| Document de référence | [brep_palier2.md](brep_palier2.md), § 8 et 9 ; bilan au § 12 |
| Fichiers | `gbs-occt/tests/tests_brep2_occt.cpp` (nouveau), `python/gbsBindBrepOps.cpp` (nouveau), `python/gbsbind.cpp`, `python/CMakeLists.txt`, `python/tests/test_intersect.py`, `python/tests/test_brep_ops.py` (nouveaux), `gbs/intersect/marching.h`, `gbs-brep/sweep.h`, `ReadMe.md`, `docs/sources/design/brep_palier2.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Clore le palier 2 :
- comparer ses opérations à celles d'OpenCascade ;
- les rendre utilisables en Python ;
- documenter le bilan.

## 2. Validation par OpenCascade (job OCCT de la CI)

Les formes sont construites par OCCT et passées à gbs par STEP. Les résultats
de gbs sont écrits par gbs, puis relus par OCCT.

| Test | Couvre | Critères |
|---|---|---|
| `booleans` | 7 paires × 3 opérations : boîtes, boîte et cylindre, boîte et sphère inclinée, cylindres croisés, boîtes coplanaires, boîtes empilées, cylindres coaxiaux | `check()` de gbs ; relu par OCCT : `BRepCheck_Analyzer` valide, volume à 2e-3 de `BRepAlgoAPI`, mêmes nombres de solides et de faces |
| `sections` | boîte, cylindre, sphère inclinée par des plans obliques | longueur à 1e-4 de `BRepAlgoAPI_Section`, wires fermés |
| `sweeps` | plaque percée construite par OCCT : extrusion en biais, révolution complète, révolution de 2 rad | comme les booléens, contre `BRepPrimAPI_MakePrism` et `MakeRevol` |

Tous passent. Les 21 booléens donnent **le même nombre de faces** qu'OCCT, ce
qui montre que la découpe et la sélection produisent la même topologie.

## 3. Corrections trouvées par la comparaison

| Défaut | Effet | Correction |
|---|---|---|
| une branche de cheminement réduite à un point (un départ sur les seams des deux surfaces à la fois) était rapportée comme contact | les cylindres d'OCCT, dont les seams se croisent sur la courbe d'intersection, étaient refusés (`NonGeneric`) | seul un point où les surfaces sont **tangentes** est un contact |
| révolution complète d'une face trouée : deux shells fermés | `sweep` exigeait un seul shell | le plus grand shell ferme le solide, les autres en sont les **cavités** (anneau à cavité torique) |

## 4. Python

Le nouveau fichier `python/gbsBindBrepOps.cpp` a deux parties.

**Sous-module `gbs.intersect`**

| Fonction | Rôle |
|---|---|
| `curves(c1, c2, range1, range2, tol, …)` | 2D ou 3D : `hits`, `overlaps` |
| `curve_surface(curve, surface, range, box, …)` | `hits`, `portions` |
| `surfaces(s1, s2, box1, box2, …, closed_forms)` | `curves` (courbe, pcurves, écart), `points`, `coincident`, `singular` |

**Fonctions ajoutées à `gbs.brep`**

| Fonction | Rôle |
|---|---|
| `classify` (face et uv, ou solide et point), `classify_face`, `interior_point` | classification |
| `section` (deux formes, ou forme et surface) | `SectionResult` |
| `split_face(model, face, [(edge, pcurve)])`, `split` (par une forme ou une surface) | découpe |
| `boolean(model, op, a, b)`, `fuse`, `cut`, `common` | `BooleanResult`, avec l'origine de chaque face |
| `make_prism`, `make_revolution(model, profile, origin, direction, angle)`, `offset_face` | balayages et faces décalées |

Les erreurs lèvent `gbs.brep.BRepError`, comme au palier 1 : configuration non
générique, profil invalide, forme morte.

## 5. Tests Python

| Test | Couvre |
|---|---|
| `test_intersect.py` | segments croisés, droite et cylindre, droites 3D avec plage ; cylindre et plan : cercle fermé posé sur les deux surfaces par ses pcurves, même paire par cheminement |
| `test_brep_ops.py` | classification (dedans, sur une face, dehors, point intérieur) ; section d'une boîte par un plan ; fusion, intersection, différence de boîtes (volumes, origines des faces) ; sphère tangente refusée ; découpe par un plan ; extrusion de la plaque percée, demi-révolution, faces décalées du cylindre |

Les 25 tests Python du BREP, de STEP, des primitives et des intersections
passent.

## 6. Documentation

- Le `ReadMe.md` montre en Python les booléens, la section et les
  intersections.
- `brep_palier2.md` passe au statut **réalisé**, avec un bilan au § 12 :
  - les PR et leurs notes ;
  - la validation par OCCT ;
  - les décisions révisées ;
  - les limites connues.

## 7. Points à valider

1. **Critères de la comparaison avec OCCT** : validité, volumes, nombres de
   solides et de faces, longueurs de section.
2. **Contact = point tangent seulement** (correction du cheminement).
3. **Cavités** des révolutions complètes de faces trouées.
4. **API Python** : sous-module `gbs.intersect`, et opérations dans `gbs.brep`.

## 8. Revue de code (a posteriori, lot 5 de la revue du noyau BREP)

Les liaisons Python des booléens rendent les erreurs comme exceptions
(`value_or_raise`) : avec le retour arrière de `boolean`, un échec laisse aussi
le modèle Python intact. Pas de défaut trouvé. La comparaison avec OpenCascade
n'a pas été relancée localement.
