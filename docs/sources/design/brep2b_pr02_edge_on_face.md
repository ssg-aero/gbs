# Palier 2 bis — PR 2 : arêtes posées sur une face

| | |
|---|---|
| PR | branche `feat/brep2b-edge-on-face` |
| Document de référence | [brep_palier2bis.md](brep_palier2bis.md), § 4 |
| Fichiers | `gbs-brep/section.h`, `gbs-brep/coincident.h`, `gbs-brep/split.h`, `gbs-brep/boolean.h`, `python/gbsBindBrepOps.cpp`, `tests/tests_brep_boolean.cpp`, `python/tests/test_brep_ops.py`, `gbs-occt/tests/tests_brep2_occt.cpp`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Les booléens acceptent qu'une arête d'un opérande soit posée sur une face de
l'autre, sans faces confondues. C'est le cas d'une boîte posée sur une arête,
de deux boîtes qui se touchent par une partie d'arête, ou d'un cylindre couché
sur un plan avec son seam dessous.

## 2. Changements

| Fichier | Changement |
|---|---|
| `section.h` | `SectionResult::edges_on_faces` : chaque arête d'une forme posée sur une face de l'autre, avec cette face. Quand les bords des deux faces se longent, les deux arêtes sont rapportées, chacune avec la face de l'autre |
| `coincident.h` | `imprint_pair` : traversée de deux arêtes (découpe, fusion de sommets), extraite de l'étape 1 de l'empreinte des faces confondues ; `merge_overlapping` : fusion d'une arête avec l'arête de bord qu'elle recouvre (étape 2) |
| `coincident.h` | `imprint_edge(m, e, f, tol, vertex_at)` : l'arête et le bord de la face sont découpés là où ils se traversent ou s'arrêtent l'un sur l'autre ; les morceaux qui longent le bord remplacent les arêtes de bord recouvertes ; les morceaux intérieurs à la face sont rendus comme coupes, avec leur pcurve sur la face |
| `split.h` | `split_at_vertex` rend la nouvelle moitié de l'arête (pour suivre les morceaux) |
| `boolean.h` | après l'empreinte des faces confondues, chaque arête posée sur une face non confondue est imprimée sur elle ; le refus « an edge lies on a face of the other operand » disparaît. Il ne reste que le cas d'une face qui n'est pas une primitive (PR 5) |
| `gbsBindBrepOps.cpp` | `SectionResult.edges_on_faces` |

Après l'empreinte, l'arête de contact est bordée par trois ou quatre faces.
L'appariement radial de la PR 1 sépare alors les deux solides.

## 3. Tests

| Test | Cas | Attendu |
|---|---|---|
| `edges_on_faces` | la boîte tournée de 45° autour de y, posée sur son arête de 4 le long du dessus de l'autre (x = 1, y de −0.5 à 3.5) : l'arête **déborde** de la face | fusion : 2 solides, 48, 13 faces (le dessus coupé en deux) ; commun vide ; différence : 24, 7 faces |
| | boîtes qui se touchent sur une partie d'arête (x = 2, y = 3, z de 1 à 4) | fusion : 2 solides ; commun vide ; différence : 24 |
| | cylindre couché sur la boîte, **seam dessous** | fusion : 2 solides ; commun vide ; différence : 24 |
| `non_generic` | la boîte tournée posée sur le dessus d'un cube à faces B-spline | `NonGeneric` (« an edge lies on a face that is not a primitive »), rien de laissé dans le modèle |
| OCCT `booleans` | `box_on_edge`, `boxes_on_edge`, `seam_on_face` × 3 opérations | valides, volumes à 2e-3, **mêmes nombres de solides et de faces** que `BRepAlgoAPI` |
| Python | boîte posée sur une arête : `edges_on_faces` non vide, fusion en deux solides | |

Les suites `brep`, `intersect`, `step`, `iges`, le binaire OCCT (29 cas) et
les tests Python passent.

## 4. Points à valider

1. **Empreinte limitée aux primitives**, comme pour les faces confondues. La
   pcurve d'un morceau sur la face vient des formes fermées. Les surfaces
   quelconques viendront avec la PR 5.
2. **Sommets posés sur une face** : ils ne découpent rien (coin d'une boîte
   sur une face), comme le prévoit le document.

## 5. Revue de code (a posteriori, lot 5 de la revue du noyau BREP)

`imprint_edge` s'arrêtait sans le dire après 1 000 pas, comme l'empreinte des
faces confondues : l'atteinte de sa borne est désormais un échec
([brep2_pr09b_coincident.md](brep2_pr09b_coincident.md) § 7). Pas d'autre défaut.
