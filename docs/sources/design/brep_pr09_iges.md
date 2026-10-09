# BREP palier 1 — PR 9 : export IGES (`gbs-io/iges_brep.h`)

| | |
|---|---|
| PR | branche `feat/brep-iges` |
| Document de référence | [brep_core.md](brep_core.md), section 7, question 10 ; notes précédentes [PR 1](brep_pr01_model.md) à [PR 8](brep_pr08_solid.md) |
| Fichiers | `gbs-io/iges_brep.h` (306 l., nouveau), `gbs-io/iges.h` (`#pragma once`, accesseur `IgesWriter::model()`), `tests/tests_brep_iges.cpp` (nouveau), `gbs-occt/tests/tests_brep_iges_occt.cpp` (nouveau, relecture par OCCT), `gbs-occt/tests/CMakeLists.txt` (libIGES, NLopt), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Écrire en IGES les faces, shells, solides et compounds du noyau natif, avec
libIGES déjà en dépendance, sous une forme que les autres systèmes relisent
comme des surfaces découpées. C'est le remplaçant natif de
`occt_utils::to_iges` pour les formes BREP.

## 2. Vue d'ensemble

```
add_brep(DLL_IGES& | IgesWriter&, model, shape, name = "", IgesExportOptions{approx_tol = 1e-5, scale = 1}) -> IgesExportReport

pour chaque face de explore<FaceId>(shape) :
  128  surface NURBS                      exacte si BSSurface / BSSurfaceRational, sinon interpolée (paramètres conservés)
  144  surface découpée                   PTS = 128, PTO = 142 du contour extérieur, PTI = 142 des trous, label = name + rang
  142  courbe sur surface (par wire)      SPTR = 128, BPTR = 102 des pcurves, CPTR = 102 des courbes 3D, PREF = 3 (les deux)
  102  courbe composite                   un segment 126 par co-arête, dans l'ordre et le sens du wire
  126  courbe NURBS                       exacte si BSCurve / BSCurveRational (tronquée à [u1, u2]), sinon interpolée ;
                                          inversée pour une co-arête Reversed ; 2D écrites avec z = 0 ;
                                          arête dégénérée : 126 de degré 1 réduite au point du sommet
arêtes de shape qui ne bordent aucune face exportée : 126 seules

IgesExportReport { faces, wires, free_edges, approximated_curves, approximated_surfaces, max_deviation }
```

## 3. Choix d'architecture et justification

### 3.1 Faces en 144, pas d'entités BREP 186

Conformément à la question 10, un shell, un solide ou un compound est écrit
comme l'ensemble de ses faces. Les entités BREP d'IGES (186 solide, 514 shell,
510 face, 508 loop, 504/502 listes) existent dans le cœur de libIGES mais ne sont
pas exposées par son API `DLL_IGES`. Les lecteurs usuels (OCCT, CATIA, NX)
reconstruisent la topologie par sewing ; le test de relecture par OCCT le
vérifie. L'orientation des `FaceUse` n'est pas transmise (une 144 n'a pas de
sens) ; le lecteur la refait en cousant.

### 3.2 Géométrie exacte quand c'est possible, paramétrage conservé sinon

- Les B-splines et NURBS gbs sont écrites **exactement**. libIGES attend, pour
  chaque pôle, les **coordonnées cartésiennes puis le poids** (`x, y, z, w`),
  alors que gbs stocke les pôles rationnels sous forme **homogène**
  (`w·x, w·y, w·z, w`) : l'export divise donc par le poids. La relecture par OCCT
  l'a établi : en passant les pôles homogènes, le cylindre rationnel relu était
  déformé (aire 15 % trop faible) et ne se recousait plus aux disques. Avec la
  conversion, le volume relu est exact à 4e-9.

  Le code d'export existant de `gbs-io/iges.h` (`add_geom`, utilisé par
  `IgesWriter::add_geometry`) passe les pôles homogènes et a donc ce défaut pour
  toute géométrie rationnelle, ainsi qu'un pas de coefficients faux pour les
  courbes 2D. Il n'est **pas** corrigé dans cette PR (hors périmètre) ; une
  correction séparée est proposée.
- Une courbe 3D est **tronquée** à l'intervalle de son arête (`trim`) plutôt que
  de compter sur les paramètres `V(0), V(1)` de l'entité 126, que tous les
  lecteurs n'honorent pas.
- Toute autre courbe (`CurveOnSurface`, `Line`, courbe analytique) ou surface
  (`SurfaceOfRevolution`, `SurfaceOffset`, analytique) est **interpolée en
  B-spline cubique à ses propres paramètres**, avec doublement de
  l'échantillonnage jusqu'à `approx_tol` (1025 points par courbe, grille 129 × 129
  au plus). Conserver le paramétrage de la surface garde les pcurves valides
  telles quelles ; le rapport donne le nombre d'approximations et l'écart
  maximal.
- La surface de révolution n'est **pas** écrite en entité 120 comme le fait
  `IgesWriter::add_geometry` : la convention de paramètres de la 120 (angle,
  génératrice) diffère selon les lecteurs, et une approximation de la
  génératrice changerait son paramétrage, donc invaliderait les pcurves.

### 3.3 Sens et ordre des segments

Une courbe composite IGES doit s'enchaîner tête-bêche. Chaque co-arête
`Reversed` est donc écrite **inversée** (`reverse()` de gbs), en 2D comme en
3D ; les segments suivent l'ordre du wire, qui est celui des PR 3 et 5
(contour extérieur direct, trous indirects en `(u, v)`).

### 3.4 Échelle et unités

`scale` multiplie les coordonnées de l'espace modèle (pôles des surfaces et des
courbes 3D, poids inchangés) et jamais les paramètres ni les pcurves. Il
remplace l'échelle fixe 1000 de `gbs-occt/export.h`. L'unité de la section
globale reste celle du modèle libIGES (`IgesWriter::model().SetUnitsFlag(...)`) :
gbs n'impose pas d'unité (question 11).

### 3.5 Pas de dépendance nouvelle

`iges_brep.h` est header-only, dans `gbs-io` à côté de `iges.h`, et n'utilise que
libIGES. `IgesWriter` gagne seulement un accesseur `model()` ; la surcharge
`add_brep(IgesWriter&, …)` évite d'alourdir `iges.h` avec `gbs-brep`.

## 4. Écarts par rapport au document de conception

| Document (§ 7) | Implémentation | Raison |
|---|---|---|
| `IgesWriter::add_geometry(model, shape, name)` | fonction libre `add_brep(writer, model, shape, name, opts)` + `IgesWriter::model()` | `iges.h` reste sans dépendance à `gbs-brep` |
| `N1 = natural_bounds ? 0 : 1` | toujours un contour extérieur (N1 = 1) | toujours valide ; l'optimisation n'apporte rien aux lecteurs |
| `CRTN = 1` si pcurve projetée | `CRTN = 0` (non spécifié) | l'origine de la pcurve n'est pas mémorisée dans le modèle |
| Noms par 406 forme 15 | label de l'entrée de répertoire (`SetLabel`), comme `iges.h` | même mécanisme que l'existant ; 8 caractères |
| 120 pour `SurfaceOfRevolution` | 128 interpolée | pcurves valides, conventions de la 120 variables (§ 3.2) |

## 5. Tests

`tests/tests_brep_iges.cpp` (3 cas, 27 assertions) — le fichier est relu par
`DLL_IGES::Read` et ses entités sont comptées dans la section « Directory
Entry » :

| Test | Couvre |
|---|---|
| `box_solid` | boîte cousue : 6 × 144, 128, 142 ; 12 × 102 ; 48 × 126 ; aucune surface approchée ; écart ≤ `approx_tol` |
| `curved_and_trimmed_faces` | sphère analytique (seule surface approchée), cylindre NURBS exact, face trouée depuis des boucles 2D, arête libre : 3 × 144, 4 × 142, 8 × 102, 31 × 126 |
| `rational_solid_and_scale` | cylindre rationnel fermé par deux disques rationnels, échelle 1000 : 3 × 144 et 128, aucune approximation de surface |

`gbs-occt/tests/tests_brep_iges_occt.cpp` (2 cas, 12 assertions, job OCCT de la CI) — relecture
par `IGESControl_Reader`, `BRepCheck_Analyzer` sur chaque face,
`BRepBuilderAPI_Sewing` d'OCCT, `BRepGProp::VolumeProperties` :

| Test | Couvre |
|---|---|
| `box_round_trip` | 6 faces valides, shell fermé, volume 1 à 1e-6 |
| `rational_cylinder_round_trip` | surface et cercles rationnels, échelle 10 : 3 faces valides, shell fermé, volume `π h × 1000` à 1e-6 (mesuré : 4e-9) |

Les volumes OCCT sont calculés par `BRepGProp::VolumeProperties` **adaptatif**
(tolérance 1e-9) : l'ordre de Gauss fixe par défaut donne 0,8 % d'erreur sur
des faces rationnelles, ce qui masquerait une erreur d'export. Les 8044
assertions de `gbs-occt_tests` passent localement.

## 6. Points à valider

1. **Conversion des pôles rationnels** homogènes → cartésiens + poids à l'écriture, et correction séparée de `gbs-io/iges.h`.
2. **Faces en 144** pour tous les shapes, sans entités BREP 186 (question 10).
3. **Géométrie non NURBS interpolée à ses propres paramètres** (cubique,
   `approx_tol` 1e-5 par défaut), y compris `SurfaceOfRevolution` (pas de 120).
4. **Courbes 3D tronquées** à l'intervalle de l'arête plutôt que bornées par
   `V(0), V(1)`.
5. **`add_brep` en fonction libre** et `IgesWriter::model()` exposé.
6. **Toujours un contour extérieur** (N1 = 1), `CRTN = 0`, label de répertoire.
7. **Orientation des faces non transmise** : le lecteur recoud.
