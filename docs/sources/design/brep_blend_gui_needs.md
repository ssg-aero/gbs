# Congés et chanfreins : manques relevés par une interface graphique

| | |
|---|---|
| PR | à venir, branche `fix/brep-blend-gui-needs` (sur master f191406) |
| Documents de référence | [brep_palier4.md](brep_palier4.md) (congés et chanfreins), [brep_palier7_blend_sections.md](brep_palier7_blend_sections.md) |
| Fichiers | `gbs-brep/sweep.h`, `gbs-brep/pipe.h` (révolution : disques plans), `gbs-brep/closure.h`, `gbs-brep/check.h`, `gbs-brep/blend.h` (pcurves dans les bornes de leur face), `gbs-brep/model.h` (`RefusalReason`, `BuildError::reason`, `BuildError::limit`), `gbs-brep/fillet.h`, `gbs-brep/blend.h` (raisons des refus, chanfrein à une face par arête), `tests/tests_brep_sweep.cpp`, `tests/tests_brep_boolean.cpp`, `tests/tests_brep_fillet.cpp`, `News.md`, cette note |
| Statut | À valider (§ 6) |

Une interface graphique de CAO bâtie sur `gbs-brep` (congés et chanfreins d'arêtes
choisies à la souris) a relevé cinq manques, numérotés ici comme dans sa liste :
G15 et G18 (une même cause), G19, G16 et G17. Chacun a son commit.

## 1. G15 et G18 : les disques d'une révolution

### 1.1 Le défaut

`make_revolution` d'un profil fermé qui touche l'axe (un rectangle dont un côté est
sur l'axe : un cylindre ; un trapèze : un tronc de cône) rendait un flanc analytique
(`CylindricalSurface`, `ConicalSurface`), mais ses disques en B-spline rationnelle
fermée non périodique (`build_revolution` du côté perpendiculaire à l'axe), avec une
couture radiale et un pôle au centre. Deux suites :

- **G15** (issue #257). Un booléen qui garde un morceau de disque que la couture traverse de bord à
  bord (cylindre R 25 en (30, 0) moins cylindre R 25 en (0, 0)) le rend en deux faces
  sur la même surface ; `unify_faces` ne les réunit pas (il ne recentre la fenêtre
  d'un paramètre périodique que pour une surface élémentaire).
- **G18.** Un congé ou un chanfrein d'un cercle du cylindre est refusé : « a contact
  crossing the seam of a face that is not periodic (a closed B-spline) ». Le même
  cylindre lu en STEP (disques plans) passe.

### 1.2 Le correctif

Le premier des deux correctifs proposés : un côté du profil perpendiculaire à l'axe
balaie une région d'un **plan** (`sweep_detail::Annulus` : rayons r1 < r2, r1 = 0 sur
l'axe ; angle balayé), et non plus une B-spline :

| Tour | r1 = 0 | r1 > 0 |
|---|---|---|
| complet | disque : un cercle | couronne : deux cercles (le petit en trou) |
| partiel | secteur : un arc, deux rayons | secteur de couronne : deux arcs, deux rayons |

Le plan a pour repère l'origine sur l'axe à la hauteur du côté, z l'axe, x vers le
côté (le repère du flanc voisin, dont le cercle commence au même point) ; ses
bords sont des `EllipseCurve` et des segments. La face est faite par
`make_face(plane, wire, holes)`. Plus de couture, plus de pôle : le booléen garde
une face, le congé et le chanfrein ne traversent plus rien.

Le second correctif proposé (un congé qui traverse la couture d'une B-spline fermée
non périodique) reste à faire pour les révolutions d'autres profils (une courbe qui
touche l'axe hors d'un côté perpendiculaire) : hors de ce lot.

### 1.3 L'exactitude des volumes

Le cercle d'un plan, paramétré par l'angle, n'a pas de pcurve B-spline exacte (la
B-spline rationnelle du cercle n'a pas le paramètre de l'angle). Deux mesures
gardent les intégrales de volume et d'aire exactes à 1e-11 près (les tests de
`full_round.h` sur des révolutions les comparent à Pappus à 1e-9) :

- les pcurves des cercles sur le plan sont interpolées à `pcurve_tol` = 1e-10 r2
  (1 025 points au plus : l'erreur de l'interpolation cubique d'un tour entier est
  d'environ 2e-11 r2) ;
- les faces planes sont construites avant les flancs : leurs arêtes ont les plus
  petits identifiants et survivent à `sew` (le survivant garde sa courbe, la pcurve
  de l'autre face est recomposée). Le flanc, dont l'arête de cercle a le même
  paramètre (l'angle, depuis le même point), garde une pcurve droite.

Sans la seconde, `sew` gardait l'iso du cylindre et recomposait la pcurve du plan à
sa propre tolérance : un écart de volume de 3e-7 (le fond d'une gorge).

### 1.4 Invariants

- `make_revolution` ne rend plus de `BSSurfaceRational` pour un côté perpendiculaire à
  l'axe ; les autres côtés gardent leurs surfaces (cylindre, cône, sphère, tore,
  révolution NURBS).
- Le nombre de faces d'une révolution ne change pas (une face par côté balayé, plus
  les deux copies du profil pour un tour partiel).
- `make_pipe` (balayage le long d'arcs) suit la même règle : il appelle
  `revolution_surface`.

## 2. G19 : une pcurve à peine hors du domaine de sa face

### 2.1 Le défaut

Sur une pièce à B-splines publique (`tests/in/step/occt_loft.stp`), le congé ou le
chanfrein d'une arête latérale, au rayon de 1 % de la boîte, était refusé par
« Surface U eval 4.000001 out of bounds [0, 4] » (sept cas sur les quatorze arêtes et
opérations possibles).

### 2.2 La cause

Le contact d'un congé, et la courbe qui le termine sur la face suivante, sont
calculés sur le **support** de leur face, la surface prolongée au-delà de ses bornes
(et au-delà des extrémités de l'arête : la pcurve couvre le congé prolongé). Là où
ils finissent sur un bord d'une face B-spline, leur pcurve dépasse ce bord de la
tolérance (de 1e-13 à 6e-6 en paramètre). Toute évaluation de la surface de la
**face** sur cette pcurve, qui refuse un paramètre hors bornes, levait l'exception :

- `check` (`Checker::same_parameter`), à la fin du congé ;
- `blend_detail::settle` (tolérance d'une arête nouvelle), dans un chanfrein ;
- le maillage du résultat (`mesh`), quand `check` était passé (un congé de 3 %).

Et la tolérance de l'arête de contact était mesurée sur le support seulement, pas sur
la face qu'elle borde.

### 2.3 Le correctif

- **Les pcurves dans les bornes de leur face.** Avant le contrôle du résultat, chaque
  pcurve d'une face qui n'est pas une primitive est rognée à la plage de son arête
  (les paramètres gardés) et ses pôles ramenés sur les bornes de la surface, quand
  ils les dépassent d'au plus 1e-5 de leur étendue (`blend_detail::onto_face_bounds`,
  comme `clamp_poles` pour une pcurve projetée) ; la tolérance de l'arête couvre
  l'écart (`settle`). Une B-spline est dans l'enveloppe convexe de ses pôles : la
  pcurve reste dans les bornes, et tout code qui évalue la face sur elle (maillage,
  contrôle, intégrales) aussi. Les primitives ne sont pas touchées : leurs
  directions ouvertes ne sont pas bornées (#260), leurs directions fermées tournent.
  Un pôle plus loin (une pcurve vraiment hors de la face) est laissé : le contrôle le
  refuse.
- `detail::value_onto_bounds(srf, uv)` (`closure.h`) : la surface au paramètre ramené
  sur ses bornes quand il les dépasse d'au plus 1e-5 de leur étendue (la règle de
  `face_volume_integral`, `onto_bounds`). `check` et `settle` l'emploient : l'écart
  d'une pcurve un peu hors du domaine est mesuré en 3D, comme tout écart, au lieu de
  lever une exception (`settle` agit avant que les pcurves soient ramenées).
- `contact_edge` mesure aussi l'écart sur la surface de la face : la tolérance de
  l'arête couvre ce que `check` mesure.

Les refus qui restent sur la pièce sont géométriques (une chaîne d'arêtes tangentes
prise en partie, des faces tangentes le long de l'arête).

### 2.4 Invariants

- Les pcurves des faces B-spline d'un solide rendu par un congé ou un chanfrein sont
  dans les bornes de leur surface (à l'arrondi près qu'une B-spline accepte).
- Non garanti ailleurs : une pcurve d'un autre opérateur peut dépasser le domaine de
  sa face (le maillage, `classify` et les intégrales le tolèrent à 1e-5 près ; `check`
  aussi désormais).

## 3. G16 : la raison d'un refus

### 3.1 Le défaut

`fillet`, `chamfer` et leurs variantes refusent presque tout avec `NonGeneric` et un
texte anglais. Une interface doit lire le texte pour dire à l'utilisateur quoi
changer.

### 3.2 Le correctif

`BuildError` reçoit deux champs, avec des valeurs par défaut (les agrégats existants
restent valides) :

```cpp
struct BuildError
{
    BuildErrc code;
    std::string message;
    std::vector<ShapeId> shapes;
    RefusalReason reason{};                                 // Unspecified hors des congés
    double limit{std::numeric_limits<double>::quiet_NaN()}; // la limite de la taille, connue
};
```

`RefusalReason` (`model.h`, avec `to_string`) :

| Raison | Ce que l'utilisateur peut changer | Refus (exemples) |
|---|---|---|
| `TooLarge` | réduire le rayon ou la distance | « a contact runs off its face », « a blend wider than a face », « no point of a face at the distance asked » |
| `Curvature` | réduire le rayon sous `limit` | « the radius exceeds a radius of curvature of a face » |
| `TangentFaces` | rien : arête entre faces tangentes | « the faces of the edge are tangent along it », sections presque plates |
| `PartialChain` | prendre toute la chaîne | « a chain of tangent edges blended in part » |
| `MixedConvexity` | séparer les arêtes | arête, onglet ou coin convexe et concave |
| `Overflow` | réduire, ou traiter l'arête voisine | débordement sur une arête vive, poche |
| `Vertex` | changer les arêtes prises à un sommet, le recul | onglets, coins, raccords de sommet, lunules |
| `Seam` | — | contact à travers une couture non traversable |
| `RadiusLaw` | changer la loi | loi non positive, sans enveloppe, pli |
| `Approximation` | relâcher les tolérances | courbe ou surface non approchée |
| `FaceNotOfEdge` | choisir une face de l'arête | « the face of d1 is not a face of the edge » (`InvalidId`) |
| `InvalidResult` | — | « the result is not a valid solid » |
| `Unsupported` | changer d'options | loi de rho avec un rayon variable, corde constante en G2 |
| `Geometry` | — | ligne des centres, section, arête à couper, exception de la géométrie |

Chaque refus de `fillet.h` et de `blend.h` pose sa raison à la source
(`fillet_detail::refuse<T>(code, reason, op, why, shapes, limit)`,
`blend_detail::refuse<T>(reason, op, why, shapes)`). Les refus des congés entre
faces (`face_blend.h`, `full_round.h`) restent `Unspecified` : hors de ce lot.

**Les textes ne changent pas** : ils restent le détail, et une interface qui les lit
encore n'est pas cassée.

### 3.3 La limite

`limit` n'est posée que là où gbs la connaît : `Curvature`, le rayon auquel le
décalage de la face (du côté de la boule, sur la boîte du support) se replie, 1 / k
pour k la plus grande courbure principale vers la boule
(`fillet_detail::curvature_limit`). Elle se lit sur la marge 1 − d k que trouve
`offset_margin` : |d| / (1 − marge), d le rayon demandé. Une seule recherche, la même
que celle qui refuse : le refus garde son coût. Exacte pour une primitive (le rayon
d'un cylindre, d'une sphère) ; pour une autre surface, aussi précise que la recherche
du point le plus courbe (échantillons puis motif, `offset_margin`). Test : la
demi-ellipse extrudée (a = 1, b = 0,3), limite b² / a à 1e-3 près, le refus en moins
de 8 recherches.

C'est une borne : un rayon plus petit peut encore être refusé pour une autre raison.
La courbure est celle de toute la boîte du support, la face prolongée au-delà de
l'arête : une limite très petite dit un point très courbe de ce support, qui n'est
pas forcément sous le contact. NaN quand elle n'est pas sous |d| (le refus d'une
primitive pour une autre cause : l'apex d'un cône) ni au-dessus de tol (aucun congé
ne la prendrait).

Pour `TooLarge`, la limite demanderait de résoudre la largeur de la face le long du
contact : non calculée (NaN).

### 3.4 Un refus rangé en `Geometry` : l'arc de plus de 300°

« an arc of more than 300 degrees between its contacts » (`arc_through`) n'est levé
que pour une section qui passe par un point imposé (`Station::via`), ce que pose
seulement `full_round.h` (le congé plein entre deux flancs qui se rapprochent loin de
la face du milieu). Un congé ou un chanfrein d'arête ne le rencontre jamais ; ce
n'est pas la tangence des faces de l'arête (les sections « nearly tangent », conique
et G2, à plus de 170°, sont `TangentFaces`). Il reste `Geometry`.

## 4. G17 : un chanfrein asymétrique à une face par arête

`chamfer(m, s, edges, f, d1, d2)` prend une face `f` que chaque arête borde. Les
arêtes sans face commune (celles d'un coin, de deux faces opposées) demandaient
plusieurs appels. Nouvelle surcharge, l'ancienne gardée :

```cpp
chamfer(m, s, std::span<const EdgeId> edges, std::span<const FaceId> faces, d1, d2, o);
chamfer(m, s, const std::vector<EdgeId> &, const std::vector<FaceId> &, d1, d2, o);
```

`faces[i]` porte d1 pour `edges[i]`. Refus : `InvalidOptions` sans une face par arête ;
`InvalidId` et `FaceNotOfEdge` pour une face qui ne borde pas son arête (le refus de
`fillet_detail::chamfer`). Les sommets, chaînes et onglets sont ceux de `chamfer`.

## 5. Tests

| Test | Ce qu'il vérifie |
|---|---|
| `tests_brep_sweep.revolutions` | couronnes planes (tour complet et demi-tour) ; un rectangle sur l'axe à 2 pi, pi, pi / 3 : un cylindre, deux plans (un cercle, ou un arc et deux rayons), volume à 1e-6 |
| `tests_brep_boolean.revolved_cylinders` | G15, la reproduction de l'issue #257 : cylindres de révolution R 25 à 30 l'un de l'autre, décalés en hauteur ; différence et intersection : une face par surface (5 et 4), volumes de la lunule |
| `tests_brep_fillet.revolved_cylinder_rims` | G18 : congé et chanfrein d'un cercle d'un cylindre de révolution (R 1 et R 5, la sonde de l'interface), puis des deux cercles ; volumes par Pappus |
| `tests_brep_fillet.bspline_faces_ending_on_a_boundary` | G19 : chaque arête de `occt_loft.stp`, congé et chanfrein à 1 % et 3 % : 14 construits chaque fois, valides, leurs pcurves dans les bornes des faces, maillés ; aucun refus « out of bounds » |
| `tests_brep_fillet.curvature_limit_on_a_bspline_face` | G16 : la limite `Curvature` sur une face B-spline (b² / a à 1e-3 près) et son coût (moins de 8 recherches de la marge) |
| `tests_brep_fillet.refusals`, `blend_refusals`, `fillet_chain` | G16 : `TooLarge`, `FaceNotOfEdge`, `Curvature` (limite 0,5 pour un cylindre R 0,5), `TangentFaces`, `PartialChain` ; `Unspecified` pour un rayon nul |
| `tests_brep_fillet.chamfer_edges` | G17 : deux arêtes sans face commune, d1 sur le dessus pour l'une, sur l'arrière pour l'autre : volume, points de contact ; refus |

## 6. Points à valider

1. G15/G18 : le plan pour un côté perpendiculaire à l'axe, plutôt qu'un congé qui
   traverse la couture d'une B-spline fermée (§ 1.2) ; ce dernier reste à faire pour
   d'autres profils.
2. L'ordre des faces avant `sew` (plans d'abord) et la tolérance de pcurve 1e-10 r2
   des cercles (§ 1.3), pour garder les volumes exacts.
3. G19 : ramener les pôles des pcurves sur les bornes des faces B-spline à la fin du
   congé, et tolérer 1e-5 de l'étendue dans `check` et `settle` (§ 2.3).
4. G16 : un champ `reason` dans `BuildError` plutôt que de nouveaux `BuildErrc`
   (les codes restent ceux des builders ; la raison ne sert qu'aux congés) ; les
   quatorze raisons et leur classement (§ 3.2).
5. G16 : la limite seulement pour `Curvature`, lue sur la marge de la recherche qui
   refuse (§ 3.3) ; l'arc de plus de 300° laissé en `Geometry` (§ 3.4).
6. G17 : `InvalidOptions` quand le nombre de faces n'est pas celui des arêtes.
