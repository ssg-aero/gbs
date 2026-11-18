# Palier 4 — PR 5 : congés et chanfreins, géométrie

| | |
|---|---|
| PR | branche `feat/brep-fillet-geometry` |
| Document de référence | [brep_palier4.md](brep_palier4.md), § 4.3 et § 4.4 |
| Fichiers | `gbs-brep/fillet.h` (nouveau), `gbs-brep/brep`, `python/gbsBindBrepOps.cpp`, `tests/tests_brep_fillet.cpp` (nouveau), `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/examples/gbs_examples.py`, `docs/sources/examples/gallery.py`, `python/tests/test_brep_ops.py`, `python/tests/test_examples.py`, `docs/sources/design/brep_palier4.md`, `docs/sources/design/brep4_pr04_offset_solid.md`, `News.md`, cette note |
| Statut | À valider (§ 9) |

## 1. Objectif

La géométrie d'un congé de rayon constant et d'un chanfrein le long d'une arête
d'un solide :

- la ligne des centres d'un congé ;
- les deux courbes de contact sur les faces de l'arête, avec leurs pcurves ;
- la surface du congé ou du chanfrein.

Le modèle n'est pas modifié. La topologie vient avec la PR 6 : faces découpées
le long des contacts, extrémités coupées par les faces voisines, couture. Elle
consommera ce résultat tel quel.

## 2. Sources

- Congé à boule roulante, ligne des centres comme intersection des surfaces
  décalées : J. R. Rossignac, A. A. G. Requicha, *Constant-radius blending in
  solid modelling*, Computers in Mechanical Engineering, 1984 ; B. K. Choi,
  S. Y. Ju, *Constant-radius blending in surface modelling*, CAD 21(4), 1989.
- Surface du congé, enveloppe de la boule (surface canal), formes exactes et
  approchées : J. Vida, R. R. Martin, T. Várady, *A survey of blending methods
  that use parametric surfaces*, CAD 26(5), 1994.
- Arc de cercle rationnel de degré 2 (un ou deux segments), *skinning* :
  L. Piegl, W. Tiller, *The NURBS Book*, 2ᵉ éd., § 7.5 et § 10.3.
- Chanfrein à deux distances mesurées depuis l'arête : comportement de CATIA
  (*Length1/Length2*) et de SolidWorks (*Distance-Distance*).

OpenCascade ne sert que de comparaison en boîte noire dans les tests
(`BRepFilletAPI_MakeFillet`, `BRepFilletAPI_MakeChamfer`).

## 3. API

```cpp
template <std::floating_point T> struct FilletOptions {
    T tol;          // tolérance du résultat (1e-6)
    T approx_tol;   // surfaces et contacts approchés : écart au raccord exact (1e-5)
    T angular_tol;  // congés approchés : sinus de l'angle de leur normale à l'exacte (1e-5)
    IntersectOptions<T> intersect;
    T extension;    // longueur au-delà des sommets d'une arête ouverte, en rayon ou en plus grande distance (1)
};

auto fillet_geometry(const Model<T> &, SolidId, EdgeId, T r, const FilletOptions<T> & = {}) -> BuildResult<EdgeBlend<T>>;
auto chamfer_geometry(const Model<T> &, SolidId, EdgeId, T d1, T d2, const FilletOptions<T> & = {}) -> BuildResult<EdgeBlend<T>>;
```

`EdgeBlend` porte :

- `surface` : cylindre, tore, plan, cône, ou B-spline (rationnelle pour un
  congé), avec `exact` et `deviation` ;
- `spine` : la ligne des centres d'un congé, l'arête (prolongée) d'un
  chanfrein ;
- `contacts[2]` (`BlendContact`) : sur la face de la première utilisation de
  l'arête, puis sur l'autre. Chacun a sa face, la surface de la face
  prolongée (`support`), la pcurve sur elle, la courbe 3D, la pcurve sur la
  surface du raccord et son écart à cette surface ;
- `convex`, `closed` ;
- `[a, b]`, le paramètre commun de toutes ces courbes, prolongement compris ;
  `t1` et `t2`, les paramètres aux sommets v1 et v2 de l'arête (`t1 > t2`
  quand la ligne des centres va dans l'autre sens).

## 4. Congé (`fillet_geometry`)

### 4.1 Côté de la boule

En 9 points de l'arête :

- n₁ et n₂ sont les normales sortantes des deux faces ;
- b₁ et b₂ sont les directions des faces qui s'éloignent de l'arête :
  la gauche du co-edge dans (u, v), le fil extérieur tournant dans le sens
  direct et les trous dans le sens indirect, portée en 3D puis rendue normale à
  l'arête.

L'arête est convexe si b₁·n₂ < 0. La boule est alors dans la matière, du côté
−n des deux faces. Elle est dans le vide (+n) le long d'une arête concave.

Des faces tangentes le long de l'arête (|b₁·n₂| ≤ 10⁻⁶) sont refusées, comme
une arête convexe par endroits et concave ailleurs.

L'angle du dièdre où va la boule, α = acos(b₁·b₂), borne l'éloignement des
contacts : r / tan(α/2) pour deux plans.

### 4.2 Ligne des centres

Les surfaces des deux faces sont prolongées (`extended_surface` de la PR 4)
puis restreintes à la bande de leurs paramètres proche de l'arête : la boîte de
sa pcurve, élargie de r / tan(α/2) + r + l'extension, divisée par la vitesse.
On les décale ensuite de r vers la boule (`offset_surface` : exact pour les
primitives, refusé là où le décalé se replie).

Leur intersection (`intersect` : formes fermées entre primitives, cheminement
sinon) donne la ligne des centres. On garde la branche la plus proche du
milieu de l'arête, entre les points les plus proches de ses sommets, prolongée
de `extension`·r au-delà. Une arête fermée prend toute la boucle depuis son
sommet. Sur une boucle, une arête ouverte prend l'arc dont le milieu est le
plus proche du sien.

### 4.3 Contacts

Le contact sur une face est le pied du centre, c − r n. Sa pcurve est celle
de la ligne des centres sur la surface décalée (mêmes paramètres, au décalage
de v d'un cône près). La courbe 3D est cette pcurve sur la surface prolongée
de la face : exacte sur la face.

Le contact doit rester sur sa face entre les sommets : il est classé
(`classify`) en 15 points intérieurs à l'arête. Les bouts sont laissés à la
PR 6, où le contact rejoint les faces d'extrémité.

### 4.4 Surface

La surface est l'enveloppe de la boule, une surface canal. Elle ne dépend que
de la ligne des centres, quelles que soient les faces :

| Ligne des centres | Surface |
|---|---|
| droite | cylindre de rayon r, d'axe la droite ; en u, l'arc des contacts |
| cercle | tore de rayon mineur r, de même repère que le cercle (u est le paramètre du cercle) ; en v, l'arc des contacts |
| autre (ellipse, B-spline cheminée) | B-spline rationnelle approchée |

La règle couvre la table du plan (plan-plan, plan-cylindre ou plan-cône
d'axe normal au plan, cylindres coaxiaux) et les cas qu'elle oubliait : plan et
sphère, cylindres parallèles.

**Stations exactes.** Chaque station de la boule est reprise par Gauss-Newton
sur les équations de la boule roulante (Choi et Ju), en partant de la ligne des
centres cheminée et de ses pcurves sur les surfaces décalées. Les inconnues
sont (U₁, V₁, U₂, V₂), et les équations O₁(U₁, V₁) = O₂(U₂, V₂) (le centre
commun aux deux décalés) et (O₁ − c)·c′ = 0 (le centre dans le plan normal à la
ligne cheminée en s). La ligne cheminée n'est juste qu'à sa tolérance divisée
par le sinus de l'angle des décalés. Là où les faces deviennent presque
tangentes, ses contacts ne sont plus à la même distance r d'elle et les arcs
construits dessus ne sont plus des cercles. Sur un té de cylindres de rayons
0,99 et 1 (r = 0,15), l'écart monte à 3·10⁻⁷ près des points presque tangents.
L'écart de normale de la surface butait alors à 10⁻⁵ et le congé était refusé.
Après Newton, il est de 5·10⁻¹⁶ partout.

La station cheminée est gardée là où Gauss-Newton ne converge pas ou s'éloigne
de plus de r / 10 (une autre racine). C'est le cas au point où les faces sont
tangentes, où le système est singulier. La ligne des centres et les pcurves des
contacts sont réinterpolées (cubique) sur les stations, avec leurs paramètres
rendus continus à travers la couture d'une surface périodique. Elles sont
contrôlées à `approx_tol` entre les stations.

**Approximation.** À n + 1 stations de la ligne des centres, on prend l'arc
de cercle exact entre les contacts : rationnel de degré 2, en un segment, ou
en deux au-delà de 120° pour garder des poids loin de 0. Les pôles homogènes
de ces arcs et leurs dérivées le long de la ligne (différences centrées)
sont interpolés par Hermite cubique (C1). La surface a donc, aux stations, la
position et le plan tangent exacts des arcs. On double n à partir de 8
jusqu'à ce que la surface soit à `approx_tol` des arcs exacts et que sa
normale soit à `angular_tol` de la leur, en trois points entre deux stations.
La normale n'est pas contrôlée là où l'arc se réduit à un point : un congé qui
s'annule là où ses faces deviennent tangentes.

Mesures sur le congé des cylindres croisés (r = 0,15, rayons 0,4 et 1), aux mêmes stations, avant la reprise par Newton :

| n | positions seules : position, normale | Hermite : position, normale |
|---|---|---|
| 32 | 5,1·10⁻⁵, 7,9·10⁻⁴ | 9,9·10⁻⁶, 1,4·10⁻⁴ |
| 64 | 3,6·10⁻⁶, 1,1·10⁻⁴ | 6,5·10⁻⁷, 1,8·10⁻⁵ |
| 128 | 3,6·10⁻⁷, 1,4·10⁻⁵ | 6,2·10⁻⁸, 3,6·10⁻⁶ |
| 256 | 5,2·10⁻⁸, 2,2·10⁻⁶ | 1,1·10⁻⁸, 1,3·10⁻⁶ |

Les deux schémas convergent au même ordre ; Hermite divise les deux écarts par 5
environ, soit un doublement de stations de moins. C'est la tolérance angulaire
qui garantit la tangence : avec la seule tolérance de position, l'arrêt se fait
à n = 64 et la normale reste à 10⁻⁴ (3,4·10⁻⁴ le long des contacts, mesuré par
les tests). Au-delà de n = 256, la normale plafonne vers 10⁻⁶ (précision de la
ligne des centres cheminée, dérivées par différences centrées), d'où
`angular_tol` = 10⁻⁵ (0,0006°). La surface retenue (n = 128) est à 3·10⁻⁸ au
plus de l'enveloppe exacte de la boule.

Sur le tore et le cylindre, les pcurves des contacts sont ajustées
(`primitive_pcurve`). Sur la B-spline, ce sont ses côtés u = 0 et u = 1, avec
leur écart à la courbe de contact.

## 5. Chanfrein (`chamfer_geometry`)

Les distances d₁ et d₂ sont mesurées depuis l'arête dans le plan qui lui est
normal. Le contact sur une face est le point de la face à la distance d du
point e(t) de l'arête, dans ce plan, du côté de la face : on résout par
Gauss-Newton (|S − e|² − d²) / 2d = 0 et (S − e)·e′/|e′| = 0.

On part du milieu de l'arête, de sa pcurve déplacée de d·b par la métrique,
puis on avance de station en station. Le déplacement de l'arête est ramené sur
le plan tangent au pied précédent (prédicteur du premier ordre), et le pas est
divisé par deux là où Gauss-Newton ne converge pas. Les pieds sont interpolés
(cubique) en une pcurve, avec n doublé jusqu'à `approx_tol` aux milieux.

L'arête est prolongée de `extension`·max(d₁, d₂) : une droite ou une B-spline
continuée par le polynôme de son dernier intervalle (`extrapolated`).

Surface réglée entre les deux contacts :

| Contacts | Surface |
|---|---|
| deux droites dans un plan | plan |
| deux cercles coaxiaux | plan (même hauteur), cylindre (même rayon), cône |
| autres | B-spline réglée, les deux contacts interpolés aux mêmes stations (degré 1 en travers) |

Les cas exacts sont reconnus sur les contacts eux-mêmes, en 33 points, à `tol`
près. Ils couvrent le chanfrein entre deux plans, et celui du bord d'un
cylindre ou d'un cône d'axe normal au plan.

**Écart au plan.** Le plan (§ 4.4) prenait le contact sur F₁ à
l'intersection de F₁ et de S₂ décalée de d₂. Cela mesure la distance à
l'autre face, pas à l'arête : les deux ne s'accordent qu'à angle droit. La
définition retenue est celle de CATIA et SolidWorks. OCCT l'applique aussi :
les deux chanfreins coïncident à 2,6·10⁻⁶ près le long d'une face B-spline
(§ 8).

## 6. Refus

- `NonGeneric` :
  - un rayon au-delà d'un rayon de courbure d'une face du côté de la boule
    (le décalé se replie) ;
  - un contact hors de sa face (un rayon ou une distance au-delà de la
    largeur de la face) ;
  - des faces tangentes le long de l'arête ;
  - une arête convexe par endroits et concave ailleurs ;
  - aucune ligne des centres près de l'arête ;
  - aucun point d'une face à la distance demandée ;
  - une approximation qui n'atteint pas ses tolérances.
- `InvalidId` : une arête qui n'est pas du solide, une arête utilisée par une
  seule face ou par plus de deux, une couture (une face l'utilise deux fois).
- `InvalidBounds` : un rayon ou une distance pas plus grand que `tol`.

## 7. Python, galerie

- `gbs.brep.fillet_geometry`, `chamfer_geometry`, `FilletOptions`,
  `EdgeBlend`, `BlendContact`.
- Galerie, section « Fillet and chamfer geometry », en 4 cas : un bloc en L
  (congé convexe et concave, chanfrein), un cylindre (tore, cône), deux
  cylindres croisés (B-spline rationnelle) et une bosse B-spline (congé et
  chanfrein approchés). Les surfaces sont dessinées prolongées au-delà des
  sommets, les contacts en rouge, la ligne des centres en noir. Images
  vérifiées.

## 8. Tests

- **`tests_brep_fillet`** (0,15 s). Chaque raccord est contrôlé en 33 × 9
  points :
  - la surface passe par les contacts ;
  - la normale de la surface est celle de la face le long d'un contact ;
  - un contact est à r de la ligne des centres, le long de la normale, du
    bon côté ;
  - la surface est à r de la ligne des centres ;
  - un chanfrein a ses contacts à d₁ et d₂ dans le plan normal.

  Les cas :
  - arête d'une boîte : cylindre, ligne des centres en (2 − r, ·, 1 − r),
    prolongée de r ;
  - arête concave d'un bloc en L : cylindre extérieur, chanfrein plan ;
  - bords d'un cylindre : tore de rayon R − r (convexe), tore de rayon R + r
    au pied d'un cylindre fusionné sur une boîte (concave), chanfrein
    conique ;
  - boule coupée par un plan : tore de rayon √((R − r)² − (h − r)²) ;
  - cylindres croisés : ligne des centres cheminée, B-spline rationnelle sous
    `approx_tol` et `angular_tol`, chanfrein réglé ;
  - prisme d'une bosse B-spline (plan et B-spline) : congé et chanfrein
    approchés ;
  - faces presque tangentes : tés de cylindres de rayons 0,99 et 0,999 sur
    un rayon 1, r = 0,05 et 0,15, sous `approx_tol` et `angular_tol` ; té de
    même rayon, où le congé s'annule (contacts confondus à 10⁻⁶ près) au
    sommet où les faces sont tangentes ;
  - coins aigus (prisme d'un triangle de 20°, 5° et 1°) : cylindre, contacts à
    r / tan(a/2) de l'arête, refusé quand ils sortent de la face ;
  - refus : rayon et distances nuls ou négatifs, rayon et chanfrein plus
    larges que la face, arête d'un autre solide, rayon au-delà de celui d'un
    cylindre, couture.
- **`tests_brep2_occt.fillet_geometry`**. Pour chaque cas, on mesure la plus
  grande distance de 17 × 9 points de la surface de gbs, sur l'arête, aux
  faces qu'OCCT génère depuis l'arête, sur les mêmes solides (ceux de gbs, lus
  par OCCT) :

  | Cas | Écart à OCCT |
  |---|---|
  | congés convexe et concave, chanfreins (bloc en L) | 3·10⁻¹⁶ |
  | congé (tore) et chanfrein (cône) d'un cylindre | 8·10⁻¹⁶ |
  | congé et chanfrein d'une bosse B-spline | 3,0·10⁻⁶ et 2,6·10⁻⁶ |
  | chanfrein de cylindres croisés | 2,9·10⁻⁶ |
  | congé de cylindres croisés | 1,1·10⁻⁴ |

  Sur les cylindres croisés, c'est OCCT qui s'écarte de la boule roulante.
  Les points de sa face sont à 8·10⁻⁵ de la distance r de la ligne des centres
  de gbs, alors que les contacts de gbs y sont à 10⁻⁶ près, le long des
  normales des faces. Sa tolérance d'approximation vaut 10⁻⁴ par défaut. Ce
  cas est donc comparé à 2·10⁻⁴, les autres cas approchés à 10⁻⁵ et les cas
  exacts à 10⁻⁷.
- **Python** : `test_brep_ops.py::test_fillet_and_chamfer_geometry` et
  `test_examples.py::test_blends`.

## 9. Points à valider

1. **Géométrie publique sans topologie** : `fillet_geometry` et
   `chamfer_geometry` ne modifient pas le modèle ; `fillet` et `chamfer`
   (PR 6) s'appuieront dessus.
2. **Côté de la boule** lu sur b₁·n₂ en 9 points ; refus des faces tangentes
   (|b₁·n₂| ≤ 10⁻⁶) et des arêtes mixtes.
3. **Surface exacte décidée par la ligne des centres** : une droite donne un
   cylindre, un cercle donne un tore, quelles que soient les faces. Cela
   élargit la table du plan (plan-sphère, cylindres parallèles). Une ellipse
   (plan oblique sur un cylindre) donne une B-spline.
4. **Congé approché par Hermite C1** des arcs exacts, avec deux tolérances :
   `approx_tol` (10⁻⁵) en position et `angular_tol` (10⁻⁵) en normale. Le plan
   ne prévoyait qu'une tolérance de position, qui laisse la tangence à 10⁻⁴ ;
   Hermite divise les écarts par 5 à nombre de stations égal.
5. **Prolongement de `extension`** (1 r ou 1 max(d₁, d₂)) au-delà des sommets
   d'une arête ouverte, pour la coupe par les faces d'extrémité de la PR 6.
6. **Chanfrein mesuré depuis l'arête dans son plan normal**, comme CATIA,
   SolidWorks et OCCT, et non par l'intersection avec le décalé de l'autre
   face, comme l'écrivait le plan (§ 5).
7. **Contact sur sa face** contrôlé en 15 points intérieurs à l'arête ; les
   bouts sont laissés à la PR 6.
8. **Comparaison à OCCT** à 10⁻⁷ pour les cas exacts, à 10⁻⁵ pour les cas
   approchés, et à 2·10⁻⁴ pour le congé des cylindres croisés, où c'est
   l'approximation d'OCCT qui s'écarte.
9. **Stations reprises par Newton** sur les équations de la boule roulante,
   la ligne des centres cheminée ne servant que d'amorce. Les faces presque
   tangentes passent ainsi.
10. **Congé qui s'annule** là où ses faces deviennent tangentes à un sommet
    (té de même rayon) : la géométrie le donne, avec des contacts confondus.
    La PR 6 devra borner une face de congé qui finit en pointe.
