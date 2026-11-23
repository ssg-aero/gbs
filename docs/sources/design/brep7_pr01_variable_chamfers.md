# Palier 7 — PR 1 : chanfreins variables

| | |
|---|---|
| PR | branche `feat/brep-variable-chamfers` |
| Document de référence | [brep_palier7_blend_sections.md](brep_palier7_blend_sections.md), § 2 et § 3 |
| Fichiers | `gbs-brep/radius_law.h`, `gbs-brep/fillet.h`, `gbs-brep/blend.h`, `tests/tests_brep_variable_chamfers.cpp` (nouveau), `gbs-occt/tests/tests_brep2_occt.cpp`, `python/gbsBindBrepOps.cpp`, `python/tests/test_brep_ops.py`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier7_blend_sections.md` (nouveau), `News.md`, cette note |
| Statut | À valider (§ 9) ; revue de code : § 10 |

## 1. Objectif

La première PR du palier 7 apporte les chanfreins dont les distances varient le
long de l'arête :

- **deux lois** : d₁ sur une face donnée, d₂ sur l'autre ;
- **une loi** pour les deux faces ;
- **distance-angle** : d₁, en loi ou constante, et un angle constant depuis la
  face de d₁.

Une loi linéaire qui change de pente fait un **pli**. Toute la topologie des
paliers 4 et 6 s'applique : extrémités libres, chaînes, boucles, onglets, coins.

## 2. Sources

- Comportement de CATIA (« Chamfer », *Length1/Length2*, *Length/Angle*) et de
  SolidWorks (*Distance-distance*, *Angle-distance*).
- OCCT (`BRepFilletAPI_MakeChamfer::Add(d1, d2, E, F)`, `AddDA`) en boîte
  noire, pour les distances constantes.

## 3. API

```cpp
template <std::floating_point T> struct BlendPoint { EdgeId edge; T at{}; T value{}; };
enum class BlendVariation { Linear, Smooth };
template <std::floating_point T> struct BlendLaw { std::vector<BlendPoint<T>> points; BlendVariation variation;
    static BlendLaw linear(EdgeId, T v1, T v2); static BlendLaw constant(EdgeId, T v); };
// alias du palier 5 : RadiusPoint = BlendPoint, RadiusVariation = BlendVariation, RadiusLaw = BlendLaw

chamfer(m, s, edges, face, BlendLaw d1, BlendLaw d2, o);   // deux lois
chamfer(m, s, edges, BlendLaw d, o);                       // la même loi des deux côtés
chamfer_angle(m, s, edges, face, d1 (valeur ou loi), angle, o);  // angle en radians, ]0, π/2[
chamfer_geometry(m, s, e, face, BlendLaw d1, BlendLaw d2, o);
chamfer_angle_geometry(m, s, e, face, d1, angle, o);
```

Le champ `radius` de `RadiusPoint` devient `value`. Les appels positionnels du
palier 5 ne changent pas. En Python, `RadiusPoint.radius` reste un nom du champ,
et `RadiusLaw`, `RadiusPoint` et `RadiusVariation` sont les mêmes classes que
`BlendLaw`, `BlendPoint` et `BlendVariation`.

## 4. Démarche

### 4.1 Lois

Les lois sont celles du palier 5 (`law_chains`, `chain_radius`) : abscisse
curviligne le long de chaque chaîne tangente, variation linéaire ou lisse,
prolongement avec la pente au-delà des points extrêmes. La fonction de la loi
sur chaque arête est la même que pour un rayon (`EdgeRadius`). Le découpage en
chaînes est désormais partagé par les congés et les chanfreins
(`law_on_edges`), et les messages de refus disent « distance » au lieu de
« rayon ».

### 4.2 Contacts

`contact_pcurve` prend une distance fonction du paramètre de l'arête, au lieu
d'une constante. Chaque pied reste le point de la face, dans le plan normal à
l'arête, à la distance d(t) du point de l'arête. Le reste est inchangé :
Gauss-Newton de station en station, interpolation, stations doublées jusqu'à la
précision des contacts. Aux points anguleux des lois, le contact est calculé
**morceau par morceau**, puis les morceaux sont joints en une B-spline (C0 à
ces points).

### 4.3 Distance-angle

Le contact 1 est à la distance d₁(t). On trace, dans le plan normal à l'arête,
la droite qui part de ce contact sous l'angle θ avec le plan tangent de la face 1
en ce point : vers l'intérieur du solide à une arête convexe, vers l'extérieur à
une arête concave. Le contact 2 est le point où elle rencontre l'autre face
(Gauss-Newton à trois inconnues, repris depuis le point trouvé le plus proche).
d₂(t) est alors la distance de ce point à l'arête, et le contact 2 se construit
comme une distance en loi. Entre deux plans à angle droit, d₂ = d₁ tan θ ; à
120°, la loi des sinus donne d₂ = d₁ sin θ / sin(120° + θ). OCCT donne le même
volume à 10⁻⁹ (§ 7).

La droite part de la tangente de la face au contact, et non de la corde vers
l'arête (revue de code, § 10) : sur une face courbe en travers de l'arête, la
corde faisait 36° au lieu de 30° (d₁ = 0,3 sur un cylindre de rayon 1).

Le prolongement et la portée des supports demandent d₂ avant de le calculer.
On l'estime par la loi des sinus entre les plans tangents des faces, au
dièdre le plus ouvert le long de l'arête, majoré de moitié pour leur courbure.

### 4.4 Surface

- **Un morceau** : primitive si les contacts le permettent (plan, cylindre,
  cône, comme au palier 4), B-spline réglée sinon.
- **Surface réglée exacte** (nouveau) : quand les deux contacts sont des
  B-splines sur des plans, la surface réglée est construite exactement, sur
  leurs pôles, après unification des degrés et des nœuds. Une pcourbe B-spline
  sur un plan est en effet une B-spline de l'espace, le plan étant affine en ses
  paramètres. C'est le cas de tous les chanfreins variables entre deux plans,
  lois lisses comprises. Avant, l'interpolation de la surface réglée laissait un
  écart de 3·10⁻⁶ pour une loi lisse ; il est maintenant nul. Une loi linéaire
  le long d'une arête droite donne ainsi un paraboloïde hyperbolique exact. Le
  chanfrein est alors marqué `exact`, comme une primitive.
- **En morceaux**, aux points anguleux : un morceau réglé par intervalle, de
  paramètres (t, en travers), exact entre les morceaux des contacts sur des
  plans. Un **pli** suit la règle à chaque point anguleux :
  segment de droite, pcourbes iso-t sur les deux morceaux. C'est la topologie
  des plis du palier 5 (`BlendPiece`, `BlendCrease`), que `blend.h` sait déjà
  traiter.

## 5. Refus

- `InvalidBounds` :
  - un angle hors de ]0, π/2[ ;
  - une distance qui ne dépasse pas la tolérance sur le chanfrein prolongé (une
    loi linéaire raide prolongée avec sa pente) ; le message nomme
    `FilletOptions::extension` ;
  - un point de loi sur une arête non chanfreinée ;
  - deux valeurs au même point.
- `InvalidId` : la face de d₁ n'est pas une face de l'arête.
- `NonGeneric` : la droite tracée sous l'angle ne rencontre pas l'autre face, en
  plus des refus du chanfrein du palier 4.

## 6. Tests

`tests_brep_variable_chamfers` (8 cas) :

| Cas | Attendu | Mesuré |
|---|---|---|
| arête de boîte, deux lois : linéaire et constante ; lisse et lisse ; linéaire avec un point au milieu | contacts à d(t) du bord ; surface exacte (`exact`) ; volume 3 − ∫ d₁ d₂ / 2 ; un morceau, ou deux et un pli | contacts à 10⁻⁶ ; écart 0 ; volumes à 10⁻⁹ ; 7 ou 8 faces |
| une loi des deux côtés | volume 3 − ∫ d² / 2 | 10⁻⁹ |
| distance-angle à 90° et à 120° | plan exact à 90° ; d₂ par la loi des sinus à 120° ; volumes exacts | 10⁻⁹ |
| distance-angle sur un demi-cylindre, d₁ sur le cylindre (courbe en travers de l'arête) | 30° depuis le plan tangent au contact, à 9 stations ; volume : triangle et segment de cercle | 10⁻⁶ |
| bord d'un cylindre, loi lisse sur la boucle | Pappus tranche par tranche : ∫ d²/2 (1 − d/3) dθ | 10⁻⁶ |
| bord d'une dalle en stade, loi lisse sur la chaîne fermée (deux droites, deux demi-cercles) | continuité d'une arête à l'autre ; volume par tranches, avec la courbure | 10⁻⁶, 10 faces |
| trois chanfreins variables au coin d'une boîte | solide valide, 9 faces, un sommet commun aux trois chanfreins ; matière enlevée entre le plus grand chanfrein seul et la somme des trois | valide |
| refus | angle à 90° et à 0 ; distance négative ; loi qui s'annule dans le prolongement (acceptée avec `extension` = 0,5) ; point sur une autre arête ; deux valeurs au même point ; mauvaise face ; angle de 89° qui sort de la face | codes et messages du § 5, modèle inchangé |

Python : volumes de deux lois, d'une loi des deux côtés et de distance-angle ;
morceaux et plis de la géométrie ; alias, et les mots-clés du palier 5
(`radius`, `r1`, `r2`). Galerie, trois cas :
`chamfer_variable_box` (pli), `chamfer_angle_prism` (120°) et
`chamfer_variable_slab` (chaîne fermée).

## 7. Comparaison à OCCT

`BRepFilletAPI_MakeChamfer` n'a pas de distance variable. Les trois cas
constants sont comparés : deux distances 0,1 / 0,3 sur une arête de boîte,
distance-angle 0,3 / 30° sur la même arête, et 0,5 / 40° entre faces à 120°.
Les deux solides sont valides, ont les mêmes faces, et leurs volumes concordent
à 10⁻⁹. L'angle d'OCCT (`AddDA`) est donc mesuré comme le nôtre, depuis la face
de la distance.

## 8. Limites

- **Lois raides.** Une loi linéaire est prolongée avec sa pente au-delà des
  sommets, sur une longueur égale à `extension` × la plus grande distance. Une
  loi qui descend vite vers un sommet peut ainsi devenir négative dans ce
  prolongement : elle est refusée, comme au palier 5. Réduire `extension` est
  alors la parade.
- **Exactitude.** La surface réglée n'est exacte que si les deux contacts sont
  des B-splines sur des plans. Sur des faces courbes, elle reste approchée sous
  `approx_tol`, comme au palier 4.
- **Distance-angle sur des faces courbes.** L'angle est mesuré dans le plan
  normal à l'arête, depuis le plan tangent de la face au contact 1. C'est, à
  notre connaissance, la définition de CATIA et d'OCCT ; le test ne la confirme
  pour OCCT qu'entre deux plans, où plan tangent et face se confondent.

## 9. Points à valider

1. **`BlendLaw`** remplace `RadiusLaw`, avec le champ `value`. `RadiusLaw`,
   `RadiusPoint` et `RadiusVariation` restent des alias, et `radius` reste un
   nom du champ en Python.
2. **Contacts en morceaux** aux points anguleux, joints C0, avec un pli le long
   de la règle.
3. **Distance-angle par intersection** de la droite avec l'autre face, à chaque
   station, l'angle pris depuis le plan tangent de la face 1 au contact. d₂
   devient une loi implicite ; l'angle est constant.
4. **Surface réglée exacte** sur les pôles quand les contacts sont des B-splines
   sur des plans (nouveau, au-delà du document de palier).
5. **Refus des lois négatives** dans le prolongement, comme au palier 5.

## 10. Revue de code

Revue du 23 novembre 2026 : 8 points, tous corrigés.

1. **Angle sur une face courbe** (important). La droite partait de la corde du
   contact 1 vers l'arête, et non de la tangente de la face. Sur un
   demi-cylindre de rayon 1, avec d₁ sur le cylindre, 30° demandés donnaient
   32,1°, 36,3° et 40,5° pour d₁ = 0,1, 0,3 et 0,5. Elle part maintenant du plan
   tangent : 30° à 10⁻⁶ près (nouveau test, § 6).
2. **Mots-clés Python du palier 5.** `RadiusPoint(radius=...)` et
   `RadiusLaw.linear(r1=..., r2=...)` sont de nouveau acceptés, à côté de
   `value`, `v1` et `v2`.
3. **Surface exacte marquée `exact`.** La surface réglée sur les pôles, d'écart
   nul, n'était pas marquée exacte : les onglets, les jonctions et les coins
   élargissaient leurs seuils de `approx_tol`. Elle l'est maintenant, et les
   morceaux entre contacts sur des plans sont aussi construits exactement.
4. **Message des lois raides.** Il nomme `FilletOptions::extension`, la parade
   du § 8.
5. **d₂ estimé pour la distance-angle.** L'estimation 2 d₁ max(1, tan θ)
   sous-estimait d₂ aux arêtes obtuses (0,57 pour 1,88 à 120° et 55°). Elle
   suit la loi des sinus entre les plans tangents (§ 4.3).
6. **Tests.** Le coin de trois chanfreins compte ses faces, vérifie le sommet
   commun et encadre la matière enlevée ; les refus vérifient leurs messages ;
   un cas de distance-angle sur une face courbe.
7. **Doc de `chamfer_geometry`** à distances constantes, de nouveau au-dessus de
   sa surcharge.
8. **Messages du palier 5.** La note `brep5_pr01_variable_geometry.md` cite le
   nouveau libellé « two values of the radius ».
