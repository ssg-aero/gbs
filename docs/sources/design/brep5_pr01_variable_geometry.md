# Palier 5 — PR 1 : congés à rayon variable, loi et géométrie

| | |
|---|---|
| PR | branche `feat/brep-variable-fillet` |
| Document de référence | [brep_palier5_variable_fillets.md](brep_palier5_variable_fillets.md), § 2 et § 3 |
| Fichiers | `gbs-brep/radius_law.h` (nouveau), `gbs-brep/fillet.h`, `tests/tests_brep_variable_fillet.cpp` (nouveau), `docs/sources/design/brep_palier5_variable_fillets.md` (nouveau), `docs/sources/design/brep4_pr07_fillet_corners.md`, `News.md`, cette note |
| Statut | À valider (§ 9) |

## 1. Objectif

La géométrie d'un congé dont le rayon suit une loi le long d'une arête :
`fillet_geometry(m, s, e, RadiusLaw, options)`. La topologie (`fillet` avec une
loi, faces des fuseaux, chaînes, coins), Python, la galerie et la comparaison à
OCCT sont la PR 2.

Le document du palier part avec cette PR (pas de PR de documentation seule). Le
statut de la note de la PR 7 du palier 4 y passe à « Validé (§ 9 ; § 10) ».

## 2. Sources

- Fritsch et Carlson, « Monotone piecewise cubic interpolation » (*SIAM Journal
  on Numerical Analysis* 17(2), 1980) ; moyenne harmonique pondérée de Fritsch
  et Butland pour les pentes.
- Choi et Ju, « Constant-radius blending in surface modelling »
  (*Computer-Aided Design* 21(4), 1989) : stations de la boule par
  Gauss-Newton.
- Do Carmo, *Differential Geometry of Curves and Surfaces*, § 3-5 : enveloppe
  d'une famille de sphères, cercle caractéristique.
- CATIA (« Variable Radius Fillet », variation linéaire ou cubique) et
  SolidWorks (« Variable Size Fillet », transition droite ou lisse) pour
  l'usage.

## 3. API

```cpp
template <std::floating_point T> struct RadiusPoint { EdgeId edge; T at{}; T radius{}; };
enum class RadiusVariation { Linear, Smooth };
template <std::floating_point T> struct RadiusLaw
{
    std::vector<RadiusPoint<T>> points;
    RadiusVariation variation{RadiusVariation::Smooth};
    static auto linear(EdgeId e, T r1, T r2) -> RadiusLaw;
};
auto fillet_geometry(const Model<T> &, SolidId, EdgeId, const RadiusLaw<T> &, const FilletOptions<T> & = {}) -> BuildResult<EdgeBlend<T>>;
```

`EdgeBlend` gagne `pieces` (`BlendPiece` : surface, intervalle, contacts sur la
surface, écart) et `lunes` (`BlendLune` : paramètre, sphère, deux arcs, leurs
pcourbes sur la sphère). Un congé sans point anguleux garde `surface` et ne les
remplit pas.

## 4. Démarche

### 4.1 Loi

`radius_detail::radius_function` (`radius_law.h`) : une `BSCfunction` de gbs,
fonction de l'abscisse curviligne.

- **Linéaire** : degré 1 par les points.
- **Lisse** : Hermite cubique, pentes de Fritsch et Carlson (nulles là où les
  pentes des intervalles voisins changent de signe, moyenne harmonique
  pondérée sinon), écrite en B-spline de degré 3 à nœuds doubles (C1).
  Monotone entre deux points. Pentes nulles aux bouts d'une loi ouverte ; une
  loi fermée prend sa pente à travers la fermeture.

Sur l'arête (`edge_radius`), les points pris sont ceux de l'arête : fraction
`at` de la longueur depuis le premier sommet. L'abscisse s(t) et son inverse
t(s) sont interpolées par les longueurs de 64 intervalles. Au-delà de l'arête,
l'abscisse continue à la vitesse de l'arête au sommet, et la loi avec sa pente
(§ 2.2 du palier : constante pour une loi lisse, droite pour une linéaire).
Une arête fermée a une loi fermée, son sommet porte un seul rayon.

`gbs::abs_curv` donne t(s) et non s(t) : sur les droites et le cercle unité les
deux coïncident, ce qui a d'abord masqué l'inversion.

### 4.2 Stations

À la station t, la boule de rayon r(t) dans le plan normal à l'arête en e(t) :
Gauss-Newton sur (u₁, v₁, u₂, v₂), les décalés `OffsetSurface` des supports à
±r(t) devant coïncider, le centre dans le plan. Le départ est la station déjà
résolue la plus proche ; à défaut, l'arête elle-même (boule de rayon nul), le
rayon grandi en huit pas. Une solution plus loin de l'arête que la portée des
contacts (une autre racine) est rejetée. Les supports sont ceux du palier 4
(faces prolongées, vue sans couture des B-splines fermées), contrôlés contre le
repli au plus grand rayon de la loi.

### 4.3 Enveloppe

Le cercle caractéristique est dans le plan normal à c′, à −r r′/‖c′‖ du
centre. Les contacts y sont, et on en tire son centre sans r′ : m = c + l t̂,
t̂ = c′/‖c′‖ par différences des stations, l la moyenne de (pᵢ − c) · t̂ sur les
deux contacts. Si |l| atteint r, une boule contient la suivante : refus. La
section est l'arc de ce cercle entre les contacts (`Station::m`, que `skinned`
prend pour centre de l'arc). La normale de l'enveloppe reste celle de la
sphère, contrôlée par `skinned` comme avant.

### 4.4 Surface

`skinned`, inchangé pour l'essentiel : arcs exacts aux stations, Hermite,
raffinement local. Une station aux nœuds des faces moins que C3
(`knot_breaks`, factorisé avec le rayon constant) et aux nœuds d'une loi
lisse. Un intervalle dont un point de contrôle sort du support (pcourbe
grossière au début du raffinement) est maintenant coupé en deux au lieu de
lever une exception.

### 4.5 Morceaux et fuseaux

Aux points anguleux d'une loi linéaire (§ 3.5 du palier), le congé est coupé :
chaque morceau est approché seul, les différences de ses stations prises dans
son intervalle, donc d'un seul côté aux bouts. Les morceaux n'ont pas le même
arc au point : ce sont les deux cercles caractéristiques. La ligne des centres
et les contacts des morceaux sont joints en une B-spline C0 (`joined`), les
pcourbes ramenées sur les périodes de la face d'un morceau au suivant. Le
fuseau a pour sphère la boule au point. Son repère met x vers le milieu du
fuseau et z en travers des contacts, ce qui tient la couture et les pôles loin
de lui. Ses arcs sont les bords des morceaux (`CurveOnSurface` sur leur iso au
point), leurs pcourbes sur la sphère interpolées sur 33 points. Sur une arête
fermée dont la loi linéaire se referme avec un angle, le fuseau est au sommet,
entre le dernier morceau et le premier.

## 5. Refus

| Situation | Code | Message |
|---|---|---|
| aucun point de la loi sur l'arête | `InvalidBounds` | `no radius on the edge` |
| deux rayons au même point, deux rayons au sommet d'une arête fermée | `InvalidBounds` | `two radii at the same point of the edge`, `... at the vertex of a closed edge` |
| point hors de l'arête (`at` hors de [0, 1]), rayon pas plus grand que `tol` | `InvalidBounds` | |
| rayon au-delà d'un rayon de courbure (au plus grand rayon de la loi) | `NonGeneric` | `the radius exceeds a radius of curvature of a face` |
| pas de boule du rayon tangente aux deux faces dans le plan normal | `NonGeneric` | `no ball of the radius touches both faces ...` |
| rayon qui varie plus vite que le centre ne se déplace | `NonGeneric` | `the radius varies faster than the centre of the ball moves (no envelope)` |
| loi prolongée au-delà d'un sommet devenue nulle ou négative | `NonGeneric` | `the radius law continued past a vertex is not positive` |
| contact hors de sa face, contact qui traverse la couture d'une B-spline fermée sans boucler | `NonGeneric` | ceux du palier 4 |

## 6. Tests

`tests/tests_brep_variable_fillet.cpp`, 6 cas. Le contrôle commun
(`expect_variable`), à 65 stations, vérifie :

- le centre dans le plan normal à l'arête ;
- les contacts à r(t) du centre, le long des normales des faces ;
- la surface tangente aux faces (sinus sous `angular_tol`) ;
- en travers, les points de la surface à r du centre et sur le cercle
  caractéristique ((q − c) · ĉ′ = −r r′/‖c′‖), sauf au ras d'un point anguleux ;
- pour chaque fuseau : sphère centrée sur la ligne des centres, au rayon du
  point ; arcs par les contacts, sur la sphère, sur leurs pcourbes ; tangence
  aux morceaux le long des arcs ; largeur non nulle.

| Test | Contenu |
|---|---|
| `radius_law` | linéaire et lisse : valeurs aux points, linéaire entre eux ; lisse C1 aux points, monotone et sans dépassement entre eux, pentes nulles aux bouts ; loi fermée C1 à la fermeture ; un seul point ; refus |
| `box_edge` | arête de boîte, loi 0,1 → 0,3 (à 40 %) → 0,15, linéaire et lisse : centres en (x, r, 1 − r) ; linéaire : deux morceaux et un fuseau ; prolongements (constant pour le lisse, de pente 0,25 pour le linéaire) |
| `constant_law` | même rayon partout : le cylindre exact du rayon constant |
| `cylinder_rim` | bord d'un cylindre (arête fermée plan-cylindre), loi lisse fermée par trois points ; linéaire : trois morceaux, trois fuseaux, dont un au sommet |
| `bspline_face` | arête du dessus du prisme d'une bosse B-spline, loi 0,05 → 0,12 → 0,08, linéaire et lisse |
| `refusals` | pas de point sur l'arête, deux rayons au même point, point hors de l'arête, rayon nul, congé plus large que les faces, deux rayons au sommet d'une arête fermée |

Les volumes de référence (intégration de la réunion des boules) viennent avec
la topologie, dans la PR 2 : sans faces, pas de solide à mesurer.

## 7. Limites

- Sur une arête fermée autour d'une B-spline fermée non périodique (pied
  d'aube), la loi part du sommet de l'arête, pas de la couture : le contact
  n'est en général pas dans les bornes de la face (refus). La PR 2 décalera le
  départ de la boucle à la couture, comme au palier 4.
- Les congés qui s'annulent (faces tangentes à un sommet) ne sont pas
  traités avec une loi.
- Le prolongement suppose l'arête parcourue à vitesse constante au-delà de ses
  sommets (vitesse au sommet) : l'abscisse n'y sert qu'à la face d'extrémité.

## 8. Suite

PR 2 : `fillet(m, s, edges, RadiusLaw)`, avec les faces des fuseaux, les chaînes
(loi sur la chaîne), les boucles, les coins au même rayon, les volumes de
référence, Python, la galerie et la comparaison à OCCT.

## 9. Points à valider

1. **Loi en `BSCfunction`** de l'abscisse (degré 1 ou 3 à nœuds doubles).
2. **Pentes nulles aux bouts d'une loi lisse** ouverte : elle rejoint le
   prolongement constant en C1.
3. **Loi linéaire prolongée avec sa pente** au-delà de ses points extrêmes
   (révision de la question 8 du palier, décidée le 19 nov. 2026).
4. **Centre du cercle caractéristique tiré des contacts** (sans r′).
5. **Morceaux et fuseaux** dans `EdgeBlend` (`pieces`, `lunes`, `surface`
   vide), plutôt qu'une seule surface.
6. **Repère de la sphère d'un fuseau** : x vers son milieu, z en travers des
   contacts.
7. **Volumes de référence dans la PR 2**, avec la topologie.
