# Noyau BREP natif — conception du palier 5 : congés à rayon variable

| | |
|---|---|
| Statut | Questions du § 8 validées (19 nov. 2026), avec le fuseau aux points anguleux d'une loi linéaire (§ 3.5) ; chaque PR valide sa note |
| Périmètre | congés dont le rayon suit une loi le long d'une arête ou d'une chaîne d'arêtes tangentes : valeurs aux sommets et en des points intermédiaires, variation linéaire ou lisse entre elles ; géométrie (enveloppe exacte des boules, surface approchée), topologie (extrémités, chaînes, coins), Python, galerie, comparaison à OpenCascade |
| Dépend de | palier 4 : géométrie des congés (`fillet.h`, PR 5 et 7), topologie (`blend.h`, PR 6 et 7) |
| Hors périmètre | chanfreins à distance variable, raccords de sommet quelconques (plus de trois arêtes, rayons différents au sommet), congés de faces (*face blend*), sections coniques ou à corde constante |

Ce document décrit *quoi* construire et *pourquoi*, comme
[brep_palier4.md](brep_palier4.md). Chaque PR a sa note d'architecture
`brep5_prNN_*.md`. C'est le premier point de la feuille de route du palier 4
(§ 9).

**Provenance.** Les algorithmes viennent de la littérature et du comportement de
CATIA et SolidWorks, la cible de qualité ; chaque note cite ses sources.
OpenCascade ne sert que de comparaison en boîte noire dans les tests.

---

## 0. Résumé des décisions proposées

| Sujet | Décision proposée | Alternative écartée |
|---|---|---|
| Loi de rayon | **points de contrôle** (sommet, ou fraction de la longueur d'une arête) et rayon en chacun ; entre eux **linéaire** ou **lisse** (cubique monotone, sans dépassement), au choix ; une **fonction B-spline de gbs** (`BSCfunction`) de l'abscisse | une fonction quelconque r(s) donnée par l'appelant (gardée en interne, pas exposée à la première version) |
| Paramètre de la loi | **abscisse curviligne** le long de l'arête, ou de la chaîne tangente entière | paramètre de la courbe de l'arête (dépend de la paramétrisation, invisible pour l'utilisateur) |
| Station de la boule | la boule de rayon r(s) **dans le plan normal à l'arête** au point d'abscisse s, tangente aux deux faces (Gauss-Newton sur les paramètres des deux contacts) | plan normal à la ligne des centres (celui de la PR 5, qui n'a pas de paramètre commun avec la loi) |
| Surface | **enveloppe exacte des boules** : en chaque station, l'arc du cercle caractéristique entre les deux contacts (plan normal à c', décalé de −r r′/‖c′‖² c′), approché par la B-spline rationnelle de la PR 5 (Hermite, raffinement local, `approx_tol`, `angular_tol`) ; elle se réduit au congé de la PR 5 quand r est constant | l'arc de centre c entre les contacts (tangent aux faces, mais pas l'enveloppe des boules : la section ne suit pas la boule qui roule) |
| Point anguleux d'une loi linéaire | la boule y a deux cercles caractéristiques : le congé est coupé en **morceaux**, raccordés par le **fuseau** de la boule entre les deux cercles (face de sphère), G1 partout ; la loi linéaire est prolongée **avec sa pente** au-delà de ses points extrêmes (pas de fuseau aux sommets) | arrondir les angles de la loi ; refuser les points intermédiaires ; section dans le plan normal (pli) |
| Cas exacts | aucun primitive nouvelle ; le **rayon constant** garde ses cylindres et tores exacts (le chemin de la PR 5) | cônes obliques exacts pour une loi linéaire entre deux plans (rationnels de degré (2, 1)) : plus tard s'il le faut |
| Topologie | **celle du palier 4**, inchangée : contacts, extrémités libres (courbe d'extrémité cheminée sur la surface), chaînes, boucles ; un coin seulement si les trois lois donnent **le même rayon au sommet** | raccord de sommet pour des rayons différents (*vertex blend*, feuille de route) |
| Variation trop rapide | **refusée** là où ‖r′‖ ≥ ‖c′‖ (les boules s'emboîtent, plus d'enveloppe) ou où le rayon dépasse un rayon de courbure d'une face | réparée |
| API | `fillet(m, s, edges, RadiusLaw, options)` en plus de la forme à rayon constant ; `fillet_geometry` avec une loi | une loi par arête seulement (pas de chaîne) |

---

## 1. Objectifs et non-objectifs

### 1.1 Objectifs

- Un congé dont le rayon varie le long d'une arête : de r₁ à r₂ entre ses
  sommets, ou par des points de contrôle intermédiaires, comme les congés à
  rayon variable de CATIA (« Variable Radius Fillet », variation linéaire ou
  cubique) et de SolidWorks (« Variable Size Fillet », transition droite ou
  lisse, points de contrôle en pourcentage de la longueur).
- La même loi le long d'une chaîne d'arêtes tangentes, continue d'une arête à
  la suivante.
- Les extrémités du palier 4 : libres (coupées par la face d'extrémité), en
  chaîne, en boucle (arête fermée), en coin (trois congés au même rayon au
  sommet).
- Toutes les faces que le palier 4 accepte : primitives, B-splines, faces
  fermées non périodiques (pied d'aube).
- Python, galerie, comparaison à `BRepFilletAPI_MakeFillet` (lois linéaires
  `Add(R1, R2, E)` et par points `SetRadius`).

### 1.2 Non-objectifs

- Chanfreins à distance variable.
- Raccords de sommet entre des rayons différents (*vertex blends*), plus de
  trois arêtes à un coin : feuille de route du palier 4, point 2.
- Sections autres que circulaires (coniques, à corde constante, à courbure
  continue) : options de CATIA et SolidWorks, plus tard.
- Loi donnée par une fonction de l'appelant : interne à la première version.

---

## 2. Loi de rayon

### 2.1 API

```cpp
struct RadiusPoint
{
    EdgeId edge;  ///< l'arête qui porte le point
    T at;         ///< fraction de la longueur de l'arête depuis son premier sommet, dans [0, 1]
    T radius;
};
enum class RadiusVariation { Linear, Smooth };
struct RadiusLaw
{
    std::vector<RadiusPoint> points;  ///< au moins un par chaîne ; aux sommets : at = 0 ou 1
    RadiusVariation variation{RadiusVariation::Smooth};
};
auto fillet(Model<T> &, SolidId, std::span<const EdgeId> edges, const RadiusLaw<T> &, const FilletOptions<T> & = {}) -> BuildResult<SolidId>;
auto fillet_geometry(const Model<T> &, SolidId, EdgeId, const RadiusLaw<T> &, const FilletOptions<T> & = {}) -> BuildResult<EdgeBlend<T>>;
```

Le rayon constant reste `fillet(m, s, edges, r)`, inchangé.

### 2.2 Le long d'une chaîne

Les arêtes données sont groupées en chaînes tangentes, comme au palier 4. Sur
chaque chaîne, l'abscisse s court d'un bout à l'autre (une chaîne fermée : son
tour, depuis un sommet). Les points de contrôle de la chaîne sont rangés par
abscisse ; deux points au même endroit (un sommet donné par ses deux arêtes)
doivent avoir le même rayon. Avant le premier et après le dernier point, la loi
continue avec sa pente : constante pour une loi lisse (pente nulle aux bouts),
droite pour une loi linéaire (qui n'a ainsi pas d'angle à ses bouts). Une arête
ou une chaîne sans point de contrôle est refusée (`InvalidBounds`).

La loi est une fonction B-spline de gbs (`BSCfunction`) de l'abscisse : de
degré 1 pour une loi linéaire, de degré 3 à nœuds doubles pour une loi lisse
(`radius_detail::radius_function`, `gbs-brep/radius_law.h`).

### 2.3 Variation

- **Linéaire** : par morceaux entre les points ; la pente saute aux points, où
  la boule qui roule balaie un fuseau de sphère (§ 3.5). C'est la variation de
  CATIA « linear » et de SolidWorks « straight transition ».
- **Lisse** : cubique d'Hermite par morceaux, pentes de Fritsch et Carlson
  (monotone) : C1, sans dépassement entre deux points de contrôle, jamais
  sous le plus petit des deux rayons voisins. Sur une chaîne fermée, la loi se
  referme avec sa pente.

---

## 3. Géométrie

### 3.1 Stations

Une **station** à l'abscisse s : e(s) le point de l'arête, t(s) sa tangente,
r = r(s). On cherche (u₁, v₁, u₂, v₂) tels que les centres des boules tangentes
aux deux faces coïncident et soient dans le plan normal à l'arête :

  S₁(u₁, v₁) + r n₁ = S₂(u₂, v₂) + r n₂ = c,  (c − e(s)) · t(s) = 0,

n₁, n₂ les normales vers la boule (le côté de la PR 5 : dans la matière le long
d'une arête convexe, dehors le long d'une concave). Quatre équations, quatre
inconnues : Gauss-Newton, comme les stations de la PR 5, depuis la station
voisine. La première station part de la station à rayon constant r(s₀) de la
PR 5 (ligne des centres par intersection des décalés). Les contacts sont
pᵢ = c − r nᵢ.

Le plan normal à l'arête donne un paramètre commun à la loi (l'abscisse le long
de l'arête) et aux stations. Le plan de la PR 5, normal à la ligne des centres,
n'en a pas : la ligne des centres n'est connue qu'une fois les rayons placés.

### 3.2 Enveloppe

La boule de centre c(s) et de rayon r(s) touche chaque face en pᵢ(s). Leur
enveloppe est la surface du congé : en chaque s, le **cercle caractéristique**,
intersection de la sphère et du plan (q − c) · c′ = −r r′, de centre
c − (r r′ / ‖c′‖²) c′ et de rayon r √(1 − r′² / ‖c′‖²). Les deux contacts y
sont : dérivant ‖c − Sᵢ‖ = r le long de s, nᵢ · c′ = r′, donc
(pᵢ − c) · c′ = −r r′. La section du congé est l'arc de ce cercle entre les
deux contacts, du côté de l'arête ; la surface est tangente aux faces le long
des contacts (l'enveloppe est tangente à chaque boule, chaque boule à chaque
face).

Pour r constant, c'est le cercle de la PR 5 (r′ = 0, plan normal à c′). Si
‖r′‖ ≥ ‖c′‖, une boule contient la suivante et il n'y a plus d'enveloppe :
refus (« the radius varies faster than the centre moves »).

c′ et r′ sont pris par différences centrées aux stations, comme les dérivées
des pôles homogènes de la PR 5.

### 3.3 Surface

La B-spline rationnelle de la PR 5 : arcs exacts (quadratiques, un ou deux
segments) aux stations, pôles homogènes interpolés le long par Hermite, ligne
des centres et contacts aussi, stations aux nœuds des faces moins que C3,
raffinement local jusqu'à la moitié de `approx_tol` et de `angular_tol`. En
plus : une station à chaque point de contrôle d'une loi linéaire, où la pente
de r saute (le congé n'y est que C0 le long). Le rayon constant garde le chemin
de la PR 5 et ses surfaces exactes (cylindres, tores).

### 3.4 Prolongement et extrémités

Au-delà des sommets d'une arête ouverte, le congé est prolongé de
`extension` fois le plus grand des rayons aux sommets, la loi continuant avec sa
pente (§ 2.2) : la face d'extrémité le coupe comme au palier 4. Un rayon qui
deviendrait nul ou négatif dans ce prolongement est refusé. Un congé qui
s'annule là où les faces deviennent tangentes reste traité comme au palier 4.

### 3.5 Points anguleux : le fuseau

Là où la pente r′ saute (point intermédiaire d'une loi linéaire, ou fermeture
d'une loi linéaire autour d'une arête fermée), la boule a deux cercles
caractéristiques, l'un pour r′ avant le point, l'autre pour r′ après. Tous deux
passent par les deux contacts, et la boule qui roule balaie entre eux un
**fuseau** de sa sphère : sur l'arête de boîte des tests (pentes 0,25 puis
−0,125), environ 0,36 rad de large.

Le congé est donc fait de **morceaux** entre les points anguleux, chacun
approché comme un congé lisse, ses stations aux bouts prises du côté du
morceau (cercle caractéristique à gauche ou à droite). Entre deux morceaux, le
**fuseau** est une face de la sphère du point (`SphericalSurface`, sa couture et
ses pôles loin du fuseau), bornée par les deux arcs, qui sont les bords des
morceaux. Les contacts et la ligne des centres courent d'un morceau à l'autre
(C0 au point). Le congé reste G1 partout : chaque morceau est tangent à la
boule le long de son arc.

Côté géométrie (`EdgeBlend`) : `pieces` (surface, intervalle, contacts sur la
surface) et `lunes` (paramètre, sphère, arcs, leurs pcourbes sur la sphère),
`surface` restant vide. La topologie (PR 2) pose les fuseaux comme les sphères
de coin.

---

## 4. Topologie

Celle du palier 4, sans changement de fond :

- **contacts** et faces des arêtes : suites de co-arêtes remplacées ;
- **extrémités libres** : la courbe d'extrémité est cheminée sur la surface du
  congé, quelle qu'elle soit ;
- **chaînes** : la loi est continue d'une arête à la suivante, les deux congés
  finissent sur le même arc au sommet (la station du sommet est commune) ;
- **boucles** : arête ou chaîne fermée, loi refermée ;
- **coins** : la sphère du palier 4 si les trois lois donnent le même rayon au
  sommet (à `tol` près) ; sinon refus (raccord de sommet, feuille de route).

---

## 5. Python

`gbs.brep.RadiusLaw`, `RadiusPoint`, `RadiusVariation` ; `brep.fillet(m, s,
edges, law)` et `brep.fillet_geometry(m, s, edge, law)`. Un raccourci pour le
cas le plus courant, de r₁ au premier sommet à r₂ au second, linéaire :
`RadiusLaw.linear(edge, r1, r2)`.

---

## 6. Validation

| Contrôle | Attendu |
|---|---|
| validité | `check` sans défaut, coques fermées |
| volumes | arête de boîte, loi linéaire et lisse : volume de référence par intégration numérique de la réunion des boules dans le coin (à 10⁻⁹) ; pour une loi constante, le volume exact du palier 4 |
| géométrie | le long des contacts : surface tangente aux faces (normales sous `angular_tol`), contacts à r(s) des centres, centres dans les plans normaux à l'arête ; en travers : points de la surface sur le cercle caractéristique (sous `approx_tol`) |
| raccords | chaîne : arcs communs aux sommets ; coin : sphère au rayon commun |
| OCCT | `BRepFilletAPI_MakeFillet` : `Add(R1, R2, E)` pour les lois linéaires sur une arête, `SetRadius` par points ; volumes à 10⁻⁵, mêmes nombres de faces, résultats de gbs relus par OCCT valides |
| galerie | arête de boîte de r₁ à r₂, chaîne d'un oblong à loi lisse, pied d'aube à rayon variable (plus grand au bord d'attaque) |

Le volume d'une loi linéaire n'est pas (1 − π/4) ∫ r² ds : la section du congé
par le plan normal à l'arête n'est pas l'arc de rayon r (le cercle
caractéristique est incliné de r′). La référence intègre donc la réunion des
boules.

---

## 7. Plan des PR

| # | PR | Contenu | Estimation |
|---|---|---|---|
| 1 | loi et géométrie | § 2, § 3 : `RadiusLaw`, stations dans le plan normal à l'arête, enveloppe, surface, morceaux et fuseaux ; `fillet_geometry` avec une loi ; tests géométriques ; ce document ([note](brep5_pr01_variable_geometry.md)) | 16 h |
| 2 | topologie, Python, OCCT | § 4, § 5, § 6 : `fillet` avec une loi (extrémités, chaînes, boucles, coins, faces des fuseaux), volumes de référence, Python, galerie, comparaison OCCT, bilan | 12 h |
| | | **Total** | **≈ 28 h** |

---

## 8. Questions à valider

1. **Loi par points de contrôle** (sommets, fractions de longueur d'arête),
   variation linéaire ou lisse ; pas de fonction de l'appelant à la première
   version. *Recommandé* : ce que CATIA et SolidWorks exposent.
2. **Abscisse curviligne** comme paramètre de la loi, le long de la chaîne
   entière. *Recommandé* : indépendante de la paramétrisation des courbes, et
   continue sur une chaîne.
3. **Variation lisse par Hermite monotone** (Fritsch-Carlson) : C1, sans
   dépassement. *Recommandé* : une spline C2 dépasserait entre des points
   inégaux, et le rayon pourrait passer sous le plus petit d'entre eux.
4. **Lisse par défaut.** *Recommandé* (défaut de SolidWorks ; CATIA a
   « cubique » par défaut).
5. **Stations dans le plan normal à l'arête.** *Recommandé* : la loi et les
   stations ont alors le même paramètre.
6. **Enveloppe exacte des boules** (cercle caractéristique) plutôt que l'arc de
   centre c. *Recommandé* : c'est la boule qui roule, et le rayon constant
   retrouve exactement la PR 5.
7. **Pas de surface exacte nouvelle** (cônes obliques des lois linéaires entre
   plans) ; le rayon constant garde les siennes. *Recommandé.*
8. **Rayon constant au-delà des points extrêmes** et dans le prolongement
   au-delà des sommets. *Recommandé.* Révisé le 19 nov. 2026 pour la loi
   linéaire : elle continue avec sa pente (§ 2.2), pour ne pas avoir de point
   anguleux aux sommets.
9. **Coin seulement au même rayon** au sommet ; sinon refus. *Recommandé* :
   les rayons différents sont un raccord de sommet, point suivant de la feuille
   de route.
10. **Variation trop rapide refusée** (‖r′‖ ≥ ‖c′‖). *Recommandé.*
11. **Deux PR** (géométrie, puis topologie et le reste). *Recommandé.*
