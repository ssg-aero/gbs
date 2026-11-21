# Palier 6 — PR 2 : surface à n côtés

| | |
|---|---|
| PR | branche `feat/nsided-surface` |
| Document de référence | [brep_palier6_vertex_blends.md](brep_palier6_vertex_blends.md), § 4.3 à § 4.6 |
| Fichiers | `gbs/nsided.h` (nouveau), `tests/tests_nsided.cpp` (nouveau), `python/gbsBindSurfaceTools.cpp`, `python/tests/test_nsided.py` (nouveau), `docs/sources/design/brep_palier6_vertex_blends.md`, `News.md`, cette note |
| Statut | À valider (§ 9) |

## 1. Objectif

La surface qui remplit un trou à n côtés, tangente aux surfaces qui l'entourent,
dans le noyau gbs (hors BREP). La PR 3 en fera la face des raccords de sommet ;
elle sert aussi à combler un trou quelconque (Python : `gbs.fill_n_sided`).

**Écart au palier.** Le § 4.4 du palier prévoyait une surface transfinie de Kato
comme forme guide, approchée sous contraintes. Mesurée, elle s'enfonce trop
(§ 4.1) : la surface retenue est celle d'**énergie de flexion minimale** sous
les contraintes des côtés. C'est le point 1 du § 9.

## 2. Sources

- K. Kato, « Generation of N-sided surface patches with holes », 1991 :
  essayée comme guide, écartée (§ 4.1).
- P. Salvi, T. Várady, A. Rockwood, « Ribbon-based transfinite surfaces »,
  2014 : rubans, domaine polygonal.
- J. A. Gregory, P. Charrot, 1984 : la compatibilité des torsions aux coins (la
  raison de relâcher les dérivées transverses aux coins, § 4.3).
- Plaques minces (énergie de flexion) : la formulation classique du lissage de
  surfaces (Greiner, Hormann et d'autres) ; ici sur une B-spline par moindres
  carrés.
- D. A. Dunavant, « High degree efficient symmetrical Gaussian quadrature rules
  for the triangle », 1985 : la règle à 6 points de degré 4.

## 3. API

```cpp
template <std::floating_point T> struct Ribbon {
    std::shared_ptr<Curve<T, 3>> curve;  // le côté, sur [t0, t1]
    std::shared_ptr<Curve<T, 3>> cross;  // X(t), tangent à la surface voisine, vers le trou
    T t0, t1;
};
auto fill_n_sided(const std::vector<Ribbon<T>> &, const NSidedOptions<T> & = {}) -> NSidedSurface<T>;
```

`NSidedSurface` : la B-spline bicubique sur [−1, 1]², le polygone régulier de
son domaine (côté i : ruban i), l'écart mesuré aux courbes (`deviation`), le
sinus de l'angle aux normales des rubans (`angle`), le nombre de carreaux par
direction. `NSidedOptions` : `approx_tol`, `angular_tol` (10⁻⁵), `max_spans`
(64), `side_weight`, `cross_weight`, `shape_tol`, `side_samples`. Python :
`gbs.Ribbon`, `gbs.NSidedOptions`, `gbs.NSidedSurface`, `gbs.fill_n_sided`.

Les côtés se suivent (la fin de l'un au début du suivant), dans le sens
trigonométrique vu du côté des normales P′ × X. La surface est à découper par
le polygone : ses côtés sont des segments du domaine (la PR 3 en fera les
pcourbes).

## 4. Démarche

### 4.1 Ce qui a été essayé

Mesuré sur des trous découpés dans des surfaces connues (paraboloïde, selle,
plan, octant de sphère), le bord image d'un polygone régulier :

| Construction | Paraboloïde, 3 à 6 côtés | Octant de sphère | Défaut |
|---|---|---|---|
| Kato, rubans linéaires (guide du palier) | s'enfonce de 0,19 (carré) à 0,33 (pentagone), sur une profondeur de 0,5 | 4,6 % du rayon | les rubans prolongent les plans tangents |
| Kato, rubans quadratiques (courbure en travers) | 2,5 à 13 % | 3,2 % | paramètres des côtés incohérents entre eux |
| énergie de flexion, normales seules | 3 % à 64 carreaux, non convergé | non convergé | les paramètres glissent le long du bord |
| **énergie, normales et dérivées transverses compatibles (retenue)** | **10⁻⁵** | **4 à 8 %** | aucun sur les quadriques ; voir § 8 |

### 4.2 Domaine

Le polygone régulier à m côtés inscrit dans le cercle unité, côté i de Vᵢ à
Vᵢ₊₁, le paramètre [t0, t1] du ruban i porté linéairement dessus.

### 4.3 Dérivées transverses compatibles

Le long du côté i, la dérivée de la surface selon la normale intérieure du côté
dans le domaine est a(s) P′(s) + b(s) X(s), a et b linéaires en s. À chaque coin,
les dérivées le long des deux côtés du domaine sont les tangentes des deux
courbes qui s'y rencontrent : la normale du côté i, combinaison des deux
directions, donne la dérivée transverse au coin, décomposée sur (P′, X). Refus
si b ≤ 0 à un coin (X ne pointe pas vers le trou).

Les torsions que demandent deux côtés à leur coin commun ne s'accordent pas en
général (la seconde dérivée croisée d'une B-spline C2 est unique) : imposées
partout, elles empêchent les tolérances au bord (octant de sphère : 4·10⁻³ en
position à 64 carreaux). Elles sont donc tenues fortement au milieu des côtés,
plus du tout aux coins (poids 16 s²(1 − s)²).

### 4.4 B-spline

Bicubique, N × N carreaux uniformes sur [−1, 1]², par moindres carrés
pondérés, chaque ligne pondérée par la racine de sa part d'un côté ou de l'aire
(des sommes indépendantes de l'échantillonnage) :

- **côtés** (poids 10⁵) : les points sur les courbes, les dérivées normales au
  ruban (n · S_u = n · S_v = 0) ;
- **dérivées transverses** (poids 10³, nul aux coins) : celles du § 4.3 ;
- **énergie de flexion** ∫ |S_uu|² + 2 |S_uv|² + |S_vv|² sur le polygone, par la
  règle de Dunavant à 6 points sur les triangles d'un éventail depuis le centre,
  chacun coupé en (2N)² ; sur le reste du carré, mille fois plus faible (elle
  tient ses coins).

Une seule matrice creuse, factorisée par `SimplicialLDLT` (équations normales).

**Quadrature.** Une première version sommait l'énergie aux points d'une grille
intérieurs au polygone : c'était l'énergie d'un polygone un peu plus petit, et
la surface s'écartait de 10⁻³ du paraboloïde qu'elle doit retrouver, l'écart ne
décroissant qu'en O(h). Sur les triangles : quelques 10⁻⁶.

### 4.5 Raffinement

N doublé à partir de 4 jusqu'à ce que :

- la surface soit à `approx_tol` des courbes et sa normale à `angular_tol` de
  celle des rubans, le long de chaque côté (64 points par carreau, 256 au
  moins) ;
- et qu'elle ne bouge plus à l'intérieur, le long de sa normale, de plus de
  `shape_tol` d'un doublement au suivant (grille de 25 × 25 points du carré
  dans le polygone ; par défaut 2·10⁻³ de la taille du trou, 10 `approx_tol` au
  moins).

Refus (`std::runtime_error`) au-delà de `max_spans`.

## 5. Refus

`std::invalid_argument` : moins de trois rubans, des côtés qui ne se suivent
pas, une direction transverse nulle ou le long de sa courbe, ou qui ne pointe
pas vers le trou à un coin ; des tolérances non positives.
`std::runtime_error` : tolérances non atteintes à `max_spans` (des rubans dont
les normales diffèrent à un coin, par exemple), système non factorisable.

## 6. Tests

`tests_nsided` (C++, 5 cas) :

| Cas | Attendu | Mesuré |
|---|---|---|
| plan, 3 à 6 côtés | la surface est le plan | à 10⁻⁵ |
| paraboloïde et selle, 3 à 6 côtés | tolérances au bord ; la quadrique retrouvée à l'intérieur (biharmonique) | à 10⁻⁴ (assertion) ; 7·10⁻⁶ à 3·10⁻⁵ mesurés, 8 carreaux, 0,02 à 0,05 s |
| octant de sphère (arcs de grands cercles) | tolérances au bord ; près de la sphère | au bord 5·10⁻⁷ ; à l'intérieur 8,3 % du rayon (dérivées transverses de longueur variable), 4,1 % (constantes, test Python) ; 64 carreaux, 2 s |
| rotation des côtés | même surface dans l'espace, le domaine tourné | à 10⁻⁴ (la grille de la B-spline ne tourne pas avec le domaine) |
| refus | deux rubans, côtés dans le désordre, X le long de la courbe, normales en désaccord à un coin | exceptions |

`python/tests/test_nsided.py` : l'octant de sphère depuis des `EllipseCurve`,
refus.

## 7. Note d'illustration

Rendus des surfaces de remplissage dans la surface d'où le trou est découpé,
en bandes zébrées : leur continuité à travers les bords montre la tangence
(hexagone dans un paraboloïde, pentagone dans une selle, octant de sphère).

## 8. Limites

- **Fidélité à la sphère.** Sur l'octant, la surface d'énergie minimale n'est
  pas la sphère : 4 à 8 % du rayon selon les longueurs des X. Le cas compte
  pour la PR 3 : quand les rayons d'un coin s'égalisent, le raccord doit tendre
  vers la sphère du coin. Deux pistes : y garder la sphère jusqu'à un écart de
  rayons, ou construire le raccord d'une famille de sphères plutôt que d'un
  remplissage générique (point 5 du § 9).
- **Coût.** L'octant demande 64 carreaux (4 489 pôles, 2 s) : la convergence est
  lente près des coins, où les angles droits viennent des angles de 60° d'un
  triangle équilatéral. Un raffinement local, ou un domaine aux angles des
  coins de l'espace (Salvi et Várady), le réduirait.
- **Équations normales.** Le poids 10⁵ des côtés limite la précision à quelques
  10⁻⁶ sur une surface exacte (le plan) ; une élimination des contraintes la
  rendrait exacte.

## 9. Points à valider

1. **Énergie de flexion minimale** à la place de la surface guide de Kato
   (écart au § 4.4 du palier), après les mesures du § 4.1.
2. **Dérivées transverses compatibles aux coins**, tenues au milieu des côtés
   et relâchées aux coins ; seules les normales y sont imposées.
3. **Domaine régulier** et B-spline découpée (sans changement).
4. **Critère d'arrêt sur la forme** relatif (2·10⁻³ de la taille du trou), les
   tolérances au bord restant absolues.
5. **Raccord presque sphérique** : choix renvoyé à la PR 3 (§ 8).
