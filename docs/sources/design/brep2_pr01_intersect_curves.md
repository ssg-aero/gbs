# Palier 2 — PR 1 : boîtes englobantes, intersection courbe-courbe

| | |
|---|---|
| PR | branche `feat/intersect-curves` |
| Document de référence | [brep_palier2.md](brep_palier2.md), § 2.1 à 2.3 (commité dans cette PR) |
| Fichiers | `gbs/intersect/box.h`, `gbs/intersect/curve_pieces.h`, `gbs/intersect/curves.h` (nouveaux), `tests/tests_intersect_curves.cpp` (nouveau), `docs/sources/design/brep_palier2.md` (nouveau), `docs/index.rst`, `News.md` |
| Statut | À valider, **avec le document de conception du palier 2** |

## 1. Objectif de la PR

Poser la première brique du palier 2 : trouver **toutes** les intersections de
deux courbes sur des plages bornées. La découpe des faces (PR 7) s'en servira
pour croiser, dans l'espace (u, v), les pcurves des arêtes de section avec
celles du bord. La classification d'un point dans une face (PR 5) s'en servira
aussi.

Le document de conception du palier 2 entre dans `master` avec cette PR :
il n'y a pas de PR de documentation seule.

## 2. API

```cpp
template <T> struct IntersectOptions { T tol = extrema_tol; T angular_tol = 1e-3; size_t max_pairs = 200000; };
template <T, dim> struct CurveCurveHit   { T t1, t2; point<T, dim> p; bool tangent; };
template <T>      struct CurveOverlap    { T t1a, t1b, t2a, t2b; };   // t2a > t2b : sens opposés
template <T, dim> struct CurveCurveResult { vector<CurveCurveHit> hits; vector<CurveOverlap> overlaps; bool complete; };

auto intersect(c1, a1, b1, c2, a2, b2, IntersectOptions = {}) -> CurveCurveResult;   // plages données
auto intersect(c1, c2, IntersectOptions = {}) -> CurveCurveResult;                   // bornes des courbes

template <T, dim> struct Box { min, max; add, inflate, overlaps(b, gap), contains(b, gap), diagonal };
```

## 3. Choix d'architecture et justification

### 3.1 Morceaux de courbe (`CurvePieces`)

| Courbe | Morceaux | Boîte | Moitiés |
|---|---|---|---|
| `BSCurve`, `BSCurveRational` | **segments de Bézier** (restriction à la plage, puis insertion des nœuds jusqu'à la multiplicité du degré) | boîte des pôles cartésiens : **garantie** (enveloppe convexe, poids positifs) | **de Casteljau** en coordonnées homogènes : exact |
| toute autre courbe (cercle générique, `CurveOnSurface`…) | un morceau échantillonné (8 intervalles) | boîte des échantillons, élargie de deux fois le plus grand écart de la courbe aux cordes : **heuristique** | nouvel échantillonnage |

Chaque morceau connaît sa **planéité** : la distance maximale de ses points
(pôles ou échantillons) à la corde de ses extrémités.

### 3.2 Algorithme

1. **Filtre.** Les paires de morceaux dont les boîtes sont plus éloignées que
   la tolérance sont écartées.
2. **Portion commune.** Si le plus petit morceau tient dans la boîte de
   l'autre, on projette cinq de ses points sur l'autre courbe. S'ils sont tous
   dans la tolérance, le morceau est une **graine de recouvrement** et n'est
   plus découpé. C'est ce qui évite l'explosion du nombre de paires le long
   d'une portion commune.
3. **Candidat.** Deux morceaux plats (écart à la corde ≤ 1e-3 × corde)
   donnent un point de départ, pris aux points les plus proches de leurs
   cordes. Il est affiné par **Gauss-Newton** sur |c1(s) − c2(t)|², sur les
   courbes entières, et gardé si la distance atteinte est sous la tolérance.
4. **Presque parallèles.** Deux morceaux plats presque parallèles (sinus de
   l'angle < 0,1) peuvent se croiser deux fois. Ils sont encore découpés
   jusqu'à 10 × la tolérance avant de donner un candidat.
5. **Découpe.** Sinon, on coupe en deux le plus grand des morceaux non plats.

Après le parcours :
- les graines qui se suivent sur les deux courbes forment les
  **recouvrements** ;
- les candidats sont triés et **dédoublonnés**. Deux solutions sont les
  mêmes si les courbes restent confondues entre elles, ce qui garde les deux
  solutions distinctes d'une courbe fermée ;
- les croisements aux bouts d'un recouvrement **l'étendent** ; ceux qu'il
  contient sont retirés ;
- une courbe fermée qui coupe l'autre en son point de départ ne compte qu'**une**
  solution.

### 3.3 Tangence ou recouvrement

C'est le seul point délicat, et il faut le valider. Une droite tangente à un
cercle de rayon R reste dans la tolérance sur une longueur d'environ
2√(2 R tol), soit 3 mm pour R = 1 m et tol = 1e-6 m. Par le seul critère de
distance, c'est une portion commune.

On les distingue par le **profil de distance** le long de la portion. Sur un
vrai recouvrement, la distance varie très peu. Sur une tangence, elle part de
presque zéro au contact et monte jusqu'à la tolérance aux bouts. Si l'écart
entre la plus grande et la plus petite distance dépasse tol / 4, la portion
est une tangence. Elle est alors remplacée par un croisement au point le plus
proche, marqué `tangent` (sinus de l'angle < `angular_tol`).

Ce critère classe aussi correctement deux courbes parallèles décalées de
tol / 2 : elles forment un recouvrement, et non une tangence.

### 3.4 Dimension quelconque

L'algorithme ne dépend pas de la dimension. Le document de conception parlait
de 2D (pcurves). En 3D, il donne les points où deux courbes se rencontrent dans
la tolérance (droites gauches : aucune solution). Cela servira à la section
pour fusionner les sommets des arêtes.

### 3.5 Dépendances

Aucune nouvelle. Pas de nlopt : Gauss-Newton 2 × 2 et la projection d'un point
par Newton sont écrits sur place. Les extrema existants (`extrema.h`) ne sont
pas touchés (question 11 du document de conception).

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| courbe-courbe en 2D | toute dimension | même code ; utile à la section en 3D |
| formule fermée pour les droites | pas de cas particulier | un segment de degré 1 est un Bézier plat : une seule paire, Newton converge en un pas |
| — | critère tangence / recouvrement (§ 3.3) | un recouvrement de √tol sinon pour toute tangence |
| — | `complete = false` au-delà de `max_pairs` | garde-fou contre une configuration pathologique |

## 5. Tests (`tests_intersect_curves`, 7 cas, 251 assertions, 50 ms)

| Test | Couvre |
|---|---|
| `boxes` | boîte vide, ajout, écart, inclusion |
| `segments` | croisement exact, parallèles, jonction en T à une extrémité, segments colinéaires dans le même sens et en sens opposé (recouvrement exact) |
| `circles_and_lines` | cercle NURBS et cercle générique (chemin échantillonné) : sécante, **tangente** (une solution marquée tangente), diamètre passant par le point de départ (deux solutions, pas trois) ; deux cercles ; arc posé sur le cercle (recouvrement) ; courbe avec elle-même |
| `bsplines` | sinusoïde cubique contre l'axe (5 solutions à kπ) et contre des horizontales ; deux sinusoïdes ; restriction à une demi-plage |
| `curves_3d` | droites concourantes et gauches, cercle incliné et une de ses cordes |
| `robustness` | cercle générique sur cercle NURBS (recouvrement sur toute la plage, à travers le point de départ) ; croisements rasants jusqu'à y = 0,9999 ; carré coupé par un cercle (8 solutions) |
| `limits` | plage restreinte, tolérance plus large, plafond de paires |

## 6. Points à valider

1. **Le document de conception du palier 2** ([brep_palier2.md](brep_palier2.md)),
   en particulier la question 1 (booléens) ; les recommandations sont
   appliquées tant qu'elles ne sont pas tranchées. Cette PR ne dépend que de
   la question 2 (emplacement `gbs/intersect/`).
2. **Bézier et de Casteljau** pour les B-splines, **échantillonnage** pour le
   reste, avec une boîte heuristique.
3. **Critère tangence / recouvrement** par le profil de distance (§ 3.3).
4. **Dimension quelconque** au lieu du 2D seul.
5. **Solutions aux extrémités** incluses, une seule solution au point de
   départ d'une courbe fermée.
