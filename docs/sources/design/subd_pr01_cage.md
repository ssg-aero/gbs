# Subdivision de Catmull–Clark — PR 1 : cage, règles, surface limite exacte

| | |
|---|---|
| PR | branche `feat/subd-catmull-clark` |
| Document de référence | [subd_catmull_clark.md](subd_catmull_clark.md), § 2, 3, 10 |
| Fichiers | `gbs/subd/cage.h`, `refine.h`, `limit.h` (nouveaux) ; `tests/tests_subd.cpp`, `tests/in/subd/` (cages types et valeurs d'OpenSubdiv) ; `scripts/subd_reference/` (`cages.py`, `reference.cpp`, `CMakeLists.txt`) ; `python/gbsBindSubd.cpp`, `python/tests/test_subd.py` ; `News.md`, `docs/index.rst` |
| Statut | À valider |

## 1. Objectif

Évaluer **exactement** la surface limite de Catmull–Clark d'une cage, en tout
point et avec ses dérivées. C'est la référence des PR suivantes : les régions
régulières (PR 2) en sont des B-splines exactes, les faces autour des sommets
extraordinaires (PR 3) des approximations mesurées sur elle.

## 2. API

```cpp
namespace gbs::subd
{
    template <std::floating_point T>
    struct Cage { points_vector<T, 3> vertices; std::vector<std::vector<std::size_t>> faces;
                  std::vector<std::array<std::size_t, 2>> creases; std::vector<std::size_t> corners; };

    enum class CageIssue { NoFace, FaceTooSmall, VertexOutOfRange, RepeatedVertex, NonFiniteVertex, NonManifoldEdge,
                           InconsistentOrientation, NonManifoldVertex, IsolatedVertex, InteriorValenceTwo,
                           CreaseNotAnEdge, CornerOutOfRange };
    class CageError : public std::invalid_argument { auto issue() const; auto index() const; };
    enum class VertexRule { Smooth, Dart, Crease, Corner };

    class Topology;  // demi-arêtes, arêtes, éventails, secteurs, règle de chaque sommet ; strict ou non
    template <std::floating_point T> auto refine(const Cage<T> &, std::size_t levels = 1) -> Cage<T>;

    template <std::floating_point T>
    class LimitSurface
    {
    public:
        explicit LimitSurface(Cage<T> cage);                       // CageError
        auto patch_count() const -> std::size_t;                   // 1 par quadrangle, n par n-gone
        auto patch(std::size_t face, std::size_t corner = 0) const -> std::size_t;
        auto patch_face(std::size_t patch) const -> std::array<std::size_t, 2>;
        auto value(std::size_t patch, T u, T v, std::size_t du = 0, std::size_t dv = 0) const -> point<T, 3>;
        auto normal(std::size_t patch, T u, T v) const -> point<T, 3>;
        auto vertex_limit(std::size_t vertex) const -> point<T, 3>;
    };
}
```

Python : `gbs.subd.Cage`, `Topology`, `refine`, `LimitSurface`, `CageError` (un
`ValueError` avec `issue` et `index`), `CageIssue`, `VertexRule`.

Le document de référence annonçait un `BuildError` : c'est une exception,
comme ailleurs dans le cœur de gbs (`std::invalid_argument`, `ValueError` en
Python) ; `BuildResult` reste au BREP.

## 3. Topologie

- Demi-arête h de la face f en position i : de `faces[f][i]` à
  `faces[f][i + 1]` ; celles d'une face sont contiguës. Les arêtes sont
  numérotées dans l'ordre où on les rencontre. Autour d'un sommet, `outgoing`
  et `ccw` parcourent ses faces dans le sens trigonométrique, depuis l'arête de
  bord sur un bord.
- **Refus** (`CageError`, avec la face, le sommet, l'arête vive ou le coin en
  cause) : tous ceux de `CageIssue`. OpenSubdiv accepte un sommet intérieur de
  deux arêtes (deux faces repliées l'une sur l'autre) et lui donne pour limite
  le sommet lui-même ; nous le refusons, la limite n'y a pas de sens.
- **Règles** (Catmark, bord « bord et coins ») : coin si marqué, sur trois
  arêtes vives ou plus, ou au bord d'une seule face ; arête vive sur deux ;
  *dart* sur une ; lisse sinon. Une arête de bord est vive.
- **Secteur** : les faces autour d'un sommet dont la limite forme un seul
  morceau lisse avec une face donnée (tout l'éventail pour un sommet lisse ou
  un dart ; entre deux arêtes vives pour une arête vive ou un coin).
- Une topologie **non stricte** prend un morceau de cage (les faces autour
  d'une face ou d'un sommet) : un sommet de sa bordure peut y être pincé ou
  paraître au bord, sa règle n'a pas de sens ; seuls servent les sommets dont
  toutes les faces sont dans le morceau.

## 4. Raffinement

Règles Catmark (point de face, d'arête lisse ou vive, de sommet lisse, dart,
arête vive, coin). La cage raffinée range les points de sommet d'abord (le
sommet v reste v), puis de face, puis d'arête. Ses faces suivent l'ordre
d'OpenSubdiv : le quadrangle du coin i d'un n-gone part du sommet i ; celui du
coin i d'un quadrangle part de son coin de plus petits (u, v) dans le parent,
dont il garde les directions (l'enfant du coin i y est en position i). Les
cages raffinées une et deux fois sont celles d'OpenSubdiv, indice pour indice.

## 5. Évaluation exacte

### 5.1 Carreaux

Un carreau par quadrangle ; n par n-gone, un par quadrangle de sa première
subdivision (les faces *ptex* d'OpenSubdiv). Sur un quadrangle, (0, 0) est son
premier sommet, u va vers le deuxième, v vers le dernier ; sur le carreau du
coin i d'un n-gone, (0, 0) est le sommet i, u va vers le milieu du côté i, v
vers celui du côté i − 1.

### 5.2 Carreau régulier

Un quadrangle dont les quatre coins sont réguliers est le carreau bicubique
uniforme de ses 16 voisins : un coin intérieur de quatre quadrangles sans
arête vive ; un coin d'arête vive ou de bord avec deux quadrangles du côté du
carreau ; un coin entre deux arêtes vives avec le carreau seul. Les sommets du
pourtour (au-delà des coins) peuvent être quelconques : les points des enfants
du carreau n'en dépendent que par des points de face. Au-delà d'un côté vif,
les points manquants sont les **points fantômes** 2B − I, réfléchis dans les
deux directions à un coin.

### 5.3 Descente

Ailleurs, on prend le carreau et les faces qui touchent ses sommets, on les
subdivise et on descend dans le quart qui contient (u, v), jusqu'à un quart
régulier : environ log₂(1 / max(u, v)) étapes près d'un sommet extraordinaire,
le même calcul que Stam (1998) fait par la base propre. Les dérivées sont
multipliées par 2 par étape et par ordre.

**Repère propre à chaque quart.** Les points de chaque quart sont pris
relativement à son premier coin et mis à l'échelle par une puissance de 2 (les
règles sont affines, la mise à l'échelle est exacte), de sorte qu'ils soient de
l'ordre de 1 ; le point et les dérivées rendus sont ramenés au repère de la
cage par `ldexp`, composante par composante. Près d'un sommet extraordinaire,
les quarts rapetissent géométriquement :

- mesurés depuis l'origine, leurs points garderaient l'arrondi des
  coordonnées de la cage : à 2⁻⁶⁰ d'un sommet de valence 6, sur une cage de
  taille 50, les normales étaient fausses de 3° ;
- à leur propre échelle, ils sortiraient des nombres flottants : à partir de
  2⁻⁵⁰⁰ d'un sommet de valence 3, S_u × S_v passait sous le plus petit double
  (`normal` levait une exception), à 2⁻¹⁰⁰⁰ S_u aussi, et S_uu débordait à
  2⁻⁵²⁰ alors qu'il vaut environ 10¹¹² (revue).

Les dérivées sont maintenant finies tant que leur valeur l'est : d'ordre k, elles
varient comme (2ᵏ λ)^profondeur, λ la valeur propre sous-dominante, ce que les
tests vérifient à 10⁻¹² près jusqu'à 2⁻¹⁰⁷⁰. La profondeur est bornée par le
type flottant (`digits − min_exponent + 64`), non par une constante.

### 5.4 Aux sommets irréguliers

**Position** : V à un coin ; (4 V + P₋ + P₊) / 6 sur une arête vive ;
(n² V + 4 Σ eᵢ + Σ dᵢ) / (n (n + 5)) à un sommet lisse (Halstead et al. 1993,
sur des quadrangles, l'anneau raffiné une fois au besoin) ; au dart, le
vecteur propre à gauche de valeur propre 1 de la matrice de subdivision de son
anneau.

**Normale** : masques tangents de Halstead et al. à un sommet lisse ; ailleurs
(dart, arête vive, coin), les deux vecteurs propres dominants de la matrice de
subdivision du secteur du carreau, quand son spectre est **propre** :
|λ₁| ≥ |λ₂| > |λ₃| et deux vecteurs propres indépendants. Son signe est celui
des normales du carreau tout près du sommet, vers lesquelles elle tend.

Une dérivée en un sommet irrégulier n'existe pas : `std::domain_error`.

### 5.5 Sommets sans plan tangent donné

Spectres (valeurs propres dominantes hors 1) des secteurs de k faces :

| k | arête vive | coin |
|---|---|---|
| 1 | ½, ¼, ¼ | ½, ½, ¼ |
| 2 | ½, ½, ¼ (régulier) | ½, ½, ½ |
| 3 | 0,580 ; ½ ; 0,410 | 0,580 ; ½ ; ½ |
| 4 | 0,611 ; ½ ; ½ | 0,611 ; ½ ; ½ ; ½ |
| 5 | 0,626 ; 0,550 ; ½ | 0,626 ; 0,550 ; ½ ; ½ |
| 6 | 0,635 ; 0,580 ; ½ | 0,635 ; 0,580 ; ½ ; ½ |

Quand λ₂ est lié à λ₃ (arête vive de 1 face ou de 4 d'un côté, coin de
2 à 4 faces, coin marqué sans arête vive), la limite n'a pas de plan
tangent, ou n'y tend qu'en 1 / profondeur (valeur propre défective : à 2⁻⁶⁰,
encore 2,6·10⁻² rad pour une arête vive de 4 faces) : `normal` lève
`std::domain_error`. À partir de 5 faces dans le secteur d'une arête vive ou d'un
coin, deux modes intérieurs dominent ceux des arêtes vives (½) : la normale
donnée est celle que les normales atteignent par l'intérieur du secteur, mais
les arêtes vives arrivent en général hors de ce plan. Ce sont les défauts connus des règles de bord standard
(Biermann, Levin, Zorin 2000) ; les corriger changerait la surface, qui ne
serait plus celle d'OpenSubdiv ni de Blender. Une cage dont les sommets de bord
sont réguliers et les sommets extraordinaires intérieurs n'en a aucun.

### 5.6 Coût

Sur le prisme hexagonal raffiné deux fois (144 quadrangles, macOS arm64,
clang 21, `-O2`) : 7 µs par point au hasard ; près d'un sommet de valence 6,
environ 10 µs par niveau de descente (56 µs à 2⁻⁵, 190 µs à 2⁻²⁰, 390 µs à
2⁻⁴⁰) ; 160 µs pour la normale au sommet. La revue a retiré 15 % (recherche
exacte des arêtes orientées dans un tableau trié au lieu d'une table de
hachage, topologie de la pièce finale gardée, plus de copie de topologie).
Restent, par niveau, deux topologies et le raffinement de toute la pièce, dont
seules les faces autour du quart choisi servent. Si l'ajustement de la PR 3 le
demande : ne raffiner que ces faces, mettre en cache les anneaux raffinés, ou
la base propre de Stam.

En Python, `value`, `normal` et `vertex_limit` relâchent le GIL ; `values` et
`normals` évaluent une liste de points d'un coup et rendent un tableau (n, 3).

## 6. Écarts avec OpenSubdiv

Constatés et vérifiés (valeurs propres en Python, descente de 80 niveaux,
normales voisines) :

| Cas | OpenSubdiv | gbs |
|---|---|---|
| position limite d'un dart | masque lisse : ses poids s'écartent jusqu'à 0,015 à 0,033 de ceux de la limite (valence 6 à 3, sur une somme de 1) | vecteur propre de 1 de l'anneau, exact |
| tangente en travers d'une arête vive de 3 faces ou plus d'un côté | formules de Biermann et al., écrites pour leurs règles modifiées : pas vecteurs propres des règles standard (résidu relatif de 0,09 à 0,21, masque à 14° à 32° du vecteur propre pour k = 3 à 6) | vecteurs propres dominants, ou refus (§ 5.5) |
| tangentes à un coin | différences des deux premières arêtes (« provisoire » dans son code) | vecteurs propres, ou refus |
| surface à moins de 2⁻¹⁰ d'un sommet irrégulier | carreaux de Gregory approchés | exacte à toute profondeur |
| dérivées sur les carreaux d'un n-gone | mesurées dans les paramètres du n-gone : 2ᵒʳᵈʳᵉ fois les nôtres | celles de `value` |
| sommet intérieur de deux arêtes | accepté, limite = le sommet | refusé |

## 7. Validation

`tests/tests_subd.cpp` : 10 cas, 39 102 assertions, 0,8 s.

| Test | Contenu |
|---|---|
| `refusals` | chaque `CageIssue` et l'élément désigné ; les cages types acceptées |
| `rules` | règles de la cage `creases` (3 darts, 2 coins, 8 arêtes vives) et de `grid` ; `hexprism` : douze sommets de valence 3, χ = 2, deux sommets de valence 6 une fois raffiné |
| `refine_as_opensubdiv` | cages raffinées une et deux fois : mêmes sommets (10⁻¹³ de la taille) et mêmes faces, indice pour indice |
| `limit_as_opensubdiv` | positions limites aux sommets du niveau 2 (sauf darts), normales où ses masques sont exacts (10⁻¹¹ rad) |
| `patches_as_opensubdiv` | surface, dérivées premières et secondes en 300 points de ses carreaux réguliers par cage (10⁻¹¹ relatif) |
| `vertex_limits_by_subdivision` | position limite de chaque sommet contre 80 subdivisions de son anneau (darts compris) |
| `normals_at_irregular_vertices` | normale au sommet contre celles du carreau à 2⁻²⁰ et 2⁻⁵⁰ (décroissance en (λ₃ / λ₂)^profondeur) ; refus aux secteurs listés au § 5.5 |
| `deep_near_extraordinary_vertices` | à 2⁻⁶⁰⁰ et 2⁻¹⁰⁷⁰ de sommets de valence 3 et 6 : dérivées premières et secondes finies, dans le rapport 2ᵏ λ d'un niveau au suivant (10⁻¹²), normale égale à celle du sommet (10⁻¹² rad) |
| `continuity` | les deux carreaux d'une arête : mêmes points, mêmes normales à travers une arête lisse |
| `errors` | dérivée à un sommet extraordinaire, hors de [0, 1]², carreau ou sommet inexistant |

Cages types (`scripts/subd_reference/cages.py`, coordonnées perturbées) :
`box` (8 sommets de valence 3), `torus` (régulier), `ngons` (triangles,
pentagone, valence 5), `grid` (bords et coins réguliers), `fans` (arêtes vives
de 3 et 4 faces au bord), `creases` (équateur vif, darts, coins, coin marqué),
`hexprism` (hexagones, sommets de valence 6 une fois raffiné). Valeurs
d'OpenSubdiv 3.7 : `scripts/subd_reference/reference.cpp`, lancé à la main
(§ 10 du document de référence), ses points tirés d'une graine fixée par le nom
de la cage ; la CI lit le JSON.

`python/tests/test_subd.py` : refus, raffinement, surface du prisme hexagonal
contre OpenSubdiv, sommets de valence 6 (position, normale, pas de dérivée),
règles, évaluation groupée.

## 8. Revue

Revue de code de la PR, sur `gbs/subd/*.h`, les tests et les liaisons Python :

| Remarque | Suite |
|---|---|
| `normal` lève une exception au-delà de 2⁻⁵⁰⁰ d'un sommet de valence 3 ; dérivées nulles, infinies ou NaN au-delà de 2⁻¹⁰⁰⁰ | quarts mis à l'échelle par puissances de 2, exposant gardé (§ 5.3) ; test à 2⁻⁶⁰⁰ et 2⁻¹⁰⁷⁰ |
| garde de profondeur fixe (2 200), insuffisante pour un `long double` | bornée par `std::numeric_limits<T>` |
| clés d'arêtes orientées (a << 32) ^ b : collisions au-delà de 2³² sommets | tableau trié des paires, comparées exactement |
| masque du dart sans contrôle du solveur | échec du solveur, valeur propre loin de 1 ou somme nulle : `std::runtime_error` |
| coins de 5 faces et plus : normale donnée, la note disait « 2 faces ou plus » | textes corrigés : refus pour 2 à 4 faces |
| coût par niveau, copies de topologie | § 5.6 |
| GIL tenu pendant les évaluations | relâché ; `values` et `normals` groupés |
| C++23 | `std::ranges::contains`, `std::unreachable` |

## 9. Hors de cette PR

Conversion en B-splines (PR 2 et 3), BREP (PR 4), interpolation et
constructeurs de cage, OBJ (PR 5), arêtes semi-vives.
