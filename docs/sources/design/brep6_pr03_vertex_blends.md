# Palier 6 — PR 3 : raccords de sommet

| | |
|---|---|
| PR | branche `feat/brep-vertex-blend-faces` |
| Document de référence | [brep_palier6_vertex_blends.md](brep_palier6_vertex_blends.md), § 4.2 à § 4.6, § 5 et § 6 |
| Fichiers | `gbs-brep/blend.h`, `gbs-brep/fillet.h`, `gbs/nsided.h`, `tests/tests_brep_vertex_blends.cpp`, `tests/tests_brep_fillet.cpp`, `tests/tests_brep_variable_fillet_solid.cpp`, `tests/tests_nsided.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `python/gbsBindBrepOps.cpp`, `python/gbsBindSurfaceTools.cpp`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier6_vertex_blends.md`, `News.md`, cette note |
| Statut | Validé (§ 9 et § 10, 23 nov. 2026) |

## 1. Objectif

Dernière PR du palier 6 : le **raccord de sommet** (*vertex blend*). C'est une
face B-spline unique, tangente aux congés et aux faces qui l'entourent, aux
sommets dont tous les congés se rencontrent sans sphère exacte. Cela couvre :

- des rayons différents au coin d'une boîte (r₁, r₂, r₃) ;
- quatre arêtes ou plus (sommet de pyramide) ;
- un rayon qui décroît vers le coin (refusé au palier 5) ;
- des arêtes convexes et concaves au même sommet (coin d'équerre) ;
- les **retraits** de l'utilisateur (`FilletOptions::setbacks`), qui donnent
  le coin « valise ».

Le palier 4 refusait tous ces cas. La sphère reste là où elle est exacte. La
PR 3 tranche aussi le point 5 de la note de la PR 2 : quand les rayons
s'égalisent, le raccord tend vers la sphère (§ 4.4). Enfin, elle corrige deux
défauts de `fill_n_sided` (§ 4.6).

## 2. Sources

- Vida, Martin et Várady, « A survey of blending methods that use parametric
  surfaces » (1994), § 6 : raccords de sommet, retraits, courbes de raccord
  (*spring curves*).
- T. Várady, A. Rockwood, « Geometric construction for setback vertex
  blending », *Computer-Aided Design* 29(6), 1997 : retraits, courbes de
  raccord, raccord à 2n côtés.
- Comportement de CATIA (« Edge Fillet », *Blend corner*, *Setback distance*) et
  de SolidWorks (*Setback fillet*) : coin sans retrait par défaut, coin valise
  avec retraits.
- Courbes de courbure d'une surface canal (les sections d'un congé) et
  compatibilité de la seconde forme fondamentale à un coin : la raison des
  retraits supplémentaires du § 4.3.
- OCCT (`BRepFilletAPI_MakeFillet`) en boîte noire.

## 3. API

```cpp
template <std::floating_point T> struct Setback { VertexId vertex; EdgeId edge; T distance{}; };
// FilletOptions<T>::setbacks : std::vector<Setback<T>>
template <std::floating_point T> struct GuideSphere { std::array<T, 3> centre{}; T radius{}; T weight{}; };
// NSidedOptions<T>::sphere : std::optional<GuideSphere<T>>
```

`fillet` (rayon constant ou loi) construit les raccords de sommet sans
fonction nouvelle. Les rayons différents passent par `RadiusLaw`, avec un point
par arête. Un retrait se mesure le long de l'arête depuis le sommet, en
longueur d'arc.

En Python : `brep.Setback(vertex, edge, distance)`, `FilletOptions.setbacks`,
`gbs.GuideSphere(centre, radius, weight)`, `NSidedOptions.sphere`.

## 4. Démarche

### 4.1 Classement

Un sommet de n arêtes toutes congées, sans chaîne qui le traverse, est un
**coin**. Sans retrait demandé, la sphère du palier 4 (§ 4.1 du palier) est
tentée. Si elle est refusée, le modèle est remis en état, et le coin devient un
raccord de sommet. La sphère est refusée dans ces cas : lignes des centres qui
ne concourent pas, rayon qui décroît vers le coin, arêtes convexes et
concaves. Avec un retrait demandé à ce sommet, c'est toujours le raccord.

### 4.2 Croisements et retraits

Autour du sommet, les congés et les faces alternent (`ring`). Sur chaque face,
les contacts des deux congés voisins se croisent en Xᵢ (Gauss-Newton). Chaque
congé s'arrête sur sa **section** à la station la plus éloignée du sommet
parmi les deux où l'un de ses contacts passe un croisement. Un retrait demandé
recule cette station : il est refusé s'il est plus court que le retrait
automatique, ou s'il dépasse l'autre bout de l'arête.

Là où les contacts de deux congés ne finissent pas au même point, une **courbe
de raccord** les relie sur la face. C'est une cubique de Hermite dans les
paramètres de la face, tangente aux deux contacts, avec des tangentes de la
longueur de la corde. Le trou a donc entre n et 2n côtés.

### 4.3 Fins proches

Sur une face, les deux contacts peuvent finir à moins de 7·10⁻³ de la taille du
trou l'un de l'autre. C'est le cas où ils finissent au croisement, ou presque.
Deux cas se présentent :

- **Rayons presque égaux, ou contacts à angle droit** (cosinus sous 0,05). Les
  deux contacts finissent au croisement : les fins y sont légèrement obliques
  au lieu d'être des sections. Cela évite un côté minuscule (§ 4.4).
- **Contacts obliques et rayons différents.** Les deux congés reculent d'un
  quart de leur rayon moyen, et une courbe de raccord arrondit le coin de la
  face. La raison est la suivante. Une section de congé est une courbe de
  courbure. Deux sections qui se rejoignent en un point imposent à la surface
  deux directions principales. Ces directions doivent être orthogonales, sauf
  si les deux courbures sont égales. Si ce n'est pas le cas, aucune surface C2
  n'est G1 le long des deux. C'était le cas au sommet de la pyramide : contacts
  à 58°, angle refusé à 7,6·10⁻⁵ à 64 carreaux.

Un retrait demandé est gardé tel quel : pas de fin au croisement sur ses faces.
S'il ne laisse entre les contacts qu'un écart non nul sous 7·10⁻³ du trou, il
est refusé (point 8 du § 9).

**Deux courbes de raccord courtes autour d'une fin** (revue de code, § 10). Une
courbe de raccord courte seule se remplit bien : 1 % du trou pour
(r, r, r(1 + ε)). Mais la fin d'un congé entre deux courbes de raccord courtes
fait plier le remplissage, ou lui fait manquer la tangence, même à 128
carreaux. C'était le cas des rayons (0,5 ; 0,51 ; 0,52), avec deux courbes de
3 % et 1,4 % du trou, et d'un retrait de 0,505 pour un rayon de 0,5, avec deux
courbes de 0,7 %. Ce congé recule donc jusqu'à ce que ses deux courbes de
raccord mesurent un dixième du trou. Si ce congé a un retrait demandé, le
retrait est refusé quand ses deux courbes resteraient sous 5 % du trou. Avec
un seuil de 5 %, les rayons (0,5 ; 0,51 ; 0,52) demandaient encore 128
carreaux (18 s) ; avec 10 %, 64 carreaux (2,3 s).

### 4.4 Rayons presque égaux : vers la sphère

C'est le point 5 du § 9 de la note de la PR 2. Avec les rayons (r, r, r(1 + ε))
au coin d'une boîte, le raccord doit tendre vers la sphère quand ε tend vers 0.
Deux obstacles se présentaient, et il fallait les deux remèdes :

1. **Côtés minuscules.** Avec les seules sections, les courbes de raccord
   mesurent environ εr. Dans le domaine régulier, chacune occupe pourtant un
   côté entier. La paramétrisation s'y pince, et la normale y devient
   imprécise : refus à ε ≤ 10⁻³ (angle 1,2·10⁻⁴). Les fins aux croisements du
   § 4.3 suppriment ces côtés.
2. **Fidélité à la sphère.** La surface d'énergie minimale s'écarte de 4 à 8 %
   du rayon de la sphère (PR 2, § 8). D'où la **sphère guide** (§ 4.6). C'est
   la sphère la plus proche des côtés du trou, par moindres carrés sur
   |p|² = 2 c · p + k. Son poids, 10⁴, est multiplié par (1 − g)². Ici g est le
   plus grand de deux rapports : l'écart des côtés à la sphère sur 5 % de son
   rayon, et le sinus de l'angle des normales des rubans aux siennes sur 0,2.
   Le poids est nul au-delà de 1, ce qui donne une transition continue.

Autres essais, écartés :

| Essai | Résultat |
|---|---|
| domaine aux côtés proportionnels aux longueurs (§ 8 de la note PR 2) | la grille uniforme ne plie pas sur un côté plus court qu'un carreau : refus pour tous les ε |
| fins obliques partout (d'un croisement à l'autre, n côtés) | continu vers la sphère, mais les fins obliques ne sont pas des courbes de courbure : refus dès ε = 0,03 et pour (0,5 ; 0,7 ; 0,9) |
| sections partout | bon dès ε = 3·10⁻³, refus en dessous |

La bascule entre fins aux croisements et sections se fait à 7·10⁻³ de la taille
du trou. Les fins obliques sont essayées jusqu'à ε = 0,02 et les sections à
partir de 3·10⁻³ : le seuil est entre les deux. Au changement de construction,
le volume saute d'environ 10⁻⁴, soit 0,3 % du volume du coin (point 3 du § 9).

### 4.5 Rubans et face

Chaque côté du trou porte un ruban :

- **Fin d'un congé**, une section ou une fin oblique, du contact sur la face
  d'avant au contact sur la face d'après. Sa direction transverse est la
  dérivée du congé le long de son arête, vers le sommet.
- **Courbe de raccord.** Sa direction transverse est dans le plan tangent de la
  face, normale à la courbe, vers le trou.

Le trou est rempli par `fill_n_sided`, jusqu'à 128 carreaux (§ 8), avec la
sphère guide. La B-spline est découpée par le polygone de son domaine.

Les arêtes gardent leurs courbes exactes : la fin sur le congé, la courbe de
raccord sur la face. L'écart à la B-spline est absorbé dans leur tolérance
(`settle`). La face est orientée contre celle du premier congé, le long de leur
arête commune, comme les fuseaux du palier 5. Les courbes de raccord rejoignent
les contacts dans la découpe de leurs faces (`face_cuts`).

### 4.6 `fill_n_sided` : deux corrections et une option

- **Invariance d'échelle (défaut de la PR 2).** Les lignes de l'ajustement
  variaient différemment avec la taille du trou h : positions en h, normales et
  dérivées transverses en 1, énergie en h². L'octant d'une sphère de rayon 0,1,
  avec une tolérance en proportion, était refusé : la forme y bougeait encore
  de 3,9 % à 64 carreaux. Toutes les lignes sont maintenant rapportées à h.
  Les rayons 0,1, 1 et 10 donnent la même surface à 10⁻⁶ près (test `scale`).
  Les poids ne changent pas à h = 1.
- **Critère de forme.** Le déplacement entre deux doublements mesure surtout
  l'erreur de la surface d'avant. Quand ces déplacements décroissent
  géométriquement, d'un rapport ρ, le reste à venir est déplacement · ρ / (1 − ρ).
  C'est le plus petit de ces deux nombres qui est comparé à `shape_tol`. Le
  critère n'est donc jamais plus strict qu'avant. Sur le coin valise, il
  accepte à 64 carreaux une surface que 128 carreaux ne déplacent que de 5·10⁻⁴.
- **Sphère guide** (`NSidedOptions::sphere`) : chaque point est tiré vers le
  pied c + r d sur la sphère, d la direction depuis le centre prise sur la
  surface d'avant. Il y a deux itérations par doublement, donc trois
  ajustements. Le poids s'annule au bord, comme le carré de la distance aux
  côtés : la tangence au bord reste celle des rubans. Sur l'octant de la sphère
  unité, l'écart à la sphère passe de 8,2 % à 3,5·10⁻⁴. Chaque coordonnée est
  tirée séparément. Une rangée le long du rayon seul (d · S = d · c + r)
  couplerait les trois coordonnées en tout point intérieur : jusqu'à 128
  carreaux, 860 Mo et 26 s au lieu de 340 Mo et 7 s (revue de code, § 10).

## 5. Refus

`NonGeneric` :

- congés et faces qui ne se referment pas autour du sommet ;
- contacts qui ne se croisent pas dans les congés prolongés ;
- retrait plus court que le croisement, ou qui dépasse l'autre bout de l'arête ;
- retrait si proche du croisement que ses courbes de raccord seraient trop
  courtes (sous 7·10⁻³ du trou, ou sous 5 % sur les deux faces de son congé) ;
- courbes de raccord qui restent courtes quel que soit le recul ;
- retrait à un sommet sans raccord de sommet ;
- trou non rempli dans les tolérances (le message de `fill_n_sided` suit).

`InvalidId` : retrait sur une arête non congée ou à un sommet qui n'est pas le
sien. `InvalidBounds` : retrait non positif, ou retrait d'un chanfrein.

## 6. Tests

`tests_brep_vertex_blends` (19 cas, dont 8 nouveaux) :

| Cas | Attendu | Mesuré |
|---|---|---|
| boîte, rayons 0,5 ; 0,7 ; 0,9, trois permutations | `check` sans défaut, 10 faces ; écart et sinus au bord ≤ 10⁻⁵ ; courbure intérieure ≤ 4/r_min ; volumes égaux à 10⁻⁵ | écart 7,6·10⁻⁷, sinus 1,7·10⁻⁶ ; courbure 3,75 (1/r_min = 2), 14,3 sur le bord même (§ 8) ; volume 26,03704 |
| rayons (r, r, r(1 + ε)), ε de 10⁻⁴ à 0,3 | G1 sous 10⁻⁵ ; volume → sphère : à 2·10⁻⁴ de la droite tangente jusqu'à ε = 10⁻², à 5·10⁻⁵ de la sphère pour ε = 10⁻⁴ | tableau ci-dessous |
| rayons presque égaux (0,5 ; 0,51 ; 0,52) | G1 ; 10 faces ; volume entre ceux des coins de rayon 0,5 et 0,52 | 26,52023 (26,53807 et 26,50128) ; 64 carreaux |
| sommet de pyramide, rayons 0,05 à 0,11 | raccord de quatre congés, G1 | 8 côtés (retraits supplémentaires du § 4.3) |
| coin valise, trois retraits de 1,2 puis deux | G1 ; 10 faces ; moins de matière enlevée avec deux retraits | 26,4536 et 26,5403 |
| retraits de 0,8 sur des congés de 0,15 à 0,25 ; retrait égal à l'automatique (0,5) | G1 le long de congés B-spline ; le retrait automatique accepté, le coin inchangé à 10⁻⁴ | sinus 9,4·10⁻⁷ ; 26,53806 |
| équerre : arête concave et deux convexes | G1 ; 12 faces | volume 3,00444 |
| rayon décroissant vers le coin (0,4 → 0,15) | G1 ; volume entre ceux des coins de rayon 0,15 et 0,4 | 26,8447 |
| refus | retraits trop courts, trop proches du croisement (0,501 et 0,52), trop longs, négatifs, mal placés, sur un chanfrein ; trou non rempli (angle de 10⁻⁹) | codes et messages du § 5, modèle inchangé |

Continuité vers la sphère (boîte d'arête 3, r = 0,5 ; sphère : 26,53807) :

| ε | 10⁻⁴ | 10⁻³ | 3·10⁻³ | 6·10⁻³ | 10⁻² | 0,1 | 0,3 |
|---|---|---|---|---|---|---|---|
| volume | 26,53803 | 26,53776 | 26,53716 | 26,53625 | 26,53513 | 26,50739 | 26,42660 |
| fins | croisements | croisements | croisements | croisements | sections | sections | sections |

Tests anciens mis à jour : le coin d'équerre (`tests_brep_fillet`) et le rayon
décroissant vers le coin (`tests_brep_variable_fillet_solid`) étaient des refus.
Ce sont maintenant des raccords. `tests_nsided` gagne `scale` et
`guide_sphere`. En Python, `test_guide_sphere_option` et
`test_fillet_setbacks_option` vérifient les options. La galerie gagne cinq cas : `vertex_blend_box`,
`vertex_blend_suitcase`, `vertex_blend_pyramid`, `vertex_blend_mixed` et
`vertex_blend_dropping`.

## 7. Comparaison à OCCT

Le coin de la boîte à trois rayons, avec `BRepFilletAPI_MakeFillet` et un rayon
par arête. Les deux solides sont valides et ont le même nombre de faces (10).
gbs enlève 0,014 de plus : 26,0370 contre 26,0513, soit 4 % environ du volume
du coin. La surface de coin d'OCCT n'est pas la nôtre ; le test vérifie l'écart
avec une marge de 0,02.

## 8. Limites

- **Couche limite au bord.** La surface n'est que G1 sur son bord : sa courbure
  saute sur une bande étroite le long des côtés. Elle vaut 14,3 sur le bord
  même du coin (0,5 ; 0,7 ; 0,9), contre 3,75 au plus à l'intérieur. Le G2
  est hors périmètre du palier.
- **Coût.** Les coins valise peuvent demander 128 carreaux : environ 400 Mo au
  pic et 5 s, autant pour un trou non rempli, refusé à 128 carreaux. Le sommet
  de pyramide passe à 64 carreaux, à 9·10⁻⁶ en angle pour une tolérance de
  10⁻⁵, donc juste.
- **Seuils empiriques.** Plusieurs seuils ont été réglés sur les cas des tests :
  - la bascule à 7·10⁻³ entre fins aux croisements et sections (§ 4.3), avec
    un saut de volume d'environ 10⁻⁴ ;
  - le recul d'un quart du rayon ;
  - le dixième du trou pour deux courbes de raccord autour d'une fin, et 5 %
    pour les retraits demandés (§ 4.3) ;
  - le poids de la sphère guide (10⁴, 5 % du rayon, sinus 0,2).
- **Retraits et paramètre.** La station d'un retrait suppose que le paramètre
  du congé suit celui de l'arête, ce qui est le cas des congés du palier 4 et 5.
- **Coins mêlés difficiles.** Des contacts qui ne se croisent pas, ou un trou
  non convexe à un coin, sont refusés.

## 9. Points à valider

1. **Sphère tentée d'abord**, le raccord sinon. Avec un retrait demandé à un
   sommet, c'est toujours le raccord.
2. **Retraits minimaux sur les sections**, et courbes de raccord de Hermite
   tangentes aux contacts, comme prévu au § 4.2 du palier.
3. **Fins aux croisements** quand les contacts finiraient à moins de 7·10⁻³ de
   la taille du trou, avec des rayons presque égaux ou des contacts à angle
   droit. Les fins y sont légèrement obliques, et le volume saute d'environ
   10⁻⁴ à la bascule. C'est l'écart au § 4.2 du palier.
4. **Recul d'un quart du rayon** là où des contacts obliques entre rayons
   différents finiraient au même point (compatibilité des courbures, § 4.3).
5. **Sphère guide** pondérée par l'accord des côtés avec une sphère : la
   réponse au point 5 de la note de la PR 2. Le raccord tend vers la sphère
   quand les rayons s'égalisent, à 4·10⁻⁵ près pour ε = 10⁻⁴.
6. **Corrections de `fill_n_sided`** : invariance d'échelle (un défaut de la PR 2),
   critère de forme par extrapolation géométrique, et jusqu'à 128 carreaux pour
   les raccords.
7. **Retrait en longueur d'arc le long de l'arête**, et refusé s'il est plus
   court que le retrait automatique (comme prévu au palier).
8. **Courbes de raccord courtes** (revue de code) : la fin d'un congé entre
   deux courbes de raccord de moins d'un dixième du trou recule jusqu'à ce
   seuil. Un retrait demandé n'est jamais ramené au croisement. Il est refusé
   s'il laisse un écart sous 7·10⁻³ du trou, ou deux sous 5 % autour de son
   congé. Le refus est immédiat, avec un message qui dit de l'allonger ou de
   l'omettre.

## 10. Revue de code

Revue du 22 novembre 2026 : 8 points, tous corrigés.

1. **Rayons presque égaux, retraits proches de l'automatique** (important).
   Les rayons (0,5 ; 0,51 ; 0,52) étaient refusés après 30 s et 850 Mo. Un
   retrait de 0,505 ou 0,51 pour un rayon de 0,5 l'était aussi : la surface
   s'y retournait. Corrigé par le recul ou le refus du § 4.3. Les rayons
   (0,5 ; 0,51 ; 0,52) et (0,5 ; 0,52 ; 0,52) aboutissent en 2,3 s et 106 Mo ;
   les retraits de 0,501 à 0,53 sont refusés aussitôt, et 0,55 aboutit.
2. **Coût de la sphère guide** (important). À 128 carreaux : 860 Mo et 26 s,
   au lieu des 400 Mo et « quelques secondes » annoncés. Les rangées couplaient
   les trois coordonnées. Chaque coordonnée est maintenant tirée vers le pied
   sur la sphère (§ 4.6) : 340 Mo et 7 s, les tests et le tableau du § 6
   inchangés à 10⁻⁵ près.
3. **Retrait proche de l'automatique ignoré.** Un retrait de 0,501 donnait le
   coin sans retrait. Il n'y a plus de fin au croisement à côté d'un retrait
   demandé : un retrait égal à l'automatique est gardé, un retrait à peine plus
   long est refusé (§ 4.3).
4. **Repli sur le raccord.** L'essai de la sphère ne se repliait que sur un
   refus. Toute exception y mène maintenant au raccord, avec le même retour
   arrière.
5. **Options en Python.** `FilletOptions.setbacks` rend un tuple : `append`
   lève une erreur au lieu de changer une copie. Les champs de `GuideSphere`
   sont en lecture seule : on affecte une nouvelle sphère à
   `NSidedOptions.sphere`.
6. **Tests des refus.** Chaque refus de retrait a son message ; un trou non
   rempli est testé. Un test couvre les retraits sur des congés B-spline.
7. **Doc de `shape_tol`.** Le texte décrit le critère extrapolé, en C++ et en
   Python.
8. **Commentaire en double** dans `tests_brep_variable_fillet_solid`.
