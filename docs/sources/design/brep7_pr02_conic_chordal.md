# Palier 7 — PR 2 : sections coniques et corde constante

| | |
|---|---|
| PR | branche `feat/brep-conic-sections` |
| Document de référence | [brep_palier7_blend_sections.md](brep_palier7_blend_sections.md), § 4.1 à § 4.3 et § 4.5 |
| Fichiers | `gbs-brep/fillet.h`, `gbs-brep/blend.h`, `tests/tests_brep_blend_sections.cpp` (nouveau), `python/gbsBindBrepOps.cpp`, `python/tests/test_brep_ops.py`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier7_blend_sections.md`, les notes du palier 6 et `brep7_pr01_variable_chamfers.md` (statuts), `News.md`, cette note |
| Statut | § 9 validé (23 nov. 2026) ; revue de code : § 10 à valider |

## 1. Objectif

- **Section conique** : un congé dont la section est la conique tangente aux
  faces, entre les contacts du congé circulaire de même rayon, de paramètre ρ.
  ρ est constant ou suit une loi.
- **Corde constante** : un congé dont les contacts restent à la même distance c
  l'un de l'autre ; le rayon suit l'angle des faces.
- **Coins** : des congés coniques ou à corde constante au même sommet font un
  raccord de sommet (palier 6), et non plus une sphère.

## 2. Sources

- L. Piegl, W. Tiller, *The NURBS Book*, § 7.5 : la conique rationnelle
  quadratique, son poids et son point d'épaule (ρ = w / (1 + w)).
- Le comportement de CATIA (« Edge Fillet », *Conic parameter*, *Chordal
  fillet*) et de SolidWorks (*Conic Rho*). OCCT n'a ni l'un ni l'autre : pas de
  comparaison.

## 3. API

```cpp
enum class SectionKind { Circular, Conic };
template <std::floating_point T> struct BlendSection { SectionKind kind; T rho; std::optional<BlendLaw<T>> rho_law;
    static circular(); static conic(T rho); static conic(const BlendLaw<T> &); };
// FilletOptions<T>::section ; EdgeBlend<T>::circular
fillet_chordal(m, s, edges, c, o);           // la section de o.section
fillet_chordal_geometry(m, s, e, c, o);
```

`fillet` (rayon constant ou loi) et `fillet_geometry` prennent la section dans
les options. En Python : `brep.SectionKind`, `brep.BlendSection`,
`FilletOptions.section`, `EdgeBlend.circular`, `brep.fillet_chordal` et
`brep.fillet_chordal_geometry`.

## 4. Démarche

### 4.1 La conique

À chaque station, le congé circulaire donne le centre, les deux contacts et
l'arc de cercle rationnel P₁, S, P₂, où S est le coin des tangentes aux faces.
La conique garde les mêmes pôles et remplace le poids cos(φ/2) de S par
w = ρ / (1 − ρ). Il n'y a donc pas de géométrie nouvelle : les contacts, la
ligne des centres et la topologie restent ceux du congé circulaire.
`conic_section` refuse un arc de plus de 170°, où S est rejeté trop loin.

### 4.2 La surface

- **Extrusion exacte.** Avec ρ constant et une ligne des centres droite (arête
  droite entre deux plans, ou cylindre parallèle), toutes les sections sont la
  même conique translatée. La surface est son extrusion, rationnelle de degré
  (2, 1). C'est vérifié : à ρ = 0,414 et à angle droit, la surface est à
  10⁻¹² du cylindre du congé circulaire.
- **Ailleurs** : la peau de Hermite du palier 4 (`skinned`), avec la conique à
  la place de l'arc. Le contrôle de la normale compare la surface à celle des
  sections exactes voisines, et non plus au rayon de la boule.
- **Points anguleux d'une loi de ρ.** La vraie surface y fait un pli : les
  sections, les mêmes de part et d'autre, balaient ensuite autrement (jusqu'à
  6,5° pour ρ linéaire 0,3 ; 0,7 ; 0,3 à r = 0,3). Une peau lisse ne le suit
  pas, et la loi était refusée, même pour un coude de 0,01 (revue de code,
  § 10). Pour un congé de rayon constant, le congé est maintenant coupé en
  morceaux aux points anguleux, une peau par morceau. Les morceaux partagent la
  section au point anguleux, un arc commun dans le solide, et la ligne des
  centres et les contacts sont joints C0 à travers eux. Avec une loi de rayon
  (corde constante comprise), une loi de ρ anguleuse est refusée.
- **Une loi de ρ** est donnée au paramètre de l'arête. Pour un congé de rayon
  constant, le paramètre de la surface est celui de la ligne des centres : la
  loi y est ramenée linéairement entre les sommets, où les deux paramètres se
  correspondent. Sans cela, ρ était lu au mauvais endroit, et le test de volume
  l'a montré (0,071 ôté au lieu de 0,047).
- **Refus** : une conique sur une loi de rayon qui change de pente. Les fuseaux
  de la boule, entre les morceaux, ne raccordent que des arcs de cercle.

### 4.3 La corde constante

À chaque paramètre t de l'arête, prolongements compris, le rayon r(t) est celui
de la boule qui, dans le plan normal à l'arête, touche les deux faces en deux
points distants de c.

- Le calcul procède par sécante sur r. Le point de départ est la formule entre
  deux plans, c / (2 sin(φ/2)), où φ est l'angle des normales des faces.
- Pour chaque r, la boule se pose par Gauss-Newton sur les quatre paramètres des
  contacts. Le centre est commun aux deux faces et se trouve dans le plan normal.
- Les boules sont **cheminées** depuis le milieu de l'arête, par pas d'un
  soixante-quatrième. Repartir du point le plus proche faisait sauter r sur une
  autre solution dans le prolongement (0,178 au lieu de 0,117), et la surface
  devenait inapprochable.
- Le rayon est calculé sur l'arête prolongée aussi loin que `fillet_law`
  prolonge le congé : max(`extension`, 1,5) fois le plus grand de c et du
  rayon, un quart de plus pour un rayon qui croîtrait ; au-delà, il est tenu
  constant. Avec `extension` × max(c, r) seulement, toute chaîne d'arêtes
  échouait, le congé étant prolongé de 1,5 r aux sommets de la chaîne (revue de
  code, § 10).

r(t) est alors la loi d'un congé à rayon variable du palier 5 (`fillet_law`,
enveloppe des boules), avec la section des options. Entre deux plans d'angle
constant, le rayon est constant et le congé est un cylindre.

### 4.4 Coins

`EdgeBlend::circular` est faux pour une section conique. La sphère d'un coin
n'est plus tentée dès qu'un congé n'est pas circulaire : le coin devient un
raccord de sommet du palier 6. Ses rubans sont les sections aux retraits, quelle
que soit leur forme.

## 5. Refus

- `InvalidBounds` : ρ hors de ]0, 1[ (constant ou en un point de sa loi), une
  corde qui ne dépasse pas la tolérance.
- `NonGeneric` :
  - une conique sur une loi de rayon qui change de pente, dans une arête ou à
    un sommet d'une chaîne ;
  - une loi de ρ anguleuse avec une loi de rayon (ou une corde constante) ;
  - un arc de plus de 170° ;
  - une boule de la corde qui ne se pose pas sur les deux faces.

## 6. Tests

`tests_brep_blend_sections` (11 cas). Les références sont calculées sans gbs :
l'aire entre la conique et l'arête par la formule de Green sur sa forme
rationnelle, et le volume par Pappus autour d'un bord circulaire.

| Cas | Attendu | Mesuré |
|---|---|---|
| ρ du cercle à angle droit | la surface du cylindre du congé circulaire ; même volume | 10⁻¹² ; 10⁻⁹ (intégration du volume) |
| ρ = 0,3 / 0,5 / 0,7 sur une arête de boîte | extrusion exacte ; volume 3 − 3 A(ρ) r² ; décroissant en ρ | 10⁻⁹ |
| aide de calcul | parabole r²/6, cercle (1 − π/4) r² | 10⁻¹² |
| loi de ρ de 0,3 à 0,7 | surface sous `approx_tol` ; volume 3 − ∫ A(ρ(s)) | 10⁻⁵ |
| loi de ρ linéaire 0,3 ; 0,7 ; 0,3 | deux morceaux sous `approx_tol` ; 8 faces ; volume 3 − ∫ A(ρ(s)) ; refusée avec une loi de rayon | écart 4,7·10⁻⁶ ; 10⁻⁵ |
| bord d'un cylindre, conique 0,6 (ligne des centres circulaire) | Pappus : π − 2π × moment de l'aire | 10⁻⁵ |
| coin de trois congés coniques ; dalle en stade conique (chaîne) | raccord de sommet sans sphère ; solides valides | 10 faces chacun |
| corde à angle droit, à 120° | contacts à c (10⁻⁷) ; le congé de rayon c / (2 sin(φ/2)) | volumes à 10⁻⁷ |
| corde au bord d'un cylindre | le congé du rayon équivalent | 10⁻⁶ |
| corde au bord d'une dalle en stade (chaîne de deux droites et deux demi-cercles) | le congé de rayon c / (2 sin 45°), solide valide | 3·10⁻⁸ |
| corde au bord d'un cylindre coupé par un plan incliné de 30° | contacts à c le long du bord elliptique, le rayon variant de plus de 0,01 | 10⁻⁶ |
| refus | ρ = 0, 1, −0,2, 1,5 ; conique sur une loi linéaire anguleuse, dans une arête et à un sommet de chaîne ; loi de ρ sur une autre arête | codes et messages du § 5, modèle inchangé |

En Python, trois vérifications : ρ du cercle donne le cercle, la parabole ôte
r²/6 et la corde donne le rayon c/√2. La galerie gagne trois cas :
`fillet_conic_box` (ρ = 0,3 et 0,7), `fillet_conic_corner` et
`fillet_chordal_prism`.

## 7. Comparaison à OCCT

Aucune : OCCT n'a ni conique ni corde constante. Le congé circulaire est
inchangé, ce que les tests d'OCCT des paliers 4 à 6 vérifient toujours.

## 8. Limites

- **Fuseaux et coniques.** Une loi de rayon linéaire à point anguleux est
  refusée avec une section conique, dans une arête comme à un sommet de chaîne
  (le fuseau de la boule ne raccorde que des arcs de cercle).
- **Loi de ρ anguleuse** : en morceaux pour un congé de rayon constant,
  refusée avec une loi de rayon.
- **ρ en loi sur un congé de rayon constant** : la correspondance entre le
  paramètre de l'arête et celui de la ligne des centres est linéaire entre les
  sommets. Elle est exacte sur une arête droite et approchée sur une arête
  courbe, ce qui déplace légèrement la loi. Sur une arête fermée, on suppose
  que la boucle part de son sommet.
- **Coût de la corde constante** : chaque rayon coûte une sécante de quelques
  Gauss-Newton à Jacobien numérique, avec un cache. Quelques dixièmes de seconde
  par arête dans les tests.
- **Bord d'un booléen** (hors de cette PR) : sur le cylindre coupé en biais, le
  bord entier n'est pas reconnu comme une chaîne, même avec un congé de rayon
  constant. Ses morceaux ne sont pas tangents à 10⁻⁶ près à leurs sommets. Une
  tâche séparée est proposée.

## 9. Points à valider

1. **Conique sur les pôles de l'arc** : même S, poids ρ / (1 − ρ), contacts du
   congé circulaire, comme validé au § 8 du palier.
2. **Extrusion exacte** autour d'une ligne des centres droite avec ρ constant ;
   peau de Hermite ailleurs, avec le contrôle de normale par sections voisines.
3. **Loi de ρ ramenée linéairement** au paramètre de la ligne des centres pour
   un congé de rayon constant (§ 8).
4. **Refus de la conique sur une loi de rayon anguleuse** (fuseaux circulaires).
5. **Corde constante par sécante et boules cheminées**, puis `fillet_law` :
   l'enveloppe des boules, la section des options possible (corde et conique
   ensemble). Le rayon est connu aussi loin que le congé est prolongé.
6. **Coins non circulaires en raccord de sommet**, la sphère n'étant tentée que
   si toutes les sections sont des arcs de cercle.
7. **Loi de ρ anguleuse en morceaux** pour un congé de rayon constant (un arc
   commun au point anguleux, la surface pliée), refusée avec une loi de rayon
   (revue de code).

## 10. Revue de code

Revue du 23 novembre 2026 : 7 points, tous corrigés.

1. **Corde constante le long d'une chaîne** (important). `fillet_chordal`
   échouait sur toute chaîne d'arêtes tangentes : le rayon de la corde n'était
   connu que sur l'arête prolongée de `extension` × max(c, r), et `fillet_law`
   prolonge le congé de 1,5 r aux sommets d'une chaîne. Le prolongement suit
   maintenant celui de `fillet_law`, et le rayon est tenu au-delà (§ 4.3). Le bord
   de la dalle en stade se construit en 0,1 s, à 3·10⁻⁸ du congé de rayon
   c/√2 (nouveau test).
2. **Loi de ρ anguleuse** (important). Elle était refusée, même pour un coude
   de 0,01 : la vraie surface y fait un pli, qu'une peau lisse ne suit pas. Le
   congé de rayon constant est maintenant coupé en morceaux aux points anguleux ;
   avec une loi de rayon, elle est refusée avec un message clair (§ 4.2).
3. **Message de la loi de ρ** : « an edge not filleted », et non plus
   « chamfered ».
4. **Jonction de congés coniques** : à un sommet de chaîne où la loi de rayon
   change de pente, le refus est explicite (le fuseau de la boule n'y serait pas
   tangent), au lieu de « do not meet along a crease ».
5. **Tests** : la corde le long d'une chaîne, la loi de ρ anguleuse, les messages
   des refus ; un `(void)A;` ôté.
6. **Python** : la docstring de `BlendSection.rho_law` dit qu'il est rendu en
   copie et se règle en affectant une loi.
7. **Note** : la base de la branche (master) et le prolongement de la corde.
