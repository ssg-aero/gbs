# Palier 2 — PR 7 : découpe des faces et des formes

| | |
|---|---|
| PR | branche `feat/brep-split` |
| Document de référence | [brep_palier2.md](brep_palier2.md), § 5 |
| Fichiers | `gbs-brep/split.h` (nouveau), `gbs-brep/brep` (agrégateur), `tests/tests_brep_split.cpp` (nouveau), `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Découper une face par des arêtes qui la traversent, puis une forme le long de
sa section par une autre forme ou par une surface. C'est la deuxième étape des
booléens (PR 9) : après la découpe, chaque face est entièrement dedans ou
dehors de l'autre solide.

## 2. API

```cpp
struct FaceCut { EdgeId edge; shared_ptr<Curve<T, 2>> pcurve; };          // pcurve sur la face, même paramètre que l'arête
struct SplitFaceResult { vector<FaceId> faces; vector<EdgeId> dangling; };
auto split_face(m, FaceId f, span<const FaceCut> cuts, T tol = …) -> BuildResult<SplitFaceResult>;

struct SplitResult { vector<pair<FaceId, vector<FaceId>>> faces; vector<EdgeId> dangling; SectionResult section; };
auto split(m, ShapeId shape, ShapeId tool, SectionOptions = {}) -> BuildResult<SplitResult>;
auto split(m, ShapeId shape, shared_ptr<Surface<T, 3>> tool, optional<array<T, 4>> box = {}, SectionOptions = {}) -> BuildResult<SplitResult>;
```

## 3. Choix d'architecture et justification

### 3.1 Découpe d'une face (`split_face`)

1. **Tolérances.** Celle de chaque coupe monte à celle de la face, et celles
   de ses sommets à celle de la coupe (règle sommet ≥ arête ≥ face du palier
   1, vérifiée par `check`).
2. **Arêtes de bord coupées.** Quand une coupe s'arrête à l'intérieur d'une
   arête du bord (projection à la tolérance), cette arête est coupée par
   `split_edge`, et le **sommet de la coupe** remplace le sommet créé. Tous les
   wires qui utilisent l'arête sont mis à jour : la face voisine reçoit les
   mêmes morceaux, et le shell reste fermé.
3. **Demi-arêtes dans l'espace (u, v).**
   - Chaque co-arête du bord donne une demi-arête, dans le sens qui laisse la
     face à gauche. Le sens est déduit de l'aire signée du wire extérieur :
     les faces qui ne suivent pas la convention du palier 1 sont acceptées.
   - Chaque coupe donne deux demi-arêtes, une par sens.
   - Les **nœuds** sont des sommets **à une position** (u, v). Les deux côtés
     d'un seam sont donc deux nœuds du même sommet. On compare les positions
     à 1e-3 près de la diagonale de la boîte : les pcurves interpolées ne se
     rejoignent pas exactement.
4. **Coupes pendantes.** Une coupe dont un bout n'atteint qu'elle-même est
   retirée (itérativement) et rapportée dans `dangling`.
5. **Boucles.**
   - Depuis chaque demi-arête, la suivante est la première demi-arête
     sortante **dans le sens horaire à partir du chemin de retour**. On tourne
     ainsi le plus à gauche possible, et la plus petite région reste à gauche.
   - Les directions sont mesurées par des cordes très courtes (1e-6 du
     paramètre), puis plus longues (1e-2, 0,2) en cas d'égalité. Deux courbes
     tangentes sont ainsi départagées par leur courbure.
6. **Faces.**
   - Les boucles d'aire positive deviennent les contours extérieurs.
   - Chaque boucle d'aire négative (trou d'origine, ou boucle fermée posée à
     l'intérieur) va dans la plus petite boucle extérieure qui la contient.
   - Les nouvelles faces **partagent la surface** de la face d'origine et
     gardent son nom.
   - Elles la remplacent dans ses shells, au même rang et dans le même sens.
     La face d'origine et ses wires sont effacés.
   - Une face que les coupes ne séparent pas est gardée telle quelle.

### 3.2 Découpe d'une forme (`split`)

On calcule la section (PR 6), on groupe ses arêtes par face de la forme, puis
on découpe chaque face. L'outil n'est pas modifié. Un solide découpé par une
surface **reste un solide**, avec ses faces subdivisées : le séparer en deux
demande des faces de fermeture, c'est-à-dire un booléen (PR 9).

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| `split_face(m, f, span<EdgeId>, tol)` | `span<FaceCut>` (arête et pcurve) | la pcurve de la coupe sur la face vient de la section |
| coupe traversant un seam : coupée et décalée d'une période | rien à faire ici | les intersections sont déjà coupées aux seams (PR 3 et 4) |
| — | sens du bord déduit de l'aire signée | faces venues de fichiers |
| — | tolérances des coupes relevées | règle sommet ≥ arête ≥ face |

## 5. Tests (`tests_brep_split`, 4 cas, 42 assertions)

`check()` passe sur toutes les formes découpées. Les volumes sont comparés à
1e-3 près, la précision de `signed_volume` sur des faces découpées.

| Test | Couvre |
|---|---|
| `plate_through_its_hole` | plaque percée coupée par un plan passant par le centre du trou : deux faces sans trou, aire conservée, `check()` |
| `solids_stay_closed` | boîte coupée par un plan oblique : quatre faces latérales coupées, dix faces, shell fermé, volume conservé ; cylindre coupé en travers (face latérale en deux, à travers le seam) et en biais (ellipse coupée au seam) |
| `holes_and_dangling_cuts` | boîte percée par un cylindre fin : dessus et dessous en un disque et une face trouée, huit faces, volume conservé ; coupe partant du bord et s'arrêtant au milieu : pendante, face gardée |
| `by_another_shape` | boîtes 0..2 et 1..3 : les trois faces de la première traversées par la seconde sont coupées, neuf faces, les deux solides valides ; erreur sur une face morte |

Les suites de section, de classification, d'écriture STEP et d'édition passent
aussi.

## 6. Points à valider

1. **Sommet de la coupe** substitué à celui de `split_edge` sur le bord,
   partagé avec la face voisine.
2. **Nœuds = sommet et position (u, v)**, pour les seams.
3. **Plus petite région à gauche** (premier sortant dans le sens horaire), et
   départage des tangences par des cordes plus longues.
4. **Trous** attribués à la plus petite boucle extérieure qui les contient.
5. **Coupes pendantes** rapportées et ignorées.
6. **Solide découpé par une surface = un solide** aux faces subdivisées ; la
   séparation viendra des booléens.
