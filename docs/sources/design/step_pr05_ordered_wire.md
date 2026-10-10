# Lecture STEP — PR 5 : wire ordonné, pôles et sens des faces

| | |
|---|---|
| PR | branche `feat/brep-ordered-wire` |
| Document de référence | [step_reader.md](step_reader.md), § 5.1 à 5.3 ; palier 1 : [brep_pr03_builders_wire.md](brep_pr03_builders_wire.md), [brep_pr05_face_wire.md](brep_pr05_face_wire.md) |
| Fichiers | `gbs-brep/builders.h`, `python/gbsBindBrep.cpp`, `tests/tests_brep_ordered_wire.cpp` (nouveau), `python/tests/test_brep.py`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Compléter les builders du palier 1 pour accepter la topologie **telle que STEP
l'écrit**, sans sewing :

1. une boucle donnée **dans l'ordre**, où l'arête de seam figure deux fois ;
2. des faces sur une sphère ou un cône dont les **pôles** n'ont pas d'arête
   (STEP n'a pas d'arête dégénérée) ;
3. le **sens** des faces : garder le sens des co-arêtes du fichier.

## 2. API

```cpp
struct OrientedEdge { EdgeId edge; Orientation orient; };

auto make_wire_ordered(m, std::span<const OrientedEdge> coedges, tol) -> BuildResult<WireId>;
auto make_face_use(m, surface, outer, holes, MakeFaceOptions) -> BuildResult<FaceUse>;
```

`make_face(m, surface, outer, holes, opts)` garde sa signature. Il gagne les
traitements du seam et des pôles (§ 3.2 et 3.3), dont `make_face_use`
partage le cœur.

Python : `brep.make_wire_ordered(m, [(edge, Orientation), …])`,
`brep.make_face_use(m, surface, outer, holes, options)`.

## 3. Choix d'architecture et justification

### 3.1 `make_wire_ordered` : l'ordre du fichier, un seam permis

`make_wire` (palier 1) réordonne des arêtes données en vrac et refuse une
arête répétée. Un `edge_loop` STEP est déjà ordonné et orienté, et une face
cylindrique complète s'écrit `[cercle bas, seam, cercle haut inversé, seam
inversé]`. Le nouveau builder :

- garde l'ordre et les sens donnés ;
- vérifie que chaque co-arête finit où commence la suivante : même sommet, ou
  sommets à moins de `max(tol, tol_a + tol_b)`, alors fusionnés comme dans
  `make_wire`. La fusion est factorisée dans `detail::merge_vertices`,
  partagée par les deux builders ;
- accepte qu'une arête serve **deux fois en sens opposés** (seam). Trois
  usages, ou deux dans le même sens, sont refusés (`DuplicateEdge`) ;
- refuse les arêtes dégénérées, insérées par la face (§ 3.3) ;
- rend un wire fermé si la dernière co-arête revient au départ, ouvert sinon.

### 3.2 Côté du seam : par les voisines, puis à l'opposé de l'autre usage

Une co-arête posée **entièrement sur le seam** d'une direction fermée
(u = u1 ou u = u2) peut prendre l'un ou l'autre côté. Au palier 1, le côté
suivait la fin de la co-arête précédente. Cela ne suffit plus quand les deux
usages d'une arête sont dans la face, ni quand la boucle commence par le seam.

Le builder procède désormais en trois temps :

1. **Calcul de toutes les pcurves**, sans indication de côté.
2. **Côté par les voisines** : une co-arête ambiguë dans la direction k prend
   le côté de la voisine (précédente ou suivante) non ambiguë qu'elle touche,
   sauf si la jonction est un pôle (k y est libre). On recommence jusqu'à ce
   que plus rien ne change. Cela suffit pour le cylindre, quel que soit le
   point de départ de la boucle, et pour le tore (deux seams : chaque
   co-arête n'est ambiguë que dans une direction).
3. **Repli** : une co-arête restée ambiguë (sphère bornée par son seam seul :
   ses seules voisines sont des pôles) va **à l'opposé de l'autre usage** de
   son arête, sinon du côté u1.

Ce choix est fait direction par direction : sur un tore, l'arête du petit
cercle est ambiguë en u et pas en v.

### 3.3 Pôles : arêtes dégénérées insérées

Après le choix des côtés, deux co-arêtes consécutives peuvent se rejoindre en
3D (même sommet) mais pas en (u, v). C'est le cas d'une sphère : le seam monte
à u = 2π, redescend à u = 0, et le pôle nord va de (2π, π/2) à (0, π/2). Le
builder insère une **co-arête dégénérée** (`make_degenerate_edge`) si les deux
conditions suivantes sont réunies :

- les deux extrémités ne diffèrent que dans une direction k ;
- la surface est singulière dans cette direction aux deux extrémités : ∂S/∂k
  nul, le critère `free_coordinate` du palier 1.

Sa pcurve est le segment entre les deux extrémités, paramétré par la
coordonnée k, comme les côtés dégénérés des faces naturelles. Ailleurs, l'écart
n'est pas comblé et le contrôle de continuité existant le refuse
(`CrossesSeam`).

Les arêtes dégénérées sont créées pendant la construction et **effacées si elle
échoue** : le modèle est rendu inchangé, comme pour la face depuis des boucles
2D du palier 1, aux pierres tombales près.

### 3.4 Sens des faces : `make_face_use`

STEP oriente la boucle extérieure dans le sens direct **vu depuis la normale de
la face**. Le builder la tourne dans le sens direct **en (u, v)**, c'est-à-dire
autour de la normale de la surface. `make_face_use` rend le `FaceUse` qui
compense :

- `Forward` si la boucle était déjà dans le sens direct en (u, v) : normale de
  la face égale à celle de la surface, `same_sense = .T.` ;
- `Reversed` si le builder l'a retournée.

Composé avec l'usage, chaque co-arête retrouve le sens du fichier. Les arêtes
partagées restent donc utilisées en sens opposés par les faces voisines, et le
shell est orienté vers l'extérieur **sans retournement**. Le test de la boîte
le vérifie : plans dont trois normales pointent vers l'intérieur, volume +1
directement.

**Cas ambigu.** Une face bornée **uniquement par des seams et des pôles**
(sphère, tore complets) peut être orientée dans les deux sens : rien ne fixe le
côté du premier seam. Une telle face n'a aucune arête commune avec une autre
face : elle forme son propre shell. `make_solid` corrige alors son sens par le
volume. Le lecteur (PR 6) pourra aussi s'appuyer sur `same_sense`.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| pcurves de seam départagées par les pcurves du fichier (§ 4.2, point 4) | départagées par les co-arêtes voisines, puis à l'opposé de l'autre usage | suffit sur les cas du plan, sans dépendre des pcurves du fichier (souvent absentes) |
| `vertex_loop` → arête dégénérée (§ 5.2) | non traité ici | une `vertex_loop` borde une face sans seam (calotte) : c'est le cas du point suivant |
| — | **limite** : une boucle qui fait le tour d'une direction fermée **sans arête de seam** (bande de cylindre bornée par deux cercles, calotte bornée par un cercle) reste refusée (`CrossesSeam`) | il faut y **insérer** un seam ; OCCT écrit toujours le seam, d'autres systèmes non. Proposition : l'ajouter à la PR 7 (découpe au seam), qui devient « seam : découpe et insertion » |

## 5. Tests (`tests_brep_ordered_wire`, 7 cas, 85 assertions)

| Test | Couvre |
|---|---|
| `make_wire_ordered` | ordre et sens gardés, wire ouvert, seam (arête deux fois), erreurs (chaîne rompue, même sens deux fois, trois usages, vide, identifiant mort, arête dégénérée) sans modifier le modèle, fusion des sommets voisins |
| `cylinder_with_seam` | NURBS exacte de cylindre (PR 3), boucle dans le sens du fichier, partant du seam, ou retournée (`Reversed`) ; seams aux deux bords ; aire 2πh ; fermé par deux disques : volume πR²h sans retournement, `check()` valide |
| `sphere_poles` | sphère bornée par son seam seul : deux arêtes dégénérées insérées, aire 2π², solide de volume 4/3 πR³ valide |
| `cone_apex` | cône jusqu'à son sommet, une arête dégénérée, solide fermé par un disque : volume πR²H/3 |
| `torus_two_seams` | tore complet borné par ses deux seams : côtés déduits des voisines, aucune arête dégénérée, volume 2π²Rr² |
| `box_face_senses` | boîte lue comme en STEP (sommets et arêtes partagés, boucles directes vues de l'extérieur), plans de normales mixtes : `FaceUse` attendus, volume +1 sans retournement, solide valide |
| `failures_leave_the_model_unchanged` | bande sans seam refusée (`CrossesSeam`) ; échec après insertion des pôles : arêtes dégénérées effacées, wires intacts |

Toutes les suites BREP passent (palier 1 inchangé), ainsi que les tests Python
du BREP (`test_ordered_wire_and_face_use`).

## 6. Points à valider

1. **`make_wire_ordered`** séparé de `make_wire`, seam permis (deux usages en
   sens opposés), fusion des sommets factorisée.
2. **Côté du seam** par les voisines, puis à l'opposé de l'autre usage ; ce
   comportement s'applique aussi à `make_face` du palier 1.
3. **Arêtes dégénérées insérées** aux pôles par le builder de face, effacées en
   cas d'échec.
4. **`make_face_use`** rend le `FaceUse` qui garde le sens du fichier ; le cas
   ambigu (face bornée par ses seams seuls) est corrigé par `make_solid`.
5. **Insertion d'un seam** (bandes et calottes sans arête de seam) ajoutée à la
   PR 7.
