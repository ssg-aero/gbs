# Palier 4 — PR 6 : congés et chanfreins, topologie

| | |
|---|---|
| PR | branche `feat/brep-fillet-topology` |
| Document de référence | [brep_palier4.md](brep_palier4.md), § 4.2 et § 4.5 |
| Fichiers | `gbs-brep/blend.h` (nouveau), `gbs-brep/fillet.h`, `gbs-brep/brep`, `python/gbsBindBrepOps.cpp`, `tests/tests_brep_fillet.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/examples/gbs_examples.py`, `docs/sources/examples/gallery.py`, `python/tests/test_brep_ops.py`, `python/tests/test_examples.py`, `docs/sources/design/brep_palier4.md`, `docs/sources/design/brep4_pr05_fillet_geometry.md`, `News.md`, cette note |
| Statut | À valider (§ 9 ; revue de code : § 10) |

## 1. Objectif

Les congés et les chanfreins d'arêtes d'un solide : un nouveau solide, dont
les faces de chaque arête sont coupées le long des contacts de la PR 5, les
faces aux extrémités coupées par le congé, et la face du congé posée entre
elles. Le solide donné n'est pas modifié.

Sont couverts :

- une ou plusieurs arêtes indépendantes, convexes ou concaves ;
- les arêtes fermées (bords d'un cylindre, pied d'un bossage) ;
- les chaînes d'arêtes tangentes, ouvertes ou fermées (le tour d'un oblong) ;
- les congés qui s'annulent là où les faces deviennent tangentes (té de même
  rayon) ;
- les congés approchés (B-splines) comme les exacts.

Les raccords de sommet (trois arêtes congées en un coin, sphère de coin) sont
la PR 7. Deux arêtes congées qui se rencontrent ailleurs que le long d'une
chaîne tangente sont refusées.

## 2. Sources

- Rossignac et Requicha, « Constant-radius blending in solid modelling »
  (*Computers in Mechanical Engineering*, 1984) : le congé comme surface entre
  deux contacts, la coupe des faces adjacentes.
- Vida, Martin et Várady, « A survey of blending methods that use parametric
  surfaces » (*Computer-Aided Design* 26(5), 1994), § 5 : les extrémités d'un
  congé (coupe par la face d'extrémité, chaînes, sommets).
- Le comportement de CATIA et SolidWorks pour l'usage : chaîne tangente
  congée d'un seul tenant, face d'extrémité coupée (convexe) ou prolongée
  (concave) par le congé. OCCT (`BRepFilletAPI_MakeFillet`,
  `BRepFilletAPI_MakeChamfer`) sert de comparaison, en boîte noire.

## 3. API

```cpp
auto fillet(Model<T> &, SolidId, std::span<const EdgeId> edges, T r, const FilletOptions<T> & = {}) -> BuildResult<SolidId>;
auto chamfer(Model<T> &, SolidId, std::span<const EdgeId> edges, T d, const FilletOptions<T> & = {}) -> BuildResult<SolidId>;
auto chamfer(Model<T> &, SolidId, std::span<const EdgeId> edges, FaceId f, T d1, T d2, const FilletOptions<T> & = {}) -> BuildResult<SolidId>;
```

Les options sont celles de la PR 5 (`tol`, `approx_tol`, `angular_tol`,
`extension`, `intersect`). Pour un chanfrein dissymétrique, d₁ est pris sur la
face `f`, qui doit border chaque arête ; le plan (§ 4.1) n'avait pas de face,
la PR 5 l'a ajoutée à `chamfer_geometry`. Python : `brep.fillet`,
`brep.chamfer` (deux formes).

Résultat : un nouveau solide. En cas d'échec, le modèle est laissé tel qu'il
était, à des entités mortes près (`rollback`).

## 4. Démarche

### 4.1 Copie, géométrie

Le solide est copié (`ShapeCopier`, qui donne aussi l'arête copiée de chaque
arête donnée). La géométrie de chaque arête est celle de la PR 5, calculée sur
la copie : surface, contacts, ligne des centres, prolongés au-delà des
sommets. La copie est ensuite modifiée en place et devient le résultat.

### 4.2 Remplacer des suites de co-arêtes plutôt que découper

Le plan (§ 4.5) prévoyait de découper F1 et F2 le long des contacts
(`split_face`) et de retirer la bande le long de l'arête ; la face
d'extrémité F3 aurait été découpée de même, son coin retiré.

Pour une arête **concave**, cela ne marche pas : la face d'extrémité ne perd
pas son coin, elle le **gagne** (le congé ajoute de la matière dans l'angle
rentrant ; la face avant d'un L s'étend jusqu'à l'arc).

La PR remplace donc, dans le fil de chaque face, la suite de co-arêtes entre
deux points de son bord par une nouvelle suite :

- sur F1 et F2, la suite entre les bouts du contact qui passe par l'arête
  congée devient le contact ;
- sur F3, la suite qui passe par le sommet (les deux morceaux d'arêtes qui y
  arrivent) devient la courbe d'extrémité.

Convexe ou concave, c'est la même opération : la face perd la bande ou le
coin, ou gagne le coin. Il n'y a pas de morceaux de faces à reconnaître après
coup.

La suite retenue est celle qui contient une arête condamnée (l'arête congée
pour F1 et F2, les morceaux d'arêtes du sommet pour F3), le reste du fil n'en
contenant aucune. Les deux sens de parcours sont essayés. Les pcurves de la
nouvelle suite sont décalées de périodes pour se raccorder en (u, v), et leurs
bouts sont contrôlés : ils tombent à 10⁻⁴ de la taille de la face près sur
ceux de la suite remplacée. Cela départage les deux côtés d'une couture.

### 4.3 Extrémités

À chaque sommet d'une arête ouverte, le congé finit de l'une de quatre façons.

**Libre** : sommet de trois arêtes, aucune autre arête congée.

1. La troisième face F3 est celle qui borde les deux autres arêtes du sommet
   (une sur F1, une sur F2). Sa surface est prolongée : pour une arête
   concave, le coin gagné est hors de la face. Il reste dans la boîte (u, v)
   de F3, entre le sommet et les deux points où les contacts coupent ses
   arêtes, si bien que la surface d'origine de F3, même une B-spline bornée à
   sa face, le porte.
2. Chaque contact rencontre F3 par Gauss-Newton sur (t, u₃, v₃) :
   C(t) = S₃(u₃, v₃). Le point doit être sur l'arête entre la face du contact
   et F3, à l'intérieur ; l'arête y est coupée (`split_at_vertex`), au sommet
   qui sera le bout du contact.
3. La **courbe d'extrémité** est l'intersection du congé et de F3 entre ces
   deux points. Elle est cheminée en travers du congé : en s ∈ [0, 1], du
   contact 0 au contact 1, le point du congé à l'abscisse s de l'arc
   (interpolation en (u, v) entre les deux contacts au même t) qui est sur F3
   (Gauss-Newton sur (t, u₃, v₃)). Les points sont interpolés (cubiques) en 3D,
   en (u, v) sur le congé et en (u, v) sur F3. Ils sont doublés de 8 à 512
   jusqu'à ce que les trois courbes soient à `tol` / 10 des points exacts
   entre eux ; les points des milieux, calculés pour ce contrôle, sont les
   nouveaux points.
4. Si un contact ne rencontre pas F3 dans le prolongement (arête oblique à F3),
   la géométrie est recalculée avec un prolongement 2, 4 puis 8 fois plus
   grand (à partir d'un rayon ou d'une distance au moins, `extension` pouvant
   être nul), tant qu'elle se laisse calculer : un échec, quel qu'il soit,
   garde la dernière géométrie.
5. Un contact qui tombe sur un sommet de l'arête entre sa face et F3 (deux
   congés qui s'y rejoignent, un arrondi complet sur une face de largeur 2r)
   est refusé comme tel, et non comme un congé plus large que la face.

**Raccord** : une autre arête congée arrive au sommet. Elle doit prolonger
l'arête tangentiellement, et les deux congés doivent finir sur le même arc :
contacts et milieu de l'arc à 10 `tol` près, plus `approx_tol` par congé
approché. Sinon c'est un coin (PR 7), refusé. Les deux faces de congé
partagent l'arc en travers au sommet. Quand le contact passe d'une face à la
suivante (oblong : plan puis demi-cylindre), son bout doit tomber sur l'arête
entre elles, qui y est coupée.

**En pointe** : les faces sont tangentes au sommet (PR 5 : les contacts se
rejoignent, le congé s'arrête là). Les deux contacts finissent au sommet, et
la face du congé y finit en pointe : une arête dégénérée y porte l'arc
réduit à un point, si ses deux bouts diffèrent en (u, v).

**Boucle** : une arête fermée. Les contacts sont des boucles, la face du
congé fait le tour avec une couture, l'arc en travers au départ, pris deux
fois (pcurves en a et en b).

### 4.4 Faces des arêtes

Les contacts d'une même face sont **chaînés** par leurs sommets communs (le
long d'une chaîne d'arêtes sur la même face). Chaque chaîne est ensuite coupée
là où elle touche le bord de la face : un sommet du fil, ou un point à moins
de 100 `tol` d'une arête non congée. C'est le cas d'une couture traversée (le
té de même rayon : le contact sur le cylindre debout passe sa couture) ou de
l'arête entre deux faces d'une chaîne.

- Chaque morceau va d'un point du bord à un autre ; ses bouts sont posés sur
  le bord (arête coupée au besoin) et il remplace la suite de co-arêtes qui
  longe ses arêtes congées.
- Une chaîne fermée qui ne touche pas le bord (le disque du dessus d'un
  cylindre, le dessus de l'oblong) remplace un fil entier fait de ses arêtes
  congées, dans le même sens de rotation en (u, v).

### 4.5 Face du congé

La face est bâtie sur la surface du congé avec un fil : contact 0, extrémité
en b, contact 1 à rebours, extrémité en a. Les pcurves sont celles des
contacts sur le congé (PR 5) et celles des extrémités. Le fil est tourné dans
le sens trigonométrique en (u, v). La face est orientée comme les faces de
l'arête : sa normale, au milieu, du côté de la somme des normales extérieures
de F1 et F2. Cela vaut pour un congé convexe comme concave, et pour un
chanfrein.

### 4.6 Nettoyage, contrôle

Les arêtes des suites remplacées qui ne servent plus sont effacées, avec leurs
sommets sans usage. Les tolérances des nouvelles arêtes sont relevées aux
écarts mesurés entre la courbe 3D et ses pcurves, et celles des sommets
couvrent les bouts. Le solide est contrôlé (`check`) ; un défaut est un refus.

### 4.7 Cônes des chanfreins

Comparée à OCCT, la PR a fait apparaître un défaut. Le chanfrein d'un bord de
cylindre donnait un cône d'angle négatif (rayon décroissant le long de l'axe).
gbs le gère, mais STEP veut un demi-angle positif, et OCCT le relisait sans
erreur avec un volume faux (6,091 au lieu de 6,216). Le cône est désormais
bâti depuis son petit cercle, l'axe retourné au besoin (`fillet.h`).

L'écrivain STEP (`gbs-io/step/geometry_writer.h`) écrit de plus tout cône
d'angle négatif, quelle que soit son origine, comme le même cône d'angle
positif : l'axe z retourné, et y avec lui, si bien que (u, v) devient
(−u, −v) ; les pcurves de ses faces sont écrites ainsi transformées. Le test
`tests_step_write.cone_of_negative_angle` relit le cône à l'angle positif, ses
arêtes et sa face sur le cône écrit.

## 5. Refus

`NonGeneric`, le modèle laissé tel qu'il était :

| Cas | Message |
|---|---|
| deux arêtes congées qui se rencontrent sans être tangentes, ou plus de deux à un sommet | vertex blends: PR 7 |
| une partie seulement d'une chaîne tangente | a chain of tangent edges blended in part |
| extrémité libre à un sommet de plus de trois arêtes | a blend ending at a vertex of more than three edges |
| contact qui ne rencontre pas F3, même prolongé 8 fois | a contact does not meet the face at an end |
| contact qui rencontre F3 au-delà de l'arête entre elles | a blend wider than a face at an end of its edge |
| contact qui rencontre F3 en un sommet de l'arête entre elles (arrondi complet) | ... at a vertex of it (blends meeting there, as a full round: not handled) |
| contact qui finit dans une face | a contact ends inside a face |
| résultat invalide | the result is not a valid solid |

S'y ajoutent les refus de la PR 5 (rayon au-delà d'une courbure, contact hors
de sa face, etc.), et `InvalidId` / `EmptyInput` pour une arête donnée deux
fois, une liste vide, une arête ou une face d'un autre solide.

## 6. Python, galerie

`brep.fillet` et `brep.chamfer` (deux formes). La galerie a une section
« Fillets and chamfers » de 6 cas : bloc en L (deux arêtes convexes, puis la
concave), cylindre (bord congé puis pied chanfreiné), bossage sur une plaque
(pied concave et bord convexe), tour d'un oblong (chaîne de quatre arêtes),
cylindres croisés (B-spline en boucle), té de même rayon (congé en pointe aux
deux bouts). Chaque cas prend de 0,01 à 0,04 s.

## 7. Tests

`tests/tests_brep_fillet.cpp`, 9 cas de plus (21 en tout). Les volumes exacts
sont comparés à l'intégrale de `signed_volume` (512 points). Celle-ci a déjà
4·10⁻⁷ d'erreur sur un cylindre, d'où une tolérance de 10⁻⁶ pour les cas
courbes.

| Test | Contenu |
|---|---|
| `fillet_box_edges` | arête de boîte (r = 0,1 et 0,3), volume 2 − 2(1 − π/4) r² à 10⁻⁹ ; deux arêtes parallèles ; les quatre arêtes verticales ; le solide donné inchangé |
| `fillet_concave_edge` | arête rentrante d'un L : 3 + (1 − π/4) r² |
| `fillet_closed_edges` | bords d'un cylindre, l'un puis les deux : Pappus, π R² h − 2π(R − x̄)A par bord |
| `fillet_chain` | tour d'un oblong, quatre arêtes en chaîne fermée : V₀ − A(2L + 2π(R − x̄)) ; la chaîne prise en partie refusée |
| `fillet_boss_and_corners` | bossage : pied concave (+2π(R + x̄)A) et bord convexe ; prisme triangulaire, extrémités coupées à 60° : V₀ − A(c − 2x̄/√3) |
| `chamfer_edges` | chanfrein de boîte, égal (d²/2) et dissymétrique (d₁d₂/2) |
| `chamfer_closed_and_concave` | bords d'un cylindre (cônes) : 2π(R − d/3) d²/2 par bord ; arête rentrante d'un L |
| `fillet_approximated` | cylindres croisés (congé et chanfrein B-splines en boucle), té presque tangent (0,99), prisme d'une bosse B-spline (extrémités libres sur des plans, congé et chanfrein) |
| `fillet_vanishing` | té de même rayon : la demi-ellipse (deux arêtes, coupées à la couture du cylindre debout) congée, en pointe aux deux bouts |
| `blend_refusals` | coin (PR 7), liste vide, arête donnée deux fois, rayon trop grand ; modèle inchangé |
| `blend_refusals_more` | extrémité libre à un sommet de quatre arêtes (dessus d'une boîte coupé en diagonale), arrondi complet (deux congés de r = largeur / 2), arête d'un autre solide, face morte, face qui ne borde pas toutes les arêtes |
| `fillet_solid_with_a_void` | boîte creuse (`make_thick`) : arête de la coque extérieure (V₀ − 2A) et arête de la cavité, concave vue du solide (V₀ + 1,8 A), volume des deux coques |

`gbs-occt/tests/tests_brep2_occt.cpp`, `fillet_solids` : chaque cas est
comparé à `BRepFilletAPI_MakeFillet` / `MakeChamfer` sur le même solide (celui
de gbs relu par OCCT). Le résultat de gbs, relu par OCCT, est valide
(`BRepCheck_Analyzer`), a le même nombre de faces et le même volume : à 10⁻⁶
près en relatif pour le L (deux arêtes convexes, l'arête concave, chanfreins),
les bords de cylindre (congés et chanfreins), le bossage et l'oblong ; à 10⁻⁵
pour les cylindres croisés (congé et chanfrein approchés), le té de rayons 0,99
et 1 (16,0996 des deux côtés ; OCCT coupe son congé en plus de faces, 8
contre 6, d'où un nombre de faces comparé au plus) et la bosse B-spline (congé
et chanfrein). Le té de même rayon, où le congé s'annule, n'est pas comparé :
le résultat d'OCCT y est invalide (`BRepCheck_Analyzer`), de volume 16,98
contre 16,19.

Python : `test_fillet_and_chamfer` (`test_brep_ops.py`), `test_fillets` sur
les 6 cas de la galerie (`test_examples.py`, volume exact quand il est connu).

## 8. Limites

- Une arête presque rasante à la face d'extrémité (19° pour la bosse de la PR
  5, dont la tangente au bout est presque dans le plan d'extrémité) demande un
  prolongement d'environ 3 r. Au-delà de 2 r, le prolongement de la B-spline
  du support ne se laisse plus approcher, et c'est refusé. Les tests utilisent
  une bosse qui part d'équerre.
- La courbe d'extrémité est toujours une B-spline interpolée à `tol` / 10,
  même quand elle est un arc exact (plan perpendiculaire à un congé
  cylindrique).
- Entre deux congés approchés qui se raccordent, l'écart (jusqu'à
  `approx_tol`) est absorbé dans les tolérances de l'arc commun et de ses
  sommets.

## 9. Points à valider

1. **API du chanfrein** : une forme à distance égale, et une forme d₁/d₂ avec
   une face de référence commune à toutes les arêtes. OCCT prend une face par
   arête ; on pourra ajouter une liste de paires (arête, face) si le besoin
   vient.
2. **Remplacement de suites de co-arêtes** au lieu de `split_face` (§ 4.2) :
   une seule opération pour le convexe et le concave, et la face d'extrémité
   qui gagne le coin d'une arête concave.
3. **Arêtes congées qui se rencontrent** : seulement le long d'une chaîne
   tangente. Une chaîne prise en partie est refusée, et non propagée
   d'elle-même comme le font CATIA et SolidWorks ; la propagation pourrait
   venir en option avec la PR 7.
4. **Extrémité libre** : seulement à un sommet de trois arêtes. Un sommet de
   plus de trois arêtes (pyramide, bout d'arête sur un sommet de plusieurs
   faces) est refusé.
5. **Prolongement agrandi d'office** (× 2, 4, 8) quand un contact n'atteint
   pas la face d'extrémité, et la limite qui en découle pour les arêtes
   rasantes (§ 8).
6. **Courbe d'extrémité interpolée** à `tol` / 10, sans reconnaissance des
   arcs exacts.
7. **Congé en pointe** : face à trois côtés, l'arc réduit à un point porté par
   une arête dégénérée.
8. **Tolérance des volumes courbes à 10⁻⁶** dans les tests, l'intégrale de
   `signed_volume` étant elle-même à 4·10⁻⁷ sur un cylindre.
9. **Cônes d'angle positif** (§ 4.7) : ceux des chanfreins bâtis ainsi, et
   l'écrivain STEP qui écrit tout cône d'angle négatif comme le même cône
   d'angle positif, ses pcurves transformées.

## 10. Revue de code

Onze points : dix corrigés dans cette PR, un sans objet.

| # | Point | Suite |
|---|---|---|
| 1 | Une arête concave dont la face d'extrémité est une B-spline bornée à sa face : la courbe d'extrémité, hors de la face, serait hors du domaine de sa surface | Sans objet : le coin gagné reste dans la boîte (u, v) de F3, entre le sommet et les points où les contacts coupent ses arêtes (§ 4.3) |
| 2 | L'écrivain STEP écrivait tel quel un cône d'angle négatif, relu par OCCT comme un autre cône ; seuls les chanfreins étaient corrigés | Tout cône d'angle négatif écrit comme le même cône d'angle positif, ses pcurves transformées (§ 4.7) ; test `cone_of_negative_angle` |
| 3 | Les essais de prolongement ne rattrapaient que les erreurs de construction : une évaluation hors bornes à 8 r faisait échouer tout le congé | Tout échec d'un essai garde la géométrie précédente |
| 4 | Les volumes des congés approchés étaient contrôlés plus largement que l'effet mesuré (0,05 pour un congé de 0,013) | Comparaison à OCCT ajoutée, à 10⁻⁵ près, pour le té de rayons 0,99 et 1 et la bosse B-spline (congé et chanfrein) ; le té de même rayon reste sans référence, OCCT y rendant un solide invalide (§ 7) |
| 5 | Plusieurs refus, et les solides à cavité, sans test | `blend_refusals_more`, `fillet_solid_with_a_void` |
| 6 | Les essais multipliaient `extension`, inutiles pour une extension nulle | À partir d'un rayon ou d'une distance au moins |
| 7 | Deux congés qui se rejoignent sur une face (arrondi complet) refusés comme « plus larges que la face » | Refusés comme tels, message distinct |
| 8 | La courbe d'extrémité recalculée en entier à chaque doublement | Les points des milieux du contrôle deviennent les nouveaux points |
| 9 | `SolidView` parcourait toutes les faces à chaque requête | Orientation et coque de chaque face retenues une fois ; arêtes par sommet, faces par arête et pcurves indexées, l'index refait après chaque découpe |
| 10 | Le point le plus proche sur une arête écrit trois fois (et une quatrième dans `fillet.h`) | `fillet_detail::edge_point` |
| 11 | La création d'une arête suivie de `settle` répétée cinq fois | `add_edge` |
