# Palier 8 — PR 3 : congé tritangent et bilan du palier

| | |
|---|---|
| PR | branche `feat/brep-full-round` (sur master e8c3114, après #249) |
| Document de référence | [brep_palier8_face_blends.md](brep_palier8_face_blends.md), § 3 (questions 8, 9 et 10 validées le 23 nov. 2026) |
| Fichiers | `gbs-brep/full_round.h` (nouveau), `gbs-brep/brep`, `gbs-brep/fillet.h` (arc par un point), `gbs-brep/blend.h` (orientation de la face du congé, message de fin libre), `tests/tests_brep_full_round.cpp` (nouveau), `gbs-occt/tests/tests_brep2_occt.cpp`, `python/gbsBindBrepOps.cpp`, `python/tests/test_brep_ops.py`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier8_face_blends.md` (bilan, § 9), `News.md`, cette note |
| Statut | À valider (§ 9) ; revue de code : § 10 |

## 1. Objectif

Le congé tritangent, comme le « Tritangent Fillet » de CATIA et le « Full round
fillet » de SolidWorks : la face du milieu M, entre deux faces F₁ et F₂, est
remplacée par le congé tangent aux trois.

```cpp
fillet_full_round(m, s, faces1, middle, faces2, options);     // Python : brep.fillet_full_round
full_round_geometry(m, s, face1, middle, face2, options);      // la géométrie seule (EdgeBlend)
```

- **Aucun rayon n'est donné** : c'est celui de la boule qui touche les trois
  faces, et il suit la géométrie le long de M (la moitié de la largeur d'une
  nervure à flancs parallèles).
- Les cas visés : le dessus d'une nervure arrondi, avec ou sans dépouille ; le
  fond d'une rainure (concave, le congé ajoute de la matière) ; le **bord de
  fuite épais d'une aube** entre ses faces B-spline.
- Le bilan du palier 8 (§ 9 de [brep_palier8_face_blends.md](brep_palier8_face_blends.md)).

## 2. Sources

- Vida, Martin et Várady, « A survey of blending methods that use parametric
  surfaces », *Computer-Aided Design* 26(5), 1994, § 5 : congés à trois
  surfaces, boule tangente aux trois.
- L'enveloppe d'une famille de sphères (surface canal) : une boule de rayon
  variable touche son enveloppe le long d'un **cercle caractéristique**, qui
  n'est pas un grand cercle (palier 5, § 3).
- Le comportement de CATIA et de SolidWorks : la face du milieu disparaît, les
  faces de bout coupent le congé, aucun rayon à saisir.
- Le code des paliers 4 et 5 (`skinned`, `free_end`, `blend_on_copy`) et la
  station à rayon inconnu de la corde constante (palier 7, § 4.3).

## 3. API

```cpp
template <std::floating_point T>
auto fillet_full_round(Model<T> &m, SolidId s, std::span<const FaceId> faces1, std::span<const FaceId> middle,
                       std::span<const FaceId> faces2, const FilletOptions<T> &o = {}) -> BuildResult<SolidId>;
// (et la surcharge sur des std::vector)

template <std::floating_point T>
auto full_round_geometry(const Model<T> &m, SolidId s, FaceId face1, FaceId middle, FaceId face2,
                         const FilletOptions<T> &o = {}) -> BuildResult<EdgeBlend<T>>;
```

- La signature de `fillet_full_round` est celle du palier (§ 0) : des ensembles
  de faces de chaque côté et au milieu. **Écart** : cette PR n'en prend qu'une
  par ensemble, et refuse les autres (§ 9, point 4).
- La **même face des deux côtés** est admise (F₁ = F₂ : une aube d'un seul
  carreau, un U balayé d'une seule courbe ; § 4.3 et § 9, point 10).
- `full_round_geometry` est ajoutée sur le modèle de `fillet_geometry` : la
  surface, la ligne des centres (`spine` ; le rayon en t est la distance de la
  ligne des centres à un contact), les contacts sur F₁ et F₂, tout sur le
  paramètre de l'arête M ∩ F₁. Elle sert aux tests et à la galerie.
- Les options sont celles des congés : `tol`, `approx_tol`, `angular_tol`,
  `extension`. La section doit être l'arc de cercle (le palier exclut les
  sections conique et G2 du tritangent) ; les retraits sont refusés.
- Python : `brep.fillet_full_round(m, s, faces1, middle, faces2, options)` et
  `brep.full_round_geometry(m, s, face1, middle, face2, options)`.

## 4. Démarche

### 4.1 Station (§ 3.1 du palier)

Le **guide** est l'arête g = M ∩ F₁. Une station se pose en g(t), dans le plan
normal à g en ce point. Les inconnues sont :

- les paramètres (u, v) des trois contacts ;
- le centre c ;
- le rayon r.

Soit 10 inconnues. Les équations sont :

- c = Sₖ(uₖ, vₖ) ∓ r Nₖ(uₖ, vₖ) pour F₁, M, F₂ (9 équations, la surface décalée
  de r vers la boule) ;
- (c − g(t)) · g′(t) = 0, le centre dans le plan de la station.

Gauss-Newton (`thick_detail::gauss_newton`), depuis la station voisine.

**Départ.** La première station part du **cercle inscrit** dans les tangentes
des trois faces, dans le plan de la station. On prend :

- le point de l'autre arête M ∩ F₂ dans ce plan (Newton sur sa projection) ;
- la largeur w de M entre les deux arêtes ;
- les directions des côtés, qui s'éloignent de M.

Le triangle formé par M et les deux côtés donne un rayon et un centre fermés.
Les contacts sont déplacés d'autant dans le plan tangent de chaque face. Le
même départ sert si une station voisine saute sur une autre racine (des
contacts qui bougent plus que le guide ne peut les porter : le même critère
que `fillet_law`).

**Côté de la boule.** La boule est dans la matière si M est **convexe** le
long de ses deux arêtes avec les côtés : le congé ôte de la matière (une
nervure). Elle est dehors si M est **concave** le long des deux : le fond d'une
rainure, le congé ajoute de la matière. Le côté est lu par `fillet_detail::setup`
sur chacune des deux arêtes. *Écart* : le palier (§ 3.1) ne prévoyait que la
boule dans la matière ; le cas concave ne coûte rien de plus (§ 9, point 3).

**Marche.** Les stations partent du milieu du guide vers ses deux sommets
(32 pas de chaque côté), puis le continuent au-delà, tant qu'une boule existe.
Le guide est prolongé de `extension` r, **plus** la distance dont l'autre arête
M ∩ F₂ dépasse ses bouts (un bout oblique, le guide du côté court : revue,
§ 10, point 4). Au-delà d'un bout où M finit en pointe, il n'y en a plus :
le congé s'arrête à la dernière station résolue.

### 4.2 Section et surface (§ 3.2 du palier)

**Section.** La section est l'arc du **cercle caractéristique**, le cercle qui
passe par les trois contacts. Il va du contact sur F₁ au contact sur F₂ en
passant par le contact P_M.

- À rayon constant, c'est un grand cercle de la boule : l'arc de la boule par
  P_M du palier.
- À rayon variable, une boule touche son enveloppe le long d'un cercle qui
  n'est pas un grand cercle. Ses trois contacts sont sur ce cercle, dont le
  centre est hors de c. *Précision au palier* : le palier disait « l'arc de la
  boule », ce qui n'est exact qu'à rayon constant (§ 9, point 2).

L'arc dépasse 180° quand les flancs se rapprochent en s'éloignant de M. Il fait
exactement 180° entre des flancs parallèles, alors que l'arc le plus court
entre deux points opposés n'est pas défini. Une fonction nouvelle de
`fillet.h`, `arc_through(c, p1, p2, via)`, en tient compte. C'est l'arc
rationnel quadratique en deux segments, de p1 à p2 par le côté de `via`, joint
en son milieu, refusé au-delà de 300°. `Station` gagne un champ facultatif
`via` : `skinned` prend alors `arc_through` au lieu de `arc`. Les congés
d'arête ne posent pas `via` et ne changent pas.

**Surface.** Les stations sont lues en 17 points ; trois cas :

- **cylindre exact** si le guide est droit, le rayon constant, et les centres
  et les contacts (dans l'espace et dans les paramètres de leurs faces) sur des
  droites le long du guide. La nervure droite, avec ou sans dépouille, est dans
  ce cas : le rayon vaut (w/2) cos α / (1 − sin α), et l'arc couvre π − 2α. Les
  faces n'ont pas à être des plans : les murs plats d'un U en B-spline de
  degré 1 donnent le même cylindre. Le guide est jugé droit sur ses points,
  quelle que soit sa courbe (un `LineCurve` ou une B-spline de degré 1 : une
  arête d'extrusion est souvent la seconde) ;
- **cône exact** si les trois faces sont des plans qui se coupent en un point A
  (une nervure qui se rétrécit) : les centres sur une droite issue de A, le
  rayon linéaire, proportionnel à la distance à A. Le cône de sommet A autour
  de cette droite, de demi-angle asin(dr/ds) (revue, § 10, point 7) ;
- sinon, la **peau de Hermite** des sections (`skinned`), sous `approx_tol` et
  `angular_tol`, avec une station là où un contact croise un nœud d'une face sous C³
  et aux nœuds du guide sous C³. C'est aussi le cas du **tore** d'une nervure
  courbe : ses arêtes, issues d'une révolution, sont des courbes sur surface
  dont le paramètre n'est pas l'angle (§ 9, point 12).

**Emboîtement** (§ 3.2 du palier). Entre deux stations voisines, si
|Δr| ≥ |Δc|, une boule contient l'autre : il n'y a pas d'enveloppe. C'est
refusé ; le contrôle est fait sur les 64 intervalles du congé.

### 4.3 Topologie (§ 3.3 du palier) : M réduite à son guide

Le palier décrit le résultat : M ôtée, F₁ et F₂ bornées par les contacts, une
fin libre sur chaque face de bout. C'est exactement ce que fait le congé d'une
arête vive entre F₁ et F₂. La PR obtient donc cette topologie par la machinerie
des congés d'arête, sur une copie :

1. **Copie** du solide.
2. **M réduite à son guide** (`collapse`) :
   - M est ôtée de sa coque ;
   - dans F₂, l'arête M ∩ F₂ est remplacée par le guide (la pcourbe : le
     segment entre les bouts de l'ancienne co-arête) ;
   - les arêtes de M avec ses faces de bout sont réduites au sommet du guide de
     leur côté, et ôtées des contours des faces de bout.
3. **Congé du guide** sur cette copie (`blend_detail::blend`), avec la géométrie
   du § 4.1–4.2 calculée sur le solide donné. Les contacts coupent F₁ et F₂ ;
   chaque face de bout coupe le congé (`free_end`) ; le résultat est contrôlé
   (`check`).
4. La copie réduite est effacée ; seul le résultat reste.

La copie réduite a la **topologie** d'un solide à arête vive, mais pas sa
géométrie : F₂ est à une largeur de M du guide. Aucune opération n'en lit la
géométrie entre F₂ et le guide : la découpe de F₂ suit le contact, qui est loin
du guide, et la fin libre coupe le congé par la face de bout prolongée. Le
contrôle final porte sur le résultat, qui est géométriquement juste.

*Écart* : le palier ne disait pas comment construire cette topologie. Écrire
une topologie à part aurait dupliqué la découpe des faces, les fins libres et
leurs contrôles (§ 9, point 1).

**Même face des deux côtés** (revue, § 10, point 3). Une aube d'un seul carreau
(un `loft` de sections ouvertes, la forme naturelle d'une aube gbs à bord de
fuite épais) a une seule face F qui touche M par ses deux bords. La station
s'en accommode (la même surface pour F₁ et F₂). La topologie non : réduire M
sur la copie mettrait le guide deux fois dans F, une couture que les congés
d'arête refusent. Donc, sur la copie seulement :

1. F est coupée (`split_face`) le long de l'isoparamétrique à mi-chemin entre
   ses deux arêtes avec M, qui doivent être isoparamétriques dans la même
   direction (`split_between`) ; sinon, refus qui dit de couper la face ;
2. le congé se fait entre les deux morceaux, comme ci-dessus ;
3. les deux morceaux sont fusionnés sur le résultat (`merge_faces` de
   `unify.h`), et les arêtes des faces de bout coupées par la découpe sont
   recousues (`join_edges`) : le résultat a les faces du solide donné.

### 4.4 Orientation de la face du congé (`blend.h`)

`blend_face` orientait la face du congé par la **somme des normales des deux
faces**, comparée à la normale du congé au milieu de sa section. Entre deux
flancs parallèles, cette somme est nulle ; au-delà de 180°, elle pointe à
l'envers. Le premier essai donnait une coque non orientable.

L'orientation est maintenant lue **à chaque contact** : la normale du congé y
est celle de la face qu'il touche (le congé lui est tangent). Pour tous les
congés d'arête, les deux règles donnent le même résultat (leur arc fait moins
de 180°) : la suite des tests le confirme (§ 6).

### 4.5 Contacts dans leurs faces

- **Sur F₁ et F₂** : le contrôle des congés d'arête (`check_on_face`). Pour F₂,
  il porte entre les stations des bouts de M ∩ F₂ (les pieds de ses sommets sur
  le guide prolongé), non entre ceux du guide : sur un bout oblique, le milieu
  du guide peut être hors de F₂ (revue, § 10, point 1).
- **Sur M** : en 17 points du guide, le contact P_M doit être dans M. Il ne peut
  en sortir que par une arête de bout, comme un congé qui dépasse un coin. Au
  milieu, il doit être dedans.

### 4.6 Refus (§ 3.4 du palier)

Avec des messages qui disent quoi faire ; le modèle reste inchangé.

`InvalidId` :
- un côté ou le milieu vide, une face morte, qui n'est pas du solide, ou
  donnée à la fois comme côté et comme milieu ;
- une face du milieu sans arête avec un côté (le palier : « M non voisine de
  F₁ ou de F₂ ») ; la même face des deux côtés, qui ne touche pas M par deux
  arêtes.

`InvalidBounds` :
- une section conique ou G2 ; des retraits ; des tolérances non positives.

`NonGeneric` :
- **aucune boule** tangente aux trois faces (Gauss-Newton sans solution, ou
  r ≤ tol) : « fillet its two edges with the sides instead » ;
- un **contact hors de sa face** : une nervure plus basse que le rayon (la
  boule touche les flancs sous leur pied), ou une boule qui sort de M ;
- **boules emboîtées** (§ 4.2) ;
- une **face de bout qui ne coupe pas** le congé : « the face at that end of
  the middle face does not cut the full round » ; une M **trop étroite** à son
  bout pour que le congé y rencontre la face de bout : « the middle face too
  narrow there: end it wider » (revue, § 10, point 6) ;
- une **face de bout tangente** à un côté là où elle le touche au bout de M
  (le bout arrondi d'une nervure) : il y faudrait un quart de sphère
  (revue, § 10, point 2 ; § 9, point 11) ;
- la même face des deux côtés dont les arêtes avec M ne sont pas
  isoparamétriques dans une même direction : « split it between them » ;
- M **convexe** le long d'un côté et **concave** le long de l'autre ;
- non traités : un côté ou un milieu de plusieurs faces, M qui fait le tour
  (le bord d'un disque), M trouée, M qui finit en pointe (ou par une arête
  dégénérée), plusieurs faces à un bout de M (une nervure en proue), M qui
  touche un côté par plusieurs arêtes.

## 5. Ce qui ne change pas

- Les congés d'arête et de faces : leur géométrie est inchangée (`via` n'est
  posé que par le tritangent). L'orientation de la face du congé (§ 4.4) donne
  le même résultat pour eux.
- Les messages d'une fin libre qui ne rencontre pas la face de bout ne changent
  que pour `fillet_full_round`. Ceux de `fillet_faces` entre faces parallèles
  nomment maintenant `fillet_full_round`.

## 6. Tests

`tests_brep_full_round` (16 cas, 0,2 s). Les volumes sont comparés aux formes
fermées, à `tol` = 10⁻⁷ (voir plus bas).

| Cas | Attendu | Mesuré |
|---|---|---|
| nervure droite, w = 0,4 et 0,8 | cylindre de rayon w/2, arc de π ; aire ôtée w²/2 − π w²/8 par unité de longueur ; solide donné intact, rien de laissé dans le modèle | rayon à 10⁻¹², volume à 2,1·10⁻¹⁰ ; 3 ms |
| nervure à dépouille α = 5°, 15°, −10° | cylindre de rayon (w/2) cos α / (1 − sin α), arc π − 2α ; aire ôtée 2 r² (tan β − β), β = π/4 − α/2 | rayon à 10⁻¹², volume à 9,4·10⁻¹¹ ; 3 ms |
| nervure droite, `tol` = 10⁻⁹ | même volume | 9,5·10⁻¹³ |
| fond de rainure (concave) ; rainure à dépouille de 10° | cylindre de rayon w/2, ou (w/2) cos α / (1 − sin α) ; deux coins ajoutés | 7,9·10⁻¹¹ et 4,3·10⁻¹¹ ; 3 ms |
| nervure à bout oblique (30°, 60°) ; `extension` = 0 | la face de bout coupe le demi-cylindre selon une ellipse ; volume des deux coins prismatiques | 6,7·10⁻¹¹ ; 3 ms |
| côtés échangés : nervure de 0,5 à bout à 60° ; nervures de 2 à bouts à 80°, 84° et 88° (dépassement de 11, 19 et 57 r) | même solide quel que soit le côté du guide ; volume en forme fermée | 2,6·10⁻¹¹ ; écart entre les deux ordres < 10⁻¹⁰ |
| nervure qui se rétrécit (0,6 → 0,2, flancs verticaux) | **cône exact** de sommet A, le point commun des trois plans ; centres sur une droite issue de A, rayon proportionnel à la distance à A | écart au cône 2,3·10⁻¹⁶ ; 2 ms |
| nervure courbe (quart de tour, flancs cylindriques) | peau approchée (tore) ; volume de Pappus | déviation 1,4·10⁻⁶ ; volume à 3,9·10⁻⁸ |
| U d'une seule face (B-spline de degré 1) fermé par un dessus plat | même face des deux côtés : coupée sur la copie, refondue ; 4 faces avant et après ; cylindre exact ; rien de laissé dans le modèle | volume à 5,5·10⁻¹¹ |
| bord de fuite d'aube (faces B-spline cubiques, épaisseur 0,08 → 0,044, vrillée, en flèche) | peau sous `approx_tol` ; contacts tangents aux côtés (sinus de l'angle des normales ≤ `angular_tol`) ; rayon lisse, entre 0,02 et 0,06 | rayon 0,044 → 0,024, différences secondes 1,7·10⁻⁵ ; déviation 1,5·10⁻⁷ ; angle 2,9·10⁻⁷ ; 20 ms |
| la même aube d'un seul carreau (intrados et extrados dans une B-spline, bord d'attaque vif à u = 1/2) | 4 faces avant et après ; le volume de l'aube de trois carreaux arrondie | écart 3,9·10⁻⁹ ; 21 ms |
| refus des faces données | `InvalidId` (côté vide, face côté et milieu, même face des deux côtés sans deux arêtes avec M, sans arête avec un côté, face d'un autre solide), `InvalidBounds` (section conique, retraits), `NonGeneric` (plusieurs faces ; convexe et concave ; plusieurs faces à un bout) ; modèle inchangé | messages vérifiés |
| refus de forme | nervure plus basse que le rayon (contact hors de sa face) ; dessus en pointe ; nervure en proue ; bord d'un disque | messages vérifiés, modèle inchangé |
| refus des bouts | bout arrondi tangent aux flancs ; dessus qui s'amincit jusqu'à 0,005 ; dessus troué ; dessus entaillé depuis un flanc (deux arêtes avec lui) ; dessus sur un carreau à arête dégénérée à la pointe | messages vérifiés, modèle inchangé |

**Tolérance des volumes.** L'écart de volume suit `tol` linéairement : 2,7·10⁻⁹
à la tolérance par défaut (`brep_default_tolerance`, 10⁻⁶), 2,1·10⁻¹⁰ à 10⁻⁷,
9,5·10⁻¹³ à 10⁻⁹. Il vient des courbes où les faces de bout coupent le congé,
des B-splines cubiques à `tol`/10 du demi-cercle exact, comme au bout de tout
congé. Les tests de volume passent `tol` = 10⁻⁷ et gardent les 10⁻⁹ du palier
(revue, § 10, point 5). La nervure courbe et l'aube restent à la tolérance par
défaut : sur une peau approchée, la courbe de bout ne peut pas suivre une
tolérance plus fine que la peau.

**OCCT** (`tests_brep2_occt.full_round_ribs`). OCCT n'a pas de congé
tritangent. Le cas équivalent : la nervure abaissée jusqu'aux contacts,
fusionnée (`BRepAlgoAPI_Fuse`) avec le cylindre tangent, ses faces d'un même
domaine fondues. Nervure droite et nervure à dépouille de 10° : même solide,
mêmes faces (10), volume à 1,0·10⁻⁹ et 8,0·10⁻¹⁰. *Écart* : le palier (§ 6)
ne prévoyait que les formes fermées pour le tritangent.

**Python** : `test_fillet_full_round` (nervure droite : volume, cylindre,
rayon de `full_round_geometry` ; dépouille de 10° : rayon ; deux refus, de
`fillet_full_round` et de `full_round_geometry` ; docstring). La galerie gagne `full_round_rib` (dépouille de 8°, volume en forme
fermée) et `full_round_blade` (bord de fuite d'aube), vérifiés par
`test_examples.py` (solide valide, maillage fermé, volume).

## 7. Illustration

Note d'illustration (montrée à la validation, hors du dépôt), rendus réels de
la branche : la nervure droite et à dépouille avant / après, le fond de
rainure, le bout oblique, la nervure qui se rétrécit, le bord de fuite d'aube
avant / après et le rayon le long de l'aube.

## 8. Estimation

| Partie | Heures |
|---|---|
| station à trois faces, section par P_M, cylindre exact, peau | 6 |
| topologie (M réduite au guide), orientation, refus | 5 |
| tests (nervures, rainure, bouts, aube), OCCT | 4 |
| Python, galerie, bilan du palier, illustration | 3 |
| **Total** | **18** |

## 9. Points à valider

1. **Topologie par réduction de M à son guide** (§ 4.3), sur une copie, puis
   la machinerie des congés d'arête (découpe des faces, fins libres, contrôle).
   *Recommandé* : le résultat est celui que décrit le palier, sans dupliquer la
   topologie des paliers 4 à 6.
2. **Section par le cercle caractéristique** des trois contacts (§ 4.2) plutôt
   que « l'arc de la boule » : identique à rayon constant, seul juste à rayon
   variable (palier 5). *Recommandé.*
3. **Côté concave** (fond de rainure) accepté en plus de la boule dans la
   matière ; M convexe d'un côté et concave de l'autre refusée. *Recommandé* :
   c'est le même code, et SolidWorks l'accepte.
4. **Un seul côté et un seul milieu** par appel (pas d'ensembles de faces), et
   pas de M qui fait le tour (bord d'un disque, bord d'attaque d'une aube
   fermée sur elle-même), ni de M qui finit en pointe ou en proue, ni de M
   trouée ou qui touche un côté par plusieurs arêtes : refusés avec « not
   handled ». C'est un écart au palier (§ 3.1, § 3.3). *Recommandé* pour
   cette PR : les cas visés (nervures, bord de fuite épais) y sont tous ; les
   ensembles demandent une chaîne de guides et le bord fermé une boucle sans
   face de bout, à reprendre dans une suite (bilan, § 9 du palier).
5. **Orientation de la face du congé par ses contacts** (§ 4.4), dans
   `blend.h`, pour tous les congés : même résultat sur les congés d'arête,
   seule juste au-delà de 180°. *Recommandé.*
6. **`arc_through`** et `Station::via` dans `fillet.h`, l'arc en deux segments
   par un point, jusqu'à 300°. *Recommandé.*
7. **Plan des stations normal au guide M ∩ F₁** (question 8 validée) : le
   paramètre du congé est celui de cette arête. Depuis la revue (§ 10, point 1),
   le résultat ne dépend plus du côté choisi comme F₁ sur un bout oblique
   (testé, à 10⁻¹⁰) ; il en dépend un peu (l'échantillonnage, pas la surface
   exacte) pour une aube. *Recommandé* tel que validé.
8. **Tolérance des volumes** : 10⁻⁹ dans les tests à `tol` = 10⁻⁷ ; l'écart
   suit `tol` (2,7·10⁻⁹ à la tolérance par défaut, 10⁻⁶). *Recommandé* : c'est
   l'écart des fins libres de tout congé.
9. **`full_round_geometry`** ajoutée à l'API, comme `fillet_geometry`, et la
   **comparaison à OCCT** par la fusion d'une nervure abaissée et d'un
   cylindre, en plus des formes fermées. *Recommandé.*
10. **Même face des deux côtés** (§ 4.3) : coupée sur la copie à mi-chemin
    entre ses deux arêtes avec M (isoparamétriques dans une même direction),
    puis refondue sur le résultat. *Recommandé* : c'est la forme naturelle d'une
    aube gbs à bord de fuite épais, sans refonte des congés d'arête ; le cas non
    isoparamétrique est refusé avec un message qui dit de couper la face.
11. **Bout arrondi refusé** (la face de bout tangente aux côtés : une nervure en
    stade), avec un message qui nomme le cas. CATIA et SolidWorks y font un
    quart de sphère. Le traiter demanderait une fin d'un genre nouveau : le
    congé s'y arrête sur la sphère de rayon r centrée au bout de sa ligne des
    centres, tangente à M et à la face de bout ; les contacts y finissent
    tangents aux arêtes entre la face de bout et les côtés, ce que `free_end`
    refuse aujourd'hui (un contact qui rencontre la face de bout à un sommet).
    Estimation : **10 à 12 h** (la sphère et sa découpe sur le modèle du coin du
    palier 4, la face de bout coupée le long d'un demi-cercle, les tests en
    forme fermée). *Recommandé* : laissé pour une suite, au bilan.
12. **Cône exact, tore approché** (§ 4.2). Le cône d'une nervure qui se rétrécit
    est exact. Le tore d'une nervure courbe reste une peau (déviation 1,4·10⁻⁶,
    volume à 3,9·10⁻⁸ de Pappus) : ses arêtes, issues d'une révolution, sont des
    courbes sur surface dont le paramètre n'est pas l'angle, et le congé est
    paramétré par le guide. Le rendre exact demande de reparamétrer les contacts
    sur l'angle (2 à 3 h). *Recommandé* : au bilan, avec les autres manques.

## 10. Revue de code

Revue de la branche locale (24 nov. 2026), avant le premier push : aucun
bloquant, trois importants, sept mineurs. Tous traités, chacun avec son test
(`tests_brep_full_round`, sauf mention). Les sondes du relecteur ont servi de
point de départ aux tests.

1. **Contact sur F₂ contrôlé au milieu du guide** (important). Sur une nervure
   courte à bout oblique, le guide du côté long avait son milieu au-delà de
   F₂ : refus « a contact runs off its face », alors que l'autre ordre passait.
   Correction : le contrôle porte entre les stations des bouts de M ∩ F₂
   (§ 4.5). Test `sides_swapped` : les deux ordres, même solide, volume en forme
   fermée.
2. **Bout arrondi** (important). Une nervure en stade était refusée avec un
   message faux (« too oblique to it, or narrower than it »). Correction : la
   face de bout tangente à un côté est reconnue d'emblée et refusée avec un
   message qui nomme le cas ; le coût de le traiter est au § 9, point 11. Test
   `refusals_of_ends`.
3. **Même face des deux côtés** (important). Une aube d'un seul carreau était
   refusée (`InvalidId`, « a face given twice »). Correction : admise (§ 4.3,
   § 9, point 10). Tests `one_face_on_both_sides` (U en B-spline de degré 1,
   cylindre exact) et `blade_of_one_patch` (nouvelle aube d'un seul carreau,
   comparée à celle de trois carreaux).
4. **Bout très oblique** (mineur). Refusé au-delà de 8 r de dépassement, la
   limite des reprises des congés et non de la géométrie. Correction : le guide
   est prolongé d'emblée du dépassement mesuré de M ∩ F₂ (§ 4.1). Tests : bouts
   à 80°, 84° et 88° (11, 19 et 57 r), acceptés, volumes à 10⁻¹¹.
5. **Tolérance par défaut** (mineur). La note la disait de 10⁻⁷ ; elle est de
   10⁻⁶. Correction de la note (§ 6), et les tests de volume passent `tol` =
   10⁻⁷ pour garder les 10⁻⁹ du palier.
6. **M qui s'amincit presque à zéro** (mineur). Refus générique de fin libre.
   Correction : un message propre au tritangent (« the middle face too narrow
   there: end it wider »). Test `refusals_of_ends`.
7. **Cône et tore approchés** (mineur). Le cône est maintenant exact (§ 4.2 ;
   test `tapered_rib`, écart au cône 2,3·10⁻¹⁶, au lieu d'une peau à 1,9·10⁻¹²).
   Le tore reste approché, justifié au § 9, point 12, et inscrit au bilan
   (test `curved_rib`, volume de Pappus à 3,9·10⁻⁸).
8. **Tests manquants** (mineur). Ajoutés : M trouée, M qui touche un côté par
   deux arêtes, carreau à arête dégénérée (refusé comme M qui finit en pointe :
   une fin faite d'arêtes dégénérées seules), retraits, côtés échangés,
   nervure courbe, rainure à dépouille, `extension` = 0 ; en Python,
   `full_round_geometry` en refus.
9. **Messages et commentaires** (mineur). Les refus de `fillet_faces` entre
   faces parallèles nomment `fillet_full_round` ; dans `test_brep_ops.py`, la
   face du dessous n'est plus appelée `plate`.
10. **C++23** (mineur). `no_ball` est une lambda `[[noreturn]]`, comme le
    refus de la découpe d'une face unique.

Le relecteur n'a trouvé de défaut ni dans l'orientation des faces de congé
(§ 4.4 : mêmes signes sur 15 points de tous les types de congés) ni dans les
fins libres sur un demi-cercle.

