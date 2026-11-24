# Palier 8 — PR 1 : congé entre deux faces

| | |
|---|---|
| PR | branche `feat/brep-face-blends` (sur master) |
| Document de référence | [brep_palier8_face_blends.md](brep_palier8_face_blends.md), § 2.1 à § 2.4 |
| Fichiers | `gbs-brep/face_blend.h` (nouveau), `gbs-brep/brep`, `tests/tests_brep_face_blends.cpp` (nouveau), `gbs-occt/tests/tests_brep2_occt.cpp`, `python/gbsBindBrepOps.cpp`, `python/tests/test_brep_ops.py`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier8_face_blends.md` (ce palier), `News.md`, cette note |
| Statut | § 9 validé (24 nov. 2026) ; revue de code : § 10, à valider |

## 1. Objectif

Un congé entre deux faces F₁ et F₂ sans arête commune, comme le « Face-Face
Fillet » de CATIA et le « Face fillet » de SolidWorks :

```cpp
fillet_faces(m, s, faces1, faces2, r, options);   // Python : brep.fillet_faces
```

- Chaque côté est une face, ou des faces jointes en tangence (le flanc d'un
  bossage en stade) : elles doivent être reliées par des arêtes le long
  desquelles leurs normales restent à `angular_tol` l'une de l'autre.
- Les faces entre les deux côtés (la **bande**) sont englouties : par-dessus un
  chanfrein, une marche, le fond d'une gorge.
- Rayon constant. Les lois et les sections du palier 7 passent par les options
  (`FilletOptions::section` est déjà lue) ; leurs tests, et les points de loi
  par position, sont pour la PR 2.

## 2. Sources

- B. K. Choi, S. Y. Ju, « Constant-radius blending in surface modelling »,
  *Computer-Aided Design* 21(4), 1989 : congés de faces par surfaces décalées.
- Le comportement de CATIA (« Face-Face Fillet ») et de SolidWorks (« Face
  fillet »). OCCT n'a pas de congé de faces : comparaison par le cas
  équivalent (§ 7).
- Le code des poches englouties (`through_pockets`, `bridge`, `fill` de
  `blend.h`), dont celui-ci reprend la démarche : copie, arête de remplacement,
  congé, contrôle.

## 3. API

```cpp
template <std::floating_point T>
BuildResult<SolidId> fillet_faces(Model<T> &m, SolidId s, std::span<const FaceId> faces1,
                                  std::span<const FaceId> faces2, T r, const FilletOptions<T> &o = {});
```

Un nouvel en-tête, `gbs-brep/face_blend.h`, inclus par `gbs-brep/brep`.
`blend.h` n'est pas modifié.

## 4. Démarche

### 4.1 La bande (`strips`)

On part de chaque face voisine de F₁ (hors des deux côtés) et l'on avance de
face en face sans traverser F₁ ni F₂. Une face rencontrée qui touche les deux
côtés est une **face de bout** (premier passage) ; sinon elle rejoint la bande.
Une candidate est gardée si chaque face des deux côtés touche la bande.

Un second passage range aussi dans la bande les faces qui touchent les deux
côtés. C'est le cas d'une bande de plusieurs faces qui fait le tour d'un côté
fermé : autour d'un stade, les chanfreins coniques touchent le plateau et le
flanc, comme les faces de bout d'une boîte. Les candidates sont essayées de la
plus petite à la plus grande ; la première qui passe le contrôle (§ 4.4) donne
le congé.

### 4.2 La copie fermée par-dessus la bande (`close_over`)

Sur une copie du solide :

1. **Faces gardées prolongées.** Les côtés et les faces de bout, prolongés de
   deux fois la taille de la bande (`extended_surface`, comme les poches).
2. **Rayons.** Les arêtes qui finissent sur la bande entre deux faces gardées :
   une face de côté et une face de bout, deux faces d'un même côté (les arêtes
   de tangence du stade), ou une couture. Chacune est prolongée au-delà de la
   bande (une droite ou un cercle sur un intervalle plus large, une B-spline par
   ses travées d'extrémité, une courbe sur surface par sa pcourbe prolongée sur
   la surface prolongée). Elle est coupée par la face de l'autre côté
   (Gauss-Newton, `meet`), la plus proche du bout de l'arête.
3. **Nœuds.** Ces points de coupe, fondus quand ils coïncident (à un bout de
   boîte, les arêtes avec le dessus et avec le flanc donnent le même coin), sont
   les sommets de l'arête virtuelle. Ils notent les faces sur lesquelles ils
   sont.
4. **Morceaux de l'arête virtuelle.** Pour chaque paire de faces (une par côté)
   qui ont des nœuds en commun : leur intersection, prolongées ; entre deux
   nœuds consécutifs, l'arc le long de la bande. Il est retenu si son milieu est
   près de la bande, à 4 fois la distance des nœuds à elle. Une courbe fermée
   avec un seul nœud fait un tour (le pied d'un bossage).
5. **Sans rayon** : refus. Une face qui fait le tour d'une période a une
   couture (`make_face` refuse une face sans couture), et cette couture finit
   sur la bande : il y a toujours un rayon. Le pied du bossage fusionné garde
   la sienne.
6. **Contours.** Dans chaque face gardée, la suite de co-arêtes le long de la
   bande (rayons et arêtes avec la bande) est remplacée par les rayons
   prolongés et les morceaux, dans l'ordre du contour (`replace_run`, qui
   recale les pcourbes sur les surfaces périodiques). Une boucle (un tour
   complet) prend le sens de la bande qu'elle remplace. Un contour fait
   seulement d'arêtes avec la bande (le trou du plateau sous le bossage) est
   remplacé par la boucle entière.
7. La bande est ôtée de la coque et effacée, avec ses arêtes, les anciens
   rayons et ses sommets. La copie doit passer `check`.

Les pcourbes des rayons prolongés sont projetées (`project_pcurve`) sur les
surfaces prolongées, une par usage pour une couture (amorcée au bout de la
co-arête). Celles des morceaux viennent de l'intersection, ou sont projetées
pour un arc qui passe l'origine d'une courbe fermée.

### 4.3 Le congé

Le congé de l'arête virtuelle sur la copie, celui de `fillet` (ses fins libres
sur les faces de bout, ses chaînes, ses sections). La convexité de l'arête
virtuelle donne le côté : convexe, il ôte (par-dessus un chanfrein) ; concave,
il ajoute (il comble une gorge).

### 4.4 Le contrôle

Celui des poches : le milieu de l'arête virtuelle classé par rapport au
résultat donne le côté ôté (dehors) ou ajouté (dedans). Tous les points de la
bande doivent être de ce côté, et aucune face de la bande ne doit couper le
résultat (`section`).

Une gorge comblée touche le résultat le long de ses arêtes avec les faces de
bout, qui restent des faces du résultat. Ces points et ces sections, sur une
face de bout, sont admis.

Si le congé a été construit mais que le contrôle échoue, le refus est celui du
débordement. La copie est effacée dans tous les cas.

### 4.5 Faces déjà voisines

Si F₁ et F₂ ont des arêtes communes, `fillet_faces` arrondit ces arêtes
(`fillet`), comme prévu au § 0 du palier. Il faut pour cela que chaque face des
deux côtés ait une arête commune avec l'autre côté et qu'aucune bande bornée
ne les sépare ailleurs. Sinon, les faces se touchent par endroits et se font
face par-dessus une bande ailleurs (un chanfrein sur une partie de l'arête
seulement) : refus, plutôt que d'arrondir les seules arêtes communes.

### 4.6 Contrôle des côtés

Avant tout le reste, `sides` vérifie que les faces sont du solide, qu'aucune
n'est des deux côtés et que les faces de chaque côté sont jointes en tangence
(§ 1). Sans ce contrôle préalable, une face donnée des deux côtés passait par
le § 4.5 : elle est voisine d'elle-même le long de ses arêtes, et toutes ses
arêtes étaient arrondies.

## 5. Refus

`NonGeneric`, le modèle inchangé :

- **débordement** : « the fillet spills over the faces between the faces given
  (overflow, not handled): a larger radius, or fillet their edges » ;
- **faces qui ne se coupent pas une fois prolongées** : flancs parallèles, ou
  deux faces que rien ne relie (deux nervures sur une plaque). Le message
  propose le congé tritangent pour le premier cas et dit le second non traité ;
- **bande qui fait le tour du solide** : toutes les candidates sont le reste du
  solide, faute d'une face de bout qui touche les deux côtés (un chanfrein qui
  finit sur un autre chanfrein : au bout, l'un touche le dessus, la face de
  bout le flanc) ;
- **faces voisines par endroits** : des arêtes communes, mais aussi une bande
  ailleurs (§ 4.5) ;
- une arête entre deux faces de bout, ou de la bande à la bande ;
- des faces données qui se touchent déjà à côté de la bande ;
- une face gardée le long de la bande en deux endroits ;
- une bande sans rayon (§ 4.2, point 5).

`InvalidId` : un côté vide, une face morte, une face qui n'est pas du solide,
une face des deux côtés, les faces d'un côté qui ne sont pas jointes en
tangence.

## 6. Tests

`tests_brep_face_blends` (9 cas, 0,45 s). Chaque cas est comparé au congé de
l'arête qui existait avant la bande, sur le même modèle.

| Cas | Attendu | Mesuré |
|---|---|---|
| arête de boîte chanfreinée (0,1), r = 0,2 / 0,3 / 0,5 | même nombre de faces ; volume 2 (1 − (1 − π/4) r²) | 10⁻⁹ |
| même cas, r = 0,13 (< 1,707 d) | refus du débordement, modèle inchangé | — |
| dessus et dessous (faces parallèles) | refus, faces qui ne se coupent pas ; modèle inchangé | — |
| dessus et face arrière (voisines) | congé de leur arête | 10⁻⁹ |
| côtés refusés : côté vide, face d'un autre solide, face des deux côtés, côté de faces non tangentes | `InvalidId`, modèle inchangé | — |
| arête chanfreinée à partir du milieu (un coin découpé) : dessus et flanc voisins par endroits | refus, modèle inchangé | — |
| chanfrein qui finit sur un autre chanfrein | refus de la bande qui fait le tour, modèle inchangé | — |
| pied de bossage chanfreiné (concave), r = 0,12 / 0,2 | bande fermée, couture du cylindre : tore ; congé du pied d'origine | 10⁻⁸ |
| pied d'un stade chanfreiné : flanc de 4 faces | chaîne fermée de 4 congés ; congé du pied d'origine | 10⁻⁷ |
| arête de boîte avec une marche (bande de 2 faces) | 7 faces ; volume de la boîte arrondie | 10⁻⁹ |
| gorge en V à fond plat, r = 0,3 | congé concave, fond comblé ; même solide que la gorge vive arrondie ; r = 0,12 : refus | 10⁻⁹ |
| flancs de deux nervures sur une plaque | refus, modèle inchangé | — |

En Python, la boîte chanfreinée (volume, refus du débordement et d'une face
des deux côtés). La galerie
gagne `fillet_faces_chamfer` et `fillet_faces_boss`.

## 7. Comparaison à OCCT

`tests_brep2_occt.fillet_faces_over_chamfers` : la boîte chanfreinée reprise
par `fillet_faces` contre le congé de l'arête de la boîte sans chanfrein par
`BRepFilletAPI_MakeFillet`. Même solide : mêmes faces, volume à 10⁻⁹.

## 8. Limites

- **Plusieurs bandes** : si les côtés se font face en deux endroits, la plus
  petite bande qui passe le contrôle est arrondie, sans refus (le document du
  palier prévoyait un refus).
- **Faces que rien ne relie** : le refus passe par l'échec des rayons, avec un
  message qui couvre les deux causes (faces parallèles, rien qui relie).
- **Courbe fermée B-spline** : un morceau d'arête virtuelle qui passerait
  l'origine d'une intersection fermée en B-spline n'est pas retenu.
- **Coût** : chaque candidate refusée coûte une copie et un congé. Le cas des
  deux nervures essaie la bande « reste du solide » avant de refuser.
- **Chanfrein qui finit sur un autre chanfrein** (les arêtes du dessus
  chanfreinées ensemble, puis un côté repris en congé) : au bout de la bande,
  deux faces touchent chacune un seul côté, aucune face de bout ne touche les
  deux. Refusé (§ 5) ; il faudrait admettre deux faces de bout, chacune
  prolongée jusqu'à l'autre.
- **Tangence des côtés** : jugée sous `angular_tol` (10⁻⁵ par défaut). Des
  carreaux importés à peine G1 demandent un `angular_tol` plus large, qui
  relâche aussi l'approximation du congé.

## 9. Points à valider

1. **En-tête séparé** `face_blend.h`, sans toucher à `blend.h`.
2. **Bande en deux passages** (faces de bout, puis bande faisant le tour), la
   plus petite candidate qui passe le contrôle.
3. **Arête virtuelle par les rayons** : les arêtes qui finissent sur la bande,
   prolongées jusqu'à la face d'en face ; les morceaux entre leurs points de
   coupe ; un nœud posé sur la courbe fermée quand il n'y a pas de rayon
   (retiré à la revue, § 10, point 4 : il y a toujours un rayon).
4. **Contrôle des poches**, avec les contacts admis sur les faces de bout pour
   une gorge comblée.
5. **Faces voisines** : le congé de leurs arêtes (si elles le sont partout,
   § 10, point 2).
6. **Plusieurs bandes** : la plus petite, sans refus (§ 8), plutôt que le refus
   prévu au palier.

## 10. Revue de code

Revue de la PR (24 nov. 2026), huit points, tous traités. Les deux premiers
ont été reproduits sur le code de la PR avant d'être corrigés.

1. **Face donnée des deux côtés** (important). Elle n'était pas refusée : le
   raccourci des faces voisines (§ 4.5) passait avant les contrôles. Une face
   est voisine d'elle-même le long de ses arêtes, si bien que
   `fillet_faces(dessus, dessus)` arrondissait les quatre arêtes du dessus.
   Correction : `sides` fait tous les contrôles `InvalidId` d'abord (§ 4.6).
2. **Faces voisines par endroits** (important). Une seule arête commune
   suffisait à arrondir ces arêtes en ignorant la bande. Exemple : le dessus
   contre {flanc, face x = 0} arrondissait l'arête dessus/x = 0 et laissait le
   chanfrein. Correction, en deux parties :
   - la tangence des faces d'un côté, annoncée mais pas vérifiée, est contrôlée
     (`InvalidId`) ;
   - le raccourci exige que chaque face ait une arête commune avec l'autre côté
     et qu'aucune bande bornée ne les sépare ailleurs. Sinon, refus. C'est le
     cas d'une arête chanfreinée sur une partie de sa longueur seulement : un
     seul dessus, un seul flanc, voisins jusqu'au chanfrein.
3. **Lien vers une page externe** retiré du document du palier.
4. **Branche sans rayon** jamais exercée. La note la disait atteinte par le pied
   d'un bossage fusionné, qui garde pourtant sa couture. Or `make_face` refuse
   une face qui fait le tour d'une période sans couture (« wire goes around the
   seam of a closed surface »), et cette couture finit toujours sur la bande.
   La branche (des nœuds posés sur les courbes fermées, une portée plus large)
   est remplacée par un refus. Les refus `InvalidId` et celui des faces
   parallèles sont maintenant testés, avec le modèle inchangé.
5. **Chanfrein qui finit sur un autre chanfrein** (§ 8). Le refus parlait de
   faces parallèles, ou d'une arête « de la bande à la bande ». Quand toutes les
   candidates sont le reste du solide, le message dit maintenant qu'aucune face
   de bout ne touche les deux côtés. Un test le vérifie.
6. **Échantillonnage d'une candidate** placé sous `try`, comme pour les poches
   (#244) : une exception écarte la candidate au lieu d'arrêter la recherche.
7. **Description de la PR** : 8 cas annoncés pour 7. Il y en a maintenant 9.
8. **C++23 et simplifications** :
   - `blend_detail::refuse`, `[[noreturn]]`, appelé directement à la place d'un
     lambda `fail` qui ne l'était pas ;
   - `std::ranges::to<std::set>()` pour les côtés et `std::views::keys` pour
     parcourir les faces prolongées ;
   - une clause redondante retirée de la boucle de `path`.

   L'attribut `[[noreturn]]` sur un lambda (P2173) n'a pas été retenu, faute
   d'être sûr de MSVC.

Sondes sans changement de comportement :
- prolongement en arc (une face de bout cylindrique) ;
- un congé de 0,1 remplacé par un congé de 0,3, convexe et concave ;
- un chanfrein dissymétrique ;
- une section conique ;
- le dessus chanfreiné d'un cylindre.

Chacune donne le même solide que le congé de l'arête d'origine.
