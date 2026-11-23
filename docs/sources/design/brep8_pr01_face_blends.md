# Palier 8 — PR 1 : congé entre deux faces

| | |
|---|---|
| PR | branche `feat/brep-face-blends` (sur master) |
| Document de référence | [brep_palier8_face_blends.md](brep_palier8_face_blends.md), § 2.1 à § 2.4 |
| Fichiers | `gbs-brep/face_blend.h` (nouveau), `gbs-brep/brep`, `tests/tests_brep_face_blends.cpp` (nouveau), `gbs-occt/tests/tests_brep2_occt.cpp`, `python/gbsBindBrepOps.cpp`, `python/tests/test_brep_ops.py`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier8_face_blends.md` (ce palier), `News.md`, cette note |
| Statut | § 9 validé (24 nov. 2026) ; § 10 après la revue |

## 1. Objectif

Un congé entre deux faces F₁ et F₂ sans arête commune, comme le « Face-Face
Fillet » de CATIA et le « Face fillet » de SolidWorks :

```cpp
fillet_faces(m, s, faces1, faces2, r, options);   // Python : brep.fillet_faces
```

- Chaque côté est une face, ou des faces jointes en tangence (le flanc d'un
  bossage en stade).
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
5. **Sans rayon** (un cylindre sans couture, après une fusion) : un nœud sur
   chaque courbe fermée où les côtés se coupent tout près de la bande.
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
(`fillet`), comme prévu au § 0 du palier.

## 5. Refus

`NonGeneric`, le modèle inchangé :

- **débordement** : « the fillet spills over the faces between the faces given
  (overflow, not handled): a larger radius, or fillet their edges » ;
- **faces qui ne se coupent pas une fois prolongées** : flancs parallèles, ou
  deux faces que rien ne relie (deux nervures sur une plaque). Le message
  propose le congé tritangent pour le premier cas et dit le second non traité ;
- une arête entre deux faces de bout, ou de la bande à la bande ;
- des faces données qui se touchent déjà à côté de la bande ;
- une face gardée le long de la bande en deux endroits.

`InvalidId` : un côté vide, une face morte, une face qui n'est pas du solide,
une face des deux côtés.

## 6. Tests

`tests_brep_face_blends` (7 cas, 0,4 s). Chaque cas est comparé au congé de
l'arête qui existait avant la bande, sur le même modèle.

| Cas | Attendu | Mesuré |
|---|---|---|
| arête de boîte chanfreinée (0,1), r = 0,2 / 0,3 / 0,5 | même nombre de faces ; volume 2 (1 − (1 − π/4) r²) | 10⁻⁹ |
| même cas, r = 0,13 (< 1,707 d) | refus du débordement, modèle inchangé | — |
| dessus et dessous (faces parallèles) | refus, faces qui ne se coupent pas | — |
| dessus et face arrière (voisines) | congé de leur arête | 10⁻⁹ |
| pied de bossage chanfreiné (concave), r = 0,12 / 0,2 | bande fermée, couture du cylindre : tore ; congé du pied d'origine | 10⁻⁸ |
| pied d'un stade chanfreiné : flanc de 4 faces | chaîne fermée de 4 congés ; congé du pied d'origine | 10⁻⁷ |
| arête de boîte avec une marche (bande de 2 faces) | 7 faces ; volume de la boîte arrondie | 10⁻⁹ |
| gorge en V à fond plat, r = 0,3 | congé concave, fond comblé ; même solide que la gorge vive arrondie ; r = 0,12 : refus | 10⁻⁹ |
| flancs de deux nervures sur une plaque | refus, modèle inchangé | — |

En Python, la boîte chanfreinée (volume, refus du débordement). La galerie
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

## 9. Points à valider

1. **En-tête séparé** `face_blend.h`, sans toucher à `blend.h`.
2. **Bande en deux passages** (faces de bout, puis bande faisant le tour), la
   plus petite candidate qui passe le contrôle.
3. **Arête virtuelle par les rayons** : les arêtes qui finissent sur la bande,
   prolongées jusqu'à la face d'en face ; les morceaux entre leurs points de
   coupe ; un nœud posé sur la courbe fermée quand il n'y a pas de rayon.
4. **Contrôle des poches**, avec les contacts admis sur les faces de bout pour
   une gorge comblée.
5. **Faces voisines** : le congé de leurs arêtes.
6. **Plusieurs bandes** : la plus petite, sans refus (§ 8), plutôt que le refus
   prévu au palier.
