# Palier 8 — PR 2 : bande partielle, lois et sections des congés de faces

| | |
|---|---|
| PR | branche `feat/brep-face-blends-2` (sur master 4f4a508, après #246) |
| Document de référence | [brep_palier8_face_blends.md](brep_palier8_face_blends.md), § 2.5 et § 2.6 (question 12 validée le 24 nov. 2026) |
| Fichiers prévus | `gbs-brep/face_blend.h`, `gbs-brep/radius_law.h`, `gbs-brep/blend.h`, `gbs-brep/fillet.h` (résolution des points par position), `tests/tests_brep_face_blends.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `python/gbsBindBrepOps.cpp`, `python/tests/test_brep_ops.py`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier8_face_blends.md`, `News.md`, cette note |
| Statut | § 9 et § 10 validés (24 nov. 2026) |

## 1. Objectif

Deux ajouts à `fillet_faces` (PR 1, #246) :

1. **La bande partielle** (§ 2.6 du palier). Les deux faces se touchent déjà sur
   une partie de leur longueur, et une bande les sépare ailleurs, **contiguë à
   leur arête commune**. Le congé engloutit la bande. C'est le manque
   d'OpenCascade, où un congé ne supprime pas de face. Trois formes :
   - une face triangulaire qui finit **en pointe** sur l'arête commune ;
   - un **chanfrein arrêté** au milieu de l'arête, fermé par une face
     d'extrémité ;
   - une bande **au milieu** de l'arête, entre deux arêtes communes.
2. **Les lois et les sections** (§ 2.5 du palier) :
   - un rayon en loi (`BlendLaw`), avec des points donnés **par position** :
     un point de l'espace et une valeur ;
   - les sections du palier 7 (conique, corde constante, G2) sur l'arête
     virtuelle, avec leurs tests.

   Les points par position servent aussi aux congés d'arête.

La gorge concave, prévue au départ pour cette PR, est déjà faite et testée dans
la PR 1.

## 2. Sources

- Le comportement de **SolidWorks** : ses congés de faces (« Face fillet »)
  sont plus avancés que ses congés d'arête. Ils engloutissent les petites faces
  entre les deux faces choisies, même là où les faces se touchent déjà. Les
  *Help points* y lèvent une ambiguïté par un point posé sur le modèle.
- **CATIA**, « Face-Face Fillet » : les deux faces choisies par l'utilisateur,
  et les faces entre elles supprimées. Les points de rayon variable y sont
  posés sur le modèle.
- Le code de la PR 1 (`strips`, `close_over`, `through_strip`) et des poches
  (`law_on_copies`, `section_on_copies` de `blend.h`).

## 3. API

Aucune fonction nouvelle pour la bande partielle :
`fillet_faces(m, s, faces1, faces2, r)` la reconnaît.

Pour les lois :

```cpp
fillet_faces(m, s, faces1, faces2, const RadiusLaw<T> &law, options);
fillet_faces_chordal(m, s, faces1, faces2, c, options);  // corde constante (§ 4.6)
BlendPoint<T>::at_point(point<T, 3> p, T value);         // un point par position
```

- `BlendPoint` gagne une position facultative : `edge` invalide et `where` le
  point.
- *Écart à l'implémentation* : `BlendPoint::at` était prévu, mais `at` est déjà
  le champ de la fraction (C++ et Python) ; la fabrique s'appelle `at_point`.
  La corde constante n'est pas une `BlendSection` mais `fillet_chordal` : son
  pendant `fillet_faces_chordal` est ajouté pour le test du § 4.6.
- Une loi mêle librement les deux formes : (arête, fraction) et position.
- Python :
  - `brep.fillet_faces(m, s, faces1, faces2, law, options)` ;
  - `brep.fillet_faces_chordal(m, s, faces1, faces2, chord, options)` ;
  - `brep.BlendPoint.at_point(point, value)`.

## 4. Démarche

### 4.1 La règle de la bande partielle

Aujourd'hui, quand les deux côtés ont une arête commune, `fillet_faces` :
- fait le congé de cette arête si aucune bande ne les sépare ailleurs ;
- sinon, refuse (« faces given that meet along edges in places only… »).

La PR remplace ce refus par la règle suivante.

**Règle.** Soit e l'arête commune de F₁ et F₂ (ou leur chaîne d'arêtes
communes). Une face G est **englobée dans la bande** si les trois conditions
suivantes sont réunies :

1. G n'est ni dans F₁ ni dans F₂ ;
2. G touche F₁ et F₂ ;
3. ses arêtes avec F₁ et F₂ **finissent sur un sommet de e**.

Une telle face est sous le prolongement de e. Exemples :
- la face en pointe, dont la pointe est ce sommet ;
- le triangle qui ferme un chanfrein arrêté.

La bande partielle est la composante connexe des faces englobées et des faces
qui les relient, sans traverser F₁, F₂ ni e. Les **faces de bout** restent celles
de la PR 1 : elles touchent F₁ et F₂ loin de e, et bornent la bande de l'autre
côté, au bout du solide.

Ce qui ne change pas :
- une face qui touche F₁ et F₂ sans finir sur un sommet de e est une face de
  bout, comme dans la PR 1 (le bout d'une boîte) ;
- sans bande contiguë à e, le comportement d'aujourd'hui reste : le congé des
  arêtes communes.

### 4.2 L'arête fondue

Sur la copie de la PR 1 :

1. L'arête virtuelle part du sommet de e où finit la bande et va jusqu'à la face
   de bout. C'est l'intersection de F₁ et F₂ prolongées, la même paire de
   surfaces que e : elle prolonge exactement la courbe de e.
2. **Fusion** : e et l'arête virtuelle deviennent une seule arête.
   - Les rayons (§ 4.2 de la PR 1) donnent le sommet du bout.
   - La courbe est l'intersection recalculée sur la réunion, depuis l'autre
     sommet de e.
   - Le sommet entre e et l'arête virtuelle est ôté, avec ses co-arêtes, et
     les contours de F₁ et F₂ sont recousus.
3. Pour une bande au milieu, e₁, l'arête virtuelle et e₂ sont fondues en une
   seule arête.
4. Sur des faces courbes, l'intersection est recalculée de bout en bout :
   - une courbe B-spline, approchée sous `tol` ;
   - un cercle ou une droite exacts quand les surfaces le permettent (plan et
     cylindre coaxial, deux plans).

La copie passe `check`, puis le congé de l'arête fondue est fait par la
machinerie de la PR 1. Le résultat est d'un seul tenant : un cylindre ou un tore
exact, sans station de chaîne au milieu.

### 4.3 Le contrôle

C'est celui de la PR 1, sur la bande partielle seule. Le congé est refusé, avec
le message du débordement, si l'une de ces conditions échoue :
- les contacts restent dans F₁ et F₂ ;
- aucune face de la bande n'est coupée par le congé ;
- les points de la bande sont du côté ôté (ou ajouté).

Exemple : dans le U de la question, pour une face en pointe de 0,12 au bout du
U, l'arc du congé est à r(√2 − 1) du coin. La face tient sous le congé dès que
r ≥ 0,3 ; à r = 0,2, le refus est celui du débordement.

### 4.4 `fillet` d'une arête

`fillet` d'une arête **ne franchit pas** une bande de lui-même (question 12).
Son refus au bout de l'arête, quand une bande partielle y est contiguë, est
complété :

« a blend wider than a face at an end of its edge: faces between its faces
there, fillet_faces swallows them »

### 4.5 Lois par position

Un point par position (p, v) est résolu en (arête, fraction) avant
l'utilisation de la loi :

- **congé d'arête** : l'arête des arêtes données dont la courbe est la plus
  proche de p (le pied par `extrema_curve_point`), et la fraction de longueur
  au pied ;
- **congé de faces** : sur l'arête virtuelle (fondue ou non) de la copie,
  après `close_over`. *Écart* : plutôt que `law_on_copies`, les points
  (arête, fraction) d'une loi de congé de faces (et d'une loi de ρ) sont
  d'abord changés en points par position, à leur place sur le solide donné,
  puis résolus sur l'arête virtuelle de la copie : une fraction d'une arête
  commune fondue n'a plus de sens sur l'arête fondue, sa position en garde un.

Refus (`InvalidBounds`, comme les autres refus des lois ; il n'y a pas de code
`InvalidArgument`) dans deux cas :
- un point à plus de `max(10 r, la taille de l'arête)` de toute arête : un
  point posé ailleurs ;
- deux points qui tombent au même pied avec des valeurs différentes.

La variation (linéaire, lisse) et l'abscisse curviligne restent celles du
palier 5.

### 4.6 Sections

`FilletOptions::section` est déjà lue par `fillet_faces` (PR 1). La PR ajoute
les tests, sur l'arête chanfreinée et sur la bande partielle :
- conique ρ ;
- corde constante (`fillet_faces_chordal` ; #248, qui en fait un cylindre exact
  entre deux plans, n'est pas encore fusionnée : le test compare au congé de
  l'arête non chanfreinée à la même corde, quelle que soit la surface) ;
- G2.

## 5. Documentation de la règle

La règle du § 4.1 doit être claire pour l'utilisateur. Elle sera écrite à trois
endroits :

- **En-tête de `fillet_faces`** : la bande, les faces de bout, la bande
  partielle (la face englobée, avec ses trois conditions), le refus du
  débordement, et ce que `fillet` d'une arête ne fait pas.
- **Docstring Python** : la même chose, en plus court, avec l'exemple du U.
- **Guide** (branche `docs/guide`, tenue par la session de documentation) :
  - une section « Which faces a face fillet swallows » avec trois figures :
    la face de bout gardée (boîte), la face en pointe englobée, le chanfrein
    arrêté englobé ;
  - la règle en une phrase ;
  - le refus du débordement, avec le rayon minimal (r(√2 − 1) contre la
    largeur de la face) ;
  - la différence avec `fillet` d'une arête.

  Le texte et les figures du guide seront transmis à cette session avec le
  commit.

Les figures sont de vrais rendus gbs, comme les notes d'illustration.

## 6. Tests

`tests_brep_face_blends` gagne les cas suivants. Chacun est comparé au congé de
l'arête vive du même solide sans la bande.

| Cas | Attendu | Mesuré |
|---|---|---|
| U, face en pointe entre paroi et fond (0,12 au bout, pointe à mi-longueur), r = 0,3 / 0,5 | mêmes faces que la référence, un cylindre d'un seul tenant, volume à 10⁻⁹ | 11 faces, 1 cylindre ; écart 0 ; 5 ms |
| même cas, r = 0,2 | refus du débordement, modèle inchangé | « spills over » ; modèle inchangé |
| U, chanfrein arrêté à mi-longueur (face d'extrémité triangulaire), r = 0,3 / 0,5 | idem | écart 0 ; 7 à 8 ms |
| U, chanfrein au milieu de l'arête (deux arêtes communes) | idem, arête fondue e₁ + virtuelle + e₂ | écart 0 ; 9 à 10 ms |
| boîte à chanfrein arrêté (refusée par la PR 1) | le congé de l'arête de la boîte | volume exact à 10⁻⁹ |
| face de bout d'une boîte (PR 1) | inchangée : gardée, prolongée | tests de la PR 1 inchangés |
| pied de bossage chanfreiné sur un demi-tour (cylindre et plan) | un tore d'un seul tenant ; volume du congé du pied non chanfreiné à 10⁻⁸ | mêmes faces ; écart 0 ; 16 à 17 ms |
| `fillet` de l'arête vive du U à face en pointe | refus, le message nomme `fillet_faces` | r = 0,2 / 0,3 / 0,5 : refusé (« …: faces between its faces there, fillet_faces swallows them ») ; idem au chanfrein arrêté |
| loi par position sur la boîte chanfreinée (rayons 0,2 → 0,4) | même solide que la loi (arête, fraction) sur la boîte sans chanfrein | linéaire : 7·10⁻¹⁴ ; lisse : 1·10⁻¹¹ ; 10 à 15 ms |
| loi par position sur un congé d'arête | même solide que la loi (arête, fraction) équivalente | écart 0 |
| point de loi loin de tout, deux points au même pied | refus `InvalidBounds` | `fillet_faces` et `fillet` : refusés, modèle inchangé |
| sections conique ρ = 0,5, corde constante, G2 sur l'arête chanfreinée | volumes du congé de l'arête non chanfreinée avec la même section (palier 7) | boîte et U à face en pointe : écarts ≤ 2·10⁻¹³ ; 5 à 17 ms |

OCCT (`tests_brep2_occt`) : le U à face en pointe repris par `fillet_faces`,
contre `BRepFilletAPI_MakeFillet` sur l'arête vive du U sans face. Même solide,
volume à 10⁻⁹ (r = 0,3 et 0,5 : passé).

Python :
- la face en pointe (volume, refus du débordement) ;
- une loi par position.

La galerie gagne `fillet_faces_tapered` (le U à face en pointe).

## 7. Illustration

Note d'illustration (montrée à la validation, hors du dépôt), rendus réels de
la branche :
- la règle : la face de bout d'une boîte gardée, le chanfrein englouti ;
- le U à face en pointe, à chanfrein arrêté et à chanfrein au milieu, avant et
  après (volume égal au congé de l'arête vive) ;
- le refus à r = 0,2, en coupe au bout du U ;
- une loi lisse 0,3 → 0,5 posée par deux points sur le U à face en pointe.

Les sections sont dans les tests (§ 6), sans figure.

## 8. Estimation

| Partie | Heures |
|---|---|
| bande partielle : règle, arête fondue, refus de `fillet`, tests, OCCT | 10 |
| lois par position (congés d'arête et de faces), sections, tests | 6 |
| Python, galerie, documentation de la règle, illustration | 2 |
| **Total** | **18** |

## 9. Points à valider

1. **La règle de la face englobée** (§ 4.1) : une face qui touche les deux
   côtés et dont les arêtes avec eux finissent sur un sommet de leur arête
   commune est dans la bande ; sinon, elle reste une face de bout.
2. **L'arête fondue** : l'arête commune et l'arête virtuelle en une seule
   arête sur la copie, recalculée de bout en bout (un congé d'un seul tenant),
   plutôt qu'une chaîne de deux arêtes.
3. **`fillet` d'une arête inchangé**, son refus complété pour nommer
   `fillet_faces`.
4. **Points de loi par position** résolus au pied le plus proche, sur les
   arêtes données (congé d'arête) ou sur l'arête virtuelle de la copie (congé
   de faces). Refus d'un point trop loin, ou de deux points au même pied.
5. **La documentation de la règle** (§ 5) : en-tête C++, docstring Python, et
   une section du guide avec trois figures, transmise à la session de
   documentation.

## 10. Revue de code

Revue de la PR (24 nov. 2026), six points, tous traités. Les trois premiers
ont été reproduits par des sondes avant d'être corrigés.

1. **Loi d'une seule valeur par position** (important). `fillet_geometry`,
   `chamfer_geometry` et la section conique d'une arête lisaient la valeur
   d'une loi constante parmi ses points tels que donnés. Un point par position
   n'a pas d'arête : la recherche passait la fin des points. Le congé de 0,2
   d'une boîte était refusé (« the radius must be larger than the
   tolerance »). Correction : la valeur est lue dans la loi posée sur l'arête
   (`fillet_detail::edge_value`).
2. **Corde constante** (important). `fillet_faces_chordal` recopiait la
   géométrie de `fillet_chordal`. La PR 248 change celle-ci (congé exact entre
   deux plans) : après les deux fusions, l'écart de volume entre les deux
   passait de 1,7e-13 à 2,8e-10, sous la tolérance de 1e-9 du test.
   Correction : une seule géométrie, `blend_detail::chordal_blend`, pour les
   deux ; le test des sections compare à 1e-11. À la fusion de master (PR 248),
   le cas des deux plans y a été reporté : `fillet_faces_chordal` le reçoit
   aussi.
3. **Arête commune de la bande à la bande** (important). Le pied d'un bossage
   chanfreiné sur un demi-tour était refusé quand le chanfrein couvre la
   couture du cylindre : le reste du pied est alors une seule arête dont les
   deux bouts sont sur la bande. Elle est maintenant sous l'arête virtuelle de
   bout en bout, sans sommet de celle-ci (l'arête virtuelle fait le tour, comme
   au pied entièrement chanfreiné). Le test du bossage couvre le demi-anneau
   des deux côtés de la couture.
4. **Refus rendu** (mineur). Parmi plusieurs bandes candidates, le refus était
   celui de la première essayée. C'était souvent une bande parasite plus
   petite : le triangle qui ferme un chanfrein arrêté, borné par ce chanfrein,
   échoue avant d'être fermé et parlait de faces parallèles. Le refus rendu
   est maintenant celui de la première bande fermée sur sa copie.
5. **Deux bandes de part et d'autre d'une arête commune** (mineur). Dans
   `fillet_faces`, le refus disait d'appeler `fillet_faces`. Il dit maintenant
   que ce cas n'est pas traité.
6. **Lien vers l'illustration** retiré du § 7 (une page hors du dépôt).
