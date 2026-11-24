# Congé de faces : bande du mauvais côté de l'arête virtuelle

| | |
|---|---|
| PR | branche `fix/brep-face-blend-undercut` (sur master dcc81fb) |
| Document de référence | [brep_palier8_face_blends.md](brep_palier8_face_blends.md), § 2.4 (contrôle d'engloutissement) ; notes [brep8_pr01_face_blends.md](brep8_pr01_face_blends.md) et [brep8_pr02_partial_strips_laws.md](brep8_pr02_partial_strips_laws.md) |
| Fichiers | `gbs-brep/face_blend.h` (`through_strip`), `gbs-brep/blend.h` (`through_pockets`), `tests/tests_brep_face_blends.cpp`, `python/tests/test_brep_ops.py`, `News.md`, cette note |
| Statut | À valider (§ 5) |

## 1. Le défaut

`fillet_faces` entre les flancs d'une nervure **en contre-dépouille** de −10°
(plus large en haut qu'en bas), sur une plaque, réussissait pour tout rayon
essayé, de 0,5 à 2 fois le rayon du congé tritangent. `check` passait, mais le
solide n'avait plus que 6 faces au lieu de 10, et 0,89 de volume au lieu de
8,85 : presque toute la plaque avait disparu. Le défaut vient des PR #246 et
#249 ; il est antérieur au congé tritangent (#251), qui a permis de le
trouver.

## 2. La cause

Deux bandes candidates séparent les flancs (§ 2.1 du palier) :

1. **le dessus de la nervure** (1 face). Les flancs, prolongés au-dessus du
   dessus, s'écartent : ils ne se rencontrent pas de ce côté. La bande est
   refusée (« faces given that do not meet once continued »), à juste titre ;
2. **la plaque** (5 faces : les deux parties du dessus de la plaque, ses deux
   côtés, son fond), bornée par les deux faces de bout. Prolongés vers le bas,
   les flancs se rencontrent sous la nervure, à z = 0,30, **dans la plaque**.
   La copie fermée sur cette bande est valide : c'est la nervure seule,
   prolongée jusqu'à son arête vive en V. Le congé de cette arête se construit.

Le contrôle d'engloutissement (§ 2.4 du palier) demandait que **tout point de
la bande soit du même côté du résultat que l'arête virtuelle**, hors du
résultat pour un congé qui ôte de la matière. C'est vrai de la plaque : la
copie, donc le résultat, n'en contient rien. Le contrôle passait.

Le palier demandait davantage : la bande « tout entière dans le volume que le
congé ôte (ou ajoute) **entre sa surface et l'arête virtuelle** ». Le code n'en
vérifiait que la moitié, le côté du résultat. Un point loin du congé, hors du
solide, le passait aussi.

La cause supposée par la revue (« les flancs prolongés se rencontrent sous la
nervure, dans la plaque ») est juste ; il faut y ajouter que c'est une **autre
bande**, la plaque, qui est retenue, et que le contrôle d'engloutissement
n'exigeait pas de la bande qu'elle soit sous l'arête virtuelle.

## 3. La correction

Dans `through_strip`, une fois la copie fermée sur la bande et **avant** le
congé : tout point de la bande doit être **dans la copie** (sur son bord admis)
si l'arête virtuelle est convexe (le congé ôtera de la matière), **hors de la
copie** si elle est concave (il en ajoutera). C'est le volume entre la bande et
l'arête virtuelle, celui que le congé doit traverser. Sinon, refus :

« faces given that meet, once continued, on the far side of the faces between
them (the flanks of a rib in undercut meet below it, not above its top): no
fillet between them over those faces (fillet_full_round replaces the face
between them) »

Le contrôle se fait avant celui de la validité de la copie, et avant le congé :
une copie fermée du mauvais côté peut être invalide, ou son congé refusé pour
une autre raison (un contact hors de sa face, à r = 2 fois le rayon du
tritangent), et le refus ne disait alors rien de la cause. Parmi les bandes
essayées, ce refus l'emporte sur ceux des bandes qui ne sont pas allées aussi
loin (§ 2.4 de la PR 1, point 4 de sa revue).

**Côté concave.** Une bande du mauvais côté d'une arête virtuelle concave ne
pouvait pas être engloutie à tort : le résultat n'y ajoute à la copie que la
matière du congé, et la bande, loin de lui, n'est pas dans le résultat. Le
contrôle d'engloutissement la refusait déjà, comme un débordement. Le cas testé,
une rainure en queue d'aronde (plus large au fond), est refusé plus tôt : la
copie fermée sur le dessus du bloc est invalide. Le message est celui de la
copie invalide, qui ne nomme pas la cause.

**Les poches** (`through_pockets` de `blend.h`, [brep_swallowed_faces.md](brep_swallowed_faces.md))
avaient le même contrôle, du seul côté du résultat. Le même complément y est
ajouté : les points de la poche doivent aussi être dans la copie remplie (hors
d'elle pour un congé concave). Aucun cas de poche fausse n'a été construit : le
pont d'une poche relie deux arêtes données voisines, et ne peut guère passer de
l'autre côté. Le complément ne change aucun test des poches.

**Bandes partielles** (PR 2) : elles passent par `through_strip` ; le même
contrôle s'y applique. Leurs tests sont inchangés.

## 4. Tests

`tests_brep_face_blends` :

| Cas | Attendu | Mesuré |
|---|---|---|
| nervure en contre-dépouille de −10°, `fillet_faces` des flancs, r de 0,5 à 2 fois le rayon du tritangent (0,252) | refus « on the far side », modèle inchangé | 6 rayons refusés (avant : 6 faces, 0,89 de volume sur 8,85) |
| nervure à dépouille de +10° (les flancs se rencontrent au-dessus du dessus), r = 0,36 et 0,4 | le dessus englouti ; le congé de l'arête au sommet de la nervure prolongée | même volume, à 10⁻⁹ |
| rainure en queue d'aronde, `fillet_faces` des parois, r = 0,2 et 0,4 | refus, modèle inchangé | refus de la copie invalide |
| tests des PR 1 et 2 (chanfrein, marche, gorge en V, bossages, bandes partielles, lois, sections) | inchangés | inchangés |

Python (`test_fillet_faces_far_side`) : le refus en contre-dépouille, le modèle
inchangé ; le dessus englouti avec une dépouille de +10°.

## 5. Point à valider

1. **Le contrôle avant le congé**, sur la copie fermée : la bande sous l'arête
   virtuelle (dans la copie pour une arête convexe, hors d'elle pour une arête
   concave), avec un refus qui nomme le cas. *Recommandé* : c'est ce que
   demandait le § 2.4 du palier ; un succès avec un solide faux n'est plus
   possible du côté convexe, et le côté concave était déjà refusé.
2. **Le même complément pour les poches**, sans cas qui le montre.
   *Recommandé* : le même raisonnement, sans coût.
3. **Le message de la queue d'aronde** reste celui de la copie invalide. Le
   rendre explicite demanderait un contrôle géométrique de la bande sans la
   copie (par exemple le côté des faces prolongées où elle se trouve).
   *Recommandé* : le laisser, le modèle restant inchangé.
