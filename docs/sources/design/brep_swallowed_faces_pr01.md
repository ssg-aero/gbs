# Congés qui engloutissent une face — PR 1 : poches, pont, débordement

| | |
|---|---|
| PR | branche `feat/brep-swallowed-faces` |
| Document de référence | [brep_swallowed_faces.md](brep_swallowed_faces.md), § 3 à § 5 |
| Fichiers | `gbs-brep/blend.h`, `tests/tests_brep_fillet.cpp`, `News.md`, `docs/sources/design/brep_swallowed_faces.md` (statut), cette note |
| Illustration | [rendus réels](https://claude.ai/artifact/MDHQ8CTKzXhi1UzBEXWyfG) |
| Statut | À valider (§ 9) |

## 1. Objectif

Première PR : les **congés** (rayon constant ou loi, arêtes convexes ou
concaves) qui passent sous une poche et l'engloutissent ; le **débordement**
détecté et refusé sous son propre message. Les chanfreins, Python, la galerie et
la comparaison à OpenCascade viennent avec la PR 2.

Le cas d'origine (§ 1 du document) : le bord d'un cylindre coupé par un plan
incliné de 30°, trois arêtes, est construit pour r ≥ 0,108 (l'engloutissement
commence à r = 0,1057), refusé comme débordement de 0,04 à 0,105, et construit en
onglets jusqu'à 0,03, comme avant.

## 2. Sources

Celles du document de référence. Rien de nouveau : la PR réutilise la marche de
la ligne des centres, l'intersection de surfaces, `replace_run`, la
classification d'un point dans un solide et la section.

## 3. API

Aucun changement : `fillet(model, solid, edges, r, options)` et
`fillet(model, solid, edges, law, options)`. Une arête donnée qui est dans une
poche engloutie disparaît sans erreur (question 3 du § 8 validée).

## 4. Démarche

### 4.1 Remplir la poche, puis arrondir (écart au § 4.1 du document)

Le document proposait de construire un seul traitement le long d'une **arête
virtuelle**, arête, pont, arête. La PR rend cette arête **réelle** sur une
copie : la poche est remplie (les deux faces continuées jusqu'à leur
intersection), puis le congé ordinaire arrondit la chaîne arête, pont, arête.
Une fois les contrôles passés, le résultat est exactement le même : la matière
ajoutée en remplissant la poche est entre ses faces et le pont, donc dans le
volume que le congé ôte (pour une arête concave, celle qui est retirée est dans
le volume qu'il ajoute). C'est le congé du solide sans poche, par exemple celui
du cylindre de hauteur 3 coupé par le même plan.

Ce choix ne touche pas à la construction des congés (fins, onglets, chaînes,
coins) : `through_pockets` est une passe avant `blend`. La chaîne garde une face
de congé par arête, pont compris : trois faces pour le bord du cylindre et pour
le méplat, et non une seule.

### 4.2 Étapes (`through_pockets`)

1. **Candidats** (`pockets`, § 3.1 du document) sur le solide donné ; aucun :
   le congé ordinaire, sans copie.
2. **Copie** du solide ; les candidats y sont retrouvés et **échantillonnés**
   (`sample`) : points intérieurs de chaque face (`interior_points`), sept points
   sur chaque arête et ses sommets, sauf v et w.
3. **Pont** (`bridge`, § 3.2) : les surfaces des deux faces continuées
   (`extended_surface`, deux fois l'étendue de la poche), leur intersection, la
   branche qui passe par v et w (à 10 tol près), l'arc du côté de la poche (le
   plus proche du centre des échantillons). Un candidat sans pont est écarté,
   comme celui dont le pont longe une arête des deux faces : c'est le « reste du
   solide », pris pour une poche en partant de l'autre bout de l'arête (vu sur
   la boîte au méplat). Les candidats retenus ne partagent aucune face.
4. **Remplissage** (`fill`, § 4.2) : les surfaces continuées deviennent celles
   des deux faces (mêmes paramètres) ; l'arête du pont remplace le parcours de
   chacune (`replace_run`) ; les faces de la poche sortent de la coque et sont
   effacées avec leurs arêtes et leurs sommets ; le solide rempli est vérifié
   (`check`).
5. **Congé** du solide rempli, sur les arêtes données hors des poches et sur
   les ponts.
6. **Contrôles** (§ 3.3) : le milieu du pont donne le côté attendu (hors du
   résultat pour un congé convexe, dedans pour un concave, jamais dessus) ;
   tous les échantillons de la poche doivent y être (`classify`), et la section
   du résultat par chaque face de la poche (celle du solide donné) doit être
   vide (`section` : ni arête, ni point, ni arête commune, complète).
7. Un échec à l'une de ces étapes (exception comprise) : retour en arrière de la
   copie, et congé du solide tel quel. Si l'étape 6 a échoué et que ce congé est
   refusé (`NonGeneric`), le refus devient **« a blend that rolls over a sharp
   edge next to it (overflow, not handled) »**, sur les arêtes de la poche.

### 4.3 Pourquoi la section en plus des points

Les points seuls laissaient passer r = 0,105 sur le cylindre : la boule dépasse
le dessus de 5·10⁻⁴ dans une petite zone que les échantillons du croissant ne
touchent pas. Le congé rempli traversait alors le croissant, et le résultat
aurait gardé de la matière au-dessus de z = 2. La section du résultat par les
faces de la poche le voit : le contrôle 2 du document (« la surface ne
rencontre aucune face de la poche ») est donc fait par intersection, comme
proposé (question 5), et les points restent pour le débordement par les
contacts.

### 4.4 Jonction de chaîne : même arc, autres paramètres (`joint`)

Au bout d'une arête et du pont, les deux congés finissent sur le même arc. Leurs
surfaces le paramètrent autrement quand l'une le met en un segment rationnel et
l'autre en deux : à r = 0,11, les arcs coïncident à 10⁻⁸ près, mais à un même
paramètre leurs points sont à 3,8·10⁻³. `joint` comparait les points à un même
paramètre et refusait. Désormais, quand les contacts se rejoignent et que
chaque point de l'arc du premier est **sur l'arc** du second (distance à la
courbe, en 15 points, sous la limite des jonctions), la pcourbe de l'arc sur le
second congé est **projetée** depuis l'arc du premier (`project_pcurve`, à
tol / 10). La distance est mesurée à l'arc, pas à la surface : deux arcs d'une
même boule dans deux plans (un congé à loi dont la ligne des centres a un coin)
ont les mêmes bouts et s'écartent au milieu, mais l'arc du premier reste près
de la surface du second ; ils gardent leur lune (`oblong_chain`, dans
`tests_brep_variable_fillet_solid`, l'a montré pendant le développement). L'arête commune a alors des pcourbes de
même paramètre sur les deux congés, comme `add_edge` l'exige. Le défaut existait
avant cette PR, mais aucune chaîne ne le rencontrait.

## 5. Refus et limites

| Cas | Comportement |
|---|---|
| le congé coupe la poche (contacts ou surface) | refus « overflow », si le congé du solide tel quel est refusé aussi |
| il la touche à peine (boule tangente) | idem : la section voit le point de contact |
| couture d'une des deux faces qui entre dans la poche (cylindre dont la couture est en x = −1) | pas de poche (le parcours de la face suit la couture) : comportement d'avant |
| poche dont le pont n'est pas sur une seule branche, B-spline fermée traversée par le pont | pas de poche |
| chanfreins | pas encore (PR 2) |

Une arête seule dont les deux bouts bordent la même poche (e = e₂) est prise en
compte par le code, mais n'est pas testée : sur le cylindre coupé, c'est le cas
de la couture dans le croissant, qui n'est pas traité.

## 6. Tests

`tests/tests_brep_fillet.cpp` :

| Test | Cas | Référence | Écart |
|---|---|---|---|
| `fillet_tilted_cut_rim` | bord du cylindre, r = 0,02 et 0,03 (onglets) | volume entre le solide et le solide moins le prisme du congé de 60° sur tout le bord | — |
| `swallow_tilted_cut_rim` | r = 0,11 ; 0,14 ; 0,2, les trois arêtes puis les deux arcs seuls | congé du cylindre de hauteur 3 coupé par le même plan (une arête fermée, pas de croissant) ; 3 faces de congé, croissant ôté | 10⁻⁶ (mesuré 3·10⁻⁸) |
| | r = 0,05 et 0,105 | refus « overflow » | — |
| `swallow_flat_on_box_edge` | boîte 2 × 2 × 2, méplat de profondeur 0,1 sur z ∈ [0,8 ; 1,2], r = 0,5 | 8 − 2 (1 − π/4) r², et le congé de la boîte sans méplat | 10⁻⁹ |
| | même méplat, loi de 0,4 à 0,6 le long de la chaîne | même loi sur l'arête de la boîte | 2·10⁻⁵ (congés approchés à 10⁻⁵) |
| | méplat de 0,3 : r = 0,2 | fins libres sur le plancher et le plafond de l'encoche, pas d'engloutissement | — |
| | méplat de 0,3 : r = 0,5 | refus « overflow » | — |
| | concave : bosse de 0,1 dans l'angle d'un L, r = 0,5 | congé de l'arête intérieure du L sans bosse | 10⁻⁹ |

Les suites `tests_brep_fillet`, `tests_brep_vertex_blends` et
`tests_brep_variable_fillet_solid` passent (34, 19 et 10 cas ; 9,0 s, 23,8 s et
10,7 s en local).

## 7. Coût

Sans poche candidate (presque toujours), la passe ne fait que parcourir les
bords des faces des arêtes données : pas de copie. Avec une poche, le congé est
construit une fois sur la copie remplie, plus une fois sur le solide tel quel
si les contrôles échouent. Les contrôles coûtent une classification par
échantillon et une section par face de la poche : les trois cas du bord du
cylindre (deux congés, plus la référence) prennent moins d'une seconde.

## 8. Écarts au document

1. Remplir la poche sur une copie au lieu de construire une arête virtuelle
   (§ 4.1) : même résultat, rien de changé dans la construction des congés ;
   une face de congé par arête de la chaîne.
2. Correction de `joint` (§ 4.4), découverte en route.
3. Couture d'une face dans la poche : non traitée (§ 5).

## 9. Points à valider

1. Le remplissage sur une copie à la place de l'arête virtuelle (§ 4.1, écart 1).
2. Les contrôles : points échantillonnés et section par les faces de la poche
   (§ 4.3).
3. La projection de l'arc à la jonction d'une chaîne (§ 4.4).
4. Le message « overflow », donné seulement quand le congé du solide tel quel
   est refusé aussi.
5. La couture dans la poche laissée hors de cette PR.
