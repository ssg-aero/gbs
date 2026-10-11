# Lecture STEP — PR 8 : comparaison avec OpenCascade

| | |
|---|---|
| PR | branche `feat/step-occt-compare` |
| Document de référence | [step_reader.md](step_reader.md), § 10 (point 2) |
| Fichiers | `gbs-occt/tests/tests_step_occt.cpp` (nouveau), `gbs-io/step/reader.h`, `gbs-brep/pcurve.h`, `tests/tests_step_read.cpp`, `tests/in/step/occt_partial_torus.stp`, `tests/in/step/occt_loft.stp` (nouveaux), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Valider le lecteur sur des fichiers **réalistes, produits par un écrivain de
référence**, sans données industrielles à versionner. Chaque forme est
construite par OCCT, écrite par `STEPControl_Writer` (AP214), puis relue par
gbs et par `STEPControl_Reader`. On compare :

- les nombres de **faces**, d'**arêtes non dégénérées** et de **sommets** ;
- les **volumes** (`BRepGProp` contre `signed_volume`, écart relatif 2e-3,
  5e-3 pour les congés et le loft) ;
- le nombre de **solides**.

Côté gbs, la lecture doit être complète : aucun échec, aucune entité non prise
en charge, `check()` valide. Les arêtes dégénérées ne comptent pas : STEP n'en
a pas, et chaque lecteur insère les siennes aux pôles.

## 2. Formes comparées (20, test `tests_step_occt` du job OCCT de la CI)

| Groupe | Formes |
|---|---|
| primitives | boîte décalée, cylindre d'axe incliné, sphère décalée, cône pointu, tronc de cône, tore, quart de cylindre, zone de sphère |
| booléens | boîte moins cylindre (trou traversant), boîte moins sphère (coin sphérique), boîte plus bossage, cylindre percé en travers de son seam |
| congés | boîte dont les 12 arêtes ont un congé (cylindres partiels, coins sphériques) |
| balayages | révolution complète et partielle d'un profil B-spline (surface de révolution), extrusion du même profil (surface d'extrusion) |
| cas durs | loft cercle-carré (surfaces B-spline), tore incliné, sphère percée par ses pôles, révolution partielle d'un cercle (tore partiel) |

## 3. Défauts trouvés et corrigés

La comparaison a révélé quatre défauts, tous corrigés ici.

### 3.1 Sphère complète écrite comme une `vertex_loop`

OCCT écrit une sphère entière comme une face bornée par une seule
`vertex_loop`, le pôle sud, sans arête. Le lecteur bornait la surface autour de
ce seul point et obtenait une calotte minuscule.

**Correction.** Une face **sans boucle d'arêtes** (aucune boucle, ou des
`vertex_loop` seulement) couvre **tout le domaine** de sa surface. C'est une
face naturelle. Le sommet de la `vertex_loop` est effacé, remplacé par ceux de
la face naturelle.

### 3.2 Seam du fichier ailleurs qu'en 0

La PR 4 plaçait la boîte d'une face qui fait le tour complet sur [0, 2π], en
supposant le seam du fichier au paramètre 0. Or un tore OCCT d'axe −z a son seam
en v à 3π/2.

**Correction.** Si la face contient une arête de seam (utilisée deux fois)
**isoparamétrique** dans la direction fermée, la boîte commence à son
paramètre (`align_with_seam_edges`).

### 3.3 Pcurve hors du rectangle paramétrique

Une arête posée sur le bord v = 0 d'une surface B-spline (loft) donne une
pcurve interpolée qui déborde d'environ 1e-13. C'est plus que `knot_eps`
(2e-14) : l'évaluation de la surface lève une exception hors bornes.

**Correction** (`project_pcurve`, palier 1). Les pôles de la pcurve sont ramenés
dans le rectangle. Une B-spline restant dans l'enveloppe convexe de ses pôles,
la pcurve ne peut plus en sortir, et l'écart à l'arête reste vérifié ensuite.

### 3.4 Boucle écrite à rebours de `same_sense`

Pour la révolution partielle d'un cercle, OCCT écrit la boucle de la face
torique à rebours de la règle STEP. Résultat : la face et le disque voisin
utilisent leur cercle commun dans le même sens, et le shell n'est pas
orientable. Le lecteur d'OCCT le répare par ShapeFix.

**Changement de règle** (révise la PR 6, § 3.3) :

- **la géométrie fixe la région** de la face : le builder tourne les boucles
  dans le sens direct en (u, v) ;
- **`same_sense` fixe sa normale** : le `FaceUse` rendu par le lecteur est
  celui de `same_sense`.

Les arêtes partagées sont alors utilisées en sens opposés dès que les normales
du fichier sont cohérentes, même si ses boucles ne le sont pas. Les désaccords
entre boucles et `same_sense` restent comptés (`sense_mismatches`). Cette
règle couvre aussi le cas ambigu de la PR 5 (face bornée par ses seuls seams) :
le sens est celui du fichier.

## 4. Tests

- `gbs-occt/tests/tests_step_occt.cpp` : 5 cas, 20 formes, 220 assertions
  (job OCCT).
- Deux fichiers écrits par OCCT sont versionnés dans `tests/in/step/` pour
  tester les corrections 3.2 à 3.4 sans OCCT (`tests_step_read.files_written_by_occt`) :
  - `occt_partial_torus.stp` (6 ko) : seam en 3π/2, un désaccord de sens,
    volume exact ;
  - `occt_loft.stp` (25 ko) : arêtes sur les bords de surfaces B-spline.

Toutes les suites BREP et STEP passent, ainsi que tout le binaire de tests
OCCT (24 cas).

## 5. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| comparaison sur primitives, booléens, congés | plus balayages (révolution, extrusion) et cas durs (loft, tore incliné, sphère percée, tore partiel) | couvrir les surfaces de révolution, d'extrusion et B-spline |
| fichiers industriels non versionnés (question 14) | deux petits fichiers **produits par OCCT** versionnés | non-régression des corrections hors du job OCCT ; ce ne sont pas des fichiers industriels |
| `FaceUse` donné par le sens des boucles (PR 5 et 6) | donné par `same_sense` | robuste aux boucles mal orientées (§ 3.4) |

## 6. Points à valider

1. **Critères de comparaison** : nombres de faces, d'arêtes non dégénérées et
   de sommets, volumes à 2e-3, `check()` valide.
2. **Face sans boucle d'arêtes = surface entière** (§ 3.1).
3. **Boîte alignée sur l'arête de seam du fichier** (§ 3.2).
4. **Pôles de pcurve ramenés dans le rectangle** (§ 3.3).
5. **`same_sense` fixe la normale de la face**, la géométrie sa région (§ 3.4).
6. **Deux fichiers OCCT versionnés** comme données de test.
