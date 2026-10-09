# BREP palier 1 — PR 8 : solide, volume signé, compound

| | |
|---|---|
| PR | branche `feat/brep-solid` |
| Document de référence | [brep_core.md](brep_core.md), section 6.4, question 15 ; notes précédentes [PR 1](brep_pr01_model.md) à [PR 7](brep_pr07_sew.md) |
| Fichiers | `gbs-brep/closure.h` (+ `face_volume_integral`, `signed_volume`), `gbs-brep/builders.h` (+ `make_solid`, `make_compound`, trois codes d'erreur), `gbs-brep/check.h` (+ `SolidNotOutward`, `VoidNotInward`), `tests/tests_brep_solid.cpp` (180 l., nouveau), `tests/tests_brep_geom.h` (`quad`, `box_corner`, `box_faces` déplacés depuis le test de sewing, taille paramétrable), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Fermer la chaîne du palier 1 : faces → sewing → **solide**. Un solide est un
shell extérieur fermé, orienté normales vers l'extérieur, et des cavités
fermées orientées normales vers la cavité. Le sens global, que le sewing
(PR 7) ne pouvait pas décider, est tranché ici par le **volume signé**. La PR
ajoute aussi `make_compound`, et le contrôle `SolidNotOutward` reporté de la
PR 6.

## 2. Vue d'ensemble

```
face_volume_integral(m, face, n = 64) = 1/3 ∬_région S · (Su × Sv) du dv      (région = contour − trous, en (u,v))
signed_volume(m, shell, n = 64)       = Σ ± face_volume_integral               (signe du FaceUse)

make_solid(m, outer, voids = {}) -> BuildResult<SolidId>
  chaque shell : vivant, fermé (is_closed), orienté (is_orientable), volume non nul
  outer : volume < 0 ⇒ tous ses FaceUse retournés ; cavité : volume > 0 ⇒ retournés
make_compound(m, shapes) -> BuildResult<CompoundId>

check : SolidNotOutward (outer de volume ≤ 0), VoidNotInward (cavité de volume ≥ 0)
```

## 3. Choix d'architecture et justification

### 3.1 Volume par le théorème de la divergence, intégré en (u, v)

`V = 1/3 ∮ S · n dA`, et pour une face paramétrée `n dA = (Su × Sv) du dv`. On
intègre par la règle du point milieu sur une grille `n × n` du rectangle
englobant le contour extérieur en `(u, v)` ; une cellule compte si son centre
est dans la région (dans le contour, hors des trous, test pair-impair sur les
polygones de la PR 5). Pour une face à bornes naturelles, la région est le
rectangle entier et la règle du point milieu est d'ordre 2 ; pour une face
découpée, l'appartenance par centre de cellule donne une erreur d'ordre 1.

Précision **mesurée** à `n = 64` (défaut) :

| Solide | Erreur relative |
|---|---|
| boîte (faces planes) | < 1e-12 (intégrande linéaire : point milieu exact) |
| tore (une face naturelle) | 1e-10 |
| sphère (une face naturelle, pôles) | 1e-4 |
| cylindre fermé par deux disques découpés | 1e-3 (1e-2 à n = 16, 3e-5 à n = 256) |

Le polygone du contour a `4n` points par co-arête, pour que son erreur reste
sous celle de la grille. Une première version le fixait à 32 points : l'erreur
du cylindre plafonnait alors à 2e-3 quel que soit `n`, ce que la mesure a
révélé. La fonction sert à **décider un signe**, pas à faire de la métrologie ;
une intégration exacte sur le contour (Green) pourra venir plus tard si on
veut des propriétés de masse.

Comparaison : OCCT oriente un solide par `BRepLib::OrientClosedSolid`, qui
classe un point à l'infini ; le volume signé est plus simple à écrire avec ce
que gbs offre et robuste sur des faces découpées (document § 6.4).

### 3.2 `make_solid` retourne des shells entiers

Un shell orienté de façon cohérente (`is_orientable`) a deux sens possibles :
on retourne **tous** ses `FaceUse` si le volume n'a pas le signe voulu. On ne
retourne jamais une face isolée, ni la géométrie. Un shell non orientable (une
face retournée, Möbius) est refusé : c'est au sewing ou à l'utilisateur de le
rendre cohérent. Le drapeau `closed` du shell est mis à vrai.

### 3.3 Cavités : stockées et orientées, inclusion non vérifiée

Conformément à la question 15, `make_solid(outer, voids)` vérifie que chaque
cavité est fermée et orientable, la tourne vers l'intérieur, mais ne vérifie pas
qu'elle est **dans** le shell extérieur ni qu'elle ne coupe pas une autre
cavité : il faut pour cela classer un point par rapport à un solide, ce qui
viendra avec les intersections du palier 2.

### 3.4 Garantie forte

Toutes les vérifications et tous les volumes sont calculés avant la première
écriture ; un échec (shell ouvert, non orientable, mort, en double, volume nul)
laisse le modèle inchangé, y compris le sens des `FaceUse` du shell extérieur
quand c'est une cavité qui pose problème.

### 3.5 `check` : orientation des solides

`check` calcule le volume signé des shells d'un solide **seulement** s'ils sont
intègres (vivants, fermés, orientables, toutes les pcurves présentes) ; sinon
les autres issues (`ShellNotClosed`, `ShellNotOrientable`, `CoEdgeWithoutPCurve`…)
sont déjà rapportées et le volume n'aurait pas de sens. Ce contrôle est
géométrique : il est désactivé par `CheckOptions::geometry = false`.

### 3.6 `make_compound`

Valide que chaque membre est vivant ; un compound peut contenir n'importe quel
shape, y compris d'autres compounds, ou rien. Pas de dédoublonnage : un
compound est une collection, l'explorateur dédoublonne au parcours.

## 4. Écarts par rapport au document de conception

| Document (§ 6.4) | Implémentation | Raison |
|---|---|---|
| `make_solid` lève si le shell n'est pas fermé ou pas orientable | `std::expected` avec `ShellNotClosed`, `ShellNotOrientable`, `ZeroVolume` | convention des builders depuis la PR 3 |
| `signed_volume(m, shell, n_u, n_v)` | un seul `n`, plus `face_volume_integral` exposée | grille carrée en `(u, v)` ; l'intégrale par face resservira (propriétés de masse) |
| `SolidNotOutward` dans `check` | idem, plus `VoidNotInward` | symétrie avec les cavités |

## 5. Tests (`tests_brep_solid`, 6 cas, 43 assertions)

| Test | Couvre |
|---|---|
| `box_from_faces_to_solid` | six faces indépendantes → `sew` → `make_solid` : volume 1 à 1e-12 après retournement éventuel, `check` valide ; la boîte faite à la main des PR 1-2 est déjà orientée vers l'extérieur |
| `inward_shell_is_turned_and_checked` | shell retourné : `check` signale `SolidNotOutward` (et rien sans géométrie) ; `make_solid` le remet à l'endroit |
| `curved_solids_volumes` | sphère (r = 2) et tore (R = 3, r = 1) d'une seule face, cylindre fermé par deux disques cousus : volumes à 1e-3 (3e-3 pour le cylindre), `check` valide |
| `hollow_box_with_a_cavity` | boîte de côté 3 avec une cavité de côté 1 : +27 et −1, total 26, `check` valide ; cavité retournée ⇒ `VoidNotInward` |
| `errors_leave_model_unchanged` | shell ouvert, shell non orientable, shell mort, shell en double, cavité ouverte avec un extérieur valide : erreurs, sens des `FaceUse` inchangés, aucun solide créé |
| `compound` | compound imbriqué solide + face libre : 7 faces, 1 solide, 2 compounds, `check` valide ; membre mort refusé ; compound vide accepté |

Les huit suites BREP passent.

## 6. Points à valider

1. **Volume signé par intégration en (u, v)** (point milieu, appartenance par
   centre de cellule) comme critère d'orientation, avec la précision mesurée
   ci-dessus ; pas de classifieur de point.
2. **Retournement de shells entiers** dans `make_solid`, jamais de faces isolées ;
   shell non orientable refusé.
3. **Cavités** orientées vers l'intérieur mais inclusion non vérifiée (question 15).
4. **`check` n'intègre le volume que sur des shells intègres.**
5. **Compound** sans dédoublonnage, vide autorisé.
