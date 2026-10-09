# BREP palier 1 — PR 6 : validation minimale (`gbs-brep/check.h`)

| | |
|---|---|
| PR | branche `feat/brep-check` |
| Document de référence | [brep_core.md](brep_core.md), sections 2.6 (invariants) et 5.4 ; notes précédentes [PR 1](brep_pr01_model.md) à [PR 5](brep_pr05_face_wire.md) |
| Fichiers | `gbs-brep/check.h` (365 l., nouveau), `gbs-brep/brep`, `tests/tests_brep_check.cpp` (278 l., nouveau), `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Donner un diagnostic de validité d'un shape : `check(model, shape)` parcourt
toutes les entités sous `shape` et liste **toutes** les violations des
invariants du modèle (§ 2.6 du document), sans s'arrêter à la première.
C'est l'équivalent réduit de `BRepCheck_Analyzer` d'OCCT. Il servira aux tests
de toutes les PR suivantes (sewing, solide, IGES), à l'utilisateur après un
import, et au lecteur STEP.

## 2. Vue d'ensemble

```
check(m, shape, CheckOptions{n_samples = 9, n_uv_polygon = 16, geometry = true}) -> CheckReport
CheckReport { std::vector<CheckEntry> entries; ok(); count(Issue); has(Issue, ShapeId); }
CheckEntry  { ShapeId shape; Issue issue; std::string detail; }

ordre de visite : solides et compounds (références) → sommets → arêtes → faces (et leurs wires) → wires libres → shells
chaque entité une seule fois (explore), modèle jamais modifié, jamais d'exception sur un modèle incohérent
```

| Niveau | Issues | Contrôle |
|---|---|---|
| référence | `DeadReference` | identifiant mort ou invalide dans une entité (non suivi ensuite) |
| sommet | `InvalidTolerance`, `InvalidPoint` | tolérance > 0 et finie, point fini |
| arête | `InvalidTolerance`, `EdgeBoundsInverted`, `EdgeWithoutCurve`, `InvalidDegenerateEdge`, `VertexOffCurve`, `ToleranceOrderViolated` | `u1 < u2` ; courbe présente sauf dégénérée ; dégénérée ⇔ pas de courbe et `v1 == v2` ; extrémités de courbe dans la boule des sommets ; `vertex.tol ≥ edge.tol` |
| wire | `WireNotChained`, `WireNotClosed` | chaînage ; fermeture si le drapeau `closed` est levé **ou** si le wire borde une face |
| face | `InvalidTolerance`, `FaceWithoutSurface`, `FaceWithoutBoundary`, `CoEdgeWithoutPCurve`, `PCurveRangeMismatch`, `SameParameterViolated`, `ToleranceOrderViolated` | surface et au moins un wire ; pcurve sur chaque co-arête, couvrant `[u1, u2]` ; `S(pcurve(t))` à `edge.tol` de `C(t)` sur 9 paramètres ; `edge.tol ≥ face.tol` |
| face, (u,v) | `OuterWireNotDirect`, `InnerWireNotIndirect`, `InnerWireOutsideOuter`, `InnerWiresOverlap` | aire signée des polygones `(u, v)` ; inclusion et chevauchement par test pair-impair |
| shell | `NonManifoldEdge`, `ShellNotOrientable`, `ShellNotClosed` | > 2 usages ; deux usages de même sens effectif ; fermeture si drapeau `closed` **ou** shell d'un solide |

## 3. Choix d'architecture et justification

### 3.1 Un rapport complet plutôt qu'un booléen

Le rapport liste chaque violation avec l'entité en cause et un détail lisible
(distance mesurée, nombre d'usages, arête concernée). On peut ainsi
tester précisément une violation (`has(Issue, ShapeId)`), compter
(`count(Issue)`), ou afficher le tout. `CheckReport` n'est pas template :
il ne contient que des identifiants et des chaînes.

### 3.2 Robuste à un modèle incohérent

`check` doit pouvoir diagnostiquer un modèle cassé, en cours de mutation
ou lu depuis un fichier douteux. Il teste donc chaque référence avant de
la suivre (`DeadReference` puis arrêt sur cette branche) et ne lève jamais
d'exception de lui-même. L'explorateur de la PR 2 ignore déjà les entités
mortes ; `check` les rapporte au niveau de l'entité qui les référence.

### 3.3 Exigences contextuelles

Un wire libre ouvert ou un shell libre ouvert sont **valides** : ce sont des
états intermédiaires normaux (wire en construction, shell avant sewing). Ils
deviennent invalides seulement si leur drapeau `closed` le prétend, ou si leur
rôle l'exige (wire bordant une face, shell bordant un solide). `check`
détermine ce rôle d'après le shape racine.

### 3.4 Contrôles géométriques par échantillonnage, désactivables

`SameParameter` est vérifié sur 9 paramètres par co-arête, l'orientation et
l'inclusion sur des polygones de 16 points par co-arête : rapide, et
suffisant pour les erreurs grossières que vise un « lite ». Une petite marge
relative (1e-9) absorbe les arrondis. `CheckOptions::geometry = false` ne garde
que les contrôles topologiques (aucune évaluation de courbe ni de surface).

### 3.5 Une seule source pour les règles de shell

`NonManifoldEdge` et `ShellNotOrientable` sont tirés de `shell_edge_uses` (PR 2),
comme `is_manifold` et `is_orientable` : impossible que `check` et les
requêtes divergent. Les arêtes dégénérées sont ignorées de la même façon.

## 4. Écarts par rapport au document de conception

| Document (§ 5.4) | Implémentation | Raison |
|---|---|---|
| `CheckReport` = paires (shape, issue) | entrées (shape, issue, détail) + `count`, `has` | diagnostic lisible, tests ciblés |
| `SolidNotOutward` | **pas encore** : ajouté en PR 8 avec `signed_volume` | le calcul de volume signé appartient à la PR du solide |
| — | `DeadReference`, `InvalidPoint`, `InvalidDegenerateEdge`, `FaceWithoutSurface`, `FaceWithoutBoundary`, `InnerWiresOverlap` | invariants du § 2.6 non listés explicitement en § 5.4 |
| `CheckReport<T>` | `CheckReport` non template | il ne contient pas de valeur de type `T` |

Hors périmètre, comme prévu : auto-intersection d'un wire en `(u, v)`,
intersection entre faces, validité des surfaces elles-mêmes.

## 5. Tests (`tests_brep_check`, 7 cas, 40 assertions)

| Test | Couvre |
|---|---|
| `builders_produce_valid_shapes` | **toutes les sorties des PR 1 à 5 passent** : boîte faite à la main, faces naturelles (cylindre analytique et NURBS, sphère, tore, cône, plan, révolution), shells fermés d'une face (sphère, tore), disque projeté, face trouée depuis des boucles 2D, wire libre ouvert, le tout dans un compound ; idem sans géométrie |
| `vertex_and_edge_issues` | tolérance nulle, point NaN, courbe absente, bornes inversées, sommet déplacé hors des courbes, arête dégénérée entre deux sommets, sommet effacé ⇒ 3 `DeadReference` |
| `tolerance_order` | arête plus tolérante que ses sommets, face plus tolérante que ses arêtes : 6 violations exactement |
| `wire_issues` | co-arêtes permutées, contour de face ouvert (et shell ouvert qui en découle), wire libre ouvert valide sauf s'il se dit fermé |
| `face_issues` | surface absente, aucun wire, pcurve absente, pcurve trop courte, pcurve décalée ⇒ `SameParameterViolated`, contour extérieur retourné (et shell non orientable qui en découle) |
| `hole_issues` | trou retourné, deux trous qui se chevauchent, trou hors du contour |
| `shell_issues` | septième face ⇒ 4 arêtes non variétés ; shell ouvert : invalide sous un solide ou s'il se dit fermé, valide sinon ; face retournée ⇒ 4 arêtes non orientables ; face, shell et membre de compound morts ; racine morte |

## 6. Points à valider

1. **Rapport complet** (toutes les violations, avec détail) plutôt qu'un arrêt
   à la première.
2. **Exigences contextuelles** : wire et shell libres ouverts valides, sauf
   drapeau `closed` ou rôle (bord de face, bord de solide).
3. **Contrôles géométriques par échantillonnage** (9 paramètres pour
   SameParameter, 16 points par co-arête pour les polygones), désactivables.
4. **`SolidNotOutward` reporté en PR 8**, avec `signed_volume`.
5. **`check` ne lève jamais d'exception**, y compris sur un modèle incohérent.
