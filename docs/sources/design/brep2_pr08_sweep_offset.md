# Palier 2 — PR 8 : extrusion, révolution, face décalée

| | |
|---|---|
| PR | branche `feat/brep-sweep-offset` |
| Document de référence | [brep_palier2.md](brep_palier2.md), § 7 ; question 9 |
| Fichiers | `gbs-brep/sweep.h`, `gbs-brep/offset.h` (nouveaux), `gbs-brep/brep`, `tests/tests_brep_sweep.cpp` (nouveau), `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Construire des formes par balayage, comme les opérations de base d'un
modeleur :
- une **extrusion** (prisme) ;
- une **révolution**, complète ou partielle ;
- une **face décalée** (offset).

Le profil peut être une arête, un wire ou une face.

## 2. API

```cpp
struct SweepOptions { T tol = brep_default_tolerance; T approx_tol = 1e-5; };
struct SweepResult  { ShapeId shape; size_t approximated; };

auto make_prism(m, ShapeId profile, point<T, 3> d, SweepOptions = {}) -> BuildResult<SweepResult>;
auto make_revolution(m, ShapeId profile, ax1<T, 3> axis, T angle = 2π, SweepOptions = {}) -> BuildResult<SweepResult>;
auto offset_face(m, FaceId f, T d, T tol = …) -> BuildResult<FaceId>;
```

Une arête donne une face, un wire un shell, une face un solide.

## 3. Choix d'architecture et justification

### 3.1 Faces latérales : primitives quand c'est possible

Chaque arête du profil donne une **face naturelle** de sa surface balayée.

| Arête | Extrusion selon d | Révolution autour de l'axis |
|---|---|---|
| segment | **plan** s'il est normal à d, sinon NURBS bilinéaire | dans un plan de l'axe : **cylindre** (parallèle), **cône** (oblique, sommet compris), NURBS pour un segment normal à l'axe (anneau) ; sur l'axe : rien |
| cercle | **cylindre** si d suit son axe, sinon NURBS | dans un plan de l'axe : **sphère** (centre sur l'axe), **tore** (hors de l'axe) |
| B-spline | NURBS exacte (`build_extrusion`) | NURBS exacte (`build_revolution`) |
| primitive non traitée | NURBS exacte de la primitive (`to_nurbs`, restreinte) | idem |
| autre courbe | B-spline interpolée, comptée dans `approximated` | idem |

Les surfaces analytiques restent analytiques. Elles sont écrites telles quelles
en STEP, et les intersections ultérieures (booléens) utilisent leurs formes
fermées.

### 3.2 Assemblage par sewing

Les faces latérales sont des faces **naturelles** : leurs bords sont les isos
de leur surface. Pour un profil face, les **faces d'extrémité** sont des
copies du profil : en place, et déplacée de d, ou tournée de l'angle. Une
révolution complète n'en a pas. Le tout est cousu par le **sewing du palier
1**, qui reconnaît les arêtes communes par leur géométrie, puis fermé par
`make_solid`.

Le document de conception prévoyait des arêtes partagées par construction,
« sans sewing ». Le sewing est plus simple : il traite les seams, les pôles et
les arêtes d'axe sans cas particuliers, et il est éprouvé (palier 1). Il coûte
un appariement de quelques arêtes.

### 3.3 Face décalée

- **Même paramétrage.** La surface décalée garde le paramétrage de la
  surface, si bien que les **pcurves sont reprises telles quelles** (§ 9.1 du
  palier 1).
- **Primitives : décalage exact.** C'est un plan translaté, un cylindre ou une
  sphère de rayon r ± d, ou un tore de petit rayon r ± d. Le bon signe est
  choisi par une vérification S′(u, v) = S(u, v) + d n en quatre points.
- **Cône.** Son décalé est un cône de rayon R + d / cos α, et son paramètre
  v se décale de −d sin α : ses pcurves sont décalées d'autant.
- **Autres surfaces.** `OffsetSurface` (S + d n), avec des dérivées premières
  par différences finies.
- **Arêtes.** Ce sont les pcurves posées sur la surface décalée
  (`CurveOnSurface`), et les sommets sont à leurs bouts. Un **pôle** reste un
  pôle si son décalé est un point (sphère). Le sommet d'un cône devient en
  revanche une arête fermée, car son décalé est un cercle.

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| shell fermé « sans sewing » | faces naturelles cousues | seams, pôles et axes sans cas particuliers |
| face décalée : `SurfaceOffset` | `OffsetSurface` avec dérivées | la `SurfaceOffset` existante n'a pas de dérivées, dont les constructeurs et `check` ont besoin |
| — | sommet de cône décalé en arête fermée | la topologie l'impose |
| balayage le long d'un chemin (question 9) | non fait | conforme à la recommandation |

## 5. Tests (`tests_brep_sweep`, 3 cas, 62 assertions, 0,9 s)

`check()` passe sur tous les résultats. Les volumes sont comparés à 2e-3 près.

| Test | Couvre |
|---|---|
| `prisms` | plaque percée extrudée de 10 mm : 7 faces (6 plans, 1 cylindre), volume = aire × 10 ; en biais (côtés NURBS), même volume ; disque du cylindre extrudé : 100π ; un cercle donne une face cylindrique, un wire un shell de 4 faces ; erreurs |
| `revolutions` | rectangle autour de z : anneau (2 cylindres) de volume 5π, demi-tour 2,5π ; disque R 5 autour d'un axe à 10 de son centre : **tore** de volume 2π²·10·25 ; triangle rectangle sur l'axe : **cône** de volume 4π ; demi-disque sur son diamètre : **sphère** de volume 32π/3 ; erreurs |
| `offset_faces` | les trois faces du cylindre décalées de 1 (primitives exactes, point décalé le long de la normale à 1e-9 près) ; sphère vers l'extérieur et vers l'intérieur (rayons 3,5 et 2) ; faces du cône (pcurves décalées, sommet devenu arête) ; face B-spline (surface décalée générique, points à 1e-12 près) ; erreur |

## 6. Points à valider

1. **Faces latérales en primitives** chaque fois que la géométrie le permet.
2. **Assemblage par le sewing** plutôt que par arêtes partagées construites.
3. **Face décalée au même paramétrage**, pcurves reprises (ou décalées pour le
   cône).
4. **`OffsetSurface`** générique, avec des dérivées par différences finies.
5. **Pas de balayage le long d'un chemin** dans ce palier.
