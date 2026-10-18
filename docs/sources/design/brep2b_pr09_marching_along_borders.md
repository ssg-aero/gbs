# Palier 2 bis — PR 9 : cheminement le long d'un bord

| | |
|---|---|
| PR | branche `fix/marching-along-borders` |
| Document de référence | [brep2b_pr08_exact_plane_pcurves.md](brep2b_pr08_exact_plane_pcurves.md), § 3 |
| Fichiers | `gbs/intersect/marching.h`, `tests/tests_intersect_marching.cpp`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

La PR 8 notait qu'une courbe d'intersection longeant le bord d'une face B-spline
ressortait du cheminement en une vingtaine de morceaux. C'est le cas du couvercle
plan d'un cylindre B-spline contre son flanc. Le test
`coincident_closed_bsplines` y passait encore 6 s.

## 2. Cause

Pour le couvercle du bas, le cercle est sur le bord v = 0 du flanc, et le
paramètre v des points prédits vaut −3e-34, sous la borne d'un arrondi.

`Marcher::correct` vérifie, avant même de corriger, que le point prédit est sur
la surface. Hors des bornes, il le rend tel quel, **sans le converger** : c'est
le point qui sert à se poser sur le bord. Le cheminement acceptait donc des
points où les deux surfaces étaient à 1e-3 ou 1e-2 l'une de l'autre. Les
branches partaient en biais et s'arrêtaient sur la boîte du plan, que le cercle
touche : 21 courbes fausses, et autant de croisements avec les bords dans la
section.

La PR 5 ramenait déjà sur la borne les dépassements d'arrondi, mais seulement
dans la boucle de Newton, pas dans ce test préalable.

## 3. Correction

Les paramètres au-delà d'une borne de la surface de moins de 1e-10 de son
étendue sont ramenés sur la borne **avant** le test préalable, comme dans la
boucle de Newton.

## 4. Effet

| | Avant | Après |
|---|---|---|
| couvercle du bas contre le flanc | 21 courbes, écarts jusqu'à 2e-2 | 1 courbe fermée |
| `coincident_closed_bsplines` | 6,1 s | 2,3 s |
| booléens de la comparaison OCCT | 2,07 s | 1,89 s |

## 5. Tests

| Test | Cas | master |
|---|---|---|
| `along_a_border` (`tests_intersect_marching`) | flanc B-spline rationnel d'un cylindre (cercle rationnel extrudé de z = 0 à 2) contre le plan z = 0 sur la boîte [−1, 1]² que le cercle touche : une courbe fermée sur les deux surfaces, longueur 2π à 1e-5 | **échoue** (courbes à 6e-3 des surfaces) |

Les 477 tests C++, le binaire OCCT (29 cas) et les 51 tests Python passent.

## 6. Points à valider

1. La marge de 1e-10 de l'étendue, la même que dans la boucle de Newton.
