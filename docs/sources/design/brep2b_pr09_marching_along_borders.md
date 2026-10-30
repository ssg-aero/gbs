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

Les paramètres au-delà d'une borne de la surface de moins de la marge d'arrondi
sont ramenés sur la borne **avant** le test préalable, comme dans la boucle de
Newton (`Marcher::snap`, une seule définition pour les deux).

La marge est celle de `inside` et `out_of` : 1e-10 (1 + |a| + |b|) pour la boîte
[a, b] du paramètre. Un point laissé au-delà d'une borne de la surface est donc
hors de la boîte, et le cheminement se pose sur son bord ; avec une marge
différente, une bande de dépassements restait acceptée sans être convergée.

## 4. Effet

| | Avant | Après |
|---|---|---|
| couvercle du bas contre le flanc | 21 courbes, écarts jusqu'à 2e-2 | 1 courbe fermée |
| `coincident_closed_bsplines` | 6,1 s | 2,3 s |
| booléens de la comparaison OCCT | 2,07 s | 1,89 s |

## 5. Tests

| Test | Cas | master |
|---|---|---|
| `along_a_border` (`tests_intersect_marching`) | flanc d'un cylindre extrudé de z = 0 à 2, contre les plans z = 0 et z = 2 (bords inférieur et supérieur) : flanc B-spline rationnel (cercle rationnel) sur la boîte [−1, 1]² que le cercle touche, une courbe de longueur 2π à 1e-5 ; sur une boîte plus large, une courbe **fermée** de longueur 2π à 1e-6 ; flanc B-spline non rationnel (cercle interpolé) : une courbe fermée de la longueur de la directrice ; plans à z = −1e-6, −1e-9 et 2 + 1e-9, au-delà de la marge : aucune courbe ; plan y = 0 par l'axe : la couture (bord en u) et la génératrice opposée, longueur totale 4 | **échoue** (courbes à 6e-3 des surfaces) |

Les 477 tests C++, le binaire OCCT (29 cas) et les 51 tests Python passent.

### 5.1 Revue de code

| Commentaire | Suite donnée |
|---|---|
| marge du recalage (1e-10 de l'étendue) différente de celle de `in_box` (1e-10 (1 + \|a\| + \|b\|)) : entre les deux, un point non convergé était encore accepté ; une étendue de 1e6 recalait à 1e-4 | une seule marge, `Marcher::margin`, celle de `in_box`, pour le recalage et `out_of` |
| recalage écrit deux fois | `Marcher::snap`, appelé avant le test préalable et après chaque pas de Newton |
| `slo` / `shi` appellent `bounds()` (virtuelle) dans la boucle intérieure | bornes des surfaces gardées dans le `Marcher` (`sb1_`, `sb2_`) |
| tests limités au bord inférieur d'un flanc rationnel | bord supérieur, bord en u, flanc non rationnel, plans au-delà de la marge, drapeau `closed` (ci-dessus) |

## 6. Points à valider

1. La marge d'arrondi, 1e-10 (1 + |a| + |b|) pour la boîte [a, b], commune au
   recalage et au test de sortie de boîte.
2. **Limite constatée** : quand la courbe touche le bord d'une boîte (le cercle
   sur la boîte [−1, 1]² du plan, en x = 1), le cheminement s'y pose, et ce
   point posé à une tangence est à quelques 1e-6 de la couture : la courbe a
   la bonne longueur mais n'est pas marquée fermée. Les points restent sur les
   deux surfaces à la tolérance.
