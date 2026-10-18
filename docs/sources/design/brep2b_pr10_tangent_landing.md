# Palier 2 bis — PR 10 : pose sur un bord touché tangentiellement

| | |
|---|---|
| PR | branche `fix/tangent-landing` |
| Document de référence | [brep2b_pr09_marching_along_borders.md](brep2b_pr09_marching_along_borders.md), point à valider 2 |
| Fichiers | `gbs/intersect/marching.h`, `tests/tests_intersect_marching.cpp`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

La PR 9 notait une limite : le cercle d'un couvercle, sur la boîte [−1, 1]² de
son disque plan, sortait du cheminement avec une longueur juste mais des bouts
écartés de 8e-6, donc sans être marqué fermé.

## 2. Cause

Le cercle touche la boîte du plan en x = 1, au point même où il traverse la
couture du flanc. Le pas qui franchit ce point sort du flanc (u < 0) ; le point
prédit, non convergé, sort aussi de la boîte du plan (x > 1). `Marcher::land` se
pose sur le premier paramètre à sortir selon l'interpolation linéaire, ici x :
il fixe x = 1 et corrige.

Mais la courbe ne fait que **toucher** x = 1 : elle y est tangente. La
correction à x fixé est mal conditionnée (une racine double), et le point posé
est à y ≈ 4,6e-6, sur les deux surfaces à 1e-11 près, mais avant la couture. Les
deux bouts de la branche sont à 8e-6 l'un de l'autre, plus que la tolérance de
fermeture.

## 3. Correction

`land` essaie les paramètres sortis dans l'ordre de leur rapport de sortie. Un
point posé où la courbe est tangente à l'iso du paramètre fixé n'est pas une
sortie : il est passé, et le paramètre suivant est essayé (ici u du flanc, posé
exactement sur la couture). S'ils sont tous tangents, le premier est gardé,
comme avant.

Le critère, `Marcher::crosses` : le sinus de l'angle entre la tangente t et
l'iso du paramètre k, |t × S_o| / |S_o| où S_o est la dérivée selon l'autre
paramètre de la même surface, est au-dessus de `angular_tol` (1e-3, le sinus
sous lequel un croisement est tangent). En un point singulier, la tangente est
inconnue : c'est un croisement. Un point posé hors des boîtes (pose hors d'une
surface) est rendu tel quel, l'appelant le refuse comme avant.

## 4. Effet

| | Avant | Après |
|---|---|---|
| cercle du couvercle sur la boîte de son disque | bouts à 8e-6, non fermé | bouts confondus (1e-12), fermé |
| `coincident_closed_bsplines` | 2,3 s | 2,2 s |

## 5. Tests

| Test | Cas | master |
|---|---|---|
| `along_a_border` (`tests_intersect_marching`) | flanc rationnel contre z = 0 et z = 2, sur la boîte [−1, 1]² que le cercle touche comme sur une boîte plus large : une courbe fermée, bouts confondus à 1e-12, longueur 2π à 1e-6 | **échoue** sur la boîte touchée (bouts à 8,2e-6, non fermée) |

Les 477 tests C++, le binaire OCCT (29 cas) et les 51 tests Python passent.

## 6. Points à valider

1. Le seuil de tangence : `angular_tol`, le sinus déjà utilisé pour dire un
   croisement tangent.
2. Une courbe tangente à tous les bords qu'elle franchit garde la première pose,
   mal conditionnée.
