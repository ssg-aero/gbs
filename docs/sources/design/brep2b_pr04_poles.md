# Palier 2 bis — PR 4 : pôles et sommets

| | |
|---|---|
| PR | branche `feat/brep2b-poles` |
| Document de référence | [brep_palier2bis.md](brep_palier2bis.md), § 6 (révisé par cette PR) |
| Fichiers | `gbs-brep/split.h`, `gbs/intersect/surface_surface.h`, `gbs-brep/classify.h`, `tests/tests_intersect_marching.cpp`, `tests/tests_brep_boolean.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `docs/sources/design/brep_palier2bis.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Les faces de sphères et de cônes sont découpées par des courbes qui passent
par leurs pôles ou leurs sommets. Cela rend possibles les booléens qui
coupent une sphère par ses pôles ou un cône par son sommet.

## 2. Ce qui manquait

Le cheminement passait déjà les pôles : depuis la PR 3, la pose sur le bord
v = ±π/2 se fait par une SVD tronquée, qui supporte le rang déficient du pôle.
Une courbe de Viviani dont la sphère a ses pôles sur la courbe a la bonne
longueur, en branches qui s'arrêtent aux pôles. Ce qui manquait est dans la
découpe et dans les pcurves.

| Fichier | Changement |
|---|---|
| `split.h` | `split_detail::split_degenerate` : une coupe qui aboutit au sommet d'une arête dégénérée (pôle, sommet de cône), en un u qui n'est pas une extrémité de son segment dans (u, v), la découpe en ce u (même sommet, même pcurve, co-arêtes insérées comme `split_edge`). Sans cela, la coupe restait pendante |
| `surface_surface.h` | `primitive_pcurve` : l'angle d'un échantillon pris au sommet d'un cône ou à un pôle est indéfini ; il est extrapolé des deux échantillons voisins. La génératrice qui finit au sommet prenait l'angle de l'autre génératrice et sa pcurve partait en diagonale |

## 3. Autres défauts trouvés par la comparaison avec OCCT

| Défaut | Effet | Correction |
|---|---|---|
| une boucle-trou pouvait être donnée à sa **jumelle**, la boucle extérieure faite des mêmes arêtes (le point testé est sur son contour) | face d'aire nulle, booléen refusé (« a face without inside ») pour une boîte dont une face contient l'axe d'un cylindre : un cas élémentaire | `split_face` exclut la jumelle |
| `classify` (point / face) rejetait des paramètres hors des bornes de la surface **d'un arrondi** (0 contre 2,4e-16 relu d'un fichier) avant de tester le bord | une génératrice posée sur le seam d'un cône n'était pas vue comme arête de bord | marge de 1e-12 de l'étendue |

## 4. Tests

| Test | Cas | Attendu |
|---|---|---|
| `through_poles` (marching) | Viviani, axe de la sphère le long du cylindre : la courbe passe par les deux pôles | non singulier, longueur à 1e-5, quatre extrémités de branche aux pôles |
| `poles` (booléens) | sphère R 3 centrée sur le flanc du cylindre R 5, axes parallèles : ses deux pôles sont sur le cylindre | commun, différences dans les deux sens et fusion aux volumes de référence (tranches : lentille de deux disques) |
| OCCT `booleans` | `viviani_poles` | valides, mêmes volumes et nombres de solides ; plus de faces qu'OCCT (§ 5) |
| | `half_cylinder` (boîte dont une face contient l'axe du cylindre), `cone_apex` (cône coupé par un plan passant par son sommet), `cone_apex_seam` (idem, seam du cône dans le plan) | valides, volumes, **mêmes nombres de solides et de faces** que `BRepAlgoAPI` |

Les suites `brep`, `intersect`, `step`, `iges`, le binaire OCCT (29 cas) et
les tests Python passent.

## 5. Écart avec OCCT

Pour `viviani_poles`, gbs rend plus de faces qu'OCCT : 4 au lieu de 2 à 3.
Dans la différence, la partie de la sphère hors du cylindre recouvre tout son
seam, d'un pôle à l'autre. Dans l'espace (u, v), elle forme deux boucles de
part et d'autre du seam, donc deux faces qui partagent l'arête du seam. Il en
va de même pour la partie du cylindre, dont les deux faces partagent aussi une
arête.

OCCT en fait une face chacune. Il coupe pourtant de la même façon le
demi-cylindre qui recouvre son seam (`half_cylinder`, mêmes nombres de faces).
Les deux représentations sont valides et ont le même volume. Ce cas est
comparé sans le nombre de faces.

## 6. Points à valider

1. La révision du § 6 du document : pas de changement du cheminement.
2. L'écart du § 5 : on l'accepte, ou on fusionne les faces de part et
   d'autre d'un seam qu'elles recouvrent entièrement (pcurves décalées d'une
   période), dans une PR ultérieure.

## 7. Revue de code (a posteriori, lot 4 de la revue du noyau BREP)

Pas de défaut trouvé dans `split_degenerate`, ni dans l'extrapolation de
l'angle aux échantillons singuliers de `primitive_pcurve`. Une sphère coupée par
des plans qui contiennent son axe, à 0, 0,7 et 1,9 rad du seam, donne 3 faces
(deux de part et d'autre du seam, comme au § 5), de volume et d'aires exacts,
et `check()` passe. Les défauts trouvés dans `split.h` sont dans
[brep2_pr07_split.md](brep2_pr07_split.md) § 7.
