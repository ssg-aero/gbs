# Palier 2 bis — PR 6 : restriction exacte, bilan

| | |
|---|---|
| PR | branche `feat/brep2b-exact-restriction` |
| Document de référence | [brep_palier2bis.md](brep_palier2bis.md), § 8 ; bilan au § 12 |
| Fichiers | `gbs/intersect/surface_surface.h`, `gbs/intersect/intersection_curve.h`, `gbs-brep/pcurve.h`, `gbs-brep/coincident.h`, `tests/tests_intersect_surfaces.cpp`, `docs/sources/design/brep_palier2bis.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Les courbes des formules fermées (droites, cercles, ellipses) sont restreintes
exactement aux boîtes de paramètres des deux faces, sans échantillonnage. La
PR clôt le palier 2 bis par son bilan.

## 2. Changements

`inside_intervals` (`surface_surface.h`), après la revue :

1. la courbe porteuse est coupée là où elle passe par un **sommet de cône ou un
   pôle de sphère** (`singular_points`) : ses pcurves y sautent ; une porteuse
   réduite à ce point (cercle de rayon nul au sommet) est un contact ponctuel ;
2. sur chaque intervalle, pour chaque surface, la pcurve de la porteuse
   (`primitive_pcurve`, paramètres périodiques déroulés) ;
3. pour chaque borne de la boîte, répétée tous les 2π pour un paramètre
   périodique, les **racines de la coordonnée moins la borne**
   (`coordinate_roots`) : pièces de Bézier de cette coordonnée (`CurvePieces`
   en dimension 1) coupées en deux jusqu'à ne plus atteindre la borne ou être
   monotones, puis dichotomie et un dernier pas de sécante (exact sur un
   segment). Un morceau qui revient vers la borne sans la franchir donne une
   racine double à son extremum (contact) ; un intervalle où la coordonnée reste
   sur la borne donne ses deux bouts. L'étendue des pcurves est la boîte de leurs
   pôles, exacte ;
4. les coupes (`struct Cut { t, seam }`) sont triées et fusionnées
   explicitement, une coupe de seam (borne d'un paramètre périodique, ou point
   singulier) le restant ;
5. les morceaux dont le milieu est dans les deux boîtes sont gardés, et rejoints
   à une coupe qui n'est pas un seam ; un tour complet est recollé à ses deux
   bouts sauf s'il y a un seam à l'origine de la porteuse ;
6. `Restriction::accurate` : faux quand la pcurve d'un intervalle qui porte un
   morceau gardé n'atteint pas sa tolérance ; `intersect` rend alors
   `complete = false`.

`PCurveFit` passe dans `gbs/intersect/intersection_curve.h` (`gbs-brep/pcurve.h`
le reprend), et `primitive_pcurve` le rend.

### 2.1 Revue de code

| Commentaire | Suite donnée |
|---|---|
| le raccord sonde à 1e-6 de la longueur, dans la marge de `on_end` : deux morceaux courts de part et d'autre d'un seam intérieur sont rejoints | plus de sonde : les coupes portent un drapeau `seam`, jamais rejointes. Test `short_parts` : sphère au seam tourné à 1 rad et bande de cylindre coaxiale autour de cet angle, dans les deux ordres ; **la version précédente échouait** (pcurve à u = −0,001, hors de la boîte) |
| l'étendue des droites de bornes vient de 65 échantillons plus une marge | boîte des pôles de la pcurve, exacte, sans marge |
| problème de racines en dimension 1 résolu par l'intersection courbe-courbe 2D | `coordinate_roots`, par Bézier en dimension 1 |
| pcurves de la porteuse recalculées sur chaque morceau | essayé : la pcurve d'un segment découpée diffère de la pcurve recalculée de quelques ulp, et le cas `tangent_slivers` d'OCCT en devenait un shell ouvert. Les morceaux gardent leur pcurve recalculée (plus serrée près d'un contact) ; le temps des booléens OCCT revient à celui de master (2,07 s contre 2,11 s, 2,38 s pour la première version de cette PR) |
| `if (!pc)` mort, écart insuffisant ignoré | `Restriction::accurate`, et `complete = false` |
| `std::unique` avec une tolérance | fusion explicite des coupes, drapeau `seam` conservé |
| tests : seam intérieur, tore, cône près du sommet | `short_parts` (seam intérieur), `restriction_torus_cone` : plan par l'axe d'un tore sur un demi-tube (deux demi-méridiens de π), ellipse d'un plan sur un cône près de son sommet, sur une partie du cône |

### 2.2 Défaut trouvé en chemin

`coincident.h` : une arête sans pcurve parce qu'elle sort de la surface
(`OffSurface`) était refusée quand son milieu touchait la face. C'est le cas
d'une arête qui ne fait que toucher la surface, comme l'arête d'une boîte tangente
à un cylindre. Elle n'est plus une coupe. Le cas des copeaux tangents est
ainsi devenu sensible, à la restriction près : un écart de 1,6e-13 sur une
extrémité y ouvrait un shell, ce que la sécante finale a corrigé. Ce point
mérite une revue à part (tolérances des sommets proches d'une tangence).

## 3. Tests

| Test | Cas | Attendu |
|---|---|---|
| `short_parts` | cercle d'un plan sur une face de cylindre large de 0,005 rad | une courbe de 0,005 |
| | même cercle sur des faces de plan en bande y ≥ 1 − 1e-8 qu'il ne fait qu'effleurer (arc de 2,8e-4 rad, 1/40 d'un ancien pas), bande tournée de kπ/8 | une courbe (deux quand la bande est à cheval sur le seam du cylindre), de longueur 2 acos(1 − 1e-8) à 1e-5 près |
| | tour complet sur le cylindre entier | une courbe fermée de longueur 2π |

L'ancien code échoue sur la première bande : il ne trouve pas la courbe.

Les suites `intersect` (dont les paires de primitives existantes), `brep`,
`step`, `iges`, le binaire OCCT (29 cas) et les tests Python passent.

## 4. Coût et précision

- Les booléens comparés à OCCT prennent le même temps que sur master (2,07 s
  contre 2,11 s pour les 66 opérations).
- Sur un croisement rasant (angle θ), l'extrémité d'un morceau glisse le long de
  la courbe de l'écart de la pcurve (tol / 10) divisé par θ. Elle reste sur la
  borne à la tolérance près.

## 5. Points à valider

1. Le bilan du palier 2 bis (§ 12 du document) et ses limites restantes.
2. La sensibilité du cas des copeaux tangents (§ 2.2), à reprendre à part.
