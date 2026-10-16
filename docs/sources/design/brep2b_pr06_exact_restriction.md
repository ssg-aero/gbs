# Palier 2 bis — PR 6 : restriction exacte, bilan

| | |
|---|---|
| PR | branche `feat/brep2b-exact-restriction` |
| Document de référence | [brep_palier2bis.md](brep_palier2bis.md), § 8 ; bilan au § 12 |
| Fichiers | `gbs/intersect/surface_surface.h`, `tests/tests_intersect_surfaces.cpp`, `docs/sources/design/brep_palier2bis.md`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Les courbes des formules fermées (droites, cercles, ellipses) sont restreintes
exactement aux boîtes de paramètres des deux faces, sans échantillonnage. La
PR clôt le palier 2 bis par son bilan.

## 2. Changements

`inside_intervals` (`surface_surface.h`) :

1. pour chaque surface, la pcurve de la courbe porteuse sur tout son domaine
   (`primitive_pcurve`, paramètres périodiques déroulés, exacte quand elle est
   linéaire) ;
2. intersection courbe-courbe de cette pcurve avec les droites des bornes de la
   boîte, répétées tous les 2π pour un paramètre périodique (toutes les
   solutions, recouvrements compris) : sorties de boîte et passages de seam ;
3. morceaux entre ces paramètres gardés quand leur milieu est dans les deux
   boîtes, rejoints là où aucun paramètre périodique ne saute ;
4. comme avant : morceaux réduits à un point écartés, tour complet coupé par une
   boîte recollé à ses deux bouts.

`primitive_pcurve` passe avant `inside_intervals`, qui s'en sert.

L'ancien échantillonnage de 512 pas pouvait manquer un morceau plus court qu'un
pas, et ne situait les sorties qu'à la précision de la dichotomie.

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

- Les booléens comparés à OCCT prennent 13 % de temps en plus (2,38 s contre
  2,11 s pour les 66 opérations) : une pcurve sur tout le domaine et quelques
  intersections courbe-droite par courbe.
- Sur un croisement rasant (angle θ), l'extrémité d'un morceau glisse le long de
  la courbe de l'écart de la pcurve (tol / 10) divisé par θ. Elle reste sur la
  borne à la tolérance près.

## 5. Points à valider

1. Le bilan du palier 2 bis (§ 12 du document) et ses limites restantes.
2. Le surcoût de 13 %.
