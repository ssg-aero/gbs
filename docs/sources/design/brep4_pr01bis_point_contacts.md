# Palier 4 — PR 1 bis : contacts ponctuels sur des surfaces B-splines

| | |
|---|---|
| PR | branche `feat/brep-bspline-point-contact` |
| Document de référence | [brep_palier4.md](brep_palier4.md), § 0 quater |
| Fichiers | `gbs/intersect/marching.h`, `tests/tests_intersect_marching.cpp`, `gbs-occt/tests/tests_brep2_occt.cpp`, `python/tests/test_examples.py`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier4.md`, `News.md`, cette note |
| Statut | À valider |

## 1. Objectif

Une sphère NURBS posée sur une boîte donnait `NonGeneric` sur les trois
booléens, alors que la même sphère exacte passait. Entre primitives, le point de
contact vient des formes fermées ; entre B-splines, il passe par le cheminement.
Un modèle converti en NURBS (un STEP d'un autre système, la comparaison à OCCT)
échouait donc là où le même modèle en primitives passait. C'était la dernière
limite connue de la galerie (`nurbs_sphere_touching`, en `xfail`).

## 2. Diagnostic

Le cheminement **trouvait** le point de contact, à 4·10⁻⁸ près : la recherche des
départs y converge (Gauss-Newton sur S₁ = S₂). Mais au point de tangence les
normales sont parallèles, la marche lancée de ce départ s'arrête aussitôt
(`End::singular`), et le résultat entier était déclaré singulier. Le booléen
refuse une intersection singulière.

Le point d'arrêt peut être de deux natures :
- un **contact isolé** : les surfaces se touchent sans se traverser (une sphère sur
  un plan, dans une cuvette) ; rien n'est singulier ;
- un **point de branchement** : deux courbes d'intersection s'y croisent (une selle
  sur son plan tangent) ; le refus est justifié.

Le commentaire de `tangent_branch` notait déjà la différence (« both are of the
same size at an isolated contact […], which is not traced »), sans en tirer de
conclusion.

## 3. Contenu

### 3.1 Certification par les secondes formes (`TangentMarcher::isolated_contact`)

Près du contact, l'écart entre les surfaces le long de la normale commune vaut
(l₁ s₁² + l₂ s₂²) / 2, où l₁ et l₂ sont les valeurs propres de II₁ − II₂, la
différence des secondes formes fondamentales dans le plan tangent commun.

- **Définie** (l₁ et l₂ de même signe) : l'écart croît dans toutes les directions,
  le contact est isolé. C'est la certification.
- **Indéfinie** : deux courbes se croisent, point de branchement ; reste singulier.
- **Presque dégénérée** (|l₁| ≤ `angular_tol` |l₂|) : c'est le domaine de
  `tangent_branch`, qui trace une courbe de tangence. Les deux tests sont
  complémentaires, avec le même seuil.
- **Nulle** (|l₂| sous le bruit des courbures, comme `direction`) : contact
  dégénéré, reste singulier.

Le point est ensuite mis sur la tangence par Gauss-Newton sur les quatre
équations de `TangentMarcher::correct` (pieds et normales parallèles). Elles sont
de rang 4 en un contact isolé : la convergence est quadratique, et le point
tombe à 10⁻¹⁵ près au lieu de 4·10⁻⁸. La forme doit rester définie au point
corrigé, et celui-ci dans les boîtes.

`direction` et `isolated_contact` partagent maintenant le calcul de II₁ − II₂
(`form_difference`) ; `direction` est inchangée.

### 3.2 Au pôle d'une surface

Une sphère posée sur un plan, son axe vertical, touche par son pôle. Là, S_u = 0 :
la seconde forme n'est pas définie dans les paramètres, et la certification
échouait.

- II est une propriété de la surface, continue à travers un pôle régulier. Une
  surface dont les paramètres dégénèrent au point (plus petite valeur singulière
  de (S_u, S_v) sous 10⁻⁶ de la plus grande) a sa forme prise à des paramètres
  rapprochés du milieu de sa boîte de 10⁻⁴ de leur étendue.
- Elle y est prise sur **sa propre normale**, orientée comme la normale commune.
  Avec la normale de l'autre surface, l'écart d'angle entre les deux normales
  (de l'ordre de 10⁻⁴) est justement ce qui porte la courbure le long de u, où
  |S_u| est du même ordre : la forme sortait nulle.
- La base du plan tangent est prise sur S₁ quand S₂ dégénère.
- Si Newton échoue au pôle (les équations utilisent la normale par les
  paramètres), le départ est gardé tel quel, à √ε de la taille des surfaces près,
  à condition que l'écart y soit dans la tolérance. Dans les cas testés, Newton
  converge.

### 3.3 Cheminement

Dans la boucle finale sur les contacts, un départ arrêté (`stopped`) et certifié
devient un point de contact non singulier, à sa position corrigée. Sans dérivées
secondes (une surface décalée, par exemple), rien ne change : il reste
singulier.

## 4. Résultats

| | Avant | Après |
|---|---|---|
| `nurbs_sphere_touching` (galerie) | `NonGeneric` sur fuse, cut, common | 7, 6 et 0 faces, comme la sphère exacte ; valide, maillage fermé |
| OCCT : `sphere_touching_nurbs` (boîte et sphère NURBS) | — | mêmes volumes, faces et solides qu'OCCT |
| OCCT : `plane_touching_nurbs_sphere` (boîte exacte, sphère NURBS) | — | idem |
| précision du point de contact | 4·10⁻⁸ | 10⁻¹⁵ |

## 5. Tests

- `tests_intersect_marching.isolated_contacts` (nouveau) : sphère NURBS sur un plan et
  sur un carreau plat B-spline, sphère NURBS dans une cuvette, cuvette dans une
  cuvette plus large. Axes de la sphère (1, 1, 1), x et z : contact hors pôle et
  au pôle. Les deux ordres des surfaces. Un point, aucune courbe, ni singulier ni
  incomplet, à 10⁻⁹ près sur les deux surfaces.
- `free_surfaces` : la sphère sur un plan par cheminement, à 10⁻⁹ près au lieu de
  10⁻³, non singulière.
- `tangency_limits` : commentaire mis à jour (le contact isolé n'est plus singulier).
- Les cas singuliers existants restent singuliers : surface sans dérivées secondes,
  tangence qui finit en se croisant (`tangency_limits`).
- Comparaison OCCT : deux paires ajoutées (§ 4).
- Python : `nurbs_sphere_touching` n'est plus une limite connue ; nouveau
  `test_nurbs_contact_like_primitives` (mêmes nombres de solides et de faces que la
  sphère exacte, pour les trois opérations).

Les 737 tests de `ctest` passent, comparaison à OCCT comprise. Les 97 tests Python
passent, sans `xfail`.

## 6. Limites

| Limite | Suite |
|---|---|
| Un contact dégénéré (II₁ − II₂ nulle : z = x⁴ + y⁴ sur son plan tangent) reste singulier | ordre supérieur de l'écart ; rare en conception |
| Un contact presque dégénéré (|l₁| ≤ `angular_tol` |l₂|) passe par le tracé d'une courbe de tangence, comme avant | — |
| Sans dérivées secondes (surface décalée), le contact reste singulier | dérivées secondes exactes de `OffsetSurface` (PR 3 du palier) |

## 7. Points à valider

1. **Certification par II₁ − II₂ définie**, avec le seuil `angular_tol` de
   `tangent_branch` pour séparer contact isolé et courbe de tangence.
2. **Au pôle, formes prises à 10⁻⁴ de l'étendue des paramètres**, sur la normale
   propre de la surface.
3. **Départ gardé** quand Newton ne peut pas être appliqué au pôle, l'écart y étant
   dans la tolérance.
