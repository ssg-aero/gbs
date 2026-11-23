# Palier 7 — PR 3 : section à courbure continue (G2) et bilan du palier

| | |
|---|---|
| PR | branche `feat/brep-g2-sections` (sur master) |
| Document de référence | [brep_palier7_blend_sections.md](brep_palier7_blend_sections.md), § 4.4 et § 4.5 ; son bilan au § 9 |
| Fichiers | `gbs-brep/fillet.h`, `gbs-brep/blend.h` (une ligne), `tests/tests_brep_g2_sections.cpp` (nouveau), `python/gbsBindBrepOps.cpp`, `python/tests/test_brep_ops.py`, `docs/sources/examples/gbs_examples.py`, `docs/sources/design/brep_palier7_blend_sections.md` (bilan), `brep7_pr02_conic_chordal.md` (statut), `News.md`, cette note |
| Statut | À valider (§ 9) ; revue de code : § 10 |

## 1. Objectif

- **Section G2** : un congé dont la section est la quartique de Bézier tangente
  aux faces et de même courbure qu'elles à ses contacts. La courbure ne saute
  plus au bord du congé : les reflets restent continus.
- **Règle C** (validée au § 8 du palier) : la forme la plus lisse, avec la
  même flèche que l'arc de rayon r. Le congé G2 enlève autant de matière au
  milieu que le congé circulaire, et ses contacts sont plus écartés.
- **Toute la topologie** des sections précédentes : rayon constant ou en loi,
  chaînes, raccords de sommet aux coins.
- **Le bilan du palier 7**.

## 2. Sources

- L. Piegl, W. Tiller, *The NURBS Book*, § 2.3 : courbure d'une Bézier à son
  extrémité, κ = (n − 1)/n · h/a².
- G. Farin, *Curves and Surfaces for CAGD*, § 23 : courbes et surfaces
  équitables (*fairing*), l'énergie ∫ κ′² ds.
- Théorème de Meusnier : la courbure d'une section plane d'une surface vaut sa
  courbure normale divisée par le cosinus de l'angle entre la normale de la
  surface et le plan.
- Le comportement de SolidWorks (« Curvature continuous ») et de CATIA (section
  G2). OCCT n'a pas de congé G2 : pas de comparaison.

## 3. API

```cpp
enum class SectionKind { Circular, Conic, G2 };
BlendSection<T>::g2();      // FilletOptions<T>::section
```

`fillet`, `fillet_geometry` (rayon constant ou loi) prennent la section dans les
options. En Python : `brep.SectionKind.G2` et `brep.BlendSection.g2()`.

## 4. Démarche

### 4.1 La section d'une station

À chaque station, la boule de rayon r donne le plan de l'arc (normale n), le
coin S_b des tangentes de l'arc à ses contacts et leur distance L à S_b. La
section G2 se construit dans ce plan (`g2_station`) :

1. **Contacts.** Sur chaque face, le point du plan dont la projection sur la
   tangente de l'arc est à m L de S_b : deux équations, Newton sur (u, v).
   Entre deux plans, c'est le point à m L de l'arête le long de la face.
2. **Tangentes et courbures.** En chaque contact, la tangente de la face dans le
   plan (vers S_b) et le vecteur courbure de sa section par le plan (Meusnier :
   la courbure normale le long de la tangente, divisée par le cosinus de la
   normale de la face au plan).
3. **Pôles.** P₁, Q₁, Q₂, Q₃, P₂ : Q₁ et Q₃ sur les tangentes, à λ du chemin
   vers leur coin S. Q₂ est à la hauteur h = 4/3 · κ · a² au-dessus de chaque
   tangente, où a = |P Q₁|. Entre deux plans, Q₂ = S.
4. **λ** minimise ∫ κ′(s)² ds (Gauss-Legendre, 8 points sur 16 intervalles).
   Brent donne λ à √ε près, puis Newton sur la dérivée de l'énergie à 5 points
   le donne à 10⁻¹³ près, et de façon lisse d'une station à l'autre. La peau
   dérive les sections le long de l'arête par différences : un λ arrêté par la
   tolérance de Brent les bruiterait. Pendant la sécante, Newton part du λ de
   la forme précédente, sans Brent.
5. **m** par la sécante, pour que le milieu B(1/2) soit aussi loin de S_b que
   celui de l'arc (r / cos(ψ/2) − r, où ψ est l'angle que l'arc balaie). Entre
   plans, la section est homothétique en m : la sécante est exacte au premier
   pas. Une sécante qui n'atteint pas la flèche à 10⁻¹⁰ près est refusée.

À angle droit : λ = 0,3525, contacts à 1,3053 r, pic de courbure de 1,234/r au
milieu, comme le prévoyait le § 4.4 du palier. Le test le recalcule sans gbs.

La station garde les centres de la boule (la ligne des centres, le paramètre)
et prend les contacts et les pôles de la section. Un arc réduit à un point
(faces tangentes) donne une section réduite au segment des contacts.

### 4.2 La surface

- **Extrusion exacte** (degré (4, 1)) quand la ligne des centres est droite et
  que les contacts de la section courent en ligne droite (u, v) entre ceux des
  bouts. C'est le cas entre deux plans, et entre un plan et un cylindre d'axe
  parallèle à l'arête. Le point du milieu est vérifié à tol/10 près ; sinon,
  c'est la peau.
- **Ailleurs** : la peau de Hermite (`skinned`) des sections exactes, de degré 4
  en travers. Son contrôle de normale compare la surface aux sections voisines,
  comme pour la conique. Les stations G2 déjà calculées sont gardées (nœuds,
  contrôles, plis).
- **Plis sur les faces B-spline.** La section G2 prend la courbure des faces à
  ses contacts. Là où un contact croise un nœud d'une face sous C3, ou la
  couture d'une face qui se ferme sans être périodique, cette courbure casse, et
  la section avec elle : le congé y fait un pli (4·10⁻⁵ rad au pied d'aube). Ces
  stations sont cherchées sur les contacts de la section elle-même (64 stations
  puis bissection), car ceux de la boule croisent les nœuds ailleurs. La peau y
  est C0 : les morceaux entre plis sont interpolés à part, avec les dérivées de
  chaque côté (différences d'ordre 2 décentrées), puis joints. Le contrôle de
  normale ne traverse pas un pli.
- **Boucle sur une face fermée non périodique** (le pied d'aube) : elle part de
  l'endroit où le contact G2 croise la couture, et non celui de la boule, pour
  que le contact reste dans les bornes de la face.
- Les contacts G2 sont plus larges que ceux de la boule : ils sont vérifiés à
  leur tour sur les faces (`check_on_face`).

### 4.3 Lois de rayon, chaînes, coins

- **Rayon en loi** (`fillet_law`) : la section G2 est construite sur chaque
  boule de l'enveloppe. Une loi anguleuse est refusée, avant le calcul des
  peaux : les fuseaux de la boule entre les morceaux ne raccordent que des arcs
  de cercle.
- **Chaînes** : les sections de part et d'autre d'un sommet de chaîne sont les
  mêmes (même plan, même rayon, mêmes faces), donc les congés se raccordent.
- **Coins** : `EdgeBlend::circular` est faux, et le coin devient un raccord de
  sommet du palier 6, G1 sur ses bords.
- **Corde constante** : refusée avec une section G2. Les contacts de la section
  ne sont pas ceux de la boule, la corde ne serait donc pas c.

### 4.4 La section de chaque arête

Comme la conique, la section G2 passe par les `SectionOf` de chaque arête :
`sections_on_edges` (`blend.h`, une ligne) et `section_of_edge` (pour
`fillet_geometry`). La première version la lisait dans les options
(`with_g2`), pour ne pas toucher à `blend.h` pendant la PR 244 ; celle-ci est
fusionnée, et la PR se fusionne sans conflit avec master. Avec les poches de la
PR 244, le congé du pont prend la section G2 comme les autres.

## 5. Refus

- `NonGeneric` (sur l'arête) :
  - un congé G2 dont la loi de rayon est anguleuse ;
  - à un sommet de chaîne où la loi de rayon change de pente, des congés G2 qui
    ne se raccordent pas (« blends that do not join there » : le plan de la
    section penche avec la loi, et les contacts G2 de part et d'autre
    diffèrent) ;
  - un congé G2 à corde constante ;
  - un arc de plus de 170° ;
  - des contacts introuvables sur les faces ;
  - des tangentes qui ne se coupent pas ;
  - une sécante sur m qui n'atteint pas la flèche de l'arc ;
  - sur une face fermée non périodique, un contact G2 qui ne croise pas la
    couture près de celui de la boule.

## 6. Tests

`tests_brep_g2_sections` (12 cas). Les références sont calculées sans gbs :
la règle C est refaite en 2D dans le test (section d'or puis Newton sur
l'énergie), et l'aire de la quartique est obtenue par la formule de Green.

| Cas | Attendu | Mesuré |
|---|---|---|
| règle C à angle droit | λ 0,3525, m 1,3053, pic 1,234/r, courbure nulle au contact | 10⁻⁴ (valeurs de la note), 10⁻¹² |
| arête de boîte, r = 0,3 | extrusion exacte de degré 4 ; contacts à m r ; courbure en travers 0 aux contacts, celle de la quartique au milieu ; flèche de l'arc ; volume 3 − 3 A | 10⁻⁹ |
| arêtes aiguë (63,4°) et obtuse (116,6°) | contacts de la règle C à ces angles ; volumes | 10⁻⁹ |
| plan et cylindre le long de l'arête | extrusion exacte ; courbure 1/R au contact sur le cylindre, 0 sur le plan | 10⁻⁷ |
| bord d'un cylindre (peau) | contacts à m r ; Pappus | 10⁻⁶ ; 10⁻⁵ |
| coin de trois congés G2 ; dalle en stade (chaîne) | raccord de sommet sans sphère ; solides valides | 10 faces chacun |
| loi de rayon lisse de 0,2 à 0,35 | solide valide, volume entre ceux des rayons extrêmes | — |
| loi de rayon sur l'arête plan-cylindre (peau) | courbure en travers 1 sur le cylindre, 0 sur le plan | 10⁻³ |
| refus | corde constante ; loi de rayon anguleuse ; chaîne de la dalle dont la loi change de pente à un sommet | messages du § 5, modèle inchangé |
| arête concave (intérieur d'un L), r = 0,3 | volume 6 + 2 A | 10⁻⁸ (mesuré 5·10⁻¹⁰) |
| pied d'aube (B-spline, nœuds C0, couture non périodique), r = 0,03 | solide valide ; matière ajoutée / celle du congé circulaire = 0,232084 / (1 − π/4) | 0,5 % (mesuré 0,07 %) |
| loi de 0,08 à 0,12 sur le bord d'un prisme à profil B-spline (nœuds C2) | matière ôtée / celle du congé circulaire de même loi, même rapport | 0,5 % (mesuré 0,02 %) |

En Python : la section G2 sur l'arête de boîte (contacts, aire 0,232084 r²,
refus de la corde). La galerie gagne `fillet_g2_box` (un congé circulaire et un
congé G2 côte à côte) et `fillet_g2_corner`.

## 7. Comparaison à OCCT

Aucune : OCCT n'a pas de congé G2. Les congés circulaires et coniques ne
changent pas, et les tests des paliers 4 à 7 le vérifient toujours.

## 8. Limites

- **Raccords de sommet G1** entre congés G2 (non-objectif du palier).
- **Plis** : sur une face B-spline dont la courbure casse (nœuds sous C3,
  couture), le congé G2 suit ce défaut de la face et n'est que C0 en travers de
  ces stations.
- **Faces courbes** : le contact est défini par sa projection sur la tangente
  de l'arc de la boule. C'est exact entre plans et naturel ailleurs, mais ce
  n'est pas une norme : CATIA et SolidWorks ne publient pas la leur.
- **Coût** : une station minimise l'énergie (Brent, puis Newton) au premier pas
  de la sécante, puis Newton seul aux suivants ; les stations sont gardées. Le
  pied d'aube prend 3 à 7 s, contre 0,1 s pour le congé circulaire ; les 12 cas
  des tests, 6 s en tout.
- **Courbure le long de l'arête** : avec un rayon en loi, rien de plus que G1
  n'est garanti le long de l'arête (§ 1.2 du palier).

## 9. Points à valider

1. **Règle C par station** : contacts par projection sur la tangente de l'arc,
   courbure des faces par Meusnier, λ de la forme la plus lisse (Brent puis
   Newton sur la dérivée à 5 points), m par la sécante pour la flèche de l'arc.
2. **Extrusion exacte de degré (4, 1)** quand la ligne des centres et les
   contacts sont droits (plans, plan et cylindre parallèle) ; peau ailleurs.
3. **Rayon en loi** avec la section G2 sur chaque boule ; loi anguleuse refusée.
4. **Corde constante refusée** avec une section G2.
5. **G2 dans les sections de chaque arête** (§ 4.4) ; plis du congé G2 sur les
   faces B-spline dont la courbure casse (§ 4.2).
6. **Bilan du palier** (§ 9 du document du palier).

## 10. Revue de code

Revue du 23 novembre 2026 : 6 points, tous traités.

1. **Faces B-spline** (important). Au pied d'aube, le congé G2 était refusé
   après 167 s (« cannot be approximated within approx_tol »). La courbure de
   l'aube casse à ses nœuds C0 et à sa couture, et la section G2 avec elle ;
   les stations n'étaient placées qu'où les contacts de la boule croisent les
   nœuds. Désormais, les plis sont cherchés sur les contacts G2, la peau y est
   C0, et la boucle part de la couture du contact G2 (§ 4.2). Le pied d'aube
   est construit (test) ; sur un prisme à profil B-spline (nœuds C2), le congé
   passe de 6,6 s à 1 s.
2. **Coût**. Les stations G2 sont gardées, l'angle des arcs n'est plus calculé
   pour G2, et Newton part du λ précédent pendant la sécante. Le bord d'un dôme
   passe de 3,7 s à 1,4 s. Une loi anguleuse est refusée avant le calcul des
   peaux (0,01 s au lieu de 23 s au pied d'aube).
3. **Sécante et messages**. La flèche est vérifiée après la sécante (refus
   sinon), avec un message propre ; les refus d'une station G2 ou d'une peau
   portent l'arête (`with_edge`).
4. **Tests** : arête concave, pied d'aube, loi sur un prisme B-spline, refus
   d'une chaîne G2 ; 12 cas (et non 10 comme annoncé).
5. **Report après la PR 244**. La section G2 passe par les sections de chaque
   arête (§ 4.4). Le message « conic fillets of a chain » n'est en fait jamais
   donné pour G2 : à un sommet de chaîne où la loi change de pente, les contacts
   G2 de part et d'autre diffèrent, et le refus est « blends that do not join
   there » (§ 5, test). Le § 8 qui l'annonçait est corrigé.
6. **Lien dans la description de la PR** : retiré.
