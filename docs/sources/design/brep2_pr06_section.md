# Palier 2 — PR 6 : section

| | |
|---|---|
| PR | branche `feat/brep-section` |
| Document de référence | [brep_palier2.md](brep_palier2.md), § 4 |
| Fichiers | `gbs-brep/section.h` (nouveau), `gbs-brep/brep` (agrégateur), `tests/tests_brep_section.cpp` (nouveau), `docs/index.rst`, `News.md` |
| Statut | À valider |

## 1. Objectif de la PR

Calculer les arêtes où les faces de deux formes se rencontrent, ou celles où
une surface coupe une forme : plan de coupe, contour d'usinage. Chaque arête
porte sa pcurve sur les deux faces. C'est la première étape de la découpe
(PR 7) et des booléens (PR 9).

## 2. API

```cpp
struct SectionOptions { T tol = brep_default_tolerance; IntersectOptions<T> intersect; };
struct SectionEdge    { EdgeId edge; FaceId face_a, face_b; shared_ptr<Curve<T, 2>> pcurve_a, pcurve_b; bool tangent; };
struct SectionResult  {
    CompoundId compound;                    // wires et arêtes de la section
    vector<WireId> wires;
    vector<SectionEdge> edges;              // arêtes nouvelles, intérieures aux deux faces
    vector<EdgeId> boundary_edges;          // arêtes des formes situées dans la section
    vector<VertexId> points;                // contacts ponctuels
    vector<pair<FaceId, FaceId>> coincident;
    bool singular, complete;
};

auto section(m, ShapeId a, ShapeId b, SectionOptions = {}) -> BuildResult<SectionResult>;
auto section(m, ShapeId a, shared_ptr<Surface<T, 3>> tool, optional<array<T, 4>> box = {}, SectionOptions = {}) -> BuildResult<SectionResult>;
```

## 3. Choix d'architecture et justification

### 3.1 Étapes

1. **Paires de faces** dont les boîtes englobantes, élargies de la
   tolérance, se recouvrent. Pour une surface de coupe, toutes les faces.
2. **Intersection** des surfaces sur la **boîte des pcurves** de chaque face
   (PR 3 et 4), pas sur les bornes des surfaces.
3. **Restriction aux faces.**
   - Chaque courbe est coupée là où ses pcurves croisent les pcurves du bord
     des deux faces (courbe-courbe, PR 1).
   - Le milieu de chaque morceau est classé dans les deux faces (PR 5) : s'il
     est **dedans** pour les deux, le morceau devient une arête.
   - Un **croisement de bord** l'emporte sur un bout de courbe à moins d'une
     tolérance : le bout peut tomber dans la marge de la boîte paramétrique.
   - Les deux bouts d'une courbe fermée sont un même point, mais pas une même
     coupure.
4. **Bord d'une face.** Un morceau dont le milieu est *sur* le bord d'une face
   suit une arête existante : c'est elle qui est la section à cet endroit. Elle
   est rapportée dans `boundary_edges`, une seule fois, sans être dupliquée.
   C'est le cas d'un plan qui passe par deux arêtes d'une boîte, ou qui
   contient l'axe d'un cylindre (le seam).
5. **Sommets.**
   - Les bouts sont confondus avec les sommets des formes, puis avec les
     sommets déjà créés, dans un rayon de **10 tol**. C'est l'écart que peut
     avoir le bout d'une courbe de cheminement interpolée.
   - La tolérance de l'arête monte jusqu'à couvrir son écart et la distance à
     ses sommets.
   - Une arête de section qui s'arrête sur le bord d'une face partage donc le
     sommet de ses voisines, et les formes ne sont pas modifiées.
6. **Wires.**
   - Les arêtes nouvelles et les arêtes de bord sont groupées par composantes
     connexes de leurs sommets, puis chaînées par `make_wire`.
   - Une composante qui se ramifie (plus de deux arêtes à un sommet) reste en
     arêtes libres dans le compound.

### 3.2 Ce qui n'est pas fait ici

- **Faces confondues** : elles sont seulement rapportées. Leur découpe
  commune est faite par la découpe et les booléens (§ 6.3 du document de
  conception).
- **Arêtes des formes** : elles ne sont pas coupées aux points de section.
  C'est le rôle de la découpe (PR 7).

## 4. Écarts par rapport au document de conception

| Document | Implémentation | Raison |
|---|---|---|
| `SectionResult` : compound et faces de chaque arête | plus `boundary_edges`, `points`, `coincident`, `singular`, `complete` | sections le long des arêtes existantes, contacts, faces confondues, diagnostics |
| fusion des sommets « à la tolérance d'opération » | rayon de 10 tol, tolérances des arêtes relevées | bouts des courbes interpolées |

## 5. Tests (`tests_brep_section`, 3 cas, 556 assertions, 0,1 s)

Chaque arête est vérifiée en 11 points sur ses deux surfaces, par ses
pcurves, à sa tolérance près, et ses sommets sont vérifiés à ses bouts.

| Test | Couvre |
|---|---|
| `box_by_planes` | boîte 2 × 3 × 4 coupée à z = 1 (4 arêtes, longueur 10, un wire fermé, boîte inchangée) ; par un plan passant par deux arêtes verticales (2 diagonales sur les sommets de la boîte, 2 arêtes de bord, rectangle fermé) ; par le plan de la face du haut (face confondue, bord de la face en 4 arêtes de bord) ; plan éloigné (rien) |
| `cylinder_and_torus` | cylindre : cercle (une arête fermée, 10π) ; ellipse oblique coupée au seam (périmètre de Ramanujan, un wire fermé) ; plan contenant l'axe (génératrice, deux diamètres, seam en arête de bord : rectangle fermé) ; tore partiel d'OCCT par un plan oblique (cheminement, wires fermés) |
| `two_solids` | boîtes 0..2 et 1..3 : hexagone gauche de longueur 6, fermé ; sphère tournée centrée sur le flanc du cylindre : une boucle fermée à travers le seam (cheminement) ; erreurs |

Les suites de classification, de sewing, de contrôle et d'écriture STEP
passent avec l'agrégateur modifié.

## 6. Limite constatée

La sphère de `sphere.stp` posée sur le flanc du cylindre **sans rotation** a
ses deux pôles sur le cylindre. La boucle de section passe par ces pôles, où le
paramétrage de la sphère dégénère, et elle est livrée en deux morceaux non
raccordés. Le test tourne la sphère. Le cheminement à travers un pôle est un
point à reprendre (PR 7 ou 9), dans l'espace des paramètres de la face qui
contient le pôle.

## 7. Points à valider

1. **Restriction par classement du milieu** des morceaux dans les deux faces.
2. **Arêtes existantes rapportées** dans `boundary_edges` plutôt que
   dupliquées.
3. **Rayon de confusion des sommets** à 10 tol, et tolérances des arêtes
   relevées.
4. **Composantes ramifiées** laissées en arêtes libres.
5. **Limite des pôles** (§ 6), à traiter plus tard.
