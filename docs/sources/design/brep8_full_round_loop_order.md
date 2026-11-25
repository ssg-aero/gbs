# Congé tritangent en boucle : l'ordre des côtés

| | |
|---|---|
| PR | branche `fix/brep-full-round-loop-order` (sur master d9b8731) |
| Document de référence | [brep8_pr04_full_round_loops.md](brep8_pr04_full_round_loops.md) (PR 4, #256) ; [brep_palier8_face_blends.md](brep_palier8_face_blends.md), § 3.5, § 9 |
| Fichiers | `gbs-brep/full_round.h`, `tests/tests_brep_full_round.cpp`, `docs/sources/design/brep_palier8_face_blends.md` (§ 9), `News.md`, cette note |
| Statut | À valider (§ 5) |

## 1. Le défaut

`fillet_full_round` sur le bord d'un disque dépendait de l'ordre des côtés :
`[bottom], [rim], [top]` construisait le demi-tore exact (volume à 5,4·10⁻¹⁰),
`[top], [rim], [bottom]` refusait « the result is not a valid solid: outer wire
not direct ». `full_round_geometry` rendait le même tore dans les deux ordres :
seule la topologie échouait. Trouvé en écrivant le guide, dans master après
la PR 4 ; ni ses tests ni sa revue de code n'avaient essayé les deux ordres sur
une boucle (seulement sur les tubes carrés et les nervures).

Rejoués dans les deux ordres sur master, les cas de boucle de la PR 4 :

| Cas | Dans l'autre ordre, sur master |
|---|---|
| bord d'un disque (le dessus d'abord) | refus « outer wire not direct » |
| bord d'une ailette | refus « outer wire not direct » |
| gorge d'un arbre (les parois échangées) | construit, mais une peau au lieu du tore exact (volume à 6,1·10⁻⁹) |
| extrémité d'un tube | même solide (4,3·10⁻¹⁰) |
| secteurs de disque, 120° et 300° | même solide |
| départ à une couture en hélice (tours 0, 0,5, 1,5) | même solide (10⁻⁹) |

## 2. Les causes

**L'orientation du guide dans F₂.** Sur la copie, M est réduite à son guide,
qui prend dans F₂ la place de l'autre arête de M (`collapse`). Son usage dans
F₂ était celui du guide par M. Un usage est relatif à l'orientation de sa face
dans la coque : il ne vaut pour F₂ que si F₂ et M y ont la même orientation.
Le fond d'un disque est retourné dans la coque, pas son bord ; donné en second
côté, il recevait le guide à l'envers, son fil tournait dans le mauvais sens
en (u, v) sur la copie, et le contact fermé qui le remplace (`blend.h`, la
boucle tournée comme le fil qu'elle remplace) aussi. Dans l'ordre qui
marchait, F₂ était le dessus, orienté comme le bord. Un guide ouvert n'y était
pas sensible : ses contacts sont placés par les sommets de la chaîne, pas par
le sens du fil (les secteurs passaient dans les deux ordres).

**Les nœuds cachés du guide** (la gorge). Les parois échangées, la boucle part
de la couture d'une paroi (une révolution B-spline fermée) : le guide devient
une `LoopCurve` qui lit l'arête au-delà de son sommet et n'a pas de nœuds à
elle. Les pcourbes des contacts, interpolées par morceaux entre les nœuds du
guide (PR 4, § 3.3), l'étaient alors d'un seul tenant, et le contrôle du tore
échouait (la vitesse d'un cercle rationnel casse aux nœuds de la révolution).
De plus, le guide était ici sur le fond de la gorge (un cylindre, sans nœuds) :
ce sont les contacts sur les parois qui croisent leurs nœuds, aux stations
mêmes (u = π/2 en t = π/2), que la recherche d'un changement de signe strict ne
voyait pas.

## 3. La correction

- **L'usage du guide par F₂** est celui de F₂ pour l'autre arête, inversé si le
  guide la parcourt en sens contraire (`q.same`), et non celui de M. Ce n'est
  pas un échange des côtés en entrée : l'échange aurait masqué le défaut pour
  le disque (le bas d'abord) sans le corriger pour un solide dont les deux
  côtés sont retournés par rapport à M.
- **Les nœuds du guide** sont lus sur la courbe de son arête (et, pour une arête
  sur une surface, là où sa pcourbe croise les nœuds de la surface), ramenés
  dans [a, b] par périodes sur une boucle ; pour les morceaux des surfaces
  exactes, les plis et les cassures de la peau.
- **Les morceaux des surfaces exactes** sont aussi coupés là où un contact
  croise un nœud de sa face sous C³, un contact posé sur le nœud à une station
  compris.

Avec, dans la même PR (même fichier) :

- **Le repli cherché autour de la station qui échoue** (mineur laissé par la
  revue de code de la PR 4). Quand la marche s'arrête, ou qu'une station
  demandée après elle échoue, le repli d'une face est cherché sur les arêtes
  des trois faces à moins de 4 rayons (plus le pas du guide) de la station qui
  échoue, au rayon de la boule résolue la plus proche, et non plus sur toutes
  les arêtes au rayon de la première boule : sur un congé à rayon variable, une
  face très courbe loin de là pouvait donner le message de repli au lieu de
  « no ball ». `fold_near` (`full_round_detail`) ; le message donne aussi le
  rayon de courbure trouvé.
- **Le bilan du § 9 du palier** : la face du milieu qui fait le tour et le tore
  exact d'une nervure courbe, faits par la PR 4, sortent de ce qui reste à
  faire ; une ligne du tableau pour la PR 4, ses écarts.

## 4. Tests

`tests_brep_full_round` : 42 cas (39 avant), 9 759 assertions.

| Cas | Attendu | Mesuré |
|---|---|---|
| bord d'un disque, les deux ordres (écrit d'abord : échoue sur master) | même tore exact, Pappus | 5,4·10⁻¹⁰ dans les deux |
| gorge, tube, secteurs, ailette, départ à la couture en hélice, côtés échangés (`expect_same_swapped`) | valide, mêmes faces, même volume | gorge 1,6·10⁻¹¹ (master : 6,1·10⁻⁹), tube 4,3·10⁻¹⁰, secteurs 6·10⁻¹⁶, ailette 2·10⁻¹⁸, couture 10⁻⁹ |
| `fold_near_the_station` : tube carré aux coins extérieurs de 0,05 | au milieu d'un côté (coins à 0,95) : rien ; près d'un coin, r = 0,15 : le côté extérieur, 0,05 ; r = 0,04 : rien | vérifié |
| refus de repli existants (ailette de 0,01, tube aux coins de 0,05 et 0,14) | inchangés | inchangés |
| disque cousu face par face (`sewn_disc_orientations`) : dessous et dessus renversés ou non, cousu dans les deux ordres, côtés dans les deux ordres (16 cas) | tore exact, même volume | 5,4·10⁻¹⁰ dans les 16 |
| cavité en disque dans un bloc (`disc_cavity`) : coque intérieure, boucle concave, les deux ordres | tore exact, matière ajoutée par Pappus | 5,4·10⁻¹⁰ dans les deux |

Pourquoi cinq passes (tests de la PR 4, revue de code, trois relectures) ont
manqué le défaut : toutes les boucles venaient de `make_prism` ou de
`make_revolution`, toujours dans le même ordre des côtés, sans jamais varier
l'orientation des faces dans la coque. La relecture de 644850e l'a vérifié
sur 64 guides ouverts (orientations de F₁, M, F₂ dans la coque, ordre de
couture, miroir, ordre des côtés : 64 construits), 16 disques cousus (8 refusés
sur master) et la cavité en disque ; les deux derniers sont devenus les tests
ci-dessus, qui font varier ces orientations. Le volume de référence y est celui
du solide d'entrée donné exactement (son cylindre s'intègre à 5·10⁻⁸ près).

Le cas d'un congé à rayon variable où l'ancienne recherche nommait à tort une
face lointaine n'a pas été construit (vu à la lecture par la revue, non
sondé) ; `fold_near_the_station` montre la recherche locale.

## 5. Point à valider

1. **L'usage du guide dans F₂ pris sur F₂ elle-même**, plutôt qu'un échange des
   côtés en entrée. *Recommandé* : c'est la cause, et l'échange ne couvrirait
   pas un solide aux deux côtés retournés.
2. **Le repli cherché autour de la station qui échoue**, à 4 rayons plus le pas,
   au rayon local. *Recommandé* : le message nomme la cause là où la marche
   s'arrête.
