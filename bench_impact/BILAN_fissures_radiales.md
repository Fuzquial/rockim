# Pourquoi nous n'avons pas (encore) l'étoile de fissures radiales

Bilan de la campagne d'impact à insert unique du 2026-08-22 (spec 005),
contre Yang, Xiang, Naderi, Wang, Aising, Ugarte, Latham — IJRMMS 191 (2025)
106125. Trois runs, arrêtés en connaissance de cause ; état au soir.

---

## 1. Ce qui était identique à l'article — et vérifié, pas supposé

Le run 3 (« fidèle », roche au 1 mm) partageait avec leur simulation :

* **la roche au chiffre près** (leur Table 4) : ρ 2731, E 57 GPa, ν 0,31,
  f_t 7 MPa, c 18,8 MPa, φ 45°, G_I 12, G_II 800, glissement 0,6 ;
* **les lois cohésives équation par équation** : amorçage Mohr-Coulomb à
  cut-off (leur éq. 1, `jointShearEnvelope = yang`), adoucissement f(D) =
  la z-curve de Munjiza (a 0,63 / b 1,8 / c 6, ∫f dD = 0,386307 — Yan
  l'avait reprise verbatim, vérifié dans `YanSoftening.hpp`), couplage
  mode I-II elliptique (leur éq. 3 = le moteur √(rn²+rs²) avec `origin`),
  DIF avec leur allocation exacte et l'exposant 0,17 de leur papier 2 ;
* **la géométrie et le maillage de roche** : cylindre Φ250 × 150, gradients
  1/2/10 mm de leur fig. 6, assemblage complet à 6 corps (piston, bit,
  insert brasé, circlip, plaque posée), insert posé sur la roche, gravité ;
* **le volume quasi élastique** comme le leur (`crushCap` relevé à 2 GPa).

La seule différence conservée — délibérément, c'est le sujet de la thèse —
est le **schéma d'insertion** : adaptatif (extrinsèque, Yan 2023) contre
intrinsèque (leurs 230 000 éléments cohésifs présents dès t = 0).

## 2. Ce qui s'est passé

| critère | run 1 (cap 150 MPa) | run 3 (cap 2 GPa) | leur Table 3 |
|---|---|---|---|
| v. d'indentation | 8,23 m/s | 7,81 | 9,40-9,85 |
| v. de rebond | (perdu) | **7,71** | 6,87-7,10 |
| rapport rebond/indentation | — | **0,99** | 0,72-0,73 |
| enfoncement | 1,85 mm | 0,90 | ~1,53 |
| fissures | halo diffus 26,7 mm | disque compact 9 mm | radiales 20-24,5 mm |
| joints rompus | 3 700 | 1 754 | — |

Aucun des deux runs ne produit l'étoile radiale de leurs fig. 11/14/16.
Le run 1 fissure large mais flou (et trop mou : le cap à 150 MPa écrêtait
la force de pointe — voir la courbe F-p, pic mesuré ~40 kN contre leur
rampe à 130). Le run 3 corrige la force (116 kN à 0,95 mm, la bonne rampe)
mais rebondit **quasi élastiquement** : 99 % de la vitesse rendue, presque
pas d'hystérésis F-p, fissuration arrêtée à 9 mm de rayon.

## 3. Le diagnostic

**À lois cohésives identiques, la différence d'insertion est une différence
de comportement effectif de la zone broyée.** Chez eux, la dissipation sous
l'insert est portée par le micro-broyage de milliers de joints intrinsèques
disponibles partout, immédiatement. Chez nous, l'adaptatif n'offre de
surface de rupture que là où le critère de face a tiré : 1 754 joints ont
existé dans tout le run 3. Moins de surface créée, moins d'énergie
consommée (17,7 J dissipés sur 52 fournis à mi-run, et la boucle se referme
ensuite), donc un rebond élastique, un enfoncement court — et pas de
poussée durable pour nourrir les radiales. La roche fine au 1 mm était
prête à les résoudre ; c'est la **loi de volume** qui ne les a pas
alimentées.

Les deux runs encadrent le vrai besoin : cap à 150 MPa = zone broyée trop
molle (enfoncement +20 %) ; volume élastique = trop dur (−40 %, restitution
totale). Ce qui manque entre les deux est une **dissipation progressive de
la zone broyée** — précisément ce que le modèle de pulvérisation de leur
papier 2 (IJRMMS 206, 2026) fournit, et que rockim possède depuis cette
nuit (`bulkDamage = yang`, WP1 : plancher exact, dissipation Y dD, 2D+3D).
Leur papier ne l'applique qu'au granite ; notre résultat suggère que le
schéma adaptatif en a besoin aussi pour le calcaire, à dose plus faible,
comme substitut du micro-broyage intrinsèque.

## 4. La suite (décidée le 2026-08-22 au soir)

Runs **rapides** de pulvérisation, montage allégé : plus de piston ni de
plaque — le corps bit+insert est lancé directement à la vitesse
d'indentation mesurée (9,5 m/s), la roche à raffinement moyen. Objectif
unique : produire l'étoile radiale et le cratère par `bulkDamage`, en
balayant sa calibration (δ₀, δ_f). Décks `impact_pulv_*.cfg`, en file
derrière le balayage λ du tunnel.

Réserve de méthode : le run 3 a été arrêté à 41 % (t = 413 µs) pour rendre
la machine — pic, indentation, rebond et faciès de charge sont acquis ;
le tri des fragments et la phase tardive ne le sont pas. Ses figures
partielles et son historique sont archivés (`out_imp_fidele`, drive).
