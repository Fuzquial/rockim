# Calibration Red Bohus — protocole (ouvert le 2026-08-15)

Calibration des paramètres de joints FDEM de rockim sur le granite **Red
Bohus**, selon la méthodologie de l'état de l'art 2025-2026 (§4). Objectif :
un jeu de paramètres **avec son incertitude**, capable de reproduire UCS,
brésilien et enveloppe triaxiale — et de **prédire** les confinements non vus.

## 1. Données expérimentales (cibles)

Source : **Dumoulin et al. (2024)**, *Geomechanics for Energy and the
Environment* **40**, 100592 ; dataset Zenodo `10.5281/zenodo.10617548`.
Granite du sud-ouest de la Suède, **60 % feldspath / 35 % quartz / 5 % biotite**
(poids), taille de grain 1-3 mm.

Extraction : `tools/extract_targets.py` → `targets/targets_redbohus.json`
(scalaires + courbes σ-ε moyennes rééchantillonnées + écarts-types).

| Cible | Valeur | Écart-type | Rôle |
|---|---|---|---|
| UCS (4 essais) | **126,6 MPa** | ±21,4 (17 %) | calibration |
| BTS (4 essais, recalculé des .ASC) | **10,27 MPa** | ±0,98 | calibration |
| q(σ₃ = 20) (3 essais) | **404,8 MPa** | ±2,8 | calibration |
| q(σ₃ = 50) (3 essais) | **599,0 MPa** | ±2,6 | calibration |
| q(σ₃ = 75) / q(σ₃ = 100) | 704,0 / 799,3 MPa | ±7,6 / ±4,2 | **prédiction pure** |
| E / ν | **77,7 GPa / 0,29** | — | sortie à reproduire |

⚠️ **Pièges documentés.** (a) Trois élasticités circulent dans l'archive pour
ce granite (52/0,25 littérature DP-DFH ; 57,3/0,17 « local moyen » PSO ;
77,66/0,29 fit des 12 branches triaxiales) — on retient la dernière et E, ν
deviennent des **sorties** du modèle, pas des entrées figées (méthode Bu 2026).
(b) La cible de traction est le **BTS mesuré** (10,3 MPa), pas le σt = 18,3 MPa
Weibull de Saadati/Shariati qui avait servi au GBM de juillet. (c) L'enveloppe
est **concave** (pente locale 13,9 → 3,8) : Mohr-Coulomb linéaire ne peut pas
la suivre partout, d'où le partage calibration / prédiction. (d) La dispersion
expérimentale de l'UCS est de 17 % — viser mieux n'a pas de sens physique ;
les triaxiaux (±0,5 %) sont les cibles exigeantes.

## 2. Méthode (état de l'art)

- **Ye et al. 2025** (IJRMMS 194, 106233) : objectifs sur la **courbe entière**
  σ-ε (pas le seul pic) + le **mode de rupture** comme objectif formel,
  optimisation **multi-objectif NSGA-II** (front de Pareto, pas de poids
  arbitraires) ; bulk élasto-plastique Mohr-Coulomb.
- **Bu et al. 2026** (IJRMMS 199, 106400, open access) : base de 3456
  simulations UDEC-BBM (UCS + confinés 10/20/30 + brésilien), comparaison
  RF / SVR / **GPR** / DNN, inversion par grid search, validation par 1485
  runs ; cibles = E, ν, UCS, BTS, c, φ.
- **Jiang et al. 2025** (Sci. Rep. 15:34923) : **328 runs suffisent**,
  analyse de corrélation pour réduire la dimension, stacking d'ensembles,
  erreurs finales 0,6 % (UCS) à 10,6 % (BTS) ; taille de bloc = L/20 ;
  **le BTS est la sortie la plus dispersée (COV > 8 %) → plusieurs graines**.

Tous calibrent des modèles à **blocs Voronoï** (BBM/GBM), pas des maillages
homogènes : c'est la structure qui produit l'enveloppe non linéaire et le bon
ratio UCS/BTS.

**Notre apport par rapport à ces trois papiers** : un **postérieur bayésien**
sur l'émulateur (incertitude + identifiabilité + corrélations entre
paramètres), et un **test de prédiction** sur des confinements non vus.

## 3. Phases

| Phase | Contenu | Budget |
|---|---|---|
| **0** | cibles (fait) ; géométries réelles ; grain L/20 ; indépendance au taux ; 3 graines pour le BTS | ~1 h |
| **1** | **criblage** mono-variable + corrélations → quels paramètres pilotent quoi | ~50 runs |
| **2** | base de données : hypercube latin sur les paramètres retenus × 4 essais | ~250 jeux, ~2 h parallélisé |
| **3** | émulateur (GPR vs RF, protocole Bu) puis **postérieur bayésien** + front de Pareto en contrôle | minutes |
| **4** | validation : 3 graines au jeu calibré + **prédiction σ₃ = 75 et 100** | ~12 runs |

## 4. Paramètres candidats

`ft_joint`, `cohesion_joint`, `frictionDeg_joint`, `Gf` (mode I),
`gfShearFactor` (mode II), `jointPenaltyFactor`, `crushCap` (cap déviatorique
du bulk — proxy de la plasticité MC de Ye), `E_bloc`.

## 5. Discipline

Schémas figés pour toute la campagne (capacités validées du solveur) :
`insertion = adaptive`, `jointSoftening = yan`, `jointShearUnload = origin`,
`contact = potential`, `gcActivation = adaptive`, `stopPeakDrop`,
`budgetAbortPct`/`budgetAbortMin`. Toute configuration de run est soumise à
validation de Fernando AVANT lancement.
