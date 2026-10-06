# Stratégie de lancement (2026-10-06, 22 h 30)

Règle : la session cloud ne lance que des calculs courts (moins de 30 min chacun à 4 fils).
Tout ce qui dépasse va sur le poste du doctorant (14 fils), lancé par lui.
Un seul calcul à la fois. Après chaque groupe : `RESULTATS.md` du groupe, commit, push.

## Cloud, dans cet ordre

| ordre | groupe | calculs | durée estimée (4 fils) | pour le chapitre |
|---|---|---|---|---|
| 1 | `yan` | 9 (en cours) | 1 à 2 h au total | r05a : rejouer août avec le binaire courant |
| 2 | suite rapide | `python3 tools/verify_suite.py --tier fast`, SANS limite de temps | 15 à 20 min | repère 50/51 (le rejeu du 06/10 a donné 49/51 : maillage T1 introuvable, à noter) |
| 3 | `calib` | 4 | 20 à 40 min | r07a : jeu de juillet avec le binaire courant |
| 4 | `diag_contact` | 3 | 25 min | r05c : origine de l'injection d'énergie par le contact (voir ci-dessous) |
| 5 | `abuaisha` | `hf_*_c`, `hf_*_m6` seulement s'ils ne sont pas déjà dans un RESULTATS | 10 à 40 min chacun | r05b : rejeu grossier de la pression de rupture |
| 6 | `impact` | `bench1_elas`, `bench1_frac` à T = 1,2e-4 s | 30 à 60 min chacun | r06 : décharge, durée de contact contre Hertz (89,4 µs), restitution |

`tunnel` (red_*) : déjà rejoué le 06/10, ne pas relancer.

### Diagnostic de l'injection d'énergie par le contact (groupe `diag_contact`)

Base : `tunnel/red_tip16.cfg` (8 min 42 s à 1 fil le 06/10 ; travail net du contact
1,94e6 J/m contre 9,57e5 J/m d'énergie cohésive). Une seule clé change par calcul :

| config | clé changée | hypothèse testée |
|---|---|---|
| `tip16_plastic` | `jointShearUnload = plastic` | la décharge `origin`, non conservative, crée l'énergie |
| `tip16_gcpen` | `gcBirth = penalty` | la naissance du contact sur un joint mort (rampe) crée l'énergie |
| `tip16_mu0` | `contactMu = 0` | le frottement de contact (ressort plus Coulomb) crée l'énergie |
| `tip16_ratchet` | `jointSecantRatchet = true` | la décharge `origin` SANS cliquet crée l'énergie (piste du chapitre coupe : +115 J/m dans les joints, supprimé par le cliquet) |

À extraire du `.log` de fin de run (bloc de bilan) : travail net du contact, énergie cohésive,
dissipation de Cundall, résidu B4, joints rompus ; comparer à la base ci-dessus.

## Poste local (14 fils), à lancer par le doctorant

| calcul | durée | pourquoi |
|---|---|---|
| tunnel de production (maillage 106 k triangles), binaire courant, `red_adapt` avec le maillage `tunnel_hs_iso.msh` | ~45 min | lever la réserve « maillage réduit trop grossier » ; refaire le bilan contact |
| St Anne prolongé à 500 µs (`configs/stanne2025_rock137_visc0_T500.cfg`) | plusieurs heures | rebond |
| St Anne sur `rock073` (251 460 tétraèdres) | plusieurs heures | convergence en maillage |
| Kuru maillage fin jusqu'au retournement (> 255 µs) | plusieurs heures | phase de décharge |
| AbuAisha de production (maille 3 mm) | 2 à 2,5 h chacun | rejouer les 7 calculs avec le binaire courant |
| Parker de production | ~11 h | lever la réserve de maillage grossier |
| Yan, compression à maille 0,75 mm | ~1 h | maillage de l'article |
| bancs `vv/` P2.2, P3.1, P4 | 2 h à 58 h | campagne V&V |
