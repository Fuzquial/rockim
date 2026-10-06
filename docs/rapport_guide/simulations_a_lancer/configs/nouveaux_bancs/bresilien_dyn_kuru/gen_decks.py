#!/usr/bin/env python3
# Écrit les quatre decks SHPB du banc B11 (brésilien dynamique, granite de Kuru,
# Saksala et al. 2013, IJRMMS 59, 128-138, doi:10.1016/j.ijrmms.2012.12.011).
# Lancer depuis n'importe où : python3 .../bresilien_dyn_kuru/gen_decks.py
#
# Impulsions incidentes : approximations linéaires par morceaux de la fig. 2a de
# l'article (essais 1, 6, 11, 16), numérisées à l'oeil sur le PDF à ±10 µs et
# ±10 MPa, puis ramenées au trapèze symétrique de rockim (shpbPulse = trapezoid :
# montée et descente linéaires de (1 - plateau)/2 * tau chacune).
# Vitesse particulaire V0 = sigma_inc / (rho_b c_b), acier rho_b = 7850 kg/m3,
# E_b = 200 GPa (valeurs SUPPOSÉES : Saksala dit seulement « high strength
# steel » ; ce sont celles de l'acier de Yang 2026, Table 1), c_b = 5047,5 m/s,
# rho_b c_b = 39,62 MPa.s/m.
# La comparaison à l'essai se fait à la vitesse de chargement MESURÉE dans le
# calcul (pente de sigma_t(t) entre 20 et 80 % du pic), contre l'ajustement de
# la fig. 5 de l'article : sigma_t = 0,017 x (GPa/s) + 21 MPa. Une erreur de
# numérisation de l'impulsion déplace donc le point le long de la droite, sans
# fausser le verdict.
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RHO_B, E_B = 7850.0, 200e9
Z_B = RHO_B * math.sqrt(E_B / RHO_B)          # impédance de l'acier [Pa.s/m]

# nom : (essai Saksala, vitesse du percuteur [m/s], sigma_inc [MPa],
#        montée [µs], palier [µs], descente [µs]) lus sur la fig. 2a
PULSES = {
    "shpb_v05": ("test 1", 5, 115.0, 120.0, 45.0, 115.0),
    "shpb_v10": ("test 6", 10, 225.0, 55.0, 150.0, 35.0),
    "shpb_v15": ("test 11", 15, 285.0, 45.0, 135.0, 40.0),
    "shpb_v20": ("test 16", 20, 400.0, 35.0, 140.0, 30.0),
}

TEMPLATE = """# ---------------------------------------------------------------------------
# NOUVEAU BANC B11 — brésilien dynamique aux barres de Hopkinson, granite de Kuru
# Essai de référence : Saksala et al. 2013 (IJRMMS 59, 128-138), {essai},
# percuteur {vstrk} m/s ; impulsion incidente {sig:.0f} MPa, montée {rise:.0f} µs,
# palier {plat:.0f} µs, descente {fall:.0f} µs (fig. 2a, lecture ±10 µs, ±10 MPa).
# Écrit par gen_decks.py (ne pas éditer à la main). Fiche : NOUVEAUX_BANCS.md, banc 2.
# Gabarit : configs_yan/shpb_complet_adaptatif.cfg (montage barre-disque-barre
# en deux phases). Déformation plane : écart au disque réel (L/D = 0,39)
# discuté dans la fiche.
# ---------------------------------------------------------------------------
mode = fdem
scenario = shpb
geometry = shpb

phases = bar rock
phase.bar.fraction  = 0.5
phase.bar.rho       = 7850
phase.bar.E         = 200e9
phase.bar.nu        = 0.29
phase.bar.ft        = 1e12
phase.bar.cohesion  = 1e12
phase.bar.Gf        = 1e6

# jeu Kuru de rockim : Yang et al. 2026, Table 1 (configs/yang2026_kuru_train1_v5.cfg)
phase.rock.fraction = 0.5
phase.rock.rho      = 2626
phase.rock.E        = 60e9
phase.rock.nu       = 0.24
phase.rock.ft       = 10.98e6
phase.rock.cohesion = 29.84e6
phase.rock.frictionDeg = 61.61
phase.rock.Gf       = 50
phase.rock.gfShearFactor = 20

rho = 2626
E = 60e9
nu = 0.24
ft = 10.98e6
cohesion = 29.84e6
frictionDeg = 61.61
Gf = 50
gfShearFactor = 20
meanTensionCapFactor = 0

# --- effet de vitesse : variante {var} ---------------------------------------
{dif}
# --- impulsion : trapèze, V0 = sigma_inc / (rho_b c_b) = {v0:.3f} m/s ---------
shpbPulse = trapezoid
shpbPulseV0 = {v0:.4f}
shpbPulseTau = {tau:.4e}
shpbPulsePlateau = {pf:.4f}
absorbFactor = 1.0

contactMu = 0.1
gcPenaltyFactor = 1.0
gcRestitution = 1.0
gcXi = 0.05
dampingLocal = 0.0
jointXi = 0.0
seed = 4211

# --- montage : barres Ø22 mm (Saksala §2.3) ; disque Ø40,8 mm (Table 1) -------
# Barre incidente portée à 1,6 m (1,2 m dans l'essai) pour que l'impulsion la
# plus longue ({taumax:.0f} µs x 5 048 m/s = {lmax:.2f} m) ne se superpose pas
# à sa réflexion sur la jauge M1 placée à 0,8 m du disque : choix numérique sans
# effet sur le chargement du disque.
shpbIncidentLength = 1.60
shpbTransmitLength = 1.20
shpbBarDiameter = 0.022
shpbDiscDiameter = 0.0408
shpbBarElemSize = 0.002
shpbDiscElemSize = 0.0005
shpbGap = 0.0
shpbMonitor1 = 0.80
shpbMonitor2 = 2.4408
gcXwindow = 0.10
frames = 30
dtFactor = 0.2
T = {T:.2e}
insertion = adaptive
insertionPenaltyFactor = 4
"""

DIF_ON = ("strainRateDIF = yang-fig2\nstrainRateTau = 1e-6\nstrainRateDIFArm = insertion")
DIF_OFF = "# aucun DIF : effet de vitesse par la seule inertie et la loi cohesive"


def main():
    taumax = max(r + p + f for (_, _, _, r, p, f) in PULSES.values())
    for name, (essai, vstrk, sig, rise, plat, fall) in PULSES.items():
        ramp = 0.5 * (rise + fall)                 # trapèze symétrique de rockim
        tau = (2 * ramp + plat) * 1e-6
        pf = plat * 1e-6 / tau
        v0 = sig * 1e6 / Z_B
        # transit 1,6 m dans l'acier (317 µs) + impulsion + 150 µs de post-rupture
        T = 1.6 / 5047.5 + tau + 150e-6
        for var, dif in (("dif", DIF_ON), ("nodif", DIF_OFF)):
            txt = TEMPLATE.format(essai=essai, vstrk=vstrk, sig=sig, rise=rise, plat=plat,
                                  fall=fall, var=var, dif=dif, v0=v0, tau=tau, pf=pf,
                                  taumax=taumax, lmax=taumax * 1e-6 * 5047.5, T=T)
            p = os.path.join(HERE, f"{name}_{var}.cfg")
            with open(p, "w", encoding="utf-8") as f:
                f.write(txt)
            print(f"{p} : V0 = {v0:.3f} m/s, tau = {tau * 1e6:.0f} us, plateau = {pf:.3f}, T = {T * 1e6:.0f} us")


if __name__ == "__main__":
    main()
