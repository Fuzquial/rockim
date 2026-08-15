# -*- coding: utf-8 -*-
"""Courbes contrainte-deformation : simulations rockim vs essais reels
Red Bohus (Dumoulin et al. 2024). UCS, bresilien et triaxial."""
import csv, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
RUNS = os.path.join(BASE, "runs")
EXP = (r"C:\Users\fuzquianoalricabi\OneDrive - Université Paris Sciences "
       r"et Lettres\Documents\phd_geothermie\FDEM\calib_cdp_bohus"
       r"\experimental_data_red_bohus.json")
TAGS = sys.argv[1:] if len(sys.argv) > 1 else ["L001", "L011", "L029", "L035"]
COL = ["C0", "C1", "C2", "C4", "C5", "C6"]


def hist(run):
    p = os.path.join(RUNS, run, "history.csv")
    if not os.path.exists(p):
        return None
    rows = [r for r in csv.DictReader(open(p))
            if all(v not in (None, "") for v in r.values())]
    if not rows:
        return None
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


fig, ax = plt.subplots(1, 3, figsize=(15, 5.2))

# ---------------------------------------------------------------- UCS
a = ax[0]
# courbes EXPERIMENTALES PROPRES (extract_curves.py, protocole des scripts de
# depouillement de la these) : deformation GLOBALE, non filtree — les jauges
# locales decrochent et donnent des courbes qui reviennent en arriere
CURV = os.path.join(BASE, "targets", "curves_redbohus.json")
d = json.load(open(CURV))
for i, (k, s) in enumerate(d["UC"].items()):
    a.plot(s["eps_global_pct"], s["stress_MPa"], color="0.5", lw=1.3, ls="--",
           label="essais réels (4)" if i == 0 else None)
for c, tag in zip(COL, TAGS):
    h = hist("%s_ucs_s4211" % tag)
    if h is None:
        continue
    a.plot(np.abs(h["epsSpec"]) * 100, np.abs(h["sigma"]) * 1e-6, color=c,
           lw=1.6, label="%s (pic %.0f MPa)" % (tag, np.abs(h["sigma"]).max() / 1e6))
a.axhline(126.6, color="C2", lw=1, ls=":")
a.set_xlim(0, 1.2); a.set_ylim(0, 260)
a.set_xlabel("déformation axiale (%)"); a.set_ylabel(r"$\sigma_1$ (MPa)")
a.set_title("Compression simple — rockim vs essais Red Bohus")
a.legend(fontsize=8); a.grid(alpha=0.3)

# ---------------------------------------------------------------- BTS
a = ax[1]
for c, tag in zip(COL, TAGS):
    h = hist("%s_bts_s4211" % tag)
    if h is None:
        continue
    a.plot(h["drive"] * 1e3, h["sigmaT"] * 1e-6, color=c, lw=1.6,
           label="%s (%.1f MPa)" % (tag, h["sigmaT"].max() / 1e6))
a.axhline(10.27, color="C2", lw=1.5, label="cible exp. 10,3 MPa")
a.set_xlabel("rapprochement des plateaux (mm)")
a.set_ylabel(r"$\sigma_t = 2P/\pi Dt$ (MPa)")
a.set_title("Brésilien")
a.legend(fontsize=8); a.grid(alpha=0.3)

# ---------------------------------------------------------------- triaxial
a = ax[2]
for c, tag in zip(COL, TAGS):
    h = hist("%s_tx20_s4211" % tag)
    if h is None:
        continue
    a.plot(np.abs(h["epsSpec"]) * 100, np.abs(h["sigma"]) * 1e-6, color=c,
           lw=1.6, label="%s (pic %.0f MPa)" % (tag, np.abs(h["sigma"]).max() / 1e6))
for i, (k, s) in enumerate(d["triaxial"].items()):
    if s["sigma3_MPa"] != 20:
        continue
    a.plot(s["eps_axial_pct"], np.array(s["q_MPa"]) + 20.0, color="0.5",
           lw=1.3, ls="--", label="essais réels (3)" if i == 0 else None)
a.axhline(424.8, color="C2", lw=1.5, label="pic exp. 424,8 MPa")
a.set_xlabel("déformation axiale (%)"); a.set_ylabel(r"$\sigma_1$ (MPa)")
a.set_xlim(0, 1.2)
a.set_title(r"Triaxial $\sigma_3$ = 20 MPa")
a.legend(fontsize=8); a.grid(alpha=0.3)

fig.suptitle("Courbes contrainte-déformation de la campagne de calibration "
             "Red Bohus (jeux de la base LHS)", fontsize=12)
fig.tight_layout()
out = os.path.join(BASE, "figures", "courbes_campagne.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
