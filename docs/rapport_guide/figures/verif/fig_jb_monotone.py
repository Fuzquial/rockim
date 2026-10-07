"""Banc d'un joint isole (scenario = jointbench, S3) : traction et cisaillement
monotones avec decharge, pour les trois branches de decharge.
Donnees : jb_monotone.csv (rejoue le 2026-10-06, prep_donnees_verif.py).
Usage : python fig_jb_monotone.py -> fig_jb_monotone.pdf
"""
import os, csv
import numpy as np
from style_verif import plt, COUL, NOM, virgule

HERE = os.path.dirname(os.path.abspath(__file__))
FT, C = 10.98e6, 29.84e6
rows = list(csv.DictReader(open(os.path.join(HERE, "jb_monotone.csv"))))

fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.6))
styles = {"plastic": "-", "origin": (0, (4, 2)), "solidity": ":"}
for law in ("plastic", "origin", "solidity"):
    for k, mode in enumerate(("traction", "cisaillement")):
        r = [x for x in rows if x["loi"] == law and x["trajet"] == mode]
        dn = np.array([float(x["dn_m"]) for x in r])
        ds = np.array([float(x["ds_m"]) for x in r])
        s = np.array([float(x["sig_Pa"]) for x in r])
        t = np.array([float(x["tau_Pa"]) for x in r])
        if mode == "traction":
            ax[0].plot(dn * 1e6, s / FT, color=COUL[law], ls=styles[law],
                       lw=1.6 if law == "solidity" else 1.1, label=NOM[law])
        else:
            ax[1].plot(ds * 1e6, t / C, color=COUL[law], ls=styles[law],
                       lw=1.6 if law == "solidity" else 1.1, label=NOM[law])
ax[0].set_xlabel(r"ouverture $\delta_n$ ($\mu$m)")
ax[0].set_ylabel(r"$\sigma / f_t$")
ax[0].set_title(r"(a) traction, décharge à 10 $\mu$m", fontsize=9)
ax[1].set_xlabel(r"glissement $\delta_s$ ($\mu$m)")
ax[1].set_ylabel(r"$\tau / c$")
ax[1].set_title(r"(b) cisaillement sous $\sigma_n = -60$ MPa, décharge à 12 $\mu$m", fontsize=9)
for a in ax:
    a.axhline(0, color="k", lw=0.4)
    a.legend(loc="upper right")
    virgule(a)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_jb_monotone.pdf"))
