"""Banc d'un joint isole, cycle ferme a pression variable (jbMode = cycle, S3 bis) :
trajet impose et boucles tau-glissement du troisieme cycle, cycle intact (D = 0)
et cycle endommage, pour les trois branches de decharge.
Donnees : jb_cycles.csv (rejoue le 2026-10-06, prep_donnees_verif.py).
Usage : python fig_jb_cycles.py -> fig_jb_cycles.pdf
"""
import os, csv
import numpy as np
from style_verif import plt, COUL, NOM, virgule

HERE = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(HERE, "jb_cycles.csv"))))


def get(regime, law):
    r = [x for x in rows if x["regime"] == regime and x["loi"] == law]
    return {k: np.array([float(x[k]) for x in r]) for k in ("dn_m", "ds_m", "sig_Pa", "tau_Pa", "phase")}


fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw=dict(width_ratios=[0.8, 1, 1]))
# (a) trajet impose, cycle intact (le cycle endommage a la meme forme, x 133 en s, x 10 en dn)
d = get("intact", "plastic")
a = ax[0]
a.plot(d["ds_m"] * 1e9, d["dn_m"] * 1e9, color="k")
for lab, (x, y), off in (("A", (7.5, -2), (0, 1.2)), ("B", (15, -11), (1.0, 0)),
                         ("C", (7.5, -20), (0, -2.6)), ("D", (0, -11), (-2.6, 0))):
    a.annotate(lab, (x + off[0], y + off[1]), ha="center", va="center", fontsize=9)
a.annotate("", xy=(11, -2), xytext=(4, -2), arrowprops=dict(arrowstyle="->", lw=0.8))
a.set_xlabel(r"glissement $\delta_s$ (nm)")
a.set_ylabel(r"ouverture $\delta_n$ (nm)")
a.set_title("(a) trajet imposé, cycle intact", fontsize=9)
a.set_xlim(-4, 19); a.set_ylim(-24, 2)
virgule(a)
styles = {"plastic": "-", "origin": (0, (4, 2)), "solidity": "-"}
for k, (regime, sx, sy, ux, uy) in enumerate((("intact", 1e9, 1e-6, "nm", "MPa"),
                                               ("endommage", 1e6, 1e-6, r"$\mu$m", "MPa"))):
    a = ax[k + 1]
    for law in ("solidity", "plastic", "origin"):
        d = get(regime, law)
        a.plot(d["ds_m"] * sx, d["tau_Pa"] * sy, color=COUL[law], ls=styles[law],
               lw=1.8 if law == "solidity" else 1.0, label=NOM[law])
    a.set_xlabel(r"glissement $\delta_s$ (%s)" % ux)
    a.set_ylabel(r"$\tau$ (%s)" % uy)
    a.set_title("(%s) %s" % ("bc"[k], "cycle intact, $D = 0$" if regime == "intact"
                               else r"cycle endommagé, $D_{\max}$ 0,47 à 0,88"), fontsize=9)
    a.legend(loc="upper left")
    virgule(a)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_jb_cycles.pdf"))
