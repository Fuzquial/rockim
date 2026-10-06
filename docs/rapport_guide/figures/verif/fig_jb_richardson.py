"""Energie nette par cycle ferme contre le pas d'echantillonnage (S3 bis) :
convergence en pas de temps et extrapolation de Richardson d'ordre 1,
W_0 = 2 W(dt/2) - W(dt). Marqueur plein : W > 0 (creation) ; creux : W < 0.
Traits horizontaux : limite de Richardson. Donnees : jb_energie_cycle.csv.
Usage : python fig_jb_richardson.py -> fig_jb_richardson.pdf
"""
import os, csv
import numpy as np
from style_verif import plt, COUL, NOM, fmt

HERE = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(HERE, "jb_energie_cycle.csv"))))
mk = {"plastic": "o", "origin": "s", "solidity": "D"}
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7))
for k, regime in enumerate(("intact", "endommage")):
    a = ax[k]
    txt = []
    for law in ("plastic", "origin", "solidity"):
        r = sorted([x for x in rows if x["regime"] == regime and x["loi"] == law],
                   key=lambda x: -float(x["dt_s"]))
        dt = np.array([float(x["dt_s"]) for x in r])
        W = np.array([float(x["W_journal_J"]) for x in r])
        W0 = 2 * W[1] - W[0]
        W0b = 2 * W[2] - W[1]
        pos = W > 0
        off = 1.06 if law == "origin" else 1.0       # origin et plastic confondus a D = 0
        a.loglog(dt * off, np.abs(W), color=COUL[law], lw=0.8)
        if pos.any():
            a.scatter(dt[pos] * off, np.abs(W[pos]), marker=mk[law], color=COUL[law], s=22,
                      zorder=3, label=NOM[law] + r", $W>0$")
        if (~pos).any():
            a.scatter(dt[~pos] * off, np.abs(W[~pos]), marker=mk[law], facecolor="white",
                      edgecolor=COUL[law], s=22, zorder=3, label=NOM[law] + r", $W<0$")
        if abs(W0) > 1e-3 * abs(W[0]):
            a.axhline(abs(W0), color=COUL[law], lw=0.7, ls=":")
            m, e = ("%.4e" % W0).split("e")
            txt.append((law, r"$W_0 = %s\times10^{%d}$ J" % (fmt(float(m), 4).replace(",", "{,}"), int(e))))
    for j, (law, t) in enumerate(txt):
        a.text(0.04, 0.60 - 0.09 * j, NOM[law] + " : " + t, transform=a.transAxes,
               fontsize=7, color=COUL[law])
    a.set_xlabel(r"pas d'échantillonnage $\Delta t$ (s)")
    a.set_ylabel(r"$|W|$ par cycle (J)")
    a.set_title("(%s) %s" % ("ab"[k], "cycle intact, $D = 0$" if regime == "intact"
                               else "cycle endommagé"), fontsize=9)
    a.legend(loc="lower right" if regime == "intact" else "center right", fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_jb_richardson.pdf"))
