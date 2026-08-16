# -*- coding: utf-8 -*-
"""Test du « facteur 2 » : la deformation simulee, telle quelle et doublee,
comparee aux essais. Verifie si l'ecart sur eps_pic peut etre une convention
de mesure (facteur 2) ou s'il est PHYSIQUE (amorcage premature)."""
import csv, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
d = json.load(open(os.path.join(BASE, "targets", "curves_redbohus.json")))


def hist(run):
    p = os.path.join(BASE, "runs", run, "history.csv")
    rows = [r for r in csv.DictReader(open(p))
            if all(v not in (None, "") for v in r.values())]
    s = np.abs(np.array([float(r["sigma"]) for r in rows])) / 1e6
    e = np.abs(np.array([float(r["epsSpec"]) for r in rows])) * 100
    ep = np.abs(np.array([float(r["epsPlaten"]) for r in rows])) * 100
    eg = np.abs(np.array([float(r["epsGauge"]) for r in rows])) * 100
    return s, e, ep, eg


fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.4))

for a, mult, ttl in ((ax[0], 1.0, "déformation simulée TELLE QUELLE"),
                     (ax[1], 2.0, "déformation simulée × 2")):
    for i, (k, s) in enumerate(d["triaxial"].items()):
        if s["sigma3_MPa"] != 20:
            continue
        a.plot(s["eps_axial_pct"], np.array(s["q_MPa"]) + 20, color="0.45",
               lw=1.4, ls="--", label="essais réels (3)" if i == 0 else None)
    for c, tag in (("C1", "E005"), ("C2", "E008")):
        sg, e, _, _ = hist("%s_tx20_s4211" % tag)
        a.plot(e * mult, sg, color=c, lw=1.7, label=tag)
    a.set_xlabel("déformation axiale (%)"); a.set_ylabel(r"$\sigma_1$ (MPa)")
    a.set_title(ttl); a.grid(alpha=0.3); a.legend(fontsize=9)
    a.set_xlim(0, 1.3); a.set_ylim(0, 460)

fig.suptitle("Triaxial $\\sigma_3$ = 20 MPa — le décalage de $\\varepsilon_{pic}$ est-il "
             "un facteur 2 de convention, ou de la physique ?\n"
             "à droite, doubler la déformation aligne le PIC mais divise la "
             "RAIDEUR par deux : l'accord élastique est détruit", fontsize=11)
fig.tight_layout()
out = os.path.join(BASE, "figures", "test_facteur2.png")
fig.savefig(out, dpi=140)
print("ecrit", out)

# --- controle chiffre : les trois extensometres du solveur ----------------
print("\nles trois mesures de deformation du solveur (au pic), en %% :")
print("%-6s %10s %10s %10s" % ("jeu", "epsPlaten", "epsSpec", "epsGauge"))
for tag in ("E005", "E008"):
    s, e, ep, eg = hist("%s_tx20_s4211" % tag)
    i = int(np.argmax(s))
    print("%-6s %10.3f %10.3f %10.3f" % (tag, ep[i], e[i], eg[i]))
epx = [np.array(v["eps_axial_pct"])[int(np.argmax(v["q_MPa"]))]
       for v in d["triaxial"].values() if v["sigma3_MPa"] == 20]
print("essais reels : %.3f %%" % np.mean(epx))
