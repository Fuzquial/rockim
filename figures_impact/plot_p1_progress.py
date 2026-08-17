# -*- coding: utf-8 -*-
"""Trace l'etat partiel du run P1 banc (1 thread) depuis history.csv."""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(here, "history_snapshot.csv"))))
# la derniere ligne peut etre tronquee (fichier ecrit en direct)
rows = [r for r in rows if all(r.get(k) not in (None, "") for k in r)]
g = lambda k: [float(r[k]) for r in rows]
t = [x * 1e6 for x in g("t")]  # µs

fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
fig.suptitle("P1 banc 842k tets — run 1 thread, état partiel (t = %.2f / 20 µs)" % t[-1])

a = ax[0][0]
a.plot(t, [z * 1e3 for z in g("grpZ")], "C0")
a.set_xlabel("t (µs)"); a.set_ylabel("grpZ insert (mm)")
a.set_title("Position verticale de l'insert")

a = ax[0][1]
a.plot(t, g("grpVz"), "C1")
a.set_xlabel("t (µs)"); a.set_ylabel("grpVz (m/s)")
a.set_title("Vitesse verticale de l'insert")

a = ax[1][0]
a.plot(t, [f * 1e-3 for f in g("grpFz")], "C2", label="grpFz (kN)")
a.plot(t, [s * 1e-6 for s in g("grpSzz")], "C3", label="grpSzz (MPa)")
a.set_xlabel("t (µs)"); a.set_title("Force de contact nette / jauge")
a.legend()

a = ax[1][1]
for k, c in [("eEl", "C0"), ("eJnt", "C1"), ("eGc", "C2"),
             ("eFric", "C4"), ("eCund", "C5"), ("eLys", "C6")]:
    a.plot(t, g(k), c, label=k)
a.set_xlabel("t (µs)"); a.set_ylabel("énergie (J)")
a.set_title("Bilan B4 (nBroken = %d)" % int(float(rows[-1]["nBroken"])))
a.legend(ncol=2, fontsize=8)

fig.tight_layout()
out = os.path.join(here, "p1_progress.png")
fig.savefig(out, dpi=150)
print(out, "|", len(rows), "pas |  t_max =", rows[-1]["t"],
      "| grpZ:", rows[0]["grpZ"], "->", rows[-1]["grpZ"],
      "| grpVz:", rows[-1]["grpVz"], "| grpFz:", rows[-1]["grpFz"],
      "| nBroken:", rows[-1]["nBroken"])
