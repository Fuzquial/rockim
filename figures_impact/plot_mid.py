# -*- coding: utf-8 -*-
"""Trace l'etat partiel du banc moyen complet (impact + rebond)."""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
rows = [r for r in csv.DictReader(open(os.path.join(here, "history_mid.csv")))
        if all(r.get(k) not in (None, "") for k in r)]
g = lambda k: [float(r[k]) for r in rows]
t = [x * 1e6 for x in g("t")]

fig, ax = plt.subplots(2, 2, figsize=(11.5, 8))
fig.suptitle("Banc moyen complet (82k tets, 18 threads) — état à t = %.1f / 120 µs" % t[-1])

a = ax[0][0]
a.plot(t, [f * 1e-3 for f in g("grpFz")], "C2")
a.axvline(6.25, color="gray", ls="--", lw=0.8)
a.annotate("contact", (6.25, 0.5), fontsize=8, color="gray")
a.set_xlabel("t (µs)"); a.set_ylabel("grpFz (kN)")
a.set_title("Force de contact nette sur l'insert")

a = ax[0][1]
a.plot(t, g("grpVz"), "C1")
a.axhline(0, color="gray", lw=0.6)
a.set_xlabel("t (µs)"); a.set_ylabel("grpVz (m/s)")
a.set_title("Vitesse de l'insert (inversion = rebond)")

a = ax[1][0]
a.plot(t, g("nBroken"), "C3", label="joints cassés")
a2 = a.twinx()
a2.plot(t, g("nFrag"), "C0", ls="--", label="fragments")
a.set_xlabel("t (µs)"); a.set_ylabel("joints cassés", color="C3")
a2.set_ylabel("fragments (corps compris)", color="C0")
a.set_title("Endommagement et fragmentation")

a = ax[1][1]
for k, c in [("eEl", "C0"), ("eJnt", "C1"), ("eGc", "C2"),
             ("eFric", "C4"), ("eCund", "C5"), ("eLys", "C6")]:
    a.plot(t, g(k), c, label=k)
a.set_xlabel("t (µs)"); a.set_ylabel("énergie (J)")
a.set_title("Bilan B4 par sous-système")
a.legend(ncol=2, fontsize=8)

fig.tight_layout()
out = os.path.join(here, "banc_mid_progress.png")
fig.savefig(out, dpi=150)
vmin = min(g("grpVz")); vlast = rows[-1]["grpVz"]
print(out, "| Fz max = %.1f kN" % (max(g("grpFz")) * 1e-3),
      "| vz min/actuel = %.2f / %s" % (vmin, vlast),
      "| nBroken =", rows[-1]["nBroken"], "| nFrag =", rows[-1]["nFrag"])
