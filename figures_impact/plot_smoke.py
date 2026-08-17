# -*- coding: utf-8 -*-
"""Trace le smoke test impact multithread (out_smoke/history.csv)."""
import csv, os, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

src = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "out_smoke", "history.csv")
rows = [r for r in csv.DictReader(open(src))
        if all(r.get(k) not in (None, "") for k in r)]
g = lambda k: [float(r[k]) for r in rows]
t = [x * 1e6 for x in g("t")]

fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
fig.suptitle("Smoke test impact bench1 (17,9k tets, gap 0,05 mm, 18 threads) — wall 158 s")

a = ax[0][0]
a.plot(t, [f * 1e-3 for f in g("grpFz")], "C2")
a.axvline(6.25, color="gray", ls="--", lw=0.8)
a.set_xlabel("t (µs)"); a.set_ylabel("grpFz (kN)")
a.set_title("Force de contact nette sur l'insert (F-δ en direct, B2)")

a = ax[0][1]
a.plot(t, g("grpVz"), "C1")
a.axvline(6.25, color="gray", ls="--", lw=0.8)
a.set_xlabel("t (µs)"); a.set_ylabel("grpVz (m/s)")
a.set_title("Vitesse verticale de l'insert")

a = ax[1][0]
a.plot(t, g("nBroken"), "C3")
a.set_xlabel("t (µs)"); a.set_ylabel("joints casses")
a.set_title("Endommagement (nBroken)")

a = ax[1][1]
for k, c in [("eEl", "C0"), ("eJnt", "C1"), ("eGc", "C2"),
             ("eFric", "C4"), ("eCund", "C5"), ("eLys", "C6")]:
    a.plot(t, g(k), c, label=k)
a.set_xlabel("t (µs)"); a.set_ylabel("energie (J)")
a.set_title("Bilan B4 par sous-systeme")
a.legend(ncol=2, fontsize=8)

fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "smoke_impact.png")
fig.savefig(out, dpi=150)
imp = next((r for r in rows if abs(float(r["grpFz"])) > 1.0), None)
print(out, "| contact des t =", (imp["t"] if imp else "jamais"),
      "| Fz max =", max(abs(f) for f in g("grpFz")),
      "| vz fin =", rows[-1]["grpVz"], "| nBroken fin =", rows[-1]["nBroken"])
