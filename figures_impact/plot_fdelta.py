# -*- coding: utf-8 -*-
"""Courbe force-penetration (F-delta) du banc moyen, charge et decharge.
delta obtenu par integration de grpVz (lisse), reference au plan de contact."""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
rows = [r for r in csv.DictReader(open(os.path.join(here, "history_mid.csv")))
        if all(r.get(k) not in (None, "") for k in r)]
t = np.array([float(r["t"]) for r in rows])
vz = np.array([float(r["grpVz"]) for r in rows])
Fz = np.array([float(r["grpFz"]) for r in rows]) * 1e-3   # kN

# position du centre par integration (evite l'escalier 1 um de grpZ)
z0 = float(rows[0]["grpZ"])
z = z0 + np.concatenate(([0.0], np.cumsum(0.5 * (vz[1:] + vz[:-1]) * np.diff(t))))
zc = 0.12 + 0.011                    # centre a l'affleurement (bloc + rayon)
delta = np.maximum(0.0, (zc - z)) * 1e3   # mm

ipk = int(np.argmax(Fz))
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(delta[:ipk + 1], Fz[:ipk + 1], "C0", lw=1.8, label="charge")
ax.plot(delta[ipk:], Fz[ipk:], "C3", lw=1.8, label="décharge (rebond)")
ax.annotate("pic : %.1f kN à δ = %.3f mm" % (Fz[ipk], delta[ipk]),
            (delta[ipk], Fz[ipk]), textcoords="offset points", xytext=(-10, 8),
            ha="right", fontsize=9)
ax.set_xlabel("pénétration δ (mm)")
ax.set_ylabel("force de contact Fz (kN)")
ax.set_title("F–δ banc moyen (82k tets, insert R = 11 mm, v = 8 m/s) — "
             "état à t = %.1f µs" % (t[-1] * 1e6))
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
out = os.path.join(here, "banc_mid_fdelta.png")
fig.savefig(out, dpi=150)
aire_ch = np.trapezoid(Fz[:ipk + 1] * 1e3, delta[:ipk + 1] * 1e-3)
aire_de = -np.trapezoid(Fz[ipk:] * 1e3, delta[ipk:] * 1e-3)
print(out, "| delta_max = %.3f mm" % delta.max(),
      "| W_charge = %.2f J, W_restitue = %.2f J (%.0f %%)"
      % (aire_ch, aire_de, 100 * aire_de / aire_ch))
