# -*- coding: utf-8 -*-
"""Les courbes du run ultra (ft = 87 MPa, l_cz = 1,0 mm) face aux runs rates
de la journee (ft = 10 MPa, l_cz = 35 mm) : force-penetration, cinematique de
l'insert, fissuration, et repartition de l'energie.

  python plot_ultra_courbes.py
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
V0, KE0 = 8.0, 2.53

RUNS = [("out_imp3d_homog", "ft = 10 MPa — $\\ell_{cz}$ = 35 mm", "0.55", "--"),
        ("out_imp3d_ultra", "ft = 87 MPa — $\\ell_{cz}$ = 1,0 mm", "C3", "-")]


def load(sub):
    p = os.path.join(ROOT, sub, "history.csv")
    if not os.path.exists(p):
        return None
    rows = [r for r in csv.DictReader(open(p))
            if all(v not in (None, "") for v in r.values())]
    if len(rows) < 30:
        return None
    g = lambda k: np.array([float(r[k]) for r in rows])
    z = g("grpZ") * 1e3
    return dict(t=g("t") * 1e6, pen=z[0] - z, fz=np.abs(g("grpFz")) * 1e-3,
                vz=g("grpVz"), nb=g("nBroken"), jnt=np.abs(g("eJnt")),
                cund=np.abs(g("eCund")), lys=np.abs(g("eLys")),
                gc=np.abs(g("eGc")), fric=np.abs(g("eFric")),
                complet=g("t")[-1] >= 0.95 * 1.2e-4)


data = [(lab, load(s), c, ls) for s, lab, c, ls in RUNS]
data = [d for d in data if d[1]]

fig, ax = plt.subplots(2, 2, figsize=(13.5, 9.0))

# (a) force - penetration
a = ax[0, 0]
for lab, d, c, ls in data:
    m = d["pen"] > 0
    a.plot(d["pen"][m], d["fz"][m], color=c, ls=ls, lw=1.9, label=lab)
a.set_xlabel("pénétration (mm)"); a.set_ylabel(r"$|F_z|$ (kN)")
a.set_title("(a) Force – pénétration", fontsize=10.5)
a.legend(fontsize=9); a.grid(alpha=0.3)

# (b) vitesse de l'insert : le rebond se lit au passage par zero
a = ax[0, 1]
a.axhline(0, color="k", lw=1.2)
a.axhspan(0, 8, color="0.93", zorder=0)
a.text(3, 1.2, "au-dessus de 0 = REBOND", fontsize=8.5, color="0.35")
for lab, d, c, ls in data:
    a.plot(d["t"], d["vz"], color=c, ls=ls, lw=1.9, label=lab)
a.set_xlabel("temps (µs)"); a.set_ylabel(r"$v_z$ de l'insert (m/s)")
a.set_title("(b) Cinématique de l'insert  ($v_0$ = −8 m/s)", fontsize=10.5)
a.legend(fontsize=9); a.grid(alpha=0.3)

# (c) fissuration cumulee
a = ax[1, 0]
for lab, d, c, ls in data:
    a.plot(d["t"], d["nb"], color=c, ls=ls, lw=1.9, label=lab)
a.set_xlabel("temps (µs)"); a.set_ylabel("joints rompus")
a.set_yscale("symlog", linthresh=10)
a.set_title("(c) Fissuration cumulée (échelle log)", fontsize=10.5)
a.legend(fontsize=9); a.grid(alpha=0.3)

# (d) ou va l'energie : fissuration contre amortissement numerique
a = ax[1, 1]
for lab, d, c, ls in data:
    a.plot(d["t"], d["jnt"], color=c, ls=ls, lw=2.1,
           label="fissuration — %s" % lab.split("—")[0].strip())
    a.plot(d["t"], d["cund"], color=c, ls=":", lw=1.5,
           label="amortissement Cundall — %s" % lab.split("—")[0].strip())
a.set_xlabel("temps (µs)"); a.set_ylabel("énergie (J)")
a.set_title("(d) Fissuration contre amortissement NUMÉRIQUE", fontsize=10.5)
a.legend(fontsize=8); a.grid(alpha=0.3)

partiel = any(not d["complet"] for _, d, _, _ in data)
fig.suptitle("Le couple (ft, Gf) est le seul verrou — impact insert/granite "
             "à 8 m/s%s" % ("   [RUN ULTRA EN COURS]" if partiel else ""),
             fontsize=13)
fig.tight_layout()
out = os.path.join(HERE, "ultra_courbes.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
for lab, d, c, ls in data:
    r = d["jnt"][-1] / d["cund"][-1] if d["cund"][-1] > 0 else float("nan")
    print("  %-34s t %.0f us | pen %.3f mm | pic %.1f kN | vz %+.2f | joints %d"
          " | fissu/amort %.2f" % (lab.replace("$\\ell_{cz}$", "l_cz"),
          d["t"][-1], d["pen"].max(), d["fz"].max(), d["vz"][-1], d["nb"][-1], r))
