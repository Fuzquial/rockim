# -*- coding: utf-8 -*-
"""Impact 3D : le run du 14/08 (sans insertion adaptative) contre les deux cas
du 17/08 — homogene et heterogene (Weibull), panoplie complete.
Trois runs, meme maillage, meme insert a 8 m/s.

  python compare_adap_homog_weib.py
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
V0 = 8.0

RUNS = [
    ("out_banc_mid",     "14/08 — insertion intrinsèque", "0.45", "--"),
    ("out_imp3d_homog",  "17/08 — adaptative, homogène",  "C0",   "-"),
    ("out_imp3d_weib",   "17/08 — adaptative + Weibull",  "C3",   "-"),
]


def load(sub):
    p = os.path.join(ROOT, sub, "history.csv")
    if not os.path.exists(p):
        return None
    rows = [r for r in csv.DictReader(open(p))
            if all(v not in (None, "") for v in r.values())]
    if not rows:
        return None
    g = lambda k: np.array([float(r[k]) for r in rows])
    return dict(t=g("t") * 1e6, z=g("grpZ") * 1e3, vz=g("grpVz"),
                fz=np.abs(g("grpFz")) * 1e-3, nb=g("nBroken"),
                nf=g("nFrag"), last=rows[-1])


data = [(lab, load(sub), c, ls) for sub, lab, c, ls in RUNS]
data = [(l, d, c, ls) for l, d, c, ls in data if d]

fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.7))

# (a) force - penetration
a = ax[0]
for lab, d, c, ls in data:
    pen = d["z"][0] - d["z"]
    a.plot(pen, d["fz"], color=c, ls=ls, lw=1.5, label=lab)
a.set_xlabel("pénétration (mm)"); a.set_ylabel(r"force $|F_z|$ (kN)")
a.set_title("(a) Force – pénétration", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3)

# (b) joints rompus
a = ax[1]
for lab, d, c, ls in data:
    a.plot(d["t"], d["nb"], color=c, ls=ls, lw=1.5, label=lab)
a.set_xlabel("temps (µs)"); a.set_ylabel("joints rompus")
a.set_title("(b) Fissuration cumulée", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3)

# (c) bilan chiffre
a = ax[2]; a.axis("off")
lignes = [("", "pénétr.\n(mm)", "pic\n(kN)", "e", "joints", "frag.")]
for lab, d, c, ls in data:
    pen = d["z"][0] - d["z"].min()
    lignes.append((lab.split(" — ")[1] if " — " in lab else lab,
                   "%.3f" % pen, "%.1f" % d["fz"].max(),
                   "%.3f" % (abs(d["vz"][-1]) / V0),
                   d["last"]["nBroken"], d["last"]["nFrag"]))
tb = a.table(cellText=[list(r) for r in lignes[1:]],
             colLabels=list(lignes[0]), loc="center", cellLoc="center")
tb.auto_set_font_size(False); tb.set_fontsize(8.5); tb.scale(1, 1.9)
tb.auto_set_column_width(list(range(6)))
for j in range(6):
    tb[(0, j)].set_facecolor("#e8e8e8")
a.set_title("(c) Bilan  (e = restitution, $v_0$ = 8 m/s)", fontsize=10)

fig.suptitle("Impact 3D insert/granite — effet de la panoplie et de "
             "l'hétérogénéité", fontsize=12.5)
fig.tight_layout()
out = os.path.join(HERE, "compare_adap_homog_weib.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
for lab, d, c, ls in data:
    print("%-38s pen %.3f mm  pic %5.1f kN  e %.3f  joints %s  frag %s"
          % (lab, d["z"][0] - d["z"].min(), d["fz"].max(),
             abs(d["vz"][-1]) / V0, d["last"]["nBroken"], d["last"]["nFrag"]))
