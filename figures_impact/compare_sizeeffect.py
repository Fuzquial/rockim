# -*- coding: utf-8 -*-
"""Impact 3D avec EFFET D'ECHELLE actif (eq. 42) : croise deux maillages
(moyen 82k / fin 259k, Zeff fige a 20 mm3) et deux materiaux (homogene /
Weibull m = 24 correle). Quatre runs, deux questions :
  - objectivite au maillage : la reponse depend-elle du maillage, l'effet
    d'echelle etant desormais comptabilise ?
  - apport du champ correle par-dessus la dispersion geometrique.

  python compare_sizeeffect.py
"""
import csv, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
V0 = 8.0

RUNS = [("out_imp3d_sz_homog_mid", "homogène — 82k",  "C0", "-"),
        ("out_imp3d_sz_homog_fin", "homogène — 259k", "C0", "--"),
        ("out_imp3d_sz_weib_mid",  "Weibull — 82k",   "C3", "-"),
        ("out_imp3d_sz_weib_fin",  "Weibull — 259k",  "C3", "--")]
# temoins sans effet d'echelle, pour situer
REF = [("out_imp3d_homog", "sans eff. échelle, m=6 — 82k", "0.6", ":")]


def load(sub):
    p = os.path.join(ROOT, sub, "history.csv")
    if not os.path.exists(p):
        return None
    rows = [r for r in csv.DictReader(open(p))
            if all(v not in (None, "") for v in r.values())]
    if len(rows) < 10:
        return None
    # Un run EN COURS a un history.csv partiel (le tampon ne se vide qu'au
    # remplissage) : l'inclure produit des chiffres faux — pic tronque,
    # restitution absurde. On n'accepte que les runs arrives au bout.
    if float(rows[-1]["t"]) < 0.98 * 1.2e-4:
        print("  (ignore %s : run incomplet, t = %.1f / 120 us)"
              % (sub, float(rows[-1]["t"]) * 1e6))
        return None
    g = lambda k: np.array([float(r[k]) for r in rows])
    return dict(t=g("t") * 1e6, z=g("grpZ") * 1e3, vz=g("grpVz"),
                fz=np.abs(g("grpFz")) * 1e-3, nb=g("nBroken"), last=rows[-1])


def logfac(sub):
    """facteur d'effet d'echelle mean/min/max lu dans le journal du run."""
    p = os.path.join(ROOT, "run_" + sub.replace("out_", "") + ".log")
    if not os.path.exists(p):
        return None
    s = open(p, encoding="utf-8", errors="ignore").read().replace("\x00", "")
    m = re.search(r"facteur mean/min/max = ([\d.eE+-]+)/([\d.eE+-]+)/"
                  r"([\d.eE+-]+)", s)
    return tuple(float(x) for x in m.groups()) if m else None


data = [(lab, load(sub), c, ls, logfac(sub)) for sub, lab, c, ls in RUNS]
data = [d for d in data if d[1]]
refs = [(lab, load(sub), c, ls, None) for sub, lab, c, ls in REF]
refs = [d for d in refs if d[1]]
if not data:
    raise SystemExit("aucun run trouve — les calculs tournent-ils encore ?")

fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8))

a = ax[0]
for lab, d, c, ls, _ in refs + data:
    a.plot(d["z"][0] - d["z"], d["fz"], color=c, ls=ls, lw=1.5, label=lab)
a.set_xlabel("pénétration (mm)"); a.set_ylabel(r"$|F_z|$ (kN)")
a.set_title("(a) Force – pénétration", fontsize=10)
a.legend(fontsize=7.5); a.grid(alpha=0.3)

a = ax[1]
for lab, d, c, ls, _ in refs + data:
    a.plot(d["t"], d["nb"], color=c, ls=ls, lw=1.5, label=lab)
a.set_xlabel("temps (µs)"); a.set_ylabel("joints rompus")
a.set_title("(b) Fissuration cumulée", fontsize=10)
a.grid(alpha=0.3)

a = ax[2]; a.axis("off")
head = ("", "pénétr.\n(mm)", "pic\n(kN)", "e", "joints", "facteur\nmoyen")
body = []
for lab, d, c, ls, f in refs + data:
    body.append((lab, "%.3f" % (d["z"][0] - d["z"].min()),
                 "%.1f" % d["fz"].max(), "%.3f" % (abs(d["vz"][-1]) / V0),
                 d["last"]["nBroken"], "%.4f" % f[0] if f else "—"))
tb = a.table(cellText=[list(r) for r in body], colLabels=list(head),
             loc="center", cellLoc="center")
tb.auto_set_font_size(False); tb.set_fontsize(8)
tb.scale(1, 1.9); tb.auto_set_column_width(list(range(6)))
for j in range(6):
    tb[(0, j)].set_facecolor("#e8e8e8")
a.set_title("(c) Bilan  ($v_0$ = 8 m/s, Zeff = 20 mm³, m = 24)", fontsize=10)

fig.suptitle("Impact 3D avec effet d'échelle statistique — objectivité au "
             "maillage et apport du champ corrélé", fontsize=12.5)
fig.tight_layout()
out = os.path.join(HERE, "compare_sizeeffect.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
for lab, d, c, ls, f in refs + data:
    print("%-30s pen %.3f mm  pic %5.1f kN  e %.3f  joints %-4s facteur %s"
          % (lab, d["z"][0] - d["z"].min(), d["fz"].max(),
             abs(d["vz"][-1]) / V0, d["last"]["nBroken"],
             ("%.4f" % f[0]) if f else "—"))
