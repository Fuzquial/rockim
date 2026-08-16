# -*- coding: utf-8 -*-
"""Controle en cours de campagne : les runs relances atteignent-ils un VRAI
pic, et l'anomalie sigma1(50) < sigma1(20) est-elle bien une troncature ?"""
import csv, glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
RUNS = os.path.join(BASE, "runs")


def hist(run):
    p = os.path.join(RUNS, run, "history.csv")
    if not os.path.exists(p):
        return None
    rows = [r for r in csv.DictReader(open(p))
            if all(v not in (None, "") for v in r.values())]
    if len(rows) < 10:
        return None
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def locked(run):
    lg = os.path.join(RUNS, run + ".log")
    if not os.path.exists(lg):
        return None
    t = open(lg, errors="replace").read()
    if "peak LOCKED" in t:
        return True
    if "NOT locked" in t:
        return False
    return None


# --- jeux termines --------------------------------------------------------
done = []
for lg in sorted(glob.glob(os.path.join(RUNS, "L*_tx20_s4211.log"))
                 + glob.glob(os.path.join(RUNS, "E*_tx20_s4211.log"))):
    tag = os.path.basename(lg).split("_")[0]
    if os.path.exists(os.path.join(RUNS, "%s_tx50_s4211.log" % tag)):
        done.append(tag)
print("jeux avec tx20 ET tx50 termines : %d" % len(done))

fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))

# (a) courbes des jeux termines
a = ax[0]
nlock = {20: [0, 0], 50: [0, 0]}
for tag in done:
    for s3, c in ((20, "C0"), (50, "C3")):
        h = hist("%s_tx%d_s4211" % (tag, s3))
        if h is None:
            continue
        lk = locked("%s_tx%d_s4211" % (tag, s3))
        nlock[s3][0 if lk else 1] += 1
        a.plot(np.abs(h["epsGauge"]) * 100, np.abs(h["sigma"]) * 1e-6,
               color=c, lw=1.0, alpha=0.55,
               ls="-" if lk else ":")
a.plot([], [], "C0-", label="$\\sigma_3$=20, pic verrouillé")
a.plot([], [], "C0:", label="$\\sigma_3$=20, NON verrouillé")
a.plot([], [], "C3-", label="$\\sigma_3$=50, pic verrouillé")
a.plot([], [], "C3:", label="$\\sigma_3$=50, NON verrouillé")
a.axvline(0.85, color="k", lw=1.4, ls=":")
a.text(0.87, 30, "ancienne\nfenêtre", fontsize=7, rotation=90)
a.set_xlabel("déformation axiale (%)"); a.set_ylabel(r"$\sigma_1$ (MPa)")
a.set_title("Courbes des jeux relancés"); a.legend(fontsize=7)
a.grid(alpha=0.3); a.set_xlim(0, 1.55)

# (b) coherence sigma1(50) > sigma1(20) ?
a = ax[1]
p20, p50, bad = [], [], []
for tag in done:
    h2, h5 = hist("%s_tx20_s4211" % tag), hist("%s_tx50_s4211" % tag)
    if h2 is None or h5 is None:
        continue
    v2 = np.abs(h2["sigma"]).max() * 1e-6
    v5 = np.abs(h5["sigma"]).max() * 1e-6
    p20.append(v2); p50.append(v5)
    if v5 < v2:
        bad.append(tag)
p20, p50 = np.array(p20), np.array(p50)
ok = p50 >= p20
a.plot(p20[ok], p50[ok], "o", color="C0", ms=6, label="cohérent")
a.plot(p20[~ok], p50[~ok], "s", color="C3", ms=7,
       label="INCOHÉRENT ($\\sigma_1$(50) < $\\sigma_1$(20))")
lim = [0, max(p20.max(), p50.max()) * 1.05]
a.plot(lim, lim, "k--", lw=0.8)
a.plot(424.8, 649.0, "*", color="C2", ms=18, label="expérimental")
a.set_xlabel(r"$\sigma_1$ à $\sigma_3$ = 20 (MPa)")
a.set_ylabel(r"$\sigma_1$ à $\sigma_3$ = 50 (MPa)")
a.set_title("Cohérence de l'enveloppe"); a.legend(fontsize=8); a.grid(alpha=0.3)

# (c) etat des pics
a = ax[2]
lab = ["$\\sigma_3$=20\nverrouillé", "$\\sigma_3$=20\nnon verr.",
       "$\\sigma_3$=50\nverrouillé", "$\\sigma_3$=50\nnon verr."]
val = [nlock[20][0], nlock[20][1], nlock[50][0], nlock[50][1]]
a.bar(range(4), val, color=["C0", "C7", "C3", "C7"])
for i, v in enumerate(val):
    a.text(i, v + 0.3, str(v), ha="center", fontsize=10)
a.set_xticks(range(4)); a.set_xticklabels(lab, fontsize=8)
a.set_ylabel("nombre de runs")
a.set_title("Les pics sont-ils atteints ?"); a.grid(alpha=0.3, axis="y")

fig.suptitle("Contrôle à mi-campagne — %d jeux relancés avec T = 6,5·10⁻³ s"
             % len(done), fontsize=12)
fig.tight_layout()
out = os.path.join(BASE, "figures", "controle_rerun.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
print("pics verrouilles : s3=20 -> %d/%d, s3=50 -> %d/%d"
      % (nlock[20][0], sum(nlock[20]), nlock[50][0], sum(nlock[50])))
if bad:
    print("jeux incoherents (s1(50) < s1(20)) :", ", ".join(bad))
