# -*- coding: utf-8 -*-
"""Enveloppe de rupture : ce que la troncature masquait.

Trois informations sur la meme figure :
 (a) les courbes sigma1(eps) aux quatre confinements, fenetre allongee,
     avec le trait vertical qui marque OU s'arretait l'ancienne fenetre ;
 (b) l'enveloppe sigma1(sigma3) simulee contre l'experimentale ;
 (c) la pente locale dq/dsigma3, qui montre la raideur excessive.
"""
import csv, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
d = json.load(open(os.path.join(BASE, "targets", "curves_redbohus.json")))
S3 = [20, 50, 75, 100]
COL = {20: "C0", 50: "C1", 75: "C2", 100: "C3"}
EPS_OLD = 0.85          # % : deformation max sous l'ancienne fenetre T = 4e-3


def hist(run):
    p = os.path.join(BASE, "runs", run, "history.csv")
    if not os.path.exists(p):
        return None
    rows = [r for r in csv.DictReader(open(p))
            if all(v not in (None, "") for v in r.values())]
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


fig = plt.figure(figsize=(15.5, 5.6))
gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1])

# ---------------------------------------------------------- (a) courbes
a = fig.add_subplot(gs[0])
for s3 in S3:
    first = True
    for k, s in d["triaxial"].items():
        if s["sigma3_MPa"] != s3:
            continue
        a.plot(s["eps_axial_pct"], np.array(s["q_MPa"]) + s3, color=COL[s3],
               lw=1.0, ls="--", alpha=0.55,
               label="essais $\\sigma_3$=%d" % s3 if first else None)
        first = False
    h = hist("CALT_tx%d_s4211" % s3)
    if h is not None:
        a.plot(np.abs(h["epsGauge"]) * 100, np.abs(h["sigma"]) * 1e-6,
               color=COL[s3], lw=2, label="rockim $\\sigma_3$=%d" % s3)
a.axvline(EPS_OLD, color="k", lw=1.6, ls=":")
a.text(EPS_OLD + 0.02, 60, "ancienne fenêtre\n(T = 4·10⁻³ s)", fontsize=8,
       rotation=90, va="bottom")
a.set_xlim(0, 1.5); a.set_ylim(0, 1300)
a.set_xlabel("déformation axiale (%)"); a.set_ylabel(r"$\sigma_1$ (MPa)")
a.set_title("Courbes aux quatre confinements — fenêtre allongée")
a.legend(fontsize=7, ncol=2); a.grid(alpha=0.3)

# ---------------------------------------------------------- (b) enveloppe
a = fig.add_subplot(gs[1])
row = [r for r in csv.DictReader(open(os.path.join(BASE, "points_results.csv")))
       if r["tag"] == "CALT"][0]
sim = [float(row["tx%d_peak_MPa" % s]) for s in S3]
exp = {20: 424.8, 50: 649.0, 75: 779.0, 100: 899.3}
a.plot([0] + S3, [126.6] + [exp[s] for s in S3], "s--", color="0.35", lw=2,
       ms=7, label="expérimental (Dumoulin 2024)")
a.plot(S3, sim, "o-", color="C3", lw=2, ms=7, label="rockim (fenêtre longue)")
a.plot([20, 50], [390.9, 625.9], "^:", color="C7", lw=1.5, ms=7,
       label="rockim (fenêtre tronquée)")
for s, v, e in zip(S3, sim, [exp[s] for s in S3]):
    a.annotate("%+.0f %%" % ((v - e) / e * 100), (s, v),
               textcoords="offset points", xytext=(6, 4), fontsize=8,
               color="C3")
a.set_xlabel(r"$\sigma_3$ (MPa)"); a.set_ylabel(r"$\sigma_1$ au pic (MPa)")
a.set_title("Enveloppe de rupture"); a.legend(fontsize=8); a.grid(alpha=0.3)

# ---------------------------------------------------------- (c) pentes
a = fig.add_subplot(gs[2])
xs = [(0 + 20) / 2, (20 + 50) / 2, (50 + 75) / 2, (75 + 100) / 2]
pe = [(exp[20] - 126.6) / 20 - 1, (exp[50] - exp[20]) / 30 - 1,
      (exp[75] - exp[50]) / 25 - 1, (exp[100] - exp[75]) / 25 - 1]
ps = [np.nan, (sim[1] - sim[0]) / 30 - 1, (sim[2] - sim[1]) / 25 - 1,
      (sim[3] - sim[2]) / 25 - 1]
a.plot(xs, pe, "s--", color="0.35", lw=2, ms=7, label="expérimental")
a.plot(xs, ps, "o-", color="C3", lw=2, ms=7, label="rockim")
a.set_xlabel(r"$\sigma_3$ (MPa)"); a.set_ylabel(r"pente locale d$q$/d$\sigma_3$")
a.set_title("Raideur de l'enveloppe"); a.legend(fontsize=8); a.grid(alpha=0.3)

fig.suptitle("Ce que la troncature masquait : l'enveloppe simulée n'est pas trop "
             "plate, elle est TROP RAIDE (et croise l'expérimentale vers "
             "$\\sigma_3$ ≈ 35 MPa)", fontsize=12)
fig.tight_layout()
out = os.path.join(BASE, "figures", "enveloppe.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
print("\n%-6s %10s %12s %9s" % ("s3", "rockim", "experimental", "ecart"))
for s, v in zip(S3, sim):
    print("%-6d %10.1f %12.1f %8.1f %%" % (s, v, exp[s], (v - exp[s]) / exp[s] * 100))
