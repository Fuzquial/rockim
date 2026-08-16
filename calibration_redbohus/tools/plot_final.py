# -*- coding: utf-8 -*-
"""FIGURE DE RESULTAT de la calibration Red Bohus : le jeu calibre confronte
aux essais, avec la deformation mesuree par l'EXTENSOMETRE (epsGauge), et la
prediction pure a sigma_3 = 50 MPa (confinement jamais utilise au calage)."""
import csv, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
TAG = "CAL1"
d = json.load(open(os.path.join(BASE, "targets", "curves_redbohus.json")))


def hist(run, ecol="epsGauge"):
    p = os.path.join(BASE, "runs", run, "history.csv")
    if not os.path.exists(p):
        return None
    rows = [r for r in csv.DictReader(open(p))
            if all(v not in (None, "") for v in r.values())]
    o = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
    return o


fig = plt.figure(figsize=(15.5, 8.4))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1])

# ---------------------------------------------------- UCS
a = fig.add_subplot(gs[0, 0])
for i, (k, s) in enumerate(d["UC"].items()):
    a.plot(s["eps_local_pct"], s["stress_local_MPa"], color="0.55", lw=1.2,
           ls="--", label="essais (jauges locales)" if i == 0 else None)
h = hist("%s_ucs_s4211" % TAG)
if h is not None:
    a.plot(np.abs(h["epsGauge"]) * 100, np.abs(h["sigma"]) * 1e-6, "C0", lw=2,
           label="rockim calibré (%.0f MPa)" % (np.abs(h["sigma"]).max() / 1e6))
a.axhline(126.6, color="C2", lw=1.2, ls=":", label="moyenne exp. 126,6 MPa")
a.set_xlim(0, 0.45); a.set_ylim(0, 180)
a.set_xlabel("déformation axiale (%)"); a.set_ylabel(r"$\sigma_1$ (MPa)")
a.set_title("Compression simple"); a.legend(fontsize=8); a.grid(alpha=0.3)

# ---------------------------------------------------- brésilien
a = fig.add_subplot(gs[0, 1])
h = hist("%s_bts_s4211" % TAG)
if h is not None:
    a.plot((h["drive"] - h["drive"][0]) * 1e3, h["sigmaT"] * 1e-6, "C0", lw=2,
           label="rockim calibré (%.1f MPa)" % (h["sigmaT"].max() / 1e6))
a.axhline(10.27, color="C2", lw=1.5, label="exp. 10,27 MPa (4 essais)")
a.axhspan(10.27 - 0.98, 10.27 + 0.98, color="C2", alpha=0.15)
a.set_xlabel("rapprochement des plateaux (mm)")
a.set_ylabel(r"$\sigma_t$ (MPa)")
a.set_title("Brésilien"); a.legend(fontsize=8); a.grid(alpha=0.3)

# ---------------------------------------------------- triaxiaux
a = fig.add_subplot(gs[0, 2])
for s3, col in ((20, "C1"), (50, "C3")):
    first = True
    for k, s in d["triaxial"].items():
        if s["sigma3_MPa"] != s3:
            continue
        a.plot(s["eps_axial_pct"], np.array(s["q_MPa"]) + s3, color="0.55",
               lw=1.1, ls="--", label="essais réels" if (first and s3 == 20) else None)
        first = False
    h = hist("%s_tx%d_s4211" % (TAG, s3))
    if h is not None:
        lab = "rockim $\\sigma_3$=%d (%.0f MPa)%s" % (
            s3, np.abs(h["sigma"]).max() / 1e6,
            "  ← PRÉDICTION" if s3 == 50 else "")
        a.plot(np.abs(h["epsGauge"]) * 100, np.abs(h["sigma"]) * 1e-6, col,
               lw=2, label=lab)
a.set_xlim(0, 1.3); a.set_ylim(0, 720)
a.set_xlabel("déformation axiale (%)"); a.set_ylabel(r"$\sigma_1$ (MPa)")
a.set_title(r"Triaxiaux $\sigma_3$ = 20 et 50 MPa"); a.legend(fontsize=8)
a.grid(alpha=0.3)

# ---------------------------------------------------- ecarts
a = fig.add_subplot(gs[1, :2])
row = [r for r in csv.DictReader(open(os.path.join(BASE, "points_results.csv")))
       if r["tag"] == TAG][0]
obs = [("UCS", float(row["ucs_peak_MPa"]), 126.6, 21.4),
       (r"$\sigma_t$ brésilien", float(row["bts_sigma_t_MPa"]), 10.27, 0.98),
       (r"$\sigma_1$ ($\sigma_3$=20)", float(row["tx20_peak_MPa"]), 424.8, 2.8),
       (r"$\sigma_1$ ($\sigma_3$=50) — prédiction", float(row["tx50_peak_MPa"]),
        649.0, 2.6)]
y = np.arange(len(obs))
err = [(s - t) / t * 100 for _, s, t, _ in obs]
cols = ["C0" if abs(e) < 15 else "C3" for e in err]
a.barh(y, err, color=cols, height=0.55)
for i, ((n, s, t, sd), e) in enumerate(zip(obs, err)):
    a.text(e + (1.5 if e > 0 else -1.5), i, "%.1f vs %.1f  (%+.1f %%)"
           % (s, t, e), va="center", ha="left" if e > 0 else "right", fontsize=9)
a.axvline(0, color="k", lw=1)
a.axvspan(-15, 15, color="C2", alpha=0.12, label="±15 % (dispersion exp. de l'UCS)")
a.set_yticks(y); a.set_yticklabels([o[0] for o in obs])
a.set_xlim(-32, 32); a.set_xlabel("écart à l'expérimental (%)")
a.set_title("Bilan du jeu calibré — runs RÉELS (pas l'émulateur)")
a.legend(fontsize=8, loc="lower right"); a.grid(alpha=0.3, axis="x")

# ---------------------------------------------------- parametres
a = fig.add_subplot(gs[1, 2]); a.axis("off")
r = json.load(open(os.path.join(BASE, "calibration_result.json")))
txt = ["JEU CALIBRÉ  (postérieur, CI95)", ""]
lbl = {"ft": "ft joint (MPa)", "cohesion": "cohésion (MPa)",
       "frictionDeg": "frottement (°)", "Gf": "Gf mode I (J/m²)",
       "gfShearFactor": "GII/GI", "crushCap": "crushCap (MPa)"}
for k, v in r["posterior"].items():
    sc = 1e-6 if k in ("ft", "cohesion", "crushCap") else 1.0
    txt.append("%-18s %7.1f   [%.1f – %.1f]"
               % (lbl[k], v[0] * sc, v[1] * sc, v[2] * sc))
txt += ["", "GBM Voronoï, grains 5 mm", "insertion adaptative + potentiel",
        "E = 77,66 GPa, ν = 0,29"]
a.text(0.02, 0.95, "\n".join(txt), va="top", family="monospace", fontsize=9,
       transform=a.transAxes)

fig.suptitle("Calibration Red Bohus — résultat final (Dumoulin et al. 2024 en "
             "référence)", fontsize=13)
fig.tight_layout()
out = os.path.join(BASE, "figures", "resultat_calibration.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
for (n, s, t, sd), e in zip(obs, err):
    print("  %-28s %8.1f   cible %7.1f   %+6.1f %%" % (n, s, t, e))
