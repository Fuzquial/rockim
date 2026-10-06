# -*- coding: utf-8 -*-
"""Figure de synthese de la calibration Red Bohus de rockim (FDEM 2D).

Trois panneaux, toutes les donnees lues dans les fichiers d'archive (aucune
valeur saisie a la main hormis les etiquettes) :
  (a) courbes experimentales sigma_1(eps_axial) de Dumoulin et al. (2024)
      aux quatre confinements (cibles) ;
  (b) enveloppe sigma_1 au pic contre sigma_3 : essais (moyenne +- ecart-type,
      sigma_1 = q + sigma_3), GBM de juillet (jeu calibre, graine 12345,
      sigma_3 REELLEMENT atteint), calibration d'aout (jeu CALT fenetre longue,
      jeu CAL1 fenetre tronquee), ajustement de Hoek-Brown des essais ;
      domaine de calage d'aout (sigma_3 <= 20 MPa vu par l'emulateur, 50 en
      controle) grise, 75 et 100 MPa = prediction pure ;
  (c) ecart relatif au pic experimental, par confinement.

Usage : python fig_calib_synthese.py [sortie.pdf]
"""
import csv
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

ROCKIM = "/home/user/rockim"
PHD = "/home/user/phd_geothermie/FDEM/rockim"
TARGETS = os.path.join(ROCKIM, "calibration_redbohus/targets/targets_redbohus.json")
CURVES = os.path.join(ROCKIM, "calibration_redbohus/targets/curves_redbohus.json")
POINTS = os.path.join(ROCKIM, "calibration_redbohus/points_results.csv")
FINAL_LOG = os.path.join(ROCKIM, "calibration_redbohus/tools/final_console.log")
JULY = os.path.join(PHD, "triax_envelope_2026-08-05.json")
PHI = os.path.join(PHD, "phijoint_sweep_2026-08-05.json")

out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "fig_calib_synthese.pdf")

plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix",
                     "font.size": 9, "axes.linewidth": 0.6,
                     "axes.unicode_minus": False})


def virgule(x, _pos):
    s = ("%g" % x).replace(".", ",").replace("-", "−")
    return s


FMT = FuncFormatter(virgule)

# ------------------------------------------------------------------ essais
T = json.load(open(TARGETS))
s3_exp = np.array([0.0] + T["triaxial"]["sigma3_MPa"])
s1_exp = np.array([T["UCS"]["mean_MPa"]]
                  + [q + s for q, s in zip(T["triaxial"]["q_peak_mean_MPa"],
                                           T["triaxial"]["sigma3_MPa"])])
sd_exp = np.array([T["UCS"]["std_MPa"]] + T["triaxial"]["q_peak_std_MPa"])
EXP = dict(zip(s3_exp, s1_exp))

# Hoek-Brown roche intacte (s = 1, a = 0,5), linearisation classique
# (Hoek & Brown 1980) : (sigma_1 - sigma_3)^2 = sigma_ci^2 + m_i sigma_ci sigma_3,
# regression lineaire ordinaire sur les 5 moyennes
def hb(s3, sci, mi):
    return s3 + sci * np.sqrt(mi * s3 / sci + 1.0)


_b, _a = np.polyfit(s3_exp, (s1_exp - s3_exp) ** 2, 1)
sci = np.sqrt(_a)
mi = _b / sci

# ------------------------------------------------- GBM de juillet (2026-07-31)
jul = [r for r in json.load(open(JULY)) if str(r.get("seed")) == "12345"]
jul_s3 = np.array([r["s3_ach_MPa"] if r.get("s3_ach_MPa") is not None
                   else r["s3_target_MPa"] for r in jul])
jul_s1 = np.array([r["peak_MPa"] for r in jul])
# point a sigma_3 = 50 du meme jeu (phi_joint = 13,4 deg) : balayage du 05/08
phi = [r for r in json.load(open(PHI)) if r["case"] == "j134_s50"]
jul_s3 = np.append(jul_s3, [r["s3_ach_MPa"] for r in phi])
jul_s1 = np.append(jul_s1, [r["peak_MPa"] for r in phi])
o = np.argsort(jul_s3)
jul_s3, jul_s1 = jul_s3[o], jul_s1[o]
jul_tgt = np.array([0, 5, 10, 20, 40, 50])  # consignes correspondantes

# ------------------------------------------------ calibration d'aout (CALT)
row = [r for r in csv.DictReader(open(POINTS)) if r["tag"] == "CALT"][0]
S3A = [20.0, 50.0, 75.0, 100.0]
calt = np.array([float(row["tx%d_peak_MPa" % s]) for s in S3A])

# jeu CAL1 (fenetre tronquee T = 4e-3 s) : journal de la passe finale
cal1 = {}
for ln in open(FINAL_LOG, encoding="utf-8", errors="replace"):
    if "CAL1" in ln and "pic=" in ln:
        tok = ln.split()
        test = tok[tok.index("CAL1") + 1]
        v = ln.split("pic=")[1].split()[0]
        try:
            cal1[test] = float(v)
        except ValueError:
            pass
cal1_s3 = np.array([0.0, 20.0, 50.0])
cal1_s1 = np.array([cal1["ucs"], cal1["tx20"], cal1["tx50"]])

# ------------------------------------------------------------------ figure
C_EXP, C_JUL, C_CALT, C_CAL1, C_HB = "0.15", "#1f6f8b", "#b03a2e", "#d68910", "0.55"
fig = plt.figure(figsize=(7.4, 6.3))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.38, wspace=0.28)

# (a) courbes experimentales
ax = fig.add_subplot(gs[0, :])
C = json.load(open(CURVES))
cols = {20: "#5d6d7e", 50: "#1f6f8b", 75: "#7d3c98", 100: "#b03a2e"}
seen = set()
for k, s in C["triaxial"].items():
    s3 = s["sigma3_MPa"]
    ax.plot(np.abs(s["eps_axial_pct"]), np.array(s["q_MPa"]) + s3,
            color=cols[s3], lw=0.9,
            label=(r"$\sigma_3 = %d$ MPa" % s3) if s3 not in seen else None)
    seen.add(s3)
first = True
for k, u in C["UC"].items():
    e = np.array(u["eps_local_pct"] if u.get("eps_local_pct") else u["eps_global_pct"], float)
    st = np.array(u["stress_local_MPa"] if u.get("stress_local_MPa") else u["stress_MPa"], float)
    ax.plot(np.abs(e), np.abs(st), color="0.2", lw=0.8, ls="--",
            label="compression simple (jauges décrochées avant le pic)" if first else None)
    first = False
ax.set_xlim(0, 1.1)
ax.set_ylim(0, 950)
ax.set_xlabel(r"déformation axiale $\varepsilon_a$ [%]")
ax.set_ylabel(r"contrainte axiale $\sigma_1$ [MPa]")
ax.set_title("(a) cibles : essais Red Bohus (Dumoulin et al., 2024)", fontsize=9.5)
ax.legend(frameon=False, fontsize=7.5, ncol=3, loc="upper left")

# (b) enveloppe
ax = fig.add_subplot(gs[1, 0])
ax.axvspan(-5, 20, color="0.88", zorder=0)
ax.axvspan(20, 50, color="0.95", zorder=0)
ax.text(-3, 1270, "calé", fontsize=7, color="0.3", va="top")
ax.text(22, 1270, "contrôle", fontsize=7, color="0.3", va="top")
ax.text(62, 1270, "prédiction pure", fontsize=7, color="0.3", va="top")
xs = np.linspace(0, 105, 200)
ax.plot(xs, hb(xs, sci, mi), color=C_HB, lw=1.0, ls=":",
        label=r"Hoek-Brown ajusté ($\sigma_{ci}$ = %s MPa ; $m_i$ = %s)"
        % (("%.0f" % sci), ("%.1f" % mi).replace(".", ",")))
ax.errorbar(s3_exp, s1_exp, yerr=sd_exp, fmt="s", color=C_EXP, ms=4.5,
            capsize=2, lw=0.8, label="essais (moyenne ± écart-type)")
ax.plot(jul_s3, jul_s1, "o-", color=C_JUL, ms=3.8, lw=1.0,
        label="GBM juillet (calé à $\\sigma_3$ = 0)")
ax.plot(S3A, calt, "D-", color=C_CALT, ms=4, lw=1.0,
        label="août CALT (fenêtre longue)")
ax.plot(cal1_s3, cal1_s1, "^--", color=C_CAL1, ms=4.5, lw=0.9,
        label="août CAL1 (fenêtre tronquée)")
ax.set_xlim(-5, 105)
ax.set_ylim(0, 1300)
ax.set_xlabel(r"confinement $\sigma_3$ [MPa]")
ax.set_ylabel(r"$\sigma_1$ au pic [MPa]")
ax.set_title("(b) enveloppe de rupture", fontsize=9.5)
axb = ax

# (c) ecarts relatifs
ax = fig.add_subplot(gs[1, 1])
ax.axvspan(-5, 20, color="0.88", zorder=0)
ax.axvspan(20, 50, color="0.95", zorder=0)
ax.axhline(0, color="0.3", lw=0.7)
rel_sd = 100 * sd_exp / s1_exp
ax.errorbar(s3_exp, 0 * s3_exp, yerr=rel_sd, fmt="none", ecolor="0.25",
            elinewidth=2.2, capsize=0, label="dispersion des essais (±1 écart-type)")
jul_err = [100 * (v - EXP[t]) / EXP[t] for t, v in zip(jul_tgt, jul_s1) if t in EXP]
jul_e3 = [s for t, s in zip(jul_tgt, jul_s3) if t in EXP]
ax.plot(jul_e3, jul_err, "o-", color=C_JUL, ms=3.8, lw=1.0, label="GBM juillet")
ax.plot(S3A, [100 * (v - EXP[s]) / EXP[s] for s, v in zip(S3A, calt)], "D-",
        color=C_CALT, ms=4, lw=1.0, label="août CALT")
ax.plot(cal1_s3, [100 * (v - EXP[s]) / EXP[s] for s, v in zip(cal1_s3, cal1_s1)],
        "^--", color=C_CAL1, ms=4.5, lw=0.9, label="août CAL1")
ax.set_xlim(-5, 105)
ax.set_ylim(-80, 55)
ax.set_xlabel(r"confinement $\sigma_3$ [MPa]")
ax.set_ylabel(r"écart au pic expérimental [%]")
ax.set_title("(c) écart relatif", fontsize=9.5)
hb_, lb_ = axb.get_legend_handles_labels()
hc_, lc_ = ax.get_legend_handles_labels()
H = hb_ + [hc_[lc_.index("dispersion des essais (±1 écart-type)")]]
L = lb_ + ["dispersion des essais (±1 écart-type)"]
fig.legend(H, L, frameon=False, fontsize=7.5, ncol=3, loc="upper center",
           bbox_to_anchor=(0.5, 0.035))

for a in fig.axes:
    a.xaxis.set_major_formatter(FMT)
    a.yaxis.set_major_formatter(FMT)
    a.tick_params(labelsize=8, width=0.6)
    a.grid(alpha=0.25, lw=0.4)

fig.savefig(out, bbox_inches="tight")
fig.savefig(out.replace(".pdf", ".png"), dpi=170, bbox_inches="tight")

print("Hoek-Brown : sigma_ci = %.1f MPa, m_i = %.2f" % (sci, mi))
print("%-6s %8s %9s %9s %9s" % ("s3", "exp s1", "juillet", "CALT", "CAL1"))
for s in s3_exp:
    j = [v for t, v in zip(jul_tgt, jul_s1) if t == s]
    print("%-6g %8.1f %9s %9s %9s" % (
        s, EXP[s], ("%.1f" % j[0]) if j else "-",
        ("%.1f" % calt[S3A.index(s)]) if s in S3A else "-",
        ("%.1f" % cal1_s1[list(cal1_s3).index(s)]) if s in cal1_s3 else "-"))
print("ecrit", out)
