# -*- coding: utf-8 -*-
"""fig_synthese.py — figure de synthese des bancs d'insertion (Computer Modern, virgule decimale).
usage : python fig_synthese.py <racine_bancs> <sortie.pdf>"""
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from matplotlib import font_manager
for f in ("lmroman10-regular.otf", "lmroman10-italic.otf", "lmroman10-bold.otf"):
    p_ = "/usr/share/texmf/fonts/opentype/public/lm/" + f
    if os.path.exists(p_): font_manager.fontManager.addfont(p_)
# Computer Modern : Latin Modern Roman (dessin de Knuth) pour le texte, fontset cm pour les maths
plt.rcParams.update({"text.usetex": False, "font.family": "Latin Modern Roman", "mathtext.fontset": "cm",
                     "font.size": 8.5, "axes.unicode_minus": True,
                     "axes.linewidth": 0.6, "lines.linewidth": 1.1, "legend.frameon": False})
R, OUT = sys.argv[1], sys.argv[2]
virg = FuncFormatter(lambda v, p: ("%g" % v).replace(".", ",").replace("-", "\u2212"))
C_AD, C_IN = "#1f5fa8", "#c4482b"

def rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))

fig, ax = plt.subplots(2, 2, figsize=(6.6, 5.6))
# (a) banc 1 : sigma_n(delta_n) mode I et tau(delta_s) mode II, b1s, normes par delta_c de la loi
a = ax[0, 0]
d = np.load(os.path.join(R, "b1s", "b1s_courbes.npz"))
ft, c, GI, GII, I = 10e6, 25e6, 70.0, 700.0, 0.386307
def fyan(D, A=0.63, B=1.8, C=6.0):
    return (1 - (A + B - 1) / (A + B) * np.exp(D * (A + C * B) / ((A + B) * (1 - A - B)))) * (A * (1 - D) + B * (1 - D) ** C)
dcI = {"lin": 2 * GI / ft, "yan": GI / (ft * I)}; dcII = {"lin": 2 * GII / c, "yan": GII / (c * I)}
for soft, ls in (("lin", "-"), ("yan", "--")):
    for ins, col in (("adapt", C_AD), ("intri", C_IN)):
        nm = "b1s_I_%s_p100_%s" % (ins, soft)
        a.plot(d[nm + "__dn"] / dcI[soft], d[nm + "__sn"] / ft, ls, color=col, lw=0.9)
        nm = "b1s_II_%s_p100_%s" % (ins, soft)
        a.plot(np.abs(d[nm + "__ds"]) / dcII[soft], np.abs(d[nm + "__tau"]) / c, ls, color=col, lw=0.9, alpha=0.45)
x = np.linspace(0, 1, 200)
a.plot(x, 1 - x, ":", color="k", lw=0.9)
a.plot(x, fyan(x), ":", color="k", lw=0.9)
a.set_xlim(0, 1.3); a.set_ylim(-0.05, 1.12)
a.set_xlabel(r"$\delta/\delta_c$ (mode I : $\delta_n$ ; mode II, pâle : $\delta_s$)")
a.set_ylabel(r"$\sigma_n/f_t$ (mode I) ; $\tau/c$ (mode II)")
a.set_title(r"(a) banc 1 : joint isolé, 100 $E/h$", fontsize=8.5)
a.plot([], [], "-", color=C_AD, label="adaptatif"); a.plot([], [], "-", color=C_IN, label="intrinsèque")
a.plot([], [], "-", color="0.4", label="adoucissement linéaire"); a.plot([], [], "--", color="0.4", label="Yan (éq. 11)")
a.plot([], [], ":", color="k", label="loi exacte")
a.legend(fontsize=6.8, loc="upper right")
# (b) banc 2 : pic et energie contre h/l_cz
b = ax[0, 1]
r2 = rows(os.path.join(R, "b2", "b2_resultats.csv")) + rows(os.path.join(R, "b2L", "b2L_resultats.csv"))
hh = {"h10": 0.5, "h5": 0.25, "h2.5": 0.125, "h1.25": 0.0625}
# adaptatif : serie b2L (T = 3 ms, barres separees) ; intrinseque : serie b2 (T = 1,25 ms)
for ins, col, serie in (("adapt", C_AD, "b2L"), ("intri", C_IN, "b2")):
    for corr, mk in (("corr", "o"), ("brut", "s")):
        pts = sorted((hh[r["name"].split("_")[1]], float(r["sig_pk_ft"]), float(r["W_GW"])) for r in r2
                     if r["name"].split("_")[0] == serie and r["name"].split("_")[2] == ins and r["name"].split("_")[3] == corr)
        if not pts: continue
        p = np.array(pts)
        off = 1.03 if corr == "brut" else 1.0
        b.plot(p[:, 0] * off, p[:, 1], "-" + mk, color=col, ms=3.5, mfc=col if corr == "corr" else "white", lw=0.9)
        b.plot(p[:, 0] * off, p[:, 2], "--" + mk, color=col, ms=3.5, mfc=col if corr == "corr" else "white", lw=0.9)
b.set_xscale("log", base=2); b.set_xticks([0.0625, 0.125, 0.25, 0.5])
b.set_xticklabels([r"1/16", r"1/8", r"1/4", r"1/2"])
b.yaxis.set_major_formatter(virg)
b.set_xlabel(r"$h/\ell_{cz}$"); b.set_ylabel(r"$\sigma_\mathrm{pic}/f_t$ (trait plein) ; $W/(G_I\,b)$ (tirets)")
b.set_title(r"(b) banc 2 : barre en traction", fontsize=8.5)
b.plot([], [], "-o", color=C_AD, ms=3, label="adaptatif, contact corrigé")
b.plot([], [], "-s", color=C_AD, mfc="none", ms=3, label="adaptatif, sans correctif")
b.plot([], [], "-o", color=C_IN, ms=3, label="intrinsèque, contact corrigé")
b.plot([], [], "-s", color=C_IN, mfc="none", ms=3, label="intrinsèque, sans correctif")
b.legend(fontsize=6.0, loc="lower center", ncol=2, columnspacing=0.8, handlelength=1.6)
b.set_ylim(0.62, 2.2)
# (c) banc 3 : fissures de la plaque entaillee
cax = ax[1, 0]
d3 = np.load(os.path.join(R, "b3", "b3_fissures.npz"))
for k, (nm, col, dx) in enumerate((("b3_h1.25_adapt", C_AD, 0.0), ("b3_h1.25_intri", C_IN, 0.045))):
    if nm + "__seg" not in d3: nm = nm.replace("h1.25", "h2.5")
    if nm + "__seg" not in d3: continue
    s = d3[nm + "__seg"] * 1e3; m = d3[nm + "__main"] * 1e3
    for q in s:
        cax.plot([q[0] + dx * 1e3, q[2] + dx * 1e3], [q[1], q[3]], "-", color="0.65", lw=0.6)
    for q in m:
        cax.plot([q[0] + dx * 1e3, q[2] + dx * 1e3], [q[1], q[3]], "-", color=col, lw=1.1)
    x0 = dx * 1e3
    cax.plot([x0, x0 + 40, x0 + 40, x0, x0, x0 + 10, x0 + 10, x0, x0],
             [0, 0, 80, 80, 40.5, 40.5, 39.5, 39.5, 0], "-", color="k", lw=0.5)
cax.set_aspect("equal"); cax.set_xlim(-2, 87); cax.set_ylim(22, 58)
cax.set_xticks([]); cax.set_yticks([])
cax.text(20, 23, "adaptatif", ha="center", va="top", fontsize=7.5, color=C_AD)
cax.text(65, 23, "intrinsèque", ha="center", va="top", fontsize=7.5, color=C_IN)
cax.set_title(r"(c) banc 3 : SENT, $h$ = 1,25 mm (zoom)", fontsize=8.5)
# (d) banc 4 : sensibilite
dax = ax[1, 1]
r4 = rows(os.path.join(R, "b4", "b4_resultats.csv"))
var = ["base", "elliptic", "volume", "max", "hold5", "tip1.6"]
lab = ["base", "elliptic", "volume", "max", "hold 5", "tip 1,6"]
for j, (hn, col) in enumerate((("h5", "0.35"), ("h2.5", "0.7"))):
    v = {r["name"].split("_", 2)[2]: r for r in r4 if r["name"].split("_")[1] == hn}
    if not v: continue
    xb = np.arange(len(var)) + (j - 0.5) * 0.36
    dax.bar(xb, [float(v[k]["n_ins_or_dmg"]) if k in v else 0 for k in var], 0.34, color=col,
            label=r"joints insérés, $h$ = %s mm" % hn[1:].replace(".", ","))
    dax.plot(xb, [float(v[k]["n_broken"]) if k in v else 0 for k in var], "k_", ms=8)
d2 = dax.twinx()
for j, (hn, mk) in enumerate((("h5", "o"), ("h2.5", "s"))):
    v = {r["name"].split("_", 2)[2]: r for r in r4 if r["name"].split("_")[1] == hn}
    if not v: continue
    d2.plot(np.arange(len(var)) + (j - 0.5) * 0.36, [float(v[k]["sig_pk_ft"]) if k in v else np.nan for k in var],
            mk, color=C_AD, ms=3.5)
d2.set_ylabel(r"$\sigma_\mathrm{pic}/f_t$ (points bleus)", color=C_AD); d2.yaxis.set_major_formatter(virg)
dax.set_xticks(range(len(var))); dax.set_xticklabels(lab, rotation=30, fontsize=7)
dax.set_ylabel("joints insérés (barres), rompus (traits)")
dax.set_title(r"(d) banc 4 : critère d'insertion (adaptatif)", fontsize=8.5)
dax.legend(fontsize=6.3, loc="center left", bbox_to_anchor=(0.0, 0.62))
for x in (a, dax):
    x.yaxis.set_major_formatter(virg)
a.xaxis.set_major_formatter(virg)
fig.tight_layout()
fig.savefig(OUT)
print("ecrit", OUT)
