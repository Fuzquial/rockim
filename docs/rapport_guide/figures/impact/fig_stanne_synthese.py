# -*- coding: utf-8 -*-
"""Planche de synthese de la replique St Anne (Yang et al. 2025), run out_stanne2025_rock137.

Donnees lues, aucune valeur inventee :
  results/data/stanne2025_rock137/history.csv      jauge, forces, vitesses, postes d'energie
  results/data/stanne2025_rock137/journal_solveur.log  masses des corps (fin du journal)
  results/fig/stanne300/evolution_cache.json        rayon de peau, extension, profondeur par trame
  output/pdf/stanne_<t>us/surface/mesures.json      rayon maximal des traces en surface
Theorie : impact 1D de barres (piston R 13,25 mm sur taillant R 15 mm, tools/make_impact_mesh.py l. 64-66).
Reperes Yang 2025 : lectures des fig. 9-10 (docs/COMPARAISON_yang2025_stanne_2026-09-14.md l. 31-57),
fin de charge 254 us (meme fichier l. 87).

Lancer depuis la racine du depot rockim :  python3 docs/rapport_guide/figures/impact/fig_stanne_synthese.py
"""
import csv, json, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
DATA = os.path.join(ROOT, "results", "data", "stanne2025_rock137")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stanne_synthese")

plt.rcParams.update({"font.size": 9.5,
                     "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.5,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.4,
                     "legend.frameon": False, "legend.fontsize": 8})
# --- Style commun du rapport : Computer Modern (CMU Serif, a defaut Latin Modern
# Roman 10), virgule decimale sur toutes les graduations, signe moins ASCII.
import glob as _glob
import matplotlib.ticker as _mticker
from matplotlib import font_manager as _fm
for _f in (_glob.glob("/usr/share/fonts/**/cmun*.[ot]tf", recursive=True)
           + _glob.glob("/usr/share/texmf/fonts/opentype/public/lm/lmroman10-*.otf")
           + _glob.glob("/usr/share/texlive/texmf-dist/fonts/opentype/public/lm/lmroman10-*.otf")):
    _fm.fontManager.addfont(_f)
_fm.fontManager.ttflist = [_e for _e in _fm.fontManager.ttflist
                           if _e.name != "Latin Modern Roman" or "lmroman10" in _e.fname]
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "Latin Modern Roman", "Computer Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
    "axes.formatter.use_locale": False,
    "pdf.fonttype": 3,
})
_sf_call = _mticker.ScalarFormatter.__call__


def _sf_virgule(self, x, pos=None):
    s = _sf_call(self, x, pos)
    return s.replace(".", "{,}") if "$" in s else s.replace(".", ",")


_mticker.ScalarFormatter.__call__ = _sf_virgule
# --- fin du style commun
VIRG = FuncFormatter(lambda x, p: ("%g" % x).replace(".", ","))
C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#555555"]

# ---- historique -------------------------------------------------------------------------
rows = list(csv.DictReader(open(os.path.join(DATA, "history.csv"))))
h = {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}
t = h["t"]; tu = t * 1e6
log = open(os.path.join(DATA, "journal_solveur.log"), errors="replace").read()
m = {k: float(v) for k, v in re.findall(r"corps '(\w+)':.*?masse = ([\d.eE+-]+) kg", log)}

# (a) jauge et theorie 1D
rho, E, v0 = 7850.0, 200e9, 10.66
c0 = np.sqrt(E / rho)
Ap, Ab = np.pi * 0.01325 ** 2, np.pi * 0.015 ** 2
s_eq = rho * c0 * v0 / 2 / 1e6                      # aires egales
s_ar = rho * c0 * v0 * Ap / (Ap + Ab) / 1e6          # aires du maillage
gauge = -h["szz_bit"] / 1e6

# (b) force mesuree et estimateurs de quantite de mouvement
def sm(y, n):
    return np.convolve(y, np.ones(n) / n, mode="same")
n = max(1, int(round(2e-6 / np.median(np.diff(t)))))
P_rig = m["piston"] * h["vz_piston"] + (m["bit"] + m["insert"] + m["circlip"]) * h["vz_bit"]
P_cor = m["piston"] * h["vz_piston"] + (m["bit"] + m["circlip"]) * h["vz_bit"] + m["insert"] * h["vz_insert"]
F_rig = sm(sm(np.gradient(P_rig, t), n), n) / 1e3
F_cor = sm(sm(np.gradient(P_cor, t), n), n) / 1e3
F_mes = h["Fc_rock_insert_z"] / 1e3

# (c) cratere et fissures par trame
ev = json.load(open(os.path.join(ROOT, "results", "fig", "stanne300", "evolution_cache.json")))
fr = sorted(ev, key=int)
fr_t = {int(r["frame"]): float(r["t"]) * 1e6 for r in csv.DictReader(open(os.path.join(DATA, "frames.csv")))}
tf = np.array([fr_t[int(k)] for k in fr])
rC = np.array([ev[k]["rCrater"] for k in fr]); rM = np.array([ev[k]["rMax"] for k in fr])
dep = np.array([ev[k]["depth"] for k in fr])
surf_t, surf_r = [], []
for k in sorted(os.listdir(os.path.join(ROOT, "output", "pdf"))):
    p = os.path.join(ROOT, "output", "pdf", k, "surface", "mesures.json")
    if k.startswith("stanne_") and os.path.isfile(p):
        d = json.load(open(p)); surf_t.append(d["time_us"]); surf_r.append(d["surface"]["max_radius_mm"])
for k in ("stanne_surface_165us", "stanne_surface_180us"):
    d = json.load(open(os.path.join(ROOT, "output", "pdf", k, "mesures.json")))
    surf_t.append(d["time_us"]); surf_r.append(d["surface"]["max_radius_mm"])
o = np.argsort(surf_t); surf_t = np.array(surf_t)[o]; surf_r = np.array(surf_r)[o]

# (d) energie
KE = 0.5 * (m["piston"] * h["vz_piston"] ** 2 + (m["bit"] + m["circlip"]) * h["vz_bit"] ** 2
            + m["insert"] * h["vz_insert"] ** 2)
T_END = 284.3

fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.6))
a = ax[0, 0]
a.plot(tu, gauge, color=C[0], label="jauge à mi-taillant (rockim)")
a.axhline(s_eq, color=C[5], ls="--", lw=1, label=(r"$\rho c v_0/2$ = %.1f MPa (aires égales)" % s_eq).replace(".", ","))
a.axhline(s_ar, color=C[1], ls=":", lw=1.4,
          label=(r"$\rho c v_0\,A_p/(A_p+A_b)$ = %.1f MPa" % s_ar).replace(".", ","))
a.set_xlim(0, 120); a.set_xlabel(r"temps [$\mu$s]"); a.set_ylabel(r"$-\sigma_{zz}$ [MPa]")
a.set_title("(a) onde incidente et théorie 1D", loc="left", fontsize=9.5)
a.legend(loc="lower center")

a = ax[0, 1]
e = tu < 297
a.plot(tu[e], F_rig[e], color=C[5], lw=0.9, ls="--", label=r"estimateur $\dot P$, train rigide")
a.plot(tu[e], F_cor[e], color=C[1], lw=1.0, label=r"estimateur $\dot P$, insert à sa vitesse")
a.plot(tu, F_mes, color=C[0], lw=1.6, label="force de contact mesurée")
a.axvline(T_END, color="#999999", lw=0.8)
a.set_xlabel(r"temps [$\mu$s]"); a.set_ylabel("force roche-insert [kN]")
a.set_title("(b) force mesurée et estimateurs", loc="left", fontsize=9.5)
a.legend(loc="lower right")

a = ax[1, 0]
a.plot(tf, rC, "o-", color=C[0], ms=3.5, label="rayon de peau du cratère")
a.plot(tf, rM, "s--", color=C[1], ms=3.5, label="extension des facettes rompues")
a.plot(surf_t, surf_r, "^:", color=C[2], ms=3.5, label="rayon max. des traces en surface")
a.plot(tf, dep, "d-.", color=C[3], ms=3.5, label="profondeur du réseau")
a.axhline(13.5, color=C[0], lw=0.8, ls=(0, (1, 2)))
a.text(118, 14.2, "cratère Yang 13,5 mm", fontsize=7.5, color="#333333")
a.axhline(32, color=C[1], lw=0.8, ls=(0, (1, 2)))
a.text(118, 32.6, "radiale Yang 32 mm", fontsize=7.5, color="#333333")
a.axvline(T_END, color="#999999", lw=0.8)
a.set_xlabel(r"temps [$\mu$s]"); a.set_ylabel("longueur [mm]")
a.set_title("(c) cratère et réseau de fissures", loc="left", fontsize=9.5)
a.set_ylim(0, 72)
a.legend(loc="upper left", ncol=1, fontsize=7.0)

a = ax[1, 1]
a.plot(tu, KE, color=C[0], label="cinétique de translation du train")
a.plot(tu, -h["eEl"], color=C[1], ls="--", label="éléments (prélevé)")
a.plot(tu, -h["eJnt"], color=C[2], ls="-.", label="joints cohésifs")
a.plot(tu, -h["eGc"], color=C[3], label="contact général (total)")
a.plot(tu, -h["eFric"], color=C[4], ls=":", lw=1.8, label="dont frottement")
a.axvline(T_END, color="#999999", lw=0.8)
a.set_xlabel(r"temps [$\mu$s]"); a.set_ylabel("énergie [J]")
a.set_title("(d) postes cumulés du bilan", loc="left", fontsize=9.5)
a.legend(loc="upper right", fontsize=7.2)

for a in ax.flat:
    a.xaxis.set_major_formatter(VIRG); a.yaxis.set_major_formatter(VIRG)
fig.tight_layout()
fig.savefig(OUT + ".pdf"); fig.savefig(OUT + ".png", dpi=200)

# ---- chiffres imprimes (pour le dossier) --------------------------------------------------
pl = (t > 40e-6) & (t < 80e-6)
i1 = np.argmax(gauge * (t < 45e-6))
print("masses", m)
print("1D : aires egales %.1f MPa, aires du maillage %.1f MPa ; jauge pic %.1f MPa a %.1f us, plateau 40-80 us moyenne %.1f mediane %.1f"
      % (s_eq, s_ar, gauge[i1], tu[i1], gauge[pl].mean(), np.median(gauge[pl])))
for T in (134e-6, 300.1e-6):
    k = t <= T
    print("fenetre 0-%.0f us : F mes max %.1f kN ; rigide %.1f (%+.0f %%) ; corrige %.1f ; impulsion mes %.3f, rigide %.3f, corrige %.3f N s"
          % (T * 1e6, F_mes[k].max(), F_rig[k].max(), 100 * (F_rig[k].max() / F_mes[k].max() - 1), F_cor[k].max(),
             np.trapezoid(h["Fc_rock_insert_z"][k], t[k]), P_rig[k][-1] - P_rig[0], P_cor[k][-1] - P_cor[0]))
k = (t > 40e-6) & (t < 295e-6)
print("ecart rms 40-295 us : rigide %.1f kN, corrige %.1f kN" % (np.sqrt(np.mean((F_rig - F_mes)[k] ** 2)),
                                                              np.sqrt(np.mean((F_cor - F_mes)[k] ** 2))))
print("fin : eEl %.2f eJnt %.2f eGc %.2f eFric %.2f ; KE corps rigides fin %.2f J" %
      (-h["eEl"][-1], -h["eJnt"][-1], -h["eGc"][-1], -h["eFric"][-1], KE[-1]))
print("surface", list(zip(surf_t.round(1), surf_r.round(2))))
