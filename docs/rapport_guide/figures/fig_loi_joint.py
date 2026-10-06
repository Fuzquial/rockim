"""Figure de principe de la loi de joint de rockim (branche parabolique de Guo,
adoucissement en z de Munjiza) et des facteurs dynamiques de Yang et al.

Aucune donnee de calcul : la figure trace les equations du rapport-guide.
Usage : python fig_loi_joint.py  (ecrit fig_loi_joint.pdf a cote du script)
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.size": 9,
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.2,
})
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

A, B, C = 0.63, 1.8, 6.0


def z(D):
    """Courbe d'adoucissement de Munjiza (eq. 2.32 de Guo)."""
    e = np.exp(D * (A + C * B) / ((A + B) * (1.0 - A - B)))
    return (1.0 - (A + B - 1.0) / (A + B) * e) * (A * (1.0 - D) + B * (1.0 - D) ** C)


fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.3))

# (a) traction normale contre ouverture, en unites reduites
dnp, dnc = 1.0, 6.0
x1 = np.linspace(-0.4, 0.0, 20)
x2 = np.linspace(0.0, dnp, 50)
x3 = np.linspace(dnp, dnc, 200)
r = x2 / dnp
ax[0].plot(x1, 2.0 * x1 / dnp, color="k")
ax[0].plot(x2, 2.0 * r - r ** 2, color="k")
ax[0].plot(x3, z((x3 - dnp) / (dnc - dnp)), color="k")
ax[0].axhline(0, color="0.6", lw=0.5)
ax[0].axvline(0, color="0.6", lw=0.5)
ax[0].set_xlabel(r"ouverture $\delta_n/\delta_{np}$")
ax[0].set_ylabel(r"$\sigma/f_t$")
ax[0].set_title("(a) branche normale", fontsize=9)
ax[0].set_ylim(-0.9, 1.3)
ax[0].text(-0.3, 1.1, "élastique", fontsize=8)
ax[0].text(2.3, 0.55, "adoucissement", fontsize=8)

# (b) courbe en z
D = np.linspace(0, 1, 200)
ax[1].plot(D, z(D), color="k")
ax[1].fill_between(D, z(D), color="0.85")
ax[1].set_xlabel(r"endommagement $D$")
ax[1].set_ylabel(r"$z(D)$")
ax[1].set_title("(b) courbe en z", fontsize=9)
ax[1].text(0.35, 0.5, r"aire $=0{,}386$", fontsize=8)

# (c) facteurs dynamiques
e = np.logspace(-6, 4, 300)
dif_c = np.minimum(1.84, 0.77 + 0.56 * e ** 0.07)
for at, ls in ((0.07, "--"), (0.17, "-")):
    dif_t = np.minimum(1.85, 0.95 + 0.41 * e ** at)
    dif_t[e > 1e2] = 1.85
    ax[2].semilogx(e, dif_t, ls, color="k", label=r"$DIF_t$, $a_t=%s$" % str(at).replace(".", "{,}"))
ax[2].semilogx(e, dif_c, color="0.5", label=r"$DIF_c$")
ax[2].set_xlabel(r"$\dot\varepsilon$ (s$^{-1}$)")
ax[2].set_ylabel("facteur")
ax[2].set_title("(c) effet de vitesse", fontsize=9)
ax[2].legend(fontsize=7, frameon=False, loc="upper left")

fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig_loi_joint.pdf")
fig.savefig(out)
print(out)
