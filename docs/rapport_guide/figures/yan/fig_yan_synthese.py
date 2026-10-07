# -*- coding: utf-8 -*-
"""Synthese rockim / Yan, Zheng & Wang (IJRMMS 169, 2023, 105439).

Aucune donnee brute de la campagne 2026-08-10/11 n'est disponible dans les
depots (history.csv et VTU restes sur le poste Windows). Cette figure ne trace
donc QUE des valeurs chiffrees consignees dans les fiches, avec leur source :

  FICHE = phd_geothermie/FDEM/rockim/FICHE_rockim.md
  CHG   = rockim/CHANGES_YAN.md
  B6    = rockim/vv/B_articles/README.md
  YAN   = texte de l'article (pdftotext de bibliographie/yan2023.pdf)

(a) enveloppe de Coulomb : critere d'entree (c = 16,4 MPa, phi = 23 deg,
    Table 1 de Yan), pics triaxiaux 3D rockim (FICHE:1067-1069), droite 2D
    rockim ajustee phi = 22,8 deg passant par l'UCS adaptatif 51,1 MPa
    (FICHE:1004-1007).
(b) ecart relatif rockim / reference : valeur publiee par Yan, ou solution
    analytique (bande pesante eq. 20, SHPB -V0/c) quand elle existe.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.6})
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
virg = FuncFormatter(lambda v, p: ("%g" % v).replace(".", ","))

# ------------------------------------------------------------------ (a) Coulomb
c_in, phi_in = 16.4, 23.0                        # Table 1 de Yan
def coulomb(s3, c, phi):
    p = np.radians(phi)
    return 2 * c * np.cos(p) / (1 - np.sin(p)) + s3 * (1 + np.sin(p)) / (1 - np.sin(p))

s3 = np.linspace(0, 42, 50)
s3_3d = np.array([0, 10, 20, 40.])               # FICHE:1067-1068
s1_3d = np.array([56.2, 80.8, 105.5, 154.5])
phi2d, ucs2d = 22.8, 51.1                        # FICHE:1004, 1006
k2d = (1 + np.sin(np.radians(phi2d))) / (1 - np.sin(np.radians(phi2d)))
phi3d, c3d = 24.9, 17.9                          # FICHE:1069-1070

fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 3.1),
                             gridspec_kw={"width_ratios": [1, 1.25]})
ax.plot(s3, coulomb(s3, c_in, phi_in), ":", color="0.35", lw=1.1,
        label="critère d'entrée\n" r"($c=16{,}4$ MPa, $\varphi=23^\circ$)")
ax.plot(s3, ucs2d + k2d * s3, "-", color="#1f77b4", lw=1.1,
        label=r"rockim 2D adaptatif ($\varphi=22{,}8^\circ$)")
ax.plot(s3, coulomb(s3, c3d, phi3d), "-", color="#b03030", lw=0.9, alpha=.7)
ax.plot(s3_3d, s1_3d, "o", color="#b03030", ms=4,
        label="rockim 3D adaptatif\n" r"($\varphi=24{,}9^\circ$, $c=17{,}9$ MPa)")
ax.set_xlabel(r"$\sigma_3$ [MPa]")
ax.set_ylabel(r"$\sigma_1$ au pic [MPa]")
ax.xaxis.set_major_formatter(virg); ax.yaxis.set_major_formatter(virg)
ax.set_ylim(40, 195)
ax.legend(frameon=False, fontsize=6.5, loc="upper left")
ax.grid(alpha=.25, lw=.4)
ax.set_title("(a) enveloppe de Coulomb", fontsize=9)

# --------------------------------------------------------- (b) ecarts a l'article
# (libelle, rockim, Yan, source rockim)
rows = [
    (u"bande pesante : $U$ tête / éq. (20)", 0.600866, 0.60087, "CHG:33 ; YAN:864"),
    (u"brésilien : pic adaptatif / 500E", 5.30 / 4.44, 1.21, "FICHE:999-1003"),
    (u"compression simple adaptative (08/2026)", 51.1, 51.0, "FICHE:1004 ; YAN fig. 19b"),
    (u"compression simple adaptative (rejeu)", 51.43, 51.0, "B6:193"),
    (u"$\\varphi$ adaptatif", 22.8, 22.87, "FICHE:1006 ; YAN:1010"),
    (u"$\\varphi$ conventionnel 100E", 21.5, 22.36, "FICHE:1007 ; YAN:1011"),
    (u"angle de rupture $\\beta$ adaptatif", 66.0, 61.87, "FICHE:1044 ; YAN:961"),
    (u"Hopkinson : pic incident / $-V_0/c$", 0.931, 0.931415, "FICHE:1038 ; CHG:37"),
]
lab = [r[0] for r in rows]
ec = np.array([100 * (r[1] / r[2] - 1) for r in rows])
y = np.arange(len(rows))[::-1]
col = ["#b03030" if abs(e) > 2 else "#1f77b4" for e in ec]
bx.barh(y, ec, color=col, height=.55)
bx.axvline(0, color="k", lw=.6)
for yi, e in zip(y, ec):
    bx.text(e + (0.3 if e >= 0 else -0.3), yi, (("%+.2f %%" % e) if abs(e) < 0.1 else ("%+.1f %%" % e)).replace(".", ",").replace("-0,00", "0,00"),
            va="center", ha="left" if e >= 0 else "right", fontsize=7)
bx.set_yticks(y); bx.set_yticklabels(lab, fontsize=7.5)
bx.set_xlim(-6.5, 9.5)
bx.xaxis.set_major_formatter(virg)
bx.set_xlabel(u"écart rockim / référence [%]")
bx.grid(axis="x", alpha=.25, lw=.4)
bx.set_title(u"(b) écarts à la référence (Yan ou analytique)", fontsize=9)

fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig_yan_synthese.pdf")
fig.savefig(out)
print("ecrit", out)
for r, e in zip(rows, ec):
    print("%-45s %8.3f %8.3f %+6.2f %%  %s" % (r[0], r[1], r[2], e, r[3]))
