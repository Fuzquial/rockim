#!/usr/bin/env python3
"""Forces de coupe d'archive (aout 2026), panneau (d) des planches de l'etude.

Recadre le panneau (d) « forces de coupe » des planches PNG archivees dans
phd_geothermie/FDEM/rockim/coupe_pdc/figures/ pour quatre runs (v3, nogc, epfl,
a2). Le titre general des planches d'archive est faux (« passe 2.216 mm, loi de
Ye » : 2,216 mm est la profondeur de l'entaille du maillage, la passe vaut
1,016 mm et la loi est celle de Heilman) ; il est donc ecarte.

    python3 fig_coupe_archives.py  ->  fig_coupe_archives.pdf
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 9})
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
from PIL import Image

SRC = "/home/user/phd_geothermie/FDEM/rockim/coupe_pdc/figures"
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = [("v3", "(a) v3, pénalité"),
        ("nogc", "(b) nogc, contact général coupé"),
        ("epfl", r"(c) epfl, $k^- = k^+(D)$"),
        ("a2", "(d) a2, epfl + écrêtage en impulsion")]
BOX = (40, 660, 860, 1385)          # panneau (d), pixels de la planche 2557 x 1395

fig, axs = plt.subplots(2, 2, figsize=(7.0, 4.6))
for ax, (run, titre) in zip(axs.flat, RUNS):
    im = Image.open(f"{SRC}/out_cut_{run}_planche.png").convert("RGB").crop(BOX)
    ax.imshow(im)
    ax.set_axis_off()
    ax.set_title(titre, fontsize=9)
fig.tight_layout(pad=0.4)
fig.savefig(os.path.join(HERE, "fig_coupe_archives.pdf"), dpi=300)
print("ecrit fig_coupe_archives.pdf")
