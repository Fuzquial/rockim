"""Style commun des figures de verification : Computer Modern, libelles francais, virgule decimale."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

plt.rcParams.update({
    "font.size": 9,
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.1,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "savefig.bbox": "tight",
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

COUL = {"plastic": "#1f4e79", "origin": "#c55a11", "solidity": "#7f7f7f"}
NOM = {"plastic": r"$\mathtt{plastic}$", "origin": r"$\mathtt{origin}$ + cliquet",
       "solidity": r"$\mathtt{solidity}$"}


def fmt(x, nd=None):
    """Nombre au format francais : virgule decimale, signe moins typographique."""
    if nd is None:
        s = f"{x:g}"
    else:
        s = f"{x:.{nd}f}"
    return s.replace(".", ",")


def virgule(ax, axes="xy"):
    f = FuncFormatter(lambda x, pos: fmt(round(x, 10)))
    if "x" in axes and ax.get_xscale() == "linear":
        ax.xaxis.set_major_formatter(f)
    if "y" in axes and ax.get_yscale() == "linear":
        ax.yaxis.set_major_formatter(f)
