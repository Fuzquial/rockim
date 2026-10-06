"""Style commun des figures de verification : STIX, libelles francais, virgule decimale."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

plt.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 9,
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.1,
    "legend.fontsize": 7.5,
    "legend.frameon": False,
    "axes.unicode_minus": True,
    "savefig.bbox": "tight",
})

COUL = {"plastic": "#1f4e79", "origin": "#c55a11", "solidity": "#7f7f7f"}
NOM = {"plastic": r"$\mathtt{plastic}$", "origin": r"$\mathtt{origin}$ + cliquet",
       "solidity": r"$\mathtt{solidity}$"}


def fmt(x, nd=None):
    """Nombre au format francais : virgule decimale, signe moins typographique."""
    if nd is None:
        s = f"{x:g}"
    else:
        s = f"{x:.{nd}f}"
    return s.replace("-", "−").replace(".", ",")


def virgule(ax, axes="xy"):
    f = FuncFormatter(lambda x, pos: fmt(round(x, 10)))
    if "x" in axes and ax.get_xscale() == "linear":
        ax.xaxis.set_major_formatter(f)
    if "y" in axes and ax.get_yscale() == "linear":
        ax.yaxis.set_major_formatter(f)
