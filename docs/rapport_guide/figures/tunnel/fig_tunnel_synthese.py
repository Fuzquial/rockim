#!/usr/bin/env python3
"""Figure de synthese du banc tunnel EDZ (Wang et al. 2024) -- donnees archivees.

Aucune valeur n'est saisie a la main : tout est relu dans les fichiers de
metriques de l'etude (rockim/tunnel_edz/*_metrics.json et
phd_geothermie/FDEM/rockim/tunnel_edz/*_metrics.txt). Seules les valeurs
publiees par Wang et al. (rayon d'EDZ, lu sur leur fig. 12f, recopie dans
rockim/tunnel_edz/README.md:102) sont reprises du README.

    python3 fig_tunnel_synthese.py  ->  fig_tunnel_synthese.pdf

(a) rayon d'EDZ en fonction de sigma_0 : rockim (max et p95) contre Wang ;
(b) facteur de pointe insertionTipFactor : part de propagation (post-
    traitement nucleation_vs_propagation.py, bilan sur 0,08-0,70 s) ;
(c) meme balayage a t = 0,41 s : joints rompus apparies et blocs ;
(d) amortissement local 0,15 contre 0,7 (maillage de fumee, T = 0,10 s).
"""
import json
import os
import re

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
# --- fin du style commun

EDZ = "/home/user/rockim/tunnel_edz"
PHD = "/home/user/phd_geothermie/FDEM/rockim/tunnel_edz"
HERE = os.path.dirname(os.path.abspath(__file__))


def virgule(x, nd=1):
    return f"{x:.{nd}f}".replace(".", ",")


def fmt_axe(ax, nd=None):
    """Decimales a virgule sur les axes."""
    from matplotlib.ticker import FuncFormatter
    def f(v, _):
        if nd is None:
            s = f"{v:g}"
        else:
            s = f"{v:.{nd}f}"
        return s.replace(".", ",")
    ax.yaxis.set_major_formatter(FuncFormatter(f))


def js(name):
    with open(os.path.join(EDZ, name)) as f:
        return json.load(f)


def txt(name):
    with open(os.path.join(PHD, name), encoding="utf-8", errors="replace") as f:
        return f.read()


# ---- (a) balayage sigma_0 ---------------------------------------------------
sig = [3, 4, 5, 6, 7]
runs = {3: "s3", 4: "s4", 5: "ref_iso", 6: "s6", 7: "s7"}
rmax = [js(runs[s] + "_metrics.json")["edz_radius_max_m"] for s in sig]
rp95 = [js(runs[s] + "_metrics.json")["edz_radius_p95_m"] for s in sig]
# Wang et al. (2024) fig. 12f, recopie dans rockim/tunnel_edz/README.md:102
readme = open(os.path.join(EDZ, "README.md"), encoding="utf-8").read()
m = re.search(r"EDZ ([\d,]+) / ([\d,]+) / ([\d,]+) / ([\d,]+) / ([\d,]+) m", readme)
wang = [float(v.replace(",", ".")) for v in m.groups()]

# ---- (b,c) facteur de pointe -----------------------------------------------
bilan = open("/home/user/rockim/BILAN_insertion_adaptative.md", encoding="utf-8").read()
m = re.search(r"\| propagation \| ([\d,]+) % \| ([\d,]+) % \| \*\*([\d,]+) %\*\* \| "
              r"([\d,]+) % \| ([\d,]+) % \|", bilan)
prop = [float(v.replace(",", ".")) for v in m.groups()]   # 1,0 1,3 1,6 2,0 | intr
facteurs = [1.0, 1.3, 1.6, 2.0]


def blocs(t, run):
    """Bloc 'run trame 20' de block_sizes.py dans un *_metrics.txt."""
    i = t.index(run + " ")
    i = t.index("trame 20", i)
    seg = t[i:i + 600]
    rompus = int(re.search(r"appariees : (\d+)", seg).group(1))
    nb = int(re.search(r"blocs :\s+(\d+)", seg).group(1))
    mono = float(re.search(r"mono-element :\s+\d+\s+\(\s*([\d.]+) %", seg).group(1))
    return rompus, nb, mono


t13 = txt("tip13_metrics.txt")
c_ref = blocs(t13, "out_tun_ref_stab")
c13 = blocs(t13, "out_tun_tip13")
c16 = blocs(txt("tip16_metrics.txt"), "out_tun_tip16")
c20 = blocs(txt("tip20_metrics.txt"), "out_tun_tip20")
cb = [c_ref, c13, c16, c20]

# ---- (d) amortissement -------------------------------------------------------
d15 = js("smoke_metrics.json")
d70 = js("smoke_d07_metrics.json")

# ---- figure -----------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.6))
gris, rouge, bleu, vert = "#4d4d4d", "#b2182b", "#2166ac", "#1b7837"

a = ax[0, 0]
a.plot(sig, wang, "s--", color=gris, label="Wang et al. (2024), fig. 12f")
a.plot(sig, rmax, "o-", color=rouge, label="rockim, rayon maximal")
a.plot(sig, rp95, "^-", color=bleu, label="rockim, quantile 95 %")
a.set_xlabel(r"contrainte in situ $\sigma_0$ [MPa]")
a.set_ylabel("rayon de l'EDZ [m]")
a.set_title(r"(a) rayon de l'EDZ, $\lambda = 1$, $T = 0{,}25$ s", fontsize=9)
a.legend(fontsize=7.5, frameon=False)
fmt_axe(a)
a.grid(alpha=.3)

a = ax[0, 1]
a.plot(facteurs, prop[:4], "o-", color=rouge, label="adaptatif + pointe relâchée")
a.axhline(prop[4], ls="--", color=bleu, label=f"intrinsèque ({virgule(prop[4])} %)")
for f, p in zip(facteurs, prop[:4]):
    a.annotate(virgule(p), (f, p), textcoords="offset points", xytext=(0, -12 if f == 1.6 else 5),
               ha="center", fontsize=7.5)
a.set_xlabel(r"facteur de pointe $k_\mathrm{tip}$")
a.set_ylabel("part de propagation [%]")
a.set_title("(b) insertions en pointe / insertions totales", fontsize=9)
a.set_ylim(35, 65)
a.set_xticks(facteurs)
a.set_xticklabels([virgule(f) for f in facteurs])
a.legend(fontsize=7.5, frameon=False, loc="lower right")
a.grid(alpha=.3)

a = ax[1, 0]
x = range(4)
w = 0.38
a.bar([i - w / 2 for i in x], [c[0] / 1e3 for c in cb], w, color=rouge,
      label="joints rompus (milliers)")
a.bar([i + w / 2 for i in x], [c[1] / 1e3 for c in cb], w, color=bleu,
      label="blocs (milliers)")
for i, c in enumerate(cb):
    a.text(i, max(c[0], c[1]) / 1e3 + 1.2, f"{virgule(c[2])} %\nmono-élém.",
           ha="center", fontsize=7)
a.set_xticks(list(x))
a.set_xticklabels([virgule(f) for f in facteurs])
a.set_xlabel(r"facteur de pointe $k_\mathrm{tip}$")
a.set_ylim(0, 52)
a.set_title(r"(c) découpage en blocs à $t = 0{,}41$ s ($r \leq 25$ m)", fontsize=9)
a.legend(fontsize=7.5, frameon=False, loc="upper left")
a.grid(alpha=.3, axis="y")

a = ax[1, 1]
cles = [("broken", "joints\nrompus"), ("crack_length_m", "longueur\nfissurée"),
        ("edz_radius_p95_m", "EDZ (quantile\n95 %)"), ("u_max_m", "déplacement\nmax.")]
ratios = [d70[k] / d15[k] for k, _ in cles]
a.bar(range(4), ratios, color=[rouge, rouge, bleu, vert])
for i, (k, _) in enumerate(cles):
    a.text(i, ratios[i] + 0.03, f"{virgule(ratios[i], 2)}", ha="center", fontsize=7.5)
a.axhline(1, color="k", lw=.8)
a.set_xticks(range(4))
a.set_xticklabels([l for _, l in cles], fontsize=7.5)
a.set_ylabel("rapport amortissement 0,7 / 0,15")
a.set_ylim(0, 1.15)
fmt_axe(a, 1)
a.set_title("(d) effet de l'amortissement local, maillage d'essai", fontsize=9)
a.grid(alpha=.3, axis="y")

fig.tight_layout()
out = os.path.join(HERE, "fig_tunnel_synthese.pdf")
fig.savefig(out)
fig.savefig(out.replace(".pdf", ".png"), dpi=150)
print("ecrit :", out)
print("Wang", wang, "rmax", [round(v, 2) for v in rmax], "p95", [round(v, 2) for v in rp95])
print("prop", prop, "blocs", cb)
print("amort.", {k: (d15[k], d70[k]) for k, _ in cles})
