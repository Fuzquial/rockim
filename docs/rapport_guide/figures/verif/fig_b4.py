"""Bilan d'energie B4 sur neuf calculs courts rejoues le 2026-10-06.
(a) residu |r| rapporte a l'echelle du bilan (max de Ec0, Ec et du flux brut
echange), verdict imprime par rockim ; (b) poste « integration » (correction
saute-mouton f^2 dt^2/2m) rapporte au travail fourni par l'outil au solide :
le bilan ferme dans tous les cas, y compris quand ce poste depasse de deux
ordres de grandeur l'apport de l'outil. Donnees : b4_residus.csv, b4_postes.csv
(prep_b4.py). Usage : python fig_b4.py -> fig_b4.pdf
"""
import os, csv
import numpy as np
from style_verif import plt, fmt

HERE = os.path.dirname(os.path.abspath(__file__))
r = list(csv.DictReader(open(os.path.join(HERE, "b4_residus.csv"))))
p = {x["cle"]: x for x in csv.DictReader(open(os.path.join(HERE, "b4_postes.csv")))}

fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw=dict(width_ratios=[1.45, 1]))
a = ax[0]
y = np.arange(len(r))[::-1]
pct = np.array([float(x["residu_pct_echelle"]) for x in r])
col = ["#c00000" if x["verdict"].startswith("CHECK") else "#1f4e79" for x in r]
VERDICT = {"OK": "conforme", "CHECK": "à vérifier", "OK (zero machine)": "zéro machine"}


def cas_lisible(c):
    """Libelle du lecteur : mode de calcul en clair, sans code de banc."""
    c = c.replace("raclage T1", "raclage").replace("(fdem3d, 2 µs)", "(FDEM 3D, 2 µs)")
    c = c.replace(" (fdem3d)", " (FDEM 3D)").replace(" (fdem)", " (FDEM 2D)")
    return c.replace(".", ",")

a.barh(y, pct, color=col, height=0.6)
a.set_xscale("log")
a.set_xlim(1e-14, 1e4)
a.set_yticks(y)
a.set_yticklabels([cas_lisible(x["cas"]) for x in r], fontsize=7)
a.axvline(1.0, color="k", lw=0.6, ls=":")
a.text(1.3, y[0] + 0.45, "seuil 1 %", fontsize=7)
for yi, x, v in zip(y, r, pct):
    a.text(v * 2, yi, VERDICT.get(x["verdict"], x["verdict"]), va="center", fontsize=6.5)
a.set_xlabel(r"$|r|$ / échelle (%)")
a.set_title("(a) résidu du bilan", fontsize=9)

b = ax[1]
keys = [("perc3d_20us", "percussion 3D"), ("t1_sig", "raclage, Signorini"), ("t1_pen", "raclage, pénalité"),
        ("perc2d", "percussion 2D")]
ratio = [float(p[k]["integration"]) / float(p[k]["outil->solide"]) for k, _ in keys]
yy = np.arange(len(keys))[::-1]
b.barh(yy, ratio, color=["#1f4e79", "#1f4e79", "#c55a11", "#c00000"], height=0.6)
b.set_xscale("log")
b.set_xlim(3e-3, 3e3)
b.set_yticks(yy)
b.set_yticklabels([l for _, l in keys], fontsize=7.5)
for yi, v in zip(yy, ratio):
    b.text(v * 1.4, yi, fmt(float(f"{v:.3g}")), va="center", fontsize=7)
b.axvline(1.0, color="k", lw=0.6, ls=":")
b.set_xlabel(r"$W_\mathrm{int\acute{e}gr} / W_\mathrm{outil \to solide}$")
b.set_title("(b) poste d'intégration", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_b4.pdf"))
