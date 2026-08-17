# -*- coding: utf-8 -*-
"""Configuration du modele de coupe PDC : geometrie, maillage GBM, outil,
conditions aux limites et heterogeneite de Weibull."""
import glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Polygon, FancyArrow

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
RUN = os.path.join(BASE, "runs", "demo_cutting")
W, H = 50.0, 30.0          # mm
DOC, RAKE, FACE = 2.0, 20.0, 13.0
TOOLX = -3.0


def arr(t, n):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, t, re.S)
    return np.fromstring(m.group(1), sep=" ") if m else None


def pts(t):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", t, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


f0 = sorted(glob.glob(os.path.join(RUN, "fdem_[0-9][0-9][0-9][0-9].vtu")))[0]
t0 = open(f0).read()
P = pts(t0) * 1e3
conn = arr(t0, "connectivity").astype(int).reshape(-1, 3)
grain = arr(t0, "grain")
stat = arr(t0, "ftScale")          # facteur de Weibull par element

fig = plt.figure(figsize=(15, 6.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.25, 1])


def cutter(ax, shown=7.0):
    """Trace le PDC : arete a (toolX, H-DOC), face de coupe inclinee du back
    rake (mesure depuis la VERTICALE). On n'affiche que la partie utile de la
    face (la face reelle fait 13 mm et sort du cadre)."""
    a = np.deg2rad(RAKE)
    tip = np.array([TOOLX, H - DOC])
    top = tip + shown * np.array([-np.sin(a), np.cos(a)])
    back = top + np.array([-3.5, 0.0])
    bot = tip + np.array([-3.5, 0.0])
    ax.add_patch(Polygon([tip, top, back, bot], closed=True,
                         facecolor="0.35", edgecolor="k", lw=1.2, zorder=5))
    ax.plot([tip[0], top[0]], [tip[1], top[1]], color="C1", lw=2.4, zorder=6)
    ax.annotate("", xy=(tip[0] + 9, tip[1]), xytext=(tip[0] + 2, tip[1]),
                arrowprops=dict(arrowstyle="-|>", lw=2.2, color="C3"), zorder=6)
    ax.text(tip[0] + 5.5, tip[1] + 1.4, "10 m/s", color="C3", fontsize=10,
            ha="center")


# --- (a) configuration ------------------------------------------------------
a = fig.add_subplot(gs[0])
a.add_patch(Polygon([[0, 0], [W, 0], [W, H], [0, H]], closed=True,
                    facecolor="#e8e4dc", edgecolor="k", lw=1.4))
a.plot([0, W], [H - DOC, H - DOC], "C0--", lw=1.4)
a.annotate("", xy=(W * 0.75, H), xytext=(W * 0.75, H - DOC),
           arrowprops=dict(arrowstyle="<->", color="C0"))
a.text(W * 0.77, H - DOC / 2, "profondeur de passe 2 mm", color="C0", fontsize=9)
cutter(a)
# conditions aux limites
a.plot([0, W], [0, 0], "k", lw=4)
for x in np.linspace(1, W - 1, 14):
    a.plot([x, x - 1.2], [0, -1.6], "k", lw=0.9)
a.text(W / 2, -3.6, "base encastrée", ha="center", fontsize=9)
for xx in (0, W):
    a.plot([xx, xx], [0, H], color="C2", lw=3, alpha=0.6)
a.text(-1.5, H / 2, "frontières absorbantes", color="C2", fontsize=9,
       rotation=90, va="center")
a.text(W / 2, H * 0.45, "granite\n50 × 30 mm\nGBM Voronoï 3 mm",
       ha="center", fontsize=10)
a.text(TOOLX - 7, H + 5.5, "face de coupe\n(back rake 20°)", fontsize=9,
       color="C1")
a.set_xlim(-11, W + 3); a.set_ylim(-6, H + 10); a.set_aspect("equal")
a.axis("off")
a.set_title("Configuration : cutter PDC, back rake 20°", fontsize=11)

# --- (b) maillage GBM -------------------------------------------------------
a = fig.add_subplot(gs[1])
rng = np.random.default_rng(3)
lut = rng.random((int(grain.max()) + 2, 3)) * 0.4 + 0.55
a.add_collection(PolyCollection([P[c] for c in conn],
                                facecolors=[lut[int(g)] for g in grain],
                                edgecolors="0.4", linewidths=0.12))
cutter(a)
a.set_xlim(-11, W + 3); a.set_ylim(-2, H + 10); a.set_aspect("equal")
a.set_xlabel("x (mm)"); a.set_ylabel("y (mm)")
a.set_title("Maillage : %d grains de Voronoï, %d éléments\n(Delaunay intra-grain "
            "1,2 mm ≈ 0,4 × grain)" % (int(grain.max()) + 1, len(conn)),
            fontsize=10)

# --- (c) heterogeneite : le facteur de Weibull vit sur les JOINTS ----------
from matplotlib.collections import LineCollection
a = fig.add_subplot(gs[2])
jf = f0.replace("fdem_", "fdem_joints_")
stat = None
if os.path.exists(jf):
    jt = open(jf).read()
    jP = pts(jt) * 1e3
    jc = arr(jt, "connectivity").astype(int)
    off = arr(jt, "offsets").astype(int)
    stat = arr(jt, "ftScale")
    segs, st = [], 0
    for i, en in enumerate(off):
        ii = jc[st:en]; st = en
        if len(ii) >= 2:
            segs.append(jP[ii[:2]])
    lc = LineCollection(segs, array=stat[:len(segs)], cmap="RdYlGn", lw=1.3)
    a.add_collection(lc)
    cb = fig.colorbar(lc, ax=a, shrink=0.75)
    cb.set_label("facteur de résistance du joint")
    ttl = ("Hétérogénéité de Weibull (m = 6)\n%d joints, facteur %.2f à %.2f"
           "\n(rouge = points faibles → amorçage)"
           % (len(segs), stat.min(), stat.max()))
else:
    ttl = "fichier de joints absent"
a.set_xlim(0, W); a.set_ylim(0, H); a.set_aspect("equal")
a.set_xlabel("x (mm)"); a.set_title(ttl, fontsize=10)

fig.suptitle("Simulation de coupe PDC en FDEM 2D — configuration du modèle",
             fontsize=13)
fig.tight_layout()
out = os.path.join(BASE, "figures", "cutting_setup.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
print("elements %d | grains %d | Weibull %.2f-%.2f"
      % (len(conn), int(grain.max()) + 1,
         stat.min() if stat is not None else 0,
         stat.max() if stat is not None else 0))
