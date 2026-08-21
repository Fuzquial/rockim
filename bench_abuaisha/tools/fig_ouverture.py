#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_ouverture.py — L'OUVERTURE DES FISSURES, reconstruite.
#
#   python bench_abuaisha/tools/fig_ouverture.py out_hfs_aniso --sh 6.8 --sv 4.6
#
# rockim N'ECRIT PAS l'ouverture : le writer ne pousse qu'UNE levre par joint
# (`lines.push_back({J.a1, J.a2})`, FdemSolver.cpp). Il faut donc retrouver la
# seconde, et c'est possible sans toucher au solveur :
#
#   * les noeuds sont dedoubles par element, n = 3e + k (verifie) ;
#   * a l'instant initial, les copies d'un meme sommet geometrique sont
#     confondues, ce qui donne un identifiant de sommet par arrondi ;
#   * chaque arete geometrique est alors portee par EXACTEMENT deux elements
#     (284 124 aretes pour 284 124 joints, plus 213 aretes de bord) ;
#   * la levre opposee d'un joint est l'arete de l'AUTRE element.
#
# L'ouverture vaut ensuite (b - a).n a chaque extremite, avec n la normale du
# segment dans la configuration COURANTE — la meme formule que celle dont le
# solveur se sert pour le volume de cavite (updateWetBoundary).
#
# Une ouverture NEGATIVE est une interpenetration : on la garde telle quelle
# dans les profils, on ne la masque pas.
# ---------------------------------------------------------------------------
import argparse
import csv
import glob
import io
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tunnel_edz"))
from plot_tunnel_fields import complete  # noqa: E402

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 10.5,
    "axes.titlesize": 11,
    "axes.linewidth": 0.8,
    "savefig.dpi": 200,
})

CX, CY, RB = 4.0, 4.0, 0.05


def vg(x, n=2):
    return (("%." + str(n) + "f") % x).replace(".", ",")


def grab(txt, nm):
    m = re.search('Name="' + nm + '"[^>]*>(.*?)</DataArray>', txt, re.S)
    return np.fromstring(m.group(1), sep=" ")


def points(txt):
    m = re.search(r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>", txt, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


def appariement(P0, JC):
    """Rend, pour chaque joint, les indices (b1, b2) de la levre opposee."""
    _, gid = np.unique(np.round(P0, 9), axis=0, return_inverse=True)
    nE = len(P0) // 3
    n = np.arange(3 * nE).reshape(nE, 3)
    u = np.concatenate([n[:, 0], n[:, 1], n[:, 2]])
    v = np.concatenate([n[:, 1], n[:, 2], n[:, 0]])
    gu, gv = gid[u], gid[v]
    lo, hi = np.minimum(gu, gv), np.maximum(gu, gv)
    key = lo.astype(np.int64) * (gid.max() + 1) + hi
    o = np.argsort(key, kind="stable")
    ks, us, vs = key[o], u[o], v[o]
    # les aretes internes viennent par paires consecutives ; on repere les
    # debuts de groupe et on ne garde que les groupes de taille 2
    start = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
    cnt = np.diff(np.r_[start, len(ks)])
    pair = start[cnt == 2]
    kpair = ks[pair]
    e0, e1 = us[pair] // 3, us[pair + 1] // 3

    jk = (np.minimum(gid[JC[:, 0]], gid[JC[:, 1]]).astype(np.int64)
          * (gid.max() + 1) + np.maximum(gid[JC[:, 0]], gid[JC[:, 1]]))
    pos = np.searchsorted(kpair, jk)
    ok = (pos < len(kpair)) & (kpair[np.clip(pos, 0, len(kpair) - 1)] == jk)
    if not ok.all():
        print("  (%d joints sans arete jumelle — ignorés)" % int((~ok).sum()))
    eA = JC[:, 0] // 3
    # le partenaire est l'entree du couple dont l'element n'est PAS celui de A
    take1 = e0[np.clip(pos, 0, len(kpair) - 1)] == eA
    pu = np.where(take1, us[pair + 1][np.clip(pos, 0, len(kpair) - 1)],
                  us[pair][np.clip(pos, 0, len(kpair) - 1)])
    pv = np.where(take1, vs[pair + 1][np.clip(pos, 0, len(kpair) - 1)],
                  vs[pair][np.clip(pos, 0, len(kpair) - 1)])
    # aligner sur les sommets geometriques de A
    same = gid[pu] == gid[JC[:, 0]]
    b1 = np.where(same, pu, pv)
    b2 = np.where(same, pv, pu)
    return b1, b2, ok


def ouverture(P, JC, b1, b2):
    A1, A2 = P[JC[:, 0]], P[JC[:, 1]]
    e = A2 - A1
    L = np.hypot(e[:, 0], e[:, 1])
    L = np.where(L < 1e-14, 1.0, L)
    nx, ny = e[:, 1] / L, -e[:, 0] / L
    w1 = (P[b1][:, 0] - A1[:, 0]) * nx + (P[b1][:, 1] - A1[:, 1]) * ny
    w2 = (P[b2][:, 0] - A2[:, 0]) * nx + (P[b2][:, 1] - A2[:, 1]) * ny
    return 0.5 * (w1 + w2), L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--sh", type=float, default=6.8)
    ap.add_argument("--sv", type=float, default=4.6)
    ap.add_argument("--dseuil", type=float, default=0.5)
    ap.add_argument("--stem", default=None)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(ROOT, a.run)

    h = list(csv.DictReader(open(os.path.join(run, "history.csv"))))
    th = np.array([float(r["t"]) for r in h]) * 1e3
    ph = np.array([float(r["hydroP"]) for r in h]) / 1e6
    Vh = np.array([float(r["hydroVol"]) for r in h])

    def kof(f):
        return int(os.path.basename(f).rsplit("_", 1)[1].split(".")[0])

    js = sorted([f for f in glob.glob(
        os.path.join(run, "fdem_joints_[0-9]*.vtu")) if complete(f)], key=kof)
    frw = list(csv.DictReader(open(os.path.join(run, "frames.csv"))))
    tfr = np.array([float(r["t"]) for r in frw]) * 1e3

    t0 = io.open(js[0], errors="ignore").read()
    P0 = points(t0)
    JC = grab(t0, "connectivity").astype(np.int64).reshape(-1, 2)
    print("appariement des lèvres…")
    b1, b2, ok = appariement(P0, JC)
    del t0

    txt = io.open(js[-1], errors="ignore").read()
    P = points(txt)
    D = grab(txt, "damage")
    w, L = ouverture(P, JC, b1, b2)
    tf = tfr[kof(js[-1])]
    jj = int(np.argmin(np.abs(th - tf)))

    m = (D > a.dseuil) & ok
    mid = 0.5 * (P[JC[m, 0]] + P[JC[m, 1]]) - [CX, CY]
    seg = (P[JC[m]] - [CX, CY]) * 1e3
    wm = w[m] * 1e6                                   # micrometres
    x = mid[:, 0] * 1e3

    fig = plt.figure(figsize=(12.6, 5.0))
    gs = fig.add_gridspec(1, 2, wspace=0.26, width_ratios=[1.0, 1.25])

    # --- (a) la carte -------------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.add_patch(plt.Circle((0.0, 0.0), RB * 1e3, fc="0.9", ec="0.5", lw=0.7,
                            zorder=1))
    vmx = float(np.percentile(np.abs(wm), 99)) if len(wm) else 1.0
    lc = LineCollection(seg, array=wm, cmap="viridis", lw=2.4, clim=(0, vmx),
                        zorder=2)
    ax.add_collection(lc)
    cb = fig.colorbar(lc, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label(r"ouverture   [$\mu$m]")
    cb.outline.set_linewidth(0.6)
    z = 1.1 * float(np.abs(seg).max()) if len(seg) else 100.0
    ax.set_xlim(-z, z)
    ax.set_ylim(-z, z)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$x$   [mm]")
    ax.set_ylabel(r"$y$   [mm]")
    ax.set_title("(a)  Ouverture le long des ailes", loc="left")

    # --- (b) le profil ------------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(x, wm, ".", ms=3.5, color="#1f4e79", alpha=0.7,
            label="joints rompus  ($D$ > " + vg(a.dseuil, 1) + ")")
    ax.axhline(0.0, color="0.5", lw=0.8)
    # la bande grise = le FORAGE projete sur l'axe des abscisses, |x| <= R.
    # Les points qui y tombent sont des joints de la PAROI, a y non nul : la
    # bande marque une projection, pas une zone sans donnee. L'etiquette va en
    # BAS du cadre — en haut elle finissait sous la legende.
    ax.axvspan(-RB * 1e3, RB * 1e3, color="0.88", zorder=0)
    ax.annotate("forage" + chr(10) + r"$\varnothing$ 100 mm", (0.0, 0.0),
                xycoords=("data", "axes fraction"), xytext=(0, 8),
                textcoords="offset points", ha="center", va="bottom",
                fontsize=9, color="0.35")

    # reference de Sneddon : fissure de demi-longueur Lc sous pression NETTE
    Lc = float(np.abs(x).max()) / 1e3
    dp = (ph[jj] - a.sv) * 1e6
    E, nu = 35.0e9, 0.27
    xs = np.linspace(-Lc, Lc, 400)
    ws = 4.0 * dp * (1 - nu * nu) / E * np.sqrt(np.maximum(Lc**2 - xs**2, 0.0))
    ax.plot(xs * 1e3, ws * 1e6, "-", lw=1.4, color="#c0392b",
            label=(r"Sneddon, $\Delta p$ = " + vg(dp / 1e6) + " MPa,  $L$ = "
                   + vg(Lc * 1e3, 0) + " mm"))
    ax.set_xlabel(r"abscisse depuis l'axe du forage   $x$   [mm]")
    ax.set_ylabel(r"ouverture   [$\mu$m]")
    ax.set_title("(b)  Profil d'ouverture", loc="left")
    ax.grid(alpha=0.22, lw=0.6)
    ax.legend(fontsize=9, loc="upper right")

    # --- recoupement : volume des fissures ---------------------------------
    vol = float((np.maximum(w[m], 0.0) * L[m]).sum())      # m2/m
    dV = float(Vh[jj] - Vh[0])
    reg = "anisotrope" if abs(a.sh - a.sv) > 1e-9 else "isotrope"
    fig.suptitle("Ouverture des fissures — état de contrainte " + reg
                 + "   ($t$ = " + vg(tf) + " ms,  $p$ = " + vg(ph[jj])
                 + " MPa)", fontsize=12.5, y=0.985)

    stem = a.stem or os.path.join(
        ROOT, "bench_abuaisha",
        "fig_ouverture_" + os.path.basename(run).replace("out_hfs_", ""))
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    fig.savefig(stem + ".png", bbox_inches="tight")
    print("écrit :", stem + ".pdf  et  .png")
    print("  %d joints rompus | ouverture max %s um, mediane %s um"
          % (int(m.sum()), vg(wm.max(), 1), vg(float(np.median(wm)), 1)))
    print("  ouverture NEGATIVE (interpénétration) : %d joints, min %s um"
          % (int((wm < 0).sum()), vg(float(wm.min()), 1)))
    print("  Sneddon au centre : %s um   (mesuré au max : %s um)"
          % (vg(float(ws.max()) * 1e6, 1), vg(wm.max(), 1)))
    print("  RECOUPEMENT volume : somme(w.L) = %.3e m2/m contre "
          "V - V0 = %.3e  (écart %+.0f %%)"
          % (vol, dV, (vol / dV - 1) * 100 if dV else float("nan")))


if __name__ == "__main__":
    main()
