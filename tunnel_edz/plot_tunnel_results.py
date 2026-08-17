#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# plot_tunnel_results.py — planche de dépouillement d'un run tunnel, calquée
# sur les figures de Wang et al. (2024) :
#
#   (a) carte des fissures classées par MODE          (leur fig. 7)
#   (b) contrainte principale MAJEURE sigma_1          (leur fig. 10a)
#   (c) cisaillement maximal tau_max = (s1 - s3)/2     (leur fig. 10b)
#   (d) module du déplacement + contour d'EDZ          (leurs fig. 9 / 11)
#   (e) nombre de fissures contre le temps             (leur fig. 8)
#
#   python tunnel_edz/plot_tunnel_results.py [out_tun_smoke] [--zoom 20]
#
# Code couleur des fissures repris de l'article : traction en vert,
# cisaillement en rouge, mixte en bleu.
# ---------------------------------------------------------------------------
import argparse
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
from make_unstructured_mesh import TUNNEL_HS  # noqa: E402
sys.path.insert(0, HERE)
from plot_tunnel_mesh import profile_xy  # noqa: E402

C_TEN, C_MIX, C_SHR = "#2E9E4F", "#2C6FB5", "#C8342B"


def vtu(path, arrays):
    with open(path) as f:
        txt = f.read()
    P = np.fromstring(re.search(
        r'<Points>.*?<DataArray[^>]*>(.*?)</DataArray>', txt, re.S).group(1),
        sep=" ").reshape(-1, 3)[:, :2]
    C = np.fromstring(re.search(
        r'Name="connectivity"[^>]*>(.*?)</DataArray>', txt, re.S).group(1),
        sep=" ").astype(int).reshape(-1, 3)
    out = {}
    for nm in arrays:
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % nm, txt, re.S)
        out[nm] = np.fromstring(m.group(1), sep=" ") if m else None
    return P, C, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_tun_smoke")
    ap.add_argument("--zoom", type=float, default=20.0)
    ap.add_argument("--qmix", type=float, default=0.5)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(HERE, "..", a.run)
    out = a.out or os.path.join(HERE, os.path.basename(a.run) + "_planche.png")

    frames = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
              if "joints" not in os.path.basename(f)]
    P0, C, _ = vtu(frames[0], [])
    P, C, S = vtu(frames[-1], ["sigmaXX", "sigmaYY", "sigmaXY"])
    W, H = P0[:, 0].max(), P0[:, 1].max()
    cx, cy = 0.5 * W, 0.5 * H
    cy0 = cy - 0.5 * TUNNEL_HS["height"]
    px, py, _ = profile_xy(cx, cy0)

    # champs de contrainte : valeurs principales du tenseur plan
    sxx, syy, sxy = S["sigmaXX"], S["sigmaYY"], S["sigmaXY"]
    tr, dif = 0.5 * (sxx + syy), 0.5 * (sxx - syy)
    rad = np.hypot(dif, sxy)
    s1, s3 = (tr + rad) / 1e6, (tr - rad) / 1e6      # MPa, tension positive
    tmax = rad / 1e6

    # deplacement nodal -> par element (moyenne des 3 sommets)
    umag = np.linalg.norm(P - P0, axis=1)
    uel = umag[C].mean(axis=1)

    T0 = mtri.Triangulation(P0[:, 0], P0[:, 1], C)   # champs sur la config. initiale

    # fissures
    J = np.genfromtxt(os.path.join(run, "fdem_final_joints.csv"), delimiter=",",
                      names=True, invalid_raise=False)
    brk = (J["tBreak"] >= 0.0) | (J["damage"] >= 0.999)
    seg = np.stack([np.c_[J["x1"][brk], J["y1"][brk]],
                    np.c_[J["x2"][brk], J["y2"][brk]]], axis=1)
    rn, rs = J["rn"][brk], J["rs"][brk]
    mx = np.maximum(rn, rs)
    q = np.where(mx > 0, np.minimum(rn, rs) / np.maximum(mx, 1e-300), 0.0)
    mixed = q >= a.qmix
    tens = (~mixed) & (rn >= rs)
    shear = (~mixed) & (rn < rs)
    xm, ym = seg[:, :, 0].mean(axis=1), seg[:, :, 1].mean(axis=1)
    redz = np.percentile(np.hypot(xm - cx, ym - cy), 95)

    Z = a.zoom
    fig = plt.figure(figsize=(16.5, 9.2))

    def frame(ax, title):
        ax.plot(px, py, color="k", lw=1.2)
        ax.set_aspect("equal")
        ax.set_xlim(cx - Z, cx + Z)
        ax.set_ylim(cy - Z, cy + Z)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("x [m]")

    ax = fig.add_subplot(2, 3, 1)
    ax.add_collection(LineCollection(seg[shear], colors=C_SHR, lw=0.7))
    ax.add_collection(LineCollection(seg[mixed], colors=C_MIX, lw=0.7))
    ax.add_collection(LineCollection(seg[tens], colors=C_TEN, lw=0.7))
    ax.add_patch(plt.Circle((cx, cy), redz, fill=False, ls="--", color="#B8860B",
                            lw=1.4))
    frame(ax, f"(a) fissures : {int(tens.sum())} traction / "
              f"{int(mixed.sum())} mixte / {int(shear.sum())} cisaillement")
    ax.set_ylabel("y [m]")
    for c, lab in ((C_TEN, "traction"), (C_MIX, "mixte"), (C_SHR, "cisaillement")):
        ax.plot([], [], color=c, lw=2, label=lab)
    ax.plot([], [], ls="--", color="#B8860B", label=f"EDZ p95 = {redz:.1f} m")
    ax.legend(loc="lower left", fontsize=7, framealpha=0.9)

    ax = fig.add_subplot(2, 3, 2)
    v = np.percentile(np.abs(s1), 99)
    h = ax.tripcolor(T0, facecolors=s1, cmap="RdBu_r", vmin=-v, vmax=v)
    frame(ax, r"(b) contrainte principale majeure $\sigma_1$ [MPa]")
    fig.colorbar(h, ax=ax, shrink=0.85)

    ax = fig.add_subplot(2, 3, 3)
    h = ax.tripcolor(T0, facecolors=tmax, cmap="inferno",
                     vmax=np.percentile(tmax, 99.5))
    frame(ax, r"(c) cisaillement maximal $\tau_{max}$ [MPa]")
    fig.colorbar(h, ax=ax, shrink=0.85)

    ax = fig.add_subplot(2, 3, 4)
    h = ax.tripcolor(T0, facecolors=uel, cmap="turbo")
    ax.add_patch(plt.Circle((cx, cy), redz, fill=False, ls="--", color="w", lw=1.4))
    frame(ax, f"(d) deplacement |u| [m] — max {umag.max():.3f} m")
    ax.set_ylabel("y [m]")
    fig.colorbar(h, ax=ax, shrink=0.85)

    ax = fig.add_subplot(2, 3, 5)
    d = np.genfromtxt(os.path.join(run, "history.csv"), delimiter=",",
                      names=True, invalid_raise=False)
    ax.plot(d["t"] * 1e3, d["nBroken"], color=C_SHR, lw=1.6)
    ax.set_xlabel("temps [ms]")
    ax.set_ylabel("joints rompus")
    ax.set_title("(e) cinetique de fissuration (leur fig. 8)", fontsize=11)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(2, 3, 6)
    ax.add_collection(LineCollection(seg[shear], colors=C_SHR, lw=0.5))
    ax.add_collection(LineCollection(seg[mixed], colors=C_MIX, lw=0.5))
    ax.add_collection(LineCollection(seg[tens], colors=C_TEN, lw=0.5))
    ax.plot(px, py, color="k", lw=1.0)
    ax.set_aspect("equal")
    ax.set_xlim(cx - 0.55 * W, cx + 0.55 * W)
    ax.set_ylim(cy - 0.55 * H, cy + 0.55 * H)
    ax.add_patch(plt.Rectangle((0, 0), W, H, fill=False, color="0.5", lw=0.8))
    ax.set_title("(f) le massif entier : l'EDZ reste loin des bords", fontsize=11)
    ax.set_xlabel("x [m]")

    fig.suptitle(f"{os.path.basename(run)} — in situ 5 MPa, lambda = 1 "
                 f"(maillage smoke : zone fine a 6 m de la paroi)", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=155)
    print("ecrit :", out)


if __name__ == "__main__":
    main()
