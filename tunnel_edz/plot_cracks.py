#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# plot_cracks.py — la carte des fissures SEULE, a deux echelles, pour juger du
# facies (leur fig. 7). Rien d'autre sur l'image : pas de champ, pas de
# maillage, juste les traces de rupture et le contour du tunnel.
#
#   python tunnel_edz/plot_cracks.py out_tun_ref_iso [--frame -1]
#                                    [--wide 15] [--tight 7]
#
# Couleurs du solveur : vert = rupture en traction (rn >= rs a l'instant de la
# rupture), rouge = cisaillement. Fonctionne sur un run EN COURS (les fissures
# sont lues dans le VTU de joints de la trame).
# ---------------------------------------------------------------------------
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
sys.path.insert(0, HERE)
from make_unstructured_mesh import TUNNEL_HS  # noqa: E402
from plot_tunnel_mesh import profile_xy  # noqa: E402
from plot_tunnel_fields import read_vtu, complete  # noqa: E402

C_TEN, C_SHR = "#1B8A3A", "#C8342B"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_tun_ref_iso")
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--wide", type=float, default=15.0)
    ap.add_argument("--tight", type=float, default=7.0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(HERE, "..", a.run)

    el = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
          if complete(f)]
    jn = [f for f in sorted(glob.glob(os.path.join(run, "fdem_joints_[0-9]*.vtu")))
          if complete(f)]
    k = min(len(el), len(jn)) - 1 if a.frame < 0 else a.frame
    tag = os.path.basename(jn[k])[12:16]
    out = a.out or os.path.join(HERE,
                                f"{os.path.basename(run)}_f{tag}_fissures.png")

    P0, _, _ = read_vtu(el[0], [])
    W, H = P0[:, 0].max(), P0[:, 1].max()
    cx, cy = 0.5 * W, 0.5 * H
    px, py, _ = profile_xy(cx, cy - 0.5 * TUNNEL_HS["height"])

    PJ, CJ, SJ = read_vtu(jn[k], ["damage", "breakMode"], ncell=2)
    bm, dmg = SJ["breakMode"], SJ["damage"]
    brk = (bm > 0) | (dmg >= 0.999)
    seg, mode = PJ[CJ[brk]], bm[brk]
    xm, ym = seg[:, :, 0].mean(axis=1), seg[:, :, 1].mean(axis=1)
    r = np.hypot(xm - cx, ym - cy)
    L = np.linalg.norm(seg[:, 1] - seg[:, 0], axis=1)
    nt, ns = int((mode == 1).sum()), int((mode == 2).sum())

    fig, ax = plt.subplots(1, 2, figsize=(17.5, 8.8))
    for i, (z, lw, ttl) in enumerate(
            [(a.wide, 0.9, f"vue d'ensemble — EDZ p95 = {np.percentile(r, 95):.1f} m"),
             (a.tight, 1.7, "zoom sur la paroi")]):
        ax[i].add_collection(LineCollection(seg[mode == 2], colors=C_SHR,
                                            lw=lw, alpha=0.85))
        ax[i].add_collection(LineCollection(seg[mode == 1], colors=C_TEN,
                                            lw=lw, alpha=0.85))
        ax[i].plot(px, py, color="k", lw=2.0)
        ax[i].set_aspect("equal")
        ax[i].set_xlim(cx - z, cx + z)
        ax[i].set_ylim(cy - z, cy + z)
        ax[i].set_title(ttl, fontsize=12)
        ax[i].set_xlabel("x [m]")
        ax[i].grid(alpha=0.15)
    ax[0].set_ylabel("y [m]")
    ax[0].plot([], [], color=C_TEN, lw=3, label=f"traction  ({nt})")
    ax[0].plot([], [], color=C_SHR, lw=3, label=f"cisaillement  ({ns})")
    ax[0].legend(loc="lower left", fontsize=10, framealpha=0.95)

    fig.suptitle(f"{os.path.basename(run)} — trame {tag} : {nt + ns} fissures, "
                 f"{100 * ns / max(nt + ns, 1):.0f} % de cisaillement, "
                 f"longueur cumulee {L.sum():.0f} m", fontsize=13)
    fig.tight_layout()
    fig.savefig(out, dpi=190)
    print("ecrit :", out)


if __name__ == "__main__":
    main()
