#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_gif_insertion.py — animation de l'INSERTION ADAPTATIVE : on suit, trame
# par trame, les joints qui naissent et ceux qui cassent.
#
#   python tunnel_edz/make_gif_insertion.py out_tun_ref_iso [--zoom 20]
#
# Trois etats par joint, lus dans le VTU de joints de chaque trame :
#   bonded = 1              -> l'arete n'existe pas encore comme joint : c'est
#                              une liaison nodale rigide. Non tracee.
#   bonded = 0, intact      -> joint INSERE et vivant : il porte de la
#                              contrainte et s'endommage. C'est la ZONE DE
#                              PROCESSUS. Trace en ambre.
#   breakMode > 0           -> joint ROMPU : vert (traction) ou rouge
#                              (cisaillement).
#
# L'interet : la couronne ambre en avant du front de rupture est exactement ce
# que l'insertion adaptative gere differemment d'un schema intrinseque, ou
# TOUS les joints existent depuis le debut. Sa largeur se compare a la longueur
# de zone cohesive l_cz = E GfI / ft^2 = 0,556 m pour ce materiau.
#
# Le panneau de droite suit deux fronts : le rayon p95 des joints inseres et
# celui des joints rompus. Leur ECART est la largeur de la zone de processus.
# ---------------------------------------------------------------------------
import argparse
import glob
import io
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
sys.path.insert(0, HERE)
from make_unstructured_mesh import TUNNEL_HS  # noqa: E402
from plot_tunnel_mesh import profile_xy  # noqa: E402
from plot_tunnel_fields import read_vtu, complete  # noqa: E402

C_INS, C_TEN, C_SHR = "#E8930C", "#1B8A3A", "#C8342B"
LCZ = 10e9 * 20.0 / (0.6e6 ** 2)          # 0,556 m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_tun_ref_iso")
    ap.add_argument("--zoom", type=float, default=20.0)
    ap.add_argument("--duree", type=int, default=650)
    ap.add_argument("--dpi", type=int, default=100)
    ap.add_argument("--sortie", default=None)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(HERE, "..", a.run)
    out = a.sortie or os.path.join(HERE,
                                   os.path.basename(run) + "_insertion.gif")

    el = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
          if complete(f)]
    jn = [f for f in sorted(glob.glob(os.path.join(run, "fdem_joints_[0-9]*.vtu")))
          if complete(f)]
    n = min(len(el), len(jn))
    P0, _, _ = read_vtu(el[0], [])
    W, H = P0[:, 0].max(), P0[:, 1].max()
    cx, cy = 0.5 * W, 0.5 * H
    px, py, _ = profile_xy(cx, cy - 0.5 * TUNNEL_HS["height"])

    try:
        meta = np.genfromtxt(os.path.join(run, "frames.csv"), delimiter=",",
                             names=True, invalid_raise=False)
        tof = {int(f): t for f, t in zip(np.atleast_1d(meta["frame"]),
                                         np.atleast_1d(meta["t"]))}
    except Exception:
        tof = {}

    # ---- premiere passe : compteurs et fronts, pour figer les axes ---------
    hist = []
    data = []
    for k in range(n):
        PJ, CJ, SJ = read_vtu(jn[k], ["bonded", "breakMode"], ncell=2)
        bd, bm = SJ["bonded"], SJ["breakMode"]
        seg = PJ[CJ]
        r = np.hypot(seg[:, :, 0].mean(axis=1) - cx,
                     seg[:, :, 1].mean(axis=1) - cy)
        ins = (bd < 0.5) & (bm == 0)          # insere, pas encore rompu
        brk = bm > 0
        f_ins = np.percentile(r[bd < 0.5], 95) if (bd < 0.5).any() else 0.0
        f_brk = np.percentile(r[brk], 95) if brk.any() else 0.0
        hist.append((tof.get(k, k) * 1e3, int((bd < 0.5).sum()), int(brk.sum()),
                     int(ins.sum()), f_ins, f_brk))
        data.append((seg, ins, brk, bm))
        print(f"  trame {k}: {(bd < 0.5).sum():6d} inseres, {brk.sum():6d} "
              f"rompus, {ins.sum():5d} actifs, fronts {f_ins:5.1f}/{f_brk:5.1f} m")
    hist = np.array(hist)

    imgs = []
    for k in range(n):
        seg, ins, brk, bm = data[k]
        fig = plt.figure(figsize=(15.5, 7.6))

        ax = fig.add_subplot(1, 2, 1)
        ax.add_collection(LineCollection(seg[ins], colors=C_INS, lw=0.8,
                                         alpha=0.9))
        ax.add_collection(LineCollection(seg[brk & (bm == 2)], colors=C_SHR,
                                         lw=0.6, alpha=0.75))
        ax.add_collection(LineCollection(seg[brk & (bm == 1)], colors=C_TEN,
                                         lw=0.6, alpha=0.75))
        ax.plot(px, py, color="k", lw=1.6)
        ax.set_aspect("equal")
        ax.set_xlim(cx - a.zoom, cx + a.zoom)
        ax.set_ylim(cy - a.zoom, cy + a.zoom)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_title(f"insertion adaptative — t = {hist[k, 0]:.0f} ms",
                     fontsize=11.5)
        ax.plot([], [], color=C_INS, lw=3,
                label=f"insere, vivant  ({int(hist[k, 3])})")
        ax.plot([], [], color=C_TEN, lw=3, label="rompu traction")
        ax.plot([], [], color=C_SHR, lw=3, label="rompu cisaillement")
        ax.legend(loc="lower left", fontsize=8.5, framealpha=0.95)

        ax = fig.add_subplot(2, 2, 2)
        ax.plot(hist[:, 0], hist[:, 1], "-", color=C_INS, lw=2, label="inseres")
        ax.plot(hist[:, 0], hist[:, 2], "-", color=C_SHR, lw=2, label="rompus")
        ax.plot(hist[:, 0], hist[:, 3], "-", color="#6A3D9A", lw=1.6,
                label="vivants (zone de processus)")
        ax.axvline(hist[k, 0], color="0.3", lw=1.2)
        ax.set_ylabel("nombre de joints")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        ax.set_title("population de joints", fontsize=11)

        ax = fig.add_subplot(2, 2, 4)
        ax.plot(hist[:, 0], hist[:, 4], "-", color=C_INS, lw=2,
                label="front d'insertion (p95)")
        ax.plot(hist[:, 0], hist[:, 5], "-", color=C_SHR, lw=2,
                label="front de rupture (p95)")
        ax.fill_between(hist[:, 0], hist[:, 5], hist[:, 4], color=C_INS,
                        alpha=0.18)
        ax.axvline(hist[k, 0], color="0.3", lw=1.2)
        d = hist[k, 4] - hist[k, 5]
        ax.set_xlabel("temps [ms]")
        ax.set_ylabel("rayon [m]")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        ax.set_title(f"fronts — ecart = {d:.2f} m "
                     f"({d / LCZ:.1f} x l_cz)", fontsize=11)

        fig.suptitle(f"{os.path.basename(run)} — naissance et rupture des "
                     f"joints cohesifs", fontsize=13)
        fig.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=a.dpi)
        plt.close(fig)
        buf.seek(0)
        imgs.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE))

    dur = [a.duree] * (len(imgs) - 1) + [4 * a.duree]
    imgs[0].save(out, save_all=True, append_images=imgs[1:], duration=dur,
                 loop=0, optimize=True)
    print(f"ecrit : {out}  ({len(imgs)} images, "
          f"{os.path.getsize(out) / 1e6:.1f} Mo)")


if __name__ == "__main__":
    main()
