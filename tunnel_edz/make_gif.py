#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_gif.py — animation d'un run tunnel : fissures, contraintes, déformations.
#
#   python tunnel_edz/make_gif.py out_tun_ref_iso [--zoom 20] [--duree 700]
#                                 [--dpi 100] [--sortie fichier.gif]
#
# Six panneaux par image, tous sur la MEME échelle de couleur du début à la
# fin — c'est la règle qui rend une animation lisible : une échelle recalculée
# à chaque image transforme la moindre fluctuation en clignotement et rend
# toute lecture d'évolution impossible. Les bornes sont donc prises sur la
# DERNIÈRE trame (la plus chargée) et imposées à toutes les autres.
#
#   (a) fissures par mode (vert traction / rouge cisaillement)
#   (b) contrainte principale majeure sigma_1
#   (c) cisaillement maximal tau_max
#   (d) module du déplacement |u|
#   (e) déformation de cisaillement gamma_max
#   (f) déformation volumique eps_v (violet = dilatance)
#
# Fonctionne sur un run EN COURS : seules les trames complètes sont lues.
# Assemblage par PIL, comme figures_impact/gif_impacted.py (imageio n'est pas
# installé sur cette machine).
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
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "tools"))
sys.path.insert(0, HERE)
from make_unstructured_mesh import TUNNEL_HS  # noqa: E402
from plot_tunnel_mesh import profile_xy  # noqa: E402
from plot_tunnel_fields import read_vtu, complete  # noqa: E402

C_TEN, C_SHR = "#1B8A3A", "#C8342B"


def fields(P0, P, C, S):
    """Contraintes principales, déformations de Green-Lagrange, déplacement."""
    sxx, syy, sxy = S["sigmaXX"], S["sigmaYY"], S["sigmaXY"]
    tr, dif = 0.5 * (sxx + syy), 0.5 * (sxx - syy)
    rad = np.hypot(dif, sxy)
    s1, tmax = (tr + rad) / 1e6, rad / 1e6
    X, x = P0[C], P[C]
    A = np.stack([X[:, 1] - X[:, 0], X[:, 2] - X[:, 0]], axis=2)
    B = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0]], axis=2)
    F = B @ np.linalg.inv(A)
    E = 0.5 * (np.transpose(F, (0, 2, 1)) @ F - np.eye(2))
    em = 0.5 * (E[:, 0, 0] + E[:, 1, 1])
    ed = 0.5 * (E[:, 0, 0] - E[:, 1, 1])
    er = np.hypot(ed, E[:, 0, 1])
    umag = np.linalg.norm(P - P0, axis=1)
    return dict(s1=s1, tmax=tmax, gmax=200.0 * er, evol=200.0 * em,
                uel=umag[C].mean(axis=1), umax=umag.max())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_tun_ref_iso")
    ap.add_argument("--zoom", type=float, default=20.0)
    ap.add_argument("--duree", type=int, default=700, help="ms par image")
    ap.add_argument("--dpi", type=int, default=100)
    ap.add_argument("--sortie", default=None)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(HERE, "..", a.run)
    out = a.sortie or os.path.join(HERE, os.path.basename(run) + "_anim.gif")

    el = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
          if complete(f)]
    jn = [f for f in sorted(glob.glob(os.path.join(run, "fdem_joints_[0-9]*.vtu")))
          if complete(f)]
    n = min(len(el), len(jn))
    if n < 2:
        raise SystemExit("moins de deux trames completes : rien a animer")
    print(f"{n} trames completes")

    P0, C, _ = read_vtu(el[0], [])
    W, H = P0[:, 0].max(), P0[:, 1].max()
    cx, cy = 0.5 * W, 0.5 * H
    px, py, _ = profile_xy(cx, cy - 0.5 * TUNNEL_HS["height"])
    T0 = mtri.Triangulation(P0[:, 0], P0[:, 1], C)

    # ---- bornes de couleur FIGEES sur la derniere trame -------------------
    Pl, _, Sl = read_vtu(el[n - 1], ["sigmaXX", "sigmaYY", "sigmaXY"])
    L = fields(P0, Pl, C, Sl)
    lim = dict(s1=np.percentile(np.abs(L["s1"]), 99),
               tmax=np.percentile(L["tmax"], 99.5),
               u=L["umax"],
               gmax=np.percentile(L["gmax"], 99.5),
               evol=np.percentile(np.abs(L["evol"]), 99))
    print("echelles figees : sigma1 +-{s1:.2f} MPa, tau {tmax:.2f}, "
          "|u| {u:.3f} m, gamma {gmax:.3f} %, eps_v +-{evol:.3f} %"
          .format(**lim))

    try:
        meta = np.genfromtxt(os.path.join(run, "frames.csv"), delimiter=",",
                             names=True, invalid_raise=False)
        tof = {int(f): t for f, t in zip(np.atleast_1d(meta["frame"]),
                                         np.atleast_1d(meta["t"]))}
    except Exception:
        tof = {}

    imgs = []
    for k in range(n):
        P, _, S = read_vtu(el[k], ["sigmaXX", "sigmaYY", "sigmaXY"])
        F = fields(P0, P, C, S)
        PJ, CJ, SJ = read_vtu(jn[k], ["damage", "breakMode"], ncell=2)
        bm = SJ["breakMode"]
        brk = (bm > 0) | (SJ["damage"] >= 0.999)
        seg, mode = PJ[CJ[brk]], bm[brk]

        fig = plt.figure(figsize=(16.5, 9.4))
        panels = [
            ("(a) fissures", None, None, None),
            (r"(b) $\sigma_1$ [MPa]", F["s1"], "RdBu_r",
             dict(vmin=-lim["s1"], vmax=lim["s1"])),
            (r"(c) $\tau_{max}$ [MPa]", F["tmax"], "inferno",
             dict(vmin=0, vmax=lim["tmax"])),
            ("(d) |u| [m]", F["uel"], "turbo", dict(vmin=0, vmax=lim["u"])),
            (r"(e) $\gamma_{max}$ [%]", F["gmax"], "magma",
             dict(vmin=0, vmax=lim["gmax"])),
            (r"(f) $\epsilon_v$ [%] (violet = dilatance)", F["evol"], "PuOr",
             dict(vmin=-lim["evol"], vmax=lim["evol"])),
        ]
        for i, (ttl, fld, cmap, kw) in enumerate(panels):
            ax = fig.add_subplot(2, 3, i + 1)
            if fld is None:
                ax.add_collection(LineCollection(seg[mode == 2], colors=C_SHR,
                                                 lw=0.7))
                ax.add_collection(LineCollection(seg[mode == 1], colors=C_TEN,
                                                 lw=0.7))
                nt, ns = int((mode == 1).sum()), int((mode == 2).sum())
                ttl += f" : {nt + ns}"
                if nt + ns:
                    ttl += f" ({100 * ns / (nt + ns):.0f} % cisaillement)"
            else:
                h = ax.tripcolor(T0, facecolors=fld, cmap=cmap, **kw)
                fig.colorbar(h, ax=ax, shrink=0.82)
            ax.plot(px, py, color="k", lw=1.2)
            ax.set_aspect("equal")
            ax.set_xlim(cx - a.zoom, cx + a.zoom)
            ax.set_ylim(cy - a.zoom, cy + a.zoom)
            ax.set_title(ttl, fontsize=10.5)
            if i >= 3:
                ax.set_xlabel("x [m]")
            if i % 3 == 0:
                ax.set_ylabel("y [m]")
        t = tof.get(k, float("nan"))
        fig.suptitle(f"{os.path.basename(run)} — trame {k}/{n - 1}   "
                     f"t = {1e3 * t:.0f} ms   |u| max = {F['umax']:.3f} m",
                     fontsize=13)
        fig.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=a.dpi)
        plt.close(fig)
        buf.seek(0)
        imgs.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE))
        print(f"  trame {k} : {int(brk.sum())} fissures")

    # la derniere image reste 4x plus longtemps : on lit l'etat final
    dur = [a.duree] * (len(imgs) - 1) + [4 * a.duree]
    imgs[0].save(out, save_all=True, append_images=imgs[1:], duration=dur,
                 loop=0, optimize=True)
    mo = os.path.getsize(out) / 1e6
    print(f"ecrit : {out}  ({len(imgs)} images, {mo:.1f} Mo)")


if __name__ == "__main__":
    main()
