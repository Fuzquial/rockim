#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# plot_tunnel_fields.py — planche complète d'un run tunnel, utilisable SUR UN
# RUN EN COURS : les fissures sont lues dans le VTU de joints de la trame (et
# non dans fdem_final_joints.csv, qui n'existe qu'à la fin).
#
#   python tunnel_edz/plot_tunnel_fields.py out_tun_ref [--frame -1] [--zoom 26]
#
# Huit panneaux :
#   (a) fissures par mode (breakMode du solveur) + contour d'EDZ
#   (b) contrainte principale majeure sigma_1
#   (c) cisaillement maximal tau_max = (s1 - s3)/2
#   (d) module du déplacement
#   (e) DÉFORMATION de cisaillement maximale gamma_max = E1 - E2
#   (f) DÉFORMATION volumique E1 + E2 (dilatance = décompression du massif)
#   (g) cinétique de fissuration
#   (h) le massif entier
#
# Les déformations viennent du gradient de transformation F de chaque triangle
# (positions de la trame 0 contre celles de la trame courante), via le tenseur
# de Green-Lagrange E = (F^T F - I)/2 — invariant par rotation, donc valable
# même là où des blocs ont tourné.
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
sys.path.insert(0, HERE)
from make_unstructured_mesh import TUNNEL_HS  # noqa: E402
from plot_tunnel_mesh import profile_xy  # noqa: E402

C_TEN, C_SHR = "#2E9E4F", "#C8342B"


def complete(path):
    """Une trame en cours d'écriture n'a pas encore sa balise de fin."""
    try:
        with open(path, "rb") as f:
            f.seek(max(0, os.path.getsize(path) - 200))
            return b"</VTKFile>" in f.read()
    except OSError:
        return False


def read_vtu(path, arrays, ncell=3, dim=2):
    """Lit un VTU rockim. dim = 2 (defaut, retro-compatible) rend x,y ;
    dim = 3 rend x,y,z.

    LE PIEGE, paye le 2026-08-18. Ce lecteur a ete ecrit pour le 2D et
    tronquait les points a [:, :2] SANS LE DIRE. Utilise tel quel sur une
    sortie 3D, il rend des tetraedres APLATIS dans le plan x-y : le calcul
    barycentrique devient singulier, et un audit d'interpenetration qui
    attrapait LinAlgError en silence a rendu « zero penetration » — un
    resultat faussement rassurant, presente comme une mesure.

    D'ou dim, EXPLICITE : un appel 3D qui oublie dim = 3 est une erreur, pas
    une approximation.
    """
    if dim not in (2, 3):
        raise ValueError("dim doit valoir 2 ou 3")
    with open(path) as f:
        txt = f.read()
    P = np.fromstring(re.search(
        r'<Points>.*?<DataArray[^>]*>(.*?)</DataArray>', txt, re.S).group(1),
        sep=" ").reshape(-1, 3)[:, :dim]
    if ncell == 4 and dim == 2:
        raise ValueError(
            "ncell = 4 (tetraedres) avec dim = 2 : les points seraient "
            "aplatis dans le plan x-y et toute mesure geometrique serait "
            "fausse. Passer dim = 3.")
    C = np.fromstring(re.search(
        r'Name="connectivity"[^>]*>(.*?)</DataArray>', txt, re.S).group(1),
        sep=" ").astype(int).reshape(-1, ncell)
    out = {}
    for nm in arrays:
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % nm, txt, re.S)
        out[nm] = np.fromstring(m.group(1), sep=" ") if m else None
    return P, C, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_tun_ref")
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--zoom", type=float, default=26.0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(HERE, "..", a.run)

    el = sorted(f for f in glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
    jn = sorted(glob.glob(os.path.join(run, "fdem_joints_[0-9]*.vtu")))
    el = [f for f in el if complete(f)]
    jn = [f for f in jn if complete(f)]
    k = min(len(el), len(jn)) - 1 if a.frame < 0 else a.frame
    tag = os.path.basename(el[k])[5:9]
    out = a.out or os.path.join(HERE,
                                f"{os.path.basename(run)}_f{tag}_champs.png")
    print(f"trame {tag} (sur {len(el)} completes)")

    P0, C, _ = read_vtu(el[0], [])
    P, C, S = read_vtu(el[k], ["sigmaXX", "sigmaYY", "sigmaXY"])
    W, H = P0[:, 0].max(), P0[:, 1].max()
    cx, cy = 0.5 * W, 0.5 * H
    px, py, _ = profile_xy(cx, cy - 0.5 * TUNNEL_HS["height"])

    # ---- contraintes principales -----------------------------------------
    sxx, syy, sxy = S["sigmaXX"], S["sigmaYY"], S["sigmaXY"]
    tr, dif = 0.5 * (sxx + syy), 0.5 * (sxx - syy)
    rad = np.hypot(dif, sxy)
    s1, tmax = (tr + rad) / 1e6, rad / 1e6

    # ---- deformations : F par triangle, Green-Lagrange -------------------
    X = P0[C]                                   # (n, 3, 2) config initiale
    x = P[C]                                    # (n, 3, 2) config courante
    A = np.stack([X[:, 1] - X[:, 0], X[:, 2] - X[:, 0]], axis=2)   # 2x2
    B = np.stack([x[:, 1] - x[:, 0], x[:, 2] - x[:, 0]], axis=2)
    F = B @ np.linalg.inv(A)
    E = 0.5 * (np.transpose(F, (0, 2, 1)) @ F - np.eye(2))
    em, ed = 0.5 * (E[:, 0, 0] + E[:, 1, 1]), 0.5 * (E[:, 0, 0] - E[:, 1, 1])
    er = np.hypot(ed, E[:, 0, 1])
    gmax = 2.0 * er                              # E1 - E2
    evol = 2.0 * em                              # E1 + E2

    umag = np.linalg.norm(P - P0, axis=1)
    uel = umag[C].mean(axis=1)
    T0 = mtri.Triangulation(P0[:, 0], P0[:, 1], C)

    # ---- fissures depuis le VTU de joints de la trame --------------------
    PJ, CJ, SJ = read_vtu(jn[k], ["damage", "breakMode"], ncell=2)
    bm, dmg = SJ["breakMode"], SJ["damage"]
    brk = (bm > 0) | (dmg >= 0.999)
    seg = PJ[CJ[brk]]
    mode = bm[brk]
    xm, ym = seg[:, :, 0].mean(axis=1), seg[:, :, 1].mean(axis=1)
    r = np.hypot(xm - cx, ym - cy)
    redz = np.percentile(r, 95) if len(r) else 0.0
    Lc = np.linalg.norm(seg[:, 1] - seg[:, 0], axis=1).sum() if len(seg) else 0

    Z = a.zoom
    fig = plt.figure(figsize=(21, 9.6))

    def frame(ax, title, ylab=False):
        ax.plot(px, py, color="k", lw=1.1)
        ax.set_aspect("equal")
        ax.set_xlim(cx - Z, cx + Z)
        ax.set_ylim(cy - Z, cy + Z)
        ax.set_title(title, fontsize=10.5)
        ax.set_xlabel("x [m]")
        if ylab:
            ax.set_ylabel("y [m]")

    ax = fig.add_subplot(2, 4, 1)
    ax.add_collection(LineCollection(seg[mode == 2], colors=C_SHR, lw=0.5))
    ax.add_collection(LineCollection(seg[mode == 1], colors=C_TEN, lw=0.5))
    if redz:
        ax.add_patch(plt.Circle((cx, cy), redz, fill=False, ls="--",
                                color="#B8860B", lw=1.5))
    ns, nt = int((mode == 2).sum()), int((mode == 1).sum())
    frame(ax, f"(a) {ns + nt} fissures : {nt} traction, {ns} cisaillement "
              f"({100 * ns / max(ns + nt, 1):.0f} %)", True)
    ax.plot([], [], color=C_TEN, lw=2, label="traction")
    ax.plot([], [], color=C_SHR, lw=2, label="cisaillement")
    ax.plot([], [], ls="--", color="#B8860B", label=f"EDZ p95 = {redz:.1f} m")
    ax.legend(loc="lower left", fontsize=7.5, framealpha=0.9)

    for i, (fld, cmap, ttl, kw) in enumerate([
            (s1, "RdBu_r", r"(b) contrainte principale majeure $\sigma_1$ [MPa]",
             dict(vmin=-np.percentile(np.abs(s1), 99),
                  vmax=np.percentile(np.abs(s1), 99))),
            (tmax, "inferno", r"(c) cisaillement maximal $\tau_{max}$ [MPa]",
             dict(vmax=np.percentile(tmax, 99.5))),
            (uel, "turbo", f"(d) deplacement |u| [m] — max {umag.max():.3f} m",
             {})]):
        ax = fig.add_subplot(2, 4, i + 2)
        h = ax.tripcolor(T0, facecolors=fld, cmap=cmap, **kw)
        frame(ax, ttl)
        fig.colorbar(h, ax=ax, shrink=0.85)

    ax = fig.add_subplot(2, 4, 5)
    h = ax.tripcolor(T0, facecolors=100 * gmax, cmap="magma",
                     vmax=np.percentile(100 * gmax, 99.5))
    frame(ax, r"(e) DEFORMATION de cisaillement $\gamma_{max}$ [%]", True)
    fig.colorbar(h, ax=ax, shrink=0.85)

    ax = fig.add_subplot(2, 4, 6)
    v = np.percentile(np.abs(100 * evol), 99)
    h = ax.tripcolor(T0, facecolors=100 * evol, cmap="PuOr", vmin=-v, vmax=v)
    # PuOr : violet = valeurs HAUTES = trace(E) > 0 = DILATANCE (la zone
    # rompue foisonne), orange = compaction.
    frame(ax, r"(f) DEFORMATION volumique $\epsilon_v$ [%] (violet = dilatance)")
    fig.colorbar(h, ax=ax, shrink=0.85)

    ax = fig.add_subplot(2, 4, 7)
    d = np.genfromtxt(os.path.join(run, "history.csv"), delimiter=",",
                      names=True, invalid_raise=False)
    ax.plot(d["t"] * 1e3, d["nBroken"], color=C_SHR, lw=1.6)
    ax.set_xlabel("temps [ms]")
    ax.set_ylabel("joints rompus")
    ax.set_title("(g) cinetique de fissuration", fontsize=10.5)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(2, 4, 8)
    ax.add_collection(LineCollection(seg[mode == 2], colors=C_SHR, lw=0.35))
    ax.add_collection(LineCollection(seg[mode == 1], colors=C_TEN, lw=0.35))
    ax.plot(px, py, color="k", lw=0.9)
    ax.add_patch(plt.Rectangle((0, 0), W, H, fill=False, color="0.5", lw=0.8))
    ax.set_aspect("equal")
    ax.set_xlim(-2, W + 2)
    ax.set_ylim(-2, H + 2)
    ax.set_title("(h) massif entier (100 x 100 m)", fontsize=10.5)
    ax.set_xlabel("x [m]")

    fig.suptitle(f"{os.path.basename(run)} — trame {tag}, "
                 f"longueur cumulee de fissures {Lc:.0f} m", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("ecrit :", out)


if __name__ == "__main__":
    main()
