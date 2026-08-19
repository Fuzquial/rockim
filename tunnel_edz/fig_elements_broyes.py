#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_elements_broyes.py — carte des ELEMENTS declares broyes, a partir du
# nombre de leurs joints ROMPUS et du mode de rupture.
#
#   python tunnel_edz/fig_elements_broyes.py out_indent2d_yan [--n 2]
#   python tunnel_edz/fig_elements_broyes.py out_indent3d_yan [--n 3]
#
# CRITERE. Un element est declare broye quand au moins n de ses joints sont
# ROMPUS (D = 1). Le seuil qui a un sens est  tous sauf un  : n = 2 pour un
# triangle (3 aretes), n = 3 pour un tetraedre (4 faces). En dessous on marque
# tout ce qui borde une fissure ; au-dessus on ne trouve que les fragments
# deja detaches, que le compteur fragment donne deja.
#
# MODE. breakMode n est renseigne qu a l instant ou D atteint 1 : le critere
# porte donc sur les joints ROMPUS, pas sur les endommages. Un element peut
# cumuler les deux modes — on classe au MAJORITAIRE, et les egalites sont
# tracees a part plutot que tranchees arbitrairement.
#
# RESERVE A LIRE AVANT D INTERPRETER. Sur les runs avec law = mc, 98 % de
# l energie part en PLASTICITE DE VOLUME et 0,6 % dans les joints (mesure du
# 2026-08-19). Cette carte montre donc le RESEAU DE FISSURES, pas la zone
# broyee, qui est dans les elements et invisible aux joints. Ne pas l appeler
#  broyage  sans cette precision.
#
# Le lien joint -> element n existe dans aucune sortie : on le reconstruit par
# GEOMETRIE sur la trame 0 (non deformee), en appariant le milieu de chaque
# joint au milieu des aretes / centres des faces des elements. L ordre des
# joints etant stable d une trame a l autre, l endommagement final s y applique.
# ---------------------------------------------------------------------------
import argparse
import itertools
import os
import re
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm", "axes.labelsize": 10, "axes.titlesize": 10.5,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
})
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

C_TEN, C_SHE, C_MIX, C_ROCK = "#C0392B", "#E8B62C", "#7D5BA6", "#EDEAE6"


def grab(path, name, ncomp=1, dtype=float):
    d = open(path, "r", errors="ignore").read()
    if name == "points":
        m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", d, re.S)
    else:
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, d, re.S)
    if m is None:
        return None
    a = np.fromstring(m.group(1), sep=" ", dtype=dtype)
    return a.reshape(-1, ncomp) if ncomp > 1 else a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_indent2d_yan")
    ap.add_argument("--n", type=int, default=None, help="joints rompus requis")
    ap.add_argument("--half", type=float, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = os.path.join(ROOT, a.run)

    els = sorted(f for f in os.listdir(run)
                 if re.match(r"fdem3?d?_\d+\.vtu$", f) or re.match(r"fdem3d_\d+\.vtu$", f))
    jts = sorted(f for f in os.listdir(run) if "joints_" in f)
    d3 = "fdem3d" in jts[0]
    nv, nper = (3, 4) if d3 else (2, 3)      # noeuds par joint, joints par element

    # ---- trame 0 : geometrie de reference, pour l appariement -------------
    P0 = grab(os.path.join(run, els[0]), "points", 3)
    C0 = grab(os.path.join(run, els[0]), "connectivity", nper if d3 else 3,
              dtype=np.int64)
    JP0 = grab(os.path.join(run, jts[0]), "points", 3)
    JC0 = grab(os.path.join(run, jts[0]), "connectivity", nv, dtype=np.int64)

    # cle = centre de l arete (2D) / de la face (3D), arrondi
    def key(pts):
        return tuple(np.round(pts.mean(axis=0), 9))

    jkey = {}
    for k in range(len(JC0)):
        jkey.setdefault(key(JP0[JC0[k]]), []).append(k)

    faces = list(itertools.combinations(range(4), 3)) if d3 \
        else [(0, 1), (1, 2), (2, 0)]
    e_of_j = defaultdict(list)
    for e in range(len(C0)):
        for fc in faces:
            k = key(P0[C0[e], :][list(fc)])
            for j in jkey.get(k, []):
                e_of_j[j].append(e)

    # ---- trame finale : endommagement et mode ------------------------------
    dm = grab(os.path.join(run, jts[-1]), "damage")
    bm = grab(os.path.join(run, jts[-1]), "breakMode")
    nT = np.zeros(len(C0), int); nS = np.zeros(len(C0), int)
    for j in np.where(dm >= 1.0)[0]:
        for e in e_of_j.get(j, []):
            if bm[j] == 1: nT[e] += 1
            elif bm[j] == 2: nS[e] += 1
    tot = nT + nS

    n = a.n if a.n is not None else (3 if d3 else 2)
    sel = tot >= n
    mode = np.where(nT[sel] > nS[sel], 0, np.where(nS[sel] > nT[sel], 1, 2))
    print("  %s : critere n >= %d sur %d joints par element" % (a.run, n, nper))
    for k in range(1, nper + 1):
        print("     n >= %d : %6d elements" % (k, int((tot >= k).sum())))
    print("     retenus : %d  (traction %d, cisaillement %d, mixte %d)"
          % (sel.sum(), (mode == 0).sum(), (mode == 1).sum(), (mode == 2).sum()))
    if not sel.any():
        print("  aucun element ne satisfait le critere"); return

    # ---- trace : projection sur (x, z) en 3D, (x, y) en 2D -----------------
    PF = grab(os.path.join(run, els[-1]), "points", 3)
    T = PF[C0]
    cen = T.mean(axis=1)
    ax_h, ax_v = (0, 2) if d3 else (0, 1)
    x0 = cen[sel][:, 0].mean()
    ys = PF[:, ax_v].max()
    half = a.half if a.half is not None else (8.0 if d3 else 14.0)
    H = half * 1e-3

    fig, ax = plt.subplots(figsize=(7.4, 5.6))
    ax.set_facecolor(C_ROCK)
    idx = np.where(sel)[0]
    if d3:   # coupe : ne garder que la tranche centrale
        keep = np.abs(cen[idx, 1] - cen[sel][:, 1].mean()) < 1.5e-3
        idx, mode = idx[keep], mode[keep]
    for mv, col, lab in ((0, C_TEN, "traction"), (1, C_SHE, "cisaillement"),
                         (2, C_MIX, "mixte")):
        ii = idx[mode == mv]
        if not len(ii):
            continue
        pol = np.stack([(T[ii][:, :, ax_h] - x0) * 1e3,
                        (T[ii][:, :, ax_v] - ys) * 1e3], axis=-1)
        ax.add_collection(PolyCollection(pol, facecolors=col, edgecolors="none",
                                         alpha=0.95, label="%s (%d)" % (lab, len(ii))))
    ax.set_xlim(-half, half); ax.set_ylim(-2 * half, 0.4 * half)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("z − z$_{surface}$ [mm]")
    ax.legend(loc="lower right", framealpha=0.95)
    ax.set_title("Éléments à ≥ %d joints rompus — %s\n"
                 "réseau de fissures, PAS la zone broyée (98 %% de l'énergie "
                 "part en plasticité de volume)" % (n, a.run), fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    out = a.out or ("fig_broyes_" + a.run.replace("out_", ""))
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, out + "." + e), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  ecrit :", out)


if __name__ == "__main__":
    main()
