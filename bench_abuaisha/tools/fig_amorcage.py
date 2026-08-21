#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_amorcage.py — POURQUOI LE COUPLAGE HYDRO EST NECESSAIRE.
#
#   python bench_abuaisha/tools/fig_amorcage.py out_f7_aniso out_f7_iso
#
# Les deux runs sont charges par une PRESSION IMPOSEE sur les faces d'origine
# du forage (confineFaces = bore) — c'est-a-dire sans que la pression puisse
# entrer dans la fissure. Le commentaire du source de rockim l'annonce :
# « faces born from cracking receive nothing ».
#
# Resultat : la fissuration S'AMORCE au bon endroit, sous la contrainte
# orthoradiale, puis S'ARRETE. Il n'y a rien pour la propager.
#
# C'est le temoin negatif du benchmark : la figure 7 d'AbuAisha montre des
# fissures developpees sur ~20 cm, impossibles a obtenir sans le fluide.
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
sys.path.insert(0, os.path.join(HERE, "..", "..", "tunnel_edz"))
from plot_tunnel_fields import read_vtu, complete  # noqa: E402

plt.rcParams.update({"font.size": 10.5, "figure.dpi": 110})
CX, CY, RB = 4.0, 4.0, 0.05


def broken(run):
    js = sorted(f for f in glob.glob(os.path.join(run, "fdem_joints_0*.vtu"))
                if complete(f))
    # ncell=2 : les cellules d'un VTU de joints sont des SEGMENTS a 2 noeuds.
    # Sans cet argument le lecteur decoupe la connectivite par paquets de 3 et
    # fabrique des segments entre noeuds sans rapport — d'ou des « fissures »
    # de 60 cm sur un maillage a 3 mm. Piege paye le 2026-08-20.
    P, C, D = read_vtu(js[-1], ["damage"], ncell=2)
    d = np.asarray(D["damage"])
    seg = np.asarray(C)
    br = np.where(d >= 1.0)[0]
    dam = np.where((d > 0.05) & (d < 1.0))[0]
    return P, seg, br, dam, len(d)


def panel(ax, run, titre, half):
    P, seg, br, dam, ntot = broken(run)
    for idx, col, lw, z in ((dam, "#e8a87c", 0.8, 2), (br, "#8c2318", 1.6, 3)):
        if len(idx) == 0:
            continue
        L = np.stack([P[seg[idx, 0], :2], P[seg[idx, 1], :2]], axis=1)
        ax.add_collection(LineCollection(L - np.array([CX, CY]), colors=col,
                                         linewidths=lw, zorder=z))
    th = np.linspace(0, 2 * np.pi, 200)
    ax.fill(RB * np.cos(th), RB * np.sin(th), color="0.35", zorder=4)
    ax.set_xlim(-half, half)
    ax.set_ylim(-half, half)
    ax.set_aspect("equal")
    ax.set_xlabel("x  [m]")
    ax.set_title(f"{titre}\n{len(br)} joints rompus sur {ntot}", fontsize=10.5)
    ax.grid(alpha=0.25)
    return len(br), ntot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("aniso")
    ap.add_argument("iso")
    ap.add_argument("--half", type=float, default=0.30)
    ap.add_argument("--out", default="f7_temoin_sans_fluide")
    a = ap.parse_args()
    root = os.path.join(HERE, "..", "..")
    ra = a.aniso if os.path.isabs(a.aniso) else os.path.join(root, a.aniso)
    ri = a.iso if os.path.isabs(a.iso) else os.path.join(root, a.iso)

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.9))
    n1, N = panel(axes[0], ri, "état de contrainte ISOTROPE\n"
                  r"$\sigma'_h=\sigma'_H=-4{,}6$ MPa", a.half)
    n2, _ = panel(axes[1], ra, "état de contrainte ANISOTROPE\n"
                  r"$\sigma'_h=-4{,}6$, $\sigma'_H=-6{,}8$ MPa", a.half)
    axes[0].set_ylabel("y  [m]")

    from matplotlib.lines import Line2D
    axes[1].legend(handles=[
        Line2D([0], [0], color="#8c2318", lw=2, label="joint rompu (D = 1)"),
        Line2D([0], [0], color="#e8a87c", lw=1.4, label="endommagé (D > 0,05)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="0.35",
               markersize=9, label="forage $\\varnothing$ 0,1 m")],
        loc="lower right", fontsize=8.5, framealpha=0.95)

    fig.suptitle("Fracturation SANS couplage hydro — la pression ne peut pas "
                 "entrer dans la fissure\n"
                 "189 487 éléments, pression de paroi portée à 20 MPa "
                 "(seuil analytique : 12 MPa)", fontsize=11.5, y=0.99)
    fig.tight_layout(rect=[0, 0.06, 1, 0.93])
    fig.text(0.5, 0.015,
             "La fissuration s'amorce au bon endroit — à la paroi, sous la contrainte orthoradiale — "
             "puis S'ARRÊTE :\nsans fluide à l'intérieur, rien ne la propage. "
             "La figure 7 d'AbuAisha montre des fissures développées sur ~20 cm.",
             ha="center", fontsize=9, style="italic", color="0.3")
    out = os.path.join(HERE, "..", a.out + ".png")
    fig.savefig(out, dpi=155)
    fig.savefig(out.replace(".png", ".pdf"))
    print(f"ecrit : {out}")
    print(f"  isotrope   : {n1} joints rompus sur {N}")
    print(f"  anisotrope : {n2} joints rompus sur {N}")


if __name__ == "__main__":
    main()
