# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig8_criteres.py — l'equivalent de la FIG. 8 de Yang et al. (IJRMMS 206,
# 2026, 106660) sur NOS donnees : la planche qui DEFINIT les criteres de
# validation multi-criteres d'un impact a insert unique.
#
#   python bench_impact/tools/fig8_criteres.py out_imperial_coulomb 8 \
#          --stem bench_impact/figures/fig8_B
#
# ⚠️ CET ARTICLE N'EST PAS CELUI QUE REPLIQUE LE DECK. Le deck reproduit
# Yang et al., IJRMMS 191 (2025) 106125, sur le CALCAIRE St Anne. Le present
# article (IJRMMS 206, 2026) est la SUITE, sur le GRANITE Kuru Grey — autres
# proprietes (E 60 GPa, ft 10,98 MPa, c 29,84 MPa contre 57 / 7,0 / 18,8).
# Les valeurs de leur fig. 8 ne sont donc PAS des cibles pour ce run ; c'est
# la DEFINITION des criteres qui est reprise ici, pas les chiffres.
#
# Leur planche porte trois blocs. On les reproduit :
#   (a) TRANSFERT D'ENERGIE, contrainte : sigma(t) a la jauge -> max stress ;
#   (b) TRANSFERT D'ENERGIE, cinematique : enfoncement(t) -> vitesse
#       d'indentation, enfoncement maximal, vitesse de rebond ;
#   (c) FISSURATION : vue de dessus -> rayon de cratere et longueur des
#       fissures radiales (imp_lib.metrics, deja ecrit pour leur fig. 8).
# ---------------------------------------------------------------------------
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from imp_lib import history, read_vtu, joints_frame, broken, metrics, CX, CY

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "font.size": 9, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42})

BLEU, ROUGE, GRIS = "#1f4e79", "#b22222", "#8a8a8a"
R_BIT = 0.015


def pente(t, y, i0, i1):
    """Vitesse par regression sur un segment (m/s, positif = enfoncement)."""
    if i1 - i0 < 3:
        return np.nan
    return -np.polyfit(t[i0:i1], y[i0:i1], 1)[0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("frame", type=int)
    ap.add_argument("--stem", default="bench_impact/figures/fig8")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.stem) or ".", exist_ok=True)

    h = history(a.run)
    t = h["t"] * 1e6                                   # us
    sig = -h["szz_bit"] / 1e6                          # MPa, compression > 0
    pen = (h["z_insert"][0] - h["z_insert"]) * 1e3     # mm, enfoncement > 0

    fig = plt.figure(figsize=(9.8, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.35], hspace=0.42,
                          wspace=0.28)

    # ---- (a) contrainte a la jauge : leur "Max stress" --------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(t, sig, lw=1.1, color=BLEU)
    k = int(np.argmax(sig))
    ax.plot(t[k], sig[k], "o", ms=4, color=ROUGE, zorder=5)
    ax.annotate("Max stress\n%.0f MPa" % sig[k], xy=(t[k], sig[k]),
                xytext=(t[k] + 14, sig[k] * 0.86), fontsize=8.2,
                color=ROUGE, arrowprops=dict(arrowstyle="->", lw=0.7,
                                             color=ROUGE))
    ax.set_xlabel(r"temps  [$\mu$s]")
    ax.set_ylabel("contrainte à la jauge  [MPa]")
    ax.set_title("(a) transfert d'énergie : contrainte", fontsize=9.4,
                 loc="left")

    # ---- (b) enfoncement : indentation / max / rebond ---------------------
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.plot(t, -pen, lw=1.1, color=BLEU)
    kmax = int(np.argmax(pen))
    fini = kmax < len(pen) - 5                     # le rebond a-t-il commence ?
    # vitesse d'indentation : la pente sur la moitie montante
    i0 = int(np.argmax(pen > 0.02 * max(pen.max(), 1e-9)))
    v_ind = pente(h["t"], h["z_insert"], i0, max(i0 + 1, kmax))
    ax2.annotate("vitesse d'indentation\n%.2f m/s" % v_ind,
                 xy=(t[(i0 + kmax) // 2], -pen[(i0 + kmax) // 2]),
                 xytext=(t[i0] - 2, -pen.max() * 0.55), fontsize=8.2,
                 color="0.25", ha="left",
                 arrowprops=dict(arrowstyle="->", lw=0.7, color="0.45"))
    if fini:
        ax2.axhline(-pen.max(), ls=(0, (5, 4)), lw=0.8, color="0.35")
        ax2.annotate("enfoncement maximal\n%.3f mm" % pen.max(),
                     xy=(t[kmax], -pen.max()), xytext=(t[kmax] * 0.55,
                                                       -pen.max() * 1.08),
                     fontsize=8.2, color=ROUGE,
                     arrowprops=dict(arrowstyle="->", lw=0.7, color=ROUGE))
        v_reb = -pente(h["t"], h["z_insert"], kmax, len(t))
        ax2.annotate("vitesse de rebond\n%.2f m/s" % v_reb,
                     xy=(t[-1], -pen[-1]), xytext=(t[-1] * 0.72,
                                                   -pen.max() * 0.45),
                     fontsize=8.2, color="0.25",
                     arrowprops=dict(arrowstyle="->", lw=0.7, color="0.45"))
    else:
        ax2.text(0.97, 0.08, "le run n'a PAS atteint l'enfoncement maximal :\n"
                 "ni le maximum ni la vitesse de rebond\nne sont encore "
                 "mesurables", transform=ax2.transAxes, ha="right",
                 va="bottom", fontsize=8, color=ROUGE, style="italic",
                 bbox=dict(fc="white", ec=ROUGE, lw=0.6, alpha=0.92, pad=3))
    ax2.set_xlabel(r"temps  [$\mu$s]")
    ax2.set_ylabel("enfoncement de l'insert  [mm]")
    ax2.set_title("(b) transfert d'énergie : cinématique", fontsize=9.4,
                  loc="left")

    # ---- (c) fissuration vue de dessus : cratere et fissures radiales -----
    ax3 = fig.add_subplot(gs[1, :])
    p, c, f = read_vtu(joints_frame(a.run, a.frame))
    ctr, nrm, mode, P = broken(p, c, f)
    m = metrics(ctr)
    tri = [np.c_[(q[:, 0] - CX) * 1e3, (q[:, 1] - CY) * 1e3] for q in P]
    for mo, col, lab in [(1, BLEU, "traction (mode I)"),
                         (2, ROUGE, "cisaillement (mode II)")]:
        sel = [q for q, k in zip(tri, mode) if k == mo]
        if sel:
            ax3.add_collection(PolyCollection(sel, facecolors=col,
                                              edgecolors="none", alpha=0.75,
                                              label="%s — %d" % (lab,
                                                                 len(sel))))
    th = np.linspace(0, 2 * np.pi, 200)
    for rr, col, lab in [(m["crater"] * 1e3, "#0b7d3e", "rayon de cratère"),
                         (m["radial"] * 1e3, "#5b2d8e",
                          "longueur des fissures radiales")]:
        ax3.plot(rr * np.cos(th), rr * np.sin(th), ls=(0, (5, 3)), lw=1.1,
                 color=col, label="%s : %.2f mm" % (lab, rr))
    lim = m["radial"] * 1e3 * 1.5 + 1
    ax3.set_xlim(-lim, lim)
    ax3.set_ylim(-lim, lim)
    ax3.set_aspect("equal")
    ax3.set_xlabel(r"$x$  [mm]")
    ax3.set_ylabel(r"$y$  [mm]")
    ax3.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8,
               frameon=True, framealpha=1, edgecolor="0.75")
    ax3.set_title("(c) fissuration, vue de dessus — les deux longueurs de "
                  "leur planche", fontsize=9.4, loc="left")
    # Les deux rayons COINCIDENT tant qu'aucune fissure radiale n'a quitte le
    # noyau : la facette la plus lointaine est encore une facette de surface.
    # Le dire, plutot que de laisser croire a un cercle manquant.
    if abs(m["crater"] - m["radial"]) < 1e-9:
        ax3.text(1.02, 0.42, "les deux rayons COÏNCIDENT :\naucune fissure "
                 "radiale n'a encore\nquitté le noyau broyé — la facette\nla "
                 "plus lointaine est de surface.", transform=ax3.transAxes,
                 fontsize=7.8, va="top", color=ROUGE, style="italic")

    fig.suptitle("Les critères de validation multi-critères (forme de la "
                 "fig. 8 de Yang et al. 2026)", fontsize=11, y=0.985)
    fig.text(0.5, 0.005, "la fig. 8 citée est celle du GRANITE Kuru Grey "
             "(IJRMMS 206, 2026) ; ce run réplique le CALCAIRE St Anne "
             "(IJRMMS 191, 2025) — on en reprend les DÉFINITIONS, pas les "
             "valeurs", ha="center", fontsize=7.6, color=ROUGE,
             style="italic")
    fig.tight_layout(rect=(0, 0.028, 1, 0.965))
    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, ext), dpi=200, bbox_inches="tight")
    print("ecrit : %s.pdf et .png" % a.stem)
    print("  max stress        : %.1f MPa a t = %.1f us" % (sig[k], t[k]))
    print("  v_indentation     : %.2f m/s" % v_ind)
    print("  enfoncement       : %.3f mm  (maximum atteint : %s)"
          % (pen.max(), "OUI" if fini else "NON"))
    print("  rayon de cratere  : %.2f mm" % (m["crater"] * 1e3))
    print("  fissures radiales : %.2f mm" % (m["radial"] * 1e3))
    print("  profondeur        : %.2f mm" % (m["depth"] * 1e3))
