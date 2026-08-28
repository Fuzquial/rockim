# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_ab_evolution.py — l'ECART A/B dans le temps : ou les deux runs se
# separent, et par quoi.
#
#   python bench_impact/tools/fig_ab_evolution.py out_imperial \
#          out_imperial_coulomb --stem bench_impact/figures/ab
#
# A = deck du 26/08 (plage de mode II divisee par la cohesion seule).
# B = A + jointShearRange = coulomb + jointFrictionScaled = 1, RIEN d'autre
# (diff verifie). Tout ecart entre les courbes ne peut donc venir que de la.
#
#   (a) RUPTURES : le compteur, avec l'instant de separation marque ;
#   (b) FROTTEMENT de contact : le poste qui doit monter si les joints morts
#       passent bien au contact frottant (cible publiee 32 J, hors d'atteinte
#       a ce stade — l'echelle le dit) ;
#   (c) MODE de rupture frame par frame : la signature, en % de cisaillement.
#
# La zone grisee marque t > fin de A : au-dela, B n'a plus de temoin.
# ---------------------------------------------------------------------------
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from imp_lib import history, read_vtu, joints_frame, frames_of, frame_times

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "font.size": 9, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42})

CA, CB = "#5b6b7a", "#b3242b"


def modes(run):
    """% de cisaillement et effectifs, frame par frame."""
    out = []
    tt = frame_times(run)
    for k in frames_of(run):
        p, c, f = read_vtu(joints_frame(run, k))
        sel = (f["tBreak"] >= 0.0) & (f["bonded"] < 0.5)
        n = int(sel.sum())
        if n == 0:
            out.append((tt[k] * 1e6, 0, 0, 0))
            continue
        bm = f["breakMode"][sel]
        out.append((tt[k] * 1e6, n, int((bm == 1).sum()),
                    int((bm == 2).sum())))
    return np.array(out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_a")
    ap.add_argument("run_b")
    ap.add_argument("--stem", default="bench_impact/figures/ab_evolution")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.stem) or ".", exist_ok=True)

    ha, hb = history(a.run_a), history(a.run_b)
    ta, tb = ha["t"] * 1e6, hb["t"] * 1e6
    t_fin_a = ta[-1]

    fig, ax = plt.subplots(1, 3, figsize=(11.2, 3.9))

    # --- (a) ruptures ------------------------------------------------------
    ax[0].plot(ta, ha["nBroken"], lw=1.3, color=CA, label="A — sans correctif")
    ax[0].plot(tb, hb["nBroken"], lw=1.3, color=CB, label="B — coulomb")
    # instant de separation : 1re rupture de B avant celle de A
    ia = np.argmax(ha["nBroken"] > 0)
    ib = np.argmax(hb["nBroken"] > 0)
    tsep = tb[ib]
    ax[0].axvline(tsep, ls=(0, (4, 3)), lw=0.8, color="0.45")
    ax[0].annotate("B casse dès %.1f " % tsep + r"$\mu$s" +
                   "\n(A : %.1f " % ta[ia] + r"$\mu$s)",
                   xy=(tsep, hb["nBroken"].max() * 0.42),
                   xytext=(tsep - 46, hb["nBroken"].max() * 0.62),
                   fontsize=7.8, color="0.3",
                   arrowprops=dict(arrowstyle="->", lw=0.6, color="0.5"))
    ax[0].set_ylabel("joints rompus")
    ax[0].set_title("(a) la rupture", fontsize=9.4, loc="left")
    ax[0].legend(fontsize=7.8, loc="upper left", frameon=False)

    # --- (b) frottement ----------------------------------------------------
    ax[1].plot(ta, -ha["eFric"], lw=1.3, color=CA)
    ax[1].plot(tb, -hb["eFric"], lw=1.3, color=CB)
    ax[1].set_ylabel("travail de frottement  [J]")
    ax[1].set_title("(b) le frottement de contact", fontsize=9.4, loc="left")
    ax[1].text(0.03, 0.94, "cible publiée : 32 J\n(l'échelle dit le chemin "
               "restant)", transform=ax[1].transAxes, fontsize=7.6,
               va="top", color="0.3", style="italic")

    # --- (c) mode, frame par frame ----------------------------------------
    ma, mb = modes(a.run_a), modes(a.run_b)
    for m, col, lab in [(ma, CA, "A"), (mb, CB, "B")]:
        ok = m[:, 1] > 0
        ax[2].plot(m[ok, 0], 100 * m[ok, 3] / m[ok, 1], "o-", ms=4, lw=1.3,
                   color=col, label=lab)
    ax[2].axhline(50, ls=(0, (2, 3)), lw=0.7, color="0.6")
    ax[2].set_ylim(-4, 104)
    ax[2].set_ylabel("ruptures en cisaillement  [%]")
    ax[2].set_title("(c) le MODE — la signature", fontsize=9.4, loc="left")
    ax[2].legend(fontsize=7.8, loc="center right", frameon=False)
    ax[2].text(0.03, 0.5, "A reste à ZÉRO :\nil fend, il ne broie pas",
               transform=ax[2].transAxes, fontsize=7.8, color=CA,
               va="center", style="italic")

    for x in ax:
        x.set_xlabel(r"temps  [$\mu$s]")
        x.axvspan(t_fin_a, max(tb[-1], t_fin_a) * 1.02, color="0.92",
                  zorder=0)
        x.spines[["top", "right"]].set_visible(False)
    ax[0].text(t_fin_a + 1.5, ax[0].get_ylim()[1] * 0.06,
               "A s'arrête ici", fontsize=7.4, color="0.45", rotation=90,
               va="bottom")

    fig.suptitle("Où les deux runs se séparent — A et B ne diffèrent que par "
                 "les deux clés du correctif", fontsize=10.6, y=1.0)
    fig.text(0.5, -0.03, "zone grisée : au-delà de la fin de A, B n'a plus de "
             "témoin — les écarts n'y sont plus des comparaisons",
             ha="center", fontsize=7.6, color="0.35", style="italic")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, ext), dpi=200, bbox_inches="tight")
    print("ecrit : %s.pdf et .png" % a.stem)
    print("  A : %d rompus a %.1f us, cisaillement final %.0f %%"
          % (ha["nBroken"][-1], ta[-1], 100 * ma[-1, 3] / max(ma[-1, 1], 1)))
    print("  B : %d rompus a %.1f us, cisaillement final %.0f %%"
          % (hb["nBroken"][-1], tb[-1], 100 * mb[-1, 3] / max(mb[-1, 1], 1)))
