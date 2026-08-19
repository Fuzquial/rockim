#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_yan_yang.py — deux planches sur les capacites ajoutees les 2026-08-18/19.
#
#   fig_dif_yang    : les deux lectures de l eq. 3 de Yang et al. (IJRMMS 191,
#                     2025) et l ATTRACTEUR que la version litterale cree dans
#                     un schema d insertion extrinseque. Le panneau du bas
#                     porte les populations de joints REELLEMENT inserees,
#                     mesurees sur les runs de non-regression.
#   fig_coupe_essais : les six essais de coupe PDC, ratio d injection et pic de
#                     force, contre la cible de Heilman et al. (ARMA 24-0238).
#
#   python tunnel_edz/fig_yan_yang.py
#
# Pas de text.usetex : les etiquettes sont en unicode direct (un escape LaTeX
# dans une chaine simple s afficherait tel quel — erreur payee le 2026-08-19).
# ---------------------------------------------------------------------------
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.labelsize": 10, "axes.titlesize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8.5,
})
HERE = os.path.dirname(os.path.abspath(__file__))
C_LIT, C_FIG, C_3D = "#C8342B", "#0B4F9E", "#1B8A3A"


def dif_t(e, n):
    e = np.asarray(e, dtype=float)
    v = np.where(e <= 5e-6, 1.0,
                 np.where(e > 1e2, 1.85,
                          0.95 + 0.41 * np.power(np.maximum(e, 1e-30), n)))
    return np.clip(v, 1.0, 1.85)


def dif_c(e):
    e = np.asarray(e, dtype=float)
    v = np.where(e <= 5e-6, 1.0,
                 np.where(e > 1e4, 1.84,
                          0.77 + 0.56 * np.power(np.maximum(e, 1e-30), 0.07)))
    return np.clip(v, 1.0, 1.84)


# populations MESUREES (min, mediane, max) du taux a l insertion, references
# figees dans verify_suite.py (2D : verify_fdem_tension ; 3D : fdem3d_tension)
POP = [
    ("3D, exposant 0,1707",              0.317,   1.5752,  57.78, C_3D),
    ("2D, exposant 0,1707, charge x3", 236.2,   315.75,  353.1,  C_FIG),
    ("2D, exposant 0,1707",             38.78,   40.216,  48.09, C_FIG),
    ("2D, exposant 0,07 (litteral)",    69.59,   99.356,  99.999, C_LIT),
]


def fig1():
    fig, (ax, bx) = plt.subplots(
        2, 1, figsize=(6.8, 5.4), sharex=True,
        gridspec_kw=dict(height_ratios=[2.6, 1.0], hspace=0.08))
    e = np.logspace(-7, 5, 4000)

    ax.semilogx(e, dif_t(e, 0.07), color=C_LIT, lw=1.9,
                label="DIF traction, exposant 0,07 — éq. 3 telle qu'imprimée")
    ax.semilogx(e, dif_t(e, 0.1707), color=C_FIG, lw=1.9,
                label="DIF traction, exposant 0,1707 — déduit de leur fig. 2b")
    ax.semilogx(e, dif_c(e), color="0.45", lw=1.2, ls="--",
                label="DIF compression, éq. 2 — exposant confirmé par leur fig. 2a")

    for x, lo, hi in ((5e-6, 1.0, float(dif_t(6e-6, 0.07))),
                      (1e2, float(dif_t(99.0, 0.07)), 1.85)):
        ax.annotate("", xy=(x, hi), xytext=(x, lo),
                    arrowprops=dict(arrowstyle="<->", color=C_LIT, lw=1.1))
        ax.text(x * 1.6, 0.5 * (lo + hi), "+%.0f %%" % (100 * (hi / lo - 1)),
                color=C_LIT, fontsize=9, va="center")

    for a in (ax, bx):
        a.axvline(1e2, color="0.72", lw=0.8, zorder=0)
    ax.set_ylabel("facteur d'amplification dynamique")
    ax.set_ylim(0.97, 1.95)
    ax.legend(loc="upper left", framealpha=0.95)
    ax.set_title("La forme imprimée de l'éq. 3 ne se raccorde à aucune de ses bornes,\n"
                 "et le saut de 1e2 s⁻¹ attire les insertions", fontsize=10)

    for k, (lab, lo, med, hi, c) in enumerate(POP):
        y = k
        bx.plot([lo, hi], [y, y], color=c, lw=6.5, solid_capstyle="butt",
                alpha=0.55)
        bx.plot([med], [y], "|", color=c, ms=15, mew=2.2)
        bx.text(lo * 0.72, y, lab, fontsize=8.2, color=c, va="center", ha="right")
    bx.set_ylim(-0.7, len(POP) - 0.3)
    bx.set_yticks([])
    bx.set_xlabel("vitesse de déformation à l'insertion  $\\dot\\varepsilon$  [s$^{-1}$]")
    bx.set_ylabel("joints insérés", fontsize=9)
    bx.text(1e2 * 1.3, len(POP) - 0.55, "discontinuité", fontsize=8, color="0.4")
    bx.set_xlim(1e-7, 1e5)
    bx.grid(axis="x", lw=0.4, color="0.9")
    bx.set_axisbelow(True)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "fig_dif_yang." + ext), dpi=200)
    plt.close(fig)
    print("ecrit : fig_dif_yang.pdf / .png")


RUNS = [
    ("v3 (réf.)",          408, 13.86, "0.40"),
    ("epfl",               110, 10.46, C_FIG),
    ("a2 = epfl + écrêt.",  57,  0.98, C_FIG),
    ("pot001",             439, 13.94, C_LIT),
    ("potxi",              446, 10.33, C_LIT),
    ("dt005",              520, 22.04, C_LIT),
]


def fig2():
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6))
    lab = [r[0] for r in RUNS]
    x = np.arange(len(RUNS))
    col = [r[3] for r in RUNS]

    ax = axes[0]
    ax.bar(x, [r[1] for r in RUNS], color=col, width=0.62)
    ax.axhline(1.0, color=C_3D, lw=1.3, ls="--")
    ax.text(-0.42, 1.25, "physique : 1", color=C_3D, fontsize=8.5,
            ha="left")
    ax.set_yscale("log")
    ax.set_ylim(0.6, 900)
    ax.set_ylabel("travail outil→solide / travail de corps rigide")
    ax.set_title("Énergie créée par le contact outil", fontsize=10)

    ax = axes[1]
    ax.bar(x, [r[2] for r in RUNS], color=col, width=0.62)
    ax.axhline(3.08, color=C_3D, lw=1.5)
    ax.text(-0.42, 4.1, "cible Heilman : 3,08", color=C_3D,
            fontsize=8.5, ha="left")
    ax.set_ylabel("pic de force outil  [MN m$^{-1}$]")
    ax.set_title("Pic de force contre la cible expérimentale", fontsize=10)

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(lab, rotation=26, ha="right", fontsize=8.2)
        ax.grid(axis="y", lw=0.4, color="0.9")
        ax.set_axisbelow(True)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "fig_coupe_essais." + ext), dpi=200)
    plt.close(fig)
    print("ecrit : fig_coupe_essais.pdf / .png")


if __name__ == "__main__":
    fig1()
    fig2()
