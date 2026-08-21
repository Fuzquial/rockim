#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_signe_hydro.py — LA PREUVE DU SIGNE.
#
#   python bench_abuaisha/tools/fig_signe_hydro.py
#
# Une cavite sous pression interne S'OUVRE. C'est le seul comportement
# admissible, et il ne demande aucun modele : le fluide pousse la roche a
# l'OPPOSE de la normale sortante du solide.
#
# Deux mesures independantes, faites sur des runs deja sur le disque :
#
#   (a) FORAGE — rayon moyen des 622 noeuds de paroi, en fonction du temps.
#       out_f7_aniso  (confiningPressure, signe de reference) : le rayon
#       REMONTE quand la pression monte.
#       out_hfp_aniso (module hydro)                          : il DESCEND.
#
#   (b) FISSURE DE PARKER — deplacement vertical des noeuds confondus au
#       centre de la discontinuite. Les deux chemins donnent la meme
#       AMPLITUDE (0,120 vs 0,121 mm) et des SIGNES OPPOSES : la fissure
#       s'ouvre d'un cote, s'interpenetre de l'autre.
#
# C'est aussi la lecon sur H3 : parker_compare.py mesurait max(y) - min(y),
# une valeur absolue. Le test etait AVEUGLE AU SIGNE et a valide une
# interpenetration comme une ouverture.
# ---------------------------------------------------------------------------
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tunnel_edz"))
from plot_tunnel_fields import read_vtu, complete  # noqa: E402

plt.rcParams.update({"font.size": 10.5, "figure.dpi": 110})
CX, CY, RB = 4.0, 4.0, 0.05
TFIN, NFRAME = 3.0e-3, 60


def wall_radius(run):
    fs = [f for f in sorted(glob.glob(os.path.join(ROOT, run, "fdem_[0-9]*.vtu")))
          if complete(f)]
    P0, _, _ = read_vtu(fs[0], [])
    r0 = np.hypot(P0[:, 0] - CX, P0[:, 1] - CY)
    w = np.where(np.abs(r0 - RB) < 2.5e-3)[0]
    t, dr = [], []
    for k, f in enumerate(fs):
        P, _, _ = read_vtu(f, [])
        r = np.hypot(P[w, 0] - CX, P[w, 1] - CY)
        t.append(k * TFIN / NFRAME * 1e3)
        dr.append((r.mean() - r0[w].mean()) * 1e6)
    return np.array(t), np.array(dr), len(w)


def parker_dy(run, half=0.75):
    fs = [f for f in sorted(glob.glob(os.path.join(ROOT, run, "fdem_[0-9]*.vtu")))
          if complete(f)]
    P0, _, _ = read_vtu(fs[0], [])
    P, _, _ = read_vtu(fs[-1], [])
    cy = 0.5 * (P0[:, 1].max() + P0[:, 1].min())
    cx = 0.5 * (P0[:, 0].max() + P0[:, 0].min())
    on = (np.abs(P0[:, 1] - cy) < 1e-9) & (np.abs(P0[:, 0] - cx) <= half + 1e-9)
    xs = np.round(P0[on, 0], 9)
    xv = np.unique(xs)
    xc = xv[np.argmin(np.abs(xv - cx))]
    g = np.where(np.abs(xs - xc) < 1e-12)[0]
    return (P[on, 1] - P0[on, 1])[g] * 1e3


def main():
    fig, (ax, axp) = plt.subplots(1, 2, figsize=(13.4, 5.2),
                                  gridspec_kw={"width_ratios": [1.5, 1.0]})

    for run, lab, col in (("out_f7_aniso", "pression imposée (signe de référence)",
                           "#1b8a3a"),
                          ("out_hfp_aniso", "module hydro", "#c8342b")):
        t, dr, nw = wall_radius(run)
        ax.plot(t, dr, "-", color=col, lw=2.0, label=lab)
    ax.axvspan(0, 0.2, color="0.85", alpha=0.6, zorder=0)
    ax.text(0.1, ax.get_ylim()[0], " excavation", fontsize=8.5, color="0.4",
            rotation=90, va="bottom")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("temps  [ms]")
    ax.set_ylabel("variation du rayon moyen de paroi  [$\\mu$m]")
    ax.set_title("(a) forage — la pression monte dans les deux cas\n"
                 "un trou sous pression interne doit S'OUVRIR")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower left", fontsize=9.5, framealpha=0.95)
    ax.annotate("il s'ouvre", (2.4, -6.5), fontsize=10, color="#1b8a3a")
    ax.annotate("il se ferme", (1.0, -22.0), fontsize=10, color="#c8342b")

    lab = ["confinement\n(référence)", "module hydro"]
    for i, run in enumerate(("out_parker_c", "out_parker_hydro_c")):
        dy = parker_dy(run)
        axp.bar(np.arange(len(dy)) + i * 0.4 - 0.2, dy, width=0.38,
                color=["#1b8a3a", "#c8342b"][i], label=lab[i])
    axp.axhline(0, color="k", lw=0.8)
    axp.set_xticks([])
    axp.set_xlabel("les 6 nœuds confondus au centre de la fissure")
    axp.set_ylabel("déplacement vertical  [mm]")
    axp.set_title("(b) fissure de Parker — même amplitude,\n"
                  "signes strictement opposés")
    axp.grid(alpha=0.3, axis="y")
    axp.legend(fontsize=9.5, framealpha=0.95)

    fig.suptitle("Le chargement hydro est appliqué avec le SIGNE OPPOSÉ à celui "
                 "du confinement", fontsize=12.5, y=0.99)
    fig.tight_layout(rect=[0, 0.06, 1, 0.93])
    fig.text(0.5, 0.012,
             "hydroForces() : « half = 0.5 * hydroP_ * L * thk_ * n » — le commentaire juste au-dessus dit "
             "« la pression pousse le solide à l'OPPOSÉ de la normale sortante ».\n"
             "confiningForces(), sur la même liste exterior_ (normale sortante garantie), écrit "
             "« half = -0.5 * p * L * thk_ * n ». Il manque le signe moins.",
             ha="center", fontsize=9, style="italic", color="0.3")
    out = os.path.join(HERE, "..", "hydro_preuve_du_signe.png")
    fig.savefig(out, dpi=155)
    fig.savefig(out.replace(".png", ".pdf"))
    print("ecrit :", out)


if __name__ == "__main__":
    main()
