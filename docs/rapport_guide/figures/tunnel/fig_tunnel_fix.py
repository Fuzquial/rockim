#!/usr/bin/env python3
"""Rejeu du tunnel EDZ (Wang et al. 2024) avec le contact corrige (b869dae),
maillage REDUIT tunnel_hs_red.msh (22 730 triangles), 4 fils, T = 0,25 s
(0,38 s pour l'intrinseque).

    python3 fig_tunnel_fix.py [DOSSIER_RUNS] [metrics_fix.json] [metrics_ref6oct.json]

Ligne du haut : joints rompus (rouge = cisaillement, bleu = traction) sur la
configuration deformee (deplacements amplifies x 5), elements en gris clair ;
adaptatif avant correction (A : cles au defaut), adaptatif apres (B :
contactCandidates = vertex, gcBirth = offset, potForceExact = true), pointe
k = 1,6 apres, intrinseque apres.
Ligne du bas : part de blocs mono-element (r <= 25 m, block_sizes.py) et
joints rompus pour toutes les variantes, 6 octobre (binaire f0209ef) / A / B.
Serie A de l'intrinseque : 1 fil (a 4 fils sur machine partagee, le calcul ne
tenait pas 40 min) ; elle redonne le 6 octobre au bit pres. Serie B de
l'intrinseque : non achevee (cout x10,5 mesure a T = 0,004 s, 1 fil).
Ecrit fig_tunnel_fix.pdf a cote du script (et un PNG de controle dans
DOSSIER_RUNS).
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fig_tunnel_rejeu as F  # style commun (Computer Modern, virgule) + analyse()  # noqa: E402

plt.rcParams.update({"font.size": 9})
SCR = "/tmp/claude-0/-home-user/95ba17a5-6321-5243-8157-c0b3e0616268/scratchpad/rejeu_tunnel_fix"
CX, CY = F.CX, F.CY
AMP = 5.0
vf = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace(".", ","))


def fr(x, n=1):
    return f"{x:.{n}f}".replace(".", ",")


def panel(a, run, lab, res):
    _, g = F.analyse(run)
    P0, P1, C = g["P0"], g["P1"], g["C"]
    Pa = P0 + AMP * (P1 - P0)
    cen = Pa[C].mean(axis=1)
    keep = (np.abs(cen[:, 0] - CX) < 23) & (np.abs(cen[:, 1] - CY) < 23)
    a.add_collection(PolyCollection(Pa[C[keep]] - [CX, CY], facecolors="0.88",
                                    edgecolors="none", rasterized=True))
    # segments rompus sur la configuration deformee amplifiee : memes points
    # (68 190 copies nodales) a la trame 0 et a la derniere trame des joints
    fj = F.frames(run, True)
    tj0, tj1 = open(fj[0]).read(), open(fj[-1]).read()
    Q0, Q1 = F.pts(tj0), F.pts(tj1)
    CJ = F.arr(tj1, "connectivity", int).reshape(-1, 2)
    brk = F.arr(tj1, "damage") >= 0.999
    Qa = Q0 + AMP * (Q1 - Q0)
    sa = Qa[CJ[brk]] - [CX, CY]
    col = np.where(g["mode"] == 2, "#cb181d", "#2171b5")
    a.add_collection(LineCollection(sa, colors=col, linewidths=0.4, rasterized=True))
    a.set_xlim(-22, 22)
    a.set_ylim(-22, 22)
    a.set_aspect("equal")
    a.set_title("%s\n%d rompus, $R_{95}$ = %s m\nblocs mono-élément %s %%" % (
        lab, res["broken"], fr(res["edz_radius_p95_m"]), fr(res["blocs_mono_pct"])),
        fontsize=7.5)
    a.xaxis.set_major_formatter(vf)
    a.yaxis.set_major_formatter(vf)
    a.set_xlabel("x [m]")


def main():
    scr = sys.argv[1] if len(sys.argv) > 1 else SCR
    mfix = json.load(open(sys.argv[2] if len(sys.argv) > 2 else os.path.join(scr, "metrics_fix.json")))
    mref = json.load(open(sys.argv[3] if len(sys.argv) > 3 else os.path.join(scr, "metrics_ref6oct.json")))
    by = {r["run"]: r for r in mfix}
    maps = [("out_red_adapt_A", "adaptatif, A (clés au défaut)"),
            ("out_red_adapt_B", "adaptatif, B (contact corrigé)"),
            ("out_red_tip16_B", "pointe $k$ = 1,6, B")]
    fig = plt.figure(figsize=(6.4, 6.0))
    gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 0.7], hspace=0.5, wspace=0.9)
    for k, (run, lab) in enumerate(maps):
        a = fig.add_subplot(gs[0, 2 * k:2 * k + 2])
        panel(a, os.path.join(scr, run), lab, by[run])
        if k:
            a.set_yticklabels([])
        else:
            a.set_ylabel("y [m]")
    # ligne du bas : parts de blocs mono-element et joints rompus
    var = [("adapt", "adaptatif"), ("tip16", "pointe 1,6"),
           ("d07", "amort. 0,7"), ("intr", "intrinsèque")]
    ref = {r["label"].replace("ref6oct_", ""): r for r in mref}
    fixA = {r["run"][len("out_red_"):-2]: r for r in mfix if r["run"].endswith("_A")}
    fixB = {r["run"][len("out_red_"):-2]: r for r in mfix if r["run"].endswith("_B")}
    series = [("6 octobre (f0209ef, 1 fil)", ref, "0.7"),
              ("A : b869dae, clés au défaut", fixA, "#6baed6"),
              ("B : contact corrigé", fixB, "#08519c")]
    x = np.arange(len(var))
    w = 0.27
    for j, key in enumerate(("blocs_mono_pct", "broken")):
        a = fig.add_subplot(gs[1, 3 * j:3 * j + 3])
        for s_, (nm, d, c) in enumerate(series):
            vals = [d[v][key] if v in d and d[v].get(key) is not None else np.nan
                    for v, _ in var]
            a.bar(x + (s_ - 1) * w, vals, w, color=c, label=nm)
        if "intr" not in fixB:
            yt = 4 if j == 0 else 300
            a.text(x[-1] + w, yt, "B non\nachevé", rotation=90, ha="center",
                   va="bottom", fontsize=6)
        a.set_xticks(x)
        a.set_xticklabels([l for _, l in var], fontsize=7, rotation=20)
        a.yaxis.set_major_formatter(vf)
        a.grid(axis="y", alpha=0.3)
        a.set_axisbelow(True)
        if j == 0:
            a.set_ylabel("blocs mono-élément ($r \\leq 25$ m) [%]")
            a.set_ylim(0, 100)
            h, l = a.get_legend_handles_labels()
        else:
            a.set_ylabel("joints rompus")
            a.yaxis.set_label_position("right")
            a.yaxis.tick_right()
    fig.legend(h, l, frameon=False, fontsize=7, ncol=3, loc="center",
               bbox_to_anchor=(0.5, 0.43))
    from matplotlib.lines import Line2D
    fig.legend([Line2D([], [], color="#cb181d"), Line2D([], [], color="#2171b5")],
               ["cisaillement", "traction"], frameon=False, fontsize=7, ncol=2,
               loc="center", bbox_to_anchor=(0.5, 0.985))
    fig.savefig(os.path.join(HERE, "fig_tunnel_fix.pdf"), dpi=220, bbox_inches="tight")
    fig.savefig(os.path.join(scr, "fig_tunnel_fix.png"), dpi=130, bbox_inches="tight")


if __name__ == "__main__":
    main()
