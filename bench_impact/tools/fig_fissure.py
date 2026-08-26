#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_fissure.py — OU la roche casse, en RENDU ELEMENTS.
#
#   python bench_impact/tools/fig_fissure.py out_imperial/fdem3d_joints_0006.vtu \
#          --stem bench_impact/fig_fissure --t-us 82.5
#
# Lit le VTU des JOINTS (ASCII, 6 noeuds : 0-1-2 cote A, 3-4-5 cote B) et
# dessine les FACETTES elles-memes, pas leurs centroides — regle maison :
# une fissure se montre en elements. Le maillage vivant reste en trait fin,
# ce qui donne l'echelle et montre que la fissure suit les faces des tetras.
#
#   (a) COUPE : les facettes dont le centroide tient dans une tranche mince
#       autour du plan y = y_axe, projetees sur (x, z). C'est une VRAIE coupe
#       de maillage, pas un nuage aplati.
#   (b) VUE EN PLAN : les memes facettes vues de dessus, sur une bande de
#       profondeur — la symetrie du facies s'y lit.
#   (c) MODE de rupture, en nombre.
#
# Trois couches, du fond vers l'avant :
#   maillage intact (gris tres clair) < zone de processus (D, colormap) <
#   facettes ROMPUES (pleines, couleur du mode, bord noir).
# ---------------------------------------------------------------------------
import argparse
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
    "font.size": 10,
})

BLEU, ROUGE, GRIS, VERT = "#1f4e79", "#b22222", "#9a9a9a", "#2e7d32"
ORANGE = "#c25e00"


def lire_vtu(path):
    txt = open(path, "r", encoding="utf-8", errors="replace").read()

    def bloc(motif):
        m = re.search(motif + r'[^>]*>(.*?)</DataArray>', txt, re.S)
        return m.group(1) if m else None

    pts = np.fromstring(bloc(r'<DataArray type="Float64" '
                             r'NumberOfComponents="3"'), sep=" ").reshape(-1, 3)
    champs = {}
    for nom in ("damage", "breakMode", "tBreak", "bonded"):
        b = bloc(r'<DataArray[^>]*Name="%s"' % nom)
        if b is not None:
            champs[nom] = np.fromstring(b, sep=" ")
    n = len(champs["damage"])
    conn = np.fromstring(bloc(r'<DataArray[^>]*Name="connectivity"'),
                         sep=" ", dtype=np.int64).reshape(n, -1)
    return pts, conn, champs


def poly(ax, tri, sel, **kw):
    """Ajoute les triangles selectionnes comme PolyCollection 2D."""
    if not np.any(sel):
        return None
    pc = PolyCollection(list(tri[sel]), **kw)
    ax.add_collection(pc)
    return pc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vtu")
    ap.add_argument("--stem", default="fig_fissure")
    ap.add_argument("--dmin", type=float, default=0.30)
    ap.add_argument("--t-us", type=float, default=None)
    ap.add_argument("--rmax", type=float, default=7.0, help="fenetre [mm]")
    ap.add_argument("--slab", type=float, default=1.2,
                    help="demi-epaisseur de la tranche de coupe [mm]")
    a = ap.parse_args()

    pts, conn, ch = lire_vtu(a.vtu)
    tri3 = pts[conn[:, :3]]                     # (n, 3, 3) triangle cote A
    cen = tri3.mean(axis=1)
    dmg = ch["damage"]
    brk = ch.get("tBreak", np.full(len(dmg), -1.0)) >= 0.0
    mode = ch.get("breakMode", np.zeros(len(dmg)))

    ref = cen[brk] if brk.sum() else cen[dmg > a.dmin]
    if len(ref) == 0:
        raise SystemExit("aucune facette rompue ni endommagee")
    x0, y0 = ref[:, 0].mean(), ref[:, 1].mean()
    zsurf = ref[:, 2].max()

    X = (cen[:, 0] - x0) * 1e3
    Y = (cen[:, 1] - y0) * 1e3
    Z = (zsurf - cen[:, 2]) * 1e3               # profondeur, > 0 vers le bas
    r = np.hypot(X, Y)
    proc = (dmg > a.dmin) & (~brk)

    nT = int((mode[brk] == 1).sum())
    nS = int((mode[brk] == 2).sum())
    nM = int(brk.sum()) - nT - nS

    fig = plt.figure(figsize=(13.6, 5.4))
    gs = fig.add_gridspec(1, 3, wspace=0.30, left=0.06, right=0.975,
                          top=0.79, bottom=0.14)
    ttl = "Ou la roche casse — replique Imperial College"
    if a.t_us is not None:
        ttl += "   ($t$ = %.1f $\\mu$s)" % a.t_us
    fig.suptitle(ttl, fontsize=13, y=0.965)
    fig.text(0.5, 0.885,
             "rendu ELEMENTS : chaque facette est une face de tetraedre. "
             "%d rompues, %d en zone de processus ($D > %.2f$), sur %d joints"
             % (brk.sum(), proc.sum(), a.dmin, len(dmg)),
             ha="center", fontsize=9.5, color="#444444", style="italic")

    # ================= (a) COUPE dans une tranche ==========================
    A = fig.add_subplot(gs[0, 0])
    slab = np.abs(Y) < a.slab
    fen = slab & (r < a.rmax * 1.6) & (Z < a.rmax * 1.6) & (Z > -1.0)
    triA = np.stack([(tri3[:, :, 0] - x0) * 1e3,
                     (zsurf - tri3[:, :, 2]) * 1e3], axis=2)

    poly(A, triA, fen & ~proc & ~brk, facecolors="none",
         edgecolors="#d8d8d8", linewidths=0.35, zorder=1)
    pc = poly(A, triA, fen & proc, array=dmg[fen & proc], cmap="YlOrBr",
              edgecolors="#b0906a", linewidths=0.3, zorder=2)
    if pc is not None:
        pc.set_clim(a.dmin, 1.0)
    for m, c in ((1, BLEU), (2, ROUGE)):
        poly(A, triA, fen & brk & (mode == m), facecolors=c,
             edgecolors="k", linewidths=0.6, zorder=4)
    poly(A, triA, fen & brk & (~np.isin(mode, (1, 2))), facecolors=ORANGE,
         edgecolors="k", linewidths=0.6, zorder=4)

    A.axhline(0, color="k", lw=0.9, zorder=5)
    A.set_xlim(-a.rmax, a.rmax)
    A.set_ylim(a.rmax, -0.6)
    A.set_aspect("equal")
    A.set_xlabel("$x$ depuis l'axe de l'insert  [mm]")
    A.set_ylabel("profondeur  [mm]")
    A.set_title("(a)  Coupe du maillage, tranche $|y| < %.1f$ mm" % a.slab,
                loc="left", fontsize=10.5)

    # ================= (b) VUE EN PLAN =====================================
    B = fig.add_subplot(gs[0, 1])
    bande = (Z > -0.6) & (Z < a.rmax)
    fenB = bande & (r < a.rmax * 1.5)
    triB = np.stack([(tri3[:, :, 0] - x0) * 1e3,
                     (tri3[:, :, 1] - y0) * 1e3], axis=2)

    poly(B, triB, fenB & ~proc & ~brk, facecolors="none",
         edgecolors="#dedede", linewidths=0.28, zorder=1)
    pc2 = poly(B, triB, fenB & proc, array=dmg[fenB & proc], cmap="YlOrBr",
               edgecolors="#b0906a", linewidths=0.25, zorder=2)
    if pc2 is not None:
        pc2.set_clim(a.dmin, 1.0)
    for m, c in ((1, BLEU), (2, ROUGE)):
        poly(B, triB, fenB & brk & (mode == m), facecolors=c,
             edgecolors="k", linewidths=0.5, zorder=4)
    poly(B, triB, fenB & brk & (~np.isin(mode, (1, 2))), facecolors=ORANGE,
         edgecolors="k", linewidths=0.5, zorder=4)

    B.set_xlim(-a.rmax, a.rmax)
    B.set_ylim(-a.rmax, a.rmax)
    B.set_aspect("equal")
    B.axhline(0, color=GRIS, lw=0.5, ls=":", zorder=5)
    B.axvline(0, color=GRIS, lw=0.5, ls=":", zorder=5)
    B.set_xlabel("$x$  [mm]")
    B.set_ylabel("$y$  [mm]")
    B.set_title("(b)  Vue de dessus, $0 < z < %.0f$ mm" % a.rmax,
                loc="left", fontsize=10.5)
    if pc2 is not None:
        cb = fig.colorbar(pc2, ax=B, fraction=0.045, pad=0.03)
        cb.set_label("endommagement $D$", fontsize=9)

    # ================= (c) MODE ============================================
    C = fig.add_subplot(gs[0, 2])
    vals = [nT, nS, nM]
    bars = C.bar(range(3), vals, color=[BLEU, ROUGE, ORANGE], width=0.6)
    for b, v in zip(bars, vals):
        C.annotate("%d" % v, (b.get_x() + b.get_width() / 2, v), ha="center",
                   va="bottom", fontsize=11, fontweight="bold")
    tot = max(1, sum(vals))
    C.set_xticks(range(3))
    C.set_xticklabels(["traction\n(mode I)", "cisaillement\n(mode II)",
                       "mixte"], fontsize=9)
    C.set_ylabel("facettes rompues")
    C.set_ylim(0, max(vals) * 1.32 + 1)
    C.grid(axis="y", lw=0.4, alpha=0.4)
    C.set_axisbelow(True)
    C.set_title("(c)  Mode de rupture : %.0f %% en traction"
                % (100.0 * nT / tot), loc="left", fontsize=10.5)
    C.annotate("Le run de reference ADAPTATIF donnait\n"
               "841 traction / 3 cisaillement (0,36 %) :\n"
               "rockim FEND, il ne broie pas.",
               (0.97, 0.93), xycoords="axes fraction", ha="right", va="top",
               fontsize=8, color="#555555", style="italic")

    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, ext), dpi=190)
    print("ecrit : %s.pdf et .png" % a.stem)
    print("  rompues %d (traction %d, cisaillement %d, mixte %d)"
          % (brk.sum(), nT, nS, nM))
    print("  zone de processus : %d facettes a D > %.2f" % (proc.sum(), a.dmin))
    print("  facettes dessinees : coupe %d, plan %d"
          % (fen.sum(), fenB.sum()))
    if brk.sum():
        print("  etendue radiale %.2f mm, profondeur max %.3f mm"
              % (r[brk].max(), Z[brk].max()))


if __name__ == "__main__":
    main()
