#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_enveloppe.py — l'ETENDUE de l'endommagement, sans le bruit du maillage.
#
#   python bench_impact/tools/fig_enveloppe.py out_imperial/fdem3d_joints_0006.vtu \
#          --stem bench_impact/fig_enveloppe --t-us 82.5 --pen 0.10
#
# fig_fissure.py dessine chaque facette ; a quelques milliers de triangles le
# maillage noie l'information. Ici on ne trace que les ENVELOPPES.
#
# METHODE — enveloppe par SECTEURS, et non enveloppe convexe : une enveloppe
# convexe lisserait des radiales en disque. On decoupe en secteurs angulaires
# (plan) ou en tranches de x (coupe) et on prend le 90e PERCENTILE dans
# chacun, PAS le maximum.
#
# DEUX PIEGES, tous deux constates le 2026-08-26 sur cette figure meme :
#   1. le rayon MAX par secteur est correle a +0,79 avec le NOMBRE de facettes
#      du secteur — l enveloppe se couvre alors de pointes qui ne sont que du
#      sur-echantillonnage. Le percentile ramene la correlation a +0,5 ;
#   2. le CENTRE pris sur les seules facettes rompues (7 a ce stade) etait
#      faux de 0,67 mm, ce qui suffisait a faire sortir le test de Rayleigh a
#      p = 0,000 (« anisotrope »). Recentre sur la zone de processus : p = 0,83,
#      donc ISOTROPE. Les deux pieges se combinaient en une belle etoile qui
#      n existait pas.
# D ou le test de Rayleigh imprime SUR la figure : l isotropie ne se juge pas
# a l oeil sur une enveloppe.
#
# Deux enveloppes superposees : zone de PROCESSUS (D > dmin, encore vivante)
# et facettes ROMPUES. L'ecart entre les deux est la marge d'avance de
# l'endommagement sur la rupture — la lecture la plus utile de la figure.
# ---------------------------------------------------------------------------
import argparse
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
    "font.size": 10,
})

BLEU, ROUGE, GRIS, VERT = "#1f4e79", "#b22222", "#9a9a9a", "#2e7d32"
ORANGE, SABLE = "#c25e00", "#e8a33d"
R_INS = 8.51            # rayon de la bille de carbure [mm]


def lire_vtu(path):
    txt = open(path, "r", encoding="utf-8", errors="replace").read()

    def bloc(motif):
        m = re.search(motif + r'[^>]*>(.*?)</DataArray>', txt, re.S)
        return m.group(1) if m else None

    pts = np.fromstring(bloc(r'<DataArray type="Float64" '
                             r'NumberOfComponents="3"'), sep=" ").reshape(-1, 3)
    ch = {}
    for nom in ("damage", "breakMode", "tBreak"):
        b = bloc(r'<DataArray[^>]*Name="%s"' % nom)
        if b is not None:
            ch[nom] = np.fromstring(b, sep=" ")
    n = len(ch["damage"])
    conn = np.fromstring(bloc(r'<DataArray[^>]*Name="connectivity"'),
                         sep=" ", dtype=np.int64).reshape(n, -1)
    return pts, conn, ch


def nsecteurs(n, par_sec=12, lo=6, hi=24):
    """Nombre de secteurs tel qu il reste ~`par_sec` facettes dans chacun.

    Le percentile ne protege du biais d echantillonnage que s il a de quoi
    travailler : a 4 facettes par secteur, le 90e percentile EST le maximum et
    l enveloppe redevient etoilee alors que le test de Rayleigh la declare
    isotrope. Le decoupage doit donc suivre la taille de l echantillon, pas
    etre fixe."""
    return int(np.clip(n // par_sec, lo, hi))


def env_polaire(x, y, nsec=36, q=90.0):
    """Rayon au q-ieme PERCENTILE par secteur -> courbe fermee (th, r).

    PAS le maximum. Mesure du 2026-08-26 : le rayon MAX par secteur est
    correle a +0,79 avec le NOMBRE de facettes du secteur — un secteur mieux
    echantillonne affiche mecaniquement un maximum plus grand, et l enveloppe
    se couvre de pointes qui n ont aucune realite physique. Au 90e percentile
    la correlation retombe a +0,5 et la courbe cesse de mentir."""
    if len(x) == 0:
        return None
    th = np.arctan2(y, x)
    r = np.hypot(x, y)
    bords = np.linspace(-np.pi, np.pi, nsec + 1)
    idx = np.clip(np.digitize(th, bords) - 1, 0, nsec - 1)
    rmax = np.zeros(nsec)
    for k in range(nsec):
        s = idx == k
        if s.any():
            rmax[k] = np.percentile(r[s], q)
    # secteurs vides : on interpole circulairement, sinon la courbe retombe a 0
    vide = rmax <= 0
    if vide.all():
        return None
    if vide.any():
        ok = ~vide
        ang = 0.5 * (bords[:-1] + bords[1:])
        rmax[vide] = np.interp(ang[vide], ang[ok], rmax[ok], period=2 * np.pi)
    ang = 0.5 * (bords[:-1] + bords[1:])
    return np.append(ang, ang[0]), np.append(rmax, rmax[0])


def env_coupe(x, z, nbin=26, xmax=7.0, q=90.0):
    """Profondeur au q-ieme PERCENTILE par tranche de x (meme raison)."""
    if len(x) == 0:
        return None
    bords = np.linspace(-xmax, xmax, nbin + 1)
    idx = np.digitize(x, bords) - 1
    xs, zs = [], []
    for k in range(nbin):
        s = idx == k
        if s.any():
            xs.append(0.5 * (bords[k] + bords[k + 1]))
            zs.append(np.percentile(z[s], q))
    return (np.array(xs), np.array(zs)) if xs else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vtu")
    ap.add_argument("--stem", default="fig_enveloppe")
    ap.add_argument("--dmin", type=float, default=0.30)
    ap.add_argument("--t-us", type=float, default=None)
    ap.add_argument("--pen", type=float, default=None,
                    help="enfoncement de l'insert [mm], pour son profil")
    ap.add_argument("--rmax", type=float, default=7.0)
    a = ap.parse_args()

    pts, conn, ch = lire_vtu(a.vtu)
    cen = pts[conn[:, :3]].mean(axis=1)
    dmg = ch["damage"]
    brk = ch.get("tBreak", np.full(len(dmg), -1.0)) >= 0.0
    mode = ch.get("breakMode", np.zeros(len(dmg)))

    # CENTRE : le barycentre de la zone de PROCESSUS, pondere par D — et non
    # celui des facettes rompues. Mesure du 2026-08-26 : avec 7 facettes
    # rompues le centre etait faux de 0,67 mm, ce qui suffisait a faire sortir
    # le test de Rayleigh a p = 0,000 (« anisotrope »). Recentre, il donne
    # p = 0,83 : le champ est ISOTROPE. Un centre mal pose fabrique une
    # anisotropie de toutes pieces.
    zone = dmg > a.dmin
    if not zone.any() and not brk.any():
        raise SystemExit("rien a montrer")
    ref = cen[zone] if zone.any() else cen[brk]
    w = dmg[zone] if zone.any() else np.ones(brk.sum())
    x0, y0 = (ref[:, :2] * w[:, None]).sum(0) / w.sum()
    zsurf = ref[:, 2].max()
    X = (cen[:, 0] - x0) * 1e3
    Y = (cen[:, 1] - y0) * 1e3
    Z = (zsurf - cen[:, 2]) * 1e3
    proc = (dmg > a.dmin) & (~brk)
    tout = proc | brk

    nT = int((mode[brk] == 1).sum())
    nS = int((mode[brk] == 2).sum())
    nM = int(brk.sum()) - nT - nS

    fig = plt.figure(figsize=(13.6, 5.2))
    gs = fig.add_gridspec(1, 3, wspace=0.30, left=0.06, right=0.975,
                          top=0.79, bottom=0.14)
    ttl = "Enveloppe de l'endommagement — replique Imperial College"
    if a.t_us is not None:
        ttl += "   ($t$ = %.1f $\\mu$s)" % a.t_us
    fig.suptitle(ttl, fontsize=13, y=0.965)
    # Test de Rayleigh sur la distribution ANGULAIRE de la zone de processus :
    # R proche de 0 = isotrope. Imprime parce qu une enveloppe se lit mal a
    # l oeil — celle-ci a deja paru etoilee alors qu elle ne l etait pas.
    thz = np.arctan2(Y[tout], X[tout])
    nz = max(1, len(thz))
    Ray = np.hypot(np.cos(thz).sum(), np.sin(thz).sum()) / nz
    pRay = float(np.exp(-nz * Ray * Ray))
    fig.text(0.5, 0.885,
             "enveloppe au 90e percentile sur %d secteurs (~12 facettes "
             "chacun ; le MAX serait biaise par l'echantillonnage)  |  "
             "Rayleigh $R$ = %.3f, $p$ = %.2f — facies %s"
             % (nsecteurs(int(tout.sum())), Ray, pRay,
                "ANISOTROPE" if pRay < 0.05 else "ISOTROPE"),
             ha="center", fontsize=9.5, color="#444444", style="italic")

    # ================= (a) COUPE ===========================================
    A = fig.add_subplot(gs[0, 0])
    ec = env_coupe(X[tout], Z[tout], nbin=nsecteurs(int(tout.sum())),
                   xmax=a.rmax)
    eb = env_coupe(X[brk], Z[brk], nbin=nsecteurs(int(brk.sum()), lo=4, hi=14),
                   xmax=a.rmax)
    if ec is not None:
        A.fill_between(ec[0], 0, ec[1], color=SABLE, alpha=0.35, zorder=2,
                       label="zone de processus ($D > %.2f$)" % a.dmin)
        A.plot(ec[0], ec[1], color="#a06a10", lw=1.6, zorder=3)
    if eb is not None:
        A.fill_between(eb[0], 0, eb[1], color=BLEU, alpha=0.45, zorder=4,
                       label="facettes rompues")
        A.plot(eb[0], eb[1], color=BLEU, lw=1.8, zorder=5)
    # profil de l'insert : bille de carbure, pointe a la profondeur `pen`
    if a.pen is not None:
        zc = a.pen - R_INS
        th = np.linspace(0, np.pi, 200)
        xi, zi = R_INS * np.cos(th), zc + R_INS * np.sin(th)
        keep = zi <= 0.02
        A.plot(xi[keep], zi[keep], color="k", lw=1.6, zorder=6)
        A.annotate("insert ($R$ = %.2f mm,\nenfonce de %.3f mm)"
                   % (R_INS, a.pen), (0, -0.25), ha="center", va="bottom",
                   fontsize=8, color="k")
    A.axhline(0, color="k", lw=1.0, zorder=6)
    A.set_xlim(-a.rmax, a.rmax)
    A.set_ylim(a.rmax, -1.1)
    A.set_aspect("equal")
    A.set_xlabel("$x$ depuis l'axe  [mm]")
    A.set_ylabel("profondeur  [mm]")
    A.legend(frameon=False, fontsize=8.5, loc="lower right")
    A.set_title("(a)  Coupe : jusqu'ou descend l'endommagement",
                loc="left", fontsize=10.5)

    # ================= (b) VUE EN PLAN =====================================
    B = fig.add_subplot(gs[0, 1])
    for sel, col, alp, lab in ((tout, SABLE, 0.35,
                                "zone de processus"),
                               (brk, BLEU, 0.45, "facettes rompues")):
        e = env_polaire(X[sel], Y[sel], nsec=nsecteurs(int(sel.sum())))
        if e is None:
            continue
        th, rr = e
        B.fill(rr * np.cos(th), rr * np.sin(th), color=col, alpha=alp,
               zorder=2 if col == SABLE else 4, label=lab)
        B.plot(rr * np.cos(th), rr * np.sin(th),
               color="#a06a10" if col == SABLE else BLEU, lw=1.7,
               zorder=3 if col == SABLE else 5)
    for rc in (2, 4, 6):
        B.add_patch(plt.Circle((0, 0), rc, fill=False, color=GRIS, lw=0.5,
                               ls=":", zorder=1))
        B.annotate("%d mm" % rc, (rc * 0.71, rc * 0.71), fontsize=7,
                   color=GRIS, ha="center", zorder=1)
    B.set_xlim(-a.rmax, a.rmax)
    B.set_ylim(-a.rmax, a.rmax)
    B.set_aspect("equal")
    B.set_xlabel("$x$  [mm]")
    B.set_ylabel("$y$  [mm]")
    B.legend(frameon=False, fontsize=8.5, loc="upper right")
    B.set_title("(b)  Vue de dessus : facies %s ($p$ = %.2f)"
                % ("ANISOTROPE" if pRay < 0.05 else "isotrope", pRay),
                loc="left", fontsize=10.5)

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

    rp = np.hypot(X[tout], Y[tout]).max() if tout.sum() else 0.0
    rb = np.hypot(X[brk], Y[brk]).max() if brk.sum() else 0.0
    fig.text(0.5, 0.028,
             "processus : rayon %.2f mm, profondeur %.2f mm  |  rupture : "
             "rayon %.2f mm, profondeur %.2f mm  |  l'endommagement precede "
             "la rupture de %.2f mm en rayon"
             % (rp, Z[tout].max() if tout.sum() else 0, rb,
                Z[brk].max() if brk.sum() else 0, rp - rb),
             ha="center", fontsize=8.5, color="#444444")

    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, ext), dpi=190)
    print("ecrit : %s.pdf et .png" % a.stem)
    print("  zone de processus : %d facettes, rayon %.2f mm, profondeur %.2f mm"
          % (proc.sum(), rp, Z[tout].max() if tout.sum() else 0))
    print("  rompues           : %d facettes, rayon %.2f mm, profondeur %.2f mm"
          % (brk.sum(), rb, Z[brk].max() if brk.sum() else 0))
    print("  modes             : traction %d, cisaillement %d, mixte %d"
          % (nT, nS, nM))


if __name__ == "__main__":
    main()
