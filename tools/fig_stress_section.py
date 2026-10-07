#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_stress_section.py — COUPE DE CONTRAINTE d'un impact fdem3d, avec le
# reseau de joints rompus par-dessus. Equivalent de la fig. 16 de Yang et al.
# (« Evolution of maximum principal stress ») superposee a leur fig. 13.
#
#   python tools/fig_stress_section.py out_dir [--frame k] [--field sigma1]
#          [--clim MIN MAX] [--lim 25] [--depth 25] [--stem results/fig/x]
#          [--smooth 0.4] [--grid 600] [--phase rock] [--plan z --z0 -2]
#
# --plan y (defaut) : coupe VERTICALE par l axe, repere (x, z sous la surface).
# --plan z --z0 D : coupe HORIZONTALE a la profondeur D mm sous la surface,
# repere (x, y) — la vue de dessus. Les deux utilisent la meme intersection
# exacte tetraedre-plan ; seule la composante testee change.
#
# --phase NOM : ne garde que les elements de cette phase (defaut `rock`). Sans
# ce filtre la coupe traverse aussi l INSERT en carbure, tres charge, dont le
# lisere rouge sombre ecrase l echelle et n a rien a voir avec la roche.
# `--phase toutes` restitue l ancien comportement.
#
# --smooth L (mm) : champ LISSE au lieu des polygones d elements. Les valeurs
# des sections (constantes par tetraedre) sont reportees sur une grille
# reguliere par moyenne ponderee gaussienne de portee L, puis tracees en aplat
# continu. Le maillage disparait, la localisation reste. L ~ la taille de
# maille (1,4 mm ici) lisse juste ce qu il faut ; L trop grand efface le
# gradient au bord du bulbe. --smooth 0 (defaut) garde les polygones exacts.
#
# REUTILISE les fonctions eprouvees de tools/fig_impact3d.py (lecture VTU,
# intersection EXACTE des tetraedres avec le plan y = YC, intersection des
# facettes de joint) : la coupe rendue est la section reelle des elements, pas
# un nuage de centroides interpole — c'etait le bug du 2026-09-03.
# Ce script n'appelle PAS son panneau de courbes, qui exige un rebond acheve
# (il plante sur un run encore en charge : IndexError dans state()).
#
# Champs disponibles dans le VTU des elements : sigma1 (contrainte principale
# maximale, traction > 0), vonMises, tauMax, pMean (pression moyenne, traction
# > 0 ; exportee seulement sous writeRuptureFields = true).
# ---------------------------------------------------------------------------
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fig_impact3d as fi                                        # noqa: E402
from fig_impact3d import rd, pts, cn, ar, slice_tets, slice_tris  # noqa: E402

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
})

DEF = {                      # (clim MPa, cmap, libelle)
    "sigma1":   ((-20.0, 20.0), "RdBu_r", r"$\sigma_1$ [MPa] (traction > 0)"),
    "pMean":    ((-200.0, 20.0), "RdBu_r", r"$p$ moyenne [MPa] (traction > 0)"),
    "vonMises": ((0.0, 140.0), "YlGnBu_r", "von Mises [MPa]"),
    "tauMax":   ((0.0, 80.0), "YlGnBu_r", r"$\tau_{max}$ [MPa]"),
}


def coupe(P, C, axe, v0, xc, yc, zsurf):
    """Intersection EXACTE des tetraedres avec le plan axe = v0 (metres).
    Meme algorithme que fig_impact3d.slice_tets, generalise a l axe z et
    rendant les deux coordonnees DU PLAN en mm : (x, z) sous `y`, (x, y)
    sous `z`. Les aretes qui changent de signe donnent les sommets du
    polygone, ordonnes par angle autour du centroide.
    """
    E6 = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    k = 1 if axe == "y" else 2
    d_all = P[C][:, :, k] - v0
    po, kp = [], []
    for i in np.where(~((d_all > 0).all(1) | (d_all < 0).all(1)))[0]:
        v = P[C[i]]
        d = v[:, k] - v0
        q = []
        for a_, b_ in E6:
            if d[a_] * d[b_] < 0:
                w = d[a_] / (d[a_] - d[b_])
                q.append(v[a_] + w * (v[b_] - v[a_]))
            elif d[a_] == 0.0:
                q.append(v[a_])
        if len(q) < 3:
            continue
        q = np.unique(np.round(np.array(q), 12), axis=0)
        if len(q) < 3:
            continue
        if axe == "y":
            uv = np.c_[(q[:, 0] - xc) * 1e3, (q[:, 2] - zsurf) * 1e3]
        else:
            uv = np.c_[(q[:, 0] - xc) * 1e3, (q[:, 1] - yc) * 1e3]
        ang = np.arctan2(uv[:, 1] - uv[:, 1].mean(), uv[:, 0] - uv[:, 0].mean())
        po.append(uv[np.argsort(ang)])
        kp.append(i)
    return po, np.array(kp, dtype=int)


def coupe_tris(P, C, axe, v0, xc, yc, zsurf):
    """Trace des facettes de joint dans le meme plan : segments, en mm."""
    k = 1 if axe == "y" else 2
    d_all = P[C][:, :, k] - v0
    sg = []
    for i in np.where(~((d_all > 0).all(1) | (d_all < 0).all(1)))[0]:
        v = P[C[i]]
        d = v[:, k] - v0
        q = []
        for a_, b_ in ((0, 1), (1, 2), (2, 0)):
            if d[a_] * d[b_] < 0:
                w = d[a_] / (d[a_] - d[b_])
                q.append(v[a_] + w * (v[b_] - v[a_]))
        if len(q) == 2:
            q = np.array(q)
            if axe == "y":
                sg.append(np.c_[(q[:, 0] - xc) * 1e3, (q[:, 2] - zsurf) * 1e3])
            else:
                sg.append(np.c_[(q[:, 0] - xc) * 1e3, (q[:, 1] - yc) * 1e3])
    return sg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--field", default="sigma1", choices=sorted(DEF))
    ap.add_argument("--clim", type=float, nargs=2, default=None)
    ap.add_argument("--lim", type=float, default=25.0)
    ap.add_argument("--depth", type=float, default=25.0)
    ap.add_argument("--xc", type=float, default=0.125)
    ap.add_argument("--yc", type=float, default=0.125)
    ap.add_argument("--zsurf", type=float, default=0.150)
    ap.add_argument("--plan", default="y", choices=("y", "z"),
                    help="y = coupe verticale par l axe ; z = vue de dessus")
    ap.add_argument("--z0", type=float, default=-2.0,
                    help="profondeur de la coupe horizontale [mm sous la surface]")
    ap.add_argument("--phase", default="rock",
                    help="phase a tracer : rock (defaut) | toutes | <nom>")
    ap.add_argument("--smooth", type=float, default=0.0,
                    help="portee du lissage gaussien [mm] ; 0 = polygones exacts")
    ap.add_argument("--grid", type=int, default=600,
                    help="points de la grille sur la largeur")
    ap.add_argument("--stem", default=None)
    ap.add_argument("--title", default=None)
    a = ap.parse_args()

    fs = sorted(glob.glob(os.path.join(a.run, "fdem3d_0*.vtu")))
    if not fs:
        raise SystemExit("aucun VTU d'elements dans " + a.run)
    k = a.frame if a.frame >= 0 else int(os.path.basename(fs[-1])[7:11])
    fe = os.path.join(a.run, "fdem3d_%04d.vtu" % k)
    fj = os.path.join(a.run, "fdem3d_joints_%04d.vtu" % k)
    clim, cmap, lab = DEF[a.field]
    if a.clim:
        clim = tuple(a.clim)
    stem = a.stem or os.path.join("results", "fig",
                                  "stress_%s_%s_%04d" % (os.path.basename(a.run.rstrip("/\\")), a.field, k))

    # slice_tets / slice_tris convertissent DEJA en mm par rapport a XC et H0,
    # les constantes du module fig_impact3d (banc de Yang : 75 mm, 100 mm). On
    # les regle sur ce deck, sinon la coupe sort hors cadre (piege du 14/09).
    fi.XC, fi.YC, fi.H0 = a.xc, a.yc, a.zsurf
    s = rd(fe)
    P, C = pts(s), cn(s, 4)
    val = ar(s, a.field)
    if val is None:
        raise SystemExit("champ %s absent du VTU (writeRuptureFields pour pMean)" % a.field)
    val = val / 1e6
    # ---- filtre de phase : la ROCHE seule par defaut ---------------------
    # Le VTU porte `phase` par element (indice de la liste `phases` du deck :
    # 0 = rock, 1 = steel, 2 = carbide ici). On restreint AVANT la coupe, donc
    # les elements du train ne sont ni traces ni pris dans le lissage.
    if a.phase != "toutes":
        ph = ar(s, "phase")
        if ph is None:
            raise SystemExit("champ `phase` absent du VTU : utiliser --phase toutes")
        idx = {"rock": 0, "steel": 1, "carbide": 2}
        want = idx.get(a.phase)
        if want is None:
            try:
                want = int(a.phase)
            except ValueError:
                raise SystemExit("--phase : rock | steel | carbide | toutes | <indice>")
        keep = np.where(ph == want)[0]
        if keep.size == 0:
            raise SystemExit("aucun element de la phase %s dans ce VTU" % a.phase)
        C, val = C[keep], val[keep]

    v0 = a.yc if a.plan == "y" else a.zsurf + a.z0 * 1e-3
    po, ei = coupe(P, C, a.plan, v0, a.xc, a.yc, a.zsurf)

    fig, A = plt.subplots(figsize=(9.2, 6.0))
    if a.smooth > 0.0:
        # ---- champ lisse -------------------------------------------------
        # Chaque section porte sa valeur (constante par tetraedre) en son
        # centroide, ponderee par son AIRE : un grand element ne pese pas
        # comme un petit. Moyenne gaussienne de portee --smooth sur une
        # grille reguliere (Nadaraya-Watson), puis aplat continu. Les points
        # de grille trop loin de toute section (hors du domaine coupe) sont
        # masques : on ne peint pas une valeur la ou il n y a pas de roche.
        cen = np.array([q.mean(axis=0) for q in po])
        are = np.array([0.5 * np.abs(np.dot(q[:, 0], np.roll(q[:, 1], -1))
                                     - np.dot(q[:, 1], np.roll(q[:, 0], -1)))
                        for q in po])
        vv = val[ei]
        nx = int(a.grid)
        nz = max(10, int(nx * (a.depth + 4.0) / (2.0 * a.lim))) if a.plan == "y" else nx
        gx = np.linspace(-a.lim, a.lim, nx)
        gz = (np.linspace(-a.depth, 4.0, nz) if a.plan == "y"
              else np.linspace(-a.lim, a.lim, nz))
        GX, GZ = np.meshgrid(gx, gz)
        num = np.zeros_like(GX)
        den = np.zeros_like(GX)
        dmin = np.full(GX.shape, 1e30)
        h2 = 2.0 * a.smooth ** 2
        for j in range(len(cen)):                 # boucle sur les sections
            d2 = (GX - cen[j, 0]) ** 2 + (GZ - cen[j, 1]) ** 2
            w = are[j] * np.exp(-d2 / h2)
            num += w * vv[j]
            den += w
            np.minimum(dmin, d2, out=dmin)
        F = np.where(den > 0, num / np.maximum(den, 1e-300), np.nan)
        F = np.ma.masked_where(~np.isfinite(F) | (dmin > (3.0 * a.smooth) ** 2), F)
        ext = ([-a.lim, a.lim, -a.depth, 4.0] if a.plan == "y"
               else [-a.lim, a.lim, -a.lim, a.lim])
        im = A.imshow(F, origin="lower", extent=ext,
                      cmap=cmap, vmin=clim[0], vmax=clim[1], interpolation="bilinear",
                      aspect="equal", zorder=1)
        cb = fig.colorbar(im, ax=A, fraction=0.046, pad=0.02)
    else:
        pc = PolyCollection(po, array=np.clip(val[ei], clim[0], clim[1]), cmap=cmap,
                            edgecolors="#8a939a", linewidths=0.15)
        pc.set_clim(*clim)
        A.add_collection(pc)
        cb = fig.colorbar(pc, ax=A, fraction=0.046, pad=0.02)
    cb.set_label(lab, fontsize=10)

    nseg = 0
    if os.path.isfile(fj):
        js = rd(fj)
        jp, jc = pts(js), cn(js, 3)
        tb = ar(js, "tBreak")
        msk = (tb >= 0.0) if tb is not None else (ar(js, "damage") >= 0.999)
        sg = coupe_tris(jp, jc[msk], a.plan, v0, a.xc, a.yc, a.zsurf)
        if sg:
            A.add_collection(LineCollection(sg, colors="#111111", linewidths=1.5,
                                            zorder=5))
            nseg = len(sg)

    if a.plan == "y":
        A.axhline(0, color="k", lw=0.8)
        A.set_xlim(-a.lim, a.lim)
        A.set_ylim(-a.depth, 4)
        A.set_ylabel("z sous la surface [mm]")
    else:
        A.set_xlim(-a.lim, a.lim)
        A.set_ylim(-a.lim, a.lim)
        A.set_ylabel("y [mm]")
    A.set_aspect("equal")
    A.set_xlabel("x [mm]")
    plan_lab = ("coupe verticale y = axe" if a.plan == "y"
                else "vue de dessus, z = %.1f mm" % a.z0)
    A.set_title(a.title or "%s — %s, trame %d ; %d segments de joints rompus"
                % (os.path.basename(a.run.rstrip("/\\")), plan_lab, k, nseg),
                fontsize=11, loc="left")
    if a.smooth > 0.0:
        A.text(0.99, 0.02, "champ lisse, portee %.1f mm ; phase %s" % (a.smooth, a.phase), ha="right",
               transform=A.transAxes, fontsize=8, color="#444444")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(stem + "." + ext, dpi=160)
    print("%s trame %d : %d elements coupes, %d segments de joints, %s dans [%.1f ; %.1f] MPa"
          % (a.run, k, len(po), nseg, a.field, np.nanmin(val[ei]), np.nanmax(val[ei])))
    print("ecrit : %s.pdf / .png" % stem)


if __name__ == "__main__":
    main()
