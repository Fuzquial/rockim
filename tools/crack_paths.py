#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# crack_paths.py — FISSURES CONNECTEES et CRATERE d'un impact fdem3d
# (campagne de correction du 13/09/2026, tache T1 ; motif : DIAGNOSTIC
# independant du 12/09 §5 — « le rayon maximal des centroides des facettes
# rompues n'est pas la longueur d'une fissure radiale identifiee »).
#
#   python tools/crack_paths.py <run> [--frame k] [--rcore 6] [--skin 1]
#                               [--radial 0.35] [--nz 0.6] [--sec 1.5]
#                               [--current] [--check] [--stem results/fig/x]
#                               [--lim 0] [--depth 0]   (0 = cadrage ajuste)
#   python tools/crack_paths.py --selftest          (mini-test falsifiant)
#
# Ce que l'outil calcule, dans l'ordre :
#   1. facettes ROMPUES = tBreak >= 0 (imp_lib.broken_mask, jamais `damage`) ;
#   2. les noeuds du VTU sont DUPLIQUES par element (4 e + k) : ils sont
#      unifies par leurs coordonnees INITIALES (trame 0), puis deux facettes
#      rompues sont ADJACENTES si elles partagent une arete (deux sommets
#      unifies). Composantes connexes de ce graphe (scipy.csgraph) ;
#   3. NOYAU = la plus grande composante dont un sommet est a r < rcore
#      (6 mm) de l'axe ; les autres composantes qui touchent r < rcore sont
#      « centrales secondaires », le reste est PERIPHERIQUE ;
#   4. BRAS du noyau = sous-composantes (meme adjacence) des facettes du
#      noyau dont le centroide est a r > rcore ; pour chaque bras et chaque
#      composante : nombre de facettes, r_min / r_max (sur les SOMMETS),
#      longueur radiale r_max − r_min, profondeur, orientation moyenne
#      ponderee par l'aire (|n.e_r| et |n_z|), azimut ;
#      orientation : RADIALE si |n.e_r| < 0,35 ET |n_z| < 0,6 (une facette
#      horizontale a aussi |n.e_r| = 0 : sans le second test elle passerait
#      pour radiale) ; HORIZONTALE si |n_z| >= 0,6 ; sinon CONIQUE /
#      circonferentielle ;
#   5. « longueur de fissure radiale » au sens de Yang = distance du centre
#      a la POINTE (sommet le plus eloigne) du plus long bras RADIAL connecte
#      au noyau ; 0 s'il n'y a aucun bras radial ;
#   6. rayon de CRATERE = frontiere exterieure des facettes de SURFACE du
#      noyau (centroide a moins de `skin` = 1 mm sous la surface) : r max
#      des sommets, et par 12 secteurs (moyenne des r max, couverture).
# Geometrie de REFERENCE (trame 0) par defaut : une fissure est une surface
# materielle, et des fragments deplaces gonflent le rayon (DIAGNOSTIC §5) ;
# `--current` prend les positions courantes. Les deux r max sont imprimes.
# Sortie : tableau stdout, <stem>_components.csv, <stem>_summary.csv, figure
# <stem>.pdf (vectoriel, Computer Modern) + .png : (a) vue de dessus, une
# couleur par composante / bras ; (b) coupe verticale dans le plan de l'axe
# et de la pointe du plus long bras radial (et non a y = 0).
# Post-traitement PUR : aucune cle, aucun effet sur le solveur.
# ---------------------------------------------------------------------------
import argparse
import csv
import io
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bench_impact", "tools"))
from imp_lib import (CX, CY, Z_SURF, broken_mask, frame_times,  # noqa: E402
                     frames_of, joints_frame, read_vtu)

Q = 1e9                       # quantification des coordonnees (1 nm) pour
                              # l'unification des noeuds dupliques


# ---------------------------------------------------------------------------
# geometrie et graphe
# ---------------------------------------------------------------------------
def unify_vertices(p0):
    """Identifiant de sommet unique par coordonnees initiales (n_noeuds,)."""
    q = np.round(p0 * Q).astype(np.int64)
    _, inv = np.unique(q, axis=0, return_inverse=True)
    return inv.ravel()


def facet_adjacency_components(tri):
    """Composantes connexes des facettes (tri : (n,3) sommets unifies) par
    ARETE partagee. Retourne (n_comp, labels)."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    n = len(tri)
    if n == 0:
        return 0, np.zeros(0, int)
    e = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [0, 2]]])
    e = np.sort(e, axis=1)
    fid = np.tile(np.arange(n), 3)
    base = int(tri.max()) + 1
    key = e[:, 0].astype(np.int64) * base + e[:, 1]
    order = np.argsort(key, kind="stable")
    ks, fs = key[order], fid[order]
    same = ks[1:] == ks[:-1]               # facettes consecutives sur la
    rows, cols = fs[:-1][same], fs[1:][same]   # meme arete : un chainage suffit
    G = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    return connected_components(G, directed=False)


def union_find_components(tri):
    """Meme partition par une methode INDEPENDANTE (union-find pur Python,
    sans scipy) — sert de controle falsifiant (--check, --selftest)."""
    parent = list(range(len(tri)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    owner = {}
    for i, t in enumerate(tri):
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[0], t[2])):
            k = (min(a, b), max(a, b))
            j = owner.get(k)
            if j is None:
                owner[k] = i
            else:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[ri] = rj
    lab = np.array([find(i) for i in range(len(tri))])
    _, lab = np.unique(lab, return_inverse=True)
    return lab.ravel()


def same_partition(la, lb):
    """Deux etiquetages definissent-ils la MEME partition ?"""
    if len(la) != len(lb):
        return False
    pa = {}
    pb = {}
    for i, (x, y) in enumerate(zip(la, lb)):
        pa.setdefault(x, []).append(i)
        pb.setdefault(y, []).append(i)
    sa = set(tuple(v) for v in pa.values())
    sb = set(tuple(v) for v in pb.values())
    return sa == sb


def facet_geometry(P, cx, cy, zsurf):
    """P (n,3,3) sommets. Centroide, normale unitaire, aire, r du centroide,
    |n.e_r|, |n_z|, r min/max des sommets, profondeurs min/max (m)."""
    c = P.mean(axis=1)
    nv = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    area = 0.5 * np.linalg.norm(nv, axis=1)
    n = nv / np.maximum(2.0 * area, 1e-30)[:, None]
    dx, dy = c[:, 0] - cx, c[:, 1] - cy
    r = np.hypot(dx, dy)
    er = np.zeros_like(c)
    ok = r > 1e-12
    er[ok, 0] = dx[ok] / r[ok]
    er[ok, 1] = dy[ok] / r[ok]
    ner = np.abs((n * er).sum(axis=1))
    nz = np.abs(n[:, 2])
    rv = np.hypot(P[:, :, 0] - cx, P[:, :, 1] - cy)
    dv = zsurf - P[:, :, 2]
    return dict(c=c, n=n, area=area, r=r, ner=ner, nz=nz,
                rvmin=rv.min(axis=1), rvmax=rv.max(axis=1),
                dmin=dv.min(axis=1), dmax=dv.max(axis=1),
                depth=zsurf - c[:, 2], az=np.arctan2(dy, dx))


def orient_class(ner, nz, thr_r, thr_z):
    if nz >= thr_z:
        return "horizontale"
    if ner < thr_r:
        return "radiale"
    return "conique"


def describe(idx, g, thr_r, thr_z):
    """Statistiques d'un groupe de facettes (indices idx dans g)."""
    w = g["area"][idx]
    W = max(float(w.sum()), 1e-30)
    ner = float((g["ner"][idx] * w).sum() / W)
    nz = float((g["nz"][idx] * w).sum() / W)
    az = g["az"][idx]
    vx = float((np.cos(az) * w).sum() / W)
    vy = float((np.sin(az) * w).sum() / W)
    rmin = float(g["rvmin"][idx].min())
    rmax = float(g["rvmax"][idx].max())
    itip = int(idx[np.argmax(g["rvmax"][idx])])
    return dict(n=int(len(idx)), area=float(W), rmin=rmin, rmax=rmax,
                length=rmax - rmin, dmin=float(g["dmin"][idx].min()),
                dmax=float(g["dmax"][idx].max()), ner=ner, nz=nz,
                orient=orient_class(ner, nz, thr_r, thr_z),
                az=float(np.degrees(np.arctan2(vy, vx))), tip=itip)


# ---------------------------------------------------------------------------
# analyse d'une population de facettes rompues
# ---------------------------------------------------------------------------
def analyse(P, tri, cx=CX, cy=CY, zsurf=Z_SURF, rcore=6e-3, skin=1e-3,
            thr_r=0.35, thr_z=0.6, nsect=12):
    """P (n,3,3) sommets des facettes rompues (geometrie de reference),
    tri (n,3) sommets unifies. Retourne un dictionnaire de resultats."""
    n = len(tri)
    g = facet_geometry(P, cx, cy, zsurf)
    ncomp, lab = facet_adjacency_components(tri)
    assert len(lab) == n
    # tri des composantes par taille decroissante -> id 0 = la plus grande
    sizes = np.bincount(lab, minlength=ncomp) if ncomp else np.zeros(0, int)
    order = np.argsort(-sizes, kind="stable")
    rank = np.empty(ncomp, int)
    rank[order] = np.arange(ncomp)
    lab = rank[lab] if n else lab
    comps = []
    for k in range(ncomp):
        idx = np.where(lab == k)[0]
        d = describe(idx, g, thr_r, thr_z)
        d["id"] = k
        d["central"] = d["rmin"] < rcore
        comps.append(d)
    central = [c for c in comps if c["central"]]
    noyau = central[0]["id"] if central else None   # la plus grande (tri)
    for c in comps:
        if c["id"] == noyau:
            c["classe"] = "noyau"
        elif c["central"]:
            c["classe"] = "centrale"
        else:
            c["classe"] = "peripherique"
    # bras du noyau : facettes du noyau a r(centroide) > rcore
    arms = []
    arm_lab = -np.ones(n, int)
    if noyau is not None:
        out = np.where((lab == noyau) & (g["r"] > rcore))[0]
        if len(out):
            na, la = facet_adjacency_components(tri[out])
            sz = np.bincount(la, minlength=na)
            for k in np.argsort(-sz, kind="stable"):
                idx = out[la == k]
                d = describe(idx, g, thr_r, thr_z)
                d["id"] = len(arms)
                arms.append(d)
                arm_lab[idx] = d["id"]
    yang = 0.0
    yang_arm = None
    for a in arms:
        if a["orient"] == "radiale" and a["rmax"] > yang:
            yang, yang_arm = a["rmax"], a["id"]
    any_arm = max((a["rmax"] for a in arms), default=0.0)
    # cratere : facettes de surface du noyau
    crater = dict(rmax=0.0, rmean_sect=0.0, rmin_sect=0.0, cover=0.0, n=0)
    if noyau is not None:
        surf = np.where((lab == noyau) & (g["depth"] < skin))[0]
        crater["n"] = int(len(surf))
        if len(surf):
            crater["rmax"] = float(g["rvmax"][surf].max())
            edges = np.linspace(-np.pi, np.pi, nsect + 1)
            az = g["az"][surf]
            per = []
            for s in range(nsect):
                m = (az >= edges[s]) & (az < edges[s + 1])
                if m.any():
                    per.append(float(g["rvmax"][surf][m].max()))
            crater["cover"] = len(per) / float(nsect)
            crater["rmean_sect"] = float(np.mean(per)) if per else 0.0
            crater["rmin_sect"] = float(np.min(per)) if per else 0.0
    return dict(n=n, ncomp=ncomp, lab=lab, comps=comps, noyau=noyau,
                arms=arms, arm_lab=arm_lab, yang=yang, yang_arm=yang_arm,
                any_arm=any_arm, crater=crater, g=g,
                rmax_all=float(g["rvmax"].max()) if n else 0.0,
                rmax_centroid=float(g["r"].max()) if n else 0.0,
                depth_all=float(g["dmax"].max()) if n else 0.0)


# ---------------------------------------------------------------------------
# lecture d'un run
# ---------------------------------------------------------------------------
def load_run(run, frame, current=False):
    ks = frames_of(run)
    k = ks[-1] if frame < 0 else frame
    tk = frame_times(run).get(k, float("nan"))
    pts, con, f = read_vtu(joints_frame(run, k))
    p0, con0, _ = read_vtu(joints_frame(run, ks[0]))
    if not np.array_equal(con, con0):
        sys.exit("connectivite differente entre la trame 0 et la trame %d" % k)
    sel = broken_mask(f)
    inv = unify_vertices(p0)
    tri = inv[con[sel]]
    P = (pts if current else p0)[con[sel]]
    Pcur = pts[con[sel]]
    return dict(k=k, t=tk, P=P, Pcur=Pcur, tri=tri, nsel=int(sel.sum()),
                mode=(f["breakMode"][sel] if "breakMode" in f
                      else np.zeros(int(sel.sum()))))


# ---------------------------------------------------------------------------
# sorties
# ---------------------------------------------------------------------------
def fmt_row(c, classe):
    return ("%4d  %-12s %6d  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %5.2f  %5.2f  "
            "%-11s %7.1f" % (c["id"], classe, c["n"], c["rmin"] * 1e3,
                             c["rmax"] * 1e3, c["length"] * 1e3,
                             c["dmin"] * 1e3, c["dmax"] * 1e3, c["ner"],
                             c["nz"], c["orient"], c["az"]))


HEAD = ("  id  classe        nfac    rmin    rmax    long   zmin    zmax   "
        "|n.er|  |nz|  orient      azimut")


def report(res, label, rcore, skin, thr_r, thr_z, nmax=25):
    print("== %s" % label)
    print("   facettes rompues (tBreak >= 0) : %d ; composantes connexes "
          "(arete partagee) : %d" % (res["n"], res["ncomp"]))
    ncl = {}
    for c in res["comps"]:
        ncl[c["classe"]] = ncl.get(c["classe"], 0) + 1
    print("   noyau : %s ; centrales secondaires (r < %.0f mm) : %d ; "
          "peripheriques : %d"
          % ("id %d, %d facettes" % (res["noyau"], res["comps"][res["noyau"]]["n"])
             if res["noyau"] is not None else "AUCUN",
             rcore * 1e3, ncl.get("centrale", 0), ncl.get("peripherique", 0)))
    print("   r max de TOUTES les facettes rompues : sommets %.2f mm, "
          "centroides %.2f mm (ancienne metrique) ; profondeur %.2f mm"
          % (res["rmax_all"] * 1e3, res["rmax_centroid"] * 1e3,
             res["depth_all"] * 1e3))
    print("   bras du noyau (facettes a r > %.0f mm) : %d ; plus long bras "
          "toutes orientations : pointe a %.2f mm"
          % (rcore * 1e3, len(res["arms"]), res["any_arm"] * 1e3))
    print("   LONGUEUR DE FISSURE RADIALE (Yang) = %.2f mm%s   [Yang 9 m/s ~10]"
          % (res["yang"] * 1e3,
             " (bras %d)" % res["yang_arm"] if res["yang_arm"] is not None
             else " (aucun bras radial connecte au noyau)"))
    cr = res["crater"]
    print("   RAYON DE CRATERE (surface du noyau, < %.0f mm sous la surface, "
          "%d facettes) : r max %.2f mm ; par 12 secteurs : moyenne %.2f, "
          "min %.2f, couverture %.0f %%   [Yang 9 m/s ~7]"
          % (skin * 1e3, cr["n"], cr["rmax"] * 1e3, cr["rmean_sect"] * 1e3,
             cr["rmin_sect"] * 1e3, 100 * cr["cover"]))
    print("   -- bras du noyau (r en mm sur les sommets, z = profondeur) --")
    print(HEAD)
    for a in res["arms"][:nmax]:
        print(fmt_row(a, "bras"))
    print("   -- composantes (les %d plus grandes sur %d) --"
          % (min(nmax, res["ncomp"]), res["ncomp"]))
    print(HEAD)
    for c in res["comps"][:nmax]:
        print(fmt_row(c, c["classe"]))
    print("   (radiale : |n.e_r| < %.2f et |n_z| < %.2f ; horizontale : "
          "|n_z| >= %.2f ; sinon conique)" % (thr_r, thr_z, thr_z))


def write_csv(res, stem, label, rcore, skin):
    with io.open(stem + "_components.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["kind", "id", "classe", "nfac", "area_mm2", "rmin_mm",
                    "rmax_mm", "length_mm", "zmin_mm", "zmax_mm", "n_er",
                    "n_z", "orient", "azimut_deg"])
        for kind, lst in (("component", res["comps"]), ("arm", res["arms"])):
            for c in lst:
                w.writerow([kind, c["id"], c.get("classe", "bras"), c["n"],
                            "%.4f" % (c["area"] * 1e6), "%.4f" % (c["rmin"] * 1e3),
                            "%.4f" % (c["rmax"] * 1e3), "%.4f" % (c["length"] * 1e3),
                            "%.4f" % (c["dmin"] * 1e3), "%.4f" % (c["dmax"] * 1e3),
                            "%.4f" % c["ner"], "%.4f" % c["nz"], c["orient"],
                            "%.1f" % c["az"]])
    cr = res["crater"]
    rows = [("label", label, "-"), ("nBroken", res["n"], "-"),
            ("nComponents", res["ncomp"], "-"),
            ("noyauId", res["noyau"] if res["noyau"] is not None else -1, "-"),
            ("noyauFacets", res["comps"][res["noyau"]]["n"]
             if res["noyau"] is not None else 0, "-"),
            ("nArms", len(res["arms"]), "-"),
            ("yangRadialLength", res["yang"], "m"),
            ("longestArmAny", res["any_arm"], "m"),
            ("craterRmax", cr["rmax"], "m"),
            ("craterRmeanSect", cr["rmean_sect"], "m"),
            ("craterRminSect", cr["rmin_sect"], "m"),
            ("craterCover", cr["cover"], "-"),
            ("rmaxAllVertices", res["rmax_all"], "m"),
            ("rmaxAllCentroids", res["rmax_centroid"], "m"),
            ("depthAll", res["depth_all"], "m"),
            ("rcore", rcore, "m"), ("skin", skin, "m")]
    with io.open(stem + "_summary.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "value", "unit"])
        for r in rows:
            w.writerow(r)


def figure(res, P, stem, title, lim, depth, sec, cx=CX, cy=CY, zsurf=Z_SURF):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.unicode_minus": False,
    })
    lab, alab, g = res["lab"], res["arm_lab"], res["g"]
    noyau = res["noyau"]
    cm10 = plt.get_cmap("tab10")
    cm20 = plt.get_cmap("tab20")
    R_INSERT, R_CRAT_Y, R_RAD_Y = 8.51, 7.0, 10.0

    def col_of(i):
        if noyau is not None and lab[i] == noyau:
            if alab[i] >= 0:
                return cm10(alab[i] % 10)
            return (0.45, 0.45, 0.45, 1.0)
        return cm20(lab[i] % 20)

    cols = np.array([col_of(i) for i in range(res["n"])]) if res["n"] else \
        np.zeros((0, 4))
    fig, (A, B) = plt.subplots(1, 2, figsize=(15.0, 7.2),
                               gridspec_kw=dict(width_ratios=[1.0, 1.15]))
    # (a) vue de dessus
    if res["n"]:
        po = (P[:, :, :2] - np.array([cx, cy])) * 1e3
        A.add_collection(PolyCollection(po, facecolors=cols, edgecolors=cols,
                                        linewidths=0.35, alpha=0.75))
    th = np.linspace(0, 2 * np.pi, 240)
    for r, c, ls, lb in ((R_INSERT, "#333", "--", "insert R 8,5"),
                         (R_CRAT_Y, "#1f4e79", ":", "Yang cratere 7"),
                         (R_RAD_Y, "#1f4e79", "-.", "Yang radiales 10")):
        A.plot(r * np.cos(th), r * np.sin(th), ls, color=c, lw=0.9, label=lb)
    cr = res["crater"]["rmax"] * 1e3
    crm = res["crater"]["rmean_sect"] * 1e3
    if cr > 0:
        A.plot(cr * np.cos(th), cr * np.sin(th), "-", color="#b22222", lw=1.3,
               label="cratere r max %.1f" % cr)
        A.plot(crm * np.cos(th), crm * np.sin(th), ":", color="#b22222", lw=1.3,
               label="cratere moyen 12 secteurs %.1f" % crm)
    azstar = None
    if res["yang_arm"] is not None:
        a = res["arms"][res["yang_arm"]]
        tp = P[a["tip"]]
        j = np.argmax(np.hypot(tp[:, 0] - cx, tp[:, 1] - cy))
        tx, ty = (tp[j, 0] - cx) * 1e3, (tp[j, 1] - cy) * 1e3
        azstar = np.arctan2(ty, tx)
        A.plot([0, tx], [0, ty], "-", color="k", lw=1.0)
        A.plot([tx], [ty], "*", color="k", ms=11,
               label="pointe radiale %.1f mm" % (res["yang"] * 1e3))
    elif res["n"]:
        i = int(np.argmax(g["rvmax"]))
        azstar = float(g["az"][i])
    for c in res["comps"][:12]:
        if c["classe"] == "noyau" or c["n"] < 3:
            continue
        i = np.where(lab == c["id"])[0]
        xc = (g["c"][i, 0].mean() - cx) * 1e3
        yc = (g["c"][i, 1].mean() - cy) * 1e3
        A.annotate("%d" % c["id"], (xc, yc), fontsize=7, color="k")
    A.set_xlim(-lim, lim); A.set_ylim(-lim, lim); A.set_aspect("equal")
    A.set_xlabel("x [mm]"); A.set_ylabel("y [mm]")
    A.set_title("(a)  vue de dessus — noyau gris, bras du noyau (tab10), "
                "composantes separees (tab20)", loc="left", fontsize=10)
    A.legend(frameon=False, fontsize=8, loc="upper right")
    # (b) coupe verticale dans le plan de l'axe et de la pointe
    if azstar is None:
        azstar = 0.0
    ca, sa = np.cos(azstar), np.sin(azstar)
    if res["n"]:
        dx = P[:, :, 0] - cx
        dy = P[:, :, 1] - cy
        s = dx * ca + dy * sa
        d = -dx * sa + dy * ca
        keep = np.abs(d.mean(axis=1)) < sec * 1e-3
        if keep.any():
            po = np.stack([s[keep] * 1e3, (P[keep, :, 2] - zsurf) * 1e3], axis=2)
            B.add_collection(PolyCollection(po, facecolors=cols[keep],
                                            edgecolors=cols[keep],
                                            linewidths=0.35, alpha=0.8))
    B.axhline(0, color="#333", lw=0.9)
    for r, c, ls in ((R_INSERT, "#333", "--"), (R_CRAT_Y, "#1f4e79", ":"),
                     (R_RAD_Y, "#1f4e79", "-.")):
        B.axvline(-r, color=c, ls=ls, lw=0.8); B.axvline(r, color=c, ls=ls, lw=0.8)
    if res["yang"] > 0:
        B.axvline(res["yang"] * 1e3, color="k", lw=1.0)
    B.set_xlim(-lim, lim); B.set_ylim(-depth, 5); B.set_aspect("equal")
    B.set_xlabel("s le long de l'azimut %.0f$^\\circ$ [mm]" % np.degrees(azstar))
    B.set_ylabel("z sous la surface [mm]")
    B.set_title("(b)  coupe verticale par l'axe et la pointe radiale, "
                "$|d| <$ %g mm" % sec, loc="left", fontsize=10)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    for ext in ("pdf", "png"):
        fig.savefig(stem + "." + ext, dpi=170)
    plt.close(fig)
    print("   ecrit : %s.pdf / .png" % stem)


# ---------------------------------------------------------------------------
# mini-test FALSIFIANT (aucun VTU : population synthetique de facettes)
# ---------------------------------------------------------------------------
def _synthetic(with_link=True):
    """Sommets (mm) et triangles ; retourne (P (n,3,3) en m, tri, attendu).
    Noyau = eventail de 4 facettes horizontales a z = -0,5 mm (r <= 2 mm)
    + deux facettes verticales de raccord (r de 0 a 6,5 mm) ; bras radial =
    bande verticale dans le plan y = 0 de x = 6,5 a 16,5 mm (z de -0,5 a
    -6 mm), 10 facettes, entierement au-dela de rcore = 6 mm ;
    peripheriques : bande radiale dans le plan x = 0, y de 14 a 18 mm (4
    facettes), une facette horizontale a r ~ 9 mm, une facette
    circonferentielle a r = 8 mm (normale e_r).
    Sans le raccord (with_link = False) la bande est une composante
    peripherique : la longueur de Yang DOIT tomber a 0."""
    V = {}

    def vid(p):
        k = tuple(np.round(p, 9))
        if k not in V:
            V[k] = len(V)
        return V[k]

    T = []

    def tri(a, b, c):
        T.append((vid(a), vid(b), vid(c)))

    # eventail (surface, z = -0,5 mm) : centre (0,0), rayon 2
    o = (0.0, 0.0, -0.5)
    ring = [(2.0, 0.0, -0.5), (0.0, 2.0, -0.5), (-2.0, 0.0, -0.5), (0.0, -2.0, -0.5)]
    for i in range(4):
        tri(o, ring[i], ring[(i + 1) % 4])
    # raccord vertical (deux facettes) : L1 partage l'arete (o, (2,0,-0.5))
    # avec l'eventail, L2 partage ((2,0,-0.5),(6.5,0,-6)) avec L1 et
    # ((6.5,0,-0.5),(6.5,0,-6)) avec la premiere facette de la bande
    if with_link:
        tri(o, (2.0, 0.0, -0.5), (6.5, 0.0, -6.0))
        tri((2.0, 0.0, -0.5), (6.5, 0.0, -0.5), (6.5, 0.0, -6.0))
    # bande radiale y = 0, x = 6,5..16,5, z = -0.5 (haut) a -6 (bas)
    xs = [6.5, 8.5, 10.5, 12.5, 14.5, 16.5]
    for x0, x1 in zip(xs[:-1], xs[1:]):
        tri((x0, 0.0, -0.5), (x1, 0.0, -0.5), (x0, 0.0, -6.0))
        tri((x1, 0.0, -0.5), (x1, 0.0, -6.0), (x0, 0.0, -6.0))
    # peripherique radiale, plan x = 0, y = 14..18, z = -1 a -4
    ys = [14.0, 16.0, 18.0]
    for y0, y1 in zip(ys[:-1], ys[1:]):
        tri((0.0, y0, -1.0), (0.0, y1, -1.0), (0.0, y0, -4.0))
        tri((0.0, y1, -1.0), (0.0, y1, -4.0), (0.0, y0, -4.0))
    # horizontale isolee a r ~ 9 mm
    tri((8.0, 3.0, -2.0), (10.0, 3.0, -2.0), (9.0, 5.0, -2.0))
    # circonferentielle a r = 8 mm dans la direction -x (normale e_x = e_r)
    tri((-8.0, -1.0, -1.0), (-8.0, 1.0, -1.0), (-8.0, 0.0, -3.0))
    pts = np.array(list(V.keys())) * 1e-3
    tri_a = np.array(T, int)
    P = pts[tri_a]
    exp = dict(ncomp=4 if with_link else 5,
               noyau_n=16 if with_link else 4,
               yang=16.5e-3 if with_link else 0.0,
               crater=2e-3, n=len(T))
    return P, tri_a, exp


def selftest():
    ok = True

    def check(name, cond, detail=""):
        nonlocal ok
        print("   [%s] %s %s" % ("OK " if cond else "ECHEC", name, detail))
        ok = ok and bool(cond)

    print("== selftest crack_paths (population synthetique, axe (0,0), "
          "surface z = 0)")
    for with_link, tag in ((True, "A : bras RELIE au noyau"),
                           (False, "B : raccord retire (variante qui DOIT "
                                   "changer le verdict)")):
        P, tri, exp = _synthetic(with_link)
        res = analyse(P, tri, cx=0.0, cy=0.0, zsurf=0.0, rcore=6e-3,
                      skin=1e-3)
        print("   -- %s : %d facettes, %d composantes, noyau %s, Yang %.6f mm, "
              "cratere %.6f mm" % (tag, res["n"], res["ncomp"],
                                   res["noyau"], res["yang"] * 1e3,
                                   res["crater"]["rmax"] * 1e3))
        check("nombre de composantes = %d" % exp["ncomp"],
              res["ncomp"] == exp["ncomp"], "(mesure %d)" % res["ncomp"])
        check("facettes du noyau = %d" % exp["noyau_n"],
              res["comps"][res["noyau"]]["n"] == exp["noyau_n"],
              "(mesure %d)" % res["comps"][res["noyau"]]["n"])
        check("longueur radiale Yang = %.3f mm a 1e-9" % (exp["yang"] * 1e3),
              abs(res["yang"] - exp["yang"]) < 1e-9,
              "(mesure %.9f mm)" % (res["yang"] * 1e3))
        check("rayon de cratere = %.3f mm a 1e-9" % (exp["crater"] * 1e3),
              abs(res["crater"]["rmax"] - exp["crater"]) < 1e-9,
              "(mesure %.9f mm)" % (res["crater"]["rmax"] * 1e3))
        check("partition = union-find independant",
              same_partition(res["lab"], union_find_components(tri)))
        check("chaque facette dans exactement une composante",
              int(np.bincount(res["lab"]).sum()) == res["n"])
        orient = {c["orient"] for c in res["comps"] if c["classe"] != "noyau"}
        check("orientations des peripheriques = {radiale, horizontale, conique}",
              orient == {"radiale", "horizontale", "conique"},
              "(mesure %s)" % sorted(orient))
        if with_link:
            check("un seul bras, radial, 10 facettes, pointe a 16,5 mm, "
                  "longueur 10 mm",
                  len(res["arms"]) == 1 and res["arms"][0]["orient"] == "radiale"
                  and res["arms"][0]["n"] == 10
                  and abs(res["arms"][0]["rmax"] - 16.5e-3) < 1e-9
                  and abs(res["arms"][0]["length"] - 10e-3) < 1e-9)
            check("r max des centroides (ancienne metrique) > Yang ? NON : "
                  "ici la peripherique a 18 mm domine",
                  res["rmax_centroid"] > res["yang"],
                  "(%.3f > %.3f mm : l'ancienne metrique prend une fissure "
                  "NON connectee)" % (res["rmax_centroid"] * 1e3,
                                     res["yang"] * 1e3))
        else:
            nper = sum(1 for c in res["comps"]
                       if c["classe"] == "peripherique" and c["orient"] == "radiale")
            check("aucun bras (le raccord manque) et DEUX composantes "
                  "peripheriques radiales (la bande + celle a 14-18 mm)",
                  len(res["arms"]) == 0 and nper == 2,
                  "(bras %d, peripheriques radiales %d)" % (len(res["arms"]), nper))
    print("== selftest : %s" % ("TOUT OK" if ok else "ECHEC"))
    return 0 if ok else 1


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="fissures connectees et cratere")
    ap.add_argument("run", nargs="?", default=None)
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--rcore", type=float, default=6.0, help="rayon du noyau, mm")
    ap.add_argument("--skin", type=float, default=1.0,
                    help="peau de surface pour le cratere, mm")
    ap.add_argument("--radial", type=float, default=0.35, help="seuil |n.e_r|")
    ap.add_argument("--nz", type=float, default=0.6, help="seuil |n_z|")
    ap.add_argument("--sec", type=float, default=1.5,
                    help="demi-epaisseur de la coupe (b), mm")
    ap.add_argument("--current", action="store_true",
                    help="geometrie courante au lieu de la reference")
    ap.add_argument("--check", action="store_true",
                    help="controle : partition = union-find independant")
    ap.add_argument("--stem", default=None)
    ap.add_argument("--title", default=None)
    ap.add_argument("--lim", type=float, default=0.0,
                    help="demi-etendue des deux panneaux, mm ; 0 = ajustee "
                         "aux donnees (au moins 12 mm : les cercles de Yang "
                         "a 7 et 10 mm restent visibles)")
    ap.add_argument("--depth", type=float, default=0.0,
                    help="profondeur du panneau (b), mm ; 0 = ajustee")
    ap.add_argument("--nofig", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    if a.run is None:
        ap.error("run manquant")
    d = load_run(a.run, a.frame, a.current)
    res = analyse(d["P"], d["tri"], rcore=a.rcore * 1e-3, skin=a.skin * 1e-3,
                  thr_r=a.radial, thr_z=a.nz)
    label = "%s  trame %d  t = %.1f us  (geometrie %s)" % (
        a.run, d["k"], d["t"] * 1e6, "courante" if a.current else "de reference")
    report(res, label, a.rcore * 1e-3, a.skin * 1e-3, a.radial, a.nz)
    if not a.current and d["nsel"]:
        rc = np.hypot(d["Pcur"][:, :, 0] - CX, d["Pcur"][:, :, 1] - CY)
        cc = d["Pcur"].mean(axis=1)
        rcc = np.hypot(cc[:, 0] - CX, cc[:, 1] - CY)
        print("   (geometrie COURANTE : r max des sommets %.2f mm, des "
              "centroides %.2f mm = l'ancienne metrique imp_lib.metrics)"
              % (rc.max() * 1e3, rcc.max() * 1e3))
    if a.check:
        same = same_partition(res["lab"], union_find_components(d["tri"]))
        print("   [%s] partition csgraph = union-find independant"
              % ("OK " if same else "ECHEC"))
        if not same:
            sys.exit(2)
    stem = a.stem or os.path.join("results", "fig",
                                  "crack_paths_" + os.path.basename(
                                      os.path.normpath(a.run)))
    os.makedirs(os.path.dirname(stem) or ".", exist_ok=True)
    write_csv(res, stem, label, a.rcore * 1e-3, a.skin * 1e-3)
    print("   ecrit : %s_components.csv / _summary.csv" % stem)
    if not a.nofig:
        # cadrage AJUSTE par defaut : 15 % de marge autour du plus grand rayon
        # mesure (facettes, cratere, cercles de Yang a 10 mm), arrondi a 2 mm.
        lim = a.lim
        if lim <= 0.0:
            r = max(res["rmax_all"] * 1e3, res["crater"]["rmax"] * 1e3, 10.0)
            lim = 2.0 * np.ceil(1.15 * r / 2.0)
        depth = a.depth
        if depth <= 0.0:
            depth = 2.0 * np.ceil(1.15 * res["depth_all"] * 1e3 / 2.0 + 0.5)
        # le r max du cratere est porte par une facette isolee sur un maillage
        # grossier : la moyenne par secteur figure a cote, jamais seule.
        title = (a.title or os.path.basename(os.path.normpath(a.run))) + \
            " — fissures connectees a $t$ = %.0f $\\mu$s : %d facettes, %d " \
            "composantes, Yang %.1f mm, cratere %.1f mm (r max) / %.1f mm " \
            "(moyenne 12 secteurs)" % (
                d["t"] * 1e6, res["n"], res["ncomp"], res["yang"] * 1e3,
                res["crater"]["rmax"] * 1e3, res["crater"]["rmean_sect"] * 1e3)
        figure(res, d["P"], stem, title, lim, depth, a.sec)


if __name__ == "__main__":
    main()
