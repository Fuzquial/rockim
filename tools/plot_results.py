#!/usr/bin/env python3
"""Quick-look plots from a rockim output directory (no ParaView needed).

usage: python3 tools/plot_results.py <out_dir> [--title "..."]

Produces, depending on what it finds in <out_dir>:
  <out_dir>/plot_field.png    damage map (FEM) or particles + broken bonds (DEM)
  <out_dir>/plot_history.png  tool force-time history
"""
import argparse
import csv
import os
import sys

import glob
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection   # chemin VTU (plot_fdem_vtu)


def vtu_array(text, name, dtype=float):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name,
                  text, re.S)
    return np.fromstring(m.group(1), sep=" ", dtype=dtype)


def vtu_points(text):
    m = re.search(r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)[:, :2]


def read_csv(path):
    with open(path) as f:
        r = csv.reader(f)
        head = next(r)
        rows = [[float(x) for x in row] for row in r if row]
    data = {h: np.array([row[i] for row in rows]) for i, h in enumerate(head)}
    return data


def plot_fem(out, ax):
    # Render the true triangulation of the last frame (the mesh is crossed
    # CST triangles; scattering element centroids, as an earlier version did,
    # produces a misleading dotted/"hexagonal" texture).
    vtus = sorted(glob.glob(os.path.join(out, "fem_[0-9]*.vtu")))
    if vtus:
        txt = open(vtus[-1]).read()
        pts = vtu_points(txt)
        conn = vtu_array(txt, "connectivity", int).reshape(-1, 3)
        dmg = vtu_array(txt, "damage")
        pc = ax.tripcolor(pts[:, 0], pts[:, 1], conn, facecolors=dmg,
                          cmap="inferno", vmin=0, vmax=1)
        plt.colorbar(pc, ax=ax, label="damage D", shrink=0.8)
        ax.set_title("FEM: damage field, deformed mesh (holes = eroded)")
        return
    d = read_csv(os.path.join(out, "fem_final_elements.csv"))
    live = d["eroded"] < 0.5
    sc = ax.scatter(d["cx"][live], d["cy"][live], c=d["damage"][live],
                    s=4, cmap="inferno", vmin=0, vmax=1, marker="s")
    ax.scatter(d["cx"][~live], d["cy"][~live], color="white", s=4, marker="s")
    plt.colorbar(sc, ax=ax, label="damage D", shrink=0.8)
    ax.set_title("FEM: damage field (white = eroded)")


EXP_JSON = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "..", "CONTINUUM", "calib_bohus_triax", "exp_qc",
    "experimental_data_red_bohus_clean.json")


def plot_exp(ax, sigma3, path=None):
    """Superpose les essais experimentaux Red Bohus a un confinement donne.

    Source : experimental_data_red_bohus_clean.json (v4, 2026-08-31), le jeu
    NETTOYE — tete q < 0 purgee, de-spike medfilt7, queue tronquee au max
    d'eps_axial (decharge machine supprimee), adoucissement CONSERVE, UCS
    ecartes. Trace chaque essai en gris fin plus leur moyenne interpolee sur
    une grille commune, pour que la dispersion se voie : a 20 MPa les trois
    essais tiennent dans +- 2,8 MPa sur le pic, soit 0,7 %.
    Unites du fichier : deviateur en MPa, deformations en microstrain.
    """
    import json
    p = path or EXP_JSON
    if not os.path.exists(p):
        print("  (donnees experimentales introuvables : %s)" % p)
        return 0
    d = json.load(open(p, encoding="utf-8"))
    runs = [v for v in d["triaxial"].values() if v["sigma3_MPa"] == sigma3]
    if not runs:
        print("  (aucun essai a sigma3 = %g MPa)" % sigma3)
        return 0
    for i, v in enumerate(runs):
        e = np.array(v["eps_axial_microstrain"]) * 1e-4      # -> %
        q = np.array(v["deviator_stress_MPa"])
        ax.plot(e, q, color="0.55", lw=0.9, zorder=1,
                label=("essais exp. (%d)" % len(runs)) if i == 0 else None)
    emax = min(np.array(v["eps_axial_microstrain"]).max() for v in runs) * 1e-4
    gr = np.linspace(0, emax, 400)
    M = np.array([np.interp(gr, np.array(v["eps_axial_microstrain"]) * 1e-4,
                            np.array(v["deviator_stress_MPa"])) for v in runs])
    ax.plot(gr, M.mean(0), color="k", lw=2.2, zorder=2,
            label="moyenne exp. ($\\sigma_3$ = %g MPa)" % sigma3)
    return len(runs)


def _delay_of(run_dir):
    """pullDelay du deck effectif : l'instant de fin de consolidation.

    Lu dans config_effective.cfg, que le solveur ecrit dans CHAQUE dossier de
    run avec toutes les cles consommees (celles du deck et celles au defaut).
    Defaut 0 si absent : la courbe part alors de sigma brut, comme avant.
    """
    p = os.path.join(run_dir, "config_effective.cfg")
    if not os.path.exists(p):
        return 0.0
    for line in open(p, encoding="utf-8", errors="replace"):
        line = line.split("#")[0]
        if "pullDelay" in line and "=" in line:
            try:
                return float(line.split("=")[1].strip())
            except ValueError:
                pass
    return 0.0


def plot_sigeps(dirs, axs, deviator=True):
    """Courbes contrainte-deformation d'essais sur eprouvette (2026-09-11).

    Ajout opt-in : `loading = platens | grips` ecrit dans history.csv les
    colonnes propres a l'essai (sigma, sigmaPeak, epsGauge, epsSpec,
    confAchieved, nInserted, nBroken, nBrokTen, nBrokShear, peakLocked) et
    AUCUN outil du depot ne les tracait — plot_force_penetration.py lit la
    force d'OUTIL, calibrate_bohus.py trace des posteriors. Se lit en cours de
    run : c'est fait pour suivre des essais qui n'ont pas encore atteint leur
    pic.

    DEVIATEUR q = sigma_ax - sigma_3 par DEFAUT (2026-09-11). C'est la
    convention du depot, etablie par la revue adverse du 2026-09-02 (memoire
    « q = sigma - sigma_3, E physique ») : la contrainte AXIALE BRUTE porte la
    part elastique du confinement lateral (mesure du 11/09 sur ces decks :
    +11,6 MPa a 40 MPa de confinement des eps = 0,02 %, soit nu du feldspath,
    la phase majoritaire — un decalage CONSTANT tant qu'aucun joint ne s'insere)
    et elle melange donc un artefact de convention a l'effet recherche. Le
    deviateur l'annule : a chargement purement elastique les quatre courbes se
    superposent, et tout ecart restant est mecanique.
    `--axial` restitue la contrainte brute.

    Deux panneaux : (1) q(eps) avec le pic marque quand il est verrouille,
    (2) population de joints — inseres, endommageants, rompus — sur la meme
    abscisse, parce que c'est la que se lit le lien microstructure/resistance.
    """
    A, B = axs
    cmap = plt.get_cmap("viridis")
    for i, d in enumerate(dirs):
        lab = os.path.basename(d.rstrip("/\\"))
        h = os.path.join(d, "history.csv")
        if not os.path.exists(h):
            print("  (pas d'historique) " + d)
            continue
        r = list(csv.DictReader(open(h)))
        if not r:
            continue
        g = lambda n: np.array([float(x[n]) for x in r]) if n in r[0] else None
        eps = g("epsGauge")
        sig = g("sigma")
        if eps is None or sig is None:
            print("  (pas de colonnes d'eprouvette) " + d)
            continue
        col = cmap(i / max(len(dirs) - 1, 1) * 0.85)
        conf = g("confAchieved")
        s3 = np.abs(conf) if conf is not None else np.zeros_like(sig)
        # DEVIATEUR : sigma MOINS SA VALEUR A LA FIN DE LA CONSOLIDATION, et
        # non sigma - sigma_3. C'est la convention du depot, lisible telle
        # quelle dans rockim_f1/calib_triax3d/courbes_2d.py:61
        #     q = sig - np.interp(DELAY, t, sig)   # offset fin de consolidation
        # Pourquoi ce n'est PAS sigma - sigma_3 : en fdem 2D le confinement
        # n'agit que sur les faces laterales et les mors bloquent les bouts
        # (FdemSolver.cpp:4881, les faces haut/bas sont exclues en TENSION).
        # L'etat de fin de consolidation n'est donc pas hydrostatique mais de
        # DEFORMATION UNIAXIALE : sigma_ax = nu sigma_3, pas sigma_3. Retrancher
        # sigma_3 fait partir la courbe a -(1-nu) sigma_3 (mesure du 11/09 :
        # -31 MPa a 40 MPa de confinement) au lieu de zero.
        # La vraie consolidation isotrope existe, mais en fem3d seulement
        # (topPressure = confiningPressure + pullDelay : les noeuds du mors
        # sont LIBRES pendant le delai ; banc falsifiant DOC:2809).
        t = g("t")
        y = sig
        if deviator:
            y = sig - np.interp(_delay_of(d), t, sig)
        cl = "%s  (sig3 = %.0f MPa)" % (lab, s3[-1] / 1e6) if conf is not None else lab
        A.plot(eps * 100, y / 1e6, color=col, lw=1.9, label=cl)
        lock = g("peakLocked")
        pk = g("sigmaPeak")
        if lock is not None and pk is not None and lock[-1] > 0.5:
            k = int(np.argmax(y))
            A.plot(eps[k] * 100, y[k] / 1e6, "o", ms=7, color=col)
        else:                       # run NON termine : marquer la fin courante
            A.plot(eps[-1] * 100, y[-1] / 1e6, "s", ms=6, color=col,
                   markerfacecolor="none")
        ins, nb = g("nInserted"), g("nBroken")
        if ins is not None:
            B.plot(eps * 100, ins, color=col, lw=1.7)
        if nb is not None:
            B.plot(eps * 100, nb, color=col, lw=1.7, ls="--")
    A.set_xlabel("deformation axiale de jauge (%)")
    A.set_ylabel((r"deviateur q = $\sigma_{ax} - \sigma_3$ (MPa)"
                  if deviator else "contrainte axiale (MPa)"))
    A.set_title(("Deviateur-deformation" if deviator else
                 "Contrainte axiale BRUTE (porte la part elastique du "
                 "confinement)")
                + "  (o = pic verrouille, carre creux = run EN COURS)",
                fontsize=10)
    A.grid(alpha=.3)
    A.legend(fontsize=8)
    B.set_xlabel("deformation axiale de jauge (%)")
    B.set_ylabel("nombre de joints")
    B.set_title("Population de joints : trait plein = inseres, "
                "tirets = rompus", fontsize=10)
    B.grid(alpha=.3)


def plot_facies(out, ax, frame=-1, dth=0.05, lw=1.6):
    """Facies de rupture : les joints ROMPUS, sur la silhouette de l'eprouvette.

    Ajout 2026-09-11. Repond a « y a-t-il LOCALISATION ? » — la question que ni
    la courbe contrainte-deformation ni les compteurs ne tranchent : un meme
    nombre de joints rompus peut former une bande nette ou un nuage diffus.

    Rouge = rompu en CISAILLEMENT (breakMode = 2), bleu = en TRACTION
    (breakMode = 1) ; gris clair = endommage sans rompre (D > dth). Le maillage
    n'est pas trace : il masque le facies. Rend (n_rompus, n_endommages).
    """
    vtus = sorted(glob.glob(os.path.join(out, "fdem_joints_[0-9]*.vtu")))
    if not vtus:
        raise SystemExit("aucun fdem_joints_XXXX.vtu dans " + out)
    txt = open(vtus[frame]).read()
    P = vtu_points(txt)
    off = vtu_array(txt, "offsets", int)
    n = int(off[0]) if len(off) else 2
    C = vtu_array(txt, "connectivity", int).reshape(-1, n)[:, :2]
    D = vtu_array(txt, "damage")
    bm = vtu_array(txt, "breakMode")
    bo = vtu_array(txt, "bonded")
    segs = P[C]
    ins = (bo < 0.5) if bo is not None else np.ones(len(C), bool)
    brk = ins & (D >= 1.0)
    dmg = ins & (D > dth) & (D < 1.0)
    # silhouette : contour de l'eprouvette, lu des bornes des noeuds
    x0, x1 = P[:, 0].min(), P[:, 0].max()
    y0, y1 = P[:, 1].min(), P[:, 1].max()
    ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                               ec="0.6", lw=0.8, zorder=1))
    if dmg.any():
        ax.add_collection(LineCollection(segs[dmg], colors="0.78",
                                         linewidths=0.6, zorder=2))
    for code, col, lab in ((2.0, "#B3121B", "cisaillement"),
                           (1.0, "#1f5fa8", "traction")):
        s = brk & (bm == code) if bm is not None else np.zeros(len(C), bool)
        if s.any():
            ax.add_collection(LineCollection(segs[s], colors=col,
                                             linewidths=lw, zorder=4,
                                             label=lab))
    other = brk & ~np.isin(bm, [1.0, 2.0]) if bm is not None else brk
    if other.any():
        ax.add_collection(LineCollection(segs[other], colors="k",
                                         linewidths=lw, zorder=3))
    ax.set_xlim(x0 - 0.002, x1 + 0.002)
    ax.set_ylim(y0 - 0.002, y1 + 0.002)
    return int(brk.sum()), int(dmg.sum())


def plot_fdem_vtu(out, ax, field="phase", frame=-1, joints=True):
    """Maillage 2D fdem d'une TRAME, lu directement du VTU (2026-09-11).

    Ajout opt-in : le chemin historique de `plot_fdem` lit les CSV de FIN de
    run, donc rien n'est visible tant que le run tourne, et il ne connait ni
    `phase` ni `grain`. Ici on rend la vraie triangulation de la trame
    demandee, coloree par un champ d'ELEMENT quelconque du VTU, avec les
    facettes de joint en surimpression — c'est ce qu'il faut pour verifier une
    microstructure GBM AVANT de lire un resultat (regle du 2026-09-07 :
    montrer le maillage avant de lancer).

    field  : 'phase', 'grain', 'vonMises', 'fragment'... (tout champ d'element)
    frame  : index de trame ; -1 = la derniere ecrite
    joints : trace les facettes de joint (gris fin), et en rouge celles qui
             separent DEUX GRAINS differents (frontieres de grain)
    """
    vtus = sorted(glob.glob(os.path.join(out, "fdem_[0-9]*.vtu")))
    if not vtus:
        raise SystemExit("aucun fdem_XXXX.vtu dans " + out)
    path = vtus[frame]
    txt = open(path).read()
    pts = vtu_points(txt)
    conn = vtu_array(txt, "connectivity", int).reshape(-1, 3)
    val = vtu_array(txt, field)
    disc = field in ("phase", "grain", "fragment")
    if disc:
        # champ DISCRET : EXACTEMENT une couleur par valeur presente. Une
        # colormap continue echantillonnee sur [vmin, vmax] affichait des
        # bandes intermediaires qui ne correspondent a aucune phase.
        from matplotlib.colors import ListedColormap, BoundaryNorm
        uniq = np.unique(val)
        idx = np.searchsorted(uniq, val).astype(float)
        n = len(uniq)
        base = plt.get_cmap("tab20" if n > 8 else "Set2")
        cmap = ListedColormap([base(i % base.N) for i in range(n)])
        norm = BoundaryNorm(np.arange(-0.5, n), n)
        pc = ax.tripcolor(pts[:, 0], pts[:, 1], conn, facecolors=idx,
                          cmap=cmap, norm=norm)
        if n <= 12:
            cb = plt.colorbar(pc, ax=ax, shrink=0.55, ticks=range(n))
            cb.ax.set_yticklabels([("%g" % u) for u in uniq])
        else:
            cb = plt.colorbar(pc, ax=ax, shrink=0.55)
        cb.set_label("%s (%d valeurs)" % (field, n))
    else:
        pc = ax.tripcolor(pts[:, 0], pts[:, 1], conn, facecolors=val,
                          cmap="viridis")
        plt.colorbar(pc, ax=ax, label=field, shrink=0.8)
    ntri = len(conn)
    njt = 0
    if joints:
        jp = path.replace("fdem_", "fdem_joints_")
        if os.path.exists(jp):
            jt = open(jp).read()
            jpts = vtu_points(jt)
            off = vtu_array(jt, "offsets", int)
            jc = vtu_array(jt, "connectivity", int)
            n = int(off[0]) if len(off) else 2
            jc = jc.reshape(-1, n)
            njt = len(jc)
            segs = [jpts[c[:2]] for c in jc]
            ax.add_collection(LineCollection(segs, colors="0.35",
                                             linewidths=0.35, zorder=3))
            # frontieres de GRAIN : les deux elements adjacents n'ont pas le
            # meme `grain`. On les retrouve geometriquement par le centre du
            # segment, seule information commune aux deux VTU.
            g = vtu_array(txt, "grain")
            cen = pts[conn].mean(axis=1)
            from scipy.spatial import cKDTree           # deja utilise ailleurs
            tree = cKDTree(cen)
            mid = np.array([s.mean(axis=0) for s in segs])
            d, nb = tree.query(mid, k=2)
            gb = g[nb[:, 0]] != g[nb[:, 1]]
            if gb.any():
                ax.add_collection(LineCollection(
                    [segs[i] for i in np.where(gb)[0]], colors="#B3121B",
                    linewidths=0.9, zorder=4))
    ax.set_title("%s · %d triangles, %d joints%s"
                 % (field, ntri, njt,
                    ", rouge = frontiere de grain" if joints else ""),
                 fontsize=10)


def plot_dem(out, ax):
    p = read_csv(os.path.join(out, "dem_final_particles.csv"))
    b = read_csv(os.path.join(out, "dem_final_bonds.csv"))
    broken = b["broken"] > 0.5
    for x1, y1, x2, y2 in zip(b["x1"][broken], b["y1"][broken],
                              b["x2"][broken], b["y2"][broken]):
        ax.plot([x1, x2], [y1, y2], color="red", lw=0.5, alpha=0.6, zorder=3)
    frag = p["fragment"]
    main = frag == 0
    ax.scatter(p["x"][main], p["y"][main], s=3, color="0.6", zorder=2)
    sc = ax.scatter(p["x"][~main], p["y"][~main], c=frag[~main], s=3,
                    cmap="turbo", zorder=2)
    ax.set_title("DEM: fragments (grey = main body) + broken bonds (red)")


def plot_dem3d(out, ax):
    """Mid-depth slice (|y - yc| < one particle diameter) of the final state."""
    p = read_csv(os.path.join(out, "dem3d_final_particles.csv"))
    yc = 0.5 * (p["y"].min() + p["y"].max())
    sel = np.abs(p["y"] - yc) < 2.05 * p["r"]
    frag = p["fragment"][sel]
    main = frag == 0
    ax.scatter(p["x"][sel][main], p["z"][sel][main], s=6, color="0.6", zorder=2)
    if (~main).any():
        ax.scatter(p["x"][sel][~main], p["z"][sel][~main], c=frag[~main], s=6,
                   cmap="turbo", zorder=3)
    ax.set_title("DEM3D: mid-depth slice, fragments (grey = main body)")


def plot_fdem(out, ax):
    e = read_csv(os.path.join(out, "fdem_final_elements.csv"))
    j = read_csv(os.path.join(out, "fdem_final_joints.csv"))
    frag = e["fragment"]
    main_b = frag == 0
    ax.scatter(e["cx"][main_b], e["cy"][main_b], s=4, color="0.75", zorder=1)
    if (~main_b).any():
        ax.scatter(e["cx"][~main_b], e["cy"][~main_b], c=frag[~main_b], s=5,
                   cmap="turbo", zorder=3)
    br = j["damage"] >= 0.99
    if br.any():
        from matplotlib.collections import LineCollection
        segs = np.stack([np.stack([j["x1"][br], j["y1"][br]], 1),
                         np.stack([j["x2"][br], j["y2"][br]], 1)], 1)
        ax.add_collection(LineCollection(segs, colors="crimson", lw=0.5, zorder=2))
    ax.set_title("FDEM: fragments (grey = main body) + broken joints (red)")


def plot_fdem3d(out, ax):
    e = read_csv(os.path.join(out, "fdem3d_final_elements.csv"))
    yc = 0.5 * (e["cy"].min() + e["cy"].max())
    hm = np.median(np.diff(np.unique(np.round(e["cy"], 5)))) * 3
    sel = np.abs(e["cy"] - yc) < hm
    frag = e["fragment"][sel]
    mb = frag == 0
    ax.scatter(e["cx"][sel][mb], e["cz"][sel][mb], s=8, color="0.75", zorder=1)
    if (~mb).any():
        ax.scatter(e["cx"][sel][~mb], e["cz"][sel][~mb], c=frag[~mb], s=9,
                   cmap="turbo", zorder=3)
    ax.set_title("FDEM3D: mid-depth slice, fragments (grey = main body)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out", nargs="+",
                    help="dossier(s) de run ; plusieurs seulement avec --sigeps")
    ap.add_argument("--facies", action="store_true",
                    help="planche de FACIES : joints rompus (rouge = "
                         "cisaillement, bleu = traction) sur la silhouette ; "
                         "accepte plusieurs dossiers")
    ap.add_argument("--labels", default=None,
                    help="etiquettes des panneaux, separees par des virgules")
    ap.add_argument("--exp", type=float, default=None,
                    help="superpose les essais Red Bohus a ce confinement "
                         "(MPa) : 20, 50, 75 ou 100")
    ap.add_argument("--axial", action="store_true",
                    help="tracer la contrainte AXIALE BRUTE au lieu du "
                         "deviateur q = sigma_ax - sigma_3 (defaut)")
    ap.add_argument("--sigeps", action="store_true",
                    help="courbes contrainte-deformation d'eprouvette "
                         "(loading = platens | grips), lisibles en cours de run")
    ap.add_argument("--title", default="")
    # ---- chemin VTU opt-in (2026-09-11) : voir plot_fdem_vtu -------------
    ap.add_argument("--field", default=None,
                    help="champ d'element d'un VTU fdem 2D (phase, grain, "
                         "vonMises...) ; sans lui, comportement historique")
    ap.add_argument("--frame", type=int, default=-1,
                    help="index de trame VTU (-1 = la derniere ecrite)")
    ap.add_argument("--no-joints", action="store_true")
    ap.add_argument("--out-png", default=None)
    ap.add_argument("--zoom", nargs=4, type=float, default=None,
                    metavar=("X0", "X1", "Y0", "Y1"),
                    help="fenetre en metres : lire la triangulation INTRA-grain")
    ap.add_argument("--figsize", nargs=2, type=float, default=None)
    a = ap.parse_args()

    if a.facies:
        N = len(a.out)
        ncol = min(N, 4); nrow = (N + ncol - 1) // ncol
        fig, axs = plt.subplots(nrow, ncol,
                                figsize=tuple(a.figsize or [3.1 * ncol, 5.6 * nrow]))
        axs = np.atleast_1d(axs).ravel()
        labs = (a.labels.split(",") if a.labels
                else [os.path.basename(d.rstrip("/\\")) for d in a.out])
        for i, d in enumerate(a.out):
            nb, nd = plot_facies(d, axs[i], a.frame)
            axs[i].set_aspect("equal"); axs[i].set_xticks([]); axs[i].set_yticks([])
            axs[i].set_title("%s — %d rompus, %d endommages"
                             % (labs[i].strip(), nb, nd), fontsize=9)
        for j in range(N, len(axs)):
            axs[j].axis("off")
        h, l = axs[0].get_legend_handles_labels()
        for i in range(1, N):
            hh, ll = axs[i].get_legend_handles_labels()
            for x, y in zip(hh, ll):
                if y not in l: h.append(x); l.append(y)
        if h:
            fig.legend(h, l, loc="lower center", ncol=2, fontsize=9)
        if a.title:
            fig.suptitle(a.title)
        p = a.out_png or "facies.png"
        fig.tight_layout(rect=[0, 0.04, 1, 0.97])
        fig.savefig(p, dpi=150)
        print("wrote", p)
        return

    if a.sigeps:
        fig, axs = plt.subplots(1, 2, figsize=tuple(a.figsize or [13.0, 5.2]))
        if a.exp is not None:
            n = plot_exp(axs[0], a.exp)
            print("  %d essai(s) experimentaux a sigma3 = %g MPa" % (n, a.exp))
        plot_sigeps(a.out, axs, deviator=not a.axial)
        if a.title:
            fig.suptitle(a.title)
        p = a.out_png or "sigeps.png"
        fig.tight_layout()
        fig.savefig(p, dpi=150)
        print("wrote", p)
        return

    out0 = a.out[0]
    if a.field:
        fs = a.figsize or ([7.0, 7.0] if a.zoom else [6.5, 10.0])
        fig, ax = plt.subplots(figsize=tuple(fs))
        plot_fdem_vtu(out0, ax, a.field, a.frame, not a.no_joints)
        if a.zoom:
            ax.set_xlim(a.zoom[0], a.zoom[1])
            ax.set_ylim(a.zoom[2], a.zoom[3])
        ax.set_aspect("equal")
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        if a.title:
            fig.suptitle(a.title)
        p = a.out_png or os.path.join(out0, "plot_%s.png" % a.field)
        fig.tight_layout()
        fig.savefig(p, dpi=150)
        print("wrote", p)
        return

    is_fem = os.path.exists(os.path.join(out0, "fem_final_elements.csv"))
    is_dem = os.path.exists(os.path.join(out0, "dem_final_particles.csv"))
    is_d3 = os.path.exists(os.path.join(out0, "dem3d_final_particles.csv"))
    is_fd = os.path.exists(os.path.join(out0, "fdem_final_elements.csv"))
    is_f3 = os.path.exists(os.path.join(out0, "fdem3d_final_elements.csv"))
    if not (is_fem or is_dem or is_d3 or is_fd or is_f3):
        sys.exit("no rockim final CSVs found in " + out0)

    fig, ax = plt.subplots(figsize=(9, 5))
    (plot_fem if is_fem else plot_dem if is_dem else
     plot_dem3d if is_d3 else plot_fdem if is_fd else plot_fdem3d)(out0, ax)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    if a.title:
        fig.suptitle(a.title)
    fig.tight_layout()
    fig.savefig(os.path.join(out0, "plot_field.png"), dpi=150)

    hist = os.path.join(out0, "history.csv")
    if os.path.exists(hist):
        h = read_csv(hist)
        if "toolFz" in h:
            fig2, ax2 = plt.subplots(figsize=(9, 4))
            ax2.plot(h["t"] * 1e3, h["toolFz"] / 1e3, label="F$_z$ on tool")
            ax2.plot(h["t"] * 1e3, h["toolFx"] / 1e3, label="F$_x$ on tool", alpha=0.8)
            ax2.set_xlabel("t [ms]")
            ax2.set_ylabel("force [kN]")
            ax2.legend()
            ax2.grid(alpha=0.3)
            ax2.set_title("Tool force history" + (" — " + a.title if a.title else ""))
            fig2.tight_layout()
            fig2.savefig(os.path.join(out0, "plot_history.png"), dpi=150)
        elif "toolFy" in h:
            fig2, ax2 = plt.subplots(figsize=(9, 4))
            ax2.plot(h["t"] * 1e3, h["toolFy"] / 1e6, label="F$_y$ on tool")
            ax2.plot(h["t"] * 1e3, h["toolFx"] / 1e6, label="F$_x$ on tool", alpha=0.8)
            ax2.set_xlabel("t [ms]")
            ax2.set_ylabel("force [MN/m]")
            ax2.legend()
            ax2.grid(alpha=0.3)
            ax2.set_title("Tool force history" + (" — " + a.title if a.title else ""))
            fig2.tight_layout()
            fig2.savefig(os.path.join(out0, "plot_history.png"), dpi=150)


if __name__ == "__main__":
    main()
