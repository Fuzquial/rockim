#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# crater_metrics.py — metriques de CRATERE d'un impact fdem3d (V2/B3).
#
#   python3 tools/crater_metrics.py <dossier_run> [options]
#
#   --frame N          frame a analyser (defaut : la derniere)
#   --center X Y       axe d'impact (defaut : centroide des joints casses,
#                      pondere par leur aire)
#   --dth D            seuil d'endommagement pour la zone endommagee (0.05)
#   --sectors N        secteurs angulaires des fissures radiales (12)
#   --skin H           epaisseur de la peau de surface [m] (defaut : 2 fois
#                      la taille mediane des joints)
#   --plot out.png     vue de dessus (joints endommages/casses, rayon,
#                      portees radiales par secteur)
#   --csv out.csv      exporte les metriques en CSV
#
# Metriques (definitions du banc Yang IJRMMS 2025, adaptees au maillage) :
#   R_crater  : rayon englobant (p95 radial) des joints CASSES dans la PEAU
#               de surface — le bol visible du cratere ;
#   R_max     : portee radiale max des joints casses, toutes profondeurs ;
#   depth     : profondeur max des joints casses sous la surface libre ;
#   A_broken  : aire cassee cumulee [m2] (les levres de fissure / 2) ;
#   V_damaged : volume des elements adjacents a un joint D > dth [m3] ;
#   V_detached: volume des fragments detaches (fragment != 0, VTU elements) ;
#   fissures radiales : par secteur angulaire, portee radiale max des joints
#               casses AU-DELA de R_crater (la longueur de fissure radiale
#               du banc), + moyenne et max.
#
# Post-traitement PUR (aucune influence sur le solveur). Lit les VTU ASCII
# de rockim (fdem3d_XXXX.vtu + fdem3d_joints_XXXX.vtu).
# ---------------------------------------------------------------------------
import argparse
import glob
import os
import re
import sys

import numpy as np


def read_vtu(path, fields):
    s = open(path).read()

    def arr(name, dtype=float):
        m = re.search(r'Name="%s"[^>]*>([^<]*)<' % name, s)
        return (np.fromstring(m.group(1), sep=" ", dtype=dtype)
                if m else None)

    pts = re.search(r'<Points>.*?format="ascii">([^<]*)<', s, re.S)
    P = np.fromstring(pts.group(1), sep=" ").reshape(-1, 3)
    conn = arr("connectivity").astype(int)
    data = {f: arr(f) for f in fields}
    return P, conn, data


def main():
    ap = argparse.ArgumentParser(description="metriques de cratere fdem3d")
    ap.add_argument("run", help="dossier de sortie du run")
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--center", type=float, nargs=2, default=None)
    ap.add_argument("--dth", type=float, default=0.05)
    ap.add_argument("--sectors", type=int, default=12)
    ap.add_argument("--skin", type=float, default=None)
    ap.add_argument("--plot", default=None)
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()

    jfiles = sorted(glob.glob(os.path.join(a.run, "fdem3d_joints_*.vtu")))
    efiles = sorted(glob.glob(os.path.join(a.run, "fdem3d_[0-9]*.vtu")))
    if not jfiles:
        sys.exit("aucun fdem3d_joints_*.vtu dans " + a.run)
    jf = jfiles[a.frame]
    ef = efiles[a.frame] if efiles else None

    Pj, cj, dj = read_vtu(jf, ["damage", "tBreak"])
    # elements (corps, fragments, volumes) — lus d'abord : la SURFACE et le
    # detache se rapportent au CORPS IMPACTE (multi-corps V1 : le VTU joints
    # melange les joints internes de tous les corps ; sans ce filtre, la
    # surface etait prise au sommet de l'INSERT et l'insert entier compte
    # comme fragment detache — mesure sur bench1)
    Pe = TE = CE = V = frag = grain = None
    if ef:
        Pe, ce, de = read_vtu(ef, ["fragment", "grain"])
        TE = ce.reshape(-1, 4)
        CE = Pe[TE].mean(axis=1)
        V = np.abs(np.einsum(
            "ij,ij->i",
            np.cross(Pe[TE[:, 1]] - Pe[TE[:, 0]], Pe[TE[:, 2]] - Pe[TE[:, 0]]),
            Pe[TE[:, 3]] - Pe[TE[:, 0]])) / 6.0
        frag = de["fragment"]
        grain = de["grain"]
    TJ = cj.reshape(-1, 3)
    CJ = Pj[TJ].mean(axis=1)                       # centroides des joints
    e1 = Pj[TJ[:, 1]] - Pj[TJ[:, 0]]
    e2 = Pj[TJ[:, 2]] - Pj[TJ[:, 0]]
    AJ = 0.5 * np.linalg.norm(np.cross(e1, e2), axis=1)
    dam = dj["damage"]
    tbk = dj["tBreak"]
    broken = tbk >= 0.0
    dmgd = dam > a.dth

    if not broken.any():
        print("[crater] aucun joint casse — pas de cratere")
        return

    # corps impacte = grain majoritaire des tets les plus proches des joints
    # casses ; surface libre = plus haut sommet de CE corps
    body = None
    if ef and grain is not None and len(np.unique(grain)) > 1:
        from scipy.spatial import cKDTree
        _, idx = cKDTree(CE).query(CJ[broken])
        gb = grain[idx].astype(int)
        body = np.bincount(gb).argmax()
        mask = grain == body
        zSurf = float(Pe[TE[mask]].reshape(-1, 3)[:, 2].max())
    else:
        zSurf = Pj[:, 2].max()
    hJ = np.sqrt(np.median(AJ))
    skin = a.skin if a.skin is not None else 2.0 * hJ

    # axe d'impact
    if a.center is not None:
        cx, cy = a.center
    else:
        w = AJ[broken]
        cx = float((CJ[broken, 0] * w).sum() / w.sum())
        cy = float((CJ[broken, 1] * w).sum() / w.sum())

    R = np.hypot(CJ[:, 0] - cx, CJ[:, 1] - cy)     # rayon cylindrique
    TH = np.arctan2(CJ[:, 1] - cy, CJ[:, 0] - cx)  # angle

    inSkin = CJ[:, 2] >= zSurf - skin
    bSkin = broken & inSkin
    rCr = float(np.percentile(R[bSkin], 95)) if bSkin.any() else 0.0
    rMax = float(R[broken].max())
    depth = float(zSurf - CJ[broken, 2].min())
    aBrk = float(AJ[broken].sum())

    # volume endommage : elements adjacents a un joint D > dth. Sans table
    # joint->elements dans le VTU, l'adjacence se prend par proximite : un
    # tet est endommage si son centroide est a moins de h_tet du centroide
    # d'un joint endommage (approximation de post-traitement, honnete pour
    # une metrique de volume ; la table exacte viendra avec les besoins V3).
    vDam = vDet = 0.0
    if ef:
        hT = np.cbrt(np.median(V))
        inBody = (grain == body) if body is not None else np.ones(len(V), bool)
        if dmgd.any():
            from scipy.spatial import cKDTree
            tree = cKDTree(CJ[dmgd])
            near = tree.query_ball_point(CE, 1.5 * hT)
            mask = np.array([len(n) > 0 for n in near]) & inBody
            vDam = float(V[mask].sum())
        # detache : fragments != fragment majoritaire DU CORPS IMPACTE (les
        # autres corps — l'insert — ne sont pas des debris)
        fb = frag[inBody].astype(int)
        fMain = np.bincount(fb).argmax()
        vDet = float(V[inBody & (frag != fMain)].sum())

    # fissures radiales par secteur : portee max des joints casses AU-DELA
    # du rayon du cratere
    edges = np.linspace(-np.pi, np.pi, a.sectors + 1)
    reach = np.zeros(a.sectors)
    reachD = np.zeros(a.sectors)           # portee de l'ENDOMMAGEMENT :
    for k in range(a.sectors):             # les bras radiaux naissent
        m = broken & (TH >= edges[k]) & (TH < edges[k + 1]) & (R > rCr)
        reach[k] = float(R[m].max() - rCr) if m.any() else 0.0
        md = dmgd & (TH >= edges[k]) & (TH < edges[k + 1]) & (R > rCr)
        reachD[k] = float(R[md].max() - rCr) if md.any() else 0.0
    nRad = int((reach > 0).sum())
    nRadD = int((reachD > 0).sum())

    print(f"[crater] frame {jf.split('_')[-1].split('.')[0]}, "
          f"{int(broken.sum())} joints casses, {int(dmgd.sum())} endommages "
          f"(D > {a.dth})")
    print(f"[crater] centre d'impact : ({cx:.4f}, {cy:.4f}) m ; surface z = "
          f"{zSurf:.4f} m"
          + (f" (corps {body})" if body is not None else "")
          + f" ; peau {skin * 1e3:.2f} mm")
    print(f"[crater] R_crater (p95 peau)   = {rCr * 1e3:.2f} mm")
    print(f"[crater] R_max (tous casses)   = {rMax * 1e3:.2f} mm")
    print(f"[crater] profondeur            = {depth * 1e3:.2f} mm")
    print(f"[crater] aire cassee           = {aBrk * 1e6:.2f} mm^2")
    if ef:
        print(f"[crater] volume endommage      = {vDam * 1e9:.1f} mm^3")
        print(f"[crater] volume detache        = {vDet * 1e9:.1f} mm^3")
    print(f"[crater] fissures radiales     : {nRad}/{a.sectors} secteurs, "
          f"portee moy {reach[reach > 0].mean() * 1e3 if nRad else 0.0:.2f} mm, "
          f"max {reach.max() * 1e3:.2f} mm")
    print(f"[crater] bras endommages (D>{a.dth}) : {nRadD}/{a.sectors} "
          f"secteurs, portee moy "
          f"{reachD[reachD > 0].mean() * 1e3 if nRadD else 0.0:.2f} mm, "
          f"max {reachD.max() * 1e3:.2f} mm")

    if a.csv:
        with open(a.csv, "w") as f:
            f.write("metric,value,unit\n")
            f.write(f"nBroken,{int(broken.sum())},-\n")
            f.write(f"rCrater,{rCr},m\nrMax,{rMax},m\ndepth,{depth},m\n")
            f.write(f"aBroken,{aBrk},m2\nvDamaged,{vDam},m3\n")
            f.write(f"vDetached,{vDet},m3\n")
            f.write(f"radialSectors,{nRad},-\n")
            f.write(f"radialReachMax,{reach.max()},m\n")
            f.write(f"damagedSectors,{nRadD},-\n")
            f.write(f"damagedReachMax,{reachD.max()},m\n")
        print("[crater] csv ->", a.csv)

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8.5, 8))
        m = dmgd & ~broken
        ax.scatter((CJ[m, 0] - cx) * 1e3, (CJ[m, 1] - cy) * 1e3,
                   s=12, c="#f5a623", alpha=0.6,
                   label=f"endommagé (D>{a.dth})")
        sc = ax.scatter((CJ[broken, 0] - cx) * 1e3, (CJ[broken, 1] - cy) * 1e3,
                        s=30, c=CJ[broken, 2] * 1e3, cmap="viridis_r",
                        edgecolors="k", linewidths=0.4, label="cassé")
        fig.colorbar(sc, ax=ax, label="z du joint cassé [mm]", shrink=0.75)
        th = np.linspace(0, 2 * np.pi, 200)
        ax.plot(rCr * 1e3 * np.cos(th), rCr * 1e3 * np.sin(th),
                "r--", lw=1.5, label=f"R_crater = {rCr * 1e3:.1f} mm")
        for k in range(a.sectors):
            if reach[k] <= 0:
                continue
            am = 0.5 * (edges[k] + edges[k + 1])
            ax.annotate("", xy=(((rCr + reach[k]) * np.cos(am)) * 1e3,
                                ((rCr + reach[k]) * np.sin(am)) * 1e3),
                        xytext=((rCr * np.cos(am)) * 1e3,
                                (rCr * np.sin(am)) * 1e3),
                        arrowprops=dict(arrowstyle="->", color="#c0392b",
                                        lw=1.6))
        ax.set_xlabel("x − x₀ [mm]")
        ax.set_ylabel("y − y₀ [mm]")
        ax.set_aspect("equal")
        ax.legend(loc="upper right", fontsize=9)
        ax.set_title("Cratère vu de dessus — joints cassés (couleur = "
                     "profondeur), flèches = fissures radiales\n"
                     f"R = {rCr * 1e3:.1f} mm, prof. {depth * 1e3:.1f} mm, "
                     f"aire cassée {aBrk * 1e6:.0f} mm², "
                     f"{nRad} secteurs radiaux")
        fig.tight_layout()
        fig.savefig(a.plot, dpi=140)
        print("[crater] figure ->", a.plot)


if __name__ == "__main__":
    main()
