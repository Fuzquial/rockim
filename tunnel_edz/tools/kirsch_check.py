#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# kirsch_check.py — LE controle quantitatif des patchs 1 et 2.
#
#   python tunnel_edz/tools/kirsch_check.py out_kirsch_lam1 \
#          [tunnel_edz/configs/verif_kirsch_lam1.cfg]
#
# Autour d'un trou circulaire de rayon a dans une plaque infinie soumise aux
# contraintes lointaines Sxx, Syy (convention TENSION POSITIVE, donc negatives
# ici), la contrainte orthoradiale en paroi vaut, exactement :
#
#     sigma_theta(theta) = (Sxx + Syy) - 2 (Sxx - Syy) cos 2theta
#
# soit, avec sigma_h et sigma_v en COMPRESSION (Sxx = -sigma_h, Syy = -sigma_v),
# une compression de (3 sigma_v - sigma_h) aux reins (theta = 0) et de
# (3 sigma_h - sigma_v) a la couronne (theta = 90 deg).
#
# Le script lit la DERNIERE trame VTU, garde le premier anneau d'elements
# autour du trou, projette leur tenseur sur la direction orthoradiale et
# compare. Tolerance par defaut 5 % : Kirsch vaut pour une plaque INFINIE,
# ici a/b = 1/10.
# ---------------------------------------------------------------------------
import glob
import os
import re
import sys

import numpy as np

TOL = 0.05


def vtu_arrays(path):
    with open(path) as f:
        txt = f.read()
    pts = re.search(r'<Points>.*?<DataArray[^>]*>(.*?)</DataArray>', txt, re.S)
    if not pts:
        raise SystemExit(f"pas de <Points> dans {path}")
    P = np.fromstring(pts.group(1), sep=" ").reshape(-1, 3)[:, :2]
    conn = re.search(r'Name="connectivity"[^>]*>(.*?)</DataArray>', txt, re.S)
    C = np.fromstring(conn.group(1), sep=" ", dtype=float).astype(int)
    C = C.reshape(-1, 3)
    out = {}
    for nm in ("sigmaXX", "sigmaYY", "sigmaXY"):
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % nm, txt, re.S)
        if not m:
            raise SystemExit(f"champ {nm} absent de {path} — le run a-t-il "
                             "ete fait avec le bon mode ?")
        out[nm] = np.fromstring(m.group(1), sep=" ")
    return P, C, out


def cfg_get(path, key, default):
    if not path or not os.path.exists(path):
        return default
    with open(path) as f:
        for line in f:
            line = line.split("#")[0]
            if "=" in line and line.split("=")[0].strip() == key:
                return float(line.split("=")[1].strip())
    return default


def main():
    run = sys.argv[1]
    cfg = sys.argv[2] if len(sys.argv) > 2 else None
    sh = cfg_get(cfg, "insituSh", 5e6)
    sv = cfg_get(cfg, "insituSv", 5e6)
    cx = cfg_get(cfg, "boreCX", 50.0)
    cy = cfg_get(cfg, "boreCY", 50.0)
    R = float(os.environ.get("KIRSCH_R", 5.0))

    frames = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
              if "joints" not in os.path.basename(f)]
    if not frames:
        raise SystemExit(f"aucune trame VTU dans {run}")
    P, C, S = vtu_arrays(frames[-1])

    ctr = P[C].mean(axis=1)
    x, y = ctr[:, 0] - cx, ctr[:, 1] - cy
    r = np.hypot(x, y)
    area = 0.5 * np.abs(
        (P[C[:, 1], 0] - P[C[:, 0], 0]) * (P[C[:, 2], 1] - P[C[:, 0], 1])
        - (P[C[:, 1], 1] - P[C[:, 0], 1]) * (P[C[:, 2], 0] - P[C[:, 0], 0]))
    # Taille d'element LOCALE (pres du trou), pas la moyenne du maillage : sur
    # un maillage gradue la moyenne est celle du champ lointain, et un anneau
    # epais de 2 h_moyen monte a r = 1,4 a, ou Kirsch predit deja 15 % de moins.
    near = r < 1.5 * R
    hLoc = np.sqrt(2.0 * np.median(area[near])) if near.any() else 0.1
    ring = (r > R) & (r < R + 3.0 * hLoc)
    if ring.sum() < 20:
        ring = (r > R) & (r < R + 6.0 * hLoc)
    th = np.arctan2(y[ring], x[ring])
    rr = r[ring]
    t = np.stack([-np.sin(th), np.cos(th)], axis=1)     # direction orthoradiale
    nn = np.stack([np.cos(th), np.sin(th)], axis=1)     # direction radiale
    sxx, syy, sxy = S["sigmaXX"][ring], S["sigmaYY"][ring], S["sigmaXY"][ring]

    def proj(a, b):
        return (a[:, 0] * (sxx * b[:, 0] + sxy * b[:, 1])
                + a[:, 1] * (sxy * b[:, 0] + syy * b[:, 1]))

    sth, srr = proj(t, t), proj(nn, nn)

    # Kirsch COMPLET (r-dependant) : on compare chaque element a la theorie a
    # SON rayon, ce qui supprime le biais d'epaisseur d'anneau.
    Sxx, Syy = -sh, -sv                                  # tension positive
    a2 = (R / rr) ** 2
    a4 = a2 * a2
    m, d = 0.5 * (Sxx + Syy), 0.5 * (Sxx - Syy)
    theo_r = m * (1.0 + a2) - d * (1.0 + 3.0 * a4) * np.cos(2 * th)
    theo_rr = m * (1.0 - a2) + d * (1.0 - 4.0 * a2 + 3.0 * a4) * np.cos(2 * th)

    print(f"--- Kirsch : {os.path.basename(run)} "
          f"(sigma_h = {sh/1e6:g}, sigma_v = {sv/1e6:g} MPa, R = {R:g} m) ---")
    print(f"  anneau : {ring.sum()} elements, r/a de {rr.min()/R:.3f} a "
          f"{rr.max()/R:.3f} (h local {hLoc:.3f} m)")
    print("  theta [deg]   sigma_theta [MPa]   Kirsch(r) [MPa]   ecart   "
          "| sigma_rr [MPa] (0 en paroi)")
    worst = 0.0
    for lo in range(0, 180, 15):
        sel = (np.degrees(th) % 180 >= lo) & (np.degrees(th) % 180 < lo + 15)
        if sel.sum() < 3:
            continue
        meas, theo = sth[sel].mean() / 1e6, theo_r[sel].mean() / 1e6
        err = abs(meas - theo) / max(abs(theo), 1e-9)
        worst = max(worst, err)
        print(f"   {lo:3d}-{lo+15:3d}        {meas:9.3f}         {theo:9.3f}   "
              f"{100*err:6.1f} %   |   {srr[sel].mean()/1e6:7.3f} "
              f"(theo {theo_rr[sel].mean()/1e6:6.3f})")
    # extrapolation a la paroi : chaque element est ramene a r = a par le
    # rapport des valeurs theoriques, ce qui donne un sigma_theta(a) mesure
    theo_wall = 2.0 * m - 4.0 * d * np.cos(2 * th)
    est = sth * theo_wall / np.where(np.abs(theo_r) > 1e-9, theo_r, 1e-9)
    for nm, sel in (("couronne (theta = 90 deg)",
                     np.abs(np.abs(np.degrees(th)) - 90) < 15),
                    ("reins (theta = 0/180 deg)",
                     (np.abs(np.degrees(th)) < 15)
                     | (np.abs(np.abs(np.degrees(th)) - 180) < 15))):
        if sel.any():
            tw = (2.0 * m - 4.0 * d * np.cos(2 * th[sel])).mean() / 1e6
            print(f"  paroi extrapolee, {nm:26s} : {est[sel].mean()/1e6:8.3f} "
                  f"MPa  (theorie {tw:8.3f})")
    print(f"  ecart maximal : {100*worst:.1f} %  "
          f"[{'PASS' if worst <= TOL else 'FAIL'} a {100*TOL:.0f} %]")
    if worst > TOL:
        print("  pistes : signe de insituS_ (patch 1) ; convention de normale "
              "du relachement (patch 2) ; run pas encore relaxe (allonger T ou "
              "monter dampingLocal) ; anneau trop epais (maillage grossier).")


if __name__ == "__main__":
    main()
