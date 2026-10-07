#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# edz_sectors.py — la FORME de l EDZ par secteur d azimut, corrigee de la paroi.
#
#   python tunnel_schisto/tools/edz_sectors.py out_A [out_B ...] [--dip 45]
#          [--cx 50 --cy 50] [--zone 6 18]
#
# POURQUOI. edz_metrics.py (p1) mesure des demi-axes HORIZONTAL / VERTICAL :
# il est AVEUGLE a une ellipse a 45 deg (mesure du 2026-09-02 : h/v = 1,01 sur
# un losange a 45 deg dont le rapport le long / en travers vaut 1,77). Ici on
# mesure, par secteur de 30 deg, la PROFONDEUR d enveloppe rompue = rayon p95
# des joints rompus du secteur MOINS le rayon de paroi du secteur (le plus
# petit rayon de joint du secteur : c est la paroi du fer a cheval, qui n est
# pas un cercle). Puis le rapport profondeur le long du litage / en travers.
#
# Geometrie de REFERENCE (trame 0), pas deformee : une enveloppe mesuree sur
# la configuration courante melangerait forme de l EDZ et convergence.
# ---------------------------------------------------------------------------
import argparse
import glob
import re
import numpy as np


def vtu(f):
    raw = open(f, "rb").read().decode("utf8", "ignore")

    def arr(name, dtype=float):
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, raw, re.S)
        return np.fromstring(m.group(1).strip(), sep=" ", dtype=dtype) if m else None
    m = re.search(r'<Points>.*?<DataArray[^>]*>(.*?)</DataArray>', raw, re.S)
    P = np.fromstring(m.group(1).strip(), sep=" ").reshape(-1, 3)[:, :2]
    return P, arr


def sectors(out, cx, cy, dip, width=30):
    fr = sorted(glob.glob(f"{out}/fdem_joints_0*.vtu"))
    P0, a0 = vtu(fr[0])
    _, a = vtu(fr[-1])
    conn = a0("connectivity", int).reshape(-1, 2)
    seg = P0[conn]
    tb = a("tBreak")
    br = tb >= 0
    c = seg.mean(1)
    rad = np.hypot(c[:, 0] - cx, c[:, 1] - cy)
    az = np.degrees(np.arctan2(c[:, 1] - cy, c[:, 0] - cx)) % 360
    rows = []
    for a0_ in range(0, 360, width):
        m = (az >= a0_) & (az < a0_ + width)
        mb = m & br
        wall = rad[m].min()
        r95 = np.percentile(rad[mb], 95) if mb.sum() > 20 else np.nan
        rows.append((a0_, wall, r95, r95 - wall, int(mb.sum())))
    depth = np.array([r[3] for r in rows])
    centers = np.array([r[0] + width / 2 for r in rows])
    if dip is not None:
        # secteurs contenant les directions du litage et sa normale
        def pick(angle):
            d = np.abs(((centers - angle) + 180) % 360 - 180)
            return depth[d <= width / 2 + 1e-9]
        along = np.nanmean(np.concatenate([pick(dip), pick(dip + 180)]))
        across = np.nanmean(np.concatenate([pick(dip + 90), pick(dip + 270)]))
    else:
        along = across = np.nan
    return rows, along, across, int(br.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--dip", type=float, default=None, help="pendage du litage [deg]")
    ap.add_argument("--cx", type=float, default=50.0)
    ap.add_argument("--cy", type=float, default=50.0)
    a = ap.parse_args()
    for run in a.runs:
        rows, along, across, nb = sectors(run, a.cx, a.cy, a.dip)
        print(f"--- {run} : {nb} joints rompus ---")
        print("  azimut  paroi[m]  p95[m]  PROFONDEUR[m]  n")
        for a0_, wall, r95, d, n in rows:
            print(f"  {a0_:4d}    {wall:6.2f}  {r95:6.2f}     {d:6.2f}     {n}")
        if a.dip is not None:
            print(f"  le long du litage ({a.dip:.0f}/{a.dip+180:.0f} deg) : {along:.2f} m ; "
                  f"en travers ({a.dip+90:.0f}/{a.dip+270:.0f}) : {across:.2f} m ; "
                  f"RAPPORT = {along/across:.2f}")
            print("  (temoin isotrope out_tun_corr_th4, memes secteurs : 1,19 ; "
                  "plancher de bruit h/v : 1,35)")


if __name__ == "__main__":
    main()
