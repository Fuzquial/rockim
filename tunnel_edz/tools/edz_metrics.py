#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# edz_metrics.py — depouillement d'un run tunnel : zone endommagee (EDZ),
# classement des fissures et convergence de la paroi.
#
#   python tunnel_edz/tools/edz_metrics.py out_tun_s5 [--cx 50 --cy 50]
#                                          [--qmix 0.5] [--json fichier.json]
#
# Entrees (aucune sortie nouvelle du solveur n'est necessaire) :
#   <run>/fdem_final_joints.csv   x1,y1,x2,y2,damage,type,breakMode,rn,rs,
#                                 tBreak,bonded
#   <run>/fdem_0000.vtu et la DERNIERE trame : positions nodales (deformees)
#                                 -> deplacements par difference
#
# Observables produits, dans l'ordre des figures de Wang et al. (2024) :
#   fig. 8/13/16  nombre de fissures par mode (traction / mixte / cisaillement)
#   fig. 12f      rayon de l'EDZ (p95 et max des distances au centre)
#   fig. 15f      demi-axes vertical et horizontal de l'EDZ (balayage lambda)
#   fig. 13b/16b  longueur cumulee et aire cumulee des fissures
#   fig. 11/14a   deplacement maximal et convergences caracteristiques
#
# CLASSEMENT DES MODES. Le solveur stampe deja `breakMode` (1 = traction si
# rn >= rs, 2 = cisaillement sinon) : c'est un partage BINAIRE. L'article
# compte trois classes ; son critere n'est pas publie. On prend donc le rapport
# q = min(rn, rs) / max(rn, rs) : mixte si q >= qmix (defaut 0,5), sinon pur.
# Les DEUX lectures sont imprimees, pour que la comparaison reste honnete.
# ---------------------------------------------------------------------------
import argparse
import glob
import json
import os
import re

import numpy as np


def read_joints(run):
    path = os.path.join(run, "fdem_final_joints.csv")
    if not os.path.exists(path):
        raise SystemExit(f"introuvable : {path}")
    d = np.genfromtxt(path, delimiter=",", names=True, invalid_raise=False)
    return d


def vtu_points(path):
    """Positions nodales (deformees) d'une trame VTU ASCII de rockim."""
    with open(path) as f:
        txt = f.read()
    m = re.search(r'<Points>.*?<DataArray[^>]*>(.*?)</DataArray>', txt,
                  re.S)
    if not m:
        raise SystemExit(f"pas de bloc <Points> dans {path}")
    a = np.fromstring(m.group(1), sep=" ")
    return a.reshape(-1, 3)[:, :2]


def displacement(run):
    frames = sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
    frames = [f for f in frames if "joints" not in os.path.basename(f)]
    if len(frames) < 2:
        return None, None, None
    P0, P1 = vtu_points(frames[0]), vtu_points(frames[-1])
    if P0.shape != P1.shape:
        print("[avert] trames de tailles differentes — deplacement ignore")
        return None, None, None
    U = P1 - P0
    return P0, np.linalg.norm(U, axis=1), U


def compute(run, cx=50.0, cy=50.0, qmix=0.5, thickness=1.0):
    """Les memes metriques, utilisables comme fonction (cf. plot_sweep.py)."""
    class A:
        pass
    a = A()
    a.run, a.cx, a.cy, a.qmix, a.thickness, a.json = (run, cx, cy, qmix,
                                                      thickness, None)
    return _run(a, quiet=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--cx", type=float, default=50.0)
    ap.add_argument("--cy", type=float, default=50.0)
    ap.add_argument("--qmix", type=float, default=0.5)
    ap.add_argument("--thickness", type=float, default=1.0)
    ap.add_argument("--json", default=None)
    _run(ap.parse_args())


def _run(a, quiet=False):
    J = read_joints(a.run)
    brk = (J["tBreak"] >= 0.0) | (J["damage"] >= 0.999)
    n = int(brk.sum())
    res = {"run": os.path.basename(a.run), "joints": int(J.size),
           "broken": n}
    if n == 0:
        if not quiet:
            print(f"{res['run']} : AUCUN joint casse sur {J.size}")
        return res

    x1, y1 = J["x1"][brk], J["y1"][brk]
    x2, y2 = J["x2"][brk], J["y2"][brk]
    L = np.hypot(x2 - x1, y2 - y1)
    xm, ym = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
    r = np.hypot(xm - a.cx, ym - a.cy)

    rn, rs = J["rn"][brk], J["rs"][brk]
    mx = np.maximum(rn, rs)
    q = np.where(mx > 0, np.minimum(rn, rs) / np.maximum(mx, 1e-300), 0.0)
    mixed = q >= a.qmix
    tens = (~mixed) & (rn >= rs)
    shear = (~mixed) & (rn < rs)
    bm = J["breakMode"][brk]

    res.update(
        tensile=int(tens.sum()), mixed=int(mixed.sum()),
        shear=int(shear.sum()),
        binary_tensile=int((bm == 1).sum()), binary_shear=int((bm == 2).sum()),
        crack_length_m=float(L.sum()),
        crack_area_m2=float(L.sum() * a.thickness),
        edz_radius_p95_m=float(np.percentile(r, 95)),
        edz_radius_max_m=float(r.max()),
        edz_halfaxis_x_p95_m=float(np.percentile(np.abs(xm - a.cx), 95)),
        edz_halfaxis_y_p95_m=float(np.percentile(np.abs(ym - a.cy), 95)))

    P0, umag, U = displacement(a.run)
    if umag is not None:
        res["u_max_m"] = float(umag.max())
        d0 = np.hypot(P0[:, 0] - a.cx, P0[:, 1] - a.cy)
        wall = d0 < 8.0                       # noeuds de la paroi et proches
        if wall.any():
            top = wall & (P0[:, 1] > a.cy + 3.0)
            side = wall & (np.abs(P0[:, 0] - a.cx) > 4.5)
            bot = wall & (P0[:, 1] < a.cy - 3.5)
            for nm, sel in (("u_crown_m", top), ("u_sidewall_m", side),
                            ("u_invert_m", bot)):
                res[nm] = float(umag[sel].max()) if sel.any() else None

    if quiet:
        return res
    print(f"--- {res['run']} : {n} joints casses / {J.size} ---")
    print(f"  modes (q >= {a.qmix}) : traction {res['tensile']}, "
          f"mixte {res['mixed']}, cisaillement {res['shear']}")
    print(f"  lecture binaire du solveur : traction {res['binary_tensile']}, "
          f"cisaillement {res['binary_shear']}")
    print(f"  EDZ : rayon p95 {res['edz_radius_p95_m']:.2f} m, "
          f"max {res['edz_radius_max_m']:.2f} m ; demi-axes "
          f"h {res['edz_halfaxis_x_p95_m']:.2f} / "
          f"v {res['edz_halfaxis_y_p95_m']:.2f} m")
    print(f"  fissures : longueur {res['crack_length_m']:.1f} m, "
          f"aire {res['crack_area_m2']:.1f} m2")
    if "u_max_m" in res:
        print(f"  deplacement max {res['u_max_m']:.4f} m "
              f"(couronne {res.get('u_crown_m') or float('nan'):.4f}, "
              f"reins {res.get('u_sidewall_m') or float('nan'):.4f}, "
              f"radier {res.get('u_invert_m') or float('nan'):.4f})")
    if a.json:
        with open(a.json, "w") as f:
            json.dump(res, f, indent=1)
        print("  ->", a.json)
    return res


if __name__ == "__main__":
    main()
