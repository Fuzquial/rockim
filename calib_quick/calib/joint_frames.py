#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# joint_frames.py — evolution des joints PAR FRAME depuis fdem_joints_XXXX.vtu :
# inseres (bonded = 0, insertion adaptative), endommages (damage > 0), D moyen
# des endommages, rompus (damage >= 1) — et le q correspondant via history.csv.
# Proxys : q_CI = q au premier frame ou des joints sont inseres/endommages
# (initiation de la microfissuration), q_CD = q ou le taux d endommagement
# par frame est maximal avant le pic (croissance instable).
#
#   python calib_quick/calib/joint_frames.py out_dir [out_dir ...]
# ---------------------------------------------------------------------------
import os
import re
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from extract import load_run  # noqa: E402


def read_joint_fields(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    out = {}
    for name in ("bonded", "damage"):
        m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
        out[name] = np.array(m.group(1).split(), dtype=float) if m else None
    return out


def frame_table(out_dir):
    frames = np.genfromtxt(os.path.join(out_dir, "frames.csv"), delimiter=",", names=True)
    r = load_run(out_dir)
    hist = np.genfromtxt(os.path.join(out_dir, "history.csv"), delimiter=",", names=True, invalid_raise=False)
    t_h = hist["t"]; ok = np.isfinite(t_h)
    rows = []
    for fr, t in zip(frames["frame"].astype(int), frames["t"]):
        p = os.path.join(out_dir, f"fdem_joints_{fr:04d}.vtu")
        if not os.path.exists(p):
            continue
        f = read_joint_fields(p)
        n = len(f["damage"]) if f["damage"] is not None else 0
        ins = int((f["bonded"] == 0).sum()) if f["bonded"] is not None else n
        dam = f["damage"] if f["damage"] is not None else np.zeros(0)
        nd = int((dam > 0).sum()); nb = int((dam >= 1).sum())
        dmean = float(dam[dam > 0].mean()) if nd else 0.0
        # q, eps au temps t (interpolation dans la courbe extraite)
        tl = t_h[ok]
        # eps/q de load_run sont indexes sur t > delay ; on interpole sur t
        t_load = tl[tl > float(r["cfg"].get("pullDelay", 0))]
        q = float(np.interp(t, t_load, r["q"])) if t > t_load[0] else 0.0
        eps = float(np.interp(t, t_load, r["eps"])) if t > t_load[0] else 0.0
        rows.append((fr, t, eps, q, n, ins, nd, dmean, nb))
    return rows, r


def proxies(rows):
    q = np.array([x[3] for x in rows]); ins = np.array([x[5] for x in rows]); nd = np.array([x[6] for x in rows])
    ipk = int(np.argmax(q)); qpk = q[ipk]
    act = np.maximum(ins, nd)
    i_ci = next((i for i in range(len(rows)) if act[i] > 0), None)
    q_ci_lo = q[i_ci - 1] if i_ci else float("nan"); q_ci_hi = q[i_ci] if i_ci is not None else float("nan")
    rate = np.diff(nd[:ipk + 1]) if ipk > 0 else np.array([0])
    i_cd = int(np.argmax(rate)) + 1 if len(rate) else None
    q_cd = q[i_cd] if i_cd is not None else float("nan")
    return dict(q_peak=qpk, q_CI_lo=q_ci_lo, q_CI_hi=q_ci_hi, CI_frac=(0.5 * (q_ci_lo + q_ci_hi) / qpk if i_ci else float("nan")),
                q_CD=q_cd, CD_frac=q_cd / qpk if i_cd is not None else float("nan"), n_act_peak=int(act[ipk]))


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    for d in sys.argv[1:]:
        rows, r = frame_table(d)
        print(f"== {d}  (sigma3 = {r['s3']:.0f} MPa)")
        print(f"{'frame':>5s} {'eps %':>6s} {'q MPa':>7s} {'joints':>6s} {'inseres':>7s} {'D>0':>6s} {'Dmoy':>5s} {'rompus':>6s}")
        for fr, t, eps, q, n, ins, nd, dmean, nb in rows:
            print(f"{fr:5d} {100*eps:6.2f} {q:7.1f} {n:6d} {ins:7d} {nd:6d} {dmean:5.2f} {nb:6d}")
        p = proxies(rows)
        print("  proxys : q_CI dans [%.0f ; %.0f] (%.2f du pic), q_CD = %.0f (%.2f du pic), actifs au pic %d" %
              (p["q_CI_lo"], p["q_CI_hi"], p["CI_frac"], p["q_CD"], p["CD_frac"], p["n_act_peak"]))


if __name__ == "__main__":
    main()
