#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# bench_compare.py — plusieurs runs d'impact (meme deck de base, variantes)
# cote a cote : physique a un instant commun + cout par pas.
#
#   python tools/bench_compare.py --t 100e-6 out_A out_B out_C ...
#
# Pour chaque run : history.csv (physique a l'instant --t, ou au dernier
# instant commun s'il est plus court) et, s'il existe, results/<nom>.log
# (dt, pas, [prof3d], potential stats, wall time, residu). Le nom du journal
# est deduit du dossier : out_X -> results/X.log ; --log peut le surcharger
# par run sous la forme dossier=journal.
# ---------------------------------------------------------------------------
import argparse
import csv
import os
import re


def load(run):
    with open(os.path.join(run, "history.csv"), newline="") as f:
        rd = csv.reader(f)
        head = next(rd)
        rows = [[float(x) for x in r] for r in rd if r]
    return head, rows


def at(rows, col, t):
    # derniere ligne dont t <= t demande
    best = None
    for r in rows:
        if r[0] <= t:
            best = r
        else:
            break
    return best


def first_positive(rows, j):
    for r in rows:
        if r[j] > 0:
            return r[0]
    return None


def logfacts(path):
    d = {}
    if not path or not os.path.exists(path):
        return d
    txt = open(path, errors="replace").read()
    m = re.search(r"dt = ([0-9.e+-]+) s, steps = (\d+)", txt)
    if m:
        d["dt_ns"] = float(m.group(1)) * 1e9
        d["steps"] = int(m.group(2))
    m = re.search(r"per step \(ms\): elem ([0-9.]+) insert ([0-9.]+) joint ([0-9.]+) "
                  r"gcontact ([0-9.]+) tool ([0-9.]+)\s+\((\d+) steps\)", txt)
    if m:
        d["ms_elem"], d["ms_ins"], d["ms_joint"], d["ms_gc"], d["ms_tool"] = \
            [float(m.group(i)) for i in range(1, 6)]
        d["ms_total"] = sum(float(m.group(i)) for i in range(1, 6))
    m = re.search(r"potential stats: (\d+) paires .*?(\d+) clip-force\), tGrid = ([0-9.]+) s, "
                  r"tLoop = ([0-9.]+) s", txt)
    if m:
        d["pairs"], d["clipF"] = int(m.group(1)), int(m.group(2))
        d["tGrid"], d["tLoop"] = float(m.group(3)), float(m.group(4))
    m = re.search(r"wall time: ([0-9.]+) s", txt)
    if m:
        d["wall_s"] = float(m.group(1))
    m = re.search(r"residu\s*:\s*([-0-9.e+]+) J \(([-0-9.e+]+) %", txt)
    if m:
        d["residu_pct"] = float(m.group(2))
    d["abort"] = "ABORT" in txt or "budget" in txt and "depasse" in txt
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--t", type=float, default=100e-6, help="instant de lecture [s]")
    ap.add_argument("--log", nargs="*", default=[], help="dossier=journal")
    a = ap.parse_args()
    logmap = dict(x.split("=", 1) for x in a.log)
    cols = ["t_us", "szz_min_MPa", "t_szz_min_us", "t_1er_rompu_us", "nBroken", "nPulv",
            "pen_mm", "vz_bit", "vz_insert", "vz_piston", "eEl", "eJnt", "eGc", "eFric",
            "bdWork", "dt_ns", "steps", "ms_total", "ms_gc", "ms_joint", "wall_s", "residu_pct"]
    table = {}
    for run in a.runs:
        head, rows = load(run)
        c = {k: i for i, k in enumerate(head)}
        r = at(rows, 0, a.t) or rows[-1]
        szz = [x[c["szz_bit"]] for x in rows] if "szz_bit" in c else [0.0]
        i = min(range(len(szz)), key=lambda k: szz[k])
        d = {
            "t_us": r[0] * 1e6,
            "szz_min_MPa": -szz[i] / 1e6,
            "t_szz_min_us": rows[i][0] * 1e6,
            "t_1er_rompu_us": (first_positive(rows, c["nBroken"]) or 0) * 1e6,
            "nBroken": r[c["nBroken"]],
            "nPulv": r[c["nPulv"]] if "nPulv" in c else float("nan"),
            "pen_mm": (rows[0][c["z_insert"]] - min(x[c["z_insert"]] for x in rows if x[0] <= r[0])) * 1e3
                      if "z_insert" in c else float("nan"),
            "vz_bit": r[c["vz_bit"]] if "vz_bit" in c else float("nan"),
            "vz_insert": r[c["vz_insert"]] if "vz_insert" in c else float("nan"),
            "vz_piston": r[c["vz_piston"]] if "vz_piston" in c else float("nan"),
            "eEl": r[c["eEl"]], "eJnt": r[c["eJnt"]], "eGc": r[c["eGc"]],
            "eFric": r[c["eFric"]],
            "bdWork": r[c["bdWork"]] if "bdWork" in c else float("nan"),
        }
        name = os.path.basename(run.rstrip("/\\"))
        log = logmap.get(run) or logmap.get(name) or os.path.join(
            "results", name[4:] + ".log" if name.startswith("out_") else name + ".log")
        d.update(logfacts(log))
        table[name] = d
    names = list(table)
    w = max(len(n) for n in names) + 2
    print("lecture a t = %.1f us (ou dernier instant du run)" % (a.t * 1e6))
    print("%-16s" % "" + "".join("%*s" % (w, n[-w + 2:]) for n in names))
    for k in cols:
        line = "%-16s" % k
        for n in names:
            v = table[n].get(k)
            if v is None:
                line += "%*s" % (w, "-")
            elif isinstance(v, float):
                line += "%*s" % (w, ("%.3g" % v) if abs(v) < 1e4 else ("%.4g" % v))
            else:
                line += "%*s" % (w, v)
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
