#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# law_matrix.py — la matrice FRACTURE x CONSERVATION du conseil du 12/09 (N4).
#
#   python tools/law_matrix.py out_A out_B ... [--t 200e-6]
#
# Une ligne par run, a l'instant --t (ou au dernier instant ecrit) :
#   fracture     : nBroken, nFrag, nPulv, enfoncement de l insert (mm)
#   conservation : poste joints (J, NEGATIF = les joints donnent de l energie),
#                  correction leapfrog (J), abort KE (t, exces), residu B4
# Les postes viennent du journal <out>.log range dans results/ (nom du dossier
# sans le prefixe out_), ou de --log.
# ---------------------------------------------------------------------------
import argparse
import csv
import os
import re
import sys


def last_row(out, t):
    with open(os.path.join(out, "history.csv"), newline="") as f:
        rd = csv.DictReader(f)
        rows = [r for r in rd if r]
    if not rows:
        return None
    if t is None:
        return rows[-1]
    best = rows[0]
    for r in rows:
        if float(r["t"]) <= t:
            best = r
        else:
            break
    return best


def posts(log):
    d = dict(joints=float("nan"), leap=float("nan"), abort="", res="")
    if not log or not os.path.isfile(log):
        return d
    txt = open(log, errors="replace").read()
    m = re.search(r"joints\s*:\s*(-?[\d.eE+-]+) J cohesif", txt)
    if m:
        d["joints"] = float(m.group(1))
    m = re.search(r"integration\s*:\s*\+?(-?[\d.eE+-]+) J", txt)
    if m:
        d["leap"] = float(m.group(1))
    m = re.search(r"ENERGY ABORT.*?a t = ([\d.eE+-]+) s.*?exces ([\d.eE+-]+) J", txt)
    if m:
        d["abort"] = "%.1f us, %.1f J" % (float(m.group(1)) * 1e6, float(m.group(2)))
    m = re.search(r"residu\s*:\s*(-?[\d.eE+-]+) J \(([\d.eE+-]+) %", txt)
    if m:
        d["res"] = "%.0e %%" % float(m.group(2))
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--t", type=float, default=None, help="instant de lecture (s)")
    a = ap.parse_args()
    print("%-34s %8s %7s %6s %6s %8s | %9s %8s %-16s %s" % (
        "run", "t(us)", "rompus", "frag", "pulv", "enf(mm)", "joints(J)", "leap(J)", "abort", "residu"))
    for out in a.runs:
        r = last_row(out, a.t)
        name = os.path.basename(out.rstrip("/\\"))
        log = os.path.join("results", name[4:] + ".log") if name.startswith("out_") else None
        p = posts(log)
        if r is None:
            print("%-34s (vide)" % name)
            continue
        zi = None
        try:
            with open(os.path.join(out, "history.csv"), newline="") as f:
                z0 = float(next(csv.DictReader(f))["z_insert"])
            zi = (z0 - float(r["z_insert"])) * 1e3
        except Exception:
            pass
        print("%-34s %8.1f %7s %6s %6s %8s | %9.2f %8.2f %-16s %s" % (
            name, float(r["t"]) * 1e6, r.get("nBroken", "?"), r.get("nFrag", "?"),
            r.get("nPulv", "?"), ("%.2f" % zi) if zi is not None else "?",
            p["joints"], p["leap"], p["abort"], p["res"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
