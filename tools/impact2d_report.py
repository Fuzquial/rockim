#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# impact2d_report.py — dépouillement d'un impact 2D (mode fdem, insert
# analytique) : restitution, enfoncement, joints rompus, fragments, énergie.
#
#   python tools/impact2d_report.py out_dir [journal.log]
#
# Restitution e = (v_y max après le pic d'enfoncement) / |v_y initiale|.
# Cible Yang 9 m/s (bit) : e = 4,65 / 5,62 = 0,83.
# ---------------------------------------------------------------------------
import csv
import re
import sys


def main():
    out = sys.argv[1]
    log = sys.argv[2] if len(sys.argv) > 2 else None
    with open(out + "/history.csv", newline="") as f:
        rd = csv.reader(f)
        h = next(rd)
        d = [[float(x) for x in r] for r in rd if r]
    c = {k: i for i, k in enumerate(h)}
    t = [r[0] for r in d]
    vy = [r[c["toolVy"]] for r in d]
    y = [r[c["toolY"]] for r in d]
    v0 = abs(vy[0])
    imin = min(range(len(y)), key=lambda i: y[i])       # enfoncement max
    pen = (y[0] - y[imin]) * 1e3
    after = range(imin, len(vy))
    imax = max(after, key=lambda i: vy[i]) if len(after) else imin
    e = vy[imax] / v0 if v0 > 0 else float("nan")
    print("=== %s : t final %.1f us, %d lignes ===" % (out, t[-1] * 1e6, len(d)))
    print("enfoncement max        : %.3f mm a %.1f us" % (pen, t[imin] * 1e6))
    print("v rebond max           : %.2f m/s a %.1f us  -> e = %.3f  [Yang : 0,83]"
          % (vy[imax], t[imax] * 1e6, e))
    print("outil : KE %.1f -> %.1f J/m ; travail outil->solide %.1f J/m"
          % (d[0][c["toolKE"]], d[-1][c["toolKE"]], d[-1][c["work"]]))
    for k in ("nBroken", "nFrag", "detachedVol", "nPulv"):
        if k in c:
            print("%-22s : %g" % (k, d[-1][c[k]]))
    print("--- energie (J/m, history) ---")
    for k in ("eEl", "eJnt", "eGc", "eFric", "eCund", "eLys", "bdWork"):
        if k in c:
            s = [r[c[k]] for r in d]
            print("   %-8s final %10.2f  min %10.2f  max %10.2f" % (k, s[-1], min(s), max(s)))
    if log:
        txt = open(log, errors="replace").read()
        for pat in (r"dt = [0-9.e+-]+ s, steps = \d+", r"breakage mode:.*",
                    r"residu\s*:.*", r"joints\s*: .*cohesif.*", r"contact\s*: .*",
                    r"elements\s*: .*", r"energy budget.*", r"ENERGY ABORT.*",
                    r"broken joints\s*:.*", r"wall time:.*"):
            for m in re.finditer(pat, txt):
                print("   " + m.group(0).strip()[:170])
    return 0


if __name__ == "__main__":
    sys.exit(main())
