# -*- coding: utf-8 -*-
"""Test n2 de la revue : histogramme de D des joints INSERES.
Un mode pique a D << 1 = tapis d insertions avortees ; bimodal = localisation."""
import glob
import io
import re
import sys

import numpy as np

for run in sys.argv[1:]:
    fs = sorted(glob.glob(run + "/fdem_joints_[0-9]*.vtu"))
    if not fs:
        print(run, ": aucun VTU joints")
        continue
    s = io.open(fs[-1], encoding="utf-8", errors="ignore").read()
    D = np.fromstring(re.search(
        r'Name="damage"[^>]*>\s*(.*?)\s*</DataArray>', s, re.S).group(1),
        sep=" ")
    n = len(D)
    q = [(D < 0.05).sum(), ((D >= 0.05) & (D < 0.5)).sum(),
         ((D >= 0.5) & (D < 0.999)).sum(), (D >= 0.999).sum()]
    print("%-22s %s  (%d joints traces)" % (run, fs[-1].split("_")[-1], n))
    print("   D < 0,05        : %6d  (%4.1f %%)  insertions quasi intactes"
          % (q[0], 100 * q[0] / n))
    print("   0,05 <= D < 0,5 : %6d  (%4.1f %%)  amorcees" % (q[1], 100 * q[1] / n))
    print("   0,5  <= D < 1   : %6d  (%4.1f %%)  en cours" % (q[2], 100 * q[2] / n))
    print("   D >= 0,999      : %6d  (%4.1f %%)  ROMPUES" % (q[3], 100 * q[3] / n))
    print("   D moyen %.3f | mediane %.3f" % (D.mean(), np.median(D)))
