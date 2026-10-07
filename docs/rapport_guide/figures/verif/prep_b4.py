"""Extrait les blocs « energy budget (V2/B4) » des journaux des calculs courts
rejoues le 2026-10-06 (dossier b4/ des calculs) et ecrit b4_residus.csv et
b4_postes.csv. Usage : python prep_b4.py <dossier des calculs rejoues>
"""
import os, re, sys, csv
RUNS = sys.argv[1] if len(sys.argv) > 1 else "."
CASES = [
    ("bar", "barre en traction, élastique (fdem3d)", "J"),
    ("rupt", "barre, rupture adaptative (fdem3d)", "J"),
    ("zl_bench1", "deux corps au repos (fdem3d)", "J"),
    ("perc3d_20us", "percussion 3D, 20 µs (fdem3d)", "J"),
    ("perc_goff", "pesanteur non comptée (fdem3d, 2 µs)", "J"),
    ("perc_gon", "pesanteur comptée (fdem3d, 2 µs)", "J"),
    ("t1_pen", "raclage T1, pénalité (fdem)", "J/m"),
    ("t1_sig", "raclage T1, Signorini (fdem)", "J/m"),
    ("perc2d", "percussion 2D, 250 µs (fdem)", "J/m"),
]
POSTES = ["elements", "joints", "contact", "Cundall", "frontieres", "outil->solide",
          "charges", "integration"]
rows, prow = [], []
for key, lab, unit in CASES:
    s = open(os.path.join(RUNS, "b4", key + ".log"), encoding="utf-8", errors="ignore").read()
    m = re.search(r"residu\s+: (\S+) J(?:/m)? \((\S+) % de l'echelle\) \[([^\]]+)\]", s)
    ke = re.search(r"energy budget \(V2/B4\): KE (\S+) -> (\S+) J", s)
    inj = re.search(r"injection outil\s+: \S+ J/m vers le solide / \S+ J/m corps rigide = ratio (\S+)", s)
    rows.append([key, lab, unit, m.group(1), m.group(2), m.group(3), ke.group(1), ke.group(2),
                 inj.group(1) if inj else ""])
    vals = {}
    for p in POSTES:
        mm = re.search(r"\]\s+%s\s*: ([-+]?\S+) J" % re.escape(p), s)
        vals[p] = mm.group(1) if mm else ""
    prow.append([key] + [vals[p] for p in POSTES])
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b4_residus.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["cle", "cas", "unite", "residu", "residu_pct_echelle", "verdict", "Ec0", "Ec_fin", "ratio_injection_outil"])
    w.writerows(rows)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "b4_postes.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["cle"] + POSTES)
    w.writerows(prow)
