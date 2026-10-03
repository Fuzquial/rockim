#!/usr/bin/env python3
"""Rang 3 : balayage jointContactPenalty x dtFactor x v0 sur le micro-banc a
deux tetraedres (tests_f2/rang3/base.cfg). Indicateur : rapport de l'enveloppe
de KE du corps B (0,5 m vz^2, vz = vitesse moyenne du corps) entre le dernier
quart du run et le deuxieme quart (apres l'endommagement initial). Un rapport
> 1 = l'amplitude d'oscillation CROIT : energie creee. Aussi : D moyen final,
joint rompu ou non, residu B4.
  python3 tests_f2/rang3/balayage.py build/rockim <dossier_sortie>"""
import csv, os, re, subprocess, sys
exe, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
base = open("tests_f2/rang3/base.cfg").read()
M = 3.09477e-07
rows = []
for v0 in (12, 13):
    for jcp in ("fixed", "adaptive"):
        for dtf in (0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 0.9):
            name = f"v{v0}_{jcp}_dt{dtf}"
            cfg = os.path.join(out, name + ".cfg")
            with open(cfg, "w") as f:
                f.write(base + f"\ngroupVel.B = 0 0 {v0}\nT = 2e-4\n"
                        f"jointContactPenalty = {jcp}\ndtFactor = {dtf}\n")
            p = subprocess.run([exe, cfg, os.path.join(out, name)],
                               capture_output=True, text=True,
                               env=dict(os.environ, OMP_NUM_THREADS="1"))
            log = p.stdout + p.stderr
            open(os.path.join(out, name + ".log"), "w").write(log)
            if p.returncode != 0:
                rows.append((v0, jcp, dtf, "ERREUR", "", "", "", ""))
                continue
            h = list(csv.DictReader(open(os.path.join(out, name, "history.csv"))))
            ke = [0.5 * M * float(r["vz_B"]) ** 2 for r in h]
            n = len(ke)
            q2 = max(ke[n // 4: n // 2]) if n > 8 else float("nan")
            q4 = max(ke[3 * n // 4:]) if n > 8 else float("nan")
            D = re.findall(r"D moyen ([0-9.eE+-]+)", log)
            br = re.findall(r"broken joints\s*:\s*(\d+)", log)
            res = re.findall(r"residu\s+: (-?[0-9.eE+-]+) J \(([0-9.eE+-]+) %", log)
            ab = "ABORT" if "ENERGY ABORT" in log else ""
            rows.append((v0, jcp, dtf, f"{q4 / q2:.3f}", D[-1] if D else "",
                         br[-1] if br else "", res[-1][1] if res else "", ab))
with open(os.path.join(out, "balayage.csv"), "w") as f:
    f.write("v0,jointContactPenalty,dtFactor,KE_q4_sur_q2,D_moyen_fin,rompus,residu_pct,abort\n")
    for r in rows:
        f.write(",".join(map(str, r)) + "\n")
for r in rows:
    print(*r, sep="\t")
