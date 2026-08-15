# -*- coding: utf-8 -*-
"""Phase A — SELECTION D'ARCHITECTURE par la mesure.

Trois architectures sur un protocole identique (memes maillages, memes
schemas, meme jeu de joints de depart) :
  A1  bulk elastique + crushCap        (ecole Bu 2026 / Jiang 2025)
  A2  bulk law = mc                    (ecole Ye 2025, loi implementee 15/08)
  A3  bulk law = dpdfh                 (carte Red Bohus de la these)

Essais : UCS, bresilien (3 GRAINES — le BTS numerique a un COV de 8-20 %,
Jiang et al. fig. 5), triaxial s3 = 20 et 50 MPa.

Criteres de selection (cf. METHODOLOGIE §3) : (a) les 4 essais cassent-ils,
(b) l'enveloppe est-elle CONCAVE, (c) le ratio UCS/BTS est-il credible
(cible experimentale 12,3), (d) cout par run.

usage : python run_phaseA.py [--dry]
"""
import csv, os, re, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
EXE = os.path.join(ROOT, "rockim_chk.exe")
CFGD = os.path.normpath(os.path.join(HERE, "..", "configs"))
RUND = os.path.normpath(os.path.join(HERE, "..", "runs"))
DRY = "--dry" in sys.argv

# --- bloc commun : materiau de depart + schemas figes de la campagne -------
COMMON = """rho = 2620
E = 77.66e9
nu = 0.29
ft = 34e6
cohesion = 13.6e6
frictionDeg = 13.4
Gf = 70
gfShearFactor = 10

insertion = adaptive
jointSoftening = yan
jointShearUnload = origin
gcActivation = adaptive
contact = potential
jointPenaltyFactor = 100
insertionPenaltyFactor = 4

dampingLocal = 0.7
jointXi = 0.0
verifyFt = false
budgetAbortPct = 5
budgetAbortMin = 0.05
"""

# --- ce qui distingue les trois architectures ------------------------------
ARCH = {
    "A1": "crushCap = 400e6\n",
    "A2": "law = mc\nmcCohesion = 13.6e6\nmcFrictionDeg = 40\nmcDilationDeg = 10\n",
    "A3": "law = dpdfh\n",
}

# --- maillage GBM : grains 5 mm (economique pour la phase A ; la production
# passera a 2,5 mm = taille reelle Bohus), Delaunay intra-grain 1,5 mm
MESH = """mesh = voronoi
grainSize = 0.005
grainSeeding = random
grainJitter = 0.5
lloydIters = 2
grainMesh = delaunay
grainElemSize = 0.0015
gbAlphaE = 1.0
gbAlphaTen = 1.0
gbAlphaCoh = 1.0
gbAlphaGf = 1.0
gbAlphaFric = 1.0
"""

UCS = """mode = fdem
scenario = tension
loading = platens
T = 4e-3
frames = 6
W = 0.05
H = 0.10
thickness = 1.0
pullV = -0.2
pullRamp = 2e-4
contactMu = 0.1
ucsStopAfterPeak = true
"""

TRIAX = UCS + """confiningPressure = {s3}e6
confiningRamp = 2e-4
confineFaces = sides
confineGaugeTime = 6e-4
pullDelay = 6e-4
"""

BTS = """mode = fdem
scenario = brazilian
geometry = disc
T = 1.1e-3
frames = 8
W = 0.0494
H = 0.0494
thickness = 1.0
brazilianLoading = platens
pullV = 0.1
pullRamp = 2e-4
platenHalfWidth = 0.0025
contactMu = 0.1
discMesh = native
discFlattenDeg = 20
elasticGaugeLo = 0.3
elasticGaugeHi = 0.8
brazilianStopAfterPeak = true
brazilianStopDelay = 6e-5
"""


def cases():
    for a in ("A1", "A2", "A3"):
        yield f"{a}_ucs", a, UCS + MESH + "seed = 4211\n"
        for sd in (4211, 4212, 4213):
            yield f"{a}_bts_s{sd}", a, BTS + MESH + f"seed = {sd}\n"
        for s3 in (20, 50):
            yield f"{a}_tx{s3}", a, TRIAX.format(s3=s3) + MESH + "seed = 4211\n"


PATS = {
    "peak_MPa": r"peak (?:axial|macro) stress\s*=\s*([-\d.eE+]+)",
    "sigma_t_MPa": r"indirect tensile strength sigma_t = 2P/\(pi D t\) = ([-\d.eE+]+)",
    "broken": r"broken joints\s*=\s*(\d+)",
    "shear_pct": r"([\d.]+) % shear",
    "s3_ach_MPa": r"achieved mean sigma_xx.*?=\s*([-\d.eE+]+)",
    "wall_s": r"wall time: ([\d.]+) s",
    "resid_pct": r"residu\s*:\s*[-\d.eE+]+ J/m \(([\d.eE+-]+) %",
    "diametral_pct": r"\((\d+) % diametral\)",
}


def main():
    os.makedirs(CFGD, exist_ok=True)
    os.makedirs(RUND, exist_ok=True)
    rows = []
    for name, arch, body in cases():
        cfg = os.path.join(CFGD, f"phaseA_{name}.cfg")
        with open(cfg, "w") as f:
            f.write(f"# Phase A — {name} ({arch})\n" + body + COMMON + ARCH[arch])
        out = os.path.join(RUND, f"phaseA_{name}")
        if DRY:
            print("[dry]", name)
            continue
        t0 = time.time()
        p = subprocess.run([EXE, cfg, out], cwd=ROOT, capture_output=True,
                           text=True, errors="replace")
        log = p.stdout + p.stderr
        with open(os.path.join(RUND, f"phaseA_{name}.log"), "w") as f:
            f.write(log)
        r = {"case": name, "arch": arch, "rc": p.returncode,
             "elapsed_s": round(time.time() - t0, 1)}
        for k, pat in PATS.items():
            m = re.search(pat, log)
            r[k] = m.group(1) if m else ""
        rows.append(r)
        print("%-14s rc=%d  pic=%-8s btsg=%-8s casse=%-5s  %5.0f s"
              % (name, p.returncode, r["peak_MPa"], r["sigma_t_MPa"],
                 r["broken"], r["elapsed_s"]), flush=True)

    if rows:
        csvp = os.path.normpath(os.path.join(HERE, "..", "phaseA_results.csv"))
        with open(csvp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print("ecrit", csvp)


if __name__ == "__main__":
    main()
