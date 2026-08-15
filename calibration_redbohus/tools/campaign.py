# -*- coding: utf-8 -*-
"""Moteur de campagne de calibration Red Bohus â€” architecture A1 (retenue
en phase A : la seule qui casse dans TOUS les regimes de confinement).

Genere les configs, lance les runs EN PARALLELE, extrait les observables.

  python campaign.py screen        # phase B : criblage mono-variable
  python campaign.py lhs   [N]     # phase C : hypercube latin
  python campaign.py points f.json # rejoue une liste de jeux (validation)

Sorties : <phase>_results.csv (un jeu de parametres par ligne, toutes les
observables) + les logs dans runs/.
"""
import concurrent.futures as cf
import csv, json, os, re, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
EXE = os.path.join(ROOT, "rockim_chk.exe")
CFGD = os.path.normpath(os.path.join(HERE, "..", "configs"))
RUND = os.path.normpath(os.path.join(HERE, "..", "runs"))
NWORK = 3                      # runs simultanes
NTHREAD = 6                    # threads OpenMP par run (18 coeurs logiques)

# --- espace des parametres (architecture A1) --------------------------------
# bornes elargies facon Bu 2026 / Jiang 2025, centrees sur ce que la phase A
# a montre : ft deja bon (BTS -3,8 %), phi tres insuffisant (enveloppe plate)
SPACE = {
    "ft":            (20e6, 60e6),      # resistance en traction des joints
    "cohesion":      (10e6, 90e6),      # cohesion des joints
    "frictionDeg":   (15.0, 60.0),      # frottement des joints â€” LE levier
    "Gf":            (30.0, 300.0),     # energie de rupture mode I
    "gfShearFactor": (1.0, 12.0),       # rapport GII/GI
    "crushCap":      (200e6, 1500e6),   # cap deviatorique du bulk
}
CENTER = {"ft": 34e6, "cohesion": 40e6, "frictionDeg": 40.0,
          "Gf": 70.0, "gfShearFactor": 10.0, "crushCap": 400e6}

MESH = """mesh = voronoi
grainSize = 0.005
grainSeeding = random
grainJitter = 0.5
lloydIters = 2
grainMesh = delaunay
grainElemSize = 0.002
gbAlphaE = 1.0
gbAlphaTen = 1.0
gbAlphaCoh = 1.0
gbAlphaGf = 1.0
gbAlphaFric = 1.0
"""
SCHEMES = """rho = 2620
E = 77.66e9
nu = 0.29
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
UCS = """mode = fdem
scenario = tension
loading = platens
T = 4e-3
frames = 2
W = 0.04
H = 0.08
thickness = 1.0
pullV = -0.2
pullRamp = 2e-4
contactMu = 0.1
ucsStopAfterPeak = true
"""
BTS = """mode = fdem
scenario = brazilian
geometry = disc
T = 1.1e-3
frames = 2
W = 0.04
H = 0.04
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


def cfg_text(p, test, seed):
    body = {"ucs": UCS, "bts": BTS}.get(test)
    if body is None:                                   # tx20, tx50, tx75...
        s3 = int(test[2:])
        body = UCS + ("confiningPressure = %de6\nconfiningRamp = 2e-4\n"
                      "confineFaces = sides\nconfineGaugeTime = 6e-4\n"
                      "pullDelay = 6e-4\n" % s3)
    mat = "".join("%s = %g\n" % (k, v) for k, v in p.items())
    return body + MESH + SCHEMES + mat + "seed = %d\n" % seed


PATS = {
    "peak_MPa": r"peak (?:axial|macro) stress\s*=\s*([-\d.eE+]+)",
    "sigma_t_MPa": r"indirect tensile strength sigma_t = 2P/\(pi D t\) = ([-\d.eE+]+)",
    "broken": r"broken joints\s*=\s*(\d+)",
    "shear_pct": r"([\d.]+) % shear",
    "diametral_pct": r"\((\d+) % diametral\)",
    "wall_s": r"wall time: ([\d.]+) s",
    "locked": r"peak (LOCKED|NOT locked)",
}


def one_run(job):
    tag, p, test, seed = job
    name = "%s_%s_s%d" % (tag, test, seed)
    cfg = os.path.join(CFGD, "c_%s.cfg" % name)
    with open(cfg, "w") as f:
        f.write(cfg_text(p, test, seed))
    env = dict(os.environ, OMP_NUM_THREADS=str(NTHREAD))
    t0 = time.time()
    try:
        r = subprocess.run([EXE, cfg, os.path.join(RUND, name)], cwd=ROOT,
                           capture_output=True, text=True, errors="replace",
                           env=env, timeout=1800)
        log = r.stdout + r.stderr
        rc = r.returncode
    except subprocess.TimeoutExpired:
        log, rc = "TIMEOUT", -9
    with open(os.path.join(RUND, name + ".log"), "w") as f:
        f.write(log)
    out = {"test": test, "seed": seed, "rc": rc,
           "elapsed_s": round(time.time() - t0, 1)}
    for k, pat in PATS.items():
        m = re.search(pat, log)
        out[k] = m.group(1) if m else ""
    return tag, out


def run_all(sets, tests, seeds=(4211,)):
    """sets : {tag: params}. Retourne {tag: {test: obs}}."""
    jobs = [(tag, p, t, s) for tag, p in sets.items()
            for t in tests for s in seeds]
    res = {tag: {} for tag in sets}
    done = 0
    with cf.ThreadPoolExecutor(max_workers=NWORK) as ex:
        for tag, out in ex.map(one_run, jobs):
            key = out["test"] if len(seeds) == 1 else "%s_s%d" % (out["test"], out["seed"])
            res[tag][key] = out
            done += 1
            print("  [%3d/%3d] %-22s %-6s pic=%-8s bts=%-8s %5.0f s"
                  % (done, len(jobs), tag, out["test"], out["peak_MPa"],
                     out["sigma_t_MPa"], out["elapsed_s"]), flush=True)
    return res


def flatten(sets, res, path):
    rows = []
    for tag, p in sets.items():
        row = {"tag": tag}
        row.update(p)
        for test, o in sorted(res[tag].items()):
            for k in ("peak_MPa", "sigma_t_MPa", "broken", "shear_pct",
                      "diametral_pct", "locked", "rc", "elapsed_s"):
                row["%s_%s" % (test, k)] = o.get(k, "")
        rows.append(row)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("ecrit", path)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "screen"
    os.makedirs(CFGD, exist_ok=True)
    os.makedirs(RUND, exist_ok=True)
    tests = ("ucs", "bts", "tx20")            # tx50 garde pour la validation

    if mode == "screen":
        sets = {"C": dict(CENTER)}
        for k, (lo, hi) in SPACE.items():
            for tag, v in (("lo", lo), ("hi", hi)):
                p = dict(CENTER); p[k] = v
                sets["%s_%s" % (k, tag)] = p
        out = os.path.normpath(os.path.join(HERE, "..", "screen_results.csv"))

    elif mode == "lhs":
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 40
        rng = np.random.default_rng(12345)
        keys = list(SPACE)
        # hypercube latin (permutations independantes, centre des strates)
        U = np.empty((n, len(keys)))
        for j in range(len(keys)):
            U[:, j] = (rng.permutation(n) + 0.5) / n
        sets = {}
        for i in range(n):
            p = {k: SPACE[k][0] + U[i, j] * (SPACE[k][1] - SPACE[k][0])
                 for j, k in enumerate(keys)}
            sets["L%03d" % i] = p
        out = os.path.normpath(os.path.join(HERE, "..", "lhs_results.csv"))

    elif mode == "points":
        sets = json.load(open(sys.argv[2]))
        tests = tuple(sys.argv[3].split(",")) if len(sys.argv) > 3 else tests
        out = os.path.normpath(os.path.join(HERE, "..", "points_results.csv"))
    else:
        raise SystemExit("mode inconnu")

    print("%d jeux x %d essais = %d runs, %d en parallele"
          % (len(sets), len(tests), len(sets) * len(tests), NWORK), flush=True)
    res = run_all(sets, tests)
    flatten(sets, res, out)


if __name__ == "__main__":
    main()

