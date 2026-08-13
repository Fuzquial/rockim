#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# verify_suite.py — runner de non-régression rockim (patron Yade : références
# en dur, tolérances PAR NATURE de test, un point d'entrée, sortie table +
# code retour). Écrit le 2026-08-11 pour verrouiller la campagne Yan avant
# toute modification du code.
#
#   python3 tools/verify_suite.py --exe build/rockim               # tier fast
#   python3 tools/verify_suite.py --exe build/rockim --tier full
#   python3 tools/verify_suite.py --exe build/rockim --only fdem3d
#   python3 tools/verify_suite.py ... --update-refs refs_linux.json
#
# Règles maison encodées ici :
#   * OMP_NUM_THREADS = 1 par défaut (certification à 1 thread) ;
#   * les contrôles à CHARGE NULLE (pullV = 1e-12) exigent 0 joint cassé —
#     le test le moins cher et le plus discriminant du dépôt ;
#   * dampWork <= 0 partout où le solveur l'imprime (un dashpot ne peut que
#     dissiper) ;
#   * les références Voronoï sont PROPRES À LA PLATEFORME (les distributions
#     std divergent entre libstdc++ et MSVC à graine égale) : le jeu embarqué
#     ici est la baseline Linux/libstdc++ du 2026-08-11.
# ---------------------------------------------------------------------------
import argparse, json, os, re, subprocess, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CFG = os.path.join(ROOT, "configs")

# ---- extracteurs (regex sur stdout) ---------------------------------------
RX = {
    "err_pct":   r"error = (-?[\d.eE+]+) %",
    "peak_mpa":  r"peak macro stress = (-?[\d.eE+]+) MPa",
    "ucs_mpa":   r"peak axial stress = (-?[\d.eE+]+) MPa",
    "shear_pct": r"([\d.eE+-]+) % shear\)",
    "broken":    r"broken joints\s*[:=]\s*(\d+)",
    "dampwork":  r"joint dashpot work: (-?[\d.eE+-]+) J",
    "yan_int":   r"int f\(D\) dD = ([\d.eE+-]+)",
    "inserted":  r"adaptive insertion: (\d+) / \d+ joints inserted",
    "ratio":     r"measured/expected ratio = ([\d.eE+-]+)",
}

# ---- définition des tests --------------------------------------------------
# check = (extracteur, référence, tolérance ABSOLUE, obligatoire?)
# ref None => seулement présence/PASS ; "PASS" => chercher [PASS] dans stdout.
# Baseline Linux g++13/libstdc++/Eigen 3.4, OMP=1, source rockim_partage_2026-08-11.
TESTS = [
    # --- tier fast : lois, 2D, charge nulle 2D (≈ 1 min) --------------------
    dict(name="selftest_saksala2011", tier="fast", selftest="selftest-saksala2011",
         checks=[]),                                    # exit 0 == PASS (throw sinon)
    dict(name="selftest_dpdfh", tier="fast", selftest="selftest-dpdfh", checks=[]),
    dict(name="yan_integral", tier="fast", cfg="verify_fdem_tension.cfg",
         over=["jointSoftening = yan", "T = 2e-6"],
         # stdout n'imprime que 6 décimales : la vérification à 1e-12 vit dans
         # tools/yan_point.cpp ; ici on verrouille la valeur imprimée
         checks=[("yan_int", 0.386307294744, 5e-7, True)]),
    dict(name="fem_bar", tier="fast", cfg="verify_fem_bar.cfg",
         checks=[("err_pct", 0.889143, 0.02, True), ("pass_tag", None, 0, True)]),
    dict(name="dem_tension", tier="fast", cfg="verify_dem_tension.cfg",
         checks=[("err_pct", -0.109232, 0.02, True), ("pass_tag", None, 0, True)]),
    dict(name="fdem_voronoi_tension", tier="fast", cfg="verify_fdem_voronoi_tension.cfg",
         checks=[("err_pct", 11.5981, 0.05, True), ("pass_tag", None, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),      # <= 0 (tol vers le +)
    dict(name="zeroload_2d_fan", tier="fast", cfg="verify_fdem_voronoi_tension.cfg",
         over=["pullV = 1e-12", "verifyFt = false"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # même contrôle sous la loi de décharge en cisaillement de l'eq. 18 : la
    # sécante à l'origine ne doit rien créer au repos non plus (A2, 2026-08-13)
    dict(name="zeroload_2d_origin", tier="fast", cfg="verify_fdem_voronoi_tension.cfg",
         over=["pullV = 1e-12", "verifyFt = false", "jointShearUnload = origin"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # le cas delaunay est LE cas historique du cliquet (158 joints cassés a
    # charge nulle avant correctif) mais coute ~2.5 min : tier full
    dict(name="zeroload_2d_delaunay", tier="full", cfg="verify_fdem_voronoi_tension.cfg",
         over=["pullV = 1e-12", "grainMesh = delaunay", "verifyFt = false"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # --- tier full : bit-repères 2D longs, adaptatif, 3D grille -------------
    dict(name="fdem_tension", tier="full", cfg="verify_fdem_tension.cfg",
         checks=[("err_pct", -1.38789, 0.01, True), ("pass_tag", None, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="fdem_tension_adaptive", tier="full", cfg="verify_fdem_tension.cfg",
         over=["insertion = adaptive"],
         checks=[("err_pct", -4.08193, 0.02, True), ("pass_tag", None, 0, True)]),
    dict(name="zeroload_2d_grid", tier="full", cfg="verify_fdem_tension.cfg",
         over=["pullV = 1e-12"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    dict(name="fdem3d_tension", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["jointXi = 0"],                          # doctrine 2026-08-05 :
         # une vérification de la LOI ne mesure pas la dissipation visqueuse
         checks=[("err_pct", -3.04802, 0.02, True), ("pass_tag", None, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="fdem3d_tension_adaptive", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["jointXi = 0", "insertion = adaptive"],
         checks=[("err_pct", -0.034309, 0.02, True), ("pass_tag", None, 0, True)]),
    dict(name="dem3d_tension", tier="full", cfg="verify_dem3d_tension.cfg",
         checks=[("err_pct", -0.0588759, 0.02, True), ("pass_tag", None, 0, True)]),
    dict(name="zeroload_3d_grid", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["pullV = 1e-12"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    dict(name="fdem3d_tension_yan", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["jointXi = 0", "jointSoftening = yan"],
         checks=[("err_pct", -1.10478, 0.05, True), ("pass_tag", None, 0, True),
                 ("yan_int", 0.386307294744, 1e-6, True)]),
    # --- A2 : décharge en cisaillement sur sécante à l'origine (eq. 18) -----
    # forme LITTÉRALE de l'article = origin + jointFrictionScaled = 1.
    # En traction (mode I dominant) l'écart au retour radial doit rester
    # marginal : c'est ce que ces deux repères verrouillent.
    dict(name="fdem_tension_origin", tier="full", cfg="verify_fdem_tension.cfg",
         over=["jointSoftening = yan", "jointFrictionScaled = 1",
               "jointShearUnload = origin"],
         checks=[("err_pct", -1.721, 0.02, True), ("pass_tag", None, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="fdem3d_tension_origin", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["jointXi = 0", "jointSoftening = yan", "jointFrictionScaled = 1",
               "jointShearUnload = origin"],
         checks=[("err_pct", -1.08522, 0.05, True), ("pass_tag", None, 0, True)]),
    dict(name="zeroload_3d_origin", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["pullV = 1e-12", "jointShearUnload = origin"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # UCS de la fig. 17 de l'article : le SEUL repère du parc dont la rupture
    # soit pilotée par le CISAILLEMENT (~47 % des joints), et donc le seul qui
    # exerce vraiment la loi de mode II — insertion adaptative + f(D) comprises.
    dict(name="ucs_yan_adaptive", tier="full", cfg="../configs_yan/ucs_adap.cfg",
         checks=[("ucs_mpa", 51.0735, 0.15, True), ("broken", 327, 0, True),
                 ("inserted", 1288, 0, True)]),
    dict(name="ucs_yan_origin", tier="full", cfg="../configs_yan/ucs_adap.cfg",
         over=["jointShearUnload = origin", "jointFrictionScaled = 1"],
         checks=[("ucs_mpa", 50.3671, 0.15, True), ("broken", 318, 0, True),
                 ("inserted", 1131, 0, True)]),
    # charge nulle sur le maillage IMPORTE (mesh = file, tets Gmsh) : certifie
    # la machinerie d'import comme les autres controles certifient les leurs
    dict(name="zeroload_3d_filemesh", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["mesh = file",
               "meshFile = " + os.path.join(ROOT, "meshes", "box3d_h45.msh"),
               "pullV = 1e-12", "T = 1e-4"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # --- tier all : Voronoï 3D + charge nulle 3D Voronoï (~40 min) ----------
    dict(name="fdem3d_voronoi_tension", tier="all", cfg="verify_fdem3d_voronoi_tension.cfg",
         over=["jointXi = 0"],
         checks=[("err_pct", 8.45551, 0.10, True), ("pass_tag", None, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="zeroload_3d_voronoi", tier="all", cfg="verify_fdem3d_voronoi_tension.cfg",
         over=["pullV = 1e-12", "verifyFt = false"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, False)]),
]

TIERS = {"fast": ["fast"], "full": ["fast", "full"], "all": ["fast", "full", "all"]}


def run_one(exe, t, outroot, env, timeout):
    out = os.path.join(outroot, t["name"])
    if "selftest" in t:
        cmd = [exe, t["selftest"], os.path.join(out + ".csv")]
    else:
        cfg = os.path.join(CFG, t["cfg"])
        if t.get("over"):
            fd, cfg2 = tempfile.mkstemp(suffix=".cfg", prefix=t["name"] + "_")
            with os.fdopen(fd, "w") as f:
                f.write(open(cfg).read())
                f.write("\n# --- verify_suite overrides ---\n")
                for line in t["over"]:
                    f.write(line + "\n")
            cfg = cfg2
        cmd = [exe, cfg, out]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                           env=env)
    except subprocess.TimeoutExpired:
        return dict(name=t["name"], ok=False, dt=time.time() - t0,
                    detail=[f"TIMEOUT après {timeout}s"])
    dt = time.time() - t0
    text = p.stdout + p.stderr
    detail, ok = [], True
    if p.returncode != 0:
        return dict(name=t["name"], ok=False, dt=dt,
                    detail=[f"exit {p.returncode}: {text.strip().splitlines()[-1] if text.strip() else '?'}"])
    for (kind, ref, tol, required) in t["checks"]:
        if kind == "pass_tag":
            n_fail = text.count("[FAIL]")
            if n_fail or "[PASS]" not in text:
                ok = False
                detail.append(f"attendu [PASS], trouvé {n_fail} [FAIL]")
            continue
        m = re.findall(RX[kind], text)
        if not m:
            if required:
                ok = False
                detail.append(f"{kind}: non trouvé dans la sortie")
            else:
                detail.append(f"{kind}: absent (toléré)")
            continue
        val = float(m[-1])
        if kind == "dampwork":                          # <= 0 exigé
            if val > tol:
                ok = False
                detail.append(f"dampWork = {val:g} > 0 : INJECTION d'énergie")
            else:
                detail.append(f"dampWork = {val:g} <= 0 ok")
            continue
        if ref is None:
            detail.append(f"{kind} = {val:g} (référence à fixer)")
            continue
        if abs(val - ref) > tol:
            ok = False
            detail.append(f"{kind} = {val:g}, attendu {ref:g} ± {tol:g}")
        else:
            detail.append(f"{kind} = {val:g} ok (ref {ref:g})")
    return dict(name=t["name"], ok=ok, dt=dt, detail=detail)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default=os.path.join(ROOT, "build", "rockim"))
    ap.add_argument("--tier", default="fast", choices=list(TIERS))
    ap.add_argument("--only", default=None, help="filtre sur le nom")
    ap.add_argument("--threads", default="1")
    ap.add_argument("--timeout", type=float, default=3600.0)
    ap.add_argument("--json", default=None, help="rapport JSON")
    args = ap.parse_args()

    env = dict(os.environ, OMP_NUM_THREADS=args.threads)
    outroot = tempfile.mkdtemp(prefix="rockim_suite_")
    sel = [t for t in TESTS if t["tier"] in TIERS[args.tier]
           and (args.only is None or args.only in t["name"])]
    print(f"[suite] {len(sel)} tests, exe = {args.exe}, "
          f"OMP_NUM_THREADS = {args.threads}, sorties {outroot}")
    results, allok = [], True
    for t in sel:
        r = run_one(args.exe, t, outroot, env, args.timeout)
        results.append(r)
        allok &= r["ok"]
        print(f"  {'PASS' if r['ok'] else 'FAIL'}  {r['name']:<28s} "
              f"{r['dt']:7.1f}s  " + "; ".join(r["detail"]))
    print(f"[suite] {'TOUT PASSE' if allok else 'ÉCHECS'} "
          f"({sum(r['ok'] for r in results)}/{len(results)})")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=1)
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
