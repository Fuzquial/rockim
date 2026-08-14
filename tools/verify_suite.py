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
    "gcact":     r"adaptive contact activation: (\d+) / \d+",
    "gcwork":    r"net work (?:injected )?by general contact: (-?[\d.eE+-]+)",
    "pot_ke":    r"pot_ke_rel = ([\d.eE+-]+)",
    "pot_mom":   r"pot_mom_rel = ([\d.eE+-]+)",
    "pot3_ke":   r"pot3_ke_rel = ([\d.eE+-]+)",
    "pot3_mom":  r"pot3_mom_rel = ([\d.eE+-]+)",
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
    # A3 : LE test decisif du contact par potentiel — collision elastique sans
    # frottement, conservation jugee sur dKE (3.7e-12 frontale, 5.1e-7 oblique)
    # et quantite de mouvement machine (3e loi exacte par construction)
    dict(name="selftest_potential2d", tier="fast", selftest="selftest-potential2d",
         checks=[("pot_ke", 0.0, 1e-5, True), ("pot_mom", 0.0, 1e-12, True),
                 ("pass_tag", None, 0, True)]),
    # ... et son miroir 3D (tet-tet, A3 phase 2) : pointe-contre-face puis
    # oblique, dKE 2e-8, quantite de mouvement machine
    dict(name="selftest_potential3d", tier="fast", selftest="selftest-potential3d",
         checks=[("pot3_ke", 0.0, 1e-5, True), ("pot3_mom", 0.0, 1e-12, True),
                 ("pass_tag", None, 0, True)]),
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
    # ... ni l'activation adaptative du contact (A1, 2026-08-13)
    dict(name="zeroload_2d_gcadaptive", tier="fast", cfg="verify_fdem_voronoi_tension.cfg",
         over=["pullV = 1e-12", "verifyFt = false", "gcActivation = adaptive"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    # ... ni le contact par potentiel (A3, 2026-08-13)
    dict(name="zeroload_2d_potential", tier="fast", cfg="verify_fdem_voronoi_tension.cfg",
         over=["pullV = 1e-12", "verifyFt = false", "contact = potential"],
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
    # --- A1 : activation adaptative du contact (gcActivation) ---------------
    # La paire percussion 2D est BIT-IDENTIQUE mode a mode (174 joints, memes
    # energies, fichiers identiques — verifie le 2026-08-13, phase debris
    # comprise, apres le correctif du cache des faces mortes). La paire SHPB
    # (multi-corps, contact necessaire des t = 0) est bit-identique a
    # T = 9e-5 ; au-dela la tempete de fragments entre dans l'enveloppe
    # chaotique (le mode full lui-meme diverge 8x plus tot sous OMP = 2).
    dict(name="percussion_2d", tier="full", cfg="fdem_percussion.cfg",
         checks=[("broken", 174, 0, True), ("dampwork", 0.0, 1e-12, True)]),
    dict(name="percussion_2d_gcadaptive", tier="full", cfg="fdem_percussion.cfg",
         over=["gcActivation = adaptive"],
         checks=[("broken", 174, 0, True), ("gcact", 21, 0, True),
                 ("dampwork", 0.0, 1e-12, True)]),
    dict(name="shpb_mini", tier="full", cfg="../configs_yan/shpb_mini.cfg",
         over=["T = 9e-5", "frames = 4"],
         checks=[("broken", 812, 0, True)]),
    dict(name="shpb_mini_gcadaptive", tier="full", cfg="../configs_yan/shpb_mini.cfg",
         over=["T = 9e-5", "frames = 4", "gcActivation = adaptive"],
         checks=[("broken", 812, 0, True), ("gcact", 122, 0, True)]),
    # --- A3 : contact par potentiel de Munjiza (contact = potential) --------
    # Conservation au niveau SOLVEUR : assemblage SHPB incassable sans
    # frottement — |gcWork| doit rester au niveau du biais du compteur
    # (mesure -2.6 J/m pour 766 J/m en jeu, 0.3 %), la ou le penalty
    # quasi-plastique par defaut dissiperait massivement.
    dict(name="shpb_elastic_potential", tier="full", cfg="../configs_yan/shpb_mini.cfg",
         over=["T = 9e-5", "frames = 4", "phase.rock.ft = 1e12",
               "phase.rock.cohesion = 1e12", "phase.rock.Gf = 1e6",
               "contactMu = 0", "contact = potential"],
         checks=[("broken", 0, 0, True), ("gcwork", 0.0, 8.0, True)]),
    # SHPB reel en potentiel : reference de plateforme (chaos de broyage
    # verrouille, comme les autres). L'onde incidente colle au penalty a
    # 3e-6 pres ; le disque casse moins (595 vs 812) — loi de contact
    # differente a l'interface, ecart PHYSIQUE assume et documente.
    # ref 595 -> 578 le 13/08 au soir : l'ordre CANONIQUE des paires (N1) a
    # change une fois pour toutes l'ordre des sommes de forces en tempete de
    # fragments (enveloppe chaotique, cf. controle OMP=2). C'est le DERNIER
    # changement d'ordre possible : l'ordre est desormais independant de
    # l'implementation de la detection.
    dict(name="shpb_mini_potential", tier="full", cfg="../configs_yan/shpb_mini.cfg",
         over=["T = 9e-5", "frames = 4", "contact = potential"],
         checks=[("broken", 578, 0, True)]),
    # Percussion 2D en potentiel : verrouille (a) la borne du residu
    # d'injection de la releve de naissance (+27 J/m mesure, garde a 60 —
    # la regression a +936 sans releve doit rester impossible) et (b) le
    # comptage debris de plateforme.
    dict(name="percussion_2d_potential", tier="full", cfg="fdem_percussion.cfg",
         over=["contact = potential"],
         checks=[("broken", 5, 0, True), ("gcwork", 0.0, 60.0, True)]),
    # A3 phase 2 (3D) : charge nulle sous le potentiel tet-tet — c'est CE
    # controle qui a attrape les slivers de tets tangents (5 joints casses au
    # repos avant les gardes plancher-de-volume + fermeture du polyedre)
    # (T reduit a 1e-4 : les slivers existent des le pas 0, la fenetre courte
    # suffit a les attraper et le controle coute ~5 min au lieu de ~11)
    dict(name="zeroload_3d_potential", tier="full", cfg="verify_fdem3d_tension.cfg",
         over=["pullV = 1e-12", "T = 1e-4", "contact = potential"],
         checks=[("broken", 0, 0, True), ("dampwork", 0.0, 1e-12, True)]),
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
    # A1 en 3D sur le maillage Gmsh importe : paire bit-identique a T = 5e-5
    # (6 joints, memes energies), et le x2.4 de mur est deja la (265 -> 111 s)
    dict(name="percussion_3d_gmsh", tier="all", cfg="fdem3d_percussion_base.cfg",
         over=["meshFile = " + os.path.join(ROOT, "meshes", "box3d_h45.msh"),
               "T = 5e-5", "frames = 4"],
         checks=[("broken", 6, 0, True)]),
    dict(name="percussion_3d_gmsh_gcadaptive", tier="all", cfg="fdem3d_percussion_base.cfg",
         over=["meshFile = " + os.path.join(ROOT, "meshes", "box3d_h45.msh"),
               "T = 5e-5", "frames = 4", "gcActivation = adaptive"],
         checks=[("broken", 6, 0, True), ("gcact", 0, 0, True)]),
    # A3 phase 2 : percussion 3D Gmsh sous le potentiel tet-tet (memes 6
    # joints que la penalite a ce T — le broyage debute a peine et les
    # debris n'ont pas encore d'influence), gcWork garde pres de zero
    dict(name="percussion_3d_gmsh_potential", tier="all", cfg="fdem3d_percussion_base.cfg",
         over=["meshFile = " + os.path.join(ROOT, "meshes", "box3d_h45.msh"),
               "T = 5e-5", "frames = 4", "contact = potential"],
         checks=[("broken", 6, 0, True), ("gcwork", 0.0, 1.0, True)]),
    # V1 : DEUX CORPS (groupes physiques Gmsh — insert sphere + bloc, aucun
    # joint inter-corps) au repos : l'insert immobile a 0.5 mm au-dessus ne
    # doit RIEN ressentir — 0 casse, travail de contact exactement nul.
    dict(name="zeroload_bench1_3d", tier="full", cfg="fdem3d_bench1_insert.cfg",
         over=["meshFile = " + os.path.join(ROOT, "meshes", "bench1_insert.msh"),
               "T = 2e-5", "frames = 2", "groupVel.insert = 0 0 0"],
         checks=[("broken", 0, 0, True), ("gcwork", 0.0, 0.0, True)]),
    # V1 : l'impact complet insert -> roche (2.53 J a -8 m/s). Reference du
    # 2026-08-14 : 12 joints (2 traction / 10 cisaillement), rebond a
    # +5.67 m/s (e = 0.71 — la restitution du potentiel, cf. percussion 2D),
    # gcWork absorbant. ~45 min : le prix du jalon multi-corps.
    dict(name="bench1_insert_impact", tier="all", cfg="fdem3d_bench1_insert.cfg",
         over=["meshFile = " + os.path.join(ROOT, "meshes", "bench1_insert.msh")],
         checks=[("broken", 12, 0, True), ("gcwork", 0.0, 1.0, True)]),
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
