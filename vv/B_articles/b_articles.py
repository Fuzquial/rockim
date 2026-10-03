#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# b_articles.py : partie III de la campagne V&V, reproductions d'articles
# deja faites (B4 a B7), consolidees et rejouables avec le binaire courant.
#
#   B4  AbuAisha et al. 2017, J. Petrol. Sci. Eng. 154 : fracturation
#       hydraulique en paroi de forage (fdem 2D, module hydro)
#   B5  Wang et al. 2024, Front. Earth Sci. 12:1517816 : EDZ d'un tunnel
#       profond, sigma0 = 5 MPa (fdem 2D, in situ + excavation)
#   B6  Yan, Zheng et Wang 2023, IJRMMS 169 : UCS de la fig. 17, insertion
#       adaptative (fdem 2D)
#   B7  Yang et al. 2025, IJRMMS 191 (St Anne, 10,66 m/s) et Yang et al.
#       2026, IJRMMS 206 (Kuru, 9 m/s) : impact d'insert unique (fdem3d)
#
# Les decks sources ont ete COPIES dans decks/ (fins de ligne et encodage
# normalises, contenu inchange). Le script les reecrit dans out/<cas>/deck.cfg
# en ne changeant que meshFile (maillages de meshes/, regeneres ici), T et
# frames quand c'est dit dans CASES, et outputDir. Aucune cle nouvelle.
#
#   python3 vv/B_articles/b_articles.py prepare      maillages + decks
#   python3 vv/B_articles/b_articles.py smoke        fumee : T court, 1 fil
#   python3 vv/B_articles/b_articles.py run [--cas B6_ucs ...] [--threads 2]
#   python3 vv/B_articles/b_articles.py analyse
#   python3 vv/run_queue.py B_articles --slots 2 --threads 2
# Les cas lourds (St Anne rock137, essais 3 AbuAisha) ne sont dans JOBS que si
# VV_B_LOURDS=1 est pose dans l'environnement.
# ---------------------------------------------------------------------------
import argparse, glob, json, os, re, shutil, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
VV = os.path.dirname(HERE)
ROOT = os.path.dirname(VV)
sys.path.insert(0, VV)
import vvcommon as C  # noqa: E402

DECKS = os.path.join(HERE, "decks")
MESH = os.path.join(HERE, "meshes")
OUT = os.path.join(HERE, "out")
OUT_SMOKE = os.path.join(OUT, "_fumee")
PHD = os.path.expanduser("~/phd_geothermie_git")
PY = sys.executable

# ---- maillages (*.msh ignores par git) : commandes des etudes d'origine -----
MESHES = {
    "hf_bore.msh": ["tunnel_edz/tools/make_circle_mesh.py", "8.0", "8.0", "0.05", "0.003",
                    "0.4", "0.3", "{out}", "1"],
    "parker_crack.msh": ["bench_abuaisha/tools/make_crack_mesh.py", "8.0", "8.0", "1.5",
                         "0.003", "0.3", "{out}", "1"],
    "parker_crack_c.msh": ["bench_abuaisha/tools/make_crack_mesh.py", "8.0", "8.0", "1.5",
                           "0.012", "0.3", "{out}", "1"],
    "tunnel_hs_iso.msh": ["tools/make_unstructured_mesh.py", "tunnelhs", "100", "100", "0.22",
                          "18", "2.0", "{out}", "1"],
    # docs/MAILLAGE_serie_2026-09-13.md : rock25 = SR 2,5 ; rock137 = SR 1,0
    # (et non 1,37 comme l'ecrit REPRODUIRE_stanne_radiales_2026-09-14.md)
    "impact_yang_train1_rock25_hxt.msh": ["tools/make_impact_mesh.py", "{out}", "1.0", "2e-5",
                                          "2.5", "gap=2e-5", "quality=hxt", "train=fixed"],
    "impact_yang_train1_rock137_hxt.msh": ["tools/make_impact_mesh.py", "{out}", "1.0", "2e-5",
                                           "1.0", "gap=2e-5", "quality=hxt", "train=fixed"],
    "impact_yang_s2.5_pose.msh": ["tools/make_impact_mesh.py", "{out}", "2.5", "gap=2e-5",
                                  "gapr=2e-5"],
}

# ---- cas ---------------------------------------------------------------------
# over : surcharges de T / frames seulement (la physique avant le pic ne depend
# pas de T en explicite). cost_h2 : cout estime a 2 fils (h), mesure par la
# fumee du 2026-10-04 (voir README) ; heavy : hors JOBS par defaut.
CASES = {
    # B4. Runs de reference du 2026-08-20 (hf_*_hydro, pompe des t = 0). Le pic
    # anisotrope est a 2,93 ms, l'isotrope a 3,18 ms : le deck isotrope d'origine
    # (T = 3,0 ms) s'arretait AVANT le pic ; l'historique archive va a 4 ms.
    "B4_aniso": dict(bench="B4", deck="b4_hf_aniso_hydro.cfg", mesh="hf_bore.msh",
                     over={"T": "3.3e-3", "frames": "11"}, kind="hydro",
                     prev_csv="FDEM/rockim/bench_abuaisha/historiques/hf_aniso.csv"),
    "B4_iso": dict(bench="B4", deck="b4_hf_iso_hydro.cfg", mesh="hf_bore.msh",
                   over={"T": "3.5e-3", "frames": "11"}, kind="hydro",
                   prev_csv="FDEM/rockim/bench_abuaisha/historiques/hf_iso.csv"),
    # protocole exact de l'article (hydroStart) : lourd, en complement
    "B4_e3_aniso": dict(bench="B4", deck="b4_e3_aniso.cfg", mesh="hf_bore.msh",
                        over={"T": "3.7e-3", "frames": "11"}, kind="hydro", heavy=True,
                        prev_csv="FDEM/rockim/bench_abuaisha/historiques/e3_aniso.csv"),
    "B4_e3_iso12": dict(bench="B4", deck="b4_e3_iso12.cfg", mesh="hf_bore.msh",
                        over={"T": "3.5e-3", "frames": "11"}, kind="hydro", heavy=True,
                        prev_csv="FDEM/rockim/bench_abuaisha/historiques/e3_iso12.csv"),
    # annexe A de l'article : seule solution fermee (Parker 1981)
    # maillage grossier (hFine = 12 mm, 125 elements sur la fissure) par defaut ;
    # le maillage de production (3 mm, ~22 h a 1 fil) est lourd
    "B4_parker_c": dict(bench="B4", deck="b4_parker.cfg", mesh="parker_crack_c.msh",
                        over={"frames": "4"}, kind="parker"),
    "B4_parker": dict(bench="B4", deck="b4_parker.cfg", mesh="parker_crack.msh",
                      over={"frames": "4"}, kind="parker", heavy=True),
    # B5. Cas de reference sigma0 = 5 MPa, lambda = 1 (out_tun_ref_iso)
    "B5_s5": dict(bench="B5", deck="b5_tunnel_ref_s5_lam1.cfg", mesh="tunnel_hs_iso.msh",
                  over={}, kind="tunnel"),
    # B6. UCS de la fig. 17 (maillage voronoi interne, graine 4211)
    "B6_ucs": dict(bench="B6", deck="b6_ucs_adap.cfg", mesh=None, over={}, kind="ucs"),
    # B7. Banc court Kuru s = 2,5 (loi plastic) du 12-13/09 ; St Anne grossier
    # rock25 lu a 134 us (ECARTS §10) ; St Anne rock137 (29 h sur 14 fils) lourd
    "B7_kuru_s25": dict(bench="B7", deck="b7_yang2026_bench_s25_plastic.cfg",
                        mesh="impact_yang_s2.5_pose.msh", over={}, kind="impact"),
    "B7_stanne_s25": dict(bench="B7", deck="b7_stanne2025_bench_s25_visc0.cfg",
                          mesh="impact_yang_train1_rock25_hxt.msh",
                          over={"T": "1.34e-4", "frames": "4"}, kind="impact"),
    "B7_stanne_137": dict(bench="B7", deck="b7_stanne2025_rock137_visc0.cfg",
                          mesh="impact_yang_train1_rock137_hxt.msh", over={}, kind="impact",
                          heavy=True),
}

# ---- references publiees et resultats anterieurs de rockim -------------------
# « paper » : valeur de l'article (ou du critere analytique qu'il utilise) ;
# « prev »  : resultat anterieur de rockim, avec binaire, plateforme et date.
REF = {
    "B4_aniso": dict(
        paper={"p_peak_MPa": 12.0}, paper_note="critere de Hubbert-Willis de leur eq. 10, "
        "3 sigma_h - sigma_H + ft = 12,0 MPa ; leur numerique ~12,5 MPa (texte), 11,69 lu sur la fig. 11b",
        prev={"p_peak_MPa": 14.993, "t_peak_ms": 2.928},
        prev_note="rockim_e1.exe (correctif de signe du 20/08), Windows/MSVC, 2026-08-20"),
    "B4_iso": dict(
        paper={"p_peak_MPa": 14.2}, paper_note="2 sigma + ft = 14,2 MPa (Kirsch isotrope)",
        prev={"p_peak_MPa": 16.078, "t_peak_ms": 3.177},
        prev_note="rockim_e1.exe, Windows/MSVC, 2026-08-20"),
    "B4_e3_aniso": dict(
        paper={"p_peak_MPa": 12.0}, paper_note="idem B4_aniso, protocole de leur §3.2 (hydroStart)",
        prev={"p_peak_MPa": 14.999, "t_peak_ms": 3.491, "p_ins_MPa": 13.324},
        prev_note="rockim_e3.exe, Windows/MSVC, 2026-08-22"),
    "B4_e3_iso12": dict(
        paper={"p_peak_MPa": 12.0}, paper_note="cible analytique 2 sigma + ft = 12,0 MPa "
        "(etat de contrainte NON publie, choisi pour egaler la cible anisotrope)",
        prev={"p_peak_MPa": 13.755, "t_peak_ms": 3.225, "p_ins_MPa": 12.596},
        prev_note="rockim_e3.exe, Windows/MSVC, 2026-08-22"),
    "B4_parker": dict(
        paper={"w0_mm": 0.0640}, paper_note="Parker (1981), leur eq. A.1 : "
        "w(0) = 2 sigma' (1 - nu^2) l / E, avec E = 45 GPa (coquille « 45 MPa » du texte)",
        prev={}, prev_note="aucune valeur posterieure au correctif de signe consignee ; "
        "les figures Parker anterieures sont d'avant le correctif"),
    "B4_parker_c": dict(
        paper={"w0_mm": 0.0640}, paper_note="idem B4_parker, maillage grossier",
        prev={}, prev_note="aucune valeur posterieure au correctif de signe consignee"),
    "B5_s5": dict(
        paper={"edz_max_m": 19.0, "u_max_m": 0.347},
        paper_note="leur fig. 12f et 11 (MultiFracS) ; hierarchie cisaillement > mixte > traction",
        prev={"edz_max_m": 17.191, "edz_p95_m": 14.486, "u_max_m": 0.389, "u_wall_mean_m": 0.069,
              "broken": 14935, "tensile": 4817, "mixed": 2371, "shear": 7747},
        prev_note="rockim_tun.exe (patchs in situ + excavation), Windows/MSVC, 2026-08-17, "
        "out_tun_ref_iso"),
    "B6_ucs": dict(
        paper={"ucs_MPa": 51.0},
        paper_note="leur fig. 17-18 : UCS ~51 MPa lu sur la figure, module d'entree 15 GPa retrouve",
        prev={"ucs_MPa": 51.0735, "broken": 327, "inserted": 1288},
        prev_note="reperes de tools/verify_suite.py (ucs_yan_adaptive, baseline Linux) ; "
        "MSVC donne 313 / 1 308 aux memes comptages"),
    "B7_kuru_s25": dict(
        paper={"sig_MPa": 160.0, "v_ind": 5.62, "p_max_mm": 1.0, "t_turn_us": 255.0},
        paper_note="Yang et al. 2026, Kuru 9 m/s (lectures de figures, +-10 %)",
        prev={"sig_MPa": 174.0, "v_ind_inst": 6.07, "p_max_mm": 0.87, "t_turn_us": 270.0,
              "broken": 32},
        prev_note="rockim_g1y10/g1y11, Windows/MSVC 10 fils, 2026-09-13 (ETAT §13.4)"),
    "B7_stanne_s25": dict(
        paper={}, paper_note="banc grossier : ne se compare pas a l'article (ECARTS §10)",
        prev={"F_max_kN": 39.5, "impulse_Ns": 1.495},
        prev_note="rockim_g1y18, Windows/MSVC, 2026-09-14, a 134 us"),
    "B7_stanne_137": dict(
        paper={"sig_MPa": 200.0, "v_ind": 8.0, "p_max_mm": 1.10, "t_turn_us": 254.0,
               "R_crater_mm": 13.5, "frag_mass_g": 2.5},
        paper_note="Yang et al. 2025 FDEM, fig. 9-10 lues (+-10 %) ; fin de charge 254 us (texte)",
        prev={"sig_MPa": 212.3, "v_ind": 7.34, "p_max_mm": 1.292, "t_turn_us": 284.3,
              "R_crater_mm": 12.0, "frag_mass_g": 4.70, "broken": 12090},
        prev_note="rockim_g1y19.exe, Windows/MSVC 14 fils, 105 690 s, 2026-09-14"),
}

# ---- criteres d'acceptation (fixes AVANT le calcul, 2026-10-04) --------------
# Critere principal d'un re-run : reproduire le resultat ANTERIEUR de rockim.
# Le maillage est regenere par Gmsh 4.15.2 (macOS) et non par la version Windows
# d'origine : sa realisation differe (hf_bore 189 403 triangles ; rock137
# 109 160 tetraedres contre 108 667), et l'effet d'une graine de maillage n'a
# jamais ete mesure sur ces bancs. Tolerances :
#   - grandeurs AVANT rupture (onde elastique, pic de pression qui precede la
#     propagation) : 2 % ;
#   - grandeurs de cinematique d'impact (enfoncement, pentes) : 5 % ;
#   - grandeurs POST-PIC ou de fissuration (rayon d'EDZ, comptages, cratere,
#     masse de fragments) : 10 a 15 %, le chaos de plateforme y est etabli
#     (UCS de Yan : 327 / 1 288 sous Linux contre 313 / 1 308 sous MSVC, soit
#     4 % ; ucs_yan_adaptive ne verrouille que le pic a 0,15 MPa).
# Critere secondaire, A POSTERIORI pour ces bancs anciens : l'ecart a l'article
# reste celui annonce (meme signe, a 2 points pres pour B4 et B6, 5 points pour
# B5 et B7). Ces ecarts ont ete constates avant toute ecriture de critere : ils
# ne valent pas validation, seulement constance.
TOL = {
    "B4_aniso": {"p_peak_MPa": 0.02, "t_peak_ms": 0.05},
    "B4_iso": {"p_peak_MPa": 0.02, "t_peak_ms": 0.05},
    "B4_e3_aniso": {"p_peak_MPa": 0.02, "p_ins_MPa": 0.02, "t_peak_ms": 0.05},
    "B4_e3_iso12": {"p_peak_MPa": 0.02, "p_ins_MPa": 0.02, "t_peak_ms": 0.05},
    "B5_s5": {"edz_max_m": 0.10, "edz_p95_m": 0.10, "u_wall_mean_m": 0.10, "broken": 0.10},
    "B6_ucs": {"ucs_MPa": 0.15 / 51.0735, "broken": 0.10, "inserted": 0.10},
    "B7_kuru_s25": {"sig_MPa": 0.02, "v_ind_inst": 0.05, "p_max_mm": 0.05, "t_turn_us": 0.05},
    "B7_stanne_s25": {"impulse_Ns": 0.03, "F_max_kN": 0.10},
    "B7_stanne_137": {"sig_MPa": 0.02, "v_ind": 0.05, "p_max_mm": 0.05, "t_turn_us": 0.05,
                      "R_crater_mm": 0.10, "frag_mass_g": 0.15, "broken": 0.15},
}
# ecart a la reference ANALYTIQUE (critere a priori, Parker n'a jamais ete
# mesure apres le correctif de signe) : 10 % (PLAN_ROBUSTESSE du 05/09, G5)
TOL_PAPER_APRIORI = {"B4_parker": {"w0_mm": 0.10}, "B4_parker_c": {"w0_mm": 0.10}}
TOL_PAPER_DRIFT_PTS = {"B4": 2.0, "B5": 5.0, "B6": 2.0, "B7": 5.0}


# ---- maillages, decks ---------------------------------------------------------
def ensure_mesh(name):
    if name is None:
        return None
    p = os.path.join(MESH, name)
    if os.path.exists(p):
        return p
    os.makedirs(MESH, exist_ok=True)
    src = os.path.join(ROOT, "meshes", name)
    if name == "impact_yang_s2.5_pose.msh" and os.path.exists(src):
        shutil.copy(src, p)
        return p
    args = [a.replace("{out}", p) for a in MESHES[name]]
    subprocess.run([PY, os.path.join(ROOT, args[0])] + args[1:], check=True, cwd=ROOT)
    if name == "tunnel_hs_iso.msh":
        drop_orphans_2d(p)
    return p


def drop_orphans_2d(path):
    """Le generateur `tunnelhs` laisse 4 noeuds de construction (centres des arcs)
    qu'aucun triangle ne reference ; la garde C3 (w20) du binaire courant refuse
    le maillage (« noeuds orphelins »), ce que rockim_tun.exe ne faisait pas.
    On retire ces noeuds et les elements ponctuels (type 15) qui les portent,
    puis on renumerote ; triangles et segments sont inchanges, dans le meme
    ordre. Le maillage d'origine est garde en <nom>_brut.msh."""
    L = open(path).read().split("\n")
    i, j = L.index("$Nodes"), L.index("$Elements")
    n, m = int(L[i + 1]), int(L[j + 1])
    nodes, elems = L[i + 2:i + 2 + n], L[j + 2:j + 2 + m]
    used = set()
    for ln in elems:
        t = ln.split()
        if t[1] == "2":
            used.update(int(x) for x in t[-3:])
    keep = [int(ln.split()[0]) for ln in nodes if int(ln.split()[0]) in used]
    if len(keep) == n:
        return 0
    ren = {k: r + 1 for r, k in enumerate(keep)}
    nn = [f"{ren[int(ln.split()[0])]} {' '.join(ln.split()[1:])}" for ln in nodes
          if int(ln.split()[0]) in ren]
    ne = []
    for ln in elems:
        t = ln.split()
        ntag = int(t[2])
        conn = t[3 + ntag:]
        if any(int(x) not in ren for x in conn):
            continue
        ne.append(" ".join(t[:3 + ntag] + [str(ren[int(x)]) for x in conn]))
    ne = [f"{r + 1} {' '.join(ln.split()[1:])}" for r, ln in enumerate(ne)]
    shutil.copy(path, path.replace(".msh", "_brut.msh"))
    with open(path, "w") as f:
        f.write("\n".join(L[:i]) + f"\n$Nodes\n{len(nn)}\n" + "\n".join(nn) +
                f"\n$EndNodes\n$Elements\n{len(ne)}\n" + "\n".join(ne) + "\n$EndElements\n")
    print(f"[maillage] {os.path.basename(path)} : {n - len(nn)} noeuds orphelins retires")
    return n - len(nn)


def deck_lines(case, over_extra=None):
    """lignes du deck source, avec meshFile et les surcharges remplacees"""
    c = CASES[case]
    over = dict(c["over"])
    over.update(over_extra or {})
    if c["mesh"]:
        over["meshFile"] = os.path.join(MESH, c["mesh"])
    seen, lines = set(), []
    for ln in open(os.path.join(DECKS, c["deck"]), encoding="utf-8"):
        ln = ln.rstrip("\n")
        m = re.match(r"\s*([A-Za-z_][\w.]*)\s*=", ln)
        if m and m.group(1) == "outputDir":
            continue
        if m and m.group(1) in over:
            k = m.group(1)
            lines.append(f"{k} = {over[k]}    # surcharge b_articles.py (source : {ln.strip()})")
            seen.add(k)
            continue
        lines.append(ln)
    for k, v in over.items():
        if k not in seen:
            lines.append(f"{k} = {v}    # ajoute par b_articles.py")
    return lines


def make_job(case, outroot=OUT, over_extra=None):
    c = CASES[case]
    ensure_mesh(c["mesh"])
    outdir = os.path.join(outroot, case)
    cfg = C.write_deck(outdir, deck_lines(case, over_extra),
                       header=f"{case} ({c['bench']}), ecrit par vv/B_articles/b_articles.py "
                              f"depuis decks/{c['deck']}")
    return C.Job(bench=c["bench"], name=case, cfg=cfg, outdir=outdir,
                 weight=COST_H2.get(case, 1.0), meta=dict(kind=c["kind"]))


def selected(include_heavy=None):
    if include_heavy is None:
        include_heavy = os.environ.get("VV_B_LOURDS", "0") == "1"
    return [k for k, c in CASES.items() if include_heavy or not c.get("heavy")]


def JOBS(args):
    return [make_job(k) for k in selected()]


# ---- fumee ----------------------------------------------------------------------
# T court, 1 trame, 1 fil : le deck demarre-t-il, toutes ses cles sont-elles lues,
# quel est le cout par pas ? Le cout du cas complet en est deduit (pas x nombre
# de pas), sans la croissance du contact en phase de fracture (a majorer, voir
# README).
SMOKE_T = {"hydro": "4e-6", "parker": "4e-6", "tunnel": "4e-4", "ucs": "2e-6", "impact": "3e-7"}
# cout estime A 2 FILS (heures), fumee du 2026-10-04 corrigee du temps de mise en
# place (runs a 3 pas), cout par pas avant fissuration, divise par 1,6 pour le
# second fil ; la phase de fracture et le contact le font croitre (facteur 2,9
# mesure sur le run St Anne rock137 de Windows). Sert aussi de poids a la file.
COST_H2 = {"B4_aniso": 2.2, "B4_iso": 2.3, "B4_e3_aniso": 2.5, "B4_e3_iso12": 2.3,
           "B4_parker_c": 0.35, "B4_parker": 11.0, "B5_s5": 0.6, "B6_ucs": 0.01,
           "B7_kuru_s25": 1.0, "B7_stanne_s25": 0.4, "B7_stanne_137": 9.0}


def smoke(cases, timeout=90):
    res = {}
    for k in cases:
        c = CASES[k]
        job = make_job(k, OUT_SMOKE, {"T": SMOKE_T[c["kind"]], "frames": "1"})
        rc, wall = C.run_job(job, threads=1, force=True) if timeout is None else \
            _run_timeout(job, timeout)
        txt = open(os.path.join(job.outdir, "rockim.log"), errors="replace").read()
        info = _smoke_info(txt)
        info.update(rc=rc, wall_s=round(wall, 2))
        Tfull = float(_deck_value(k, "T"))
        if info.get("dt") and info.get("steps") and info.get("solver_s"):
            per_step = info["solver_s"] / info["steps"]
            nfull = Tfull / info["dt"]
            info.update(T_full=Tfull, steps_full=int(nfull), s_per_step_1t=per_step,
                        cost_h_1t=per_step * nfull / 3600.0)
        res[k] = info
        print(f"  {k:15s} rc = {rc}  {wall:6.1f} s  dt = {info.get('dt')}  pas = {info.get('steps')}"
              f"  cles non lues = {info.get('unread')}  cout 1 fil ~ {info.get('cost_h_1t', float('nan')):.2f} h",
              flush=True)
    os.makedirs(OUT_SMOKE, exist_ok=True)
    old = {}
    p = os.path.join(OUT_SMOKE, "smoke.json")
    if os.path.exists(p):
        old = json.load(open(p))
    old.update(res)
    json.dump(old, open(p, "w"), indent=1)
    return res


def _run_timeout(job, timeout):
    env = dict(os.environ, OMP_NUM_THREADS="1")
    t0 = time.time()
    with open(os.path.join(job.outdir, "rockim.log"), "w") as log:
        try:
            rc = subprocess.run([C.EXE_DEFAULT, job.cfg], stdout=log, stderr=subprocess.STDOUT,
                                cwd=ROOT, env=env, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = "timeout"
    return rc, time.time() - t0


def _smoke_info(txt):
    def f(p, cast=float):
        m = re.search(p, txt)
        return cast(m.group(1).rstrip(",")) if m else None
    return dict(
        dt=f(r"dt\s*=\s*([\d.eE+-]+)\s*s"),
        steps=f(r"steps\s*=\s*(\d+)", int),
        solver_s=f(r"wall time:\s*([\d.eE+-]+)\s*s"),
        consumed=f(r"cles : (\d+) consommees", int),
        unread=f(r"(\d+) du deck non lues", int),
        errors=[l for l in txt.splitlines() if re.search(r"ERREUR|ERROR|inconnue|unknown key", l, re.I)][:5],
        warnings=len(re.findall(r"WARNING|AVERTISSEMENT", txt)),
        elems=f(r"(\d+) tets", int) or f(r"\] (\d+) elements,", int),
        joints=f(r"(\d+) joints", int),
    )


def _deck_value(case, key):
    for ln in deck_lines(case):
        m = re.match(rf"\s*{re.escape(key)}\s*=\s*([^#\s]+)", ln)
        if m:
            return m.group(1)
    return None




# ---- depouillement ----------------------------------------------------------
def m_hydro(outdir):
    h = C.read_history(outdir)
    p = h["hydroP"] / 1e6
    i = int(np.argmax(p))
    r = dict(p_peak_MPa=float(p[i]), t_peak_ms=float(h["t"][i] * 1e3))
    if "nInserted" in h.dtype.names and h["nInserted"].max() > 0:
        r["p_ins_MPa"] = float(p[np.argmax(h["nInserted"] > 0)])
    if h["nBroken"].max() > 0:
        r["p_break1_MPa"] = float(p[np.argmax(h["nBroken"] > 0)])
    r["peak_inside_run"] = bool(i < len(p) - 5)
    return r


def m_parker(outdir, l=0.75, cx=4.0, cy=4.0):
    """ouverture SIGNEE au centre : y(levre du haut) - y(levre du bas), chaque
    copie de noeud etant rattachee au cote de son element (en FDEM chaque
    element a ses propres noeuds). Pas de max - min : un controle de signe qui
    passe par une norme ne controle pas le signe (VALIDATION_hydro §2)."""
    sys.path.insert(0, os.path.join(ROOT, "tunnel_edz"))
    from plot_tunnel_fields import read_vtu, complete
    el = [f for f in sorted(glob.glob(os.path.join(outdir, "fdem_[0-9]*.vtu"))) if complete(f)]
    P0, Cn, _ = read_vtu(el[0], [])
    P, _, _ = read_vtu(el[-1], [])
    side = np.zeros(len(P0))
    cyel = P0[Cn].mean(axis=1)[:, 1]
    for k in range(3):
        side[Cn[:, k]] = np.sign(cyel - cy)
    on = (np.abs(P0[:, 1] - cy) < 1e-9) & (np.abs(P0[:, 0] - cx) <= l + 1e-9)
    xs = np.round(P0[on, 0] - cx, 9)
    ys, ss = P[on, 1], side[on]
    xa, wa = [], []
    for xv in np.unique(xs):
        g = xs == xv
        up, lo = ys[g & (ss > 0)], ys[g & (ss < 0)]
        if len(up) and len(lo):
            xa.append(xv)
            wa.append(up.mean() - lo.mean())
    xa, wa = np.array(xa), np.array(wa)
    w_an = lambda x: 2 * 2e6 * (1 - 0.2 ** 2) / 45e9 * np.sqrt(np.maximum(l * l - x * x, 0))
    i0 = int(np.argmin(np.abs(xa)))
    sel = np.abs(xa) < 0.9 * l
    return dict(w0_mm=float(wa[i0] * 1e3), w0_parker_mm=float(w_an(0.0) * 1e3),
                err_mean_pct=float(np.mean((wa[sel] - w_an(xa[sel])) / w_an(xa[sel])) * 100),
                n_pairs=int(len(xa)), profile=dict(x=xa.tolist(), w=wa.tolist()))


def m_tunnel(outdir):
    js = os.path.join(outdir, "edz_metrics.json")
    subprocess.run([PY, os.path.join(ROOT, "tunnel_edz/tools/edz_metrics.py"), outdir,
                    "--json", js], check=True, cwd=ROOT, stdout=subprocess.DEVNULL)
    e = json.load(open(js))
    sys.path.insert(0, os.path.join(ROOT, "tunnel_edz/tools"))
    import wall_convergence as wc
    cv = wc.convergence(outdir)
    return dict(edz_max_m=e["edz_radius_max_m"], edz_p95_m=e["edz_radius_p95_m"],
                u_max_m=e["u_max_m"], u_wall_mean_m=float(cv[0]) if cv else None,
                broken=e["broken"], tensile=e["tensile"], mixed=e["mixed"], shear=e["shear"],
                hierarchy=">".join(k for k, _ in sorted(
                    [("cisaillement", e["shear"]), ("mixte", e["mixed"]), ("traction", e["tensile"])],
                    key=lambda kv: -kv[1])))


def m_ucs(outdir):
    txt = open(os.path.join(outdir, "rockim.log"), errors="replace").read()
    g = lambda p, c=float: (c(re.search(p, txt).group(1)) if re.search(p, txt) else None)
    return dict(ucs_MPa=g(r"peak axial stress = (-?[\d.eE+]+) MPa"),
                broken=g(r"broken joints\s*[:=]\s*(\d+)", int),
                inserted=g(r"adaptive insertion: (\d+) / \d+ joints inserted", int),
                gcfric=g(r"dont frottement (-?[\d.eE+-]+) J"),
                E_app_GPa=_e_app(outdir))


def _e_app(outdir):
    """module apparent : pente sigma / epsGauge (extensometre) entre 25 et 50 %
    du pic, branche montante"""
    h = C.read_history(outdir)
    s, e = h["sigma"], h["epsGauge"]
    ip = int(np.argmax(s))
    sel = (s[:ip] > 0.25 * s[ip]) & (s[:ip] < 0.50 * s[ip])
    if sel.sum() < 3:
        return None
    return float(np.polyfit(e[:ip][sel], s[:ip][sel], 1)[0] / 1e9)


def m_impact(outdir):
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import yang_estimators as ye
    h = ye.load_history(outdir)
    k = ye.kinetics(h)
    b = k["bodies"].get("bit") or {}
    r = {}
    if k["gauge"] and k["gauge"]["ok"]:
        r["sig_MPa"] = k["gauge"]["sig"] / 1e6
        r["sig_max_MPa"] = k["gauge"]["sig_max"] / 1e6
    if b:
        r["p_max_mm"] = b["ind"]["p_max"] * 1e3
        r["t_turn_us"] = b["ind"]["t_pmax"] * 1e6
        if b["ind"]["ok"]:
            r["v_ind"] = b["ind"]["v"]
        if b["inst"]:
            r["v_ind_inst"] = abs(float(np.min(h["vz_bit"])))
    t = h["t"]
    fz = [c for c in h if c.startswith("Fc_rock_insert") and c.endswith("z")]
    if fz:
        F = np.abs(h[fz[0]])
        r["F_max_kN"] = float(F.max() / 1e3)
        r["impulse_Ns"] = float(np.sum(0.5 * (F[1:] + F[:-1]) * np.diff(t)))
    for c in ("nBroken",):
        if c in h:
            r["broken"] = int(h[c][-1])
    r["t_end_us"] = float(t[-1] * 1e6)
    if "stanne" in os.path.basename(outdir.rstrip("/")):
        csvp = os.path.join(outdir, "crater_metrics.csv")
        subprocess.run([PY, os.path.join(ROOT, "tools/crater_metrics.py"), outdir, "--sectors", "16",
                        "--csv", csvp], cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if os.path.exists(csvp):
            cm = {ln.split(",")[0]: float(ln.split(",")[1]) for ln in open(csvp).read().split("\n")[1:] if ln}
            r["R_crater_mm"] = cm.get("rCrater", np.nan) * 1e3
            r["R_max_mm"] = cm.get("rMax", np.nan) * 1e3
            r["frag_mass_g"] = cm.get("vBrushed", np.nan) * RHO_STANNE * 1e3
    return r


RHO_STANNE = 2731.0   # Yang et al. 2025, Table 4


METRICS = dict(hydro=m_hydro, parker=m_parker, tunnel=m_tunnel, ucs=m_ucs, impact=m_impact)


def gap(a, b):
    return None if (a is None or b in (None, 0)) else (a - b) / b


def analyse(cases=None):
    cases = cases or list(CASES)
    res = {}
    print("[B] cas  grandeur : article | rockim anterieur | rockim courant | ecart courant/anterieur")
    for k in cases:
        c = CASES[k]
        out = os.path.join(OUT, k)
        ref = REF[k]
        entry = dict(bench=c["bench"], deck=c["deck"], mesh=c["mesh"], paper=ref["paper"],
                     paper_note=ref["paper_note"], prev=ref["prev"], prev_note=ref["prev_note"],
                     cur=None, status="non lance")
        if C.is_done(out):
            lg = C.parse_log(out)
            try:
                cur = METRICS[c["kind"]](out)
            except Exception as ex:   # sortie incomplete : on le dit
                cur = {"erreur": repr(ex)}
            cur["wall_s"] = lg["wall"]
            entry["cur"] = cur
            checks = []
            for q, tol in TOL.get(k, {}).items():
                if q in ref["prev"] and cur.get(q) is not None:
                    g = gap(cur[q], ref["prev"][q])
                    checks.append(dict(q=q, prev=ref["prev"][q], cur=cur[q], gap=g, tol=tol,
                                       ok=abs(g) <= tol, kind="reproduction"))
            for q, tol in TOL_PAPER_APRIORI.get(k, {}).items():
                if cur.get(q) is not None:
                    g = gap(cur[q], ref["paper"][q])
                    checks.append(dict(q=q, paper=ref["paper"][q], cur=cur[q], gap=g, tol=tol,
                                       ok=abs(g) <= tol, kind="article, a priori"))
            for q, pv in ref["paper"].items():
                if q in ref["prev"] and cur.get(q) is not None:
                    gp, gc = gap(ref["prev"][q], pv), gap(cur[q], pv)
                    d = abs(gc - gp) * 100
                    checks.append(dict(q=q, paper=pv, gap_prev=gp, gap_cur=gc, drift_pts=d,
                                       tol_pts=TOL_PAPER_DRIFT_PTS[c["bench"]],
                                       ok=d <= TOL_PAPER_DRIFT_PTS[c["bench"]] and gp * gc >= 0,
                                       kind="ecart a l'article, a posteriori"))
            entry["checks"] = checks
            entry["status"] = "PASSE" if checks and all(x["ok"] for x in checks) else \
                ("ECHEC" if checks else "sans critere")
            for x in checks:
                print(f"  {k:14s} {x['q']:12s} {x['kind']:32s} {'PASSE' if x['ok'] else 'ECHEC'}  "
                      + json.dumps({kk: (round(v, 4) if isinstance(v, float) else v)
                                    for kk, v in x.items() if kk not in ('q', 'kind', 'ok')}))
        else:
            print(f"  {k:14s} non lance ou non fini")
        res[k] = entry
    p = os.path.join(HERE, "resultats.json")
    old = json.load(open(p)) if os.path.exists(p) else {}
    cases_all = old.get("cases", {})
    cases_all.update(res)
    old.update(dict(cases=cases_all, date=time.strftime("%Y-%m-%d"), exe=C.EXE_DEFAULT,
                    smoke=json.load(open(os.path.join(OUT_SMOKE, "smoke.json")))
                    if os.path.exists(os.path.join(OUT_SMOKE, "smoke.json")) else None))
    json.dump(old, open(p, "w"), indent=1, default=float)
    figures(res)
    return res


def figures(res):
    plt = C.plot_style()
    # 1) B4 : p(t) anterieure (historiques archives) et courante
    fig, ax = plt.subplots(figsize=(6.0, 3.2))
    any_curve = False
    for k, col in (("B4_aniso", "C0"), ("B4_iso", "C3"), ("B4_e3_aniso", "C2"), ("B4_e3_iso12", "C1")):
        pc = os.path.join(PHD, CASES[k]["prev_csv"])
        if os.path.exists(pc):
            d = np.genfromtxt(pc, delimiter=",", names=True)
            ax.plot(d["t"] * 1e3, d["hydroP"] / 1e6, color=col, lw=0.8, ls="--",
                    label=f"{k}, antérieur")
            any_curve = True
        out = os.path.join(OUT, k)
        if os.path.exists(os.path.join(out, "history.csv")):
            d = C.read_history(out)
            ax.plot(d["t"] * 1e3, d["hydroP"] / 1e6, color=col, lw=1.2, label=f"{k}, courant")
    for v, lab in ((12.0, "12,0 MPa (éq. 10)"), (14.2, "14,2 MPa")):
        ax.axhline(v, color="0.5", lw=0.6, ls=":")
        ax.text(0.05, v + 0.15, lab, fontsize=7, color="0.4")
    ax.set_xlabel("t (ms)")
    ax.set_ylabel("pression du fluide (MPa)")
    if any_curve:
        ax.legend(fontsize=6.5, frameon=False, loc="center left")
    fig.tight_layout()
    C.savefig(fig, os.path.join(HERE, "fig_b4_pression"))
    plt.close(fig)
    # 2) ecarts a l'article : anterieur et courant, par grandeur
    rows = []
    for k, e in res.items():
        for q, pv in e["paper"].items():
            gp = gap(e["prev"].get(q), pv)
            gc = gap((e["cur"] or {}).get(q), pv)
            if gp is not None or gc is not None:
                rows.append((f"{k}\n{q}", gp, gc))
    if rows:
        fig, ax = plt.subplots(figsize=(6.4, 0.28 * len(rows) + 1.0))
        y = np.arange(len(rows))
        ax.barh(y + 0.2, [100 * (r[1] if r[1] is not None else np.nan) for r in rows], 0.38,
                color="0.65", label="rockim antérieur")
        ax.barh(y - 0.2, [100 * (r[2] if r[2] is not None else np.nan) for r in rows], 0.38,
                color="C0", label="rockim courant")
        ax.set_yticks(y)
        ax.set_yticklabels([r[0].replace("\n", " ") for r in rows], fontsize=6.5)
        ax.axvline(0, color="k", lw=0.6)
        ax.set_xlabel("écart à la valeur publiée (%)")
        ax.invert_yaxis()
        ax.legend(fontsize=7, frameon=False)
        fig.tight_layout()
        C.savefig(fig, os.path.join(HERE, "fig_b_ecarts"))
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["prepare", "smoke", "run", "analyse", "all"])
    ap.add_argument("--cas", nargs="+", default=None)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--lourds", action="store_true", help="inclut les cas lourds")
    ap.add_argument("--timeout", type=float, default=90.0, help="fumee : borne par cas (s)")
    a = ap.parse_args()
    cases = a.cas or selected(a.lourds or None)
    if a.action == "prepare":
        for k in cases:
            j = make_job(k)
            print(f"  {k:15s} -> {os.path.relpath(j.cfg, ROOT)}")
    if a.action == "smoke":
        smoke(a.cas or list(CASES), a.timeout)
    if a.action in ("run", "all"):
        for k in cases:
            j = make_job(k)
            rc, wall = C.run_job(j, threads=a.threads)
            print(f"  {k:15s} rc = {rc}  {wall:8.1f} s", flush=True)
    if a.action in ("analyse", "all"):
        analyse(a.cas)


if __name__ == "__main__":
    main()
