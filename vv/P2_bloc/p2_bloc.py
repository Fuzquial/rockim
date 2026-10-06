#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# P2.1 : bloc glissant jusqu'a l'arret (validation physique du frottement).
#
# Cas de Xiang, Munjiza, Latham et Guises (Eng. Comput. 26, 2009, §3) repris en
# 3D par Fukuda et al. (RMRE 53, 2020, fig. 9). Un cube de cote a pose sur une
# dalle fixe recoit une vitesse horizontale initiale v0 ; la pesanteur g (selon
# -z dans fdem3d) et le frottement de Coulomb mu le freinent. Solution exacte
# tant que le bloc glisse sans basculer (mu < a/a = 1) :
#   x(t) = v0 t - mu g t^2 / 2   pour t <= ts = v0 / (mu g)
#   L    = v0^2 / (2 mu g)        distance d'arret
#
# Deux corps nommes du maillage (`base`, `bloc`), chacun continu
# (groupContinuum.<corps> = true : aucun joint), sans liaison entre eux : le
# seul lien est le contact general. Mesures dans history.csv par les charges
# nulles force.bloc = 0 0 0 (deplacement moyen des copies du bloc, U_bloc_*).
#
#   python3 vv/P2_bloc/p2_bloc.py all --exe build_nofma/rockim
# ---------------------------------------------------------------------------
import argparse, json, math, os, re, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, "out")

# ---- probleme (fixe AVANT le calcul) ---------------------------------------
A_BLOC = 0.05                          # cube de 50 mm (Xiang : carre de 0,05 m)
LB, WB, HB = 0.60, 0.10, 0.02          # dalle 600 x 100 x 20 mm
H_MESH = 0.0125                        # taille de maille
RHO, E, NU = 2650.0, 1e9, 0.25         # E = 1 GPa comme Xiang (pas de temps raisonnable)
MU, G = 0.5, 9.81
V0S = [0.5, 1.0, 2.0]                  # m/s : L = 25,5 / 102 / 408 mm
CONTACTS = ["penalty", "potential"]

CRIT = dict(L_err=0.02, x_l2=0.02, tip_rot=0.01)   # 2 % sur L et sur x(t)


def L_exact(v0):
    return v0 ** 2 / (2 * MU * G)


def x_exact(t, v0):
    ts = v0 / (MU * G)
    tt = np.minimum(t, ts)
    return v0 * tt - 0.5 * MU * G * tt ** 2


def mesh_path():
    return os.path.join(HERE, "meshes", f"bloc_dalle_h{H_MESH * 1e3:g}.msh")


def ensure_mesh():
    p = mesh_path()
    if os.path.exists(p):
        return p
    import gmsh
    os.makedirs(os.path.dirname(p), exist_ok=True)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("bloc")
    base = gmsh.model.occ.addBox(0, 0, 0, LB, WB, HB)
    x0 = 0.05                                          # bloc a 50 mm du bord
    bloc = gmsh.model.occ.addBox(x0, (WB - A_BLOC) / 2, HB, A_BLOC, A_BLOC, A_BLOC)
    gmsh.model.occ.synchronize()                       # PAS de fragment : interfaces non conformes
    gmsh.model.addPhysicalGroup(3, [base], name="base")
    gmsh.model.addPhysicalGroup(3, [bloc], name="bloc")
    gmsh.option.setNumber("Mesh.MeshSizeMax", H_MESH)
    gmsh.option.setNumber("Mesh.MeshSizeMin", H_MESH / 2)
    gmsh.option.setNumber("Mesh.RandomSeed", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(3)
    gmsh.write(p)
    gmsh.finalize()
    return p


def tag(v0, contact):
    return f"{contact}_v{v0:g}".replace(".", "p")


def write_deck(v0, contact, outdir):
    T = 1.25 * v0 / (MU * G)
    lines = [f"# P2.1 bloc glissant, v0 = {v0} m/s, contact = {contact} (ecrit par p2_bloc.py)",
             "mode = fdem3d", "scenario = loads", "mesh = file", f"meshFile = {mesh_path()}",
             f"T = {T:.6g}", "frames = 4", "historyFlush = true",
             f"rho = {RHO}", f"E = {E}", f"nu = {NU}",
             "ft = 1e12", "cohesion = 1e12", "frictionDeg = 30", "Gf = 1e6",
             "insertion = adaptive",
             "groupContinuum.base = true", "groupContinuum.bloc = true",
             f"contact = {contact}", f"contactMu = {MU}",
             f"gravity = {G}", "energyBodyForces = on",
             "fix.base = all",
             f"groupVel.bloc = {v0} 0 0",
             "force.bloc = 0 0 0",
             f"outputDir = {outdir}"]
    cfg = os.path.join(outdir, "deck.cfg")
    open(cfg, "w").write("\n".join(lines) + "\n")
    return cfg


def run(exe, v0, contact):
    od = os.path.join(OUT, tag(v0, contact))
    os.makedirs(od, exist_ok=True)
    cfg = write_deck(v0, contact, od)
    t0 = time.time()
    with open(os.path.join(od, "rockim.log"), "w") as log:
        rc = subprocess.run([exe, cfg], stdout=log, stderr=subprocess.STDOUT, cwd=ROOT).returncode
    print(f"  {contact:9s} v0 = {v0:4g} m/s  rc = {rc}  {time.time() - t0:7.1f} s", flush=True)


def JOBS(args=None):
    """runs du banc pour vv/run_queue.py (contact par potentiel seul : la
    penalite a echoue au premier run, 0,5 m/s, garde pour la trace)"""
    sys.path.insert(0, os.path.dirname(HERE))
    import vvcommon as C
    ensure_mesh()
    jobs = []
    for v0 in V0S:
        od = os.path.join(OUT, tag(v0, "potential"))
        os.makedirs(od, exist_ok=True)
        jobs.append(C.Job("P2.1", tag(v0, "potential"), write_deck(v0, "potential", od), od,
                          weight=v0))
    return jobs


def analyse():
    res = {}
    for contact in CONTACTS:
        for v0 in V0S:
            od = os.path.join(OUT, tag(v0, contact))
            lg = os.path.join(od, "rockim.log")
            if not os.path.exists(lg) or "wall time" not in open(lg).read():
                continue
            log = open(lg).read()
            h = np.genfromtxt(os.path.join(od, "history.csv"), delimiter=",", names=True)
            t, x = h["t"], h["U_bloc_x"]
            z = h["U_bloc_z"]
            xe = x_exact(t, v0)
            Lsim, Le = float(x[-1]), L_exact(v0)
            m = re.search(r"residu\s*:\s*\S+ J \((\S+) % de l'echelle\)", log)
            r = dict(L_sim=Lsim, L_exact=Le, L_err=Lsim / Le - 1,
                     x_l2=float(np.sqrt(np.sum((x - xe) ** 2) / np.sum(xe ** 2))),
                     z_min=float(z.min()), z_max=float(z.max()),
                     budget_pct=float(m.group(1)) if m else float("nan"),
                     wall=float(re.search(r"wall time: (\S+) s", log).group(1)))
            r["passe"] = abs(r["L_err"]) <= CRIT["L_err"] and r["x_l2"] <= CRIT["x_l2"]
            res[f"{contact}_v{v0:g}"] = r
            print(f"{contact:9s} v0 {v0:4g}  L {Lsim * 1e3:8.2f} mm (exact {Le * 1e3:8.2f})  "
                  f"err {100 * r['L_err']:+6.2f} %  x_L2 {100 * r['x_l2']:5.2f} %  "
                  f"z [{r['z_min'] * 1e6:.1f}, {r['z_max'] * 1e6:.1f}] um  "
                  f"B4 {r['budget_pct']:.1e} %  {'PASSE' if r['passe'] else 'ECHEC'}")
    json.dump(dict(reference=dict(mu=MU, g=G, a=A_BLOC, E=E, rho=RHO, criteres=CRIT),
                   runs=res), open(os.path.join(HERE, "resultats.json"), "w"), indent=1)
    figure()
    return res


def figure():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX Two Text"],
                         "mathtext.fontset": "stix", "font.size": 9})
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    for k, v0 in enumerate(V0S):
        tt = np.linspace(0, 1.25 * v0 / (MU * G), 400)
        ax.plot(tt, x_exact(tt, v0) * 1e3, color=f"C{k}", lw=2.2, alpha=0.35)
        for contact, ls in zip(CONTACTS, ["--", ":"]):
            od = os.path.join(OUT, tag(v0, contact))
            p = os.path.join(od, "history.csv")
            if os.path.exists(p):
                h = np.genfromtxt(p, delimiter=",", names=True)
                ax.plot(h["t"], h["U_bloc_x"] * 1e3, color=f"C{k}", lw=0.9, ls=ls,
                        label=f"v0 = {v0:g} m/s, {contact}")
    ax.set_xlabel("temps (s)")
    ax.set_ylabel("déplacement du bloc (mm)")
    ax.set_title(r"Bloc glissant : trait large exact $x(t) = v_0 t - \mu g t^2/2$", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_bloc.pdf"))
    fig.savefig(os.path.join(HERE, "fig_bloc.png"), dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["run", "analyse", "all"])
    ap.add_argument("--exe", default=os.path.join(ROOT, "build_nofma", "rockim"))
    ap.add_argument("--v0", nargs="+", type=float, default=None)
    ap.add_argument("--contacts", nargs="+", default=None)
    a = ap.parse_args()
    global V0S, CONTACTS
    V0S = a.v0 or V0S
    CONTACTS = a.contacts or CONTACTS
    if a.action in ("run", "all"):
        ensure_mesh()
        for c in CONTACTS:
            for v in V0S:
                run(os.path.abspath(a.exe), v, c)
    if a.action in ("analyse", "all"):
        analyse()


if __name__ == "__main__":
    main()
