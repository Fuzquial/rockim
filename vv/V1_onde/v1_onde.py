#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# V1 : propagation d'une onde de compression dans une barre (verification).
#
# Barre W x W x L, axe z. Faces laterales sur rouleaux (fix.xmin = x, ...) :
# deformation uniaxiale, donc le probleme est EXACTEMENT unidimensionnel pour
# le continu, sans dispersion geometrique (Pochhammer). Celerite exacte :
#   c = sqrt(M / rho),  M = E (1 - nu) / ((1 + nu)(1 - 2 nu))  (module oedometrique)
# Chargement : vitesse imposee V(t) sur `bottom` (z = 0), trapeze lineaire par
# morceaux (amplitude.bottom). Bout `top` (z = L) libre.
#
# Solution de d'Alembert avec reflexions (bout impose = deplacement impose,
# reflexion de signe oppose ; bout libre = reflexion de meme signe) :
#   u(z, t) = sum_n (-1)^n [ D(t - (2nL + z)/c) + D(t - (2(n+1)L - z)/c) ]
#   D(t) = integrale de V de 0 a t (nulle pour t < 0)
#   u_z(0, t) = sum_n (-1)^n [ -V(t - 2nL/c) + V(t - 2(n+1)L/c) ] / c
#   reaction de l'appui sur le solide : RF_z = -M u_z(0, t) A
#
# Grandeurs comparees (sorties de history.csv, scenario = loads) :
#   U_s1_z, U_s2_z, U_s3_z, U_top_z : deplacement axial de 3 stations et du
#     bout libre ; la position EXACTE du noeud retenu par point.<g> est lue
#     dans le journal du solveur (ligne « point.sK : sommet ... repere solveur »)
#   RF_bottom_z : reaction de l'appui charge
# Metriques :
#   e_u  = ||u_sim - u_ex||_2 / ||u_ex||_2 sur [0, T], par station
#   e_F  = idem pour la reaction
#   c_mes : (z3 - z1) / (t3 - t1), t_k = instant ou u franchit la moitie du
#     deplacement final de l'impulsion (D_inf / 2), interpole lineairement
#   F_plateau : moyenne de RF sur le plateau de l'impulsion, contre rho c v0 A
#   ordre observe : pente de log(e_u(s2)) en fonction de log(h)
#
#   python3 vv/V1_onde/v1_onde.py run      --exe build_nofma/rockim [--variants ...] [--h ...]
#   python3 vv/V1_onde/v1_onde.py analyse  [--variants ...]
#   python3 vv/V1_onde/v1_onde.py all      --exe build_nofma/rockim
# Lance depuis la racine rockim. Sorties : vv/V1_onde/out/<variante>_h<h>/,
# vv/V1_onde/resultats.json, vv/V1_onde/fig_*.pdf
# ---------------------------------------------------------------------------
import argparse, json, math, os, re, subprocess, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, "out")

# ---- probleme (fixe AVANT le calcul) ---------------------------------------
W, L = 0.008, 0.4                      # section 8 x 8 mm, longueur 400 mm
RHO, E, NU = 2650.0, 60e9, 0.25
V0 = 0.1                               # m/s : sigma = rho c v0 = 1,38 MPa, elastique
AMP = [(0.0, 0.0), (4e-6, 1.0), (16e-6, 1.0), (20e-6, 0.0)]   # trapeze 4/12/4 us
T = 2.5e-4                             # ~3,3 allers-retours
STATIONS = {"s1": 0.1, "s2": 0.2, "s3": 0.3}
M = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU))
C = math.sqrt(M / RHO)
A = W * W

# ---- variantes ---------------------------------------------------------------
# h en mm. fdem3d coute ~50 fois fem3d par tet : h = 1 mm y est reserve a un
# run long, lance a part.
VARIANTS = {
    "fem3d": dict(h=[4, 2.83, 2, 1.41, 1], keys=["mode = fem3d", "law = elastic"]),
    # adaptatif sans insertion (ft >> sigma) : doit reproduire le continu
    "fdem3d_adaptive": dict(h=[4, 2.83, 2], keys=[
        "mode = fdem3d", "insertion = adaptive",
        "ft = 50e6", "cohesion = 100e6", "frictionDeg = 40", "Gf = 70"]),
    # intrinseque : joints de penalite partout, souplesse ajoutee
    "fdem3d_intrinsic_pf20": dict(h=[4, 2.83, 2], keys=[
        "mode = fdem3d", "insertion = intrinsic", "jointPenaltyFactor = 20",
        "ft = 50e6", "cohesion = 100e6", "frictionDeg = 40", "Gf = 70"]),
}

# ---- criteres d'acceptation (fixes AVANT le calcul) ------------------------
# appliques au maillage le plus fin de chaque variante continue
CRIT = dict(e_u_s2=0.02, c_err=0.01, F_err=0.02, budget_pct=1e-6, order_min=1.0)


def amp(t):
    ts, vs = zip(*AMP)
    return np.interp(t, ts, vs, left=0.0, right=vs[-1])


def V(t):
    return V0 * amp(np.asarray(t))


_TF = np.linspace(0.0, T + 2 * L / C, 400001)
_DF = np.concatenate([[0.0], np.cumsum(0.5 * (V(_TF[1:]) + V(_TF[:-1])) * np.diff(_TF))])


def D(t):
    t = np.asarray(t, dtype=float)
    return np.where(t <= 0.0, 0.0, np.interp(t, _TF, _DF))


def nrefl(t):
    return int(np.max(t) * C / (2 * L)) + 2


def u_exact(z, t):
    u = np.zeros_like(t, dtype=float)
    for n in range(nrefl(t)):
        s = (-1) ** n
        u += s * (D(t - (2 * n * L + z) / C) + D(t - (2 * (n + 1) * L - z) / C))
    return u


def rf_exact(t):
    uz = np.zeros_like(t, dtype=float)
    for n in range(nrefl(t)):
        s = (-1) ** n
        tn, tn1 = t - 2 * n * L / C, t - 2 * (n + 1) * L / C
        uz += s * (-np.where(tn > 0, V(tn), 0.0) + np.where(tn1 > 0, V(tn1), 0.0)) / C
    return -M * uz * A


# ---- maillages, decks, runs -------------------------------------------------
def tag(h):
    return f"{h:g}".replace(".", "p")


def mesh_path(h):
    return os.path.join(HERE, "meshes", f"bar_h{tag(h)}.msh")


def ensure_mesh(h):
    p = mesh_path(h)
    if not os.path.exists(p):
        os.makedirs(os.path.dirname(p), exist_ok=True)
        subprocess.run([sys.executable, os.path.join(ROOT, "tools", "make_unstructured_mesh.py"),
                        "box3dbc", str(W), str(W), str(L), str(h * 1e-3), p, "1"],
                       check=True, cwd=ROOT)
    return p


def write_deck(var, h, outdir):
    lines = [f"# V1 onde, variante {var}, h = {h} mm (ecrit par v1_onde.py)",
             *VARIANTS[var]["keys"],
             "scenario = loads", "mesh = file", f"meshFile = {mesh_path(h)}",
             f"T = {T}", "frames = 1", f"rho = {RHO}", f"E = {E}", f"nu = {NU}",
             "fix.xmin = x", "fix.xmax = x", "fix.ymin = y", "fix.ymax = y",
             f"velocity.bottom = free free {V0}",
             "amplitude.bottom = " + " ".join(f"{t:g} {a:g}" for t, a in AMP)]
    for g, z in STATIONS.items():
        lines += [f"point.{g} = {W / 2} {W / 2} {z}", f"force.{g} = 0 0 0"]
    lines += ["force.top = 0 0 0", f"outputDir = {outdir}"]
    cfg = os.path.join(outdir, "deck.cfg")
    with open(cfg, "w") as f:
        f.write("\n".join(lines) + "\n")
    return cfg


def run(exe, var, h):
    outdir = os.path.join(OUT, f"{var}_h{tag(h)}")
    os.makedirs(outdir, exist_ok=True)
    ensure_mesh(h)
    cfg = write_deck(var, h, outdir)
    t0 = time.time()
    with open(os.path.join(outdir, "rockim.log"), "w") as log:
        rc = subprocess.run([exe, cfg], stdout=log, stderr=subprocess.STDOUT, cwd=ROOT).returncode
    print(f"  {var:24s} h = {h:5g} mm  rc = {rc}  {time.time() - t0:7.1f} s", flush=True)
    return rc


# ---- depouillement ----------------------------------------------------------
def read_run(outdir):
    log = open(os.path.join(outdir, "rockim.log")).read()
    zs = {}
    for g in STATIONS:
        m = re.search(rf"point\.{g} : sommet \d+ a \S+ m du point demande \(\s*(\S+)\s+(\S+)\s+(\S+),", log)
        zs[g] = float(m.group(3))
    m = re.search(r"residu\s*:\s*\S+ J \((\S+) % de l'echelle\)", log)
    budget = float(m.group(1)) if m else float("nan")
    m = re.search(r"dt = (\S+) s", log)
    dt = float(m.group(1).rstrip(",")) if m else float("nan")
    m = re.search(r"(\d+) tets", log)
    ntet = int(m.group(1)) if m else 0
    m = re.search(r"wall time: (\S+) s", log)
    wall = float(m.group(1)) if m else float("nan")
    data = np.genfromtxt(os.path.join(outdir, "history.csv"), delimiter=",", names=True)
    return dict(z=zs, budget=budget, dt=dt, ntet=ntet, wall=wall, h=data)


def l2rel(a, b):
    return float(np.sqrt(np.sum((a - b) ** 2) / np.sum(b ** 2)))


def crossing(t, u, level):
    i = int(np.argmax(u >= level))
    if i == 0:
        return float("nan")
    return float(t[i - 1] + (level - u[i - 1]) * (t[i] - t[i - 1]) / (u[i] - u[i - 1]))


def analyse_run(r):
    h = r["h"]
    t = h["t"]
    res = dict(ntet=r["ntet"], dt=r["dt"], wall=r["wall"], budget_pct=r["budget"], z=r["z"])
    sig = {}
    for g in list(STATIONS) + ["top"]:
        z = r["z"][g] if g in STATIONS else L
        ue = u_exact(z, t)
        us = h[f"U_{g}_z"]
        res[f"e_u_{g}"] = l2rel(us, ue)
        sig[g] = (us, ue)
    Fe, Fs = rf_exact(t), h["RF_bottom_z"]
    res["e_F"] = l2rel(Fs, Fe)
    # premier passage, avant tout retour de reflexion
    half = 0.5 * D(np.array([1.0]))[0]
    tk = {g: crossing(t, h[f"U_{g}_z"], half) for g in STATIONS}
    res["c_mes"] = (r["z"]["s3"] - r["z"]["s1"]) / (tk["s3"] - tk["s1"])
    res["c_err"] = res["c_mes"] / C - 1.0
    w = (t > AMP[1][0] + 1e-6) & (t < AMP[2][0] - 1e-6)
    Fp = RHO * C * V0 * A
    res["F_plateau"] = float(np.mean(Fs[w]))
    res["F_err"] = res["F_plateau"] / Fp - 1.0
    return res, sig, (t, Fs, Fe)


def order(hs, es):
    hs, es = np.log(np.asarray(hs)), np.log(np.asarray(es))
    return float(np.polyfit(hs, es, 1)[0]) if len(hs) >= 2 else float("nan")


def figures(results, sigs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX Two Text"],
                         "mathtext.fontset": "stix", "font.size": 10})
    # 1. signaux du maillage le plus fin de chaque variante
    for var, runs in sigs.items():
        hmin = min(runs)
        sig, (t, Fs, Fe) = runs[hmin]
        fig, ax = plt.subplots(2, 1, figsize=(6.3, 5.2), sharex=True)
        for g, col in zip(list(STATIONS) + ["top"], ["C0", "C1", "C2", "C3"]):
            us, ue = sig[g]
            ax[0].plot(t * 1e6, ue * 1e6, color=col, lw=1.6, alpha=0.45)
            ax[0].plot(t * 1e6, us * 1e6, color=col, lw=0.8, ls="--", label=g)
        ax[0].set_ylabel("déplacement axial (µm)")
        ax[0].legend(ncol=4, fontsize=8, frameon=False)
        ax[1].plot(t * 1e6, Fe, color="k", lw=1.6, alpha=0.45, label="exact")
        ax[1].plot(t * 1e6, Fs, color="C3", lw=0.8, ls="--", label="rockim")
        ax[1].set_ylabel("réaction de l'appui (N)")
        ax[1].set_xlabel("temps (µs)")
        ax[1].legend(fontsize=8, frameon=False)
        fig.suptitle(f"V1, {var}, h = {hmin:g} mm : trait plein exact, tirets rockim", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, f"fig_signaux_{var}.pdf"))
        plt.close(fig)
    # 2. convergence
    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    for var, res in results.items():
        hs = sorted(res["runs"], key=float)
        ax.loglog([float(h) for h in hs], [res["runs"][h]["e_u_s2"] for h in hs], "o-", label=var)
    hh = np.array([1.0, 4.0])
    ax.loglog(hh, 0.01 * (hh / 1.0) ** 2, "k:", lw=0.8, label="pente 2")
    ax.loglog(hh, 0.01 * (hh / 1.0) ** 1, "k--", lw=0.8, label="pente 1")
    ax.set_xlabel("taille de maille h (mm)")
    ax.set_ylabel(r"erreur $L^2$ relative en s2")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_convergence.pdf"))
    plt.close(fig)


def analyse(variants):
    results, sigs = {}, {}
    for var in variants:
        runs, sv = {}, {}
        for h in VARIANTS[var]["h"]:
            od = os.path.join(OUT, f"{var}_h{tag(h)}")
            lg = os.path.join(od, "rockim.log")
            if not os.path.exists(lg) or "wall time" not in open(lg).read():
                continue                      # run absent ou inacheve
            res, sig, F = analyse_run(read_run(od))
            runs[f"{h:g}"] = res
            sv[h] = (sig, F)
        if not runs:
            continue
        hs = sorted(runs, key=float)
        fin = runs[hs[0]]
        verdict = {
            "e_u_s2": fin["e_u_s2"] <= CRIT["e_u_s2"],
            "c_err": abs(fin["c_err"]) <= CRIT["c_err"],
            "F_err": abs(fin["F_err"]) <= CRIT["F_err"],
            "budget": all(abs(r["budget_pct"]) <= CRIT["budget_pct"] for r in runs.values()),
        }
        ordre = order([float(h) for h in hs], [runs[h]["e_u_s2"] for h in hs])
        if len(hs) >= 3:
            verdict["ordre"] = ordre >= CRIT["order_min"]
        results[var] = dict(runs=runs, ordre=ordre, verdict=verdict,
                            passe=all(verdict.values()))
        sigs[var] = sv
    ref = dict(c=C, M=M, F_plateau=RHO * C * V0 * A, L=L, W=W, rho=RHO, E=E, nu=NU,
               v0=V0, amp=AMP, T=T, criteres=CRIT)
    with open(os.path.join(HERE, "resultats.json"), "w") as f:
        json.dump(dict(reference=ref, variantes=results), f, indent=1)
    # table console
    print(f"\nreference : c = {C:.2f} m/s, F_plateau = {RHO * C * V0 * A:.4f} N")
    print(f"{'variante':24s} {'h':>5s} {'ntet':>7s} {'e_u s1':>8s} {'e_u s2':>8s} {'e_u s3':>8s} "
          f"{'e_u top':>8s} {'e_F':>8s} {'c err':>8s} {'F err':>8s} {'B4 %':>9s} {'mur s':>7s}")
    for var, res in results.items():
        for h in sorted(res["runs"], key=lambda x: -float(x)):
            r = res["runs"][h]
            print(f"{var:24s} {h:>5s} {r['ntet']:7d} {r['e_u_s1']:8.4f} {r['e_u_s2']:8.4f} "
                  f"{r['e_u_s3']:8.4f} {r['e_u_top']:8.4f} {r['e_F']:8.4f} {r['c_err']:+8.4f} "
                  f"{r['F_err']:+8.4f} {r['budget_pct']:9.1e} {r['wall']:7.1f}")
        print(f"{'':24s} ordre observe {res['ordre']:.2f} ; verdict "
              f"{'PASSE' if res['passe'] else 'ECHEC'} {res['verdict']}")
    figures(results, sigs)
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["run", "analyse", "all"])
    ap.add_argument("--exe", default=os.path.join(ROOT, "build_nofma", "rockim"))
    ap.add_argument("--variants", nargs="+", default=list(VARIANTS))
    ap.add_argument("--h", nargs="+", type=float, default=None, help="restreint les h")
    a = ap.parse_args()
    if a.h:
        for v in a.variants:
            VARIANTS[v]["h"] = [h for h in VARIANTS[v]["h"] if h in a.h] or a.h
    if a.action in ("run", "all"):
        exe = os.path.abspath(a.exe)
        print(f"[V1] exe = {exe}, c = {C:.2f} m/s, 2L/c = {2 * L / C * 1e6:.2f} us")
        for v in a.variants:
            for h in VARIANTS[v]["h"]:
                run(exe, v, h)
    if a.action in ("analyse", "all"):
        analyse(a.variants)


if __name__ == "__main__":
    main()
