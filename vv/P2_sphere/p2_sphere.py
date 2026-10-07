#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# P2.2 : sphere elastique sur un plan (contact de Hertz, impact et indentation).
#
# Essai du type Fukuda et al. (RMRE 53, 2020, fig. 8) complete par la solution
# de Hertz (Johnson, Contact Mechanics, 1985, §4.2 et §11.4). Une sphere
# elastique de rayon R, lancee normalement a v0 sur un bloc elastique fixe par
# sa base, sans frottement (contactMu = 0) ni pesanteur.
#
# Solution exacte (impact quasi statique de Hertz, bloc = demi-espace) :
#   E*  = 1 / ((1 - nu^2)/E + (1 - nu^2)/E) = E / (2 (1 - nu^2))
#   K   = 4/3 E* sqrt(R)                 F = K delta^(3/2)
#   m delta'' = -K delta^(3/2)           (masse du bloc infinie)
#   delta_max = (15 m v0^2 / (16 E* sqrt(R)))^(2/5)
#   F_max     = K delta_max^(3/2)
#   t_c       = 2,9432 delta_max / v0
#   a_max     = sqrt(R delta_max)        rayon de contact maximal
#   e         = 1 (restitution ; la part rayonnee en ondes elastiques dans un
#               demi-espace vaut lambda ~ 1,04 (v0/c0)^(3/5), Hunter 1957,
#               Johnson §11.4, c0 = sqrt(E/rho) : documentee, pas exigee)
# Indentation quasi statique (variante qs) : sphere rigide (vitesse imposee a
# tous ses noeuds selon z), E*_qs = E / (1 - nu^2), F = 4/3 E*_qs sqrt(R) d^(3/2).
#
# Modele quart (symetrie) : quart de sphere et quart de bloc, plans x = 0 et
# y = 0 tenus par fix.symx = x et fix.symy = y, base du bloc fix.bottom = all.
# Deux corps nommes (`block`, `sphere`), continus (groupContinuum), lies par le
# seul contact general. Mesures : force de contact exercee par le bloc sur la
# sphere (contactForcePairs = block:sphere -> Fc_block_sphere_z), centroide et
# vitesse moyenne de la sphere (trackGroup = sphere -> grpZ, grpVz), bilan B4
# et energies par corps du journal. Toutes les forces et masses du quart sont
# multipliees par 4 avant comparaison.
#
# Raffinement : taille hc = a_max / NA dans une boule de rayon 1,5 a_max autour
# du point de contact (les deux corps), croissance lineaire jusqu'a HFAR.
#
#   python3 vv/P2_sphere/p2_sphere.py prepare            # maillages + decks
#   python3 vv/P2_sphere/p2_sphere.py run   [--variants ..] [--na ..] [--threads 2]
#   python3 vv/P2_sphere/p2_sphere.py analyse
#   python3 vv/P2_sphere/p2_sphere.py all
#   python3 vv/P2_sphere/p2_sphere.py smoke [--tfrac 0.6]  # essai court, out/_smoke
#   python3 vv/run_queue.py P2_sphere --slots 4 --threads 2
# ---------------------------------------------------------------------------
import argparse, json, math, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
VV = os.path.dirname(HERE)
ROOT = os.path.dirname(VV)
OUT = os.path.join(HERE, "out")
sys.path.insert(0, VV)
import vvcommon as C  # noqa: E402

BENCH = "P2.2"

# ---- probleme (fixe AVANT le calcul) ---------------------------------------
# Materiau unique pour les deux corps. E = 1 GPa (comme P2.1) : seul le rapport
# rho v0^2 / E fixe a_max / R (0,11 ici) ; un module de roche (50 GPa) a v0 = 1 m/s
# donnerait a_max / R = 0,03 et un maillage 4 fois plus fin au contact.
RHO, E, NU = 2650.0, 1.0e9, 0.25
R = 0.020                       # rayon de la sphere [m]
V0 = 1.0                        # vitesse d'impact [m/s]
GAP = 2.0e-5                    # jeu initial sphere-bloc [m] (20 us de vol)
LB, HB = 3 * R, 3 * R           # quart de bloc LB x LB x HB (bloc entier 6R x 6R x 3R)
HFAR = 4.0e-3                   # taille loin du contact [m]
RFINE = 1.5                     # rayon de la zone fine, en multiples de a_max
GROWTH = 0.3                    # pente de croissance de la taille (dh/dr)
VQS = 0.2                       # vitesse d'indentation quasi statique [m/s]
TRAMP_QS = 2.0e-4               # montee en cosinus de la vitesse qs [s]
FTHR = 0.01                     # seuil de force (fraction de F_max exact) pour t_c force
FWIN = 0.05                     # largeur du filtre de force, en multiples de t_c (ajout
                                # apres l'essai court : oscillations de la force brute)

ESTAR = E / (2 * (1 - NU ** 2))
ESTAR_QS = E / (1 - NU ** 2)
MASS = RHO * 4.0 / 3.0 * math.pi * R ** 3
KH = 4.0 / 3.0 * ESTAR * math.sqrt(R)
DMAX = (15 * MASS * V0 ** 2 / (16 * ESTAR * math.sqrt(R))) ** 0.4
FMAX = KH * DMAX ** 1.5
TC = 2.9432 * DMAX / V0
AMAX = math.sqrt(R * DMAX)
C0 = math.sqrt(E / RHO)
HUNTER = 1.04 * (V0 / C0) ** 0.6

# ---- variantes ---------------------------------------------------------------
# na = a_max / hc (nombre de mailles sur le rayon de contact maximal)
VARIANTS = {
    # variante principale : potentiel de Munjiza, penalite par defaut (p = E)
    "pot_pf1": dict(na=[2, 4, 8, 16], kind="impact",
                    keys=["contact = potential", "potPenaltyFactor = 1"]),
    # penalite du potentiel x10 : part de la souplesse de contact dans l'ecart
    "pot_pf10": dict(na=[4, 8], kind="impact",
                     keys=["contact = potential", "potPenaltyFactor = 10"]),
    # pour memoire : contact de penalite noeud-face (echoue P2.1), restitution 1
    "pen": dict(na=[2], kind="impact",
                keys=["contact = penalty", "gcRestitution = 1.0"]),
    # indentation quasi statique, sphere rigide (vitesse imposee)
    "qs_pf1": dict(na=[2, 4], kind="qs",
                   keys=["contact = potential", "potPenaltyFactor = 1"]),
}

# ---- criteres d'acceptation (fixes AVANT le calcul, 2026-10-03) ------------
# appliques au maillage le plus fin de chaque variante d'impact au potentiel ;
# la variante `pen` et les mailles grossieres sont documentees sans verdict.
# Les SEUILS n'ont pas bouge depuis l'essai court (2026-10-04, a/h = 2,
# T = 0,6 t_c) ; seules les DEFINITIONS ont ete precisees apres lui : F_max lu
# sur la force filtree (FWIN), t_c lu sur le retour a zero de delta.
CRIT = dict(
    Fmax_err=0.05,        # |F_max / F_max_Hertz - 1| (force filtree, exact filtre de meme)
    tc_err=0.05,          # |t_c / t_c_Hertz - 1| (delta repasse par zero)
    dmax_err=0.05,        # |delta_max / delta_max_Hertz - 1|
    e_min=0.90,           # restitution documentee ; critere lache (ondes, maillage)
    budget_pct=1e-6,      # residu B4 du journal, % de l'echelle
    closure=0.01,         # |KE0 - (KE fin des corps + elastique stocke)| / KE0
    qs_F_l2=0.05,         # indentation : ecart L2 de F(d) sur d in [0,2 ; 1] d_max
    monotone=True,        # l'erreur sur F_max decroit quand on raffine
)


# ---- solution exacte ---------------------------------------------------------
def hertz_curve(n=20000):
    """integration de m d'' = -K d^(3/2) depuis d = 0, d' = v0 (RK4), jusqu'a
    la separation ; renvoie t (depuis le premier contact), d, F"""
    dt = 1.2 * TC / n
    t, d, v = [0.0], [0.0], [V0]

    def acc(x):
        return -KH * max(x, 0.0) ** 1.5 / MASS

    x, u = 0.0, V0
    for i in range(n):
        k1x, k1v = u, acc(x)
        k2x, k2v = u + 0.5 * dt * k1v, acc(x + 0.5 * dt * k1x)
        k3x, k3v = u + 0.5 * dt * k2v, acc(x + 0.5 * dt * k2x)
        k4x, k4v = u + dt * k3v, acc(x + dt * k3x)
        x += dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        u += dt / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
        t.append((i + 1) * dt); d.append(x); v.append(u)
        if x < 0:
            break
    t, d = np.array(t), np.array(d)
    return t, np.maximum(d, 0.0), KH * np.maximum(d, 0.0) ** 1.5


def smooth(t, y, win=None):
    """moyenne glissante centree de largeur win (s) sur une grille uniforme"""
    win = FWIN * TC if win is None else win
    tg = np.linspace(t[0], t[-1], max(len(t), 2000))
    yg = np.interp(tg, t, y)
    n = max(1, int(round(win / (tg[1] - tg[0]))))
    k = np.ones(n) / n
    ys = np.convolve(np.pad(yg, (n // 2, n - 1 - n // 2), mode="edge"), k, mode="valid")
    return np.interp(t, tg, ys)


def contact_window(t, F, thr):
    """instants de montee et de descente de F au seuil thr (interpolation)"""
    on = np.where(F >= thr)[0]
    if len(on) == 0:
        return float("nan"), float("nan")
    i0, i1 = on[0], on[-1]

    def cross(a, b):
        return t[a] + (thr - F[a]) * (t[b] - t[a]) / (F[b] - F[a])

    t0 = cross(i0 - 1, i0) if i0 > 0 else t[i0]
    t1 = cross(i1, i1 + 1) if i1 + 1 < len(t) else float("nan")
    return float(t0), float(t1)


_TH, _DH, _FH = hertz_curve()
_T0H, _T1H = contact_window(_TH, _FH, FTHR * FMAX)
TC_THR = _T1H - _T0H              # duree exacte au meme seuil que la mesure


# ---- maillage ----------------------------------------------------------------
def hc_of(na):
    return AMAX / na


def mesh_path(na):
    return os.path.join(HERE, "meshes", f"sphere_plan_na{na}.msh")


def ensure_mesh(na):
    p = mesh_path(na)
    if os.path.exists(p):
        return p
    import gmsh
    os.makedirs(os.path.dirname(p), exist_ok=True)
    hc = hc_of(na)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("sphere_plan")
    occ = gmsh.model.occ
    blk = occ.addBox(0, 0, -HB, LB, LB, HB)
    sph = occ.addSphere(0, 0, GAP + R, R)
    cut = occ.addBox(0, 0, GAP, 1.01 * R, 1.01 * R, 2.02 * R)
    quart, _ = occ.intersect([(3, sph)], [(3, cut)])
    occ.synchronize()                  # PAS de fragment : deux corps disjoints
    sq = quart[0][1]
    gmsh.model.addPhysicalGroup(3, [blk], name="block")
    gmsh.model.addPhysicalGroup(3, [sq], name="sphere")
    eps = 1e-6                         # tolerance OCC des boites englobantes
    symx, symy, bot = [], [], []
    for dim, tg in gmsh.model.getEntities(2):
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(dim, tg)
        if abs(x0) < eps and abs(x1) < eps:
            symx.append(tg)
        elif abs(y0) < eps and abs(y1) < eps:
            symy.append(tg)
        elif abs(z0 + HB) < eps and abs(z1 + HB) < eps:
            bot.append(tg)
    assert len(symx) == 2 and len(symy) == 2 and len(bot) == 1, (symx, symy, bot)
    gmsh.model.addPhysicalGroup(2, symx, name="symx")
    gmsh.model.addPhysicalGroup(2, symy, name="symy")
    gmsh.model.addPhysicalGroup(2, bot, name="bottom")
    # champ de taille : distance au point de contact (0, 0, GAP/2)
    r1 = RFINE * AMAX
    r2 = r1 + (HFAR - hc) / GROWTH
    f = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(
        f, "F", "%.9g + (%.9g) * Max(0, Min(1, (Sqrt(x^2 + y^2 + (z - %.9g)^2) - %.9g) / %.9g))"
        % (hc, HFAR - hc, GAP / 2, r1, r2 - r1))
    gmsh.model.mesh.field.setAsBackgroundMesh(f)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", hc)
    gmsh.option.setNumber("Mesh.MeshSizeMax", HFAR)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.option.setNumber("Mesh.RandomSeed", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(3)
    gmsh.write(p)
    gmsh.finalize()
    return p


# ---- decks et jobs -----------------------------------------------------------
def tag(var, na):
    return f"{var}_na{na}"


def t_end(kind):
    if kind == "impact":
        return GAP / V0 + 1.4 * TC         # contact + rebond libre (~0,4 t_c)
    return (GAP + DMAX) / VQS + 0.5 * TRAMP_QS


def deck_lines(var, na, T):
    v = VARIANTS[var]
    lines = ["mode = fdem3d", "scenario = loads", "mesh = file", f"meshFile = {mesh_path(na)}",
             f"T = {T:.6g}", "frames = 4", "historyFlush = true",
             f"rho = {RHO}", f"E = {E}", f"nu = {NU}",
             "ft = 1e12", "cohesion = 1e12", "frictionDeg = 30", "Gf = 1e6",
             "insertion = adaptive",
             "groupContinuum.block = true", "groupContinuum.sphere = true",
             *v["keys"], "contactMu = 0",
             "fix.bottom = all", "fix.symx = x", "fix.symy = y",
             "trackGroup = sphere", "contactForcePairs = block:sphere"]
    if v["kind"] == "impact":
        lines += [f"groupVel.sphere = 0 0 {-V0}"]
    else:
        lines += [f"velocity.sphere = free free {-VQS}", f"amplitude.sphere = ramp {TRAMP_QS}"]
    return lines


def make_job(var, na, outdir=None, T=None):
    v = VARIANTS[var]
    od = outdir or os.path.join(OUT, tag(var, na))
    ensure_mesh(na)
    T = T if T is not None else t_end(v["kind"])
    cfg = C.write_deck(od, deck_lines(var, na, T),
                       header=f"{BENCH} sphere sur plan, {var}, a/hc = {na} (ecrit par p2_sphere.py)")
    w = na ** 4 * (T / t_end("impact")) * (2.0 if "pf10" in var else 1.0)
    return C.Job(BENCH, tag(var, na), cfg, od, weight=w, meta=dict(var=var, na=na, T=T))


def JOBS(args=None):
    vs = getattr(args, "variants", None) or list(VARIANTS)
    nas = getattr(args, "na", None)
    jobs = []
    for var in vs:
        for na in VARIANTS[var]["na"]:
            if nas and na not in nas:
                continue
            jobs.append(make_job(var, na))
    return jobs


# ---- depouillement -----------------------------------------------------------
def parse_bodies(txt):
    out = {}
    for m in re.finditer(r"corps '(\w+)': KE = (\S+) J, vz moyenne = (\S+) m/s, masse = (\S+) kg", txt):
        out[m.group(1)] = dict(KE=float(m.group(2)), vz=float(m.group(3).rstrip(",")),
                               m=float(m.group(4)))
    m = re.search(r"stocke elastique (\S+) J\)", txt)
    eel = float(m.group(1)) if m else float("nan")
    m = re.search(r"energy budget \(V2/B4\): KE (\S+) -> (\S+) J", txt)
    ke0 = float(m.group(1)) if m else float("nan")
    return out, eel, ke0


def analyse_run(od, var):
    lg = C.parse_log(od)
    h = C.read_history(od)
    kind = VARIANTS[var]["kind"]
    t = h["t"]
    F = 4.0 * h["Fc_block_sphere_z"]               # quart -> sphere entiere, N
    z = h["grpZ"]
    res = dict(ntet=lg["ntet"], dt=lg["dt"], steps=lg["steps"], wall=lg["wall"],
               budget_pct=lg["budget_pct"], done=lg["done"], T=float(t[-1]))
    if kind == "qs":
        d = (z[0] - z) - GAP                        # enfoncement impose
        Fe = 4.0 / 3.0 * ESTAR_QS * math.sqrt(R) * np.maximum(d, 0.0) ** 1.5
        w = (d >= 0.2 * DMAX) & (d <= DMAX)
        res["qs_F_l2"] = float(np.sqrt(np.sum((F[w] - Fe[w]) ** 2) / np.sum(Fe[w] ** 2))) if w.any() else float("nan")
        res["qs_d_reached"] = float(d.max())
        return res, dict(t=t, d=d, F=F, Fe=Fe)
    d = (z[0] - z) - GAP                            # rapprochement des points lointains
    tr = t - GAP / V0                               # temps depuis le premier contact exact
    vz = h["grpVz"]
    Fe = np.interp(tr, _TH, _FH, left=0.0, right=0.0)
    # force : brute (pic, documente) et filtree par moyenne glissante de largeur
    # FWIN t_c (les modes propres de la sphere et du bloc, non amortis, font
    # osciller la force autour de la courbe quasi statique de Hertz) ; la
    # solution exacte passe par le MEME filtre avant comparaison
    Ff, Fef = smooth(tr, F), smooth(tr, Fe)
    res["Fmax_brut"] = float(F.max())
    res["Fmax_brut_err"] = res["Fmax_brut"] / FMAX - 1
    res["Fmax"] = float(Ff.max())
    res["Fmax_err"] = res["Fmax"] / float(Fef.max()) - 1
    res["dmax"] = float(d.max())
    res["dmax_err"] = res["dmax"] / DMAX - 1
    # duree de contact cinematique : delta repasse par zero (centroide de la
    # sphere, insensible aux vibrations) ; exact 2,9432 delta_max / v0
    i_on = np.where(d > 0)[0]
    t_up = float("nan")
    if len(i_on) and i_on[-1] + 1 < len(d):
        k = i_on[-1]
        t_up = float(tr[k] + d[k] * (tr[k + 1] - tr[k]) / (d[k] - d[k + 1]))
    res["tc"] = t_up                                 # depuis le premier contact exact
    res["tc_err"] = res["tc"] / TC - 1
    # duree au seuil de force (filtree), documentee
    t0, t1 = contact_window(tr, Ff, FTHR * FMAX)
    e0, e1 = contact_window(tr, Fef, FTHR * FMAX)
    res["tc_force"] = t1 - t0
    res["tc_force_err"] = res["tc_force"] / (e1 - e0) - 1
    res["impulse"] = float(np.trapezoid(F, t))       # = m (1 + e) v0 si le contact est fini
    sep = bool(np.isfinite(t_up) and tr[-1] > t_up)
    res["separe"] = sep
    res["e"] = float(vz[-1] / V0) if sep else float("nan")
    res["e_impulse"] = res["impulse"] / (MASS * V0) - 1
    res["F_l2"] = float(np.sqrt(np.sum((Ff - Fef) ** 2) / np.sum(Fef ** 2)))
    res["F_l2_brut"] = float(np.sqrt(np.sum((F - Fe) ** 2) / np.sum(Fe ** 2)))
    bodies, eel, ke0 = parse_bodies(lg["text"])
    if "sphere" in bodies and np.isfinite(ke0) and ke0 > 0:
        ms = bodies["sphere"]["m"]
        ke_tr = 0.5 * ms * vz[-1] ** 2
        res["KE0_quart"] = ke0
        res["KE0_ref_quart"] = 0.25 * 0.5 * MASS * V0 ** 2
        res["e_KE"] = float(math.sqrt(ke_tr / ke0)) if sep else float("nan")
        res["part_translation"] = ke_tr / ke0
        res["part_vib_sphere"] = (bodies["sphere"]["KE"] - ke_tr) / ke0
        res["part_KE_bloc"] = bodies.get("block", {}).get("KE", float("nan")) / ke0
        res["part_elastique"] = eel / ke0
        res["closure"] = ((bodies["sphere"]["KE"] + bodies.get("block", {}).get("KE", 0.0) + eel) / ke0 - 1
                          if sep else float("nan"))   # potentiel de contact nul apres separation
    return res, dict(t=tr, d=d, F=F, Ff=Ff, Fe=Fe, vz=vz)


def verdict(var, runs):
    """verdict au maillage le plus fin de la variante"""
    if not runs:
        return {}
    nas = sorted(runs, key=int)
    fin = runs[nas[-1]]
    if VARIANTS[var]["kind"] == "qs":
        return {"qs_F_l2": fin["qs_F_l2"] <= CRIT["qs_F_l2"]}
    v = {
        "Fmax": abs(fin["Fmax_err"]) <= CRIT["Fmax_err"],
        "tc": abs(fin["tc_err"]) <= CRIT["tc_err"],
        "dmax": abs(fin["dmax_err"]) <= CRIT["dmax_err"],
        "e": bool(fin["e"] >= CRIT["e_min"]),
        "budget": all(abs(r["budget_pct"]) <= CRIT["budget_pct"] for r in runs.values()
                      if r["budget_pct"] is not None),
        "closure": abs(fin.get("closure", float("nan"))) <= CRIT["closure"],
    }
    if len(nas) >= 2:
        errs = [abs(runs[n]["Fmax_err"]) for n in nas]
        v["monotone"] = all(b <= a for a, b in zip(errs, errs[1:]))
    return v


def reference():
    return dict(R=R, v0=V0, rho=RHO, E=E, nu=NU, gap=GAP, Estar=ESTAR, Estar_qs=ESTAR_QS,
                m=MASS, K=KH, delta_max=DMAX, F_max=FMAX, t_c=TC, t_c_seuil=TC_THR,
                seuil_force=FTHR, a_max=AMAX, c0=C0, hunter_lambda=HUNTER,
                e_hunter=math.sqrt(1 - HUNTER), criteres=CRIT)


def print_reference():
    print(f"[{BENCH}] Hertz : m = {MASS * 1e3:.4f} g, E* = {ESTAR / 1e9:.4f} GPa, "
          f"delta_max = {DMAX * 1e6:.3f} um, F_max = {FMAX:.3f} N, t_c = {TC * 1e6:.2f} us "
          f"(au seuil {FTHR:g} F_max : {TC_THR * 1e6:.2f} us), a_max = {AMAX * 1e3:.4f} mm, "
          f"a/R = {AMAX / R:.4f}, Hunter lambda = {HUNTER:.4f}")


def print_run(name, r):
    if "Fmax" in r:
        print(f"  {name:18s} ntet {r['ntet']}  dt {r['dt']:.3e}  F_max filtre {r['Fmax']:9.3f} N ({100 * r['Fmax_err']:+6.2f} %)"
              f"  brut {r['Fmax_brut']:9.3f} N ({100 * r['Fmax_brut_err']:+6.2f} %)"
              f"  d_max {r['dmax'] * 1e6:8.3f} um ({100 * r['dmax_err']:+6.2f} %)"
              f"  t_c {r['tc'] * 1e6:8.2f} us ({100 * r['tc_err']:+6.2f} %)  e {r['e']:.4f}"
              f"  F_L2 {100 * r['F_l2']:5.2f} %  B4 {r['budget_pct']} %  bilan {100 * r.get('closure', float('nan')):+.3f} %"
              f"  mur {r['wall']} s")
    else:
        print(f"  {name:18s} ntet {r['ntet']}  qs F_L2 {100 * r['qs_F_l2']:.2f} %  d atteint "
              f"{r['qs_d_reached'] * 1e6:.2f} um  mur {r['wall']} s")


def analyse(variants=None):
    variants = variants or list(VARIANTS)
    print_reference()
    results, sigs = {}, {}
    for var in variants:
        runs, sv = {}, {}
        for na in VARIANTS[var]["na"]:
            od = os.path.join(OUT, tag(var, na))
            if not C.is_done(od):
                continue
            r, s = analyse_run(od, var)
            runs[str(na)] = r
            sv[na] = s
            print_run(tag(var, na), r)
        if not runs:
            continue
        ver = verdict(var, runs)
        jug = var.startswith("pot") or var.startswith("qs")
        results[var] = dict(runs=runs, verdict=ver, juge=jug,
                            passe=all(ver.values()) if jug else None)
        sigs[var] = sv
        print(f"  {var} : {'PASSE' if results[var]['passe'] else ('ECHEC' if jug else 'documente')} {ver}")
    with open(os.path.join(HERE, "resultats.json"), "w") as f:
        json.dump(dict(banc=BENCH, reference=reference(), variantes=results), f, indent=1)
    figures(sigs, results)
    return results


def figures(sigs, results):
    if not sigs:
        return
    plt = C.plot_style()
    th = _TH * 1e6
    # 1. F(t) et F(delta) de l'impact
    imp = [v for v in sigs if VARIANTS[v]["kind"] == "impact"]
    if imp:
        fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.2))
        ax[0].plot(th, _FH, color="k", lw=2.2, alpha=0.35, label="Hertz")
        ax[1].plot(_DH * 1e6, _FH, color="k", lw=2.2, alpha=0.35, label="Hertz")
        k = 0
        for var in imp:
            for na, s in sorted(sigs[var].items()):
                lab = f"{var}, a/h = {na}"
                ax[0].plot(s["t"] * 1e6, s["F"], lw=0.4, color=f"C{k}", alpha=0.35)
                ax[0].plot(s["t"] * 1e6, s["Ff"], lw=0.9, color=f"C{k}", label=lab)
                ax[1].plot(s["d"] * 1e6, s["Ff"], lw=0.9, color=f"C{k}", label=lab)
                k += 1
        ax[0].set_xlabel("temps depuis le premier contact (µs)")
        ax[0].set_ylabel("force de contact (N)")
        ax[1].set_xlabel(r"rapprochement $\delta$ (µm)")
        ax[1].set_ylabel("force de contact filtrée (N)")
        ax[0].set_xlim(-0.1 * TC * 1e6, 1.4 * TC * 1e6)
        ax[1].set_xlim(left=-0.05 * DMAX * 1e6)
        ax[0].legend(fontsize=6, frameon=False)
        fig.suptitle(r"P2.2 impact de Hertz (quart $\times$ 4) : trait large exact, trait fin brut, trait moyen filtré", fontsize=9)
        fig.tight_layout()
        C.savefig(fig, os.path.join(HERE, "fig_impact"))
        plt.close(fig)
        # 2. convergence des ecarts
        fig, ax = plt.subplots(figsize=(4.6, 3.4))
        for var in imp:
            runs = results[var]["runs"]
            nas = sorted(runs, key=int)
            if not nas:
                continue
            x = [int(n) for n in nas]
            ax.loglog(x, [max(abs(runs[n]["Fmax_err"]), 1e-5) for n in nas], "o-", label=f"{var} : F_max")
            ax.loglog(x, [max(abs(runs[n]["tc_err"]), 1e-5) for n in nas], "s--", label=f"{var} : t_c")
        ax.axhline(CRIT["Fmax_err"], color="k", ls=":", lw=0.8, label="critère 5 %")
        ax.set_xlabel(r"$a_{max}/h_c$")
        ax.set_ylabel("écart relatif à Hertz")
        ax.legend(fontsize=6, frameon=False)
        fig.tight_layout()
        C.savefig(fig, os.path.join(HERE, "fig_convergence"))
        plt.close(fig)
    qs = [v for v in sigs if VARIANTS[v]["kind"] == "qs"]
    if qs:
        fig, ax = plt.subplots(figsize=(4.6, 3.4))
        dd = np.linspace(0, DMAX, 200)
        ax.plot(dd * 1e6, 4 / 3 * ESTAR_QS * math.sqrt(R) * dd ** 1.5, color="k", lw=2.2, alpha=0.35,
                label="Hertz, sphère rigide")
        for var in qs:
            for na, s in sorted(sigs[var].items()):
                ax.plot(s["d"] * 1e6, s["F"], lw=0.8, label=f"{var}, a/h = {na}")
        ax.set_xlabel("enfoncement (µm)")
        ax.set_ylabel("force (N)")
        ax.set_xlim(0, 1.05 * DMAX * 1e6)
        ax.legend(fontsize=7, frameon=False)
        fig.tight_layout()
        C.savefig(fig, os.path.join(HERE, "fig_indentation"))
        plt.close(fig)


# ---- essai court (smoke) -----------------------------------------------------
def smoke(exe, tfrac, var="pot_pf1", na=2, rerun=True):
    """run COURT dans out/_smoke (jamais compte comme run de campagne)"""
    od = os.path.join(OUT, "_smoke", tag(var, na))
    T = GAP / V0 + tfrac * TC
    job = make_job(var, na, outdir=od, T=T)
    print_reference()
    rc, wall = C.run_job(job, exe=exe, threads=1, force=rerun)
    print(f"  smoke {job.name} T = {T * 1e6:.1f} us  rc = {rc}  {wall:.1f} s")
    if rc == 0:
        r, _ = analyse_run(od, var)
        print_run(job.name, r)
        json.dump(dict(reference=reference(), run=r, T=T), open(os.path.join(od, "smoke.json"), "w"),
                  indent=1, default=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["prepare", "run", "analyse", "all", "smoke"])
    ap.add_argument("--exe", default=C.EXE_DEFAULT)
    ap.add_argument("--variants", nargs="+", default=None)
    ap.add_argument("--na", nargs="+", type=int, default=None)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--tfrac", type=float, default=0.6, help="smoke : T = gap/v0 + tfrac t_c")
    ap.add_argument("--smoke-var", default="pot_pf1")
    ap.add_argument("--no-rerun", action="store_true", help="smoke : re-analyse sans relancer")
    a = ap.parse_args()
    exe = os.path.abspath(a.exe)
    if a.action == "smoke":
        smoke(exe, a.tfrac, a.smoke_var, (a.na or [2])[0], rerun=not a.no_rerun)
        return
    jobs = JOBS(a)
    if a.action == "prepare":
        print_reference()
        for j in jobs:
            print(f"  {j.name:18s} {j.cfg}")
    if a.action in ("run", "all"):
        for j in sorted(jobs, key=lambda j: j.weight):
            rc, wall = C.run_job(j, exe=exe, threads=a.threads)
            print(f"  {j.name:18s} rc = {rc}  {wall:8.1f} s", flush=True)
    if a.action in ("analyse", "all"):
        analyse(a.variants)


if __name__ == "__main__":
    main()
