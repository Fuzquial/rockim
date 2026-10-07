#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# P3.1 : fissure pressurisee, sensibilite au maillage (Guo 2014, §2.4).
#
# Domaine 120 x 120 mm (xy), epaisseur 20 mm (z), fissure centrale traversante
# de longueur 2a = 40 mm dans le plan y = 60 mm, chargee par une pression
# interne P(t) = 1e10 t (Pa) sur ses deux levres. Demi-modele x >= 0 (plan de
# symetrie yz) avec rouleau fix.sym = x. Materiau de la table 2.3 de Guo :
# rho 2340, E 26 GPa, nu 0,2, ft 3 MPa, c 15 MPa, phi 30 deg, Gf 10 et 50 N/m,
# frottement 0,6. Maillages tetraedriques STRUCTURES (24 tets par cube : un
# point au centre de chaque face et du cube), h = 20, 10, 5, 2,5 (1,25) mm :
# 432, 3 456, 27 648, 221 184 (1 769 472) elements, les nombres de la table 2.4.
#
# Modele de la fissure : noeuds DUPLIQUES sur la partie fissuree du plan
# y = 60 mm (x < a), comme le plugin Crack de Gmsh. rockim ne soude pas les
# noeuds coincidents (buildMeshFile), donc les faces des deux levres sont des
# faces EXTERIEURES distinctes du meme corps, sans joint entre elles : chaque
# levre est un groupe physique de surface (lip_lo, lip_hi) qui recoit
# pressure.<g> (suiveuse, normale sortante de SON tetraedre, donc la pression
# ouvre la fissure). Le front (x = a) n'est pas duplique.
#
# Metriques :
#   P_onset : pression a l'amorcage = 1e10 x (premier instant ou nBroken >= 1
#             dans history.csv) ; recoupe par min(tBreak) des joints du VTU final
#   t_reach : premier joint rompu touchant le bord x = W (VTU final, tBreak)
#   t_prop  : t_reach - t_onset (temps de propagation de Guo, fig. 2.30)
#   t_split : dernier joint rompu du plan y = 60 mm au-dela du front (plan
#             entierement separe), si tout le ligament a rompu
#   souplesse elastique : integrale de l'ouverture de la bouche (x = 0,
#             U_mouth_hi_y - U_mouth_lo_y) sur t in [5, 60] us (P <= 0,6 MPa),
#             rapportee au meme maillage en fem3d (continu exact)
#
#   python3 vv/P3_fissure_guo/p3_fissure_guo.py prepare            # maillages + decks
#   python3 vv/P3_fissure_guo/p3_fissure_guo.py run --threads 2 [--var ...] [--h ...] [--gf ...]
#   python3 vv/P3_fissure_guo/p3_fissure_guo.py analyse
#   python3 vv/run_queue.py P3_fissure_guo --slots 4 --threads 2
# ---------------------------------------------------------------------------
import argparse, json, math, os, re, sys, warnings
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
VV = os.path.dirname(HERE)
ROOT = os.path.dirname(VV)
OUT = os.path.join(HERE, "out")
MESHES = os.path.join(HERE, "meshes")
sys.path.insert(0, VV)
import vvcommon as C  # noqa: E402

BENCH = "P3.1"

# ---- probleme (fixe AVANT le calcul) ---------------------------------------
W, HY, TZ = 0.060, 0.120, 0.020        # demi-modele : x in [0, 60] mm
A_CRACK = 0.020                        # demi-longueur de fissure (fig. 2.17 : 2a = 40 mm)
YC = 0.060                             # plan de la fissure
RHO, E, NU = 2340.0, 26e9, 0.2
FT, COH, PHI = 3e6, 15e6, 30.0
MU = 0.6
RATE = 1e10                            # Pa/s, P = 1e10 t (eq. 2.61)
GFS = [10.0, 50.0]
HS = [20.0, 10.0, 5.0, 2.5]            # mm ; 1,25 mm avec --fine (cout, voir README)
H_FINE = 1.25
T_END = {10.0: 5.0e-4, 50.0: 8.0e-4}   # s : P_max 5 et 8 MPa (> limite de ligament 6 MPa)
T_ELAS = 1.0e-4                        # fem3d : phase elastique seule (P <= 1 MPa)
OPEN_WIN = (5e-6, 6e-5)                # fenetre de la mesure de souplesse (s)

COMMON = [f"rho = {RHO}", f"E = {E}", f"nu = {NU}",
          "scenario = loads", "mesh = file", "historyFlush = true"]
FRACT = [f"ft = {FT}", f"cohesion = {COH}", f"frictionDeg = {PHI}",
         "gfShearFactor = 1", f"contactMu = {MU}"]

# ---- variantes ---------------------------------------------------------------
# intr_pf100 : variante principale. Loi de joint par defaut de rockim
#   (elastique lineaire, adoucissement lineaire), insertion intrinseque,
#   pf = 100 : biais de celerite predit -0,6 % par la loi de P1.1
#   (c_eff/c = (1 + 1,24/pf)^-1/2), marge sur le seuil pf >= 62 (1 %) parce que
#   alpha = 1,24 a ete mesure sur des maillages non structures ou h = hmin est
#   commande par les tets aplatis ; ici tous les tets sont congruents.
# guo_pf5 : la loi de Solidity selon la these (eq. 2.25-2.31) : branche
#   parabolique, z-curve de Munjiza, delta_c = 3 Gf/f depuis zero, longueur
#   d'arete, joint rompu remis au contact. p0 = 10 E (borne haute de l'eq.
#   2.28) donne pf = p0/(2E) = 5 ; Guo ne publie pas la valeur de §2.4.
# adapt : insertion adaptative (aucun joint initial, continuum exact).
# elas_fem3d : continu exact sur le meme maillage, phase elastique seule.
VARIANTS = {
    "intr_pf100": dict(mode="fdem3d", pf=100, keys=[
        "mode = fdem3d", "insertion = intrinsic", "jointPenaltyFactor = 100",
        "jointXi = 0.01"]),
    "guo_pf5": dict(mode="fdem3d", pf=5, keys=[
        "mode = fdem3d", "insertion = intrinsic", "jointPenaltyFactor = 5",
        "jointPenaltyLength = edge", "jointElastic = parabolic",
        "jointSoftening = munjiza", "jointDeltaC = guo", "jointDeath = damage",
        "jointXi = 0"]),
    "adapt": dict(mode="fdem3d", pf=4, keys=[
        "mode = fdem3d", "insertion = adaptive", "jointXi = 0.01"]),
    "elas_fem3d": dict(mode="fem3d", pf=0, elastic=True, keys=[
        "mode = fem3d", "law = elastic"]),
}
MAIN = "intr_pf100"

# ---- criteres d'acceptation (fixes AVANT le premier calcul, 2026-10-04) -----
# C1 convergence : |P(2,5) / P(5) - 1| <= 5 % par Gf (variante principale)
# C2 Guo : |P / P_Guo - 1| <= 15 % pour guo_pf5 a chaque h <= 10 mm
#    (rapporte, sans verdict, pour les autres variantes)
# C3 physique, maillage le plus fin, variante principale :
#    Gf = 10 : P dans [P_LEFM largeur finie ; P_cohesif lineaire plaque infinie]
#    Gf = 50 : P dans [P_LEFM largeur finie ; limite de ligament ft (W - a)/a]
# C4 la fissure atteint le bord x = W avant T pour tous les runs fdem3d
# C5 souplesse elastique de intr_pf100 a 2 % du fem3d sur le meme maillage
# C6 residu du bilan B4 <= 0,1 % de l'echelle
CRIT = dict(conv=0.05, guo=0.15, compliance=0.02, budget_pct=0.1)


def htag(h):
    return f"{h:g}".replace(".", "p")


def gtag(gf):
    return f"{gf:g}"


def tag(var, h, gf=None):
    return f"{var}_h{htag(h)}" + ("" if gf is None else f"_gf{gtag(gf)}")


# ---- maillage structure a fissure dupliquee --------------------------------
def mesh_path(h):
    return os.path.join(MESHES, f"guo_h{htag(h)}.msh")


def ntets(h):
    return int(round(W / (h * 1e-3))) * int(round(HY / (h * 1e-3))) * int(round(TZ / (h * 1e-3))) * 24


def ensure_mesh(h):
    """24 tets par cube (centre du cube, centre de face, deux coins d'une
    arete de la face), noeuds dupliques sur la levre superieure x < a."""
    p = mesh_path(h)
    if os.path.exists(p):
        return p
    os.makedirs(MESHES, exist_ok=True)
    hm = h * 1e-3
    nx, ny, nz = (int(round(L / hm)) for L in (W, HY, TZ))
    jc, ia = int(round(YC / hm)), int(round(A_CRACK / hm))
    assert abs(nx * hm - W) < 1e-12 and abs(jc * hm - YC) < 1e-12 and abs(ia * hm - A_CRACK) < 1e-12
    xyz, nid = [], 0

    def block(shape, fx, fy, fz):
        nonlocal nid
        I, J, K = np.meshgrid(*(np.arange(s) for s in shape), indexing="ij")
        ids = np.arange(nid, nid + I.size).reshape(shape)
        xyz.append(np.stack([fx(I), fy(J), fz(K)], -1).reshape(-1, 3) * hm)
        nid += I.size
        return ids

    G = block((nx + 1, ny + 1, nz + 1), lambda i: i, lambda j: j, lambda k: k)
    FX = block((nx + 1, ny, nz), lambda i: i, lambda j: j + .5, lambda k: k + .5)
    FY = block((nx, ny + 1, nz), lambda i: i + .5, lambda j: j, lambda k: k + .5)
    FZ = block((nx, ny, nz + 1), lambda i: i + .5, lambda j: j + .5, lambda k: k)
    CC = block((nx, ny, nz), lambda i: i + .5, lambda j: j + .5, lambda k: k + .5)
    xyz = np.concatenate(xyz)
    # duplication des noeuds de la partie fissuree du plan y = YC (x < a)
    src = np.concatenate([G[:ia, jc, :].ravel(), FY[:ia, jc, :].ravel()])
    dup = {int(n): nid + q for q, n in enumerate(src)}
    xyz = np.vstack([xyz, xyz[src]])
    nid += len(src)
    mapping = np.arange(nid)
    for k, v in dup.items():
        mapping[k] = v

    I, J, K = (a.ravel() for a in np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij"))
    cen = CC[I, J, K]
    upper = J == jc                                  # cubes juste au-dessus du plan
    g = lambda di, dj, dk: G[I + di, J + dj, K + dk]
    faces = {  # nom -> (centre, 4 coins en cycle)
        "xm": (FX[I, J, K], [g(0, 0, 0), g(0, 1, 0), g(0, 1, 1), g(0, 0, 1)]),
        "xp": (FX[I + 1, J, K], [g(1, 0, 0), g(1, 1, 0), g(1, 1, 1), g(1, 0, 1)]),
        "ym": (FY[I, J, K], [g(0, 0, 0), g(1, 0, 0), g(1, 0, 1), g(0, 0, 1)]),
        "yp": (FY[I, J + 1, K], [g(0, 1, 0), g(1, 1, 0), g(1, 1, 1), g(0, 1, 1)]),
        "zm": (FZ[I, J, K], [g(0, 0, 0), g(1, 0, 0), g(1, 1, 0), g(0, 1, 0)]),
        "zp": (FZ[I, J, K + 1], [g(0, 0, 1), g(1, 0, 1), g(1, 1, 1), g(0, 1, 1)]),
    }
    tets, tris = [], {"sym": [], "lip_lo": [], "lip_hi": []}
    for nm, (fc, q) in faces.items():
        for m in range(4):
            t = np.stack([cen, fc, q[m], q[(m + 1) % 4]], 1)
            t[upper] = mapping[t[upper]]
            tets.append(t)
            tri = t[:, 1:]
            if nm == "xm":
                tris["sym"].append(tri[I == 0])
            if nm == "yp":
                tris["lip_lo"].append(tri[(J == jc - 1) & (I < ia)])
            if nm == "ym":
                tris["lip_hi"].append(tri[(J == jc) & (I < ia)])
    tets = np.concatenate(tets)
    assert len(tets) == ntets(h)
    phys = {"sym": 1, "lip_lo": 2, "lip_hi": 3}
    # bouche de la fissure (x = 0, y = YC) : sommets de chaque levre, groupes de
    # dimension 0 pour la mesure d'ouverture (force.<g> = 0 0 0)
    # (sommet a mi-epaisseur s'il existe, sinon les deux faces z = 0 et z = TZ,
    # symetriques) : un seul sommet evite la ponderation par copies de fdem3d
    km = [nz // 2] if nz % 2 == 0 else list(range(nz + 1))
    mouth = {"mouth_lo": G[0, jc, km], "mouth_hi": mapping[G[0, jc, km]]}
    with open(p, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n6\n")
        for nm, k in phys.items():
            f.write(f'2 {k} "{nm}"\n')
        for q, nm in enumerate(mouth):
            f.write(f'0 {5 + q} "{nm}"\n')
        f.write('3 10 "roche"\n$EndPhysicalNames\n')
        f.write(f"$Nodes\n{len(xyz)}\n")
        f.write("".join(f"{i + 1} {x:.10g} {y:.10g} {z:.10g}\n" for i, (x, y, z) in enumerate(xyz)))
        f.write("$EndNodes\n")
        trs = [(phys[nm], np.concatenate(v)) for nm, v in tris.items()]
        ne = sum(len(t) for _, t in trs) + sum(len(v) for v in mouth.values()) + len(tets)
        f.write(f"$Elements\n{ne}\n")
        eid = 1
        for q, (nm, vv) in enumerate(mouth.items()):
            for v in vv + 1:
                f.write(f"{eid} 15 2 {5 + q} {5 + q} {v}\n")
                eid += 1
        for k, tt in trs:
            for a, b, c in tt + 1:
                f.write(f"{eid} 2 2 {k} {k} {a} {b} {c}\n")
                eid += 1
        lines = []
        for a, b, c, d in tets + 1:
            lines.append(f"{eid} 4 2 10 10 {a} {b} {c} {d}\n")
            eid += 1
        f.write("".join(lines))
        f.write("$EndElements\n")
    print(f"[P3] maillage h = {h:g} mm : {len(tets)} tets, {len(xyz)} noeuds "
          f"({len(dup)} dupliques), {p}")
    return p


# ---- decks et Jobs ------------------------------------------------------------
def deck_lines(var, h, gf, T=None):
    V = VARIANTS[var]
    if T is None:
        T = T_ELAS if V.get("elastic") else T_END[gf]
    lines = list(V["keys"]) + COMMON + [f"meshFile = {mesh_path(h)}", f"T = {T:.6g}",
                                        "frames = 1"]
    if V["mode"] == "fdem3d":
        lines += FRACT + [f"Gf = {gf:g}"]
    # P(t) = RATE t : valeur nominale a T, amplitude lineaire 0 -> 1 sur [0, T]
    lines += ["fix.sym = x", "force.mouth_lo = 0 0 0", "force.mouth_hi = 0 0 0",
              f"pressure.lip_lo = {RATE * T:.6g}", f"amplitude.lip_lo = 0 0 {T:.6g} 1",
              f"pressure.lip_hi = {RATE * T:.6g}", f"amplitude.lip_hi = 0 0 {T:.6g} 1"]
    return lines


def plan(variants=None, hs=None, gfs=None):
    out = []
    for var in variants or [MAIN, "guo_pf5", "adapt", "elas_fem3d"]:
        V = VARIANTS[var]
        for h in hs or HS:
            for gf in ([None] if V.get("elastic") else (gfs or GFS)):
                out.append((var, h, gf))
    return out


def weight(var, h, gf):
    V = VARIANTS[var]
    T = T_ELAS if V.get("elastic") else T_END[gf]
    pf = max(V["pf"], 1)
    k = 0.02 if V["mode"] == "fem3d" else (1.1 if var == "adapt" else 1.0)
    return k * ntets(h) * (T / 5e-4) * (20.0 / h) * math.sqrt(pf / 100.0) / 1e3


def make_job(var, h, gf, outroot=OUT, T=None):
    ensure_mesh(h)
    od = os.path.join(outroot, tag(var, h, gf))
    cfg = C.write_deck(od, deck_lines(var, h, gf, T),
                       header=f"P3.1 fissure pressurisee (Guo 2014 §2.4), {var}, h = {h:g} mm"
                              + ("" if gf is None else f", Gf = {gf:g} N/m")
                              + " (ecrit par p3_fissure_guo.py)")
    return C.Job(BENCH, tag(var, h, gf), cfg, od, weight(var, h, gf),
                 dict(var=var, h=h, gf=gf))


def JOBS(args=None):
    """campagne complete hors 1,25 mm (run_queue.py)"""
    return [make_job(v, h, g) for v, h, g in plan()]


# ---- references physiques ----------------------------------------------------
def tada_F(alpha):
    """facteur de largeur finie d'une fissure centrale (Tada, plaque de
    largeur 2b, alpha = a/b, precision 0,1 %)"""
    return (1 - 0.025 * alpha ** 2 + 0.06 * alpha ** 4) * math.sqrt(1 / math.cos(math.pi * alpha / 2))


def _kernel(x, s1, s2, c):
    from scipy.integrate import quad
    xa = math.sqrt(max(c * c - x * x, 0.0))

    def g(s):
        r = math.sqrt(max(c * c - s * s, 0.0))
        d = r - xa
        return 0.0 if d == 0 else math.log(abs((r + xa) / d))
    pts = [x] if s1 < x < s2 else None
    return quad(g, s1, s2, points=pts, limit=200)[0]


def zcurve(D):
    A, B, Cc = 0.63, 1.8, 6.0
    D = min(max(D, 0.0), 1.0)
    return (1 - (A + B - 1) / (A + B) * math.exp(D * (A + Cc * B) / ((A + B) * (1 - A - B))))\
        * (A * (1 - D) + B * (1 - D) ** Cc)


def cohesive_pc(gf, law="lin", N=40):
    """pression critique d'une fissure cohesive 2a pressurisee dans une plaque
    INFINIE en deformation plane, loi rigide-adoucissante : instant ou
    l'ouverture en x = a atteint delta_c. Ouverture de Sneddon par fonction de
    Green (fissure fictive de demi-longueur c, contraintes cohesives sur
    [a, c], facteur K nul en c), collocation a N panneaux constants.
    law : 'dug' (Dugdale, sigma = ft, delta_c = Gf/ft), 'lin' (lineaire,
    delta_c = 2 Gf/ft), 'z' (z-curve de Munjiza, delta_c = 3 Gf/ft de Guo)."""
    from scipy.optimize import brentq, fsolve
    Ep = E / (1 - NU ** 2)
    k = 4 / (math.pi * Ep)
    a = A_CRACK
    dc = {"dug": gf / FT, "lin": 2 * gf / FT, "z": 3 * gf / FT}[law]
    f = {"lin": lambda d: max(1 - d / dc, 0.0), "z": lambda d: zcurve(d / dc)}.get(law)

    def solve(c):
        s = a + (c - a) * np.linspace(0, 1, N + 1)
        xm = 0.5 * (s[1:] + s[:-1])
        dA = np.arcsin(np.clip(s[1:] / c, 0, 1)) - np.arcsin(np.clip(s[:-1] / c, 0, 1))
        if law == "dug":
            P, sig = FT * dA.sum() / math.asin(a / c), np.full(N, FT)
        else:
            I0 = np.array([_kernel(x, 0, a, c) for x in xm])
            Ij = np.array([[_kernel(x, s[j], s[j + 1], c) for j in range(N)] for x in xm])
            M = np.zeros((N + 1, N + 1)); b = np.zeros(N + 1)
            M[:N, 0] = FT * k * I0 / (2 * gf / FT)
            M[:N, 1:] = -FT * k * Ij / (2 * gf / FT) + np.eye(N)
            b[:N] = FT
            M[N, 0], M[N, 1:] = math.asin(a / c), -dA
            x0 = np.linalg.solve(M, b)
            if law == "z":
                def res(u):
                    P, sg = u[0], u[1:]
                    d = k * (P * I0 - Ij @ sg)
                    return np.concatenate([sg - FT * np.array([f(q) for q in d]),
                                           [P * math.asin(a / c) - dA @ sg]])
                with warnings.catch_warnings():   # loin de la racine fsolve peut stagner
                    warnings.simplefilter("ignore")
                    x0 = fsolve(res, x0, xtol=1e-10)
            P, sig = x0[0], x0[1:]
        da = k * (P * _kernel(a, 0, a, c) - sum(sig[j] * _kernel(a, s[j], s[j + 1], c) for j in range(N)))
        return P, da
    cs = a * np.exp(np.linspace(np.log(1.01), np.log(30), 50))
    prev = cs[0]
    for c in cs:
        if solve(c)[1] > dc:
            break
        prev = c
    c = brentq(lambda c: solve(c)[1] - dc, prev, c, xtol=1e-7)
    return float(solve(c)[0]), float(c)


def references():
    Ep = E / (1 - NU ** 2)
    a, b = A_CRACK, W
    F = tada_F(a / b)
    ref = dict(E_prime=Ep, a=a, b=b, F_tada=F,
               P_ligament=FT * (b - a) / a,
               note="plaque infinie et deformation plane sauf mention ; "
                    "P_ligament = ft (W - a)/a, equilibre du demi-domaine a ligament entierement a ft")
    for gf in GFS:
        K = math.sqrt(Ep * gf)
        r = dict(K_Ic=K, K_Ic_contraintes_planes=math.sqrt(E * gf), l_ch=E * gf / FT ** 2,
                 P_LEFM_inf=K / math.sqrt(math.pi * a), P_LEFM_fw=K / (F * math.sqrt(math.pi * a)))
        for law in ("dug", "lin", "z"):
            P, c = cohesive_pc(gf, law)
            r[f"P_coh_{law}"], r[f"c_coh_{law}"] = P, c
        ref[f"Gf{gtag(gf)}"] = r
    return ref


# ---- depouillement -----------------------------------------------------------
def read_vtu_joints(path):
    """points, triangles et champs de cellule d'un fdem3d_joints_XXXX.vtu ascii"""
    txt = open(path).read()

    def arr(name, dtype=float):
        m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
        return np.array(m.group(1).split(), dtype=dtype) if m else None
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", txt, re.S)
    pts = np.array(m.group(1).split(), float).reshape(-1, 3)
    tri = arr("connectivity", int).reshape(-1, 3)
    return pts, tri, arr("tBreak"), arr("damage")


def last_joint_vtu(od):
    fs = sorted(f for f in os.listdir(od) if f.startswith("fdem3d_joints_") and f.endswith(".vtu"))
    return os.path.join(od, fs[-1]) if fs else None


def analyse_run(od, var, h, gf):
    L = C.parse_log(od)
    r = dict(var=var, h=h, gf=gf, done=L["done"], wall=L["wall"], dt=L["dt"], steps=L["steps"],
             ntet=L["ntet"], budget_pct=L["budget_pct"], broken=L["broken"])
    H = C.read_history(od)
    t = H["t"]
    P = RATE * t
    # souplesse : integrale en temps de l'ouverture de la bouche (x = 0) sur
    # [5, 60] us (P <= 0,6 MPa, avant tout endommagement de la pointe, qui
    # commence vers 0,6 a 1 MPa en intrinseque). Le rapport a fem3d sur la
    # meme fenetre elimine les oscillations dynamiques, communes aux deux.
    dop = H["U_mouth_hi_y"] - H["U_mouth_lo_y"]
    tg = np.linspace(OPEN_WIN[0], OPEN_WIN[1], 400)
    if t[-1] >= OPEN_WIN[1]:
        r["compliance"] = float(np.trapezoid(np.interp(tg, t, dop), tg))   # m s
    if VARIANTS[var].get("elastic"):
        return r
    nb = H["nBroken"]
    i = int(np.argmax(nb >= 1)) if (nb >= 1).any() else None
    if i is not None:
        r["t_onset_hist"] = float(t[i])
        r["t_onset_hist_lo"] = float(t[i - 1]) if i > 0 else 0.0
    vp = last_joint_vtu(od)
    if vp:
        pts, tri, tb, dmg = read_vtu_joints(vp)
        # positions de reference : trame 0 (meme liste de joints jt_, figee a
        # la construction en intrinseque comme en adaptatif) ; la trame finale
        # porte les fragments deplaces apres la separation
        v0 = os.path.join(od, "fdem3d_joints_0000.vtu")
        if os.path.exists(v0) and v0 != vp:
            p0, t0_, _, _ = read_vtu_joints(v0)
            if len(t0_) == len(tri) and np.array_equal(t0_, tri):
                pts = p0
        X = pts[tri]                                   # (n, 3, 3)
        cx, cy = X[:, :, 0].mean(1), X[:, :, 1].mean(1)
        hm = h * 1e-3
        brk = tb >= 0
        band = np.abs(cy - YC) < 0.5 * hm
        plane = np.all(np.abs(X[:, :, 1] - YC) < 0.05 * hm, axis=1) & (cx > A_CRACK)
        edge = band & (X[:, :, 0].max(1) > W - 1e-3 * hm)
        if brk.any():
            j = np.argmin(np.where(brk, tb, np.inf))
            r["t_onset_vtu"] = float(tb[j])
            r["onset_xyz"] = [float(v) for v in X[j].mean(0)]
            r["n_broken_off_band"] = int((brk & ~band).sum())
        if (brk & edge).any():
            r["t_reach"] = float(tb[brk & edge].min())
        if plane.any():
            r["ligament_broken_frac"] = float(brk[plane].mean())
            if brk[plane].all():
                r["t_split"] = float(tb[plane].max())
    to = r.get("t_onset_vtu", r.get("t_onset_hist"))
    if to is not None:
        r["P_onset"] = RATE * to
        if "t_reach" in r:
            r["t_prop"] = r["t_reach"] - to
    return r


def collect(outroot=OUT):
    runs = {}
    for var, h, gf in plan(hs=HS + [H_FINE]):
        od = os.path.join(outroot, tag(var, h, gf))
        if C.is_done(od):
            runs[tag(var, h, gf)] = analyse_run(od, var, h, gf)
    return runs


def verdicts(runs, ref, guo):
    out = {}
    get = lambda v, h, g: runs.get(tag(v, h, g), {})
    for gf in GFS:
        g = gtag(gf)
        # C1 convergence entre les deux maillages les plus fins disponibles
        hs = sorted([h for h in HS + [H_FINE] if "P_onset" in get(MAIN, h, gf)])
        if len(hs) >= 2:
            p1, p2 = get(MAIN, hs[0], gf)["P_onset"], get(MAIN, hs[1], gf)["P_onset"]
            out[f"C1_conv_Gf{g}"] = dict(h=[hs[0], hs[1]], ecart=p1 / p2 - 1,
                                         passe=abs(p1 / p2 - 1) <= CRIT["conv"])
            # C3 physique au plus fin
            P = p1
            lo = ref[f"Gf{g}"]["P_LEFM_fw"]
            hi = ref[f"Gf{g}"]["P_coh_lin"] if gf == 10.0 else ref["P_ligament"]
            out[f"C3_phys_Gf{g}"] = dict(h=hs[0], P=P, borne_basse=lo, borne_haute=hi,
                                         passe=lo <= P <= hi)
        # C2 Guo
        for var in VARIANTS:
            if VARIANTS[var].get("elastic"):
                continue
            for h in HS + [H_FINE]:
                rr = get(var, h, gf)
                if "P_onset" not in rr:
                    continue
                pg = guo["fracture_load_MPa"][f"Gf{g}"][f"{h:g}"]["valeur"] * 1e6
                e = rr["P_onset"] / pg - 1
                d = dict(ecart=e)
                if var == "guo_pf5" and h <= 10:
                    d["passe"] = abs(e) <= CRIT["guo"]
                out[f"C2_guo_{var}_h{htag(h)}_Gf{g}"] = d
    # C4 atteinte du bord, C5 souplesse, C6 bilan
    for k, rr in runs.items():
        if not VARIANTS[rr["var"]].get("elastic"):
            out[f"C4_bord_{k}"] = dict(passe="t_reach" in rr)
        if rr.get("budget_pct") is not None:
            out[f"C6_bilan_{k}"] = dict(val=rr["budget_pct"],
                                       passe=abs(rr["budget_pct"]) <= CRIT["budget_pct"])
    for h in HS + [H_FINE]:
        e = get("elas_fem3d", h, None).get("compliance")
        for var in (MAIN, "guo_pf5", "adapt"):
            for gf in GFS:
                c = get(var, h, gf).get("compliance")
                if e and c:
                    d = dict(ecart=c / e - 1)
                    if var == MAIN:
                        d["passe"] = abs(c / e - 1) <= CRIT["compliance"]
                    out[f"C5_souplesse_{var}_h{htag(h)}_Gf{gtag(gf)}"] = d
    return out


def figures(runs, ref, guo):
    plt = C.plot_style()
    cols = {MAIN: "C0", "guo_pf5": "C3", "adapt": "C2"}
    mk = {10.0: "s", 50.0: "o"}
    for what, ylabel, fn in (("P_onset", "pression d'amorçage (MPa)", "fig_charge_rupture"),
                             ("t_prop", "temps de propagation (µs)", "fig_temps_propagation")):
        fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=(what == "t_prop"))
        for ax, gf in zip(axs, GFS):
            g = gtag(gf)
            key = "fracture_load_MPa" if what == "P_onset" else "propagation_time_us"
            G = guo[key][f"Gf{g}"]
            hh = sorted(float(x) for x in G)
            ax.plot([x / 20 for x in hh], [G[f"{x:g}"]["valeur"] for x in hh], "k" + mk[gf],
                    mfc="none", label="Guo 2014 (fig. 2.29-2.30)")
            for var, col in cols.items():
                pts = sorted((h, runs[tag(var, h, gf)][what]) for h in HS + [H_FINE]
                             if what in runs.get(tag(var, h, gf), {}))
                if pts:
                    sc = 1e-6 if what == "P_onset" else 1e6
                    ax.plot([p[0] / 20 for p in pts], [p[1] * sc for p in pts], mk[gf] + "-",
                            color=col, label=f"rockim {var}")
            if what == "P_onset":
                R = ref[f"Gf{g}"]
                for val, ls, lab in ((R["P_LEFM_fw"], ":", "LEFM, largeur finie"),
                                     (R["P_LEFM_inf"], "-.", "LEFM, plaque infinie"),
                                     (R["P_coh_lin"], "--", "cohésif linéaire, plaque infinie"),
                                     (ref["P_ligament"], "-", "limite du ligament")):
                    if val < 9e6:
                        ax.axhline(val * 1e-6, color="0.5", ls=ls, lw=0.8, label=lab)
            ax.set_xlabel("taille de maille normalisée h/20 mm")
            ax.set_title(f"Gf = {g} N/m", fontsize=9)
            ax.set_xlim(0, 1.1)
        axs[0].set_ylabel(ylabel)
        from matplotlib.ticker import FuncFormatter
        virg = FuncFormatter(lambda v, _: f"{v:g}".replace(".", ","))
        for ax in axs:
            ax.xaxis.set_major_formatter(virg)
            ax.yaxis.set_major_formatter(virg)
            ax.legend(fontsize=6, frameon=False)
        fig.tight_layout()
        C.savefig(fig, os.path.join(HERE, fn))
        plt.close(fig)


def analyse(outroot=OUT, write=True):
    guo = json.load(open(os.path.join(HERE, "guo_reference.json")))
    ref = references()
    runs = collect(outroot)
    ver = verdicts(runs, ref, guo)
    print(f"\n[P3] references : P_ligament = {ref['P_ligament'] / 1e6:.2f} MPa, F_Tada = {ref['F_tada']:.4f}")
    for gf in GFS:
        R = ref[f"Gf{gtag(gf)}"]
        print(f"  Gf = {gf:g} : K_Ic = {R['K_Ic'] / 1e6:.4f} MPa m^1/2, l_ch = {R['l_ch'] * 1e3:.1f} mm, "
              f"P_LEFM inf/fw = {R['P_LEFM_inf'] / 1e6:.3f}/{R['P_LEFM_fw'] / 1e6:.3f} MPa, "
              f"P_coh dug/lin/z = {R['P_coh_dug'] / 1e6:.3f}/{R['P_coh_lin'] / 1e6:.3f}/"
              f"{R['P_coh_z'] / 1e6:.3f} MPa (c = {R['c_coh_lin'] * 1e3:.1f} mm lin)")
    print(f"\n{'run':28s} {'ntet':>8s} {'dt':>9s} {'P_onset':>8s} {'t_prop':>7s} {'lig':>5s} "
          f"{'hors':>5s} {'B4 %':>9s} {'mur s':>8s}")
    for k, r in sorted(runs.items()):
        print(f"{k:28s} {r['ntet'] or 0:8d} {r['dt'] or float('nan'):9.2e} "
              f"{r.get('P_onset', float('nan')) / 1e6:8.3f} {r.get('t_prop', float('nan')) * 1e6:7.1f} "
              f"{r.get('ligament_broken_frac', float('nan')):5.2f} {r.get('n_broken_off_band', 0):5d} "
              f"{r['budget_pct'] if r['budget_pct'] is not None else float('nan'):9.1e} "
              f"{r['wall'] or float('nan'):8.1f}")
    for k, v in ver.items():
        print(f"  {k:44s} {v}")
    if write:
        json.dump(dict(reference=ref, criteres=CRIT, runs=runs, verdicts=ver,
                       passe=all(v.get("passe", True) for v in ver.values())),
                  open(os.path.join(HERE, "resultats.json"), "w"), indent=1)
        figures(runs, ref, guo)
    return runs, ver


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["prepare", "run", "analyse", "all", "smoke"])
    ap.add_argument("--exe", default=C.EXE_DEFAULT)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--var", nargs="+", default=None)
    ap.add_argument("--h", nargs="+", type=float, default=None)
    ap.add_argument("--gf", nargs="+", type=float, default=None)
    ap.add_argument("--fine", action="store_true", help="ajoute h = 1,25 mm")
    ap.add_argument("--T", type=float, default=None, help="duree (smoke)")
    a = ap.parse_args()
    hs = a.h or (HS + ([H_FINE] if a.fine else []))
    if a.action == "smoke":
        # essai de fumee : maillage grossier, duree courte, sorties dans out/smoke/
        root = os.path.join(OUT, "smoke")
        jobs = [make_job(v, h, g, root, a.T) for v, h, g in plan(a.var, a.h or [20.0], a.gf or [10.0])]
        for j in jobs:
            rc, wall = C.run_job(j, os.path.abspath(a.exe), a.threads or 1, force=True)
            print(f"  {j.name:28s} rc = {rc}  {wall:7.1f} s", flush=True)
        analyse(root, write=False)
        return
    jobs = [make_job(v, h, g) for v, h, g in plan(a.var, hs, a.gf)]
    if a.action == "prepare":
        for j in jobs:
            print(f"  {j.name:28s} poids {j.weight:9.1f}  {j.cfg}")
    if a.action in ("run", "all"):
        for j in sorted(jobs, key=lambda j: j.weight):
            rc, wall = C.run_job(j, os.path.abspath(a.exe), a.threads)
            print(f"  {j.name:28s} rc = {rc}  {wall:7.1f} s", flush=True)
    if a.action in ("analyse", "all"):
        analyse()


if __name__ == "__main__":
    main()
