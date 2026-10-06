#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# P4 : essais de laboratoire de Guo 2014 (these, chapitre 3) en fdem3d.
#
#   P4.1  flexion trois points : poutre 320 x 40 x 20 mm entre trois platines
#         d'acier mobiles (Guo §3.3, tableau 3.1, fig. 3.3 et 3.7)
#   P4.2  essai bresilien : disque D = 40 mm, t = 15 mm entre deux platines
#         planes (Guo §3.4, tableau 3.2, fig. 3.9, eq. 3.1 a 3.3)
#
# fdem3d n'a pas de scenario `brazilian` (2D seulement, FdemSolver) ni de
# flexion : les deux essais sont montes en `scenario = loads` avec des CORPS
# nommes du maillage (§5.21) : l'eprouvette (rupture, insertion adaptative ou
# intrinseque) et des platines d'acier continues (groupContinuum.<corps> =
# true) pilotees en bloc par `velocity.<platine> = 0 -Vy 0` et une rampe
# lineaire `amplitude.<platine>` (Guo fig. 3.2). Le seul lien entre corps est
# le contact general par potentiel (contact = potential ; la penalite echoue
# en glissement aux interfaces non conformes, docs/VV_campagne.md).
# Frottement : contactMu.roche = 0,6 (levres de fissure), contactMu.acier =
# 0,1 ; la regle de paire prend le minimum, donc 0,1 entre acier et roche.
#
# Mesures (history.csv) : RF_<platine>_y (reaction de la liaison = force de
# contact transmise), Fc_<platine>_<eprouvette>_y (contactForcePairs, controle),
# U_<point>_y des points de l'axe neutre de la poutre (force nulle).
#
#   python3 vv/P4_essais_guo/p4_essais_guo.py prepare        # maillages + decks
#   python3 vv/P4_essais_guo/p4_essais_guo.py run [--threads 2] [--only motif]
#   python3 vv/P4_essais_guo/p4_essais_guo.py analyse
#   python3 vv/P4_essais_guo/p4_essais_guo.py all
#   python3 vv/P4_essais_guo/p4_essais_guo.py all --smoke    # essai a blanc
#   python3 vv/run_queue.py P4_essais_guo --slots 4 --threads 2
# ---------------------------------------------------------------------------
import argparse, json, math, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
VV = os.path.dirname(HERE)
ROOT = os.path.dirname(VV)
sys.path.insert(0, VV)
import vvcommon as C  # noqa: E402

OUT = os.path.join(HERE, "out")
MESH = os.path.join(HERE, "meshes")

# ===========================================================================
# 1. Donnees de Guo (tableaux 3.1 et 3.2) ; voir guo_reference.json
# ===========================================================================
ROCK41 = dict(rho=2340.0, E=26e9, nu=0.2, ft=3e6, cohesion=10.5e6, frictionDeg=30.0, Gf=30.0)
ROCK42 = dict(rho=2340.0, E=26e9, nu=0.2, ft=3e6, cohesion=15e6, frictionDeg=30.0, Gf=50.0)
STEEL = dict(rho=7850.0, E=200e9, nu=0.28)
MU_CRACK, MU_PLATEN = 0.6, 0.1

# ---- P4.1 geometrie (repere du maillage = repere solveur, boite a l'origine)
LB, HB, BB = 0.320, 0.040, 0.020        # longueur x, hauteur y, epaisseur z
WP, HP = 0.005, 0.005                   # platine : 5 mm en x, 5 mm en y, toute l'epaisseur
Y0 = HP                                 # la poutre repose sur les appuis : y in [HP, HP + HB]
SPAN = LB - WP                          # entre-axe des appuis : 315 mm (Guo ~ 314 mm, lu sur fig. 3.1)
XS = (WP / 2, LB - WP / 2)              # axes des appuis
XC = LB / 2

# ---- P4.2 geometrie
DD, TD = 0.040, 0.015                   # diametre, epaisseur (z)
RD = DD / 2
WPL, HPL = 0.024, 0.005                 # platine 24 x 5 x 15 mm
CX, CY = RD, HPL + RD                   # centre du disque (repere solveur)

# ===========================================================================
# 2. References physiques
# ===========================================================================
def beam_w(xi, L, E, nu, b, h):
    """fleche par unite de force (m/N) d'une poutre sur appuis simples chargee
    au milieu, Euler-Bernoulli + Timoshenko (kappa = 5/6), xi = distance a
    l'axe d'appui (xi < 0 : porte-a-faux, rotation de corps rigide)"""
    I, A = b * h ** 3 / 12, b * h
    G, kap = E / (2 * (1 + nu)), 5.0 / 6.0
    xi = np.asarray(xi, dtype=float)
    xin = np.clip(xi, 0.0, L / 2)
    w = xin * (3 * L ** 2 - 4 * xin ** 2) / (48 * E * I) + xin / (2 * kap * G * A)
    th0 = L ** 2 / (16 * E * I) + 1 / (2 * kap * G * A)
    return np.where(xi < 0, th0 * xi, w)


def beam_ref(m=ROCK41):
    E, nu, ft, Gf = m["E"], m["nu"], m["ft"], m["Gf"]
    I = BB * HB ** 3 / 12
    kb = 48 * E * I / SPAN ** 3
    G = E / (2 * (1 + nu))
    ks = 1 / (1 / kb + SPAN / (4 * (5 / 6) * G * BB * HB))
    lch = E * Gf / ft ** 2
    # premiere frequence propre de flexion (appuis simples, Euler-Bernoulli)
    f1 = (math.pi / 2) / SPAN ** 2 * math.sqrt(E * I / (m["rho"] * BB * HB))
    F_ft = ft * 2 * BB * HB ** 2 / (3 * SPAN)      # charge qui porte la contrainte de flexion a ft
    return dict(I=I, K_EB=kb, K_timo=ks, shear_corr=kb / ks - 1, lch=lch, h_over_lch=HB / lch,
                f1=f1, F_at_ft=F_ft, MOR_band=MOR_BAND)


def mor(F):
    return 3 * F * SPAN / (2 * BB * HB ** 2)


def hondros(rho, F, alpha, R=RD, t=TD):
    """contraintes de Hondros sur le diametre charge (traction > 0) :
    sigma_xx (perpendiculaire au diametre) et sigma_yy, charge F repartie sur
    un arc 2 alpha"""
    p = F / (2 * alpha * R * t)
    r2 = np.asarray(rho, dtype=float) ** 2
    A = (1 - r2) * math.sin(2 * alpha) / (1 - 2 * r2 * math.cos(2 * alpha) + r2 ** 2)
    B = np.arctan((1 + r2) / (1 - r2) * math.tan(alpha))
    return 2 * p / math.pi * (A - B), -2 * p / math.pi * (A + B)


def contact_E(m):
    return 1 / ((1 - m["nu"] ** 2) / m["E"] + (1 - STEEL["nu"] ** 2) / STEEL["E"])


def hertz_b(F, m=ROCK42, R=RD, t=TD):
    """demi-largeur de contact cylindre-plan (Hertz 2D, par unite d'epaisseur)"""
    return math.sqrt(4 * (F / t) * R / (math.pi * contact_E(m)))


def disc_closure(F, m=ROCK42, R=RD, t=TD):
    """rapprochement des platines (estimation, Johnson 1985 §5.6) : sert a
    choisir T, jamais de critere"""
    b = hertz_b(F, m, R, t)
    return 2 * (F / t) / (math.pi * contact_E(m)) * (2 * math.log(4 * R / b) - 1)


def brazil_ref(m=ROCK42):
    F_ft = math.pi * DD * TD * m["ft"] / 2            # charge pour 2F/(pi D t) = ft
    c = math.sqrt(m["E"] / m["rho"])
    return dict(F_at_ft=F_ft, b_at_Fft=hertz_b(F_ft), closure_at_Fft=disc_closure(F_ft),
                t_transit=DD / c, lch=m["E"] * m["Gf"] / m["ft"] ** 2)


# ===========================================================================
# 3. Criteres d'acceptation (fixes le 2026-10-04, AVANT tout calcul complet)
# ===========================================================================
MOR_BAND = (1.3, 2.1)     # MOR/ft attendu a h/l_ch = 0,46 (fissure cohesive, Petersson 1981, Hillerborg 1983)
CRIT = dict(
    # P4.1
    K41_err_fin=0.03,        # raideur elastique F/delta contre Timoshenko, maillage 2,5 mm, adaptatif
    K41_err_gros=0.06,       # idem au maillage 5 mm (rigidite de flexion des tets lineaires)
    MOR_band=MOR_BAND,       # MOR / ft dans la bande de l'effet d'echelle cohesif
    rate41=0.03,             # pics des deux vitesses les plus lentes a 3 % (convergence en vitesse)
    mesh41=0.10,             # pic 5 mm contre 2,5 mm a 10 % (Guo : -7 %)
    guo41=0.15,              # pic contre Guo (299 N a 2,5 mm, 0,005 m/s) a 15 %
    init41_dx=2.0,           # amorcage a moins de 2 h du milieu, en fibre inferieure (y < HP + h)
    # P4.2
    hondros_centre=0.05,     # sigma1 au centre contre Hondros a F = 0,5 F_pic : 5 %
    hondros_L2=0.10,         # profil sigma1(y), |y| <= 0,8 R, ecart L2 relatif : 10 %
    fbt_band=(0.85, 1.15),   # f_bt / ft pour un disque sans defaut et un critere de Rankine
    rate42=0.03,             # deux vitesses les plus lentes a 3 %
    mesh42=0.10,             # pic 2 mm contre 1,2 mm a 10 %
    guo42=0.15,              # f_bt contre Guo 1,96 MPa a 15 % (echec ATTENDU, voir README)
    init42_r=0.5,            # premiere facette rompue a |y - yc| <= 0,5 R (amorcage central)
    budget_pct=1.0,          # residu du bilan B4 < 1 % de l'echelle
)

# ===========================================================================
# 4. Variantes (cout estime en tet-pas, voir README)
# ===========================================================================
# vitesse principale 0,02 m/s : Guo trouve un pic a +1,5 % de la valeur
# convergee (fig. 3.3) pour P4.1 et identique a 0,01 m/s (fig. 3.9) pour P4.2.
V41 = [
    dict(id="ad_h5_v0p05", h=5.0, vy=0.05, ins="adaptive"),
    dict(id="ad_h5_v0p02", h=5.0, vy=0.02, ins="adaptive"),
    dict(id="ad_h5_v0p01", h=5.0, vy=0.01, ins="adaptive"),
    dict(id="ad_h2p5_v0p02", h=2.5, vy=0.02, ins="adaptive"),
    dict(id="in20_h5_v0p02", h=5.0, vy=0.02, ins="intrinsic", pf=20),
]
V42 = [
    dict(id="ad_h1p2_v0p1", h=1.2, vy=0.1, ins="adaptive"),
    dict(id="ad_h1p2_v0p05", h=1.2, vy=0.05, ins="adaptive"),
    dict(id="ad_h1p2_v0p02", h=1.2, vy=0.02, ins="adaptive"),
    dict(id="ad_h2_v0p02", h=2.0, vy=0.02, ins="adaptive"),
    dict(id="in20_h2_v0p02", h=2.0, vy=0.02, ins="intrinsic", pf=20),
]
# essai a blanc : maillages grossiers, vitesse forte, T court (chaine seule)
SMOKE41 = [dict(id="smoke_h10_v1", h=10.0, vy=1.0, ins="adaptive", ramp=1e-5, T=6e-5)]
SMOKE42 = [dict(id="smoke_h4_v1", h=4.0, vy=1.0, ins="adaptive", ramp=5e-6, T=4e-5)]
RAMP41, RAMP42 = 1e-3, 2e-4                 # rampes de Guo (fig. 3.2 ; §3.4.1)
PLATEN_H41 = 5.0                            # maille des platines de flexion (mm)
PLATEN_H42 = 2.0


def T41(v):
    if "T" in v:
        return v["T"]
    d = 2.0 * 1.7 * beam_ref()["F_at_ft"] / beam_ref()["K_timo"]   # 2 x fleche au pic attendu
    if v["ins"] == "intrinsic":
        d *= 1.5
    return RAMP41 / 2 + d / (2 * v["vy"])


def T42(v):
    if "T" in v:
        return v["T"]
    d = 2.5 * disc_closure(brazil_ref()["F_at_ft"])
    if v["ins"] == "intrinsic":
        d *= 1.5
    return RAMP42 / 2 + d / (2 * v["vy"])


# ===========================================================================
# 5. Maillages
# ===========================================================================
def _box24(x0, y0, z0, nx, ny, nz, h):
    """parallelepipede structure : cubes de cote h coupes en 24 tets (sommets,
    centres des faces, centre du cube), le decoupage de Guo (24 x 2 048 =
    49 152 tets a 5 mm, 24 x 16 384 = 393 216 a 2,5 mm)"""
    idx, pts, tets = {}, [], []

    def node(i2, j2, k2):                 # coordonnees en demi-mailles
        key = (i2, j2, k2)
        if key not in idx:
            idx[key] = len(pts)
            pts.append((x0 + i2 * h / 2, y0 + j2 * h / 2, z0 + k2 * h / 2))
        return idx[key]

    faces = [((0, 0, 0), (0, 2, 0), (0, 2, 2), (0, 0, 2)),     # x-
             ((2, 0, 0), (2, 0, 2), (2, 2, 2), (2, 2, 0)),     # x+
             ((0, 0, 0), (0, 0, 2), (2, 0, 2), (2, 0, 0)),     # y-
             ((0, 2, 0), (2, 2, 0), (2, 2, 2), (0, 2, 2)),     # y+
             ((0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0)),     # z-
             ((0, 0, 2), (0, 2, 2), (2, 2, 2), (2, 0, 2))]     # z+
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                o = (2 * i, 2 * j, 2 * k)
                B = node(o[0] + 1, o[1] + 1, o[2] + 1)
                for f in faces:
                    fc = tuple(o[d] + sum(c[d] for c in f) // 4 for d in range(3))
                    Fc = node(*fc)
                    cs = [node(o[0] + c[0], o[1] + c[1], o[2] + c[2]) for c in f]
                    for a in range(4):
                        tets.append([cs[a], cs[(a + 1) % 4], Fc, B])
    P = np.array(pts)
    T = np.array(tets)
    v = np.einsum("ij,ij->i", np.cross(P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]]),
                  P[T[:, 3]] - P[T[:, 0]])
    T[v < 0] = T[v < 0][:, [0, 2, 1, 3]]
    return P, T


def write_msh(path, bodies):
    """bodies = [(nom, P, T)] -> MSH 2.2 ASCII, un volume physique par corps,
    noeuds NON partages entre corps (interfaces de contact)"""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$PhysicalNames\n%d\n" % len(bodies))
        for k, (nm, _, _) in enumerate(bodies):
            f.write(f'3 {k + 1} "{nm}"\n')
        f.write("$EndPhysicalNames\n$Nodes\n%d\n" % sum(len(P) for _, P, _ in bodies))
        off, offs = 0, []
        for _, P, _ in bodies:
            offs.append(off)
            for q, p in enumerate(P):
                f.write(f"{off + q + 1} {p[0]:.9g} {p[1]:.9g} {p[2]:.9g}\n")
            off += len(P)
        f.write("$EndNodes\n$Elements\n%d\n" % sum(len(T) for _, _, T in bodies))
        e = 0
        for k, (_, _, T) in enumerate(bodies):
            for t in T:
                e += 1
                n = " ".join(str(offs[k] + int(q) + 1) for q in t)
                f.write(f"{e} 4 2 {k + 1} {k + 1} {n}\n")
        f.write("$EndElements\n")


def mesh41(h_mm):
    p = os.path.join(MESH, f"flexion_h{tag(h_mm)}.msh")
    if os.path.exists(p):
        return p
    h, hp = h_mm * 1e-3, PLATEN_H41 * 1e-3
    nx, ny, nz = (round(LB / h), round(HB / h), round(BB / h))
    beam = _box24(0, Y0, 0, nx, ny, nz, h)
    npz = round(BB / hp)
    pg = _box24(0, 0, 0, 1, 1, npz, hp)
    pd = _box24(LB - WP, 0, 0, 1, 1, npz, hp)
    pc = _box24(XC - WP / 2, Y0 + HB, 0, 1, 1, npz, hp)
    write_msh(p, [("poutre",) + beam, ("appui_g",) + pg, ("appui_d",) + pd, ("poincon",) + pc])
    return p


def mesh42(h_mm):
    p = os.path.join(MESH, f"bresilien_h{tag(h_mm)}.msh")
    if os.path.exists(p):
        return p
    import gmsh
    os.makedirs(MESH, exist_ok=True)
    h, hp = h_mm * 1e-3, PLATEN_H42 * 1e-3
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("bresilien")
    o = gmsh.model.occ
    # disque par extrusion d'un cercle en quatre arcs : sommets exacts aux
    # generatrices haute et basse (contact initial sans jeu avec les platines)
    c = o.addPoint(CX, CY, 0)
    q = [o.addPoint(CX + RD * math.cos(a), CY + RD * math.sin(a), 0)
         for a in (0, math.pi / 2, math.pi, 3 * math.pi / 2)]
    arcs = [o.addCircleArc(q[k], c, q[(k + 1) % 4]) for k in range(4)]
    s = o.addPlaneSurface([o.addCurveLoop(arcs)])
    ext = o.extrude([(2, s)], 0, 0, TD)
    disc = [e[1] for e in ext if e[0] == 3][0]
    pb = o.addBox(CX - WPL / 2, 0, 0, WPL, HPL, TD)
    ph = o.addBox(CX - WPL / 2, HPL + DD, 0, WPL, HPL, TD)
    o.remove([(0, c)])
    o.synchronize()                       # pas de fragment : corps separes
    for vol, nm in ((disc, "disque"), (pb, "plateau_b"), (ph, "plateau_h")):
        gmsh.model.addPhysicalGroup(3, [vol], name=nm)
        pts = gmsh.model.getBoundary([(3, vol)], recursive=True)
        gmsh.model.mesh.setSize(pts, h if nm == "disque" else hp)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.option.setNumber("Mesh.RandomSeed", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.model.mesh.generate(3)
    gmsh.write(p)
    gmsh.finalize()
    return p


# ===========================================================================
# 6. Decks et Job
# ===========================================================================
def tag(x):
    return f"{x:g}".replace(".", "p")


def material_lines(m):
    return [f"rho = {m['rho']}", f"E = {m['E']:g}", f"nu = {m['nu']}", f"ft = {m['ft']:g}",
            f"cohesion = {m['cohesion']:g}", f"frictionDeg = {m['frictionDeg']}", f"Gf = {m['Gf']}",
            "phases = roche acier", "phase.roche.fraction = 0.95", "phase.acier.fraction = 0.05",
            f"phase.acier.rho = {STEEL['rho']}", f"phase.acier.E = {STEEL['E']:g}",
            f"phase.acier.nu = {STEEL['nu']}", "phase.acier.ft = 1e12", "phase.acier.cohesion = 1e12",
            "contact = potential", f"contactMu = {MU_CRACK}", f"contactMu.roche = {MU_CRACK}",
            f"contactMu.acier = {MU_PLATEN}"]


def insertion_lines(v):
    if v["ins"] == "intrinsic":
        return ["insertion = intrinsic", f"jointPenaltyFactor = {v['pf']}"]
    return ["insertion = adaptive"]


def deck41(v, od):
    ramp = v.get("ramp", RAMP41)
    amp = f"0 0 {ramp:g} 1"
    yN = Y0 + HB / 2
    lines = ["mode = fdem3d", "scenario = loads", "mesh = file", f"meshFile = {mesh41(v['h'])}",
             f"T = {T41(v):.6g}", f"frames = {v.get('frames', 4)}", "historyFlush = true",
             *material_lines(ROCK41), *insertion_lines(v),
             "groupPhase.poutre = roche", "groupPhase.appui_g = acier",
             "groupPhase.appui_d = acier", "groupPhase.poincon = acier",
             "groupContinuum.appui_g = true", "groupContinuum.appui_d = true",
             "groupContinuum.poincon = true",
             f"velocity.poincon = 0 {-v['vy']} 0", f"amplitude.poincon = {amp}",
             f"velocity.appui_g = 0 {v['vy']} 0", f"amplitude.appui_g = {amp}",
             f"velocity.appui_d = 0 {v['vy']} 0", f"amplitude.appui_d = {amp}",
             f"point.mil = {XC} {yN} {BB / 2}", "force.mil = 0 0 0",
             f"point.sg = {XS[0]} {yN} {BB / 2}", "force.sg = 0 0 0",
             f"point.sd = {XS[1]} {yN} {BB / 2}", "force.sd = 0 0 0",
             f"point.fib = {XC} {Y0} {BB / 2}", "force.fib = 0 0 0",
             "contactForcePairs = poincon:poutre appui_g:poutre appui_d:poutre",
             "writeRuptureFields = true"]
    return C.write_deck(od, lines, header=f"P4.1 flexion trois points, {v['id']} (ecrit par p4_essais_guo.py)")


def deck42(v, od):
    ramp = v.get("ramp", RAMP42)
    amp = f"0 0 {ramp:g} 1"
    lines = ["mode = fdem3d", "scenario = loads", "mesh = file", f"meshFile = {mesh42(v['h'])}",
             f"T = {T42(v):.6g}", f"frames = {v.get('frames', 10)}", "historyFlush = true",
             *material_lines(ROCK42), *insertion_lines(v),
             "groupPhase.disque = roche", "groupPhase.plateau_h = acier", "groupPhase.plateau_b = acier",
             "groupContinuum.plateau_h = true", "groupContinuum.plateau_b = true",
             f"velocity.plateau_h = 0 {-v['vy']} 0", f"amplitude.plateau_h = {amp}",
             f"velocity.plateau_b = 0 {v['vy']} 0", f"amplitude.plateau_b = {amp}",
             "contactForcePairs = plateau_h:disque plateau_b:disque",
             "writeRuptureFields = true"]
    return C.write_deck(od, lines, header=f"P4.2 bresilien, {v['id']} (ecrit par p4_essais_guo.py)")


def _weight41(v):
    nt = 24 * round(LB / (v["h"] * 1e-3)) * round(HB / (v["h"] * 1e-3)) * round(BB / (v["h"] * 1e-3))
    return nt * T41(v) / (v["h"] * 1e-3) * (2.2 if v["ins"] == "intrinsic" else 1.0)


def _weight42(v):
    nt = 6.0 * math.pi * RD ** 2 * TD / (v["h"] * 1e-3) ** 3 * 1.3
    return nt * T42(v) / (v["h"] * 1e-3) * (2.2 if v["ins"] == "intrinsic" else 1.0)


def JOBS(args=None):
    smoke = bool(getattr(args, "smoke", False)) if args is not None else False
    only = getattr(args, "only", None) if args is not None else None
    l41, l42 = (SMOKE41, SMOKE42) if smoke else (V41, V42)
    jobs = []
    for v in l41:
        od = os.path.join(OUT, f"p41_{v['id']}")
        if only and only not in od:
            continue
        jobs.append(C.Job("P4.1", f"p41_{v['id']}", deck41(v, od), od, _weight41(v), dict(v)))
    for v in l42:
        od = os.path.join(OUT, f"p42_{v['id']}")
        if only and only not in od:
            continue
        jobs.append(C.Job("P4.2", f"p42_{v['id']}", deck42(v, od), od, _weight42(v), dict(v)))
    return jobs


# ===========================================================================
# 7. Depouillement
# ===========================================================================
def _vtu_arrays(path, names):
    """lecteur minimal des VTU ASCII de rockim"""
    txt = open(path).read()
    out = {}
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", txt, re.S)
    out["points"] = np.array(m.group(1).split(), float).reshape(-1, 3)
    m = re.search(r'Name="connectivity"[^>]*>(.*?)</DataArray>', txt, re.S)
    out["conn"] = np.array(m.group(1).split(), int)
    m = re.search(r'Name="offsets"[^>]*>(.*?)</DataArray>', txt, re.S)
    out["offsets"] = np.array(m.group(1).split(), int)
    for nm in names:
        m = re.search(rf'<CellData.*?Name="{nm}"[^>]*>(.*?)</DataArray>', txt, re.S)
        if m:
            out[nm] = np.array(m.group(1).split(), float)
    return out


def _centroids(a, nper):
    c = a["conn"].reshape(-1, nper)
    return a["points"][c].mean(axis=1)


def _translation(log):
    m = re.search(r"maillage translate de \((\S+), (\S+), (\S+)\) m", log)
    return np.array([float(m.group(k)) for k in (1, 2, 3)]) if m else np.zeros(3)


def _point_pos(log, g):
    m = re.search(rf"point\.{g} : sommet \d+ a \S+ m du point demande \(\s*(\S+)\s+(\S+)\s+(\S+),", log)
    return np.array([float(m.group(k)) for k in (1, 2, 3)]) if m else None


def _first_breaks(od, n=5):
    fs = sorted(f for f in os.listdir(od) if f.startswith("fdem3d_joints_") and f.endswith(".vtu"))
    if not fs:
        return None
    a = _vtu_arrays(os.path.join(od, fs[-1]), ["tBreak"])
    if "tBreak" not in a:
        return None
    tb = a["tBreak"]
    ok = np.where(tb >= 0)[0]
    if len(ok) == 0:
        return dict(n_broken=0)
    cen = _centroids(a, 3)
    sel = ok[np.argsort(tb[ok])[:n]]
    return dict(n_broken=int(len(ok)), t_first=float(tb[sel[0]]), xyz_first=cen[sel].mean(axis=0).tolist())


def _peak(t, F):
    i = int(np.argmax(F))
    post = F[i:]
    return i, float(F[i]), bool(len(post) > 3 and post.min() < 0.7 * F[i])


def analyse41(od, v):
    r = C.parse_log(od)
    log = r.pop("text")
    h = C.read_history(od)
    t = h["t"]
    # force de contact du poincon sur la poutre (la reaction RF de la liaison
    # contient en plus l'inertie de la platine pendant la rampe)
    F = -h["Fc_poincon_poutre_y"]
    Fs = h["Fc_appui_g_poutre_y"] + h["Fc_appui_d_poutre_y"]
    Fc = -h["RF_poincon_y"]
    d = -(h["U_mil_y"] - 0.5 * (h["U_sg_y"] + h["U_sd_y"]))
    i, Fmax, dropped = _peak(t, F)
    w = (F[:i + 1] >= 0.1 * Fmax) & (F[:i + 1] <= 0.4 * Fmax)
    K = float(np.polyfit(d[:i + 1][w], F[:i + 1][w], 1)[0]) if w.sum() >= 3 else float("nan")
    # reference evaluee aux positions reelles des sommets retenus
    pos = {g: _point_pos(log, g) for g in ("mil", "sg", "sd")}
    tr = _translation(log)
    if all(p is not None for p in pos.values()):
        xs = [pos["sg"][0] - tr[0] - XS[0], XS[1] - (pos["sd"][0] - tr[0])]
        xm = abs(pos["mil"][0] - tr[0] - XC)
        wm = beam_w(SPAN / 2 - xm, SPAN, ROCK41["E"], ROCK41["nu"], BB, HB)
        ws = 0.5 * (beam_w(xs[0], SPAN, ROCK41["E"], ROCK41["nu"], BB, HB)
                    + beam_w(xs[1], SPAN, ROCK41["E"], ROCK41["nu"], BB, HB))
        K_th = float(1 / (wm - ws))
    else:
        K_th = beam_ref()["K_timo"]
    br = beam_ref()
    hm = v["h"] * 1e-3
    fb = _first_breaks(od)
    init_ok = None
    if fb and fb.get("n_broken", 0) > 0:
        x0, y0 = fb["xyz_first"][0] - tr[0], fb["xyz_first"][1] - tr[1]
        init_ok = bool(abs(x0 - XC) <= CRIT["init41_dx"] * hm and y0 <= Y0 + hm)
    res = dict(id=v["id"], h_mm=v["h"], vy=v["vy"], insertion=v["ins"], ntet=r["ntet"], dt=r["dt"],
               steps=r["steps"], wall=r["wall"], budget_pct=r["budget_pct"], broken=r["broken"],
               F_peak=Fmax, t_peak=float(t[i]), delta_peak=float(d[i]), peak_reached=dropped,
               K_sim=K, K_th=K_th, K_err=K / K_th - 1, MOR=mor(Fmax), MOR_over_ft=mor(Fmax) / ROCK41["ft"],
               F_support_over_F=float(Fs[i] / Fmax) if Fmax else None,
               RF_over_Fc=float(Fc[i] / Fmax) if Fmax else None,
               t_peak_over_T1=float(t[i] * br["f1"]), first_breaks=fb, init_bottom_mid=init_ok)
    return res, dict(t=t, F=F, d=d)


def analyse42(od, v):
    r = C.parse_log(od)
    log = r.pop("text")
    h = C.read_history(od)
    t = h["t"]
    F = 0.5 * (-h["Fc_plateau_h_disque_y"] + h["Fc_plateau_b_disque_y"])   # Guo eq. 3.1
    eps = ((h["U_plateau_b_y"] - h["U_plateau_h_y"])) / DD           # Guo eq. 3.2
    i, Fmax, dropped = _peak(t, F)
    fbt = 2 * Fmax / (math.pi * DD * TD)
    tr = _translation(log)
    # Hondros a la frame la plus proche de F = 0,5 F_pic avant le pic
    hond = None
    fr = os.path.join(od, "frames.csv")
    if os.path.exists(fr) and Fmax > 0:
        fd = np.genfromtxt(fr, delimiter=",", names=True)
        ft_ = np.atleast_1d(fd["t"])
        Ffr = np.interp(ft_, t, F)
        cand = [k for k in range(len(ft_)) if ft_[k] <= t[i] and Ffr[k] > 0.2 * Fmax]
        if cand:
            k = min(cand, key=lambda q: abs(Ffr[q] - 0.5 * Fmax))
            a = _vtu_arrays(os.path.join(od, f"fdem3d_{int(np.atleast_1d(fd['frame'])[k]):04d}.vtu"),
                            ["sigma1", "phase"])
            cen = _centroids(a, 4) - tr
            hm = v["h"] * 1e-3
            disc = a["phase"] == 0 if "phase" in a else np.ones(len(cen), bool)
            band = disc & (np.abs(cen[:, 0] - CX) < 1.0 * hm)
            yy = cen[band, 1] - CY
            s1 = a["sigma1"][band]
            Fk = float(Ffr[k])
            alpha = hertz_b(Fk) / RD
            bins = np.linspace(-RD, RD, 17)
            yc, sm = [], []
            for b0, b1 in zip(bins[:-1], bins[1:]):
                q = (yy >= b0) & (yy < b1)
                if q.sum() >= 3:
                    yc.append(0.5 * (b0 + b1)); sm.append(float(np.mean(s1[q])))
            yc, sm = np.array(yc), np.array(sm)
            sx = hondros(yc / RD, Fk, alpha)[0]
            inner = np.abs(yc) <= 0.8 * RD
            q0 = np.abs(yy) <= 0.15 * RD
            s_c = float(np.mean(s1[q0])) if q0.sum() else float("nan")
            s_c_th = float(hondros(np.array([0.0]), Fk, alpha)[0][0])
            hond = dict(frame=int(np.atleast_1d(fd["frame"])[k]), t=float(ft_[k]), F=Fk, alpha=alpha,
                        sigma_centre=s_c, sigma_centre_th=s_c_th, err_centre=s_c / s_c_th - 1,
                        L2=float(np.sqrt(np.sum((sm[inner] - sx[inner]) ** 2) / np.sum(sx[inner] ** 2))),
                        y=yc.tolist(), s1=sm.tolist(), s_th=sx.tolist(), n_el=int(band.sum()))
    fb = _first_breaks(od)
    init_ok = None
    if fb and fb.get("n_broken", 0) > 0:
        init_ok = bool(abs(fb["xyz_first"][1] - tr[1] - CY) <= CRIT["init42_r"] * RD)
    res = dict(id=v["id"], h_mm=v["h"], vy=v["vy"], insertion=v["ins"], ntet=r["ntet"], dt=r["dt"],
               steps=r["steps"], wall=r["wall"], budget_pct=r["budget_pct"], broken=r["broken"],
               F_peak=Fmax, t_peak=float(t[i]), eps_peak=float(eps[i]), peak_reached=dropped,
               f_bt=fbt, fbt_over_ft=fbt / ROCK42["ft"], fbt_over_guo=fbt / GUO["P42"]["f_bt"] - 1,
               t_peak_over_transit=float(t[i] / brazil_ref()["t_transit"]),
               hondros=hond, first_breaks=fb, init_centre=init_ok)
    return res, dict(t=t, F=F, eps=eps)


GUO = json.load(open(os.path.join(HERE, "guo_reference.json")))["valeurs_cles"]


def verdicts(R41, R42):
    v = {}
    ad41 = {k: r for k, r in R41.items() if r["insertion"] == "adaptive" and not k.startswith("smoke")}
    fin = [r for r in ad41.values() if r["h_mm"] == 2.5]
    gros = [r for r in ad41.values() if r["h_mm"] == 5.0]
    if fin:
        v["P4.1 raideur 2,5 mm"] = abs(fin[0]["K_err"]) <= CRIT["K41_err_fin"]
    if gros:
        v["P4.1 raideur 5 mm"] = all(abs(r["K_err"]) <= CRIT["K41_err_gros"] for r in gros)
    allr = list(ad41.values())
    if allr:
        v["P4.1 MOR/ft"] = all(MOR_BAND[0] <= r["MOR_over_ft"] <= MOR_BAND[1] for r in allr)
        v["P4.1 amorcage"] = all(r["init_bottom_mid"] for r in allr)
    sl = sorted([r for r in gros], key=lambda r: r["vy"])
    if len(sl) >= 2:
        v["P4.1 vitesse"] = abs(sl[1]["F_peak"] / sl[0]["F_peak"] - 1) <= CRIT["rate41"]
    same = [r for r in gros if r["vy"] == 0.02]
    if fin and same:
        v["P4.1 maillage"] = abs(same[0]["F_peak"] / fin[0]["F_peak"] - 1) <= CRIT["mesh41"]
    ref = fin[0] if fin else (same[0] if same else None)
    if ref:
        v["P4.1 Guo (B1)"] = abs(ref["F_peak"] / GUO["P41"]["F_peak_h2p5"] - 1) <= CRIT["guo41"]
    ad42 = {k: r for k, r in R42.items() if r["insertion"] == "adaptive" and not k.startswith("smoke")}
    f42 = sorted([r for r in ad42.values() if r["h_mm"] == 1.2], key=lambda r: r["vy"])
    for r in ad42.values():
        if r["hondros"]:
            v.setdefault("P4.2 Hondros centre", True)
            v["P4.2 Hondros centre"] &= abs(r["hondros"]["err_centre"]) <= CRIT["hondros_centre"]
            v.setdefault("P4.2 Hondros profil", True)
            v["P4.2 Hondros profil"] &= r["hondros"]["L2"] <= CRIT["hondros_L2"]
    if ad42:
        b = CRIT["fbt_band"]
        v["P4.2 f_bt/ft"] = all(b[0] <= r["fbt_over_ft"] <= b[1] for r in ad42.values())
        v["P4.2 amorcage central"] = all(r["init_centre"] for r in ad42.values())
    if len(f42) >= 2:
        v["P4.2 vitesse"] = abs(f42[1]["F_peak"] / f42[0]["F_peak"] - 1) <= CRIT["rate42"]
    g42 = [r for r in ad42.values() if r["h_mm"] == 2.0]
    if f42 and g42:
        ref = [r for r in f42 if r["vy"] == g42[0]["vy"]]
        if ref:
            v["P4.2 maillage"] = abs(g42[0]["F_peak"] / ref[0]["F_peak"] - 1) <= CRIT["mesh42"]
    if f42:
        v["P4.2 Guo (B1)"] = abs(f42[0]["fbt_over_guo"]) <= CRIT["guo42"]
    return v


def figures(S41, S42, R41, R42):
    plt = C.plot_style()
    if S41:
        fig, ax = plt.subplots(figsize=(6.0, 3.8))
        for k, (name, s) in enumerate(sorted(S41.items())):
            ax.plot(s["d"] * 1e3, s["F"], lw=0.9, color=f"C{k}", label=name.replace(".", ","))
        dd = np.linspace(0, 0.1, 50)
        K = beam_ref()["K_timo"]
        ax.plot(dd, K * dd * 1e-3, "k--", lw=0.8, label="Timoshenko, élastique")
        ax.axhline(GUO["P41"]["F_peak_h2p5"], color="0.5", ls=":", lw=0.8, label="pic de Guo (2,5 mm)")
        ax.axhline(beam_ref()["F_at_ft"], color="0.3", ls="-.", lw=0.6, label="MOR = ft")
        ax.set_xlabel("flèche relative de l'axe neutre (mm)")
        ax.set_ylabel("force du poinçon F (N)")
        ax.legend(fontsize=7, frameon=False)
        fig.tight_layout()
        C.savefig(fig, os.path.join(HERE, "fig_flexion"))
        plt.close(fig)
    if S42:
        fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.4))
        for k, (name, s) in enumerate(sorted(S42.items())):
            ax[0].plot(s["eps"] * 100, s["F"], lw=0.9, color=f"C{k}", label=name.replace(".", ","))
        ax[0].axhline(GUO["P42"]["F_peak_v0p01"], color="0.5", ls=":", lw=0.8, label="pic de Guo (0,01 m/s)")
        ax[0].axhline(brazil_ref()["F_at_ft"], color="0.3", ls="-.", lw=0.6, label="2F/(πDt) = ft")
        ax[0].set_xlabel(r"déformation verticale $\varepsilon_{yy}$ (%)")
        ax[0].set_ylabel("force F (N)")
        ax[0].legend(fontsize=6, frameon=False)
        for k, (name, r) in enumerate(sorted(R42.items())):
            hd = r.get("hondros")
            if not hd:
                continue
            y = np.array(hd["y"]) / RD
            ax[1].plot(y, np.array(hd["s1"]) / 1e6, "o", ms=3, color=f"C{k}", label=name.replace(".", ","))
            yy = np.linspace(-0.95, 0.95, 200)
            ax[1].plot(yy, hondros(yy, hd["F"], hd["alpha"])[0] / 1e6, "-", lw=0.8, color=f"C{k}")
        ax[1].set_xlabel("y / R sur le diamètre chargé")
        ax[1].set_ylabel(r"$\sigma_1$ (MPa) ; trait : Hondros")
        ax[1].legend(fontsize=6, frameon=False)
        fig.tight_layout()
        C.savefig(fig, os.path.join(HERE, "fig_bresilien"))
        plt.close(fig)


def analyse(smoke=False):
    R41, R42, S41, S42 = {}, {}, {}, {}
    for v in (SMOKE41 if smoke else V41):
        od = os.path.join(OUT, f"p41_{v['id']}")
        if C.is_done(od):
            R41[v["id"]], S41[v["id"]] = analyse41(od, v)
    for v in (SMOKE42 if smoke else V42):
        od = os.path.join(OUT, f"p42_{v['id']}")
        if C.is_done(od):
            R42[v["id"]], S42[v["id"]] = analyse42(od, v)
    ver = verdicts(R41, R42)
    for k, r in R41.items():
        print(f"P4.1 {k:16s} ntet {r['ntet']}  dt {r['dt']:.3g}  pas {r['steps']}  {r['wall']:.0f} s  "
              f"F_pic {r['F_peak']:.1f} N  K {r['K_sim']:.4g} (th {r['K_th']:.4g}, {100 * r['K_err']:+.2f} %)  "
              f"MOR/ft {r['MOR_over_ft']:.3f}  pic atteint {r['peak_reached']}  B4 {r['budget_pct']}")
    for k, r in R42.items():
        hd = r["hondros"] or {}
        print(f"P4.2 {k:16s} ntet {r['ntet']}  dt {r['dt']:.3g}  pas {r['steps']}  {r['wall']:.0f} s  "
              f"F_pic {r['F_peak']:.1f} N  f_bt {r['f_bt'] / 1e6:.3f} MPa (/ft {r['fbt_over_ft']:.3f})  "
              f"Hondros centre {hd.get('err_centre', float('nan')):+.3f} L2 {hd.get('L2', float('nan')):.3f}  "
              f"pic atteint {r['peak_reached']}")
    for k, ok in ver.items():
        print(f"  {k:28s} {'PASSE' if ok else 'ECHEC'}")
    out = dict(reference=dict(P41=beam_ref(), P42=brazil_ref(), criteres=CRIT, guo=GUO),
               P41=R41, P42=R42, verdicts=ver, smoke=smoke)
    fn = "resultats_smoke.json" if smoke else "resultats.json"
    json.dump(out, open(os.path.join(HERE, fn), "w"), indent=1, default=float)
    figures(S41, S42, R41, R42)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["prepare", "run", "analyse", "all", "ref"])
    ap.add_argument("--exe", default=C.EXE_DEFAULT)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    if a.action == "ref":
        print(json.dumps(dict(P41=beam_ref(), P42=brazil_ref(),
                              T41={v["id"]: T41(v) for v in V41}, T42={v["id"]: T42(v) for v in V42}),
                         indent=1, default=float))
        return
    jobs = JOBS(a)
    if a.action == "prepare":
        for j in jobs:
            print(f"  {j.name:24s} {j.cfg}  poids {j.weight:.3g}")
    if a.action in ("run", "all"):
        for j in sorted(jobs, key=lambda j: -j.weight):
            rc, w = C.run_job(j, os.path.abspath(a.exe), a.threads)
            print(f"  {j.name:24s} rc = {rc}  {w:7.1f} s", flush=True)
    if a.action in ("analyse", "all"):
        analyse(a.smoke)


if __name__ == "__main__":
    main()
