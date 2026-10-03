#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# P1.2 + P1.3 : elasticite lineaire de l'element en 3D, quasi statique.
#
# P1.3, patch test (cube biaxial, Lei et Rougier 2016, cas 1, reproduit en
#   forme) : cube L x L x L, traction morte -PX sur xmax, -PY sur ymax,
#   rouleaux sur xmin (x), ymin (y), bottom (z), top libre. Champ exact
#   homogene (tets lineaires : represente exactement) :
#     sxx = -PX, syy = -PY, szz = 0
#     exx = (sxx - nu syy)/E, eyy = (syy - nu sxx)/E, ezz = -nu (sxx + syy)/E
#     u = (exx x, eyy y, ezz z)   (origine au coin c000, tenu par les rouleaux)
#   Metrique : erreur nodale max |u - u_ex| / max |u_ex| sur TOUS les noeuds,
#   lue dans les VTU (coordonnees a 12 chiffres : resolution ~5e-13 m), plus
#   les deplacements moyens de groupes U_xmax_x, U_ymax_y, U_top_z de
#   history.csv (6 chiffres seulement).
#
# P1.2, Kirsch en 3D : quart de plaque [0,W]x[0,W]x[0,t] percee d'un quart de
#   trou de rayon a centre a l'origine (symetries : rouleaux sur xmin et ymin),
#   deformation plane par rouleaux en z sur bottom et top, traction morte
#   SX sur xmax et LAMBDA*SX sur ymax (compression : SX < 0). Solution de
#   Kirsch (Jaeger et Cook, §10.4), en polaire pour une traction uniaxiale s
#   selon x (theta mesure depuis x) :
#     srr = s/2 (1 - a2/r2) + s/2 (1 - 4 a2/r2 + 3 a4/r4) cos 2t
#     stt = s/2 (1 + a2/r2) - s/2 (1 + 3 a4/r4) cos 2t
#     srt = -s/2 (1 + 2 a2/r2 - 3 a4/r4) sin 2t
#     ur  = s/(4G) [ r ((k-1)/2 + cos 2t) + a2/r (1 + (1+k) cos 2t) - a4/r3 cos 2t ]
#     ut  = -s/(4G) [ r + (1-k) a2/r + a4/r3 ] sin 2t,     k = 3 - 4 nu
#   la traction selon y s'obtient par t -> t - pi/2 et superposition.
#   Concentration en (0, a) : sxx = (3 - LAMBDA) SX.
#   Metriques : fem3d, sondes `probes` (contrainte constante par tet) le long
#   de l'axe y (theta = 90 deg) et de l'axe x ; fem3d et fdem3d, deplacement
#   des sommets de paroi (point.wA en (a,0), point.wB en (0,a), force nulle).
#
# Quasi statique : scenario = loads, amplitude `ramp`, dampingLocal (Cundall),
# duree en nombre de traversees L/c de l'onde P.
#
#   python3 vv/P1_elastique/p1_elastique.py prepare
#   python3 vv/P1_elastique/p1_elastique.py run     [--only patch_fem3d ...] [--threads 2]
#   python3 vv/P1_elastique/p1_elastique.py analyse
#   python3 vv/P1_elastique/p1_elastique.py all
#   python3 vv/run_queue.py P1_elastique --slots 4 --threads 2
# Option --smoke : maillage le plus grossier de chaque variante, sorties dans
# out_smoke/ (non lues par analyse sauf --smoke).
# ---------------------------------------------------------------------------
import argparse, json, math, os, re, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
VVDIR = os.path.dirname(HERE)
ROOT = os.path.dirname(VVDIR)
sys.path.insert(0, VVDIR)
import vvcommon as C  # noqa: E402

OUT = os.path.join(HERE, "out")
MESHES = os.path.join(HERE, "meshes")

# ---- materiau et chargement (fixes AVANT le calcul) ------------------------
RHO, E, NU = 2650.0, 60e9, 0.25
G = E / (2 * (1 + NU))
KAPPA = 3 - 4 * NU                         # deformation plane
CP = math.sqrt(E * (1 - NU) / ((1 + NU) * (1 - 2 * NU)) / RHO)   # onde P, 5 212 m/s
DAMP = 0.5                                 # dampingLocal (Cundall)
# fdem3d : insertion adaptative sans rupture (ft >> contraintes) = continu
# exact (V1, P1.1). NB : `dtUpdate = inserted` (doc §5.1, pas initial sans les
# joints lies) est refusee par build_nofma/rockim (cle inconnue) : pas de temps
# par defaut, ~5 fois plus petit qu'en fem3d sur ces maillages.
FDEM_KEYS = ["insertion = adaptive", "ft = 50e6", "cohesion = 100e6",
             "frictionDeg = 40", "Gf = 70"]

# P1.3 cube biaxial
L_CUBE = 1.0
PX, PY = 1e6, 0.5e6                        # Pa, compression (eps ~ 1,5e-5)
H_CUBE = [0.25, 0.125, 0.0625]             # L/h = 4, 8, 16

# P1.2 Kirsch
A_HOLE = 0.1                               # rayon du trou
W_PLATE = 2.0                              # W / a = 20
SX = -1e6                                  # traction lointaine selon x (compression)
LAMBDAS = [0.0, 0.5]
S_KIRSCH = [1.0, 0.5, 0.25]                # facteur d'echelle des tailles
HW0 = A_HOLE / 10                          # taille a la paroi pour s = 1 (a/10, a/20, a/40)
GROW = 1.5                                 # h(r) = hw (1 + GROW (r - a)/a)
HCAP0 = 0.25                               # plafond de taille pour s = 1
NLAYERS = 2                                # couches en z, epaisseur t = 2 hw

# duree : montee en cosinus sur N_RAMP traversees L/c, total N_T traversees.
# Mesure du smoke test (maillages grossiers, dampingLocal = 0,5) : U a 1e-5
# pres de sa valeur finale apres 40 traversees (cube) et 25 (Kirsch, L = W).
N_RAMP = 5
N_T = dict(patch=60, kirsch=40)

# ---- variantes ---------------------------------------------------------------
VARIANTS = {
    "patch_fem3d":  dict(kind="patch", mode="fem3d", sizes=H_CUBE),
    "patch_fdem3d": dict(kind="patch", mode="fdem3d", sizes=H_CUBE),
    "kirsch_fem3d_l0":   dict(kind="kirsch", mode="fem3d", lam=0.0, sizes=S_KIRSCH),
    "kirsch_fem3d_l0p5": dict(kind="kirsch", mode="fem3d", lam=0.5, sizes=S_KIRSCH),
    # fdem3d coute ~50 fois fem3d par tet : les deux maillages les plus grossiers
    "kirsch_fdem3d_l0":  dict(kind="kirsch", mode="fdem3d", lam=0.0, sizes=S_KIRSCH[:2]),
}

# ---- criteres d'acceptation (fixes AVANT le calcul, 2026-10-03) ------------
CRIT = dict(
    # P1.3 : chaque maillage, chaque mode
    patch_nodal=1e-6,        # max |u - u_ex| / max |u_ex| sur tous les noeuds (VTU)
    patch_group=2e-5,        # U de groupe (history.csv, 6 chiffres significatifs)
    # P1.2 : maillage le plus fin de chaque variante
    kt_centroid=0.02,        # contrainte du tet de paroi contre l'exact a SON centroide
    kt_raw=0.05,             # sxx(0, a) du tet de paroi contre (3 - lambda) SX
    line_centroid=0.03,      # L2 relative de sxx le long de theta = 90 deg, r in [a, 4a], exact aux centroides
    u_wall=0.01,             # deplacement des sommets de paroi wA, wB
    order_min=0.8,           # ordre observe de l'erreur L2 brute de sxx (theorique 1, tets a contrainte constante)
    stationary=1e-4,         # variation relative de U_wB sur les 10 derniers % de T (etat statique atteint)
)


# ---- solutions exactes -------------------------------------------------------
def patch_exact(xyz):
    sxx, syy = -PX, -PY
    exx = (sxx - NU * syy) / E
    eyy = (syy - NU * sxx) / E
    ezz = -NU * (sxx + syy) / E
    return xyz * np.array([exx, eyy, ezz])


def _kirsch_uni(r, th, s, a=A_HOLE):
    c2, s2 = np.cos(2 * th), np.sin(2 * th)
    q2, q4 = a * a / (r * r), a ** 4 / r ** 4
    srr = s / 2 * (1 - q2) + s / 2 * (1 - 4 * q2 + 3 * q4) * c2
    stt = s / 2 * (1 + q2) - s / 2 * (1 + 3 * q4) * c2
    srt = -s / 2 * (1 + 2 * q2 - 3 * q4) * s2
    ur = s / (4 * G) * (r * ((KAPPA - 1) / 2 + c2) + a * a / r * (1 + (1 + KAPPA) * c2)
                        - a ** 4 / r ** 3 * c2)
    ut = -s / (4 * G) * (r + (1 - KAPPA) * a * a / r + a ** 4 / r ** 3) * s2
    return srr, stt, srt, ur, ut


def kirsch(x, y, lam):
    """contraintes cartesiennes (sxx, syy, sxy) et deplacements (ux, uy), plan"""
    x, y = np.asarray(x, float), np.asarray(y, float)
    r, th = np.hypot(x, y), np.arctan2(y, x)
    a1 = _kirsch_uni(r, th, SX)
    a2 = _kirsch_uni(r, th - np.pi / 2, lam * SX)
    srr, stt, srt, ur, ut = (p + q for p, q in zip(a1, a2))
    c, s = np.cos(th), np.sin(th)
    sxx = srr * c * c + stt * s * s - 2 * srt * s * c
    syy = srr * s * s + stt * c * c + 2 * srt * s * c
    sxy = (srr - stt) * s * c + srt * (c * c - s * s)
    return sxx, syy, sxy, ur * c - ut * s, ur * s + ut * c


# ---- maillages -----------------------------------------------------------------
def tag(x):
    return f"{x:g}".replace(".", "p")


def cube_mesh(h):
    p = os.path.join(MESHES, f"cube_h{tag(h)}.msh")
    if not os.path.exists(p):
        os.makedirs(MESHES, exist_ok=True)
        subprocess.run([sys.executable, os.path.join(ROOT, "tools", "make_unstructured_mesh.py"),
                        "box3dbc", str(L_CUBE), str(L_CUBE), str(L_CUBE), str(h), p, "1"],
                       check=True, cwd=ROOT)
    return p


def kirsch_geom(s):
    hw = s * HW0
    return dict(hw=hw, hcap=s * HCAP0, t=NLAYERS * hw)


def kirsch_mesh(s):
    p = os.path.join(MESHES, f"kirsch_s{tag(s)}.msh")
    if os.path.exists(p):
        return p
    import gmsh
    os.makedirs(MESHES, exist_ok=True)
    g = kirsch_geom(s)
    a, W, t = A_HOLE, W_PLATE, g["t"]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("kirsch")
    occ = gmsh.model.occ
    rect = occ.addRectangle(0, 0, 0, W, W)
    disk = occ.addDisk(0, 0, 0, a, a)
    surf, _ = occ.cut([(2, rect)], [(2, disk)])
    occ.synchronize()
    ext = occ.extrude(surf, 0, 0, t, numElements=[NLAYERS], recombine=False)
    occ.synchronize()
    vols = [tg for d, tg in gmsh.model.getEntities(3)]
    gmsh.model.addPhysicalGroup(3, vols, name="solid")
    tol = 1e-6 * W
    groups = {}
    for _, sf in gmsh.model.getEntities(2):
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, sf)
        if z1 < tol:
            nm = "bottom"
        elif z0 > t - tol:
            nm = "top"
        elif x1 < tol:
            nm = "xmin"
        elif y1 < tol:
            nm = "ymin"
        elif x0 > W - tol:
            nm = "xmax"
        elif y0 > W - tol:
            nm = "ymax"
        else:
            nm = "hole"
        groups.setdefault(nm, []).append(sf)
    for nm, sfs in groups.items():
        gmsh.model.addPhysicalGroup(2, sfs, name=nm)
    f = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(
        f, "F", f"Min({g['hcap']}, {g['hw']}*(1+{GROW}*(Sqrt(x*x+y*y)-{a})/{a}))")
    gmsh.model.mesh.field.setAsBackgroundMesh(f)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 5)
    gmsh.option.setNumber("Mesh.RandomSeed", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(3)
    gmsh.write(p)
    gmsh.finalize()
    return p


# ---- decks -------------------------------------------------------------------
def outdir(var, size, smoke=False):
    return os.path.join(HERE, "out_smoke" if smoke else "out", f"{var}_{tag(size)}")


def mode_keys(mode):
    return ["mode = fem3d", "law = elastic"] if mode == "fem3d" else ["mode = fdem3d", *FDEM_KEYS]


def times(Lchar, kind, tscale=1.0):
    tr = Lchar / CP
    return N_RAMP * tr, N_T[kind] * tr * tscale


def probe_points(lam, g):
    """sondes : axe y (theta = 90 deg) r in [a, 4a], resserrees a la paroi ;
    axe x (theta = 0) idem, moins dense"""
    a, z = A_HOLE, g["t"] / 4
    eps = 1e-6
    ry = a + 3 * a * (np.arange(25) / 24.0) ** 2
    rx = a + 3 * a * (np.arange(10) / 9.0) ** 2
    ry[0] = rx[0] = a * (1 + 1e-6)
    pts = [(eps, r, z) for r in ry] + [(r, eps, z) for r in rx]
    return pts, len(ry)


def deck(var, size, od, tscale=1.0):
    v = VARIANTS[var]
    lines = [f"# P1 elastique, {var}, taille {size:g} (ecrit par p1_elastique.py)",
             *mode_keys(v["mode"]), "scenario = loads", "mesh = file",
             f"rho = {RHO}", f"E = {E}", f"nu = {NU}", f"dampingLocal = {DAMP}",
             "frames = 2", "historyFlush = true"]
    if v["kind"] == "patch":
        tr, T = times(L_CUBE, "patch", tscale)
        lines += [f"meshFile = {cube_mesh(size)}", f"T = {T:.6g}",
                  "fix.xmin = x", "fix.ymin = y", "fix.bottom = z",
                  f"traction.xmax = {-PX:g} 0 0", f"traction.ymax = 0 {-PY:g} 0",
                  f"amplitude.xmax = ramp {tr:.6g}", f"amplitude.ymax = ramp {tr:.6g}",
                  "force.top = 0 0 0"]
    else:
        g = kirsch_geom(size)
        tr, T = times(W_PLATE, "kirsch", tscale)
        lam = v["lam"]
        a, zc = A_HOLE, g["t"] / 2
        lines += [f"meshFile = {kirsch_mesh(size)}", f"T = {T:.6g}",
                  "fix.xmin = x", "fix.ymin = y", "fix.bottom = z", "fix.top = z",
                  f"traction.xmax = {SX:g} 0 0", f"amplitude.xmax = ramp {tr:.6g}"]
        if lam != 0:
            lines += [f"traction.ymax = 0 {lam * SX:g} 0", f"amplitude.ymax = ramp {tr:.6g}"]
        for nm, (x, y) in dict(wA=(a, 0), wB=(0, a), mA=(2 * a, 0), mB=(0, 2 * a)).items():
            lines += [f"point.{nm} = {x:.9g} {y:.9g} {zc:.9g}", f"force.{nm} = 0 0 0"]
        if v["mode"] == "fem3d":
            pts, _ = probe_points(lam, g)
            lines.append("probes = " + " ; ".join(f"{x:.9g},{y:.9g},{z:.9g}" for x, y, z in pts))
    return C.write_deck(od, lines)


def weight(var, size):
    """cout estime en s a 1 fil (smoke test : cube h = 0,25 m 0,16 s fem3d, 4,6 s
    fdem3d ; Kirsch s = 1, 40 traversees, ~30 s fem3d, ~2 950 s fdem3d) ;
    cube ~ h^-4 (tets h^-3, pas h^-1), Kirsch ~ s^-3 (tets s^-2, pas s^-1)"""
    v = VARIANTS[var]
    if v["kind"] == "patch":
        return (0.16 if v["mode"] == "fem3d" else 4.6) * (0.25 / size) ** 4
    return (30.0 if v["mode"] == "fem3d" else 2950.0) * (1.0 / size) ** 3


def JOBS(args=None, smoke=False, only=None, tscale=1.0):
    jobs = []
    for var, v in VARIANTS.items():
        if only and var not in only:
            continue
        sizes = [v["sizes"][0]] if smoke else v["sizes"]
        for s in sizes:
            od = outdir(var, s, smoke)
            cfg = deck(var, s, od, tscale)
            jobs.append(C.Job("P1", f"{var}_{tag(s)}", cfg, od, weight(var, s),
                              dict(var=var, size=s)))
    return jobs


# ---- lecture des sorties ---------------------------------------------------
def vtu_points(path):
    txt = open(path).read()
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", txt, re.S)
    return np.array(m.group(1).split(), dtype=float).reshape(-1, 3)


def vtu_frames(od):
    fs = sorted(f for f in os.listdir(od) if re.match(r"(fem3d|fdem3d)_\d+\.vtu$", f))
    return [os.path.join(od, f) for f in fs]


def log_points(log):
    out = {}
    for m in re.finditer(r"point\.(\w+) : sommet \d+ a \S+ m du point demande \(\s*(\S+)\s+(\S+)\s+(\S+),", log):
        out[m.group(1)] = np.array([float(m.group(i)) for i in (2, 3, 4)])
    return out


def log_probes(log):
    cen = {}
    for m in re.finditer(r"probe (\d+) : \(.*?\) m -> element (\d+) \(centroide (\S+), (\S+), (\S+), lc", log):
        cen[int(m.group(1))] = np.array([float(m.group(i)) for i in (3, 4, 5)])
    return cen


def log_translation(log):
    m = re.search(r"translat\w*[^\n]*?\(\s*(\S+)\s*,?\s+(\S+)\s*,?\s+(\S+)\s*\)", log)
    return m.group(0) if m else None


def stationarity(t, u):
    n = len(t)
    i0 = int(0.9 * n)
    ref = max(abs(u[-1]), 1e-300)
    return float(np.max(np.abs(u[i0:] - u[-1])) / ref)


def analyse_patch(od, mode):
    lg = C.parse_log(od)
    fr = vtu_frames(od)
    X0, X1 = vtu_points(fr[0]), vtu_points(fr[-1])
    U = X1 - X0
    Ue = patch_exact(X0 - X0.min(axis=0))
    nodal = float(np.max(np.linalg.norm(U - Ue, axis=1)) / np.max(np.linalg.norm(Ue, axis=1)))
    h = C.read_history(od)
    ex = patch_exact(np.array([L_CUBE, L_CUBE, L_CUBE]))
    grp = dict(xmax_x=float(h["U_xmax_x"][-1] / ex[0] - 1),
               ymax_y=float(h["U_ymax_y"][-1] / ex[1] - 1),
               top_z=float(h["U_top_z"][-1] / ex[2] - 1))
    return dict(ntet=lg["ntet"], dt=lg["dt"], wall=lg["wall"], budget_pct=lg["budget_pct"],
                nnodes=int(len(X0)), nodal=nodal, group=grp,
                stationary=stationarity(h["t"], h["U_xmax_x"]),
                passe=bool(nodal <= CRIT["patch_nodal"]
                           and max(abs(x) for x in grp.values()) <= CRIT["patch_group"]))


def analyse_kirsch(od, var, size):
    v = VARIANTS[var]
    lam = v["lam"]
    lg = C.parse_log(od)
    log = lg["text"]
    h = C.read_history(od)
    pts = log_points(log)
    r = dict(ntet=lg["ntet"], dt=lg["dt"], wall=lg["wall"], budget_pct=lg["budget_pct"],
             hw=kirsch_geom(size)["hw"])
    _, _, _, uxA, _ = kirsch(pts["wA"][0], pts["wA"][1], lam)
    _, _, _, _, uyB = kirsch(pts["wB"][0], pts["wB"][1], lam)
    r["u_wA_err"] = float(h["U_wA_x"][-1] / uxA - 1)
    r["u_wB_err"] = float(h["U_wB_y"][-1] / uyB - 1)
    _, _, _, uxm, _ = kirsch(pts["mA"][0], pts["mA"][1], lam)
    _, _, _, _, uym = kirsch(pts["mB"][0], pts["mB"][1], lam)
    r["u_mA_err"] = float(h["U_mA_x"][-1] / uxm - 1)
    r["u_mB_err"] = float(h["U_mB_y"][-1] / uym - 1)
    r["stationary"] = stationarity(h["t"], h["U_wB_y"])
    r["uz_mid_max"] = float(max(abs(h["U_wA_z"][-1]), abs(h["U_wB_z"][-1])))
    if v["mode"] == "fem3d":
        g = kirsch_geom(size)
        P, ny = probe_points(lam, g)
        P = np.array(P)
        cen = log_probes(log)
        pr = np.genfromtxt(os.path.join(od, "probes.csv"), delimiter=",", names=True)
        last = pr[-1]
        sxx = np.array([last[f"p{k + 1}_sxx"] for k in range(len(P))])
        syy = np.array([last[f"p{k + 1}_syy"] for k in range(len(P))])
        Cc = np.array([cen[k + 1] for k in range(len(P))])
        ex_p = kirsch(P[:, 0], P[:, 1], lam)
        ex_c = kirsch(Cc[:, 0], Cc[:, 1], lam)
        iy = slice(0, ny)
        kt_ex = (3 - lam) * SX
        r["kt_exact"] = 3 - lam
        r["kt_h"] = float(sxx[0] / SX)
        r["kt_raw_err"] = float(sxx[0] / kt_ex - 1)
        r["kt_centroid_err"] = float(sxx[0] / ex_c[0][0] - 1)
        r["wall_centroid_r"] = float(np.hypot(*Cc[0, :2]))
        r["line_raw"] = float(np.linalg.norm(sxx[iy] - ex_p[0][iy]) / np.linalg.norm(ex_p[0][iy]))
        r["line_centroid"] = float(np.linalg.norm(sxx[iy] - ex_c[0][iy]) / np.linalg.norm(ex_c[0][iy]))
        ix = slice(ny, len(P))
        r["xaxis_syy_wall"] = float(syy[ny])
        r["xaxis_syy_wall_exact_c"] = float(ex_c[1][ny])
        r["profile"] = dict(r=np.hypot(P[iy, 0], P[iy, 1]).tolist(), sxx=sxx[iy].tolist(),
                            r_c=np.hypot(Cc[iy, 0], Cc[iy, 1]).tolist(),
                            sxx_ex_c=ex_c[0][iy].tolist())
    return r


def order(hs, es):
    hs, es = np.log(np.asarray(hs, float)), np.log(np.abs(np.asarray(es, float)))
    return float(np.polyfit(hs, es, 1)[0]) if len(hs) >= 2 else float("nan")


def analyse(smoke=False):
    res = {}
    for var, v in VARIANTS.items():
        runs = {}
        for s in v["sizes"]:
            od = outdir(var, s, smoke)
            if not C.is_done(od):
                continue
            runs[f"{s:g}"] = (analyse_patch(od, v["mode"]) if v["kind"] == "patch"
                              else analyse_kirsch(od, var, s))
        if not runs:
            continue
        ent = dict(runs=runs)
        if v["kind"] == "patch":
            ent["passe"] = all(r["passe"] for r in runs.values())
        else:
            fin = runs[min(runs, key=float)]
            vd = dict(u_wall=max(abs(fin["u_wA_err"]), abs(fin["u_wB_err"])) <= CRIT["u_wall"],
                      stationary=fin["stationary"] <= CRIT["stationary"])
            if v["mode"] == "fem3d":
                vd["kt_centroid"] = abs(fin["kt_centroid_err"]) <= CRIT["kt_centroid"]
                vd["kt_raw"] = abs(fin["kt_raw_err"]) <= CRIT["kt_raw"]
                vd["line_centroid"] = fin["line_centroid"] <= CRIT["line_centroid"]
                if len(runs) >= 3:
                    hs = sorted(runs, key=float)
                    ent["ordre_line_raw"] = order([runs[k]["hw"] for k in hs],
                                                  [runs[k]["line_raw"] for k in hs])
                    ent["ordre_u_wB"] = order([runs[k]["hw"] for k in hs],
                                              [runs[k]["u_wB_err"] for k in hs])
                    vd["ordre"] = ent["ordre_line_raw"] >= CRIT["order_min"]
            ent["verdict"] = vd
            ent["passe"] = all(vd.values())
        res[var] = ent
    ref = dict(rho=RHO, E=E, nu=NU, cp=CP, damping=DAMP, L_cube=L_CUBE, PX=PX, PY=PY,
               a=A_HOLE, W=W_PLATE, SX=SX, n_ramp=N_RAMP, n_T=N_T, criteres=CRIT)
    name = "resultats_smoke.json" if smoke else "resultats.json"
    with open(os.path.join(HERE, name), "w") as f:
        json.dump(dict(reference=ref, variantes=res), f, indent=1)
    report(res)
    if not smoke:
        figures(res)
    return res


def report(res):
    for var, ent in res.items():
        print(f"\n{var}  ->  {'PASSE' if ent['passe'] else 'ECHEC'}")
        for s, r in sorted(ent["runs"].items(), key=lambda kv: -float(kv[0])):
            if "nodal" in r:
                print(f"  h {s:>7s}  ntet {r['ntet']}  nodal {r['nodal']:.2e}  groupes "
                      + " ".join(f"{k} {x:+.1e}" for k, x in r["group"].items())
                      + f"  stat {r['stationary']:.1e}  {r['wall']} s")
            else:
                line = (f"  s {s:>5s}  ntet {r['ntet']}  u_wA {r['u_wA_err']:+.4f}  u_wB {r['u_wB_err']:+.4f}"
                        f"  u_mA {r['u_mA_err']:+.4f}  u_mB {r['u_mB_err']:+.4f}  stat {r['stationary']:.1e}")
                if "kt_h" in r:
                    line += (f"  Kt {r['kt_h']:.4f} (exact {r['kt_exact']:.2f}, brut {r['kt_raw_err']:+.4f},"
                             f" centroide {r['kt_centroid_err']:+.4f})  L2 brut {r['line_raw']:.4f}"
                             f"  L2 centr {r['line_centroid']:.4f}")
                print(line + f"  {r['wall']} s")
        if "verdict" in ent:
            print(f"  verdict {ent['verdict']}  ordres "
                  f"{ent.get('ordre_line_raw', float('nan')):.2f} / {ent.get('ordre_u_wB', float('nan')):.2f}")


def figures(res):
    plt = C.plot_style()
    # 1. patch test : erreur nodale contre h
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    for var in ("patch_fem3d", "patch_fdem3d"):
        if var in res:
            rr = res[var]["runs"]
            hs = sorted(rr, key=float)
            ax.semilogy([float(h) for h in hs], [rr[h]["nodal"] for h in hs], "o-", label=var)
    ax.axhline(CRIT["patch_nodal"], color="k", ls="--", lw=0.8, label="critère")
    ax.set_xlabel("taille de maille h (m)")
    ax.set_ylabel("erreur nodale max relative")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    C.savefig(fig, os.path.join(HERE, "fig_patch"))
    plt.close(fig)
    # 2. Kirsch : profil sxx le long de theta = 90 deg (maillage le plus fin)
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for k, var in enumerate(("kirsch_fem3d_l0", "kirsch_fem3d_l0p5")):
        if var not in res:
            continue
        lam = VARIANTS[var]["lam"]
        rr = np.linspace(A_HOLE, 4 * A_HOLE, 300)
        ax.plot(rr / A_HOLE, kirsch(1e-9 * rr, rr, lam)[0] / SX, color=f"C{k}", lw=2, alpha=0.35)
        for s, r in res[var]["runs"].items():
            if "profile" in r and float(s) == min(map(float, res[var]["runs"])):
                p = r["profile"]
                ax.plot(np.array(p["r_c"]) / A_HOLE, np.array(p["sxx"]) / SX, "o", ms=3,
                        color=f"C{k}", label=f"λ = {lam:g}, s = {s}")
    ax.set_xlabel("r / a  (θ = 90°)")
    ax.set_ylabel(r"$\sigma_{xx} / \sigma_x^\infty$")
    ax.set_title("Kirsch : trait large exact, points rockim (centroïdes)", fontsize=9)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    C.savefig(fig, os.path.join(HERE, "fig_kirsch_profil"))
    plt.close(fig)
    # 3. Kirsch : convergence
    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    for var, ent in res.items():
        if not var.startswith("kirsch"):
            continue
        rr = ent["runs"]
        hs = sorted(rr, key=float)
        x = [rr[h]["hw"] / A_HOLE for h in hs]
        if "line_raw" in rr[hs[0]]:
            ax.loglog(x, [rr[h]["line_raw"] for h in hs], "o-", label=f"{var}, L2 σxx")
            ax.loglog(x, [abs(rr[h]["kt_raw_err"]) for h in hs], "s--", label=f"{var}, |Kt − Kt_ex|/Kt_ex")
        ax.loglog(x, [abs(rr[h]["u_wB_err"]) for h in hs], "^:", label=f"{var}, u paroi")
    ax.set_xlabel("h paroi / a")
    ax.set_ylabel("erreur relative")
    ax.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    C.savefig(fig, os.path.join(HERE, "fig_kirsch_convergence"))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["prepare", "run", "analyse", "all"])
    ap.add_argument("--exe", default=C.EXE_DEFAULT)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--only", nargs="+", default=None)
    ap.add_argument("--smoke", action="store_true", help="maillage grossier seul, out_smoke/")
    ap.add_argument("--tscale", type=float, default=1.0, help="T multiplie (essai court)")
    a = ap.parse_args()
    if a.action in ("prepare", "run", "all"):
        jobs = JOBS(None, smoke=a.smoke, only=a.only, tscale=a.tscale)
        print(f"[P1] {len(jobs)} decks ecrits")
    if a.action in ("run", "all"):
        for j in sorted(jobs, key=lambda j: j.weight):
            rc, wall = C.run_job(j, os.path.abspath(a.exe), a.threads)
            print(f"  {j.name:30s} rc = {rc}  {wall:7.1f} s", flush=True)
    if a.action in ("analyse", "all"):
        analyse(smoke=a.smoke)


if __name__ == "__main__":
    main()
