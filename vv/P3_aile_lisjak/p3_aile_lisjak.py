#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# P3 aile : amorcage de fissures en aile depuis une fissure inclinee sous
# compression uniaxiale (benchmark B2 et validation physique P3.3).
#
# Reference : Lisjak (2013), these Univ. Toronto, §3.4.3, pp. 36-39 (Y-Geo,
# FDEM 2D). Eprouvette 50 x 100 mm, fissure preexistante de 5 mm a 45 deg au
# centre, compression par deux platines rigides frottantes (k = 0,1).
# Y-Geo : sigma_1c = 7,0 MPa, angle de branchement 64 deg ; MLER : 5,1 a
# 14,6 MPa (Table 3.3) ; contrainte circonferentielle max : 70,5 deg.
#
# Modelisation rockim (mode fdem, deformation plane, thickness = 1) :
#   - maillage Delaunay gmsh (Mesh.Algorithm = 5) de la boite, la fissure etant
#     une LIGNE INCLUSE (embed) : ses aretes sont exactement sur le segment ;
#   - fissure = preBrokenJoints (§5.11) sur ce segment, tolerance 1e-6 m et
#     5 deg : seules les aretes de la ligne incluse naissent rompues (D = 1),
#     frottement residuel jointResidualMu = tan(35 deg) = 0,70 ;
#   - scenario = tension, loading = platens, pullV < 0 (fermeture totale) ;
#     Lisjak : deux platines a 0,05 m/s en sens opposes = fermeture 0,1 m/s
#     (coherent avec sigma = E' eps = 6,1 MPa a t = 1,87 ms, evenements 1-2) ;
#   - frottement platine-echantillon : contactMu = 0,1 (la meme cle regit le
#     contact general des fissures neuves : ecart assume, voir README) ;
#   - amortissement : Lisjak eq. 3.2 C = mu I, mu = 7,4e3 kg/(m s) = 2 mu_c
#     (eq. 3.21), convention de Munjiza sigma_v = mu D. rockim ecrit 2 mu D
#     (DOCUMENTATION §5.4 quinquies, meme lecture pour Solidity) : on pose
#     donc bulkViscosity = mu / 2 = 3,7e3 Pa s (corrige le 2026-10-04 a la
#     relecture, avant tout calcul de campagne). dampingLocal = 0 (Y-Geo n'a pas de Cundall), plafond de
#     traction moyenne coupe (meanTensionCapFactor = 0, regle de replique).
# Detection (fixee AVANT le calcul) :
#   amorcage = premier joint NON preexistant rompu (tBreak > 0 dans
#   fdem_final_joints.csv) dont le milieu est a moins de R_TIP = 1,5 h d'une
#   pointe de la fissure ; sigma_1c = sigma(t) de history.csv (reaction de la
#   platine superieure / W) interpolee a cet instant. Indicateur secondaire :
#   premier joint entre en endommagement pres d'une pointe (tYield du
#   catalogue microsismique fdem_seismic.csv), borne basse.
#   angle de branchement = angle SIGNE (trigonometrique) entre la tangente
#   sortante de la fissure a la pointe et le vecteur pointe -> extremite
#   lointaine du premier joint rompu (positif = cote en traction attendu) ;
#   trajet = amas connexe de joints rompus issu de chaque pointe en fin de
#   calcul (longueur de corde, orientation par rapport a sigma_1).
#
#   python3 vv/P3_aile_lisjak/p3_aile_lisjak.py prepare|run|analyse|all [--smoke]
#   python3 vv/run_queue.py P3_aile_lisjak --slots 4 --threads 2
# ---------------------------------------------------------------------------
import argparse, csv, json, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import vvcommon as C          # noqa: E402
import lefm_aile as LE        # noqa: E402

OUT = os.path.join(HERE, "out")
MESHES = os.path.join(HERE, "meshes")
BENCH = "P3.aile"

# ---- probleme (fixe AVANT le calcul) ---------------------------------------
W, H = 0.050, 0.100
XC, YC = 0.5 * W, 0.5 * H
FLAW_L, GAMMA = 5e-3, 45.0                     # 2c, inclinaison sur sigma_1 (vertical)
_d = 0.5 * FLAW_L / math.sqrt(2.0)
TIP_A = np.array([XC - _d, YC - _d])          # pointe basse-gauche
TIP_B = np.array([XC + _d, YC + _d])          # pointe haute-droite
MAT = ["rho = 2300", "E = 3e9", "nu = 0.29", "ft = 3e6", "cohesion = 15e6",
       "frictionDeg = 35", "Gf = 2.0", "gfShearFactor = 5"]       # G_IIc = 10 J/m2
CLOSURE_V = 0.1                               # m/s, fermeture totale des platines
T_REF = 3.5e-3                                # ~1,8 x l'amorcage Y-Geo (1,87 ms), << 9 ms
COMMON = [
    "mode = fdem", "scenario = tension", "loading = platens", "mesh = file",
    "thickness = 1.0", "pullRamp = 1e-4", "verifyFt = false",
    "jointSoftening = yan",                    # z-curve de Munjiza (eq. 3.4, A B C = 0,63 1,8 6,0)
    "jointPenaltyFactor = 20", "jointXi = 0.01",
    "contactMu = 0.1",                         # frottement platine (Table 3.1)
    "jointResidualMu = 0.70",                  # tan(phi_f = 35 deg), fissure preexistante
    "preBrokenTol = 1e-6", "preBrokenAngleDeg = 5",
    "bulkViscosity = 3.7e3",                   # mu/2 : Lisjak mu = 7,4e3 = 2 mu_c, rockim ecrit 2 mu D
    "dampingLocal = 0", "meanTensionCapFactor = 0",
    "microseismic = true", "writeJointMode = true",
    "gcSurfaceRefresh = eager",                # recommande par le solveur avec preBrokenJoints
]

# variantes : h en mm, fermeture en m/s, cles supplementaires
VARIANTS = {
    "h0p7":       dict(h=0.7, v=CLOSURE_V, T=T_REF, ref=True),             # Lisjak
    "h1p0":       dict(h=1.0, v=CLOSURE_V, T=T_REF),                       # plus grossier
    "h0p5":       dict(h=0.5, v=CLOSURE_V, T=T_REF),                       # plus fin
    "h1p0_v0p05": dict(h=1.0, v=0.5 * CLOSURE_V, T=2 * T_REF),             # vitesse / 2
    "h1p0_pf2p5": dict(h=1.0, v=CLOSURE_V, T=T_REF, keys=["jointPenaltyFactor = 2.5"]),
    # Lisjak pf = 15 GPa = 5 E : equivalent rockim 5 si la raideur Y-Geo vaut pf/h,
    # 2,5 si elle vaut pf/(2h) ; les deux encadrent (relecture 2026-10-04)
    "h1p0_pf5":   dict(h=1.0, v=CLOSURE_V, T=T_REF, keys=["jointPenaltyFactor = 5"]),
}
SMOKE = {"smoke_h2p5": dict(h=2.5, v=CLOSURE_V, T=2.5e-3, frames=5),
         "smoke_h2_v0p5": dict(h=2.0, v=0.5, T=1.2e-3, frames=6)}
# taille gmsh = GMSH_SIZE x h : calee pour que l'ARETE MOYENNE vaille h
# (mesure : taille 0,7 mm -> arete moyenne 0,643 mm). Lisjak : arete moyenne
# 0,70 mm et 21 600 triangles ; ici ~24 000 a h = 0,7 mm.
GMSH_SIZE = 1.09

R_TIP = 1.5          # rayon de detection autour d'une pointe, en multiples de h

# ---- criteres d'acceptation (fixes AVANT le calcul, 2026-10-04) -------------
# appliques a la variante de reference h0p7
CRIT = dict(
    lefm_min=5.1, lefm_max=14.6,              # P3.3 : dans la fourchette MLER publiee
    ygeo=7.0, ygeo_tol=0.15,                  # B2 : a 15 % de Y-Geo
    theta_refs=[64.0, 70.53], theta_tol=10.0, # B2 / P3.3 : angle a 10 deg de 64 ou 70,5
    two_tips=True,                            # les deux pointes amorcent
    mesh_tol=0.15,                            # sensibilite : |s1c(h) / s1c(0,7) - 1| <= 15 %
    rate_tol=0.10,                            # quasi-statisme : vitesse / 2 change s1c de <= 10 %
)


# ---- maillage -------------------------------------------------------------
def tag(h):
    return f"{h:g}".replace(".", "p")


def mesh_path(h):
    return os.path.join(MESHES, f"aile_h{tag(h)}.msh")


def ensure_mesh(h_mm, seed=1):
    p = mesh_path(h_mm)
    if os.path.exists(p):
        return p
    import gmsh
    os.makedirs(MESHES, exist_ok=True)
    h = h_mm * 1e-3 * GMSH_SIZE
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("aile")
    s = gmsh.model.occ.addRectangle(0, 0, 0, W, H)
    a = gmsh.model.occ.addPoint(TIP_A[0], TIP_A[1], 0)
    b = gmsh.model.occ.addPoint(TIP_B[0], TIP_B[1], 0)
    ln = gmsh.model.occ.addLine(a, b)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.embed(1, [ln], 2, s)       # aretes conformes a la fissure
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.RandomSeed", seed)
    gmsh.option.setNumber("Mesh.Algorithm", 5)  # Delaunay (regle maison, jamais 6)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.model.addPhysicalGroup(2, [s], name="roche")
    gmsh.model.mesh.generate(2)
    gmsh.write(p)
    _, _, conn = gmsh.model.mesh.getElements(2)
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    idx = {t: i for i, t in enumerate(tags)}
    P = np.array(xyz).reshape(-1, 3)[[idx[t] for t in conn[0]], :2].reshape(-1, 3, 2)
    e = np.linalg.norm(P[:, [1, 2, 0]] - P, axis=2)
    print(f"[P3 aile] maillage h = {h_mm} mm : {len(conn[0]) // 3} triangles, arete moyenne "
          f"{e.mean() * 1e3:.3f} mm -> {p}")
    gmsh.finalize()
    return p


# ---- decks et Job -----------------------------------------------------------
def outdir(name):
    return os.path.join(OUT, name)


def deck_lines(name, var):
    h = var["h"]
    lines = list(COMMON) + MAT + [
        f"meshFile = {mesh_path(h)}",
        f"T = {var['T']:g}", f"frames = {var.get('frames', 35)}",
        f"pullV = {-var['v']:g}",
        f"preBrokenJoints = {TIP_A[0]:.9g} {TIP_A[1]:.9g} {TIP_B[0]:.9g} {TIP_B[1]:.9g}",
    ]
    for k in var.get("keys", []):              # surcharge : la derniere ligne gagne ?
        key = k.split("=")[0].strip()
        lines = [l for l in lines if l.split("=")[0].strip() != key] + [k]
    return lines


def make_job(name, var):
    ensure_mesh(var["h"])
    od = outdir(name)
    cfg = C.write_deck(od, deck_lines(name, var),
                       header=f"P3 aile Lisjak, variante {name} (ecrit par p3_aile_lisjak.py)")
    # cout ~ N_tri x N_pas ~ h^-2 x h^-2 x T : le pas est commande par la borne
    # diffusive rho h^2 / (4 mu) de la viscosite (journal du smoke), pas par h/c
    w = (0.7 / var["h"]) ** 4 * var["T"] / T_REF
    return C.Job(BENCH, name, cfg, od, weight=w, meta=dict(var))


def JOBS(args=None):
    return [make_job(n, v) for n, v in VARIANTS.items()]


# ---- depouillement ----------------------------------------------------------
def read_csv(path):
    with open(path) as f:
        r = csv.DictReader(f)
        rows = list(r)
    return {k: np.array([float(x[k]) for x in rows]) for k in (rows[0].keys() if rows else [])}


def unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def angle_deg(u, v):
    return math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(unit(u), unit(v)))))))


def kink_deg(w, out_dir):
    """angle de branchement SIGNE de la tangente sortante a l'aile, positif dans
    le sens trigonometrique. Pour cette fissure (bas-gauche -> haut-droite) et
    sigma_1 vertical, le cote en traction est le sens trigonometrique aux DEUX
    pointes : MLER attend +70,5 deg ; une valeur negative = mauvais cote."""
    cr = out_dir[0] * w[1] - out_dir[1] * w[0]
    return math.degrees(math.atan2(cr, float(np.dot(out_dir, w))))


def wing_cluster(J, tip, h):
    """amas connexe de joints rompus (non preexistants) issu d'une pointe :
    deux joints sont relies s'ils partagent une extremite (a 0,2 h pres)"""
    P1 = np.c_[J["x1"], J["y1"]]
    P2 = np.c_[J["x2"], J["y2"]]
    tol = 0.2 * h
    near = lambda p, q: np.linalg.norm(p - q, axis=-1) < tol
    seeds = np.where(near(P1, tip) | near(P2, tip))[0]
    seen, stack = set(seeds.tolist()), list(seeds)
    while stack:
        i = stack.pop()
        for q in (P1[i], P2[i]):
            nb = np.where(near(P1, q) | near(P2, q))[0]
            for j in nb:
                if j not in seen:
                    seen.add(j)
                    stack.append(j)
    return sorted(seen)


def analyse_run(name, var):
    od = outdir(name)
    lg = C.parse_log(od)
    hist = C.read_history(od)
    t, sig = hist["t"], hist["sigma"]
    J = read_csv(os.path.join(od, "fdem_final_joints.csv"))
    h = var["h"] * 1e-3
    mid = np.c_[0.5 * (J["x1"] + J["x2"]), 0.5 * (J["y1"] + J["y2"])]
    tb = J["tBreak"]
    # preexistants : tBreak = 0 (§5.11) ; rompus en cours de calcul : tBreak > 0
    pre = (tb == 0.0) & (J["bonded"] == 0)
    brk = tb > 0.0
    res = dict(h_mm=var["h"], v_fermeture=var["v"], T=var["T"], wall=lg["wall"], dt=lg["dt"],
               steps=lg["steps"], done=lg["done"], n_pre=int(pre.sum()), n_rompus=int(brk.sum()),
               sigma_max_MPa=float(sig.max() / 1e6))
    if brk.any():
        i0 = int(np.argmin(np.where(brk, tb, np.inf)))
        res["premiere_rupture"] = dict(t=float(tb[i0]), x=float(mid[i0, 0]), y=float(mid[i0, 1]),
                                       sigma_MPa=float(np.interp(tb[i0], t, sig) / 1e6))
    tips = {}
    for lab, tip, out_dir in (("A_bas", TIP_A, unit(TIP_A - TIP_B)), ("B_haut", TIP_B, unit(TIP_B - TIP_A))):
        d = np.linalg.norm(mid - tip, axis=1)
        sel = np.where(brk & (d < R_TIP * h))[0]
        r = dict(amorce=bool(len(sel)))
        if len(sel):
            k = sel[np.argmin(tb[sel])]
            p1, p2 = np.array([J["x1"][k], J["y1"][k]]), np.array([J["x2"][k], J["y2"][k]])
            far = p1 if np.linalg.norm(p1 - tip) > np.linalg.norm(p2 - tip) else p2
            r.update(t=float(tb[k]), sigma_MPa=float(np.interp(tb[k], t, sig) / 1e6),
                     mode=int(J["breakMode"][k]) if "breakMode" in J else None,
                     theta_deg=kink_deg(far - tip, out_dir),
                     angle_sigma1_deg=angle_deg(far - tip, np.array([0.0, 1.0 if far[1] >= tip[1] else -1.0])),
                     d_pointe_mm=float(d[k] * 1e3))
            # trajet en fin de calcul
            sub = {kk: vv[brk] for kk, vv in J.items()}
            cl = wing_cluster(sub, tip, h)
            if cl:
                pts = np.r_[np.c_[sub["x1"][cl], sub["y1"][cl]], np.c_[sub["x2"][cl], sub["y2"][cl]]]
                dd = np.linalg.norm(pts - tip, axis=1)
                fp = pts[np.argmax(dd)]
                r.update(n_joints_aile=len(cl), corde_mm=float(dd.max() * 1e3),
                         corde_angle_sigma1_deg=angle_deg(fp - tip, np.array([0.0, 1.0 if fp[1] >= tip[1] else -1.0])),
                         corde_theta_deg=kink_deg(fp - tip, out_dir))
        tips[lab] = r
    res["pointes"] = tips
    # indicateur secondaire : entree en endommagement (catalogue microsismique)
    sp = os.path.join(od, "fdem_seismic.csv")
    if os.path.exists(sp):
        S = read_csv(sp)
        if S:
            ds = np.minimum(np.hypot(S["x"] - TIP_A[0], S["y"] - TIP_A[1]),
                            np.hypot(S["x"] - TIP_B[0], S["y"] - TIP_B[1]))
            m = ds < R_TIP * h
            if m.any():
                ty = float(S["tYield"][m].min())
                res["endommagement_pointe"] = dict(t=ty, sigma_MPa=float(np.interp(ty, t, sig) / 1e6))
    am = [r for r in tips.values() if r["amorce"]]
    if am:
        f = min(am, key=lambda r: r["t"])
        res["sigma1c_MPa"] = f["sigma_MPa"]
        res["t_amorcage"] = f["t"]
        res["theta_deg"] = float(np.mean([r["theta_deg"] for r in am]))
    return res, (t, sig), J, (pre, brk)


def verdict(runs):
    v = {}
    ref = runs.get("h0p7")
    if ref and "sigma1c_MPa" in ref:
        s = ref["sigma1c_MPa"]
        v["P3.3_fourchette_MLER"] = CRIT["lefm_min"] <= s <= CRIT["lefm_max"]
        v["B2_ecart_YGeo"] = abs(s / CRIT["ygeo"] - 1) <= CRIT["ygeo_tol"]
        th = [r["theta_deg"] for r in ref["pointes"].values() if r["amorce"]]
        v["angle"] = bool(th) and all(min(abs(x - tr) for tr in CRIT["theta_refs"]) <= CRIT["theta_tol"] for x in th)
        v["deux_pointes"] = all(r["amorce"] for r in ref["pointes"].values())
        for n in ("h1p0", "h0p5"):
            if n in runs and "sigma1c_MPa" in runs[n]:
                v[f"maillage_{n}"] = abs(runs[n]["sigma1c_MPa"] / s - 1) <= CRIT["mesh_tol"]
    a, b = runs.get("h1p0"), runs.get("h1p0_v0p05")
    if a and b and "sigma1c_MPa" in a and "sigma1c_MPa" in b:
        v["vitesse"] = abs(b["sigma1c_MPa"] / a["sigma1c_MPa"] - 1) <= CRIT["rate_tol"]
    return v


def figures(data):
    plt = C.plot_style()
    lisjak = json.load(open(os.path.join(HERE, "lisjak_reference.json")))
    pub = lisjak["table_3_3"]
    # 1. contrainte axiale et amorcage
    fig, ax = plt.subplots(figsize=(6.3, 3.8))
    ax.axhspan(CRIT["lefm_min"], CRIT["lefm_max"], color="0.9", label="fourchette MLER (Table 3.3)")
    ax.axhline(pub["ygeo_fdem"], color="k", ls="--", lw=0.8, label="Y-Geo, 7,0 MPa")
    for i, (n, (res, (t, s), J, _)) in enumerate(data.items()):
        col = f"C{i}"
        ax.plot(t * 1e3, s / 1e6, color=col, lw=0.9, label=n.replace("p", ","))
        if "sigma1c_MPa" in res:
            ax.plot(res["t_amorcage"] * 1e3, res["sigma1c_MPa"], "o", color=col, ms=5)
    ax.set_xlabel("temps (ms)")
    ax.set_ylabel(r"contrainte axiale $\sigma_1$ (MPa)")
    ax.legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout()
    C.savefig(fig, os.path.join(HERE, "fig_sigma_amorcage"))
    plt.close(fig)
    # 2. faciès des ailes, zoom au centre, par variante
    n = len(data)
    fig, axs = plt.subplots(1, n, figsize=(2.4 * n + 0.6, 3.4), squeeze=False)
    for ax, (name, (res, _, J, (pre, brk))) in zip(axs[0], data.items()):
        tb = J["tBreak"]
        for k in np.where(pre)[0]:
            ax.plot([J["x1"][k] * 1e3, J["x2"][k] * 1e3], [J["y1"][k] * 1e3, J["y2"][k] * 1e3], "k-", lw=1.6)
        if brk.any():
            tmax = tb[brk].max()
            for k in np.where(brk)[0]:
                c = plt.cm.viridis(tb[k] / tmax)
                ax.plot([J["x1"][k] * 1e3, J["x2"][k] * 1e3], [J["y1"][k] * 1e3, J["y2"][k] * 1e3], "-", color=c, lw=1.0)
        ax.set_xlim(XC * 1e3 - 10, XC * 1e3 + 10)
        ax.set_ylim(YC * 1e3 - 12, YC * 1e3 + 12)
        ax.set_aspect("equal")
        ax.set_title(name.replace("p", ","), fontsize=8)
        ax.set_xlabel("x (mm)")
    axs[0][0].set_ylabel("y (mm)")
    fig.suptitle("noir : fissure ; couleur : instant de rupture", fontsize=8)
    fig.tight_layout()
    C.savefig(fig, os.path.join(HERE, "fig_ailes"))
    plt.close(fig)


def analyse(names=None, extra=None):
    allv = dict(VARIANTS)
    allv.update(SMOKE)
    names = names or list(allv)
    data, runs = {}, {}
    for n in names:
        if not C.is_done(outdir(n)):
            continue
        res, ts, J, masks = analyse_run(n, allv[n])
        data[n], runs[n] = (res, ts, J, masks), res
    if not os.path.exists(os.path.join(HERE, "lisjak_reference.json")):
        json.dump(LE.reference(), open(os.path.join(HERE, "lisjak_reference.json"), "w"), indent=1)
    v = verdict(runs)
    json.dump(dict(criteres=CRIT, variantes=runs, verdict=v,
                   passe=bool(v) and all(v.values())),
              open(os.path.join(HERE, "resultats.json"), "w"), indent=1, ensure_ascii=False)
    print(f"\n{'variante':12s} {'h':>4s} {'pas':>8s} {'mur s':>7s} {'rompus':>6s} {'s1c MPa':>8s} "
          f"{'t_am ms':>8s} {'theta A':>8s} {'theta B':>8s} {'s_endo':>7s}")
    for n, r in runs.items():
        th = {k: (f"{p['theta_deg']:.1f}" if p["amorce"] else "-") for k, p in r["pointes"].items()}
        print(f"{n:12s} {r['h_mm']:4g} {r['steps'] or 0:8d} {r['wall'] or 0:7.1f} {r['n_rompus']:6d} "
              f"{r.get('sigma1c_MPa', float('nan')):8.2f} {1e3 * r.get('t_amorcage', float('nan')):8.3f} "
              f"{th['A_bas']:>8s} {th['B_haut']:>8s} "
              f"{r.get('endommagement_pointe', {}).get('sigma_MPa', float('nan')):7.2f}")
    print("verdict :", v if v else "(variante de reference absente)")
    if data:
        figures(data)
    return runs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["prepare", "run", "analyse", "all"])
    ap.add_argument("--exe", default=C.EXE_DEFAULT)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--variants", nargs="+", default=None)
    ap.add_argument("--smoke", action="store_true", help="variante d'essai h = 2,5 mm, T = 2,5 ms")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    table = dict(SMOKE) if a.smoke else dict(VARIANTS)
    if a.variants:
        allv = dict(VARIANTS); allv.update(SMOKE)
        table = {n: allv[n] for n in a.variants}
    if a.action in ("prepare", "run", "all"):
        jobs = [make_job(n, v) for n, v in table.items()]
        print(f"[P3 aile] {len(jobs)} decks ecrits dans {OUT}")
        if a.action in ("run", "all"):
            for j in jobs:
                rc, wall = C.run_job(j, os.path.abspath(a.exe), a.threads, a.force)
                print(f"  {j.name:14s} rc = {rc}  {wall:7.1f} s", flush=True)
    if a.action in ("analyse", "all"):
        analyse(list(table))


if __name__ == "__main__":
    main()
