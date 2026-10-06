#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# B_impact_coupe : partie III de la campagne V&V (reproduction d'articles),
# bancs d'IMPACT et de COUPE.
#
#   B8  impact : Saksala (IJNAMG 35, 2011), bouton unique pousse par une
#       impulsion de contrainte dans la tige (200 MPa, 110 us), boutons
#       cylindrique et hemispherique, granite-marbre de la Table I. Code contre
#       code : la loi de Saksala est portee telle quelle dans rockim
#       (law = saksala2011, defauts = Table I). Modele : fem3d, quart de bloc,
#       contact de Signorini, source a impedance (toolPulse*).
#   B9  coupe  : Heilman et al. (ARMA 24-0238, 2024, HOSS), cutter PDC 3D sur
#       granite Utah FORGE, garde 20 deg, 10 m/s, passes 0,508 / 1,016 /
#       1,524 mm. Code contre code (force de coupe, fig. 3B). Reference
#       secondaire, experimentale : LeBaron et al. (GRC Trans. 47, 2023),
#       raclage du granitoide FORGE aux memes passes et a la meme garde
#       (force moyenne, MSE). Modele : fdem3d, scenario = shear,
#       toolShape = pdc, toolContact = signorini.
#
#   python3 vv/B_impact_coupe/b_impact_coupe.py all            # campagne complete
#   python3 vv/B_impact_coupe/b_impact_coupe.py all --smoke    # essai de fumee
#   python3 vv/run_queue.py B_impact_coupe --slots 2 --threads 2
#
# Regles : cles rockim existantes seulement, src/ intouche, criteres ecrits
# ci-dessous AVANT le premier calcul (2026-10-04).
# ---------------------------------------------------------------------------
import argparse, json, math, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
VV = os.path.dirname(HERE)
ROOT = os.path.dirname(VV)
sys.path.insert(0, VV)
import vvcommon as C  # noqa: E402

REF = json.load(open(os.path.join(HERE, "reference_publiee.json")))
MESHES = os.path.join(HERE, "meshes")

# ===========================================================================
# CRITERES D'ACCEPTATION (fixes le 2026-10-04, avant tout calcul)
# ---------------------------------------------------------------------------
# B8. Ecart relatif a la valeur publiee. Budget d'ecart : lecture +-3 a 6 %,
# axisymetrique CST (Saksala) contre tetraedres 3D en quart, maillage que
# Saksala declare lui-meme non objectif (fig. 13 : courbes « considerablement »
# differentes entre deux maillages), celerite de la tige et forme exacte du
# trapeze non publiees, rayon du bouton hemispherique ambigu.
CRIT_B8 = dict(F_max=0.15, u_max=0.20, u_res=0.25)
# B9 (Heilman, code contre code). Pic de force filtre (fenetre 0,05 mm de
# course) a +-25 % : lecture +-0,5 kN, maille 0,5 mm contre 0,25 mm, chanfrein
# non cote, amortissement de HOSS non publie. Position du pic a +-0,5 mm.
# Ordre des pics F(0,508) < F(1,016) < F(1,524) comme sur la fig. 3B.
# Signature : la force retombe sous 20 % du pic avant 3,0 mm de course.
CRIT_B9 = dict(F_pic=0.25, x_pic=0.5e-3, ordre=True, chute_frac=0.20, chute_x=3.0e-3)
# B9 (LeBaron, essai) : force moyenne sur la course [0,2 ; 2,0] mm dans la bande
# publiee +-1 ecart-type. Aucun ajustement en cas d'echec : la vitesse de
# l'essai (25 mm/s) est 400 fois plus faible que celle du calcul (10 m/s).
CRIT_B9_LEB = dict(bande_sd=1.0, fenetre=(0.2e-3, 2.0e-3))

# ===========================================================================
# B8 : Saksala 2011
# ---------------------------------------------------------------------------
SK = REF["saksala2011"]
SK_W, SK_H = 0.10, 0.15                    # quart : rayon de domaine 100 mm (fig. 5a), hauteur 150 mm
SK_ROD_A = SK["chargement"]["A_rod"]        # 166,7 mm2 = 1 000 mm2 / 6 boutons
SK_SIG = SK["chargement"]["sigma_A"]        # 200 MPa
SK_CB = SK["chargement"]["c_bar"]["valeur"]  # 5 189 m/s (hypothese acier)
SK_RHOB = SK["chargement"]["rho_bit"]
SK_MBIT = SK["chargement"]["m_bit"]
SK_TABLE = "0:0 1e-05:1 0.0001:1 0.00011:0"  # t_rise 10 us, palier jusqu'a 100 us, descente 10 us
SK_T = 2.2e-4                               # couvre la decharge (u_res stabilise)
SK_CASES = {                                # nom -> (toolShape, toolRadius)
    "sk11_cyl": ("flat", 5.0e-3),           # poincon plat de diametre 10 mm
    "sk11_hemi_R10": ("sphere", 10.0e-3),   # cote R = 10 mm de la fig. 5d
    "sk11_hemi_R5": ("sphere", 5.0e-3),     # sensibilite : vraie demi-sphere de diametre 10 mm
}
SK_MESH_FULL = dict(hFine=0.6e-3, rFine=8e-3, hFar=6e-3, rFar=60e-3)
SK_MESH_SMOKE = dict(hFine=2.5e-3, rFine=6e-3, hFar=15e-3, rFar=50e-3)


def sk_mesh(m):
    tag = f"sk11_q_h{m['hFine'] * 1e3:g}".replace(".", "p")
    p = os.path.join(MESHES, tag + ".msh")
    if os.path.exists(p):
        return p
    import gmsh
    os.makedirs(MESHES, exist_ok=True)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("sk11")
    v = gmsh.model.occ.addBox(0, 0, 0, SK_W, SK_W, SK_H)
    gmsh.model.occ.synchronize()
    gmsh.model.addPhysicalGroup(3, [v], name="rock")
    # sommet du coin (0, 0, H) : point d'impact du quart, deja dans la geometrie
    # (pas de point libre, donc pas de noeud orphelin sous le bouton, cf. §5.20)
    pt = [t for d, t in gmsh.model.getEntitiesInBoundingBox(-1e-6, -1e-6, SK_H - 1e-6,
                                                            1e-6, 1e-6, SK_H + 1e-6, 0)]
    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "PointsList", pt)
    ft = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(ft, "InField", fd)
    gmsh.model.mesh.field.setNumber(ft, "SizeMin", m["hFine"])
    gmsh.model.mesh.field.setNumber(ft, "SizeMax", m["hFar"])
    gmsh.model.mesh.field.setNumber(ft, "DistMin", m["rFine"])
    gmsh.model.mesh.field.setNumber(ft, "DistMax", m["rFar"])
    gmsh.model.mesh.field.setAsBackgroundMesh(ft)
    for k in ("MeshSizeExtendFromBoundary", "MeshSizeFromPoints", "MeshSizeFromCurvature"):
        gmsh.option.setNumber("Mesh." + k, 0)
    gmsh.option.setNumber("Mesh.RandomSeed", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    gmsh.write(p)
    gmsh.finalize()
    return p


def sk_deck(name, outdir, smoke):
    shape, R = SK_CASES[name]
    mesh = sk_mesh(SK_MESH_SMOKE if smoke else SK_MESH_FULL)
    F_inc = SK_ROD_A * SK_SIG / 4.0                 # quart
    Z = SK_RHOB * SK_CB * SK_ROD_A / 4.0
    lines = [f"# B8 Saksala 2011, {name} (ecrit par b_impact_coupe.py){' FUMEE' if smoke else ''}",
             "mode = fem3d", "law = saksala2011", "scenario = percussion",
             "mesh = file", f"meshFile = {mesh}", "quarterModel = true",
             f"T = {SK_T:g}", f"frames = {2 if smoke else 6}",
             # Table I ; les parametres sk* gardent leurs defauts = Table I (§5.5)
             "rho = 2600", "E = 60e9", "nu = 0.2", "ft = 13e6", "cohesion = 37.5e6",
             "frictionDeg = 30",
             f"toolShape = {shape}", f"toolRadius = {R:g}",
             "toolX = 0", "toolY = 0", "toolGap = 1e-6",
             f"toolMass = {SK_MBIT / 4.0:g}", "impactSpeed = 0",
             f"toolPulseForce = {F_inc:.6g}", f"toolPulseTable = {SK_TABLE}",
             f"toolPulseImpedance = {Z:.6g}",
             "toolContact = signorini", "toolSignoriniRelax = 0",
             "contactMu = 0",                      # Saksala : seul le ddl vertical du bouton
             "absorbing = all", "absorbSpringFactor = 0",
             "dampingLocal = 0"]
    return C.write_deck(outdir, lines)


# ===========================================================================
# B9 : Heilman 2024 (+ LeBaron 2023)
# ---------------------------------------------------------------------------
HE = REF["heilman2024"]
CUT_W, CUT_D, CUT_H = 0.040, 0.030, 0.020
CUT_NOTCH, CUT_JEU = 0.003, 0.000284        # face verticale de l'entaille a x = 3 mm
CUT_TOOLX0 = 0.0025                          # arete 0,5 mm avant la face
CUT_V = 10.0
CUT_STROKE = 3.0e-3                          # course utile (fig. 3B : 0 a 3 mm)
CUT_DEPTHS = [0.508e-3, 1.016e-3, 1.524e-3]
CUT_MESH_FULL = dict(hFine=0.5e-3, hFar=4e-3)
CUT_MESH_SMOKE = dict(hFine=1.2e-3, hFar=6e-3)


def cut_tag(d, br):
    return f"he24_br{br:g}_d{d * 1e3:.3f}".replace(".", "p")


def cut_mesh(d, m):
    tag = f"cut3d_d{d * 1e6:.0f}_h{m['hFine'] * 1e3:g}".replace(".", "p")
    p = os.path.join(MESHES, tag + ".msh")
    if not os.path.exists(p):
        os.makedirs(MESHES, exist_ok=True)
        args = [CUT_W, CUT_D, CUT_H, d, CUT_JEU, CUT_NOTCH, m["hFine"], 0.001, 0.010,
                0.006, 0.003, m["hFar"]]
        subprocess.run([sys.executable, os.path.join(ROOT, "tools", "make_cut3d_mesh.py")]
                       + [f"{a:g}" for a in args] + [p, "1"], check=True,
                       stdout=subprocess.DEVNULL)
    return p


def cut_deck(d, br, outdir, smoke):
    mesh = cut_mesh(d, CUT_MESH_SMOKE if smoke else CUT_MESH_FULL)
    T = 3.0e-5 if smoke else (CUT_NOTCH - CUT_TOOLX0 + CUT_STROKE) / CUT_V
    lines = [f"# B9 Heilman 2024, garde {br} deg, passe {d * 1e3:.3f} mm, {CUT_V} m/s "
             f"(ecrit par b_impact_coupe.py){' FUMEE' if smoke else ''}",
             "# deck derive de configs/cut3d_heilman.cfg (memes cles, memes choix commentes) ;",
             "# differences voulues : dampingLocal = 0 (regle des bancs de force, §5.6 ter, T0b) et",
             "# jointSecantRatchet = on (origin seul est non conservatif, avertissement du solveur, §5.4)",
             "mode = fdem3d", "scenario = shear", f"T = {T:.6g}",
             f"frames = {2 if smoke else 30}", "historyFlush = true",
             "mesh = file", f"meshFile = {mesh}",
             "rho = 2610", "E = 48.26e9", "nu = 0.23", "ft = 10.62e6", "cohesion = 44.815e6",
             "frictionDeg = 38.66", "Gf = 152.9", "gfShearFactor = 4.220", "crushCap = 1e12",
             "toolShape = pdc", "cutterDia = 0.013", "cutterThick = 0.0025", "chamferLen = 0.0",
             f"backRakeDeg = {-br:g}", "cutterFloor = true",
             f"cutDepth = {d:g}", f"cutSpeed = {CUT_V:g}",
             f"toolX = {CUT_TOOLX0:g}", f"toolY = {CUT_D / 2:g}", "contactMu = 0.8",
             "toolContact = signorini", "toolSignoriniRelax = 0.0",
             "insertion = adaptive", "jointSoftening = yan", "jointShearUnload = origin",
             "jointFrictionScaled = 1", "jointSecantRatchet = on", "contact = potential", "gcActivation = adaptive",
             "insertionPenaltyFactor = 4", "dtFactor = 0.15", "dampingLocal = 0",
             "jointXi = 0.05", "budgetAbortPct = 5", "budgetAbortMin = 1.0"]
    return C.write_deck(outdir, lines)


# ===========================================================================
# file de runs
# ---------------------------------------------------------------------------
def cases(smoke=False, etendu=False):
    """(banc, nom, fabrique de deck, poids)"""
    out = []
    sk = ["sk11_cyl", "sk11_hemi_R10"] if smoke else list(SK_CASES)
    for n in sk:
        out.append(("B8", n, lambda od, n=n: sk_deck(n, od, smoke), 1.0))
    brs = [20, 10] if etendu else [20]
    deps = [CUT_DEPTHS[1]] if smoke else CUT_DEPTHS
    for br in brs:
        for d in deps:
            out.append(("B9", cut_tag(d, br), lambda od, d=d, br=br: cut_deck(d, br, od, smoke),
                        20.0 * d / 1e-3))
    return out


def outroot(smoke):
    return os.path.join(HERE, "out_smoke" if smoke else "out")


def JOBS(args=None):
    smoke = bool(getattr(args, "smoke", False))
    etendu = bool(getattr(args, "etendu", False))
    jobs = []
    for bench, name, mk, w in cases(smoke, etendu):
        od = os.path.join(outroot(smoke), name)
        os.makedirs(od, exist_ok=True)
        jobs.append(C.Job(bench, name, mk(od), od, weight=w))
    return jobs


# ===========================================================================
# depouillement
# ---------------------------------------------------------------------------
def movavg_t(t, y, win):
    """moyenne glissante sur une fenetre de duree win (echantillonnage quelconque)"""
    c = np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])
    out = np.empty_like(y)
    for i, ti in enumerate(t):
        a, b = ti - win / 2, ti + win / 2
        ia, ib = np.searchsorted(t, [a, b])
        ia, ib = max(ia, 0), min(ib, len(t) - 1)
        out[i] = (c[ib] - c[ia]) / (t[ib] - t[ia]) if t[ib] > t[ia] else y[i]
    return out


def analyse_b8(smoke):
    res = {}
    for name in SK_CASES:
        od = os.path.join(outroot(smoke), name)
        if not C.is_done(od):
            continue
        lg = C.parse_log(od)
        h = C.read_history(od)
        t, z, fz = h["t"], h["toolZ"], h["toolFz"]
        u = (z[0] - z) - 1e-6                     # enfoncement sous la surface initiale
        F = 4.0 * np.abs(fz)                      # quart -> bouton entier
        Ff = movavg_t(t, F, 2e-6)
        i = int(np.argmax(Ff))
        key = "cyl" if "cyl" in name else "hemi"
        r = dict(F_max=float(Ff[i]), F_max_brut=float(F.max()), u_at_Fmax=float(u[i]),
                 u_max=float(u.max()), u_fin=float(u[-1]), t_fin=float(t[-1]),
                 wall=lg["wall"], dt=lg["dt"], ntet=lg["ntet"])
        # enfoncement residuel : u quand la force est retombee sous 2 % du pic apres le pic
        after = np.where((np.arange(len(F)) > i) & (Ff < 0.02 * Ff[i]))[0]
        r["u_res"] = float(u[after[0]]) if len(after) else float("nan")
        ref = SK["resultats"][key]
        ecarts, ok = {}, True
        for q in ("F_max", "u_max", "u_res"):
            e = r[q] / ref[q]["valeur"] - 1 if np.isfinite(r[q]) else float("nan")
            ecarts[q] = e
            ok &= bool(np.isfinite(e) and abs(e) <= CRIT_B8[q])
        r["ecarts"], r["passe"] = ecarts, ok and not smoke
        res[name] = r
        print(f"[B8] {name:15s} Fmax {r['F_max'] / 1e3:7.2f} kN ({100 * ecarts['F_max']:+6.1f} %)  "
              f"umax {r['u_max'] * 1e3:6.3f} mm ({100 * ecarts['u_max']:+6.1f} %)  "
              f"ures {r['u_res'] * 1e3:6.3f} mm  t_fin {r['t_fin'] * 1e6:.0f} us  "
              f"{'FUMEE' if smoke else ('PASSE' if r['passe'] else 'ECHEC')}")
    return res


def course(x, F):
    """course de l'arete depuis le premier contact outil-roche (force non nulle).
    Avec la garde de 20 deg la face de coupe penche en avant (b < 0, §5.6 quater) :
    elle touche l'arete superieure de l'entaille 0,37 mm (passe 1,016) avant que
    l'arete n'atteigne la face verticale ; Heilman trace la force depuis le debut
    de la coupe, d'ou l'origine au premier contact."""
    nz = np.where(F > 0)[0]
    return x - (x[nz[0]] if len(nz) else CUT_NOTCH)


def analyse_b9(smoke):
    res = {}
    leb = REF["lebaron2023"]["resultats"]
    for br in (20, 10):
        for d in CUT_DEPTHS:
            name = cut_tag(d, br)
            od = os.path.join(outroot(smoke), name)
            if not C.is_done(od):
                continue
            lg = C.parse_log(od)
            h = C.read_history(od)
            t, x, fx, W = h["t"], h["toolX"], h["toolFx"], h["work"]
            Fraw = np.abs(fx)
            s = course(x, Fraw)                       # course depuis le premier contact
            # force par le travail : F = dW/ds sur une fenetre de 0,05 mm de course,
            # insensible au pic d'impulsion du premier toucher (Signorini : r/dt)
            Wc = np.abs(W)
            Fw = np.full_like(s, np.nan)
            for k in range(len(s)):
                ia = np.searchsorted(s, s[k] - 0.025e-3)
                ib = min(np.searchsorted(s, s[k] + 0.025e-3), len(s) - 1)
                if s[ib] > s[ia]:
                    Fw[k] = (Wc[ib] - Wc[ia]) / (s[ib] - s[ia])
            m = (s > 0.05e-3) & np.isfinite(Fw)
            r = dict(course_max=float(s.max()), wall=lg["wall"], dt=lg["dt"], ntet=lg["ntet"],
                     broken=lg["broken"], budget_pct=lg["budget_pct"],
                     F_brut_max=float(Fraw.max()))
            if m.any():
                i = np.where(m)[0][np.argmax(Fw[m])]
                r.update(F_pic=float(Fw[i]), x_pic=float(s[i]))
                after = np.where(m & (s > s[i]) & (Fw < CRIT_B9["chute_frac"] * Fw[i]))[0]
                r["x_chute"] = float(s[after[0]]) if len(after) else float("nan")
                fen = m & (s >= CRIT_B9_LEB["fenetre"][0]) & (s <= CRIT_B9_LEB["fenetre"][1])
                r["F_moy"] = float(np.nanmean(Fw[fen])) if fen.any() else float("nan")
            key = f"d{d * 1e3:.3f}"
            refs = HE["resultats"].get(f"br{br}_v10", {}).get(key)
            if refs and "F_pic" in r:
                r["ecart_F_pic"] = r["F_pic"] / refs["F_pic"] - 1
                r["ecart_x_pic"] = r["x_pic"] - refs["x_pic"]
                r["passe_heilman"] = bool(abs(r["ecart_F_pic"]) <= CRIT_B9["F_pic"]
                                          and abs(r["ecart_x_pic"]) <= CRIT_B9["x_pic"]
                                          and np.isfinite(r["x_chute"])
                                          and r["x_chute"] <= CRIT_B9["chute_x"]) and not smoke
            if br == 20 and key in leb and np.isfinite(r.get("F_moy", np.nan)):
                L = leb[key]
                r["ecart_F_moy_leb"] = r["F_moy"] / L["F_moy"] - 1
                r["passe_lebaron"] = bool(abs(r["F_moy"] - L["F_moy"])
                                          <= CRIT_B9_LEB["bande_sd"] * L["F_sd"]) and not smoke
            res[name] = r
            print(f"[B9] {name:20s} course {r['course_max'] * 1e3:5.2f} mm  "
                  f"F_pic {r.get('F_pic', float('nan')) / 1e3:7.2f} kN a {r.get('x_pic', float('nan')) * 1e3:5.2f} mm  "
                  f"F_moy {r.get('F_moy', float('nan')) / 1e3:6.2f} kN  brut max {r['F_brut_max'] / 1e3:7.2f} kN  "
                  f"rompus {r['broken']}  B4 {r['budget_pct']}  "
                  f"{'FUMEE' if smoke else ''}")
    # ordre des pics (fig. 3B)
    p = [res.get(cut_tag(d, 20), {}).get("F_pic") for d in CUT_DEPTHS]
    ordre = None if any(v is None for v in p) else bool(p[0] < p[1] < p[2])
    return res, ordre


def figures(r8, r9, smoke):
    plt = C.plot_style()
    sfx = "_smoke" if smoke else ""
    # B8 : force-enfoncement
    fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.8))
    for ax, key, names in ((axs[0], "cyl", ["sk11_cyl"]), (axs[1], "hemi", ["sk11_hemi_R10", "sk11_hemi_R5"])):
        ref = SK["resultats"][key]
        ax.plot([ref["u_max"]["valeur"] * 1e3], [ref["F_max"]["valeur"] / 1e3], "ks", ms=5,
                label="Saksala 2011 (pic)")
        ax.plot([ref["u_res"]["valeur"] * 1e3], [0], "k^", ms=5, label="Saksala 2011 (résiduel, lu)")
        for k, n in enumerate(names):
            od = os.path.join(outroot(smoke), n)
            if C.is_done(od):
                h = C.read_history(od)
                ax.plot(((h["toolZ"][0] - h["toolZ"]) - 1e-6) * 1e3, 4 * np.abs(h["toolFz"]) / 1e3,
                        color=f"C{k}", lw=0.9, label=f"rockim {n.replace('sk11_', '')}")
        ax.set_xlabel("enfoncement du bouton (mm)")
        ax.set_ylabel("force (kN)")
        ax.set_title("bouton cylindrique" if key == "cyl" else "bouton hémisphérique", fontsize=9)
        ax.legend(fontsize=6.5, frameon=False)
    fig.tight_layout()
    C.savefig(fig, os.path.join(HERE, "fig_b8_force_enfoncement" + sfx))
    # B9 : force de coupe-course
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    for k, d in enumerate(CUT_DEPTHS):
        ref = HE["resultats"]["br20_v10"][f"d{d * 1e3:.3f}"]
        ax.plot([ref["x_pic"] * 1e3], [ref["F_pic"] / 1e3], "s", color=f"C{k}", ms=5)
        od = os.path.join(outroot(smoke), cut_tag(d, 20))
        if C.is_done(od):
            h = C.read_history(od)
            ax.plot(course(h["toolX"], np.abs(h["toolFx"])) * 1e3, np.abs(h["toolFx"]) / 1e3, color=f"C{k}",
                    lw=0.7, label=f"rockim, passe {d * 1e3:.3f} mm")
    ax.set_xlim(-0.5, 3.0)
    ax.set_xlabel("course depuis le premier contact (mm)")
    ax.set_ylabel("force de coupe (kN)")
    ax.set_title("carrés : pics de Heilman et al. 2024, fig. 3B (lus)", fontsize=8)
    ax.legend(fontsize=6.5, frameon=False)
    fig.tight_layout()
    C.savefig(fig, os.path.join(HERE, "fig_b9_force_coupe" + sfx))


def analyse(smoke=False):
    r8 = analyse_b8(smoke)
    r9, ordre = analyse_b9(smoke)
    out = dict(date="2026-10-04", smoke=smoke,
               criteres=dict(B8=CRIT_B8, B9=CRIT_B9, B9_lebaron=CRIT_B9_LEB),
               B8=r8, B9=r9, B9_ordre_des_pics=ordre)
    fn = "resultats_smoke.json" if smoke else "resultats.json"
    json.dump(out, open(os.path.join(HERE, fn), "w"), indent=1, default=float)
    figures(r8, r9, smoke)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["prepare", "run", "analyse", "all"])
    ap.add_argument("--exe", default=C.EXE_DEFAULT)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--smoke", action="store_true", help="maillages grossiers, T court")
    ap.add_argument("--etendu", action="store_true", help="ajoute la garde 10 deg (fig. 2B)")
    ap.add_argument("--only", nargs="+", default=None, help="restreint aux runs nommes")
    a = ap.parse_args()
    if a.action in ("prepare", "run", "all"):
        jobs = JOBS(a)
        if a.only:
            jobs = [j for j in jobs if j.name in a.only]
        for j in jobs:
            print(f"[prepare] {j.bench} {j.name} -> {j.cfg}")
        if a.action in ("run", "all"):
            for j in jobs:
                rc, wall = C.run_job(j, os.path.abspath(a.exe), a.threads)
                print(f"[run] {j.bench} {j.name:22s} rc = {rc}  {wall:7.1f} s", flush=True)
    if a.action in ("analyse", "all"):
        analyse(a.smoke)


if __name__ == "__main__":
    main()
