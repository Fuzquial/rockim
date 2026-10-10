# -*- coding: utf-8 -*-
"""analyse_frac.py — depouillement des bancs 2, 3, 4 (fissuration d'une barre / plaque SENT).

usage : python analyse_frac.py <dossier_banc> [W H h_par_nom]
Lit, pour chaque <nom>.cfg : <nom>.log (resume), <nom>.dur, <nom>_out/history.csv, la PREMIERE
et la DERNIERE trame fdem_joints_*.vtu (positions initiales ; damage, tInsert, tBreak, bonded).
Topologie : les copies de noeud sont identifiees par leurs coordonnees initiales (arrondi 1e-9 m).
- rompu : tBreak >= 0 ; insere (adaptatif) : bonded == 0 ; endommage : damage > 0.
- N_fiss : composantes connexes des joints rompus (partage d'un sommet initial).
- fissure principale : la plus grande composante (nombre de joints).
- hors fissure : joints inseres (adaptatif) ou endommages (intrinseque) dont le milieu est a plus
  de 2h du milieu du joint le plus proche de la fissure principale.
- propagation (adaptation de tunnel_edz/tools/nucleation_vs_propagation.py, instants exacts au lieu
  des trames) : un joint rompu est « propagation » s'il partage un sommet avec un joint rompu
  strictement avant lui ; le tout premier est exclu du rapport. Meme chose a l'insertion (tInsert).
Ecrit <banc>_resultats.csv et <banc>_fissures.npz (segments initiaux des joints rompus).
"""
import glob, io, math, os, re, sys
import numpy as np

def arr(s, name):
    m = re.search(r'Name="%s"[^>]*>\s*(.*?)\s*</DataArray>' % name, s, re.S)
    return np.array(m.group(1).split(), float) if m else None

def jvtu(path):
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    P = np.array(re.search(r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>', s, re.S)
                 .group(1).split(), float).reshape(-1, 3)[:, :2]
    con = arr(s, "connectivity").astype(int).reshape(-1, 2)
    return P, con, {k: arr(s, k) for k in ("damage", "tInsert", "tBreak", "bonded")}

def summary(log):
    s = open(log, errors="ignore").read()
    def g(r):
        m = re.search(r, s)
        return float(m.group(1)) if m else float("nan")
    return dict(Wj=g(r"joints\s+: ([-0-9.e+]+) J/m cohesif"), Wc=g(r"contact\s+: ([-0-9.e+]+) J/m"),
                gcnet=g(r"net work injected by general contact: ([-0-9.e+]+)"),
                resid=g(r"residu\s+: ([-0-9.e+]+) J/m"), wall=g(r"wall time: ([0-9.e+-]+)"),
                dt=g(r"dt = ([0-9.e+-]+) s"), nj=g(r"elements, ([0-9]+) joints"),
                nel=g(r"\] ([0-9]+) elements,"), ins=g(r"adaptive insertion: ([0-9]+) /"))

def comps(n1, n2, idx):
    """composantes connexes des joints idx (union-find sur les sommets)"""
    par = {}
    def f(x):
        while par.setdefault(x, x) != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for j in idx:
        a, b = f(n1[j]), f(n2[j])
        if a != b: par[a] = b
    lab = np.array([f(n1[j]) for j in idx])
    return lab

def prop_share(t, n1, n2, mask):
    idx = np.nonzero(mask)[0]
    if len(idx) < 2: return float("nan"), 0
    order = idx[np.argsort(t[idx], kind="stable")]
    first_t = {}
    for j in order:
        for n in (n1[j], n2[j]):
            first_t.setdefault(n, t[j])
    prop = 0
    for j in order[1:]:
        if min(first_t[n1[j]], first_t[n2[j]]) < t[j] - 1e-12: prop += 1
    return prop / (len(order) - 1), len(order)

def analyse(run, name, W, h, notch=None):
    fs = sorted(glob.glob(run + "/fdem_joints_[0-9]*.vtu"))
    P0, con, _ = jvtu(fs[0]); _, _, F = jvtu(fs[-1])
    _, node = np.unique(np.round(P0 / 1e-9).astype(np.int64), axis=0, return_inverse=True)
    node = node.ravel()
    n1, n2 = node[con[:, 0]], node[con[:, 1]]
    A, B = P0[con[:, 0]], P0[con[:, 1]]
    mid = 0.5 * (A + B); Lj = np.linalg.norm(B - A, axis=1)
    D, tI, tB, bd = F["damage"], F["tInsert"], F["tBreak"], F["bonded"]
    broken = tB >= 0
    adaptive = "adapt" in name or name.startswith("b4")
    act = (bd < 0.5) if adaptive else (D > 0)
    r = dict(name=name, njoints=len(con), n_ins_or_dmg=int(act.sum()), n_broken=int(broken.sum()))
    ib = np.nonzero(broken)[0]
    if len(ib):
        lab = comps(n1, n2, ib)
        u, cnt = np.unique(lab, return_counts=True)
        r["N_fiss"] = len(u)
        r["N_fiss_ge3"] = int((cnt >= 3).sum())
        # fusion des composantes distantes de moins de 1,5 h (fissure qui saute un element :
        # le motif de rupture est le meme, la topologie par sommets partages la coupe en deux)
        grp = {k: k for k in u}
        def g(x):
            while grp[x] != x: x = grp[x]
            return x
        pts = {k: np.r_[A[ib[lab == k]], B[ib[lab == k]]] for k in u}
        for i1, k1 in enumerate(u):
            for k2 in u[i1 + 1:]:
                if np.min(np.linalg.norm(pts[k1][:, None] - pts[k2][None], axis=2)) < 1.5 * h:
                    a1, a2 = g(k1), g(k2)
                    if a1 != a2: grp[a1] = a2
        lab = np.array([g(k) for k in lab])
        u, cnt = np.unique(lab, return_counts=True)
        r["N_fiss_fus"] = len(u)
        main = ib[lab == u[np.argmax(cnt)]]
        r["main_n"] = len(main); r["main_frac"] = len(main) / len(ib)
        xs = np.r_[A[main, 0], B[main, 0]]
        r["main_xext_W"] = (xs.max() - xs.min()) / W
        r["main_L_mm"] = Lj[main].sum() * 1e3
        r["main_ymean_mm"] = mid[main, 1].mean() * 1e3
        r["main_ystd_mm"] = mid[main, 1].std() * 1e3
        ia = np.nonzero(act)[0]
        d = np.min(np.linalg.norm(mid[ia, None, :] - mid[None, main, :], axis=2), axis=1)
        r["off_frac"] = float((d > 2 * h).mean())
        r["prop_break"], _ = prop_share(tB, n1, n2, broken)
        if adaptive:
            r["prop_insert"], _ = prop_share(tI, n1, n2, act)
        else:
            r["prop_insert"] = float("nan")
        if notch is not None:
            tip = np.array(notch)
            j0 = ib[np.argmin(tB[ib])]
            r["first_dist_tip_mm"] = np.linalg.norm(mid[j0] - tip) * 1e3
            dtip = np.min(np.linalg.norm(np.r_[A[main], B[main]] - tip, axis=1))
            r["main_dist_tip_mm"] = dtip * 1e3
            r["main_xmax_mm"] = xs.max() * 1e3
            dev = np.abs(mid[main, 1] - tip[1])
            r["dev_mean_mm"] = dev.mean() * 1e3; r["dev_max_mm"] = dev.max() * 1e3
            r["L_lig"] = r["main_L_mm"] / ((W - tip[0]) * 1e3)
            r["broken_off_main"] = 1 - r["main_frac"]
    seg = np.c_[A[ib], B[ib]]
    segm = np.c_[A[main], B[main]] if len(ib) else np.zeros((0, 4))
    return r, seg, segm

def main():
    d = os.path.abspath(sys.argv[1]); banc = os.path.basename(d)
    rows, npz = [], {}
    for cfg in sorted(glob.glob(d + "/*.cfg")):
        name = os.path.basename(cfg)[:-4]; run = cfg[:-4] + "_out"
        cache = cfg[:-4] + "_depouille.npz"      # cache : les .vtu sont supprimes apres depouillement
        if not glob.glob(run + "/fdem_joints_[0-9]*.vtu"):
            if os.path.exists(cache):
                z = np.load(cache, allow_pickle=True)
                rows.append(z["row"].item())
                for k in ("seg", "main", "t", "sig"):
                    npz[name + "__" + k] = z[k]
            continue
        m = re.search(r"_h([0-9.]+)_", name); h = float(m.group(1)) * 1e-3
        W = 0.04 if banc == "b3" else 0.02
        notch = (0.01, 0.04) if banc == "b3" else None
        r, seg, segm = analyse(run, name, W, h, notch)
        S = summary(cfg[:-4] + ".log"); r.update(S)
        hh = np.genfromtxt(run + "/history.csv", delimiter=",", names=True)
        r["sig_pk_MPa"] = np.nanmax(hh["sigma"]) / 1e6
        r["sig_pk_ft"] = r["sig_pk_MPa"] / 10.0
        r["sig_end_ft"] = float(hh["sigma"][-1]) / 1e7      # contrainte residuelle en fin de calcul
        r["W_GW"] = S["Wj"] / (40.0 * (0.03 if banc == "b3" else W) * 1.0)   # b3 : ligament 30 mm
        dur = open(cfg[:-4] + ".dur").read() if os.path.exists(cfg[:-4] + ".dur") else ""
        r["dur_s"] = float(re.search(r"dur=([0-9.]+)", dur).group(1)) if "dur=" in dur else float("nan")
        r["rc"] = re.search(r"rc=(\d+)", dur).group(1) if "rc=" in dur else "?"
        rows.append(r)
        npz[name + "__seg"] = seg; npz[name + "__main"] = segm
        npz[name + "__t"] = hh["t"]; npz[name + "__sig"] = hh["sigma"]
        if os.path.exists(cfg[:-4] + ".dur"):   # cache seulement pour un calcul termine
            np.savez_compressed(cache, row=np.array(r, dtype=object), seg=seg, main=segm, t=hh["t"], sig=hh["sigma"])
    keys = []
    for r in rows:
        for k in r:
            if k not in keys: keys.append(k)
    with open(os.path.join(d, banc + "_resultats.csv"), "w") as f:
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join(("%.6g" % r[k]) if isinstance(r.get(k), float) else str(r.get(k, "")) for k in keys) + "\n")
    np.savez_compressed(os.path.join(d, banc + "_fissures.npz"), **npz)
    show = ["sig_pk_ft", "sig_end_ft", "W_GW", "n_ins_or_dmg", "n_broken", "N_fiss", "N_fiss_fus", "N_fiss_ge3", "main_frac", "main_xext_W",
            "off_frac", "prop_break", "prop_insert", "Wc", "gcnet", "dur_s"]
    if banc == "b3":
        show += ["first_dist_tip_mm", "main_dist_tip_mm", "main_xmax_mm", "dev_mean_mm", "dev_max_mm", "L_lig"]
    for r in rows:
        print(r["name"], " ".join("%s=%s" % (k, ("%.4g" % r[k]) if isinstance(r.get(k), (float, int)) else r.get(k)) for k in show))

if __name__ == "__main__":
    main()
