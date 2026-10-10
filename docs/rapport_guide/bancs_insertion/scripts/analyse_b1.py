# -*- coding: utf-8 -*-
"""analyse_b1.py — depouille le banc 1 (joint isole). Ecrit b1_resultats.csv et b1_courbes.npz.

Ouverture du joint cible : aux deux extremites du joint, difference des positions moyennes des
copies de noeud des elements au-dessus et au-dessous de la ligne du joint (trames fdem_*.vtu),
projetee sur n (normale, vers le haut) et t (tangente). Contrainte : sigma = -gripFy/(W t) de
history.csv (signee, traction > 0), projetee sur le joint (champ uniaxial uniforme) :
sigma_n = sigma cos^2 theta, tau = sigma sin theta cos theta.
"""
import glob, io, math, os, re, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
B1 = os.path.join(ROOT, sys.argv[1] if len(sys.argv) > 1 else "b1")
W, a, t = 0.010, 0.040, 1.0
FT, C, GI = 10e6, 25e6, 70.0
GII = 10 * GI

def pts(path):
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    P = np.array(re.search(r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>', s, re.S)
                 .group(1).split(), float).reshape(-1, 3)[:, :2]
    return P

def summary(log):
    s = open(log).read()
    g = lambda r: float(re.search(r, s).group(1)) if re.search(r, s) else float("nan")
    return dict(Wj=g(r"joints\s+: ([-0-9.e+]+) J/m cohesif"),
                wall=g(r"wall time: ([0-9.e+-]+)"),
                dt=g(r"dt = ([0-9.e+-]+) s"))

def analyse(name):
    run = os.path.join(B1, name + "_out")
    mode = name.split("_")[1]
    th = 0.0 if mode == "I" else 45.0  # I : horizontal ; II, mix, mixc : 45 deg
    thr = math.radians(th)
    n = np.array([-math.sin(thr), math.cos(thr)]); tv = np.array([math.cos(thr), math.sin(thr)])
    ends = [np.array([0.0, a]), np.array([W, a + W * math.tan(thr)])]
    fr = sorted(glob.glob(run + "/fdem_[0-9]*.vtu"))
    P0 = pts(fr[0])
    cen = P0.reshape(-1, 3, 2).mean(1)
    # au-dessus de la ligne du joint ?
    above_el = ((cen - ends[0]) @ n) > 0
    above = np.repeat(above_el, 3)
    sel = []
    for e in ends:
        m = np.linalg.norm(P0 - e, axis=1) < 1e-9
        sel.append((m & above, m & ~above))
    tf = np.loadtxt(run + "/frames.csv", delimiter=",", skiprows=1)[:, 1]
    dn, ds = [], []
    for f in fr:
        P = pts(f)
        d = np.mean([P[u].mean(0) - P[l].mean(0) - (P0[u].mean(0) - P0[l].mean(0)) for u, l in sel], axis=0)
        dn.append(d @ n); ds.append(d @ tv)
    dn, ds = np.array(dn), np.array(ds)
    h = np.genfromtxt(run + "/history.csv", delimiter=",", names=True)
    sig = -h["gripFy"] / (W * t)
    sn_h = sig * math.cos(thr) ** 2; tau_h = sig * math.sin(thr) * math.cos(thr)
    sn = np.interp(tf, h["t"], sn_h); tau = np.interp(tf, h["t"], tau_h)
    S = summary(os.path.join(B1, name + ".log"))
    L = W / math.cos(thr)
    r = dict(name=name, mode=mode, ins=name.split("_")[2], pen=name.split("_")[3][1:],
             soft=name.split("_")[4], wall=S["wall"], dt=S["dt"])
    r["sn_pk_ft"] = sn_h.max() / FT if mode != "II" else float("nan")
    r["tau_pk_c"] = np.abs(tau_h).max() / C if mode != "I" else float("nan")
    # sigma_n au joint a l instant du pic de |sigma|
    r["sn_at_pk_ft"] = sn_h[np.argmax(np.abs(sig))] / FT
    r["W"] = S["Wj"]
    r["W_GI"] = S["Wj"] / (GI * L * t)
    r["W_GII"] = S["Wj"] / (GII * L * t)
    # ouverture de rupture : premier instant ou |contrainte| < 1 % du pic apres le pic
    sref = sn if mode == "I" else (np.abs(tau) if mode == "II" else sn)
    ipk = int(np.argmax(sref))
    post = np.nonzero(sref[ipk:] < 0.01 * sref[ipk])[0]
    r["dn_rupt_um"] = dn[ipk + post[0]] * 1e6 if len(post) else float("nan")
    r["ds_rupt_um"] = abs(ds[ipk + post[0]]) * 1e6 if len(post) else float("nan")
    r["dn_pic_um"] = dn[ipk] * 1e6; r["ds_pic_um"] = abs(ds[ipk]) * 1e6
    return r, dict(tf=tf, dn=dn, ds=ds, sn=sn, tau=tau)

def main():
    names = sorted(os.path.basename(f)[:-4] for f in glob.glob(B1 + "/b1*_*.cfg"))
    rows, curves = [], {}
    for nm in names:
        r, cv = analyse(nm)
        rows.append(r)
        for k, v in cv.items():
            curves[nm + "__" + k] = v
    keys = list(rows[0].keys())
    with open(os.path.join(B1, os.path.basename(B1) + "_resultats.csv"), "w") as f:
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join(("%.6g" % r[k]) if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
    np.savez_compressed(os.path.join(B1, os.path.basename(B1) + "_courbes.npz"), **curves)
    for r in rows:
        print("%-24s pk(sn/ft)=%.4f pk(tau/c)=%.4f sn@pk=%.4f W=%.5f W/GIL=%.4f W/GIIL=%.4f "
              "dn_pk=%.2f dn_r=%.2f ds_pk=%.2f ds_r=%.2f um wall=%.2f"
              % (r["name"], r["sn_pk_ft"], r["tau_pk_c"], r["sn_at_pk_ft"], r["W"], r["W_GI"],
                 r["W_GII"], r["dn_pic_um"], r["dn_rupt_um"], r["ds_pic_um"], r["ds_rupt_um"], r["wall"]))

if __name__ == "__main__":
    main()
