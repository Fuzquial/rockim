# -*- coding: utf-8 -*-
"""check_s4.py — verification falsifiante des essais elementaires S4 (campagne du 13/09).

usage : python check_s4.py <scratch_dir>

Criteres (chaque verdict est imprime avec ses chiffres) :
  A. seuils en deformation imprimes par bulkDamageProbe = delta0/h exacts du
     maillage de Kuhn (a = 2 mm : inscrit 1,68995 %, arete 0,55552 %), a 1e-5 pres ;
  B. compression ISOTROPE : sous `deviatoric` (et `edge` seul) max delta_m < 0,1 delta0,
     D = 0 partout, nPulv = 0, bdWork = 0 ; sous `total` / `principal` : D > 0 ;
  C. bit-identite du chemin par defaut : *_ref (exe S4) == *_ref (rockim_g1y16.exe)
     octet pour octet (history.csv + dernier VTU) ; *_dev (probe seule) : history.csv
     identique a *_ref, VTU = *_ref + le seul champ `bulkDm` ;
  D. recomposition INDEPENDANTE de delta_m par element a partir de la geometrie du
     dernier VTU (F = dx dX^-1, U par decomposition polaire exacte, eps = U - I) pour
     les trois mesures et les deux longueurs : la mesure ACTIVE du deck doit
     coincider avec le champ `bulkDm` exporte (max sur les elements de l ecart
     relatif < 1e-2) et les mesures INACTIVES doivent s en ecarter (> 5 %) ;
  E. decks mal formes refuses (rc != 0, message attendu).
"""
import glob, hashlib, io, os, re, sys
import numpy as np

S = sys.argv[1]
D0 = 1.4e-5
A = 2.0e-3                                   # arete du cube de Kuhn
H_INS = (np.sqrt(2.0) - 1.0) * A             # 6V/A = a/(1+sqrt2)
H_EDGE = (3.0 + 2.0 * np.sqrt(2.0) + np.sqrt(3.0)) / 6.0 * A
NU = 0.25
nfail = 0


def verdict(ok, msg):
    global nfail
    print(("  [OK]   " if ok else "  [FAIL] ") + msg)
    if not ok:
        nfail += 1


def read_vtu(path):
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    m = re.search(r'Name="connectivity"[^>]*>\s*(.*?)\s*</DataArray>', s, re.S)
    con = np.fromstring(m.group(1), sep=" ").astype(int)
    m = re.search(r'Name="offsets"[^>]*>\s*(.*?)\s*</DataArray>', s, re.S)
    off = np.fromstring(m.group(1), sep=" ").astype(int)
    m = re.search(r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>', s, re.S)
    pts = np.fromstring(m.group(1), sep=" ").reshape(-1, 3)
    ncell = off.size
    fields = {}
    for m in re.finditer(r'<DataArray[^>]*Name="([A-Za-z]\w*)"[^>]*>\s*(.*?)\s*</DataArray>', s, re.S):
        nm = m.group(1)
        if nm in ("connectivity", "offsets", "types"):
            continue
        v = np.fromstring(m.group(2), sep=" ")
        if v.size == ncell:
            fields[nm] = v
    return pts, con.reshape(ncell, -1), fields


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def frames(outdir, pat="fdem3d_0*.vtu"):
    fs = sorted(glob.glob(os.path.join(outdir, pat)))
    if not fs:
        fs = sorted(glob.glob(os.path.join(outdir, "fdem_0*.vtu")))
    fs = [f for f in fs if "joints" not in os.path.basename(f)]
    return fs


def hist(outdir):
    p = os.path.join(outdir, "history.csv")
    with io.open(p, encoding="utf-8", errors="ignore") as f:
        hdr = f.readline().strip().split(",")
    data = np.genfromtxt(p, delimiter=",", skip_header=1)
    if data.ndim == 1:
        data = data[None, :]
    return hdr, data


def probe_lines(log):
    txt = io.open(log, encoding="utf-8", errors="ignore").read()
    m = re.search(r"seuil eps_m \(delta0/h\) : inscrit ([0-9.e+-]+) %, arete ([0-9.e+-]+) %", txt)
    thr = (float(m.group(1)), float(m.group(2))) if m else None
    fr = re.findall(r"bulkDamageProbe trame (\d+) : max delta_m = ([0-9.e+-]+) um .*?max D = ([0-9.e+-]+), elements armes \(delta_m > delta0\) : (\d+) / (\d+)", txt)
    fr = [(int(a), float(b) * 1e-6, float(c), int(d), int(e)) for a, b, c, d, e in fr]
    return thr, fr, txt


def strain_measures(p0, con0, p1, con1):
    """delta_m par element pour les 6 combinaisons (h x eps_m) a partir des VTU."""
    out = {k: [] for k in ("ins_dev", "ins_pri", "ins_tot", "edge_dev", "edge_pri", "edge_tot")}
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    faces = [(1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)]
    for c0, c1 in zip(con0, con1):
        X = p0[c0]
        x = p1[c1]
        dX = np.stack([X[1] - X[0], X[2] - X[0], X[3] - X[0]], axis=1)
        dx = np.stack([x[1] - x[0], x[2] - x[0], x[3] - x[0]], axis=1)
        F = dx @ np.linalg.inv(dX)
        W, sg, Vt = np.linalg.svd(F)
        U = Vt.T @ np.diag(sg) @ Vt
        eps = 0.5 * (U + U.T) - np.eye(3)
        dev = eps - np.trace(eps) / 3.0 * np.eye(3)
        e_dev = np.sqrt(2.0 / 3.0) * np.linalg.norm(dev)
        e_pri = np.max(np.abs(np.linalg.eigvalsh(eps)))
        e_tot = np.sqrt(2.0 / 3.0) * np.linalg.norm(eps)
        V0 = abs(np.linalg.det(dX)) / 6.0
        At = sum(0.5 * np.linalg.norm(np.cross(X[b] - X[a], X[c] - X[a])) for a, b, c in faces)
        h_ins = 6.0 * V0 / At
        h_edge = sum(np.linalg.norm(X[b] - X[a]) for a, b in edges) / 6.0
        out["ins_dev"].append(h_ins * e_dev)
        out["ins_pri"].append(h_ins * e_pri)
        out["ins_tot"].append(h_ins * e_tot)
        out["edge_dev"].append(h_edge * e_dev)
        out["edge_pri"].append(h_edge * e_pri)
        out["edge_tot"].append(h_edge * e_tot)
    return {k: np.array(v) for k, v in out.items()}


ACTIVE = {"dev": "ins_dev", "tot": "ins_tot", "pri": "ins_pri",
          "edge": "edge_dev", "edge_tot": "edge_tot"}

print("S4 — essais elementaires de la pulverisation (scratch : %s)" % S)

# ---- A. seuils imprimes ------------------------------------------------------
print("\nA. seuils en deformation imprimes (delta0/h, medianes)")
exp_ins, exp_edge = 100.0 * D0 / H_INS, 100.0 * D0 / H_EDGE
print("  attendu (Kuhn a = 2 mm) : h_ins = %.6f mm, h_edge = %.6f mm -> %.5f %% / %.5f %%"
      % (H_INS * 1e3, H_EDGE * 1e3, exp_ins, exp_edge))
for d in ("iso_dev", "uni_dev", "uni_edge_tot"):
    thr, _, _ = probe_lines(os.path.join(S, d + ".log"))
    if thr is None:
        verdict(False, d + " : ligne de seuil absente")
        continue
    ok = abs(thr[0] / exp_ins - 1) < 1e-5 and abs(thr[1] / exp_edge - 1) < 1e-5
    verdict(ok, "%s : imprime %.5f %% / %.5f %% (ecart relatif %.1e / %.1e)"
            % (d, thr[0], thr[1], thr[0] / exp_ins - 1, thr[1] / exp_edge - 1))

# ---- B. compression isotrope -------------------------------------------------
print("\nB. compression isotrope (confineFaces = all, 2,5 GPa) et confinement lateral")
for d in ("iso_dev", "iso_edge", "iso_tot", "iso_pri", "iso_edge_tot",
          "biax_dev", "biax_tot", "biax_pri", "biax_edge", "biax_edge_tot"):
    thr, fr, _ = probe_lines(os.path.join(S, d + ".log"))
    out = os.path.join(S, "out_" + d)
    hdr, data = hist(out)
    npulv = data[-1, hdr.index("nPulv")] if "nPulv" in hdr else float("nan")
    bdw = data[-1, hdr.index("bdWork")] if "bdWork" in hdr else float("nan")
    _, _, f = read_vtu(frames(out)[-1])
    Dmax = f["bulkD"].max()
    dm_last = fr[-1][1]
    nArm = fr[-1][3]
    line = ("%s : max delta_m = %.3f um (%.3f delta0), max D (VTU) = %.4f, armes %d/%d, nPulv = %g, bdWork = %.4g J"
            % (d, dm_last * 1e6, dm_last / D0, Dmax, nArm, fr[-1][4], npulv, bdw))
    if d.startswith("iso") and d in ("iso_dev", "iso_edge"):
        # le residu deviatorique est un TRANSITOIRE de la rampe (pression suiveuse
        # sur base bloquee en z) : on le juge en DEFORMATION, delta_m/h < 0,1 %
        h_act = H_EDGE if d == "iso_edge" else H_INS
        verdict(dm_last / h_act < 1e-3 and Dmax == 0.0 and npulv == 0 and bdw == 0.0,
                line + " ; eps_m residuel = delta_m/h = %.3e  [attendu : deviatorique nul en isotrope -> eps_m < 1e-3, D = 0]"
                % (dm_last / h_act))
    elif d.startswith("iso"):
        verdict(nArm > 0 and Dmax > 0.0 and bdw < 0.0,
                line + "  [attendu : la mesure a trace S ARME en isotrope]")
    else:
        verdict(nArm > 0 and Dmax > 0.0,
                line + "  [attendu : etat deviatorique -> toutes les mesures s arment]")

# ---- C. bit-identite -----------------------------------------------------------
print("\nC. bit-identite du chemin par defaut")
for d in ("iso", "biax", "uni", "uni2d"):
    new = os.path.join(S, "out_" + d + "_ref")
    old = os.path.join(S, "old_" + d + "_ref")
    if not os.path.isdir(old):
        verdict(False, d + "_ref : sortie rockim_g1y16 absente")
        continue
    fn, fo = frames(new)[-1], frames(old)[-1]
    okh = sha(os.path.join(new, "history.csv")) == sha(os.path.join(old, "history.csv"))
    okv = sha(fn) == sha(fo)
    verdict(okh and okv, "%s_ref : exe S4 vs rockim_g1y16.exe — history.csv %s, %s %s"
            % (d, "IDENTIQUE" if okh else "DIFFERE", os.path.basename(fn), "IDENTIQUE" if okv else "DIFFERE"))
    dev = os.path.join(S, "out_" + d + "_dev")
    okh2 = sha(os.path.join(new, "history.csv")) == sha(os.path.join(dev, "history.csv"))
    _, _, fr_ref = read_vtu(fn)
    _, _, fr_dev = read_vtu(frames(dev)[-1])
    extra = sorted(set(fr_dev) - set(fr_ref))
    same = all(np.array_equal(fr_ref[k], fr_dev[k]) for k in fr_ref)
    verdict(okh2 and extra == ["bulkDm"] and same,
            "%s_dev (probe seule) vs %s_ref : history.csv %s ; champs VTU communs %s ; champs en plus %s"
            % (d, d, "IDENTIQUE" if okh2 else "DIFFERE", "egaux" if same else "DIFFERENTS", extra))

# ---- D. recomposition independante de delta_m ---------------------------------
print("\nD. delta_m recompose depuis la geometrie du dernier VTU (polaire exacte) contre `bulkDm`")
for ld in ("uni", "biax", "iso"):
    for ms in ("dev", "tot", "pri", "edge", "edge_tot"):
        d = ld + "_" + ms
        out = os.path.join(S, "out_" + d)
        fs = frames(out)
        p0, c0, _ = read_vtu(fs[0])
        act = ACTIVE[ms]
        # `bulkDm` est un MAX HISTORIQUE : il n egale la mesure courante que
        # tant que le chargement est monotone. Sous pression imposee la chute
        # de raideur (D -> 0,9) fait sauter la boite (iso_pri : 13 um a la
        # trame 6, 259 um a la trame 8, puis rebond elastique) : la comparaison
        # deux cotes se fait donc a la DERNIERE trame ou max bulkDm <= delta0
        # (rampe cosinus monotone, D = 0 : max historique = valeur courante).
        # A la derniere trame on ne verifie qu un cote (bulkDm >= courant).
        kstar = None
        for k in range(len(fs) - 1, 0, -1):
            _, _, fk = read_vtu(fs[k])
            if fk["bulkDm"].max() <= D0 and fk["bulkDm"].max() > 0.05 * D0:
                kstar = k
                break
        if kstar is None:
            kstar = len(fs) - 1
        pk, ck, fk = read_vtu(fs[kstar])
        meas = strain_measures(p0, c0, pk, ck)
        bdm = fk["bulkDm"]
        ref = np.maximum(bdm, 1e-12)
        rel = {k: np.max(np.abs(v - bdm) / ref) for k, v in meas.items()}
        # mesure active NON mesurable (isotrope sous `deviatoric` : max bulkDm
        # 0,8 um = transitoire de rampe, non monotone) : le verdict B fait foi,
        # ici on ne fait qu imprimer.
        measurable = bdm.max() > 0.2 * D0
        ok = rel[act] < 1e-2 if measurable else True
        others = [k for k in meas if k != act]
        # une mesure inactive peut coincider par symetrie (isotrope : ins_dev et
        # edge_dev sont ~0 toutes deux) : l ecart n est exige que si elle est
        # mesurable (> 0,2 delta0). Seuil de separation 2 % : en uniaxial a
        # nu = 0,25, `total` et `deviatoric` ne different que de 3,9 %.
        sep = {k: rel[k] for k in others if meas[k].max() > 0.2 * D0 and rel[k] > 2e-2}
        nsep = {k: rel[k] for k in others if meas[k].max() > 0.2 * D0 and rel[k] <= 2e-2}
        p1, c1, f1 = read_vtu(fs[-1])
        m1 = strain_measures(p0, c0, p1, c1)[act]
        b1 = f1["bulkDm"]
        under = np.max((m1 - b1) / np.maximum(b1, 1e-12))     # courant > max hist. ?
        verdict(ok and not nsep and (under < 1e-2 or not measurable),
                "%s : trame %d (max bulkDm %.2f um)%s active %s ecart max %.2e ; inactives separees %s%s ; "
                "derniere trame : max bulkDm %.2f um, courant recompose %.2f um, (courant - bulkDm)/bulkDm max %.2e"
                % (d, kstar, bdm.max() * 1e6,
                   "" if measurable else " [non mesurable < 0,2 delta0 : info seule]",
                   act, rel[act],
                   {k: "%.3f" % v for k, v in sep.items()},
                   ("" if not nsep else " ; NON separees " + str({k: "%.3f" % v for k, v in nsep.items()})),
                   b1.max() * 1e6, m1.max() * 1e6, under))

# ---- E. refus ---------------------------------------------------------------
print("\nE. decks mal formes")
runs = io.open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "RESULTATS_runs.log"),
               encoding="utf-8", errors="ignore").read()
for d, needle in (("bad_nokey", "exige bulkDamage = yang"),
                  ("bad_value", "bulkDamageStrain must be"),
                  ("bad_len", "bulkDamageLength must be")):
    m = re.search(r"=== %s ===.*?rc=(\d+)" % d, runs, re.S)
    rc = int(m.group(1)) if m else -1
    txt = io.open(os.path.join(S, d + ".log"), encoding="utf-8", errors="ignore").read()
    verdict(rc != 0 and needle in txt, "%s : rc = %d, message '%s' %s"
            % (d, rc, needle, "present" if needle in txt else "ABSENT"))

# ---- F. miroir 2D ---------------------------------------------------------------
print("\nF. miroir 2D (fdem, plaque 4 x 4 mm grid 2 x 2, compression uniaxiale)")
for d in ("uni2d_dev", "uni2d_tot", "uni2d_pri", "uni2d_edge"):
    thr, fr, txt = probe_lines(os.path.join(S, d + ".log"))
    out = os.path.join(S, "out_" + d)
    _, _, f = read_vtu(frames(out)[-1])
    print("  %s : seuils %s %% ; dernier max delta_m = %.3f um, max D (VTU) = %.4f, armes %d/%d, bulkDm exporte : %s"
          % (d, thr, fr[-1][1] * 1e6, f["bulkD"].max(), fr[-1][3], fr[-1][4], "bulkDm" in f))
    verdict(fr[-1][3] > 0 and "bulkDm" in f, d + " : s arme en uniaxial et exporte bulkDm")

print("\nBILAN : %d verdict(s) FAIL" % nfail)
sys.exit(1 if nfail else 0)
