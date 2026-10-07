# -*- coding: utf-8 -*-
"""check_s3.py - verifie le mini-test S3 (banc de joint cinematique) contre son schema FALSIFIANT.

Lit les journaux [JOINTBENCH] des neuf decks, les deux refus, les deux ancres bitid, et recalcule
independamment l aire sous sigma(dn) depuis jointbench.csv (trapezes numpy, meme donnee, autre
code) pour recouper le chiffre imprime par le solveur.

usage : python check_s3.py <scratch_dir>   (defaut : le dossier s3 du scratchpad)
Ecrit RESULTATS_check_s3.txt a cote de ce script. Code de sortie 0 si tout le schema tient.
"""
import io, os, re, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
S = sys.argv[1] if len(sys.argv) > 1 else (
    "C:/Users/FUZQUI~1/AppData/Local/Temp/claude/C--Users-fuzquianoalricabi-simulations/"
    "01c1f07a-4d93-418a-851e-98d14203ea97/scratchpad/s3")

DECKS = ["ten_solidity", "ten_plastic", "ten_origin", "sh_solidity", "sh_plastic", "sh_origin",
         "tilt_sol", "tilt_pl_maj", "tilt_pl_any"]
OUT = []


def say(s=""):
    print(s)
    OUT.append(s)


def num(s):
    return float(s)


def parse(log):
    """Extrait les chiffres des lignes [JOINTBENCH] du resume."""
    d = {"ok": False}
    txt = io.open(log, encoding="utf-8", errors="replace").read()
    d["wall"] = num(m.group(1)) if (m := re.search(r"wall time: ([0-9.eE+-]+)", txt)) else None
    d["ok"] = d["wall"] is not None
    m = re.search(r"\(iii\) aire sous sigma\(dn\).*?W_I = ([0-9.eE+-]+) J/m2 = ([0-9.eE+-]+) Gf ; attendu par la loi codee ([0-9.eE+-]+) J/m2 = ([0-9.eE+-]+) Gf(?: -> ecart ([0-9.eE+-]+) % \[(\w+))?", txt)
    if m:
        d["wI"], d["wI_Gf"], d["wIexp"], d["wIexp_Gf"] = map(num, m.groups()[:4])
        d["wI_err"] = num(m.group(5)) if m.group(5) else None
        d["wI_verdict"] = m.group(6)
    m = re.search(r"\(iii\) aire sous tau\(ds\).*?W_II = ([0-9.eE+-]+) J/m2 = ([0-9.eE+-]+) GfII ; attendu a pression normale nulle ([0-9.eE+-]+) J/m2 = ([0-9.eE+-]+) GfII", txt)
    if m:
        d["wII"], d["wII_GfII"], d["wIIexp"], d["wIIexp_GfII"] = map(num, m.groups())
    m = re.search(r"\(i\) retrace : max \|(\w+)_recharge - \w+_charge\| = ([0-9.eE+-]+) Pa = ([0-9.eE+-]+) (\w+) \((\d+) echantillons apparies, (\d+) sans jumeau\) ; decharge : ([0-9.eE+-]+) Pa = ([0-9.eE+-]+) \w+ \((\d+) apparies\) -> \[(\w+)", txt)
    if m:
        d["i_var"] = m.group(1); d["i_devR"] = num(m.group(2)); d["i_devR_rel"] = num(m.group(3))
        d["i_ref"] = m.group(4); d["i_nR"] = int(m.group(5)); d["i_nRnm"] = int(m.group(6))
        d["i_devU"] = num(m.group(7)); d["i_devU_rel"] = num(m.group(8)); d["i_nU"] = int(m.group(9))
        d["i_verdict"] = m.group(10)
    m = re.search(r"\(ii\) (\w+) residuel\(le\) a traction nulle sur la decharge, par point \[m\] :((?: \S+){3}) ; glissement plastique interne \|slip\| :((?: \S+){3})", txt)
    if m:
        d["ii_var"] = m.group(1)
        # "<=x" = pas de changement de signe jusqu au retournement : residuel borne par |x|min
        d["ii_res"] = [num(v[2:]) if v.startswith("<=") else (None if v.startswith("non") else num(v))
                       for v in m.group(2).split()]
        d["ii_slip"] = [num(v) for v in m.group(3).split()]
    m = re.search(r"\(iv\) rupture des points \(D >= 1\) a t =((?: \S+){3}) s ; joint rompu \(tBreak\) a t = (\S+) s(?: \(nFail = (\d+) a cet instant\))?(?: -> (.*?) \[(\w[\w ]*)\])?", txt)
    if m:
        d["iv_tk"] = [None if v == "jamais" else num(v) for v in m.group(1).split()]
        d["iv_tB"] = None if m.group(2) == "jamais" else num(m.group(2))
        d["iv_nFail"] = int(m.group(3)) if m.group(3) else None
        d["iv_text"] = m.group(4); d["iv_verdict"] = m.group(5)
    m = re.search(r"dt = ([0-9.eE+-]+) s \(echantillonnage", txt)
    d["dt"] = num(m.group(1)) if m else None
    m = re.search(r"joint : pj = ([0-9.eE+-]+) Pa/m, ft = ([0-9.eE+-]+) Pa, c = ([0-9.eE+-]+) Pa, dnE = ([0-9.eE+-]+) m, dnF = ([0-9.eE+-]+) m, sE = ([0-9.eE+-]+) m, slipF = ([0-9.eE+-]+) m, Gf = ([0-9.eE+-]+) J/m2, GfII = ([0-9.eE+-]+)", txt)
    if m:
        d["pj"], d["ft"], d["c"], d["dnE"], d["dnF"], d["sE"], d["slipF"], d["Gf"], d["GfII"] = map(num, m.groups())
    return d


def area_from_csv(path, upto_break=True):
    """Aire sous sigma(dn) et tau(ds) par point, moyenne des trois, jusqu a la rupture (trapezes)."""
    a = np.genfromtxt(path, delimiter=",", names=True)
    n = len(a)
    if upto_break:
        ib = np.nonzero(a["broken"] > 0.5)[0]
        n = int(ib[0]) + 1 if len(ib) else n
    wI = wII = 0.0
    for k in range(3):
        dn, sg = a["dn%d" % k][:n], a["sig%d" % k][:n]
        ds, ta = a["ds%d" % k][:n], a["tau%d" % k][:n]
        wI += np.sum(0.5 * (sg[1:] + sg[:-1]) * np.diff(dn))
        wII += np.sum(0.5 * (ta[1:] + ta[:-1]) * np.diff(ds))
    # coherence broken -> nFail >= 2 (majority) : compte des lignes qui violent
    viol = int(np.sum((a["broken"] > 0.5) & (a["nFail"] < 2)))
    return wI / 3.0, wII / 3.0, len(a), viol


def main():
    fails = []
    say("S3 - banc de joint cinematique : verification du schema falsifiant")
    say("scratch : %s" % S)
    R = {}
    for d in DECKS:
        log = os.path.join(S, d + ".log")
        if not os.path.isfile(log):
            say("%-13s : journal absent" % d); fails.append(d + " absent"); continue
        r = parse(log)
        R[d] = r
        csv = os.path.join(S, "out_" + d, "jointbench.csv")
        if os.path.isfile(csv):
            r["csv_wI"], r["csv_wII"], r["csv_n"], r["csv_viol"] = area_from_csv(csv)
        say("")
        say("== %s (wall %.2f s, %s lignes csv, dt %s s)" % (d, r.get("wall") or -1, r.get("csv_n"), r.get("dt")))
        if "wI" in r:
            say("  (iii) W_I = %.4f J/m2 = %.5f Gf (attendu %.5f Gf) ecart %s %% [%s] ; recalcul csv %.4f J/m2 (ecart solveur/csv %.2e)"
                % (r["wI"], r["wI_Gf"], r["wIexp_Gf"], r.get("wI_err"), r.get("wI_verdict"),
                   r.get("csv_wI", float("nan")), abs(r["wI"] - r.get("csv_wI", r["wI"])) / max(r["wI"], 1e-30)))
            say("  (iii) W_II = %.4f J/m2 = %.5f GfII (attendu a sigma_n = 0 : %.5f GfII) ; recalcul csv %.4f"
                % (r["wII"], r["wII_GfII"], r["wIIexp_GfII"], r.get("csv_wII", float("nan"))))
        if "i_verdict" in r:
            say("  (i)   retrace %s : recharge %.3e %s (%d apparies, %d sans jumeau), decharge %.3e %s -> %s"
                % (r["i_var"], r["i_devR_rel"], r["i_ref"], r["i_nR"], r["i_nRnm"], r["i_devU_rel"], r["i_ref"], r["i_verdict"]))
        if "ii_res" in r:
            say("  (ii)  %s residuel a traction nulle [m] : %s ; |slip| interne : %s"
                % (r["ii_var"], ["%.3e" % v if v is not None else "non-atteint" for v in r["ii_res"]],
                   ["%.3e" % v for v in r["ii_slip"]]))
        if "iv_tB" in r:
            say("  (iv)  points rompus a t = %s s ; joint a t = %s s (nFail %s) -> %s [%s]"
                % (["%.4e" % v if v is not None else "jamais" for v in r["iv_tk"]],
                   "%.4e" % r["iv_tB"] if r["iv_tB"] is not None else "jamais", r.get("iv_nFail"), r.get("iv_text"), r.get("iv_verdict")))
        if r.get("csv_viol"):
            say("  !! %d lignes broken = 1 avec nFail < 2" % r["csv_viol"])

    say("")
    say("== schema falsifiant (ce que chaque deck DOIT donner)")

    def expect(cond, label):
        say("  [%s] %s" % ("OK " if cond else "KO ", label))
        if not cond:
            fails.append(label)

    def g(d, k, default=None):
        return R.get(d, {}).get(k, default)

    dlt_ten = 2.0e-5 / 4000.0
    dlt_sh = 3.0e-5 / 4000.0
    # (i) retrace : solidity PASS, plastic/origin FAIL
    expect(g("ten_solidity", "i_verdict") == "PASS", "ten_solidity (i) retrace PASS (loi sans memoire)")
    expect(g("ten_plastic", "i_verdict") == "FAIL", "ten_plastic  (i) retrace FAIL (secante de decharge eq. 17) - variante qui DOIT echouer")
    expect(g("ten_origin", "i_verdict") == "FAIL", "ten_origin   (i) retrace FAIL (mode I identique a plastic)")
    expect(g("sh_solidity", "i_verdict") == "PASS", "sh_solidity  (i) retrace PASS")
    expect(g("sh_plastic", "i_verdict") == "FAIL", "sh_plastic   (i) retrace FAIL (retour radial)")
    expect(g("sh_origin", "i_verdict") == "FAIL", "sh_origin    (i) retrace FAIL contre la charge (secante a l origine)")
    # (iii) aire : plastic/origin = Gf, solidity = 1,159 Gf, a 2 %
    for d, ratio in (("ten_solidity", 1.1593), ("ten_plastic", 1.0), ("ten_origin", 1.0)):
        w = g(d, "wI_Gf")
        expect(w is not None and abs(w / ratio - 1.0) <= 0.02 and g(d, "wI_verdict") == "PASS",
               "%-12s (iii) W_I/Gf = %s, attendu %.3f (+/- 2 %%)" % (d, w, ratio))
        cw = g(d, "csv_wI"); sw = g(d, "wI")
        expect(cw is not None and sw is not None and abs(cw - sw) <= 1e-6 * max(abs(sw), 1e-30),
               "%-12s (iii) recalcul independant depuis jointbench.csv = chiffre du solveur (1e-6)" % d)
    # (ii) residuel : plastic > 0, solidity/origin = 0
    res = g("sh_plastic", "ii_res") or []
    expect(len(res) == 3 and all(v is not None and v > 1.0e-6 for v in res),
           "sh_plastic   (ii) glissement residuel > 1 um sur les trois points (plasticite conservee)")
    for d in ("sh_solidity", "sh_origin"):
        res = g(d, "ii_res") or []
        expect(len(res) == 3 and all(v is not None and abs(v) <= 1e-3 * dlt_sh for v in res),
               "%-12s (ii) glissement residuel nul (< 1e-3 pas de chemin = %.1e m)" % (d, 1e-3 * dlt_sh))
    # (iv) regle 2/3 : majority -> mort au 2e point, instants distincts ; any -> au premier
    for d in ("tilt_sol", "tilt_pl_maj"):
        # sous majority le joint meurt au 2e point et n est plus evalue : le
        # 3e point reste "jamais" - c est precisement la regle nfail > 1
        tk = [v for v in (g(d, "iv_tk") or []) if v is not None]
        ok = len(tk) >= 2 and g(d, "iv_tB") is not None
        if ok:
            ts = sorted(tk); dt = g(d, "dt") or 1e-9
            ok = (ts[1] - ts[0] > 0.5 * dt) and abs(g(d, "iv_tB") - ts[1]) <= 0.5 * dt
        expect(ok and g(d, "iv_verdict") == "PASS",
               "%-12s (iv) instants distincts, joint mort au 2e point, 3e point jamais rompu [%s] (points %s)"
               % (d, g(d, "iv_verdict"), g(d, "iv_tk")))
    tk = g("tilt_pl_any", "iv_tk") or []
    ok = len(tk) == 3 and all(v is not None for v in tk) and g("tilt_pl_any", "iv_tB") is not None
    if ok:
        dt = g("tilt_pl_any", "dt") or 1e-9
        ok = abs(g("tilt_pl_any", "iv_tB") - min(tk)) <= 0.5 * dt
    expect(ok and g("tilt_pl_any", "iv_verdict") == "conforme", "tilt_pl_any  (iv) jointFailRule = any : mort au PREMIER point - variante qui DOIT differer")
    tsol = g("tilt_sol", "iv_tk"); tany = g("tilt_pl_any", "iv_tB"); tmaj = g("tilt_pl_maj", "iv_tB")
    if tany is not None and tmaj is not None:
        expect(tmaj > tany, "tilt : tBreak(majority) = %.4e s > tBreak(any) = %.4e s sur le meme trajet plastic" % (tmaj, tany))
    # refus
    for d, msg in (("bad_key_percussion", "jointbench"), ("bad_mesh", "mesh")):
        log = os.path.join(S, d + ".log")
        txt = io.open(log, encoding="utf-8", errors="replace").read() if os.path.isfile(log) else ""
        expect(("error" in txt.lower()) and (msg in txt) and ("wall time" not in txt),
               "%-18s refuse a l initialisation (message contenant '%s')" % (d, msg))
    # bitid
    # ancres jouees avec l exe FINAL (rockim_s3final.exe) ; bitid_k / bitid_y = premier exe (17acb6b5)
    for tag, name in (("bitid_kf", "fdem3d_kuru9_court"), ("bitid_yf", "fdem3d_yang_v2_court")):
        log = os.path.join(S, tag + ".log")
        txt = io.open(log, encoding="utf-8", errors="replace").read() if os.path.isfile(log) else ""
        m = re.search(r"(IDENTIQUE|DIFF\w*|ECHEC)\s+%s\s+([0-9.]+) s" % name, txt)
        expect(m is not None and m.group(1) == "IDENTIQUE",
               "ancre bitid %s : %s" % (name, (m.group(1) + " " + m.group(2) + " s") if m else "absente"))
    say("")
    say("BILAN : %d verification(s) en echec" % len(fails) if fails else "BILAN : schema falsifiant entierement tenu")
    with io.open(os.path.join(HERE, "RESULTATS_check_s3.txt"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(OUT) + "\n")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
