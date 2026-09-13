# -*- coding: utf-8 -*-
"""check_s1.py — verification falsifiante des mini-tests S1 (campagne 13/09).

usage : python check_s1.py <scratch_dir>
Lit les sorties out_* et bitid_*, imprime chiffres et verdicts.
"""
import glob, hashlib, io, os, re, sys
import numpy as np

S = sys.argv[1]


def read_vtu(path):
    s = io.open(path, encoding="utf-8", errors="ignore").read()
    m = re.search(r'Name="connectivity"[^>]*>\s*(.*?)\s*</DataArray>', s, re.S)
    con = np.fromstring(m.group(1), sep=" ")
    m = re.search(r'Name="offsets"[^>]*>\s*(.*?)\s*</DataArray>', s, re.S)
    off = np.fromstring(m.group(1), sep=" ")
    ncell = off.size
    fields = {}
    for m in re.finditer(r'<DataArray[^>]*Name="([A-Za-z]\w*)"[^>]*>\s*(.*?)\s*</DataArray>', s, re.S):
        nm = m.group(1)
        if nm in ("connectivity", "offsets", "types"):
            continue
        v = np.fromstring(m.group(2), sep=" ")
        if v.size == ncell:
            fields[nm] = v
    return ncell, fields


def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def last(pattern):
    fs = sorted(glob.glob(pattern))
    return fs[-1] if fs else None


def find_out(name):
    """dossier de sortie : out_<name> ou bitid_<x>/<deck>."""
    d = os.path.join(S, "out_" + name)
    if os.path.isdir(d):
        return d
    for d in glob.glob(os.path.join(S, "bitid_*", "*")):
        if os.path.isdir(d) and os.path.basename(d) == name:
            return d
    return None


def joints(d):
    f = last(os.path.join(d, "fdem*_joints_*.vtu"))
    return f, read_vtu(f)[1] if f else (None, None)


def elems(d):
    f = last(os.path.join(d, "fdem3d_[0-9]*.vtu")) or last(os.path.join(d, "fdem_[0-9]*.vtu"))
    return f, read_vtu(f)[1] if f else (None, None)


def cap_lines(log):
    if not os.path.isfile(log):
        return []
    return [l.strip() for l in io.open(log, encoding="utf-8", errors="ignore")
            if "meanTensionCap" in l]


def report_rupt(name):
    d = find_out(name)
    print("\n--- %s (%s)" % (name, d))
    if d is None:
        print("  ABSENT")
        return None
    fj, J = joints(d)
    fe, E = elems(d)
    print("  joints VTU :", os.path.basename(fj), "champs :", sorted(J.keys()))
    print("  elems  VTU :", os.path.basename(fe), "champs :", sorted(E.keys()))
    tb = J["tBreak"]
    nBrk = int((tb >= 0).sum())
    print("  facettes : %d ; rompues (tBreak>=0) : %d ; damage>=0.999 : %d ; bmode 1/2 : %d/%d"
          % (tb.size, nBrk, int((J["damage"] >= 0.999).sum()),
             int((J["breakMode"] == 1).sum()), int((J["breakMode"] == 2).sum())))
    if "dead" in J:
        dead = J["dead"] > 0.5
        viol = int((dead & (tb < 0)).sum())
        om = J["openMax"]
        print("  dead : %d ; dead ET tBreak<0 (doit etre 0) : %d  -> %s"
              % (int(dead.sum()), viol, "OK" if viol == 0 else "ECHEC"))
        print("  openMax : max %.3e m ; rompues ouvertes (tBreak>=0 & openMax>0) : %d / %d ; rompues fermees : %d"
              % (om.max(), int(((tb >= 0) & (om > 0)).sum()), nBrk, int(((tb >= 0) & (om <= 0)).sum())))
        print("  openMax>0 sur facettes NON rompues : %d (ouverture elastique, attendu > 0 possible)"
              % int(((tb < 0) & (om > 0)).sum()))
    if "pMean" in E:
        pm = E["pMean"]
        print("  pMean : min %.3e Pa, max %.3e Pa, >0 : %d / %d" % (pm.min(), pm.max(), int((pm > 0).sum()), pm.size))
    return d, J, E


def compare(nameA, nameB, label):
    dA, dB = find_out(nameA), find_out(nameB)
    print("\n=== %s : %s vs %s" % (label, nameA, nameB))
    if dA is None or dB is None:
        print("  ABSENT", dA, dB)
        return
    for f in ("history.csv", "frames.csv"):
        a, b = os.path.join(dA, f), os.path.join(dB, f)
        same = sha(a) == sha(b)
        print("  %-12s sha256 %s" % (f, "IDENTIQUE" if same else "DIFFERENT"))
    _, JA = joints(dA)
    _, JB = joints(dB)
    for k in ("damage", "tBreak", "breakMode"):
        if k in JA and k in JB:
            nd = int((JA[k] != JB[k]).sum())
            print("  joints.%-10s : %d valeurs differentes / %d" % (k, nd, JA[k].size))
    if "failMode" in JA and "failMode" in JB:
        nd = int((JA["failMode"] != JB["failMode"]).sum())
        print("  joints.failMode  : %d valeurs differentes" % nd)
    _, EA = elems(dA)
    _, EB = elems(dB)
    for k in ("vonMises", "sigma1", "sigmaXX", "bulkD"):
        if k in EA and k in EB:
            nd = int((EA[k] != EB[k]).sum())
            print("  elems.%-10s : %d valeurs differentes / %d" % (k, nd, EA[k].size))


for nm in sys.argv[2:]:
    pass

print("##### 2D (fdem) #####")
for nm in ("u_rupt_slipRef", "u_co_slipF", "u_co_slipRef", "u2_ref", "u2_slipRef", "u2_co_slipF", "u2_co_slipRef"):
    report_rupt(nm)
compare("u2_ref", "u2_slipRef", "2D UCS 6 ms : slipRef sans coulomb vs ref (0 etiquette changee attendu)")
compare("u2_co_slipF", "u2_co_slipRef", "2D UCS 6 ms : slipF vs slipRef sous coulomb")
compare("fdem_ucs_yan_adaptive_court", "u_rupt_slipRef", "cles armees (sans coulomb) vs reference bitid")
compare("u_co_slipF", "u_co_slipRef", "jointBreakModeRef slipF vs slipRef sous coulomb")
print("\n=== compteur du cap 2D")
for nm, log in (("bitid ref", os.path.join(S, "bitid_u", "fdem_ucs_yan_adaptive_court.log")),
                ("u_cap", os.path.join(S, "u_cap.log")), ("u_rupt_slipRef", os.path.join(S, "u_rupt_slipRef.log"))):
    ls = cap_lines(log)
    if not ls:
        ls = cap_lines(last(os.path.join(S, "bitid_u", "*", "*.log")) or "") if nm == "bitid ref" else ls
    print("  [%s] %d lignes" % (nm, len(ls)))
    for l in ls[:6]:
        print("     ", l)

print("\n##### 3D (fdem3d) #####")
for nm in ("k9_rupt", "k9_pl_slipF", "k9_pl_slipRef", "h_pl_slipF", "h_pl_slipRef", "h_ori_slipRef", "v_slipRef", "v_co_slipF", "v_co_slipRef"):
    report_rupt(nm)
compare("fdem3d_cut3d_heilman_court", "h_ori_slipRef", "heilman origin : cles armees vs ref bitid")
compare("h_pl_slipF", "h_pl_slipRef", "heilman plastic+coulomb : slipF vs slipRef")
compare("fdem3d_visc_yan_3d", "v_slipRef", "visc_yan tension : slipRef sans coulomb vs ref bitid")
compare("v_co_slipF", "v_co_slipRef", "visc_yan tension coulomb : slipF vs slipRef")
compare("fdem3d_kuru9_court", "k9_rupt", "writeRuptureFields vs reference bitid")
compare("k9_pl_slipF", "k9_pl_slipRef", "jointBreakModeRef slipF vs slipRef (plastic + coulomb)")
print("\n=== compteur du cap 3D")
for nm, log in (("bitid ref", os.path.join(S, "bitid_k", "fdem3d_kuru9_court.log")),
                ("k9_cap", os.path.join(S, "k9_cap.log")), ("k9_cap2", os.path.join(S, "k9_cap2.log")), ("k9_rupt", os.path.join(S, "k9_rupt.log"))):
    ls = cap_lines(log)
    if not ls and nm == "bitid ref":
        ls = cap_lines(last(os.path.join(S, "bitid_k", "*", "*.log")) or "")
    print("  [%s] %d lignes" % (nm, len(ls)))
    for l in ls[:6]:
        print("     ", l)
