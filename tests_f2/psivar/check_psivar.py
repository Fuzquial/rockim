#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# check_psivar.py — BANC COURT de la dilatance variable psi(pbar) du DP-DFH
# (capacite portee de la branche insertion-pointe, commit c2c2d42).
#
#   python tests_f2/psivar/check_psivar.py --exe build/rockim.exe
#   python tests_f2/psivar/check_psivar.py --exe build/rockim.exe \
#          --ref  ../rockim_f2/rockim_f2w21.exe      # + preuve bit-identite
#
# Deux points materiels DP-DFH (aucun maillage, ~1 s chacun) :
#   TEMOIN  mp_dpdfh_psivar_off.cfg  dfhPsiVar absent -> psi = dfhPsiDeg = 15 deg
#   ESSAI   mp_dpdfh_psivar_on.cfg   dfhPsiVar = 1    -> psi = clamp(160,345 -
#                                    0,213793 pbar[MPa], 0, 51,7)
#
# CRITERES (falsifiants, verdict OK / ECHEC imprime) :
#   A. bit-identite : avec --ref, la trace du TEMOIN doit etre RIGOUREUSEMENT
#      celle du binaire d avant le portage (la cle absente ne change rien).
#   B. la cle AGIT : a sigma3 = 0 les deux traces doivent differer.
#   C. psi DECROIT avec p : psi(p) coupe les 15 deg fixes du temoin a
#      pbar = (160,345 - 15)/0,213793 = 679,8 MPa. L ecart de deformation
#      volumique essai - temoin doit donc etre POSITIF (essai plus dilatant)
#      tant que le pbar au seuil reste sous 679,8 MPa, et NEGATIF au-dela.
#      Avec beta = 51,7 deg et dcoh = 153,3 MPa le seuil vaut analytiquement
#      pbar = 1,7297 sigma3 + 88,4 MPa, soit le basculement entre
#      sigma3 = 300 MPa (pbar 607) et sigma3 = 400 MPa (pbar 780).
#      La VARIANTE QUI DOIT ECHOUER : rejouer l ESSAI avec dfhPsiVar = 0 —
#      l ecart tombe alors a zero partout (option --falsify).
# ---------------------------------------------------------------------------
import argparse, csv, collections, hashlib, io, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
PSI0, KPSI, PSIMAX, PSIFIX = 160.345, 0.213793, 51.7, 15.0
PCROSS = (PSI0 - PSIFIX) / KPSI            # 679,8 MPa


def run(exe, cfg, out):
    env = dict(os.environ, OMP_NUM_THREADS="2")
    r = subprocess.run([exe, "matpoint", cfg, out], cwd=ROOT, env=env,
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("echec du run %s :\n%s\n%s" % (cfg, r.stdout, r.stderr))
    return out


def read(p):
    d = collections.OrderedDict()
    for r in csv.DictReader(open(p)):
        d.setdefault(float(r["sigma3"]), []).append(r)
    return d


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default=os.path.join("build", "rockim.exe"))
    ap.add_argument("--ref", default=None, help="binaire d avant le portage (critere A)")
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--falsify", action="store_true",
                    help="rejoue l ESSAI avec dfhPsiVar = 0 : l ecart doit"
                         " tomber a zero partout, donc le critere B doit ECHOUER")
    a = ap.parse_args()
    exe = a.exe if os.path.isabs(a.exe) else os.path.join(ROOT, a.exe)
    tmp = tempfile.mkdtemp(prefix="rockim_psivar_")
    off = os.path.join(HERE, "mp_dpdfh_psivar_off.cfg")
    on = os.path.join(HERE, "mp_dpdfh_psivar_on.cfg")
    fOff = run(exe, off, os.path.join(tmp, "off.csv"))
    fOn = run(exe, on, os.path.join(tmp, "on.csv"))
    ok = True

    print("[psivar] sorties dans %s" % tmp)
    if a.ref:
        ref = a.ref if os.path.isabs(a.ref) else os.path.join(ROOT, a.ref)
        fRef = run(ref, off, os.path.join(tmp, "off_ref.csv"))
        same = sha(fOff) == sha(fRef)
        ok &= same
        print("[A] bit-identite du TEMOIN vs %s : %s\n    %s"
              % (os.path.basename(ref), "OK" if same else "ECHEC", sha(fOff)))
    else:
        print("[A] bit-identite : non evaluee (passer --ref <exe d avant>)")

    dOff, dOn = read(fOff), read(fOn)
    print("\n sigma3 | pbar seuil |  psi(p)  |  eps_vol TEMOIN |  eps_vol ESSAI  |"
          "   ecart   | attendu")
    first = True
    for s3 in dOff:
        A, B = dOff[s3], dOn[s3]
        # pbar au SEUIL = celui du premier pas ou les deux traces divergent
        pb = None
        for ra, rb in zip(A, B):
            if abs(float(ra["eps_vol"]) - float(rb["eps_vol"])) > 1e-12:
                pb = -(float(rb["sig_ax"]) + 2 * float(rb["sig_lat"])) / 3.0 / 1e6
                break
        psi = None if pb is None else max(0.0, min(PSIMAX, PSI0 - KPSI * pb))
        ea, eb = float(A[-1]["eps_vol"]), float(B[-1]["eps_vol"])
        d = eb - ea
        if first:
            ok &= abs(d) > 1e-9                       # critere B
            print("[B] la cle AGIT a sigma3 = 0 : %s (ecart %+.3e)"
                  % ("OK" if abs(d) > 1e-9 else "ECHEC (cle INERTE)", d))
            print("\n sigma3 | pbar seuil |  psi(p)  |  eps_vol TEMOIN |  eps_vol ESSAI  |"
                  "   ecart   | attendu")
            first = False
        if pb is None:
            att, good = "aucune plasticite", True
        elif pb < PCROSS:
            att, good = "psi > 15 -> ecart > 0", d > 0
        else:
            att, good = "psi < 15 -> ecart < 0", d < 0
        ok &= good
        print("%7.0f | %10s | %8s | %15.4e | %15.4e | %+9.2e | %-22s %s"
              % (s3 / 1e6, "n/a" if pb is None else "%.0f MPa" % pb,
                 "n/a" if psi is None else "%.1f deg" % psi, ea, eb, d, att,
                 "OK" if good else "ECHEC"))
    print("\n[C] psi decroissant avec p (changement de signe de l ecart a "
          "pbar = %.1f MPa) : %s" % (PCROSS, "OK" if ok else "ECHEC"))
    if a.falsify:
        # VARIANTE QUI DOIT ECHOUER (annoncee dans l en-tete de ce fichier et
        # dans mp_dpdfh_psivar_on.cfg) : le MEME deck d essai, dfhPsiVar = 0.
        # Si psi(p) est bien pilotee par cette seule cle, la trace redevient
        # RIGOUREUSEMENT celle du temoin — le critere B (« la cle AGIT ») doit
        # donc ECHOUER sur cette variante. Un banc dont la variante negative
        # passerait ne prouverait rien.
        txt = io.open(on, encoding="utf-8").read()
        txt2 = re.sub(r"(?m)^(\s*dfhPsiVar\s*=\s*)1\b", r"\g<1>0", txt)
        assert txt2 != txt, "dfhPsiVar = 1 introuvable dans mp_dpdfh_psivar_on.cfg"
        fz = os.path.join(tmp, "mp_dpdfh_psivar_falsify.cfg")
        io.open(fz, "w", encoding="utf-8", newline="\n").write(txt2)
        fF = run(exe, fz, os.path.join(tmp, "falsify.csv"))
        dF = read(fF)
        worst = max(abs(float(dF[s3][-1]["eps_vol"])
                        - float(dOff[s3][-1]["eps_vol"])) for s3 in dOff)
        same = sha(fOff) == sha(fF)
        ok &= same
        print("\n[FALSIFY] meme deck d essai avec dfhPsiVar = 0 : ecart max "
              "%+.3e, trace %s au TEMOIN — le critere B echoue comme il doit : %s"
              % (worst, "IDENTIQUE" if same else "DIFFERENTE",
                 "OK" if same else "ECHEC (la cle n est pas le seul pilote)"))

    print("\n[psivar] VERDICT : %s" % ("TOUS LES CRITERES PASSES" if ok else "ECHEC"))
    if not a.keep:
        import shutil; shutil.rmtree(tmp, ignore_errors=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
