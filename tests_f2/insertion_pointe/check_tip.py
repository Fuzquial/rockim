#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# check_tip.py — BANC COURT des deux capacites d insertion portees de la
# branche insertion-pointe : `insertionTipFactor` (enveloppe relachee en
# pointe de fissure) et `insertion = none` (continuum pur).
#
#   python tests_f2/insertion_pointe/check_tip.py --exe build/rockim.exe
#
# Quatre runs 2D courts (compression uniaxiale entre plateaux, maillage
# voronoi 2 399 elements, arret automatique apres le pic ; ~2 a 3 min chacun
# a OMP_NUM_THREADS = 2). CRITERES falsifiants :
#
#   A. DEFAUT NEUTRE     ins_tip_ref (cle absente) et ins_tip_1
#      (insertionTipFactor = 1) doivent rendre un history.csv BIT-IDENTIQUE.
#      C est ce qui autorise a porter la capacite sans changer l ancre.
#   B. LA CLE AGIT       ins_tip_16 (facteur 1,6) doit rendre un history.csv
#      DIFFERENT, et imprimer le compte propagations / nucleations.
#   C. L AMORCAGE EST INTACT   c est la raison d etre du choix « relacher la
#      pointe » plutot que « penaliser la nucleation » : la resistance
#      macroscopique mesuree (pic de contrainte axiale) doit rester la MEME a
#      1e-3 MPa pres. La variante qui DOIT echouer serait une penalisation de
#      la nucleation : elle deplacerait le pic (mesure du 2026-08-24 : +60 %).
#   D. CONTINUUM PUR     ins_none doit annoncer insertion = none, n inserer
#      AUCUN joint, n en rompre aucun, et tourner a un pas de temps PLUS
#      GRAND que l adaptatif (le ressort de penalite des joints sort du
#      budget de stabilite).
# ---------------------------------------------------------------------------
import argparse, hashlib, os, re, subprocess, sys, tempfile, shutil

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
DECKS = ["ins_tip_ref", "ins_tip_1", "ins_tip_16", "ins_none"]


def run(exe, deck, outdir, threads):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    cfg = os.path.join("tests_f2", "insertion_pointe", deck + ".cfg")
    r = subprocess.run([exe, cfg, outdir], cwd=ROOT, env=env,
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("echec du run %s (rc %d) :\n%s\n%s"
                         % (deck, r.returncode, r.stdout[-2000:], r.stderr[-2000:]))
    return r.stdout


def num(log, rx, cast=float):
    m = re.search(rx, log)
    return cast(m.group(1)) if m else None


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default=os.path.join("build", "rockim.exe"))
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    exe = a.exe if os.path.isabs(a.exe) else os.path.join(ROOT, a.exe)
    tmp = tempfile.mkdtemp(prefix="rockim_tip_")
    print("[tip] sorties dans %s" % tmp)
    log, hh = {}, {}
    for d in DECKS:
        o = os.path.join(tmp, d)
        log[d] = run(exe, d, o, a.threads)
        hh[d] = sha(os.path.join(o, "history.csv"))
        ins = num(log[d], r"adaptive insertion: (\d+) / \d+ joints inserted", int)
        brk = num(log[d], r"joints inserted \([^)]*\), (\d+) fully broken", int)
        pk = num(log[d], r"peak axial stress = ([0-9.eE+-]+) MPa")
        dt = num(log[d], r"dt = ([0-9.eE+-]+) s")
        print("  %-12s dt %.5e s  inseres %-5s rompus %-5s  pic %s MPa  %s"
              % (d, dt, ins, brk, pk, hh[d][:16]))
    ok = True

    same = hh["ins_tip_ref"] == hh["ins_tip_1"]
    ok &= same
    print("\n[A] defaut NEUTRE (cle absente == insertionTipFactor = 1) : %s"
          % ("OK" if same else "ECHEC"))

    diff = hh["ins_tip_16"] != hh["ins_tip_ref"]
    cnt = re.search(r"insertions en POINTE : (\d+) propagations / (\d+) nucleations",
                    log["ins_tip_16"])
    ok &= diff and cnt is not None
    print("[B] la cle AGIT (facteur 1,6) : %s%s"
          % ("OK" if diff else "ECHEC (trace identique = cle inerte)",
             "" if cnt is None else " — %s propagations / %s nucleations"
                                   % (cnt.group(1), cnt.group(2))))
    if cnt is None:
        print("    ECHEC : le compte propagations / nucleations n a pas ete imprime")

    p0 = num(log["ins_tip_ref"], r"peak axial stress = ([0-9.eE+-]+) MPa")
    p1 = num(log["ins_tip_16"], r"peak axial stress = ([0-9.eE+-]+) MPa")
    intact = p0 is not None and p1 is not None and abs(p1 - p0) < 1e-3
    ok &= intact
    print("[C] amorcage INTACT (pic %.4f -> %.4f MPa) : %s"
          % (p0 or -1, p1 or -1, "OK" if intact else "ECHEC"))

    L = log["ins_none"]
    dtN = num(L, r"dt = ([0-9.eE+-]+) s")
    dtA = num(log["ins_tip_ref"], r"dt = ([0-9.eE+-]+) s")
    annonce = "insertion = none : continuum pur" in L
    jamais = "joints inserted" not in L
    plusGrand = dtN is not None and dtA is not None and dtN > dtA
    good = annonce and jamais and plusGrand
    ok &= good
    print("[D] continuum pur (annonce %s, aucun joint insere %s, dt x%.3f %s) : %s"
          % (annonce, jamais, (dtN / dtA) if (dtN and dtA) else 0.0, plusGrand,
             "OK" if good else "ECHEC"))

    print("\n[tip] VERDICT : %s" % ("TOUS LES CRITERES PASSES" if ok else "ECHEC"))
    if not a.keep:
        shutil.rmtree(tmp, ignore_errors=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
