#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# t4_check_instr.py — verificateur FALSIFIANT de l instrumentation des decks de
# conformite (T4, 2e passe de la campagne du 13/09). Il ne lance rien : il lit un
# dossier de fumee produit par tools/deck_smoke.py et repond par des PASS/FAIL
# mesures sur les fichiers.
#
#   python tools/t4_check_instr.py tests_f2/campagne13/T4/<stem>_<tag> [...] \
#          [--expect-fc rock:insert,piston:bit,plate:bit] [--expect-fields] \
#          [--negative] [--piston-mass 1.057] [--bit-mass 1.288] [--tcontact 1.88e-6]
#
# Controles (chacun imprime la valeur mesuree) :
#   M  masses par corps lues au journal (et ecart aux masses attendues du train
#      T3 a s = 1 quand --piston-mass / --bit-mass sont donnes) ;
#   C  colonnes Fc_<a>_<b>_x/y/z de history.csv : presence, max |F|, et pour
#      rock:insert la NULLITE EXACTE tant que l onde n a pas atteint la roche
#      (a quelques microsecondes l insert est encore immobile) ;
#   B  Fc_piston_bit_z NON nulle apres l instant de contact attendu (jeu / v) ;
#   V  champs `dead` et `openMax` dans le dernier VTU des joints, `pMean` dans
#      celui des elements (writeRuptureFields) ;
#   N  --negative : les MEMES controles doivent ECHOUER (aucune colonne Fc_*,
#      aucun champ dead/openMax/pMean) — le temoin sans les cles.
#
# Code de retour 0 si tous les controles demandes sont conformes, 1 sinon.
# ---------------------------------------------------------------------------
import argparse
import csv
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RE_MASS = re.compile(r"corps '(\w+)': .*masse = ([0-9.eE+-]+) kg")
RE_DA = re.compile(r'<DataArray[^>]*Name="([^"]+)"')


def read_log(run_dir):
    logs = sorted(glob.glob(os.path.join(run_dir, "*.log")))
    if not logs:
        return ""
    with open(logs[0], encoding="utf-8", errors="replace") as f:
        return f.read()


def vtu_names(path, chunk=1 << 20):
    """Noms de TOUS les DataArray du fichier. Les donnees etant ecrites en ASCII
    a la suite de chaque balise, les noms sont disperses dans tout le fichier
    (un VTU de 43 772 joints fait 7 Mo) : lire seulement l en-tete en manquerait
    — c est l erreur qui a fait echouer le premier jet de ce verificateur."""
    names, tail = set(), ""
    with open(path, encoding="utf-8", errors="replace") as f:
        while True:
            blk = f.read(chunk)
            if not blk:
                break
            buf = tail + blk
            names |= set(RE_DA.findall(buf))
            tail = buf[-4096:]          # une balise coupee par la frontiere
    return names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--expect-fc", default="")
    ap.add_argument("--expect-fields", action="store_true")
    ap.add_argument("--negative", action="store_true")
    ap.add_argument("--piston-mass", type=float, default=None)
    ap.add_argument("--bit-mass", type=float, default=None)
    ap.add_argument("--tcontact", type=float, default=None)
    a = ap.parse_args()
    pairs = [p for p in a.expect_fc.split(",") if p]
    bad = 0

    for d in a.dirs:
        dd = os.path.abspath(os.path.join(ROOT, d))
        run = os.path.join(dd, "run")
        print("=" * 78)
        print(os.path.basename(dd))
        log = read_log(dd)
        masses = {k: float(v) for k, v in RE_MASS.findall(log)}
        if masses:
            print("  M masses (kg) : " + ", ".join(
                "%s %.5g" % (k, v) for k, v in masses.items()))
        for nom, att in (("piston", a.piston_mass), ("bit", a.bit_mass)):
            if att is not None:
                got = masses.get(nom)
                if got is None:
                    print("  M %-6s ABSENTE du journal  [FAIL]" % nom); bad += 1
                else:
                    rel = abs(got - att) / att
                    verdict = "PASS" if rel < 0.01 else "FAIL"
                    if verdict == "FAIL":
                        bad += 1
                    print("  M %-6s %.5f kg contre %.5f attendue (%.2f %%)  [%s]"
                          % (nom, got, att, 100 * rel, verdict))

        hist = os.path.join(run, "history.csv")
        cols, rows = [], []
        if os.path.exists(hist):
            with open(hist, newline="", encoding="utf-8", errors="replace") as f:
                rd = csv.reader(f)
                cols = next(rd)
                cols = [c.strip() for c in cols]
                for r in rd:
                    if len(r) == len(cols):
                        rows.append(r)
            print("  C history.csv : %d colonnes, %d lignes ; colonnes Fc_* : %s"
                  % (len(cols), len(rows),
                     ", ".join(c for c in cols if c.startswith("Fc_")) or "(aucune)"))
        else:
            print("  C history.csv ABSENT")

        def col(name):
            return cols.index(name) if name in cols else -1

        def series(name):
            i = col(name)
            return [float(r[i]) for r in rows] if i >= 0 else None

        if a.negative:
            fc = [c for c in cols if c.startswith("Fc_")]
            v = "PASS" if not fc else "FAIL"
            if v == "FAIL":
                bad += 1
            print("  N aucune colonne Fc_* attendue : %d trouvee(s)  [%s]" % (len(fc), v))
        for p in pairs:
            aa, bb = p.split(":")
            names = ["Fc_%s_%s_%s" % (aa, bb, c) for c in "xyz"]
            if any(col(n) < 0 for n in names):
                print("  C %-18s colonnes ABSENTES  [FAIL]" % p); bad += 1
                continue
            mx = {}
            for n in names:
                s = series(n)
                mx[n] = max(abs(x) for x in s) if s else float("nan")
            print("  C %-18s max |Fx| %.4g  |Fy| %.4g  |Fz| %.4g N"
                  % (p, mx[names[0]], mx[names[1]], mx[names[2]]))
            if aa == "rock" or bb == "rock":
                nz = sum(1 for n in names for x in (series(n) or []) if x != 0.0)
                v = "PASS" if nz == 0 else "FAIL"
                if v == "FAIL":
                    bad += 1
                print("     -> nullite EXACTE de la paire roche sur %d lignes : "
                      "%d valeur(s) non nulle(s)  [%s]" % (len(rows), nz, v))
            if a.tcontact is not None and aa == "piston" and bb == "bit":
                t = series("t") or series("time")
                sz = series(names[2])
                if t is None or sz is None:
                    print("     -> colonne t absente  [FAIL]"); bad += 1
                else:
                    av = [abs(z) for tt, z in zip(t, sz) if tt < a.tcontact]
                    ap_ = [abs(z) for tt, z in zip(t, sz) if tt >= a.tcontact]
                    m_av = max(av) if av else 0.0
                    m_ap = max(ap_) if ap_ else 0.0
                    v = "PASS" if (m_av == 0.0 and m_ap > 0.0) else "FAIL"
                    if v == "FAIL":
                        bad += 1
                    print("     -> B max |Fz| avant %.3g s = %.4g N, apres = %.4g N  [%s]"
                          % (a.tcontact, m_av, m_ap, v))

        jv = sorted(glob.glob(os.path.join(run, "*joints_*.vtu")))
        ev = sorted(p for p in glob.glob(os.path.join(run, "fdem*_[0-9][0-9][0-9][0-9].vtu"))
                    if "joints" not in os.path.basename(p))
        for tag, files, want in (("joints", jv, {"dead", "openMax"}),
                                 ("elements", ev, {"pMean"})):
            if not files:
                print("  V VTU %s ABSENT" % tag)
                continue
            names = vtu_names(files[-1])
            got = want & names
            if a.negative:
                v = "PASS" if not got else "FAIL"
                print("  N VTU %-8s %s : champs de S1 attendus ABSENTS, trouves %s  [%s]"
                      % (tag, os.path.basename(files[-1]), sorted(got) or "(aucun)", v))
            elif a.expect_fields:
                v = "PASS" if got == want else "FAIL"
                print("  V VTU %-8s %s : %s  [%s]"
                      % (tag, os.path.basename(files[-1]), sorted(got) or "(aucun)", v))
            else:
                v = "PASS"
                print("  V VTU %-8s %s : %s" % (tag, os.path.basename(files[-1]), sorted(got)))
            if v == "FAIL":
                bad += 1

    print("=" * 78)
    print("controles non conformes : %d" % bad)
    sys.exit(0 if bad == 0 else 1)


if __name__ == "__main__":
    main()
