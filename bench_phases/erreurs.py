# -*- coding: utf-8 -*-
"""erreurs.py — banc des CAS FAUTIFS du materiau par phase en fem3d.

Regle du depot : une cle inconnue, ou lue et sans effet, doit lever une ERREUR
EXPLICATIVE — jamais etre ignoree en silence. Une cle qui serait acceptee,
validee et inerte est le motif interdit numero 1 : le deck a l'air de dire
quelque chose et le calcul ne l'entend pas.

Chaque cas ci-dessous DOIT echouer, et le message DOIT contenir le fragment
attendu (c'est le message qui est teste, pas seulement le code de retour :
un refus qui envoie chercher au mauvais endroit coute une demi-journee).

usage : python bench_phases/erreurs.py [chemin/vers/rockim.exe]
"""
import os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

BASE_FILE = """mode = fem3d
mesh = file
meshFile = bench_phases/bar2.msh
law = elastic
scenario = tension
T = 1e-6
frames = 1
pullV = 0.01
E = 60e9
nu = 0.25
rho = 2650
ft = 10e6
cohesion = 30e6
frictionDeg = 40
Gf = 60
"""

BASE_GRID = BASE_FILE.replace("mesh = file\nmeshFile = bench_phases/bar2.msh\n",
                              "mesh = grid\nW = 0.01\nD = 0.01\nH = 0.01\n"
                              "nx = 4\nny = 4\nnz = 4\n")

DEUX = ("phases = dur mou\n"
        "phase.dur.fraction = 0.5\nphase.dur.E = 83.1e9\n"
        "phase.mou.fraction = 0.5\nphase.mou.E = 29.3e9\n")

CAS = [
    ("phases sur mesh = grid",
     BASE_GRID + DEUX, "AUCUNE source de phase"),
    ("loi par phase (phase.<nom>.law)",
     BASE_FILE + DEUX + "phase.dur.law = cdp\n", "SEULE loi de volume"),
    ("cle de LOI ecrite par phase",
     BASE_FILE + DEUX + "phase.dur.erodeD = 0.9\n",
     "n'est pas une propriete de phase"),
    ("nom de phase inconnu",
     BASE_FILE + DEUX + "phase.durr.E = 70e9\n", "n'est pas declaree"),
    ("propriete de JOINT (gbAlphaTen)",
     BASE_FILE + DEUX + "gbAlphaTen = 0.3\n", "propriete de JOINT"),
    ("propriete de JOINT (gb.<a>.<b>.ft)",
     BASE_FILE + DEUX + "gb.dur.mou.ft = 2e6\n", "propriete de JOINT"),
    ("brasage de groupes (groupBond)",
     BASE_FILE + DEUX + "groupBond.dur.mou = 1\n", "deja SOUDEE"),
    ("groupe physique sans phase",
     BASE_FILE + "phases = quartz biotite\n"
     "phase.quartz.fraction = 0.5\nphase.biotite.fraction = 0.5\n",
     "n'a pas de phase"),
    ("phases sans groupes nommes",
     BASE_FILE.replace("bar2.msh", "bar1.msh") + DEUX,
     "exige des groupes physiques nommes"),
    ("Weibull x phases sans phaseWeibull",
     BASE_FILE + DEUX + "matWeibullM = 3\n", "se COMPOSENT"),
    ("valeur non physique (E < 0)",
     BASE_FILE + DEUX.replace("phase.mou.E = 29.3e9", "phase.mou.E = -1"),
     "E must be > 0"),
    ("valeur non physique (nu >= 0.5)",
     BASE_FILE + DEUX + "phase.mou.nu = 0.6\n", "nu must be in"),
    ("fraction manquante",
     BASE_FILE + "phases = dur mou\nphase.dur.fraction = 0.5\n"
     "phase.dur.E = 83.1e9\nphase.mou.E = 29.3e9\n", "fraction"),
    # Sans la cle `phases`, PhaseSet nomme sa fiche unique « rock » : une
    # fiche phase.rock.* serait acceptee, consommee et INERTE. Trou ferme.
    ("fiche de phase sans cle `phases`",
     BASE_FILE + "phase.rock.E = 70e9\n", "aucune phase n'est declaree"),
    ("groupPhase hors mesh = file",
     BASE_FILE.replace("mesh = file\nmeshFile = bench_phases/bar2.msh\n",
                       "mesh = grid\nW = 0.01\nD = 0.01\nH = 0.01\n"
                       "nx = 4\nny = 4\nnz = 4\n")
     + "groupPhase.dur = mou\n", "n'a de sens qu'avec mesh = file"),
    ("phaseWeibull sans objet",
     BASE_FILE + DEUX + "phaseWeibull = true\n", "n'autoriserait rien"),
]


def main(exe):
    tmp = tempfile.mkdtemp(prefix="rockim_err_")
    ok = True
    for name, deck, want in CAS:
        p = os.path.join(tmp, "deck.cfg")
        with open(p, "w", newline="\n") as f:
            f.write(deck)
        r = subprocess.run([exe, p, os.path.join(tmp, "out")],
                           cwd=ROOT, capture_output=True, text=True)
        out = (r.stdout or "") + (r.stderr or "")
        good = r.returncode != 0 and want in out
        if not good:
            ok = False
        print("[%s] %-38s %s" % ("OK  " if good else "ECHEC", name,
                                 "refuse avec le bon message"
                                 if good else
                                 "code %d, message : %s"
                                 % (r.returncode,
                                    " / ".join(out.strip().splitlines()[-2:])
                                    or "(aucun)")))
    print("\nVERDICT : %s" % ("les %d cas fautifs sont refuses avec un message "
                              "explicatif" % len(CAS) if ok
                              else "AU MOINS UN CAS FAUTIF PASSE OU MENT"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1
                  else os.path.join(ROOT, "build", "rockim.exe")))
