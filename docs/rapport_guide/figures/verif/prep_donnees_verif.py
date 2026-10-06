"""Prepare les CSV des figures de verification (chapitre 3 du rapport-guide).

Les donnees viennent de calculs REJOUES le 2026-10-06 (Linux, g++, un fil,
binaire compile hors depot depuis le commit 170786c) :
  - bancs de joint : tests_f2/campagne13/S3/*.cfg et S3bis/*.cfg, sorties
    jointbench.csv et journaux ;
  - auto-tests de contact : rockim selftest-potential2d | potential3d |
    potvolume3d | toolcontact <csv> ;
  - chaine T0b : copie de la lambda `bar` de src/FdemSolver.cpp (cas C8/C9)
    appelant le vrai noyau include/rockim/ToolSignorini.hpp, balayee en N
    (t0b_sweep.cpp) ; N = 100 et 400 redonnent au chiffre pres le selftest ;
  - bilan B4 : residus lus dans les journaux de huit calculs courts.
Usage : python prep_donnees_verif.py <dossier des calculs rejoues>
Ecrit les CSV a cote de ce script. Aucune donnee n'est inventee : toute
valeur ecrite est lue dans un fichier de sortie de rockim.
"""
import os, re, sys, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else "."
A0 = 4.33013e-07          # aire de la facette du banc (journal [JOINTBENCH])


def jb(name):
    return np.genfromtxt(os.path.join(RUNS, "jb", "out_" + name, "jointbench.csv"),
                         delimiter=",", names=True)


def avg(d, key):
    return (d[key + "0"] + d[key + "1"] + d[key + "2"]) / 3.0


# 1. traction et cisaillement monotones avec decharge (S3)
rows = []
for law in ("plastic", "origin", "solidity"):
    for mode, pre in (("traction", "ten"), ("cisaillement", "sh")):
        d = jb(f"{pre}_{law}")
        dead = np.where(d["dead"] > 0)[0]
        n = dead[0] + 1 if len(dead) else len(d)
        step = max(1, n // 1500)
        for i in range(0, n, step):
            rows.append([law, mode, d["t"][i], d["phase"][i], avg(d, "dn")[i],
                         avg(d, "ds")[i], avg(d, "sig")[i], avg(d, "tau")[i],
                         avg(d, "D")[i]])
with open(os.path.join(HERE, "jb_monotone.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["loi", "trajet", "t_s", "phase", "dn_m", "ds_m", "sig_Pa", "tau_Pa", "D"])
    w.writerows(rows)

# 2. cycles fermes (S3bis) : dernier cycle, et energie par cycle
rows, erows = [], []
for regime, pre in (("intact", "cyce"), ("endommage", "cyc")):
    for law in ("plastic", "origin", "solidity"):
        for suf, k in (("", 1), ("_dt2", 2), ("_dt4", 4)):
            name = f"{pre}_{law}{suf}"
            log = open(os.path.join(RUNS, "jb", name + ".log")).read()
            dt = float(re.search(r"cycle\(s\) de \d+ pas \(dt = (\S+) s\)", log).group(1))
            wlast = float(re.search(r"W net du DERNIER cycle = (\S+) J", log).group(1))
            ampl = float(re.search(r"amplitude du travail dans le cycle = (\S+) J", log).group(1))
            dmax = float(re.search(r"D max sur le cycle = (\S+)", log).group(1))
            d = jb(name)
            # troisieme cycle : lignes ou deux cycles sont acheves, precedees du
            # dernier echantillon du deuxieme (= point de depart du troisieme)
            i = np.where(d["cyc"] == 2)[0]
            i = np.r_[i[0] - 1, i]
            dn, ds = d["dn0"][i], d["ds0"][i]
            sig, tau = avg(d, "sig")[i], avg(d, "tau")[i]
            wjw = d["jw"][i[-1]] - d["jw"][i[0]]
            # travail des tractions du joint sur le tetraedre B (W > 0 : creation)
            wg = -A0 * np.sum(sig[:-1] * np.diff(dn) + tau[:-1] * np.diff(ds))
            wt = -A0 * np.sum(0.5 * (sig[1:] + sig[:-1]) * np.diff(dn)
                              + 0.5 * (tau[1:] + tau[:-1]) * np.diff(ds))
            erows.append([regime, law, k, dt, wlast, wjw, wg, wt, ampl, dmax])
            if k == 1:
                step = max(1, len(i) // 2000)
                for j in range(0, len(i), step):
                    rows.append([regime, law, dn[j], ds[j], sig[j], tau[j], d["phase"][i[j]]])
with open(os.path.join(HERE, "jb_cycles.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["regime", "loi", "dn_m", "ds_m", "sig_Pa", "tau_Pa", "phase"])
    w.writerows(rows)
with open(os.path.join(HERE, "jb_energie_cycle.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["regime", "loi", "division_dt", "dt_s", "W_journal_J", "W_jw_J",
                "W_rectangle_J", "W_trapeze_J", "amplitude_J", "Dmax"])
    w.writerows(erows)

# 3. collisions par potentiel
rows = []
for test, m in (("potentiel_2D", 3.0), ("potentiel_3D", 4.0), ("volume_3D", 4.0)):
    fn = {"potentiel_2D": "potential2d", "potentiel_3D": "potential3d",
          "volume_3D": "potvolume3d"}[test]
    d = np.genfromtxt(os.path.join(RUNS, "selftest", fn + ".csv"), delimiter=",", names=True)
    vA = d["vA"] if "vA" in d.dtype.names else d["vAx"]
    vB = d["vB"] if "vB" in d.dtype.names else d["vBx"]
    for r, a, b in zip(d, vA, vB):
        rows.append([test, int(r["phase"]), r["t"], a, b, m * (a + b), r["W"], r["KE"]])
with open(os.path.join(HERE, "collision_potentiel.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["test", "phase", "t", "vAx", "vBx", "Px", "W_contact", "Ec"])
    w.writerows(rows)

# 4. bilan des collisions : lignes [POT]/[POT3] des journaux des auto-tests
out = []
for test, fn in (("potentiel_2D", "potential2d"), ("potentiel_3D", "potential3d"),
                 ("volume_3D", "potvolume3d")):
    s = open(os.path.join(RUNS, "selftest", fn + ".log")).read()
    for ph, m in zip(("frontale", "oblique"), re.finditer(
            r"\|W_contact\|/KE0 = (\S+), \|dKE\|/KE0 = (\S+), \|dP\|/\|P0\| = (\S+)", s)):
        out.append([test, ph, m.group(1).rstrip(","), m.group(2).rstrip(","), m.group(3)])
with open(os.path.join(HERE, "collision_bilan.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["test", "phase", "W_contact_sur_Ec0", "dEc_sur_Ec0", "dP_sur_P0"])
    w.writerows(out)
