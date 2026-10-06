"""Auto-tests de contact par potentiel : collision de deux corps rigides sans
frottement ni amortissement (selftest-potential2d, -potential3d, -potvolume3d).
(a) vitesses en x du choc frontal ; (b) energie cinetique rapportee a Ec0 ;
(c) ecarts relatifs finaux d'energie cinetique et de quantite de mouvement.
Donnees : collision_potentiel.csv (traces CSV des auto-tests, echantillonnees
tous les 0,02 s) et collision_bilan.csv (lignes [POT] des journaux),
rejoues le 2026-10-06. Usage : python fig_collision.py -> fig_collision.pdf
"""
import os, csv
import numpy as np
from style_verif import plt, virgule

HERE = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(HERE, "collision_potentiel.csv"))))
bil = list(csv.DictReader(open(os.path.join(HERE, "collision_bilan.csv"))))
TESTS = [("potentiel_2D", "potentiel 2D", "#1f4e79"), ("potentiel_3D", "potentiel 3D", "#c55a11"),
         ("volume_3D", "volume 3D", "#548235")]


def get(test, ph):
    r = [x for x in rows if x["test"] == test and x["phase"] == ph]
    return {k: np.array([float(x[k]) for x in r]) for k in ("t", "vAx", "vBx", "Px", "W_contact", "Ec")}


fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
for test, lab, c in TESTS:
    d = get(test, "1")
    ax[0].plot(d["t"], d["vAx"], color=c, lw=1.0, label=lab)
    ax[0].plot(d["t"], d["vBx"], color=c, lw=1.0, ls="--")
    for ph, ls in (("1", "-"), ("2", ":")):
        d = get(test, ph)
        ax[1].plot(d["t"], d["Ec"] / d["Ec"][0], color=c, ls=ls, lw=1.1,
                   label=lab + (" frontal" if ph == "1" else " oblique"))
ax[0].set_xlabel("temps (unités du test)")
ax[0].set_ylabel(r"vitesse $v_x$")
ax[0].set_title("(a) choc frontal : A (trait), B (tirets)", fontsize=9)
ax[0].legend(loc="center right")
ax[1].set_xlabel("temps (unités du test)")
ax[1].set_ylabel(r"$E_c / E_{c0}$")
ax[1].set_title("(b) énergie cinétique", fontsize=9)
ax[1].legend(loc="lower right", fontsize=6.5)
for a in ax[:2]:
    virgule(a)
# (c) ecarts finaux
a = ax[2]
x = np.arange(6)
NOMT = {"potentiel_2D": "pot. 2D", "potentiel_3D": "pot. 3D", "volume_3D": "vol. 3D"}
lab = [f"{NOMT[b['test']]}, {'frontal' if b['phase'] == 'frontale' else 'oblique'}" for b in bil]
dE = [float(b["dEc_sur_Ec0"]) for b in bil]
dP = [float(b["dP_sur_P0"]) for b in bil]
a.bar(x - 0.18, dE, 0.36, color="#1f4e79", label=r"$|\Delta E_c|/E_{c0}$")
a.bar(x + 0.18, dP, 0.36, color="#bfbfbf", label=r"$|\Delta \mathbf{P}|/|\mathbf{P}_0|$")
a.set_yscale("log")
a.set_ylim(1e-17, 1e-3)
a.set_xticks(x)
a.set_xticklabels(lab, fontsize=6.5, rotation=40, ha="right", rotation_mode="anchor")
a.set_title("(c) écarts relatifs en fin d'essai", fontsize=9)
a.legend(loc="upper right", fontsize=7, ncol=1)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_collision.pdf"))
