"""Banc T0b : chaine de N masses (barre de Saint-Venant, granite Utah FORGE,
L = 0,1 m, c = 4 300 m/s) frappee a 10 m/s par le mur rigide de Signorini (e = 0
au noeud). (a) histoire a N = 100 ; (b) vitesse de sortie / 2v et duree de
contact / (2L/c) contre N ; (c) perte d'energie contre N, attendue 1/N.
Donnees : t0b_balayage_N.csv et t0b_histoire_N100.csv, produits le 2026-10-06
par t0b_sweep.cpp (copie de la lambda `bar` du selftest, vrai noyau
toolsig::impulse) ; N = 100 et 400 redonnent les chiffres du selftest.
Usage : python fig_t0b_chaine.py -> fig_t0b_chaine.pdf
"""
import os, csv
import numpy as np
from style_verif import plt, virgule

HERE = os.path.dirname(os.path.abspath(__file__))
b = list(csv.DictReader(open(os.path.join(HERE, "t0b_balayage_N.csv"))))
h = np.genfromtxt(os.path.join(HERE, "t0b_histoire_N100.csv"), delimiter=",", names=True)
T2 = 46.5111  # 2L/c en microsecondes

fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5))
a = ax[0]
a.plot(h["t_us"], h["v_moy"], color="#1f4e79", label="vitesse moyenne")
a.plot(h["t_us"], h["v_noeudN"], color="#c55a11", lw=0.8, label="extrémité libre")
a.axhline(20, color="k", lw=0.5, ls=":")
a.axvline(T2, color="k", lw=0.5, ls=":")
a.set_xlim(0, 58)
a.set_xlabel(r"temps ($\mu$s)")
a.set_ylabel("vitesse (m/s)")
a.set_title(r"(a) $N = 100$, mur à 10 m/s", fontsize=9)
a.text(T2 + 2, 2, r"$2L/c$", fontsize=8)
a.legend(loc="upper left", fontsize=7)
virgule(a)

for amort, mk, lab in (("0.00", "o", "sans amortissement"), ("0.05", "s", r"Cundall 0,05")):
    r = [x for x in b if x["amort"] == amort]
    N = np.array([int(x["N"]) for x in r])
    v = np.array([float(x["v_sur_2v"]) for x in r])
    tc = np.array([float(x["t_sur_2Lc"]) for x in r])
    ax[1].semilogx(N, v, marker=mk, color="#1f4e79", ms=4,
                   ls="-" if amort == "0.00" else "--", label=r"$v_s/2v$, " + lab)
    if amort == "0.00":
        ax[1].semilogx(N, tc, marker="^", color="#c55a11", ms=4, label=r"$t_c/(2L/c)$, " + lab)
        loss = np.array([float(x["perte_rel"]) for x in r])
        ax[2].loglog(N, loss, "o", color="#1f4e79", ms=4, label="perte mesurée")
        ax[2].loglog(N, 1.0 / N, color="k", lw=0.8, label=r"$1/N$")
        ax[2].loglog(N, 1 - v, "s", color="#c55a11", ms=4, label=r"$1 - v_s/2v$")
ax[1].axhline(1, color="k", lw=0.5, ls=":")
ax[1].set_xlabel(r"nombre de masses $N$")
ax[1].set_ylabel("rapport à la théorie")
ax[1].set_title("(b) sortie et durée de contact", fontsize=9)
ax[1].legend(loc="lower right", fontsize=6.5)
virgule(ax[1], "y")
ax[2].set_xlabel(r"nombre de masses $N$")
ax[2].set_ylabel(r"perte rapportée à $\frac{1}{2}Mv^2$")
ax[2].set_title("(c) perte d'énergie", fontsize=9)
ax[2].legend(loc="lower left", fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig_t0b_chaine.pdf"))
