#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_exo_sweep.py — la signature d'obscuration : N(ppoint), rockim contre le
# banc 6 Abaqus.
#
#   python exo_tunnel/tools/fig_exo_sweep.py
# ---------------------------------------------------------------------------
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "axes.unicode_minus": False,
})

# ppoint [MPa/us], N mesure dans rockim (methode angulaire, couronne 30-55 mm)
PDOT = np.array([0.83, 8.3, 25.0, 83.0, 250.0])
NROC = np.array([5, 10, 12, 20, 20], float)
# N_span publie (RESULTATS_bench6.md, Abaqus + vumat_hole.f)
NPUB = np.array([4, 6, 7, 10, 15], float)

fig, ax = plt.subplots(figsize=(7.4, 5.0))
ax.loglog(PDOT, NROC, "o-", color="#b3202f", lw=1.8, ms=7,
          label="rockim — FEM pur + law = dpdfh")
ax.loglog(PDOT, NPUB, "s--", color="#444", lw=1.4, ms=6,
          label="Abaqus + vumat_hole.f (banc 6)")

# pentes ajustees sur le regime dynamique (les 4 points > 1 MPa/us)
for y, c, nom in ((NROC, "#b3202f", "rockim"), (NPUB, "#444", "banc 6")):
    b, a = np.polyfit(np.log(PDOT[1:]), np.log(y[1:]), 1)
    ax.plot(PDOT[1:], np.exp(a) * PDOT[1:] ** b, ":", color=c, lw=1.0)
    print("%s : exposant dynamique b = %.3f" % (nom, b))

ax.set_xlabel(r"$\dot p$ [MPa/$\mu$s]")
ax.set_ylabel("nombre d'armes radiales")
ax.set_title("Signature d'obscuration : le nombre de fissures radiales\n"
             "croît avec la vitesse de chargement", fontsize=12)
ax.grid(alpha=0.3, which="both")
ax.legend(frameon=False, fontsize=9, loc="upper left")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig("exo_tunnel/fig_exo_sweep." + ext, dpi=165)
print("ecrit : exo_tunnel/fig_exo_sweep")
