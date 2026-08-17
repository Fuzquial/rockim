# -*- coding: utf-8 -*-
"""Les deux MODES de rupture, en couleurs : BLEU = mode I (ouverture,
traction), ROUGE = mode II (glissement, cisaillement).

Le mode est attribue par rockim a l'instant ou D atteint 1, en comparant
l'avancement en ouverture rn = (dn-dnE)/(dnF-dnE) a celui en glissement
rs = |slip|/slipF (eq. 16 de Yan) : bmode = 1 si rn >= rs, sinon 2.

Rappel de vocabulaire : il n'y a PAS de fissure « de compression ». Un joint
cohesif s'ouvre (mode I) ou glisse (mode II) ; la compression le FERME et
renforce sa resistance au glissement via le terme de Coulomb c + sigma_n tan(phi).
Ce qu'on appelle rupture en compression dans une roche est ce mode II sous
contrainte normale compressive — sous l'indenteur, precisement.

  python plot_modes_rupture.py [run]        (defaut out_imp3d_ultra)
"""
import os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
RUN = sys.argv[1] if len(sys.argv) > 1 else "out_imp3d_ultra"
A_HERTZ = 1.45          # rayon de contact de Hertz a delta = 0,19 mm [mm]
BLEU, ROUGE = "#1f6fb4", "#c0392b"

D = os.path.join(ROOT, RUN)
f = sorted(x for x in os.listdir(D)
           if re.fullmatch(r"fdem3d_joints_\d{4}\.vtu", x))[-1]
s = open(os.path.join(D, f), encoding="utf-8", errors="ignore").read()


def arr(n, d=float):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n, s, re.S)
    return None if not m else np.fromstring(m.group(1), sep=" ", dtype=d)


P = np.fromstring(re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>",
                            s, re.S).group(1), sep=" ").reshape(-1, 3) * 1e3
conn, off = arr("connectivity", int), arr("offsets", int)
dmg, bm = arr("damage"), arr("breakMode")

idx = np.where(dmg >= 1.0)[0]
tri = np.array([P[conn[(0 if i == 0 else off[i - 1]):off[i]][:3]] for i in idx])
mode = bm[idx]
cen = tri.mean(axis=1)
r = np.hypot(cen[:, 0] - 60, cen[:, 1] - 60)
nI, nII = int((mode == 1).sum()), int((mode == 2).sum())

# cadrage sur l'etendue reelle des fissures
rad = max(3.0, 1.25 * r.max())
dep = max(3.0, 1.25 * (120 - cen[:, 2]).max())

fig, ax = plt.subplots(1, 3, figsize=(15.5, 5.2),
                       gridspec_kw={"width_ratios": [1, 1, 0.85]})

for k, (cols, xl, yl, xlim, ylim, titre) in enumerate((
        ([0, 1], "x (mm)", "y (mm)", (60 - rad, 60 + rad), (60 - rad, 60 + rad),
         "(a) Vue de dessus"),
        ([0, 2], "x (mm)", "z (mm)", (60 - rad, 60 + rad), (120 - dep, 121),
         "(b) Coupe verticale"))):
    a = ax[k]
    # mode II dessine EN DERNIER : c'est le moins nombreux et le plus confine
    for v, col in ((1, BLEU), (2, ROUGE)):
        m = mode == v
        if m.any():
            a.add_collection(PolyCollection(tri[m][:, :, cols], facecolors=col,
                                            edgecolors="none", alpha=0.75))
    if k == 0:
        a.add_patch(plt.Circle((60, 60), A_HERTZ, fill=False, ec="k", lw=1.6,
                              ls="--"))
        a.text(60, 60 + A_HERTZ + 0.12 * rad, "contact de Hertz  a = %.2f mm"
               % A_HERTZ, fontsize=8.5, ha="center")
    else:
        a.axhline(120, color="0.4", lw=1.1)
        a.text(60 - 0.95 * rad, 120.3, "surface du bloc", fontsize=8, color="0.4")
    a.set_xlim(*xlim); a.set_ylim(*ylim); a.set_aspect("equal")
    a.set_xlabel(xl); a.set_ylabel(yl)
    a.set_title(titre, fontsize=10.5)
    if k == 0:
        a.legend(handles=[Line2D([], [], marker="s", ls="", ms=9, color=BLEU,
                                 label="mode I — ouverture / traction (%d)" % nI),
                          Line2D([], [], marker="s", ls="", ms=9, color=ROUGE,
                                 label="mode II — glissement / cisaillement (%d)" % nII)],
                 fontsize=8.5, loc="upper right")

# (c) profil radial : la separation des deux modes, quantifiee
a = ax[2]
bins = np.linspace(0, min(rad, 6.0), 13)
for v, col, lab in ((1, BLEU, "mode I (traction)"),
                    (2, ROUGE, "mode II (cisaillement)")):
    h, _ = np.histogram(r[mode == v], bins=bins)
    a.step(0.5 * (bins[1:] + bins[:-1]), h, where="mid", color=col, lw=2.0,
           label=lab)
    a.fill_between(0.5 * (bins[1:] + bins[:-1]), h, step="mid", color=col,
                   alpha=0.25)
a.axvline(A_HERTZ, color="k", lw=1.5, ls="--")
a.text(A_HERTZ * 1.06, a.get_ylim()[1] * 0.92, "a = %.2f mm" % A_HERTZ,
       fontsize=8.5)
a.set_xlabel("distance à l'axe d'impact (mm)")
a.set_ylabel("nombre de joints rompus")
a.set_title("(c) Le cisaillement est confiné sous le contact", fontsize=10.5)
a.legend(fontsize=8.5); a.grid(alpha=0.3)

fig.suptitle("Modes de rupture — %s : %d joints rompus, %.0f %% en traction, "
             "%.0f %% en cisaillement"
             % (RUN.replace("out_imp3d_", ""), nI + nII,
                100 * nI / (nI + nII), 100 * nII / (nI + nII)), fontsize=12.5)
fig.tight_layout()
out = os.path.join(HERE, "modes_rupture_%s.png" % RUN.replace("out_imp3d_", ""))
fig.savefig(out, dpi=140)
print("ecrit", out)
print("  mode I  : %4d | rayon median %.2f mm, profondeur mediane %.2f mm"
      % (nI, np.median(r[mode == 1]), np.median(120 - cen[mode == 1][:, 2])))
print("  mode II : %4d | rayon median %.2f mm, profondeur mediane %.2f mm"
      % (nII, np.median(r[mode == 2]), np.median(120 - cen[mode == 2][:, 2])))
print("  part de mode II sous le contact (r < a) : %.1f %%"
      % (100 * (mode[r < A_HERTZ] == 2).mean()))
print("  part de mode II au-dela de 3 mm         : %.1f %%"
      % (100 * (mode[r > 3.0] == 2).mean() if (r > 3.0).any() else 0.0))
