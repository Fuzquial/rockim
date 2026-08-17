# -*- coding: utf-8 -*-
"""ETAPE 1 — la reponse de CONTACT du banc converge-t-elle au maillage ?

Quatre maillages, impact elastique pur (ft = 1e12, aucune rupture possible),
fenetre 40 us : on compare la RAIDEUR, c'est-a-dire la force a penetration
donnee sur la branche de charge — grandeur convergente, contrairement au pic
dont l'instant se deplace d'un maillage a l'autre.

  python plot_convergence_contact.py
"""
import csv, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))

# (dossier, libelle, taille caracteristique dans la zone de contact [mm], couleur)
RUNS = [("out_elast_p1_banc_mid", "82k uniforme",      2.76, "C0"),
        ("out_elast_p1_grad_15",  "gradué 1,5 mm",     1.50, "C1"),
        ("out_elast_p1_banc_fin", "259k uniforme",     1.88, "C2"),
        ("out_elast_p1_grad_7",   "gradué 0,7 mm",     0.70, "C3")]
DELTAS = (0.05, 0.10, 0.15, 0.20)      # penetrations de mesure [mm]
TWIN = 4e-5                            # fenetre des runs elastiques [s]


def load(sub):
    p = os.path.join(ROOT, sub, "history.csv")
    if not os.path.exists(p):
        return None
    rows = [r for r in csv.DictReader(open(p))
            if all(v not in (None, "") for v in r.values())]
    if len(rows) < 50:
        return None
    # Un run EN COURS donne une branche de charge tronquee, donc une raideur
    # apparente fausse : on exige d'avoir atteint la fenetre T = 4e-5.
    if float(rows[-1]["t"]) < 0.95 * TWIN:
        print("  (en cours : %s, t = %.1f / %.0f us)"
              % (sub, float(rows[-1]["t"]) * 1e6, TWIN * 1e6))
        return None
    g = lambda k: np.array([float(r[k]) for r in rows])
    z = g("grpZ") * 1e3
    return dict(t=g("t") * 1e6, pen=z[0] - z, fz=np.abs(g("grpFz")) * 1e-3,
                nb=int(rows[-1]["nBroken"]), tend=g("t")[-1] * 1e6)


def nel(sub):
    p = os.path.join(ROOT, "run_" + sub.replace("out_", "") + ".log")
    if not os.path.exists(p):
        return None
    s = open(p, encoding="utf-8", errors="ignore").read().replace("\x00", "")
    m = re.search(r"(\d+) tets, (\d+) joints", s)
    return int(m.group(1)) if m else None


data = []
for sub, lab, h, c in RUNS:
    d = load(sub)
    if d is None:
        print("  (absent ou incomplet : %s)" % sub)
        continue
    if d["nb"]:
        print("  ATTENTION %s : %d joints rompus — le run n'est PAS elastique !"
              % (sub, d["nb"]))
    d.update(lab=lab, h=h, c=c, n=nel(sub))
    data.append(d)
if not data:
    raise SystemExit("aucun run exploitable")

fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8))

# (a) branches de charge superposees
a = ax[0]
for d in data:
    m = d["pen"] > 0
    a.plot(d["pen"][m], d["fz"][m], color=d["c"], lw=1.6,
           label="%s (%d k tets)" % (d["lab"], (d["n"] or 0) / 1000))
for x in DELTAS:
    a.axvline(x, color="0.85", lw=0.8, zorder=0)
a.set_xlabel("pénétration (mm)"); a.set_ylabel(r"$|F_z|$ (kN)")
a.set_title("(a) Branche de charge, élastique pur", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3)

# (b) convergence : F(delta fixe) en fonction de la taille d'element
a = ax[1]
tab = {}
for x in DELTAS:
    hs, fs = [], []
    for d in sorted(data, key=lambda q: -q["h"]):
        if d["pen"].max() < x:
            continue
        f = float(np.interp(x, d["pen"][d["pen"] > -1e-9], d["fz"][d["pen"] > -1e-9]))
        hs.append(d["h"]); fs.append(f)
    if len(hs) >= 2:
        a.plot(hs, fs, "o-", lw=1.5, ms=6, label=r"$\delta$ = %.2f mm" % x)
        tab[x] = (hs, fs)
a.set_xlabel("taille d'élément dans la zone de contact (mm)")
a.set_ylabel(r"$|F_z|$ à pénétration fixée (kN)")
a.set_title("(b) Converge-t-elle quand h → 0 ?", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3)
a.invert_xaxis()

# (c) verdict chiffre
a = ax[2]; a.axis("off")
lines = [("δ (mm)", "F au plus\ngrossier", "F au plus\nfin", "écart", "verdict")]
for x, (hs, fs) in tab.items():
    ecart = 100 * (fs[-1] / fs[0] - 1)
    verdict = ("convergé" if abs(ecart) < 5 else
               "limite" if abs(ecart) < 15 else "NON convergé")
    lines.append(("%.2f" % x, "%.1f" % fs[0], "%.1f" % fs[-1],
                  "%+.1f %%" % ecart, verdict))
tb = a.table(cellText=[list(r) for r in lines[1:]], colLabels=list(lines[0]),
             loc="center", cellLoc="center")
tb.auto_set_font_size(False); tb.set_fontsize(9); tb.scale(1, 2.0)
tb.auto_set_column_width(list(range(5)))
for j in range(5):
    tb[(0, j)].set_facecolor("#e8e8e8")
a.set_title("(c) Verdict  (élément 2,76 → 0,70 mm)", fontsize=10)

fig.suptitle("Convergence de la réponse de contact — impact élastique pur, "
             "insert R = 11 mm à 8 m/s", fontsize=12.5)
fig.tight_layout()
out = os.path.join(HERE, "convergence_contact.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
for d in sorted(data, key=lambda q: -q["h"]):
    print("  %-16s h = %.2f mm  %7d tets  t final %.1f us  pen max %.3f mm  "
          "joints rompus %d" % (d["lab"], d["h"], d["n"] or 0, d["tend"],
                                d["pen"].max(), d["nb"]))
