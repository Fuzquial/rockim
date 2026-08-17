# -*- coding: utf-8 -*-
"""Les branches de charge face a la solution ANALYTIQUE de Hertz — le seul juge
exterieur disponible.

    F(delta) = 4/3 E* sqrt(R) delta^(3/2)
    a(delta) = sqrt(R delta)          p0 = 3F/(2 pi a^2)

avec 1/E* = (1-nu1^2)/E1 + (1-nu2^2)/E2. Hypotheses : demi-espace et
quasi-statique — licites ici tant que delta << R et que a reste petit devant
le bloc, ce qui est le cas sur la branche de charge.

  python plot_hertz.py
"""
import csv, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
TWIN = 4e-5
E1, n1, E2, n2, Rs = 50e9, 0.25, 600e9, 0.22, 0.011
ES = 1 / ((1 - n1**2) / E1 + (1 - n2**2) / E2)

RUNS = [("p1_banc_mid", "82k uniforme — roche 4,74 mm",  4.74, "C0", "-"),
        ("p1_banc_fin", "259k uniforme — roche 3,29 mm", 3.29, "C2", "-"),
        ("p1_grad_15",  "gradué — roche 1,69 mm",        1.69, "C1", "--"),
        ("p1_grad_7",   "gradué — roche 0,80 mm",        0.80, "C3", "--")]


def hertzF(d):
    return (4 / 3) * ES * np.sqrt(Rs) * np.maximum(d, 0) ** 1.5


def hertzP0(d):
    F, a = hertzF(d), np.sqrt(Rs * np.maximum(d, 1e-12))
    return 3 * F / (2 * np.pi * a**2)


def load(run):
    D = os.path.join(ROOT, "out_elast_" + run)
    p = os.path.join(D, "history.csv")
    if not os.path.exists(p):
        return None
    rows = [r for r in csv.DictReader(open(p))
            if all(v not in (None, "") for v in r.values())]
    if len(rows) < 50 or float(rows[-1]["t"]) < 0.95 * TWIN:
        return None
    g = lambda k: np.array([float(r[k]) for r in rows])
    z = g("grpZ") * 1e3
    return z[0] - z, np.abs(g("grpFz")) * 1e-3, D


def peak_vm(D, pen, target=0.19):
    """von Mises max dans la roche, a la frame la plus proche de `target`."""
    fr = {int(r["frame"]): float(r["t"])
          for r in csv.DictReader(open(os.path.join(D, "frames.csv")))}
    hp = os.path.join(D, "history.csv")
    rows = [r for r in csv.DictReader(open(hp))
            if all(v not in (None, "") for v in r.values())]
    t = np.array([float(r["t"]) for r in rows])
    best, bd = None, 1e9
    for k, tk in fr.items():
        pk = float(np.interp(tk, t, pen))
        if abs(pk - target) < bd:
            best, bd = k, abs(pk - target)
    s = open(os.path.join(D, "fdem3d_%04d.vtu" % best), encoding="utf-8",
             errors="ignore").read()
    a = lambda n, d=float: np.fromstring(
        re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % n,
                  s, re.S).group(1), sep=" ", dtype=d)
    vm, ph = a("vonMises") * 1e-6, a("phase", int)
    return vm[ph == 0].max()


data = []
for run, lab, hr, c, ls in RUNS:
    d = load(run)
    if d is None:
        print("  (en cours ou absent : %s)" % run)
        continue
    data.append((run, lab, hr, c, ls) + d)

fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.9))

# (a) courbes vs Hertz
a = ax[0]
dd = np.linspace(0, 0.27, 300)
a.plot(dd, hertzF(dd * 1e-3) * 1e-3, "k", lw=2.4, label="Hertz (analytique)")
for run, lab, hr, c, ls, pen, fz, D in data:
    m = pen > 0
    a.plot(pen[m], fz[m], color=c, ls=ls, lw=1.6, label=lab)
a.set_xlabel("pénétration δ (mm)"); a.set_ylabel(r"$|F_z|$ (kN)")
a.set_title("(a) Branches de charge contre Hertz", fontsize=10)
a.legend(fontsize=7.5); a.grid(alpha=0.3); a.set_xlim(0, 0.27)

# (b) rapport mesure / Hertz
a = ax[1]
a.axhline(1.0, color="k", lw=2.0)
a.axhspan(0.95, 1.05, color="0.9", zorder=0)
for run, lab, hr, c, ls, pen, fz, D in data:
    m = (pen > 0.06) & (pen < 0.26)
    a.plot(pen[m], fz[m] / (hertzF(pen[m] * 1e-3) * 1e-3), color=c, ls=ls,
           lw=1.7, label=lab.split("—")[1].strip())
a.set_xlabel("pénétration δ (mm)"); a.set_ylabel("F mesurée / F Hertz")
a.set_ylim(0, 2.0)
a.set_title("(b) Écart à l'analytique\n(bande grise = ±5 %)", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3)

# (c) le pic de contrainte, lui, ne suit pas
a = ax[2]
hs, ps = [], []
for run, lab, hr, c, ls, pen, fz, D in data:
    try:
        p = peak_vm(D, pen)
    except Exception:
        continue
    hs.append(hr); ps.append(p)
    a.plot([hr], [p], "o", color=c, ms=9)
    a.annotate("%.0f" % p, (hr, p), xytext=(6, -3),
               textcoords="offset points", fontsize=8, color=c)
if hs:
    o = np.argsort(hs)[::-1]
    a.plot(np.array(hs)[o], np.array(ps)[o], "0.5", lw=1.2, zorder=0)
p0 = hertzP0(0.19e-3) * 1e-6
a.axhline(p0, color="k", lw=2.2, ls="--")
a.text(max(hs) if hs else 4, p0 * 1.04, "Hertz : $p_0$ = %.0f MPa" % p0,
       fontsize=9)
a.set_xlabel("taille d'élément dans la roche (mm)")
a.set_ylabel("von Mises max dans la roche (MPa)")
a.set_title("(c) Le PIC local, lui, est loin\n(δ ≈ 0,19 mm)", fontsize=10)
a.invert_xaxis(); a.grid(alpha=0.3); a.set_yscale("log")

fig.suptitle("Le contact fonctionne — la réponse globale suit Hertz ; c'est le "
             "pic local qui n'est pas résolu", fontsize=12.5)
fig.tight_layout()
out = os.path.join(HERE, "hertz_compare.png")
fig.savefig(out, dpi=140)
print("ecrit", out)
print("  E* = %.2f GPa | Hertz a delta=0,19 mm : F = %.1f kN, a = %.2f mm, "
      "p0 = %.0f MPa" % (ES * 1e-9, hertzF(0.19e-3) * 1e-3,
                         np.sqrt(Rs * 0.19e-3) * 1e3, p0))
for run, lab, hr, c, ls, pen, fz, D in data:
    r15 = float(np.interp(0.15, pen, fz)) / (hertzF(0.15e-3) * 1e-3)
    print("  %-32s F/Hertz a 0,15 mm = %.2f" % (lab, r15))
