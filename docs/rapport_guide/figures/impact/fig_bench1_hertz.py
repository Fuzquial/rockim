# -*- coding: utf-8 -*-
"""Impact court d'un insert spherique maille (banc configs/fdem3d_bench1_insert.cfg), rejoue le
2026-10-06 avec le binaire f0209ef, contre la solution de Hertz.

Deux calculs de 50 us, OMP_NUM_THREADS = 1 (configs a cote des sorties, dans le scratchpad de la
session) : roche elastique (ft = c = 1e12 Pa, sans amortissement local) et roche fissurable (deck
d'origine). Seule modification de montage : le jeu insert-roche du maillage passe de 0,5 mm a 20 um
(translation rigide des noeuds de l'insert) pour ne pas payer 60 us de vol libre.

Usage : python3 fig_bench1_hertz.py <dossier_runs>   (contient out_elas/ et out_frac/ et le .msh)
"""
import csv, math, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

RUNS = sys.argv[1]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bench1_hertz")
plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9.5,
                     "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.5,
                     "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.4,
                     "legend.frameon": False, "legend.fontsize": 8})
VIRG = FuncFormatter(lambda x, p: ("%g" % x).replace(".", ","))
C = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#555555"]

# ---- geometrie : masse et rayon de l'insert lus dans le maillage --------------------------
L = open(os.path.join(RUNS, "bench1_insert_gap20um.msh")).read().split("\n")
i0 = L.index("$Nodes"); nn = int(L[i0 + 1])
X = {int(l.split()[0]): np.array([float(v) for v in l.split()[1:4]]) for l in L[i0 + 2:i0 + 2 + nn]}
e0 = L.index("$Elements"); ne = int(L[e0 + 1])
vol, zmin_i, zmax_r, ins_nodes = 0.0, 1e9, -1e9, set()
for l in L[e0 + 2:e0 + 2 + ne]:
    p = [int(v) for v in l.split()]
    if p[1] != 4:
        continue
    nt = p[2]; v = p[3 + nt:]
    if p[3] == 2:
        a, b, c, d = (X[k] for k in v)
        vol += abs(np.dot(b - a, np.cross(c - a, d - a))) / 6
        ins_nodes.update(v)
    else:
        zmax_r = max(zmax_r, max(X[k][2] for k in v))
P = np.array([X[k] for k in ins_nodes])
zmin_i = P[:, 2].min(); R = (P[:, 2].max() - zmin_i) / 2
GAP = zmin_i - zmax_r
RHO_I, E_I, NU_I = 14500.0, 600e9, 0.22
E_R, NU_R = 50e9, 0.25
V0 = 8.0
m = RHO_I * vol
Es = 1.0 / ((1 - NU_R ** 2) / E_R + (1 - NU_I ** 2) / E_I)
K = 4.0 / 3.0 * Es * math.sqrt(R)
dmax = (15 * m * V0 ** 2 / (16 * Es * math.sqrt(R))) ** 0.4
Fmax = K * dmax ** 1.5
tc = 2.9432 * dmax / V0
amax = math.sqrt(R * dmax)
p0 = 3 * Fmax / (2 * math.pi * amax ** 2)

# courbe de Hertz en temps (RK4 sur m d'' = -K d^1.5, d(0) = 0, d'(0) = V0)
def hertz(T, n=20000):
    dt = T / n; d, v = 0.0, V0; out = [(0.0, 0.0)]
    acc = lambda x: -K * max(x, 0) ** 1.5 / m
    for i in range(n):
        k1 = (v, acc(d)); k2 = (v + dt / 2 * k1[1], acc(d + dt / 2 * k1[0]))
        k3 = (v + dt / 2 * k2[1], acc(d + dt / 2 * k2[0])); k4 = (v + dt * k3[1], acc(d + dt * k3[0]))
        d += dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); v += dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        out.append(((i + 1) * dt, d))
    o = np.array(out); return o[:, 0], o[:, 1], K * np.clip(o[:, 1], 0, None) ** 1.5


def load(run):
    rows = [r for r in csv.DictReader(open(os.path.join(RUNS, run, "history.csv"))) if r and r.get("eLys")]
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def sm(y, n):
    return np.convolve(y, np.ones(n) / n, mode="same")


runs = {}
for name in ("out_elas", "out_frac"):
    if os.path.isfile(os.path.join(RUNS, name, "history.csv")):
        runs[name] = load(name)

t0 = GAP / V0
th, dh, Fh = hertz(min(tc, 60e-6))
print("insert : R %.3f mm, volume %.1f mm3, masse %.4f kg, jeu %.1f um" % (R * 1e3, vol * 1e9, m, GAP * 1e6))
print("Hertz : E* %.2f GPa, dmax %.1f um, Fmax %.2f kN, tc %.1f us, amax %.3f mm, p0 %.0f MPa"
      % (Es / 1e9, dmax * 1e6, Fmax / 1e3, tc * 1e6, amax * 1e3, p0 / 1e6))

fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.7))
lab = {"out_elas": "roche élastique", "out_frac": "roche fissurable"}
col = {"out_elas": C[0], "out_frac": C[1]}
for name, h in runs.items():
    t = h["t"]; z = h["grpZ"]
    d = (z[0] - z) - GAP                       # rapprochement apres fermeture du jeu
    F = h["Fc_rock_insert_z"]
    n = max(1, int(round(1e-6 / np.median(np.diff(t)))))
    Fest = sm(sm(m * np.gradient(h["grpVz"], t), n), n)
    ax[0].plot((t - t0) * 1e6, F / 1e3, color=col[name], label=lab[name])
    if name == "out_frac":
        ax[0].plot((t - t0) * 1e6, Fest / 1e3, color=C[4], ls="--", lw=0.9, label=r"$m\,\dot v$, roche fissurable")
    ax[1].plot(d * 1e3, F / 1e3, color=col[name], label=lab[name])
    ax[2].plot((t - t0) * 1e6, h["nBroken"], color=col[name], label=lab[name])
    k = d > 0.2 * dmax
    i = np.argmax(F)
    tt = (t - t0)
    Fh_t = np.interp(tt, th, Fh)
    print("%s : t fin %.1f us, F max %.2f kN a d = %.1f um (t - t0 = %.1f us) ; d max %.1f um ; F/Fhertz(t) a la fin %.3f ; "
          "nBroken fin %d ; impulsion mesuree %.4f N s, m dv %.4f N s ; ecart rms F - m dv/dt %.2f kN"
          % (name, t[-1] * 1e6, F.max() / 1e3, d[i] * 1e6, tt[i] * 1e6, d.max() * 1e6,
             F[-1] / max(Fh_t[-1], 1e-9), h["nBroken"][-1], np.trapezoid(F, t), m * (h["grpVz"][-1] - h["grpVz"][0]),
             np.sqrt(np.mean((Fest - F)[t > t0] ** 2)) / 1e3))
    on = np.where(F > 0.01 * Fmax)[0]
    if len(on):
        t_on, t_off = t[on[0]] - t0, t[on[-1]] - t0
        fin = "contact NON termine a la fin du calcul (F = %.2f kN)" % (F[-1] / 1e3) if on[-1] == len(F) - 1 else \
            "fin du contact a %.1f us, duree %.1f us (Hertz %.1f)" % (t_off * 1e6, (t_off - t_on) * 1e6, tc * 1e6)
        print("   debut du contact (F > 1 %% Fmax Hertz) a %.1f us apres t0 ; %s ; vz centroide fin %.3f m/s ; Hertz : F max a t_c/2 = %.1f us"
              % (t_on * 1e6, fin, h["grpVz"][-1], tc / 2 * 1e6))
    iload = np.arange(np.argmax(d) + 1)
    for dd in (0.25, 0.5, 0.75, 1.0):
        if dd * dmax > d.max():
            print("   d/dmax %.2f : non atteint (d max %.1f um)" % (dd, d.max() * 1e6)); continue
        j = iload[np.argmin(abs(d[iload] - dd * dmax))]
        print("   d/dmax %.2f : F = %.2f kN, Hertz %.2f kN, rapport %.3f" % (dd, F[j] / 1e3, K * (dd * dmax) ** 1.5 / 1e3,
                                                                       F[j] / (K * (dd * dmax) ** 1.5)))
    eb = {k: -h[k][-1] for k in ("eEl", "eJnt", "eGc", "eFric")}
    print("   postes fin :", {k: round(v, 3) for k, v in eb.items()})

ax[0].plot(th * 1e6, Fh / 1e3, color="black", lw=2.2, alpha=0.35, label="Hertz")
dd = np.linspace(0, dmax, 200)
ax[1].plot(dd * 1e3, K * dd ** 1.5 / 1e3, color="black", lw=2.2, alpha=0.35, label=r"Hertz $\frac{4}{3}E^*\sqrt{R}\,\delta^{3/2}$")
ax[0].set_xlabel(r"temps depuis le contact [$\mu$s]"); ax[0].set_ylabel(r"$F_z$ roche-insert [kN]")
ax[1].set_xlabel(r"rapprochement $\delta$ [mm]"); ax[1].set_ylabel(r"$F_z$ [kN]")
ax[2].set_xlabel(r"temps depuis le contact [$\mu$s]"); ax[2].set_ylabel("joints rompus")
ax[0].set_title("(a) force contre temps", loc="left", fontsize=9.5)
ax[1].set_title("(b) force contre rapprochement", loc="left", fontsize=9.5)
ax[2].set_title("(c) rupture", loc="left", fontsize=9.5)
ax[0].legend(loc="upper left", fontsize=7); ax[1].legend(loc="upper left", fontsize=7)
for a in ax:
    a.xaxis.set_major_formatter(VIRG); a.yaxis.set_major_formatter(VIRG)
fig.text(0.99, 0.005, "rejoué le 2026-10-06 (binaire f0209ef)", ha="right", fontsize=7, color="#555555")
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(OUT + ".pdf"); fig.savefig(OUT + ".png", dpi=200)
