# -*- coding: utf-8 -*-
"""Cas 1 — impact 3D homogene, panoplie complete (out_imp3d_homog).
Force-penetration, coupe verticale en von Mises, et carte des joints INSERES
par le critere adaptatif face a l'empreinte de contact (la question du run).

  python plot_imp3d_homog.py [nom_du_run]
"""
import csv, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
RUN = sys.argv[1] if len(sys.argv) > 1 else "out_imp3d_homog"
D = os.path.join(ROOT, RUN)
V0, ZTOP = 8.0, 0.12          # vitesse d'impact, sommet du bloc (m)


_CACHE = {}


def _txt(path):
    if path not in _CACHE:
        _CACHE[path] = open(path, encoding="utf-8", errors="ignore").read()
    return _CACHE[path]


def darray(path, name, ncomp=1, dtype=float):
    """Extrait un DataArray nomme d'un .vtu ASCII."""
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name,
                  _txt(path), re.S)
    if not m:
        return None
    v = np.fromstring(m.group(1), sep=" ", dtype=dtype)
    return v.reshape(-1, ncomp) if ncomp > 1 else v


def points(path):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>",
                  _txt(path), re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)


# ---- historique ----------------------------------------------------------
rows = [r for r in csv.DictReader(open(os.path.join(D, "history.csv")))
        if all(v not in (None, "") for v in r.values())]
g = lambda k: np.array([float(r[k]) for r in rows])
t = g("t") * 1e6
pen = (g("grpZ")[0] - g("grpZ")) * 1e3
fz = np.abs(g("grpFz")) * 1e-3
nb = g("nBroken")

# ---- derniere frame ------------------------------------------------------
fs = sorted(f for f in os.listdir(D) if re.fullmatch(r"fdem3d_\d{4}\.vtu", f))
# La frame a montrer est celle du PIC de force, pas la derniere : a t = T le
# rebond est fini et le bloc est decharge (contraintes ~ 0, image trompeuse).
ftimes = {int(r["frame"]): float(r["t"]) * 1e6
          for r in csv.DictReader(open(os.path.join(D, "frames.csv")))}
tpeak = t[int(np.argmax(fz))]
kbest = min(ftimes, key=lambda k: abs(ftimes[k] - tpeak))
last = os.path.join(D, "fdem3d_%04d.vtu" % kbest)
print("frame du pic : %04d a t = %.1f us (pic de force a %.1f us)"
      % (kbest, ftimes[kbest], tpeak))
pts = points(last)
conn = darray(last, "connectivity", dtype=int).reshape(-1, 4)
vm = darray(last, "vonMises")
phase = darray(last, "phase", dtype=int)

jl = last.replace("fdem3d_", "fdem3d_joints_")
jconn = darray(jl, "connectivity", dtype=int)
joff = darray(jl, "offsets", dtype=int)
bonded = darray(jl, "bonded", dtype=int)
dmg = darray(jl, "damage")
tbr = darray(jl, "tBreak")

cen = pts[conn].mean(axis=1)
rock = phase == 0

# Un tetra projete par UNE SEULE de ses faces ne pave pas le plan : mesure
# faite le 2026-08-17, 96,8 % de couverture seulement, les 3,2 % restants
# apparaissant en trous blancs. On projette donc les QUATRE faces (100 %).
FACES = ([0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3])


def tri_proj(idx, cols, vals=None):
    """Les 4 faces de chaque tetra idx, projetees sur les axes cols (en mm)."""
    polys, out = [], []
    for n, k in enumerate(idx):
        q = pts[conn[k]][:, cols] * 1e3
        for f in FACES:
            polys.append(q[f])
            if vals is not None:
                out.append(vals[n])
    return polys, np.array(out)

fig = plt.figure(figsize=(15.5, 8.4))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15])

# --- (a) force - penetration ---------------------------------------------
a = fig.add_subplot(gs[0, 0])
a.plot(pen, fz, "C0", lw=1.5)
i = int(np.argmax(fz))
a.plot(pen[i], fz[i], "o", color="C3", ms=6)
a.annotate("pic %.1f kN\nà %.3f mm" % (fz[i], pen[i]), (pen[i], fz[i]),
           xytext=(-72, -26), textcoords="offset points", fontsize=8, color="C3")
a.set_xlabel("pénétration (mm)"); a.set_ylabel(r"$|F_z|$ (kN)")
a.set_title("(a) Force – pénétration", fontsize=10)
a.grid(alpha=0.3)

# --- (b) force et penetration dans le temps ------------------------------
a = fig.add_subplot(gs[0, 1])
a.plot(t, fz, "C0", lw=1.3, label=r"$|F_z|$")
a.set_xlabel("temps (µs)"); a.set_ylabel(r"$|F_z|$ (kN)", color="C0")
a2 = a.twinx()
a2.plot(t, pen, "C2", lw=1.3)
a2.axhline(0, color="0.7", lw=0.8)
a2.set_ylabel("pénétration (mm)", color="C2")
a.set_title("(b) Contact, charge–décharge, rebond", fontsize=10)
a.grid(alpha=0.3)

# --- (c) joints rompus dans le temps -------------------------------------
a = fig.add_subplot(gs[0, 2])
a.step(t, nb, "C3", lw=1.6, where="post")
a.set_xlabel("temps (µs)"); a.set_ylabel("joints rompus")
a.set_title("(c) Fissuration : %d joints sur %d" % (nb[-1], len(bonded)),
            fontsize=10)
a.grid(alpha=0.3)

# --- (d) coupe verticale, von Mises --------------------------------------
a = fig.add_subplot(gs[1, 0])
band = rock & (np.abs(cen[:, 1] - 0.06) < 0.003)
kb = np.where(band)[0]
polys, vals = tri_proj(kb, [0, 2], vm[kb] * 1e-6)
pc = a.add_collection(PolyCollection(polys, array=vals,
                                     cmap="inferno", edgecolors="none"))
pc.set_clim(0, 40.0)   # MPa — 4 x ft, echelle commune (d)/(e)
cb = fig.colorbar(pc, ax=a, shrink=0.85); cb.set_label("von Mises (MPa)")
a.set_xlim(0, 120); a.set_ylim(60, 122); a.set_aspect("equal")
a.set_xlabel("x (mm)"); a.set_ylabel("z (mm)")
a.set_title("(d) Coupe y = 60 mm à t = %.0f µs (pic)\n"
            "roche : $f_t$ = 10 MPa" % ftimes[kbest], fontsize=10)

# --- (e) vue de dessus : von Mises en surface ----------------------------
a = fig.add_subplot(gs[1, 1])
surf = rock & (cen[:, 2] > ZTOP - 0.006)
ks = np.where(surf)[0]
polys, vals = tri_proj(ks, [0, 1], vm[ks] * 1e-6)
pc = a.add_collection(PolyCollection(polys, array=vals,
                                     cmap="inferno", edgecolors="none"))
pc.set_clim(0, 40.0)   # MPa — 4 x ft, echelle commune (d)/(e)
cb = fig.colorbar(pc, ax=a, shrink=0.85); cb.set_label("von Mises (MPa)")
a.set_xlim(0, 120); a.set_ylim(0, 120); a.set_aspect("equal")
a.set_xlabel("x (mm)"); a.set_ylabel("y (mm)")
a.set_title("(e) Vue de dessus, 6 mm sous la surface\n"
            "à t = %.0f µs" % ftimes[kbest], fontsize=10)

# --- (f) LA question : ou le critere a-t-il insere ? ---------------------
a = fig.add_subplot(gs[1, 2])
jc = [jconn[(0 if i == 0 else joff[i - 1]):joff[i]] for i in range(len(joff))]
ins = np.where(bonded == 0)[0]
brk = np.where(tbr >= 0)[0]
jcen = np.array([pts[jc[i]].mean(axis=0) for i in ins]) if len(ins) else np.empty((0, 3))
bcen = np.array([pts[jc[i]].mean(axis=0) for i in brk]) if len(brk) else np.empty((0, 3))
a.add_collection(PolyCollection(tri_proj(ks, [0, 1])[0],
                                facecolors="#ece8e0", edgecolors="0.8",
                                linewidths=0.1))
if len(jcen):
    a.scatter(jcen[:, 0] * 1e3, jcen[:, 1] * 1e3, s=18, c="C0",
              label="joints insérés (%d)" % len(ins), zorder=5)
if len(bcen):
    a.scatter(bcen[:, 0] * 1e3, bcen[:, 1] * 1e3, s=70, marker="x", c="crimson",
              lw=2, label="rompus (%d)" % len(brk), zorder=6)
a.add_patch(plt.Circle((60, 60), 11, fill=False, ec="C2", lw=1.6, ls="--"))
a.text(60, 74, "insert R = 11 mm", color="C2", fontsize=8, ha="center")
a.set_xlim(20, 100); a.set_ylim(20, 100); a.set_aspect("equal")
a.set_xlabel("x (mm)"); a.set_ylabel("y (mm)")
a.set_title("(f) Où le critère adaptatif a-t-il inséré ?", fontsize=10)
a.legend(fontsize=8, loc="upper right")

TITRE = {"out_imp3d_homog": "homogène",
         "out_imp3d_weib": "hétérogène (Weibull m = 6, facteurs 0,17–1,63)"}
fig.suptitle("Impact 3D %s, panoplie complète — bloc 120³ mm, insert "
             "R = 11 mm à 8 m/s, 82k tets, 6 min de calcul"
             % TITRE.get(RUN, RUN), fontsize=12.5)
fig.tight_layout()
out = os.path.join(HERE, "imp3d_%s.png" % RUN.replace("out_imp3d_", ""))
fig.savefig(out, dpi=140)
print("ecrit", out)
print("joints inseres : %d / %d (%.4f %%) | rompus : %d | damage max : %.3f"
      % (len(ins), len(bonded), 100 * len(ins) / len(bonded), len(brk),
         np.nanmax(dmg) if dmg is not None else -1))
