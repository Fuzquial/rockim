#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_controle_run.py - PLANCHE DE CONTROLE D'UN RUN HYDRO EN COURS.
#
#   python bench_abuaisha/tools/fig_controle_run.py out_hfs_aniso [--sh 6.8]
#
# Quatre lectures independantes, toutes falsifiables, sur un run qui n'est pas
# fini. Elles ne mesurent pas un resultat : elles verifient que la physique
# part dans le bon sens.
#
#   (a) p(t) contre la loi de compressibilite (leur eq. 6) evaluee sur le
#       volume MESURE. Si le module hydro est juste, les deux se superposent.
#   (b) le VOLUME de cavite. Il doit DESCENDRE tant que p < contrainte in situ
#       (le trou se referme apres l'excavation), puis REMONTER. Une cavite qui
#       ne remonte jamais = signe de chargement inverse.
#   (c) contrainte principale majeure autour du forage, rendu par ELEMENTS.
#   (d) marge a la rupture sigma_1 / ft. Les lobes doivent etre alignes sur
#       sigma'_H (axe x) : c'est la direction ou Kirsch predit l'amorcage,
#       sigma_theta(0) = (sH + sh) - 2(sH - sh) - p.
# ---------------------------------------------------------------------------
import argparse, csv, glob, io, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "tunnel_edz"))
from plot_tunnel_fields import read_vtu, complete  # noqa: E402

for f in ("CMU Serif", "Latin Modern Roman", "DejaVu Serif"):
    if any(f in x.name for x in matplotlib.font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = f
        break
plt.rcParams.update({"font.size": 10, "mathtext.fontset": "cm", "figure.dpi": 120})

CX, CY, RB = 4.0, 4.0, 0.05
KF, RHO0, Q, FT = 2.2e9, 1000.0, 0.02, 5.0e6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--sh", type=float, default=6.8, help="sigma'_H [MPa]")
    ap.add_argument("--sv", type=float, default=4.6, help="sigma'_h [MPa]")
    ap.add_argument("--zoom", type=float, default=0.13)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    run = os.path.join(ROOT, a.run) if not os.path.isabs(a.run) else a.run

    h = list(csv.DictReader(open(os.path.join(run, "history.csv"))))
    t = np.array([float(r["t"]) for r in h]) * 1e3
    p = np.array([float(r["hydroP"]) for r in h]) / 1e6
    V = np.array([float(r["hydroVol"]) for r in h])
    nb = np.array([int(r["nBroken"]) for r in h])
    V0 = V[0]
    pcib = -a.sh + 3.0 * a.sv + FT / 1e6      # leur eq. 10, compression < 0

    fig = plt.figure(figsize=(11.5, 8.2))
    gs = fig.add_gridspec(2, 2, hspace=0.30, wspace=0.24)

    # --- (a) pression -------------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    # eq. 6 : m est la masse TOTALE, remplissage initial compris. Le
    # solveur part de hydroMass_ = V0 rho0 ; l'oublier fait un log(0).
    m = RHO0 * (V0 + Q * (t * 1e-3))
    pth = KF * np.log(np.maximum(m / (V * RHO0), 1e-300)) / 1e6
    ax.plot(t, p, lw=2.0, label="mesure")
    ax.plot(t, pth, "--", lw=1.2, color="crimson",
            label=r"$K_f\ln(m/V\rho_0)$ (leur eq. 6)")
    ax.axhline(pcib, color="0.4", ls=":", lw=1.0)
    ax.text(t[1], pcib + 0.3, "cible %.1f MPa (leur eq. 10)" % pcib,
            color="0.3", fontsize=8.5)
    ip = int(np.argmax(p))
    if nb[-1] > 0:
        ax.plot(t[ip], p[ip], "o", ms=6, mfc="none", mec="k", mew=1.4)
        ax.annotate("pic %.2f MPa" % p[ip], (t[ip], p[ip]),
                    textcoords="offset points", xytext=(-64, 6), fontsize=9)
    ax.set_xlabel("temps [ms]"); ax.set_ylabel("pression de puits [MPa]")
    ax.set_title("(a) la pompe et la compressibilite", loc="left")
    ax.legend(fontsize=8.5, loc="lower right"); ax.grid(alpha=0.25)

    # --- (b) volume ---------------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(t, (V / V0 - 1.0) * 100.0, lw=2.0, color="teal")
    ax.axhline(0.0, color="0.4", lw=0.8)
    ic = np.argmin(np.abs(p - a.sv))
    ax.axvline(t[ic], color="0.6", ls=":", lw=1.0)
    ax.text(t[ic], ax.get_ylim()[1] * 0.72,
            r"  $p = \sigma'_h$" + " = %.1f MPa" % a.sv, fontsize=8.5, color="0.3")
    ax.set_xlabel("temps [ms]")
    ax.set_ylabel(r"variation du volume de cavite  $V/V_0-1$  [%]")
    ax.set_title("(b) la cavite se referme, puis S'OUVRE", loc="left")
    ax.grid(alpha=0.25)

    # --- champs -------------------------------------------------------------
    fs = [f for f in sorted(glob.glob(os.path.join(run, "fdem_[0-9]*.vtu")))
          if complete(f)]
    P, C, F = read_vtu(fs[-1], ["sigmaXX", "sigmaYY", "sigmaXY"])
    k = int(os.path.basename(fs[-1]).split("_")[1].split(".")[0])
    fr = list(csv.DictReader(open(os.path.join(run, "frames.csv"))))
    tf = float(fr[k]["t"]) * 1e3
    sxx, syy, sxy = F["sigmaXX"], F["sigmaYY"], F["sigmaXY"]
    hm, r = 0.5 * (sxx + syy), np.hypot(0.5 * (sxx - syy), sxy)
    s1 = hm + r                                   # principale MAJEURE
    cen = P[C].mean(axis=1)
    sel = np.where((np.abs(cen[:, 0] - CX) < a.zoom) &
                   (np.abs(cen[:, 1] - CY) < a.zoom))[0]
    verts = [(P[C[i]] - [CX, CY]) * 1e3 for i in sel]

    # --- la FISSURE : joints endommages, rendu par SEGMENTS ----------------
    jf = [f for f in sorted(glob.glob(os.path.join(run, "fdem_joints_[0-9]*.vtu")))
          if complete(f)]
    JD = JS = None
    if jf:
        txt = io.open(jf[-1], errors="ignore").read()

        def jarr(nm):
            m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % nm, txt, re.S)
            return np.fromstring(m.group(1), sep=" ") if m else None
        JP = np.fromstring(re.search(
            r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>", txt, re.S).group(1),
            sep=" ").reshape(-1, 3)[:, :2]
        JC = np.fromstring(re.search(
            r'Name="connectivity"[^>]*>(.*?)</DataArray>', txt, re.S).group(1),
            sep=" ").astype(int).reshape(-1, 2)
        JD = jarr("damage")
        m = JD > 0.01
        if m.any():
            JS = (JP[JC[m]] - [CX, CY]) * 1e3
            JD = JD[m]
        else:
            JS = None

    for j, (val, ttl, cm, lab) in enumerate((
            (s1[sel] / 1e6, r"(c) contrainte principale majeure $\sigma_1$",
             "RdBu_r", r"$\sigma_1$  [MPa]"),
            (np.clip(s1[sel] / FT, 0, 1),
             r"(d) marge a la rupture  $\sigma_1/f_t$", "inferno", None))):
        if j == 1 and JS is not None:
            ax = fig.add_subplot(gs[1, 1])
            lc = LineCollection(JS, array=JD, cmap="inferno_r", lw=2.0,
                                clim=(0, 1))
            ax.add_collection(lc)
            th = np.linspace(0, 2 * np.pi, 200)
            ax.plot(RB * 1e3 * np.cos(th), RB * 1e3 * np.sin(th), "-",
                    color="0.5", lw=0.9)
            ext = max(70.0, 1.25 * np.abs(JS).max())
            ax.set_xlim(-ext, ext); ax.set_ylim(-ext, ext)
            ax.set_aspect("equal")
            ax.set_xlabel(r"x  [mm]   ($\sigma'_H$ = %.1f MPa)" % a.sh)
            ax.set_ylabel(r"y  [mm]   ($\sigma'_h$ = %.1f MPa)" % a.sv)
            ax.set_title(r"(d) LA FISSURE : joints endommages   $t$ = %.2f ms"
                         % tf, loc="left")
            ax.grid(alpha=0.2)
            cb = fig.colorbar(lc, ax=ax, fraction=0.046, pad=0.03)
            cb.set_label("endommagement du joint  $D$")
            msg = ("%d joints rompus" + chr(10) +
                   "%d faces mouillees (105 au depart)") % (
                       nb[-1], int(float(h[-1]["hydroNWet"])))
            ax.text(0.03, 0.03, msg,
                    transform=ax.transAxes, fontsize=9, va="bottom",
                    bbox=dict(fc="w", ec="0.7", alpha=0.85, pad=3))
            continue
        ax = fig.add_subplot(gs[1, j])
        v = np.max(np.abs(val)) if j == 0 else 1.0
        pc = PolyCollection(verts, array=val, cmap=cm, lw=0,
                            clim=(-v, v) if j == 0 else (0, 1))
        ax.add_collection(pc)
        th = np.linspace(0, 2 * np.pi, 200)
        ax.plot(RB * 1e3 * np.cos(th), RB * 1e3 * np.sin(th), "k-", lw=0.8)
        ax.set_xlim(-a.zoom * 1e3, a.zoom * 1e3)
        ax.set_ylim(-a.zoom * 1e3, a.zoom * 1e3)
        ax.set_aspect("equal")
        ax.set_xlabel(r"x  [mm]   ($\sigma'_H$ = %.1f MPa)" % a.sh)
        ax.set_ylabel(r"y  [mm]   ($\sigma'_h$ = %.1f MPa)" % a.sv)
        ax.set_title(ttl + "   $t$ = %.2f ms" % tf, loc="left")
        cb = fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.03)
        if lab:
            cb.set_label(lab)
        else:
            cb.set_label(r"$\sigma_1/f_t$   (1 = amorcage)")

    fig.suptitle("%s  -  t = %.2f / 4,00 ms,  p = %.2f MPa,  %d joints rompus"
                 % (os.path.basename(run), t[-1], p[-1], nb[-1]),
                 fontsize=11.5, y=0.985)
    out = a.out or os.path.join(run + "_controle.png")
    fig.savefig(out, bbox_inches="tight")
    print("ecrit :", out)
    print("  t = %.3f ms, p = %.2f MPa, V/V0-1 = %+.4f %%, nBroken = %d"
          % (t[-1], p[-1], (V[-1] / V0 - 1) * 100, nb[-1]))
    print("  sigma_1 max dans la fenetre : %.2f MPa  (ft = %.1f)"
          % (s1[sel].max() / 1e6, FT / 1e6))
    i = sel[np.argmax(s1[sel])]
    ang = np.degrees(np.arctan2(cen[i, 1] - CY, cen[i, 0] - CX))
    print("  ... atteint a %.0f deg de l'axe x, r = %.1f mm"
          % (ang, np.hypot(cen[i, 0] - CX, cen[i, 1] - CY) * 1e3))


if __name__ == "__main__":
    main()
