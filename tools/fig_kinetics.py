#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# fig_kinetics.py — la CINETIQUE de l'impact contre l'article (Yang et al. 2026,
# fig. 8 schema + fig. 9 a 9 m/s) : contrainte de jauge, vitesse du bit,
# enfoncement de l'insert, contre le temps. Les reperes de Yang a 9 m/s :
# contrainte max ~160 MPa (onde 1D rho c v/2 = 178) ; vitesse d'indentation
# 5,62 m/s ; enfoncement max ~1,0 mm (exp. 0,6-1,1) ; retournement ~255-300 us ;
# rebond 4,65 m/s (lu apres 450 us) ; impulsion du piston ~100 us (2L/c).
#
#   python tools/fig_kinetics.py out_A[:label] ... [--stem results/fig/kin]
# ---------------------------------------------------------------------------
import argparse
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
})

YANG = dict(sig=160.0, vind=5.62, depth=1.0, tturn=255.0, vreb=4.65, pulse=100.0)


def load(run):
    with open(os.path.join(run, "history.csv"), newline="") as f:
        r = [x for x in csv.DictReader(f) if x]
    t = np.array([float(x["t"]) for x in r]) * 1e6
    vz = np.array([float(x["vz_bit"]) for x in r])
    vp = np.array([float(x["vz_piston"]) for x in r])
    zi = np.array([float(x["z_insert"]) for x in r])
    szz = np.array([float(x["szz_bit"]) for x in r]) / 1e6
    return t, vz, vp, (zi[0] - zi) * 1e3, -szz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--stem", default="results/fig/kinetics")
    a = ap.parse_args()
    fig, AX = plt.subplots(1, 3, figsize=(16, 5.2))
    cols = ["#1f4e79", "#b22222", "#2e7d32", "#e08a00"]
    rows = []
    for i, spec in enumerate(a.runs):
        run, _, lab = spec.partition(":")
        lab = lab or os.path.basename(run)
        t, vz, vp, p, sig = load(run)
        c = cols[i % len(cols)]
        AX[0].plot(t, sig, color=c, lw=1.3, label=lab)
        AX[1].plot(t, -vz, color=c, lw=1.4, label=lab + " (bit)")
        AX[1].plot(t, -vp, color=c, lw=0.8, ls=":", label=lab + " (piston)")
        AX[2].plot(t, p, color=c, lw=1.4, label=lab)
        i0 = int(np.argmax(-vz))
        it = int(np.argmax(p))
        # duree de l impulsion : jauge > 50 % du max
        above = np.where(sig > 0.5 * sig.max())[0]
        pulse = (t[above[-1]] - t[above[0]]) if len(above) > 1 else float("nan")
        rows.append((lab, sig.max(), t[int(np.argmax(sig))], pulse, -vz[i0], t[i0],
                     p.max(), t[it], -vz[-1], t[-1]))
    # reperes Yang
    AX[0].axhline(YANG["sig"], color="k", ls="--", lw=0.8, label="Yang 9 m/s : ~160 MPa")
    AX[1].axhline(YANG["vind"], color="k", ls="--", lw=0.8, label="Yang : indentation 5,62 m/s")
    AX[1].axhline(-YANG["vreb"], color="k", ls="-.", lw=0.8, label="Yang : rebond 4,65 m/s")
    AX[1].axvline(YANG["tturn"], color="#888", ls=":", lw=0.8, label="Yang : retournement ~255 µs")
    AX[2].axhline(YANG["depth"], color="k", ls="--", lw=0.8, label="Yang : ~1,0 mm (exp. 0,6-1,1)")
    AX[2].axvline(YANG["tturn"], color="#888", ls=":", lw=0.8)
    AX[0].set_title("(a)  Contrainte de jauge à mi-bit (leur fig. 8-9a)", loc="left", fontsize=11)
    AX[0].set_xlabel(r"temps [$\mu$s]"); AX[0].set_ylabel(r"$-\sigma_{zz}$ [MPa]")
    AX[1].set_title("(b)  Vitesse du bit et du piston, vers le bas > 0", loc="left", fontsize=11)
    AX[1].set_xlabel(r"temps [$\mu$s]"); AX[1].set_ylabel("$v$ [m/s]")
    AX[2].set_title("(c)  Enfoncement de l'insert (leur fig. 8-9b)", loc="left", fontsize=11)
    AX[2].set_xlabel(r"temps [$\mu$s]"); AX[2].set_ylabel("$p$ [mm]")
    for A in AX:
        A.axhline(0, color="k", lw=0.5)
        A.legend(frameon=False, fontsize=7.5)
    fig.suptitle("Cinétique de l'impact contre Yang et al. 2026 (Kuru Grey, piston 9 m/s)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=160)
    print("%-32s %8s %7s %8s %8s %7s %7s %7s %8s %6s" % (
        "run", "sig_max", "t_sig", "impuls.", "v_ind", "t_vind", "p_max", "t_pmax", "v_fin", "t_fin"))
    for r in rows:
        print("%-32s %6.0f MPa %5.0f us %6.0f us %6.2f m/s %5.0f us %5.2f mm %5.0f us %6.2f m/s %5.0f us"
              % r)
    print("Yang 9 m/s                       ~160 MPa    ~40 us   ~100 us   5.62 m/s  ~100 us  ~1.0 mm  ~255 us  -4.65 m/s (rebond, apres 450 us)")
    print("ecrit : %s.pdf / .png" % a.stem)


if __name__ == "__main__":
    main()
