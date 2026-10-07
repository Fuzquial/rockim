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
#
# Campagne du 13/09 (T2, COMPLEMENT_YANG §5) : les vitesses de Yang 2025 §4.1
# sont des PENTES des portions lineaires du deplacement (indentation : 10-90 %
# de l'enfoncement ; rebond : portion remontante apres le retournement) et la
# contrainte de reference est le pic de la PREMIERE onde a mi-bit. Les deux
# estimateurs (pente / max instantane) sont imprimes cote a cote avec leurs
# fenetres (tools/yang_estimators.py) ; la fenetre de pente est hachuree sur
# (c) et la droite ajustee tracee, la 1re onde est hachuree sur (a). Les
# colonnes historiques du tableau (max instantanes) sont conservees.
# ---------------------------------------------------------------------------
import argparse
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import yang_estimators as ye  # noqa: E402

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


def _fmt_slope(d):
    """'6.858 m/s [64.8-169.4 us]' ou 'n.m.' (non mesurable) pour le tableau."""
    return ("%6.3f m/s [%5.1f-%5.1f us]" % (d["v"], d["t0"] * 1e6, d["t1"] * 1e6)
            if d["ok"] else "  n.m.  [%s]" % d["reason"][:34])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--stem", default="results/fig/kinetics")
    ap.add_argument("--body", default="insert",
                    help="corps dont le deplacement donne la pente 10-90 %% (insert | bit)")
    a = ap.parse_args()
    fig, AX = plt.subplots(1, 3, figsize=(16, 5.2))
    cols = ["#1f4e79", "#b22222", "#2e7d32", "#e08a00"]
    rows = []
    yang_rows = []
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
        # --- estimateurs de Yang (T2) : pentes + 1re onde, fenetres tracees ---
        res = ye.kinetics(ye.load_history(run), bodies=(a.body, "bit", "insert"))
        yang_rows.append((lab, res))
        d = res["bodies"].get(a.body)
        if d is not None:
            ind, reb = d["ind"], d["reb"]
            tt = res["t"] * 1e6
            pp = d["p"] * 1e3
            if ind["ok"]:
                AX[2].axvspan(ind["t0"] * 1e6, ind["t1"] * 1e6, color=c, alpha=0.08, lw=0)
                sl = ind["v"] * 1e-3            # mm/us
                tw = tt[ind["i0"]:ind["i1"] + 1]
                pw = pp[ind["i0"]:ind["i1"] + 1]
                b = pw.mean() - sl * tw.mean()
                AX[2].plot(tw, sl * tw + b, color=c, lw=0.9, ls="--",
                           label="%s : pente 10-90 %% (%s) = %.2f m/s" % (lab, a.body, ind["v"]))
                AX[1].axhline(ind["v"], color=c, lw=0.8, ls="--",
                              label="%s : pente 10-90 %% (%s) = %.2f m/s" % (lab, a.body, ind["v"]))
            if reb["ok"]:
                AX[2].axvspan(reb["t0"] * 1e6, reb["t1"] * 1e6, color=c, alpha=0.08, lw=0,
                              hatch="//")
                sl = -reb["v"] * 1e-3
                tw = tt[reb["i0"]:reb["i1"] + 1]
                pw = pp[reb["i0"]:reb["i1"] + 1]
                b = pw.mean() - sl * tw.mean()
                AX[2].plot(tw, sl * tw + b, color=c, lw=0.9, ls="-.",
                           label="%s : pente de rebond (%s) = %.2f m/s" % (lab, a.body, reb["v"]))
                AX[1].axhline(-reb["v"], color=c, lw=0.8, ls="-.",
                              label="%s : pente de rebond (%s) = %.2f m/s" % (lab, a.body, reb["v"]))
        g = res["gauge"]
        if g is not None and g["ok"]:
            AX[0].axvspan(g["t_on"] * 1e6, g["t_off"] * 1e6, color=c, alpha=0.08, lw=0)
            AX[0].plot([g["t"] * 1e6], [g["sig"] / 1e6], "o", color=c, ms=4,
                       label="%s : pic 1re onde %.0f MPa" % (lab, g["sig"] / 1e6))
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
    for A, loc in zip(AX, ("center right", "best", "lower right")):
        A.axhline(0, color="k", lw=0.5)
        A.legend(frameon=False, fontsize=7.5, loc=loc)
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
    # --- T2 : les estimateurs de Yang 2025 par. 4.1 (pentes, 1re onde) cote a
    # cote avec les max instantanes historiques, fenetres comprises ----------
    print()
    print("Estimateurs de Yang (pente des portions lineaires du deplacement, pic de la 1re onde)"
          " contre les max instantanes ; n.m. = non mesurable dans l'enregistrement")
    print("%-24s %-6s %-30s %-20s %-30s %-24s" % (
        "run", "corps", "v_ind pente 10-90 % [fenetre]", "v_ind max [date]",
        "v_reb pente [fenetre]", "v_reb max [date]"))
    for lab, res in yang_rows:
        for b, d in res["bodies"].items():
            inst = d["inst"]
            s_iind = ("%6.3f m/s [%5.1f us]" % (inst["v_ind"], inst["t_ind"] * 1e6)
                      if inst else "-")
            s_ireb = ("%6.3f m/s [%5.1f us]" % (inst["v_reb"], inst["t_reb"] * 1e6)
                      if inst and inst["ok_reb"] else
                      "  n.m.  [vz fin %+.2f m/s]" % (inst["v_end"] if inst else float("nan")))
            print("%-24s %-6s %-30s %-20s %-30s %-24s" % (
                lab[:24], b, _fmt_slope(d["ind"]), s_iind, _fmt_slope(d["reb"]), s_ireb))
        g = res["gauge"]
        if g is not None and g["ok"]:
            print("%-24s jauge  pic 1re onde %6.2f MPa a %5.1f us [onde %5.1f-%5.1f us%s] ; "
                  "max global %6.2f MPa a %5.1f us%s" % (
                      lab[:24], g["sig"] / 1e6, g["t"] * 1e6, g["t_on"] * 1e6, g["t_off"] * 1e6,
                      "" if g["finished"] else ", non retombee",
                      g["sig_max"] / 1e6, g["t_max"] * 1e6,
                      " (= 1re onde)" if g["t_max"] == g["t"] else " (onde ULTERIEURE)"))
    print("%-24s %-6s %-30s %-20s %-30s %-24s" % (
        "Yang 9 m/s (par. 4.1)", "bit", " 5.62 m/s [pente lineaire]", "-",
        " 4.65 m/s [pente, apres 450 us]", "-"))
    print("ecrit : %s.pdf / .png" % a.stem)


if __name__ == "__main__":
    main()
