#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_fp.py — la courbe FORCE-PENETRATION de l'impact (leur fig. 7b).
#
#   python bench_impact/tools/fig_fp.py out_imp_stanne --stem bench_impact/fig_fp_stanne
#
# DEUX definitions de F, parce qu'elles ne racontent pas la meme chose :
#   (a) la force BIT-ROCHE, reconstruite par la quantite de mouvement du
#       corps bit+insert : F = m dv/dt - m g. Elle n'est valide qu'une fois
#       le piston SEPARE (le trait plein) — avant, l'acceleration melange la
#       pousse du piston et la reaction de la roche (pointille). C'est
#       l'analogue de leur fig. 7b ;
#   (b) la force de SECTION a mi-bit (la jauge, F = -sigma_zz A) : elle voit
#       d'abord l'onde incidente du piston (~140 kN AVANT le contact roche),
#       puis la force transmise en regime d'indentation. C'est leur fig. 9a
#       reportee en fonction de p.
# La penetration p est le deplacement du bit moins le jeu initial de 0,2 mm.
# ---------------------------------------------------------------------------
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from imp_lib import history

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.unicode_minus": False,
})

R_BIT, GAP = 0.015, 2.0e-4
RHO_S, RHO_C = 7850.0, 15250.0
R_INS, R_SHANK, H_INS = 0.00851, 0.00794, 0.0232


def masse_bit_insert():
    vb = np.pi * R_BIT ** 2 * (0.265 - H_INS)
    # insert : sphere + fut - recouvrement (calcul par tranches)
    z = np.linspace(-R_INS, H_INS - R_INS, 4000)
    dz = z[1] - z[0]
    r_sph = np.sqrt(np.maximum(R_INS ** 2 - z ** 2, 0.0))
    r_fut = np.where(z > 0, R_SHANK, 0.0)
    r = np.maximum(r_sph, r_fut)
    vi = float(np.pi * (r ** 2).sum() * dz)
    return RHO_S * vb + RHO_C * vi


def lisse(y, w):
    k = np.ones(w) / w
    return np.convolve(y, k, mode="same")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default="fig_fp")
    a = ap.parse_args()

    h = history(a.run)
    t = h["t"]
    p = (h["z_bit"][0] - h["z_bit"] - GAP) * 1e3
    m = masse_bit_insert()

    # (a) force de contact INSERT-ROCHE, exacte par Newton sur LES DEUX
    # corps : le piston ne touche que le bit, donc
    #   m_p dv_p/dt = F(bit->piston) - m_p g
    #   m_bi dv_bi/dt = F(roche->bit) + F(piston->bit) - m_bi g
    # et la somme elimine le contact interne (3e loi), ONDES COMPRISES —
    # les moyennes de corps rendent Newton exact pour les forces EXTERNES :
    #   F_roche = m_bi (dv_bi/dt + g) + m_p (dv_p/dt + g).
    # Controle integre : avant le contact insert-roche, la somme doit valoir
    # zero — c'est le segment p < 0 de la courbe.
    dt = float(np.median(np.diff(t)))
    w = max(3, int(10e-6 / dt) | 1)
    m_p = np.pi * 0.01325 ** 2 * 0.260 * RHO_S
    vb = lisse(h["vz_bit"], w)
    vp = lisse(h["vz_piston"], w)
    F = (m * (np.gradient(vb, t) + 9.81)
         + m_p * (np.gradient(vp, t) + 9.81)) / 1e3     # kN, compression > 0
    F = lisse(F, w)
    w2 = max(3, int(30e-6 / dt) | 1)
    Fm = lisse(F, w2)
    cut = len(t) - 2 * max(w, w2)
    sep = int(np.argmax(t > 3.5e-5))    # bord de convolution du depart

    # (b) jauge
    Fg = -h["szz_bit"] * np.pi * R_BIT ** 2 / 1e3

    fig, ax = plt.subplots(1, 2, figsize=(12.6, 5.4))
    A = ax[0]
    seg = np.array([p[sep:cut], F[sep:cut]]).T.reshape(-1, 1, 2)
    seg = np.concatenate([seg[:-1], seg[1:]], axis=1)
    lc = LineCollection(seg, cmap="viridis", array=t[sep:cut - 1] * 1e6,
                        linewidths=1.0, alpha=0.55)
    A.add_collection(lc)
    A.plot(p[sep:cut], Fm[sep:cut], color="#1f4e79", lw=2.2,
           label="moyenne (30 $\mu$s)")
    cb = fig.colorbar(lc, ax=A, pad=0.02)
    cb.set_label(r"temps [$\mu$s]")
    k = sep + int(np.argmax(Fm[sep:cut]))
    A.annotate("pic moyen %.0f kN à p = %.2f mm" % (Fm[k], p[k]), (p[k], Fm[k]),
               textcoords="offset points", xytext=(-8, 10), fontsize=10, ha="right")
    A.axhline(0, color="k", lw=0.5)
    A.axvline(0, color="#999", lw=0.6, ls=":")
    A.set_xlim(-0.25, max(p.max() * 1.12, 1.0))
    A.set_ylim(min(-5, F[sep:cut].min() * 1.15), F[sep:cut].max() * 1.25)
    A.set_xlabel("enfoncement du bit  [mm]")
    A.set_ylabel("force bit-roche  [kN]")
    A.set_title("(a)  $F = m\\,\\dot v$ du corps bit+insert "
                "(leur fig. 7b)", loc="left", fontsize=11)
    A.legend(frameon=False, fontsize=9, loc="lower right")

    B = ax[1]
    B.plot(p, Fg, color="#1f4e79", lw=1.1)
    B.axhline(0, color="k", lw=0.5)
    B.axvline(0, color="#999", lw=0.6, ls=":")
    B.set_xlabel("enfoncement du bit  [mm]")
    B.set_ylabel("force de section à mi-bit  [kN]")
    B.set_title("(b)  Jauge : onde incidente puis transmission",
                loc="left", fontsize=11)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(a.stem + "." + ext, dpi=170)
    print("écrit : %s | m(bit+insert) = %.3f kg | pic contact %.0f kN à "
          "p = %.2f mm | p_max %.2f mm"
          % (a.stem, m, Fm[k], p[k], p.max()))


if __name__ == "__main__":
    main()
