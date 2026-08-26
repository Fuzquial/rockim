#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_fp_imperial.py — la courbe FORCE-PENETRATION de la replique (leur fig. 7b).
#
#   python bench_impact/tools/fig_fp_imperial.py out_imperial \
#          --stem bench_impact/fig_fp_imperial
#
# DEUX mesures qui ne portent PAS sur la meme grandeur — erreur de la
# premiere version de ce script, corrigee le 2026-08-27. La jauge donne la
# force de SECTION a mi-bit ; Newton donne la reaction de la ROCHE. Elles ne
# coincident qu a l equilibre quasi statique de l outil, or le transit de
# l onde dans le bit dure 48 us pour un evenement de ~100 : on n y est
# JAMAIS. Leur ecart n est donc pas un defaut, c est la dynamique de l onde.
# La jauge est l analogue de leur fig. 9a, Newton celui de leur fig. 7b.
#
#   (a) NEWTON sur tout ce qui n'est pas la roche. Le piston ne touche que
#       l'outil, donc la 3e loi elimine le contact interne :
#           F_roche = m_outil (dv_outil/dt + g) + m_piston (dv_piston/dt + g)
#       Exact ondes comprises, parce que les vitesses sont des MOYENNES DE
#       CORPS. Controle integre : avant le contact roche la somme doit valoir
#       zero — c'est le segment a penetration nulle.
#   (b) LA JAUGE a mi-bit : F = -sigma_zz A. Mesure directe, aucune masse en
#       jeu, mais elle est 24 cm AU-DESSUS de la pointe : elle voit l'onde
#       incidente avant le contact roche et porte les reflexions.
#
# ABSCISSE : la penetration de l'INSERT (z_insert), jamais celle du bit. Le
# bit se fait comprimer par le piston a son sommet pendant que l'insert n'a
# pas bouge — confondre les deux fait lire 0,10 mm de penetration la ou il y
# en a 0,004 (constate le 2026-08-26).
#
# La masse de l'outil est VERIFIEE par conservation de la quantite de
# mouvement dans la fenetre 25-60 us (apres l'impact du piston, avant le
# contact roche) : 1,2182 kg +- 0,0004, soit le bit seul — l'insert n'y bouge
# pas encore. Voir --m-outil.
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
    "font.size": 10,
})

BLEU, ROUGE, GRIS, VERT = "#1f4e79", "#b22222", "#9a9a9a", "#2e7d32"
R_BIT = 0.015
M_PISTON = 0.962595          # lu au datacheck
M_BIT = 1.21859              # lu au datacheck, confirme par la quantite de
                             # mouvement : 1,2182 +- 0,0004 kg sur 25-60 us
M_INSERT = 0.0641428         # lu au datacheck


def lisse(y, w):
    return np.convolve(y, np.ones(w) / w, mode="same")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--stem", default="fig_fp_imperial")
    ap.add_argument("--m-bit", type=float, default=M_BIT)
    ap.add_argument("--m-insert", type=float, default=M_INSERT)
    ap.add_argument("--m-piston", type=float, default=M_PISTON)
    ap.add_argument("--lisse-us", type=float, default=3.0)
    a = ap.parse_args()

    h = history(a.run)
    t = h["t"]
    tus = t * 1e6
    pen = (h["z_insert"][0] - h["z_insert"]) * 1e3          # mm
    dt = float(np.median(np.diff(t)))
    w = max(3, int(a.lisse_us * 1e-6 / dt) | 1)

    # SIGNE — piege verifie le 2026-08-27. L axe z est VERS LE HAUT, donc on
    # garde vz TEL QUEL. Pour un corps soumis a la gravite et a la reaction
    # de la roche :   m a_z = F_roche - m g   =>   F_roche = m (a_z + g).
    # Prendre la vitesse en « descente positive » tout en gardant le +g inverse
    # le terme d acceleration et fabrique une force de -65 kN pendant la phase
    # ou l insert accelere : c est la premiere version de cette figure.
    # CHAQUE CORPS AVEC SA PROPRE ACCELERATION. Erreur du 2026-08-27 :
    # appliquer la masse bit+insert a l acceleration du BIT SEUL. L insert
    # passe de 0 a 11 m/s en 12 us, soit ~9e5 m/s^2 ; sa contribution vaut
    # 0,064 x 9e5 = 58 kN a elle seule, l ordre EXACT du desaccord constate
    # avec la jauge. Tant que l onde n a pas fini de traverser le bit
    # (transit 48 us), bit et insert ont des histoires tres differentes et ne
    # peuvent pas etre traites comme un corps unique.
    vb = lisse(h["vz_bit"], w)
    vi = lisse(h["vz_insert"], w)
    vp = lisse(h["vz_piston"], w)
    F = (a.m_bit * (np.gradient(vb, t) + 9.81)
         + a.m_insert * (np.gradient(vi, t) + 9.81)
         + a.m_piston * (np.gradient(vp, t) + 9.81)) / 1e3  # kN, roche > 0
    F = lisse(F, w)
    Fg = -h["szz_bit"] * np.pi * R_BIT ** 2 / 1e3           # kN

    # Les bords de la convolution `same` sont FAUX (elle y moyenne du vide) :
    # une reconstruction par gradient y produisait -150 kN. On les coupe.
    lo, hi = 2 * w, len(t) - 2 * w
    if hi <= lo:
        raise SystemExit("serie trop courte pour ce lissage")
    sl = slice(lo, hi)

    fig = plt.figure(figsize=(13.6, 4.9))
    gs = fig.add_gridspec(1, 3, wspace=0.30, left=0.06, right=0.975,
                          top=0.80, bottom=0.15)
    fig.suptitle("Force-penetration — replique Imperial College", fontsize=13,
                 y=0.965)
    fig.text(0.5, 0.895,
             "ETAT A MI-VOL : $t$ = %.1f $\\mu$s, penetration %.3f mm sur "
             "~1,53 attendus (%.0f %%). Leur fig. 7b culmine vers 130 kN."
             % (tus[-1], pen.max(), 100 * pen.max() / 1.53),
             ha="center", fontsize=9.5, color="#444444", style="italic")

    # ---- (a) Newton -------------------------------------------------------
    A = fig.add_subplot(gs[0, 0])
    seg = np.array([pen[sl], F[sl]]).T.reshape(-1, 1, 2)
    seg = np.concatenate([seg[:-1], seg[1:]], axis=1)
    lc = LineCollection(seg, cmap="viridis", array=tus[sl][:-1],
                        linewidths=1.6)
    A.add_collection(lc)
    cb = fig.colorbar(lc, ax=A, pad=0.02)
    cb.set_label(r"temps  [$\mu$s]", fontsize=9)
    k = lo + int(np.argmax(F[sl]))
    A.annotate("pic %.0f kN\na %.3f mm" % (F[k], pen[k]), (pen[k], F[k]),
               xytext=(-8, 8), textcoords="offset points", ha="right",
               fontsize=9, color=ROUGE)
    A.axhline(0, color="k", lw=0.6)
    A.axvline(0, color=GRIS, lw=0.6, ls=":")
    A.set_xlim(-0.02, max(pen.max() * 1.15, 0.05))
    A.set_ylim(min(-5, F[sl].min() * 1.2), max(F[sl].max() * 1.25, 10))
    A.set_xlabel("penetration de l'insert  [mm]")
    A.set_ylabel("force outil-roche  [kN]")
    A.set_title("(a)  Reaction de la ROCHE (leur fig. 7b)\n"
                "par Newton sur chaque corps", loc="left", fontsize=10.5)

    # ---- (b) la jauge -----------------------------------------------------
    B = fig.add_subplot(gs[0, 1])
    B.plot(pen[sl], Fg[sl], color=BLEU, lw=1.3)
    B.axhline(0, color="k", lw=0.6)
    B.axvline(0, color=GRIS, lw=0.6, ls=":")
    B.set_xlabel("penetration de l'insert  [mm]")
    B.set_ylabel("force de section a mi-bit  [kN]")
    B.set_title("(b)  Force de SECTION a mi-bit\n(leur fig. 9a)",
                loc="left", fontsize=10.5)
    B.annotate("la jauge est 24 cm AU-DESSUS\nde la pointe : elle voit l'onde\n"
               "avant que la roche soit touchee",
               (0.96, 0.06), xycoords="axes fraction", ha="right",
               fontsize=7.5, color="#555555", style="italic")

    # ---- (c) les deux dans le temps ---------------------------------------
    C = fig.add_subplot(gs[0, 2])
    C.plot(tus[sl], F[sl], color=ROUGE, lw=1.5, label="Newton (outil+piston)")
    C.plot(tus[sl], Fg[sl], color=BLEU, lw=1.2, alpha=0.8, label="jauge mi-bit")
    C.axhline(0, color="k", lw=0.6)
    C.set_xlabel(r"temps  [$\mu$s]")
    C.set_ylabel("force  [kN]")
    C.legend(frameon=False, fontsize=8.5, loc="upper left")
    C2 = C.twinx()
    C2.plot(tus, pen, color=VERT, lw=1.2, ls="--")
    C2.set_ylabel("penetration de l'insert  [mm]", color=VERT)
    C2.tick_params(axis="y", labelcolor=VERT)
    C.set_title("(c)  Deux grandeurs DIFFERENTES : section\n"
                "a mi-bit contre reaction de la roche",
                loc="left", fontsize=10.5)

    fig.text(0.5, 0.028,
             "bit %.4f + insert %.4f + piston %.4f kg, CHACUN avec sa propre "
             "acceleration (masse du bit confirmee a 0,03 %% par la quantite "
             "de mouvement sur 25-60 $\\mu$s)  |  lissage %.0f $\\mu$s  |  "
             "circlip et platine NON comptes"
             % (a.m_bit, a.m_insert, a.m_piston, a.lisse_us),
             ha="center", fontsize=8.5, color="#444444")

    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (a.stem, ext), dpi=190)
    print("ecrit : %s.pdf et .png" % a.stem)
    print("  penetration max   : %.4f mm" % pen.max())
    print("  pic Newton        : %.1f kN a %.4f mm" % (F[k], pen[k]))
    print("  pic jauge         : %.1f kN" % Fg[sl].max())
    print("  cible publiee     : ~130 kN a ~1,53 mm")


if __name__ == "__main__":
    main()
