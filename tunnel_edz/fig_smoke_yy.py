#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# fig_smoke_yy.py — planche des deux smokes de percussion 3D du 2026-08-19
# (out_smoke_yy et out_smoke_yy_iso), qui valident le portage 3D de la
# viscosite de Yan et du DIF de Yang.
#
#   python tunnel_edz/fig_smoke_yy.py
#
# LIMITE ASSUMEE : la part visqueuse n a PAS de colonne dans history.csv (le
# poste eVisc n a pas ete ajoute, decision D5 du plan de portage). Sa valeur
# n est donc connue qu en fin de run, par le resume. Le panneau (c) porte ce
# chiffre en annotation, pas en courbe.
# ---------------------------------------------------------------------------
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "axes.labelsize": 9.5, "axes.titlesize": 9.5,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5, "legend.fontsize": 8,
})
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def load(run):
    return np.genfromtxt(os.path.join(ROOT, run, "history.csv"),
                         delimiter=",", names=True, invalid_raise=False)


def main():
    a = load("out_smoke_yy")          # viscousInInsertion = 1 (defaut)
    b = load("out_smoke_yy_iso")      # viscousInInsertion = 0
    t = a["t"] * 1e6                  # us

    fig, ax = plt.subplots(1, 3, figsize=(10.2, 3.3))

    # ---- (a) signature de percussion -----------------------------------
    p = ax[0]
    p.plot(t, a["grpFz"] / 1e3, color="#0B4F9E", lw=1.3)
    p.set_xlabel("temps [µs]")
    p.set_ylabel("force axiale sur l'insert  [kN]", color="#0B4F9E")
    p.tick_params(axis="y", labelcolor="#0B4F9E")
    q = p.twinx()
    q.plot(t, -a["grpVz"], color="#C8342B", lw=1.3, ls="--")
    q.set_ylabel("vitesse de l'insert  [m s$^{-1}$]", color="#C8342B")
    q.tick_params(axis="y", labelcolor="#C8342B")
    p.set_title("(a) l'insert décélère sous la force de contact", fontsize=9.5)

    # ---- (b) postes du bilan, cumules ----------------------------------
    p = ax[1]
    for nom, cle, c, ls in (("éléments", "eEl", "#0B4F9E", "-"),
                            ("joints", "eJnt", "#1B8A3A", "-"),
                            ("contact", "eGc", "#8A5A00", "-"),
                            ("frottement", "eFric", "#8A5A00", ":"),
                            ("frontières", "eLys", "0.45", "-"),
                            ("Cundall", "eCund", "#C8342B", "-")):
        p.plot(t, -a[cle], color=c, ls=ls, lw=1.4, label=nom)
    p.axhline(0.0, color="0.8", lw=0.6)
    p.set_xlabel("temps [µs]")
    p.set_ylabel("énergie prélevée à l'énergie cinétique  [J]")
    p.legend(loc="upper left", framealpha=0.95, ncol=2)
    p.set_title("(b) Cundall reste plat à zéro : la viscosité\na remplacé l'amortisseur numérique",
                fontsize=9.5)

    # ---- (c) effet de viscousInInsertion, decomposition de fin de run ---
    # La difference n est PAS dans le travail preleve (identique au 4e chiffre)
    # mais dans la part que le bilan attribue a l elastique STOCKE : isole, sigG
    # redevient purement elastique et uEl cesse de compter 2 mu D comme stocke.
    p = ax[2]
    PREL = [0.404187, 0.404196]
    STOC = [0.368899, 0.324958]
    VISQ = [0.0358876, 0.0359059]
    xb = np.arange(2)
    dis = [PREL[k] - STOC[k] for k in range(2)]
    p.bar(xb, STOC, width=0.5, color="#0B4F9E", label="stocké élastique (uEl)")
    p.bar(xb, dis, width=0.5, bottom=STOC, color="#C8342B",
          label="reste du poste (dissipé)")
    for k in range(2):
        p.plot([xb[k] - 0.25, xb[k] + 0.25], [PREL[k] - VISQ[k]] * 2,
               color="k", lw=1.4, ls=":")
    p.text(1.28, PREL[1] - VISQ[1], "  visqueux mesuré\n  0,0359 J",
           fontsize=7.4, va="center")
    p.set_xticks(xb)
    p.set_xticklabels(["viscousInInsertion\n= 1 (fidèle à Yan)",
                       "viscousInInsertion\n= 0 (isolé)"], fontsize=8)
    p.set_ylabel("poste éléments en fin de run  [J]")
    p.set_ylim(0, 0.52)
    p.legend(loc="upper center", framealpha=0.95, fontsize=7.6)
    p.set_title("(c) le prélevé est identique ; ce qui bouge,\nc'est ce que le bilan appelle stocké",
                fontsize=9.5)

    for p in ax[:2]:
        p.grid(lw=0.4, color="0.92")
        p.set_axisbelow(True)
    for p in (ax[2],):
        p.grid(axis="y", lw=0.4, color="0.92")
        p.set_axisbelow(True)
    fig.suptitle("Smokes de percussion 3D — validation du portage de la viscosité de Yan "
                 "et du DIF de Yang (2026-08-19)", fontsize=10.5, y=1.02)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "fig_smoke_yy." + ext), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)
    print("ecrit : fig_smoke_yy.pdf / .png")


def fig_fp():
    """Force-penetration. L outil de tools/plot_force_penetration.py suppose un
    outil RIGIDE (toolY) ; ici l insert est un corps MAILLE pilote par groupe,
    donc la penetration se lit sur grpZ. Origine prise au PREMIER CONTACT
    (premier pas ou |grpFz| depasse 1 N), pas au debut du run : l insert part
    avec un jeu, et compter depuis t = 0 decalerait toute la courbe."""
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    for run, lab, c, ls in (("out_smoke_yy", "viscousInInsertion = 1", "#0B4F9E", "-"),
                            ("out_smoke_yy_iso", "viscousInInsertion = 0", "#C8342B", "--")):
        h = load(run)
        F = h["grpFz"] / 1e3                       # kN
        i0 = int(np.argmax(np.abs(h["grpFz"]) > 1.0))
        d = (h["grpZ"][i0] - h["grpZ"]) * 1e6      # um
        m = np.arange(len(F)) >= i0
        ax.plot(d[m], F[m], color=c, ls=ls, lw=1.5, label=lab)
        W = np.trapezoid(F[m] * 1e3, d[m] * 1e-6)
        print("  %-16s : delta_max %6.1f um, F_max %6.2f kN, aire %6.4f J"
              % (run, d[m].max(), F[m].max(), W))
    ax.set_xlabel(u"pénétration de l'insert  δ  [µm]")
    ax.set_ylabel("force axiale  [kN]")
    ax.set_title(u"Courbe force–pénétration, percussion 3D à 8 m s⁻¹\n"
                 u"(smoke : l'impact n'est pas terminé, pas de décharge)",
                 fontsize=9.5)
    ax.legend(loc="upper left", framealpha=0.95)
    ax.grid(lw=0.4, color="0.92")
    ax.set_axisbelow(True)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(HERE, "fig_smoke_fp." + ext), dpi=200)
    plt.close(fig)
    print("ecrit : fig_smoke_fp.pdf / .png")


if __name__ == "__main__":
    main()
    fig_fp()
