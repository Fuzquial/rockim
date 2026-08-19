#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# plot_sweep.py — les courbes du balayage in situ contre celles de l'article.
#
#   python tunnel_edz/plot_sweep.py [--out fichier.png]
#
# Quatre panneaux, calqués sur leurs figures :
#   (a) rayon de l'EDZ contre sigma_0            (leur fig. 12f)
#   (b) déplacement maximal contre sigma_0       (leur fig. 14a)
#   (c) nombre de fissures par mode              (leur fig. 13a)
#   (d) longueur cumulée de fissures             (leur fig. 13b)
#
# Deux estimateurs du rayon d'EDZ sont tracés, et c'est délibéré : le p95 des
# distances au centre et le maximum. L'article trace un cercle enveloppe sur
# ses figures 12a-e sans dire comment il est construit — montrer les deux
# évite de choisir après coup celui qui arrange.
#
# Les cas absents (runs pas encore terminés) sont simplement sautés.
# ---------------------------------------------------------------------------
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "tools"))
from edz_metrics import compute  # noqa: E402
from wall_convergence import convergence  # noqa: E402

# (sigma0 [MPa], dossier de sortie)
CASES = [(3.0, "out_tun_s3"), (4.0, "out_tun_s4"), (5.0, "out_tun_ref_iso"),
         (6.0, "out_tun_s6"), (7.0, "out_tun_s7")]

# Valeurs publiées — Wang et al. (2024), fig. 12f et 14a
ART_S = np.array([3.0, 4.0, 5.0, 6.0, 7.0])
ART_EDZ = np.array([11.1, 16.0, 19.0, 22.0, 22.8])
ART_U = np.array([0.244, 0.278, 0.347, 0.393, 0.451])

C_ART, C_P95, C_MAX = "#444444", "#2C6FB5", "#C8342B"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "sweep_insitu.png"))
    a = ap.parse_args()

    s, edz95, edzmx, umax, nt, nm, ns, lg = ([] for _ in range(8))
    cmean, cp90 = [], []
    for sig, run in CASES:
        path = os.path.join(HERE, "..", run)
        if not os.path.exists(os.path.join(path, "fdem_final_joints.csv")):
            print(f"  (saute {run} : pas encore termine)")
            continue
        r = compute(path)
        cv = convergence(path)
        s.append(sig)
        edz95.append(r["edz_radius_p95_m"])
        edzmx.append(r["edz_radius_max_m"])
        umax.append(r.get("u_max_m", np.nan))
        cmean.append(cv[0]); cp90.append(cv[1])
        nt.append(r["tensile"]); nm.append(r["mixed"]); ns.append(r["shear"])
        lg.append(r["crack_length_m"])
        print(f"  {sig:g} MPa : EDZ p95 {edz95[-1]:5.2f} / max {edzmx[-1]:5.2f} m,"
              f" u max {umax[-1]:.3f} / paroi moy {cmean[-1]:.4f} m, "
              f"{nt[-1]+nm[-1]+ns[-1]} fissures")
    s = np.array(s)
    if len(s) < 2:
        raise SystemExit("moins de deux cas termines")

    fig, ax = plt.subplots(2, 2, figsize=(13.5, 9.6))

    p = ax[0, 0]
    p.plot(ART_S, ART_EDZ, "s--", color=C_ART, lw=1.6, ms=7,
           label="Wang et al. 2024")
    p.plot(s, edzmx, "o-", color=C_MAX, lw=2, ms=7, label="rockim — max")
    p.plot(s, edz95, "^-", color=C_P95, lw=2, ms=7, label="rockim — p95")
    p.set_xlabel(r"contrainte in situ $\sigma_0$ [MPa]")
    p.set_ylabel("rayon de l'EDZ [m]")
    p.set_title("(a) etendue de la zone endommagee", fontsize=11)
    p.legend(fontsize=9)
    p.grid(alpha=0.3)

    p = ax[0, 1]
    p.plot(ART_S, ART_U, "s--", color=C_ART, lw=1.6, ms=7,
           label="Wang et al. 2024 (max de champ)")
    p.plot(s, umax, "o-", color=C_MAX, lw=2, ms=7,
           label="rockim — max de champ (fragile)")
    p.plot(s, cp90, "^-", color=C_P95, lw=2, ms=7, label="rockim — paroi p90")
    p.plot(s, cmean, "v-", color="#6A3D9A", lw=2.4, ms=7,
           label="rockim — paroi, moyenne (robuste)")
    p.axvline(3.94, color="#B8860B", ls=":", lw=1.6)
    p.text(3.99, 0.02, "UCS mesuree 3,94 MPa", rotation=90, fontsize=8,
           color="#B8860B", va="bottom")
    p.set_xlabel(r"contrainte in situ $\sigma_0$ [MPa]")
    p.set_ylabel("convergence [m]")
    p.set_title("(b) convergence de la paroi — trois estimateurs", fontsize=11)
    p.legend(fontsize=8)
    p.grid(alpha=0.3)

    p = ax[1, 0]
    w = 0.25
    p.bar(s - w, nt, w, color="#2E9E4F", label="traction")
    p.bar(s, nm, w, color="#2C6FB5", label="mixte")
    p.bar(s + w, ns, w, color="#C8342B", label="cisaillement")
    p.plot(s, np.array(nt) + np.array(nm) + np.array(ns), "ko-", lw=1.5, ms=5,
           label="total")
    p.set_xlabel(r"$\sigma_0$ [MPa]")
    p.set_ylabel("nombre de fissures")
    p.set_title("(c) repartition par mode de rupture", fontsize=11)
    p.legend(fontsize=9)
    p.grid(alpha=0.3, axis="y")

    p = ax[1, 1]
    p.plot(s, lg, "o-", color="#6A3D9A", lw=2, ms=7)
    p.set_xlabel(r"$\sigma_0$ [MPa]")
    p.set_ylabel("longueur cumulee de fissures [m]")
    p.set_title("(d) longueur totale de fissuration", fontsize=11)
    p.grid(alpha=0.3)

    fig.suptitle("Balayage de la contrainte in situ, lambda = 1 — maillage "
                 "isotrope, 106 298 triangles", fontsize=13)
    fig.tight_layout()
    fig.savefig(a.out, dpi=165)
    print("ecrit :", a.out)


if __name__ == "__main__":
    main()
