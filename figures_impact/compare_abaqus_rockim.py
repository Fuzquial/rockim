# -*- coding: utf-8 -*-
"""Comparaison de proprete des F-delta : FDEM Abaqus 2D (indenter_voronoi_rod,
tige+Lysmer, 2026-06) vs rockim 3D (banc moyen, contact potentiel, 2026-08-14).
Bruit mesure = RMS de (F - mediane glissante 21 pts) / pic, sur la phase de charge."""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))

# --- Abaqus : indenter_voronoi_rod (2D, tige V=11, unites mm-t-s-MPa-N) ---
d = np.load(r"C:\Users\fuzquianoalricabi\simulations\FDEM\percussion_2d_rod"
            r"\indenter_voronoi_rod_data.npz", allow_pickle=True)
cfn = d["h_Step_1_ElementSetPIBATCH_CFN2 on surface ASSEMBLY_GRANITE-1_SURF_TOP_IMPACT"]
u2 = d["h_Step_1_NodePART_INSERT_1_1_U2"]
tA, FA = cfn[:, 0], cfn[:, 1]          # N (mm-t-s -> N)
uA = np.interp(tA, u2[:, 0], u2[:, 1])  # mm
dA = -(uA - uA[0])                      # enfoncement positif [mm]
FA = np.abs(FA) * 1e-3                  # kN

# --- rockim : banc moyen 3D ---
rows = [r for r in csv.DictReader(open(os.path.join(here, "history_mid.csv")))
        if all(r.get(k) not in (None, "") for k in r)]
tR = np.array([float(r["t"]) for r in rows])
vz = np.array([float(r["grpVz"]) for r in rows])
FR = np.array([float(r["grpFz"]) for r in rows]) * 1e-3
z0 = float(rows[0]["grpZ"])
z = z0 + np.concatenate(([0.0], np.cumsum(0.5 * (vz[1:] + vz[:-1]) * np.diff(tR))))
dR = np.maximum(0.0, (0.131 - z)) * 1e3

def bruit(F):
    n = 21
    med = np.array([np.median(F[max(0, i - n // 2):i + n // 2 + 1])
                    for i in range(len(F))])
    hp = F - med
    return np.sqrt(np.mean(hp**2)) / F.max() * 100.0, med

iA = FA > 0.02 * FA.max()
bA, medA = bruit(FA[iA])
bR, _ = bruit(FR[FR > 0.02 * FR.max()])

fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.6), sharey=False)
a = ax[0]
a.plot(dA[iA], FA[iA], color="C7", lw=0.5, alpha=0.8, label="brut")
a.plot(dA[iA], medA, color="C0", lw=1.8, label="médiane glissante")
a.set_xlabel("enfoncement (mm)"); a.set_ylabel("force de contact (kN)")
a.set_title("Abaqus FDEM 2D — indenter_voronoi_rod (06/2026)\n"
            "pénalité + cohésifs, tige V = 11 m/s — bruit RMS %.1f %% du pic" % bA)
a.legend(fontsize=8); a.grid(alpha=0.3)

a = ax[1]
a.plot(dR, FR, color="C2", lw=1.2, label="brut (aucun filtrage)")
a.set_xlabel("pénétration δ (mm)")
a.set_title("rockim FDEM 3D — banc moyen (08/2026)\n"
            "contact potentiel Munjiza + adaptatif, insert 8 m/s — bruit RMS %.1f %% du pic" % bR)
a.legend(fontsize=8); a.grid(alpha=0.3)
fig.suptitle("Propreté des F–δ : pénalité Abaqus vs potentiel rockim (cas différents — comparer le BRUIT, pas les amplitudes)")
fig.tight_layout()
out = os.path.join(here, "compare_fdelta_abaqus_rockim.png")
fig.savefig(out, dpi=150)
print(out, "| bruit Abaqus = %.2f %% | bruit rockim = %.2f %% | ratio = %.0f x" % (bA, bR, bA / bR))
