#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# plot_impact3d.py — planche d'un run d'indentation 3D, par COUPE MEDIANE.
#
#   python tunnel_edz/plot_impact3d.py out_indent3d_ye [--slab 0.001]
#                                      [--zoom 0.012] [--toolr 0.004]
#
# On garde les tetraedres dont le centre est a moins de `slab` du plan
# y = D/2, et on les trace dans le plan (x, z) — un nuage de centres colore
# par le champ, la taille du marqueur suivant la taille de l'element. C'est
# volontairement rustique : une vraie coupe conforme demanderait de decouper
# les tetraedres, ce qui n'apporterait rien pour un controle.
#
# Champs disponibles en 3D : sigma1 (contrainte principale majeure), tauMax,
# vonMises, fragment. sigma1/tauMax viennent de la sortie ajoutee au solveur
# 3D le 2026-08-17.
# ---------------------------------------------------------------------------
import argparse
import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from plot_tunnel_fields import complete  # noqa: E402


def read3d(path, arrays, ncell=4):
    with open(path) as f:
        txt = f.read()
    P = np.fromstring(re.search(r'<Points>.*?<DataArray[^>]*>(.*?)</DataArray>',
                                txt, re.S).group(1), sep=" ").reshape(-1, 3)
    C = np.fromstring(re.search(r'Name="connectivity"[^>]*>(.*?)</DataArray>',
                                txt, re.S).group(1),
                      sep=" ").astype(int).reshape(-1, ncell)
    out = {}
    for nm in arrays:
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % nm, txt, re.S)
        out[nm] = np.fromstring(m.group(1), sep=" ") if m else None
    return P, C, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default="out_indent3d_ye")
    ap.add_argument("--slab", type=float, default=0.001)
    ap.add_argument("--zoom", type=float, default=0.012)
    ap.add_argument("--toolr", type=float, default=0.004)
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--titre", default=None)
    a = ap.parse_args()
    run = a.run if os.path.isabs(a.run) else os.path.join(HERE, "..", a.run)
    out = os.path.join(HERE, os.path.basename(run) + "_planche.png")

    el = [f for f in sorted(glob.glob(os.path.join(run, "fdem3d_[0-9]*.vtu")))
          if complete(f)]
    jn = [f for f in sorted(glob.glob(os.path.join(run, "fdem3d_joints_[0-9]*.vtu")))
          if complete(f)]
    k = min(len(el), len(jn)) - 1 if a.frame < 0 else a.frame
    print(f"trame {k} (sur {len(el)} completes)")

    P0, C, _ = read3d(el[0], [])
    P, _, S = read3d(el[k], ["sigma1", "tauMax", "vonMises"])
    W, D, H = P0[:, 0].max(), P0[:, 1].max(), P0[:, 2].max()
    cx, cy, cz = 0.5 * W, 0.5 * D, H

    ctr0 = P0[C].mean(axis=1)
    slab = np.abs(ctr0[:, 1] - cy) < a.slab
    umag = np.linalg.norm(P - P0, axis=1)
    uel = umag[C].mean(axis=1)
    vol = np.abs(np.einsum("ij,ij->i", P0[C[:, 1]] - P0[C[:, 0]],
                           np.cross(P0[C[:, 2]] - P0[C[:, 0]],
                                    P0[C[:, 3]] - P0[C[:, 0]]))) / 6.0
    ms = 90.0 * (vol ** (1 / 3)) / (vol ** (1 / 3)).max()

    PJ, CJ, SJ = read3d(jn[k], ["breakMode"], ncell=3)
    bm = SJ["breakMode"]
    jc = PJ[CJ].mean(axis=1)
    brk = (bm > 0) & (np.abs(jc[:, 1] - cy) < 2.5 * a.slab)

    near = slab & (np.hypot(ctr0[:, 0] - cx, ctr0[:, 2] - cz) < 0.6 * a.zoom)
    def lim(f, q=98.0):
        return np.percentile(np.abs(f[near]), q) if near.any() else 1.0

    fig = plt.figure(figsize=(16.5, 9.2))
    panels = [(None, None, None, "(a) fissures dans la coupe"),
              (S["sigma1"] / 1e6, "RdBu_r", "s1", r"(b) $\sigma_1$ [MPa]"),
              (S["tauMax"] / 1e6, "inferno", "t", r"(c) $\tau_{max}$ [MPa]"),
              (S["vonMises"] / 1e6, "viridis", "t", "(d) von Mises [MPa]"),
              (uel * 1e3, "turbo", "t", f"(e) |u| [mm] — max {umag.max()*1e3:.3f}")]
    for i, (fld, cmap, kind, ttl) in enumerate(panels):
        ax = fig.add_subplot(2, 3, i + 1)
        if fld is None:
            for m, c, lab in ((1, "#1B8A3A", "traction"), (2, "#C8342B", "cisaillement")):
                s = brk & (bm == m)
                ax.plot(jc[s, 0], jc[s, 2], ".", color=c, ms=4, label=f"{lab} ({int(s.sum())})")
            ax.legend(fontsize=8, loc="lower left")
        else:
            v = lim(fld)
            kw = dict(vmin=-v, vmax=v) if kind == "s1" else dict(vmin=0, vmax=v)
            h = ax.scatter(ctr0[slab, 0], ctr0[slab, 2], c=fld[slab], s=ms[slab],
                           cmap=cmap, marker="o", linewidths=0, **kw)
            fig.colorbar(h, ax=ax, shrink=0.82)
        if a.toolr > 0:
            th = np.linspace(np.pi, 2 * np.pi, 120)
            ax.plot(cx + a.toolr * np.cos(th), cz + a.toolr + a.toolr * np.sin(th),
                    color="#0B4F9E", lw=2.2)
        ax.set_aspect("equal")
        ax.set_xlim(cx - a.zoom, cx + a.zoom)
        ax.set_ylim(cz - 1.6 * a.zoom, cz + 0.15 * a.zoom)
        ax.set_title(ttl, fontsize=10.5)
        ax.set_xlabel("x [m]")
        if i % 3 == 0:
            ax.set_ylabel("z [m]")

    ax = fig.add_subplot(2, 3, 6)
    d = np.genfromtxt(os.path.join(run, "history.csv"), delimiter=",",
                      names=True, invalid_raise=False)
    t, vz, f = d["t"], d["toolVz"], np.abs(d["toolFz"]) / 1e3
    dlt = -np.concatenate([[0.], np.cumsum(0.5*(vz[1:]+vz[:-1])*np.diff(t))]) * 1e3
    ax.plot(dlt, f, color="#0B4F9E", lw=1.8)
    ax2 = ax.twinx()
    ax2.plot(dlt, d["nBroken"], color="#C8342B", lw=1.5)
    ax2.set_ylabel("joints rompus", color="#C8342B")
    ax.set_xlabel("enfoncement du bouton [mm]")
    ax.set_ylabel("force outil |Fz| [kN]", color="#0B4F9E")
    ax.set_title("(f) force-penetration", fontsize=10.5)
    ax.grid(alpha=0.3)
    print(f"  pic {f.max():.2f} kN a {dlt[np.argmax(f)]:.3f} mm, "
          f"enfoncement max {dlt.max():.3f} mm, {int(d['nBroken'][-1])} fissures")

    fig.suptitle(a.titre or f"{os.path.basename(run)} — coupe mediane "
                 f"(|y - D/2| < {a.slab*1e3:g} mm)", fontsize=12.5)
    fig.tight_layout()
    fig.savefig(out, dpi=155)
    print("ecrit :", out)


if __name__ == "__main__":
    main()
