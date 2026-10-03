#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# V1 : figures du maillage et de la propagation (lecture directe des VTU ASCII
# de rockim, sans VTK ni ParaView).
#
#   python3 vv/V1_onde/viz_v1.py [--exe build_nofma/rockim] [--h 2] [--frames 25]
#
# Lance un run fem3d de la barre V1 avec `frames` trames, puis trace :
#   fig_maillage.pdf     : face y = 0 de la barre entiere et vue 3D des 40 premiers mm
#   fig_propagation.pdf  : face y = 0 coloree par la contrainte axiale a plusieurs
#                          instants (pression du VTU, signe change : compression < 0)
#   fig_profils.pdf      : deplacement axial u(z) des noeuds contre la solution exacte
# ---------------------------------------------------------------------------
import argparse, os, re, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v1_onde as V  # noqa: E402

ROOT = V.ROOT


def read_vtu(path):
    txt = open(path).read()

    def arr(name, dtype=float):
        m = re.search(rf'<DataArray[^>]*Name="{name}"[^>]*>(.*?)</DataArray>', txt, re.S)
        return np.array(m.group(1).split(), dtype=dtype)

    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", txt, re.S)
    pts = np.array(m.group(1).split(), dtype=float).reshape(-1, 3)
    conn = arr("connectivity", int).reshape(-1, 4)
    out = dict(pts=pts, tets=conn)
    for nm in ("pressure", "velocity"):
        try:
            out[nm] = arr(nm)
        except AttributeError:
            pass
    if "velocity" in out:
        out["velocity"] = out["velocity"].reshape(-1, 3)
    return out


def face_y0(pts0, tets, tol=1e-9):
    """triangles de bord sur y = 0, avec l'indice du tet porteur"""
    tris, owner = [], []
    for k, t in enumerate(tets):
        on = [i for i in t if abs(pts0[i, 1]) < tol]
        if len(on) == 3:
            tris.append(on)
            owner.append(k)
    return np.array(tris), np.array(owner)


def run_viz(exe, h, frames):
    od = os.path.join(V.OUT, f"viz_fem3d_h{V.tag(h)}")
    os.makedirs(od, exist_ok=True)
    V.ensure_mesh(h)
    cfg = V.write_deck("fem3d", h, od)
    txt = open(cfg).read().replace("frames = 1", f"frames = {frames}")
    open(cfg, "w").write(txt)
    with open(os.path.join(od, "rockim.log"), "w") as log:
        subprocess.run([exe, cfg], stdout=log, stderr=subprocess.STDOUT, cwd=ROOT, check=True)
    return od


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default=os.path.join(ROOT, "build_nofma", "rockim"))
    ap.add_argument("--h", type=float, default=2.0)
    ap.add_argument("--frames", type=int, default=25)
    ap.add_argument("--norun", action="store_true")
    a = ap.parse_args()
    od = os.path.join(V.OUT, f"viz_fem3d_h{V.tag(a.h)}")
    if not a.norun:
        od = run_viz(os.path.abspath(a.exe), a.h, a.frames)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIX Two Text"],
                         "mathtext.fontset": "stix", "font.size": 9})

    fr = np.genfromtxt(os.path.join(od, "frames.csv"), delimiter=",", names=True)
    files = sorted(f for f in os.listdir(od) if re.match(r"fem3d_\d+\.vtu", f))
    v0 = read_vtu(os.path.join(od, files[0]))
    X0, tets = v0["pts"], v0["tets"]
    tris, owner = face_y0(X0, tets)
    zmm, xmm = X0[:, 2] * 1e3, X0[:, 0] * 1e3

    # ---- 1. maillage -------------------------------------------------------
    fig = plt.figure(figsize=(7.0, 4.6))
    ax = fig.add_axes([0.06, 0.70, 0.90, 0.22])
    polys = [np.c_[zmm[t], xmm[t]] for t in tris]
    ax.add_collection(PolyCollection(polys, facecolor="0.92", edgecolor="0.25", lw=0.15))
    for g, z in V.STATIONS.items():
        ax.axvline(z * 1e3, color="C3", lw=0.8, ls="--")
        ax.text(z * 1e3, 9.5, g, color="C3", ha="center", fontsize=8)
    ax.text(-2, 4, "bottom\n(vitesse\nimposée)", ha="right", va="center", fontsize=7)
    ax.text(402, 4, "top\n(libre)", ha="left", va="center", fontsize=7)
    ax.set_xlim(-30, 430); ax.set_ylim(-1, 12)
    ax.set_aspect("equal"); ax.set_xlabel("z (mm)"); ax.set_ylabel("x (mm)")
    ax.set_title(f"Face y = 0 de la barre, h = {a.h:g} mm, {len(tets)} tétraèdres", fontsize=9)
    ax3 = fig.add_axes([0.05, 0.02, 0.9, 0.58], projection="3d")
    sel = [t for t in tets if X0[t, 2].max() <= 0.04]
    faces = []
    for t in sel:
        for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            faces.append(X0[t[list(f)]] * 1e3)
    ax3.add_collection3d(Poly3DCollection(faces, facecolor=(0.75, 0.82, 0.92, 0.25),
                                          edgecolor="0.2", lw=0.12))
    ax3.set_xlim(0, 8); ax3.set_ylim(0, 8); ax3.set_zlim(0, 40)
    ax3.set_box_aspect((8, 8, 40)); ax3.view_init(elev=18, azim=-58)
    ax3.set_xlabel("x"); ax3.set_ylabel("y"); ax3.set_zlabel("z (mm)")
    ax3.set_title("Vue 3D des 40 premiers millimètres (arêtes des tétraèdres)", fontsize=9)
    fig.savefig(os.path.join(HERE, "fig_maillage.pdf"), dpi=200)
    fig.savefig(os.path.join(HERE, "fig_maillage.png"), dpi=200)
    plt.close(fig)

    # ---- 2. propagation ----------------------------------------------------
    want = [10e-6, 40e-6, 70e-6, 90e-6, 120e-6, 160e-6, 200e-6]
    idx = [int(np.argmin(abs(fr["t"] - w))) for w in want]
    smax = V.RHO * V.C * V.V0 * 2 / 1e6           # 2 rho c v0 : reflexion
    fig, axs = plt.subplots(len(idx), 1, figsize=(7.0, 1.0 * len(idx) + 0.8), sharex=True)
    for ax, i in zip(axs, idx):
        d = read_vtu(os.path.join(od, files[i]))
        # le champ VTU « pressure » est tr(sigma)/3, POSITIF EN TRACTION
        # (Fem3dSolver.cpp, pm = sig.trace() / 3) ; en deformation uniaxiale
        # sigma_xx = sigma_yy = nu/(1-nu) sigma_zz, d'ou sigma_zz = 3 pm (1-nu)/(1+nu)
        sz = d["pressure"][owner] * 3 * (1 - V.NU) / (1 + V.NU) / 1e6   # sigma_zz (MPa)
        pc = PolyCollection(polys, array=sz, cmap="RdBu_r", edgecolor="face", lw=0.05)
        pc.set_clim(-smax, smax)
        ax.add_collection(pc)
        ax.set_xlim(0, 400); ax.set_ylim(0, 8); ax.set_aspect(4)
        ax.set_yticks([]); ax.set_ylabel(f"{fr['t'][i] * 1e6:.0f} µs", rotation=0,
                                         ha="right", va="center")
    axs[-1].set_xlabel("z (mm)")
    cb = fig.colorbar(pc, ax=axs, fraction=0.03, pad=0.01)
    cb.set_label(r"$\sigma_{zz}$ (MPa) : compression en bleu, traction en rouge")
    fig.suptitle("Propagation de l'impulsion : contrainte axiale sur la face y = 0", fontsize=9)
    fig.savefig(os.path.join(HERE, "fig_propagation.pdf"), dpi=200)
    fig.savefig(os.path.join(HERE, "fig_propagation.png"), dpi=200)
    plt.close(fig)

    # ---- 3. profils u(z) ---------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.3, 3.6))
    zz = np.linspace(0, V.L, 801)
    for k, i in enumerate(idx):
        d = read_vtu(os.path.join(od, files[i]))
        uz = (d["pts"][:, 2] - X0[:, 2]) * 1e6
        t = fr["t"][i]
        ue = np.array([V.u_exact(z, np.array([t]))[0] for z in zz]) * 1e6
        col = f"C{k}"
        ax.plot(zz * 1e3, ue, color=col, lw=1.8, alpha=0.4)
        ax.plot(zmm, uz, ".", color=col, ms=0.8, label=f"{t * 1e6:.0f} µs")
    ax.set_xlabel("z (mm)"); ax.set_ylabel("déplacement axial (µm)")
    ax.legend(fontsize=7, ncol=4, frameon=False, markerscale=8)
    ax.set_title("Déplacement des nœuds (points) contre solution exacte (trait large)", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fig_profils.pdf"), dpi=200)
    fig.savefig(os.path.join(HERE, "fig_profils.png"), dpi=200)
    plt.close(fig)
    print("figures : fig_maillage, fig_propagation, fig_profils (.pdf et .png) dans", HERE)


if __name__ == "__main__":
    main()
