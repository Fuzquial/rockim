#!/usr/bin/env python3
"""fig_loads.py — figures des charges et CL par groupes (fdem3d, 2026-10-03).

  python3 tools/fig_loads.py MAILLAGE.msh RUN_ELASTIQUE RUN_TENSION RUN_LOADS SORTIE/

  1. fig_loads_maillage   : le maillage et ses groupes physiques nommes
  2. fig_loads_elastique  : banc A (traction 5 MPa) — u_z(z) des noeuds contre
                            la solution exacte, histoires U_top_z et RF_bottom_z
  3. fig_loads_rupture    : banc de rupture — contrainte d appui contre le temps,
                            scenario tension et scenario loads superposes,
                            joints inseres et rompus a la fin
PDF vectoriel (Computer Modern via mathtext) + PNG.
"""
import csv, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "font.serif": ["DejaVu Serif"], "font.size": 9,
                     "axes.titlesize": 9, "legend.fontsize": 8})

E_MOD, SIG, L, A = 50e9, 5e6, 0.04, 0.02 * 0.02


def read_msh(path):
    names, nodes, tris, pts = {}, {}, [], []
    with open(path) as f:
        lines = f.read().split("\n")
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if s == "$PhysicalNames":
            n = int(lines[i + 1])
            for k in range(n):
                d, pid, nm = lines[i + 2 + k].split(maxsplit=2)
                names[int(pid)] = (int(d), nm.strip('"'))
            i += n + 2
        elif s == "$Nodes":
            n = int(lines[i + 1])
            for k in range(n):
                t = lines[i + 2 + k].split()
                nodes[int(t[0])] = np.array(list(map(float, t[1:4])))
            i += n + 2
        elif s == "$Elements":
            n = int(lines[i + 1])
            for k in range(n):
                t = list(map(int, lines[i + 2 + k].split()))
                typ, nt = t[1], t[2]
                phys, conn = t[3], t[3 + nt:]
                if typ == 2:
                    tris.append((phys, conn))
                elif typ == 15:
                    pts.append((phys, conn[0]))
            i += n + 2
        else:
            i += 1
    return names, nodes, tris, pts


def read_vtu(path):
    txt = open(path).read()
    out = {}
    for m in re.finditer(r'<DataArray([^>]*)>(.*?)</DataArray>', txt, re.S):
        attrs, body = m.group(1), m.group(2)
        nm = re.search(r'Name="([^"]*)"', attrs)
        nc = re.search(r'NumberOfComponents="(\d+)"', attrs)
        key = nm.group(1) if nm else "Points"
        arr = np.array(body.split(), dtype=float)
        if nc and int(nc.group(1)) > 1:
            arr = arr.reshape(-1, int(nc.group(1)))
        out.setdefault(key, arr)
    return out


def hist(path):
    rows = list(csv.DictReader(open(path)))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def save(fig, out, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, name + "." + ext), dpi=200,
                    bbox_inches="tight")
    plt.close(fig)


def fig_mesh(msh, out):
    names, nodes, tris, pts = read_msh(msh)
    col = {"top": "#d62728", "bottom": "#1f77b4", "xmin": "#cfd8dc",
           "xmax": "#b0bec5", "ymin": "#eceff1", "ymax": "#90a4ae"}
    fig = plt.figure(figsize=(6.6, 3.6))
    for k, (elev, azim) in enumerate([(22, -55), (-22, 125)]):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        for nm in ("xmin", "xmax", "ymin", "ymax", "bottom", "top"):
            P = [np.array([nodes[c] for c in conn]) * 1e3
                 for phys, conn in tris if names.get(phys, (0, ""))[1] == nm]
            ax.add_collection3d(Poly3DCollection(
                P, facecolor=col[nm], edgecolor="k", linewidths=0.15,
                alpha=0.95))
        for phys, nid in pts:
            nm = names[phys][1]
            x = nodes[nid] * 1e3
            ax.scatter(*x, color="k", s=6)
            if k == 0 and nm in ("c000", "c100", "c111"):
                off = {"c000": [-1, -6, -2], "c100": [1, -3, -3],
                       "c111": [1, 1, 2]}[nm]
                ax.text(*(x + np.array(off)), nm, fontsize=7)
        ax.set_xlim(0, 20); ax.set_ylim(0, 20); ax.set_zlim(0, 40)
        ax.set_box_aspect((1, 1, 2))
        ax.view_init(elev, azim)
        if k == 0:
            ax.set_xlabel("$x$ [mm]", labelpad=-6)
            ax.set_ylabel("$y$ [mm]", labelpad=-6)
            ax.set_zlabel("$z$ [mm]", labelpad=-6)
            ax.tick_params(pad=-2, labelsize=7)
        else:
            ax.set_axis_off()
        fig.text(0.31 if k == 0 else 0.71, 0.86,
                 "vue de dessus" if k == 0 else "vue de dessous",
                 ha="center")
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=col["top"], label="top : traction / vitesse"),
                        Patch(color=col["bottom"], label="bottom : appui"),
                        Patch(color=col["xmax"], label="faces laterales (libres)")],
               loc="lower center", ncol=3, frameon=False)
    fig.suptitle("meshes/loads_bar_h4.msh : 20 x 20 x 40 mm, 1 662 tetraedres, "
                 "groupes physiques Gmsh nommes (points noirs : coins c000...c111)",
                 y=0.98)
    save(fig, out, "fig_loads_maillage")


def fig_elastic(runA, out):
    v0 = read_vtu(os.path.join(runA, "fdem3d_0000.vtu"))
    last = sorted(f for f in os.listdir(runA)
                  if re.match(r"fdem3d_\d+\.vtu", f))[-1]
    v1 = read_vtu(os.path.join(runA, last))
    X0, X1 = v0["Points"], v1["Points"]
    uz = X1[:, 2] - X0[:, 2]
    h = hist(os.path.join(runA, "history.csv"))
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
    z = np.linspace(0, L, 50)
    ax[0].plot(X0[:, 2] * 1e3, uz * 1e6, ".", ms=1.5, color="0.45",
               label="noeuds (6 648 copies)")
    ax[0].plot(z * 1e3, SIG / E_MOD * z * 1e6, "r-", lw=1,
               label=r"exact $\sigma z/E$")
    ax[0].set_xlabel("$z$ [mm]"); ax[0].set_ylabel(r"$u_z$ [$\mu$m]")
    ax[0].legend(frameon=False, loc="upper left"); ax[0].set_title("fin du run, t = 0,3 ms")
    t = h["t"] * 1e6
    ax[1].plot(t, h["U_top_z"] * 1e6, "k-", lw=1, label=r"$U_{top,z}$ calcule")
    ax[1].axhline(SIG * L / E_MOD * 1e6, color="r", lw=0.8, ls="--",
                  label=r"exact $\sigma L/E$ = 4,000")
    ax[1].set_xlabel(r"$t$ [$\mu$s]"); ax[1].set_ylabel(r"$U_z$ [$\mu$m]")
    ax[1].legend(frameon=False, loc="lower right")
    ax[1].set_title(r"fin : 4,00122 $\mu$m (+0,03 %)")
    ax[2].plot(t, h["F_top_z"], "r-", lw=1, label=r"charge $F_{top,z}$")
    ax[2].plot(t, -h["RF_bottom_z"], "k-", lw=0.8,
               label=r"$-$reaction $RF_{bottom,z}$")
    ax[2].set_xlabel(r"$t$ [$\mu$s]"); ax[2].set_ylabel("force [N]")
    ax[2].legend(frameon=False, loc="lower right")
    ax[2].set_title("fin : 2 000 N contre 2 000,49 N")
    fig.suptitle("Banc A : barre elastique, traction.top = 5 MPa en ramp 0,1 ms, "
                 "fix.bottom = z (insertion adaptative, aucun joint insere)",
                 y=1.03)
    fig.tight_layout()
    save(fig, out, "fig_loads_elastique")


def fig_rupture(runT, runL, out):
    ht, hl = hist(os.path.join(runT, "history.csv")), hist(os.path.join(runL, "history.csv"))
    last = sorted(f for f in os.listdir(runL)
                  if re.match(r"fdem3d_joints_\d+\.vtu", f))[-1]
    J = read_vtu(os.path.join(runL, last))
    P, conn = J["Points"], J["connectivity"].astype(int).reshape(-1, 3)
    # regle de la base : rompu = tBreak >= 0 (jamais damage >= 0,999, qui
    # compte aussi les joints morts en compression : 380 faux positifs mesures)
    tb, bonded = J["tBreak"], J["bonded"]
    fig = plt.figure(figsize=(7.2, 3.0))
    ax = fig.add_subplot(1, 3, (1, 2))
    ax.plot(ht["t"] * 1e6, np.abs(ht["gripFz"]) / A / 1e6, color="0.6", lw=2.2,
            label="scenario = tension (mors cables)")
    ax.plot(hl["t"] * 1e6, hl["RF_top_z"] / A / 1e6, "k-", lw=0.8,
            label="scenario = loads (fix.bottom, velocity.top)")
    ax.set_xlabel(r"$t$ [$\mu$s]"); ax.set_ylabel(r"$F_{appui}/A$ [MPa]")
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("pic 12,27 MPa des deux cotes ; ecart max 0,1 % sur 2 061 lignes")
    a3 = fig.add_subplot(1, 3, 3, projection="3d")
    ins = bonded < 0.5
    br = tb >= 0.0
    tri = P[conn] * 1e3
    D = J["damage"]
    cmap = plt.get_cmap("YlOrRd")
    lo = ins & (D < 0.5)
    hi = ins & (D >= 0.5)
    a3.add_collection3d(Poly3DCollection(tri[lo], facecolor="#ffe0b2",
                                         edgecolor="none", alpha=0.25))
    a3.add_collection3d(Poly3DCollection(tri[hi], facecolor=cmap(0.4 + 0.6 * D[hi]),
                                         edgecolor="k", linewidths=0.15))
    for x in (0, 20):
        for y in (0, 20):
            a3.plot([x, x], [y, y], [0, 40], color="0.5", lw=0.4)
    for z in (0, 40):
        a3.plot([0, 20, 20, 0, 0], [0, 0, 20, 20, 0], [z] * 5, color="0.5", lw=0.4)
    a3.set_xlim(0, 20); a3.set_ylim(0, 20); a3.set_zlim(0, 40)
    a3.set_box_aspect((1, 1, 2)); a3.view_init(18, -55)
    a3.tick_params(pad=-2, labelsize=7)
    a3.set_title("%d joints inseres, %d a D >= 0,5 (fonces)\n"
                 "%d a D >= 0,999, %d rompus (tBreak >= 0)"
                 % (int(ins.sum()), int(hi.sum()), int((ins & (D >= 0.999)).sum()),
                    int(br.sum())), fontsize=8)
    fig.tight_layout()
    save(fig, out, "fig_loads_rupture")


if __name__ == "__main__":
    msh, ra, rt, rl, out = sys.argv[1:6]
    os.makedirs(out, exist_ok=True)
    fig_mesh(msh, out)
    fig_elastic(ra, out)
    fig_rupture(rt, rl, out)
