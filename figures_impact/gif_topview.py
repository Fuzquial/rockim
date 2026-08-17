# -*- coding: utf-8 -*-
"""Vue de dessus en rendu elements (style figures percussion de la these) :
surface superieure de la roche en triangles colores von Mises, fragment
detache en rouge, joints casses en rouge sombre. GIF sur toutes les frames."""
import io, os, re, sys
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm as mcm

here = os.path.dirname(os.path.abspath(__file__))
snap = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "out_banc_mid")
out_gif = sys.argv[2] if len(sys.argv) > 2 else os.path.join(here, "banc_mid_topview.gif")

def arr(text, name):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name, text, re.S)
    return np.fromstring(m.group(1), sep=" ")

def points(text):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)

hist = np.genfromtxt(os.path.join(snap, "history.csv"), delimiter=",", names=True, invalid_raise=False)
frames_meta = np.genfromtxt(os.path.join(snap, "frames.csv"), delimiter=",", names=True, invalid_raise=False)

FACES = [(0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 3, 1), (1, 2, 3, 0)]  # face + sommet oppose

vtus = sorted(f for f in os.listdir(snap) if re.fullmatch(r"fdem3d_\d{4}\.vtu", f))
imgs = []
for k, fn in enumerate(vtus):
    txt = open(os.path.join(snap, fn)).read()
    P = points(txt)
    conn = arr(txt, "connectivity").astype(int).reshape(-1, 4)
    phase = arr(txt, "phase")
    frag = arr(txt, "fragment")
    vm = arr(txt, "vonMises") * 1e-6

    jfn = os.path.join(snap, fn.replace("fdem3d_", "fdem3d_joints_"))
    jtxt = open(jfn).read()
    jP = points(jtxt)
    jconn = arr(jtxt, "connectivity").astype(int)
    joff = arr(jtxt, "offsets").astype(int)
    jtb = arr(jtxt, "tBreak")      # -1 = intact, sinon instant de rupture
    jdam = arr(jtxt, "damage")     # endommagement continu du joint

    tk = float(np.atleast_1d(frames_meta["t"])[k])
    Fz = float(np.interp(tk, hist["t"], hist["grpFz"]))
    vz = float(np.interp(tk, hist["t"], hist["grpVz"]))

    rock = phase < 0.5
    ids, counts = np.unique(frag[rock], return_counts=True)
    main = ids[np.argmax(counts)]

    # faces exterieures de la roche (tuple trie vu une seule fois)
    cnt = Counter()
    for e in np.where(rock)[0]:
        for fa in FACES:
            cnt[tuple(sorted(conn[e][list(fa[:3])]))] += 1
    polys, cols = [], []
    cmap = mcm.get_cmap("viridis")
    for e in np.where(rock)[0]:
        for fa in FACES:
            tri = conn[e][list(fa[:3])]
            key = tuple(sorted(tri))
            if cnt[key] != 1:
                continue
            A, B, C = P[tri[0]], P[tri[1]], P[tri[2]]
            n = np.cross(B - A, C - A)
            if np.dot(n, P[conn[e][fa[3]]] - A) > 0:   # orienter vers l'exterieur
                n = -n
            nz = n[2] / (np.linalg.norm(n) + 1e-30)
            detached = frag[e] != main
            if nz < 0.05 and not detached:
                continue                                # on ne garde que le dessus
            polys.append([A, B, C])
            cols.append((0.85, 0.1, 0.1) if detached
                        else cmap(min(vm[e] / 40.0, 1.0)))
    # joints rompus (tBreak >= 0) et fortement endommages (damage > 0.5)
    jpolys, jdpolys = [], []
    start = 0
    for ci, end in enumerate(joff):
        idx = jconn[start:end]
        start = end
        if len(idx) < 3:
            continue
        if jtb[ci] >= 0.0:
            jpolys.append([jP[i] for i in idx[:3]])
        elif jdam[ci] > 0.5:
            jdpolys.append([jP[i] for i in idx[:3]])

    fig = plt.figure(figsize=(7.2, 6.6))
    ax = fig.add_subplot(111, projection="3d")
    pc = Poly3DCollection(polys, facecolors=cols, edgecolors="k", linewidths=0.05)
    ax.add_collection3d(pc)
    if jdpolys:
        ax.add_collection3d(Poly3DCollection(jdpolys, facecolors=(1.0, 0.6, 0.0),
                                             edgecolors="none"))
    if jpolys:
        ax.add_collection3d(Poly3DCollection(jpolys, facecolors=(0.55, 0.0, 0.0),
                                             edgecolors="none"))
    ax.set_xlim(0, 0.12); ax.set_ylim(0, 0.12); ax.set_zlim(0.10, 0.135)
    ax.set_box_aspect((1, 1, 0.3))
    ax.view_init(elev=90, azim=-90)
    ax.set_axis_off()
    ax.set_title("Vue de dessus — t = %.0f µs   Fz = %.1f kN   vz = %+.2f m/s\n"
                 "(insert masqué ; von Mises 0-40 MPa ; fragment détaché rouge ; "
                 "joints rompus rouge sombre, endommagés > 0,5 orange)"
                 % (tk * 1e6, Fz * 1e-3, vz), fontsize=9)
    m = mcm.ScalarMappable(cmap=cmap); m.set_clim(0, 40)
    cb = fig.colorbar(m, ax=ax, shrink=0.6); cb.set_label("von Mises (MPa)")
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=100); plt.close(fig)
    buf.seek(0); imgs.append(Image := __import__("PIL.Image", fromlist=["Image"]).open(buf).convert("P"))

imgs[0].save(out_gif, save_all=True, append_images=imgs[1:], duration=800, loop=0)
print("ecrit", out_gif, "(%d frames)" % len(imgs))
