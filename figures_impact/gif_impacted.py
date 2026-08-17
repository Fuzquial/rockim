# -*- coding: utf-8 -*-
"""GIF vue oblique des SEULS elements impactes (style 'damaged elements only') :
element impacte = tet adjacent a un joint endommage (D > 0.05) ou appartenant a
un fragment detache. Couleur = endommagement du joint voisin (Reds), fragments
en rouge sombre. Variante --envelope : + enveloppe du maillage (aretes vives
du bloc extraites de la frame 0 par dedoublonnage geometrique des noeuds FDEM
dupliques, silhouette translucide de l'insert), cadrage bloc entier.
Sans --envelope : cadrage zoome sur la zone d'impact."""
import io, os, re, sys
from collections import Counter, defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from matplotlib import colormaps
from PIL import Image

here = os.path.dirname(os.path.abspath(__file__))
snap = os.path.join(here, "..", "out_banc_mid")
ENV = "--envelope" in sys.argv
out_gif = os.path.join(here, "banc_mid_impacted_env.gif" if ENV else "banc_mid_impacted.gif")

def arr(text, name):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name, text, re.S)
    return np.fromstring(m.group(1), sep=" ")

def points(text):
    m = re.search(r"<Points>\s*<DataArray[^>]*>(.*?)</DataArray>", text, re.S)
    return np.fromstring(m.group(1), sep=" ").reshape(-1, 3)

hist = np.genfromtxt(os.path.join(snap, "history.csv"), delimiter=",", names=True, invalid_raise=False)
frames_meta = np.genfromtxt(os.path.join(snap, "frames.csv"), delimiter=",", names=True, invalid_raise=False)
FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
cmap = colormaps["Reds"]

# ---- topologie sur la frame 0 (noeuds FDEM dupliques -> groupes geometriques)
txt0 = open(os.path.join(snap, "fdem3d_0000.vtu")).read()
P0 = points(txt0)
conn = arr(txt0, "connectivity").astype(int).reshape(-1, 4)
phase = arr(txt0, "phase")
gid = {}
node_gid = np.empty(len(P0), dtype=int)
for i, p in enumerate(np.round(P0 * 1e7).astype(np.int64)):
    node_gid[i] = gid.setdefault(tuple(p), len(gid))
tets_g = node_gid[conn]

env_edges, env_ins_faces = [], []
if ENV:
    cnt = Counter()
    for e in range(len(conn)):
        for fa in FACES:
            cnt[tuple(sorted(tets_g[e][list(fa)]))] += 1
    edge_normals = defaultdict(list)
    edge_nodes = {}
    for e in range(len(conn)):
        for fa in FACES:
            tri_g = tuple(sorted(tets_g[e][list(fa)]))
            if cnt[tri_g] != 1:
                continue
            tri = conn[e][list(fa)]
            if phase[e] >= 0.5:
                env_ins_faces.append(tri)
                continue
            A, B, C = P0[tri[0]], P0[tri[1]], P0[tri[2]]
            n = np.cross(B - A, C - A)
            n /= (np.linalg.norm(n) + 1e-30)
            for a, b in ((0, 1), (1, 2), (0, 2)):
                key = tuple(sorted((tets_g[e][fa[a]], tets_g[e][fa[b]])))
                edge_normals[key].append(n)
                edge_nodes[key] = (tri[a], tri[b])
    for key, ns in edge_normals.items():
        if len(ns) == 2 and abs(np.dot(ns[0], ns[1])) < 0.7:
            env_edges.append(edge_nodes[key])

vtus = sorted(f for f in os.listdir(snap) if re.fullmatch(r"fdem3d_\d{4}\.vtu", f))
imgs = []
for k, fn in enumerate(vtus):
    txt = open(os.path.join(snap, fn)).read()
    P = points(txt)
    frag = arr(txt, "fragment")
    cen = P[conn].mean(axis=1)

    jtxt = open(os.path.join(snap, fn.replace("fdem3d_", "fdem3d_joints_"))).read()
    jP = points(jtxt)
    jconn = arr(jtxt, "connectivity").astype(int)
    joff = arr(jtxt, "offsets").astype(int)
    jdam = arr(jtxt, "damage")

    tk = float(np.atleast_1d(frames_meta["t"])[k])
    Fz = float(np.interp(tk, hist["t"], hist["grpFz"]))
    nB = int(np.interp(tk, hist["t"], hist["nBroken"]))

    rock = phase < 0.5
    ids, counts = np.unique(frag[rock], return_counts=True)
    main = ids[np.argmax(counts)]
    detached = rock & (frag != main)

    edam = np.zeros(len(conn))
    start = 0
    for ci, end in enumerate(joff):
        idx = jconn[start:end]; start = end
        if jdam[ci] <= 0.05 or len(idx) < 3:
            continue
        jc = jP[idx[:3]].mean(axis=0)
        m = np.einsum("ij,ij->i", cen - jc, cen - jc) < 0.004**2
        edam[m] = np.maximum(edam[m], jdam[ci])
    hit = rock & ((edam > 0.05) | detached)

    polys, cols = [], []
    for e in np.where(hit)[0]:
        c = ((0.45, 0.02, 0.02) if detached[e]
             else cmap(0.25 + 0.75 * min(edam[e], 1.0)))
        for fa in FACES:
            polys.append(P[conn[e][list(fa)]])
            cols.append(c)

    fig = plt.figure(figsize=(7.6, 6.6))
    ax = fig.add_subplot(111, projection="3d")
    if polys:
        ax.add_collection3d(Poly3DCollection(polys, facecolors=cols,
                                             edgecolors="k", linewidths=0.15))
    if ENV:
        if env_edges:
            segs = [[P[a], P[b]] for a, b in env_edges]
            ax.add_collection3d(Line3DCollection(segs, colors="0.5", linewidths=0.9))
        if env_ins_faces:
            ax.add_collection3d(Poly3DCollection([P[list(t)] for t in env_ins_faces],
                                                 facecolors=(0.55, 0.55, 0.62, 0.12),
                                                 edgecolors="none"))
        ax.set_xlim(0.0, 0.12); ax.set_ylim(0.0, 0.12); ax.set_zlim(0.02, 0.135)
        ax.set_box_aspect((1, 1, 0.95))
    else:
        ax.set_xlim(0.035, 0.09); ax.set_ylim(0.035, 0.09); ax.set_zlim(0.098, 0.125)
        ax.set_box_aspect((1, 1, 0.5))
    ax.view_init(elev=28, azim=-55)
    ax.set_axis_off()
    ax.set_title("Éléments impactés%s — t = %.0f µs   Fz = %.1f kN   %d joints cassés"
                 % (" + enveloppe" if ENV else " seuls", tk * 1e6, Fz * 1e-3, nB),
                 fontsize=10)
    m = plt.cm.ScalarMappable(cmap=cmap); m.set_clim(0, 1)
    cb = fig.colorbar(m, ax=ax, shrink=0.55); cb.set_label("endommagement du joint voisin")
    fig.tight_layout()
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=105); plt.close(fig)
    buf.seek(0); imgs.append(Image.open(buf).convert("P"))

imgs[0].save(out_gif, save_all=True, append_images=imgs[1:], duration=800, loop=0)
print("ecrit", out_gif, "(%d frames)" % len(imgs))
