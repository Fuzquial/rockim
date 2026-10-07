import sys, math, numpy as np
sys.path.insert(0, "calib_quick")
from plot_tess import read_vtu
def intra_metrics(d):
    pts, tri, grain, phase = read_vtu(d + "/fdem_0000.vtu")
    # aretes INTRA-grain : les deux triangles incidents dans le meme grain ; cle geometrique
    key = {}
    for t, g in zip(tri, grain):
        for k in range(3):
            P, Q = pts[t[k]], pts[t[(k + 1) % 3]]
            kk = tuple(sorted((tuple(np.round(P, 7)), tuple(np.round(Q, 7)))))
            key.setdefault(kk, []).append(g)
    th = []
    for kk, gs in key.items():
        if len(gs) == 2 and gs[0] == gs[1]:
            (x1, y1), (x2, y2) = kk
            th.append(math.atan2(y2 - y1, x2 - x1) % math.pi)
    th = np.array(th)
    hist, _ = np.histogram(th, bins=36, range=(0, math.pi))
    r6 = abs(np.exp(6j * th).mean())
    u = pts[tri[:, 1]] - pts[tri[:, 0]]; v = pts[tri[:, 2]] - pts[tri[:, 0]]
    a, b, c = (np.linalg.norm(pts[tri[:,1]]-pts[tri[:,2]],axis=1), np.linalg.norm(pts[tri[:,2]]-pts[tri[:,0]],axis=1), np.linalg.norm(pts[tri[:,0]]-pts[tri[:,1]],axis=1))
    ang = lambda x, y, z: np.degrees(np.arccos(np.clip((y*y+z*z-x*x)/(2*y*z), -1, 1)))
    amin = np.minimum(np.minimum(ang(a,b,c), ang(b,c,a)), ang(c,a,b))
    return len(th), hist.max()/max(1,hist.min()), r6, float(np.median(amin)), len(tri)
for name, d in (("reseau (defaut)", "calib_quick/_bit_n"), ("grainMeshRandom", "calib_quick/_gr_a")):
    n, pc, r6, am, nt = intra_metrics(d)
    print(f"{name:18s} {nt:5d} tri  aretes intra {n:5d}  pic/creux {pc:5.2f}  R6 {r6:6.3f}  angle min median {am:5.1f} deg")
