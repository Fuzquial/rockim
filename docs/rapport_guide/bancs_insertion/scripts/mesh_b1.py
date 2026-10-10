# -*- coding: utf-8 -*-
"""mesh_b1.py — mini-maillage du banc 1 : un joint traversant une bande, incline de theta.

Bande W x H (mm -> m). Le joint cible relie (0, a) a (W, a + W tan theta). Chaque
moitie (quadrilatere) est decoupee en 4 triangles autour de son barycentre :
8 triangles, 1 joint cible + 6 joints internes (aretes coin-centre) presque verticaux.
Sous contrainte uniaxiale sigma_yy (nu = 0, mors laterallement libres) le champ est
UNIFORME et exact pour des CST : le joint cible voit sigma_n = sigma cos^2 theta,
tau = sigma sin theta cos theta ; on imprime le rapport de charge des aretes parasites.
usage : python mesh_b1.py theta_deg out.msh [W_mm a_mm b_mm]
"""
import math, sys

def write_msh(path, P, T):
    with open(path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n%d\n" % len(P))
        for i, (x, y) in enumerate(P):
            f.write("%d %.12g %.12g 0\n" % (i + 1, x, y))
        f.write("$EndNodes\n$Elements\n%d\n" % len(T))
        for k, t in enumerate(T):
            f.write("%d 2 2 1 1 %d %d %d\n" % (k + 1, t[0] + 1, t[1] + 1, t[2] + 1))
        f.write("$EndElements\n")

def build(theta, W=0.010, a=0.040, b=0.040):
    yR = a + W * math.tan(math.radians(theta))
    H = yR + b
    P = [(0, 0), (W, 0), (W, yR), (0, a), (W, H), (0, H)]
    def quad(i, j, k, l):
        cx = sum(P[n][0] for n in (i, j, k, l)) / 4; cy = sum(P[n][1] for n in (i, j, k, l)) / 4
        P.append((cx, cy)); c = len(P) - 1
        return [(i, j, c), (j, k, c), (k, l, c), (l, i, c)]
    T = quad(0, 1, 2, 3) + quad(3, 2, 4, 5)
    return P, T, H

def ratios(P, T, theta, ft, c, sig):
    """rapports sigma_n/ft et |tau|/c de toutes les aretes internes sous sigma_yy = sig"""
    from collections import Counter
    E = Counter()
    for t in T:
        for u in range(3):
            e = tuple(sorted((t[u], t[(u + 1) % 3]))); E[e] += 1
    out = []
    for (i, j), n in E.items():
        if n != 2: continue
        dx, dy = P[j][0] - P[i][0], P[j][1] - P[i][1]; L = math.hypot(dx, dy)
        nx, ny = -dy / L, dx / L
        sn = sig * ny * ny; tau = abs(sig * nx * ny)
        out.append(((i, j), sn / ft, tau / c))
    return out

if __name__ == "__main__":
    th = float(sys.argv[1]); out = sys.argv[2]
    P, T, H = build(th)
    write_msh(out, P, T)
    print("ecrit", out, "H =", H)
    for e, r1, r2 in ratios(P, T, th, 1.0, 1.0, 1.0):
        print("arete", e, "sigma_n/sigma = %.3f  |tau|/sigma = %.3f" % (r1, r2))
