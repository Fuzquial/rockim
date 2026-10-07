# -*- coding: utf-8 -*-
"""Le MAILLAGE lui-meme : les tetraedres et leurs aretes, sur la peau de l'eprouvette,
colores par phase minerale. Panneau de droite = agrandissement pour voir les elements.
    python fig_maillage.py [out_U_gbm]
"""
from __future__ import print_function
import sys, os, re
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.collections import PolyCollection

plt.rcParams.update({'font.family': 'serif', 'font.serif': ['CMU Serif', 'DejaVu Serif'],
                     'mathtext.fontset': 'cm', 'font.size': 9, 'pdf.fonttype': 42})
PHASES = ['feldspath', 'quartz', 'biotite']
COL = ['#d9c8a9', '#8fb8de', '#4a3b2a']


def arrays(path):
    s = open(path, 'r', errors='ignore').read()
    out = {}
    for m in re.finditer(r'<DataArray([^>]*)>(.*?)</DataArray>', s, re.S):
        nm = re.search(r'Name="([^"]+)"', m.group(1))
        if nm:
            out[nm.group(1)] = np.fromstring(m.group(2).strip(), sep=' ')
    p = re.search(r'<Points>\s*<DataArray[^>]*>(.*?)</DataArray>', s, re.S)
    return np.fromstring(p.group(1).strip(), sep=' ').reshape(-1, 3), out


def peau(conn):
    """Faces qui n'appartiennent qu'a UN tetraedre = la peau. -> (faces, tet proprietaire)"""
    fi = [(1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)]
    F = np.concatenate([conn[:, list(f)] for f in fi])
    own = np.tile(np.arange(len(conn)), 4)
    key = np.sort(F, axis=1)
    o = np.lexsort((key[:, 2], key[:, 1], key[:, 0]))
    key, F, own = key[o], F[o], own[o]
    same = np.zeros(len(key), bool)
    same[:-1] |= np.all(key[1:] == key[:-1], axis=1)
    same[1:] |= np.all(key[1:] == key[:-1], axis=1)
    return F[~same], own[~same]


def main(outdir='out_U_gbm'):
    pts, A = arrays(os.path.join(outdir, 'fem3d_0000.vtu'))
    conn = A['connectivity'].astype(np.int64).reshape(-1, 4)
    ph = A['phase'].astype(int)
    F, own = peau(conn)
    P = pts * 1e3
    print('%d tets, %d faces de peau' % (len(conn), len(F)))

    fig = plt.figure(figsize=(8.8, 4.6), constrained_layout=True)
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    tri = P[F]
    pc = Poly3DCollection(tri, facecolors=[COL[p] for p in ph[own]],
                          edgecolors='0.25', linewidths=0.10)
    ax.add_collection3d(pc)
    ax.set_xlim(0, 16); ax.set_ylim(0, 16); ax.set_zlim(0, 32)
    ax.set_box_aspect((1, 1, 2)); ax.view_init(elev=18, azim=-58)
    ax.set_xlabel('x [mm]', labelpad=-6); ax.set_ylabel('y [mm]', labelpad=-6)
    ax.set_zlabel('z [mm]', labelpad=-4)
    ax.tick_params(labelsize=7, pad=-2)
    ax.set_title(u'peau de l\'éprouvette : %d faces' % len(F), fontsize=9)

    # --- agrandissement : la face y = 0, a plat, sur 10 x 10 mm ---
    ax2 = fig.add_subplot(1, 2, 2)
    n = np.cross(P[F[:, 1]] - P[F[:, 0]], P[F[:, 2]] - P[F[:, 0]])
    face_y0 = (np.abs(P[F].mean(axis=1)[:, 1]) < 1e-6)
    sel = face_y0 & (P[F].mean(axis=1)[:, 0] < 10.2) & (P[F].mean(axis=1)[:, 2] < 10.2)
    poly = [P[f][:, [0, 2]] for f in F[sel]]
    ax2.add_collection(PolyCollection(poly, facecolors=[COL[p] for p in ph[own][sel]],
                                      edgecolors='0.15', linewidths=0.35))
    ax2.set_xlim(0, 10); ax2.set_ylim(0, 10); ax2.set_aspect('equal')
    ax2.set_xlabel('x [mm]'); ax2.set_ylabel('z [mm]')
    ax2.set_title(u'agrandissement de la face $y = 0$ : les triangles sont les faces\n'
                  u'des tétraèdres ; les grains sont les plages de même couleur', fontsize=9)
    h = [plt.Line2D([], [], marker='s', ls='', mfc=c, mec='0.3', ms=8) for c in COL]
    ax2.legend(h, PHASES, fontsize=8, loc='lower right', frameon=True, framealpha=0.9)
    fig.suptitle(u'Maillage GBM en éléments finis — %d grains de Voronoï, %d tétraèdres, '
                 u'$l_c$ médian 0,58 mm' % (len(np.unique(A['grain'])), len(conn)), fontsize=9)
    os.makedirs('figures', exist_ok=True)
    for e in ('pdf', 'png'):
        fig.savefig('figures/maillage_gbm.' + e, dpi=220)
    print('-> figures/maillage_gbm.pdf et .png')


if __name__ == '__main__':
    main(*(sys.argv[1:] or []))
