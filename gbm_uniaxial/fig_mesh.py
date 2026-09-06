# -*- coding: utf-8 -*-
"""Vues de la microstructure GBM en elements finis : coupe mediane et face laterale,
colorees par phase minerale puis par module d'Young. Lecture de la premiere image
(geometrie NON deformee).      python fig_mesh.py [out_U_gbm]
"""
from __future__ import print_function
import sys, os, re, base64, zlib
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

plt.rcParams.update({'font.family': 'serif', 'font.serif': ['CMU Serif', 'DejaVu Serif'],
                     'mathtext.fontset': 'cm', 'font.size': 9, 'pdf.fonttype': 42})
PHASES = ['feldspath', 'quartz', 'biotite']
COL = ['#d9c8a9', '#8fb8de', '#4a3b2a']          # feldspath clair, quartz bleu, biotite sombre


def arrays(path):
    s = open(path, 'r', errors='ignore').read()
    out = {}
    for m in re.finditer(r'<DataArray([^>]*)>(.*?)</DataArray>', s, re.S):
        att, body = m.group(1), m.group(2)
        nm = re.search(r'Name="([^"]+)"', att)
        if not nm:
            continue
        out[nm.group(1)] = np.fromstring(body.strip(), sep=' ')
    p = re.search(r'<Points>\s*<DataArray[^>]*>(.*?)</DataArray>', s, re.S)
    pts = np.fromstring(p.group(1).strip(), sep=' ').reshape(-1, 3)
    return pts, out


def main(outdir='out_U_gbm'):
    pts, A = arrays(os.path.join(outdir, 'fem3d_0000.vtu'))
    conn = A['connectivity'].astype(np.int64).reshape(-1, 4)
    cen = pts[conn].mean(axis=1) * 1e3                     # mm
    ph = A['phase'].astype(int)
    E = A['matE'] / 1e9
    W, H = pts[:, 0].max() * 1e3, pts[:, 2].max() * 1e3
    vol = np.abs(np.einsum('ij,ij->i', np.cross(pts[conn[:, 1]] - pts[conn[:, 0]],
                                                pts[conn[:, 2]] - pts[conn[:, 0]]),
                           pts[conn[:, 3]] - pts[conn[:, 0]])) / 6.0
    lc = (vol ** (1.0 / 3.0)) * 1e3

    cmap = ListedColormap(COL); norm = BoundaryNorm([-.5, .5, 1.5, 2.5], 3)
    fig, ax = plt.subplots(1, 3, figsize=(8.6, 3.6), constrained_layout=True)

    sl = np.abs(cen[:, 1] - W / 2) < 0.5                    # coupe mediane, slab 1 mm
    s = (lc[sl] * 7.0) ** 2
    ax[0].scatter(cen[sl, 0], cen[sl, 2], c=ph[sl], cmap=cmap, norm=norm, s=s,
                  marker='s', linewidths=0, rasterized=True)
    ax[0].set_title(u'coupe médiane — phase', fontsize=9)

    ax[1].scatter(cen[sl, 0], cen[sl, 2], c=E[sl], cmap='viridis', s=s, marker='s',
                  linewidths=0, vmin=25, vmax=85, rasterized=True)
    ax[1].set_title(u"coupe médiane — module $E$ [GPa]", fontsize=9)

    fa = np.abs(cen[:, 1] - W) < 0.8                        # face laterale
    s2 = (lc[fa] * 7.0) ** 2
    m = ax[2].scatter(cen[fa, 0], cen[fa, 2], c=E[fa], cmap='viridis', s=s2, marker='s',
                      linewidths=0, vmin=25, vmax=85, rasterized=True)
    ax[2].set_title(u'face latérale — module $E$', fontsize=9)
    fig.colorbar(m, ax=ax[2], fraction=0.05, pad=0.02, label='GPa')

    for a in ax:
        a.set_aspect('equal'); a.set_xlim(0, W); a.set_ylim(0, H)
        a.set_xlabel('x [mm]'); a.set_xticks([0, 8, 16])
    ax[0].set_ylabel('z [mm]')
    for a in ax[1:]:
        a.set_yticklabels([])
    h = [plt.Line2D([], [], marker='s', ls='', mfc=c, mec='none', ms=7) for c in COL]
    n = np.bincount(ph, minlength=3)
    ax[0].legend(h, ['%s  %.0f %%' % (p, 100.0 * np.sum(vol[ph == i]) / vol.sum())
                     for i, p in enumerate(PHASES)], fontsize=7, loc='upper center',
                 bbox_to_anchor=(0.5, -0.13), ncol=1, frameon=False)
    fig.suptitle(u'Microstructure GBM en éléments finis — %d grains, %d tétraèdres, '
                 u'éprouvette 16 x 16 x 32 mm' % (len(np.unique(A['grain'])), len(conn)),
                 fontsize=9)
    os.makedirs('figures', exist_ok=True)
    for e in ('pdf', 'png'):
        fig.savefig('figures/mesh_gbm.' + e, dpi=200)
    print('-> figures/mesh_gbm.pdf et .png ;  %d grains, %d tets, lc median %.2f mm'
          % (len(np.unique(A['grain'])), len(conn), np.median(lc)))


if __name__ == '__main__':
    main(*(sys.argv[1:] or []))
