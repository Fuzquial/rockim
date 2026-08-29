#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# pack_run.py — comprime une sortie de run rockim (history.csv + frames VTU)
# en UN fichier .npz transportable, pour que le depouillement puisse se faire
# ailleurs que sur la machine de calcul.
#
#   python tools/pack_run.py out_pulv_coulomb bench_impact/donnees/P1.npz
#   python tools/pack_run.py out_pulv_coulomb P1.npz --rayon 0.07
#
# CE QUI EST GARDE (et pourquoi) :
#   * history.csv en entier (float32) : F-p, vitesses, rebond, pulverisation,
#     canaux d energie — 3 figures sur 4 n ont besoin QUE de ca ;
#   * par frame, les ELEMENTS DE ROCHE proches de l impact : positions,
#     bulkD, vonMises — de quoi refaire coupes, cartes et film ;
#   * les JOINTS de la DERNIERE frame seulement : tBreak est cumulatif, donc
#     la derniere frame porte TOUTE la chronologie de fissuration ;
#   * l outil (acier + carbure) : sa position par frame, en enveloppe.
# Le reste (tets d outil, champs redondants) est jete : ~355 Mo -> ~10 Mo.
#
# Aucune dependance exotique : numpy seul, le XML VTU est lu a la main.
# ---------------------------------------------------------------------------
import sys, os, glob, csv
import xml.etree.ElementTree as ET
import numpy as np


def lire_vtu(f):
    """Retourne (points, {nom: tableau}) d un .vtu ascii ecrit par rockim."""
    root = ET.parse(f).getroot()
    piece = next(root.iter('Piece'))
    pts = np.fromstring(piece.find('Points').find('DataArray').text.replace('\n', ' '),
                        sep=' ').reshape(-1, 3)
    a = {}
    for da in root.iter('DataArray'):
        n = da.get('Name')
        if n and da.text:
            a[n] = np.fromstring(da.text.replace('\n', ' '), sep=' ')
    return pts, a


def centres(pts, A, nsom):
    """Barycentres des cellules a nsom sommets (4 = tet, 3 = triangle)."""
    conn = A['connectivity'].astype(np.int64)
    off = A['offsets'].astype(np.int64)
    st = np.concatenate(([0], off[:-1]))
    idx = np.stack([conn[st + k] for k in range(nsom)], 1)
    return pts[idx].mean(axis=1), idx


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else 'out_pulv_coulomb'
    dst = sys.argv[2] if len(sys.argv) > 2 else 'run.npz'
    rayon = 0.07
    if '--rayon' in sys.argv:
        rayon = float(sys.argv[sys.argv.index('--rayon') + 1])

    if not os.path.isdir(src):
        sys.exit(f"dossier introuvable : {src}")
    out = {}

    # ---- history.csv ------------------------------------------------------
    hpath = os.path.join(src, 'history.csv')
    if os.path.exists(hpath):
        with open(hpath) as fh:
            rd = csv.reader(fh)
            cols = next(rd)
            data = np.array([[float(v) for v in r] for r in rd if len(r) == len(cols)],
                            dtype=np.float32)
        out['hist_cols'] = np.array(cols)
        out['hist'] = data
        print(f"history.csv : {data.shape[0]} lignes x {data.shape[1]} colonnes")
    else:
        print("ATTENTION : pas de history.csv")

    # ---- frames elements --------------------------------------------------
    ef = sorted(glob.glob(os.path.join(src, 'fdem3d_0*.vtu')))
    ef = [f for f in ef if 'joints' not in os.path.basename(f)]
    if not ef:
        sys.exit("aucune frame d elements trouvee")

    # le masque est fige sur la PREMIERE frame (configuration de reference) :
    # roche uniquement (phase 0), dans un cylindre autour de l axe d impact.
    pts0, A0 = lire_vtu(ef[0])
    c0, idx0 = centres(pts0, A0, 4)
    ph = A0.get('phase', np.zeros(len(c0)))
    cx, cy = c0[:, 0].mean(), c0[:, 1].mean()          # axe = centre du bloc
    roche = ph < 0.5
    if roche.sum() == 0:                                # secours : pas de phase
        roche = np.ones(len(c0), bool)
    r2 = (c0[:, 0] - cx) ** 2 + (c0[:, 1] - cy) ** 2
    zmax = c0[roche][:, 2].max()
    garde = roche & (r2 <= rayon ** 2) & (c0[:, 2] >= zmax - 3.0 * rayon)
    outil = ~roche
    print(f"elements : {garde.sum()} gardes / {len(c0)} "
          f"(roche {roche.sum()}, outil {outil.sum()}), axe ({cx:.4f}, {cy:.4f})")

    P, BD, VM, TOOLZ, TVTX = [], [], [], [], []
    for f in ef:
        try:
            p, A = lire_vtu(f)
        except Exception as e:
            print(f"  frame illisible ignoree : {os.path.basename(f)} ({e})")
            continue
        _, idx = centres(p, A, 4)
        P.append(p[idx[garde]].astype(np.float32))          # (nel, 4, 3)
        BD.append(A.get('bulkD', np.zeros(len(idx)))[garde].astype(np.float32))
        VM.append(A.get('vonMises', np.zeros(len(idx)))[garde].astype(np.float32))
        zo = p[idx[outil]][:, :, 2] if outil.any() else np.zeros((1, 4))
        TOOLZ.append([zo.min(), zo.max()])
    out['el_xyz'] = np.array(P)            # (nframe, nel, 4, 3)
    out['el_bulkD'] = np.array(BD)
    out['el_vonMises'] = np.array(VM)
    out['tool_z'] = np.array(TOOLZ, dtype=np.float32)
    out['axe'] = np.array([cx, cy, zmax], dtype=np.float32)

    # ---- joints : la DERNIERE frame porte toute la chronologie -------------
    jf = sorted(glob.glob(os.path.join(src, 'fdem3d_joints_0*.vtu')))
    if jf:
        pj, Aj = lire_vtu(jf[-1])
        cj, idxj = centres(pj, Aj, 3)
        tb = Aj.get('tBreak', np.full(len(cj), -1.0))
        rompu = tb >= 0
        out['jt_xyz'] = pj[idxj[rompu]].astype(np.float32)   # (njt, 3, 3)
        out['jt_tBreak'] = tb[rompu].astype(np.float32)
        out['jt_mode'] = Aj.get('breakMode', np.zeros(len(cj)))[rompu].astype(np.int8)
        out['jt_frame'] = np.array([len(jf) - 1])
        print(f"joints rompus : {int(rompu.sum())} / {len(cj)} "
              f"(derniere frame, tBreak porte la chronologie)")

    np.savez_compressed(dst, **out)
    mo = os.path.getsize(dst) / 1e6
    print(f"ecrit : {dst}  ({mo:.1f} Mo, {len(ef)} frames)")


if __name__ == '__main__':
    main()
