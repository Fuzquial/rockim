#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# export_abaqus.py — exporte un maillage rockim (frame 0) vers Abaqus .inp,
# pour la validation croisée rockim <-> Abaqus/Explicit+VUMAT sur maillage
# IDENTIQUE (memes noeuds, memes tets C3D4).
#
#   python export_abaqus.py <run_dir | frame0.vtu> <sortie.inp> [options]
#
# Options :
#   --units mm|si   mm (defaut) = convention these mm-t-s-MPa : coordonnees
#                   multipliees par 1000 ; si = metres inchanges.
#   --no-field      ne pas exporter le champ ftScale meme s'il existe.
#
# Entree : le VTU de la frame 0 (fem3d_0000.vtu ou fdem3d_0000.vtu).
#   * fem3d : noeuds PARTAGES — export direct.
#   * fdem3d : noeuds DUPLIQUES par tet — a la frame 0 les copies sont
#     co-localisees, l'exporteur les SOUDE par coordonnees (le continuum
#     geometrique, y compris un maillage Voronoi/GBM, devient un maillage
#     C3D4 classique ; les joints cohesifs ne sont PAS exportes — cote
#     Abaqus la rupture vient de la loi/VUMAT, c'est le principe meme de la
#     comparaison implicite/explicite).
#
# Sorties :
#   <sortie.inp>            *NODE, *ELEMENT C3D4, NSET NALL/NBOTTOM/NTOP,
#                           ELSET EALL (+ ELSET_PHASE<i> si 'phase' present)
#   <sortie>_field.csv      elem, cx, cy, cz, ftScale   (si champ present)
#   <sortie>_field.inp      *INITIAL CONDITIONS, TYPE=FIELD, VAR=1 (valeurs
#                           nodales = moyenne des elements incidents) — a
#                           inclure dans le deck pour le mecanisme FIELD des
#                           VUMAT (saksala2011 SDV15/16, DP-DFH FIELD).
#
# Autocontroles : volume total des tets reordonnes (jacobien positif,
# convention C3D4) compare au volume VTU ; comptes noeuds/elements.
# Rappel these : point decimal partout, unites mm-t-s-MPa (locale FR : ce
# script n'ecrit que des points).
# ---------------------------------------------------------------------------
import os
import re
import sys

import numpy as np


def read_vtu(path):
    txt = open(path).read()
    m = re.search(r"<Points>.*?<DataArray[^>]*>(.*?)</DataArray>", txt, re.S)
    pts = np.fromstring(m.group(1), sep=" ").reshape(-1, 3)
    conn = np.fromstring(
        re.search(r'Name="connectivity"[^>]*>(.*?)</DataArray>', txt, re.S)
        .group(1), sep=" ").astype(int)
    offs = np.fromstring(
        re.search(r'Name="offsets"[^>]*>(.*?)</DataArray>', txt, re.S)
        .group(1), sep=" ").astype(int)
    if not np.all(np.diff(np.concatenate([[0], offs])) == 4):
        raise SystemExit("export_abaqus: le VTU n'est pas un maillage de "
                         "tets (offsets != 4) — exporter la frame 0 d'un "
                         "run fem3d/fdem3d")
    tets = conn.reshape(-1, 4)
    cell = {}
    for name in ("ftScale", "phase", "grain"):
        mm = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
        if mm:
            cell[name] = np.fromstring(mm.group(1), sep=" ")
    return pts, tets, cell


def weld(pts, tets):
    """Soude les noeuds co-localises (fdem3d frame 0). Tolerance : 1e-9 de
    la plus grande dimension — les copies sont exactement co-localisees a
    la frame 0, la tolerance ne fait qu'absorber l'ecriture ASCII."""
    scale = max(pts.max(axis=0) - pts.min(axis=0))
    q = np.round(pts / (1e-9 * scale)).astype(np.int64)
    _, first, inv = np.unique(q, axis=0, return_index=True,
                              return_inverse=True)
    return pts[first], inv[tets]


def main(argv):
    args = [a for a in argv if not a.startswith("--")]
    units = "mm"
    for a in argv:
        if a.startswith("--units"):
            units = a.split("=", 1)[1] if "=" in a else "mm"
    noField = "--no-field" in argv
    if len(args) != 2:
        raise SystemExit(__doc__ or "usage: export_abaqus.py <run|vtu> <inp>")
    src, out = args
    if os.path.isdir(src):
        for base in ("fem3d_0000.vtu", "fdem3d_0000.vtu"):
            p = os.path.join(src, base)
            if os.path.exists(p):
                src = p
                break
        else:
            raise SystemExit("export_abaqus: aucune frame 0 (fem3d_0000.vtu "
                             "/ fdem3d_0000.vtu) dans " + src)

    pts, tets, cell = read_vtu(src)
    n0, e0 = len(pts), len(tets)
    duplicated = len(pts) == 4 * len(tets)     # fdem3d: un noeud par coin
    if duplicated:
        pts, tets = weld(pts, tets)

    # noeuds references seulement, renumerotes 1-based
    used = np.unique(tets)
    remap = -np.ones(len(pts), dtype=int)
    remap[used] = np.arange(len(used))
    pts = pts[used]
    tets = remap[tets]

    # jacobien positif dans la convention C3D4 (base 1-2-3 vue de 4)
    a = pts[tets[:, 1]] - pts[tets[:, 0]]
    b = pts[tets[:, 2]] - pts[tets[:, 0]]
    c = pts[tets[:, 3]] - pts[tets[:, 0]]
    v6 = np.einsum("ij,ij->i", np.cross(a, b), c)
    flip = v6 < 0
    tets[flip] = tets[flip][:, [0, 2, 1, 3]]
    volSI = np.abs(v6).sum() / 6.0

    scale = 1000.0 if units == "mm" else 1.0
    P = pts * scale

    zmin, zmax = P[:, 2].min(), P[:, 2].max()
    tol = 1e-6 * (zmax - zmin)
    nbot = np.where(P[:, 2] < zmin + tol)[0]
    ntop = np.where(P[:, 2] > zmax - tol)[0]

    with open(out, "w") as f:
        f.write("*HEADING\n")
        f.write("rockim mesh export: %s\n" % os.path.basename(src))
        f.write("** units: %s (%s), %d nodes, %d C3D4%s\n"
                % ("mm-t-s-MPa" if units == "mm" else "SI", units,
                   len(P), len(tets),
                   ", fdem3d welded" if duplicated else ""))
        f.write("*NODE\n")
        for i, p in enumerate(P, 1):
            f.write("%d, %.9g, %.9g, %.9g\n" % (i, p[0], p[1], p[2]))
        f.write("*ELEMENT, TYPE=C3D4, ELSET=EALL\n")
        for e, t in enumerate(tets, 1):
            f.write("%d, %d, %d, %d, %d\n"
                    % (e, t[0] + 1, t[1] + 1, t[2] + 1, t[3] + 1))
        for name, idx in (("NBOTTOM", nbot), ("NTOP", ntop)):
            f.write("*NSET, NSET=%s\n" % name)
            for k in range(0, len(idx), 8):
                f.write(", ".join(str(i + 1) for i in idx[k:k + 8]) + "\n")
        f.write("*NSET, NSET=NALL, GENERATE\n1, %d, 1\n" % len(P))
        if "phase" in cell:
            ph = cell["phase"].astype(int)
            for p in np.unique(ph):
                f.write("*ELSET, ELSET=EPHASE%d\n" % p)
                ids = np.where(ph == p)[0] + 1
                for k in range(0, len(ids), 8):
                    f.write(", ".join(map(str, ids[k:k + 8])) + "\n")

    msg = ["export_abaqus: %d noeuds, %d C3D4 -> %s" % (len(P), len(tets), out),
           "  volume tets = %.9g m^3 (VTU %d pts / %d tets%s)"
           % (volSI, n0, e0, ", soudes" if duplicated else "")]

    fts = cell.get("ftScale")
    if fts is not None and not noField and not np.allclose(fts, fts[0]):
        cen = P[tets].mean(axis=1)
        csv = os.path.splitext(out)[0] + "_field.csv"
        with open(csv, "w") as f:
            f.write("elem,cx,cy,cz,ftScale\n")
            for e in range(len(tets)):
                f.write("%d,%.9g,%.9g,%.9g,%.9g\n"
                        % (e + 1, cen[e, 0], cen[e, 1], cen[e, 2], fts[e]))
        # champ nodal = moyenne des elements incidents (mecanisme FIELD :
        # Abaqus interpole les valeurs NODALES vers les points materiels)
        acc = np.zeros(len(P))
        cnt = np.zeros(len(P))
        for e in range(len(tets)):
            for nId in tets[e]:
                acc[nId] += fts[e]
                cnt[nId] += 1
        nodal = acc / np.maximum(cnt, 1)
        finp = os.path.splitext(out)[0] + "_field.inp"
        with open(finp, "w") as f:
            f.write("** rockim ftScale -> Abaqus predefined field VAR=1\n")
            f.write("*INITIAL CONDITIONS, TYPE=FIELD, VARIABLE=1\n")
            for i, v in enumerate(nodal, 1):
                f.write("%d, %.9g\n" % (i, v))
        msg.append("  champ ftScale: %s + %s (min/moy/max %.3f/%.3f/%.3f)"
                   % (os.path.basename(csv), os.path.basename(finp),
                      fts.min(), fts.mean(), fts.max()))
    print("\n".join(msg))


if __name__ == "__main__":
    main(sys.argv[1:])
