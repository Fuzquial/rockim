# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# imp_lib.py — lecture et depouillement des impacts a insert unique
# (spec 005, WP5). Partage par fig_impact.py et gif_impact.py.
#
# Les VTU de rockim sont en ASCII (VtkWriter) : un parseur minimal suffit,
# et il ne depend d aucune bibliotheque VTK.
# ---------------------------------------------------------------------------
import csv
import glob
import io
import os
import re

import numpy as np

Z_SURF = 0.15                # sommet de la roche (maillage translate a lo=0)
CX = CY = 0.125              # axe de l impact : le solveur translate le
                             # maillage a lo = 0, l axe passe au centre


def read_vtu(path):
    """Points, connectivite triangles, et champs par cellule d un VTU ascii."""
    s = io.open(path, encoding="utf-8", errors="ignore").read()

    def arr(name, dtype=float):
        m = re.search(r'Name="%s"[^>]*>\s*(.*?)\s*</DataArray>' % name, s,
                      re.S)
        if not m:
            return None
        return np.fromstring(m.group(1), sep=" ", dtype=float).astype(dtype)

    pts = arr("Points") if 'Name="Points"' in s else None
    if pts is None:                        # points sans Name : premier bloc
        m = re.search(r'<Points>.*?<DataArray[^>]*>\s*(.*?)\s*</DataArray>',
                      s, re.S)
        pts = np.fromstring(m.group(1), sep=" ")
    pts = pts.reshape(-1, 3)
    con = arr("connectivity", int).reshape(-1, 3)
    fields = {}
    for m in re.finditer(r'<DataArray[^>]*Name="([A-Za-z]\w*)"[^>]*>\s*'
                         r'(.*?)\s*</DataArray>', s, re.S):
        nm = m.group(1)
        if nm in ("connectivity", "offsets", "types"):
            continue
        v = np.fromstring(m.group(2), sep=" ")
        if v.size == con.shape[0]:
            fields[nm] = v
    return pts, con, fields


def joints_frame(run, k):
    return os.path.join(run, "fdem3d_joints_%04d.vtu" % k)


def frames_of(run):
    fs = sorted(glob.glob(os.path.join(run, "fdem3d_joints_[0-9]*.vtu")))
    return [int(re.search(r"_(\d+)\.vtu$", f).group(1)) for f in fs]


def frame_times(run):
    t = {}
    with io.open(os.path.join(run, "frames.csv")) as f:
        for row in csv.DictReader(f):
            t[int(row["frame"])] = float(row["t"])
    return t


def history(run):
    r = list(csv.DictReader(io.open(os.path.join(run, "history.csv"))))
    return {k: np.array([float(x[k]) for x in r]) for k in r[0]}


def broken(pts, con, f):
    """Faces ROMPUES (D = 1, plus bondees) : sommets (n,3,3), centroides,
    normales unitaires et mode de rupture (1 = traction, 2 = cisaillement)."""
    sel = (f["damage"] >= 0.999) & (f["bonded"] < 0.5)
    P = pts[con[sel]]
    c = P.mean(axis=1)
    n = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    n = n / np.maximum(np.linalg.norm(n, axis=1), 1e-30)[:, None]
    mode = f["breakMode"][sel] if "breakMode" in f else np.zeros(len(c))
    return c, n, mode, P


def metrics(c):
    """Les metriques morphologiques de leur fig. 8, sur les centroides c."""
    if len(c) == 0:
        return dict(radial=0.0, crater=0.0, depth=0.0, n=0)
    r = np.hypot(c[:, 0] - CX, c[:, 1] - CY)
    zrel = Z_SURF - c[:, 2]                       # profondeur sous la surface
    surf = zrel < 0.003                           # les 3 mm sous la surface
    return dict(
        radial=float(r.max()),                     # leur "radial crack length"
        crater=float(r[surf].max()) if surf.any() else 0.0,
        depth=float(zrel.max()),
        n=int(len(c)))
