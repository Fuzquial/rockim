"""Convertit un dossier de run 2D (VTU ASCII + history.csv) en cache binaire.

Prototype jetable de la spec 007 (jalon J1), commun aux prototypes web et QML.
Le cache est un dossier de fichiers bruts little-endian plus un meta.json :
le navigateur les charge par fetch() en ArrayBuffer, numpy par np.fromfile.

    python cache_run.py <dossier_run> [dossier_cache]

Fichiers écrits :
    meta.json           dimensions, temps des frames, bornes des champs, colonnes
    pts.bin             float32 [F, nVert, 2]   positions déformées (m)
    <champ>.bin         float32 [F, nTri]       sigmaXX, sigmaYY, sigmaXY, vonMises, epsXX
    phase.bin           float32 [nTri]          phase minérale (constante)
    jseg.bin            uint32  [nJoint, 2]     sommets de chaque joint
    jmode.bin           uint8   [F, nJoint]     breakMode (0 intact, 1 traction, 2 cisaillement, 4 pré-rompu)
    hist.bin            float32 [nRow, 2]       (epsGauge en %, q en MPa)
"""
import json
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

CHAMPS = ["sigmaXX", "sigmaYY", "sigmaXY", "vonMises", "epsXX"]


def _lire_frame(args):
    """Lit une frame (éléments + joints). Exécuté dans un processus fils."""
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy as v2n

    run, i = args
    out = {}
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(os.path.join(run, "fdem_%04d.vtu" % i))
    r.Update()
    g = r.GetOutput()
    out["pts"] = v2n(g.GetPoints().GetData())[:, :2].astype(np.float32)
    conn = v2n(g.GetCells().GetConnectivityArray())
    # Le solveur écrit des nœuds propres à chaque triangle : sommet k -> triangle k // 3.
    if not np.array_equal(conn, np.arange(conn.size)):
        raise ValueError("connectivité non triviale dans fdem_%04d.vtu" % i)
    cd = g.GetCellData()
    for c in CHAMPS + (["phase"] if i == 0 else []):
        out[c] = v2n(cd.GetArray(c)).astype(np.float32)

    rj = vtk.vtkXMLUnstructuredGridReader()
    rj.SetFileName(os.path.join(run, "fdem_joints_%04d.vtu" % i))
    rj.Update()
    gj = rj.GetOutput()
    out["jmode"] = v2n(gj.GetCellData().GetArray("breakMode")).astype(np.uint8)
    if i == 0:
        out["jseg"] = v2n(gj.GetCells().GetConnectivityArray()).astype(np.uint32).reshape(-1, 2)
    return i, out


def _cle_cfg(run, cle, defaut=None):
    with open(os.path.join(run, "config_effective.cfg"), encoding="utf-8", errors="replace") as f:
        for ligne in f:
            m = re.match(r"\s*%s\s*=\s*([^\s#]+)" % re.escape(cle), ligne)
            if m:
                return m.group(1)
    return defaut


def _historique(run):
    """Courbe q-epsilon, convention du dépôt : q = sigma - sigma(fin de consolidation)."""
    with open(os.path.join(run, "history.csv"), encoding="utf-8") as f:
        tete = f.readline().strip().split(",")
    h = np.loadtxt(os.path.join(run, "history.csv"), delimiter=",", skiprows=1, ndmin=2)
    col = {n: h[:, k] for k, n in enumerate(tete)}
    delai = float(_cle_cfg(run, "pullDelay", "0"))
    q = col["sigma"] - np.interp(delai, col["t"], col["sigma"])
    return np.column_stack([col["epsGauge"] * 100.0, q / 1e6]).astype(np.float32), col["t"]


def convertir(run, cache):
    t0 = time.perf_counter()
    nf = len([f for f in os.listdir(run) if re.fullmatch(r"fdem_\d{4}\.vtu", f)])
    with ProcessPoolExecutor(max_workers=min(nf, os.cpu_count() or 4)) as ex:
        frames = dict(ex.map(_lire_frame, [(run, i) for i in range(nf)]))
    os.makedirs(cache, exist_ok=True)

    def ecrire(nom, tab):
        np.ascontiguousarray(tab).tofile(os.path.join(cache, nom + ".bin"))

    ecrire("pts", np.stack([frames[i]["pts"] for i in range(nf)]))
    bornes = {}
    for c in CHAMPS:
        a = np.stack([frames[i][c] for i in range(nf)])
        ecrire(c, a)
        # Bornes robustes (percentiles 1-99) : une poignée d'éléments isolés ne doit
        # pas écraser l'échelle de couleur.
        bornes[c] = [float(np.percentile(a, 1)), float(np.percentile(a, 99))]
    ecrire("phase", frames[0]["phase"])
    ecrire("jseg", frames[0]["jseg"])
    ecrire("jmode", np.stack([frames[i]["jmode"] for i in range(nf)]))
    hist, th = _historique(run)
    ecrire("hist", hist)

    tf = np.loadtxt(os.path.join(run, "frames.csv"), delimiter=",", skiprows=1, ndmin=2)[:, 1]
    # Pour chaque frame, l'indice de la ligne d'historique la plus proche en temps :
    # c'est ce qui synchronise le curseur de la courbe avec la vue des champs.
    lien = [int(np.argmin(np.abs(th - t))) for t in tf]
    p = frames[0]["pts"]
    meta = {
        "run": os.path.basename(os.path.normpath(run)),
        "nFrames": nf,
        "temps": [float(t) for t in tf],
        "nVert": int(p.shape[0]),
        "nTri": int(p.shape[0] // 3),
        "nJoint": int(frames[0]["jseg"].shape[0]),
        "nHist": int(hist.shape[0]),
        "frameVersHist": lien,
        "champs": CHAMPS,
        "bornes": bornes,
        "bbox": [float(p[:, 0].min()), float(p[:, 1].min()), float(p[:, 0].max()), float(p[:, 1].max())],
        "sigma3_MPa": float(_cle_cfg(run, "confiningPressure", "0")) / 1e6,
        "conversion_s": round(time.perf_counter() - t0, 2),
    }
    with open(os.path.join(cache, "meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=1)
    return meta


if __name__ == "__main__":
    run = sys.argv[1]
    ici = os.path.dirname(os.path.abspath(__file__))
    cache = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        ici, "..", "_cache", os.path.basename(os.path.normpath(run)))
    m = convertir(run, cache)
    print("%s : %d frames, %d triangles, %d joints, conversion %.2f s -> %s"
          % (m["run"], m["nFrames"], m["nTri"], m["nJoint"], m["conversion_s"], os.path.abspath(cache)))
