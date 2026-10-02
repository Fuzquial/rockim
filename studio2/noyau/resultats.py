"""Lecture d'un dossier de run : historique (y compris en cours d'écriture) et frames.

Conventions du dépôt (etude_triax_hetero/depouille.py, tools/plot_results.py) :
la contrainte et la déformation sont prises en valeur absolue, et le déviateur
est q = sigma - sigma(fin de consolidation), PAS sigma - sigma3. En fdem 2D la
fin de consolidation est un état de déformation uniaxiale (sigma_ax = nu
sigma3) ; retrancher sigma3 ferait partir la courbe à -(1 - nu) sigma3.

Le cache binaire (convertir) est celui validé au jalon J1 : toutes les frames
en tableaux bruts little-endian, lus tels quels par le navigateur.
"""
import json
import os
import re
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from . import cfg

CHAMPS = ["sigmaXX", "sigmaYY", "sigmaXY", "vonMises", "epsXX"]


# ---------------------------------------------------------------- historique
def lire_historique(dossier):
    """Colonnes de history.csv en tableaux numpy. Tolère un run en cours :
    une dernière ligne incomplète (écriture en cours) est écartée."""
    p = os.path.join(dossier, "history.csv")
    if not os.path.exists(p):
        return None
    with open(p, encoding="utf-8", errors="replace") as f:
        lignes = f.read().splitlines()
    if len(lignes) < 2:
        return None
    tete = lignes[0].split(",")
    rangs = []
    for l in lignes[1:]:
        v = l.split(",")
        if len(v) != len(tete):
            continue
        try:
            rangs.append([float(x) for x in v])
        except ValueError:
            continue
    if not rangs:
        return None
    a = np.array(rangs)
    return {n: a[:, k] for k, n in enumerate(tete)}


def dernier_temps(dossier):
    """Temps simulé de la dernière ligne complète, sans relire tout le fichier."""
    p = os.path.join(dossier, "history.csv")
    try:
        with open(p, "rb") as f:
            f.seek(0, os.SEEK_END)
            n = f.tell()
            f.seek(max(0, n - 4096))
            lignes = f.read().decode("utf-8", "replace").splitlines()
    except OSError:
        return None
    for l in reversed(lignes[:-1] if len(lignes) > 1 else lignes):
        try:
            return float(l.split(",")[0])
        except ValueError:
            continue
    return None


def delai_consolidation(dossier):
    p = os.path.join(dossier, "config_effective.cfg")
    if not os.path.exists(p):
        return 0.0
    try:
        return float(cfg.lire(p).get("pullDelay", "0").split()[0])
    except ValueError:
        return 0.0


def courbe(h, delai=0.0):
    """(epsilon axial en %, q en MPa) selon la convention du dépôt."""
    sig = np.abs(h["sigma"])
    s0 = np.interp(delai, h["t"], sig) if delai > 0 else 0.0
    return 100.0 * np.abs(h["epsGauge"]), (sig - s0) / 1e6


# ---------------------------------------------------------------- frames et cache
def frames_disponibles(dossier):
    return sorted(int(m.group(1)) for f in os.listdir(dossier)
                  for m in [re.fullmatch(r"fdem_(\d{4})\.vtu", f)] if m)


def _lire_frame(args):
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
    for c in CHAMPS + ["phase"]:
        a = cd.GetArray(c)
        if a is not None:
            out[c] = v2n(a).astype(np.float32)
    rj = vtk.vtkXMLUnstructuredGridReader()
    rj.SetFileName(os.path.join(run, "fdem_joints_%04d.vtu" % i))
    rj.Update()
    gj = rj.GetOutput()
    out["jmode"] = v2n(gj.GetCellData().GetArray("breakMode")).astype(np.uint8)
    out["jseg"] = v2n(gj.GetCells().GetConnectivityArray()).astype(np.uint32).reshape(-1, 2)
    return i, out


def convertir(run, cache):
    """Run -> cache binaire (meta.json + .bin). Voir specs/007-studio-v2/j1."""
    t0 = time.perf_counter()
    idx = frames_disponibles(run)
    with ProcessPoolExecutor(max_workers=max(1, min(len(idx), os.cpu_count() or 4))) as ex:
        F = dict(ex.map(_lire_frame, [(run, i) for i in idx]))
    F = [F[i] for i in idx]
    os.makedirs(cache, exist_ok=True)

    def ecrire(nom, tab):
        np.ascontiguousarray(tab).tofile(os.path.join(cache, nom + ".bin"))

    ecrire("pts", np.stack([f["pts"] for f in F]))
    bornes, champs = {}, [c for c in CHAMPS if c in F[0]]
    for c in champs:
        a = np.stack([f[c] for f in F])
        ecrire(c, a)
        bornes[c] = [float(np.percentile(a, 1)), float(np.percentile(a, 99))]
    ecrire("phase", F[0].get("phase", np.zeros(len(F[0]["pts"]) // 3, np.float32)))
    ecrire("jseg", F[0]["jseg"])
    ecrire("jmode", np.stack([f["jmode"] for f in F]))

    h = lire_historique(run)
    eps, q = courbe(h, delai_consolidation(run))
    ecrire("hist", np.column_stack([eps, q]).astype(np.float32))
    fc = os.path.join(run, "frames.csv")
    tf = (np.loadtxt(fc, delimiter=",", skiprows=1, ndmin=2)[:, 1] if os.path.exists(fc)
          else np.linspace(0, h["t"][-1], len(F)))[:len(F)]
    p = F[0]["pts"]
    d = cfg.lire(os.path.join(run, "config_effective.cfg")) if os.path.exists(
        os.path.join(run, "config_effective.cfg")) else {}
    meta = dict(run=os.path.basename(os.path.normpath(run)), nFrames=len(F),
                temps=[float(t) for t in tf], nVert=int(len(p)), nTri=int(len(p) // 3),
                nJoint=int(len(F[0]["jseg"])), nHist=int(len(q)),
                frameVersHist=[int(np.argmin(np.abs(h["t"] - t))) for t in tf],
                champs=champs, bornes=bornes,
                bbox=[float(p[:, 0].min()), float(p[:, 1].min()), float(p[:, 0].max()), float(p[:, 1].max())],
                sigma3_MPa=float(d.get("confiningPressure", "0").split()[0]) / 1e6,
                conversion_s=round(time.perf_counter() - t0, 2))
    with open(os.path.join(cache, "meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=1)
    return meta
