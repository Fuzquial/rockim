"""Vérification de l'extension writeStrainFields (spec 007, J3).

Trois runs du même deck (F7_disc_gbm_P020, tronqué), lancés en parallèle par la
file du noyau, un fil chacun (OMP_NUM_THREADS = 1 : bit-identique au série) :
  A  rockim_j3ref.exe (sources d'avant J3), sans la clé
  B  rockim_j3.exe,                          sans la clé
  C  rockim_j3.exe,                          avec writeStrainFields = true
Critères :
  1. A et B : TOUS les fichiers de sortie identiques octet à octet (défaut intact) ;
  2. A et C : history.csv, joints et csv finaux identiques ; VTU d'éléments identiques
     sur tous les tableaux communs (la clé n'a changé aucune trajectoire) ;
  3. C : déformation physiquement cohérente (voir controles_physiques).

    python specs/007-studio-v2/j3/verif_j3.py [--lancer]
Sans --lancer : écrit les decks et imprime le plan, ne lance rien.
"""
import copy
import filecmp
import os
import sys
import time

import numpy as np

ICI = os.path.dirname(os.path.abspath(__file__))
G1 = os.path.abspath(os.path.join(ICI, "..", "..", ".."))
sys.path[:0] = [os.path.join(G1, "studio2"), os.path.join(G1, "studio2", "tests")]

from campagne_triax_hetero import essai_du_cas, gen_decks   # noqa: E402
from noyau.file import File                                 # noqa: E402

T, FRAMES = 8.0e-4, 4
TRAVAIL = os.path.join(G1, "studio2", "_travail", "j3")


def essais():
    base = essai_du_cas(gen_decks(), "F7_disc_gbm_P020", "F7_disc_gbm", 20, 4211)
    base.sorties.T, base.sorties.frames = T, FRAMES
    sans = copy.deepcopy(base)
    sans.nom = "F7_j3_sans"
    avec = copy.deepcopy(base)
    avec.nom, avec.sorties.champs_deformation = "F7_j3_avec", True
    return sans, avec


def files():
    sans, avec = essais()
    plan = [("A", "rockim_j3ref.exe", sans), ("B", "rockim_j3.exe", sans), ("C", "rockim_j3.exe", avec)]
    out = []
    for nom, exe, e in plan:
        f = File(os.path.join(TRAVAIL, nom), [os.path.join(G1, exe)], jobs=1, fils=1, cwd=G1)
        if not f.travaux:
            f.ajouter(e)
        out.append((nom, exe, f))
    return out


def lire_vtu(p):
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy as v2n
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(p)
    r.Update()
    g = r.GetOutput()
    d = {"points": v2n(g.GetPoints().GetData())}
    for src in (g.GetCellData(), g.GetPointData()):
        for i in range(src.GetNumberOfArrays()):
            d[src.GetArrayName(i)] = v2n(src.GetArray(i))
    return d


def comparer(F):
    oa, ob, oc = (f.travaux[0]["out"] for _, _, f in F)
    noms = sorted(os.listdir(oa))
    print("\n1. A (reference) contre B (j3, sans la cle) : %d fichiers" % len(noms))
    diff_ab = [n for n in noms if not filecmp.cmp(os.path.join(oa, n), os.path.join(ob, n), shallow=False)]
    print("   differents :", diff_ab or "AUCUN")
    print("\n2. A contre C (j3, avec la cle)")
    for n in noms:
        a, c = os.path.join(oa, n), os.path.join(oc, n)
        if n.startswith("fdem_") and n.endswith(".vtu") and "joints" not in n:
            va, vc = lire_vtu(a), lire_vtu(c)
            communs = [k for k in va if k in vc]
            ecarts = [k for k in communs if not np.array_equal(va[k], vc[k])]
            print("   %-22s tableaux communs %d, differents %s, ajoutes %s"
                  % (n, len(communs), ecarts or "aucun", sorted(set(vc) - set(va))))
        elif n != "config_effective.cfg":
            print("   %-22s %s" % (n, "identique" if filecmp.cmp(a, c, shallow=False) else "DIFFERENT"))
    return oc


def controles_physiques(oc):
    """La déformation écrite doit dire la même chose que ce qu'on sait déjà."""
    import glob
    der = sorted(glob.glob(os.path.join(oc, "fdem_[0-9]*.vtu")))[-1]
    v = lire_vtu(der)
    p = v["points"][:, :2]
    c = p.reshape(-1, 3, 2).mean(1)
    print("\n3. controles physiques sur", os.path.basename(der))
    # (a) déplacement = position - position initiale
    v0 = lire_vtu(os.path.join(oc, "fdem_0000.vtu"))
    print("   (a) |displacement - (x - x0)| max = %.3g m" % np.abs(v["displacement"][:, :2] - (p - v0["points"][:, :2])).max())
    # (b) avant rupture les rotations sont petites : strainXX ~ epsXX (co-roté)
    print("   (b) |strainXX - epsXX| max = %.3g (|epsXX| max = %.3g)"
          % (np.abs(v["strainXX"] - v["epsXX"]).max(), np.abs(v["epsXX"]).max()))
    # (c) moyenne de strainYY entre les jauges ~ déformation de l'extensomètre (history)
    import csv
    h = list(csv.DictReader(open(os.path.join(oc, "history.csv"))))
    zone = (c[:, 1] > 0.018) & (c[:, 1] < 0.054)
    print("   (c) moyenne strainYY entre jauges = %.4g ; epsGauge (history) = %.4g"
          % (v["strainYY"][zone].mean(), float(h[-1]["epsGauge"])))
    # (d) sous confinement seul (avant pullDelay) : strainXX < 0 au coeur, tenseur symétrique
    print("   (d) moyenne strainXX = %.4g, strainXY = %.3g" % (v["strainXX"].mean(), v["strainXY"].mean()))


def main():
    F = files()
    print("Plan : 3 runs du deck F7_disc_gbm_P020, T = %g s, frames = %d, 1 fil chacun" % (T, FRAMES))
    for nom, exe, f in F:
        t = f.travaux[0]
        print("  %s  %-18s %s  -> %s" % (nom, exe, os.path.basename(t["deck"]), t["out"]))
    if "--lancer" not in sys.argv:
        print("\n(aucun lancement : ajouter --lancer)")
        return
    t0 = time.time()
    while any(t["etat"] in ("en_attente", "en_cours") for _, _, f in F for t in f.travaux):
        for _, _, f in F:
            f.pas()
        time.sleep(2)
    for nom, _, f in F:
        t = f.travaux[0]
        print("%s : %s (code %s) en %.0f s" % (nom, t["etat"], t["code"], (t["fin"] or 0) - (t["debut"] or 0)))
    print("duree totale %.0f s" % (time.time() - t0))
    controles_physiques(comparer(F))


if __name__ == "__main__":
    main()
