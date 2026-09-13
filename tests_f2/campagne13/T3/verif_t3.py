#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# verif_t3.py — verification INDEPENDANTE (2e passe T3, campagne du 13/09) de
# la serie de maillages a train fige : docs/MAILLAGE_serie_2026-09-13.md.
# Ne genere aucun maillage, ne lance aucun run : ne LIT que les .msh du depot.
#
#   python tests_f2/campagne13/T3/verif_t3.py
#
# Critere FALSIFIANT de chaque bloc (un echec => code de retour 1) :
#  A. train identique sur les trois membres de la serie : memes coordonnees
#     (a 0 ULP), memes tetras (connectivite canonique), memes volumes ;
#  B. train DIFFERENT entre les maillages du chemin par defaut a SR != s
#     (impact_yang_s1_pose.msh et impact_yang_s2.5_pose.msh) — c'est le defaut
#     que T3 corrige : si ce test passait « identique », le motif de la tache
#     n'existerait pas ;
#  C. roche de rock137 == roche de impact_yang_s1_pose_hxt05.msh (le run s = 1
#     actuel), a 0 ULP sur les coordonnees et sur la connectivite ;
#  D. masses par corps reproduites a 1e-6 kg pres contre les valeurs ecrites
#     dans le rapport, et contre les masses IMPRIMEES par le solveur dans
#     results/yang_bench_s25_v3P.log et results/yang_adaptive_smoke.log ;
#  E. arete moyenne mediane dans la boule R 12,5 mm et h inscrit min de la
#     roche, contre les valeurs du rapport (tolerance 1e-3 mm) ;
#  F. format des fichiers fusionnes : CRLF, elements de type 4 seuls, aucun
#     noeud orphelin, six tags physiques, aucun octet de controle < 32 hors
#     tab/LF/CR dans les outils Python edites par T3.
# ---------------------------------------------------------------------------
import os
import re
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import mesh_quality as mq                                        # noqa: E402

SERIE = {
    "rock25": "meshes/impact_yang_train1_rock25_hxt.msh",
    "rock137": "meshes/impact_yang_train1_rock137_hxt.msh",
    "rock073": "meshes/impact_yang_train1_rock073_hxt.msh",
}
DEFAUT = {
    "s1_pose": "meshes/impact_yang_s1_pose.msh",
    "s25_pose": "meshes/impact_yang_s2.5_pose.msh",
    "s1_hxt05": "meshes/impact_yang_s1_pose_hxt05.msh",
}
TRAIN = ("insert", "bit", "piston", "circlip", "plate")
# valeurs ECRITES dans docs/MAILLAGE_serie_2026-09-13.md (1re passe) : ce
# script les re-mesure et refuse tout ecart au-dela de la tolerance.
ATTENDU_MASSE = {"insert": 0.064584, "bit": 1.288297, "piston": 1.058054,
                 "circlip": 0.014237, "plate": 0.190551}
ATTENDU_ROCHE = {                     # masse, arete med. boule, h min (mm)
    "rock25": (19.247244, 3.594, 0.7654),
    "rock137": (19.320000, 1.425, 0.2700),
    "rock073": (19.327251, 1.033, 0.2000),
}
ATTENDU_NTET = {"rock25": (24010, 6221), "rock137": (108667, 90878),
                "rock073": (251460, 233671)}
# masses imprimees par le solveur (grep dans les journaux du depot)
LOGS = {"results/yang_bench_s25_v3P.log": None, "results/yang_adaptive_smoke.log": None}

ok_global = True


def verdict(nom, ok, txt):
    global ok_global
    ok_global = ok_global and ok
    print("[%s] %-52s %s" % ("OK " if ok else "ECHEC", nom, txt))


def charge(rel):
    nodes, tets, phys, names = mq.read_msh(os.path.join(ROOT, rel))
    h, em, cen, V = mq.inscribed(nodes, tets)
    par = {}
    for tag in sorted(set(phys.tolist())):
        nm = names.get(tag, str(tag))
        m = phys == tag
        par[nm] = dict(n=int(m.sum()), tets=tets[m], V=float(V[m].sum()),
                       h=h[m], em=em[m], cen=cen[m])
    return nodes, par


def signature(nodes, bloc):
    """(coordonnees triees, connectivite canonique triee) d'un corps : deux
    maillages ont la meme signature si et seulement si ils ont les memes
    noeuds et les memes tetras, quelle que soit la numerotation."""
    t = bloc["tets"]
    xyz = nodes[np.unique(t.ravel())]
    o = np.lexsort((xyz[:, 2], xyz[:, 1], xyz[:, 0]))
    xyz = xyz[o]
    # connectivite en coordonnees : chaque tetra = ses 4 sommets tries
    P = nodes[t].reshape(-1, 3)
    key = np.round(P / 1e-12).astype(np.int64).reshape(len(t), 4, 3)
    ks = np.array([sorted(map(tuple, k)) for k in key]).reshape(len(t), 12)
    ks = ks[np.lexsort(ks.T[::-1])]
    return xyz, ks


def main():
    print("=" * 100)
    print("A. train identique sur les trois membres de la serie (critere falsifiant)")
    ref_nodes, ref = charge(SERIE["rock137"])
    sig_ref = {nm: signature(ref_nodes, ref[nm]) for nm in TRAIN}
    charges = {"rock137": (ref_nodes, ref)}
    for nom in ("rock25", "rock073"):
        nodes, par = charge(SERIE[nom])
        charges[nom] = (nodes, par)
        for nm in TRAIN:
            xyz, ks = signature(nodes, par[nm])
            same = (xyz.shape == sig_ref[nm][0].shape and np.array_equal(xyz, sig_ref[nm][0])
                    and ks.shape == sig_ref[nm][1].shape and np.array_equal(ks, sig_ref[nm][1]))
            dv = abs(par[nm]["V"] - ref[nm]["V"])
            verdict("%s / rock137 : %s" % (nom, nm), same and dv == 0.0,
                    "N %d = %d, noeuds %d, dV %.3e m3" % (par[nm]["n"], ref[nm]["n"], len(xyz), dv))

    print("=" * 100)
    print("B. contre-epreuve : le chemin par DEFAUT fait varier le train avec SR (DOIT differer)")
    n1, p1 = charge(DEFAUT["s1_pose"])
    n25, p25 = charge(DEFAUT["s25_pose"])
    diff = [nm for nm in TRAIN if p1[nm]["n"] != p25[nm]["n"]]
    verdict("defaut s = 1 vs s = 2,5 : train different", len(diff) == len(TRAIN),
            "N piston %d vs %d, insert %d vs %d, bit %d vs %d ; corps differents : %d/5"
            % (p1["piston"]["n"], p25["piston"]["n"], p1["insert"]["n"], p25["insert"]["n"],
               p1["bit"]["n"], p25["bit"]["n"], len(diff)))
    mp1 = p1["piston"]["V"] * 7850.0
    mp25 = p25["piston"]["V"] * 7850.0
    verdict("piston du defaut : ecart s = 2,5 / s = 1", abs(100 * (mp25 / mp1 - 1) + 26.5) < 0.6,
            "%.6f vs %.6f kg = %+.2f %% (diagnostic §4 : -26,5 %%)" % (mp25, mp1, 100 * (mp25 / mp1 - 1)))

    print("=" * 100)
    print("C. roche de rock137 == roche de impact_yang_s1_pose_hxt05.msh")
    nh, ph = charge(DEFAUT["s1_hxt05"])
    xyzA, ksA = signature(charges["rock137"][0], charges["rock137"][1]["rock"])
    xyzB, ksB = signature(nh, ph["rock"])
    same = (xyzA.shape == xyzB.shape and np.array_equal(xyzA, xyzB)
            and ksA.shape == ksB.shape and np.array_equal(ksA, ksB))
    verdict("roche identique (0 ULP)", same, "N %d / %d tetras, %d / %d noeuds"
            % (charges["rock137"][1]["rock"]["n"], ph["rock"]["n"], len(xyzA), len(xyzB)))
    dtrain = [nm for nm in TRAIN if ph[nm]["n"] != ref[nm]["n"]]
    verdict("train de hxt05 != train de la serie (attendu : neuf)", len(dtrain) > 0,
            "corps differents %d/5 (insert %d vs %d, piston %d vs %d)"
            % (len(dtrain), ph["insert"]["n"], ref["insert"]["n"], ph["piston"]["n"], ref["piston"]["n"]))

    print("=" * 100)
    print("D. masses par corps (V x rho) contre le rapport et contre les journaux du solveur")
    for nom in ("rock25", "rock137", "rock073"):
        nodes, par = charges[nom]
        for nm in TRAIN:
            m = par[nm]["V"] * mq.RHO_DEFAUT[nm]
            verdict("%s : masse %s" % (nom, nm), abs(m - ATTENDU_MASSE[nm]) < 1e-6,
                    "%.6f kg (rapport %.6f)" % (m, ATTENDU_MASSE[nm]))
        mr = par["rock"]["V"] * mq.RHO_DEFAUT["rock"]
        verdict("%s : masse roche" % nom, abs(mr - ATTENDU_ROCHE[nom][0]) < 1e-3,
                "%.6f kg (rapport %.6f)" % (mr, ATTENDU_ROCHE[nom][0]))
        ntot = sum(par[k]["n"] for k in par)
        verdict("%s : nombre de tetras" % nom,
                (ntot, par["rock"]["n"]) == ATTENDU_NTET[nom],
                "total %d (rapport %d), roche %d (rapport %d)"
                % (ntot, ATTENDU_NTET[nom][0], par["rock"]["n"], ATTENDU_NTET[nom][1]))

    # masses imprimees par le solveur dans les journaux du depot
    rex = re.compile(r"corps '(\w+)': .*masse = ([0-9.eE+-]+) kg")
    for rel, cible in (("results/yang_bench_s25_v3P.log", "s25_pose"),
                       ("results/yang_adaptive_smoke.log", "s1_pose")):
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            verdict("journal %s" % rel, False, "absent")
            continue
        with open(p, errors="replace") as f:
            txt = f.read()
        lu = {k: float(v) for k, v in rex.findall(txt)}
        _, pm = charge(DEFAUT[cible])
        pires = []
        for nm, mv in lu.items():
            if nm in pm:
                mo = pm[nm]["V"] * mq.RHO_DEFAUT[nm]
                pires.append((abs(mo / mv - 1.0), nm, mo, mv))
        pires.sort(reverse=True)
        verdict("solveur vs outil : %s (%d corps)" % (os.path.basename(rel), len(pires)),
                bool(pires) and pires[0][0] < 2e-5,
                "pire ecart %s %.6f vs %.6f = %.2e" % (pires[0][1], pires[0][2], pires[0][3], pires[0][0])
                if pires else "aucun corps lu")

    print("=" * 100)
    print("E. boule R 12,5 mm (roche) : arete moyenne mediane et h inscrit min")
    for nom in ("rock25", "rock137", "rock073"):
        nodes, par = charges[nom]
        r = np.linalg.norm(par["rock"]["cen"], axis=1)
        m = r < 0.0125
        eb = float(np.median(par["rock"]["em"][m])) * 1e3
        hmin = float(par["rock"]["h"].min()) * 1e3
        att = ATTENDU_ROCHE[nom]
        verdict("%s : boule N %d" % (nom, int(m.sum())), abs(eb - att[1]) < 1e-3 and abs(hmin - att[2]) < 1e-3,
                "arete med. %.4f mm (rapport %.3f), h min roche %.4f mm (rapport %.4f)"
                % (eb, att[1], hmin, att[2]))
    # insert : le plus petit h du train, celui qui pourrait commander dt
    hins = float(ref["insert"]["h"].min()) * 1e3
    verdict("train : h inscrit min de l'insert", abs(hins - 0.2528) < 1e-3, "%.4f mm (rapport 0,2528)" % hins)

    print("=" * 100)
    print("F. format des fichiers fusionnes et hygiene des sources")
    for nom, rel in SERIE.items():
        p = os.path.join(ROOT, rel)
        with open(p, "rb") as f:
            raw = f.read()
        crlf = raw.count(b"\r\n")
        lf = raw.count(b"\n")
        txt = raw.decode("ascii", "replace")
        types = set()
        for line in txt.split("\n"):
            q = line.split()
            if len(q) > 8 and q[0].isdigit() and q[1].isdigit() and q[2] == "2":
                types.add(q[1])
        nodes, par = charges[nom]
        utilises = set()
        for nm in par:
            utilises.update(par[nm]["tets"].ravel().tolist())
        nnodes = int(re.search(r"\$Nodes\s+(\d+)", txt).group(1))
        nphys = int(re.search(r"\$PhysicalNames\s+(\d+)", txt).group(1))
        verdict("%s : format" % nom,
                crlf == lf and types == {"4"} and len(utilises) == nnodes and nphys == 6,
                "CRLF %d/%d, types {%s}, noeuds %d tous references, %d tags physiques"
                % (crlf, lf, ",".join(sorted(types)), nnodes, nphys))
    for rel in ("tools/make_impact_mesh.py", "tools/mesh_quality.py",
                "tests_f2/campagne13/T3/verif_t3.py"):
        with open(os.path.join(ROOT, rel), "rb") as f:
            raw = f.read()
        bad = [b for b in raw if b < 32 and b not in (9, 10, 13)]
        verdict("%s : aucun octet de controle" % rel, not bad, "%d octets < 32 hors tab/LF/CR" % len(bad))

    print("=" * 100)
    print("VERDICT GLOBAL : %s" % ("TOUT OK" if ok_global else "AU MOINS UN ECHEC"))
    return 0 if ok_global else 1


if __name__ == "__main__":
    sys.exit(main())
