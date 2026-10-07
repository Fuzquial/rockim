#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# gbm_neper_like.py — generateur de microstructure GBM CYLINDRIQUE pour rockim
# (`mode = fem3d`, `mesh = file`, volumes physiques nommes = phases minerales).
#
# Remplace la tessellation INTERNE de rockim, qui souffre de deux defauts
# mesures le 2026-09-06 :
#   (a) ses eclats de bord ont une plus petite arete qui varie d'un facteur 90
#       d'un tirage a l'autre a taille de grain constante — donc un pas de
#       temps critique qui varie de 0,077 a 6,7 ns, ce qui rend une campagne
#       incomparable d'un germe a l'autre ;
#   (b) son maillage intragranulaire par defaut est un EVENTAIL depuis le
#       centre du grain : les aretes RAYONNENT, ce que la litterature FDEM
#       proscrit pour tout essai dont la reponse est un TRAJET de fissure
#       (le trajet suit les aretes, donc le rayon).
#
# Ici : germes en DISQUE DE POISSON (pas de cellule mince), cellules obtenues
# comme INTERSECTION DE DEMI-ESPACES (donc jamais ouvertes, voir §2), decoupe
# par le cylindre incluse dans la meme intersection, complexe de faces PARTAGE
# (donc maillage conforme sans booleen), maillage Delaunay NON STRUCTURE de
# Gmsh (aucun eventail), export MSH 2.2 ASCII verifie.
#
# Unites : SI (metres), comme les decks rockim (W = 0.016 ...).
#
#   python gbm_uniaxial/gbm_neper_like.py --diametre 8e-3 --hauteur 16e-3 \
#          --grain 4e-3 --maille 0.8e-3 --germe 1 --sortie cyl_test.msh
#
# ---------------------------------------------------------------------------
# METHODE, ET LES PIEGES QU'ELLE TRAITE
#
# §1 GERMES — DISQUE DE POISSON.  Un tirage uniforme met des germes a distance
#     arbitrairement petite : la face bissectrice du couple est alors minuscule
#     et ses aretes aussi, ce qui est exactement le defaut (a). On impose donc
#     une distance minimale r entre germes (echantillonnage par disque de
#     Poisson, grille de hachage 3x3x3, r relaxe de 7 % quand le tirage sature)
#     avec r = --minDist x (V/N)^(1/3), N fixe par la taille de grain visee.
#
# §2 BORD — CELLULES OUVERTES.  scipy.spatial.Voronoi sur des germes non
#     periodiques rend des regions OUVERTES au bord (sommet d'indice -1) :
#     elles n'ont pas de volume, donc pas de solide. Traitement retenu, celui
#     de Neper : on n'utilise JAMAIS les regions de Voronoi. Une cellule de
#     Voronoi est par definition l'intersection des demi-espaces
#         2 (q - p) . x <= |q|^2 - |p|^2      pour tout autre germe q ,
#     et on AJOUTE a cette liste les demi-espaces du domaine (les M facettes
#     du cylindre et les deux plans z = 0 et z = H). L'intersection est alors
#     bornee par construction : AUCUNE cellule ne peut etre ouverte, et la
#     decoupe par l'eprouvette est faite dans la meme operation (scipy
#     HalfspaceIntersection). L'option --miroirs ajoute en plus les germes
#     images (reflexion radiale sur la paroi, reflexion sur z = 0 et z = H) :
#     les grains de bord sont alors termines par une face plane au lieu d'etre
#     tronques ; c'est un CHOIX de microstructure, pas une necessite numerique.
#
# §3 VOISINS — CERTIFICAT.  On ne prend pas tous les germes mais les K plus
#     proches (KD-tree). C'est exact sous certificat : si d_K est la distance
#     au K-ieme voisin retenu et rho le rayon de la cellule obtenue, alors tout
#     germe omis est a distance >= d_K, et sa bissectrice ne peut couper la
#     cellule des que d_K >= 2 rho (pour x dans la cellule, |x - q| >= d_K -
#     rho >= rho >= |x - p|). On VERIFIE ce certificat cellule par cellule et
#     on double K tant qu'il n'est pas satisfait.
#
# §4 PAROI — FACETTAGE.  Le cylindre est approche par M plans. La distance des
#     plans a l'axe est prise EGALISATRICE D'AIRE, a = R sqrt(pi / (M tan(pi/M))) :
#     le polygone a exactement l'aire du disque, donc le volume de reference
#     est exact et le controle du §5 n'est pas biaise. M est choisi pour que
#     l'ecart radial maximal (fleche vers l'interieur a - R... et bosse
#     a/cos(pi/M) - R vers l'exterieur) reste sous --tolParoi x maille.
#
# §5 CONFORMITE — NI TROU NI RECOUVREMENT.  Les cellules etant des
#     intersections de demi-espaces sur la MEME liste de plans, deux cellules
#     voisines partagent exactement le plan bissecteur : le pavage est un
#     partitionnement exact. On le VERIFIE trois fois :
#       - somme des volumes de cellules / volume du cylindre facette ;
#       - chaque face interne est vue par EXACTEMENT deux cellules, chaque
#         arete d'une cellule par exactement deux de ses faces (fermeture) ;
#       - somme des volumes de tetraedres apres maillage.
#     Le maillage est conforme sans booleen OCC parce que chaque face est
#     construite UNE SEULE FOIS et referencee par ses deux cellules (les
#     sommets sont fusionnes globalement, cf. --fusion). Aucun `fragment`
#     n'est donc necessaire — c'est plus sur et sans commune mesure plus
#     rapide que de laisser OCC recoller des centaines de polyedres.
#
# §6 REGULARISATION.  Les sommets plus proches que --fusion x maille sont
#     fusionnes (les aretes plus courtes que la maille ne sont de toute facon
#     pas representables : elles ne produisent que des tetraedres plats, donc
#     un pas de temps effondre). Si la fusion casse la fermeture d'une
#     cellule, la tolerance est divisee par deux et on recommence, jusqu'a 0
#     ou la construction est exacte. On imprime la plus petite arete avant et
#     apres.
# ---------------------------------------------------------------------------

import argparse
import json
import math
import os
import sys
import time

import numpy as np
from scipy.spatial import HalfspaceIntersection, cKDTree


# ===========================================================================
# 1. GERMES — disque de Poisson dans le cylindre
# ===========================================================================
def germes_poisson(R, H, dGrain, rng, minDist=0.80, marge=0.0, relax=0.93,
                   echecsMax=2000):
    """Germes a distance minimale r dans le cylindre d'axe (R, R), hauteur H.

    N est fixe par la taille de grain VISEE (diametre de la sphere de meme
    volume) : N = V / (pi/6 dGrain^3). r part de minDist x (V/N)^(1/3) et est
    relaxe de `relax` quand le tirage sature, de sorte que N germes sont
    toujours places (le generateur est donc reproductible en NOMBRE).
    """
    V = math.pi * R * R * H
    N = max(1, int(round(V / (math.pi / 6.0 * dGrain ** 3))))
    r = minDist * (V / N) ** (1.0 / 3.0)
    m = marge * r
    pts = []
    grille = {}

    def cle(p, c):
        return (int(math.floor(p[0] / c)), int(math.floor(p[1] / c)),
                int(math.floor(p[2] / c)))

    def reconstruire(c):
        g = {}
        for i, p in enumerate(pts):
            g.setdefault(cle(p, c), []).append(i)
        return g

    echecs = 0
    while len(pts) < N:
        u = rng.random()
        rad = (R - m) * math.sqrt(u)
        th = rng.random() * 2.0 * math.pi
        z = m + rng.random() * max(H - 2.0 * m, 1e-12)
        p = np.array([R + rad * math.cos(th), R + rad * math.sin(th), z])
        k0 = cle(p, r)
        ok = True
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grille.get((k0[0] + dx, k0[1] + dy, k0[2] + dz), ()):
                        if np.linalg.norm(p - pts[j]) < r:
                            ok = False
                            break
                    if not ok:
                        break
                if not ok:
                    break
            if not ok:
                break
        if ok:
            grille.setdefault(k0, []).append(len(pts))
            pts.append(p)
            echecs = 0
        else:
            echecs += 1
            if echecs > echecsMax:
                r *= relax
                m = marge * r
                grille = reconstruire(r)
                echecs = 0
    return np.array(pts), N, r


# ===========================================================================
# 2. PLANS DU DOMAINE — cylindre facette (§4)
# ===========================================================================
def facettes_cylindre(R, H, tol):
    """M facettes, distance a l'axe egalisatrice d'aire ; M tel que l'ecart
    radial max (interieur et exterieur) soit sous `tol`. Renvoie
    (M, a, ecartInt, ecartExt)."""
    for M in range(8, 1025):
        a = R * math.sqrt(math.pi / (M * math.tan(math.pi / M)))
        eInt = R - a                       # milieu de corde, vers l'interieur
        eExt = a / math.cos(math.pi / M) - R   # sommet, vers l'exterieur
        if max(abs(eInt), abs(eExt)) <= tol:
            return M, a, eInt, eExt
    M = 1024
    a = R * math.sqrt(math.pi / (M * math.tan(math.pi / M)))
    return M, a, R - a, a / math.cos(math.pi / M) - R


def plans_domaine(R, H, M, a):
    """Demi-espaces A.x + b <= 0 du cylindre facette (axe en (R, R))."""
    A, b, nom = [], [], []
    for j in range(M):
        th = 2.0 * math.pi * (j + 0.5) / M
        n = np.array([math.cos(th), math.sin(th), 0.0])
        A.append(n)
        b.append(-(a + n[0] * R + n[1] * R))
        nom.append(("paroi", j))
    A.append(np.array([0.0, 0.0, -1.0])); b.append(0.0); nom.append(("bas", 0))
    A.append(np.array([0.0, 0.0, 1.0])); b.append(-H); nom.append(("haut", 0))
    return np.array(A), np.array(b), nom


def germes_miroirs(P, R, H):
    """Images des germes : reflexion radiale sur la paroi + sur z = 0 et z = H
    (§2, option --miroirs). Ne generent aucune cellule, seulement des plans."""
    c = np.array([R, R, 0.0])
    d = P - c
    d[:, 2] = 0.0
    rho = np.linalg.norm(d[:, :2], axis=1)
    rho = np.maximum(rho, 1e-12)
    fac = (2.0 * R / rho - 1.0)[:, None]
    lat = np.column_stack([c[0] + fac[:, 0] * d[:, 0],
                           c[1] + fac[:, 0] * d[:, 1], P[:, 2]])
    bas = P.copy(); bas[:, 2] = -P[:, 2]
    haut = P.copy(); haut[:, 2] = 2.0 * H - P[:, 2]
    return np.vstack([lat, bas, haut])


# ===========================================================================
# 3. CELLULES — intersection de demi-espaces + certificat (§2, §3)
# ===========================================================================
def cellules(P, Pall, Adom, bdom, nomdom, K0=48, S=1.0):
    """Pour chaque germe : sommets, faces (clef -> sommets), volume.

    Une face est identifiee par sa clef : ('g', min(i,j), max(i,j)) pour une
    face interne, ('d', type, indice) pour une facette du domaine. Les deux
    cellules d'une face interne portent la MEME clef : c'est ce qui rend le
    complexe partage (§5).
    """
    # Qhull juge la faisabilite du point interieur en ABSOLU : en metres, un
    # germe a 5 um d'un plan de paroi est refuse (QH6023). On travaille donc
    # en longueurs REDUITES (x <- S x, S = 1/R) et on redimensionne ensuite ;
    # les plans unitaires deviennent n.y <= S d, donc b <- S b.
    P = P * S
    Pall = Pall * S
    bdom = bdom * S
    tree = cKDTree(Pall)
    n, Nall = len(P), len(Pall)
    out = []
    Kmax = 0
    for i in range(n):
        p = P[i]
        K = K0
        while True:
            kk = int(min(K + 1, Nall))
            dist, idx = tree.query(p, k=kk)
            dist = np.atleast_1d(dist); idx = np.atleast_1d(idx)
            keep = idx != i
            vois = idx[keep]
            dvois = dist[keep]
            Q = Pall[vois]
            Ab = 2.0 * (Q - p)
            bb = (p @ p) - np.einsum('ij,ij->i', Q, Q)
            A = np.vstack([Ab, Adom])
            b = np.concatenate([bb, bdom])
            hs = HalfspaceIntersection(np.column_stack([A, b]), p)
            Vx = hs.intersections
            rho = float(np.max(np.linalg.norm(Vx - p, axis=1)))
            dK = float(dvois.max()) if len(dvois) else float('inf')
            if dK >= 2.0 * rho - 1e-15 or kk >= Nall:
                Kmax = max(Kmax, kk)
                break
            K *= 2
        # faces : halfspace -> sommets (dual_facets[k] = plans du sommet k)
        nv = len(vois)
        parPlan = {}
        for k, planes in enumerate(hs.dual_facets):
            for s in planes:
                parPlan.setdefault(int(s), []).append(k)
        faces = {}
        for s, ks in parPlan.items():
            if len(ks) < 3:
                continue
            if s < nv:
                j = int(vois[s])
                clef = ('g', min(i, j), max(i, j)) if j < len(P) else ('m', i, j)
                nrm = A[s] / np.linalg.norm(A[s])
            else:
                t, q = nomdom[s - nv]
                clef = ('d', t, q)
                nrm = A[s] / np.linalg.norm(A[s])
            faces[clef] = (np.array(ks, dtype=int), nrm)
        out.append(dict(seed=p / S, V=Vx / S, faces=faces))
        if (i + 1) % 200 == 0:
            print(f"    ... {i + 1}/{n} cellules", flush=True)
    return out, Kmax


# ===========================================================================
# 4. COMPLEXE PARTAGE — fusion des sommets, faces uniques, fermeture (§5, §6)
# ===========================================================================
def _ordonne(poly, nrm):
    """Ordonne les sommets d'une face plane par angle autour de son centre."""
    c = poly.mean(axis=0)
    e1 = np.array([1.0, 0.0, 0.0])
    if abs(nrm @ e1) > 0.9:
        e1 = np.array([0.0, 1.0, 0.0])
    u = np.cross(nrm, e1); u /= np.linalg.norm(u)
    v = np.cross(nrm, u)
    d = poly - c
    return np.argsort(np.arctan2(d @ v, d @ u))


def complexe(cells, tolFusion, echelle):
    """Fusionne les sommets, construit les faces uniques, verifie la fermeture.

    Renvoie (X, faces, cellFaces, diag) ou faces[clef] = liste ordonnee
    d'indices de sommets, cellFaces[i] = [(clef, signe), ...].
    """
    brut, off = [], []
    for c in cells:
        off.append(len(brut))
        brut.extend(c['V'])
    brut = np.array(brut)
    tol = max(tolFusion, 1e-12 * echelle)
    # union-find sur les paires a distance < tol
    par = list(range(len(brut)))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for (aa, bb) in cKDTree(brut).query_pairs(tol):
        ra, rb = find(aa), find(bb)
        if ra != rb:
            par[ra] = rb
    grp = {}
    for k in range(len(brut)):
        grp.setdefault(find(k), []).append(k)
    ids = {}
    X, diam = [], 0.0
    for g, membres in grp.items():
        pts = brut[membres]
        X.append(pts.mean(axis=0))
        if len(membres) > 1:
            diam = max(diam, float(np.max(np.ptp(pts, axis=0))))
        for m in membres:
            ids[m] = len(X) - 1
    X = np.array(X)

    # faces : union des deux cotes (robuste a un sommet vu d'un seul cote)
    brutFaces, nrmFaces = {}, {}
    for i, c in enumerate(cells):
        for clef, (ks, nrm) in c['faces'].items():
            s = brutFaces.setdefault(clef, set())
            for k in ks:
                s.add(ids[off[i] + int(k)])
            nrmFaces.setdefault(clef, nrm)
    faces = {}
    for clef, s in brutFaces.items():
        idx = np.array(sorted(s), dtype=int)
        if len(idx) < 3:
            continue
        o = _ordonne(X[idx], nrmFaces[clef])
        faces[clef] = idx[o]

    # fermeture : chaque arete d'une cellule vue par exactement 2 de ses faces
    cellFaces, ok = [], True
    for i, c in enumerate(cells):
        mesFaces = [k for k in c['faces'] if k in faces]
        aretes = {}
        for clef in mesFaces:
            pl = faces[clef]
            for q in range(len(pl)):
                e = (min(pl[q], pl[(q + 1) % len(pl)]),
                     max(pl[q], pl[(q + 1) % len(pl)]))
                aretes[e] = aretes.get(e, 0) + 1
        if any(v != 2 for v in aretes.values()):
            ok = False
        signe = []
        for clef in mesFaces:
            pl = X[faces[clef]]
            nn = np.cross(pl[1] - pl[0], pl[2] - pl[0])
            s = 1 if nn @ (pl.mean(axis=0) - c['seed']) > 0 else -1
            signe.append((clef, s))
        cellFaces.append(signe)

    # chaque face interne vue par 2 cellules exactement
    compte = {}
    for cf in cellFaces:
        for clef, _ in cf:
            compte[clef] = compte.get(clef, 0) + 1
    for clef, n in compte.items():
        if clef[0] == 'g' and n != 2:
            ok = False
        if clef[0] == 'd' and n != 1:
            ok = False

    return X, faces, cellFaces, dict(ok=ok, diamCluster=diam, tol=tol,
                                     nBrut=len(brut), nFusion=len(X))


def volume_cellules(cells):
    """Volume exact de chaque cellule (somme des tetraedres germe-face)."""
    V = []
    for c in cells:
        tot = 0.0
        for clef, (ks, nrm) in c['faces'].items():
            poly = c['V'][ks]
            o = _ordonne(poly, nrm)
            poly = poly[o]
            c0 = poly.mean(axis=0)
            for q in range(len(poly)):
                a = poly[q] - c['seed']
                b = poly[(q + 1) % len(poly)] - c['seed']
                d = c0 - c['seed']
                tot += abs(np.dot(np.cross(a, b), d)) / 6.0
        V.append(tot)
    return np.array(V)


def arete_min(X, faces):
    m = float('inf')
    for pl in faces.values():
        for q in range(len(pl)):
            m = min(m, float(np.linalg.norm(X[pl[q]] - X[pl[(q + 1) % len(pl)]])))
    return m


# ===========================================================================
# 5. PHASES — tirage par deficit de volume
# ===========================================================================
def phases_par_grain(Vcell, noms, cibles, rng, mode="deficit"):
    n = len(Vcell)
    ph = np.zeros(n, dtype=int)
    if mode == "multinomial":
        ph = rng.choice(len(noms), size=n, p=np.array(cibles) / sum(cibles))
    else:
        Vtot = float(Vcell.sum())
        acc = np.zeros(len(noms))
        for i in rng.permutation(n):
            k = int(np.argmax(np.array(cibles) * Vtot - acc))
            ph[i] = k
            acc[k] += Vcell[i]
    reel = np.array([Vcell[ph == k].sum() for k in range(len(noms))])
    return ph, reel / reel.sum()


# ===========================================================================
# 6. GMSH — geometrie partagee + maillage Delaunay non structure
# ===========================================================================
def maille(X, faces, cellFaces, phGrain, noms, h, out, algo3d=1,
           optNetgen=True, planTol=None, verbeux=False):
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if verbeux else 0)
    gmsh.model.add("gbm_cyl")
    geo = gmsh.model.geo

    planTol = planTol if planTol is not None else 1e-3 * h
    used = set()
    for pl in faces.values():
        used.update(int(k) for k in pl)
    ptTag = {}
    for k in sorted(used):
        ptTag[k] = geo.addPoint(float(X[k][0]), float(X[k][1]), float(X[k][2]), h)

    lig = {}

    def ligne(a, b):
        key = (min(a, b), max(a, b))
        if key not in lig:
            lig[key] = geo.addLine(ptTag[key[0]], ptTag[key[1]])
        return lig[key] if a < b else -lig[key]

    surfFace, nSplit = {}, 0
    for clef, pl in faces.items():
        P = X[pl]
        c = P.mean(axis=0)
        nn = np.cross(P[1] - P[0], P[2] - P[0])
        nn = nn / max(np.linalg.norm(nn), 1e-30)
        dev = float(np.max(np.abs((P - c) @ nn)))
        if dev <= planTol:
            bd = [ligne(int(pl[q]), int(pl[(q + 1) % len(pl)]))
                  for q in range(len(pl))]
            cl = geo.addCurveLoop(bd)
            surfFace[clef] = [geo.addPlaneSurface([cl])]
        else:                      # face non plane apres fusion : eventail
            nSplit += 1            # partage (meme centre des deux cotes)
            pc = geo.addPoint(float(c[0]), float(c[1]), float(c[2]), h)
            ptTag[('c', clef)] = pc
            ss = []
            radial = {}

            def rad(a):
                if a not in radial:
                    radial[a] = geo.addLine(pc, ptTag[a])
                return radial[a]

            for q in range(len(pl)):
                a, b = int(pl[q]), int(pl[(q + 1) % len(pl)])
                cl = geo.addCurveLoop([rad(a), ligne(a, b), -rad(b)])
                ss.append(geo.addPlaneSurface([cl]))
            surfFace[clef] = ss

    volTag = []
    for cf in cellFaces:
        ss = []
        for clef, s in cf:
            for t in surfFace[clef]:
                ss.append(s * t)
        sl = geo.addSurfaceLoop(ss)
        volTag.append(geo.addVolume([sl]))
    geo.synchronize()

    for k, nom in enumerate(noms):
        vv = [volTag[i] for i in range(len(volTag)) if phGrain[i] == k]
        if vv:
            gmsh.model.addPhysicalGroup(3, vv, k + 1)
            gmsh.model.setPhysicalName(3, k + 1, nom)

    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 5)          # Delaunay 2D
    gmsh.option.setNumber("Mesh.Algorithm3D", algo3d)   # 1 = Delaunay 3D
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1 if optNetgen else 0)
    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.option.setNumber("Mesh.SaveAll", 0)

    # tetraedres par groupe physique (pour l'ecriture de secours + le rapport)
    nt, coords, _ = gmsh.model.mesh.getNodes()
    xyz = np.array(coords).reshape(-1, 3)
    idx = {int(t): i for i, t in enumerate(nt)}
    tets, tphys = [], []
    for k, nom in enumerate(noms):
        try:
            ents = gmsh.model.getEntitiesForPhysicalGroup(3, k + 1)
        except Exception:
            continue
        for e in ents:
            ty, et, en = gmsh.model.mesh.getElements(3, int(e))
            for t, tags, conn in zip(ty, et, en):
                if t != 4:
                    continue
                c = np.array(conn, dtype=np.int64).reshape(-1, 4)
                for row in c:
                    tets.append([idx[int(z)] for z in row])
                    tphys.append(k + 1)
    tets = np.array(tets, dtype=np.int64)
    tphys = np.array(tphys, dtype=np.int64)

    gmsh.write(out)
    gmsh.finalize()
    return xyz, tets, tphys, nSplit


def ecrit_msh22(path, xyz, tets, tphys, noms):
    """Ecriture de secours : MSH 2.2 ASCII, tetraedres SEULS, aucun noeud
    orphelin (le lecteur de rockim refuse un noeud jamais reference)."""
    used = np.unique(tets)
    ren = {int(u): i + 1 for i, u in enumerate(used)}
    with open(path, "w", encoding="ascii", newline="\n") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write(f"$PhysicalNames\n{len(noms)}\n")
        for k, nom in enumerate(noms):
            f.write(f'3 {k + 1} "{nom}"\n')
        f.write("$EndPhysicalNames\n")
        f.write(f"$Nodes\n{len(used)}\n")
        for i, u in enumerate(used):
            p = xyz[int(u)]
            f.write(f"{i + 1} {p[0]:.17g} {p[1]:.17g} {p[2]:.17g}\n")
        f.write("$EndNodes\n")
        f.write(f"$Elements\n{len(tets)}\n")
        for e in range(len(tets)):
            a, b, c, d = (ren[int(z)] for z in tets[e])
            f.write(f"{e + 1} 4 2 {tphys[e]} {tphys[e]} {a} {b} {c} {d}\n")
        f.write("$EndElements\n")


def verifie_msh(path):
    """Relit le fichier comme le fait rockim : version, tags physiques
    declares, aucun noeud orphelin."""
    nomsPhys, nodes, tets, tags = {}, set(), 0, set()
    ver, refs = None, set()
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        it = iter(f.read().split("\n"))
        for line in it:
            if line.startswith("$MeshFormat"):
                ver = float(next(it).split()[0])
            elif line.startswith("$PhysicalNames"):
                n = int(next(it))
                for _ in range(n):
                    t = next(it).split(None, 2)
                    if int(t[0]) == 3:
                        nomsPhys[int(t[1])] = t[2].strip().strip('"')
            elif line.startswith("$Nodes"):
                n = int(next(it))
                for _ in range(n):
                    nodes.add(int(next(it).split()[0]))
            elif line.startswith("$Elements"):
                n = int(next(it))
                for _ in range(n):
                    t = next(it).split()
                    ty, ntg = int(t[1]), int(t[2])
                    if ty == 4:
                        tets += 1
                        tags.add(int(t[3]))
                        refs.update(int(z) for z in t[3 + ntg:3 + ntg + 4])
    pb = []
    if ver is None or not (2.0 <= ver < 3.0):
        pb.append(f"version MSH {ver} (2.2 attendu)")
    if tets == 0:
        pb.append("aucun tetraedre")
    manque = tags - set(nomsPhys)
    if manque:
        pb.append(f"tags physiques non declares : {sorted(manque)}")
    orph = nodes - refs
    if orph:
        pb.append(f"{len(orph)} noeud(s) orphelin(s) (rockim les refuse)")
    return dict(ok=not pb, pb=pb, ver=ver, tets=tets, noeuds=len(nodes),
                noms=nomsPhys, tags=sorted(tags))


# ===========================================================================
# 7. RAPPORT DE QUALITE
# ===========================================================================
def qualite_tets(xyz, tets):
    p = xyz[tets]
    v0, v1, v2, v3 = p[:, 0], p[:, 1], p[:, 2], p[:, 3]
    V = np.abs(np.einsum('ij,ij->i', np.cross(v1 - v0, v2 - v0), v3 - v0)) / 6.0
    A = np.zeros(len(tets))
    for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
        a, b, c = p[:, f[0]], p[:, f[1]], p[:, f[2]]
        A += 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    lc = np.cbrt(np.maximum(V, 1e-300))
    rap = (6.0 * V / np.maximum(A, 1e-300)) / lc
    ar = []
    for (i, j) in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
        ar.append(np.linalg.norm(p[:, i] - p[:, j], axis=1))
    ar = np.array(ar)
    return V, lc, rap, ar.min()


# ===========================================================================
# 8. MAIN
# ===========================================================================
def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    ap = argparse.ArgumentParser(
        description="Microstructure GBM cylindrique (Voronoi decoupe + Gmsh) "
                    "pour rockim fem3d mesh = file.")
    ap.add_argument("--diametre", type=float, default=16e-3, help="D [m]")
    ap.add_argument("--hauteur", type=float, default=32e-3, help="H [m]")
    ap.add_argument("--grain", type=float, default=2.5e-3,
                    help="diametre equivalent moyen de grain [m]")
    ap.add_argument("--maille", type=float, default=0.6e-3, help="taille cible [m]")
    ap.add_argument("--germe", type=int, default=20260906, help="germe aleatoire")
    ap.add_argument("--sortie", default=None, help="fichier .msh")
    ap.add_argument("--phases", nargs="+",
                    default=["feldspath", "quartz", "biotite"])
    ap.add_argument("--fractions", nargs="+", type=float,
                    default=[0.62, 0.31, 0.07])
    ap.add_argument("--tirage", choices=["deficit", "multinomial"],
                    default="deficit")
    ap.add_argument("--minDist", type=float, default=0.80,
                    help="distance mini de Poisson / (V/N)^(1/3)")
    ap.add_argument("--margeBord", type=float, default=0.0,
                    help="marge germe-paroi, en fraction de la distance mini")
    ap.add_argument("--facettes", type=int, default=0,
                    help="M plans de paroi (0 = auto par --tolParoi)")
    ap.add_argument("--tolParoi", type=float, default=0.02,
                    help="ecart radial max tolere, en fraction de --maille")
    ap.add_argument("--fusion", type=float, default=0.05,
                    help="tolerance de fusion des sommets, fraction de --maille")
    ap.add_argument("--voisins", type=int, default=48, help="K initial")
    ap.add_argument("--miroirs", action="store_true",
                    help="ajoute les germes images (grains de bord a face plane)")
    ap.add_argument("--algo3d", type=int, default=1,
                    help="Mesh.Algorithm3D (1 = Delaunay, 10 = HXT)")
    ap.add_argument("--sansNetgen", action="store_true")
    ap.add_argument("--sansMaillage", action="store_true",
                    help="geometrie et controles seulement")
    ap.add_argument("--rapport", default=None, help="rapport JSON optionnel")
    ap.add_argument("--verbeux", action="store_true")
    a = ap.parse_args()

    R, H, h = 0.5 * a.diametre, a.hauteur, a.maille
    if len(a.phases) != len(a.fractions):
        sys.exit("--phases et --fractions doivent avoir la meme longueur")
    fr = np.array(a.fractions, dtype=float)
    fr = fr / fr.sum()
    out = a.sortie or os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "meshes",
        f"gbm_cyl_D{a.diametre*1e3:g}_H{a.hauteur*1e3:g}"
        f"_g{a.grain*1e3:g}_h{a.maille*1e3:g}_s{a.germe}.msh")
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    rng = np.random.default_rng(a.germe)
    t0 = time.time()

    print("=" * 74)
    print(f"GBM cylindrique  D = {a.diametre*1e3:g} mm, H = {a.hauteur*1e3:g} mm, "
          f"grain {a.grain*1e3:g} mm, maille {a.maille*1e3:g} mm, germe {a.germe}")
    print("=" * 74)

    # --- §1 germes -------------------------------------------------------
    P, N, rmin = germes_poisson(R, H, a.grain, rng, a.minDist, a.margeBord)
    dmin = float(cKDTree(P).query(P, k=2)[0][:, 1].min())
    print(f"[1] germes : {len(P)} (vise {N}), disque de Poisson r = "
          f"{rmin*1e3:.3f} mm, plus proche couple reel {dmin*1e3:.3f} mm "
          f"({dmin/a.grain:.2f} x grain)")

    # --- §4 paroi --------------------------------------------------------
    if a.facettes > 0:
        M = a.facettes
        aa = R * math.sqrt(math.pi / (M * math.tan(math.pi / M)))
        eInt, eExt = R - aa, aa / math.cos(math.pi / M) - R
    else:
        M, aa, eInt, eExt = facettes_cylindre(R, H, a.tolParoi * h)
    print(f"[2] paroi : {M} plans, rayon egalisateur d'aire a = {aa*1e3:.6f} mm ; "
          f"ecart radial {eInt*1e6:+.2f} um (fleche) / {eExt*1e6:+.2f} um "
          f"(sommet), soit {max(abs(eInt),abs(eExt))/h*100:.2f} % de la maille")
    Adom, bdom, nomdom = plans_domaine(R, H, M, aa)
    Pall = P if not a.miroirs else np.vstack([P, germes_miroirs(P, R, H)])
    print(f"    germes images : {'oui' if a.miroirs else 'non'} "
          f"({len(Pall) - len(P)} plans supplementaires)")

    # --- §2/§3 cellules --------------------------------------------------
    cells, Kmax = cellules(P, Pall, Adom, bdom, nomdom, a.voisins, 1.0 / R)
    Vc = volume_cellules(cells)
    Vcyl = M * aa * aa * math.tan(math.pi / M) * H      # = pi R^2 H a 1e-15 pres
    Vnom = math.pi * R * R * H
    print(f"[3] cellules : {len(cells)}, certificat des voisins satisfait "
          f"(K max = {Kmax}) ; aucune cellule ouverte par construction")
    print(f"    somme des volumes {Vc.sum()*1e9:.6f} mm3 / cylindre facette "
          f"{Vcyl*1e9:.6f} mm3  -> ecart {abs(Vc.sum()/Vcyl-1)*100:.3e} % "
          f"(cylindre exact {Vnom*1e9:.6f} mm3)")
    if abs(Vc.sum() / Vcyl - 1.0) > 1e-6:
        print("    *** ATTENTION : trou ou recouvrement, la partition n'est "
              "pas exacte")

    # --- §6 complexe partage --------------------------------------------
    tol = a.fusion * h
    while True:
        X, faces, cellFaces, dg = complexe(cells, tol, R)
        if dg["ok"] or tol <= 0.0:
            break
        print(f"    fusion {tol*1e6:.2f} um casse la fermeture -> "
              f"{tol*0.5*1e6:.2f} um")
        tol = 0.0 if tol < 1e-4 * h else 0.5 * tol
    X0, faces0, _, _ = complexe(cells, 0.0, R)
    aMin0, aMin = arete_min(X0, faces0), arete_min(X, faces)
    nInt = sum(1 for k in faces if k[0] in ("g", "m"))
    print(f"[4] complexe : {len(X)} sommets ({dg['nBrut']} bruts fusionnes a "
          f"{dg['tol']*1e6:.2f} um), {len(faces)} faces dont {nInt} internes")
    print(f"    fermeture (2 faces par arete, 2 cellules par face interne) : "
          f"{'OK' if dg['ok'] else 'ECHEC'}")
    print(f"    plus petite arete {aMin0*1e6:.2f} um -> {aMin*1e6:.2f} um apres "
          f"regularisation ({aMin/h*100:.1f} % de la maille)")

    # --- §5 phases -------------------------------------------------------
    ph, reel = phases_par_grain(Vc, a.phases, fr, rng, a.tirage)
    print("[5] phases (fraction VOLUMIQUE realisee / visee) :")
    for k, nom in enumerate(a.phases):
        print(f"    {nom:<12s} {100*reel[k]:6.2f} % / {100*fr[k]:6.2f} %  "
              f"({int((ph == k).sum())} grains)")
    deq = np.cbrt(6.0 * Vc / math.pi)
    print(f"    grains : d_eq median {np.median(deq)*1e3:.3f} mm, "
          f"[min {deq.min()*1e3:.3f} ; p10 {np.percentile(deq,10)*1e3:.3f} ; "
          f"p90 {np.percentile(deq,90)*1e3:.3f} ; max {deq.max()*1e3:.3f}] mm")

    if a.sansMaillage:
        print(f"[6] --sansMaillage : arret avant Gmsh ({time.time()-t0:.1f} s)")
        return

    # --- §6 maillage -----------------------------------------------------
    print(f"[6] maillage Gmsh (Delaunay non structure, h = {h*1e3:g} mm) ...")
    xyz, tets, tphys, nSplit = maille(X, faces, cellFaces, ph, a.phases, h,
                                      out, a.algo3d, not a.sansNetgen,
                                      verbeux=a.verbeux)
    if nSplit:
        print(f"    {nSplit} face(s) non planes apres fusion -> triangulees "
              "(partagees, donc toujours conformes)")
    V, lc, rap, arMin = qualite_tets(xyz, tets)
    print(f"    {len(tets)} tetraedres, {len(np.unique(tets))} noeuds")
    print(f"    somme des volumes de tets {V.sum()*1e9:.6f} mm3 / cellules "
          f"{Vc.sum()*1e9:.6f} mm3 -> ecart {abs(V.sum()/Vc.sum()-1)*100:.3e} %")
    print(f"[7] qualite : lc median {np.median(lc)*1e3:.4f} mm "
          f"(vise {h*1e3:g}) ; diametre inscrit / lc : median "
          f"{np.median(rap):.3f}, min {rap.min():.3f}  "
          f"(tetraedre REGULIER = 0,833 ; sous 0,1 la bande de fissuration "
          f"est fausse)")
    print(f"    p1 = {np.percentile(rap,1):.3f}, "
          f"tets sous 0,1 : {int((rap < 0.1).sum())} "
          f"({100*(rap < 0.1).mean():.3f} %)")
    print(f"    plus petite arete de tet {arMin*1e6:.2f} um "
          f"({arMin/h*100:.1f} % de la maille)")
    vph = np.array([V[tphys == k + 1].sum() for k in range(len(a.phases))])
    print("    fractions volumiques MAILLEES / visees :")
    for k, nom in enumerate(a.phases):
        print(f"      {nom:<12s} {100*vph[k]/vph.sum():6.2f} % / "
              f"{100*fr[k]:6.2f} %")

    # --- verification du fichier ----------------------------------------
    v = verifie_msh(out)
    if not v["ok"]:
        print(f"    export Gmsh non conforme ({'; '.join(v['pb'])}) "
              "-> reecriture directe")
        ecrit_msh22(out, xyz, tets, tphys, a.phases)
        v = verifie_msh(out)
    print(f"[8] fichier : {out}")
    print(f"    MSH {v['ver']}, {v['noeuds']} noeuds, {v['tets']} tets, "
          f"groupes {v['noms']} -> {'CONFORME' if v['ok'] else 'PROBLEME: ' + '; '.join(v['pb'])}")
    print(f"    deck : mode = fem3d / mesh = file / meshFile = {os.path.basename(out)} ; "
          "phases = " + " ".join(a.phases))
    print(f"[9] duree totale {time.time()-t0:.1f} s")

    if a.rapport:
        json.dump(dict(
            D=a.diametre, H=a.hauteur, grain=a.grain, maille=h, germe=a.germe,
            grains=len(cells), facettesParoi=M, ecartParoi=[eInt, eExt],
            volCellules=float(Vc.sum()), volCylindre=Vcyl, volTets=float(V.sum()),
            fractionsVisees=fr.tolist(), fractionsRealisees=reel.tolist(),
            fractionsMaillees=(vph / vph.sum()).tolist(),
            dEq=dict(min=float(deq.min()), med=float(np.median(deq)),
                     max=float(deq.max())),
            tets=int(len(tets)), lcMed=float(np.median(lc)),
            rapMed=float(np.median(rap)), rapMin=float(rap.min()),
            areteMinGeom=float(aMin), areteMinTet=float(arMin),
            fermeture=bool(dg["ok"]), fichier=out, msh=v["ok"]),
            open(a.rapport, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
        print(f"    rapport JSON : {a.rapport}")


if __name__ == "__main__":
    main()
