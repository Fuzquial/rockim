#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_tunnel_bedded_mesh.py — BRIQUE 3 de Lisjak (these Toronto 2013, ch. 5,
# fig. 5.8b) : le maillage du tunnel en fer a cheval avec des ARETES ALIGNEES
# SUR LE LITAGE dans la zone raffinee.
#
#   python tunnel_schisto/make_tunnel_bedded_mesh.py W H hFine rFine hFar \
#          dipDeg t rBed out.msh [seed] [algo]
#   ex :   python tunnel_schisto/make_tunnel_bedded_mesh.py 100 100 0.12 18 1.8 \
#          45 0.35 24 meshes/tunnel_hs_bed45.msh 1
#
# POURQUOI. La loi cohesive directionnelle (bedding*) donne une resistance
# faible aux joints presque paralleles au litage — mais si le maillage n a
# aucune arete CONTINUE le long du litage, il n existe aucun CHEMIN faible :
# la delamination doit zigzaguer sur des aretes mal orientees et la
# resistance au cisaillement est gonflee. Verbatim de la source : « the mesh
# topology must combine a random triangulation for the intra-layer material
# (i.e., matrix) together with crack elements preferably aligned along the
# bedding planes ». Regle de discretisation constante aux deux echelles de
# Lisjak : h ~ t/3 (labo : t = 1,0 mm, h = 0,3 mm ; tunnel de 3 m : t = 10 cm,
# h = 3 cm). Par homothetie pour le tunnel de 11 m : t ~ 35 cm, h ~ 12 cm.
#
# COMMENT (v6, 2026-09-02). La famille de cordes paralleles (pendage dipDeg,
# espacement t) du disque r < rBed est DECOUPEE EN PYTHON : chaque corde est
# echantillonnee (pas h/4), un point est utilisable s il est dans la roche ET
# a plus de 0,8 h de la paroi (gmsh.model.isInside + getClosestPoint), les
# runs utilisables sont bornes par bissection (1e-4 m) puis raccourcis de
# 0,5 h, et les troncons sont PLONGES dans la face roche par mesh.embed.
# Gmsh y pose des aretes conformes ; ailleurs, triangulation aleatoire.
# Voir l HISTORIQUE v1-v6 dans le corps du script : le fragment OCC ne fait
# pas ce travail, et un retrait fixe laisse des slivers le long de la paroi.
#
# QUALITE. Le pas de temps de rockim est fixe par le PLUS PETIT element
# (m_min) ET, en mesh = file, par pj = 4E/hmin GLOBAL : un seul sliver
# penalise tout le modele. Cible : diametre inscrit minimal >= celui du
# maillage isotrope de reference (73 mm a h = 0,2). v6 : 76 mm.
#
# CE QUE LE LECTEUR DE ROCKIM EN FAIT. rockim lit les triangles de type 2 et
# JETTE les tags (src/FdemSolver.cpp, `in >> tag`) : aucun grain ni phase
# parasite. Verifie le 2026-09-01.
#
# CONTROLES imprimes : nombre de triangles, h median dans la zone litee, et
# — le chiffre qui compte — la CONTINUITE : longueur des aretes exactement
# alignees (< 1 deg) rapportee a la longueur des cordes dans la roche. Ne
# PAS juger sur une tolerance angulaire large : la triangulation contrainte
# appauvrit son voisinage dans la direction de la ligne et masque tout.
# ---------------------------------------------------------------------------
import os
import sys

import gmsh
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "rockim", "rockim_p1", "tools"))
from make_unstructured_mesh import build_tunnel_hs, TUNNEL_HS  # noqa: E402


def main():
    if len(sys.argv) < 10:
        raise SystemExit(__doc__ or "usage: voir l en-tete du script")
    W, H, hFine, rFine, hFar, dip, t, rBed = map(float, sys.argv[1:9])
    out = sys.argv[9]
    seed = int(sys.argv[10]) if len(sys.argv) > 10 else 1
    algo = int(sys.argv[11]) if len(sys.argv) > 11 else 5
    if t <= 0 or rBed <= 0:
        raise SystemExit("t et rBed doivent etre > 0")
    if hFine > t / 2.5:
        print(f"[bedded] WARNING: hFine = {hFine} > t/2.5 = {t/2.5:.3f} : Lisjak "
              f"tient h ~ t/3, sinon les couches ne sont pas resolues")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("bedded")
    gmsh.option.setNumber("Mesh.RandomSeed", seed)
    occ = gmsh.model.occ

    cx = 0.5 * W
    cy0 = 0.5 * H - 0.5 * TUNNEL_HS["height"]
    cyc = 0.5 * H                                 # centre du disque lite
    plate = occ.addRectangle(0, 0, 0, W, H)
    loop, g = build_tunnel_hs(cx, cy0)
    cav = occ.addPlaneSurface([loop])
    rock, _ = occ.cut([(2, plate)], [(2, cav)])
    occ.synchronize()

    # ---- la famille de cordes, PRE-DECOUPEES par la cavite ------------------
    # HISTORIQUE (2026-09-01/02), pour ne pas y retomber :
    #   v1  fragment par un disque + cordes : les extremites sur le cercle
    #       ressortaient sans face (tolerance OCC) -> supprimees -> isotrope.
    #   v2  fragment par les cordes seules : le General Fuse DECOUPE a la
    #       paroi mais n ATTACHE pas ; mesh.embed des troncons regle ca.
    #   v3  MAIS le fragment cree un SOMMET sur la paroi a chaque extremite de
    #       troncon ; quand il tombe a quelques mm d un noeud regulier de la
    #       paroi, Gmsh fabrique un sliver. Mesure : 17 triangles sous 40 mm,
    #       d_inscrit min 20 mm au lieu de 73 -> et comme pj = 4E/hmin GLOBAL
    #       en mesh = file, la penalite de TOUS les joints etait gonflee x3,6.
    #       Cout : dt divise par 2 pour tout le modele.
    # v4 (ici) : plus de fragment du tout. Chaque corde est decoupee EN
    # PYTHON par bissection sur gmsh.model.isInside (precision 1e-4 m), et
    # chaque troncon est RACCOURCI de 0,5 hFine a ses extremites de paroi :
    # il s arrete DANS la roche, aucun sommet n est impose sur le contour, la
    # paroi garde son maillage regulier, et le dernier segment vers la paroi
    # est une arete ordinaire. Puis mesh.embed direct.
    faces = [tg for _, tg in gmsh.model.getEntities(2)]
    if len(faces) != 1:
        raise SystemExit(f"bedded : {len(faces)} faces apres cut, 1 attendue")
    rockTag = faces[0]
    inside = lambda p: bool(gmsh.model.isInside(2, rockTag, [p[0], p[1], 0.0]))
    d = np.array([np.cos(np.radians(dip)), np.sin(np.radians(dip))])   # litage
    n = np.array([-d[1], d[0]])                                        # normale
    kmax = int(np.floor(rBed / t))
    ds = 0.25 * hFine                     # pas d echantillonnage le long
    # v5 : le RETRAIT depend de l angle theta entre la corde et la tangente
    # de la paroi. Mesure v4 (retrait fixe 0,5 h) : 109 slivers, TOUS a
    # 3-10 cm de la paroi — dans le vide en coin entre la fin de corde et le
    # contour, ou Gmsh doit loger des elements plus petits que h. Pour qu un
    # element entier y tienne, il faut un ecart PERPENDICULAIRE >= h, donc
    # un retrait le long de la corde l = h / sin(theta), borne a [0,75 h ; 4 h].
    # Les courbes de paroi servent a lire theta (point le plus proche +
    # derivee parametrique).
    wallTags = []
    for dim_, tag_ in gmsh.model.getEntities(1):
        x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(dim_, tag_)
        if (x0 > cx - g["halfSpan"] - 1e-6 and x1 < cx + g["halfSpan"] + 1e-6
                and y0 > cy0 - 1e-6 and y1 < cy0 + g["height"] + 1e-6):
            wallTags.append(tag_)
    if not wallTags:
        raise SystemExit("bedded : contour de la cavite introuvable")

    def trim_at(p):
        """retrait le long de la corde au point de paroi p"""
        best = (1e30, None, None)
        for tg in wallTags:
            cp, u = gmsh.model.getClosestPoint(1, tg, [p[0], p[1], 0.0])
            dist = np.hypot(cp[0] - p[0], cp[1] - p[1])
            if dist < best[0]:
                best = (dist, tg, u[0])
        der = gmsh.model.getDerivative(1, best[1], [best[2]])
        tan = np.array([der[0], der[1]])
        tan /= max(np.linalg.norm(tan), 1e-30)
        sinth = abs(d[0] * tan[1] - d[1] * tan[0])       # |d x tangente|
        return float(np.clip(hFine / max(sinth, 0.25), 0.75 * hFine, 4.0 * hFine))

    def near_wall(p):
        """distance du point p a la paroi (grande valeur loin du tunnel)"""
        if np.hypot(p[0] - cx, p[1] - cyc) > 9.0:
            return 1e9                     # la paroi est a r < 6 m du centre
        best = 1e30
        for tg in wallTags:
            cp, _ = gmsh.model.getClosestPoint(1, tg, [p[0], p[1], 0.0])
            best = min(best, np.hypot(cp[0] - p[0], cp[1] - p[1]))
        return best

    usable = lambda p: inside(p) and near_wall(p) >= 0.8 * hFine

    keep, chordLenRock = [], []
    nChords = 0
    for k in range(-kmax, kmax + 1):
        off = k * t
        half = np.sqrt(max(rBed * rBed - off * off, 0.0))
        if half < 0.5 * t:
            continue
        nChords += 1
        c0 = np.array([cx, cyc]) + off * n
        s = np.arange(-half, half + ds, ds)
        pts = c0[None, :] + s[:, None] * d[None, :]
        # v6 : un point de corde est « utilisable » s il est dans la roche ET
        # a plus de 0,8 h de la paroi. Mesure v5 : les 6 slivers restants
        # n etaient PAS a des extremites (retrait plus large : maillage
        # identique) mais la ou une corde LONGE le radier a 6-7 cm, presque
        # parallele — entre la corde et la paroi il ne tient qu un element
        # d un tiers de h. Exclure la bande de 0,8 h garantit un element
        # entier entre les deux, quelle que soit l incidence.
        flags = np.array([inside(p) and near_wall(p) >= 0.8 * hFine
                          for p in pts])
        # runs de points DANS la roche -> troncons ; bissection aux bords
        i = 0
        while i < len(s):
            if not flags[i]:
                i += 1
                continue
            j = i
            while j + 1 < len(s) and flags[j + 1]:
                j += 1
            sa, sb = s[i], s[j]
            # bissection sur le predicat COMPLET (roche ET a >= 0,8 h de la
            # paroi), puis retrait fixe de 0,5 h : la bande garantit deja
            # l ecart perpendiculaire, le retrait ne sert qu a eloigner le
            # sommet terminal de la frontiere de bande.
            if i > 0:                      # bord gauche : bissection
                lo, hi = s[i - 1], s[i]
                for _ in range(20):
                    mid = 0.5 * (lo + hi)
                    if usable(c0 + mid * d): hi = mid
                    else: lo = mid
                sa = hi + 0.5 * hFine
            if j + 1 < len(s):             # bord droit : bissection
                lo, hi = s[j], s[j + 1]
                for _ in range(20):
                    mid = 0.5 * (lo + hi)
                    if usable(c0 + mid * d): lo = mid
                    else: hi = mid
                sb = lo - 0.5 * hFine
            if sb - sa > 2.0 * hFine:      # troncon utile
                p1 = c0 + sa * d
                p2 = c0 + sb * d
                a = occ.addPoint(p1[0], p1[1], 0.0)
                b = occ.addPoint(p2[0], p2[1], 0.0)
                keep.append(occ.addLine(a, b))
                chordLenRock.append(float(sb - sa))
            i = j + 1
    occ.synchronize()
    print(f"[bedded] {nChords} cordes -> {len(keep)} troncons dans la roche, "
          f"retrait a la paroi h/sin(theta) borne [0,75 h ; 4 h], plonges (mesh.embed)")
    if not keep:
        raise SystemExit("bedded : aucun troncon de litage dans la roche")
    gmsh.model.mesh.embed(1, keep, 2, rockTag)

    # ---- paroi de la cavite (pour le champ de distance) --------------------
    wall = []
    for dim, tag in gmsh.model.getEntities(1):
        x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(dim, tag)
        if (x0 > cx - g["halfSpan"] - 1e-6 and x1 < cx + g["halfSpan"] + 1e-6
                and y0 > cy0 - 1e-6 and y1 < cy0 + g["height"] + 1e-6):
            wall.append(tag)
    if not wall:
        raise SystemExit("bedded : contour de la cavite introuvable")

    # ---- champ de taille : identique a make_unstructured_mesh tunnelhs ----
    fd = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(fd, "CurvesList", wall)
    gmsh.model.mesh.field.setNumber(fd, "Sampling", 400)
    fth = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(fth, "InField", fd)
    gmsh.model.mesh.field.setNumber(fth, "SizeMin", hFine)
    gmsh.model.mesh.field.setNumber(fth, "SizeMax", hFar)
    gmsh.model.mesh.field.setNumber(fth, "DistMin", rFine)
    gmsh.model.mesh.field.setNumber(fth, "DistMax", rFine + 0.8 * rFine)
    gmsh.model.mesh.field.setAsBackgroundMesh(fth)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeMin", hFine)
    gmsh.option.setNumber("Mesh.MeshSizeMax", hFar)
    gmsh.option.setNumber("Mesh.Algorithm", algo)
    gmsh.logger.start()
    gmsh.model.mesh.generate(2)
    # Controle que le plongement a PRIS : des noeuds doivent exister sur les
    # troncons plonges. Zero = Gmsh les a ignores en silence (vu le
    # 2026-09-01 dans une version anterieure de ce script).
    nEmb = sum(len(gmsh.model.mesh.getNodes(1, tg)[0]) for tg in keep)
    warns = [l for l in gmsh.logger.get()
             if any(k in l for k in ("mbed", "arning", "rror"))]
    gmsh.logger.stop()
    print(f"[bedded] noeuds poses sur les troncons plonges : {nEmb}"
          + ("  <<< ZERO : plongement ignore" if nEmb == 0 else ""))
    for w in warns[:8]:
        print("[bedded]   gmsh :", w[:150])
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(out)

    # ---- controles ---------------------------------------------------------
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    idx = {tg: i for i, tg in enumerate(tags)}
    P = np.array(xyz).reshape(-1, 3)[:, :2]
    _, _, conn = gmsh.model.mesh.getElements(2)
    T = np.array([[idx[v] for v in tri] for tri in np.array(conn[0]).reshape(-1, 3)])
    gmsh.finalize()

    E = np.vstack([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
    E.sort(axis=1)
    E = np.unique(E, axis=0)
    mid = 0.5 * (P[E[:, 0]] + P[E[:, 1]])
    v = P[E[:, 1]] - P[E[:, 0]]
    L = np.linalg.norm(v, axis=1)
    inBed = np.hypot(mid[:, 0] - cx, mid[:, 1] - cyc) < rBed
    cosg = np.abs(v[:, 0] * d[0] + v[:, 1] * d[1]) / np.maximum(L, 1e-30)
    gam = np.degrees(np.arccos(np.clip(cosg, 0.0, 1.0)))
    # METRIQUE. Une premiere version comptait les aretes a moins de 10 deg du
    # litage ; elle ne bougeait pas (10,2 % lite contre 12,1 % isotrope)
    # alors que 1614 aretes etaient EXACTEMENT alignees : la triangulation
    # contrainte par une ligne APPAUVRIT son voisinage dans la direction de
    # la ligne (aretes a +-60 deg), ce qui compense. On compte donc les
    # aretes a moins de 1 deg — celles qui SONT les cordes — et on rapporte
    # leur longueur a la longueur ideale des cordes dans la roche : c est la
    # CONTINUITE, meme definition que pour weakPlanes dans le solveur.
    exact = inBed & (gam < 1.0)
    lAligned = L[exact].sum()
    lIdeal = sum(gmsh_len for gmsh_len in chordLenRock)
    cont = lAligned / max(lIdeal, 1e-30)
    print(f"[bedded] {len(T)} triangles, {len(P)} noeuds -> {out}")
    print(f"[bedded] zone litee r < {rBed} m : {inBed.sum()} aretes, h median "
          f"{np.median(L[inBed])*1000:.0f} mm, litage {dip} deg, t = {t} m "
          f"({2*kmax+1} cordes)")
    print(f"[bedded] aretes EXACTEMENT alignees (< 1 deg) : {exact.sum()} "
          f"({100*np.mean(gam[inBed] < 1.0):.1f} % ; isotrope ~0,7 %), "
          f"longueur {lAligned:.0f} m pour {lIdeal:.0f} m de cordes dans la "
          f"roche -> CONTINUITE = {cont:.3f}")
    if cont < 0.9:
        print("[bedded] WARNING: continuite < 0,9 — une partie des cordes n a "
              "pas ete plongee ; verifier le journal gmsh ci-dessus")


if __name__ == "__main__":
    main()
