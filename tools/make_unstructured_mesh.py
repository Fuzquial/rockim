#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_unstructured_mesh.py — genere un maillage simplexe NON STRUCTURE
# uniforme (le maillage "a la Yan et al. 2023") pour rockim `mesh = file`.
#
#   python3 tools/make_unstructured_mesh.py box3d W D H h out.msh [seed]
#   python3 tools/make_unstructured_mesh.py box2d W H   h out.msh [seed]
#   python3 tools/make_unstructured_mesh.py bench1 W D H R gap h hIns out.msh [seed]
#   python3 tools/make_unstructured_mesh.py tunnelhs W H hFine rFine hFar out.msh [seed]
#
# bench1 (V1) : DEUX corps — bloc W x D x H (volume physique "rock") et
# INSERT spherique de rayon R centre au-dessus (volume physique "insert"),
# separes de `gap`, mailles a h (roche) et hIns (insert). C'est le premier
# jalon de la trajectoire Solidity : l'impact insert-roche en maille.
#
# Sortie : Gmsh MSH 2.2 ASCII (tets type 4 en 3D, triangles type 2 en 2D).
# Gmsh est LE mailleur de la pratique FDEM (Akantu, OpenFDEM, litterature) ;
# pip install gmsh. Algorithmes : Delaunay + optimisation Netgen (qualite —
# le min du diametre inscrit pilote la CFL de rockim). Le seed perturbe le
# point de depart du frontal 2D pour varier la realisation.
# ---------------------------------------------------------------------------
import math
import sys
import gmsh


# ---------------------------------------------------------------------------
# Profil en FER A CHEVAL du tunnel de Hutou Beishan — Wang et al., Front. Earth
# Sci. 12:1517816 (2024), fig. 6a. Le dessin ne cote pas la construction, il
# cote des LONGUEURS : 11 m de portee, 8,85 m de hauteur totale, 5,55 m pour la
# voute, 1,6 m de piedroit, 1,73 m de conge d'about, 4,25 / 1,225 m au radier,
# plus un rayon "R555" et deux angles (17,64 deg, 111,45 deg).
#
# RECONSTRUCTION RETENUE (a verifier si le plan d'origine devient disponible) :
#   - voute      : arc R = 5,55 m centre sur l'axe a y = 8,85 - 5,55 = 3,30 m,
#                  de naissance a naissance (la portee vaut donc 11,1 m, cote
#                  11 m au dessin) ;
#   - piedroits  : verticaux de y = 1,73 a y = 3,30 m -> 1,57 m (cote 1,6) ;
#   - conges     : arc R3 = demi-portee - 4,25 = 1,30 m (cote 1,225), TANGENT
#                  au piedroit en (5,55 ; 1,73) ;
#   - radier     : arc R4 resolu pour etre TANGENT aux conges et horizontal sur
#                  l'axe (radier en cuvette, point bas y = 0).
# Toutes les longueurs cotees sont restituees a quelques cm pres et le contour
# est G1 partout : aucun angle vif, donc aucune singularite de contrainte
# parasite a l'about du radier — ce qui compte davantage, pour une EDZ, que la
# construction exacte a trois centres. Les deux ANGLES cotes ne sont pas
# utilises : ils suggerent un dessin a trois centres un peu plus fin.
# ---------------------------------------------------------------------------
TUNNEL_HS = dict(halfSpan=5.55, height=8.85, rCrown=5.55, yWall=1.73, xInv=4.25)


def tunnel_hs_profile(P=TUNNEL_HS):
    """Geometrie derivee du profil : rayons, centres, points de tangence."""
    halfSpan, height = P["halfSpan"], P["height"]
    ySpring = height - P["rCrown"]                    # naissance de la voute
    r3 = halfSpan - P["xInv"]                         # conge d'about
    k = P["yWall"] - r3
    if k <= 0.0:
        raise SystemExit("profil tunnel : yWall doit depasser le rayon de conge")
    # tangence interne conge / radier : |C3 - C4| = R4 - R3, avec
    # C3 = (xInv, yWall) et C4 = (0, R4)  ->  R4 ferme en une ligne
    r4 = (P["xInv"] ** 2 - k * k) / (2.0 * k) + P["yWall"]
    c3 = (P["xInv"], P["yWall"])
    c4 = (0.0, r4)
    d = math.hypot(c3[0] - c4[0], c3[1] - c4[1])
    ux, uy = (c3[0] - c4[0]) / d, (c3[1] - c4[1]) / d
    tang = (c4[0] + r4 * ux, c4[1] + r4 * uy)         # radier <-> conge
    return dict(ySpring=ySpring, r3=r3, r4=r4, c3=c3, c4=c4, tang=tang, **P)


def build_tunnel_hs(cx, cy0, P=TUNNEL_HS):
    """Boucle fermee du contour, invert au niveau cy0, axe en x = cx.

    Renvoie (tag de la boucle, geometrie derivee).
    """
    g = tunnel_hs_profile(P)
    occ = gmsh.model.occ
    hs, ht, ys, yw = g["halfSpan"], g["height"], g["ySpring"], g["yWall"]
    tx, ty = g["tang"]

    def pt(x, y):
        return occ.addPoint(cx + x, cy0 + y, 0.0)

    # points du contour (droite -> haut -> gauche) et centres des arcs
    pA = pt(0.0, 0.0)                       # point bas du radier (sur l'axe)
    pBr, pBl = pt(tx, ty), pt(-tx, ty)      # radier <-> conge
    pCr, pCl = pt(hs, yw), pt(-hs, yw)      # conge <-> piedroit
    pDr, pDl = pt(hs, ys), pt(-hs, ys)      # piedroit <-> voute (naissance)
    pE = pt(0.0, ht)                        # clef de voute
    o4 = pt(0.0, g["r4"])
    o3r, o3l = pt(g["c3"][0], g["c3"][1]), pt(-g["c3"][0], g["c3"][1])
    o1 = pt(0.0, ys)
    arcs = [occ.addCircleArc(pA, o4, pBr),
            occ.addCircleArc(pBr, o3r, pCr),
            occ.addLine(pCr, pDr),
            occ.addCircleArc(pDr, o1, pE),   # demi-voute droite (90 deg)
            occ.addCircleArc(pE, o1, pDl),   # demi-voute gauche
            occ.addLine(pDl, pCl),
            occ.addCircleArc(pCl, o3l, pBl),
            occ.addCircleArc(pBl, o4, pA)]
    return occ.addCurveLoop(arcs), g


def main():
    kind = sys.argv[1]
    if kind == "box3d":
        W, D, H, h = map(float, sys.argv[2:6])
        out = sys.argv[6]
        seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    elif kind == "box2d":
        W, H, h = map(float, sys.argv[2:5])
        out = sys.argv[5]
        seed = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    elif kind == "bench1":
        W, D, H, R, gap, h, hIns = map(float, sys.argv[2:9])
        out = sys.argv[9]
        seed = int(sys.argv[10]) if len(sys.argv) > 10 else 1
    elif kind == "bench1g":
        # bench1 GRADUE : meme geometrie, mais la roche est fine dans un
        # CYLINDRE sous le point d'impact et grossiere ailleurs. Impose par le
        # diagnostic du 2026-08-17 : resoudre le contact de Hertz demande
        # ~15 elements sur le rayon de contact (dx < 0.5 mm), ce qu'un
        # raffinement uniforme du bloc 120^3 mettrait hors d'atteinte.
        #   hFin  : taille dans la zone fine
        #   rFin  : rayon du cylindre fin (autour de l'axe d'impact)
        #   dFin  : profondeur du cylindre fin sous la surface
        W, D, H, R, gap, h, hIns, hFin, rFin, dFin = map(float, sys.argv[2:12])
        out = sys.argv[12]
        seed = int(sys.argv[13]) if len(sys.argv) > 13 else 1
    elif kind == "tunnel":
        # plaque W x H percee d'un trou circulaire R au centre (2D) —
        # cavite pressurisee par confineFaces = bore
        W, H, R, h = map(float, sys.argv[2:6])
        out = sys.argv[6]
        seed = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    elif kind == "tunnelhs":
        # massif W x H perce du tunnel en FER A CHEVAL (Wang et al. 2024),
        # maillage GRADUE : hFine jusqu'a rFine de la PAROI, hFar au loin.
        W, H, hFine, rFine, hFar = map(float, sys.argv[2:7])
        out = sys.argv[7]
        seed = int(sys.argv[8]) if len(sys.argv) > 8 else 1
        # ALGORITHME 2D — mesure du 2026-08-17, a ne pas remettre a 6 :
        # `Mesh.Algorithm = 6` (frontal-Delaunay) pave en triangles QUASI
        # EQUILATERAUX des que la taille est localement constante. Mesure sur
        # la zone fine de tunnel_hs.msh : mediane du plus petit angle 60,0 deg,
        # 95,7 % des triangles au-dessus de 50 deg, histogramme d'orientation
        # periodique a 60 deg avec un rapport pic/creux de 19 a 30. Les
        # fissures n'ont alors que trois directions disponibles et suivent des
        # droites sur plusieurs metres — condition d'invalidite pour un calcul
        # de facies. L'algorithme 5 (Delaunay) ramene pic/creux a 1,16.
        # Le lissage n'y est pour rien (sans lissage : 32,2) et bruiter le
        # champ de taille AGGRAVE le cas de l'algo 6 (R6 0,167 -> 0,552).
        algo2d = int(sys.argv[9]) if len(sys.argv) > 9 else 5
        h = hFine
    else:
        raise SystemExit("usage: make_unstructured_mesh.py box3d W D H h out.msh [seed]\n"
                         "       make_unstructured_mesh.py box2d W H h out.msh [seed]\n"
                         "       make_unstructured_mesh.py bench1 W D H R gap h hIns out.msh [seed]\n"
                         "       make_unstructured_mesh.py tunnel W H R h out.msh [seed]\n"
                         "       make_unstructured_mesh.py tunnelhs W H hFine rFine hFar out.msh [seed]")

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("rockim")
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.RandomSeed", seed)
    gmsh.option.setNumber("Mesh.Algorithm", 6)        # frontal-Delaunay 2D
    gmsh.option.setNumber("Mesh.Optimize", 1)
    if kind == "box3d":
        gmsh.model.occ.addBox(0, 0, 0, W, D, H)
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)  # Delaunay 3D
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        gmsh.model.mesh.generate(3)
    elif kind == "bench1":
        rock = gmsh.model.occ.addBox(0, 0, 0, W, D, H)
        ins = gmsh.model.occ.addSphere(0.5 * W, 0.5 * D, H + gap + R, R)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(3, [rock], name="rock")
        gmsh.model.addPhysicalGroup(3, [ins], name="insert")
        # taille locale : fine sur l'insert, h partout ailleurs
        gmsh.option.setNumber("Mesh.MeshSizeMin", min(h, hIns))
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        insPts = gmsh.model.getBoundary([(3, ins)], recursive=True)
        for dim, tag in insPts:
            if dim == 0:
                gmsh.model.mesh.setSize([(0, tag)], hIns)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        gmsh.model.mesh.generate(3)
    elif kind == "bench1g":
        rock = gmsh.model.occ.addBox(0, 0, 0, W, D, H)
        ins = gmsh.model.occ.addSphere(0.5 * W, 0.5 * D, H + gap + R, R)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(3, [rock], name="rock")
        gmsh.model.addPhysicalGroup(3, [ins], name="insert")
        # Champ de taille : hFin dans le cylindre (axe d'impact, rayon rFin,
        # profondeur dFin sous la surface), h au loin, transition lissee sur
        # une longueur egale a rFin pour eviter un saut de taille brutal (qui
        # degrade la qualite et donc le dt).
        cx, cy, zTop = 0.5 * W, 0.5 * D, H
        f = gmsh.model.mesh.field.add("MathEval")
        # r = distance a l'axe ; p = profondeur sous la surface
        expr = ("%.9g + (%.9g - %.9g) * max(0, min(1, (sqrt((x-%.9g)^2 + "
                "(y-%.9g)^2)/%.9g - 1)))" % (hFin, h, hFin, cx, cy, rFin))
        exprZ = ("%.9g + (%.9g - %.9g) * max(0, min(1, ((%.9g - z)/%.9g - 1)))"
                 % (hFin, h, hFin, zTop, dFin))
        # L'insert est ENTIEREMENT au-dessus de z = H : on lui impose hIns dans
        # le champ lui-meme. Sans cela, MeshSizeFromPoints = 0 rend inerte le
        # setSize() pose sur ses sommets, la sphere tombe dans la zone fine du
        # champ (elle est sur l'axe d'impact) et se fait mailler a hFin —
        # constate le 2026-08-17 : 62 752 tets dans l'insert au lieu de ~1 400,
        # soit 28 % du maillage pour 0,3 % du volume.
        step = "max(0, min(1, (z - %.9g)/%.9g))" % (zTop, 1e-4)
        rockF = "Max(%s, %s)" % (expr, exprZ)
        gmsh.model.mesh.field.setString(
            f, "F", "(1 - %s)*(%s) + (%s)*%.9g" % (step, rockF, step, hIns))
        gmsh.model.mesh.field.setAsBackgroundMesh(f)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.MeshSizeMin", min(hFin, hIns))
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        # l'insert garde sa propre finesse (il porte le contact)
        for dim, tag in gmsh.model.getBoundary([(3, ins)], recursive=True):
            if dim == 0:
                gmsh.model.mesh.setSize([(0, tag)], hIns)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        gmsh.model.mesh.generate(3)
    elif kind == "tunnel":
        plate = gmsh.model.occ.addRectangle(0, 0, 0, W, H)
        hole = gmsh.model.occ.addDisk(0.5 * W, 0.5 * H, 0, R, R)
        gmsh.model.occ.cut([(2, plate)], [(2, hole)])
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(2)
    elif kind == "tunnelhs":
        cx = 0.5 * W
        cy0 = 0.5 * H - 0.5 * TUNNEL_HS["height"]     # section centree
        plate = gmsh.model.occ.addRectangle(0, 0, 0, W, H)
        loop, g = build_tunnel_hs(cx, cy0)
        cav = gmsh.model.occ.addPlaneSurface([loop])
        gmsh.model.occ.cut([(2, plate)], [(2, cav)])
        gmsh.model.occ.synchronize()
        # Champ de taille : distance a la PAROI (les courbes du contour de la
        # cavite, retrouvees par boite englobante apres la coupe, qui renumerote
        # tout), puis seuil hFine -> hFar. La transition s'etale sur 0,8 rFine :
        # une marche de taille degrade la qualite, donc le dt (lecon bench1g).
        wall = []
        for dim, tag in gmsh.model.getEntities(1):
            x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(dim, tag)
            if (x0 > cx - g["halfSpan"] - 1e-6 and x1 < cx + g["halfSpan"] + 1e-6
                    and y0 > cy0 - 1e-6 and y1 < cy0 + g["height"] + 1e-6):
                wall.append(tag)
        if not wall:
            raise SystemExit("tunnelhs : contour de la cavite introuvable")
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
        gmsh.option.setNumber("Mesh.Algorithm", algo2d)   # cf. commentaire plus haut
        gmsh.model.mesh.generate(2)
        print(f"[tunnel] Mesh.Algorithm = {algo2d}"
              f"{' (Delaunay, isotrope)' if algo2d == 5 else ''}"
              f"{' (frontal-Delaunay : PAVAGE QUASI STRUCTURE, cf. commentaire)' if algo2d == 6 else ''}")
        print(f"[tunnel] portee {2*g['halfSpan']:.2f} m, hauteur {g['height']:.2f} m, "
              f"R voute {g['rCrown']:.2f} m (centre y = {g['ySpring']:.2f}), "
              f"conge R {g['r3']:.2f} m, radier R {g['r4']:.2f} m "
              f"(tangence x = {g['tang'][0]:.2f}, releve {g['tang'][1]:.2f} m)")
        print(f"[tunnel] section centree en ({cx:.1f}, {0.5*H:.1f}) m, "
              f"radier a y = {cy0:.3f} m, {len(wall)} courbes de paroi")
    else:
        gmsh.model.occ.addRectangle(0, 0, 0, W, H)
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(out)
    # bilan qualite : diametre inscrit min/med (ce qui pilote le dt rockim)
    import numpy as np
    if kind in ("box3d", "bench1"):
        _, _, conn = gmsh.model.mesh.getElements(3)
        tets = np.array(conn[0], dtype=int).reshape(-1, 4)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        idx = {t: i for i, t in enumerate(tags)}
        P = np.array(xyz).reshape(-1, 3)[[idx[t] for t in tets.flatten()]]
        P = P.reshape(-1, 4, 3)
        V = np.abs(np.einsum("ij,ij->i",
                             np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]),
                             P[:, 3] - P[:, 0])) / 6.0
        F = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
        A = sum(0.5 * np.linalg.norm(
                np.cross(P[:, f[1]] - P[:, f[0]], P[:, f[2]] - P[:, f[0]]),
                axis=1) for f in F)
        hin = 6.0 * V / A
        print(f"[mesh] {len(tets)} tets, h inscrit min/med/max = "
              f"{hin.min()*1e3:.3f}/{np.median(hin)*1e3:.3f}/{hin.max()*1e3:.3f} mm")
    elif kind == "tunnelhs":
        # meme bilan en 2D : c'est le diametre inscrit MIN qui fixe le dt, et
        # le pire angle qui dit si le mailleur a lache des slivers.
        _, _, conn = gmsh.model.mesh.getElements(2)
        tri = np.array(conn[0], dtype=int).reshape(-1, 3)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        idx = {t: i for i, t in enumerate(tags)}
        P = np.array(xyz).reshape(-1, 3)[[idx[t] for t in tri.flatten()], :2]
        P = P.reshape(-1, 3, 2)
        e = np.stack([np.linalg.norm(P[:, 2] - P[:, 1], axis=1),
                      np.linalg.norm(P[:, 0] - P[:, 2], axis=1),
                      np.linalg.norm(P[:, 1] - P[:, 0], axis=1)], axis=1)
        u, v = P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]   # np.cross 2D : retire en numpy 2
        Ar = 0.5 * np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0])
        hin = 4.0 * Ar / e.sum(axis=1)             # diametre du cercle inscrit
        ang = np.degrees(np.arccos(np.clip(
            (e[:, [1, 2, 0]] ** 2 + e[:, [2, 0, 1]] ** 2 - e ** 2)
            / (2 * e[:, [1, 2, 0]] * e[:, [2, 0, 1]]), -1, 1)))
        worst = ang.min()
        print(f"[mesh] {len(tri)} triangles, h inscrit min/med/max = "
              f"{hin.min():.4f}/{np.median(hin):.4f}/{hin.max():.4f} m, "
              f"pire angle {worst:.1f} deg, "
              f"{100.0*np.mean(ang.min(axis=1) < 20.0):.2f} % sous 20 deg")
        # Regle de dimensionnement maison (2026-08-17) : 2 dx < l_cz = E Gf/ft^2.
        # Materiau de l'article : E = 10 GPa, Gf = 20 J/m^2, ft = 0,6 MPa.
        lcz = 10e9 * 20.0 / (0.6e6 ** 2)
        print(f"[mesh] l_cz (E Gf/ft^2, materiau Wang 2024) = {lcz:.3f} m -> "
              f"dx admissible <= {0.5*lcz:.3f} m ; hFine = {hFine:.3f} m "
              f"[{'OK' if hFine <= 0.5*lcz else 'TROP GROSSIER'}]")
    gmsh.finalize()


if __name__ == "__main__":
    main()
