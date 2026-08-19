#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_yang_bench.py — maillage EQUIVALENT au banc d impact de Yang, Xiang,
# Naderi, Wang, Aising, Ugarte & Latham (IJRMMS 191, 2025, 106125), section 3.1.
#
#   python tunnel_edz/tools/make_yang_bench.py meshes/yang_bench.msh [seed]
#
# LEUR GRADATION, mot pour mot (leur section 3.1) :
#   « Within the rock model, a hemispherical volume with a fine mesh and a
#     diameter of 25 mm is assigned a mesh size of 1.0 mm. The mesh size
#     gradually increases from 25 mm diameter to 50 mm diameter, where it
#     becomes 2.0 mm. Beyond the 50 mm diameter, the mesh size gradually
#     coarsens outwards, reaching 10 mm at the boundary. »
#
# TROIS niveaux, donc, avec deux transitions DOUCES. Le mode bench1g du
# generateur commun n en a que DEUX (une zone fine, un lointain), ce qui
# produisait un saut de 2,3 a 13,4 mm en une marche a r = 25 mm — mesure du
# 2026-08-19. Or leurs fissures radiales font 20 a 25 mm (leur Table 3) :
# elles tombaient donc dans des elements de 13 mm et n auraient jamais pu etre
# resolues. D ou ce generateur dedie.
#
# EQUIVALENCE ASSUMEE, a ecrire partout ou ce maillage sert :
#   * leur eprouvette est un cylindre Ø250 x 150 ; ici un BLOC de meme cote.
#     Leur propre texte dit que la forme est indifferente tant que le bord est
#     a plus de trois longueurs de fissure du point d impact.
#   * leur insert hemispherique R 8,51 mm est soude a un taillant Ø30 x 265 ;
#     ici une SPHERE libre de meme rayon. La masse du taillant est rendue par
#     la densite de la phase insert (voir la config), donc la quantite de
#     mouvement et l energie sont les leurs — mais PAS la dynamique d onde
#     dans le taillant (aller-retour 106 us, comparable a l evenement).
#   * l insert est maille a 1,5 mm et non 0,7 : il reste elastique et ne sert
#     qu a pousser. A 0,7 mm il pesait 31 000 tetraedres sur 101 000, soit
#     30 % du maillage la ou il ne se passe rien.
# ---------------------------------------------------------------------------
import sys

import gmsh

# --- leur geometrie et leur gradation ---------------------------------------
W = 0.160          # cote du bloc [m] : leur bord est a >3 longueurs de
                   # fissure (24,5 mm) du point d impact ; 80 mm de demi-cote
                   # tient ce critere pour un quart des elements
HR = 0.100         # hauteur de roche [m]
RI = 0.00851       # rayon de l insert [m]       (leur R = 8,51 mm)
GAP = 1.0e-4       # jeu initial insert / surface [m]
H1 = 0.0010        # taille fine                 (leur 1,0 mm)
R1 = 0.0125        # ... jusqu a r = 12,5 mm     (leur Ø25)
H2 = 0.0020        # taille intermediaire        (leur 2,0 mm)
R2 = 0.0250        # ... a r = 25 mm             (leur Ø50)
H3 = 0.0080        # taille au bord (leur 10 mm, resserre avec la boite)
PENTE = 0.10       # dh/dr : leurs DEUX transitions ont la meme pente,
                   # (2-1)/(25-12,5) = (10-2)/(125-25) = 0,08 — la gradation
                   # se ramene donc a UNE droite. On la resserre a 0,10.
HINS = 0.0015      # insert (voir l en-tete)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "meshes/yang_bench.msh"
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.RandomSeed", seed / 100.0)

    # roche : sommet a z = 0, descend en -z
    rock = gmsh.model.occ.addBox(-0.5 * W, -0.5 * W, -HR, W, W, HR)
    # insert : sphere posee au-dessus du point d impact, separee de GAP
    ins = gmsh.model.occ.addSphere(0.0, 0.0, RI + GAP, RI)
    gmsh.model.occ.synchronize()

    gmsh.model.addPhysicalGroup(3, [rock], name="rock")
    gmsh.model.addPhysicalGroup(3, [ins], name="insert")

    # ---- champ de taille : UNE formule, aucune entite ---------------------
    # La gradation de Yang se ramene a une droite : h(r) = min(H3, max(H1,
    # pente*r)). Le champ Distance de la version precedente s appuyait sur un
    # POINT isole que gmsh a ignore — resultat, 4,5 M de tetraedres a taille
    # uniforme (mesure du 2026-08-19). MathEval n a besoin d aucune entite.
    fm = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(
        fm, "F", "min(%g, max(%g, %g*sqrt(x^2+y^2+z^2)))" % (H3, H1, PENTE))
    fr = gmsh.model.mesh.field.add("Restrict")
    gmsh.model.mesh.field.setNumber(fr, "InField", fm)
    gmsh.model.mesh.field.setNumbers(fr, "VolumesList", [rock])

    fi = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(fi, "F", "%g" % HINS)
    fri = gmsh.model.mesh.field.add("Restrict")
    gmsh.model.mesh.field.setNumber(fri, "InField", fi)
    gmsh.model.mesh.field.setNumbers(fri, "VolumesList", [ins])

    fmin = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", [fr, fri])
    gmsh.model.mesh.field.setAsBackgroundMesh(fmin)

    for k in ("Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature",
              "Mesh.MeshSizeExtendFromBoundary"):
        gmsh.option.setNumber(k, 0)
    gmsh.option.setNumber("Mesh.Algorithm", 5)     # Delaunay isotrope (peau)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)   # Delaunay (volume)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    gmsh.model.mesh.generate(3)
    gmsh.write(out)

    # ---- controle imprime : la gradation est-elle celle demandee ? ---------
    import numpy as np
    nt, ntg, npar = gmsh.model.mesh.getElements(3)
    tags = np.concatenate([np.array(x).reshape(-1, 4) for x in npar])
    ncoord = gmsh.model.mesh.getNodes()[1].reshape(-1, 3)
    nid = gmsh.model.mesh.getNodes()[0]
    idx = {int(t): k for k, t in enumerate(nid)}
    P = np.array([ncoord[idx[int(t)]] for t in tags.ravel()]).reshape(-1, 4, 3)
    cen = P.mean(axis=1)
    e = np.zeros(len(P))
    for i in range(4):
        for j in range(i + 1, 4):
            e += np.linalg.norm(P[:, i] - P[:, j], axis=1)
    e /= 6.0
    r = np.hypot(cen[:, 0], cen[:, 1])
    haut = cen[:, 2] > 0.0
    print("\n[maillage] %d tetraedres (roche %d, insert %d)"
          % (len(P), int((~haut).sum()), int(haut.sum())))
    for lo, hi in [(0, 12.5), (12.5, 25), (25, 50), (50, 125)]:
        m = (~haut) & (r >= lo * 1e-3) & (r < hi * 1e-3)
        if m.sum():
            print("  roche r %3.0f-%3.0f mm : %7d tet, arete med %5.2f mm"
                  % (lo, hi, m.sum(), np.median(e[m]) * 1e3))
    if haut.any():
        print("  insert            : %7d tet, arete med %5.2f mm"
              % (haut.sum(), np.median(e[haut]) * 1e3))
    gmsh.finalize()
    print("  ecrit :", out)


if __name__ == "__main__":
    main()
