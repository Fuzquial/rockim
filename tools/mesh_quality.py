#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# mesh_quality.py — les SLIVERS d'un maillage .msh v2.2 : diametre inscrit
# h = 6 V / (somme des aires des faces) par tetra (la mesure que rockim
# appelle « h inscrit » et qui commande le pas de temps), par corps physique.
#   python tools/mesh_quality.py meshes/a.msh [meshes/b.msh ...] [--worst 5]
# Imprime : N, h min, quantiles 0,1 % / 1 %, nombre de tetras sous 0,3 / 0,5 mm,
# et les pires elements (position, h, arete moyenne, corps).
# ---------------------------------------------------------------------------
import sys
import numpy as np


def read_msh(path):
    names, nodes, tets, phys = {}, None, [], []
    with open(path, errors="replace") as f:
        lines = f.read().split(chr(10))
    i = 0
    while i < len(lines):
        L = lines[i].strip()
        if L == "$PhysicalNames":
            n = int(lines[i + 1])
            for k in range(n):
                d, tag, nm = lines[i + 2 + k].split(maxsplit=2)
                names[int(tag)] = nm.strip().strip('"')
            i += n + 3
        elif L == "$Nodes":
            n = int(lines[i + 1])
            arr = np.array(" ".join(lines[i + 2:i + 2 + n]).split(), float).reshape(n, 4)
            nodes = np.zeros((int(arr[:, 0].max()) + 1, 3))
            nodes[arr[:, 0].astype(int)] = arr[:, 1:]
            i += n + 3
        elif L == "$Elements":
            n = int(lines[i + 1])
            for k in range(n):
                p = lines[i + 2 + k].split()
                if p[1] == "4":
                    nt = int(p[2])
                    phys.append(int(p[3]))
                    tets.append([int(x) for x in p[3 + nt:7 + nt]])
            i += n + 3
        else:
            i += 1
    return nodes, np.array(tets), np.array(phys), names


def inscribed(nodes, tets):
    P = nodes[tets]                                   # (n, 4, 3)
    a, b, c, d = P[:, 0], P[:, 1], P[:, 2], P[:, 3]
    V = np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))) / 6.0
    A = (np.linalg.norm(np.cross(b - a, c - a), axis=1) + np.linalg.norm(np.cross(b - a, d - a), axis=1)
         + np.linalg.norm(np.cross(c - a, d - a), axis=1) + np.linalg.norm(np.cross(c - b, d - b), axis=1)) / 2.0
    h = 6.0 * V / A
    E = np.stack([np.linalg.norm(b - a, axis=1), np.linalg.norm(c - a, axis=1), np.linalg.norm(d - a, axis=1),
                  np.linalg.norm(c - b, axis=1), np.linalg.norm(d - b, axis=1), np.linalg.norm(d - c, axis=1)], 1)
    return h, E.mean(1), P.mean(1), V


def report(path, worst=5):
    nodes, tets, phys, names = read_msh(path)
    h, em, cen, V = inscribed(nodes, tets)
    print("== %s : %d tetras, h min %.4f mm, q0.1%% %.4f mm, q1%% %.4f mm, mediane %.3f mm ; < 0,3 mm : %d ; < 0,5 mm : %d"
          % (path, len(h), h.min() * 1e3, np.quantile(h, 1e-3) * 1e3, np.quantile(h, 1e-2) * 1e3,
             np.median(h) * 1e3, (h < 3e-4).sum(), (h < 5e-4).sum()))
    for tag in sorted(set(phys.tolist())):
        m = phys == tag
        hh = h[m]
        print("   %-8s N %7d  h min %.4f mm  q1%% %.4f  mediane %.3f  arete med. %.3f mm  h/arete min %.3f  < 0,3 mm : %d"
              % (names.get(tag, str(tag)), m.sum(), hh.min() * 1e3, np.quantile(hh, 1e-2) * 1e3,
                 np.median(hh) * 1e3, np.median(em[m]) * 1e3, (hh / em[m]).min(), (hh < 3e-4).sum()))
    idx = np.argsort(h)[:worst]
    for j in idx:
        print("   pire : h %.4f mm  arete moy. %.3f mm  h/arete %.3f  corps %-7s  x %+.4f y %+.4f z %+.4f  (r %.4f)"
              % (h[j] * 1e3, em[j] * 1e3, h[j] / em[j], names.get(int(phys[j]), "?"),
                 cen[j, 0], cen[j, 1], cen[j, 2], np.hypot(cen[j, 0], cen[j, 1])))
    return h


# ---------------------------------------------------------------------------
# Options du 13/09 (T3, campagne de correction) — sortie par defaut inchangee :
#   --ball R    [m] statistiques des tetras dont le CENTROIDE est a r < R de
#               l'origine, par corps : N, mediane de la moyenne des 6 aretes,
#               h min, mediane de h. La « demi-boule R 12,5 mm » du diagnostic
#               du 12/09 (14 722 tetras, 1,372 mm sur impact_yang_s1_pose.msh).
#   --masses    volume et masse par corps (V x rho). rho par corps = phases des
#               decks yang2026_* : rock 2626, bit/piston/plate 7850 (acier),
#               insert/circlip 15250 (carbure) ; surcharge : --rho nom=val.
#               Le solveur imprime la meme masse (somme des masses nodales
#               condensees) dans son resume par corps : c'est le controle.
#   --yang      avec --masses : volumes ANALYTIQUES des corps tels que dessines
#               par make_impact_mesh.py (cylindres, hemisphere + fut, anneau,
#               plaque percee), perte de facettisation, masses publiees par
#               Yang 2026 (piston 1,173 kg, bit 1,509 kg) et densite corrigee
#               rho x m_Yang / m_maille pour les atteindre.
# ---------------------------------------------------------------------------
RHO_DEFAUT = {"rock": 2626.0, "bit": 7850.0, "piston": 7850.0, "plate": 7850.0,
              "insert": 15250.0, "circlip": 15250.0}
M_YANG = {"piston": 1.173, "bit": 1.509}     # ICL p. 11 (Yang et al. 2026)


def analytic_volumes():
    """Volumes des corps de make_impact_mesh.py (memes constantes, en m)."""
    R_ROCK, H_ROCK = 0.125, 0.150
    R_INS, R_SHANK, H_INS = 0.00851, 0.00794, 0.0232
    R_BIT, L_BIT = 0.015, 0.265 - 0.0232
    R_PIS, L_PIS = 0.01325, 0.260
    R_CLIP, H_CLIP = 0.018, 0.003
    PL_X, PL_Y, PL_H, R_HOLE = 0.1199, 0.040, 0.006, 0.0155
    pi = np.pi
    # insert = sphere R_INS fusionnee avec un fut R_SHANK qui part du CENTRE de
    # la sphere : union = sphere + fut - (fut inclus dans l'hemisphere haut)
    h1 = np.sqrt(R_INS ** 2 - R_SHANK ** 2)
    v_int = pi * R_SHANK ** 2 * h1 + pi * (R_INS ** 2 * (R_INS - h1) - (R_INS ** 3 - h1 ** 3) / 3.0)
    v_ins = 4.0 / 3.0 * pi * R_INS ** 3 + pi * R_SHANK ** 2 * (H_INS - R_INS) - v_int
    return {"rock": pi * R_ROCK ** 2 * H_ROCK,
            "insert": v_ins,
            "bit": pi * R_BIT ** 2 * L_BIT,
            "piston": pi * R_PIS ** 2 * L_PIS,
            "circlip": pi * (R_CLIP ** 2 - R_BIT ** 2) * H_CLIP,
            "plate": PL_X * PL_Y * PL_H - pi * R_HOLE ** 2 * PL_H}


def ball_report(path, R):
    nodes, tets, phys, names = read_msh(path)
    h, em, cen, V = inscribed(nodes, tets)
    r = np.linalg.norm(cen, axis=1)
    print("== %s : boule r < %.1f mm (centroides)" % (path, R * 1e3))
    for tag in sorted(set(phys.tolist())):
        m = (phys == tag) & (r < R)
        if not m.any():
            continue
        print("   %-8s N %7d  arete moy. mediane %.4f mm  q10 %.4f  q90 %.4f  h min %.4f mm  h mediane %.4f mm"
              % (names.get(tag, str(tag)), m.sum(), np.median(em[m]) * 1e3, np.quantile(em[m], 0.1) * 1e3,
                 np.quantile(em[m], 0.9) * 1e3, h[m].min() * 1e3, np.median(h[m]) * 1e3))


def masses_report(path, rho, yang=False):
    nodes, tets, phys, names = read_msh(path)
    h, em, cen, V = inscribed(nodes, tets)
    ana = analytic_volumes() if yang else {}
    print("== %s : masses par corps (V x rho)" % path)
    tot = 0.0
    out = {}
    for tag in sorted(set(phys.tolist())):
        nm = names.get(tag, str(tag))
        m = phys == tag
        vol = V[m].sum()
        rh = rho.get(nm, float("nan"))
        mass = vol * rh
        tot += mass
        out[nm] = (m.sum(), vol, mass)
        line = "   %-8s N %7d  V %.6e m3  rho %7.0f  masse %.6f kg" % (nm, m.sum(), vol, rh, mass)
        if yang and nm in ana:
            line += "  | V analytique %.6e  ecart maillage %+.2f %%" % (ana[nm], 100.0 * (vol / ana[nm] - 1.0))
        print(line)
    print("   total    masse %.6f kg" % tot)
    if yang:
        for nm, my in M_YANG.items():
            if nm not in out:
                continue
            n, vol, mass = out[nm]
            print("   Yang %-7s %.3f kg : maille %.6f (%+.2f %%), analytique %.6f (%+.2f %%) ; rho corrigee = %.0f x %.3f / %.6f = %.1f kg/m3"
                  % (nm, my, mass, 100.0 * (mass / my - 1.0), ana[nm] * rho[nm], 100.0 * (ana[nm] * rho[nm] / my - 1.0),
                     rho[nm], my, mass, rho[nm] * my / mass))
        if all(k in out for k in ("bit", "insert", "circlip")):
            mb = out["bit"][2] + out["insert"][2] + out["circlip"][2]
            print("   Yang bit 1.509 kg lu comme bit+insert+circlip : maille %.6f (%+.2f %%) ; rho corrigee du seul corps bit = %.1f kg/m3"
                  % (mb, 100.0 * (mb / 1.509 - 1.0), rho["bit"] * (1.509 - out["insert"][2] - out["circlip"][2]) / out["bit"][2]))
    return out


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    worst = 5
    if "--worst" in sys.argv:
        worst = int(sys.argv[sys.argv.index("--worst") + 1])
        args = [a for a in args if a != str(worst)]
    ball = None
    if "--ball" in sys.argv:
        ball = float(sys.argv[sys.argv.index("--ball") + 1])
        args = [a for a in args if a != sys.argv[sys.argv.index("--ball") + 1]]
    rho = dict(RHO_DEFAUT)
    for i, a in enumerate(sys.argv):
        if a == "--rho":
            k, v = sys.argv[i + 1].split("=")
            rho[k] = float(v)
            args = [x for x in args if x != sys.argv[i + 1]]
    masses = "--masses" in sys.argv
    yang = "--yang" in sys.argv
    for p in args:
        if not masses and ball is None:
            report(p, worst)
        if ball is not None:
            ball_report(p, ball)
        if masses:
            masses_report(p, rho, yang)
