#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# make_masses_decks.py — construit les trois decks de LECTURE (2 us) qui
# mesurent la correction de masse proposee par T3 (docs/MAILLAGE_serie §4).
# Aucun deck de configs/ n'est modifie : les trois sortent dans ce dossier.
#
#   python tests_f2/campagne13/T3/make_masses_decks.py
#
# Temoin  : configs/yang2026_bench_s25_v3P.cfg, maillage de la serie (rock25,
#           train fige) — le solveur doit imprimer les masses de l'outil
#           tools/mesh_quality.py --masses (piston 1,058054, bit 1,288297).
# _rho    : deux phases acier de plus (steelPiston 8702,8 ; steelBit 8714,5)
#           + groupPhase.piston / .bit — le solveur doit imprimer 1,173 et
#           1,4302 kg (bit + insert + circlip = 1,509 = Yang 2026, p. 11).
# _rhoE   : idem avec E scale du meme facteur (221,7 et 222,0 GPa) : c = sqrt(E/rho)
#           inchange, donc le pas de temps de l'acier inchange.
# Critere FALSIFIANT : si `groupPhase.<corps> = <phase nouvelle>` etait lu et
# INERTE (piege n. 1 du depot), les trois decks imprimeraient les MEMES masses.
# ---------------------------------------------------------------------------
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
ICI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import mesh_quality as mq                                        # noqa: E402

TEMOIN = "configs/yang2026_bench_s25_v3P.cfg"
MESH = "meshes/impact_yang_train1_rock25_hxt.msh"
M_YANG_PISTON, M_YANG_BIT = 1.173, 1.509     # Yang et al. 2026, p. 11


def masses(rel):
    nodes, tets, phys, names = mq.read_msh(os.path.join(ROOT, rel))
    h, em, cen, V = mq.inscribed(nodes, tets)
    out = {}
    for tag in sorted(set(phys.tolist())):
        nm = names.get(tag, str(tag))
        out[nm] = float(V[phys == tag].sum()) * mq.RHO_DEFAUT[nm]
    return out


def bloc(rho_p, rho_b, e_p, e_b):
    return "\n".join([
        "",
        "# --- T3 (13/09) : correction de masse du train par DENSITE PAR CORPS ---",
        "# Les masses publiees (piston 1,173 kg, bit 1,509 kg) ne sont pas celles du",
        "# dessin V1 du generateur (facettisation -6,0 / -4,0 %, longueurs -4,0 / -11,1 %,",
        "# docs/MAILLAGE_serie_2026-09-13.md §4). Deux phases acier de plus portent la",
        "# densite corrigee ; la geometrie et le maillage ne changent pas.",
        "phases = rock steel steelPiston steelBit carbide",
        "phase.steelPiston.fraction = 0.0140",
        "phase.steelPiston.rho = %.1f" % rho_p,
        "phase.steelPiston.E   = %.4g" % e_p,
        "phase.steelPiston.nu  = 0.29",
        "phase.steelPiston.ft  = 1e12",
        "phase.steelPiston.cohesion = 1e12",
        "phase.steelBit.fraction = 0.0170",
        "phase.steelBit.rho = %.1f" % rho_b,
        "phase.steelBit.E   = %.4g" % e_b,
        "phase.steelBit.nu  = 0.29",
        "phase.steelBit.ft  = 1e12",
        "phase.steelBit.cohesion = 1e12",
        "groupPhase.piston = steelPiston",
        "groupPhase.bit    = steelBit",
        "",
    ])


def ecrire(nom, txt, extra):
    # une seule ligne `phases = ...` dans le deck final : celle du bloc T3
    lignes = []
    for l in txt.splitlines():
        s = l.strip()
        if s.startswith("meshFile"):
            lignes.append("meshFile = %s" % MESH)
            continue
        if extra and s.startswith("phases ="):
            lignes.append("# [T3] phases remplace plus bas : " + s)
            continue
        # une seule affectation par corps : les deux lignes d'origine sont
        # commentees, pour ne pas dependre de la regle « premiere ou derniere
        # occurrence gagne » du lecteur de deck.
        if extra and (s.startswith("groupPhase.piston") or s.startswith("groupPhase.bit ")
                      or s.startswith("groupPhase.bit=")):
            lignes.append("# [T3] remplace plus bas : " + s)
            continue
        lignes.append(l)
    p = os.path.join(ICI, nom)
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lignes) + (extra or "") + "\n")
    return p


def main():
    m = masses(MESH)
    fp = M_YANG_PISTON / m["piston"]
    fb = (M_YANG_BIT - m["insert"] - m["circlip"]) / m["bit"]
    rho_p, rho_b = 7850.0 * fp, 7850.0 * fb
    print("maillage %s" % MESH)
    print("  masses maillees  : piston %.6f  bit %.6f  insert %.6f  circlip %.6f  plate %.6f  rock %.6f"
          % (m["piston"], m["bit"], m["insert"], m["circlip"], m["plate"], m["rock"]))
    print("  facteurs         : piston x %.6f -> rho %.1f ; bit x %.6f -> rho %.1f" % (fp, rho_p, fb, rho_b))
    print("  ATTENDU au solveur : piston %.4f  bit %.4f  bit+insert+circlip %.6f (Yang %.3f)"
          % (m["piston"] * fp, m["bit"] * fb, m["bit"] * fb + m["insert"] + m["circlip"], M_YANG_BIT))
    with open(os.path.join(ROOT, TEMOIN), encoding="utf-8") as f:
        txt = f.read()
    p0 = ecrire("t3_masses_temoin.cfg", txt, None)
    p1 = ecrire("t3_masses_rho.cfg", txt, bloc(rho_p, rho_b, 200e9, 200e9))
    p2 = ecrire("t3_masses_rhoE.cfg", txt, bloc(rho_p, rho_b, 200e9 * fp, 200e9 * fb))
    for p in (p0, p1, p2):
        print("  ecrit : %s" % os.path.relpath(p, ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
