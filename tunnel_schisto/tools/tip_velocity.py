#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# tip_velocity.py — S7, la vitesse de pointe de fissure, a comparer a c_R.
#
#   python tunnel_schisto/tools/tip_velocity.py out_X [--E 10e9 --nu 0.26 --rho 2650]
#
# POURQUOI. He & Hutchinson reservent leur critere aux fissures LENTES ;
# Chalivendra & Rosakis 2008 montrent qu une fissure rapide est plus facilement
# piegee par une interface. Dans un code explicite, le relachement d excavation
# en 0,08 s peut lancer des fissures a une fraction notable de c_R.
#
# MESURE. Pour chaque joint rompu, la vitesse d arrivee = distance au joint
# rompu VOISIN (partageant un noeud) le plus recent qui l a precede, divisee
# par l ecart de leurs tBreak. On rapporte la distribution (p50, p90, p99, max)
# et la fraction des ruptures qui depassent 0,3 c_R. Le champ tBreak est
# ecrit par defaut : aucun run a refaire.
# ---------------------------------------------------------------------------
import argparse
import glob
import re
from collections import defaultdict
import numpy as np


def vtu(f):
    raw = open(f, "rb").read().decode("utf8", "ignore")

    def arr(name, dtype=float):
        m = re.search(r'Name="%s"[^>]*>(.*?)</DataArray>' % name, raw, re.S)
        return np.fromstring(m.group(1).strip(), sep=" ", dtype=dtype) if m else None
    m = re.search(r'<Points>.*?<DataArray[^>]*>(.*?)</DataArray>', raw, re.S)
    P = np.fromstring(m.group(1).strip(), sep=" ").reshape(-1, 3)[:, :2]
    return P, arr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--E", type=float, default=10e9)
    ap.add_argument("--nu", type=float, default=0.26)
    ap.add_argument("--rho", type=float, default=2650.0)
    a = ap.parse_args()
    fr = sorted(glob.glob(f"{a.run}/fdem_joints_0*.vtu"))
    P0, a0 = vtu(fr[0]); _, ar = vtu(fr[-1])
    conn = a0("connectivity", int).reshape(-1, 2); tb = ar("tBreak")
    G = a.E / (2 * (1 + a.nu)); cS = np.sqrt(G / a.rho)
    cR = cS * (0.862 + 1.14 * a.nu) / (1 + a.nu)          # Viktorov
    seg = P0[conn]; c = seg.mean(1)
    br = np.where(tb >= 0)[0]
    node_j = defaultdict(list)
    for j in br:
        node_j[conn[j, 0]].append(j); node_j[conn[j, 1]].append(j)
    vel = []
    for j in br:
        best = None
        for nd in (conn[j, 0], conn[j, 1]):
            for k in node_j[nd]:
                if k == j or tb[k] >= tb[j]:
                    continue
                if best is None or tb[k] > tb[best]:
                    best = k
        if best is not None and tb[j] - tb[best] > 0:
            vel.append(np.linalg.norm(c[j] - c[best]) / (tb[j] - tb[best]))
    vel = np.array(vel)
    print(f"--- {a.run} : {len(br)} rompus, {len(vel)} vitesses d arrivee ---")
    print(f"  c_S = {cS:.0f} m/s, c_R = {cR:.0f} m/s (Viktorov)")
    for q in (50, 90, 99):
        print(f"  p{q} = {np.percentile(vel, q):8.1f} m/s  ({np.percentile(vel, q)/cR:5.2f} c_R)")
    print(f"  max = {vel.max():8.1f} m/s  ({vel.max()/cR:5.2f} c_R)")
    print(f"  fraction des ruptures a > 0,3 c_R : {100*np.mean(vel > 0.3*cR):.1f} %")
    print("  (note : une vitesse 'infinie' entre deux joints rompus au meme pas est exclue ; un")
    print("   front qui traverse un element en un pas est borne par h/dt, non par c_R)")


if __name__ == "__main__":
    main()
