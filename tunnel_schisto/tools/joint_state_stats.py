#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# joint_state_stats.py — S8, le diagnostic de Renshaw-Pollard aux plans
# delamines. Exige un run avec writeJointState = true (rockim_f2j+).
#
#   python tunnel_schisto/tools/joint_state_stats.py out_X --dip 45 [--frame N]
#
# QUESTION. Les plans de litage qui BORNENT l EDZ sont-ils ouverts, fermes et
# glissants, ou fermes et bloques a la fin du run ? Un plan ouvert ou glissant
# ne transmet aucune traction : l arret d une fissure qui y arrive est
# MECANIQUE (Renshaw & Pollard 1995 : pas de traversee sans blocage en
# frottement). Un plan ferme-bloque qui arrete quand meme, c est le rapport
# d energies (He & Hutchinson) qui bloque — et la revue dit que 0,057 le
# garantit.
#
# MESURES.
#  1. Etat de contact (contactState) des joints rompus PARALLELES au litage
#     (< 15 deg), a la FRONTIERE de l EDZ (10 % les plus externes par secteur)
#     contre l interieur.
#  2. Contrainte normale sigN et cisaillement tauS moyens par etat.
#  3. EVENEMENTS DE TRAVERSEE : joints rompus EN TRAVERS (> 60 deg) dont les
#     deux extremites touchent chacune un joint rompu parallele DIFFERENT —
#     c est-a-dire une fissure continue de part et d autre d un plan.
# ---------------------------------------------------------------------------
import argparse
import glob
import re
from collections import defaultdict
import numpy as np

STATE = {0: "ouvert", 1: "ferme-glissant", 2: "ferme-bloque", 3: "cohesif",
         4: "mort (contact general)", 5: "lie (non insere)"}


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
    ap.add_argument("--dip", type=float, required=True)
    ap.add_argument("--frame", type=int, default=-1)
    ap.add_argument("--cx", type=float, default=50.0)
    ap.add_argument("--cy", type=float, default=50.0)
    a = ap.parse_args()
    fr = sorted(glob.glob(f"{a.run}/fdem_joints_0*.vtu"))
    f = fr[a.frame]
    P0, a0 = vtu(fr[0])
    P, ar = vtu(f)
    conn = a0("connectivity", int).reshape(-1, 2)
    st = ar("contactState")
    if st is None:
        raise SystemExit("pas de champ contactState : relancer avec writeJointState = true")
    sig = ar("sigN"); tau = ar("tauS"); tb = ar("tBreak"); br = tb >= 0
    seg = P0[conn]; v = seg[:, 1] - seg[:, 0]
    ang = np.degrees(np.arctan2(v[:, 1], v[:, 0])) % 180
    dev = np.abs(((ang - a.dip) + 90) % 180 - 90)
    para = dev < 15; cross = dev > 60
    c = seg.mean(1); rad = np.hypot(c[:, 0] - a.cx, c[:, 1] - a.cy)
    az = np.degrees(np.arctan2(c[:, 1] - a.cy, c[:, 0] - a.cx)) % 360
    bnd = np.zeros_like(br)
    for a0_ in range(0, 360, 10):
        m = br & (az >= a0_) & (az < a0_ + 10)
        if m.sum() > 30:
            bnd[m & (rad >= np.percentile(rad[m], 90))] = True
    print(f"--- {a.run}, trame {f[-8:-4]} : {int(br.sum())} rompus, "
          f"{int((br & para).sum())} paralleles au litage, {int((br & cross).sum())} en travers ---")
    for label, msk in (("PARALLELES a la FRONTIERE", br & para & bnd),
                       ("paralleles a l interieur", br & para & ~bnd),
                       ("en travers (tous)", br & cross)):
        n = int(msk.sum())
        if n == 0:
            print(f"  {label}: aucun"); continue
        print(f"  {label} ({n}) :")
        for s in sorted(set(st[msk].astype(int))):
            mm = msk & (st == s)
            print(f"     {STATE.get(s, s):24s} {100*mm.sum()/n:5.1f} %   sigN moyen {sig[mm].mean()/1e6:7.3f} MPa, |tauS| moyen {np.abs(tau[mm]).mean()/1e6:6.3f} MPa")
    # ---- evenements de traversee ---------------------------------------
    node_par = defaultdict(set)              # noeud -> ids de joints paralleles rompus
    for j in np.where(br & para)[0]:
        node_par[conn[j, 0]].add(j); node_par[conn[j, 1]].add(j)
    nCross = 0
    for j in np.where(br & cross)[0]:
        A = node_par.get(conn[j, 0], set()); B = node_par.get(conn[j, 1], set())
        if A and B and (A - B):                # deux plans DIFFERENTS aux deux bouts
            nCross += 1
    print(f"  EVENEMENTS DE TRAVERSEE (joint rompu en travers reliant deux plans rompus distincts) : {nCross}")
    print("  Lecture : frontiere majoritairement OUVERTE ou GLISSANTE -> arret mecanique (Renshaw-Pollard) ;")
    print("            frontiere FERMEE-BLOQUEE et zero traversee -> arret energetique (He-Hutchinson, rapport 0,057).")


if __name__ == "__main__":
    main()
