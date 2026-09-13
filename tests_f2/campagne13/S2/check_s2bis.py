#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""check_s2bis.py — complement FALSIFIANT de S2 (campagne du 13/09).

check_s2.py a verifie la cle contactForcePairs sur le banc s = 2,5 de 30 us
(C1-C6). Deux points du cadrage S2 y restaient sans chiffre :

  * le canal insert -> roche n'y est jamais NON nul (a 30 us l'onde n'a pas
    traverse le bit : C1 mesure 0 exact). Une colonne toujours nulle ne prouve
    pas qu'elle est branchee. -> C7, sur le banc s25_fcrock (train lance a
    -20 m/s, T = 4 us) : la MEME colonne doit devenir non nulle, du bon signe,
    antisymetrique, et recouper grpFz (trackGroup = insert).
  * la comparaison a l'estimateur « quantite de mouvement » de tools/fig_fp.py
    (F_r = m_p dv_p/dt + m_train dv_bit/dt, « la roche est la seule force
    exterieure au systeme piston + train »). -> C8 : sur le banc 30 us, cet
    estimateur est confronte a la somme des forces EXTERIEURES reellement
    mesurees par les colonnes Fc (plate->bit, rock->bit, rock->insert) et au
    poids. C9 : la meme quantite sur le temoin 200 us out_yang_bench_s25_v3P,
    ou les colonnes Fc n'existent pas (run anterieur a la cle).

Usage :
  python check_s2bis.py <out_s25_fcrock> <out_s25_fc> <out_temoin_200us>
"""
import csv
import os
import sys

import numpy as np

G = 9.81
# masses lues dans results/yang_bench_s25_v3P.log (memes corps, meme maillage
# s = 2,5 : le journal du run les reimprime, cf. RESULTATS_s25_fc.log)
M = {"rock": 19.2472, "insert": 0.0631335, "bit": 1.0802,
     "piston": 0.776709, "circlip": 0.0142247, "plate": 0.191708}

ok_all = True


def verdict(name, ok, msg):
    global ok_all
    ok_all = ok_all and ok
    print("%-4s %s : %s" % (name, "OK   " if ok else "ECHEC", msg))


def load(run):
    with open(os.path.join(run, "history.csv"), newline="") as f:
        rows = [x for x in csv.DictReader(f) if x and x.get("t")]
    out = {}
    for k in rows[0]:
        if k is None:
            continue
        try:
            out[k] = np.array([float(r[k]) for r in rows])
        except (TypeError, ValueError):
            pass
    return out


def masses_from_log(path):
    """relit les masses par corps dans un journal de run (si dispo)"""
    import re
    if not path or not os.path.isfile(path):
        return {}
    txt = open(path, errors="replace").read()
    return {m.group(1): float(m.group(2)) for m in
            re.finditer(r"corps '(\w+)':.*?masse = ([\d.eE+-]+) kg", txt)}


def smooth(y, t, tau):
    """moyenne glissante de largeur tau — identique a tools/fig_fp.py"""
    if len(t) < 5:
        return y
    dt = np.median(np.diff(t))
    n = max(1, int(round(tau / dt)))
    if n <= 1:
        return y
    return np.convolve(y, np.ones(n) / n, mode="same")


# --------------------------------------------------------------------------
# C7 : le canal insert -> roche est VIVANT (banc s25_fcrock)
# --------------------------------------------------------------------------
def c7(run):
    H = load(run)
    t = H["t"]
    fir = np.c_[H["Fc_insert_rock_x"], H["Fc_insert_rock_y"],
                H["Fc_insert_rock_z"]]
    fri = np.c_[H["Fc_rock_insert_x"], H["Fc_rock_insert_y"],
                H["Fc_rock_insert_z"]]
    nz = np.nonzero(fir[:, 2])[0]
    nneg = int(np.sum(fir[:, 2] < 0))
    if nz.size == 0:
        verdict("C7a", False, "Fc_insert_rock_z reste nul sur %d lignes — le "
                "canal insert/roche n est PAS branche" % len(t))
        return H
    i0, im = nz[0], int(np.argmax(np.abs(fir[:, 2])))
    verdict("C7a", fir[i0, 2] < 0 and np.abs(fir[im, 2]) > 1.0
            and nneg == nz.size,
            "premiere force non nulle a t = %.4g us (ligne %d) : "
            "Fc_insert_rock_z = %+.4g N ; max |Fc_insert_rock_z| = %.4g N a "
            "t = %.4g us ; lignes non nulles %d, dont %d negatives (l insert "
            "POUSSE la roche vers le bas : signe attendu negatif partout ; le "
            "meme canal donne 0 exact sur le banc 30 us, C1)"
            % (t[i0] * 1e6, i0, fir[i0, 2], np.abs(fir[im, 2]),
               t[im] * 1e6, nz.size, nneg))
    asym = float(np.max(np.abs(fir + fri)))
    verdict("C7b", asym == 0.0,
            "max |Fc_insert_rock + Fc_rock_insert| = %g N (attendu 0 exact, "
            "action = reaction sur les trois composantes)" % asym)
    part = ["Fc_rock_insert_z", "Fc_bit_insert_z", "Fc_circlip_insert_z",
            "Fc_plate_insert_z", "Fc_piston_insert_z"]
    part = [c for c in part if c in H]
    ssum = np.sum([H[c] for c in part], axis=0)
    d = np.abs(ssum - H["grpFz"])
    tol = np.maximum(1.0, 1e-6 * np.abs(H["grpFz"]) * 10)   # arrondi 6 chiffres
    ig = int(np.argmax(np.abs(H["grpFz"])))
    verdict("C7c", bool(np.all(d <= tol)),
            "somme %s == grpFz (trackGroup = insert) : ecart absolu max %.3g N "
            "sur %d lignes ; la somme n est pas triviale : max |grpFz| = %.4g N "
            "a t = %.4g us, somme au meme pas = %.4g N"
            % (part, float(np.max(d)), len(t), H["grpFz"][ig], t[ig] * 1e6,
               ssum[ig]))
    return H


# --------------------------------------------------------------------------
# C8 : l estimateur « quantite de mouvement » de fig_fp.py contre les Fc
# --------------------------------------------------------------------------
def c8(run, log=None):
    H = load(run)
    t = H["t"]
    m = dict(M)
    m.update(masses_from_log(log))
    mp = m["piston"]
    mt = m["bit"] + m["insert"] + m["circlip"]
    # (1) estimateur publie dans tools/fig_fp.py : tout le train suit le bit
    dP_fig = mp * np.gradient(H["vz_piston"], t) + mt * np.gradient(H["vz_bit"], t)
    Fr_fig = smooth(smooth(dP_fig, t, 2e-6), t, 2e-6)
    # (2) meme chose SANS approximation de corps rigide : m_g vz_g est la
    # quantite de mouvement EXACTE du corps g (vz_g est la moyenne ponderee
    # par les masses nodales, cf. historyRow). Le circlip n est pas suivi.
    P_ex = (mp * H["vz_piston"] + m["bit"] * H["vz_bit"]
            + m["insert"] * H["vz_insert"])
    Fr_ex = smooth(smooth(np.gradient(P_ex, t), t, 2e-6), t, 2e-6)
    # forces EXTERIEURES mesurees sur le systeme {piston, bit, insert, circlip}
    ext = np.zeros_like(t)
    used = []
    for c in ("Fc_plate_bit_z", "Fc_rock_bit_z", "Fc_rock_insert_z",
              "Fc_plate_insert_z", "Fc_rock_piston_z"):
        if c in H:
            ext = ext + H[c]
            used.append(c)
    poids = -(mp + mt) * G
    tot = ext + poids
    r_fig = Fr_fig - tot
    r_ex = Fr_ex - tot
    ech = max(float(np.max(np.abs(H["Fc_piston_bit_z"]))), 1.0)
    detail = ", ".join("%s max %.4g N" % (c, float(np.max(np.abs(H[c]))))
                       for c in used)
    verdict("C8a", len(used) >= 3,
            "forces EXTERIEURES au systeme {piston, bit, insert, circlip}, "
            "mesurees colonne par colonne : %s ; max |somme| = %.4g N ; poids "
            "= %.3g N. La roche ne porte encore RIEN a 30 us (C1) : tout vient "
            "du contact plaque/bit — l estimateur de fig_fp.py l attribue a la "
            "roche" % (detail, float(np.max(np.abs(ext))), poids))
    verdict("C8b", float(np.max(np.abs(r_ex))) < 0.05 * ech,
            "estimateur EXACT (m_g vz_g par corps) : max |F_r - F_ext| = %.4g N, "
            "RMS %.4g N ; echelle = max |Fc_piston_bit_z| = %.4g N "
            "(%.2f %% et %.2f %%)" % (float(np.max(np.abs(r_ex))),
                                      float(np.sqrt(np.mean(r_ex ** 2))), ech,
                                      100 * np.max(np.abs(r_ex)) / ech,
                                      100 * np.sqrt(np.mean(r_ex ** 2)) / ech))
    verdict("C8c", True,
            "estimateur de fig_fp.py (m_train dv_bit/dt) : max |F_r - F_ext| = "
            "%.4g N, RMS %.4g N, soit %.1f %% et %.1f %% de l echelle — biais "
            "de l hypothese « train rigide » (mesure, pas un verdict)"
            % (float(np.max(np.abs(r_fig))),
               float(np.sqrt(np.mean(r_fig ** 2))),
               100 * np.max(np.abs(r_fig)) / ech,
               100 * np.sqrt(np.mean(r_fig ** 2)) / ech))
    # controle FALSIFIANT : si le systeme est FAUX (piston oublie, la poussee
    # du piston devient exterieure), la fermeture C8b DOIT s effondrer.
    P_faux = m["bit"] * H["vz_bit"] + m["insert"] * H["vz_insert"]
    Fr_faux = smooth(smooth(np.gradient(P_faux, t), t, 2e-6), t, 2e-6)
    r_faux = Fr_faux - tot
    rap = float(np.max(np.abs(r_faux))) / max(float(np.max(np.abs(r_ex))), 1e-30)
    verdict("C8d", rap > 10.0,
            "controle falsifiant (systeme FAUX : piston exclu, sa poussee "
            "devient exterieure) : max |F_r - F_ext| = %.4g N, soit %.0f x le "
            "residu du systeme correct (doit exceder 10 x)"
            % (float(np.max(np.abs(r_faux))), rap))
    return Fr_fig, Fr_ex


# --------------------------------------------------------------------------
# C9 : la valeur a comparer plus tard, sur le temoin 200 us
# --------------------------------------------------------------------------
def c9(run, log=None):
    H = load(run)
    t = H["t"]
    m = dict(M)
    m.update(masses_from_log(log))
    mp = m["piston"]
    mt = m["bit"] + m["insert"] + m["circlip"]
    dP = mp * np.gradient(H["vz_piston"], t) + mt * np.gradient(H["vz_bit"], t)
    Fr = smooth(smooth(dP, t, 2e-6), t, 2e-6)
    P_ex = (mp * H["vz_piston"] + m["bit"] * H["vz_bit"]
            + m["insert"] * H["vz_insert"])
    Fr_ex = smooth(smooth(np.gradient(P_ex, t), t, 2e-6), t, 2e-6)
    i = int(np.argmax(Fr))
    j = int(np.argmax(Fr_ex))
    has = [c for c in H if c.startswith("Fc_")]
    verdict("C9", not has,
            "temoin %s : %d lignes jusqu a %.1f us ; estimateur fig_fp.py "
            "F_r max = %.4g N a t = %.1f us ; estimateur exact (3 corps) "
            "F_r max = %.4g N a t = %.1f us ; colonnes Fc dans ce run : %s "
            "(run anterieur a la cle S2 : la comparaison directe exige un "
            "rejeu >= 50 us, hors budget de 5 min)"
            % (os.path.basename(run.rstrip("/\\")), len(t), t[-1] * 1e6,
               Fr[i], t[i] * 1e6, Fr_ex[j], t[j] * 1e6, has or "aucune"))


def main():
    fcrock, fc, temoin = sys.argv[1:4]
    print("--- C7 : canal insert/roche (%s) ---" % fcrock)
    c7(fcrock)
    print("--- C8 : estimateur q. de mouvement contre les Fc (%s) ---" % fc)
    c8(fc, os.environ.get("LOG_FC"))
    print("--- C9 : temoin 200 us (%s) ---" % temoin)
    c9(temoin, "results/yang_bench_s25_v3P.log")
    print("BILAN S2bis :", "TOUS OK" if ok_all else "AU MOINS UN ECHEC")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
