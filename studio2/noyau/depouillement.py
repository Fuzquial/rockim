"""Grandeurs d'un essai triaxial terminé, avec les définitions de la campagne.

Mêmes définitions que etude_triax_hetero/depouille.py, réécrites en numpy ; le
test tests/test_depouillement.py vérifie l'égalité sur tous les runs présents.

  sigma3_atteint  pression latérale RÉELLEMENT établie au cœur : à lire avant le pic ;
  q_pic           déviateur au pic (écart à la fin de consolidation) ;
  E_secante       pente moindres carrés entre 30 et 60 % du pic, extensomètre central ;
  eps_pic         déformation au pic ;
  chute_post_pic  q à eps_pic + 0,2 % rapporté à q_pic : la fragilité ;
  modes           joints rompus en traction / en cisaillement ;
  intergranulaire part des joints rompus PENDANT l'essai qui sont des joints de grain.
"""
import csv
import math
import os
import re

import numpy as np

from .resultats import courbe, deformation_axiale, delai_consolidation, est_bresilien, lire_historique


def mesures(h, delai=0.0):
    if h is None or len(h["t"]) == 0:
        return None
    t = h["t"]
    sb = np.abs(h["sigma"])
    s0 = np.interp(delai, t, sb) if delai > 0 else 0.0
    sig = sb - s0
    eps = deformation_axiale(h)
    ipk = int(np.argmax(sig))
    spk = sig[ipk]
    conf = np.max(np.abs(h["confAchieved"])) if "confAchieved" in h else 0.0

    sel = np.zeros(len(t), bool)
    sel[:ipk + 1] = True
    sel &= (sig >= 0.30 * spk) & (sig <= 0.60 * spk) & (eps > 0.0)
    E = float("nan")
    if sel.sum() >= 2:
        x, y = eps[sel], sig[sel]
        n = len(x)
        den = n * (x * x).sum() - x.sum() ** 2
        if abs(den) > 0:
            E = (n * (x * y).sum() - x.sum() * y.sum()) / den

    post = np.nonzero(eps[ipk:] >= eps[ipk] + 0.002)[0]
    chute = sig[ipk + post[0]] / spk if len(post) and spk > 0 else float("nan")

    def fin(c):
        return int(h[c][-1]) if c in h else 0

    sf = mediane_glissante(sig, DEMI_FENETRE)
    jf = int(np.argmax(sf))
    return dict(sigma3_atteint_MPa=conf / 1e6, q_pic_MPa=spk / 1e6, sigma_pic_MPa=(spk + s0) / 1e6,
                E_secante_GPa=E / 1e9, eps_pic_pct=100.0 * eps[ipk], t_pic_ms=1e3 * t[ipk],
                chute_post_pic=chute, pic_verrouille=fin("peakLocked"), n_rompus=fin("nBroken"),
                n_tension=fin("nBrokTen"), n_cisaillement=fin("nBrokShear"), n_insere=fin("nInserted"),
                n_fragments=fin("nFrag"), t_fin_ms=1e3 * t[-1],
                q_pic_filtre_MPa=sf[jf] / 1e6, eps_pic_filtre_pct=100.0 * eps[jf])


# Pic FILTRÉ (spec 007 §2.6, mesure du 2026-10-02) : sur le banc de 10 grains, des pointes
# isolées de contrainte (16 lignes sur 1 103 en UCS) faisaient passer le maximum brut de 48 à
# 75 MPa. La médiane glissante sur ±7 lignes d'historique (sur ~2 000) les écarte sans
# déplacer le pic de l'enveloppe (stable de ±3 à ±15 lignes). Le maximum brut reste
# rendu (q_pic_MPa), pour l'identité avec depouille.py.
DEMI_FENETRE = 7


def mediane_glissante(y, k):
    if len(y) == 0:
        return y
    p = np.pad(y, k, mode="edge")
    return np.median(np.lib.stride_tricks.sliding_window_view(p, 2 * k + 1), axis=1)


def joints(chemin):
    """Joints rompus PENDANT l'essai (tBreak > 0) : les pré-fissures démarrent à
    damage = 1 avec tBreak = 0 et sont comptées à part (prerompus_csv)."""
    if not os.path.exists(chemin):
        return {}
    with open(chemin, encoding="utf-8", errors="replace") as f:
        r = list(csv.DictReader(f))
    if not r:
        return {}
    col = lambda k: np.array([float(x[k]) for x in r])
    dmg, tb, ty = col("damage"), col("tBreak"), col("type").astype(int)
    pre = (dmg >= 0.999) & (tb <= 0.0)
    ro = (tb > 0.0) & ~pre
    ty = np.minimum(ty[ro], 2)
    out = dict(rompus=int(ro.sum()), intra=int((ty == 0).sum()), homo=int((ty == 1).sum()),
               hetero=int((ty == 2).sum()), prerompus_csv=int(pre.sum()))
    if out["rompus"]:
        dx = np.abs(col("x2") - col("x1"))[ro]
        dy = np.abs(col("y2") - col("y1"))[ro]
        ang = np.sort(np.degrees(np.arctan2(dy, dx)))
        out["part_intergranulaire"] = (out["homo"] + out["hetero"]) / out["rompus"]
        out["angle_median_deg"] = float(ang[len(ang) // 2])
    return out


def diagnostics(chemin_log):
    """Ce que le solveur imprime à l'initialisation et qu'aucun CSV ne porte."""
    if not os.path.exists(chemin_log):
        return {}
    txt = open(chemin_log, encoding="utf-8", errors="replace").read()
    d = {}
    m = re.search(r"population pre-endommagee : (\d+) / (\d+)", txt)
    if m:
        d["prerompus"] = int(m.group(1))
        d["prerompus_pct"] = 100.0 * int(m.group(1)) / int(m.group(2))
    m = re.search(r"; (\d+) avec une extremite deja scindee", txt)
    if m:
        d["prerompus_libres"] = int(m.group(1))
    m = re.search(r"-> (\d+) / (\d+) joints affaiblis", txt)
    if m:
        d["joints_affaiblis"] = int(m.group(1))
    m = re.search(r"(\d+) elements, (\d+) joints", txt)
    if m:
        d["elements"], d["joints"] = int(m.group(1)), int(m.group(2))
    m = re.search(r"voronoi: (\d+) grains, \d+ phase\(s\), (\d+) grain-boundary joints, hmin = ([\d.e+-]+)", txt)
    if m:
        d["grains"], d["joints_de_grain"], d["hmin_m"] = int(m.group(1)), int(m.group(2)), float(m.group(3))
    else:
        m = re.search(r"voronoi: (\d+) grains", txt)
        if m:
            d["grains"] = int(m.group(1))
    # aires de phase réalisées contre visées (fin du journal)
    m = re.search(r"grains: \d+, phases: (.*)", txt)
    if m:
        d["phases_realisees"] = [{"nom": n, "realise_pct": float(r), "cible_pct": float(c)}
                                 for n, r, c in re.findall(r"(\w+) ([\d.]+)% \(target ([\d.]+)%\)", m.group(1))]
    if "WARNING" in txt:
        d["avertissements"] = txt.count("WARNING")
    return d


def mesures_bresilien(h):
    """Brésilien : sigma_t = 2P/(pi D t) écrit par le solveur à chaque ligne d'historique.
    Le pic brut, le pic filtré (médiane glissante, comme en triaxial), la force et la
    fermeture des plateaux au pic."""
    if h is None or len(h["t"]) == 0:
        return None
    ferm, st = courbe(h)
    st = st * 1e6
    i = int(np.argmax(st))
    sf = mediane_glissante(st, DEMI_FENETRE)
    j = int(np.argmax(sf))

    def fin(c):
        return int(h[c][-1]) if c in h else 0

    return dict(essai="bresilien", sigma_t_pic_MPa=st[i] / 1e6, sigma_t_pic_filtre_MPa=sf[j] / 1e6,
                P_pic_N=float(np.abs(h["P"][i])), fermeture_pic_mm=float(ferm[i]), t_pic_ms=1e3 * h["t"][i],
                pic_verrouille=fin("peakLocked"), n_rompus=fin("nBroken"), n_fragments=fin("nFrag"),
                t_fin_ms=1e3 * h["t"][-1])


def diagnostics_bresilien(chemin_log):
    """Les verdicts du banc brésilien, imprimés par le solveur en fin de run
    (FdemSolver.cpp, résumé « brazilian ») :
      jauge élastique de BANDE : sigma_xx au centre / sigma_t sur sigma_t dans [0,3 ; 0,8] ft,
        attendue à 1 (solution fermée du disque), PASS dans [0,85 ; 1,25] ;
      diamétralité : part des joints rompus près de l'axe de charge ;
      sigma_t au pic et son rapport à ft."""
    if not os.path.exists(chemin_log):
        return {}
    txt = open(chemin_log, encoding="utf-8", errors="replace").read()
    d = {}
    m = re.search(r"-> ratio ([-\d.eE+]+)\s+\[(PASS|FAIL)\]\s+\(band 0\.85-1\.25; mean sigma_yy = ([-\d.eE+]+) MPa, "
                  r"sigma_yy/sigma_xx = ([-\d.eE+]+)\)", txt)
    if m:
        d.update(jauge_elastique=float(m.group(1)), jauge_elastique_verdict=m.group(2),
                 syy_sur_sxx=float(m.group(4)))
    elif "ELASTIC-BAND gauge: NOT measured" in txt:
        d["jauge_elastique_verdict"] = "non mesurée"
    m = re.search(r"sigma_t = 2P/\(pi D t\) = ([-\d.eE+]+) MPa", txt)
    if m:
        d["sigma_t_solveur_MPa"] = float(m.group(1))
    m = re.search(r"ratio to the bulk ft \(([-\d.eE+]+) MPa\) = ([-\d.eE+]+)", txt)
    if m:
        d["sigma_t_sur_ft"] = float(m.group(2))
    m = re.search(r"\(([\d.]+) % diametral\)", txt)
    if m:
        d["diametral_pct"] = float(m.group(1))
    m = re.search(r"peak force P = ([-\d.eE+]+) N", txt)
    if m:
        d["P_pic_solveur_N"] = float(m.group(1))
    d["pic_verrouille_solveur"] = "peak LOCKED" in txt
    m = re.search(r"early stop at t = [-\d.eE+]+ s \((\d+) / (\d+) steps\)", txt)
    if m:
        d["pas_effectues"], d["pas_plafond"] = int(m.group(1)), int(m.group(2))
    return d


def synthese(dossier_run, chemin_log=None):
    """Toutes les grandeurs d'un run, prêtes pour le tableau comparatif."""
    h = lire_historique(dossier_run)
    if est_bresilien(h):
        m = mesures_bresilien(h)
        m.update(joints(os.path.join(dossier_run, "fdem_final_joints.csv")))
        if chemin_log:
            m.update(diagnostics(chemin_log))
            m.update(diagnostics_bresilien(chemin_log))
        return m
    m = mesures(h, delai_consolidation(dossier_run))
    if m is None:
        return None
    m.update(joints(os.path.join(dossier_run, "fdem_final_joints.csv")))
    if chemin_log:
        m.update(diagnostics(chemin_log))
    return m


def mode_dominant(m):
    return "traction" if m["n_tension"] > m["n_cisaillement"] else "cisaillement"


def est_nan(x):
    return isinstance(x, float) and math.isnan(x)
