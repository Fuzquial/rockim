#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# yang_estimators.py — les estimateurs de Yang et al. (2025, §4.1 et fig. 8)
# pour la cinetique d'un impact rockim fdem3d, lus dans history.csv.
#
# Campagne de correction du 13/09/2026, tache T2 (docs/CAMPAGNE_correction_
# 2026-09-13.md ; motif : COMPLEMENT_YANG §5). Yang mesure :
#   * la vitesse d'indentation = PENTE de la portion lineaire de la courbe
#     deplacement-temps du bit (ici : entre 10 % et 90 % de l'enfoncement
#     maximal, avant le retournement) ;
#   * la vitesse de rebond = PENTE de la portion lineaire remontante apres le
#     retournement (ici : entre 10 % et 90 % du deplacement recupere dans
#     l'enregistrement) ;
#   * la contrainte de reference = PIC DE LA PREMIERE ONDE a la jauge a mi-bit,
#     pas le plus grand pic de toute la simulation.
# Les extrema instantanes (min de vz_<corps>, max de -szz_bit) restent
# calcules et imprimes COTE A COTE : ce sont deux estimateurs differents de la
# meme grandeur, et l'ecart entre eux est une information (ECARTS §5 : 6,86
# contre 7,37 m/s sur le s = 1).
#
# Conventions rockim : z vers le haut, le bit descend (vz < 0) ; enfoncement
# p = z(0) - z(t) > 0 vers le bas ; szz_bit < 0 en compression, on travaille
# sur sig = -szz_bit > 0.
#
# Aucune cle solveur : module de post-traitement pur (numpy seul), partage par
# tools/fig_kinetics.py et tools/yang_report.py. Test falsifiant :
# tools/test_yang_estimators.py.
# ---------------------------------------------------------------------------
import csv
import io
import math

import numpy as np

NAN = float("nan")


# ---------------------------------------------------------------- lecture --
def load_history(run_or_path):
    """history.csv -> dict colonne -> np.ndarray (float). Accepte le dossier
    du run ou le chemin du CSV."""
    path = run_or_path
    if not path.lower().endswith(".csv"):
        path = path.rstrip("/\\") + "/history.csv"
    with io.open(path, newline="") as f:
        rows = [r for r in csv.DictReader(f) if r]
    if not rows:
        return {}
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def penetration(z):
    """Enfoncement p = z(0) - z(t), positif vers le bas (m)."""
    z = np.asarray(z, dtype=float)
    return z[0] - z


# ------------------------------------------------------------- utilitaires --
def _linfit(t, y):
    """Droite y = a t + b aux moindres carres ; renvoie (a, b, r2, rms)."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(t) < 2:
        return NAN, NAN, NAN, NAN
    tm = t.mean()
    ym = y.mean()
    dt = t - tm
    den = float(np.dot(dt, dt))
    if den <= 0.0:
        return NAN, NAN, NAN, NAN
    a = float(np.dot(dt, y - ym)) / den
    b = ym - a * tm
    res = y - (a * t + b)
    ss_res = float(np.dot(res, res))
    ss_tot = float(np.dot(y - ym, y - ym))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else NAN
    rms = math.sqrt(ss_res / len(t))
    return a, b, r2, rms


def _first_index(mask, start=0):
    """Premier indice >= start ou mask est vrai, -1 sinon."""
    idx = np.nonzero(mask[start:])[0]
    return int(idx[0]) + start if len(idx) else -1


# ------------------------------------------------------- vitesse (pentes) --
def indentation_slope(t, p, frac=(0.10, 0.90)):
    """Vitesse d'indentation = pente de p(t) entre frac[0] et frac[1] de
    l'enfoncement maximal, avant le retournement (Yang 2025 §4.1).

    Renvoie un dict : v (m/s, > 0 vers le bas), t0, t1 (s) bornes de la
    fenetre, i0, i1 indices, r2, rms (m) du residu, p_max (m), t_pmax (s),
    n points, ok (bool), reason (str)."""
    t = np.asarray(t, dtype=float)
    p = np.asarray(p, dtype=float)
    out = dict(v=NAN, t0=NAN, t1=NAN, i0=-1, i1=-1, r2=NAN, rms=NAN,
               p_max=NAN, t_pmax=NAN, n=0, ok=False, reason="",
               frac=tuple(frac))
    if len(t) < 3:
        out["reason"] = "moins de 3 points"
        return out
    i_turn = int(np.argmax(p))
    p_max = float(p[i_turn])
    out["p_max"], out["t_pmax"] = p_max, float(t[i_turn])
    if not p_max > 0.0:
        out["reason"] = "enfoncement max nul ou negatif"
        return out
    head = p[: i_turn + 1]
    i0 = _first_index(head >= frac[0] * p_max)
    i1 = _first_index(head >= frac[1] * p_max)
    if i0 < 0 or i1 < 0 or i1 - i0 < 2:
        out["reason"] = "fenetre %g-%g %% trop courte (%d points)" % (
            100 * frac[0], 100 * frac[1], max(0, i1 - i0 + 1))
        return out
    a, b, r2, rms = _linfit(t[i0:i1 + 1], p[i0:i1 + 1])
    out.update(v=a, t0=float(t[i0]), t1=float(t[i1]), i0=i0, i1=i1, r2=r2,
               rms=rms, n=i1 - i0 + 1, ok=True)
    return out


def rebound_slope(t, p, frac=(0.10, 0.90), min_amp=0.05):
    """Vitesse de rebond = pente de la portion lineaire REMONTANTE de p(t)
    apres le retournement (Yang 2025 §4.1, lue apres 450 us dans leur run).

    Fenetre : deplacement recupere u = p_max - p(t) entre frac[0] et frac[1]
    de l'amplitude recuperee A = p_max - p(fin d'enregistrement). Le rebond est
    declare NON MESURABLE si le retournement n'est pas dans l'enregistrement
    (p_max au dernier point) ou si A < min_amp * p_max : on ne fabrique pas
    une pente sur un retournement inacheve. Si l'enregistrement s'arrete en
    plein rebond, la fenetre est relative a l'amplitude ENREGISTREE : lire r2
    et comparer a la vitesse instantanee finale.

    Renvoie : v (m/s, > 0 vers le haut), t0, t1, i0, i1, r2, rms, amp (m),
    t_turn (s), ok, reason."""
    t = np.asarray(t, dtype=float)
    p = np.asarray(p, dtype=float)
    out = dict(v=NAN, t0=NAN, t1=NAN, i0=-1, i1=-1, r2=NAN, rms=NAN, amp=NAN,
               t_turn=NAN, ok=False, reason="", frac=tuple(frac),
               min_amp=min_amp)
    if len(t) < 3:
        out["reason"] = "moins de 3 points"
        return out
    i_turn = int(np.argmax(p))
    p_max = float(p[i_turn])
    out["t_turn"] = float(t[i_turn])
    if i_turn >= len(t) - 1:
        out["reason"] = "pas de retournement dans l'enregistrement (p_max au dernier point)"
        return out
    amp = p_max - float(p[-1])
    out["amp"] = amp
    if not p_max > 0.0 or amp < min_amp * p_max:
        out["reason"] = ("retournement inacheve : %.1f %% de p_max recuperes "
                         "(< %g %%)" % (100.0 * amp / p_max if p_max > 0 else 0.0,
                                        100.0 * min_amp))
        return out
    u = p_max - p
    i0 = _first_index(u >= frac[0] * amp, i_turn)
    i1 = _first_index(u >= frac[1] * amp, i_turn)
    if i0 < 0 or i1 < 0 or i1 - i0 < 2:
        out["reason"] = "fenetre %g-%g %% trop courte (%d points)" % (
            100 * frac[0], 100 * frac[1], max(0, i1 - i0 + 1))
        return out
    a, b, r2, rms = _linfit(t[i0:i1 + 1], p[i0:i1 + 1])
    out.update(v=-a, t0=float(t[i0]), t1=float(t[i1]), i0=i0, i1=i1, r2=r2,
               rms=rms, n=i1 - i0 + 1, ok=True)
    return out


# ------------------------------------------------ vitesse (max instantane) --
def instantaneous(t, vz):
    """L'ancien estimateur : extrema instantanes de la vitesse moyenne du
    corps (rockim vz < 0 vers le bas). v_ind = -min(vz) ; v_reb = max(vz)
    APRES l'instant du minimum (> 0 seulement si le corps remonte)."""
    t = np.asarray(t, dtype=float)
    vz = np.asarray(vz, dtype=float)
    out = dict(v_ind=NAN, t_ind=NAN, v_reb=NAN, t_reb=NAN, v_end=NAN,
               ok_reb=False)
    if len(t) == 0:
        return out
    i = int(np.argmin(vz))
    j = i + int(np.argmax(vz[i:]))
    out.update(v_ind=-float(vz[i]), t_ind=float(t[i]), v_reb=float(vz[j]),
               t_reb=float(t[j]), v_end=float(vz[-1]), ok_reb=bool(vz[j] > 0.0))
    return out


# -------------------------------------------------- contrainte de reference --
def first_wave_peak(t, sig, onset_frac=0.05, end_frac=0.10):
    """Pic de la PREMIERE onde d'un signal de jauge sig(t) >= 0 (compression
    positive). Debut = premier point ou sig > onset_frac * max global ;
    fin = premier point suivant ou sig retombe sous end_frac * (max courant
    depuis le debut) ; pic = max sur [debut, fin]. Le max GLOBAL et sa date
    sont renvoyes a cote (l'ancien estimateur).

    Renvoie : sig (pic 1re onde), t (s), t_on, t_off, i_on, i_off,
    finished (bool : la 1re onde est retombee avant la fin), sig_max, t_max,
    ok, reason."""
    t = np.asarray(t, dtype=float)
    sig = np.asarray(sig, dtype=float)
    out = dict(sig=NAN, t=NAN, t_on=NAN, t_off=NAN, i_on=-1, i_off=-1,
               finished=False, sig_max=NAN, t_max=NAN, ok=False, reason="",
               onset_frac=onset_frac, end_frac=end_frac)
    if len(t) == 0:
        out["reason"] = "signal vide"
        return out
    ig = int(np.argmax(sig))
    smax = float(sig[ig])
    out["sig_max"], out["t_max"] = smax, float(t[ig])
    if not smax > 0.0:
        out["reason"] = "aucune compression a la jauge"
        return out
    i_on = _first_index(sig > onset_frac * smax)
    if i_on < 0:
        out["reason"] = "pas de depart d'onde"
        return out
    run_max = np.maximum.accumulate(sig[i_on:])
    below = sig[i_on:] < end_frac * run_max
    below[0] = False
    k = _first_index(below)
    if k < 0:
        i_off = len(t) - 1
        finished = False
    else:
        i_off = i_on + k
        finished = True
    ip = i_on + int(np.argmax(sig[i_on:i_off + 1]))
    out.update(sig=float(sig[ip]), t=float(t[ip]), t_on=float(t[i_on]),
               t_off=float(t[i_off]), i_on=i_on, i_off=i_off,
               finished=finished, ok=True)
    return out


# ------------------------------------------------------------ depouillement --
def kinetics(hist, bodies=("insert", "bit"), gauge="szz_bit",
             frac=(0.10, 0.90), min_amp=0.05):
    """Les deux estimateurs pour chaque corps suivi (colonnes z_<corps>,
    vz_<corps>) et pour la jauge. Renvoie un dict : t, bodies -> {ind, reb,
    inst, p}, gauge -> first_wave_peak."""
    t = hist["t"]
    res = dict(t=t, bodies={}, gauge=None)
    for b in bodies:
        if "z_" + b not in hist:
            continue
        p = penetration(hist["z_" + b])
        vz = hist.get("vz_" + b)
        res["bodies"][b] = dict(
            p=p,
            ind=indentation_slope(t, p, frac),
            reb=rebound_slope(t, p, frac, min_amp),
            inst=instantaneous(t, vz) if vz is not None else None)
    if gauge in hist:
        res["gauge"] = first_wave_peak(t, -hist[gauge])
    return res


def _us(x):
    return x * 1e6


# ------------------------------------------- T2 deuxieme passe (13/09 soir) --
# Les trois fonctions qui suivent ne changent RIEN au chemin par defaut :
# elles repondent a la question « la pente 10-90 % est-elle une mesure ou un
# artefact de fenetre ? » et « le pic retenu est-il bien celui de la premiere
# onde, au sens des temps de transit du train ? ». Elles sont appelees par
# l'option --robust du module et par la section F de
# tools/test_yang_estimators.py ; sans l'option, la sortie est inchangee.

def slope_theilsen(t, y, max_points=300):
    """Pente robuste (Theil-Sen) : mediane des pentes de toutes les paires de
    points, sur un sous-echantillonnage a max_points points au plus (le nombre
    de paires croit en n^2). Estimateur INDEPENDANT des moindres carres :
    s'ils divergent, la fenetre n'est pas lineaire ou un point aberrant tire la
    droite. Renvoie (pente, nombre de paires)."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(t)
    if n < 2:
        return NAN, 0
    step = max(1, n // int(max_points))
    ti = t[::step]
    yi = y[::step]
    m = len(ti)
    sl = []
    for a in range(m):
        dt = ti[a + 1:] - ti[a]
        dy = yi[a + 1:] - yi[a]
        ok = dt > 0.0
        if np.any(ok):
            sl.append(dy[ok] / dt[ok])
    if not sl:
        return NAN, 0
    allsl = np.concatenate(sl)
    return float(np.median(allsl)), int(len(allsl))


DEFAULT_FRACS = ((0.05, 0.95), (0.10, 0.90), (0.20, 0.80), (0.30, 0.70),
                 (0.40, 0.60))


def window_sensitivity(t, p, fracs=DEFAULT_FRACS, ref=(0.10, 0.90)):
    """Pente d'indentation pour plusieurs fenetres (a, b) de l'enfoncement.
    Renvoie dict(rows=[{frac, v, t0, t1, n, r2, dev}], v_ref, spread_pct,
    spread_pct_wide) : `dev` = ecart relatif a la fenetre de reference (%),
    `spread_pct` = ecart max sur les fenetres LARGES (celles dont la demi-plage
    est >= 0,2, c.-a-d. 5-95, 10-90, 20-80), `spread_pct_wide` = ecart max sur
    TOUTES les fenetres demandees. Une pente qui bouge de plus de quelques
    pourcents entre 5-95 % et 20-80 % n'est pas une droite : p(t) accelere
    puis decelere, et le chiffre publie doit etre donne avec sa fenetre."""
    t = np.asarray(t, dtype=float)
    p = np.asarray(p, dtype=float)
    d0 = indentation_slope(t, p, ref)
    out = dict(rows=[], v_ref=d0["v"], ref=tuple(ref), spread_pct=NAN,
               spread_pct_wide=NAN)
    dev_wide, dev_all = [], []
    for fr in fracs:
        d = indentation_slope(t, p, fr)
        dev = (100.0 * abs(d["v"] - d0["v"]) / abs(d0["v"])
               if d["ok"] and d0["ok"] and d0["v"] != 0.0 else NAN)
        out["rows"].append(dict(frac=tuple(fr), v=d["v"], t0=d["t0"], t1=d["t1"],
                                n=d.get("n", 0), r2=d["r2"], dev=dev,
                                ok=d["ok"], reason=d["reason"]))
        if d["ok"] and dev == dev:
            dev_all.append(dev)
            if (fr[1] - fr[0]) >= 0.6 - 1e-12:
                dev_wide.append(dev)
    if dev_wide:
        out["spread_pct"] = max(dev_wide)
    if dev_all:
        out["spread_pct_wide"] = max(dev_all)
    return out


def wave_transit(z_gauge=(0.28, 0.31), z_bit=(0.17322, 0.41502),
                 L_piston=0.26, gap=2e-5, v0=9.0, E=200e9, rho=7850.0,
                 nu=0.29):
    """Temps de transit 1D du train de frappe, pour juger la fenetre de
    premiere onde de `first_wave_peak`. Geometrie par defaut = celle que le
    solveur imprime pour les maillages yang2026 (`bit : z = [0.17322,
    0.41502] m`, `gauge bit : ... z = [0.28, 0.31] m`, piston 0,26 m, jeu
    0,02 mm, 9 m/s) et phase acier du deck (200 GPa, 7850, 0,29) ; tout est
    surchargeable (voir transit_from_log).

    Renvoie : c_bar = sqrt(E/rho) (onde de barre, celle qui porte le palier),
    c_p = sqrt(E(1-nu)/(rho(1+nu)(1-2nu))) (dilatation, front le plus rapide
    d'un milieu massif), t_contact = jeu / v0, t_arr_bar / t_arr_p (arrivee au
    CENTRE de la bande de jauge), t_first_p (arrivee au bord HAUT de la bande,
    a c_p : la borne physique la plus precoce d'un signal a la jauge),
    t_refl (retour au centre de la bande de l'onde reflechie par le bas du
    bit, a c_bar : au-dela, le signal de la jauge est une superposition),
    T_pulse = 2 L_piston / c_bar (duree du creneau incident), sig_1d =
    rho c_bar v0 / 2 (contrainte d'impact 1D de deux barres de meme
    impedance)."""
    c_bar = math.sqrt(E / rho)
    c_p = math.sqrt(E * (1.0 - nu) / (rho * (1.0 + nu) * (1.0 - 2.0 * nu)))
    zg = 0.5 * (z_gauge[0] + z_gauge[1])
    d_top = z_bit[1] - zg
    d_top_edge = z_bit[1] - z_gauge[1]
    d_bot = zg - z_bit[0]
    t_c = gap / v0 if v0 > 0 else 0.0
    return dict(c_bar=c_bar, c_p=c_p, t_contact=t_c,
                t_arr_bar=t_c + d_top / c_bar, t_arr_p=t_c + d_top / c_p,
                t_first_p=t_c + d_top_edge / c_p,
                t_refl=t_c + d_top / c_bar + 2.0 * d_bot / c_bar,
                T_pulse=2.0 * L_piston / c_bar,
                sig_1d=rho * c_bar * v0 / 2.0, z_gauge=tuple(z_gauge),
                z_bit=tuple(z_bit), L_piston=L_piston, gap=gap, v0=v0,
                E=E, rho=rho, nu=nu)


def transit_from_log(path, **kw):
    """Meme chose, mais la geometrie est LUE dans le journal du solveur (lignes
    `[FDEM3D]   bit : z = [a, b] m`, `[FDEM3D]   piston : z = [a, b] m`,
    `[FDEM3D] gauge bit : n tets dans z = [z0, z1] m`, `groupVel.piston =
    (0, 0, -v)`) au lieu d'etre supposee. Les cles E, rho, nu (phase acier du
    deck) restent des arguments. Renvoie (transit, source) ou source est un
    dict des valeurs effectivement trouvees ; les absentes gardent le defaut."""
    import re
    src = {}
    try:
        txt = io.open(path, errors="replace").read()
    except Exception as exc:                                   # noqa: BLE001
        return wave_transit(**kw), dict(erreur=str(exc))
    m = re.search(r"bit : z = \[([-0-9.e+]+), ([-0-9.e+]+)\] m", txt)
    if m:
        src["z_bit"] = (float(m.group(1)), float(m.group(2)))
    m = re.search(r"piston : z = \[([-0-9.e+]+), ([-0-9.e+]+)\] m", txt)
    if m:
        z0, z1 = float(m.group(1)), float(m.group(2))
        src["L_piston"] = z1 - z0
        if "z_bit" in src:
            src["gap"] = max(0.0, z0 - src["z_bit"][1])
    m = re.search(r"gauge \w+ : \d+ tets dans z = \[([-0-9.e+]+), ([-0-9.e+]+)\] m", txt)
    if m:
        src["z_gauge"] = (float(m.group(1)), float(m.group(2)))
    m = re.search(r"groupVel\.piston = \(0, 0, (-?[0-9.e+]+)\) m/s", txt)
    if m:
        src["v0"] = abs(float(m.group(1)))
    args = dict(src)
    args.update(kw)
    return wave_transit(**args), src


def format_robustness(res, transit=None, yang_vind=5.62, yang_sig=160.0,
                      fracs=DEFAULT_FRACS):
    """Bloc texte du controle croise : pente robuste, sensibilite a la fenetre,
    estimateurs naifs, et coherence de la premiere onde avec les transits."""
    L = ["--- controle croise T2 (pente robuste, fenetres, transits) ---"]
    t = res["t"]
    for b, d in res["bodies"].items():
        ind, p = d["ind"], d["p"]
        if not ind["ok"]:
            L.append("%-7s pente non mesurable (%s)" % (b, ind["reason"]))
            continue
        ts, npair = slope_theilsen(t[ind["i0"]:ind["i1"] + 1],
                                   p[ind["i0"]:ind["i1"] + 1])
        dev_ts = 100.0 * abs(ts - ind["v"]) / abs(ind["v"])
        L.append("%-7s pente MC %.4f m/s | Theil-Sen %.4f m/s sur la MEME fenetre "
                 "(%d paires, ecart %.2f %%)" % (b, ind["v"], ts, npair, dev_ts))
        ws = window_sensitivity(t, p, fracs)
        for r in ws["rows"]:
            L.append("        fenetre %2.0f-%2.0f %% : %7.4f m/s  [%6.2f-%6.2f us, "
                     "%4d pts, r2 %.5f]  ecart %5.2f %%"
                     % (100 * r["frac"][0], 100 * r["frac"][1], r["v"],
                        _us(r["t0"]), _us(r["t1"]), r["n"], r["r2"], r["dev"]))
        L.append("        dispersion des fenetres LARGES (>= 60 %% de p_max) %.2f %% ; "
                 "toutes fenetres %.2f %%" % (ws["spread_pct"], ws["spread_pct_wide"]))
        naive = ind["p_max"] / ind["t_pmax"] if ind["t_pmax"] > 0 else NAN
        L.append("        estimateurs naifs : p_max/t_pmax %.4f m/s (ecart %.1f %%)"
                 % (naive, 100.0 * abs(naive - ind["v"]) / abs(ind["v"])))
        if d["inst"] is not None:
            L.append("                            max instantane %.4f m/s (ecart %.1f %%)"
                     % (d["inst"]["v_ind"],
                        100.0 * abs(d["inst"]["v_ind"] - ind["v"]) / abs(ind["v"])))
        if yang_vind == yang_vind:
            L.append("        Yang 5,62 m/s : ecart %+.1f %% a la pente 10-90 %%, "
                     "%+.1f %% a la fenetre la plus FAVORABLE (%.4f m/s)"
                     % (100.0 * (ind["v"] - yang_vind) / yang_vind,
                        100.0 * (min(r["v"] for r in ws["rows"] if r["ok"]) - yang_vind)
                        / yang_vind, min(r["v"] for r in ws["rows"] if r["ok"])))
    g = res["gauge"]
    if transit is not None and g is not None and g["ok"]:
        tr = transit
        L.append("jauge   transits 1D (c_barre %.0f m/s, c_P %.0f m/s, contact a %.2f us) :"
                 % (tr["c_bar"], tr["c_p"], _us(tr["t_contact"])))
        L.append("        premier signal possible a la bande %.2f us (bord haut, c_P) ; "
                 "arrivee au centre %.2f us (c_barre)"
                 % (_us(tr["t_first_p"]), _us(tr["t_arr_bar"])))
        L.append("        retour de la reflexion du bas du bit %.2f us ; creneau de "
                 "piston 2L/c %.1f us" % (_us(tr["t_refl"]), _us(tr["T_pulse"])))
        L.append("        mesure : depart %.2f us, pic %.2f us, fin de fenetre %.2f us "
                 "(duree %.1f us)" % (_us(g["t_on"]), _us(g["t"]), _us(g["t_off"]),
                                      _us(g["t_off"] - g["t_on"])))
        ok_peak = tr["t_arr_bar"] <= g["t"] <= tr["t_refl"]
        L.append("        pic dans [arrivee, retour de reflexion] : %s  -> le pic de "
                 "reference %s contamine par la reflexion du bit"
                 % ("OUI" if ok_peak else "NON", "n'est PAS" if ok_peak else "PEUT etre"))
        L.append("        pic %.2f MPa contre l'impact 1D rho c v/2 = %.2f MPa "
                 "(ecart %+.1f %%) ; Yang ~%.0f MPa (ecart %+.1f %%)"
                 % (g["sig"] / 1e6, tr["sig_1d"] / 1e6,
                    100.0 * (g["sig"] - tr["sig_1d"]) / tr["sig_1d"],
                    yang_sig, 100.0 * (g["sig"] / 1e6 - yang_sig) / yang_sig))
        if g["t_off"] - g["t_on"] < tr["T_pulse"]:
            L.append("        NB : la fenetre mesuree (%.1f us) est plus COURTE que le "
                     "creneau incident (%.1f us) : la jauge est dechargee avant la fin "
                     "du creneau (reflexions), le critere de retombee ne delimite donc "
                     "pas le creneau complet." % (_us(g["t_off"] - g["t_on"]),
                                                  _us(tr["T_pulse"])))
    return "\n".join(L)


def format_report(res, label="", yang=None):
    """Tableau texte ASCII : les deux estimateurs cote a cote avec leurs
    fenetres. yang = dict(vind, vreb, sig) pour rappeler les reperes."""
    L = []
    if label:
        L.append("--- estimateurs de Yang, pente 10-90 %% | max instantane : %s ---" % label)
    for b, d in res["bodies"].items():
        ind, reb, inst = d["ind"], d["reb"], d["inst"]
        s_ind = ("%6.3f m/s  [%5.1f-%5.1f us, %d pts, r2 %.4f, rms %.1e mm]" % (
            ind["v"], _us(ind["t0"]), _us(ind["t1"]), ind["n"], ind["r2"],
            ind["rms"] * 1e3) if ind["ok"] else "n.m. (%s)" % ind["reason"])
        s_reb = ("%6.3f m/s  [%5.1f-%5.1f us, %d pts, r2 %.4f, %.0f %% de p_max recuperes]"
                 % (reb["v"], _us(reb["t0"]), _us(reb["t1"]), reb["n"], reb["r2"],
                    100.0 * reb["amp"] / ind["p_max"])
                 if reb["ok"] else "n.m. (%s)" % reb["reason"])
        if inst is not None:
            s_iind = "%6.3f m/s  [%5.1f us]" % (inst["v_ind"], _us(inst["t_ind"]))
            s_ireb = ("%6.3f m/s  [%5.1f us]" % (inst["v_reb"], _us(inst["t_reb"]))
                      if inst["ok_reb"] else
                      "n.m. (vz final %+.3f m/s, encore descendant)" % inst["v_end"])
        else:
            s_iind = s_ireb = "colonne vz absente"
        L.append("%-7s v_ind pente %s" % (b, s_ind))
        L.append("%-7s v_ind max   %s" % ("", s_iind))
        L.append("%-7s v_reb pente %s" % ("", s_reb))
        L.append("%-7s v_reb max   %s" % ("", s_ireb))
        L.append("%-7s p_max %.4f mm a %.1f us" % ("", ind["p_max"] * 1e3, _us(ind["t_pmax"])))
    g = res["gauge"]
    if g is not None and g["ok"]:
        L.append("jauge   pic 1re onde %7.2f MPa a %5.1f us  [onde %5.1f-%5.1f us%s]"
                 % (g["sig"] / 1e6, _us(g["t"]), _us(g["t_on"]), _us(g["t_off"]),
                    "" if g["finished"] else ", NON retombee a la fin"))
        L.append("jauge   max global   %7.2f MPa a %5.1f us  %s"
                 % (g["sig_max"] / 1e6, _us(g["t_max"]),
                    "(= pic de la 1re onde)" if g["t_max"] == g["t"]
                    else "(onde ULTERIEURE : ce n'est pas la reference de Yang)"))
    elif g is not None:
        L.append("jauge   n.m. (%s)" % g["reason"])
    if yang:
        L.append("Yang 9 m/s : v_ind %.2f m/s (pente, par. 4.1), v_reb %.2f m/s (pente apres "
                 "450 us), jauge ~%.0f MPa (pic de la 1re onde a mi-bit)"
                 % (yang.get("vind", NAN), yang.get("vreb", NAN), yang.get("sig", NAN)))
    return "\n".join(L)


def guess_log(run):
    """Journal probable d'un run rockim : out_<tag> -> results/<tag>.log."""
    import os
    base = os.path.basename(run.rstrip("/\\"))
    if base.startswith("out_"):
        base = base[4:]
    return os.path.join("results", base + ".log")


if __name__ == "__main__":
    import os
    import sys
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    flags = [a for a in sys.argv[1:] if a.startswith("--")]
    robust = "--robust" in flags
    if not args:
        print("usage: yang_estimators.py [--robust] out_dir [out_dir ...]\n"
              "  --robust : ajoute le controle croise (pente de Theil-Sen, "
              "sensibilite a la fenetre, transits 1D de la premiere onde)")
        sys.exit(1)
    for run in args:
        h = load_history(run)
        res = kinetics(h)
        print(format_report(res, label=run,
                            yang=dict(vind=5.62, vreb=4.65, sig=160.0)))
        if robust:
            log = guess_log(run)
            tr, src = transit_from_log(log)
            print("(geometrie lue dans %s : %s)"
                  % (log if os.path.isfile(log) else log + " ABSENT -> defauts",
                     ", ".join("%s=%s" % (k, v) for k, v in sorted(src.items()))
                     or "aucune"))
            print(format_robustness(res, transit=tr))
