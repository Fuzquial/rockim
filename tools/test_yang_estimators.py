#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ---------------------------------------------------------------------------
# test_yang_estimators.py — mini-test FALSIFIANT de tools/yang_estimators.py
# (campagne du 13/09, tache T2). Aucun solveur : signaux synthetiques a pentes
# connues, puis les deux runs reels si presents.
#
#   python tools/test_yang_estimators.py [out_yang2026_v3] [out_yang_bench_s25_v3P]
#
# Ce qui DOIT tenir (a 1e-9 relatif) :
#   A. deplacement p(t) = rampe quadratique 0-20 us, puis lineaire a 6 m/s,
#      retournement parabolique 200-220 us, remontee lineaire a 4 m/s :
#      pente 10-90 % = 6 m/s EXACTEMENT ; pente de rebond = 4 m/s EXACTEMENT ;
#      fenetres entierement dans les portions lineaires.
#   B. vitesse instantanee = dp/dt + pic gaussien de +1,4 m/s a 100 us :
#      max instantane = 7,4 m/s (et non 6) -> les deux estimateurs
#      DIFFERENT sur le meme signal ; c'est l'ecart 6,86 / 7,37 du s = 1.
#   C. jauge = demi-sinus 160 MPa (20-120 us) puis demi-sinus 200 MPa
#      (300-400 us) : pic de la 1re onde = 160 ; max global = 200.
#      Variante naive (max global) : 200 != 160 -> DOIT ECHOUER.
#   D. enregistrement tronque a 210 us (avant le retournement a 212 us) et a
#      216 us (retournement inacheve, 0,3 % de p_max recuperes) : le rebond
#      DOIT etre declare non mesurable ; tronque a 300 us : 4 m/s exactement.
#   E. runs reels : les chiffres publies du s = 1 (6,86 / 7,37 / 175,95 MPa).
# Ajouts de la deuxieme passe (13/09 au soir) — le chiffre de pente est-il une
# mesure ou un artefact d'ajustement ?
#   F. pour chaque run donne : F1 la pente de Theil-Sen (mediane des pentes de
#      paires, estimateur independant) egale les moindres carres a 3 % sur la
#      MEME fenetre ; F2 la pente varie de moins de 6 % entre les fenetres
#      larges 5-95 / 10-90 / 20-80 ; F3 et F4 les deux estimateurs naifs
#      (vitesse moyenne p_max/t_pmax, max instantane) DOIVENT differer de plus
#      de 15 % et 10 % (sinon la definition de Yang ne changerait rien) ;
#      F5 le pic retenu tombe entre l'arrivee de l'onde a la bande de jauge et
#      le retour de la reflexion du bas du bit (transits 1D lus dans le
#      journal du run) ; F6 ce pic vaut l'impact 1D rho c v / 2 a 10 % pres.
#   G. les criteres F2 et F5 PEUVENT echouer : un deplacement bilineaire
#      (3 puis 9 m/s) viole F2, une jauge muette avant 80 us viole F5, et les
#      transits du train sont verifies contre le calcul a la main.
# Sortie : PASS/FAIL par critere et code retour 0/1.
# ---------------------------------------------------------------------------
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import yang_estimators as ye  # noqa: E402

TOL = 1e-9
fails = []


def check(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    print("  [%s] %s %s" % (tag, name, detail))
    if not cond:
        fails.append(name)


def synthetic(t_end=600e-6, dt=1e-7):
    """p(t) en m, v(t) = dp/dt en m/s, C1 par morceaux, v_lin = 6, v_reb = 4."""
    t = np.arange(0.0, t_end + 0.5 * dt, dt)
    t1, t2, t3 = 20e-6, 200e-6, 220e-6
    v0 = 6.0
    a1 = v0 / t1                     # rampe quadratique 0..t1
    a2 = (v0 + 4.0) / (t3 - t2)      # deceleration constante t2..t3 (6 -> -4)
    p = np.empty_like(t)
    v = np.empty_like(t)
    m = t < t1
    p[m] = 0.5 * a1 * t[m] ** 2
    v[m] = a1 * t[m]
    p1 = 0.5 * a1 * t1 ** 2
    m = (t >= t1) & (t < t2)
    p[m] = p1 + v0 * (t[m] - t1)
    v[m] = v0
    p2 = p1 + v0 * (t2 - t1)
    m = (t >= t2) & (t < t3)
    tau = t[m] - t2
    p[m] = p2 + v0 * tau - 0.5 * a2 * tau ** 2
    v[m] = v0 - a2 * tau
    p3 = p2 + v0 * (t3 - t2) - 0.5 * a2 * (t3 - t2) ** 2
    m = t >= t3
    p[m] = p3 - 4.0 * (t[m] - t3)
    v[m] = -4.0
    t_turn = t2 + v0 / a2
    p_max = p2 + 0.5 * v0 ** 2 / a2
    return t, p, v, dict(t1=t1, t2=t2, t3=t3, t_turn=t_turn, p_max=p_max, p1=p1, p3=p3)


def main():
    print("=== T2 : test falsifiant des estimateurs de Yang ===")
    t, p, v, g = synthetic()
    print("synthetique : p_max %.4f mm a %.1f us ; lineaire 6 m/s sur %.0f-%.0f us ; "
          "rebond -4 m/s apres %.0f us" % (g["p_max"] * 1e3, g["t_turn"] * 1e6,
                                          g["t1"] * 1e6, g["t2"] * 1e6, g["t3"] * 1e6))

    # ---- A : pentes exactes ------------------------------------------------
    ind = ye.indentation_slope(t, p)
    check("A1 pente 10-90 % = 6 m/s a 1e-9", ind["ok"] and abs(ind["v"] - 6.0) < TOL * 6.0,
          "mesure %.12f m/s, fenetre %.1f-%.1f us, r2 %.12f" % (
              ind["v"], ind["t0"] * 1e6, ind["t1"] * 1e6, ind["r2"]))
    check("A2 fenetre dans la portion lineaire", ind["ok"] and ind["t0"] >= g["t1"] - 1e-12
          and ind["t1"] <= g["t2"] + 1e-12,
          "[%.1f, %.1f] us dans [%.0f, %.0f]" % (ind["t0"] * 1e6, ind["t1"] * 1e6,
                                                 g["t1"] * 1e6, g["t2"] * 1e6))
    check("A3 p_max et t_turn retrouves", ind["ok"] and abs(ind["p_max"] - g["p_max"]) < 1e-12
          and abs(ind["t_pmax"] - g["t_turn"]) < 1.1e-7,
          "p_max %.6f mm a %.2f us (attendu %.6f mm a %.2f us)" % (
              ind["p_max"] * 1e3, ind["t_pmax"] * 1e6, g["p_max"] * 1e3, g["t_turn"] * 1e6))
    reb = ye.rebound_slope(t, p)
    check("A4 pente de rebond = 4 m/s a 1e-9", reb["ok"] and abs(reb["v"] - 4.0) < TOL * 4.0,
          "mesure %.12f m/s, fenetre %.1f-%.1f us, r2 %.12f" % (
              reb["v"], reb["t0"] * 1e6, reb["t1"] * 1e6, reb["r2"]))
    check("A5 fenetre de rebond apres t3", reb["ok"] and reb["t0"] >= g["t3"] - 1e-12,
          "t0 %.1f us >= %.0f us" % (reb["t0"] * 1e6, g["t3"] * 1e6))

    # ---- B : max instantane != pente ----------------------------------------
    vz = -(v + 1.4 * np.exp(-0.5 * ((t - 100e-6) / 5e-6) ** 2))   # rockim : vz < 0 vers le bas
    inst = ye.instantaneous(t, vz)
    check("B1 max instantane = 7,4 m/s", abs(inst["v_ind"] - 7.4) < 1e-9,
          "mesure %.12f m/s a %.1f us" % (inst["v_ind"], inst["t_ind"] * 1e6))
    check("B2 les deux estimateurs DIFFERENT sur le meme signal",
          abs(inst["v_ind"] - ind["v"]) > 1.0,
          "max %.3f contre pente %.3f m/s" % (inst["v_ind"], ind["v"]))
    check("B3 rebond instantane = 4 m/s (max de vz apres le min)",
          inst["ok_reb"] and abs(inst["v_reb"] - 4.0) < 1e-9,
          "mesure %.12f m/s" % inst["v_reb"])

    # ---- C : premiere onde vs max global -----------------------------------
    sig = np.zeros_like(t)
    m = (t >= 20e-6) & (t <= 120e-6)
    sig[m] = 160e6 * np.sin(np.pi * (t[m] - 20e-6) / 100e-6)
    m = (t >= 300e-6) & (t <= 400e-6)
    sig[m] = 200e6 * np.sin(np.pi * (t[m] - 300e-6) / 100e-6)
    fw = ye.first_wave_peak(t, sig)
    check("C1 pic de la 1re onde = 160 MPa", fw["ok"] and abs(fw["sig"] - 160e6) < 1e-9 * 160e6,
          "mesure %.6f MPa a %.1f us, onde %.1f-%.1f us, retombee %s" % (
              fw["sig"] / 1e6, fw["t"] * 1e6, fw["t_on"] * 1e6, fw["t_off"] * 1e6, fw["finished"]))
    check("C2 max global = 200 MPa (onde ulterieure)", abs(fw["sig_max"] - 200e6) < 1e-9 * 200e6
          and fw["t_max"] > 300e-6, "mesure %.6f MPa a %.1f us" % (fw["sig_max"] / 1e6, fw["t_max"] * 1e6))
    naive_ok = abs(fw["sig_max"] - 160e6) < 1e-9 * 160e6
    check("C3 variante NAIVE (max global) prise pour la reference : DOIT ECHOUER",
          not naive_ok, "max global %.0f != 160 MPa : l'ancien estimateur se trompe d'onde ici"
          % (fw["sig_max"] / 1e6))
    check("C4 1re onde bornee avant la 2e", fw["t_off"] < 300e-6,
          "t_off %.1f us < 300" % (fw["t_off"] * 1e6))

    # ---- D : enregistrements tronques ---------------------------------------
    for t_cut, expect_ok, lab in ((210e-6, False, "avant le retournement"),
                                  (216e-6, False, "retournement inacheve"),
                                  (300e-6, True, "en plein rebond lineaire")):
        m = t <= t_cut + 1e-12
        r = ye.rebound_slope(t[m], p[m])
        if expect_ok:
            check("D tronque a %.0f us (%s) : 4 m/s a 1e-9" % (t_cut * 1e6, lab),
                  r["ok"] and abs(r["v"] - 4.0) < TOL * 4.0,
                  "mesure %.12f m/s, %.1f %% de p_max recuperes" % (r["v"], 100 * r["amp"] / g["p_max"]))
        else:
            check("D tronque a %.0f us (%s) : NON mesurable" % (t_cut * 1e6, lab),
                  not r["ok"] and r["v"] != r["v"], "raison : %s" % r["reason"])
        i = ye.indentation_slope(t[m], p[m])
        check("D tronque a %.0f us : v_ind pente inchangee = 6 m/s" % (t_cut * 1e6),
              i["ok"] and abs(i["v"] - 6.0) < TOL * 6.0, "mesure %.12f m/s" % i["v"])

    # ---- E : runs reels (si presents) ---------------------------------------
    for run in sys.argv[1:]:
        if not os.path.isfile(os.path.join(run, "history.csv")):
            print("  (run %s absent : saute)" % run)
            continue
        h = ye.load_history(run)
        res = ye.kinetics(h)
        print(ye.format_report(res, label=run, yang=dict(vind=5.62, vreb=4.65, sig=160.0)))
        if os.path.basename(run.rstrip("/\\")) == "out_yang2026_v3":
            b = res["bodies"]
            check("E1 s = 1 : v_ind pente insert = 6,86 +- 0,01 m/s (ECARTS par. 5)",
                  abs(b["insert"]["ind"]["v"] - 6.86) < 0.01,
                  "mesure %.3f m/s [%.1f-%.1f us]" % (b["insert"]["ind"]["v"],
                                                       b["insert"]["ind"]["t0"] * 1e6,
                                                       b["insert"]["ind"]["t1"] * 1e6))
            check("E2 s = 1 : v_ind max instantane bit = 7,37 +- 0,01 m/s",
                  abs(b["bit"]["inst"]["v_ind"] - 7.37) < 0.01,
                  "mesure %.3f m/s" % b["bit"]["inst"]["v_ind"])
            check("E3 s = 1 : rebond non mesurable (run arrete a 183 us, bit descendant)",
                  not b["bit"]["reb"]["ok"] and not b["bit"]["inst"]["ok_reb"],
                  "raison : %s" % b["bit"]["reb"]["reason"])
            check("E4 s = 1 : pic 1re onde = max global = 175,95 MPa",
                  abs(res["gauge"]["sig"] / 1e6 - 175.95) < 0.01 and res["gauge"]["t"] == res["gauge"]["t_max"],
                  "mesure %.2f MPa a %.1f us" % (res["gauge"]["sig"] / 1e6, res["gauge"]["t"] * 1e6))

        # ---- F : controle croise (T2 deuxieme passe, 13/09 soir) ------------
        # La pente 10-90 % est-elle une MESURE ou un artefact de fenetre et
        # d'ajustement ? Trois questions, chacune avec un chiffre :
        #   F1 un estimateur independant (Theil-Sen, mediane des pentes de
        #      paires) sur la MEME fenetre doit donner la meme pente ;
        #   F2 la pente doit peu bouger entre les fenetres LARGES (5-95,
        #      10-90, 20-80) : sinon p(t) n'a pas de portion lineaire ;
        #   F3-F4 les deux estimateurs naifs (vitesse moyenne p_max/t_pmax et
        #      max instantane) doivent DIFFERER franchement : c'est la raison
        #      d'etre de la definition de Yang ;
        #   F5-F6 le pic retenu doit tomber entre l'arrivee de l'onde a la
        #      jauge et le retour de la reflexion du bas du bit (sinon ce
        #      n'est plus « la premiere onde »), et valoir l'impact 1D
        #      rho c v / 2 a 10 % pres.
        log = ye.guess_log(run)
        tr, src = ye.transit_from_log(log)
        print("  (transits : geometrie %s de %s : %s)"
              % ("LUE" if src and "erreur" not in src else "PAR DEFAUT", log,
                 ", ".join("%s=%s" % (k, v) for k, v in sorted(src.items())) or "aucune"))
        print(ye.format_robustness(res, transit=tr))
        for body, d in res["bodies"].items():
            ind = d["ind"]
            if not ind["ok"]:
                continue
            ts, npair = ye.slope_theilsen(res["t"][ind["i0"]:ind["i1"] + 1],
                                          d["p"][ind["i0"]:ind["i1"] + 1])
            dev = 100.0 * abs(ts - ind["v"]) / abs(ind["v"])
            check("F1 %s %s : Theil-Sen == moindres carres a 3 %%" % (run, body),
                  dev < 3.0, "MC %.4f, TS %.4f m/s (%d paires), ecart %.2f %%"
                  % (ind["v"], ts, npair, dev))
            ws = ye.window_sensitivity(res["t"], d["p"])
            check("F2 %s %s : dispersion des fenetres larges < 6 %%" % (run, body),
                  ws["spread_pct"] < 6.0,
                  "%.2f %% (5-95 %.4f, 10-90 %.4f, 20-80 %.4f m/s) ; toutes fenetres %.2f %%"
                  % (ws["spread_pct"], ws["rows"][0]["v"], ws["rows"][1]["v"],
                     ws["rows"][2]["v"], ws["spread_pct_wide"]))
            naive = ind["p_max"] / ind["t_pmax"]
            dn = 100.0 * abs(naive - ind["v"]) / abs(ind["v"])
            check("F3 %s %s : la vitesse moyenne p_max/t_pmax DOIT differer (> 15 %%)"
                  % (run, body), dn > 15.0,
                  "%.4f contre %.4f m/s, ecart %.1f %%" % (naive, ind["v"], dn))
            if d["inst"] is not None:
                di = 100.0 * abs(d["inst"]["v_ind"] - ind["v"]) / abs(ind["v"])
                check("F4 %s %s : le max instantane DOIT differer (> 10 %%)"
                      % (run, body), di > 10.0,
                      "%.4f contre %.4f m/s, ecart %.1f %%"
                      % (d["inst"]["v_ind"], ind["v"], di))
        g = res["gauge"]
        if g is not None and g["ok"]:
            check("F5 %s : le pic est entre l'arrivee et le retour de la reflexion" % run,
                  tr["t_arr_bar"] <= g["t"] <= tr["t_refl"],
                  "pic %.2f us dans [%.2f, %.2f] us (depart mesure %.2f us, premier "
                  "signal possible %.2f us a c_P)"
                  % (g["t"] * 1e6, tr["t_arr_bar"] * 1e6, tr["t_refl"] * 1e6,
                     g["t_on"] * 1e6, tr["t_first_p"] * 1e6))
            dev1d = 100.0 * abs(g["sig"] - tr["sig_1d"]) / tr["sig_1d"]
            check("F6 %s : pic = impact 1D rho c v/2 a 10 %% pres" % run, dev1d < 10.0,
                  "%.2f contre %.2f MPa, ecart %.1f %%"
                  % (g["sig"] / 1e6, tr["sig_1d"] / 1e6, dev1d))

    # ---- G : les deux criteres de la section F PEUVENT echouer --------------
    # Un test qui ne peut pas echouer ne mesure rien : on construit les deux
    # contre-exemples et on exige qu'ils VIOLENT F2 et F5.
    tg = np.arange(0.0, 200.1e-6, 1e-7)
    pg = np.where(tg <= 100e-6, 3.0 * tg, 3.0 * 100e-6 + 9.0 * (tg - 100e-6))
    wsg = ye.window_sensitivity(tg, pg)
    check("G1 deplacement BILINEAIRE (3 puis 9 m/s) : F2 DOIT echouer",
          wsg["spread_pct"] >= 6.0,
          "dispersion large %.2f %% (5-95 %.3f, 10-90 %.3f, 20-80 %.3f, 40-60 %.3f m/s)"
          % (wsg["spread_pct"], wsg["rows"][0]["v"], wsg["rows"][1]["v"],
             wsg["rows"][2]["v"], wsg["rows"][4]["v"]))
    trg = ye.wave_transit()
    # jauge qui ne voit RIEN avant 80 us puis une onde large : le pic retenu
    # (120 us) tombe apres le retour de la reflexion (74,25 us), donc ce n'est
    # plus « la premiere onde » mais une superposition -> F5 doit tomber.
    sg = np.zeros_like(tg)
    m = tg >= 80e-6
    sg[m] = 180e6 * np.sin(np.pi * (tg[m] - 80e-6) / 160e-6)
    fwg = ye.first_wave_peak(tg, sg)
    check("G2 jauge dont le pic tombe APRES le retour de reflexion : F5 DOIT echouer",
          not (trg["t_arr_bar"] <= fwg["t"] <= trg["t_refl"]),
          "pic retenu %.1f us hors de [%.1f, %.1f] us"
          % (fwg["t"] * 1e6, trg["t_arr_bar"] * 1e6, trg["t_refl"] * 1e6))
    check("G3 transits : c_barre, arrivee et retour de reflexion conformes a la "
          "geometrie du train", abs(trg["c_bar"] - 5047.5) < 0.5
          and abs(trg["t_arr_bar"] * 1e6 - 26.00) < 0.02
          and abs(trg["t_refl"] * 1e6 - 74.25) < 0.02
          and abs(trg["sig_1d"] / 1e6 - 178.30) < 0.02,
          "c_barre %.1f m/s, arrivee %.2f us, reflexion %.2f us, rho c v/2 %.2f MPa"
          % (trg["c_bar"], trg["t_arr_bar"] * 1e6, trg["t_refl"] * 1e6,
             trg["sig_1d"] / 1e6))

    print("=== %d echec(s) : %s ===" % (len(fails), ", ".join(fails) if fails else "aucun"))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
