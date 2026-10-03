#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# lefm_aile.py : modeles MLER (LEFM) d'amorcage d'une fissure en aile a partir
# d'une fissure frottante inclinee sous compression (Lisjak 2013, these, §3.4.3,
# Tables 3.2 et 3.3). Fichier de calcul pur (pas de JOBS) : importe par
# p3_aile_lisjak.py, executable seul pour regenerer lisjak_reference.json.
#
# Convention : compression POSITIVE. Fissure de longueur 2c, inclinee de gamma
# sur sigma1, frottement mu ; aile droite de longueur l, inclinee de theta sur
# le plan de la fissure (angle de branchement) ; sigma2 lateral.
#   sigma_n = (s1 + s2)/2 - (s1 - s2)/2 cos(2 gamma)
#   tau     = (s1 - s2)/2 sin(2 gamma)
#   tau*    = tau - mu sigma_n           (cisaillement moteur, fissure fermee)
#   F       = 2 c tau*                   (force de coin equivalente, par m)
#   sigma_N = s1 sin^2(theta - gamma) + s2 cos^2(theta - gamma)
#             (contrainte normale distante sur le plan de l'aile)
# Amorcage : K_I(sigma1) = K_Ic (eq. 3.22), K_I lineaire en sigma1.
#
# Statut des formules (Lisjak ne les donne pas ; references originales absentes
# du depot) :
#   cotterell_rice : VERIFIEE, reproduit 5,1 MPa (5,07) ;
#   horii_nemat_nasser : VERIFIEE, reproduit 14,6 MPa (14,57) ;
#   zaitsev, lehner_kachanov : RECONSTRUCTIONS CANDIDATES, NON VERIFIEES sur
#     les textes originaux (voir README) ; servent seulement d'ordre de grandeur.
#     Le critere P3.3 utilise les valeurs PUBLIEES de la Table 3.3.
# ---------------------------------------------------------------------------
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))

# Table 3.1 et 3.2 (these, pp. 37 et 39 imprimees)
P = dict(E=3e9, nu=0.29, rho=2300.0, GIc=2.0, GIIc=10.0, ft=3e6, c_coh=15e6,
         phi_i=35.0, phi_f=35.0, mu_damp=7.4e3, pn=30e9, pt=3e9, pf=15e9,
         k_platen=0.1, c=2.5e-3, gamma=45.0, mu=0.7, theta=64.0, l=0.7e-3,
         s2=0.0)


def stresses(s1, s2, gamma):
    g = math.radians(gamma)
    sn = 0.5 * (s1 + s2) - 0.5 * (s1 - s2) * math.cos(2 * g)
    tau = 0.5 * (s1 - s2) * math.sin(2 * g)
    return sn, tau


def tau_star(s1, p):
    sn, tau = stresses(s1, p["s2"], p["gamma"])
    return tau - p["mu"] * sn


def sigma_N(s1, p):
    a = math.radians(p["theta"] - p["gamma"])
    return s1 * math.sin(a) ** 2 + p["s2"] * math.cos(a) ** 2


def kI_cotterell_rice(s1, p):
    """branchement infinitesimal (Cotterell et Rice 1980, premier ordre) :
    K_II = tau* sqrt(pi c), K_I = 0 (fissure fermee) ;
    k_I(theta) = 3/4 [sin(theta/2) + sin(3 theta/2)] K_II"""
    th = math.radians(p["theta"])
    KII = tau_star(s1, p) * math.sqrt(math.pi * p["c"])
    return 0.75 * (math.sin(th / 2) + math.sin(1.5 * th)) * KII


def kI_horii_nemat_nasser(s1, p):
    """aile droite chargee par la force de coin F sin(theta) en son centre,
    longueur effective l + l*, l* = 0,27 c (Horii et Nemat-Nasser 1986) :
    K_I = F sin(theta) / sqrt(pi (l + l*)) - sigma_N sqrt(pi l)"""
    th = math.radians(p["theta"])
    F = 2 * p["c"] * tau_star(s1, p)
    return (F * math.sin(th) / math.sqrt(math.pi * (p["l"] + 0.27 * p["c"]))
            - sigma_N(s1, p) * math.sqrt(math.pi * p["l"]))


def kI_zaitsev_candidat(s1, p):
    """CANDIDAT non verifie : coin de glissement moyen b = pi tau* c / E'
    ouvrant l'aile a sa racine, K_I = E' b sin(theta) / (2 sqrt(pi l))
    = pi tau* c sin(theta) / (2 sqrt(pi l)). Donne 6,82 MPa contre 6,9 publie."""
    th = math.radians(p["theta"])
    return math.pi * tau_star(s1, p) * p["c"] * math.sin(th) / (2 * math.sqrt(math.pi * p["l"]))


def kI_lehner_kachanov_candidat(s1, p):
    """CANDIDAT non verifie (ajustement de forme, pas une derivation) :
    K_I = F / sqrt(pi l) - sigma_N sqrt(pi (l + 0,27 c)). Donne 8,53 MPa
    contre 8,5 publie ; la coincidence ne prouve pas la formule."""
    F = 2 * p["c"] * tau_star(s1, p)
    return (F / math.sqrt(math.pi * p["l"])
            - sigma_N(s1, p) * math.sqrt(math.pi * (p["l"] + 0.27 * p["c"])))


MODELS = {
    "cotterell_rice_1980": (kI_cotterell_rice, 5.1, "verifiee"),
    "zaitsev_1985": (kI_zaitsev_candidat, 6.9, "candidate non verifiee"),
    "horii_nemat_nasser_1986": (kI_horii_nemat_nasser, 14.6, "verifiee"),
    "lehner_kachanov_1996": (kI_lehner_kachanov_candidat, 8.5, "candidate non verifiee"),
}


def KIc(p, which="E"):
    """eq. 3.23 : G_Ic = K_Ic^2 / E (Lisjak) ; variante E' = E/(1 - nu^2)"""
    E = p["E"] if which == "E" else p["E"] / (1 - p["nu"] ** 2)
    return math.sqrt(E * p["GIc"])


def sigma1c(model, p, K):
    f = MODELS[model][0]
    return K / f(1.0, p)          # K_I lineaire en sigma1 (s2 = 0)


def theta_max_hoop():
    """angle de branchement de la contrainte circonferentielle maximale
    (Erdogan et Sih) en mode II pur : K_II (3 cos(theta) - 1) = 0, soit
    theta = arccos(1/3) = 70,53 deg"""
    return math.degrees(math.acos(1.0 / 3.0))


def reference():
    p = dict(P)
    K_E, K_Ep = KIc(p, "E"), KIc(p, "Ep")
    out = {}
    for m, (f, pub, stat) in MODELS.items():
        out[m] = dict(publie_MPa=pub, statut_formule=stat,
                      recalcule_KIc_0077_MPa=round(sigma1c(m, p, 0.077e6) / 1e6, 3),
                      recalcule_KIc_E_MPa=round(sigma1c(m, p, K_E) / 1e6, 3),
                      recalcule_KIc_Eprime_MPa=round(sigma1c(m, p, K_Ep) / 1e6, 3),
                      formule=f.__doc__.strip().replace("\n    ", " "))
    mu_c = 2 * 0.7e-3 * math.sqrt(p["rho"] * p["E"])
    return dict(
        source="Lisjak Bradley A. (2013), PhD thesis, University of Toronto, chap. 3, §3.4.3 "
               "« Seismic analysis of AE », pp. 36-39 imprimees (PDF + 20)",
        geometrie=dict(L_x_mm=50, L_y_mm=100, fissure_mm=5, inclinaison_deg=45,
                       position="centre", maillage="Delaunay h = 0,70 mm, 21 600 triangles",
                       dt_s=5e-9, platines="deux platines rigides a 0,05 m/s en sens opposes (p. 36)",
                       page=36),
        table_3_1=dict(page=37, rho=2300, E_GPa=3, nu=0.29, phi_i_deg=35, cohesion_MPa=15,
                       phi_f_deg=35, ft_MPa=3, GIc_J_m2=2.0, GIIc_J_m2=10,
                       mu_amort_kg_m_s=7.4e3, pn_GPa_m=30, pt_GPa_par_m=3, pf_GPa=15,
                       frottement_platine=0.1),
        amortissement=dict(page=36, eq="3.21 mu_c = 2 h sqrt(rho E)", h_m=0.7e-3,
                           mu_c_recalcule=round(mu_c, 1), mu_c_publie=3.7e3,
                           mu_utilise="2 mu_c = 7,4e3 kg/(m s)"),
        table_3_2=dict(page=39, deux_c_mm=5, gamma_deg=45, mu=0.7, theta_deg=64, l_mm=0.7,
                       sigma2_MPa=0),
        table_3_3=dict(page=39, cotterell_rice_1980=5.1, zaitsev_1985=6.9,
                       horii_nemat_nasser_1986=14.6, lehner_kachanov_1996=8.5, ygeo_fdem=7.0),
        resultats_ygeo=dict(sigma1c_MPa=7.0, page=37, theta_branchement_deg=64,
                            sigma1u_MPa=[33, 32], note_sigma1u="33 MPa dans le texte p. 37, 32 MPa "
                            "dans la legende de la fig. 3.4c", ratio=0.21,
                            t_evenements_1_2_ms=[[1.86, 1.87], [1.89, 1.90]],
                            t_rupture_ms=9, detection="sigma1c lu sur les reactions nodales des "
                            "platines a l'amorcage des ailes (evenements EA 1 et 2, rupture "
                            "d'element de fissure)"),
        KIc=dict(eq="3.23 G_Ic = K_Ic^2 / E", publie_MPa_sqrt_m=0.077,
                 recalcule_E=round(K_E / 1e6, 5), recalcule_Eprime=round(K_Ep / 1e6, 5)),
        theta_contrainte_circonferentielle_max_deg=round(theta_max_hoop(), 2),
        lefm=out)


if __name__ == "__main__":
    ref = reference()
    with open(os.path.join(HERE, "lisjak_reference.json"), "w") as f:
        json.dump(ref, f, indent=1, ensure_ascii=False)
    for m, r in ref["lefm"].items():
        print(f"{m:26s} publie {r['publie_MPa']:5.1f}  K=0,077 {r['recalcule_KIc_0077_MPa']:6.2f}  recalcule (E) {r['recalcule_KIc_E_MPa']:6.2f}"
              f"  (E') {r['recalcule_KIc_Eprime_MPa']:6.2f}  [{r['statut_formule']}]")
    print("K_Ic", ref["KIc"], "theta_max", ref["theta_contrainte_circonferentielle_max_deg"])
