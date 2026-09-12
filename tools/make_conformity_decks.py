#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# make_conformity_decks.py — decks de conformite de la campagne du 13/09
# (T4, docs/CAMPAGNE_correction_2026-09-13.md). Genere, SANS lancer :
#
#   (1) configs/stanne2025_bench_s25_visc.cfg
#       = configs/stanne2025_bench_s25_visc0.cfg (ecrit a la main, la source)
#         + bulkViscosity = 2000 (exploratoire, ARMA 24-0952 « mass damping
#         4000 » lu comme eta D de Guo eq. 2.6 = 2 mu D, convention du banc C).
#
#   (2) configs/yang2026_bench_s25_v4_<X>.cfg, X dans A B B1 B2 C D
#       = le deck v3 du banc X (corps IDENTIQUE, octet pour octet) + un en-tete
#         + les trois cles de conformite en fin de deck :
#           meanTensionCapFactor = 0     (cap cache 3 ft neutralise, DIAGNOSTIC §4)
#           jointBreakModeRef = slipRef  (etiquette rn/rs de S1, sortie seule)
#           writeRuptureFields = true    (dead / openMax / pMean, S1(a), sortie
#                                         seule, necessaires a T1 ; history.csv
#                                         intact, contrairement a contactForcePairs
#                                         qui reste hors de la serie v4)
#       L en-tete porte la liste EXACTE des cles qui different du temoin A
#       (configs/yang2026_bench_s25_v3P_300.cfg), calculee ici par diff des
#       paires cle = valeur (commentaires ignores), pas recopiee a la main.
#
#   python tools/make_conformity_decks.py            # ecrit les 7 decks
#   python tools/make_conformity_decks.py --check    # ne fait que verifier
#                                                    # (code 1 si un deck differe)
#
# Regles : aucun run, aucune cle nouvelle dans le solveur ; les decks sources
# ne sont jamais modifies ; un deck source qui porterait deja l une des deux
# cles est REFUSE (le doublon serait silencieux dans Config).
# ---------------------------------------------------------------------------
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CFG = os.path.join(ROOT, "configs")

KEY_RE = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_.]*)\s*=\s*(.*?)\s*(#.*)?$")

# banc -> (deck source v3, description courte)
SERIE_V4 = [
    ("A",  "yang2026_bench_s25_v3P_300.cfg",
     "le TEMOIN A (v3-P : plastic + coulomb + ratchet, 300 us)"),
    ("B",  "yang2026_bench_s25_solidity.cfg",
     "le banc B (loi solidity + penalite edge 25 + contact de Solidity)"),
    ("B1", "yang2026_bench_s25_law.cfg",
     "le banc B1 (temoin + la SEULE loi jointShearUnload = solidity)"),
    ("B2", "yang2026_bench_s25_pen.cfg",
     "le banc B2 (temoin + la SEULE penalite edge / 25)"),
    ("C",  "yang2026_bench_s25_solidity_visc.cfg",
     "le banc C (banc B + bulkViscosity = 2000)"),
    ("D",  "yang2026_bench_s25_v3P_contact.cfg",
     "le banc D (temoin + le SEUL contact de Solidity : potPenaltyFactor 0.25, gcBirth penalty)"),
]
TEMOIN = "yang2026_bench_s25_v3P_300.cfg"
# Les cles de conformite ajoutees en fin de chaque deck v4. TOUTES sont soit une
# neutralisation d un garde-fou maison (meanTensionCapFactor), soit une SORTIE
# (jointBreakModeRef n est qu une etiquette, writeRuptureFields n ajoute que des
# champs VTU) : aucune ne change une force.
# `contactForcePairs` (S2) n est deliberement PAS dans cette liste : elle ajoute
# des colonnes a history.csv, ce qui casserait le controle « history.csv de la v4
# identique a celui de la v3 » qui est justement le test de bit-identite de S1 sur
# un vrai deck (critere (b) de l en-tete). Les decks St Anne, qui n ont pas de
# jumeau v3 a comparer, la portent (leur critere premier EST la reaction).
V4_KEYS = [("meanTensionCapFactor", "0"),
           ("jointBreakModeRef", "slipRef"),
           ("writeRuptureFields", "true")]

STANNE_SRC = "stanne2025_bench_s25_visc0.cfg"
STANNE_DST = "stanne2025_bench_s25_visc.cfg"


def read(path):
    with open(path, encoding="utf-8", newline="") as f:
        return f.read()


def keys_of(text):
    """Paires cle -> valeur (derniere occurrence), commentaires ignores."""
    out = {}
    for line in text.splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        m = KEY_RE.match(line)
        if m:
            out[m.group(1)] = m.group(2).strip()
    return out


def diff_keys(src, ref):
    """Lignes 'cle  valeur_src  (temoin : valeur_ref | absente)' pour les cles
    de src qui different du temoin, puis les cles du temoin absentes de src."""
    rows = []
    for k, v in src.items():
        if k not in ref:
            rows.append("%-24s %-22s (absente du temoin)" % (k, v))
        elif ref[k] != v:
            rows.append("%-24s %-22s (temoin : %s)" % (k, v, ref[k]))
    for k, v in ref.items():
        if k not in src:
            rows.append("%-24s %-22s (RETIREE ; temoin : %s)" % (k, "-", v))
    return rows


def v4_text(banc, src_name, desc, src_text, ref_keys):
    src_keys = keys_of(src_text)
    for k, _ in V4_KEYS:
        if k in src_keys:
            raise SystemExit("REFUS : %s porte deja la cle %s" % (src_name, k))
    rows = diff_keys(src_keys, ref_keys)
    rows += ["%-24s %-22s (AJOUTEE ici)" % (k, v) for k, v in V4_KEYS]
    dst_name = "yang2026_bench_s25_v4_%s.cfg" % banc
    head = [
        "# ---------------------------------------------------------------------------",
        "# %s — SERIE v4 (T4, campagne de correction du 13/09" % dst_name,
        "# apres-midi, docs/CAMPAGNE_correction_2026-09-13.md) : %s" % desc,
        "# + les TROIS cles de conformite, ajoutees en FIN de deck :",
        "#   meanTensionCapFactor = 0     cap cache sur la pression moyenne en traction",
        "#                                (defaut 3 ft = 32,9 MPa, actif dans TOUS les runs",
        "#                                v3, DIAGNOSTIC §4 / ECARTS §5) neutralise ;",
        "#   jointBreakModeRef = slipRef  etiquette rn/rs a la rupture normalisee par la",
        "#                                plage du moteur (S1, SORTIE SEULE, aucune force",
        "#                                changee) — EXIGE rockim_g1y17.exe : sous g1y16",
        "#                                la cle est inconnue et le deck est REFUSE",
        "#                                (unknownKeys = error), c est voulu ;",
        "#   writeRuptureFields = true    champs VTU `dead` (0/1) et `openMax` (m) sur les",
        "#                                joints, `pMean` (Pa) sur les elements (S1(a),",
        "#                                SORTIE SEULE) — necessaires a T1",
        "#                                (tools/crack_paths.py) pour separer une facette",
        "#                                rompue mais FERMEE d une fissure ouverte ; ne",
        "#                                touche PAS history.csv.",
        "# PAS de `contactForcePairs` (S2) ici : elle ajouterait des colonnes a history.csv",
        "# et casserait le controle (b) ci-dessous. Les decks St Anne la portent.",
        "# Source : configs/%s (corps recopie octet pour octet par" % src_name,
        "# tools/make_conformity_decks.py ; meme maillage s = 2,5 historique que la serie",
        "# v3 pour que la SEULE difference v3 -> v4 soit ces deux cles).",
        "#",
        "# CLES QUI DIFFERENT DU TEMOIN A (configs/%s), liste EXACTE" % TEMOIN,
        "# calculee par diff des paires cle = valeur (les commentaires ne comptent pas) :",
    ]
    if rows:
        head += ["#   " + r for r in rows]
    else:
        head += ["#   (aucune)"]
    head += [
        "#",
        "# COUT ESTIME : docs/DECKS_conformite_2026-09-13.md (dt mesure par le lancement",
        "# de 2 us, tools/deck_smoke.py ; 77,7 ms/pas machine partagee 14 fils d apres",
        "# results/yang_bench_s25_v3P.log, ~17 ms/pas seul).",
        "# CRITERE DE SUCCES : celui du banc source (temoin A : plateau de reaction",
        "# >= 45 kN avant p = 0,5 mm, v_bit max <= 6,3 m/s, nPulv > 0, traces de joints",
        "# au-dela de r = 9 mm) ; EN PLUS (a) le cap est ETEINT ici, donc le compteur",
        "# d ecretage de S1 n a rien a compter : savoir si le cap mordait dans la v3",
        "# demande de rejouer un deck v3 (cap = 3) AVEC writeRuptureFields = true, ce que",
        "# la v4 ne fait pas — la v4 mesure l effet du cap par la DIFFERENCE v3/v4, pas",
        "# par un compteur ; (b) proportions traction/cisaillement a la rupture a comparer",
        "# a la v3 : l etiquette seule change, pas les forces — `history.csv` doit etre",
        "# identique a celui de la v3 hors les colonnes de recensement nBrokTen/nBrokShear",
        "# si le cap n a jamais mordu (c est le test de bit-identite de S1 sur un vrai",
        "# deck ; les VTU, eux, portent les champs de S1(a) en plus).",
        "# ---------------------------------------------------------------------------",
    ]
    tail = [
        "",
        "# ---------------------------------------------------------------------------",
        "# SERIE v4 (T4, 13/09 apres-midi) — les trois cles de conformite",
        "# ---------------------------------------------------------------------------",
    ] + ["%s = %s" % (k, v) for k, v in V4_KEYS]
    body = src_text if src_text.endswith("\n") else src_text + "\n"
    return dst_name, "\n".join(head) + "\n" + body + "\n".join(tail) + "\n"


def stanne_visc_text(src_text):
    keys = keys_of(src_text)
    if "bulkViscosity" in keys:
        raise SystemExit("REFUS : %s porte deja bulkViscosity" % STANNE_SRC)
    old_head = "# stanne2025_bench_s25_visc0.cfg — CAS DE CONFORMITE ST ANNE"
    if old_head not in src_text:
        raise SystemExit("en-tete de %s non reconnu" % STANNE_SRC)
    t = src_text.replace(old_head,
                         "# stanne2025_bench_s25_visc.cfg — CAS DE CONFORMITE ST ANNE", 1)
    t = t.replace("# amortissement. Impact a insert unique",
                  "# amortissement REMPLACEE par : variante AVEC bulkViscosity = 2000\n"
                  "# (EXPLORATOIRE, genere par tools/make_conformity_decks.py depuis\n"
                  "# _visc0 : seule la cle bulkViscosity est ajoutee, en fin de section\n"
                  "# AMORTISSEMENT). Impact a insert unique", 1)
    marker = "dampingLocal = 0\n"
    if t.count(marker) != 1:
        raise SystemExit("marqueur dampingLocal introuvable ou multiple")
    add = (marker +
           "# EXPLORATOIRE (T4) : ARMA 24-0952 « Mass Damping Coefficient 4000 » lu comme\n"
           "# eta = 4000 Pa.s de Guo 2014 eq. 2.6 (T = ... + eta D), soit 2000 dans la\n"
           "# convention 2 mu D de rockim — la convention du banc C. Ni l unite ni\n"
           "# l operateur ne sont publies (COMPLEMENT §3) : un amortissement de masse\n"
           "# -alpha M v freinerait aussi la translation rigide, ce qu une viscosite de\n"
           "# volume ne fait pas. Viscosite GLOBALE (acier et carbure compris).\n"
           "bulkViscosity = 2000\n")
    t = t.replace(marker, add, 1)
    # la ligne de lancement de l en-tete
    t = t.replace("stanne2025_bench_s25_visc0.cfg out_stanne_s25_visc0",
                  "stanne2025_bench_s25_visc.cfg out_stanne_s25_visc", 1)
    return t


def main():
    check = "--check" in sys.argv
    ref_keys = keys_of(read(os.path.join(CFG, TEMOIN)))
    todo = []
    for banc, src_name, desc in SERIE_V4:
        src_text = read(os.path.join(CFG, src_name))
        todo.append(v4_text(banc, src_name, desc, src_text, ref_keys))
    todo.append((STANNE_DST, stanne_visc_text(read(os.path.join(CFG, STANNE_SRC)))))
    bad = 0
    for name, text in todo:
        path = os.path.join(CFG, name)
        if any(ord(c) < 32 and c not in "\t\n\r" for c in text):
            raise SystemExit("caractere de controle dans %s" % name)
        if check:
            same = os.path.exists(path) and read(path) == text
            print("%-40s %s" % (name, "identique" if same else "DIFFERE / absent"))
            bad += 0 if same else 1
        else:
            with open(path, "w", encoding="utf-8", newline="\n") as f:
                f.write(text)
            print("ecrit %s (%d lignes)" % (name, text.count("\n")))
    if check and bad:
        sys.exit(1)


if __name__ == "__main__":
    main()
