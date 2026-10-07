#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""scan_decks.py — balayage STATIQUE des decks (*.cfg) selon la regle C1 des cles (w21, 2026-09-05 ; w22 : lecteurs).

Sans lancer rockim : on ne connait pas les cles CONSOMMEES, on applique la regle a partir du registre
(tools/keys_by_mode.json, genere par gen_keys_by_mode.py : `readers` = lecteurs de chaque cle) et de la table
des obsoletes (tools/obsolete_keys.json). Pour chaque cle d'un deck :
  * obsolete                      -> « cle obsolete X : remplacee par Z » (corrigee en place avec --fix-obsolete) ;
  * garde PRE-INIT de main.cpp    -> refusee avant meme l'audit : hydro* hors fdem (C6) ; thermal hors fdem ;
                                     beddingDip hors fdem ; law hors fem3d/fdem/fdem3d ; phases hors fdem/fdem3d ;
                                     mesh = voronoi hors fdem/fdem3d ; mesh = file hors fdem/fdem3d/fem3d ; mode inconnu
                                     (w22 : la relecture w21 avait releve tg_3d.cfg et uniax.cfg declares « sans faute ») ;
  * cle du registre dont le mode courant ou le code partage (`shared`) est LECTEUR -> legitime (une lecture
                                     conditionnelle n'est pas jugeable ici) ;
  * cle du registre dont AUCUN lecteur n'est le mode courant -> « cle X sans effet en mode Y : cles des modes Z1, Z2 »
                                     (w22, D1 : jointSoftening en fem3d, bond* hors dem, kpFactor hors fem...) ;
  * famille dynamique (phase.<nom>.E, groupBond.<A>.<B>...) -> non jugeable statiquement, comptee a part ;
  * sinon                         -> faute de frappe (Levenshtein <= 2 sur le registre, suggestion) ou inconnue de rockim.
Un fichier SANS cle de solveur (ni mode, mesh, scenario, W, H, T, nx, ny, frames : deck de point materiel
`rockim matpoint <cfg>`, carte materiau `law = cdp` + constantes, deck temporaire de selftest) n'est pas un deck
de solveur : il est compte a part et non juge (il ne passe jamais par l'audit de main.cpp). Ce que le juge statique ne voit PAS :
les gardes de maillage (C3 : noeud orphelin, tet plat — il faut le .msh), les lectures conditionnelles, les
familles dynamiques ; le verdict definitif est celui du binaire.
Le rapport (markdown) liste TOUTES les fautes ; seules les obsoletes sont corrigees (--fix-obsolete) : le reste
est pour decision de Fernando.

usage : python tools/scan_decks.py [--roots D1 D2 ...] [--fix-obsolete] [--out tools/scan_decks_<date>.md]
        (defaut : racine rockim_f2 + CONTINUUM/calib_bohus_triax/cdp_rockim ; dossiers out_*, orig/, selftest_* ignores)
"""
import argparse, io, json, os, re, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
CDP = os.path.normpath(os.path.join(ROOT, "..", "..", "CONTINUUM", "calib_bohus_triax", "cdp_rockim"))
SKIP_DIRS = re.compile(r"^(out_|__pycache__|\.git|orig$|selftest_|adverse$)")   # orig = copies avant edition ; selftest_* / adverse = bancs falsifiants (fautes VOULUES)


def levenshtein(a, b, maxd=2):
    if abs(len(a) - len(b)) > maxd:
        return maxd + 1
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        if min(cur) > maxd:
            return maxd + 1
        prev = cur
    return prev[-1]


def read_deck(path):
    """[(cle, valeur, ligne)] dans l'ordre du fichier (derniere definition gardee, comme Config::load)."""
    kv = {}
    with io.open(path, encoding="utf-8", errors="ignore") as f:
        for n, line in enumerate(f, 1):
            line = line.split("#", 1)[0].strip()
            if not line or "=" not in line:
                continue
            k, v = line.split("=", 1)
            kv[k.strip()] = (v.strip(), n)
    return kv


def load_registry():
    with io.open(os.path.join(HERE, "keys_by_mode.json"), encoding="utf-8") as f:
        reg = json.load(f)
    mode_of = {}
    for m, ks in reg["modes"].items():
        for k in ks:
            mode_of[k] = m
    known = set(reg["common"]) | set(mode_of)
    readers = {k: set(v) for k, v in reg.get("readers", {}).items()}
    known |= set(readers)
    dyn = tuple(reg.get("prefixes_dynamic", []))
    obs = {}
    p = os.path.join(HERE, "obsolete_keys.json")
    if os.path.isfile(p):
        with io.open(p, encoding="utf-8") as f:
            obs = {k: v["new"] for k, v in json.load(f)["obsolete"].items()}
    return reg, mode_of, known, dyn, obs, readers


MODES = ("fem", "fem3d", "dem", "dem3d", "fdem", "fdem3d")


SOLVER_KEYS = ("mode", "mesh", "meshFile", "scenario", "W", "H", "T", "nx", "ny", "frames")


def is_matpoint(kv):
    """deck qui n'est PAS un deck de solveur : point materiel (rockim matpoint, cles mp*), carte materiau
    (law + constantes seulement), deck temporaire de selftest — aucune cle de solveur (mode, mesh, scenario,
    W, H, T, nx, ny, frames). Ces fichiers ne passent jamais par l'audit de main.cpp."""
    return not any(k in kv for k in SOLVER_KEYS)


def preinit_guard(k, v, mode):
    """garde de main.cpp qui refuse AVANT l'audit (message court) ou None"""
    if k.startswith("hydro") and mode != "fdem":
        return "hydro* n'est implemente que pour mode = fdem (garde main.cpp, C6)"
    if k == "thermal" and mode != "fdem":
        return "'thermal' n'est implemente que pour mode = fdem (garde main.cpp)"
    if k == "beddingDip" and mode != "fdem":
        return "'bedding*' n'est implemente que pour mode = fdem (garde main.cpp)"
    if k == "law" and mode not in ("fem3d", "fdem", "fdem3d"):
        return "'law' n'est implemente que pour mode = fem3d | fdem | fdem3d (garde main.cpp)"
    if k == "phases" and mode not in ("fdem", "fdem3d"):
        return "'phases' n'est implemente que pour mode = fdem | fdem3d (garde main.cpp)"
    if k == "mesh" and v == "voronoi" and mode not in ("fdem", "fdem3d"):
        return "mesh = voronoi n'est implemente que pour mode = fdem | fdem3d (garde main.cpp)"
    if k == "mesh" and v == "file" and mode not in ("fdem", "fdem3d", "fem3d"):
        return "mesh = file n'est implemente que pour mode = fdem | fdem3d | fem3d (garde main.cpp)"
    if k == "mesh" and v not in ("grid", "voronoi", "file"):
        return "mesh = '%s' inconnu (grid | voronoi | file) (garde main.cpp)" % v
    if k == "mode" and v not in MODES:
        return "mode = '%s' inconnu (fem | fem3d | dem | dem3d | fdem | fdem3d)" % v
    return None


def modes_text(rd):
    rd = sorted(rd)
    if len(rd) == 1:
        return "cle du mode %s seulement" % rd[0]
    return "cles des modes %s" % ", ".join(rd)


def judge(kv, mode_of, known, dyn, obs, readers):
    """-> (fautes [(cle, ligne, kind, message)], n_dynamiques, mode)"""
    mode = kv.get("mode", ("fem", 0))[0]
    faults, ndyn = [], 0
    for k, (v, n) in kv.items():
        g = preinit_guard(k, v, mode)
        if k in obs:
            faults.append((k, n, "obsolete", "remplacee par '%s'" % obs[k]))
        elif g:
            faults.append((k, n, "garde", g))
        elif k in readers:
            rd = readers[k]
            if mode in rd or "shared" in rd:
                continue                                    # lecteur du mode / code partage : legitime
            faults.append((k, n, "mode", "sans effet en mode %s : %s" % (mode, modes_text(rd))))
        elif k in known:
            continue
        elif k.startswith(dyn):
            ndyn += 1
        else:
            cands = sorted((levenshtein(k, c), c) for c in known if c != k)
            cands = [c for d, c in cands if d <= 2] or \
                    [c for c in known if c.lower() == k.lower()]
            if cands:
                sug = " ou ".join("'%s'%s" % (c, " (%s)" % modes_text(readers[c])
                                              if c in readers and mode not in readers[c] and "shared" not in readers[c]
                                              else "")
                                  for c in cands[:3])
                faults.append((k, n, "typo", "vouliez-vous dire %s ?" % sug))
            else:
                faults.append((k, n, "unknown", "inconnue de rockim"))
    return faults, ndyn, mode


def fix_obsolete(path, kv, obs):
    """Remplace en place `old =` par `new =` (debut de ligne, espaces preserves) ; rend les remplacements."""
    done = []
    with io.open(path, encoding="utf-8", errors="ignore", newline="") as f:
        text = f.read()
    for old, new in obs.items():
        if old not in kv:
            continue
        rx = re.compile(r"(?m)^(\s*)%s(\s*=)" % re.escape(old))
        text, n = rx.subn(r"\g<1>%s\g<2>" % new, text)
        if n:
            done.append((old, new, n))
    if done:
        with io.open(path, "w", encoding="utf-8", newline="") as f:
            f.write(text)
    return done


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--roots", nargs="*", default=None)
    ap.add_argument("--fix-obsolete", action="store_true", dest="fix")
    ap.add_argument("--out", default=os.path.join(HERE, "scan_decks_%s.md" % time.strftime("%Y-%m-%d")))
    args = ap.parse_args()
    roots = args.roots or [ROOT, CDP]
    reg, mode_of, known, dyn, obs, readers = load_registry()

    decks = []
    for root in roots:
        for d, dirs, fns in os.walk(root):
            dirs[:] = [x for x in dirs if not SKIP_DIRS.match(x)]
            for fn in sorted(fns):
                if fn.endswith(".cfg"):
                    decks.append(os.path.join(d, fn))
    decks.sort()

    def rel(p):
        for root, tag in ((ROOT, ""), (CDP, "CDP/")):
            try:
                r = os.path.relpath(p, root)
                if not r.startswith(".."):
                    return tag + r.replace(os.sep, "/")
            except ValueError:
                pass
        return p.replace(os.sep, "/")

    rows, fixed, by_key, kinds, ndyn_tot, nbad, matpoint = [], [], {}, {}, 0, 0, []
    for p in decks:
        kv = read_deck(p)
        if is_matpoint(kv):
            matpoint.append(rel(p))
            continue
        faults, ndyn, mode = judge(kv, mode_of, known, dyn, obs, readers)
        ndyn_tot += ndyn
        if args.fix and any(k in obs for k in kv):
            for old, new, n in fix_obsolete(p, kv, obs):
                fixed.append((rel(p), old, new, kv[old][1]))
        if faults:
            nbad += 1
        rows.append((rel(p), mode, faults, ndyn))
        for k, n, kind, msg in faults:
            kinds[kind] = kinds.get(kind, 0) + 1
            by_key.setdefault((kind, k, msg), []).append("%s:%d" % (rel(p), n))

    now = time.strftime("%Y-%m-%d %H:%M:%S")
    L = []
    L.append("# Balayage statique des decks — regle C1 des cles (w21, lecteurs w22) — %s" % now)
    L.append("")
    L.append("Script : `tools/scan_decks.py` (analyse STATIQUE : registre `tools/keys_by_mode.json` du %s, %d cles connues avec "
             "leurs lecteurs, %d propres a un seul mode ; obsoletes `tools/obsolete_keys.json` : %s). Racines : %s. Dossiers `out_*`, "
             "`orig/` (copies avant edition) et `selftest_*` (bancs falsifiants, fautes voulues) ignores. "
             "Aucun run de rockim : les lectures conditionnelles ne sont pas jugeables ici (regle w22 : cle du registre dont le mode "
             "courant ou le code partage est lecteur = legitime ; cle qu'aucun lecteur du mode ne lit = « sans effet en mode »). "
             "Les gardes PRE-INIT de `main.cpp` (thermal/bedding/law/phases/mesh/hydro hors mode) sont modelisees ; les gardes de "
             "MAILLAGE (C3 : noeud orphelin, tet plat) ne le sont pas (il faudrait lire le .msh). %d fichier(s) sans cle de solveur "
             "(point materiel `rockim matpoint`, cartes materiau, decks temporaires de selftest) comptes a part, non juges."
             % (reg["_meta"]["generated"], len(known), len(mode_of),
                ", ".join("`%s` -> `%s`" % kv for kv in sorted(obs.items())) or "aucune",
                ", ".join("`%s`" % r.replace(os.sep, "/") for r in roots), len(matpoint)))
    L.append("")
    L.append("## Bilan")
    L.append("")
    L.append("| decks de solveur lus | sans faute | avec faute(s) | cles fautives | obsoletes | autre mode (lecteurs) | gardes pre-init main.cpp | fautes de frappe | inconnues | cles dynamiques (non jugees) | fichiers sans cle de solveur (non juges) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    L.append("| %d | %d | %d | %d | %d | %d | %d | %d | %d | %d | %d |" % (
        len(rows), len(rows) - nbad, nbad, sum(kinds.values()), kinds.get("obsolete", 0), kinds.get("mode", 0),
        kinds.get("garde", 0), kinds.get("typo", 0), kinds.get("unknown", 0), ndyn_tot, len(matpoint)))
    L.append("")
    L.append("Lecture : un deck « avec faute » serait REFUSE par `rockim_f2w22.exe` (`unknownKeys = error`, defaut) — sauf si la cle "
             "fautive est en realite lue par un getter que le registre ne voit pas (cle construite) ; le verdict definitif est "
             "celui du binaire. Une cle dynamique (`phase.<nom>.E`, `groupBond.<A>.<B>`, `gauge.<groupe>`...) n'est jugee qu'au run.")
    L.append("")
    if matpoint:
        L.append("Fichiers sans cle de solveur (matpoint / cartes / temporaires, non juges) : " + ", ".join("`%s`" % m for m in matpoint))
        L.append("")
    if args.fix:
        L.append("## Cles obsoletes CORRIGEES en place (%d remplacements)" % len(fixed))
        L.append("")
        if fixed:
            L.append("| deck | ligne | avant | apres |")
            L.append("|---|---|---|---|")
            for p, old, new, n in fixed:
                L.append("| `%s` | %d | `%s` | `%s` |" % (p, n, old, new))
        else:
            L.append("aucune")
        L.append("")
    L.append("## Fautes par cle (toutes ; a decider par Fernando, sauf les obsoletes deja corrigees)")
    L.append("")
    L.append("| type | cle | diagnostic | decks | ou (deck:ligne) |")
    L.append("|---|---|---|---|---|")
    order = {"obsolete": 0, "mode": 1, "garde": 2, "typo": 3, "unknown": 4}
    for (kind, k, msg), where in sorted(by_key.items(), key=lambda x: (order[x[0][0]], -len(x[1]), x[0][1])):
        w = ", ".join("`%s`" % x for x in where[:12]) + (" … (+%d)" % (len(where) - 12) if len(where) > 12 else "")
        L.append("| %s | `%s` | %s | %d | %s |" % (kind, k, msg.replace("|", "\\|"), len(where), w))
    L.append("")
    L.append("## Par deck (decks avec au moins une faute : %d)" % nbad)
    L.append("")
    L.append("| deck | mode | fautes |")
    L.append("|---|---|---|")
    for p, mode, faults, ndyn in rows:
        if not faults:
            continue
        L.append("| `%s` | %s | %s |" % (p, mode, " ; ".join(
            "`%s` (l.%d, %s : %s)" % (k, n, kind, msg.replace("|", "\\|")) for k, n, kind, msg in faults)))
    L.append("")
    with io.open(args.out, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(L) + "\n")
    print("%d decks, %d avec faute(s), %d cles fautives (%s), %d dynamiques non jugees -> %s" % (
        len(decks), nbad, sum(kinds.values()),
        ", ".join("%s %d" % kv for kv in sorted(kinds.items())), ndyn_tot, args.out))
    if args.fix:
        print("%d cles obsoletes corrigees en place" % len(fixed))
    return 0


if __name__ == "__main__":
    sys.exit(main())
