# -*- coding: utf-8 -*-
"""gen_keys_by_mode.py — registre des cles de deck PAR MODE (C6, plan de robustesse 2026-09-05, w20).

METHODE (documentee ici et dans DOCUMENTATION_rockim.md §8 « Gardes des entrees ») :
  1. On lit les getters de configuration dans src/*.cpp et include/rockim/*.hpp (sauf *.w0, Guards.hpp,
     KeysByMode.hpp) : getd / getb / gets / geti / reqs / reqd / has avec une cle LITTERALE complete,
     c'est-a-dire `getX("cle"` suivi de `,` ou `)`. Une cle construite (`"groupBond." + nom`) n'est pas
     un litteral complet : elle n'entre PAS au registre (jamais refusee). keysWithPrefix("p") donne un
     prefixe commun (jamais refuse).
  2. Chaque fichier a un proprietaire : FemSolver -> fem, Fem3dSolver -> fem3d, DemSolver -> dem,
     Dem3dSolver -> dem3d, FdemSolver -> fdem, Fdem3dSolver -> fdem3d ; tout autre fichier (main,
     Material, MatLaw, Tool*, PotentialContact, RandomField, Tessellation*, VtkWriter, Yan*, Yang*...)
     est du code PARTAGE.
  3. Une cle lue par le code partage, ou par PLUSIEURS solveurs, est « commune » : jamais refusee.
     Une cle lue par UN SEUL solveur est « du mode X » : un deck d'un autre mode qui la porte est refuse
     (« cle X sans effet en mode Y ») par keysbymode::check() dans main.cpp.
  4. En cas de doute (cle listee dans ALWAYS_COMMON ci-dessous), la cle est commune.
  5. C1 (w21, 2026-09-05) : le meme registre sert de liste des cles CONNUES du code (kKnown = communes +
     propres a un mode) pour la suggestion de Levenshtein ; les familles DYNAMIQUES (`"phase." + nom`, `gb.`,
     `contactMu.`, `groupBond.`, `groupPhase.`, `groupVel.`, `gauge.` : litteral `"xxx." +` dans le code)
     sont relevees comme prefixes (kDynamicPrefix) ; la table des cles OBSOLETES vient de
     tools/obsolete_keys.json (kObsolete : « cle obsolete X : remplacee par Z »).
  6. w22 (2026-09-05, relecture D1) : les LECTEURS de chaque cle sont exportes (kReaders : cle -> "fdem fdem3d",
     "shared"...). La regle de KeyGuard.hpp pour une cle du deck NON consommee devient : legitime si et seulement si
     le mode courant est un de ses lecteurs (lecture conditionnelle : capP0 si dprCap = true) ou si le code partage
     la lit ("shared", ou ALWAYS_COMMON) ; sinon « cle X sans effet en mode Y : cles des modes Z1, Z2 ». Une cle
     lue par >= 2 solveurs mais par AUCUN du mode courant (jointSoftening, insertion, bulkDamage en fem3d ;
     bond*, packing hors dem/dem3d ; kpFactor hors fem/fem3d) n'est donc plus « commune = legitime » (trou w21).
Sorties : tools/keys_by_mode.json (registre lisible) et include/rockim/KeysByMode.hpp (table compilee,
NE PAS EDITER A LA MAIN). Le registre doit etre regenere a chaque nouvelle cle : le build ne le fait pas.

usage : python tools/gen_keys_by_mode.py                  # regenere json + hpp
        python tools/gen_keys_by_mode.py --scan-decks D1 D2 ...   # essai a blanc : decks (*.cfg) qui
                                                                  # seraient refuses, sans rien ecrire
"""
import io, json, os, re, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
OWNER = {"FemSolver": "fem", "Fem3dSolver": "fem3d", "DemSolver": "dem", "Dem3dSolver": "dem3d",
         "FdemSolver": "fdem", "Fdem3dSolver": "fdem3d"}
# Cles a garder communes quoi qu'en dise le comptage (doute = commune, jamais refusee a tort).
ALWAYS_COMMON = {"mode", "mesh", "meshFile", "outputDir", "frames", "historyFlush", "nanCheckEvery",
                 "seed", "T", "W", "H", "D", "nx", "ny", "nz", "scenario", "dtFactor", "cfl"}
RX_GET = re.compile(r'\b(?:getd|getb|gets|geti|reqs|reqd|has)\s*\(\s*"([A-Za-z0-9_.]+)"\s*[,)]')
RX_PFX = re.compile(r'\bkeysWithPrefix\s*\(\s*"([A-Za-z0-9_.]+)"\s*\)')
# famille dynamique : litteral "xxx." concatene a un nom (phase.<nom>.E, groupBond.<A>.<B>...)
RX_DYN = re.compile(r'"([A-Za-z0-9_]+\.)"\s*\+')
OBSOLETE_JSON = os.path.join(HERE, "obsolete_keys.json")


def load_obsolete():
    if not os.path.isfile(OBSOLETE_JSON):
        return {}
    with io.open(OBSOLETE_JSON, encoding="utf-8") as f:
        return json.load(f).get("obsolete", {})


def scan():
    files = []
    for d in ("src", os.path.join("include", "rockim")):
        for fn in sorted(os.listdir(os.path.join(ROOT, d))):
            if not (fn.endswith(".cpp") or fn.endswith(".hpp")):
                continue
            if fn in ("Guards.hpp", "KeysByMode.hpp"):
                continue
            files.append(os.path.join(d, fn))
    readers = {}                                   # cle -> set(proprietaires)
    prefixes = set()
    dynamic = set()
    for rel in files:
        base = os.path.basename(rel).split(".")[0]
        owner = OWNER.get(base, "shared")
        txt = io.open(os.path.join(ROOT, rel), encoding="utf-8", errors="ignore").read()
        for k in RX_GET.findall(txt):
            readers.setdefault(k, set()).add(owner)
        for p in RX_PFX.findall(txt):
            prefixes.add(p)
        for p in RX_DYN.findall(txt):
            dynamic.add(p)
    modes = {m: [] for m in sorted(set(OWNER.values()))}
    common = []
    for k in sorted(readers):
        own = readers[k]
        if k in ALWAYS_COMMON or "shared" in own or len(own) > 1:
            common.append(k)
        else:
            modes[next(iter(own))].append(k)
    return files, readers, modes, common, sorted(prefixes), sorted(dynamic)


def write(files, readers, modes, common, prefixes, dynamic):
    obsolete = load_obsolete()
    reg = {
        "_meta": {
            "generated": time.strftime("%Y-%m-%d %H:%M:%S"),
            "generator": "tools/gen_keys_by_mode.py (voir la docstring : methode)",
            "rule": "cle lue par UN SEUL solveur = cle de ce mode, refusee dans les autres modes ; "
                    "lue par le code partage ou par >= 2 solveurs = commune, jamais refusee ; "
                    "cle absente du registre (construite dynamiquement) = jamais refusee",
            "files_scanned": files,
            "always_common": sorted(ALWAYS_COMMON),
            "n_keys": len(readers),
            "n_common": len(common),
            "n_mode_specific": sum(len(v) for v in modes.values()),
            "n_obsolete": len(obsolete),
            "rule_c1": "C1 (w21) : cle du deck NON consommee ET hors registre (ni commune, ni du mode, ni "
                       "famille dynamique lue) = faute ; cle du registre non consommee (lecture "
                       "conditionnelle) = legitime ; cle obsolete = « remplacee par » ; cle d'un autre "
                       "mode = « sans effet en mode » ; source des obsoletes : tools/obsolete_keys.json",
            "rule_w22": "w22 (relecture D1) : une cle du deck NON consommee est legitime ssi le mode courant "
                        "figure dans ses lecteurs (readers) ou si 'shared' y figure (code partage, ALWAYS_COMMON) ; "
                        "sinon « cle X sans effet en mode Y : cles des modes Z1, Z2 » (kReaders dans KeysByMode.hpp)",
        },
        "modes": modes,
        "common": common,
        "prefixes_common": prefixes,
        "prefixes_dynamic": dynamic,
        "obsolete": {k: v.get("new") for k, v in sorted(obsolete.items())},
        # lecteurs ; une cle ALWAYS_COMMON porte aussi 'shared' (doute = commune, jamais refusee)
        "readers": {k: sorted(v | ({"shared"} if k in ALWAYS_COMMON else set()))
                    for k, v in sorted(readers.items())},
    }
    with io.open(os.path.join(HERE, "keys_by_mode.json"), "w", encoding="utf-8", newline="\n") as f:
        json.dump(reg, f, indent=1, ensure_ascii=False, sort_keys=False)
        f.write("\n")
    lines = [
        "#pragma once",
        "// KeysByMode.hpp — GENERE par tools/gen_keys_by_mode.py le %s : NE PAS EDITER." % reg["_meta"]["generated"],
        "// Registre lisible : tools/keys_by_mode.json (methode dans la docstring du generateur).",
        "// Table des cles PROPRES A UN MODE (lues par un seul solveur) ; toute autre cle est commune.",
        "#include <stdexcept>",
        "#include <string>",
        '#include "rockim/Config.hpp"',
        "",
        "namespace rockim {",
        "namespace keysbymode {",
        "",
        "struct Entry { const char* key; const char* mode; };",
        "",
        "inline const Entry kTable[] = {",
    ]
    for m in sorted(modes):
        for k in modes[m]:
            lines.append('    {"%s", "%s"},' % (k, m))
    lines += [
        '    {nullptr, nullptr}',
        "};",
        "",
        "// C1 (w21) : TOUTES les cles connues du code (communes + propres a un mode),",
        "// ordre lexicographique — « presente dans le registre = legitime meme si non",
        "// consommee » (lecture conditionnelle) et candidates de la suggestion.",
        "inline const char* const kKnown[] = {",
    ]
    for k in sorted(readers):
        lines.append('    "%s",' % k)
    lines += [
        "    nullptr",
        "};",
        "",
        "// C1 (w21) : familles de cles construites dynamiquement (\"phase.\" + nom...) :",
        "// une cle de ces familles non lue = nom de phase/groupe inconnu.",
        "inline const char* const kDynamicPrefix[] = {",
    ]
    for p in dynamic:
        lines.append('    "%s",' % p)
    lines += [
        "    nullptr",
        "};",
        "",
        "// C1 (w21) : prefixes valides par le code lui-meme (keysWithPrefix).",
        "inline const char* const kCommonPrefix[] = {",
    ]
    for p in prefixes:
        lines.append('    "%s",' % p)
    lines += [
        "    nullptr",
        "};",
        "",
        "// C1 (w21) : cles OBSOLETES -> nom courant (tools/obsolete_keys.json).",
        "struct Obsolete { const char* old; const char* now; };",
        "inline const Obsolete kObsolete[] = {",
    ]
    for k in sorted(obsolete):
        lines.append('    {"%s", "%s"},' % (k, obsolete[k].get("new", "")))
    lines += [
        "    {nullptr, nullptr}",
        "};",
        "",
        "// w22 (relecture D1) : LECTEURS de chaque cle connue (proprietaires des fichiers",
        "// qui la lisent : fem fem3d dem dem3d fdem fdem3d, ou shared = code partage /",
        "// ALWAYS_COMMON), separes par un espace. KeyGuard : une cle du deck NON consommee",
        "// est legitime ssi le mode courant ou shared figure ici ; sinon « sans effet en",
        "// mode Y : cles des modes Z1, Z2 ».",
        "struct Readers { const char* key; const char* readers; };",
        "inline const Readers kReaders[] = {",
    ]
    for k in sorted(readers):
        rd = sorted(readers[k] | ({"shared"} if k in ALWAYS_COMMON else set()))
        lines.append('    {"%s", "%s"},' % (k, " ".join(rd)))
    lines += [
        "    {nullptr, nullptr}",
        "};",
        "",
        "// C6 (w20) : refuse toute cle du deck qui n'est lue que par le solveur d'un",
        "// AUTRE mode (« cle X sans effet en mode Y »). Les cles communes et les cles",
        "// construites dynamiquement (absentes de la table) passent toujours.",
        "inline void check(const Config& cfg, const std::string& mode) {",
        "    for (const auto& k : cfg.keys())",
        "        for (const Entry* e = kTable; e->key; ++e)",
        "            if (k == e->key && mode != e->mode)",
        "                throw std::runtime_error(\"cle '\" + k + \"' sans effet en mode \" + mode",
        "                    + \" : cle du mode \" + e->mode + \" seulement (registre \"",
        "                    \"tools/keys_by_mode.json) ; retirez-la du deck ou changez de mode\");",
        "}",
        "",
        "} // namespace keysbymode",
        "} // namespace rockim",
        "",
    ]
    with io.open(os.path.join(ROOT, "include", "rockim", "KeysByMode.hpp"), "w", encoding="utf-8",
                 newline="\n") as f:
        f.write("\n".join(lines))
    print("registre : %d cles lues, %d communes, %d propres a un mode (%s), %d prefixes communs, "
          "%d familles dynamiques, %d obsoletes" % (
        len(readers), len(common), reg["_meta"]["n_mode_specific"],
        ", ".join("%s %d" % (m, len(v)) for m, v in sorted(modes.items())), len(prefixes),
        len(dynamic), len(obsolete)))


def read_cfg(p):
    kv = {}
    for line in io.open(p, encoding="utf-8", errors="ignore"):
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            kv[k.strip()] = v.strip()
    return kv


def scan_decks(dirs, modes):
    table = {}
    for m, ks in modes.items():
        for k in ks:
            table[k] = m
    n = nbad = 0
    for d in dirs:
        for root, _, fns in os.walk(d):
            for fn in sorted(fns):
                if not fn.endswith(".cfg"):
                    continue
                p = os.path.join(root, fn)
                kv = read_cfg(p)
                mode = kv.get("mode", "fem")
                bad = [(k, table[k]) for k in kv if k in table and table[k] != mode]
                hyd = [k for k in kv if k.startswith("hydro")] if mode != "fdem" else []
                n += 1
                if bad or hyd:
                    nbad += 1
                    print("REFUSE  %-70s mode=%-6s %s%s" % (
                        os.path.relpath(p, ROOT), mode,
                        " ".join("%s(%s)" % b for b in bad),
                        (" hydro:" + ",".join(hyd)) if hyd else ""))
    print("%d decks lus, %d seraient refuses" % (n, nbad))
    return nbad


if __name__ == "__main__":
    files, readers, modes, common, prefixes, dynamic = scan()
    if len(sys.argv) > 2 and sys.argv[1] == "--scan-decks":
        sys.exit(1 if scan_decks(sys.argv[2:], modes) else 0)
    write(files, readers, modes, common, prefixes, dynamic)
