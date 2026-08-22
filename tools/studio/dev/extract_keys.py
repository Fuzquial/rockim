#!/usr/bin/env python3
"""extract_keys.py — extrait le registre des clés de config depuis src/*.cpp.

La vérité du parseur est l'ensemble des sites d'appel cfg.getd/geti/getb/gets
(clé + défaut) et cfg.reqd/reqs (clé requise). Ce script les balaie, les
regroupe par clé et écrit tools/studio/rockim_studio/model/keys_extracted.json :

    {"jointXi": {"type": "float", "required": false,
                 "sites": [{"file": "FdemSolver.cpp", "default": "0.05"},
                           {"file": "Fdem3dSolver.cpp", "default": "0.05"}]},
     ...}

Le type est déduit de la méthode (getd/reqd -> float, geti -> int,
getb -> bool, gets/reqs -> str) ; un conflit de type entre sites est signalé.
Les clés DYNAMIQUES (construites à l'exécution : groupBond.<A>.<B>,
groupVel.<nom>, gauge.<nom>, phase.<nom>.<prop>, groupPhase.<nom>) ne sont
pas des littéraux dans le C++ ; elles sont déclarées ici comme préfixes.

Usage :  python3 tools/studio/dev/extract_keys.py [--check]
  --check : ne réécrit rien ; échoue (code 1) si le JSON commité diffère de
            l'extraction courante — c'est le garde anti-dérive (risque R3 de
            la spec 006) à brancher dans la suite de vérification.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SOURCES = (ROOT / "src", ROOT / "include/rockim")
OUT = ROOT / "tools/studio/rockim_studio/model/keys_extracted.json"

METHOD_TYPE = {"getd": "float", "reqd": "float", "geti": "int",
               "getb": "bool", "gets": "str", "reqs": "str"}

# Préfixes de familles de clés construites dynamiquement (DOCUMENTATION §5.1,
# §5.2, §5.4) — le balayage des littéraux ne peut pas les voir.
DYNAMIC_PREFIXES = {
    "phase.": "propriété de phase minérale : phase.<nom>.<prop>",
    "groupPhase.": "matériau d'un corps : groupPhase.<nom> = <phase>",
    "groupVel.": "vitesse initiale d'un corps : groupVel.<nom> = vx vy vz",
    "groupBond.": "liaison entre corps : groupBond.<A>.<B> = joints",
    "gauge.": "jauge en tranche : gauge.<nom> = \"z0 z1\"",
}

CALL_RE = re.compile(
    r'\.(getd|geti|getb|gets|reqd|reqs)\(\s*"([A-Za-z0-9_.]+)"'
    r'(?:\s*,\s*([^)]*?))?\)')


def extract():
    keys = {}
    files = [p for d in SOURCES for pat in ("*.cpp", "*.hpp")
             for p in sorted(d.glob(pat))]
    for cpp in files:
        text = cpp.read_text(encoding="utf-8", errors="replace")
        for m in CALL_RE.finditer(text):
            method, key, default = m.group(1), m.group(2), m.group(3)
            entry = keys.setdefault(key, {"type": None, "required": False,
                                          "sites": []})
            ktype = METHOD_TYPE[method]
            if entry["type"] is None:
                entry["type"] = ktype
            elif entry["type"] != ktype:
                # gets() relisant une clé numérique (dispatch) : str l'emporte
                # seulement si aucun site numérique n'existe ; sinon on garde
                # le numérique et on note le conflit.
                if ktype != "str" and entry["type"] == "str":
                    entry["type"] = ktype
                entry.setdefault("type_conflicts", []).append(
                    {"file": cpp.name, "method": method})
            if method in ("reqd", "reqs"):
                entry["required"] = True
            site = {"file": cpp.name}
            if default is not None:
                site["default"] = default.strip()
            entry["sites"].append(site)
    return {"keys": dict(sorted(keys.items())),
            "dynamic_prefixes": DYNAMIC_PREFIXES}


def main():
    data = extract()
    payload = json.dumps(data, indent=1, ensure_ascii=False,
                         sort_keys=True) + "\n"
    if "--check" in sys.argv:
        if not OUT.exists():
            print(f"ECHEC : {OUT} absent — lancer extract_keys.py sans --check")
            return 1
        if OUT.read_text(encoding="utf-8") != payload:
            print("ECHEC : keys_extracted.json ne correspond plus a src/*.cpp "
                  "— relancer tools/studio/dev/extract_keys.py et committer")
            return 1
        print(f"OK : registre a jour ({len(data['keys'])} cles)")
        return 0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(payload, encoding="utf-8")
    n_conf = sum(1 for v in data["keys"].values() if "type_conflicts" in v)
    print(f"{len(data['keys'])} cles extraites -> {OUT.relative_to(ROOT)}"
          f" ({n_conf} conflits de type notes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
