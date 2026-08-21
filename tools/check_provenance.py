#!/usr/bin/env python3
"""Verifie que chaque extrait de code cite dans un .tex correspond au depot.

Pour chaque bloc

    \\provenance{src/Fichier.cpp:1234 --- libelle}
    \\begin{lstlisting}
    ...code...
    \\end{lstlisting}

le script relit le fichier a la reference git indiquee et controle, ligne a ligne :

  ABSENTE  la ligne citee n'existe nulle part dans le fichier
           -> texte ajoute, traduit ou reformule dans le .tex
  DECALEE  la ligne existe mais loin du numero annonce
           -> renvoi obsolete (mauvaise branche, code deplace)
  OK       la ligne existe a +/- tolerance du numero annonce

Une ligne reduite a "..." est un saut explicite et n'est pas controlee.

    python check_provenance.py solveur_rockim.tex --repo <chemin> --ref article-exact
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
import unicodedata

PROV = re.compile(r"\\provenance\{(.+?)\}", re.S)
LOC = re.compile(r"([\w./\-]+\.(?:cpp|hpp|h|py|f)):(\d+)")
LST = re.compile(r"\\begin\{lstlisting\}(?:\[.*?\])?\n(.*?)\\end\{lstlisting\}", re.S)
SRC = re.compile(r"\\src\{([\w./\-]+\.(?:cpp|hpp|h|py|f)):(\d+)(?:-(\d+))?\}")


def norm(s: str) -> str:
    """Compare a l'identique modulo espaces, tirets et accents."""
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = s.replace("\u2014", "-").replace("\u2013", "-").replace("\u2019", "'")
    return re.sub(r"\s+", " ", s).strip()


def git_file(repo: str, ref: str, path: str) -> list[str] | None:
    try:
        out = subprocess.run(
            ["git", "-C", repo, "show", f"{ref}:{path}"],
            capture_output=True, check=True, text=True, errors="replace",
        ).stdout
    except subprocess.CalledProcessError:
        return None
    return out.splitlines()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tex")
    ap.add_argument("--repo", required=True, help="racine du depot git")
    ap.add_argument("--ref", required=True, help="commit, tag ou branche faisant foi")
    ap.add_argument("--tol", type=int, default=3,
                    help="ecart tolere en lignes (defaut 3)")
    ap.add_argument("--quiet", action="store_true", help="ne lister que les defauts")
    args = ap.parse_args()

    tex = open(args.tex, encoding="utf-8").read()
    cache: dict[str, list[str] | None] = {}

    def load(path: str) -> list[str] | None:
        if path not in cache:
            got = git_file(args.repo, args.ref, path)
            if got is None and not path.startswith(("src/", "include/")):
                for pre in ("src/", "include/rockim/", "include/"):
                    got = git_file(args.repo, args.ref, pre + path)
                    if got is not None:
                        break
            cache[path] = got
        return cache[path]

    # --- 1. blocs provenance + listing -----------------------------------
    blocs = []
    for m in PROV.finditer(tex):
        lst = LST.search(tex, m.end())
        if lst:
            blocs.append((m.group(1).strip(), lst.group(1),
                          tex[:m.start()].count("\n") + 1))

    n_ok = n_abs = n_dec = n_skip = 0
    print(f"=== {args.tex}")
    print(f"=== reference : {args.ref}\n")

    for prov, code, texline in blocs:
        loc = LOC.search(prov)
        libelle = prov.split("---")[-1].strip() if "---" in prov else prov
        if not loc:
            n_skip += 1
            if not args.quiet:
                print(f"[.tex:{texline}] SANS RENVOI  {prov[:60]}")
            continue
        path, cited = loc.group(1), int(loc.group(2))
        src = load(path)
        if src is None:
            print(f"[.tex:{texline}] FICHIER INTROUVABLE  {path} @ {args.ref}")
            continue

        index: dict[str, list[int]] = {}
        for i, line in enumerate(src, 1):
            index.setdefault(norm(line), []).append(i)

        span = len(code.splitlines())
        defauts = []
        for raw in code.splitlines():
            if not raw.strip() or raw.strip() == "...":
                continue
            hits = index.get(norm(raw))
            if not hits:
                defauts.append(("ABSENTE", raw.strip(), None))
            else:
                best = min(hits, key=lambda h: abs(h - cited))
                if abs(best - cited) > args.tol + span:
                    defauts.append(("DECALEE", raw.strip(), best))

        if not defauts:
            n_ok += 1
            if not args.quiet:
                print(f"[.tex:{texline}] OK  {path}:{cited}  {libelle}")
        else:
            print(f"[.tex:{texline}] {path}:{cited}  {libelle}")
            for kind, txt, at in defauts:
                if kind == "ABSENTE":
                    n_abs += 1
                    print(f"      ABSENTE  {txt[:76]}")
                else:
                    n_dec += 1
                    print(f"      DECALEE  ligne reelle {at} : {txt[:60]}")

    # --- 2. renvois \src{} isoles ----------------------------------------
    bad_src = []
    for m in SRC.finditer(tex):
        path, a = m.group(1), int(m.group(2))
        src = load(path)
        if src is None:
            bad_src.append((path, a, "fichier introuvable"))
        elif a > len(src):
            bad_src.append((path, a, f"hors fichier ({len(src)} lignes)"))

    print(f"\n--- blocs : {n_ok} conformes, {n_abs} lignes absentes, "
          f"{n_dec} lignes decalees, {n_skip} sans renvoi")
    if bad_src:
        print(f"--- renvois src hors bornes : {len(bad_src)}")
        for p, a, why in bad_src:
            print(f"      {p}:{a} - {why}")
    return 1 if (n_abs or n_dec or bad_src) else 0


if __name__ == "__main__":
    sys.exit(main())
