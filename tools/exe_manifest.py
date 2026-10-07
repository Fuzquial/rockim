#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# exe_manifest.py — empreinte des binaires existants de rockim_f2 (mesure B6,
# plan de robustesse 2026-09-05) : pour chaque rockim_f2*.exe de la racine,
# SHA-256, taille, date, et la ligne du tableau « Binaires » de
# etude_lois_fem/JOURNAL.md qui le decrit (ou « non documente »).
#
#   python tools/exe_manifest.py                 # ecrit tools/exe_manifest.json
#   python tools/exe_manifest.py --out autre.json
#
# Le tableau du JOURNAL a une ligne a motif (`rockim_f2w1*.exe` = les liens
# partiels w1a..w1d) : `*` est lu comme « lettres seulement » (pas de chiffre),
# pour ne pas attribuer cette ligne a w10..w19. Le champ `journal_mentions`
# donne en plus les numeros de ligne ou le nom complet de l'exe apparait
# ailleurs dans le JOURNAL (aide a la lecture, pas une documentation).
# Les binaires sont LUS, jamais modifies.
# ---------------------------------------------------------------------------
import argparse, glob, hashlib, json, os, re, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
JOURNAL = os.path.join(ROOT, "etude_lois_fem", "JOURNAL.md")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def journal_rows(path):
    """Lignes du tableau des binaires : (numero, motif regex, texte)."""
    rows = []
    if not os.path.isfile(path):
        return rows
    rx = re.compile(r"^\|\s*\**`(rockim_f2[^`]*\.exe)`")
    with open(path, encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            m = rx.match(line)
            if m:
                pat = re.escape(m.group(1)).replace(r"\*", r"[a-z_]*")
                rows.append((i, re.compile("^" + pat + "$"), line.rstrip("\n")))
    return rows


def journal_mentions(path, name):
    out = []
    if not os.path.isfile(path):
        return out
    with open(path, encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            if name in line:
                out.append(i)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "exe_manifest.json"))
    ap.add_argument("--journal", default=JOURNAL)
    args = ap.parse_args()
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

    rows = journal_rows(args.journal)
    exes = sorted(glob.glob(os.path.join(ROOT, "rockim_f2*.exe")), key=os.path.getmtime)
    entries = []
    n_doc = 0
    for p in exes:
        name = os.path.basename(p)
        st = os.stat(p)
        row = next(((i, txt) for (i, rx, txt) in rows if rx.match(name)), None)
        mentions = [i for i in journal_mentions(args.journal, name) if not row or i != row[0]]
        e = dict(name=name, sha256=sha256_file(p), size=st.st_size,
                 mtime=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(st.st_mtime)),
                 journal_line=row[0] if row else None,
                 journal_row=row[1] if row else "non documente",
                 journal_mentions=mentions)
        n_doc += bool(row)
        entries.append(e)
        print(f"{name:<28s} {e['sha256'][:16]}  {st.st_size:>9d}  {e['mtime']}  "
              f"{'JOURNAL l.' + str(row[0]) if row else 'non documente'}"
              f"{'  (mentions ' + ','.join(map(str, mentions[:6])) + ('...' if len(mentions) > 6 else '') + ')' if mentions else ''}")
    out = dict(_meta=dict(date=time.strftime("%Y-%m-%d %H:%M:%S"), root=ROOT,
                          journal=os.path.relpath(args.journal, ROOT),
                          n_exe=len(entries), n_documentes=n_doc,
                          note="SHA-256 du fichier exe ; journal_row = ligne du tableau « Binaires » de "
                               "etude_lois_fem/JOURNAL.md dont le motif (`*` = lettres seulement) "
                               "correspond au nom, sinon « non documente » ; journal_mentions = "
                               "autres lignes du JOURNAL citant le nom complet"),
               exe=entries)
    with open(args.out, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(f"[manifest] {len(entries)} exe, {n_doc} documentes dans le tableau du JOURNAL -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
