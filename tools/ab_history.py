#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# ab_history.py — deux history.csv du MEME deck, deux binaires (ou deux
# nombres de fils) : sont-ils bit-identiques sur leurs instants communs ?
#
#   python tools/ab_history.py out_A out_B
#
# Les deux runs peuvent avoir une cadence d ecriture differente (frames, T) :
# les lignes sont appariees par leur t ECRIT (meme dt, meme pas -> meme
# chaine), pas par leur rang. Verdict : IDENTIQUE (aucune cellule ne differe
# sur les instants communs) ou DIFFERENT (premier instant, pire colonne).
# Sert au controle d une parallelisation : le resultat ne doit dependre ni du
# nombre de fils ni du decoupage.
# ---------------------------------------------------------------------------
import csv
import sys


def load(path):
    with open(path + "/history.csv", newline="") as f:
        rd = csv.reader(f)
        head = next(rd)
        rows = [r for r in rd if r]
    return head, rows


def main():
    if len(sys.argv) < 3:
        print("usage: ab_history.py out_A out_B")
        return 2
    ha, ra = load(sys.argv[1])
    hb, rb = load(sys.argv[2])
    if ha != hb:
        print("en-tetes differents :\n  A %s\n  B %s" % (ha, hb))
        return 1
    tb = {r[0]: r for r in rb}
    pairs = [(r, tb[r[0]]) for r in ra if r[0] in tb]
    print("A %d lignes, B %d lignes, %d instants communs (dernier t commun %s s)"
          % (len(ra), len(rb), len(pairs), pairs[-1][0][0] if pairs else "-"))
    if not pairs:
        print("aucun instant commun : cadences incompatibles")
        return 1
    first = None
    worst = {}
    for a, b in pairs:
        for j, k in enumerate(ha):
            if a[j] == b[j]:
                continue
            try:
                d = abs(float(a[j]) - float(b[j]))
            except ValueError:
                d = float("inf")
            if k not in worst or d > worst[k][0]:
                worst[k] = (d, a[0], a[j], b[j])
            if first is None:
                first = (a[0], k)
    if first is None:
        print("VERDICT : IDENTIQUE au caractere pres sur %d instants communs"
              % len(pairs))
        return 0
    print("VERDICT : DIFFERENT — premiere divergence a t = %s s, colonne %s"
          % first)
    for k, (d, t, a, b) in sorted(worst.items(), key=lambda kv: -kv[1][0])[:12]:
        print("  %-14s ecart max %.3e (t = %s : A %s | B %s)" % (k, d, t, a, b))
    return 1


if __name__ == "__main__":
    sys.exit(main())
