#!/usr/bin/env python3
"""Depouillement du rejeu tunnel avec le contact corrige (b869dae).

    python3 extract.py RUNDIR:LOG:TBLOC:etiquette ... > metrics_fix.json

Reutilise les outils du depot, sans les modifier :
  - analyse() de docs/rapport_guide/figures/tunnel/fig_tunnel_rejeu.py
    (rompus, modes, R_EDZ, demi-axes, longueur, paroi : memes definitions
    que rejeu_metrics.json du 6 octobre) ;
  - tunnel_edz/tools/block_sizes.py --t TBLOC --r 25 ;
  - tunnel_edz/tools/nucleation_vs_propagation.py ;
  - le resume du journal (bilan V2/B4, naissances, duree).
"""
import json
import os
import re
import subprocess
import sys

ROOT = "/home/user/rockim"
sys.path.insert(0, os.path.join(ROOT, "docs/rapport_guide/figures/tunnel"))
import fig_tunnel_rejeu as F  # noqa: E402

TOOLS = os.path.join(ROOT, "tunnel_edz/tools")


def num(pat, s, cast=float):
    m = re.search(pat, s)
    return cast(m.group(1)) if m else None


def blocks(run, t):
    out = subprocess.run([sys.executable, os.path.join(TOOLS, "block_sizes.py"),
                          run, "--t", str(t), "--r", "25"],
                         capture_output=True, text=True).stdout
    big = re.search(r"5 plus gros blocs \[m2\] : (.*)", out)
    return dict(
        bloc_trame_t=num(r"\(t = ([0-9.]+) s\)", out),
        blocs=num(r"blocs : +(\d+)", out, int),
        blocs_mono=num(r"blocs mono-element : +(\d+)", out, int),
        blocs_mono_pct=num(r"blocs mono-element : +\d+ +\( *([0-9.]+) %\)", out),
        blocs_5_plus_gros_m2=[float(x) for x in big.group(1).split()] if big else None,
        blocs_hors_massif_moy_m2=num(r"hors massif : aire moyenne +([0-9.]+)", out),
        blocs_txt=out)


def propagation(run):
    out = subprocess.run([sys.executable, os.path.join(TOOLS, "nucleation_vs_propagation.py"),
                          run], capture_output=True, text=True).stdout
    m = re.search(r"BILAN : (\d+) propagations / (\d+) nucleations -> ([0-9.]+) %", out)
    return dict(n_prop=int(m.group(1)), n_nuc=int(m.group(2)),
                propagation_pct=float(m.group(3))) if m else {}


def journal(log):
    s = open(log, errors="ignore").read()
    g = r"(-?[0-9.]+(?:e[+-]?\d+)?)"
    return dict(
        contact_travail_net_Jm=num(r"net work injected by general contact: " + g, s),
        energie_cohesive_Jm=num(r"joints +: " + g + " J/m cohesif", s),
        frottement_Jm=num(r"dont frottement " + g, s),
        cundall_Jm=num(r"Cundall +: " + g, s),
        residu_B4_pct=num(r"residu +: .*?\(" + g + " % de l'echelle\)", s),
        naissances_paires=num(r"contact, naissances : (\d+) paires", s, int),
        energie_naissance_Jm=num(r"energie de recouvrement a la naissance " + g, s),
        naissance_statut=num(r"contact, naissances : .*? J/m \((\w+)", s, str),
        renaissances_paires=num(r"renaissances \(gcBirth = offset\) : (\d+) paires", s, int),
        energie_renaissance_Jm=num(r"energie neutralisee " + g, s),
        candidats_ajoutes=num(r"contactCandidates = vertex, (\d+) elements-pas", s, int),
        ke_fin_Jm=num(r"block kinetic energy at end: " + g, s),
        dt=num(r"dt = " + g + " s", s),
        duree_s=num(r"wall time: " + g + " s", s),
        avertissements=[l for l in s.splitlines() if "AVERTISSEMENT" in l or "WARNING" in l])


def main():
    out = []
    for a in sys.argv[1:]:
        run, log, tb, lab = a.split(":", 3)
        res, _ = F.analyse(run)
        res["label"] = lab
        res.update(journal(log))
        b = blocks(run, float(tb))
        b.pop("blocs_txt")
        res.update(b)
        res.update(propagation(run))
        out.append(res)
        print(json.dumps(res, ensure_ascii=False), file=sys.stderr)
    json.dump(out, sys.stdout, indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
