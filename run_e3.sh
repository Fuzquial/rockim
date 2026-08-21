#!/bin/bash
# ESSAI 3 — protocole article : injection apres excavation (hydroStart).
# Les deux cas EN SERIE, isotrope puis anisotrope : un seul job a la fois.
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"

run () {           # $1 = config, $2 = dossier de sortie, $3 = journal, $4 = libelle
    t0=$(date +%s)
    echo "=== $4 : debut $(date +%H:%M:%S) ===" > "$3"
    ./rockim_e3.exe "$1" "$2" >> "$3" 2>&1
    rc=$?
    echo "=== TERMINE (rc=$rc) en $(( $(date +%s) - t0 )) s ===" >> "$3"
    return $rc
}

run bench_abuaisha/configs/e3_iso12.cfg  out_e3_iso12 bench_abuaisha/run_e3_iso.log   "essai 3 isotrope"
run bench_abuaisha/configs/e3_aniso.cfg  out_e3_aniso bench_abuaisha/run_e3_aniso.log "essai 3 anisotrope"
echo "=== CHAINE ESSAI 3 — protocole article : injection apres excavation (hydroStart).
