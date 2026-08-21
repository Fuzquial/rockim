#!/bin/bash
# ESSAI 2 — suppression de l incubation a zone cohesive seche (hydroWetDamage = 0).
# Les deux cas EN SERIE, isotrope puis anisotrope : un seul job a la fois.
cd "C:/Users/fuzquianoalricabi/simulations/FDEM/rockim/rockim_p1"

run () {           # $1 = config, $2 = dossier de sortie, $3 = journal, $4 = libelle
    t0=$(date +%s)
    echo "=== $4 : debut $(date +%H:%M:%S) ===" > "$3"
    ./rockim_e2.exe "$1" "$2" >> "$3" 2>&1
    rc=$?
    echo "=== TERMINE (rc=$rc) en $(( $(date +%s) - t0 )) s ===" >> "$3"
    return $rc
}

run bench_abuaisha/configs/e2_iso12.cfg  out_e2_iso12 bench_abuaisha/run_e2_iso.log   "essai 2 isotrope"
run bench_abuaisha/configs/e2_aniso.cfg  out_e2_aniso bench_abuaisha/run_e2_aniso.log "essai 2 anisotrope"
echo "=== CHAINE ESSAI 2 TERMINEE $(date +%H:%M:%S) ==="
