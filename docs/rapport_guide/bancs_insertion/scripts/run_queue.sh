#!/bin/bash
# run_queue.sh liste_de_decks.txt : lance chaque deck (chemin absolu) avec au plus 2 calculs
# en parallele (NPAR, defaut 2), OMP_NUM_THREADS=2, nice -n 10, 20 min max (timeout 1200 s). Le dossier de
# sortie est <deck sans .cfg>_out ; le journal <deck sans .cfg>.log ; la duree dans .dur.
BIN=/home/user/rockim/build_contact/rockim
run_one() {
  d="$1"; b="${d%.cfg}"
  s=$(date +%s.%N)
  OMP_NUM_THREADS=2 timeout 1200 nice -n 10 $BIN "$d" "${b}_out" > "${b}.log" 2>&1
  rc=$?
  e=$(date +%s.%N)
  echo "rc=$rc dur=$(echo "$e - $s" | bc)" > "${b}.dur"
}
export -f run_one; export BIN
cat "$1" | xargs -P ${NPAR:-2} -I{} bash -c 'run_one {}'
