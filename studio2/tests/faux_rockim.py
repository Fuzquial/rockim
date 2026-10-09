"""Faux solveur pour tester la file sans lancer g1 : `python faux_rockim.py <deck> <out>`.

Lit T dans le deck, écrit un history.csv ligne à ligne pendant FAUX_DUREE secondes
(défaut 0,6), imprime un résumé comme le solveur, et sort avec le code FAUX_CODE
(défaut 0). Note l'heure de début et de fin dans <out>/faux.json pour vérifier
le parallélisme.
"""
import json
import os
import sys
import time

deck, out = sys.argv[1], sys.argv[2]
T = 1e-2
for l in open(deck, encoding="utf-8"):
    l = l.split("#")[0]
    if l.strip().startswith("T ") or l.strip().startswith("T="):
        T = float(l.split("=")[1])
duree = float(os.environ.get("FAUX_DUREE", "0.6"))
os.makedirs(out, exist_ok=True)
debut = time.time()
print("[FDEM] 12995 elements, 19339 joints, 38985 nodes", flush=True)
with open(os.path.join(out, "history.csv"), "w", encoding="utf-8") as f:
    f.write("t,gripFy,sigma,epsGauge\n")
    n = 20
    for i in range(n + 1):
        f.write("%g,0,%g,%g\n" % (T * i / n, -1e6 * i, -1e-4 * i))
        f.flush()
        time.sleep(duree / n)
code = int(os.environ.get("FAUX_CODE", "0"))
if code == 0:
    print("[FDEM] ---- summary ----", flush=True)
json.dump({"debut": debut, "fin": time.time()}, open(os.path.join(out, "faux.json"), "w"))
sys.exit(code)
