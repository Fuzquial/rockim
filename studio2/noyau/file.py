"""File de calculs persistante, exécutée en parallèle.

Reprend la logique de etude_triax_hetero/lancer_lot.py :
  * N jobs x M fils (mesure de LANCER.md : le 2D plafonne vers 8 fils par job,
    4 x 4 vaut 2 x 8) ;
  * un run dont le dossier a déjà son history.csv n'est PAS relancé ;
  * stdout du solveur redirigé vers un journal par run (le résumé 2D y est écrit).
Et ajoute ce qu'un script ne fait pas :
  * l'état de la file est écrit sur disque (file.json) à chaque changement ;
  * fermer l'application ne tue pas les calculs : à la réouverture, un run dont
    le processus vit encore est rattaché par son PID, sinon son état final est
    déduit de son dossier (ERROR.txt, résumé du journal) ;
  * relancer un run déjà fait RENOMME l'ancien dossier, il ne le supprime jamais.

La file n'a pas de fil d'exécution propre : l'appelant (le serveur de
l'interface) appelle pas() périodiquement. Chaque appel est court (pas de
lecture de gros fichiers), pour ne jamais bloquer l'interface.
"""
import ctypes
import json
import os
import signal
import subprocess
import sys
import time
import uuid

from .essai import Essai
from .resultats import dernier_temps

ETATS_FINAUX = ("fini", "echec", "arrete", "interrompu", "deja_fait")


def pid_vivant(pid):
    if not pid:
        return False
    if sys.platform == "win32":
        k = ctypes.windll.kernel32
        h = k.OpenProcess(0x1000, False, int(pid))      # PROCESS_QUERY_LIMITED_INFORMATION
        if not h:
            return False
        code = ctypes.c_ulong()
        ok = k.GetExitCodeProcess(h, ctypes.byref(code))
        k.CloseHandle(h)
        return bool(ok) and code.value == 259            # STILL_ACTIVE
    try:
        os.kill(int(pid), 0)
        return True
    except OSError:
        return False


class File:
    def __init__(self, dossier, commande, jobs=4, fils=4, cwd=None):
        """commande : liste de base, ex. ["C:/.../rockim_g1y19.exe"] ; le deck et le
        dossier de sortie sont ajoutés à la fin, comme dans lancer_lot.py."""
        self.dossier = os.path.abspath(dossier)
        self.commande = list(commande)
        self.jobs, self.fils = jobs, fils
        self.cwd = cwd
        self.procs = {}                                  # id -> Popen (runs lancés par CET objet)
        self.journaux = {}
        for d in ("decks", "out", "logs"):
            os.makedirs(os.path.join(self.dossier, d), exist_ok=True)
        self.chemin = os.path.join(self.dossier, "file.json")
        self.travaux = []
        if os.path.exists(self.chemin):
            with open(self.chemin, encoding="utf-8") as f:
                self.travaux = json.load(f)["travaux"]
            self._rattacher()

    # ------------------------------------------------------------ persistance
    def _sauver(self):
        tmp = self.chemin + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump({"version": 1, "travaux": self.travaux}, f, indent=1)
        os.replace(tmp, self.chemin)

    def _rattacher(self):
        for t in self.travaux:
            if t["etat"] == "en_cours" and not pid_vivant(t.get("pid")):
                t["etat"] = self._etat_deduit(t)
                t["fin"] = t.get("fin") or time.time()
        self._sauver()

    def _etat_deduit(self, t):
        """État final d'un run dont on n'a pas vu le code de retour."""
        if os.path.exists(os.path.join(t["out"], "ERROR.txt")):
            return "echec"
        try:
            log = open(t["log"], encoding="utf-8", errors="replace").read()
        except OSError:
            log = ""
        return "fini" if "---- summary ----" in log else "interrompu"

    # ------------------------------------------------------------ édition
    def trouver(self, ident):
        for t in self.travaux:
            if t["id"] == ident:
                return t
        raise KeyError(ident)

    def ajouter(self, essai: Essai):
        if any(t["nom"] == essai.nom for t in self.travaux if t["etat"] not in ETATS_FINAUX):
            raise ValueError("un travail « %s » est déjà dans la file" % essai.nom)
        deck = os.path.join(self.dossier, "decks", essai.nom + ".cfg")
        with open(deck, "w", encoding="utf-8", newline="\n") as f:
            f.write(essai.vers_cfg())
        t = dict(id=uuid.uuid4().hex[:8], nom=essai.nom, essai=essai.vers_dict(), deck=deck,
                 out=os.path.join(self.dossier, "out", essai.nom),
                 log=os.path.join(self.dossier, "logs", essai.nom + ".log"),
                 etat="en_attente", pid=None, code=None, ajout=time.time(), debut=None, fin=None)
        self.travaux.append(t)
        self._sauver()
        return t

    def retirer(self, ident):
        t = self.trouver(ident)
        if t["etat"] == "en_cours":
            raise ValueError("arrêter le run avant de le retirer")
        self.travaux.remove(t)
        self._sauver()

    def deplacer(self, ident, delta):
        """Change la place d'un travail dans la file (delta = -1 : plus tôt)."""
        t = self.trouver(ident)
        i = self.travaux.index(t)
        j = max(0, min(len(self.travaux) - 1, i + delta))
        self.travaux.insert(j, self.travaux.pop(i))
        self._sauver()

    def relancer(self, ident):
        """Remet un travail terminé en attente. L'ancien dossier est RENOMMÉ, jamais supprimé."""
        t = self.trouver(ident)
        if t["etat"] == "en_cours":
            raise ValueError("le run tourne encore")
        if os.path.exists(t["out"]):
            os.rename(t["out"], t["out"] + "__ancien_" + time.strftime("%Y%m%d-%H%M%S"))
        t.update(etat="en_attente", pid=None, code=None, debut=None, fin=None)
        self._sauver()

    # ------------------------------------------------------------ exécution
    def en_cours(self):
        return [t for t in self.travaux if t["etat"] == "en_cours"]

    def arreter(self, ident):
        t = self.trouver(ident)
        if t["etat"] == "en_attente":
            t["etat"] = "arrete"
        elif t["etat"] == "en_cours":
            p = self.procs.get(ident)
            if p is not None:
                p.terminate()
                p.wait(timeout=30)
            elif pid_vivant(t["pid"]):
                os.kill(int(t["pid"]), signal.SIGTERM)   # TerminateProcess sous Windows
            self._fermer_journal(ident)
            t.update(etat="arrete", fin=time.time())
        self._sauver()

    def _fermer_journal(self, ident):
        j = self.journaux.pop(ident, None)
        if j:
            j.close()
        self.procs.pop(ident, None)

    def _lancer(self, t):
        if os.path.exists(os.path.join(t["out"], "history.csv")):
            t.update(etat="deja_fait", fin=time.time())
            return
        env = dict(os.environ, OMP_NUM_THREADS=str(self.fils))
        log = open(t["log"], "w", encoding="utf-8", errors="replace")
        drapeaux = 0
        if sys.platform == "win32":
            # Nouveau groupe de processus et pas de console : le calcul survit à la
            # fermeture de l'interface et ne reçoit pas ses Ctrl+C.
            drapeaux = subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.CREATE_NO_WINDOW
        p = subprocess.Popen(self.commande + [t["deck"], t["out"]], cwd=self.cwd, env=env,
                             stdout=log, stderr=subprocess.STDOUT, creationflags=drapeaux)
        self.procs[t["id"]], self.journaux[t["id"]] = p, log
        t.update(etat="en_cours", pid=p.pid, debut=time.time())

    def pas(self):
        """Récolte les runs terminés et lance les suivants. Rend la liste des changements."""
        changes = []
        for t in self.en_cours():
            p = self.procs.get(t["id"])
            if p is not None:
                code = p.poll()
                if code is None:
                    continue
                self._fermer_journal(t["id"])
                t.update(etat="fini" if code == 0 else "echec", code=code, fin=time.time())
            elif not pid_vivant(t["pid"]):           # rattaché après réouverture
                t.update(etat=self._etat_deduit(t), fin=time.time())
            else:
                continue
            changes.append(t["id"])
        for t in self.travaux:
            if len(self.en_cours()) >= self.jobs:
                break
            if t["etat"] == "en_attente":
                self._lancer(t)
                changes.append(t["id"])
        if changes:
            self._sauver()
        return changes

    def avancement(self, ident):
        """Fraction du temps simulé atteinte (0-1), lue sur la dernière ligne d'historique."""
        t = self.trouver(ident)
        if t["etat"] in ("fini", "deja_fait"):
            return 1.0
        tf = dernier_temps(t["out"])
        T = t["essai"]["sorties"]["T"]
        return min(1.0, tf / T) if tf is not None and T > 0 else 0.0

    def attendre(self, delai=None, periode=0.2):
        """Fait tourner la file jusqu'à ce qu'elle soit vide (scripts et tests)."""
        t0 = time.time()
        while True:
            self.pas()
            if not any(t["etat"] in ("en_attente", "en_cours") for t in self.travaux):
                return
            if delai is not None and time.time() - t0 > delai:
                raise TimeoutError("file non terminée après %g s" % delai)
            time.sleep(periode)
