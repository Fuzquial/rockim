"""Serveur local de l'interface (spec 007, J4) : le noyau exposé en API JSON.

Bibliothèque standard seulement. Une fois lancé :

    python studio2/app/serveur.py [--espace DOSSIER] [--port 8770]
    puis ouvrir http://localhost:8770

L'espace de travail contient la file (file.json, decks/, out/, logs/), les caches
d'affichage (caches/) et les réglages (reglages.json). Les runs des dossiers de
campagne déclarés dans les réglages (« bibliothèques ») sont lus sur place, jamais
modifiés : leurs caches vont dans l'espace.

Un fil de fond fait avancer la file toutes les secondes. Les conversions de runs en
cache tournent dans un autre fil, une à la fois. Aucun appel d'API ne fait de
travail long sur le fil qui répond.
"""
import argparse
import copy
import json
import os
import re
import sys
import threading
import time
import traceback
from concurrent.futures import ThreadPoolExecutor
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, unquote, urlparse

ICI = os.path.dirname(os.path.abspath(__file__))
STUDIO2 = os.path.dirname(ICI)
G1 = os.path.dirname(STUDIO2)
sys.path.insert(0, STUDIO2)

from noyau import depouillement, formulaire, maillage, materiaux, resultats, validation   # noqa: E402
from noyau.essai import Essai                                         # noqa: E402
from noyau.file import ETATS_FINAUX, File                             # noqa: E402

REGLAGES_DEFAUT = {
    "exe": "rockim_j3.exe",
    "pause": False,
    "jobs": 4,
    "fils": 4,
    "bibliotheques": [{"nom": "etude_triax_hetero",
                       "out": os.path.join(G1, "etude_triax_hetero", "out"),
                       "logs": os.path.join(G1, "etude_triax_hetero", "logs")}],
}
SEP = "~"                                    # identifiant de run : <source>~<nom>


class Studio:
    def __init__(self, espace, commande=None, periode=1.0):
        self.espace = os.path.abspath(espace)
        os.makedirs(os.path.join(self.espace, "caches"), exist_ok=True)
        self.chemin_reglages = os.path.join(self.espace, "reglages.json")
        self.reglages = copy.deepcopy(REGLAGES_DEFAUT)
        if os.path.exists(self.chemin_reglages):
            self.reglages.update(json.load(open(self.chemin_reglages, encoding="utf-8")))
        self.commande_forcee = commande          # tests : faux solveur
        self.verrou = threading.RLock()
        self.file = File(self.espace, self.commande(), self.reglages["jobs"], self.reglages["fils"], cwd=G1)
        self.conversions = {}                    # id -> "en_cours" | "erreur: ..."
        self.pool = ThreadPoolExecutor(max_workers=1)
        self.arret = threading.Event()
        self.periode = periode
        self.fil = threading.Thread(target=self._boucle, daemon=True)
        self.fil.start()

    # ------------------------------------------------------------ réglages
    def commande(self):
        if self.commande_forcee:
            return list(self.commande_forcee)
        exe = self.reglages["exe"]
        return [exe if os.path.isabs(exe) else os.path.join(G1, exe)]

    def sauver_reglages(self, nouveaux):
        with self.verrou:
            for k in ("exe", "jobs", "fils", "bibliotheques", "pause"):
                if k in nouveaux:
                    self.reglages[k] = nouveaux[k]
            json.dump(self.reglages, open(self.chemin_reglages, "w", encoding="utf-8"), indent=1)
            self.file.commande = self.commande()
            self.file.jobs, self.file.fils = int(self.reglages["jobs"]), int(self.reglages["fils"])
        return self.reglages

    def exes(self):
        return sorted((f for f in os.listdir(G1) if re.fullmatch(r"rockim_.*\.exe", f)),
                      key=lambda f: os.path.getmtime(os.path.join(G1, f)), reverse=True)

    # ------------------------------------------------------------ boucle de fond
    def _boucle(self):
        while not self.arret.is_set():
            try:
                with self.verrou:
                    self.file.pas(lancer=not self.reglages.get("pause"))
            except Exception:
                traceback.print_exc()
            self.arret.wait(self.periode)

    def fermer(self):
        self.arret.set()
        self.fil.join(timeout=5)
        self.pool.shutdown(wait=False)

    # ------------------------------------------------------------ runs
    def runs(self):
        """Tous les runs visibles : ceux de l'espace, puis ceux des bibliothèques."""
        out = []
        with self.verrou:
            for t in self.file.travaux:
                out.append(dict(id="espace" + SEP + t["nom"], nom=t["nom"], source="espace",
                                etat=t["etat"], dossier=t["out"], log=t["log"],
                                sigma3_MPa=t["essai"]["chargement"]["sigma3_MPa"]))
        for i, b in enumerate(self.reglages["bibliotheques"]):
            if not os.path.isdir(b["out"]):
                continue
            for nom in sorted(os.listdir(b["out"]), key=tri_naturel):
                d = os.path.join(b["out"], nom)
                if nom.startswith("_") or not os.path.exists(os.path.join(d, "history.csv")):
                    continue
                # Un run arrêté avant sa première ligne d'historique n'a rien à montrer.
                vide = resultats.dernier_temps(d) is None
                out.append(dict(id="b%d" % i + SEP + nom, nom=nom, source=b["nom"], etat="vide" if vide else "fini",
                                dossier=d, log=os.path.join(b["logs"], nom + ".log"), sigma3_MPa=None))
        return out

    def run(self, ident):
        for r in self.runs():
            if r["id"] == ident:
                return r
        raise KeyError(ident)

    def synthese(self, ident):
        r = self.run(ident)
        m = depouillement.synthese(r["dossier"], r["log"])
        return propre(m)

    def historique(self, ident):
        r = self.run(ident)
        h = resultats.lire_historique(r["dossier"])
        if h is None:
            return {"eps": [], "q": [], "t": []}
        e, q = resultats.courbe(h, resultats.delai_consolidation(r["dossier"]))
        return {"eps": propre(list(e)), "q": propre(list(q)), "t": list(h["t"] * 1e3)}

    # ------------------------------------------------------------ caches d'affichage
    def dossier_cache(self, ident):
        return os.path.join(self.espace, "caches", ident.replace(SEP, "__"))

    def etat_cache(self, ident):
        r = self.run(ident)
        c = self.dossier_cache(ident)
        meta = os.path.join(c, "meta.json")
        if self.conversions.get(ident) == "en_cours":
            return {"etat": "en_cours"}
        if os.path.exists(meta):
            frames = resultats.frames_disponibles(r["dossier"])
            m = json.load(open(meta, encoding="utf-8"))
            if m["nFrames"] == len(frames):
                return {"etat": "pret", "meta": m}
        if str(self.conversions.get(ident, "")).startswith("erreur"):
            return {"etat": self.conversions[ident]}
        return {"etat": "absent"}

    def demander_cache(self, ident):
        e = self.etat_cache(ident)
        if e["etat"] in ("pret", "en_cours"):
            return e
        r = self.run(ident)
        if not resultats.frames_disponibles(r["dossier"]):
            return {"etat": "erreur: aucune frame VTU"}
        self.conversions[ident] = "en_cours"

        def tache():
            try:
                resultats.convertir(r["dossier"], self.dossier_cache(ident))
                self.conversions.pop(ident, None)
            except Exception as ex:
                self.conversions[ident] = "erreur: %s" % ex
        self.pool.submit(tache)
        return {"etat": "en_cours"}

    # ------------------------------------------------------------ file
    def etat_file(self):
        with self.verrou:
            out = []
            for t in self.file.travaux:
                d = {k: t[k] for k in ("id", "nom", "etat", "code", "debut", "fin", "ajout")}
                d["avancement"] = self.file.avancement(t["id"]) if t["etat"] in ("en_cours", "fini", "deja_fait") else 0.0
                d["sigma3_MPa"] = t["essai"]["chargement"]["sigma3_MPa"]
                d["resume"] = resume_essai(t["essai"])
                d["lot"] = t.get("lot")
                out.append(d)
            return {"travaux": out, "jobs": self.file.jobs, "fils": self.file.fils,
                    "exe": os.path.basename(self.file.commande[-1]), "pause": bool(self.reglages.get("pause"))}

    def action_file(self, ident, action):
        with self.verrou:
            if action == "arreter":
                self.file.arreter(ident)
            elif action == "relancer":
                self.file.relancer(ident)
            elif action == "retirer":
                self.file.retirer(ident)
            elif action == "monter":
                self.file.deplacer(ident, -1)
            elif action == "descendre":
                self.file.deplacer(ident, +1)
            else:
                raise ValueError("action inconnue : %s" % action)
        return self.etat_file()

    def journal(self, ident, n=200):
        with self.verrou:
            t = self.file.trouver(ident)
        try:
            lignes = open(t["log"], encoding="utf-8", errors="replace").read().splitlines()
        except OSError:
            lignes = []
        return {"lignes": lignes[-n:], "total": len(lignes)}

    # ------------------------------------------------------------ essais
    def verifier(self, d):
        e = Essai.depuis_dict(d)
        est = maillage.estimation_voronoi(e) if e.maillage.type == "voronoi" else {}
        return {"avis": [a.__dict__ for a in validation.verifier(e)], "estimation": est, "deck": e.vers_cfg()}

    def construire(self, choix):
        """Choix de l'écran Essai -> essai, vérifications, estimations, deck."""
        e = formulaire.essai_depuis_choix(choix)
        est = maillage.estimation_voronoi(e) if e.maillage.type == "voronoi" else {}
        cout = formulaire.estimation_cout(e, fils=int(self.reglages["fils"]))
        return {"essai": e.vers_dict(), "avis": [a.__dict__ for a in validation.verifier(e)],
                "estimation": dict(est, **cout), "deck": e.vers_cfg(),
                "segments_prerompus": [list(x) for x in e.segments_prerompus()]}

    def formulaire(self):
        """Valeurs par défaut et préréglages, pour construire l'écran."""
        niveaux = {}
        for n in materiaux.NIVEAUX:
            m = materiaux.materiau(n)
            niveaux[n] = {"materiau": m.__dict__, "phases": [p.__dict__ for p in materiaux.phases(n)],
                          "alpha": materiaux.NIVEAUX[n]["gbAlpha"]}
        return {"defaut": formulaire.CHOIX_DEFAUT, "niveaux": niveaux}

    def ajouter(self, d, sigma3_liste=None, lot=None):
        e = Essai.depuis_dict(d)
        if validation.erreurs(validation.verifier(e)):
            raise ValueError("l'essai a des erreurs : corriger avant d'ajouter")
        essais = [e]
        if sigma3_liste:
            essais = []
            for s in sigma3_liste:
                x = copy.deepcopy(e)
                x.chargement.sigma3_MPa = float(s)
                x.nom = "%s_P%03d" % (e.nom, int(round(float(s))))
                essais.append(x)
        ajoutes = []
        with self.verrou:
            for x in essais:
                t = self.file.ajouter(x)
                if lot:
                    t["lot"] = lot
                ajoutes.append(t["id"])
            self.file._sauver()
        return {"ajoutes": ajoutes}


def tri_naturel(nom):
    """F2 avant F10 : les nombres sont comparés comme des nombres."""
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", nom)]


def resume_essai(d):
    m, c = d["maillage"], d["chargement"]
    s = "Voronoï" if m["type"] == "voronoi" else "Gmsh"
    s += ", %d phase%s" % (max(1, len(d["phases"])), "s" if len(d["phases"]) > 1 else "") if m["type"] == "voronoi" else ""
    s += ", σ₃ = %g MPa" % c["sigma3_MPa"] if c["type_essai"] == "triaxial" else ", traction directe"
    if d["discontinuites"]["fraction_diffuse"]:
        s += ", %g %% pré-rompus" % (100 * d["discontinuites"]["fraction_diffuse"])
    return s


def propre(x):
    """JSON strict : NaN et infinis deviennent null."""
    if isinstance(x, dict):
        return {k: propre(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [propre(v) for v in x]
    try:
        f = float(x)
        if f != f or f in (float("inf"), float("-inf")):
            return None
        return x if isinstance(x, (int, bool)) else f
    except (TypeError, ValueError):
        return x


# ---------------------------------------------------------------- HTTP
TYPES = {".js": "text/javascript; charset=utf-8", ".css": "text/css; charset=utf-8",
         ".html": "text/html; charset=utf-8", ".json": "application/json", ".bin": "application/octet-stream",
         ".svg": "image/svg+xml", ".png": "image/png"}


def fabrique(studio):
    class Gestionnaire(SimpleHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def _json(self, obj, code=200):
            corps = json.dumps(obj, ensure_ascii=False).encode("utf-8")
            self.send_response(code)
            self.send_header("Content-Type", "application/json; charset=utf-8")
            self.send_header("Content-Length", str(len(corps)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(corps)

        def _fichier(self, chemin):
            if not os.path.isfile(chemin):
                return self._json({"erreur": "introuvable"}, 404)
            corps = open(chemin, "rb").read()
            self.send_response(200)
            self.send_header("Content-Type", TYPES.get(os.path.splitext(chemin)[1], "application/octet-stream"))
            self.send_header("Content-Length", str(len(corps)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(corps)

        def _corps(self):
            n = int(self.headers.get("Content-Length", 0))
            return json.loads(self.rfile.read(n) or b"{}")

        def _route(self, methode):
            u = urlparse(self.path)
            p = [unquote(x) for x in u.path.strip("/").split("/") if x]
            q = parse_qs(u.query)
            try:
                if methode == "GET" and (not p or p[0] in ("js", "css") or p[-1].endswith((".html", ".js", ".css", ".svg", ".png"))):
                    rel = os.path.normpath(os.path.join(*p)) if p else "index.html"
                    if rel.startswith(".."):
                        return self._json({"erreur": "interdit"}, 403)
                    return self._fichier(os.path.join(ICI, "static", rel))
                if p[0] == "cache" and len(p) == 3:
                    return self._fichier(os.path.join(studio.dossier_cache(p[1]), os.path.basename(p[2])))
                if p[0] != "api":
                    return self._json({"erreur": "introuvable"}, 404)
                a = p[1:]
                if methode == "GET":
                    if a == ["etat"]:
                        return self._json({"reglages": studio.reglages, "exes": studio.exes()})
                    if a == ["runs"]:
                        return self._json([{k: r[k] for k in ("id", "nom", "source", "etat", "sigma3_MPa")} for r in studio.runs()])
                    if len(a) == 3 and a[0] == "runs" and a[2] == "synthese":
                        return self._json(studio.synthese(a[1]))
                    if len(a) == 3 and a[0] == "runs" and a[2] == "historique":
                        return self._json(studio.historique(a[1]))
                    if len(a) == 3 and a[0] == "runs" and a[2] == "cache":
                        return self._json(studio.etat_cache(a[1]))
                    if a == ["file"]:
                        return self._json(studio.etat_file())
                    if a == ["formulaire"]:
                        return self._json(propre(studio.formulaire()))
                    if len(a) == 3 and a[0] == "file" and a[2] == "journal":
                        return self._json(studio.journal(a[1], int(q.get("n", ["200"])[0])))
                if methode == "POST":
                    if len(a) == 3 and a[0] == "runs" and a[2] == "cache":
                        return self._json(studio.demander_cache(a[1]))
                    if len(a) == 3 and a[0] == "file":
                        return self._json(studio.action_file(a[1], a[2]))
                    if a == ["reglages"]:
                        return self._json(studio.sauver_reglages(self._corps()))
                    if a == ["essai", "construire"]:
                        return self._json(propre(studio.construire(self._corps())))
                    if a == ["essai", "verifier"]:
                        return self._json(studio.verifier(self._corps()))
                    if a == ["essai", "ajouter"]:
                        c = self._corps()
                        return self._json(studio.ajouter(c["essai"], c.get("sigma3_liste"), c.get("lot")))
                return self._json({"erreur": "route inconnue : %s %s" % (methode, u.path)}, 404)
            except KeyError as ex:
                return self._json({"erreur": "introuvable : %s" % ex}, 404)
            except ValueError as ex:
                return self._json({"erreur": str(ex)}, 400)
            except Exception as ex:
                traceback.print_exc()
                return self._json({"erreur": "%s: %s" % (type(ex).__name__, ex)}, 500)

        def do_GET(self):
            self._route("GET")

        def do_POST(self):
            self._route("POST")

    return Gestionnaire


def demarrer(espace, port=8770, commande=None, periode=1.0):
    studio = Studio(espace, commande=commande, periode=periode)
    serveur = ThreadingHTTPServer(("127.0.0.1", port), fabrique(studio))
    serveur.daemon_threads = True
    return studio, serveur


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--espace", default=os.path.join(STUDIO2, "espace"))
    ap.add_argument("--port", type=int, default=8770)
    ap.add_argument("--faux-solveur", action="store_true",
                    help="démonstration et tests : la file lance tests/faux_rockim.py au lieu de g1")
    a = ap.parse_args()
    faux = [sys.executable, os.path.join(STUDIO2, "tests", "faux_rockim.py")] if a.faux_solveur else None
    studio, serveur = demarrer(a.espace, a.port, commande=faux)
    print("Rockim : http://localhost:%d  (espace %s)" % (a.port, studio.espace), flush=True)
    try:
        serveur.serve_forever()
    finally:
        studio.fermer()
