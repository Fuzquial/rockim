"""Prototype C (spec 007, jalon J1) : serveur local du prototype web.

Bibliothèque standard seulement. Sert l'interface (index.html, app.js, style.css)
et les caches binaires de ../_cache tels quels :
le navigateur reçoit les .bin bruts et les place directement en mémoire GPU.

    python serveur.py [port]          puis ouvrir http://localhost:8765
"""
import json
import os
import sys
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

ICI = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.normpath(os.path.join(ICI, "..", "_cache"))

TYPES = {".js": "text/javascript", ".css": "text/css", ".html": "text/html; charset=utf-8",
         ".json": "application/json", ".bin": "application/octet-stream"}


class Gestionnaire(SimpleHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def _envoyer(self, chemin):
        if not os.path.isfile(chemin):
            self.send_error(404)
            return
        with open(chemin, "rb") as f:
            corps = f.read()
        self.send_response(200)
        self.send_header("Content-Type", TYPES.get(os.path.splitext(chemin)[1], "application/octet-stream"))
        self.send_header("Content-Length", str(len(corps)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(corps)

    def do_GET(self):
        p = self.path.split("?")[0]
        if p == "/api/runs":
            runs = sorted(d for d in os.listdir(CACHE) if os.path.isfile(os.path.join(CACHE, d, "meta.json")))
            corps = json.dumps(runs).encode()
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(corps)))
            self.end_headers()
            self.wfile.write(corps)
        elif p.startswith("/cache/"):
            rel = os.path.normpath(p[len("/cache/"):])
            if rel.startswith(".."):
                self.send_error(403)
                return
            self._envoyer(os.path.join(CACHE, rel))
        else:
            self._envoyer(os.path.join(ICI, "static", os.path.basename(p) if p != "/" else "index.html"))

    def do_POST(self):
        # Le mode mesure du navigateur dépose ici ses chiffres (N1-N6).
        if self.path != "/api/mesures":
            self.send_error(404)
            return
        corps = self.rfile.read(int(self.headers.get("Content-Length", 0)))
        with open(os.path.join(ICI, "mesures_web.json"), "wb") as f:
            f.write(corps)
        self.send_response(204)
        self.end_headers()


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8765
    print("Prototype web : http://localhost:%d" % port)
    ThreadingHTTPServer(("127.0.0.1", port), Gestionnaire).serve_forever()
