"""Prototype B (jetable) de la spec 007 : écran « Résultats » d'un triaxial 2D en PySide6 + Qt Quick.

    python main.py [run] [--mesure] [--sansvsync]

Rendu des champs par QtQuick3D (QQuick3DGeometry, caméra orthographique) : les sommets sont
envoyés au GPU en un seul tampon entrelacé (position xyz + couleur rgba, float32). Changer de
frame ou de champ = recalculer ce tampon en numpy vectorisé et le remplacer.

--mesure : déroule seul N1, N3, N4, N5, N6 puis écrit mesures_qml.json et capture_qml.png.
"""
import time

T0 = time.perf_counter()  # N1 : origine des temps, avant tout import lourd

import json
import os
import sys
import threading

import numpy as np
import PySide6

# Python >= 3.8 ne cherche plus les DLL dans PATH : sans cela les greffons QML (qtquick2plugin.dll,
# etc.) ne trouvent pas leurs dépendances Qt6*.dll rangées dans le dossier de PySide6.
_DLL = os.add_dll_directory(os.path.dirname(PySide6.__file__)) if os.name == "nt" else None

from PySide6.QtCore import (QByteArray, QObject, QPointF, Qt, QTimer, QUrl, Property, Signal,
                            Slot)
from PySide6.QtGui import QColor, QGuiApplication, QSurfaceFormat, QVector3D
from PySide6.QtQml import QmlElement, QQmlApplicationEngine
from PySide6.QtQuick import QQuickWindow
from PySide6.QtQuick3D import QQuick3DGeometry
from PySide6.QtGraphs import QLineSeries  # noqa: F401  (type Python de la série reçue du QML)
from PySide6.QtQuickControls2 import QQuickStyle

T_IMPORTS = time.perf_counter()  # fin des imports (N1 détaillé)
ICI = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.normpath(os.path.join(ICI, "..", "_cache"))
CHAMPS = ["sigmaXX", "sigmaYY", "sigmaXY", "vonMises", "epsXX", "phase"]
NB_BANDES = 12          # palette à bandes discrètes, comme Abaqus
MM = 1000.0             # positions du cache en m, affichées en mm

QML_IMPORT_NAME = "Proto"
QML_IMPORT_MAJOR_VERSION = 1


# --------------------------------------------------------------------------------------------
# Géométries GPU
# --------------------------------------------------------------------------------------------
class _GeometrieCouleur(QQuick3DGeometry):
    """Tampon entrelacé position (3 float) + couleur (4 float), 28 octets par sommet."""

    STRIDE = 7 * 4

    def __init__(self, primitive, parent=None):
        super().__init__(parent)
        self.setStride(self.STRIDE)
        self.setPrimitiveType(primitive)
        A = QQuick3DGeometry.Attribute
        self.addAttribute(A.Semantic.PositionSemantic, 0, A.ComponentType.F32Type)
        self.addAttribute(A.Semantic.ColorSemantic, 12, A.ComponentType.F32Type)

    def pousser(self, tampon, bmin, bmax):
        """tampon : ndarray float32 [n, 7] contigu."""
        self.setVertexData(QByteArray(tampon.tobytes()))
        self.setBounds(QVector3D(bmin[0], bmin[1], -1.0), QVector3D(bmax[0], bmax[1], 1.0))
        self.update()


@QmlElement
class GeometrieChamps(_GeometrieCouleur):
    def __init__(self, parent=None):
        super().__init__(QQuick3DGeometry.PrimitiveType.Triangles, parent)


@QmlElement
class GeometrieJoints(QQuick3DGeometry):
    """Rubans des joints : position (3) + autre extrémité (2) + côté (2) + couleur (4) = 44 octets.
    L'épaisseur en pixels est appliquée dans joints.vert (CustomMaterial) : les lignes D3D11 ne
    dépassent pas 1 px."""

    STRIDE = 11 * 4

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setStride(self.STRIDE)
        self.setPrimitiveType(QQuick3DGeometry.PrimitiveType.Triangles)
        A = QQuick3DGeometry.Attribute
        self.addAttribute(A.Semantic.PositionSemantic, 0, A.ComponentType.F32Type)
        self.addAttribute(A.Semantic.TexCoord0Semantic, 12, A.ComponentType.F32Type)
        self.addAttribute(A.Semantic.TexCoord1Semantic, 20, A.ComponentType.F32Type)
        self.addAttribute(A.Semantic.ColorSemantic, 28, A.ComponentType.F32Type)

    pousser = _GeometrieCouleur.pousser


# --------------------------------------------------------------------------------------------
# Données d'un run (cache binaire)
# --------------------------------------------------------------------------------------------
class Run:
    def __init__(self, nom):
        d = os.path.join(CACHE, nom)
        with open(os.path.join(d, "meta.json"), encoding="utf-8") as f:
            m = json.load(f)
        self.nom, self.meta = nom, m
        F, nV, nT, nJ = m["nFrames"], m["nVert"], m["nTri"], m["nJoint"]
        lire = lambda n, t: np.fromfile(os.path.join(d, n), dtype=t)
        self.pts = lire("pts.bin", "<f4").reshape(F, nV, 2) * np.float32(MM)
        self.champs = {c: lire(c + ".bin", "<f4").reshape(F, nT) for c in CHAMPS[:-1]}
        self.champs["phase"] = np.broadcast_to(lire("phase.bin", "<f4"), (F, nT))
        self.jseg = lire("jseg.bin", "<u4").reshape(nJ, 2)
        self.jmode = lire("jmode.bin", "u1").reshape(F, nJ)
        self.hist = lire("hist.bin", "<f4").reshape(-1, 2)
        self.temps = np.asarray(m["temps"])
        self.f2h = np.asarray(m["frameVersHist"], dtype=int)
        self.bmin = self.pts.min(axis=(0, 1))
        self.bmax = self.pts.max(axis=(0, 1))
        # Bornes d'affichage : MPa pour les contraintes, % pour epsXX, entiers pour la phase.
        b = dict(m["bornes"])
        ph = self.champs["phase"][0]
        b["phase"] = [float(ph.min()) - 0.5, float(ph.max()) + 0.5]
        self.bornes = {c: [v * self.echelle(c) for v in b[c]] for c in CHAMPS}

    @staticmethod
    def echelle(champ):
        return 1e-6 if champ.startswith("sigma") or champ == "vonMises" else (
            100.0 if champ == "epsXX" else 1.0)

    @staticmethod
    def unite(champ):
        return {"epsXX": "%", "phase": ""}.get(champ, "MPa")


def palette_bandes(stops, n):
    """Arc-en-ciel à n bandes interpolé entre les couleurs d'ancrage (bleu -> rouge)."""
    s = np.array([[c.redF(), c.greenF(), c.blueF()] for c in stops], dtype=np.float32)
    if len(s) == n:
        return np.hstack([s, np.ones((n, 1), np.float32)])
    x = np.linspace(0.0, 1.0, len(s))
    t = (np.arange(n) + 0.5) / n
    rgb = np.stack([np.interp(t, x, s[:, k]) for k in range(3)], axis=1)
    return np.hstack([rgb, np.ones((n, 1), np.float32)]).astype(np.float32)


# --------------------------------------------------------------------------------------------
# Contrôleur exposé au QML
# --------------------------------------------------------------------------------------------
class Studio(QObject):
    runChange = Signal()
    frameChange = Signal()
    champChange = Signal()
    modesChange = Signal()
    imageAffichee = Signal()

    def __init__(self, run_initial):
        super().__init__()
        self._runs = sorted(d for d in os.listdir(CACHE)
                            if os.path.isfile(os.path.join(CACHE, d, "meta.json")))
        self._nomInitial = run_initial
        self.run = None
        self._frame = 0
        self._champ = "sigmaYY"
        self._modes = {1: True, 2: True, 4: True}
        self.geoC = self.geoJ = self.serie = None
        self.lut = None
        self.coulModes = {}
        self.dureeChargement = None
        self._runLu.connect(self._installer, Qt.QueuedConnection)

    # -- branchement depuis le QML (Component.onCompleted) --
    @Slot(QObject, QObject, QObject, "QVariantList", QColor, QColor, QColor)
    def demarrer(self, geoC, geoJ, serie, arcEnCiel, cTraction, cCisail, cPreRompu):
        self.geoC, self.geoJ, self.serie = geoC, geoJ, serie
        self.lut = palette_bandes([QColor(c) for c in arcEnCiel], NB_BANDES)
        self.coulModes = {m: np.array([c.redF(), c.greenF(), c.blueF(), 1.0], np.float32)
                          for m, c in ((1, cTraction), (2, cCisail), (4, cPreRompu))}
        self.chargerRun(self._nomInitial, True)

    @Slot(str)
    def chargerRun(self, nom, synchrone=False):
        """Lecture des .bin dans un fil secondaire (numpy relâche le GIL pendant les E/S),
        puis installation des tampons GPU sur le fil UI. Au démarrage : synchrone, pour que
        la première image montre déjà le run."""
        self._t0Chargement = time.perf_counter()
        if synchrone:
            self._installer(Run(nom))
            return
        threading.Thread(target=lambda: self._runLu.emit(Run(nom)), daemon=True).start()

    _runLu = Signal(object)

    def _installer(self, r):
        self.run = r
        self._frame = min(self._frame, r.meta["nFrames"] - 1)
        # tampon de positions + couleurs alloué une fois par run
        self._bufC = np.zeros((r.meta["nVert"], 7), np.float32)
        self._remplirCourbe()
        self._majTampons()
        self.dureeChargement = time.perf_counter() - self._t0Chargement
        self.runChange.emit()
        self.frameChange.emit()
        self.champChange.emit()
        self.modesChange.emit()

    def _remplirCourbe(self):
        h = self.run.hist
        self.serie.replace([QPointF(float(x), float(y)) for x, y in h])

    # -- mise à jour GPU : tout en numpy vectorisé --
    def _majTampons(self):
        r, f = self.run, self._frame
        lo, hi = r.bornes[self._champ]
        v = r.champs[self._champ][f] * np.float32(r.echelle(self._champ))
        idx = np.clip(((v - lo) / (hi - lo) * NB_BANDES).astype(np.int32), 0, NB_BANDES - 1)
        b = self._bufC
        b[:, 0:2] = r.pts[f]
        b[:, 3:7] = np.repeat(self.lut[idx], 3, axis=0)   # sommet k -> triangle k // 3
        self.geoC.pousser(b, r.bmin, r.bmax)

        modes = r.jmode[f]
        sel = np.zeros(modes.shape, bool)
        for m, vis in self._modes.items():
            if vis:
                sel |= modes == m
        seg = r.jseg[sel]
        n = len(seg)
        bj = np.zeros((max(6 * n, 3), 11), np.float32)
        if n:
            P = r.pts[f]
            a, b = P[seg[:, 0]], P[seg[:, 1]]
            # 6 sommets par segment : (A-, A+, B-) et (B-, A+, B+) ; le côté est inversé
            # aux sommets B, car le vertex shader y mesure la direction B -> A.
            propre = np.stack([a, a, b, b, a, b], axis=1)
            autre = np.stack([b, b, a, a, b, a], axis=1)
            cote = np.array([-1, 1, 1, 1, 1, -1], np.float32)
            v = bj[:6 * n].reshape(n, 6, 11)
            v[:, :, 0:2] = propre
            v[:, :, 2] = 0.05                        # devant les triangles
            v[:, :, 3:5] = autre
            v[:, :, 5] = cote
            table = np.zeros((5, 4), np.float32)
            for m, c in self.coulModes.items():
                table[m] = c
            v[:, :, 7:11] = table[modes[sel]][:, None, :]
        self.geoJ.pousser(bj, r.bmin, r.bmax)

    # -- propriétés --
    @Property("QStringList", constant=True)
    def runs(self):
        return self._runs

    @Property(str, notify=runChange)
    def nomRun(self):
        return self.run.nom if self.run else ""

    @Property(int, notify=runChange)
    def nFrames(self):
        return self.run.meta["nFrames"] if self.run else 1

    @Property(float, notify=runChange)
    def sigma3(self):
        return float(self.run.meta.get("sigma3_MPa", 0)) if self.run else 0.0

    @Property(int, notify=runChange)
    def nTri(self):
        return self.run.meta["nTri"] if self.run else 0

    @Property(int, notify=runChange)
    def nJoint(self):
        return self.run.meta["nJoint"] if self.run else 0

    @Property("QVariantList", notify=runChange)
    def boite(self):
        """[xmin, ymin, xmax, ymax] en mm sur toutes les frames."""
        if not self.run:
            return [0, 0, 1, 1]
        return [float(self.run.bmin[0]), float(self.run.bmin[1]),
                float(self.run.bmax[0]), float(self.run.bmax[1])]

    @Property("QVariantList", notify=runChange)
    def courbeBornes(self):
        """[epsMax, qMin, qMax] pour les axes de la courbe."""
        if not self.run:
            return [1, 0, 1]
        h = self.run.hist
        return [float(h[:, 0].max()), float(min(0.0, h[:, 1].min())), float(h[:, 1].max())]

    @Property(float, notify=runChange)
    def qPic(self):
        return float(self.run.hist[:, 1].max()) if self.run else 0.0

    @Property(float, notify=runChange)
    def epsPic(self):
        return float(self.run.hist[np.argmax(self.run.hist[:, 1]), 0]) if self.run else 0.0

    def _get_frame(self):
        return self._frame

    def _set_frame(self, f):
        f = int(max(0, min(f, self.nFrames - 1)))
        if f != self._frame:
            self._frame = f
            self._majTampons()
            self.frameChange.emit()

    frame = Property(int, _get_frame, _set_frame, notify=frameChange)

    @Property(float, notify=frameChange)
    def tempsMs(self):
        return float(self.run.temps[self._frame] * 1e3) if self.run else 0.0

    @Property(float, notify=frameChange)
    def epsCourant(self):
        return float(self.run.hist[self.run.f2h[self._frame], 0]) if self.run else 0.0

    @Property(float, notify=frameChange)
    def qCourant(self):
        return float(self.run.hist[self.run.f2h[self._frame], 1]) if self.run else 0.0

    @Property("QVariantList", notify=frameChange)
    def comptesModes(self):
        """Nombre de joints [traction, cisaillement, pré-rompu] à la frame courante."""
        if not self.run:
            return [0, 0, 0]
        c = np.bincount(self.run.jmode[self._frame], minlength=5)
        return [int(c[1]), int(c[2]), int(c[4])]

    def _get_champ(self):
        return self._champ

    def _set_champ(self, c):
        if c in CHAMPS and c != self._champ:
            self._champ = c
            self._majTampons()
            self.champChange.emit()

    champ = Property(str, _get_champ, _set_champ, notify=champChange)

    @Property("QVariantList", notify=champChange)
    def echelleCouleurs(self):
        """[min, max, unité] du champ courant."""
        lo, hi = self.run.bornes[self._champ] if self.run else (0, 1)
        return [lo, hi, Run.unite(self._champ)]

    @Property(int, constant=True)
    def nbBandes(self):
        return NB_BANDES

    @Property("QVariantList", notify=modesChange)
    def modesVisibles(self):
        return [self._modes[1], self._modes[2], self._modes[4]]

    @Slot(int, bool)
    def montrerMode(self, mode, visible):
        self._modes[mode] = visible
        self._majTampons()
        self.modesChange.emit()

    @Slot(float, result=int)
    def frameProche(self, eps):
        """Frame dont le point d'historique est le plus proche de eps (axe x de la courbe)."""
        e = self.run.hist[self.run.f2h, 0]
        return int(np.argmin(np.abs(e - eps)))


# --------------------------------------------------------------------------------------------
# Instrumentation (--mesure)
# --------------------------------------------------------------------------------------------
class Mesure(QObject):
    """Automate piloté par frameSwapped : N1, N3, N4, N5, N6, capture."""

    def __init__(self, app, fenetre, studio, racine, sansvsync):
        super().__init__()
        self.app, self.w, self.s, self.racine = app, fenetre, studio, racine
        self.res = {"vsync": not sansvsync, "dpr": fenetre.devicePixelRatio()}
        self.etape = "N1"
        self.attente = None          # (nom, t0) : mesure en cours jusqu'à la prochaine image
        self.synchro = False         # vrai si une synchronisation a eu lieu après le changement
        self.n4 = []
        self.phaseN6 = "repos_apres_demarrage"
        self.n6 = {}
        self.tic = None
        self._n3Attente = None
        studio.runChange.connect(self._runInstalle)
        fenetre.frameSwapped.connect(self._image, Qt.QueuedConnection)
        fenetre.afterSynchronizing.connect(self._synchro, Qt.DirectConnection)
        # N6 : sonde sur le fil UI, l'écart entre deux ticks à 1 ms = plus long blocage
        self.sonde = QTimer(self)
        self.sonde.setTimerType(Qt.PreciseTimer)
        self.sonde.setInterval(1)
        self.sonde.timeout.connect(self._tic)

    def _tic(self):
        t = time.perf_counter()
        if self.tic is not None:
            dt = (t - self.tic) * 1e3
            self.n6[self.phaseN6] = max(self.n6.get(self.phaseN6, 0.0), dt)
        self.tic = t

    def _runInstalle(self):
        if self._n3Attente:
            self.synchro = False
            self.attente, self._n3Attente = self._n3Attente, None

    def _synchro(self):
        if self.attente:
            self.synchro = True

    def _armer(self, nom):
        self.synchro = False
        self.attente = (nom, time.perf_counter())

    def _image(self):
        t = time.perf_counter()
        if self.etape == "N1":
            self.res["N1_demarrage_ms"] = (t - T0) * 1e3
            self.res["N1_detail_ms"] = {"imports": (T_IMPORTS - T0) * 1e3,
                                        "qml_charge": (T_QML - T0) * 1e3,
                                        "premiere_image": (t - T0) * 1e3}
            self.res["N3_ouverture_initiale_ms"] = self.s.dureeChargement * 1e3
            self.etape = "pause"
            self.sonde.start()
            QTimer.singleShot(400, self._n3)
            return
        if self.attente and self.synchro:
            nom, t0 = self.attente
            self.attente = None
            self._fini(nom, (t - t0) * 1e3)
        if self.etape == "N5":
            self._n5_pas(t)

    # -- N3 : changer de run puis revenir au run de référence --
    def _n3(self):
        self.phaseN6 = "N3"
        autre = next((r for r in self.s.runs if r != self.s.nomRun), self.s.nomRun)
        self._suite = [autre, self.res.get("run", self.s.nomRun)]
        self.res["run"] = self.s.nomRun
        self._n3_suivant()

    def _n3_suivant(self):
        if not self._suite:
            QTimer.singleShot(200, self._n4_debut)
            return
        nom = self._suite.pop(0)
        # la mesure court du clic jusqu'à l'image qui suit l'installation du run
        self._n3Attente = ("N3:" + nom, time.perf_counter())
        self.s.chargerRun(nom)
        self.w.update()

    # -- N4 : balayage des frames --
    def _n4_debut(self):
        self.phaseN6 = "N4"
        self.s.frame = 0
        self._fsuiv = 0
        QTimer.singleShot(100, self._n4_suivant)

    def _n4_suivant(self):
        self._fsuiv += 1
        if self._fsuiv >= self.s.nFrames:
            QTimer.singleShot(200, self._n5_debut)
            return
        self._armer("N4")
        self.s.frame = self._fsuiv

    def _fini(self, nom, dt_ms):
        if nom.startswith("N3:"):
            self.res.setdefault("N3_rechargement", []).append(
                {"run": nom[3:], "lecture_et_tampons_ms": self.s.dureeChargement * 1e3,
                 "jusqu_image_ms": dt_ms})
            QTimer.singleShot(150, self._n3_suivant)
        elif nom == "N4":
            self.n4.append(dt_ms)
            QTimer.singleShot(30, self._n4_suivant)     # repos : aucune image en vol
        elif nom == "capture":
            QTimer.singleShot(50, self._capture)

    # -- N5 : zoom et déplacement animés pendant 3 s --
    def _n5_debut(self):
        self.phaseN6 = "N5"
        self.vue = self.racine.findChild(QObject, "vueChamps")
        self.etape = "N5"
        self.n5_t0 = time.perf_counter()
        self.n5_images = 0
        self.vue.setProperty("zoom", 1.0)
        self.w.update()

    def _n5_pas(self, t):
        ecoule = t - self.n5_t0
        if ecoule >= 3.0:
            self.res["N5_images_par_s"] = self.n5_images / ecoule
            self.etape = "capture"
            QTimer.singleShot(0, self._capture_prep)
            return
        self.n5_images += 1
        # zoom 1x -> 4x -> 1x et cercle de rayon 10 mm autour du centre
        a = 2 * np.pi * ecoule / 3.0
        self.vue.setProperty("zoom", float(1.0 + 1.5 * (1 - np.cos(a))))
        self.vue.setProperty("decalX", float(10 * np.sin(a)))
        self.vue.setProperty("decalY", float(10 * np.sin(2 * a)))

    # -- capture à la frame 20 --
    def _capture_prep(self):
        self.phaseN6 = "capture"
        self.vue.setProperty("zoom", 1.0)
        self.vue.setProperty("decalX", 0.0)
        self.vue.setProperty("decalY", 0.0)
        self._armer("capture")
        self.s.frame = min(20, self.s.nFrames - 1)
        self.w.update()

    def _capture(self):
        # QQuickWindow.grabWindow() fait planter le processus avec la boucle de rendu threadée
        # (D3D11 + View3D) : on passe par grabToImage, asynchrone, sur l'élément racine.
        self.sonde.stop()
        self._grab = self.w.contentItem().grabToImage()
        self._grab.ready.connect(self._ecrire)

    def _ecrire(self):
        img = self._grab.image()
        img.save(os.path.join(ICI, "capture_qml.png"))
        n4 = np.array(self.n4)
        self.res.update({
            "N4_frame_moy_ms": float(n4.mean()), "N4_frame_max_ms": float(n4.max()),
            "N4_frames_ms": [round(x, 2) for x in n4],
            "N6_max_boucle_ui_ms": {k: round(v, 2) for k, v in self.n6.items()},
            "capture": "capture_qml.png", "capture_taille": [img.width(), img.height()],
        })
        nom = "mesures_qml_sansvsync.json" if not self.res["vsync"] else "mesures_qml.json"
        with open(os.path.join(ICI, nom), "w", encoding="utf-8") as f:
            json.dump(self.res, f, indent=1, ensure_ascii=False)
        print(json.dumps({k: v for k, v in self.res.items() if k != "N4_frames_ms"},
                         indent=1, ensure_ascii=False))
        self.app.quit()


# --------------------------------------------------------------------------------------------
def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    mesure = "--mesure" in sys.argv
    sansvsync = "--sansvsync" in sys.argv
    run = args[0] if args else "F7_disc_gbm_P020"

    if sansvsync:
        fmt = QSurfaceFormat.defaultFormat()
        fmt.setSwapInterval(0)
        QSurfaceFormat.setDefaultFormat(fmt)
    QQuickStyle.setStyle("Basic")
    app = QGuiApplication(sys.argv)
    app.setApplicationName("rockim · Résultats (prototype QML)")

    studio = Studio(run)
    engine = QQmlApplicationEngine()
    engine.rootContext().setContextProperty("studio", studio)
    engine.load(QUrl.fromLocalFile(os.path.join(ICI, "Main.qml")))
    if not engine.rootObjects():
        sys.exit(1)
    global T_QML
    T_QML = time.perf_counter()
    fenetre = engine.rootObjects()[0]
    if mesure:
        fenetre.resize(1600, 1000)
        app._mesure = Mesure(app, fenetre, studio, fenetre, sansvsync)
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
