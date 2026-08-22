"""app.py — la fenêtre principale de rockim-studio (spec 006, WP0.1/0.4/0.5).

Disposition M0 (la scène 3D PyVista arrive en M1 au centre) :
  gauche  = arbre du modèle ; droite = propriétés du groupe sélectionné ;
  centre  = courbe live du run ; bas = console (journal + validation).
"""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QAction, QKeySequence
from PySide6.QtWidgets import (QDockWidget, QFileDialog, QInputDialog,
                               QMainWindow, QMessageBox, QSpinBox, QToolBar,
                               QWidget)

from . import modules
from .controller import Controller
from .model.cfg_io import CfgFile
from .run.monitor import HistoryMonitor
from .run.runner import Runner
from .views.console import Console
from .views.plots import LivePlot
from .views.props import PropertyPanel
from .views.tree import ModelTree


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("rockim-studio")
        self.settings = QSettings("rockim", "studio")

        self.ctrl = Controller(parent=self)
        self.runner = Runner(parent=self)
        self.monitor = HistoryMonitor(parent=self)

        # centre : scène 3D + courbes
        from PySide6.QtWidgets import QTabWidget

        from .views.geometry import GeometryPanel
        from .views.scene import SceneView
        self.geometry = GeometryPanel(self.ctrl)
        self.scene = SceneView()
        self.plot = LivePlot()
        self.center = QTabWidget()
        self.center.addTab(self.geometry, "Géométrie")
        self.center.addTab(self.scene, "Résultats 3D")
        self.center.addTab(self.plot, "Courbes")
        self.setCentralWidget(self.center)
        self.center.currentChanged.connect(
            lambda i: self.console.append_log(
                "onglet actif : %s" % self.center.tabText(i))
            if hasattr(self, "console") else None)
        self.geometry.mesh_ready.connect(self._mesh_ready)
        self.scene.point_picked.connect(self._on_point_picked)

        # docks
        self.tree = ModelTree(self.ctrl)
        d_tree = self._dock("Modèle", self.tree, Qt.LeftDockWidgetArea)
        self.props = PropertyPanel(self.ctrl)
        d_props = self._dock("Propriétés", self.props, Qt.RightDockWidgetArea)
        self.console = Console(self.ctrl)
        d_console = self._dock("Console", self.console,
                               Qt.BottomDockWidgetArea)
        self._docks = (d_tree, d_props, d_console)

        self.tree.group_selected.connect(self.props.show_group)
        self.ctrl.model_reset.connect(self._refresh_title)
        self.ctrl.key_changed.connect(lambda _k: self._refresh_title())

        self.runner.output.connect(self.console.append_log)
        self.runner.started.connect(self._run_started)
        self.runner.finished.connect(self._run_finished)
        self.monitor.header_ready.connect(self.plot.set_header)
        self.monitor.rows_added.connect(self.plot.add_rows)

        self._build_actions()
        self._restore_state()
        # ---- MODULE metier (Fernando 2026-08-22) : choisi au demarrage,
        # memorise ; l'interface ne montre que ce qui appartient au module.
        mk = self.settings.value("module", "")
        if not mk:
            mk = self._choose_module(initial=True)
        self.set_module(modules.by_key(mk), startup=True)
        self.ctrl.new()

    # --- actions et barre d'outils ---------------------------------------
    def _build_actions(self):
        bar = QToolBar("Principal")
        bar.setObjectName("mainToolbar")
        self.addToolBar(bar)
        self._toolbar = bar
        menu_f = self.menuBar().addMenu("&Fichier")
        menu_e = self.menuBar().addMenu("&Édition")
        menu_v = self.menuBar().addMenu("&Affichage")
        menu_r = self.menuBar().addMenu("&Calcul")
        # panneaux fermables -> TOUJOURS re-ouvrables (bug signale par
        # Fernando : un dock ferme etait perdu, faute de menu Affichage)
        for d in self._docks:
            menu_v.addAction(d.toggleViewAction())
        menu_v.addSeparator()
        ra = QAction("&Réinitialiser la disposition", self)
        ra.triggered.connect(self._reset_layout)
        menu_v.addAction(ra)

        def act(text, slot, seq=None, menus=(), toolbar=False):
            a = QAction(text, self)
            if seq:
                a.setShortcut(QKeySequence(seq))
            a.triggered.connect(slot)
            for m in menus:
                m.addAction(a)
            if toolbar:
                bar.addAction(a)
            return a

        act("Changer de &module…", self._change_module, None, (menu_f,))
        menu_f.addSeparator()
        act("&Nouveau (fdem)", self._new, "Ctrl+N", (menu_f,))
        self._build_templates_menu(menu_f)
        act("&Ouvrir un cfg…", self._open, "Ctrl+O", (menu_f,), True)
        act("&Enregistrer", self._save, "Ctrl+S", (menu_f,), True)
        act("Enregistrer &sous…", self._save_as, "Ctrl+Shift+S", (menu_f,))
        menu_f.addSeparator()
        act("&Quitter", self.close, "Ctrl+Q", (menu_f,))

        for i, seq in ((0, "Ctrl+1"), (1, "Ctrl+2"), (2, "Ctrl+3")):
            a = QAction(self)
            a.setShortcut(QKeySequence(seq))
            a.triggered.connect(
                lambda _c=False, k=i: self.center.setCurrentIndex(k))
            self.addAction(a)
        act("&Annuler", self.ctrl.undo, "Ctrl+Z", (menu_e,))
        act("&Rétablir", self.ctrl.redo, "Ctrl+Shift+Z", (menu_e,))

        bar.addSeparator()
        act("&Lancer", self._launch, "F5", (menu_r,), True)
        act("&Arrêter", self._stop, "Shift+F5", (menu_r,), True)
        act("&Exécutable rockim…", self._pick_exe, None, (menu_r,))
        act("Ouvrir un &dossier de résultats…", self._open_results,
            "Ctrl+R", (menu_r,), True)
        act("&Comparer avec un run…", self._compare_run, None, (menu_r,))
        self._sweep_menu = menu_r.addMenu("&Balayages du module")
        menu_r.addSeparator()

        self.threads = QSpinBox()
        self.threads.setRange(0, 128)
        self.threads.setToolTip("OMP_NUM_THREADS (0 = environnement)")
        self.threads.setValue(int(self.settings.value("threads", 0)))
        bar.addWidget(self.threads)

        self.statusBar().showMessage("prêt")
        self.ctrl.validation_changed.connect(self._status_validation)

    # Gabarits FDEM : les configs de référence du dépôt, par essai. Chemins
    # relatifs à la racine du repo (détectée depuis ce fichier).
    _TEMPLATES = [
        # — les trois filières de production —
        ("Tunnel EDZ pressurisé (rapide)", "configs/tunnel_bore_fast.cfg"),
        ("Tunnel EDZ pressurisé (production)", "configs/tunnel_bore.cfg"),
        ("Tunnel EDZ Weibull", "configs/tunnel_bore_weib.cfg"),
        None,
        ("Impact 3D smoke (bench1 réduit)", "configs/smoke_impact.cfg"),
        ("Impact banc St Anne s1,5 (spec 005)",
         "bench_impact/configs/impact_stanne_s15.cfg"),
        ("Impact St Anne FIDÈLE (tout comme eux sauf adaptatif)",
         "bench_impact/configs/impact_stanne_fidele_s15.cfg"),
        None,
        ("Hydro-frac Abu-Aisha ISO (grossier, hydro=on)",
         "bench_abuaisha/configs/hf_iso_hydro_c.cfg"),
        ("Hydro-frac Abu-Aisha ISO (production)",
         "bench_abuaisha/configs/hf_iso_hydro.cfg"),
        ("Hydro-frac Abu-Aisha ANISO (production)",
         "bench_abuaisha/configs/hf_aniso_hydro.cfg"),
        None,
        # — essais de laboratoire et divers —
        ("Percussion 2D (insert disque)", "configs/fdem_percussion.cfg"),
        ("Percussion 3D (insert sphère)", "configs/fdem3d_percussion.cfg"),
        ("Percussion 3D GBM Voronoï",
         "configs/fdem3d_voronoi_percussion.cfg"),
        ("UCS Bohus (platines)", "configs/cal_ucs_bohus.cfg"),
        ("Brésilien Bohus (disque)", "configs/cal_bts_bohus.cfg"),
        ("Triaxial Bohus GBM", "configs/triax_bohus_gbm.cfg"),
        ("Coupe 2D (couteau PDC)", "configs/fdem_shear.cfg"),
        ("Vérification traction FDEM", "configs/verify_fdem_tension.cfg"),
    ]

    def _build_templates_menu(self, menu_parent):
        self._tpl_menu = menu_parent.addMenu("Nouveau depuis un &modèle FDEM")
        self._fill_templates_menu()

    def _fill_templates_menu(self):
        root = Path(__file__).resolve().parents[3]
        self._tpl_menu.clear()
        mod = getattr(self, "module", None)
        entries = (list(mod.templates) if mod and mod.templates
                   else self._TEMPLATES)
        for entry in entries:
            if entry is None:
                self._tpl_menu.addSeparator()
                continue
            label, rel = entry
            path = root / rel
            a = QAction(label, self)
            a.setEnabled(path.exists())
            a.triggered.connect(
                lambda _c=False, p=path, lbl=label: self._from_template(
                    p, lbl))
            self._tpl_menu.addAction(a)

    # --- modules metier ----------------------------------------------------
    def _choose_module(self, initial=False):
        from PySide6.QtWidgets import QInputDialog
        items = ["%s — %s" % (m.label, m.doc) for m in modules.MODULES]
        cur = 0
        saved = self.settings.value("module", "")
        for i, m in enumerate(modules.MODULES):
            if m.key == saved:
                cur = i
        text, ok = QInputDialog.getItem(
            self, "Choisir un module",
            "Le studio ne montrera que les gabarits, groupes de clés et\n"
            "actions de ce module (le mode Expert montre tout) :",
            items, cur, False)
        if not ok and initial:
            return "expert"
        if not ok:
            return self.module.key
        return modules.MODULES[items.index(text)].key

    def _change_module(self):
        self.set_module(modules.by_key(self._choose_module()))

    def set_module(self, mod, startup=False):
        self.module = mod
        self.settings.setValue("module", mod.key)
        self.tree.set_allowed(set(mod.groups) if mod.groups else None)
        if hasattr(self, "_tpl_menu"):
            self._fill_templates_menu()
        if hasattr(self, "_sweep_menu"):
            self._sweep_menu.clear()
            self._sweep_menu.setEnabled(bool(mod.sweeps))
            for sw in mod.sweeps:
                a = QAction(sw.name, self)
                a.triggered.connect(lambda _c=False, s=sw: self._do_sweep(s))
                self._sweep_menu.addAction(a)
        self.setWindowTitle("rockim-studio — module %s" % mod.label)
        if not startup:
            self.statusBar().showMessage("module : %s" % mod.label)

    def _do_sweep(self, sweep):
        """Genere un deck par valeur (copie du cas courant + cles du
        balayage) et, sur demande, met les runs en FILE (le runner est
        deja sequentiel : un seul job a la fois)."""
        from PySide6.QtWidgets import QCheckBox, QInputDialog
        vals, ok = QInputDialog.getText(
            self, sweep.name, sweep.doc + "\n\nValeurs :",
            text=sweep.values)
        if not ok or not vals.strip():
            return
        try:
            values = [float(v) for v in vals.split()]
        except ValueError:
            QMessageBox.warning(self, "Balayage", "Valeurs illisibles.")
            return
        base, ok = QInputDialog.getText(
            self, sweep.name, "Prefixe des decks et des out_dir :",
            text=self.settings.value("sweepBase", "sweep"))
        if not ok or not base:
            return
        self.settings.setValue("sweepBase", base)
        lancer = QMessageBox.question(
            self, sweep.name,
            "Mettre les %d runs en file maintenant ?\n(Non = generer "
            "seulement les decks)" % len(values),
            QMessageBox.Yes | QMessageBox.No) == QMessageBox.Yes
        model = self.ctrl.model
        root = Path(model.source_dir or Path.cwd())
        made = []
        for tag, pairs in modules.sweep_cases(dict(model.cfg.pairs),
                                              sweep, values):
            copy = CfgFile(pairs=pairs, comments=list(model.cfg.comments))
            mesh = model.mesh_file_path()
            if mesh is not None:
                copy.pairs["meshFile"] = str(mesh)
            cfg_path = root / ("%s_%s.cfg" % (base, tag))
            copy.write(cfg_path, header="balayage %s — rockim-studio"
                       % sweep.name)
            made.append((cfg_path, root / ("out_%s_%s" % (base, tag))))
        self.console.append_log("balayage : %d decks écrits (%s)"
                                % (len(made), root))
        if lancer:
            exe = self.settings.value("exe", "")
            if not exe or not Path(exe).exists():
                self._pick_exe()
                exe = self.settings.value("exe", "")
                if not exe:
                    return
            self.runner.exe = exe
            self.runner.threads = self.threads.value()
            for cfg_path, out_dir in made:
                self.runner.launch(cfg_path, out_dir)
            self.console.append_log(
                "balayage : %d runs en file (sequentiels)" % len(made))

    def _from_template(self, path, label):
        if self._confirm_discard():
            self.ctrl.open_template(path, label)

    # --- fichier ----------------------------------------------------------
    def _confirm_discard(self) -> bool:
        if not self.ctrl.model.dirty:
            return True
        r = QMessageBox.question(self, "Modifications non enregistrées",
                                 "Abandonner les modifications en cours ?")
        return r == QMessageBox.Yes

    def _new(self):
        if self._confirm_discard():
            self.ctrl.new()

    def _open(self):
        if not self._confirm_discard():
            return
        start = self.settings.value("lastDir", str(Path.cwd()))
        path, _f = QFileDialog.getOpenFileName(
            self, "Ouvrir une configuration", start, "Config rockim (*.cfg)")
        if path:
            self.settings.setValue("lastDir", str(Path(path).parent))
            self.ctrl.open(path)

    def _save(self):
        if self.ctrl.model.path is None:
            self._save_as()
        else:
            self.ctrl.save()
            self._refresh_title()

    def _save_as(self):
        start = self.settings.value("lastDir", str(Path.cwd()))
        path, _f = QFileDialog.getSaveFileName(
            self, "Enregistrer la configuration", start,
            "Config rockim (*.cfg)")
        if path:
            self.settings.setValue("lastDir", str(Path(path).parent))
            self.ctrl.save(path)
            self._refresh_title()

    # --- calcul -----------------------------------------------------------
    def _pick_exe(self):
        path, _f = QFileDialog.getOpenFileName(
            self, "Exécutable rockim",
            self.settings.value("exe", str(Path.cwd())),
            "Exécutables rockim (rockim* rockim*.exe);;Tous (*)")
        # rockim*N : la convention du depot est UN exe dedie par
        # chantier (rockim_i4.exe, rockim_e3.exe...) — le filtre
        # initial ne montrait que "rockim.exe" et masquait tout le parc
        if path:
            self.settings.setValue("exe", path)

    def _launch(self):
        issues = self.ctrl.model.validate()
        errors = [i for i in issues if i[0] == "erreur"]
        if errors:
            QMessageBox.warning(
                self, "Validation",
                "Erreurs bloquantes :\n" + "\n".join(
                    f"• {k} : {m}" for _l, k, m in errors[:12]))
            return
        exe = self.settings.value("exe", "")
        if not exe or not Path(exe).exists():
            self._pick_exe()
            exe = self.settings.value("exe", "")
            if not exe:
                return
        out, ok = QInputDialog.getText(
            self, "Dossier de sortie", "out_dir :",
            text=self.settings.value("lastOut", "out_studio"))
        if not ok or not out:
            return
        self.settings.setValue("lastOut", out)
        self.settings.setValue("threads", self.threads.value())
        self.launch_case(out, exe=exe, threads=self.threads.value())

    def launch_case(self, out_dir: str | Path, exe: str = "",
                    threads: int = 0):
        """Lancement NON interactif (utilisé par l'UI et les tests) : copie
        la config dans out_dir avec le meshFile ABSOLUTISÉ (les decks du
        dépôt utilisent des chemins relatifs à la racine, qui casseraient
        depuis la copie), puis lance."""
        model = self.ctrl.model
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        cfg_path = out_dir / "studio.cfg"
        copy = CfgFile(pairs=dict(model.cfg.pairs),
                       comments=list(model.cfg.comments))
        mesh = model.mesh_file_path()
        if mesh is not None:
            copy.pairs["meshFile"] = str(mesh)
        copy.write(cfg_path, header="écrit par rockim-studio")
        self.runner.exe = exe or self.settings.value("exe", "")
        self.runner.threads = threads
        self.runner.launch(cfg_path, out_dir)

    def _stop(self):
        self.runner.stop()

    def _run_started(self, out_dir: str):
        self.console.append_log(f"=== run lancé -> {out_dir} ===")
        self.plot.reset()
        self.monitor.watch(out_dir)
        self.statusBar().showMessage(f"run en cours : {out_dir}")

    def _run_finished(self, code: int, out_dir: str):
        self.monitor.stop()
        verdict = "OK" if code == 0 else f"code retour {code}"
        self.console.append_log(f"=== run terminé ({verdict}) ===")
        summary = Path(out_dir) / "summary.txt"
        if summary.exists():
            self.console.append_log(summary.read_text(
                encoding="utf-8", errors="replace"))
        self.statusBar().showMessage(f"terminé : {out_dir} ({verdict})")
        if code == 0:
            self.scene.load(out_dir)
            self.center.setCurrentWidget(self.scene)

    def _open_results(self):
        start = self.settings.value("lastOut", str(Path.cwd()))
        path = QFileDialog.getExistingDirectory(
            self, "Dossier de résultats (out_*)", start)
        if path:
            self.load_results(path)

    def load_results(self, path: str):
        self.scene.load(path)
        self.plot.load_csv(Path(path) / "history.csv")
        if self.scene.series is not None and len(self.scene.series):
            self.plot.attach_series(self.scene.series)
        self.center.setCurrentWidget(self.scene)
        self.console.append_log(f"résultats chargés : {path}")

    def _on_point_picked(self, x: float, y: float, z: float):
        self.plot.set_probe_point(x, y, z)
        self.center.setCurrentWidget(self.plot)

    def _mesh_ready(self, result: dict):
        self.scene.show_mesh(result["vtu"])
        self.center.setCurrentWidget(self.scene)

    def _compare_run(self):
        start = self.settings.value("lastOut", str(Path.cwd()))
        path = QFileDialog.getExistingDirectory(
            self, "Run de référence à superposer (out_*)", start)
        if not path:
            self.plot.set_reference(None)
            return
        header, data = self.plot.read_csv(Path(path) / "history.csv")
        if header is None:
            self.console.append_log(f"pas de history.csv dans {path}")
            return
        self.plot.set_reference(Path(path).name, header, data)
        self.center.setCurrentWidget(self.plot)

    # --- divers -----------------------------------------------------------
    def _dock(self, title: str, widget: QWidget, area) -> QDockWidget:
        d = QDockWidget(title, self)
        d.setObjectName(title)
        d.setWidget(widget)
        self.addDockWidget(area, d)
        return d

    def _refresh_title(self):
        name = self.ctrl.model.path.name if self.ctrl.model.path \
            else "sans titre"
        star = " *" if self.ctrl.model.dirty else ""
        self.setWindowTitle(f"rockim-studio — {name}{star} "
                            f"[{self.ctrl.model.mode}]")

    def _status_validation(self, issues: list):
        n_err = sum(1 for i in issues if i[0] == "erreur")
        if n_err:
            self.statusBar().showMessage(
                f"{n_err} erreur(s) de validation — voir la console")
        elif "erreur" in self.statusBar().currentMessage():
            self.statusBar().showMessage("validation OK")

    def _restore_state(self):
        geo = self.settings.value("geometry")
        if geo:
            self.restoreGeometry(geo)
        state = self.settings.value("windowState")
        if state:
            self.restoreState(state)
        else:
            self._first_layout = True

    def showEvent(self, event):
        super().showEvent(event)
        # resizeDocks n'a d'effet qu'une fois la fenêtre réellement layoutée
        if getattr(self, "_first_layout", False):
            self._first_layout = False
            d_tree, d_props, d_console = self._docks
            self.resizeDocks([d_tree, d_props], [300, 420], Qt.Horizontal)
            self.resizeDocks([d_console], [220], Qt.Vertical)

    def closeEvent(self, event):
        if self.runner.busy:
            r = QMessageBox.question(self, "Run en cours",
                                     "Un run tourne encore. L'arrêter et "
                                     "quitter ?")
            if r != QMessageBox.Yes:
                event.ignore()
                return
            self.runner.stop()
        if not self._confirm_discard():
            event.ignore()
            return
        self.settings.setValue("geometry", self.saveGeometry())
        self.settings.setValue("windowState", self.saveState())
        super().closeEvent(event)

    def _reset_layout(self):
        """Re-ouvre tous les panneaux et rend les largeurs par defaut."""
        for d in self._docks:
            d.show()
        if hasattr(self, "_toolbar"):
            self._toolbar.show()
        d_tree, d_props, d_console = self._docks
        self.resizeDocks([d_tree, d_props], [300, 420], Qt.Horizontal)
        self.resizeDocks([d_console], [170], Qt.Vertical)
        self.statusBar().showMessage("disposition réinitialisée")
