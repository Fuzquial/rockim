"""Tests WP0.1/0.4 : RockimModel (sans Qt) + fumée GUI si PySide6 présent."""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools/studio"))

from rockim_studio.model.rockim_model import RockimModel   # noqa: E402


def test_model_open_edit_save(tmp_path=None):
    m = RockimModel()
    m.open(ROOT / "configs/fdem_percussion.cfg")
    assert m.mode == "fdem"
    assert not m.dirty
    # défaut du mode courant quand la clé n'est pas posée
    if not m.is_explicit("jointXi"):
        assert m.value("jointXi") == "0.05"
    m.set_value("jointXi", "0.01")
    assert m.dirty and m.value("jointXi") == "0.01"
    out = Path(tmp_path or "/tmp") / "studio_model_test.cfg"
    m.save(out)
    m2 = RockimModel()
    m2.open(out)
    assert m2.cfg.pairs == m.cfg.pairs


def test_new_case_is_fdem_and_templates_are_copies():
    m = RockimModel()
    m.new()
    assert m.mode == "fdem" and not m.dirty
    m.open_template(ROOT / "configs/cal_ucs_bohus.cfg")
    assert m.path is None and m.dirty          # jamais d'écrasement de la réf
    assert m.mode in ("fdem", "fdem3d")
    try:
        m.save()
        raise RuntimeError("save() sans chemin aurait dû échouer")
    except ValueError:
        pass


def test_groups_follow_mode():
    m = RockimModel()
    m.set_value("mode", "fdem3d")
    g3 = m.groups()
    joints = {k.name for k in g3.get("Joints", [])}
    assert "jointXi" in joints
    # une clé explicite hors portée reste visible
    m.set_value("mode", "fem")
    m.set_value("jointXi", "0.02")
    gf = m.groups()
    assert any(k.name == "jointXi" for ks in gf.values() for k in ks)


def test_validation_catches_house_traps():
    m = RockimModel()
    m.set_value("jointXi", "0,05")          # virgule FR
    m.set_value("nu", "0.7")                # hors bornes
    m.set_value("mode", "fdme")             # énumération
    m.set_value("cleImaginaire", "1")       # inconnue
    levels = {(lvl, key) for lvl, key, _m in m.validate()}
    assert ("erreur", "jointXi") in levels
    assert ("erreur", "nu") in levels
    assert ("erreur", "mode") in levels
    assert ("alerte", "cleImaginaire") in levels
    # familles dynamiques : pas des inconnues
    m.set_value("groupVel.insert", "0 0 -8")
    assert not any(k == "groupVel.insert" for _l, k, _m in m.validate())


def test_gui_smoke_if_pyside_available():
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print("  (PySide6 absent : fumée GUI sautée)")
        return
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from rockim_studio.app import MainWindow
    app = QApplication.instance() or QApplication([])
    win = MainWindow()
    win.ctrl.open(ROOT / "configs/fdem_percussion.cfg")
    win.ctrl.set_key("jointXi", "0.01")
    assert win.ctrl.model.dirty
    win.ctrl.undo()
    assert not win.ctrl.model.is_explicit("jointXi") \
        or win.ctrl.model.value("jointXi") != "0.01"
    win.props.show_group("Joints")
    win.tree.rebuild()
    app.processEvents()
    win.ctrl.model.dirty = False    # pas de dialogue à la fermeture
    win.close()


def _main():
    import tempfile
    fails = 0
    with tempfile.TemporaryDirectory() as td:
        for name, fn in sorted(globals().items()):
            if not name.startswith("test_"):
                continue
            try:
                fn(td) if "tmp_path" in fn.__code__.co_varnames else fn()
                print(f"PASS {name}")
            except AssertionError as e:
                print(f"FAIL {name}: {e}")
                fails += 1
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    _main()
