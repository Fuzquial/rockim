"""Tests des MODULES metier (stdlib + Qt offscreen pour le smoke)."""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
ROOT = HERE.parents[2]

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")

from rockim_studio import modules  # noqa: E402


def test_modules_sane():
    keys = [m.key for m in modules.MODULES]
    assert len(keys) == len(set(keys)), "cles de module dupliquees"
    assert keys[-1] == "expert"
    for m in modules.MODULES:
        for label, rel in m.templates:
            assert (ROOT / rel).exists(), f"{m.key}: gabarit absent {rel}"
    print("PASS test_modules_sane")


def test_groups_exist_in_registry():
    from rockim_studio.model.registry import Registry
    reg = Registry.load()
    known = {k.group for k in reg.curated.values()} if hasattr(reg, "curated") \
        else set()
    if not known:                       # API alternative : balayer les cles
        known = {getattr(k, "group", "") for k in
                 getattr(reg, "keys", {}).values()}
    known |= {"Corps et groupes", "Autres", "Inconnues"}
    for m in modules.MODULES:
        for g in m.groups:
            assert g in known, f"{m.key}: groupe inconnu du registre '{g}'"
    print("PASS test_groups_exist_in_registry")


def test_sweep_cases_lambda():
    tun = modules.by_key("tunnel")
    sw = tun.sweeps[0]
    cases = modules.sweep_cases({"E": "10e9", "insituSh": "5e6",
                                 "insituSv": "5e6"}, sw, [0.5, 1.5])
    assert [t for t, _ in cases] == ["0p5", "1p5"]
    p05 = dict(cases[0][1])
    assert p05["insituSh"] == "5e6" and abs(float(p05["insituSv"]) - 1e7) < 1
    p15 = dict(cases[1][1])
    assert abs(float(p15["insituSv"]) - 5e6 / 1.5) < 1
    assert p05["E"] == "10e9", "les cles du cas de base doivent survivre"
    print("PASS test_sweep_cases_lambda")


def test_gui_module_smoke():
    try:
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print("SKIP test_gui_module_smoke (PySide6 absent)")
        return
    from rockim_studio.app import MainWindow
    app = QApplication.instance() or QApplication([])
    w = MainWindow()
    tun = modules.by_key("tunnel")
    w.set_module(tun)
    assert w.module.key == "tunnel"
    assert w.tree.allowed == set(tun.groups)
    assert w._sweep_menu.isEnabled()
    w.set_module(modules.by_key("expert"))
    assert w.tree.allowed is None
    assert not w._sweep_menu.isEnabled()
    w.close()
    print("PASS test_gui_module_smoke")


if __name__ == "__main__":
    test_modules_sane()
    test_groups_exist_in_registry()
    test_sweep_cases_lambda()
    test_gui_module_smoke()
