"""Tests WP0.2/WP0.3 (spec 006) : registre des clés et round-trip cfg.

Lancement :  python3 -m pytest tools/studio/tests -q
(ou python3 tools/studio/tests/test_registry_roundtrip.py sans pytest).
Aucune dépendance Qt/VTK — stdlib seulement.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools/studio"))

from rockim_studio.model.cfg_io import CfgFile            # noqa: E402
from rockim_studio.model.registry import Registry          # noqa: E402

CONFIGS = sorted((ROOT / "configs").glob("*.cfg"))


def test_extraction_up_to_date():
    """Garde anti-dérive R3 : le JSON commité colle à src/*.cpp."""
    r = subprocess.run([sys.executable,
                        str(ROOT / "tools/studio/dev/extract_keys.py"),
                        "--check"], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_registry_loads_and_is_sane():
    reg = Registry.load()
    assert len(reg.keys) > 200, "extraction suspecte (< 200 clés)"
    assert reg.phantom_curated() == [], (
        "entrées CURATED sans clé dans le C++ : "
        f"{reg.phantom_curated()}")
    # quelques clés pivots, avec leur type attendu
    assert reg.keys["jointXi"].type == "float"
    assert reg.keys["nx"].type == "int"
    assert reg.keys["historyFlush"].type == "bool"
    assert reg.keys["mode"].type == "str"
    # portée : jointXi est lue par les deux solveurs FDEM
    assert {"fdem", "fdem3d"} <= set(reg.keys["jointXi"].scope)


def test_roundtrip_all_configs():
    """cfg -> modèle -> cfg -> modèle : mêmes paires, pour les 104 configs."""
    assert CONFIGS, "configs/ introuvable"
    for path in CONFIGS:
        first = CfgFile.parse(path)
        second = CfgFile().parse_text(first.dumps())
        assert second.pairs == first.pairs, f"round-trip cassé : {path.name}"


def test_all_config_keys_known_to_registry():
    """Toute clé utilisée par une config du dépôt existe dans le registre
    (ou appartient à une famille dynamique). Détecte les clés mortes des
    configs comme les trous de l'extracteur."""
    reg = Registry.load()
    unknown = {}
    for path in CONFIGS:
        bad = reg.unknown_keys(CfgFile.parse(path).pairs)
        if bad:
            unknown[path.name] = bad
    assert not unknown, f"clés inconnues du registre : {unknown}"


def _main():
    fails = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as e:
                print(f"FAIL {name}: {e}")
                fails += 1
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    _main()
