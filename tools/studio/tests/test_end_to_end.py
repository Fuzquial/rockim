"""Test bout-en-bout M0 (spec 006) : Runner + HistoryMonitor sur un VRAI run.

Nécessite PySide6 (offscreen) et l'exécutable build/rockim ; sinon le test
se saute proprement. Cas court : verify_fdem_tension.cfg (la référence de
la suite du solveur) lancé PAR LE RUNNER DU STUDIO, courbe suivie PAR LE
MONITOR, verdict lu dans summary.txt.
"""
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tools/studio"))

EXE = ROOT / "build/rockim"


def test_runner_monitor_end_to_end():
    try:
        from PySide6.QtCore import QCoreApplication, QEventLoop, QTimer
    except ImportError:
        print("  (PySide6 absent : bout-en-bout sauté)")
        return
    if not EXE.exists():
        print("  (build/rockim absent : bout-en-bout sauté)")
        return
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from rockim_studio.model.rockim_model import RockimModel
    from rockim_studio.run.monitor import HistoryMonitor
    from rockim_studio.run.runner import Runner

    app = QCoreApplication.instance() or QCoreApplication([])
    with tempfile.TemporaryDirectory() as td:
        out_dir = Path(td) / "out_e2e"
        # le studio écrit toujours SA copie du cfg dans le dossier de sortie
        model = RockimModel()
        model.open(ROOT / "configs/verify_fdem_tension.cfg")
        assert model.validate() == [], model.validate()
        out_dir.mkdir(parents=True)
        cfg_path = out_dir / "studio.cfg"
        model.cfg.write(cfg_path, header="écrit par rockim-studio (test e2e)")

        runner = Runner()
        runner.exe = str(EXE)
        runner.threads = 2
        monitor = HistoryMonitor(interval_ms=200)

        state = {"rows": 0, "code": None, "log": []}
        monitor.rows_added.connect(
            lambda rows: state.__setitem__("rows", state["rows"] + len(rows)))
        runner.output.connect(lambda t: state["log"].append(t))
        runner.started.connect(lambda out: monitor.watch(out))

        loop = QEventLoop()

        def done(code, _out):
            monitor.stop()
            state["code"] = code
            loop.quit()

        runner.finished.connect(done)
        QTimer.singleShot(180000, loop.quit)   # garde-fou 3 min
        runner.launch(cfg_path, out_dir)
        loop.exec()

        assert state["code"] == 0, \
            f"run non abouti (code {state['code']}) : {state['log'][-8:]}"
        assert state["rows"] > 10, \
            f"monitor n'a vu que {state['rows']} lignes de history.csv"
        # le verdict PASS/FAIL de la vérification arrive sur stdout — c'est
        # exactement ce que la console du studio affiche
        log = "\n".join(state["log"])
        assert "[PASS]" in log and "[FAIL]" not in log, \
            f"verdict solveur manquant ou FAIL :\n{log[-800:]}"
        assert (out_dir / "history.csv").exists()
        print(f"  run OK ([PASS] au journal), "
              f"{state['rows']} lignes suivies en live")


if __name__ == "__main__":
    try:
        test_runner_monitor_end_to_end()
        print("PASS test_runner_monitor_end_to_end")
    except AssertionError as e:
        print(f"FAIL test_runner_monitor_end_to_end: {e}")
        sys.exit(1)
