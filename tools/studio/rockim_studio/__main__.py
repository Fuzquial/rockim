"""Point d'entrée :  python -m rockim_studio  (depuis tools/studio/)."""
import sys

from PySide6.QtWidgets import QApplication

from .app import MainWindow


def main():
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setApplicationName("rockim-studio")
    app.setOrganizationName("rockim")
    win = MainWindow()
    win.resize(1400, 900)
    win.show()
    if len(sys.argv) > 1:
        win.ctrl.open(sys.argv[1])
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
