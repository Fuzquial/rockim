# En-têtes de rockim

`include/rockim/` contient les déclarations des solveurs et les noyaux partagés écrits en fonctions
pures (loi cohésive `JointTsl.hpp`, DIF `YangDif.hpp`, adoucissement `YanSoftening.hpp`, contact par
potentiel `PotentialContact.hpp`, contact outil `ToolSignorini.hpp`, cutter `ToolPdc3d.hpp`).

La carte du code, le rôle de chaque en-tête et le déroulé d'un pas de temps sont dans
[`src/README.md`](../src/README.md).

`KeysByMode.hpp` est généré par `tools/gen_keys_by_mode.py` : ne pas l'éditer à la main.
