# rockim-studio — WP0.2/WP0.3 (spec 006)

Premier livrable du chantier GUI (specs/006-interface-graphique/spec.md) :

* `dev/extract_keys.py` — extrait le registre des clés depuis les sites
  d'appel `getd/geti/getb/gets/reqd/reqs` de `src/*.cpp` + `include/rockim/*.hpp`
  → `rockim_studio/model/keys_extracted.json` (256 clés). `--check` = garde
  anti-dérive (échoue si le JSON commité ne colle plus au C++).
* `rockim_studio/model/registry.py` — le registre (types, défauts par mode,
  portée, familles dynamiques `phase.*`/`group*`/`gauge.*`) + métadonnées
  d'ergonomie CURATED (groupes UI, docs, bornes, énumérations).
* `rockim_studio/model/cfg_io.py` — lecture/écriture `.cfg`, sémantique
  identique à `Config::load`, round-trip garanti.
* `tests/test_registry_roundtrip.py` — 4 tests, stdlib seulement :
  `python3 tools/studio/tests/test_registry_roundtrip.py`
  (round-trip vérifié sur les 104 configs du dépôt ; toute clé de config
  inconnue du registre fait échouer la suite).

Après toute modification des clés dans le C++ : relancer
`python3 tools/studio/dev/extract_keys.py` et committer le JSON.
