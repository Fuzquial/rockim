# Tout récupérer sur un poste neuf — 13/09/2026

Deux dépôts portent le travail. Selon ce que vous voulez faire, vous n'avez pas besoin des deux.

| Dépôt | Ce qu'il porte | À télécharger | Accès |
|---|---|---|---|
| `Fuzquial/rockim`, branche `g1` | le code, les decks, les outils, les figures, les modes d'emploi | **112 Mo** | **public** |
| `Fuzquial/phd_geothermie`, branche `master` | toute la thèse : Abaqus, FDEM, vault Obsidian, Overleaf, et le miroir `FDEM/rockim_g1` | **1,0 Go** (≈ 3,2 Go sur disque après extraction) | privé |

⚠️ `Fuzquial/rockim` est **public**, pas privé : le solveur, les decks et les figures sont visibles de
tous. Si ce n'était pas voulu, c'est à changer dans les réglages du dépôt sur GitHub (Settings →
General → Danger Zone → Change visibility).

## 1. Les outils, une fois par poste

Sur un PC Windows, dans PowerShell :

```
winget install --id Git.Git -e
winget install --id GitHub.cli -e
gh auth login
```

`gh auth login` n'est nécessaire que pour `phd_geothermie`, qui est privé. Choisir GitHub.com, HTTPS,
puis l'authentification par navigateur. Sur macOS, `brew install git gh` fait la même chose.

## 2. Pour écrire le rapport : le dépôt rockim suffit

C'est le cas le plus fréquent, et le plus léger. 112 Mo, aucune authentification.

```
git clone --branch g1 https://github.com/Fuzquial/rockim.git
cd rockim
```

Vous y trouvez :

- **220 figures** dont **90 PDF vectoriels en Computer Modern**, directement inclusibles en LaTeX,
  dans `results/fig/` et `output/pdf/` — les planches de St Anne à 192 µs sont dans
  `results/fig/stanne192/` et `output/pdf/stanne_*_192us/` ;
- **22 documents** dans `docs/`, dont `ECARTS_guo2014_rockim_2026-09-13.md` qui porte toutes les
  mesures chiffrées, et `REPRODUIRE_stanne_radiales_2026-09-14.md` qui dit comment rejouer le calcul ;
- les decks dans `configs/`, les scripts de figures dans `tools/`, les bancs dans `tests_f2/` ;
- `rockim_g1y19.exe`, le binaire Windows qui a produit les fissures radiales.

## 3. Pour tout avoir : l'archive de thèse

```
gh repo clone Fuzquial/phd_geothermie
```

1,0 Go à télécharger. Compter dix à trente minutes selon la ligne. L'archive contient en plus du
dépôt rockim : les campagnes Abaqus (`CONTINUUM/`, 3241 images dont 209 animations), le vault Obsidian
à la racine (ouvrir le dossier comme coffre dans Obsidian), le manuscrit (`OVERLEAF/`), et le miroir
`FDEM/rockim_g1/` sans les sorties lourdes.

**Version allégée**, si vous ne voulez qu'une partie :

```
git clone --filter=blob:none --sparse https://github.com/Fuzquial/phd_geothermie.git
cd phd_geothermie
git sparse-checkout set CERVEAU.md MOC ETUDES JOURNAL FDEM/rockim_g1/docs FDEM/rockim_g1/results/fig
```

Git ne télécharge alors que les fichiers des dossiers listés, et va chercher les autres à la demande.
Ajouter un dossier plus tard : relancer `git sparse-checkout set` avec la liste complétée.

## 4. Pour relancer des calculs

Il faut recompiler : le `.exe` versionné est un binaire Windows, et les sources font foi.

**Sur Windows**, avec Visual Studio 2022 Community installé (CMake et Ninja viennent avec, rien à
installer de plus) :

```
powershell -ExecutionPolicy Bypass -File tools\build.ps1 -Jobs 8
```

Le résultat est `build\rockim.exe`. Contrôle de santé : `python tools\bitid.py --exe build\rockim.exe --list`.

**Sur macOS**, la marche à suivre est dans `docs/INSTALLER_macOS_2026-09-14.md` : clang, `brew install
cmake libomp`, puis CMake. Attention, Apple clang ne fournit pas OpenMP, sans `libomp` le code tourne
en série et quatre à dix fois plus lentement.

**Le maillage n'est pas dans le dépôt** (5,5 Mo) mais se régénère à l'identique, le mailleur fixe sa
graine :

```
python tools/make_impact_mesh.py meshes/impact_yang_train1_rock137_hxt.msh 1.0 2e-5 1.37 gap=2e-5 quality=hxt train=fixed
```

Python demande `numpy`, `matplotlib`, `scipy` et `gmsh` (`python -m pip install ...`).

## 5. Ce qui n'est dans aucun des deux dépôts

Les sorties de calcul : les dossiers `out_*/` et leurs trames VTU, six gigaoctets pour la seule
branche `g1`. Elles ne sont pas perdues pour autant, elles se rejouent depuis les decks, et toutes les
mesures qu'on en a tirées sont consignées dans `docs/`. Les maillages `.msh` et les binaires
intermédiaires sont exclus pour la même raison.

Si vous voulez malgré tout les trames brutes, elles ne vivent que sur le poste de calcul : il faut les
copier à la main, par disque ou par OneDrive.
