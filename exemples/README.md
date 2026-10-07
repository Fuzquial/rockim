# exemples/ — exemples prêts à lancer

Chaque sous-dossier est un exemple autonome : un maillage, un deck commenté, un script de figure
et la figure du run de référence. Ils servent de point de départ pour écrire un nouveau calcul.

| Exemple | Mode | Objet | Contenu |
|---|---|---|---|
| [`barre_encastree/`](barre_encastree/README.md) | `fdem3d` | barre de granite 100 × 20 × 20 mm encastrée à gauche, pression de 5 MPa à droite ; le plus petit montage écrit avec les charges et conditions aux limites par groupes | `barre.msh` (3 673 tétraèdres), `barre.cfg`, `fig_barre.py`, `barre_resultat.png` |

![Barre encastrée, résultat du run de référence](barre_encastree/barre_resultat.png)

Lancement, depuis la racine du dépôt :

```sh
./build/rockim exemples/barre_encastree/barre.cfg exemples/barre_encastree/out
python3 exemples/barre_encastree/fig_barre.py exemples/barre_encastree/out
```

Pour d'autres cas, la suite de vérification (`configs/verify_*.cfg`, voir
[`configs/README.md`](../configs/README.md)) donne des decks courts à réponse connue.
