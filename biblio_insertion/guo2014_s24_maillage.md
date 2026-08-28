# Guo 2014, §2.4-2.5 — sensibilité à la taille de maille du modèle de rupture 3D

*Fiche du 2026-08-28. Source : thèse Guo (Imperial College), §2.4 « Mesh size
sensitivity » + §2.5 « Computational efficiency », p. 79-100, lues intégralement
sur PDF fourni par Fernando.* **[V]** = lu dans ce PDF ; **[D]** = dérivation
personnelle explicite.

## 1. Le dispositif **[V]**

Fissure préexistante de 40 mm au centre d'un domaine 120×120×20 mm, pression
interne en rampe P = 1e10·t Pa/s sur les lèvres (éq. 2.61). Demi-modèle avec
rouleaux sur le plan de symétrie — légitime ICI : une seule fissure plane de
mode I DANS le plan de symétrie (rien à voir avec nos motifs de fissuration,
où le quart-modèle reste banni). Matériau type béton/roche tendre :
ρ 2340, E 26 GPa, ν 0,2, ft 3 MPa, c 15 MPa, φ 30°, μ 0,6 (Table 2.3).
Deux énergies : G_f = 10 et 50 N/m. Cinq maillages STRUCTURÉS h = 20/10/5/2,5/
1,25 mm (N = 432 → 1 769 472, Table 2.4) — structuré À DESSEIN : « to
eliminate the influence of different mesh orientations », donc pour isoler
l'effet de TAILLE. La sensibilité d'ORIENTATION est discutée à part (fin de
§2.4.3) : les fissures suivent les frontières d'éléments ; assez petit devant
le domaine, l'écart local ne fausse pas le trajet global. Nos deux règles
maison sont donc complémentaires, pas contradictoires : étude de taille =
structuré admissible ; étude de MOTIF = non-structuré obligatoire
(FICHE 2026-08-06).

## 2. Les résultats **[V]**

Zone plastique en tête de fissure (mesurée du vrai front au max de σ_yy) :

| h (mm) | G_f = 10 : zone plastique | verdict Guo |
|---|---|---|
| 20 | = 1 élément | pilotée par la maille, « only qualitative » |
| 10 | ~2 éléments, floue | pas encore libre de la maille |
| 5 | ~3 éléments | **libre de la maille** (longueur constante de 5 à 1,25) |
| 2,5 | ~6 éléments | idem, champ « nuageux » bien défini |
| 1,25 | ~12 éléments | idem |

À G_f = 50 : h = 20 et 10 pilotés par la maille ; h = 5 : la zone plastique ne
tient plus dans le domaine (120 mm trop court) ; 2,5 et 1,25 corrects.

Charge de rupture (fig. 2.29) : CROÎT avec h — un maillage grossier renforce
structurellement. G_f = 50 : 4,90 MPa (h=20) → 3,09 (h=1,25), écart 1,81 MPa ;
G_f = 10 : 2,63 → 2,10, écart 0,53. Temps de propagation (fig. 2.30) : croît
quand h décroît. **Ni l'un ni l'autre ne converge, même à h = 1,25 mm** — Guo
le dit sans détour : « the smallest element size in this series is not small
enough for a convergent solution », arrêt pour cause de CPU (code SÉRIE :
2,62 s/pas à 1,77 M d'éléments, Table 2.5, Xeon E5-2680 de 2014).

Règle énoncée (§2.4.3) : maillage « fin » = h vaut une FRACTION (un tiers à
un quart) de la taille de la zone plastique en tête de fissure.

## 3. La règle chiffrée pour nos matériaux **[D]**

Guo ne calcule pas la longueur théorique ; en la prenant à la Hillerborg,
ℓ_ch = E·G_f/ft², ses cas donnent ℓ = 28,9 mm (G_f = 10) et 144 mm (G_f = 50 —
supérieure au domaine de 120 mm, ce qui explique exactement son constat à
h = 5). Le seuil « libre de la maille » observé (h = 5 à G_f = 10) correspond
à h ≈ ℓ/6, et sa règle prudente h ≤ ℓ/3 à ℓ/4 est cohérente. Appliquée chez
nous (mode I, le mode II n'est jamais contraignant : ℓ_II = E·G_II/c² =
129 mm pour St Anne) :

| matériau | ℓ_ch = E·G_f/ft² | h ≤ ℓ/3-ℓ/4 | nos maillages |
|---|---|---|---|
| St Anne (E 57, G_f 12, ft 7) | **14,0 mm** | 3,5-4,7 mm | Bgrad fin 0,81 (ℓ/17) OK ; Crebond fin 1,48 (ℓ/9) OK ; **au-delà de r = 15 mm, 3-9,5 mm : HORS règle** — radiales lointaines sous-résolues (déjà assumé au deck) |
| Kuru Grey (E 60, G_I 50, ft 10,98) | **24,9 mm** | 6-8 mm | leur zone fine 1 mm = ℓ/25 ; large |
| Guo test G_f = 10 | 28,9 mm | 7-10 mm | son h = 5 conforme |

Le garde-fou crack-band de rockim (`MatLaw.cpp:1304-1314`, lois dpr/saksala :
throw si le plus gros élément dépasse E·G_f/ft²) est exactement ce critère,
pris à la borne ℓ/1 — il attrape l'inadmissible, pas l'imprécis. Les joints
n'ont AUCUN garde équivalent : la discipline reste à l'auteur du deck (cette
fiche + FICHE 2026-08-06).

## 4. Conséquences pratiques

1. **Comparer à h FIXE.** Charge de rupture et chronologie ne convergent pas
   même à ℓ/23 : toute calibration (calibration_redbohus, benchs) n'a de sens
   qu'à maillage constant entre variantes — jamais « calibré à h1, prédit à
   h2 ».
2. **La charge structurelle n'est pas ft.** Ses charges encadrent les 3 MPa
   du matériau selon G_f (au-dessus à 50, en-dessous à 10) : ne jamais lire
   une résistance de structure comme une résistance de matériau.
3. **Crebond/radiales** : le verdict quantitatif au-delà de r = 15 mm est
   directionnel (comptage, azimuts), pas métrique (longueurs). C'était déjà
   la réserve du deck ; elle a maintenant sa source.
4. Le gros run (fin 0,7-1 mm sur St Anne, ℓ/14-ℓ/20) est confortablement
   dans la règle partout où les fissures comptent.
