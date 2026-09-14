# rockim contre Yang et al. 2025, calcaire St Anne à 10,66 m/s — 14/09/2026

Confrontation critère par critère du run `out_stanne2025_rock137` aux résultats publiés dans
*Multi-criteria validation of hi-fidelity numerical model of impact breakage*, Yang, Xiang, Naderi,
Wang, Aising, Ugarte, Latham, **Int. J. Rock Mech. Min. Sci. 191 (2025) 106125**, DOI
10.1016/j.ijrmms.2025.106125.

C'est l'article qu'il faut citer ici, et pas celui de 2026 sur le granite de Kuru. Les figures de
`fig_kinetics.py` portent encore des repères « Yang 9 m/s » qui viennent du Kuru : **ils ne
s'appliquent pas à ce run** et sont à ignorer sur les planches de St Anne.

## 1. D'où viennent les chiffres publiés

Trois provenances, de qualité décroissante, et il faut les distinguer :

| Provenance | Fiabilité | Ce qu'elle donne |
|---|---|---|
| Texte de l'article | exacte | chronologie de fissuration, §5.1 p. 9-10 |
| Table 4 p. 6 | exacte | paramètres matériau, déjà dans notre deck |
| **Figures 9 et 10, lues sur le nuage de points** | **±10 % au mieux** | les sept critères en fonction de la vitesse de piston |

Les sept critères ne sont **publiés que sous forme de graphiques** en fonction de la vitesse du
piston. Aucune table de valeurs. Les chiffres ci-dessous marqués « ≈ » sont donc des lectures au
point le plus proche de 10,66 m/s sur les croix bleues (calcaire, FDEM), avec une incertitude que je
ne sais pas réduire sans les données des auteurs.

## 2. Transfert d'énergie (leur figure 9)

| Critère | Yang FDEM ≈ | Yang expérience ≈ | rockim | Écart à Yang FDEM |
|---|---:|---:|---:|---:|
| Contrainte maximale dans le taillant | 200 MPa | 200 MPa | **212,3 MPa** | +6 % |
| Vitesse d'indentation du taillant | 8,0 m/s | 7,9 m/s | **7,34 m/s** | −8 % |
| Enfoncement maximal | 1,10 mm | 1,05 mm | **1,292 mm** | +17 % |
| Vitesse de rebond du taillant | 5,6 m/s | 4,2 m/s | **non mesurable** | — |

**Contrainte et vitesse d'indentation : accord.** Les deux tombent dans la barre de lecture des
figures. La contrainte est prise comme eux, pic de la **première** onde à mi-taillant (212,3 MPa à
34,3 µs), pas le maximum global. La vitesse d'indentation est prise comme eux, **pente** de la portion
linéaire du déplacement entre 10 et 90 %, pas un maximum instantané : le maximum instantané vaut
8,68 m/s pour le taillant et 12,0 m/s pour l'insert, et se comparerait à tort.

**Enfoncement : 17 % de trop.** C'est le premier écart réel. Il va dans le même sens que la masse de
fragments ci-dessous, un cratère plus ouvert que le leur.

**Rebond : hors de portée de ce run.** Le calcul s'arrête à 300 µs. Yang lit sa pente de rebond sur la
remontée établie, et sa figure 8 montre un enregistrement qui court jusqu'à 750 µs. À 298 µs notre
contact porte encore 80 kN et la vitesse oscille entre 1,3 et 2,5 m/s : ce n'est pas une vitesse de
sortie. **Ce critère reste à mesurer**, il faut prolonger à 450 µs au moins.

## 3. Fissuration (leur figure 10)

| Critère | Yang FDEM ≈ | Yang expérience ≈ | rockim | Écart |
|---|---:|---:|---:|---:|
| Longueur de fissure radiale | 32 mm | 31 mm | **29,6 mm** (rayon maximal rompu) | −8 % |
| | | | *17,7 mm* (portée radiale par secteur) | *−45 %* |
| Rayon de cratère | 13,5 mm | 13 mm | **11,8 mm** | −12 % |
| Masse de fragments | 2,5 g | 1,2 g | **4,67 g** | +87 % |

**Rayon de cratère : accord à 12 %**, dans la dispersion de leurs propres points.

**Longueur radiale : dépend de la définition, et c'est un problème.** Ils définissent la longueur
radiale comme *la plus grande distance du centre de la roche à la pointe d'une fissure radiale*.
Nous produisons deux grandeurs candidates, et l'écart entre elles est d'un facteur 1,7 :

- le **rayon maximal des facettes rompues**, 29,6 mm, qui suit leur définition à la lettre mais
  englobe les fissures latérales ;
- la **portée radiale par secteur**, 17,7 mm au maximum sur seize secteurs, qui isole les fissures
  classées radiales mais mesure une portée locale et non la plus grande.

Le rapport longueur radiale sur rayon de cratère vaut **2,4 chez eux et 2,5 chez nous** avec la
première grandeur, ce qui plaide pour elle. Mais l'accord à 8 % qui en résulte ne doit pas être
présenté comme acquis tant que la convention n'est pas tranchée sur leurs images.

**Masse de fragments : 87 % de trop, et c'est l'écart le plus net.** Nous brossons 1 711 mm³, soit
4,67 g avec leur masse volumique de 2 731 kg/m³, contre environ 2,5 g pour leur FDEM et 1,2 g mesuré
au laboratoire. Trois causes possibles, non départagées :
le critère de brossage n'est pas réglé comme le leur ; le cratère est réellement plus ouvert, ce que
l'enfoncement en excès de 17 % corrobore ; ou notre volume détaché compte des éléments qu'un brossage
physique ne détacherait pas. Les auteurs signalent eux-mêmes que leur FDEM dépasse l'expérience de
0,85 g environ à haute vitesse, par recompactage de la poudre sous l'insert.

## 4. Chronologie (leur §5.1, texte, donc exact)

| Événement | Yang | rockim | Écart |
|---|---:|---:|---:|
| Zone broyée en place | avant 130 µs | premières ruptures à 60,9 µs | — |
| Fin de la phase de charge | **254 µs** | **284,3 µs** | **+12 %** |
| Enfoncement maximal | — | 283,1 µs | — |
| Arrêt des fissures médianes | 291 µs | non atteint | — |
| Arrêt des radiales | 388 µs | non atteint | — |
| Arrêt des latérales | 482 µs | non atteint | — |

**La fin de charge est notre meilleur repère temporel : 284,3 contre 254 µs, soit 12 % de retard.**
C'est cohérent avec l'enfoncement plus grand : on charge plus longtemps et plus profond.

⚠️ **Comment ce chiffre a failli être faux.** La vitesse de l'insert croise zéro une première fois à
260,3 µs, puis repart vers le bas et l'enfoncement reprend jusqu'à 283 µs. Le premier croisement est
une **pause sous charge maintenue**, pas un retournement : la force de contact ne retombe pas, elle
reste entre 86 et 90 kN de part et d'autre. Un croisement de vitesse ne suffit pas, il faut vérifier
la décharge. `tools/fig_retournement.py` porte désormais la force de contact sur la figure et marque
les pauses d'une croix rouge.

## 5. Ce que la comparaison établit, et ce qu'elle ne peut pas établir

**Établi.** Sur les quatre critères mesurables ici — contrainte, vitesse d'indentation, rayon de
cratère, date de fin de charge — rockim tombe entre 6 et 12 % de Yang. C'est du même ordre que la
dispersion de leurs propres points expérimentaux et que leur coefficient de variation annoncé de 10 %
sur la fissuration. La morphologie est la bonne : radiales dans les seize secteurs, réseau à 62 % de
faces sub-verticales, 47 % de ruptures en traction.

**Non établi, par ordre de gravité.**

1. **La vitesse de rebond n'est pas mesurée** : le run s'arrête à 300 µs au lieu des 450 nécessaires.
   C'est le critère le plus sensible au contact, donc celui qui manque le plus.
2. **La masse de fragments est à 87 %** et l'enfoncement à 17 %. Les deux vont dans le même sens, un
   cratère trop ouvert. À instruire ensemble, pas séparément.
3. **La longueur radiale dépend d'une convention non tranchée**, avec un facteur 1,7 entre nos deux
   définitions candidates.
4. **Les valeurs publiées sont lues sur des graphiques** à ±10 %. Tous les accords annoncés ci-dessus
   sont à prendre avec cette barre. Les auteurs annoncent leurs données disponibles sur demande : une
   demande ciblée des sept critères en table vaudrait mieux que toute relecture de figure.
5. **Une seule résolution de maillage.** Notre cœur est à 1,42 mm d'arête médiane, le leur à 1 mm
   annoncé. `rock073` (1,03 mm, 251 460 tétraèdres) reste à jouer pour savoir ce qui, dans les écarts
   ci-dessus, est du maillage.

## 6. Le run est terminé (12 h 05) — valeurs finales et bilan d'énergie

29 h 21 de calcul (105 690 s) sur 14 fils, 22 trames, arrêt propre à 300,001 µs.
Rien n'a été enchaîné derrière : la file de nuit a été désarmée à la demande.

**Le bilan d'énergie ferme, et c'est le résultat le plus solide du run.**

| | |
|---|---|
| Énergie cinétique, début → fin | 60,12 → 3,28 J |
| Résidu du budget B4 | −9,8e-8 J, soit **1,5e-7 %** de l'échelle |

Après une semaine passée à traquer des lois de joint qui créaient de l'énergie, celle-ci n'en crée
pas, sur 300 µs et 12 090 facettes rompues. C'est la validation que les bancs de cycle fermé
annonçaient.

**Les postes de dissipation**, à confronter à `yang2024energy` :

| Poste | Valeur |
|---|---:|
| Rupture (Gc) | 34,5 J |
| Frottement | 27,8 J |
| Joints (élastique + adoucissement) | 14,5 J |
| Élastique stocké | 18,1 J |

Le frottement pèse enfin quelque chose : 27,8 J, contre les **0,66 J** mesurés sur les anciens bancs
Kuru face aux 32 J publiés. C'est l'écart que la mémoire du projet identifiait comme le déficit
principal, et il n'existe plus sur ce montage. ⚠️ Comparaison à prendre avec précaution : ce sont deux
roches, deux vitesses et deux decks différents, et le chiffre de 32 J est celui du Kuru.

**Valeurs finales des critères** (trame 21, 300 µs), qui ne changent pas les conclusions du §2 et §3 :

| | 285 µs | 300 µs |
|---|---:|---:|
| Facettes rompues | 11 914 | 12 090 |
| Fragments | 3 904 | 4 089 |
| Rayon de cratère, peau | 11,84 mm | **12,00 mm** |
| Rayon maximal rompu | 29,56 mm | 29,60 mm |
| Volume détaché | 2 139 mm³ | 2 153 mm³ |
| Masse brossée | 4,67 g | **4,70 g** |
| Radiales, portée moyenne | 8,79 mm | 8,66 mm |

**Un fait nouveau dans les quinze dernières microsecondes : la profondeur du réseau passe de 36,9 à
40,9 mm**, +11 %, alors que tout le reste est figé. Une fissure descend pendant la décharge. C'est
cohérent avec Yang, qui fait cesser les médianes à 291 µs, donc après la fin de charge. À vérifier
sur les coupes avant d'en faire un résultat.

**Le rebond n'est toujours pas mesurable** : à 300 µs l'insert repart à −0,33 m/s pendant que le
taillant monte à 1,47 m/s et le piston à 1,07. Le train se décompose, il n'est pas sorti. Les énergies
cinétiques finales valent 1,84 J pour le taillant et 1,36 J pour le piston, contre 60,1 J au départ.

## 7. Le montage, pour mémoire

Deck `configs/stanne2025_rock137_visc0.cfg`, binaire `rockim_g1y19.exe`, maillage
`impact_yang_train1_rock137_hxt.msh` (108 667 tétraèdres, cœur 14 030 à 1,42 mm).
Matériau Table 4 : ρ 2731, E 57 GPa, ν 0,31, ft 7,0 MPa, c 18,8 MPa, tan φ 1,0, GI 12 J/m²,
GII 800 J/m². Pénalité `jointPenaltyLength = edge` avec le facteur 26,316 = p0/(2E) et les 3 000 GPa
de l'ARMA 24-0952 pour St Anne. Loi `plastic` + plage coulomb + `jointSecantRatchet`.
**Pulvérisation désactivée** : c'est le modèle 2026 du granite, absent de l'article 2025.
Sans viscosité. Forces de contact mesurées par `contactForcePairs`.

Reproduction : `docs/REPRODUIRE_stanne_radiales_2026-09-14.md`.
Figures : `results/fig/stanne285/` et `output/pdf/stanne_285us/`, index dans
`docs/FIGURES_pour_le_rapport_2026-09-13.md`.
