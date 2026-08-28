# Bilan du run « Crebond » — dépistage 3D du cycle charge–rebond

**Date : 28 août 2026** · Cas : impact sphère R = 8,51 mm, 0,53 kg à 9,5 m/s (23,92 J) sur bloc 100×100×55 mm, St Anne, maillage gradué 54 010 tets · Verdict : **VALIDÉ PARTIELLEMENT** (périmètre restreint, voir §3)

## 1. Contexte

Le run Crebond (deck du 2026-08-28) devait couvrir pour la première fois un cycle de contact complet — charge, rebroussement, décharge, séparation — afin de dépister le mécanisme décrit par Yang et al. (IJRMMS 206, 2026) pour St Anne : radiales amorcées en charge, écaillage pendant le rebond. Les screenings antérieurs (A, Bg, Bgrad) s'arrêtaient tous avant le rebroussement. Le dimensionnement reposait sur une ODE à loi de contact F = 15,5·p^0,63 kN calée sur Bg, prédisant rebroussement à ~310 µs (p_max 1,76 mm, F_max 22 kN), séparation à ~445 µs et restitution e ≈ 0,37, d'où T = 550 µs.

Le run a été interrompu à t = 204,86 µs (37,2 % de T) par un redémarrage du conteneur de calcul, après environ 8 h. L'examen des pièces est formel : dernière ligne de history.csv complète et saine, vitesse de l'outil strictement monotone, forces bornées, détecteur de NaN silencieux, journal tronqué après le jalon 30 % sans aucun message d'erreur. L'arrêt est exogène ; aucune divergence numérique.

## 2. Méthode de dépouillement

Quatre analyses indépendantes ont été menées sur les mêmes pièces (history.csv, 753 lignes de données ; 11 sorties de champs jusqu'à 196,4 µs ; deck ; journal ; code source) : cinématique de l'outil, morphologie du dommage, conformité aux prédictions du deck, santé numérique. Leurs conclusions ont ensuite été soumises à une contre-expertise contradictoire chargée de les réfuter. Les contradictions chiffrées résiduelles ont été tranchées par retour aux pièces brutes ; celles qui ne pouvaient l'être sont déclarées au §6.

## 3. Verdict en trois volets

**(i) Phase de charge.** La dynamique côté outil est établie et irréprochable : impulsion ∫F dt = 2,016 N·s contre mΔv = 2,017 N·s (0,05 %) ; travail 15,323 J égal à la perte d'énergie cinétique de l'outil à 8·10⁻⁶ près ; forces latérales ≤ 11 % de Fz. La signature mesurée — montée quasi linéaire (~12,6 kN/mm) jusqu'à p ≈ 1,25 mm puis plateau à 16,0 ± 0,8 kN pendant au moins 65 µs, pic brut 18,47 kN à 153,4 µs — réfute la loi de dimensionnement (voir ii). En revanche, **l'identification mécanistique (broyage cisaillant de type coulomb) n'est pas validée** : la contre-expertise a établi, sans réponse à ce jour, que le champ de rupture est suspect d'artefact de sur-pulvérisation (6,49 cm³ détachés à 37 % du run, soit ~6 fois le cratère attendu par le deck ; fragments médians de 0,355 mm³, soit un tétraèdre ; énergie spécifique implicite de 2,36 MJ/m³, vérifiée sur la colonne specificEnergy = work/detachedVol, environ deux ordres de grandeur sous le forage percussif en roche dure) et que la fraction de 78,6 % de ruptures en cisaillement dans le noyau est une signature nécessaire mais non discriminante sous poinçon confiné (marqueur purement diagnostique, classement par ratio d'endommagement). Ce volet est **suspendu** jusqu'à instruction de l'artefact au banc dédié.

**(ii) Objectifs de cycle complet.** Intégralement **non testés** : aucune donnée au-delà de 204,86 µs, donc ni rebroussement, ni décharge, ni séparation, ni rebond, ni écaillage. Mais la fenêtre couverte suffit à réfuter trois prédictions du deck : la loi F(p) (force mesurée 20-35 % sous 15,5·p^0,63, exposant effectif 0,86-1,00 puis saturation) ; p_max = 1,76 mm (l'arrêt dans les 0,183 mm restants exigerait 47 kN moyens, 2,5 fois le maximum jamais observé) ; le rebroussement à 310 µs (il exigerait ~29 kN moyens). Extrapolations concordantes depuis l'état mesuré : rebroussement 340-410 µs, p_max 1,96-2,15 mm, séparation ~480-540 µs. **T = 550 µs était sous-dimensionné même sans l'interruption** ; un successeur doit viser T ≥ 650-700 µs.

**(iii) Chaîne de dépouillement.** Opérationnelle pour la cinématique de l'outil (bilans exacts ; un artefact du lissage zéro-paddé de rapport.py, qui fabriquait une fausse chute de force de −5,6 kN en fin de série, a été identifié et corrigé) et pour la reproductibilité morphologique (volume détaché retrouvé au mm³ près par deux méthodes indépendantes). **Non opérationnelle pour le bilan d'énergie** : les termes toolWork_, biasW_, bcWork_ et l'énergie cinétique de la roche ne sont pas exportés dans history.csv, laissant un résidu de +8,32 J (57 % de l'apport au dernier frame) attribué par simple analogie au biais d'intégration — infermable sur pièces pour un run interrompu. Le poste eFric = −55 J a été tranché sur code : sous-ensemble tangentiel de eGc (net du contact interne : −3,5 J), non un canal additionnel ; mais les flux internes de ±51-55 J (3,4 fois l'apport) rendent les canaux eGc/eFric individuellement inexploitables. Enfin, les champs s'arrêtent à 196,4 µs alors que history court jusqu'à 204,9 µs : les 8,5 dernières µs (1 518 ruptures) échappent à toute analyse spatiale.

## 4. Tableau des critères

| Critère | Mesure | Attendu | Statut |
|---|---|---|---|
| Stabilité numérique, cause d'arrêt | dt constant 3,784 ns ; fin de fichier saine ; log sans erreur | Arrêt exogène, pas de divergence | PASS |
| Hypothèses d'entrée | KE0 = 23,9163 J exact ; 54 010 tets ; absorbantes actives ; contact 11,99 µs | Conformes au deck | PASS |
| Dynamique de charge côté outil | Bilans 0,05 % ; plateau 15,99 ± 0,84 kN (140 µs→kill) ; F_max brut 18,47 kN à 153,4 µs | Bilans fermés, charge axisymétrique | PASS |
| Loi F = 15,5·p^0,63 kN | Mesuré 20-35 % sous la loi ; exposant 0,86-1,00 puis saturation ~16 kN | Loi suivie | **FAIL** |
| p_max ≈ 1,76 mm | p = 1,577 mm au kill avec 8,59 J restants ; arrêt exigerait 47 kN moyens ; extrapolé 1,96-2,15 mm | 1,76 mm | **FAIL** |
| Rebroussement ≈ 310 µs | Exigerait ~28,7 kN moyens (1,55× le max observé) ; extrapolé 340-410 µs | 310 µs | **FAIL** |
| F_max ≈ 22 kN | Max observé 18,47 kN puis plateau ; pic du cycle (au voisinage de p_max) non atteint ; 16-26 kN possibles | 22 kN | NON TESTÉ |
| Broyage cisaillant type coulomb | 78,6 % mode II dans le noyau (gradient 25→90 %) mais marqueur non discriminant et champ suspect d'artefact (cf. §5-6) | Signature coulomb validée | **AMBIGU** |
| Forme conique du noyau | Bol quasi hémisphérique R = 14,4 mm (rms 0,31 mm) contre cône (rms 0,63 mm) | Cône | AMBIGU |
| Radiales amorcées en charge | Structure bipolaire 55°/255° au-dessus du bruit, mais orientation des plans isotrope (28,5 % vs 30,3 %), pas de bras discrets | Amorces radiales | AMBIGU |
| « L'écaillage attend la décharge » | 3,2 cm³ éjectés de r = 10-30 mm dès la charge — mais éjection elle-même suspecte d'artefact | Éjection en décharge seulement | AMBIGU |
| Cycle complet (décharge, séparation, e ≈ 0,37, écaillage en rebond) | 0 µs couverte au-delà de 37,2 % | Cycle observé (objet du cas) | NON TESTÉ |
| Confinement / pollution de bord | Ruptures à ≥ 22,6 mm des faces libres ; eLys = 0,18 J (1,15 %) ; réserve : halo à 28,5 mm > zone fine 15 mm | Phénomènes confinés | PASS |
| Bilan d'énergie fermable sur pièces | Résidu +8,32 J ; termes non exportés ; canaux internes ±51-55 J inexploitables | Bilan fermé | **FAIL** |
| Chaîne cinématique et morphologique | Bilans exacts ; artefact de lissage corrigé ; volumes reproduits au recalcul ; champs en retard de 8,5 µs sur history | Chaîne fiable | PASS |

## 5. Chiffres clés (avec méthodes)

- État au kill : t = 204,86 µs ; p = 1,577 mm (= 0,06351 − toolZ, référence contact) ; v = −5,694 m/s ; F = 15,74 kN ; énergie cinétique restante de l'outil 8,59 J.
- Bilans outil : ∫F dt (trapèzes) = 2,016 N·s vs mΔv = 2,017 N·s ; work = 15,323 J = KE0 − KE_fin exactement.
- Loi de contact mesurée : branche montante (0,15-1,24 mm) F ≈ 12,6·p^1,00 kN (fit log-log) ; fenêtre étendue au plateau : F ≈ 12,0·p^0,86 kN ; saturation à ~16 kN au-delà de p ≈ 1,25 mm ; pic brut 18,47 kN.
- Extrapolations depuis l'état mesuré (intégration de m·dv/dt = −F(p) sous trois lois encadrantes) : rebroussement 340-410 µs, p_max 1,96-2,15 mm.
- Ruptures : 32 712 / 105 996 joints (30,9 %) au kill, taux encore en accélération ; 61,4 % en cisaillement global, 78,6 % dans le noyau (r < 8 mm, prof. < 10 mm) ; enveloppe hémisphérique R = 14,4 mm ; halo jusqu'à r = 28,5 mm.
- Volumes (frame 10, 196,4 µs) : broyé 11 559 mm³ (voxels 1,5 mm, positions courantes ; 7,2-13,5 cm³ selon convention) ; détaché 6 492 mm³ (history et recalcul indépendant par champ fragment, identiques) ; « bourrelet » 1 807 mm³ dont 73 % d'éjecta en vol ; creux net 70 mm³ ≈ calotte de l'indenteur (59 mm³).
- Énergie spécifique implicite : 2,36 MJ/m³ (colonne specificEnergy, vérifiée = work/detachedVol) contre ~150-500 MJ/m³ en forage percussif réel. **Annotation post-contre-expertise** : ce rapprochement mélange deux définitions — la MSE du forage divise par le volume EXCAVÉ, pas par le volume détaché (fragments séparés par joints rompus, majoritairement encore en place dans le bol compacté). Rapportée au creux net (70 mm³), la même énergie donne 219 MJ/m³, dans la fourchette réelle. Le « deux ordres de grandeur » est donc un artefact de définition, pas une mesure de sur-tendreté ; en revanche les deux autres pièces de la suspicion (granulométrie à l'échelle de la maille, volume détaché 6× le cratère attendu) tiennent telles quelles, et la colonne specificEnergy reste inexploitable en l'état.
- Canaux d'énergie au kill (J) : eEl −4,86 ; eJnt −12,77 ; eGc −3,51 (dont part tangentielle eFric −55,03 et part normale +51,52) ; eLys −0,176 ; résidu du bilan partiel +8,32 J au frame 10.

## 6. Réserves et contradictions

**Tranchées sur pièces.** Pénétration au kill : 1,577 mm (le 1,68 mm cité par l'analyse morphologique est le déplacement depuis z0, gap de 0,1 mm compris). Premier contact : première force non nulle (124 N) à t = 11,9865 µs ; le « 10,6 µs » cité ailleurs est l'instant géométrique théorique (10,53 µs) ; l'écart de ~14 µm est un retard d'activation du contact adaptatif, et non une « anticipation » comme écrit dans un rapport. Exposant de F(p) : 1,00 contre 0,86 selon que le fit inclut ou non le plateau — les deux sont exacts dans leur fenêtre, la convention doit désormais être déclarée. F_max « lissé » : dépend de la méthode (17,1-17,7 kN) ; seul le brut (18,47 kN) est citable sans convention. eFric : l'inclusion dans eGc est tranchée sur code ; la conséquence retenue est l'inexploitabilité des canaux, non un bogue du solveur.

**Non tranchées.** (1) La part physique du champ de rupture : la sur-pulvérisation (granulométrie à l'échelle de la maille, volumes détachés, énergie spécifique) peut contaminer le noyau broyé, la chronologie d'éjection et la souplesse même de la réponse F(p) ; renvoyée au banc dédié. (2) La fermeture exacte du bilan d'énergie : le résidu de +8,32 J est compatible en ordre de grandeur avec le biais d'intégration mesuré sur un autre run, mais invérifiable ici. (3) Le F_max du cycle (16-26 kN selon le comportement tardif du broyat). (4) L'origine de la structure azimutale bipolaire (anisotropie naissante ou biais de maillage).

## 7. Décisions et suites

1. **Aucun verdict mécanistique avant instruction de l'artefact de sur-pulvérisation** au banc dédié (recoupement avec les bancs du même jour) : granulométrie découplée de la maille, énergie spécifique plausible, volume détaché compatible avec le cratère.
2. **Pas de relance à l'identique, ni de transfert de la preuve à un run héritier non corrigé.** Précision de contexte : le run « P1 » évoqué au dépouillement existe (bench_impact/configs/impact_pulv_coulomb.cfg, en cours sur la machine du doctorant) mais son objet est le mécanisme de couplage pulvérisation→contact (WP6) sur la géométrie insert + taillant — il ne teste PAS le cycle charge-rebond et ne peut donc pas porter cette preuve : la contre-expertise a raison sur le fond. Le successeur de Crebond est un run distinct : tout successeur doit porter T ≥ 650-700 µs, une loi de dimensionnement recalée sur le mesuré (~12·p^0,9 kN avec saturation), et un réexamen de la zone fine (radiales attendues 10-17 mm à la limite de r = 15 mm, halo déjà à 28,5 mm ; marge frontière réelle ~1,75× et non 3×).
3. **Avant tout gros run** : exporter toolWork_, biasW_, bcWork_ et une énergie cinétique de la roche dans history ; mettre en place une reprise sur checkpoint (l'interruption a coûté ~8 h).
4. **Corrections de chaîne** : padding miroir dans le lissage de rapport.py ; conventions figées (pénétration, instant de contact, fenêtres de fit, méthode de V_broyé).
5. **Ne pas exploiter en livrable** : canaux eGc/eFric isolés, volumes détachés absolus, vitesses d'éjection individuelles, specificEnergy comme grandeur physique.
6. **Acquis à verser aux fiches** : signature de charge (linéaire puis plateau ~16 kN, type capacité portante avec broyage entretenu), réfutation chiffrée du dimensionnement ODE hérité de Bg, et le fait que Bgrad annonçait déjà l'écart (−22 % sous la loi) — le glissement était prévisible avant le lancement.