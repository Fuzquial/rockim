**Diagnostic indépendant de la reproduction d’ICL dans Rockim — 12 septembre 2026**

**Complément documentaire après lecture du corpus Yang :** voir `COMPLEMENT_YANG_sources_et_corrections_2026-09-12.md`. La comparaison 7,37/5,62 m/s ci-dessous doit être nuancée : Yang 2025 mesure une pente linéaire du déplacement, alors que le chiffre Rockim est un maximum instantané. Ce rapprochement est un signal, pas une preuve autonome d’erreur de chargement. La pénalité 3000 GPa retrouvée dans ARMA concerne St Anne, pas une valeur publiée pour le Kuru. Le complément propose un impact St Anne sans pulvérisation comme premier cas discriminant.

Le répertoire examiné est `FDEM/rockim_g1`. Conclusion : les données ne démontrent pas une incapacité intrinsèque du FDEM 3D à localiser. Elles montrent une reproduction encore non validée, avec une loi de joint hybride, une pulvérisation pratiquement inactive, un chargement différent de la référence et des erreurs de post-traitement. Raffiner ou diminuer Gf au hasard ne constitue pas une correction démontrée.

Périmètre : code source présent, configuration effectivement consommée et résultats de `out_yang2026_v3` (maillage s=1, historique disponible jusqu’à 182,925 µs, VTU jusqu’à 179,988 µs), quelques configurations et journaux s=2,5, article ICL et chapitre 2 de Guo. Les commentaires des configurations et les notes antérieures sont traités comme des affirmations à vérifier. La présence de nouvelles options dans g1y16 ne prouve pas qu’elles étaient actives dans le calcul v3. Aucun calcul complet supplémentaire n’a été exécuté ; aucun solveur ni calcul existant n’a été modifié.

**1. Ce qui est effectivement calculé**

Mesures recalculées depuis les CSV et VTU, sans reprendre les nombres des notes :

| Grandeur | Résultat Rockim | Référence / signification |
|---|---:|---|
| Vitesse maximale descendante, groupe `bit` | 7,369 m/s | ICL p. 11 : 5,62 m/s à piston 9 m/s ; vérifier la même définition du corps et de la vitesse |
| Compression maximale à la jauge du bit | 175,95 MPa | Un bon pic d’onde ne suffit pas à valider toute la cinématique |
| À environ 100 µs | 623 joints rompus ; `bulkD=0` partout dans la roche | Rupture des interfaces sans pulvérisation des tétraèdres |
| À 179,988 µs | 3 357 joints rompus ; seulement 10 tétraèdres avec `bulkD>0`, maximum 0,25883 | Aucun tétraèdre n’a atteint 0,9 à cette trame |
| À la dernière ligne CSV, 182,925 µs | 3 530 joints rompus ; 2 éléments pulvérisés ; bit encore descendant à 5,730 m/s | État incomplet, pas de rebond final mesurable |
| Maillage total | 120 185 tétraèdres | ICL : 230 788, §3.1 |
| Dans la demi-boule de rayon 12,5 mm | 14 722 tétraèdres ; médiane de la moyenne des 6 arêtes : 1,372 mm | La consigne gmsh « 1 mm » n’est pas une mesure des arêtes effectivement produites |
| Même zone, diamètre inscrit médian | 0,51035 mm | C’est cette longueur qui pilote actuellement la pulvérisation |

La légère différence entre les nombres du CSV à 179,936 µs et du VTU à 179,988 µs vient des dates d’échantillonnage différentes.

Les coupes exactes à 100 et 180 µs confirment un réseau dense sous l’insert, avec quelques extensions périphériques. À 180 µs, les centroïdes des facettes effectivement rompues atteignent un rayon de 9,16 mm et une profondeur de 10,14 mm ; 95 % sont à moins de 6,67 mm de l’axe. **Le rayon maximal de cette population n’est pas la longueur d’une fissure radiale identifiée.** Une coupe à y=0 ne remplace pas non plus une coupe suivant la plus longue fissure, comme dans ICL.

L’article distingue les joints rompus (fissures mésoscopiques) et `D` dans les tétraèdres (pulvérisation). Sa figure 14 à 9 m/s comporte un noyau compact, des fissures périphériques limitées, et une évolution jusqu’à 543 µs ; le retournement est à 255 µs. La figure 15 à 11 m/s présente davantage d’éjection et de fissures périphériques. Il ne faut ni attendre un grand cône médian dans ce granite, ni déclarer le résultat conforme à partir de son seul diamètre. ICL décrit justement l’absence de la grande zone cisaillée profonde caractéristique d’autres roches.

**2. La loi active n’est pas la loi complète de Guo**

`out_yang2026_v3/config_effective.cfg` contient `jointShearUnload = plastic`, `jointShearRange = coulomb`, `jointElastic = parabolic`, `jointSecantRatchet = on`.

Dans `src/Fdem3dSolver.cpp:4868`, le cisaillement utilise :

```text
tau_trial = pj * (s - s_plastique)
tau = projection sur tau_lim
s_plastique += (tau_trial - tau) / pj
D_modeII piloté par |s_plastique| / slipRef
```

Guo, équations 2.25–2.33, pp. 67–70, prescrit une relation traction–déplacement avec une branche parabolique avant le pic, puis un adoucissement dépendant des déplacements normal et tangentiel. La substitution du cisaillement à l’ouverture est explicitement indiquée p. 70.

Même avant l’endommagement, l’écart est vérifiable : avec `sp=fs/pj` et `s=sp/2`, la branche active donne `tau=fs/2`, la parabole donne `tau=3fs/4`. La pente initiale tangentielle vaut donc `pj` dans cette branche, contre `2pj` pour la parabole. Le choix `jointElastic=parabolic` ne rend pas le retour plastique tangentiel parabolique. La variable d’histoire et la décharge sont également différentes.

**Conséquence certaine :** copier ft, c, GI et GII de l’article dans cette loi ne reproduit pas automatiquement son comportement. **Hypothèse mécanique à tester :** les interfaces peuvent accommoder trop de déplacement et perdre leur portance avant que les tétraèdres développent les déformations de pulvérisation ; l’insert pénètre dans un réseau de joints cassés au lieu de produire la séquence noyau pulvérisé / redistribution / fissures périphériques. Les données sont compatibles avec cette hypothèse, mais ne suffisent pas à attribuer tout l’écart à ce seul mécanisme.

La branche `jointShearUnload=solidity` existe maintenant dans le source (`:4580` environ). Elle mérite une comparaison contrôlée. Son nom ne vaut pas validation contre le binaire des auteurs. L’absence de maximum historique explicite dans l’équation 2.33 ne démontre pas, à elle seule, toute la gestion de décharge du logiciel des auteurs. Et Guo p. 69 dit clairement qu’un joint effectivement rompu quitte la loi cohésive : parler d’un joint qui « guérit » doit être limité à la branche pré-rupture, pas aux fissures déjà créées.

**3. La pulvérisation est une approximation encore non identifiée**

Le code (`src/Fdem3dSolver.cpp:4042–4078`) calcule :

```text
h = 6 V0 / somme_des_aires_des_faces
dm = h * sqrt(2/3) * ||dev(epsilon_corotationnelle)||
kappa = maximum historique de dm
D = min(0.9, df*(kappa-d0)/(kappa*(df-d0))) au-delà de d0
sigma *= Cd*(1-D)
```

Deux faits indépendants :

- En compression isotrope, `dev(epsilon)=0` : ce mécanisme ne s’active pas, quelle que soit la contraction volumique. C’est une propriété du modèle codé, pas un défaut de l’affichage.
- Pour le diamètre inscrit médian mesuré, le seuil d’amorçage correspond à environ 2,74 % de déformation équivalente déviatorique. Avec la formule codée, `D=0.9` est atteint à `kappa≈106,46 µm`, soit environ 20,86 % pour cette taille. Ce seuil n’est pas 0,4 mm : le plafonnement à 0,9 intervient avant. La formule rationnelle de D n’est pas une croissance linéaire de D avec le déplacement.

À 100 µs, le déplacement effectif courant maximal reconstruit dans la zone centrale est environ 10,79 µm, inférieur aux 14 µm d’amorçage. À 180 µs, il atteint seulement 14,40 µm ; quatre éléments de cette zone dépassent alors le seuil courant. L’historique exporté de D peut naturellement concerner davantage d’éléments que le critère courant.

Le contrôle utilise la décomposition polaire exacte des positions exportées ; le solveur utilise trois itérations de Higham. Il s’agit d’une vérification indépendante des ordres de grandeur, pas d’une reconstruction bit à bit des contraintes internes.

ICL §2.2, équations 3–4 et figure 4, ne fournit pas une définition opérationnelle assez complète de `delta_m`, de la mesure de déformation et de Cd pour affirmer que cette implémentation est identique. La notation de l’équation 4 et son explication comportent en outre une ambiguïté entre déplacement courant, maximal et final. L’équation 3 publiée ne se réduit pas de façon univoque au multiplicateur adimensionnel `Cd*(1-D)` utilisé ici.

**Correction nécessaire pour une reproduction stricte :** obtenir la définition algorithmique exacte de `delta_m`, sa longueur de référence, la loi contrainte–déformation après amorçage, Cd et epsilon_d. Tester alors un élément en compression isotrope, compression uniaxiale et cisaillement avant de relancer l’impact. Diminuer arbitrairement `bulkDamageDelta0` peut activer D, mais constitue une recalibration, pas une réparation démontrée.

**4. Chargement et maillage sont aussi en cause**

Les volumes initiaux donnent : piston 1,05703 kg ; corps acier `bit` 1,28830 kg ; ensemble bit+insert+circlip 1,36713 kg. ICL p. 11 indique une masse de bit de 1,509 kg. Il faut vérifier le périmètre exact de cette masse. L’énergie cinétique initiale du piston maillé à 9 m/s vaut 42,81 J ; ce n’est pas l’énergie que reçoit le bit.

Le journal du jumeau s=2,5 donne un piston de 0,776709 kg : **26,5 % plus léger que le maillage s=1**. La discrétisation grossière modifie aussi les masses et la géométrie des corps. Comparer ces deux runs comme une étude de résolution de la roche uniquement est donc invalide. Pour une convergence utile : garder le train métallique, ses masses et son maillage fixes ; raffiner seulement la roche et contrôler aussi la surface de contact de l’insert.

Le choix de longueur de pénalité est un autre écart : le chemin par défaut des maillages fichiers utilise `hmin_`, alors que Guo p. 67 utilise la moyenne des arêtes de la facette. `jointPenaltyLength=edge` corrige cette convention dans g1y16. Avec la parabole, la relation est `p0 = 2*jointPenaltyFactor*E`. Ainsi le facteur 25 correspond à `p0=50E`, qui n’est pas la plage E–10E donnée par Guo ; cela doit être justifié par la référence de l’impact ou exploré comme sensibilité, pas présenté comme une valeur publiée dans la Table 1 d’ICL.

Enfin, le calcul actif garde `bulkModel=corotational` et le plafonnement caché `meanTensionCapFactor=3` : la pression moyenne en traction est écrêtée à 32,94 MPa dans la roche (`:4038`). ICL utilise une loi néo-hookéenne avant dommage. Le cap doit être neutralisé (`meanTensionCapFactor=0`) dans un essai de conformité et son effet mesuré. Ce n’est pas une preuve qu’il gouverne le run actuel. Dans la zone centrale des trames examinées, le minimum de det F est voisin de 0,96 : les affirmations antérieures de det F≈0,5–0,7 ne sont pas vérifiées pour cette zone à ces instants.

**5. Deux défauts de lecture des fissures sont établis**

`bench_impact/tools/imp_lib.py:77` sélectionne `damage>=0.999` et `bonded<0.5`. Or `damage` est le maximum des trois points ; sous `majority`, une facette peut avoir ce maximum égal à 1 sans être rompue. `bonded` n’est pas le drapeau `dead`. Le vrai événement de rupture est enregistré dans `tBreak` (`src/Fdem3dSolver.cpp:4916`).

À 100 µs, ce filtre annonce 711 facettes contre 623 ruptures ; à 180 µs, 3 737 contre 3 357, soit 380 faux positifs. La correction appliquée dans le script indépendant est `tBreak>=0`. Pour distinguer une fissure rompue mais fermée d’une fissure ouverte, il faut aussi exporter l’ouverture et/ou `dead` ; la rupture seule ne donne pas l’ouverture géométrique.

De plus, `imp_lib.py:94–95` assimile rayon maximal des centroïdes à longueur de fissure et rayon de cratère. Ce sont des extensions spatiales, pas les métriques expérimentales : des fragments déplacés peuvent gonfler le rayon, un réseau diffus peut atteindre 10 mm sans former une radiale connectée, et une facette cassée en surface ne signifie pas que la matière a été retirée. Identifier des chemins de fissures connectés et un cratère après retrait des fragments.

Réserve sur les couleurs : le moteur plastique de dommage emploie `slipRef`, qui dépend de la pression (`:4879`), mais la classification exportée à la rupture renormalise par `J.slipF` (`:4934`). Les composantes normales et tangentielles sont aussi prises comme maxima sur la facette. L’étiquette traction/cisaillement peut donc différer de celle du point qui a réellement provoqué la rupture. Corriger cette instrumentation avant une comparaison quantitative rouge/jaune avec ICL.

**6. Solution proposée, dans l’ordre**

1. **Rendre les mesures fiables.** Conserver `tBreak>=0`, séparer D volumique, rupture et ouverture, mesurer les fissures connectées. Exporter les forces de contact piston/bit et insert/roche, et les énergies par corps. `toolFz=0` est normal avec `toolShape=none` : ce champ n’est pas la force de l’insert maillé. Ne pas calibrer une réaction sur cette colonne ni assimiler sans vérification `m_bit*dv/dt` à la seule réaction de la roche pendant l’interaction avec le piston.
2. **Valider une loi de joint unique sur un petit banc 3D.** Traction monotone et décharge, cisaillement à plusieurs pressions, trajet à pression variable et relais joint/contact. Vérifier courbes, aire dissipée, ouverture de rupture et passage à deux points sur trois. Comparer la branche actuelle à `solidity` avec le même maillage et la même pénalité. Un essai de traction/flexion entaillé doit produire une fissure connectée avec une réponse convergente avant de parler de capacité ou d’incapacité générale de localisation.
3. **Séparer les sensibilités.** Comparer indépendamment longueur de pénalité, loi de cisaillement, naissance du contact (`gcBirth=ramp` versus `penalty`), puis raideur de contact. Le banc `yang2026_bench_s25_solidity.cfg` change simultanément loi, longueur et facteur de pénalité, raideur et naissance du contact. Il peut tester un paquet de changements, mais ne permet pas d’attribuer son résultat à la seule mémoire de la loi.
4. **Identifier la pulvérisation.** Obtenir les données manquantes des auteurs, ou assumer explicitement une autre loi régularisée et la calibrer sur les essais. Vérifier l’amorçage en compression confinée et le travail dissipé. Un recalcul avec `bulkModel=neohookean` est une comparaison utile, pas une solution suffisante à un dommage presque inactif.
5. **Recaler le train de frappe puis raffiner.** Contrôler masses, impulsions, jauge, vitesse incidente effective et profondeur avec le même périmètre du bit que les auteurs. Conserver les corps métalliques identiques entre niveaux de maille. Étudier au moins trois résolutions de roche ; ne pas considérer 1,372 mm comme « 37 % trop grossier » au sens strict sans connaître la même statistique du maillage des auteurs.
6. **Terminer l’impact.** Pour comparer la figure 14 entière, aller au-delà de 543 µs, idéalement jusqu’à disparition du contact et rebond stabilisé. Pour la masse de débris, ajouter ensuite le tri approprié. Le résultat s=1 actuellement disponible s’arrête avant le retournement de référence. Le prolonger ne résoudra toutefois pas à lui seul l’absence de pulvérisation au début de la charge.

La bit-identité des tests existants garantit une absence de changement sur ces cas, pas la justesse du modèle d’ICL. Le garde-fou `KE<=KE_initiale+travail_exterieur` est utile, mais insuffisant pour prouver la thermodynamique de toutes les branches : il faut aussi des bilans de stockage/dissipation sur les trajets constitutifs.

**Verdict : priorité à la cohérence joint–volume–contact et au transfert de charge.** Le maillage et la durée comptent, mais il n’y a pas de base pour promettre que changer seulement une taille de maille ou un Gf donnera les figures. Le réseau diffus actuel est réel ; il ne prouve pas que l’architecture FDEM 3D est incapable de fissurer. Une reproduction exacte reste conditionnée aux détails non publiés du modèle de pulvérisation et à la validation des corrections par des essais contrôlés.

Sources scientifiques : ICL.pdf, pp. 4–6, 9–12, figures 14–18 ; Guo, chapitre 2, pp. 67–70 et §2.4. Article : [Yang et al., 2026](https://doi.org/10.1016/j.ijrmms.2026.106660).

Fichiers indépendants de vérification : `mesures_ICL_2026-09-12.json`, `coupes_ICL_verifiees_2026-09-12.png` dans ce dossier ; script reproductible `simulations/tmp/icl_audit/analyse.py`. Les images PDF sources sont des intermédiaires d’examen, pas des résultats de simulation.
