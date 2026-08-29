# Guo 2014, §2.3.4-2.3.5 — contact, frottement, intégration explicite, pas de temps

*Fiche du 2026-08-28. Source : thèse Guo (Imperial College), p. 75-79, lues
intégralement sur PDF fourni par Fernando.* **[V]** = lu dans ce PDF ;
**[D]** = dérivation/mesure personnelle.

## 1. §2.3.4 — contact et frottement **[V]**

- Détection : **NBS** (Munjiza & Andrews 1998), mémoire et CPU linéaires en N
  (éq. 2.47-2.48). Interaction : **pénalité par potentiel** (Munjiza &
  Andrews 2000), force distribuée f_contact = ΣΣ ∫ (grad φ_c − grad φ_t) dV
  (éq. 2.49), paires contacteur/cible de TÉTRAÈDRES, forces réparties aux
  nœuds. La rugosité microscopique des surfaces est explicitement ignorée.
- **Frottement — LE point qui clôt une recherche** : la thèse n'y consacre
  que TROIS PHRASES : « Sliding friction is also considered as a type of
  contact. A Coulomb friction law was implemented into the three-dimensional
  FEMDEM code by Dr Jiansheng Xiang. Sliding [...] will occur when the
  tangential contact force f_tan is greater than μN. » AUCUN détail
  d'algorithme tangentiel : pas de ressort k_t, pas de 2/7, pas de
  régularisation. L'implémentation est attribuée à Xiang, non publiée ici.
  > **MISE A JOUR DU 2026-08-29 (LOT 2c) — LE TROU EST COMBLE, ET AILLEURS.**
> L'algorithme tangentiel EST publie, non pas dans la these ni dans les articles
> d'impact, mais dans Xiang, Latham & Farsi (2017), « Algorithms and Capabilities
> of Solidity », eq. (4)-(5) p. 4 : f_t = -k_t delta_t - eta v_t en regime
> adherent, bascule a f_t = -mu f_n des que f_t depasse mu f_n. Il y a donc bien
> un RESSORT TANGENTIEL et un amortisseur visqueux, comme dans rockim. Le
> chapitre attribue l'implementation a « Xiang et al (2009) », ce qui boucle avec
> l'attribution a Xiang de la these. La conclusion ci-dessous — « la source de
> verite est le code Solidity lui-meme » — est donc CADUQUE : la source de verite
> est une publication. Voir
> [la fiche du lot 2c](2026-08-29_lot2c_frottement_tangentiel.md).

**Conséquence** : la vérification du relais joint→contact et du
  k_t = 2/7 de rockim ne peut PAS se faire contre la thèse — la source de
  vérité est le code Solidity lui-même (Y3Dfd.c, déjà cloné dans le
  conteneur le 2026-08-27). Inutile de chercher d'autres pages de thèse
  pour ça.

## 2. §2.3.5 — équations et intégration **[V]**

- Bilan nodal : m·v̇ + f_int = f_ext, avec f_ext = f_joint + f_contact +
  f_load (éq. 2.50-2.51).
- **Intégration : Euler explicite SEMI-IMPLICITE (symplectique)** —
  éq. 2.53-2.55 : v̇ = (f_ext − f_int)/m ; v_{t+1} = v_t + v̇Δt ;
  x_{t+1} = x_t + **v_{t+1}**Δt (la position avance avec la vitesse NEUVE).
  **Aucun amortissement numérique dans l'intégrateur** — le point que la
  citation partielle du §2.3.5.2 laissait à sceller est scellé : toute
  dissipation vient des joints, du frottement et des amortisseurs
  explicites, jamais du schéma. rockim utilise EXACTEMENT le même schéma,
  pour les nœuds comme pour l'outil rigide
  (`Tool3::integrate : v += (dt/m)F ; x += dt·v`).

## 3. Pas de temps (éq. 2.56-2.60) **[V]**

- FEM : Δt_FEM ~ 0,1·h·√(ρ/E), h = plus petite ARÊTE (un dixième du transit
  d'onde dans l'élément, Kolsky).
- DEM : Δt_DEM ~ (π/5)·√(m/k) avec k = E·h ; pour un tet régulier
  (m = ρh³/6√2) → Δt_DEM ≈ 0,2·h·√(ρ/E).
- Δt = min des deux (éq. 2.58) → **c'est la condition FEM qui gouverne**.

**[D]** Contrôle sur nos runs : grad_b (h_min = 0,459 mm, St Anne) donnerait
Δt_Guo = 0,1×4,59e-4×√(2731/57e9) ≈ 1,0e-8 s ; rockim (dtFactor 0,15 sur sa
propre formule interne) a pris 1,93e-9 s — un facteur ~5 SOUS la règle de
Guo. Notre marge de stabilité est donc large ; si un jour le coût CPU des
gros runs devient bloquant, il y a de la marge documentée pour relever
dtFactor — après un test de bilan d'énergie, jamais à l'aveugle.

## 4. Ce que ça change à la liste de partage

Le §2.3.5.2 est acquis en entier (plus besoin), et la section frottement
n'existe pas dans la thèse au-delà des trois phrases ci-dessus : la
vérification tangentielle de WP6 se fera contre Y3Dfd.c. Reste UN SEUL item
de thèse en attente : la table des matières des chapitres 5+ (validations
dynamiques éventuelles à transformer en benchmarks).
