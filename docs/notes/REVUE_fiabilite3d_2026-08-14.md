# Fiabilité 3D — revue, diagnostic et plan d'attaque (2026-08-14)

*Objet : l'instabilité latente du 3D en phase débris (chasse gelée le
2026-08-07, FICHE §instabilité). Commande de Fernando : diagnostic d'abord,
revue de littérature et des codes similaires, puis plan d'attaque.*

## 1. Le dossier interne (ce qu'on sait déjà)

**Symptôme** (percussion 3D homogène, grille de Kuhn, T doublé 2e-4 → 4e-4) :
**2,4 MJ créés pour 16 J incidents**, cascade auto-entretenue entre 240 et
400 µs SANS charge extérieure. Les DEUX schémas explosent (l'insertion
adaptative est disculpée). Autopsie : des nœuds ISOLÉS du cratère à
250-300 m/s dès 133 µs (p99 du bloc : 0,5 m/s), population contagieuse,
la pompe vit dans la **phase débris** : petits fragments entre contact
général, joints comprimés survivants et frottement tanh.

**Ce qui a changé depuis le gel** (audit 08-11 → 13) :
- borne MOOSE `cd ≤ m_eff/dt` AJOUTÉE au dashpot de joint 3D (le suspect n°1)
  + dashpot bilatéral sur joint intact + compteur `dampWork_` à verdict ;
- CFL sur le min RÉEL des diamètres inscrits (le nominal surestimait ×5,3
  sur Kuhn jitté) — nécessaire mais démontré insuffisant ;
- **B4 (bilan par sous-système) construit** — c'est précisément l'instrument
  que le plan gelé réclamait ;
- contact par potentiel de Munjiza (conservatif) disponible — au moment de
  l'explosion, le contact général était en pénalité.

**⚠️ Le cas explosif n'a JAMAIS été rejoué depuis ces correctifs.**

## 2. Revue externe (littérature + codes)

**La pièce maîtresse — arXiv:2511.14323 (nov. 2025), « Stability of
Extrinsic Cohesive-Zone Model with Penalty-Based Contact in Explicit
Dynamic Fragmentation »** : analyse de stabilité de NOTRE configuration
exacte (CZM extrinsèque + pénalité + fragmentation). Trois sources
d'instabilité :
1. raideur cohésive divergente à endommagement → 0 (seuil d̃(Δt) sous
   lequel l'explicite est instable) ;
2. **LA DOMINANTE : la transition cohésif ↔ contact** — la discontinuité
   entre raideur cohésive k⁺ et pénalité de contact k⁻ à séparation nulle
   accumule des erreurs d'énergie à chaque passage répété, « même sous les
   conditions de stabilité standard », jusqu'à dérive énergétique et
   fragmentation artificielle. C'est mot pour mot notre pompe de phase
   débris (joints comprimés survivants qui alternent avec le contact) ;
3. l'adoucissement progressif est bénin (perte algorithmique compensée par
   la dissipation de rupture).

Remèdes qu'ils évaluent : Δt ≲ 0,2·Δt_c pendant les transitions ;
**pénalité adaptative k⁻(d) = k⁺(d)** (continuité de raideur, au prix
d'interpénétration) ; surveillance des « points chauds » ; critère de
Gershgorin (le CFL standard IGNORE les raideurs de contact) ; et leur
recommandation finale : schémas **non lisses par impulsions**
(arXiv:2606.01355, Newmark semi-explicite non lisse) — la pénalité pure
n'est « pas viable » pour une fragmentation conservative.

**Corollaire pour rockim** : notre contact général au POTENTIEL est déjà
conservatif (ΔKE/KE₀ ~ 1e-12) — le maillon faible restant est côté
**joints** : dashpot compressif + pénalité du joint rompu + la transition
joint↔contact.

**Pratiques des codes de production** :
- *Abaqus (doc Explicit)* : l'amortissement RÉDUIT le pas stable (un dt
  calculé sans amortissement est trop grand avec) ; les nœuds à petite
  masse pilotés par des éléments sans masse sont un danger identifié.
- *Yade (doc formulation)* : amortissement de Cundall borné par
  construction (non visqueux) ; facteur de sécurité systématique sur
  Δt_cr ; « le Δt requis peut être inférieur à Δt_cr » à hautes vitesses.
- *MOOSE* : borne dashpot cd ≤ m/dt (déjà adoptée).
- *OpenFDEM (github.com/OpenFDEM-geomechanics), Y-HFDEM, HOSS, Irazu* :
  mêmes fondements Munjiza ; la sensibilité aux paramètres de contrôle
  (pénalités, amortissement) est documentée comme un enjeu de premier
  ordre (Springer 2023, FDEM GPU).

## 3. Diagnostic proposé (hiérarchisé, falsifiable)

Hypothèse causale : **le cycle {joint comprimé survivant → force de
joint/dashpot → éjection du nœud → contact général → re-fermeture}** sur
des tets minces de Kuhn à petite masse nodale, où trois mécanismes peuvent
pomper :
- H1 (transition, favori depuis 2511.14323) : discontinuité de raideur
  joint↔contact à chaque alternance — gain net d'énergie par cycle ;
- H2 (dt local) : le dt global (ressorts de joints + CFL) ne voit PAS la
  raideur de contact effective locale sur un nœud quasi isolé — Gershgorin
  local violé pendant les impacts de débris ;
- H3 (dashpot résiduel) : la borne MOOSE utilise un m_eff qui peut rester
  trop optimiste pour un nœud d'éventail quasi isolé.
La borne MOOSE posée en 08-11/13 a peut-être déjà tué une partie de la
pompe : à VÉRIFIER avant toute nouvelle modification.

## 4. Plan d'attaque (chaque étape a un critère de réfutation)

- **E0 — Rejouer et instrumenter** (l'étape gelée, maintenant outillée) :
  le cas explosif intrinsèque T = 4e-4 avec le code ACTUEL (borne MOOSE +
  CFL réelle) et le budget B4 par famille + log des hotspots (nœud le plus
  rapide, famille de force dominante à l'instant du kick). Verdicts
  possibles : (a) plus d'explosion → la borne MOOSE était la pompe, on
  passe direct à la validation E3 ; (b) explosion → LE compteur B4 qui
  monte en positif désigne H1/H2/H3.
- **E1 — Correctif ciblé selon E0** (opt-in, bit-identique par défaut) :
  H1 → pénalité de joint rompu continue k⁻(d) = k⁺(d) (éq. 2511.14323) ;
  H2 → critère de Gershgorin LOCAL incluant les raideurs de contact
  actives → warning puis abort propre ; H3 → m_eff = masse nodale réelle
  min de la paire.
- **E2 — Garde-fous systémiques** (défense en profondeur, tous opt-in) :
  moniteur d'énergie runtime (résidu B4 > seuil → arrêt propre + dump =
  l'« energy sanity abort » des codes de production) ; plancher de
  masse/fusion des micro-fragments OU suppression comptabilisée au budget
  (pratique débris standard) ; option Δt réduit en présence de transitions
  fréquentes (compteur de passages joint↔contact).
- **E3 — Validation** : rebond complet de la percussion 3D homogène
  (le cas déclencheur) avec budget < 1 % et e physique, les DEUX schémas,
  + repère de suite dédié (« rebond 3D ») pour que la pathologie de phase
  tardive ne redevienne jamais invisible.

Règle de travail : chaque run de ce plan est soumis à validation de
Fernando (config exacte + coût) AVANT lancement.

## Sources

- arXiv:2511.14323 — Stability of Extrinsic CZM with Penalty-Based Contact
  in Explicit Dynamic Fragmentation (2025)
- arXiv:2606.01355 — Semi-explicit nonsmooth Newmark for robust unilateral
  contact in dynamic fragmentation (2026)
- Docs Abaqus Explicit §6.3.3 (stabilité, amortissement, petites masses)
- Yade — DEM Formulation (amortissement de Cundall, facteurs de sécurité)
- OpenFDEM (github.com/OpenFDEM-geomechanics), openfdem.com ; Y-HFDEM ;
  HOSS (LANL) ; Irazu (Geomechanica)
- Springer G&G 2023 — Influence of model control parameters on fracture
  characteristics of GPU parallel FDEM
- FICHE_rockim.md §instabilité latente (autopsie interne 2026-08-07)
