# Constitution rockim

*La discipline maison du solveur FDEM rockim, éprouvée sur les chantiers
A1-A3 / N1 / V1 / V2 (2026-08), formalisée en loi de projet. Tout plan,
toute tâche et toute implémentation issus de spec-kit se vérifient contre
ce document.*

## Core Principles

### I. Bit-neutralité des défauts (NON NÉGOCIABLE)
Toute capacité nouvelle est **opt-in** derrière une clé de configuration ;
toute configuration existante produit des trajectoires **bit-identiques**
(à 1 thread) avant/après le changement. Un changement qui altère les défauts
exige : (a) une décision explicite de Fernando, (b) le recalage documenté des
références de la suite, (c) la mention « dernier changement de ce type » si
la classe de changement est close (ex. ordre canonique des paires). Les
compteurs et l'instrumentation sont des lectures pures : aucun flottant de la
physique ne change.

### II. Le contrôle à charge nulle
Après toute modification touchant joints, contact ou maillage : un run à
charge nulle (`pullV = 1e-12` ou corps au repos) doit donner **0 joint cassé**
et un travail de contact **exactement nul**. C'est le test le plus
discriminant du dépôt (3 bugs attrapés sur le seul chantier A3) ; il est
non négociable et s'exécute AVANT toute mesure de performance.

### III. 2D et 3D de front
Toute capacité physique naît dans les deux solveurs (FdemSolver,
Fdem3dSolver) dans le même chantier. Exception admise : l'objet n'existe pas
dans l'une des dimensions (ex. multi-corps par groupes physiques, 3D-only) —
l'exception est alors documentée. L'instrumentation diagnostique temporaire
peut être mono-dimension si le problème diagnostiqué l'est.

### IV. La mesure fait foi
Aucune optimisation sans mesure avant/après, **séquentielle, machine idle,
OMP_NUM_THREADS = 1**. Les critères d'acceptation chiffrés sont posés AVANT
le chantier. Un diagnostic s'établit aux **compteurs** (issues par branche,
temps par section), pas à l'intuition — le chantier N1 a réfuté deux
diagnostics successifs par la mesure avant de trouver le poste dominant.
Un critère manqué se déclare NON VALIDÉ avec ses raisons mesurées ; il n'est
jamais maquillé en succès partiel.

### V. La suite de non-régression comme contrat
Chaque capacité laisse des repères dans `tools/verify_suite.py` (références
en dur, tolérances par nature de test, tiers fast/full/all). La suite fast
tourne après chaque étape de chantier ; full/all avant commit quand le
chantier touche leurs domaines. Les références sont PAR PLATEFORME
(libstdc++ vs MSVC divergent à graine égale) ; la certification est à
1 thread. Un repère qui casse arrête le chantier jusqu'à explication.

### VI. Le bilan d'énergie comme juge physique
`dampWork ≤ 0` partout (un amortisseur ne peut que dissiper) ; le budget
d'énergie par sous-système (V2/B4) doit boucler à **< 1 % du flux brut** sur
tout scénario d'impact ; tout canal de force nouveau est instrumenté
(travail cumulé) dès sa naissance. Une injection d'énergie inexpliquée est
un bug bloquant, pas un bruit.

### VII. Documentation au fil de l'eau
`DOCUMENTATION_rockim.md` (référence des clés et sorties — le code fait
foi), `ROADMAP_rockim.md` et le plan actif sont mis à jour **dans le même
commit** que le code. Les leçons chèrement acquises (pièges, règles maison)
sont consignées avec leur date et leur mesure. Les messages de commit
racontent le POURQUOI et les chiffres, pas seulement le quoi.

## Contraintes techniques

- **C++17, Eigen headers-only, OpenMP optionnel** : le code compile et
  tourne en sériel sans OpenMP ; MSVC et g++ sont supportés à parts égales.
  Toute dépendance nouvelle (bibliothèque, runtime, toolchain — CUDA/MPI
  compris) est une décision d'architecture explicite, jamais un effet de
  bord d'implémentation.
- **SI partout** (m, s, Pa, kg) ; en 2D les grandeurs sont par mètre
  d'épaisseur.
- **Reproductibilité** : graine par binaire ; à N threads fixé le résultat
  est déterministe ; les comparaisons fines se font à threads égaux.
- **Un exécutable unique** piloté par fichier de configuration texte
  (`clé = valeur`) ; pas de format de config nouveau sans décision.

## Workflow de développement

- Tout chantier d'effort moyen ou plus suit le flux spec-kit :
  `constitution → specify → clarify → plan → tasks → implement`, avec les
  artefacts sous `specs/`. Les chantiers courts (< 0,5 j) peuvent s'en
  dispenser mais respectent les principes I-VII.
- Les décisions **non bit-neutres** (changement de méthode numérique,
  d'architecture, de dépendances) sont présentées à Fernando avec options
  chiffrées AVANT implémentation.
- Commits par chantier propre sur la branche `article-exact`, bundle livré ;
  les runs de mesure chronométrés ne partagent jamais la machine avec un
  build ou un autre run.

## Governance

La constitution prime sur toute autre pratique. Amendement = commit dédié
avec version incrémentée (sémantique : MAJEUR = principe retiré/redéfini,
MINEUR = principe ou section ajouté, PATCH = clarification), date, et
justification. Les plans et tâches spec-kit citent la constitution dans
leurs « constitution checks » ; une violation détectée à la revue arrête
l'implémentation.

**Version**: 1.0.0 | **Ratified**: 2026-08-14 | **Last Amended**: 2026-08-14
