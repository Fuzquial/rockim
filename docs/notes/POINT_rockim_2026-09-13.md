# Le point sur rockim — 13 septembre 2026 (demande de Fernando, 11 h 15)

**Ce que rockim sait faire aujourd'hui.** Un FDEM 3D complet pour l'impact d'un insert sur du granite :
tétraèdres et joints cohésifs intrinsèques entre toutes les facettes, la loi de Solidity ingrédient par
ingrédient (enveloppe de Mohr-Coulomb avec cut-off, parabole élastique, adoucissement en z-curve, plages
3G/ft et 3G_II/f_s, quadrature à trois points, mort à deux points sur trois, DIF de Yang), endommagement de
volume pour la pulvérisation, contact par potentiel de Munjiza entre fragments avec frottement par
matériau, corps d'acier et de carbure continus, train de frappe complet à quatre corps, insertion adaptative
en option, un bilan d'énergie poste par poste avec une borne physique qui arrête tout run qui crée de
l'énergie, une bit-identité vérifiée à chaque changement. Le run s = 1 du 13/09 (120 185 tétraèdres) tourne
sur une loi qui ne crée plus d'énergie et fait à 160 µs ce que Yang montre à 168 µs : un disque broyé de
9 mm, des amorces radiales au bord, aucune radiale longue — leur cas à 9 m/s aussi.

**Ce qui ne manque pas.** Les équations : écrites, comparées à leur code source, les deux erreurs de
transcription corrigées (sécante suivant la pression, σ_n divisée par deux). Les garde-fous : résidu,
borne physique, ancres. Les outils de lecture : sept critères, coupes exactes, force–pénétration,
cinétique. Le maillage : le générateur reproduit leur géométrie (rendu 37 % plus grossier au cœur, réglable
par l'échelle). Les données d'Aising en base.

**Ce qui manque, par ordre d'importance.**
1. La tenue du lit broyé sous l'insert : 20–45 kN de freinage contre ~57 kN impliqués par Yang. Loi
   (joint dissipatif et à glissement résiduel contre joint réversible qui guérit ; mort après 4 µm sous
   400 MPa) ou paramètres non publiés (pénalité de joint, de contact). Deux bancs de 2 h tranchent.
2. La pulvérisation : 0 élément contre ~360 ; chez eux le volume s'écrase parce que le joint tient. Et
   δ_m = h·ε n'est pas objective au maillage.
3. L'énergie de l'essai : piston maillé 42,8 J à 9 m/s, leur bit en reçoit 23,84 ; piston 36 % trop lourd.
4. Rebond et fragments : 600–750 µs et la brosse de tri ; le run s'arrête à 400.
5. Un maillage à leur taille (0,7 / 1,0 mm), qui triple le coût.

**Pourquoi nos runs sont plus longs.** Eux : 230 788 éléments, dt 2,5 ns, 5 h. Nous : 120 185, dt 1,0 ns,
15 h. Le pas de temps : slivers de gmsh (1,0 ns contre 4,7 ns pour les éléments seuls), 98 % de la raideur
qui commande le pas vient des ressorts de pénalité (pf = 20 E/h). Le coût par pas : contact et joints
étaient 80 % du pas sans profiter des fils ; contact parallélisé le 11/09, fusion des forces de joint le
13/09 (×6 sur ce poste) ; le régime des débris multiplie encore le pas par trois. Reste ~×5 par élément et
par pas, moitié dt, moitié implémentation (audit C : CSR des groupes, cache du contact, AVX2).

**Pourquoi pas de radiales longues.** Chez eux non plus à 9 m/s : radiales de 9–10 mm depuis le centre
pour un cratère de 7 mm, nées en charge, visibles à 168 µs ; latérales et copeaux seulement à 11 m/s. Nous
à 160 µs : amorces de 2–3 mm au bord d'un disque de 8,7 mm, la même image. Ce qui les borne chez nous :
le lit broyé trop mou pousse moins en coin (moins de traction circonférentielle), la maille de 1,4 mm (six
facettes en zigzag pour 10 mm), et à 2,5 mm elles n'existaient pas.

**Corrigé cette semaine, définitivement.** Deck v1 (onze écarts), loi `origin` (moitié de l'énergie du
bit créée en 20 µs, cuvette de 33 mm prise pour un cône), résidu aveugle, runs 2D sur amortissements,
σ_n/2 sous parabole, acier et carbure endommageables, contact roche/roche à la raideur du carbure, fusion
série des joints.

**La question ouverte.** Notre joint est conforme à la thermodynamique, le leur pas tout à fait ; la seule
différence de physique qui reste est ce que porte un joint comprimé avant de mourir — mesurable en deux
bancs de deux heures (`origin` + coulomb + ratchet ; leur loi mot pour mot sous la borne).
