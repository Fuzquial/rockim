#pragma once
// ---------------------------------------------------------------------------
// Fem3dSolver — 3D explicit continuum FEM on Kuhn tetrahedra with PLUGGABLE
// constitutive laws (MatLaw: elastic | dpr | saksala). The 3D sibling of the
// 2D fem module and the sandbox counterpart of an Abaqus/Explicit + VUMAT
// percussion run: shared-node mesh (no cohesive joints — fracture is smeared
// damage + element EROSION), co-rotational linear tets, rigid spherical
// tool, Lysmer quiet boundaries, percussion and tension/compression
// scenarios (tension with pullV < 0 IS the uniaxial compression test, with
// pullRamp and gripLateralFree available as in the 2D FDEM).
//
// Verification hooks: the DP cone's uniaxial compressive strength and the
// Perzyna viscous overstress have closed forms (MatLaw::sigmaCdp, README) —
// the tension/compression scenario prints measured vs expected.
// ---------------------------------------------------------------------------
#include <array>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "rockim/Config.hpp"
#include "rockim/MatLaw.hpp"
#include "rockim/Material.hpp"
#include "rockim/Solver.hpp"

namespace rockim {

class Fem3dSolver : public Solver {
public:
    Fem3dSolver(const Config& cfg, std::string outDir);

    void init() override;
    void step() override;
    void writeFrame(int frame) override;
    void historyHeader(std::ostream&) const override;
    void historyRow(std::ostream&) const override;
    void finalize() override;

private:
    enum Flag { FREE = 0, FIXED = 1, PRESCRIBED = 2 };
    enum class Scenario { PERCUSSION, SHEAR, TENSION };

    struct Elem {
        std::array<int, 4> n;                  // SHARED node ids
        Eigen::Matrix<double, 3, 4> dN;
        double V0 = 0.0, lc = 0.0;             // volume, size V0^(1/3)
        double svm = 0.0, pm = 0.0, szz = 0.0; // outputs
        // etude briques (2026-09-03) : contrainte laterale moyenne (jauge de
        // confinement), det F (soupape geometrique), energie elastique
        // nominale (energie effacee a l'erosion), canal d'erosion
        double slat = 0.0, J = 1.0, wEl = 0.0;
        char erodedBy = 0;                     // 0 vivant, 1 loi, 2 det F
        // essai triaxial continu (triaxStats, 2026-09-04) : composantes
        // normales et deformations de Biot pour la jauge du tiers central
        double sxx = 0.0, syy = 0.0, ezz = 0.0, ev = 0.0;
        // viscosite de volume (bulkViscosity, 2026-09-05) : J = det F du pas
        // PRECEDENT, memoire du taux volumique (J_n+1 - J_n)/(dt J_n+1) ;
        // lue et mise a jour seulement quand la cle est active
        double Jbv = 1.0;
        // ---- materiau PAR PHASE (2026-09-06, opt-in, defaut bit-identique) --
        // phase : indice dans PhaseSet (0 = la fiche materiau globale quand la
        // cle `phases` est absente, donc chemin historique inchange) ; grain :
        // cellule de Voronoi (mesh = voronoi) ou groupe physique (mesh = file).
        // Miroir litteral de Fdem3dSolver.hpp l. 100-101. Ces deux indices sont
        // poses DANS la boucle de construction des elements (jamais dans une
        // seconde passe indexee) : finishMesh conserve l'ordre des tets
        // d'entree mais rien ne le garantirait apres un futur tri.
        int phase = 0, grain = 0;
        MatState st;
    };
    struct BFace { int elem; std::array<int, 3> n; };
    struct Tool3 {
        double mass = 0.5, radius = 0.015;
        bool free = true;
        bool flat = false;   // flat-ended cylindrical punch (axis z); x is
                             // then the CENTER OF THE BOTTOM FACE, radius
                             // the punch radius — the 3D lift of the 2D
                             // FLAT tool (percussion only)
        // toolShape = blade (shear, 2026-09-03) : lame rigide plane sur toute
        // l'epaisseur, x = POINTE (arete), face de coupe a `rake` de la
        // verticale (vers l'arriere), face de depouille a `clear` de
        // l'horizontale, hauteur de face `hFace`
        bool blade = false;
        double rake = 0.0, clear_ = 0.0, hFace = 0.02;
        // toolLockXY (2026-09-05, w17, opt-in, defaut false = bit-identique) :
        // l'outil rigide est bloque en x, y — v_x = v_y = 0 apres la mise a
        // jour par F, donc x et y restent EXACTEMENT a leur valeur initiale
        // (x += dt * 0). F.x, F.y sont toujours calcules par toolContact() et
        // comptes (toolFx de l'historique, peakF_) : c'est la reaction du
        // noeud de reference bloque (*Boundary RPB1 1,2 des decks Abaqus).
        bool lockXY = false;
        Eigen::Vector3d x{0, 0, 0}, v{0, 0, 0}, F{0, 0, 0};
        void integrate(double dt) {
            if (free) v += (dt / mass) * F;
            if (lockXY) { v.x() = 0.0; v.y() = 0.0; }
            x += dt * v;
        }
        double ke() const { return 0.5 * mass * v.squaredNorm(); }
    };

    void buildMesh();
    void buildMeshFile();                        // mesh = file (Gmsh MSH 2.2)
    void buildMeshVoronoi();                     // mesh = voronoi (grains + phases)
    void finishMesh(const std::vector<std::array<int, 4>>& tets);
    void checkMassAudit(const char* tag);        // sum(m) == sum_p rho_p V_p
    void phaseKeyGuards();                       // cles de JOINT refusees en fem3d
    void checkFinite();                          // C4 (w20) : NaN/Inf reel
    void placeTool();
    void setupBoundaries();
    void setupConfinement();
    void confiningForces();
    double achievedConfinement() const;
    void refreshActiveNodes();
    void computeStableDt();
    void elementForces();
    void toolContact();
    void integrate();
    double craterVol() const;
    void setupProbes();                          // probes = x,y,z ; ... (opt-in)
    void probesHeader(std::ostream&) const;
    void probesRow(std::ostream&) const;

    Config cfg_;
    std::string out_;
    Material mat_;
    // ---- materiau par phase (2026-09-06) -----------------------------------
    // phases_ contient TOUJOURS au moins une fiche : sans la cle `phases`,
    // PhaseSet::from rend {Material::from(cfg)} = mat_ champ par champ, donc
    // phases_.mat[0] est le MEME objet que mat_ et toute substitution
    // mat_.X -> phases_.mat[0].X rend le meme double, bit pour bit. C'est le
    // pivot de la bit-identite : un SEUL chemin de code, jamais deux.
    PhaseSet phases_;
    // Une INSTANCE DE LOI PAR PHASE, meme cle `law` pour toutes. Impose par le
    // code, pas par gout : MatLaw fige lam_, G_, K_, adp_, kdp_ dans son
    // CONSTRUCTEUR a partir de la fiche recue (MatLaw.hpp) et stress() ne
    // recoit PAS de materiau — un objet unique ne peut pas porter trois
    // modules d'Young, et un simple champ `phase` sur l'element ne changerait
    // STRICTEMENT RIEN a l'elasticite. Les options de loi (erodeD, dfh*, cdp*,
    // meridian...) sont lues dans le Config GLOBAL par MatLaw::make : elles
    // restent donc communes a toutes les phases (documente, refuse si ecrit
    // par phase).
    std::vector<std::unique_ptr<MatLaw>> laws_;
    // pointeur NON PROPRIETAIRE sur laws_[0] : les sites qui n'appellent que
    // name() / sigmaCdp() / viscousOverstress() restent ecrits comme avant
    // (diff minimal). stress() est const et tout l'etat mutable vit dans
    // MatState : un vecteur de lois partage en lecture par N fils est sur.
    const MatLaw* law_ = nullptr;
    // tables CHAUDES par phase (modele Fdem3dSolver.cpp l. 410-418) : rho pour
    // la masse condensee et la viscosite de volume, rho c_d pour sa part
    // lineaire. cP() contient un sqrt : on ne l'appelle pas par element et par
    // pas.
    std::vector<double> rhoP_, rhoCdP_;
    // phase / grain par tet, remplis par le front-end de maillage AVANT
    // finishMesh. VIDES = comportement historique (tout en phase 0).
    std::vector<int> tetPhase_, tetGrain_;
    std::vector<std::string> groupName_;         // mesh = file : volumes physiques
    int nGrains_ = 0;
    bool meshVoronoi_ = false;                   // chemin mesh = voronoi actif
    // La cle `phases` est-elle POSEE dans le deck ? A distinguer de
    // phases_.n() > 1 (correction du 2026-09-06) : avec `phases = dur`, une
    // seule fiche est declaree et phases_.n() vaut 1 — toutes les gardes
    // conditionnees par n() > 1 etaient alors muettes, alors que le deck
    // AFFIRME decrire une microstructure. C'est l'intention du deck qui doit
    // armer les gardes, pas le nombre de fiches qu'il a fini par produire.
    bool phasesDeclared_ = false;
    // phaseWeibull (opt-in) : autorise explicitement la COMPOSITION du champ
    // de Weibull (matWeibullM, qui ecrit ftScale par element) avec les phases.
    // Refusee par defaut : les deux heterogeneites se multiplient et un run
    // qui localise ne serait plus attribuable — or le constat qui motive tout
    // ce chantier est justement qu'un Weibull sur une raideur UNIFORME ne
    // pouvait pas localiser.
    bool phaseWeibull_ = false;
    Scenario scen_ = Scenario::PERCUSSION;

    double W_ = 0.1, D_ = 0.1, H_ = 0.08;
    int nx_ = 24, ny_ = 24, nz_ = 18;
    double hmin_ = 1e-3, lcMax_ = 1e-3;

    std::vector<Eigen::Vector3d> X0_, u_, v_, f_;
    std::vector<double> m_;
    std::vector<int> flag_;
    std::vector<Elem> el_;
    std::vector<BFace> exterior_;
    std::vector<Eigen::Vector3d> cAbs_, kAbs_;

    Tool3 tool_;
    double toolKE0_ = 0.0;
    // suivi lateral de l'outil (2026-09-05, w17 ; console seulement, resume) :
    // position initiale, max |x - x0|, |y - y0| et max |F_x|, |F_y| sur tous
    // les pas — la derive laterale du quart de bloc (contre-expertise
    // perc3d/DECHARGE_abaqus_vs_rockim.md § 6 : -0,18 mm a 320 us) etait
    // invisible hors de history.csv (toolX) et frames.csv (toolY par image).
    Eigen::Vector3d toolX0_{0, 0, 0};
    double toolDxMax_ = 0.0, toolDyMax_ = 0.0, toolFxMax_ = 0.0, toolFyMax_ = 0.0;
    double kp_ = 0.0, muC_ = 0.5, xiC_ = 0.05, vReg_ = 1e-3;
    // contactPenaltyFactor (2026-09-04, opt-in, 1 = bit-identique) : kp =
    // facteur x E hmin. Sur un Delaunay hmin est un sliver, E hmin est alors
    // 4-5 fois plus souple que la raideur nodale de la roche (interpenetration
    // de l'ordre de l'indentation au pole de l'insert). dt reste borne par
    // 2 sqrt(m_min/kp) (computeStableDt) ; omega_p/(2/dt) imprime, avertissement
    // au-dela de 0,5 du rapport COMBINE sqrt((dt/CFL)^2 + (omega_p dt/2)^2)
    // (revue : la raideur d'element du noeud de contact compte aussi).
    // kpFactor (cle du mode fem 2D) accepte comme alias en fem3d (revue).
    double kpFac_ = 1.0;
    // vtkCap (revue du cap cdp, opt-in, false = VTU bit-identiques) : ajoute
    // aux VTU les champs cellulaires capPc (MatState::pc, Pa) et epsVpl
    // (-tr(eps_pl), compaction > 0) — la zone broyee par le cap (cdp ou dpr)
    // n'etait visible dans aucun champ.
    bool vtkCap_ = false;
    // tensionDamage = fixed (2026-09-05, opt-in) : champs VTU dFix1..3 et
    // colonnes probes dFix1..3 (validite de la cle : MatLaw::make)
    bool fixedDam_ = false;
    // toolContact = signorini (2026-09-04, port de la branche A1 de fdem3d,
    // chantier du 2026-09-02 ; defaut penalty = bit-identique) : contact
    // outil en IMPULSION (CD-Lagrange, ToolSignorini.hpp) — par noeud, avec
    // H = 1/m_i, v* = v + (dt/m) f, condition de Signorini sur le gap predit,
    // impulsion normale r_n = m (relax pen/dt - v_n), cap de Coulomb sur
    // l'impulsion tangentielle ; report en force r/dt pour integrate(). Sans
    // raideur ni interpenetration cumulative : le banc de Hertz est le juge.
    // Sphere et poincon plat ; la lame garde la penalite (geometrie a deux
    // faces, non portee).
    bool toolSig_ = false;
    double toolSigRelax_ = 0.0;
    // toolPulseForce / toolPulseTable (2026-09-05, opt-in, defaut absent =
    // bit-identique) : l'insert est POUSSE par une impulsion de force imposee,
    // comme le bouton des decks Abaqus de la these (tige + *Cload amplitude
    // Pulse) — F_z(t) = -toolPulseForce x a(t), a(t) lineaire par morceaux
    // dans la table "t:a t:a ..." (0 hors de la table), appliquee a la masse
    // de l'outil APRES la comptabilite du contact (toolFz de l'historique et
    // peakF_ restent la force de CONTACT seule). Avec impactSpeed = 0 et
    // toolGap petit : le meme protocole que P_cdpQ_v11. Compteurs :
    // pulseImp_ = int F dt, pulseWork_ = int F v dt (resume).
    double pulseF_ = 0.0, pulseWork_ = 0.0, pulseImp_ = 0.0;
    std::vector<double> pulseT_, pulseA_;
    // toolPulseImpedance (N s/m, defaut 0 = force pure) : SOURCE A IMPEDANCE
    // — la tige est une ligne caracteristique d'impedance Z = rho c A ; l'onde
    // incidente F_inc(t) = toolPulseForce x a(t) donne au bouton la force
    // F = 2 F_inc − Z v_bouton (bouton bloque : 2 F_inc ; bouton libre :
    // v = 2 F_inc / Z, F = 0). Sans Z, une force imposee pousse l'insert a
    // travers la roche (verifie : 8 mm et 13 m/s sur la carte historique).
    // Dans les decks Abaqus a *Cload au sommet d'une tige absorbee par des
    // dashpots Z_d : F_inc = Cload x Z/(Z + Z_d).
    double pulseZ_ = 0.0;
    double damping_ = 0.05, pullV_ = 0.05, pullRamp_ = 0.0;
    bool gripFree_ = false;

    long stepCount_ = 0, nEroded_ = 0;
    long nanEvery_ = 256;                        // nanCheckEvery (C4, w20)
    std::vector<long> gmshNodeId_;               // mesh = file : id Gmsh par noeud (C3)
    double work_ = 0.0, peakF_ = 0.0, sigmaPeak_ = 0.0;

    // ---- etude briques (2026-09-03), tout opt-in, defaut bit-identique ----
    // confinement suiveur (copie de Fdem3dSolver) + pression de dessus (boue)
    double confP_ = 0.0, confRamp_ = 30e-6, topP_ = 0.0, bottomP_ = 0.0;
    std::vector<int> confFaces_, topFaces_, bottomFaces_;   // indices dans exterior_
    double confWork_ = 0.0, confAchieved_ = 0.0, confGaugeT_ = 0.0;
    bool confLatched_ = false;
    double toolDelay_ = 0.0;                     // outil gele avant t
    bool symY_ = false;                          // tranche plane : v_y = 0
    // quarterModel (2026-09-04, quart de bloc pour la percussion confinee) :
    // les plans x = 0 et y = 0 sont des plans de SYMETRIE (v_x = 0 sur x = 0,
    // v_y = 0 sur y = 0, rien ailleurs — ce n'est pas la deformation plane de
    // symmetryY) ; ni pression suiveuse ni Lysmer sur ces deux faces ; les
    // faces x = W et y = D restent des bords ordinaires. Insert sur l'arete
    // (0, 0), masse d'outil a diviser par 4 dans le deck.
    bool symQ_ = false;
    // bottomFree (2026-09-05, opt-in) : avec absorbing = none | sides, le fond
    // n'est plus encastre — bloc flottant tenu par son inertie, comme le bloc
    // des decks Abaqus de percussion sans condition au fond. Defaut false =
    // bit-identique (fond FIXED).
    bool bottomFree_ = false;
    // mecanisme de suppression : masque de noeuds actifs + soupape det F
    bool activeNodes_ = false, activeDirty_ = false;
    std::vector<char> active_;
    double erodeDetMin_ = 0.0, erodeStrainMax_ = 0.0;
    long nErodedLaw_ = 0, nErodedGeo_ = 0;
    double vErodedLaw_ = 0.0, vErodedGeo_ = 0.0, eRemoved_ = 0.0;
    // observables de champ (fieldStats) : colonnes ajoutees en fin de ligne
    bool stats_ = false;
    Eigen::Vector3d gripF_{0, 0, 0};
    // mid-third gauge: the grip force carries the Cundall-damping drag of
    // the flowing column (measured +11 % at damping 0.7), the mid-specimen
    // stress does not — verifications read here
    std::vector<int> midEl_;
    double sigMid_ = 0.0, sigMidPeak_ = 0.0;
    double sigMidSum_ = 0.0;               // plateau average (last quarter):
    long sigMidN_ = 0;                     // the end-state snapshot of a
                                           // ringing signal is mesh-fragile

    // ---- essai triaxial continu (2026-09-04), opt-in, defaut bit-identique --
    // pullDelay : avant t = pullDelay les noeuds PRESCRIBED du mors sont
    // LIBRES (charges par topPressure) — consolidation isotrope exacte, puis
    // phase deviatoire avec la rampe pullRamp comptee depuis pullDelay ;
    // gripSection : section reelle du mors (defaut W D) ; triaxStats :
    // colonnes sigZZmid, sigXXmid, sigYYmid, epsAxMid, epsVolMid, epsAxGrip
    double pullDelay_ = 0.0, gripSection_ = 0.0;
    bool triax_ = false;

    // ---- viscosite de volume a la Abaqus/Explicit (2026-09-05), opt-in ----
    // bulkViscosity = b1 b2 (defaut absent, ou "0 0" = inactif : bit-identique).
    // Par element, avec edot = taux de deformation volumique (J_n+1 - J_n)/(dt J_n+1),
    // L_e = lc = V0^(1/3), c_d = cP() non endommage :
    //   p1 = b1 rho c_d L_e edot            (lineaire, tout signe)
    //   p2 = rho (b2 L_e edot)^2            (quadratique, compression edot < 0 seulement)
    // contrainte visqueuse sig_bv = q I, q = p1 - [edot < 0] p2 (q du signe de edot :
    // en compression rapide elle AJOUTE de la compression), ajoutee a la contrainte
    // de la loi pour les FORCES INTERNES seulement (svm/pm/szz/VTU = loi, comme
    // Abaqus). Dissipation wBulk = sum V0 q edot dt >= 0 (colonne de history.csv
    // et resume, seulement si actif). dt : facteur sqrt(1 + b1^2) - b1 sur la CFL.
    bool bv_ = false;
    double bvB1_ = 0.0, bvB2_ = 0.0, wBulk_ = 0.0;

    // ---- kinematics = biot | hencky (2026-09-05, w16), opt-in ---------------
    // biot (defaut, cle absente) : chemin historique bit-identique — F = R U par
    // iteration, eps = sym(R^T F) - I (Biot), forces f_a = -V0 (R sigma_c) dN_a/dX
    // (couple conjugue Biot/Biot, reference). hencky : decomposition spectrale de
    // C = F^T F (Fem3dKinematics.hpp), R = F U^-1 exacte, eps = ln U (deformation
    // logarithmique TOTALE : = ln U d'Abaqus/Explicit pour un trajet coaxial, le
    // taux integre d'Abaqus n'est pas porte), la loi rend sigma_c = Cauchy
    // co-rotationnelle, P = J R sigma_c U^-1 (forces sur la configuration
    // courante). Inchanges : masse lumpee sur V0, Lysmer, contact, lc = V0^(1/3),
    // CFL sur la geometrie initiale, soupapes (erodeStrainMax lit |eps| = |ln U|
    // en hencky), sorties svm/pm/szz/slat/jauges = contrainte de la loi (Cauchy
    // en hencky), wEl = 0.5 sigma_c : (eps - epsP) par unite de volume de
    // reference (nominal). Viscosite de volume : q ajoute a sigma_c avant P dans
    // les deux cinematiques.
    bool hencky_ = false;

    // ---- sondes de point materiel `probes = x,y,z ; x,y,z ; ...` (2026-09-05,
    // w18, opt-in ; cle absente = bit-identique, AUCUN test par element) -------
    // Metres, repere du bloc recale [0,W]x[0,D]x[0,H] (dessus z = H ; mesh =
    // file est translate a l'origine a la lecture). A l'init : tetraedre qui
    // contient le point (coordonnees barycentriques par les gradients dN deja
    // stockes, tolerance 1e-9) — sinon le plus proche par centroide, avec
    // WARNING (point hors du bloc ou dans un trou). Plusieurs sondes peuvent
    // partager un element (un CRENEAU par element distinct, probeSlot_[e] >= 0
    // seulement pour les elements sondes : dans elementForces chaque element
    // est traite par un seul thread, aucune course). A chaque ligne
    // d'historique (meme cadence que history.csv) : une ligne de
    // <outputDir>/probes.csv — t, puis par sonde l'element, la contrainte de la
    // LOI (sigma_c : co-rotationnelle, avant viscosite de volume, celle des
    // forces), la deformation passee a la loi (Biot, ou ln U en hencky), p =
    // -tr/3 (compression > 0), q = sqrt(3/2 s:s), les trois principales
    // decroissantes, tr(eps_pl) et eps_pl equivalente sqrt(2/3 dev:dev), epvEq
    // de la loi, d_t = st.D, d_c = st.Dc, pc, det F, drapeau eroded, puis si cdp
    // epsTpl/epsCpl/d, si dpdfh Dv1..3. Element erode : sigma = 0, eps = derniere
    // valeur calculee (soupape det F : celle du pas precedent), drapeau 1.
    struct Probe {
        Eigen::Vector3d x{0, 0, 0};              // point demande (m)
        int elem = -1;                           // tetraedre retenu
        int slot = -1;                           // creneau dans probeSig_/probeEps_
        bool inside = true;                      // false = plus proche par centroide
        double dist = 0.0;                       // distance point - centroide (m)
    };
    std::vector<Probe> probes_;
    std::vector<int> probeSlot_;                 // par element : -1 ou creneau
    std::vector<Eigen::Matrix3d> probeSig_, probeEps_;   // par creneau
    std::unique_ptr<std::ofstream> probesOut_;   // <outputDir>/probes.csv

    // OpenMP scratch (shared-node scatter)
    std::vector<std::vector<Eigen::Vector3d>> fTL_;
    std::vector<std::vector<char>> seenTL_;
    std::vector<std::vector<int>> touchedTL_;
};

} // namespace rockim
