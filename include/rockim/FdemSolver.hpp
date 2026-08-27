#pragma once
// ---------------------------------------------------------------------------
// FdemSolver — 2D combined finite-discrete element method (FDEM) in the
// architecture of Munjiza (Y2D lineage): the v2 mode sketched in the README.
//
//   * Every CST owns its three nodes (node duplication). Elements are linear
//     elastic, CO-ROTATIONAL (2D polar decomposition), so detached fragments
//     can tumble without spurious straining.
//   * 4-node cohesive joint elements are inserted on EVERY interior edge from
//     the start (intrinsic approach). Glued behaviour through a penalty
//     p = jointPenaltyFactor * E / h; fracture through a damage-plasticity
//     cohesive law: mode I softens from ft over 2 Gf/ft of opening, mode II
//     from c over 2 Gf,II/c of frictional slip, with Coulomb friction
//     tan(phi) * (-sigma_n) active throughout compression and surviving as
//     the residual strength of fully broken joints. Damage is the shared
//     scalar max of the two drivers (simple mixed-mode coupling).
//   * Once a joint is fully broken and its faces have moved apart, it goes
//     DEAD and the freed faces enter a general node-to-edge penalty contact
//     (cell grid over active edges) so debris interacts with the crater and
//     with other debris. The rigid tool (disc, free or prescribed) presses on
//     any node, as in the FEM module.
//   * Lysmer + Deeks-Randolph viscous-spring quiet boundaries on the
//     exterior, fragments as connected components of live joints, same
//     outputs and history format as the other modes.
//
// Two lab-scale test capabilities sit on top of that core:
//
//   * scenario = brazilian — the INDIRECT TENSION test. geometry = disc cuts
//     a disc of diameter min(W, H) out of the mesh (grid or voronoi alike:
//     triangles whose centroid falls outside are dropped, so the rim is
//     ragged at ELEMENT scale, not at grain scale), and the disc is
//     compressed between two rigid flat platens — top prescribed at -pullV
//     with the pullRamp cosine rise, bottom fixed — in the same penalty +
//     dashpot + Coulomb contact as the flat tool. Platens rather than pinned
//     node rows on purpose: the contact arc must GROW with the load, which is
//     what keeps the centre (not the loading points) the most tensile place.
//     The reported strength is the ISRM formula sigma_t = 2 P / (pi D t).
//
//   * confiningPressure — a follower pressure on the exterior faces (current
//     configuration, lumped L/2 per node), ramped by confiningRamp. Applied
//     to the ORIGINAL exterior only: fluid presses on the specimen surface,
//     not inside freshly opened cracks (pressurizing crack faces would drive
//     them apart and is a different physical problem). Combined with
//     scenario = tension and pullV < 0 this is the TRIAXIAL test; combined
//     with percussion it is the confined-impact case. The achieved mean
//     lateral stress is measured and reported against the target at the end.
//
// Deviations from Munjiza's book, stated plainly: linear cohesive softening
// instead of his z(D) curve, 2-point (node-pair) joint integration,
// node-edge penalty contact instead of distributed potential contact, and a
// structured cross-diagonal mesh (with optional jitter) instead of an
// unstructured one — crack paths follow mesh edges.
// ---------------------------------------------------------------------------
#include <array>
#include <cstdint>
#include <unordered_map>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "rockim/Config.hpp"
#include "rockim/MatLaw.hpp"
#include "rockim/Material.hpp"
#include "rockim/Solver.hpp"
#include "rockim/Tool.hpp"
#include "rockim/YanSoftening.hpp"

namespace rockim {

class FdemSolver : public Solver {
public:
    FdemSolver(const Config& cfg, std::string outDir);

    void init() override;
    void step() override;
    void writeFrame(int frame) override;
    void historyHeader(std::ostream&) const override;
    void historyRow(std::ostream&) const override;
    // Essai 0 : compteurs d etat des interfaces pour
    // l historique (nInserted, nDamaging). Voir .cpp.
    void countInserted(long& nIns, long& nDam) const;
    // Early stop of the brazilian test once the post-peak load drop has been
    // seen (brazilianStopAfterPeak). Off by default: the run length stays T.
    bool finished() const override;
    void finalize() override;

private:
    // ROLLERX: lateral roller — the x displacement is held at zero, the y
    // dof stays free. Needed by the confined-strip verification of Yan et
    // al. (their section 3.1), whose flanks are "fixed in the normal
    // direction" only. Set by `lateralRollers = true`.
    // DRIVEX: SHPB velocity boundary — v_x is prescribed by the pulse
    // history shpbVel(t), the y dof stays FREE (a bar end is not laterally
    // clamped; clamping it would radiate a spurious shear wave).
    enum Flag { FREE = 0, FIXED = 1, PRESCRIBED = 2, ROLLERX = 3,
                DRIVEX = 4 };
    enum class Scenario { PERCUSSION, SHEAR, TENSION, BRAZILIAN, SHPB };

    struct Elem {
        std::array<int, 3> n;              // its own three nodes
        Eigen::Matrix<double, 2, 3> dN;    // reference shape-fn gradients
        double A0;                         // reference area
        double svm = 0.0;                  // von Mises (output)
        double exx = 0.0;                  // co-rotated axial strain
                                           // (output; SHPB gauges)
        double sxx = 0.0, syy = 0.0;       // global stress (gauges: confinement,
        double sxy = 0.0;                  // brazilian centre; full in-plane
                                           // tensor for the insertion criterion
        // Taux de deformation EQUIVALENT lisse [1/s] : max des |valeurs
        // propres| du tenseur taux de deformation D = sym(Fdot F^-1), passe
        // dans un filtre exponentiel de constante strainRateTau. Sert au DIF
        // de Yang et al. 2025. La mesure est la principale MAXIMALE et non un
        // equivalent de von Mises parce que les courbes DIF de la litterature
        // sont mesurees en essai UNIAXIAL : en uniaxial la principale max EST
        // le taux axial, sans facteur de convention. Vaut 0 tant que ni
        // bulkViscosity ni strainRateDIF ne sont armes (branche non calculee).
        double edot = 0.0;
        // WP1 pulverisation (bulkDamage = yang) : endommagement D et max
        // historique de delta_m = h_e * eps_vm. 0 si la cle est absente.
        double bdD = 0.0;
        double bdDm = 0.0;
        int phase = 0;                     // mineral phase (index in phases_)
        int grain = 0;                     // grain id (voronoi) / 0 (grid)
        MatState st;                       // only used when law_ is set
    };

    // 4-node cohesive joint between elements eA and eB. Node pairs
    // (a1, b1) and (a2, b2) are co-located initially; (a1 -> a2) is the CCW
    // traversal of the edge in element A, so the outward normal of A is
    // n = normalize( (Q - P).y, -(Q - P).x ).
    // Cohesive properties live PER JOINT: intra-grain joints carry the phase
    // material, grain-boundary joints the attenuated mean of the two phases
    // (type 0 = intra-grain, 1 = homophase boundary, 2 = heterophase).
    struct Joint {
        int eA, eB;
        int a1, a2, b1, b2;
        double L0;
        double D = 0.0;                    // shared cohesive damage [0..1]
        double slip[2] = {0.0, 0.0};       // jointShearUnload = plastic : le
                                           // glissement plastique cumule par
                                           // le retour radial.
                                           // = origin : l'ORIGINE FIGEE de la
                                           // secante de l'eq. 18, stampee par
                                           // activateJoint() (0 en intrinseque),
                                           // jamais mise a jour ensuite.
        double omax[2] = {0.0, 0.0};       // largest opening ever reached
                                           // (jointSoftening = yan: the omax
                                           // of eq. 17, per integration point)
        double smax[2] = {0.0, 0.0};       // largest SLIDING ever reached, the
                                           // s_max of eq. 18 (jointShearUnload
                                           // = origin only; unused otherwise)
        // Force NORMALE nette que le joint transmettait a l instant EXACT de
        // sa mort, en N par metre d epaisseur (negatif = compression). Sortie
        // de mesure seulement : c est la charge que le relais au contact doit
        // reprendre. Voir jointDeath dans ce header.
        double fDeath = 0.0;
        bool dead = false;                 // faces released to general contact
        double tBreak = -1.0;              // first time D reached 1
        int type = 0;
        double pj = 0.0;                   // penalty per area [Pa/m]
        double ft = 0.0, coh = 0.0;        // strengths [Pa]
        double Gf = 0.0, GfII = 0.0;       // fracture energies [J/m^2]
        double dnE = 0.0, dnF = 0.0;       // mode I elastic / final opening
        double slipF = 0.0;                // mode II softening slip
        double tanPhi = 0.0;               // friction
        double stat = 1.0;                 // Weibull strength factor (output)
        // ---- DIF de Yang et al. 2025, GELE a l instant de l insertion ----
        // Sortie seulement : les facteurs ont deja ete appliques a ft, coh,
        // Gf et GfII lors de activateJoint(). edotIns est le taux de
        // deformation moyen des deux elements a cet instant.
        double difT = 1.0, difC = 1.0, edotIns = 0.0;
        // En schema INTRINSEQUE (difIntrinsic_), le meme gel a lieu au
        // franchissement de l enveloppe et non a l insertion : ce drapeau
        // garantit qu il n a lieu QU UNE FOIS. Inerte en adaptatif.
        bool difStamped = false;
        // ---- mode de rupture (Yan et al. fig. 18 / 20) --------------------
        // Renseigne UNE FOIS, a l'instant ou D atteint 1, en comparant les
        // deux moteurs de l'eq. 16 : rn = (dn - dnE)/(dnF - dnE) (ouverture)
        // et rs = |slip|/slipF (glissement). bmode = 1 si rn >= rs (fissure
        // de TRACTION), 2 sinon (fissure de CISAILLEMENT). 0 = intact.
        // Sortie seulement : aucune force n'en depend.
        int bmode = 0;
        double rnB = 0.0, rsB = 0.0;       // les deux moteurs a la rupture
        // ---- adaptive insertion (Yan, Zheng & Wang, IJRMMS 169, 2023) ----
        // bonded = the joint does not exist yet: the co-located node copies
        // are RIGIDLY BOUND (exact shared-node FEM, zero artificial
        // compliance) and processJoint skips it. It is activated — "inserted"
        // — the first time the edge-averaged traction reaches the strength
        // envelope. dn0 is the opening offset stamped at activation so the
        // newborn joint reproduces the transmitted normal traction at zero
        // geometric opening (stress continuity: the "time discontinuity"
        // guard of the article, transposed to this joint law).
        bool bonded = false;
        double dn0 = 0.0;
        // ---- failure mode, recorded ONCE when D first reaches 1 ------------
        // The article reads its fig. 12 patterns as "tensile failure or mixed
        // tensile-shear mode". Telling the two apart needs the partition of the
        // eq. 16 damage driver AT the instant of breakage, which is not
        // recoverable afterwards (D saturates at 1 while the opening keeps
        // growing). failMode = rn^2 / (rn^2 + rs^2) in [0, 1]: 1 = pure mode I
        // (opening), 0 = pure mode II (sliding). -1 while the joint is intact.
        // Written to the joint VTU as `failMode` when writeJointMode = true.
        double failMode = -1.0;
        // ---- jointFailRule = majority (Solidity Y3Dfd.c, Sigma_tau) -------
        // Chez eux chaque point d integration porte SON z et SA traction, et
        // le joint ne meurt qu au-dela d UN point rompu (`nfail>1`, l. 1175).
        // Le joint 2D n a que deux points : la regle y demande donc les DEUX.
        // Inerte sous `any` (defaut) — le scalaire partage D pilote seul.
        double Dk[2] = {0.0, 0.0};
        // ---- strainRateDIFArm = continuous --------------------------------
        // Leur dpeftdif est recalcule a chaque pas (l. 1448-1456), jamais
        // fige : on garde donc les proprietes de BASE intactes et la loi
        // reconstruit ft = ftB*DIF a chaque pas, sans jamais composer.
        double ftB = 0.0, cohB = 0.0, GfB = 0.0, GfIIB = 0.0;
        bool baseSnapped = false;
    };

    struct BEdge {                         // active contact edge (elem face)
        int elem, na, nb;
    };

    // ---- brazilian loading --------------------------------------------------
    // brazilianLoading = platens (default) | traction.
    //
    // `platens` is what the established FDEM codes do (Y-Geo: Mahabadi et al.
    // 2010; the FEM-DEM BTS/UCS studies; calibration procedure of Tatone &
    // Grasselli, IJRMMS 75 (2015) 56-72): two rigid platens closing on the disc
    // at a constant velocity of order 0.05 m/s. The essential detail, and the
    // one this solver got wrong at first, is that BOTH platens move INWARD.
    // With a fixed lower platen the disc drifts, and whatever load the support
    // fails to take is silently carried by the Cundall damping — measured on
    // the first attempt: 1.39 MN read on top against 11 kN underneath, with the
    // centre of the disc never seeing the diametral field at all. Closing
    // symmetrically keeps the specimen still and makes the damping forces of
    // the two halves cancel.
    //
    // `traction` is the alternative kept for comparison: a prescribed pressure
    // growing at loadRate over two opposite rim arcs of half-angle loadArcDeg
    // (the ISRM curved-jaw idealisation). It is self-equilibrated by
    // construction, adds no contact stiffness to the stable time step and the
    // applied force is known exactly — but it is load-controlled, so it cannot
    // follow the post-peak.
    struct Platen {
        double y = 0.0, v = 0.0, halfW = 0.0, xc = 0.0;
        int sign = -1;                     // -1 above (pushes down), +1 below
        Eigen::Vector2d F{0.0, 0.0};       // force ON the platen
    };
    struct LoadArc {
        std::vector<BEdge> edge;           // rim faces carrying the traction
        double length = 0.0;               // total loaded length [m]
        Eigen::Vector2d F{0.0, 0.0};       // force applied this step [N]
    };

    // setup
    void buildMesh();                      // dispatch: grid | voronoi | file
    void buildMeshVoronoi();
    // mesh = file — import d'un maillage triangulaire NON STRUCTURE (Gmsh
    // MSH 2.2 ASCII, elements type 2) : le maillage "a la Yan et al.",
    // simplexes uniformes sans structure de grains (box uniquement).
    void buildMeshFile();
    void buildMeshDisc();                  // discMesh = native: exact rim
    void buildMeshShpb();                  // geometry = shpb: bar-disc-bar
    void setupShpbGauges();                // monitor-point element lists
    double shpbVel(double t) const;        // prescribed pulse [m/s]
    void shpbGaugeRead();                  // area-averaged eps_xx
    // cuts the disc out of the meshed box AND fits the resulting rim to the
    // circle (vpos is edited in place — the staircase rim is not a detail:
    // its flat facet at the pole destroys the diametral field)
    void cutDisc(std::vector<Eigen::Vector2d>& vpos,
                 std::vector<std::array<int, 3>>& tris,
                 std::vector<int>& triGrain) const;
    void buildFromTriangles(const std::vector<Eigen::Vector2d>& vpos,
                            const std::vector<std::array<int, 3>>& tris,
                            const std::vector<int>& triGrain,
                            const std::vector<int>& grainPhase);
    void assignJointProps();
    void applyJointStatistics();
    void applyJointSizeEffect(double mWeib);  // eq. 42 : ft *= (Zeff/V_J)^(1/m)
    void placeTool();
    // adaptive insertion machinery (insertion = adaptive)
    void buildBindingTables();
    void rebindVertex(int v);
    void insertionSweep();
    void activateJoint(int jI, double sig, double tau);
    // ---- DIF de Yang : application des facteurs a UN joint ---------------
    // Extrait de activateJoint() le 2026-08-25 pour etre partage avec le gel
    // a l AMORCAGE du schema intrinseque (difIntrinsic_). Arithmetique
    // inchangee : les reperes dif_yang_* de la suite verrouillent la
    // bit-identite du chemin adaptatif.
    void stampDif(Joint& J, double er);
    // dnE / dnF / slipF depuis les resistances COURANTES du joint — un seul
    // site, partage par assignJointProps, applyJointStatistics, stampDif et
    // refreshDif (voir jointDeltaC).
    void setJointLengths(Joint& J);
    void snapBase(Joint& J);              // proprietes de base, une fois
    void refreshDif(Joint& J, double er); // strainRateDIFArm = continuous
    void setupBrazilianLoad();
    void initPlatens(double xc, double hw);  // geometry + tributary weights
    void setupBoundaries();
    void setupConfinement();
    void setupExcavation();                // etude tunnel : relachement paroi
    void excavationForces();
    double excavRelief() const;            // facteur de relachement rel(t)
    void computeStableDt();

    // stepping
    void elementForces();
    void bodyForces();                     // gravity (config key `gravity`)
    void jointForces();
    void rebuildContactEdges();
    void generalContact();
    void toolContact();
    void brazilianForces();
    void platenForces();                   // shared by brazilian and UCS
    void confiningForces();
    void integrate();
    void computeFragments();

    // reporting
    double achievedConfinement() const;    // mean sigma_xx in the core [Pa]
    void discCentreStress(double& sxx, double& syy) const;

    Config cfg_;
    std::string out_;
    Material mat_;
    // Optional constitutive law for the BULK elements (law = dpr | saksala |
    // saksala2011 | dpdfh). Absent by default, in which case elements stay
    // linear elastic with the crush cap and every earlier result is
    // reproducible bit for bit. With a law set, this is the COUPLED
    // configuration of the thesis' vumat_fdem_coupled: a dissipative bulk
    // PLUS discrete cohesive joints. Plane strain is imposed by handing the
    // 3D law a Biot strain whose out-of-plane component is exactly zero.
    std::unique_ptr<MatLaw> law_;
    PhaseSet phases_;                      // per-phase materials (>= 1)
    Scenario scen_ = Scenario::PERCUSSION;

    double W_ = 0.2, H_ = 0.2, thk_ = 1.0;
    int nx_ = 64, ny_ = 64;
    double hmin_ = 1e-3;
    bool voronoi_ = false;
    int nGrains_ = 1;

    // disc geometry (geometry = disc; mandatory for scenario = brazilian).
    // discFlat_ is the FULL loading angle 2*alpha of the FLATTENED brazilian
    // disc: two horizontal chords are cut at |y - yc| = R cos(alpha), giving a
    // true flat bearing surface of width 2 R sin(alpha). 0 = plain disc.
    bool disc_ = false;
    double discR_ = 0.0, discFlat_ = 0.0;
    Eigen::Vector2d discC_{0.0, 0.0};

    // joint law (per-joint properties live in Joint; xiJ_ is the shared
    // compressive dashpot ratio)
    double xiJ_ = 0.05;

    // ---- post-peak softening shape (jointSoftening) ------------------------
    // linear (default, unchanged): ft -> 0 over dnF - dnE = 2 Gf / ft.
    // yan: the exponential reduction factor f(D) of Yan et al. (IJRMMS 169,
    // 2023) eq. 11, with the critical opening/slip calibrated so the area
    // under the softening branch is still exactly GfI and GfII (eq. 13, 15):
    //     dnF - dnE = Gf / (ft yanI_),  slipF = GfII / (c yanI_).
    // yanFricScaled_ selects whether the Coulomb term of the shear cap is
    // multiplied by f(D) as well: eq. 10 of the article does multiply it,
    // which leaves a fully broken joint with ZERO shear resistance under
    // compression. rockim keeps the residual friction by default (as the
    // linear law does), because a crushed joint sliding under compression
    // must stay a frictional contact; set jointFrictionScaled = true for
    // the literal form of eq. 10.
    bool yanSoft_ = false;
    bool yanFricScaled_ = false;
    // ---- jointResidualMu : le frottement RESIDUEL du joint rompu ---------
    // Coefficient de frottement vers lequel le joint GLISSE quand il
    // s endommage, par la meme f(D) que la cohesion : mu_eff = muRes +
    // (tan(frictionDeg) - muRes) f(D). A D = 0 c est le frottement de PIC, a
    // D = 1 c est muRes. C est la distinction de Y-Geo (AbuAisha et al. 2015,
    // eq. 7.5 : angle de frottement de FRACTURE phi_f distinct de l angle
    // interne), et l equivalent de ce que Solidity obtient en remettant le
    // joint rompu au contact et a son glissement (0,6 calcaire / 0,18
    // granite). rockim gardait jusqu ici le frottement de PIC a vie.
    // < 0 = non pose = comportement historique, bit-identique.
    // Generalise jointFrictionScaled : muRes = tan(frictionDeg) reproduit le
    // defaut, muRes = 0 reproduit jointFrictionScaled = 1. Les deux cles ne
    // peuvent donc pas etre posees ensemble (erreur de config).
    double muRes_ = -1.0;
    yan::Params yanP_;
    double yanI_ = 1.0;                    // int_0^1 f(D) dD

    // ---- decharge en CISAILLEMENT (jointShearUnload) -----------------------
    // plastic (defaut, inchange) : plasticite a retour — la decharge suit la
    //   secante de PENALITE et le glissement plastique est conserve.
    // origin : l'eq. 18 de Yan et al., symetrique de l'eq. 17 du mode I — la
    //   decharge et la recharge suivent la secante a l'ORIGINE passant par
    //   (s_max, tau_env(s_max)). C'est la forme litterale de l'article.
    //   Le glissement au pic est le s_p de Munjiza, (c + tan(phi) |sigma_n|)/p,
    //   evalue sur l'enveloppe NON endommagee et sur la part geometrique pj*dn
    //   de la contrainte normale : l'endommagement de mode II demarre donc au
    //   moment ou l'enveloppe de Mohr-Coulomb est atteinte, exactement comme
    //   dans le retour radial — les deux modes coincident en charge monotone
    //   et ne different qu'a la decharge.
    // ⚠️ l'eq. 18 place TOUT le cap dans la secante, frottement de Coulomb
    // compris ; combinee a jointFrictionScaled = 0 (defaut rockim, frottement
    // residuel non degrade) elle rend le glissement frottant reversible. La
    // configuration litterale de l'article est donc origin + fricScaled = 1
    // (c'est ce que pose configs_yan/article_exact_base.cfg).
    bool shearOrigin_ = false;

    // jointShearRange = cohesion | coulomb (defaut cohesion, bit-identique).
    // `coulomb` : la PLAGE d'adoucissement de mode II est divisee a chaque pas
    // par la resistance de Mohr-Coulomb a la pression courante
    // fs = c + tan(phi) |sigma_n| (compression seule), plancher 2 sE — la
    // convention du code Solidity (Y3Dfd.c l. 1110-1126 : st = max(2 sp,
    // 3 GfII/dpefs) avec dpefs friction-inclus, recalcule par point et par
    // pas). Sous fort confinement le glissement critique s'effondre de
    // 3 GfII/c (~1e2 um aux cartes d'impact) a quelques microns : c'est ce
    // qui autorise la rupture en cisaillement des joints COMPRIMES (noyau
    // broye d'indentation). `cohesion` = plage figee 3 GfII/c, l'historique.
    bool shearRangeCoulomb_ = false;

    // ---- adaptive insertion (insertion = adaptive) --------------------------
    // No cohesive joint exists at t = 0: every interior edge starts BONDED and
    // its co-located node copies move as ONE node (groups below), which is
    // bit-equal to the shared-node FEM of the article — no artificial bulk
    // compliance, no penalty in the stable-dt budget beyond the (softer)
    // activation penalty. Each step, insertionSweep() averages the two
    // neighbouring element stress tensors on every bonded edge and activates
    // the joint when sigma_n >= ft or |tau| >= fs (Mohr-Coulomb, eq. 7-8 of
    // the article). Node "splitting" (Fig. 7) falls out of re-running the
    // union-find at the two endpoint vertices: copies bind per connected
    // component of the element fan over still-bonded edges.
    bool adaptive_ = false;
    long nInserted_ = 0;
    std::vector<std::vector<int>> copiesOfVert_;   // vertex -> node copies
    std::vector<std::vector<int>> jointsOfVert_;   // vertex -> incident joints
    std::vector<std::vector<std::vector<int>>> grpsOfVert_; // vertex -> groups
    int nVert_ = 0;
    // per-phase element data (indexed by Elem::phase)
    std::vector<Eigen::Matrix3d> DmP_;
    std::vector<double> nuP_, crushCapP_, ftP_, rhoP_;
    // Constantes de Lame par phase, pour bulkModel = neohookean seulement
    // (la branche co-rotationnelle passe par DmP_, inchangee).
    std::vector<double> lamP_, mu2P_;
    std::vector<double> hEl_;              // per-element inscribed size (4A/P)

    // contact
    double kp_ = 0.0, muC_ = 0.5, xiC_ = 0.05, vReg_ = 1e-3;
    double kpGC_ = 0.0, xiGC_ = 0.4, gcWork_ = 0.0, relax_ = 1.0, gcRest_ = 0.2;
    std::unordered_map<uint64_t, double> pen0_;

    double damping_ = 0.02, pullV_ = 0.05;
    double bulkVisc_ = 0.0;                // 2*mu*D de l'eq. 6 (Yan) [Pa.s]
    // ---- DIF dependant du taux de deformation (Yang et al. 2025) --------
    // Leurs eq. 2-3 : un facteur d amplification dynamique applique a la
    // resistance en traction et a l energie de rupture mode I (DIF_traction),
    // a la cohesion et a l energie mode II (DIF_compression). L angle de
    // frottement N EST PAS touche (ils citent Zhao : influence negligeable).
    //
    // GELE A L INSERTION. Le schema extrinseque de Yan insere le joint a
    // l instant precis ou le materiau atteint son enveloppe : c est donc
    // l instant dont le taux de deformation gouverne la rupture, et figer le
    // DIF la evite de faire varier une resistance dont depend deja
    // l historique d endommagement. Note : puisque ft ET Gf recoivent le MEME
    // facteur, l ouverture critique dnF = dnE + kI Gf/ft est quasi invariante
    // (seule sa part elastique dnE = ft/pj suit le DIF), donc le compteur
    // d endommagement D = f(ouverture) n est pas corrompu.
    //
    // Exige insertion = adaptive : sans instant d insertion, pas de DIF.
    // ---- enveloppe de cisaillement : Yan (eq. 8) ou Yang (eq. 1) --------
    // Voir include/rockim/YangDif.hpp. Defaut yan = comportement historique,
    // bit-identique. La forme de Yang affaiblit le cisaillement dans les
    // zones TENDUES (elle y fait decroitre le terme de frottement au lieu de
    // l annuler), ce qui deplace le partage traction / cisaillement — donc le
    // facies. Les deux formes sont rigoureusement identiques en compression.
    bool yangEnv_ = false;
    // Facteur du cap de traction moyenne (pm <= f * ft). Historiquement 3, en
    // dur et sans echappatoire en 3D. Ce cap N EST DANS AUCUN des deux
    // articles : il faut pouvoir le desarmer pour faire tourner le modele de
    // quelqu un d autre. <= 0 le desarme. Defaut 3 = inchange.
    double mtCap_ = 3.0;
    bool difOn_ = false;
    // ---- DIF en schema INTRINSEQUE : le gel a l AMORCAGE (2026-08-25) -----
    // Vrai quand strainRateDIF est arme ET insertion = intrinsic. Le meme
    // balayage d enveloppe que l insertion adaptative est alors execute, mais
    // au franchissement il ne CREE pas le joint (il existe deja) : il se
    // contente d y STAMPER le DIF. Les deux schemas partagent donc le critere
    // exact et ne different que par sa consequence — c est le temoin recherche
    // pour la comparaison a Yang et al. (joints intrinseques AVEC DIF).
    // Combinaison auparavant REFUSEE par une exception : aucune configuration
    // valide ne change de comportement (principe I).
    // ---- bulkModel : la loi de VOLUME (2026-08-25) ---------------------
    // `corotational` (defaut, historique) : Biot + P = R sigma.
    // `neohookean` : Guo (these Imperial 2014, eq. 2.6),
    //   T = (mu/J)(B - I) + (lambda/J) ln(J) I,   B = F F^T,
    // avec l assemblage EXACT P = J T F^-T = R sigma cof(U). La loi est
    // hyperelastique (W = mu/2 (tr B - 3) - mu ln J + lambda/2 (ln J)^2,
    // le neo-hookeen compressible de Simo-Hughes) et redonne l elasticite
    // lineaire au premier ordre avec les MEMES lambda et mu.
    // ATTENTION : le facteur d ecart a la forme co-rotationnelle vaut
    // J^(-2/3) en 3D mais J^(-1/2) en DEFORMATION PLANE. On n ecrit donc
    // JAMAIS l exposant en dur : la forme generique cof(U) = J U^-1 est
    // correcte dans les deux dimensions.
    // ---- Loi de joint : les deux dernieres conventions de Guo -------------
    // jointElastic = linear (defaut) | parabolic : la branche elastique du
    // joint. Guo eq. 2.31 la pose PARABOLIQUE, sigma = ft(2r - r^2) avec
    // r = dn/dnE, ce qui donne une tangente NULLE au pic (transition douce
    // vers l adoucissement) et une pente initiale 2 pj, valable des DEUX
    // cotes de dn = 0 (loi C1 a l origine). rockim posait une droite de
    // pente pj, avec un coude au pic.
    // ---- jointQuadrature : les points d integration du joint -------------
    // `vertex` (defaut) : aux NOEUDS. C est la regle de Newton-Cotes nodale,
    // celle des elements cohesifs d Abaqus, retenue contre les oscillations
    // parasites du champ de traction a forte penalite (Schellekens & de
    // Borst). `midedge` : aux MILIEUX D ARETES, poids 1/3 — la regle de Guo
    // (Table 2.2), exacte a l ordre 2 au lieu de 1.
    // Les deux coincident EXACTEMENT en chargement uniforme ; elles different
    // de 50 a 200 % la ou l ouverture a un gradient a travers la facette,
    // c est-a-dire au front de fissure. Ce n est donc PAS un raffinement de
    // second ordre : c est un effet nul la ou on le mesure d habitude et fort
    // la ou la fissure se decide.
    bool midEdge_ = false;
    bool paraElastic_ = false;
    // jointDeltaC = exact (defaut) | guo : l ouverture critique. Guo
    // eq. 2.30 pose delta_c = 3 Gf/f mesure depuis ZERO, en approchant
    // l integrale de la z-curve par 1/3 la ou elle vaut 0,386307. Son
    // modele dissipe donc 3/0,386307 = 1,159 fois son Gf nominal. C est
    // SA convention, et ses Gf publies ont ete calibres avec : il faut la
    // reproduire pour retrouver ses chiffres.
    bool guoDeltaC_ = false;
    // jointDeltaC = solidity : ce que fait REELLEMENT leur code, qui n est
    // pas ce qu ecrit la these. Y3Dfd.c lignes 1098-1099 et 1125-1126 :
    //     op = 2 el ft / p0                   (l ouverture au pic)
    //     ot = MAXIM(2 op, 3 Gf/ft)           (la PLAGE d adoucissement)
    // la rupture est donc a op + ot, pas a 3 Gf/ft : la convention `guo`
    // oublie et l offset op et le plancher 2 op. Le plancher mord quand le
    // maillage est fin devant Gf/ft — exactement le regime d un impact.
    bool solidityDeltaC_ = false;
    bool neoHooke_ = false;
    bool difIntrinsic_ = false;
    // strainRateDIFArm = continuous : DIF recalcule a CHAQUE pas depuis les
    // proprietes de base, comme leur dpeftdif. Exclusif de insertion/envelope.
    bool difContinuous_ = false;
    // jointFailRule = majority : mort du joint seulement au-dela d UN point
    // d integration a D >= 1 (leur `nfail>1`), donc les DEUX points en 2D.
    bool majorityFail_ = false;
    // ---- gcBirth : la naissance d un contact sur un joint qui vient de mourir
    // `ramp` (defaut, historique) : la force part de ZERO et monte sur
    // gcBirthTau (releve de naissance par AIRE, voir PotHist::aRef).
    // `penalty` : leur solution (Y3Did.c l. 915-964) — la penalite de la paire
    // est calee sur fn_joint/fn_contact, bornee, pour que la force soit
    // CONTINUE. Le facteur persiste, la raideur tangentielle le suit.
    bool birthPenalty_ = false;
    double birthPenMin_ = 0.01, birthPenMax_ = 3.0;   // leurs deux bornes
    long nBirthScaled_ = 0;                           // mesure : paires calees
    double birthScaleSum_ = 0.0;                      // ... et facteur moyen
    // ---- strainRateFilter : le taux qui alimente le DIF ------------------
    // `exponential` (defaut) : passe-bas de constante strainRateTau.
    // `none` : le taux BRUT, ce que fait leur code (Y3Dfd.c l. 1448).
    bool srFilterOff_ = false;
    long nDifStamped_ = 0;                 // joints ayant recu le gel
    long nDifLate_ = 0;                    // dont DEJA endommages au gel
    // ---- jointDeath : QUAND le joint passe la main au contact -------------
    // `separation` (defaut, historique) : le joint ne meurt qu une fois
    // FRANCHEMENT ouvert (dnMax > 3 dnF). En compression il ne meurt donc
    // JAMAIS, et l algorithme de contact — qui porte le glissement contactMu,
    // le 0,6 de leur Table 4 — ne prend jamais le relais sous l insert. Le
    // commentaire du site de mort explique le choix : un joint broye qui
    // glisse en compression reste VIVANT et sert de contact frottant de ses
    // propres levres. C est defendable, a une chose pres — il frotte alors a
    // tan(frictionDeg), le frottement de PIC, et non au glissement residuel.
    // `damage` : mort des que D >= 1, quel que soit le signe de l ouverture,
    // la regle de Guo (these Imperial 2014, §2.3.3). Le risque documente est
    // la pompe a energie (le contact materialise 1/2 k pen^2 sur des levres
    // interpenetrees) ; la releve de naissance pen0_ ajoutee depuis devrait
    // la neutraliser — c est ce que la mesure fDeath doit trancher.
    bool deathOnDamage_ = false;
    double srTau_ = 0.0;                   // constante du filtre de taux [s]
    double srRelax_ = 0.0;                 // exp(-dt/srTau_), pose apres dt
    // Exposant de DIF_traction. 0,07 = transcription LITTERALE de leur eq. 3.
    // 0,1707 = valeur deduite de LEUR PROPRE figure 2(b), qui rend la loi
    // continue a ses deux bornes (voir la note dans FdemSolver.cpp).
    double difExpT_ = 0.07;
    double gravity_ = 0.0;                 // body-force acceleration, -y [m/s^2]
    double pullRamp_ = 0.0;                // grip velocity rise time [s]
    bool gripFree_ = false;                // frictionless tension grips

    std::vector<Eigen::Vector2d> X0_, u_, v_, f_;
    std::vector<double> m_;
    std::vector<int> flag_, elemOf_, vOf_; // node flag, owning element, origin vertex
    std::vector<Elem> el_;
    std::vector<Joint> jt_;
    std::vector<int> fragId_;              // per element
    std::vector<char> extEdge_;            // per joint? exterior edges kept apart
    std::vector<BEdge> exterior_;          // domain-boundary faces

    // quiet boundaries (per node)
    std::vector<double> cAbsX_, cAbsY_, kAbsX_, kAbsY_;

    // general contact machinery
    std::vector<BEdge> act_;               // active edges (exterior + freed)
    std::vector<int> actNodes_;

    // ---- activation adaptative du contact (gcActivation = adaptive) --------
    // Fukuda et al. (GPGPU-FDEM) : act_ ne contient que les faces qui PEUVENT
    // toucher, au lieu de tout l'exterieur balaye a chaque pas. Une face du
    // pool (l'exterieur, fenetre par gcXwindow) s'active — DEFINITIVEMENT,
    // l'activation est monotone — par l'une des trois regles :
    //   (C) peau endommagee : son element porte un joint casse (D >= 1) ou
    //       mort — le voisinage immediat d'une fissure, active d'office ;
    //   (A) autre corps : sa sphere englobante approche a moins de gcActMargin
    //       cellules d'une face d'une AUTRE composante connexe (union-find sur
    //       les joints porteurs : bonded ou vivants avec D < 1, recalcule
    //       quand nBroken_ change) — c'est ce qui arme le SHPB multi-corps des
    //       t = 0 et les debris detaches qui retombent ;
    //   (B) voisinage d'un contact reel : idem, contre une face active dont
    //       les noeuds ont deja PORTE une force de contact (estampille
    //       lastTouch_) — un debris qui racle propage l'activation a la
    //       surface qu'il aborde, une face active mais inerte ne propage rien
    //       (sans ce filtre, l'adjacence a gap nul contaminait tout le bord de
    //       proche en proche).
    // Le balayage est periodique, cadence par v_max : entre deux balayages
    // deux surfaces ne peuvent se rapprocher de plus de 2 v_max dt par pas,
    // la cadence est bornee pour que l'approche reste sous la demi-marge
    // (facteur de securite 2). Il se declenche aussi des que nBroken_ change
    // (aligne sur la cadence du rebuild, % 8). Cout : O(N) tous les
    // gcActEvery pas — mesure nul devant un seul pas de detection.
    // Approximation ASSUMEE (la meme que Fukuda) : un continuum INTACT ne
    // peut pas se replier sur lui-meme. gcActivation = full (defaut,
    // bit-identique) reste la reference.
    bool gcAdaptive_ = false;
    double gcActMargin_ = 2.0;             // marge d'activation [cellules]
    long gcActEvery_ = 64;                 // cadence max du balayage [pas]
    std::vector<BEdge> pool_;              // exterieur fenetre (fige a l'init)
    std::vector<char> extOn_;              // drapeau actif, par face du pool
    std::vector<char> elemDam_;            // regle C : element au bord casse
                                           // + un anneau par sommet
    std::vector<std::vector<int>> vElems_; // sommet -> elements (statique)
    std::vector<int> bodyOf_;              // composante connexe par element
    std::vector<long> lastTouch_;          // dernier pas ou le NOEUD a porte
                                           // une force de contact (-1 jamais)
    long bodyStamp_ = -1;                  // nBroken_ au dernier union-find
    int nBodies_ = 1;
    long nextSweep_ = 0;                   // stepCount_ du prochain balayage
    long sweepBroken_ = -1;                // nBroken_ au dernier balayage
    bool poolBuilt_ = false;
    bool haveDead_ = false;                // des faces liberees existent
    long nActivated_ = 0;                  // compteur (resume de fin de run)
    std::vector<int> actPool_;             // act idx -> pool idx (-1 liberee)
    std::vector<BEdge> deadList_;          // cache des faces liberees, ne se
                                           // rafraichit QUE sur le declencheur
                                           // historique (timing = mode full)
    std::vector<long> poolTouch_;          // premier pas avec force (diag)
    std::vector<long> actStep_;            // pas d'activation (diag)
    void activationSweep();

    // ---- contact par POTENTIEL de Munjiza (contact = potential) — A3 -------
    // La forme des eq. 2-5 de Yan et al. 2023 : force normale distribuee
    // = p ∮ (phi_A - phi_B) n dG sur le bord du RECOUVREMENT des deux
    // triangles (conservative par construction, 3e loi de Newton machine —
    // voir PotentialContact.hpp), frottement tangentiel INCREMENTAL a ressort
    // avec cap de Coulomb (eq. 4-5). Paires d'ELEMENTS (pas noeud-arete),
    // detection O(N) type NBS : binning AABB en grille de cellules. Une paire
    // liee par un joint VIVANT est exclue (le joint porte son interaction,
    // contact du joint casse compris) ; un joint MORT rend la paire au
    // contact general. Compose avec gcActivation et gcXwindow (le jeu actif
    // fournit les elements candidats). Le contact OUTIL analytique reste en
    // penalite (l'outil n'est pas maille — leve par B1, physical groups).
    bool contactPot_ = false;
    double potP_ = 0.0;                    // penalite normale [Pa.m = N/m.phi]
    double potKt_ = 0.0;                   // raideur tangentielle (eq. 4-5)
    // Amortissement NORMAL du contact potentiel. OPT-IN, defaut 0 = aucun
    // terme calcule, trajectoires bit-identiques (le selftest-potential2d
    // teste la CONSERVATION, il doit rester exact).
    //
    // Pourquoi il manquait. La branche PENALITE se donne deja xiGC (0,8)
    // et une restitution gcRest (0,2), avec la justification ecrite en
    // toutes lettres a l'initialisation : a la raideur pleine E*t le
    // contact de debris  pompe de l'energie dans la zone broyee . Le mode
    // POTENTIEL ignorait les deux : il est conservatif par construction
    // (voir PotentialContact.hpp), donc un nuage de fragments n'y perd
    // jamais rien et ne se pose jamais. Mesure du 2026-08-18 (coupe PDC
    // out_cut_v3) : 202 fragments apparaissent en 10 us et 37 kJ/m sont
    // dissipes pendant que l'outil n'a fourni que 60 J/m.
    double potXi_ = 0.0;                   // fraction critique, masse reduite
    struct PotHist {                       // historique par paire d'elements
        Eigen::Vector2d Ft{0.0, 0.0};      // force tangentielle sur l'el. MIN
        long step = -1000;                 // dernier pas vu (peremption du Ft)
        double aRef = 0.0;                 // releve de naissance par AIRE de
                                           // recouvrement — le pen0_ exact du
                                           // potentiel. A la premiere vue,
                                           // aRef = aire courante ; ensuite il
                                           // decroit en exp(-t/gcBirthTau) et
                                           // la force est multipliee par
                                           // max(0, 1 - aRef/aire) :
                                           //  * paire nee en APPROCHE reelle :
                                           //    aRef ~ 0 des la naissance, la
                                           //    conservation est intacte ;
                                           //  * paire nee EN RECOUVREMENT (un
                                           //    joint meurt comprime dans la
                                           //    zone broyee) : seule la
                                           //    NOUVELLE approche resiste, et
                                           //    la sortie sous aRef est LIBRE
                                           //    — asymetrie du signe ABSORBANT.
                                           // Une rampe TEMPORELLE fait
                                           // l'inverse (approche bradee,
                                           // sortie plein tarif) et INJECTE
                                           // sur chaque premier cycle —
                                           // mesure : +936 J/m sans releve,
                                           // +179 avec rampe temporelle,
                                           // asymetrie d'etat requise.
                                           // Jamais reinitialise (parite
                                           // pen0_) : une paire retrouvee
                                           // reprend sa force pleine.
        // ---- gcBirth = penalty (Solidity Y3Did.c l. 915-964) --------------
        // Facteur de penalite PROPRE A CETTE PAIRE, fige au pas de naissance
        // pour que la force du contact naissant egale celle du joint mourant.
        // C est leur d1pepe[icoup]. < 0 = pas encore ne. Inerte sous `ramp`.
        double penScale = -1.0;
    };
    std::unordered_map<uint64_t, PotHist> potFt_;
    std::unordered_map<uint64_t, int> jointOfPair_;   // (eMin,eMax) -> joint
    void potentialContact();
    // OpenMP scratch: per-thread force buffers for the joint scatter
    std::vector<std::vector<Eigen::Vector2d>> fTL_;
    std::vector<std::vector<char>> seenTL_;
    std::vector<std::vector<int>> touchedTL_;
    long actStamp_ = -1;                   // estampille du dernier rebuild :
                                           // nBroken_ en mode legacy, nDead_
                                           // en mode eager (voir gcEager_)
    // --- A' : rafraichissement des surfaces liberees ------------------------
    // gcSurfaceRefresh = legacy (defaut, bit-identique) | eager
    //
    // En legacy le cache deadList_ n'est rafraichi que si nBroken_ CHANGE et a
    // un pas multiple de 8. Or nBroken_ compte les CASSES (D -> 1) tandis que
    // J.dead marque les SEPARATIONS (dnMax > 3 dnF) : deux evenements
    // decouples dans le temps. Un joint qui se separe alors que plus rien ne
    // casse ne voit donc JAMAIS ses faces declarees — latence non bornee,
    // interpenetration libre, puis impulsion maximale a la declaration.
    // Mesure du 2026-08-18 (coupe PDC, bissection RKM_NOGC) : ce canal porte
    // un facteur 6 a 8 sur toutes les observables (vitesse nodale 2544 -> 314
    // m/s, fragments 368 -> 45, injection outil 77 -> 12 kJ/m).
    // En eager l'estampille suit nDead_ et la grille % 8 saute : les faces
    // entrent au pas ou leur joint meurt.
    // ATTENTION : eager rompt volontairement l'identite au bit pres entre
    // modes full et adaptive que le cache legacy preserve.
    bool gcEager_ = false;
    long nDead_ = 0;                       // joints separes (J.dead), cumul

    // --- PENALITE DE CONTACT ADAPTATIVE (EPFL) ------------------------------
    // jointContactPenalty = fixed (defaut, bit-identique) | adaptive
    //
    // SOURCE : T. Ghesquiere-Dierickx, J.-F. Molinari & G. Anciaux, « Stability
    // of Extrinsic Cohesive-Zone Model with Penalty-Based Contact in Explicit
    // Dynamic Fragmentation Simulations », arXiv:2511.14323v1 [physics.comp-ph],
    // 18 novembre 2025, EPFL. Leur section 4, « Adaptive Contact Penalty » :
    //
    //     « we replace the fixed contact penalty k- with a damage-dependent
    //       variant following the cohesive stiffness, i.e. k- = k+(d),
    //       thereby restoring continuity at delta = 0 and suppressing any
    //       energy errors. »
    //
    // LE MECANISME VISE. En compression le joint applique la penalite PLEINE
    // J.pj, alors qu'en traction il applique la secante endommagee (1-D) pj.
    // Le saut de raideur a dn = 0 est la discontinuite que leur article
    // designe comme « the most critical source of instability » : chaque
    // traversee de dn = 0 injecte ou dissipe un increment, et les traversees
    // repetees derivent. Mesure du 2026-08-18 sur la coupe PDC : couper le
    // contact general divise la vitesse nodale par 8, mais rafraichir les
    // surfaces PLUS TOT l'aggrave d'un facteur 2 a 3 — l'injection suit le
    // NOMBRE de basculements, pas leur retard. C'est leur diagnostic.
    //
    // TRANSCRIPTION EXACTE. Leur k+(d) est la secante du modele d'endommagement.
    // Dans rockim la traction d'essai vaut tr = (1 - D) pj dn, donc la secante
    // est (1 - D) pj — et a D = 0 elle vaut pj, si bien qu'a l'etat neuf il n'y
    // a AUCUNE discontinuite, ce qui est la propriete recherchee.
    //
    // CE QUE LES AUTEURS EN DISENT EUX-MEMES, et qu'il faut lire avant d'y voir
    // un remede : « non-physical interpenetration becomes more pronounced at
    // higher cohesive damage levels, since the contact stiffness decays with
    // increasing damage. Consequently, fragment statistics, although
    // numerically stable, cannot be regarded as physically accurate. [...] it
    // cannot be viewed as a viable long-term remedy. » C'est un INSTRUMENT DE
    // DIAGNOSTIC, pas une correction de production.
    //
    // LE PAS DE TEMPS RESTE CELUI DE LA PENALITE PLEINE : computeStableDt()
    // budgete J.pj, non la secante. C'est deliberement leur protocole
    // (« with dt_c,G still computed using the originally fixed k- »).
    bool jcAdaptive_ = false;

    // --- A : ECRETAGE EN IMPULSION DU CONTACT OUTIL --------------------------
    // toolImpulseCap = kappa >= 0 ; 0 (defaut) desarme, comportement historique.
    // La force outil sur un noeud est bornee par
    //
    //     |Fc| <= kappa * 2 * |v_outil| * m_i / dt
    //
    // c'est-a-dire que l'outil ne peut pas changer la vitesse d'un noeud de
    // plus de kappa * 2 v_outil en un pas.
    //
    // POURQUOI CETTE BORNE. Un corps rigide de masse INFINIE anime de v
    // communique au plus 2v a une masse ponctuelle, en choc parfaitement
    // elastique. Un outil a vitesse IMPOSEE est exactement ce corps : il ne
    // decelere jamais, donc c'est un reservoir d'energie infini. Tout ce qui
    // depasse 2v est numerique par construction. kappa = 1 est la borne
    // stricte ; kappa > 1 la relache si un cas l'exige.
    //
    // CE QU'IL CORRIGE, MESURE LE 2026-08-18 SUR LA COUPE PDC. Le solveur
    // imprime deja les deux nombres qui exposent la pompe :
    //   outil->solide  77 286 J/m   (ce que l'outil injecte dans le solide)
    //   tool work out     189 J/m   (ce qu'il depense en corps rigide)
    // soit un rapport de 408. L'ecretage geometrique existant plafonne la
    // force a kp_ * 0,6 h = 7,24 MN/m, ce qui donne encore 377 m/s par pas sur
    // une masse nodale de 2,46e-5 kg/m — 992 m/s sur la plus legere.
    //
    // CE QU'IL NE CORRIGE PAS : ni le contact general, ni les joints, ni
    // l'interpenetration. L'audit du 2026-08-18 montre d'ailleurs que celle-ci
    // est bornee a 0,6 h partout, donc sous controle.
    //
    // TEST D'ACCEPTATION, binaire : aucun noeud au-dessus de 2 v_outil, et le
    // rapport injection / travail de corps rigide qui retombe vers 1.
    double toolVCap_ = 0.0;

    // --- A1 : CONTACT OUTIL EN CONDITION DE VITESSE (CD-Lagrange) ----------
    // toolContact = penalty (defaut, bit-identique) | signorini
    //
    // SOURCES. Fekak, Brun, Gravouil et al. (2017), schema CD-Lagrange ;
    // expose complet dans Dureisseix, Greco, Brun et al., « Explicit dynamics
    // and non-smooth interface behaviors », JTCAM (2024), leur Algorithme 1 ;
    // et Ghesquiere-Dierickx, Anciaux, Acary & Molinari, arXiv:2606.01355
    // (Nonsmooth Newmark), qui applique la meme idee a la fragmentation avec
    // modele cohesif extrinseque — exactement notre cas.
    //
    // CE QUE CA REMPLACE. La voie penalite borne la PENETRATION (0,6 h) et ne
    // borne pas l'IMPULSION, qui est la grandeur qui lance les noeuds. Mesure
    // du 2026-08-18 (coupe PDC) : l'outil injecte 77 286 J/m dans le solide
    // pour 189 J/m de travail de corps rigide — rapport 408 — et lance des
    // noeuds a 2 544 m/s contre une borne physique de 2 v_outil = 20 m/s,
    // soit 127 fois trop. L'ecretage toolImpulseCap pose le meme jour plafonne
    // l'increment PAR PAS, alors que la borne physique porte sur l'impulsion
    // de TOUTE la collision : un noeud en contact soutenu prend 2 v a chaque
    // pas et seize pas donnent les 311 m/s mesures.
    //
    // LA FORMULATION. Le contact n'est plus une force penalisee mais une
    // relation IMPULSION / SAUT DE VITESSE (lemme de viabilite de Moreau) :
    //     si g > 0  alors r = 0 ;  sinon  0 <= v ⊥ r >= 0
    // avec la dynamique reduite v = v_libre + H r, H = L M^-1 L^T l'operateur
    // de Delassus.
    //
    // POURQUOI C'EST EXPLICITE ICI. L'article le dit : « due to the properties
    // of the lumped mass matrix M, this leads to a Delassus operator which is
    // diagonal, definite positive and spherical per node ». La masse de rockim
    // EST diagonale et l'outil est un obstacle RIGIDE : L est local a chaque
    // noeud, H vaut 1/m_i, et l'impulsion se calcule en FORME FERMEE noeud par
    // noeud. Aucun systeme a resoudre, le solveur reste matrix-free.
    //
    // CONSEQUENCE SUR LE PAS DE TEMPS. kp_ sort du budget de computeStableDt()
    // en mode signorini : la borne ressort-masse perd son terme de contact
    // outil. C'est le gain structurel du schema — le pas cesse de dependre
    // d'une raideur de penalite arbitraire.
    //
    // CE QU'IL FAUT ACCEPTER. Une interpenetration RESIDUELLE subsiste : la
    // condition porte sur la vitesse, pas sur le deplacement. Les auteurs la
    // quantifient par eta = -min(g)/max|g|, 0,43 % sur leur cas de reference.
    // toolSignoriniRelax > 0 en resorbe une fraction par pas (rattrapage facon
    // Baumgarte) ; 0 = condition de vitesse PURE, le defaut, celui de
    // l'article.
    //
    // LIMITE CONNUE ET ASSUMEE : la vitesse libre est formee avec f_ tel qu'il
    // est a l'appel de toolContact(). En percussion et en coupe, c'est le
    // total (l'outil est la derniere famille de forces). Si un confinement est
    // arme, sa contribution — rampee et lente — n'y figure pas ; l'erreur est
    // du second ordre mais elle existe.
    bool toolSig_ = false;
    double toolSigRelax_ = 0.0;

    // --- BALAI NUMERIQUE : retrait des fragments ----------------------------
    // SOURCE : X. Yang, J. Xiang, S. Naderi, Y. Wang, J. Aising, I. Ugarte &
    // J.-P. Latham, « Multi-criteria validation of hi-fidelity numerical model
    // of impact breakage: towards next generation percussion drill
    // simulation », Int. J. Rock Mech. Min. Sci. 191 (2025) 106125, sec. 2.
    // Donnees experimentales de Mines Paris-PSL (Aising, Gerbaud, Sellami).
    //
    // LE DEFAUT CORRIGE. Un fragment identifie par la seule CONNECTIVITE peut
    // etre topologiquement detache et physiquement prisonnier : « some
    // simulated fragments may be completely surrounded by failed joint
    // elements but still mechanically constrained by adjacent intact rock
    // blocks ». computeFragments() declare detache TOUT ce qui n'est pas la
    // plus grosse composante — c'est exactement ce critere naif, et il
    // SUREVALUE le volume detache d'une quantite jamais mesuree.
    //
    // LE PRINCIPE. On reproduit le geste experimental : apres l'impact, un
    // petit pinceau balaie les debris. On donne donc a tous les candidats une
    // vitesse initiale et une acceleration « anti-gravite », on integre un
    // court instant, et CE QUI BOUGE est du debris, CE QUI RESTE COINCE n'en
    // est pas. Ce ne sont ni un seuil ni un critere qui tranchent, c'est le
    // contact et le frottement.
    //
    // L'OBJET DE REFERENCE, l'idee centrale. On ne compare pas a un
    // deplacement absolu mais a celui d'une PARTICULE LIBRE fictive soumise
    // aux memes vitesse initiale et acceleration, donc purement balistique :
    //     d_ref = v0 t + 1/2 a t^2
    // Un fragment est retire si son deplacement depasse beta * d_ref. La
    // question posee n'est pas « a-t-il bouge ? » mais « a-t-il bouge aussi
    // librement que ce que rien ne retient ? ».
    //
    // LEURS VALEURS, reprises en defaut : v0 = 2,5 mm/s, a = 10 g = 98,1 m/s2,
    // duree 0,5 ms, beta = 0,8. Soit d_ref = 13,5 um et un seuil de 10,8 um.
    // beta = 0,8 est justifie par un PLATEAU : leur figure 4 balaie 0,2 a 1,0
    // et le classement cesse de changer au-dela de 0,8.
    //
    // POURQUOI CE SEUIL ETAIT INATTEIGNABLE HIER. 10,8 um, contre un pas de
    // lecture des VTU de 10 um avant le correctif de precision du 2026-08-18.
    // Le critere tombait pile sur la quantification.
    //
    // CONTRAINTE SUR L'ACCELERATION : « the anti-gravity force should not
    // cause any new cracks ». Verification pour le granite de l'article :
    // rho a h = 2626 x 98,1 x 0,15 = 38,6 kPa contre ft = 10,98 MPa, marge
    // d'un facteur 284. La condition est tres largement satisfaite — ce qui
    // laisse de la place pour accelerer le balayage si son cout derange.
    //
    // COMPTABILITE. Le balai INJECTE de l'energie deliberement. Son travail va
    // dans un poste SEPARE (brushWork_), affiche a part et signale comme tel :
    // apres la journee du 2026-08-18 passee a traquer des pompes logees dans
    // des canaux comptabilises, un balai numerique qui alimenterait sumW en
    // silence serait exactement l'erreur a ne pas commettre.
    //
    // RIEN N'EST SUPPRIME : « the rock fragments are not physically deleted ».
    // C'est une CLASSIFICATION. Le run reste physique, seule la masse de
    // fragments rapportee change.
    //
    // Cles : fragBrushStart (0 = desarme), fragBrushV0, fragBrushAccel,
    //        fragBrushBeta, fragBrushDirX, fragBrushDirY.
    // --- ETAPE 1 DE L'ALGORITHME : « After completing the impact simulation »
    // Yang et al. 2025, sec. 2.3, etape 1. Le balai suppose que le chargement
    // est TERMINE : il mesure un deplacement et l'attribue a l'anti-gravite,
    // ce qui n'est valide que si rien d'autre ne bouge les fragments.
    //
    // Leur geometrie satisfait cette condition gratuitement : le taillant
    // frappe, rebondit, s'en va. Une coupe continue n'a pas d'« apres » — il
    // faut le FABRIQUER. D'ou cette cle, distincte de fragBrushStart : on
    // arrete l'outil a toolStop, on laisse l'amortissement eteindre les
    // vitesses, puis on arme le balai plus tard.
    //
    // MESURE DU 2026-08-18 QUI L'IMPOSE : en confondant les deux instants,
    // trois fragments volaient encore a plus de 50 m/s a l'armement, dont un
    // a 108,5 m/s VERS LE BAS. En 140 us de balayage il a traverse 18 mm d'un
    // bloc qui en fait 20 — soit 15,2 mm de balistique pure, le compte y est.
    // Il a ete classe « libre » sans que le balai y soit pour rien.
    //
    // L'ECHELLE DE VITESSE DU BALAI vaut v0 + a t_balayage. Toute vitesse
    // residuelle du meme ordre noie le signal : c'est ce que armBrush()
    // imprime desormais, pour que le journal dise lui-meme si l'etape 1 est
    // respectee.
    double toolStop_ = 0.0;                // 0 = l'outil ne s'arrete jamais
    bool toolStopped_ = false;
    double brushStart_ = 0.0;              // instant d'armement (0 = desarme)
    double brushV0_ = 2.5e-3;              // m/s
    double brushA_ = 98.1;                 // m/s2 (10 g)
    double brushBeta_ = 0.8;
    Eigen::Vector2d brushDir_ = Eigen::Vector2d(0.0, 1.0);
    bool brushZeroV_ = false;              // remplacer v au lieu d'ajouter
    bool brushArmed_ = false;
    double brushT0_ = 0.0;                 // instant reel de l'armement
    double brushWork_ = 0.0;               // poste d'energie SEPARE
    std::vector<char> brushCand_;          // par ELEMENT : candidat au retrait
    std::vector<int> brushFrag_;           // fragId_ fige a l'armement
    int brushNFrag_ = 0;
    std::vector<Eigen::Vector2d> brushU0_; // deplacement nodal de reference
    void armBrush();
    void brushReport();
    double cell_ = 0.0;
    Eigen::Vector2d gmin_;
    int gx_ = 1, gy_ = 1;
    std::vector<std::vector<int>> grid_;   // edge ids per cell

    Tool tool_;
    double toolKE0_ = 0.0;

    // UCS / triaxial loaded through PLATENS instead of clamped node rows.
    // Clamping a node row imposes v = 0 over zero thickness, which creates a
    // boundary layer the mesh resolves better and better as it is refined: on
    // the unstructured intra-grain mesh, 100 % of the first fifty joints broke
    // within 2.5 mm of the fixed grip, at two wave transits, and the specimen
    // read six times ft. Y-Geo and the FEM-DEM calibration literature load UCS
    // through rigid platens in frictional contact, which has no clamped layer
    // at all and lets the specimen expand laterally. loading = grips restores
    // the old behaviour.
    bool tensionPlatens_ = false;
    bool brazPlatens_ = true;              // platens (default) vs traction arcs
    double kpPlaten_ = 0.0;                // platen contact penalty [N/m]
    // Per-node weight of the platen contact = tributary rim length / mean.
    // Without it the penalty is applied once per DUPLICATED node, so a corner
    // where k elements meet gets k springs and the bearing pressure is
    // distributed by mesh VALENCE (3 to 8 on a voronoi rim) instead of by
    // geometry. Zero for nodes that are not on the exterior, which also kills
    // the spurious contact of interior duplicates. platenTributary = false
    // restores the old per-node behaviour for comparison.
    bool platenTrib_ = true;
    std::vector<double> platenW_;
    // participation ratio (sum f)^2 / sum f^2 of the bearing nodes = effective
    // number of nodes carrying the load; latched with the elastic gauge
    double gPR_ = 0.0;
    long gNC_ = 0;
    // Diffuse SUB-CRITICAL damage at the moment of peak load. The intrinsic
    // penalty makes the elastic range of a joint microscopic (mode I opens
    // elastically over ft/pj, mode II slips past coh/pj — a few nanometres at
    // pj = 20 E/h), so joint vibration can ratchet a damage field across the
    // whole specimen long before any crack exists. If that field is what eats
    // the brazilian strength, its extent must correlate with the deficit.
    double gDfrac_ = 0.0, gDmean_ = 0.0;
    // Work done by the joint dashpot. A dashpot can only DISSIPATE, so this
    // must stay <= 0. It is the direct detector of the rectifier failure mode
    // that a one-sided viscous term produces, and it costs one multiply per
    // integration point.
    double dampWork_ = 0.0;
    // ---- V2/B4 : bilan d'energie par sous-systeme (miroir du 3D) ----------
    // Instrumentation PURE : KE(t) - KE(0) = elWork_ + jointWork_ + gcWork_
    // + cundWork_ + lysWork_ + toolWork_ + bcWork_ + residu O(dt).
    double elWork_ = 0.0;      // forces internes des elements
    // Part VISQUEUSE (2 mu D) du travail des elements, comptee a part. Ce
    // n est PAS un poste supplementaire du bilan : elle est deja DANS
    // elWork_ (la contrainte visqueuse produit les memes forces nodales).
    // On l isole parce que sans elle le poste  elements  melange l energie
    // elastique STOCKEE et la dissipation visqueuse, et le budget devient
    // illisible des que bulkViscosity est arme.
    double viscWork_ = 0.0;
    // WP1 pulverisation — miroir 2D des membres du 3D. bdWork_ est une
    // ENERGIE (somme des Y dD), ventilee dans elWork_.
    bool bdOn_ = false;
    double bdD0_ = 1.4e-5, bdDf_ = 4.0e-4, bdDmax_ = 0.9, bdCd_ = 1.0;
    double bdWork_ = 0.0;
    long nPulv_ = 0;
    double jointWork_ = 0.0;   // tractions des joints (visqueux INCLUS)
    double cundWork_ = 0.0;    // damping local de Cundall (<= 0)
    double lysWork_ = 0.0;     // frontieres absorbantes (ressort+amortisseur)
    double gcFricWork_ = 0.0;  // part tangentielle du contact general
    double toolWork_ = 0.0;    // outil rigide -> solide
    double bcWork_ = 0.0;      // platines/grips (PRESCRIBED) -> solide
    double confWork_ = 0.0;    // pression de confinement/bore -> solide
    // E2 : moniteur d'energie runtime (opt-in budgetAbortPct, 0 = off)
    double eAbortPct_ = -1.0;  // -1 = pas encore lu dans la config
    double eAbortMin_ = 0.0;   // plancher absolu [J/m]
    bool eAbort_ = false;
    void checkEnergyAbort();
    double biasW_ = 0.0;       // correction leapfrog EXACTE : les compteurs
                               // lisent v- ; le theoreme discret veut
                               // (v- + v+)/2 -> ecart = f_tot^2 dt^2 / 2m par
                               // noeud et par pas, accumule ici (>= 0)
    double keInit_ = -1.0;     // KE au premier pas (< 0 = pas encore prise)
    Platen plTop_, plBot_;                 // displacement control
    double platenV_ = 0.0;                 // CLOSURE rate; each platen v/2
    LoadArc arcTop_, arcBot_;              // load control
    double loadRate_ = 0.0;                // [Pa/s] on the loaded arcs
    double brazP_ = 0.0;                   // current arc pressure [Pa]
    double sigmaT_ = 0.0, sigmaTpeak_ = 0.0;
    // The BTS is the peak of the load curve at FAILURE. Past that the platens
    // keep closing and crush the two halves against each other, which drives
    // the force far above it (measured: 1113 MPa of nominal sigma_t against a
    // 10 MPa material). The peak is therefore LOCKED at the first genuine
    // post-cracking load drop.
    bool peakLocked_ = false;
    // End-of-test marker. peakLocked_ already latches the strength at the
    // first genuine post-cracking load drop; past that instant the platens are
    // no longer running an indirect tension test, they are crushing the two
    // halves against each other, and the load rises again for a purely
    // kinematic reason. brazilianStopAfterPeak = true ends the run there (plus
    // brazilianStopDelay seconds, so the drop itself is on the curve);
    // whatever the setting, history.csv carries the peakLocked flag so the
    // post-processing can truncate the curve at the same instant.
    bool brazStop_ = false;
    double brazStopDelay_ = 0.0;
    // ---- metrologie de la compression uniaxiale/triaxiale (platens) -------
    // Le meme marqueur de fin d'essai que le bresilien, transpose a la
    // compression : ucsStopAfterPeak = true arrete le run ucsStopDelay
    // secondes apres la premiere chute franche (30 % du pic) post-fissuration.
    // Par defaut false : la duree reste T, comportement d'origine inchange.
    // pullDelay : retard du demarrage du chargement axial (s). 0 = inchange.
    double pullDelay_ = 0.0;
    bool ucsStop_ = false;
    double ucsStopDelay_ = 0.0;
    bool peakLockedU_ = false;
    double tLockedU_ = -1.0;
    // Deformation axiale. TROIS mesures, portees par history.csv :
    //   epsPlaten : fermeture des plateaux / ecartement initial. C'est la
    //               deformation MACHINE : elle contient la compliance du
    //               contact penalite plateau-eprouvette.
    //   epsSpec   : (uy moyen de la rangee y = 0) - (uy moyen de y = H), sur H.
    //   epsGauge  : extensometre interieur entre deux bandes de nœuds situees
    //               a gaugeLoFrac et gaugeHiFrac de la hauteur — la mesure
    //               affranchie a la fois du contact et des effets de bord.
    double gap0_ = 0.0;                    // ecartement initial des plateaux
    double gLoFrac_ = 0.25, gHiFrac_ = 0.75;
    std::vector<int> gLoNodes_, gHiNodes_, topNodes_, botNodes_;
    double gLoY_ = 0.0, gHiY_ = 0.0;       // y0 moyen effectif des deux bandes
    void setupStrainGauge();
    void gaugeStrain(double& epsPlaten, double& epsSpec,
                     double& epsGauge) const;
    double tLocked_ = -1.0;                // time peakLocked_ turned true
    double tPeak_ = 0.0;
    long nBrokenAtPeak_ = 0;
    // elastic gauge at the disc centre, latched on the last step BEFORE the
    // first joint breaks: the parameter-free check of the whole platen+disc
    // setup against the classical solution (sxx = +2P/(pi D t), syy = -3 sxx)
    double gXX_ = 0.0, gYY_ = 0.0, gSigT_ = 0.0;
    // Second, EARLY gauge. The gauge above reads the last step before the
    // first joint breaks, which on an adaptive run is already deep in the
    // diffuse-damage regime: the core has cracked micro-mechanically, its
    // secant stiffness has dropped, and the ratio it reports mixes the platen
    // setup (what the gauge is meant to check) with the constitutive softening
    // (what it is not). This one latches the FIRST time the nominal sigma_t
    // crosses elasticGaugeFrac * ft while the disc is still fully intact, i.e.
    // in the linear range where sigma_xx = +2P/(pi D t) actually holds.
    // It is ACCUMULATED over a band of load rather than read at a single
    // step: the platen contact is a penalty spring and the load it transmits
    // rings at the contact frequency, so any single-step ratio carries +-12 %
    // of that ringing (measured on bd_yan_adaptatif). Averaging sigma_xx and
    // sigma_t separately over the band and dividing the two means is the same
    // estimator with the ringing integrated out.
    double eSumXX_ = 0.0, eSumYY_ = 0.0, eSumSig_ = 0.0;
    long eN_ = 0;
    double eGaugeLo_ = 0.3, eGaugeHi_ = 0.8;   // band, in units of ft
    double eDfrac_ = 0.0;                      // damage at the end of the band

    // ---- SHPB assembly (scenario = shpb, geometry = shpb) ------------------
    // Three DISJOINT bodies meshed in one pass — incident bar, brazilian
    // disc, transmission bar — that talk to each other only through the
    // existing node-to-edge general contact (friction contactMu = the k_c of
    // the article). Bar and rock properties come from two mineral PHASES
    // (phases = bar rock), each body being stamped as one "grain", so the
    // per-phase element tables and the per-joint property assignment work
    // unchanged. Nothing here runs unless geometry = shpb.
    bool shpb_ = false;
    double shpbLIB_ = 2.0, shpbLTB_ = 1.5, shpbD_ = 0.05;   // [m]
    double shpbDrock_ = 0.05, shpbGap_ = 0.0;               // [m]
    double hBar_ = 5e-3, hDisc_ = 7.5e-4;                   // mesh sizes
    bool shpbNoDisc_ = false;              // bar-only wave verification
    int shpbPulse_ = 0;                    // 0 = half-sine, 1 = trapezoid
    double shpbV0_ = 5.2, shpbTau_ = 2.2e-4, shpbPlateau_ = 0.5;
    double shpbM1_ = 0.0, shpbM2_ = 0.0;   // monitor abscissae [m]
    double shpbGaugeW_ = 0.0;              // gauge half-length [m]
    std::vector<int> monEl1_, monEl2_;     // gauge element lists
    std::vector<double> monA1_, monA2_;    // and their areas
    double monArea1_ = 0.0, monArea2_ = 0.0;
    double epsM1_ = 0.0, epsM2_ = 0.0;     // axial strain at the monitors
    double vDrive_ = 0.0;                  // pulse velocity this step
    double xEndAbs_ = 0.0;                 // abscissa of the viscous end
    // eq. 21 of the article reads sigma = -2 rho cp vn, i.e. TWICE the
    // classical Lysmer dashpot rho c v that rockim implements. Keep 1.0
    // (impedance-matched, zero reflection) by default and let the config ask
    // for the literal 2.0 — the two are compared in the report.
    double absFac_ = 1.0;
    // contact-grid controls: a 3.55 m x 0.05 m assembly meshed at 0.75 mm
    // makes the default cell (2 hmin) allocate millions of buckets per step
    double gcCell_ = 0.0;                  // 0 = 2 hmin (unchanged)
    bool gcBoxMesh_ = false;               // AABB of the mesh instead of W/H
    double gcXmin_ = -1e300, gcXmax_ = 1e300;  // active-face x window

    // confining pressure (0 = off)
    double confP_ = 0.0, confRamp_ = 0.0, confL0_ = 0.0;
    std::vector<BEdge> confEdges_;         // ORIGINAL exterior faces only

    // ---- COUPLAGE HYDRO-MECANIQUE (spec 004) — hydro = off | on ----------
    // Modele d'AbuAisha, Eaton, Priest & Wong (JPSE 154, 2017, 100-113), le
    // module HF du code Y-Geo. Hypothese centrale, qu'ils enoncent en toutes
    // lettres (leur §2.2) : « the fluid pressure inside the cavities/fractures
    // is assumed CONSTANT, i.e. no flow due to hydraulic gradients, i.e.
    // similar to an INVISCID flow ». Pas de loi cubique, pas de gradient, pas
    // de leak-off, pas de pas de temps hydraulique propre.
    //
    // C'est exactement ce qui manque a confineFaces = bore, dont le commentaire
    // avoue le defaut : « faces born from cracking receive nothing ». Ici, au
    // contraire, LA PRESSION SUIT LA FISSURE.
    //
    // Trois grandeurs, et une seule variable d'etat :
    //   * la FRONTIERE MOUILLEE — les faces reliees a la source par un chemin
    //     de joints rompus (recherche de composante connexe, leur module 2) ;
    //   * le VOLUME de cavite par le theoreme de Green, leur eq. 4 :
    //         V = 1/2 ∮ x dy - y dx
    //     Discretisee, c'est la formule du lacet SOMMEE FACE PAR FACE, donc
    //     INDEPENDANTE DE L'ORDRE DE PARCOURS : aucun contour a assembler ;
    //   * la PRESSION par compressibilite lineaire, leur eq. 6 :
    //         p = p0 + Kf log( m / (V rho_f0) )
    //     ou m, la masse injectee, est LA variable d'etat.
    // La pression charge ensuite les faces mouillees en force suiveuse.
    //
    // (!) LA COQUILLE DE LEUR EQ. 7, ET CE QU'ELLE A COUTE. Elle s'ecrit
    //     F_p12 = -(p/2) [y2-y1 ; x2-x1]
    // et ce vecteur n'est PAS orthogonal au segment (produit scalaire avec
    // (dx, dy) = 2 dx dy). Le premier lecteur en a conclu « coquille de
    // rotation OU de signe », a voulu appliquer « la forme physique », et a
    // supprime le moins de tete. C'etait l'inverse : LE MOINS ETAIT JUSTE, la
    // coquille est dans la SECONDE COMPOSANTE seule. La forme correcte est
    //     F_p12 = -(p/2) [y2-y1 ; -(x2-x1)] = -(p/2) L n,   n = (dy,-dx)/L
    // soit exactement ce que fait confiningForces(). Deux elements de
    // l'article le tranchent : (a) leur section 3.1 pose « negative sign to
    // compressive stresses », donc p > 0 impose sigma = -p I et la traction
    // t = sigma.n = -p n ; (b) Y-Geo triangule en CCW (Munjiza 2004 ;
    // Mahabadi 2012 ; Lisjak 2014a), donc leur (1 -> 2) est le parcours CCW
    // et (dy, -dx) est bien sortant, comme ici.
    //
    // Le signe inverse a vecu du 19 au 20/08 : le fluide SERRAIT la cavite,
    // le forage faisait un breakout aligne sur sigma'_h et rompait a 6,6 MPa
    // au lieu de 12. Corrige le 2026-08-20 dans hydroForces().
    //
    // LA LECON DE CONTROLE. Le garde-fou annonce ici — H3, le pont Parker —
    // EXISTAIT et a ete passe : il n'a rien vu, parce que parker_compare.py
    // mesurait max(y) - min(y), une valeur ABSOLUE, aveugle au signe. Un
    // controle de signe qui passe par une norme ne controle pas le signe.
    // (Au passage : le label « V3 » qui figurait ici est mort. La spec 004 ne
    // definit plus que H1-H5 pour le fluide et F1-F8 pour les figures ; V1-V8
    // ne survit que dans son tableau d'effort.)
    //
    // UNE SEULE CAVITE pour l'instant (decision F. Uzquiano du 19/08) : leur
    // modele suppose une source unique. Les structures sont dimensionnees pour
    // en accueillir plusieurs (fusion/scission) sans reecriture.
    bool hydroOn_ = false;
    double fluidK_ = 2.2e9;                // K_f [Pa]
    double fluidRho0_ = 1000.0;            // rho_f0 [kg/m3]
    double hydroP0_ = 0.0;                 // p_0 [Pa]
    // ESSAI 2 : seuil d endommagement de conduction du fluide.
    // 1.0 (defaut) = seules les interfaces rompues conduisent.
    double wetDmin_ = 1.0;                 // cle hydroWetDamage
    bool hydroRateMode_ = true;            // rate | pressure
    double hydroRate_ = 0.0;               // [m3/s par m d'epaisseur]
    double hydroPimp_ = 0.0;               // pression imposee [Pa]
    double hydroRamp_ = 0.0;               // rampe cosinus de la pompe [s]
    // ESSAI 3 : instant de demarrage de la pompe [s]. 0 (defaut) =
    // demarrage a t = 0, bit-identique. > 0 = protocole de l article
    // (injection apres le pas geostatique et l excavation).
    double hydroStart_ = 0.0;              // cle hydroStart
    double hydroMass_ = 0.0;               // LA variable d'etat [kg/m]
    double hydroP_ = 0.0, hydroVol_ = 0.0, hydroVol0_ = 0.0;
    double hydroWork_ = 0.0;               // poste SEPARE du bilan B4
    long hydroNWet_ = 0;                   // faces mouillees (sortie)
    long wetStamp_ = -1;                   // nBroken_ au dernier mouillage
    double hydroClose_ = 0.0;              // defaut de FERMETURE du lacet de
                                           // la SOURCE, rapporte a la longueur
                                           // MOYENNE D'UNE FACE (et non au
                                           // perimetre : les noeuds FDEM sont
                                           // dedoubles, l'anneau ne peut pas
                                           // fermer au zero machine). Un lacet
                                           // ne mesure une aire que ferme :
                                           // au-dela de 0,05 le volume de la
                                           // cavite source n'a pas de sens.
    bool hydroCloseWarned_ = false;
    double hydroVolCrack_ = 0.0;           // part FISSURES du volume
    std::vector<char> wetJoint_;           // joint mouille (par joint)
    std::vector<int> wetJointIdx_;         // ... et leur LISTE. Le volume
                                           // est recalcule a CHAQUE pas :
                                           // parcourir les 284 124 joints
                                           // pour n'en trouver aucun coutait
                                           // 50x le temps du run (mesure du
                                           // 2026-08-20). On n'itere que sur
                                           // la liste, vide avant fissuration.
    std::vector<BEdge> hydroSrc_;          // faces SOURCE (le forage)
    std::vector<BEdge> wetEdges_;          // frontiere mouillee courante
    void setupHydro();
    void updateWetBoundary();              // connexite + volume de Green
    void hydroForces();

    // ---- contrainte in situ (etude tunnel EDZ, Wang et al. 2024) ----------
    // Tenseur de contrainte initial GLOBAL, en convention rockim TENSION
    // POSITIVE (les cles, elles, sont des pressions positives : insituSh =
    // 5e6 signifie 5 MPa de COMPRESSION horizontale). Ajoute a la contrainte
    // de chaque element dans elementForces() : sigma_total = sigma0 + D:eps.
    bool hasInsitu_ = false;
    Eigen::Matrix2d insituS_ = Eigen::Matrix2d::Zero();
    // ---- excavation : relachement de la traction de paroi -----------------
    bool excRelease_ = false;
    double excStart_ = 0.0, excRamp_ = 0.0;
    std::vector<BEdge> excEdges_;          // faces ORIGINALES de la cavite
    // the confinement is checked at the END OF ITS RAMP, before the axial load
    // builds: later on the core stress is the triaxial state, not the cell
    // pressure, and comparing it to the target would be meaningless
    bool confLatched_ = false;
    double confAchieved_ = 0.0;

    long stepCount_ = 0;
    double work_ = 0.0, peakF_ = 0.0;
    long nBroken_ = 0;
    int nFrag_ = 1;
    double detachedVol_ = 0.0;
    Eigen::Vector2d gripF_{0.0, 0.0};
    double sigmaPeak_ = 0.0;
};

} // namespace rockim
