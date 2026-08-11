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
        double slip[2] = {0.0, 0.0};       // frictional slip per point
        double omax[2] = {0.0, 0.0};       // largest opening ever reached
                                           // (jointSoftening = yan: the omax
                                           // of eq. 17, per integration point)
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
    void buildMesh();                      // dispatch: grid | voronoi
    void buildMeshVoronoi();
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
    void placeTool();
    // adaptive insertion machinery (insertion = adaptive)
    void buildBindingTables();
    void rebindVertex(int v);
    void insertionSweep();
    void activateJoint(int jI, double sig, double tau);
    void setupBrazilianLoad();
    void initPlatens(double xc, double hw);  // geometry + tributary weights
    void setupBoundaries();
    void setupConfinement();
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
    yan::Params yanP_;
    double yanI_ = 1.0;                    // int_0^1 f(D) dD

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
    std::vector<double> hEl_;              // per-element inscribed size (4A/P)

    // contact
    double kp_ = 0.0, muC_ = 0.5, xiC_ = 0.05, vReg_ = 1e-3;
    double kpGC_ = 0.0, xiGC_ = 0.4, gcWork_ = 0.0, relax_ = 1.0, gcRest_ = 0.2;
    std::unordered_map<uint64_t, double> pen0_;

    double damping_ = 0.02, pullV_ = 0.05;
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
    // OpenMP scratch: per-thread force buffers for the joint scatter
    std::vector<std::vector<Eigen::Vector2d>> fTL_;
    std::vector<std::vector<char>> seenTL_;
    std::vector<std::vector<int>> touchedTL_;
    long actStamp_ = -1;                   // nBroken at last rebuild
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
