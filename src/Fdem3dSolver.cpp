// ---------------------------------------------------------------------------
// Fdem3dSolver — 3D Munjiza-style FDEM. Derivations in comments; the model
// mirrors the (verified) 2D FdemSolver, lifted to Kuhn tetrahedra, with the
// same two mesh front-ends (structured grid | Voronoi grains + phases).
// ---------------------------------------------------------------------------
#include "rockim/Fdem3dSolver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <stdexcept>

#include "rockim/PotentialContact.hpp"
#include "rockim/RandomField.hpp"
#include "rockim/Tessellation3.hpp"
#include "rockim/YangDif.hpp"
#include "rockim/VtkWriter.hpp"

#include <chrono>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace rockim {

namespace {
// ROCKIM_PROF=1 — per-phase step profile, the 2D FdemSolver pattern in 3D.
struct F3Prof {
    double tEl = 0, tIn = 0, tJt = 0, tGc = 0, tTc = 0;
    long n = 0;
    bool on = std::getenv("ROCKIM_PROF") != nullptr;
    ~F3Prof() {
        if (!on || n == 0) return;
        std::fprintf(stderr,
                     "[prof3d] per step (ms): elem %.2f insert %.2f joint %.2f "
                     "gcontact %.2f tool %.2f  (%ld steps)\n",
                     1e3 * tEl / n, 1e3 * tIn / n, 1e3 * tJt / n,
                     1e3 * tGc / n, 1e3 * tTc / n, n);
    }
} f3Prof;
double f3now() {
    return std::chrono::duration<double>(
               std::chrono::steady_clock::now().time_since_epoch()).count();
}
} // namespace

Fdem3dSolver::Fdem3dSolver(const Config& cfg, std::string outDir)
    : cfg_(cfg), out_(std::move(outDir)) {}

void Fdem3dSolver::init() {
    // ---- viscosite de Yan (eq. 6) et DIF de Yang (eq. 2-3) --------------
    // Portes du 2D le 2026-08-19. La garde de parite qui les refusait ici
    // est levee. Opt-in strict : cles absentes = branches non calculees.
    viscIns_ = cfg_.geti("viscousInInsertion", 1) != 0;
    {
        std::string sd = cfg_.gets("strainRateDIF", "off");
        if (sd != "off" && sd != "yang" && sd != "yang-fig2")
            throw std::runtime_error("strainRateDIF must be off | yang | "
                                     "yang-fig2 (Yang et al., IJRMMS 191, "
                                     "2025 : yang = leur eq. 3 litterale, "
                                     "exposant 0,07 et loi discontinue ; "
                                     "yang-fig2 = exposant 0,1707 deduit de "
                                     "leur figure 2b, loi continue)");
        difOn_ = sd != "off";
        difExpT_ = sd == "yang-fig2" ? 0.1707 : 0.07;
    }
    if (difOn_ && cfg_.gets("insertion", "intrinsic") != "adaptive")
        throw std::runtime_error("strainRateDIF = yang exige insertion = "
                                 "adaptive : le DIF est fige a l instant "
                                 "d insertion du joint, et le schema "
                                 "intrinseque n en a pas");
    // jointShearEnvelope / meanTensionCapFactor : voir le header.
    {
        std::string se = cfg_.gets("jointShearEnvelope", "yan");
        if (se != "yan" && se != "yang")
            throw std::runtime_error("jointShearEnvelope must be yan | yang "
                                     "(yan = son eq. 8, le frottement tombe a "
                                     "zero des la traction ; yang = son eq. 1, "
                                     "il decroit jusqu au cut-off en ft)");
        yangEnv_ = se == "yang";
    }
    // ---- WP1 : pulverisation (Yang et al. 2026, IJRMMS 206, eq. 3-4) ----
    // Degradation de raideur des tetraedres : sigma = Cd (1 - D) sigma_b
    // au-dela du seuil, D lineaire en deplacement effectif delta_m = h_e *
    // eps_vm (Camanho), irreversible, plafonne a Dmax. delta_m0 / delta_mf
    // en METRES : la calibration de l article (0,014 / 0,4, elements de
    // 1 mm) se lit 1.4e-5 / 4.0e-4. ADDITION (principe VIII) : le crushCap
    // reste tel quel — le deck granite le neutralise (crushCap = 1e12).
    {
        std::string bd = cfg_.gets("bulkDamage", "off");
        if (bd != "off" && bd != "yang")
            throw std::runtime_error("bulkDamage must be off | yang "
                                     "(Yang et al. 2026, eq. 3-4 : "
                                     "pulverisation par degradation de "
                                     "raideur des elements volumiques)");
        bdOn_ = bd == "yang";
    }
    if (bdOn_) {
        bdD0_   = cfg_.getd("bulkDamageDelta0", 1.4e-5);
        bdDf_   = cfg_.getd("bulkDamageDeltaF", 4.0e-4);
        bdDmax_ = cfg_.getd("bulkDamageDmax", 0.9);
        bdCd_   = cfg_.getd("bulkDamageCd", 1.0);
        if (!(bdD0_ > 0.0) || !(bdDf_ > bdD0_))
            throw std::runtime_error("bulkDamage : il faut 0 < "
                                     "bulkDamageDelta0 < bulkDamageDeltaF "
                                     "[m] (deplacements effectifs h_e*eps)");
        if (!(bdDmax_ > 0.0) || bdDmax_ > 1.0 || !(bdCd_ > 0.0))
            throw std::runtime_error("bulkDamage : bulkDamageDmax dans "
                                     "]0, 1] et bulkDamageCd > 0");
        if (cfg_.gets("law", "elastic") != "elastic")
            throw std::runtime_error("bulkDamage = yang exige law = elastic "
                                     "(la degradation s applique a la "
                                     "branche co-rotationnelle elastique, "
                                     "pas a une loi MatLaw)");
    }
    mtCap_ = cfg_.getd("meanTensionCapFactor", 3.0);
    srTau_ = cfg_.getd("strainRateTau", 1.0e-6);
    if (difOn_ && !(srTau_ > 0.0))
        throw std::runtime_error("strainRateTau must be > 0 [s]");
    mat_ = Material::from(cfg_);
    phases_ = PhaseSet::from(cfg_);

    std::string sc = cfg_.gets("scenario", "percussion");
    if      (sc == "percussion") scen_ = Scenario::PERCUSSION;
    else if (sc == "shear")      scen_ = Scenario::SHEAR;
    else if (sc == "tension")    scen_ = Scenario::TENSION;
    else throw std::runtime_error("fdem3d scenario must be percussion | shear | tension");

    W_ = cfg_.getd("W", 0.08);
    D_ = cfg_.getd("D", 0.08);
    H_ = cfg_.getd("H", 0.06);
    nx_ = cfg_.geti("nx", 20);
    ny_ = cfg_.geti("ny", 20);
    nz_ = cfg_.geti("nz", 15);
    T_ = cfg_.getd("T", 2e-4);
    damping_ = cfg_.getd("dampingLocal", scen_ == Scenario::TENSION ? 0.7 : 0.05);

    buildMesh();

    // voronoi contact cell: 2 x median element size, floored so the ±1-cell
    // sweep still covers the deepest allowed penetration (0.6 x the largest
    // element) — hmin (the thinnest sliver) would explode the cell count
    if (voronoi_) {
        std::vector<double> h = hEl_;
        std::nth_element(h.begin(), h.begin() + h.size() / 2, h.end());
        double hMed = h[h.size() / 2];
        double hMax = *std::max_element(hEl_.begin(), hEl_.end());
        cellV_ = std::max(2.0 * hMed, 0.61 * hMax);
    }

    // ---- per-phase element tables (elementForces hot loop) ------------------
    lamP_.clear(); mu2P_.clear(); crushCapP_.clear(); ftP_.clear(); rhoP_.clear();
    for (const Material& m : phases_.mat) {
        lamP_.push_back(m.E * m.nu / ((1 + m.nu) * (1 - 2 * m.nu)));
        mu2P_.push_back(m.E / (1 + m.nu));                 // 2G
        crushCapP_.push_back(cfg_.getd("crushCap", 8.0 * m.cohesion));
        ftP_.push_back(m.ft);
        rhoP_.push_back(m.rho);
    }

    // ---- cohesive joint law (PER JOINT, as in 2D) ---------------------------
    // ---- optional bulk constitutive law -------------------------------------
    if (cfg_.has("law")) {
        if (phases_.n() > 1)
            throw std::runtime_error("'law' (bulk constitutive law) is a SINGLE "
                "material model: it cannot be combined with mineral 'phases'.");
        double lcMax = 0.0;
        for (double h : hEl_) lcMax = std::max(lcMax, h);
        // ---- SURCHARGES DE VOLUME (2026-08-19) ---------------------------
        // Le bloc materiau sert A LA FOIS a la loi de volume et aux joints,
        // et les deux ne veulent pas les memes nombres : dans une
        // architecture « a la Ye », ft et cohesion sont des MICRO-parametres
        // de joint (25 et 60 MPa sur indent3d_grad) tandis que la roche a ses
        // valeurs macro (11,4 et 42,8 MPa pour le granite de Kuru). Jusqu'ici
        // seule la loi mc avait ses surcharges (mcCohesion, mcFrictionDeg) ;
        // dpr, saksala et saksala2011 heritaient donc en silence des valeurs
        // de JOINT — ce qui rend impossible de changer de loi a joints
        // constants, c'est-a-dire l'experience meme qu'on veut faire.
        // Absentes, ces cles ne changent rien : comportement bit-identique.
        Material mBulk = mat_;
        mBulk.ft       = cfg_.getd("bulkFt", mBulk.ft);
        mBulk.cohesion = cfg_.getd("bulkCohesion", mBulk.cohesion);
        mBulk.phiDeg   = cfg_.getd("bulkFrictionDeg", mBulk.phiDeg);
        mBulk.Gf       = cfg_.getd("bulkGf", mBulk.Gf);
        if (cfg_.has("bulkFt") || cfg_.has("bulkCohesion")
            || cfg_.has("bulkFrictionDeg") || cfg_.has("bulkGf")) {
            PhaseSet::validate(mBulk, "loi de volume (cles bulk*)");
            std::cout << "[MATLAW] loi de volume dissociee des joints : ft = "
                      << mBulk.ft << " Pa, cohesion = " << mBulk.cohesion
                      << " Pa, phi = " << mBulk.phiDeg << " deg, Gf = "
                      << mBulk.Gf << " J/m2 (les joints gardent ft = "
                      << mat_.ft << ", c = " << mat_.cohesion << ")\n";
        }
        law_ = MatLaw::make(cfg_.gets("law", "elastic"), mBulk, cfg_, lcMax);
        // C2 (audit 2026-08-11, corrige 2026-08-15) : le centroide INITIAL de
        // chaque element doit etre renseigne AVANT le premier appel a la loi.
        // dpdfh en tire ses trois seuils de Weibull par hash spatial 64 bits
        // (exactement comme la VUMAT hashe coordMp) : sans lui, tous les
        // elements partagent le meme hash — heterogeneite DFH morte.
        for (auto& e : el_)
            e.st.x0 = 0.25 * (X0_[e.n[0]] + X0_[e.n[1]]
                              + X0_[e.n[2]] + X0_[e.n[3]]);
        std::cout << "[FDEM3D] bulk law = " << law_->name() << ", "
                  << el_.size() << " tets, lc max = " << lcMax << " m\n";
    }

    // insertion = adaptive : extrinsic scheme of Yan et al. (IJRMMS 169,
    // 2023), ported from the 2D solver. Read BEFORE assignJointProps, which
    // selects the (softer) activation penalty in that mode.
    {
        std::string ins = cfg_.gets("insertion", "intrinsic");
        if (ins != "intrinsic" && ins != "adaptive")
            throw std::runtime_error("insertion must be intrinsic | adaptive "
                                     "(got '" + ins + "')");
        adaptive_ = ins == "adaptive";
    }
    // jointSoftening = linear (default, unchanged) | yan — the exponential
    // reduction factor f(D) of Yan et al. eq. 11, ported from the 2D solver.
    // Read BEFORE assignJointProps(): it sets the critical opening/slip from
    // the fracture energies through I = int_0^1 f(D) dD.
    {
        std::string js = cfg_.gets("jointSoftening", "linear");
        if (js != "linear" && js != "yan")
            throw std::runtime_error("jointSoftening must be linear | yan "
                                     "(got '" + js + "')");
        yanSoft_ = js == "yan";
    }
    if (yanSoft_) {
        yanP_.a = cfg_.getd("yanA", 0.63);
        yanP_.b = cfg_.getd("yanB", 1.8);
        yanP_.c = cfg_.getd("yanC", 6.0);
        yanFricScaled_ = cfg_.geti("jointFrictionScaled", 0) != 0;
        yanI_ = yan::integralFD(yanP_, cfg_.geti("yanQuadN", 4096));
        if (!(yanI_ > 1e-6))
            throw std::runtime_error("jointSoftening = yan: int f(D) dD is "
                                     "not positive — check yanA/yanB/yanC");
        std::cout << "[FDEM3D] joint softening: Yan et al. f(D), a = "
                  << yanP_.a << ", b = " << yanP_.b << ", c = " << yanP_.c
                  << ", int f(D) dD = " << yanI_ << "\n";
    }
    // jointShearUnload = plastic (defaut, inchange) | origin — eq. 18 de Yan
    // et al., miroir exact du 2D (voir FdemSolver.cpp).
    {
        std::string ju = cfg_.gets("jointShearUnload", "plastic");
        if (ju != "plastic" && ju != "origin")
            throw std::runtime_error("jointShearUnload must be plastic | origin "
                                     "(got '" + ju + "')");
        shearOrigin_ = ju == "origin";
    }
    if (shearOrigin_) {
        std::cout << "[FDEM3D] shear unloading: origin secant (Yan eq. 18)\n";
        if (!yanFricScaled_)
            std::cout << "[FDEM3D] WARNING: jointShearUnload = origin with "
                         "jointFrictionScaled = 0 — the Coulomb term rides the "
                         "origin secant, so frictional sliding is REVERSIBLE "
                         "(no hysteresis loop). Literal article form is "
                         "origin + jointFrictionScaled = 1.\n";
    }
    assignJointProps();
    if (adaptive_) {
        for (auto& J : jt_) J.bonded = true;
        buildBindingTables();
        std::cout << "[FDEM3D] adaptive insertion: " << jt_.size()
                  << " bonded faces, activation penalty "
                  << cfg_.getd("insertionPenaltyFactor", 4.0) << " E/h\n";
    }
    xiJ_ = cfg_.getd("jointXi", 0.05);

    kp_ = phases_.maxE() * hmin_;                      // tool contact [N/m]
    // ---- A1 : loi de contact de l'outil, miroir du 2D --------------------
    {
        std::string tc = cfg_.gets("toolContact", "penalty");
        if (tc != "penalty" && tc != "signorini")
            throw std::runtime_error("toolContact must be penalty | signorini");
        toolSig_ = (tc == "signorini");
        toolSigRelax_ = cfg_.getd("toolSignoriniRelax", 0.0);
        if (!(toolSigRelax_ >= 0.0 && toolSigRelax_ <= 1.0))
            throw std::runtime_error("toolSignoriniRelax must be in [0, 1]");
        if (toolSig_)
            std::cout << "[FDEM3D] contact outil : SIGNORINI en vitesse "
                         "(CD-Lagrange, Delassus diagonal par noeud) — la "
                         "penalite kp = " << kp_ << " N/m ne s'applique plus "
                         "et SORT du budget du pas de temps ; rattrapage "
                         "d'interpenetration = " << toolSigRelax_ << "\n";
    }
    muC_ = cfg_.getd("contactMu", 0.5);
    xiC_ = cfg_.getd("contactXi", 0.05);
    vReg_ = cfg_.getd("contactVreg", 1e-3);
    kpGC_ = cfg_.getd("gcPenaltyFactor", 0.01) * phases_.maxE() * hmin_;
    xiGC_ = cfg_.getd("gcXi", 0.8);
    gcRest_ = cfg_.getd("gcRestitution", 0.2);
    // contact = penalty (defaut, inchange) | potential — A3 phase 2 : le
    // contact par POTENTIEL de Munjiza en tet-tet, miroir exact du 2D.
    {
        std::string cm = cfg_.gets("contact", "penalty");
        if (cm != "penalty" && cm != "potential")
            throw std::runtime_error("contact must be penalty | potential "
                                     "(got '" + cm + "')");
        contactPot_ = cm == "potential";
    }
    if (contactPot_) {
        potP_ = cfg_.getd("potPenaltyFactor", 1.0) * phases_.maxE();
        jcAdaptive_ = cfg_.gets("jointContactPenalty", "fixed") == "adaptive";
        if (jcAdaptive_)
            std::cout << "[FDEM3D] jointContactPenalty = adaptive : k- = "
                         "k+(D) = (1-D) pj (EPFL arXiv:2511.14323 sec. 4)\n";
        potKt_ = cfg_.getd("potTangentFactor", 1.0) * phases_.maxE() * hmin_;
        std::cout << "[FDEM3D] contact: Munjiza potential (eq. 2-5, tet-tet)"
                     ", p = " << potP_ << " Pa, kt = " << potKt_ << " N/m\n";
    }
    // gcActivation = full (defaut, inchange) | adaptive — miroir exact du 2D
    // (FdemSolver.hpp) : seules les faces qui PEUVENT toucher sont balayees.
    {
        std::string ga = cfg_.gets("gcActivation", "full");
        if (ga != "full" && ga != "adaptive")
            throw std::runtime_error("gcActivation must be full | adaptive "
                                     "(got '" + ga + "')");
        gcAdaptive_ = ga == "adaptive";
    }
    if (gcAdaptive_) {
        gcActMargin_ = cfg_.getd("gcActMargin", 2.0);
        gcActEvery_ = cfg_.geti("gcActEvery", 64);
        if (!(gcActMargin_ > 0.0))
            throw std::runtime_error("gcActMargin must be > 0 [cells]");
        if (gcActEvery_ < 1)
            throw std::runtime_error("gcActEvery must be >= 1 [steps]");
        std::cout << "[FDEM3D] contact activation: adaptive (Fukuda) — margin "
                  << gcActMargin_ << " cells, sweep <= every "
                  << gcActEvery_ << " steps\n";
    }

    placeTool();
    setupBoundaries();
    setupConfinement();
    // ---- viscosite de Yan : mu par element ------------------------------
    // mu = xi h sqrt(E rho) : xi est le taux d amortissement a l echelle de
    // la MAILLE. xi = 2 est exactement l amortissement critique de Munjiza
    // mu = 2 h sqrt(E rho) — applique aux chiffres de la Table 1 de Yan
    // (h = 0,75 mm, E = 15 GPa, rho = 1704) il redonne 7583 Pa.s contre les
    // 7600 publies, a 0,2 %. Leur viscosite EST le critique de Munjiza, ce
    // que leur article ne dit pas (il renvoie a Tatone & Grasselli sans
    // formule ni etude de sensibilite).
    {
        double muLit = cfg_.getd("bulkViscosity", 0.0);
        double xiV   = cfg_.getd("bulkViscosityXi", 0.0);
        bool graded  = cfg_.geti("bulkViscosityGraded", 0) != 0;
        if (muLit < 0.0 || xiV < 0.0)
            throw std::runtime_error("bulkViscosity / bulkViscosityXi must "
                                     "be >= 0");
        if (muLit > 0.0 && xiV > 0.0)
            throw std::runtime_error("bulkViscosity et bulkViscosityXi sont "
                                     "exclusives : la seconde CALCULE la "
                                     "premiere a partir du maillage");
        viscOn_ = muLit > 0.0 || xiV > 0.0;
        if (viscOn_) {
            muEl_.assign(el_.size(), muLit);
            if (xiV > 0.0) {
                std::vector<double> srt(el_.size());
                for (std::size_t eI = 0; eI < el_.size(); ++eI) {
                    double E = phases_.mat[el_[eI].phase].E;
                    double rho = rhoP_[el_[eI].phase];
                    srt[eI] = xiV * hEl_[eI] * std::sqrt(E * rho);
                }
                if (graded) muEl_ = srt;
                else {
                    std::vector<double> s2 = srt;
                    std::sort(s2.begin(), s2.end());
                    double med = s2.empty() ? 0.0 : s2[s2.size() / 2];
                    muEl_.assign(el_.size(), med);
                }
            }
            double lo = 1e300, hi = 0.0;
            for (double m : muEl_) { lo = std::min(lo, m); hi = std::max(hi, m); }
            std::cout << "[FDEM3D] viscosite newtonienne (Yan eq. 6) : mu "
                      << (graded ? "GRADUE par element, " : "GLOBAL, ")
                      << lo << " a " << hi << " Pa.s";
            if (xiV > 0.0)
                std::cout << " (bulkViscosityXi = " << xiV << ", soit "
                          << 0.5 * xiV << " x le critique de Munjiza ; Yan "
                             "et al. Table 1 = 2,00 = critique)";
            std::cout << "\n";
            if (!graded && xiV > 0.0)
                std::cout << "[FDEM3D]   NOTE : mu GLOBAL sur un maillage "
                             "gradue paie le pas de temps du plus FIN tetra "
                             "partout. bulkViscosityGraded = 1 cale mu sur "
                             "chaque element et rend la borne diffusive "
                             "independante de la gradation.\n";
            if (!viscIns_)
                std::cout << "[FDEM3D]   viscousInInsertion = 0 : la "
                             "contrainte visqueuse n entre PAS dans sigG, "
                             "donc ni dans le critere d insertion ni dans la "
                             "jauge de confinement. Ecart assume avec Yan, "
                             "dont l eq. 6 EST la contrainte lue par son "
                             "eq. 7.\n";
        }
    }
    computeStableDt();
    relax_ = std::exp(-dt_ / cfg_.getd("gcBirthTau", 1e-6));
    srRelax_ = srTau_ > 0.0 ? std::exp(-dt_ / srTau_) : 0.0;
    if (difOn_) {
        std::cout << "[FDEM3D] strainRateDIF = "
                  << (difExpT_ > 0.1 ? "yang-fig2" : "yang")
                  << " (Yang et al., IJRMMS 191, 2025, eq. 2-3) : ft et Gf "
                     "recoivent DIF_traction, cohesion et GfII recoivent "
                     "DIF_compression, l angle de frottement est inchange. "
                     "Facteurs FIGES a l insertion.\n"
                  << "[FDEM3D]   taux = principale max de sym(Fdot F^-1), "
                     "filtre a strainRateTau = " << srTau_ << " s ("
                  << (srTau_ / dt_) << " pas)\n";
        if (difExpT_ < 0.1)
            std::cout << "[FDEM3D]   AVERTISSEMENT : exposant litteral 0,07 "
                         "— la loi saute de 12 % en 5e-6 /s et de 22 % en "
                         "1e2 /s. En insertion extrinseque ce saut est un "
                         "ATTRACTEUR : la population inseree s empile juste "
                         "sous 1e2 (mesure 2D du 2026-08-18). "
                         "strainRateDIF = yang-fig2 l evite.\n";
        if (viscOn_ && viscIns_)
            std::cout << "[FDEM3D]   AVERTISSEMENT : la viscosite est armee "
                         "et entre dans sigG. Le taux de deformation agit "
                         "alors DEUX FOIS : il gonfle la contrainte d essai "
                         "ET le seuil. Fidele a Yan, mais ce n est le modele "
                         "de personne. viscousInInsertion = 0 les separe.\n";
    }

    if (scen_ == Scenario::TENSION) pullV_ = cfg_.getd("pullV", 0.05);
    pullDelay_ = cfg_.getd("pullDelay", 0.0);
    gripFree_ = cfg_.getb("gripLateralFree", false);
    pullRamp_ = cfg_.getd("pullRamp", 0.0);
    fragId_.assign(el_.size(), 0);
    toolKE0_ = tool_.ke();

    // Un maillage STRUCTURE n'est pas un choix esthetique : il biaise les
    // trajets de fissures ET peut diverger en phase debris (2D mesure
    // 2026-08-06 ; 3D confirme 2026-08-11 — grille de Kuhn intrinseque en
    // cascade, gcWork +584 J pour 16 J incidents). Les scenarios d'impact et
    // de coupe doivent partir du maillage desordonne.
    if (!voronoi_ && scen_ != Scenario::TENSION)
        std::cout << "[FDEM3D] WARNING: mesh = grid + scenario de "
                     "fissuration — un maillage structure biaise les trajets "
                     "et peut diverger en phase debris (FICHE 2026-08-06/11). "
                     "Utiliser mesh = voronoi + grainSeeding = random pour "
                     "tout resultat ou la casse compte.\n";
    std::cout << "[FDEM3D] " << el_.size() << " tets, " << jt_.size()
              << " joints, " << X0_.size() << " nodes, dt = " << dt_
              << " s, steps = " << (long)std::ceil(T_ / dt_) << "\n";
    if (voronoi_) {
        long nGB = 0;
        for (const auto& J : jt_) if (J.type > 0) ++nGB;
        std::cout << "[FDEM3D] voronoi: " << nGrains_ << " grains, "
                  << phases_.n() << " phase(s), " << nGB
                  << " grain-boundary joints, hmin = " << hmin_ << " m\n";
    }
}

// ---------------------------------------------------------------------------
// Mesh generation. Two front-ends share one topology builder, as in 2D:
//   mesh = grid    — structured hex grid split into the 6 Kuhn tetrahedra
//                    per cell (identical corner ordering makes the induced
//                    face diagonals compatible across cells), optional
//                    interior-vertex jitter (one implicit "grain").
//   mesh = voronoi — 3D Voronoi grain structure (Tessellation3) with
//                    per-grain mineral phases: the GBM mode.
// Every interior triangular face gets a cohesive joint; boundary faces feed
// the quiet boundaries and the general contact.
// ---------------------------------------------------------------------------
void Fdem3dSolver::buildMesh() {
    std::string mesh = cfg_.gets("mesh", "grid");
    if (mesh != "grid" && mesh != "voronoi" && mesh != "file")
        throw std::runtime_error("mesh must be grid | voronoi | file (got '"
                                 + mesh + "')");
    voronoi_ = mesh == "voronoi";
    if (mesh == "file") {
        // Plusieurs phases sont admises si (et seulement si) le maillage
        // porte des groupes physiques nommes : chaque corps prend la phase
        // homonyme (ou groupPhase.<nom>). La verification se fait dans
        // buildMeshFile, une fois les groupes connus.
        buildMeshFile();
        return;
    }
    if (!voronoi_ && phases_.n() > 1)
        throw std::runtime_error("'phases' declares "
            + std::to_string(phases_.n()) + " minerals but mesh = grid would "
            "silently use only the first: set mesh = voronoi (or drop the "
            "phases key)");
    if (voronoi_) { buildMeshVoronoi(); return; }

    double dx = W_ / nx_, dy = D_ / ny_, dz = H_ / nz_;
    hmin_ = std::min({dx, dy, dz});
    double jit = cfg_.getd("meshJitter", 0.0) * 0.5 * hmin_;
    std::mt19937 rng(cfg_.geti("seed", 12345));
    std::uniform_real_distribution<double> U(-jit, jit);

    int vnx = nx_ + 1, vny = ny_ + 1, vnz = nz_ + 1;
    std::vector<Eigen::Vector3d> vc((std::size_t)vnx * vny * vnz);
    auto vid = [&](int i, int j, int k) { return (k * vny + j) * vnx + i; };
    for (int k = 0; k < vnz; ++k)
        for (int j = 0; j < vny; ++j)
            for (int i = 0; i < vnx; ++i) {
                Eigen::Vector3d p(i * dx, j * dy, k * dz);
                if (jit > 0 && i > 0 && i < nx_ && j > 0 && j < ny_
                    && k > 0 && k < nz_)
                    p += Eigen::Vector3d(U(rng), U(rng), U(rng));
                vc[vid(i, j, k)] = p;
            }

    std::vector<std::array<int, 4>> tets;
    tets.reserve((std::size_t)6 * nx_ * ny_ * nz_);
    const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                             {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (int k = 0; k < nz_; ++k)
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                int c[2][2][2];
                for (int a = 0; a < 2; ++a)
                    for (int b = 0; b < 2; ++b)
                        for (int cph = 0; cph < 2; ++cph)
                            c[a][b][cph] = vid(i + a, j + b, k + cph);
                for (const auto& p : perms) {
                    int s[3] = {0, 0, 0};
                    std::array<int, 4> vv;
                    vv[0] = c[0][0][0];
                    s[p[0]] = 1; vv[1] = c[s[0]][s[1]][s[2]];
                    s[p[1]] = 1; vv[2] = c[s[0]][s[1]][s[2]];
                    vv[3] = c[1][1][1];
                    tets.push_back(vv);
                }
            }
    std::vector<int> tetGrain(tets.size(), 0);         // one implicit grain
    nGrains_ = 1;
    buildFromTets(vc, tets, tetGrain, {0});
}

// ---------------------------------------------------------------------------
// mesh = file — maillage tetraedrique non structure importe (Gmsh MSH 2.2
// ASCII, $Nodes + $Elements de type 4). C'est le maillage "a la Yan et al."
// (IJRMMS 169, 2023, figs 8/9b/11) : simplexes uniformes desordonnes, sans
// structure de grains — la pratique des codes etablis (Gmsh chez Akantu,
// OpenFDEM, la litterature FDEM ; cf. veille P1). Generer par exemple avec
// tools/make_unstructured_mesh.py, puis :
//     mesh = file
//     meshFile = maillage.msh
// Le maillage est TRANSLATE pour que sa boite englobante parte de l'origine,
// et W/D/H sont REDEFINIS depuis cette boite (les frontieres, l'outil et les
// jauges s'appuient sur les plans x/y/z = 0 et W/D/H). Ids de noeuds non
// contigus acceptes ; l'orientation des tets est reparee par buildFromTets.
// Le chemin de dimensionnement non uniforme (h local par element, grille de
// contact creuse) est celui du voronoi : voronoi_ = true, un seul grain.
// ---------------------------------------------------------------------------
void Fdem3dSolver::buildMeshFile() {
    std::string path = cfg_.reqs("meshFile");
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("meshFile: cannot open '" + path + "'");
    std::string line;
    std::map<long, int> id2idx;
    std::vector<Eigen::Vector3d> vpos;
    std::vector<std::array<int, 4>> tets;
    std::vector<long> tetPhys;             // tag physique par tet (0 = aucun)
    std::map<long, std::string> physVol;   // id physique (dim 3) -> nom
    bool sawFormat = false;
    while (std::getline(in, line)) {
        if (line.rfind("$MeshFormat", 0) == 0) {
            std::getline(in, line);
            double ver = std::atof(line.c_str());
            if (ver < 2.0 || ver >= 3.0)
                throw std::runtime_error("meshFile: MSH version "
                    + std::to_string(ver) + " unsupported — export ASCII 2.2 "
                    "(Gmsh: Mesh.MshFileVersion = 2.2)");
            sawFormat = true;
        } else if (line.rfind("$PhysicalNames", 0) == 0) {
            long n = 0;
            in >> n;
            for (long k = 0; k < n; ++k) {
                int dim; long id; std::string nm;
                in >> dim >> id;
                std::getline(in, nm);
                auto q0 = nm.find('"');
                auto q1 = nm.rfind('"');
                if (q0 != std::string::npos && q1 > q0)
                    nm = nm.substr(q0 + 1, q1 - q0 - 1);
                if (dim == 3) physVol[id] = nm;        // surfaces : plus tard
            }
        } else if (line.rfind("$Nodes", 0) == 0) {
            long n = 0;
            in >> n;
            for (long k = 0; k < n; ++k) {
                long id; double x, y, z;
                in >> id >> x >> y >> z;
                id2idx[id] = (int)vpos.size();
                vpos.push_back({x, y, z});
            }
        } else if (line.rfind("$Elements", 0) == 0) {
            long n = 0;
            in >> n;
            for (long k = 0; k < n; ++k) {
                long id; int type, ntags;
                in >> id >> type >> ntags;
                long tag, phys = 0;
                for (int t = 0; t < ntags; ++t) {
                    in >> tag;
                    if (t == 0) phys = tag;            // 1er tag = physique
                }
                int nn = type == 15 ? 1 : type == 1 ? 2 : type == 2 ? 3
                       : type == 4 ? 4 : -1;
                if (nn < 0)
                    throw std::runtime_error("meshFile: element type "
                        + std::to_string(type) + " unsupported (tets only; "
                        "export a pure tetrahedral mesh)");
                std::array<int, 4> vv{};
                for (int q = 0; q < nn; ++q) {
                    long nid; in >> nid;
                    if (nn == 4) {
                        auto it = id2idx.find(nid);
                        if (it == id2idx.end())
                            throw std::runtime_error("meshFile: element "
                                "references unknown node id");
                        vv[q] = it->second;
                    }
                }
                if (nn == 4) {                         // points/lines/tris:
                    tets.push_back(vv);                // boundary — skipped
                    tetPhys.push_back(phys);
                }
            }
        }
    }
    if (!sawFormat || vpos.empty() || tets.empty())
        throw std::runtime_error("meshFile: no tetrahedra found in '" + path
                                 + "' (need ASCII MSH 2.2 with type-4 "
                                 "elements)");
    // translate to the origin and take the box from the bounding box
    Eigen::Vector3d lo = vpos[0], hi = vpos[0];
    for (const auto& p : vpos) { lo = lo.cwiseMin(p); hi = hi.cwiseMax(p); }
    for (auto& p : vpos) p -= lo;
    W_ = hi.x() - lo.x();
    D_ = hi.y() - lo.y();
    H_ = hi.z() - lo.z();
    if (!(W_ > 0 && D_ > 0 && H_ > 0))
        throw std::runtime_error("meshFile: degenerate bounding box");
    // ---- physical groups (V1) : volumes physiques -> groupes/corps ---------
    // Sans $PhysicalNames : un seul groupe, comportement inchange. Avec :
    // un groupe par volume physique, materiau = phase du meme nom (ou
    // groupPhase.<nom>), et buildFromTets ne cree AUCUN joint entre groupes.
    groupName_.clear();
    std::map<long, int> phys2grp;
    if (!physVol.empty()) {
        for (const auto& [pid, nm] : physVol) {
            phys2grp[pid] = (int)groupName_.size();
            groupName_.push_back(nm);
        }
    } else {
        groupName_.push_back("all");
    }
    nGroups_ = (int)groupName_.size();
    if (phases_.n() > 1 && nGroups_ <= 1)
        throw std::runtime_error("'phases' avec mesh = file exige des groupes "
            "physiques nommes ($PhysicalNames, dim 3) : sans groupes le "
            "maillage est un seul corps et une seule phase s'applique");
    // ---- WP2 : paires de corps LIES (groupBond.<A>.<B> = joints) --------
    // L interface conforme entre deux volumes physiques nommes recoit des
    // joints cohesifs (brasage insert/bit) au lieu d etre remise au contact.
    // La cle accepte les deux ordres de noms. Absente : rien ne change.
    gbond_.assign((std::size_t)nGroups_ * nGroups_, 0);
    for (int g1 = 0; g1 < nGroups_; ++g1)
        for (int g2 = g1 + 1; g2 < nGroups_; ++g2) {
            std::string k1 = "groupBond." + groupName_[g1] + "."
                           + groupName_[g2];
            std::string k2 = "groupBond." + groupName_[g2] + "."
                           + groupName_[g1];
            std::string bv = cfg_.gets(k1, cfg_.gets(k2, ""));
            if (bv.empty()) continue;
            if (bv != "joints")
                throw std::runtime_error(k1 + " : seule la valeur 'joints' "
                    "est supportee (interface cohesive conforme)");
            gbond_[(std::size_t)g1 * nGroups_ + g2] = 1;
            gbond_[(std::size_t)g2 * nGroups_ + g1] = 1;
            std::cout << "[FDEM3D] groupBond: " << groupName_[g1] << " <-> "
                      << groupName_[g2] << " — interface JOINTE (joints de "
                      "frontiere, proprietes GBM moyennes x gb*)\n";
        }
    tetGroupTmp_.assign(tets.size(), 0);
    if (!physVol.empty())
        for (std::size_t k = 0; k < tets.size(); ++k) {
            auto it = phys2grp.find(tetPhys[k]);
            if (it == phys2grp.end())
                throw std::runtime_error("meshFile: un tet porte le tag "
                    "physique " + std::to_string(tetPhys[k]) + " sans volume "
                    "physique declare — nommer TOUS les volumes ou aucun");
            tetGroupTmp_[k] = it->second;
        }
    // phase par groupe : groupPhase.<nom> sinon la phase homonyme sinon 0
    std::vector<int> grpPhase(nGroups_, 0);
    for (int g = 0; g < nGroups_; ++g) {
        std::string want = cfg_.gets("groupPhase." + groupName_[g],
                                     groupName_[g]);
        int ph = -1;
        for (int p = 0; p < phases_.n(); ++p)
            if (phases_.name[p] == want) ph = p;
        if (ph < 0) {
            ph = 0;
            if (phases_.n() > 1)
                std::cout << "[FDEM3D] WARNING: groupe '" << groupName_[g]
                          << "' sans phase homonyme ni groupPhase — phase 0 ("
                          << phases_.name[0] << ")\n";
        }
        grpPhase[g] = ph;
    }
    std::vector<int> tetGrain = tetGroupTmp_;          // grain = groupe (VTU)
    nGrains_ = nGroups_;
    voronoi_ = true;             // non-uniform sizing paths (local h, sparse
                                 // contact grid) — a file mesh is not a grid
    std::cout << "[FDEM3D] mesh = file: '" << path << "' — "
              << vpos.size() << " nodes, " << tets.size()
              << " tets, box " << W_ << " x " << D_ << " x " << H_ << " m\n";
    if (nGroups_ > 1) {
        std::cout << "[FDEM3D] physical groups: " << nGroups_ << " corps —";
        for (int g = 0; g < nGroups_; ++g)
            std::cout << " " << groupName_[g] << " (phase "
                      << phases_.name[grpPhase[g]] << ")";
        std::cout << "\n";
    }
    buildFromTets(vpos, tets, tetGrain, grpPhase);
    // groupe par element (ordre des elements = ordre des tets)
    elemGroup_.assign(el_.size(), 0);
    for (std::size_t e = 0; e < el_.size(); ++e)
        elemGroup_[e] = tetGroupTmp_[e];
    {   // WP2 : bilan des joints d interface entre corps lies
        long nb = 0;
        for (const auto& J : jt_)
            if (elemGroup_[J.eA] != elemGroup_[J.eB]) ++nb;
        if (nb > 0)
            std::cout << "[FDEM3D] groupBond: " << nb
                      << " joints d interface crees entre corps lies\n";
    }
    // vitesse initiale par groupe : groupVel.<nom> = "vx vy vz"
    for (int g = 0; g < nGroups_; ++g) {
        std::string vs = cfg_.gets("groupVel." + groupName_[g], "");
        if (vs.empty()) continue;
        std::istringstream iss(vs);
        double vx = 0, vy = 0, vz = 0;
        if (!(iss >> vx >> vy >> vz))
            throw std::runtime_error("groupVel." + groupName_[g]
                                     + " : attendu \"vx vy vz\"");
        long nset = 0;
        for (std::size_t e = 0; e < el_.size(); ++e)
            if (elemGroup_[e] == g)
                for (int a = 0; a < 4; ++a) {
                    v_[el_[e].n[a]] = {vx, vy, vz};
                    ++nset;
                }
        std::cout << "[FDEM3D] groupVel." << groupName_[g] << " = (" << vx
                  << ", " << vy << ", " << vz << ") m/s sur " << nset
                  << " noeuds\n";
    }
    // groupe suivi dans history.csv (colonnes grpZ, grpVz)
    {
        std::string tg = cfg_.gets("trackGroup", "");
        if (!tg.empty()) {
            for (int g = 0; g < nGroups_; ++g)
                if (groupName_[g] == tg) trackGroup_ = g;
            if (trackGroup_ < 0)
                throw std::runtime_error("trackGroup: groupe '" + tg
                                         + "' inconnu");
        }
    }
    // ---- WP3 : cinematique par corps + jauges en tranche ---------------
    // trackGroups = "nom1 nom2 ..." : z et vz massiques de CHAQUE corps
    // liste (le bit ET le piston de la spec 005). gauge.<nom> = "z0 z1" :
    // sigma_zz moyen en volume dans la tranche [z0, z1] du corps <nom> —
    // la jauge de leur fig. 8. La tranche est figee en configuration de
    // REFERENCE (les deplacements d un impact sont de l ordre du mm pour
    // des tranches de cm : biais negligeable, et la liste d elements se
    // fige a l init, cout nul en course).
    {
        std::string tl = cfg_.gets("trackGroups", "");
        std::istringstream iss(tl);
        std::string nm;
        while (iss >> nm) {
            int gg = -1;
            for (int g = 0; g < nGroups_; ++g)
                if (groupName_[g] == nm) gg = g;
            if (gg < 0)
                throw std::runtime_error("trackGroups: groupe '" + nm
                                         + "' inconnu");
            trkGrps_.push_back(gg);
        }
    }
    for (int g = 0; g < nGroups_; ++g) {
        std::string gs = cfg_.gets("gauge." + groupName_[g], "");
        if (gs.empty()) continue;
        std::istringstream iss(gs);
        double z0 = 0, z1 = 0;
        if (!(iss >> z0 >> z1) || !(z1 > z0))
            throw std::runtime_error("gauge." + groupName_[g]
                                     + " : attendu \"z0 z1\" avec z1 > z0");
        Gauge3 gg;
        gg.grp = g;
        gg.z0 = z0;
        gg.z1 = z1;
        for (std::size_t e = 0; e < el_.size(); ++e) {
            if (elemGroup_[e] != g) continue;
            Eigen::Vector3d c = Eigen::Vector3d::Zero();
            for (int a = 0; a < 4; ++a) c += X0_[el_[e].n[a]];
            c /= 4.0;
            if (c.z() >= z0 && c.z() <= z1) gg.elems.push_back((int)e);
        }
        if (gg.elems.empty())
            throw std::runtime_error("gauge." + groupName_[g]
                                     + " : aucune element dans la tranche");
        std::cout << "[FDEM3D] gauge " << groupName_[g] << " : "
                  << gg.elems.size() << " tets dans z = [" << z0 << ", "
                  << z1 << "] m\n";
        gauges_.push_back(std::move(gg));
    }
    hmin_ = 1e30;
    for (double h : hEl_) hmin_ = std::min(hmin_, h);
}

void Fdem3dSolver::buildMeshVoronoi() {
    double d = cfg_.reqd("grainSize");
    double jit = cfg_.getd("grainJitter", 0.5);
    int lloyd = cfg_.geti("lloydIters", 2);
    double mf = cfg_.getd("vertexMergeFrac", 0.12);
    int refine = cfg_.geti("refineLevels", 0);
    std::string seeding = cfg_.gets("grainSeeding", "hex");
    if (seeding != "hex" && seeding != "random")
        throw std::runtime_error("grainSeeding must be hex | random (got '"
                                 + seeding + "')");
    std::mt19937 rng(cfg_.geti("seed", 12345));

    Tessellation3 T = Tessellation3::build(W_, D_, H_, d, jit, lloyd, mf,
                                           refine, phases_.fraction, rng,
                                           seeding == "random");
    nGrains_ = T.nGrains;

    std::vector<std::array<int, 4>> tets;
    std::vector<int> tetGrain;
    tets.reserve(T.tet.size());
    tetGrain.reserve(T.tet.size());
    for (const auto& t : T.tet) {
        tets.push_back(t.v);
        tetGrain.push_back(t.grain);
    }
    buildFromTets(T.vtx, tets, tetGrain, T.phaseOfGrain);

    // the length scale of the contact machinery and the CFL is now set by
    // the smallest element the tessellation produced (inscribed size 6V/A)
    hmin_ = 1e30;
    for (double h : hEl_) hmin_ = std::min(hmin_, h);
}

// ---------------------------------------------------------------------------
// Shared topology builder: per-tet node quadruples, cohesive joints on the
// doubly-shared virtual faces (outward-oriented, node pairs matched by
// virtual id), exterior faces, masses, tension grips. Virtual vertex ids
// come from the caller: two tets are joined iff they share a vertex triple,
// inside grains and across them alike.
// ---------------------------------------------------------------------------
void Fdem3dSolver::buildFromTets(const std::vector<Eigen::Vector3d>& vpos,
                                 const std::vector<std::array<int, 4>>& tets,
                                 const std::vector<int>& tetGrain,
                                 const std::vector<int>& grainPhase) {
    // face registry: sorted virtual triple -> owners (elem + 3 (virt, real))
    struct Owner { int elem; std::array<int, 3> vv, nn; };
    std::map<std::array<int, 3>, std::vector<Owner>> faces;
    el_.reserve(tets.size());
    hEl_.reserve(tets.size());

    for (std::size_t tId = 0; tId < tets.size(); ++tId) {
        std::array<int, 4> vv = tets[tId];
        int base = (int)X0_.size();
        Elem e;
        for (int q = 0; q < 4; ++q) {
            X0_.push_back(vpos[vv[q]]);
            elemOf_.push_back((int)el_.size());
            vOf_.push_back(vv[q]);
            e.n[q] = base + q;
        }
        // ensure positive volume: swap two vertices if needed
        Eigen::Matrix3d J;
        J.col(0) = X0_[e.n[1]] - X0_[e.n[0]];
        J.col(1) = X0_[e.n[2]] - X0_[e.n[0]];
        J.col(2) = X0_[e.n[3]] - X0_[e.n[0]];
        double det = J.determinant();
        if (det < 0) {
            std::swap(e.n[2], e.n[3]);
            std::swap(vv[2], vv[3]);
            J.col(1) = X0_[e.n[2]] - X0_[e.n[0]];
            J.col(2) = X0_[e.n[3]] - X0_[e.n[0]];
            det = -det;
        }
        e.V0 = det / 6.0;
        if (e.V0 <= 0) throw std::runtime_error("degenerate tet");
        Eigen::Matrix3d Jinv = J.inverse();
        // grad N_a: N0 = 1 - xi - eta - zeta, N1 = xi, ...
        e.dN.col(1) = Jinv.row(0);
        e.dN.col(2) = Jinv.row(1);
        e.dN.col(3) = Jinv.row(2);
        e.dN.col(0) = -(e.dN.col(1) + e.dN.col(2) + e.dN.col(3));
        e.grain = tetGrain.empty() ? 0 : tetGrain[tId];
        e.phase = grainPhase.empty() ? 0 : grainPhase[e.grain];
        // inscribed-sphere diameter 6V/A — the local length scale of the
        // joint penalty and the contact caps on the voronoi mesh
        double Atot = 0.0;
        const int fi4[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
        for (const auto& fi : fi4) {
            Eigen::Vector3d A = X0_[e.n[fi[0]]], B = X0_[e.n[fi[1]]],
                            C = X0_[e.n[fi[2]]];
            Atot += 0.5 * (B - A).cross(C - A).norm();
        }
        hEl_.push_back(6.0 * e.V0 / Atot);
        int id = (int)el_.size();
        el_.push_back(e);

        // four faces, ordered OUTWARD (right-hand rule away from the
        // opposite vertex)
        const int faceIdx[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
        for (const auto& fi : faceIdx) {
            std::array<int, 3> fvv = {vv[fi[0]], vv[fi[1]], vv[fi[2]]};
            std::array<int, 3> fnn = {e.n[fi[0]], e.n[fi[1]], e.n[fi[2]]};
            // verify outwardness against the tet centroid
            Eigen::Vector3d A = X0_[fnn[0]], B = X0_[fnn[1]], C = X0_[fnn[2]];
            Eigen::Vector3d cen = 0.25 * (X0_[e.n[0]] + X0_[e.n[1]]
                                          + X0_[e.n[2]] + X0_[e.n[3]]);
            Eigen::Vector3d nrm = (B - A).cross(C - A);
            if (nrm.dot((A + B + C) / 3.0 - cen) < 0) {
                std::swap(fvv[1], fvv[2]);
                std::swap(fnn[1], fnn[2]);
            }
            std::array<int, 3> key = fvv;
            std::sort(key.begin(), key.end());
            faces[key].push_back({id, fvv, fnn});
        }
    }

    for (auto& [key, lst] : faces) {
        if (lst.size() == 2) {
            // V1 — physical groups : AUCUN joint entre deux groupes ; les
            // deux faces deviennent exterieures et l'interface est portee
            // par le contact general (penalite ou potentiel).
            // WP2 : SAUF si la paire est declaree liee (groupBond.<A>.<B> =
            // joints) — l interface conforme recoit alors des joints
            // cohesifs ordinaires, et rebindVertex liera ses noeuds comme
            // partout ailleurs (insertion adaptative comprise).
            if (!tetGroupTmp_.empty()
                && tetGroupTmp_[lst[0].elem] != tetGroupTmp_[lst[1].elem]
                && !(!gbond_.empty()
                     && gbond_[(std::size_t)tetGroupTmp_[lst[0].elem]
                               * nGroups_ + tetGroupTmp_[lst[1].elem]])) {
                exterior_.push_back({lst[0].elem, lst[0].nn});
                exterior_.push_back({lst[1].elem, lst[1].nn});
                continue;
            }
            Joint J;
            J.eA = lst[0].elem;
            J.eB = lst[1].elem;
            J.a = lst[0].nn;                 // outward-ordered in A
            for (int q = 0; q < 3; ++q)      // pair B nodes by virtual id
                for (int r = 0; r < 3; ++r)
                    if (lst[1].vv[r] == lst[0].vv[q]) J.b[q] = lst[1].nn[r];
            Eigen::Vector3d A = X0_[J.a[0]], B = X0_[J.a[1]], C = X0_[J.a[2]];
            J.A0 = 0.5 * (B - A).cross(C - A).norm();
            for (auto& sl : J.slip) sl.setZero();
            jt_.push_back(J);
        } else if (lst.size() == 1) {
            exterior_.push_back({lst[0].elem, lst[0].nn});
        } else {
            throw std::runtime_error("mesh topology error: a face is shared "
                                     "by more than two tets (vertex weld "
                                     "produced a non-manifold mesh)");
        }
    }

    u_.assign(X0_.size(), Eigen::Vector3d::Zero());
    v_.assign(X0_.size(), Eigen::Vector3d::Zero());
    f_.assign(X0_.size(), Eigen::Vector3d::Zero());
    m_.assign(X0_.size(), 0.0);
    flag_.assign(X0_.size(), FREE);
    for (const auto& e : el_)
        for (int a = 0; a < 4; ++a)
            m_[e.n[a]] += phases_.mat[e.phase].rho * e.V0 / 4.0;

    if (scen_ == Scenario::TENSION) {
        for (int i = 0; i < (int)X0_.size(); ++i) {
            if (X0_[i].z() < 1e-9)           flag_[i] = FIXED;
            else if (X0_[i].z() > H_ - 1e-9) flag_[i] = PRESCRIBED;
        }
    }
}

// ---------------------------------------------------------------------------
// Per-joint cohesive properties — the 2D GBM rule verbatim. Intra-grain
// joints (and every joint of the grid mesh) carry the bulk material of
// their phase. Grain-boundary joints take the MEAN of the two neighbouring
// phases times the alpha attenuation factors; heterophase boundaries get
// the extra heteroFactor on the strength-like properties. The penalty uses
// the local element size for the voronoi mesh (elements are not uniform)
// and the global hmin for the grid mesh (bit-compatible with the pre-GBM
// behaviour).
// ---------------------------------------------------------------------------
void Fdem3dSolver::assignJointProps() {
    // Adaptive insertion: the penalty only serves ACTIVATED joints as
    // unloading/contact stiffness (bonded faces are kinematic), so it can be
    // much softer than the intrinsic factor — that is where the dt gain of
    // Yan et al. comes from. Same convention as the 2D solver.
    double pf = adaptive_ ? cfg_.getd("insertionPenaltyFactor", 4.0)
                          : cfg_.getd("jointPenaltyFactor", 20.0);
    for (auto& J : jt_) {
        const Elem& A = el_[J.eA];
        const Elem& B = el_[J.eB];
        const Material& mA = phases_.mat[A.phase];
        const Material& mB = phases_.mat[B.phase];
        double E, ft, coh, Gf, GfII, phiDeg;
        if (!voronoi_ || A.grain == B.grain) {
            J.type = 0;                                // intra-grain: bulk
            E = mA.E; ft = mA.ft; coh = mA.cohesion;
            Gf = mA.Gf; GfII = mA.gfShearFactor * mA.Gf;
            phiDeg = mA.phiDeg;
        } else {
            bool hetero = A.phase != B.phase;
            J.type = hetero ? 2 : 1;
            double s = hetero ? phases_.heteroFactor : 1.0;
            E    = phases_.aE   * 0.5 * (mA.E + mB.E);
            ft   = s * phases_.aTen * 0.5 * (mA.ft + mB.ft);
            coh  = s * phases_.aCoh * 0.5 * (mA.cohesion + mB.cohesion);
            Gf   = s * phases_.aGf  * 0.5 * (mA.Gf + mB.Gf);
            GfII = s * phases_.aGf  * 0.5 * (mA.gfShearFactor * mA.Gf
                                             + mB.gfShearFactor * mB.Gf);
            phiDeg = phases_.aFric * 0.5 * (mA.phiDeg + mB.phiDeg);
        }
        // defense in depth: the attenuated MEANS must stay physical (ft = 0
        // would make dnF infinite, the envelope NaN, and the joint
        // unbreakable; tan(phi >= 90 deg) would flip the friction sign)
        if (!(ft > 0.0) || !(coh > 0.0) || !(E > 0.0)
            || !(phiDeg >= 0.0 && phiDeg < 89.0))
            throw std::runtime_error("assignJointProps: attenuated joint "
                "properties left the physical range (ft/coh/E must be > 0, "
                "friction angle in [0, 89) deg) — check gb* factors");
        double h = voronoi_ ? 0.5 * (hEl_[J.eA] + hEl_[J.eB]) : hmin_;
        J.pj = pf * E / h;
        J.ft = ft;
        J.coh = coh;
        J.Gf = Gf;
        J.GfII = GfII;
        J.dnE = ft / J.pj;
        // Critical opening / slip: the softening branch must enclose exactly
        // GfI (resp. GfII) — linear branch has area ft w / 2, the f(D) branch
        // has area ft w I with I = int f(D) dD (eq. 13/15), as in 2D.
        double kI = yanSoft_ ? 1.0 / yanI_ : 2.0;
        J.dnF = J.dnE + kI * Gf / ft;                  // mode-I critical opening
        J.slipF = kI * GfII / coh;                     // mode-II softening slip
        J.tanPhi = std::tan(phiDeg * M_PI / 180.0);
    }
    applyJointStatistics();
}

// ---------------------------------------------------------------------------
// Statistical joint strengths — the 2D implicit-DFH bridge in 3D. With
// jointWeibullM = m > 1, every joint's ft and cohesion are multiplied by a
// Weibull(m) factor of MEAN 1; fracture energies are NOT scaled and the
// openings dnF/slipF are recomputed. strengthCorrLength = 0 draws
// independently per joint; > 0 samples ONE Gaussian random field
// (RandomField3) through the Gaussian copula at the joint face centroid —
// the field lives in SPACE, independent of the mesh (fieldSeed controls it
// separately from the mesh seed), so two meshes see the same weak zones.
// Anisotropy: strengthCorrLengthB across the texture plane and
// strengthCorrAngleDeg tilting that plane about the y-axis (foliation).
// ---------------------------------------------------------------------------
void Fdem3dSolver::applyJointStatistics() {
    const bool wGf =
        cfg_.gets("weibullScope", "strength") == "strengthGf";
    double m = cfg_.getd("jointWeibullM", 0.0);
    // Effet d'echelle statistique de Weibull, MEME CONVENTION que les VUMAT
    // DP-DFH d'Abaqus (VUMATS/dfh/vumat_kstdfh_psivar.f) :
    //     sig_k = sigw * (Zeff/V_el)^(1/m),   V_el = charLength^3
    // Ici le "point materiel" d'un joint est la face entre deux tetras : son
    // volume represente est la moyenne des deux volumes adjacents. Le rapport
    // DP-DFH (eq. 42, §13.1) impose ce recalage pour toute etude d'objectivite
    // STRUCTURALE : la variante a V_el fige (vumat_psivar_rc99_veff1.f, qui ne
    // differe que par V_el = 1 mm3) est reservee aux controles au point
    // materiel et desactive l'effet d'echelle.
    // OPT-IN : sans jointSizeEffect, ce bloc est saute et les resultats sont
    // inchanges au bit pres.
    bool szOn = cfg_.getb("jointSizeEffect", false);
    if (m <= 0.0 && !szOn) return;                     // deterministic joints
    if (m > 0.0 && m <= 1.0)
        throw std::runtime_error("jointWeibullM must be > 1 (typical rock "
                                 "values 5-30)");
    if (m <= 0.0) for (auto& J : jt_) J.stat = 1.0;    // taille seule
    double ell = cfg_.getd("strengthCorrLength", 0.0);
    unsigned fseed = (unsigned)cfg_.geti("fieldSeed",
                                         cfg_.geti("seed", 12345) + 777);
    double gam = std::tgamma(1.0 + 1.0 / m);
    auto weib = [&](double u) {
        u = std::clamp(u, 1e-12, 1.0 - 1e-12);
        return std::pow(-std::log(1.0 - u), 1.0 / m) / gam;
    };

    double xmin = 1e300, xmax = 0.0, xsum = 0.0;
    if (m > 0.0) {
    if (ell > 0.0) {
        double ellB = cfg_.getd("strengthCorrLengthB", ell);
        double ang = cfg_.getd("strengthCorrAngleDeg", 0.0);
        RandomField3 F(W_, D_, H_, ell, ellB, ang, fseed);
        for (auto& J : jt_) {
            Eigen::Vector3d mid = (X0_[J.a[0]] + X0_[J.a[1]] + X0_[J.a[2]])
                                  / 3.0;
            double g = F(mid);
            J.stat = weib(0.5 * std::erfc(-g / std::sqrt(2.0)));   // Phi(g)
        }
    } else {
        std::mt19937 rng(fseed);
        std::uniform_real_distribution<double> U(0.0, 1.0);
        for (auto& J : jt_) J.stat = weib(U(rng));
    }
    }
    if (szOn) applyJointSizeEffect(m);
    for (auto& J : jt_) {
        J.ft *= J.stat;
        J.coh *= J.stat;
        // weibullScope — miroir exact du 2D, voir MatLaw.hpp
        if (wGf) { J.Gf *= J.stat; J.GfII *= J.stat; }
        J.dnE = J.ft / J.pj;
        double kI = yanSoft_ ? 1.0 / yanI_ : 2.0;
        J.dnF = J.dnE + kI * J.Gf / J.ft;
        J.slipF = kI * J.GfII / J.coh;
        xmin = std::min(xmin, J.stat);
        xmax = std::max(xmax, J.stat);
        xsum += J.stat;
    }
    if (m > 0.0)
        std::cout << "[FDEM3D] joint strength statistics: Weibull m = " << m
                  << (ell > 0.0 ? " correlated, ell = " + std::to_string(ell)
                                : std::string(" independent per joint"))
                  << ", factor mean/min/max = " << xsum / jt_.size() << "/"
                  << xmin << "/" << xmax << "\n";
    if (szOn)
        std::cout << "[FDEM3D] TOTAL ft factor (Weibull x taille) "
                     "mean/min/max = " << xsum / jt_.size() << "/"
                  << xmin << "/" << xmax << "\n";
}

// ---------------------------------------------------------------------------
// Effet d'echelle statistique : ft <- ft * (Zeff/V_J)^(1/m), la formule des
// VUMAT DP-DFH (sig_k = sigw*(Zeff/V_el)^(1/m)) transposee au joint cohesif.
// V_J = moyenne des volumes des deux tetras adjacents : le volume de matiere
// que le lien echantillonne. Gf n'est PAS recale (energie de fissuration =
// propriete du materiau, pas une resistance au pic gouvernee par le maillon
// faible) ; dnF/slipF sont donc recalcules par l'appelant apres coup, ce qui
// raccourcit la branche adoucissante quand ft monte — consequence assumee.
// Le facteur est replie dans J.stat, donc le champ ftScale des VTU montre le
// facteur TOTAL sans changer l'ecriture.
// ---------------------------------------------------------------------------
void Fdem3dSolver::applyJointSizeEffect(double mWeib) {
    double mS = cfg_.getd("jointSizeEffectM", mWeib);
    if (!(mS > 1.0))
        throw std::runtime_error("jointSizeEffect: exposant invalide — poser "
            "jointSizeEffectM > 1 (ou jointWeibullM). C'est le meme m que la "
            "dispersion : eq. 42 du rapport DP-DFH");
    // Zeff : volume de REFERENCE auquel ft/cohesion de la config sont
    // declares. Doit etre une constante PHYSIQUE, jamais deduite du maillage —
    // un Zeff qui suivrait la moyenne du maillage rendrait le facteur moyen
    // egal a 1 et desactiverait l'effet d'echelle (le piege du §13.1).
    // Defaut 1e-9 m3 = 1 mm3, le Zeff par defaut des VUMAT (echelle de
    // l'indentation, cf. Table 3 du rapport).
    double Zeff = cfg_.getd("jointZeff", 1e-9);
    if (!(Zeff > 0.0))
        throw std::runtime_error("jointZeff doit etre > 0 [m^3]");
    double cap = cfg_.getd("jointSizeEffectClamp", 5.0);
    if (!(cap >= 1.0))
        throw std::runtime_error("jointSizeEffectClamp doit etre >= 1");
    double fmin = 1e300, fmax = 0.0, fsum = 0.0, vmin = 1e300, vmax = 0.0;
    std::size_t nclip = 0;
    for (auto& J : jt_) {
        double Vj = 0.5 * (el_[J.eA].V0 + el_[J.eB].V0);
        double f = 1.0;
        if (Vj > 0.0) f = std::pow(Zeff / Vj, 1.0 / mS);
        if (f > cap)            { f = cap;       ++nclip; }
        else if (f < 1.0 / cap) { f = 1.0 / cap; ++nclip; }
        J.stat *= f;
        fmin = std::min(fmin, f); fmax = std::max(fmax, f); fsum += f;
        vmin = std::min(vmin, Vj); vmax = std::max(vmax, Vj);
    }
    std::cout << "[FDEM3D] effet d'echelle (Zeff/V_J)^(1/m) : Zeff = " << Zeff
              << " m^3, m = " << mS << ", V_J min/max = " << vmin << "/"
              << vmax << " m^3, facteur mean/min/max = " << fsum / jt_.size()
              << "/" << fmin << "/" << fmax << "\n";
    if (nclip)
        std::cout << "[FDEM3D] WARNING: " << nclip << " joints bornes a "
                  << cap << "x (jointSizeEffectClamp) — maillage tres "
                     "heterogene ou Zeff mal choisi\n";
}

// ===========================================================================
// Adaptive insertion — 3D port of the FdemSolver machinery (Yan, Zheng &
// Wang, IJRMMS 169, 2023). Nodes are already duplicated per tet, so "shared"
// is enforced kinematically: the co-located copies of an original vertex are
// bound into groups that integrate as one node. Binding groups are the
// connected components of the TET FAN around the vertex, two tets being
// connected when the face between them is still bonded; activating a joint
// re-runs the union-find at the face's three vertices, which reproduces the
// progressive node splitting of the 2D fig. 7 in 3D (a crack-front vertex
// stays whole while the fan is still connected around it).
// ===========================================================================
void Fdem3dSolver::buildBindingTables() {
    nVert_ = 0;
    for (int v : vOf_) nVert_ = std::max(nVert_, v + 1);
    copiesOfVert_.assign(nVert_, {});
    for (int i = 0; i < (int)X0_.size(); ++i)
        copiesOfVert_[vOf_[i]].push_back(i);
    jointsOfVert_.assign(nVert_, {});
    for (int jI = 0; jI < (int)jt_.size(); ++jI)
        for (int k = 0; k < 3; ++k)
            jointsOfVert_[vOf_[jt_[jI].a[k]]].push_back(jI);
    grpsOfVert_.assign(nVert_, {});
    for (int v = 0; v < nVert_; ++v) rebindVertex(v);
}

void Fdem3dSolver::rebindVertex(int v) {
    const auto& copies = copiesOfVert_[v];
    auto& grps = grpsOfVert_[v];
    grps.clear();
    if (copies.empty()) return;
    std::vector<int> elems;
    elems.reserve(copies.size());
    for (int i : copies) elems.push_back(elemOf_[i]);
    auto local = [&](int e) {
        for (int k = 0; k < (int)elems.size(); ++k)
            if (elems[k] == e) return k;
        return -1;
    };
    std::vector<int> par(elems.size());
    for (int k = 0; k < (int)par.size(); ++k) par[k] = k;
    std::function<int(int)> find = [&](int x) {
        while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; }
        return x;
    };
    for (int jI : jointsOfVert_[v]) {
        const Joint& J = jt_[jI];
        if (!J.bonded) continue;
        int a = local(J.eA), b = local(J.eB);
        if (a < 0 || b < 0) continue;
        par[find(a)] = find(b);
    }
    std::vector<int> root(copies.size());
    for (int c = 0; c < (int)copies.size(); ++c) root[c] = find(c);
    std::vector<char> done(copies.size(), 0);
    for (int c = 0; c < (int)copies.size(); ++c) {
        if (done[c]) continue;
        std::vector<int> g;
        for (int d = c; d < (int)copies.size(); ++d)
            if (!done[d] && root[d] == root[c]) { done[d] = 1; g.push_back(copies[d]); }
        grps.push_back(std::move(g));
    }
}

// Insertion criterion on every bonded face: average of the two tets' GLOBAL
// stress tensors projected on the current face frame; activate when
// sigma_n >= ft or |tau| >= fs with fs = c - sigma_n tan(phi) in compression.
void Fdem3dSolver::insertionSweep() {
    struct Hit { int jI; double sig, fs; Eigen::Vector3d tauV; };
    std::vector<Hit> hits;
    auto testJoint = [&](int jI, std::vector<Hit>& out) {
        const Joint& J = jt_[jI];
        if (!J.bonded) return;
        Eigen::Vector3d P1 = 0.5 * ((X0_[J.a[0]] + u_[J.a[0]]) + (X0_[J.b[0]] + u_[J.b[0]]));
        Eigen::Vector3d P2 = 0.5 * ((X0_[J.a[1]] + u_[J.a[1]]) + (X0_[J.b[1]] + u_[J.b[1]]));
        Eigen::Vector3d P3 = 0.5 * ((X0_[J.a[2]] + u_[J.a[2]]) + (X0_[J.b[2]] + u_[J.b[2]]));
        Eigen::Vector3d nr = (P2 - P1).cross(P3 - P1);
        double nn = nr.norm();
        if (nn < 1e-18) return;
        Eigen::Vector3d n = nr / nn;
        Eigen::Matrix3d S = 0.5 * (el_[J.eA].sigG + el_[J.eB].sigG);
        Eigen::Vector3d tvec = S * n;
        double sig = n.dot(tvec);
        Eigen::Vector3d tauV = tvec - sig * n;
        // DIF de Yang : le critere doit voir la resistance DYNAMIQUE, sinon
        // le joint s insere au seuil statique et le facteur applique ensuite
        // serait sans effet sur l INSTANT d insertion. Le terme de frottement
        // -sig tan(phi) n est PAS amplifie (leur choix, source Zhao). Le taux
        // est la moyenne non ponderee des deux tetras, par coherence avec la
        // moyenne des contraintes 0,5 (sigG_A + sigG_B) deja en place.
        double dT = 1.0, dC = 1.0;
        if (difOn_) {
            double er = 0.5 * (el_[J.eA].edot + el_[J.eB].edot);
            dT = rockim::difTensionYang(er, difExpT_);
            dC = rockim::difCompressionYang(er);
        }
        double fs = dC * J.coh
                  + J.tanPhi * rockim::mcFrictionTerm(sig, J.ft, yangEnv_);
        if (fs < 0.0) fs = 0.0;
        if (sig >= dT * J.ft || tauV.norm() >= fs)
            out.push_back({jI, sig, fs, tauV});
    };
#ifdef _OPENMP
    #pragma omp parallel
    {
        std::vector<Hit> mine;
        #pragma omp for schedule(static) nowait
        for (int jI = 0; jI < (int)jt_.size(); ++jI) testJoint(jI, mine);
        #pragma omp critical
        hits.insert(hits.end(), mine.begin(), mine.end());
    }
#else
    for (int jI = 0; jI < (int)jt_.size(); ++jI) testJoint(jI, hits);
#endif
    if (hits.empty()) return;
    std::sort(hits.begin(), hits.end(),
              [](const Hit& x, const Hit& y) { return x.jI < y.jI; });
    for (const Hit& h : hits) activateJoint(h.jI, h.sig, h.tauV, h.fs);
}

// Stress continuity at insertion, as in 2D: opening offset so the elastic
// branch reproduces min(sigma, ft) at zero geometric opening, and a vector
// slip offset so the trial shear pj*(dt3 - slip) equals the transmitted
// tangential traction (clamped to the current Coulomb cap) at dt3 = 0.
void Fdem3dSolver::activateJoint(int jI, double sig,
                                 const Eigen::Vector3d& tauV, double fsNow) {
    Joint& J = jt_[jI];
    if (!J.bonded) return;
    J.bonded = false;
    // ---- DIF de Yang et al. 2025, FIGE ICI ------------------------------
    // Applique AVANT les decalages de continuite de contrainte ci-dessous,
    // qui lisent J.ft, J.coh et J.pj. Se COMPOSE avec le facteur statistique
    // deja porte par ft et coh (multiplicatif). fsNow arrive deja calcule
    // avec dC par insertionSweep, donc il reste coherent avec le J.coh neuf.
    // Les ouvertures critiques sont recalculees : dnE = ft/pj suit le DIF,
    // tandis que kI Gf/ft ne bouge pas puisque ft et Gf recoivent le MEME
    // facteur — dnF est donc quasi invariante et le compteur
    // d endommagement reste coherent.
    if (difOn_) {
        double er = 0.5 * (el_[J.eA].edot + el_[J.eB].edot);
        double dT = rockim::difTensionYang(er, difExpT_);
        double dC = rockim::difCompressionYang(er);
        J.ft   *= dT;
        J.Gf   *= dT;
        J.coh  *= dC;
        J.GfII *= dC;
        J.difT = dT; J.difC = dC; J.edotIns = er;
        double kI = yanSoft_ ? 1.0 / yanI_ : 2.0;
        J.dnE   = J.ft / J.pj;
        J.dnF   = J.dnE + kI * J.Gf / J.ft;
        J.slipF = kI * J.GfII / J.coh;
    }
    J.dn0 = std::min(sig, J.ft) / J.pj;
    Eigen::Vector3d tau0 = tauV;
    double tn = tau0.norm();
    if (tn > fsNow && tn > 0.0) tau0 *= fsNow / tn;
    for (int k = 0; k < 3; ++k) J.slip[k] = -tau0 / J.pj;
    ++nInserted_;
    rebindVertex(vOf_[J.a[0]]);
    rebindVertex(vOf_[J.a[1]]);
    rebindVertex(vOf_[J.a[2]]);
}

void Fdem3dSolver::placeTool() {
    if (scen_ == Scenario::TENSION) return;
    tool_.mass   = cfg_.getd("toolMass", 0.5);
    tool_.radius = cfg_.getd("toolRadius", 0.015);
    // ---- force volumique et tri des fragments : voir Fdem3dSolver.hpp -----
    gravity_ = cfg_.getd("gravity", 0.0);
    toolStop_ = cfg_.getd("toolStop", 0.0);
    brushStart_ = cfg_.getd("fragBrushStart", 0.0);
    if (toolStop_ > 0.0 && brushStart_ > 0.0 && brushStart_ <= toolStop_)
        throw std::runtime_error("fragBrushStart doit etre APRES toolStop : "
                                 "c'est l'intervalle de repos qui rend le tri "
                                 "interpretable (Yang et al. 2025, etape 1)");
    if (brushStart_ > 0.0) {
        brushV0_   = cfg_.getd("fragBrushV0", 2.5e-3);
        brushA_    = cfg_.getd("fragBrushAccel", 98.1);
        brushBeta_ = cfg_.getd("fragBrushBeta", 0.8);
        brushZeroV_ = cfg_.gets("fragBrushZeroV", "false") == "true";
        brushDir_  = Eigen::Vector3d(cfg_.getd("fragBrushDirX", 0.0),
                                     cfg_.getd("fragBrushDirY", 0.0),
                                     cfg_.getd("fragBrushDirZ", 1.0));
        if (brushDir_.norm() < 1e-12)
            throw std::runtime_error("fragBrushDir est nul : donner une "
                                     "direction (opposee a l'impact)");
        brushDir_.normalize();
        if (brushBeta_ <= 0.0)
            throw std::runtime_error("fragBrushBeta doit etre > 0");
    }
    toolVCap_ = cfg_.getd("toolImpulseCap", 0.0);
    if (toolVCap_ > 0.0)
        std::cout << "[FDEM3D] toolImpulseCap = " << toolVCap_
                  << " : |Fc| <= kappa * 2 v_outil * m / dt (borne du choc "
                     "elastique contre une masse infinie)\n";
    double gap = cfg_.getd("toolGap", 1e-4);
    // toolShape = sphere | flat ('disc' accepted as the 2D synonym) | none.
    // none (V1) : PAS d'outil analytique — l'outil est un corps MAILLE
    // (physical group + groupVel), et toolContact / tool_.integrate sont
    // court-circuites DUR (lecon du disque fantome brezilien 2D).
    std::string sh = cfg_.gets("toolShape", "sphere");
    if (sh == "none") {
        toolNone_ = true;
        tool_.free = false;
        tool_.x = {1e9, 1e9, 1e9};         // hors de portee de tout
        tool_.v.setZero();
        return;
    }
    if (sh != "sphere" && sh != "disc" && sh != "flat")
        throw std::runtime_error("toolShape must be sphere | flat | none (3D)");
    tool_.flat = sh == "flat" && scen_ == Scenario::PERCUSSION;
    if (scen_ == Scenario::PERCUSSION) {
        tool_.free = true;
        double zTip = tool_.flat ? H_ + gap : H_ + tool_.radius + gap;
        tool_.x = {cfg_.getd("toolX", 0.5 * W_), cfg_.getd("toolY", 0.5 * D_),
                   zTip};
        tool_.v = {0.0, 0.0, -cfg_.getd("impactSpeed", 8.0)};
    } else {
        tool_.free = false;
        double depth = cfg_.getd("cutDepth", 0.004);
        tool_.x = {cfg_.getd("toolX", -tool_.radius - gap), 0.5 * D_,
                   H_ - depth + tool_.radius};
        tool_.v = {cfg_.getd("cutSpeed", 10.0), 0.0, 0.0};
    }
}

void Fdem3dSolver::setupBoundaries() {
    cAbs_.assign(X0_.size(), Eigen::Vector3d::Zero());
    kAbs_.assign(X0_.size(), Eigen::Vector3d::Zero());
    if (scen_ == Scenario::TENSION) return;

    std::string ab = cfg_.gets("absorbing", "none");
    if (ab != "none" && ab != "sides" && ab != "all")
        throw std::runtime_error("absorbing must be none | sides | all");

    double sF = cfg_.getd("absorbSpringFactor", 1.0);
    double Rx = cfg_.getd("absorbSpringR", 0.5 * W_);
    double Ry = cfg_.getd("absorbSpringR", 0.5 * D_);
    double Rz = cfg_.getd("absorbSpringR", H_);
    double tol = 1e-9;

    for (const auto& bf : exterior_) {
        // impedances of the LOCAL phase: with mineral phases the truncated
        // continuum behind each boundary face is the one of the face's grain
        const Material& mp = phases_.mat[el_[bf.elem].phase];
        double G = mp.G();
        Eigen::Vector3d A = X0_[bf.n[0]], B = X0_[bf.n[1]], C = X0_[bf.n[2]];
        double At3 = 0.5 * (B - A).cross(C - A).norm() / 3.0;  // per node
        bool xlo = A.x() < tol && B.x() < tol && C.x() < tol;
        bool xhi = A.x() > W_ - tol && B.x() > W_ - tol && C.x() > W_ - tol;
        bool ylo = A.y() < tol && B.y() < tol && C.y() < tol;
        bool yhi = A.y() > D_ - tol && B.y() > D_ - tol && C.y() > D_ - tol;
        bool zlo = A.z() < tol && B.z() < tol && C.z() < tol;
        int nAxis = xlo || xhi ? 0 : (ylo || yhi ? 1 : (zlo ? 2 : -1));
        if (nAxis < 0) continue;
        double R = nAxis == 0 ? Rx : (nAxis == 1 ? Ry : Rz);
        bool lateral = nAxis != 2;
        for (int nid : bf.n) {
            if (lateral && ab != "none") {
                for (int a = 0; a < 3; ++a) {
                    double c = (a == nAxis ? mp.cP() : mp.cS());
                    cAbs_[nid](a) += mp.rho * c * At3;
                    kAbs_[nid](a) += sF * G / (a == nAxis ? R : 2.0 * R) * At3;
                }
            }
            if (nAxis == 2) {
                if (ab == "all") {
                    for (int a = 0; a < 3; ++a) {
                        double c = (a == 2 ? mp.cP() : mp.cS());
                        cAbs_[nid](a) += mp.rho * c * At3;
                        kAbs_[nid](a) += sF * G / (a == 2 ? R : 2.0 * R) * At3;
                    }
                } else {                       // percussion AND shear: the
                    flag_[nid] = FIXED;        // block needs its support
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Confinement triaxial 3D — portage direct du 2D : la pression est une charge
// SUIVEUSE sur les faces exterieures LATERALES d'origine (une cellule
// triaxiale a une membrane sur le fut, pas sur les plateaux), recalculee sur
// la geometrie courante, lumped A/3 par noeud, montee en rampe cosinus.
// Comme en 2D, les faces creees par la fissuration ne recoivent RIEN (la
// membrane presse la surface de l'eprouvette, pas l'interieur des fissures).
// ---------------------------------------------------------------------------
void Fdem3dSolver::setupConfinement() {
    confP_ = cfg_.getd("confiningPressure", 0.0);
    confRamp_ = cfg_.getd("confiningRamp", 2e-4);
    if (confP_ == 0.0) return;
    std::string cf = cfg_.gets("confineFaces", "lateral");
    if (cf != "lateral" && cf != "all")
        throw std::runtime_error("confineFaces must be lateral | all");
    for (const auto& bf : exterior_) {
        Eigen::Vector3d A = X0_[bf.n[0]], B = X0_[bf.n[1]], C = X0_[bf.n[2]];
        Eigen::Vector3d n = (B - A).cross(C - A);
        double nn = n.norm();
        if (nn < 1e-18) continue;
        if (cf == "lateral" && std::abs(n.z() / nn) > 0.5) continue;
        confFaces_.push_back(bf);
    }
    std::cout << "[FDEM3D] confinement: " << confP_ / 1e6 << " MPa sur "
              << confFaces_.size() << " faces laterales d'origine, rampe "
              << confRamp_ << " s\n";
}

void Fdem3dSolver::confiningForces() {
    if (confP_ == 0.0) return;
    double p = confP_;
    if (confRamp_ > 0.0 && t_ < confRamp_)
        p *= 0.5 * (1.0 - std::cos(M_PI * t_ / confRamp_));
    for (const auto& bf : confFaces_) {
        Eigen::Vector3d A = X0_[bf.n[0]] + u_[bf.n[0]];
        Eigen::Vector3d B = X0_[bf.n[1]] + u_[bf.n[1]];
        Eigen::Vector3d C = X0_[bf.n[2]] + u_[bf.n[2]];
        // aire orientee courante : les faces exterieures sont construites
        // sortantes, la traction -p n s'applique vers l'interieur
        Eigen::Vector3d S = 0.5 * (B - A).cross(C - A);
        Eigen::Vector3d F = -p * S / 3.0;              // par noeud
        f_[bf.n[0]] += F;
        f_[bf.n[1]] += F;
        f_[bf.n[2]] += F;
        // V2/B4 : travail de la pression sur le solide (compteur en v-,
        // la correction leapfrog globale couvre le demi-pas comme ailleurs)
        confWork_ += dt_ * (F.dot(v_[bf.n[0]]) + F.dot(v_[bf.n[1]])
                            + F.dot(v_[bf.n[2]]));
    }
}

// contrainte laterale moyenne (sxx+syy)/2 dans le coeur (tiers central) —
// la jauge falsifiable du confinement, lue apres l'equilibrage
double Fdem3dSolver::achievedConfinement() const {
    double s = 0.0, V = 0.0;
    for (const auto& e : el_) {
        Eigen::Vector3d c = 0.25 * (X0_[e.n[0]] + X0_[e.n[1]]
                                    + X0_[e.n[2]] + X0_[e.n[3]]);
        if (std::abs(c.x() - 0.5 * W_) > 0.25 * W_
            || std::abs(c.y() - 0.5 * D_) > 0.25 * D_
            || std::abs(c.z() - 0.5 * H_) > 0.25 * H_) continue;
        s += 0.5 * (e.sigG(0, 0) + e.sigG(1, 1)) * e.V0;
        V += e.V0;
    }
    return V > 0.0 ? s / V : 0.0;
}

void Fdem3dSolver::computeStableDt() {
    std::vector<double> K(X0_.size());
    for (std::size_t i = 0; i < X0_.size(); ++i) {
        const Elem& e = el_[elemOf_[i]];
        double h = voronoi_ ? hEl_[elemOf_[i]] : hmin_;
        K[i] = 4.0 * phases_.mat[e.phase].E * h;
    }
    for (const auto& J : jt_) {
        double k = J.pj * J.A0 / 3.0;
        for (int q = 0; q < 3; ++q) { K[J.a[q]] += k; K[J.b[q]] += k; }
    }
    double nExtra = cfg_.getd("extraContacts", 2.0);
    double dtMin = 1e30;
    for (std::size_t i = 0; i < X0_.size(); ++i)
        // A1 : en Signorini l'outil n'a plus de raideur (condition sur la
        // VITESSE) — kp_ sort du budget. Gain structurel : le pas cesse de
        // dependre d'une penalite arbitraire.
        dtMin = std::min(dtMin, 2.0 * std::sqrt(m_[i] / (K[i]
                    + (toolSig_ ? 0.0 : nExtra * kp_))));
    // CFL from the TRUE minimum inscribed diameter, not the nominal grid
    // pitch. The Kuhn split makes tets thin BY CONSTRUCTION (median 6V/A is
    // ~0.4x the cell edge, jitter takes the worst to ~0.2x), so the nominal
    // hmin overestimates the stable dt by up to 5x. The intrinsic runs never
    // saw it: their pf = 20 joint springs dominated dtMin and provided the
    // safety margin by accident. insertion = adaptive removed those springs
    // from the bonded phase and the first homogeneous grid impact EXPLODED
    // (2.4 MJ of block KE from a 16 J impact) — measured 2026-08-07.
    double hCfl = hmin_;
    for (double h : hEl_) hCfl = std::min(hCfl, h);
    double cfl = hCfl / phases_.maxCp();
    // ---- borne DIFFUSIVE du terme visqueux 2 mu D (eq. 6 de Yan) --------
    // La viscosite longitudinale effective vaut 2 mu, donc la diffusivite de
    // quantite de mouvement est nu = 2 mu / rho et le critere explicite
    // dt <= h^2/(2 nu) donne dt <= rho h^2 / (4 mu). ELEMENT PAR ELEMENT :
    // ecrire cette borne avec hmin_ serait faux — en mesh = grid, hmin_ reste
    // le pas NOMINAL et Kuhn donne hEl = 0,414 a, donc la borne sortirait
    // 5,8 fois trop permissive, instable en silence. C est exactement le bug
    // CFL du 2026-08-07 (2,4 MJ d energie de bloc pour 16 J incidents).
    double dtVis = 1e30;
    if (viscOn_)
        for (std::size_t eI = 0; eI < el_.size(); ++eI)
            if (muEl_[eI] > 0.0)
                dtVis = std::min(dtVis, rhoP_[el_[eI].phase] * hEl_[eI]
                                        * hEl_[eI] / (4.0 * muEl_[eI]));
    double dtSpr = std::min(dtMin, cfl);
    dt_ = cfg_.getd("dtFactor", 0.15) * std::min(dtSpr, dtVis);
    if (viscOn_)
        std::cout << "[FDEM3D] borne diffusive rho h^2/(4 mu) = " << dtVis
                  << " s contre " << dtSpr << " s pour les ressorts"
                  << (dtVis < dtSpr
                      ? "  <-- C EST ELLE QUI COMMANDE LE PAS" : "")
                  << " (facteur " << (dtSpr / std::min(dtSpr, dtVis))
                  << " sur le cout du run)\n";
}

// ===========================================================================

void Fdem3dSolver::step() {
    for (auto& fi : f_) fi.setZero();
    tool_.F.setZero();
    // tri des fragments : armement a l'instant demande (voir Fdem3dSolver.hpp)
    if (brushStart_ > 0.0 && !brushArmed_ && t_ >= brushStart_) armBrush();

    if (f3Prof.on) {
        double t0 = f3now(); elementForces(); bodyForces();
        double t1 = f3now(); if (adaptive_) insertionSweep();
        double t2 = f3now(); jointForces();
        double t3 = f3now(); generalContact();
        double t4 = f3now(); toolContact();
        double t5 = f3now();
        f3Prof.tEl += t1 - t0; f3Prof.tIn += t2 - t1; f3Prof.tJt += t3 - t2;
        f3Prof.tGc += t4 - t3; f3Prof.tTc += t5 - t4;
        ++f3Prof.n;
    } else {
        elementForces();
        bodyForces();                      // gravite + anti-gravite du tri
        if (adaptive_) insertionSweep();   // before jointForces: a joint born
                                           // this step carries traction now
        jointForces();
        generalContact();
        toolContact();
    }
    confiningForces();                     // no-op si confiningPressure = 0
    if (confP_ > 0.0 && !confLatched_
        && t_ >= std::max(cfg_.getd("confineGaugeTime", 3.0 * confRamp_),
                          20.0 * dt_)) {
        confAchieved_ = achievedConfinement();
        confLatched_ = true;
    }

    if (scen_ == Scenario::TENSION) {
        gripF_.setZero();
        for (int i = 0; i < (int)X0_.size(); ++i)
            if (flag_[i] == PRESCRIBED) gripF_ += f_[i];
        double sigNow = std::abs(gripF_.z()) / (W_ * D_);
        sigmaPeak_ = std::max(sigmaPeak_, sigNow);
        // Arret post-rupture (opt-in stopPeakDrop, 0 = off) : quand la
        // contrainte est retombee sous (1 - drop) x pic, l'essai est fini —
        // le post-pic profond (bande entiere en contact de levres) coute le
        // prix D0 sans rien apporter a la mesure. Garde : pic significatif
        // (> 1 MPa) pour ne pas declencher sur le bruit d'avant-charge.
        if (stopDrop_ < 0.0) stopDrop_ = cfg_.getd("stopPeakDrop", 0.0);
        if (stopDrop_ > 0.0 && sigmaPeak_ > 1e6
            && sigNow < (1.0 - stopDrop_) * sigmaPeak_)
            peakStop_ = true;
    } else {
        peakF_ = std::max(peakF_, tool_.F.norm());
        work_ += -tool_.F.dot(tool_.v) * dt_;
    }

    integrate();
    t_ += dt_;
    if ((++stepCount_ & 1023) == 0) {
        // ---- E5 (2026-08-19), miroir du 2D : u_[0] peut etre un noeud FIXED,
        // donc toujours fini — le detecteur etait AVEUGLE. Echantillonnage de
        // tout le maillage a pas constant (~256 noeuds), decalage tournant.
        bool bad = !std::isfinite(work_);
        const std::size_t nN = X0_.size();
        const std::size_t stride = (nN > 256) ? nN / 256 : 1;
        const std::size_t off = (std::size_t)((stepCount_ >> 10) % (long)stride);
        for (std::size_t i = off; i < nN && !bad; i += stride)
            if (!std::isfinite(u_[i].x()) || !std::isfinite(u_[i].y())
                || !std::isfinite(u_[i].z()))
                bad = true;
        if (bad)
            throw std::runtime_error("FDEM3D instability (NaN)");
        checkEnergyAbort();                // opt-in (budgetAbortPct), E2
    }
}

// ---------------------------------------------------------------------------
// E2 (fiabilite) : moniteur d'energie runtime — l'« energy sanity abort » des
// codes de production. Opt-in par budgetAbortPct (0 = off, defaut) : si le
// residu B4 courant depasse budgetAbortPct % de l'echelle (meme definition
// que le resume), on arrete PROPREMENT via finished() : derniere frame,
// derniere ligne d'history et summary sont ecrits — un run qui diverge laisse
// ainsi son autopsie au lieu de 2,4 MJ de debris (cas gele du 2026-08-07).
// ---------------------------------------------------------------------------
void Fdem3dSolver::checkEnergyAbort() {
    if (eAbortPct_ < 0.0) {                // lecture paresseuse, une fois
        eAbortPct_ = cfg_.getd("budgetAbortPct", 0.0);
        // plancher ABSOLU [J] : en quasi-statique le flux part de zero et le
        // critere relatif declenche sur le bruit de fermeture des rampes
        // (mesure : 1,2 mJ = 7,4 % d'une echelle de 8 mJ, hotspot 0,05 m/s)
        eAbortMin_ = cfg_.getd("budgetAbortMin", 0.0);
    }
    if (eAbortPct_ <= 0.0 || eAbort_) return;
    double ke = 0.0;
    for (std::size_t i = 0; i < X0_.size(); ++i)
        ke += 0.5 * m_[i] * v_[i].squaredNorm();
    double sumW = elWork_ + jointWork_ + gcWork_ + cundWork_ + lysWork_
                + toolWork_ + bcWork_ + confWork_ + biasW_;
    double gross = std::abs(elWork_) + std::abs(jointWork_)
                 + std::abs(gcWork_) + std::abs(cundWork_)
                 + std::abs(lysWork_) + std::abs(toolWork_)
                 + std::abs(bcWork_) + std::abs(confWork_);
    double scale = std::max({keInit_, ke, gross, 1e-30});
    if (scale < 1e-12) return;             // charge nulle : pas de verdict
    double resid = (ke - keInit_) - sumW;
    if (std::abs(resid) <= 0.01 * eAbortPct_ * scale
        || std::abs(resid) <= eAbortMin_) return;
    int iw = 0; double vw = 0.0;           // hotspot : le noeud le plus rapide
    for (std::size_t i = 0; i < X0_.size(); ++i) {
        double vn = v_[i].squaredNorm();
        if (vn > vw) { vw = vn; iw = (int)i; }
    }
    std::cout << "[FDEM3D] ENERGY ABORT (budgetAbortPct = " << eAbortPct_
              << ") a t = " << t_ << " s : residu B4 " << resid << " J = "
              << 100.0 * std::abs(resid) / scale
              << " % de l'echelle. Hotspot : noeud " << iw << ", |v| = "
              << std::sqrt(vw) << " m/s, position ("
              << (X0_[iw] + u_[iw]).transpose() << ")\n";
    eAbort_ = true;
}

bool Fdem3dSolver::finished() const { return eAbort_ || peakStop_; }

// Co-rotational linear tet. R from F by Higham iteration R <- (R + R^-T)/2,
// warm-started from F scaled to unit Frobenius norm of a rotation. The
// Lame constants, crush cap and mean-tension cap are those of the
// element's PHASE (single-phase tables reduce to the pre-GBM arithmetic).
// Node duplication makes the loop embarrassingly parallel, as in 2D:
// every element writes only its OWN four nodes.
void Fdem3dSolver::elementForces() {
    double wEl = 0.0;                      // V2/B4 : travail des forces
                                           // internes ce pas (compteur pur)
    double wVi = 0.0;                      // dont part VISQUEUSE
    double wBd = 0.0;                      // dont part PULVERISATION (WP1)
    long nPv = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:wEl,wVi,wBd,nPv)
#endif
    for (int eI = 0; eI < (int)el_.size(); ++eI) {
        Elem& e = el_[eI];
        double lam = lamP_[e.phase];
        double mu2 = mu2P_[e.phase];
        Eigen::Matrix3d F = Eigen::Matrix3d::Zero();
        for (int a = 0; a < 4; ++a)
            F += (X0_[e.n[a]] + u_[e.n[a]]) * e.dN.col(a).transpose();
        Eigen::Matrix3d R;
        double det = F.determinant();
        if (det > 1e-9) {
            R = F * std::sqrt(3.0) / F.norm();
            for (int it = 0; it < 3; ++it)
                R = 0.5 * (R + R.inverse().transpose());
        } else R.setIdentity();               // degenerate: crush cap handles it
        Eigen::Matrix3d Ub = R.transpose() * F;
        Eigen::Matrix3d eps = 0.5 * (Ub + Ub.transpose()) - Eigen::Matrix3d::Identity();
        double tr = eps.trace();
        Eigen::Matrix3d sig;
        if (law_) sig = law_->stress(eps, e.st, dt_, hEl_[eI]);
        else      sig = lam * tr * Eigen::Matrix3d::Identity() + mu2 * eps;
        double pm = sig.trace() / 3.0;
        Eigen::Matrix3d dev = sig - pm * Eigen::Matrix3d::Identity();
        e.svm = std::sqrt(1.5) * dev.norm();
        if (!law_ && e.svm > crushCapP_[e.phase]) {   // deviatoric crush cap
            dev *= crushCapP_[e.phase] / e.svm;
            e.svm = crushCapP_[e.phase];
        }
        if (mtCap_ > 0.0 && pm > mtCap_ * ftP_[e.phase])
            pm = mtCap_ * ftP_[e.phase];              // mean-tension cap
        sig = dev + pm * Eigen::Matrix3d::Identity();
        // ---- WP1 : pulverisation (Yang et al. 2026, eq. 3-4) ------------
        // delta_m = h_e * eps_vm (deformation equivalente deviatorique) ;
        // D = delta_f (dm_max - delta_0) / (dm_max (delta_f - delta_0)),
        // irreversible via le max historique, plafonne a Dmax ; puis
        // sigma <- Cd (1 - D) sigma. Applique APRES les caps : le deck
        // choisit (principe VIII) — le granite de l article neutralise le
        // crushCap (1e12) et laisse ce modele seul degrader.
        if (bdOn_ && !law_) {
            Eigen::Matrix3d ed = eps
                - (eps.trace() / 3.0) * Eigen::Matrix3d::Identity();
            double dm = hEl_[eI] * std::sqrt(2.0 / 3.0) * ed.norm();
            if (dm > e.bdDm) e.bdDm = dm;
            if (e.bdDm > bdD0_) {
                double D = bdDf_ * (e.bdDm - bdD0_)
                         / (e.bdDm * (bdDf_ - bdD0_));
                if (D > bdDmax_) D = bdDmax_;
                // Dissipation d endommagement : Phi = Y dD, avec le taux de
                // restitution Y = 1/2 Cd eps : C : eps = 1/2 Cd sig_b : eps
                // (psi = 1/2 Cd (1-D) eps:C:eps, Y = -dpsi/dD). NULLE quand
                // D n evolue pas — un endommagement fige est elastique.
                if (D > e.bdD)
                    wBd -= 0.5 * bdCd_
                         * (sig.array() * eps.array()).sum()
                         * e.V0 * (D - e.bdD);
                e.bdD = D;
                double k = bdCd_ * (1.0 - D);
                sig *= k;
                e.svm *= k;
                if (D >= bdDmax_) ++nPv;
            }
        }
        Eigen::Matrix3d P = R * sig;
        e.sigG = P * R.transpose();        // global Cauchy (insertion sweep)
        // ---- viscosite newtonienne 2 mu D (Yan eq. 6) et taux pour le DIF
        // Ce n est PAS une  viscosite de volume  au sens zeta tr(D) I : le
        // terme agit sur le tenseur COMPLET, trace comprise. En 3D sig est
        // deja une Matrix3d : un seul geste, sinon sigma_zz, sigma_xz et
        // sigma_yz resteraient sans terme visqueux — contrainte ajoutee non
        // objective, et dissipation reelle differente de la ventilation.
        if ((viscOn_ || difOn_) && det > 1e-9) {
            Eigen::Matrix3d Fd = Eigen::Matrix3d::Zero();
            for (int a = 0; a < 4; ++a)
                Fd += v_[e.n[a]] * e.dN.col(a).transpose();
            Eigen::Matrix3d Lv = Fd * F.inverse();
            Eigen::Matrix3d Dr = 0.5 * (Lv + Lv.transpose());
            Eigen::Matrix3d Dc = R.transpose() * Dr * R;   // co-rote
            if (difOn_)
                e.edot = srRelax_ * e.edot
                       + (1.0 - srRelax_) * rockim::maxAbsEigSym3(Dc);
            if (viscOn_) {
                double mue = muEl_[eI];
                sig += 2.0 * mue * Dc;
                // Puissance 2 mu D:D, >= 0 par construction, comptee
                // NEGATIVEMENT comme tout ce qui quitte l energie cinetique.
                // Deja incluse dans wEl : VENTILATION, pas un poste de plus.
                // Nuance : 2 mu D:D est une densite par volume COURANT et V0
                // est le volume de REFERENCE — sous un insert ou det F tombe
                // a 0,5-0,7, cette ligne SOUS-ESTIME la dissipation reelle.
                wVi -= 2.0 * mue * Dc.squaredNorm() * e.V0;
                P = R * sig;                   // forces AVEC le visqueux
                if (viscIns_) e.sigG = P * R.transpose();
            }
        }
        for (int a = 0; a < 4; ++a) {
            Eigen::Vector3d fe = e.V0 * (P * e.dN.col(a));
            f_[e.n[a]] -= fe;
            wEl -= fe.dot(v_[e.n[a]]);     // V2/B4 : travail (lecture pure)
        }
    }
    elWork_ += wEl * dt_;
    viscWork_ += wVi * dt_;                // ventilation, incluse dans elWork_
    bdWork_ += wBd;                        // WP1 : deja une ENERGIE (Y dD)
    nPulv_ = nPv;
}

// Triangular cohesive joints, 3 node-pair integration points (A0/3 each).
// Same constitutive scheme as 2D: damage envelope in mode I, damage-plastic
// friction in mode II with a VECTOR slip return mapping in the face plane,
// shared scalar damage, death by clear separation only. All properties are
// PER JOINT (GBM). Parallel as in 2D: per-joint state is private, the
// force scatter goes through per-thread buffers reduced in thread order.
void Fdem3dSolver::jointForces() {
    auto processJoint = [&](Joint& J, auto&& addF, long& nb, double& dampW,
                            double& jw) {
        if (J.dead || J.bonded) return;    // bonded: node binding carries it
        Eigen::Vector3d P1 = 0.5 * ((X0_[J.a[0]] + u_[J.a[0]]) + (X0_[J.b[0]] + u_[J.b[0]]));
        Eigen::Vector3d P2 = 0.5 * ((X0_[J.a[1]] + u_[J.a[1]]) + (X0_[J.b[1]] + u_[J.b[1]]));
        Eigen::Vector3d P3 = 0.5 * ((X0_[J.a[2]] + u_[J.a[2]]) + (X0_[J.b[2]] + u_[J.b[2]]));
        Eigen::Vector3d nr = (P2 - P1).cross(P3 - P1);
        double nn = nr.norm();
        if (nn < 1e-18) return;
        Eigen::Vector3d n = nr / nn;                   // outward from A

        double At = J.A0 / 3.0;
        double dnMax = -1e30;
        double rsMaxO = 0.0;               // moteur de mode II du pas courant
                                           // (jointShearUnload = origin)

        for (int k = 0; k < 3; ++k) {
            int ia = J.a[k], ib = J.b[k];
            Eigen::Vector3d delta = (X0_[ib] + u_[ib]) - (X0_[ia] + u_[ia]);
            // dn0: adaptive-insertion opening offset (0 for intrinsic
            // joints) — a joint born under tension starts AT the envelope
            // peak, stress continuity as in the 2D solver
            double dn = delta.dot(n) + J.dn0;
            Eigen::Vector3d dt3 = delta - delta.dot(n) * n;
            dnMax = std::max(dnMax, dn);

            // ---- eq. 18 (jointShearUnload = origin) : miroir exact du 2D ----
            // s_max est mis a jour AVANT la traction normale : en mode
            // `origin` le moteur de mode II de l'eq. 16 en depend, et ce
            // moteur alimente D, donc f(D), donc l'enveloppe normale.
            // sEff = glissement VECTORIEL depuis l'origine figee J.slip[k]
            // (0 en intrinseque, -tau0/pj stampee a l'insertion : continuite
            // en cisaillement). sE = s_p de Munjiza, evalue sur l'enveloppe
            // de Mohr-Coulomb NON endommagee et sur la part geometrique pj*dn
            // de la contrainte normale — non circulaire, independant de D.
            Eigen::Vector3d sEff = Eigen::Vector3d::Zero();
            double smx = 0.0, rsO = 0.0;
            if (shearOrigin_) {
                J.slip[k] -= J.slip[k].dot(n) * n;     // origine dans le plan
                sEff = dt3 - J.slip[k];
                double sm = sEff.norm();
                if (sm > J.smax[k]) J.smax[k] = sm;
                smx = J.smax[k];
                double sE = (J.coh + J.tanPhi
                             * rockim::mcFrictionTerm(J.pj * dn, J.ft,
                                                      yangEnv_)) / J.pj;
                if (sE < 0.0) sE = 0.0;
                rsO = (J.slipF > 0.0 && smx > sE) ? (smx - sE) / J.slipF : 0.0;
                rsMaxO = std::max(rsMaxO, rsO);
            }

            // ================= normal traction =========================
            // The 2026-08-05 2D refonte, ported verbatim (FdemSolver.cpp):
            //  (1) the ELASTIC/cohesive part is computed FIRST — a viscous
            //      term must never be able to break a joint;
            //  (2) an INTACT cohesive joint is BILATERAL (it is a bond, it
            //      may pull); clipping belongs to a BROKEN joint acting as a
            //      contact, and then it applies to the RESULTANT;
            //  (3) the dashpot coefficient carries the MOOSE hard bound
            //      cd <= m_eff/dt (past it the dashpot reverses the approach
            //      velocity within one step instead of damping it) — the
            //      bound the 3D port had LOST (suspect #1 of the frozen
            //      debris-phase instability, 2026-08-07).
            double sigEl;                              // elastic / cohesive
            if (yanSoft_) {
                // exponential softening of Yan et al., eq. 11-13/16-17, as
                // in 2D: rn measured beyond the elastic branch dnE, mixed
                // driver D = sqrt(rn^2 + rs^2), origin-secant unloading.
                double ot = J.dnF - J.dnE;
                double rn = (ot > 0.0 && dn > J.dnE) ? (dn - J.dnE) / ot : 0.0;
                double rs = shearOrigin_
                    ? rsO
                    : ((J.slipF > 0.0) ? J.slip[k].norm() / J.slipF : 0.0);
                double Dnow = std::sqrt(rn * rn + rs * rs);
                if (Dnow > J.D) J.D = std::min(1.0, Dnow);
                double fdY = yan::fD(J.D, yanP_);
                if (dn >= 0.0) {
                    if (dn > J.omax[k]) J.omax[k] = dn;
                    double om = J.omax[k];
                    double sMax = std::min(J.pj * om, fdY * J.ft);
                    sigEl = (om > 1e-30) ? sMax * dn / om : 0.0;
                } else {
                    sigEl = (jcAdaptive_ ? (1.0 - J.D) * J.pj : J.pj) * dn;
                }
            } else if (dn >= 0.0) {
                double env = (dn <= J.dnE) ? J.pj * dn
                           : (dn >= J.dnF) ? 0.0
                           : J.ft * (J.dnF - dn) / (J.dnF - J.dnE);
                double tr2 = (1.0 - J.D) * J.pj * dn;
                if (tr2 > env) {
                    sigEl = env;
                    if (dn > 1e-30) {
                        double Dn = 1.0 - env / (J.pj * dn);
                        if (Dn > J.D) J.D = std::min(1.0, Dn);
                    }
                } else sigEl = tr2;
            } else {
                sigEl = (jcAdaptive_ ? (1.0 - J.D) * J.pj : J.pj) * dn;
            }

            double sig = sigEl;
            if (xiJ_ > 0.0) {
                Eigen::Vector3d vrel = v_[ib] - v_[ia];
                double meff = 0.5 * std::min(m_[ia], m_[ib]);
                double cd = 2.0 * xiJ_ * std::sqrt(J.pj * At * meff);
                // hard bound: past m_eff/dt the dashpot reverses the
                // approach velocity within one step (MOOSE), as in 2D
                cd = std::min(cd, meff / dt_);
                // SIGN: the traction acts on B as -sig*At*n, so a term that
                // OPPOSES the opening rate vrel.n carries a PLUS sign (the
                // historical minus was anti-damping — 2D lesson 2026-08-05)
                double sigV = cd * vrel.dot(n) / At;
                if (J.D < 1.0) {
                    sig = sigEl + sigV;                // bilateral bond
                } else if (dn < 0.0) {
                    sig = std::min(0.0, sigEl + sigV); // contact: clip the SUM
                }
                // power of the viscous part: F_B.v_B + F_A.v_A
                dampW -= (sig - sigEl) * At * vrel.dot(n) * dt_;
            }

            // yan: cohesion scaled by f(D); Coulomb term unscaled by default
            // so a crushed joint keeps residual friction (jointFrictionScaled
            // = 1 for the literal eq. 10), exactly as in 2D
            double fdS = yanSoft_ ? yan::fD(J.D, yanP_) : 0.0;
            double coh = yanSoft_ ? fdS * J.coh : (1.0 - J.D) * J.coh;
            double muS = (yanSoft_ && yanFricScaled_) ? fdS : 1.0;
            double tauLim = coh + muS * J.tanPhi
                          * rockim::mcFrictionTerm(sig, J.ft, yangEnv_);
            if (tauLim < 0.0) tauLim = 0.0;
            Eigen::Vector3d tau;
            if (shearOrigin_) {
                // ---- eq. 18 : secante a l'origine, version vectorielle ----
                // enveloppe = min(branche elastique, cap) au glissement
                // maximal ; l'etat courant est lu sur la droite origine ->
                // (s_max, tau_env). La DIRECTION est celle de sEff, comme le
                // retour radial prend celle de tauTr.
                double tauEnv = std::min(J.pj * smx, tauLim);
                tau.setZero();
                if (smx > 1e-30) tau = (tauEnv / smx) * sEff;
                // en `yan`, D a deja ete mis a jour par l'eq. 16 au-dessus
                if (!yanSoft_ && rsO > J.D) J.D = std::min(1.0, rsO);
            } else {
                J.slip[k] -= J.slip[k].dot(n) * n;     // keep slip in-plane
                Eigen::Vector3d tauTr = J.pj * (dt3 - J.slip[k]);
                double tn = tauTr.norm();
                tau = tauTr;
                if (tn > tauLim && tn > 0) {
                    tau *= tauLim / tn;
                    J.slip[k] += (tauTr - tau) / J.pj; // vector return mapping
                    double Dt2;
                    if (yanSoft_) {
                        double ot = J.dnF - J.dnE;
                        double rn = (ot > 0.0 && dn > J.dnE) ? (dn - J.dnE) / ot
                                                             : 0.0;
                        double rs = J.slip[k].norm() / J.slipF;
                        Dt2 = std::sqrt(rn * rn + rs * rs);   // eq. 16 / 14
                    } else {
                        Dt2 = J.slip[k].norm() / J.slipF;
                    }
                    if (Dt2 > J.D) J.D = std::min(1.0, Dt2);
                }
            }

            Eigen::Vector3d trac = (sig * n + tau) * At;
            addF(ib, -trac);
            addF(ia, trac);
            // V2/B4 : travail TOTAL des tractions de joint (visqueux inclus,
            // deja isole dans dampW) — lecture pure des vitesses
            jw += trac.dot(v_[ia] - v_[ib]) * dt_;
        }

        if (J.D >= 1.0) {
            if (J.tBreak < 0) {
                J.tBreak = t_; ++nb;
                // partition of the damage driver at the breaking instant
                // (ported from 2D): rn = mode I (opening), rs = mode II
                // (sliding). Output only — no force depends on it.
                double otF = J.dnF - J.dnE;
                double rnF = (otF > 0.0 && dnMax > J.dnE)
                           ? (dnMax - J.dnE) / otF : 0.0;
                double rsF;
                if (shearOrigin_) {            // pas de glissement plastique
                    rsF = rsMaxO;
                } else {
                    double sMx = 0.0;
                    for (int q = 0; q < 3; ++q)
                        sMx = std::max(sMx, J.slip[q].norm());
                    rsF = (J.slipF > 0.0) ? sMx / J.slipF : 0.0;
                }
                double den = rnF * rnF + rsF * rsF;
                J.failMode = (den > 1e-300) ? (rnF * rnF) / den : 1.0;
                J.rnB = rnF;
                J.rsB = rsF;
                J.bmode = (rnF >= rsF) ? 1 : 2;  // 1 traction, 2 cisaillement
            }
            if (dnMax > 3.0 * J.dnF) J.dead = true;    // separation only
        }
    };

#ifdef _OPENMP
    int nT = omp_get_max_threads();
    if (nT > 1) {
        if ((int)fTL_.size() != nT) {
            fTL_.assign(nT, std::vector<Eigen::Vector3d>(
                                X0_.size(), Eigen::Vector3d::Zero()));
            seenTL_.assign(nT, std::vector<char>(X0_.size(), 0));
            touchedTL_.assign(nT, {});
        }
        std::vector<long> nbT(nT, 0);
        std::vector<double> dwT(nT, 0.0), jwT(nT, 0.0);
#pragma omp parallel
        {
            int t = omp_get_thread_num();
            auto& fb = fTL_[t];
            auto& seen = seenTL_[t];
            auto& tl = touchedTL_[t];
            tl.clear();
            long nb = 0;
            double dw = 0.0, jw = 0.0;
            auto addF = [&](int i, const Eigen::Vector3d& v3) {
                if (!seen[i]) { seen[i] = 1; tl.push_back(i); }
                fb[i] += v3;
            };
#pragma omp for schedule(static)
            for (int jI = 0; jI < (int)jt_.size(); ++jI)
                processJoint(jt_[jI], addF, nb, dw, jw);
            nbT[t] = nb;
            dwT[t] = dw;
            jwT[t] = jw;
        }
        for (int t = 0; t < nT; ++t) {
            for (int i : touchedTL_[t]) {
                f_[i] += fTL_[t][i];
                fTL_[t][i].setZero();
                seenTL_[t][i] = 0;
            }
            nBroken_ += nbT[t];
            dampWork_ += dwT[t];
            jointWork_ += jwT[t];
        }
        return;
    }
#endif
    long nb1 = 0;                        // 1 thread: bit-identical to serial
    double dw1 = 0.0, jw1 = 0.0;
    auto addF1 = [&](int i, const Eigen::Vector3d& v3) { f_[i] += v3; };
    for (auto& J : jt_) processJoint(J, addF1, nb1, dw1, jw1);
    nBroken_ += nb1;
    dampWork_ += dw1;
    jointWork_ += jw1;
}

void Fdem3dSolver::rebuildContactFaces() {
    if (!poolBuilt_) {                     // le pool est FIXE (pas de fenetre)
        pool_ = exterior_;
        extOn_.assign(pool_.size(), gcAdaptive_ ? 0 : 1);
        poolBuilt_ = true;
    }
    act_.clear();
    if (gcAdaptive_) {                    // sous-ensemble ACTIF, ordre du pool
        for (std::size_t k = 0; k < pool_.size(); ++k)
            if (extOn_[k]) act_.push_back(pool_[k]);
    } else {
        act_.insert(act_.end(), pool_.begin(), pool_.end());
    }
    // faces liberees : depuis le CACHE, rafraichi seulement sur le
    // declencheur historique — voir le commentaire du 2D (un joint peut
    // mourir sans casser ; avancer son entree ferait diverger les modes)
    act_.insert(act_.end(), deadList_.begin(), deadList_.end());
    haveDead_ = !deadList_.empty();
    std::vector<char> inAct(X0_.size(), 0);
    for (const auto& bf : act_)
        for (int nid : bf.n) inAct[nid] = 1;
    actNodes_.clear();
    for (int i = 0; i < (int)X0_.size(); ++i)
        if (inAct[i]) actNodes_.push_back(i);
}

// General node-triangle contact with every 2D-learned safeguard: clipped
// grid box, co-location exclusion by origin vertex, birth-gap relief and a
// quasi-plastic soft normal law. On the voronoi mesh, faces can span
// several grid cells, so they are binned into ALL cells their AABB covers
// (with per-node pair dedup) and the deep-penetration cap uses the local
// element size — the 2D general-contact note, applied in 3D. The grid mesh
// keeps the original midpoint binning, bit-identical.
// gcActivation = adaptive — le balayage d'activation, miroir exact du 2D
// (regles C/A/B, cadence par v_max — voir FdemSolver.hpp / .cpp).
void Fdem3dSolver::activationSweep() {
    if (bodyStamp_ != nBroken_) {
        std::vector<int> uf(el_.size());
        for (int e2 = 0; e2 < (int)el_.size(); ++e2) uf[e2] = e2;
        auto find = [&](int x) {
            while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
            return x;
        };
        elemDam_.assign(el_.size(), 0);
        for (const auto& J : jt_) {
            if (J.bonded || (!J.dead && J.D < 1.0)) {
                int a = find(J.eA), b = find(J.eB);
                if (a != b) uf[a] = b;
            } else {
                elemDam_[J.eA] = elemDam_[J.eB] = 1;   // regle C
            }
        }
        bodyOf_.resize(el_.size());
        nBodies_ = 0;
        for (int e2 = 0; e2 < (int)el_.size(); ++e2) {
            bodyOf_[e2] = find(e2);
            if (uf[e2] == e2) ++nBodies_;
        }
        // regle C etendue d'un ANNEAU par sommet (voir le 2D : le coin du
        // bourrelet de cratere se deforme avant de rien casser lui-meme)
        if (vElems_.empty()) {
            int nV = 0;
            for (int v : vOf_) nV = std::max(nV, v + 1);
            vElems_.assign(nV, {});
            for (int nd = 0; nd < (int)vOf_.size(); ++nd)
                vElems_[vOf_[nd]].push_back(elemOf_[nd]);
        }
        std::vector<char> ring(el_.size(), 0);
        for (int e2 = 0; e2 < (int)el_.size(); ++e2)
            if (elemDam_[e2])
                for (int a = 0; a < 4; ++a)
                    for (int e3 : vElems_[vOf_[el_[e2].n[a]]]) ring[e3] = 1;
        for (int e2 = 0; e2 < (int)el_.size(); ++e2)
            if (ring[e2]) elemDam_[e2] = 1;
        bodyStamp_ = nBroken_;
    }
    if (nBodies_ <= 1 && nActivated_ == 0 && !haveDead_) {
        nextSweep_ = stepCount_ + gcActEvery_;
        return;
    }
    double cl = voronoi_ ? cellV_ : 2.0 * hmin_;       // = cell_ de detect()
    double M = gcActMargin_ * cl;
    struct SFace {
        Eigen::Vector3d c;
        double r;
        int body, idx;
        char on, touched;
    };
    static std::vector<SFace> sf;
    sf.clear();
    auto push = [&](const std::array<int, 3>& n3, int elem, char on, int idx) {
        Eigen::Vector3d A = X0_[n3[0]] + u_[n3[0]];
        Eigen::Vector3d B = X0_[n3[1]] + u_[n3[1]];
        Eigen::Vector3d C = X0_[n3[2]] + u_[n3[2]];
        Eigen::Vector3d c = (A + B + C) / 3.0;
        double r = std::sqrt(std::max({(A - c).squaredNorm(),
                                       (B - c).squaredNorm(),
                                       (C - c).squaredNorm()}));
        char tch = on && (lastTouch_[n3[0]] >= 0 || lastTouch_[n3[1]] >= 0
                          || lastTouch_[n3[2]] >= 0);
        sf.push_back({c, r, bodyOf_[elem], idx, on, tch});
    };
    for (std::size_t k = 0; k < pool_.size(); ++k)
        push(pool_[k].n, pool_[k].elem, extOn_[k], (int)k);
    for (const auto& J : jt_)
        if (J.dead) {
            push(J.a, J.eA, 1, -1);
            push({J.b[0], J.b[2], J.b[1]}, J.eB, 1, -1);
        }
    double rmax = 0.0;
    for (const auto& s : sf) rmax = std::max(rmax, s.r);
    double cs = M + 2.0 * rmax + 1e-300;
    static std::unordered_map<uint64_t, std::vector<int>> hg;
    hg.clear();
    auto keyIJK = [&](long long ix, long long iy, long long iz) {
        return (uint64_t)(ix & 0x1FFFFF) | ((uint64_t)(iy & 0x1FFFFF) << 21)
               | ((uint64_t)(iz & 0x1FFFFF) << 42);
    };
    auto cellOf = [&](const Eigen::Vector3d& p, long long& ix, long long& iy,
                      long long& iz) {
        ix = (long long)std::floor(p.x() / cs);
        iy = (long long)std::floor(p.y() / cs);
        iz = (long long)std::floor(p.z() / cs);
    };
    for (int q = 0; q < (int)sf.size(); ++q) {
        long long ix, iy, iz;
        cellOf(sf[q].c, ix, iy, iz);
        hg[keyIJK(ix, iy, iz)].push_back(q);
    }
    bool changed = false;
    for (int q = 0; q < (int)sf.size(); ++q) {
        SFace& f = sf[q];
        if (f.on || f.idx < 0) continue;
        if (elemDam_[pool_[f.idx].elem]) {             // regle C
            extOn_[f.idx] = 1;
            ++nActivated_;
            changed = true;
            continue;
        }
        long long ix, iy, iz;
        cellOf(f.c, ix, iy, iz);
        bool hit = false;
        for (long long dz = -1; dz <= 1 && !hit; ++dz)
            for (long long dy = -1; dy <= 1 && !hit; ++dy)
                for (long long dx = -1; dx <= 1 && !hit; ++dx) {
                    auto it = hg.find(keyIJK(ix + dx, iy + dy, iz + dz));
                    if (it == hg.end()) continue;
                    for (int g : it->second) {
                        if (g == q) continue;
                        const SFace& s = sf[g];
                        if (!(s.body != f.body || s.touched)) continue;
                        double gap = (s.c - f.c).norm() - s.r - f.r;
                        if (gap < M) { hit = true; break; }
                    }
                }
        if (hit) {
            extOn_[f.idx] = 1;
            ++nActivated_;
            changed = true;
        }
    }
    if (changed) rebuildContactFaces();    // recompose (cache mort INCHANGE)
    double v2 = 0.0;
    for (const auto& v : v_) v2 = std::max(v2, v.squaredNorm());
    double closing = 2.0 * std::sqrt(v2) * dt_;
    long kn = gcActEvery_;
    if (closing > 1e-300)
        kn = std::max((long)1,
                      std::min(gcActEvery_, (long)(0.5 * M / closing)));
    nextSweep_ = stepCount_ + kn;
}

// ---------------------------------------------------------------------------
// contact = potential — A3 phase 2, miroir exact du 2D (voir FdemSolver.cpp)
// en paires de TETS : detection O(N) par binning AABB, exclusion des paires
// liees par un joint vivant, force du polyedre de recouvrement
// (PotentialContact.hpp / pot3), releve de naissance par VOLUME, frottement
// incremental vectoriel. Serie et deterministe.
// ---------------------------------------------------------------------------
void Fdem3dSolver::potentialContact() {
    if (jointOfPair_.empty() && !jt_.empty()) {
        jointOfPair_.reserve(2 * jt_.size());
        for (int j = 0; j < (int)jt_.size(); ++j) {
            uint64_t a = (uint64_t)std::min(jt_[j].eA, jt_[j].eB);
            uint64_t b = (uint64_t)std::max(jt_[j].eA, jt_[j].eB);
            jointOfPair_[(a << 32) | b] = j;
        }
    }
    // ---- (1) elements uniques du jeu actif -------------------------------
    static std::vector<long> emark;
    static std::vector<int> elems;
    static long epoch = 0;
    if (emark.size() != el_.size()) emark.assign(el_.size(), -1);
    ++epoch;
    elems.clear();
    for (const auto& bf : act_)
        if (emark[bf.elem] != epoch) {
            emark[bf.elem] = epoch;
            elems.push_back(bf.elem);
        }
    if (elems.size() < 2) return;
    auto potTic = std::chrono::steady_clock::now();   // diagnostic (resume)

    // ---- (2) AABB courantes + grille DENSE a seaux REUTILISES (N1) -------
    // Miroir du 2D : la grille de hachage reallouait chaque seau a chaque
    // pas (la detection dominait le cout, x2,6 mesure sur la percussion 3D).
    // Boite fixe comme la grille de faces, seaux clear() sans liberation.
    double cl = voronoi_ ? cellV_ : 2.0 * hmin_;
    Eigen::Vector3d ebLo(-0.5 * W_, -0.5 * D_, -0.5 * H_);
    Eigen::Vector3d ebHi(1.5 * W_, 1.5 * D_, 2.0 * H_);
    Eigen::Vector3d egMin = ebLo - Eigen::Vector3d::Constant(cl);
    Eigen::Vector3d espan = ebHi - egMin + Eigen::Vector3d::Constant(cl);
    int egx = std::max(1, int(espan.x() / cl) + 1);
    int egy = std::max(1, int(espan.y() / cl) + 1);
    int egz = std::max(1, int(espan.z() / cl) + 1);
    auto ecid = [&](int cx, int cy, int cz) {
        return ((std::size_t)cz * egy + cy) * egx + cx;
    };
    static std::vector<Eigen::Vector3d> elo, ehi;
    elo.resize(elems.size());
    ehi.resize(elems.size());
    static std::vector<char> einb;
    einb.resize(elems.size());
    static std::vector<std::vector<int>> eg;
    {
        std::size_t nC = (std::size_t)egx * egy * egz;
        if (eg.size() != nC) eg.assign(nC, {});
        else for (auto& c : eg) c.clear();
    }
    for (int q = 0; q < (int)elems.size(); ++q) {
        const Elem& E = el_[elems[q]];
        Eigen::Vector3d lo = X0_[E.n[0]] + u_[E.n[0]], hi = lo;
        for (int a = 1; a < 4; ++a) {
            Eigen::Vector3d p = X0_[E.n[a]] + u_[E.n[a]];
            lo = lo.cwiseMin(p);
            hi = hi.cwiseMax(p);
        }
        elo[q] = lo;
        ehi[q] = hi;
        einb[q] = (hi.array() > ebLo.array()).all()
                  && (lo.array() < ebHi.array()).all();
        if (!einb[q]) continue;
        int x0 = std::clamp(int((lo.x() - egMin.x()) / cl), 0, egx - 1);
        int x1 = std::clamp(int((hi.x() - egMin.x()) / cl), 0, egx - 1);
        int y0 = std::clamp(int((lo.y() - egMin.y()) / cl), 0, egy - 1);
        int y1 = std::clamp(int((hi.y() - egMin.y()) / cl), 0, egy - 1);
        int z0 = std::clamp(int((lo.z() - egMin.z()) / cl), 0, egz - 1);
        int z1 = std::clamp(int((hi.z() - egMin.z()) / cl), 0, egz - 1);
        for (int cz = z0; cz <= z1; ++cz)
            for (int cy = y0; cy <= y1; ++cy)
                for (int cx = x0; cx <= x1; ++cx)
                    eg[ecid(cx, cy, cz)].push_back(q);
    }

    // ---- (3) paires candidates, en ordre CANONIQUE (voir 2D) -------------
    static std::vector<int> pstamp;
    if (pstamp.size() != elems.size()) pstamp.assign(elems.size(), -1);
    else std::fill(pstamp.begin(), pstamp.end(), -1);
    static std::vector<uint64_t> pairs;
    pairs.clear();
    for (int q = 0; q < (int)elems.size(); ++q) {
        if (!einb[q]) continue;
        int x0 = std::clamp(int((elo[q].x() - egMin.x()) / cl), 0, egx - 1);
        int x1 = std::clamp(int((ehi[q].x() - egMin.x()) / cl), 0, egx - 1);
        int y0 = std::clamp(int((elo[q].y() - egMin.y()) / cl), 0, egy - 1);
        int y1 = std::clamp(int((ehi[q].y() - egMin.y()) / cl), 0, egy - 1);
        int z0 = std::clamp(int((elo[q].z() - egMin.z()) / cl), 0, egz - 1);
        int z1 = std::clamp(int((ehi[q].z() - egMin.z()) / cl), 0, egz - 1);
        for (int cz = z0; cz <= z1; ++cz)
            for (int cy = y0; cy <= y1; ++cy)
                for (int cx = x0; cx <= x1; ++cx)
                    for (int r : eg[ecid(cx, cy, cz)]) {
                        if (r <= q || pstamp[r] == q) continue;
                        pstamp[r] = q;
                        if ((elo[r].array() > ehi[q].array()).any()
                            || (elo[q].array() > ehi[r].array()).any())
                            continue;
                        uint64_t a = (uint64_t)std::min(elems[q], elems[r]);
                        uint64_t b = (uint64_t)std::max(elems[q], elems[r]);
                        pairs.push_back((a << 32) | b);
                    }
    }
    std::sort(pairs.begin(), pairs.end());
    {
        auto t1 = std::chrono::steady_clock::now();
        potStats_.tGrid += std::chrono::duration<double>(t1 - potTic).count();
        potTic = t1;
    }
    potStats_.pairs += pairs.size();

    for (uint64_t pk : pairs) {
        {
            {
                        int eLo = (int)(pk >> 32);
                        int eHi = (int)(pk & 0xFFFFFFFFu);
                        auto itJ = jointOfPair_.find(pk);
                        if (itJ != jointOfPair_.end()
                            && !jt_[itJ->second].dead) {
                            ++potStats_.joint;
                            continue;      // le joint vivant porte la paire
                        }
                        const Elem& EA = el_[eLo];
                        const Elem& EB = el_[eHi];
                        pot3::V3 pa[4], pb[4];
                        for (int k = 0; k < 4; ++k) {
                            pa[k] = X0_[EA.n[k]] + u_[EA.n[k]];
                            pb[k] = X0_[EB.n[k]] + u_[EB.n[k]];
                        }
                        // pre-filtre d'axe separateur avec cache par paire :
                        // en regime etabli un voisin tangent coute UN test de
                        // plan au lieu d'un clip complet (bit-neutre — voir
                        // PotentialContact.hpp)
                        auto& H = potFt_[pk];
                        {
                            const int h0 = H.sepAxis;
                            if (pot3::separated(pa, pb, H.sepAxis)) {
                                if (H.sepAxis == h0) ++potStats_.sepHint;
                                else if (H.sepAxis < 8) ++potStats_.sepFace;
                                else ++potStats_.sepEdge;
                                continue;
                            }
                        }
                        pot3::PairForce3 R;
                        if (!pot3::pairForce(pa, pb, potP_, R)) {
                            ++potStats_.clipMiss;
                            continue;
                        }
                        ++potStats_.clipHit;
                        // ---- releve de naissance par VOLUME (cf. 2D) -----
                        // vRef < 0 = premiere fois en recouvrement (l'entree
                        // peut preexister via le cache d'axe separateur)
                        if (H.vRef < 0.0) H.vRef = R.vol;
                        else H.vRef *= relax_;
                        double sc = std::max(0.0, 1.0 - H.vRef / R.vol);
                        R.F *= sc;
                        for (int k = 0; k < 4; ++k) {
                            R.fA[k] *= sc;
                            R.fB[k] *= sc;
                        }
                        double w = 0.0;
                        for (int k = 0; k < 4; ++k) {
                            f_[EA.n[k]] += R.fA[k];
                            f_[EB.n[k]] += R.fB[k];
                            w += R.fA[k].dot(v_[EA.n[k]])
                               + R.fB[k].dot(v_[EB.n[k]]);
                        }
                        gcWork_ += w * dt_;
                        if (trackGroup_ >= 0) {        // V2/B2
                            if (elemGroup_[eLo] == trackGroup_) grpF_ += R.F;
                            if (elemGroup_[eHi] == trackGroup_) grpF_ -= R.F;
                        }
                        // barycentriques du centroide : estampilles PONDEREES
                        // de la regle B + frottement (voir le 2D — estampiller
                        // tous les noeuds sur-propageait l'activation : 96 %
                        // des faces activees mesurees sur la percussion longue)
                        pot3::Bary4 bA, bB;
                        bA.set(pa[0], pa[1], pa[2], pa[3]);
                        bB.set(pb[0], pb[1], pb[2], pb[3]);
                        double la[4], lb[4];
                        bA.lam(R.cen, la);
                        bB.lam(R.cen, lb);
                        if (gcAdaptive_)
                            for (int k = 0; k < 4; ++k) {
                                if (la[k] > 0.1)
                                    lastTouch_[EA.n[k]] = stepCount_;
                                if (lb[k] > 0.1)
                                    lastTouch_[EB.n[k]] = stepCount_;
                            }
                        // ---- frottement incremental (eq. 4-5) ------------
                        if (muC_ > 0.0 && potKt_ > 0.0) {
                            double Fn = R.F.norm();
                            if (Fn > 1e-300) {
                                Eigen::Vector3d vA =
                                    la[0] * v_[EA.n[0]] + la[1] * v_[EA.n[1]]
                                    + la[2] * v_[EA.n[2]]
                                    + la[3] * v_[EA.n[3]];
                                Eigen::Vector3d vB =
                                    lb[0] * v_[EB.n[0]] + lb[1] * v_[EB.n[1]]
                                    + lb[2] * v_[EB.n[2]]
                                    + lb[3] * v_[EB.n[3]];
                                Eigen::Vector3d vrel = vA - vB;
                                Eigen::Vector3d nh = R.F / Fn;
                                if (H.step < stepCount_ - 1) H.Ft.setZero();
                                Eigen::Vector3d Ft =
                                    H.Ft - H.Ft.dot(nh) * nh;
                                Eigen::Vector3d vt =
                                    vrel - vrel.dot(nh) * nh;
                                Ft -= potKt_ * dt_ * vt;
                                double cap = muC_ * Fn;
                                double Ftn = Ft.norm();
                                if (Ftn > cap && Ftn > 0.0)
                                    Ft *= cap / Ftn;
                                H.Ft = Ft;
                                H.step = stepCount_;
                                for (int k = 0; k < 4; ++k) {
                                    f_[EA.n[k]] += la[k] * Ft;
                                    f_[EB.n[k]] -= lb[k] * Ft;
                                }
                                gcWork_ += Ft.dot(vrel) * dt_;
                                gcFricWork_ += Ft.dot(vrel) * dt_;  // V2/B4
                                if (trackGroup_ >= 0) {  // V2/B2
                                    if (elemGroup_[eLo] == trackGroup_)
                                        grpF_ += Ft;
                                    if (elemGroup_[eHi] == trackGroup_)
                                        grpF_ -= Ft;
                                }
                            } else {
                                H.step = stepCount_;
                            }
                        } else {
                            H.step = stepCount_;
                        }
                    }
                }
    }
    potStats_.tLoop += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - potTic).count();
}

void Fdem3dSolver::generalContact() {
    grpF_.setZero();                       // V2/B2 : force du pas courant
    if (gcAdaptive_) {
        if (!poolBuilt_) rebuildContactFaces();
        if (lastTouch_.empty()) lastTouch_.assign(X0_.size(), -1);
        if (stepCount_ >= nextSweep_
            || (sweepBroken_ != nBroken_
                && (sweepBroken_ < 0 || stepCount_ % 8 == 0))) {
            activationSweep();             // peut recomposer (cache inchange)
            sweepBroken_ = nBroken_;
        }
    }
    if (actStamp_ != nBroken_ && (actStamp_ < 0 || stepCount_ % 8 == 0)) {
        deadList_.clear();                 // le declencheur historique : seul
        for (const auto& J : jt_)          // endroit ou le cache se rafraichit
            if (J.dead) {
                deadList_.push_back({J.eA, J.a});
                deadList_.push_back({J.eB, {J.b[0], J.b[2], J.b[1]}});
            }
        rebuildContactFaces();
        actStamp_ = nBroken_;
    }
    if (act_.empty()) return;
    if (contactPot_) {                     // A3 : contact par potentiel —
        potentialContact();                // meme jeu actif, autre physique
        return;
    }

    // Cell size: 2 hmin on the grid mesh (unchanged); on the voronoi mesh
    // hmin is the THINNEST sliver, and (L/hmin)^3 dense cells reallocated
    // every step cost ~0.5 GB/step on the percussion demo — the voronoi
    // path uses 2 x median element size and a SPARSE hash grid instead
    // (allocation proportional to active faces, not to the box volume).
    cell_ = voronoi_ ? cellV_ : 2.0 * hmin_;
    Eigen::Vector3d boxLo(-0.5 * W_, -0.5 * D_, -0.5 * H_);
    Eigen::Vector3d boxHi(1.5 * W_, 1.5 * D_, 2.0 * H_);
    std::vector<Eigen::Vector3d> cen(act_.size());
    std::vector<char> inBox(act_.size(), 0);
    for (std::size_t k = 0; k < act_.size(); ++k) {
        cen[k] = ((X0_[act_[k].n[0]] + u_[act_[k].n[0]])
                  + (X0_[act_[k].n[1]] + u_[act_[k].n[1]])
                  + (X0_[act_[k].n[2]] + u_[act_[k].n[2]])) / 3.0;
        inBox[k] = (cen[k].array() > boxLo.array()).all()
                   && (cen[k].array() < boxHi.array()).all();
    }
    gmin_ = boxLo - Eigen::Vector3d::Constant(cell_);
    Eigen::Vector3d span = boxHi - gmin_ + Eigen::Vector3d::Constant(cell_);
    gx_ = std::max(1, int(span.x() / cell_) + 1);
    gy_ = std::max(1, int(span.y() / cell_) + 1);
    gz_ = std::max(1, int(span.z() / cell_) + 1);
    auto cidx = [&](int cx, int cy, int cz) {
        return ((uint64_t)cz * gy_ + cy) * gx_ + cx;
    };
    if (voronoi_) {
        gridV_.clear();
        for (std::size_t k = 0; k < act_.size(); ++k) {
            if (!inBox[k]) continue;
            // AABB binning: a voronoi face can span several cells
            Eigen::Vector3d lo = X0_[act_[k].n[0]] + u_[act_[k].n[0]];
            Eigen::Vector3d hi = lo;
            for (int q = 1; q < 3; ++q) {
                Eigen::Vector3d p = X0_[act_[k].n[q]] + u_[act_[k].n[q]];
                lo = lo.cwiseMin(p);
                hi = hi.cwiseMax(p);
            }
            int x0 = std::clamp(int((lo.x() - gmin_.x()) / cell_), 0, gx_ - 1);
            int y0 = std::clamp(int((lo.y() - gmin_.y()) / cell_), 0, gy_ - 1);
            int z0 = std::clamp(int((lo.z() - gmin_.z()) / cell_), 0, gz_ - 1);
            int x1 = std::clamp(int((hi.x() - gmin_.x()) / cell_), 0, gx_ - 1);
            int y1 = std::clamp(int((hi.y() - gmin_.y()) / cell_), 0, gy_ - 1);
            int z1 = std::clamp(int((hi.z() - gmin_.z()) / cell_), 0, gz_ - 1);
            for (int cz = z0; cz <= z1; ++cz)
                for (int cy = y0; cy <= y1; ++cy)
                    for (int cx = x0; cx <= x1; ++cx)
                        gridV_[cidx(cx, cy, cz)].push_back((int)k);
        }
    } else {
        // reuse the buckets instead of destroying them: assign() frees every
        // inner vector every step, clear() keeps their capacity (bit-neutral)
        std::size_t nCells = (std::size_t)gx_ * gy_ * gz_;
        if (grid_.size() != nCells) grid_.assign(nCells, {});
        else for (auto& c : grid_) c.clear();
        for (std::size_t k = 0; k < act_.size(); ++k) {
            if (!inBox[k]) continue;
            int cx = std::clamp(int((cen[k].x() - gmin_.x()) / cell_), 0, gx_ - 1);
            int cy = std::clamp(int((cen[k].y() - gmin_.y()) / cell_), 0, gy_ - 1);
            int cz = std::clamp(int((cen[k].z() - gmin_.z()) / cell_), 0, gz_ - 1);
            grid_[cidx(cx, cy, cz)].push_back((int)k);
        }
    }

    double cap = 0.6 * hmin_;

    // ---- geometry of the ACTIVE FACES, computed ONCE per step --------------
    // The sweep below used to rebuild A/B/C, the face normal (a cross product
    // AND a sqrt) and the barycentric denominators for EVERY (node, face)
    // pair — i.e. dozens of times per face on a dense exterior, where the
    // face has not moved between two of those evaluations. Hoisting it out is
    // arithmetically NEUTRAL (same values, same pair list, same order) and
    // removes one sqrt plus ~10 dot products per candidate pair.
    struct FGeo {
        Eigen::Vector3d A, v0, v1, nrm, cen;
        double d00, d01, d11, den, capLoc, rad2;
        bool ok;
    };
    static std::vector<FGeo> fgeo;
    fgeo.resize(act_.size());
    for (std::size_t k = 0; k < act_.size(); ++k) {
        FGeo& g = fgeo[k];
        g.ok = false;
        if (!inBox[k]) continue;                  // never binned = never a candidate
        const BFace& bf = act_[k];
        g.A = X0_[bf.n[0]] + u_[bf.n[0]];
        Eigen::Vector3d B = X0_[bf.n[1]] + u_[bf.n[1]];
        Eigen::Vector3d C = X0_[bf.n[2]] + u_[bf.n[2]];
        Eigen::Vector3d nr = (B - g.A).cross(C - g.A);
        double a2 = nr.norm();
        if (a2 < 1e-18) continue;
        g.nrm = nr / a2;                          // outward of bf.elem
        g.v0 = B - g.A;
        g.v1 = C - g.A;
        g.d00 = g.v0.dot(g.v0);
        g.d01 = g.v0.dot(g.v1);
        g.d11 = g.v1.dot(g.v1);
        g.den = g.d00 * g.d11 - g.d01 * g.d01;
        if (g.den < 1e-24) continue;
        g.capLoc = voronoi_ ? 0.6 * hEl_[bf.elem] : cap;
        g.cen = cen[k];
        // EXACT culling sphere. An accepted contact projects INSIDE the
        // triangle (hence within rMax of its centroid) at a depth <= capLoc,
        // so an accepted node always lies within rMax + capLoc of the
        // centroid. Rejecting the others early therefore removes NO pair —
        // it is a bound, not an approximation.
        double rMax = std::max({(g.A - g.cen).norm(), (B - g.cen).norm(),
                                (C - g.cen).norm()});
        double rc = rMax + g.capLoc;
        g.rad2 = rc * rc;
        g.ok = true;
    }

    // Detection is the expensive part (grid sweep + geometry) and is PURE:
    // it parallelizes over nodes into per-thread candidate lists, exactly as
    // the 2D solver does (FdemSolver::generalContact). The DELICATE part —
    // birth-gap bookkeeping (pen0_), damping, force application and the
    // net-work meter — stays SERIAL, walking the candidates in thread order.
    // schedule(static) hands each thread a contiguous chunk of actNodes_, so
    // the concatenated candidate order equals the serial sweep order: the
    // result is BIT-IDENTICAL to the serial path for any thread count.
    struct CPair {
        int i, k;
        double pen, wa, wb, wc;
        Eigen::Vector3d nrm;
    };
    auto detect = [&](int i, std::vector<CPair>& outC,
                      std::vector<int>& stamp) {
        Eigen::Vector3d p = X0_[i] + u_[i];
        if ((p.array() <= boxLo.array()).any()
            || (p.array() >= boxHi.array()).any()) return;
        int ci = std::clamp(int((p.x() - gmin_.x()) / cell_), 0, gx_ - 1);
        int cj = std::clamp(int((p.y() - gmin_.y()) / cell_), 0, gy_ - 1);
        int ck = std::clamp(int((p.z() - gmin_.z()) / cell_), 0, gz_ - 1);
        int myElem = elemOf_[i];
        // no per-node clear of the dedup structure: the stamp below is keyed
        // on the node id, which is unique within one sweep
        for (int dk = -1; dk <= 1; ++dk)
        for (int dj = -1; dj <= 1; ++dj)
        for (int di = -1; di <= 1; ++di) {
            int cx = ci + di, cy = cj + dj, cz = ck + dk;
            if (cx < 0 || cy < 0 || cz < 0 || cx >= gx_ || cy >= gy_ || cz >= gz_)
                continue;
            const std::vector<int>* lst;
            if (voronoi_) {
                auto itc = gridV_.find(cidx(cx, cy, cz));
                if (itc == gridV_.end()) continue;
                lst = &itc->second;
            } else {
                lst = &grid_[cidx(cx, cy, cz)];
            }
            for (int k : *lst) {
                const BFace& bf = act_[k];
                if (bf.elem == myElem) continue;
                if (voronoi_) {
                    // AABB binning can list one face in several stencil
                    // cells: dedup per NODE. A stamp array is O(1) where the
                    // former linear scan was O(candidates) per candidate,
                    // i.e. quadratic on the crowded cells of a fine mesh.
                    if (stamp[k] == i) continue;
                    stamp[k] = i;
                }
                const FGeo& g = fgeo[k];
                if (!g.ok) continue;
                if ((p - g.cen).squaredNorm() > g.rad2) continue;   // exact bound
                bool colocated = false;
                for (int q = 0; q < 3; ++q)
                    if (vOf_[i] == vOf_[bf.n[q]]) {
                        Eigen::Vector3d Pq = X0_[bf.n[q]] + u_[bf.n[q]];
                        if ((p - Pq).norm() < 0.25 * hmin_) colocated = true;
                    }
                if (colocated) continue;
                double d = (p - g.A).dot(g.nrm);
                if (d >= 0.0 || d < -g.capLoc) continue;
                // barycentric inside test of the projection
                Eigen::Vector3d q = p - d * g.nrm;
                Eigen::Vector3d v2 = q - g.A;
                double d20 = v2.dot(g.v0), d21 = v2.dot(g.v1);
                double wb = (g.d11 * d20 - g.d01 * d21) / g.den;
                double wc = (g.d00 * d21 - g.d01 * d20) / g.den;
                double wa = 1.0 - wb - wc;
                if (wa < 0 || wb < 0 || wc < 0) continue;
                outC.push_back({i, k, -d, wa, wb, wc, g.nrm});
            }
        }
    };

    static std::vector<std::vector<CPair>> cpTL;       // per-thread lists
    bool par = false;
#ifdef _OPENMP
    // fork/join costs more than the sweep while few faces are active: go
    // parallel only when the debris population makes it genuinely heavy
    par = omp_get_max_threads() > 1 && actNodes_.size() >= 4096;
#endif
    if (par) {
#ifdef _OPENMP
        int nT = omp_get_max_threads();
        if ((int)cpTL.size() < nT) cpTL.assign(nT, {});
#pragma omp parallel
        {
            int t = omp_get_thread_num();
            cpTL[t].clear();
            std::vector<int> stamp(act_.size(), -1);   // per-node dedup stamp
#pragma omp for schedule(static)
            for (int a = 0; a < (int)actNodes_.size(); ++a)
                detect(actNodes_[a], cpTL[t], stamp);
        }
#endif
    } else {
        if (cpTL.empty()) cpTL.assign(1, {});
        for (auto& lst : cpTL) lst.clear();
        std::vector<int> stamp(act_.size(), -1);
        for (int i : actNodes_) detect(i, cpTL[0], stamp);
    }

    for (auto& lst : cpTL)
        for (const CPair& cp : lst) {
            int i = cp.i;
            const BFace& bf = act_[cp.k];
            const Eigen::Vector3d& nrm = cp.nrm;
            double pen = cp.pen;
            uint64_t pkey = (uint64_t(uint32_t(i)) << 40)
                            ^ (uint64_t(uint32_t(bf.n[0])) << 20)
                            ^ uint32_t(bf.n[1]);
            auto [it0, isNew] = pen0_.try_emplace(pkey, pen);
            if (!isNew) it0->second *= relax_;
            pen = std::max(0.0, pen - it0->second);
            if (pen <= 0.0) continue;
            Eigen::Vector3d vFace = cp.wa * v_[bf.n[0]] + cp.wb * v_[bf.n[1]]
                                    + cp.wc * v_[bf.n[2]];
            Eigen::Vector3d vrel = v_[i] - vFace;
            double vn = vrel.dot(nrm);
            double cdmp = 2.0 * xiGC_ * std::sqrt(kpGC_ * m_[i]);
            double fn = kpGC_ * pen * (vn < 0.0 ? 1.0 : gcRest_) - cdmp * vn;
            if (fn < 0) fn = 0;
            Eigen::Vector3d vt = vrel - vn * nrm;
            double vtn = vt.norm();
            Eigen::Vector3d ftv = Eigen::Vector3d::Zero();
            if (vtn > 0)
                ftv = -muC_ * fn * std::tanh(vtn / vReg_) * vt / vtn;
            Eigen::Vector3d Fc = fn * nrm + ftv;
            double capF = 20.0 * m_[i] / dt_;
            double Fn2 = Fc.norm();
            if (Fn2 > capF) Fc *= capF / Fn2;
            gcWork_ += Fc.dot(vrel) * dt_;
            // V2/B4 : part tangentielle (ftv est normal a nrm, le cap scale)
            gcFricWork_ += (Fc - Fc.dot(nrm) * nrm).dot(vrel) * dt_;
            f_[i] += Fc;
            f_[bf.n[0]] -= cp.wa * Fc;
            f_[bf.n[1]] -= cp.wb * Fc;
            f_[bf.n[2]] -= cp.wc * Fc;
            if (trackGroup_ >= 0) {        // V2/B2 (lecture pure)
                if (elemGroup_[elemOf_[i]] == trackGroup_) grpF_ += Fc;
                if (elemGroup_[bf.elem] == trackGroup_) grpF_ -= Fc;
            }
            if (gcAdaptive_) {                         // source de la regle B
                lastTouch_[i] = stepCount_;
                lastTouch_[bf.n[0]] = stepCount_;
                lastTouch_[bf.n[1]] = stepCount_;
                lastTouch_[bf.n[2]] = stepCount_;
            }
        }
}

// Rigid tool (sphere or flat punch) against every node. Per-node force
// writes are race-free; the tool reaction is reduced through per-thread
// partial sums added in thread order (deterministic per thread count).
void Fdem3dSolver::toolContact() {
    if (scen_ == Scenario::TENSION || toolNone_) return;
    // TRI : l'outil est RETIRE pendant la classification, sinon il continue de
    // fracturer pendant qu'on classe. C'est aussi le geste experimental.
    if (brushArmed_) return;
    // ETAPE 1 de Yang : arret de l'outil AVANT l'armement, pour menager un
    // intervalle de repos (voir Fdem3dSolver.hpp).
    if (toolStop_ > 0.0 && t_ >= toolStop_) {
        if (!toolStopped_) {
            toolStopped_ = true;
            std::cout << "[FDEM3D] OUTIL ARRETE a t = " << t_ << " s. Repos "
                         "jusqu'a l'armement du tri (t = " << brushStart_
                      << " s), soit " << (brushStart_ - t_) << " s.\n";
        }
        return;
    }
    auto nodeFc = [&](int i, Eigen::Vector3d& Fc) {
        Eigen::Vector3d p = X0_[i] + u_[i];
        if (tool_.flat) {
            // flat-ended punch: vertical contact only against the bottom
            // face (sharp edge, as in 2D)
            double rx = p.x() - tool_.x.x(), ry = p.y() - tool_.x.y();
            if (rx * rx + ry * ry > tool_.radius * tool_.radius) return false;
            double pen = p.z() - tool_.x.z();
            if (pen <= 0) return false;
            double c = 2.0 * xiC_ * std::sqrt(kp_ * m_[i]);
            double fn = kp_ * pen + c * (v_[i].z() - tool_.v.z());
            if (fn < 0) fn = 0;
            Eigen::Vector3d vt(v_[i].x() - tool_.v.x(),
                               v_[i].y() - tool_.v.y(), 0.0);
            double vtn = vt.norm();
            Eigen::Vector3d ftv = Eigen::Vector3d::Zero();
            if (vtn > 0) ftv = -muC_ * fn * std::tanh(vtn / vReg_) * vt / vtn;
            Fc = ftv + Eigen::Vector3d(0.0, 0.0, -fn);
        } else {
            Eigen::Vector3d d = p - tool_.x;
            double dist = d.norm();
            if (dist >= tool_.radius || dist < 1e-14) return false;
            Eigen::Vector3d n = d / dist;
            double pen = tool_.radius - dist;
            Eigen::Vector3d vrel = v_[i] - tool_.v;
            double c = 2.0 * xiC_ * std::sqrt(kp_ * m_[i]);
            double fn = kp_ * pen - c * vrel.dot(n);
            if (fn < 0) fn = 0;
            Eigen::Vector3d vt = vrel - vrel.dot(n) * n;
            double vtn = vt.norm();
            Eigen::Vector3d ftv = Eigen::Vector3d::Zero();
            if (vtn > 0) ftv = -muC_ * fn * std::tanh(vtn / vReg_) * vt / vtn;
            Fc = fn * n + ftv;
        }
        // --- ECRETAGE EN IMPULSION (porte du 2D) : voir Fdem3dSolver.hpp ----
        // Borne PHYSIQUE du choc elastique contre une masse infinie : un outil
        // a vitesse imposee ne decelere jamais, c'est un reservoir d'energie
        // infini. On ecrete le VECTEUR complet, frottement compris, car c'est
        // l'impulsion TOTALE qui accelere le noeud.
        if (toolVCap_ > 0.0) {
            double vtool = tool_.v.norm();
            if (vtool > 1e-12) {
                double fmax = toolVCap_ * 2.0 * vtool * m_[i] / dt_;
                double f = Fc.norm();
                if (f > fmax) Fc *= fmax / f;
            }
        }
        return true;
    };

    // -----------------------------------------------------------------------
    // A1 — variante SIGNORINI EN VITESSE, miroir exact du 2D.
    // Geometrie DUPLIQUEE a dessein : la voie penalite ne doit pas bouger.
    // Par noeud, avec H = 1/m_i (masse concentree + obstacle rigide) :
    //   v* = v_i + (dt/m_i) f_i ;  vn = (v* - v_outil).n
    //   gap predit g+ = -pen + dt vn >= 0  ->  r = 0 (le noeud se separe seul)
    //   sinon r_n = m_i (relax*pen/dt - vn) > 0, et frottement par cap de
    //   Coulomb sur l'impulsion tangentielle |r_t| <= mu r_n.
    // Report en force r/dt pour que integrate() (v += dt/m f) produise le saut
    // de vitesse voulu et que les compteurs d'energie restent comparables.
    // -----------------------------------------------------------------------
    auto nodeSig = [&](int i, Eigen::Vector3d& Fc) {
        Eigen::Vector3d p = X0_[i] + u_[i];
        Eigen::Vector3d n;
        double pen;
        if (tool_.flat) {
            double rx = p.x() - tool_.x.x(), ry = p.y() - tool_.x.y();
            if (rx * rx + ry * ry > tool_.radius * tool_.radius) return false;
            pen = p.z() - tool_.x.z();
            if (pen <= 0.0) return false;
            n = Eigen::Vector3d(0.0, 0.0, -1.0);   // repousse vers le bas
        } else {
            Eigen::Vector3d d = p - tool_.x;
            double dist = d.norm();
            if (dist >= tool_.radius || dist < 1e-14) return false;
            n = d / dist;
            pen = tool_.radius - dist;
        }
        Eigen::Vector3d vFree = v_[i] + (dt_ / m_[i]) * f_[i];
        Eigen::Vector3d vrel = vFree - tool_.v;
        double vn = vrel.dot(n);
        if (-pen + dt_ * vn >= 0.0) return false;
        double rn = m_[i] * (toolSigRelax_ * pen / dt_ - vn);
        if (rn <= 0.0) return false;
        // frottement : impulsion de collage, ecretee par le cap de Coulomb
        Eigen::Vector3d vt = vrel - vn * n;
        double vtn = vt.norm();
        Eigen::Vector3d rt = Eigen::Vector3d::Zero();
        if (vtn > 1e-14) {
            double rts = m_[i] * vtn;                  // colle exactement
            double cap = muC_ * rn;
            if (rts > cap) rts = cap;
            rt = -rts * vt / vtn;
        }
        Fc = (rn / dt_) * n + rt / dt_;
        return true;
    };
#ifdef _OPENMP
    int nT = omp_get_max_threads();
    std::vector<Eigen::Vector3d> FT(nT, Eigen::Vector3d::Zero());
    double tw = 0.0;                       // V2/B4 : travail outil -> solide
#pragma omp parallel reduction(+:tw)
    {
        int t = omp_get_thread_num();
        Eigen::Vector3d Floc = Eigen::Vector3d::Zero();
#pragma omp for schedule(static)
        for (int i = 0; i < (int)X0_.size(); ++i) {
            Eigen::Vector3d Fc;
            if (!(toolSig_ ? nodeSig(i, Fc) : nodeFc(i, Fc))) continue;
            f_[i] += Fc;
            tw += Fc.dot(v_[i]) * dt_;
            Floc -= Fc;
        }
        FT[t] = Floc;
    }
    for (const auto& F : FT) tool_.F += F;
    toolWork_ += tw;
#else
    for (int i = 0; i < (int)X0_.size(); ++i) {
        Eigen::Vector3d Fc;
        if (!(toolSig_ ? nodeSig(i, Fc) : nodeFc(i, Fc))) continue;
        f_[i] += Fc;
        toolWork_ += Fc.dot(v_[i]) * dt_;  // V2/B4 : outil -> solide
        tool_.F -= Fc;
    }
#endif
}

void Fdem3dSolver::integrate() {
    // V2/B4 : KE du solide au premier pas (les vitesses initiales — groupVel,
    // percussion — sont encore intactes avant la premiere integration)
    if (keInit_ < 0.0) {
        double ke0 = 0.0;
        for (std::size_t i = 0; i < X0_.size(); ++i)
            ke0 += 0.5 * m_[i] * v_[i].squaredNorm();
        keInit_ = ke0;
    }
    double cw = 0.0, lw = 0.0, bw = 0.0, bias = 0.0;  // V2/B4 compteurs
    if (adaptive_) {
        // Bound groups integrate as ONE node (sum of forces and masses,
        // damping and quiet-boundary terms on the sums), exactly as in the
        // 2D solver. Copies of a group share flags (same position) and stay
        // bit-identical: groups only ever split, never merge.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:cw,lw,bw,bias)
#endif
        for (int vv = 0; vv < nVert_; ++vv) {
            for (const auto& g : grpsOfVert_[vv]) {
                int i0 = g[0];
                if (flag_[i0] == FIXED) {
                    if (gripFree_) {
                        double Fx = 0.0, Fy = 0.0, M = 0.0;
                        for (int i : g) { Fx += f_[i].x(); Fy += f_[i].y(); M += m_[i]; }
                        double vx = v_[i0].x() + (dt_ / M) * Fx;
                        double vy = v_[i0].y() + (dt_ / M) * Fy;
                        for (int i : g) {
                            v_[i] = {vx, vy, 0.0};
                            u_[i] += dt_ * v_[i];
                        }
                    } else for (int i : g) v_[i].setZero();
                    continue;
                }
                if (flag_[i0] == PRESCRIBED) {
                    double vx = 0.0, vy = 0.0;
                    if (gripFree_) {
                        double Fx = 0.0, Fy = 0.0, M = 0.0;
                        for (int i : g) { Fx += f_[i].x(); Fy += f_[i].y(); M += m_[i]; }
                        vx = v_[i0].x() + (dt_ / M) * Fx;
                        vy = v_[i0].y() + (dt_ / M) * Fy;
                    }
                    double tv = t_ - pullDelay_;   // triaxial : sigma3 d'abord
                    double vg = tv <= 0.0 ? 0.0 : pullV_;
                    if (vg != 0.0 && pullRamp_ > 0.0 && tv < pullRamp_)
                        vg *= 0.5 * (1.0 - std::cos(M_PI * tv / pullRamp_));
                    for (int i : g) {
                        // V2/B4 (corrige 2026-08-14, formule Sierra) :
                        // travail de la LIAISON R = m.a - f en vitesse
                        // moyenne — les familles comptent deja f.v ici
                        double vgOld = v_[i].z();
                        double Rz = m_[i] * (vg - vgOld) / dt_ - f_[i].z();
                        bw += Rz * 0.5 * (vg + vgOld) * dt_;
                        v_[i] = {vx, vy, vg};
                        u_[i] += dt_ * v_[i];
                    }
                    continue;
                }
                Eigen::Vector3d F = Eigen::Vector3d::Zero();
                Eigen::Vector3d cS = Eigen::Vector3d::Zero();
                double M = 0.0;
                for (int i : g) {
                    F += f_[i];
                    for (int a = 0; a < 3; ++a)
                        if (kAbs_[i](a) > 0) {
                            double fk = kAbs_[i](a) * u_[i](a);
                            F(a) -= fk;
                            lw -= fk * v_[i0](a) * dt_;   // V2/B4 ressort
                        }
                    cS += cAbs_[i];
                    M += m_[i];
                }
                if (damping_ > 0)
                    for (int a = 0; a < 3; ++a) {
                        double fd = damping_ * std::abs(F(a))
                                * (v_[i0](a) > 0 ? 1.0 : (v_[i0](a) < 0 ? -1.0 : 0.0));
                        F(a) -= fd;
                        cw -= fd * v_[i0](a) * dt_;       // V2/B4 Cundall
                    }
                bias += F.squaredNorm() * dt_ * dt_ / (2.0 * M);
                Eigen::Vector3d vn = v_[i0] + (dt_ / M) * F;
                for (int a = 0; a < 3; ++a)
                    if (cS(a) > 0) {
                        vn(a) /= 1.0 + dt_ * cS(a) / M;
                        lw -= cS(a) * vn(a) * vn(a) * dt_;  // V2/B4 amortisseur
                    }
                for (int i : g) {
                    v_[i] = vn;
                    u_[i] += dt_ * vn;
                }
            }
        }
        cundWork_ += cw;                   // V2/B4
        lysWork_ += lw;
        bcWork_ += bw;
        biasW_ += bias;
        if (scen_ != Scenario::TENSION && !toolNone_) tool_.integrate(dt_);
        return;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:cw,lw,bw,bias)
#endif
    for (int i = 0; i < (int)X0_.size(); ++i) {
        if (flag_[i] == FIXED) {
            // gripFree: frictionless grips — only the axial dof is held
            // (2D lesson: fully clamped grips decide where the specimen
            // breaks through the Saint-Venant corner concentration)
            if (gripFree_) {
                v_[i].x() += (dt_ / m_[i]) * f_[i].x();
                v_[i].y() += (dt_ / m_[i]) * f_[i].y();
                v_[i].z() = 0.0;
                u_[i] += dt_ * v_[i];
            } else v_[i].setZero();
            continue;
        }
        if (flag_[i] == PRESCRIBED) {
            if (gripFree_) {
                v_[i].x() += (dt_ / m_[i]) * f_[i].x();
                v_[i].y() += (dt_ / m_[i]) * f_[i].y();
            } else { v_[i].x() = 0.0; v_[i].y() = 0.0; }
            // pullRamp: cosine rise of the grip velocity — a stepped grip
            // launches a transient that unzips the first joint row under
            // the grip whatever the strength map says (2D lesson)
            double tv = t_ - pullDelay_;       // triaxial : sigma3 d'abord
            double vg = tv <= 0.0 ? 0.0 : pullV_;
            if (vg != 0.0 && pullRamp_ > 0.0 && tv < pullRamp_)
                vg *= 0.5 * (1.0 - std::cos(M_PI * tv / pullRamp_));
            {   // V2/B4 (corrige 2026-08-14, formule Sierra) : travail de la
                // LIAISON R = m.a - f en vitesse moyenne — les familles
                // comptent deja f.v sur ce noeud (l'ancien +f.vg avait le
                // signe inverse ET ignorait l'inertie de la rampe)
                double vgOld = v_[i].z();
                double Rz = m_[i] * (vg - vgOld) / dt_ - f_[i].z();
                bw += Rz * 0.5 * (vg + vgOld) * dt_;
            }
            v_[i].z() = vg;
            u_[i] += dt_ * v_[i];
            continue;
        }
        for (int a = 0; a < 3; ++a) {
            if (kAbs_[i](a) > 0) {
                double fk = kAbs_[i](a) * u_[i](a);
                f_[i](a) -= fk;
                lw -= fk * v_[i](a) * dt_;     // V2/B4 : ressort frontiere
            }
            if (damping_ > 0) {
                double fd = damping_ * std::abs(f_[i](a))
                            * (v_[i](a) > 0 ? 1.0 : (v_[i](a) < 0 ? -1.0 : 0.0));
                f_[i](a) -= fd;
                cw -= fd * v_[i](a) * dt_;     // V2/B4 : Cundall (<= 0)
            }
        }
        bias += f_[i].squaredNorm() * dt_ * dt_ / (2.0 * m_[i]);
        v_[i] += (dt_ / m_[i]) * f_[i];
        for (int a = 0; a < 3; ++a)
            if (cAbs_[i](a) > 0) {
                v_[i](a) /= 1.0 + dt_ * cAbs_[i](a) / m_[i];
                lw -= cAbs_[i](a) * v_[i](a) * v_[i](a) * dt_;  // amortisseur
            }
        u_[i] += dt_ * v_[i];
    }
    cundWork_ += cw;                       // V2/B4
    lysWork_ += lw;
    bcWork_ += bw;
    biasW_ += bias;
    if (scen_ != Scenario::TENSION && !toolNone_) tool_.integrate(dt_);
}

// ---------------------------------------------------------------------------
// Force volumique 3D : gravite + anti-gravite du tri des fragments.
// Absente (ou nulle) laisse tout modele existant bit-identique.
// ---------------------------------------------------------------------------
void Fdem3dSolver::bodyForces() {
    if (gravity_ > 0.0) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int eI = 0; eI < (int)el_.size(); ++eI) {
            const Elem& e = el_[eI];
            double w = phases_.mat[e.phase].rho * e.V0 * gravity_ / 4.0;
            for (int a = 0; a < 4; ++a) f_[e.n[a]].z() -= w;
        }
    }
    // anti-gravite sur les SEULS candidats. Serie : le travail s'accumule dans
    // un scalaire, et la phase de tri est courte devant le run.
    if (!brushArmed_) return;
    double bw = 0.0;
    for (int eI = 0; eI < (int)el_.size(); ++eI) {
        if (!brushCand_[eI]) continue;
        const Elem& e = el_[eI];
        double w = phases_.mat[e.phase].rho * e.V0 * brushA_ / 4.0;
        for (int a = 0; a < 4; ++a) {
            Eigen::Vector3d F = w * brushDir_;
            f_[e.n[a]] += F;
            bw += F.dot(v_[e.n[a]]) * dt_;
        }
    }
    brushWork_ += bw;                      // POSTE SEPARE, jamais dans sumW
}

// ---------------------------------------------------------------------------
// Armement du tri : fige les candidats et l'etat de reference.
// ---------------------------------------------------------------------------
void Fdem3dSolver::armBrush() {
    computeFragments();
    brushFrag_ = fragId_;
    brushNFrag_ = nFrag_;
    brushCand_.assign(el_.size(), 0);
    long nc = 0;
    for (std::size_t e = 0; e < el_.size(); ++e)
        if (fragId_[e] != 0) { brushCand_[e] = 1; ++nc; }
    brushU0_ = u_;
    std::vector<char> touched(X0_.size(), 0);
    for (std::size_t e = 0; e < el_.size(); ++e)
        if (brushCand_[e])
            for (int a = 0; a < 4; ++a) touched[el_[e].n[a]] = 1;
    double vres = 0.0, vresAll = 0.0;
    for (std::size_t i = 0; i < X0_.size(); ++i) {
        double vi = v_[i].norm();
        vresAll = std::max(vresAll, vi);
        if (touched[i]) vres = std::max(vres, vi);
    }
    long nn = 0;
    for (std::size_t i = 0; i < X0_.size(); ++i)
        if (touched[i]) {
            if (brushZeroV_) v_[i]  = brushV0_ * brushDir_;
            else             v_[i] += brushV0_ * brushDir_;
            ++nn;
        }
    brushArmed_ = true;
    brushT0_ = t_;
    std::cout << "[FDEM3D] TRI DES FRAGMENTS arme a t = " << t_ << " s : "
              << nc << " / " << el_.size() << " tetraedres candidats ("
              << brushNFrag_ - 1 << " fragments hors corps principal), " << nn
              << " noeuds\n"
              << "[FDEM3D]   v0 = " << brushV0_ << " m/s, a = " << brushA_
              << " m/s2, direction (" << brushDir_.x() << ", " << brushDir_.y()
              << ", " << brushDir_.z() << "), beta = " << brushBeta_ << "\n"
              << "[FDEM3D]   CLASSIFICATION, pas de physique : travail compte "
                 "a part, hors bilan B4\n"
              << "[FDEM3D]   REPOS (etape 1) : vitesse residuelle max des "
                 "candidats " << vres << " m/s, du bloc entier " << vresAll
              << " m/s\n";
    if (vres > 100.0 * brushV0_)
        std::cout << "[FDEM3D]   >>> ETAPE 1 VIOLEE : les candidats bougent "
                     "encore a " << vres << " m/s. Le classement mesurera le "
                     "mouvement RESIDUEL, pas la reponse a l'anti-gravite. "
                     "Allonger le repos entre toolStop et fragBrushStart, ou "
                     "poser fragBrushZeroV = true.\n";
}

// ---------------------------------------------------------------------------
// Verdict : deplacement de chaque fragment contre celui d'une particule LIBRE.
// ---------------------------------------------------------------------------
void Fdem3dSolver::brushReport() {
    if (!brushArmed_) return;
    double tb = t_ - brushT0_;
    double dref = brushV0_ * tb + 0.5 * brushA_ * tb * tb;
    double seuil = brushBeta_ * dref;
    std::vector<double> sum(brushNFrag_, 0.0), vol(brushNFrag_, 0.0);
    std::vector<long> cnt(brushNFrag_, 0);
    for (std::size_t e = 0; e < el_.size(); ++e) {
        int f = brushFrag_[e];
        double d = 0.0;
        for (int a = 0; a < 4; ++a)
            d += (u_[el_[e].n[a]] - brushU0_[el_[e].n[a]]).norm() / 4.0;
        sum[f] += d; ++cnt[f];
        vol[f] += el_[e].V0;
    }
    double vLibre = 0.0, vCoince = 0.0;
    long nLibre = 0, nCoince = 0;
    for (int f = 1; f < brushNFrag_; ++f) {
        if (cnt[f] == 0) continue;
        if (sum[f] / cnt[f] > seuil) { vLibre += vol[f]; ++nLibre; }
        else                         { vCoince += vol[f]; ++nCoince; }
    }
    double vTot = vLibre + vCoince;
    std::cout << "[FDEM3D] TRI, verdict apres " << tb << " s :\n"
              << "[FDEM3D]   particule libre de reference : d_ref = " << dref
              << " m, seuil = beta d_ref = " << seuil << " m\n"
              << "[FDEM3D]   fragments LIBRES  : " << nLibre << " (" << vLibre
              << " m^3)\n"
              << "[FDEM3D]   fragments COINCES : " << nCoince << " (" << vCoince
              << " m^3)\n"
              << "[FDEM3D]   volume detache CORRIGE " << vLibre << " contre "
              << vTot << " par le critere naif : "
              << (vTot > 0 ? 100.0 * (vTot - vLibre) / vTot : 0.0)
              << " % de SUREVALUATION\n"
              << "[FDEM3D]   travail du tri : " << brushWork_
              << " J (hors bilan B4)\n";
}

void Fdem3dSolver::computeFragments() {
    std::vector<std::vector<int>> adj(el_.size());
    for (const auto& J : jt_)
        if (J.D < 1.0) {
            adj[J.eA].push_back(J.eB);
            adj[J.eB].push_back(J.eA);
        }
    std::fill(fragId_.begin(), fragId_.end(), -1);
    int nid = 0;
    std::vector<std::pair<int, double>> mass;
    for (int s = 0; s < (int)el_.size(); ++s) {
        if (fragId_[s] >= 0) continue;
        double mm = 0;
        std::queue<int> qu;
        qu.push(s);
        fragId_[s] = nid;
        while (!qu.empty()) {
            int e = qu.front(); qu.pop();
            mm += rhoP_[el_[e].phase] * el_[e].V0;
            for (int w : adj[e])
                if (fragId_[w] < 0) { fragId_[w] = nid; qu.push(w); }
        }
        mass.push_back({nid, mm});
        ++nid;
    }
    std::sort(mass.begin(), mass.end(),
              [](auto& a, auto& b) { return a.second > b.second; });
    std::vector<int> rank(nid);
    for (int k = 0; k < nid; ++k) rank[mass[k].first] = k;
    for (auto& fid : fragId_) fid = rank[fid];
    nFrag_ = nid;
    // E7 (2026-08-19) : le volume detache ne compte plus les CORPS ETRANGERS.
    // Un maillage multi-corps (groupes physiques Gmsh) fait de l'insert une
    // composante connexe distincte de la roche : il etait donc compte comme
    // « fragment detache » a chaque run — 5 446 mm^3, soit tout son volume,
    // dans le resume et dans l'energie specifique work_/detachedVol_. On ne
    // retient desormais que les elements appartenant au MEME groupe physique
    // que le plus gros fragment, c'est-a-dire a la roche. En mono-corps
    // (elemGroup_ vide) le comportement est strictement inchange.
    int rockGroup = -1;
    if (!elemGroup_.empty())
        for (int e = 0; e < (int)el_.size(); ++e)
            if (fragId_[e] == 0) { rockGroup = elemGroup_[e]; break; }
    double vDet = 0;                                   // volume, phase-neutral
    long nForeign = 0;
    for (int e = 0; e < (int)el_.size(); ++e) {
        if (fragId_[e] == 0) continue;                 // le corps principal
        if (rockGroup >= 0 && elemGroup_[e] != rockGroup) { ++nForeign; continue; }
        vDet += el_[e].V0;
    }
    if (nForeign > 0)
        std::cout << "[FDEM3D] volume detache : " << nForeign
                  << " elements ecartes (corps etrangers a la roche : outil "
                     "maille, platines)\n";
    detachedVol_ = vDet;
}

// ===========================================================================

void Fdem3dSolver::writeFrame(int frame) {
    computeFragments();

    std::vector<Eigen::Vector3d> pts(X0_.size()), vel(X0_.size());
    for (std::size_t i = 0; i < X0_.size(); ++i) {
        pts[i] = X0_[i] + u_[i];
        vel[i] = v_[i];
    }
    std::vector<std::array<int, 4>> tets(el_.size());
    std::vector<double> svm(el_.size()), frag(el_.size());
    std::vector<double> phs(el_.size()), grn(el_.size());
    // Le 3D n'ecrivait QUE vonMises : impossible de distinguer ce qui tire de
    // ce qui cisaille alors que les deux modes de rupture le font. On sort donc
    // les deux moteurs, decomposition spectrale du Cauchy stocke sigG :
    //   sigma1  = contrainte principale MAJEURE, > 0 = TRACTION (moteur mode I)
    //   tauMax  = (sigma1 - sigma3)/2, cisaillement maximal  (moteur mode II)
    // Purement sortie : aucun effet sur le calcul (bit-neutre), le cout d'une
    // decomposition 3x3 n'etant paye qu'aux ecritures de frame.
    std::vector<double> sg1(el_.size()), tmx(el_.size());
    for (std::size_t e = 0; e < el_.size(); ++e) {
        tets[e] = el_[e].n;
        svm[e] = el_[e].svm;
        frag[e] = fragId_[e];
        phs[e] = el_[e].phase;
        grn[e] = el_[e].grain;
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(el_[e].sigG);
        const Eigen::Vector3d& ev = es.eigenvalues();   // triees croissant
        sg1[e] = ev(2);
        tmx[e] = 0.5 * (ev(2) - ev(0));
    }
    char name[64];
    std::snprintf(name, sizeof(name), "/fdem3d_%04d.vtu", frame);
    if (bdOn_) {
        std::vector<double> bdv(el_.size());
        for (std::size_t e = 0; e < el_.size(); ++e) bdv[e] = el_[e].bdD;
        vtk::writeTetMesh(out_ + name, pts, tets,
                          {{"vonMises", &svm}, {"sigma1", &sg1},
                           {"tauMax", &tmx}, {"fragment", &frag},
                           {"phase", &phs}, {"grain", &grn},
                           {"bulkD", &bdv}},
                          {{"velocity", &vel}});
    } else
    vtk::writeTetMesh(out_ + name, pts, tets,
                      {{"vonMises", &svm}, {"sigma1", &sg1},
                       {"tauMax", &tmx}, {"fragment", &frag},
                       {"phase", &phs}, {"grain", &grn}},
                      {{"velocity", &vel}});

    std::vector<std::array<int, 3>> tris;
    std::vector<double> Dj, tb, Tp, Fs, Bd, Fm, Bm, Dt, Ed;
    for (const auto& J : jt_) {
        tris.push_back(J.a);
        Dj.push_back(J.D);
        tb.push_back(J.tBreak);
        Tp.push_back(J.type);
        Fs.push_back(J.stat);
        Bd.push_back(J.bonded ? 1.0 : 0.0);
        Fm.push_back(J.failMode);
        Bm.push_back(J.bmode);
        // E8 (2026-08-19) : le DIF gele a l'insertion et le taux de
        // deformation qui l'a produit sont desormais EXPORTES. Sans eux la
        // population de joints inseree n'est pas auditable — or c'est
        // exactement la qu'un attracteur s'etait loge le 18/08 (l'exposant
        // litteral de Yang empilant les insertions juste sous 1e2 /s, mediane
        // 99,36). On ne peut pas surveiller ce qu'on n'ecrit pas.
        Dt.push_back(J.difT);
        Ed.push_back(J.edotIns);
    }
    std::snprintf(name, sizeof(name), "/fdem3d_joints_%04d.vtu", frame);
    vtk::ScalarField jf{
        {"damage", &Dj}, {"tBreak", &tb}, {"type", &Tp},
        {"ftScale", &Fs}, {"bonded", &Bd}, {"breakMode", &Bm}};
    if (difOn_) { jf["difT"] = &Dt; jf["edotIns"] = &Ed; }
    if (cfg_.getb("writeJointMode", false)) jf["failMode"] = &Fm;
    vtk::writeTriangles3(out_ + name, pts, tris, jf);

    std::ofstream fm(out_ + "/frames.csv",
                     frame == 0 ? std::ios::trunc : std::ios::app);
    if (frame == 0) fm << "frame,t,toolX,toolY,toolZ\n";
    fm << frame << "," << t_ << "," << tool_.x.x() << "," << tool_.x.y()
       << "," << tool_.x.z() << "\n";
}

void Fdem3dSolver::historyHeader(std::ostream& os) const {
    if (scen_ == Scenario::TENSION) {
        os << "t,gripFz,sigma,sigmaPeak,nBroken";
        for (int g : trkGrps_)                         // WP3 : par corps
            os << ",z_" << groupName_[g] << ",vz_" << groupName_[g];
        for (const auto& gg : gauges_)                 // WP3 : jauges
            os << ",szz_" << groupName_[gg.grp];
        os << "\n";
        return;
    }
    os << "t,toolFx,toolFy,toolFz,toolX,toolY,toolZ,toolVx,toolVy,toolVz,"
          "work,toolKE,nBroken,nFrag,detachedVol,specificEnergy";
    if (trackGroup_ >= 0)                              // corps suivi (V1+V2)
        os << ",grpZ,grpVz,grpFx,grpFy,grpFz,grpSzz";
    for (int g : trkGrps_)                             // WP3 : par corps
        os << ",z_" << groupName_[g] << ",vz_" << groupName_[g];
    for (const auto& gg : gauges_)                     // WP3 : jauges
        os << ",szz_" << groupName_[gg.grp];
    // V2/B4 : travaux cumules par famille (signes : negatif = preleve)
    os << ",eEl,eJnt,eGc,eFric,eCund,eLys";
    if (bdOn_) os << ",nPulv,bdWork";
    os << "\n";
}

void Fdem3dSolver::historyRow(std::ostream& os) const {
    if (scen_ == Scenario::TENSION) {
        os << t_ << "," << gripF_.z() << ","
           << std::abs(gripF_.z()) / (W_ * D_) << ","
           << sigmaPeak_ << "," << nBroken_;
        for (int g : trkGrps_) {                       // WP3 : par corps
            double mz = 0.0, mvz = 0.0, mm = 0.0;
            for (std::size_t e = 0; e < el_.size(); ++e) {
                if (elemGroup_[e] != g) continue;
                for (int a = 0; a < 4; ++a) {
                    int i = el_[e].n[a];
                    mm += m_[i];
                    mz += m_[i] * (X0_[i].z() + u_[i].z());
                    mvz += m_[i] * v_[i].z();
                }
            }
            os << "," << (mm > 0 ? mz / mm : 0.0)
               << "," << (mm > 0 ? mvz / mm : 0.0);
        }
        for (const auto& gg : gauges_) {               // WP3 : jauges
            double vs = 0.0, sz = 0.0;
            for (int e : gg.elems) {
                vs += el_[e].V0;
                sz += el_[e].V0 * el_[e].sigG(2, 2);
            }
            os << "," << (vs > 0 ? sz / vs : 0.0);
        }
        os << "\n";
        return;
    }
    double Es = detachedVol_ > 0 ? work_ / detachedVol_ : 0.0;
    os << t_ << "," << tool_.F.x() << "," << tool_.F.y() << "," << tool_.F.z()
       << "," << tool_.x.x() << "," << tool_.x.y() << "," << tool_.x.z()
       << "," << tool_.v.x() << "," << tool_.v.y() << "," << tool_.v.z()
       << "," << work_ << "," << tool_.ke() << "," << nBroken_ << ","
       << nFrag_ << "," << detachedVol_ << "," << Es;
    if (trackGroup_ >= 0) {
        // corps suivi : z du centroide massique et vz moyenne (ponderation
        // par les masses nodales — les noeuds appartiennent a un seul corps)
        double mz = 0.0, mvz = 0.0, mm = 0.0;
        for (std::size_t e = 0; e < el_.size(); ++e) {
            if (elemGroup_[e] != trackGroup_) continue;
            for (int a = 0; a < 4; ++a) {
                int i = el_[e].n[a];
                mm += m_[i];
                mz += m_[i] * (X0_[i].z() + u_[i].z());
                mvz += m_[i] * v_[i].z();
            }
        }
        os << "," << (mm > 0 ? mz / mm : 0.0)
           << "," << (mm > 0 ? mvz / mm : 0.0);
        // V2/B2 : force de contact nette du pas + jauge sigma_zz volumique
        double vsum = 0.0, szz = 0.0;
        for (std::size_t e = 0; e < el_.size(); ++e) {
            if (elemGroup_[e] != trackGroup_) continue;
            vsum += el_[e].V0;
            szz += el_[e].V0 * el_[e].sigG(2, 2);
        }
        os << "," << grpF_.x() << "," << grpF_.y() << "," << grpF_.z()
           << "," << (vsum > 0 ? szz / vsum : 0.0);
    }
    for (int g : trkGrps_) {                           // WP3 : par corps
        double mz = 0.0, mvz = 0.0, mm = 0.0;
        for (std::size_t e = 0; e < el_.size(); ++e) {
            if (elemGroup_[e] != g) continue;
            for (int a = 0; a < 4; ++a) {
                int i = el_[e].n[a];
                mm += m_[i];
                mz += m_[i] * (X0_[i].z() + u_[i].z());
                mvz += m_[i] * v_[i].z();
            }
        }
        os << "," << (mm > 0 ? mz / mm : 0.0)
           << "," << (mm > 0 ? mvz / mm : 0.0);
    }
    for (const auto& gg : gauges_) {                   // WP3 : jauges
        double vs = 0.0, sz = 0.0;
        for (int e : gg.elems) {
            vs += el_[e].V0;
            sz += el_[e].V0 * el_[e].sigG(2, 2);
        }
        os << "," << (vs > 0 ? sz / vs : 0.0);
    }
    os << "," << elWork_ << "," << jointWork_ << "," << gcWork_ << ","
       << gcFricWork_ << "," << cundWork_ << "," << lysWork_;   // V2/B4
    if (bdOn_) os << "," << nPulv_ << "," << bdWork_;
    os << "\n";
}

void Fdem3dSolver::finalize() {
    computeFragments();

    std::ofstream fe(out_ + "/fdem3d_final_elements.csv");
    fe << "cx,cy,cz,fragment,phase,grain\n";
    for (std::size_t e = 0; e < el_.size(); ++e) {
        Eigen::Vector3d c = Eigen::Vector3d::Zero();
        for (int a = 0; a < 4; ++a) c += X0_[el_[e].n[a]] + u_[el_[e].n[a]];
        c /= 4.0;
        fe << c.x() << "," << c.y() << "," << c.z() << "," << fragId_[e]
           << "," << el_[e].phase << "," << el_[e].grain << "\n";
    }

    double keBlock = 0.0;
    for (std::size_t i = 0; i < X0_.size(); ++i)
        keBlock += 0.5 * m_[i] * v_[i].squaredNorm();

    std::cout << "\n[FDEM3D] ---- summary ----\n"
              << "[FDEM3D] block kinetic energy at end: " << keBlock << " J\n"
              << "[FDEM3D] net work by general contact: " << gcWork_ << " J\n";
    if (contactPot_ && potStats_.pairs > 0) {
        const auto& S = potStats_;
        std::cout << "[FDEM3D] potential stats: " << S.pairs << " paires ("
                  << S.joint << " joint vivant, " << S.sepHint << " sep-hint, "
                  << S.sepFace << " sep-face, " << S.sepEdge << " sep-arete, "
                  << S.clipMiss << " clip-vide, " << S.clipHit
                  << " clip-force), tGrid = " << S.tGrid
                  << " s, tLoop = " << S.tLoop << " s\n";
    }
    // A dashpot can only DISSIPATE. A positive figure here means the viscous
    // branch is injecting energy — the rectifier failure mode — and every
    // number the run produced is suspect (2D lesson, now measured in 3D too).
    std::cout << "[FDEM3D] joint dashpot work: " << dampWork_ << " J  ["
              << (dampWork_ <= 0.0 ? "OK, dissipative"
                                   : "FAIL - the dashpot INJECTED energy")
              << "]\n";
    {   // ---- V2/B4 : bilan d'energie par sous-systeme -------------------
        // Theoreme travail-energie sur les noeuds : KE(t) - KE(0) = somme
        // des travaux par famille de forces + residu O(dt) (les compteurs
        // lisent v au moment de l'application, leapfrog). U_el se lit sur
        // le Cauchy stocke sigG (invariants par rotation, compliance
        // isotrope de la PHASE — exact en elastique, approx sous law/caps).
        double uEl = 0.0;
        for (const auto& e : el_) {
            double mu2 = mu2P_[e.phase], lam = lamP_[e.phase];
            double ss = e.sigG.squaredNorm(), trs = e.sigG.trace();
            uEl += e.V0 * (ss - lam * trs * trs / (3.0 * lam + mu2))
                   / (2.0 * mu2);
        }
        double uSpr = 0.0;                 // ressorts absorbants (stocke)
        for (std::size_t i = 0; i < X0_.size(); ++i)
            for (int a = 0; a < 3; ++a)
                if (kAbs_[i](a) > 0)
                    uSpr += 0.5 * kAbs_[i](a) * u_[i](a) * u_[i](a);
        double sumW = elWork_ + jointWork_ + gcWork_ + cundWork_ + lysWork_
                    + toolWork_ + bcWork_ + confWork_ + biasW_;
        double dKE = keBlock - keInit_;
        double resid = dKE - sumW;
        // echelle du verdict : le flux BRUT echange (la somme signee est ~0
        // precisement quand le bilan boucle — premier faux CHECK mesure sur
        // la percussion outil rigide : residu 0.014 J sur 1.53 J injectes,
        // soit 0.9 %, juge a 147 % de la somme signee)
        double gross = std::abs(elWork_) + std::abs(jointWork_)
                     + std::abs(gcWork_) + std::abs(cundWork_)
                     + std::abs(lysWork_) + std::abs(toolWork_)
                     + std::abs(bcWork_) + std::abs(confWork_);
        double scale = std::max({keInit_, keBlock, gross, 1e-30});
        // a charge nulle l'echelle est elle-meme un zero machine : le ratio
        // de deux zeros n'a pas de sens, le verdict se rend sur l'absolu
        bool zeroCase = scale < 1e-12;
        std::cout << "[FDEM3D] energy budget (V2/B4): KE " << keInit_
                  << " -> " << keBlock << " J\n"
                  << "[FDEM3D]   elements     : " << -elWork_
                  << " J preleve (stocke elastique " << uEl << " J)\n";
        if (viscOn_)
            std::cout << "[FDEM3D]      dont visqueux (2 mu D) : "
                      << -viscWork_ << " J dissipes, soit "
                      << (std::abs(elWork_) > 1e-300
                          ? 100.0 * viscWork_ / elWork_ : 0.0)
                      << " % du poste elements  "
                      << (viscWork_ <= 0.0
                          ? "[OK, dissipatif]"
                          : "[FAIL - le terme visqueux a INJECTE de l energie]")
                      << ". VENTILATION : deja comptee ci-dessus, pas un "
                         "poste de plus. Sous-estimee la ou det F < 1 (2 mu "
                         "D:D est par volume COURANT, V0 est de reference).\n";
        if (bdOn_)
            std::cout << "[FDEM3D]      dont pulverisation (bulkDamage) : "
                      << -bdWork_ << " J dissipes, " << nPulv_
                      << " elements a D = Dmax. VENTILATION : deja comptee "
                         "dans le poste elements.\n";
        std::cout << "[FDEM3D]   joints       : " << -(jointWork_ - dampWork_)
                  << " J cohesif (fissuration + stocke), dashpot "
                  << -dampWork_ << " J\n"
                  << "[FDEM3D]   contact      : " << -gcWork_
                  << " J (dont frottement " << -gcFricWork_ << " J)\n"
                  << "[FDEM3D]   Cundall      : " << -cundWork_ << " J\n"
                  << "[FDEM3D]   frontieres   : " << -lysWork_
                  << " J (dont stocke ressorts " << uSpr << " J)\n"
                  << "[FDEM3D]   outil->solide: " << toolWork_
                  << " J, platines: " << bcWork_ << " J\n";
        if (confP_ > 0.0)                  // sortie inchangee si pas confine
            std::cout << "[FDEM3D]   confinement  : " << confWork_
                      << " J (pression suiveuse -> solide)\n";
        std::cout << "[FDEM3D]   integration  : +" << biasW_
                  << " J (correction leapfrog f^2 dt^2/2m)\n"
                  << "[FDEM3D]   residu       : " << resid << " J ("
                  << 100.0 * std::abs(resid) / scale << " % de l'echelle) ["
                  << (zeroCase ? "OK (zero machine)"
                      : std::abs(resid) <= 0.01 * scale ? "OK" : "CHECK")
                  << "]\n";
    }
    {   // failure-mode census (breakMode: 1 = tensile, 2 = shear)
        long nTm = 0, nSm = 0;
        for (const auto& J : jt_) {
            if (J.bmode == 1) ++nTm;
            else if (J.bmode == 2) ++nSm;
        }
        if (nTm + nSm > 0)
            std::cout << "[FDEM3D] breakage mode: " << nTm << " tensile, "
                      << nSm << " shear ("
                      << 100.0 * nSm / (double)(nTm + nSm) << " % shear)\n";
    }
    if (adaptive_)
        std::cout << "[FDEM3D] adaptive insertion: " << nInserted_ << " / "
                  << jt_.size() << " joints inserted ("
                  << (jt_.empty() ? 0.0 : 100.0 * nInserted_ / jt_.size())
                  << " %), " << nBroken_ << " fully broken\n";
    if (difOn_) {
        std::vector<double> er, dtv;
        for (const auto& J : jt_)
            if (!J.bonded) { er.push_back(J.edotIns); dtv.push_back(J.difT); }
        if (er.empty()) {
            std::cout << "[FDEM3D] strainRateDIF : aucun joint insere — le "
                         "DIF n a jamais ete evalue.\n";
        } else {
            std::sort(er.begin(), er.end());
            std::sort(dtv.begin(), dtv.end());
            std::size_t n = er.size();
            std::cout << "[FDEM3D] strainRateDIF (sur " << n
                      << " joints inseres) : taux de deformation a "
                         "l insertion, mediane " << er[n / 2] << " /s, min "
                      << er.front() << ", max " << er.back()
                      << " ; DIF_traction median " << dtv[n / 2]
                      << " (min " << dtv.front() << ", max " << dtv.back()
                      << ")\n";
            std::size_t sat = 0;
            for (double x : er) if (x > 1.0e2) ++sat;
            std::cout << "[FDEM3D]   " << (100.0 * sat / n)
                      << " % des insertions sont AU PLATEAU de DIF_traction "
                         "(edot > 1e2 /s) : sur cette part le facteur ne "
                         "discrimine plus rien.\n";
        }
    }
    if (nGroups_ > 1) {                    // V1 : bilan par corps
        for (int g = 0; g < nGroups_; ++g) {
            double ke = 0.0, mvz = 0.0, mm = 0.0;
            for (std::size_t e = 0; e < el_.size(); ++e) {
                if (elemGroup_[e] != g) continue;
                for (int a = 0; a < 4; ++a) {
                    int i = el_[e].n[a];
                    ke += 0.5 * m_[i] * v_[i].squaredNorm();
                    mm += m_[i];
                    mvz += m_[i] * v_[i].z();
                }
            }
            std::cout << "[FDEM3D] corps '" << groupName_[g] << "': KE = "
                      << ke << " J, vz moyenne = "
                      << (mm > 0 ? mvz / mm : 0.0) << " m/s, masse = "
                      << mm << " kg\n";
        }
    }
    if (gcAdaptive_)
        std::cout << "[FDEM3D] adaptive contact activation: " << nActivated_
                  << " / " << pool_.size() << " exterior faces activated ("
                  << (pool_.empty() ? 0.0
                                    : 100.0 * nActivated_ / pool_.size())
                  << " %)\n";
    if (confP_ > 0.0)
        std::cout << "[FDEM3D] confinement: vise " << -confP_ / 1e6
                  << " MPa lateral, atteint (coeur, apres equilibrage) "
                  << confAchieved_ / 1e6 << " MPa ("
                  << 100.0 * std::abs(confAchieved_ / -confP_) << " %)\n";

    if (voronoi_) {
        // grain/phase bookkeeping: achieved volume fractions and the
        // inter/intra-granular split of the broken joints — the observable
        // a GBM exists to produce.
        std::vector<double> vPh(phases_.n(), 0.0);
        double vTot = 0.0;
        for (const auto& e : el_) { vPh[e.phase] += e.V0; vTot += e.V0; }
        std::cout << "[FDEM3D] grains: " << nGrains_ << ", phases:";
        for (int p = 0; p < phases_.n(); ++p)
            std::cout << " " << phases_.name[p] << " "
                      << 100.0 * vPh[p] / vTot << "%"
                      << " (target " << 100.0 * phases_.fraction[p] << "%)";
        std::cout << "\n";
        long nt[3] = {0, 0, 0}, nb[3] = {0, 0, 0};
        for (const auto& J : jt_) {
            ++nt[J.type];
            if (J.D >= 1.0) ++nb[J.type];
        }
        std::cout << "[FDEM3D] joints intra/homo/hetero: " << nt[0] << "/"
                  << nt[1] << "/" << nt[2] << ", broken: " << nb[0] << "/"
                  << nb[1] << "/" << nb[2];
        long nbTot = nb[0] + nb[1] + nb[2];
        if (nbTot > 0)
            std::cout << "  (intergranular fraction "
                      << 100.0 * (nb[1] + nb[2]) / (double)nbTot << " %)";
        std::cout << "\n";
    }

    if (scen_ == Scenario::TENSION) {
        if (!cfg_.getb("verifyFt", true)) {
            std::cout << "[FDEM3D] tension result (GB-controlled, verifyFt=off):\n"
                      << "[FDEM3D]   peak macro stress = " << sigmaPeak_ / 1e6
                      << " MPa = " << sigmaPeak_ / mat_.ft << " x bulk ft"
                      << " (weak-boundary reference gbAlphaTen = "
                      << phases_.aTen << ")\n"
                      << "[FDEM3D]   broken joints = " << nBroken_ << " / "
                      << jt_.size() << "\n";
            return;
        }
        double err = 100.0 * (sigmaPeak_ - mat_.ft) / mat_.ft;
        // On the voronoi mesh the crack must follow a tortuous joint path
        // whose facets are inclined to the loading axis (each sees
        // sigma cos^2 theta), so the macroscopic peak sits AT or somewhat
        // ABOVE ft: the band is widened on that side, as in 2D. The Kuhn
        // grid keeps complete horizontal joint planes and must hit ft.
        double hi = voronoi_ ? 25.0 : 5.0;
        bool pass = err > -5.0 && err < hi;
        std::cout << "[FDEM3D] tension verification ("
                  << (voronoi_ ? "voronoi grains"
                               : "full horizontal joint planes exist") << "):\n"
                  << "[FDEM3D]   peak macro stress = " << sigmaPeak_ / 1e6
                  << " MPa, expected ft = " << mat_.ft / 1e6
                  << " MPa, error = " << err << " %  ["
                  << (pass ? "PASS" : "FAIL") << "]\n"
                  << "[FDEM3D]   broken joints = " << nBroken_ << " / "
                  << jt_.size() << "\n";
        return;
    }
    double Es = detachedVol_ > 0 ? work_ / detachedVol_ : 0.0;
    std::cout << "[FDEM3D] peak tool force   : " << peakF_ << " N\n"
              << "[FDEM3D] tool work output  : " << work_ << " J";
    if (tool_.free)
        std::cout << "  (tool KE loss: " << toolKE0_ - tool_.ke() << " J)";
    std::cout << "\n[FDEM3D] broken joints     : " << nBroken_ << " / " << jt_.size()
              << "\n[FDEM3D] fragments         : " << nFrag_
              << " (detached vol " << detachedVol_ << " m^3)"
              << "\n[FDEM3D] specific energy   : " << Es << " J/m^3\n";
    brushReport();          // no-op si le tri n'a pas ete arme
}

// ---------------------------------------------------------------------------
// selftest-potential3d — le test decisif du portage 3D (A3 phase 2).
// Deux TETS rigides (centre + rotation par Rodrigues, equations d'Euler en
// repere monde), lies uniquement par le contact par potentiel tet-tet, sans
// frottement. Frontale symetrique (transfert exact attendu) puis oblique
// (couple non nul, le polyedre de coupe tourne et se deforme). Compteur de
// travail de convention solveur (f.v avant le kick) — la conservation se
// juge sur dKE, comme en 2D.
// ---------------------------------------------------------------------------
int potentialSelftest3d(const std::string& csvPath) {
    using V3 = Eigen::Vector3d;
    using M3 = Eigen::Matrix3d;
    std::ofstream csv(csvPath);
    csv << "phase,t,xA,vAx,xB,vBx,vol,W,KE\n";

    struct Rigid {
        V3 c, v;
        M3 R;
        V3 om;
        V3 r0[4];
        double m = 4.0;                    // 4 masses ponctuelles unite
        M3 Ib;                             // inertie en repere corps
        void setTet(const V3& A, const V3& B, const V3& C, const V3& D) {
            c = (A + B + C + D) / 4.0;
            r0[0] = A - c;
            r0[1] = B - c;
            r0[2] = C - c;
            r0[3] = D - c;
            R.setIdentity();
            om.setZero();
            v.setZero();
            Ib.setZero();
            for (int k = 0; k < 4; ++k)
                Ib += r0[k].squaredNorm() * M3::Identity()
                      - r0[k] * r0[k].transpose();
        }
        void pos(V3 p[4]) const {
            for (int k = 0; k < 4; ++k) p[k] = c + R * r0[k];
        }
        V3 vel(const V3& x) const { return v + om.cross(x - c); }
        M3 Iw() const { return R * Ib * R.transpose(); }
    };
    auto rodrigues = [](const V3& w, double dt) {
        double th = w.norm() * dt;
        if (th < 1e-30) return M3(M3::Identity());
        V3 k = w.normalized();
        M3 K;
        K << 0, -k.z(), k.y(), k.z(), 0, -k.x(), -k.y(), k.x(), 0;
        return M3(M3::Identity() + std::sin(th) * K
                  + (1.0 - std::cos(th)) * K * K);
    };

    const double p = 1.0e3, dt = 1.0e-4, v0 = 1.0;
    int fails = 0;
    double worstW = 0.0, worstKE = 0.0, worstP = 0.0;

    for (int phase = 1; phase <= 2; ++phase) {
        Rigid A, B;
        // POINTE-contre-FACE, la configuration generique du contact FDEM :
        // la BASE de A (a x = -0.10, tournee vers +x) recoit l'APEX de B
        // (a x = -0.08, pointant vers -x), les deux centroides sur l'axe de
        // vol. Le tip-contre-tip de deux tets a ete essaye et ecarte : la
        // repulsion y est quasi nulle (phi -> 0 aux sommets), les deux cones
        // broutent ~10 000 pas en position degeneree du clip et l'erreur
        // d'integration s'accumule (mesure : dKE 5e-3, contre 1e-8 ici).
        A.setTet(V3(-0.10, -0.45, -0.26), V3(-0.10, 0.00, 0.52),
                 V3(-0.10, 0.45, -0.26), V3(-1.15, 0.00, 0.00));
        double yoff = (phase == 2) ? 0.18 : 0.0;
        B.setTet(V3(1.00, -0.45 + yoff, -0.26), V3(1.00, 0.00 + yoff, 0.52),
                 V3(1.00, 0.45 + yoff, -0.26), V3(-0.08, yoff, 0.00));
        A.v = V3(v0, 0, 0);

        double KE0 = 0.5 * A.m * A.v.squaredNorm();
        V3 P0 = A.m * A.v + B.m * B.v;
        double W = 0.0, volMax = 0.0;
        long nTouch = 0;

        const long nSteps = (long)(3.0 / dt);
        for (long s = 0; s < nSteps; ++s) {
            V3 pa[4], pb[4];
            A.pos(pa);
            B.pos(pb);
            pot3::PairForce3 Rp;
            V3 FA = V3::Zero(), FB = V3::Zero();
            V3 tA = V3::Zero(), tB = V3::Zero();
            if (pot3::pairForce(pa, pb, p, Rp)) {
                ++nTouch;
                volMax = std::max(volMax, Rp.vol);
                for (int k = 0; k < 4; ++k) {
                    W += (Rp.fA[k].dot(A.vel(pa[k]))
                          + Rp.fB[k].dot(B.vel(pb[k]))) * dt;
                    FA += Rp.fA[k];
                    FB += Rp.fB[k];
                    tA += (pa[k] - A.c).cross(Rp.fA[k]);
                    tB += (pb[k] - B.c).cross(Rp.fB[k]);
                }
            }
            auto kick = [&](Rigid& X, const V3& F, const V3& tq) {
                X.v += F / X.m * dt;
                M3 I = X.Iw();
                X.om += I.inverse() * (tq - X.om.cross(I * X.om)) * dt;
                X.c += X.v * dt;
                X.R = rodrigues(X.om, dt) * X.R;
            };
            kick(A, FA, tA);
            kick(B, FB, tB);
            if (s % 200 == 0) {
                double KE = 0.5 * A.m * A.v.squaredNorm()
                          + 0.5 * B.m * B.v.squaredNorm()
                          + 0.5 * A.om.dot(A.Iw() * A.om)
                          + 0.5 * B.om.dot(B.Iw() * B.om);
                csv << phase << "," << s * dt << "," << A.c.x() << ","
                    << A.v.x() << "," << B.c.x() << "," << B.v.x() << ","
                    << volMax << "," << W << "," << KE << "\n";
            }
        }

        double KE1 = 0.5 * A.m * A.v.squaredNorm()
                   + 0.5 * B.m * B.v.squaredNorm()
                   + 0.5 * A.om.dot(A.Iw() * A.om)
                   + 0.5 * B.om.dot(B.Iw() * B.om);
        V3 P1 = A.m * A.v + B.m * B.v;
        double wRel = std::abs(W) / KE0;
        double keRel = std::abs(KE1 - KE0) / KE0;
        double pRel = (P1 - P0).norm() / P0.norm();
        worstW = std::max(worstW, wRel);
        worstKE = std::max(worstKE, keRel);
        worstP = std::max(worstP, pRel);
        std::cout << "[POT3] phase " << phase
                  << (phase == 1 ? " (frontale)" : " (oblique) ")
                  << ": contact " << nTouch << " pas, volume max " << volMax
                  << "\n[POT3]   vA_fin = (" << A.v.x() << ", " << A.v.y()
                  << ", " << A.v.z() << "), vB_fin = (" << B.v.x() << ", "
                  << B.v.y() << ", " << B.v.z() << "), |omB| = "
                  << B.om.norm()
                  << "\n[POT3]   |W_contact|/KE0 = " << wRel
                  << ", |dKE|/KE0 = " << keRel
                  << ", |dP|/|P0| = " << pRel << "\n";
        if (phase == 1) {
            if (std::abs(B.v.x() - v0) > 0.02 || std::abs(A.v.x()) > 0.02)
                ++fails;
        }
    }
    std::cout << "pot3_work_rel = " << worstW << "\n"
              << "pot3_ke_rel = " << worstKE << "\n"
              << "pot3_mom_rel = " << worstP << "\n";
    bool ok = fails == 0 && worstW < 5e-3 && worstKE < 1e-4
              && worstP < 1e-10;
    std::cout << (ok ? "[PASS]" : "[FAIL]")
              << " selftest-potential3d : contact conservatif de Munjiza en "
                 "3D (conservation jugee sur dKE ; biais O(dt) du compteur "
                 "documente en 2D)\n";
    return ok ? 0 : 1;
}

} // namespace rockim
