// ---------------------------------------------------------------------------
// FdemSolver — 2D combined finite-discrete element method (Munjiza-style).
// See the header for the model overview; comments here carry the derivations.
// ---------------------------------------------------------------------------
#include "rockim/FdemSolver.hpp"

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
#include "rockim/YangDif.hpp"
#include "rockim/RandomField.hpp"
#include "rockim/Tessellation.hpp"
#include "rockim/VtkWriter.hpp"

#include <chrono>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace rockim {

namespace {
struct FProf {
    double tEl = 0, tJt = 0, tGc = 0, tTc = 0, tIn = 0;
    long n = 0;
    bool on = std::getenv("ROCKIM_PROF") != nullptr;
    ~FProf() {
        if (!on || n == 0) return;
        std::fprintf(stderr,
                     "[prof] per step (ms): elem %.2f insert %.2f joint %.2f "
                     "gcontact %.2f tool %.2f  (%ld steps)\n",
                     1e3 * tEl / n, 1e3 * tIn / n, 1e3 * tJt / n,
                     1e3 * tGc / n, 1e3 * tTc / n, n);
    }
} fProf;
double fnow() {
    return std::chrono::duration<double>(
               std::chrono::steady_clock::now().time_since_epoch()).count();
}

// ---------------------------------------------------------------------------
// Facteurs d amplification dynamique de Yang, Xiang, Naderi, Wang, Aising,
// Ugarte & Latham (IJRMMS 191, 2025, 106125), leurs eq. 2 et 3, transcrites
// LITTERALEMENT. edot est un taux de deformation en 1/s.
//
// ATTENTION : la forme publiee de DIF_traction est DISCONTINUE a ses deux
// bornes (elle vaut 1,124 juste au-dessus de 5e-6 contre 1 en dessous, et
// 1,516 juste en dessous de 1e2 contre 1,85 au-dessus). On la transcrit telle
// quelle — c est le modele des auteurs — et le solveur imprime ces sauts au
// demarrage pour qu ils ne soient jamais decouverts en aval.
// Les deux facteurs DIF de Yang et al. et la mesure du taux principal 3x3
// vivent desormais dans include/rockim/YangDif.hpp : le solveur 3D les
// utilise aussi, et il ne doit exister qu UNE transcription de l article
// dans le depot. Les expressions y sont inchangees ; le commentaire sur la
// coquille de l exposant 0,07 et sur l attracteur en 1e2 /s y a suivi.
} // namespace

using rockim::difTensionYang;
using rockim::difCompressionYang;

FdemSolver::FdemSolver(const Config& cfg, std::string outDir)
    : cfg_(cfg), out_(std::move(outDir)) {}

void FdemSolver::init() {
    mat_ = Material::from(cfg_);
    phases_ = PhaseSet::from(cfg_);

    std::string sc = cfg_.gets("scenario", "percussion");
    if      (sc == "percussion") scen_ = Scenario::PERCUSSION;
    else if (sc == "shear")      scen_ = Scenario::SHEAR;
    else if (sc == "tension")    scen_ = Scenario::TENSION;
    else if (sc == "brazilian")  scen_ = Scenario::BRAZILIAN;
    else if (sc == "shpb")       scen_ = Scenario::SHPB;
    else throw std::runtime_error("fdem scenario must be percussion | shear | "
                                  "tension | brazilian | shpb");

    W_ = cfg_.getd("W", 0.2);
    H_ = cfg_.getd("H", 0.2);
    nx_ = cfg_.geti("nx", 64);
    ny_ = cfg_.geti("ny", 64);
    thk_ = cfg_.getd("thickness", 1.0);
    T_ = cfg_.getd("T", 2.5e-4);
    bool quasiStatic = scen_ == Scenario::TENSION || scen_ == Scenario::BRAZILIAN;
    // SHPB default 0: Cundall local damping is a non-physical force
    // proportional to |f| that ATTENUATES a travelling elastic wave, which
    // is exactly the quantity this test measures. Any non-zero value here
    // would show up as the "cohesive-element attenuation" of fig. 23.
    damping_ = cfg_.getd("dampingLocal",
                         scen_ == Scenario::SHPB ? 0.0
                                                 : (quasiStatic ? 0.7 : 0.02));
    // Viscosite de volume 2*mu*D de l'eq. 6 de Yan et al. (le terme de taux
    // DANS la contrainte — leur amortissement est VISQUEUX, Table 1 :
    // mu = 7,6e3 Pa.s, la ou rockim n'avait que le Cundall nodal). Opt-in :
    // 0 (defaut) = strictement bit-identique (le terme n'est pas calcule).
    bulkVisc_ = cfg_.getd("bulkViscosity", 0.0);
    if (bulkVisc_ < 0.0)
        throw std::runtime_error("bulkViscosity must be >= 0 [Pa.s] (eq. 6 "
                                 "de Yan et al. : terme 2*mu*D dans la "
                                 "contrainte)");
    // strainRateDIF = off (defaut) | yang — voir FdemSolver.hpp.
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
    // ---- WP1 : pulverisation (Yang et al. 2026) — portage 2D ------------
    // Miroir exact du 3D (principe III) : memes cles, meme loi. Voir
    // Fdem3dSolver pour la doc detaillee. delta_m = h_e * eps_vm avec
    // eps_zz = 0 (deformation plane) dans le deviateur.
    {
        std::string bd = cfg_.gets("bulkDamage", "off");
        if (bd != "off" && bd != "yang")
            throw std::runtime_error("bulkDamage must be off | yang "
                                     "(Yang et al. 2026, eq. 3-4)");
        bdOn_ = bd == "yang";
    }
    if (bdOn_) {
        bdD0_   = cfg_.getd("bulkDamageDelta0", 1.4e-5);
        bdDf_   = cfg_.getd("bulkDamageDeltaF", 4.0e-4);
        bdDmax_ = cfg_.getd("bulkDamageDmax", 0.9);
        bdCd_   = cfg_.getd("bulkDamageCd", 1.0);
        if (!(bdD0_ > 0.0) || !(bdDf_ > bdD0_))
            throw std::runtime_error("bulkDamage : il faut 0 < "
                                     "bulkDamageDelta0 < bulkDamageDeltaF [m]");
        if (!(bdDmax_ > 0.0) || bdDmax_ > 1.0 || !(bdCd_ > 0.0))
            throw std::runtime_error("bulkDamage : bulkDamageDmax dans "
                                     "]0, 1] et bulkDamageCd > 0");
        if (cfg_.gets("law", "elastic") != "elastic")
            throw std::runtime_error("bulkDamage = yang exige law = elastic");
    }
    mtCap_ = cfg_.getd("meanTensionCapFactor", 3.0);
    srTau_ = cfg_.getd("strainRateTau", 1.0e-6);
    if (difOn_ && !(srTau_ > 0.0))
        throw std::runtime_error("strainRateTau must be > 0 [s] (filtre du "
                                 "taux de deformation : le taux brut par "
                                 "element est trop bruite pour figer un DIF "
                                 "dessus)");
    // Body force. Absent (or 0) leaves every existing model bit-identical:
    // bodyForces() returns immediately and no other code path is touched.
    gravity_ = cfg_.getd("gravity", 0.0);
    if (gravity_ < 0.0)
        throw std::runtime_error("gravity is a magnitude in m/s^2 (it acts along -y): use a positive value");

    // ---- balai numerique (Yang et al. 2025) : voir FdemSolver.hpp ---------
    // etape 1 : arret de l'outil, DISTINCT de l'armement du balai
    toolStop_ = cfg_.getd("toolStop", 0.0);
    brushStart_ = cfg_.getd("fragBrushStart", 0.0);
    if (toolStop_ > 0.0 && brushStart_ > 0.0 && brushStart_ <= toolStop_)
        throw std::runtime_error("fragBrushStart doit etre APRES toolStop : "
                                 "c'est l'intervalle de repos qui rend le "
                                 "balai interpretable (Yang et al. 2025, "
                                 "etape 1 — « after completing the impact »)");
    if (brushStart_ > 0.0) {
        brushV0_   = cfg_.getd("fragBrushV0", 2.5e-3);
        brushA_    = cfg_.getd("fragBrushAccel", 98.1);
        brushBeta_ = cfg_.getd("fragBrushBeta", 0.8);
        brushZeroV_ = cfg_.gets("fragBrushZeroV", "false") == "true";
        brushDir_  = Eigen::Vector2d(cfg_.getd("fragBrushDirX", 0.0),
                                     cfg_.getd("fragBrushDirY", 1.0));
        if (brushDir_.norm() < 1e-12)
            throw std::runtime_error("fragBrushDir is null: give a direction");
        brushDir_.normalize();
        if (brushBeta_ <= 0.0)
            throw std::runtime_error("fragBrushBeta must be > 0 (0,8 chez Yan"
                                     "g et al. 2025, justifie par un plateau)");
    }

    // ---- contrainte in situ (etude tunnel EDZ) ---------------------------
    // insituSh / insituSv : contraintes principales horizontale et verticale
    // du massif, en PRESSION POSITIVE (5e6 = 5 MPa de compression), meme
    // convention que confiningPressure. Le coefficient de pression laterale
    // de l'article est lambda = insituSh / insituSv. Absentes (ou nulles) :
    // aucun chemin de code n'est touche, trajectoires bit-identiques.
    {
        double sh = cfg_.getd("insituSh", 0.0);
        double sv = cfg_.getd("insituSv", 0.0);
        double sxy0 = cfg_.getd("insituSxy", 0.0);
        if (sh < 0.0 || sv < 0.0)
            throw std::runtime_error("insituSh / insituSv sont des PRESSIONS "
                                     "positives (compression) : utiliser 5e6 "
                                     "pour 5 MPa de compression");
        if (sh > 0.0 || sv > 0.0 || sxy0 != 0.0) {
            hasInsitu_ = true;
            insituS_ << -sh, sxy0, sxy0, -sv;
            std::cout << "[FDEM] in situ : sigma_h = " << sh / 1e6
                      << " MPa, sigma_v = " << sv / 1e6
                      << " MPa (compression), lambda = "
                      << (sv > 0.0 ? sh / sv : 0.0) << "\n";
            // NB : adaptive_ n'est lu que plus bas dans init(), on relit donc
            // la cle elle-meme.
            if (cfg_.gets("insertion", "intrinsic") != "adaptive")
                std::cout << "[FDEM] WARNING: contrainte in situ avec "
                             "insertion = intrinsic — les joints intrinseques "
                             "ne transmettent sigma0 qu'apres s'etre fermes de "
                             "sigma0/pj, donc l'etat initial n'est PAS un "
                             "equilibre exact (transitoire parasite, controle "
                             "a charge nulle inexact). insertion = adaptive "
                             "rend l'etat initial exact : les aretes non "
                             "encore inserees sont des noeuds partages.\n";
            if (cfg_.getd("crushCap", 8.0 * cfg_.getd("cohesion", 0.0)) < 1e11)
                std::cout << "[FDEM] NOTE: le cap deviatorique (crushCap) ne "
                             "voit PAS la contrainte in situ (e.svm est "
                             "calcule avant l'ajout). Sous champ anisotrope il "
                             "peut mordre pour une raison purement numerique : "
                             "poser crushCap tres grand si le bulk doit rester "
                             "elastique.\n";
        }
    }

    // ---- geometry: box (default) or disc (mandatory for the brazilian) -----
    std::string geo = cfg_.gets("geometry",
                                scen_ == Scenario::BRAZILIAN ? "disc"
                              : scen_ == Scenario::SHPB      ? "shpb" : "box");
    if (geo != "box" && geo != "disc" && geo != "shpb")
        throw std::runtime_error("geometry must be box | disc | shpb (got '"
                                 + geo + "')");
    disc_ = geo == "disc";
    shpb_ = geo == "shpb";
    if ((scen_ == Scenario::SHPB) != shpb_)
        throw std::runtime_error("scenario = shpb and geometry = shpb go "
                                 "together (the SHPB assembly IS the geometry)");
    if (shpb_) {
        // ---- SHPB assembly (fig. 22 of Yan, Zheng & Wang 2023) -------------
        shpbLIB_   = cfg_.getd("shpbIncidentLength", 2.00);
        shpbLTB_   = cfg_.getd("shpbTransmitLength", 1.50);
        shpbD_     = cfg_.getd("shpbBarDiameter", 0.05);
        shpbDrock_ = cfg_.getd("shpbDiscDiameter", 0.050);
        shpbGap_   = cfg_.getd("shpbGap", 0.0);
        hBar_      = cfg_.getd("shpbBarElemSize", 5e-3);
        hDisc_     = cfg_.getd("shpbDiscElemSize", 7.5e-4);
        shpbNoDisc_ = cfg_.getb("shpbNoDisc", false);
        // NOTE: the SHPB disc is left ROUND, as in fig. 22 (flat bar ends
        // against a circular disc). discFlattenDeg is a brazilian-platen
        // key and is deliberately NOT honoured here.
        for (auto [v, nm] : {std::pair<double, const char*>{shpbLIB_, "shpbIncidentLength"},
                             {shpbLTB_, "shpbTransmitLength"},
                             {shpbD_, "shpbBarDiameter"},
                             {shpbDrock_, "shpbDiscDiameter"},
                             {hBar_, "shpbBarElemSize"},
                             {hDisc_, "shpbDiscElemSize"}})
            if (!(v > 0.0))
                throw std::runtime_error(std::string("shpb: ") + nm
                                         + " must be > 0");
        if (shpbGap_ < 0.0)
            throw std::runtime_error("shpbGap must be >= 0 (0 = faces flush)");
        std::string pl = cfg_.gets("shpbPulse", "halfsine");
        if (pl != "halfsine" && pl != "trapezoid")
            throw std::runtime_error("shpbPulse must be halfsine | trapezoid");
        shpbPulse_ = pl == "trapezoid" ? 1 : 0;
        shpbV0_  = cfg_.getd("shpbPulseV0", 5.2);
        shpbTau_ = cfg_.getd("shpbPulseTau", 2.2e-4);
        shpbPlateau_ = cfg_.getd("shpbPulsePlateau", 0.5);
        if (!(shpbV0_ > 0.0) || !(shpbTau_ > 0.0))
            throw std::runtime_error("shpbPulseV0 / shpbPulseTau must be > 0");
        if (!(shpbPlateau_ > 0.0 && shpbPlateau_ < 1.0))
            throw std::runtime_error("shpbPulsePlateau (plateau fraction of "
                                     "the trapezoid) must be in (0, 1)");
        absFac_ = cfg_.getd("absorbFactor", 1.0);
        if (!(absFac_ > 0.0))
            throw std::runtime_error("absorbFactor must be > 0 (1 = classical "
                                     "Lysmer rho c v, 2 = eq. 21 as printed)");
        double xDiscL = shpbLIB_ + shpbGap_;
        discR_ = 0.5 * shpbDrock_;
        discC_ = {xDiscL + discR_, 0.5 * shpbD_};
        W_ = shpbLIB_ + shpbLTB_ + shpbDrock_ + 2.0 * shpbGap_;
        H_ = shpbD_;
        xEndAbs_ = W_;
        shpbM1_ = cfg_.getd("shpbMonitor1", shpbLIB_ - 1.00);
        shpbM2_ = cfg_.getd("shpbMonitor2",
                            xDiscL + shpbDrock_ + shpbGap_ + 1.00);
        shpbGaugeW_ = cfg_.getd("shpbGaugeHalfLength", 2.0 * hBar_);
        // contact cell sized on the DISC elements (the small ones): a cell
        // of 8 bar elements put every interface face in one bucket and made
        // the detection O(n_faces^2) — measured 64 ms/step of the 80 ms.
        gcCell_ = cfg_.getd("gcCell", 2.0 * hDisc_);
        gcBoxMesh_ = cfg_.getb("gcBoxMesh", true);   // defaut PROPRE au SHPB
        double win = cfg_.getd("gcXwindow", 0.10);
        if (win > 0.0) {
            gcXmin_ = discC_.x() - discR_ - win;
            gcXmax_ = discC_.x() + discR_ + win;
        }
        brazStop_ = false;
    }
    if (scen_ == Scenario::BRAZILIAN && !disc_)
        throw std::runtime_error("scenario = brazilian needs geometry = disc "
                                 "(the indirect tension test is a DISC "
                                 "compressed diametrically)");
    if (disc_) {
        discR_ = 0.5 * std::min(W_, H_);
        discC_ = {0.5 * W_, 0.5 * H_};
        // FLATTENED brazilian disc: the standard answer to the failure mode
        // this solver measured on the plain disc — with a round rim the first
        // joints break AT THE CONTACT (measured: first break at r/R = 0.98,
        // 9 deg from the loaded pole, the centre only eighth), so the peak load
        // is a contact-crushing load and not a tensile strength at all. Cutting
        // two chords gives a real flat bearing, the load spreads over it, and
        // the crack initiates at the centre as the test requires. Literature
        // loading angles are 20-30 deg.
        discFlat_ = cfg_.getd("discFlattenDeg", 0.0);
        // End-of-test handling of the brazilian (see FdemSolver.hpp). Default
        // false = the run length is T, exactly as before.
        brazStop_ = cfg_.getb("brazilianStopAfterPeak", false);
        eGaugeLo_ = cfg_.getd("elasticGaugeLo", 0.3);
        eGaugeHi_ = cfg_.getd("elasticGaugeHi", 0.8);
        if (!(eGaugeHi_ > eGaugeLo_ && eGaugeLo_ > 0.0))
            throw std::runtime_error("elasticGaugeLo/Hi must satisfy "
                                     "0 < Lo < Hi (band in units of ft)");
        brazStopDelay_ = cfg_.getd("brazilianStopDelay", 5e-5);
        if (discFlat_ < 0.0 || discFlat_ > 60.0)
            throw std::runtime_error("discFlattenDeg (FULL loading angle of the "
                                     "flattened brazilian disc) must be in "
                                     "[0, 60]");
        if (cfg_.gets("absorbing", "none") != "none")
            throw std::runtime_error("geometry = disc has no box faces to put "
                                     "Lysmer boundaries on: set absorbing = none");
    }

    // loading mode of the uniaxial test — read BEFORE buildMesh, which is
    // where the grip flags are stamped
    if (scen_ == Scenario::TENSION) {
        std::string ld = cfg_.gets("loading", "grips");
        if (ld != "grips" && ld != "platens")
            throw std::runtime_error("loading must be grips | platens (got '"
                                     + ld + "')");
        tensionPlatens_ = ld == "platens";
        if (tensionPlatens_ && cfg_.getd("pullV", 0.05) > 0.0)
            throw std::runtime_error("loading = platens can only PUSH: use "
                                     "pullV < 0 (compression). Direct tension "
                                     "needs glued grips.");
        // ---- metrologie de la compression (voir FdemSolver.hpp) -----------
        ucsStop_ = cfg_.getb("ucsStopAfterPeak", false);
        ucsStopDelay_ = cfg_.getd("ucsStopDelay", 5e-5);
        gLoFrac_ = cfg_.getd("gaugeLoFrac", 0.25);
        gHiFrac_ = cfg_.getd("gaugeHiFrac", 0.75);
        if (!(gHiFrac_ > gLoFrac_ && gLoFrac_ >= 0.0 && gHiFrac_ <= 1.0))
            throw std::runtime_error("gaugeLoFrac/gaugeHiFrac must satisfy "
                                     "0 <= Lo < Hi <= 1 (extensometer bands, "
                                     "fractions of the specimen height)");
    }

    buildMesh();

    // ---- per-phase element tables (elementForces hot loop) ------------------
    DmP_.clear(); nuP_.clear(); crushCapP_.clear(); ftP_.clear(); rhoP_.clear();
    for (const Material& m : phases_.mat) {
        DmP_.push_back(m.Dmat());
        nuP_.push_back(m.nu);
        crushCapP_.push_back(cfg_.getd("crushCap", 8.0 * m.cohesion));
        ftP_.push_back(m.ft);
        rhoP_.push_back(m.rho);
    }

    // ---- optional bulk constitutive law -------------------------------------
    if (cfg_.has("law")) {
        if (phases_.n() > 1)
            throw std::runtime_error("'law' (bulk constitutive law) is a SINGLE "
                "material model: it cannot be combined with mineral 'phases'. "
                "Drop one of the two.");
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
        // C2 (audit 2026-08-11, corrige 2026-08-15) : centroide INITIAL de
        // l'element, requis par le hash spatial des tirages de Weibull de
        // dpdfh (z = 0 en 2D). Sans lui tous les elements tirent la meme
        // resistance et l'heterogeneite DFH est morte.
        for (auto& e : el_) {
            Eigen::Vector2d c = (X0_[e.n[0]] + X0_[e.n[1]] + X0_[e.n[2]]) / 3.0;
            e.st.x0 = Eigen::Vector3d(c.x(), c.y(), 0.0);
        }
        std::cout << "[FDEM] bulk law = " << law_->name()
                  << " (plane strain), elements " << el_.size()
                  << ", lc max = " << lcMax << " m\n";
    }

    // ---- cohesive joint law (PER JOINT) -------------------------------------
    // Intrinsic penalty: p = factor * E / h. The glued assembly then has the
    // series compliance 1/E_eff ~ 1/E + O(1)/(p h): factor 20 keeps the
    // artificial softening at the few-percent level while costing only
    // sqrt(20/100) of the dt a factor-100 penalty would. Intra-grain joints
    // carry the phase material; grain-boundary joints the attenuated mean of
    // the two phases (assignJointProps).
    //
    // insertion = adaptive switches to the EXTRINSIC scheme of Yan, Zheng &
    // Wang (IJRMMS 169, 2023, 105439): no joint exists at t = 0 (bonded,
    // handled by rigid node binding = exact shared-node FEM), each is
    // activated when the edge-averaged traction reaches the envelope.
    {
        std::string ins = cfg_.gets("insertion", "intrinsic");
        if (ins != "intrinsic" && ins != "adaptive")
            throw std::runtime_error("insertion must be intrinsic | adaptive "
                                     "(got '" + ins + "')");
        adaptive_ = ins == "adaptive";
    }
    // jointSoftening = linear (default, unchanged) | yan — the exponential
    // reduction factor f(D) of the article (its eq. 11), see YanSoftening.hpp.
    // Read BEFORE assignJointProps(): it sets the critical opening/slip from
    // the fracture energies through I = int_0^1 f(D) dD.
    {
        std::string js = cfg_.gets("jointSoftening", "linear");
        if (js != "linear" && js != "yan" && js != "munjiza")
            throw std::runtime_error("jointSoftening must be linear | yan "
                                     "| munjiza (got '" + js + "') — munjiza "
                                     "est un ALIAS de yan : Yan et al. 2023 "
                                     "ont repris la z-curve f(D) de Munjiza "
                                     "2004 (a=0,63, b=1,8, c=6), celle de "
                                     "Y-Geo/Solidity");
        yanSoft_ = js == "yan" || js == "munjiza";
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
        std::cout << "[FDEM] joint softening: Yan et al. f(D), a = "
                  << yanP_.a << ", b = " << yanP_.b << ", c = " << yanP_.c
                  << ", int f(D) dD = " << yanI_ << "\n";
    }
    // jointShearUnload = plastic (defaut, inchange) | origin — l'eq. 18 de
    // Yan et al. : decharge ET recharge en cisaillement sur la SECANTE A
    // L'ORIGINE passant par (s_max, tau_env(s_max)), symetrique exact de
    // l'eq. 17 du mode I. Voir l'en-tete pour la mise en garde.
    {
        std::string ju = cfg_.gets("jointShearUnload", "plastic");
        if (ju != "plastic" && ju != "origin")
            throw std::runtime_error("jointShearUnload must be plastic | origin "
                                     "(got '" + ju + "')");
        shearOrigin_ = ju == "origin";
    }
    {
        std::string sr = cfg_.gets("jointShearRange", "cohesion");
        if (sr != "cohesion" && sr != "coulomb")
            throw std::runtime_error("jointShearRange must be cohesion | "
                                     "coulomb (cohesion = plage 3 GfII/c "
                                     "figee, defaut historique ; coulomb = "
                                     "plage divisee par fs = c + tan(phi)"
                                     "|sigma_n| a la pression courante, la "
                                     "convention Solidity Y3Dfd.c l. 1126)");
        shearRangeCoulomb_ = sr == "coulomb";
        if (shearRangeCoulomb_ && !shearOrigin_)
            throw std::runtime_error("jointShearRange = coulomb exige "
                                     "jointShearUnload = origin (le moteur "
                                     "(smax - sE)/plage est celui de la "
                                     "branche origin)");
        if (shearRangeCoulomb_)
            std::cout << "[FDEM] jointShearRange = coulomb : plage de mode II"
                         " divisee par fs(sigma_n) a chaque pas (plancher 2 "
                         "sE)\n";
    }
    if (shearOrigin_) {
        std::cout << "[FDEM] shear unloading: origin secant (Yan eq. 18)\n";
        if (!yanFricScaled_)
            std::cout << "[FDEM] WARNING: jointShearUnload = origin with "
                         "jointFrictionScaled = 0 — the Coulomb term rides the "
                         "origin secant, so frictional sliding is REVERSIBLE "
                         "(no hysteresis loop). Literal article form is "
                         "origin + jointFrictionScaled = 1.\n";
    }
    assignJointProps();
    if (adaptive_) {
        for (auto& J : jt_) J.bonded = true;
        buildBindingTables();
        std::cout << "[FDEM] adaptive insertion: " << jt_.size()
                  << " bonded edges, activation penalty "
                  << cfg_.getd("insertionPenaltyFactor", 4.0) << " E/h\n";
    }
    xiJ_ = cfg_.getd("jointXi", 0.05);

    // ---- E9 (2026-08-20) : gcCell et gcBoxMesh etaient INERTES hors SHPB.
    // Ils n'etaient lus que dans la branche du montage barre-disque-barre, si
    // bien que tout autre cas subissait la grille par defaut : cellule 2 hmin
    // sur le domaine AGRANDI DE MOITIE. Mesure du 20/08 sur le forage
    // d'AbuAisha (8 x 8 m, hmin 3,86 mm) : 6,5 Go de grille et un run qui
    // stagne. Ces cles sont maintenant lues pour TOUS les scenarios ; les
    // defauts hors SHPB reproduisent exactement le comportement anterieur,
    // donc rien ne bouge tant qu'on ne les pose pas.
    if (scen_ != Scenario::SHPB) {
        gcCell_ = cfg_.getd("gcCell", 0.0);        // 0 = 2 hmin, inchange
        gcBoxMesh_ = cfg_.getb("gcBoxMesh", false);// false = boite large, idem
    }
    kp_ = phases_.maxE() * thk_;                       // tool contact penalty
    // ---- A1 : loi de contact de l'outil (voir FdemSolver.hpp) -------------
    {
        std::string tc = cfg_.gets("toolContact", "penalty");
        if (tc != "penalty" && tc != "signorini")
            throw std::runtime_error("toolContact must be penalty | signorini");
        toolSig_ = (tc == "signorini");
        toolSigRelax_ = cfg_.getd("toolSignoriniRelax", 0.0);
        if (!(toolSigRelax_ >= 0.0 && toolSigRelax_ <= 1.0))
            throw std::runtime_error("toolSignoriniRelax must be in [0, 1]");
        if (toolSig_)
            std::cout << "[FDEM] contact outil : SIGNORINI en vitesse "
                         "(CD-Lagrange, Delassus diagonal par noeud) — la "
                         "penalite kp = " << kp_ << " N/m ne s'applique plus "
                         "et SORT du budget du pas de temps ; rattrapage "
                         "d'interpenetration = " << toolSigRelax_
                      << " (0 = condition de vitesse pure)\n";
    }
    // General (debris) contact: node-segment penalty on ROTATING faces is a
    // follower force — at full E*t stiffness it pumps energy into the
    // crushed zone. Debris does not need bulk stiffness: soften and damp.
    kpGC_ = cfg_.getd("gcPenaltyFactor", 0.01) * phases_.maxE() * thk_;
    xiGC_ = cfg_.getd("gcXi", 0.8);
    gcRest_ = cfg_.getd("gcRestitution", 0.2);
    // REPARATION (2026-08-28) : defaut bruyant desormais — miroir du 3D.
    std::cout << "[FDEM] gcRestitution = " << gcRest_
              << (cfg_.has("gcRestitution") ? " (deck)" : " (DEFAUT)")
              << " : detente normale du contact general a ce facteur — a "
                 "figer au deck des que e ou l ejection est une metrique\n";
    // A' : voir FdemSolver.hpp pour la justification. Opt-in, defaut legacy =
    // bit-identique.
    gcEager_ = cfg_.gets("gcSurfaceRefresh", "legacy") == "eager";
    // EPFL (arXiv:2511.14323) sec. 4 : k- = k+(d). Voir FdemSolver.hpp.
    jcAdaptive_ = cfg_.gets("jointContactPenalty", "fixed") == "adaptive";
    if (jcAdaptive_)
        std::cout << "[FDEM] jointContactPenalty = adaptive : k- = k+(D) = "
                     "(1-D) pj — continuite en dn = 0 (EPFL arXiv:2511.14323, "
                     "sec. 4). DIAGNOSTIC : leurs auteurs previennent que "
                     "l'interpenetration croit avec D et que les statistiques "
                     "de fragments cessent d'etre physiques.\n";
    if (gcEager_)
        std::cout << "[FDEM] gcSurfaceRefresh = eager : les faces liberees "
                     "entrent au pas ou leur joint se separe (estampille "
                     "nDead_, pas de grille % 8)\n";
    // contact = penalty (defaut, inchange) | potential — A3 : le contact
    // general par POTENTIEL de Munjiza (eq. 2-5 de l'article), paires
    // d'elements, conservatif, voir PotentialContact.hpp et l'en-tete.
    {
        std::string cm = cfg_.gets("contact", "penalty");
        if (cm != "penalty" && cm != "potential")
            throw std::runtime_error("contact must be penalty | potential "
                                     "(got '" + cm + "')");
        contactPot_ = cm == "potential";
    }
    if (contactPot_) {
        potP_ = cfg_.getd("potPenaltyFactor", 1.0) * phases_.maxE() * thk_;
        potKt_ = cfg_.getd("potTangentFactor", 1.0) * phases_.maxE() * thk_;
        potXi_ = cfg_.getd("potXi", 0.0);
        if (potXi_ < 0.0)
            throw std::runtime_error("potXi is a critical-damping fraction "
                                     ">= 0 (0 = none = the conservative "
                                     "Munjiza potential as published)");
        if (potXi_ > 0.0)
            std::cout << "[FDEM] contact potentiel : amortissement normal "
                         "potXi = " << potXi_ << " (masse reduite des deux "
                         "elements). Le potentiel nu est CONSERVATIF : sans "
                         "ce terme un nuage de debris ne se pose jamais. "
                         "Strictement dissipatif, compte dans gcWork_.\n";
        std::cout << "[FDEM] contact: Munjiza potential (eq. 2-5), p = "
                  << potP_ << " N/m, kt = " << potKt_ << " N/m\n";
    }
    // gcActivation = full (defaut, inchange) | adaptive — voir l'en-tete :
    // seules les faces qui PEUVENT toucher sont balayees (Fukuda et al.).
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
        std::cout << "[FDEM] contact activation: adaptive (Fukuda) — margin "
                  << gcActMargin_ << " cells, sweep <= every "
                  << gcActEvery_ << " steps\n";
    }
    relax_ = 1.0;   // set after dt is known (init order): see below
    muC_ = cfg_.getd("contactMu", 0.5);
    // ---- WP6 : mu de contact residuel post-pulverisation (miroir du 3D,
    // voir Fdem3dSolver.cpp et specs/005-impact-insert-yang/WP6_...) -------
    muCRes_ = cfg_.getd("contactResidualMu", -1.0);
    if (muCRes_ >= 0.0) {
        if (!bdOn_)
            throw std::runtime_error(
                "contactResidualMu exige bulkDamage = yang : sans source "
                "d endommagement volumique aucun element n est jamais "
                "pulverise, et la cle serait lue mais INOPERANTE (regle "
                "E3/E6 : un reglage qui ne fait rien est interdit)");
        if (muCRes_ > muC_)
            std::cout << "\n[FDEM] *** AVERTISSEMENT *** contactResidualMu"
                         " = " << muCRes_ << " > contactMu = " << muC_
                      << " : le frottement AUGMENTE a la pulverisation. "
                         "C est l inverse du modele de Yang et al. 2026 "
                         "(0,18 residuel contre 0,6 intact). Voulu ?\n\n";
        std::cout << "[FDEM] contactResidualMu = " << muCRes_
                  << " : toute interaction de contact (outil ou general) "
                     "impliquant un element pulverise (bulkDamage : D = "
                  << bdDmax_ << ") bascule de contactMu = " << muC_
                  << " a ce mu residuel (sliding friction post-rupture, "
                     "Yang et al. 2026)\n";
    }
    xiC_ = cfg_.getd("contactXi", 0.05);
    vReg_ = cfg_.getd("contactVreg", 1e-3);

    placeTool();
    setupBoundaries();
    if (scen_ == Scenario::SHPB) setupShpbGauges();
    setupConfinement();
    setupHydro();
    setupExcavation();                     // no-op si excavRelease = false
    setupBrazilianLoad();                  // before computeStableDt: it sets
                                           // kpPlaten_, which enters the budget
    setupStrainGauge();                    // after the platens: needs plTop_.y
    // ---- bulkViscosityXi : donner la viscosite en TAUX D AMORTISSEMENT ---
    // Yan et al. publient mu = 7,6e3 Pa.s pour E = 15 GPa, rho = 1704 et un
    // maillage h = 0,75 mm. Transporter ce chiffre tel quel vers un autre
    // materiau ou un autre maillage n a pas de sens : ce qui est invariant,
    // c est le taux d amortissement VU PAR UN ELEMENT.
    //
    // Pour sigma = E eps + 2 mu epsdot, le facteur de perte a la pulsation
    // propre de l element omega = c/h vaut 2 mu omega / E, soit
    //     xi = mu / (h sqrt(E rho))   et donc   mu = xi h sqrt(E rho).
    //
    // CE QUE VAUT xi = 2. L amortissement CRITIQUE de Munjiza / Y-Geo est
    // mu_crit = 2 h sqrt(E rho), soit exactement xi = 2. Et ce n est pas une
    // analogie : applique aux chiffres de la Table 1 de Yan (h = 0,75 mm,
    // E = 15 GPa, rho = 1704), mu_crit = 7583 Pa.s contre les 7600 publies —
    // 0,2 % d ecart. Leur viscosite EST l amortissement critique de Munjiza,
    // ce que leur article ne dit pas (il renvoie a Tatone & Grasselli sans
    // aucune formule ni etude de sensibilite). bulkViscosityXi = 2 reproduit
    // donc leur choix sur n importe quel maillage et n importe quelle roche.
    //
    // POURQUOI PAS LE TEMPS DE RELAXATION. Conserver tau = 2 mu / E donnerait,
    // pour un granite E = 48,26 GPa maille a 0,238 mm, mu = 24 450 Pa.s soit
    // xi = 9,2 : quatre fois et demie l amortissement critique, et la borne
    // diffusive du pas de temps deviendrait franchement contraignante. La
    // grandeur transportable est le taux d amortissement A L ECHELLE DE LA
    // MAILLE, pas un temps de relaxation materiau — parce que ce terme est un
    // filtre de bruit de maillage, pas une viscosite mesuree sur la roche.
    {
        double xiV = cfg_.getd("bulkViscosityXi", 0.0);
        if (xiV < 0.0)
            throw std::runtime_error("bulkViscosityXi must be >= 0 "
                                     "(taux d amortissement a l echelle de "
                                     "l element ; Yan et al. valent 2,00)");
        if (xiV > 0.0 && bulkVisc_ > 0.0)
            throw std::runtime_error("bulkViscosity et bulkViscosityXi sont "
                                     "exclusives : la seconde CALCULE la "
                                     "premiere a partir du maillage");
        if (xiV > 0.0) {
            std::vector<double> mus(el_.size());
            for (std::size_t eI = 0; eI < el_.size(); ++eI)
                mus[eI] = xiV * hEl_[eI]
                        * std::sqrt(phases_.mat[el_[eI].phase].E
                                    * rhoP_[el_[eI].phase]);
            std::vector<double> srt = mus;
            std::sort(srt.begin(), srt.end());
            bulkVisc_ = srt.empty() ? 0.0 : srt[srt.size() / 2];
            double xLo = 1e30, xHi = 0.0;
            for (std::size_t eI = 0; eI < el_.size(); ++eI) {
                double x = mus[eI] > 0.0 ? xiV * bulkVisc_ / mus[eI] : 0.0;
                xLo = std::min(xLo, x); xHi = std::max(xHi, x);
            }
            std::cout << "[FDEM] bulkViscosityXi = " << xiV << " (soit "
                      << 0.5 * xiV << " x l amortissement CRITIQUE de "
                         "Munjiza mu = 2 h sqrt(E rho) ; Yan et al. Table 1 "
                         "= 2,00 = critique) -> mu = "
                      << bulkVisc_ << " Pa.s, pose sur la maille MEDIANE. "
                         "Sur ce maillage le xi effectif va de " << xLo
                      << " (grosses mailles) a " << xHi << " (fines)"
                      << (xHi > 1.05 * xLo
                          ? " : le maillage est gradue, l amortissement"
                            " effectif l est donc aussi." : " (uniforme).")
                      << "\n";
        }
    }
    computeStableDt();

    relax_ = std::exp(-dt_ / cfg_.getd("gcBirthTau", 1e-6));
    srRelax_ = srTau_ > 0.0 ? std::exp(-dt_ / srTau_) : 0.0;
    if (difOn_) {
        std::cout << "[FDEM] strainRateDIF = yang (Yang et al., IJRMMS 191, "
                     "2025, eq. 2-3) : ft et Gf recoivent DIF_traction, "
                     "cohesion et GfII recoivent DIF_compression, l angle de "
                     "frottement est inchange. Facteurs FIGES a l insertion "
                     "du joint.\n"
                  << "[FDEM]   taux de deformation = principale max de "
                     "sym(Fdot F^-1), filtre a strainRateTau = " << srTau_
                  << " s (" << (srTau_ / dt_) << " pas)\n"
                  << "[FDEM]   exposant de DIF_traction = " << difExpT_
                  << (difExpT_ > 0.1 ? " (deduit de LEUR figure 2b : la loi "
                                       "se raccorde alors a ses deux bornes)"
                                     : " (transcription LITTERALE de leur "
                                       "eq. 3)")
                  << " : la loi va de "
                  << difTensionYang(6.0e-6, difExpT_) << " (juste au-dessus "
                     "de 5e-6 /s, contre 1 en dessous) a "
                  << difTensionYang(99.0, difExpT_) << " (juste sous 1e2 /s, "
                     "contre 1,85 au-dessus)\n";
        if (difExpT_ < 0.1)
            std::cout << "[FDEM]   AVERTISSEMENT : avec l exposant litteral "
                         "0,07 ces deux paires ne se raccordent PAS (sauts de "
                         "12 % et 22 %). Un joint dont le taux traverse 1e2 /s "
                         "voit sa resistance sauter de 22 % sans physique. "
                         "strainRateDIF = yang-fig2 evite cela.\n";
        if (bulkVisc_ > 0.0)
            std::cout << "[FDEM]   AVERTISSEMENT : bulkViscosity est arme EN "
                         "MEME TEMPS. La contrainte visqueuse 2 mu D entre "
                         "dans sigG, donc dans le critere d insertion : le "
                         "taux de deformation agit alors DEUX FOIS, une fois "
                         "en gonflant la contrainte d essai et une fois en "
                         "gonflant le seuil. C est fidele a Yan (leur eq. 6 "
                         "EST la contrainte que lit leur eq. 7) mais ce n est "
                         "pas le modele de Yang. A separer avant toute "
                         "calibration.\n";
    }
    if (scen_ == Scenario::TENSION || scen_ == Scenario::BRAZILIAN)
        pullV_ = cfg_.getd("pullV", 0.05);
    gripFree_ = cfg_.getb("gripLateralFree", false);
    pullRamp_ = cfg_.getd("pullRamp", 0.0);
    pullDelay_ = cfg_.getd("pullDelay", 0.0);
    if (pullDelay_ < 0.0)
        throw std::runtime_error("pullDelay must be >= 0 s (it delays the "
                                 "start of the axial loading so a confining "
                                 "pressure can equilibrate first)");
    fragId_.assign(el_.size(), 0);
    toolKE0_ = tool_.ke();

    // Meme garde-fou qu'en 3D : la grille en croix diverge en phase debris
    // (mesure 2026-08-06 : 379 joints, Ec > energie livree) la ou le Voronoi
    // desordonne est propre. Condition de validite, pas d'esthetique.
    if (!voronoi_ && (scen_ == Scenario::PERCUSSION || scen_ == Scenario::SHEAR))
        std::cout << "[FDEM] WARNING: mesh = grid + scenario de fissuration — "
                     "un maillage structure biaise les trajets et peut "
                     "diverger en phase debris (FICHE 2026-08-06). Utiliser "
                     "mesh = voronoi + grainSeeding = random pour tout "
                     "resultat ou la casse compte.\n";
    std::cout << "[FDEM] " << el_.size() << " elements, " << jt_.size()
              << " joints, " << X0_.size() << " nodes, dt = " << dt_
              << " s, steps = " << (long)std::ceil(T_ / dt_) << "\n";
    if (disc_)
        std::cout << "[FDEM] disc geometry: diameter " << 2.0 * discR_
                  << " m, centre (" << discC_.x() << ", " << discC_.y() << ")\n";
    if (voronoi_) {
        long nGB = 0;
        for (const auto& J : jt_) if (J.type > 0) ++nGB;
        std::cout << "[FDEM] voronoi: " << nGrains_ << " grains, "
                  << phases_.n() << " phase(s), " << nGB
                  << " grain-boundary joints, hmin = " << hmin_ << " m\n";
    }
}

// ---------------------------------------------------------------------------
// Mesh generation. Two front-ends share one topology builder:
//   mesh = grid    — structured cross-diagonal mesh with optional corner
//                    jitter (the original layout; one implicit "grain").
//   mesh = voronoi — Voronoi grain structure (Tessellation) with per-grain
//                    mineral phases: the GBM mode. Crack paths then follow
//                    grain boundaries and intra-grain fans instead of the
//                    lattice directions.
// Every interior edge gets a 4-node cohesive joint; exterior faces are kept
// for the quiet boundaries and the general contact.
// ---------------------------------------------------------------------------
void FdemSolver::buildMesh() {
    std::string mesh = cfg_.gets("mesh", "grid");
    if (mesh != "grid" && mesh != "voronoi" && mesh != "file")
        throw std::runtime_error("mesh must be grid | voronoi | file (got '"
                                 + mesh + "')");
    voronoi_ = mesh == "voronoi";
    if (mesh == "file") {
        if (shpb_ || disc_)
            throw std::runtime_error("mesh = file is implemented for the BOX "
                "geometry (percussion / shear / tension); the disc and shpb "
                "assemblies build their own meshes");
        if (phases_.n() > 1)
            throw std::runtime_error("'phases' needs the grain machinery: "
                "mesh = file imports a single-material unstructured mesh");
        buildMeshFile();
        return;
    }
    if (shpb_) { buildMeshShpb(); return; }            // three-body assembly
    if (!voronoi_ && phases_.n() > 1)
        throw std::runtime_error("'phases' declares "
            + std::to_string(phases_.n()) + " minerals but mesh = grid would "
            "silently use only the first: set mesh = voronoi (or drop the "
            "phases key)");
    // discMesh = cut (default: disc cut out of the box mesh, rim ragged at
    // element scale) | native (the disc is MESHED as a disc: boundary ring
    // exactly on the circle). cut is kept as the default so every earlier
    // result stays reproducible bit for bit.
    std::string dm = cfg_.gets("discMesh", "cut");
    if (dm != "cut" && dm != "native")
        throw std::runtime_error("discMesh must be cut | native (got '" + dm
                                 + "')");
    if (dm == "native") {
        if (!disc_)
            throw std::runtime_error("discMesh = native needs geometry = disc");
        voronoi_ = true;                   // grain machinery (GBM joints) on
        buildMeshDisc();
        return;
    }
    if (voronoi_) { buildMeshVoronoi(); return; }

    double dx = W_ / nx_, dy = H_ / ny_;
    hmin_ = std::min(dx, dy);
    double jit = cfg_.getd("meshJitter", 0.0) * 0.5 * hmin_;
    std::mt19937 rng(cfg_.geti("seed", 12345));
    std::uniform_real_distribution<double> U(-jit, jit);

    int vnx = nx_ + 1, vny = ny_ + 1;
    std::vector<Eigen::Vector2d> vpos(vnx * vny + nx_ * ny_);
    for (int j = 0; j < vny; ++j)
        for (int i = 0; i < vnx; ++i) {
            Eigen::Vector2d p(i * dx, j * dy);
            if (jit > 0 && i > 0 && i < nx_ && j > 0 && j < ny_)
                p += Eigen::Vector2d(U(rng), U(rng));
            vpos[j * vnx + i] = p;
        }
    int centerBase = vnx * vny;
    auto vid = [&](int i, int j) { return j * vnx + i; };
    for (int j = 0; j < ny_; ++j)
        for (int i = 0; i < nx_; ++i)
            vpos[centerBase + j * nx_ + i] =
                0.25 * (vpos[vid(i, j)] + vpos[vid(i + 1, j)]
                        + vpos[vid(i + 1, j + 1)] + vpos[vid(i, j + 1)]);

    std::vector<std::array<int, 3>> tris;
    tris.reserve((std::size_t)4 * nx_ * ny_);
    for (int j = 0; j < ny_; ++j)
        for (int i = 0; i < nx_; ++i) {
            int c00 = vid(i, j), c10 = vid(i + 1, j);
            int c11 = vid(i + 1, j + 1), c01 = vid(i, j + 1);
            int cc = centerBase + j * nx_ + i;
            tris.push_back({c00, c10, cc});
            tris.push_back({c10, c11, cc});
            tris.push_back({c11, c01, cc});
            tris.push_back({c01, c00, cc});
        }
    std::vector<int> triGrain(tris.size(), 0);         // one implicit grain
    nGrains_ = 1;
    cutDisc(vpos, tris, triGrain);
    buildFromTriangles(vpos, tris, triGrain, {0});
}


// ---------------------------------------------------------------------------
// geometry = shpb — the SHPB assembly of Yan, Zheng & Wang (2023) fig. 22:
// incident bar (length LIB) | brazilian disc (D_rock) | transmission bar
// (LTB), all of height D_bar, meshed in ONE pass as THREE DISJOINT bodies.
//
// Why disjoint bodies and not one welded mesh: in the real test the bars only
// PRESS on the rock. The interfaces must transmit compression and friction
// (k_c = 0.577 of table 2) and must be able to separate when the reflected
// wave unloads them, which is exactly what the existing node-to-edge general
// contact does. Welding them would make the assembly one continuum: no
// reflection at the impedance jump, no separation, and the transmitted wave
// would be wrong from the first microsecond.
//
// Why the two materials come from mineral PHASES: the phase machinery already
// carries a full Material per phase through the element tables (DmP_, rhoP_),
// the mass lumping, the joint properties and the boundary impedances. One
// "grain" per body plus a phase per grain therefore gives the bars E = 240 GPa,
// rho = 7700, nu = 0.01 and the disc E = 100.8 GPa, rho = 2800, nu = 0.297
// with no new material plumbing at all. The bodies being disjoint, no
// grain-boundary joint is ever created between two phases.
//
// The bars are meshed with the structured cross-diagonal grid at
// shpbBarElemSize (5 mm by default: 40 elements per pulse wavelength
// c*tau = 1.2 m, and the bars only have to carry a 1D elastic wave), the disc
// with the same ring + hex-fill + Delaunay generator as discMesh = native at
// shpbDiscElemSize (0.75 mm, the article's mesh size). Keeping the bars coarse
// is what makes the run affordable: the stable dt is fixed by the FINEST
// element, which is in the disc, and refining the bars would only multiply the
// element count at constant dt.
// ---------------------------------------------------------------------------
void FdemSolver::buildMeshShpb() {
    std::mt19937 rng(cfg_.geti("seed", 12345));
    std::vector<Eigen::Vector2d> vpos;
    std::vector<std::array<int, 3>> tris;
    std::vector<int> triGrain;

    // ---- a rectangular body on the cross-diagonal grid ---------------------
    auto addBar = [&](double x0, double x1, double h, int grain) {
        int nxb = std::max(1, (int)std::llround((x1 - x0) / h));
        int nyb = std::max(1, (int)std::llround(H_ / h));
        double dx = (x1 - x0) / nxb, dy = H_ / nyb;
        int base = (int)vpos.size();
        int vnx = nxb + 1;
        for (int j = 0; j <= nyb; ++j)
            for (int i = 0; i <= nxb; ++i)
                vpos.push_back({x0 + i * dx, j * dy});
        int cBase = (int)vpos.size();
        for (int j = 0; j < nyb; ++j)
            for (int i = 0; i < nxb; ++i)
                vpos.push_back({x0 + (i + 0.5) * dx, (j + 0.5) * dy});
        auto vid = [&](int i, int j) { return base + j * vnx + i; };
        for (int j = 0; j < nyb; ++j)
            for (int i = 0; i < nxb; ++i) {
                int c00 = vid(i, j), c10 = vid(i + 1, j);
                int c11 = vid(i + 1, j + 1), c01 = vid(i, j + 1);
                int cc = cBase + j * nxb + i;
                tris.push_back({c00, c10, cc});
                tris.push_back({c10, c11, cc});
                tris.push_back({c11, c01, cc});
                tris.push_back({c01, c00, cc});
                for (int k = 0; k < 4; ++k) triGrain.push_back(grain);
            }
    };

    // ---- the disc: boundary ring on the circle + jittered hex fill ---------
    auto addDisc = [&](double h, int grain) {
        double R = discR_;
        const Eigen::Vector2d C = discC_;
        int base = (int)vpos.size();
        std::vector<Eigen::Vector2d> p;
        // radial dither of 1e-7 R breaks the exact cocircularity of the ring
        // (Bowyer-Watson's inCircle test is degenerate on it) — same guard as
        // buildMeshDisc, and three orders below the chord sagitta.
        std::uniform_real_distribution<double> Ud(0.2, 1.0);
        int nR = std::max(8, (int)std::ceil(2.0 * M_PI * R / h));
        for (int k = 0; k < nR; ++k) {
            double t = 2.0 * M_PI * k / nR;
            double s = 1.0 - 1e-7 * Ud(rng);
            p.push_back(C + R * s * Eigen::Vector2d(std::cos(t), std::sin(t)));
        }
        double jit = cfg_.getd("meshJitter", 0.25);
        std::uniform_real_distribution<double> U(-jit * h, jit * h);
        double dy = h * std::sqrt(3.0) / 2.0;
        int j = 0;
        for (double y = C.y() - R + 0.7 * h; y <= C.y() + R - 0.7 * h;
             y += dy, ++j)
            for (double x = C.x() - R + ((j & 1) ? 0.5 * h : 0.0);
                 x <= C.x() + R; x += h) {
                Eigen::Vector2d q(x + U(rng), y + U(rng));
                if ((q - C).norm() <= R - 0.7 * h) p.push_back(q);
            }
        std::vector<std::array<int, 3>> dt = delaunayCCW(p);
        if (dt.empty())
            throw std::runtime_error("geometry = shpb: disc triangulation "
                                     "failed");
        // ---- Laplacian smoothing of the INTERIOR points ---------------------
        // The stable time step is set by the SMALLEST inscribed diameter
        // anywhere in the assembly, and the ring/fill transition of a jittered
        // hex fill routinely leaves slivers: measured on the first build,
        // h_min = 0.199 mm for a nominal 0.75 mm mesh, i.e. dt divided by ~13
        // for a handful of degenerate elements. A few Laplacian sweeps (each
        // interior vertex to the mean of its Delaunay neighbours, ring frozen
        // so the disc stays a disc) regularise them at negligible cost. The
        // topology is NOT re-triangulated: smoothing a Delaunay mesh of a
        // convex domain cannot invert a triangle here, and the areas are
        // checked below.
        int nSm = cfg_.geti("shpbDiscSmooth", 8);
        for (int it = 0; it < nSm; ++it) {
            std::vector<Eigen::Vector2d> acc(p.size(), Eigen::Vector2d::Zero());
            std::vector<int> cnt(p.size(), 0);
            for (const auto& t : dt)
                for (int k = 0; k < 3; ++k) {
                    int a = t[k], b = t[(k + 1) % 3];
                    acc[a] += p[b]; ++cnt[a];
                    acc[b] += p[a]; ++cnt[b];
                }
            for (std::size_t k = (std::size_t)nR; k < p.size(); ++k)
                if (cnt[k] > 0) p[k] = 0.7 * (acc[k] / cnt[k]) + 0.3 * p[k];
        }
        if (nSm > 0) {
            dt = delaunayCCW(p);                       // re-triangulate once
            if (dt.empty())
                throw std::runtime_error("geometry = shpb: re-triangulation "
                                         "after smoothing failed");
        }
        for (const auto& q : p) vpos.push_back(q);
        for (const auto& t : dt) {
            tris.push_back({base + t[0], base + t[1], base + t[2]});
            triGrain.push_back(grain);
        }
    };

    // phase lookup by NAME so the order of the `phases` line does not matter
    auto phaseIdx = [&](const char* key, const char* def) {
        std::string want = cfg_.gets(key, def);
        for (int i = 0; i < phases_.n(); ++i)
            if (phases_.name[i] == want) return i;
        if (phases_.n() == 1) return 0;                // single-material model
        throw std::runtime_error(std::string("geometry = shpb: no phase named '")
            + want + "' in the `phases` list (set " + key + ")");
    };
    int pBar  = phaseIdx("shpbBarPhase", "bar");
    int pRock = phaseIdx("shpbRockPhase", "rock");

    double xDiscL = shpbLIB_ + shpbGap_;
    addBar(0.0, shpbLIB_, hBar_, 0);                   // incident bar
    std::vector<int> grainPhase{pBar};
    if (!shpbNoDisc_) {
        addDisc(hDisc_, 1);
        addBar(xDiscL + shpbDrock_ + shpbGap_, W_, hBar_, 2);
        grainPhase.push_back(pRock);
        grainPhase.push_back(pBar);
        nGrains_ = 3;
    } else {
        // bar-only wave-propagation verification: one bar, the pulse at its
        // left end and the viscous boundary at its right end
        xEndAbs_ = shpbLIB_;
        nGrains_ = 1;
    }
    voronoi_ = true;                                   // per-element h, phases
    buildFromTriangles(vpos, tris, triGrain, grainPhase);
    hmin_ = 1e30;
    for (double h : hEl_) hmin_ = std::min(hmin_, h);
    std::cout << "[FDEM] shpb assembly: incident bar " << shpbLIB_
              << " m, disc " << (shpbNoDisc_ ? 0.0 : 2.0 * discR_)
              << " m, transmission bar " << (shpbNoDisc_ ? 0.0 : shpbLTB_)
              << " m, height " << H_ << " m, " << el_.size()
              << " elements (bar h = " << hBar_ << " m, disc h = " << hDisc_
              << " m), hmin = " << hmin_ << " m\n";
}

// ---------------------------------------------------------------------------
// Monitor points 1 and 2 (fig. 22). A strain gauge is not a point: it averages
// over its own length. Here the reading is the AREA-WEIGHTED mean of the
// co-rotated axial strain eps_xx over every element of the bar whose centroid
// lies within shpbGaugeHalfLength of the monitor abscissa — i.e. a gauge of
// 2 x that length, 4 bar elements by default. Averaging over the full bar
// height also removes the lateral (Pochhammer-Chree) ringing that a
// single-element read would show.
// ---------------------------------------------------------------------------
void FdemSolver::setupShpbGauges() {
    monEl1_.clear(); monEl2_.clear();
    monA1_.clear(); monA2_.clear();
    monArea1_ = monArea2_ = 0.0;
    for (int e = 0; e < (int)el_.size(); ++e) {
        const Elem& E = el_[e];
        double cx = (X0_[E.n[0]].x() + X0_[E.n[1]].x() + X0_[E.n[2]].x()) / 3.0;
        if (E.grain == 0 && std::abs(cx - shpbM1_) <= shpbGaugeW_) {
            monEl1_.push_back(e); monA1_.push_back(E.A0); monArea1_ += E.A0;
        }
        if (E.grain == 2 && std::abs(cx - shpbM2_) <= shpbGaugeW_) {
            monEl2_.push_back(e); monA2_.push_back(E.A0); monArea2_ += E.A0;
        }
    }
    if (monEl1_.empty())
        throw std::runtime_error("geometry = shpb: monitor point 1 caught no "
                                 "element — check shpbMonitor1 / "
                                 "shpbGaugeHalfLength");
    if (!shpbNoDisc_ && monEl2_.empty())
        throw std::runtime_error("geometry = shpb: monitor point 2 caught no "
                                 "element — check shpbMonitor2");
    std::cout << "[FDEM] shpb gauges: monitor 1 at x = " << shpbM1_ << " m ("
              << monEl1_.size() << " elements), monitor 2 at x = " << shpbM2_
              << " m (" << monEl2_.size() << " elements), half-length "
              << shpbGaugeW_ << " m\n";
}

// Prescribed velocity history on the struck face (fig. 22, bottom inset).
// The article fits experimental data it does not publish; the inset shows a
// single smooth compressive pulse peaking at 5.2 m/s about 0.09 ms after the
// start and returning to zero at ~0.22 ms. Two shapes are offered:
//   halfsine  v(t) = V0 sin(pi t / tau)                    (default)
//   trapezoid V0 with a linear rise/fall of (1-plateau)/2 tau each.
// Both are zero outside [0, tau]: the bar is then free and the reflected wave
// travels back through a quiet boundary condition (v prescribed = 0 would
// clamp it and reflect a second time).
double FdemSolver::shpbVel(double t) const {
    if (t <= 0.0 || t >= shpbTau_) return 0.0;
    if (shpbPulse_ == 0)
        return shpbV0_ * std::sin(M_PI * t / shpbTau_);
    double rise = 0.5 * (1.0 - shpbPlateau_) * shpbTau_;
    if (t < rise) return shpbV0_ * t / rise;
    if (t > shpbTau_ - rise) return shpbV0_ * (shpbTau_ - t) / rise;
    return shpbV0_;
}

void FdemSolver::shpbGaugeRead() {
    double s1 = 0.0, s2 = 0.0;
    for (std::size_t k = 0; k < monEl1_.size(); ++k)
        s1 += el_[monEl1_[k]].exx * monA1_[k];
    for (std::size_t k = 0; k < monEl2_.size(); ++k)
        s2 += el_[monEl2_[k]].exx * monA2_[k];
    epsM1_ = monArea1_ > 0.0 ? s1 / monArea1_ : 0.0;
    epsM2_ = monArea2_ > 0.0 ? s2 / monArea2_ : 0.0;
}

// ---------------------------------------------------------------------------
// discMesh = native — the disc is meshed AS a disc, the way Y-Geo/Irazu (and
// Yan et al.'s fig. 11) build their brazilian models, instead of being cut
// out of a box mesh:
//   1. boundary ring EXACTLY on the circle (and on the two flattening chords
//      when discFlattenDeg > 0), spaced at the element size, so the rim error
//      is the chord sagitta h^2/(8R) — micrometres — instead of the
//      half-element staircase of the cut;
//   2. interior fill on a jittered hex lattice, kept clear of the ring by
//      0.7 h so no boundary sliver can form;
//   3. Bowyer-Watson Delaunay of ring + fill. The domain is CONVEX, so every
//      triangle of the triangulation belongs to it: nothing is dropped, no
//      flapping-element pass is needed;
//   4. grain ids by nearest Voronoi seed of the centroid (jittered hex or
//      Poisson at grainSize), so the GBM joint classification — and the
//      per-grain phase draw — work exactly as with the box tessellation.
// The two earlier attempts document why this exists: cutting by centroid
// leaves a ragged rim (the user's objection), and projecting the cut rim
// onto the circle crushed elements into slivers (dt / 8.6, rejected
// 2026-08-05).
// ---------------------------------------------------------------------------
void FdemSolver::buildMeshDisc() {
    double h = cfg_.reqd("grainElemSize");
    double d = cfg_.reqd("grainSize");
    double R = discR_;
    const Eigen::Vector2d C = discC_;
    std::mt19937 rng(cfg_.geti("seed", 12345));

    // flattened disc: chords at |y - yc| = R cos(alpha), alpha = discFlat_/2
    double ca = 1.0;                                   // cos(alpha)
    if (discFlat_ > 0.0) ca = std::cos(0.5 * discFlat_ * M_PI / 180.0);
    double yCut = R * ca;
    auto inside = [&](const Eigen::Vector2d& p, double margin) {
        return (p - C).norm() <= R - margin
            && std::abs(p.y() - C.y()) <= yCut - margin;
    };

    // ---- 1. boundary ring ---------------------------------------------------
    // Every ring point lies on ONE circle, so every ring quadruple is exactly
    // cocircular — the degenerate case of the strict inCircle test, where
    // Bowyer-Watson's flip decisions ride on floating-point noise (measured:
    // the triangulation span, 18+ CPU-minutes without terminating). A
    // deterministic RADIAL dither of 1e-7 R (nanometres — three orders below
    // the chord sagitta this mesher exists to achieve) breaks the exact
    // cocircularity while leaving the rim physically on the circle.
    std::uniform_real_distribution<double> Ud(0.2, 1.0);
    auto dither = [&]() { return 1.0 - 1e-7 * Ud(rng); };
    std::vector<Eigen::Vector2d> vpos;
    if (discFlat_ > 0.0) {
        double a = 0.5 * discFlat_ * M_PI / 180.0;
        double sa = std::sin(a);
        // four arcs boundaries: theta in (a-90 -> 90-a) right, (90+a -> 270-a)
        // left, plus the two chords y = yc +- yCut, x in [-R sa, R sa]
        auto arc = [&](double t0, double t1) {
            int n = std::max(2, (int)std::ceil(R * (t1 - t0) / h));
            for (int k = 0; k < n; ++k) {              // t1 excluded: next
                double t = t0 + (t1 - t0) * k / n;     // segment starts there
                vpos.push_back(C + R * dither()
                               * Eigen::Vector2d(std::cos(t), std::sin(t)));
            }
        };
        auto chord = [&](double y, double x0, double x1) {
            int n = std::max(2, (int)std::ceil(std::abs(x1 - x0) / h));
            for (int k = 0; k < n; ++k)
                vpos.push_back({C.x() + x0 + (x1 - x0) * k / n, C.y() + y});
        };
        arc(a - M_PI / 2, M_PI / 2 - a);               // right arc, CCW
        chord(+yCut, +R * sa, -R * sa);                // top chord, CCW
        arc(M_PI / 2 + a, 3 * M_PI / 2 - a);           // left arc
        chord(-yCut, -R * sa, +R * sa);                // bottom chord
    } else {
        int n = std::max(8, (int)std::ceil(2.0 * M_PI * R / h));
        for (int k = 0; k < n; ++k) {
            double t = 2.0 * M_PI * k / n;
            vpos.push_back(C + R * dither()
                           * Eigen::Vector2d(std::cos(t), std::sin(t)));
        }
    }
    const int nRing = (int)vpos.size();

    // ---- 2. interior fill ---------------------------------------------------
    double jit = cfg_.getd("meshJitter", 0.25);        // fraction of h
    std::uniform_real_distribution<double> U(-jit * h, jit * h);
    double dy = h * std::sqrt(3.0) / 2.0;
    int j = 0;
    for (double y = C.y() - R + 0.7 * h; y <= C.y() + R - 0.7 * h; y += dy, ++j)
        for (double x = C.x() - R + ((j & 1) ? 0.5 * h : 0.0);
             x <= C.x() + R; x += h) {
            Eigen::Vector2d p(x + U(rng), y + U(rng));
            if (inside(p, 0.7 * h)) vpos.push_back(p);
        }

    // ---- 3. Delaunay --------------------------------------------------------
    std::vector<std::array<int, 3>> tris = delaunayCCW(vpos);
    if (tris.empty())
        throw std::runtime_error("discMesh = native: triangulation failed");

    // ---- 4. grains by nearest seed ------------------------------------------
    std::vector<Eigen::Vector2d> seeds;
    {
        std::uniform_real_distribution<double> Us(-0.5 * d * 0.5, 0.5 * d * 0.5);
        double sy = d * std::sqrt(3.0) / 2.0;
        int r = 0;
        for (double y = C.y() - R; y <= C.y() + R; y += sy, ++r)
            for (double x = C.x() - R + ((r & 1) ? 0.5 * d : 0.0);
                 x <= C.x() + R; x += d) {
                Eigen::Vector2d s(x + Us(rng), y + Us(rng));
                if ((s - C).norm() <= R) seeds.push_back(s);
            }
        if (seeds.empty()) seeds.push_back(C);
    }
    std::vector<int> triGrain(tris.size(), 0);
    for (std::size_t t = 0; t < tris.size(); ++t) {
        Eigen::Vector2d cen = (vpos[tris[t][0]] + vpos[tris[t][1]]
                               + vpos[tris[t][2]]) / 3.0;
        double best = 1e300;
        for (int s = 0; s < (int)seeds.size(); ++s) {
            double q = (cen - seeds[s]).squaredNorm();
            if (q < best) { best = q; triGrain[t] = s; }
        }
    }
    nGrains_ = (int)seeds.size();

    // per-grain phase draw, area-greedy like Tessellation: single phase = all 0
    std::vector<int> grainPhase(nGrains_, 0);
    if (phases_.n() > 1) {
        std::vector<int> order(nGrains_);
        for (int g = 0; g < nGrains_; ++g) order[g] = g;
        std::shuffle(order.begin(), order.end(), rng);
        std::vector<double> deficit = phases_.fraction;
        for (int g : order) {
            int bestP = 0;
            for (int p = 1; p < phases_.n(); ++p)
                if (deficit[p] > deficit[bestP]) bestP = p;
            grainPhase[g] = bestP;
            deficit[bestP] -= 1.0 / nGrains_;
        }
    }

    hmin_ = 1e300;
    for (const auto& t : tris) {
        const auto &A = vpos[t[0]], &B = vpos[t[1]], &Cc = vpos[t[2]];
        double det = (B.x() - A.x()) * (Cc.y() - A.y())
                   - (Cc.x() - A.x()) * (B.y() - A.y());
        double per = (B - A).norm() + (Cc - B).norm() + (A - Cc).norm();
        if (det > 0) hmin_ = std::min(hmin_, 2.0 * det / per);
    }
    std::cout << "[FDEM] native disc mesh: " << nRing << " rim nodes on the "
              << (discFlat_ > 0 ? "flattened circle" : "circle") << ", "
              << tris.size() << " triangles, " << nGrains_
              << " grains, hmin = " << hmin_ << " m (rim sagitta "
              << h * h / (8.0 * R) << " m)\n";
    buildFromTriangles(vpos, tris, triGrain, grainPhase);
}

// ---------------------------------------------------------------------------
// geometry = disc: keep only the triangles whose CENTROID lies inside the disc
// of radius min(W, H)/2 centred in the box. Filtering triangles (rather than
// clipping them) keeps every element well shaped — the price is a rim that is
// ragged at element scale. That is the right trade here: the brazilian
// observable is where the crack runs and what the platen reads, neither of
// which is decided by a half-element of rim roughness, whereas a clipped
// sliver at the rim would set the stable time step for the whole run.
//
// Grain ids are left untouched (gaps are harmless: phaseOfGrain is indexed by
// id), only the count reported to the user is refreshed.
// ---------------------------------------------------------------------------
void FdemSolver::cutDisc(std::vector<Eigen::Vector2d>& vpos,
                         std::vector<std::array<int, 3>>& tris,
                         std::vector<int>& triGrain) const {
    if (!disc_) return;
    std::vector<std::array<int, 3>> keptT;
    std::vector<int> keptG;
    keptT.reserve(tris.size());
    keptG.reserve(tris.size());
    double R2 = discR_ * discR_;
    double yCut = discFlat_ > 0.0
                  ? discR_ * std::cos(0.5 * discFlat_ * M_PI / 180.0)
                  : 1e300;                             // no chord
    for (std::size_t k = 0; k < tris.size(); ++k) {
        Eigen::Vector2d c = (vpos[tris[k][0]] + vpos[tris[k][1]]
                             + vpos[tris[k][2]]) / 3.0;
        if ((c - discC_).squaredNorm() > R2) continue;
        if (std::abs(c.y() - discC_.y()) > yCut) continue;   // flattening
        keptT.push_back(tris[k]);
        keptG.push_back(triGrain.empty() ? 0 : triGrain[k]);
    }
    if (keptT.empty())
        throw std::runtime_error("geometry = disc removed every element: check "
                                 "W/H (the disc is inscribed in the box)");
    tris.swap(keptT);
    triGrain.swap(keptG);
    (void)vpos;

    // ---- drop the flapping elements -----------------------------------------
    // A staircase rim leaves triangles hanging by a SINGLE edge: in the
    // cross-diagonal grid the four triangles of a cell have four different
    // centroids, so a rim cell routinely keeps one of them while its three
    // siblings go. Such an element is a flap on one cohesive joint, and that
    // joint fails under any load at all — measured on the 20x20 disc: the first
    // joint broke at an arc pressure of 0.70 MPa against ft = 10 MPa, and the
    // debris it freed sent the general contact (and the run time) through the
    // roof. Anything held by fewer than two joints is removed, iteratively
    // because removing a flap can orphan its neighbour.
    long dropped = 0;
    for (int pass = 0;; ++pass) {
        std::map<std::pair<int, int>, int> use;
        for (const auto& t : tris)
            for (int k = 0; k < 3; ++k) {
                auto key = std::minmax(t[k], t[(k + 1) % 3]);
                ++use[{key.first, key.second}];
            }
        std::vector<std::array<int, 3>> okT;
        std::vector<int> okG;
        for (std::size_t k = 0; k < tris.size(); ++k) {
            int shared = 0;
            for (int e = 0; e < 3; ++e) {
                auto key = std::minmax(tris[k][e], tris[k][(e + 1) % 3]);
                if (use[{key.first, key.second}] == 2) ++shared;
            }
            if (shared >= 2) { okT.push_back(tris[k]); okG.push_back(triGrain[k]); }
        }
        if (okT.size() == tris.size()) break;
        dropped += (long)(tris.size() - okT.size());
        tris.swap(okT);
        triGrain.swap(okG);
        if (tris.empty())
            throw std::runtime_error("geometry = disc: the rim cleanup removed "
                                     "every element — mesh far too coarse for "
                                     "the disc");
        if (pass > 50) break;                          // pathological mesh guard
    }
    // "held by 2 joints" does NOT imply "attached to the specimen": a cluster
    // of three mutually-joined triangles satisfies it and floats free. Keep the
    // LARGEST connected component only, otherwise those islands are counted as
    // fragments from frame 0 and pollute the breakage census.
    std::map<std::pair<int, int>, std::vector<int>> owner;
    for (std::size_t k = 0; k < tris.size(); ++k)
        for (int e = 0; e < 3; ++e) {
            auto key = std::minmax(tris[k][e], tris[k][(e + 1) % 3]);
            owner[{key.first, key.second}].push_back((int)k);
        }
    std::vector<int> comp(tris.size(), -1);
    int nComp = 0;
    for (std::size_t s = 0; s < tris.size(); ++s) {
        if (comp[s] >= 0) continue;
        std::queue<int> qu;
        qu.push((int)s);
        comp[s] = nComp;
        while (!qu.empty()) {
            int a = qu.front(); qu.pop();
            for (int e = 0; e < 3; ++e) {
                auto key = std::minmax(tris[a][e], tris[a][(e + 1) % 3]);
                for (int b : owner[{key.first, key.second}])
                    if (comp[b] < 0) { comp[b] = nComp; qu.push(b); }
            }
        }
        ++nComp;
    }
    if (nComp > 1) {
        std::vector<int> size(nComp, 0);
        for (int c : comp) ++size[c];
        int best = (int)(std::max_element(size.begin(), size.end())
                         - size.begin());
        std::vector<std::array<int, 3>> mainT;
        std::vector<int> mainG;
        for (std::size_t k = 0; k < tris.size(); ++k)
            if (comp[k] == best) {
                mainT.push_back(tris[k]);
                mainG.push_back(triGrain[k]);
            }
        dropped += (long)(tris.size() - mainT.size());
        tris.swap(mainT);
        triGrain.swap(mainG);
    }
    std::cout << "[FDEM] disc cut: " << tris.size() << " elements kept, "
              << dropped << " removed (flapping or detached islands; "
              << nComp << " components before cleanup)\n";
    // The rim is left as the STAIRCASE the cut produces. Projecting the
    // boundary vertices radially onto the circle was tried and dropped: it does
    // give a smooth disc, but stretching those rim triangles collapses their
    // inscribed size, and the joint penalty (20 E/h) is inversely proportional
    // to it — the stable time step fell by 8.6x on the 20x20 case for a
    // cosmetic gain. The brazilian load does not need a smooth rim anyway: it
    // is applied as a traction on the horizontal facets inside the load arc
    // (setupBrazilianLoad), and element-scale roughness of the rim does not
    // reach the centre of the disc, where the measurement is made.
}

// ---------------------------------------------------------------------------
// mesh = file — maillage triangulaire non structure importe (Gmsh MSH 2.2
// ASCII, elements type 2) : le maillage "a la Yan et al." en 2D. Meme
// contrat que le 3D : translation a l'origine, W/H relus de la boite
// englobante, orientation CCW reparee ici (buildFromTriangles la refuse),
// chemin de dimensionnement non uniforme (voronoi_ = true, un grain).
// ---------------------------------------------------------------------------
void FdemSolver::buildMeshFile() {
    std::string path = cfg_.reqs("meshFile");
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("meshFile: cannot open '" + path + "'");
    std::string line;
    std::map<long, int> id2idx;
    std::vector<Eigen::Vector2d> vpos;
    std::vector<std::array<int, 3>> tris;
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
        } else if (line.rfind("$Nodes", 0) == 0) {
            long n = 0;
            in >> n;
            for (long k = 0; k < n; ++k) {
                long id; double x, y, z;
                in >> id >> x >> y >> z;
                id2idx[id] = (int)vpos.size();
                vpos.push_back({x, y});
            }
        } else if (line.rfind("$Elements", 0) == 0) {
            long n = 0;
            in >> n;
            for (long k = 0; k < n; ++k) {
                long id; int type, ntags;
                in >> id >> type >> ntags;
                long tag;
                for (int t = 0; t < ntags; ++t) in >> tag;
                int nn = type == 15 ? 1 : type == 1 ? 2 : type == 2 ? 3
                       : type == 4 ? 4 : -1;
                if (nn < 0)
                    throw std::runtime_error("meshFile: element type "
                        + std::to_string(type) + " unsupported (2D: "
                        "triangles only)");
                if (nn == 4)
                    throw std::runtime_error("meshFile: tetrahedra found — "
                        "this is a 3D mesh; use mode = fdem3d");
                std::array<int, 3> vv{};
                for (int q = 0; q < nn; ++q) {
                    long nid; in >> nid;
                    if (nn == 3) {
                        auto it = id2idx.find(nid);
                        if (it == id2idx.end())
                            throw std::runtime_error("meshFile: element "
                                "references unknown node id");
                        vv[q] = it->second;
                    }
                }
                if (nn == 3) tris.push_back(vv);       // points/lines skipped
            }
        }
    }
    if (!sawFormat || vpos.empty() || tris.empty())
        throw std::runtime_error("meshFile: no triangles found in '" + path
                                 + "' (need ASCII MSH 2.2 with type-2 "
                                 "elements)");
    Eigen::Vector2d lo = vpos[0], hi = vpos[0];
    for (const auto& p : vpos) { lo = lo.cwiseMin(p); hi = hi.cwiseMax(p); }
    for (auto& p : vpos) p -= lo;
    W_ = hi.x() - lo.x();
    H_ = hi.y() - lo.y();
    if (!(W_ > 0 && H_ > 0))
        throw std::runtime_error("meshFile: degenerate bounding box");
    for (auto& t : tris) {                             // force CCW
        const auto &A = vpos[t[0]], &B = vpos[t[1]], &C = vpos[t[2]];
        double det = (B.x() - A.x()) * (C.y() - A.y())
                   - (C.x() - A.x()) * (B.y() - A.y());
        if (det < 0) std::swap(t[1], t[2]);
    }
    std::vector<int> triGrain(tris.size(), 0);
    nGrains_ = 1;
    voronoi_ = true;             // non-uniform sizing paths — not a grid
    std::cout << "[FDEM] mesh = file: '" << path << "' — " << vpos.size()
              << " nodes, " << tris.size() << " triangles, box " << W_
              << " x " << H_ << " m\n";
    buildFromTriangles(vpos, tris, triGrain, {0});
    hmin_ = 1e30;
    for (double h : hEl_) hmin_ = std::min(hmin_, h);
}

void FdemSolver::buildMeshVoronoi() {
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

    // grainMesh = fan (default, unchanged) | delaunay. The delaunay front-end
    // is what the grain-based FDEM literature does: an unstructured mesh INSIDE
    // every grain at about 0.18 x the grain diameter (Y-Geo / Irazu GBM), so
    // that transgranular cracking is not confined to the spokes of a centroid
    // fan. It is opt-in so every earlier voronoi result stays reproducible.
    std::string gm = cfg_.gets("grainMesh", "fan");
    if (gm != "fan" && gm != "delaunay")
        throw std::runtime_error("grainMesh must be fan | delaunay (got '"
                                 + gm + "')");
    double gh = cfg_.getd("grainElemSize", 0.0);
    Tessellation T = Tessellation::build(W_, H_, d, jit, lloyd, mf, refine,
                                         phases_.fraction, rng,
                                         seeding == "random",
                                         gm == "delaunay", gh);
    nGrains_ = T.nGrains;

    std::vector<std::array<int, 3>> tris;
    std::vector<int> triGrain;
    tris.reserve(T.tri.size());
    triGrain.reserve(T.tri.size());
    for (const auto& t : T.tri) {
        tris.push_back(t.v);
        triGrain.push_back(t.grain);
    }
    cutDisc(T.vtx, tris, triGrain);
    if (disc_) {                                       // refresh the census
        std::vector<char> seen(nGrains_, 0);
        for (int g : triGrain) if (g >= 0 && g < nGrains_) seen[g] = 1;
        nGrains_ = (int)std::count(seen.begin(), seen.end(), (char)1);
    }
    buildFromTriangles(T.vtx, tris, triGrain, T.phaseOfGrain);

    // the length scale of the contact machinery and the CFL is now set by
    // the smallest element the tessellation produced (inscribed size 4A/P)
    hmin_ = 1e30;
    for (double h : hEl_) hmin_ = std::min(hmin_, h);
}

// ---------------------------------------------------------------------------
// Shared topology builder: per-element node triplets, cohesive joints on the
// doubly-shared virtual edges (with the CCW/outward orientation fixes that
// the energy-pump hunt made mandatory), exterior faces, masses, tension
// grips. Virtual vertex ids come from the caller: two triangles are joined
// iff they share two virtual ids, inside grains and across them alike.
// ---------------------------------------------------------------------------
void FdemSolver::buildFromTriangles(const std::vector<Eigen::Vector2d>& vpos,
                                    const std::vector<std::array<int, 3>>& tris,
                                    const std::vector<int>& triGrain,
                                    const std::vector<int>& grainPhase) {
    std::map<std::pair<int, int>, std::vector<std::array<int, 3>>> edges;
    el_.reserve(tris.size());
    hEl_.reserve(tris.size());

    for (std::size_t tId = 0; tId < tris.size(); ++tId) {
        int va = tris[tId][0], vb = tris[tId][1], vs = tris[tId][2];
        Elem e;
        int base = (int)X0_.size();
        for (int v : {va, vb, vs}) {
            X0_.push_back(vpos[v]);
            elemOf_.push_back((int)el_.size());
            vOf_.push_back(v);
        }
        e.n = {base, base + 1, base + 2};
        const auto& A = X0_[e.n[0]];
        const auto& B = X0_[e.n[1]];
        const auto& C = X0_[e.n[2]];
        double det = (B.x() - A.x()) * (C.y() - A.y())
                   - (C.x() - A.x()) * (B.y() - A.y());
        e.A0 = 0.5 * det;
        if (e.A0 <= 0) throw std::runtime_error("inverted element in mesh gen");
        // dN(:,a) = grad N_a in the reference configuration
        e.dN.col(0) = Eigen::Vector2d(B.y() - C.y(), C.x() - B.x()) / det;
        e.dN.col(1) = Eigen::Vector2d(C.y() - A.y(), A.x() - C.x()) / det;
        e.dN.col(2) = Eigen::Vector2d(A.y() - B.y(), B.x() - A.x()) / det;
        e.grain = triGrain.empty() ? 0 : triGrain[tId];
        e.phase = grainPhase.empty() ? 0 : grainPhase[e.grain];
        double per = (B - A).norm() + (C - B).norm() + (A - C).norm();
        hEl_.push_back(4.0 * e.A0 / per);              // inscribed diameter
        int id = (int)el_.size();
        el_.push_back(e);
        int vv[4] = {va, vb, vs, va};
        int nn[4] = {e.n[0], e.n[1], e.n[2], e.n[0]};
        for (int k = 0; k < 3; ++k) {
            auto key = std::minmax(vv[k], vv[k + 1]);
            edges[{key.first, key.second}].push_back(
                {id, /*P*/ vv[k] == key.first ? nn[k] : nn[k + 1],
                     /*Q*/ vv[k] == key.first ? nn[k + 1] : nn[k]});
        }
    }

    // joints on doubly-shared edges, exterior list on the rest
    for (auto& [key, lst] : edges) {
        if (lst.size() == 2) {
            Joint J;
            J.eA = lst[0][0];
            J.eB = lst[1][0];
            // recover A's CCW traversal (P -> Q as stored may be either
            // orientation; the stored pair is keyed on the sorted virtual ids,
            // so find which of A's stored nodes is first in its own CCW walk)
            J.a1 = lst[0][1]; J.a2 = lst[0][2];
            J.b1 = lst[1][1]; J.b2 = lst[1][2];
            // ensure (a1 -> a2) is CCW in element A: the stored order is the
            // sorted-virtual order; flip if needed using the outward test
            const Elem& EA = el_[J.eA];
            Eigen::Vector2d cenA = (X0_[EA.n[0]] + X0_[EA.n[1]] + X0_[EA.n[2]]) / 3.0;
            Eigen::Vector2d P = X0_[J.a1], Q = X0_[J.a2];
            Eigen::Vector2d nrm(( Q - P ).y(), -( Q - P ).x());
            if (nrm.dot(0.5 * (P + Q) - cenA) < 0) {   // normal must leave A
                std::swap(J.a1, J.a2);
                std::swap(J.b1, J.b2);
            }
            J.L0 = (X0_[J.a2] - X0_[J.a1]).norm();
            jt_.push_back(J);
        } else if (lst.size() == 1) {
            // Same orientation fix as for joints: the stored pair is in
            // sorted-virtual-id order, but the outward-normal formula
            // n = ((Q-P).y, -(P-Q).x)... assumes (P -> Q) is the CCW
            // traversal of the owning element. Half the exterior edges were
            // flipped, turning them into permanent phantom pushers for any
            // node on their inside — the energy pump found by bisection.
            BEdge be{lst[0][0], lst[0][1], lst[0][2]};
            const Elem& E = el_[be.elem];
            Eigen::Vector2d cen = (X0_[E.n[0]] + X0_[E.n[1]] + X0_[E.n[2]]) / 3.0;
            Eigen::Vector2d P = X0_[be.na], Q = X0_[be.nb];
            Eigen::Vector2d nrm((Q - P).y(), -(Q - P).x());
            if (nrm.dot(0.5 * (P + Q) - cen) < 0) std::swap(be.na, be.nb);
            exterior_.push_back(be);
        } else {
            throw std::runtime_error("mesh topology error: an edge is shared "
                                     "by more than two triangles (vertex weld "
                                     "produced a non-manifold mesh)");
        }
    }

    u_.assign(X0_.size(), Eigen::Vector2d::Zero());
    v_.assign(X0_.size(), Eigen::Vector2d::Zero());
    f_.assign(X0_.size(), Eigen::Vector2d::Zero());
    m_.assign(X0_.size(), 0.0);
    flag_.assign(X0_.size(), FREE);
    for (const auto& e : el_)
        for (int a = 0; a < 3; ++a)
            m_[e.n[a]] += phases_.mat[e.phase].rho * e.A0 * thk_ / 3.0;

    // Clamped grip rows ONLY when the test is grip-driven. With platens the
    // specimen is held by contact alone — no zero-thickness clamped layer, so
    // no boundary layer for a refined mesh to resolve and fail in.
    if (scen_ == Scenario::TENSION && !tensionPlatens_) {
        for (int i = 0; i < (int)X0_.size(); ++i) {
            if (X0_[i].y() < 1e-9)      flag_[i] = FIXED;
            else if (X0_[i].y() > H_ - 1e-9) flag_[i] = PRESCRIBED;
        }
    }

    // Lateral rollers on the two flanks of the box: u_x = 0, u_y free.
    // Off by default. setupBoundaries() runs later and may promote the
    // bottom row to FIXED, which is what the confined-strip test wants.
    if (cfg_.getb("lateralRollers", false)) {
        if (disc_)
            throw std::runtime_error("lateralRollers needs geometry = box "
                                     "(a disc has no flanks)");
        for (int i = 0; i < (int)X0_.size(); ++i)
            if (flag_[i] == FREE
                && (X0_[i].x() < 1e-9 || X0_[i].x() > W_ - 1e-9))
                flag_[i] = ROLLERX;
    }
}

// ---------------------------------------------------------------------------
// Per-joint cohesive properties. Intra-grain joints (and every joint of the
// grid mesh) carry the bulk material of their phase. Grain-boundary joints
// take the MEAN of the two neighbouring phases times the alpha attenuation
// factors; heterophase boundaries get the extra heteroFactor on the
// strength-like properties — the classification the GBM literature uses.
// The penalty uses the local element size for the voronoi mesh (elements are
// not uniform) and the global hmin for the grid mesh (bit-compatible with
// the pre-GBM behaviour).
// ---------------------------------------------------------------------------
void FdemSolver::assignJointProps() {
    // Adaptive insertion: the penalty never glues an intact continuum (bonded
    // edges are handled kinematically), it only serves the ACTIVATED joints as
    // unloading/contact stiffness. It can therefore be much softer than the
    // intrinsic factor without softening the specimen — that is the whole
    // point of the article (its fig. 9-10: accuracy AND dt). Default 4 E/h.
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
        // defense in depth: PhaseSet::validate already rejects zero-strength
        // materials, but the attenuated MEANS must stay physical too (ft = 0
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
        // Critical opening / slip. The softening branch must enclose exactly
        // GfI (resp. GfII), the elastic part excluded — that is what fixes
        // the factor: linear branch of peak ft and width w has area ft w / 2,
        // the f(D) branch of eq. 11 has area ft w I with I = int_0^1 f(D) dD
        // (eq. 13 and 15 of the article).
        double kI = yanSoft_ ? 1.0 / yanI_ : 2.0;
        J.dnF = J.dnE + kI * Gf / ft;                  // mode-I critical opening
        J.slipF = kI * GfII / coh;                     // mode-II critical slip
        J.tanPhi = std::tan(phiDeg * M_PI / 180.0);
    }
    applyJointStatistics();
}

// ---------------------------------------------------------------------------
// Statistical joint strengths — the implicit-DFH idea transplanted onto the
// joint network. With jointWeibullM = m > 1, every joint's ft and cohesion
// are multiplied by a Weibull(m) factor of MEAN 1 (scale 1/Gamma(1+1/m)),
// so the calibrated deterministic strengths stay the ensemble mean.
// Fracture energies are NOT scaled (defect statistics affect the strength,
// the toughness stays a material property); the openings dnF/slipF are
// recomputed accordingly.
//
// strengthCorrLength selects the spatial structure:
//   = 0 : independent draw per joint (the analogue of the per-element
//         Weibull draw of smeared DFH implementations — statistics converge
//         with the mesh, the crack MAP does not);
//   > 0 : the factors sample ONE Gaussian random field (RandomField) with
//         that correlation length through the Gaussian copula
//         u = Phi(g(x_mid)) — the field lives in SPACE, independent of the
//         mesh, so two different meshes see the same weak zones and the
//         crack map becomes reproducible. fieldSeed (default seed + 777)
//         controls the field independently of the mesh seed.
// ---------------------------------------------------------------------------
void FdemSolver::applyJointStatistics() {
    const bool wGf =
        cfg_.gets("weibullScope", "strength") == "strengthGf";
    double m = cfg_.getd("jointWeibullM", 0.0);
    // Effet d'echelle statistique de Weibull — voir le commentaire detaille de
    // Fdem3dSolver::applyJointSizeEffect. Meme convention que les VUMAT DP-DFH
    // (sig_k = sigw*(Zeff/V_el)^(1/m)). OPT-IN : sans jointSizeEffect, chemin
    // inchange et resultats identiques au bit pres.
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
        // anisotropic option: a second length across the first (bands) and
        // an orientation — the foliation-like texture knob
        double ellB = cfg_.getd("strengthCorrLengthB", ell);
        double ang = cfg_.getd("strengthCorrAngleDeg", 0.0);
        RandomField F(W_, H_, ell, ellB, ang, fseed);
        for (auto& J : jt_) {
            Eigen::Vector2d mid = 0.5 * (X0_[J.a1] + X0_[J.a2]);
            double g = F(mid);
            double u = 0.5 * std::erfc(-g / std::sqrt(2.0));   // Phi(g)
            J.stat = weib(u);
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
        // ---- weibullScope (2026-08-19) : la TENACITE suit-elle la
        // resistance ? Voir MatLaw.hpp pour les deux conventions. Defaut
        // `strength` = comportement historique, bit-identique.
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
        std::cout << "[FDEM] joint strength statistics: Weibull m = " << m
                  << (ell > 0.0 ? " correlated, ell = " + std::to_string(ell)
                                : std::string(" independent per joint"))
                  << ", factor mean/min/max = " << xsum / jt_.size() << "/"
                  << xmin << "/" << xmax << "\n";
    if (szOn)
        std::cout << "[FDEM] TOTAL ft factor (Weibull x taille) "
                     "mean/min/max = " << xsum / jt_.size() << "/"
                  << xmin << "/" << xmax << "\n";
}

// ---------------------------------------------------------------------------
// Effet d'echelle statistique 2D : ft <- ft * (Zeff/V_J)^(1/m).
// ATTENTION a l'epaisseur : en 2D le volume represente vaut aire x thickness,
// et `thickness` vaut souvent 1 m (valeur conventionnelle, pas physique). Zeff
// doit donc etre declare de façon COHERENTE avec l'epaisseur employee — sinon
// le facteur est correct en tendance mais faux en niveau. Le log imprime V_J
// pour permettre de le verifier d'un coup d'oeil.
// ---------------------------------------------------------------------------
void FdemSolver::applyJointSizeEffect(double mWeib) {
    double mS = cfg_.getd("jointSizeEffectM", mWeib);
    if (!(mS > 1.0))
        throw std::runtime_error("jointSizeEffect: exposant invalide — poser "
            "jointSizeEffectM > 1 (ou jointWeibullM)");
    double Zeff = cfg_.getd("jointZeff", 1e-9);        // 1 mm^3, defaut VUMAT
    if (!(Zeff > 0.0))
        throw std::runtime_error("jointZeff doit etre > 0 [m^3]");
    double cap = cfg_.getd("jointSizeEffectClamp", 5.0);
    if (!(cap >= 1.0))
        throw std::runtime_error("jointSizeEffectClamp doit etre >= 1");
    double fmin = 1e300, fmax = 0.0, fsum = 0.0, vmin = 1e300, vmax = 0.0;
    std::size_t nclip = 0;
    for (auto& J : jt_) {
        double Vj = 0.5 * (el_[J.eA].A0 + el_[J.eB].A0) * thk_;
        double f = 1.0;
        if (Vj > 0.0) f = std::pow(Zeff / Vj, 1.0 / mS);
        if (f > cap)            { f = cap;       ++nclip; }
        else if (f < 1.0 / cap) { f = 1.0 / cap; ++nclip; }
        J.stat *= f;
        fmin = std::min(fmin, f); fmax = std::max(fmax, f); fsum += f;
        vmin = std::min(vmin, Vj); vmax = std::max(vmax, Vj);
    }
    std::cout << "[FDEM] effet d'echelle (Zeff/V_J)^(1/m) : Zeff = " << Zeff
              << " m^3, m = " << mS << ", V_J = aire x thickness ("
              << thk_ << " m) min/max = " << vmin << "/" << vmax
              << " m^3, facteur mean/min/max = " << fsum / jt_.size()
              << "/" << fmin << "/" << fmax << "\n";
    if (nclip)
        std::cout << "[FDEM] WARNING: " << nclip << " joints bornes a " << cap
                  << "x (jointSizeEffectClamp) — maillage tres heterogene, "
                     "Zeff mal choisi, ou thickness non physique\n";
}

// ===========================================================================
// Adaptive insertion (Yan, Zheng & Wang, IJRMMS 169, 2023, 105439)
//
// The article starts from a SHARED-NODE mesh and splits nodes when a cohesive
// element is inserted (its fig. 7). This solver's data layout is the exact
// dual: nodes are ALREADY duplicated per element, so "shared" is enforced
// kinematically — the co-located copies of an original vertex are bound into
// groups that integrate as one node (sum of forces, sum of masses, common
// velocity). Binding groups are the connected components of the element fan
// around the vertex, two elements being connected when the edge between them
// is still BONDED. Activating a joint and re-running the union-find at its
// two endpoint vertices reproduces the progressive splitting of fig. 7: a
// crack-tip vertex stays whole (the fan is still connected around the tip),
// a traversed vertex splits into exactly the fan components.
// ===========================================================================
void FdemSolver::buildBindingTables() {
    nVert_ = 0;
    for (int v : vOf_) nVert_ = std::max(nVert_, v + 1);
    copiesOfVert_.assign(nVert_, {});
    for (int i = 0; i < (int)X0_.size(); ++i)
        copiesOfVert_[vOf_[i]].push_back(i);
    jointsOfVert_.assign(nVert_, {});
    for (int jI = 0; jI < (int)jt_.size(); ++jI) {
        jointsOfVert_[vOf_[jt_[jI].a1]].push_back(jI);
        jointsOfVert_[vOf_[jt_[jI].a2]].push_back(jI);
    }
    grpsOfVert_.assign(nVert_, {});
    for (int v = 0; v < nVert_; ++v) rebindVertex(v);
}

void FdemSolver::rebindVertex(int v) {
    const auto& copies = copiesOfVert_[v];
    auto& grps = grpsOfVert_[v];
    grps.clear();
    if (copies.empty()) return;
    // local union-find over the elements of the fan
    std::vector<int> elems;
    elems.reserve(copies.size());
    for (int i : copies) elems.push_back(elemOf_[i]);
    auto local = [&](int e) {
        for (int k = 0; k < (int)elems.size(); ++k)
            if (elems[k] == e) return k;
        return -1;                                     // not in this fan
    };
    std::vector<int> par(elems.size());
    for (int k = 0; k < (int)par.size(); ++k) par[k] = k;
    std::function<int(int)> find = [&](int x) {
        while (par[x] != x) { par[x] = par[par[x]]; x = par[x]; }
        return x;
    };
    for (int jI : jointsOfVert_[v]) {
        const Joint& J = jt_[jI];
        if (!J.bonded) continue;                       // inserted: edge is cut
        int a = local(J.eA), b = local(J.eB);
        if (a < 0 || b < 0) continue;
        par[find(a)] = find(b);
    }
    // groups of copies keyed by the component of their element
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

// The insertion criterion of the article, eq. 7-8: the traction on each
// bonded edge is the average of the two neighbouring element stress tensors
// projected on the edge frame; the joint is activated when sigma_n >= ft or
// |tau| >= fs, with fs = c - sigma_n tan(phi) in compression, c otherwise.
// The sweep is O(bonded edges) per step and runs before jointForces so a
// newborn joint carries traction the very step it is inserted.
void FdemSolver::insertionSweep() {
    struct Hit { int jI; double sig, tau; };
    std::vector<Hit> hits;
#ifdef _OPENMP
    #pragma omp parallel
    {
        std::vector<Hit> mine;
        #pragma omp for schedule(static) nowait
        for (int jI = 0; jI < (int)jt_.size(); ++jI) {
            const Joint& J = jt_[jI];
            if (!J.bonded) continue;
            Eigen::Vector2d P = 0.5 * (X0_[J.a1] + u_[J.a1] + X0_[J.b1] + u_[J.b1]);
            Eigen::Vector2d Q = 0.5 * (X0_[J.a2] + u_[J.a2] + X0_[J.b2] + u_[J.b2]);
            Eigen::Vector2d ed = Q - P;
            double L = ed.norm();
            if (L < 1e-14) continue;
            Eigen::Vector2d e = ed / L;
            Eigen::Vector2d n(e.y(), -e.x());
            const Elem& A = el_[J.eA];
            const Elem& B = el_[J.eB];
            double sxx = 0.5 * (A.sxx + B.sxx);
            double syy = 0.5 * (A.syy + B.syy);
            double sxy = 0.5 * (A.sxy + B.sxy);
            double sig = n.x() * (sxx * n.x() + sxy * n.y())
                       + n.y() * (sxy * n.x() + syy * n.y());
            double tau = e.x() * (sxx * n.x() + sxy * n.y())
                       + e.y() * (sxy * n.x() + syy * n.y());
            // DIF de Yang et al. : le critere doit voir la resistance
            // DYNAMIQUE, sinon le joint s insere au seuil statique et le
            // facteur applique ensuite serait sans effet sur l instant
            // d insertion. Le terme de frottement -sig tan(phi) n est PAS
            // amplifie (leur choix, source Zhao).
            double dT = 1.0, dC = 1.0;
            if (difOn_) {
                double er = 0.5 * (A.edot + B.edot);
                dT = difTensionYang(er, difExpT_);
                dC = difCompressionYang(er);
            }
            double fs = dC * J.coh
                      + J.tanPhi * rockim::mcFrictionTerm(sig, J.ft, yangEnv_);
            if (fs < 0.0) fs = 0.0;
            if (sig >= dT * J.ft || std::abs(tau) >= fs)
                mine.push_back({jI, sig, tau});
        }
        #pragma omp critical
        hits.insert(hits.end(), mine.begin(), mine.end());
    }
#else
    for (int jI = 0; jI < (int)jt_.size(); ++jI) {
        const Joint& J = jt_[jI];
        if (!J.bonded) continue;
        Eigen::Vector2d P = 0.5 * (X0_[J.a1] + u_[J.a1] + X0_[J.b1] + u_[J.b1]);
        Eigen::Vector2d Q = 0.5 * (X0_[J.a2] + u_[J.a2] + X0_[J.b2] + u_[J.b2]);
        Eigen::Vector2d ed = Q - P;
        double L = ed.norm();
        if (L < 1e-14) continue;
        Eigen::Vector2d e = ed / L;
        Eigen::Vector2d n(e.y(), -e.x());
        const Elem& A = el_[J.eA];
        const Elem& B = el_[J.eB];
        double sxx = 0.5 * (A.sxx + B.sxx);
        double syy = 0.5 * (A.syy + B.syy);
        double sxy = 0.5 * (A.sxy + B.sxy);
        double sig = n.x() * (sxx * n.x() + sxy * n.y())
                   + n.y() * (sxy * n.x() + syy * n.y());
        double tau = e.x() * (sxx * n.x() + sxy * n.y())
                   + e.y() * (sxy * n.x() + syy * n.y());
        double dT = 1.0, dC = 1.0;             // DIF de Yang (voir plus haut)
        if (difOn_) {
            double er = 0.5 * (A.edot + B.edot);
            dT = difTensionYang(er, difExpT_);
            dC = difCompressionYang(er);
        }
        double fs = dC * J.coh
                  + J.tanPhi * rockim::mcFrictionTerm(sig, J.ft, yangEnv_);
        if (fs < 0.0) fs = 0.0;
        if (sig >= dT * J.ft || std::abs(tau) >= fs)
            hits.push_back({jI, sig, tau});
    }
#endif
    if (hits.empty()) return;
    // deterministic activation order whatever the thread count
    std::sort(hits.begin(), hits.end(),
              [](const Hit& x, const Hit& y) { return x.jI < y.jI; });
    for (const Hit& h : hits) activateJoint(h.jI, h.sig, h.tau);
}

// Stress continuity at insertion (the article's guard against the classical
// "time discontinuity" of extrinsic CZM, its section 2.4-2.5, transposed to
// this joint law): the newborn joint must transmit at zero geometric opening
// exactly the traction the continuum was carrying.
//   * normal — opening offset dn0 = min(sig, ft)/pj: the effective opening
//     dn0 puts the elastic branch at sig (tension-triggered edges, where
//     sig >= ft up to the one-step overshoot, start exactly at the envelope
//     peak; shear-triggered edges keep their sub-critical or compressive
//     normal traction, the article's case 2);
//   * shear — plastic-slip offset so the trial traction pj*(dtg - slip)
//     equals the transmitted tau (clamped to the current Coulomb cap, the
//     article's f(D)*fs) at dtg = 0.
// Then the union-find is re-run at the two endpoint vertices: fig. 7 for
// free, including the third re-split of its node 2.
void FdemSolver::activateJoint(int jI, double sig, double tau) {
    Joint& J = jt_[jI];
    if (!J.bonded) return;
    J.bonded = false;
    // ---- DIF de Yang et al. 2025, FIGE ICI ------------------------------
    // Applique AVANT les decalages de continuite de contrainte ci-dessous,
    // qui lisent J.ft, J.coh et J.pj. Se COMPOSE avec le facteur de Weibull
    // J.stat deja porte par ft et coh (multiplicatif, il ne le remplace pas).
    // L angle de frottement est inchange. Les ouvertures critiques sont
    // recalculees : dnE = ft/pj suit le DIF, tandis que kI Gf/ft ne bouge pas
    // puisque ft et Gf recoivent le meme facteur — dnF est donc quasi
    // invariante et le compteur d endommagement reste coherent.
    if (difOn_) {
        double er = 0.5 * (el_[J.eA].edot + el_[J.eB].edot);
        double dT = difTensionYang(er, difExpT_);
        double dC = difCompressionYang(er);
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
    double fsNow = J.coh
                 + J.tanPhi * rockim::mcFrictionTerm(sig, J.ft, yangEnv_);
    if (fsNow < 0.0) fsNow = 0.0;
    double tau0 = std::clamp(tau, -fsNow, fsNow);
    J.slip[0] = J.slip[1] = -tau0 / J.pj;
    ++nInserted_;
    rebindVertex(vOf_[J.a1]);
    rebindVertex(vOf_[J.a2]);
}

void FdemSolver::placeTool() {
    if (scen_ == Scenario::TENSION || scen_ == Scenario::BRAZILIAN
        || scen_ == Scenario::SHPB) return;        // no rigid tool
    tool_.mass = cfg_.getd("toolMass", 5.0);
    double gap = cfg_.getd("toolGap", 1e-4);
    std::string sh = cfg_.gets("toolShape", "disc");
    tool_.shape = (sh == "flat") ? Tool::Shape::FLAT
                : (sh == "pdc")  ? Tool::Shape::PDC : Tool::Shape::DISC;
    if (sh != "flat" && sh != "pdc" && sh != "disc")
        throw std::runtime_error("toolShape must be disc | flat | pdc (2D)");
    tool_.width  = cfg_.getd("toolWidth", 0.02);
    tool_.radius = cfg_.getd("toolRadius", 0.015);

    if (scen_ == Scenario::PERCUSSION) {
        tool_.motion = Tool::Motion::FREE;
        double vImp = cfg_.getd("impactSpeed", 8.0);
        double xc = cfg_.getd("toolX", 0.5 * W_);
        // ---- E2 (2026-08-19) : toolShape est desormais HONORE en percussion.
        // Avant ce correctif la valeur lue etait ECRASEE par DISC trois lignes
        // plus bas : un run « poincon plat » et un run « disque » de meme rayon
        // etaient bit-identiques, et c'est ainsi que le defaut a ete decouvert
        // (18/08). Le 3D honorait deja la cle. Audit prealable du 19/08 : les
        // DIX configs 2D-percussion du depot posent toolShape = disc, donc ce
        // correctif ne change aucun resultat existant — il rend seulement le
        // poincon plat 2D atteignable.
        if (sh == "pdc")
            throw std::runtime_error(
                "toolShape = pdc n'a pas de sens en scenario = percussion : "
                "un cutter PDC est un outil de COUPE (scenario = shear). "
                "Utiliser disc (bouton) ou flat (poincon plat).");
        // Placement selon la geometrie : `x` est le CENTRE pour un disque, le
        // milieu de la FACE INFERIEURE pour un poincon plat (cf. Tool.hpp).
        double yTop = (tool_.shape == Tool::Shape::FLAT) ? H_ + gap
                                                         : H_ + tool_.radius + gap;
        tool_.x = {xc, yTop};
        tool_.v = {0.0, -vImp};
        std::cout << "[FDEM] outil percussif : "
                  << (tool_.shape == Tool::Shape::FLAT
                          ? "poincon plat, largeur " : "disque, rayon ")
                  << (tool_.shape == Tool::Shape::FLAT ? tool_.width
                                                       : tool_.radius)
                  << " m, masse " << tool_.mass << " kg/m, vitesse d'impact "
                  << vImp << " m/s, jeu initial " << gap << " m (contact a t = "
                  << gap / vImp << " s)\n";
        // REPARATION (2026-08-28, decision F. Uzquiano) : la cle n etait lue
        // que dans la branche PDC du scenario shear — en percussion elle
        // etait ignoree EN SILENCE (imp2d_panoplie posait 1.0 en croyant
        // l armer : son controle ne controlait rien). L ecretage lui-meme
        // (site nodeFc) existait deja ; seul le READ manquait.
        toolVCap_ = cfg_.getd("toolImpulseCap", 0.0);
        if (toolVCap_ > 0.0)
            std::cout << "[FDEM] toolImpulseCap = " << toolVCap_
                      << " : |Fc| <= kappa * 2 v_outil * m / dt par noeud "
                         "(percussion — actif depuis le 2026-08-28)\n";
        // (5) meme visibilite du piege jointDeath que le 3D
        // (mpka9c : jointDeath n existe pas sur cette branche — seule la
        // semantique separation existe, la notice est inconditionnelle)
        std::cout << "[FDEM] jointDeath = separation (defaut) : sous "
                         "l indenteur un joint ecroui en compression ne meurt "
                         "jamais, le relais contact roche/roche ne s engage "
                         "pas. Levier : jointDeath = damage.\n";
    } else {                                           // SHEAR: lateral cut
        tool_.motion = Tool::Motion::PRESCRIBED;
        double depth = cfg_.getd("cutDepth", 0.004);
        double vCut  = cfg_.getd("cutSpeed", 10.0);
        if (sh == "pdc") {
            // PDC cutter: `x` IS the cutting edge, placed at the depth of cut
            // and started clear of the specimen so it engages progressively.
            tool_.shape = Tool::Shape::PDC;
            tool_.rakeDeg = cfg_.getd("backRakeDeg", 20.0);
            tool_.faceLen = cfg_.getd("cutterLen", 0.013);
            tool_.chamLen = cfg_.getd("chamferLen", 0.0);
            tool_.chamDeg = cfg_.getd("chamferDeg", 45.0);
            // ---- E3 (2026-08-19) : le chanfrein est LU, STOCKE, et n'entre
            // dans AUCUN test de contact du cutter. Dix configs de la campagne
            // de coupe le posent, en croyant modeler une arete chanfreinee.
            // Le rendre operant changerait leurs resultats ; la constitution
            // (I) l'interdit sans decision explicite. En attendant, il ne sera
            // plus SILENCIEUX : un reglage qui ne fait rien est pire que son
            // absence, et il doit au moins le dire.
            if (tool_.chamLen > 0.0)
                std::cout << "\n[FDEM] *** AVERTISSEMENT *** chamferLen = "
                          << tool_.chamLen << " m est lu mais **NON "
                             "IMPLEMENTE** : le chanfrein n'intervient dans "
                             "aucun test de contact du cutter. La geometrie "
                             "reellement simulee est un coin SANS chanfrein. "
                             "Retirer la cle, ou implementer FR-008 (spec "
                             "003-cutter-pdc-3d).\n\n";
            tool_.thick   = cfg_.getd("cutterThick", 0.0);
            // A : ecretage en impulsion — voir FdemSolver.hpp. Lu ICI plutot
            // qu'avec les autres cles de contact pour rester a cote de la
            // geometrie d'outil qu'il borne.
            toolVCap_ = cfg_.getd("toolImpulseCap", 0.0);
            // ---- E6 (2026-08-19) : la trace MENTAIT. Elle lisait
            // tool_.v.norm() alors que tool_.v n'est assigne qu'a la fin de
            // cette fonction : le journal affichait invariablement « 0 m/s »
            // pendant que le plafond agissait bel et bien. Une trace qui ment
            // est exactement le piege rencontre avec toolShape. On lit donc
            // vCut, qui EST la vitesse imposee de l'outil.
            if (toolVCap_ > 0.0)
                std::cout << "[FDEM] toolImpulseCap = " << toolVCap_
                          << " : |Fc| <= " << toolVCap_
                          << " * 2 v_outil * m / dt — l'outil ne peut pas "
                             "changer la vitesse d'un noeud de plus de "
                          << toolVCap_ * 2.0 * std::abs(vCut) << " m/s par pas\n";
            if (!(tool_.rakeDeg > -60.0 && tool_.rakeDeg < 60.0))
                throw std::runtime_error("backRakeDeg must be in (-60, 60)");
            tool_.x = {cfg_.getd("toolX", -0.002), H_ - depth};
            std::cout << "[FDEM] PDC cutter: thickness "
                      << (tool_.thick > 0.0 ? tool_.thick : -1.0)
                      << " m (negative = unbounded wedge), back rake "
                      << tool_.rakeDeg
                      << " deg, face " << tool_.faceLen << " m, depth of cut "
                      << depth << " m, speed " << vCut << " m/s\n";
        } else {
            tool_.shape = Tool::Shape::DISC;
            tool_.x = {cfg_.getd("toolX", -tool_.radius - gap),
                       H_ - depth + tool_.radius};
        }
        tool_.v = {vCut, 0.0};
    }
}

// ---------------------------------------------------------------------------
// Encastrement / quiet boundaries on the exterior faces, mirroring the FEM
// module: absorbing = none | sides | all; 'all' replaces the fixed bottom by
// the viscous-spring (Lysmer + Deeks-Randolph) support.
// ---------------------------------------------------------------------------
void FdemSolver::setupBoundaries() {
    cAbsX_.assign(X0_.size(), 0.0);
    cAbsY_.assign(X0_.size(), 0.0);
    kAbsX_.assign(X0_.size(), 0.0);
    kAbsY_.assign(X0_.size(), 0.0);
    // ---- SHPB: driven face + viscous end (eq. 21-22 of the article) --------
    if (scen_ == Scenario::SHPB) {
        // (a) struck face x = 0 of the incident bar: v_x = shpbVel(t), v_y free
        long nDrive = 0;
        for (int i = 0; i < (int)X0_.size(); ++i)
            if (X0_[i].x() < 1e-9) { flag_[i] = DRIVEX; ++nDrive; }
        // (b) viscous boundary at the far end of the LAST bar. rockim's Lysmer
        //     dashpot is the classical, impedance-matched rho c v per unit
        //     area; the article's eq. 21 prints sigma = -2 rho cp vn, i.e.
        //     TWICE that. Twice the matched impedance is NOT transparent (the
        //     reflection coefficient of a dashpot z_d against a bar of
        //     impedance z is (z_d - z)/(z_d + z), so z_d = 2 z reflects +1/3 of
        //     the incident amplitude, sign-inverted); the standard
        //     Lysmer-Kuhlemeyer factor is 1. absorbFactor exposes the choice:
        //     default 1 (matched, what a "viscous boundary that absorbs
        //     outgoing waves" must do), absorbFactor = 2 reproduces eq. 21 as
        //     printed. The two are compared in the report.
        long nAbs = 0;
        for (const auto& be : exterior_) {
            const Material& mp = phases_.mat[el_[be.elem].phase];
            double zP = absFac_ * mp.rho * mp.cP() * thk_;
            double zS = absFac_ * mp.rho * mp.cS() * thk_;
            Eigen::Vector2d Pp = X0_[be.na], Qq = X0_[be.nb];
            if (!(Pp.x() > xEndAbs_ - 1e-9 && Qq.x() > xEndAbs_ - 1e-9))
                continue;
            double L2 = 0.5 * (Qq - Pp).norm();
            for (int nid : {be.na, be.nb}) {
                cAbsX_[nid] += zP * L2;                // normal = x here
                cAbsY_[nid] += zS * L2;
                ++nAbs;
            }
        }
        if (nDrive == 0 || nAbs == 0)
            throw std::runtime_error("scenario = shpb: the driven face or the "
                                     "viscous end caught no node");
        const Material& mb = phases_.mat[el_[0].phase];
        std::cout << "[FDEM] shpb: " << nDrive << " driven nodes at x = 0, "
                     "viscous end at x = " << xEndAbs_ << " m (absorbFactor = "
                  << absFac_ << " x rho c, cP = " << mb.cP() << " m/s, cS = "
                  << mb.cS() << " m/s, cBar = " << mb.cBar() << " m/s)\n";
        return;
    }
    // the brazilian disc stands on its platens: no encastred face, no Lysmer
    if (scen_ == Scenario::TENSION || scen_ == Scenario::BRAZILIAN) return;

    std::string ab = cfg_.gets("absorbing", "none");
    if (ab != "none" && ab != "sides" && ab != "all")
        throw std::runtime_error("absorbing must be none | sides | all");

    double sF = cfg_.getd("absorbSpringFactor", 1.0);
    double Rside = cfg_.getd("absorbSpringR", 0.5 * W_);
    double Rbot  = cfg_.getd("absorbSpringR", H_);
    double tol = 1e-9;

    for (const auto& be : exterior_) {
        // impedances of the LOCAL phase: with mineral phases the truncated
        // continuum behind each boundary face is the one of the face's grain
        const Material& mp = phases_.mat[el_[be.elem].phase];
        double zP = mp.rho * mp.cP() * thk_;
        double zS = mp.rho * mp.cS() * thk_;
        double G  = mp.G();
        Eigen::Vector2d P = X0_[be.na], Q = X0_[be.nb];
        double L2 = 0.5 * (Q - P).norm();
        bool left  = P.x() < tol && Q.x() < tol;
        bool right = P.x() > W_ - tol && Q.x() > W_ - tol;
        bool bot   = P.y() < tol && Q.y() < tol;
        for (int nid : {be.na, be.nb}) {
            if ((left || right) && ab != "none") {
                cAbsX_[nid] += zP * L2;
                cAbsY_[nid] += zS * L2;
                kAbsX_[nid] += sF * G / Rside * L2 * thk_;
                kAbsY_[nid] += sF * G / (2.0 * Rside) * L2 * thk_;
            }
            if (bot) {
                if (ab == "all") {
                    cAbsY_[nid] += zP * L2;
                    cAbsX_[nid] += zS * L2;
                    kAbsY_[nid] += sF * G / Rbot * L2 * thk_;
                    kAbsX_[nid] += sF * G / (2.0 * Rbot) * L2 * thk_;
                } else if (scen_ == Scenario::PERCUSSION) {
                    flag_[nid] = FIXED;                // encastred bottom
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Brazilian load arcs. The two arcs are selected from the ORIGINAL exterior
// faces by the angular position of their midpoint (within loadArcDeg of the
// top and bottom poles) AND by the orientation of their outward normal: only
// the predominantly VERTICAL facets are loaded. That second test matters
// because cutting the disc out of a meshed box leaves a staircase rim whose
// vertical risers would otherwise take a horizontal traction and shear the rim
// for nothing; the horizontal treads are the bearing surface, and their total
// length is exactly the bearing width the classical solution assumes.
//
// Loading is by PRESSURE RATE (loadRate, Pa/s) rather than by displacement:
// the arcs are self-equilibrated, so there is no support to react against and
// nothing constrains the specimen's rigid-body motion. The default rate takes
// the arc pressure to 20 x ft over the run, which brackets failure for the
// usual bearing widths (sigma_t = 2P/(pi D t) with P = p x bearing width).
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Platen geometry and per-node bearing weights, shared by the brazilian and by
// the platen-loaded uniaxial test. The contact planes start FLUSH with the
// extreme nodes of the strip: no gap (the platen would accelerate into the
// specimen) and no initial penetration (which would release k pen^2/2 of energy
// created from nothing).
// ---------------------------------------------------------------------------
void FdemSolver::initPlatens(double xc, double hw) {
    // Platen contact stiffness. It MUST be set here and not in the caller: it
    // lived in the brazilian branch alone, so the uniaxial test ran with
    // kpPlaten_ = 0 and the platens transmitted nothing at all (peak 0.0 MPa,
    // zero broken joints — a silent "ok" that looked like a result).
    kpPlaten_ = cfg_.getd("platenPenaltyFactor", 1.0) * phases_.maxE() * thk_;
    if (!(kpPlaten_ > 0.0))
        throw std::runtime_error("platenPenaltyFactor must be > 0");
    // platen geometry FIRST: the tributary weights below test against
    // plTop_.xc, and computing them before it is set silently measures the
    // strip |x| <= hw at the LEFT EDGE of the box instead of the bearing.
    // That mistake cost four runs of a sweep: on a flattened disc
    // (hw = R sin alpha = 5.2 mm) the two strips do not even overlap, so
    // every weight was zero and the platens transmitted nothing.
    plTop_.xc = plBot_.xc = xc;
    plTop_.halfW = plBot_.halfW = hw;
    plTop_.sign = -1;                                  // above, presses down
    plBot_.sign = +1;                                  // below, presses up

    // tributary rim length per node, normalized to its mean over the
    // bearing: turns a per-duplicated-node spring into a per-unit-length
    // pressure, which is what a platen actually applies
    platenTrib_ = cfg_.getb("platenTributary", true);
    platenW_.assign(X0_.size(), 0.0);
    if (platenTrib_) {
        for (const auto& be : exterior_) {
            double L = (X0_[be.nb] - X0_[be.na]).norm();
            for (int nid : {be.na, be.nb}) {
                if (std::abs(X0_[nid].x() - plTop_.xc) > hw) continue;
                platenW_[nid] += 0.5 * L;
            }
        }
        double sum = 0.0;
        long n = 0;
        for (double w : platenW_) if (w > 0.0) { sum += w; ++n; }
        if (n == 0)
            throw std::runtime_error("brazilian: no exterior node under the "
                                     "platens — platenHalfWidth too small?");
        double mean = sum / n;
        for (double& w : platenW_) w /= mean;
    } else {
        std::fill(platenW_.begin(), platenW_.end(), 1.0);
    }
    double yHi = -1e300, yLo = 1e300;
    for (const auto& X : X0_) {
        if (std::abs(X.x() - plTop_.xc) > hw) continue;
        yHi = std::max(yHi, X.y());
        yLo = std::min(yLo, X.y());
    }
    if (yHi < yLo)
        throw std::runtime_error("brazilian: no node under the platens");
    // start flush with the extreme nodes: no gap (the platen would
    // accelerate into the disc) and no initial penetration (which would
    // release k pen^2/2 of energy created from nothing)
    plTop_.y = yHi;
    plBot_.y = yLo;
}

void FdemSolver::setupBrazilianLoad() {
    if (scen_ == Scenario::TENSION && tensionPlatens_) {
        // UCS / triaxial through platens: they span the FULL specimen width
        initPlatens(0.5 * W_, 0.5 * W_ + 1e-9);
        platenV_ = std::abs(cfg_.getd("pullV", 0.05));     // closure rate
        std::cout << "[FDEM] uniaxial loaded by PLATENS closing at "
                  << platenV_ << " m/s total (no clamped grip row)\n";
        return;
    }

    if (scen_ != Scenario::BRAZILIAN) return;
    std::string how = cfg_.gets("brazilianLoading", "platens");
    if (how != "platens" && how != "traction")
        throw std::runtime_error("brazilianLoading must be platens | traction "
                                 "(got '" + how + "')");
    brazPlatens_ = how == "platens";

    if (brazPlatens_) {
        // Two rigid platens closing on the disc, as in Y-Geo and the FEM-DEM
        // BTS literature. platenHalfWidth defaults to the full radius (flat
        // platens); the contact arc then grows with the load, which is what
        // keeps the CENTRE the most tensile point of the disc.
        // on a flattened disc the natural bearing is the flat itself
        double hwDef = discFlat_ > 0.0
                       ? discR_ * std::sin(0.5 * discFlat_ * M_PI / 180.0)
                       : discR_;
        double hw = cfg_.getd("platenHalfWidth", hwDef);
        if (!(hw > 0.0)) throw std::runtime_error("platenHalfWidth must be > 0");
        // Platen contact stiffness, separate from the tool's. Real brazilian
        // rigs put a cardboard or steel cushion between jaw and rock exactly to
        // spread the contact; numerically the same job is done by softening the
        // penalty, which lets more rim nodes engage instead of one asperity
        // taking the whole load.
        initPlatens(discC_.x(), hw);
        platenV_ = std::abs(cfg_.getd("pullV", 0.05));      // CLOSURE rate
        std::cout << "[FDEM] brazilian: flattening 2*alpha = " << discFlat_
                  << " deg, platen penalty " << kpPlaten_ / (phases_.maxE() * thk_)
                  << " x E*t\n";
        std::cout << "[FDEM] brazilian: two rigid platens closing at "
                  << platenV_ << " m/s total (" << 0.5 * platenV_
                  << " m/s each, both inward), half-width " << hw << " m\n";
        return;
    }

    double halfAng = cfg_.getd("loadArcDeg", 7.5) * M_PI / 180.0;
    if (!(halfAng > 0.0 && halfAng < 0.5 * M_PI))
        throw std::runtime_error("loadArcDeg must be in (0, 90)");
    loadRate_ = cfg_.getd("loadRate", 20.0 * mat_.ft / T_);
    if (!(loadRate_ > 0.0))
        throw std::runtime_error("loadRate must be > 0 [Pa/s]");

    for (const auto& be : exterior_) {
        Eigen::Vector2d P = X0_[be.na], Q = X0_[be.nb];
        Eigen::Vector2d mid = 0.5 * (P + Q);
        Eigen::Vector2d d = Q - P;
        double L = d.norm();
        if (L < 1e-14) continue;
        Eigen::Vector2d n(d.y() / L, -d.x() / L);      // outward
        if (std::abs(n.y()) < 0.5) continue;           // staircase riser
        Eigen::Vector2d r = mid - discC_;
        double rn = r.norm();
        if (rn < 1e-14) continue;
        // angle away from the vertical: |cos| of the angle between r and +/-y
        double cosToPole = std::abs(r.y()) / rn;
        if (cosToPole < std::cos(halfAng)) continue;
        if (r.y() > 0.0 && n.y() > 0.0) {
            arcTop_.edge.push_back(be);
            arcTop_.length += L;
        } else if (r.y() < 0.0 && n.y() < 0.0) {
            arcBot_.edge.push_back(be);
            arcBot_.length += L;
        }
    }
    if (arcTop_.edge.empty() || arcBot_.edge.empty())
        throw std::runtime_error("brazilian: no loadable rim facet found in one "
                                 "of the arcs — widen loadArcDeg or refine the "
                                 "mesh");
    std::cout << "[FDEM] brazilian load arcs: +/-" << cfg_.getd("loadArcDeg", 7.5)
              << " deg, " << arcTop_.edge.size() << " / " << arcBot_.edge.size()
              << " faces, bearing width " << arcTop_.length << " / "
              << arcBot_.length << " m (" << 100.0 * arcTop_.length
              / (2.0 * discR_) << " % of D), rate " << loadRate_ / 1e6
              << " MPa/s\n";
    if (std::abs(arcTop_.length - arcBot_.length)
        > 0.2 * arcTop_.length)
        std::cout << "[FDEM] WARNING: the two bearing widths differ by more "
                     "than 20 %: the load pair is not balanced and the disc "
                     "will drift. Refine the mesh or widen loadArcDeg.\n";
}

// ---------------------------------------------------------------------------
// Extensometre de la compression uniaxiale/triaxiale (Yan et al. fig. 19b).
// Trois mesures de la deformation axiale, ecrites cote a cote dans
// history.csv pour que la comparaison au module d'entree soit verifiable :
//   epsPlaten : la fermeture des plateaux rapportee a leur ecartement initial
//               — la deformation MACHINE, qui contient la compliance du
//               contact penalite ;
//   epsSpec   : le raccourcissement des deux faces de l'eprouvette ;
//   epsGauge  : un extensometre interieur entre deux bandes de nœuds a
//               gaugeLoFrac et gaugeHiFrac de la hauteur, donc affranchi du
//               contact ET des effets de bord.
// Aucune force n'en depend : sortie seulement.
// ---------------------------------------------------------------------------
void FdemSolver::setupStrainGauge() {
    if (!(scen_ == Scenario::TENSION && tensionPlatens_)) return;
    gap0_ = plTop_.y - plBot_.y;
    if (!(gap0_ > 0.0))
        throw std::runtime_error("platen gap is not positive — initPlatens "
                                 "did not find the bearing rows");
    // les deux bandes de l'extensometre : une demi-maille de part et d'autre
    // des cotes y = gLoFrac*H et y = gHiFrac*H
    double band = 0.5 * hmin_ + 1e-12;
    double yLo = gLoFrac_ * H_, yHi = gHiFrac_ * H_;
    double sLo = 0.0, sHi = 0.0;
    for (int i = 0; i < (int)X0_.size(); ++i) {
        double y = X0_[i].y();
        if (std::abs(y - yLo) <= band) { gLoNodes_.push_back(i); sLo += y; }
        if (std::abs(y - yHi) <= band) { gHiNodes_.push_back(i); sHi += y; }
        if (y < 1e-9) botNodes_.push_back(i);
        if (y > H_ - 1e-9) topNodes_.push_back(i);
    }
    if (gLoNodes_.empty() || gHiNodes_.empty())
        throw std::runtime_error("strain gauge: no node found on one of the "
                                 "extensometer bands — widen the band or move "
                                 "gaugeLoFrac / gaugeHiFrac");
    gLoY_ = sLo / gLoNodes_.size();
    gHiY_ = sHi / gHiNodes_.size();
    std::cout << "[FDEM] axial strain gauge: platen gap " << gap0_
              << " m; extensometer " << gHiY_ - gLoY_ << " m between "
              << gLoNodes_.size() << " and " << gHiNodes_.size()
              << " nodes at y = " << gLoY_ << " / " << gHiY_ << " m\n";
}

void FdemSolver::gaugeStrain(double& epsPlaten, double& epsSpec,
                             double& epsGauge) const {
    epsPlaten = epsSpec = epsGauge = 0.0;
    if (gap0_ <= 0.0) return;
    epsPlaten = (gap0_ - (plTop_.y - plBot_.y)) / gap0_;
    auto meanUy = [&](const std::vector<int>& ns) {
        double s = 0.0;
        for (int i : ns) s += u_[i].y();
        return ns.empty() ? 0.0 : s / ns.size();
    };
    if (!topNodes_.empty() && !botNodes_.empty())
        epsSpec = (meanUy(botNodes_) - meanUy(topNodes_)) / H_;
    double L0 = gHiY_ - gLoY_;
    if (L0 > 0.0)
        epsGauge = (meanUy(gLoNodes_) - meanUy(gHiNodes_)) / L0;
}

// ---------------------------------------------------------------------------
// Confining pressure. The loaded set is fixed once, at init, from the ORIGINAL
// exterior faces:
//   box  + confineFaces = sides : the two lateral faces (x = 0, x = W) — the
//                                 triaxial cell, and the confined-percussion
//                                 case;
//   box  + confineFaces = all   : every exterior face except the ones a
//                                 boundary condition already owns (the tension
//                                 grips; the struck top face in percussion —
//                                 the tool is there, not a fluid);
//   disc                        : the whole perimeter minus the platen strips;
//   confineFaces = bore         : ONLY the original exterior faces whose
//                                 midpoint lies within boreSelectR of
//                                 (boreCX, boreCY) — a pressurized cavity
//                                 (tunnel / borehole). The same follower load
//                                 pushes the cavity wall INTO the solid, and
//                                 faces born from cracking receive nothing,
//                                 so the fluid never enters the fissures.
// ---------------------------------------------------------------------------
void FdemSolver::setupConfinement() {
    confP_ = cfg_.getd("confiningPressure", 0.0);
    confRamp_ = cfg_.getd("confiningRamp", 0.0);
    if (confP_ == 0.0) return;
    if (confP_ < 0.0)
        throw std::runtime_error("confiningPressure must be >= 0 (it is a "
                                 "PRESSURE: positive squeezes the specimen)");
    std::string which = cfg_.gets("confineFaces", "sides");
    if (which != "sides" && which != "all" && which != "bore")
        throw std::runtime_error("confineFaces must be sides | all | bore "
                                 "(got '" + which + "')");
    bool all = which == "all";
    bool bore = which == "bore";
    double tol = 1e-9;
    double hw = cfg_.getd("platenHalfWidth", discR_);
    double bcx = cfg_.getd("boreCX", 0.5 * W_);
    double bcy = cfg_.getd("boreCY", 0.5 * H_);
    double bsr = cfg_.getd("boreSelectR", 0.0);
    if (bore && bsr <= 0.0)
        throw std::runtime_error("confineFaces = bore requires boreSelectR > 0 "
                                 "(faces whose midpoint is within that radius "
                                 "of boreCX/boreCY are pressurized)");

    for (const auto& be : exterior_) {
        Eigen::Vector2d P = X0_[be.na], Q = X0_[be.nb];
        if (bore) {
            Eigen::Vector2d mid = 0.5 * (P + Q);
            double dx = mid.x() - bcx, dy = mid.y() - bcy;
            if (dx * dx + dy * dy <= bsr * bsr) confEdges_.push_back(be);
            continue;
        }
        if (disc_) {
            // skip the faces the platens will press on: a fluid pressure and a
            // rigid platen on the same facet is a double load
            double ym = 0.5 * (P.y() + Q.y());
            bool underPlaten = std::abs(0.5 * (P.x() + Q.x()) - discC_.x()) <= hw
                               && (ym > discC_.y() + 0.9 * discR_
                                   || ym < discC_.y() - 0.9 * discR_);
            if (!underPlaten) confEdges_.push_back(be);
            continue;
        }
        bool left  = P.x() < tol && Q.x() < tol;
        bool right = P.x() > W_ - tol && Q.x() > W_ - tol;
        bool bot   = P.y() < tol && Q.y() < tol;
        bool top   = P.y() > H_ - tol && Q.y() > H_ - tol;
        if (left || right) { confEdges_.push_back(be); continue; }
        if (!all) continue;
        if (scen_ == Scenario::TENSION && (bot || top)) continue;  // grips
        if (scen_ != Scenario::TENSION && top) continue;           // tool face
        if (bot && (flag_[be.na] == FIXED || flag_[be.nb] == FIXED)) continue;
        confEdges_.push_back(be);
    }
    if (confEdges_.empty())
        throw std::runtime_error("confiningPressure is set but no exterior face "
                                 "qualifies — check confineFaces / geometry");
    for (const auto& be : confEdges_)
        confL0_ += (X0_[be.nb] - X0_[be.na]).norm();

    std::cout << "[FDEM] confinement: " << confP_ / 1e6 << " MPa on "
              << confEdges_.size() << " exterior faces (" << which
              << "), total length " << confL0_ << " m, ramp " << confRamp_
              << " s\n";
    if (confRamp_ <= 0.0)
        std::cout << "[FDEM] WARNING: confiningRamp = 0 applies the pressure as "
                     "a STEP, which launches a wave through the specimen "
                     "(same failure mode as a stepped grip velocity — see "
                     "pullRamp). Ramp over several wave transits.\n";
    for (const auto& be : confEdges_)
        if (kAbsX_[be.na] > 0 || kAbsY_[be.na] > 0) {
            std::cout << "[FDEM] WARNING: confined faces also carry Lysmer "
                         "springs (absorbing): the springs pull the surface "
                         "back toward u = 0 and FIGHT the confinement. Read the "
                         "achieved lateral stress in the summary before "
                         "trusting the target.\n";
            break;
        }
}

// ---------------------------------------------------------------------------
// Stable time step. Per-node stiffness sums (joints dominate through the
// intrinsic penalty):  elements ~ 2 E t per node, each joint point
// p * (L/2) * t on its two nodes, plus a budget of general contacts at kp.
//   dt = dtFactor * min_i 2 sqrt(m_i / K_i),   also capped by the mesh CFL.
// ---------------------------------------------------------------------------
void FdemSolver::computeStableDt() {
    std::vector<double> K(X0_.size());
    for (std::size_t i = 0; i < X0_.size(); ++i)
        K[i] = 2.0 * phases_.mat[el_[elemOf_[i]].phase].E * thk_;
    for (const auto& J : jt_) {
        double k = J.pj * 0.5 * J.L0 * thk_;
        K[J.a1] += k; K[J.a2] += k; K[J.b1] += k; K[J.b2] += k;
    }
    double nExtra = cfg_.getd("extraContacts", 2.0);
    // the platen penalty is a spring on the bearing nodes exactly like the
    // tool's: budget the STIFFER of the two, or a platenPenaltyFactor > 1
    // would silently control stability
    // A1 : en contact de Signorini l'outil n'a plus de raideur — il impose une
    // condition sur la VITESSE. kp_ sort donc du budget, et c'est le gain
    // structurel du schema : le pas de temps cesse de dependre d'une penalite
    // arbitraire. Les platines, elles, restent en penalite.
    double kContact = toolSig_ ? kpPlaten_ : std::max(kp_, kpPlaten_);
    // SHPB: the bar/rock interfaces are the ONLY load path, so the general
    // contact penalty is a first-class spring and must enter the budget (it
    // does not for the other scenarios, where general contact only handles
    // debris at gcPenaltyFactor = 0.01 E t and never controls stability).
    if (scen_ == Scenario::SHPB) kContact = std::max(kContact, kpGC_);
    // Contact de debris REELLEMENT employe. En mode potentiel c'est potP_ /
    // potKt_ ; kpGC_ n'y est jamais evalue (generalContact() derive vers
    // potentialContact() avant de le lire), donc le budget ci-dessus ne
    // parlait pas de la bonne raideur. A potPenaltyFactor = 1 la valeur vaut
    // exactement kp_ et ce max ne change rien — mais il SUIT desormais la
    // cle dans les deux sens.
    if (contactPot_) kContact = std::max(kContact, std::max(potP_, potKt_));
    double dtMin = 1e30;
    for (std::size_t i = 0; i < X0_.size(); ++i)
        dtMin = std::min(dtMin,
                         2.0 * std::sqrt(m_[i] / (K[i] + nExtra * kContact)));
    // CFL from the TRUE minimum inscribed diameter (4A/P), not the nominal
    // grid pitch: the cross-diagonal split makes triangles with h ~ 0.41 dx,
    // so nominal hmin overestimates the ceiling ~2.4x. Masked in intrinsic
    // mode by the joint-spring term dominating dtMin; exposed by
    // insertion = adaptive (see the 3D solver, where the same defect blew up
    // the first homogeneous grid impact, 2026-08-07).
    double hCfl = hmin_;
    for (double h : hEl_) hCfl = std::min(hCfl, h);
    double cfl = hCfl / phases_.maxCp();
    // ---- borne DIFFUSIVE du terme visqueux 2 mu D (eq. 6 de Yan) --------
    // Un terme de Kelvin-Voigt ajoute sa PROPRE contrainte de stabilite au
    // schema explicite, que le budget elastique ci-dessus ignore. Pour une
    // barre 1D de longueur h, la force visqueuse equivaut a un amortisseur
    // c = 2 mu A / h sur une masse nodale m = rho A h / 2, et le critere de
    // la difference centree tend vers dt <= 2 m / c = rho h^2 / (2 mu) quand
    // l amortissement domine. Sans cette borne, monter bulkViscosity finit
    // par faire diverger le calcul SANS que le pas de temps ne bouge — le
    // piege exact que la cle est censee eviter.
    // Forme SERREE : la viscosite longitudinale effective est 2 mu, donc la
    // diffusivite de quantite de mouvement vaut nu = 2 mu / rho et le critere
    // explicite classique dt <= h^2/(2 nu) donne dt <= rho h^2 / (4 mu).
    // (Le raisonnement par amortisseur equivalent 2m/c donne rho h^2/(2 mu),
    // deux fois plus permissif : on retient le plus severe des deux.)
    double dtVis = 1e30;
    if (bulkVisc_ > 0.0)
        for (int eI = 0; eI < (int)el_.size(); ++eI)
            dtVis = std::min(dtVis, rhoP_[el_[eI].phase] * hEl_[eI] * hEl_[eI]
                                    / (4.0 * bulkVisc_));
    dt_ = cfg_.getd("dtFactor", 0.2) * std::min(std::min(dtMin, cfl), dtVis);
    if (bulkVisc_ > 0.0)
        std::cout << "[FDEM] viscosite newtonienne (Yan eq. 6) : mu = "
                  << bulkVisc_ << " Pa.s ; borne diffusive rho h^2/(4 mu) = "
                  << dtVis << " s contre " << std::min(dtMin, cfl)
                  << " s pour les ressorts"
                  << (dtVis < std::min(dtMin, cfl)
                      ? "  <-- C EST ELLE QUI COMMANDE LE PAS" : "")
                  << "\n";
}

// ===========================================================================
// Time stepping
// ===========================================================================

void FdemSolver::step() {
    for (auto& fi : f_) fi.setZero();
    tool_.resetForce();
    // balai numerique : armement a l'instant demande (voir FdemSolver.hpp)
    if (brushStart_ > 0.0 && !brushArmed_ && t_ >= brushStart_) armBrush();
    // SHPB: the pulse velocity of THIS step, read by integrate() (DRIVEX)
    if (scen_ == Scenario::SHPB) vDrive_ = shpbVel(t_);

    if (fProf.on) {
        double t0 = fnow(); elementForces(); bodyForces();
        double t05 = fnow(); if (adaptive_) insertionSweep();
        double t1 = fnow(); jointForces();
        double t2 = fnow(); generalContact();
        double t3 = fnow(); toolContact();
        double t4 = fnow();
        fProf.tEl += t05 - t0; fProf.tIn += t1 - t05; fProf.tJt += t2 - t1;
        fProf.tGc += t3 - t2; fProf.tTc += t4 - t3;
        ++fProf.n;
    } else {
        elementForces();
        bodyForces();                      // gravity (no-op when gravity = 0)
        if (adaptive_) insertionSweep();   // before jointForces: a joint born
                                           // this step carries traction now
        jointForces();
        generalContact();
        toolContact();
    }
    if (scen_ == Scenario::SHPB) shpbGaugeRead();   // monitor points 1 and 2
    brazilianForces();                     // no-op outside the brazilian
    if (scen_ == Scenario::TENSION && tensionPlatens_) platenForces();
    confiningForces();                     // no-op when confiningPressure = 0
    excavationForces();                    // no-op si excavRelease = false
    hydroForces();                         // no-op si hydro = off (spec 004)
    // gauge the confinement AFTER the ramp has had time to equilibrate through
    // the specimen: read exactly at the end of the ramp and the interior is
    // still catching up (measured 55 % of the target that way, 3 ramp times
    // later it is there). Keep it well before the axial load matters.
    if (confP_ > 0.0 && !confLatched_
        && t_ >= std::max(cfg_.getd("confineGaugeTime", 3.0 * confRamp_),
                          20.0 * dt_)) {
        confAchieved_ = achievedConfinement();
        confLatched_ = true;
    }

    if (scen_ == Scenario::TENSION) {
        if (tensionPlatens_) {
            // the platen reaction IS the axial load; no grip row exists
            gripF_ = plTop_.F;
        } else {
            gripF_.setZero();
            for (int i = 0; i < (int)X0_.size(); ++i)
                if (flag_[i] == PRESCRIBED) gripF_ += f_[i];
        }
        double sigma = std::abs(gripF_.y()) / (W_ * thk_);
        // Verrouillage du pic, transpose du bresilien : tant que la premiere
        // chute franche post-fissuration n'a pas eu lieu, l'essai est une
        // mesure de resistance valable ; apres, les plateaux ecrasent les
        // fragments et la contrainte nominale n'a plus de sens. Le seuil est
        // 30 % du pic (l'article coupe ses courbes de la fig. 19b sur la
        // branche descendante). Sortie seulement : aucune force n'en depend.
        if (!peakLockedU_) {
            sigmaPeak_ = std::max(sigmaPeak_, sigma);
            if (nBroken_ > 0 && sigmaPeak_ > 0.0 && sigma < 0.3 * sigmaPeak_) {
                peakLockedU_ = true;
                tLockedU_ = t_;
            }
        }
    } else if (scen_ == Scenario::BRAZILIAN) {
        // ISRM indirect tensile strength from the platen reaction:
        //   sigma_t = 2 P / (pi D t)
        // the platen reaction and the applied traction play the same role: the
        // vertical force squeezing the disc along its diameter
        double P = brazPlatens_ ? std::abs(plTop_.F.y())
                                : std::abs(arcTop_.F.y());
        // Flattened brazilian disc: sigma_t = k * 2P/(pi D t), the correction
        // k of Wang et al. (IJRMMS 41, 2004) as tabulated by Yu et al. (EFM
        // 247, 2021): 0.9644 at 2*alpha = 20 deg, 0.9205 at 30 deg. Linear in
        // between, 1 for the plain disc. Wang's validity criterion is
        // 2*alpha >= 20 deg — below that the opening stress is no longer
        // maximal at the centre.
        double kFBD = 1.0;
        if (discFlat_ > 0.0)
            kFBD = discFlat_ <= 20.0
                   ? 1.0 + (0.9644 - 1.0) * discFlat_ / 20.0
                   : 0.9644 + (0.9205 - 0.9644) * (discFlat_ - 20.0) / 10.0;
        sigmaT_ = kFBD * 2.0 * P / (M_PI * 2.0 * discR_ * thk_);
        if (!peakLocked_) {
            if (sigmaT_ > sigmaTpeak_) {
                sigmaTpeak_ = sigmaT_;
                peakF_ = P;
                tPeak_ = t_;
                nBrokenAtPeak_ = nBroken_;
            }
            // a drop to 70 % of the peak once cracking has started is the
            // failure event; everything after it is the platens grinding the
            // halves together
            if (nBroken_ > 0 && sigmaTpeak_ > 0.0
                && sigmaT_ < 0.7 * sigmaTpeak_) {
                peakLocked_ = true;
                tLocked_ = t_;                 // end-of-test marker
            }
        }
        // EARLY elastic gauge, accumulated while the nominal sigma_t is
        // inside [eGaugeLo_, eGaugeHi_] x ft and the disc is still intact.
        // That band is the genuinely linear part of the loading: the centre
        // there carries sigma_xx = 2P/(pi D t) to better than a percent, while
        // above it the core starts to soften and the closed-form solution stops
        // describing the disc (measured: the ratio holds 0.99 up to
        // sigma_t ~ 1.2 x ft, then falls away as the centre unloads).
        if (nBroken_ == 0 && sigmaT_ >= eGaugeLo_ * mat_.ft
            && sigmaT_ <= eGaugeHi_ * mat_.ft) {
            double sx = 0.0, sy = 0.0;
            discCentreStress(sx, sy);
            eSumXX_ += sx;
            eSumYY_ += sy;
            eSumSig_ += sigmaT_;
            ++eN_;
            eDfrac_ = gDfrac_;
        }
        // latch the elastic gauge while the disc is still intact; the LAST
        // such step is the one just before first breakage, i.e. the highest
        // load at which the closed-form solution still applies
        if (nBroken_ == 0) {
            discCentreStress(gXX_, gYY_);
            gSigT_ = sigmaT_;
            long nD = 0;
            double sD = 0.0;
            for (const auto& J : jt_) {
                sD += J.D;
                if (J.D > 0.01) ++nD;
            }
            gDfrac_ = jt_.empty() ? 0.0 : (double)nD / jt_.size();
            gDmean_ = jt_.empty() ? 0.0 : sD / jt_.size();
        }
    } else {
        peakF_ = std::max(peakF_, tool_.F.norm());
        work_ += -tool_.F.dot(tool_.v) * dt_;
    }

    integrate();
    t_ += dt_;

    if ((++stepCount_ & 1023) == 0) {                  // cheap stability guard
        checkEnergyAbort();                // opt-in (budgetAbortPct), E2
        // ---- E5 (2026-08-19) : la garde testait u_[0].x(), or le noeud 0
        // peut etre FIXED — donc rigoureusement nul, donc toujours fini. Le
        // detecteur etait AVEUGLE a une divergence qui n'aurait pas touche ce
        // noeud precis. On echantillonne desormais tout le maillage a pas
        // constant (~256 noeuds), avec un decalage qui tourne d'un controle a
        // l'autre pour couvrir l'integralite des noeuds au fil du run. Cout :
        // 256 tests tous les 1024 pas. Lecture PURE, aucun flottant ne change.
        bool bad = !std::isfinite(work_);
        const std::size_t nN = X0_.size();
        const std::size_t stride = (nN > 256) ? nN / 256 : 1;
        const std::size_t off = (std::size_t)((stepCount_ >> 10) % (long)stride);
        for (std::size_t i = off; i < nN && !bad; i += stride)
            if (!std::isfinite(u_[i].x()) || !std::isfinite(u_[i].y()))
                bad = true;
        if (bad)
            throw std::runtime_error("FDEM instability (NaN) — reduce dtFactor");
    }
}

// ---------------------------------------------------------------------------
// Co-rotational CST internal forces. F = sum_a x_a (grad N_a)^T; 2D polar
// decomposition in closed form: with c1 = F00 + F11, c2 = F10 - F01,
//   R = [[c1, -c2], [c2, c1]] / sqrt(c1^2 + c2^2).
// Biot strain eps = sym(R^T F) - I, sigma = D eps in the co-rotated frame,
// nodal forces f_a -= A0 t (R sigma) grad N_a. Exact for arbitrary rigid
// rotations, small elastic strains — which is the regime of flying rock
// fragments.
// ---------------------------------------------------------------------------
void FdemSolver::elementForces() {
    // node duplication makes this loop embarrassingly parallel: every
    // element writes only its OWN three nodes
    double wEl = 0.0;                      // V2/B4 : travail de ce pas
    double wVi = 0.0;                      // dont part VISQUEUSE (ventilation)
    double wBd = 0.0;                      // dont PULVERISATION (WP1, energie)
    long nPv = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:wEl,wVi,wBd,nPv)
#endif
    for (int eI = 0; eI < (int)el_.size(); ++eI) {
        Elem& e = el_[eI];
        const Eigen::Matrix3d& Dm = DmP_[e.phase];
        Eigen::Matrix2d F = Eigen::Matrix2d::Zero();
        for (int a = 0; a < 3; ++a) {
            Eigen::Vector2d x = X0_[e.n[a]] + u_[e.n[a]];
            F += x * e.dN.col(a).transpose();
        }
        double c1 = F(0, 0) + F(1, 1);
        double c2 = F(1, 0) - F(0, 1);
        double hyp = std::sqrt(c1 * c1 + c2 * c2);
        Eigen::Matrix2d R;
        if (hyp > 1e-14) R << c1 / hyp, -c2 / hyp, c2 / hyp, c1 / hyp;
        else             R.setIdentity();
        Eigen::Matrix2d U = R.transpose() * F;
        Eigen::Vector3d eps(U(0, 0) - 1.0, U(1, 1) - 1.0, U(0, 1) + U(1, 0));
        Eigen::Vector3d s;
        if (law_) {
            // PLANE STRAIN is exactly eps_zz = 0, so the 3D law can be used
            // verbatim: hand it the co-rotated Biot strain with a zero
            // out-of-plane component and keep the in-plane stresses. The law
            // returns sigma_zz itself (it is a reaction, not an input).
            Eigen::Matrix3d E3 = Eigen::Matrix3d::Zero();
            E3(0, 0) = eps(0);
            E3(1, 1) = eps(1);
            E3(0, 1) = E3(1, 0) = 0.5 * eps(2);        // tensorial shear
            Eigen::Matrix3d S3 = law_->stress(E3, e.st, dt_, hEl_[eI]);
            s << S3(0, 0), S3(1, 1), S3(0, 1);
        } else {
            s = Dm * eps;
        }
        double szz = nuP_[e.phase] * (s(0) + s(1));
        double pm = (s(0) + s(1) + szz) / 3.0;
        e.svm = std::sqrt(1.5 * ((s(0) - pm) * (s(0) - pm)
                                 + (s(1) - pm) * (s(1) - pm)
                                 + (szz - pm) * (szz - pm) + 2.0 * s(2) * s(2)));
        // Crush cap (elastic-perfectly-plastic ceiling on the deviator, plus
        // a cap on mean tension). FDEM has no erosion: crushed elements under
        // the tool otherwise reach O(1) strain where the geometric stiffness
        // is no longer covered by the linear stable dt, and the explicit
        // update pumps energy (measured: 1e6 J from a 148 J impact). Bounding
        // the deviatoric stress bounds the storable elastic energy and turns
        // over-compression into plastic dissipation. The cap is far above
        // every legitimate wave amplitude in the demos.
        double cap = law_ ? 1e300 : crushCapP_[e.phase];
        if (e.svm > cap) {
            double fscale = cap / e.svm;
            s(0) = pm + (s(0) - pm) * fscale;
            s(1) = pm + (s(1) - pm) * fscale;
            s(2) *= fscale;
            szz = pm + (szz - pm) * fscale;
            e.svm = cap;
        }
        if (!law_ && mtCap_ > 0.0 && pm > mtCap_ * ftP_[e.phase]) {
            double shift = pm - mtCap_ * ftP_[e.phase];
            s(0) -= shift; s(1) -= shift;
        }
        // ---- WP1 : pulverisation (Yang et al. 2026, eq. 3-4), miroir 2D --
        // delta_m = h_e * eps_vm en DEFORMATION PLANE : le deviateur inclut
        // eps_zz = 0 ; eps(2) est le cisaillement d INGENIEUR (gamma), la
        // composante tensorielle vaut gamma/2. Dissipation Y dD avec
        // Y = 1/2 Cd s:eps = 1/2 Cd (s0 e0 + s1 e1 + s2 g) — le produit de
        // Voigt rend la double contraction exacte. Applique APRES les caps,
        // AVANT la viscosite, comme en 3D.
        if (bdOn_ && !law_) {
            double m3 = (eps(0) + eps(1)) / 3.0;       // tr/3, eps_zz = 0
            double d2 = (eps(0) - m3) * (eps(0) - m3)
                      + (eps(1) - m3) * (eps(1) - m3) + m3 * m3
                      + 0.5 * eps(2) * eps(2);         // 2 (gamma/2)^2
            double dm = hEl_[eI] * std::sqrt(2.0 / 3.0 * d2);
            if (dm > e.bdDm) e.bdDm = dm;
            if (e.bdDm > bdD0_) {
                double D = bdDf_ * (e.bdDm - bdD0_)
                         / (e.bdDm * (bdDf_ - bdD0_));
                if (D > bdDmax_) D = bdDmax_;
                if (D > e.bdD)
                    wBd -= 0.5 * bdCd_
                         * (s(0) * eps(0) + s(1) * eps(1) + s(2) * eps(2))
                         * e.A0 * thk_ * (D - e.bdD);
                e.bdD = D;
                double k = bdCd_ * (1.0 - D);
                s *= k;
                e.svm *= k;
                if (D >= bdDmax_) ++nPv;
            }
        }
        // ---- viscosite NEWTONIENNE ISOTROPE (eq. 6 de Yan : + 2 mu D) -----
        // NB de vocabulaire : ce n est PAS une  viscosite de volume  au sens
        // zeta tr(D) I (le bulk viscosity d Abaqus). Le terme agit sur le
        // tenseur COMPLET, trace comprise, donc il amortit aussi bien le
        // volumique que le deviatorique. Confondre les deux fausserait toute
        // calibration.
        // D = sym(L), L = Fdot F^-1 (taux de deformation en configuration
        // courante), tourne dans le repere co-rote et ajoute a la contrainte
        // AVANT la rotation retour — l'equivalent exact du 2*mu*D de leur
        // forme, dissipatif par construction (puissance 2 mu D:D >= 0).
        // bulkViscosity = 0 (defaut) : branche non executee, bit-identique.
        if (bulkVisc_ > 0.0 || difOn_) {
            Eigen::Matrix2d Fd = Eigen::Matrix2d::Zero();
            for (int a = 0; a < 3; ++a)
                Fd += v_[e.n[a]] * e.dN.col(a).transpose();
            // Garde : F est inverse ici, ce que ne fait aucun autre chemin
            // du solveur 2D. Un element retourne ou ecrase (det F <= 0) rend
            // l inverse infinie et propagerait des NaN dans TOUTE la
            // contrainte. On saute l element pour ce pas — il garde son
            // taux filtre precedent, ce qui est le comportement sur.
            double detF = F.determinant();
            if (!(detF > 1e-12)) continue;
            Eigen::Matrix2d L = Fd * F.inverse();
            Eigen::Matrix2d Dr = 0.5 * (L + L.transpose());
            Eigen::Matrix2d Dc = R.transpose() * Dr * R;   // co-rote
            if (bulkVisc_ > 0.0) {
                s(0) += 2.0 * bulkVisc_ * Dc(0, 0);
                s(1) += 2.0 * bulkVisc_ * Dc(1, 1);
                s(2) += 2.0 * bulkVisc_ * Dc(0, 1);
                // VENTILATION de la dissipation visqueuse : la puissance
                // 2 mu D:D par unite de volume est >= 0 par construction, donc
                // comptee NEGATIVEMENT comme tout ce qui quitte l energie
                // cinetique. Elle est DEJA dans wEl (la contrainte visqueuse
                // produit les memes forces nodales) : c est une ventilation,
                // pas un poste supplementaire du bilan. En deformation plane
                // D_zz = 0, donc la norme de Frobenius 2D EST D:D.
                wVi -= 2.0 * bulkVisc_ * Dc.squaredNorm() * e.A0 * thk_;
            }
            if (difOn_) {
                // principale MAXIMALE en valeur absolue du tenseur taux, puis
                // filtre exponentiel de constante srTau_ (voir le header : le
                // taux brut par element est trop bruite pour figer un DIF).
                double trD = Dc(0, 0) + Dc(1, 1);
                double haD = 0.5 * (Dc(0, 0) - Dc(1, 1));
                double rad = std::sqrt(haD * haD + Dc(0, 1) * Dc(0, 1));
                double lm1 = 0.5 * trD + rad;
                double lm2 = 0.5 * trD - rad;
                double er  = std::max(std::abs(lm1), std::abs(lm2));
                e.edot = srRelax_ * e.edot + (1.0 - srRelax_) * er;
            }
        }
        Eigen::Matrix2d sig;
        sig << s(0), s(2), s(2), s(1);
        Eigen::Matrix2d P = R * sig;                   // rotate back
        // ---- contrainte in situ ---------------------------------------
        // sigma_global_total = R sig R^T + sigma0, et la force interne
        // s'ecrit avec P = sigma_global * R : il suffit donc d'ajouter
        // sigma0 * R ici. sigG (calcule juste en dessous par P * R^T) porte
        // alors la contrainte TOTALE, et c'est elle que lisent la sortie
        // VTU, la jauge achievedConfinement() et — surtout — le critere
        // d'insertion adaptative de insertionSweep().
        if (hasInsitu_) P += insituS_ * R;
        // sigma_xx in the GLOBAL frame (output only — the confinement gauge);
        // P * R^T is the Cauchy stress rotated out of the co-rotated frame
        Eigen::Matrix2d sigG = P * R.transpose();
        e.sxx = sigG(0, 0);
        e.syy = sigG(1, 1);
        e.sxy = 0.5 * (sigG(0, 1) + sigG(1, 0));
        e.exx = eps(0);                            // axial strain (gauges)
        for (int a = 0; a < 3; ++a) {
            Eigen::Vector2d fe = e.A0 * thk_ * (P * e.dN.col(a));
            f_[e.n[a]] -= fe;
            wEl -= fe.dot(v_[e.n[a]]);     // V2/B4 : lecture pure
        }
    }
    elWork_ += wEl * dt_;
    viscWork_ += wVi * dt_;                // ventilation, incluse dans elWork_
    bdWork_ += wBd;                        // WP1 : deja une ENERGIE (Y dD)
    nPulv_ = nPv;
}

// ---------------------------------------------------------------------------
// Body force (gravity). Each element hands rho * A0 * thk * g to its three
// OWN nodes, one third each — the same lumping the mass matrix uses
// (buildMesh: m_[n] += rho * A0 * thk / 3), so the nodal weight is exactly
// m_[n] * g however the nodes are duplicated or bound. It acts along -y and
// is accumulated into f_ with the SAME sign convention as every other
// external force (f_ is the total nodal force entering v += dt/m * f).
//
// Parallelised exactly like elementForces(): node duplication means every
// element writes only to its own three nodes, so the element loop has no
// write conflict. Works unchanged in adaptive mode — integrate() sums f_ and
// m_ over each binding group, and both sums scale identically.
//
// gravity = 0 (the default) short-circuits: no existing model changes.
// ---------------------------------------------------------------------------
void FdemSolver::bodyForces() {
    if (gravity_ > 0.0) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int eI = 0; eI < (int)el_.size(); ++eI) {
            const Elem& e = el_[eI];
            double w = rhoP_[e.phase] * e.A0 * thk_ * gravity_ / 3.0;
            for (int a = 0; a < 3; ++a) f_[e.n[a]].y() -= w;
        }
    }
    // ---- balai numerique : anti-gravite sur les SEULS candidats -----------
    // Serie et non parallele : le travail s'accumule dans un scalaire, et la
    // phase de balayage represente quelques pour cent du run. Voir
    // FdemSolver.hpp pour la justification de la comptabilite separee.
    if (!brushArmed_) return;
    double bw = 0.0;
    for (int eI = 0; eI < (int)el_.size(); ++eI) {
        if (!brushCand_[eI]) continue;
        const Elem& e = el_[eI];
        double w = rhoP_[e.phase] * e.A0 * thk_ * brushA_ / 3.0;
        for (int a = 0; a < 3; ++a) {
            Eigen::Vector2d F = w * brushDir_;
            f_[e.n[a]] += F;
            bw += F.dot(v_[e.n[a]]) * dt_;
        }
    }
    brushWork_ += bw;                      // POSTE SEPARE, jamais dans sumW
}

// ---------------------------------------------------------------------------
// Armement du balai : fige les candidats et l'etat de reference.
// ---------------------------------------------------------------------------
void FdemSolver::armBrush() {
    computeFragments();                    // composantes connexes a cet instant
    brushFrag_ = fragId_;                  // fige : le classement se fera sur
    brushNFrag_ = nFrag_;                  // CES composantes, pas sur d'autres
    brushCand_.assign(el_.size(), 0);
    long nc = 0;
    for (std::size_t e = 0; e < el_.size(); ++e)
        if (fragId_[e] != 0) { brushCand_[e] = 1; ++nc; }
    brushU0_ = u_;                         // deplacement nodal de reference
    // vitesse initiale communiquee aux seuls noeuds des candidats
    std::vector<char> touched(X0_.size(), 0);
    for (std::size_t e = 0; e < el_.size(); ++e)
        if (brushCand_[e])
            for (int a = 0; a < 3; ++a) touched[el_[e].n[a]] = 1;
    // fragBrushZeroV : REMPLACER la vitesse au lieu de l'ajouter.
    //
    // L'article AJOUTE v0 a l'etat existant, ce qui suppose un etat au repos
    // (leur etape 1, « after completing the impact simulation »). Sur une
    // COUPE CONTINUE cet etat n'existe pas : mesure du 2026-08-18, apres
    // 244 us de repos outil arrete, la vitesse mediane des fragments est
    // encore de 5,4 m/s et la maximale ne decroit pas du tout — un fragment
    // libre qui vole n'a rien pour le freiner, et l'amortissement de Cundall
    // agit sur la deformation, pas sur un mouvement de corps rigide.
    //
    // Or l'echelle propre du balai est bornee par le critere « pas de
    // nouvelle fissure » : a 20 000 m/s2 sur un bloc de 20 mm, rho a h = 1,0
    // MPa contre ft = 10,6 MPa, et a t = 6 m/s. Les deux sont du MEME ORDRE
    // que le residuel : aucun choix de parametres ne les separe.
    //
    // Remettre les vitesses a zero rend la mesure exactement ce qu'elle doit
    // etre — la reponse a l'anti-gravite SEULE. C'est licite parce que la
    // phase de balayage est une CLASSIFICATION : rien de ce qui la suit n'est
    // de la physique, et le bilan B4 de l'impact est deja clos.
    long nn = 0;
    for (std::size_t i = 0; i < X0_.size(); ++i)
        if (touched[i]) {
            if (brushZeroV_) v_[i]  = brushV0_ * brushDir_;
            else             v_[i] += brushV0_ * brushDir_;
            ++nn;
        }
    // ---- CONTROLE DE REPOS (etape 1) -------------------------------------
    // L'echelle de vitesse du balai vaut v0 + a t_balayage : c'est la vitesse
    // qu'il communique lui-meme. Toute vitesse residuelle du meme ordre noie
    // le signal, et le classement mesure alors le mouvement d'AVANT, pas la
    // reponse a l'anti-gravite. On l'imprime pour que le journal tranche.
    double vres = 0.0, vresAll = 0.0;
    for (std::size_t i = 0; i < X0_.size(); ++i) {
        double vi = v_[i].norm();
        vresAll = std::max(vresAll, vi);
        if (touched[i]) vres = std::max(vres, vi);
    }
    brushArmed_ = true;
    brushT0_ = t_;
    std::cout << "[FDEM] BALAI arme a t = " << t_ << " s : " << nc << " / "
              << el_.size() << " elements candidats (" << brushNFrag_ - 1
              << " fragments hors corps principal), " << nn << " noeuds\n"
              << "[FDEM]   v0 = " << brushV0_ << " m/s, a = " << brushA_
              << " m/s2, direction (" << brushDir_.x() << ", " << brushDir_.y()
              << "), beta = " << brushBeta_ << "\n"
              << "[FDEM]   ATTENTION : phase de CLASSIFICATION, pas de "
                 "physique. Son travail est compte a part (brushWork_) et "
                 "n'entre pas dans le bilan B4.\n"
              << "[FDEM]   REPOS (etape 1) : vitesse residuelle max des "
                 "candidats " << vres << " m/s, du bloc entier " << vresAll
              << " m/s\n";
    // l'echelle propre du balai n'est connue qu'avec la duree de balayage ;
    // on la borne ici par la vitesse acquise en un temps caracteristique egal
    // a v0/a majore : on se contente de comparer a v0 + a * (duree deja
    // ecoulee est nulle), donc on signale par rapport a v0 seul, plus la
    // consigne d'usage : vres doit rester tres petit devant a * t_balayage.
    if (vres > 100.0 * brushV0_)
        std::cout << "[FDEM]   >>> ETAPE 1 VIOLEE : les candidats bougent "
                     "encore a " << vres << " m/s, soit " << vres / brushV0_
                  << " fois la vitesse initiale du balai. Le classement "
                     "mesurera le mouvement RESIDUEL, pas la reponse a "
                     "l'anti-gravite. Allonger le repos entre toolStop et "
                     "fragBrushStart.\n";
}

// ---------------------------------------------------------------------------
// Verdict : deplacement de chaque fragment contre celui d'une particule LIBRE.
// ---------------------------------------------------------------------------
void FdemSolver::brushReport() {
    if (!brushArmed_) return;
    double tb = t_ - brushT0_;
    double dref = brushV0_ * tb + 0.5 * brushA_ * tb * tb;
    double seuil = brushBeta_ * dref;
    // deplacement moyen par fragment, mesure DEPUIS l'armement
    std::vector<double> sum(brushNFrag_, 0.0), vol(brushNFrag_, 0.0);
    std::vector<long> cnt(brushNFrag_, 0);
    for (std::size_t e = 0; e < el_.size(); ++e) {
        int f = brushFrag_[e];
        double d = 0.0;
        for (int a = 0; a < 3; ++a)
            d += (u_[el_[e].n[a]] - brushU0_[el_[e].n[a]]).norm() / 3.0;
        sum[f] += d; ++cnt[f];
        vol[f] += el_[e].A0 * thk_;
    }
    double vLibre = 0.0, vCoince = 0.0;
    long nLibre = 0, nCoince = 0;
    for (int f = 1; f < brushNFrag_; ++f) {          // 0 = corps principal
        if (cnt[f] == 0) continue;
        if (sum[f] / cnt[f] > seuil) { vLibre += vol[f]; ++nLibre; }
        else                         { vCoince += vol[f]; ++nCoince; }
    }
    double vTot = vLibre + vCoince;
    std::cout << "[FDEM] BALAI, verdict apres " << tb << " s :\n"
              << "[FDEM]   particule libre de reference : d_ref = " << dref
              << " m, seuil = beta d_ref = " << seuil << " m\n"
              << "[FDEM]   fragments LIBRES  : " << nLibre << " (" << vLibre
              << " m^3/m)\n"
              << "[FDEM]   fragments COINCES : " << nCoince << " (" << vCoince
              << " m^3/m)\n"
              << "[FDEM]   volume detache CORRIGE " << vLibre << " contre "
              << vTot << " par le critere naif : "
              << (vTot > 0 ? 100.0 * (vTot - vLibre) / vTot : 0.0)
              << " % de SUREVALUATION\n"
              << "[FDEM]   travail du balai : " << brushWork_
              << " J/m (hors bilan B4)\n";
}

// ---------------------------------------------------------------------------
// Cohesive joints, 2-point (node-pair) integration. Per point:
//   delta   = u(B side) - u(A side); n = outward normal of A (from current
//             edge midpoints), e = edge tangent (A's CCW direction).
//   Opening dn = delta.n; tangential dtg = delta.e.
// Mode I (dn > 0):   sigma = min( (1-D) p dn, envelope(dn) ), envelope
//   linear ft -> 0 over [dnE, dnF]; the cap binding updates D.
// jointSoftening = yan replaces that envelope by the exponential reduction
//   factor f(D) ft of Yan et al. (IJRMMS 169, 2023) eq. 11, with D the
//   mixed-mode driver of eq. 16 and the origin-secant unloading of eq. 17.
// Compression:       sigma = p dn (+ small dashpot), no damage.
// Mode II:           trial tau = p (dtg - slip), capped by
//   (1-D) c + tanPhi * max(0, -sigma); on the cap, slip flows (return
//   mapping) and D grows with |slip| / slipF. Friction survives at D = 1.
// Forces: traction integrated with the trapezoid rule (L/2 per point),
//   equal and opposite on the two faces.
// A fully broken joint whose faces have clearly separated or slid by half an
// edge length goes DEAD: it stops carrying anything and its faces are handed
// to the general contact.
// ---------------------------------------------------------------------------
void FdemSolver::jointForces() {
    static const bool noTau = std::getenv("RKM_NOTAU") != nullptr;

    // Per-joint state updates (D, slip, tBreak) are private to the joint, so
    // the loop parallelizes over joints — but the force scatter hits nodes
    // shared by the up-to-three joints of an element, so each thread
    // accumulates into its own buffer (touched-list + per-thread seen flags)
    // and the buffers are reduced SERIALLY IN THREAD ORDER: deterministic
    // for a fixed thread count.
    auto processJoint = [&](Joint& J, auto&& addF, long& nb, long& nd,
                            double& dampW, double& jw) {
        if (J.dead || J.bonded) return;    // bonded: node binding carries it

        Eigen::Vector2d pA1 = X0_[J.a1] + u_[J.a1], pA2 = X0_[J.a2] + u_[J.a2];
        Eigen::Vector2d pB1 = X0_[J.b1] + u_[J.b1], pB2 = X0_[J.b2] + u_[J.b2];
        Eigen::Vector2d P = 0.5 * (pA1 + pB1), Q = 0.5 * (pA2 + pB2);
        Eigen::Vector2d ed = Q - P;
        double L = ed.norm();
        if (L < 1e-14) return;
        Eigen::Vector2d e = ed / L;
        Eigen::Vector2d n(e.y(), -e.x());              // outward from A

        double Ltrib = 0.5 * J.L0 * thk_;
        double dnMax = -1e30;
        double rsMaxO = 0.0;               // moteur de mode II du pas courant
                                           // (jointShearUnload = origin)

        const int ia[2] = {J.a1, J.a2};
        const int ib[2] = {J.b1, J.b2};
        for (int k = 0; k < 2; ++k) {
            Eigen::Vector2d delta = (X0_[ib[k]] + u_[ib[k]])
                                  - (X0_[ia[k]] + u_[ia[k]]);
            // dn0: adaptive-insertion opening offset (0 for intrinsic joints).
            // A joint born under tension starts AT the envelope peak, so the
            // traction it hands the fresh crack faces equals the traction the
            // bonded edge was transmitting — no release spike at insertion.
            double dn = delta.dot(n) + J.dn0;
            double dtg = delta.dot(e);
            dnMax = std::max(dnMax, dn);

            // ---- eq. 18 (jointShearUnload = origin seulement) ---------------
            // s_max, le plus grand glissement JAMAIS atteint, est mis a jour
            // ICI — avant la traction normale — parce qu'en mode `origin` le
            // moteur de mode II de l'eq. 16 en depend, et que ce moteur
            // alimente D, donc f(D), donc l'enveloppe NORMALE evaluee juste
            // en dessous.
            //   * sEff : le glissement compte depuis l'origine STAMPEE. En
            //     intrinseque J.slip vaut 0 et sEff = dtg ; en insertion
            //     adaptative activateJoint() y a ecrit -tau0/pj, si bien que
            //     le joint naissant transmet tau0 a dtg = 0 (continuite en
            //     cisaillement, symetrique de dn0 en mode I). En mode `origin`
            //     J.slip n'evolue plus : ce n'est plus un glissement plastique
            //     mais l'origine figee de la secante.
            //   * sE : le s_p de Munjiza, glissement au PIC = penalite contre
            //     l'enveloppe de Mohr-Coulomb NON endommagee, exact symetrique
            //     de dnE = ft/pj en mode I. La contrainte normale y est prise
            //     sur sa part geometrique pj*dn (regle (1) de la loi : un
            //     terme visqueux ne fixe jamais une resistance) — ce qui rend
            //     le calcul NON CIRCULAIRE, sE ne dependant pas de D.
            //     Sous confinement la branche elastique s'elargit, exactement
            //     comme le cap : l'endommagement de mode II demarre quand
            //     l'enveloppe de Mohr-Coulomb est atteinte, c'est-a-dire au
            //     meme instant que dans le retour radial — d'ou l'accord des
            //     deux modes en charge monotone.
            double smx = 0.0, sEff = 0.0, rsO = 0.0;
            if (shearOrigin_) {
                sEff = dtg - J.slip[k];
                double sm = std::abs(sEff);
                if (sm > J.smax[k]) J.smax[k] = sm;
                smx = J.smax[k];
                double sE = (J.coh + J.tanPhi
                             * rockim::mcFrictionTerm(J.pj * dn, J.ft,
                                                      yangEnv_)) / J.pj;
                if (sE < 0.0) sE = 0.0;
                double den = J.slipF;
                if (shearRangeCoulomb_) {
                    // Y3Dfd.c l. 1110-1126 : st = max(2 sp, 3 GfII/dpefs),
                    // dpefs = c + tan(phi)|sigma_n| en compression (en
                    // traction dpefs = c, leur clamp), reevalue a chaque pas.
                    // Ici : meme convention kI que J.slipF, donc on scale par
                    // c/fs au lieu de recalculer, et le plancher est 2 sE.
                    double fs = J.coh
                              + J.tanPhi * std::max(0.0, -(J.pj * dn));
                    if (fs > J.coh)
                        den = std::max(2.0 * sE, J.slipF * (J.coh / fs));
                }
                rsO = (den > 0.0 && smx > sE) ? (smx - sE) / den : 0.0;
                rsMaxO = std::max(rsMaxO, rsO);
            }

            // ================= normal traction =========================
            // Structured as the established codes do (Kratos DEM_KDEM, Yade
            // ViscoelasticPM/HertzMindlin, LAMMPS granular), after measuring
            // that this solver's original formulation — a one-sided dashpot
            // with a constant coefficient whose VISCOUS TERM ALONE was
            // clipped — has no equivalent anywhere. Two rules carry the fix:
            //
            //  (1) The failure criterion is evaluated on the ELASTIC part
            //      only. A viscous term must never be able to break a joint.
            //      (Kratos builds sigma on the elastic part and calls
            //      CalculateNormalForces BEFORE CalculateViscoDamping.)
            //  (2) An INTACT cohesive joint is BILATERAL: it is a bond, it is
            //      meant to pull, and clipping it at zero is what turned the
            //      dashpot into a rectifier — pushing while the faces close,
            //      resisting nothing while they separate, so every vibration
            //      cycle injected momentum. Clipping belongs to a BROKEN
            //      joint acting as a contact, and then it applies to the
            //      RESULTANT, with the excess handed back to the viscous part
            //      so the dissipation ledger stays exact (Yade HertzMindlin
            //      does `Fn = 0; Fvisc += normTemp`). LAMMPS goes further and
            //      makes the clip flag a FATAL ERROR on a cohesive normal
            //      model — the rule this code was missing.
            //
            // Measured before/after, under STRICTLY ZERO load on the
            // unstructured intra-grain mesh: 158 broken joints and a 57.5 MPa
            // phantom grip reaction, against 0 joints and 1.9e-10 MPa.
            double sigEl;                              // elastic / cohesive
            if (yanSoft_) {
                // ---- exponential softening of Yan et al., eq. 9/11-13/16-17 --
                // The article inserts an EXTRINSIC element that carries ft at
                // zero opening; rockim keeps an elastic branch of width
                // dnE = ft/pj (needed by the intrinsic mode, and the adaptive
                // insertion offset dn0 lands a tension-born joint exactly at
                // dnE, i.e. at the peak). The article's opening o is therefore
                // the opening BEYOND the elastic branch, and the critical
                // opening ot = dnF - dnE is the one calibrated on GfI.
                double ot = J.dnF - J.dnE;
                double rn = (ot > 0.0 && dn > J.dnE) ? (dn - J.dnE) / ot : 0.0;
                // rs : moteur de mode II. `plastic` le lit sur le glissement
                // PLASTIQUE cumule par le retour radial ; `origin` sur
                // (s_max - s_p)/(s_t - s_p), la forme litterale de l'eq. 14.
                double rs = shearOrigin_
                    ? rsO
                    : ((J.slipF > 0.0) ? std::abs(J.slip[k]) / J.slipF : 0.0);
                // eq. 16, which degenerates into eq. 12 (pure tension) and
                // eq. 14 (pure shear) when the other driver vanishes
                double Dnow = std::sqrt(rn * rn + rs * rs);
                // irreversibility (Fukuda et al. rule quoted by the article):
                // J.D IS Dmax, it never decreases
                if (Dnow > J.D) J.D = std::min(1.0, Dnow);
                double fdY = yan::fD(J.D, yanP_);
                if (dn >= 0.0) {
                    // envelope: min(elastic branch, f(D) ft) — the min is what
                    // makes shear damage cut the tensile strength too, and it
                    // is continuous at dn = dnE where pj dnE = ft
                    if (dn > J.omax[k]) J.omax[k] = dn;
                    double om = J.omax[k];
                    double sMax = std::min(J.pj * om, fdY * J.ft);
                    // eq. 17: below omax the element unloads/reloads on the
                    // secant to the origin. At dn = omax this returns sMax, so
                    // loading and unloading agree on the envelope.
                    sigEl = (om > 1e-30) ? sMax * dn / om : 0.0;
                } else {
                    // k- = k+(D) en mode adaptatif : la MEME secante des deux
                    // cotes de dn = 0, donc plus de saut de raideur.
                    // (EPFL arXiv:2511.14323 sec. 4 ; voir FdemSolver.hpp)
                    sigEl = (jcAdaptive_ ? (1.0 - J.D) * J.pj : J.pj) * dn;
                }
            } else if (dn >= 0.0) {
                double env = (dn <= J.dnE) ? J.pj * dn
                           : (dn >= J.dnF) ? 0.0
                           : J.ft * (J.dnF - dn) / (J.dnF - J.dnE);
                double tr = (1.0 - J.D) * J.pj * dn;
                if (tr > env) {
                    sigEl = env;
                    if (dn > 1e-30) {
                        double Dn = 1.0 - env / (J.pj * dn);
                        if (Dn > J.D) J.D = std::min(1.0, Dn);
                    }
                } else sigEl = tr;
            } else {
                // idem : k- = k+(D) = (1-D) pj (EPFL arXiv:2511.14323 sec. 4)
                sigEl = (jcAdaptive_ ? (1.0 - J.D) * J.pj : J.pj) * dn;
            }

            double sig = sigEl;
            if (xiJ_ > 0.0) {
                Eigen::Vector2d vrel = v_[ib[k]] - v_[ia[k]];
                double meff = 0.5 * std::min(m_[ia[k]], m_[ib[k]]);
                double cd = 2.0 * xiJ_ * std::sqrt(J.pj * Ltrib * meff);
                // hard bound: past m_eff/dt the dashpot reverses the approach
                // velocity within one step instead of damping it (MOOSE)
                cd = std::min(cd, meff / dt_);
                // SIGN. The traction acts on B as -sig*Ltrib*n (it pulls B
                // back toward A), so a term that OPPOSES the opening rate
                // vrel.n must carry a PLUS sign. The original code had
                // `- cd * vrel.dot(n)`, whose power is +cd (vrel.n)^2: the
                // joint dashpot has been ANTI-damping since the code was
                // written. The one-sided clip hid it in tension and the
                // fade-out patch hid it near closure; making the bond
                // bilateral exposed it at once — every joint of the specimen
                // broke under zero load.
                double sigV = cd * vrel.dot(n) / Ltrib;
                if (J.D < 1.0) {
                    sig = sigEl + sigV;                // bilateral bond
                } else if (dn < 0.0) {
                    sig = std::min(0.0, sigEl + sigV); // contact: clip the SUM
                }
                // power of the viscous part: F_B.v_B + F_A.v_A
                dampW -= (sig - sigEl) * Ltrib * vrel.dot(n) * dt_;
            }

            // --- tangential traction (damage-plastic, frictional) ---
            // yan: eq. 10, the cohesion is scaled by f(D) instead of (1 - D).
            // The Coulomb term is left unscaled by default so a crushed joint
            // keeps residual friction (jointFrictionScaled = 1 for the literal
            // eq. 10, where f(D) multiplies the whole cap).
            double fdS = yanSoft_ ? yan::fD(J.D, yanP_) : 0.0;
            double coh = yanSoft_ ? fdS * J.coh : (1.0 - J.D) * J.coh;
            double muS = (yanSoft_ && yanFricScaled_) ? fdS : 1.0;
            double tauLim = coh + muS * J.tanPhi
                          * rockim::mcFrictionTerm(sig, J.ft, yangEnv_);
            if (tauLim < 0.0) tauLim = 0.0;
            double tau;
            if (shearOrigin_) {
                // ---- eq. 18 : sécante à l'origine ------------------------
                // Symetrique EXACT de l'eq. 17 ecrite plus haut pour le mode
                // I : enveloppe = min(branche elastique, cap), evaluee AU
                // glissement maximal s_max ; l'etat courant est lu sur la
                // droite qui joint l'origine a ce point. En charge monotone
                // |dtg| = s_max et l'on retombe sur l'enveloppe ; en decharge
                // la traction retourne a zero avec le glissement, sans
                // conserver de glissement plastique.
                double tauEnv = std::min(J.pj * smx, tauLim);
                tau = (noTau || smx <= 1e-30) ? 0.0 : tauEnv * sEff / smx;
                // en `yan`, D a deja ete mis a jour par l'eq. 16 au-dessus
                if (!yanSoft_ && rsO > J.D) J.D = std::min(1.0, rsO);
            } else {
                // ---- retour radial (defaut, inchange) --------------------
                // Plasticite a retour, comme dans la loi lineaire : sur le cap
                // cela reproduit tau = f(D) c de l'eq. 10 en glissement
                // monotone, et la decharge suit la secante de PENALITE plutot
                // que la secante a l'origine de l'eq. 18 — ecart assume, leve
                // par jointShearUnload = origin.
                double tauTr = noTau ? 0.0 : J.pj * (dtg - J.slip[k]);
                tau = std::clamp(tauTr, -tauLim, tauLim);
                if (tau != tauTr) {
                    J.slip[k] += (tauTr - tau) / J.pj;     // return mapping
                    double Dt;
                    if (yanSoft_) {
                        double ot = J.dnF - J.dnE;
                        double rn = (ot > 0.0 && dn > J.dnE) ? (dn - J.dnE) / ot
                                                             : 0.0;
                        double rs = std::abs(J.slip[k]) / J.slipF;
                        Dt = std::sqrt(rn * rn + rs * rs);   // eq. 16 / 14
                    } else {
                        Dt = std::abs(J.slip[k]) / J.slipF;
                    }
                    if (Dt > J.D) J.D = std::min(1.0, Dt);
                }
            }

            Eigen::Vector2d trac = (sig * n + tau * e) * Ltrib;
            addF(ib[k], -trac);                        // pull B back toward A
            addF(ia[k], trac);
            // V2/B4 : travail TOTAL des tractions (visqueux inclus, isole
            // dans dampW) — lecture pure des vitesses
            jw += trac.dot(v_[ia[k]] - v_[ib[k]]) * dt_;
        }

        if (J.D >= 1.0) {
            if (J.tBreak < 0) {
                J.tBreak = t_; ++nb;
                // partition of the eq. 16 driver at the breaking instant:
                // rn = normal (mode I) ratio, rs = sliding (mode II) ratio.
                // Recomputed here from the same state processJoint used above,
                // so no extra bookkeeping is carried per step.
                double otF = J.dnF - J.dnE;
                double rnF = (otF > 0.0 && dnMax > J.dnE)
                           ? (dnMax - J.dnE) / otF : 0.0;
                double rsF;
                if (shearOrigin_) {           // pas de glissement plastique :
                                              // le moteur est celui de l'eq. 14
                    rsF = rsMaxO;
                } else {
                    double sMx = std::max(std::abs(J.slip[0]),
                                          std::abs(J.slip[1]));
                    rsF = (J.slipF > 0.0) ? sMx / J.slipF : 0.0;
                }
                double den = rnF * rnF + rsF * rsF;
                J.failMode = (den > 1e-300) ? (rnF * rnF) / den : 1.0;
                J.rnB = rnF;
                J.rsB = rsF;
                J.bmode = (rnF >= rsF) ? 1 : 2;  // 1 traction, 2 cisaillement
            }
            // Release to general contact ONLY on clear separation. A crushed
            // joint sliding under compression must stay alive as the paired
            // frictional contact of its own faces: killing it by slip hands
            // interpenetrated faces to the general contact, whose penalty
            // then releases 1/2 k pen^2 of energy created from nothing (the
            // pump isolated by the net-work meter).
            // ++nd tire EXACTEMENT une fois par joint : au pas suivant le
            // garde `if (J.dead || J.bonded) return;` en tete de lambda coupe
            // court. C'est ce compteur qui estampille le rebuild en mode eager.
            if (dnMax > 3.0 * J.dnF) { J.dead = true; ++nd; }
        }
    };

#ifdef _OPENMP
    int nT = omp_get_max_threads();
    if (nT == 1) {                       // bit-identical to the serial build
        long nb1 = 0, nd1 = 0;
        double dw1 = 0.0, jw1 = 0.0;
        auto addF1 = [&](int i, const Eigen::Vector2d& v) { f_[i] += v; };
        for (auto& J : jt_) processJoint(J, addF1, nb1, nd1, dw1, jw1);
        nBroken_ += nb1;
        nDead_ += nd1;
        dampWork_ += dw1;
        jointWork_ += jw1;
        return;
    }
    if ((int)fTL_.size() != nT) {
        fTL_.assign(nT, std::vector<Eigen::Vector2d>(
                            X0_.size(), Eigen::Vector2d::Zero()));
        seenTL_.assign(nT, std::vector<char>(X0_.size(), 0));
        touchedTL_.assign(nT, {});
    }
    std::vector<long> nbT(nT, 0), ndT(nT, 0);
    std::vector<double> dwT(nT, 0.0), jwT(nT, 0.0);
#pragma omp parallel
    {
        int t = omp_get_thread_num();
        auto& fb = fTL_[t];
        auto& seen = seenTL_[t];
        auto& tl = touchedTL_[t];
        tl.clear();
        long nb = 0, nd = 0;
        double dw = 0.0, jw = 0.0;
        auto addF = [&](int i, const Eigen::Vector2d& v) {
            if (!seen[i]) { seen[i] = 1; tl.push_back(i); }
            fb[i] += v;
        };
#pragma omp for schedule(static)
        for (int jI = 0; jI < (int)jt_.size(); ++jI)
            processJoint(jt_[jI], addF, nb, nd, dw, jw);
        nbT[t] = nb;
        ndT[t] = nd;
        dwT[t] = dw;
        jwT[t] = jw;
    }
    for (int t = 0; t < nT; ++t) {                     // deterministic order
        for (int i : touchedTL_[t]) {
            f_[i] += fTL_[t][i];
            fTL_[t][i].setZero();
            seenTL_[t][i] = 0;
        }
        nBroken_ += nbT[t];
        nDead_ += ndT[t];
        dampWork_ += dwT[t];
        jointWork_ += jwT[t];
    }
#else
    long nb = 0, nd = 0;
    double dw = 0.0, jw = 0.0;
    auto addF = [&](int i, const Eigen::Vector2d& v) { f_[i] += v; };
    for (auto& J : jt_) processJoint(J, addF, nb, nd, dw, jw);
    nBroken_ += nb;
    nDead_ += nd;
    dampWork_ += dw;
    jointWork_ += jw;
#endif
}

// ---------------------------------------------------------------------------
// General node-edge contact between released faces: exterior faces plus the
// faces of every broken (D >= 1) joint. A cell grid over active edges is
// rebuilt whenever the broken count changes; nodes of active edges are tested
// against nearby edges of OTHER elements (skipping the partner across a
// still-live joint, whose interaction the joint itself handles). Penetration
// = signed distance along the edge's outward normal, capped at half the mesh
// size; penalty + tanh-regularised Coulomb friction, reaction distributed to
// the edge nodes by the projection weights.
// ---------------------------------------------------------------------------
void FdemSolver::rebuildContactEdges() {
    if (!poolBuilt_) {
        // gcXwindow: on the SHPB assembly all but a few hundred of the
        // ~40 000 exterior faces are bar flanks metres away from anything
        // they could ever touch. Restricting the ACTIVE set to a window
        // around the specimen is what makes the contact sweep cost
        // O(interface) instead of O(assembly). The window tests X0_ (initial
        // positions), so the pool is CONSTANT: built once, kept.
        pool_ = exterior_;
        if (gcXmin_ > -1e299) {
            std::vector<BEdge> keep;
            for (const auto& be : pool_) {
                double xm = 0.5 * (X0_[be.na].x() + X0_[be.nb].x());
                if (xm >= gcXmin_ && xm <= gcXmax_) keep.push_back(be);
            }
            pool_.swap(keep);
        }
        extOn_.assign(pool_.size(), gcAdaptive_ ? 0 : 1);
        poolBuilt_ = true;
    }
    act_.clear();
    actPool_.clear();                      // act idx -> pool idx (-1 liberee)
    if (gcAdaptive_) {                    // sous-ensemble ACTIF, ordre du pool
        for (std::size_t k = 0; k < pool_.size(); ++k)
            if (extOn_[k]) {
                act_.push_back(pool_[k]);
                actPool_.push_back((int)k);
            }
    } else {
        act_.insert(act_.end(), pool_.begin(), pool_.end());
        for (std::size_t k = 0; k < pool_.size(); ++k)
            actPool_.push_back((int)k);
    }
    // Les faces liberees viennent du CACHE deadList_, rafraichi UNIQUEMENT
    // par le declencheur historique (nBroken_ change, pas % 8) dans
    // generalContact(). Un joint peut mourir SANS casser (separation tardive,
    // dnMax > 3 dnF longtemps apres D = 1) : le mode full n'integre alors ses
    // faces qu'au rebuild suivant la casse SUIVANTE. Les recompositions
    // declenchees par l'activation adaptative reutilisent le cache tel quel,
    // sinon elles avanceraient cette entree de quelques pas et les deux modes
    // divergeraient au bit pres (mesure : c'etait LA source de l'ecart).
    act_.insert(act_.end(), deadList_.begin(), deadList_.end());
    actPool_.insert(actPool_.end(), deadList_.size(), -1);
    haveDead_ = !deadList_.empty();
    if (poolTouch_.empty()) poolTouch_.assign(pool_.size(), -1);
    std::vector<char> inAct(X0_.size(), 0);
    for (const auto& be : act_) inAct[be.na] = inAct[be.nb] = 1;
    actNodes_.clear();
    for (int i = 0; i < (int)X0_.size(); ++i)
        if (inAct[i]) actNodes_.push_back(i);
}

// ---------------------------------------------------------------------------
// gcActivation = adaptive — le balayage d'activation. Regles C / A / B de
// l'en-tete ; monotone ; serie (deterministe). Voir FdemSolver.hpp pour la
// justification de chaque regle et de la cadence.
// ---------------------------------------------------------------------------
void FdemSolver::activationSweep() {
    // ---- (1) composantes connexes et peau endommagee — seulement quand la
    // topologie a change (un joint ne change de statut qu'en cassant, et
    // nBroken_ est incremente a l'instant meme ou D atteint 1)
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
        // regle C etendue d'un ANNEAU : tout element partageant un SOMMET
        // avec un element au bord casse est peau endommagee aussi. C'est le
        // coin du bourrelet de cratere : son element n'a encore rien casse,
        // mais la surface s'y deforme deja et son premier recouvrement peut
        // preceder tout contact porte — mesure sur la percussion 2D (les deux
        // faces du bord divergeaient au meme pas quelle que soit la marge).
        if (vElems_.empty()) {             // topologie STATIQUE : bati une fois
            int nV = 0;
            for (int v : vOf_) nV = std::max(nV, v + 1);
            vElems_.assign(nV, {});
            for (int nd = 0; nd < (int)vOf_.size(); ++nd)
                vElems_[vOf_[nd]].push_back(elemOf_[nd]);
        }
        std::vector<char> ring(el_.size(), 0);
        for (int e2 = 0; e2 < (int)el_.size(); ++e2)
            if (elemDam_[e2])
                for (int a = 0; a < 3; ++a)
                    for (int e3 : vElems_[vOf_[el_[e2].n[a]]]) ring[e3] = 1;
        for (int e2 = 0; e2 < (int)el_.size(); ++e2)
            if (ring[e2]) elemDam_[e2] = 1;
        bodyStamp_ = nBroken_;
    }
    // ---- (2) un seul corps, rien d'active, rien de casse : rien ne peut
    // toucher (l'approximation assumee : un continuum intact ne se replie
    // pas sur lui-meme). Le cas ecrasamment majoritaire d'un run de
    // traction/UCS avant le pic — et TOUT l'avant-endommagement en percussion.
    if (nBodies_ <= 1 && nActivated_ == 0 && !haveDead_) {
        nextSweep_ = stepCount_ + gcActEvery_;
        return;
    }
    // ---- (3) sources et cibles ------------------------------------------
    double cl = gcCell_ > 0.0 ? gcCell_ : 2.0 * hmin_; // = cell_ de detect()
    double M = gcActMargin_ * cl;
    struct SFace {
        Eigen::Vector2d c;
        double r;
        int body, idx;                     // idx dans le pool, -1 = liberee
        char on, touched;
    };
    static std::vector<SFace> sf;
    sf.clear();
    auto push = [&](int na, int nb, int elem, char on, int idx) {
        Eigen::Vector2d P = X0_[na] + u_[na], Q = X0_[nb] + u_[nb];
        char tch = on && (lastTouch_[na] >= 0 || lastTouch_[nb] >= 0);
        sf.push_back({0.5 * (P + Q), 0.5 * (Q - P).norm(),
                      bodyOf_[elem], idx, on, tch});
    };
    for (std::size_t k = 0; k < pool_.size(); ++k)
        push(pool_[k].na, pool_[k].nb, pool_[k].elem, extOn_[k], (int)k);
    for (const auto& J : jt_)
        if (J.dead) {
            push(J.a1, J.a2, J.eA, 1, -1);
            push(J.b2, J.b1, J.eB, 1, -1);
        }
    double rmax = 0.0;
    for (const auto& s : sf) rmax = std::max(rmax, s.r);
    // ---- (4) grille de hachage, pas = M + 2 rmax : le pochoir 3x3 couvre
    // toute paire dont l'ecart des spheres englobantes est < M
    double cs = M + 2.0 * rmax + 1e-300;
    static std::unordered_map<uint64_t, std::vector<int>> hg;
    hg.clear();
    auto key = [&](const Eigen::Vector2d& p) {
        long long ix = (long long)std::floor(p.x() / cs);
        long long iy = (long long)std::floor(p.y() / cs);
        return (uint64_t(uint32_t(ix)) | (uint64_t(uint32_t(iy)) << 32));
    };
    auto keyIJ = [&](long long ix, long long iy) {
        return (uint64_t(uint32_t(ix)) | (uint64_t(uint32_t(iy)) << 32));
    };
    for (int q = 0; q < (int)sf.size(); ++q) hg[key(sf[q].c)].push_back(q);
    bool changed = false;
    for (int q = 0; q < (int)sf.size(); ++q) {
        SFace& f = sf[q];
        if (f.on || f.idx < 0) continue;               // cibles : pool inactif
        if (elemDam_[pool_[f.idx].elem]) {             // regle C
            extOn_[f.idx] = 1;
            ++nActivated_;
            changed = true;
            continue;
        }
        long long ix = (long long)std::floor(f.c.x() / cs);
        long long iy = (long long)std::floor(f.c.y() / cs);
        bool hit = false;
        for (long long dy = -1; dy <= 1 && !hit; ++dy)
            for (long long dx = -1; dx <= 1 && !hit; ++dx) {
                auto it = hg.find(keyIJ(ix + dx, iy + dy));
                if (it == hg.end()) continue;
                for (int g : it->second) {
                    if (g == q) continue;
                    const SFace& s = sf[g];
                    // regle A : autre corps. regle B : deja en contact.
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
    if (changed) {
        rebuildContactEdges();             // recompose (cache mort INCHANGE)
        if (actStep_.empty()) actStep_.assign(pool_.size(), -1);
        for (std::size_t k = 0; k < pool_.size(); ++k)
            if (extOn_[k] && actStep_[k] < 0) actStep_[k] = stepCount_;
    }
    // ---- (5) cadence : entre deux balayages, l'approche par pas est bornee
    // par 2 v_max dt ; on la garde sous la demi-marge (securite 2)
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
// contact = potential — A3. Paires d'ELEMENTS candidates depuis le jeu actif
// (compose donc avec gcActivation et gcXwindow), detection O(N) type NBS
// (binning AABB en grille de hachage, paires uniques par stamp), exclusion
// des paires liees par un joint VIVANT (le joint porte leur interaction),
// puis la force distribuee de Munjiza calculee par PotentialContact.hpp et
// le frottement incremental de l'eq. 4-5. Serie et deterministe.
// ---------------------------------------------------------------------------
void FdemSolver::potentialContact() {
    // table (eMin,eMax) -> joint, batie UNE fois : les paires sont figees a
    // la construction du maillage, seul l'etat dead evolue
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
    for (const auto& be : act_)
        if (emark[be.elem] != epoch) {
            emark[be.elem] = epoch;
            elems.push_back(be.elem);
        }
    if (elems.size() < 2) return;

    // ---- (2) AABB courantes + grille DENSE a seaux REUTILISES (N1) -------
    // La grille de hachage reallouait chaque seau a chaque pas (mesure : la
    // detection dominait le cout du potentiel, x2,6 sur la percussion 3D).
    // Grille dense fenetree comme celle des faces : boite fixe autour du
    // domaine, seaux clear() sans liberation, elements hors boite ignores
    // (meme regle que le contact noeud-face — un debris parti n'est plus
    // resolu).
    double cl = gcCell_ > 0.0 ? gcCell_ : 2.0 * hmin_;
    Eigen::Vector2d ebLo(-0.5 * W_, -0.5 * H_);
    Eigen::Vector2d ebHi(1.5 * W_, 2.0 * H_);
    Eigen::Vector2d egMin = ebLo - Eigen::Vector2d(cl, cl);
    Eigen::Vector2d espan = ebHi - egMin + Eigen::Vector2d(cl, cl);
    int egx = std::max(1, int(espan.x() / cl) + 1);
    int egy = std::max(1, int(espan.y() / cl) + 1);
    static std::vector<Eigen::Vector2d> elo, ehi;
    elo.resize(elems.size());
    ehi.resize(elems.size());
    static std::vector<char> einb;
    einb.resize(elems.size());
    static std::vector<std::vector<int>> eg;
    {
        std::size_t nC = (std::size_t)egx * egy;
        if (eg.size() != nC) eg.assign(nC, {});
        else for (auto& c : eg) c.clear();
    }
    for (int q = 0; q < (int)elems.size(); ++q) {
        const Elem& E = el_[elems[q]];
        Eigen::Vector2d lo = X0_[E.n[0]] + u_[E.n[0]], hi = lo;
        for (int a = 1; a < 3; ++a) {
            Eigen::Vector2d p = X0_[E.n[a]] + u_[E.n[a]];
            lo = lo.cwiseMin(p);
            hi = hi.cwiseMax(p);
        }
        elo[q] = lo;
        ehi[q] = hi;
        einb[q] = (hi.x() > ebLo.x() && lo.x() < ebHi.x()
                   && hi.y() > ebLo.y() && lo.y() < ebHi.y());
        if (!einb[q]) continue;
        int x0 = std::clamp(int((lo.x() - egMin.x()) / cl), 0, egx - 1);
        int x1 = std::clamp(int((hi.x() - egMin.x()) / cl), 0, egx - 1);
        int y0 = std::clamp(int((lo.y() - egMin.y()) / cl), 0, egy - 1);
        int y1 = std::clamp(int((hi.y() - egMin.y()) / cl), 0, egy - 1);
        for (int cy = y0; cy <= y1; ++cy)
            for (int cx = x0; cx <= x1; ++cx)
                eg[(std::size_t)cy * egx + cx].push_back(q);
    }

    // ---- (3) paires candidates, en ordre CANONIQUE -----------------------
    // Les paires sont collectees puis TRIEES par (eLo, eHi) : l'ordre de
    // traitement — donc l'ordre des sommes de forces — devient independant
    // de l'implementation de la detection. Toute optimisation future de la
    // grille (ou un portage GPU) est alors bit-neutre par construction.
    static std::vector<int> pstamp;        // dedup de paire par q courant
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
        for (int cy = y0; cy <= y1; ++cy)
            for (int cx = x0; cx <= x1; ++cx)
                for (int r : eg[(std::size_t)cy * egx + cx]) {
                    if (r <= q || pstamp[r] == q) continue;
                    pstamp[r] = q;
                    if (elo[r].x() > ehi[q].x() || elo[q].x() > ehi[r].x()
                        || elo[r].y() > ehi[q].y() || elo[q].y() > ehi[r].y())
                        continue;
                    uint64_t a = (uint64_t)std::min(elems[q], elems[r]);
                    uint64_t b = (uint64_t)std::max(elems[q], elems[r]);
                    pairs.push_back((a << 32) | b);
                }
    }
    std::sort(pairs.begin(), pairs.end());

    for (uint64_t pk : pairs) {
        {
            int eLo = (int)(pk >> 32);
            int eHi = (int)(pk & 0xFFFFFFFFu);
                    // exclusion : joint VIVANT (bonded ou non mort) — le
                    // joint porte l'interaction de sa propre paire, contact
                    // du joint casse-mais-vivant compris. Un joint MORT rend
                    // la paire au contact general.
                    auto itJ = jointOfPair_.find(pk);
                    if (itJ != jointOfPair_.end() && !jt_[itJ->second].dead)
                        continue;
                    const Elem& EA = el_[eLo];
                    const Elem& EB = el_[eHi];
                    pot::V2 pa[3], pb[3];
                    for (int k = 0; k < 3; ++k) {
                        pa[k] = X0_[EA.n[k]] + u_[EA.n[k]];
                        pb[k] = X0_[EB.n[k]] + u_[EB.n[k]];
                    }
                    pot::PairForce R;
                    if (!pot::pairForce(pa, pb, potP_, R)) continue;
                    // ---- releve de naissance par AIRE (voir PotHist::aRef),
                    // constante de temps gcBirthTau (1 us), semantique pen0_
                    auto [itH, isNewH] = potFt_.try_emplace(pk);
                    PotHist& H = itH->second;
                    if (isNewH) H.aRef = R.area;
                    else H.aRef *= relax_;
                    double sc = std::max(0.0, 1.0 - H.aRef / R.area);
                    R.F *= sc;
                    for (int k = 0; k < 3; ++k) {
                        R.fA[k] *= sc;
                        R.fB[k] *= sc;
                    }
                    // ---- forces nodales normales + compteur de travail ---
                    double w = 0.0;
                    for (int k = 0; k < 3; ++k) {
                        f_[EA.n[k]] += R.fA[k];
                        f_[EB.n[k]] += R.fB[k];
                        w += R.fA[k].dot(v_[EA.n[k]])
                           + R.fB[k].dot(v_[EB.n[k]]);
                    }
                    gcWork_ += w * dt_;
                    // barycentriques du centroide du recouvrement : servent
                    // aux estampilles PONDEREES de la regle B et au frottement
                    pot::Bary bA, bB;
                    bA.set(pa[0], pa[1], pa[2]);
                    bB.set(pb[0], pb[1], pb[2]);
                    double la[3], lb[3];
                    bA.lam(R.cen, la);
                    bB.lam(R.cen, lb);
                    // estampilles regle B (A1) PONDEREES : seuls les noeuds
                    // qui FONT FACE au contact (poids barycentrique > 0,1)
                    // deviennent sources d'activation. L'estampillage de tous
                    // les noeuds des deux elements sur-propageait la regle B
                    // (mesure : 96 % des faces activees sur la percussion 3D
                    // longue, contre 4 % en penalite noeud-face).
                    if (gcAdaptive_)
                        for (int k = 0; k < 3; ++k) {
                            if (la[k] > 0.1) lastTouch_[EA.n[k]] = stepCount_;
                            if (lb[k] > 0.1) lastTouch_[EB.n[k]] = stepCount_;
                        }
                    // ---- amortissement normal OPT-IN (potXi) -------------
                    // Amortisseur visqueux sur la vitesse relative NORMALE au
                    // centroide du recouvrement, calibre sur la masse reduite
                    // des deux elements : c = 2 xi sqrt(potP_ mR). Ecrete du
                    // cote SEPARATION pour que la resultante normale ne
                    // devienne jamais tractive (deux fragments ne doivent pas
                    // pouvoir se coller). Voir FdemSolver.hpp pour le pourquoi.
                    if (potXi_ > 0.0) {
                        double FnD = R.F.norm();
                        if (FnD > 1e-300) {
                            Eigen::Vector2d nhD = R.F / FnD;
                            Eigen::Vector2d vAD =
                                la[0] * v_[EA.n[0]] + la[1] * v_[EA.n[1]]
                                + la[2] * v_[EA.n[2]];
                            Eigen::Vector2d vBD =
                                lb[0] * v_[EB.n[0]] + lb[1] * v_[EB.n[1]]
                                + lb[2] * v_[EB.n[2]];
                            Eigen::Vector2d vrelD = vAD - vBD;
                            double mAD = rhoP_[EA.phase] * EA.A0 * thk_;
                            double mBD = rhoP_[EB.phase] * EB.A0 * thk_;
                            double mRD = mAD * mBD / (mAD + mBD);
                            double cD = 2.0 * potXi_ * std::sqrt(potP_ * mRD);
                            double fdD = -cD * vrelD.dot(nhD);
                            if (fdD < -FnD) fdD = -FnD;   // jamais tractif
                            Eigen::Vector2d FdD = fdD * nhD;
                            for (int k = 0; k < 3; ++k) {
                                f_[EA.n[k]] += la[k] * FdD;
                                f_[EB.n[k]] -= lb[k] * FdD;
                            }
                            gcWork_ += FdD.dot(vrelD) * dt_;   // V2/B4
                        }
                    }
                    // ---- frottement incremental (eq. 4-5) ----------------
                    // ressort tangentiel a HISTOIRE : Ft += -kt (vrel.t) dt,
                    // plafonne au cap de Coulomb mu |Fn| (glissement). La
                    // direction normale est celle de la resultante du
                    // potentiel ; l'historique est tourne dans le plan
                    // tangent courant a chaque pas.
                    if (muC_ > 0.0 && potKt_ > 0.0) {
                        double Fn = R.F.norm();
                        if (Fn > 1e-300) {
                            Eigen::Vector2d vA =
                                la[0] * v_[EA.n[0]] + la[1] * v_[EA.n[1]]
                                + la[2] * v_[EA.n[2]];
                            Eigen::Vector2d vB =
                                lb[0] * v_[EB.n[0]] + lb[1] * v_[EB.n[1]]
                                + lb[2] * v_[EB.n[2]];
                            Eigen::Vector2d vrel = vA - vB;
                            Eigen::Vector2d nh = R.F / Fn;
                            Eigen::Vector2d th(-nh.y(), nh.x());
                            if (H.step < stepCount_ - 1) H.Ft.setZero();
                            Eigen::Vector2d Ft =
                                H.Ft - H.Ft.dot(nh) * nh;      // rotation
                            Ft -= potKt_ * (vrel.dot(th) * dt_) * th;
                            double cap = ctcMu(eLo, eHi) * Fn;  // WP6
                            double Ftn = Ft.norm();
                            if (Ftn > cap && Ftn > 0.0)
                                Ft *= cap / Ftn;               // glissement
                            H.Ft = Ft;
                            H.step = stepCount_;
                            for (int k = 0; k < 3; ++k) {
                                f_[EA.n[k]] += la[k] * Ft;
                                f_[EB.n[k]] -= lb[k] * Ft;
                            }
                            gcWork_ += Ft.dot(vrel) * dt_;
                            gcFricWork_ += Ft.dot(vrel) * dt_;  // V2/B4
                        }
                    }
        }
    }
}

void FdemSolver::generalContact() {
    if (std::getenv("RKM_NOGC")) return;               // bisection switch
    if (gcAdaptive_) {
        if (!poolBuilt_) rebuildContactEdges();        // le pool avant tout
        if (lastTouch_.empty()) lastTouch_.assign(X0_.size(), -1);
        // le balayage se declenche a sa cadence propre ET des que nBroken_
        // change, aligne sur la cadence % 8 du rebuild : les faces liberees
        // par les joints frais morts deviennent sources au meme pas que leur
        // entree dans act_
        if (stepCount_ >= nextSweep_
            || (sweepBroken_ != nBroken_
                && (sweepBroken_ < 0 || stepCount_ % 8 == 0))) {
            activationSweep();             // peut recomposer act_ (cache mort
            sweepBroken_ = nBroken_;       // inchange : timing = mode full)
        }
    }
    // A' : en eager l'estampille suit les SEPARATIONS (nDead_) et non les
    // casses, et la grille % 8 saute — les faces entrent au pas ou leur joint
    // meurt. En legacy, comportement historique inchange au bit pres.
    const long gcStamp = gcEager_ ? nDead_ : nBroken_;
    if (actStamp_ != gcStamp
        && (actStamp_ < 0 || gcEager_ || stepCount_ % 8 == 0)) {
        deadList_.clear();                 // LE declencheur historique : seul
        for (const auto& J : jt_)          // endroit ou le cache se rafraichit
            if (J.dead) {
                deadList_.push_back({J.eA, J.a1, J.a2});
                deadList_.push_back({J.eB, J.b2, J.b1});  // CCW oppose dans B
            }
        rebuildContactEdges();
        actStamp_ = gcStamp;
    }
    if (act_.empty()) return;
    if (contactPot_) {                     // A3 : contact par potentiel —
        potentialContact();                // meme jeu actif, autre physique
        return;
    }

    // Grid over current edge midpoints, CLIPPED to a fixed box around the
    // domain: without the clip, one fast ejected fragment stretches the
    // bounding box and the per-step grid allocation explodes (this was a
    // measured 20x slowdown). Debris outside the box no longer needs
    // resolved contact and is simply skipped.
    cell_ = gcCell_ > 0.0 ? gcCell_ : 2.0 * hmin_;
    // Default box: the domain blown up by half its size, which is what a
    // percussion crater needs. On the 3.55 m x 0.05 m SHPB assembly that box is
    // 5.3 m x 0.15 m and, binned at 2 hmin = 1.5 mm, allocates 3.5e5 cells
    // EVERY step for 3 contacting faces. gcBoxMesh keeps the box tight around
    // the mesh instead (the bars cannot go anywhere: their far ends are driven
    // and damped).
    Eigen::Vector2d boxLo(-0.5 * W_, -0.5 * H_);
    Eigen::Vector2d boxHi(1.5 * W_, 2.0 * H_);
    if (gcBoxMesh_) {
        // AABB of the ACTIVE faces (already restricted to the gcXwindow), not
        // of the whole assembly: the grid then has O(window/cell) buckets
        // instead of O(assembly/cell), and every node's 3x3 stencil holds a
        // handful of candidates rather than the entire interface.
        Eigen::Vector2d lo(1e300, 1e300), hi(-1e300, -1e300);
        for (const auto& be : act_)
            for (int nid : {be.na, be.nb}) {
                Eigen::Vector2d q = X0_[nid] + u_[nid];
                lo = lo.cwiseMin(q);
                hi = hi.cwiseMax(q);
            }
        if (hi.x() > lo.x()) {
            Eigen::Vector2d pad(4.0 * cell_, 4.0 * cell_);
            boxLo = lo - pad;
            boxHi = hi + pad;
        }
    }
    std::vector<Eigen::Vector2d> mid(act_.size());
    std::vector<char> inBox(act_.size(), 0);
    for (std::size_t k = 0; k < act_.size(); ++k) {
        Eigen::Vector2d P = X0_[act_[k].na] + u_[act_[k].na];
        Eigen::Vector2d Q = X0_[act_[k].nb] + u_[act_[k].nb];
        mid[k] = 0.5 * (P + Q);
        inBox[k] = (mid[k].x() > boxLo.x() && mid[k].x() < boxHi.x()
                    && mid[k].y() > boxLo.y() && mid[k].y() < boxHi.y());
    }
    gmin_ = boxLo - Eigen::Vector2d(cell_, cell_);
    Eigen::Vector2d span = boxHi - gmin_ + Eigen::Vector2d(cell_, cell_);
    gx_ = std::max(1, int(span.x() / cell_) + 1);
    gy_ = std::max(1, int(span.y() / cell_) + 1);
    // reuse the buckets instead of destroying them: assign() frees every
    // inner vector every step, clear() keeps their capacity (bit-neutral)
    {
        std::size_t nCells = (std::size_t)gx_ * gy_;
        if (grid_.size() != nCells) grid_.assign(nCells, {});
        else for (auto& c : grid_) c.clear();
    }
    for (std::size_t k = 0; k < act_.size(); ++k) {
        if (!inBox[k]) continue;
        // Bin the edge into EVERY cell its AABB covers, not only the
        // midpoint cell: voronoi faces can be several hmin long, and
        // midpoint-only binning left contacts near the ends of long faces
        // undetected (outside the node's 3x3 stencil).
        Eigen::Vector2d P = X0_[act_[k].na] + u_[act_[k].na];
        Eigen::Vector2d Q = X0_[act_[k].nb] + u_[act_[k].nb];
        int cx0 = std::clamp(int((std::min(P.x(), Q.x()) - gmin_.x()) / cell_), 0, gx_ - 1);
        int cx1 = std::clamp(int((std::max(P.x(), Q.x()) - gmin_.x()) / cell_), 0, gx_ - 1);
        int cy0 = std::clamp(int((std::min(P.y(), Q.y()) - gmin_.y()) / cell_), 0, gy_ - 1);
        int cy1 = std::clamp(int((std::max(P.y(), Q.y()) - gmin_.y()) / cell_), 0, gy_ - 1);
        for (int cy = cy0; cy <= cy1; ++cy)
            for (int cx = cx0; cx <= cx1; ++cx)
                grid_[(std::size_t)cy * gx_ + cx].push_back((int)k);
    }

    double cap = 0.6 * hmin_;                          // deep-pen sanity cap

    // Detection is the expensive part (grid sweep + geometry) and is pure:
    // it parallelizes over nodes into per-thread candidate lists. The
    // DELICATE part — birth-gap bookkeeping (pen0_), damping, the force
    // application and the net-work meter — stays SERIAL, walking the
    // candidates in thread order: deterministic for a fixed thread count,
    // and none of the energy-pump safeguards touches shared state from two
    // threads.
    // ---- geometry of the active edges, computed ONCE per step --------------
    // Same hoisting as the 3D solver: the sweep rebuilt P/Q, the edge length
    // (a sqrt), the outward normal and the local cap for EVERY (node, edge)
    // pair. Arithmetically neutral — same values, same pair list, same order.
    struct EGeo {
        Eigen::Vector2d P, Q, ed, e, nrm, mid;   // Q stored, NOT recomputed as
        double L2, capk, rad2;                   // P + ed (not bit-identical)
        bool ok;
    };
    static std::vector<EGeo> egeo;
    egeo.resize(act_.size());
    for (std::size_t k = 0; k < act_.size(); ++k) {
        EGeo& g = egeo[k];
        g.ok = false;
        if (!inBox[k]) continue;                  // never binned = never a candidate
        const BEdge& be = act_[k];
        g.P = X0_[be.na] + u_[be.na];
        g.Q = X0_[be.nb] + u_[be.nb];
        g.ed = g.Q - g.P;
        g.L2 = g.ed.squaredNorm();
        if (g.L2 < 1e-20) continue;
        double L = std::sqrt(g.L2);
        g.e = g.ed / L;
        g.nrm = Eigen::Vector2d(g.e.y(), -g.e.x());   // outward of be.elem
        g.capk = voronoi_ ? 0.6 * hEl_[be.elem] : cap;
        g.mid = 0.5 * (g.P + g.Q);
        // exact culling disc: an accepted contact projects ONTO the segment
        // (within L/2 of its midpoint) at a depth <= capk
        double rc = 0.5 * L + g.capk;
        g.rad2 = rc * rc;
        g.ok = true;
    }

    struct CPair {
        int i, k;
        double s, d;
        Eigen::Vector2d e, nrm;
    };
    auto detect = [&](int i, std::vector<CPair>& outC,
                      std::vector<int>& stamp) {
        Eigen::Vector2d p = X0_[i] + u_[i];
        if (p.x() <= boxLo.x() || p.x() >= boxHi.x()
            || p.y() <= boxLo.y() || p.y() >= boxHi.y()) return;
        int ci = std::clamp(int((p.x() - gmin_.x()) / cell_), 0, gx_ - 1);
        int cj = std::clamp(int((p.y() - gmin_.y()) / cell_), 0, gy_ - 1);
        int myElem = elemOf_[i];
        // no per-node clear: the dedup stamp below is keyed on the node id,
        // which is unique within one sweep
        for (int dj = -1; dj <= 1; ++dj)
            for (int di = -1; di <= 1; ++di) {
                int cx = ci + di, cy = cj + dj;
                if (cx < 0 || cy < 0 || cx >= gx_ || cy >= gy_) continue;
                for (int k : grid_[(std::size_t)cy * gx_ + cx]) {
                    const BEdge& be = act_[k];
                    if (be.elem == myElem) continue;
                    // AABB binning can list one edge in several stencil
                    // cells: dedup per NODE. The stamp array is PER THREAD
                    // (a shared one caused cache-line ping-pong and DOUBLED
                    // the contact cost at 18 threads) and O(1), where the
                    // former local linear scan was O(candidates) each.
                    if (stamp[k] == i) continue;
                    stamp[k] = i;
                    const EGeo& g = egeo[k];
                    if (!g.ok) continue;
                    if ((p - g.mid).squaredNorm() > g.rad2) continue;  // exact bound
                    // Nodes born from the same virtual vertex as an edge end
                    // are CO-LOCATED with it while the mesh hangs together:
                    // micron-scale compressive interpenetration otherwise
                    // turns every such corner into a phantom follower force —
                    // the energy pump isolated by bisection. The pair is
                    // re-admitted once it has genuinely moved apart.
                    if (vOf_[i] == vOf_[be.na]
                        && (p - g.P).norm() < 0.25 * hmin_) continue;
                    if (vOf_[i] == vOf_[be.nb]
                        && (p - g.Q).norm() < 0.25 * hmin_) continue;
                    double s = (p - g.P).dot(g.ed) / g.L2;
                    if (s < 0.0 || s > 1.0) continue;
                    double d = (p - g.P).dot(g.nrm);
                    if (d >= 0.0 || d < -g.capk) continue; // outside / too deep
                    outC.push_back({i, k, s, d, g.e, g.nrm});
                }
            }
    };

    static std::vector<std::vector<CPair>> cpTL;       // per-thread lists
    // fork/join costs more than the sweep while few faces are active
    // (early in a run only the exterior is): go parallel only when the
    // debris population makes the sweep genuinely heavy
    bool par = false;
#ifdef _OPENMP
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
            const BEdge& be = act_[cp.k];
            const Eigen::Vector2d& e = cp.e;
            const Eigen::Vector2d& nrm = cp.nrm;
            double s = cp.s;
            double pen = -cp.d;
                    // Initial-penetration relief: a pair discovered already
                    // overlapping (rebuild latency, freed faces) must not be
                    // loaded on its full overlap at birth. The reference
                    // penetration is recorded on first sight and relaxed on
                    // a ~1 us time scale, so only NEW approach is resisted.
                    uint64_t pkey = (uint64_t(uint32_t(i)) << 40)
                                    ^ (uint64_t(uint32_t(be.na)) << 20)
                                    ^ uint32_t(be.nb);
                    auto [it0, isNew] = pen0_.try_emplace(pkey, pen);
                    if (!isNew) it0->second *= relax_;
                    pen = std::max(0.0, pen - it0->second);
                    if (pen <= 0.0) continue;
                    Eigen::Vector2d vEdge = (1.0 - s) * v_[be.na] + s * v_[be.nb];
                    Eigen::Vector2d vrel = v_[i] - vEdge;
                    double cdmp = 2.0 * xiGC_ * std::sqrt(kpGC_ * m_[i]);
                    double vn = vrel.dot(nrm);
                    // Quasi-plastic normal contact: the spring pushes at full
                    // stiffness only while APPROACHING. On release it returns
                    // a fraction gcRest of that force. The geometric spring on
                    // rotating segments is a follower force (pen is geometric,
                    // the work is kinematic): full-stiffness release is what
                    // kept pumping energy (net-work meter). Crushed-rock
                    // contacts have restitution ~0.1-0.3 anyway.
                    double fn = kpGC_ * pen * (vn < 0.0 ? 1.0 : gcRest_)
                                - cdmp * vn;
                    if (fn < 0) fn = 0;
                    double ftg = -ctcMu(elemOf_[i], be.elem)   // WP6
                                 * fn * std::tanh(vrel.dot(e) / vReg_);
                    Eigen::Vector2d Fc = fn * nrm + ftg * e;
                    // impulse cap: no single general contact may change this
                    // node's velocity by more than vCap in one step (deep
                    // captures by fast-sweeping faces otherwise act as guns)
                    double capF = 20.0 * m_[i] / dt_;
                    double Fn2 = Fc.norm();
                    if (Fn2 > capF) Fc *= capF / Fn2;
            gcWork_ += Fc.dot(vrel) * dt_;             // net energy meter
            // V2/B4 : part tangentielle (ftg*e est normal a nrm, cap scale)
            gcFricWork_ += (Fc - Fc.dot(nrm) * nrm).dot(vrel) * dt_;
            f_[i] += Fc;
            f_[be.na] -= (1.0 - s) * Fc;
            f_[be.nb] -= s * Fc;
            {   // premiere force portee par une face du pool (diagnostic
                // RKM_GCLOG + source de la regle B en mode adaptatif)
                int pi = actPool_[cp.k];
                if (pi >= 0 && poolTouch_[pi] < 0) poolTouch_[pi] = stepCount_;
            }
            if (gcAdaptive_) {                         // source de la regle B
                lastTouch_[i] = stepCount_;
                lastTouch_[be.na] = stepCount_;
                lastTouch_[be.nb] = stepCount_;
            }
        }
}

// Rigid tool (disc or flat) against every node — FemSolver's contact law.
// Per-node force writes are race-free; the tool reaction is reduced through
// per-thread partial sums added in thread order (deterministic for a fixed
// thread count).
void FdemSolver::toolContact() {
    // BALAI : l'outil est RETIRE pendant la phase de classification. Sans
    // cela il continuerait a couper et fabriquerait de nouveaux fragments
    // pendant qu'on essaie de classer les anciens. C'est aussi ce que fait
    // l'experience : on balaie apres l'impact, outil ecarte.
    if (brushArmed_) return;
    // ETAPE 1 : arret de l'outil AVANT l'armement, pour menager un intervalle
    // de repos. Voir FdemSolver.hpp — sans lui, des fragments encore lances a
    // 100 m/s traversent le bloc pendant le balayage.
    if (toolStop_ > 0.0 && t_ >= toolStop_) {
        if (!toolStopped_) {
            toolStopped_ = true;
            // REPARATION (2026-08-28) : arret REEL — voir le miroir 3D.
            tool_.v.setZero();
            std::cout << "[FDEM] OUTIL ARRETE (v = 0) a t = " << t_ << " s. Repos "
                         "jusqu'a l'armement du balai (t = " << brushStart_
                      << " s), soit " << (brushStart_ - t_) << " s pour que "
                         "l'amortissement eteigne les vitesses residuelles.\n";
        }
        return;
    }
    // BRAZILIAN must be excluded here as well as TENSION. placeTool() returns
    // early for both, so tool_ keeps its STRUCT DEFAULTS — a disc of radius
    // 10 mm sitting at the origin — and this routine would press that phantom
    // obstacle into the specimen at kp_ = E*t. It bit exactly once and hard:
    // the disc centre sits at (W/2, H/2), so the clearance between the origin
    // and the rim is 0.2071*W, i.e. 10.36 mm for a 50 mm specimen against a
    // 10 mm tool. The Red Bohus disc (W = 50 mm) therefore had 82 joints break
    // at t = 4e-7 s in the LOWER-LEFT quadrant at 40-50 deg under zero applied
    // load, while the 60 mm verification disc (clearance 12.43 mm) was clean.
    if (scen_ == Scenario::TENSION || scen_ == Scenario::BRAZILIAN
        || scen_ == Scenario::SHPB) return;
    auto nodeFc = [&](int i, Eigen::Vector2d& Fc) {
        Eigen::Vector2d p = X0_[i] + u_[i];
        if (tool_.shape == Tool::Shape::PDC) {
            // Rigid wedge: the cutter is the half-space BEHIND the rake face
            // and ABOVE the cutting edge. A node is in contact when it is on
            // the wrong side of the rake face, still within the face extent,
            // and above the edge — so the rock below the depth of cut is
            // untouched, exactly as a real cutter leaves it.
            Eigen::Vector2d n = tool_.rakeNormal();
            Eigen::Vector2d tdir = tool_.rakeDir();
            Eigen::Vector2d rel = p - tool_.x;          // from the edge
            double dn = rel.dot(n);                     // >0 = in front, free
            double dt2 = rel.dot(tdir);                 // along the face
            if (dn >= 0.0) return false;
            // BACK FACE (cutterThick > 0): the cutter has a finite thickness,
            // so a node more than `thick` behind the rake face is OUTSIDE the
            // tool -- it is in the void the cutter left behind. Without this
            // bound the wedge is a half-space and traps its own cut floor and
            // chip, which the capped penalty then pumps to km/s.
            if (tool_.thick > 0.0 && dn < -tool_.thick) return false;
            if (dt2 < 0.0 || dt2 > tool_.faceLen) return false;
            // DEEP-PENETRATION CAP on the local element size, exactly as the
            // general contact does. Capping on the face length instead (6.5 mm
            // on a 13 mm cutter) let a node of the chip end up millimetres
            // inside the wedge and produced isolated 32 MN/m spikes on an
            // otherwise 0.5 MN/m curve — the penalty force is linear in the
            // penetration, so a node that slips deep dominates everything.
            double pen = -dn;
            double capPen = 0.6 * hEl_[elemOf_[i]];
            if (pen > capPen) pen = capPen;
            Eigen::Vector2d vrel = v_[i] - tool_.v;
            double c = 2.0 * xiC_ * std::sqrt(kp_ * m_[i]);
            double fn = kp_ * pen - c * vrel.dot(n);
            if (fn < 0) fn = 0;
            double ftg = -ctcMu(elemOf_[i]) * fn               // WP6
                         * std::tanh(vrel.dot(tdir) / vReg_);
            Fc = fn * n + ftg * tdir;
        } else if (tool_.shape == Tool::Shape::FLAT) {
            if (std::abs(p.x() - tool_.x.x()) > 0.5 * tool_.width) return false;
            double pen = p.y() - tool_.x.y();
            if (pen <= 0) return false;
            double c = 2.0 * xiC_ * std::sqrt(kp_ * m_[i]);
            double fn = kp_ * pen + c * (v_[i].y() - tool_.v.y());
            if (fn < 0) fn = 0;
            double ftg = -ctcMu(elemOf_[i]) * fn               // WP6
                         * std::tanh((v_[i].x() - tool_.v.x()) / vReg_);
            Fc = {ftg, -fn};
        } else {
            Eigen::Vector2d d = p - tool_.x;
            double dist = d.norm();
            if (dist >= tool_.radius || dist < 1e-14) return false;
            Eigen::Vector2d n = d / dist;
            Eigen::Vector2d tdir(-n.y(), n.x());
            double pen = tool_.radius - dist;
            Eigen::Vector2d vrel = v_[i] - tool_.v;
            double c = 2.0 * xiC_ * std::sqrt(kp_ * m_[i]);
            double fn = kp_ * pen - c * vrel.dot(n);
            if (fn < 0) fn = 0;
            double ftg = -ctcMu(elemOf_[i]) * fn               // WP6
                         * std::tanh(vrel.dot(tdir) / vReg_);
            Fc = fn * n + ftg * tdir;
        }
        // --- A : ECRETAGE EN IMPULSION (voir FdemSolver.hpp) ----------------
        // L'ecretage historique borne la PENETRATION (0,6 h) ; il ne borne pas
        // l'IMPULSION. Mesure du 2026-08-18 (coupe PDC) : la force ainsi
        // plafonnee, 7,24 MN/m, communique encore 377 m/s en UN pas a un noeud
        // de 2,46e-5 kg/m. On borne donc la grandeur qui lance les noeuds.
        // La borne est PHYSIQUE, pas un reglage : un corps rigide de masse
        // infinie anime de v ne peut, en choc parfaitement elastique,
        // communiquer plus de 2v. On ecrete le VECTEUR complet, frottement
        // compris, car c'est l'impulsion totale qui accelere le noeud.
        if (toolVCap_ > 0.0) {
            double vt = tool_.v.norm();
            if (vt > 1e-12) {                  // outil a l'arret : pas de borne
                double fmax = toolVCap_ * 2.0 * vt * m_[i] / dt_;
                double f = Fc.norm();
                if (f > fmax) Fc *= fmax / f;
            }
        }
        return true;
    };

    // -----------------------------------------------------------------------
    // A1 — variante SIGNORINI EN VITESSE (toolContact = signorini).
    //
    // La geometrie est DUPLIQUEE a dessein plutot que factorisee avec nodeFc :
    // la voie penalite doit rester bit-identique, et un refactor la toucherait.
    //
    // Algorithme, par noeud (Delassus diagonal, H = 1/m_i) :
    //   1. vitesse LIBRE : v* = v_i + (dt/m_i) f_i, avec f_i le total des
    //      autres forces deja assemblees ;
    //   2. vitesse relative normale a l'outil : vn = (v* - v_outil).n ;
    //   3. gap predit : g+ = -pen + dt vn. Si g+ >= 0 le noeud se separe tout
    //      seul, l'impulsion est NULLE (condition de Signorini) ;
    //   4. sinon impulsion repulsive r_n = m_i (v_cible - vn) >= 0, avec
    //      v_cible = relax * pen / dt (0 = pas de rattrapage : on annule
    //      l'approche, on ne resorbe pas la penetration existante) ;
    //   5. frottement : impulsion de COLLAGE r_t = -m_i vt, ecretee par le
    //      cap de Coulomb |r_t| <= mu r_n. Pas de regularisation en tanh :
    //      le cap porte sur l'impulsion, il n'a pas besoin d'etre lisse.
    //   6. report en FORCE, r/dt, pour que integrate() (v += dt/m f) produise
    //      exactement le saut de vitesse voulu et que les compteurs d'energie
    //      restent comparables a la voie penalite.
    //
    // La borne physique est ainsi respectee PAR CONSTRUCTION : l'impulsion ne
    // peut qu'annuler l'approche, donc un noeud ne peut jamais repartir a plus
    // de 2 v_outil. Aucun reglage, aucune raideur.
    // -----------------------------------------------------------------------
    auto nodeSig = [&](int i, Eigen::Vector2d& Fc) {
        Eigen::Vector2d p = X0_[i] + u_[i];
        Eigen::Vector2d n, tdir;
        double pen;
        if (tool_.shape == Tool::Shape::PDC) {
            n = tool_.rakeNormal();
            tdir = tool_.rakeDir();
            Eigen::Vector2d rel = p - tool_.x;
            double dn = rel.dot(n), dt2 = rel.dot(tdir);
            if (dn >= 0.0) return false;
            if (tool_.thick > 0.0 && dn < -tool_.thick) return false;
            if (dt2 < 0.0 || dt2 > tool_.faceLen) return false;
            pen = -dn;
        } else if (tool_.shape == Tool::Shape::FLAT) {
            if (std::abs(p.x() - tool_.x.x()) > 0.5 * tool_.width) return false;
            pen = p.y() - tool_.x.y();
            if (pen <= 0.0) return false;
            n = {0.0, -1.0};                   // repousse le noeud vers le bas
            tdir = {1.0, 0.0};
        } else {
            Eigen::Vector2d d = p - tool_.x;
            double dist = d.norm();
            if (dist >= tool_.radius || dist < 1e-14) return false;
            n = d / dist;
            tdir = {-n.y(), n.x()};
            pen = tool_.radius - dist;
        }
        // (1) vitesse libre — f_[i] porte deja toutes les autres forces
        Eigen::Vector2d vFree = v_[i] + (dt_ / m_[i]) * f_[i];
        Eigen::Vector2d vrel = vFree - tool_.v;
        double vn = vrel.dot(n);
        // (3) le noeud se separe-t-il de lui-meme au pas suivant ?
        if (-pen + dt_ * vn >= 0.0) return false;
        // (4) impulsion normale, positive par construction
        double vTarget = toolSigRelax_ * pen / dt_;
        double rn = m_[i] * (vTarget - vn);
        if (rn <= 0.0) return false;
        // (5) frottement de Coulomb sur l'IMPULSION
        double vt = vrel.dot(tdir);
        double rt = -m_[i] * vt;
        double cap = ctcMu(elemOf_[i]) * rn;                   // WP6
        if (rt > cap) rt = cap;
        else if (rt < -cap) rt = -cap;
        // (6) report en force
        Fc = (rn / dt_) * n + (rt / dt_) * tdir;
        return true;
    };
#ifdef _OPENMP
    int nT = omp_get_max_threads();
    std::vector<Eigen::Vector2d> FT(nT, Eigen::Vector2d::Zero());
    double tw = 0.0;                       // V2/B4 : travail outil -> solide
#pragma omp parallel reduction(+:tw)
    {
        int t = omp_get_thread_num();
        Eigen::Vector2d Floc = Eigen::Vector2d::Zero();
#pragma omp for schedule(static)
        for (int i = 0; i < (int)X0_.size(); ++i) {
            Eigen::Vector2d Fc;
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
        Eigen::Vector2d Fc;
        if (!(toolSig_ ? nodeSig(i, Fc) : nodeFc(i, Fc))) continue;
        f_[i] += Fc;
        toolWork_ += Fc.dot(v_[i]) * dt_;  // V2/B4
        tool_.F -= Fc;
    }
#endif
}

// ---------------------------------------------------------------------------
// Brazilian traction pair. Same follower-load machinery as the confinement
// (current face, lumped L/2 per node), driven by a pressure that grows at
// loadRate. The two arcs are loaded with the SAME pressure, so the resultant
// pair is self-balancing to within the mesh's top/bottom asymmetry — reported
// as arcTop_.F + arcBot_.F, which the history carries so drift is visible
// rather than silent.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Rigid platen pair in frictional penalty contact. Shared by the brazilian and
// by the UCS/triaxial test when loading = platens: in both cases the specimen
// is held by CONTACT, never clamped, which is how Y-Geo and the FEM-DEM
// calibration literature load laboratory specimens. Both platens close
// symmetrically so the specimen does not drift and the Cundall damping forces
// of the two halves cancel.
// ---------------------------------------------------------------------------
void FdemSolver::platenForces() {
    // ramp resolved here so the dashpot sees the velocity the platen
    // actually has this step; integrate() reuses it to advance the planes
    double vg = 0.5 * platenV_;                    // each platen, inward
    // pullDelay : le chargement axial ne demarre qu'a t = pullDelay. En
    // triaxial la pression de confinement doit etre etablie et equilibree
    // AVANT que les plateaux ne se ferment, sinon la courbe contrainte-
    // deformation demarre sur un etat qui n'est ni uniaxial ni confine.
    // Defaut 0 : tous les modeles existants sont bit-identiques.
    double tl = t_ - pullDelay_;
    if (tl <= 0.0) vg = 0.0;
    else if (pullRamp_ > 0.0 && tl < pullRamp_)
        vg *= 0.5 * (1.0 - std::cos(M_PI * tl / pullRamp_));
    plTop_.v = -vg;                                // downward
    plBot_.v = +vg;                                // upward
    plTop_.F.setZero();
    plBot_.F.setZero();
    double sF = 0.0, sF2 = 0.0;                    // participation ratio
    long nC = 0;
    double pw = 0.0;                               // V2/B4 : platine -> solide
    for (Platen* pl : {&plTop_, &plBot_}) {
        for (int i = 0; i < (int)X0_.size(); ++i) {
            if (platenW_[i] <= 0.0) continue;      // not a bearing node
            Eigen::Vector2d p = X0_[i] + u_[i];
            if (std::abs(p.x() - pl->xc) > pl->halfW) continue;
            double pen = -pl->sign * (p.y() - pl->y);
            if (pen <= 0.0) continue;
            double kNode = kpPlaten_ * platenW_[i];
            double c = 2.0 * xiC_ * std::sqrt(kNode * m_[i]);
            // d(pen)/dt = -(v_node - v_platen) for BOTH platens, so the
            // dashpot opposes -vrel in both cases
            double vrel = v_[i].y() - pl->v;
            double fn = kNode * pen - c * vrel;
            if (fn < 0.0) fn = 0.0;                // no adhesion on release
            // WP6 : volontairement PAS de ctcMu ici — le plateau est une
            // frontiere de machine, pas un support de fragments.
            double ftg = -muC_ * fn * std::tanh((v_[i].x()) / vReg_);
            Eigen::Vector2d Fc(ftg, pl->sign * fn);   // force ON the node
            f_[i] += Fc;
            pl->F -= Fc;                           // reaction ON the platen
            pw += Fc.dot(v_[i]);                   // V2/B4 (2026-08-15) :
            // les plateaux RIGIDES poussent par contact penalise — leur
            // travail sur le solide n'etait compte NULLE PART (4e trou
            // designe par le moniteur E2, apres confinement 2D/3D et mors).
            // Meme convention que les autres familles : (force appliquee au
            // noeud) . v-, le biais O(dt) est couvert par biasW_.
            if (pl->sign < 0) { sF += fn; sF2 += fn * fn; ++nC; }
        }
    }
    bcWork_ += pw * dt_;                           // poste « platines/grips »
    // How many nodes actually carry the bearing? The participation ratio
    // (sum f)^2 / sum f^2 is that number: 1 if a single asperity takes
    // everything, nC if the load is uniform. Latched with the gauge.
    if (nBroken_ == 0 && sF2 > 0.0) {
        gPR_ = sF * sF / sF2;
        gNC_ = nC;
    }
}

// ---------------------------------------------------------------------------
// End of the brazilian test. peakLocked_ marks the first genuine post-cracking
// load drop: past it the platens keep closing on two separated halves and the
// nominal sigma_t rises again for a kinematic reason (measured on
// bd_yan_adaptatif: 1.30 MPa at the bottom of the drop, 2.10 MPa afterwards,
// unchanged with the general contact switched off, so it is the crushing of
// the halves and not a contact artefact). The published brazilian curve stops
// there. brazStopDelay_ keeps enough steps after the lock for the drop itself
// to be on the curve.
// E2 (fiabilite) : moniteur d'energie runtime, opt-in par budgetAbortPct —
// meme definition de residu/echelle que le resume B4, arret PROPRE via
// finished() (frame + history + summary ecrits). Voir Fdem3dSolver pour le
// commentaire complet.
void FdemSolver::checkEnergyAbort() {
    if (eAbortPct_ < 0.0) {
        eAbortPct_ = cfg_.getd("budgetAbortPct", 0.0);
        eAbortMin_ = cfg_.getd("budgetAbortMin", 0.0);   // plancher absolu [J/m]
    }
    if (eAbortPct_ <= 0.0 || eAbort_) return;
    double ke = 0.0;
    for (std::size_t i = 0; i < X0_.size(); ++i)
        ke += 0.5 * m_[i] * v_[i].squaredNorm();
    double sumW = elWork_ + jointWork_ + gcWork_ + cundWork_ + lysWork_
                + toolWork_ + bcWork_ + confWork_ + hydroWork_ + biasW_;
    double gross = std::abs(elWork_) + std::abs(jointWork_)
                 + std::abs(gcWork_) + std::abs(cundWork_)
                 + std::abs(lysWork_) + std::abs(toolWork_)
                 + std::abs(bcWork_) + std::abs(confWork_)
                 + std::abs(hydroWork_);
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
    std::cout << "[FDEM] ENERGY ABORT (budgetAbortPct = " << eAbortPct_
              << ") a t = " << t_ << " s : residu B4 " << resid << " J/m = "
              << 100.0 * std::abs(resid) / scale
              << " % de l'echelle. Hotspot : noeud " << iw << ", |v| = "
              << std::sqrt(vw) << " m/s\n";
    eAbort_ = true;
}

bool FdemSolver::finished() const {
    if (eAbort_) return true;              // E2 : moniteur d'energie
    // meme mecanisme cote compression (ucsStopAfterPeak), voir FdemSolver.hpp
    if (ucsStop_ && scen_ == Scenario::TENSION && tensionPlatens_)
        return peakLockedU_ && tLockedU_ >= 0.0
               && t_ >= tLockedU_ + ucsStopDelay_;
    if (!brazStop_ || scen_ != Scenario::BRAZILIAN) return false;
    return peakLocked_ && tLocked_ >= 0.0 && t_ >= tLocked_ + brazStopDelay_;
}

void FdemSolver::brazilianForces() {
    if (scen_ != Scenario::BRAZILIAN) return;

    if (brazPlatens_) { platenForces(); return; }
    brazP_ = loadRate_ * t_;
    for (LoadArc* arc : {&arcTop_, &arcBot_}) {
        arc->F.setZero();
        for (const auto& be : arc->edge) {
            Eigen::Vector2d P = X0_[be.na] + u_[be.na];
            Eigen::Vector2d Q = X0_[be.nb] + u_[be.nb];
            Eigen::Vector2d d = Q - P;
            double L = d.norm();
            if (L < 1e-14) continue;
            Eigen::Vector2d n(d.y() / L, -d.x() / L);   // outward unit normal
            Eigen::Vector2d Fe = -brazP_ * L * thk_ * n;   // inward
            f_[be.na] += 0.5 * Fe;
            f_[be.nb] += 0.5 * Fe;
            arc->F += Fe;
        }
    }
}

// ---------------------------------------------------------------------------
// Confining pressure as a FOLLOWER load: the traction -p n is recomputed on
// the CURRENT face (position and orientation both move), lumped L/2 per node.
// (P -> Q) is the CCW traversal of the owning element, so ((Q-P).y, -(Q-P).x)
// points OUT of the specimen — the same convention buildFromTriangles fixed
// for the exterior faces, and the reason this loop can trust it.
// ---------------------------------------------------------------------------
void FdemSolver::confiningForces() {
    if (confP_ == 0.0) return;
    double p = confP_;
    if (confRamp_ > 0.0 && t_ < confRamp_)
        p *= 0.5 * (1.0 - std::cos(M_PI * t_ / confRamp_));
    for (const auto& be : confEdges_) {
        Eigen::Vector2d P = X0_[be.na] + u_[be.na];
        Eigen::Vector2d Q = X0_[be.nb] + u_[be.nb];
        Eigen::Vector2d d = Q - P;
        double L = d.norm();
        if (L < 1e-14) continue;
        Eigen::Vector2d n(d.y() / L, -d.x() / L);       // outward unit normal
        Eigen::Vector2d half = -0.5 * p * L * thk_ * n; // inward, half per node
        f_[be.na] += half;
        f_[be.nb] += half;
        // V2/B4 : travail de la pression sur le solide (compteur en v-)
        confWork_ += dt_ * (half.dot(v_[be.na]) + half.dot(v_[be.nb]));
    }
}

// ===========================================================================
// COUPLAGE HYDRO-MECANIQUE (spec 004) — modele d'AbuAisha et al. 2017, Y-Geo
// ===========================================================================

// Selection des faces SOURCE : les memes que confineFaces = bore, mais pour un
// role oppose. Le confinement presse une membrane sur l'exterieur d'origine ;
// ici, le forage est la BOUCHE par ou le fluide entre, et tout ce qui s'y
// connecte par une fissure deviendra mouille.
void FdemSolver::setupHydro() {
    hydroOn_ = cfg_.getb("hydro", false);
    if (!hydroOn_) return;
    if (scen_ == Scenario::BRAZILIAN || scen_ == Scenario::SHPB)
        throw std::runtime_error("hydro = on n'a pas de sens en scenario "
                                 "brazilian ou shpb");
    fluidK_    = cfg_.getd("fluidBulk", 2.2e9);
    fluidRho0_ = cfg_.getd("fluidDensity", 1000.0);
    // ESSAI 2 (2026-08-20). Seuil d endommagement a partir duquel une
    // interface conduit le fluide. Defaut 1.0 = comportement historique,
    // seules les interfaces ROMPUES sont mouillees ; le chemin est alors
    // bit-identique. Poser 0.0 mouille toute interface inseree des que
    // son endommagement demarre, ce qui supprime l incubation a sec.
    wetDmin_ = cfg_.getd("hydroWetDamage", 1.0);
    if (wetDmin_ < 0.0 || wetDmin_ > 1.0)
        throw std::runtime_error("hydroWetDamage doit etre dans [0, 1] : "
                                 "1 = seules les interfaces rompues "
                                 "conduisent (defaut), 0 = toute "
                                 "interface inseree conduit");
    hydroP0_   = cfg_.getd("hydroP0", 0.0);
    hydroRamp_ = cfg_.getd("hydroRamp", 0.0);
    // ESSAI 3 (2026-08-21). Protocole d AbuAisha 2017, section 3.2 : le pas
    // geostatique et l excavation sont ENTIEREMENT anterieurs a l injection.
    // Avant hydroStart la pompe est a l arret, p = p0, et la masse de fluide
    // est re-basee en continu sur le volume COURANT de la cavite, si bien que
    // l injection demarre de l etat excave et converge — pas du volume de
    // maillage a t = 0. Defaut 0 = comportement historique, bit-identique.
    hydroStart_ = cfg_.getd("hydroStart", 0.0);
    if (hydroStart_ < 0.0)
        throw std::runtime_error("hydroStart doit etre >= 0 [s] : retard du\n"
                                 "  demarrage de la pompe (l excavation doit\n"
                                 "  etre bouclee avant cet instant)");
    if (!(fluidK_ > 0.0) || !(fluidRho0_ > 0.0))
        throw std::runtime_error("fluidBulk et fluidDensity doivent etre > 0");
    std::string inj = cfg_.gets("hydroInjection", "rate");
    if (inj != "rate" && inj != "pressure")
        throw std::runtime_error("hydroInjection must be rate | pressure");
    hydroRateMode_ = (inj == "rate");
    hydroRate_ = cfg_.getd("hydroRate", 0.0);
    hydroPimp_ = cfg_.getd("hydroPressure", 0.0);

    // --- faces source. Deux modes : un disque autour de (boreCX, boreCY),
    // comme le confinement ; ou TOUTES les faces exterieures d'origine, ce qui
    // sert au controle de Parker ou la « source » est la discontinuite
    // elle-meme, dedoublee par le greffon Crack de gmsh (ses deux levres sont
    // des faces exterieures confondues a t = 0).
    std::string src = cfg_.gets("hydroSource", "bore");
    if (src != "bore" && src != "all")
        throw std::runtime_error("hydroSource must be bore | all");
    double bcx = cfg_.getd("boreCX", 0.5 * W_);
    double bcy = cfg_.getd("boreCY", 0.5 * H_);
    double bsr = cfg_.getd("boreSelectR", 0.0);
    if (src == "bore" && bsr <= 0.0)
        throw std::runtime_error("hydroSource = bore exige boreSelectR > 0 "
                                 "(rayon de selection autour de boreCX/boreCY)");
    for (const auto& be : exterior_) {
        if (src == "all") { hydroSrc_.push_back(be); continue; }
        Eigen::Vector2d mid = 0.5 * (X0_[be.na] + X0_[be.nb]);
        double dx = mid.x() - bcx, dy = mid.y() - bcy;
        if (dx * dx + dy * dy <= bsr * bsr) hydroSrc_.push_back(be);
    }
    if (hydroSrc_.empty())
        throw std::runtime_error("hydro : aucune face source selectionnee — "
                                 "verifier boreCX/boreCY/boreSelectR");
    updateWetBoundary();
    hydroVol0_ = hydroVol_;
    // ---- VOLUME INITIAL NUL : ce n'est pas toujours une erreur ------------
    // Une discontinuite d'epaisseur nulle — la fissure de Parker, dont les
    // deux levres sont CONFONDUES a t = 0 — enferme rigoureusement zero aire.
    // Le lacet somme donc a zero, et c'est geometriquement juste.
    //   * en pression imposee, cela n'a aucune consequence : la pression ne
    //     depend pas du volume, qui n'est plus qu'une sortie ;
    //   * en debit impose, la pression passe par log(m / V rho_0) : un volume
    //     nul rend le modele de compressibilite INDEFINI. Il faut alors une
    //     cavite physique de volume fini — un forage — ou une ouverture
    //     initiale non nulle a la Lisjak.
    // C'est d'ailleurs pour cela que l'annexe A d'AbuAisha impose une pression
    // UNIFORME et ne fait pas tourner sa pompe sur ce cas.
    if (!(hydroVol0_ > 0.0)) {
        if (hydroRateMode_)
            throw std::runtime_error(
                "hydro : volume de cavite initial nul ou negatif ("
                + std::to_string(hydroVol0_) + " m3/m) en injection a DEBIT "
                "impose — la pression p0 + Kf log(m / V rho0) y est indefinie. "
                "Il faut une cavite de volume fini (forage) ou une ouverture "
                "initiale. En pression imposee, ce cas est licite.");
        std::cout << "[HYDRO] NOTE : volume de cavite initial nul ("
                  << hydroVol0_ << " m3/m). Normal pour une discontinuite "
                     "d'epaisseur nulle dont les levres sont confondues ; sans "
                     "consequence en pression imposee.\n";
    }
    // masse initiale telle que p(t=0) = p0, par inversion de leur eq. 6
    hydroMass_ = hydroVol0_ * fluidRho0_;
    hydroP_ = hydroP0_;
    std::cout << "[HYDRO] couplage hydro-mecanique ACTIF (AbuAisha 2017, "
                 "fluide NON VISQUEUX : pression uniforme dans la cavite)\n"
              << "[HYDRO]   source = " << src << " : " << hydroSrc_.size()
              << " faces, volume initial " << hydroVol0_ << " m3/m\n"
              << "[HYDRO]   K_f = " << fluidK_ << " Pa, rho_f0 = " << fluidRho0_
              << " kg/m3, p_0 = " << hydroP0_ << " Pa\n"
              << "[HYDRO]   pompe = " << inj << " ("
              << (hydroRateMode_ ? hydroRate_ : hydroPimp_)
              << (hydroRateMode_ ? " m3/s/m)" : " Pa)")
              << (hydroRamp_ > 0.0 ? ", rampe " : "")
              << (hydroRamp_ > 0.0 ? std::to_string(hydroRamp_) + " s" : "")
              << (hydroStart_ > 0.0 ? ", demarrage a " : "")
              << (hydroStart_ > 0.0 ? std::to_string(hydroStart_) + " s" : "")
              << "\n";
}

// Frontiere mouillee : les faces reliees a la source par un chemin de joints
// ROMPUS, plus le volume de cavite qu'elles enferment.
//
// C'est leur module 2, « tracks the newly created wet boundaries by checking
// their connection with the initial source of fluid ». En pratique une
// recherche de composante connexe sur le graphe des SOMMETS (au sens de
// l'union-find, pas des noeuds dupliques), dont les aretes sont les joints
// rompus. La machinerie de sommets existe deja : vOf_ et nVert_.
void FdemSolver::updateWetBoundary() {
    // (1) propagation de la mouillabilite, seulement si la topologie a change.
    //
    // ESSAI 2 (2026-08-20). Le timbre historique est nBroken_ : il suffit tant
    // que seule une interface ROMPUE conduit, puisque la topologie mouillee ne
    // peut alors changer qu a une rupture. Des que wetDmin_ < 1 cette equivalence
    // tombe — une interface se met a conduire quand son endommagement franchit
    // le seuil, sans rupture — et un timbre sur nBroken_ laisserait le front
    // mouille fige entre deux ruptures, ce qui viderait l essai de son objet.
    // On timbre alors sur le nombre d interfaces conductrices. La passe est
    // O(nJoints) mais n a lieu que dans ce mode.
    long wetKey = nBroken_;
    if (wetDmin_ < 1.0) {
        long nc = 0;
        for (const auto& J : jt_)
            if (!J.bonded && J.D >= wetDmin_) ++nc;
        wetKey = nc;
    }
    if (wetStamp_ != wetKey || wetEdges_.empty()) {
        std::vector<char> wetV(std::max(1, nVert_), 0);
        for (const auto& be : hydroSrc_) {
            wetV[vOf_[be.na]] = 1;
            wetV[vOf_[be.nb]] = 1;
        }
        // relaxation jusqu'a stabilite : une fissure longue demande autant de
        // passes que d'aretes en enfilade. Recalcule seulement quand nBroken_
        // change, donc rarement — mesurer avant d'optimiser.
        bool changed = true;
        while (changed) {
            changed = false;
            for (const auto& J : jt_) {
                // ESSAI 2 : seuil reglable. Avec wetDmin_ = 1 (defaut)
                // seule une interface ROMPUE conduit ; abaisser le seuil
                // fait conduire les interfaces en adoucissement.
                if (J.bonded || J.D < wetDmin_) continue;  // etanche
                int v1 = vOf_[J.a1], v2 = vOf_[J.a2];
                if (wetV[v1] && !wetV[v2]) { wetV[v2] = 1; changed = true; }
                else if (wetV[v2] && !wetV[v1]) { wetV[v1] = 1; changed = true; }
            }
        }
        wetEdges_ = hydroSrc_;
        wetJoint_.assign(jt_.size(), 0);
        wetJointIdx_.clear();
        for (std::size_t jI = 0; jI < jt_.size(); ++jI) {
            const Joint& J = jt_[jI];
            if (J.bonded || J.D < wetDmin_) continue;      // ESSAI 2
            if (!(wetV[vOf_[J.a1]] && wetV[vOf_[J.a2]])) continue;
            wetJoint_[jI] = 1;
            wetJointIdx_.push_back((int)jI);
            // les deux levres, orientees chacune vers l'exterieur de SON
            // element — meme convention que le cache deadList_ du contact
            wetEdges_.push_back({J.eA, J.a1, J.a2});
            wetEdges_.push_back({J.eB, J.b2, J.b1});
        }
        wetStamp_ = wetKey;
        hydroNWet_ = (long)wetEdges_.size();
        // --- diagnostic H5, temporaire : pourquoi la cavite ne grandit pas ---
        static int dbg = 0;
        if (std::getenv("RKM_HYDRODBG") && dbg < 6) {
            ++dbg;
            long nD1 = 0, nOne = 0, nBoth = 0;
            for (const auto& J : jt_) {
                if (J.bonded || J.D < 1.0) continue;
                ++nD1;
                bool w1 = wetV[vOf_[J.a1]], w2 = wetV[vOf_[J.a2]];
                if (w1 || w2) ++nOne;
                if (w1 && w2) ++nBoth;
            }
            std::cout << "[H5] t=" << t_ << " nBroken=" << nBroken_
                      << " joints D>=1 : " << nD1
                      << " | au moins un sommet mouille : " << nOne
                      << " | les deux : " << nBoth
                      << " | faces mouillees : " << wetEdges_.size()
                      << " (nVert=" << nVert_ << ")\n";
            // ou sont les sommets ENSEMENCES, et ou sont les joints rompus ?
            double rs0 = 1e30, rs1 = 0.0;
            for (const auto& be : hydroSrc_) {
                Eigen::Vector2d m = 0.5 * (X0_[be.na] + X0_[be.nb]);
                double r = (m - Eigen::Vector2d(4.0, 4.0)).norm();
                rs0 = std::min(rs0, r); rs1 = std::max(rs1, r);
            }
            std::cout << "[H5]   faces source : r dans [" << rs0 << ", "
                      << rs1 << "] m, premiers sommets ensemences : ";
            int shown = 0;
            for (const auto& be : hydroSrc_) {
                if (shown++ >= 3) break;
                std::cout << vOf_[be.na] << "/" << vOf_[be.nb] << " ";
            }
            std::cout << "| sommets des joints rompus : ";
            shown = 0;
            for (std::size_t q = 0; q < jt_.size() && shown < 4; ++q) {
                if (jt_[q].bonded || jt_[q].D < 1.0) continue;
                ++shown;
                Eigen::Vector2d m = 0.5 * (X0_[jt_[q].a1] + X0_[jt_[q].a2]);
                std::cout << vOf_[jt_[q].a1] << "/" << vOf_[jt_[q].a2]
                          << "(r=" << (m - Eigen::Vector2d(4.0, 4.0)).norm()
                          << ") ";
            }
            std::cout << "\n";
        }
    }
    // (2) volume par le theoreme de Green, sur la configuration COURANTE.
    //
    // ⚠️ LECON DU 2026-08-20, payee par un run divergent. « La somme du lacet
    // face par face est independante de l'ordre » n'est vraie que pour un
    // contour FERME. Chaque terme x_P y_Q - x_Q y_P depend de l'ORIGINE du
    // repere ; seule la somme sur un contour ferme est intrinseque. Le forage
    // d'AbuAisha etant centre en (4, 4), chaque terme vaut ~|r|^2 = 32 m^2
    // contre une aire reelle de 7,8e-3 m^2 : le moindre defaut de fermeture
    // produisait une erreur ENORME. Tant que seul le forage etait mouille le
    // contour fermait, d'ou le succes de H2 a 5,6e-07 ; des que les levres de
    // fissure s'ajoutaient, le volume DECROISSAIT au lieu de croitre, et comme
    // K_f / V vaut 2,8e11 Pa/m^2, 1,9 % d'erreur de volume faisaient 42 MPa
    // d'erreur de pression. La pression montait a la rupture au lieu de
    // chuter.
    //
    // Deux protections, dans cet ordre :
    //   (a) referencer le lacet au CENTROIDE de la frontiere mouillee. Pour un
    //       contour ferme cela ne change rien ; pour un contour imparfait cela
    //       ramene l'erreur de l'ordre de |r|^2 a l'ordre de la cavite ;
    //   (b) CONTROLER la fermeture : la somme vectorielle des segments
    //       orientes doit etre nulle. Si elle ne l'est pas, le contour n'est
    //       pas la frontiere d'un domaine et son aire n'a AUCUN sens — il faut
    //       le dire, pas rendre un nombre.
    // DECOMPOSITION LOCALE (2026-08-20). Le lacet global est abandonne pour le
    // volume : il n'est intrinseque que sur un contour ferme, et un contour de
    // fissure ouverte est trop fragile pour porter une grandeur amplifiee par
    // K_f / V = 2,8e11 Pa/m^2. On somme donc DEUX contributions independantes :
    //
    //   (a) la CAVITE SOURCE, par le lacet sur ses seules faces d'origine —
    //       contour ferme par construction, et verifie a 5,6e-07 pres contre
    //       l'aire exacte du polygone (controle H2) ;
    //   (b) chaque FISSURE mouillee, par son aire propre L * (ouverture
    //       moyenne), calculee LOCALEMENT sur le joint.
    //
    // C'est exactement la forme qu'emploie Lisjak et al. 2017 pour le volume
    // de cavite, V = somme_j L_j (a_0 + a_1)/2. Elle est insensible a
    // l'origine du repere, insensible a la fermeture, et chaque terme est
    // positif par construction — donc une fissure qui s'ouvre ne peut PAS
    // faire decroitre le volume, ce que l'ancienne formule autorisait.
    Eigen::Vector2d cen = Eigen::Vector2d::Zero();
    for (const auto& be : hydroSrc_)
        cen += 0.5 * (X0_[be.na] + u_[be.na] + X0_[be.nb] + u_[be.nb]);
    if (!hydroSrc_.empty()) cen /= (double)hydroSrc_.size();
    double A2 = 0.0, perim = 0.0;
    Eigen::Vector2d closure = Eigen::Vector2d::Zero();
    for (const auto& be : hydroSrc_) {                 // (a) la source seule
        Eigen::Vector2d P = X0_[be.na] + u_[be.na] - cen;
        Eigen::Vector2d Q = X0_[be.nb] + u_[be.nb] - cen;
        A2 += P.x() * Q.y() - Q.x() * P.y();
        closure += (Q - P);
        perim += (Q - P).norm();
    }
    double vol = -0.5 * A2 * thk_;
    // (b) les fissures mouillees, aire propre de chaque levre ecartee
    double vCrack = 0.0;
    for (int jI : wetJointIdx_) {                      // liste, pas balayage
        const Joint& J = jt_[jI];
        Eigen::Vector2d P = 0.5 * (X0_[J.a1] + u_[J.a1] + X0_[J.b1] + u_[J.b1]);
        Eigen::Vector2d Q = 0.5 * (X0_[J.a2] + u_[J.a2] + X0_[J.b2] + u_[J.b2]);
        Eigen::Vector2d e = Q - P;
        double L = e.norm();
        if (L < 1e-14) continue;
        Eigen::Vector2d n(e.y() / L, -e.x() / L);
        double a0 = (X0_[J.b1] + u_[J.b1] - X0_[J.a1] - u_[J.a1]).dot(n);
        double a1 = (X0_[J.b2] + u_[J.b2] - X0_[J.a2] - u_[J.a2]).dot(n);
        // une levre refermee ne rend pas de volume : on ne compte que l'ouvert
        vCrack += L * 0.5 * (std::max(0.0, a0) + std::max(0.0, a1)) * thk_;
    }
    hydroVol_ = vol + vCrack;
    hydroVolCrack_ = vCrack;
    // Le defaut de fermeture du lacet de la SOURCE (partie (a) ci-dessus).
    //
    // (!) CORRECTIF DU 2026-08-20 — L'ECHELLE DE REFERENCE ETAIT FAUSSE, et
    // l'alarme criait au loup sur des runs sains : run_hf_diag2 et
    // run_hfp_aniso la declenchaient avec les SEULES faces du forage, un
    // anneau ferme par construction, sans une seule levre de fissure.
    //
    // La raison : en FDEM les noeuds sont INTEGRALEMENT DEDOUBLES — 3 par
    // triangle, le log le dit (41 511 noeuds pour 13 837 triangles). Deux
    // faces consecutives de la paroi appartiennent a deux elements differents
    // et ne partagent AUCUN indice de noeud : la somme des segments orientes
    // ne mesure pas un defaut de fermeture topologique, elle mesure la somme
    // des OUVERTURES ELASTIQUES des joints de l'anneau. Elle ne peut donc pas
    // descendre au zero machine, et le seuil de 1e-6 rapporte au perimetre
    // etait SOUS le bruit physique (~1e-8 m par joint sous quelques MPa).
    //
    // La bonne echelle est la longueur MOYENNE D'UNE FACE : un contour
    // reellement ouvert laisse un trou de l'ordre d'une face (rapport ~ 1),
    // la respiration elastique des joints donne ~ 3e-5. Seuil a 0,05 : plus
    // d'un ordre de grandeur de marge des deux cotes.
    double lFace = (!hydroSrc_.empty()) ? perim / (double)hydroSrc_.size() : 0.0;
    hydroClose_ = (lFace > 0.0) ? closure.norm() / lFace : 0.0;
    if (hydroClose_ > 0.05 && !hydroCloseWarned_) {
        hydroCloseWarned_ = true;
        std::cout << "\n[HYDRO] *** CONTOUR SOURCE NON FERME *** defaut de "
                     "fermeture " << hydroClose_ << " fois la longueur d'une "
                     "face, a t = " << t_ << " s, " << hydroSrc_.size()
                  << " faces source.\n"
                     "[HYDRO] Le lacet ne mesure alors PAS une aire : le volume "
                     "de la cavite SOURCE, donc la pression, sont sans valeur. "
                     "Cause probable : la selection boreSelectR ne prend pas "
                     "tout le pourtour, ou l'anneau est coupe par une face "
                     "manquante. (Les levres de fissure ne passent PAS par ce "
                     "lacet : elles sont comptees localement, terme (b).)\n\n";
    }
}

// Pompe, pression, et chargement des levres.
void FdemSolver::hydroForces() {
    if (!hydroOn_) return;
    updateWetBoundary();
    const double tInj = t_ - hydroStart_;   // horloge de l injection (essai 3)
    if (tInj < 0.0) {
        // Pompe a l arret : phase geostatique + excavation du protocole de
        // l article. La masse SUIT le volume courant pour que p reste a p0
        // exactement — sans ce re-basage, la convergence de paroi pendant
        // l excavation comprimerait un fluide que personne n injecte encore.
        hydroP_ = hydroP0_;
        hydroMass_ = hydroVol_ * fluidRho0_;
    } else {
    double ramp = 1.0;
    if (hydroRamp_ > 0.0 && tInj < hydroRamp_)
        ramp = 0.5 * (1.0 - std::cos(M_PI * tInj / hydroRamp_));
    if (hydroRateMode_) {
        hydroMass_ += ramp * hydroRate_ * fluidRho0_ * dt_;
        // leur eq. 6 ; le plancher evite un log(<=0) si la cavite se ferme
        double r = hydroMass_ / std::max(hydroVol_ * fluidRho0_, 1e-300);
        hydroP_ = (r > 1e-300) ? hydroP0_ + fluidK_ * std::log(r) : hydroP0_;
    } else {
        hydroP_ = ramp * hydroPimp_;
        hydroMass_ = hydroVol_ * fluidRho0_
                   * std::exp((hydroP_ - hydroP0_) / fluidK_);  // pour la sortie
    }
    }
    // chargement : la pression pousse le solide a l'OPPOSE de la normale
    // sortante, moitie par noeud — meme forme que confiningForces(). Voir
    // l'en-tete au sujet de la coquille de leur eq. 7.
    //
    // (!) CORRECTIF DU 2026-08-20 — LE SIGNE ETAIT INVERSE. Cette boucle
    // appliquait +p n au lieu de -p n : le fluide SERRAIT la cavite au lieu
    // de l'ouvrir. Le forage produisait un breakout aligne sur sigma'_h et
    // rompait a 6,6 MPa au lieu des 12 MPa de leur eq. 10.
    //
    // La cause est une mauvaise lecture de leur eq. 7. Elle s'ecrit
    //     F_p12 = -(p/2) [y2-y1 ; x2-x1]
    // dont le vecteur n'est PAS orthogonal au segment (produit scalaire avec
    // (dx, dy) = 2 dx dy). La coquille est dans la SECONDE COMPOSANTE seule,
    // pas dans le signe de tete : la forme juste est
    //     F_p12 = -(p/2) [y2-y1 ; -(x2-x1)] = -(p/2) L n,   n = (dy,-dx)/L
    // Deux elements de l'article le tranchent : (a) leur section 3.1 pose
    // « negative sign to compressive stresses », donc une pression p > 0
    // impose sigma = -p I et la traction t = sigma.n = -p n ; (b) Y-Geo
    // triangule en CCW (Munjiza 2004 ; Mahabadi 2012 ; Lisjak 2014a), donc
    // leur (1 -> 2) est le parcours CCW et (dy, -dx) est bien sortant, comme
    // ici. Le moins de tete de leur eq. 7 etait donc PHYSIQUE ; c'est lui
    // qu'on avait supprime en croyant corriger la coquille.
    for (const auto& be : wetEdges_) {
        Eigen::Vector2d P = X0_[be.na] + u_[be.na];
        Eigen::Vector2d Q = X0_[be.nb] + u_[be.nb];
        Eigen::Vector2d d = Q - P;
        double L = d.norm();
        if (L < 1e-14) continue;
        Eigen::Vector2d n(d.y() / L, -d.x() / L);      // sortante du solide
        Eigen::Vector2d half = -0.5 * hydroP_ * L * thk_ * n;  // vers l'INTERIEUR du solide
        f_[be.na] += half;
        f_[be.nb] += half;
        hydroWork_ += dt_ * (half.dot(v_[be.na]) + half.dot(v_[be.nb]));
    }
}

// ---------------------------------------------------------------------------
// EXCAVATION (etude tunnel EDZ, Wang et al. 2024). Le massif est maille AVEC
// sa cavite et pre-contraint (insituSh / insituSv). Les faces de la paroi
// portent donc a t = 0 une traction DESEQUILIBREE : celle que la roche excavee
// exercait, t = sigma0 . n. On la retablit exactement (rel = 1 : equilibre
// parfait, rien ne bouge), puis on la fait decroitre jusqu'a zero. C'est la
// methode convergence-confinement ; elle remplace le "core modulus reduction"
// de l'article (meme etat initial, meme etat final) sans avoir a faire varier
// un module en cours de calcul.
//
// La traction est evaluee FACE PAR FACE comme sigma0 . n : exact pour un champ
// anisotrope (lambda != 1), ce qu'une pression scalaire ne saurait pas faire.
// Convention de normale et de signe identiques a confiningForces() : n sortant
// du solide, et sigma0 . n avec sigma0 = -p I redonne exactement -p n.
// ---------------------------------------------------------------------------
void FdemSolver::setupExcavation() {
    excRelease_ = cfg_.getb("excavRelease", false);
    if (!excRelease_) return;
    if (!hasInsitu_)
        throw std::runtime_error("excavRelease exige une contrainte in situ "
                                 "(insituSh / insituSv) : sans elle il n'y a "
                                 "rien a relacher");
    excStart_ = cfg_.getd("excavStart", 0.0);
    excRamp_ = cfg_.getd("excavRamp", 0.0);
    if (excRamp_ <= 0.0)
        throw std::runtime_error("excavRamp doit etre > 0 : un relachement "
                                 "instantane lance une onde de choc dans le "
                                 "massif (meme piege que confiningRamp = 0). "
                                 "Compter plusieurs transits d'onde.");
    double bcx = cfg_.getd("boreCX", 0.5 * W_);
    double bcy = cfg_.getd("boreCY", 0.5 * H_);
    double bsr = cfg_.getd("boreSelectR", 0.0);
    if (bsr <= 0.0)
        throw std::runtime_error("excavRelease exige boreSelectR > 0 : rayon "
                                 "de selection des faces de paroi autour de "
                                 "(boreCX, boreCY)");
    for (const auto& be : exterior_) {
        Eigen::Vector2d M = 0.5 * (X0_[be.na] + X0_[be.nb]);
        double dx = M.x() - bcx, dy = M.y() - bcy;
        if (dx * dx + dy * dy <= bsr * bsr) excEdges_.push_back(be);
    }
    if (excEdges_.empty())
        throw std::runtime_error("excavRelease : aucune face de cavite "
                                 "selectionnee — verifier boreCX/boreCY/"
                                 "boreSelectR contre la geometrie du maillage");
    double len = 0.0;
    for (const auto& be : excEdges_)
        len += (X0_[be.nb] - X0_[be.na]).norm();
    std::cout << "[FDEM] excavation : " << excEdges_.size()
              << " faces de paroi, perimetre " << len << " m, relachement de "
              << excStart_ << " s a " << (excStart_ + excRamp_) << " s\n";
}

double FdemSolver::excavRelief() const {
    if (!excRelease_) return 0.0;
    if (t_ <= excStart_) return 1.0;
    double x = (t_ - excStart_) / excRamp_;
    if (x >= 1.0) return 0.0;
    return 0.5 * (1.0 + std::cos(M_PI * x));
}

void FdemSolver::excavationForces() {
    if (!excRelease_) return;
    double rel = excavRelief();
    if (rel <= 0.0) return;
    for (const auto& be : excEdges_) {
        Eigen::Vector2d P = X0_[be.na] + u_[be.na];
        Eigen::Vector2d Q = X0_[be.nb] + u_[be.nb];
        Eigen::Vector2d d = Q - P;
        double L = d.norm();
        if (L < 1e-14) continue;
        Eigen::Vector2d n(d.y() / L, -d.x() / L);      // sortante du solide
        // traction que la roche excavee exercait sur la paroi : t = sigma0 . n
        Eigen::Vector2d half = 0.5 * rel * L * thk_ * (insituS_ * n);
        f_[be.na] += half;
        f_[be.nb] += half;
        // B4 : verse dans le canal confinement (deja cable au budget)
        confWork_ += dt_ * (half.dot(v_[be.na]) + half.dot(v_[be.nb]));
    }
}

// Mean sigma_xx over the CORE of the specimen (middle half in both
// directions), i.e. away from the loaded faces where the pressure has not yet
// diffused into a uniform state. For a working confinement this must come out
// at -confiningPressure.
double FdemSolver::achievedConfinement() const {
    double sum = 0.0, area = 0.0;
    for (const auto& e : el_) {
        Eigen::Vector2d c = (X0_[e.n[0]] + X0_[e.n[1]] + X0_[e.n[2]]) / 3.0;
        if (std::abs(c.x() - 0.5 * W_) > 0.25 * W_) continue;
        if (std::abs(c.y() - 0.5 * H_) > 0.25 * H_) continue;
        sum += e.sxx * e.A0;
        area += e.A0;
    }
    return area > 0.0 ? sum / area : 0.0;
}

// Area-averaged stress over the core of the disc (r <= 0.15 R), the zone the
// classical solution describes as uniform enough to call "the centre".
void FdemSolver::discCentreStress(double& sxx, double& syy) const {
    double sx = 0.0, sy = 0.0, area = 0.0;
    double r2 = 0.15 * discR_ * 0.15 * discR_;
    for (const auto& e : el_) {
        Eigen::Vector2d c = (X0_[e.n[0]] + X0_[e.n[1]] + X0_[e.n[2]]) / 3.0;
        if ((c - discC_).squaredNorm() > r2) continue;
        sx += e.sxx * e.A0;
        sy += e.syy * e.A0;
        area += e.A0;
    }
    sxx = area > 0.0 ? sx / area : 0.0;
    syy = area > 0.0 ? sy / area : 0.0;
}

void FdemSolver::integrate() {
    // V2/B4 : KE au premier pas (vitesses initiales intactes), compteurs purs
    if (keInit_ < 0.0) {
        double ke0 = 0.0;
        for (std::size_t i = 0; i < X0_.size(); ++i)
            ke0 += 0.5 * m_[i] * v_[i].squaredNorm();
        keInit_ = ke0;
    }
    double cw = 0.0, lw = 0.0, bw = 0.0, bias = 0.0;
    if (adaptive_) {
        // Bound groups integrate as ONE node: forces and masses summed,
        // Cundall damping and the quiet-boundary terms applied to the sums,
        // the common velocity written back to every copy. For a singleton
        // group this reduces exactly to the per-node path below. Groups never
        // span vertices, so the vertex loop parallelizes cleanly.
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:cw,lw,bw,bias)
#endif
        for (int vv = 0; vv < nVert_; ++vv) {
            for (const auto& g : grpsOfVert_[vv]) {
                int i0 = g[0];                         // copies share flags
                if (flag_[i0] == FIXED) {
                    if (gripFree_) {
                        double F = 0.0, M = 0.0;
                        for (int i : g) { F += f_[i].x(); M += m_[i]; }
                        double vx = v_[i0].x() + (dt_ / M) * F;
                        for (int i : g) {
                            v_[i] = {vx, 0.0};
                            u_[i] += dt_ * v_[i];
                        }
                    } else for (int i : g) v_[i].setZero();
                    continue;
                }
                if (flag_[i0] == DRIVEX) {         // SHPB struck face
                    double F = 0.0, M = 0.0;
                    for (int i : g) { F += f_[i].y(); M += m_[i]; }
                    double vy = v_[i0].y() + (dt_ / M) * F;
                    for (int i : g) {
                        v_[i] = {vDrive_, vy};
                        u_[i] += dt_ * v_[i];
                    }
                    continue;
                }
                if (flag_[i0] == PRESCRIBED) {
                    double vg = pullV_;
                    double tl = t_ - pullDelay_;   // voir platenForces()
                    if (tl <= 0.0) vg = 0.0;
                    else if (pullRamp_ > 0.0 && tl < pullRamp_)
                        vg *= 0.5 * (1.0 - std::cos(M_PI * tl / pullRamp_));
                    double vx = 0.0;
                    if (gripFree_) {
                        double F = 0.0, M = 0.0;
                        for (int i : g) { F += f_[i].x(); M += m_[i]; }
                        vx = v_[i0].x() + (dt_ / M) * F;
                    }
                    for (int i : g) {
                        // V2/B4 (corrige 2026-08-14, formule Sierra) : les
                        // familles comptent deja f.v sur ce noeud ; le mors
                        // ajoute le travail de la LIAISON R = m.a - f, integre
                        // en vitesse moyenne (trapeze, coherent leapfrog)
                        double vgOld = v_[i].y();
                        double Ry = m_[i] * (vg - vgOld) / dt_ - f_[i].y();
                        bw += Ry * 0.5 * (vg + vgOld) * dt_;
                        v_[i] = {vx, vg};
                        u_[i] += dt_ * v_[i];
                    }
                    continue;
                }
                Eigen::Vector2d F = Eigen::Vector2d::Zero();
                double M = 0.0, cX = 0.0, cY = 0.0;
                for (int i : g) {
                    F += f_[i];
                    if (kAbsX_[i] > 0) {
                        double fk = kAbsX_[i] * u_[i].x();
                        F.x() -= fk;
                        lw -= fk * v_[i0].x() * dt_;   // V2/B4 ressort
                    }
                    if (kAbsY_[i] > 0) {
                        double fk = kAbsY_[i] * u_[i].y();
                        F.y() -= fk;
                        lw -= fk * v_[i0].y() * dt_;
                    }
                    M += m_[i];
                    cX += cAbsX_[i];
                    cY += cAbsY_[i];
                }
                if (damping_ > 0) {
                    double fdx = damping_ * std::abs(F.x())
                             * (v_[i0].x() > 0 ? 1.0 : (v_[i0].x() < 0 ? -1.0 : 0.0));
                    F.x() -= fdx;
                    cw -= fdx * v_[i0].x() * dt_;      // V2/B4 Cundall
                    double fdy = damping_ * std::abs(F.y())
                             * (v_[i0].y() > 0 ? 1.0 : (v_[i0].y() < 0 ? -1.0 : 0.0));
                    F.y() -= fdy;
                    cw -= fdy * v_[i0].y() * dt_;
                }
                bias += F.squaredNorm() * dt_ * dt_ / (2.0 * M);
                Eigen::Vector2d vn = v_[i0] + (dt_ / M) * F;
                if (cX > 0) {
                    vn.x() /= 1.0 + dt_ * cX / M;
                    lw -= cX * vn.x() * vn.x() * dt_;  // V2/B4 amortisseur
                }
                if (cY > 0) {
                    vn.y() /= 1.0 + dt_ * cY / M;
                    lw -= cY * vn.y() * vn.y() * dt_;
                }
                if (flag_[i0] == ROLLERX) vn.x() = 0.0;
                for (int i : g) {
                    if (flag_[i0] == ROLLERX) u_[i].x() = 0.0;
                    v_[i] = vn;
                    u_[i] += dt_ * vn;
                }
            }
        }
        cundWork_ += cw;                   // V2/B4
        lysWork_ += lw;
        bcWork_ += bw;
        biasW_ += bias;
        if (scen_ == Scenario::BRAZILIAN || (scen_ == Scenario::TENSION
                                             && tensionPlatens_)) {
            if (brazPlatens_ || tensionPlatens_) {
                plTop_.y += dt_ * plTop_.v;
                plBot_.y += dt_ * plBot_.v;
            }
        } else if (scen_ != Scenario::TENSION) {
            tool_.integrate(dt_);
        }
        return;
    }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(+:cw,lw,bw,bias)
#endif
    for (int i = 0; i < (int)X0_.size(); ++i) {
        if (flag_[i] == FIXED) {
            // gripFree: frictionless grips — only the axial dof is held, the
            // lateral one follows the forces. Removes the Saint-Venant
            // corner concentration of fully clamped grips, which otherwise
            // decides where a tension specimen breaks (it drowned the
            // random-strength-field localization in the wbm experiment).
            if (gripFree_) {
                v_[i].x() += (dt_ / m_[i]) * f_[i].x();
                v_[i].y() = 0.0;
                u_[i] += dt_ * v_[i];
            } else v_[i].setZero();
            continue;
        }
        if (flag_[i] == DRIVEX) {                      // SHPB struck face
            // v_x is imposed by the pulse; the LATERAL dof integrates freely
            // (a bar end is not clamped: clamping it would radiate a shear
            // wave from the driven face at every step)
            v_[i].y() += (dt_ / m_[i]) * f_[i].y();
            v_[i].x() = vDrive_;
            u_[i] += dt_ * v_[i];
            continue;
        }
        if (flag_[i] == PRESCRIBED) {
            if (gripFree_) v_[i].x() += (dt_ / m_[i]) * f_[i].x();
            else v_[i].x() = 0.0;
            // pullRamp: smooth (cosine) rise of the grip velocity over the
            // given time instead of a step. A stepped grip launches a stress
            // transient that localizes failure at the FIRST joint row under
            // the grip whatever the strength map says (measured: straight
            // "unzipping" of the top row); ramping over many wave transits
            // makes the loading quasi-static so the specimen breaks where it
            // is WEAK, not where it is pulled.
            double vg = pullV_;
            if (pullRamp_ > 0.0 && t_ < pullRamp_)
                vg *= 0.5 * (1.0 - std::cos(M_PI * t_ / pullRamp_));
            {   // V2/B4 (corrige 2026-08-14, formule Sierra) : travail de la
                // LIAISON R = m.a - f en vitesse moyenne — les familles
                // comptent deja f.v sur ce noeud (l'ancien +f.vg avait le
                // signe inverse ET ignorait l'inertie de la rampe)
                double vgOld = v_[i].y();
                double Ry = m_[i] * (vg - vgOld) / dt_ - f_[i].y();
                bw += Ry * 0.5 * (vg + vgOld) * dt_;
            }
            v_[i].y() = vg;
            u_[i] += dt_ * v_[i];
            continue;
        }
        if (kAbsX_[i] > 0 || kAbsY_[i] > 0) {          // boundary springs
            double fkx = kAbsX_[i] * u_[i].x();
            double fky = kAbsY_[i] * u_[i].y();
            f_[i].x() -= fkx;
            f_[i].y() -= fky;
            lw -= (fkx * v_[i].x() + fky * v_[i].y()) * dt_;  // V2/B4
        }
        if (damping_ > 0) {                            // Cundall local damping
            double fdx = damping_ * std::abs(f_[i].x())
                         * (v_[i].x() > 0 ? 1.0 : (v_[i].x() < 0 ? -1.0 : 0.0));
            f_[i].x() -= fdx;
            double fdy = damping_ * std::abs(f_[i].y())
                         * (v_[i].y() > 0 ? 1.0 : (v_[i].y() < 0 ? -1.0 : 0.0));
            f_[i].y() -= fdy;
            cw -= (fdx * v_[i].x() + fdy * v_[i].y()) * dt_;  // V2/B4
        }
        bias += f_[i].squaredNorm() * dt_ * dt_ / (2.0 * m_[i]);
        v_[i] += (dt_ / m_[i]) * f_[i];
        if (cAbsX_[i] > 0) {
            v_[i].x() /= 1.0 + dt_ * cAbsX_[i] / m_[i];
            lw -= cAbsX_[i] * v_[i].x() * v_[i].x() * dt_;    // V2/B4
        }
        if (cAbsY_[i] > 0) {
            v_[i].y() /= 1.0 + dt_ * cAbsY_[i] / m_[i];
            lw -= cAbsY_[i] * v_[i].y() * v_[i].y() * dt_;
        }
        if (flag_[i] == ROLLERX) { v_[i].x() = 0.0; u_[i].x() = 0.0; }
        u_[i] += dt_ * v_[i];
    }
    cundWork_ += cw;                       // V2/B4
    lysWork_ += lw;
    bcWork_ += bw;
    biasW_ += bias;
    if (scen_ == Scenario::BRAZILIAN || (scen_ == Scenario::TENSION
                                         && tensionPlatens_)) {
        if (brazPlatens_ || tensionPlatens_) {         // kinematic planes
            plTop_.y += dt_ * plTop_.v;                // v ramped in
            plBot_.y += dt_ * plBot_.v;                // platenForces()
        }
    } else if (scen_ != Scenario::TENSION) {
        tool_.integrate(dt_);
    }
}

// Fragments: connected components of elements over still-cohesive joints.
void FdemSolver::computeFragments() {
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
            mm += rhoP_[el_[e].phase] * el_[e].A0 * thk_;
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
    double vDet = 0;                                   // volume, phase-neutral
    for (int e = 0; e < (int)el_.size(); ++e)
        if (fragId_[e] != 0) vDet += el_[e].A0 * thk_;
    detachedVol_ = vDet;
}

// ===========================================================================
// Output
// ===========================================================================

void FdemSolver::writeFrame(int frame) {
    computeFragments();

    std::vector<Eigen::Vector2d> pts(X0_.size());
    std::vector<Eigen::Vector2d> vel(X0_.size());
    for (std::size_t i = 0; i < X0_.size(); ++i) {
        pts[i] = X0_[i] + u_[i];
        vel[i] = v_[i];
    }
    std::vector<std::array<int, 3>> tris(el_.size());
    std::vector<double> svm(el_.size()), frag(el_.size());
    std::vector<double> phs(el_.size()), grn(el_.size());
    std::vector<double> sxx(el_.size()), syy(el_.size());
    // sigmaXY : le cisaillement manquait dans les .vtu du FDEM 2D alors que
    // el_[e].sxy est deja calcule par elementForces (l'insertion adaptative
    // s'en sert). Sans lui, ni les contraintes principales ni l'orientation
    // des bandes de cisaillement ne sont lisibles en post-traitement.
    std::vector<double> sxy(el_.size()), exx(el_.size());
    for (std::size_t e = 0; e < el_.size(); ++e) {
        tris[e] = el_[e].n;
        svm[e] = el_[e].svm;
        frag[e] = fragId_[e];
        phs[e] = el_[e].phase;
        grn[e] = el_[e].grain;
        sxx[e] = el_[e].sxx;
        syy[e] = el_[e].syy;
        sxy[e] = el_[e].sxy;                           // shear (was missing)
        exx[e] = el_[e].exx;                           // axial strain
    }
    char name[64];
    std::snprintf(name, sizeof(name), "/fdem_%04d.vtu", frame);
    if (bdOn_) {
        std::vector<double> bdv(el_.size());
        for (std::size_t e = 0; e < el_.size(); ++e) bdv[e] = el_[e].bdD;
        vtk::writeTriMesh(out_ + name, pts, tris,
                          {{"vonMises", &svm}, {"fragment", &frag},
                           {"phase", &phs}, {"grain", &grn},
                           {"sigmaXX", &sxx}, {"sigmaYY", &syy},
                           {"sigmaXY", &sxy}, {"epsXX", &exx},
                           {"bulkD", &bdv}},
                          {{"velocity", &vel}});
    } else
    vtk::writeTriMesh(out_ + name, pts, tris,
                      {{"vonMises", &svm}, {"fragment", &frag},
                       {"phase", &phs}, {"grain", &grn},
                       {"sigmaXX", &sxx}, {"sigmaYY", &syy},
                       // sigmaXY completes the in-plane tensor the solver
                       // already carries per element (Elem::sxy, set in
                       // elementStress): without it a post-processor cannot
                       // form the principal stresses or the shear field.
                       {"sigmaXY", &sxy}, {"epsXX", &exx}},
                      {{"velocity", &vel}});

    std::vector<std::array<int, 2>> lines;
    std::vector<double> Dj, tb, Tp, Fs, Bd, Fm, Bm, Dt, Ed;
    for (const auto& J : jt_) {
        lines.push_back({J.a1, J.a2});
        Dj.push_back(J.D);
        tb.push_back(J.tBreak);
        Tp.push_back(J.type);
        Fs.push_back(J.stat);
        Bd.push_back(J.bonded ? 1.0 : 0.0);
        Fm.push_back(J.failMode);
        Bm.push_back(J.bmode);
        Dt.push_back(J.difT);          // E8 : voir Fdem3dSolver::writeFrame
        Ed.push_back(J.edotIns);
    }
    std::snprintf(name, sizeof(name), "/fdem_joints_%04d.vtu", frame);
    vtk::ScalarField jf{
        {"damage", &Dj}, {"tBreak", &tb}, {"type", &Tp},
        {"ftScale", &Fs}, {"bonded", &Bd}, {"breakMode", &Bm}};
    if (difOn_) { jf["difT"] = &Dt; jf["edotIns"] = &Ed; }
    if (cfg_.getb("writeJointMode", false)) jf["failMode"] = &Fm;
    vtk::writeLines(out_ + name, pts, lines, jf);

    std::ofstream fm(out_ + "/frames.csv",
                     frame == 0 ? std::ios::trunc : std::ios::app);
    if (frame == 0) fm << "frame,t,toolX,toolY\n";
    // the brazilian has no tool: the platen plane goes in the toolY column so
    // make_gif.py and the GUI keep reading a meaningful loading position
    bool noTool = scen_ == Scenario::BRAZILIAN || scen_ == Scenario::SHPB;
    double tx = noTool ? discC_.x() : tool_.x.x();
    double ty = noTool ? discC_.y() + discR_ : tool_.x.y();
    fm << frame << "," << t_ << "," << tx << "," << ty << "\n";
}

void FdemSolver::historyHeader(std::ostream& os) const {
    if (scen_ == Scenario::SHPB) {
        // epsM1 / epsM2 = area-averaged axial strain at monitor points 1 and 2
        // (fig. 23); vDrive = the prescribed pulse; sxxC/syyC = the stress at
        // the centre of the disc (fig. 25b reads syyC); nInserted = joints
        // activated so far (adaptive only).
        os << "t,vDrive,epsM1,epsM2,sxxC,syyC,nBroken,nFrag,nInserted\n";
        return;
    }
    // (ancienne ligne unique supprimee : elle court-circuitait le bloc
    // ci-dessous — doublon laisse par la fusion a trois voies)
    if (scen_ == Scenario::TENSION) {
        // compression par plateaux : on ajoute la metrologie de l'essai —
         // trois deformations axiales (machine, faces, extensometre), le
        // decompte des fissures par MODE (traction / cisaillement), la
        // pression laterale effectivement atteinte et le marqueur de fin
        // d'essai peakLocked. Le montage a grips garde EXACTEMENT ses
        // anciennes colonnes.
        if (tensionPlatens_)
            os << "t,gripFy,sigma,sigmaPeak,nBroken,epsPlaten,epsSpec,"
                  "epsGauge,nBrokTen,nBrokShear,nFrag,confAchieved,"
                  "peakLocked";
        else
            os << "t,gripFy,sigma,sigmaPeak,nBroken";
        if (hydroOn_) os << ",hydroP,hydroVol,hydroMass,hydroNWet,eHydro";
        // Essai 0 : sans ces deux compteurs, la pression d insertion
        // n est connue qu a l espacement des trames pres. nBroken ne
        // compte que D >= 1, soit un evenement TARDIF : sur le bench
        // AbuAisha l insertion precede la premiere rupture de ~300 us.
        if (adaptive_) os << ",nInserted,nDamaging";
        if (bdOn_) os << ",nPulv,bdWork";
        os << "\n";
        return;
    }
    if (scen_ == Scenario::BRAZILIAN) {
        // P = vertical force applied on the top arc, Pbot on the bottom one;
        // their SUM is the unbalanced force that makes the disc drift and must
        // stay small next to P (mesh asymmetry between the two arcs)
        // sxxC/syyC = area-averaged stress over the core of the disc, the
        // quantity the elastic gauge reads; carrying them at every history row
        // (instead of only at the last intact step) is what lets the gauge be
        // read in the genuinely elastic part of the loading. peakLocked = the
        // end-of-test marker: 0 while the test is a valid indirect tension
        // measurement, 1 once the post-peak drop has been seen and the platens
        // are only crushing the halves together. Truncate the published curve
        // at the first row carrying 1.
        os << "t,P,Pbot,drive,sigmaT,sigmaTpeak,nBroken,nFrag,sxxC,syyC,"
              "peakLocked\n";
        return;
    }
    os << "t,toolFx,toolFy,toolX,toolY,toolVx,toolVy,work,toolKE,"
          "nBroken,nFrag,detachedVol,specificEnergy"
          ",eEl,eJnt,eGc,eFric,eCund,eLys";      // V2/B4
    // FR-010 : une cavite qui se remplit sans qu on puisse la regarder
    // est ingouvernable. hydroP est LA courbe de leur fig. 11.
    if (hydroOn_) os << ",hydroP,hydroVol,hydroMass,hydroNWet,eHydro";
    if (adaptive_) os << ",nInserted,nDamaging";
    if (bdOn_) os << ",nPulv,bdWork";
    os << "\n";
}

// Essai 0 (2026-08-20). Deux compteurs d etat des interfaces, pour
// l historique. nInserted = interfaces effectivement creees par le
// critere extrinseque (bonded = false) ; nDamaging = interfaces en
// cours d adoucissement, 0 < D < 1, celles qui travaillent SANS
// recevoir la pression de fluide (updateWetBoundary les exclut).
// Balayage O(nJoints) appele une fois par ligne d historique, soit
// ~2000 fois par run : negligeable devant le meme balayage effectue
// a chaque pas par jointForces().
void FdemSolver::countInserted(long& nIns, long& nDam) const {
    nIns = 0; nDam = 0;
    for (const auto& J : jt_) {
        if (J.bonded) continue;
        ++nIns;
        if (J.D > 0.0 && J.D < 1.0) ++nDam;
    }
}

void FdemSolver::historyRow(std::ostream& os) const {
    if (scen_ == Scenario::SHPB) {
        double sx = 0.0, sy = 0.0;
        if (!shpbNoDisc_) discCentreStress(sx, sy);
        os << t_ << "," << vDrive_ << "," << epsM1_ << "," << epsM2_ << ","
           << sx << "," << sy << "," << nBroken_ << "," << nFrag_ << ","
           << nInserted_ << "\n";
        return;
    }
    if (scen_ == Scenario::TENSION) {
        os << t_ << "," << gripF_.y() << ","
           << std::abs(gripF_.y()) / (W_ * thk_) << ","
           << sigmaPeak_ << "," << nBroken_;
        if (tensionPlatens_) {
            double ep = 0.0, es = 0.0, eg = 0.0;
            gaugeStrain(ep, es, eg);
            long nT = 0, nS = 0;
            for (const auto& J : jt_) {
                if (J.bmode == 1) ++nT;
                else if (J.bmode == 2) ++nS;
            }
            os << "," << ep << "," << es << "," << eg << "," << nT << ","
               << nS << "," << nFrag_ << "," << confAchieved_ << ","
               << (peakLockedU_ ? 1 : 0);
        }
        if (hydroOn_)
            os << "," << hydroP_ << "," << hydroVol_ << "," << hydroMass_
               << "," << hydroNWet_ << "," << hydroWork_;
        if (adaptive_) { long ni, nd; countInserted(ni, nd);
                         os << "," << ni << "," << nd; }
        if (bdOn_) os << "," << nPulv_ << "," << bdWork_;
        os << "\n";
        return;
    }
    if (scen_ == Scenario::BRAZILIAN) {
        double ft = brazPlatens_ ? plTop_.F.y() : arcTop_.F.y();
        double fb = brazPlatens_ ? plBot_.F.y() : arcBot_.F.y();
        double drive = brazPlatens_ ? plTop_.y : brazP_;
        double sxxC = 0.0, syyC = 0.0;
        discCentreStress(sxxC, syyC);
        os << t_ << "," << ft << "," << fb << "," << drive << ","
           << sigmaT_ << "," << sigmaTpeak_ << "," << nBroken_ << ","
           << nFrag_ << "," << sxxC << "," << syyC << ","
           << (peakLocked_ ? 1 : 0) << "\n";
        return;
    }
    double Es = detachedVol_ > 0 ? work_ / detachedVol_ : 0.0;
    os << t_ << "," << tool_.F.x() << "," << tool_.F.y() << ","
       << tool_.x.x() << "," << tool_.x.y() << ","
       << tool_.v.x() << "," << tool_.v.y() << ","
       << work_ << "," << tool_.ke() << "," << nBroken_ << "," << nFrag_ << ","
       << detachedVol_ << "," << Es
       << "," << elWork_ << "," << jointWork_ << "," << gcWork_ << ","
       << gcFricWork_ << "," << cundWork_ << "," << lysWork_;   // B4
    if (hydroOn_)
        os << "," << hydroP_ << "," << hydroVol_ << "," << hydroMass_
           << "," << hydroNWet_ << "," << hydroWork_;
    if (adaptive_) { long ni, nd; countInserted(ni, nd);
                     os << "," << ni << "," << nd; }
    if (bdOn_) os << "," << nPulv_ << "," << bdWork_;
    os << "\n";
}

void FdemSolver::finalize() {
    computeFragments();

    std::ofstream fe(out_ + "/fdem_final_elements.csv");
    fe << "cx,cy,fragment,phase,grain\n";
    for (std::size_t e = 0; e < el_.size(); ++e) {
        Eigen::Vector2d c = Eigen::Vector2d::Zero();
        for (int a = 0; a < 3; ++a) c += (X0_[el_[e].n[a]] + u_[el_[e].n[a]]);
        c /= 3.0;
        fe << c.x() << "," << c.y() << "," << fragId_[e] << ","
           << el_[e].phase << "," << el_[e].grain << "\n";
    }
    // Body-force verification output (Yan et al. section 3.1). Only written
    // when a body force is active, so no existing run gains a file.
    if (gravity_ > 0.0) {
        std::ofstream fg(out_ + "/fdem_nodal_displacement.csv");
        fg << "x0,y0,ux,uy\n";
        for (std::size_t i = 0; i < X0_.size(); ++i)
            fg << X0_[i].x() << "," << X0_[i].y() << ","
               << u_[i].x() << "," << u_[i].y() << "\n";
        // Mean settlement of the TOP row (y = H): the quantity the closed-form
        // solution of a self-weighted confined column predicts.
        double s = 0.0; long n = 0;
        for (std::size_t i = 0; i < X0_.size(); ++i)
            if (X0_[i].y() > H_ - 1e-9) { s += u_[i].y(); ++n; }
        std::cout << "[FDEM] gravity g = " << gravity_ << " m/s^2, mean top "
                     "displacement uy = " << (n ? s / n : 0.0) << " m over "
                  << n << " nodes\n";
    }

    std::ofstream fj(out_ + "/fdem_final_joints.csv");
    // breakMode / rn / rs ajoutes en QUEUE de ligne : les colonnes
    // existantes gardent leur place et leur ordre.
    fj << "x1,y1,x2,y2,damage,type,breakMode,rn,rs,tBreak,bonded\n";
    for (const auto& J : jt_) {
        Eigen::Vector2d P = X0_[J.a1] + u_[J.a1], Q = X0_[J.a2] + u_[J.a2];
        fj << P.x() << "," << P.y() << "," << Q.x() << "," << Q.y() << ","
           << J.D << "," << J.type << "," << J.bmode << "," << J.rnB << ","
           << J.rsB << "," << J.tBreak << "," << (J.bonded ? 1 : 0) << "\n";
    }

    double keBlock = 0.0;
    for (std::size_t i = 0; i < X0_.size(); ++i)
        keBlock += 0.5 * m_[i] * v_[i].squaredNorm();
    std::cout << "\n[FDEM] ---- summary ----\n"
              << "[FDEM] block kinetic energy at end: " << keBlock << " J/m\n"
              << "[FDEM] net work injected by general contact: " << gcWork_
              << " J/m\n";
    if (contactPot_)
        std::cout << "[FDEM] (contact = potential : champ conservatif — un "
                     "petit residu positif est le biais O(dt) du compteur "
                     "plus la releve de naissance, pas une pathologie ; en "
                     "mode penalty tout positif est une injection)\n";
    // A dashpot can only DISSIPATE. A positive figure here means the viscous
    // branch is injecting energy — the rectifier failure mode — and every
    // number the run produced is suspect. One multiply per integration point.
    std::cout << "[FDEM] joint dashpot work: " << dampWork_ << " J/m  ["
              << (dampWork_ <= 0.0 ? "OK, dissipative"
                                   : "FAIL - the dashpot INJECTED energy")
              << "]\n";
    {   // ---- V2/B4 : bilan d'energie par sous-systeme (miroir du 3D) -----
        // KE(t) - KE(0) = somme des travaux par famille + residu O(dt).
        // U_el en invariants isotropes (plane strain, szz = nu (sxx+syy)) ;
        // E de la phase retrouve de Dm(2,2) = G — exact en elastique, approx
        // sous law/caps. Toutes les grandeurs 2D sont PAR METRE d'epaisseur
        // (thk_ inclus), comme le reste du resume.
        double uEl = 0.0;
        for (const auto& e : el_) {
            double nu = nuP_[e.phase];
            double Eph = 2.0 * (1.0 + nu) * DmP_[e.phase](2, 2);
            double szz = nu * (e.sxx + e.syy);
            double ss = e.sxx * e.sxx + e.syy * e.syy + szz * szz
                      + 2.0 * e.sxy * e.sxy;
            double trs = e.sxx + e.syy + szz;
            uEl += e.A0 * thk_ * ((1.0 + nu) * ss - nu * trs * trs)
                   / (2.0 * Eph);
        }
        double uSpr = 0.0;                 // ressorts absorbants (stocke)
        for (std::size_t i = 0; i < X0_.size(); ++i) {
            if (kAbsX_[i] > 0) uSpr += 0.5 * kAbsX_[i] * u_[i].x() * u_[i].x();
            if (kAbsY_[i] > 0) uSpr += 0.5 * kAbsY_[i] * u_[i].y() * u_[i].y();
        }
        double sumW = elWork_ + jointWork_ + gcWork_ + cundWork_ + lysWork_
                    + toolWork_ + bcWork_ + confWork_ + hydroWork_ + biasW_;
        double dKE = keBlock - keInit_;
        double resid = dKE - sumW;
        // echelle du verdict : le flux BRUT echange (la somme signee est ~0
        // precisement quand le bilan boucle — premier faux CHECK mesure sur
        // la percussion outil rigide : residu 0.014 J sur 1.53 J injectes,
        // soit 0.9 %, juge a 147 % de la somme signee)
        double gross = std::abs(elWork_) + std::abs(jointWork_)
                     + std::abs(gcWork_) + std::abs(cundWork_)
                     + std::abs(lysWork_) + std::abs(toolWork_)
                     + std::abs(bcWork_) + std::abs(confWork_)
                     + std::abs(hydroWork_);
        double scale = std::max({keInit_, keBlock, gross, 1e-30});
        bool zeroCase = scale < 1e-12;     // charge nulle : verdict en absolu
        std::cout << "[FDEM] energy budget (V2/B4): KE " << keInit_ << " -> "
                  << keBlock << " J/m\n"
                  << "[FDEM]   elements     : " << -elWork_
                  << " J/m preleve (stocke elastique " << uEl << " J/m)\n";
        if (bulkVisc_ > 0.0)
            std::cout << "[FDEM]      dont visqueux (2 mu D) : "
                      << -viscWork_ << " J/m dissipes, soit "
                      << (std::abs(elWork_) > 1e-300
                          ? 100.0 * viscWork_ / elWork_ : 0.0)
                      << " % du poste elements  "
                      << (viscWork_ <= 0.0
                          ? "[OK, dissipatif]"
                          : "[FAIL - le terme visqueux a INJECTE de l energie]")
                      << ". VENTILATION : deja comptee dans la ligne "
                         "ci-dessus, pas un poste de plus.\n";
        if (bdOn_)
            std::cout << "[FDEM]      dont pulverisation (bulkDamage) : "
                      << -bdWork_ << " J/m dissipes, " << nPulv_
                      << " elements a D = Dmax. VENTILATION : deja "
                         "comptee dans le poste elements." << std::endl;
        std::cout << "[FDEM]   joints       : " << -(jointWork_ - dampWork_)
                  << " J/m cohesif (fissuration + stocke), dashpot "
                  << -dampWork_ << " J/m\n"
                  << "[FDEM]   contact      : " << -gcWork_
                  << " J/m (dont frottement " << -gcFricWork_ << " J/m)\n"
                  << "[FDEM]   Cundall      : " << -cundWork_ << " J/m\n"
                  << "[FDEM]   frontieres   : " << -lysWork_
                  << " J/m (dont stocke ressorts " << uSpr << " J/m)\n"
                  << "[FDEM]   outil->solide: " << toolWork_
                  << " J/m, platines/grips: " << bcWork_ << " J/m\n";
        if (confP_ > 0.0)                  // sortie inchangee si pas confine
            std::cout << "[FDEM]   confinement  : " << confWork_
                      << " J/m (pression suiveuse -> solide)\n";
        std::cout << "[FDEM]   integration  : +" << biasW_
                  << " J/m (correction leapfrog f^2 dt^2/2m)\n"
                  << "[FDEM]   residu       : " << resid << " J/m ("
                  << 100.0 * std::abs(resid) / scale << " % de l'echelle) ["
                  << (zeroCase ? "OK (zero machine)"
                      : std::abs(resid) <= 0.01 * scale ? "OK" : "CHECK")
                  << "]\n";
    }
    if (adaptive_)
        std::cout << "[FDEM] adaptive insertion: " << nInserted_ << " / "
                  << jt_.size() << " joints inserted ("
                  << (jt_.empty() ? 0.0 : 100.0 * nInserted_ / jt_.size())
                  << " %), " << nBroken_ << " fully broken\n";
    if (difOn_) {
        // Sur les joints REELLEMENT inseres : c est la seule population sur
        // laquelle le DIF a ete evalue.
        std::vector<double> er, dt_v;
        for (const auto& J : jt_)
            if (!J.bonded) { er.push_back(J.edotIns); dt_v.push_back(J.difT); }
        if (er.empty()) {
            std::cout << "[FDEM] strainRateDIF : aucun joint insere — le DIF "
                         "n a jamais ete evalue.\n";
        } else {
            std::sort(er.begin(), er.end());
            std::sort(dt_v.begin(), dt_v.end());
            std::size_t n = er.size();
            std::cout << "[FDEM] strainRateDIF (sur " << n
                      << " joints inseres) : taux de deformation a "
                         "l insertion, mediane " << er[n / 2]
                      << " /s, min " << er.front() << ", max " << er.back()
                      << " ; DIF_traction median " << dt_v[n / 2]
                      << " (min " << dt_v.front() << ", max " << dt_v.back()
                      << ")\n";
            std::size_t sat = 0;
            for (double x : er) if (x > 1.0e2) ++sat;
            std::cout << "[FDEM]   " << (100.0 * sat / n)
                      << " % des insertions sont AU PLATEAU de DIF_traction "
                         "(edot > 1e2 /s) : sur cette part le facteur ne "
                         "discrimine plus rien, il agit comme une simple "
                         "multiplication de ft par 1,85.\n";
        }
    }
    if (gcAdaptive_)
        std::cout << "[FDEM] adaptive contact activation: " << nActivated_
                  << " / " << pool_.size() << " exterior faces activated ("
                  << (pool_.empty() ? 0.0
                                    : 100.0 * nActivated_ / pool_.size())
                  << " %)\n";
    if (std::getenv("RKM_GCLOG")) {        // diagnostic : premiere force et
        for (std::size_t k = 0; k < poolTouch_.size(); ++k) {  // activation
            long ac = (k < actStep_.size()) ? actStep_[k] : -1;
            if (poolTouch_[k] >= 0 || ac >= 0)
                std::cout << "[GCLOG] pool " << k << " touch " << poolTouch_[k]
                          << " act " << ac << "\n";
        }
    }
    if (scen_ == Scenario::BRAZILIAN
        || (scen_ == Scenario::TENSION && tensionPlatens_)) {
        // Free quasi-static check: two platens must read the same force, and
        // their imbalance IS the inertia the specimen is carrying. Nobody
        // formalises it as a criterion, every platen-loaded code computes it.
        double a = std::abs(plTop_.F.y()), b = std::abs(plBot_.F.y());
        double mn = 0.5 * (a + b);
        std::cout << "[FDEM] platen balance (|Ftop|-|Fbot|)/mean = "
                  << (mn > 0.0 ? 100.0 * (a - b) / mn : 0.0)
                  << " % (quasi-static if a few %)\n";
    }

    if (voronoi_) {
        // grain/phase bookkeeping: achieved area fractions and the
        // inter/intra-granular split of the broken joints — the observable
        // a GBM exists to produce.
        std::vector<double> aPh(phases_.n(), 0.0);
        double aTot = 0.0;
        for (const auto& e : el_) { aPh[e.phase] += e.A0; aTot += e.A0; }
        std::cout << "[FDEM] grains: " << nGrains_ << ", phases:";
        for (int p = 0; p < phases_.n(); ++p)
            std::cout << " " << phases_.name[p] << " "
                      << 100.0 * aPh[p] / aTot << "%"
                      << " (target " << 100.0 * phases_.fraction[p] << "%)";
        std::cout << "\n";
        long nt[3] = {0, 0, 0}, nb[3] = {0, 0, 0};
        for (const auto& J : jt_) {
            ++nt[J.type];
            if (J.D >= 1.0) ++nb[J.type];
        }
        std::cout << "[FDEM] joints intra/homo/hetero: " << nt[0] << "/"
                  << nt[1] << "/" << nt[2] << ", broken: " << nb[0] << "/"
                  << nb[1] << "/" << nb[2];
        long nbTot = nb[0] + nb[1] + nb[2];
        if (nbTot > 0)
            std::cout << "  (intergranular fraction "
                      << 100.0 * (nb[1] + nb[2]) / (double)nbTot << " %)";
        std::cout << "\n";
    }

    if (confP_ > 0.0) {
        // Falsifiable check of the follower load: away from the loaded faces
        // the specimen must sit at sigma_xx = -p. A large gap means the
        // pressure never diffused (ramp too fast / run too short) or something
        // is holding the surface back (Lysmer springs, grips).
        double err = 100.0 * (confAchieved_ + confP_) / confP_;
        std::cout << "[FDEM] confinement: target " << confP_ / 1e6
                  << " MPa, achieved mean sigma_xx in the core at the end of "
                     "the ramp = " << confAchieved_ / 1e6 << " MPa (" << err
                  << " %); at the END of the run (axial load included) "
                  << achievedConfinement() / 1e6 << " MPa\n";
    }

    if (scen_ == Scenario::BRAZILIAN) {
        // Where did it break? A valid brazilian test fails on the LOAD
        // DIAMETER: the crack runs vertically through the centre. Failure that
        // wanders off the axis (or crushes under a platen) is not an indirect
        // tension measurement, so the diametral fraction is reported next to
        // the strength and must be read WITH it.
        long nbTot = 0, nbAxis = 0;
        double band = cfg_.getd("diametralBand", 0.15) * discR_;
        double sumAbs = 0.0;
        for (const auto& J : jt_) {
            if (J.D < 1.0) continue;
            Eigen::Vector2d mid = 0.5 * (X0_[J.a1] + X0_[J.a2]);
            double dx = std::abs(mid.x() - discC_.x());
            ++nbTot;
            sumAbs += dx / discR_;
            if (dx <= band) ++nbAxis;
        }
        // Parameter-free verification of the setup: on an elastic disc under a
        // diametral load the centre carries sigma_xx = +2P/(pi D t) (exactly
        // what the ISRM formula reports) and sigma_yy = -6P/(pi D t). The
        // RATIO -3 involves no material constant, no geometry and not even the
        // load: it checks the platens, the disc cut and the stress recovery in
        // one number. Finite flat platens (instead of the line load of the
        // closed form) spread the contact and soften both slightly.
        double rXX = gSigT_ != 0.0 ? gXX_ / gSigT_ : 0.0;
        bool passXX = rXX > 0.85 && rXX < 1.25;
        // The ELASTIC-BAND gauge is the one to read: it checks the platens,
        // the disc cut and the stress recovery where the closed-form solution
        // actually applies. The single-step gauge above is kept for continuity
        // but it is read at the last intact step, which on an adaptive run is
        // already past the onset of core softening, so a low ratio there says
        // the CENTRE has yielded, not that the setup is wrong.
        if (eN_ > 0) {
            double mXX = eSumXX_ / eN_, mYY = eSumYY_ / eN_,
                   mSig = eSumSig_ / eN_;
            double eR = mSig != 0.0 ? mXX / mSig : 0.0;
            std::cout << "[FDEM] brazilian ELASTIC-BAND gauge (sigma_t in ["
                      << eGaugeLo_ << ", " << eGaugeHi_ << "] x ft, disc "
                         "intact, " << eN_ << " history steps, "
                      << 100.0 * eDfrac_ << " % of joints damaged at the top of "
                         "the band):\n"
                      << "[FDEM]   mean sigma_t = " << mSig / 1e6
                      << " MPa, mean centre sigma_xx = " << mXX / 1e6
                      << " MPa -> ratio " << eR << "  ["
                      << (eR > 0.85 && eR < 1.25 ? "PASS" : "FAIL")
                      << "]  (band 0.85-1.25; mean sigma_yy = " << mYY / 1e6
                      << " MPa, sigma_yy/sigma_xx = "
                      << (mXX != 0.0 ? mYY / mXX : 0.0) << ")\n";
        } else {
            std::cout << "[FDEM] brazilian ELASTIC-BAND gauge: NOT measured — "
                         "the disc broke before sigma_t reached "
                      << eGaugeLo_ << " x ft, or the run stopped below it\n";
        }
        std::cout << "[FDEM] brazilian elastic gauge (last step before first "
                     "breakage, sigma_t = " << gSigT_ / 1e6 << " MPa):\n"
                  << "[FDEM]   centre sigma_xx = " << gXX_ / 1e6
                  << " MPa, expected +" << gSigT_ / 1e6 << " MPa -> ratio "
                  << rXX << "  [" << (passXX ? "PASS" : "FAIL")
                  << "]  (this IS what sigma_t reports, band 0.85-1.25)\n"
                  << "[FDEM]   centre sigma_yy = " << gYY_ / 1e6
                  << " MPa, sigma_yy/sigma_xx = "
                  << (gXX_ != 0.0 ? gYY_ / gXX_ : 0.0)
                  << " (informative only: -3 is the LINE-load value; a finite "
                     "bearing cuts the centre compression much more than the "
                     "tension, so a flattened disc reads well above it)\n";
        if (brazPlatens_)
            std::cout << "[FDEM]   bearing: " << gNC_ << " nodes in contact, "
                         "participation ratio " << gPR_
                      << " (= effective number carrying the load; 1 means a "
                         "single asperity takes everything)\n";
        std::cout << "[FDEM]   sub-critical damage at peak: "
                  << 100.0 * gDfrac_ << " % of joints above D = 0.01, mean D = "
                  << gDmean_ << " (diffuse ratcheting of the intrinsic penalty "
                     "— if this eats the strength it must track the deficit)\n";
        if (std::getenv("ROCKIM_BRAZ_DEBUG")) {         // load-path profile
            const int NB = 10;
            std::vector<double> sx(NB, 0.0), sy(NB, 0.0), ar(NB, 0.0);
            for (const auto& e : el_) {
                Eigen::Vector2d c = (X0_[e.n[0]] + X0_[e.n[1]] + X0_[e.n[2]]) / 3.0;
                if (std::abs(c.x() - discC_.x()) > 0.1 * discR_) continue;
                int b = (int)((c.y() - (discC_.y() - discR_)) / (2.0 * discR_) * NB);
                b = std::clamp(b, 0, NB - 1);
                sx[b] += e.sxx * e.A0; sy[b] += e.syy * e.A0; ar[b] += e.A0;
            }
            Eigen::Vector2d fNet = Eigen::Vector2d::Zero();
            for (int i = 0; i < (int)X0_.size(); ++i) fNet += f_[i];
            Eigen::Vector2d Ft = brazPlatens_ ? plTop_.F : arcTop_.F;
            Eigen::Vector2d Fb = brazPlatens_ ? plBot_.F : arcBot_.F;
            std::cout << "[DBG] Ftop = (" << Ft.x() << ", " << Ft.y()
                      << ") N, Fbot = (" << Fb.x() << ", " << Fb.y()
                      << ") N, net nodal = (" << fNet.x() << ", " << fNet.y()
                      << ") N\n";
            for (int b = 0; b < NB; ++b)
                if (ar[b] > 0)
                    std::cout << "[DBG] band " << b << " sxx = "
                              << sx[b] / ar[b] / 1e6 << " MPa, syy = "
                              << sy[b] / ar[b] / 1e6 << " MPa\n";
        }
        std::cout << "[FDEM] brazilian (D = " << 2.0 * discR_ << " m, t = "
                  << thk_ << " m, bearing width " << arcTop_.length << " m):\n"
                  << "[FDEM]   peak force P = " << peakF_ << " N at t = "
                  << tPeak_ << " s (" << nBrokenAtPeak_ << " joints already "
                  << "broken; peak "
                  << (peakLocked_ ? "LOCKED at the post-failure load drop"
                                  : "NOT locked — no load drop seen, the run "
                                    "may have stopped before failure") << ")\n"
                  << "[FDEM]   indirect tensile strength sigma_t = 2P/(pi D t) = "
                  << sigmaTpeak_ / 1e6 << " MPa\n"
                  << "[FDEM]   ratio to the bulk ft (" << mat_.ft / 1e6
                  << " MPa) = " << sigmaTpeak_ / mat_.ft << "\n";
        if (nbTot > 0)
            std::cout << "[FDEM]   crack location: " << nbAxis << " / " << nbTot
                      << " broken joints within " << 100.0 * band / discR_
                      << " % of R from the load axis ("
                      << 100.0 * nbAxis / (double)nbTot
                      << " % diametral), mean |x-xc|/R = " << sumAbs / nbTot
                      << "\n";
        else
            std::cout << "[FDEM]   NO joint broke: the disc never failed — "
                         "lengthen T or raise pullV\n";
        std::cout << "[FDEM]   fragments = " << nFrag_ << "\n";
        return;
    }

    if (scen_ == Scenario::TENSION) {
        if (tensionPlatens_) {
            // COMPRESSION through platens (UCS / triaxial): the specimen fails
            // by shear-band formation and the peak is the compressive
            // strength — comparing it against the tensile ft would print a
            // nonsense FAIL. Report the raw figures.
            long nT = 0, nS = 0;
            for (const auto& J : jt_) {
                if (J.bmode == 1) ++nT;
                else if (J.bmode == 2) ++nS;
            }
            double ep = 0.0, es = 0.0, eg = 0.0;
            gaugeStrain(ep, es, eg);
            std::cout << "[FDEM] uniaxial/triaxial compression (platens):\n"
                      << "[FDEM]   peak axial stress = " << sigmaPeak_ / 1e6
                      << " MPa (" << sigmaPeak_ / mat_.ft
                      << " x ft, cohesion " << mat_.cohesion / 1e6 << " MPa)"
                      << (peakLockedU_ ? "  [peak LOCKED at the post-failure "
                                         "load drop]"
                                      : "  [peak NOT locked: no load drop seen "
                                        "— the run may have stopped before "
                                        "failure]") << "\n"
                      << "[FDEM]   broken joints = " << nBroken_ << " / "
                      << jt_.size() << " (mode: " << nT << " tensile, " << nS
                      << " shear";
            if (nT + nS > 0)
                std::cout << ", " << 100.0 * nS / (double)(nT + nS)
                          << " % shear";
            std::cout << ")\n"
                      << "[FDEM]   final axial strain: platen " << ep
                      << ", faces " << es << ", extensometer " << eg << "\n"
                      << "[FDEM]   fragments = " << nFrag_ << "\n";
            if (confP_ > 0.0)
                std::cout << "[FDEM]   confinement: target " << confP_ / 1e6
                          << " MPa, achieved (gauged in the elastic range) "
                          << confAchieved_ / 1e6 << " MPa\n";
            return;
        }
        if (!cfg_.getb("verifyFt", true)) {
            // demo configs with deliberately weakened boundaries: the
            // macroscopic strength is MEANT to sit below the bulk ft, so a
            // PASS/FAIL against ft would be noise. Report the ratio instead
            // (reference: a flat weak-boundary path fails at gbAlphaTen*ft,
            // tortuosity raises it some percent above that).
            std::cout << "[FDEM] tension result (GB-controlled, verifyFt=off):\n"
                      << "[FDEM]   peak macro stress = " << sigmaPeak_ / 1e6
                      << " MPa = " << sigmaPeak_ / mat_.ft << " x bulk ft"
                      << " (weak-boundary reference gbAlphaTen = "
                      << phases_.aTen << ")\n"
                      << "[FDEM]   broken joints = " << nBroken_ << " / "
                      << jt_.size() << "\n";
            return;
        }
        double err = 100.0 * (sigmaPeak_ - mat_.ft) / mat_.ft;
        // On the voronoi mesh the crack must follow a tortuous joint path
        // whose facets are inclined to the loading axis, so the macroscopic
        // peak sits ABOVE ft (each inclined facet sees sigma*cos^2). The
        // verification band is widened on that side and the excess is
        // expected physics, not an error of the joint law.
        double hi = voronoi_ ? 25.0 : 5.0;
        bool pass = err > -5.0 && err < hi;
        std::cout << "[FDEM] tension verification ("
                  << (voronoi_ ? "voronoi grains" : "uniform strip") << "):\n"
                  << "[FDEM]   peak macro stress = " << sigmaPeak_ / 1e6
                  << " MPa, expected ft = " << mat_.ft / 1e6
                  << " MPa, error = " << err << " %  ["
                  << (pass ? "PASS" : "FAIL") << "]\n"
                  << "[FDEM]   broken joints = " << nBroken_ << " / " << jt_.size()
                  << "\n";
        return;
    }
    double Es = detachedVol_ > 0 ? work_ / detachedVol_ : 0.0;
    std::cout << "[FDEM] peak tool force   : " << peakF_ << " N/m\n"
              << "[FDEM] tool work output  : " << work_ << " J/m";
    if (tool_.motion == Tool::Motion::FREE)
        std::cout << "  (tool KE loss: " << toolKE0_ - tool_.ke() << " J/m)";
    std::cout << "\n[FDEM] broken joints     : " << nBroken_ << " / " << jt_.size()
              << "\n[FDEM] fragments         : " << nFrag_
              << " (detached vol " << detachedVol_ << " m^3/m)"
              << "\n[FDEM] specific energy   : " << Es << " J/m^3\n";
    if (muCRes_ >= 0.0) {                                        // WP6
        std::cout << "[FDEM] contact residuel  : " << nCtcPulv_
                  << " evaluations au mu pulverise (" << muCRes_ << ")";
        if (tCtcPulv0_ >= 0.0)
            std::cout << ", premier engagement a t = " << tCtcPulv0_ << " s";
        else
            std::cout << " — JAMAIS engage (aucun element pulverise n a "
                         "touche un contact)";
        std::cout << "\n";
    }
    brushReport();          // no-op si le balai n'a pas ete arme
}

// ---------------------------------------------------------------------------
// selftest-potential2d — le test decisif du chantier A3.
//
// Deux corps RIGIDES triangulaires (3 ddl chacun : centre + rotation), lies
// uniquement par le contact par potentiel, sans frottement. Saute-mouton
// identique au solveur (forces puis v puis x), compteur de travail identique
// (somme des f.v AVANT le kick). Deux phases :
//   1. frontale symetrique  — collision elastique de masses egales : A doit
//      s'arreter, B repartir a v0 ;
//   2. oblique (B decale)   — couple non nul, le clip est asymetrique et
//      tourne : la conservation doit tenir aussi en rotation.
// Verdicts : |W_contact| / KE0 et |dKE| / KE0 au niveau de l'erreur du
// saute-mouton (<< 1), quantite de mouvement conservee MACHINE (la 3e loi
// est exacte par construction dans pairForce). Le contact penalite
// quasi-plastique du solveur dissipe ~80 % d'un rebond PAR CONSTRUCTION :
// c'est l'ecart categorique que ce test verrouille.
// ---------------------------------------------------------------------------
int potentialSelftest(const std::string& csvPath) {
    using V2 = Eigen::Vector2d;
    std::ofstream csv(csvPath);
    csv << "phase,t,xA,vA,xB,vB,area,W,KE\n";

    struct Rigid {
        V2 c, v;                            // centre, vitesse
        double th = 0.0, om = 0.0;          // rotation, vitesse angulaire
        V2 r0[3];                           // sommets dans le repere du corps
        double m = 3.0, I = 0.0;            // 3 masses ponctuelles unite
        void pos(V2 p[3]) const {
            double cs = std::cos(th), sn = std::sin(th);
            for (int k = 0; k < 3; ++k)
                p[k] = c + V2(cs * r0[k].x() - sn * r0[k].y(),
                              sn * r0[k].x() + cs * r0[k].y());
        }
        V2 vel(const V2& x) const {         // vitesse rigide au point x
            V2 r = x - c;
            return v + om * V2(-r.y(), r.x());
        }
        void setTri(const V2& A, const V2& B, const V2& C) {
            c = (A + B + C) / 3.0;
            r0[0] = A - c;
            r0[1] = B - c;
            r0[2] = C - c;
            I = r0[0].squaredNorm() + r0[1].squaredNorm()
              + r0[2].squaredNorm();
        }
    };

    const double p = 1.0e3;                 // penalite du potentiel
    const double dt = 1.0e-4;
    const double v0 = 1.0;
    int fails = 0;
    double worstW = 0.0, worstKE = 0.0, worstP = 0.0;

    for (int phase = 1; phase <= 2; ++phase) {
        Rigid A, B;
        // triangle unite pointe vers +x, et son miroir pointe vers -x
        A.setTri(V2(-1.10, -0.50), V2(-0.10, 0.00), V2(-1.10, 0.50));
        double yoff = (phase == 2) ? 0.22 : 0.0;   // oblique : couple non nul
        B.setTri(V2(1.10, -0.50 + yoff), V2(1.10, 0.50 + yoff),
                 V2(0.10, 0.00 + yoff));
        A.v = V2(v0, 0.0);
        B.v = V2(0.0, 0.0);

        double KE0 = 0.5 * A.m * A.v.squaredNorm();
        V2 P0 = A.m * A.v + B.m * B.v;
        double W = 0.0, areaMax = 0.0;
        long nTouch = 0;

        const long nSteps = (long)(3.0 / dt);
        for (long s = 0; s < nSteps; ++s) {
            V2 pa[3], pb[3];
            A.pos(pa);
            B.pos(pb);
            pot::PairForce R;
            V2 FA = V2::Zero(), FB = V2::Zero();
            double tA = 0.0, tB = 0.0;
            if (pot::pairForce(pa, pb, p, R)) {
                ++nTouch;
                areaMax = std::max(areaMax, R.area);
                for (int k = 0; k < 3; ++k) {
                    // compteur de travail du solveur : f.v AVANT le kick
                    W += (R.fA[k].dot(A.vel(pa[k]))
                          + R.fB[k].dot(B.vel(pb[k]))) * dt;
                    FA += R.fA[k];
                    FB += R.fB[k];
                    tA += pot::cross2(pa[k] - A.c, R.fA[k]);
                    tB += pot::cross2(pb[k] - B.c, R.fB[k]);
                }
            }
            A.v += FA / A.m * dt;
            A.om += tA / A.I * dt;
            B.v += FB / B.m * dt;
            B.om += tB / B.I * dt;
            A.c += A.v * dt;
            A.th += A.om * dt;
            B.c += B.v * dt;
            B.th += B.om * dt;
            if (s % 200 == 0) {
                double KE = 0.5 * A.m * A.v.squaredNorm()
                          + 0.5 * B.m * B.v.squaredNorm()
                          + 0.5 * A.I * A.om * A.om
                          + 0.5 * B.I * B.om * B.om;
                csv << phase << "," << s * dt << "," << A.c.x() << ","
                    << A.v.x() << "," << B.c.x() << "," << B.v.x() << ","
                    << (s == 0 ? 0.0 : areaMax) << "," << W << "," << KE
                    << "\n";
            }
        }

        double KE1 = 0.5 * A.m * A.v.squaredNorm()
                   + 0.5 * B.m * B.v.squaredNorm()
                   + 0.5 * A.I * A.om * A.om + 0.5 * B.I * B.om * B.om;
        V2 P1 = A.m * A.v + B.m * B.v;
        double wRel = std::abs(W) / KE0;
        double keRel = std::abs(KE1 - KE0) / KE0;
        double pRel = (P1 - P0).norm() / P0.norm();
        worstW = std::max(worstW, wRel);
        worstKE = std::max(worstKE, keRel);
        worstP = std::max(worstP, pRel);
        std::cout << "[POT] phase " << phase
                  << (phase == 1 ? " (frontale)" : " (oblique) ")
                  << ": contact " << nTouch << " pas, aire max " << areaMax
                  << "\n[POT]   vA_fin = (" << A.v.x() << ", " << A.v.y()
                  << "), vB_fin = (" << B.v.x() << ", " << B.v.y()
                  << "), omB = " << B.om
                  << "\n[POT]   |W_contact|/KE0 = " << wRel
                  << ", |dKE|/KE0 = " << keRel
                  << ", |dP|/|P0| = " << pRel << "\n";
        if (phase == 1) {
            // collision elastique de masses egales : transfert quasi total
            if (std::abs(B.v.x() - v0) > 0.02 || std::abs(A.v.x()) > 0.02)
                ++fails;
        }
    }

    // Seuils. LA preuve de conservation est |dKE|/KE0 (mesure : 3.7e-12 en
    // frontale, 5.1e-7 en oblique — l'erreur du saute-mouton en rotation) et
    // la quantite de mouvement MACHINE (3e loi exacte de pairForce). Le
    // compteur de travail lit v AVANT le kick (la convention gcWork_ du
    // solveur) : sur une force conservative il porte un biais systematique
    // POSITIF de Sum |F|^2 dt^2 / 2m — un artefact O(dt) de la mesure, pas
    // de la physique (mesure : ~8e-4 de KE0 ici, a comparer aux ~80 % que le
    // contact penalite quasi-plastique dissipe PAR CONSTRUCTION sur le meme
    // rebond). Ce biais existera aussi dans le gcWork_ des runs en mode
    // potential : un petit POSITIF n'y est pas une pathologie, contrairement
    // au mode penalty ou tout positif est une injection.
    std::cout << "pot_work_rel = " << worstW << "\n"
              << "pot_ke_rel = " << worstKE << "\n"
              << "pot_mom_rel = " << worstP << "\n";
    bool ok = fails == 0 && worstW < 5e-3 && worstKE < 1e-5
              && worstP < 1e-12;
    std::cout << (ok ? "[PASS]" : "[FAIL]")
              << " selftest-potential2d : contact conservatif de Munjiza "
                 "(la conservation est jugee sur dKE ; le compteur de "
                 "travail porte un biais O(dt) documente)\n";
    return ok ? 0 : 1;
}

} // namespace rockim
