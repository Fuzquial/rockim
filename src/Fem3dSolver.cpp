// ---------------------------------------------------------------------------
// Fem3dSolver — 3D continuum FEM with pluggable laws. Mesh/boundaries/tool
// mirror the (verified) Fdem3dSolver; the bulk physics is delegated to
// MatLaw and fracture is smeared damage + erosion, as in the 2D fem module.
// ---------------------------------------------------------------------------
#include "rockim/Fem3dSolver.hpp"
#include "rockim/Guards.hpp"
#include "rockim/Fem3dKinematics.hpp"
#include "rockim/ToolSignorini.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>

#include "rockim/RandomField.hpp"
#include "rockim/VtkWriter.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace rockim {

Fem3dSolver::Fem3dSolver(const Config& cfg, std::string outDir)
    : cfg_(cfg), out_(std::move(outDir)) {}

void Fem3dSolver::init() {
    mat_ = Material::from(cfg_);
    PhaseSet::validate(mat_, "global");

    std::string sc = cfg_.gets("scenario", "percussion");
    if      (sc == "percussion") scen_ = Scenario::PERCUSSION;
    else if (sc == "shear")      scen_ = Scenario::SHEAR;
    else if (sc == "tension")    scen_ = Scenario::TENSION;
    else throw std::runtime_error("fem3d scenario must be percussion | "
                                  "shear | tension (tension with pullV < 0 "
                                  "is the uniaxial compression test)");

    W_ = cfg_.getd("W", 0.1);
    D_ = cfg_.getd("D", 0.1);
    H_ = cfg_.getd("H", 0.08);
    nx_ = cfg_.geti("nx", 24);
    ny_ = cfg_.geti("ny", 24);
    nz_ = cfg_.geti("nz", 18);
    T_ = cfg_.getd("T", 2e-4);
    damping_ = cfg_.getd("dampingLocal",
                         scen_ == Scenario::TENSION ? 0.7 : 0.05);
    nanEvery_ = cfg_.geti("nanCheckEvery", 256);       // C4 (w20), 0 = off

    // mesh = grid (defaut, tets de Kuhn miroites + jitter : STRUCTURE) ou
    // mesh = file (Gmsh MSH 2.2 non structure, la regle de la these pour tout
    // essai dont la reponse est un trajet de fissure ou une energie de bande)
    {
        std::string mesh = cfg_.gets("mesh", "grid");
        if (mesh == "file") buildMeshFile();
        else buildMesh();
    }
    law_ = MatLaw::make(cfg_.gets("law", "dpr"), mat_, cfg_, lcMax_);
    fixedDam_ = cfg_.gets("tensionDamage", "scalar") == "fixed";   // validee par make

    // initial element centroids: dpdfh seeds its deterministic Weibull
    // draws from a spatial hash of these coordinates (the VUMAT's coordMp)
    for (auto& e : el_)
        e.st.x0 = 0.25 * (X0_[e.n[0]] + X0_[e.n[1]]
                          + X0_[e.n[2]] + X0_[e.n[3]]);

    // Per-element Weibull strength heterogeneity (matWeibullM > 1): a
    // mean-1 factor per element scales the local strengths (wired to
    // saksala2011's SDV15/16 mechanism — the sandbox version of the
    // paper's FIELD random strengths, which are what turns the diffuse
    // damage patch into DISCRETE crack bands: with homogeneous strength
    // nothing selects individual cracks). strengthCorrLength = 0 draws
    // independently per element (the paper's choice); > 0 samples ONE
    // Gaussian random field of that correlation length through the
    // Gaussian copula (RandomField3, same keys and semantics as the 2D
    // FDEM joint statistics): the field lives in space, independent of the
    // mesh, so two meshes see the same weak zones.
    double mW = cfg_.getd("matWeibullM", 0.0);
    if (mW > 0.0) {
        if (mW <= 1.0)
            throw std::runtime_error("matWeibullM must be > 1 (the paper "
                                     "uses 3)");
        unsigned fseed = (unsigned)cfg_.geti("fieldSeed",
                                             cfg_.geti("seed", 12345) + 777);
        double gam = std::tgamma(1.0 + 1.0 / mW);
        auto weib = [&](double u) {
            u = std::clamp(u, 1e-12, 1.0 - 1e-12);
            return std::pow(-std::log(1.0 - u), 1.0 / mW) / gam;
        };
        double ell = cfg_.getd("strengthCorrLength", 0.0);
        double mn = 1e300, mx = 0.0, sum = 0.0;
        if (ell > 0.0) {
            double ellB = cfg_.getd("strengthCorrLengthB", ell);
            double ang = cfg_.getd("strengthCorrAngleDeg", 0.0);
            RandomField3 F(W_, D_, H_, ell, ellB, ang, fseed);
            for (auto& e : el_) {
                Eigen::Vector3d cen = 0.25 * (X0_[e.n[0]] + X0_[e.n[1]]
                                              + X0_[e.n[2]] + X0_[e.n[3]]);
                double g = F(cen);
                e.st.ftScale = weib(0.5 * std::erfc(-g / std::sqrt(2.0)));
            }
        } else {
            std::mt19937 rng(fseed);
            std::uniform_real_distribution<double> U(0.0, 1.0);
            for (auto& e : el_) e.st.ftScale = weib(U(rng));
        }
        for (const auto& e : el_) {
            mn = std::min(mn, e.st.ftScale);
            mx = std::max(mx, e.st.ftScale);
            sum += e.st.ftScale;
        }
        std::cout << "[FEM3D] Weibull strengths: m = " << mW
                  << (ell > 0.0 ? " correlated, ell = " + std::to_string(ell)
                                : std::string(" independent per element"))
                  << ", factor mean/min/max = " << sum / el_.size() << "/"
                  << mn << "/" << mx << "\n";
    }

    kp_ = mat_.E * hmin_;                              // tool penalty [N/m]
    // toolContact = penalty | signorini (port du 2026-09-04 par la session
    // parallele, ToolSignorini.hpp ; defaut penalty = bit-identique)
    {
        std::string tc = cfg_.gets("toolContact", "penalty");
        if (tc != "penalty" && tc != "signorini")
            throw std::runtime_error("toolContact must be penalty | signorini");
        toolSig_ = (tc == "signorini");
        toolSigRelax_ = cfg_.getd("toolSignoriniRelax", 0.0);
        if (toolSigRelax_ < 0.0 || toolSigRelax_ > 1.0)
            throw std::runtime_error("toolSignoriniRelax must be in [0, 1]");
        if (toolSig_)
            std::cout << "[FEM3D] contact outil : SIGNORINI en vitesse (impulsion, "
                         "relax " << toolSigRelax_ << ") ; la penalite kp n'est "
                         "pas utilisee\n";
    }
    // contactPenaltyFactor (2026-09-04, percussion quart de bloc) : sur un
    // Delaunay hmin est un SLIVER (0,11 mm pour des elements de 0,5 mm), donc
    // kp = E hmin = 8,6e6 N/m par noeud, 4 a 5 fois plus souple que la raideur
    // nodale de la roche (interpenetration 0,07 mm en moyenne, 0,43 mm au
    // pole = l'indentation). Facteur multiplicatif opt-in ; a 1 (defaut)
    // aucune operation n'est appliquee (bit-identique). La stabilite est
    // preservee par construction : computeStableDt borne dt par
    // 2 sqrt(m_min/kp) (pulsation de penalite omega_p = sqrt(kp/m_min)), le
    // rapport omega_p/(2/dt) est imprime a l'init et un avertissement est
    // emis si le rapport COMBINE sqrt((omega_el dt)^2 + (omega_p dt)^2)/2
    // depasse 0,5 (revue : un noeud de contact porte aussi la raideur
    // d'element, omega_el = 2 c_P/hmin ; le ressort seul vaut au plus dtFactor).
    // REVUE : kpFactor (la cle du mode fem 2D, FemSolver.cpp) est acceptee
    // comme ALIAS en fem3d — auparavant ignoree en silence ; les deux cles
    // avec des valeurs differentes sont refusees. Aucun deck fem3d du depot ne
    // porte kpFactor : bit-identique.
    // toolPulseForce / toolPulseTable : impulsion de force imposee sur l'outil
    // (voir Fem3dSolver.hpp). Table "t:a t:a ..." en secondes / amplitude.
    pulseF_ = cfg_.getd("toolPulseForce", 0.0);
    if (pulseF_ < 0.0)
        throw std::runtime_error("toolPulseForce must be >= 0 (N, applied along -z)");
    if (pulseF_ > 0.0) {
        std::istringstream ss(cfg_.gets("toolPulseTable", ""));
        std::string tok;
        while (ss >> tok) {
            auto c = tok.find(':');
            if (c == std::string::npos)
                throw std::runtime_error("toolPulseTable: expected t:a pairs, got '" + tok + "'");
            double tt = std::stod(tok.substr(0, c)), aa = std::stod(tok.substr(c + 1));
            if (!pulseT_.empty() && tt <= pulseT_.back())
                throw std::runtime_error("toolPulseTable: times must increase");
            pulseT_.push_back(tt);
            pulseA_.push_back(aa);
        }
        if (pulseT_.size() < 2)
            throw std::runtime_error("toolPulseForce > 0 requires toolPulseTable with >= 2 pairs t:a");
        pulseZ_ = cfg_.getd("toolPulseImpedance", 0.0);
        if (pulseZ_ < 0.0)
            throw std::runtime_error("toolPulseImpedance must be >= 0 (N s/m, rod rho c A)");
        if (pulseZ_ > 0.0)
            std::cout << "[FEM3D]   source a impedance : F = 2 F_inc - Z v, Z = " << pulseZ_
                      << " N s/m (F_inc = toolPulseForce x a(t))\n";
        std::cout << "[FEM3D] impulsion de force sur l'outil : " << pulseF_
                  << " N x table de " << pulseT_.size() << " points ("
                  << pulseT_.front() << " -> " << pulseT_.back() << " s)\n";
    } else if (cfg_.has("toolPulseTable")) {
        throw std::runtime_error("toolPulseTable given without toolPulseForce > 0 "
                                 "(nothing is applied in silence)");
    }
    kpFac_ = cfg_.getd("contactPenaltyFactor", 1.0);
    if (cfg_.has("kpFactor")) {
        const double kf = cfg_.getd("kpFactor", 1.0);
        if (cfg_.has("contactPenaltyFactor") && kf != kpFac_)
            throw std::runtime_error("kpFactor and contactPenaltyFactor both given "
                                     "with different values (fem3d : kpFactor is an "
                                     "alias of contactPenaltyFactor = factor x E hmin)");
        kpFac_ = kf;
    }
    if (!(kpFac_ > 0.0) || !std::isfinite(kpFac_))
        throw std::runtime_error("contactPenaltyFactor > 0 expected (1 = E hmin, "
                                 "the historical tool penalty)");
    if (kpFac_ != 1.0) kp_ *= kpFac_;
    vtkCap_ = cfg_.getb("vtkCap", false);              // champs VTU capPc / epsVpl (revue)
    muC_ = cfg_.getd("contactMu", 0.5);
    xiC_ = cfg_.getd("contactXi", 0.05);
    vReg_ = cfg_.getd("contactVreg", 1e-3);

    // ---- etude briques (2026-09-03) : cles opt-in, defaut bit-identique ----
    toolDelay_ = cfg_.getd("toolDelay", 0.0);          // outil gele avant t
    symY_ = cfg_.getb("symmetryY", false);             // tranche plane v_y = 0
    symQ_ = cfg_.getb("quarterModel", false);          // quart de bloc : plans x = 0 et y = 0 symetriques
    bottomFree_ = cfg_.getb("bottomFree", false);      // fond libre (absorbing = none | sides) : bloc flottant
    activeNodes_ = cfg_.getb("activeNodes", false);    // masque des orphelins
    erodeDetMin_ = cfg_.getd("erodeDetMin", 0.0);      // soupape det F (0 = off)
    erodeStrainMax_ = cfg_.getd("erodeStrainMax", 0.0); // plafond |eps| (0 = off)
    stats_ = cfg_.getb("fieldStats", false);           // colonnes de champ
    if (toolDelay_ < 0.0 || erodeDetMin_ < 0.0 || erodeDetMin_ >= 1.0
        || erodeStrainMax_ < 0.0)
        throw std::runtime_error("toolDelay >= 0, 0 <= erodeDetMin < 1 and "
                                 "erodeStrainMax >= 0 (0 = off) expected");
    if (activeNodes_) active_.assign(X0_.size(), 1);

    // ---- bulkViscosity = b1 b2 (2026-09-05) : viscosite de volume a la
    // Abaqus/Explicit (Analysis User's Guide, « Bulk viscosity » : defauts
    // Abaqus 0,06 et 1,2). Cle absente OU "0 0" : rien n'est ajoute, aucune
    // colonne, dt inchange (bit-identique). Voir Fem3dSolver.hpp et doc §5.20.
    if (cfg_.has("bulkViscosity")) {
        std::istringstream ss(cfg_.gets("bulkViscosity", ""));
        std::string s1, s2, s3;
        if (!(ss >> s1 >> s2) || (ss >> s3))
            throw std::runtime_error("bulkViscosity expects exactly two numbers "
                                     "'b1 b2' (Abaqus defaults : 0.06 1.2)");
        auto num = [](const std::string& s, const char* what) {
            std::size_t used = 0;
            double d = 0.0;
            try { d = std::stod(s, &used); } catch (const std::exception&) { used = 0; }
            if (used != s.size() || !std::isfinite(d) || d < 0.0)
                throw std::runtime_error(std::string("bulkViscosity: ") + what
                                         + " must be a finite number >= 0, got '" + s + "'");
            return d;
        };
        bvB1_ = num(s1, "b1 (linear)");
        bvB2_ = num(s2, "b2 (quadratic)");
        bv_ = bvB1_ > 0.0 || bvB2_ > 0.0;
        if (!bv_)
            std::cout << "[FEM3D] bulkViscosity = 0 0 : viscosite de volume INACTIVE "
                         "(chemin bit-identique, pas de colonne wBulk)\n";
    }

    // ---- kinematics = biot | hencky (2026-09-05, w16) : mesure de deformation
    // et conjugue. Cle absente ou "biot" : chemin historique (Biot, forces
    // sur la reference, bit-identique). "hencky" : eps = ln U (spectrale de
    // C = F^T F), sigma_c = Cauchy co-rotationnelle, P = J R sigma_c U^-1.
    // Voir Fem3dKinematics.hpp et doc §5.20 ¶ « Cinematique hencky ».
    {
        const std::string kin = cfg_.gets("kinematics", "biot");
        if (kin != "biot" && kin != "hencky")
            throw std::runtime_error("kinematics must be biot | hencky "
                                     "(fem3d ; biot = default, unchanged path)");
        hencky_ = kin == "hencky";
        if (hencky_)
            std::cout << "[FEM3D] cinematique HENCKY : eps = ln U (spectrale de "
                         "F^T F), contrainte de la loi = Cauchy co-rotationnelle, "
                         "forces internes sur la configuration courante "
                         "(P = J R sigma U^-1) ; masse, Lysmer, contact, lc, "
                         "CFL sur la geometrie initiale inchanges\n";
    }

    placeTool();
    setupBoundaries();
    setupConfinement();
    computeStableDt();
    setupProbes();                         // probes = ... (opt-in ; rien sans la cle)
    if (bv_) {
        const double xi = bvB1_;
        std::cout << "[FEM3D] viscosite de volume (bulkViscosity) : b1 = " << bvB1_
                  << " (lineaire, p1 = b1 rho c_d L_e edot), b2 = " << bvB2_
                  << " (quadratique, p2 = rho (b2 L_e edot)^2, compression seule) ; "
                     "c_d = " << mat_.cP() << " m/s, L_e = lc = V0^(1/3) (max "
                  << lcMax_ << " m) ; dt CFL x (sqrt(1 + b1^2) - b1) = "
                  << std::sqrt(1.0 + xi * xi) - xi << " (part lineaire seule, "
                     "la part quadratique -b2^2 L_e edot/c_d n'est pas appliquee)\n";
    }
    // ---- stabilite du contact de penalite (imprimee a chaque init) --------
    // Schema centre : stable si omega dt < 2 (amortissement xi compris :
    // < 2 (sqrt(1 + xi^2) - xi)). Le rapport omega_p/(2/dt) = omega_p dt/2
    // vaut au plus dtFactor par construction de computeStableDt. REVUE : un
    // noeud de contact porte AUSSI la raideur d'element (omega_el = 2 c_P/hmin,
    // omega_el dt/2 = dt/CFL) ; la borne de Rayleigh sur la somme donne
    // omega_tot <= sqrt(omega_el^2 + omega_p^2), et c'est le rapport COMBINE
    // sqrt((dt/CFL)^2 + (omega_p dt/2)^2) qui declenche l'avertissement (pas
    // d'exception : dtFactor = 0,7 est un choix documente sur les maillages
    // Gmsh ; a dtFactor <= 0,35 le combine reste < 0,5 quel que soit kp). En
    // mode signorini la penalite kp n'est pas utilisee : rien n'est imprime.
    if (!toolSig_) {
        double mMin = 1e300;
        for (std::size_t i = 0; i < X0_.size(); ++i)
            if (m_[i] > 0.0) mMin = std::min(mMin, m_[i]);
        const double omegaP = std::sqrt(kp_ / mMin);
        const double ratio = omegaP * dt_ / 2.0;
        const double cflDt = hmin_ / mat_.cP();
        const double ratioEl = dt_ / cflDt;                       // omega_el dt / 2
        const double ratioTot = std::sqrt(ratio * ratio + ratioEl * ratioEl);
        std::cout << "[FEM3D] tool penalty kp = " << kp_ << " N/m (contactPenaltyFactor "
                  << kpFac_ << " x E hmin, hmin = " << hmin_ << " m), m_min = " << mMin
                  << " kg, omega_p = sqrt(kp/m_min) = " << omegaP << " rad/s, omega_p/(2/dt) = "
                  << ratio << " (dt = " << dt_ << " s, CFL " << cflDt
                  << " s, 2 sqrt(m_min/kp) = " << 2.0 / omegaP << " s) ; omega_el/(2/dt) = dt/CFL = "
                  << ratioEl << ", rapport combine sqrt(omega_el^2 + omega_p^2)/(2/dt) = "
                  << ratioTot << "\n";
        if (ratioTot > 0.5)
            std::cout << "[FEM3D] WARNING: combined contact ratio sqrt(omega_el^2 + omega_p^2)/(2/dt) = "
                      << ratioTot << " > 0.5 (spring alone " << ratio << ", element alone "
                      << ratioEl << ") — the contact node is close to the central-difference "
                         "stability limit (1, or " << std::sqrt(1.0 + xiC_ * xiC_) - xiC_
                      << " with contactXi = " << xiC_ << ") : lower dtFactor or "
                         "contactPenaltyFactor\n";
    }

    if (scen_ == Scenario::TENSION) pullV_ = cfg_.getd("pullV", 0.05);
    pullRamp_ = cfg_.getd("pullRamp", 0.0);
    gripFree_ = cfg_.getb("gripLateralFree", false);
    // ---- essai triaxial continu (2026-09-04), opt-in, defaut bit-identique --
    pullDelay_ = cfg_.getd("pullDelay", 0.0);
    gripSection_ = cfg_.getd("gripSection", W_ * D_);
    triax_ = cfg_.getb("triaxStats", false);
    if (pullDelay_ < 0.0 || !(gripSection_ > 0.0))
        throw std::runtime_error("pullDelay >= 0 and gripSection > 0 [m^2] "
                                 "expected");
    if (pullDelay_ > 0.0 && scen_ != Scenario::TENSION)
        throw std::runtime_error("pullDelay only applies to scenario = tension");
    if (pullDelay_ > 0.0)
        std::cout << "[FEM3D] pullDelay " << pullDelay_ << " s : mors superieur "
                     "LIBRE (pression de dessus " << topP_ / 1e6 << " MPa) puis "
                     "prescrit depuis sa position courante, rampe " << pullRamp_
                  << " s ; section de mors " << gripSection_ << " m^2\n";
    toolKE0_ = tool_.ke();
    toolX0_ = tool_.x;                     // suivi lateral (resume, console seulement)

    std::cout << "[FEM3D] law = " << law_->name() << ", " << el_.size()
              << " tets, " << X0_.size() << " nodes, dt = " << dt_
              << " s, steps = " << (long)std::ceil(T_ / dt_) << "\n";
    if (law_->name() != "elastic")
        std::cout << "[FEM3D] DP uniaxial compressive strength (analytic) = "
                  << law_->sigmaCdp() / 1e6 << " MPa\n";
}

// ---------------------------------------------------------------------------
// Kuhn tet mesh with SHARED nodes: 6 tets per hex cell, compatible face
// diagonals, optional interior jitter. Only the exterior faces are kept
// (quiet boundaries); interior faces need no bookkeeping — the continuum is
// glued by the shared nodes themselves.
// ---------------------------------------------------------------------------
void Fem3dSolver::buildMesh() {
    double dx = W_ / nx_, dy = D_ / ny_, dz = H_ / nz_;
    hmin_ = std::min({dx, dy, dz});
    double jit = cfg_.getd("meshJitter", 0.0) * 0.5 * hmin_;
    std::mt19937 rng(cfg_.geti("seed", 12345));
    std::uniform_real_distribution<double> U(-jit, jit);

    int vnx = nx_ + 1, vny = ny_ + 1, vnz = nz_ + 1;
    X0_.resize((std::size_t)vnx * vny * vnz);
    auto vid = [&](int i, int j, int k) { return (k * vny + j) * vnx + i; };
    for (int k = 0; k < vnz; ++k)
        for (int j = 0; j < vny; ++j)
            for (int i = 0; i < vnx; ++i) {
                Eigen::Vector3d p(i * dx, j * dy, k * dz);
                if (jit > 0 && i > 0 && i < nx_ && j > 0 && j < ny_
                    && k > 0 && k < nz_)
                    p += Eigen::Vector3d(U(rng), U(rng), U(rng));
                X0_[vid(i, j, k)] = p;
            }

    std::map<std::array<int, 3>, std::pair<int, std::array<int, 3>>> faces;
    auto addTet = [&](std::array<int, 4> nn) {
        Elem e;
        e.n = nn;
        Eigen::Matrix3d J;
        J.col(0) = X0_[e.n[1]] - X0_[e.n[0]];
        J.col(1) = X0_[e.n[2]] - X0_[e.n[0]];
        J.col(2) = X0_[e.n[3]] - X0_[e.n[0]];
        double det = J.determinant();
        if (det < 0) {
            std::swap(e.n[2], e.n[3]);
            J.col(1) = X0_[e.n[2]] - X0_[e.n[0]];
            J.col(2) = X0_[e.n[3]] - X0_[e.n[0]];
            det = -det;
        }
        e.V0 = det / 6.0;
        if (e.V0 <= 0)                         // C3 (w20) : erreur NOMMEE
            guards::degenerateError("fem3d mesh = grid", (long)el_.size(),
                                    e.n, X0_, guards::kNoIds, e.V0);
        e.lc = std::cbrt(e.V0);
        lcMax_ = std::max(lcMax_, e.lc);
        Eigen::Matrix3d Jinv = J.inverse();
        e.dN.col(1) = Jinv.row(0);
        e.dN.col(2) = Jinv.row(1);
        e.dN.col(3) = Jinv.row(2);
        e.dN.col(0) = -(e.dN.col(1) + e.dN.col(2) + e.dN.col(3));
        int id = (int)el_.size();
        el_.push_back(e);

        const int faceIdx[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
        for (const auto& fi : faceIdx) {
            std::array<int, 3> fn = {e.n[fi[0]], e.n[fi[1]], e.n[fi[2]]};
            std::array<int, 3> key = fn;
            std::sort(key.begin(), key.end());
            auto it = faces.find(key);
            if (it == faces.end()) faces[key] = {id, fn};
            else it->second.first = -1;                // interior: two owners
        }
    };

    // geometry = cylinder carves the structured grid: cells whose centroid
    // lies outside the radius min(W, D)/2 around the vertical axis are
    // dropped. The carved faces become exterior automatically through the
    // face registry; the staircase boundary is documented (keep it far
    // from the impact).
    bool cyl = cfg_.gets("geometry", "box") == "cylinder";
    double Rcyl = 0.5 * std::min(W_, D_);
    double cxc = 0.5 * W_, cyc = 0.5 * D_;

    // meshMirror (default true): MIRRORED Kuhn split — each cell is
    // reflected according to the parity of its (i, j, k) indices, which is
    // face-compatible by construction and alternates the diagonal
    // directions in a checkerboard. The plain split threads ONE global
    // family of diagonals through the mesh and crack patterns visibly
    // snap to it (measured on the cylinder percussion); the mirrored
    // split spreads them over 8 alternating orientations. Set
    // meshMirror = false to recover the previous mesh exactly.
    bool mirror = cfg_.getb("meshMirror", true);

    lcMax_ = 0.0;
    const int perms[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2},
                             {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (int k = 0; k < nz_; ++k)
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) {
                if (cyl) {
                    double xc = (i + 0.5) * dx - cxc;
                    double yc = (j + 0.5) * dy - cyc;
                    if (xc * xc + yc * yc > Rcyl * Rcyl) continue;
                }
                int fx = mirror ? (i & 1) : 0;
                int fy = mirror ? (j & 1) : 0;
                int fz = mirror ? (k & 1) : 0;
                int c[2][2][2];
                for (int a = 0; a < 2; ++a)
                    for (int b = 0; b < 2; ++b)
                        for (int cc = 0; cc < 2; ++cc)
                            c[a][b][cc] = vid(i + (fx ? 1 - a : a),
                                              j + (fy ? 1 - b : b),
                                              k + (fz ? 1 - cc : cc));
                for (const auto& p : perms) {
                    int s[3] = {0, 0, 0};
                    std::array<int, 4> nn;
                    nn[0] = c[0][0][0];
                    s[p[0]] = 1; nn[1] = c[s[0]][s[1]][s[2]];
                    s[p[1]] = 1; nn[2] = c[s[0]][s[1]][s[2]];
                    nn[3] = c[1][1][1];
                    addTet(nn);
                }
            }

    for (auto& [key, fo] : faces) {
        if (fo.first < 0) continue;                    // interior
        // outward orientation vs owning tet centroid
        const Elem& e = el_[fo.first];
        std::array<int, 3>& fn = fo.second;
        Eigen::Vector3d A = X0_[fn[0]], B = X0_[fn[1]], C = X0_[fn[2]];
        Eigen::Vector3d cen = 0.25 * (X0_[e.n[0]] + X0_[e.n[1]]
                                      + X0_[e.n[2]] + X0_[e.n[3]]);
        if ((B - A).cross(C - A).dot((A + B + C) / 3.0 - cen) < 0)
            std::swap(fn[1], fn[2]);
        exterior_.push_back({fo.first, fn});
    }

    u_.assign(X0_.size(), Eigen::Vector3d::Zero());
    v_.assign(X0_.size(), Eigen::Vector3d::Zero());
    f_.assign(X0_.size(), Eigen::Vector3d::Zero());
    m_.assign(X0_.size(), 0.0);
    flag_.assign(X0_.size(), FREE);
    for (const auto& e : el_)
        for (int a = 0; a < 4; ++a)
            m_[e.n[a]] += mat_.rho * e.V0 / 4.0;
    // carved geometries leave unused grid nodes: pin them (zero mass would
    // otherwise divide the integrator)
    {   // C3 (w20) : les seuls noeuds sans element toleres sont ceux de la
        // grille hors du cylindre taille (geometry = cylinder) — comptes et
        // imprimes ; sinon noeud orphelin, sliver ou masse nulle = erreur
        std::vector<double> vols;
        std::vector<std::array<int, 4>> conn;
        vols.reserve(el_.size());
        conn.reserve(el_.size());
        for (const auto& e : el_) { vols.push_back(e.V0); conn.push_back(e.n); }
        std::vector<char> used = guards::referencedMask(X0_.size(), conn);
        std::size_t nUnused = 0;
        for (char c : used) if (!c) ++nUnused;
        if (nUnused > 0 && !cyl) guards::checkOrphans(X0_, conn, guards::kNoIds);
        guards::checkDegenerate("fem3d mesh = grid", vols, conn, X0_, guards::kNoIds);
        guards::checkMasses("fem3d mesh = grid", m_, X0_, used, guards::kNoIds);
        if (nUnused > 0)
            std::cout << "[FEM3D] geometry = cylinder : " << nUnused
                      << " noeuds de grille hors du cylindre (sans element, "
                         "masse nulle) epingles FIXED\n";
    }
    for (std::size_t i = 0; i < X0_.size(); ++i)
        if (m_[i] <= 0.0) flag_[i] = FIXED;

    if (scen_ == Scenario::TENSION) {
        for (int i = 0; i < (int)X0_.size(); ++i) {
            if (X0_[i].z() < 1e-9)           flag_[i] = FIXED;
            else if (X0_[i].z() > H_ - 1e-9) flag_[i] = PRESCRIBED;
        }
        for (int e = 0; e < (int)el_.size(); ++e) {
            double zc = 0.0;
            for (int a = 0; a < 4; ++a) zc += X0_[el_[e].n[a]].z();
            zc /= 4.0;
            if (zc > H_ / 3.0 && zc < 2.0 * H_ / 3.0) midEl_.push_back(e);
        }
    }
}

// ---------------------------------------------------------------------------
// mesh = file (2026-09-04) : maillage tetraedrique NON STRUCTURE importe
// (Gmsh MSH 2.2 ASCII, $Nodes + elements de type 4), meme lecteur que
// Fdem3dSolver::buildMeshFile. Le maillage est translate a l'origine et
// W/D/H sont redefinis depuis la boite englobante (frontieres, outil et
// jauges s'appuient sur les plans x/y/z = 0 et W/D/H). hmin = plus petit
// diametre inscrit 6V/A (pilote la CFL et la penalite de contact) ; lc par
// element = V0^(1/3) comme sur la grille. Chemin grille inchange.
// ---------------------------------------------------------------------------
void Fem3dSolver::buildMeshFile() {
    std::string path = cfg_.reqs("meshFile");
    std::ifstream in(path);
    if (!in) throw std::runtime_error("meshFile: cannot open '" + path + "'");
    std::string line;
    std::map<long, int> id2idx;
    std::vector<Eigen::Vector3d> vpos;
    std::vector<std::array<int, 4>> tets;
    bool sawFormat = false;
    while (std::getline(in, line)) {
        if (line.rfind("$MeshFormat", 0) == 0) {
            std::getline(in, line);
            double ver = std::atof(line.c_str());
            if (ver < 2.0 || ver >= 3.0)
                throw std::runtime_error("meshFile: MSH version "
                    + std::to_string(ver) + " unsupported — export ASCII 2.2");
            sawFormat = true;
        } else if (line.rfind("$Nodes", 0) == 0) {
            long n = 0; in >> n;
            for (long k = 0; k < n; ++k) {
                long id; double x, y, z;
                in >> id >> x >> y >> z;
                id2idx[id] = (int)vpos.size();
                gmshNodeId_.push_back(id);
                vpos.push_back({x, y, z});
            }
        } else if (line.rfind("$Elements", 0) == 0) {
            long n = 0; in >> n;
            for (long k = 0; k < n; ++k) {
                long id; int type, ntags;
                in >> id >> type >> ntags;
                long tag;
                for (int t = 0; t < ntags; ++t) in >> tag;
                int nn = type == 15 ? 1 : type == 1 ? 2 : type == 2 ? 3
                       : type == 4 ? 4 : -1;
                if (nn < 0)
                    throw std::runtime_error("meshFile: element type "
                        + std::to_string(type) + " unsupported (tets only)");
                std::array<int, 4> vv{};
                for (int q = 0; q < nn; ++q) {
                    long nid; in >> nid;
                    if (nn == 4) {
                        auto it = id2idx.find(nid);
                        if (it == id2idx.end())
                            throw std::runtime_error("meshFile: element "
                                + std::to_string(id) + " reference le noeud "
                                "inconnu id " + std::to_string(nid));
                        vv[q] = it->second;
                    }
                }
                if (nn == 4) tets.push_back(vv);
            }
        }
    }
    if (!sawFormat || vpos.empty() || tets.empty())
        throw std::runtime_error("meshFile: no tetrahedra in '" + path + "'");
    // C3 (w20) : noeud jamais reference (point du champ de taille Gmsh...) =
    // erreur nommee, AVANT la translation (coordonnees du fichier)
    guards::checkOrphans(vpos, tets, gmshNodeId_);
    Eigen::Vector3d lo = vpos[0], hi = vpos[0];
    for (const auto& p : vpos) { lo = lo.cwiseMin(p); hi = hi.cwiseMax(p); }
    for (auto& p : vpos) p -= lo;
    W_ = hi.x() - lo.x(); D_ = hi.y() - lo.y(); H_ = hi.z() - lo.z();
    if (!(W_ > 0 && D_ > 0 && H_ > 0))
        throw std::runtime_error("meshFile: degenerate bounding box");
    X0_ = vpos;
    finishMesh(tets);
    std::cout << "[FEM3D] mesh = file: '" << path << "' — " << X0_.size()
              << " nodes, " << el_.size() << " tets, box " << W_ << " x " << D_
              << " x " << H_ << " m, h_min (diametre inscrit) = " << hmin_
              << " m, lc max = " << lcMax_ << " m\n";
}

// elements, registre de faces, masses et drapeaux a partir d'une liste de
// tets sur X0_ (utilise par mesh = file ; la grille garde son propre chemin)
void Fem3dSolver::finishMesh(const std::vector<std::array<int, 4>>& tetsIn) {
    std::map<std::array<int, 3>, std::pair<int, std::array<int, 3>>> faces;
    lcMax_ = 0.0;
    hmin_ = 1e30;
    el_.clear();
    el_.reserve(tetsIn.size());
    for (auto nn : tetsIn) {
        Elem e;
        e.n = nn;
        Eigen::Matrix3d J;
        J.col(0) = X0_[e.n[1]] - X0_[e.n[0]];
        J.col(1) = X0_[e.n[2]] - X0_[e.n[0]];
        J.col(2) = X0_[e.n[3]] - X0_[e.n[0]];
        double det = J.determinant();
        if (det < 0) {
            std::swap(e.n[2], e.n[3]);
            J.col(1) = X0_[e.n[2]] - X0_[e.n[0]];
            J.col(2) = X0_[e.n[3]] - X0_[e.n[0]];
            det = -det;
        }
        e.V0 = det / 6.0;
        if (e.V0 <= 0)                         // C3 (w20) : erreur NOMMEE
            guards::degenerateError("fem3d mesh = file", (long)el_.size(),
                                    e.n, X0_, gmshNodeId_, e.V0);
        e.lc = std::cbrt(e.V0);
        lcMax_ = std::max(lcMax_, e.lc);
        // diametre inscrit 6V/A
        double A = 0.0;
        const int faceIdx[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
        for (const auto& fi : faceIdx) {
            Eigen::Vector3d a = X0_[e.n[fi[0]]], b = X0_[e.n[fi[1]]],
                            c = X0_[e.n[fi[2]]];
            A += 0.5 * (b - a).cross(c - a).norm();
        }
        hmin_ = std::min(hmin_, 6.0 * e.V0 / A);
        Eigen::Matrix3d Jinv = J.inverse();
        e.dN.col(1) = Jinv.row(0);
        e.dN.col(2) = Jinv.row(1);
        e.dN.col(3) = Jinv.row(2);
        e.dN.col(0) = -(e.dN.col(1) + e.dN.col(2) + e.dN.col(3));
        int id = (int)el_.size();
        el_.push_back(e);
        for (const auto& fi : faceIdx) {
            std::array<int, 3> fn = {e.n[fi[0]], e.n[fi[1]], e.n[fi[2]]};
            std::array<int, 3> key = fn;
            std::sort(key.begin(), key.end());
            auto it = faces.find(key);
            if (it == faces.end()) faces[key] = {id, fn};
            else it->second.first = -1;
        }
    }
    exterior_.clear();
    for (auto& [key, fo] : faces) {
        if (fo.first < 0) continue;
        const Elem& e = el_[fo.first];
        std::array<int, 3>& fn = fo.second;
        Eigen::Vector3d A = X0_[fn[0]], B = X0_[fn[1]], C = X0_[fn[2]];
        Eigen::Vector3d cen = 0.25 * (X0_[e.n[0]] + X0_[e.n[1]]
                                      + X0_[e.n[2]] + X0_[e.n[3]]);
        if ((B - A).cross(C - A).dot((A + B + C) / 3.0 - cen) < 0)
            std::swap(fn[1], fn[2]);
        exterior_.push_back({fo.first, fn});
    }
    u_.assign(X0_.size(), Eigen::Vector3d::Zero());
    v_.assign(X0_.size(), Eigen::Vector3d::Zero());
    f_.assign(X0_.size(), Eigen::Vector3d::Zero());
    m_.assign(X0_.size(), 0.0);
    flag_.assign(X0_.size(), FREE);
    for (const auto& e : el_)
        for (int a = 0; a < 4; ++a) m_[e.n[a]] += mat_.rho * e.V0 / 4.0;
    {   // C3 (w20) : sliver (< 1e-6 x mediane) et masse nulle = erreurs
        // nommees ; plus jamais d'epinglage FIXED silencieux (broche fantome)
        std::vector<double> vols;
        std::vector<std::array<int, 4>> conn;
        vols.reserve(el_.size());
        conn.reserve(el_.size());
        for (const auto& e : el_) { vols.push_back(e.V0); conn.push_back(e.n); }
        guards::checkDegenerate("fem3d mesh = file", vols, conn, X0_, gmshNodeId_);
        guards::checkMasses("fem3d mesh = file", m_, X0_, guards::kNoMask,
                            gmshNodeId_);
    }
    if (scen_ == Scenario::TENSION) {
        for (int i = 0; i < (int)X0_.size(); ++i) {
            if (X0_[i].z() < 1e-9)           flag_[i] = FIXED;
            else if (X0_[i].z() > H_ - 1e-9) flag_[i] = PRESCRIBED;
        }
        for (int e = 0; e < (int)el_.size(); ++e) {
            double zc = 0.0;
            for (int a = 0; a < 4; ++a) zc += X0_[el_[e].n[a]].z();
            zc /= 4.0;
            if (zc > H_ / 3.0 && zc < 2.0 * H_ / 3.0) midEl_.push_back(e);
        }
    }
}

void Fem3dSolver::placeTool() {
    if (scen_ == Scenario::TENSION) return;
    tool_.mass   = cfg_.getd("toolMass", 0.5);
    tool_.radius = cfg_.getd("toolRadius", 0.015);
    double gap = cfg_.getd("toolGap", 1e-4);
    // toolShape = sphere | flat ('disc' accepted as the 2D-key synonym of
    // sphere). The flat punch is a vertical cylinder of radius toolRadius
    // whose bottom face indents the top surface — the 3D lift of the 2D
    // FLAT tool. Lateral cutting needs a full 3D normal, so shear forces
    // the sphere, exactly as 2D shear forces the disc.
    std::string sh = cfg_.gets("toolShape", "sphere");
    if (sh != "sphere" && sh != "disc" && sh != "flat" && sh != "blade")
        throw std::runtime_error("toolShape must be sphere | flat | blade "
                                 "(3D ; blade = shear only)");
    tool_.flat = sh == "flat";
    if (scen_ == Scenario::PERCUSSION) {
        if (sh == "blade")
            throw std::runtime_error("toolShape = blade is a SHEAR tool");
        tool_.free = true;
        double zTip = tool_.flat ? H_ + gap : H_ + tool_.radius + gap;
        tool_.x = {cfg_.getd("toolX", 0.5 * W_), cfg_.getd("toolY", 0.5 * D_),
                   zTip};
        tool_.v = {0.0, 0.0, -cfg_.getd("impactSpeed", 8.0)};
    } else {                     // SHEAR: dragged spherical cutter, as fdem3d
        tool_.free = false;
        tool_.flat = false;
        double depth = cfg_.getd("cutDepth", 0.004);
        tool_.blade = sh == "blade";
        if (tool_.blade) {
            // lame rigide plane (2026-09-03) : x = la POINTE, a la profondeur
            // de coupe, hors bloc par defaut (x = -gap) ; face de coupe vers
            // le haut et l'arriere a backRakeDeg de la verticale, face de
            // depouille vers l'arriere a clearanceDeg de l'horizontale
            tool_.rake = cfg_.getd("backRakeDeg", 20.0) * M_PI / 180.0;
            tool_.clear_ = cfg_.getd("clearanceDeg", 10.0) * M_PI / 180.0;
            tool_.hFace = cfg_.getd("bladeHeight", 0.02);
            if (tool_.rake < 0.0 || tool_.rake > 0.5 * M_PI
                || tool_.clear_ <= 0.0 || tool_.clear_ > 0.5 * M_PI
                || tool_.hFace <= 0.0)
                throw std::runtime_error("blade: 0 <= backRakeDeg < 90, "
                                         "0 < clearanceDeg < 90, bladeHeight > 0");
            tool_.x = {cfg_.getd("toolX", -gap), 0.5 * D_, H_ - depth};
        } else {
            tool_.x = {cfg_.getd("toolX", -tool_.radius - gap), 0.5 * D_,
                       H_ - depth + tool_.radius};
        }
        tool_.v = {cfg_.getd("cutSpeed", 10.0), 0.0, 0.0};
    }
    // ---- toolLockXY (2026-09-05, w17, opt-in ; defaut absent = bit-identique)
    // Percussion : l'outil rigide garde v_x = v_y = 0 et x = toolX, y = toolY
    // a tous les pas (Tool3::integrate) ; la force de contact laterale est
    // calculee et comptee (toolFx, peakF_) mais ne le deplace pas — le noeud
    // de reference du bouton bloque en 1,2 des decks Abaqus. Motif : dans le
    // quart de bloc l'insert pose sur l'arete (0, 0) DERIVE vers les plans
    // de symetrie (le contact discretise n'est pas exactement symetrique et
    // rien n'equilibre F_x, F_y : -0,18 mm a 320 us, DECHARGE § 6).
    // Shear : l'outil est pilote en deplacement (v_x = cutSpeed, lame ou
    // sphere trainee) — la cle n'y a pas de sens et est IGNOREE (avertie).
    if (cfg_.getb("toolLockXY", false)) {
        if (scen_ == Scenario::PERCUSSION) {
            tool_.lockXY = true;
            std::cout << "[FEM3D] toolLockXY : outil rigide BLOQUE en x, y (v_x = v_y = 0, "
                         "x = " << tool_.x.x() << " m, y = " << tool_.x.y()
                      << " m fixes a tous les pas) ; la force de contact laterale est "
                         "calculee et comptee (toolFx de history.csv) mais ne deplace "
                         "pas l'outil — equivalent de *Boundary RPB1 1,2 des decks Abaqus\n";
        } else {
            std::cout << "[FEM3D] WARNING: toolLockXY = true IGNORE en scenario shear "
                         "(lame / sphere trainee : l'outil est pilote en deplacement, "
                         "v_x = cutSpeed) — la cle ne s'applique qu'a la percussion\n";
        }
    }
}

// ---------------------------------------------------------------------------
// Confinement suiveur (2026-09-03) — copie de Fdem3dSolver::setupConfinement /
// confiningForces / achievedConfinement sur le registre exterior_ : pression
// sur la geometrie COURANTE des faces exterieures d'origine, lumped A/3 par
// noeud, rampe cosinus. confineFaces = lateral (|n_z| <= 0,5, defaut) | all
// (+ fond : n'agit que si le fond n'est pas encastre, absorbing = all).
// topPressure : pression de boue sur la face z = H (les faces sous l'outil
// aussi — approximation documentee). Les faces des elements erodes ne
// recoivent rien ; les faces y d'une tranche symmetryY non plus.
// Tout est inerte a confiningPressure = topPressure = 0 (bit-identique).
// ---------------------------------------------------------------------------
void Fem3dSolver::setupConfinement() {
    confP_ = cfg_.getd("confiningPressure", 0.0);
    topP_ = cfg_.getd("topPressure", 0.0);
    bottomP_ = cfg_.getd("bottomPressure", 0.0);
    confRamp_ = cfg_.getd("confiningRamp", 30e-6);
    if (confP_ == 0.0 && topP_ == 0.0 && bottomP_ == 0.0) return;
    if (confP_ < 0.0 || topP_ < 0.0 || bottomP_ < 0.0)
        throw std::runtime_error("confiningPressure / topPressure / "
                                 "bottomPressure are POSITIVE pressures");
    std::string cf = cfg_.gets("confineFaces", "lateral");
    if (cf != "lateral" && cf != "all")
        throw std::runtime_error("confineFaces must be lateral | all");
    if (cfg_.gets("absorbing", "none") != "none"
        && cfg_.getd("absorbSpringFactor", 1.0) > 0.0)
        std::cout << "[FEM3D] WARNING: confinement + ressorts de Lysmer "
                     "(absorbSpringFactor > 0) — les ressorts rappellent les "
                     "faces vers u = 0 et avalent une part de la pression ; "
                     "poser absorbSpringFactor = 0\n";
    double tol = 1e-9;
    for (int k = 0; k < (int)exterior_.size(); ++k) {
        const auto& bf = exterior_[k];
        Eigen::Vector3d A = X0_[bf.n[0]], B = X0_[bf.n[1]], C = X0_[bf.n[2]];
        Eigen::Vector3d n = (B - A).cross(C - A);
        double nn = n.norm();
        if (nn < 1e-18) continue;
        double nz = n.z() / nn;
        bool top = nz > 0.5;
        bool bottom = nz < -0.5;
        bool yface = std::abs(n.y() / nn) > 0.5;
        if (top) {
            if (topP_ > 0.0) topFaces_.push_back(k);
            continue;
        }
        if (symY_ && yface) continue;
        // quart de bloc : pas de pression sur les plans de symetrie x = 0 et
        // y = 0 (les faces x = W et y = D la recoivent normalement)
        if (symQ_ && ((A.x() < tol && B.x() < tol && C.x() < tol)
                      || (A.y() < tol && B.y() < tol && C.y() < tol))) continue;
        // bottomPressure : appui sous le bloc (equilibre d'une pression de
        // dessus quand le fond est un amortisseur et non un encastrement)
        if (bottom && bottomP_ > 0.0 && flag_[bf.n[0]] != FIXED)
            bottomFaces_.push_back(k);
        if (bottom && cf != "all") continue;
        if (bottom && A.z() < tol && B.z() < tol && C.z() < tol
            && flag_[bf.n[0]] == FIXED)
            std::cout << "[FEM3D] WARNING: confineFaces = all sur un fond "
                         "ENCASTRE : la pression du fond est sans effet\n";
        if (confP_ > 0.0) confFaces_.push_back(k);
    }
    confGaugeT_ = cfg_.getd("confineGaugeTime", 2.0 * confRamp_);
    std::cout << "[FEM3D] confinement: " << confP_ / 1e6 << " MPa sur "
              << confFaces_.size() << " faces (" << cf << "), pression de "
                 "dessus " << topP_ / 1e6 << " MPa sur " << topFaces_.size()
              << " faces, rampe " << confRamp_ << " s, jauge a "
              << confGaugeT_ << " s, outil gele jusqu'a " << toolDelay_
              << " s\n";
    if (toolDelay_ < confGaugeT_ && scen_ != Scenario::TENSION)
        std::cout << "[FEM3D] WARNING: l'outil frappe avant la jauge de "
                     "confinement (toolDelay < confineGaugeTime)\n";
}

void Fem3dSolver::confiningForces() {
    if (confFaces_.empty() && topFaces_.empty() && bottomFaces_.empty()) return;
    double ramp = 1.0;
    if (confRamp_ > 0.0 && t_ < confRamp_)
        ramp = 0.5 * (1.0 - std::cos(M_PI * t_ / confRamp_));
    auto apply = [&](const std::vector<int>& faces, double p) {
        for (int k : faces) {
            const auto& bf = exterior_[k];
            if (el_[bf.elem].st.eroded) continue;      // face decharge
            Eigen::Vector3d A = X0_[bf.n[0]] + u_[bf.n[0]];
            Eigen::Vector3d B = X0_[bf.n[1]] + u_[bf.n[1]];
            Eigen::Vector3d C = X0_[bf.n[2]] + u_[bf.n[2]];
            Eigen::Vector3d S = 0.5 * (B - A).cross(C - A);   // sortante
            Eigen::Vector3d F = -p * S / 3.0;                 // par noeud
            f_[bf.n[0]] += F;
            f_[bf.n[1]] += F;
            f_[bf.n[2]] += F;
            confWork_ += dt_ * (F.dot(v_[bf.n[0]]) + F.dot(v_[bf.n[1]])
                                + F.dot(v_[bf.n[2]]));
        }
    };
    apply(confFaces_, confP_ * ramp);
    apply(topFaces_, topP_ * ramp);
    apply(bottomFaces_, bottomP_ * ramp);
}

// contrainte laterale moyenne (sxx + syy)/2 dans le coeur (moitie centrale)
double Fem3dSolver::achievedConfinement() const {
    double s = 0.0, V = 0.0;
    for (const auto& e : el_) {
        if (e.st.eroded) continue;
        Eigen::Vector3d c = 0.25 * (X0_[e.n[0]] + X0_[e.n[1]]
                                    + X0_[e.n[2]] + X0_[e.n[3]]);
        if (std::abs(c.x() - 0.5 * W_) > 0.25 * W_
            || std::abs(c.y() - 0.5 * D_) > 0.25 * D_
            || std::abs(c.z() - 0.5 * H_) > 0.25 * H_) continue;
        s += e.slat * e.V0;
        V += e.V0;
    }
    return V > 0.0 ? s / V : 0.0;
}

// masque des noeuds encore portes par au moins un element vivant (copie du
// fem 2D) : un noeud orphelin ne voit plus l'outil (fantome de Saksala 2023)
void Fem3dSolver::refreshActiveNodes() {
    std::fill(active_.begin(), active_.end(), 0);
    for (const auto& e : el_)
        if (!e.st.eroded)
            for (int a = 0; a < 4; ++a) active_[e.n[a]] = 1;
    activeDirty_ = false;
}

void Fem3dSolver::setupBoundaries() {
    cAbs_.assign(X0_.size(), Eigen::Vector3d::Zero());
    kAbs_.assign(X0_.size(), Eigen::Vector3d::Zero());
    if (scen_ == Scenario::TENSION) return;

    std::string ab = cfg_.gets("absorbing", "none");
    if (ab != "none" && ab != "sides" && ab != "all")
        throw std::runtime_error("absorbing must be none | sides | all");

    double G  = mat_.G();
    double sF = cfg_.getd("absorbSpringFactor", 1.0);
    double Rx = cfg_.getd("absorbSpringR", 0.5 * W_);
    double Ry = cfg_.getd("absorbSpringR", 0.5 * D_);
    double Rz = cfg_.getd("absorbSpringR", H_);
    double tol = 1e-9;

    if (cfg_.gets("geometry", "box") == "cylinder") {
        // Curved lateral surface: face normals are oblique, so the
        // diagonal Lysmer is applied with per-component weighting by the
        // geometric normal (|n_a| picks cP, the rest cS) — the standard
        // cheap approximation; a few percent of obliquely incident energy
        // reflects, as on the box faces. Top disc (z = H) stays free, the
        // bottom is absorbed with 'all' or FIXED otherwise.
        double Rc = 0.5 * std::min(W_, D_);
        for (const auto& bf : exterior_) {
            Eigen::Vector3d A = X0_[bf.n[0]], B = X0_[bf.n[1]],
                            C = X0_[bf.n[2]];
            Eigen::Vector3d cen = (A + B + C) / 3.0;
            if (cen.z() > H_ - tol) continue;              // impact surface
            bool bottom = A.z() < tol && B.z() < tol && C.z() < tol;
            if (bottom && ab != "all") {
                for (int nid : bf.n) flag_[nid] = FIXED;
                continue;
            }
            if (!bottom && ab == "none") continue;
            Eigen::Vector3d nrm = (B - A).cross(C - A);
            double A2 = nrm.norm();
            if (A2 < 1e-20) continue;
            double At3 = 0.5 * A2 / 3.0;
            nrm /= A2;
            double R = bottom ? Rz : Rc;
            for (int nid : bf.n)
                for (int a = 0; a < 3; ++a) {
                    double w = std::abs(nrm(a));
                    double c = mat_.cP() * w + mat_.cS() * (1.0 - w);
                    cAbs_[nid](a) += mat_.rho * c * At3;
                    kAbs_[nid](a) += sF * G
                                     / (w * R + (1.0 - w) * 2.0 * R) * At3;
                }
        }
        return;
    }

    for (const auto& bf : exterior_) {
        Eigen::Vector3d A = X0_[bf.n[0]], B = X0_[bf.n[1]], C = X0_[bf.n[2]];
        double At3 = 0.5 * (B - A).cross(C - A).norm() / 3.0;
        bool xlo = A.x() < tol && B.x() < tol && C.x() < tol;
        bool xhi = A.x() > W_ - tol && B.x() > W_ - tol && C.x() > W_ - tol;
        bool ylo = A.y() < tol && B.y() < tol && C.y() < tol;
        bool yhi = A.y() > D_ - tol && B.y() > D_ - tol && C.y() > D_ - tol;
        bool zlo = A.z() < tol && B.z() < tol && C.z() < tol;
        int nAxis = xlo || xhi ? 0 : (ylo || yhi ? 1 : (zlo ? 2 : -1));
        if (nAxis < 0) continue;
        if (symY_ && nAxis == 1) continue;     // plans de symetrie : ni Lysmer
        if (symQ_ && ((nAxis == 0 && xlo) || (nAxis == 1 && ylo))) continue;   // quart de bloc
        double R = nAxis == 0 ? Rx : (nAxis == 1 ? Ry : Rz);
        bool lateral = nAxis != 2;
        for (int nid : bf.n) {
            if (lateral && ab != "none") {
                for (int a = 0; a < 3; ++a) {
                    double c = (a == nAxis ? mat_.cP() : mat_.cS());
                    cAbs_[nid](a) += mat_.rho * c * At3;
                    kAbs_[nid](a) += sF * G / (a == nAxis ? R : 2.0 * R) * At3;
                }
            }
            if (nAxis == 2) {
                if (ab == "all") {
                    for (int a = 0; a < 3; ++a) {
                        double c = (a == 2 ? mat_.cP() : mat_.cS());
                        cAbs_[nid](a) += mat_.rho * c * At3;
                        kAbs_[nid](a) += sF * G / (a == 2 ? R : 2.0 * R) * At3;
                    }
                } else if (!bottomFree_) {     // percussion AND shear: the
                    flag_[nid] = FIXED;        // block needs its support
                }                              // bottomFree : fond libre (bloc flottant, comme un
                                               // deck Abaqus sans condition au fond ; opt-in)
            }
        }
    }
}

void Fem3dSolver::computeStableDt() {
    double cfl = hmin_ / mat_.cP();
    // viscosite de volume (bulkViscosity, 2026-09-05) : Abaqus/Explicit reduit
    // le pas stable de l'element du facteur sqrt(1 + xi^2) - xi, xi = b1 -
    // b2^2 L_e edot/c_d ; on applique la part LINEAIRE xi = b1 (le terme
    // quadratique depend du taux courant : non porte, documente). Ne touche
    // pas la borne du ressort de contact. Inactif (bv_ = false) : cfl intact.
    if (bv_) cfl *= std::sqrt(1.0 + bvB1_ * bvB1_) - bvB1_;
    double dtTool = 1e30;
    for (std::size_t i = 0; i < X0_.size(); ++i)
        if (m_[i] > 0.0)                   // carved grids leave unused nodes
            dtTool = std::min(dtTool, 2.0 * std::sqrt(m_[i] / kp_));
    dt_ = cfg_.getd("dtFactor", 0.3) * std::min(cfl, dtTool);
}

// ===========================================================================

void Fem3dSolver::step() {
    for (auto& fi : f_) fi.setZero();
    tool_.F.setZero();

    elementForces();
    toolContact();
    confiningForces();                     // no-op sans confinement
    if ((confP_ > 0.0 || topP_ > 0.0) && !confLatched_ && t_ >= confGaugeT_) {
        confAchieved_ = achievedConfinement();
        confLatched_ = true;
    }

    if (scen_ == Scenario::TENSION) {
        gripF_.setZero();
        for (int i = 0; i < (int)X0_.size(); ++i)
            if (flag_[i] == PRESCRIBED) gripF_ += f_[i];
        sigmaPeak_ = std::max(sigmaPeak_,
                              std::abs(gripF_.z()) / gripSection_);
        double s = 0.0, vv = 0.0;
        for (int e : midEl_)
            if (!el_[e].st.eroded) {
                s += el_[e].szz * el_[e].V0;
                vv += el_[e].V0;
            }
        sigMid_ = vv > 0 ? s / vv : 0.0;
        sigMidPeak_ = std::max(sigMidPeak_, std::abs(sigMid_));
        if (t_ > 0.75 * T_) {
            sigMidSum_ += std::abs(sigMid_);
            ++sigMidN_;
        }
    } else {
        peakF_ = std::max(peakF_, tool_.F.norm());
        work_ += -tool_.F.dot(tool_.v) * dt_;
        toolFxMax_ = std::max(toolFxMax_, std::abs(tool_.F.x()));   // suivi lateral
        toolFyMax_ = std::max(toolFyMax_, std::abs(tool_.F.y()));   // (resume)
    }

    if (pulseF_ > 0.0 && scen_ != Scenario::TENSION) {
        // impulsion de force imposee : appliquee a la masse de l'outil pour
        // l'integration seulement, la force de contact enregistree (tool_.F,
        // toolFz, peakF_, work_) reste celle du contact
        double a = 0.0;
        if (t_ >= pulseT_.front() && t_ <= pulseT_.back()) {
            std::size_t k = 1;
            while (k < pulseT_.size() && pulseT_[k] < t_) ++k;
            double w = (t_ - pulseT_[k - 1]) / (pulseT_[k] - pulseT_[k - 1]);
            a = pulseA_[k - 1] + w * (pulseA_[k] - pulseA_[k - 1]);
        }
        // force pure : F = -F0 a(t) ; source a impedance : F = -(2 F_inc - Z v_down)
        // avec v_down = -v_z, soit F_z = -2 F0 a - Z v_z (la tige peut tirer :
        // bouton lie a la tige, comme la contrainte *Equation des decks)
        double Fp = pulseZ_ > 0.0 ? -2.0 * pulseF_ * a - pulseZ_ * tool_.v.z()
                                  : -pulseF_ * a;       // vers -z
        pulseImp_ += -Fp * dt_;
        pulseWork_ += Fp * tool_.v.z() * dt_;           // F.v, > 0 quand il pousse (Fp et v_z < 0)
        Eigen::Vector3d Fsave = tool_.F;
        tool_.F.z() += Fp;
        integrate();
        tool_.F = Fsave;
    } else {
        integrate();
    }
    if (scen_ != Scenario::TENSION) {      // suivi lateral de l'outil (resume)
        toolDxMax_ = std::max(toolDxMax_, std::abs(tool_.x.x() - toolX0_.x()));
        toolDyMax_ = std::max(toolDyMax_, std::abs(tool_.x.y() - toolX0_.y()));
    }
    t_ += dt_;
    ++stepCount_;
    // C4 (w20) : detecteur REEL (toutes les composantes de u, v, f de tous
    // les noeuds) tous les nanCheckEvery pas — remplace le test aveugle sur
    // u_[0] (noeud 0 possiblement FIXED, donc toujours fini)
    if (nanEvery_ > 0 && stepCount_ % nanEvery_ == 0) checkFinite();
}

void Fem3dSolver::checkFinite() {
    static const char* const names[3] = {"u", "v", "f"};
    guards::checkFinite("FEM3D", stepCount_, t_, X0_.size(), 3, names,
        [&](std::size_t i, int k) -> const Eigen::Vector3d& {
            return k == 0 ? u_[i] : k == 1 ? v_[i] : f_[i];
        },
        [&](std::size_t i) -> const Eigen::Vector3d& { return X0_[i]; },
        [&](std::size_t i) -> long {
            for (std::size_t e = 0; e < el_.size(); ++e)
                for (int a = 0; a < 4; ++a)
                    if ((std::size_t)el_[e].n[a] == i) return (long)e;
            return -1;
        },
        {{"work", work_}, {"sigmaPeak", sigmaPeak_},
         {"|toolF|", tool_.F.norm()}, {"|toolV|", tool_.v.norm()}});
}

// Co-rotational tets; the law does the physics. Shared nodes make the
// scatter race under OpenMP: per-thread buffers reduced in thread order
// (deterministic per thread count; 1 thread = serial arithmetic).
// relance hors region OpenMP l'exception capturee par un thread (revue cdp
// 2026-09-04) : indice d'element, taille et facteur de Weibull en tete
static void lawError(long el, double lc, double ftScale, const std::string& what) {
    throw std::runtime_error("[fem3d] constitutive law threw in element "
                             + std::to_string(el) + " (lc = " + std::to_string(lc)
                             + " m, ftScale = " + std::to_string(ftScale) + ") : "
                             + what);
}

void Fem3dSolver::elementForces() {
    long nEro = 0, nEroGeo = 0;
    double vEroLaw = 0.0, vEroGeo = 0.0, eRem = 0.0;
    const bool needWel = stats_ || erodeDetMin_ > 0.0 || erodeStrainMax_ > 0.0;
    // Revue cdp du 2026-09-04 : une exception levee par la loi (cdp « no
    // bracket / return mapping failed », dpr « crack-band limit » a ftScale
    // > 1...) DANS la region OpenMP appelait std::terminate — le processus
    // mourait sans message (stdout non vide, code 127), le catch de main.cpp
    // ne voyant rien. Chaque thread capture la premiere exception de ses
    // elements (message, indice, lc, ftScale) dans son accumulateur ; elle
    // est relancee apres la barriere. Sans exception : aucun chemin nouveau.
    struct EroAcc {
        long law = 0, geo = 0; double vLaw = 0, vGeo = 0, eRem = 0;
        long errEl = -1; double errLc = 0.0, errFt = 1.0; std::string err;
        double wBv = 0.0;                      // dissipation de la viscosite de volume (J)
    };
    // viscosite de volume : rho c_d (non endommage) et rho, lus une fois
    const double bvRhoC = bv_ ? mat_.rho * mat_.cP() : 0.0;
    const double bvRho = mat_.rho;
    // sondes de point materiel (probes, 2026-09-05, w18) : sans la cle prb =
    // false et aucun test n'est fait par element ; avec, un element sonde
    // (probeSlot_ >= 0) depose la contrainte de la loi et la deformation qu'elle
    // a recue dans son creneau (un seul thread par element : pas de course).
    // Erosion : sigma = 0, eps = derniere valeur calculee (nullptr = inchangee).
    const bool prb = !probes_.empty();
    auto probeRec = [&](const Elem& e, const Eigen::Matrix3d* sig,
                        const Eigen::Matrix3d* eps) {
        const int sl = probeSlot_[(std::size_t)(&e - el_.data())];
        if (sl < 0) return;
        if (sig) probeSig_[sl] = *sig; else probeSig_[sl].setZero();
        if (eps) probeEps_[sl] = *eps;
    };
    auto processElem = [&](Elem& e, auto&& addF, EroAcc& acc) {
        if (e.st.eroded || acc.errEl >= 0) return;
        Eigen::Matrix3d F = Eigen::Matrix3d::Zero();
        for (int a = 0; a < 4; ++a)
            F += (X0_[e.n[a]] + u_[e.n[a]]) * e.dN.col(a).transpose();
        Eigen::Matrix3d R;
        double det = F.determinant();
        // soupape GEOMETRIQUE (erodeDetMin, 2026-09-03) : un element ecrase
        // sous det F < seuil est retire et compte a part — un mecanisme
        // numerique, pas une physique ; l'energie elastique effacee est
        // comptee (eRemoved). Inerte a erodeDetMin = 0.
        if (erodeDetMin_ > 0.0 && det < erodeDetMin_) {
            e.st.eroded = true;
            e.erodedBy = 2;
            ++acc.geo;
            acc.vGeo += e.V0;
            acc.eRem += e.wEl * e.V0;
            e.J = det;
            e.svm = 0.0; e.pm = 0.0; e.szz = 0.0; e.slat = 0.0; e.wEl = 0.0;
            if (prb) probeRec(e, nullptr, nullptr);
            return;
        }
        // cinematique hencky (kinematics = hencky, 2026-09-05, opt-in) :
        // R polaire exacte, eps = ln U et U^-1 par la decomposition spectrale
        // de F^T F (Fem3dKinematics.hpp) ; hk = false (cle absente, ou
        // element degenere : det F <= 1e-9, lambda_min <= 1e-9) -> chemin
        // historique ci-dessous, inchange
        Eigen::Matrix3d eps, Uinv;
        bool hk = false;
        if (hencky_) hk = fem3dkin::hencky(F, det, R, eps, Uinv);
        if (!hk) {
        if (det > 1e-9) {
            R = F * std::sqrt(3.0) / F.norm();
            for (int it = 0; it < 3; ++it)
                R = 0.5 * (R + R.inverse().transpose());
        } else R.setIdentity();
        Eigen::Matrix3d Ub = R.transpose() * F;
        eps = 0.5 * (Ub + Ub.transpose())
              - Eigen::Matrix3d::Identity();
        }
        // plafond de deformation (erodeStrainMax, le strainCap du fem 2D) :
        // un element ETIRE au-dela du seuil (norme de Frobenius de la
        // deformation de Biot, >= la plus grande principale) est retire et
        // compte avec la soupape — sans lui, les noeuds orphelins d'une
        // premiere erosion etirent leurs voisins en cascade (mini-bloc B3)
        if (erodeStrainMax_ > 0.0 && eps.norm() > erodeStrainMax_) {
            e.st.eroded = true;
            e.erodedBy = 2;
            ++acc.geo;
            acc.vGeo += e.V0;
            acc.eRem += e.wEl * e.V0;
            e.J = det;
            e.svm = 0.0; e.pm = 0.0; e.szz = 0.0; e.slat = 0.0; e.wEl = 0.0;
            if (prb) probeRec(e, nullptr, &eps);
            return;
        }
        const double wPrev = e.wEl;
        Eigen::Matrix3d sig;
        try {
            sig = law_->stress(eps, e.st, dt_, e.lc);
        } catch (const std::exception& ex) {
            acc.errEl = (long)(&e - el_.data());
            acc.errLc = e.lc;
            acc.errFt = e.st.ftScale;
            acc.err = ex.what();
            e.J = det;
            return;
        }
        if (e.st.eroded) {
            ++acc.law;
            acc.vLaw += e.V0;
            acc.eRem += wPrev * e.V0;
            e.erodedBy = 1;
            e.svm = 0.0; e.pm = 0.0; e.szz = 0.0; e.slat = 0.0; e.wEl = 0.0;
            e.J = det;
            if (prb) probeRec(e, nullptr, &eps);
            return;
        }
        if (prb) probeRec(e, &sig, &eps);      // contrainte de la LOI, avant viscosite de volume
        double pm = sig.trace() / 3.0;
        e.pm = pm;
        e.szz = sig(2, 2);
        e.svm = std::sqrt(1.5)
                * (sig - pm * Eigen::Matrix3d::Identity()).norm();
        // jauge laterale : (sxx + syy)/2, ou sxx seul en tranche symmetryY
        // (la deformation plane y impose syy = nu (sxx + szz), pas -P)
        e.slat = symY_ ? sig(0, 0) : 0.5 * (sig(0, 0) + sig(1, 1));
        if (triax_) {                          // jauge triaxiale (Biot ; hencky : Cauchy, ln U)
            e.sxx = sig(0, 0); e.syy = sig(1, 1);
            e.ezz = eps(2, 2); e.ev = eps.trace();
        }
        e.J = det;
        if (needWel)                           // energie elastique nominale
            e.wEl = 0.5 * sig.cwiseProduct(eps - e.st.epsP).sum();
        if (bv_) {
            // viscosite de volume (bulkViscosity, 2026-09-05) : taux
            // volumique de l'element edot = (J_n+1 - J_n)/(dt J_n+1) (J = det F
            // courant, memoire e.Jbv), L_e = lc ; contrainte visqueuse
            // isotrope q I ajoutee a la contrainte de la LOI pour les forces
            // internes seulement (les sorties svm/pm/szz/slat/wEl ci-dessus
            // restent la contrainte materiau). q est du signe de edot : la
            // pression visqueuse s'oppose au taux volumique et la dissipation
            // q edot est >= 0 par construction ; le terme quadratique
            // n'agit qu'en compression volumique (edot < 0).
            const double edot = (det - e.Jbv) / (dt_ * det);
            e.Jbv = det;
            double q = bvB1_ * bvRhoC * e.lc * edot;
            if (edot < 0.0) {
                const double b2Le = bvB2_ * e.lc * edot;
                q -= bvRho * b2Le * b2Le;
            }
            acc.wBv += e.V0 * q * edot * dt_;
            sig(0, 0) += q; sig(1, 1) += q; sig(2, 2) += q;
        }
        Eigen::Matrix3d P = R * sig;
        // hencky : P = J R sigma_c U^-1 (sigma_c = Cauchy co-rotationnelle,
        // viscosite de volume comprise), forces sur la configuration
        // courante avec les gradients de REFERENCE ; biot : P = R sigma_c
        if (hk) P = det * (P * Uinv);
        for (int a = 0; a < 4; ++a)
            addF(e.n[a], -e.V0 * (P * e.dN.col(a)));
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
        std::vector<EroAcc> eroT(nT);
#pragma omp parallel
        {
            int t = omp_get_thread_num();
            auto& fb = fTL_[t];
            auto& seen = seenTL_[t];
            auto& tl = touchedTL_[t];
            tl.clear();
            EroAcc acc;
            auto addF = [&](int i, const Eigen::Vector3d& v3) {
                if (!seen[i]) { seen[i] = 1; tl.push_back(i); }
                fb[i] += v3;
            };
#pragma omp for schedule(static)
            for (int eI = 0; eI < (int)el_.size(); ++eI)
                processElem(el_[eI], addF, acc);
            eroT[t] = acc;
        }
        for (int t = 0; t < nT; ++t) {
            for (int i : touchedTL_[t]) {
                f_[i] += fTL_[t][i];
                fTL_[t][i].setZero();
                seenTL_[t][i] = 0;
            }
            nEro += eroT[t].law;
            nEroGeo += eroT[t].geo;
            vEroLaw += eroT[t].vLaw;
            vEroGeo += eroT[t].vGeo;
            eRem += eroT[t].eRem;
            if (bv_) wBulk_ += eroT[t].wBv;    // ordre des threads : deterministe
        }
        for (int t = 0; t < nT; ++t)
            if (eroT[t].errEl >= 0) lawError(eroT[t].errEl, eroT[t].errLc,
                                             eroT[t].errFt, eroT[t].err);
    } else
#endif
    {
        auto addF = [&](int i, const Eigen::Vector3d& v3) { f_[i] += v3; };
        EroAcc acc;
        for (auto& e : el_) processElem(e, addF, acc);
        nEro = acc.law; nEroGeo = acc.geo;
        vEroLaw = acc.vLaw; vEroGeo = acc.vGeo; eRem = acc.eRem;
        if (bv_) wBulk_ += acc.wBv;
        if (acc.errEl >= 0) lawError(acc.errEl, acc.errLc, acc.errFt, acc.err);
    }
    nEroded_ += nEro + nEroGeo;            // total (loi + soupape)
    nErodedLaw_ += nEro;
    nErodedGeo_ += nEroGeo;
    vErodedLaw_ += vEroLaw;
    vErodedGeo_ += vEroGeo;
    eRemoved_ += eRem;
    if (activeNodes_ && (nEro + nEroGeo) > 0) activeDirty_ = true;
}

void Fem3dSolver::toolContact() {
    if (scen_ == Scenario::TENSION) return;
    if (t_ < toolDelay_) return;               // outil gele (confinement)
    if (activeNodes_ && activeDirty_) refreshActiveNodes();
    for (int i = 0; i < (int)X0_.size(); ++i) {
        if (activeNodes_ && !active_[i]) continue;   // noeud orphelin
        // noeud de masse nulle (noeud 0-D d'un maillage Gmsh jamais reference
        // par un tetraedre, epingle FIXED a la lecture) : il ne porte aucune
        // matiere, il ne doit recevoir aucune force de contact. Avant ce garde
        // (2026-09-05) la penalite lui appliquait kp x penetration — une
        // BROCHE FANTOME sous le pole de l'insert quand le point du champ de
        // taille coincide avec le point d'impact — et la voie Signorini
        // divisait par sa masse (NaN au premier contact). Sans noeud orphelin
        // le resultat est inchange.
        if (m_[i] <= 0.0) continue;
        Eigen::Vector3d p = X0_[i] + u_[i];
        Eigen::Vector3d Fc;
        if (tool_.blade) {
            // lame rigide plane : deux faces planes issues de la pointe.
            // Face de coupe : direction d_r = (-sin a, 0, cos a) (haut et
            // arriere), normale vers la roche n_r = (cos a, 0, sin a).
            // Face de depouille : direction d_c = (-cos b, 0, sin b)
            // (arriere et montante), normale vers la roche n_c =
            // (-sin b, 0, -cos b). Un noeud derriere une face (g < 0) dans
            // l'emprise de cette face est en penetration ; la plus petite
            // penetration gagne pres de la pointe.
            Eigen::Vector3d r = p - tool_.x;
            double sa = std::sin(tool_.rake), ca = std::cos(tool_.rake);
            double sb = std::sin(tool_.clear_), cb = std::cos(tool_.clear_);
            Eigen::Vector3d nr(ca, 0.0, sa), dr(-sa, 0.0, ca);
            Eigen::Vector3d nc(-sb, 0.0, -cb), dc(-cb, 0.0, sb);
            double gr = r.dot(nr), sr = r.dot(dr);
            double gc = r.dot(nc), sc = r.dot(dc);
            // le corps de la lame = derriere le plan de coupe ET derriere le
            // plan de depouille (l'intersection des deux demi-espaces), borne
            // par la hauteur des faces ; un noeud derriere le seul plan de
            // coupe est dans la roche deja coupee, pas dans l'outil
            // le corps est le coin convexe engendre par dr et dc depuis la
            // pointe : coordonnees sr, sc >= 0 (bornees par la hauteur des
            // faces). Face qui repousse : la face de COUPE pour tout noeud
            // au-dessus du niveau de la pointe, la face de DEPOUILLE en
            // dessous — la regle « face la plus proche » envoyait sous
            // l'outil les noeuds pres de la pointe (l'outil enjambait la
            // matiere, B4 a depouille 60 deg).
            bool inside = gr < 0.0 && gc < 0.0 && sr >= 0.0 && sc >= 0.0
                          && sr <= tool_.hFace && sc <= tool_.hFace;
            if (!inside) continue;
            // face qui repousse : la face de COUPE pour un noeud au-dessus du
            // niveau de la pointe, la face de DEPOUILLE en dessous. La regle
            // « face la plus proche » (bissectrice), essayee sur conseil de la
            // revue, envoie sous l'outil les noeuds pres de la pointe et fait
            // perdre le contact a depouille 60 deg (banc B4, w7) ; la regle par
            // hauteur donne Fz/Fx = tan(rake) exact (w6). Limite connue : un
            // noeud entrant par en dessous et remontant au-dessus de la pointe
            // change de face brutalement (rare : la matiere sous l'outil est
            // enlevee par l'erosion).
            Eigen::Vector3d n;
            double pen;
            if (r.z() >= 0.0) { n = nr; pen = -gr; }
            else              { n = nc; pen = -gc; }
            Eigen::Vector3d vrel = v_[i] - tool_.v;
            double c = 2.0 * xiC_ * std::sqrt(kp_ * m_[i]);
            double fn = kp_ * pen - c * vrel.dot(n);
            if (fn < 0) fn = 0;
            Eigen::Vector3d vt = vrel - vrel.dot(n) * n;
            double vtn = vt.norm();
            Eigen::Vector3d ftv = Eigen::Vector3d::Zero();
            if (vtn > 0) ftv = -muC_ * fn * std::tanh(vtn / vReg_) * vt / vtn;
            Fc = fn * n + ftv;
        } else if (tool_.flat) {
            // flat-ended punch: vertical contact only against the bottom
            // face (sharp edge, as in 2D — nodes outside the radius are
            // untouched)
            double rx = p.x() - tool_.x.x(), ry = p.y() - tool_.x.y();
            if (rx * rx + ry * ry > tool_.radius * tool_.radius) continue;
            double pen = p.z() - tool_.x.z();
            if (pen <= 0) continue;
            double c = 2.0 * xiC_ * std::sqrt(kp_ * m_[i]);
            double fn = kp_ * pen + c * (v_[i].z() - tool_.v.z());
            if (fn < 0) fn = 0;
            Eigen::Vector3d vt(v_[i].x() - tool_.v.x(),
                               v_[i].y() - tool_.v.y(), 0.0);
            double vtn = vt.norm();
            Eigen::Vector3d ftv = Eigen::Vector3d::Zero();
            if (vtn > 0) ftv = -muC_ * fn * std::tanh(vtn / vReg_) * vt / vtn;
            Fc = ftv + Eigen::Vector3d(0.0, 0.0, -fn);
        } else if (toolSig_) {
            // sphere, voie SIGNORINI (port de Fdem3dSolver::nodeSig) : la
            // geometrie est dupliquee a dessein, la voie penalite ne bouge pas.
            // Vitesse libre = apres les forces d'elements (toolContact est
            // appele apres elementForces, avant confiningForces : la pression
            // suiveuse des faces laterales n'entre pas dans v*, elle est loin
            // du pole).
            Eigen::Vector3d d = p - tool_.x;
            double dist = d.norm();
            if (dist >= tool_.radius || dist < 1e-14) continue;
            Eigen::Vector3d n = d / dist;
            double pen = tool_.radius - dist;
            Eigen::Vector3d vFree = v_[i] + (dt_ / m_[i]) * f_[i];
            Eigen::Vector3d vrel = vFree - tool_.v;
            double vn = vrel.dot(n);
            Eigen::Vector3d vt = vrel - vn * n;
            double vtn = vt.norm();
            toolsig::Impulse R = toolsig::impulse(pen, vn, vtn, m_[i], dt_, muC_,
                                                  toolSigRelax_);
            if (!R.active) continue;
            Fc = (R.rn / dt_) * n;
            if (vtn > 1e-14) Fc += (R.rt / dt_) * (vt / vtn);
        } else {
            Eigen::Vector3d d = p - tool_.x;
            double dist = d.norm();
            if (dist >= tool_.radius || dist < 1e-14) continue;
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
        f_[i] += Fc;
        tool_.F -= Fc;
    }
}

void Fem3dSolver::integrate() {
    for (int i = 0; i < (int)X0_.size(); ++i) {
        if (flag_[i] == FIXED) {
            if (gripFree_) {                   // hold z only, lateral free
                v_[i].x() += (dt_ / m_[i]) * f_[i].x();
                v_[i].y() += (dt_ / m_[i]) * f_[i].y();
                v_[i].z() = 0.0;
                u_[i] += dt_ * v_[i];
            } else v_[i].setZero();
            continue;
        }
        // pullDelay (2026-09-04) : avant t = pullDelay le noeud du mors est
        // traite comme LIBRE (il tombe dans la branche generale ci-dessous,
        // avec la pression de dessus et l'amortissement) ; a pullDelay = 0
        // la condition est fausse des t = 0 et le chemin est inchange
        if (flag_[i] == PRESCRIBED && !(t_ < pullDelay_)) {
            if (gripFree_) {
                v_[i].x() += (dt_ / m_[i]) * f_[i].x();
                v_[i].y() += (dt_ / m_[i]) * f_[i].y();
            } else { v_[i].x() = 0.0; v_[i].y() = 0.0; }
            double vg = pullV_;
            const double tp = t_ - pullDelay_;     // rampe comptee depuis pullDelay
            if (pullRamp_ > 0.0 && tp < pullRamp_)
                vg *= 0.5 * (1.0 - std::cos(M_PI * tp / pullRamp_));
            v_[i].z() = vg;
            u_[i] += dt_ * v_[i];
            continue;
        }
        for (int a = 0; a < 3; ++a) {
            if (kAbs_[i](a) > 0) f_[i](a) -= kAbs_[i](a) * u_[i](a);
            if (damping_ > 0)
                f_[i](a) -= damping_ * std::abs(f_[i](a))
                            * (v_[i](a) > 0 ? 1.0 : (v_[i](a) < 0 ? -1.0 : 0.0));
        }
        v_[i] += (dt_ / m_[i]) * f_[i];
        for (int a = 0; a < 3; ++a)
            if (cAbs_[i](a) > 0)
                v_[i](a) /= 1.0 + dt_ * cAbs_[i](a) / m_[i];
        if (symY_) v_[i].y() = 0.0;            // tranche plane
        if (symQ_) {                           // quart de bloc : plans de symetrie
            if (X0_[i].x() < 1e-9) v_[i].x() = 0.0;
            if (X0_[i].y() < 1e-9) v_[i].y() = 0.0;
        }
        u_[i] += dt_ * v_[i];
    }
    if (scen_ != Scenario::TENSION && t_ >= toolDelay_) tool_.integrate(dt_);
}

double Fem3dSolver::craterVol() const {
    double v = 0.0;
    for (const auto& e : el_)
        if (e.st.eroded) v += e.V0;
    return v;
}

// ===========================================================================

void Fem3dSolver::writeFrame(int frame) {
    std::vector<Eigen::Vector3d> pts(X0_.size()), vel(X0_.size());
    for (std::size_t i = 0; i < X0_.size(); ++i) {
        pts[i] = X0_[i] + u_[i];
        vel[i] = v_[i];
    }
    std::vector<std::array<int, 4>> tets(el_.size());
    std::vector<double> svm(el_.size()), pm(el_.size()), dmg(el_.size()),
        ero(el_.size()), epv(el_.size()), kdp(el_.size()), fts(el_.size());
    for (std::size_t e = 0; e < el_.size(); ++e) {
        tets[e] = el_[e].n;
        svm[e] = el_[e].svm;
        pm[e] = el_[e].pm;
        dmg[e] = el_[e].st.D;
        ero[e] = el_[e].st.eroded ? 1.0 : 0.0;
        epv[e] = el_[e].st.epvEq;
        kdp[e] = el_[e].st.kappa;
        fts[e] = el_[e].st.ftScale;
    }
    char name[64];
    std::snprintf(name, sizeof(name), "/fem3d_%04d.vtu", frame);
    if (fixedDam_) {
        // tensionDamage = fixed (2026-09-05) : les trois endommagements
        // directionnels d_i du repere fige s'AJOUTENT (dFix1..3) aux champs du
        // mode courant (stats / vtkCap honores) ; « damage » reste max_i d_i.
        // Branche jamais atteinte sans la cle : appels d'origine intacts.
        std::vector<double> d1(el_.size()), d2(el_.size()), d3(el_.size());
        for (std::size_t e = 0; e < el_.size(); ++e) {
            d1[e] = el_[e].st.fcm.d[0];
            d2[e] = el_[e].st.fcm.d[1];
            d3[e] = el_[e].st.fcm.d[2];
        }
        vtk::ScalarField sf{{"vonMises", &svm}, {"pressure", &pm},
                            {"damage", &dmg}, {"eroded", &ero}, {"epvEq", &epv},
                            {"kapDP", &kdp}, {"ftScale", &fts},
                            {"dFix1", &d1}, {"dFix2", &d2}, {"dFix3", &d3}};
        std::vector<double> cpc, evp, omc, detF, eby, slat;
        if (vtkCap_) {
            cpc.resize(el_.size()); evp.resize(el_.size());
            for (std::size_t e = 0; e < el_.size(); ++e) {
                cpc[e] = el_[e].st.pc;
                evp[e] = -el_[e].st.epsP.trace();
            }
            sf["capPc"] = &cpc; sf["epsVpl"] = &evp;
        }
        if (stats_) {
            omc.resize(el_.size()); detF.resize(el_.size());
            eby.resize(el_.size()); slat.resize(el_.size());
            for (std::size_t e = 0; e < el_.size(); ++e) {
                omc[e] = el_[e].st.Dc;
                detF[e] = el_[e].J;
                eby[e] = el_[e].erodedBy == 2 ? 5.0 : (double)el_[e].st.eroCode;
                slat[e] = el_[e].slat;
            }
            sf["omegaC"] = &omc; sf["detF"] = &detF;
            sf["erodedBy"] = &eby; sf["sigLat"] = &slat;
        }
        vtk::writeTetMesh(out_ + name, pts, tets, sf, {{"velocity", &vel}});
    } else
    if (vtkCap_) {
        // vtkCap = true (revue du cap cdp, 2026-09-04 nuit) : la zone
        // compactee etait INVISIBLE (pc et eps_v^pl ne sont dans aucun champ,
        // seul le resume max cap pc les lit). Deux champs cellulaires
        // s'AJOUTENT a la liste du mode courant : capPc = MatState::pc [Pa]
        // (0 sans cap : dpr sans capP0, cdp sans cdpCap) et epsVpl =
        // -tr(eps_pl) (compaction volumique plastique > 0, dilatance < 0).
        // Valable pour cdp ET dpr/saksala. Sans la cle : appels d'origine intacts.
        std::vector<double> cpc(el_.size()), evp(el_.size());
        for (std::size_t e = 0; e < el_.size(); ++e) {
            cpc[e] = el_[e].st.pc;
            evp[e] = -el_[e].st.epsP.trace();
        }
        vtk::ScalarField sf{{"vonMises", &svm}, {"pressure", &pm},
                            {"damage", &dmg}, {"eroded", &ero}, {"epvEq", &epv},
                            {"kapDP", &kdp}, {"ftScale", &fts},
                            {"capPc", &cpc}, {"epsVpl", &evp}};
        std::vector<double> omc, detF, eby, slat;
        if (stats_) {
            omc.resize(el_.size()); detF.resize(el_.size());
            eby.resize(el_.size()); slat.resize(el_.size());
            for (std::size_t e = 0; e < el_.size(); ++e) {
                omc[e] = el_[e].st.Dc;
                detF[e] = el_[e].J;
                eby[e] = el_[e].erodedBy == 2 ? 5.0 : (double)el_[e].st.eroCode;
                slat[e] = el_[e].slat;
            }
            sf["omegaC"] = &omc; sf["detF"] = &detF;
            sf["erodedBy"] = &eby; sf["sigLat"] = &slat;
        }
        vtk::writeTetMesh(out_ + name, pts, tets, sf, {{"velocity", &vel}});
    } else
    if (stats_) {
        // champs de l'etude briques (fichier VTU inchange sans fieldStats)
        std::vector<double> omc(el_.size()), detF(el_.size()), eby(el_.size()),
            slat(el_.size());
        for (std::size_t e = 0; e < el_.size(); ++e) {
            omc[e] = el_[e].st.Dc;
            detF[e] = el_[e].J;
            // 0 vivant, 1 spall, 2 broyage, 3 omega_c, 4 dfh, 5 soupape det F
            eby[e] = el_[e].erodedBy == 2 ? 5.0 : (double)el_[e].st.eroCode;
            slat[e] = el_[e].slat;
        }
        vtk::writeTetMesh(out_ + name, pts, tets,
                          {{"vonMises", &svm}, {"pressure", &pm},
                           {"damage", &dmg}, {"eroded", &ero}, {"epvEq", &epv},
                           {"kapDP", &kdp}, {"ftScale", &fts},
                           {"omegaC", &omc}, {"detF", &detF},
                           {"erodedBy", &eby}, {"sigLat", &slat}},
                          {{"velocity", &vel}});
    } else
    vtk::writeTetMesh(out_ + name, pts, tets,
                      {{"vonMises", &svm}, {"pressure", &pm},
                       {"damage", &dmg}, {"eroded", &ero}, {"epvEq", &epv},
                       {"kapDP", &kdp}, {"ftScale", &fts}},
                      {{"velocity", &vel}});

    std::ofstream fm(out_ + "/frames.csv",
                     frame == 0 ? std::ios::trunc : std::ios::app);
    if (frame == 0) fm << "frame,t,toolX,toolY,toolZ\n";
    fm << frame << "," << t_ << "," << tool_.x.x() << "," << tool_.x.y()
       << "," << tool_.x.z() << "\n";
}

// Colonnes de l'etude briques (2026-09-03), ajoutees EN FIN DE LIGNE et
// SEULEMENT si la cle correspondante est posee : le fichier history.csv
// reste bit-identique sinon (extracteurs positionnels).
static void extraHeader(std::ostream& os, bool conf, bool stats) {
    if (conf) os << ",confP,confAch,confWork";
    if (stats) os << ",keBlock,wPlas,wDamT,wDamC,eRemoved,V_D09,V_wc05,"
                     "V_eroLaw,V_eroGeo,detFmin,pMin,slatMean,nEroLaw,nEroGeo,"
                     "nEroSpall,nEroCrush,nEroWc";
}

void Fem3dSolver::historyHeader(std::ostream& os) const {
    const bool conf = confP_ > 0.0 || topP_ > 0.0;
    if (scen_ == Scenario::TENSION) {
        os << "t,sigma,sigmaPeak,nEroded";
        extraHeader(os, conf, stats_);
        if (triax_) os << ",sigZZmid,sigXXmid,sigYYmid,epsAxMid,epsVolMid,epsAxGrip";
        if (bv_) os << ",wBulk";               // viscosite de volume : derniere colonne
        os << "\n";
        return;
    }
    // toolFx/toolX appended at the END so percussion post-processing that
    // reads columns by position keeps working; they are the cutting force
    // and advance of the shear scenario
    os << "t,toolFz,toolZ,toolVz,work,toolKE,nEroded,craterVol,"
          "toolFx,toolX";
    extraHeader(os, conf, stats_);
    if (triax_) os << ",sigZZmid,sigXXmid,sigYYmid,epsAxMid,epsVolMid,epsAxGrip";
    if (bv_) os << ",wBulk";                   // viscosite de volume : derniere colonne
    os << "\n";
}

void Fem3dSolver::historyRow(std::ostream& os) const {
    const bool conf = confP_ > 0.0 || topP_ > 0.0;
    if (scen_ == Scenario::TENSION) {
        os << t_ << "," << std::abs(gripF_.z()) / gripSection_ << ","
           << sigmaPeak_ << "," << nEroded_;
    } else {
        os << t_ << "," << tool_.F.z() << "," << tool_.x.z() << ","
           << tool_.v.z() << "," << work_ << "," << tool_.ke() << ","
           << nEroded_ << "," << craterVol() << ","
           << tool_.F.x() << "," << tool_.x.x();
    }
    if (conf) {
        double ramp = 1.0;
        if (confRamp_ > 0.0 && t_ < confRamp_)
            ramp = 0.5 * (1.0 - std::cos(M_PI * t_ / confRamp_));
        os << "," << confP_ * ramp << "," << confAchieved_ << "," << confWork_;
    }
    if (stats_) {
        double ke = 0.0;
        for (std::size_t i = 0; i < X0_.size(); ++i)
            ke += 0.5 * m_[i] * v_[i].squaredNorm();
        double wP = 0.0, wT = 0.0, wC = 0.0, vD = 0.0, vC = 0.0;
        double jmin = 1e300, pmin = 0.0, sl = 0.0, vAlive = 0.0;
        long nSp = 0, nCr = 0, nWc = 0;
        for (const auto& e : el_) {
            wP += e.st.wPlas * e.V0;
            wT += e.st.wDamT * e.V0;
            wC += e.st.wDamC * e.V0;
            if (e.st.eroded) {
                if (e.st.eroCode == 1) ++nSp;
                else if (e.st.eroCode == 2) ++nCr;
                else if (e.st.eroCode == 3) ++nWc;
                continue;
            }
            if (e.st.D >= 0.9) vD += e.V0;
            if (e.st.Dc >= 0.5) vC += e.V0;
            jmin = std::min(jmin, e.J);
            pmin = std::min(pmin, e.pm);
            sl += e.slat * e.V0;
            vAlive += e.V0;
        }
        os << "," << ke << "," << wP << "," << wT << "," << wC << ","
           << eRemoved_ << "," << vD << "," << vC << "," << vErodedLaw_
           << "," << vErodedGeo_ << "," << (jmin < 1e299 ? jmin : 0.0)
           << "," << pmin << "," << (vAlive > 0 ? sl / vAlive : 0.0) << ","
           << nErodedLaw_ << "," << nErodedGeo_ << "," << nSp << "," << nCr
           << "," << nWc;
    }
    if (triax_) {
        // jauge du tiers central (moyennes ponderees par V0, elements
        // vivants) + deformation du mors u_z/H (moyenne des noeuds prescrits)
        double szz = 0.0, sxx = 0.0, syy = 0.0, ezz = 0.0, ev = 0.0, V = 0.0;
        for (int e : midEl_) {
            const Elem& el = el_[e];
            if (el.st.eroded) continue;
            szz += el.szz * el.V0; sxx += el.sxx * el.V0; syy += el.syy * el.V0;
            ezz += el.ezz * el.V0; ev += el.ev * el.V0; V += el.V0;
        }
        double uz = 0.0;
        long n = 0;
        for (int i = 0; i < (int)X0_.size(); ++i)
            if (flag_[i] == PRESCRIBED) { uz += u_[i].z(); ++n; }
        const double iv = V > 0.0 ? 1.0 / V : 0.0;
        os << "," << szz * iv << "," << sxx * iv << "," << syy * iv << ","
           << ezz * iv << "," << ev * iv << ","
           << (n > 0 ? uz / (double)n / H_ : 0.0);
    }
    if (bv_) os << "," << wBulk_;              // dissipation cumulee de la viscosite de volume (J)
    os << "\n";
    // sondes de point materiel : une ligne de probes.csv a la meme cadence
    // (history.csv lui-meme n'est pas touche ; rien sans la cle probes)
    if (probesOut_) probesRow(*probesOut_);
}

// ---------------------------------------------------------------------------
// Sondes de point materiel `probes = x,y,z ; x,y,z ; ...` (2026-09-05, w18,
// opt-in ; voir Fem3dSolver.hpp). Analyse stricte de la cle (trois nombres par
// point, separes par des virgules, points separes par des point-virgules ;
// virgule decimale = erreur), localisation barycentrique a l'init, en-tete
// explicite p<k>_<champ>, une ligne par ligne d'historique.
// ---------------------------------------------------------------------------
void Fem3dSolver::setupProbes() {
    if (!cfg_.has("probes")) return;
    const std::string raw = cfg_.gets("probes", "");
    auto num = [&](std::string s, int k, int comp) {
        auto a = s.find_first_not_of(" \t"), b = s.find_last_not_of(" \t");
        s = a == std::string::npos ? "" : s.substr(a, b - a + 1);
        std::size_t used = 0;
        double d = 0.0;
        try { d = std::stod(s, &used); } catch (const std::exception&) { used = 0; }
        if (s.empty() || used != s.size() || !std::isfinite(d))
            throw std::runtime_error("probes: point " + std::to_string(k)
                                     + ", component " + std::to_string(comp + 1)
                                     + " is not a number: '" + s + "' (expected "
                                     "'x,y,z ; x,y,z ; ...' in metres, decimal POINT)");
        return d;
    };
    {
        std::stringstream ss(raw);
        std::string pt;
        int k = 0;
        while (std::getline(ss, pt, ';')) {
            if (pt.find_first_not_of(" \t") == std::string::npos) continue;
            ++k;
            std::vector<std::string> comp;
            std::stringstream sc(pt);
            std::string c;
            while (std::getline(sc, c, ',')) comp.push_back(c);
            if (comp.size() != 3)
                throw std::runtime_error("probes: point " + std::to_string(k)
                                         + " has " + std::to_string(comp.size())
                                         + " components, expected 3 (x,y,z in metres)");
            Probe p;
            for (int i = 0; i < 3; ++i) p.x(i) = num(comp[i], k, i);
            probes_.push_back(p);
        }
    }
    if (probes_.empty())
        throw std::runtime_error("probes: key present but no point given "
                                 "(expected 'x,y,z ; x,y,z ; ...' in metres)");

    // localisation : coordonnees barycentriques lam_a (a = 1..3) = dN_a . (x -
    // X0[n0]) — les colonnes de dN sont les lignes de J^-1 —, lam_0 = 1 - somme ;
    // dedans si toutes >= -tol. Premier element trouve (point sur une face
    // partagee : l'un des deux, deterministe = ordre du maillage). Sinon le
    // centroide le plus proche, avec avertissement.
    const double tol = 1e-9;
    probeSlot_.assign(el_.size(), -1);
    int nSlots = 0;
    for (std::size_t k = 0; k < probes_.size(); ++k) {
        Probe& p = probes_[k];
        int best = -1;
        double bestD = 1e300;
        bool found = false;
        for (std::size_t ei = 0; ei < el_.size(); ++ei) {
            const Elem& e = el_[ei];
            const Eigen::Vector3d r = p.x - X0_[e.n[0]];
            const double l1 = e.dN.col(1).dot(r), l2 = e.dN.col(2).dot(r),
                         l3 = e.dN.col(3).dot(r), l0 = 1.0 - l1 - l2 - l3;
            if (l0 >= -tol && l1 >= -tol && l2 >= -tol && l3 >= -tol) {
                best = (int)ei;
                found = true;
                break;
            }
            const double d = (p.x - e.st.x0).norm();
            if (d < bestD) { bestD = d; best = (int)ei; }
        }
        if (best < 0) throw std::runtime_error("probes: empty mesh");
        p.inside = found;                      // false : aucun tetraedre ne contient le point
        p.dist = (p.x - el_[best].st.x0).norm();
        p.elem = best;
        if (probeSlot_[best] < 0) probeSlot_[best] = nSlots++;
        p.slot = probeSlot_[best];
        const bool outBox = p.x.x() < -tol || p.x.x() > W_ + tol || p.x.y() < -tol
                            || p.x.y() > D_ + tol || p.x.z() < -tol || p.x.z() > H_ + tol;
        std::cout << "[FEM3D] probe " << k + 1 << " : (" << p.x.x() << ", " << p.x.y()
                  << ", " << p.x.z() << ") m -> element " << p.elem << " (centroide "
                  << el_[best].st.x0.x() << ", " << el_[best].st.x0.y() << ", "
                  << el_[best].st.x0.z() << ", lc = " << el_[best].lc << " m, distance "
                  << p.dist << " m)" << (p.inside ? "" : " [PLUS PROCHE, point hors maillage]")
                  << "\n";
        if (!p.inside)
            std::cout << "[FEM3D] WARNING: probe " << k + 1 << " lies in no tetrahedron"
                      << (outBox ? " (OUTSIDE the block [0,W]x[0,D]x[0,H] = [0,"
                                   + std::to_string(W_) + "]x[0," + std::to_string(D_)
                                   + "]x[0," + std::to_string(H_) + "] m)"
                                 : " (inside the box but in no element : hole or face gap)")
                      << " — nearest centroid taken, element " << p.elem << " at "
                      << p.dist << " m\n";
    }
    for (std::size_t k = 0; k < probes_.size(); ++k)
        for (std::size_t j = 0; j < k; ++j)
            if (probes_[j].elem == probes_[k].elem) {
                std::cout << "[FEM3D] probes " << j + 1 << " and " << k + 1
                          << " share element " << probes_[k].elem
                          << " (identical columns)\n";
                break;
            }
    probeSig_.assign(nSlots, Eigen::Matrix3d::Zero());
    probeEps_.assign(nSlots, Eigen::Matrix3d::Zero());
    probesOut_ = std::make_unique<std::ofstream>(out_ + "/probes.csv");
    if (!*probesOut_)
        throw std::runtime_error("probes: cannot open '" + out_ + "/probes.csv'");
    probesOut_->precision(9);
    probesHeader(*probesOut_);
    probesOut_->flush();
    std::cout << "[FEM3D] probes : " << probes_.size() << " sonde(s) dans " << nSlots
              << " element(s) distinct(s) -> " << out_ << "/probes.csv (une ligne par "
                 "ligne d'historique ; contrainte de la loi avant viscosite de volume, "
                 "deformation " << (hencky_ ? "ln U" : "de Biot") << ", repere "
                 "co-rotationnel ; cisaillements tensoriels)"
              << (law_->name() == "cdp" ? " + colonnes cdp epsTpl/epsCpl/d" : "")
              << (law_->name() == "dpdfh" ? " + colonnes dpdfh Dv1..3" : "") << "\n";
}

void Fem3dSolver::probesHeader(std::ostream& os) const {
    static const char* base[] = {
        "elem", "sxx", "syy", "szz", "sxy", "syz", "sxz",
        "exx", "eyy", "ezz", "exy", "eyz", "exz",
        "p", "q", "s1", "s2", "s3", "trEpl", "eplEq", "epvEq",
        "dt", "dc", "pc", "detF", "eroded"};
    const std::string law = law_->name();
    os << "t";
    for (std::size_t k = 0; k < probes_.size(); ++k) {
        const std::string pre = ",p" + std::to_string(k + 1) + "_";
        for (const char* c : base) os << pre << c;
        if (law == "cdp") os << pre << "epsTpl" << pre << "epsCpl" << pre << "dcdp";
        if (law == "dpdfh") os << pre << "Dv1" << pre << "Dv2" << pre << "Dv3";
        if (fixedDam_) os << pre << "dFix1" << pre << "dFix2" << pre << "dFix3";
    }
    os << "\n";
}

void Fem3dSolver::probesRow(std::ostream& os) const {
    const std::string law = law_->name();
    os << t_;
    for (const Probe& p : probes_) {
        const Elem& e = el_[p.elem];
        const Eigen::Matrix3d& s = probeSig_[p.slot];
        const Eigen::Matrix3d& ep = probeEps_[p.slot];
        const double pm = s.trace() / 3.0;
        const Eigen::Matrix3d dev = s - pm * Eigen::Matrix3d::Identity();
        const double q = std::sqrt(1.5) * dev.norm();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(s, Eigen::EigenvaluesOnly);
        Eigen::Vector3d pr = es.info() == Eigen::Success ? es.eigenvalues()
                                                         : Eigen::Vector3d::Zero();
        // eps_pl equivalente sqrt(2/3 dev(eps_pl):dev(eps_pl))
        const Eigen::Matrix3d epd = e.st.epsP - (e.st.epsP.trace() / 3.0)
                                                    * Eigen::Matrix3d::Identity();
        os << "," << p.elem
           << "," << s(0, 0) << "," << s(1, 1) << "," << s(2, 2)
           << "," << s(0, 1) << "," << s(1, 2) << "," << s(0, 2)
           << "," << ep(0, 0) << "," << ep(1, 1) << "," << ep(2, 2)
           << "," << ep(0, 1) << "," << ep(1, 2) << "," << ep(0, 2)
           << "," << -pm << "," << q
           << "," << pr(2) << "," << pr(1) << "," << pr(0)        // decroissantes
           << "," << e.st.epsP.trace() << "," << std::sqrt(2.0 / 3.0) * epd.norm()
           << "," << e.st.epvEq
           << "," << e.st.D << "," << e.st.Dc << "," << e.st.pc
           << "," << e.J << "," << (e.st.eroded ? 1 : 0);
        if (law == "cdp")
            os << "," << e.st.cdp.epsTpl << "," << e.st.cdp.epsCpl << "," << e.st.cdp.d;
        if (law == "dpdfh")
            os << "," << e.st.dfh.Dv[0] << "," << e.st.dfh.Dv[1] << "," << e.st.dfh.Dv[2];
        if (fixedDam_)
            os << "," << e.st.fcm.d[0] << "," << e.st.fcm.d[1] << "," << e.st.fcm.d[2];
    }
    os << "\n";
    os.flush();
}

void Fem3dSolver::finalize() {
    if (nanEvery_ > 0) checkFinite();          // C4 (w20) : dernier controle
    double keBlock = 0.0;
    for (std::size_t i = 0; i < X0_.size(); ++i)
        keBlock += 0.5 * m_[i] * v_[i].squaredNorm();
    double pcMax = 0.0, epvMax = 0.0, pMin = 0.0;
    for (const auto& e : el_) {
        pcMax = std::max(pcMax, e.st.pc);
        epvMax = std::max(epvMax, e.st.epvEq);
        pMin = std::min(pMin, e.pm);
    }
    std::cout << "\n[FEM3D] ---- summary (law = " << law_->name() << ") ----\n"
              << "[FEM3D] block kinetic energy at end: " << keBlock << " J\n"
              << "[FEM3D] eroded elements: " << nEroded_ << " / " << el_.size()
              << " (crater vol " << craterVol() << " m^3)\n"
              << "[FEM3D] max equiv viscoplastic strain: " << epvMax
              << ", min pressure: " << pMin / 1e6 << " MPa, max cap pc: "
              << pcMax / 1e6 << " MPa\n";
    // etude briques : lignes supplementaires seulement si une cle est posee
    if (confP_ > 0.0 || topP_ > 0.0)
        std::cout << "[FEM3D] confinement: vise " << confP_ / 1e6
                  << " MPa (dessus " << topP_ / 1e6 << "), jauge coeur "
                  << confAchieved_ / 1e6 << " MPa ("
                  << (confP_ > 0.0 ? 100.0 * std::abs(confAchieved_ / -confP_) : 0.0)
                  << " % a t = " << confGaugeT_ << " s), travail "
                  << confWork_ << " J\n";
    if (erodeDetMin_ > 0.0 || erodeStrainMax_ > 0.0 || activeNodes_ || stats_) {
        double wP = 0.0, wT = 0.0, wC = 0.0;
        for (const auto& e : el_) {
            wP += e.st.wPlas * e.V0;
            wT += e.st.wDamT * e.V0;
            wC += e.st.wDamC * e.V0;
        }
        std::cout << "[FEM3D] erosion: loi " << nErodedLaw_ << " el. ("
                  << vErodedLaw_ * 1e9 << " mm^3), soupape geometrique (det F < "
                  << erodeDetMin_ << " ou |eps| > " << erodeStrainMax_ << ") : "
                  << nErodedGeo_ << " el. (" << vErodedGeo_ * 1e9
                  << " mm^3), energie elastique effacee " << eRemoved_ << " J\n"
                  << "[FEM3D] dissipation: plastique " << wP
                  << " J, endommagement traction " << wT
                  << " J, compression " << wC << " J\n";
    }

    if (bv_)
        std::cout << "[FEM3D] viscosite de volume (b1 = " << bvB1_ << ", b2 = " << bvB2_
                  << ") : dissipation wBulk = " << wBulk_ << " J (>= 0 attendu)\n";

    if (scen_ == Scenario::TENSION) {
        bool comp = pullV_ < 0.0;
        // verifications read the MID-SPECIMEN stress gauge: the grip force
        // additionally carries the Cundall-damping drag of the flowing
        // column (+11 % measured at damping 0.7) and would fail the sharp
        // bands for reasons that have nothing to do with the law
        std::cout << "[FEM3D] peak |sigma| grip / mid-specimen = "
                  << sigmaPeak_ / 1e6 << " / " << sigMidPeak_ / 1e6
                  << " MPa (" << (comp ? "compression" : "tension") << ")\n";
        if (law_->name() == "dpr") {
            double ref = comp ? law_->sigmaCdp() : mat_.ft;
            double err = 100.0 * (sigMidPeak_ - ref) / ref;
            bool pass = std::abs(err) < 5.0;
            std::cout << "[FEM3D]   reference ("
                      << (comp ? "DP cone, analytic" : "ft") << ") = "
                      << ref / 1e6 << " MPa, deviation = " << err
                      << " %  [" << (pass ? "PASS" : "FAIL")
                      << "] (band 5 %, mid gauge)\n";
        } else if (law_->name() == "saksala" && comp) {
            double epdot = std::abs(pullV_) / H_;
            double pred = law_->viscousOverstress(epdot);
            double plateau = sigMidN_ > 0 ? sigMidSum_ / sigMidN_
                                          : std::abs(sigMid_);
            double meas = plateau - law_->sigmaCdp();
            double r = pred > 0 ? meas / pred : 0.0;
            bool pass = r > 0.75 && r < 1.25;
            std::cout << "[FEM3D]   viscous overstress (mid gauge, last-"
                         "quarter average): measured " << meas / 1e6
                      << " MPa, Perzyna predicts " << pred / 1e6
                      << " MPa at epdot = " << epdot << " /s  -> ratio "
                      << r << "  [" << (pass ? "PASS" : "FAIL")
                      << "] (band 25 %)\n";
        }
        return;
    }
    if (pulseF_ > 0.0)
        std::cout << "[FEM3D] impulsion de force sur l'outil : int F dt = " << pulseImp_
                  << " N s, travail int F v dt = " << pulseWork_ << " J, KE outil finale "
                  << tool_.ke() << " J\n";
    std::cout << "[FEM3D] tool lateral    : x = " << tool_.x.x() << ", y = "
              << tool_.x.y() << " m (depart " << toolX0_.x() << ", " << toolX0_.y()
              << "), max |dx| = " << toolDxMax_ << ", max |dy| = " << toolDyMax_
              << " m, max |Fx| = " << toolFxMax_ << ", max |Fy| = " << toolFyMax_
              << " N" << (tool_.lockXY ? " [toolLockXY actif : x, y fixes]" : "")
              << "\n";
    std::cout << "[FEM3D] peak tool force : " << peakF_ << " N\n"
              << "[FEM3D] tool work       : " << work_ << " J  (tool KE loss: "
              << toolKE0_ - tool_.ke() << " J)\n";
    double cv = craterVol();
    if (cv > 0)
        std::cout << "[FEM3D] specific energy : " << work_ / cv << " J/m^3\n";
}

} // namespace rockim
