#include "rockim/DemSolver.hpp"
#include "rockim/VtkWriter.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <random>

namespace rockim {

static inline uint64_t pairKey(int i, int j) {
    if (i > j) std::swap(i, j);
    return (uint64_t(uint32_t(i)) << 32) | uint32_t(j);
}

// ===========================================================================
// Construction / setup
// ===========================================================================

DemSolver::DemSolver(const Config& cfg, std::string outDir)
    : cfg_(cfg), out_(std::move(outDir)) {}

void DemSolver::init() {
    mat_ = Material::from(cfg_);
    W_ = cfg_.getd("W", W_);
    H_ = cfg_.getd("H", H_);
    r_ = cfg_.getd("particleRadius", r_);
    T_ = cfg_.getd("T", 250e-6);
    packing_ = cfg_.gets("packing", "hex");

    std::string s = cfg_.gets("scenario", "percussion");
    if      (s == "percussion") scen_ = Scenario::PERCUSSION;
    else if (s == "shear")      scen_ = Scenario::SHEAR;
    else if (s == "tension")    scen_ = Scenario::TENSION;
    else throw std::runtime_error("DemSolver: unknown scenario '" + s + "'");

    // ---- micro parameters -------------------------------------------------
    // Contact: linear springs. kn ~ E * t makes the grain-scale contact
    // compliance consistent with the elastic modulus (standard 2D heuristic;
    // quantitative macro elasticity requires calibration, see README).
    knC_ = cfg_.getd("knFactor", 1.0) * mat_.E * thk_;
    ksRatio_ = cfg_.getd("ksRatio", ksRatio_);
    ksC_ = ksRatio_ * knC_;
    mu_  = cfg_.getd("contactMu", mu_);
    xiC_ = cfg_.getd("contactXi", xiC_);

    // Parallel bonds: per-area stiffness kbn = E_bond / L0 (PFC-style bond
    // modulus), shear = ksRatio * normal. Strengths default to the macro
    // strengths (exactly recoverable on an aligned square lattice, cf. the
    // tension verification; on a hex lattice they act as micro-parameters).
    tanPhiB_ = std::tan(cfg_.getd("bondFrictionDeg", mat_.phiDeg) * M_PI / 180.0);
    lambda_  = cfg_.getd("bondRadiusFactor", 1.0);
    damping_ = cfg_.getd("dampingLocal", scen_ == Scenario::TENSION ? 0.7 : 0.02);

    bottomWall_ = cfg_.getb("bottomWall", scen_ != Scenario::TENSION);
    sideWalls_  = cfg_.getb("sideWalls", false);

    buildPacking();
    buildBonds();
    placeTool();
    computeStableDt();

    // ------------------------------------------------------------------
    // Lysmer-Kuhlemeyer quiet boundaries for the bonded assembly: viscous
    // dashpots on the outermost particle layer, with the continuum
    // impedances (rho c_p normal to the face, rho c_s tangential) lumped
    // over each boundary particle's tributary length ~ 2r. 'sides' damps
    // the lateral columns (the bottom wall stays), 'all' also damps the
    // bottom row and removes the rigid bottom wall it replaces. Applied
    // implicitly in integrate(), so no extra time-step restriction.
    // ------------------------------------------------------------------
    cAbsX_.assign(p_.size(), 0.0);
    cAbsY_.assign(p_.size(), 0.0);
    kAbsX_.assign(p_.size(), 0.0);
    kAbsY_.assign(p_.size(), 0.0);
    xAnchor_.resize(p_.size());
    for (std::size_t i = 0; i < p_.size(); ++i) xAnchor_[i] = p_[i].x;
    std::string ab = cfg_.gets("absorbing", "none");
    if (scen_ != Scenario::TENSION && ab != "none") {
        if (ab != "sides" && ab != "all")
            throw std::runtime_error("absorbing must be none | sides | all");
        if (ab == "all") bottomWall_ = false;
        // Viscous-spring (Deeks-Randolph) quiet boundary: Lysmer dashpots for
        // the waves + a soft spring to the initial position restoring the
        // static stiffness ~ G/R of the truncated half-space (a pure dashpot
        // has none, letting the block bend as if free-floating under the slow
        // part of the load). Set absorbSpringFactor = 0 for pure Lysmer.
        double layer = cfg_.getd("absorbLayer", 2.2) * r_;
        double G  = mat_.E / (2.0 * (1.0 + mat_.nu));
        double sF = cfg_.getd("absorbSpringFactor", 1.0);
        double Rside = cfg_.getd("absorbSpringR", 0.5 * W_);
        double Rbot  = cfg_.getd("absorbSpringR", H_);
        for (std::size_t i = 0; i < p_.size(); ++i) {
            double At = 2.0 * p_[i].r * thk_;
            double zP = mat_.rho * mat_.cP() * At;
            double zS = mat_.rho * mat_.cS() * At;
            double x = p_[i].x.x(), y = p_[i].x.y();
            if (x < layer || x > W_ - layer) {
                cAbsX_[i] += zP;                cAbsY_[i] += zS;
                kAbsX_[i] += sF * G / Rside * At;
                kAbsY_[i] += sF * G / (2.0 * Rside) * At;
            }
            if (ab == "all" && y < layer) {
                cAbsY_[i] += zP;                cAbsX_[i] += zS;
                kAbsY_[i] += sF * G / Rbot * At;
                kAbsX_[i] += sF * G / (2.0 * Rbot) * At;
            }
        }
    }

    if (scen_ == Scenario::TENSION) pullV_ = cfg_.getd("pullV", 0.05);

    fragId_.assign(p_.size(), 0);
    toolKE0_ = tool_.ke();

    std::cout << "[DEM] " << p_.size() << " particles, " << b_.size()
              << " bonds, dt = " << dt_ << " s, steps = "
              << (long)std::ceil(T_ / dt_) << "\n";
}

void DemSolver::buildPacking() {
    // Regular lattices: hexagonal (isotropic-ish, 6 neighbours) for impact /
    // cutting, square (4 neighbours, load-aligned) for the tension
    // verification. Nearest-neighbour distance is exactly 2r in both.
    double m = mat_.rho * M_PI * r_ * r_ * thk_;
    double I = 0.5 * m * r_ * r_;
    double tiny = 1e-9;

    if (packing_ == "hex") {
        double dy = r_ * std::sqrt(3.0);
        int row = 0;
        for (double y = r_; y <= H_ - r_ + tiny; y += dy, ++row) {
            double x0 = r_ + (row % 2) * r_;
            for (double x = x0; x <= W_ - r_ + tiny; x += 2.0 * r_)
                p_.push_back({{x, y}, {0, 0}, {0, 0}, 0, 0, r_, m, I, FREE});
        }
    } else {  // square
        for (double y = r_; y <= H_ - r_ + tiny; y += 2.0 * r_)
            for (double x = r_; x <= W_ - r_ + tiny; x += 2.0 * r_)
                p_.push_back({{x, y}, {0, 0}, {0, 0}, 0, 0, r_, m, I, FREE});
    }

    if (scen_ == Scenario::TENSION) {
        // grips: bottom row fixed, top row pulled at constant velocity
        double yTop = 0;
        for (const auto& q : p_) yTop = std::max(yTop, q.x.y());
        for (auto& q : p_) {
            if (q.x.y() < 2.0 * r_)          q.flag = FIXED;
            else if (q.x.y() > yTop - r_)    q.flag = PRESCRIBED;
        }
    }
}

void DemSolver::buildBonds() {
    // Bond every pair closer than 1.05*(ri+rj): exactly the lattice
    // nearest neighbours. Parallel-bond cross-section in 2D:
    //   Rb = lambda * min(ri, rj),  A = 2 Rb t,  I = (2/3) t Rb^3
    // (rectangular beam of width 2Rb and thickness t).
    double scatter = cfg_.getd("bondStrengthScatter", 0.0);
    double scBase  = cfg_.getd("bondTensile", mat_.ft);
    double tcBase  = cfg_.getd("bondCohesion", mat_.cohesion);
    std::mt19937 rng(cfg_.geti("seed", 42));
    std::uniform_real_distribution<double> U(-scatter, scatter);

    rebuildGrid();
    for (int i = 0; i < (int)p_.size(); ++i) {
        int ci = std::clamp(int((p_[i].x.x() - gmin_.x()) / cell_), 0, gx_ - 1);
        int cj = std::clamp(int((p_[i].x.y() - gmin_.y()) / cell_), 0, gy_ - 1);
        for (int dj = -1; dj <= 1; ++dj) for (int di = -1; di <= 1; ++di) {
            int cx = ci + di, cy = cj + dj;
            if (cx < 0 || cy < 0 || cx >= gx_ || cy >= gy_) continue;
            for (int j = head_[cy * gx_ + cx]; j >= 0; j = nxt_[j]) {
                if (j <= i) continue;
                double L = (p_[j].x - p_[i].x).norm();
                if (L > 1.05 * (p_[i].r + p_[j].r)) continue;
                Bond bo;
                bo.i = i; bo.j = j;
                bo.Rb = lambda_ * std::min(p_[i].r, p_[j].r);
                bo.A  = 2.0 * bo.Rb * thk_;
                bo.Ib = (2.0 / 3.0) * thk_ * bo.Rb * bo.Rb * bo.Rb;
                double kbn = mat_.E / L;              // bond modulus / length
                bo.knA = kbn * bo.A;
                bo.ksA = ksRatio_ * kbn * bo.A;
                bo.knI = kbn * bo.Ib;
                double sf = 1.0 + U(rng);
                bo.sc  = scBase * sf;
                bo.tc0 = tcBase * sf;
                b_.push_back(bo);
            }
        }
    }

    // Intact-bond count per particle: gates the absorbing dashpots (a fully
    // detached fragment is debris, not truncated continuum) and is cheaper
    // than rescanning the bond list.
    nIntact_.assign(p_.size(), 0);
    for (const auto& bo : b_) { ++nIntact_[bo.i]; ++nIntact_[bo.j]; }
}

void DemSolver::placeTool() {
    if (scen_ == Scenario::TENSION) return;
    tool_.mass = cfg_.getd("toolMass", 5.0);
    double gap = cfg_.getd("toolGap", 1e-4);
    std::string sh = cfg_.gets("toolShape", "disc");
    tool_.shape = (sh == "flat") ? Tool::Shape::FLAT : Tool::Shape::DISC;
    tool_.width  = cfg_.getd("toolWidth", 0.02);
    tool_.radius = cfg_.getd("toolRadius", 0.015);

    double yTop = 0;
    for (const auto& q : p_) yTop = std::max(yTop, q.x.y() + q.r);

    if (scen_ == Scenario::PERCUSSION) {
        tool_.motion = Tool::Motion::FREE;
        double vImp = cfg_.getd("impactSpeed", 15.0);
        double xc = cfg_.getd("toolX", 0.5 * W_);
        if (tool_.shape == Tool::Shape::FLAT) tool_.x = {xc, yTop + gap};
        else                                  tool_.x = {xc, yTop + tool_.radius + gap};
        tool_.v = {0.0, -vImp};
    } else {  // SHEAR
        tool_.motion = Tool::Motion::PRESCRIBED;
        tool_.shape  = Tool::Shape::DISC;
        double depth = cfg_.getd("cutDepth", 0.004);
        double vCut  = cfg_.getd("cutSpeed", 10.0);
        tool_.x = {-tool_.radius - gap, yTop - depth + tool_.radius};
        tool_.v = {vCut, 0.0};
    }
}

// ---------------------------------------------------------------------------
// Stable explicit time step (velocity-Verlet on a mass-spring network).
//
// A single mass-spring gives dt_crit = 2 sqrt(m/k) = 2/omega. In a dense
// assembly each particle carries several springs (contacts + bond axial +
// bond shear) and rotational DOFs, so the practical rule is
//     dt = f * sqrt( m_min / k_trans,max )   and
//     dt = f * sqrt( I_min / k_rot,max )
// with k_rot ~ kn*Ib + ksA*(L/2)^2 (bending + shear-force lever arm) and a
// safety fraction f ~ 0.15-0.2 (Cundall's recommendation for BPMs).
// ---------------------------------------------------------------------------
void DemSolver::computeStableDt() {
    // Per-particle stiffness sums: dt_crit = 2/omega_max with
    // omega_max^2 ~ (sum of spring stiffnesses attached to a particle)/m.
    // Bonds contribute (knA + ksA) translationally and knI + ksA (L/2)^2
    // rotationally. On top of the initial bonds, the crushed zone under the
    // tool develops many extra frictional contacts, so a budget of
    // `extraContacts` (default 8) linear contacts per particle is added to
    // the bound. dtFactor (default 0.2) is then a true fraction of critical.
    std::vector<double> kT(p_.size(), 0.0), kR(p_.size(), 0.0);
    for (const auto& bo : b_) {
        double L = 2.0 * r_;
        double kt = bo.knA + bo.ksA;
        double kr = bo.knI + bo.ksA * 0.25 * L * L;
        kT[bo.i] += kt; kT[bo.j] += kt;
        kR[bo.i] += kr; kR[bo.j] += kr;
    }
    double nExtra = cfg_.getd("extraContacts", 8.0);
    double dtMin = 1e30;
    for (std::size_t i = 0; i < p_.size(); ++i) {
        double kt = kT[i] + nExtra * (knC_ + ksC_);
        double kr = kR[i] + nExtra * ksC_ * p_[i].r * p_[i].r;
        dtMin = std::min(dtMin, 2.0 * std::sqrt(p_[i].m / kt));
        dtMin = std::min(dtMin, 2.0 * std::sqrt(p_[i].I / kr));
    }
    double f = cfg_.getd("dtFactor", 0.2);
    dt_ = f * dtMin;
}

void DemSolver::rebuildGrid() {
    Eigen::Vector2d lo(1e30, 1e30), hi(-1e30, -1e30);
    for (const auto& q : p_) {
        lo = lo.cwiseMin(q.x);
        hi = hi.cwiseMax(q.x);
    }
    cell_ = 2.2 * r_;
    gmin_ = lo - Eigen::Vector2d(cell_, cell_);
    Eigen::Vector2d span = hi - gmin_ + Eigen::Vector2d(cell_, cell_);
    gx_ = std::max(1, int(span.x() / cell_) + 1);
    gy_ = std::max(1, int(span.y() / cell_) + 1);
    head_.assign((std::size_t)gx_ * gy_, -1);
    nxt_.assign(p_.size(), -1);
    for (int i = 0; i < (int)p_.size(); ++i) {
        int cx = std::clamp(int((p_[i].x.x() - gmin_.x()) / cell_), 0, gx_ - 1);
        int cy = std::clamp(int((p_[i].x.y() - gmin_.y()) / cell_), 0, gy_ - 1);
        int c = cy * gx_ + cx;
        nxt_[i] = head_[c];
        head_[c] = i;
    }
}

// ===========================================================================
// Time stepping
// ===========================================================================

void DemSolver::step() {
    for (auto& q : p_) { q.f.setZero(); q.tq = 0; }
    tool_.resetForce();

    bondForces();
    rebuildGrid();
    contactForces();
    wallAndToolForces();

    if (scen_ == Scenario::TENSION) {
        gripF_.setZero();
        for (const auto& q : p_)
            if (q.flag == PRESCRIBED) gripF_ += q.f;   // reaction on the grip
        double sigma = std::abs(gripF_.y()) / (W_ * thk_);
        sigmaPeak_ = std::max(sigmaPeak_, sigma);
    }
    peakF_ = std::max(peakF_, tool_.F.norm());
    // Energy delivered by the tool (thrust/cutting work). Tool-side
    // bookkeeping deliberately includes the interface friction and dashpot
    // dissipation: that is part of the drilling energy input.
    if (scen_ != Scenario::TENSION) work_ += -tool_.F.dot(tool_.v) * dt_;

    integrate();
    t_ += dt_;
}

// ---------------------------------------------------------------------------
// Parallel bonds (Potyondy & Cundall's BPM, 2D).
//
// Each bond is a short beam glued between two particles. It carries an axial
// force Fn (TENSION POSITIVE), a shear force Fs and a bending moment Mb,
// updated INCREMENTALLY from the relative velocities:
//     dFn = + kn A  (v_rel . n) dt          (stretch -> tension grows)
//     dFs = - ks A  v_t,rel     dt          (v_t,rel includes rotation terms)
//     dMb = - kn I  (w_j - w_i) dt          (bending resisted via kn, PFC2D)
// with v_t,rel = (v_j - v_i).t - (w_i + w_j) L/2 the tangential slip of the
// contact mid-point.
//
// Beam theory gives the extreme-fibre stresses at the bond periphery:
//     sigma_max = Fn/A + |Mb| Rb / I        (max tensile stress)
//     tau_max   = |Fs| / A
// Breakage — the tensile AND shear criteria required by the spec:
//     sigma_max > sigma_c            -> tensile bond failure
//     tau_max   > c_b + tan(phi_b) * max(0, -Fn/A)   -> shear failure
// (shear strength = cohesion + friction * compressive normal stress,
// Mohr-Coulomb-like). A broken bond simply stops transmitting bond loads;
// the pair automatically falls back to the frictional contact below, so
// cracks and fragments EMERGE from the topology of broken bonds.
// ---------------------------------------------------------------------------
void DemSolver::bondForces() {
    for (auto& bo : b_) {
        if (bo.broken) continue;
        Part& pi = p_[bo.i];
        Part& pj = p_[bo.j];

        Eigen::Vector2d d = pj.x - pi.x;
        double L = d.norm();
        Eigen::Vector2d n = d / L;
        Eigen::Vector2d tdir(-n.y(), n.x());

        Eigen::Vector2d vrel = pj.v - pi.v;
        double vn = vrel.dot(n);
        double vt = vrel.dot(tdir) - (pi.w + pj.w) * 0.5 * L;
        double wrel = pj.w - pi.w;

        bo.Fn += bo.knA * vn * dt_;
        bo.Fs += -bo.ksA * vt * dt_;
        bo.Mb += -bo.knI * wrel * dt_;

        double sigMax = bo.Fn / bo.A + std::abs(bo.Mb) * bo.Rb / bo.Ib;
        double tau    = std::abs(bo.Fs) / bo.A;
        double tauC   = bo.tc0 + tanPhiB_ * std::max(0.0, -bo.Fn / bo.A);
        if (sigMax > bo.sc || tau > tauC) {
            bo.broken = true;
            bo.tBreak = t_;
            ++nBroken_;
            --nIntact_[bo.i];
            --nIntact_[bo.j];
            continue;
        }

        // Apply. Tension (Fn > 0) pulls the particles together.
        Eigen::Vector2d Fj = -bo.Fn * n + bo.Fs * tdir;
        pj.f += Fj;
        pi.f -= Fj;
        // Tangential force couple (applied at the mid-point) + bending moment
        pi.tq += -0.5 * L * bo.Fs - bo.Mb;
        pj.tq += -0.5 * L * bo.Fs + bo.Mb;
    }
}

// ---------------------------------------------------------------------------
// Grain-to-grain contact (acts in PARALLEL with bonds, and alone once a bond
// has broken): linear normal spring with viscous damping (compression only),
// tangential spring with history capped by Coulomb friction mu*Fn.
// ---------------------------------------------------------------------------
void DemSolver::contactForces() {
    tangNew_.clear();
    for (int i = 0; i < (int)p_.size(); ++i) {
        int ci = std::clamp(int((p_[i].x.x() - gmin_.x()) / cell_), 0, gx_ - 1);
        int cj = std::clamp(int((p_[i].x.y() - gmin_.y()) / cell_), 0, gy_ - 1);
        for (int dj = -1; dj <= 1; ++dj) for (int di = -1; di <= 1; ++di) {
            int cx = ci + di, cy = cj + dj;
            if (cx < 0 || cy < 0 || cx >= gx_ || cy >= gy_) continue;
            for (int j = head_[cy * gx_ + cx]; j >= 0; j = nxt_[j]) {
                if (j <= i) continue;
                Part& pi = p_[i];
                Part& pj = p_[j];
                Eigen::Vector2d d = pj.x - pi.x;
                double dist2 = d.squaredNorm();
                double rsum = pi.r + pj.r;
                if (dist2 >= rsum * rsum) continue;
                double dist = std::sqrt(dist2);
                double delta = rsum - dist;
                Eigen::Vector2d n = d / dist;
                Eigen::Vector2d tdir(-n.y(), n.x());

                double li = pi.r - 0.5 * delta;   // lever arms to contact point
                double lj = pj.r - 0.5 * delta;
                Eigen::Vector2d vrel = pj.v - pi.v;
                double vn = vrel.dot(n);
                double vt = vrel.dot(tdir) - (pi.w * li + pj.w * lj);

                double meff = pi.m * pj.m / (pi.m + pj.m);
                double cn = 2.0 * xiC_ * std::sqrt(knC_ * meff);
                double fn = knC_ * delta - cn * vn;
                if (fn < 0) fn = 0;

                uint64_t key = pairKey(i, j);
                auto it = tang_.find(key);
                double fs = (it == tang_.end()) ? 0.0 : it->second;
                fs += -ksC_ * vt * dt_;
                double fmax = mu_ * fn;
                fs = std::clamp(fs, -fmax, fmax);
                tangNew_[key] = fs;

                Eigen::Vector2d Fj = fn * n + fs * tdir;
                pj.f += Fj;
                pi.f -= Fj;
                pi.tq += -li * fs;
                pj.tq += -lj * fs;
            }
        }
    }
    tang_.swap(tangNew_);   // springs of separated contacts are dropped
}

void DemSolver::toolParticleContact(Part& q) {
    Eigen::Vector2d n;      // outward from the tool, i.e. push direction
    double delta;
    Eigen::Vector2d vTool = tool_.v;

    if (tool_.shape == Tool::Shape::DISC) {
        Eigen::Vector2d d = q.x - tool_.x;
        double dist = d.norm();
        double rsum = q.r + tool_.radius;
        if (dist >= rsum || dist < 1e-14) return;
        n = d / dist;
        delta = rsum - dist;
    } else {  // FLAT bottom face pushing down
        if (std::abs(q.x.x() - tool_.x.x()) > 0.5 * tool_.width + q.r) return;
        delta = q.x.y() + q.r - tool_.x.y();
        if (delta <= 0) return;
        n = {0.0, -1.0};
    }
    Eigen::Vector2d tdir(-n.y(), n.x());
    Eigen::Vector2d vrel = q.v - vTool;
    double cn = 2.0 * xiC_ * std::sqrt(knC_ * q.m);
    double fn = knC_ * delta - cn * vrel.dot(n);
    if (fn < 0) fn = 0;
    double vt = vrel.dot(tdir) - q.w * (q.r - 0.5 * delta);
    double fs = -mu_ * fn * std::tanh(vt / 1e-3);   // regularized Coulomb

    Eigen::Vector2d F = fn * n + fs * tdir;
    q.f += F;
    // Torque: tangential force fs*t applied at the contact point
    // r_c = -(r - delta/2) n  ->  tau = cross(r_c, fs t) = -(r - delta/2) fs
    q.tq += -(q.r - 0.5 * delta) * fs;

    tool_.F -= F;
}

void DemSolver::wallAndToolForces() {
    for (auto& q : p_) {
        if (bottomWall_) {
            double delta = q.r - q.x.y();
            if (delta > 0) {
                double cn = 2.0 * xiC_ * std::sqrt(knC_ * q.m);
                double fn = knC_ * delta - cn * q.v.y();
                if (fn < 0) fn = 0;
                double vt = q.v.x() + q.w * (q.r - 0.5 * delta);  // slip of bottom point
                double fs = -mu_ * fn * std::tanh(vt / 1e-3);
                q.f += Eigen::Vector2d(fs, fn);
                // r_c = -(r - delta/2) j_hat, F_t = (fs, 0)
                //   -> tau = cross(r_c, F_t) = +(r - delta/2) fs
                q.tq += (q.r - 0.5 * delta) * fs;
            }
        }
        if (sideWalls_) {
            double dL = q.r - q.x.x();
            if (dL > 0) q.f.x() += knC_ * dL;
            double dR = q.r - (W_ - q.x.x());
            if (dR > 0) q.f.x() -= knC_ * dR;
        }
        if (scen_ != Scenario::TENSION) toolParticleContact(q);
    }
}

// Leapfrog kick-drift, with Cundall local non-viscous damping on free
// particles:  F <- F - alpha |F| sign(v) (componentwise). Used heavily
// (0.7) for the quasi-static tension test, lightly (~0.02) for dynamics.
void DemSolver::integrate() {
    for (auto& q : p_) {
        if (q.flag == FIXED) { q.v.setZero(); q.w = 0; continue; }
        if (q.flag == PRESCRIBED) {
            q.v = {0.0, pullV_};
            q.w = 0;
            q.x += dt_ * q.v;
            continue;
        }
        std::size_t ib = (std::size_t)(&q - p_.data());
        if (nIntact_[ib] > 0 && (kAbsX_[ib] > 0 || kAbsY_[ib] > 0)) {
            // boundary spring of the viscous-spring quiet boundary
            q.f.x() -= kAbsX_[ib] * (q.x.x() - xAnchor_[ib].x());
            q.f.y() -= kAbsY_[ib] * (q.x.y() - xAnchor_[ib].y());
        }
        if (damping_ > 0) {
            for (int k = 0; k < 2; ++k)
                q.f(k) -= damping_ * std::abs(q.f(k)) * (q.v(k) > 0 ? 1.0 : (q.v(k) < 0 ? -1.0 : 0.0));
            q.tq -= damping_ * std::abs(q.tq) * (q.w > 0 ? 1.0 : (q.w < 0 ? -1.0 : 0.0));
        }
        q.v += (dt_ / q.m) * q.f;
        q.w += (dt_ / q.I) * q.tq;
        // Lysmer dashpot (quiet boundaries), implicit in v — see init().
        // Only while the particle is still bonded to the block: the quiet
        // boundary stands in for the truncated elastic continuum, and must
        // not drag on detached debris flying through the layer.
        if (nIntact_[ib] > 0) {
            if (cAbsX_[ib] > 0) q.v.x() /= 1.0 + dt_ * cAbsX_[ib] / q.m;
            if (cAbsY_[ib] > 0) q.v.y() /= 1.0 + dt_ * cAbsY_[ib] / q.m;
        }
        q.x += dt_ * q.v;
    }
    if (scen_ != Scenario::TENSION) tool_.integrate(dt_);
}

// ===========================================================================
// Fragments: connected components of the INTACT bond network (union-find).
// Fragment 0 = the largest (the remaining specimen); everything else counts
// as detached / comminuted material for the specific-energy estimate.
// ===========================================================================

int DemSolver::uf(std::vector<int>& parent, int a) const {
    while (parent[a] != a) {
        parent[a] = parent[parent[a]];
        a = parent[a];
    }
    return a;
}

void DemSolver::computeFragments() {
    std::vector<int> parent(p_.size());
    for (int i = 0; i < (int)p_.size(); ++i) parent[i] = i;
    for (const auto& bo : b_)
        if (!bo.broken) {
            int ra = uf(parent, bo.i);
            int rb = uf(parent, bo.j);
            if (ra != rb) parent[ra] = rb;
        }

    std::map<int, int> count;
    for (int i = 0; i < (int)p_.size(); ++i) ++count[uf(parent, i)];

    std::vector<std::pair<int, int>> order;   // (size, root), largest first
    for (auto& [root, c] : count) order.push_back({c, root});
    std::sort(order.rbegin(), order.rend());

    std::map<int, int> remap;
    for (int k = 0; k < (int)order.size(); ++k) remap[order[k].second] = k;

    fragId_.resize(p_.size());
    detachedVol_ = 0;
    for (int i = 0; i < (int)p_.size(); ++i) {
        fragId_[i] = remap[uf(parent, i)];
        if (fragId_[i] != 0) detachedVol_ += M_PI * p_[i].r * p_[i].r * thk_;
    }
    nFrag_ = (int)order.size();
}

// ===========================================================================
// Output
// ===========================================================================

void DemSolver::writeFrame(int frame) {
    computeFragments();

    std::vector<Eigen::Vector2d> pts(p_.size());
    std::vector<double> rad(p_.size()), frag(p_.size()), spd(p_.size());
    std::vector<Eigen::Vector2d> vel(p_.size());
    for (std::size_t i = 0; i < p_.size(); ++i) {
        pts[i] = p_[i].x;
        rad[i] = p_[i].r;
        frag[i] = fragId_[i];
        spd[i] = p_[i].v.norm();
        vel[i] = p_[i].v;
    }

    char name[64];
    std::snprintf(name, sizeof(name), "/dem_particles_%04d.vtu", frame);
    vtk::writeParticles(out_ + name, pts,
                        {{"radius", &rad}, {"fragment", &frag}, {"speed", &spd}},
                        {{"velocity", &vel}});

    // bond network: state 0 = intact, 1 = broken (the crack network)
    std::vector<std::array<int, 2>> lines;
    std::vector<double> state, tb;
    for (const auto& bo : b_) {
        lines.push_back({bo.i, bo.j});
        state.push_back(bo.broken ? 1.0 : 0.0);
        tb.push_back(bo.tBreak);
    }
    std::snprintf(name, sizeof(name), "/dem_bonds_%04d.vtu", frame);
    vtk::writeLines(out_ + name, pts, lines, {{"state", &state}, {"tBreak", &tb}});

    // frame index -> (time, tool pose) for make_gif.py
    std::ofstream fm(out_ + "/frames.csv",
                     frame == 0 ? std::ios::trunc : std::ios::app);
    if (frame == 0) fm << "frame,t,toolX,toolY\n";
    fm << frame << "," << t_ << "," << tool_.x.x() << "," << tool_.x.y() << "\n";
}

void DemSolver::historyHeader(std::ostream& os) const {
    if (scen_ == Scenario::TENSION) { os << "t,gripFy,sigma,sigmaPeak,nBroken\n"; return; }
    os << "t,toolFx,toolFy,toolX,toolY,toolVx,toolVy,work,toolKE,"
          "nBroken,nFrag,detachedVol,specificEnergy\n";
}

void DemSolver::historyRow(std::ostream& os) const {
    if (scen_ == Scenario::TENSION) {
        os << t_ << "," << gripF_.y() << "," << std::abs(gripF_.y()) / (W_ * thk_)
           << "," << sigmaPeak_ << "," << nBroken_ << "\n";
        return;
    }
    double Es = detachedVol_ > 0 ? work_ / detachedVol_ : 0.0;
    os << t_ << "," << tool_.F.x() << "," << tool_.F.y() << ","
       << tool_.x.x() << "," << tool_.x.y() << ","
       << tool_.v.x() << "," << tool_.v.y() << ","
       << work_ << "," << tool_.ke() << ","
       << nBroken_ << "," << nFrag_ << "," << detachedVol_ << "," << Es << "\n";
}

void DemSolver::finalize() {
    computeFragments();

    // fragment-size distribution
    std::map<int, double> mass;
    std::map<int, int> np;
    for (int i = 0; i < (int)p_.size(); ++i) {
        mass[fragId_[i]] += p_[i].m;
        np[fragId_[i]] += 1;
    }
    std::ofstream ff(out_ + "/dem_fragments.csv");
    ff << "fragment,nParticles,mass\n";
    for (auto& [id, mm] : mass) ff << id << "," << np[id] << "," << mm << "\n";

    // convenience CSV snapshots for quick plotting
    std::ofstream fp(out_ + "/dem_final_particles.csv");
    fp << "x,y,r,fragment\n";
    for (int i = 0; i < (int)p_.size(); ++i)
        fp << p_[i].x.x() << "," << p_[i].x.y() << "," << p_[i].r << "," << fragId_[i] << "\n";
    std::ofstream fb(out_ + "/dem_final_bonds.csv");
    fb << "x1,y1,x2,y2,broken,tBreak\n";
    for (const auto& bo : b_)
        fb << p_[bo.i].x.x() << "," << p_[bo.i].x.y() << ","
           << p_[bo.j].x.x() << "," << p_[bo.j].x.y() << ","
           << (bo.broken ? 1 : 0) << "," << bo.tBreak << "\n";

    std::cout << "\n[DEM] ---- summary ----\n";
    if (scen_ == Scenario::TENSION) {
        double expected = cfg_.getd("bondTensile", mat_.ft) * lambda_;
        double err = 100.0 * (sigmaPeak_ - expected) / expected;
        std::cout << "[DEM] tension verification (square lattice, load-aligned):\n"
                  << "[DEM]   peak macro stress = " << sigmaPeak_ / 1e6 << " MPa, "
                  << "expected ~ lambda * sigma_c = " << expected / 1e6 << " MPa, "
                  << "error = " << err << " %  "
                  << (std::abs(err) < 15.0 ? "[PASS]" : "[FAIL]") << "\n"
                  << "[DEM]   broken bonds = " << nBroken_ << "\n";
        return;
    }
    double Es = detachedVol_ > 0 ? work_ / detachedVol_ : 0.0;
    std::cout << "[DEM] peak tool force   : " << peakF_ << " N/m\n"
              << "[DEM] tool work output  : " << work_ << " J/m";
    if (tool_.motion == Tool::Motion::FREE)
        std::cout << "  (tool KE loss: " << toolKE0_ - tool_.ke() << " J/m)";
    std::cout << "\n"
              << "[DEM] broken bonds      : " << nBroken_ << " / " << b_.size() << "\n"
              << "[DEM] fragments         : " << nFrag_ << " (detached vol "
              << detachedVol_ << " m^3/m)\n"
              << "[DEM] specific energy   : " << Es << " J/m^3\n";
}

} // namespace rockim
