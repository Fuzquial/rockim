// ---------------------------------------------------------------------------
// Dem3dSolver — 3D bonded-particle model. See the header for the model
// overview; the comments here focus on the derivations (bond kinematics,
// failure criteria, torque bookkeeping, stable time step).
// ---------------------------------------------------------------------------
#include "rockim/Dem3dSolver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <random>

#include "rockim/VtkWriter.hpp"

#include <chrono>
#include <cstdlib>

namespace rockim {

// crude phase profiler, enabled with ROCKIM_PROF=1 (report at destruction)
namespace {
struct Prof {
    double tBond = 0, tGrid = 0, tCont = 0, tTool = 0, tInt = 0;
    long n = 0;
    bool on = std::getenv("ROCKIM_PROF") != nullptr;
    ~Prof() {
        if (!on || n == 0) return;
        std::fprintf(stderr,
                     "[prof] per step (ms): bond %.2f grid %.2f contact %.2f "
                     "tool %.2f integ %.2f  (%ld steps)\n",
                     1e3 * tBond / n, 1e3 * tGrid / n, 1e3 * tCont / n,
                     1e3 * tTool / n, 1e3 * tInt / n, n);
    }
} gProf;
double now() {
    return std::chrono::duration<double>(
               std::chrono::steady_clock::now().time_since_epoch()).count();
}
} // namespace

Dem3dSolver::Dem3dSolver(const Config& cfg, std::string outDir)
    : cfg_(cfg), out_(std::move(outDir)) {}

void Dem3dSolver::init() {
    mat_ = Material::from(cfg_);

    std::string sc = cfg_.gets("scenario", "percussion");
    if      (sc == "percussion") scen_ = Scenario::PERCUSSION;
    else if (sc == "shear")      scen_ = Scenario::SHEAR;
    else if (sc == "tension")    scen_ = Scenario::TENSION;
    else throw std::runtime_error("dem3d scenario must be percussion | shear | tension");

    W_ = cfg_.getd("W", 0.1);
    D_ = cfg_.getd("D", 0.1);
    H_ = cfg_.getd("H", 0.08);
    r_ = cfg_.getd("particleRadius", 1.5e-3);
    packing_ = cfg_.gets("packing", scen_ == Scenario::TENSION ? "cubic" : "hcp");
    T_ = cfg_.getd("T", 2e-4);

    // Contact stiffness heuristic (same spirit as 2D's kn = E*t): a linear
    // spring with k ~ E * (grain size) so grain-scale contact compliance is
    // consistent with the elastic modulus. Quantitative macro elasticity
    // requires calibration either way.
    knC_ = cfg_.getd("knFactor", 1.0) * mat_.E * 2.0 * r_;
    ksRatio_ = cfg_.getd("ksRatio", ksRatio_);
    ksC_ = ksRatio_ * knC_;
    mu_  = cfg_.getd("contactMu", mu_);
    xiC_ = cfg_.getd("contactXi", xiC_);

    tanPhiB_ = std::tan(cfg_.getd("bondFrictionDeg", mat_.phiDeg) * M_PI / 180.0);
    lambda_  = cfg_.getd("bondRadiusFactor", 1.0);
    damping_ = cfg_.getd("dampingLocal", scen_ == Scenario::TENSION ? 0.7 : 0.02);
    bottomWall_ = cfg_.getb("bottomWall", scen_ != Scenario::TENSION);

    buildPacking();
    buildBonds();
    placeTool();
    computeStableDt();
    setupAbsorbing();

    if (scen_ == Scenario::TENSION) pullV_ = cfg_.getd("pullV", 0.05);
    pullRamp_ = cfg_.getd("pullRamp", 0.0);
    gripFree_ = cfg_.getb("gripLateralFree", false);

    fragId_.assign(p_.size(), 0);
    toolKE0_ = tool_.ke();

    std::cout << "[DEM3D] " << p_.size() << " particles, " << b_.size()
              << " bonds, dt = " << dt_ << " s, steps = "
              << (long)std::ceil(T_ / dt_) << "\n";
}

// ---------------------------------------------------------------------------
// Packings. HCP (ABAB stacking of hexagonal layers): every sphere touches 12
// neighbours at exactly 2r — 6 in its layer, 3 below, 3 above. Layer spacing
// dz = 2r sqrt(2/3); B-layers are shifted by (r, r/sqrt(3)) to sit in the
// hollows of A-layers. Simple cubic (6 neighbours at 2r, load-aligned) is
// used for the tension verification: only the z-bonds carry load, so the
// macroscopic peak stress must equal
//     sigma_peak = (A_bond / A_cell) * sigma_c = (pi/4) lambda^2 sigma_c,
// a sharp parameter-free check of the 3D bond implementation.
// ---------------------------------------------------------------------------
void Dem3dSolver::buildPacking() {
    double m = mat_.rho * (4.0 / 3.0) * M_PI * r_ * r_ * r_;
    double I = 0.4 * m * r_ * r_;             // solid sphere: (2/5) m r^2
    double tiny = 1e-9;

    if (packing_ == "hcp") {
        double dy = r_ * std::sqrt(3.0);
        double dz = 2.0 * r_ * std::sqrt(2.0 / 3.0);
        int layer = 0;
        for (double z = r_; z <= H_ - r_ + tiny; z += dz, ++layer) {
            double ox = (layer % 2) * r_;
            double oy = (layer % 2) * r_ / std::sqrt(3.0);
            int row = 0;
            for (double y = r_ + oy; y <= D_ - r_ + tiny; y += dy, ++row) {
                double x0 = r_ + ox + (row % 2) * r_;
                for (double x = x0; x <= W_ - r_ + tiny; x += 2.0 * r_)
                    p_.push_back({{x, y, z}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
                                  r_, m, I, FREE});
            }
        }
    } else {  // simple cubic
        for (double z = r_; z <= H_ - r_ + tiny; z += 2.0 * r_)
            for (double y = r_; y <= D_ - r_ + tiny; y += 2.0 * r_)
                for (double x = r_; x <= W_ - r_ + tiny; x += 2.0 * r_)
                    p_.push_back({{x, y, z}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
                                  r_, m, I, FREE});
    }

    // Sort particles along a coarse spatial key: neighbours in space become
    // neighbours in memory, which roughly halves the cache-miss cost of the
    // cell-list sweeps (indices are assigned after the sort, so bonds and
    // flags are unaffected).
    double ck = 2.2 * r_;
    std::sort(p_.begin(), p_.end(), [ck](const Part& a, const Part& b) {
        auto key = [ck](const Part& q) {
            return std::array<long, 3>{lround(q.x.z() / ck),
                                       lround(q.x.y() / ck),
                                       lround(q.x.x() / ck)};
        };
        return key(a) < key(b);
    });

    if (scen_ == Scenario::TENSION) {
        double zTop = 0;
        for (const auto& q : p_) zTop = std::max(zTop, q.x.z());
        for (auto& q : p_) {
            if (q.x.z() < 2.0 * r_)       q.flag = FIXED;
            else if (q.x.z() > zTop - r_) q.flag = PRESCRIBED;
        }
    }
}

// ---------------------------------------------------------------------------
// Parallel bonds. Cross-section of a 3D parallel bond of radius Rb:
//   A = pi Rb^2,  I = pi Rb^4 / 4 (bending),  J = pi Rb^4 / 2 (torsion).
// Per-area bond modulus kbn = E / L0 (PFC convention), shear = ksRatio * kbn.
// ---------------------------------------------------------------------------
void Dem3dSolver::buildBonds() {
    double scatter = cfg_.getd("bondStrengthScatter", 0.0);
    double scBase  = cfg_.getd("bondTensile", mat_.ft);
    double tcBase  = cfg_.getd("bondCohesion", mat_.cohesion);
    std::mt19937 rng(cfg_.geti("seed", 12345));
    std::uniform_real_distribution<double> U(-scatter, scatter);

    rebuildGrid();
    for (int i = 0; i < (int)p_.size(); ++i) {
        int ci = std::clamp(int((p_[i].x.x() - gmin_.x()) / cell_), 0, gx_ - 1);
        int cj = std::clamp(int((p_[i].x.y() - gmin_.y()) / cell_), 0, gy_ - 1);
        int ck = std::clamp(int((p_[i].x.z() - gmin_.z()) / cell_), 0, gz_ - 1);
        for (int dk = -1; dk <= 1; ++dk)
        for (int dj = -1; dj <= 1; ++dj)
        for (int di = -1; di <= 1; ++di) {
            int cx = ci + di, cy = cj + dj, cz = ck + dk;
            if (cx < 0 || cy < 0 || cz < 0 || cx >= gx_ || cy >= gy_ || cz >= gz_)
                continue;
            for (int j = head_[(cz * gy_ + cy) * gx_ + cx]; j >= 0; j = nxt_[j]) {
                if (j <= i) continue;
                double L = (p_[j].x - p_[i].x).norm();
                if (L > 1.05 * (p_[i].r + p_[j].r)) continue;
                Bond bo;
                bo.i = i; bo.j = j;
                bo.Rb = lambda_ * std::min(p_[i].r, p_[j].r);
                bo.A  = M_PI * bo.Rb * bo.Rb;
                bo.Ib = 0.25 * M_PI * std::pow(bo.Rb, 4);
                bo.Jb = 2.0 * bo.Ib;
                double kbn = mat_.E / L;
                bo.knA = kbn * bo.A;
                bo.ksA = ksRatio_ * kbn * bo.A;
                bo.knI = kbn * bo.Ib;
                bo.ksJ = ksRatio_ * kbn * bo.Jb;
                double sf = 1.0 + U(rng);
                bo.sc  = scBase * sf;
                bo.tc0 = tcBase * sf;
                b_.push_back(bo);
            }
        }
    }

    nIntact_.assign(p_.size(), 0);
    for (const auto& bo : b_) { ++nIntact_[bo.i]; ++nIntact_[bo.j]; }
}

void Dem3dSolver::placeTool() {
    if (scen_ == Scenario::TENSION) return;
    tool_.mass   = cfg_.getd("toolMass", 0.5);
    tool_.radius = cfg_.getd("toolRadius", 0.015);
    double gap = cfg_.getd("toolGap", 1e-4);
    // toolShape = sphere | flat ('disc' accepted as the 2D synonym); shear
    // forces the sphere, exactly as 2D shear forces the disc
    std::string sh = cfg_.gets("toolShape", "sphere");
    if (sh != "sphere" && sh != "disc" && sh != "flat")
        throw std::runtime_error("toolShape must be sphere | flat (3D)");
    tool_.flat = sh == "flat" && scen_ == Scenario::PERCUSSION;

    double zTop = 0;
    for (const auto& q : p_) zTop = std::max(zTop, q.x.z() + q.r);

    if (scen_ == Scenario::PERCUSSION) {
        tool_.free = true;
        double vImp = cfg_.getd("impactSpeed", 8.0);
        double zTip = tool_.flat ? zTop + gap : zTop + tool_.radius + gap;
        tool_.x = {cfg_.getd("toolX", 0.5 * W_), cfg_.getd("toolY", 0.5 * D_),
                   zTip};
        tool_.v = {0.0, 0.0, -vImp};
    } else {  // SHEAR: spherical cutter dragged along +x at a depth of cut
        tool_.free = false;
        double depth = cfg_.getd("cutDepth", 0.004);
        double vCut  = cfg_.getd("cutSpeed", 10.0);
        tool_.x = {cfg_.getd("toolX", -tool_.radius - gap), 0.5 * D_,
                   zTop - depth + tool_.radius};
        tool_.v = {vCut, 0.0, 0.0};
    }
}

// ---------------------------------------------------------------------------
// Quiet boundaries: Lysmer dashpots (rho c_p normal to the face, rho c_s in
// the two tangential directions) + Deeks-Randolph springs (k_n = G/R,
// k_t = G/2R) lumped over each boundary particle's tributary face area,
// taken as the close-packed per-particle area 2*sqrt(3) r^2. Applied on the
// four lateral faces ('sides') and additionally on the bottom ('all', which
// then replaces the rigid bottom wall). Gated on bondedness in integrate():
// the boundary stands in for the truncated continuum, not for flying debris.
// ---------------------------------------------------------------------------
void Dem3dSolver::setupAbsorbing() {
    cAbs_.assign(p_.size(), Eigen::Vector3d::Zero());
    kAbs_.assign(p_.size(), Eigen::Vector3d::Zero());
    xAnchor_.resize(p_.size());
    for (std::size_t i = 0; i < p_.size(); ++i) xAnchor_[i] = p_[i].x;

    std::string ab = cfg_.gets("absorbing", "none");
    if (scen_ == Scenario::TENSION || ab == "none") return;
    if (ab != "sides" && ab != "all")
        throw std::runtime_error("absorbing must be none | sides | all");
    if (ab == "all") bottomWall_ = false;

    double layer = cfg_.getd("absorbLayer", 2.2) * r_;
    double G  = mat_.E / (2.0 * (1.0 + mat_.nu));
    double sF = cfg_.getd("absorbSpringFactor", 1.0);
    double Rx = cfg_.getd("absorbSpringR", 0.5 * W_);
    double Ry = cfg_.getd("absorbSpringR", 0.5 * D_);
    double Rz = cfg_.getd("absorbSpringR", H_);

    for (std::size_t i = 0; i < p_.size(); ++i) {
        double At = 2.0 * std::sqrt(3.0) * p_[i].r * p_[i].r;
        double zP = mat_.rho * mat_.cP() * At;
        double zS = mat_.rho * mat_.cS() * At;
        const auto& x = p_[i].x;

        auto face = [&](int nAxis, double R) {   // face with outward normal nAxis
            for (int a = 0; a < 3; ++a) {
                if (a == nAxis) {
                    cAbs_[i](a) += zP;
                    kAbs_[i](a) += sF * G / R * At;
                } else {
                    cAbs_[i](a) += zS;
                    kAbs_[i](a) += sF * G / (2.0 * R) * At;
                }
            }
        };
        if (x.x() < layer || x.x() > W_ - layer) face(0, Rx);
        if (x.y() < layer || x.y() > D_ - layer) face(1, Ry);
        if (ab == "all" && x.z() < layer)        face(2, Rz);
    }
}

// ---------------------------------------------------------------------------
// Stable time step: same philosophy as 2D — per-particle sums of ALL attached
// stiffnesses. Translational: bonds (knA + ksA) + a budget of extra frictional
// contacts. Rotational: bond bending knI + torsion ksJ + the shear-force
// lever arm ksA (L/2)^2, plus the contact lever ksC r^2 budget.
//     dt = dtFactor * min_i ( 2 sqrt(m_i/K_t,i), 2 sqrt(I_i/K_r,i) )
// ---------------------------------------------------------------------------
void Dem3dSolver::computeStableDt() {
    std::vector<double> kT(p_.size(), 0.0), kR(p_.size(), 0.0);
    for (const auto& bo : b_) {
        double L = 2.0 * r_;
        double kt = bo.knA + bo.ksA;
        double kr = bo.knI + bo.ksJ + bo.ksA * 0.25 * L * L;
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
    dt_ = cfg_.getd("dtFactor", 0.2) * dtMin;
}

void Dem3dSolver::rebuildGrid() {
    Eigen::Vector3d lo(1e30, 1e30, 1e30), hi(-1e30, -1e30, -1e30);
    for (const auto& q : p_) {
        lo = lo.cwiseMin(q.x);
        hi = hi.cwiseMax(q.x);
    }
    cell_ = 2.2 * r_;
    gmin_ = lo - Eigen::Vector3d::Constant(cell_);
    Eigen::Vector3d span = hi - gmin_ + Eigen::Vector3d::Constant(cell_);
    gx_ = std::max(1, int(span.x() / cell_) + 1);
    gy_ = std::max(1, int(span.y() / cell_) + 1);
    gz_ = std::max(1, int(span.z() / cell_) + 1);
    head_.assign((std::size_t)gx_ * gy_ * gz_, -1);
    nxt_.assign(p_.size(), -1);
    for (int i = 0; i < (int)p_.size(); ++i) {
        int cx = std::clamp(int((p_[i].x.x() - gmin_.x()) / cell_), 0, gx_ - 1);
        int cy = std::clamp(int((p_[i].x.y() - gmin_.y()) / cell_), 0, gy_ - 1);
        int cz = std::clamp(int((p_[i].x.z() - gmin_.z()) / cell_), 0, gz_ - 1);
        int c = (cz * gy_ + cy) * gx_ + cx;
        nxt_[i] = head_[c];
        head_[c] = i;
    }
}

// ===========================================================================
// Time stepping
// ===========================================================================

void Dem3dSolver::step() {
    for (auto& q : p_) { q.f.setZero(); q.tq.setZero(); }
    tool_.F.setZero();

    if (gProf.on) {
        double t0 = now(); bondForces();
        double t1 = now(); rebuildGrid();
        double t2 = now(); contactForces();
        double t3 = now(); toolAndWallForces();
        double t4 = now();
        gProf.tBond += t1 - t0; gProf.tGrid += t2 - t1;
        gProf.tCont += t3 - t2; gProf.tTool += t4 - t3;
        ++gProf.n;
    } else {
        bondForces();
        rebuildGrid();
        contactForces();
        toolAndWallForces();
    }

    if (scen_ == Scenario::TENSION) {
        gripF_.setZero();
        for (const auto& q : p_)
            if (q.flag == PRESCRIBED) gripF_ += q.f;
        double sigma = std::abs(gripF_.z()) / (W_ * D_);
        sigmaPeak_ = std::max(sigmaPeak_, sigma);
    }
    peakF_ = std::max(peakF_, tool_.F.norm());
    if (scen_ != Scenario::TENSION) work_ += -tool_.F.dot(tool_.v) * dt_;

    integrate();
    t_ += dt_;
}

// ---------------------------------------------------------------------------
// Parallel-bond kinematics (contact point c at the bond midpoint):
//   v_rel = v(c on j) - v(c on i) = (vj - vi) - (L/2)(wi + wj) x n
//   normal rate  vn = v_rel . n          -> dFn = +knA vn dt (tension > 0)
//   shear rate   vt = v_rel - vn n       -> dFs = -ksA vt dt
//   twist rate   (wj - wi) . n           -> dMt = -ksJ  twist dt
//   bending rate (wj - wi) perp part     -> dMb = -knI  bend  dt
// Stored Fs and Mb live in the plane perpendicular to n; they are kept there
// by projection each step (exact for the small per-step bond rotations).
// Failure (beam theory on the bond periphery):
//   sigma_max = Fn/A + |Mb| Rb / I   > sigma_c            (tension)
//   tau_max   = |Fs|/A + |Mt| Rb / J > c + tan(phi) max(0, -Fn/A)  (shear)
// Application (energy-consistent, mirrors the verified 2D scheme):
//   F_j = -Fn n + Fs,  F_i = -F_j
//   both particles receive the shear couple  -(L/2) n x Fs,
//   and the moment M = Mb + Mt n acts as  tq_i -= M, tq_j += M.
// ---------------------------------------------------------------------------
void Dem3dSolver::bondForces() {
    for (auto& bo : b_) {
        if (bo.broken) continue;
        Part& pi = p_[bo.i];
        Part& pj = p_[bo.j];

        Eigen::Vector3d d = pj.x - pi.x;
        double L = d.norm();
        Eigen::Vector3d n = d / L;

        // keep stored perpendicular quantities perpendicular to the new axis
        bo.Fs -= bo.Fs.dot(n) * n;
        bo.Mb -= bo.Mb.dot(n) * n;

        Eigen::Vector3d vrel = (pj.v - pi.v) - 0.5 * L * (pi.w + pj.w).cross(n);
        double vn = vrel.dot(n);
        Eigen::Vector3d vt = vrel - vn * n;
        Eigen::Vector3d wrel = pj.w - pi.w;
        double twist = wrel.dot(n);
        Eigen::Vector3d bend = wrel - twist * n;

        bo.Fn += bo.knA * vn * dt_;
        bo.Fs += -bo.ksA * vt * dt_;
        bo.Mt += -bo.ksJ * twist * dt_;
        bo.Mb += -bo.knI * bend * dt_;

        double sigMax = bo.Fn / bo.A + bo.Mb.norm() * bo.Rb / bo.Ib;
        double tau    = bo.Fs.norm() / bo.A + std::abs(bo.Mt) * bo.Rb / bo.Jb;
        double tauC   = bo.tc0 + tanPhiB_ * std::max(0.0, -bo.Fn / bo.A);
        if (sigMax > bo.sc || tau > tauC) {
            bo.broken = true;
            bo.tBreak = t_;
            ++nBroken_;
            --nIntact_[bo.i];
            --nIntact_[bo.j];
            continue;
        }

        Eigen::Vector3d Fj = -bo.Fn * n + bo.Fs;
        pj.f += Fj;
        pi.f -= Fj;
        Eigen::Vector3d couple = -0.5 * L * n.cross(bo.Fs);
        Eigen::Vector3d M = bo.Mb + bo.Mt * n;
        pi.tq += couple - M;
        pj.tq += couple + M;
    }
}

// ---------------------------------------------------------------------------
// Frictional contact (in parallel with bonds; alone once broken): linear
// normal spring + viscous damping (compression only), vector tangential
// spring with history (rotated into the contact plane each step), Coulomb cap
// mu*fn. Torques use the lever arms to the actual contact point.
// ---------------------------------------------------------------------------
void Dem3dSolver::contactForces() {
    tangNew_.clear();
    for (int i = 0; i < (int)p_.size(); ++i) {
        int ci = std::clamp(int((p_[i].x.x() - gmin_.x()) / cell_), 0, gx_ - 1);
        int cj = std::clamp(int((p_[i].x.y() - gmin_.y()) / cell_), 0, gy_ - 1);
        int ck = std::clamp(int((p_[i].x.z() - gmin_.z()) / cell_), 0, gz_ - 1);
        // Half stencil: the 13 "forward" neighbour cells take every j, the
        // home cell takes j > i — each unordered pair is visited exactly once
        // instead of twice, halving the linked-list traffic.
        for (int dk = -1; dk <= 1; ++dk)
        for (int dj = -1; dj <= 1; ++dj)
        for (int di = -1; di <= 1; ++di) {
            bool home = (dk == 0 && dj == 0 && di == 0);
            if (dk < 0 || (dk == 0 && (dj < 0 || (dj == 0 && di < 0))))
                continue;                      // backward half: already visited
            int cx = ci + di, cy = cj + dj, cz = ck + dk;
            if (cx < 0 || cy < 0 || cz < 0 || cx >= gx_ || cy >= gy_ || cz >= gz_)
                continue;
            for (int j = head_[(cz * gy_ + cy) * gx_ + cx]; j >= 0; j = nxt_[j]) {
                if (home && j <= i) continue;
                Part& pi = p_[i];
                Part& pj = p_[j];
                Eigen::Vector3d d = pj.x - pi.x;
                double dist2 = d.squaredNorm();
                double rsum = pi.r + pj.r;
                if (dist2 >= rsum * rsum) continue;
                double dist = std::sqrt(dist2);
                double delta = rsum - dist;
                // Resting HCP neighbours sit at exactly 2r: floating-point
                // noise puts half of them "in contact" with delta ~ 1e-16 m,
                // churning the tangential-spring map for zero physics. Skip
                // sub-nanometre overlaps (real contacts are micrometres+).
                if (delta <= 1e-12) continue;
                Eigen::Vector3d n = d / dist;

                double li = pi.r - 0.5 * delta;
                double lj = pj.r - 0.5 * delta;
                Eigen::Vector3d vrel =
                    (pj.v - pi.v) - (li * pi.w + lj * pj.w).cross(n);
                double vn = vrel.dot(n);
                Eigen::Vector3d vt = vrel - vn * n;

                double meff = pi.m * pj.m / (pi.m + pj.m);
                double cn = 2.0 * xiC_ * std::sqrt(knC_ * meff);
                double fn = knC_ * delta - cn * vn;
                if (fn < 0) fn = 0;

                uint64_t key = pairKey(i, j);
                auto it = tang_.find(key);
                Eigen::Vector3d fs =
                    (it == tang_.end()) ? Eigen::Vector3d::Zero() : it->second;
                fs -= fs.dot(n) * n;              // rotate into contact plane
                fs += -ksC_ * vt * dt_;
                double fmax = mu_ * fn;
                double fsn = fs.norm();
                if (fsn > fmax && fsn > 0) fs *= fmax / fsn;
                tangNew_[key] = fs;

                Eigen::Vector3d Fj = fn * n + fs;
                pj.f += Fj;
                pi.f -= Fj;
                Eigen::Vector3d nxfs = n.cross(fs);
                pi.tq += -li * nxfs;
                pj.tq += -lj * nxfs;
            }
        }
    }
    tang_.swap(tangNew_);
}

// Rigid spherical tool + optional rigid bottom wall (z = 0), same contact law.
void Dem3dSolver::toolAndWallForces() {
    toolTangNew_.clear();
    for (int i = 0; i < (int)p_.size(); ++i) {
        Part& q = p_[i];

        if (bottomWall_) {
            double delta = q.r - q.x.z();
            if (delta > 0) {
                double cn = 2.0 * xiC_ * std::sqrt(knC_ * q.m);
                double fn = knC_ * delta - cn * q.v.z();
                if (fn < 0) fn = 0;
                double l = q.r - 0.5 * delta;
                Eigen::Vector3d nW(0, 0, 1);   // wall normal (into the block)
                Eigen::Vector3d vt = q.v - q.v.z() * nW - (l * q.w).cross(nW);
                Eigen::Vector3d fs = Eigen::Vector3d::Zero();
                double vtn = vt.norm();
                if (vtn > 0)
                    fs = -mu_ * fn * std::tanh(vtn / 1e-3) * vt / vtn;
                q.f += fn * nW + fs;
                q.tq += -l * nW.cross(fs);
            }
        }

        if (scen_ == Scenario::TENSION) continue;

        Eigen::Vector3d n;                     // push direction on the particle
        double delta;
        if (tool_.flat) {                      // flat bottom face pushing down
            double rx = q.x.x() - tool_.x.x(), ry = q.x.y() - tool_.x.y();
            if (std::sqrt(rx * rx + ry * ry) > tool_.radius + q.r) continue;
            delta = q.x.z() + q.r - tool_.x.z();
            if (delta <= 0) continue;
            n = {0.0, 0.0, -1.0};
        } else {
            Eigen::Vector3d d = q.x - tool_.x;
            double dist = d.norm();
            delta = (q.r + tool_.radius) - dist;
            if (delta <= 0) continue;
            n = d / dist;                      // from tool centre to particle
        }

        Eigen::Vector3d vrel = q.v - tool_.v;  // tool does not spin
        double l = q.r - 0.5 * delta;
        vrel -= (l * q.w).cross(-n);           // particle surface point velocity
        double vn = vrel.dot(n);
        Eigen::Vector3d vt = vrel - vn * n;

        double cn = 2.0 * xiC_ * std::sqrt(knC_ * q.m);
        double fn = knC_ * delta - cn * vn;
        if (fn < 0) fn = 0;

        auto it = toolTang_.find(i);
        Eigen::Vector3d fs =
            (it == toolTang_.end()) ? Eigen::Vector3d::Zero() : it->second;
        fs -= fs.dot(n) * n;
        fs += -ksC_ * vt * dt_;
        double fmax = mu_ * fn;
        double fsn = fs.norm();
        if (fsn > fmax && fsn > 0) fs *= fmax / fsn;
        toolTangNew_[i] = fs;

        Eigen::Vector3d F = fn * n + fs;       // on the particle
        q.f += F;
        q.tq += -l * (-n).cross(fs);           // lever from particle centre
        tool_.F -= F;
    }
    toolTang_.swap(toolTangNew_);
}

void Dem3dSolver::integrate() {
    for (auto& q : p_) {
        if (q.flag == FIXED) {
            // gripFree: frictionless grips — hold the axial dof only (see
            // the 2D FDEM comment: fully clamped grips concentrate stress
            // at the corners and decide where the specimen breaks)
            if (gripFree_) {
                q.v.x() += (dt_ / q.m) * q.f.x();
                q.v.y() += (dt_ / q.m) * q.f.y();
                q.v.z() = 0.0;
                q.w.setZero();
                q.x += dt_ * q.v;
            } else { q.v.setZero(); q.w.setZero(); }
            continue;
        }
        if (q.flag == PRESCRIBED) {
            // pullRamp: cosine rise of the grip velocity — a stepped grip
            // launches a transient that breaks the first row under the
            // grip whatever the strength map says (2D lesson)
            double vg = pullV_;
            if (pullRamp_ > 0.0 && t_ < pullRamp_)
                vg *= 0.5 * (1.0 - std::cos(M_PI * t_ / pullRamp_));
            if (gripFree_) {
                q.v.x() += (dt_ / q.m) * q.f.x();
                q.v.y() += (dt_ / q.m) * q.f.y();
            } else { q.v.x() = 0.0; q.v.y() = 0.0; }
            q.v.z() = vg;
            q.w.setZero();
            q.x += dt_ * q.v;
            continue;
        }
        std::size_t ib = (std::size_t)(&q - p_.data());
        if (nIntact_[ib] > 0) {                // viscous-spring quiet boundary
            for (int a = 0; a < 3; ++a)
                if (kAbs_[ib](a) > 0)
                    q.f(a) -= kAbs_[ib](a) * (q.x(a) - xAnchor_[ib](a));
        }
        if (damping_ > 0) {                    // Cundall local damping
            for (int a = 0; a < 3; ++a) {
                q.f(a)  -= damping_ * std::abs(q.f(a)) *
                           (q.v(a) > 0 ? 1.0 : (q.v(a) < 0 ? -1.0 : 0.0));
                q.tq(a) -= damping_ * std::abs(q.tq(a)) *
                           (q.w(a) > 0 ? 1.0 : (q.w(a) < 0 ? -1.0 : 0.0));
            }
        }
        q.v += (dt_ / q.m) * q.f;
        q.w += (dt_ / q.I) * q.tq;
        if (nIntact_[ib] > 0) {                // implicit Lysmer dashpot
            for (int a = 0; a < 3; ++a)
                if (cAbs_[ib](a) > 0)
                    q.v(a) /= 1.0 + dt_ * cAbs_[ib](a) / q.m;
        }
        q.x += dt_ * q.v;
    }
    if (scen_ != Scenario::TENSION) tool_.integrate(dt_);
}

// ---------------------------------------------------------------------------
// Fragments: connected components of the intact bond network (BFS).
// ---------------------------------------------------------------------------
void Dem3dSolver::computeFragments() {
    std::vector<std::vector<int>> adj(p_.size());
    for (const auto& bo : b_)
        if (!bo.broken) {
            adj[bo.i].push_back(bo.j);
            adj[bo.j].push_back(bo.i);
        }
    std::fill(fragId_.begin(), fragId_.end(), -1);
    int nid = 0;
    std::vector<std::pair<int, double>> sizes;   // (id, mass)
    for (int s = 0; s < (int)p_.size(); ++s) {
        if (fragId_[s] >= 0) continue;
        double mass = 0;
        std::queue<int> qu;
        qu.push(s);
        fragId_[s] = nid;
        while (!qu.empty()) {
            int u = qu.front(); qu.pop();
            mass += p_[u].m;
            for (int v : adj[u])
                if (fragId_[v] < 0) { fragId_[v] = nid; qu.push(v); }
        }
        sizes.push_back({nid, mass});
        ++nid;
    }
    // relabel so fragment 0 is the most massive (the "main body")
    std::sort(sizes.begin(), sizes.end(),
              [](auto& a, auto& b) { return a.second > b.second; });
    std::vector<int> rank(nid);
    for (int k = 0; k < nid; ++k) rank[sizes[k].first] = k;
    for (auto& f : fragId_) f = rank[f];
    nFrag_ = nid;
    double mDet = 0;
    for (int i = 0; i < (int)p_.size(); ++i)
        if (fragId_[i] != 0) mDet += p_[i].m;
    detachedVol_ = mDet / mat_.rho;
}

// ===========================================================================
// Output
// ===========================================================================

void Dem3dSolver::writeFrame(int frame) {
    computeFragments();

    std::vector<Eigen::Vector3d> pts(p_.size()), vel(p_.size());
    std::vector<double> rad(p_.size()), frag(p_.size()), spd(p_.size());
    for (std::size_t i = 0; i < p_.size(); ++i) {
        pts[i] = p_[i].x;
        rad[i] = p_[i].r;
        frag[i] = fragId_[i];
        spd[i] = p_[i].v.norm();
        vel[i] = p_[i].v;
    }
    char name[64];
    std::snprintf(name, sizeof(name), "/dem3d_particles_%04d.vtu", frame);
    vtk::writeParticles(out_ + name, pts,
                        {{"radius", &rad}, {"fragment", &frag}, {"speed", &spd}},
                        {{"velocity", &vel}});

    std::vector<std::array<int, 2>> lines;
    std::vector<double> state, tb;
    for (const auto& bo : b_) {
        lines.push_back({bo.i, bo.j});
        state.push_back(bo.broken ? 1.0 : 0.0);
        tb.push_back(bo.tBreak);
    }
    std::snprintf(name, sizeof(name), "/dem3d_bonds_%04d.vtu", frame);
    vtk::writeLines(out_ + name, pts, lines, {{"state", &state}, {"tBreak", &tb}});

    std::ofstream fm(out_ + "/frames.csv",
                     frame == 0 ? std::ios::trunc : std::ios::app);
    if (frame == 0) fm << "frame,t,toolX,toolY,toolZ\n";
    fm << frame << "," << t_ << "," << tool_.x.x() << "," << tool_.x.y()
       << "," << tool_.x.z() << "\n";
}

void Dem3dSolver::historyHeader(std::ostream& os) const {
    if (scen_ == Scenario::TENSION) { os << "t,gripFz,sigma,sigmaPeak,nBroken\n"; return; }
    os << "t,toolFx,toolFy,toolFz,toolX,toolY,toolZ,toolVx,toolVy,toolVz,"
          "work,toolKE,nBroken,nFrag,detachedVol,specificEnergy\n";
}

void Dem3dSolver::historyRow(std::ostream& os) const {
    if (scen_ == Scenario::TENSION) {
        os << t_ << "," << gripF_.z() << ","
           << std::abs(gripF_.z()) / (W_ * D_) << ","
           << sigmaPeak_ << "," << nBroken_ << "\n";
        return;
    }
    double Es = detachedVol_ > 0 ? work_ / detachedVol_ : 0.0;
    os << t_ << "," << tool_.F.x() << "," << tool_.F.y() << "," << tool_.F.z()
       << "," << tool_.x.x() << "," << tool_.x.y() << "," << tool_.x.z()
       << "," << tool_.v.x() << "," << tool_.v.y() << "," << tool_.v.z()
       << "," << work_ << "," << tool_.ke() << "," << nBroken_ << ","
       << nFrag_ << "," << detachedVol_ << "," << Es << "\n";
}

void Dem3dSolver::finalize() {
    computeFragments();

    std::map<int, double> mass;
    std::map<int, int> np;
    for (int i = 0; i < (int)p_.size(); ++i) {
        mass[fragId_[i]] += p_[i].m;
        np[fragId_[i]] += 1;
    }
    std::ofstream ff(out_ + "/dem3d_fragments.csv");
    ff << "fragment,nParticles,mass\n";
    for (auto& [id, mm] : mass) ff << id << "," << np[id] << "," << mm << "\n";

    std::ofstream fp(out_ + "/dem3d_final_particles.csv");
    fp << "x,y,z,r,fragment\n";
    for (int i = 0; i < (int)p_.size(); ++i)
        fp << p_[i].x.x() << "," << p_[i].x.y() << "," << p_[i].x.z() << ","
           << p_[i].r << "," << fragId_[i] << "\n";

    std::cout << "\n[DEM3D] ---- summary ----\n";
    if (scen_ == Scenario::TENSION) {
        double expected = cfg_.getd("bondTensile", mat_.ft)
                          * lambda_ * lambda_ * M_PI / 4.0;
        double err = 100.0 * (sigmaPeak_ - expected) / expected;
        bool pass = std::abs(err) < 5.0;
        std::cout << "[DEM3D] tension verification (simple cubic, load-aligned):\n"
                  << "[DEM3D]   peak macro stress = " << sigmaPeak_ / 1e6
                  << " MPa, expected (pi/4) lambda^2 sigma_c = " << expected / 1e6
                  << " MPa, error = " << err << " %  ["
                  << (pass ? "PASS" : "FAIL") << "]\n"
                  << "[DEM3D]   broken bonds = " << nBroken_ << "\n";
        return;
    }
    double Es = detachedVol_ > 0 ? work_ / detachedVol_ : 0.0;
    std::cout << "[DEM3D] peak tool force   : " << peakF_ << " N\n"
              << "[DEM3D] tool work output  : " << work_ << " J";
    if (tool_.free)
        std::cout << "  (tool KE loss: " << toolKE0_ - tool_.ke() << " J)";
    std::cout << "\n[DEM3D] broken bonds      : " << nBroken_ << " / " << b_.size()
              << "\n[DEM3D] fragments         : " << nFrag_
              << " (detached vol " << detachedVol_ << " m^3)"
              << "\n[DEM3D] specific energy   : " << Es << " J/m^3\n";
}

} // namespace rockim
