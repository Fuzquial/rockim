#include "rockim/FemSolver.hpp"
#include "rockim/VtkWriter.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

namespace rockim {

// ===========================================================================
// Construction / setup
// ===========================================================================

FemSolver::FemSolver(const Config& cfg, std::string outDir)
    : cfg_(cfg), out_(std::move(outDir)) {}

void FemSolver::init() {
    mat_ = Material::from(cfg_);
    W_ = cfg_.getd("W", W_);
    H_ = cfg_.getd("H", H_);
    nx_ = cfg_.geti("nx", nx_);
    ny_ = cfg_.geti("ny", ny_);
    T_  = cfg_.getd("T", 200e-6);

    std::string s = cfg_.gets("scenario", "percussion");
    if      (s == "percussion") scen_ = Scenario::PERCUSSION;
    else if (s == "shear")      scen_ = Scenario::SHEAR;
    else if (s == "bar_wave")   scen_ = Scenario::BAR_WAVE;
    else throw std::runtime_error("FemSolver: unknown scenario '" + s + "'");

    damageOn_  = cfg_.getb("damage", scen_ != Scenario::BAR_WAVE);
    erodeD_    = cfg_.getd("erodeD", erodeD_);
    strainCap_ = cfg_.getd("strainCap", strainCap_);

    Dm_ = mat_.Dmat();
    mat_.dpParams(dpAlpha_, dpK_);
    // Damage initiation thresholds expressed as equivalent strains:
    //   tension: kappa0_t = ft / E        (Rankine, strain at first cracking)
    //   shear  : kappa0_s = k_DP / (2G)   (DP stress measure mapped to strain)
    kappa0T_ = mat_.ft / mat_.E;
    kappa0S_ = dpK_ / (2.0 * mat_.G());

    buildMesh();
    placeTool();

    // Penalty contact stiffness. kp ~ E * t gives a contact compliance of the
    // same order as one element row, i.e. stiff enough not to pollute the
    // response, soft enough not to dominate the stable time step (it is
    // accounted for in computeStableDt() anyway).
    kp_  = cfg_.getd("kpFactor", 1.0) * mat_.E * thk_;
    xiC_ = cfg_.getd("contactXi", xiC_);
    muC_ = cfg_.getd("contactMu", muC_);

    computeStableDt();
    toolKE0_ = tool_.ke();

    if (scen_ == Scenario::BAR_WAVE) {
        barV0_ = cfg_.getd("barV0", 1.0);
        double frac = cfg_.getd("barGaugeFrac", 0.8);
        gaugeX_ = frac * W_;
        double dx = W_ / nx_;
        for (int i = 0; i < (int)X0_.size(); ++i) {
            if (X0_[i].x() < 0.5 * dx) barBcNodes_.push_back(i);
            if (std::abs(X0_[i].x() - gaugeX_) < 0.51 * dx) gaugeNodes_.push_back(i);
        }
    }

    std::cout << "[FEM] " << el_.size() << " CST elements, " << X0_.size()
              << " nodes, dt = " << dt_ << " s, steps = "
              << (long)std::ceil(T_ / dt_) << "\n";
}

void FemSolver::buildMesh() {
    // Structured grid, each quad split into two CSTs with alternating
    // diagonals ("union-jack"-ish) to reduce directional mesh bias in the
    // crack patterns.
    int nnx = nx_ + 1, nny = ny_ + 1;
    double dx = W_ / nx_, dy = H_ / ny_;

    X0_.resize((std::size_t)nnx * nny);
    for (int j = 0; j < nny; ++j)
        for (int i = 0; i < nnx; ++i)
            X0_[(std::size_t)j * nnx + i] = {i * dx, j * dy};

    u_.assign(X0_.size(), Eigen::Vector2d::Zero());
    v_.assign(X0_.size(), Eigen::Vector2d::Zero());
    f_.assign(X0_.size(), Eigen::Vector2d::Zero());
    m_.assign(X0_.size(), 0.0);
    fix_.assign(X0_.size(), 0);
    active_.assign(X0_.size(), 1);

    auto nid = [nnx](int i, int j) { return j * nnx + i; };

    for (int j = 0; j < ny_; ++j) {
        for (int i = 0; i < nx_; ++i) {
            int a = nid(i, j), b = nid(i + 1, j), c = nid(i + 1, j + 1), d = nid(i, j + 1);
            if ((i + j) % 2 == 0) {           // diagonal a-c
                el_.push_back({{a, b, c}});
                el_.push_back({{a, c, d}});
            } else {                          // diagonal b-d
                el_.push_back({{a, b, d}});
                el_.push_back({{b, c, d}});
            }
        }
    }

    // Geometry, B matrices, lumped mass (row-sum: A*rho*t/3 per node)
    for (auto& e : el_) {
        const auto& p1 = X0_[e.n[0]];
        const auto& p2 = X0_[e.n[1]];
        const auto& p3 = X0_[e.n[2]];
        double x1 = p1.x(), y1 = p1.y(), x2 = p2.x(), y2 = p2.y(), x3 = p3.x(), y3 = p3.y();
        double det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
        e.A = 0.5 * det;                       // CCW ordering guarantees det > 0
        double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
        double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;
        double inv = 1.0 / det;                // = 1/(2A)
        e.B.setZero();
        e.B(0, 0) = b1 * inv; e.B(0, 2) = b2 * inv; e.B(0, 4) = b3 * inv;
        e.B(1, 1) = c1 * inv; e.B(1, 3) = c2 * inv; e.B(1, 5) = c3 * inv;
        e.B(2, 0) = c1 * inv; e.B(2, 1) = b1 * inv;
        e.B(2, 2) = c2 * inv; e.B(2, 3) = b2 * inv;
        e.B(2, 4) = c3 * inv; e.B(2, 5) = b3 * inv;

        double l12 = std::hypot(x2 - x1, y2 - y1);
        double l23 = std::hypot(x3 - x2, y3 - y2);
        double l31 = std::hypot(x1 - x3, y1 - y3);
        double lmax = std::max({l12, l23, l31});
        e.hMin = 2.0 * e.A / lmax;             // smallest altitude
        e.lc   = std::sqrt(2.0 * e.A);         // crack-band width ~ grid spacing

        double mNode = mat_.rho * e.A * thk_ / 3.0;
        for (int k = 0; k < 3; ++k) m_[e.n[k]] += mNode;
    }

    // Boundary conditions: bottom fully fixed (bedrock support) except for
    // the bar-wave test, which needs a free-free bar with a velocity BC.
    if (scen_ != Scenario::BAR_WAVE) {
        for (int i = 0; i < (int)X0_.size(); ++i)
            if (X0_[i].y() < 1e-12) fix_[i] = 3;
        if (cfg_.getb("fixSides", false))
            for (int i = 0; i < (int)X0_.size(); ++i)
                if (X0_[i].x() < 1e-12 || X0_[i].x() > W_ - 1e-12) fix_[i] |= 1;
    }

    // ------------------------------------------------------------------
    // Lysmer-Kuhlemeyer absorbing (quiet) boundaries.
    // A plane wave hitting a boundary of normal n is absorbed exactly (and
    // oblique incidence approximately) by viscous tractions matching the 1D
    // impedances:  t_n = rho c_p v_n,  t_t = rho c_s v_t.  Lumped on each
    // boundary node over its tributary edge length, this gives one dashpot
    // per direction. 'absorbing = sides' treats the lateral faces (bottom
    // stays fixed), 'absorbing = all' also replaces the bottom support by
    // dashpots. Applied implicitly in the velocity update (see integrate()),
    // so it is unconditionally stable and never controls the time step.
    // ------------------------------------------------------------------
    cAbsX_.assign(X0_.size(), 0.0);
    cAbsY_.assign(X0_.size(), 0.0);
    kAbsX_.assign(X0_.size(), 0.0);
    kAbsY_.assign(X0_.size(), 0.0);
    std::string ab = cfg_.gets("absorbing", "none");
    if (scen_ != Scenario::BAR_WAVE && ab != "none") {
        if (ab != "sides" && ab != "all")
            throw std::runtime_error("absorbing must be none | sides | all");
        // Pure viscous (Lysmer) boundaries have ZERO static stiffness: under
        // the slow part of the loading the block behaves as if free-floating,
        // which puts spurious bending tension at the bottom centre. The
        // viscous-SPRING boundary (Deeks & Randolph 1994) restores the static
        // support of the truncated half-space with a spring in parallel:
        //   k_n = G/R, k_t = G/(2R) per unit boundary length,
        // R being the spreading distance from the loaded zone to that face.
        // The spring is orders softer than the medium, so wave absorption is
        // essentially unaffected.
        double zP = mat_.rho * mat_.cP() * thk_;   // impedance per unit length
        double zS = mat_.rho * mat_.cS() * thk_;
        double G  = mat_.E / (2.0 * (1.0 + mat_.nu));
        double sF = cfg_.getd("absorbSpringFactor", 1.0);   // 0 = pure Lysmer
        double Rside = cfg_.getd("absorbSpringR", 0.5 * W_);
        double Rbot  = cfg_.getd("absorbSpringR", H_);
        for (int j = 0; j < nny; ++j) {            // lateral faces, n = +-x
            double L = dy * ((j == 0 || j == nny - 1) ? 0.5 : 1.0);
            for (int i : {nid(0, j), nid(nnx - 1, j)}) {
                cAbsX_[i] += zP * L;
                cAbsY_[i] += zS * L;
                kAbsX_[i] += sF * G / Rside * L * thk_;
                kAbsY_[i] += sF * G / (2.0 * Rside) * L * thk_;
            }
        }
        if (ab == "all") {                         // bottom face, n = -y
            for (int i = 0; i < nnx; ++i) {
                double L = dx * ((i == 0 || i == nnx - 1) ? 0.5 : 1.0);
                cAbsY_[nid(i, 0)] += zP * L;
                cAbsX_[nid(i, 0)] += zS * L;
                kAbsY_[nid(i, 0)] += sF * G / Rbot * L * thk_;
                kAbsX_[nid(i, 0)] += sF * G / (2.0 * Rbot) * L * thk_;
            }
            for (int i = 0; i < (int)X0_.size(); ++i)
                if (X0_[i].y() < 1e-12) fix_[i] = 0;   // dashpots replace the support
        }
    }
}

void FemSolver::placeTool() {
    if (scen_ == Scenario::BAR_WAVE) return;
    tool_.mass = cfg_.getd("toolMass", 5.0);
    double gap = cfg_.getd("toolGap", 1e-4);

    std::string sh = cfg_.gets("toolShape", scen_ == Scenario::SHEAR ? "disc" : "flat");
    tool_.shape = (sh == "disc") ? Tool::Shape::DISC : Tool::Shape::FLAT;
    tool_.width  = cfg_.getd("toolWidth", 0.02);
    tool_.radius = cfg_.getd("toolRadius", 0.01);

    if (scen_ == Scenario::PERCUSSION) {
        tool_.motion = Tool::Motion::FREE;
        double vImp = cfg_.getd("impactSpeed", 15.0);
        double xc = cfg_.getd("toolX", 0.5 * W_);
        if (tool_.shape == Tool::Shape::FLAT) tool_.x = {xc, H_ + gap};
        else                                  tool_.x = {xc, H_ + tool_.radius + gap};
        tool_.v = {0.0, -vImp};
    } else {  // SHEAR: prescribed lateral pass at a fixed depth of cut
        tool_.motion = Tool::Motion::PRESCRIBED;
        tool_.shape  = Tool::Shape::DISC;   // lateral cutting needs a 2D normal
        double depth = cfg_.getd("cutDepth", 0.003);
        double vCut  = cfg_.getd("cutSpeed", 10.0);
        tool_.x = {-tool_.radius - gap, H_ - depth + tool_.radius};
        tool_.v = {vCut, 0.0};
    }
}

// ---------------------------------------------------------------------------
// Stable explicit time step.
//
// Central difference / velocity-Verlet is conditionally stable:
//     dt <= dt_crit = 2 / omega_max.
// For a lumped-mass CST mesh, omega_max is bounded element-wise and the usual
// CFL estimate is
//     dt_crit ~ min_e ( h_min(e) / c_p ),
// where h_min is the smallest element altitude and c_p the plane-strain
// dilatational wave speed (the fastest wave). The penalty contact adds a
// spring k_p on single nodes, i.e. an extra frequency sqrt(k_p/m_node), so
// its own limit 2*sqrt(m_min/k_p) is taken into account too. A safety factor
// (default 0.7) covers the softening branch of the damage law and the
// contact nonlinearity.
// ---------------------------------------------------------------------------
void FemSolver::computeStableDt() {
    double cP = mat_.cP();
    double hMin = 1e30;
    for (const auto& e : el_) hMin = std::min(hMin, e.hMin);
    double dtMesh = hMin / cP;

    double mMin = 1e30;
    for (double mm : m_) mMin = std::min(mMin, mm);
    double dtContact = (kp_ > 0) ? 2.0 * std::sqrt(mMin / kp_) : 1e30;

    double cfl = cfg_.getd("cfl", 0.7);
    dt_ = cfl * std::min(dtMesh, dtContact);
}

// ===========================================================================
// Time stepping
// ===========================================================================

void FemSolver::step() {
    for (auto& fv : f_) fv.setZero();

    if (scen_ != Scenario::BAR_WAVE) applyContact();
    internalForcesAndDamage();
    if (activeDirty_) refreshActiveNodes();
    integrate();

    t_ += dt_;
}

// ---------------------------------------------------------------------------
// Rigid tool <-> mesh contact, node-to-rigid-surface penalty:
//     f_n = k_p * penetration + c * v_rel_n,   c = 2 xi sqrt(k_p m_node)
// clamped to compression only (no adhesion), plus regularized Coulomb
// friction f_t = -mu |f_n| tanh(v_slip / v_reg).
// ---------------------------------------------------------------------------
void FemSolver::applyContact() {
    tool_.resetForce();
    for (int i = 0; i < (int)X0_.size(); ++i) {
        if (!active_[i]) continue;             // orphan (fully eroded) nodes fly free
        Eigen::Vector2d p = X0_[i] + u_[i];
        Eigen::Vector2d Fc = Eigen::Vector2d::Zero();

        if (tool_.shape == Tool::Shape::FLAT) {
            if (std::abs(p.x() - tool_.x.x()) > 0.5 * tool_.width) continue;
            double pen = p.y() - tool_.x.y();  // node above the bottom face
            if (pen <= 0) continue;
            double vreln = v_[i].y() - tool_.v.y();
            double c = 2.0 * xiC_ * std::sqrt(kp_ * m_[i]);
            double fn = kp_ * pen + c * vreln;
            if (fn < 0) fn = 0;
            double vslip = v_[i].x() - tool_.v.x();
            double ftang = -muC_ * fn * std::tanh(vslip / vReg_);
            Fc = {ftang, -fn};                 // tool pushes the node down
        } else {                               // DISC
            Eigen::Vector2d d = p - tool_.x;
            double dist = d.norm();
            if (dist >= tool_.radius || dist < 1e-14) continue;
            Eigen::Vector2d n = d / dist;      // outward normal (pushes node out)
            Eigen::Vector2d tdir(-n.y(), n.x());
            double pen = tool_.radius - dist;
            Eigen::Vector2d vrel = v_[i] - tool_.v;
            double c = 2.0 * xiC_ * std::sqrt(kp_ * m_[i]);
            double fn = kp_ * pen - c * vrel.dot(n);
            if (fn < 0) fn = 0;
            double ftang = -muC_ * fn * std::tanh(vrel.dot(tdir) / vReg_);
            Fc = fn * n + ftang * tdir;
        }

        f_[i] += Fc;
        tool_.F -= Fc;                         // Newton's third law
    }
    peakF_ = std::max(peakF_, tool_.F.norm());
    // Energy delivered by the tool (thrust/cutting work). Tool-side
    // bookkeeping deliberately includes the interface friction and dashpot
    // dissipation: that is part of the drilling energy input.
    work_ += -tool_.F.dot(tool_.v) * dt_;
}

// ---------------------------------------------------------------------------
// Internal forces + damage.
//
// Small-strain, total-strain format:  eps = B u_e,  effective stress
// s_eff = D eps (with s_zz = nu (s_xx + s_yy) in plane strain), nominal
// stress s = (1 - D) s_eff, and  f_int = B^T s A t.
//
// Damage model (engineering-grade continuum damage, see README limitations):
//  * Tension (Rankine): driver  e_t = max principal effective stress / E,
//    threshold kappa0_t = ft/E.
//  * Shear (Drucker-Prager): driver e_s = (sqrt(J2) + alpha I1)/(2G),
//    threshold kappa0_s = k/(2G). Pure hydrostatic compression never damages
//    (alpha I1 << 0 keeps the driver below threshold), which is the point of
//    using DP instead of a von Mises-type surface for rock.
//  * Both use the exponential softening law
//        D = 1 - (kappa0/kappa) exp( -(kappa - kappa0) / eps_f )
//    with eps_f tied to the fracture energy and the element size through the
//    crack-band scaling (Bazant-Oh):  eps_f = G_f / (f * l_c). This keeps the
//    dissipated energy per unit crack area ~ G_f independent of the mesh.
//  * D = max(D_t, D_s), monotonic. Element erosion at D >= erodeD (or at a
//    strain cap, to kill inverted junk) removes the element => cracks and
//    material removal appear as bands/craters of eroded elements.
// ---------------------------------------------------------------------------
void FemSolver::internalForcesAndDamage() {
    for (auto& e : el_) {
        if (e.eroded) continue;

        Eigen::Matrix<double, 6, 1> ue;
        for (int k = 0; k < 3; ++k) {
            ue(2 * k)     = u_[e.n[k]].x();
            ue(2 * k + 1) = u_[e.n[k]].y();
        }
        Eigen::Vector3d eps = e.B * ue;

        // strain cap: erode grossly distorted elements before they misbehave
        if (std::abs(eps(0)) > strainCap_ || std::abs(eps(1)) > strainCap_ ||
            std::abs(eps(2)) > 2.0 * strainCap_) {
            erode(e);
            continue;
        }

        Eigen::Vector3d sEff = Dm_ * eps;
        double szz = mat_.nu * (sEff(0) + sEff(1));   // plane strain

        if (damageOn_) {
            updateDamage(e, sEff, szz);
            if (e.D >= erodeD_) { erode(e); continue; }
        }

        Eigen::Vector3d sig = (1.0 - e.D) * sEff;

        // stored for output
        double sm = (sig(0) + sig(1) + (1.0 - e.D) * szz) / 3.0;
        double dxx = sig(0) - sm, dyy = sig(1) - sm, dzz = (1.0 - e.D) * szz - sm;
        e.smean = sm;
        e.svm = std::sqrt(1.5 * (dxx * dxx + dyy * dyy + dzz * dzz) + 3.0 * sig(2) * sig(2));

        Eigen::Matrix<double, 6, 1> fe = e.B.transpose() * sig * (e.A * thk_);
        for (int k = 0; k < 3; ++k) {
            f_[e.n[k]].x() -= fe(2 * k);
            f_[e.n[k]].y() -= fe(2 * k + 1);
        }
    }
}

void FemSolver::updateDamage(Elem& e, const Eigen::Vector3d& sEff, double szz) {
    double sxx = sEff(0), syy = sEff(1), txy = sEff(2);

    // In-plane principal + out-of-plane
    double ctr = 0.5 * (sxx + syy);
    double R   = std::sqrt(0.25 * (sxx - syy) * (sxx - syy) + txy * txy);
    double s1  = std::max(ctr + R, szz);

    // Invariants (tension positive), z included
    double I1 = sxx + syy + szz;
    double sm = I1 / 3.0;
    double dxx = sxx - sm, dyy = syy - sm, dzz = szz - sm;
    double J2 = 0.5 * (dxx * dxx + dyy * dyy + dzz * dzz) + txy * txy;

    // --- tension (Rankine) ---
    double et = s1 / mat_.E;
    if (et > e.kappaT) e.kappaT = et;
    double Dt = 0.0;
    if (e.kappaT > kappa0T_) {
        double efT = std::max(mat_.Gf / (mat_.ft * e.lc), 1.5 * kappa0T_);
        Dt = 1.0 - (kappa0T_ / e.kappaT) * std::exp(-(e.kappaT - kappa0T_) / efT);
    }

    // --- shear (Drucker-Prager) ---
    double q = std::sqrt(J2) + dpAlpha_ * I1;
    double es = q / (2.0 * mat_.G());
    if (es > e.kappaS) e.kappaS = es;
    double Ds = 0.0;
    if (e.kappaS > kappa0S_) {
        double efS = std::max(mat_.gfShearFactor * mat_.Gf / (dpK_ * e.lc), 1.5 * kappa0S_);
        Ds = 1.0 - (kappa0S_ / e.kappaS) * std::exp(-(e.kappaS - kappa0S_) / efS);
    }

    double Dn = std::min(std::max(Dt, Ds), 0.999);
    if (Dn > e.D) e.D = Dn;
}

void FemSolver::erode(Elem& e) {
    e.eroded = true;
    e.D = 1.0;
    erodedVol_ += e.A * thk_;
    ++nEroded_;
    activeDirty_ = true;
}

void FemSolver::refreshActiveNodes() {
    std::fill(active_.begin(), active_.end(), 0);
    for (const auto& e : el_)
        if (!e.eroded)
            for (int k = 0; k < 3; ++k) active_[e.n[k]] = 1;
    activeDirty_ = false;
}

// Central-difference / leapfrog kick-drift:
//   v^{n+1/2} = v^{n-1/2} + dt a^n ;  u^{n+1} = u^n + dt v^{n+1/2}
void FemSolver::integrate() {
    for (int i = 0; i < (int)X0_.size(); ++i) {
        // boundary spring (viscous-spring quiet boundary), anchored at X0
        if (active_[i] && (kAbsX_[i] > 0 || kAbsY_[i] > 0)) {
            f_[i].x() -= kAbsX_[i] * u_[i].x();
            f_[i].y() -= kAbsY_[i] * u_[i].y();
        }
        Eigen::Vector2d a = f_[i] / m_[i];
        if (!(fix_[i] & 1)) v_[i].x() += dt_ * a.x(); else v_[i].x() = 0;
        if (!(fix_[i] & 2)) v_[i].y() += dt_ * a.y(); else v_[i].y() = 0;
        // Lysmer dashpot, implicit in v:  m dv/dt = F - c v  discretized as
        //   v <- (v + dt F/m) / (1 + dt c/m)
        // exact impedance at low frequency, unconditionally stable.
        if (cAbsX_[i] > 0) v_[i].x() /= 1.0 + dt_ * cAbsX_[i] / m_[i];
        if (cAbsY_[i] > 0) v_[i].y() /= 1.0 + dt_ * cAbsY_[i] / m_[i];
    }
    if (scen_ == Scenario::BAR_WAVE)
        for (int i : barBcNodes_) v_[i] = {barV0_, 0.0};   // velocity step at x=0

    for (int i = 0; i < (int)X0_.size(); ++i) u_[i] += dt_ * v_[i];

    if (scen_ != Scenario::BAR_WAVE) {
        tool_.integrate(dt_);
        if (tool_.motion == Tool::Motion::FREE && tool_.v.y() > 0 &&
            cfg_.getb("stopOnRebound", false)) { /* optional early stop hook */ }
    }

    if (scen_ == Scenario::BAR_WAVE && tArrive_ < 0)
        for (int i : gaugeNodes_)
            if (std::abs(v_[i].x()) > 0.1 * std::abs(barV0_)) { tArrive_ = t_; break; }
}

// ===========================================================================
// Output
// ===========================================================================

void FemSolver::writeFrame(int frame) {
    std::vector<Eigen::Vector2d> pts(X0_.size());
    for (std::size_t i = 0; i < X0_.size(); ++i) pts[i] = X0_[i] + u_[i];

    std::vector<std::array<int, 3>> tris;
    std::vector<double> dmg, svm, smean;
    for (const auto& e : el_) {
        if (e.eroded) continue;                // eroded elements = open cracks/crater
        tris.push_back(e.n);
        dmg.push_back(e.D);
        svm.push_back(e.svm);
        smean.push_back(e.smean);
    }

    char name[64];
    std::snprintf(name, sizeof(name), "/fem_%04d.vtu", frame);
    vtk::writeTriMesh(out_ + name, pts, tris,
                      {{"damage", &dmg}, {"vonMises", &svm}, {"meanStress", &smean}},
                      {{"velocity", &v_}, {"displacement", &u_}});

    // frame index -> (time, tool pose): lets post-processing (make_gif.py)
    // place the rigid tool without re-deriving the frame stride.
    std::ofstream fm(out_ + "/frames.csv",
                     frame == 0 ? std::ios::trunc : std::ios::app);
    if (frame == 0) fm << "frame,t,toolX,toolY\n";
    fm << frame << "," << t_ << "," << tool_.x.x() << "," << tool_.x.y() << "\n";
}

void FemSolver::historyHeader(std::ostream& os) const {
    if (scen_ == Scenario::BAR_WAVE) { os << "t,gaugeVx\n"; return; }
    os << "t,toolFx,toolFy,toolX,toolY,toolVx,toolVy,work,toolKE,"
          "erodedVol,nEroded,specificEnergy\n";
}

void FemSolver::historyRow(std::ostream& os) const {
    if (scen_ == Scenario::BAR_WAVE) {
        double vg = 0;
        for (int i : gaugeNodes_) vg = std::max(vg, std::abs(v_[i].x()));
        os << t_ << "," << vg << "\n";
        return;
    }
    double Es = erodedVol_ > 0 ? work_ / erodedVol_ : 0.0;
    os << t_ << "," << tool_.F.x() << "," << tool_.F.y() << ","
       << tool_.x.x() << "," << tool_.x.y() << ","
       << tool_.v.x() << "," << tool_.v.y() << ","
       << work_ << "," << tool_.ke() << ","
       << erodedVol_ << "," << nEroded_ << "," << Es << "\n";
}

void FemSolver::finalize() {
    // Convenience CSV snapshot of the final damage field (easy to plot)
    std::ofstream fe(out_ + "/fem_final_elements.csv");
    fe << "cx,cy,damage,eroded\n";
    for (const auto& e : el_) {
        Eigen::Vector2d c = Eigen::Vector2d::Zero();
        for (int k = 0; k < 3; ++k) c += (X0_[e.n[k]] + u_[e.n[k]]) / 3.0;
        fe << c.x() << "," << c.y() << "," << e.D << "," << (e.eroded ? 1 : 0) << "\n";
    }

    std::cout << "\n[FEM] ---- summary ----\n";
    if (scen_ == Scenario::BAR_WAVE) {
        double cTh = mat_.cBar();
        if (tArrive_ > 0) {
            double cMe = gaugeX_ / tArrive_;
            double err = 100.0 * (cMe - cTh) / cTh;
            std::cout << "[FEM] bar-wave verification: measured c = " << cMe
                      << " m/s, theory sqrt(E/rho) = " << cTh
                      << " m/s, error = " << err << " %  "
                      << (std::abs(err) < 3.0 ? "[PASS]" : "[FAIL]") << "\n";
        } else {
            std::cout << "[FEM] bar-wave verification: wave did not reach the gauge [FAIL]\n";
        }
        return;
    }
    double Es = erodedVol_ > 0 ? work_ / erodedVol_ : 0.0;
    std::cout << "[FEM] peak tool force   : " << peakF_ << " N/m\n"
              << "[FEM] tool work output  : " << work_ << " J/m";
    if (tool_.motion == Tool::Motion::FREE)
        std::cout << "  (tool KE loss: " << toolKE0_ - tool_.ke() << " J/m)";
    std::cout << "\n"
              << "[FEM] eroded volume     : " << erodedVol_ << " m^3/m ("
              << nEroded_ << " elements)\n"
              << "[FEM] specific energy   : " << Es << " J/m^3\n";
}

} // namespace rockim
