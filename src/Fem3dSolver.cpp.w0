// ---------------------------------------------------------------------------
// Fem3dSolver — 3D continuum FEM with pluggable laws. Mesh/boundaries/tool
// mirror the (verified) Fdem3dSolver; the bulk physics is delegated to
// MatLaw and fracture is smeared damage + erosion, as in the 2D fem module.
// ---------------------------------------------------------------------------
#include "rockim/Fem3dSolver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
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

    buildMesh();
    law_ = MatLaw::make(cfg_.gets("law", "dpr"), mat_, cfg_, lcMax_);

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
    muC_ = cfg_.getd("contactMu", 0.5);
    xiC_ = cfg_.getd("contactXi", 0.05);
    vReg_ = cfg_.getd("contactVreg", 1e-3);

    placeTool();
    setupBoundaries();
    computeStableDt();

    if (scen_ == Scenario::TENSION) pullV_ = cfg_.getd("pullV", 0.05);
    pullRamp_ = cfg_.getd("pullRamp", 0.0);
    gripFree_ = cfg_.getb("gripLateralFree", false);
    toolKE0_ = tool_.ke();

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
        if (e.V0 <= 0) throw std::runtime_error("degenerate tet");
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
    if (sh != "sphere" && sh != "disc" && sh != "flat")
        throw std::runtime_error("toolShape must be sphere | flat (3D)");
    tool_.flat = sh == "flat";
    if (scen_ == Scenario::PERCUSSION) {
        tool_.free = true;
        double zTip = tool_.flat ? H_ + gap : H_ + tool_.radius + gap;
        tool_.x = {cfg_.getd("toolX", 0.5 * W_), cfg_.getd("toolY", 0.5 * D_),
                   zTip};
        tool_.v = {0.0, 0.0, -cfg_.getd("impactSpeed", 8.0)};
    } else {                     // SHEAR: dragged spherical cutter, as fdem3d
        tool_.free = false;
        tool_.flat = false;
        double depth = cfg_.getd("cutDepth", 0.004);
        tool_.x = {cfg_.getd("toolX", -tool_.radius - gap), 0.5 * D_,
                   H_ - depth + tool_.radius};
        tool_.v = {cfg_.getd("cutSpeed", 10.0), 0.0, 0.0};
    }
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
                } else {                       // percussion AND shear: the
                    flag_[nid] = FIXED;        // block needs its support
                }
            }
        }
    }
}

void Fem3dSolver::computeStableDt() {
    double cfl = hmin_ / mat_.cP();
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

    if (scen_ == Scenario::TENSION) {
        gripF_.setZero();
        for (int i = 0; i < (int)X0_.size(); ++i)
            if (flag_[i] == PRESCRIBED) gripF_ += f_[i];
        sigmaPeak_ = std::max(sigmaPeak_,
                              std::abs(gripF_.z()) / (W_ * D_));
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
    }

    integrate();
    t_ += dt_;
    if ((++stepCount_ & 1023) == 0)
        if (!std::isfinite(work_ + sigmaPeak_) || !std::isfinite(u_[0].x()))
            throw std::runtime_error("FEM3D instability (NaN) — reduce "
                                     "dtFactor");
}

// Co-rotational tets; the law does the physics. Shared nodes make the
// scatter race under OpenMP: per-thread buffers reduced in thread order
// (deterministic per thread count; 1 thread = serial arithmetic).
void Fem3dSolver::elementForces() {
    long nEro = 0;
    auto processElem = [&](Elem& e, auto&& addF, long& ero) {
        if (e.st.eroded) return;
        Eigen::Matrix3d F = Eigen::Matrix3d::Zero();
        for (int a = 0; a < 4; ++a)
            F += (X0_[e.n[a]] + u_[e.n[a]]) * e.dN.col(a).transpose();
        Eigen::Matrix3d R;
        double det = F.determinant();
        if (det > 1e-9) {
            R = F * std::sqrt(3.0) / F.norm();
            for (int it = 0; it < 3; ++it)
                R = 0.5 * (R + R.inverse().transpose());
        } else R.setIdentity();
        Eigen::Matrix3d Ub = R.transpose() * F;
        Eigen::Matrix3d eps = 0.5 * (Ub + Ub.transpose())
                              - Eigen::Matrix3d::Identity();
        Eigen::Matrix3d sig = law_->stress(eps, e.st, dt_, e.lc);
        if (e.st.eroded) { ++ero; e.svm = 0.0; e.pm = 0.0; e.szz = 0.0; return; }
        double pm = sig.trace() / 3.0;
        e.pm = pm;
        e.szz = sig(2, 2);
        e.svm = std::sqrt(1.5)
                * (sig - pm * Eigen::Matrix3d::Identity()).norm();
        Eigen::Matrix3d P = R * sig;
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
        std::vector<long> eroT(nT, 0);
#pragma omp parallel
        {
            int t = omp_get_thread_num();
            auto& fb = fTL_[t];
            auto& seen = seenTL_[t];
            auto& tl = touchedTL_[t];
            tl.clear();
            long ero = 0;
            auto addF = [&](int i, const Eigen::Vector3d& v3) {
                if (!seen[i]) { seen[i] = 1; tl.push_back(i); }
                fb[i] += v3;
            };
#pragma omp for schedule(static)
            for (int eI = 0; eI < (int)el_.size(); ++eI)
                processElem(el_[eI], addF, ero);
            eroT[t] = ero;
        }
        for (int t = 0; t < nT; ++t) {
            for (int i : touchedTL_[t]) {
                f_[i] += fTL_[t][i];
                fTL_[t][i].setZero();
                seenTL_[t][i] = 0;
            }
            nEro += eroT[t];
        }
        nEroded_ += nEro;
        return;
    }
#endif
    auto addF = [&](int i, const Eigen::Vector3d& v3) { f_[i] += v3; };
    for (auto& e : el_) processElem(e, addF, nEro);
    nEroded_ += nEro;
}

void Fem3dSolver::toolContact() {
    if (scen_ == Scenario::TENSION) return;
    for (int i = 0; i < (int)X0_.size(); ++i) {
        Eigen::Vector3d p = X0_[i] + u_[i];
        Eigen::Vector3d Fc;
        if (tool_.flat) {
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
        if (flag_[i] == PRESCRIBED) {
            if (gripFree_) {
                v_[i].x() += (dt_ / m_[i]) * f_[i].x();
                v_[i].y() += (dt_ / m_[i]) * f_[i].y();
            } else { v_[i].x() = 0.0; v_[i].y() = 0.0; }
            double vg = pullV_;
            if (pullRamp_ > 0.0 && t_ < pullRamp_)
                vg *= 0.5 * (1.0 - std::cos(M_PI * t_ / pullRamp_));
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
        u_[i] += dt_ * v_[i];
    }
    if (scen_ != Scenario::TENSION) tool_.integrate(dt_);
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

void Fem3dSolver::historyHeader(std::ostream& os) const {
    if (scen_ == Scenario::TENSION) {
        os << "t,sigma,sigmaPeak,nEroded\n";
        return;
    }
    // toolFx/toolX appended at the END so percussion post-processing that
    // reads columns by position keeps working; they are the cutting force
    // and advance of the shear scenario
    os << "t,toolFz,toolZ,toolVz,work,toolKE,nEroded,craterVol,"
          "toolFx,toolX\n";
}

void Fem3dSolver::historyRow(std::ostream& os) const {
    if (scen_ == Scenario::TENSION) {
        os << t_ << "," << std::abs(gripF_.z()) / (W_ * D_) << ","
           << sigmaPeak_ << "," << nEroded_ << "\n";
        return;
    }
    os << t_ << "," << tool_.F.z() << "," << tool_.x.z() << ","
       << tool_.v.z() << "," << work_ << "," << tool_.ke() << ","
       << nEroded_ << "," << craterVol() << ","
       << tool_.F.x() << "," << tool_.x.x() << "\n";
}

void Fem3dSolver::finalize() {
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
    std::cout << "[FEM3D] peak tool force : " << peakF_ << " N\n"
              << "[FEM3D] tool work       : " << work_ << " J  (tool KE loss: "
              << toolKE0_ - tool_.ke() << " J)\n";
    double cv = craterVol();
    if (cv > 0)
        std::cout << "[FEM3D] specific energy : " << work_ / cv << " J/m^3\n";
}

} // namespace rockim
