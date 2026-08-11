#pragma once
// ---------------------------------------------------------------------------
// RandomField: stationary Gaussian random field g(x) ~ N(0,1) on a rectangle
// with a Gaussian correlation kernel of length ell, built as kernel-smoothed
// lattice noise (moving-average construction):
//
//   g(x) = sum_k w_k(x) eta_k / sqrt( sum_k w_k(x)^2 ),
//   w_k(x) = exp( -|x - x_k|^2 / (2 ell^2) ),   eta_k iid N(0,1)
//
// on a noise lattice of spacing ell/2 extended 3 ell beyond the domain (so
// boundary points see a full kernel). The normalization makes the marginal
// EXACTLY N(0,1) at every point; the correlation between two points falls
// off over ~ell (Gaussian-kernel moving average: correlation length of the
// field is ell up to a factor sqrt(2), which is the right order — the GBM
// use case sets ell to the grain size, not to a precise covariance model).
//
// Deterministic for a given (domain, ell, seed) and query point — the field
// exists independently of any mesh, which is the whole point: two DIFFERENT
// meshes sampling the SAME field see the same weak zones. This is the
// sandbox version of the thesis' correlated-sigma_w idea (inject the defect
// statistics as a spatial field instead of an independent draw per element).
// ---------------------------------------------------------------------------
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>

namespace rockim {

// The correlation may be ANISOTROPIC and ROTATED: lengths (ellA, ellB) along
// axes rotated by angleDeg. ellB << ellA gives band-like weak zones along
// the ellA direction — the sandbox representation of a foliation or any
// oriented texture. Isotropic = ellB == ellA (angle irrelevant).
class RandomField {
public:
    RandomField(double W, double H, double ellA, double ellB, double angleDeg,
                unsigned seed)
        : ea_(ellA), eb_(ellB) {
        if (!(ellA > 0.0) || !(ellB > 0.0))
            throw std::runtime_error("RandomField: correlation lengths must "
                                     "be > 0");
        ca_ = std::cos(angleDeg * M_PI / 180.0);
        sa_ = std::sin(angleDeg * M_PI / 180.0);
        double emin = std::min(ellA, ellB), emax = std::max(ellA, ellB);
        h_ = 0.5 * emin;
        R_ = (int)std::ceil(3.0 * emax / h_);
        double margin = 3.0 * emax;
        x0_ = -margin;
        y0_ = -margin;
        nx_ = (int)std::ceil((W + 2.0 * margin) / h_) + 2;
        ny_ = (int)std::ceil((H + 2.0 * margin) / h_) + 2;
        std::mt19937 rng(seed);
        std::normal_distribution<double> N(0.0, 1.0);
        noise_.resize((std::size_t)nx_ * ny_);
        for (double& v : noise_) v = N(rng);
    }

    RandomField(double W, double H, double ell, unsigned seed)
        : RandomField(W, H, ell, ell, 0.0, seed) {}

    double operator()(const Eigen::Vector2d& p) const {
        int ci = (int)std::floor((p.x() - x0_) / h_);
        int cj = (int)std::floor((p.y() - y0_) / h_);
        double s = 0.0, s2 = 0.0;
        for (int j = cj - R_; j <= cj + R_; ++j) {
            if (j < 0 || j >= ny_) continue;
            for (int i = ci - R_; i <= ci + R_; ++i) {
                if (i < 0 || i >= nx_) continue;
                double dx = p.x() - (x0_ + i * h_);
                double dy = p.y() - (y0_ + j * h_);
                double da = ( ca_ * dx + sa_ * dy) / ea_;   // along ellA
                double db = (-sa_ * dx + ca_ * dy) / eb_;   // across
                double w = std::exp(-0.5 * (da * da + db * db));
                s += w * noise_[(std::size_t)j * nx_ + i];
                s2 += w * w;
            }
        }
        return s2 > 0.0 ? s / std::sqrt(s2) : 0.0;
    }

private:
    double ea_, eb_, ca_, sa_, h_, x0_, y0_;
    int nx_, ny_, R_;
    std::vector<double> noise_;
};

// The 3D lift of RandomField, same moving-average construction on a 3D
// noise lattice — the marginal is exactly N(0,1) at every point and the
// field exists independently of any mesh (same fieldSeed = same weak zones
// for two different meshes), which is the property the correlated-sigma_w
// experiments need.
//
// Anisotropy models a FOLIATION: two in-plane correlation lengths ellA
// (isotropic within the texture plane) and one length ellB ACROSS it, the
// plane being tilted by angleDeg about the y-axis (rotation in the x-z
// plane, the 3D reading of the 2D strengthCorrAngleDeg). angle = 0 keeps
// the texture plane horizontal: ellB << ellA then gives weak BANDS normal
// to z — the 3D analogue of the 2D banded field. Isotropic = ellB == ellA.
class RandomField3 {
public:
    RandomField3(double W, double D, double H, double ellA, double ellB,
                 double angleDeg, unsigned seed)
        : ea_(ellA), eb_(ellB) {
        if (!(ellA > 0.0) || !(ellB > 0.0))
            throw std::runtime_error("RandomField3: correlation lengths must "
                                     "be > 0");
        ca_ = std::cos(angleDeg * M_PI / 180.0);
        sa_ = std::sin(angleDeg * M_PI / 180.0);
        double emin = std::min(ellA, ellB), emax = std::max(ellA, ellB);
        h_ = 0.5 * emin;
        R_ = (int)std::ceil(3.0 * emax / h_);
        double margin = 3.0 * emax;
        x0_ = -margin; y0_ = -margin; z0_ = -margin;
        nx_ = (int)std::ceil((W + 2.0 * margin) / h_) + 2;
        ny_ = (int)std::ceil((D + 2.0 * margin) / h_) + 2;
        nz_ = (int)std::ceil((H + 2.0 * margin) / h_) + 2;
        // the noise lattice is O((L/ell)^3): a correlation length far below
        // the domain size exhausts memory silently — fail with the fix
        if ((double)nx_ * ny_ * nz_ > 2.68e8)
            throw std::runtime_error("RandomField3: noise lattice exceeds "
                                     "2.7e8 points — increase "
                                     "strengthCorrLength (or shrink the "
                                     "domain)");
        std::mt19937 rng(seed);
        std::normal_distribution<double> N(0.0, 1.0);
        noise_.resize((std::size_t)nx_ * ny_ * nz_);
        for (double& v : noise_) v = N(rng);
    }

    RandomField3(double W, double D, double H, double ell, unsigned seed)
        : RandomField3(W, D, H, ell, ell, 0.0, seed) {}

    double operator()(const Eigen::Vector3d& p) const {
        int ci = (int)std::floor((p.x() - x0_) / h_);
        int cj = (int)std::floor((p.y() - y0_) / h_);
        int ck = (int)std::floor((p.z() - z0_) / h_);
        double s = 0.0, s2 = 0.0;
        for (int k = ck - R_; k <= ck + R_; ++k) {
            if (k < 0 || k >= nz_) continue;
            for (int j = cj - R_; j <= cj + R_; ++j) {
                if (j < 0 || j >= ny_) continue;
                for (int i = ci - R_; i <= ci + R_; ++i) {
                    if (i < 0 || i >= nx_) continue;
                    double dx = p.x() - (x0_ + i * h_);
                    double dy = p.y() - (y0_ + j * h_);
                    double dz = p.z() - (z0_ + k * h_);
                    // texture frame: e1 = x rotated in the x-z plane,
                    // e2 = y (in-plane), en = tilted z (across the plane)
                    double d1 = ( ca_ * dx + sa_ * dz) / ea_;
                    double d2 = dy / ea_;
                    double dn = (-sa_ * dx + ca_ * dz) / eb_;
                    double w = std::exp(-0.5 * (d1 * d1 + d2 * d2 + dn * dn));
                    s += w * noise_[((std::size_t)k * ny_ + j) * nx_ + i];
                    s2 += w * w;
                }
            }
        }
        return s2 > 0.0 ? s / std::sqrt(s2) : 0.0;
    }

private:
    double ea_, eb_, ca_, sa_, h_, x0_, y0_, z0_;
    int nx_, ny_, nz_, R_;
    std::vector<double> noise_;
};

} // namespace rockim
