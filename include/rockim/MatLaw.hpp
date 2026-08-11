#pragma once
// ---------------------------------------------------------------------------
// MatLaw: pluggable constitutive laws for the 3D continuum FEM (mode fem3d).
// One integration point per tet: the solver hands the law the CO-ROTATED
// total (Biot) strain and the element size, the law returns the Cauchy
// stress in the same frame and updates its internal state. Sign convention:
// tension positive, p = tr(sigma)/3 (tension positive).
//
// Laws (config key `law = ...`):
//
//  * elastic  — linear elasticity (lambda, mu from E, nu). The control case.
//
//  * dpr      — rate-INdependent elastoplastic damage:
//               Drucker-Prager cone F = sqrt(J2) + alpha I1 - k (alpha, k
//               matched to Mohr-Coulomb in TRIAXIAL COMPRESSION from
//               cohesion/frictionDeg), deviatoric radial return (psi = 0);
//               Rankine tensile damage: when the max principal EFFECTIVE
//               stress exceeds ft, a scalar damage D grows with exponential
//               softening regularized by the crack band (Gf smeared over the
//               element size lc), nominal stress = (1-D) sigma_eff.
//               The 3D sibling of the 2D fem module's DP/Rankine + damage.
//
//  * saksala2011 — FAITHFUL PORT of the thesis' vumat_saksala_2011.f90
//               (Saksala, IJNAMG 35 (2011) 1483-1505): DP cone with
//               non-associated potential (skBetaDP), modified Rankine
//               (associated), parabolic cap (associated) joined to the
//               cone, viscoplastic CONSISTENCY law with two viscosities
//               (skSdp compression, skSmr tension), confinement-dependent
//               cohesion degradation (skNd, confinement frozen at the
//               elastic trial state), logarithmic volumetric cap hardening
//               (skDcap, skWcap), exponential tensile damage driven by the
//               positive norm of the viscoplastic strain increment (skAt,
//               skBetaT), unilateral crack closure, Koiter corner rule.
//               Deliberately NO erosion, NO fracture-energy regularization,
//               NO compressive damage — as in the published law and the
//               VUMAT. Verified by trace superposition against the Fortran
//               reference on the three loading paths of its test harness.
//
//  * dpdfh    — FAITHFUL PORT of the thesis' VUMATS/dfh/vumat_kstdfh.f
//               (DP-DFH Bohus, Shariati/Saadati/Hild arXiv:2201.01870;
//               obscuration: Forquin & Hild 2010): Drucker-Prager
//               plasticity on the EFFECTIVE stress (f = q - p tan(beta) -
//               d, non-associated psi, radial return + apex), gated to
//               COMPRESSION only (OPTION-3: in tension the principal
//               stress rises elastically to the DFH initiation threshold);
//               ANISOTROPIC local DFH damage: three directional damages
//               D_i in a frame FROZEN at first initiation (Euler ZYX),
//               initiation thresholds = three sorted Weibull draws
//               sig_k = sigw (Zeff/V_el)^(1/m), V_el = lc^3, drawn
//               DETERMINISTICALLY from a 64-bit spatial hash of the
//               element's initial centroid (micrometre resolution,
//               xorshift64 — bit-compatible with the VUMAT's kst_seed);
//               growth by the cube-root obscuration integrator
//               D_i = 1 - exp(-x^3), dx = (S lam)^(1/3) k c dt with
//               lam = (sigma_i/sigw)^m / Zeff floored at 1/V_el and
//               c = sqrt(E/rho); nominal stress: tensile principal
//               components scaled by (1-D_i), shears by min(f_i, f_j)
//               (closed crack = full shear transfer). NO regularization,
//               deletion OFF by default (dfhDeld = 1e9), exactly like the
//               VUMAT (percussion practice: damaged elements stay as the
//               rubble bed). ftScale multiplies the Weibull scale (the
//               solver's FIELD mechanism — sandbox extension, neutral at
//               ftScale = 1). Verified by trace superposition against the
//               ifx-compiled Fortran reference on four loading paths.
//
//  * saksala  — rate-DEPENDENT damage-viscoplasticity in the spirit of
//               Saksala's model for percussive drilling: the SAME DP cone
//               but with a PERZYNA viscoplastic return (linear overstress,
//               viscosity saksalaEta [Pa s]: consistency F = eta dlam/dt,
//               closed form dlam = F_tr / (G + eta/dt)), a compression CAP
//               on the pressure (crush pressure capP0, linear hardening
//               capH: volumetric plastic return with pc += H |deps_v^p|)
//               and the same crack-band Rankine tensile damage.
//               SIMPLIFICATIONS stated plainly: Perzyna overstress instead
//               of Saksala's consistency formulation, scalar tensile damage
//               only (no separate compressive damage variable), no
//               viscosity on the damage side, zero dilatancy. The rate
//               effect on strength is exact and verified: steady uniaxial
//               overstress = sqrt(3) eta epdot / (1/sqrt(3) - alpha).
//
// The crack band imposes the classical limit lc <= E Gf / ft^2 (snapback):
// the law throws at init if the mesh is too coarse for (Gf, ft).
// ---------------------------------------------------------------------------
#include <memory>
#include <string>

#include <Eigen/Dense>

#include "rockim/Config.hpp"
#include "rockim/Material.hpp"

namespace rockim {

struct MatState {
    Eigen::Matrix3d epsP = Eigen::Matrix3d::Zero();    // plastic strain
    double D = 0.0;                                    // tensile damage
    double kappa = 0.0;                                // damage driving strain
    double epvEq = 0.0;                                // equiv. viscopl. strain
    double pc = 0.0;                                   // current cap pressure
    bool eroded = false;

    // per-element strength factor set by the SOLVER before the first law
    // call (Weibull heterogeneity — the sandbox version of the VUMAT's
    // FIELD mechanism); 1 = homogeneous
    double ftScale = 1.0;

    // initial element centroid, set by the solver after meshing (dpdfh
    // seeds its Weibull draws from a spatial hash of these coordinates,
    // exactly like the VUMAT hashes coordMp)
    Eigen::Vector3d x0 = Eigen::Vector3d::Zero();

    // dpdfh sub-state (mirrors the VUMAT's 16-SDV layout; Voigt-6 order
    // 11 22 33 12 23 13, tensorial shears)
    struct Dfh {
        bool seeded = false;
        double snom[6] = {0, 0, 0, 0, 0, 0};   // last NOMINAL stress
        double epsPrev[6] = {0, 0, 0, 0, 0, 0};
        double Dv[3] = {0, 0, 0};              // SDV 4-6
        double eul[3] = {0, 0, 0};             // SDV 7-9 (Euler ZYX, rad)
        double sc[3] = {0, 0, 0};              // SDV 10-12 (sorted draws)
        double ti[3] = {0, 0, 0};              // SDV 13-15 (init times)
        double peeq = 0.0;                     // SDV 3
        double smaxh = 0.0;                    // SDV 16
        double t = 0.0;                        // accumulated total time
        bool dead = false;                     // SDV 1 (STATUS)
    } dfh;

    // saksala2011 sub-state (Voigt-6, order 11 22 33 12 23 13, tensorial
    // shears — mirrors the VUMAT's SDV layout, incl. the local strengths
    // SDV15/16)
    struct Sk11 {
        bool init = false;
        double sbar[6] = {0, 0, 0, 0, 0, 0};           // effective stress
        double epsPrev[6] = {0, 0, 0, 0, 0, 0};        // last total strain
        double kapDP = 0.0, kapMR = 0.0, eqvt = 0.0, epsv = 0.0, pp = 0.0;
        double ftLoc = 0.0, c0Loc = 0.0;               // local strengths
        int active = 0, failed = 0;
    } sk;
};

class MatLaw {
public:
    virtual ~MatLaw() = default;
    // eps: co-rotated Biot strain; lc: element size (V^(1/3)) for the crack
    // band; dt: time step (rate laws). Returns nominal Cauchy stress.
    virtual Eigen::Matrix3d stress(const Eigen::Matrix3d& eps, MatState& s,
                                   double dt, double lc) const = 0;
    virtual std::string name() const = 0;
    double cP() const { return mat_.cP(); }

    // uniaxial compressive yield of the DP cone (analytic, for verification)
    double sigmaCdp() const { return kdp_ / (1.0 / std::sqrt(3.0) - adp_); }
    // steady uniaxial viscous overstress at axial strain rate epdot
    // (0 for rate-independent laws)
    virtual double viscousOverstress(double) const { return 0.0; }

    static std::unique_ptr<MatLaw> make(const std::string& kind,
                                        const Material& m, const Config& c,
                                        double lcMax);

protected:
    explicit MatLaw(const Material& m) : mat_(m) {
        lam_ = m.E * m.nu / ((1.0 + m.nu) * (1.0 - 2.0 * m.nu));
        G_ = m.G();
        K_ = m.K();
        double sf = std::sin(m.phiDeg * M_PI / 180.0);
        double cf = std::cos(m.phiDeg * M_PI / 180.0);
        adp_ = 2.0 * sf / (std::sqrt(3.0) * (3.0 - sf));
        kdp_ = 6.0 * mat_.cohesion * cf / (std::sqrt(3.0) * (3.0 - sf));
    }
    Eigen::Matrix3d elastic(const Eigen::Matrix3d& eps) const {
        return lam_ * eps.trace() * Eigen::Matrix3d::Identity()
               + 2.0 * G_ * eps;
    }
    Material mat_;
    double lam_ = 0.0, G_ = 0.0, K_ = 0.0;
    double adp_ = 0.0, kdp_ = 0.0;                     // DP cone (MC triax)
};

// Replays the three material-point loading paths of the Fortran reference
// harness (VUMATS/saksala/test_saksala_2011.f90) through the C++ port and
// writes the same CSV trace, for quantitative superposition against the
// ifx-compiled reference. Returns 0.
int saksala2011Selftest(const std::string& csvPath);

// Same idea for the DP-DFH port: replays the four material-point paths of
// VUMATS/dfh/test_kstdfh.f90 (two tension rates, deviatoric compression,
// oblique tension with shears) and writes the same CSV (stresses in MPa).
int dpdfhSelftest(const std::string& csvPath);

} // namespace rockim
