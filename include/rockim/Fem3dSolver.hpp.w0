#pragma once
// ---------------------------------------------------------------------------
// Fem3dSolver — 3D explicit continuum FEM on Kuhn tetrahedra with PLUGGABLE
// constitutive laws (MatLaw: elastic | dpr | saksala). The 3D sibling of the
// 2D fem module and the sandbox counterpart of an Abaqus/Explicit + VUMAT
// percussion run: shared-node mesh (no cohesive joints — fracture is smeared
// damage + element EROSION), co-rotational linear tets, rigid spherical
// tool, Lysmer quiet boundaries, percussion and tension/compression
// scenarios (tension with pullV < 0 IS the uniaxial compression test, with
// pullRamp and gripLateralFree available as in the 2D FDEM).
//
// Verification hooks: the DP cone's uniaxial compressive strength and the
// Perzyna viscous overstress have closed forms (MatLaw::sigmaCdp, README) —
// the tension/compression scenario prints measured vs expected.
// ---------------------------------------------------------------------------
#include <array>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "rockim/Config.hpp"
#include "rockim/MatLaw.hpp"
#include "rockim/Material.hpp"
#include "rockim/Solver.hpp"

namespace rockim {

class Fem3dSolver : public Solver {
public:
    Fem3dSolver(const Config& cfg, std::string outDir);

    void init() override;
    void step() override;
    void writeFrame(int frame) override;
    void historyHeader(std::ostream&) const override;
    void historyRow(std::ostream&) const override;
    void finalize() override;

private:
    enum Flag { FREE = 0, FIXED = 1, PRESCRIBED = 2 };
    enum class Scenario { PERCUSSION, SHEAR, TENSION };

    struct Elem {
        std::array<int, 4> n;                  // SHARED node ids
        Eigen::Matrix<double, 3, 4> dN;
        double V0 = 0.0, lc = 0.0;             // volume, size V0^(1/3)
        double svm = 0.0, pm = 0.0, szz = 0.0; // outputs
        MatState st;
    };
    struct BFace { int elem; std::array<int, 3> n; };
    struct Tool3 {
        double mass = 0.5, radius = 0.015;
        bool free = true;
        bool flat = false;   // flat-ended cylindrical punch (axis z); x is
                             // then the CENTER OF THE BOTTOM FACE, radius
                             // the punch radius — the 3D lift of the 2D
                             // FLAT tool (percussion only)
        Eigen::Vector3d x{0, 0, 0}, v{0, 0, 0}, F{0, 0, 0};
        void integrate(double dt) {
            if (free) v += (dt / mass) * F;
            x += dt * v;
        }
        double ke() const { return 0.5 * mass * v.squaredNorm(); }
    };

    void buildMesh();
    void placeTool();
    void setupBoundaries();
    void computeStableDt();
    void elementForces();
    void toolContact();
    void integrate();
    double craterVol() const;

    Config cfg_;
    std::string out_;
    Material mat_;
    std::unique_ptr<MatLaw> law_;
    Scenario scen_ = Scenario::PERCUSSION;

    double W_ = 0.1, D_ = 0.1, H_ = 0.08;
    int nx_ = 24, ny_ = 24, nz_ = 18;
    double hmin_ = 1e-3, lcMax_ = 1e-3;

    std::vector<Eigen::Vector3d> X0_, u_, v_, f_;
    std::vector<double> m_;
    std::vector<int> flag_;
    std::vector<Elem> el_;
    std::vector<BFace> exterior_;
    std::vector<Eigen::Vector3d> cAbs_, kAbs_;

    Tool3 tool_;
    double toolKE0_ = 0.0;
    double kp_ = 0.0, muC_ = 0.5, xiC_ = 0.05, vReg_ = 1e-3;
    double damping_ = 0.05, pullV_ = 0.05, pullRamp_ = 0.0;
    bool gripFree_ = false;

    long stepCount_ = 0, nEroded_ = 0;
    double work_ = 0.0, peakF_ = 0.0, sigmaPeak_ = 0.0;
    Eigen::Vector3d gripF_{0, 0, 0};
    // mid-third gauge: the grip force carries the Cundall-damping drag of
    // the flowing column (measured +11 % at damping 0.7), the mid-specimen
    // stress does not — verifications read here
    std::vector<int> midEl_;
    double sigMid_ = 0.0, sigMidPeak_ = 0.0;
    double sigMidSum_ = 0.0;               // plateau average (last quarter):
    long sigMidN_ = 0;                     // the end-state snapshot of a
                                           // ringing signal is mesh-fragile

    // OpenMP scratch (shared-node scatter)
    std::vector<std::vector<Eigen::Vector3d>> fTL_;
    std::vector<std::vector<char>> seenTL_;
    std::vector<std::vector<int>> touchedTL_;
};

} // namespace rockim
