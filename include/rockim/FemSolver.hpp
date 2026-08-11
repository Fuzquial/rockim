#pragma once
// ---------------------------------------------------------------------------
// FemSolver: explicit-dynamics FEM in plane strain with constant-strain
// triangles (CST), lumped mass, an isotropic damage model built on a
// Drucker-Prager shear surface + Rankine tension cutoff, and element erosion.
// See FemSolver.cpp for the full derivations (time step, damage, contact).
// ---------------------------------------------------------------------------
#include <array>
#include <cstdint>
#include <string>
#include <vector>
#include <Eigen/Dense>

#include "rockim/Config.hpp"
#include "rockim/Material.hpp"
#include "rockim/Solver.hpp"
#include "rockim/Tool.hpp"

namespace rockim {

class FemSolver : public Solver {
public:
    FemSolver(const Config& cfg, std::string outDir);

    void init() override;
    void step() override;
    void writeFrame(int frame) override;
    void historyHeader(std::ostream&) const override;
    void historyRow(std::ostream&) const override;
    void finalize() override;

private:
    enum class Scenario { PERCUSSION, SHEAR, BAR_WAVE };

    struct Elem {
        std::array<int, 3> n{};              // node indices (CCW)
        double A = 0;                        // reference area
        double lc = 0;                       // characteristic length (crack band)
        double hMin = 0;                     // smallest altitude (time-step length)
        Eigen::Matrix<double, 3, 6> B;       // constant strain-displacement matrix
        double kappaT = 0, kappaS = 0;       // damage history variables
        double D = 0;                        // scalar damage
        bool eroded = false;
        double svm = 0, smean = 0;           // stored for output (nominal stress)
    };

    // --- setup -----------------------------------------------------------
    void buildMesh();
    void placeTool();
    void computeStableDt();

    // --- per-step pieces -------------------------------------------------
    void applyContact();
    void internalForcesAndDamage();
    void integrate();
    void updateDamage(Elem& e, const Eigen::Vector3d& sEff, double szz);
    void erode(Elem& e);
    void refreshActiveNodes();

    // --- data ------------------------------------------------------------
    Config cfg_;
    std::string out_;
    Material mat_;
    Tool tool_;
    Scenario scen_ = Scenario::PERCUSSION;

    double W_ = 0.2, H_ = 0.1, thk_ = 1.0;
    int nx_ = 96, ny_ = 48;

    std::vector<Eigen::Vector2d> X0_, u_, v_, f_;
    std::vector<double> m_;
    std::vector<double> cAbsX_, cAbsY_;  // Lysmer dashpot coeffs per node [N s/m]
    std::vector<double> kAbsX_, kAbsY_;  // boundary spring coeffs per node [N/m]
    std::vector<uint8_t> fix_;       // bit0: ux fixed, bit1: uy fixed
    std::vector<uint8_t> active_;    // attached to >= 1 intact element
    std::vector<Elem> el_;
    bool activeDirty_ = false;

    // material-derived constants
    Eigen::Matrix3d Dm_;
    double dpAlpha_ = 0, dpK_ = 0;
    double kappa0T_ = 0, kappa0S_ = 0;

    // damage / erosion controls
    bool damageOn_ = true;
    double erodeD_ = 0.98;
    double strainCap_ = 0.15;

    // contact (penalty)
    double kp_ = 0, xiC_ = 0.05, muC_ = 0.3, vReg_ = 1e-3;

    // diagnostics
    double work_ = 0;               // energy delivered by the tool (tool-side)
    double toolKE0_ = 0;
    double erodedVol_ = 0;
    long   nEroded_ = 0;
    double peakF_ = 0;

    // bar-wave verification
    std::vector<int> barBcNodes_, gaugeNodes_;
    double barV0_ = 1.0, gaugeX_ = 0, tArrive_ = -1.0;
};

} // namespace rockim
