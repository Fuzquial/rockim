#pragma once
// ---------------------------------------------------------------------------
// DemSolver: 2D bonded-particle model (BPM). Circular particles, linear
// contact law with Coulomb friction, parallel bonds carrying force AND
// moment between bonded neighbours, bond breakage by tensile and shear
// stress criteria. Fracture and fragmentation emerge from bond breakage.
// See DemSolver.cpp for the derivations (bond mechanics, breakage, dt).
// ---------------------------------------------------------------------------
#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>
#include <Eigen/Dense>

#include "rockim/Config.hpp"
#include "rockim/Material.hpp"
#include "rockim/Solver.hpp"
#include "rockim/Tool.hpp"

namespace rockim {

class DemSolver : public Solver {
public:
    DemSolver(const Config& cfg, std::string outDir);

    void init() override;
    void step() override;
    void writeFrame(int frame) override;
    void historyHeader(std::ostream&) const override;
    void historyRow(std::ostream&) const override;
    void finalize() override;

private:
    enum class Scenario { PERCUSSION, SHEAR, TENSION };
    enum : uint8_t { FREE = 0, FIXED = 1, PRESCRIBED = 2 };

    struct Part {
        Eigen::Vector2d x, v, f;
        double w = 0, tq = 0;        // angular velocity, torque (z)
        double r = 0, m = 0, I = 0;
        uint8_t flag = FREE;
    };

    struct Bond {
        int i = 0, j = 0;
        double Fn = 0;               // axial force, TENSION POSITIVE
        double Fs = 0;               // tangential force on particle j
        double Mb = 0;               // bending moment on particle j
        double A = 0, Ib = 0, Rb = 0;// cross-section, inertia, half-width
        double knA = 0, ksA = 0, knI = 0; // precomputed stiffnesses
        double sc = 0, tc0 = 0;      // tensile strength, cohesion (with scatter)
        bool broken = false;
        double tBreak = -1.0;
    };

    // --- setup -----------------------------------------------------------
    void buildPacking();
    void buildBonds();
    void placeTool();
    void computeStableDt();

    // --- per-step pieces -------------------------------------------------
    void bondForces();
    void contactForces();
    void wallAndToolForces();
    void integrate();
    void rebuildGrid();
    void toolParticleContact(Part& p);

    // --- fragments -------------------------------------------------------
    void computeFragments();          // union-find over intact bonds
    int  uf(std::vector<int>& parent, int a) const;

    // --- data ------------------------------------------------------------
    Config cfg_;
    std::string out_;
    Material mat_;
    Tool tool_;
    Scenario scen_ = Scenario::PERCUSSION;

    double W_ = 0.2, H_ = 0.1, thk_ = 1.0;
    double r_ = 1.25e-3;
    std::string packing_ = "hex";

    std::vector<Part> p_;
    std::vector<Bond> b_;
    long nBroken_ = 0;

    // frictional contact tangential springs (history), keyed by pair id
    std::unordered_map<uint64_t, double> tang_, tangNew_;

    // contact / bond micro-parameters
    double knC_ = 0, ksC_ = 0, mu_ = 0.5, xiC_ = 0.1;
    double kbn_ = 0;                  // bond modulus per length: E / L0 (per bond, via knA)
    double ksRatio_ = 0.4;
    double tanPhiB_ = 0;
    double lambda_ = 1.0;             // bond radius factor: Rb = lambda*min(ri,rj)
    double damping_ = 0.0;            // Cundall local non-viscous damping

    // walls
    bool bottomWall_ = true, sideWalls_ = false;

    // cell list
    double cell_ = 0;
    int gx_ = 0, gy_ = 0;
    Eigen::Vector2d gmin_;
    std::vector<int> head_, nxt_;

    // tension-test grips
    double pullV_ = 0.05;
    Eigen::Vector2d gripF_{0, 0};
    double sigmaPeak_ = 0;

    // fragments / diagnostics (refreshed at output frames)
    std::vector<int> fragId_;
    std::vector<int> nIntact_;           // intact bonds per particle (dashpot gate)
    std::vector<double> cAbsX_, cAbsY_;  // Lysmer dashpot coeffs per particle
    std::vector<double> kAbsX_, kAbsY_;  // boundary spring coeffs per particle
    std::vector<Eigen::Vector2d> xAnchor_;
    int nFrag_ = 1;
    double detachedVol_ = 0;
    double work_ = 0, peakF_ = 0, toolKE0_ = 0;
};

} // namespace rockim
