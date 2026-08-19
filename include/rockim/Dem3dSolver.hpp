#pragma once
// ---------------------------------------------------------------------------
// Dem3dSolver — 3D bonded-particle model (spheres + parallel bonds).
//
// Direct extension of the 2D DemSolver to three dimensions:
//   * equal spheres on an HCP lattice (12 bonded neighbours) for impact /
//     cutting, or a simple-cubic lattice (6 neighbours, load-aligned) for the
//     tension verification;
//   * PFC3D-style parallel bonds carrying a normal force, a shear-force
//     vector, a twisting moment and a bending-moment vector; brittle failure
//     by maximum tensile fibre stress or shear stress (with frictional
//     strengthening under compression);
//   * broken pairs fall back to frictional contact (linear normal spring +
//     vector tangential spring with history, Coulomb cap);
//   * rigid spherical tool (free percussive impactor or prescribed cutter);
//   * Lysmer + Deeks-Randolph viscous-spring quiet boundaries on the four
//     lateral faces and the bottom, gated on bondedness;
//   * fragments from connected components of the intact bond network.
//
// Scope notes (documented in the README): regular lattice (calibration and a
// disordered packing are required for quantitative work), and the stored bond
// shear/bending vectors are kept perpendicular to the bond axis by projection
// each step — exact for the small per-step rotations of an explicit run.
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

namespace rockim {

class Dem3dSolver : public Solver {
public:
    Dem3dSolver(const Config& cfg, std::string outDir);

    void init() override;
    void step() override;
    void writeFrame(int frame) override;
    void historyHeader(std::ostream&) const override;
    void historyRow(std::ostream&) const override;
    void finalize() override;

private:
    enum Flag { FREE = 0, FIXED = 1, PRESCRIBED = 2 };
    enum class Scenario { PERCUSSION, SHEAR, TENSION };

    struct Part {
        Eigen::Vector3d x, v, w, f, tq;
        double r, m, I;
        int flag;
    };

    struct Bond {
        int i, j;
        double Rb, A, Ib, Jb;         // radius, area, bending I, polar J
        double knA, ksA, knI, ksJ;    // stiffnesses (force / moment per unit)
        double sc, tc0;               // tensile strength, cohesion
        double Fn = 0.0;              // normal force (tension positive)
        Eigen::Vector3d Fs{0, 0, 0};  // shear force (perp. to bond axis)
        double Mt = 0.0;              // twisting moment (about bond axis)
        Eigen::Vector3d Mb{0, 0, 0};  // bending moment (perp. to bond axis)
        bool broken = false;
        double tBreak = -1.0;
    };

    struct Tool3 {
        bool free = true;
        bool flat = false;   // flat-ended cylindrical punch (axis z); x is
                             // then the center of the bottom face — the 3D
                             // lift of the 2D FLAT tool (percussion only)
        double mass = 0.5, radius = 0.015;
        Eigen::Vector3d x{0, 0, 0}, v{0, 0, 0}, F{0, 0, 0};
        void integrate(double dt) { if (free) v += (dt / mass) * F; x += dt * v; }
        double ke() const { return 0.5 * mass * v.squaredNorm(); }
    };

    // setup
    void buildPacking();
    void buildBonds();
    void placeTool();
    void setupAbsorbing();
    void computeStableDt();

    // stepping
    void bondForces();
    void rebuildGrid();
    void contactForces();
    void toolAndWallForces();
    void integrate();
    void computeFragments();

    // E4 (2026-08-19) : la cle est desormais TRIEE. Sans tri, (i,j) et (j,i)
    // donnent deux cles differentes pour la MEME paire. La boucle en cellules
    // liees ne garantit l'ordre i<j que dans la cellule d'origine (garde
    // `home && j <= i`) : des qu'une particule change de cellule, les roles
    // s'inversent, la cle change, et l'historique du ressort tangentiel est
    // perdu EN SILENCE — le contact repart d'un glissement nul. Le tri rend la
    // cle invariante par permutation.
    // ⚠️ Ce correctif n'est PAS bit-neutre : le repere dem3d_tension de la
    // suite doit etre re-mesure. Il l'est a raison — l'ancienne valeur etait
    // celle d'un historique tangentiel amnesique.
    static uint64_t pairKey(int i, int j) {
        int a = (i < j) ? i : j, b = (i < j) ? j : i;
        return (uint64_t(uint32_t(a)) << 32) | uint32_t(b);
    }

    Config cfg_;
    std::string out_;
    Material mat_;
    Scenario scen_ = Scenario::PERCUSSION;

    double W_ = 0.1, D_ = 0.1, H_ = 0.08;   // x, y, z extents
    double r_ = 1.5e-3;
    std::string packing_ = "hcp";

    double knC_ = 0.0, ksC_ = 0.0, ksRatio_ = 0.4;
    double mu_ = 0.5, xiC_ = 0.05;
    double tanPhiB_ = 0.8, lambda_ = 1.0, damping_ = 0.02;
    bool bottomWall_ = true;
    double pullV_ = 0.05;
    double pullRamp_ = 0.0;               // grip velocity rise time [s]
    bool gripFree_ = false;               // frictionless tension grips

    std::vector<Part> p_;
    std::vector<Bond> b_;
    std::vector<int> nIntact_;
    std::vector<int> fragId_;

    // quiet boundaries
    std::vector<Eigen::Vector3d> cAbs_, kAbs_, xAnchor_;

    // contact history (particle-particle and particle-tool tangential springs)
    std::unordered_map<uint64_t, Eigen::Vector3d> tang_, tangNew_;
    std::unordered_map<int, Eigen::Vector3d> toolTang_, toolTangNew_;

    // cell grid
    double cell_ = 0.0;
    Eigen::Vector3d gmin_;
    int gx_ = 1, gy_ = 1, gz_ = 1;
    std::vector<int> head_, nxt_;

    Tool3 tool_;
    double toolKE0_ = 0.0;

    // results
    double work_ = 0.0, peakF_ = 0.0;
    long nBroken_ = 0;
    int nFrag_ = 1;
    double detachedVol_ = 0.0;
    Eigen::Vector3d gripF_{0, 0, 0};
    double sigmaPeak_ = 0.0;
};

} // namespace rockim
