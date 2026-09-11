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
#include "rockim/Tessellation.hpp"
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
    enum class Scenario { PERCUSSION, SHEAR, BAR_WAVE, TENSION };

    struct Elem {
        std::array<int, 3> n{};              // node indices (CCW)
        double A = 0;                        // reference area
        double lc = 0;                       // characteristic length (crack band)
        double hMin = 0;                     // smallest altitude (time-step length)
        // GBM 2D (2026-09-07) : grain et phase du triangle. En mesh = grid
        // (defaut) les deux valent 0 et phases_ n'a qu'une fiche = mat_, donc
        // tous les chemins indexes ci-dessous rendent exactement les anciennes
        // constantes globales.
        int phase = 0;
        int grain = 0;
        // B6 (2026-09-06) : initialiseur par defaut. Sans lui, les quatre
        // el_.push_back({{a, b, c}}) de FemSolver.cpp laissaient B non
        // initialise et clang levait « missing field 'B' initializer ».
        Eigen::Matrix<double, 3, 6> B =
            Eigen::Matrix<double, 3, 6>::Zero();  // constant strain-displacement matrix
        double kappaT = 0, kappaS = 0;       // damage history variables
        double D = 0;                        // scalar damage
        bool eroded = false;
        double svm = 0, smean = 0;           // stored for output (nominal stress)
        double syy = 0;                      // contrainte axiale (jauge de traction)
    };

    // --- setup -----------------------------------------------------------
    void buildMesh();                    // grille structuree (historique)
    void buildMeshVoronoi();             // GBM 2D (2026-09-07), opt-in
    void finishMesh();                   // geometrie / B / masses / gardes
    void applyBoundaryConditions();      // appuis, mors, frontieres absorbantes
    void phaseKeyGuards() const;         // cles de JOINT refusees en fem
    void placeTool();
    void computeStableDt();

    // --- per-step pieces -------------------------------------------------
    void applyContact();
    void internalForcesAndDamage();
    void integrate();
    void checkFinite();                  // C4 (w20) : NaN/Inf reel
    void updateDamage(Elem& e, const Eigen::Vector3d& sEff, double szz);
    void erode(Elem& e);
    void refreshActiveNodes();

    // --- data ------------------------------------------------------------
    Config cfg_;
    std::string out_;
    Material mat_;
    PhaseSet phases_;                // >= 1 fiche ; sans `phases` = {mat_}
    bool phasesDeclared_ = false;
    bool voronoi_ = false;           // mesh = voronoi
    int  nGrains_ = 1;
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

    // material-derived constants, UNE ENTREE PAR PHASE (2026-09-07). Sans la
    // cle `phases` il n'y en a qu'une, egale champ pour champ aux anciennes
    // constantes globales : le chemin grille est bit-identique.
    std::vector<Eigen::Matrix3d> DmP_;
    std::vector<double> dpAlphaP_, dpKP_, kappa0TP_, kappa0SP_;

    // ---- amortissement local de Cundall (2026-09-07, opt-in) -------------
    // Porte depuis fdem / fem3d, ou il existait deja. Force non visqueuse
    // -alpha |f| sign(v) par composante : sans dimension, elle ne freine que
    // ce qui bouge et s'annule a l'equilibre. Le FEM 2D n'en avait AUCUN, ce
    // qui rendait toute comparaison FEM / FEMDEM « a conditions egales »
    // fausse par construction (mesure du 2026-09-07 : sur le run de reference
    // l'amortissement dissipe 3,39 J/m contre 1,38 J/m de travail cohesif —
    // il travaille PLUS que la fissuration).
    // Defaut 0 sur les scenarios historiques = chemin bit-identique ; 0,7 en
    // scenario = tension, qui est ne le meme jour, pour s'aligner sur fdem et
    // fem3d en quasi-statique.
    double damping_ = 0.0;

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
    long nanEvery_ = 256, nanStep_ = 0;  // nanCheckEvery (C4, w20)

    // bar-wave verification
    std::vector<int> barBcNodes_, gaugeNodes_;
    double barV0_ = 1.0, gaugeX_ = 0, tArrive_ = -1.0;

    // ---- scenario = tension (2026-09-07) --------------------------------
    // Miroir 2D du scenario TENSION de fem3d : la rangee du bas est encastree
    // (ou tenue en y seulement si gripLateralFree), celle du haut avance a
    // pullV. pullV < 0 = COMPRESSION uniaxiale. La contrainte macroscopique
    // est la reaction des mors divisee par gripSection = W x thk.
    std::vector<uint8_t> grip_;      // 0 libre, 1 mors bas, 2 mors haut
    std::vector<int> midEl_;         // tiers central, jauge de contrainte
    double pullV_ = 0.05, pullRamp_ = 0.0, gripSection_ = 0.0;
    bool   gripFree_ = false;
    Eigen::Vector2d gripF_ = Eigen::Vector2d::Zero();
    double sigmaPeak_ = 0.0, sigMid_ = 0.0;
};

} // namespace rockim
