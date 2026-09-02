#pragma once
// ---------------------------------------------------------------------------
// Material: single definition of the rock, shared by FEM and DEM solvers.
//
// FEM uses (E, nu, rho) + strength (ft, cohesion, phi) + fracture energies.
// DEM derives micro (bond/contact) parameters from the same macro inputs;
// bond strengths can be overridden for calibration.
// Sign convention everywhere: tension positive.
// ---------------------------------------------------------------------------
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "rockim/Config.hpp"

namespace rockim {

struct Material {
    // Elasticity / inertia
    double rho = 2650.0;   // [kg/m^3]
    double E   = 50e9;     // Young's modulus [Pa]
    double nu  = 0.25;     // Poisson ratio

    // Strength (macro)
    double ft        = 10e6;  // uniaxial tensile strength [Pa]
    double cohesion  = 25e6;  // Mohr-Coulomb cohesion [Pa]
    double phiDeg    = 40.0;  // internal friction angle [deg]

    // Fracture energies (FEM crack-band regularization)
    double Gf        = 70.0;  // mode-I fracture energy [J/m^2]
    double gfShearFactor = 10.0; // G_f,II = factor * G_f,I (mode II >> mode I in rock)

    double G() const { return E / (2.0 * (1.0 + nu)); }
    double K() const { return E / (3.0 * (1.0 - 2.0 * nu)); }

    // Plane-strain P-wave (dilatational) speed: the fastest wave in the mesh,
    // which controls the FEM critical time step.
    //   c_p = sqrt( E(1-nu) / ((1+nu)(1-2nu) rho) )
    double cP() const {
        return std::sqrt(E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu) * rho));
    }

    // 1D bar speed sqrt(E/rho) — used by the bar-wave verification (with nu=0
    // the plane-strain modulus reduces to E, so the model must recover this).
    double cBar() const { return std::sqrt(E / rho); }

    // Shear wave speed sqrt(G/rho) — tangential impedance of the
    // Lysmer-Kuhlemeyer absorbing boundaries.
    double cS() const { return std::sqrt(E / (2.0 * (1.0 + nu)) / rho); }

    // Drucker-Prager parameters matched to Mohr-Coulomb (c, phi) under PLANE
    // STRAIN (Chen & Han):  F = sqrt(J2) + alpha*I1 - k <= 0
    //   alpha = tan(phi) / sqrt(9 + 12 tan^2(phi))
    //   k     = 3 c      / sqrt(9 + 12 tan^2(phi))
    void dpParams(double& alpha, double& k) const {
        double t = std::tan(phiDeg * M_PI / 180.0);
        double den = std::sqrt(9.0 + 12.0 * t * t);
        alpha = t / den;
        k     = 3.0 * cohesion / den;
    }

    // Plane-strain elastic matrix for (eps_xx, eps_yy, gamma_xy) -> (s_xx, s_yy, t_xy)
    Eigen::Matrix3d Dmat() const {
        double f = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        Eigen::Matrix3d D;
        D << f * (1.0 - nu), f * nu,          0.0,
             f * nu,         f * (1.0 - nu),  0.0,
             0.0,            0.0,             f * (1.0 - 2.0 * nu) / 2.0;
        return D;
    }

    static Material from(const Config& c) {
        Material m;
        m.rho = c.getd("rho", m.rho);
        m.E   = c.getd("E",   m.E);
        m.nu  = c.getd("nu",  m.nu);
        m.ft        = c.getd("ft", m.ft);
        m.cohesion  = c.getd("cohesion", m.cohesion);
        m.phiDeg    = c.getd("frictionDeg", m.phiDeg);
        m.Gf        = c.getd("Gf", m.Gf);
        m.gfShearFactor = c.getd("gfShearFactor", m.gfShearFactor);
        return m;
    }
};

// ---------------------------------------------------------------------------
// PhaseSet: one Material per mineral phase, plus the grain-boundary joint
// rule. Config syntax (all per-phase keys default to the global material, so
// only the differences need spelling out; fraction is required):
//
//   phases = quartz feldspar biotite
//   phase.quartz.fraction  = 0.33
//   phase.quartz.E         = 94e9
//   phase.quartz.ft        = 13e6
//   ...
//
// Grain-boundary joints take the MEAN of the two neighbouring phases times an
// attenuation factor per property (the alpha coefficients of the GBM
// literature); heterophase boundaries (different minerals) get one extra
// multiplier on the strength-like properties:
//
//   gbAlphaTen / gbAlphaCoh / gbAlphaGf / gbAlphaE / gbAlphaFric   (default 1)
//   gbHeteroFactor                                                  (default 1)
//
// Without a 'phases' key the set degenerates to the single global material
// and every joint keeps the bulk properties — the pre-GBM behaviour.
// ---------------------------------------------------------------------------
struct PhaseSet {
    std::vector<Material>    mat;       // >= 1 entries
    std::vector<double>      fraction;  // normalized to sum 1
    std::vector<std::string> name;

    double aTen = 1.0, aCoh = 1.0, aGf = 1.0, aE = 1.0, aFric = 1.0;
    double heteroFactor = 1.0;

    // ---- proprietes de joint PAR PAIRE DE PHASES (2026-09-01) -------------
    // Les alpha ci-dessus donnent UN facteur pour tout le reseau ; la
    // litterature GBM tabule au contraire chaque paire (Aboayanah et al.
    // RMRE 2024, Table 2 : six paires Bt-Bt, Fsp-Fsp, Qz-Qz, Bt-Fsp, Bt-Qz,
    // Qz-Fsp, dont les GfI vont de 1,07 a 907 J/m2 — deux ordres de
    // grandeur que deux boutons ne peuvent pas representer).
    // Syntaxe (opt-in, l'ordre des deux phases est indifferent) :
    //   gb.feldspar.quartz.ft            = 2.0e6
    //   gb.feldspar.quartz.cohesion      = 32e6
    //   gb.feldspar.quartz.Gf            = 300
    //   gb.feldspar.quartz.gfShearFactor = 4.83
    //   gb.feldspar.quartz.frictionDeg   = 39.0
    //   gb.feldspar.quartz.E             = 71.5e9
    // Une valeur posee REMPLACE le resultat de la regle alpha/hetero pour
    // cette paire (elle ne s'y multiplie pas) ; une valeur absente laisse
    // la regle alpha agir. Sans aucune cle gb.<a>.<b>.*, chemin et
    // resultats strictement inchanges.
    struct GbPair {
        double ft = -1.0, coh = -1.0, Gf = -1.0, gfs = -1.0;
        double phiDeg = -1.0, E = -1.0;
        bool any = false;
    };
    std::vector<GbPair> pairGb;            // n x n, symetrique ; vide = aucun

    const GbPair* gbOf(int a, int b) const {
        if (pairGb.empty() || a < 0 || b < 0 || a >= n() || b >= n())
            return nullptr;
        const GbPair& g = pairGb[(std::size_t)a * n() + b];
        return g.any ? &g : nullptr;
    }

    std::string gbOverrideSummary() const {
        if (pairGb.empty()) return {};
        std::string s;
        for (int i = 0; i < n(); ++i)
            for (int j = i; j < n(); ++j) {
                const GbPair* g = gbOf(i, j);
                if (!g) continue;
                s += "  " + name[i] + "-" + name[j] + " :";
                if (g->ft > 0)  s += " ft " + std::to_string(g->ft / 1e6) + " MPa";
                if (g->coh > 0) s += " c " + std::to_string(g->coh / 1e6) + " MPa";
                if (g->Gf > 0)  s += " GfI " + std::to_string(g->Gf);
                if (g->gfs > 0) s += " x" + std::to_string(g->gfs);
                if (g->phiDeg > 0) s += " phi " + std::to_string(g->phiDeg);
                if (g->E > 0)   s += " E " + std::to_string(g->E / 1e9) + " GPa";
                s += "\n";
            }
        return s;
    }

    int    n() const { return (int)mat.size(); }
    double maxE() const {
        double e = 0.0;
        for (const auto& m : mat) e = std::max(e, m.E);
        return e;
    }
    double maxCp() const {
        double c = 0.0;
        for (const auto& m : mat) c = std::max(c, m.cP());
        return c;
    }

    // A zero or negative strength/stiffness silently INVERTS the joint
    // semantics downstream (ft = 0 makes dnF infinite, the cohesive
    // envelope NaN, and the joint unbreakable instead of strengthless), so
    // material validation is strict: model pre-cracked boundaries with a
    // SMALL positive alpha (e.g. 1e-3), not zero.
    static void validate(const Material& m, const std::string& who) {
        auto bad = [&](const std::string& what) {
            throw std::runtime_error("Material (" + who + "): " + what);
        };
        if (!(m.E > 0.0))   bad("E must be > 0");
        if (!(m.rho > 0.0)) bad("rho must be > 0");
        if (!(m.nu >= 0.0 && m.nu < 0.5)) bad("nu must be in [0, 0.5)");
        if (!(m.ft > 0.0))       bad("ft must be > 0");
        if (!(m.cohesion > 0.0)) bad("cohesion must be > 0");
        if (!(m.Gf > 0.0))       bad("Gf must be > 0");
        if (!(m.gfShearFactor > 0.0)) bad("gfShearFactor must be > 0");
        if (!(m.phiDeg >= 0.0 && m.phiDeg < 89.0))
            bad("frictionDeg must be in [0, 89)");
    }

    static PhaseSet from(const Config& c) {
        PhaseSet ps;
        Material base = Material::from(c);
        ps.aTen  = c.getd("gbAlphaTen", 1.0);
        ps.aCoh  = c.getd("gbAlphaCoh", 1.0);
        ps.aGf   = c.getd("gbAlphaGf", 1.0);
        ps.aE    = c.getd("gbAlphaE", 1.0);
        ps.aFric = c.getd("gbAlphaFric", 1.0);
        ps.heteroFactor = c.getd("gbHeteroFactor", 1.0);
        for (auto [v, nm] : {std::pair<double, const char*>{ps.aTen, "gbAlphaTen"},
                             {ps.aCoh, "gbAlphaCoh"}, {ps.aGf, "gbAlphaGf"},
                             {ps.aE, "gbAlphaE"}, {ps.aFric, "gbAlphaFric"},
                             {ps.heteroFactor, "gbHeteroFactor"}})
            if (!(v > 0.0))
                throw std::runtime_error(std::string("PhaseSet: ") + nm
                    + " must be > 0 (use a small value like 1e-3 for "
                      "pre-cracked boundaries, not 0)");

        std::string list = c.gets("phases", "");
        if (list.empty()) {
            validate(base, "global");
            ps.mat = {base};
            ps.fraction = {1.0};
            ps.name = {"rock"};
            return ps;
        }
        std::istringstream ss(list);
        std::string nm;
        double fsum = 0.0;
        while (ss >> nm) {
            std::string k = "phase." + nm + ".";
            Material m = base;
            m.rho = c.getd(k + "rho", m.rho);
            m.E   = c.getd(k + "E",   m.E);
            m.nu  = c.getd(k + "nu",  m.nu);
            m.ft       = c.getd(k + "ft", m.ft);
            m.cohesion = c.getd(k + "cohesion", m.cohesion);
            m.phiDeg   = c.getd(k + "frictionDeg", m.phiDeg);
            m.Gf       = c.getd(k + "Gf", m.Gf);
            m.gfShearFactor = c.getd(k + "gfShearFactor", m.gfShearFactor);
            validate(m, "phase " + nm);
            double f = c.reqd(k + "fraction");
            if (f <= 0.0)
                throw std::runtime_error("PhaseSet: fraction of phase '" + nm
                                         + "' must be positive");
            ps.mat.push_back(m);
            ps.fraction.push_back(f);
            ps.name.push_back(nm);
            fsum += f;
        }
        if (ps.mat.empty())
            throw std::runtime_error("PhaseSet: 'phases' key present but empty");
        for (double& f : ps.fraction) f /= fsum;

        // ---- surcharges par paire de phases (voir GbPair ci-dessus) ------
        // Lecture directe par cle : pas besoin d'enumerer le fichier, et
        // les deux ordres d'ecriture sont acceptes. Une seule cle presente
        // suffit a armer la paire ; les autres proprietes restent sous la
        // regle alpha.
        const int np = ps.n();
        std::vector<GbPair> pg((std::size_t)np * np);
        bool anyPair = false;
        auto rd = [&](const std::string& a, const std::string& b,
                      const char* prop) {
            double v = c.getd("gb." + a + "." + b + "." + prop, -1.0);
            if (v < 0.0) v = c.getd("gb." + b + "." + a + "." + prop, -1.0);
            return v;
        };
        for (int i = 0; i < np; ++i)
            for (int j = i; j < np; ++j) {
                GbPair g;
                g.ft     = rd(ps.name[i], ps.name[j], "ft");
                g.coh    = rd(ps.name[i], ps.name[j], "cohesion");
                g.Gf     = rd(ps.name[i], ps.name[j], "Gf");
                g.gfs    = rd(ps.name[i], ps.name[j], "gfShearFactor");
                g.phiDeg = rd(ps.name[i], ps.name[j], "frictionDeg");
                g.E      = rd(ps.name[i], ps.name[j], "E");
                g.any = g.ft > 0 || g.coh > 0 || g.Gf > 0 || g.gfs > 0
                        || g.phiDeg > 0 || g.E > 0;
                if (!g.any) continue;
                if (g.phiDeg > 0 && !(g.phiDeg < 89.0))
                    throw std::runtime_error(
                        "PhaseSet: gb." + ps.name[i] + "." + ps.name[j]
                        + ".frictionDeg doit etre dans [0, 89)");
                if (g.gfs > 0 && !(g.Gf > 0))
                    throw std::runtime_error(
                        "PhaseSet: gb." + ps.name[i] + "." + ps.name[j]
                        + ".gfShearFactor sans .Gf — poser les deux (GfII = "
                          "Gf x gfShearFactor)");
                anyPair = true;
                pg[(std::size_t)i * np + j] = g;
                pg[(std::size_t)j * np + i] = g;
            }
        if (anyPair) ps.pairGb.swap(pg);
        return ps;
    }
};

} // namespace rockim
