// ---------------------------------------------------------------------------
// rockim — 2D rock impact & shear-failure simulator (FEM + DEM).
//
//   usage: rockim <config.cfg> [output_dir]
//
// The driver only talks to the abstract Solver interface: it loads the shared
// problem definition, builds the requested solver (mode = fem | dem), runs the
// explicit time loop, and writes history.csv + VTK frames. A future FDEM
// solver drops in here with one extra line in the factory.
// ---------------------------------------------------------------------------
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>

#include "rockim/Config.hpp"
#include "rockim/Dem3dSolver.hpp"
#include "rockim/DemSolver.hpp"
#include "rockim/Fdem3dSolver.hpp"
#include "rockim/FdemSolver.hpp"
#include "rockim/Fem3dSolver.hpp"
#include "rockim/FemSolver.hpp"
#include "rockim/PotentialContact.hpp"
#include "rockim/Solver.hpp"

using namespace rockim;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: rockim <config.cfg> [output_dir]\n"
                     "       rockim selftest-saksala2011 [out.csv]\n";
        return 1;
    }

    try {
        if (std::string(argv[1]) == "selftest-saksala2011") {
            std::string csv = argc > 2 ? argv[2] : "rockim_saksala.csv";
            saksala2011Selftest(csv);
            std::cout << "[rockim] saksala2011 selftest traces written to "
                      << csv << "\n";
            return 0;
        }
        if (std::string(argv[1]) == "selftest-mc") {
            std::string csv = argc > 2 ? argv[2] : "rockim_mc.csv";
            int rc = mcSelftest(csv);
            std::cout << "[rockim] mc selftest traces written to " << csv
                      << "\n";
            return rc;
        }
        if (std::string(argv[1]) == "selftest-dpdfh") {
            std::string csv = argc > 2 ? argv[2] : "rockim_dpdfh.csv";
            dpdfhSelftest(csv);
            std::cout << "[rockim] dpdfh selftest traces written to "
                      << csv << "\n";
            return 0;
        }
        if (std::string(argv[1]) == "selftest-potential2d") {
            std::string csv = argc > 2 ? argv[2] : "rockim_potential2d.csv";
            int rc = potentialSelftest(csv);
            std::cout << "[rockim] potential2d selftest traces written to "
                      << csv << "\n";
            return rc;
        }
        if (std::string(argv[1]) == "selftest-potential3d") {
            std::string csv = argc > 2 ? argv[2] : "rockim_potential3d.csv";
            int rc = potentialSelftest3d(csv);
            std::cout << "[rockim] potential3d selftest traces written to "
                      << csv << "\n";
            return rc;
        }
        Config cfg = Config::load(argv[1]);
        std::string out = (argc > 2) ? argv[2] : cfg.gets("outputDir", "out");
        std::filesystem::create_directories(out);

        std::string mode = cfg.gets("mode", "fem");
        std::string mesh = cfg.gets("mesh", "grid");
        if (mesh != "grid" && mesh != "voronoi" && mesh != "file")
            throw std::runtime_error("unknown mesh '" + mesh
                                     + "' (grid | voronoi | file)");
        if (mesh == "voronoi" && mode != "fdem" && mode != "fdem3d")
            throw std::runtime_error("mesh = voronoi (grains + phases) is only "
                                     "implemented for mode = fdem | fdem3d");
        if (mesh == "file" && mode != "fdem" && mode != "fdem3d")
            throw std::runtime_error("mesh = file (unstructured import) is "
                                     "only implemented for mode = fdem | "
                                     "fdem3d");
        if (cfg.has("phases") && mode != "fdem" && mode != "fdem3d")
            throw std::runtime_error("'phases' (mineral phases) is only "
                                     "implemented for mode = fdem | fdem3d");
        if (cfg.has("law") && mode != "fem3d" && mode != "fdem"
            && mode != "fdem3d")
            throw std::runtime_error("'law' (bulk constitutive law) is "
                                     "implemented for mode = fem3d (3D "
                                     "continuum) and mode = fdem (2D bulk in "
                                     "PLANE STRAIN, coupled with the cohesive "
                                     "joints — the vumat_fdem_coupled "
                                     "configuration) and mode = fdem3d (3D "
                                     "bulk, no plane-strain trick needed)");
        std::unique_ptr<Solver> solver;
        if      (mode == "fem") solver = std::make_unique<FemSolver>(cfg, out);
        else if (mode == "fem3d") solver = std::make_unique<Fem3dSolver>(cfg, out);
        else if (mode == "dem") solver = std::make_unique<DemSolver>(cfg, out);
        else if (mode == "dem3d") solver = std::make_unique<Dem3dSolver>(cfg, out);
        else if (mode == "fdem") solver = std::make_unique<FdemSolver>(cfg, out);
        else if (mode == "fdem3d") solver = std::make_unique<Fdem3dSolver>(cfg, out);
        else throw std::runtime_error("unknown mode '" + mode + "' (fem | fem3d | dem | dem3d | fdem | fdem3d)");

        solver->init();

        long nSteps = (long)std::ceil(solver->duration() / solver->dt());
        int  nFrames = cfg.geti("frames", 50);
        long outEvery  = std::max(1L, nSteps / std::max(1, nFrames));
        long histEvery = std::max(1L, nSteps / 2000);

        std::ofstream hist(out + "/history.csv");
        solver->historyHeader(hist);
        // history.csv est vide jusqu'a la fin si on laisse l'OS bufferiser :
        // impossible de suivre un run en cours, et un run TUE laisse une
        // derniere ligne tronquee au milieu du tampon (constate le 2026-08-14
        // sur out_banc_mid : 26 colonnes au lieu de 28, terminees par ",-").
        // On vide donc apres CHAQUE ligne : histEvery borne le nombre de
        // lignes a ~2000 sur tout le run, le cout est negligeable, et le
        // fichier se termine toujours sur une ligne complete. Purement I/O :
        // aucun effet sur le calcul (bit-neutre par construction).
        bool histFlush = cfg.getb("historyFlush", true);
        auto histRow = [&] {
            solver->historyRow(hist);
            if (histFlush) hist.flush();
        };
        hist.flush();

        auto t0 = std::chrono::steady_clock::now();
        int frame = 0;
        long nextPct = 10;
        for (long i = 0; i < nSteps; ++i) {
            if (i % outEvery == 0) solver->writeFrame(frame++);
            if (i % histEvery == 0) histRow();
            solver->step();
            if (solver->finished()) {
                std::cout << "\n[rockim] solver requested an early stop at t = "
                          << solver->time() << " s (" << i + 1 << " / " << nSteps
                          << " steps)\n";
                solver->writeFrame(frame++);
                histRow();
                break;
            }
            if (100 * (i + 1) / nSteps >= nextPct) {
                std::cout << "  " << nextPct << "%" << std::flush
                          << (nextPct == 100 ? "\n" : " ");
                nextPct += 10;
            }
        }
        solver->writeFrame(frame);
        histRow();
        solver->finalize();

        auto t1 = std::chrono::steady_clock::now();
        std::cout << "[rockim] wall time: "
                  << std::chrono::duration<double>(t1 - t0).count() << " s, output in '"
                  << out << "'\n";
    } catch (const std::exception& e) {
        std::cerr << "[rockim] error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
