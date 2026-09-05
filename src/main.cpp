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
#include <sstream>

#include "rockim/Config.hpp"
#include "rockim/Dem3dSolver.hpp"
#include "rockim/DemSolver.hpp"
#include "rockim/Fdem3dSolver.hpp"
#include "rockim/FdemSolver.hpp"
#include "rockim/Fem3dSolver.hpp"
#include "rockim/FemSolver.hpp"
#include "rockim/Guards.hpp"
#include "rockim/KeyGuard.hpp"
#include "rockim/KeysByMode.hpp"
#include "rockim/PotentialContact.hpp"
#include "rockim/Solver.hpp"
#include "rockim/ToolSignorini.hpp"
#include "rockim/ToolPdc3d.hpp"

using namespace rockim;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: rockim <config.cfg> [output_dir]\n"
                     "       rockim selftest-saksala2011 [out.csv]\n"
                     "       rockim selftest-cdp [out.csv]\n"
                     "       rockim selftest-fixed [out.csv]\n"
                     "       rockim matpoint <cfg> [out.csv]\n";
        return 1;
    }

    std::string outDir;                    // connu apres la lecture du deck
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
        if (std::string(argv[1]) == "selftest-triax") {
            std::string csv = argc > 2 ? argv[2] : "rockim_triax.csv";
            int rc = triaxSelftest(csv);
            std::cout << "[rockim] triax selftest traces written to " << csv
                      << "\n";
            return rc;
        }
        // CDP (2026-09-04) : banc point-materiel et pilote generique
        if (std::string(argv[1]) == "selftest-cdp") {
            std::string csv = argc > 2 ? argv[2] : "rockim_cdp.csv";
            int rc = cdpSelftest(csv);
            std::cout << "[rockim] cdp selftest traces written to " << csv
                      << "\n";
            return rc;
        }
        // tensionDamage = fixed (2026-09-05) : bancs falsifiants (a)-(e)
        if (std::string(argv[1]) == "selftest-fixed") {
            std::string csv = argc > 2 ? argv[2] : "rockim_fixed.csv";
            int rc = fixedCrackSelftest(csv);
            std::cout << "[rockim] fixed-crack selftest traces written to " << csv
                      << "\n";
            return rc;
        }
        if (std::string(argv[1]) == "matpoint") {
            if (argc < 3)
                throw std::runtime_error("usage: rockim matpoint <cfg> [out.csv]");
            Config mp = Config::load(argv[2]);
            std::string csv = argc > 3 ? argv[3] : mp.gets("mpOut", "matpoint.csv");
            int rc = matpointDrive(mp, csv);
            std::cout << "[rockim] matpoint trace written to " << csv << "\n";
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
        // T0 (2026-09-02) : contact outil de Signorini, en forme fermee.
        // Sans maillage ni simulation — voir ToolSignorini.hpp.
        if (std::string(argv[1]) == "selftest-toolcontact") {
            std::string csv = argc > 2 ? argv[2] : "rockim_toolcontact.csv";
            int rc = toolSignoriniSelftest(csv);
            std::cout << "[rockim] toolcontact selftest traces written to "
                      << csv << "\n";
            return rc;
        }
        // Cutter PDC 3D (2026-09-03) : geometrie en forme fermee, sans
        // maillage — voir ToolPdc3d.hpp et src/ToolPdc3d.cpp (G1..G5).
        if (std::string(argv[1]) == "selftest-pdc3d") {
            std::string csv = argc > 2 ? argv[2] : "rockim_pdc3d.csv";
            int rc = pdc3dSelftest(csv);
            std::cout << "[rockim] pdc3d selftest traces written to " << csv << "\n";
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
        outDir = out;

        std::string mode = cfg.gets("mode", "fem");
        std::string mesh = cfg.gets("mesh", "grid");
        if (mesh != "grid" && mesh != "voronoi" && mesh != "file")
            throw std::runtime_error("unknown mesh '" + mesh
                                     + "' (grid | voronoi | file)");
        if (mesh == "voronoi" && mode != "fdem" && mode != "fdem3d")
            throw std::runtime_error("mesh = voronoi (grains + phases) is only "
                                     "implemented for mode = fdem | fdem3d");
        if (mesh == "file" && mode != "fdem" && mode != "fdem3d"
            && mode != "fem3d")
            throw std::runtime_error("mesh = file (unstructured import) is "
                                     "only implemented for mode = fdem | "
                                     "fdem3d | fem3d");
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
        // Le couplage thermo-mecanique n existe qu en 2D `fdem`. Le lecteur
        // de configuration ignorant les cles inconnues EN SILENCE, un deck 3D
        // portant `thermal = on` tournerait sans thermique et sans un mot —
        // exactement le resultat faux-mais-plausible que le depot proscrit.
        if (cfg.has("thermal") && mode != "fdem")
            throw std::runtime_error("'thermal' (choc thermique de paroi) "
                                     "n'est implemente que pour mode = fdem "
                                     "(2D). Le deck demande mode = " + mode);
        // Idem pour le litage de Lisjak (bedding*) : 2D `fdem` seulement.
        if (cfg.has("beddingDip") && mode != "fdem")
            throw std::runtime_error("'bedding*' (schistosite pervasive de "
                                     "Lisjak) n'est implemente que pour "
                                     "mode = fdem (2D). Le deck demande "
                                     "mode = " + mode);
        // C6 (w20) : idem pour le couplage hydro-mecanique (AbuAisha) —
        // `hydro` et toute cle `hydro*` : 2D `fdem` seulement. Un deck 3D
        // avec hydro = on tournait sans fluide et sans un mot (05/09).
        // (w21 : keysWithPrefix MARQUE les cles rendues comme consommees —
        // on ne l'appelle donc que hors fdem, pour qu'en fdem une cle hydro*
        // mal orthographiee reste visible par l'audit C1.)
        if (mode != "fdem") {
            auto hk = cfg.keysWithPrefix("hydro");
            if (!hk.empty())
                throw std::runtime_error("'" + hk.front() + "' (couplage "
                    "hydro-mecanique, cavite + fissures) n'est implemente "
                    "que pour mode = fdem (2D). Le deck demande mode = "
                    + mode + " ; retirez les cles hydro*");
        }
        // C6 (w20) : registre des cles PAR MODE (tools/keys_by_mode.json ->
        // include/rockim/KeysByMode.hpp, genere par tools/gen_keys_by_mode.py) :
        // une cle lue par UN SEUL solveur, presente dans un deck d'un autre
        // mode, est refusee (« cle X sans effet en mode Y »). Les cles lues
        // par plusieurs solveurs ou par le code partage sont communes et ne
        // sont jamais refusees. w21 (C1) : ce controle est FONDU dans l'audit
        // des cles consommees fait apres init() (KeyGuard.hpp) — toutes les
        // cles fautives (obsolete, autre mode, faute de frappe, inconnue)
        // sont listees d'un coup et unknownKeys = warn s'y applique aussi ;
        // keysbymode::check() reste disponible mais n'est plus appele ici.
        std::unique_ptr<Solver> solver;
        if      (mode == "fem") solver = std::make_unique<FemSolver>(cfg, out);
        else if (mode == "fem3d") solver = std::make_unique<Fem3dSolver>(cfg, out);
        else if (mode == "dem") solver = std::make_unique<DemSolver>(cfg, out);
        else if (mode == "dem3d") solver = std::make_unique<Dem3dSolver>(cfg, out);
        else if (mode == "fdem") solver = std::make_unique<FdemSolver>(cfg, out);
        else if (mode == "fdem3d") solver = std::make_unique<Fdem3dSolver>(cfg, out);
        else throw std::runtime_error("unknown mode '" + mode + "' (fem | fem3d | dem | dem3d | fdem | fdem3d)");

        solver->init();

        int  nFrames = cfg.geti("frames", 50);
        bool histFlush = cfg.getb("historyFlush", true);
        // C1 (w21, decision de Fernando du 2026-09-05 20:00) : toute cle du
        // deck qu'aucun getter n'a consommee pendant l'initialisation et que
        // ni le solveur du mode courant ni le code partage ne lisent (registre
        // des lecteurs kReaders, w22) est une ERREUR nommee — obsolete (->
        // nouveau nom), d'un autre mode (-> les modes ou elle agit), faute de
        // frappe (suggestion Levenshtein <= 2), ou inconnue de rockim — toutes
        // listees d'un coup, code 1 ; unknownKeys = warn : avertissement et le
        // run continue. Regle detaillee dans KeyGuard.hpp et DOC §8.11.
        // C7 : <outputDir>/config_effective.cfg = les cles consommees et leur
        // valeur effective (deck ou defaut).
        {
            const std::string policy = cfg.gets("unknownKeys", "error");
            auto found = keyguard::enforce(cfg, mode, policy);
            auto eff = cfg.effective();
            std::size_t nDeck = 0;
            for (const auto& e : eff) if (e.fromDeck) ++nDeck;
            std::ostringstream stamp;
            stamp << "exe " << argv[0] << " ; unknownKeys = " << policy;
            keyguard::writeEffective(cfg, out + "/config_effective.cfg", mode, found,
                                     stamp.str());
            std::cout << "[rockim] cles : " << eff.size() << " consommees ("
                      << nDeck << " du deck, " << (eff.size() - nDeck)
                      << " au defaut), " << cfg.size() - nDeck
                      << " du deck non lues" << (found.empty() ? "" : " dont "
                      + std::to_string(found.size()) + " fautives")
                      << " -> " << out << "/config_effective.cfg\n";
            // w22 : fin du suivi — les getters appeles dans step() (ex.
            // confineGaugeTime a chaque pas) lisent la table sans verrou ni
            // chaine de defaut, comme avant w21. Lecture pure.
            cfg.seal();
        }

        long nSteps = (long)std::ceil(solver->duration() / solver->dt());
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
    } catch (const NanError& e) {
        // C4 (w20) : NaN/Inf detecte par le garde reel — pas, temps, noeud,
        // element voisin dans le message ; trace ecrite dans ERROR.txt, code 3
        std::cerr << "[rockim] error: " << e.what() << "\n";
        if (!outDir.empty()) {
            std::ofstream ef(outDir + "/ERROR.txt");
            ef << "[rockim] error: " << e.what() << "\n";
            std::cerr << "[rockim] trace : " << outDir << "/ERROR.txt\n";
        }
        return 3;
    } catch (const std::exception& e) {
        std::cerr << "[rockim] error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
