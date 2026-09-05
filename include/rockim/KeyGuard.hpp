#pragma once
// ---------------------------------------------------------------------------
// KeyGuard.hpp — C1 « la cle inconnue est une ERREUR » et C7 config_effective
// (plan de robustesse du 2026-09-05, chantier C ; decision de Fernando 20:00 ;
// binaire w21, corrige w22 apres relecture : D1 lecteurs, D2 rejeu).
//
// Appele par main.cpp APRES solver->init() (toutes les lectures conditionnelles
// de l'initialisation ont eu lieu ; Config a marque chaque cle CONSOMMEE par
// gets/getd/geti/getb/reqs/reqd/has/keysWithPrefix — Config.hpp).
//
// REGLE, pour chaque cle du deck JAMAIS consommee :
//   1. nom dans la table des OBSOLETES (tools/obsolete_keys.json -> kObsolete)
//        -> « cle obsolete X (ligne N du deck) : remplacee par Z » ;
//   2. cle du registre (kReaders, w22 : les LECTEURS de chaque cle = solveurs
//      dont le fichier la lit, ou `shared` = code partage / ALWAYS_COMMON) :
//        - le mode courant OU `shared` est lecteur -> LEGITIME, rien. C'est la
//          lecture CONDITIONNELLE : capP0 n'est lu que si dprCap = true,
//          hydroStart que si hydro = on, jointResidualMu que sous certains
//          adoucissements ; gravity n'est lu qu'apres init (placeTool). Une cle
//          que le solveur du mode lit dans un chemin du code existe et agit :
//          la refuser serait un faux refus. (Le « reglage inerte » d'une telle
//          cle reste possible sans message : prix assume de zero faux refus.)
//        - sinon -> « cle X sans effet en mode Y : cle du mode Z seulement »
//          (un lecteur) ou « … : cles des modes Z1, Z2 » (plusieurs). w21
//          traitait toute cle a >= 2 lecteurs comme « commune = legitime » :
//          jointSoftening, insertion, bulkDamage, fragBrushV0 (fdem + fdem3d)
//          passaient en fem3d sans un mot — trou D1 de la relecture, ferme ici.
//   3. cle d'une famille DYNAMIQUE (kDynamicPrefix : phase.<nom>.E,
//      groupBond.<A>.<B>, gauge.<groupe>...) non lue -> « cle X (famille F) jamais
//      lue : nom de phase / groupe inconnu ? » ;
//   4. sinon, suggestion de Levenshtein <= 2 sur l'union (registre + cles
//      consommees) -> « cle inconnue X (ligne N du deck) : vouliez-vous dire Y ? »
//      (si Y n'agit pas dans ce mode, on dit dans lesquels) ;
//   5. sinon « cle X (ligne N du deck) inconnue de rockim ».
// TOUTES les cles fautives sont listees d'un coup, puis :
//   unknownKeys = error (DEFAUT)  -> runtime_error, code de retour 1 ;
//   unknownKeys = warn            -> memes messages en [rockim] WARNING, le run
//                                    continue (vieux decks).
// Lecture PURE : aucun flottant du calcul n'est touche (bit-identite).
//
// C7 : writeEffective() ecrit <outputDir>/config_effective.cfg :
//   - lignes ACTIVES = toutes les cles du deck sauf les fautives (consommees :
//     « # deck ligne n » ; non consommees mais legitimes : « # deck ligne n ;
//     non lue a l'initialisation ») ;
//   - defauts du code consommes = lignes COMMENTEES « # cle = valeur   (defaut) »
//     (w22, D2 : ecrits actifs, ils devenaient des cles POSEES et declenchaient
//     les gardes « satellite orpheline » — tensionShearRetention en fem3d,
//     dampingLocalAfter en fdem — au rejeu) ;
//   tri stable par cle. Le fichier est un deck REJOUABLE : ses cles actives sont
//   celles du deck d'origine moins les refusees -> memes resultats (verifie par
//   le banc selftest_cles partie C : history.csv bit-identique sur 6 modes).
// ---------------------------------------------------------------------------
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "rockim/Config.hpp"
#include "rockim/KeysByMode.hpp"

namespace rockim {
namespace keyguard {

struct Finding {
    std::string key;
    int line = 0;
    std::string kind;      // obsolete | mode | dynamic | typo | unknown
    std::string msg;
};

// distance de Levenshtein (edition), bornee : rend > maxD des que depasse
inline int levenshtein(const std::string& a, const std::string& b, int maxD) {
    const int n = (int)a.size(), m = (int)b.size();
    if (std::abs(n - m) > maxD) return maxD + 1;
    std::vector<int> prev(m + 1), cur(m + 1);
    for (int j = 0; j <= m; ++j) prev[j] = j;
    for (int i = 1; i <= n; ++i) {
        cur[0] = i;
        int rowMin = cur[0];
        for (int j = 1; j <= m; ++j) {
            int c = (a[i - 1] == b[j - 1]) ? 0 : 1;
            cur[j] = std::min({prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + c});
            rowMin = std::min(rowMin, cur[j]);
        }
        if (rowMin > maxD) return maxD + 1;
        std::swap(prev, cur);
    }
    return prev[m];
}

inline std::string lower(std::string s) {
    for (auto& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

inline bool isKnown(const std::string& key) {
    for (const char* const* k = keysbymode::kKnown; *k; ++k)
        if (key == *k) return true;
    return false;
}

inline const char* obsoleteOf(const std::string& key) {
    for (const keysbymode::Obsolete* o = keysbymode::kObsolete; o->old; ++o)
        if (key == o->old) return o->now;
    return nullptr;
}

inline const char* dynamicFamily(const std::string& key) {
    for (const char* const* p = keysbymode::kDynamicPrefix; *p; ++p)
        if (key.compare(0, std::string(*p).size(), *p) == 0) return *p;
    return nullptr;
}

// Lecteurs d'une cle du registre (w22), decoupes ("fdem fdem3d" -> {fdem,
// fdem3d}) ; vide si la cle est hors registre.
inline std::vector<std::string> readersOf(const std::string& key) {
    std::vector<std::string> out;
    for (const keysbymode::Readers* r = keysbymode::kReaders; r->key; ++r)
        if (key == r->key) {
            std::istringstream is(r->readers);
            std::string w;
            while (is >> w) out.push_back(w);
            break;
        }
    return out;
}

// La cle agit-elle dans `mode` ? (mode lecteur, ou code partage)
inline bool actsIn(const std::vector<std::string>& readers, const std::string& mode) {
    for (const auto& r : readers)
        if (r == mode || r == "shared") return true;
    return false;
}

// « cle du mode Z seulement » / « cles des modes Z1, Z2 »
inline std::string modesText(const std::vector<std::string>& readers) {
    std::ostringstream os;
    if (readers.size() == 1) os << "cle du mode " << readers[0] << " seulement";
    else {
        os << "cles des modes ";
        for (std::size_t i = 0; i < readers.size(); ++i) os << (i ? ", " : "") << readers[i];
    }
    return os.str();
}

// Suggestions (distance <= 2, casse ignoree en premier) parmi `cands` ;
// au plus 3, les plus proches d'abord, ordre lexicographique a egalite.
inline std::vector<std::string> suggest(const std::string& key,
                                        const std::set<std::string>& cands) {
    std::vector<std::pair<int, std::string>> hits;
    const std::string lk = lower(key);
    for (const auto& c : cands) {
        if (c == key) continue;
        int d = (lower(c) == lk) ? 0 : levenshtein(key, c, 2);
        if (d <= 2) hits.emplace_back(d, c);
    }
    std::sort(hits.begin(), hits.end());
    std::vector<std::string> out;
    for (const auto& h : hits) {
        if (out.size() >= 3) break;
        out.push_back(h.second);
    }
    return out;
}

// Audit du deck : toutes les cles fautives, ordre du deck (ligne croissante).
inline std::vector<Finding> audit(const Config& cfg, const std::string& mode) {
    std::set<std::string> cands;
    for (const char* const* k = keysbymode::kKnown; *k; ++k) cands.insert(*k);
    for (const auto& k : cfg.consumedKeys()) cands.insert(k);

    std::vector<Finding> out;
    for (const auto& key : cfg.unusedKeys()) {
        Finding f;
        f.key = key;
        f.line = cfg.lineOf(key);
        std::ostringstream os;
        const std::string where = " (ligne " + std::to_string(f.line) + " du deck)";
        const auto readers = readersOf(key);
        if (const char* now = obsoleteOf(key)) {
            f.kind = "obsolete";
            os << "cle obsolete '" << key << "'" << where << " : remplacee par '"
               << now << "' (tools/obsolete_keys.json)";
        } else if (!readers.empty()) {
            if (actsIn(readers, mode)) continue;   // lue sous condition / apres init : legitime
            f.kind = "mode";
            os << "cle '" << key << "' sans effet en mode " << mode << " : "
               << modesText(readers) << where
               << " (registre tools/keys_by_mode.json) ; retirez-la du deck ou changez de mode";
        } else if (const char* fam = dynamicFamily(key)) {
            f.kind = "dynamic";
            os << "cle '" << key << "'" << where << " jamais lue : famille dynamique '"
               << fam << "<nom>...' — le nom de phase / groupe qu'elle porte n'est "
                  "declare nulle part (phases, groupes du maillage) ?";
        } else {
            auto sug = suggest(key, cands);
            if (!sug.empty()) {
                f.kind = "typo";
                os << "cle inconnue '" << key << "'" << where << " : vouliez-vous dire ";
                for (std::size_t i = 0; i < sug.size(); ++i) {
                    os << (i ? " ou " : "") << "'" << sug[i] << "'";
                    const auto sr = readersOf(sug[i]);
                    if (!sr.empty() && !actsIn(sr, mode)) os << " (" << modesText(sr) << ")";
                }
                os << " ?";
            } else {
                f.kind = "unknown";
                os << "cle '" << key << "'" << where << " inconnue de rockim (aucun getter "
                      "du code ne la lit, aucune cle connue a moins de 2 caracteres)";
            }
        }
        f.msg = os.str();
        out.push_back(f);
    }
    std::stable_sort(out.begin(), out.end(),
                     [](const Finding& a, const Finding& b) { return a.line < b.line; });
    return out;
}

// Applique la politique `unknownKeys` (error | warn). Rend les fautes (pour
// config_effective.cfg) ; leve runtime_error en mode error s'il y en a.
inline std::vector<Finding> enforce(const Config& cfg, const std::string& mode,
                                    const std::string& policy) {
    if (policy != "error" && policy != "warn")
        throw std::runtime_error("unknownKeys = '" + policy + "' : valeurs admises error "
                                 "(defaut : cle inconnue = erreur, decision du 2026-09-05) | "
                                 "warn (avertissement, le run continue)");
    auto found = audit(cfg, mode);
    if (found.empty()) return found;
    std::ostringstream os;
    os << found.size() << " cle" << (found.size() > 1 ? "s" : "") << " du deck '"
       << cfg.path() << "' refusee" << (found.size() > 1 ? "s" : "")
       << " (unknownKeys = " << policy << ")";
    if (policy == "error") {
        os << " ; poser unknownKeys = warn pour continuer avec un avertissement :";
        for (const auto& f : found) os << "\n  - " << f.msg;
        throw std::runtime_error(os.str());
    }
    std::cerr << "[rockim] WARNING: " << os.str() << " — le run continue :\n";
    for (const auto& f : found) std::cerr << "[rockim] WARNING:   - " << f.msg << "\n";
    return found;
}

// C7 : <outputDir>/config_effective.cfg (format : en-tete de ce fichier)
inline void writeEffective(const Config& cfg, const std::string& path,
                           const std::string& mode, const std::vector<Finding>& found,
                           const std::string& stamp) {
    std::ofstream f(path);
    if (!f) return;                       // I/O seulement : jamais bloquant
    std::set<std::string> faulty;
    for (const auto& x : found) faulty.insert(x.key);
    const auto deck = cfg.deck();         // cles du deck, tri par cle
    const auto eff = cfg.effective();     // consommees (deck + defauts), tri par cle
    struct Row { std::string key, text; };
    std::vector<Row> rows;
    std::size_t nDeck = 0, nUnread = 0, nDef = 0;
    for (const auto& k : deck) {
        if (faulty.count(k.key)) continue;
        Row r;
        r.key = k.key;
        r.text = k.key + " = " + k.value;
        if (k.consumed) { ++nDeck; r.text += "\t# deck ligne " + std::to_string(k.line); }
        else {
            ++nUnread;
            r.text += "\t# deck ligne " + std::to_string(k.line)
                    + " ; non lue a l'initialisation (lecture conditionnelle ou apres init)";
        }
        rows.push_back(r);
    }
    for (const auto& e : eff) {
        if (e.fromDeck) continue;
        ++nDef;
        Row r;
        r.key = e.key;
        r.text = "# " + e.key + " = " + e.value + "\t(defaut)";
        if (e.defaults.size() > 1) {
            r.text += " ; autres defauts lus :";
            for (std::size_t i = 1; i < e.defaults.size(); ++i) r.text += " " + e.defaults[i];
        }
        rows.push_back(r);
    }
    std::stable_sort(rows.begin(), rows.end(),
                     [](const Row& a, const Row& b) { return a.key < b.key; });
    f << "# config_effective.cfg — rockim (C7, w22) : cles du deck retenues + defauts du code consommes\n"
      << "# deck : " << cfg.path() << " ; mode : " << mode << " ; " << stamp << "\n"
      << "# " << (nDeck + nDef) << " cles consommees a l'initialisation : " << nDeck
      << " du deck (# deck ligne n), " << nDef << " au defaut du code (lignes commentees « # cle = valeur (defaut) »)"
      << " ; " << nUnread << " cle" << (nUnread > 1 ? "s" : "") << " du deck non lue"
      << (nUnread > 1 ? "s" : "") << " a l'initialisation mais legitime" << (nUnread > 1 ? "s" : "")
      << " (gardee" << (nUnread > 1 ? "s" : "") << " active" << (nUnread > 1 ? "s" : "") << ") ; tri par cle.\n"
      << "# Deck REJOUABLE : lignes actives = cles du deck d'origine moins les refusees ; les defauts sont\n"
      << "# commentes pour ne pas devenir des cles posees (gardes « satellite orpheline » du code).\n";
    if (!found.empty()) {
        f << "# cles du deck REFUSEES / AVERTIES (" << found.size() << "), retirees des lignes actives :\n";
        for (const auto& x : found) f << "#   - " << x.msg << "\n";
    }
    f << "\n";
    // alignement des commentaires : colonne commune
    std::size_t w = 0;
    for (const auto& r : rows) {
        auto t = r.text.find('\t');
        if (t != std::string::npos) w = std::max(w, t);
    }
    for (const auto& r : rows) {
        auto t = r.text.find('\t');
        if (t == std::string::npos) { f << r.text << "\n"; continue; }
        f << r.text.substr(0, t) << std::string(w + 4 - t, ' ') << r.text.substr(t + 1) << "\n";
    }
}

} // namespace keyguard
} // namespace rockim
