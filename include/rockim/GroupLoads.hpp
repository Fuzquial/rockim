#pragma once
// ---------------------------------------------------------------------------
// GroupLoads.hpp — briques COMMUNES des chargements et conditions aux limites
// par groupes (scenario = loads, cles fix./velocity./traction./pressure./
// force./amplitude./point./box.), partagees par fdem3d (Fdem3dLoads.cpp) et
// fem3d (Fem3dLoads.cpp) depuis le 2026-10-03.
//
// Seulement ce qui ne depend pas du solveur : lecture stricte des nombres
// (virgule decimale refusee), table d'amplitude (lecture et evaluation),
// axes de fix.<g>, vitesse de velocity.<g>. Les messages d'erreur sont ceux
// de la premiere version de Fdem3dLoads.cpp, au caractere pres : un deck
// reste interchangeable entre les deux modes (a la ligne `mode` pres).
// Les groupes eux-memes (copies de noeuds en fdem3d, noeuds partages en
// fem3d) restent dans chaque solveur.
// ---------------------------------------------------------------------------
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace rockim {
namespace grouploads {

// prefixes des cles, dans l'ordre de lecture
inline const char* const kPrefixes[] = {"fix.", "velocity.", "traction.",
                                        "pressure.", "force.", "amplitude.",
                                        "point.", "box."};

inline std::vector<std::string> tokens(const std::string& s) {
    std::istringstream iss(s);
    std::vector<std::string> out;
    std::string t;
    while (iss >> t) out.push_back(t);
    return out;
}

// nombre strict : la virgule decimale et les suffixes parasites sont refuses
// (meme regle que Config::getd), avec le nom de la cle dans le message
inline double num(const std::string& key, const std::string& t) {
    std::size_t pos = 0;
    double v = 0.0;
    bool ok = true;
    try { v = std::stod(t, &pos); } catch (...) { ok = false; }
    if (!ok || pos != t.size() || !std::isfinite(v))
        throw std::runtime_error(key + " : valeur numerique invalide '" + t
                                 + "' (point decimal obligatoire)");
    return v;
}

inline std::vector<double> nums(const std::string& key, const std::string& s,
                                std::size_t n) {
    auto tk = tokens(s);
    if (tk.size() != n)
        throw std::runtime_error(key + " : attendu " + std::to_string(n)
                                 + " valeurs, lu '" + s + "'");
    std::vector<double> v;
    for (const auto& t : tk) v.push_back(num(key, t));
    return v;
}

// amplitude.<g> = t0 a0 t1 a1 ... | ramp T
struct Amp {
    std::vector<double> t, a;
    bool smooth = false;                       // ramp T : montee en cosinus
};

inline Amp parseAmp(const std::string& k, const std::string& value) {
    auto tk = tokens(value);
    Amp A;
    if (!tk.empty() && tk[0] == "ramp") {
        if (tk.size() != 2)
            throw std::runtime_error(k + " : attendu 'ramp T' (T > 0 [s])");
        double T = num(k, tk[1]);
        if (!(T > 0.0))
            throw std::runtime_error(k + " : ramp T exige T > 0 [s]");
        A.smooth = true;
        A.t = {0.0, T};
        A.a = {0.0, 1.0};
    } else {
        if (tk.size() < 2 || tk.size() % 2 != 0)
            throw std::runtime_error(k + " : attendu 't0 a0 t1 a1 ...' "
                "(paires temps [s], facteur) ou 'ramp T'");
        for (std::size_t q = 0; q < tk.size(); q += 2) {
            A.t.push_back(num(k, tk[q]));
            A.a.push_back(num(k, tk[q + 1]));
            if (A.t.size() > 1 && !(A.t.back() > A.t[A.t.size() - 2]))
                throw std::runtime_error(k + " : les temps doivent etre "
                                         "strictement croissants");
        }
    }
    return A;
}

// facteur d'amplitude a l'instant t (constant hors de la table)
inline double ampAt(const std::vector<double>& at, const std::vector<double>& aa,
                    bool smooth, double t) {
    if (smooth) {
        const double T = at[1];
        if (t <= 0.0) return 0.0;
        if (t >= T) return 1.0;
        return 0.5 * (1.0 - std::cos(M_PI * t / T));
    }
    if (t <= at.front()) return aa.front();
    if (t >= at.back()) return aa.back();
    std::size_t j = 1;
    while (at[j] < t) ++j;
    const double s = (t - at[j - 1]) / (at[j] - at[j - 1]);
    return aa[j - 1] + s * (aa[j] - aa[j - 1]);
}

// fix.<g> = x y z | all (xy, "x z"... acceptes) -> masque de bits x=1 y=2 z=4
inline int parseFixMask(const std::string& k, const std::string& value) {
    int mask = 0;
    for (const auto& t : tokens(value)) {
        if (t == "all" || t == "xyz") mask |= 7;
        else if (t.find_first_not_of("xyz") == std::string::npos)
            for (char c : t) mask |= 1 << (c - 'x');
        else
            throw std::runtime_error(k + " : axes attendus parmi x y z "
                                     "all, lu '" + t + "'");
    }
    if (mask == 0)
        throw std::runtime_error(k + " : aucun axe (x y z | all)");
    return mask;
}

// velocity.<g> = vx vy vz ('free' = axe libre) -> masque et vitesse
inline int parseVelocity(const std::string& k, const std::string& value,
                         Eigen::Vector3d& v) {
    auto tk = tokens(value);
    if (tk.size() != 3)
        throw std::runtime_error(k + " : attendu 'vx vy vz' [m/s], 'free' "
                                 "pour un axe libre");
    int mask = 0;
    v = Eigen::Vector3d::Zero();
    for (int a = 0; a < 3; ++a) {
        if (tk[a] == "free") continue;
        v(a) = num(k, tk[a]);
        mask |= 1 << a;
    }
    if (mask == 0)
        throw std::runtime_error(k + " : les trois axes sont 'free' — la "
                                 "cle n'imposerait rien");
    return mask;
}

}  // namespace grouploads
}  // namespace rockim
