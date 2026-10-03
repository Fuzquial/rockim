// ---------------------------------------------------------------------------
// Fdem3dLoads — chargements et conditions aux limites PAR GROUPES (fdem3d).
//
// Ajout du 2026-10-03. Le pendant des tables /YD/YDB/ du code public de
// Solidity (Y3Drd.c l. 1221-1260, Y3Dsd.c l. 80-102 et 238 : vitesse imposee
// par axe, force nodale, acceleration, rattachees a chaque noeud par un
// numero de propriete ecrit par le preprocesseur), mais adresse par GROUPES
// NOMMES et module par une amplitude en temps. Cles et regles : bloc
// « Chargements et conditions aux limites PAR GROUPES » de Fdem3dSolver.hpp
// et DOCUMENTATION_rockim.md §5.21.
//
// Principe VIII : toutes les cles sont opt-in. Sans aucune d'elles,
// setupUserLoads() rend la main avant d'ecrire quoi que ce soit, uOn_ reste
// faux, et aucune autre instruction du solveur ne change (bit-identite).
// ---------------------------------------------------------------------------
#include "rockim/Fdem3dSolver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace rockim {

namespace {

std::vector<std::string> tokens(const std::string& s) {
    std::istringstream iss(s);
    std::vector<std::string> out;
    std::string t;
    while (iss >> t) out.push_back(t);
    return out;
}

// nombre strict : la virgule decimale et les suffixes parasites sont refuses
// (meme regle que Config::getd), avec le nom de la cle dans le message
double num(const std::string& key, const std::string& t) {
    std::size_t pos = 0;
    double v = 0.0;
    bool ok = true;
    try { v = std::stod(t, &pos); } catch (...) { ok = false; }
    if (!ok || pos != t.size() || !std::isfinite(v))
        throw std::runtime_error(key + " : valeur numerique invalide '" + t
                                 + "' (point decimal obligatoire)");
    return v;
}

std::vector<double> nums(const std::string& key, const std::string& s,
                         std::size_t n) {
    auto tk = tokens(s);
    if (tk.size() != n)
        throw std::runtime_error(key + " : attendu " + std::to_string(n)
                                 + " valeurs, lu '" + s + "'");
    std::vector<double> v;
    for (const auto& t : tk) v.push_back(num(key, t));
    return v;
}

std::array<int, 3> faceKey(const std::vector<int>& vOf,
                           const std::array<int, 3>& n) {
    std::array<int, 3> k = {vOf[n[0]], vOf[n[1]], vOf[n[2]]};
    std::sort(k.begin(), k.end());
    return k;
}

}  // namespace

double Fdem3dSolver::userAmp(int k, double t) const {
    const UAmp& A = uAmp_[k];
    if (A.smooth) {
        const double T = A.t[1];
        if (t <= 0.0) return 0.0;
        if (t >= T) return 1.0;
        return 0.5 * (1.0 - std::cos(M_PI * t / T));
    }
    if (t <= A.t.front()) return A.a.front();
    if (t >= A.t.back()) return A.a.back();
    std::size_t j = 1;
    while (A.t[j] < t) ++j;
    const double s = (t - A.t[j - 1]) / (A.t[j] - A.t[j - 1]);
    return A.a[j - 1] + s * (A.a[j] - A.a[j - 1]);
}

void Fdem3dSolver::setupUserLoads() {
    static const char* const kPre[] = {"fix.", "velocity.", "traction.",
                                       "pressure.", "force.", "amplitude.",
                                       "point.", "box."};
    std::map<std::string, std::vector<std::string>> keys;   // prefixe -> cles
    bool any = false;
    for (const char* p : kPre) {
        keys[p] = cfg_.keysWithPrefix(p);
        if (!keys[p].empty()) any = true;
    }
    if (!any) return;                      // aucune cle : rien ne change
    if (jbOn_)
        throw std::runtime_error("scenario = jointbench construit son propre "
            "montage : les cles fix./velocity./force./... y sont refusees");

    // ---- inventaire des noms de groupes disponibles ----------------------
    // sommet virtuel -> copies, et position de reference d'un sommet
    int nV = 0;
    for (int v : vOf_) nV = std::max(nV, v + 1);
    std::vector<std::vector<int>> copiesV(nV);
    for (int i = 0; i < (int)X0_.size(); ++i) copiesV[vOf_[i]].push_back(i);
    // faces exterieures d'origine par triplet de sommets virtuels
    std::map<std::array<int, 3>, std::vector<int>> extOf;
    for (int f = 0; f < (int)exterior_.size(); ++f)
        extOf[faceKey(vOf_, exterior_[f].n)].push_back(f);

    std::map<std::string, std::string> source;   // nom -> origine (messages)
    auto declare = [&](const std::string& nm, const std::string& src) {
        auto it = source.find(nm);
        if (it != source.end())
            throw std::runtime_error("groupe '" + nm + "' defini deux fois ("
                + it->second + " et " + src + ") : renommer l'un des deux");
        source[nm] = src;
    };
    for (const auto& G : mshLow_)
        declare(G.name, "groupe physique Gmsh de dimension "
                        + std::to_string(G.dim));
    if (!elemGroup_.empty() || !groupName_.empty())
        for (const auto& nm : groupName_) declare(nm, "corps (volume physique)");
    for (const auto& k : keys["point."]) declare(k.substr(6), k);
    for (const auto& k : keys["box."]) declare(k.substr(4), k);

    std::map<std::string, int> gIdx;       // nom -> uGrp_ (construit a la demande)
    auto group = [&](const std::string& key, const std::string& nm) -> int {
        auto it = gIdx.find(nm);
        if (it != gIdx.end()) return it->second;
        if (!source.count(nm)) {
            std::string avail;
            for (const auto& s : source) avail += " " + s.first;
            throw std::runtime_error(key + " : groupe '" + nm + "' inconnu. "
                "Groupes disponibles :" + (avail.empty() ? " (aucun)" : avail)
                + ". Les groupes viennent des $PhysicalNames du maillage "
                  "(dim 0, 1, 2 ou volumes), de point.<g> ou de box.<g>.");
        }
        UGroup G;
        G.name = nm;
        auto addFace = [&](int f) {
            G.faces.push_back(exterior_[f]);
            for (int q = 0; q < 3; ++q) G.verts.push_back(vOf_[exterior_[f].n[q]]);
        };
        bool done = false;
        for (const auto& M : mshLow_) {
            if (M.name != nm) continue;
            G.dim = M.dim;
            G.verts = M.verts;
            for (const auto& tri : M.tris) {
                std::array<int, 3> k = tri;
                std::sort(k.begin(), k.end());
                auto e = extOf.find(k);
                if (e == extOf.end()) { ++G.missing; continue; }
                if (e->second.size() != 1) {
                    G.missing += 1;          // interface de deux corps : ambigu
                    continue;
                }
                G.faces.push_back(exterior_[e->second[0]]);
            }
            done = true;
        }
        if (!done)
            for (int g = 0; g < (int)groupName_.size() && !done; ++g) {
                if (groupName_[g] != nm) continue;
                G.dim = 3;
                for (std::size_t e = 0; e < el_.size(); ++e) {
                    if (!elemGroup_.empty() && elemGroup_[e] != g) continue;
                    for (int a = 0; a < 4; ++a) G.copies.push_back(el_[e].n[a]);
                }
                for (int i : G.copies) G.verts.push_back(vOf_[i]);
                done = true;
            }
        if (!done && cfg_.has("point." + nm)) {
            const std::string pk = "point." + nm;
            auto x = nums(pk, cfg_.gets(pk, ""), 3);
            Eigen::Vector3d P(x[0], x[1], x[2]);
            int best = -1;
            double dBest = 1e300;
            for (int v = 0; v < nV; ++v) {
                if (copiesV[v].empty()) continue;
                double d = (X0_[copiesV[v][0]] - P).squaredNorm();
                if (d < dBest) { dBest = d; best = v; }
            }
            G.dim = 0;
            G.verts = {best};
            std::cout << "[FDEM3D] " << pk << " : sommet " << best << " a "
                      << std::sqrt(dBest) << " m du point demande ("
                      << X0_[copiesV[best][0]].transpose()
                      << ", repere solveur)\n";
            done = true;
        }
        if (!done && cfg_.has("box." + nm)) {
            const std::string bk = "box." + nm;
            auto b = nums(bk, cfg_.gets(bk, ""), 6);
            Eigen::Vector3d lo(std::min(b[0], b[3]), std::min(b[1], b[4]),
                               std::min(b[2], b[5]));
            Eigen::Vector3d hi(std::max(b[0], b[3]), std::max(b[1], b[4]),
                               std::max(b[2], b[5]));
            const double tol = 1e-9 * std::max({W_, D_, H_, 1e-12});
            auto inside = [&](const Eigen::Vector3d& X) {
                return (X.array() >= lo.array() - tol).all()
                    && (X.array() <= hi.array() + tol).all();
            };
            G.dim = 2;
            for (int f = 0; f < (int)exterior_.size(); ++f) {
                const auto& n = exterior_[f].n;
                if (inside(X0_[n[0]]) && inside(X0_[n[1]]) && inside(X0_[n[2]]))
                    addFace(f);
            }
            if (G.faces.empty())
                throw std::runtime_error(bk + " : aucune face exterieure dans "
                    "la boite (repere SOLVEUR : le maillage est translate a "
                    "l'origine, boite [0, W] x [0, D] x [0, H])");
            done = true;
        }
        std::sort(G.verts.begin(), G.verts.end());
        G.verts.erase(std::unique(G.verts.begin(), G.verts.end()), G.verts.end());
        if (G.copies.empty())
            for (int v : G.verts)
                for (int i : copiesV[v]) G.copies.push_back(i);
        std::sort(G.copies.begin(), G.copies.end());
        G.copies.erase(std::unique(G.copies.begin(), G.copies.end()),
                       G.copies.end());
        if (G.verts.empty())
            throw std::runtime_error(key + " : le groupe '" + nm + "' n'a "
                "aucun sommet (groupe physique vide ?)");
        int id = (int)uGrp_.size();
        uGrp_.push_back(std::move(G));
        uGrpAmp_.push_back(-1);
        gIdx[nm] = id;
        return id;
    };

    // ---- amplitudes -------------------------------------------------------
    std::map<std::string, int> ampOf;
    for (const auto& k : keys["amplitude."]) {
        const std::string nm = k.substr(10);
        auto tk = tokens(cfg_.gets(k, ""));
        UAmp A;
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
        ampOf[nm] = (int)uAmp_.size();
        uAmp_.push_back(A);
    }

    // ---- conditions aux limites : fix. et velocity. ------------------------
    const std::size_t nN = X0_.size();
    uMask_.assign(nN, 0);
    uVel_.assign(nN, Eigen::Vector3d::Zero());
    uAmpOf_.assign(nN, {-1, -1, -1});
    uR_.assign(nN, Eigen::Vector3d::Zero());
    auto addBc = [&](const std::string& key, int g, int mask,
                     const Eigen::Vector3d& v) {
        UBc B;
        B.grp = g;
        B.mask = mask;
        B.v = v;
        B.amp = uGrpAmp_[g];
        long clash = 0;
        for (int i : uGrp_[g].copies) {
            if (flag_[i] != FREE)
                throw std::runtime_error(key + " : le groupe '" + uGrp_[g].name
                    + "' touche des noeuds deja tenus par le montage du "
                      "scenario (fond encastre en percussion/coupe, mors en "
                      "traction). Poser scenario = loads pour decrire tout le "
                      "montage par groupes.");
            for (int a = 0; a < 3; ++a) {
                if (!((mask >> a) & 1)) continue;
                if ((uMask_[i] >> a) & 1) {
                    if (uVel_[i](a) != v(a) || uAmpOf_[i][a] != B.amp) ++clash;
                    continue;
                }
                uMask_[i] |= (unsigned char)(1u << a);
                uVel_[i](a) = v(a);
                uAmpOf_[i][a] = B.amp;
            }
        }
        if (clash > 0)
            throw std::runtime_error(key + " : " + std::to_string(clash)
                + " ddl deja imposes a une AUTRE valeur par un autre groupe "
                  "(arete ou coin commun) — rendre les conditions compatibles");
        uBc_.push_back(B);
    };
    auto groupFor = [&](const std::string& key, const std::string& nm) {
        int g = group(key, nm);
        auto a = ampOf.find(nm);
        if (a != ampOf.end()) uGrpAmp_[g] = a->second;
        return g;
    };
    for (const auto& k : keys["fix."]) {
        const std::string nm = k.substr(4);
        if (cfg_.has("velocity." + nm))
            throw std::runtime_error(k + " et velocity." + nm + " portent sur "
                "le meme groupe : un seul des deux (fix = velocity 0)");
        int mask = 0;
        for (const auto& t : tokens(cfg_.gets(k, ""))) {
            if (t == "all" || t == "xyz") mask |= 7;
            else if (t.find_first_not_of("xyz") == std::string::npos)
                for (char c : t) mask |= 1 << (c - 'x');
            else
                throw std::runtime_error(k + " : axes attendus parmi x y z "
                                         "all, lu '" + t + "'");
        }
        if (mask == 0)
            throw std::runtime_error(k + " : aucun axe (x y z | all)");
        int g = groupFor(k, nm);
        addBc(k, g, mask, Eigen::Vector3d::Zero());
    }
    for (const auto& k : keys["velocity."]) {
        const std::string nm = k.substr(9);
        auto tk = tokens(cfg_.gets(k, ""));
        if (tk.size() != 3)
            throw std::runtime_error(k + " : attendu 'vx vy vz' [m/s], 'free' "
                                     "pour un axe libre");
        int mask = 0;
        Eigen::Vector3d v = Eigen::Vector3d::Zero();
        for (int a = 0; a < 3; ++a) {
            if (tk[a] == "free") continue;
            v(a) = num(k, tk[a]);
            mask |= 1 << a;
        }
        if (mask == 0)
            throw std::runtime_error(k + " : les trois axes sont 'free' — la "
                                     "cle n'imposerait rien");
        int g = groupFor(k, nm);
        addBc(k, g, mask, v);
    }

    // ---- charges : traction., pressure., force. ---------------------------
    auto faceArea0 = [&](const BFace& f) {
        const Eigen::Vector3d A = X0_[f.n[0]], B = X0_[f.n[1]], C = X0_[f.n[2]];
        return 0.5 * (B - A).cross(C - A).norm();
    };
    auto needFaces = [&](const std::string& key, int g) {
        const UGroup& G = uGrp_[g];
        if (G.dim != 2)
            throw std::runtime_error(key + " : '" + G.name + "' n'est pas une "
                "SURFACE (dimension " + std::to_string(G.dim) + ") — traction "
                "et pression exigent un groupe de dimension 2 ou box.<g> ; "
                "pour un point, une courbe ou un corps, utiliser force.<g>");
        if (G.missing > 0)
            throw std::runtime_error(key + " : " + std::to_string(G.missing)
                + " triangle(s) de '" + G.name + "' ne sont pas des faces "
                  "EXTERIEURES du maillage (face interne a un corps, ou "
                  "interface ambigue entre deux corps non lies)");
        if (G.faces.empty())
            throw std::runtime_error(key + " : '" + G.name + "' n'a aucun "
                                     "triangle");
    };
    for (const auto& k : keys["traction."]) {
        const std::string nm = k.substr(9);
        auto x = nums(k, cfg_.gets(k, ""), 3);
        ULoad L;
        L.grp = groupFor(k, nm);
        needFaces(k, L.grp);
        L.kind = 0;
        L.val = {x[0], x[1], x[2]};
        L.amp = uGrpAmp_[L.grp];
        for (const auto& f : uGrp_[L.grp].faces) L.area0 += faceArea0(f);
        uLoad_.push_back(L);
    }
    for (const auto& k : keys["pressure."]) {
        const std::string nm = k.substr(9);
        auto x = nums(k, cfg_.gets(k, ""), 1);
        ULoad L;
        L.grp = groupFor(k, nm);
        needFaces(k, L.grp);
        L.kind = 1;
        L.val = {x[0], 0.0, 0.0};
        L.amp = uGrpAmp_[L.grp];
        for (const auto& f : uGrp_[L.grp].faces) L.area0 += faceArea0(f);
        uLoad_.push_back(L);
    }
    for (const auto& k : keys["force."]) {
        const std::string nm = k.substr(6);
        auto x = nums(k, cfg_.gets(k, ""), 3);
        ULoad L;
        L.grp = groupFor(k, nm);
        L.kind = 2;
        L.val = {x[0], x[1], x[2]};
        L.amp = uGrpAmp_[L.grp];
        const UGroup& G = uGrp_[L.grp];
        if (G.dim == 2) {
            needFaces(k, L.grp);           // repartie au prorata des aires
            for (const auto& f : G.faces) L.area0 += faceArea0(f);
        } else if (G.dim == 3) {           // corps : au prorata des masses
            double M = 0.0;
            for (int i : G.copies) M += m_[i];
            for (int i : G.copies) L.nodeW.push_back({i, m_[i] / M});
        } else {                           // points, courbes : parts egales
            const double share = 1.0 / (double)G.verts.size();
            for (int v : G.verts) {        // puis masses des copies du sommet
                double M = 0.0;
                for (int i : copiesV[v]) M += m_[i];
                for (int i : copiesV[v])
                    L.nodeW.push_back({i, share * m_[i] / M});
            }
        }
        uLoad_.push_back(L);
    }
    // amplitude.<g> qui ne sert a rien : refusee (piege maison n. 1)
    for (const auto& [nm, k] : ampOf) {
        (void)k;
        if (!gIdx.count(nm))
            throw std::runtime_error("amplitude." + nm + " : aucune cle "
                "fix./velocity./traction./pressure./force. sur ce groupe");
    }
    for (const auto& k : keys["point."])
        if (!gIdx.count(k.substr(6)))
            throw std::runtime_error(k + " : groupe defini mais jamais charge "
                                     "ni tenu");
    for (const auto& k : keys["box."])
        if (!gIdx.count(k.substr(4)))
            throw std::runtime_error(k + " : groupe defini mais jamais charge "
                                     "ni tenu");
    uOn_ = true;

    // ---- compte rendu -----------------------------------------------------
    static const char* const kDim[] = {"point", "courbe", "surface", "corps"};
    std::cout << "[FDEM3D] charges et CL par groupes : " << uGrp_.size()
              << " groupe(s), " << uBc_.size() << " liaison(s), "
              << uLoad_.size() << " charge(s)";
    if (meshOrigin_.squaredNorm() > 0.0)
        std::cout << " — repere solveur = fichier - ("
                  << meshOrigin_.transpose() << ")";
    std::cout << "\n";
    for (int g = 0; g < (int)uGrp_.size(); ++g) {
        const UGroup& G = uGrp_[g];
        std::cout << "[FDEM3D]   " << G.name << " (" << kDim[G.dim] << ") : "
                  << G.verts.size() << " sommets, " << G.copies.size()
                  << " copies";
        if (G.dim == 2) {
            double A = 0.0;
            for (const auto& f : G.faces) A += faceArea0(f);
            std::cout << ", " << G.faces.size() << " faces, aire " << A << " m2";
        }
        if (uGrpAmp_[g] >= 0) std::cout << ", amplitude." << G.name;
        std::cout << "\n";
    }
    if (adaptive_)
        std::cout << "[FDEM3D]   insertion adaptative : traction/pression sur "
                     "les copies du tetra de surface, force de sommet au "
                     "prorata des masses des copies, vitesse imposee a toutes "
                     "les copies\n";
}

void Fdem3dSolver::userLoadForces() {
    double w = 0.0;
    for (auto& L : uLoad_) {
        const double a = L.amp < 0 ? 1.0 : userAmp(L.amp, t_);
        L.Fnow.setZero();
        if (a == 0.0) continue;
        const UGroup& G = uGrp_[L.grp];
        if (L.kind == 2 && G.dim != 2) {
            const Eigen::Vector3d F = a * L.val;
            for (const auto& [i, wi] : L.nodeW) {
                const Eigen::Vector3d Fi = wi * F;
                f_[i] += Fi;
                w += Fi.dot(v_[i]);
                L.Fnow += Fi;
            }
            continue;
        }
        for (const auto& bf : G.faces) {
            Eigen::Vector3d F;             // par noeud
            if (L.kind == 1) {             // pression suiveuse, aire courante
                const Eigen::Vector3d A = X0_[bf.n[0]] + u_[bf.n[0]];
                const Eigen::Vector3d B = X0_[bf.n[1]] + u_[bf.n[1]];
                const Eigen::Vector3d C = X0_[bf.n[2]] + u_[bf.n[2]];
                F = -a * L.val.x() * 0.5 * (B - A).cross(C - A) / 3.0;
            } else {                       // charge morte, aire de reference
                const Eigen::Vector3d A = X0_[bf.n[0]], B = X0_[bf.n[1]],
                                      C = X0_[bf.n[2]];
                const double Af = 0.5 * (B - A).cross(C - A).norm();
                F = L.kind == 0 ? Eigen::Vector3d(a * L.val * Af / 3.0)
                                : Eigen::Vector3d(a * L.val * (Af / L.area0) / 3.0);
            }
            for (int q = 0; q < 3; ++q) {
                f_[bf.n[q]] += F;
                w += F.dot(v_[bf.n[q]]);
            }
            L.Fnow += 3.0 * F;
        }
    }
    uLoadW_ += dt_ * w;                    // meme convention que confWork_
}

void Fdem3dSolver::userHistoryHeader(std::ostream& os) const {
    for (int g = 0; g < (int)uGrp_.size(); ++g) {
        const std::string& n = uGrp_[g].name;
        os << ",U_" << n << "_x,U_" << n << "_y,U_" << n << "_z";
        for (const auto& B : uBc_)
            if (B.grp == g) {
                os << ",RF_" << n << "_x,RF_" << n << "_y,RF_" << n << "_z";
                break;
            }
        for (const auto& L : uLoad_)
            if (L.grp == g) {
                os << ",F_" << n << "_x,F_" << n << "_y,F_" << n << "_z";
                break;
            }
    }
    os << ",eLoad,eBc";
}

void Fdem3dSolver::userHistoryRow(std::ostream& os) const {
    for (int g = 0; g < (int)uGrp_.size(); ++g) {
        const UGroup& G = uGrp_[g];
        Eigen::Vector3d U = Eigen::Vector3d::Zero();
        for (int i : G.copies) U += u_[i];
        U /= (double)G.copies.size();
        os << "," << U.x() << "," << U.y() << "," << U.z();
        for (const auto& B : uBc_)
            if (B.grp == g) {
                // reaction du pas : sur les axes tenus par CE groupe
                Eigen::Vector3d R = Eigen::Vector3d::Zero();
                for (int i : G.copies)
                    for (int a = 0; a < 3; ++a)
                        if ((B.mask >> a) & 1) R(a) += uR_[i](a);
                os << "," << R.x() << "," << R.y() << "," << R.z();
                break;
            }
        Eigen::Vector3d F = Eigen::Vector3d::Zero();
        bool has = false;
        for (const auto& L : uLoad_)
            if (L.grp == g) { F += L.Fnow; has = true; }
        if (has) os << "," << F.x() << "," << F.y() << "," << F.z();
    }
    os << "," << uLoadW_ << "," << uBcW_;
}

// fin de run : etat final de chaque groupe (lu par tools/verify_suite.py)
void Fdem3dSolver::userSummary() const {
    for (int g = 0; g < (int)uGrp_.size(); ++g) {
        const UGroup& G = uGrp_[g];
        Eigen::Vector3d U = Eigen::Vector3d::Zero();
        for (int i : G.copies) U += u_[i];
        U /= (double)G.copies.size();
        std::cout << "[FDEM3D]   groupe " << G.name << " : U = " << U.x()
                  << " " << U.y() << " " << U.z() << " m";
        for (const auto& B : uBc_)
            if (B.grp == g) {
                Eigen::Vector3d R = Eigen::Vector3d::Zero();
                for (int i : G.copies)
                    for (int a = 0; a < 3; ++a)
                        if ((B.mask >> a) & 1) R(a) += uR_[i](a);
                std::cout << " ; RF = " << R.x() << " " << R.y() << " "
                          << R.z() << " N";
                break;
            }
        Eigen::Vector3d F = Eigen::Vector3d::Zero();
        bool has = false;
        for (const auto& L : uLoad_)
            if (L.grp == g) { F += L.Fnow; has = true; }
        if (has)
            std::cout << " ; F = " << F.x() << " " << F.y() << " " << F.z()
                      << " N";
        std::cout << "\n";
    }
}

}  // namespace rockim
