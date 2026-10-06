// ---------------------------------------------------------------------------
// Fem3dLoads — chargements et conditions aux limites PAR GROUPES (fem3d).
//
// Portage du 2026-10-03 de src/Fdem3dLoads.cpp (meme jour) : memes cles, meme
// semantique, memes messages d'erreur, memes colonnes de history.csv et meme
// ligne de resume, pour qu'un deck soit interchangeable entre fdem3d et fem3d
// a la ligne `mode` pres. Les briques sans etat (nombres stricts, amplitudes,
// axes) sont COMMUNES : include/rockim/GroupLoads.hpp.
//
// Difference de structure : fem3d est a NOEUDS PARTAGES, sans copies ni
// insertion — un sommet du maillage = un noeud, une face exterieure = trois
// noeuds. Une face dont l'element est erode ne recoit plus de charge
// (meme regle que confiningForces).
//
// Principe VIII : toutes les cles sont opt-in. Sans aucune d'elles,
// setupUserLoads() rend la main avant d'ecrire quoi que ce soit, uOn_ reste
// faux, et aucune autre instruction du solveur ne change (bit-identite).
// ---------------------------------------------------------------------------
#include "rockim/Fem3dSolver.hpp"
#include "rockim/GroupLoads.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>

namespace rockim {

namespace {

using grouploads::nums;

std::array<int, 3> sortedKey(const std::array<int, 3>& n) {
    std::array<int, 3> k = n;
    std::sort(k.begin(), k.end());
    return k;
}

}  // namespace

double Fem3dSolver::userAmp(int k, double t) const {
    const UAmp& A = uAmp_[k];
    return grouploads::ampAt(A.t, A.a, A.smooth, t);
}

// puissance des forces nodales courantes sur les vitesses du pas (v_ vit en
// t - dt/2) : sert a ventiler le bilan (elements, puis outil) — uOn_ seul
double Fem3dSolver::userPower() const {
    double p = 0.0;
    for (std::size_t i = 0; i < X0_.size(); ++i) p += f_[i].dot(v_[i]);
    return p;
}

// Integration d'un noeud LIBRE quand les cles par groupes sont actives.
// Axes libres : exactement les formules de Fem3dSolver::integrate (ressort de
// champ lointain relatif a uConf_, Cundall, mise a jour de v, Lysmer par
// division implicite, locWork_, lysWork_, plans de symetrie). Axes imposes
// (fix./velocity.) : ni ressort, ni Cundall, ni Lysmer ; v = vitesse imposee
// a mi-pas (t + dt/2, la vitesse du schema saute-mouton), reaction
// R = m (v_imp - v)/dt - f (meme formule que fdem3d), travail de liaison
// R (v_imp + v)/2 dt. Comptabilite du bilan en plus, sans toucher aux
// valeurs calculees.
void Fem3dSolver::integrateUserNode(int i) {
    const unsigned mk = uMask_[i];
    const double tMid = t_ + 0.5 * dt_;
    const Eigen::Vector3d vOld = v_[i];
    double fLoc[3] = {0.0, 0.0, 0.0};
    for (int a = 0; a < 3; ++a) {
        if ((mk >> a) & 1u) continue;      // axe impose : ni ressort ni Cundall
        if (kAbs_[i](a) > 0) {
            const double fk = kAbs_[i](a) * (u_[i](a) - uConf_[i](a));
            f_[i](a) -= fk;
            uSprW_ -= fk * vOld(a) * dt_;
        }
        if (damping_ > 0) {
            fLoc[a] = -damping_ * std::abs(f_[i](a))
                      * (v_[i](a) > 0 ? 1.0 : (v_[i](a) < 0 ? -1.0 : 0.0));
            f_[i](a) += fLoc[a];
            uCundW_ += fLoc[a] * vOld(a) * dt_;
        }
    }
    for (int a = 0; a < 3; ++a) {
        if ((mk >> a) & 1u) {
            const int k = uAmpOf_[i][a];
            const double vt = k < 0 ? uVel_[i](a) : uVel_[i](a) * userAmp(k, tMid);
            const double vo = v_[i](a);
            const double Ri = m_[i] * (vt - vo) / dt_ - f_[i](a);
            uR_[i](a) = Ri;
            uBcW_ += Ri * 0.5 * (vt + vo) * dt_;
            // dKE = (R + f)(vt + vo)/2 dt : la part f vo dt est deja dans les
            // postes (puissance sur v_(n-1/2)), le reste est ici
            uBias_ += f_[i](a) * 0.5 * (vt - vo) * dt_;
            v_[i](a) = vt;
            continue;
        }
        uBias_ += f_[i](a) * f_[i](a) * dt_ * dt_ / (2.0 * m_[i]);
        v_[i](a) += (dt_ / m_[i]) * f_[i](a);
        if (cAbs_[i](a) > 0) {
            const double v1 = v_[i](a);
            v_[i](a) /= 1.0 + dt_ * cAbs_[i](a) / m_[i];
            lysWork_ += cAbs_[i](a) * v_[i](a) * v_[i](a) * dt_;
            uLysW_ -= 0.5 * m_[i] * (v1 * v1 - v_[i](a) * v_[i](a));
        }
        if (damping_ > 0) locWork_ -= fLoc[a] * dt_ * v_[i](a);
    }
    // plans de symetrie (refuses en scenario = loads, ou le bilan est
    // imprime ; conserves ailleurs pour ne rien changer au montage)
    if (symY_) v_[i].y() = 0.0;
    if (symQ_) {
        if (X0_[i].x() < 1e-9) v_[i].x() = 0.0;
        if (X0_[i].y() < 1e-9) v_[i].y() = 0.0;
    }
    u_[i] += dt_ * v_[i];
}

void Fem3dSolver::setupUserLoads() {
    std::map<std::string, std::vector<std::string>> keys;   // prefixe -> cles
    bool any = false;
    for (const char* p : grouploads::kPrefixes) {
        keys[p] = cfg_.keysWithPrefix(p);
        if (!keys[p].empty()) any = true;
    }
    if (!any) return;                      // aucune cle : rien ne change
    if (scen_ == Scenario::LOADS && (symY_ || symQ_))
        throw std::runtime_error("scenario = loads : symmetryY / quarterModel "
            "annulent des vitesses hors du bilan des liaisons — decrire le plan "
            "de symetrie par un groupe : fix.<g> = y (ou x)");

    // ---- inventaire des noms de groupes disponibles ----------------------
    const int nN = (int)X0_.size();
    std::map<std::array<int, 3>, std::vector<int>> extOf;   // faces exterieures
    for (int f = 0; f < (int)exterior_.size(); ++f)
        extOf[sortedKey(exterior_[f].n)].push_back(f);

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
        bool done = false;
        for (const auto& M : mshLow_) {
            if (M.name != nm) continue;
            G.dim = M.dim;
            G.nodes = M.verts;
            for (const auto& tri : M.tris) {
                auto e = extOf.find(sortedKey(tri));
                if (e == extOf.end()) { ++G.missing; continue; }
                if (e->second.size() != 1) {
                    G.missing += 1;          // face vue deux fois : ambigu
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
                for (const auto& e : el_) {
                    if (e.grain != g) continue;
                    for (int a = 0; a < 4; ++a) G.nodes.push_back(e.n[a]);
                }
                done = true;
            }
        if (!done && cfg_.has("point." + nm)) {
            const std::string pk = "point." + nm;
            auto x = nums(pk, cfg_.gets(pk, ""), 3);
            Eigen::Vector3d P(x[0], x[1], x[2]);
            int best = -1;
            double dBest = 1e300;
            for (int i = 0; i < nN; ++i) {
                if (!(m_[i] > 0.0)) continue;    // noeud sans matiere
                double d = (X0_[i] - P).squaredNorm();
                if (d < dBest) { dBest = d; best = i; }
            }
            G.dim = 0;
            G.nodes = {best};
            std::cout << "[FEM3D] " << pk << " : sommet " << best << " a "
                      << std::sqrt(dBest) << " m du point demande ("
                      << X0_[best].transpose() << ", repere solveur)\n";
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
            for (const auto& bf : exterior_) {
                const auto& n = bf.n;
                if (inside(X0_[n[0]]) && inside(X0_[n[1]]) && inside(X0_[n[2]])) {
                    G.faces.push_back(bf);
                    for (int q = 0; q < 3; ++q) G.nodes.push_back(n[q]);
                }
            }
            if (G.faces.empty())
                throw std::runtime_error(bk + " : aucune face exterieure dans "
                    "la boite (repere SOLVEUR : le maillage est translate a "
                    "l'origine, boite [0, W] x [0, D] x [0, H])");
            done = true;
        }
        std::sort(G.nodes.begin(), G.nodes.end());
        G.nodes.erase(std::unique(G.nodes.begin(), G.nodes.end()),
                      G.nodes.end());
        if (G.nodes.empty())
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
        const grouploads::Amp P = grouploads::parseAmp(k, cfg_.gets(k, ""));
        UAmp A;
        A.t = P.t;
        A.a = P.a;
        A.smooth = P.smooth;
        ampOf[nm] = (int)uAmp_.size();
        uAmp_.push_back(A);
    }

    // ---- conditions aux limites : fix. et velocity. ------------------------
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
        for (int i : uGrp_[g].nodes) {
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
        const int mask = grouploads::parseFixMask(k, cfg_.gets(k, ""));
        int g = groupFor(k, nm);
        addBc(k, g, mask, Eigen::Vector3d::Zero());
    }
    for (const auto& k : keys["velocity."]) {
        const std::string nm = k.substr(9);
        Eigen::Vector3d v;
        const int mask = grouploads::parseVelocity(k, cfg_.gets(k, ""), v);
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
            for (int i : G.nodes) M += m_[i];
            for (int i : G.nodes) L.nodeW.push_back({i, m_[i] / M});
        } else {                           // points, courbes : parts egales
            const double share = 1.0 / (double)G.nodes.size();
            for (int i : G.nodes) L.nodeW.push_back({i, share});
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
    uKe0_ = 0.0;
    for (int i = 0; i < nN; ++i) uKe0_ += 0.5 * m_[i] * v_[i].squaredNorm();

    // ---- compte rendu -----------------------------------------------------
    static const char* const kDim[] = {"point", "courbe", "surface", "corps"};
    std::cout << "[FEM3D] charges et CL par groupes : " << uGrp_.size()
              << " groupe(s), " << uBc_.size() << " liaison(s), "
              << uLoad_.size() << " charge(s)";
    if (meshOrigin_.squaredNorm() > 0.0)
        std::cout << " — repere solveur = fichier - ("
                  << meshOrigin_.transpose() << ")";
    std::cout << "\n";
    for (int g = 0; g < (int)uGrp_.size(); ++g) {
        const UGroup& G = uGrp_[g];
        std::cout << "[FEM3D]   " << G.name << " (" << kDim[G.dim] << ") : "
                  << G.nodes.size() << " noeuds";
        if (G.dim == 2) {
            double A = 0.0;
            for (const auto& f : G.faces) A += faceArea0(f);
            std::cout << ", " << G.faces.size() << " faces, aire " << A << " m2";
        }
        if (uGrpAmp_[g] >= 0) std::cout << ", amplitude." << G.name;
        std::cout << "\n";
    }
    if (scen_ != Scenario::LOADS)
        std::cout << "[FEM3D]   scenario avec montage cable : bilan d'energie "
                     "par postes imprime en scenario = loads seulement\n";
}

void Fem3dSolver::userLoadForces() {
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
            if (el_[bf.elem].st.eroded) continue;      // face decharge
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

void Fem3dSolver::userHistoryHeader(std::ostream& os) const {
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

void Fem3dSolver::userHistoryRow(std::ostream& os) const {
    for (int g = 0; g < (int)uGrp_.size(); ++g) {
        const UGroup& G = uGrp_[g];
        Eigen::Vector3d U = Eigen::Vector3d::Zero();
        for (int i : G.nodes) U += u_[i];
        U /= (double)G.nodes.size();
        os << "," << U.x() << "," << U.y() << "," << U.z();
        for (const auto& B : uBc_)
            if (B.grp == g) {
                // reaction du pas : sur les axes tenus par CE groupe
                Eigen::Vector3d R = Eigen::Vector3d::Zero();
                for (int i : G.nodes)
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
void Fem3dSolver::userSummary() const {
    for (int g = 0; g < (int)uGrp_.size(); ++g) {
        const UGroup& G = uGrp_[g];
        Eigen::Vector3d U = Eigen::Vector3d::Zero();
        for (int i : G.nodes) U += u_[i];
        U /= (double)G.nodes.size();
        std::cout << "[FEM3D]   groupe " << G.name << " : U = " << U.x()
                  << " " << U.y() << " " << U.z() << " m";
        for (const auto& B : uBc_)
            if (B.grp == g) {
                Eigen::Vector3d R = Eigen::Vector3d::Zero();
                for (int i : G.nodes)
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

// Bilan d'energie par postes (scenario = loads). Tous les travaux sont
// cumules au sens f . v_(n-1/2) dt, si bien que, pas par pas et noeud par
// noeud, dKE = somme des postes + correction saute-mouton (f^2 dt^2/2m sur
// les axes libres, f (v_imp - v)/2 dt sur les axes imposes) : le residu est
// un zero d'arrondi quand la comptabilite est complete. Lysmer est compte
// EXACTEMENT (variation d'energie cinetique de la division implicite).
// Hors scenario = loads, seule la ligne charges/liaisons est imprimee.
void Fem3dSolver::userEnergySummary() const {
    std::cout << "[FEM3D] charges et CL par groupes :\n";
    if (scen_ != Scenario::LOADS) {
        std::cout << "[FEM3D]   charges      : " << uLoadW_
                  << " J (force./traction./pressure. -> solide), "
                     "liaisons " << uBcW_
                  << " J (fix./velocity. -> solide)\n";
        userSummary();
        return;
    }
    double ke = 0.0;
    for (std::size_t i = 0; i < X0_.size(); ++i)
        ke += 0.5 * m_[i] * v_[i].squaredNorm();
    const double sumW = uElW_ + uToolW_ + confWork_ + uLoadW_ + uBcW_
                      + uSprW_ + uCundW_ + uLysW_ + uBias_;
    const double resid = (ke - uKe0_) - sumW;
    const double gross = std::abs(uElW_) + std::abs(uToolW_)
                       + std::abs(confWork_) + std::abs(uLoadW_)
                       + std::abs(uBcW_) + std::abs(uSprW_)
                       + std::abs(uCundW_) + std::abs(uLysW_);
    const double scale = std::max({uKe0_, ke, gross, 1e-30});
    const bool zeroCase = scale < 1e-12;
    std::cout << "[FEM3D] energy budget (charges par groupes) : KE " << uKe0_
              << " -> " << ke << " J\n"
              << "[FEM3D]   elements     : " << -uElW_ << " J preleve\n"
              << "[FEM3D]   Cundall      : " << -uCundW_ << " J\n"
              << "[FEM3D]   frontieres   : " << -uLysW_
              << " J (Lysmer), ressorts " << -uSprW_ << " J\n";
    if (confP_ > 0.0 || topP_ > 0.0 || bottomP_ > 0.0)
        std::cout << "[FEM3D]   confinement  : " << confWork_
                  << " J (pression suiveuse -> solide)\n";
    std::cout << "[FEM3D]   charges      : " << uLoadW_
              << " J (force./traction./pressure. -> solide), "
                 "liaisons " << uBcW_
              << " J (fix./velocity. -> solide)\n";
    userSummary();
    std::cout << "[FEM3D]   integration  : +" << uBias_
              << " J (correction leapfrog f^2 dt^2/2m)\n"
              << "[FEM3D]   residu       : " << resid << " J ("
              << 100.0 * std::abs(resid) / scale << " % de l'echelle) ["
              << (zeroCase ? "OK (zero machine)"
                  : std::abs(resid) <= 0.01 * scale ? "OK" : "CHECK")
              << "]\n";
}

}  // namespace rockim
