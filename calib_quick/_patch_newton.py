import io, os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
H="include/rockim/Tessellation.hpp"; C="src/Tessellation.cpp"
h=io.open(H,encoding="utf-8").read(); c=io.open(C,encoding="utf-8").read()
def rep(s,a,b,name):
    assert s.count(a)==1,(name,s.count(a)); return s.replace(a,b,1)
# ---- voronoiCells : sortie optionnelle du Laplacien (voisin par arete) ----
c=rep(c,"""std::vector<std::vector<Vec2>> voronoiCells(const std::vector<Vec2>& seeds,
                                            double W, double H,
                                            const std::vector<double>* wgt = nullptr) {
    const int N = (int)seeds.size();
    const bool lag = (wgt != nullptr);
""","""// lap (optionnel, mode Laguerre) : pour chaque cellule i, la liste des
// (j, l_ij / (2 d_ij)) de ses aretes partagees — les coefficients du
// Laplacien dA_i/dw_j. Le voisin d une arete est identifie par le test de
// puissance a son milieu m : |m-p_j|^2 - w_j = |m-p_i|^2 - w_i ; un mur n a
// pas d egal.
std::vector<std::vector<Vec2>> voronoiCells(const std::vector<Vec2>& seeds,
                                            double W, double H,
                                            const std::vector<double>* wgt = nullptr,
                                            std::vector<std::vector<std::pair<int, double>>>* lap = nullptr) {
    const int N = (int)seeds.size();
    const bool lag = (wgt != nullptr);
    if (lag && lap) lap->assign((std::size_t)N, {});
""","lap sig")
c=rep(c,"""        cells[i] = std::move(poly);
    }
    return cells;
}
""","""        if (lag && lap && poly.size() >= 3) {
            std::vector<std::pair<int, double>>& adj = (*lap)[i];
            const double tol = 1e-10 * (W * W + H * H);
            for (std::size_t k = 0; k < poly.size(); ++k) {
                const Vec2& P = poly[k];
                const Vec2& Q = poly[(k + 1) % poly.size()];
                const Vec2 m = 0.5 * (P + Q);
                const double ell = (Q - P).norm();
                const double powI = (m - seeds[i]).squaredNorm() - (*wgt)[i];
                int best = -1;
                double bestPow = 1e300;
                for (const auto& cd : cand) {
                    const int j = cd.second;
                    const double pw = (m - seeds[j]).squaredNorm() - (*wgt)[j];
                    if (pw < bestPow) { bestPow = pw; best = j; }
                }
                if (best >= 0 && std::abs(bestPow - powI) < tol) {
                    const double dij = (seeds[best] - seeds[i]).norm();
                    adj.push_back({best, ell / (2.0 * dij)});
                }
            }
        }
        cells[i] = std::move(poly);
    }
    return cells;
}
""","lap fill")
# ---- build : poids a aires prescrites par Newton amorti ----
start = "    // ---- 2. Voronoi cells (+ Lloyd relaxation) -------------------------------\n"
end = "        if (it >= lloyd) break;\n"
i0 = c.index(start); i1 = c.index(end, i0) + len(end)
assert c.count(start) == 1
new = '''    // ---- 2. Voronoi cells (+ Lloyd relaxation) -------------------------------
    // POLYDISPERSITE : diagramme de Laguerre a AIRES PRESCRITES. Pour des
    // graines fixees il existe des poids, uniques a une constante pres,
    // realisant exactement toute famille d aires positives de somme W H
    // (Aurenhammer, Hoffmann & Aronov 1998) ; on les obtient par Newton
    // amorti sur le dual semi-discret (Kitagawa, Merigot & Thibert 2019),
    // depart w = 0 (Voronoi : aucune cellule vide), amortissement gardant
    // min A >= eps0 — c est la construction de Bourne, Kok, Roper & Spanjer
    // 2020 pour des grains de volumes donnes. Aires cibles A_i ~ s_i^2.
    // (Des poids FIXES w_i = (kappa s_i)^2, essayes d abord, ne rendent que
    // 40 % de la dispersion : les interstices du Poisson-disc se partagent
    // au perimetre et nivellent les aires — replica calib_quick/_lag_*.py.)
    //   dA_i/dw_i = sum_j l_ij / (2 d_ij),  dA_i/dw_j = -l_ij / (2 d_ij)
    //   (Laplacien L, singulier sur les constantes -> L + 11^T/N, resolu par
    //   gradient conjugue) ; Newton : L delta = A_cible - A.
    std::vector<double> wgt, aTgt;
    std::vector<std::vector<std::pair<int, double>>> lap;
    int newtonIt = 0;
    double newtonErr = 0.0;
    auto polyArea = [](const std::vector<Vec2>& p) {
        double a = 0.0;
        for (std::size_t k = 0; k < p.size(); ++k) {
            const Vec2& P = p[k];
            const Vec2& Q = p[(k + 1) % p.size()];
            a += P.x() * Q.y() - Q.x() * P.y();
        }
        return 0.5 * a;
    };
    auto setTargets = [&]() {
        double sum2 = 0.0;
        for (double v : sp) sum2 += v * v;
        aTgt.resize(seeds.size());
        for (std::size_t i = 0; i < seeds.size(); ++i)
            aTgt[i] = W * H * sp[i] * sp[i] / sum2;
    };
    auto solveWeights = [&]() {
        const std::size_t N = seeds.size();
        std::vector<std::vector<Vec2>> cl = voronoiCells(seeds, W, H, &wgt, &lap);
        std::vector<double> A(N);
        for (std::size_t i = 0; i < N; ++i) A[i] = cl[i].size() >= 3 ? polyArea(cl[i]) : 0.0;
        double eps0 = 1e300;
        for (std::size_t i = 0; i < N; ++i) eps0 = std::min(eps0, std::min(A[i], aTgt[i]));
        eps0 *= 0.5;
        std::vector<double> g(N), delta(N), r(N), p(N), q(N), w2(N), A2(N);
        auto dot = [](const std::vector<double>& a, const std::vector<double>& b) {
            double s = 0.0;
            for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
            return s;
        };
        auto matvec = [&](const std::vector<double>& x, std::vector<double>& y) {
            double mean = 0.0;
            for (double v : x) mean += v;
            mean /= (double)N;
            std::fill(y.begin(), y.end(), mean);          // + 11^T x / N
            for (std::size_t i = 0; i < N; ++i)            // symetrise : c/2 des deux cotes
                for (const auto& e : lap[i]) {
                    const double h = 0.5 * e.second * (x[i] - x[(std::size_t)e.first]);
                    y[i] += h;
                    y[(std::size_t)e.first] -= h;
                }
        };
        for (int k = 0; k < 40; ++k) {
            double err = 0.0;
            for (std::size_t i = 0; i < N; ++i) {
                g[i] = aTgt[i] - A[i];
                err = std::max(err, std::abs(g[i]) / aTgt[i]);
            }
            newtonIt = k;
            newtonErr = err;
            if (err < 1e-3) break;
            std::fill(delta.begin(), delta.end(), 0.0);
            r = g;
            p = r;
            double rr = dot(r, r);
            const double gg = rr;
            for (int cg = 0; cg < 1000 && rr > 1e-24 * gg; ++cg) {
                matvec(p, q);
                const double alpha = rr / dot(p, q);
                for (std::size_t i = 0; i < N; ++i) { delta[i] += alpha * p[i]; r[i] -= alpha * q[i]; }
                const double rr2 = dot(r, r);
                for (std::size_t i = 0; i < N; ++i) p[i] = r[i] + (rr2 / rr) * p[i];
                rr = rr2;
            }
            double tau = 1.0;
            std::vector<std::vector<Vec2>> cl2;
            for (;;) {
                for (std::size_t i = 0; i < N; ++i) w2[i] = wgt[i] + tau * delta[i];
                cl2 = voronoiCells(seeds, W, H, &w2, &lap);
                double amin = 1e300;
                for (std::size_t i = 0; i < N; ++i) {
                    A2[i] = cl2[i].size() >= 3 ? polyArea(cl2[i]) : 0.0;
                    amin = std::min(amin, A2[i]);
                }
                if (amin >= eps0 || tau < 1e-4) break;
                tau *= 0.5;
            }
            wgt.swap(w2);
            A.swap(A2);
            cl.swap(cl2);
        }
    };
    if (sizeSpread > 0.0) {
        setTargets();
        wgt.assign(seeds.size(), 0.0);
    }
    int nRedundant = 0;
    std::vector<std::vector<Vec2>> cells;
    for (int it = 0;; ++it) {
        if (!wgt.empty()) solveWeights();
        cells = voronoiCells(seeds, W, H, wgt.empty() ? nullptr : &wgt);
        if (!wgt.empty()) {
            // securite : une cellule encore vide apres Newton (amortissement
            // epuise) est retiree et les poids sont resolus a nouveau
            std::vector<int> keep;
            for (std::size_t i = 0; i < cells.size(); ++i)
                if (cells[i].size() >= 3) keep.push_back((int)i);
            if (keep.size() < cells.size()) {
                std::vector<Vec2> s2;
                std::vector<double> p2;
                for (int i : keep) { s2.push_back(seeds[i]); p2.push_back(sp[i]); }
                nRedundant += (int)(cells.size() - keep.size());
                seeds.swap(s2);
                sp.swap(p2);
                setTargets();
                wgt.assign(seeds.size(), 0.0);
                --it;                  // recalcul sans consommer une iteration
                continue;
            }
        }
        if (it >= lloyd) break;
'''
c = c[:i0] + new + c[i1:]
c=rep(c,"""    if (nRedundant > 0)
        std::printf("[tess] laguerre: %d graine(s) redondante(s) retiree(s) "
""","""    if (sizeSpread > 0.0)
        std::printf("[tess] laguerre: aires prescrites (sigma_ln demande %.3g) "
                    "atteintes a %.3f %% pres en %d iteration(s) de Newton "
                    "(dernier cycle)\\n", sizeSpread, 100.0 * newtonErr, newtonIt);
    if (nRedundant > 0)
        std::printf("[tess] laguerre: %d graine(s) redondante(s) retiree(s) "
""","printf")
# ---- hpp : commentaire ----
h=rep(h,"""    //            cellules sont celles du diagramme de LAGUERRE de poids
    //            (0,315 s_i)^2 : la graine lourde recoit la grande cellule (un
    //            Voronoi ordinaire moyennerait les tailles des voisines et
    //            ecraserait la dispersion : 0,5 demande -> 0,16 realise).
""","""    //            cellules sont celles du diagramme de LAGUERRE dont les poids
    //            sont resolus (Newton amorti, Kitagawa-Merigot-Thibert 2019)
    //            pour que l aire de chaque cellule soit EXACTEMENT ~ s_i^2
    //            (Bourne et al. 2020). Un Voronoi ordinaire moyenne les
    //            tailles des voisines (0,5 demande -> 0,16 realise) ; des
    //            poids fixes (kappa s_i)^2 n en rendent que 40 %.
""","hpp")
io.open(H,"w",encoding="utf-8").write(h); io.open(C,"w",encoding="utf-8").write(c)
print("newton : 5 greffes appliquees (voronoiCells lap x2, build x2, hpp)")
