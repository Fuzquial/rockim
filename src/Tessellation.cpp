// ---------------------------------------------------------------------------
// Tessellation — Voronoi grains + phases. See the header for the pipeline.
// ---------------------------------------------------------------------------
#include "rockim/Tessellation.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace rockim {

namespace {

using Vec2 = Eigen::Vector2d;

// Keep the part of a CCW polygon with n . p <= b (Sutherland-Hodgman step).
// Clipping a convex CCW polygon by a half-plane preserves convexity and CCW.
std::vector<Vec2> clipHalfPlane(const std::vector<Vec2>& poly,
                                const Vec2& n, double b) {
    std::vector<Vec2> out;
    out.reserve(poly.size() + 2);
    const int k = (int)poly.size();
    for (int i = 0; i < k; ++i) {
        const Vec2& P = poly[i];
        const Vec2& Q = poly[(i + 1) % k];
        double dP = n.dot(P) - b, dQ = n.dot(Q) - b;
        bool inP = dP <= 0.0, inQ = dQ <= 0.0;
        if (inP) out.push_back(P);
        if (inP != inQ) out.push_back(P + (dP / (dP - dQ)) * (Q - P));
    }
    return out;
}

double polyArea(const std::vector<Vec2>& p) {
    double a = 0.0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        const Vec2& A = p[i];
        const Vec2& B = p[(i + 1) % p.size()];
        a += A.x() * B.y() - B.x() * A.y();
    }
    return 0.5 * a;
}

Vec2 polyCentroid(const std::vector<Vec2>& p) {
    double a = 0.0;
    Vec2 c = Vec2::Zero();
    for (std::size_t i = 0; i < p.size(); ++i) {
        const Vec2& A = p[i];
        const Vec2& B = p[(i + 1) % p.size()];
        double w = A.x() * B.y() - B.x() * A.y();
        a += w;
        c += w * (A + B);
    }
    if (std::abs(a) < 1e-30) return p.empty() ? Vec2::Zero() : p[0];
    return c / (3.0 * a);
}

// Voronoi cells of all seeds inside the rectangle. Neighbours are visited in
// order of increasing distance with the exact early-out: once
// 0.5 * |sj - si| exceeds the current cell circumradius (max vertex distance
// from the seed), no farther bisector can cut the cell.
// wgt (2026-09-02, optionnel) : diagramme de LAGUERRE (cellules de puissance).
// Avec des poids w_i, la frontiere entre i et j n est plus la mediatrice mais
// la ligne radicale { x : |x-p_i|^2 - w_i = |x-p_j|^2 - w_j }, a distance
// t_i = (d^2 + w_i - w_j) / (2 d) de p_i : la graine la plus lourde recoit la
// plus grande cellule. C est la construction standard des microstructures
// polydisperses (Neper, Quey & Renversade 2018 ; Falco et al. 2017). Sans
// poids (nullptr) le chemin historique est execute a l identique.
//   - coupe possible seulement si t_i < R, soit d < R + sqrt(R^2 + w_j - w_i)
//     <= 2 R + sqrt(max(0, w_max - w_i)) : borne de l early-out et de la
//     couverture en anneaux ;
//   - une graine dont la cellule de puissance est VIDE est redondante (cachee
//     par des voisines plus lourdes) : la cellule est rendue vide et l appelant
//     la retire.
// lap (optionnel, mode Laguerre) : pour chaque cellule i, la liste des
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
    double wMax = 0.0;
    if (lag) for (double v : *wgt) wMax = std::max(wMax, v);
    std::vector<std::vector<Vec2>> cells(N);
    const std::vector<Vec2> rect = {
        {0.0, 0.0}, {W, 0.0}, {W, H}, {0.0, H}};

    // distance-sorted neighbour order per seed, via a coarse uniform grid
    double cs = std::sqrt(W * H / std::max(1, N));      // ~1 seed per cell
    int gx = std::max(1, (int)(W / cs)), gy = std::max(1, (int)(H / cs));
    std::vector<std::vector<int>> grid((std::size_t)gx * gy);
    auto cellOf = [&](const Vec2& p, int& cx, int& cy) {
        cx = std::clamp((int)(p.x() / W * gx), 0, gx - 1);
        cy = std::clamp((int)(p.y() / H * gy), 0, gy - 1);
    };
    for (int i = 0; i < N; ++i) {
        int cx, cy;
        cellOf(seeds[i], cx, cy);
        grid[(std::size_t)cy * gx + cx].push_back(i);
    }

    // Coverage guarantee of a ring block: any seed OUTSIDE the block
    // [ci +- ring, cj +- ring] differs by at least ring cells on one axis,
    // i.e. lies at Euclidean distance >= ring * min(cell width, height).
    // A bisector can only cut the cell if its seed is closer than twice the
    // cell circumradius R, so the block is PROVEN sufficient only once
    // 2 R <= ring * csMin — checking the early-out on the candidates alone
    // is NOT enough (a cutting seed can be absent from the block: this
    // produced overlapping cells and a +0.23 % area excess on the
    // percussion domain before the fix).
    double csx = W / gx, csy = H / gy;
    double csMin = std::min(csx, csy);
    std::vector<std::pair<double, int>> cand;
    for (int i = 0; i < N; ++i) {
        int ci, cj;
        cellOf(seeds[i], ci, cj);
        std::vector<Vec2> poly;
        const double dwi = lag ? std::sqrt(std::max(0.0, wMax - (*wgt)[i])) : 0.0;

        for (int ring = 1;; ++ring) {
            cand.clear();
            for (int dj = -ring; dj <= ring; ++dj)
                for (int di = -ring; di <= ring; ++di) {
                    int cx = ci + di, cy = cj + dj;
                    if (cx < 0 || cy < 0 || cx >= gx || cy >= gy) continue;
                    for (int j : grid[(std::size_t)cy * gx + cx])
                        if (j != i)
                            cand.push_back({(seeds[j] - seeds[i]).squaredNorm(), j});
                }
            std::sort(cand.begin(), cand.end());

            poly = rect;                     // re-clip from scratch (cheap)
            double r2max = 0.0;
            for (const Vec2& v : poly)
                r2max = std::max(r2max, (v - seeds[i]).squaredNorm());
            for (std::size_t k = 0; k < cand.size(); ++k) {
                if (!lag) {
                    if (0.25 * cand[k].first > r2max) break;   // sorted: none closer
                } else if (std::sqrt(cand[k].first) > 2.0 * std::sqrt(r2max) + dwi) {
                    break;                                     // borne de Laguerre
                }
                int j = cand[k].second;
                Vec2 n = seeds[j] - seeds[i];
                double b = lag ? n.dot(seeds[i]) + 0.5 * (cand[k].first + (*wgt)[i] - (*wgt)[j])
                               : n.dot(0.5 * (seeds[i] + seeds[j]));
                poly = clipHalfPlane(poly, n, b);
                if (poly.size() < 3) {
                    if (lag) { poly.clear(); break; }          // graine redondante
                    throw std::runtime_error("Tessellation: seed cell vanished "
                                             "(coincident seeds?)");
                }
                r2max = 0.0;
                for (const Vec2& v : poly)
                    r2max = std::max(r2max, (v - seeds[i]).squaredNorm());
            }
            bool coveredAll = (ci - ring < 0 && cj - ring < 0
                               && ci + ring >= gx && cj + ring >= gy);
            if (lag && poly.empty()) break;
            double rcov = ring * csMin;
            if (lag) {
                if (coveredAll || 2.0 * std::sqrt(r2max) + dwi <= rcov) break;
            } else if (coveredAll || 4.0 * r2max <= rcov * rcov) break;
        }
        if (lag && lap && poly.size() >= 3) {
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

// ---------------------------------------------------------------------------
// Bowyer-Watson Delaunay triangulation of a planar point set, returning CCW
// index triples. Used to mesh the INSIDE of a Voronoi grain the way the
// grain-based FDEM literature does. The cells are convex (half-space clipping),
// and the Delaunay triangulation of a point set covers exactly its convex hull,
// so no constrained triangulation is needed: the boundary points are on the
// hull by construction. The caller checks the area anyway, because the
// short-edge contraction can leave a cell marginally non-convex.
//
// Deterministic: insertion order is the caller's point order, no randomness.
//
// Namespace-scope (declared in Tessellation.hpp) since the native disc mesher
// of FdemSolver uses it too; the rest of this file stays internal.
// ---------------------------------------------------------------------------
} // namespace (internal helpers above; delaunayCCW is exported)
std::vector<std::array<int, 3>> delaunayCCW(const std::vector<Vec2>& p) {
    const int n = (int)p.size();
    if (n < 3) return {};
    Vec2 lo = p[0], hi = p[0];
    for (const Vec2& q : p) {
        lo = lo.cwiseMin(q);
        hi = hi.cwiseMax(q);
    }
    Vec2 c = 0.5 * (lo + hi);
    double R = std::max((hi - lo).norm(), 1e-12);
    std::vector<Vec2> pt(p);
    pt.push_back(Vec2(c.x() - 20.0 * R, c.y() - R));       // super-triangle
    pt.push_back(Vec2(c.x() + 20.0 * R, c.y() - R));
    pt.push_back(Vec2(c.x(), c.y() + 20.0 * R));

    auto ccw = [&](int a, int b, int d) {
        return (pt[b].x() - pt[a].x()) * (pt[d].y() - pt[a].y())
             - (pt[d].x() - pt[a].x()) * (pt[b].y() - pt[a].y());
    };
    // strictly-inside test of the circumcircle of the CCW triangle (a,b,d)
    auto inCircle = [&](int a, int b, int d, int e) {
        double ax = pt[a].x() - pt[e].x(), ay = pt[a].y() - pt[e].y();
        double bx = pt[b].x() - pt[e].x(), by = pt[b].y() - pt[e].y();
        double cx = pt[d].x() - pt[e].x(), cy = pt[d].y() - pt[e].y();
        double det = (ax * ax + ay * ay) * (bx * cy - cx * by)
                   - (bx * bx + by * by) * (ax * cy - cx * ay)
                   + (cx * cx + cy * cy) * (ax * by - bx * ay);
        return det > 1e-18 * R * R * R * R;
    };

    std::vector<std::array<int, 3>> tri{{n, n + 1, n + 2}};
    std::vector<std::array<int, 2>> cav;
    std::vector<std::array<int, 3>> keep;
    for (int i = 0; i < n; ++i) {
        cav.clear();
        keep.clear();
        for (const auto& t : tri) {
            if (inCircle(t[0], t[1], t[2], i)) {
                cav.push_back({t[0], t[1]});
                cav.push_back({t[1], t[2]});
                cav.push_back({t[2], t[0]});
            } else {
                keep.push_back(t);
            }
        }
        // cavity boundary = the edges appearing exactly once
        tri.swap(keep);
        for (std::size_t a = 0; a < cav.size(); ++a) {
            bool shared = false;
            for (std::size_t b = 0; b < cav.size() && !shared; ++b)
                if (a != b && cav[a][0] == cav[b][1] && cav[a][1] == cav[b][0])
                    shared = true;
            if (shared) continue;
            int u = cav[a][0], v = cav[a][1];
            if (ccw(u, v, i) > 0.0) tri.push_back({u, v, i});
            else                    tri.push_back({v, u, i});
        }
    }
    std::vector<std::array<int, 3>> out;
    out.reserve(tri.size());
    for (const auto& t : tri)
        if (t[0] < n && t[1] < n && t[2] < n && ccw(t[0], t[1], t[2]) > 0.0)
            out.push_back(t);
    return out;
}

Tessellation Tessellation::build(double W, double H, double targetD,
                                 double jitter, int lloyd, double mergeFrac,
                                 int refine,
                                 const std::vector<double>& phaseFraction,
                                 std::mt19937& rng, bool randomSeeds,
                                 bool useDelaunay, double elemSize,
                                 double sizeSpread,
                                 const std::vector<double>& phaseSize,
                                 bool randomInterior) {
    if (targetD <= 0.0 || targetD > 0.5 * std::min(W, H))
        throw std::runtime_error("Tessellation: grainSize must be positive and "
                                 "much smaller than the domain");
    if (refine < 0 || refine > 4)
        throw std::runtime_error("Tessellation: refineLevels must be in 0..4 "
                                 "(each level multiplies the element count "
                                 "by 4)");

    // ---- 1. seeds ------------------------------------------------------------
    // hex cell area = sqrt(3)/2 s^2 matched to the disc area pi/4 d^2
    double s = targetD * std::sqrt(M_PI / (2.0 * std::sqrt(3.0)));
    double sy = s * std::sqrt(3.0) / 2.0;
    std::vector<Vec2> seeds;
    std::vector<double> sp;      // espacement propre a chaque graine (polydispersite)
    double eps = 1e-6 * std::min(W, H);
    if (randomSeeds && sizeSpread > 0.0) {
        // POLYDISPERSITE (2026-09-02) : addition sequentielle TRIEE. Les N
        // espacements s_i = s L_i (L_i log-normale de moyenne 1, ecart-type
        // de ln = sizeSpread) sont tires d abord, puis places du plus GRAND
        // au plus petit, chaque graine acceptee a distance
        // >= 0,35 (s_i + s_j) de toutes les precedentes. Tirer la taille a
        // chaque essai (Poisson-disc ordinaire) rejette surtout les grandes
        // graines, qui exigent le plus de place : l ecart-type realise
        // s effondrait (0,5 demande -> 0,20 obtenu). N est divise par
        // E[L^2] = exp(sigma^2) pour conserver la taille moyenne.
        long target = std::lround(W * H / (std::sqrt(3.0) / 2.0 * s * s
                                           * std::exp(sizeSpread * sizeSpread)));
        target = std::max(2L, target);
        std::normal_distribution<double> N01(0.0, 1.0);
        std::vector<double> want((std::size_t)target);
        for (double& v : want)
            v = s * std::exp(sizeSpread * N01(rng) - 0.5 * sizeSpread * sizeSpread);
        std::sort(want.begin(), want.end(), std::greater<double>());
        const double cs = 0.7 * want.front();     // >= toute distance d acceptation
        int gx = std::max(1, (int)(W / cs)), gy = std::max(1, (int)(H / cs));
        std::vector<std::vector<int>> acc((std::size_t)gx * gy);
        std::uniform_real_distribution<double> Ux(eps, W - eps);
        std::uniform_real_distribution<double> Uy(eps, H - eps);
        int skipped = 0;
        for (double spi : want) {
            bool placed = false;
            for (int tr = 0; tr < 4000 && !placed; ++tr) {
                Vec2 p(Ux(rng), Uy(rng));
                int ci = std::clamp((int)(p.x() / W * gx), 0, gx - 1);
                int cj = std::clamp((int)(p.y() / H * gy), 0, gy - 1);
                bool ok = true;
                for (int dj = -1; dj <= 1 && ok; ++dj)
                    for (int di = -1; di <= 1 && ok; ++di) {
                        int cx = ci + di, cy = cj + dj;
                        if (cx < 0 || cy < 0 || cx >= gx || cy >= gy) continue;
                        for (int q : acc[(std::size_t)cy * gx + cx])
                            if ((seeds[q] - p).norm() < 0.35 * (spi + sp[q])) { ok = false; break; }
                    }
                if (!ok) continue;
                acc[(std::size_t)cj * gx + ci].push_back((int)seeds.size());
                seeds.push_back(p);
                sp.push_back(spi);
                placed = true;
            }
            if (!placed) ++skipped;
        }
        if (skipped > 0)
            std::printf("[tess] polydispersite: %d graine(s) sur %ld non placee(s) "
                        "(domaine sature) - l ecart-type realise en tient compte\n",
                        skipped, target);
    } else if (randomSeeds) {
        // Poisson-disc dart throwing: uniform candidates, accepted when at
        // least dmin from every accepted seed. Unlike the jittered hex
        // lattice (whose Voronoi keeps three preferred boundary
        // orientations even after Lloyd), the resulting boundary
        // orientation distribution is isotropic.
        double dmin = 0.7 * s;
        long target = std::lround(W * H / (std::sqrt(3.0) / 2.0 * s * s));
        double cs = dmin;
        int gx = std::max(1, (int)(W / cs)), gy = std::max(1, (int)(H / cs));
        std::vector<std::vector<int>> acc((std::size_t)gx * gy);
        std::uniform_real_distribution<double> Ux(eps, W - eps);
        std::uniform_real_distribution<double> Uy(eps, H - eps);
        long attempts = 0, maxAttempts = 400 * std::max(1L, target);
        while ((long)seeds.size() < target && attempts < maxAttempts) {
            ++attempts;
            Vec2 p(Ux(rng), Uy(rng));
            int ci = std::clamp((int)(p.x() / W * gx), 0, gx - 1);
            int cj = std::clamp((int)(p.y() / H * gy), 0, gy - 1);
            bool ok = true;
            for (int dj = -1; dj <= 1 && ok; ++dj)
                for (int di = -1; di <= 1 && ok; ++di) {
                    int cx = ci + di, cy = cj + dj;
                    if (cx < 0 || cy < 0 || cx >= gx || cy >= gy) continue;
                    for (int q : acc[(std::size_t)cy * gx + cx])
                        if ((seeds[q] - p).norm() < dmin) { ok = false; break; }
                }
            if (!ok) continue;
            acc[(std::size_t)cj * gx + ci].push_back((int)seeds.size());
            seeds.push_back(p);
        }
    } else {
        if (sizeSpread > 0.0)
            throw std::runtime_error("Tessellation: grainSizeSpread exige "
                                     "grainSeeding = random (le reseau "
                                     "hexagonal n a qu une taille)");
        std::uniform_real_distribution<double> U(-1.0, 1.0);
        for (int j = 0;; ++j) {
            double y = (j + 0.5) * sy;
            if (y >= H) break;
            double off = (j % 2 == 0) ? 0.25 * s : -0.25 * s;
            for (int i = 0;; ++i) {
                double x = (i + 0.5) * s + off;
                if (x >= W) break;
                Vec2 p(x + jitter * 0.5 * s * U(rng),
                       y + jitter * 0.5 * sy * U(rng));
                p.x() = std::clamp(p.x(), eps, W - eps);
                p.y() = std::clamp(p.y(), eps, H - eps);
                seeds.push_back(p);
            }
        }
    }
    if (seeds.size() < 2)
        throw std::runtime_error("Tessellation: fewer than 2 grains — "
                                 "reduce grainSize");

    // ---- 2. Voronoi cells (+ Lloyd relaxation) -------------------------------
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
        const bool dbg = std::getenv("ROCKIM_TESS_DEBUG") != nullptr;
        std::vector<std::vector<Vec2>> cl = voronoiCells(seeds, W, H, &wgt, &lap);
        std::vector<double> A(N);
        bool anyEmpty = false;
        for (std::size_t i = 0; i < N; ++i) {
            A[i] = cl[i].size() >= 3 ? polyArea(cl[i]) : 0.0;
            if (A[i] <= 0.0) anyEmpty = true;
        }
        if (anyEmpty) {            // depart de Newton amorti : Voronoi, aucune cellule vide
            std::fill(wgt.begin(), wgt.end(), 0.0);
            cl = voronoiCells(seeds, W, H, &wgt, &lap);
            for (std::size_t i = 0; i < N; ++i) A[i] = cl[i].size() >= 3 ? polyArea(cl[i]) : 0.0;
        }
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
            int cgIt = 0;
            for (int cg = 0; cg < 1000 && rr > 1e-24 * gg; ++cg, ++cgIt) {
                matvec(p, q);
                const double pq = dot(p, q);
                if (!(pq > 0.0)) break;                    // matrice degeneree : on garde delta courant
                const double alpha = rr / pq;
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
            if (dbg) {
                double dmax = 0.0, amin = 1e300;
                std::size_t nEdges = 0;
                for (std::size_t i = 0; i < N; ++i) { dmax = std::max(dmax, std::abs(delta[i])); amin = std::min(amin, A2[i]); nEdges += lap[i].size(); }
                std::printf("[tess-dbg] newton %d : err %.3f  cg %d  aretes %zu  |delta|max %.3e  tau %.4f  amin/eps0 %.3f\n",
                            k, err, cgIt, nEdges, dmax, tau, amin / eps0);
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
        for (std::size_t i = 0; i < seeds.size(); ++i) {
            Vec2 c = polyCentroid(cells[i]);
            seeds[i].x() = std::clamp(c.x(), eps, W - eps);
            seeds[i].y() = std::clamp(c.y(), eps, H - eps);
        }
    }

    if (sizeSpread > 0.0)
        std::printf("[tess] laguerre: aires prescrites (sigma_ln demande %.3g) "
                    "atteintes a %.3f %% pres en %d iteration(s) de Newton "
                    "(dernier cycle)\n", sizeSpread, 100.0 * newtonErr, newtonIt);
    if (nRedundant > 0)
        std::printf("[tess] laguerre: %d graine(s) redondante(s) retiree(s) "
                    "(cellule de puissance vide)\n", nRedundant);

    Tessellation T;
    T.nGrains = (int)cells.size();

    // ---- 3. short-EDGE contraction (the sliver guard) ------------------------
    // Naive spatial clustering of vertices is WRONG here: it can merge two
    // vertices that are close in space but not connected by a diagram edge
    // (e.g. an interior vertex and a wall vertex separated by a thin
    // boundary cell). Their incident cells then sweep across the cell in
    // between and the tiling gains double-covered area (measured: +0.23 %
    // on the percussion domain). The safe operation is EDGE CONTRACTION:
    // only vertices connected by a diagram edge shorter than tol are merged
    // (transitively, via union-find), every incident cell is updated
    // consistently, and the merged position is the component mean snapped
    // back onto any domain wall a member lay on. The strict area check
    // below remains as the backstop.
    double tol = std::max(mergeFrac * targetD, 1e-12);
    double epsF = 1e-9 * std::min(W, H);               // fp-copy identification

    // 3a. exact dedup of the per-cell vertex copies into diagram vertices
    std::vector<Vec2> vFine;
    std::unordered_map<long long, std::vector<int>> hash;
    auto hkey = [&](int ix, int iy) {
        return ((long long)ix << 32) ^ (unsigned long long)(unsigned int)iy;
    };
    auto fineId = [&](const Vec2& p) {
        int ix = (int)std::floor(p.x() / (4.0 * epsF));
        int iy = (int)std::floor(p.y() / (4.0 * epsF));
        for (int dj = -1; dj <= 1; ++dj)
            for (int di = -1; di <= 1; ++di) {
                auto it = hash.find(hkey(ix + di, iy + dj));
                if (it == hash.end()) continue;
                for (int id : it->second)
                    if ((vFine[id] - p).norm() < epsF) return id;
            }
        int id = (int)vFine.size();
        vFine.push_back(p);
        hash[hkey(ix, iy)].push_back(id);
        return id;
    };
    std::vector<std::vector<int>> cellFine(cells.size());
    for (std::size_t g = 0; g < cells.size(); ++g) {
        std::vector<int>& ids = cellFine[g];
        for (const Vec2& v : cells[g]) {
            int id = fineId(v);
            if (ids.empty() || ids.back() != id) ids.push_back(id);
        }
        while (ids.size() > 1 && ids.back() == ids.front()) ids.pop_back();
        if (ids.size() < 3)
            throw std::runtime_error("Tessellation: degenerate cell polygon");
    }

    // 3b. union-find over short diagram edges
    std::vector<int> parent(vFine.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::function<int(int)> find = [&](int a) {
        while (parent[a] != a) { parent[a] = parent[parent[a]]; a = parent[a]; }
        return a;
    };
    for (const auto& ids : cellFine)
        for (std::size_t k = 0; k < ids.size(); ++k) {
            int a = ids[k], b = ids[(k + 1) % ids.size()];
            if ((vFine[a] - vFine[b]).norm() < tol) {
                int ra = find(a), rb = find(b);
                if (ra != rb) parent[std::max(ra, rb)] = std::min(ra, rb);
            }
        }

    // 3c. component positions: mean of members, snapped onto any wall a
    // member lies on (wall membership detected at fp precision — clipping
    // puts boundary vertices exactly on the walls)
    std::unordered_map<int, int> rootToNew;
    std::vector<Vec2> sum;
    std::vector<int> cnt, flags;
    auto wallOf = [&](const Vec2& p) {
        return (p.x() < epsF ? 1 : 0) | (p.x() > W - epsF ? 2 : 0)
             | (p.y() < epsF ? 4 : 0) | (p.y() > H - epsF ? 8 : 0);
    };
    for (int v = 0; v < (int)vFine.size(); ++v) {
        int r = find(v);
        auto [it, isNew] = rootToNew.try_emplace(r, (int)sum.size());
        if (isNew) { sum.push_back(Vec2::Zero()); cnt.push_back(0); flags.push_back(0); }
        int id = it->second;
        sum[id] += vFine[v];
        cnt[id] += 1;
        flags[id] |= wallOf(vFine[v]);
    }
    T.vtx.resize(sum.size());
    for (std::size_t id = 0; id < sum.size(); ++id) {
        Vec2 p = sum[id] / cnt[id];
        if (flags[id] & 1) p.x() = 0.0;
        if (flags[id] & 2) p.x() = W;
        if (flags[id] & 4) p.y() = 0.0;
        if (flags[id] & 8) p.y() = H;
        T.vtx[id] = p;
    }

    // 3d. rebuild the cell polygons on the contracted vertices
    std::vector<std::vector<int>> cellIds(cells.size());
    for (std::size_t g = 0; g < cells.size(); ++g) {
        std::vector<int>& ids = cellIds[g];
        for (int f : cellFine[g]) {
            int id = rootToNew.at(find(f));
            if (ids.empty() || ids.back() != id) ids.push_back(id);
        }
        while (ids.size() > 1 && ids.back() == ids.front()) ids.pop_back();
        if (ids.size() < 3)
            throw std::runtime_error("Tessellation: a grain degenerated during "
                                     "edge contraction — reduce vertexMergeFrac");
        // a repeated non-consecutive vertex would pinch the polygon (and
        // silently corrupt areas): reject loudly, this indicates a needle-
        // shaped cell that contraction cannot fix at this tolerance
        std::vector<int> sorted = ids;
        std::sort(sorted.begin(), sorted.end());
        if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
            throw std::runtime_error("Tessellation: pinched cell after edge "
                                     "contraction — reduce vertexMergeFrac");
    }

    // ---- 4. phase assignment (area-greedy on shuffled grains) ----------------
    std::vector<double> frac = phaseFraction;
    if (frac.empty()) frac = {1.0};
    double fsum = std::accumulate(frac.begin(), frac.end(), 0.0);
    if (fsum <= 0.0) throw std::runtime_error("Tessellation: phase fractions "
                                              "must sum to a positive value");
    for (double& f : frac) f /= fsum;

    T.grainArea.resize(cells.size());
    double Atot = 0.0;
    for (std::size_t g = 0; g < cells.size(); ++g) {
        std::vector<Vec2> poly;
        for (int id : cellIds[g]) poly.push_back(T.vtx[id]);
        T.grainArea[g] = polyArea(poly);
        Atot += T.grainArea[g];
    }
    if (std::abs(Atot - W * H) > 1e-6 * W * H) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
                      "Tessellation: cells do not tile the domain after "
                      "welding (sum %.9g vs %.9g, mismatch %.3g = %.2e rel)",
                      Atot, W * H, Atot - W * H, (Atot - W * H) / (W * H));
        throw std::runtime_error(msg);
    }

    std::vector<int> order(cells.size());
    std::iota(order.begin(), order.end(), 0);
    std::vector<double> assigned(frac.size(), 0.0);
    T.phaseOfGrain.assign(cells.size(), 0);
    bool anyPhaseSize = false;
    for (double v : phaseSize) if (v > 0.0) anyPhaseSize = true;
    if (!anyPhaseSize) {
        // chemin historique : grains melanges, glouton sur le deficit d aire
        std::shuffle(order.begin(), order.end(), rng);
        for (int g : order) {
            int best = 0;
            double bestDef = -1e300;
            for (std::size_t p = 0; p < frac.size(); ++p) {
                double def = frac[p] * Atot - assigned[p];
                if (def > bestDef) { bestDef = def; best = (int)p; }
            }
            T.phaseOfGrain[g] = best;
            assigned[best] += T.grainArea[g];
        }
    } else {
        // TAILLE PAR PHASE (2026-09-02) : grains pris du plus grand au plus
        // petit ; chaque phase a une affinite log-normale (sigma 0,35) a sa
        // taille cible, multipliee par son deficit RELATIF d aire. Le deficit
        // est en aire, donc les fractions globales sont respectees quelle que
        // soit la distribution des tailles ; l affinite ne fait que choisir
        // QUELS grains vont a quelle phase. Une phase saturee (deficit <= 0)
        // ne recrute plus tant qu une autre a du deficit.
        if (phaseSize.size() != frac.size())
            throw std::runtime_error("Tessellation: phaseSize doit avoir une "
                                     "entree par phase");
        std::sort(order.begin(), order.end(), [&](int a, int b) {
            if (T.grainArea[a] != T.grainArea[b]) return T.grainArea[a] > T.grainArea[b];
            return a < b;
        });
        for (int g : order) {
            const double dg = 2.0 * std::sqrt(T.grainArea[g] / M_PI);
            int best = -1;
            double bestScore = -1e300;
            for (std::size_t p = 0; p < frac.size(); ++p) {
                const double def = frac[p] * Atot - assigned[p];
                if (def <= 0.0) continue;
                double aff = 1.0;
                if (phaseSize[p] > 0.0) {
                    const double z = std::log(dg / phaseSize[p]) / 0.35;
                    aff = std::exp(-0.5 * z * z);
                }
                const double score = aff * def / (frac[p] * Atot);
                if (score > bestScore) { bestScore = score; best = (int)p; }
            }
            if (best < 0) {                        // toutes saturees : max deficit
                double bestDef = -1e300;
                for (std::size_t p = 0; p < frac.size(); ++p) {
                    const double def = frac[p] * Atot - assigned[p];
                    if (def > bestDef) { bestDef = def; best = (int)p; }
                }
            }
            T.phaseOfGrain[g] = best;
            assigned[best] += T.grainArea[g];
        }
    }

    // ---- 5. intra-grain meshing ---------------------------------------------
    // Two front-ends. `fan` is the original: one triangle per cell edge, all
    // sharing the centroid. It is cheap and star-shaped, but it fixes the
    // intra-grain topology to a HUB — a transgranular crack can only run along
    // a spoke. `delaunay` is what the grain-based FDEM literature actually does
    // (Y-Geo / Irazu GBM studies mesh the interior of every grain with an
    // unstructured Delaunay mesh at roughly 0.18 x the grain diameter, i.e.
    // about six elements across a grain), precisely so that intergranular AND
    // transgranular cracking are both representable. Boundary edges are
    // subdivided from the edge's own geometry, keyed on its sorted vertex ids,
    // so the two cells sharing it produce IDENTICAL points and the cohesive
    // joints stay paired.
    if (!useDelaunay) {
        for (std::size_t g = 0; g < cells.size(); ++g) {
            const std::vector<int>& ids = cellIds[g];
            std::vector<Vec2> poly;
            for (int id : ids) poly.push_back(T.vtx[id]);
            Vec2 cen = polyCentroid(poly);
            int cenId = (int)T.vtx.size();
            T.vtx.push_back(cen);
            for (std::size_t k = 0; k < ids.size(); ++k) {
                int a = ids[k], b = ids[(k + 1) % ids.size()];
                const Vec2& A = T.vtx[a];
                const Vec2& B = T.vtx[b];
                double ar = 0.5 * ((A.x() - cen.x()) * (B.y() - cen.y())
                                   - (B.x() - cen.x()) * (A.y() - cen.y()));
                if (ar <= 0.0)
                    throw std::runtime_error("Tessellation: non-CCW fan triangle "
                                             "(cell not star-shaped after welding) "
                                             "— reduce vertexMergeFrac");
                T.tri.push_back({{cenId, a, b}, (int)g});
            }
        }
    } else {
        double h = elemSize > 0.0 ? elemSize : 0.18 * targetD;
        // shared subdivision of every Voronoi edge
        std::map<std::pair<int, int>, std::vector<int>> split;
        auto edgePoints = [&](int a, int b) -> const std::vector<int>& {
            auto key = std::minmax(a, b);
            auto it = split.find({key.first, key.second});
            if (it != split.end()) return it->second;
            // COPIES, not references: the loop below push_back()s into T.vtx,
            // which reallocates and would leave A and B dangling. That bug
            // segfaulted the mesher as soon as a domain was large enough for
            // the reallocation to land inside the loop.
            const Vec2 A = T.vtx[key.first];
            const Vec2 B = T.vtx[key.second];
            double L = (B - A).norm();
            int n = std::max(1, (int)std::llround(L / h));
            std::vector<int> ids;                      // interior points only
            for (int k = 1; k < n; ++k) {
                ids.push_back((int)T.vtx.size());
                T.vtx.push_back(A + (B - A) * ((double)k / n));
            }
            return split.emplace(std::make_pair(key.first, key.second),
                                 std::move(ids)).first->second;
        };

        long nFallback = 0;
        for (std::size_t g = 0; g < cells.size(); ++g) {
            const std::vector<int>& ids = cellIds[g];
            std::vector<Vec2> poly;
            for (int id : ids) poly.push_back(T.vtx[id]);

            // CCW boundary loop with the shared subdivisions inserted
            std::vector<int> loop;
            for (std::size_t k = 0; k < ids.size(); ++k) {
                int a = ids[k], b = ids[(k + 1) % ids.size()];
                loop.push_back(a);
                const std::vector<int>& mid = edgePoints(a, b);
                if (a < b)
                    for (int m : mid) loop.push_back(m);
                else
                    for (auto it = mid.rbegin(); it != mid.rend(); ++it)
                        loop.push_back(*it);
            }

            // interior points on a triangular lattice, kept one margin away
            // from every edge so the Delaunay does not produce rim slivers
            std::vector<int> inner;
            Vec2 cen = polyCentroid(poly);
            // test commun : dans le polygone ET a plus de 0,55 h de chaque arete
            auto insideMargin = [&](const Vec2& p) {
                for (std::size_t k = 0; k < poly.size(); ++k) {
                    const Vec2& A = poly[k];
                    const Vec2& B = poly[(k + 1) % poly.size()];
                    Vec2 e = B - A;
                    double Le = e.norm();
                    if (Le < 1e-15) continue;
                    double d = (e.x() * (p.y() - A.y())
                                - e.y() * (p.x() - A.x())) / Le;
                    if (d < 0.55 * h) return false;     // outside or too close
                }
                return true;
            };
            if (!randomInterior) {
            double rmax = 0.0;
            for (const Vec2& p : poly) rmax = std::max(rmax, (p - cen).norm());
            int nr = (int)std::ceil(rmax / h) + 1;
            double dy = h * std::sqrt(3.0) / 2.0;
            for (int j = -nr; j <= nr; ++j)
                for (int i = -nr; i <= nr; ++i) {
                    Vec2 p(cen.x() + (i + 0.5 * (j & 1)) * h, cen.y() + j * dy);
                    bool ok = true;
                    for (std::size_t k = 0; k < poly.size() && ok; ++k) {
                        const Vec2& A = poly[k];
                        const Vec2& B = poly[(k + 1) % poly.size()];
                        Vec2 e = B - A;
                        double Le = e.norm();
                        if (Le < 1e-15) continue;
                        double d = (e.x() * (p.y() - A.y())
                                    - e.y() * (p.x() - A.x())) / Le;
                        if (d < 0.55 * h) ok = false;  // outside or too close
                    }
                    if (ok) { inner.push_back((int)T.vtx.size()); T.vtx.push_back(p); }
                }
            } else {
                // POINTS INTERIEURS ALEATOIRES (2026-09-02, opt-in
                // grainMeshRandom) : Poisson-disc dans le polygone, distance
                // >= 0,75 h entre points (densite equivalente au reseau a la
                // saturation du lancer), meme marge 0,55 h aux aretes. Le
                // reseau triangulaire pave les grains en triangles quasi
                // equilateraux ALIGNES — le defaut du frontal de Gmsh, banni.
                double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
                for (const Vec2& p : poly) {
                    xmin = std::min(xmin, p.x()); xmax = std::max(xmax, p.x());
                    ymin = std::min(ymin, p.y()); ymax = std::max(ymax, p.y());
                }
                const double aP = std::abs(polyArea(poly));
                const long target = std::lround(aP / (std::sqrt(3.0) / 2.0 * h * h));
                const double dmin = 0.75 * h;
                std::uniform_real_distribution<double> Ux(xmin, xmax), Uy(ymin, ymax);
                std::vector<Vec2> acc;
                long attempts = 0;
                const long maxAttempts = 80 * std::max(1L, target);
                while ((long)acc.size() < target && attempts < maxAttempts) {
                    ++attempts;
                    Vec2 p(Ux(rng), Uy(rng));
                    if (!insideMargin(p)) continue;
                    bool ok = true;
                    for (const Vec2& q : acc)
                        if ((q - p).norm() < dmin) { ok = false; break; }
                    if (!ok) continue;
                    acc.push_back(p);
                }
                for (const Vec2& p : acc) {
                    inner.push_back((int)T.vtx.size());
                    T.vtx.push_back(p);
                }
            }

            std::vector<int> all(loop);
            all.insert(all.end(), inner.begin(), inner.end());
            std::vector<Vec2> pts;
            pts.reserve(all.size());
            for (int id : all) pts.push_back(T.vtx[id]);

            std::vector<std::array<int, 3>> tri = delaunayCCW(pts);
            // accept only if the triangulation reproduces the cell area: a
            // cell left slightly NON-convex by the edge contraction would have
            // its Delaunay hull spill outside, and the domain tiling check
            // downstream would fail with a far less helpful message
            double aSum = 0.0;
            for (const auto& t : tri) {
                const Vec2& A = pts[t[0]];
                const Vec2& B = pts[t[1]];
                const Vec2& C = pts[t[2]];
                aSum += 0.5 * ((B.x() - A.x()) * (C.y() - A.y())
                               - (C.x() - A.x()) * (B.y() - A.y()));
            }
            double aCell = polyArea(poly);
            if (!tri.empty() && std::abs(aSum - aCell) <= 1e-9 * aCell) {
                for (const auto& t : tri)
                    T.tri.push_back({{all[t[0]], all[t[1]], all[t[2]]}, (int)g});
            } else {
                // fall back to a fan over the SUBDIVIDED boundary: still
                // conforming with the neighbours, just star-shaped here
                ++nFallback;
                int cenId = (int)T.vtx.size();
                T.vtx.push_back(cen);
                for (std::size_t k = 0; k < loop.size(); ++k) {
                    int a = loop[k], b = loop[(k + 1) % loop.size()];
                    const Vec2& A = T.vtx[a];
                    const Vec2& B = T.vtx[b];
                    double ar = 0.5 * ((A.x() - cen.x()) * (B.y() - cen.y())
                                       - (B.x() - cen.x()) * (A.y() - cen.y()));
                    if (ar <= 0.0)
                        throw std::runtime_error("Tessellation: cell is neither "
                            "Delaunay-meshable nor star-shaped — reduce "
                            "vertexMergeFrac");
                    T.tri.push_back({{cenId, a, b}, (int)g});
                }
            }
        }
        std::printf("[tess] delaunay grain mesh: target element %.4g m "
                    "(%.2f x grain), %d cells, %ld fell back to fan%s\n",
                    h, h / targetD, (int)cells.size(), nFallback,
                    randomInterior ? " — points interieurs ALEATOIRES (Poisson-disc)" : "");
    }

    // ---- 6. conforming refinement -------------------------------------------
    for (int r = 0; r < refine; ++r) {
        std::map<std::pair<int, int>, int> mid;
        auto midpoint = [&](int a, int b) {
            auto key = std::minmax(a, b);
            auto it = mid.find({key.first, key.second});
            if (it != mid.end()) return it->second;
            int id = (int)T.vtx.size();
            T.vtx.push_back(0.5 * (T.vtx[a] + T.vtx[b]));
            mid[{key.first, key.second}] = id;
            return id;
        };
        std::vector<Tri> fine;
        fine.reserve(4 * T.tri.size());
        for (const Tri& t : T.tri) {
            int a = t.v[0], b = t.v[1], c = t.v[2];
            int ab = midpoint(a, b), bc = midpoint(b, c), ca = midpoint(c, a);
            fine.push_back({{a, ab, ca}, t.grain});
            fine.push_back({{ab, b, bc}, t.grain});
            fine.push_back({{ca, bc, c}, t.grain});
            fine.push_back({{ab, bc, ca}, t.grain});
        }
        T.tri = std::move(fine);
    }

    // ---- final validity check ------------------------------------------------
    double Asum = 0.0;
    for (const Tri& t : T.tri) {
        const Vec2& A = T.vtx[t.v[0]];
        const Vec2& B = T.vtx[t.v[1]];
        const Vec2& C = T.vtx[t.v[2]];
        double ar = 0.5 * ((B.x() - A.x()) * (C.y() - A.y())
                           - (C.x() - A.x()) * (B.y() - A.y()));
        if (ar <= 0.0)
            throw std::runtime_error("Tessellation: inverted triangle after "
                                     "refinement");
        Asum += ar;
    }
    if (std::abs(Asum - W * H) > 1e-6 * W * H)
        throw std::runtime_error("Tessellation: triangulation does not tile "
                                 "the domain");
    return T;
}

} // namespace rockim
