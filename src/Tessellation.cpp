// ---------------------------------------------------------------------------
// Tessellation — Voronoi grains + phases. See the header for the pipeline.
// ---------------------------------------------------------------------------
#include "rockim/Tessellation.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
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
std::vector<std::vector<Vec2>> voronoiCells(const std::vector<Vec2>& seeds,
                                            double W, double H) {
    const int N = (int)seeds.size();
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
                if (0.25 * cand[k].first > r2max) break;   // sorted: none closer
                int j = cand[k].second;
                Vec2 n = seeds[j] - seeds[i];
                double b = n.dot(0.5 * (seeds[i] + seeds[j]));
                poly = clipHalfPlane(poly, n, b);
                if (poly.size() < 3)
                    throw std::runtime_error("Tessellation: seed cell vanished "
                                             "(coincident seeds?)");
                r2max = 0.0;
                for (const Vec2& v : poly)
                    r2max = std::max(r2max, (v - seeds[i]).squaredNorm());
            }
            bool coveredAll = (ci - ring < 0 && cj - ring < 0
                               && ci + ring >= gx && cj + ring >= gy);
            double rcov = ring * csMin;
            if (coveredAll || 4.0 * r2max <= rcov * rcov) break;
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
                                 bool useDelaunay, double elemSize) {
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
    double eps = 1e-6 * std::min(W, H);
    if (randomSeeds) {
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
    std::vector<std::vector<Vec2>> cells;
    for (int it = 0;; ++it) {
        cells = voronoiCells(seeds, W, H);
        if (it >= lloyd) break;
        for (std::size_t i = 0; i < seeds.size(); ++i) {
            Vec2 c = polyCentroid(cells[i]);
            seeds[i].x() = std::clamp(c.x(), eps, W - eps);
            seeds[i].y() = std::clamp(c.y(), eps, H - eps);
        }
    }

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
    std::shuffle(order.begin(), order.end(), rng);
    std::vector<double> assigned(frac.size(), 0.0);
    T.phaseOfGrain.assign(cells.size(), 0);
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
                    "(%.2f x grain), %d cells, %ld fell back to fan\n",
                    h, h / targetD, (int)cells.size(), nFallback);
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
