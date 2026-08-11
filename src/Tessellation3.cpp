// ---------------------------------------------------------------------------
// Tessellation3 — 3D Voronoi grains + phases. See the header for the
// pipeline; comments here carry the geometry. Every hard-won 2D lesson
// (ring coverage proof, edge contraction instead of spatial clustering,
// strict tiling check as backstop) is applied verbatim in 3D.
// ---------------------------------------------------------------------------
#include "rockim/Tessellation3.hpp"

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

using Vec3 = Eigen::Vector3d;

// Convex polyhedron as vertex list + outward-oriented face loops (vertex
// indices, ordered CCW seen from outside). Kept convex by construction:
// we only ever intersect with half-spaces.
struct Poly {
    std::vector<Vec3> v;
    std::vector<std::vector<int>> f;
};

Poly makeBox(double W, double D, double H) {
    Poly P;
    P.v = {{0, 0, 0}, {W, 0, 0}, {W, D, 0}, {0, D, 0},
           {0, 0, H}, {W, 0, H}, {W, D, H}, {0, D, H}};
    // outward loops (CCW from outside)
    P.f = {{0, 3, 2, 1},    // z = 0
           {4, 5, 6, 7},    // z = H
           {0, 1, 5, 4},    // y = 0
           {2, 3, 7, 6},    // y = D
           {0, 4, 7, 3},    // x = 0
           {1, 2, 6, 5}};   // x = W
    return P;
}

// Keep the part with n . p <= b. Face loops are clipped Sutherland-Hodgman
// style (order preserved, so outward orientation survives); each cut edge
// contributes ONE crossing vertex, keyed on the sorted old edge pair so the
// two faces sharing the edge reference the same new vertex. The cap face
// (on the plane) is the convex polygon of the crossing vertices, ordered by
// angle about its centroid and oriented along +n (the cell's outward
// direction there).
Poly clipHalfSpace(const Poly& P, const Vec3& n, double b) {
    const int nv = (int)P.v.size();
    std::vector<double> d(nv);
    bool anyIn = false, anyOut = false;
    for (int i = 0; i < nv; ++i) {
        d[i] = n.dot(P.v[i]) - b;
        (d[i] <= 0.0 ? anyIn : anyOut) = true;
    }
    if (!anyOut) return P;                        // untouched
    if (!anyIn)
        throw std::runtime_error("Tessellation3: seed cell vanished "
                                 "(coincident seeds?)");

    Poly Q;
    std::vector<int> keep(nv, -1);
    auto keptId = [&](int i) {
        if (keep[i] < 0) { keep[i] = (int)Q.v.size(); Q.v.push_back(P.v[i]); }
        return keep[i];
    };
    std::map<std::pair<int, int>, int> crossId;
    auto crossing = [&](int a, int bb) {
        auto key = std::minmax(a, bb);
        auto it = crossId.find({key.first, key.second});
        if (it != crossId.end()) return it->second;
        double t = d[a] / (d[a] - d[bb]);
        int id = (int)Q.v.size();
        Q.v.push_back(P.v[a] + t * (P.v[bb] - P.v[a]));
        crossId[{key.first, key.second}] = id;
        return id;
    };

    std::vector<int> capIds;
    for (const auto& loop : P.f) {
        std::vector<int> out;
        const int k = (int)loop.size();
        for (int i = 0; i < k; ++i) {
            int a = loop[i], bb = loop[(i + 1) % k];
            bool inA = d[a] <= 0.0, inB = d[bb] <= 0.0;
            if (inA) out.push_back(keptId(a));
            if (inA != inB) {
                int c = crossing(a, bb);
                out.push_back(c);
                capIds.push_back(c);
            }
        }
        // drop consecutive duplicates a grazing cut can produce
        std::vector<int> clean;
        for (int id : out)
            if (clean.empty() || clean.back() != id) clean.push_back(id);
        while (clean.size() > 1 && clean.back() == clean.front())
            clean.pop_back();
        if (clean.size() >= 3) Q.f.push_back(std::move(clean));
    }

    // cap face: unique crossings, angle-ordered about the centroid in the
    // cutting plane, oriented so the loop normal points along +n
    std::sort(capIds.begin(), capIds.end());
    capIds.erase(std::unique(capIds.begin(), capIds.end()), capIds.end());
    if (capIds.size() >= 3) {
        Vec3 cen = Vec3::Zero();
        for (int id : capIds) cen += Q.v[id];
        cen /= (double)capIds.size();
        Vec3 nn = n.normalized();
        Vec3 u = nn.cross(std::abs(nn.x()) < 0.9 ? Vec3(1, 0, 0)
                                                 : Vec3(0, 1, 0)).normalized();
        Vec3 w = nn.cross(u);
        std::sort(capIds.begin(), capIds.end(), [&](int a, int bb) {
            Vec3 pa = Q.v[a] - cen, pb = Q.v[bb] - cen;
            return std::atan2(pa.dot(w), pa.dot(u))
                 < std::atan2(pb.dot(w), pb.dot(u));
        });
        Vec3 pn = Vec3::Zero();                    // polygon normal (Newell)
        for (std::size_t i = 0; i < capIds.size(); ++i)
            pn += Q.v[capIds[i]].cross(Q.v[capIds[(i + 1) % capIds.size()]]);
        if (pn.dot(nn) < 0.0) std::reverse(capIds.begin(), capIds.end());
        Q.f.push_back(std::move(capIds));
    }
    return Q;
}

// Signed volume and centroid: face-fan tets with the vertex mean as apex
// (outward loops make every signed volume positive for a convex cell).
void polyVolCentroid(const Poly& P, double& vol, Vec3& cen) {
    Vec3 apex = Vec3::Zero();
    for (const Vec3& p : P.v) apex += p;
    apex /= (double)P.v.size();
    vol = 0.0;
    cen = Vec3::Zero();
    for (const auto& loop : P.f)
        for (std::size_t i = 1; i + 1 < loop.size(); ++i) {
            const Vec3& a = P.v[loop[0]];
            const Vec3& b = P.v[loop[i]];
            const Vec3& c = P.v[loop[i + 1]];
            double v6 = (a - apex).dot((b - apex).cross(c - apex));
            vol += v6;
            cen += v6 * (apex + a + b + c);
        }
    vol /= 6.0;
    cen = vol > 1e-30 ? cen / (24.0 * vol) : apex;
}

// Voronoi cells of all seeds inside the box, by half-space clipping with the
// distance-sorted early-out and the ring coverage proof (2 R_cell <= ring *
// grid cell — checking the early-out on the candidates alone is NOT enough,
// the 2D +0.23 % overlap bug).
std::vector<Poly> voronoiCells3(const std::vector<Vec3>& seeds,
                                double W, double D, double H) {
    const int N = (int)seeds.size();
    std::vector<Poly> cells(N);

    double cs = std::cbrt(W * D * H / std::max(1, N));  // ~1 seed per cell
    int gx = std::max(1, (int)(W / cs));
    int gy = std::max(1, (int)(D / cs));
    int gz = std::max(1, (int)(H / cs));
    std::vector<std::vector<int>> grid((std::size_t)gx * gy * gz);
    auto cellOf = [&](const Vec3& p, int& cx, int& cy, int& cz) {
        cx = std::clamp((int)(p.x() / W * gx), 0, gx - 1);
        cy = std::clamp((int)(p.y() / D * gy), 0, gy - 1);
        cz = std::clamp((int)(p.z() / H * gz), 0, gz - 1);
    };
    for (int i = 0; i < N; ++i) {
        int cx, cy, cz;
        cellOf(seeds[i], cx, cy, cz);
        grid[((std::size_t)cz * gy + cy) * gx + cx].push_back(i);
    }
    double csMin = std::min({W / gx, D / gy, H / gz});

    std::vector<std::pair<double, int>> cand;
    for (int i = 0; i < N; ++i) {
        int ci, cj, ck;
        cellOf(seeds[i], ci, cj, ck);
        Poly poly;

        for (int ring = 1;; ++ring) {
            cand.clear();
            for (int dk = -ring; dk <= ring; ++dk)
                for (int dj = -ring; dj <= ring; ++dj)
                    for (int di = -ring; di <= ring; ++di) {
                        int cx = ci + di, cy = cj + dj, cz = ck + dk;
                        if (cx < 0 || cy < 0 || cz < 0
                            || cx >= gx || cy >= gy || cz >= gz) continue;
                        for (int j : grid[((std::size_t)cz * gy + cy) * gx + cx])
                            if (j != i)
                                cand.push_back(
                                    {(seeds[j] - seeds[i]).squaredNorm(), j});
                    }
            std::sort(cand.begin(), cand.end());

            poly = makeBox(W, D, H);           // re-clip from scratch (cheap)
            double r2max = 0.0;
            for (const Vec3& v : poly.v)
                r2max = std::max(r2max, (v - seeds[i]).squaredNorm());
            for (std::size_t k = 0; k < cand.size(); ++k) {
                if (0.25 * cand[k].first > r2max) break;  // sorted: none closer
                int j = cand[k].second;
                Vec3 n = seeds[j] - seeds[i];
                double b = n.dot(0.5 * (seeds[i] + seeds[j]));
                poly = clipHalfSpace(poly, n, b);
                if (poly.f.size() < 4)
                    throw std::runtime_error("Tessellation3: seed cell "
                                             "degenerated (coincident "
                                             "seeds?)");
                r2max = 0.0;
                for (const Vec3& v : poly.v)
                    r2max = std::max(r2max, (v - seeds[i]).squaredNorm());
            }
            bool coveredAll = (ci - ring < 0 && cj - ring < 0 && ck - ring < 0
                               && ci + ring >= gx && cj + ring >= gy
                               && ck + ring >= gz);
            double rcov = ring * csMin;
            if (coveredAll || 4.0 * r2max <= rcov * rcov) break;
        }
        cells[i] = std::move(poly);
    }
    return cells;
}

} // namespace

Tessellation3 Tessellation3::build(double W, double D, double H,
                                   double targetD, double jitter, int lloyd,
                                   double mergeFrac, int refine,
                                   const std::vector<double>& phaseFraction,
                                   std::mt19937& rng, bool randomSeeds) {
    double Lmin = std::min({W, D, H});
    if (targetD <= 0.0 || targetD > 0.5 * Lmin)
        throw std::runtime_error("Tessellation3: grainSize must be positive "
                                 "and much smaller than the domain");
    if (refine < 0 || refine > 2)
        throw std::runtime_error("Tessellation3: refineLevels must be in 0..2 "
                                 "in 3D (each level multiplies the element "
                                 "count by 8)");

    // ---- 1. seeds ------------------------------------------------------------
    // HCP per-seed volume a^3/sqrt(2) matched to the sphere volume pi/6 d^3
    double a = targetD * std::cbrt(M_PI * std::sqrt(2.0) / 6.0);
    double dyRow = a * std::sqrt(3.0) / 2.0;
    double dzLay = a * std::sqrt(2.0 / 3.0);
    std::vector<Vec3> seeds;
    double eps = 1e-6 * Lmin;
    if (randomSeeds) {
        // Poisson-disc dart throwing, isotropic boundary orientations
        double dmin = 0.7 * a;
        long target = std::lround(W * D * H / (a * a * a / std::sqrt(2.0)));
        double cs = dmin;
        int gx = std::max(1, (int)(W / cs));
        int gy = std::max(1, (int)(D / cs));
        int gz = std::max(1, (int)(H / cs));
        std::vector<std::vector<int>> acc((std::size_t)gx * gy * gz);
        std::uniform_real_distribution<double> Ux(eps, W - eps);
        std::uniform_real_distribution<double> Uy(eps, D - eps);
        std::uniform_real_distribution<double> Uz(eps, H - eps);
        long attempts = 0, maxAttempts = 400 * std::max(1L, target);
        while ((long)seeds.size() < target && attempts < maxAttempts) {
            ++attempts;
            Vec3 p(Ux(rng), Uy(rng), Uz(rng));
            int ci = std::clamp((int)(p.x() / W * gx), 0, gx - 1);
            int cj = std::clamp((int)(p.y() / D * gy), 0, gy - 1);
            int ck = std::clamp((int)(p.z() / H * gz), 0, gz - 1);
            bool ok = true;
            for (int dk = -1; dk <= 1 && ok; ++dk)
                for (int dj = -1; dj <= 1 && ok; ++dj)
                    for (int di = -1; di <= 1 && ok; ++di) {
                        int cx = ci + di, cy = cj + dj, cz = ck + dk;
                        if (cx < 0 || cy < 0 || cz < 0
                            || cx >= gx || cy >= gy || cz >= gz) continue;
                        for (int q : acc[((std::size_t)cz * gy + cy) * gx + cx])
                            if ((seeds[q] - p).norm() < dmin) {
                                ok = false;
                                break;
                            }
                    }
            if (!ok) continue;
            acc[((std::size_t)ck * gy + cj) * gx + ci]
                .push_back((int)seeds.size());
            seeds.push_back(p);
        }
    } else {
        // jittered HCP lattice (ABAB), the 3D reading of the jittered hex
        std::uniform_real_distribution<double> U(-1.0, 1.0);
        int layer = 0;
        for (double z = 0.5 * dzLay; z < H; z += dzLay, ++layer) {
            double ox = (layer % 2) * 0.5 * a;
            double oy = (layer % 2) * a / (2.0 * std::sqrt(3.0));
            int row = 0;
            for (double y = 0.5 * dyRow + oy; y < D; y += dyRow, ++row) {
                double x0 = 0.5 * a + ox + (row % 2) * 0.5 * a;
                for (double x = x0; x < W; x += a) {
                    Vec3 p(x + jitter * 0.5 * a * U(rng),
                           y + jitter * 0.5 * dyRow * U(rng),
                           z + jitter * 0.5 * dzLay * U(rng));
                    p.x() = std::clamp(p.x(), eps, W - eps);
                    p.y() = std::clamp(p.y(), eps, D - eps);
                    p.z() = std::clamp(p.z(), eps, H - eps);
                    seeds.push_back(p);
                }
            }
        }
    }
    if (seeds.size() < 2)
        throw std::runtime_error("Tessellation3: fewer than 2 grains — "
                                 "reduce grainSize");

    // ---- 2. Voronoi cells (+ Lloyd relaxation) -------------------------------
    std::vector<Poly> cells;
    for (int it = 0;; ++it) {
        cells = voronoiCells3(seeds, W, D, H);
        if (it >= lloyd) break;
        for (std::size_t i = 0; i < seeds.size(); ++i) {
            double vol;
            Vec3 c;
            polyVolCentroid(cells[i], vol, c);
            seeds[i].x() = std::clamp(c.x(), eps, W - eps);
            seeds[i].y() = std::clamp(c.y(), eps, D - eps);
            seeds[i].z() = std::clamp(c.z(), eps, H - eps);
        }
    }

    Tessellation3 T;
    T.nGrains = (int)cells.size();

    // ---- 3. short-EDGE contraction (the sliver guard) ------------------------
    // 3a. exact dedup of the per-cell vertex copies into diagram vertices:
    // each cell computed its vertices independently (different clip chains),
    // so the SAME diagram vertex exists as several fp copies a hair apart.
    double tol = std::max(mergeFrac * targetD, 1e-12);
    double epsF = 1e-9 * Lmin;

    std::vector<Vec3> vFine;
    std::unordered_map<long long, std::vector<int>> hash;
    auto hkey = [&](int ix, int iy, int iz) {
        return ((long long)ix * 73856093LL) ^ ((long long)iy * 19349663LL)
             ^ ((long long)iz * 83492791LL);
    };
    auto fineId = [&](const Vec3& p) {
        int ix = (int)std::floor(p.x() / (4.0 * epsF));
        int iy = (int)std::floor(p.y() / (4.0 * epsF));
        int iz = (int)std::floor(p.z() / (4.0 * epsF));
        for (int dk = -1; dk <= 1; ++dk)
            for (int dj = -1; dj <= 1; ++dj)
                for (int di = -1; di <= 1; ++di) {
                    auto it = hash.find(hkey(ix + di, iy + dj, iz + dk));
                    if (it == hash.end()) continue;
                    for (int id : it->second)
                        if ((vFine[id] - p).norm() < epsF) return id;
                }
        int id = (int)vFine.size();
        vFine.push_back(p);
        hash[hkey(ix, iy, iz)].push_back(id);
        return id;
    };
    // per-cell faces on fine ids
    std::vector<std::vector<std::vector<int>>> cellFine(cells.size());
    for (std::size_t g = 0; g < cells.size(); ++g) {
        for (const auto& loop : cells[g].f) {
            std::vector<int> ids;
            for (int lv : loop) {
                int id = fineId(cells[g].v[lv]);
                if (ids.empty() || ids.back() != id) ids.push_back(id);
            }
            while (ids.size() > 1 && ids.back() == ids.front()) ids.pop_back();
            if (ids.size() >= 3) cellFine[g].push_back(std::move(ids));
        }
        if (cellFine[g].size() < 4)
            throw std::runtime_error("Tessellation3: degenerate cell "
                                     "polyhedron");
    }

    // 3b. union-find over short diagram edges (edges = consecutive pairs in
    // the face loops)
    std::vector<int> parent(vFine.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::function<int(int)> find = [&](int x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    for (const auto& fl : cellFine)
        for (const auto& ids : fl)
            for (std::size_t k = 0; k < ids.size(); ++k) {
                int p = ids[k], q = ids[(k + 1) % ids.size()];
                if ((vFine[p] - vFine[q]).norm() < tol) {
                    int rp = find(p), rq = find(q);
                    if (rp != rq) parent[std::max(rp, rq)] = std::min(rp, rq);
                }
            }

    // 3c. component positions: mean of members, snapped onto any of the six
    // box walls a member lies on
    std::unordered_map<int, int> rootToNew;
    std::vector<Vec3> sum;
    std::vector<int> cnt, flags;
    auto wallOf = [&](const Vec3& p) {
        return (p.x() < epsF ? 1 : 0) | (p.x() > W - epsF ? 2 : 0)
             | (p.y() < epsF ? 4 : 0) | (p.y() > D - epsF ? 8 : 0)
             | (p.z() < epsF ? 16 : 0) | (p.z() > H - epsF ? 32 : 0);
    };
    for (int vI = 0; vI < (int)vFine.size(); ++vI) {
        int r = find(vI);
        auto [it, isNew] = rootToNew.try_emplace(r, (int)sum.size());
        if (isNew) {
            sum.push_back(Vec3::Zero());
            cnt.push_back(0);
            flags.push_back(0);
        }
        int id = it->second;
        sum[id] += vFine[vI];
        cnt[id] += 1;
        flags[id] |= wallOf(vFine[vI]);
    }
    T.vtx.resize(sum.size());
    for (std::size_t id = 0; id < sum.size(); ++id) {
        Vec3 p = sum[id] / cnt[id];
        if (flags[id] & 1)  p.x() = 0.0;
        if (flags[id] & 2)  p.x() = W;
        if (flags[id] & 4)  p.y() = 0.0;
        if (flags[id] & 8)  p.y() = D;
        if (flags[id] & 16) p.z() = 0.0;
        if (flags[id] & 32) p.z() = H;
        T.vtx[id] = p;
    }

    // 3d. rebuild the face loops on the contracted vertices; faces that
    // collapse below a triangle are gone (the neighbours close over them),
    // pinched loops (a repeated non-consecutive vertex) are rejected loudly
    std::vector<std::vector<std::vector<int>>> cellFaces(cells.size());
    for (std::size_t g = 0; g < cells.size(); ++g) {
        for (const auto& fl : cellFine[g]) {
            std::vector<int> ids;
            for (int fid : fl) {
                int id = rootToNew.at(find(fid));
                if (ids.empty() || ids.back() != id) ids.push_back(id);
            }
            while (ids.size() > 1 && ids.back() == ids.front()) ids.pop_back();
            if (ids.size() < 3) continue;
            std::vector<int> sorted = ids;
            std::sort(sorted.begin(), sorted.end());
            if (std::adjacent_find(sorted.begin(), sorted.end())
                != sorted.end())
                throw std::runtime_error("Tessellation3: pinched face after "
                                         "edge contraction — reduce "
                                         "vertexMergeFrac");
            cellFaces[g].push_back(std::move(ids));
        }
        if (cellFaces[g].size() < 4)
            throw std::runtime_error("Tessellation3: a grain degenerated "
                                     "during edge contraction — reduce "
                                     "vertexMergeFrac");
    }

    // ---- 4. tet meshing ------------------------------------------------------
    // Face centroids are SHARED between the two cells of a face (keyed on
    // the sorted vertex set), so both sides fan-triangulate identically and
    // the FDEM face registry pairs their tets into joints. Cell centroid +
    // face-centroid fan cones every face triangle into a tet.
    // Contraction can leave a shared face slightly NON-PLANAR, and a
    // non-planar polygon has no unique volume attribution: fanning it from
    // different apexes (as the two adjacent cells would) books different
    // volumes and the tiling check false-alarms at ~1e-4 relative. The
    // face-centroid triangulation below IS shared by both cells, so all
    // volume bookkeeping (grainVol, tiling check) runs on the actual tets.
    std::map<std::vector<int>, int> faceCen;
    T.grainVol.assign(cells.size(), 0.0);
    double volSum = 0.0;
    for (std::size_t g = 0; g < cells.size(); ++g) {
        // cell centroid estimate from the face loops (fan volumes; the
        // slight non-planarity bias is irrelevant for a fan apex)
        Vec3 apex = Vec3::Zero();
        int na = 0;
        for (const auto& ids : cellFaces[g])
            for (int id : ids) { apex += T.vtx[id]; ++na; }
        apex /= (double)na;
        double vol = 0.0;
        Vec3 cen = Vec3::Zero();
        for (const auto& ids : cellFaces[g])
            for (std::size_t i = 1; i + 1 < ids.size(); ++i) {
                const Vec3& A = T.vtx[ids[0]];
                const Vec3& B = T.vtx[ids[i]];
                const Vec3& C = T.vtx[ids[i + 1]];
                double v6 = (A - apex).dot((B - apex).cross(C - apex));
                vol += v6;
                cen += v6 * (apex + A + B + C);
            }
        vol /= 6.0;
        if (vol <= 0.0)
            throw std::runtime_error("Tessellation3: non-positive cell "
                                     "volume after contraction — reduce "
                                     "vertexMergeFrac");
        cen /= 24.0 * vol;

        int cenId = (int)T.vtx.size();
        T.vtx.push_back(cen);
        for (const auto& ids : cellFaces[g]) {
            std::vector<int> key = ids;
            std::sort(key.begin(), key.end());
            auto [it, isNew] = faceCen.try_emplace(key, (int)T.vtx.size());
            if (isNew) {
                Vec3 fc = Vec3::Zero();
                for (int id : ids) fc += T.vtx[id];
                T.vtx.push_back(fc / (double)ids.size());
            }
            int fcId = it->second;
            for (std::size_t k = 0; k < ids.size(); ++k) {
                int p = ids[k], q = ids[(k + 1) % ids.size()];
                T.tet.push_back({{cenId, fcId, p, q}, (int)g});
                const Vec3& A = T.vtx[cenId];
                const Vec3& B = T.vtx[fcId];
                const Vec3& C = T.vtx[p];
                const Vec3& E = T.vtx[q];
                double tv = std::abs((B - A).dot((C - A).cross(E - A))) / 6.0;
                T.grainVol[g] += tv;
            }
        }
        volSum += T.grainVol[g];
    }
    // Tolerance note: after edge contraction the shared faces are slightly
    // non-planar and a handful of fan tets can locally overlap their
    // neighbours by |v| bookkeeping — a sub-1e-4 residue at the recommended
    // vertexMergeFrac (0.25), against the 2.3e-3 overlap the 2D clipping
    // bug produced. 1e-4 still catches real tiling errors by two orders.
    if (std::abs(volSum - W * D * H) > 1e-4 * W * D * H) {
        char msg[256];
        std::snprintf(msg, sizeof(msg),
                      "Tessellation3: cells do not tile the domain after "
                      "welding (sum %.9g vs %.9g, mismatch %.3g = %.2e rel)",
                      volSum, W * D * H, volSum - W * D * H,
                      (volSum - W * D * H) / (W * D * H));
        throw std::runtime_error(msg);
    }

    // ---- 5. phase assignment (volume-greedy on shuffled grains) --------------
    std::vector<double> frac = phaseFraction;
    if (frac.empty()) frac = {1.0};
    double fsum = std::accumulate(frac.begin(), frac.end(), 0.0);
    if (fsum <= 0.0)
        throw std::runtime_error("Tessellation3: phase fractions must sum "
                                 "to a positive value");
    for (double& f : frac) f /= fsum;

    std::vector<int> order(cells.size());
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rng);
    std::vector<double> assigned(frac.size(), 0.0);
    T.phaseOfGrain.assign(cells.size(), 0);
    for (int g : order) {
        int best = 0;
        double bestDef = -1e300;
        for (std::size_t p = 0; p < frac.size(); ++p) {
            double def = frac[p] * volSum - assigned[p];
            if (def > bestDef) { bestDef = def; best = (int)p; }
        }
        T.phaseOfGrain[g] = best;
        assigned[best] += T.grainVol[g];
    }

    // ---- 6. conforming 1:8 refinement ----------------------------------------
    // Corner tets at the 4 vertices + the interior octahedron cut into 4
    // tets around the (m01, m23) diagonal. Edge midpoints are keyed on the
    // sorted vertex pair, so faces split 1:4 identically on both sides of
    // every grain boundary and the joints stay paired.
    for (int r = 0; r < refine; ++r) {
        std::map<std::pair<int, int>, int> mid;
        auto midpoint = [&](int p, int q) {
            auto key = std::minmax(p, q);
            auto it = mid.find({key.first, key.second});
            if (it != mid.end()) return it->second;
            int id = (int)T.vtx.size();
            T.vtx.push_back(0.5 * (T.vtx[p] + T.vtx[q]));
            mid[{key.first, key.second}] = id;
            return id;
        };
        std::vector<Tet> fine;
        fine.reserve(8 * T.tet.size());
        for (const Tet& t : T.tet) {
            int v0 = t.v[0], v1 = t.v[1], v2 = t.v[2], v3 = t.v[3];
            int m01 = midpoint(v0, v1), m02 = midpoint(v0, v2);
            int m03 = midpoint(v0, v3), m12 = midpoint(v1, v2);
            int m13 = midpoint(v1, v3), m23 = midpoint(v2, v3);
            fine.push_back({{v0, m01, m02, m03}, t.grain});
            fine.push_back({{v1, m01, m12, m13}, t.grain});
            fine.push_back({{v2, m02, m12, m23}, t.grain});
            fine.push_back({{v3, m03, m13, m23}, t.grain});
            // octahedron m02-m12-m13-m03 equator around the m01-m23 axis
            fine.push_back({{m01, m23, m02, m12}, t.grain});
            fine.push_back({{m01, m23, m12, m13}, t.grain});
            fine.push_back({{m01, m23, m13, m03}, t.grain});
            fine.push_back({{m01, m23, m03, m02}, t.grain});
        }
        T.tet = std::move(fine);
    }

    // ---- final validity check ------------------------------------------------
    // (orientation is normalized by the solver's tet builder; here the
    // magnitude check guards degenerate slivers and the tiling backstop)
    double vSum = 0.0;
    for (const Tet& t : T.tet) {
        const Vec3& A = T.vtx[t.v[0]];
        const Vec3& B = T.vtx[t.v[1]];
        const Vec3& C = T.vtx[t.v[2]];
        const Vec3& E = T.vtx[t.v[3]];
        double v = std::abs((B - A).dot((C - A).cross(E - A))) / 6.0;
        if (v <= 0.0)
            throw std::runtime_error("Tessellation3: degenerate tet after "
                                     "meshing/refinement — reduce "
                                     "vertexMergeFrac");
        vSum += v;
    }
    if (std::abs(vSum - W * D * H) > 1e-4 * W * D * H)
        throw std::runtime_error("Tessellation3: tet mesh does not tile "
                                 "the domain");
    return T;
}

} // namespace rockim
