#pragma once
// ---------------------------------------------------------------------------
// Tessellation: Voronoi grain structure of the rectangle [0,W] x [0,H],
// fan-triangulated into CCW constant-strain triangles, with per-grain phase
// labels. This is the geometric front-end of the GBM (grain-based model)
// mode of the FDEM solver: triangles carry their grain id, so the solver can
// classify every cohesive joint as intra-grain or grain-boundary.
//
// Pipeline:
//   1. Seeds on a hexagonal lattice sized for the target grain diameter,
//      jittered uniformly (grainJitter in [0,1], fraction of the half
//      spacing), then optionally relaxed by Lloyd iterations (seed -> cell
//      centroid), which rounds the cells without erasing the disorder.
//   2. Voronoi cell of each seed by half-plane clipping of the domain
//      rectangle against the bisectors of its neighbours, visited in order
//      of increasing distance with the exact early-out (a bisector at
//      half-distance beyond the cell's circumradius cannot cut it). The
//      neighbour rings are enlarged until the coverage PROOF holds
//      (2 R_cell <= ring * grid cell): the early-out alone cannot certify
//      completeness, since a cutting seed may be absent from the ring.
//   3. Short-EDGE contraction: Voronoi diagrams of jittered lattices
//      contain very short edges; edges shorter than mergeFrac * targetD
//      are contracted transitively (union-find), the merged vertex sits at
//      the component mean and is snapped back onto any domain wall a
//      member lay on. This is the sliver guard: without it the fan
//      triangulation produces degenerate CSTs and the explicit time step
//      collapses. (Contraction, NOT spatial clustering: merging vertices
//      that are close but not edge-connected makes neighbouring cells
//      overlap — the strict tiling check at the end guards this.)
//   4. Fan triangulation of each (convex) cell from its centroid; optional
//      conforming refinement (each level splits every triangle in 4 through
//      the edge midpoints, shared across grains so joints stay paired).
//   5. Phase assignment per grain by area-greedy matching of the target
//      fractions (each grain goes to the phase with the largest remaining
//      area deficit, in shuffled order), so achieved area fractions match
//      the requested ones to within about one grain.
//
// Vertex ids returned here are the "virtual" ids the FDEM mesh builder keys
// its edge map on: two triangles share a joint iff they reference the same
// two virtual vertices, across grains as well as inside them.
// ---------------------------------------------------------------------------
#include <array>
#include <random>
#include <vector>

#include <Eigen/Dense>

namespace rockim {

// Bowyer-Watson Delaunay triangulation of a point set: CCW triangles covering
// the convex hull, deterministic in the caller's point order. Shared with the
// native disc mesher of the brazilian test (a disc is convex, so every
// returned triangle belongs to the domain and no constrained triangulation is
// needed).
std::vector<std::array<int, 3>> delaunayCCW(
    const std::vector<Eigen::Vector2d>& p);

struct Tessellation {
    struct Tri {
        std::array<int, 3> v;   // virtual vertex ids, CCW
        int grain;
    };

    std::vector<Eigen::Vector2d> vtx;   // welded virtual vertices
    std::vector<Tri> tri;
    std::vector<int> phaseOfGrain;      // one entry per grain
    std::vector<double> grainArea;      // polygon area per grain
    int nGrains = 0;

    // W, H     : domain size [m]
    // targetD  : target mean grain diameter [m]
    // jitter   : seed jitter, fraction of the half lattice spacing [0..1]
    //            (hex seeding only)
    // lloyd    : Lloyd relaxation iterations (>= 0)
    // mergeFrac: edge-contraction tolerance, fraction of targetD
    // refine   : conforming 4-way refinement levels (0..4)
    // phaseFraction: target area fractions (empty or size 1 = single phase)
    // randomSeeds: false = jittered hex lattice (compact, but boundary
    //            orientations stay quantized near the three lattice
    //            directions); true = Poisson-disc dart throwing (isotropic
    //            boundary orientations, slightly wider grain-size spread)
    // useDelaunay: false = fan triangulation from the cell centroid (the
    //            original; cheap, but the intra-grain topology is a HUB and a
    //            transgranular crack can only follow a spoke); true = the
    //            unstructured Delaunay intra-grain mesh of the grain-based
    //            FDEM literature (Y-Geo / Irazu GBM studies use about
    //            0.18 x the grain diameter, i.e. six elements across a grain,
    //            so that transgranular cracking is representable at all)
    // elemSize : target element size [m] for the delaunay front-end
    //            (<= 0 selects the literature default 0.18 * targetD)
    static Tessellation build(double W, double H, double targetD,
                              double jitter, int lloyd, double mergeFrac,
                              int refine,
                              const std::vector<double>& phaseFraction,
                              std::mt19937& rng, bool randomSeeds = false,
                              bool useDelaunay = false, double elemSize = 0.0);
};

} // namespace rockim
