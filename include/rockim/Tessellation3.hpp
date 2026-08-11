#pragma once
// ---------------------------------------------------------------------------
// Tessellation3: Voronoi grain structure of the box [0,W] x [0,D] x [0,H],
// meshed into tetrahedra, with per-grain phase labels — the 3D lift of the
// 2D Tessellation and the geometric front-end of the 3D GBM mode of the
// FDEM solver (tets carry their grain id, so the solver can classify every
// triangular cohesive joint as intra-grain or grain-boundary).
//
// Pipeline (each step mirrors its 2D counterpart; the 2D bring-up lessons —
// ring coverage proof, edge contraction instead of spatial clustering,
// strict tiling check — are built in from the start):
//   1. Seeds on an HCP lattice sized for the target grain diameter,
//      jittered (fraction of the half spacing), or Poisson-disc dart
//      throwing (grainSeeding = random: isotropic boundary orientations),
//      then optionally relaxed by Lloyd iterations (seed -> cell centroid).
//   2. Voronoi cell of each seed by half-SPACE clipping of the domain box
//      against the bisector planes of its neighbours, visited in order of
//      increasing distance with the exact early-out; the neighbour rings
//      grow until the coverage PROOF holds (2 R_cell <= ring * grid cell).
//      Cells are convex polyhedra stored as vertex-index face loops
//      (outward-oriented); the cap face cut by each bisector is assembled
//      from the face intersection segments (convexity makes it a convex
//      polygon, ordered by angle about its centroid).
//   3. Short-EDGE contraction (union-find over diagram edges, component
//      mean position, wall-snapped on the six box walls): the sliver
//      guard. Contraction, NOT spatial clustering — the 2D lesson that
//      welding by proximity makes neighbouring cells overlap holds
//      verbatim in 3D; the strict volume check below is the backstop.
//   4. Tet meshing: each face is fan-triangulated from its OWN centroid
//      vertex (shared between the two cells of the face, so joints stay
//      paired), and each face triangle is coned to the cell centroid.
//      Optional conforming refinement: each level splits every tet 1:8
//      through the edge midpoints (4 corner tets + the interior octahedron
//      cut into 4 around one diagonal); face triangles split 1:4 with
//      shared midpoints, so grain-boundary faces stay paired.
//   5. Phase assignment per grain by volume-greedy matching of the target
//      fractions, exactly as in 2D.
//
// Vertex ids returned here are the "virtual" ids the FDEM mesh builder
// keys its face map on: two tets share a joint iff they reference the same
// vertex triple, across grains as well as inside them.
// ---------------------------------------------------------------------------
#include <array>
#include <random>
#include <vector>

#include <Eigen/Dense>

namespace rockim {

struct Tessellation3 {
    struct Tet {
        std::array<int, 4> v;   // virtual vertex ids
        int grain;
    };

    std::vector<Eigen::Vector3d> vtx;   // welded virtual vertices
    std::vector<Tet> tet;
    std::vector<int> phaseOfGrain;      // one entry per grain
    std::vector<double> grainVol;       // polyhedron volume per grain
    int nGrains = 0;

    // W, D, H  : domain size [m]
    // targetD  : target mean grain diameter [m]
    // jitter   : seed jitter, fraction of the half lattice spacing [0..1]
    //            (hcp seeding only)
    // lloyd    : Lloyd relaxation iterations (>= 0)
    // mergeFrac: edge-contraction tolerance, fraction of targetD
    // refine   : conforming 1:8 refinement levels (0..2 — each level
    //            multiplies the element count by 8)
    // phaseFraction: target volume fractions (empty or size 1 = single phase)
    // randomSeeds: false = jittered HCP lattice (compact, but boundary
    //            orientations keep the lattice preferences); true =
    //            Poisson-disc dart throwing (isotropic orientations,
    //            slightly wider grain-size spread)
    static Tessellation3 build(double W, double D, double H, double targetD,
                               double jitter, int lloyd, double mergeFrac,
                               int refine,
                               const std::vector<double>& phaseFraction,
                               std::mt19937& rng, bool randomSeeds = false);
};

} // namespace rockim
