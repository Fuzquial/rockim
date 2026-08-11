#pragma once
// ---------------------------------------------------------------------------
// Tool: rigid impactor / cutter, shared by both solvers.
//
// Shapes:
//   FLAT — flat-ended punch. x = midpoint of the *bottom face*; occupies
//          [x.x - width/2, x.x + width/2] x [x.y, +inf). Vertical contact only
//          (percussive flat punch).
//   DISC — circular tool (button insert / cutter). x = center, radius R.
//          Full 2D normal, so it works for both indentation and lateral cutting.
//
// Motion:
//   FREE       — rigid mass driven by contact forces (percussive impact).
//   PRESCRIBED — kinematic path at constant velocity (cutting / shearing);
//                the recorded force is then the reaction on the tool.
//
// 2D plane model: all masses/forces are per unit out-of-plane thickness.
// ---------------------------------------------------------------------------
#include <cmath>

#include <Eigen/Dense>

namespace rockim {

struct Tool {
    // PDC — polycrystalline diamond compact cutter, the tool of the drilling
    // literature (Heilman et al., ARMA 24-0238). In 2D it reduces to a rigid
    // WEDGE: a rake face through the cutting edge, inclined by the BACK RAKE
    // angle from the vertical so that its top trails behind, plus an optional
    // chamfer at the edge. `x` is the CUTTING EDGE; the cutter occupies the
    // half-space behind the rake face and above the edge, so intact rock below
    // the cut depth is never touched. Contrast with DISC, whose round profile
    // has no rake angle at all and cannot reproduce the rake dependence that
    // dominates the cutting force.
    enum class Shape  { FLAT, DISC, PDC };
    enum class Motion { FREE, PRESCRIBED };

    Shape  shape  = Shape::DISC;
    Motion motion = Motion::FREE;

    double mass   = 5.0;    // [kg per m thickness]
    double width  = 0.02;   // FLAT: face width [m]
    double radius = 0.01;   // DISC: radius [m]
    // PDC geometry
    double rakeDeg = 20.0;  // back rake angle from the vertical [deg]
    double faceLen = 0.013; // rake face extent from the edge [m]
    double chamLen = 0.0;   // chamfer length at the edge (0 = none) [m]
    double chamDeg = 45.0;  // chamfer angle from the rake face [deg]

    // Outward normal of the rake face (points INTO the rock, i.e. forward and
    // slightly up for a positive back rake) and the along-face direction.
    Eigen::Vector2d rakeNormal() const {
        double b = rakeDeg * M_PI / 180.0;
        return {std::cos(b), std::sin(b)};
    }
    Eigen::Vector2d rakeDir() const {          // edge -> top of the face
        double b = rakeDeg * M_PI / 180.0;
        return {-std::sin(b), std::cos(b)};
    }

    Eigen::Vector2d x{0.0, 0.0};  // position (see shape docs)
    Eigen::Vector2d v{0.0, 0.0};
    Eigen::Vector2d F{0.0, 0.0};  // total contact force ON the tool (this step)

    void resetForce() { F.setZero(); }

    void integrate(double dt) {
        if (motion == Motion::FREE) v += (dt / mass) * F;  // gravity neglected at impact time scales
        x += dt * v;
    }

    double ke() const { return 0.5 * mass * v.squaredNorm(); }
};

} // namespace rockim
