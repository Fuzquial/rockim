#pragma once
// ---------------------------------------------------------------------------
// Solver: common interface for the continuous (FEM) and discrete (DEM) modes.
//
// The driver (main.cpp) only sees this interface: it builds a solver from the
// shared problem definition (Config + Material + Tool), advances it in time,
// and asks it to write frames / history rows. A future v3 FDEM solver (v2 of
// the physics) implements the same interface and reuses the same setup and
// output layers — see README "Extending to FDEM".
// ---------------------------------------------------------------------------
#include <ostream>
#include <string>

namespace rockim {

class Solver {
public:
    virtual ~Solver() = default;

    virtual void init() = 0;                       // build mesh/packing, compute dt
    virtual void step() = 0;                       // one explicit time step
    virtual void writeFrame(int frame) = 0;        // VTK snapshot
    virtual void historyHeader(std::ostream&) const = 0;
    virtual void historyRow(std::ostream&) const = 0;
    virtual void finalize() = 0;                   // summary + verification checks

    // Optional early-stop hook. The driver polls it after every step and ends
    // the time loop when it returns true, then writes the last frame, the last
    // history row and the summary exactly as it would at t = T. The default is
    // false, so every existing solver runs to T with no change of behaviour.
    virtual bool finished() const { return false; }

    double time()     const { return t_; }
    double dt()       const { return dt_; }
    double duration() const { return T_; }

protected:
    double t_  = 0.0;
    double dt_ = 0.0;
    double T_  = 0.0;
};

} // namespace rockim
