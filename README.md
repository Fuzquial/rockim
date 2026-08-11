# rockim — 2D rock impact & shear-cutting simulator (FEM + DEM)

A small, self-contained C++17 code that simulates the dynamic response and
failure of a rectangular rock specimen under two loading scenarios relevant to
percussive/rotary drilling:

- **percussion** — a free rigid impactor (given mass and velocity) strikes the
  top surface;
- **shear** — a rigid cutter is dragged laterally at a fixed depth of cut and
  imposed velocity.

Two independent solvers live behind one common `Solver` interface:

- **FEM** — explicit dynamics on constant-strain triangles (plane strain),
  Drucker–Prager + Rankine failure with scalar damage, crack-band softening
  and element erosion;
- **DEM** — bonded circular particles (parallel bonds carrying normal force,
  shear force and moment), tensile/shear bond breakage, frictional contact
  between unbonded particles, emergent fragmentation.

Both solvers share the configuration format, the tool models, the output
formats and the time-stepping philosophy, which is exactly the layout you want
if the next step is a combined FDEM code (see *Roadmap*).

---

## Build

Requirements: CMake ≥ 3.16, a C++17 compiler, Eigen3 (header-only; a system
install under `/usr/include/eigen3` is found automatically).

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```

This produces a single executable `rockim`.

## Run

```bash
./rockim ../configs/fem_percussion.cfg  out_femp
./rockim ../configs/dem_percussion.cfg  out_demp
./rockim ../configs/fem_shear.cfg       out_fems
./rockim ../configs/dem_shear.cfg       out_dems

# verification cases
./rockim ../configs/verify_fem_bar.cfg      out_vbar
./rockim ../configs/verify_dem_tension.cfg  out_vdt
```

The first argument is a plain-text `key = value` config file (see `configs/`
for commented examples of every option), the optional second argument is the
output directory. Quick-look figures:

```bash
python3 ../tools/plot_results.py out_demp   # writes plot_field.png + plot_history.png
```

---

## Theory notes

### Explicit time step and stability

Both solvers use central-difference (leapfrog) explicit integration, which is
conditionally stable: `dt <= 2/omega_max`, with `omega_max` the highest
natural frequency of the discrete system.

- **FEM.** `omega_max` is bounded through the smallest element: the code uses
  `dt = CFL * h_min / c_p`, where `h_min` is the minimum element inscribed
  size and `c_p = sqrt(E(1-nu) / (rho(1+nu)(1-2nu)))` is the plane-strain
  dilatational wave speed (the fastest signal in the mesh). The contact
  penalty stiffness adds a local spring `k_p` on surface nodes; the code also
  checks `dt <= CFL * 2 sqrt(m_node/k_p)` so the penalty never controls
  stability silently. Default `CFL = 0.7`.

- **Absorbing boundaries (Lysmer–Kuhlemeyer).** With `absorbing = sides`
  (lateral faces) or `absorbing = all` (lateral faces + bottom, which then
  replaces the fixed support / rigid wall), outgoing waves are absorbed by
  viscous dashpots matched to the continuum impedances: normal traction
  `rho c_p v_n`, tangential `rho c_s v_t`, lumped over each boundary node's
  (FEM) or boundary particle's (DEM) tributary length. The dashpot is applied
  *implicitly* in the velocity update, `v <- (v + dt F/m)/(1 + dt c/m)`, so it
  is unconditionally stable and never affects the time step. In the DEM the
  dashpots act only while a particle still carries at least one intact bond:
  the quiet boundary stands in for the truncated elastic continuum and must
  not drag on detached debris flying through the layer. Effect on the demos
  (both percussion cases use the *same* insert: a free rigid disc, R = 15 mm,
  5 kg, 8 m/s, so the two modes are directly comparable): the spurious
  reflected-wave damage disappears (FEM percussion: elements with D > 0.1
  drop from 1127 to 82, erosion unchanged; DEM percussion: 1113 -> 716
  broken bonds and the large boundary-driven fragments vanish), while the
  contact-force peak is unchanged. Default is
  `none` (the verification cases need the original boundary conditions).

  With the identical insert, the two modes still answer differently — peak
  0.55 MN/m, ~1 mm penetration, no rebound within the window, 143 J/m
  absorbed (FEM) versus peak 2.7 MN/m, 0.22 mm penetration, restitution
  ~0.79, 61 J/m absorbed (DEM, H = 0.2 block). That gap is model form, not a
  bug: the FEM damage law caps the contact pressure and keeps dissipating in
  the crushed pocket, while the uncalibrated hex-lattice DEM is much stronger
  under the confined contact than its nominal bond strengths suggest. Closing
  it is precisely the DEM calibration exercise described under *Honest
  limitations*.

  A quiet boundary absorbs *waves*; it cannot stand in for the missing rock
  around a *crack* tip, whose driving field is static. The percussion demo
  documents this: on an H = 0.1 block the median crack (which the bond break
  times show running top-down at ~3 km/s) reached the bottom face and picked
  up extra breaks there, while on an H = 0.2 block the same impact arrests
  the crack on its own at 0.101 m depth and leaves the bottom half pristine
  (322 vs 544 broken bonds). Hence the rule of thumb used by the shipped
  configs: keep every quiet boundary at least one expected crack length away
  from the process zone, and treat any fracture touching a boundary as a
  domain-size red flag, not a result.

- **DEM.** The classical estimate `dt = f * 2 sqrt(m/k)` with a single bond
  stiffness is *not* safe once a particle accumulates many bonds plus several
  frictional contacts (exactly what happens in the crushed zone under the
  tool). The code therefore sums, per particle, the translational stiffness of
  all attached bonds (`k_n A + k_s A`) plus a budget of extra contacts, and
  the rotational stiffness (`k_n I_b + k_s A (L/2)^2`), and takes
  `dt = f * min_i ( 2 sqrt(m_i / K_trans,i), 2 sqrt(I_i / K_rot,i) )` with
  `f = dtFactor = 0.2` by default. This fixed an energy-pumping instability
  observed with the naive estimate (tool bouncing back with more energy than
  it brought).

### Failure representation

- **FEM.** Effective stress is checked against a Rankine cut-off (tension,
  `ft`) and a Drucker–Prager cone matched to Mohr–Coulomb in plane strain
  (cohesion `c`, friction angle `phi`). Violation drives a scalar damage
  `D in [0,1]` with exponential softening whose slope is regularized by the
  crack-band method: the fracture energy `Gf` (and `gfShearFactor * Gf` in
  shear) is smeared over the element size `l_c = sqrt(2A)`, so the dissipated
  energy per unit crack area is mesh-independent. Nominal stress is
  `(1-D) sigma_eff`; elements with `D >= 0.98` (or a strain cap) are eroded —
  removed from the force computation — which is how macro-cracks, craters and
  chips appear.

- **DEM.** Particles are glued by parallel bonds: incremental springs in the
  normal direction, tangential direction and in rotation. A bond breaks when
  the maximum tensile stress at the bond periphery
  `sigma_max = F_n/A + |M_b| R_b / I_b` exceeds `sigma_c`, or when the shear
  stress `|F_s|/A` exceeds `c_bar + tan(phi_b) * sigma_n` (compression
  raising the shear resistance). Broken pairs fall back to frictional
  Hertz-like (linear here) contact with Coulomb friction and a stored
  tangential spring. Cracks are simply the set of broken bonds; fragments are
  connected components of the surviving bond network, extracted at the end of
  the run.

---

## Outputs (per run directory)

| file | content |
|---|---|
| `history.csv` | time, tool force (N/m), tool position/velocity, work, kinetic energy, broken bonds / eroded elements, damage volume |
| `frame_XXXX.vtu` | VTK unstructured-grid frames (ParaView-ready): FEM = mesh with `damage`, `eroded`, velocity; DEM = points with radius, velocity, fragment id, plus line cells for intact bonds |
| `fem_final_elements.csv` / `dem_final_particles.csv`, `dem_final_bonds.csv` | final-state dumps for scripting |
| `dem_fragments.csv` | fragment id, particle count, mass, volume per fragment |
| `frames.csv` | frame index → time and tool pose (used by `make_gif.py`) |
| `summary.txt` | the end-of-run summary (peak force, work, specific energy, …) |

ParaView tips: for DEM `.vtu` frames, apply **Glyph → Sphere**, scale by the
`radius` array (factor 2 = diameter) to see the packing; color by
`fragment` or `speed`. Bonds are rendered as line cells — color them by
`broken` on the final frame, or just threshold them away.

`tools/plot_results.py` produces the same views without ParaView, and
`tools/make_gif.py <config.cfg> <run_dir> [out.gif]` renders an animated GIF
of the whole run with the rigid tool drawn (DEM: particles coloured by speed +
broken bonds in red; FEM: damage field with eroded elements as holes).

---

## Verification

| case | measured | expected | error |
|---|---|---|---|
| FEM bar wave speed (`nu = 0`, `verify_fem_bar.cfg`) | 4382.3 m/s | `sqrt(E/rho)` = 4343.7 m/s | **+0.89 %** |
| DEM direct tension, load-aligned square lattice (`verify_dem_tension.cfg`) | 9.99 MPa peak | calibrated `sigma_c` = 10 MPa | **−0.11 %** |

The tension case uses a square lattice on purpose: only load-aligned bonds
carry stress, so the macroscopic strength must equal the bond strength — a
sharp, parameter-free check of the bond force/failure implementation. On the
hexagonal lattices used for the demos, the macro strength is an emergent
property (lattice geometry, disorder) and requires calibration, as in any DEM.

Energy bookkeeping is done on the tool side (`W = -∫ F_tool · v_tool dt`),
so for a free impactor the reported work matches its kinetic-energy loss to
within a few percent in all demos.

---

## The 3D DEM (`mode = dem3d`)

The bonded-particle model also exists in three dimensions, behind the same
`Solver` interface and config format. What it is:

* **Packing.** Equal spheres on an HCP lattice (ABAB stacking): every sphere
  is bonded to its 12 neighbours at exactly 2r — 6 in its layer, 3 below,
  3 above. The tension verification instead uses a load-aligned simple-cubic
  lattice so the expected macroscopic strength is known in closed form.
* **Parallel bonds, full 3D.** Each bond of radius Rb = lambda*r carries a
  normal force (tension positive), a shear-force *vector* in the plane
  perpendicular to the bond axis, a twisting moment about the axis and a
  bending-moment *vector*: A = pi Rb^2, I = pi Rb^4/4, J = 2I, per-area
  modulus E/L0 and shear ratio `ksRatio`. Failure by maximum tensile fibre
  stress sigma = Fn/A + |Mb| Rb/I > sigma_c, or by shear
  tau = |Fs|/A + |Mt| Rb/J exceeding c + tan(phi) * (compressive normal
  stress). The stored perpendicular vectors are kept perpendicular by
  projection each step — exact for the small per-step rotations of an
  explicit run. Force and moment application mirrors the energy-consistent
  2D scheme (shear couple -(L/2) n x Fs on both particles).
* **Contacts.** Broken pairs (and any overlapping pair) interact through a
  linear normal spring with viscous damping and a *vector* tangential spring
  with history, rotated into the contact plane each step and capped by
  Coulomb friction; torques use the true lever arms. Sub-nanometre overlaps
  are ignored (resting HCP neighbours sit at exactly 2r, and floating-point
  noise would otherwise churn the tangential-spring map for zero physics).
* **Tool.** A rigid sphere: free flight for percussion (button insert),
  prescribed path for shear (dragged cutter at a depth of cut). Percussion
  alternatively takes `toolShape = flat`: a flat-ended cylindrical punch of
  radius `toolRadius` (axis z), the 3D lift of the 2D FLAT tool.
* **Tension extras.** `pullRamp` (cosine rise of the grip velocity) and
  `gripLateralFree` (frictionless grips) work as in the 2D FDEM.
* **Quiet boundaries.** The same Lysmer + Deeks-Randolph viscous-spring
  treatment as in 2D, applied on the four lateral faces and the bottom
  (`absorbing = all`), lumped over the close-packed tributary area
  2 sqrt(3) r^2 and gated on bondedness.
* **Outputs.** `dem3d_particles_XXXX.vtu` / `dem3d_bonds_XXXX.vtu` (open in
  ParaView and add a Glyph > Sphere scaled by the `radius` array for the
  true 3D view; threshold the bonds on `state` for the crack surface),
  `dem3d_final_particles.csv`, `dem3d_fragments.csv`, `history.csv` in true
  SI units (N, J, m^3). `make_gif.py` and `plot_results.py` render a
  mid-depth slice for a quick look without ParaView.

**3D verification.** Direct tension on the load-aligned simple-cubic lattice:
only the z-bonds carry load, so the macroscopic peak stress must equal
(A_bond/A_cell) sigma_c = (pi/4) lambda^2 sigma_c = 7.854 MPa with the
defaults. Measured: 7.853 MPa, error -0.02 % (`verify_dem3d_tension.cfg`,
PASS printed by `finalize()`).

**3D-specific limitations.** The regular HCP lattice needs the same
calibration caveat as the 2D hex lattice, only more so (12-fold coordination,
preferred crack planes); a disordered packing generator is the natural
upgrade. The per-step cost is ~50 ms for 39k particles / 225k bonds
single-threaded (bond sweep and contact sweep in roughly equal parts), so the
shipped demo runs in ~7 minutes; OpenMP over particles is the obvious next
step. For reference, that demo (R = 15 mm sphere, 0.5 kg, 8 m/s, 16 J on a
0.1 x 0.1 x 0.08 m block) gives a 93.6 kN force peak at 93 us, 0.55 mm of
penetration, restitution 0.96, 2632 broken bonds and 96 small fragments
around the indent — the stiff, quasi-elastic answer expected from an
uncalibrated close-packed lattice, consistent with the 2D comparison above.

## Honest limitations

- 2D plane strain: out-of-plane relief and 3D fragment shapes are absent;
  2D wave trapping tends to over-fragment compared to 3D.
- FEM: small-strain kinematics; the DP/Rankine + scalar-damage coupling is a
  deliberately simple model (no plastic flow rule, no dilatancy); erosion
  removes mass and its pattern keeps some mesh bias despite crack-band
  regularization.
- DEM: quantitative macro-properties (UCS, tensile strength, friction) of a
  hexagonal packing must be *calibrated* against lab tests before any
  engineering claim; there is no re-bonding, no moisture/rate effects.
- The reported "detached volume" counts every fragment other than the
  largest, so a single through-going crack that splits the block in two makes
  it jump — read it together with the fragment-size table.
- Rigid tools only; no tool wear, no confining pressure (both are natural
  extensions).
- The Lysmer boundaries absorb normally-incident waves exactly and oblique
  incidence only approximately; a few percent of the radiated energy still
  reflects (perfect absorption would need PML-type layers).

## The combined FDEM (`mode = fdem`) — the v2, delivered

The Munjiza-architecture FDEM sketched in the roadmap below now exists as a
fourth mode, and the tension verification plus a three-way percussion
comparison back it up.

**The model.** Every CST owns its three nodes (node duplication); elements
are linear elastic and CO-ROTATIONAL (closed-form 2D polar decomposition), so
detached fragments tumble without spurious straining, with an
elastic-perfectly-plastic ceiling on the deviator (`crushCap`, default
8*cohesion) bounding the energy storable by crushed elements. 4-node cohesive
joints sit on every interior edge from the start (intrinsic approach,
penalty `jointPenaltyFactor`*E/h, default 20): mode I softens linearly from
ft over an opening of 2 Gf/ft, mode II from c over 2 Gf_II/c of frictional
slip (return mapping), friction tan(phi)*(-sigma_n) acts throughout
compression and survives as the residual strength of broken joints; damage is
the shared max of the two drivers. A broken joint stays alive as the paired
frictional contact of its own faces and is released to the general contact
ONLY on clear separation. Freed faces and the exterior enter a node-edge
penalty contact on a cell grid with three safeguards that the debugging
session below made necessary; fragments are connected components of cohesive
joints; the viscous-spring quiet boundaries are shared with the other modes.

**Verification.** Direct tension on a uniform strip: every joint carries the
same stress, so the macroscopic peak must equal ft exactly — the joint
penalty softens the modulus, never the strength. Measured 9.93 MPa for
ft = 10 MPa (-0.7 %), clean single-plane crack (`verify_fdem_tension.cfg`).

**A debugging story worth keeping.** The first percussion run pulverised the
whole block: 21k broken joints from a 148 J/m impact, 1e6 J/m of kinetic
energy — an energy pump. The hunt, by bisection with kill-switches and a
net-work meter on the general contact, isolated three real defects, each
fixed and each a known failure mode of penalty-contact FDEM:
(1) half the exterior edges had inverted outward normals (they skipped the
CCW re-orientation the joints got), turning them into permanent phantom
pushers; (2) joints crushed in shear died while still compressed, handing
interpenetrated faces to the general contact, whose penalty then released
1/2 k pen^2 of energy created from nothing — joints now die by clear
separation only, and the contact got an initial-penetration relief
(birth-gap, LS-DYNA-style); (3) the node-segment penalty spring on fast
ROTATING debris faces is a follower force (the penetration is geometric, the
work is kinematic) and pumps at full stiffness — debris contact is now soft
(0.01 E t), heavily damped (xi 0.8) and quasi-plastic (restitution 0.2 on
release, physically sensible for crushed rock). After the fixes the tool
energy balance closes to 0.03 % and the block's residual kinetic energy is
70 J/m for 139 J/m delivered.

**Three-way percussion comparison** (same insert everywhere: free rigid
disc R = 15 mm, 5 kg, 8 m/s, 0.2 x 0.2 m block, viscous-spring boundaries):

| mode | peak force | penetration | rebound | absorbed | breakage |
|------|-----------|-------------|---------|----------|----------|
| FEM (damage) | 0.55 MN/m | ~1 mm | none in window | 143 J/m | 22 eroded elems |
| FDEM | 1.23 MN/m | 0.39 mm | 2.9 m/s (e~0.36) | 139 J/m | 253 joints, 160 fragments |
| DEM (BPM) | 2.68 MN/m | 0.22 mm | 6.3 m/s (e~0.79) | 61 J/m | 322 bonds |

The FDEM lands between the two, as its construction predicts: continuum
elasticity like the FEM (so no lattice over-strength), discrete fracture
surfaces with an explicit Gf (so stiffer and less dissipative than the
smeared damage + erosion of the FEM). Specific energy of the FDEM crater:
0.33 MJ/m^3.

**FDEM-specific limitations.** Structured cross-diagonal mesh (jittered by
`meshJitter`) biases crack paths to mesh edges — an unstructured mesher is
the natural upgrade; 2-point joint integration and linear softening instead
of Munjiza's z(D); node-edge penalty contact instead of potential contact;
the quasi-plastic debris contact trades exact restitution for robustness;
`gcWork` printed at the end is a debugging meter (relative-velocity work),
judge energy sanity on the tool balance and block KE it reports alongside.
Finally, the intrinsic penalty makes the elastic joint threshold microscopic
(ft/p ~ 30 nm here), so high-frequency joint vibrations ratchet a diffuse
sub-critical damage field across the block (most joints end with
0.1 < D < 1); the D = 1 crack set, the fragment count and the energy balance
are the robust outputs, and an extrinsic (insert-on-criterion) variant would
remove the artefact at the cost of topology changes at run time.

## Voronoi grains + mineral phases (`mesh = voronoi`) — the GBM mode

The FDEM's structured cross-diagonal mesh has been joined by a grain-based
front-end: a Voronoi tessellation of the block whose cells are grains, whose
grains carry mineral phases, and whose cohesive joints know whether they sit
*inside* a grain or *on a grain boundary*. Crack paths then follow the
boundary network instead of the lattice directions — the upgrade the
FDEM-specific limitations paragraph called for, plus the mineralogy.

**Geometry.** Seeds on a hexagonal lattice sized for `grainSize` (target mean
diameter), jittered by `grainJitter`, relaxed by `lloydIters` Lloyd
iterations; Voronoi cells by half-plane clipping with the exact
sorted-by-distance early-out; **short-edge contraction** with tolerance
`vertexMergeFrac * grainSize` collapses the short Voronoi edges that would
otherwise produce sliver triangles and kill the stable time step. Cells are
fan-triangulated from their centroids; `refineLevels` conforming 4-way
refinements set the intra-grain crack resolution. Everything is reproducible
from `seed`, and a strict tiling check (sum of cell areas = W x H to 1e-6)
guards the whole chain.

Two geometry bugs were caught by that check during bring-up, and both are
worth keeping in mind for anyone touching the code: (1) the neighbour-ring
early-out is not a completeness proof — a cutting seed can be absent from
the ring block entirely, which left some cells unclipped and OVERLAPPING
(+0.23 % area on the percussion domain); the rings now grow until
2 R_cell <= ring x grid-cell, the actual coverage guarantee. (2) welding
vertices by spatial proximity is wrong — an interior vertex and a wall
vertex separated by a thin boundary cell got merged, the cells swept across
the cell in between, and the boundary dented; the weld is now a transitive
short-EDGE contraction (union-find over diagram edges, component mean
position, wall-snapped), which only merges vertices whose incident cells
all follow the move.

**Phases.** `phases = quartz feldspar biotite` declares the set; each phase
overrides the global material through `phase.<name>.<prop>` keys (E, nu, rho,
ft, cohesion, frictionDeg, Gf, gfShearFactor) and must give a `fraction`.
Grains are assigned by area-greedy matching, so achieved area fractions track
the targets to within about one grain (measured: < 0.1 %). Elements carry
their phase in the element loop (per-phase D-matrix, crush cap, mass), the
quiet boundaries use the impedances of their local phase, and the stable-dt
scan uses each node's own phase stiffness and each joint's own penalty.

**Joints.** Intra-grain joints (type 0) carry the bulk properties of their
phase. Grain-boundary joints take the *mean of the two neighbouring phases*
times attenuation factors — `gbAlphaTen`, `gbAlphaCoh`, `gbAlphaGf`,
`gbAlphaE`, `gbAlphaFric` (the alpha coefficients of the GBM calibration
literature; all default 1) — and heterophase boundaries (type 2, different
minerals) get one extra `gbHeteroFactor` on the strength-like properties.
The summary reports the joint census and, after a run, the
**intergranular fraction** of the broken joints — the observable a GBM
exists to produce. The `.vtu` frames carry `phase` and `grain` per element
and `type` per joint for ParaView.

**Verification.**

| case | measured | expected | verdict |
|---|---|---|---|
| single-phase voronoi tension, all alphas = 1 (`verify_fdem_voronoi_tension.cfg`) | 11.20 MPa peak (seeds 7 / 42: 10.91 / 10.85) | ft = 10 MPa + path tortuosity (each inclined facet sees sigma cos^2 theta, so the peak must sit AT or somewhat ABOVE ft) | **+12.0 %, PASS** (band −5…+25 %; +8.5…+12 % across three seeds) |
| two-phase, weak boundaries gbAlphaTen = 0.4 (`fdem_voronoi_tension_gb.cfg`) | 4.44 MPa = 0.44 ft | ~0.4 ft (flat weak path) + tortuosity | **consistent**, and the crack is **100 % intergranular** (15/15 broken joints on boundaries) |

The second line is the GBM signature: strength and crack path are controlled
by the boundary network, not by the bulk — grid FDEM cannot represent that
at all. (`verifyFt = false` in that config replaces the PASS/FAIL-against-ft
line, which would be noise there, by the ratio report.)

**Percussion demo** (`fdem_voronoi_percussion.cfg`): the same rigid disc as
every other percussion demo (R = 15 mm, 5 kg, 8 m/s) on a three-phase
granite-like microstructure (33 % quartz / 59 % feldspar / 8 % biotite,
textbook contrasts, Poisson-disc seeding, boundaries at half strength,
heterophase 0.8). Shipped-seed numbers: 509 grains / 11 220 elements
(refine 1), achieved phase fractions within 0.1 % of target, peak force
2.25 MN/m, 160 J/m absorbed with the tool energy balance closing to
0.01 %, 40 broken joints (17.5 % intergranular), 19 fragments,
~4.5 min single-threaded. Color by `grain` or `phase` and threshold the
joint `damage` in ParaView for the faciès. All demo numbers are
UNCALIBRATED illustrations (textbook phase properties, guessed alphas) —
calibration against lab data is the whole subject of the GBM literature.

**Seeding and isotropy.** `grainSeeding = hex` (default) places seeds on a
jittered hexagonal lattice: compact, uniform grain sizes — but measurably
ANISOTROPIC: at the shipped defaults the grain-boundary orientations stay
quantized near the three lattice directions (~30/90/150 deg; the two
60-degree sectors around 60/120 deg carry < 1 % of the boundary length
against 33 % for an isotropic Voronoi), because Lloyd relaxation
re-converges toward the honeycomb. Crack paths inherit that bias — the very
artefact the GBM mode is meant to remove. `grainSeeding = random`
(Poisson-disc dart throwing, used by the shipped percussion demo) gives
isotropic boundary orientations at the price of a wider grain-size spread.
Use `random` for any result where crack directions matter; `hex` remains
useful for controlled parameter studies.

**Statistical joint strengths (`jointWeibullM`) — the implicit-DFH bridge.**
With `jointWeibullM = m` (> 1), every joint's ft and cohesion are multiplied
by a Weibull(m) factor of mean 1; fracture energies stay deterministic and
the softening openings are recomputed. `strengthCorrLength` picks the
spatial structure: `0` draws independently per joint (the analogue of the
per-element Weibull draw in smeared DFH codes — statistics converge with
the mesh, the crack MAP does not), while `> 0` samples ONE Gaussian random
field of that correlation length through the Gaussian copula
(`RandomField.hpp`, moving-average construction, exactly N(0,1) marginals).
The field lives in space, independent of the mesh — `fieldSeed` controls it
separately from the mesh `seed` — so two different meshes see the same weak
zones: this is the sandbox for the "correlated sigma_w as initial
condition" idea, whose whole point is a reproducible crack map at
continuum-level cost. The joints `.vtu` carries `ftScale` for ParaView.
The crack-map experiment (same field, several meshes — do they break in
the same place?) needed two loading fixes before the field could speak,
both now in the demo config (`fdem_voronoi_tension_weibull.cfg`):
(1) `gripLateralFree = true` — with fully clamped grips the Saint-Venant
corner concentration decides where the specimen breaks; (2) `pullRamp` —
a STEPPED grip velocity launches a transient that unzips the first joint
row under the grip whatever the strength map says (a straight, unphysical
crack at the very top); a smooth cosine ramp over ~10 wave transits makes
the loading quasi-static. With both fixes the demonstration lands:
THREE different meshes (seeds 12345/111/222) on the SAME field
(fieldSeed 555, weakest band at mid-height) all break through that band
(mean crack y = 0.065/0.064/0.069 on H = 0.12, spread below the
correlation length) at consistent peaks (4.6-5.1 MPa), while the
independent-draw control breaks much STRONGER (7.4 MPa: uncorrelated
draws offer no continuous weak path — the correlated field is what
creates band-like weakness). A reproducible crack map from a
mesh-independent strength field is precisely the property the
"correlated sigma_w as initial condition" idea needs — demonstrated here
at sandbox scale. The unzipping artefact itself is a faithful
loading-regime effect (the wave beats the heterogeneity), worth keeping
in mind for any velocity-driven test on heterogeneous specimens.

**OpenMP.** The FDEM hot loops are parallel when built with OpenMP
(`/openmp`, or CMake finds it automatically): element forces (node
duplication makes them race-free), joint forces (per-thread accumulation,
reduced in thread order), tool contact (ordered partial sums) and the
general-contact SEARCH (per-thread candidate lists above a 4096-active-node
threshold — below it fork/join costs more than the sweep; the birth-gap
bookkeeping and force application stay serial by design). Guarantees:
OMP_NUM_THREADS=1 is BIT-IDENTICAL to the serial build; a fixed thread
count is deterministic run-to-run; different thread counts differ at
floating-point associativity level, which the chaotic fragmentation phase
amplifies into realization-level spread (same class as a seed change).
Measured on the shipped voronoi percussion demo (18 logical hybrid cores):
267 s serial -> 180 s, a modest 1.5x — the remaining serial parts (contact
application, fragment bookkeeping, output) and the hybrid-core stragglers
bound it; per-step hot loops profile at ~1.5 ms (ROCKIM_PROF=1 to see).

**Bayesian calibration bench (`tools/bayes_bench.py`).** An end-to-end
SYNTHETIC test of the emulator-based calibration strategy (LHS design ->
rockim runs -> one GP per observable -> Metropolis-Hastings), on the
two-phase GB tension case with a hidden truth theta* = (gbAlphaTen,
gbAlphaCoh) and mesh-seed scatter as the experimental noise. The physics
makes the expected verdict sharp — the mode-I tension peak identifies
alphaTen and barely constrains alphaCoh — so the bench checks BOTH that
the pipeline recovers the truth and that it exposes non-identifiability
(posterior ~ prior on the weak component) instead of hiding it.
Self-contained numpy/matplotlib (scipy optional). This is the cheap dress
rehearsal of the SHPB-style calibration campaign: validate the machinery
on runs that cost seconds before spending the real FDEM budget.
Shipped-seed verdict (15 rockim runs, ~5 min): gbAlphaTen recovered —
median 0.428, CI95 [0.396, 0.458] around the hidden 0.45,
posterior/prior sd 0.07; gbAlphaCoh correctly flagged NOT identified
(CI95 spans the prior, ratio 1.06); R-hat <= 1.007, GP leave-one-out
coverage 0.92-1.00. PASS.

**General-contact note.** Voronoi faces can be several `hmin` long, so the
debris contact now bins every active edge into ALL grid cells its AABB
covers (with node-edge pair dedup) instead of the midpoint cell only, and
the deep-penetration cap uses the local element size on the voronoi mesh.
For the grid mesh the detected pair set is unchanged, but the summation
order shifts at bit level, which the chaotic fragmentation phase amplifies:
grid percussion debris statistics move within realization-to-realization
noise (measured: far less than a seed change moves them; the force peak is
bit-identical).

**Config validation.** The GBM keys are strictly validated: unknown `mesh`
values, `phases` with `mesh = grid` or with non-FDEM modes, zero or negative
strengths/stiffnesses (a zero-strength joint would become UNBREAKABLE
through the cohesive-envelope arithmetic — model pre-cracked boundaries
with a small alpha like 1e-3, never 0), friction angles >= 89 deg,
`refineLevels` outside 0..4, and numeric values with trailing garbage
(a French-locale decimal comma reads as garbage now instead of silently
truncating: `0,5` used to parse as 0.0) all abort at startup with a named
key. Reproducibility from `seed` is guaranteed per binary; across standard
libraries (MSVC vs libstdc++) `std::shuffle` and the distributions may
legitimately produce different draws for the same seed.

**GBM-specific limitations, stated plainly.** Fan triangulation fixes the
intra-grain mesh topology (a centroid hub); grain-shape statistics are those
of a relaxed Voronoi, not of a real granite texture (no elongated micas);
phases are elastic contrasts + strength contrasts only (no cleavage
anisotropy inside grains); the mean-times-alpha boundary rule is the
simplest of the literature's options. `mesh = voronoi` is wired to the
FDEM solvers only (`fem`/`dem` modes reject it) — the 3D lift
(Tessellation3 + tet front-end) now exists in `fdem3d`, see below.

## The 3D continuum FEM with pluggable laws (`mode = fem3d`)

The sandbox counterpart of an Abaqus/Explicit + VUMAT percussion run:
shared-node Kuhn tetrahedra, co-rotational, rigid spherical impactor,
Lysmer quiet boundaries — and the CONSTITUTIVE LAW as a config switch
(`law = elastic | dpr | saksala | saksala2011 | dpdfh`), so the same 3D
impact can be replayed under different rock models. `dpdfh` is the
line-by-line port of the thesis' own DP-DFH VUMAT (see the MatLaw header
for the card): anisotropic DFH obscuration damage with deterministic
per-element Weibull draws, verified at 4.7e-12 against the ifx-compiled
Fortran reference (`rockim selftest-dpdfh`, four paths incl. an oblique
frozen frame; the spatial-hash draws are bit-identical). Fracture is smeared damage + element erosion
(no cohesive joints): the FEM way of making a crater. Scenarios:
percussion (sphere, or `toolShape = flat` for the flat-ended punch),
shear (dragged spherical cutter at `cutDepth`/`cutSpeed`, the 3D lift of
the 2D fem shear — `fem3d_shear.cfg`), tension/compression.

* **`dpr`** — rate-independent elastoplastic damage: Drucker-Prager cone
  (alpha, k matched to Mohr-Coulomb in triaxial compression from
  cohesion/frictionDeg, deviatoric return, psi = 0) + Rankine tensile
  damage with crack-band exponential softening (Gf over the element size;
  init aborts beyond the snapback limit E Gf/ft^2).
* **`saksala`** — rate-DEPENDENT damage-viscoplasticity in the spirit of
  Saksala's percussive-drilling model: same cone but PERZYNA overstress
  (viscosity `saksalaEta`, closed-form return dlam = F/(G + eta/dt)), a
  compression CAP on the pressure (`capP0`, hardening `capH`) for the
  confined crush under the tool, same tensile damage. Simplifications
  stated plainly: Perzyna instead of the consistency formulation, one
  scalar (tensile) damage variable, no viscosity on damage, no dilatancy.
* **Unilateral damage, learned the hard way**: a scalar (1-D) on the full
  stress also kills the COMPRESSIVE bearing capacity — the damaged surface
  collapsed and the tool TUNNELED through the block at constant velocity.
  Damage now degrades only the tensile principal components (spectral
  split), and erosion removes only fully damaged elements in net TENSION
  (spalled chips) or over-crushed ones (`erodeEpv` on the viscoplastic
  strain): a damaged element under compression stays — it is the rubble
  bed. This mirrors percussion-VUMAT practice (keep damaged elements
  alive; deletion reserved for spall).

**Verification** (mid-specimen stress gauge — the grip force additionally
carries the Cundall-damping drag of the flowing column, measured +11 % at
damping 0.7, which has nothing to do with the law):

| case | measured | expected | verdict |
|---|---|---|---|
| uniaxial compression, dpr (`verify_fem3d_dp.cfg`) | 107.35 MPa | DP analytic k/(1/sqrt3 - alpha) = 107.23 MPa | **+0.11 %, PASS** |
| direct tension, dpr (`verify_fem3d_tension.cfg`) | 10.000 MPa | ft = 10 MPa | **-0.002 %, PASS** |
| Perzyna overstress at 8.3 /s (`verify_fem3d_rate1.cfg`) | 2.83 MPa | 2.75 MPa closed form | **ratio 1.03, PASS** |
| same at 16.7 /s (`verify_fem3d_rate2.cfg`) | 5.63 MPa | 5.50 MPa | **ratio 1.02**, rate2/rate1 = **1.99: linear in rate** |

Two measurement lessons are baked into those configs: the overstress gauge
averages the mid-specimen stress over the LAST QUARTER of the run (the
end-state snapshot of a ringing signal moves with the mesh), and
`dampingLocal = 0.1` — heavy Cundall damping (0.7) acts as a body force
along the flowing column and biased the measured overstress by up to
+/- 35 % depending on the mesh.

**Mesh disorder (`meshMirror`, `meshJitter`).** The plain Kuhn split
threads ONE global family of diagonals through the mesh, and without
regularization crack patterns snap to it (measured on the cylinder
percussion: star-shaped damage arms along the diagonals). `meshMirror`
(default ON) uses the MIRRORED Kuhn split — cells reflected by index
parity, face-compatible by construction — which alternates the diagonal
directions checkerboard-wise; combined with `meshJitter` (0.45 in the
cylinder demo) the damage pattern becomes an irregular isotropic disc
with no preferred directions. Distinct radial macro-cracks would need
finer meshes and/or Weibull strength fields — `matWeibullM` (the VUMAT's
FIELD mechanism) draws a mean-1 strength factor per element, either
independently (the paper's choice) or from ONE correlated random field
(`strengthCorrLength`/`strengthCorrLengthB`/`strengthCorrAngleDeg`/
`fieldSeed`, RandomField3 — same keys and semantics as the FDEM joint
statistics). Set `meshMirror = false` for the previous mesh, bit-exactly.

**Three-law percussion comparison** (shipped `fem3d_percussion_*.cfg`:
sphere R = 15 mm, 0.5 kg, 8 m/s, 0.1 x 0.1 x 0.08 m, 62k tets, ~60 s each):

| law | peak force | penetration | restitution | absorbed | eroded |
|---|---|---|---|---|---|
| elastic | 65.2 kN | 0.64 mm | 0.99 | 0.4 J | 0 |
| dpr | 54.1 kN | 0.84 mm | 0.74 | 7.3 J | 106 |
| saksala | 55.3 kN | 0.86 mm | 0.67 | 8.8 J | 58 (cap pc 250 -> 973 MPa under the tool) |
| saksala2011 (Table I material) | 50.9 kN | 0.95 mm | 0.41 | 13.3 J | 0 by design (cap pc 1040 -> 2176 MPa) |

The faithful law is the most dissipative of all and its F-delta curve has a
different grammar: an EARLY knee and a long low plateau (strong cohesion
softening h_DP = -10 GPa, where dpr/saksala are perfectly plastic), then
the stiff cap-controlled rise with its characteristic oscillations, and
the largest hysteresis loop.

The ordering is the physical one: yielding caps the force below the
elastic peak; the rate-dependent cap law dissipates the most (confined
crush + viscous work) and rebounds the least. All demo parameters are
plausible orders of magnitude, UNCALIBRATED.

**The faithful Saksala (2011) law (`law = saksala2011`).** Beyond the
simplified `saksala`, this is a LINE-BY-LINE PORT of the thesis'
`vumat_saksala_2011.f90` (Saksala, IJNAMG 35 (2011) 1483-1505): DP cone
with non-associated potential, modified Rankine (associated), parabolic
cap joined to the cone, viscoplastic CONSISTENCY law with two viscosities
(compression `skSdp`, tension `skSmr`), confinement-dependent cohesion
degradation (`skNd`, confinement frozen at the elastic trial state),
logarithmic volumetric cap hardening, exponential tensile damage driven by
the positive norm of the viscoplastic strain increment, unilateral crack
closure, Koiter corner rule. Deliberately NO erosion and NO
fracture-energy regularization — like the published law (mesh-sensitive by
design). All `sk*` parameters default to Table I of the paper (converted
to SI); E, nu, phi, ft, cohesion come from the shared material block.
**Verification by trace superposition**: `rockim selftest-saksala2011`
replays the three material-point loading paths of the VUMAT's own test
harness (tension/modified-Rankine + damage, uniaxial compression/DP,
hydrostatic compression/cap hardening) and writes the same CSV as the
ifx-compiled Fortran reference — measured agreement **8e-14 relative** on
all 551 trace points (stresses, damage, hardening variables, cap
pressure). The shipped percussion demo `fem3d_percussion_saksala2011.cfg`
runs Table I material under the standard insert.

**Abaqus cross-validation exporter (`tools/export_abaqus.py`).** Writes the
frame-0 mesh of a fem3d run (or a fdem3d run, duplicated nodes WELDED — a
Voronoi/GBM mesh becomes a plain C3D4 mesh with one ELSET per mineral
phase) as an Abaqus .inp in the thesis' mm-t-s-MPa units, with
bottom/top NSETs and, when a `ftScale` field exists (matWeibullM), the
element-centroid CSV plus a nodal `*INITIAL CONDITIONS, TYPE=FIELD,
VARIABLE=1` include — the SAME correlated strength field sampled by both
codes, which is what the rockim <-> Abaqus/VUMAT comparison on identical
meshes (and the correlated-sigma_w structural validation) needs.
Self-check: total volume of the reordered (positive-Jacobian) tets.

**GUI (`tools/rockim_gui.py`).** A tkinter front-end (stdlib + matplotlib
only): browse and edit the shipped configs, launch runs with a live log
(exe picker, output dir, OMP thread count, stop button), one-click
verification suite with a PASS/FAIL summary, and a results pane that plots
the force-penetration curve, every history.csv column, or the mid-plane
slice (fem3d tets and 2D fdem alike) for any out_* directory.

**fem3d limitations.** Small-strain co-rotational kinematics (the crush
zone reaches ~100 % equivalent plastic strain — inside erodeEpv territory
by design); structured tet mesh only (no 3D Voronoi/GBM here: grains +
phases live in the FDEM solvers, 2D and 3D); a single material (no
`phases`); erosion loses mass, and its pattern keeps mesh bias despite
the crack band — the same honest caveats as the 2D fem module, in 3D.

## The 3D FDEM (`mode = fdem3d`)

The 2D FDEM lifted to tetrahedra, with every 2D-learned contact safeguard
built in from the first line: structured hex grid split into 6 Kuhn
tetrahedra per cell (compatible face diagonals across cells, optional
jitter), per-tet node duplication, CO-ROTATIONAL linear tets (R from F by 3
Higham iterations, Biot strain, deviatoric crush cap + mean-tension cap),
triangular 6-node cohesive joints on every interior face with 3 node-pair
integration points (mode I ft/Gf envelope, mode II vector slip return
mapping with c and Gf_II, shared damage, death by clear separation only),
general node-triangle contact (clipped-grid, co-location exclusion by origin
vertex, birth-gap relief, quasi-plastic soft normal law), rigid spherical
tool (or `toolShape = flat`, the flat-ended punch, in percussion),
viscous-spring quiet boundaries on the five truncated faces, and the
tension extras `pullRamp` + `gripLateralFree` of the 2D FDEM.

**The 3D GBM (`mesh = voronoi`) — grains, phases and joint statistics in
3D.** The whole 2D grain-based front-end exists in 3D behind the same
config keys. `Tessellation3` builds a 3D Voronoi tessellation of the box
(seeds on a jittered HCP lattice or by Poisson-disc dart throwing
[`grainSeeding = hex | random`], Lloyd relaxation, half-SPACE clipping
with the ring coverage PROOF of the 2D bring-up, transitive short-EDGE
contraction with wall snapping on the six box walls, strict volume-tiling
check), then meshes every cell by coning its faces — each fan-triangulated
from its own SHARED centroid vertex, so grain-boundary faces triangulate
identically on both sides and pair into joints — to the cell centroid.
`refineLevels` (0..2 in 3D) splits every tet 1:8 conformingly. Grains
carry mineral `phases` (volume-greedy assignment); joints are classified
intra-grain / homophase / heterophase and take the mean-times-alpha rule
(`gbAlphaTen/Coh/Gf/E/Fric`, `gbHeteroFactor`) exactly as in 2D; the
summary reports the achieved volume fractions and the intergranular
fraction of the breakage; the `.vtu` frames carry `phase`, `grain`,
joint `type` and `ftScale`. `jointWeibullM` statistical joint strengths
work too, independent per joint or sampling ONE correlated Gaussian
random field (`RandomField3`; `strengthCorrLength`, anisotropic
`strengthCorrLengthB` across a texture plane tilted by
`strengthCorrAngleDeg` about the y-axis — the 3D foliation; `fieldSeed`
independent of the mesh seed, so different meshes see the same weak
zones). Two 3D-specific notes, learned bringing it up: (1) contraction
can leave a shared face slightly NON-PLANAR, so all volume bookkeeping
runs on the shared face-centroid tets (fanning a non-planar polygon from
the two cells' different apexes books different volumes and false-alarms
the tiling check); (2) the thin fan tets of elongated Voronoi faces set
the stable dt — `vertexMergeFrac = 0.25` is the practical 3D setting
(0.12 in 2D; beyond ~0.3 the distortion trips the tiling check), and the
general contact bins voronoi faces into ALL grid cells their AABB covers
with the deep-penetration cap on the local element size, as in 2D.

**3D verification.** The Kuhn subdivision keeps complete planes of
triangular joints normal to z at every cell interface, all carrying exactly
the macroscopic stress, so the tension peak must equal ft. Measured
9.74 MPa for 10 (-2.6 %), and the specimen broke on exactly one such plane:
200 broken joints = one full 10x10-cell interface (`verify_fdem3d_tension`).
On the voronoi mesh the crack must instead follow a tortuous grain-boundary
surface whose inclined facets see sigma cos^2 theta, so the peak must sit
AT or somewhat ABOVE ft (band -5..+25 %, as in 2D): measured 10.89 MPa for
10 (+8.8 %), with 67 % of the broken joints on grain boundaries — the GBM
signature in 3D (`verify_fdem3d_voronoi_tension`, 122 grains, ~17 min
single-threaded — the joint penalty of the thinnest fan tets sets dt).

**Cost.** ~35 ms/step single-threaded for 36k tets / 70k joints / 144k
nodes; the shipped percussion demo (sphere R = 15 mm, 0.5 kg, 8 m/s on a
0.08 x 0.08 x 0.06 m block) runs in ~15 minutes. The element, joint, tool
and integration loops are OpenMP-parallel with the 2D guarantees
(OMP_NUM_THREADS=1 bit-identical to the serial build, a fixed thread
count deterministic run-to-run); the general-contact sweep stays serial. ParaView reads the
`fdem3d_*.vtu` tet mesh and joint triangles directly (threshold the joints
on `damage` for the crack surface); `make_gif.py`/`plot_results.py` render
mid-depth slices.

## Roadmap notes (what v1 was, and the Munjiza mapping)

**What v1 is, in the literature's terms.** The `fem` module is an explicit
continuum-damage code with crack-band regularisation (Bazant & Oh 1983) —
the *smeared* side of fracture. The `dem`/`dem3d` modules are bonded-particle
models with parallel bonds in the sense of Potyondy & Cundall (2004), i.e.
PFC-style BPMs. **Neither is Munjiza's FDEM**, and the difference is
structural, not cosmetic:

* *Where elasticity lives.* BPM: in the bond network, so the macroscopic
  E, nu emerge from the lattice and must be calibrated (and a regular
  lattice is anisotropic at short wavelengths). FDEM: inside deformable
  finite elements meshing every body, so E and nu are direct inputs.
* *How rock breaks.* BPM: point bonds fail instantly at a strength
  criterion — the fracture energy is whatever the lattice happens to
  dissipate. FDEM: cohesive joint elements are inserted between **all**
  element faces from the start and soften following the combined single and
  smeared crack model (Munjiza, Andrews & White 1999), so mode I is driven
  by ft *and* Gf, mode II by c and phi, and cracks are real new surfaces
  with an explicit fracture energy.
* *How bodies touch.* This code: point penalty contacts between spheres
  (and the analytic tool). FDEM: distributed potential contact forces
  between overlapping element pairs (Munjiza & Andrews 2000), penalty-based
  but smooth and resultant-consistent for arbitrary polygons/polyhedra.
* *How contacts are found.* The uniform cell grid used here is the one
  place v1 is close in spirit to Munjiza's linear NBS detection
  (Munjiza & Andrews 1998).

**The v2 path is exactly the Munjiza architecture.** Keep the CST meshes of
the `fem` module; insert 4-node cohesive joints on every internal edge with
penalty stiffness ~10-100 E/h (intrinsic insertion adds artificial
compliance, and dt drops with the joint stiffness); drive their softening
with the ft/Gf and c/phi laws already coded in `Material`; let broken joints
release element faces as new discrete boundaries; and reuse the DEM contact
machinery (node-edge penalty first, potential contact later) plus the
existing fragment bookkeeping. Crack paths then follow mesh edges, so use
unstructured, randomised meshes rather than the current structured
cross-diagonal one. The `Solver` interface, `Config`/`Material`/`Tool`
layer and output formats need no change — only a new `FdemSolver`.

**Reading list.** Munjiza, Owen & Bicanic (1995); Munjiza & Andrews (1998,
NBS contact detection); Munjiza, Andrews & White (1999, combined
single/smeared crack); Munjiza, *The Combined Finite-Discrete Element
Method*, Wiley (2004); Potyondy & Cundall (2004, the BPM that v1's DEM
actually is); Lisjak & Grasselli (2014, FDEM review for rock mechanics).
Reference implementations: Y2D/Y3D, Solidity, Irazu, HOSS, OpenFDEM.
