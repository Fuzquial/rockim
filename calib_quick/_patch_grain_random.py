import io, os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
H="include/rockim/Tessellation.hpp"; C="src/Tessellation.cpp"; F="src/FdemSolver.cpp"
h=io.open(H,encoding="utf-8").read(); c=io.open(C,encoding="utf-8").read(); f=io.open(F,encoding="utf-8").read()
def rep(s,a,b,name):
    assert s.count(a)==1,(name,s.count(a)); return s.replace(a,b,1)
# ---- hpp : signature ----
h=rep(h,"""                              double sizeSpread = 0.0,
                              const std::vector<double>& phaseSize = {});
""","""                              double sizeSpread = 0.0,
                              const std::vector<double>& phaseSize = {},
                              bool randomInterior = false);
    // randomInterior : POINTS INTERIEURS ALEATOIRES du Delaunay intra-grain
    //            (2026-09-02, opt-in grainMeshRandom). Le reseau triangulaire
    //            historique pave chaque grain en triangles quasi equilateraux
    //            alignes — le meme defaut que l algorithme frontal de Gmsh,
    //            banni pour l eprouvette (trois directions de fissure imposees
    //            a l interieur des grains). Poisson-disc dans le polygone :
    //            distance >= 0,75 h entre points, marge 0,55 h aux aretes.
    //            false = comportement historique, bit-identique.
""","hpp")
# ---- cpp : signature de la definition ----
c=rep(c,"""                                 double sizeSpread,
                                 const std::vector<double>& phaseSize) {""",
"""                                 double sizeSpread,
                                 const std::vector<double>& phaseSize,
                                 bool randomInterior) {""","cpp sig")
# ---- cpp : points interieurs ----
c=rep(c,"""            // interior points on a triangular lattice, kept one margin away
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
""","""            // interior points on a triangular lattice, kept one margin away
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
""","cpp interior")
c=rep(c,"""        std::printf("[tess] delaunay grain mesh: target element %.4g m "
                    "(%.2f x grain), %d cells, %ld fell back to fan\\n",
                    h, h / targetD, (int)cells.size(), nFallback);
""","""        std::printf("[tess] delaunay grain mesh: target element %.4g m "
                    "(%.2f x grain), %d cells, %ld fell back to fan%s\\n",
                    h, h / targetD, (int)cells.size(), nFallback,
                    randomInterior ? " — points interieurs ALEATOIRES (Poisson-disc)" : "");
""","cpp printf")
# ---- FdemSolver : cle + passage ----
f=rep(f,"""    Tessellation T = Tessellation::build(W_, H_, d, jit, lloyd, mf, refine,
                                         phases_.fraction, rng,
                                         seeding == "random",
                                         gm == "delaunay", gh,
                                         gSpread, pSize);
""","""    // grainMeshRandom (2026-09-02, opt-in) : points interieurs aleatoires du
    // Delaunay intra-grain (le reseau triangulaire par defaut est structure).
    const bool gRandom = cfg_.getb("grainMeshRandom", false);
    if (gRandom && gm != "delaunay")
        throw std::runtime_error("grainMeshRandom exige grainMesh = delaunay");
    Tessellation T = Tessellation::build(W_, H_, d, jit, lloyd, mf, refine,
                                         phases_.fraction, rng,
                                         seeding == "random",
                                         gm == "delaunay", gh,
                                         gSpread, pSize, gRandom);
""","fdem")
io.open(H,"w",encoding="utf-8").write(h); io.open(C,"w",encoding="utf-8").write(c); io.open(F,"w",encoding="utf-8").write(f)
print("grainMeshRandom : 5 greffes appliquees (hpp 1, Tessellation.cpp 3, FdemSolver.cpp 1)")
