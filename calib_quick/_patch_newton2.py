import io, os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
C="src/Tessellation.cpp"
c=io.open(C,encoding="utf-8").read()
def rep(s,a,b,name):
    assert s.count(a)==1,(name,s.count(a)); return s.replace(a,b,1)
# demarrage : Voronoi (w = 0) si une cellule est vide avec les poids courants (KMT exige un depart sans cellule vide)
c=rep(c,"""        const std::size_t N = seeds.size();
        std::vector<std::vector<Vec2>> cl = voronoiCells(seeds, W, H, &wgt, &lap);
        std::vector<double> A(N);
        for (std::size_t i = 0; i < N; ++i) A[i] = cl[i].size() >= 3 ? polyArea(cl[i]) : 0.0;
""","""        const std::size_t N = seeds.size();
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
""","warm start")
c=rep(c,"""            std::fill(delta.begin(), delta.end(), 0.0);
            r = g;
            p = r;
            double rr = dot(r, r);
            const double gg = rr;
            for (int cg = 0; cg < 1000 && rr > 1e-24 * gg; ++cg) {
                matvec(p, q);
                const double alpha = rr / dot(p, q);
""","""            std::fill(delta.begin(), delta.end(), 0.0);
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
""","cg guard")
c=rep(c,"""            wgt.swap(w2);
            A.swap(A2);
            cl.swap(cl2);
        }
    };
""","""            if (dbg) {
                double dmax = 0.0, amin = 1e300;
                std::size_t nEdges = 0;
                for (std::size_t i = 0; i < N; ++i) { dmax = std::max(dmax, std::abs(delta[i])); amin = std::min(amin, A2[i]); nEdges += lap[i].size(); }
                std::printf("[tess-dbg] newton %d : err %.3f  cg %d  aretes %zu  |delta|max %.3e  tau %.4f  amin/eps0 %.3f\\n",
                            k, err, cgIt, nEdges, dmax, tau, amin / eps0);
            }
            wgt.swap(w2);
            A.swap(A2);
            cl.swap(cl2);
        }
    };
""","dbg")
if "#include <cstdlib>" not in c:
    c = c.replace("#include <cstdio>\n", "#include <cstdio>\n#include <cstdlib>\n", 1)
io.open(C,"w",encoding="utf-8").write(c)
print("warm start + garde CG + debug appliques")
