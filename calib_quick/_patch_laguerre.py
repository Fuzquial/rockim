import io, os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
H="include/rockim/Tessellation.hpp"; C="src/Tessellation.cpp"; F="src/FdemSolver.cpp"
h=io.open(H,encoding="utf-8").read(); c=io.open(C,encoding="utf-8").read(); f=io.open(F,encoding="utf-8").read()
def rep(s,a,b,name):
    assert s.count(a)==1,(name,s.count(a)); return s.replace(a,b,1)
# ---- A1 : signature de voronoiCells + poids ----
c=rep(c,"""std::vector<std::vector<Vec2>> voronoiCells(const std::vector<Vec2>& seeds,
                                            double W, double H) {
    const int N = (int)seeds.size();
""","""// wgt (2026-09-02, optionnel) : diagramme de LAGUERRE (cellules de puissance).
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
std::vector<std::vector<Vec2>> voronoiCells(const std::vector<Vec2>& seeds,
                                            double W, double H,
                                            const std::vector<double>* wgt = nullptr) {
    const int N = (int)seeds.size();
    const bool lag = (wgt != nullptr);
    double wMax = 0.0;
    if (lag) for (double v : *wgt) wMax = std::max(wMax, v);
""","A1")
# ---- A2 : portee par graine ----
c=rep(c,"""    for (int i = 0; i < N; ++i) {
        int ci, cj;
        cellOf(seeds[i], ci, cj);
        std::vector<Vec2> poly;
""","""    for (int i = 0; i < N; ++i) {
        int ci, cj;
        cellOf(seeds[i], ci, cj);
        std::vector<Vec2> poly;
        const double dwi = lag ? std::sqrt(std::max(0.0, wMax - (*wgt)[i])) : 0.0;
""","A2")
# ---- A3 : early-out, ligne radicale, cellule vide ----
c=rep(c,"""                if (0.25 * cand[k].first > r2max) break;   // sorted: none closer
                int j = cand[k].second;
                Vec2 n = seeds[j] - seeds[i];
                double b = n.dot(0.5 * (seeds[i] + seeds[j]));
                poly = clipHalfPlane(poly, n, b);
                if (poly.size() < 3)
                    throw std::runtime_error("Tessellation: seed cell vanished "
                                             "(coincident seeds?)");
""","""                if (!lag) {
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
""","A3")
# ---- A4 : couverture en anneaux ----
c=rep(c,"""            double rcov = ring * csMin;
            if (coveredAll || 4.0 * r2max <= rcov * rcov) break;
""","""            if (lag && poly.empty()) break;
            double rcov = ring * csMin;
            if (lag) {
                if (coveredAll || 2.0 * std::sqrt(r2max) + dwi <= rcov) break;
            } else if (coveredAll || 4.0 * r2max <= rcov * rcov) break;
""","A4")
# ---- A5 : portee de sp ----
c=rep(c,"""    std::vector<Vec2> seeds;
    double eps = 1e-6 * std::min(W, H);
    if (randomSeeds) {
""","""    std::vector<Vec2> seeds;
    std::vector<double> sp;      // espacement propre a chaque graine (polydispersite)
    double eps = 1e-6 * std::min(W, H);
    if (randomSeeds) {
""","A5a")
c=rep(c,"""        std::vector<double> sp;
        std::normal_distribution<double> N01(0.0, 1.0);
""","""        std::normal_distribution<double> N01(0.0, 1.0);
""","A5b")
# ---- A6 : nombre de graines vise ----
c=rep(c,"""        long target = std::lround(W * H / (std::sqrt(3.0) / 2.0 * s * s));
""","""        long target = std::lround(W * H / (std::sqrt(3.0) / 2.0 * s * s));
        if (sizeSpread > 0.0)      // aire moyenne x E[L^2] = exp(sigma^2)
            target = std::lround(W * H / (std::sqrt(3.0) / 2.0 * s * s
                                          * std::exp(sizeSpread * sizeSpread)));
""","A6")
# ---- A7 : poids de Laguerre + retrait des graines redondantes ----
c=rep(c,"""    // ---- 2. Voronoi cells (+ Lloyd relaxation) -------------------------------
    std::vector<std::vector<Vec2>> cells;
    for (int it = 0;; ++it) {
        cells = voronoiCells(seeds, W, H);
        if (it >= lloyd) break;
""","""    // ---- 2. Voronoi cells (+ Lloyd relaxation) -------------------------------
    // Polydispersite : poids de Laguerre w_i = (kappa s_i)^2. Avec la bande
    // d exclusion d_ij >= 0,35 (s_i + s_j), kappa <= 0,35 garantit
    // |w_i - w_j| < d_ij^2, donc la ligne radicale passe ENTRE les deux graines
    // (pas de redondance par paire) ; kappa = 0,315 laisse une marge. La
    // relaxation de Lloyd deplace les graines vers les centroides de leurs
    // cellules de puissance, les poids restent : les tailles sont conservees.
    std::vector<double> wgt;
    if (sizeSpread > 0.0) {
        const double kappa = 0.315;
        wgt.resize(seeds.size());
        for (std::size_t i = 0; i < seeds.size(); ++i)
            wgt[i] = (kappa * sp[i]) * (kappa * sp[i]);
    }
    int nRedundant = 0;
    std::vector<std::vector<Vec2>> cells;
    for (int it = 0;; ++it) {
        cells = voronoiCells(seeds, W, H, wgt.empty() ? nullptr : &wgt);
        if (!wgt.empty()) {
            std::vector<int> keep;
            for (std::size_t i = 0; i < cells.size(); ++i)
                if (cells[i].size() >= 3) keep.push_back((int)i);
            if (keep.size() < cells.size()) {
                std::vector<Vec2> s2; std::vector<double> w2, p2;
                for (int i : keep) { s2.push_back(seeds[i]); w2.push_back(wgt[i]); p2.push_back(sp[i]); }
                nRedundant += (int)(cells.size() - keep.size());
                seeds.swap(s2); wgt.swap(w2); sp.swap(p2);
                --it;                  // recalcul sans consommer une iteration
                continue;
            }
        }
        if (it >= lloyd) break;
""","A7")
c=rep(c,"""    Tessellation T;
    T.nGrains = (int)cells.size();
""","""    if (nRedundant > 0)
        std::printf("[tess] laguerre: %d graine(s) redondante(s) retiree(s) "
                    "(cellule de puissance vide)\\n", nRedundant);

    Tessellation T;
    T.nGrains = (int)cells.size();
""","A8")
# ---- hpp : commentaire ----
h=rep(h,"""    // sizeSpread : POLYDISPERSITE (2026-09-02, opt-in). Ecart-type de
    //            ln(taille) : chaque graine de Poisson porte son propre
    //            espacement s_i = s L_i, L_i log-normale de moyenne 1, et deux
    //            graines s acceptent a distance >= 0,35 (s_i + s_j). Moyenne
    //            conservee, cellules petites et grandes melangees. Exige
    //            randomSeeds. 0 = comportement historique, bit-identique.
""","""    // sizeSpread : POLYDISPERSITE (2026-09-02, opt-in). Ecart-type de
    //            ln(taille) : chaque graine de Poisson porte son propre
    //            espacement s_i = s L_i, L_i log-normale de moyenne 1, deux
    //            graines s acceptent a distance >= 0,35 (s_i + s_j), et les
    //            cellules sont celles du diagramme de LAGUERRE de poids
    //            (0,315 s_i)^2 : la graine lourde recoit la grande cellule (un
    //            Voronoi ordinaire moyennerait les tailles des voisines et
    //            ecraserait la dispersion : 0,5 demande -> 0,16 realise).
    //            Moyenne conservee, nombre de grains / exp(sigma^2). Exige
    //            randomSeeds. 0 = comportement historique, bit-identique.
""","hpp")
# ---- FdemSolver : journal ----
f=rep(f,"""                     "ln(d_eq) REALISE = " << sdLn << " (demande " << gSpread
                  << ", Lloyd " << lloyd << " iterations homogeneise)" << std::endl;
""","""                     "ln(d_eq) REALISE = " << sdLn << " (demande " << gSpread
                  << ", diagramme de Laguerre, Lloyd " << lloyd << ")" << std::endl;
""","fdem 1")
f=rep(f,"""        if (gSpread > 0.0 && lloyd > 1)
            std::cout << "[FDEM] WARNING: lloydIters = " << lloyd << " avec grainSizeSpread : "
                         "chaque iteration de Lloyd rapproche les cellules d une taille "
                         "unique — comparer l ecart-type REALISE a la demande, et "
                         "reduire lloydIters si l ecart est trop grand" << std::endl;
""","""        if (gSpread > 0.0 && sdLn < 0.6 * gSpread)
            std::cout << "[FDEM] WARNING: ecart-type de ln(d_eq) realise (" << sdLn
                      << ") < 60 % de la demande (" << gSpread << ") : domaine trop "
                         "petit pour la queue de la distribution, ou contraction des "
                         "aretes (vertexMergeFrac) trop forte" << std::endl;
""","fdem 2")
io.open(H,"w",encoding="utf-8").write(h); io.open(C,"w",encoding="utf-8").write(c); io.open(F,"w",encoding="utf-8").write(f)
print("laguerre : 11 greffes appliquees (Tessellation.cpp 9, hpp 1, FdemSolver.cpp 2)")
