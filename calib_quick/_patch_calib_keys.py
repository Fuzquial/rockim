import io, os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
H="include/rockim/FdemSolver.hpp"; C="src/FdemSolver.cpp"; M="src/MatLaw.cpp"
h=io.open(H,encoding="utf-8").read(); c=io.open(C,encoding="utf-8").read(); m=io.open(M,encoding="utf-8").read()
def rep(s,a,b,name):
    assert s.count(a)==1,(name,s.count(a)); return s.replace(a,b,1)
# ---- hpp : membres ----
h=rep(h,"""    bool ucsStop_ = false;
    double ucsStopDelay_ = 0.0;
""","""    bool ucsStop_ = false;
    double ucsStopDelay_ = 0.0;
    // ---- ajouts calibration (2026-09-02, opt-in, DOCUMENTATION 5.17) --------
    // gripsStopAfterPeak : arret apres la chute post-pic en montage a MORS
    // (miroir de ucsStopAfterPeak, qui ne vaut qu en plateaux) ; stopPeakDrop :
    // fraction de chute sous le pic qui verrouille (defaut historique : sigma
    // < 0,3 pic, soit 70 % de chute — inatteignable sous confinement) ;
    // historyStrains : colonnes epsAx, epsLat, epsVol (faces de la boite).
    bool gripsStop_ = false;
    double gripsStopDelay_ = 0.0;
    double lockDrop_ = -1.0;
    bool stopKeysRead_ = false;
    bool histStrains_ = false;
    std::vector<int> hsTop_, hsBot_, hsLeft_, hsRight_;
    void setupHistoryStrains();
    void historyStrains(double& eAx, double& eLat) const;
""","hpp")
# ---- cpp : verrouillage du pic (fraction configurable) ----
c=rep(c,"""        if (!peakLockedU_) {
            sigmaPeak_ = std::max(sigmaPeak_, sigma);
            if (nBroken_ > 0 && sigmaPeak_ > 0.0 && sigma < 0.3 * sigmaPeak_) {
""","""        if (!stopKeysRead_) {              // lecture unique (opt-in, 2026-09-02)
            stopKeysRead_ = true;
            gripsStop_ = cfg_.getb("gripsStopAfterPeak", false);
            gripsStopDelay_ = cfg_.getd("gripsStopDelay", 0.0);
            lockDrop_ = cfg_.getd("stopPeakDrop", -1.0);
            if (lockDrop_ >= 0.0 && !(lockDrop_ > 0.0 && lockDrop_ < 1.0))
                throw std::runtime_error("stopPeakDrop doit etre dans ]0 ; 1[ "
                                         "(fraction de chute sous le pic)");
        }
        // lockFrac = 0,3 (historique : chute de 70 %) ou 1 - stopPeakDrop
        const double lockFrac = (lockDrop_ > 0.0) ? (1.0 - lockDrop_) : 0.3;
        if (!peakLockedU_) {
            sigmaPeak_ = std::max(sigmaPeak_, sigma);
            if (nBroken_ > 0 && sigmaPeak_ > 0.0 && sigma < lockFrac * sigmaPeak_) {
""","lock")
# ---- cpp : finished() ----
c=rep(c,"""    if (ucsStop_ && scen_ == Scenario::TENSION && tensionPlatens_)
        return peakLockedU_ && tLockedU_ >= 0.0
               && t_ >= tLockedU_ + ucsStopDelay_;
""","""    if (gripsStop_ && scen_ == Scenario::TENSION && !tensionPlatens_)
        return peakLockedU_ && tLockedU_ >= 0.0
               && t_ >= tLockedU_ + gripsStopDelay_;
    if (ucsStop_ && scen_ == Scenario::TENSION && tensionPlatens_)
        return peakLockedU_ && tLockedU_ >= 0.0
               && t_ >= tLockedU_ + ucsStopDelay_;
""","finished")
# ---- cpp : appel de setupHistoryStrains ----
c=rep(c,"""    setupStrainGauge();                    // after the platens: needs plTop_.y
""","""    setupStrainGauge();                    // after the platens: needs plTop_.y
    if (cfg_.getb("historyStrains", false)) setupHistoryStrains();
""","init call")
# ---- cpp : fonctions ----
c=rep(c,"""    double L0 = gHiY_ - gLoY_;
    if (L0 > 0.0)
        epsGauge = (meanUy(gLoNodes_) - meanUy(gHiNodes_)) / L0;
}
""","""    double L0 = gHiY_ - gLoY_;
    if (L0 > 0.0)
        epsGauge = (meanUy(gLoNodes_) - meanUy(gHiNodes_)) / L0;
}

// ---------------------------------------------------------------------------
// historyStrains (2026-09-02, opt-in) : deformations moyennes des quatre faces
// de la boite pour l historique — epsAx = (u_y haut - u_y bas)/H, epsLat =
// (u_x droite - u_x gauche)/W, epsVol = epsAx + epsLat (deformation plane,
// eps_zz = 0), convention TRACTION POSITIVE. Sortie seulement. C est ce qui
// permet les seuils sigma_ci / sigma_cd par la methode SBM (inversion de la
// deformation volumique) avec le MEME operateur que sur l essai. Les noeuds
// sont dupliques par element (maillage cohesif) : la moyenne sur toutes les
// copies d une face pondere chaque noeud par ses elements incidents.
// ---------------------------------------------------------------------------
void FdemSolver::setupHistoryStrains() {
    histStrains_ = true;
    hsTop_.clear(); hsBot_.clear(); hsLeft_.clear(); hsRight_.clear();
    for (int i = 0; i < (int)X0_.size(); ++i) {
        const double x = X0_[i].x(), y = X0_[i].y();
        if (y < 1e-9) hsBot_.push_back(i);
        if (y > H_ - 1e-9) hsTop_.push_back(i);
        if (x < 1e-9) hsLeft_.push_back(i);
        if (x > W_ - 1e-9) hsRight_.push_back(i);
    }
    if (hsTop_.empty() || hsBot_.empty() || hsLeft_.empty() || hsRight_.empty())
        throw std::runtime_error("historyStrains : une face de la boite n a aucun "
                                 "noeud (geometrie non rectangulaire ?)");
    std::cout << "[FDEM] historyStrains : faces haut/bas/gauche/droite = "
              << hsTop_.size() << "/" << hsBot_.size() << "/" << hsLeft_.size()
              << "/" << hsRight_.size()
              << " noeuds -> colonnes epsAx, epsLat, epsVol (traction > 0)"
              << std::endl;
}

void FdemSolver::historyStrains(double& eAx, double& eLat) const {
    auto meanU = [&](const std::vector<int>& ns, int comp) {
        double s = 0.0;
        for (int i : ns) s += (comp == 0) ? u_[i].x() : u_[i].y();
        return ns.empty() ? 0.0 : s / (double)ns.size();
    };
    eAx = (meanU(hsTop_, 1) - meanU(hsBot_, 1)) / H_;
    eLat = (meanU(hsRight_, 0) - meanU(hsLeft_, 0)) / W_;
}
""","functions")
# ---- cpp : historique (en-tete + ligne) ----
c=rep(c,"""        if (adaptive_) os << ",nInserted,nDamaging";
        if (bdOn_) os << ",nPulv,bdWork";
        os << "\\n";
        return;
""","""        if (adaptive_) os << ",nInserted,nDamaging";
        if (bdOn_) os << ",nPulv,bdWork";
        if (histStrains_) os << ",epsAx,epsLat,epsVol";
        os << "\\n";
        return;
""","header")
c=rep(c,"""        if (adaptive_) { long ni, nd; countInserted(ni, nd);
                         os << "," << ni << "," << nd; }
        if (bdOn_) os << "," << nPulv_ << "," << bdWork_;
        os << "\\n";
        return;
""","""        if (adaptive_) { long ni, nd; countInserted(ni, nd);
                         os << "," << ni << "," << nd; }
        if (bdOn_) os << "," << nPulv_ << "," << bdWork_;
        if (histStrains_) {
            double ea = 0.0, el = 0.0;
            historyStrains(ea, el);
            os << "," << ea << "," << el << "," << (ea + el);
        }
        os << "\\n";
        return;
""","row")
# ---- cpp : weibullScope = lcz ----
c=rep(c,"""    const bool wGf =
        cfg_.gets("weibullScope", "strength") == "strengthGf";
""","""    const bool wGf =
        cfg_.gets("weibullScope", "strength") == "strengthGf";
    // lcz (2026-09-02) : Gf et GfII suivent stat^2, donc l_cz = E Gf / ft^2 est
    // CONSTANT joint a joint — disperser la resistance sans changer la
    // ductilite (la regle du balayage calib_quick : Gf ~ ft^2)
    const bool wLcz = cfg_.gets("weibullScope", "strength") == "lcz";
""","weibull a")
c=rep(c,"""        if (wGf) { J.Gf *= J.stat; J.GfII *= J.stat; }
""","""        if (wGf) { J.Gf *= J.stat; J.GfII *= J.stat; }
        else if (wLcz) { J.Gf *= J.stat * J.stat; J.GfII *= J.stat * J.stat; }
""","weibull b")
m=rep(m,"""    if (wsc != "strength" && wsc != "strengthGf")
        throw std::runtime_error("weibullScope must be strength | strengthGf");
""","""    if (wsc != "strength" && wsc != "strengthGf" && wsc != "lcz")
        throw std::runtime_error("weibullScope must be strength | strengthGf | lcz");
""","matlaw")
# ---- cpp : avertissement jointPenaltyFactor inerte ----
c=rep(c,"""        std::cout << "[FDEM] adaptive insertion: " << jt_.size()
                  << " bonded edges, activation penalty "
                  << cfg_.getd("insertionPenaltyFactor", 4.0) << " E/h\\n";
""","""        std::cout << "[FDEM] adaptive insertion: " << jt_.size()
                  << " bonded edges, activation penalty "
                  << cfg_.getd("insertionPenaltyFactor", 4.0) << " E/h\\n";
        if (cfg_.has("jointPenaltyFactor"))
            std::cout << "[FDEM] WARNING: jointPenaltyFactor est INERTE en insertion = "
                         "adaptive : la penalite des joints inseres est "
                         "insertionPenaltyFactor (" << cfg_.getd("insertionPenaltyFactor", 4.0)
                      << " E/h)" << std::endl;
""","warning")
io.open(H,"w",encoding="utf-8").write(h); io.open(C,"w",encoding="utf-8").write(c); io.open(M,"w",encoding="utf-8").write(m)
print("cles de calibration : 10 greffes appliquees (hpp 1, FdemSolver.cpp 8, MatLaw.cpp 1)")
