import io, os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
F="src/FdemSolver.cpp"
f=io.open(F,encoding="utf-8").read()
def rep(s,a,b,name):
    assert s.count(a)==1,(name,s.count(a)); return s.replace(a,b,1)
f=rep(f,"""        std::vector<double> aSum(phases_.n(), 0.0), dSum(phases_.n(), 0.0),
                            d2Sum(phases_.n(), 0.0);
""","""        std::vector<double> aSum(phases_.n(), 0.0), dSum(phases_.n(), 0.0),
                            d2Sum(phases_.n(), 0.0), adSum(phases_.n(), 0.0);
""","decl")
f=rep(f,"""            aSum[p] += a; dSum[p] += dg; d2Sum[p] += dg * dg; ++nG[p]; aTot += a;
""","""            aSum[p] += a; dSum[p] += dg; d2Sum[p] += dg * dg; ++nG[p]; aTot += a;
            adSum[p] += a * dg;                       // moyenne PONDEREE PAR L AIRE
""","acc")
f=rep(f,"""                      << " %), " << nG[p] << " grains, d_eq = " << 1000.0 * m << " +- "
                      << 1000.0 * sd << " mm" << (pSize.size() && pSize[p] > 0.0
""","""                      << " %), " << nG[p] << " grains, d_eq = " << 1000.0 * m << " +- "
                      << 1000.0 * sd << " mm en nombre, "
                      << (aSum[p] > 0.0 ? 1000.0 * adSum[p] / aSum[p] : 0.0)
                      << " mm pondere par l aire" << (pSize.size() && pSize[p] > 0.0
""","print")
io.open(F,"w",encoding="utf-8").write(f)
print("journal : moyenne ponderee par l aire ajoutee")
