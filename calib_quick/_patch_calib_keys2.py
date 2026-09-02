import io, os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
C="src/FdemSolver.cpp"
c=io.open(C,encoding="utf-8").read()
a = """    const bool wLcz = cfg_.gets("weibullScope", "strength") == "lcz";
"""
b = """    const bool wLcz = cfg_.gets("weibullScope", "strength") == "lcz";
    // validation ici aussi : MatLaw::make ne la fait que sous `law`, une faute
    // de frappe passait donc en silence avec un bulk elastique
    {
        const std::string wsc = cfg_.gets("weibullScope", "strength");
        if (cfg_.getd("jointWeibullM", 0.0) > 0.0 && wsc != "strength"
            && wsc != "strengthGf" && wsc != "lcz")
            throw std::runtime_error("weibullScope must be strength | strengthGf | lcz");
    }
"""
assert c.count(a) == 1
c = c.replace(a, b, 1)
io.open(C,"w",encoding="utf-8").write(c)
print("validation weibullScope ajoutee dans applyJointStatistics")
