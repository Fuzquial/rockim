import io, os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
C="src/FdemSolver.cpp"
c=io.open(C,encoding="utf-8").read()
a = """    if (mesh == "file") {
        if (shpb_ || disc_)
            throw std::runtime_error("mesh = file is implemented for the BOX "
                "geometry (percussion / shear / tension); the disc and shpb "
                "assemblies build their own meshes");
"""
b = """    if (mesh == "file") {
        // 2026-09-02 : le DISQUE accepte aussi un maillage fichier (bresilien
        // sur le meme maillage Delaunay-Gmsh que l eprouvette triaxiale, pour
        // que le BTS soit un observable de calibration coherent). Le fichier
        // EST le disque (meplats compris) : aucune decoupe, discR_/discC_
        // viennent des cles W/H lues avant (diametre = W). shpb reste exclu.
        if (shpb_)
            throw std::runtime_error("mesh = file is implemented for the BOX "
                "and DISC geometries; the shpb assembly builds its own mesh");
        if (disc_)
            std::cout << "[FDEM] geometry = disc + mesh = file : le fichier est "
                         "pris tel quel comme disque (diametre " << 2.0 * discR_
                      << " m, centre " << discC_.x() << ", " << discC_.y() << ")"
                      << std::endl;
"""
assert c.count(a) == 1
io.open(C,"w",encoding="utf-8").write(c.replace(a, b, 1))
print("garde mesh = file + disc relachee")
