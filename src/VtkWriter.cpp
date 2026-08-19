#include "rockim/VtkWriter.hpp"
#include <fstream>
#include <iomanip>
#include <stdexcept>

// ---------------------------------------------------------------------------
// PRECISION D'ECRITURE DES COORDONNEES. std::ofstream applique par defaut 6
// chiffres significatifs. Sur un domaine metrique c'est un PAS DE LECTURE de
// 10 micrometres, et tout deplacement plus petit est ecrase a zero.
//
// Mesure du 2026-08-18, benchmark Parker (AbuAisha et al. 2017, annexe A) :
// bloc de 8 m, fissure sous pression dont l'ouverture theorique vaut 0,128 mm,
// soit TREIZE pas de quantification. Les deplacements du VTU ne prenaient plus
// que les valeurs 1e-5 et 1,414e-5 m, les levres lisaient exactement zero, et
// l'ecart de -6,25 % annonce sur le maillage grossier tenait entierement dans
// un cran de graduation — indecidable entre 0 et -14 %.
//
// C'est un defaut de SORTIE, pas de calcul : le solveur travaille en double.
// Il mord des que le rapport taille du domaine / effet mesure depasse ~1e5 —
// donc sur toute etude de petits deplacements sur grand domaine. L'etude
// tunnel y echappait (on y mesurait des metres d'EDZ), la coupe PDC aussi
// (bloc de 40 mm). Meme famille que le defaut toolY, qui rendait la courbe
// force-penetration en escalier.
//
// 12 chiffres portent la resolution a ~1e-9 m sur un domaine metrique, soit
// cinq ordres de grandeur sous l'effet cherche. Le fichier grossit d'environ
// un tiers. AUCUN effet sur la physique.
// ---------------------------------------------------------------------------
static constexpr int kCoordDigits = 12;

namespace rockim::vtk {

static void openVtu(std::ofstream& out, const std::string& path,
                    std::size_t nPts, std::size_t nCells) {
    out.open(path);
    if (!out) throw std::runtime_error("VtkWriter: cannot open '" + path + "'");
    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "<UnstructuredGrid>\n"
        << "<Piece NumberOfPoints=\"" << nPts << "\" NumberOfCells=\"" << nCells << "\">\n";
}

static void writePoints(std::ofstream& out, const std::vector<Eigen::Vector2d>& pts) {
    out << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    out << std::setprecision(kCoordDigits);
    for (const auto& p : pts) out << p.x() << " " << p.y() << " 0\n";
    out << "</DataArray>\n</Points>\n";
}

static void writeCellScalars(std::ofstream& out, const ScalarField& fields) {
    if (fields.empty()) return;
    out << "<CellData>\n";
    for (const auto& [name, vec] : fields) {
        out << "<DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        for (double v : *vec) out << v << "\n";
        out << "</DataArray>\n";
    }
    out << "</CellData>\n";
}

static void writePointData(std::ofstream& out, const ScalarField& scalars,
                           const VectorField& vectors) {
    if (scalars.empty() && vectors.empty()) return;
    out << "<PointData>\n";
    for (const auto& [name, vec] : scalars) {
        out << "<DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        for (double v : *vec) out << v << "\n";
        out << "</DataArray>\n";
    }
    for (const auto& [name, vec] : vectors) {
        out << "<DataArray type=\"Float64\" Name=\"" << name
            << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (const auto& v : *vec) out << v.x() << " " << v.y() << " 0\n";
        out << "</DataArray>\n";
    }
    out << "</PointData>\n";
}

static void closeVtu(std::ofstream& out) {
    out << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

void writeTriMesh(const std::string& path,
                  const std::vector<Eigen::Vector2d>& points,
                  const std::vector<std::array<int, 3>>& tris,
                  const ScalarField& cellScalars,
                  const VectorField& pointVectors) {
    std::ofstream out;
    openVtu(out, path, points.size(), tris.size());
    writePoints(out, points);
    out << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& t : tris) out << t[0] << " " << t[1] << " " << t[2] << "\n";
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 1; i <= tris.size(); ++i) out << 3 * i << "\n";
    out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < tris.size(); ++i) out << "5\n";  // VTK_TRIANGLE
    out << "</DataArray>\n</Cells>\n";
    writeCellScalars(out, cellScalars);
    writePointData(out, {}, pointVectors);
    closeVtu(out);
}

void writeParticles(const std::string& path,
                    const std::vector<Eigen::Vector2d>& points,
                    const ScalarField& pointScalars,
                    const VectorField& pointVectors) {
    std::ofstream out;
    openVtu(out, path, points.size(), points.size());
    writePoints(out, points);
    out << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < points.size(); ++i) out << i << "\n";
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 1; i <= points.size(); ++i) out << i << "\n";
    out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < points.size(); ++i) out << "1\n";  // VTK_VERTEX
    out << "</DataArray>\n</Cells>\n";
    writePointData(out, pointScalars, pointVectors);
    closeVtu(out);
}

void writeLines(const std::string& path,
                const std::vector<Eigen::Vector2d>& points,
                const std::vector<std::array<int, 2>>& lines,
                const ScalarField& cellScalars) {
    std::ofstream out;
    openVtu(out, path, points.size(), lines.size());
    writePoints(out, points);
    out << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& l : lines) out << l[0] << " " << l[1] << "\n";
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 1; i <= lines.size(); ++i) out << 2 * i << "\n";
    out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < lines.size(); ++i) out << "3\n";  // VTK_LINE
    out << "</DataArray>\n</Cells>\n";
    writeCellScalars(out, cellScalars);
    closeVtu(out);
}

// ---------------------------------------------------------------------------
// 3D variants
// ---------------------------------------------------------------------------
static void writePoints3(std::ofstream& out, const std::vector<Eigen::Vector3d>& pts) {
    out << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    out << std::setprecision(kCoordDigits);
    for (const auto& p : pts) out << p.x() << " " << p.y() << " " << p.z() << "\n";
    out << "</DataArray>\n</Points>\n";
}

static void writePointData3(std::ofstream& out, const ScalarField& sc, const VectorField3& vec) {
    out << "<PointData>\n";
    for (const auto& [name, v] : sc) {
        out << "<DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        for (double x : *v) out << x << "\n";
        out << "</DataArray>\n";
    }
    for (const auto& [name, v] : vec) {
        out << "<DataArray type=\"Float64\" Name=\"" << name
            << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (const auto& q : *v) out << q.x() << " " << q.y() << " " << q.z() << "\n";
        out << "</DataArray>\n";
    }
    out << "</PointData>\n";
}

void writeParticles(const std::string& path,
                    const std::vector<Eigen::Vector3d>& points,
                    const ScalarField& pointScalars,
                    const VectorField3& pointVectors) {
    std::ofstream out;
    openVtu(out, path, points.size(), points.size());
    writePoints3(out, points);
    out << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < points.size(); ++i) out << i << "\n";
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 1; i <= points.size(); ++i) out << i << "\n";
    out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < points.size(); ++i) out << "1\n";  // VTK_VERTEX
    out << "</DataArray>\n</Cells>\n";
    writePointData3(out, pointScalars, pointVectors);
    closeVtu(out);
}

void writeLines(const std::string& path,
                const std::vector<Eigen::Vector3d>& points,
                const std::vector<std::array<int, 2>>& lines,
                const ScalarField& cellScalars) {
    std::ofstream out;
    openVtu(out, path, points.size(), lines.size());
    writePoints3(out, points);
    out << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& l : lines) out << l[0] << " " << l[1] << "\n";
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 1; i <= lines.size(); ++i) out << 2 * i << "\n";
    out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < lines.size(); ++i) out << "3\n";  // VTK_LINE
    out << "</DataArray>\n</Cells>\n";
    writeCellScalars(out, cellScalars);
    closeVtu(out);
}

void writeTetMesh(const std::string& path,
                  const std::vector<Eigen::Vector3d>& points,
                  const std::vector<std::array<int, 4>>& tets,
                  const ScalarField& cellScalars,
                  const VectorField3& pointVectors) {
    std::ofstream out;
    openVtu(out, path, points.size(), tets.size());
    writePoints3(out, points);
    out << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& c : tets) out << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << "\n";
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 1; i <= tets.size(); ++i) out << 4 * i << "\n";
    out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < tets.size(); ++i) out << "10\n";  // VTK_TETRA
    out << "</DataArray>\n</Cells>\n";
    writeCellScalars(out, cellScalars);
    writePointData3(out, {}, pointVectors);
    closeVtu(out);
}

void writeTriangles3(const std::string& path,
                     const std::vector<Eigen::Vector3d>& points,
                     const std::vector<std::array<int, 3>>& tris,
                     const ScalarField& cellScalars) {
    std::ofstream out;
    openVtu(out, path, points.size(), tris.size());
    writePoints3(out, points);
    out << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& c : tris) out << c[0] << " " << c[1] << " " << c[2] << "\n";
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 1; i <= tris.size(); ++i) out << 3 * i << "\n";
    out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < tris.size(); ++i) out << "5\n";   // VTK_TRIANGLE
    out << "</DataArray>\n</Cells>\n";
    writeCellScalars(out, cellScalars);
    closeVtu(out);
}

} // namespace rockim::vtk
