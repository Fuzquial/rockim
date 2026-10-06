#include "rockim/VtkWriter.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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

// ---------------------------------------------------------------------------
// (2026-10-03, performances) FORMATAGE PARALLELE, OCTETS IDENTIQUES.
// L ecriture ASCII d une trame (159 k noeuds, 40 k tetras) coutait ~2 s de
// formatage de flottants en SERIE, soit l equivalent de 30 a 90 pas de calcul
// par trame. Les lignes d un tableau sont maintenant formatees par tranches
// dans des ostringstream qui recopient TOUS les reglages du flux de sortie
// (copyfmt : precision, drapeaux, largeur, remplissage, locale), puis les
// tranches sont ecrites dans l ordre : le fichier est le meme octet pour octet
// quel que soit le nombre de fils. Les petits tableaux restent en serie.
// ---------------------------------------------------------------------------
template <class Fmt>
static void emitRows(std::ofstream& out, std::size_t n, Fmt&& fmt) {
    constexpr std::size_t kSerialMax = 8192;
    if (n < kSerialMax) {
        for (std::size_t i = 0; i < n; ++i) fmt(static_cast<std::ostream&>(out), i);
        return;
    }
    constexpr int kParts = 64;
    std::vector<std::string> part(kParts);
#pragma omp parallel for schedule(dynamic, 1)
    for (int c = 0; c < kParts; ++c) {
        std::ostringstream os;
        os.copyfmt(out);
        os.exceptions(std::ios::goodbit);
        const std::size_t i0 = n * (std::size_t)c / kParts;
        const std::size_t i1 = n * (std::size_t)(c + 1) / kParts;
        for (std::size_t i = i0; i < i1; ++i) fmt(static_cast<std::ostream&>(os), i);
        part[c] = os.str();
    }
    for (const auto& sp : part) out.write(sp.data(), (std::streamsize)sp.size());
}

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
    emitRows(out, pts.size(), [&](std::ostream& o, std::size_t i) {
        o << pts[i].x() << " " << pts[i].y() << " 0\n";
    });
    out << "</DataArray>\n</Points>\n";
}

static void writeCellScalars(std::ofstream& out, const ScalarField& fields) {
    if (fields.empty()) return;
    out << "<CellData>\n";
    for (const auto& [name, vec] : fields) {
        out << "<DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        const auto& vv = *vec;
        emitRows(out, vv.size(), [&](std::ostream& o, std::size_t i) { o << vv[i] << "\n"; });
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
        const auto& vv = *vec;
        emitRows(out, vv.size(), [&](std::ostream& o, std::size_t i) { o << vv[i] << "\n"; });
        out << "</DataArray>\n";
    }
    for (const auto& [name, vec] : vectors) {
        out << "<DataArray type=\"Float64\" Name=\"" << name
            << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        const auto& vv = *vec;
        emitRows(out, vv.size(), [&](std::ostream& o, std::size_t i) {
            o << vv[i].x() << " " << vv[i].y() << " 0\n";
        });
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
    emitRows(out, tris.size(), [&](std::ostream& o, std::size_t i) {
        o << tris[i][0] << " " << tris[i][1] << " " << tris[i][2] << "\n";
    });
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    emitRows(out, tris.size(), [&](std::ostream& o, std::size_t i) { o << 3 * (i + 1) << "\n"; });
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
    emitRows(out, points.size(), [&](std::ostream& o, std::size_t i) { o << i << "\n"; });
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    emitRows(out, points.size(), [&](std::ostream& o, std::size_t i) { o << i + 1 << "\n"; });
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
    emitRows(out, lines.size(), [&](std::ostream& o, std::size_t i) {
        o << lines[i][0] << " " << lines[i][1] << "\n";
    });
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    emitRows(out, lines.size(), [&](std::ostream& o, std::size_t i) { o << 2 * (i + 1) << "\n"; });
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
    emitRows(out, pts.size(), [&](std::ostream& o, std::size_t i) {
        o << pts[i].x() << " " << pts[i].y() << " " << pts[i].z() << "\n";
    });
    out << "</DataArray>\n</Points>\n";
}

static void writePointData3(std::ofstream& out, const ScalarField& sc, const VectorField3& vec) {
    out << "<PointData>\n";
    for (const auto& [name, v] : sc) {
        out << "<DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        const auto& vv = *v;
        emitRows(out, vv.size(), [&](std::ostream& o, std::size_t i) { o << vv[i] << "\n"; });
        out << "</DataArray>\n";
    }
    for (const auto& [name, v] : vec) {
        out << "<DataArray type=\"Float64\" Name=\"" << name
            << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        const auto& vv = *v;
        emitRows(out, vv.size(), [&](std::ostream& o, std::size_t i) {
            o << vv[i].x() << " " << vv[i].y() << " " << vv[i].z() << "\n";
        });
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
    emitRows(out, points.size(), [&](std::ostream& o, std::size_t i) { o << i << "\n"; });
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    emitRows(out, points.size(), [&](std::ostream& o, std::size_t i) { o << i + 1 << "\n"; });
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
    emitRows(out, lines.size(), [&](std::ostream& o, std::size_t i) {
        o << lines[i][0] << " " << lines[i][1] << "\n";
    });
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    emitRows(out, lines.size(), [&](std::ostream& o, std::size_t i) { o << 2 * (i + 1) << "\n"; });
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
    emitRows(out, tets.size(), [&](std::ostream& o, std::size_t i) {
        o << tets[i][0] << " " << tets[i][1] << " " << tets[i][2] << " " << tets[i][3] << "\n";
    });
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    emitRows(out, tets.size(), [&](std::ostream& o, std::size_t i) { o << 4 * (i + 1) << "\n"; });
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
    emitRows(out, tris.size(), [&](std::ostream& o, std::size_t i) {
        o << tris[i][0] << " " << tris[i][1] << " " << tris[i][2] << "\n";
    });
    out << "</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    emitRows(out, tris.size(), [&](std::ostream& o, std::size_t i) { o << 3 * (i + 1) << "\n"; });
    out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < tris.size(); ++i) out << "5\n";   // VTK_TRIANGLE
    out << "</DataArray>\n</Cells>\n";
    writeCellScalars(out, cellScalars);
    closeVtu(out);
}

} // namespace rockim::vtk
