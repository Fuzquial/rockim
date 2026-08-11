#pragma once
// ---------------------------------------------------------------------------
// Minimal ASCII VTK XML (.vtu) writers, readable directly in ParaView.
// Three flavours:
//   - triangle meshes (FEM): cell scalars (damage...) + point vectors (vel...)
//   - particle clouds (DEM): VTK_VERTEX cells with scalars (radius, fragment...)
//   - bond networks (DEM):   VTK_LINE cells with scalars (state, break time)
// ---------------------------------------------------------------------------
#include <array>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Dense>

namespace rockim::vtk {

using ScalarField  = std::map<std::string, const std::vector<double>*>;
using VectorField  = std::map<std::string, const std::vector<Eigen::Vector2d>*>;
using VectorField3 = std::map<std::string, const std::vector<Eigen::Vector3d>*>;

void writeTriMesh(const std::string& path,
                  const std::vector<Eigen::Vector2d>& points,
                  const std::vector<std::array<int, 3>>& tris,
                  const ScalarField& cellScalars,
                  const VectorField& pointVectors);

void writeParticles(const std::string& path,
                    const std::vector<Eigen::Vector2d>& points,
                    const ScalarField& pointScalars,
                    const VectorField& pointVectors);

void writeLines(const std::string& path,
                const std::vector<Eigen::Vector2d>& points,
                const std::vector<std::array<int, 2>>& lines,
                const ScalarField& cellScalars);

// 3D meshes
void writeTetMesh(const std::string& path,
                  const std::vector<Eigen::Vector3d>& points,
                  const std::vector<std::array<int, 4>>& tets,
                  const ScalarField& cellScalars,
                  const VectorField3& pointVectors);

void writeTriangles3(const std::string& path,
                     const std::vector<Eigen::Vector3d>& points,
                     const std::vector<std::array<int, 3>>& tris,
                     const ScalarField& cellScalars);

// 3D variants (true z coordinate, 3-component vectors)
void writeParticles(const std::string& path,
                    const std::vector<Eigen::Vector3d>& points,
                    const ScalarField& pointScalars,
                    const VectorField3& pointVectors);

void writeLines(const std::string& path,
                const std::vector<Eigen::Vector3d>& points,
                const std::vector<std::array<int, 2>>& lines,
                const ScalarField& cellScalars);

} // namespace rockim::vtk
