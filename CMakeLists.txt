cmake_minimum_required(VERSION 3.16)
project(rockim CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Eigen : le paquet système s'il existe, sinon téléchargement automatique des
# en-têtes (header-only, aucune étape de build).
# CORRECTIF 2026-08-30 — deux casses réparées, mesurées sur macOS ARM :
#  1. `find_package(Eigen3 3.3)` REJETTE Eigen >= 4 : le fichier de version
#     d'Eigen déclare une compatibilité même-version-majeure, donc une machine
#     avec Eigen 5.x (conda-forge le sert par défaut) partait télécharger —
#     alors que le code compile sans warning contre 3.4.0 ET 5.0.1 (vérifié).
#     On accepte donc toute version trouvée, avec un garde-fou >= 3.3 seul.
#  2. `FetchContent_Populate(<nom déclaré>)` est déprécié (CMP0169) et casse
#     sous CMake >= 4 en mode strict ; la forme DIRECTE ci-dessous ne l'est
#     pas et reproduit l'intention d'origine (en-têtes seuls, pas de
#     add_subdirectory du CMake d'Eigen).
find_package(Eigen3 QUIET NO_MODULE)
if(TARGET Eigen3::Eigen AND DEFINED Eigen3_VERSION AND Eigen3_VERSION VERSION_LESS 3.3)
  message(STATUS "rockim: Eigen ${Eigen3_VERSION} trouvé mais < 3.3 — téléchargement des en-têtes 3.4.0")
  set(_rockim_fetch_eigen ON)
elseif(NOT TARGET Eigen3::Eigen)
  message(STATUS "rockim: Eigen introuvable sur le système — téléchargement des en-têtes 3.4.0 "
                 "(réseau requis ; hors ligne : installer Eigen puis relancer cmake, "
                 "ou passer -DEigen3_DIR=<prefix>/share/eigen3/cmake)")
  set(_rockim_fetch_eigen ON)
else()
  message(STATUS "rockim: Eigen ${Eigen3_VERSION} (système)")
endif()
if(_rockim_fetch_eigen)
  include(FetchContent)
  FetchContent_Populate(eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    SOURCE_DIR   "${CMAKE_BINARY_DIR}/_deps/eigen-src"
    BINARY_DIR   "${CMAKE_BINARY_DIR}/_deps/eigen-bin"
    SUBBUILD_DIR "${CMAKE_BINARY_DIR}/_deps/eigen-subbuild")
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_BINARY_DIR}/_deps/eigen-src")
endif()

add_executable(rockim
  src/main.cpp
  src/Config.cpp
  src/FemSolver.cpp
  src/FdemSolver.cpp
  src/Fdem3dSolver.cpp
  src/DemSolver.cpp
  src/Dem3dSolver.cpp
  src/Fem3dSolver.cpp
  src/MatLaw.cpp
  src/Tessellation.cpp
  src/Tessellation3.cpp
  src/VtkWriter.cpp
)

target_include_directories(rockim PRIVATE include)
target_link_libraries(rockim PRIVATE Eigen3::Eigen)

# OpenMP: parallel element/joint/integrate loops of the FDEM (optional —
# the code builds and runs serial without it). Results are deterministic
# for a FIXED thread count (ordered reduction), not across thread counts.
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
  target_link_libraries(rockim PRIVATE OpenMP::OpenMP_CXX)
endif()

if(MSVC)
  # M_PI et amis, et neutraliser les macros min/max de windows.h
  target_compile_definitions(rockim PRIVATE _USE_MATH_DEFINES NOMINMAX)
  target_compile_options(rockim PRIVATE /O2 /W3)
else()
  target_compile_options(rockim PRIVATE -O3 -Wall -Wextra)
endif()
