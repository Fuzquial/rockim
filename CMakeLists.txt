cmake_minimum_required(VERSION 3.16)
project(rockim CXX)

# ---------------------------------------------------------------------------
# rockim_g0 — configuration CMake de NAISSANCE (2026-09-05, mesure D1a du
# PLAN_ROBUSTESSE). Decision (9) : g0 nait avec CMake + Ninja, pas avec
# build_f2.cmd.
#
# CONTRAINTE : les drapeaux MSVC doivent etre STRICTEMENT ceux de
# `build_f2.cmd` de rockim_f2, faute de quoi l'exe ne reproduit pas l'ancre
# de bit-identite `tools/bitid_refs.json` (prise sur rockim_f2w18.exe) :
#
#   cl /nologo /std:c++17 /O2 /EHsc /openmp /D_USE_MATH_DEFINES /DNOMINMAX
#      /I include /I ..\rockim\eigen-3.4.0 src\*.cpp /Fe:rockim.exe
#
# CMake ajoute par defaut, sous MSVC, des drapeaux que `cl` seul n'a pas :
#   * CMAKE_CXX_FLAGS          = /DWIN32 /D_WINDOWS /W3 /GR /EHsc
#   * CMAKE_CXX_FLAGS_RELEASE  = /O2 /Ob2 /DNDEBUG   <-- /DNDEBUG desarme les
#                                                        assert() du code
#   * bibliotheque d'execution = /MD (CRT dynamique) alors que `cl` sans
#     drapeau prend /MT (CRT statique — verifie sur rockim_f2w21.exe, qui
#     n'importe que KERNEL32.dll et VCOMP140.DLL)
# Les trois sont donc remis a la main ci-dessous. /W4 /WX viendront au lot
# D1b, separement, et chaque lot y sera verifie bit-identique par B6.
# ---------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

if(MSVC)
  # Drapeaux exactement equivalents a build_f2.cmd (rien de plus).
  set(CMAKE_CXX_FLAGS                "/nologo /EHsc" CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE        "/O2"           CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/O2 /Zi"       CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG          "/Od /Zi"       CACHE STRING "" FORCE)
  # CRT statique, comme `cl` sans drapeau (CMP0091 = NEW des 3.15).
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

# ---------------------------------------------------------------------------
# Eigen (en-tetes seules). Priorite a la copie locale ..\rockim\eigen-3.4.0,
# celle que build_f2.cmd utilise : c'est la seule qui garantit la reprise de
# l'ancre. A defaut, paquet systeme, puis telechargement (CI Linux, macOS).
# ---------------------------------------------------------------------------
set(ROCKIM_EIGEN_LOCAL "${CMAKE_CURRENT_SOURCE_DIR}/../rockim/eigen-3.4.0"
    CACHE PATH "Copie locale des en-tetes Eigen (celle de build_f2.cmd)")
if(EXISTS "${ROCKIM_EIGEN_LOCAL}/Eigen/Dense")
  message(STATUS "rockim: Eigen 3.4.0 (copie locale ${ROCKIM_EIGEN_LOCAL})")
  add_library(rockim_eigen INTERFACE)
  target_include_directories(rockim_eigen INTERFACE "${ROCKIM_EIGEN_LOCAL}")
else()
  find_package(Eigen3 QUIET NO_MODULE)
  if(TARGET Eigen3::Eigen AND NOT (DEFINED Eigen3_VERSION AND Eigen3_VERSION VERSION_LESS 3.3))
    message(STATUS "rockim: Eigen ${Eigen3_VERSION} (systeme)")
    add_library(rockim_eigen INTERFACE)
    target_link_libraries(rockim_eigen INTERFACE Eigen3::Eigen)
  else()
    message(STATUS "rockim: Eigen introuvable — telechargement des en-tetes 3.4.0 (reseau requis)")
    include(FetchContent)
    FetchContent_Populate(eigen
      URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE
      SOURCE_DIR   "${CMAKE_BINARY_DIR}/_deps/eigen-src"
      BINARY_DIR   "${CMAKE_BINARY_DIR}/_deps/eigen-bin"
      SUBBUILD_DIR "${CMAKE_BINARY_DIR}/_deps/eigen-subbuild")
    add_library(rockim_eigen INTERFACE)
    target_include_directories(rockim_eigen INTERFACE "${CMAKE_BINARY_DIR}/_deps/eigen-src")
  endif()
endif()

# ---------------------------------------------------------------------------
# Cible. L'ordre des sources reprend celui que `cl src\*.cpp` recoit sur ce
# disque (ordre NTFS, lisible dans build_w21.log) : il ne change pas le
# resultat numerique, il rend seulement les deux chaines comparables.
# ---------------------------------------------------------------------------
add_executable(rockim
  src/Config.cpp
  src/Dem3dSolver.cpp
  src/DemSolver.cpp
  src/Fdem3dSolver.cpp
  src/FdemSolver.cpp
  src/Fem3dSolver.cpp
  src/FemSolver.cpp
  src/main.cpp
  src/MatLaw.cpp
  src/Tessellation.cpp
  src/Tessellation3.cpp
  src/ThermoBench.cpp
  src/ToolPdc3d.cpp
  src/VtkWriter.cpp
)

target_include_directories(rockim PRIVATE include)
target_link_libraries(rockim PRIVATE rockim_eigen)

# OpenMP : boucles elements/joints/integration du FDEM (facultatif — le code
# compile et tourne en serie sans). Les resultats sont deterministes a nombre
# de fils FIXE (reduction ordonnee), pas d'un nombre de fils a l'autre : c'est
# pourquoi tools/bitid.py impose OMP_NUM_THREADS.
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
  target_link_libraries(rockim PRIVATE OpenMP::OpenMP_CXX)
else()
  message(WARNING "rockim: OpenMP introuvable — build SERIE, les hachages bitid differeront")
endif()

if(MSVC)
  # M_PI et amis ; neutraliser les macros min/max de windows.h.
  target_compile_definitions(rockim PRIVATE _USE_MATH_DEFINES NOMINMAX)
else()
  target_compile_options(rockim PRIVATE -O3 -Wall -Wextra)
endif()
