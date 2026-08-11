cmake_minimum_required(VERSION 3.16)
project(rockim CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Eigen: use the system package if present, otherwise fetch the headers
# automatically (Eigen is header-only, so this needs no build step).
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  include(FetchContent)
  FetchContent_Declare(eigen
    URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz)
  FetchContent_GetProperties(eigen)
  if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
  endif()
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${eigen_SOURCE_DIR}")
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
