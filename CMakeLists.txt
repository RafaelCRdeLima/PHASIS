cmake_minimum_required(VERSION 3.16)
project(phasis LANGUAGES CXX)

# =====================================================================
# PHASIS -- transporte de neutrinos UHE ao longo de geodesicas nulas.
#
# Sem dependencias externas. Apenas C++17 e a STL.
#
# O subprojeto dipole/ tem build proprio (Makefile) e NAO e construido
# aqui: o acoplamento entre os dois e por ARQUIVO (as tabelas de secao
# de choque), nao por link.
# =====================================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_compile_options(-O2 -Wall -Wextra)

add_library(phasis_core
  src/integrate.cpp
  src/metric.cpp
  src/density.cpp
  src/cross_section.cpp
  src/trace.cpp
)
target_include_directories(phasis_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(trace src/trace_main.cpp)
target_link_libraries(trace PRIVATE phasis_core)

enable_testing()
add_executable(test_phase1 tests/test_phase1.cpp)
target_link_libraries(test_phase1 PRIVATE phasis_core)
add_test(NAME phase1 COMMAND test_phase1)
