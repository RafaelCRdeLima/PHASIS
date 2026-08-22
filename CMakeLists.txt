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

find_package(OpenMP)


add_library(phasis_core
  src/shell.cpp
  src/cascade.cpp
  src/ode.cpp
  src/sweep.cpp
  src/integrate.cpp
  src/metric.cpp
  src/density.cpp
  src/cross_section.cpp
  src/trace.cpp
)
target_include_directories(phasis_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(OpenMP_CXX_FOUND)
  target_link_libraries(phasis_core PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(trace src/trace_main.cpp)
target_link_libraries(trace PRIVATE phasis_core)

enable_testing()
add_executable(test_phase1 tests/test_phase1.cpp)
target_link_libraries(test_phase1 PRIVATE phasis_core)
add_test(NAME phase1 COMMAND test_phase1)

add_executable(test_phase2 tests/test_phase2.cpp)
target_link_libraries(test_phase2 PRIVATE phasis_core)
add_test(NAME phase2 COMMAND test_phase2)

add_executable(test_phase3 tests/test_phase3.cpp)
target_link_libraries(test_phase3 PRIVATE phasis_core)
add_test(NAME phase3 COMMAND test_phase3)

add_executable(test_phase4 tests/test_phase4.cpp)
target_link_libraries(test_phase4 PRIVATE phasis_core)
add_test(NAME phase4 COMMAND test_phase4)

add_executable(test_phase5 tests/test_phase5.cpp)
target_link_libraries(test_phase5 PRIVATE phasis_core)
add_test(NAME phase5 COMMAND test_phase5)
