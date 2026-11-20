cmake_minimum_required(VERSION 3.24)

project(multiax VERSION 0.1.0 LANGUAGES CXX
        DESCRIPTION "Multiax-Invariant fatigue assessment (Bolchoun 2018) - data-oriented C++ core")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# When built by scikit-build-core (pip install), only the Python module is needed.
if(SKBUILD)
  set(_mx_default_extras OFF)
else()
  set(_mx_default_extras ON)
endif()
option(MULTIAX_BUILD_PYTHON "Build the Python module multiax._core" ON)
option(MULTIAX_BUILD_TESTS "Build the C++ unit tests (doctest)" ${_mx_default_extras})
option(MULTIAX_BUILD_BENCH "Build the benchmarks (nanobench)" ${_mx_default_extras})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(MultiaxOptions)

add_subdirectory(multiax_core)

if(MULTIAX_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests/cpp)
endif()

if(MULTIAX_BUILD_BENCH)
  add_subdirectory(bench)
endif()

if(MULTIAX_BUILD_PYTHON)
  add_subdirectory(python)
endif()
