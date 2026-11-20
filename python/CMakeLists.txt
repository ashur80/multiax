# Python extension multiax._core (pybind11). Built by scikit-build-core via `pip install -e .`
# or as part of a normal CMake build (then importable from build/<preset>/python).

find_package(Python 3.10 COMPONENTS Interpreter Development.Module REQUIRED)
if(NOT pybind11_DIR)
  execute_process(COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
                  OUTPUT_VARIABLE pybind11_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core MODULE src/module.cpp)
target_link_libraries(_core PRIVATE multiax::core)
multiax_copy_runtime(_core)

if(SKBUILD)
  install(TARGETS _core LIBRARY DESTINATION multiax RUNTIME DESTINATION multiax)
  if(MINGW AND EXISTS "${MULTIAX_LIBOMP_DLL}")
    install(FILES "${MULTIAX_LIBOMP_DLL}" DESTINATION multiax)
  endif()
endif()
