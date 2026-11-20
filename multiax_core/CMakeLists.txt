add_library(multiax_core STATIC
  src/build_info.cpp
)
add_library(multiax::core ALIAS multiax_core)

target_include_directories(multiax_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_compile_definitions(multiax_core PRIVATE MULTIAX_VERSION="${PROJECT_VERSION}")
target_link_libraries(multiax_core PUBLIC multiax::options)
set_target_properties(multiax_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
