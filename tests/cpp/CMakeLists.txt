add_executable(multiax_tests
  test_main.cpp
  test_core.cpp
  test_woehler.cpp
)
target_include_directories(multiax_tests SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/third_party/doctest)
target_link_libraries(multiax_tests PRIVATE multiax::core)
multiax_copy_runtime(multiax_tests)

# one CTest entry per doctest test case suite; '--force-colors=0' keeps logs readable
add_test(NAME multiax_tests COMMAND multiax_tests --force-colors=0)
