add_executable(multiax_bench
  bench_main.cpp
)
target_include_directories(multiax_bench SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/third_party/nanobench)
target_link_libraries(multiax_bench PRIVATE multiax::core)
multiax_copy_runtime(multiax_bench)
