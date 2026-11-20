# Common compile/link settings for all multiax targets.

add_library(multiax_options INTERFACE)
add_library(multiax::options ALIAS multiax_options)

target_compile_features(multiax_options INTERFACE cxx_std_20)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(multiax_options INTERFACE
    -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion)
  if(MINGW)
    # link libc++/libunwind statically: the only runtime DLL left is libomp.dll
    target_link_options(multiax_options INTERFACE
      -static-libgcc -static-libstdc++ -Wl,-Bstatic -lc++ -lunwind -Wl,-Bdynamic)
  endif()
elseif(MSVC)
  target_compile_options(multiax_options INTERFACE /W4 /permissive-)
endif()

find_package(OpenMP REQUIRED)
target_link_libraries(multiax_options INTERFACE OpenMP::OpenMP_CXX)

# Copy libomp.dll next to a target's output (tests, benchmarks).
function(multiax_copy_runtime target)
  if(MINGW AND EXISTS "${MULTIAX_LIBOMP_DLL}")
    add_custom_command(TARGET ${target} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different "${MULTIAX_LIBOMP_DLL}" "$<TARGET_FILE_DIR:${target}>")
  endif()
endfunction()
