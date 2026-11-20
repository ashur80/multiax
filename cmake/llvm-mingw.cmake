# Toolchain file for llvm-mingw (clang + libc++ + libomp, UCRT), see CLAUDE.md.
#
# The compiler is searched in this order:
#   1. environment variable LLVM_MINGW_ROOT (folder containing bin/clang++.exe)
#   2. clang++ on PATH (winget adds the llvm-mingw bin folder to the user PATH)
#   3. the winget package folder %LOCALAPPDATA%/Microsoft/WinGet/Packages/MartinStorsjo.LLVM-MinGW.UCRT*/llvm-mingw-*

set(_mx_hints "")
if(DEFINED ENV{LLVM_MINGW_ROOT})
  list(APPEND _mx_hints "$ENV{LLVM_MINGW_ROOT}/bin")
endif()
if(DEFINED ENV{LOCALAPPDATA})
  file(TO_CMAKE_PATH "$ENV{LOCALAPPDATA}" _mx_lad)
  file(GLOB _mx_winget "${_mx_lad}/Microsoft/WinGet/Packages/MartinStorsjo.LLVM-MinGW.UCRT*/llvm-mingw-*/bin")
  list(APPEND _mx_hints ${_mx_winget})
endif()

find_program(MX_CLANGXX NAMES clang++ HINTS ${_mx_hints} REQUIRED)
find_program(MX_CLANG NAMES clang HINTS ${_mx_hints} REQUIRED)
set(CMAKE_C_COMPILER "${MX_CLANG}" CACHE FILEPATH "" FORCE)
set(CMAKE_CXX_COMPILER "${MX_CLANGXX}" CACHE FILEPATH "" FORCE)

get_filename_component(_mx_bin "${MX_CLANGXX}" DIRECTORY)
get_filename_component(LLVM_MINGW_ROOT "${_mx_bin}" DIRECTORY CACHE)

# libomp.dll must lie next to the Python extension / test executables (no static libomp available)
set(MULTIAX_LIBOMP_DLL "${LLVM_MINGW_ROOT}/x86_64-w64-mingw32/bin/libomp.dll" CACHE FILEPATH "OpenMP runtime DLL")
