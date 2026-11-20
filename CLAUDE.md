# CLAUDE.md

## Project

`multiax` is a C++ library for fatigue life assessment under multiaxial loading
with constant and variable amplitudes. It implements the method **Multiax-Invariant**.

- `multiax_core/` – static core library: headers in `include/multiax/<module>/`, sources in `src/`.
- `python/` – pybind11 bindings (`src/module.cpp` → `multiax._core`) and the Python package `python/multiax/`.
- `tests/cpp/` (doctest), `tests/python/` (pytest), `bench/` (nanobench), `third_party/` (doctest, nanobench; vendored).
- `glm/` – third-party math library (vendored, not used in milestone I). Do not modify.
- `docs/learning/` – guided learning path (exercises E1–E6). **The exercises are implemented by the author**: don't write solutions for `TODO(E<n>)` code unless asked; provide scaffolding, tests, hints and reviews.
- `docs/` – method documentation; `docs/sources/` – original literature.
- `notebooks/` – Python reference (package `mxpy`, Jupyter notebooks) recomputing Appendix D and Section 6.3; test oracle for the C++ code (`notebooks/results/*.csv`). Edit `notebooks/src/*.py`, then `python notebooks/build_notebooks.py --execute`. See `notebooks/README.md` for results and open questions.

### Build
- Toolchain: llvm-mingw (clang + libc++ + libomp, UCRT; winget `MartinStorsjo.LLVM-MinGW.UCRT`) + Ninja + CMake ≥ 3.24, C++20, OpenMP.
  Toolchain file `cmake/llvm-mingw.cmake` locates clang (env `LLVM_MINGW_ROOT`, PATH, winget folder).
- `powershell -File tools/check.ps1` builds everything (preset `release`), runs ctest, `pip install -e .` and pytest.
- Presets: `cmake --preset release|debug`, `cmake --build --preset release`, `ctest --preset release`.
- Python: `pip install --user --no-build-isolation -e .` (scikit-build-core; rebuilds on import). libc++ is linked statically;
  `libomp.dll` is installed/copied next to the extension and test executables.
- CUDA (later): nvcc 12.8 works with the VS2017 host compiler (C++17 only, tested on the RTX 2070).

### Code conventions (follow existing code)
- Namespace `MultiaxCore`.
- Data-oriented design: SoA, free functions over `std::span`, no allocation in hot loops; kernels marked `MX_HOST_DEVICE` (future CUDA) use no std containers/exceptions.
- The core is **specimen-independent**: no hard-coded specimen data; S-N curves, transfer matrices, load channels are inputs (plane stress sx, sy, txy at hot spots).
- Numeric templates constrained with `std::floating_point TNum`; default alias uses `MXScalar` (`multiax/core/types.hpp`, `double` when `MX_USES_DOUBLES`), e.g. `using WoehlerCurve = WoehlerCurveT<MXScalar>;`.
- Use `MXMachineEpsilon` for tolerance checks.
- Input parameters prefixed `i_`; include guards `__Name_hpp__`.

## Primary source of truth

**`docs/sources/N-515832.pdf`** – A. Bolchoun: *Eine Methode zur Festigkeitsbeurteilung von
laserstrahlgeschweißten Magnesium-Verbindungen unter mehrachsigen Beanspruchungen mit konstanten
und variablen Amplituden*, Dissertation TU Darmstadt 2018, LBF-Bericht FB-253 (ISBN 978-3-8396-1393-1).

This document is the primary source of truth for all formulas, symbols, algorithms and
terminology in this repository.

Rules:
- Implementations follow the dissertation. Cite the equation number in a code comment,
  e.g. `// Eq. (D.53)` or `// Eq. (5.12)`.
- If code, docs or other literature contradict the dissertation, the dissertation wins.
- Suspected misprints in the dissertation are listed in
  `docs/multiax-invariant-appendix-d.md` (section "Suspected misprints"). Use the corrected form
  given there; never silently "fix" anything else – ask first.
- The "Clarifications by the author" in that section fill gaps in the source and are binding.
- Start with the summary `docs/multiax-invariant-appendix-d.md`, then verify details in the PDF.
  **Page offset:** PDF page = printed page + 5 (e.g. Appendix D starts on printed p. 165 = PDF p. 170).
- Key locations (printed pages): Ch. 5 Multiax-Invariant p. 68–83, Appendix B (rainflow, maximum
  variance, nonproportionality) p. 149–158, Appendix C (equivalent stress amplitude as integral)
  p. 159–164, Appendix D (step-by-step procedure, plane stress) p. 165–179, symbols p. 180.
- If something is unclear or not covered by the source, ask instead of guessing.
