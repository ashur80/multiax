# Learning path: data-oriented, parallel C++ for Multiax-Invariant

Milestone I of the C++ core is built as a guided learning path:
- The author implements the core pieces (E1–E6) as exercises.
- Claude provides the scaffold, the tests against the Python reference (`notebooks/mxpy`), and a review of each exercise.

Each exercise has its own file `E<n>_*.md`, created when the exercise starts. It contains:
- **Learning goal** – which concept (DOD, parallelism, …) is practised
- **Background** – method context with equation numbers from the dissertation
- **Task** – interface stub in the code, marked `// TODO(E<n>)`
- **Tests** – doctest/pytest. They are skipped at first (`doctest::skip`) and are switched on when the exercise starts, so they are red first
- **Benchmark** – a nanobench case in `bench/`, results in `benchmarks.md`
- **Review** – Claude reviews the solution (correctness, layout, performance)

| # | Exercise | Focus | Oracle |
|---|---|---|---|
| E1 | S-N curves: `WoehlerCurveT::LofN/NofL` (2.1), then `SNCurveSet` as SoA + interpolation (5.12) | AoS vs SoA, branch-free evaluation, vectorisation | `mxpy/woehler.py` |
| E2 | Cumulative integrals and interval sums → covariances per cycle (D.71–D.73) | flat arrays + index arithmetic instead of object lists | `mxpy/stats.py` |
| E3 | Standard 4-point rainflow count (B.2) | stack on a flat array, no allocation in the hot path | `mxpy/rainflow.py` |
| E4 | Modified rainflow count with time intervals (Sec. 5.5), result as CSR | interval set, CSR layout, edge cases (ties, zero-length intervals) | `mxpy/rainflow.py`, notebook 02 |
| E5 | Rainflow in parallel across many histories (OpenMP) | data races, false sharing, thread-local buffers + prefix sums | E3/E4 serial |
| E6 | One long sequence in parallel: chunks + residue merge | parallel algorithm design, proof of equivalence, scaling | E3/E4 serial |
| (E7) | F_PV kernel (D.64) per cycle (optional) | SoA over cycles, `MX_HOST_DEVICE` kernels | `mxpy/multiax_invariant.py` |

Later: designing the rainflow count on the GPU (CUDA, chunk-parallel per block + residue merge).

## Workflow
```
powershell -File tools/check.ps1              # build + C++ tests + pip install -e . + pytest
powershell -File tools/check.ps1 -Bench       # plus benchmarks
build/release/tests/cpp/multiax_tests.exe -tc="E1*"   # run only the tests of one exercise
```
