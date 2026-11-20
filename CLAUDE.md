# CLAUDE.md

## Project

`multiax` is a C++ library for fatigue life assessment under multiaxial loading
with constant and variable amplitudes. It implements the method **Multiax-Invariant**.

- `multiax_core/` – core library (`incl/` headers, `src/` sources). DLL export macro: `MULTIAXCORE_DLL` → `MULTIAXCORE_DLLEXPORTS` (`MultiaxCoreConfig.hpp`).
- `multiax/` – application / driver (`src/main.cpp`), solution file `multiax/multiax.sln`.
- `glm/` – third-party math library (vendored). Do not modify.
- `docs/` – method documentation; `docs/sources/` – original literature.

Build: Visual Studio / MSBuild via `multiax/multiax.sln` (C++20, uses `<concepts>`).

### Code conventions (follow existing code)
- Namespace `MultiaxCore`.
- Numeric templates constrained with `std::floating_point TNum`; default alias uses `MXScalar` (`mxtypes.hpp`, `double` when `MX_USES_DOUBLES`), e.g. `using WoehlerCurve = WoehlerCurveT<MXScalar>;`.
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
