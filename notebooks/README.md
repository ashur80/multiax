# Python reference: Appendix D and Section 6.3 recomputed

This folder recomputes the derivations of Appendix D and the results of Section 6.3 (Multiax-Invariant and the code interaction equations) of the dissertation `docs/sources/N-515832.pdf`.
It serves as an independent reference and as the test oracle for the C++ implementation.

| Notebook | Content |
|---|---|
| `01_appendix_d_derivations.ipynb` | All derivation steps of Appendix D, done symbolically (sympy) and checked numerically. Finds exactly the misprints listed in `docs/multiax-invariant-appendix-d.md` §8 |
| `02_rainflow_and_gauss_sequence.ipynb` | Modified rainflow count (Sec. 5.5, Fig. 5.5) and the synthetic Gauss sequence (Sec. 3.4). Validated against Fig. 5.5, the standard rainflow count and Tab. 3.12–3.14 |
| `03_section_6_3_5_multiax_invariant.ipynb` | Tab. 6.5 and Fig. 6.45–6.56 (Multiax-Invariant, constant and variable amplitudes, σ_e^(I)/σ_e^(II)) |
| `04_section_6_3_4_codes.ipynb` | Fig. 6.27–6.44 (IIW, FKM, Eurocode 3) |

The shared code is the package `mxpy/`; every function cites its equation number. The data files are
`data/tests_appendix_a.csv` (Tab. A.2–A.19) and `data/parameters.json` (Tab. 3.8, 3.11, 6.1, k*, m, D_zul, L_TF).
Results per test (oracle for C++) are written to `results/*.csv`.

## Running

```
pip install -r notebooks/requirements.txt
python notebooks/build_notebooks.py --execute     # builds notebooks/*.ipynb from notebooks/src/*.py and executes them (~5 min)
```
Edit the sources in `notebooks/src/*.py` (`# %%` cell format), not the generated `.ipynb`.

## Results

**Reproduced exactly (hard checks via `assert`):**
- Every closed form in Appendix D. Misprints found: D.7, D.16, D.30 (integrand), D.51, 5.22.
- Table 6.5: σ_e^(I), σ_e^(II) and F_PV (4430.90 / 3133.12) to all printed digits.
- The modified rainflow count gives the same cycles as the standard count, and the intervals of the Fig. 5.5 example match the figure.

**Reproduced well (constant amplitudes):** Multiax-Invariant and IIW agree with the statements and figures of Sections 6.3.5 and 6.3.4. Examples:
- δ = 0° is non-conservative with σ_e^(I) and corrected with σ_e^(II).
- δ = 90° is close to the axial results.
- δ = 45° is conservative.

**Deviations / questions to the author:**
1. **f_I in Tab. 6.5:** the full search of (D.57) gives f_I(45°) = 0.50 (table: 0.48). The value 0.48 results only if the search is restricted to ψ = 0, but then f_I(0°) = 0.478 (table: 0.50). The life is not affected, because the axial curve is used for f_I ≤ 0.5 in both cases.
2. **Variable amplitudes (MI):** the recomputed lives are less conservative than in Fig. 6.51–6.56, by a factor of about 1.2 (axial) to 2.
   - About 15 % of this comes from the variance-based cycle amplitude for the synthetic sequence, about 5–10 % from the amplitude classification.
   - The rest is unexplained. Is the original sequence available, and how exactly was d_C computed per cycle?
3. **FKM:** with σ∥ = σ_y and the axial curve as bearable amplitude, FKM is clearly more conservative than IIW for nonproportional loading. Without the σ∥ term the results agree with IIW and with the statement "IIW and FKM very similar". Which bearable amplitude was used for σ∥, or was σ∥ neglected?
4. **EC3, δ = 90°:** (2.41) gives the same N for δ = 0° and 90°, but in Fig. 6.33 the δ = 90° points lie clearly higher. What additional element was used for δ = 90°?
5. **Tab. 3.14, torsion 196 Nm:** D_TF = 0.0903 is not consistent with D_tat = 0.53 (N̄_exp/L_TF = 10.47). The recomputed D_TF ≈ 0.046 gives D_tat ≈ 0.48. Possibly a misprint.
6. **Synthetic sequence:** amplitudes on the upper limits of 32 classes (step height in Fig. 3.10 ≈ 1/32). D_TF of Tab. 3.14 is reproduced to about 0.93 (median).
7. **Runouts in the boxplots (preliminary):** DL are excluded, hDL are counted. Switch: `RUNOUT_POLICY` in notebooks 03/04.
