"""mxpy - Python reference implementation of the Multiax-Invariant method and the code
interaction equations, used to recompute Appendix D and Section 6.3 of

    A. Bolchoun, Dissertation TU Darmstadt 2018, LBF-Bericht FB-253
    (docs/sources/N-515832.pdf - primary source of truth).

Equation numbers in comments refer to the dissertation. Corrected misprints are documented
in docs/multiax-invariant-appendix-d.md, Section 8.

Stress components are always ordered (sigma_x, tau_xy, sigma_y), matching the covariance
numbering c11 = Var[sigma_x], c22 = Var[tau_xy], c33 = Var[sigma_y] of Eq. (D.50).
"""
