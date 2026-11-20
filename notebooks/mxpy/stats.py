"""Mean, variance and covariance over [0,T] or over finite sets of disjoint intervals.

Eqs. (D.24)-(D.26) and (D.68)-(D.73). Sampled signals are integrated with the trapezoidal
rule; integrals over arbitrary intervals use the cumulative integral, linearly interpolated
at the interval bounds.
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid


def cumint(t, y):
    """Cumulative trapezoidal integral of y(t), same length as t (starts at 0)."""
    return cumulative_trapezoid(y, t, initial=0.0)


def interval_integrals(t, cum, a, b, cid, ncycles):
    """Sum over intervals [a_j, b_j] (belonging to cycle cid_j) of the integral given by cum.

    cum may be 1D (one channel) or 2D (n_samples, n_channels). Returns (ncycles,) or (ncycles, n_channels).
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if cum.ndim == 1:
        d = np.interp(b, t, cum) - np.interp(a, t, cum)
        out = np.zeros(ncycles)
        np.add.at(out, cid, d)
        return out
    out = np.zeros((ncycles, cum.shape[1]))
    for c in range(cum.shape[1]):
        d = np.interp(b, t, cum[:, c]) - np.interp(a, t, cum[:, c])
        np.add.at(out[:, c], cid, d)
    return out


def cycle_covariances(t, X, a, b, cid, ncycles):
    """Mean vectors and covariance matrices of the channels X (n_samples, 3) per cycle.

    E[f] = (1/L) sum_i int_{a_i}^{b_i} f dt  (D.71); Var/Cov from E (D.72), (D.73).
    Returns (means (n,3), C (n,3,3), L (n,)).
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    L = np.zeros(ncycles)
    np.add.at(L, cid, b - a)
    m = X.shape[1]
    cum1 = np.column_stack([cumint(t, X[:, i]) for i in range(m)])
    pairs = [(i, j) for i in range(m) for j in range(i, m)]
    cum2 = np.column_stack([cumint(t, X[:, i] * X[:, j]) for i, j in pairs])
    I1 = interval_integrals(t, cum1, a, b, cid, ncycles)
    I2 = interval_integrals(t, cum2, a, b, cid, ncycles)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = I1 / L[:, None]
        C = np.zeros((ncycles, m, m))
        for p, (i, j) in enumerate(pairs):
            C[:, i, j] = C[:, j, i] = I2[:, p] / L - mean[:, i] * mean[:, j]
    return mean, C, L


def covariance_full(t, X):
    """Mean and covariance of channels X over the whole sampled interval [t0, tN]."""
    T = t[-1] - t[0]
    mean = np.array([np.trapezoid(X[:, i], t) for i in range(X.shape[1])]) / T
    Xc = X - mean
    C = np.array([[np.trapezoid(Xc[:, i] * Xc[:, j], t) / T for j in range(X.shape[1])]
                  for i in range(X.shape[1])])
    return mean, C
