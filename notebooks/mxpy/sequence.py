"""Synthetic Gauss load sequence (Section 3.4) and two-channel local stress histories.

Section 3.4: L_TF = 50441 cycles, Gauss-distributed amplitudes (Fig. 3.10), R = -1.
The cumulative spectrum of Fig. 3.10 is H(x) = L_TF^(1 - x^2) with x = L_a / L_a,max
(standard Gauss spectrum; check: H(0.6) ~ 1e3, H(0.8) ~ 50).

Construction of the reversal sequence: one enclosing half cycle -1 -> +1 that contains all
other cycles as symmetric excursions (+x_i, -x_i) in random order. With this nesting the
4-point rainflow count reproduces the prescribed amplitudes exactly (checked in notebook 02),
while the order of the cycles is random.

Between reversal points the signal is a cosine half wave, Eq. (5.25).
A phase shift delta = 90 deg is realised as a time shift of half a reversal interval
(Sec. 5.3 end, Sec. 3.6; confirmed by the author).
"""
import numpy as np


def gauss_amplitudes(L_TF):
    """Amplitudes x_i (descending) with H(x_i) = i, i = 1..L_TF-1 (x = 0 for i = L_TF is dropped)."""
    i = np.arange(1, L_TF)
    return np.sqrt(1.0 - np.log(i) / np.log(L_TF))


def gauss_reversals(L_TF=50441, seed=1, n_classes=32):
    """Reversal point sequence p_0..p_M (p_0 = p_M = 0), normalised to max |p| = 1.

    n_classes: amplitudes are set to the upper limit of n_classes equal amplitude classes
    (Fig. 3.10 is a staircase with a step height of about 1/32). None -> continuous amplitudes.
    With 32 classes the damage sums D_TF of Tab. 3.14 are reproduced within about 10 %
    (notebook 02); with continuous amplitudes they are about 20 % lower.
    """
    x = gauss_amplitudes(L_TF)
    if n_classes:
        x = np.ceil(x * n_classes) / n_classes
    small = x[1:].copy()
    rng = np.random.default_rng(seed)
    rng.shuffle(small)
    inner = np.empty(2 * len(small))
    inner[0::2] = small
    inner[1::2] = -small
    return np.concatenate([[0.0, -1.0], inner, [1.0, 0.0]])


def cosine_interp(p, t):
    """Signal value and time derivative between reversal points p_i at integer times, Eq. (5.25).

    Outside [0, M] the signal is p_0 resp. p_M (both 0 for the sequences used here).
    """
    p = np.asarray(p, float)
    M = len(p) - 1
    t = np.asarray(t, float)
    inside = (t >= 0) & (t <= M)
    tc = np.clip(t, 0, M)
    i = np.minimum(np.floor(tc).astype(int), M - 1)
    frac = tc - i
    a, b = p[i], p[i + 1]
    f = 0.5 * (a - b) * np.cos(frac * np.pi) + 0.5 * (a + b)
    df = -0.5 * (a - b) * np.pi * np.sin(frac * np.pi)
    f = np.where(inside, f, np.where(t < 0, p[0], p[-1]))
    df = np.where(inside, df, 0.0)
    return f, df


def sample_times(M, n_per=16, extra=0.5):
    """Uniform sampling of [0, M + extra] with n_per samples per reversal interval."""
    n = int(round((M + extra) * n_per))
    return np.arange(n + 1) / n_per


def local_channels(p, F, Mt, c_x, c_y, c_xy, shift=0.0, n_per=16):
    """Local plane-stress history (sx, txy, sy) and derivatives for a two-channel test.

    sx = c_x F f(t), sy = c_y F f(t), txy = c_xy Mt f(t - shift); f from Eq. (5.25).
    Returns t, X (n,3), dX (n,3) with columns (sx, txy, sy).
    """
    M = len(p) - 1
    t = sample_times(M, n_per, extra=shift if shift > 0 else 0.0)
    f, df = cosine_interp(p, t)
    g, dg = cosine_interp(p, t - shift)
    X = np.column_stack([c_x * F * f, c_xy * Mt * g, c_y * F * f])
    dX = np.column_stack([c_x * F * df, c_xy * Mt * dg, c_y * F * df])
    return t, X, dX
