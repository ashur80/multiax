"""Bilinear S-N curves, S-N curve interpolation and damage-equivalent amplitudes."""
from dataclasses import dataclass, replace
import numpy as np


@dataclass(frozen=True)
class SN:
    """Bilinear S-N curve, Eq. (2.1): slope k above and k2 below the knee (L_k, N_k).

    k2 is k* (constant amplitudes) or the Haibach slope k' = 2k - m (variable amplitudes).
    """
    k: float
    L_k: float
    N_k: float
    k2: float

    def N(self, L):
        L = np.asarray(L, float)
        with np.errstate(divide='ignore'):
            return np.where(L >= self.L_k,
                            self.N_k * (self.L_k / L) ** self.k,
                            self.N_k * (self.L_k / L) ** self.k2)

    def L(self, N):
        N = np.asarray(N, float)
        return np.where(N <= self.N_k,
                        self.L_k * (self.N_k / N) ** (1.0 / self.k),
                        self.L_k * (self.N_k / N) ** (1.0 / self.k2))

    def scaled(self, s):
        """Same curve with the load axis multiplied by s (e.g. kN -> MPa)."""
        return replace(self, L_k=self.L_k * s)

    def with_k2(self, k2):
        return replace(self, k2=k2)


def haibach(k, m=2.0):
    """Haibach slope after the knee k' = 2k - m, Eq. (2.23); m = 2 for welded joints."""
    return 2.0 * k - m


def interpolate_params(ax: SN, tors: SN, f_I, interpolate_k2=False):
    """S-N curve interpolation, Eq. (5.12), vectorised over f_I.

    X = 2 (X_ax (1 - f_I) + X_tors (f_I - 1/2)) for X in {k, L_k, N_k} (N_k linearly, Eq. (D.21)).
    f_I < 1/2 -> axial curve; f_I > 1 -> torsion curve (Sec. 5.2, confirmed by the author).
    k2: k* is not interpolated (Weg A); k' is interpolated in Weg B (Eqs. (D.78), (D.84)).
    Returns arrays (k, L_k, N_k, k2).
    """
    f = np.clip(np.asarray(f_I, float), 0.5, 1.0)
    a, b = 2.0 * (1.0 - f), 2.0 * (f - 0.5)
    k = a * ax.k + b * tors.k
    L_k = a * ax.L_k + b * tors.L_k
    N_k = a * ax.N_k + b * tors.N_k
    k2 = a * ax.k2 + b * tors.k2 if interpolate_k2 else np.full_like(f, ax.k2)
    return k, L_k, N_k, k2


def N_of(L, k, L_k, N_k, k2):
    """Vectorised Eq. (2.1) for per-element curve parameters."""
    L = np.asarray(L, float)
    with np.errstate(divide='ignore'):
        return np.where(L >= L_k, N_k * (L_k / L) ** k, N_k * (L_k / L) ** k2)


def L_eq(amplitudes, counts, sn: SN, D_zul=1.0):
    """Damage-equivalent constant amplitude, Eq. (2.25) (with k1 = sn.k, k2 = sn.k2).

    L_eq = [ (1/D_zul) (sum_{L_i >= L_k} n_i L_i^k1 + L_k^(k1-k2) sum_{L_i < L_k} n_i L_i^k2) / sum n_i ]^(1/k1)
    """
    L = np.asarray(amplitudes, float)
    n = np.broadcast_to(np.asarray(counts, float), L.shape)
    above = L >= sn.L_k
    s = np.sum(n[above] * L[above] ** sn.k) + sn.L_k ** (sn.k - sn.k2) * np.sum(n[~above] * L[~above] ** sn.k2)
    return (s / (D_zul * np.sum(n))) ** (1.0 / sn.k)
