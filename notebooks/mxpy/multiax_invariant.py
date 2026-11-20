"""Multiax-Invariant: equivalent stress amplitudes, interpolation factor, nonproportionality
measure F_PV, damage function d_PV(F_PV) and the life calculation (Weg A1/A2/B, Appendix D).
"""
import numpy as np
from .stress import w_normal, quad, ShearHull
from .woehler import SN, interpolate_params, N_of, haibach

# ----------------------------------------------------------------------------------------
# Covariances of sinusoidal histories
# ----------------------------------------------------------------------------------------

def cov_sinusoids(amps, phases):
    """Covariance matrix of x_i(t) = a_i sin(t + phi_i): C_ij = a_i a_j cos(phi_i - phi_j) / 2."""
    a = np.asarray(amps, float)
    p = np.asarray(phases, float)
    return 0.5 * np.outer(a, a) * np.cos(p[:, None] - p[None, :])


# ----------------------------------------------------------------------------------------
# Equivalent stress amplitudes
# ----------------------------------------------------------------------------------------

def sigma_e_II(C):
    """sigma_e,a^(II) = sqrt(4/(15 pi) (3 c11 + 2 c13 + 3 c33 + 4 c22)), Eq. (D.53)."""
    C = np.asarray(C, float)
    c11, c22, c33, c13 = C[..., 0, 0], C[..., 1, 1], C[..., 2, 2], C[..., 0, 2]
    return np.sqrt(np.maximum(4.0 / (15.0 * np.pi) * (3 * c11 + 2 * c13 + 3 * c33 + 4 * c22), 0.0))


_NTH = 2880
_TH = np.arange(_NTH) * np.pi / _NTH
_WN = w_normal(_TH, 0.0)


def var_normal_profile(C):
    """g(theta) = Var[sigma_n(theta, psi=0)]; Var[sigma_n(theta, psi)] = cos^4(psi) g(theta) (D.52)."""
    return np.maximum(quad(_WN, C), 0.0)


def _chunked(fun, C, chunk=2000):
    C = np.asarray(C, float)
    if C.ndim == 2:
        return fun(C)
    return np.concatenate([np.atleast_1d(fun(C[s:s + chunk])) for s in range(0, len(C), chunk)])


def sigma_e_I(C):
    """sigma_e,a^(I) = 1/pi^2 int int sqrt(2 Var sigma_n) cos(psi) dtheta dpsi, Eq. (5.4).

    With Var sigma_n = cos^4(psi) g(theta), the psi integral gives int cos^3 = 4/3:
    sigma_e,a^(I) = 4/(3 pi^2) int_0^pi sqrt(2 g(theta)) dtheta (periodic trapezoidal rule).
    """
    f = lambda c: 4.0 / (3.0 * np.pi ** 2) * np.pi * np.mean(np.sqrt(2.0 * var_normal_profile(c)), axis=-1)
    return _chunked(f, C)


def sigma_a_max(C):
    """sigma_a,max = sqrt(2 max Var sigma_n), Eq. (D.54) (maximum at psi = 0)."""
    return _chunked(lambda c: np.sqrt(2.0 * np.max(var_normal_profile(c), axis=-1)), C)


_HULL = None


def shear_hull(step_deg=1.0):
    global _HULL
    if _HULL is None or step_deg != getattr(_HULL, '_step', None):
        _HULL = ShearHull(step_deg)
        _HULL._step = step_deg
    return _HULL


def tau_a_max(C, hull=None):
    """tau_a,max = sqrt(2 max Var tau_{theta,psi,xi}), Eq. (D.57)."""
    hull = hull or shear_hull()
    v = hull.max_var(C)[0]
    return np.sqrt(2.0 * np.maximum(v, 0.0))


def f_I(C, hull=None, eps=1e-12):
    """Interpolation factor f_I = tau_a,max / sigma_a,max, Eq. (5.11)."""
    s = sigma_a_max(C)
    t = tau_a_max(C, hull)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(s > eps, t / np.maximum(s, eps), 0.5)


# ----------------------------------------------------------------------------------------
# Nonproportionality measure F_PV
# ----------------------------------------------------------------------------------------

def C_PV(nu):
    """C_PV,nu, Eq. (5.23)/(D.40)."""
    return (np.pi * (nu - 1) + 4 * (1 + nu) * np.sqrt(nu / (1 + nu) ** 2)
            + 2 * (1 - nu) * np.arccos((nu - 1) / (nu + 1))) / 8.0


def _abs_cos_integral(A, R):
    """int_0^pi |A + R cos(2 theta - phi)| dtheta (closed form, independent of phi)."""
    A = np.abs(np.asarray(A, float))
    R = np.abs(np.asarray(R, float))
    out = np.pi * A
    m = A < R
    Am, Rm = A[m], R[m]
    out[m] = 2.0 * np.sqrt(Rm ** 2 - Am ** 2) + 2.0 * Am * np.arcsin(Am / Rm)
    return out


def fpv_rate(X0, dX):
    """Integrand of F_PV per unit time for mean-free channels X0 = (sx, txy, sy) and dX, Eq. (D.64):

    dF_PV/dt = 1/8 int_0^pi |A + B cos 2th + C sin 2th| dth with
      A = (dsx - dsy) txy - (sx - sy) dtxy,  B = (dsx + dsy) txy - (sx + sy) dtxy,  C = dsx sy - sx dsy.
    The theta integral is evaluated in closed form with R = sqrt(B^2 + C^2).
    """
    sx, txy, sy = X0[..., 0], X0[..., 1], X0[..., 2]
    dsx, dtxy, dsy = dX[..., 0], dX[..., 1], dX[..., 2]
    A = (dsx - dsy) * txy - (sx - sy) * dtxy
    B = (dsx + dsy) * txy - (sx + sy) * dtxy
    Cc = dsx * sy - sx * dsy
    return _abs_cos_integral(A, np.hypot(B, Cc)) / 8.0


def F_PV_sinusoids(amps, phases, n=4096):
    """F_PV over one period [0, 2 pi] of x_i = a_i sin(t + phi_i), Eq. (D.64)."""
    t = np.arange(n) * 2 * np.pi / n
    a = np.asarray(amps, float)
    p = np.asarray(phases, float)
    X = a * np.sin(t[:, None] + p)
    dX = a * np.cos(t[:, None] + p)
    return 2 * np.pi * np.mean(fpv_rate(X, dX))


def F_PV_cycles(t, X, dX, rf, means):
    """F_PV per rainflow cycle, Eq. (D.64), integrated over the cycle's intervals with the
    channels made mean-free with the cycle means (tau_n,0, Eqs. (D.60)-(D.62)). Vectorised."""
    a, b, cid = rf.a, rf.b, rf.cid
    i0 = np.searchsorted(t, a, side='right')
    i1 = np.searchsorted(t, b, side='left')
    cnt = np.maximum(i1 - i0, 0)
    seg = cnt + 2                                   # interior samples + both end points
    k = np.repeat(np.arange(len(a)), seg)           # interval id of every position
    start = np.cumsum(seg) - seg
    j = np.arange(seg.sum()) - start[k]             # position inside the interval
    is_a = j == 0
    is_b = j == seg[k] - 1
    inner = ~(is_a | is_b)
    tt = np.empty(len(k))
    tt[is_a] = a[k[is_a]]
    tt[is_b] = b[k[is_b]]
    idx = i0[k[inner]] + j[inner] - 1
    tt[inner] = t[idx]
    XX = np.empty((len(k), 3))
    dXX = np.empty((len(k), 3))
    for c in range(3):
        XX[:, c] = np.interp(tt, t, X[:, c])
        dXX[:, c] = np.interp(tt, t, dX[:, c])
        XX[inner, c] = X[idx, c]
        dXX[inner, c] = dX[idx, c]
    rate = fpv_rate(XX - means[cid[k]], dXX)
    same = k[1:] == k[:-1]
    contrib = 0.5 * (rate[1:] + rate[:-1]) * (tt[1:] - tt[:-1])
    out = np.zeros(rf.n_cycles)
    np.add.at(out, cid[k[:-1][same]], contrib[same])
    return out


# ----------------------------------------------------------------------------------------
# Life calculation
# ----------------------------------------------------------------------------------------

class MIModel:
    """Multiax-Invariant life model for the tube-tube specimens of the dissertation.

    measure: 'II' -> sigma_e,a^(II) (Eq. (5.5)/(D.53)), 'I' -> sigma_e,a^(I) (Eq. (5.4)).
    amplitudes: 'CA' -> slope k* after the knee, not interpolated (Weg A);
                'VA' -> Haibach slope k' = 2k - m, interpolated (Weg B, Eqs. (D.76)-(D.78), (D.84)).
    """

    def __init__(self, params, measure='II', amplitudes='CA', notch='inner'):
        self.params = params
        self.measure = measure
        self.amplitudes = amplitudes
        tf = params['transfer_factors_tab_6_1'][notch]
        self.c_x, self.c_y, self.c_xy = tf['c_x'], tf['c_y'], tf['c_xy']
        self.r = params['Mt_over_F']
        self.se = sigma_e_II if measure == 'II' else sigma_e_I
        wl = params['woehler_normed_tab_3_8']
        m = params['haibach_m']
        kstar = params['k_star']

        def k2(k):
            return kstar if amplitudes == 'CA' else haibach(k, m)

        # unit equivalent amplitudes: axial per kN, torsion per Nm, combined per kN (M = r F)
        self.unit = {
            'ax': self.se(self.cov(1.0, 0.0, 0.0)),
            'tors': self.se(self.cov(0.0, 1.0, 0.0)),
            'd90': self.se(self.cov(1.0, self.r, 90.0)),
        }
        self.sn_ax = SN(wl['ax']['k'], wl['ax']['L_k'] * self.unit['ax'], wl['ax']['N_k'], k2(wl['ax']['k']))
        self.sn_tors = SN(wl['tors']['k'], wl['tors']['L_k'] * self.unit['tors'], wl['tors']['N_k'], k2(wl['tors']['k']))
        self.sn_d90 = SN(wl['d90']['k'], wl['d90']['L_k'] * self.unit['d90'], wl['d90']['N_k'], k2(wl['d90']['k']))
        # calibration load case delta = 90 deg (Weg A2, step 2), per kN with M = r F
        C90 = self.cov(1.0, self.r, 90.0)
        self.fpv90_unit = self.F_PV_CA(1.0, self.r, 90.0)
        self.fI90 = float(f_I(C90))
        self.interp_k2 = amplitudes == 'VA'

    # local stresses -------------------------------------------------------------------
    def amps_phases(self, F, Mt, delta_deg):
        return ([self.c_x * F, self.c_xy * Mt, self.c_y * F],
                [0.0, -np.deg2rad(delta_deg or 0.0), 0.0])

    def cov(self, F, Mt, delta_deg):
        return cov_sinusoids(*self.amps_phases(F, Mt, delta_deg))

    def F_PV_CA(self, F, Mt, delta_deg):
        return F_PV_sinusoids(*self.amps_phases(F, Mt, delta_deg))

    # S-N curves -------------------------------------------------------------------------
    def N_interp(self, se, fI):
        """Life from the interpolated S-N curve (Eq. (5.12)) at equivalent amplitude se."""
        k, L_k, N_k, k2 = interpolate_params(self.sn_ax, self.sn_tors, fI, self.interp_k2)
        return N_of(se, k, L_k, N_k, k2)

    def dpv2(self, F_PV):
        """d_PV^2(F_PV) = d_exp^2 - d_int^2 at the calibration amplitude, Eqs. (D.44)-(D.47).

        sigma_e,cal(F_PV) = sigma_e,90^unit sqrt(F_PV / F_PV,90^unit) is the general form of (D.44)
        for the actual calibration load (sigma_e ~ F, F_PV ~ F^2). Negative values correspond to an
        imaginary d_PV (life extension); for Al and Mg the source expects positive values.
        """
        F_PV = np.asarray(F_PV, float)
        s = self.unit['d90'] * np.sqrt(np.maximum(F_PV, 0.0) / self.fpv90_unit)
        with np.errstate(divide='ignore'):
            d_exp = np.where(s > 0, 1.0 / self.sn_d90.N(np.where(s > 0, s, 1.0)), 0.0)
            d_int = np.where(s > 0, 1.0 / self.N_interp(np.where(s > 0, s, 1.0), self.fI90), 0.0)
        return d_exp ** 2 - d_int ** 2

    def damage(self, se, fI, F_PV):
        """Damage of one cycle: sqrt(d_sigma^2 + d_PV^2), Eqs. (5.3), (D.66)."""
        d_s = 1.0 / self.N_interp(se, fI)
        return np.sqrt(np.maximum(d_s ** 2 + self.dpv2(F_PV), 0.0))

    def life_CA(self, F, Mt, delta_deg):
        """Weg A1/A2: life under constant amplitudes, Eq. (D.67). Returns (N, se, fI, F_PV)."""
        C = self.cov(F, Mt, delta_deg)
        se = float(self.se(C))
        fI = float(f_I(C))
        fpv = self.F_PV_CA(F, Mt, delta_deg) if delta_deg not in (None, 0, 0.0) else 0.0
        return 1.0 / float(self.damage(se, fI, fpv)), se, fI, fpv

    def reference_amplitude(self, N):
        """Equivalent amplitude on the axial reference S-N curve at life N (y-axis of the S-N plots)."""
        return self.sn_ax.L(N)
