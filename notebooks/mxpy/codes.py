"""Interaction equations of the design codes (Section 2.8, application in Section 6.3.4).

Local stresses: sigma_x (perpendicular to the weld), sigma_y (parallel), tau_xy at the inner
root notch. Bearable amplitudes are taken from the axial resp. torsion S-N curves at the
life N (normed curves, Tab. 3.8, in local stresses). N is found iteratively.

FKM: sigma_perp = sigma_x, sigma_par = sigma_y, tau = tau_xy; the bearable amplitude for
sigma_par is taken from the axial curve. The 45 deg case is nonproportional for IIW and FKM
(confirmed by the author).
"""
import numpy as np
from scipy.optimize import brentq
from .woehler import SN, haibach


def local_curves(params, amplitudes='CA', notch='inner'):
    """Axial and torsion S-N curves in local stresses (sigma_x resp. tau_xy)."""
    wl = params['woehler_normed_tab_3_8']
    tf = params['transfer_factors_tab_6_1'][notch]
    m, kstar = params['haibach_m'], params['k_star']

    def k2(k):
        return kstar if amplitudes == 'CA' else haibach(k, m)

    ax = SN(wl['ax']['k'], wl['ax']['L_k'] * tf['c_x'], wl['ax']['N_k'], k2(wl['ax']['k']))
    tors = SN(wl['tors']['k'], wl['tors']['L_k'] * tf['c_xy'], wl['tors']['N_k'], k2(wl['tors']['k']))
    return ax, tors


def _solve_N(lhs, rhs=1.0, logN_min=0.0, logN_max=16.0):
    """Solve lhs(N) = rhs for N; lhs increases with N."""
    g = lambda x: lhs(10.0 ** x) - rhs
    if g(logN_min) >= 0:
        return 10.0 ** logN_min
    if g(logN_max) <= 0:
        return 10.0 ** logN_max
    return 10.0 ** brentq(g, logN_min, logN_max, xtol=1e-10)


def life_iiw(s, tau, ax: SN, tors: SN, proportional=True):
    """Modified Gough-Pollard, Eq. (2.34): (s/s_ax(N))^2 + (tau/tau_t(N))^2 = CV, CV = 1 or 0.5."""
    CV = 1.0 if proportional else 0.5
    return _solve_N(lambda N: (s / ax.L(N)) ** 2 + (tau / tors.L(N)) ** 2, CV)


def life_fkm(s_perp, s_par, tau, ax: SN, tors: SN, proportional=True):
    """FKM: Eq. (2.37) (proportional) or Eq. (2.38) (nonproportional) with a = amplitude / bearable amplitude."""
    def lhs(N):
        ap, aq, at = s_perp / ax.L(N), s_par / ax.L(N), tau / tors.L(N)
        if proportional:
            return 0.5 * (abs(ap + aq) + np.sqrt((ap - aq) ** 2 + 4 * at ** 2))
        return ap + aq + at
    return _solve_N(lhs, 1.0)


def life_ec3(s, tau, ax: SN, tors: SN):
    """Eurocode 3, Eq. (2.41): (ds/ds_c(N))^3 + (dtau/dtau_c(N))^5 = 1 (ratios of ranges = ratios of amplitudes)."""
    return _solve_N(lambda N: (s / ax.L(N)) ** 3 + (tau / tors.L(N)) ** 5, 1.0)


def life(code, s_x, s_y, tau, ax, tors, proportional):
    if code == 'IIW':
        return life_iiw(s_x, tau, ax, tors, proportional)
    if code == 'FKM':
        return life_fkm(s_x, s_y, tau, ax, tors, proportional)
    if code == 'EC3':
        return life_ec3(s_x, tau, ax, tors)
    raise ValueError(code)
