"""Cutting-plane geometry and linear stress projections for plane stress.

All projections are linear in (sigma_x, tau_xy, sigma_y), so variances are quadratic forms
w^T C w with the 3x3 covariance matrix C (Eq. (D.50), ordering sx, txy, sy).
"""
import numpy as np
from scipy.spatial import ConvexHull


def normal(theta, psi):
    """Unit normal of a cutting plane, Eq. (D.5)."""
    theta, psi = np.broadcast_arrays(np.asarray(theta, float), np.asarray(psi, float))
    return np.stack([np.cos(theta) * np.cos(psi),
                     np.sin(theta) * np.cos(psi),
                     np.sin(psi)], axis=-1)


def direction(theta, psi, xi):
    """In-plane unit direction r = r1 cos(xi) + r2 sin(xi), Eqs. (D.33), (D.34)."""
    theta, psi, xi = np.broadcast_arrays(*(np.asarray(a, float) for a in (theta, psi, xi)))
    r1 = np.stack([-np.sin(theta), np.cos(theta), np.zeros_like(theta)], axis=-1)
    r2 = np.stack([-np.cos(theta) * np.sin(psi),
                   -np.sin(theta) * np.sin(psi),
                   np.cos(psi)], axis=-1)
    return r1 * np.cos(xi)[..., None] + r2 * np.sin(xi)[..., None]


def w_normal(theta, psi=0.0):
    """Coefficients w with sigma_n = w . (sx, txy, sy).

    Eq. (D.51) corrected: sigma_n = cos^2(psi) (sx cos^2 th + sy sin^2 th + txy sin 2th).
    """
    theta, psi = np.broadcast_arrays(np.asarray(theta, float), np.asarray(psi, float))
    c2 = np.cos(psi) ** 2
    return np.stack([c2 * np.cos(theta) ** 2,
                     c2 * np.sin(2 * theta),
                     c2 * np.sin(theta) ** 2], axis=-1)


def w_shear(theta, psi, xi):
    """Coefficients w with tau_{theta,psi,xi} = w . (sx, txy, sy), Eq. (D.35)/(D.55).

    Since r is orthogonal to n, tau = r^T sigma n.
    """
    n = normal(theta, psi)
    r = direction(theta, psi, xi)
    return np.stack([r[..., 0] * n[..., 0],
                     r[..., 0] * n[..., 1] + r[..., 1] * n[..., 0],
                     r[..., 1] * n[..., 1]], axis=-1)


def quad(w, C):
    """Quadratic form w^T C w, broadcasting over leading axes of w (directions) and C (cases).

    w: (..., m, 3) or (m, 3); C: (n, 3, 3) or (3, 3). Returns (n, m) or (m,).
    """
    C = np.asarray(C, float)
    if C.ndim == 2:
        return np.einsum('mi,ij,mj->m', w, C, w)
    return np.einsum('mi,nij,mj->nm', w, C, w)


class ShearHull:
    """Extreme points of the set {w_shear(theta, psi, xi)} for maximising Var[tau].

    Var[tau] = w^T C w is convex in w, so its maximum over the direction set is attained at a
    vertex of the convex hull of {w} U {-w}. Evaluating only the hull vertices makes the
    maximum-variance search (Eqs. (D.57), (D.82)) cheap for many cases at once.
    """

    def __init__(self, step_deg=1.0):
        st = np.deg2rad(step_deg)
        th = np.arange(0.0, np.pi, st)
        ps = np.arange(-np.pi / 2, np.pi / 2 + 1e-12, st)
        xi = np.arange(0.0, np.pi, st)
        TH, PS, XI = np.meshgrid(th, ps, xi, indexing='ij')
        W = w_shear(TH.ravel(), PS.ravel(), XI.ravel())
        P = np.vstack([W, -W])
        hull = ConvexHull(P)
        v = hull.vertices
        self.W = P[v]
        n = W.shape[0]
        idx = v % n
        self.theta = TH.ravel()[idx]
        self.psi = PS.ravel()[idx]
        self.xi = XI.ravel()[idx]

    def max_var(self, C, chunk=4000):
        """Max over directions of Var[tau] for C of shape (3,3) or (n,3,3).

        Returns (max_var, theta, psi, xi) of the maximising direction (grid resolution).
        """
        C = np.asarray(C, float)
        single = C.ndim == 2
        if single:
            C = C[None]
        out_v = np.empty(len(C))
        out_i = np.empty(len(C), dtype=int)
        for s in range(0, len(C), chunk):
            V = quad(self.W, C[s:s + chunk])
            out_i[s:s + chunk] = np.argmax(V, axis=1)
            out_v[s:s + chunk] = V[np.arange(V.shape[0]), out_i[s:s + chunk]]
        res = (out_v, self.theta[out_i], self.psi[out_i], self.xi[out_i])
        if single:
            return tuple(np.asarray(r)[0] for r in res)
        return res
