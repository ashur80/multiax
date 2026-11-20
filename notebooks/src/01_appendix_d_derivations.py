# %% [markdown]
# # 01 – Appendix D: derivations recomputed (sympy + numerical checks)
#
# Source: A. Bolchoun, Dissertation TU Darmstadt 2018, Appendix D (printed p. 165–179 = PDF p. 170–184).
# Every derivation step of Appendix D is recomputed symbolically with sympy. Each closed form is then
# cross-checked numerically (brute-force integration at random parameter values).
#
# Suspected misprints (documented in `docs/multiax-invariant-appendix-d.md`, Section 8) are shown explicitly:
# **D.7, D.16, D.30, D.51, 5.22** (and D.3, notation only).
#
# Conventions: plane stress, components ordered $(\sigma_x, \tau_{xy}, \sigma_y)$, normal $\mathbf n(\theta,\psi)$ per (D.5).

# %%
import os, sys
sys.path.insert(0, os.getcwd())  # notebooks/ (contains the mxpy package)
import numpy as np
import sympy as sp
from scipy import integrate
from IPython.display import display, Markdown

th, ps, xi, t = sp.symbols('theta psi xi t', real=True)
nu, sI, sII = sp.symbols('nu sigma_I sigma_II', real=True)
h, sig, tau = sp.symbols('h sigma tau', positive=True)
results = []  # (equation, printed, derived, status)


def record(eq, printed, derived, ok, note=''):
    results.append((eq, str(printed), str(derived), 'OK' if ok else 'MISPRINT', note))
    print(f'{eq}: {"OK" if ok else "MISPRINT"}  {note}')


n = sp.Matrix([sp.cos(th) * sp.cos(ps), sp.sin(th) * sp.cos(ps), sp.sin(ps)])        # (D.5)
r1 = sp.Matrix([-sp.sin(th), sp.cos(th), 0])                                           # (D.33)
r2 = sp.Matrix([-sp.cos(th) * sp.sin(ps), -sp.sin(th) * sp.sin(ps), sp.cos(ps)])
r = r1 * sp.cos(xi) + r2 * sp.sin(xi)                                                  # (D.34)
assert sp.simplify(n.dot(n) - 1) == 0 and sp.simplify(r.dot(r) - 1) == 0 and sp.simplify(n.dot(r)) == 0
print('n, r1, r2 orthonormal: OK')


def amp_integral(amp2):
    """(1/pi) sqrt( int_{-pi/2}^{pi/2} int_0^pi amp^2 cos(psi) dtheta dpsi ), Eq. (D.3)."""
    I = sp.integrate(sp.integrate(sp.expand(amp2 * sp.cos(ps)), (th, 0, sp.pi)), (ps, -sp.pi / 2, sp.pi / 2))
    return sp.simplify(sp.sqrt(sp.simplify(I)) / sp.pi)

# %% [markdown]
# ## Weg A1, step 1: base S-N curves → equivalent stress amplitude (D.6)–(D.13)
# Pure axial load $\boldsymbol\sigma_{ax} = \mathrm{diag}(\sigma, \nu\sigma, 0)\sin t$ (D.4).

# %%
S_ax = sp.diag(sig, nu * sig, 0) * sp.sin(t)
sn_ax = sp.simplify((n.T * S_ax * n)[0])
display(sp.Eq(sp.Symbol('sigma_{theta,psi,ax}'), sp.factor(sn_ax)))
amp_ax = sig * sp.cos(ps) ** 2 * (sp.cos(th) ** 2 + nu * sp.sin(th) ** 2)              # amplitude, corrected (D.7)
assert sp.simplify(sn_ax - amp_ax * sp.sin(t)) == 0
printed_D7 = sig * sp.cos(th) ** 2 * (sp.cos(th) ** 2 + nu * sp.sin(th) ** 2)
record('D.7', 'sigma cos^2(theta)(cos^2 th + nu sin^2 th)', 'sigma cos^2(psi)(cos^2 th + nu sin^2 th)',
       sp.simplify(printed_D7 - amp_ax) == 0, 'cos^2(theta) must be cos^2(psi), follows from (D.6)')

se_ax = amp_integral(amp_ax ** 2)
s_ax_printed = sp.sqrt(sp.Rational(2, 15) / sp.pi * (3 * nu ** 2 + 2 * nu + 3))
display(sp.Eq(sp.Symbol('sigma_e,ax^(II)'), se_ax))
ok = sp.simplify(se_ax / sig - s_ax_printed) == 0 or abs(float((se_ax / sig - s_ax_printed).subs({sig: 1, nu: 0.3}))) < 1e-14
record('D.8/D.9', s_ax_printed, se_ax / sig, ok, 's_ax')

# torsion
S_t = sp.Matrix([[0, tau, 0], [tau, 0, 0], [0, 0, 0]]) * sp.sin(t)
sn_t = sp.simplify((n.T * S_t * n)[0])
assert sp.simplify(sn_t - tau * sp.cos(ps) ** 2 * sp.sin(2 * th) * sp.sin(t)) == 0   # (D.10)
se_t = amp_integral((tau * sp.cos(ps) ** 2 * sp.sin(2 * th)) ** 2)
record('D.12/D.13', sp.sqrt(sp.Rational(8, 15) / sp.pi), sp.simplify(se_t / tau),
       sp.simplify(se_t / tau - sp.sqrt(sp.Rational(8, 15) / sp.pi)) == 0, 's_tors')

# %% [markdown]
# ## Weg A1, step 2: proportional load in principal axes (D.14)–(D.16), (D.22)

# %%
amp_p = sp.cos(ps) ** 2 * (sI * sp.cos(th) ** 2 + sII * sp.sin(th) ** 2)              # (D.15)
se_p = amp_integral(amp_p ** 2)
correct = sp.sqrt(sp.Rational(2, 15) / sp.pi * (3 * sI ** 2 + 2 * sI * sII + 3 * sII ** 2))   # (D.22)
printed_D16 = sp.sqrt(sp.Rational(2, 15) / sp.pi * (3 * sI + 2 * sI * sII + 3 * sII))
display(sp.Eq(sp.Symbol('sigma_e,a^(II)'), sp.simplify(se_p)))
record('D.22', correct, se_p, sp.simplify(se_p ** 2 - correct ** 2) == 0)
record('D.16', printed_D16, correct, sp.simplify(printed_D16 ** 2 - correct ** 2) == 0,
       'squares missing on sigma_I, sigma_II (dimensionally inconsistent)')

# (D.17) tau_a,max = max(sigma_I/2, (sigma_I - sigma_II)/2) vs. numerical max over all planes/directions
from mxpy.stress import w_shear, quad
TH, PS, XI = np.meshgrid(np.linspace(0, np.pi, 91), np.linspace(-np.pi / 2, np.pi / 2, 91), np.linspace(0, np.pi, 91), indexing='ij')
Wg = w_shear(TH.ravel(), PS.ravel(), XI.ravel())
for a, b in [(1.0, 0.3), (1.0, -0.5), (1.0, -1.0), (1.0, 1.0)]:
    C = 0.5 * np.diag([a * a, 0.0, b * b]); C[0, 2] = C[2, 0] = 0.5 * a * b
    num = np.sqrt(2 * quad(Wg, C).max())
    assert abs(num - max(a / 2, (a - b) / 2)) < 1e-3, (a, b, num)
print('D.17 tau_a,max = max(sI/2, (sI-sII)/2): OK (numerical check)')

# %% [markdown]
# ## Weg A2, step 2a: 90° calibration load (D.27)–(D.31)
# $\boldsymbol\sigma_{\delta=90°} = \begin{pmatrix}\sigma\sin t & h\sigma\cos t\\ h\sigma\cos t & \nu\sigma\sin t\end{pmatrix}$.

# %%
S90 = sp.Matrix([[sig * sp.sin(t), h * sig * sp.cos(t), 0], [h * sig * sp.cos(t), nu * sig * sp.sin(t), 0], [0, 0, 0]])
sn90 = sp.expand((n.T * S90 * n)[0])
a_sin = sp.simplify(sn90.coeff(sp.sin(t)))
a_cos = sp.simplify(sn90.coeff(sp.cos(t)))
assert sp.simplify(a_sin - sig * sp.cos(ps) ** 2 * (sp.cos(th) ** 2 + nu * sp.sin(th) ** 2)) == 0      # (D.28)
assert sp.simplify(a_cos - sig * sp.cos(ps) ** 2 * h * sp.sin(2 * th)) == 0
amp90_sq = a_sin ** 2 + a_cos ** 2                                                    # (D.29) squared
se90 = amp_integral(amp90_sq)
s90_printed = sp.sqrt(sp.Rational(2, 15) / sp.pi * (3 + 4 * h ** 2 + 2 * nu + 3 * nu ** 2))
record('D.30/D.31', s90_printed, sp.simplify(se90 / sig), sp.simplify((se90 / sig) ** 2 - s90_printed ** 2) == 0, 'final result')
# printed integrand of (D.30): cos^5 psi ((cos^2 th + nu sin^2 th) + (h sin 2th)^2)  (first term not squared)
I_printed = sp.integrate(sp.integrate(sp.cos(ps) ** 5 * ((sp.cos(th) ** 2 + nu * sp.sin(th) ** 2) + (h * sp.sin(2 * th)) ** 2),
                                      (th, 0, sp.pi)), (ps, -sp.pi / 2, sp.pi / 2))
record('D.30 integrand', 'cos^5psi((c^2+nu s^2) + (h sin2th)^2)', 'cos^5psi((c^2+nu s^2)^2 + (h sin2th)^2)',
       sp.simplify(sp.sqrt(I_printed) / sp.pi - s90_printed) == 0, 'printed integrand does not give the printed result')

# %% [markdown]
# ## Weg A2, step 2b: shear stress amplitude (D.35), (D.36)

# %%
tau90 = sp.expand((r.T * S90 * n)[0])
b_sin = tau90.coeff(sp.sin(t))
b_cos = tau90.coeff(sp.cos(t))
amp_tau_sq = b_sin ** 2 + b_cos ** 2
printed_D36_sq = (sig ** 2 / 4 * sp.cos(ps) ** 2 *
                  (((nu - 1) * sp.cos(xi) * sp.sin(2 * th) + ((nu - 1) * sp.cos(2 * th) - (1 + nu)) * sp.sin(xi) * sp.sin(ps)) ** 2
                   + 4 * h ** 2 * (sp.cos(xi) * sp.cos(2 * th) - sp.sin(xi) * sp.sin(ps) * sp.sin(2 * th)) ** 2))
f_diff = sp.lambdify((th, ps, xi, nu, h, sig), amp_tau_sq - printed_D36_sq, 'numpy')
rng = np.random.default_rng(0)
vals = rng.uniform(-3, 3, size=(6, 5000))
err = np.max(np.abs(f_diff(*vals)))
record('D.36', 'printed', 'derived from (D.35)', err < 1e-10, f'max abs diff {err:.1e} (random points)')

# %% [markdown]
# ## Weg A2, step 2d: nonproportionality measure (D.39)–(D.42), constant $C_{PV,\nu}$ (5.22), (5.23), (D.40)

# %%
sx, txy = sig * sp.sin(t), h * sig * sp.cos(t)
integrand = sp.simplify(sp.diff(sx, t) * txy - sx * sp.diff(txy, t))
assert sp.simplify(integrand - h * sig ** 2) == 0                       # |dsx tau - sx dtau| = h sigma^2  (D.41)
record('D.41', '2 pi h sigma^2', sp.integrate(integrand, (t, 0, 2 * sp.pi)), True)

# C_PV,nu = 1/8 int_0^pi |1 - nu + (1 + nu) cos 2 theta| dtheta  (theta integral of (5.22) with the psi integral = 1/2)
a_, b_ = 1 - nu, 1 + nu
th0 = sp.acos(-a_ / b_) / 2                                          # zero of the integrand in (0, pi/2)
f_ = a_ + b_ * sp.cos(2 * th)
I_abs = 2 * (sp.integrate(f_, (th, 0, th0)) - sp.integrate(f_, (th, th0, sp.pi / 2)))
C_derived = sp.simplify(I_abs / 8)
C_printed = (sp.pi * (nu - 1) + 4 * (1 + nu) * sp.sqrt(nu / (1 + nu) ** 2) + 2 * (1 - nu) * sp.acos((nu - 1) / (nu + 1))) / 8
diff = max(abs(float((C_derived - C_printed).subs(nu, v))) for v in (0.05, 0.3, 0.3262, 0.5, 0.9, 1.0))
record('D.40/5.23', C_printed, C_derived, diff < 1e-12, f'max diff over nu samples {diff:.1e}')

# psi weighting in (5.22): with the area weight cos(psi) of (5.17) the psi integral is 1/2 and gives (5.23)
w_with = sp.integrate(sp.cos(ps) ** 3 * sp.Abs(sp.sin(ps)), (ps, -sp.pi / 2, sp.pi / 2))
w_without = sp.integrate(sp.cos(ps) ** 2 * sp.Abs(sp.sin(ps)), (ps, -sp.pi / 2, sp.pi / 2))
print('psi integral with cos(psi) weight:', w_with, ' without:', w_without)
# (5.22) prefactor 1/4: C = 1/4 * w * int|..|dtheta must equal 1/8 int|..|dtheta  ->  w = 1/2
record('5.22', '|cos^2 psi sin psi| (no cos psi weight)', '|cos^2 psi sin psi| cos psi',
       sp.Rational(1, 4) * w_without == sp.Rational(1, 8), f'w = {w_with} (with weight) required for (5.23); printed form gives {w_without}')

# %% [markdown]
# ## Weg A2, step 2e: equivalent amplitude as a function of $F_{PV}$ (D.43), (D.44)

# %%
F = sp.symbols('F_PV', positive=True)
C_ = sp.symbols('C_PV', positive=True)
sxa = sp.symbols('sigma_xa', positive=True)
sol = sp.solve(sp.Eq(F, 2 * sp.pi * h * sxa ** 2 * C_), sxa)
sxa_F = [s for s in sol if s.is_positive is not False][0]
se_F = sp.simplify(sxa_F * s90_printed)
printed_D44 = sp.sqrt(F * (3 + 4 * h ** 2 + 2 * nu + 3 * nu ** 2) / (15 * h * C_)) / sp.pi
chk = float((se_F - printed_D44).subs({F: 1234.5, h: 0.5, nu: 0.33, C_: 0.37}))
record('D.44', printed_D44, se_F, abs(chk) < 1e-12)

# %% [markdown]
# ## Weg A2, step 3: general plane stress (D.48)–(D.53)

# %%
Sx, Sy, Txy = sp.symbols('s_x s_y t_xy', real=True)
S_gen = sp.Matrix([[Sx, Txy, 0], [Txy, Sy, 0], [0, 0, 0]])
sn_gen = sp.simplify((n.T * S_gen * n)[0])
corrected_D51 = sp.cos(ps) ** 2 * (Sx * sp.cos(th) ** 2 + Sy * sp.sin(th) ** 2 + Txy * sp.sin(2 * th))
printed_D51 = sp.cos(ps) * (Sx * sp.cos(th) ** 2 + Sy * sp.sin(th) ** 2 + Txy * sp.sin(2 * th))
assert sp.simplify(sn_gen - corrected_D51) == 0
record('D.51', 'cos(psi)(...)', 'cos^2(psi)(...)', sp.simplify(sn_gen - printed_D51) == 0,
       'n^T sigma n has n_x^2, n_y^2, n_x n_y ~ cos^2 psi')

c11, c22, c33, c12, c13, c23 = sp.symbols('c11 c22 c33 c12 c13 c23', real=True)
Cm = sp.Matrix([[c11, c12, c13], [c12, c22, c23], [c13, c23, c33]])        # ordering (sx, txy, sy)
w = sp.Matrix([sp.cos(ps) ** 2 * sp.cos(th) ** 2, sp.cos(ps) ** 2 * sp.sin(2 * th), sp.cos(ps) ** 2 * sp.sin(th) ** 2])
var_sn = sp.expand((w.T * Cm * w)[0])
printed_D52 = sp.cos(ps) ** 4 * (c11 * sp.cos(th) ** 4 + c22 * sp.sin(2 * th) ** 2 + c33 * sp.sin(th) ** 4
                                  + c12 * 4 * sp.cos(th) ** 3 * sp.sin(th) + c13 * 2 * sp.cos(th) ** 2 * sp.sin(th) ** 2
                                  + c23 * 4 * sp.cos(th) * sp.sin(th) ** 3)
record('D.52', 'printed', 'w^T C w', sp.simplify(sp.expand_trig(var_sn - printed_D52)) == 0)
se_gen = amp_integral(2 * var_sn)
printed_D53 = sp.sqrt(sp.Rational(4, 15) / sp.pi * (3 * c11 + 2 * c13 + 3 * c33 + 4 * c22))
record('D.53', printed_D53, se_gen, sp.simplify(se_gen ** 2 - printed_D53 ** 2) == 0)

# numerical brute-force check of (D.53) for random covariance matrices
from mxpy.multiax_invariant import sigma_e_II, sigma_e_I
for k in range(3):
    A = rng.normal(size=(3, 3)); Cn = A @ A.T
    g = lambda p_, t_: 2 * quad(np.array([[np.cos(p_) ** 2 * np.cos(t_) ** 2, np.cos(p_) ** 2 * np.sin(2 * t_),
                                           np.cos(p_) ** 2 * np.sin(t_) ** 2]]), Cn)[0] * np.cos(p_)
    I, _ = integrate.dblquad(g, 0, np.pi, -np.pi / 2, np.pi / 2, epsabs=1e-12, epsrel=1e-12)
    assert abs(np.sqrt(I) / np.pi - sigma_e_II(Cn)) < 1e-8 * sigma_e_II(Cn)
    # sigma_e^(I) (5.4) by brute force vs. the 1D reduction used in mxpy
    gI = lambda p_, t_: np.sqrt(np.maximum(2 * quad(np.array([[np.cos(p_) ** 2 * np.cos(t_) ** 2, np.cos(p_) ** 2 * np.sin(2 * t_),
                                                              np.cos(p_) ** 2 * np.sin(t_) ** 2]]), Cn)[0], 0)) * np.cos(p_)
    II, _ = integrate.dblquad(gI, 0, np.pi, -np.pi / 2, np.pi / 2, epsabs=1e-10, epsrel=1e-10)
    assert abs(II / np.pi ** 2 - sigma_e_I(Cn)) < 1e-6 * sigma_e_I(Cn), (II / np.pi ** 2, sigma_e_I(Cn))
print('D.53 and (5.4): closed form / 1D reduction == brute-force double integral: OK')

# %% [markdown]
# ## Weg A2, step 4: shear stress in direction ξ and its variance (D.55), (D.56)

# %%
tau_gen = sp.expand((r.T * S_gen * n)[0])
printed_D55 = sp.cos(ps) / 2 * (-Sx * (sp.cos(xi) * sp.sin(2 * th) + 2 * sp.cos(th) ** 2 * sp.sin(xi) * sp.sin(ps))
                                + Sy * (sp.cos(xi) * sp.sin(2 * th) - 2 * sp.sin(xi) * sp.sin(th) ** 2 * sp.sin(ps))
                                + 2 * Txy * (sp.cos(xi) * sp.cos(2 * th) - sp.sin(xi) * sp.sin(2 * th) * sp.sin(ps)))
record('D.55', 'printed', 'r^T sigma n', sp.simplify(sp.expand_trig(tau_gen - printed_D55)) == 0)

wt = sp.Matrix([tau_gen.coeff(Sx), tau_gen.coeff(Txy), tau_gen.coeff(Sy)])
var_tau = (wt.T * Cm * wt)[0]
cs, sn_ = sp.cos, sp.sin
printed_D56 = cs(ps) ** 2 * (c11 * cs(th) ** 2 * (cs(xi) * sn_(th) + cs(th) * sn_(xi) * sn_(ps)) ** 2
                             + c22 * (cs(xi) * cs(2 * th) - sn_(xi) * sn_(2 * th) * sn_(ps)) ** 2
                             + c33 * sn_(th) ** 2 * (cs(xi) * cs(th) - sn_(xi) * sn_(th) * sn_(ps)) ** 2
                             + c12 * (-sp.Rational(1, 2) * cs(xi) ** 2 * sn_(4 * th) + cs(th) ** 2 * sn_(ps) * ((1 - 2 * cs(2 * th)) * sn_(2 * xi) + 2 * sn_(2 * th) * sn_(xi) ** 2 * sn_(ps)))
                             - sp.Rational(1, 4) * c13 * (2 * cs(xi) ** 2 * sn_(2 * th) ** 2 + sn_(ps) * (sn_(2 * xi) * sn_(4 * th) - 2 * sn_(xi) ** 2 * sn_(2 * th) ** 2 * sn_(ps)))
                             + c23 * (sp.Rational(1, 2) * cs(xi) ** 2 * sn_(4 * th) + sn_(th) ** 2 * sn_(ps) * (-sn_(2 * xi) * (1 + 2 * cs(2 * th)) + 2 * sn_(2 * th) * sn_(xi) ** 2 * sn_(ps))))
f56 = sp.lambdify((th, ps, xi, c11, c22, c33, c12, c13, c23), var_tau - printed_D56, 'numpy')
vals = rng.uniform(-3, 3, size=(9, 20000))
err = np.max(np.abs(f56(*vals)))
record('D.56', 'printed', 'w_tau^T C w_tau', err < 1e-10, f'max abs diff {err:.1e} (random points, each c_ij term)')

# %% [markdown]
# ## Weg A2, step 5: swept area and $F_{PV}$ for general plane stress (D.59)–(D.64)
# Check of (D.63): $\tfrac12\|\boldsymbol\tau_{\mathbf n,0}\times\dot{\boldsymbol\tau}_{\mathbf n,0}\| = \tfrac14\cos^2\psi\,|\sin\psi|\,|Q|$,
# with $Q = (\dot\sigma_x-\dot\sigma_y)\tau_{xy}-(\sigma_x-\sigma_y)\dot\tau_{xy}+((\dot\sigma_x+\dot\sigma_y)\tau_{xy}-(\sigma_x+\sigma_y)\dot\tau_{xy})\cos2\theta+(\dot\sigma_x\sigma_y-\sigma_x\dot\sigma_y)\sin2\theta$.
# Integrating $\psi$ out with weight $\cos\psi$ gives the factor $\tfrac18$ of (D.64).

# %%
dSx, dSy, dTxy = sp.symbols('ds_x ds_y dt_xy', real=True)
dS_gen = sp.Matrix([[dSx, dTxy, 0], [dTxy, dSy, 0], [0, 0, 0]])
tv = S_gen * n - (n.T * S_gen * n)[0] * n
dtv = dS_gen * n - (n.T * dS_gen * n)[0] * n
cross_norm = sp.sqrt(sum(c ** 2 for c in tv.cross(dtv))) / 2
Q = ((dSx - dSy) * Txy - (Sx - Sy) * dTxy + ((dSx + dSy) * Txy - (Sx + Sy) * dTxy) * sp.cos(2 * th)
     + (dSx * Sy - Sx * dSy) * sp.sin(2 * th))
printed_D63 = sp.Rational(1, 4) * sp.cos(ps) ** 2 * sp.Abs(sp.sin(ps)) * sp.Abs(Q)
f63 = sp.lambdify((th, ps, Sx, Sy, Txy, dSx, dSy, dTxy), cross_norm - printed_D63, 'numpy')
vals = rng.uniform(-3, 3, size=(8, 20000)); vals[1] = rng.uniform(-np.pi / 2, np.pi / 2, 20000)
err = np.max(np.abs(f63(*vals)))
record('D.63', 'printed', '|tau x dtau|/2', err < 1e-10, f'max abs diff {err:.1e}')
wpsi = sp.integrate(sp.Rational(1, 4) * sp.cos(ps) ** 3 * sp.Abs(sp.sin(ps)), (ps, -sp.pi / 2, sp.pi / 2))
record('D.64', '1/8', wpsi, wpsi == sp.Rational(1, 8), 'psi integral of (D.63) with weight cos(psi)')
# reduction to (5.21)/(5.24) for sigma_y = nu sigma_x
Q_nu = sp.factor(sp.expand(Q.subs({Sy: nu * Sx, dSy: nu * dSx})))
display(sp.Eq(sp.Symbol('Q|_{s_y = nu s_x}'), Q_nu))
assert sp.simplify(Q_nu - (dSx * Txy - Sx * dTxy) * (1 - nu + (1 + nu) * sp.cos(2 * th))) == 0
print('(D.64) reduces to (5.21)/(5.24) for sigma_y = nu sigma_x: OK')

# %% [markdown]
# ## mxpy implementation vs. closed forms
# `mxpy.multiax_invariant.fpv_rate` evaluates the θ-integral of (D.64) in closed form.
# Here it is checked against direct numerical θ-integration.

# %%
from mxpy.multiax_invariant import fpv_rate, C_PV
for k in range(200):
    X0 = rng.normal(size=3); dX = rng.normal(size=3)
    qf = lambda tt: abs((dX[0] - dX[2]) * X0[1] - (X0[0] - X0[2]) * dX[1]
                        + ((dX[0] + dX[2]) * X0[1] - (X0[0] + X0[2]) * dX[1]) * np.cos(2 * tt)
                        + (dX[0] * X0[2] - X0[0] * dX[2]) * np.sin(2 * tt))
    A_ = (dX[0] - dX[2]) * X0[1] - (X0[0] - X0[2]) * dX[1]
    B_ = (dX[0] + dX[2]) * X0[1] - (X0[0] + X0[2]) * dX[1]
    C2 = dX[0] * X0[2] - X0[0] * dX[2]
    R_, phi = np.hypot(B_, C2), np.arctan2(C2, B_)
    kinks = []
    if abs(A_) < R_:  # zeros of A + R cos(2 theta - phi) -> break points for quad
        u0 = np.arccos(-A_ / R_)
        kinks = sorted(((phi + sgn * u0) / 2) % np.pi for sgn in (1, -1))
    Inum = integrate.quad(qf, 0, np.pi, points=kinks or None, limit=400, epsabs=1e-14, epsrel=1e-13)[0] / 8
    assert abs(Inum - fpv_rate(X0[None], dX[None])[0]) < 1e-10 * max(1.0, Inum)
print('closed-form theta integral of (D.64): OK')
for v in (0.1, 0.3262, 0.7):
    t0 = 0.5 * np.arccos(-(1 - v) / (1 + v))
    Iv = integrate.quad(lambda tt: abs(1 - v + (1 + v) * np.cos(2 * tt)), 0, np.pi, points=[t0, np.pi - t0], epsabs=1e-14)[0]
    assert abs(C_PV(v) - Iv / 8) < 1e-12
print('C_PV,nu (mxpy) == 1/8 int |1 - nu + (1 + nu) cos 2 theta|: OK')

# %% [markdown]
# ## Summary: printed vs. derived
# This list matches Section 8 of `docs/multiax-invariant-appendix-d.md`.

# %%
import pandas as pd
summary = pd.DataFrame(results, columns=['equation', 'printed', 'derived', 'status', 'note'])
display(summary[['equation', 'status', 'note']])
expected_misprints = {'D.7', 'D.16', 'D.30 integrand', '5.22', 'D.51'}
found = set(summary.loc[summary.status == 'MISPRINT', 'equation'])
assert found == expected_misprints, found
print('Misprints found:', sorted(found), '-> consistent with docs Section 8')
