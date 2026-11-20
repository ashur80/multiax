# Multiax-Invariant – Fatigue Life Calculation for Plane Stress (Summary of Appendix D)

**Source (primary source of truth):** A. Bolchoun, Dissertation TU Darmstadt 2018, LBF-Bericht FB-253,
`docs/sources/N-515832.pdf`, Appendix D "Lebensdauerberechnung mit dem Verfahren Multiax-Invariant im
Falle eines ebenen Spannungszustandes", printed p. 165–179 (= PDF p. 170–184).
Equations referenced from Chapter 5 (printed p. 68–83) are included where Appendix D relies on them.

Equation numbers in parentheses refer to the dissertation. Symbols follow the dissertation.
Formulas marked ⚠ deviate from the printed text (suspected misprints, see [Section 8](#8-suspected-misprints-and-open-points)).

---

## 1. Scope and conventions

- **Plane stress**, $T$-periodic (constant amplitude) or finite (variable amplitude) local stress history (D.1):

$$
\boldsymbol\sigma(t)=\begin{pmatrix}\sigma_x(t)&\tau_{xy}(t)&0\\ \tau_{xy}(t)&\sigma_y(t)&0\\ 0&0&0\end{pmatrix}
$$

- **No mean stress assessment** in Appendix D ($g_m = 0$, $\bar{\boldsymbol\sigma}=0$, (5.2)). A proposal for mean stresses is in Sec. 5.7 (not part of this procedure).
- **Damage combination** (5.3), (D.66):
  $D_{gesamt}=\sqrt{D_{\sigma_{e,a}}^2 \pm D_{PV}^2}$.
  Use "+" if phase shift **shortens** life (ductile behaviour) and "−" if it **extends** life (brittle). Appendix D writes only "+" and allows $d_{PV}$ to be purely imaginary for the "−" case.
- **Cutting plane normal** (D.5), with $\theta\in[0,\pi)$ and $\psi\in(-\tfrac\pi2,\tfrac\pi2]$:

$$
\mathbf n=\begin{pmatrix}\cos\theta\cos\psi\\ \sin\theta\cos\psi\\ \sin\psi\end{pmatrix}
$$

- **In-plane orthonormal basis and direction** (D.33), (D.34), with $\xi\in[0,\pi)$:

$$
\mathbf r_1=\begin{pmatrix}-\sin\theta\\ \cos\theta\\ 0\end{pmatrix},\quad
\mathbf r_2=\begin{pmatrix}-\cos\theta\sin\psi\\ -\sin\theta\sin\psi\\ \cos\psi\end{pmatrix},\quad
\mathbf r_{\theta,\psi,\xi}=\mathbf r_1\cos\xi+\mathbf r_2\sin\xi
$$

- Normal stress in a plane: $\sigma_{\theta,\psi}(t)=\mathbf n^T\boldsymbol\sigma(t)\mathbf n$.
  Shear stress vector: $\boldsymbol\tau_{\theta,\psi}(t)=\boldsymbol\sigma\mathbf n-(\mathbf n^T\boldsymbol\sigma\mathbf n)\mathbf n$.
  Shear stress in direction $\xi$: $\tau_{\theta,\psi,\xi}(t)=\boldsymbol\tau_{\theta,\psi}^T\mathbf r_{\theta,\psi,\xi}$ (D.35).

---

## 2. Building blocks

### 2.1 Statistics and amplitude definition

Over $[0,T]$ (D.24)–(D.26):

$$
E[f]=\frac1T\int_0^T f\,dt,\qquad
\mathrm{Var}[f]=E\big[(f-E[f])^2\big],\qquad
\mathrm{Cov}[f,g]=E\big[(f-E[f])(g-E[g])\big]
$$

**Amplitude of an arbitrary signal:** $f_a=\sqrt{2\,\mathrm{Var}[f]}$. This is exact for a sine wave (B.9).

**Variable amplitudes (Weg B):** the same definitions are used over a finite set of disjoint intervals
$\Omega=\{[a_1,b_1],\dots,[a_n,b_n]\}$ with total length $L=\sum(b_i-a_i)$ (D.68)–(D.73):

$$
E[f]=\frac1L\sum_{i=1}^n\int_{a_i}^{b_i}f(t)\,dt
$$

Var and Cov follow from $E$.

**Component covariances** (D.50), (D.81):

$$
c_{11}=\mathrm{Var}[\sigma_x],\; c_{22}=\mathrm{Var}[\tau_{xy}],\; c_{33}=\mathrm{Var}[\sigma_y],\;
c_{12}=\mathrm{Cov}[\sigma_x,\tau_{xy}],\; c_{13}=\mathrm{Cov}[\sigma_x,\sigma_y],\; c_{23}=\mathrm{Cov}[\tau_{xy},\sigma_y]
$$

### 2.2 Equivalent stress amplitude $\sigma_{e,a}^{(II)}$ (5.5), (D.3)

$$
\sigma_{e,a}^{(II)}=\sqrt{\frac1{\pi^2}\iint_{\mathbf n}\sigma_a^2(\mathbf n)\,d\mathbf n}
=\frac1\pi\sqrt{\int_{-\pi/2}^{\pi/2}\int_0^\pi 2\,\mathrm{Var}\big[\sigma_{\theta,\psi}(t)\big]\cos\psi\,d\theta\,d\psi}\quad⚠
$$

Normal stress in a plane (D.51) ⚠ and its variance (D.52):

$$
\sigma_{\theta,\psi}(t)=\cos^2\psi\,\big(\sigma_x\cos^2\theta+\sigma_y\sin^2\theta+\tau_{xy}\sin2\theta\big)
$$

$$
\mathrm{Var}[\sigma_{\theta,\psi}]=\cos^4\psi\big(c_{11}\cos^4\theta+c_{22}\sin^22\theta+c_{33}\sin^4\theta
+4c_{12}\cos^3\theta\sin\theta+2c_{13}\cos^2\theta\sin^2\theta+4c_{23}\cos\theta\sin^3\theta\big)
$$

**Closed form** (D.53). This is the formula to implement:

$$
\boxed{\sigma_{e,a}^{(II)}=\sqrt{\frac{4}{15\pi}\big(3c_{11}+2c_{13}+3c_{33}+4c_{22}\big)}}
$$

Special case for proportional loading in principal axes, $\sigma_I g(t)$, $\sigma_{II} g(t)$, with $|\sigma_I|\ge|\sigma_{II}|$ and $\sigma_I>0$ (D.16) ⚠, (D.22):
$\sigma_{e,a}^{(II)}=\sqrt{\tfrac{2}{15\pi}(3\sigma_I^2+2\sigma_I\sigma_{II}+3\sigma_{II}^2)}$.

### 2.3 Maximum amplitudes and interpolation factor $f_I$ (maximum variance method)

$$
\sigma_{a,max}=\sqrt{2\max_{\theta,\psi}\mathrm{Var}[\sigma_{\theta,\psi}(t)]}\;\;(D.54),\qquad
\tau_{a,max}=\sqrt{2\max_{\theta,\psi,\xi}\mathrm{Var}[\tau_{\theta,\psi,\xi}(t)]}\;\;(D.57)
$$

$$
f_I=\frac{\tau_{a,max}}{\sigma_{a,max}}\qquad(5.11)
$$

Shear stress in direction $\xi$ for the general plane stress history (D.55):

$$
\tau_{\theta,\psi,\xi}(t)=\tfrac12\cos\psi\Big[-\sigma_x(\cos\xi\sin2\theta+2\cos^2\theta\sin\xi\sin\psi)
+\sigma_y(\cos\xi\sin2\theta-2\sin\xi\sin^2\theta\sin\psi)
+2\tau_{xy}(\cos\xi\cos2\theta-\sin\xi\sin2\theta\sin\psi)\Big]
$$

Its variance (D.56), identical to (D.80):

$$
\begin{aligned}
\mathrm{Var}[\tau_{\theta,\psi,\xi}]=\cos^2\psi\Big[\,&c_{11}\cos^2\theta(\cos\xi\sin\theta+\cos\theta\sin\xi\sin\psi)^2
+c_{22}(\cos\xi\cos2\theta-\sin\xi\sin2\theta\sin\psi)^2\\
&+c_{33}\sin^2\theta(\cos\xi\cos\theta-\sin\xi\sin\theta\sin\psi)^2\\
&+c_{12}\Big(-\tfrac12\cos^2\xi\sin4\theta+\cos^2\theta\sin\psi\big((1-2\cos2\theta)\sin2\xi+2\sin2\theta\sin^2\xi\sin\psi\big)\Big)\\
&-\tfrac14c_{13}\Big(2\cos^2\xi\sin^22\theta+\sin\psi\big(\sin2\xi\sin4\theta-2\sin^2\xi\sin^22\theta\sin\psi\big)\Big)\\
&+c_{23}\Big(\tfrac12\cos^2\xi\sin4\theta+\sin^2\theta\sin\psi\big(-\sin2\xi(1+2\cos2\theta)+2\sin2\theta\sin^2\xi\sin\psi\big)\Big)\Big]
\end{aligned}
$$

> **Implementation hint:** (D.52) and (D.56) are quadratic forms. $\sigma_{\theta,\psi}$ and $\tau_{\theta,\psi,\xi}$ are linear in
> $(\sigma_x,\tau_{xy},\sigma_y)$ with a coefficient vector $\mathbf w(\theta,\psi[,\xi])$, so $\mathrm{Var}=\mathbf w^T\mathbf C\,\mathbf w$
> with the 3×3 covariance matrix $\mathbf C=(c_{ij})$. This is equivalent to the expanded form above and less error-prone.
> I checked the $c_{11}$ term of (D.56) against (D.55).
> The maxima over $\theta,\psi,\xi$ are found numerically.

### 2.4 Woehler curve interpolation (5.12)

The base S-N curves are pure axial ($ax$) and pure torsion ($tors$). Their parameters (Table D.1) are the slopes $k_{ax},k_{tors}$,
the knee-point amplitudes $\sigma_{x,ax,A},\tau_{xy,tors,A}$, the knee-point cycles $N_{ax,A},N_{tors,A}$, and the slope after the knee $k^*$.

For $\tfrac12\le f_I\le1$, each parameter $X\in\{k,\ \sigma^{(II)}_{e,A},\ N_A\}$ (and $k'$ in Weg B) is interpolated **linearly** (D.19)–(D.21), (D.58):

$$
X_{f_I}=2\Big(X_{ax}(1-f_I)+X_{tors}\big(f_I-\tfrac12\big)\Big)
$$

- $f_I<\tfrac12$: use the axial curve (D.58).
- $f_I=\tfrac12$ gives the axial curve and $f_I=1$ the torsion curve.
- $f_I>1$: use the torsion curve (Sec. 5.2, confirmed by the author). It does not occur in practice.
- In Weg A, $k^*$ is not interpolated. In Weg B, $k'$ is interpolated (D.78), (D.84).

### 2.5 Converting the base curves to $\sigma_{e,a}^{(II)}$

All S-N curves must be expressed in terms of $\sigma_{e,a}^{(II)}$. Only the knee-point amplitudes change, by a scaling factor. Slopes and $N_A$ stay the same.

| Load case | Tensor | Scaling factor $s$ | Eq. |
|---|---|---|---|
| Pure axial | $\sigma_x=\sigma_{x,a}\sin t,\ \sigma_y=\nu\sigma_x$ (D.4) | $s_{ax}=\sqrt{\tfrac{2}{15\pi}(3\nu^2+2\nu+3)}$ | (D.8), (D.9) |
| Pure torsion | $\tau_{xy}=\tau_{xy,a}\sin t$ (D.4) | $s_{tors}=\sqrt{\tfrac{8}{15\pi}}$ | (D.12) |
| Phase shift 90° | $\sigma_x=\sigma_{x,a}\sin t,\ \sigma_y=\nu\sigma_x,\ \tau_{xy}=h\sigma_{x,a}\cos t$ (D.27) | $s_{\delta=90°}=\sqrt{\tfrac{2}{15\pi}(3+4h^2+2\nu+3\nu^2)}$ | (D.30), (D.31) |

$\sigma^{(II)}_{e,ax,A}=s_{ax}\,\sigma_{x,ax,A}$, $\sigma^{(II)}_{e,tors,A}=s_{tors}\,\tau_{xy,tors,A}$ (D.13), and
$\sigma^{(II)}_{e,\delta=90°,A}=s_{\delta=90°}\,\sigma_{x,\delta=90°,A}$. Here $h=\tau_{xy,\delta=90°,a}/\sigma_{x,\delta=90°,a}$.

### 2.6 Nonproportionality measure $F_{PV}$

Mean-free shear stress vector (5.15), (D.60)–(D.62):
$\boldsymbol\tau_{\mathbf n,0}(t)=\boldsymbol\tau_{\theta,\psi}(t)-\bar{\boldsymbol\tau}_{\theta,\psi}$, where $\bar{\boldsymbol\tau}_{\theta,\psi}=\bar{\boldsymbol\sigma}\mathbf n-(\mathbf n^T\bar{\boldsymbol\sigma}\mathbf n)\mathbf n$.
The area swept per $dt$ is $dA_{\mathbf n}=\tfrac12\|\boldsymbol\tau_{\mathbf n,0}\times\dot{\boldsymbol\tau}_{\mathbf n,0}\|\,dt$ (5.16).

**General definition** (5.18):

$$
F_{PV}=\frac12\int_{t_1}^{t_2}\int_{-\pi/2}^{\pi/2}\int_0^\pi\|\boldsymbol\tau_{\mathbf n,0}\times\dot{\boldsymbol\tau}_{\mathbf n,0}\|\cos\psi\,d\theta\,d\psi\,dt
$$

**General plane stress, $\psi$ integrated analytically** (D.63), (D.64). Use mean-free components here, see [Section 8](#8-suspected-misprints-and-open-points):

$$
F_{PV}=\frac18\int_{t_1}^{t_2}\int_0^\pi\Big|(\dot\sigma_x-\dot\sigma_y)\tau_{xy}-(\sigma_x-\sigma_y)\dot\tau_{xy}
+\big((\dot\sigma_x+\dot\sigma_y)\tau_{xy}-(\sigma_x+\sigma_y)\dot\tau_{xy}\big)\cos2\theta
+(\dot\sigma_x\sigma_y-\sigma_x\dot\sigma_y)\sin2\theta\Big|\,d\theta\,dt
$$

**Special case $\sigma_y=\nu\sigma_x$** (5.19), (5.23), (5.24), (D.40):

$$
F_{PV}=C_{PV,\nu}\int_{t_1}^{t_2}\big|\dot\sigma^0_x\tau^0_{xy}-\sigma^0_x\dot\tau^0_{xy}\big|\,dt,\qquad
C_{PV,\nu}=\frac18\Big(\pi(\nu-1)+4(1+\nu)\sqrt{\tfrac{\nu}{(1+\nu)^2}}+2(1-\nu)\arccos\tfrac{\nu-1}{\nu+1}\Big)
$$

I checked analytically that $C_{PV,\nu}=\tfrac18\int_0^\pi|1-\nu+(1+\nu)\cos2\theta|\,d\theta$ (e.g. $C_{PV,1}=\tfrac12$).

For the 90° calibration load over one period $[0,2\pi]$ (D.41), (D.42): $F_{PV,\delta=90°}=2\pi h\,\sigma^2_{x,\delta=90°,a}\,C_{PV,\nu}$.

### 2.7 Damage function $d_{PV}(F_{PV})$ (D.44)–(D.47)

For the calibration load (D.27), $\sigma^{(II)}_{e,a}$ and $F_{PV}$ are linked one-to-one. Eliminating $\sigma_{x,\delta=90°,a}$ gives (D.44):

$$
\sigma^{(II)}_{e,\delta=90°,a}(F_{PV})=\frac1\pi\sqrt{F_{PV}\,\frac{3+4h^2+2\nu+3\nu^2}{15\,h\,C_{PV,\nu}}}
$$

Two damages are defined at the same amplitude:
- $d_{exp}$ comes from the **experimental** 90° curve ($k_{\delta=90°}$, $\sigma^{(II)}_{e,\delta=90°,A}$, $N_{\delta=90°,A}$, $k^*$).
- $d_{int}$ comes from the **interpolated** curve with $f_{I,\delta=90°}$ (D.38). It ignores the phase shift.

Then (D.45)–(D.47):

$$
d_{PV}(F_{PV})=\sqrt{d^2_{exp}\big(\sigma^{(II)}_{e,\delta=90°,a}(F_{PV})\big)-d^2_{int}\big(\sigma^{(II)}_{e,\delta=90°,a}(F_{PV})\big)}
$$

The result is real if the phase shift shortens life. If it extends life, $d_{PV}$ is imaginary, which is equivalent to (D.46) and the "−" sign in (5.3).
By construction, (D.47) reproduces the experimental 90° result exactly.

### 2.8 Damage from an S-N curve

$d=1/N(\sigma^{(II)}_{e,a})$ (5.13), using the (interpolated) curve.
Above the knee the slope is $k$; below it, $k^*$ (constant amplitude) or the Haibach slope $k'$ (variable amplitude).

---

## 3. Weg A1 – proportional loading, constant amplitude (D.1.1)

Input: $\boldsymbol\sigma(t)=\mathrm{diag}(\sigma_I,\sigma_{II},0)\,g(t)$ with $|\sigma_I|\ge|\sigma_{II}|$ and $\sigma_I>0$. Here $g$ is $T$-periodic with mean 0 and amplitude 1 (D.2), (D.14).

1. **Convert the base curves** to $\sigma^{(II)}_{e,a}$ using $s_{ax}$ and $s_{tors}$ ([2.5](#25-converting-the-base-curves-to-sigma_eaii)).
2. **Equivalent amplitude:** $\sigma^{(II)}_{e,a}=\sqrt{\tfrac{2}{15\pi}(3\sigma_I^2+2\sigma_I\sigma_{II}+3\sigma_{II}^2)}$ (D.16) ⚠, (D.22).
3. **Load state:** $\sigma_{a,max}=\sigma_I$ and $\tau_{a,max}=\max\big(\tfrac{\sigma_I}2,\tfrac{\sigma_I-\sigma_{II}}2\big)$ (D.17), giving $f_I=\tau_{a,max}/\sigma_{a,max}$.
   In plane stress, $f_I=\tfrac12$ for uniaxial or biaxial tension and $f_I=1$ for pure shear ($\sigma_{II}=-\sigma_I$).
4. **Interpolated curve** for $\tfrac12\le f_I\le1$ (D.18)–(D.21). Use the axial curve for $f_I\le\tfrac12$.
5. **Life:** $N$ from $\sigma^{(II)}_{e,a}$ and the interpolated curve.

> A1 is the special case of A2 with $F_{PV}=0$. The general formulas (D.53), (D.54), (D.57) give the same results.

---

## 4. Weg A2 – nonproportional (phase-shifted) loading, constant amplitude (D.1.2)

Input: a $T$-periodic general plane stress history (D.23) without mean stress assessment.
Additional data: the experimental S-N curve for 90° phase shift, load (D.27), parameters $h,\nu,k_{\delta=90°},\sigma_{x,\delta=90°,A},N_{\delta=90°,A},k^*$.

1. **Convert the base curves** as in A1, step 1.
2. **Calibrate the additional damage** $d_{PV}(F_{PV})$ from the 90° test:
   - 2a. Convert the 90° curve: $\sigma^{(II)}_{e,\delta=90°,A}=s_{\delta=90°}\sigma_{x,\delta=90°,A}$ (D.28)–(D.31).
   - 2b. Find $f_{I,\delta=90°}$. The normal stress amplitude in a plane is
     $\sigma_{\delta=90°,a}(\theta,\psi)=\sigma_{x,a}\cos^2\psi\sqrt{(\cos^2\theta+\nu\sin^2\theta)^2+(h\sin2\theta)^2}$ (D.29); maximise it (D.32).
     The shear stress amplitude is (D.36):
     $\tau_{\delta=90°,a}(\theta,\psi,\xi)=\sigma_{x,a}\tfrac12|\cos\psi|\sqrt{[(\nu-1)\cos\xi\sin2\theta+((\nu-1)\cos2\theta-(1+\nu))\sin\xi\sin\psi]^2+4h^2[\cos\xi\cos2\theta-\sin\xi\sin\psi\sin2\theta]^2}$; maximise it numerically (D.37).
     Both maxima scale with $\sigma_{x,a}$, so $f_{I,\delta=90°}$ depends only on $\nu$ and $h$.
   - 2c. Interpolate the curve with $f_{I,\delta=90°}$ (D.38).
   - 2d. Compute $F_{PV,\delta=90°}=2\pi h\sigma^2_{x,\delta=90°,a}C_{PV,\nu}$ (D.39)–(D.42).
   - 2e. Build $d_{PV}(F_{PV})$ (D.44)–(D.47) as in [2.7](#27-damage-function-d_pvf_pv-d44d47).
3. **Equivalent amplitude** of the load to assess: compute $c_{ij}$ (D.50), then $\sigma^{(II)}_{e,a}$ from (D.53).
4. **Interpolation and amplitude damage:** $\sigma_{a,max}$ (D.54) and $\tau_{a,max}$ (D.57) give $f_I$. Use the axial curve if $f_I<\tfrac12$, otherwise interpolate (D.58).
   Then $d_{\sigma^{(II)}_{e,a}}=1/N(\sigma^{(II)}_{e,a})$.
5. **Nonproportional damage:** compute $F_{PV}$ over one period, using (D.59)–(D.64) in general or (5.24) if $\sigma_y=\nu\sigma_x$.
   Then $d_{PV}=d_{PV}(F_{PV})$ (D.65).
6. **Total:** $d_{gesamt}=\sqrt{d^2_{\sigma^{(II)}_{e,a}}+d^2_{PV}}$ (D.66) and $N=1/d_{gesamt}$ (D.67).

---

## 5. Weg B – variable amplitudes (D.2)

Input: a finite plane stress history (D.74) on $[0,T]$ with several reversal points per component.
Statistics are taken over interval sets ([2.1](#21-statistics-and-amplitude-definition)).

1. **Base curves:** Palmgren-Miner modified by Haibach, with slope after the knee $k'=2k-m$ (2.23).
   Appendix D uses $m=2$ (D.75), giving $k'_{ax}=2k_{ax}-2$ and $k'_{tors}=2k_{tors}-2$ (D.76).
   Convert to $\sigma^{(II)}_{e,a}$ as in A1, step 1.
2. **$d_{PV}(F_{PV})$** as in A2, step 2, with $k'_{\delta=90°}=2k_{\delta=90°}-2$ (D.77).
   The interpolated curve also interpolates $k'$ (D.78). Skip this step if the history is proportional.
3. **Critical plane and direction:** maximise $\mathrm{Var}[\tau_{\theta,\psi,\xi}(t)]$ over $[0,T]$ (D.79)–(D.82). This gives $\theta_{max},\psi_{max},\xi_{max}$.
   The scalar signal $\tau_{\theta_{max},\psi_{max},\xi_{max}}(t)$ drives the rainflow count.
4. **Modified rainflow count** (Sec. 5.5) of that scalar signal. Each cycle $C$ is a set of disjoint time intervals $\Omega_C$ (D.83).
   $N_{TF}$ is the number of identified cycles.
5. **Damage per cycle** $d_C$: apply Weg A1 (proportional) or Weg A2 (nonproportional) to the full tensor restricted to $\Omega_C$.
   The distinction comes from $F_{PV}$ itself: $F_{PV}^C=0$ gives $d_{PV}=0$, which is Weg A1.
   - $E$, Var and Cov are taken over $\Omega_C$ (D.71)–(D.73).
   - $F_{PV}$ is integrated over $\Omega_C$.
   - The interpolation also includes $k'_{f_I}=2(k'_{ax}(1-f_I)+k'_{tors}(f_I-\tfrac12))$ (D.84).
6. **Total damage:** $D_{TF}=\sum_C d_C$ (D.85). The calculated life is (D.86), (2.21):

$$
N_{rech}=\frac{N_{TF}}{D_{TF}}\cdot\frac{D_{zul}}{D_{theor}=1}
$$

### 5.1 Modified rainflow algorithm (Sec. 5.5, Fig. 5.5 flowchart, printed p. 77–80)

Unlike standard rainflow, each cycle also records the **time intervals** it consists of, so the full tensor history per cycle is known.

**Simplifying assumptions:**
- Reversal points $p_0,\dots,p_N$ occur at integer times $0,\dots,N$.
- $p_0=p_N=0$.
- Between neighbouring reversal points the signal is a cosine half-wave (5.25):

$$
f_i(t)=\frac{p_i-p_{i+1}}2\cos\big((t-i)\pi\big)+\frac{p_i+p_{i+1}}2
$$

**Data structures:**
- `IntervalCache`: intervals $(a,b)$ not yet assigned to a cycle.
- `PointCache`: reversal points with their times (the residue).
- `CyclesCache`: counted cycles $((p_{von},p_{bis}),Int_1,\dots,Int_k)$.

**Algorithm:**
```
for each reversal point (p_i, t_i):
    push (p_i, t_i) to PointCache
    move the end of the last interval in IntervalCache to t_i
    while PointCache has >= 4 points and its last four points (p_{i-3}, p_{i-2}, p_{i-1}, p_i) satisfy
          p_{i-3} >= p_{i-1} > p_{i-2} >= p_i   or   p_{i-3} <= p_{i-1} < p_{i-2} <= p_i:
        remove (p_{i-2}, t_{i-2}) and (p_{i-1}, t_{i-1}) from PointCache
        in the IntervalCache intervals between t_{i-1} and t_i, find t* with f(t*) = p_{i-2}
        if t_{i-2} and t* lie in the same interval:
            cut [t_{i-2}, t*] out of it (both remaining parts stay in IntervalCache)
        else:
            remove all intervals between them from IntervalCache
            trim the intervals that contain t_{i-2} and t*
        append cycle ((p_{i-2}, p_{i-1}), removed/cut intervals) to CyclesCache
output: CyclesCache; residue = IntervalCache + PointCache
```

**Residue:** the residue (`PointCache` + `IntervalCache`) is appended to itself and the count is repeated, so the residual half cycles are counted as full cycles. This is an author clarification, see [Section 8](#8-suspected-misprints-and-open-points).

---

## 6. Input data

| Data | Symbols | Used in |
|---|---|---|
| Axial base S-N curve (local stresses) | $k_{ax},\ \sigma_{x,ax,A},\ N_{ax,A}$ | A1, A2, B |
| Torsion base S-N curve | $k_{tors},\ \tau_{xy,tors,A},\ N_{tors,A}$ | A1, A2, B |
| Slope after knee (constant amplitude) | $k^*$ | A1, A2 |
| Haibach exponent | $m$ (=2) → $k'=2k-m$ | B |
| 90° phase-shift S-N curve | $k_{\delta=90°},\ \sigma_{x,\delta=90°,A},\ N_{\delta=90°,A},\ h$ | A2, B (nonproportional) |
| Poisson ratio | $\nu$ | all |
| Sign in (5.3) | life shortening (+) or extension (−); Al, Mg: shortening (+) | A2, B |
| Allowable damage sum | $D_{zul}$ | B |

---

## 7. Implementation notes

- **Closed form:** $\sigma^{(II)}_{e,a}$ (D.53) from the $c_{ij}$, the scaling factors $s$, $C_{PV,\nu}$, and the interpolation (5.12).
- **Numerical:**
  - time integrals ($E$, Var, Cov, $F_{PV}$) over $[0,T]$ or $\Omega_C$
  - the $\theta$-integral in $F_{PV}$ (D.64)
  - maximisation over $\theta,\psi$ (D.54) and $\theta,\psi,\xi$ (D.57), (D.82)
  - finding $t^*$ in the rainflow count

  The integrands are smooth except where the absolute value has kinks, so a grid search followed by local refinement is adequate for the maxima.
- **Time derivatives** ($\dot\sigma$) are needed for $F_{PV}$. With the cosine interpolation (5.25) they are analytic between reversal points.
- **N_A interpolation is linear** in $N$ (D.21), not logarithmic. Implement it as written.
- **Mapping to `MultiaxCore::WoehlerCurveT(k, kStar, lK, NK)`:** $lK=\sigma^{(II)}_{e,A}$ (already scaled with $s$), $NK=N_A$, and `kStar` $=k^*$ (Weg A) or $k'$ (Weg B).
  An interpolated curve is a new `WoehlerCurveT` built from the interpolated parameters.

---

## 8. Suspected misprints and open points

The formulas in this summary use the corrected forms. Each item below was checked against neighbouring equations or recomputed.

| Eq. | Printed p. (PDF p.) | Printed | Corrected / interpretation | Reason |
|---|---|---|---|---|
| (D.3), last line | 165 (170) | $2\,\mathrm{Var}\big(\sigma_a^2(\theta,\psi)\big)$ | $2\,\mathrm{Var}\big[\sigma_{\theta,\psi}(t)\big]$ | $\sigma_a^2=2\mathrm{Var}[\sigma_{\theta,\psi}(t)]$; cf. (5.5), (D.53) |
| (D.7) | 166 (171) | $\sigma_{x,ax,a}\cos^2\theta(\cos^2\theta+\nu\sin^2\theta)$ | $\sigma_{x,ax,a}\cos^2\psi(\cos^2\theta+\nu\sin^2\theta)$ | follows from (D.6); (D.8) integrates $\cos^5\psi$ |
| (D.16), result | 168 (173) | $3\sigma_I+2\sigma_I\sigma_{II}+3\sigma_{II}$ | $3\sigma_I^2+2\sigma_I\sigma_{II}+3\sigma_{II}^2$ | dimensionally required; matches (D.22) and recomputation |
| (D.30), integrand | 170 (175) | $\cos^5\psi\,(\cos^2\theta+\nu\sin^2\theta)+(h\sin2\theta)^2$; line 1 without $\cos\psi$ | $\cos^5\psi\big[(\cos^2\theta+\nu\sin^2\theta)^2+(h\sin2\theta)^2\big]$; weight $\cos\psi$ in line 1 | the printed final result is only obtained with the corrected integrand |
| (D.51) | 173 (178) | $\cos\psi(\dots)$ | $\cos^2\psi(\dots)$ | $\mathbf n^T\boldsymbol\sigma\mathbf n$ has $n_x^2,n_y^2,n_xn_y\propto\cos^2\psi$; (D.52) has $\cos^4\psi$ |
| (5.22) | 76 (81) | $\lvert\cos^2\psi\sin\psi\rvert$ without $\cos\psi$ weight | $\lvert\cos^2\psi\sin\psi\rvert\cos\psi$ | (5.23) is only obtained with the area weight $\cos\psi$ from (5.17) (checked: $C_{PV,1}=\tfrac12$) |

**Clarifications by the author** (A. Bolchoun, 2026-11-20). These points are not explicit in the source; treat them as binding:
1. **$f_I>1$:** use the torsion curve (Sec. 5.2). It does not occur in practice.
2. **(D.63)/(D.64):** use the **mean-free** components $\sigma^0_x,\sigma^0_y,\tau^0_{xy}$, as defined by $\boldsymbol\tau_{\mathbf n,0}$ (5.15), (D.62).
3. **Proportional or nonproportional (Weg B, step 5):** the distinction follows from the nonproportionality factor $F_{PV}$ itself.
   - $F_{PV}=0$ gives $d_{PV}=0$, and the per-cycle evaluation reduces to Weg A1.
   - For a globally proportional history, step 2 is skipped.
   - No separate classification criterion is needed.
4. **Residue of the modified rainflow count:** append the residue to itself and repeat the count. The half cycles in the residue are thus counted as full cycles.
5. **Sign in (5.3):** a material property. For aluminium and magnesium, nonproportional loading **shortens** life, so use "+" and real $d_{PV}$.
