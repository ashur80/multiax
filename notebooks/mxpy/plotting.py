"""Plots in the style of Section 6.3: S-N plot, N_exp-N_calc plot, boxplot (Fig. 6.1)."""
import numpy as np
import matplotlib.pyplot as plt
from .data import LABELS, COLORS, MARKERS, LOADCASES


def quantile_ms07(x, q):
    """q-quantile per Eq. (6.3): x_ceil(qn) if qn is not integer, else (x_qn + x_qn+1)/2 (1-based)."""
    x = np.sort(np.asarray(x, float))
    n = len(x)
    qn = q * n
    if abs(qn - round(qn)) < 1e-9:
        j = int(round(qn))
        if j <= 0:
            return x[0]
        if j >= n:
            return x[-1]
        return 0.5 * (x[j - 1] + x[j])
    return x[int(np.ceil(qn)) - 1]


def box_stats(ratios):
    r = np.asarray(ratios, float)
    return {q: quantile_ms07(r, q) for q in (0.1, 0.25, 0.5, 0.75, 0.9)} | {'mean': float(np.mean(r)), 'n': len(r)}


def sn_plot(ax, df, y, ref_curve, T=1.29, N_range=(3e4, 1.2e7), ylabel='Ersatzspannungsamplitude in MPa',
            title=None):
    """S-N plot: tests at (N_B, y) with reference axial curve and scatter band T (P_S = 10/50/90 %)."""
    N = np.logspace(np.log10(N_range[0]), np.log10(N_range[1]), 300)
    L = ref_curve.L(N)
    ax.loglog(N, L, 'k-', lw=1.2, label='Referenzwöhlerlinie (axial)')
    ax.loglog(N, L * np.sqrt(T), 'k-.', lw=0.8)
    ax.loglog(N, L / np.sqrt(T), 'k-.', lw=0.8, label=f'Streuband $T$ = 1:{T}')
    for lc in LOADCASES:
        m = (df['loadcase'] == lc).to_numpy()
        if not m.any():
            continue
        ro = df['runout'].to_numpy() & m
        fa = m & ~ro
        ax.scatter(df['N_B'][fa], y[fa], marker=MARKERS[lc], facecolor='none', edgecolor=COLORS[lc], label=LABELS[lc])
        if ro.any():
            ax.scatter(df['N_B'][ro], y[ro], marker='>', color=COLORS[lc])
    ax.set_xlabel('Bruchschwingspielzahl $N_B$')
    ax.set_ylabel(ylabel)
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=7, loc='upper right')
    if title:
        ax.set_title(title, fontsize=9)


def nn_plot(ax, df, N_calc, lim=(1e3, 1e9), title=None):
    """N_exp - N_calc plot with scatter bands T_N = 1:4 and 1:9."""
    x = np.array(lim)
    ax.loglog(x, x, 'k-', lw=1)
    for f, ls, lab in ((2.0, '--', '$T_N$ = 1:4'), (3.0, ':', '$T_N$ = 1:9')):
        ax.loglog(x, x * f, 'k' + ls, lw=0.8, label=lab)
        ax.loglog(x, x / f, 'k' + ls, lw=0.8)
    for lc in LOADCASES:
        m = (df['loadcase'] == lc).to_numpy()
        if not m.any():
            continue
        ro = df['runout'].to_numpy() & m
        fa = m & ~ro
        ax.scatter(df['N_B'][fa], N_calc[fa], marker=MARKERS[lc], facecolor='none', edgecolor=COLORS[lc], label=LABELS[lc])
        if ro.any():
            ax.scatter(df['N_B'][ro], N_calc[ro], marker='>', color=COLORS[lc])
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel('Experimentelle Lebensdauer $N_{exp}$')
    ax.set_ylabel('Rechnerische Lebensdauer $N_{rechn}$')
    ax.text(0.05, 0.9, 'nichtkonservativ', transform=ax.transAxes, fontsize=8)
    ax.text(0.65, 0.1, 'konservativ', transform=ax.transAxes, fontsize=8)
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=7, loc='lower right')
    if title:
        ax.set_title(title, fontsize=9)


def boxplot(ax, df, N_calc, mask, groups=None, title=None):
    """Boxplot of N_calc/N_exp as in Fig. 6.1: box 0.25-0.75, whiskers 0.1-0.9, median (solid), mean (dashed).

    Returns dict of quantile statistics per group.
    """
    ratio = np.asarray(N_calc, float) / df['N_B'].to_numpy()
    groups = groups or [lc for lc in LOADCASES if (df['loadcase'] == lc).any()] + ['all']
    stats = {}
    for i, g in enumerate(groups):
        m = mask & ((df['loadcase'] == g).to_numpy() if g != 'all' else True)
        if not np.any(m):
            continue
        st = box_stats(ratio[m])
        stats[g] = st
        yv = len(groups) - i
        ax.plot([st[0.1], st[0.9]], [yv, yv], color='gray', lw=1)
        ax.plot([st[0.1]] * 2, [yv - 0.15, yv + 0.15], color='gray')
        ax.plot([st[0.9]] * 2, [yv - 0.15, yv + 0.15], color='gray')
        ax.add_patch(plt.Rectangle((st[0.25], yv - 0.3), st[0.75] - st[0.25], 0.6, color=COLORS[g]))
        ax.plot([st[0.5]] * 2, [yv - 0.3, yv + 0.3], 'w-', lw=1.5)
        ax.plot([st['mean']] * 2, [yv - 0.3, yv + 0.3], 'w--', lw=1.5)
    ax.axvline(1.0, color='k', lw=0.8)
    ax.set_xscale('log')
    ax.set_yticks(range(1, len(groups) + 1))
    ax.set_yticklabels([LABELS[g] for g in reversed(groups)], fontsize=8)
    ax.set_xlabel('Quotient $N_{rechn}/N_{exp}$')
    ax.grid(True, which='both', axis='x', alpha=0.3)
    if title:
        ax.set_title(title, fontsize=9)
    return stats


def gassner_plot(ax, df_va, lines, Mt_over_F=28.0, title=None):
    """Gassner plot (Fig. 6.36, 6.51, 6.54): tests (F_a vs N_B) and lines {label: (N, F, color, ls)}."""
    for lab, (N, F, color, ls) in lines.items():
        ax.loglog(N, F, color=color, ls=ls, label=lab)
    for lc in ('d0', 'd90'):
        m = (df_va['loadcase'] == lc).to_numpy()
        ro = df_va['runout'].to_numpy() & m
        fa = m & ~ro
        ax.scatter(df_va['N_B'][fa], df_va['F_a'][fa], marker=MARKERS[lc], facecolor='none', edgecolor=COLORS[lc])
        ax.scatter(df_va['N_B'][ro], df_va['F_a'][ro], marker='>', color=COLORS[lc])
    ax.set_xlim(1e5, 3e7)
    ax.set_ylim(1.5, 7.5)
    ax.set_xlabel('Bruchschwingspielzahl $N_B$')
    ax.set_ylabel(r'Axiale Kraftamplitude $\bar F_a$ in kN')
    sec = ax.secondary_yaxis('right', functions=(lambda F: F * Mt_over_F, lambda M: M / Mt_over_F))
    sec.set_ylabel(r'Torsionsmomentenamplitude $\bar M_{t,a}$ in Nm')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=7)
    if title:
        ax.set_title(title, fontsize=9)
