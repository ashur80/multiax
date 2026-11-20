"""Rainflow counting.

* rainflow_values: standard 4-point rainflow count (Appendix B.2) on a reversal sequence.
* ModifiedRainflow: modified rainflow count of Sec. 5.5 (flowchart Fig. 5.5), which also
  records the time intervals each cycle consists of (IntervalCache / PointCache / CyclesCache).

Closing condition for the last four points (p_{i-3}, p_{i-2}, p_{i-1}, p_i), Fig. 5.5:
    p_{i-3} >= p_{i-1} > p_{i-2} >= p_i   or   p_{i-3} <= p_{i-1} < p_{i-2} <= p_i
The closed cycle is (p_{i-2}, p_{i-1}).

Residue (author clarification): the residue is appended to itself and the count is repeated,
so the half cycles of the residue are counted as full cycles.
"""
from bisect import bisect_left, bisect_right
import numpy as np


def closes(p3, p2, p1, p0):
    return (p3 >= p1 > p2 >= p0) or (p3 <= p1 < p2 <= p0)


def to_reversals(p):
    """Remove repeated values and intermediate (non-extremal) points from a sequence."""
    out = []
    for x in p:
        if out and x == out[-1]:
            continue
        if len(out) >= 2 and (out[-1] - out[-2]) * (x - out[-1]) > 0:
            out[-1] = x
        else:
            out.append(x)
    return out


def rainflow_values(p, close_residue=True):
    """Standard 4-point rainflow count. Returns list of cycles (p_from, p_to) and the residue."""
    stack, cycles = [], []
    for x in to_reversals(p):
        stack.append(x)
        while len(stack) >= 4 and closes(*stack[-4:]):
            cycles.append((stack[-3], stack[-2]))
            del stack[-3:-1]
    residue = list(stack)
    if close_residue and len(residue) > 2:
        st2 = []
        for x in to_reversals(residue + residue):
            st2.append(x)
            while len(st2) >= 4 and closes(*st2[-4:]):
                cycles.append((st2[-3], st2[-2]))
                del st2[-3:-1]
    return cycles, residue


def reversal_indices(y, rtol=1e-10):
    """Indices of the reversal points (local extrema incl. both ends) of a sampled signal.

    Differences below rtol * max|y| are treated as plateaus (avoids spurious reversals from
    rounding at the junctions of residual intervals).
    """
    y = np.asarray(y, float)
    dy = np.diff(y)
    eps = rtol * (np.max(np.abs(y)) if len(y) else 0.0)
    d = np.where(np.abs(dy) <= eps, 0.0, np.sign(dy))
    nz = np.flatnonzero(d)
    if len(nz) == 0:
        return np.array([0, len(y) - 1])
    # forward-fill zero slopes so that plateaus belong to the preceding slope
    idx = np.maximum.accumulate(np.where(d != 0, np.arange(len(d)), 0))
    d = d[idx]
    d[:nz[0]] = d[nz[0]]
    change = np.flatnonzero(d[1:] != d[:-1]) + 1
    return np.concatenate([[0], change, [len(y) - 1]])


class ModifiedRainflow:
    """Modified rainflow count (Sec. 5.5) of a sampled scalar signal y(t).

    After run(), the result is stored in
      p_from, p_to : arrays (n_cycles,) - reversal values of each cycle
      a, b, cid    : flat arrays of the time intervals [a, b] belonging to cycle cid
      pass_id      : 1 = counted in the first pass, 2 = counted from the doubled residue
    """

    def __init__(self, t, y):
        self.t = np.asarray(t, float)
        self.y = np.asarray(y, float)

    # -- helpers ---------------------------------------------------------------------
    @staticmethod
    def _crossing(t, y, i1, i0, level):
        """Time t* in [t_i1, t_i0] with y(t*) = level (y is monotonic between reversals)."""
        ts, ys = t[i1:i0 + 1], y[i1:i0 + 1]
        if ys[-1] < ys[0]:
            ts, ys = ts[::-1], ys[::-1]
        return float(np.interp(level, ys, ts))

    @staticmethod
    def _cut(starts, ends, a, b):
        """Remove [a, b] from the interval set (starts, ends); return the removed pieces."""
        j0 = bisect_right(ends, a)
        j1 = bisect_left(starts, b)
        last = len(starts) - 1
        removed, ns, ne = [], [], []
        for j in range(j0, j1):
            s, e = starts[j], ends[j]
            ps, pe = max(s, a), min(e, b)
            if pe > ps:
                removed.append((ps, pe))
            if s < a:
                ns.append(s)
                ne.append(a)
            # keep the remainder; the last (growing) interval is kept even with zero length,
            # otherwise the next reversal would extend an older interval over counted parts
            if e > b or (j == last and e >= b):
                ns.append(b)
                ne.append(e)
        starts[j0:j1] = ns
        ends[j0:j1] = ne
        return removed

    def _pass(self, t, y):
        ridx = reversal_indices(y)
        starts, ends = [t[ridx[0]]], [t[ridx[0]]]  # IntervalCache
        P = []                                      # PointCache (sample indices)
        cycles = []                                 # CyclesCache
        for i in ridx:
            P.append(i)
            ends[-1] = t[i]
            while len(P) >= 4:
                i3, i2, i1, i0 = P[-4:]
                if not closes(y[i3], y[i2], y[i1], y[i0]):
                    break
                tstar = self._crossing(t, y, i1, i0, y[i2])
                removed = self._cut(starts, ends, t[i2], tstar)
                cycles.append((y[i2], y[i1], removed))
                del P[-3:-1]
        intervals = [(s, e) for s, e in zip(starts, ends) if e > s]
        return cycles, P, intervals

    def _residue_signal(self, intervals):
        """Concatenate the residual intervals twice into a virtual signal (v, yv) and a piece table."""
        t, y = self.t, self.y
        vs, ys, pieces = [], [], []   # pieces: (v_start, v_end, t_start)
        v0 = 0.0
        for _rep in range(2):
            for (s, e) in intervals:
                inner = np.flatnonzero((t > s) & (t < e))
                tt = np.concatenate([[s], t[inner], [e]])
                yy = np.concatenate([[np.interp(s, t, y)], y[inner], [np.interp(e, t, y)]])
                vv = v0 + (tt - s)
                if vs and abs(yy[0] - ys[-1][-1]) <= 1e-10 * max(1.0, abs(yy[0])):
                    yy[0] = ys[-1][-1]   # identical value at the junction (cut at the same level)
                if vs and vv[0] <= vs[-1][-1]:
                    vv = vv + 1e-9  # keep virtual time strictly increasing at junctions
                vs.append(vv)
                ys.append(yy)
                pieces.append((vv[0], vv[-1], s))
                v0 = vv[-1]
        return np.concatenate(vs), np.concatenate(ys), pieces

    @staticmethod
    def _map_back(removed, pieces):
        out = []
        pv0 = np.array([p[0] for p in pieces])
        pv1 = np.array([p[1] for p in pieces])
        pt0 = np.array([p[2] for p in pieces])
        for (va, vb) in removed:
            j = np.flatnonzero((pv1 > va) & (pv0 < vb))
            for jj in j:
                ps, pe = max(va, pv0[jj]), min(vb, pv1[jj])
                if pe > ps:
                    out.append((pt0[jj] + ps - pv0[jj], pt0[jj] + pe - pv0[jj]))
        return out

    # -- main ------------------------------------------------------------------------
    def run(self):
        cycles1, P, intervals = self._pass(self.t, self.y)
        self.residue_points = [(self.y[i], self.t[i]) for i in P]
        self.residue_intervals = intervals
        cycles2 = []
        if len(P) > 2:
            v, yv, pieces = self._residue_signal(intervals)
            c2, _, _ = self._pass(v, yv)
            cycles2 = [(pf, pt, self._map_back(rem, pieces)) for (pf, pt, rem) in c2]
        allc = [(c, 1) for c in cycles1] + [(c, 2) for c in cycles2]
        self.p_from = np.array([c[0] for c, _ in allc])
        self.p_to = np.array([c[1] for c, _ in allc])
        self.pass_id = np.array([pid for _, pid in allc])
        a, b, cid = [], [], []
        for k, (c, _) in enumerate(allc):
            for (s, e) in c[2]:
                a.append(s)
                b.append(e)
                cid.append(k)
        self.a, self.b, self.cid = np.array(a), np.array(b), np.array(cid, dtype=int)
        self.n_cycles = len(allc)
        return self

    def cycle_sample_slices(self):
        """For each cycle, list of (i_start, i_end) sample index ranges strictly inside its intervals."""
        i0 = np.searchsorted(self.t, self.a, side='right')
        i1 = np.searchsorted(self.t, self.b, side='left')
        out = [[] for _ in range(self.n_cycles)]
        for k in range(len(self.a)):
            out[self.cid[k]].append((k, i0[k], i1[k]))
        return out
