"""Reference ROI methods used as benchmark baselines.

* :func:`slice_traces` - fixed log-spaced m/z bins split at scan gaps (an
  "m/z slice" / binned XIC).
* :func:`greedy_roi` - the classic centWave-style ROI builder: scan by scan,
  each centroid joins the nearest open ROI within a fixed ppm window (one
  point per scan), otherwise opens a new ROI.
* :func:`rect_roi` - the original macha ``rectroi``: alternately split by
  m/z gaps and by scan gaps until nothing changes.
"""

from __future__ import annotations

import numba as nb
import numpy as np

from ..run import MSRun


def _split_by_scan_gaps(group, scan, max_gap):
    order = np.lexsort((scan, group))
    g, s = group[order], scan[order]
    new = np.ones(order.size, bool)
    new[1:] = (g[1:] != g[:-1]) | (s[1:] - s[:-1] > max_gap + 1)
    lab = np.cumsum(new) - 1
    out = np.empty_like(lab)
    out[order] = lab
    return out


def slice_traces(run: MSRun, ppm: float = 10.0, max_gap: int = 2) -> np.ndarray:
    bins = np.floor(np.log(run.mz) / np.log1p(ppm * 1e-6)).astype(np.int64)
    return _split_by_scan_gaps(bins, run.scan_of_point.astype(np.int64), max_gap)


@nb.njit(cache=True)
def _greedy(mz, offsets, ppm, max_gap):
    n = mz.shape[0]
    label = np.full(n, -1, np.int64)
    roi_sum = np.zeros(n)
    roi_n = np.zeros(n)
    roi_last = np.zeros(n, np.int64)
    active = np.empty(n, np.int64)
    n_active = 0
    n_roi = 0
    means = np.empty(n)
    for s in range(offsets.shape[0] - 1):
        w = 0
        for i in range(n_active):
            r = active[i]
            if s - roi_last[r] - 1 <= max_gap:
                active[w] = r
                w += 1
        n_active = w
        for i in range(n_active):
            means[i] = roi_sum[active[i]] / roi_n[active[i]]
        # sort active by mean (insertion sort)
        for i in range(1, n_active):
            km = means[i]
            kr = active[i]
            j = i - 1
            while j >= 0 and means[j] > km:
                means[j + 1] = means[j]
                active[j + 1] = active[j]
                j -= 1
            means[j + 1] = km
            active[j + 1] = kr
        n_before = n_active
        for p in range(offsets[s], offsets[s + 1]):
            x = mz[p]
            tol = x * ppm * 1e-6
            # binary search nearest
            lo, hi = 0, n_before
            while lo < hi:
                mid = (lo + hi) // 2
                if means[mid] < x:
                    lo = mid + 1
                else:
                    hi = mid
            best = -1
            bd = tol
            for j in (lo - 1, lo):
                if 0 <= j < n_before:
                    r = active[j]
                    d = abs(means[j] - x)
                    if d <= bd and roi_last[r] != s:
                        bd = d
                        best = r
            if best < 0:
                best = n_roi
                n_roi += 1
                roi_sum[best] = 0.0
                roi_n[best] = 0.0
                active[n_active] = best
                n_active += 1
            roi_sum[best] += x
            roi_n[best] += 1
            roi_last[best] = s
            label[p] = best
    return label


def greedy_roi(run: MSRun, ppm: float = 5.0, max_gap: int = 2) -> np.ndarray:
    return _greedy(run.mz, run.offsets, ppm, max_gap)


def rect_roi(run: MSRun, ppm: float = 8.0, max_gap: int = 4, max_iter: int = 50) -> np.ndarray:
    mz = run.mz
    scan = run.scan_of_point.astype(np.int64)
    g = np.zeros(mz.size, np.int64)
    prev = -1
    for _ in range(max_iter):
        order = np.lexsort((mz, g))
        gm, m = g[order], mz[order]
        new = np.ones(order.size, bool)
        new[1:] = (gm[1:] != gm[:-1]) | ((m[1:] - m[:-1]) / m[1:] * 1e6 > ppm)
        lab = np.cumsum(new) - 1
        g = np.empty_like(lab)
        g[order] = lab
        g = _split_by_scan_gaps(g, scan, max_gap)
        n = int(g.max()) + 1
        if n == prev:
            break
        prev = n
    return g
