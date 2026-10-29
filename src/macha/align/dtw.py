"""Banded dynamic time warping (numba)."""

from __future__ import annotations

import numba as nb
import numpy as np


@nb.njit(cache=True)
def dtw_path(a, b, band):
    """DTW between 1-D series with a Sakoe-Chiba band (in samples, after length scaling).

    Returns the warping path as two index arrays and the total cost.
    """
    n, m = a.shape[0], b.shape[0]
    INF = 1e300
    D = np.full((n + 1, m + 1), INF)
    D[0, 0] = 0.0
    ratio = m / n
    for i in range(1, n + 1):
        c = int((i - 1) * ratio)
        lo = max(1, c - band)
        hi = min(m, c + band + 1)
        for j in range(lo, hi + 1):
            cost = abs(a[i - 1] - b[j - 1])
            best = D[i - 1, j - 1]
            if D[i - 1, j] < best:
                best = D[i - 1, j]
            if D[i, j - 1] < best:
                best = D[i, j - 1]
            D[i, j] = cost + best
    # backtrack
    pi = np.empty(n + m, np.int64)
    pj = np.empty(n + m, np.int64)
    i, j = n, m
    k = 0
    while i > 0 and j > 0:
        pi[k] = i - 1
        pj[k] = j - 1
        k += 1
        d = D[i - 1, j - 1]
        u = D[i - 1, j]
        lft = D[i, j - 1]
        if d <= u and d <= lft:
            i -= 1
            j -= 1
        elif u <= lft:
            i -= 1
        else:
            j -= 1
    return pi[:k][::-1], pj[:k][::-1], D[n, m]


def _prep(x):
    x = np.sqrt(np.maximum(np.asarray(x, float), 0.0))
    return (x - x.mean()) / (x.std() + 1e-12)


def dtw_warp(rt_query, y_query, rt_ref, y_ref, band_seconds: float = 60.0, smooth: int = 25):
    """Warp function mapping query retention times onto the reference time axis."""
    a, b = _prep(y_query), _prep(y_ref)
    step = float(np.median(np.diff(rt_ref)))
    pi, pj, _ = dtw_path(a, b, int(band_seconds / step))
    tr = rt_ref[pj]
    # average reference time per query index, smooth, enforce monotonicity
    uq, inv = np.unique(pi, return_inverse=True)
    mean_tr = np.bincount(inv, weights=tr) / np.bincount(inv)
    if smooth > 1 and mean_tr.size > smooth:
        k = np.ones(smooth) / smooth
        pad = np.pad(mean_tr, smooth // 2, mode="edge")
        mean_tr = np.convolve(pad, k, mode="valid")[: uq.size]
    mean_tr = np.maximum.accumulate(mean_tr)
    xq = rt_query[uq]

    def f(t):
        return np.interp(t, xq, mean_tr)

    return f
