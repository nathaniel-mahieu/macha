"""Noise in an ROI's own m/z window.

Whether a match is a coincidence depends on how busy the m/z window it was accepted in is, not on a run-wide or
wide-tile density: a window next to a persistent background ion is occupied in almost every scan, one in an empty
region almost never, and a tile of a few Da averages the two away.  So noise is measured per ROI, in the window the
relation search itself accepts a partner in at the noise floor (``gate`` standard deviations of a noise-level
centroid's m/z error, the ROI's own scatter if wider, plus the systematic term), over the ROI's own scans.
Centroids of the ROI itself are excluded; everything else in the window (unassigned centroids, other traces) is
interference.

Per window: the share of scans holding at least one interfering centroid (the per-scan chance of a spurious
partner), the log-intensity mean and sd of those centroids (the local noise level) and the strongest of them.
"""

from __future__ import annotations

import numba as nb
import numpy as np
import pandas as pd


@nb.njit(cache=True)
def _window_noise(mz, intensity, label, offsets, q_mz, q_hw, q_s0, q_s1, q_ex):
    """Columns: scans, scans hit, hits, mean log intensity of hits, sd log intensity of hits, max intensity."""
    n = q_mz.shape[0]
    n_scans = offsets.shape[0] - 1
    out = np.zeros((n, 6))
    for q in range(n):
        lo_mz = q_mz[q] - q_hw[q]
        hi_mz = q_mz[q] + q_hw[q]
        ns = 0
        nh = 0
        nhit = 0
        s1 = 0.0
        s2 = 0.0
        mx = 0.0
        for s in range(max(q_s0[q], 0), min(q_s1[q], n_scans - 1) + 1):
            ns += 1
            a, b = offsets[s], offsets[s + 1]
            lo, hi = a, b
            while lo < hi:  # first centroid >= lo_mz
                mid = (lo + hi) // 2
                if mz[mid] < lo_mz:
                    lo = mid + 1
                else:
                    hi = mid
            j = lo
            hit = False
            while j < b and mz[j] <= hi_mz:
                if label[j] < 0 or label[j] != q_ex[q]:
                    hit = True
                    nhit += 1
                    li = np.log(intensity[j])
                    s1 += li
                    s2 += li * li
                    if intensity[j] > mx:
                        mx = intensity[j]
                j += 1
            if hit:
                nh += 1
        out[q, 0] = ns
        out[q, 1] = nh
        out[q, 2] = nhit
        if nhit > 0:
            m = s1 / nhit
            out[q, 3] = m
            out[q, 4] = np.sqrt(max(s2 / nhit - m * m, 0.0))
            out[q, 5] = mx
        else:
            out[q, 3] = np.nan
            out[q, 4] = np.nan
            out[q, 5] = np.nan
    return out


def floor_halfwidth_ppm(traces, gate: float = 4.0, sys_ppm: float = 0.3, floor_q: float = 0.01) -> float:
    """Half-width (ppm) inside which the relation search accepts a noise-floor partner."""
    run = traces.run
    sig = float(traces.precision.sigma_ppm(np.quantile(run.intensity, floor_q))) if run.n_points else 5.0
    return gate * np.sqrt(sig * sig + sys_ppm * sys_ppm)


def window_noise(traces, mz, s0, s1, halfwidth_ppm, exclude=None) -> pd.DataFrame:
    """Interference in m/z windows ``mz +/- halfwidth_ppm`` over scans ``s0..s1``.

    ``exclude`` gives, per window, the trace whose own centroids are not interference (-1 for none).
    """
    run = traces.run
    mz = np.atleast_1d(np.asarray(mz, float))
    n = mz.size
    hw = np.broadcast_to(np.asarray(halfwidth_ppm, float), (n,)) * 1e-6 * mz
    ex = np.full(n, -1, np.int64) if exclude is None else np.asarray(exclude, np.int64)
    S = _window_noise(run.mz, run.intensity.astype(np.float64), traces.label.astype(np.int64),
                      run.offsets.astype(np.int64), mz, np.ascontiguousarray(hw),
                      np.broadcast_to(np.asarray(s0, np.int64), (n,)).copy(),
                      np.broadcast_to(np.asarray(s1, np.int64), (n,)).copy(), ex)
    out = pd.DataFrame(S, columns=["n_scans", "hit_scans", "hits", "noise_log_mean", "noise_log_sd", "noise_max"])
    for c in ("n_scans", "hit_scans", "hits"):
        out[c] = out[c].astype(np.int64)
    out["hit_frac"] = out["hit_scans"] / np.maximum(out["n_scans"], 1)
    return out


def roi_noise(traces, gate: float = 4.0, sys_ppm: float = 0.3, pad: int = 20) -> pd.DataFrame:
    """Interference in every ROI's own window: its m/z +/- the noise-floor acceptance width (or its own scatter
    if wider), over its scans +/- ``pad``."""
    t = traces.table
    base = floor_halfwidth_ppm(traces, gate, sys_ppm)
    own = gate * np.sqrt(t["mz_sd_ppm"].fillna(0).to_numpy() ** 2 + sys_ppm * sys_ppm)
    hw = np.maximum(base, own)
    out = window_noise(traces, t["mz"].to_numpy(), t["scan_start"].to_numpy() - pad, t["scan_end"].to_numpy() + pad,
                       hw, np.arange(len(t)))
    out.index = t.index
    out.insert(0, "halfwidth_ppm", hw)
    # apex over the window's typical interference; with no interference, over the run's low-intensity floor
    floor = np.log(np.quantile(traces.run.intensity, 0.05)) if traces.run.n_points else 0.0
    lvl = out["noise_log_mean"].fillna(floor).to_numpy()
    out["snr"] = t["apex_intensity"].to_numpy() / np.exp(lvl)
    return out
