"""Censoring-aware baseline and noise estimation for traces.

Orbitrap centroid data has no zeros: intensities below the instrument's
noise threshold are simply not reported.  A missing point is therefore an
*upper bound* (``<= detection floor``), not a zero and not "the last value"
(the original macha rolled the last observation forward).

Model: asymmetric Whittaker smoothing (arPLS, Baek et al. 2015) with a
second-difference penalty, solved in O(n) with a banded LDL^T solver in
numba.  Observed points get arPLS weights (points above the baseline are
down-weighted); censored points pull the baseline down only when it rises
above the detection floor.
"""

from __future__ import annotations

from dataclasses import dataclass

import numba as nb
import numpy as np

from .run import MSRun


# ----------------------------------------------------------------- detection floor
@dataclass
class DetectionFloor:
    """Lowest reportable intensity as a function of (scan, m/z), from a coarse grid."""

    scan_edges: np.ndarray
    mz_edges: np.ndarray
    grid: np.ndarray  # (n_scan_bins, n_mz_bins)

    def __call__(self, scan, mz):
        si = np.clip(np.searchsorted(self.scan_edges, scan, "right") - 1, 0, self.grid.shape[0] - 1)
        mi = np.clip(np.searchsorted(self.mz_edges, mz, "right") - 1, 0, self.grid.shape[1] - 1)
        return self.grid[si, mi]


def estimate_floor(run: MSRun, n_mz_bins: int = 12, scans_per_bin: int = 50, q: float = 0.01) -> DetectionFloor:
    """Low quantile of reported intensities per (scan block, m/z block)."""
    scan = run.scan_of_point
    mz_edges = np.quantile(run.mz, np.linspace(0, 1, n_mz_bins + 1)) if run.n_points else np.array([0, 1])
    mz_edges[0], mz_edges[-1] = -np.inf, np.inf
    scan_edges = np.arange(0, run.n_scans + scans_per_bin, scans_per_bin)
    si = np.clip(np.searchsorted(scan_edges, scan, "right") - 1, 0, len(scan_edges) - 2)
    mi = np.clip(np.searchsorted(mz_edges, run.mz, "right") - 1, 0, n_mz_bins - 1)
    key = si * n_mz_bins + mi
    order = np.lexsort((run.intensity, key))
    ks = key[order]
    grid = np.full((len(scan_edges) - 1) * n_mz_bins, np.nan)
    starts = np.searchsorted(ks, np.arange(grid.size), "left")
    stops = np.searchsorted(ks, np.arange(grid.size), "right")
    cnt = stops - starts
    ok = cnt > 0
    pos = starts[ok] + np.minimum((q * cnt[ok]).astype(np.int64), cnt[ok] - 1)
    grid[ok] = run.intensity[order][pos]
    grid = grid.reshape(len(scan_edges) - 1, n_mz_bins)
    fill = np.nanmedian(grid) if np.isfinite(grid).any() else 1.0
    grid = np.where(np.isfinite(grid), grid, fill)
    return DetectionFloor(scan_edges, mz_edges, grid)


# ----------------------------------------------------------------- Whittaker / arPLS
@nb.njit(cache=True)
def whittaker_solve(y, w, lam):
    """Solve (W + lam D2'D2) z = W y with a pentadiagonal LDL^T factorisation."""
    n = y.shape[0]
    if n < 3:
        return y.copy()
    # banded matrix: diagonal d0, first off-diagonal d1, second off-diagonal d2
    d0 = np.empty(n)
    d1 = np.empty(n)
    d2 = np.empty(n)
    for i in range(n):
        # squared column norm of D2: rows i-2 (coef 1), i-1 (-2), i (1) when within [0, n-3]
        c = 0.0
        if 0 <= i - 2 <= n - 3:
            c += 1.0
        if 0 <= i - 1 <= n - 3:
            c += 4.0
        if i <= n - 3:
            c += 1.0
        d0[i] = w[i] + lam * c
    for i in range(n - 1):
        c = -2.0 if (i == 0 or i == n - 2) else -4.0
        d1[i] = lam * c
    for i in range(n - 2):
        d2[i] = lam
    # LDL^T for symmetric pentadiagonal
    D = np.empty(n)
    L1 = np.zeros(n)
    L2 = np.zeros(n)
    for i in range(n):
        D[i] = d0[i]
        if i >= 1:
            D[i] -= L1[i - 1] ** 2 * D[i - 1]
        if i >= 2:
            D[i] -= L2[i - 2] ** 2 * D[i - 2]
        if i < n - 1:
            v = d1[i]
            if i >= 1:
                v -= L1[i - 1] * L2[i - 1] * D[i - 1]
            L1[i] = v / D[i]
        if i < n - 2:
            L2[i] = d2[i] / D[i]
    b = w * y
    z = np.empty(n)
    for i in range(n):
        v = b[i]
        if i >= 1:
            v -= L1[i - 1] * z[i - 1]
        if i >= 2:
            v -= L2[i - 2] * z[i - 2]
        z[i] = v
    for i in range(n):
        z[i] /= D[i]
    for i in range(n - 1, -1, -1):
        v = z[i]
        if i < n - 1:
            v -= L1[i] * z[i + 1]
        if i < n - 2:
            v -= L2[i] * z[i + 2]
        z[i] = v
    return z


@nb.njit(cache=True)
def _ndtri(p):
    """Inverse standard normal CDF (Acklam's rational approximation, |rel err| < 1.2e-9)."""
    a = (-3.969683028665376e01, 2.209460984245205e02, -2.759285104469687e02, 1.383577518672690e02,
         -3.066479806614716e01, 2.506628277459239e00)
    b = (-5.447609879822406e01, 1.615858368580409e02, -1.556989798598866e02, 6.680131188771972e01,
         -1.328068155288572e01)
    c = (-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e00, -2.549732539343734e00,
         4.374664141464968e00, 2.938163982698783e00)
    d = (7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00, 3.754408661907416e00)
    if p < 0.02425:
        q = np.sqrt(-2 * np.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    if p > 1 - 0.02425:
        q = np.sqrt(-2 * np.log(1 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)


@nb.njit(cache=True)
def tobit_fill(observed, floor, half_window=15):
    """Censored-value imputation from the local fraction of observed scans.

    If the baseline is b with noise sd ~ floor/3, a scan is reported with
    probability p = 1 - Phi((floor - b) / sd), so b = floor - sd * Phi^-1(1 - p).
    Mostly-censored stretches therefore sit near zero and occasionally observed
    stretches just under the floor.
    """
    n = observed.shape[0]
    cs = np.zeros(n + 1)
    for i in range(n):
        cs[i + 1] = cs[i] + (1.0 if observed[i] else 0.0)
    out = np.empty(n)
    for i in range(n):
        a = max(i - half_window, 0)
        b = min(i + half_window + 1, n)
        p = (cs[b] - cs[a]) / (b - a)
        p = min(max(p, 0.002), 0.998)
        v = floor[i] - floor[i] / 3.0 * _ndtri(1.0 - p)
        out[i] = min(max(v, 0.0), floor[i])
    return out


@nb.njit(cache=True)
def censored_arpls(y, observed, floor, lam, n_iter, censored_weight=1.0, fill=-1.0):
    """Baseline of ``y`` where ``observed[i]`` is False for censored points (value <= floor[i]).

    ``fill < 0`` imputes censored points with :func:`tobit_fill`; otherwise with
    ``fill * floor``.
    """
    n = y.shape[0]
    fv = tobit_fill(observed, floor) if fill < 0 else fill * floor
    yy = np.where(observed, y, fv)
    w = np.where(observed, 1.0, censored_weight)
    z = yy.copy()
    for it in range(n_iter):
        z = whittaker_solve(yy, w, lam)
        # arPLS statistics from negative residuals of observed points
        s1 = 0.0
        s2 = 0.0
        cnt = 0
        for i in range(n):
            if observed[i]:
                d = y[i] - z[i]
                if d < 0:
                    s1 += d
                    s2 += d * d
                    cnt += 1
        if cnt > 1:
            m = s1 / cnt
            sd = np.sqrt(max(s2 / cnt - m * m, 1e-12))
        else:
            m = 0.0
            sd = 0.0
        new_w = np.empty(n)
        for i in range(n):
            if observed[i]:
                if sd > 0:
                    arg = 2.0 * (y[i] - z[i] - (2.0 * sd - m)) / sd
                    arg = min(max(arg, -50.0), 50.0)
                    new_w[i] = 1.0 / (1.0 + np.exp(arg))
                else:
                    new_w[i] = 1.0 if y[i] <= z[i] else 0.01
            else:
                # upper bound: matters only when the baseline exceeds the floor
                new_w[i] = 1.0 if z[i] > floor[i] else censored_weight
        if fill < 0:
            # Tobit fraction from baseline-like observations only (peak points do not count)
            fv = tobit_fill(observed & (new_w > 0.5), floor)
        for i in range(n):
            if not observed[i]:
                yy[i] = floor[i] if z[i] > floor[i] else fv[i]
        w = new_w
    return z


@nb.njit(cache=True)
def _noise_from_residuals(y, observed, z, floor):
    """Robust noise sd from first differences of observed points near the baseline."""
    n = y.shape[0]
    diffs = np.empty(n)
    k = 0
    for i in range(1, n):
        if observed[i] and observed[i - 1]:
            if y[i] - z[i] < 3 * floor[i] and y[i - 1] - z[i - 1] < 3 * floor[i]:
                diffs[k] = abs(y[i] - y[i - 1])
                k += 1
    if k < 5:
        return -1.0
    med = np.median(diffs[:k])
    return 1.4826 * med / np.sqrt(2.0)


@nb.njit(cache=True)
def _batch(intensity_flat, offs, observed_flat, floor_flat, lam, n_iter, censored_weight, noise_default, fill):
    nt = offs.shape[0] - 1
    base = np.empty(intensity_flat.shape[0])
    noise = np.empty(nt)
    for t in range(nt):
        a, b = offs[t], offs[t + 1]
        y = intensity_flat[a:b]
        ob = observed_flat[a:b]
        fl = floor_flat[a:b]
        z = censored_arpls(y, ob, fl, lam, n_iter, censored_weight, fill)
        for i in range(b - a):
            if z[i] < 0:
                z[i] = 0.0
        base[a:b] = z
        nz = _noise_from_residuals(y, ob, z, fl)
        noise[t] = nz if nz > 0 else noise_default[t]
    return base, noise


@dataclass
class BaselineResult:
    """Dense, padded per-trace arrays in CSR layout (``offsets`` into the flat arrays)."""

    offsets: np.ndarray
    scans: np.ndarray
    intensity: np.ndarray  # NaN where censored
    baseline: np.ndarray
    floor: np.ndarray
    noise: np.ndarray  # per trace
    point: np.ndarray  # centroid index per dense position, -1 where censored

    def trace(self, i):
        a, b = self.offsets[i], self.offsets[i + 1]
        return dict(scans=self.scans[a:b], intensity=self.intensity[a:b], baseline=self.baseline[a:b],
                    floor=self.floor[a:b], noise=self.noise[i], point=self.point[a:b])


def dense_traces(traces, pad: int = 10):
    """Padded dense intensity arrays for all traces (CSR). Missing scans -> NaN."""
    run = traces.run
    tab = traces.table
    s0 = np.maximum(tab["scan_start"].to_numpy() - pad, 0)
    s1 = np.minimum(tab["scan_end"].to_numpy() + pad, run.n_scans - 1)
    lens = s1 - s0 + 1
    offs = np.concatenate([[0], np.cumsum(lens)])
    tid = np.repeat(np.arange(len(tab)), lens)
    scans = s0[tid] + np.arange(offs[-1]) - offs[:-1][tid]
    inten = np.full(offs[-1], np.nan)
    point = np.full(offs[-1], -1, np.int64)
    idx = traces.order
    lab = traces.label[idx]
    pos = offs[lab] + run.scan_of_point[idx] - s0[lab]
    inten[pos] = run.intensity[idx]
    point[pos] = idx
    return offs, scans, inten, point


def trace_baselines(traces, lam: float = 1e4, n_iter: int = 8, pad: int = 10, floor: DetectionFloor | None = None,
                    censored_weight: float = 1.0, censored_fill: float = -1.0) -> BaselineResult:
    """Baseline and noise for every trace of a :class:`~macha.traces.Traces`."""
    run = traces.run
    floor = floor or estimate_floor(run)
    offs, scans, inten, point = dense_traces(traces, pad)
    mz = np.repeat(traces.table["mz"].to_numpy(), np.diff(offs))
    fl = floor(scans, mz).astype(np.float64)
    observed = np.isfinite(inten)
    y = np.where(observed, inten, 0.0)
    noise_default = floor(traces.table["apex_scan"].to_numpy(), traces.table["mz"].to_numpy()) / 3.0
    base, noise = _batch(y, offs, observed, fl, lam, n_iter, censored_weight, noise_default.astype(np.float64),
                         censored_fill)
    noise = np.maximum(noise, noise_default)
    return BaselineResult(offs, scans, inten, base, fl, noise, point)


# ----------------------------------------------------------------- reference methods
def rolling_min_baseline(y, window: int = 31, smooth: int = 15):
    """Morphological reference: rolling minimum then moving average (NaN -> 0)."""
    from scipy.ndimage import minimum_filter1d, uniform_filter1d

    yy = np.where(np.isfinite(y), y, 0.0)
    return uniform_filter1d(minimum_filter1d(yy, window, mode="nearest"), smooth, mode="nearest")


def als_baseline(y, lam: float = 1e4, p: float = 0.01, n_iter: int = 10, fill: str = "zero"):
    """Classic Eilers ALS; ``fill='zero'`` or ``'locf'`` (roll last observation, as original macha)."""
    yy = np.asarray(y, dtype=np.float64).copy()
    miss = ~np.isfinite(yy)
    if fill == "locf":
        idx = np.where(~miss, np.arange(yy.size), 0)
        np.maximum.accumulate(idx, out=idx)
        yy = yy[idx]
        yy[~np.isfinite(yy)] = 0.0
    else:
        yy[miss] = 0.0
    w = np.ones(yy.size)
    z = yy
    for _ in range(n_iter):
        z = whittaker_solve(yy, w, lam)
        w = np.where(yy > z, p, 1 - p)
    return z
