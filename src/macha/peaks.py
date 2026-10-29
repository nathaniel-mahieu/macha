"""Chromatographic peak detection on traces.

Detection uses 0-dimensional topological persistence of the smoothed trace:
every local maximum is born at its height and dies at the saddle where it
merges with a taller one, so persistence equals prominence and is computed
for all maxima in one O(n log n) union-find sweep.  Peaks are kept when both
their prominence and height above the baseline exceed multiples of the
trace noise and enough non-censored points support them.  Boundaries are the
valleys between kept peaks, trimmed to where the signal returns to the
baseline.

:func:`fit_emg_mixture` optionally refines overlapping peaks with a sum of
exponentially modified Gaussians (analytic Jacobian, scipy least squares).
:func:`sspp_peaks` ports the original macha scale-space picker and
:func:`cwt_peaks` wraps scipy's CWT detector, both for benchmarking.
"""

from __future__ import annotations

from dataclasses import dataclass

import numba as nb
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.special import erfcx

from .baseline import BaselineResult, trace_baselines


@dataclass
class PeakParams:
    smooth_scans: float = 1.0
    min_snr: float = 5.0
    min_prominence: float = 3.0  # in noise units
    min_points: int = 4  # observed centroids inside the peak
    min_fwhm_scans: float = 1.5
    edge_snr: float = 1.0  # boundaries stop where smoothed signal <= baseline + edge_snr * noise


# ------------------------------------------------------------------ numba core
@nb.njit(cache=True)
def gaussian_smooth(y, sigma):
    n = y.shape[0]
    if sigma <= 0:
        return y.copy()
    h = int(np.ceil(3 * sigma))
    k = np.exp(-0.5 * (np.arange(-h, h + 1) / sigma) ** 2)
    out = np.empty(n)
    for i in range(n):
        s = 0.0
        ws = 0.0
        for j in range(-h, h + 1):
            q = i + j
            if 0 <= q < n:
                s += k[j + h] * y[q]
                ws += k[j + h]
        out[i] = s / ws
    return out


@nb.njit(cache=True)
def persistence(y):
    """Return (peak index, persistence) for every local maximum of ``y``.

    The global maximum gets persistence ``max - min``.
    """
    n = y.shape[0]
    order = np.argsort(-y)
    comp = np.full(n, -1, np.int64)  # union-find parent per index (-1 = not yet added)
    peak_of = np.empty(n, np.int64)
    pers = np.zeros(n)
    is_peak = np.zeros(n, np.bool_)
    for q in range(n):
        i = order[q]
        left = comp[i - 1] if i > 0 else -1
        right = comp[i + 1] if i < n - 1 else -1
        if left < 0 and right < 0:
            comp[i] = i
            peak_of[i] = i
            is_peak[i] = True
            continue
        roots = np.empty(2, np.int64)
        nr = 0
        for nbh in (i - 1, i + 1):
            if 0 <= nbh < n and comp[nbh] >= 0:
                r = nbh
                while comp[r] != r:
                    comp[r] = comp[comp[r]]
                    r = comp[r]
                if nr == 0 or roots[0] != r:
                    roots[nr] = r
                    nr += 1
        if nr == 1:
            comp[i] = roots[0]
        else:
            ra, rb = roots[0], roots[1]
            pa, pb = peak_of[ra], peak_of[rb]
            if y[pa] < y[pb]:
                ra, rb = rb, ra
                pa, pb = pb, pa
            pers[pb] = y[pb] - y[i]  # younger (lower) peak dies at this saddle
            comp[rb] = ra
            comp[i] = ra
    gmax = order[0]
    pers[gmax] = y[gmax] - y[order[n - 1]]
    idx = np.flatnonzero(is_peak)
    return idx, pers[idx]


@nb.njit(cache=True)
def _detect_one(y_raw, observed, base, floor, noise, sigma, min_snr, min_prom, min_points, min_fwhm, edge_snr):
    n = y_raw.shape[0]
    y = np.where(observed, y_raw, np.minimum(0.5 * floor, base))
    ys = gaussian_smooth(y, sigma)
    pk, pr = persistence(ys)
    # candidate filter
    keep = np.zeros(pk.shape[0], np.bool_)
    for q in range(pk.shape[0]):
        i = pk[q]
        if pr[q] >= min_prom * noise and ys[i] - base[i] >= min_snr * noise:
            keep[q] = True
    cand = np.sort(pk[keep])
    m = cand.shape[0]
    out = np.zeros((m, 12))
    k = 0
    for q in range(m):
        i = cand[q]
        # valley bounds to neighbouring candidates
        lo = 0
        if q > 0:
            prev = cand[q - 1]
            lo = prev
            for j in range(prev, i):
                if ys[j] < ys[lo]:
                    lo = j
        hi = n - 1
        if q < m - 1:
            nxt = cand[q + 1]
            hi = i
            for j in range(i, nxt + 1):
                if ys[j] < ys[hi]:
                    hi = j
        # trim to baseline return
        left = i
        while left > lo and ys[left - 1] - base[left - 1] > edge_snr * noise:
            left -= 1
        if left > lo:
            left -= 1
        right = i
        while right < hi and ys[right + 1] - base[right + 1] > edge_snr * noise:
            right += 1
        if right < hi:
            right += 1
        nobs = 0
        for j in range(left, right + 1):
            if observed[j]:
                nobs += 1
        h = ys[i] - base[i]
        # FWHM by linear interpolation on the smoothed trace
        half = base[i] + 0.5 * h
        a = i
        while a > left and ys[a] > half:
            a -= 1
        fl = a + (half - ys[a]) / (ys[a + 1] - ys[a]) if ys[a + 1] != ys[a] and ys[a] <= half else float(a)
        b = i
        while b < right and ys[b] > half:
            b += 1
        fr = b - (half - ys[b]) / (ys[b - 1] - ys[b]) if ys[b - 1] != ys[b] and ys[b] <= half else float(b)
        fwhm = fr - fl
        if nobs < min_points or fwhm < min_fwhm:
            continue
        # parabolic apex refinement
        off = 0.0
        if 0 < i < n - 1:
            den = ys[i - 1] - 2 * ys[i] + ys[i + 1]
            if den < 0:
                off = 0.5 * (ys[i - 1] - ys[i + 1]) / den
        area = 0.0
        for j in range(left, right + 1):
            area += max(y[j] - base[j], 0.0)
        apex_left = i - fl
        apex_right = fr - i
        out[k, 0] = i
        out[k, 1] = i + off
        out[k, 2] = left
        out[k, 3] = right
        out[k, 4] = h
        out[k, 5] = area
        out[k, 6] = fwhm
        out[k, 7] = pr[np.searchsorted(pk, i)]
        out[k, 8] = nobs
        out[k, 9] = base[i]
        out[k, 10] = apex_right / apex_left if apex_left > 0 else 1.0
        out[k, 11] = fl
        k += 1
    return out[:k]


def detect_peaks_1d(y, baseline=None, noise=None, floor=None, params: PeakParams | None = None):
    """Detect peaks in a single dense trace (NaN = censored). Returns a structured dict per peak."""
    p = params or PeakParams()
    y = np.asarray(y, dtype=np.float64)
    observed = np.isfinite(y)
    yy = np.where(observed, y, 0.0)
    base = np.zeros_like(yy) if baseline is None else np.asarray(baseline, dtype=np.float64)
    fl = np.full_like(yy, np.nanmin(y) if observed.any() else 1.0) if floor is None else np.asarray(floor, float)
    nz = float(noise) if noise is not None else float(max(np.nanmin(y) / 3.0, 1e-9)) if observed.any() else 1.0
    res = _detect_one(yy, observed, base, fl, nz, p.smooth_scans, p.min_snr, p.min_prominence, p.min_points,
                      p.min_fwhm_scans, p.edge_snr)
    cols = ["apex", "apex_frac", "start", "end", "height", "area", "fwhm", "prominence", "n_obs", "baseline",
            "asymmetry", "fwhm_left"]
    return pd.DataFrame(res, columns=cols)


@nb.njit(cache=True)
def _detect_all(y_flat, obs_flat, base_flat, floor_flat, offs, noise, sigma, min_snr, min_prom, min_points,
                min_fwhm, edge_snr):
    nt = offs.shape[0] - 1
    chunks = []
    tids = []
    for t in range(nt):
        a, b = offs[t], offs[t + 1]
        r = _detect_one(y_flat[a:b], obs_flat[a:b], base_flat[a:b], floor_flat[a:b], noise[t], sigma, min_snr,
                        min_prom, min_points, min_fwhm, edge_snr)
        if r.shape[0]:
            chunks.append(r)
            tids.append(t)
    total = 0
    for c in chunks:
        total += c.shape[0]
    out = np.empty((total, 13))
    k = 0
    for ci in range(len(chunks)):
        c = chunks[ci]
        for r in range(c.shape[0]):
            out[k, 0] = tids[ci]
            out[k, 1:] = c[r]
            k += 1
    return out


def detect_peaks(traces, baselines: BaselineResult | None = None, params: PeakParams | None = None) -> pd.DataFrame:
    """Detect chromatographic peaks in every trace; returns the per-sample feature table."""
    p = params or PeakParams()
    bl = baselines or trace_baselines(traces)
    run = traces.run
    obs = np.isfinite(bl.intensity)
    y = np.where(obs, bl.intensity, 0.0)
    res = _detect_all(y, obs, bl.baseline, bl.floor, bl.offsets, bl.noise, p.smooth_scans, p.min_snr,
                      p.min_prominence, p.min_points, p.min_fwhm_scans, p.edge_snr)
    cols = ["trace", "apex", "apex_frac", "start", "end", "height", "area", "fwhm", "prominence", "n_obs",
            "baseline", "asymmetry", "fwhm_left"]
    df = pd.DataFrame(res, columns=cols)
    if df.empty:
        return df
    t = df["trace"].to_numpy().astype(np.int64)
    a = bl.offsets[t]
    s_apex = bl.scans[a + df["apex"].to_numpy().astype(np.int64)]
    rt = run.rt
    frac = df["apex_frac"].to_numpy() - df["apex"].to_numpy()
    s_float = s_apex + frac
    df["scan_apex"] = s_apex
    df["rt"] = np.interp(s_float, np.arange(run.n_scans), rt)
    df["scan_start"] = bl.scans[a + df["start"].to_numpy().astype(np.int64)]
    df["scan_end"] = bl.scans[a + df["end"].to_numpy().astype(np.int64)]
    df["rt_start"] = rt[df["scan_start"].to_numpy()]
    df["rt_end"] = rt[df["scan_end"].to_numpy()]
    cyc = run.cycle_time
    df["fwhm_s"] = df["fwhm"] * cyc
    df["area"] = df["area"] * cyc
    df["noise"] = bl.noise[t]
    df["snr"] = df["height"] / df["noise"]
    # m/z: inverse-variance weighted mean of observed centroids within the FWHM core, using the
    # trace detector's calibrated precision model; mz_se is the standard error of that mean
    mz = np.empty(len(df))
    mz_se = np.empty(len(df))
    lo = a + np.floor(df["fwhm_left"].to_numpy()).astype(np.int64)
    width = np.maximum(np.ceil(df["fwhm"].to_numpy()).astype(np.int64), 1)
    starts = df["start"].to_numpy().astype(np.int64)
    ends = df["end"].to_numpy().astype(np.int64)
    pm = traces.precision
    for r in range(len(df)):
        pts = bl.point[lo[r]: lo[r] + width[r] + 1]
        pts = pts[pts >= 0]
        if pts.size == 0:
            pts = bl.point[a[r] + starts[r]: a[r] + ends[r] + 1]
            pts = pts[pts >= 0]
        if pts.size == 0:
            mz[r] = traces.table["mz"].iat[int(t[r])]
            mz_se[r] = traces.table["mz_se"].iat[int(t[r])]
            continue
        m = run.mz[pts]
        sig = pm.sigma_ppm(run.intensity[pts]) * 1e-6 * m
        w = 1.0 / sig**2
        mz[r] = np.sum(m * w) / np.sum(w)
        mz_se[r] = np.sqrt(1.0 / np.sum(w))
    df["mz"] = mz
    df["mz_se"] = mz_se
    df["trace"] = t
    df = df.drop(columns=["apex", "apex_frac", "start", "end", "fwhm_left"])
    df = df.sort_values(["mz", "rt"]).reset_index(drop=True)
    df.index.name = "feature"
    lead = ["mz", "rt", "rt_start", "rt_end", "height", "area", "snr", "fwhm_s", "trace"]
    return df[lead + [c for c in df.columns if c not in lead]]


# ------------------------------------------------------------------ EMG mixture fit
def _emg_and_grad(t, h, mu, sig, tau):
    """Height-parametrised EMG (h = density * scale) and its partial derivatives."""
    x = t - mu
    z = (sig / tau - x / sig) / np.sqrt(2)
    g = np.exp(-0.5 * (x / sig) ** 2)
    ex = erfcx(np.clip(z, -25, None))
    f = g * ex / (2 * tau)
    # derivative of erfcx: d/dz erfcx(z) = 2 z erfcx(z) - 2/sqrt(pi)
    dex = 2 * z * ex - 2 / np.sqrt(np.pi)
    dz_dmu = 1 / (np.sqrt(2) * sig)
    dz_dsig = (1 / tau + x / sig**2) / np.sqrt(2)
    dz_dtau = -sig / (tau**2 * np.sqrt(2))
    df_dmu = g * (x / sig**2 * ex + dex * dz_dmu) / (2 * tau)
    df_dsig = g * (x**2 / sig**3 * ex + dex * dz_dsig) / (2 * tau)
    df_dtau = g * (dex * dz_dtau) / (2 * tau) - f / tau
    return h * f, f, h * df_dmu, h * df_dsig, h * df_dtau


def fit_emg_mixture(t, y, init, baseline=None, weights=None, bounds_scale=3.0):
    """Least-squares fit of a sum of EMGs.

    ``init`` is a list of (height, mu, sigma, tau) initial guesses where
    height is the apex height.  Returns a DataFrame with fitted parameters,
    apex rt, apex height and area per component, plus the fitted curve.
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    base = np.zeros_like(y) if baseline is None else np.asarray(baseline, float)
    w = np.ones_like(y) if weights is None else np.asarray(weights, float)
    ok = np.isfinite(y)
    t, y, base, w = t[ok], y[ok], base[ok], w[ok]
    from .synth.lc import emg as emg_density, emg_mode

    x0, lo, hi = [], [], []
    span = t.max() - t.min() if t.size else 1.0
    for h, mu, sig, tau in init:
        dens_apex = emg_density(emg_mode(mu, sig, tau), mu, sig, tau)
        amp = h / max(float(dens_apex), 1e-12)
        x0 += [amp, mu, sig, tau]
        lo += [0.0, mu - bounds_scale * sig - span * 0.05, sig / bounds_scale, 1e-3]
        hi += [np.inf, mu + bounds_scale * sig + span * 0.05, sig * bounds_scale, max(tau, sig) * 10]
    x0 = np.clip(np.array(x0), lo, hi)

    def resid(x):
        model = base.copy()
        for k in range(len(init)):
            a, mu, sig, tau = x[4 * k: 4 * k + 4]
            model += _emg_and_grad(t, a, mu, sig, tau)[0]
        return w * (model - y)

    def jac(x):
        J = np.empty((t.size, x.size))
        for k in range(len(init)):
            a, mu, sig, tau = x[4 * k: 4 * k + 4]
            _, f, dmu, dsig, dtau = _emg_and_grad(t, a, mu, sig, tau)
            J[:, 4 * k] = f * w
            J[:, 4 * k + 1] = dmu * w
            J[:, 4 * k + 2] = dsig * w
            J[:, 4 * k + 3] = dtau * w
        return J

    sol = least_squares(resid, x0, jac=jac, bounds=(lo, hi), method="trf", x_scale="jac", max_nfev=200)
    rows = []
    for k in range(len(init)):
        a, mu, sig, tau = sol.x[4 * k: 4 * k + 4]
        mode = float(emg_mode(mu, sig, tau))
        rows.append(dict(amplitude=a, mu=mu, sigma=sig, tau=tau, rt=mode,
                         height=a * float(emg_density(mode, mu, sig, tau)), area=a))
    rss = float(np.sum(sol.fun**2))
    out = pd.DataFrame(rows)
    out.attrs.update(rss=rss, n=int(t.size), bic=t.size * np.log(rss / max(t.size, 1) + 1e-300)
                     + 4 * len(init) * np.log(max(t.size, 2)), success=bool(sol.success))
    return out


# ------------------------------------------------------------------ reference detectors
def sspp_peaks(y, baseline, noise, S=range(3, 9), maxdist=4, min_snr=5.0):
    """Port of the original macha scale-space peak picking (``sspp`` + SN filter)."""
    from scipy.signal.windows import hamming

    v0 = np.where(np.isfinite(y), y, 0.0)
    pad = max(S)
    v = np.concatenate([np.zeros(pad), v0, np.zeros(pad)])
    C = np.zeros(v.size)
    O = None
    vorig = v.copy()
    for i, s in enumerate(sorted(S)):
        W = hamming(s)
        W /= W.sum()
        vs = np.convolve(v, W, mode="same")
        v = vs
        P = np.flatnonzero((v[1:-1] > v[:-2]) & (v[1:-1] >= v[2:])) + 1
        if P.size == 0:
            break
        dC = np.zeros(v.size)
        if O is None:
            dC[P] = v[P]
            O = P
        else:
            near = np.abs(O[:, None] - P[None, :]) < maxdist
            ints = np.where(near, vorig[O][:, None], -np.inf)
            best = np.argmax(ints, axis=0)
            hit = near[best, np.arange(P.size)]
            dC[O[best[hit]]] = v[P[hit]] * s**1.4
        C += dC
    C = C[pad:-pad]
    idx = np.flatnonzero(C > 0)
    b = np.asarray(baseline, float)
    snr = (v0[idx] - b[idx]) / noise
    return idx[snr >= min_snr]


def cwt_peaks(y, baseline, noise, widths=np.arange(1, 12), min_snr=5.0):
    from scipy.signal import find_peaks_cwt

    v0 = np.where(np.isfinite(y), y, 0.0)
    idx = np.asarray(find_peaks_cwt(v0, widths), dtype=np.int64)
    if idx.size == 0:
        return idx
    b = np.asarray(baseline, float)
    return idx[(v0[idx] - b[idx]) / noise >= min_snr]
