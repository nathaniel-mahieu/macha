"""Trace (ROI) detection: the public API.

A trace is an m/z-continuous series of centroids, at most one per scan.
Unlike an m/z slice (fixed window) or a greedy ROI builder, traces are found
by multi-target tracking (see :mod:`macha.traces.tracking`) whose m/z
tolerance adapts to each centroid's intensity and to how well the trace's
m/z is already known.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numba as nb
import numpy as np
import pandas as pd
from scipy.optimize import minimize

from .._kernels import csr_from_labels
from ..run import MSRun
from .tracking import kf_loglik, solve_links, trace_stats, track, track_kf


@dataclass
class PrecisionModel:
    """Centroid m/z precision: ``sigma_ppm(I) = sqrt(a^2 + (b * I^-c)^2)``."""

    a: float = 1.0
    b: float = 100.0
    c: float = 0.5

    def sigma_ppm(self, intensity):
        i = np.maximum(np.asarray(intensity, dtype=np.float64), 1e-9)
        return np.sqrt(self.a**2 + (self.b * i ** (-self.c)) ** 2)

    @classmethod
    def initial(cls, run: MSRun, floor_ppm: float = 1.0, low_ppm: float = 5.0) -> "PrecisionModel":
        i5 = float(np.quantile(run.intensity, 0.05)) if run.n_points else 1.0
        return cls(a=floor_ppm, b=low_ppm * i5**0.5, c=0.5)

    @classmethod
    def fit_from_traces(cls, run: MSRun, label: np.ndarray, start: "PrecisionModel | None" = None,
                        max_pairs: int = 400_000, seed: int = 0) -> "PrecisionModel":
        """Maximum-likelihood fit on differences of consecutive centroids in traces.

        Differences cancel slow calibration drift, so the floor ``a`` reflects
        scan-to-scan precision only.  Residuals beyond 5 sigma are trimmed.
        """
        order, offs = csr_from_labels(label)
        scan = run.scan_of_point
        i1, i2 = _consecutive_pairs(order, offs, scan)
        if i1.size < 50:
            return start or cls()
        if i1.size > max_pairs:
            sel = np.random.default_rng(seed).choice(i1.size, max_pairs, replace=False)
            i1, i2 = i1[sel], i2[sel]
        mz1, mz2 = run.mz[i1], run.mz[i2]
        d = (mz2 - mz1) / (0.5 * (mz1 + mz2)) * 1e6
        I1 = run.intensity[i1].astype(np.float64)
        I2 = run.intensity[i2].astype(np.float64)
        m0 = start or cls()
        x0 = np.array([np.log(m0.a), np.log(m0.b), m0.c])
        keep = np.ones(d.size, bool)

        def nll(x):
            a, b, c = np.exp(x[0]), np.exp(x[1]), x[2]
            v = 2 * a * a + (b * I1[keep] ** -c) ** 2 + (b * I2[keep] ** -c) ** 2
            return 0.5 * np.sum(np.log(v) + d[keep] ** 2 / v)

        for _ in range(3):
            res = minimize(nll, x0, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-3, maxiter=2000))
            x0 = res.x
            a, b, c = np.exp(x0[0]), np.exp(x0[1]), x0[2]
            v = 2 * a * a + (b * I1**-c) ** 2 + (b * I2**-c) ** 2
            keep = d**2 < 25 * v
        return cls(a=float(np.exp(x0[0])), b=float(np.exp(x0[1])), c=float(x0[2]))


@dataclass
class TrackParams:
    max_gap: int = 2  # missing scans bridged while tracking
    gate: float = 4.0  # m/z gate in combined standard deviations
    max_ppm: float = 25.0  # absolute cap of the search window
    birth_cost: float = 8.0  # NLL of starting a new trace instead of linking
    gap_cost: float = 0.75  # NLL per missing scan
    intensity_sigma: float = 0.9  # sd of log intensity ratio between consecutive scans
    intensity_cap: float = 3.0
    motion: str = "velocity"  # "velocity": constant-velocity Kalman per trace; "level": forgetting mean
    q_ppm: float | None = None  # velocity model: acceleration noise (ppm per scan^1.5); learned if None
    v0_ppm: float = 0.2  # velocity model: prior sd of a new trace's drift, ppm per scan
    q_init: float = 0.05  # velocity model: generous first-pass value before calibration
    forget: float = 0.97  # level model: per-scan decay of a trace's m/z memory
    drift_ppm: float = 0.15  # level model: per-scan random-walk allowance for a trace's m/z
    close_gaps: int = 8  # second pass: rejoin trace ends/starts up to this many missing scans
    min_points: int = 5
    calibrate: bool = True
    precision: PrecisionModel | None = None


@nb.njit(cache=True)
def _consecutive_pairs(order, offs, scan):
    n = 0
    for t in range(offs.shape[0] - 1):
        n += max(offs[t + 1] - offs[t] - 1, 0)
    i1 = np.empty(n, np.int64)
    i2 = np.empty(n, np.int64)
    k = 0
    for t in range(offs.shape[0] - 1):
        for q in range(offs[t], offs[t + 1] - 1):
            a, b = order[q], order[q + 1]
            if scan[b] - scan[a] == 1:
                i1[k] = a
                i2[k] = b
                k += 1
    return i1[:k], i2[:k]


@nb.njit(cache=True)
def _gap_candidates(mean, se, s0, s1, first_int, last_int, max_gap_lo, max_gap_hi, gate, drift_ppm,
                    birth_cost, gap_cost, int_sigma, int_cap):
    nt = mean.shape[0]
    o = np.argsort(mean)
    ca = []
    cb = []
    cg = []
    for ii in range(nt):
        A = o[ii]
        tol = mean[A] * 30e-6
        for jj in range(ii + 1, nt):
            B = o[jj]
            if mean[B] - mean[A] > tol:
                break
            for (E, S_) in ((A, B), (B, A)):
                gap = s0[S_] - s1[E] - 1
                if gap < max_gap_lo or gap > max_gap_hi:
                    continue
                v = se[E] ** 2 + se[S_] ** 2 + (drift_ppm * 1e-6 * mean[E]) ** 2 * (gap + 1)
                d2 = (mean[E] - mean[S_]) ** 2 / v
                if d2 > gate * gate:
                    continue
                r = np.log(first_int[S_] / last_int[E])
                ci = min(0.5 * r * r / (int_sigma * int_sigma * (gap + 1)), int_cap)
                cost = 0.5 * d2 + ci + gap_cost * gap
                if cost < birth_cost:
                    ca.append(E)
                    cb.append(S_)
                    cg.append(birth_cost - cost)
    return np.array(ca, np.int64), np.array(cb, np.int64), np.array(cg, np.float64)


@dataclass
class Traces:
    run: MSRun
    label: np.ndarray
    order: np.ndarray
    offsets: np.ndarray
    table: pd.DataFrame
    precision: PrecisionModel
    params: TrackParams
    timings: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.table)

    def points(self, i: int) -> np.ndarray:
        return self.order[self.offsets[i] : self.offsets[i + 1]]

    def eic(self, i: int, pad: int = 0):
        """Dense trace: (scans, rt, intensity, mz) with NaN where no centroid was reported."""
        idx = self.points(i)
        sc = self.run.scan_of_point[idx]
        s0 = max(int(sc[0]) - pad, 0)
        s1 = min(int(sc[-1]) + pad, self.run.n_scans - 1)
        scans = np.arange(s0, s1 + 1)
        inten = np.full(scans.size, np.nan)
        mz = np.full(scans.size, np.nan)
        inten[sc - s0] = self.run.intensity[idx]
        mz[sc - s0] = self.run.mz[idx]
        return scans, self.run.rt[scans], inten, mz

    def filter(self, min_points: int) -> "Traces":
        keep = self.table["n_points"].to_numpy() >= min_points
        return _finalise(self.run, self.label, self.precision, self.params, keep_mask=keep,
                         timings=self.timings)


def _finalise(run, label, precision, params, keep_mask=None, timings=None) -> Traces:
    """Relabel densely (sorted by m/z) and compute the trace table."""
    label = _compact(label)
    n_raw = int(label.max()) + 1 if label.size else 0
    counts = np.bincount(label[label >= 0], minlength=n_raw)
    keep = counts >= params.min_points if keep_mask is None else keep_mask
    order, offs = csr_from_labels(label, n_raw)
    scan = run.scan_of_point
    pm = precision
    stats = trace_stats(label, order, offs, scan, run.mz, run.intensity, pm.a, pm.b, pm.c) if n_raw else None
    if n_raw == 0:
        empty = pd.DataFrame(columns=["mz", "mz_se", "mz_sd_ppm", "scan_start", "scan_end", "n_points",
                                      "apex_scan", "apex_intensity", "total_intensity", "rt_start",
                                      "rt_end", "apex_rt"])
        return Traces(run, label, np.empty(0, np.int64), np.zeros(1, np.int64), empty, pm, params,
                      timings or {})
    mean = stats[0]
    keep_ids = np.flatnonzero(keep)
    keep_ids = keep_ids[np.argsort(mean[keep_ids], kind="stable")]
    remap = np.full(n_raw, -1, np.int64)
    remap[keep_ids] = np.arange(keep_ids.size)
    new_label = np.where(label >= 0, remap[np.maximum(label, 0)], -1)
    new_order, new_offs = csr_from_labels(new_label, keep_ids.size)
    # sort each trace's points by scan (they are, since labels were assigned scan-ordered,
    # but gap-closing unions can interleave): stable argsort by (label, scan)
    o2 = np.lexsort((scan[new_order], new_label[new_order]))
    new_order = new_order[o2]
    cols = ["mz", "mz_se", "mz_sd_ppm", "scan_start", "scan_end", "n_points", "apex_scan",
            "apex_intensity", "total_intensity", "first_intensity", "last_intensity"]
    table = pd.DataFrame({c: v[keep_ids] for c, v in zip(cols, stats)})
    table = table.drop(columns=["first_intensity", "last_intensity"])
    table["rt_start"] = run.rt[table["scan_start"].to_numpy()]
    table["rt_end"] = run.rt[table["scan_end"].to_numpy()]
    table["apex_rt"] = run.rt[table["apex_scan"].to_numpy()]
    table.index.name = "trace"
    return Traces(run, new_label, new_order, new_offs, table, pm, params, timings or {})


def _compact(label) -> np.ndarray:
    label = np.asarray(label, dtype=np.int64)
    out = np.full(label.shape, -1, np.int64)
    ok = label >= 0
    if ok.any():
        _, out[ok] = np.unique(label[ok], return_inverse=True)
    return out


def _track_once(run: MSRun, p: TrackParams, pm: PrecisionModel, q: float | None = None) -> np.ndarray:
    it = run.intensity.astype(np.float64)
    if p.motion == "velocity":
        return track_kf(run.mz, it, run.offsets, pm.a, pm.b, pm.c, p.max_gap, p.gate, p.max_ppm, p.birth_cost,
                        p.gap_cost, p.intensity_sigma, p.intensity_cap, p.q_init if q is None else q, p.v0_ppm)
    return track(run.mz, it, run.offsets, pm.a, pm.b, pm.c, p.max_gap, p.gate, p.max_ppm, p.birth_cost,
                 p.gap_cost, p.intensity_sigma, p.intensity_cap, p.forget, p.drift_ppm)


def fit_motion(run: MSRun, label: np.ndarray, pm: PrecisionModel, v0_ppm: float = 0.2,
               grid=(0.0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4)) -> tuple[float, list]:
    """Acceleration noise maximising the (capped) innovation likelihood of the given traces."""
    order, offs = csr_from_labels(label)
    it = run.intensity.astype(np.float64)
    lls = [kf_loglik(order, offs, run.scan_of_point, run.mz, it, pm.a, pm.b, pm.c, q, v0_ppm, 12.5)[0]
           for q in grid]
    return float(grid[int(np.argmax(lls))]), lls


def close_gaps(run: MSRun, label: np.ndarray, p: TrackParams, pm: PrecisionModel) -> np.ndarray:
    """Rejoin trace segments separated by short dropouts (an LAP between ends and starts)."""
    label = np.asarray(label, np.int64)
    n_raw = int(label.max()) + 1
    order, offs = csr_from_labels(label, n_raw)
    st = trace_stats(label, order, offs, run.scan_of_point, run.mz, run.intensity, pm.a, pm.b, pm.c)
    mean, se, _, s0, s1, npt, _, _, _, fi, li = st
    ca, cb, cg = _gap_candidates(mean, se, s0, s1, fi, li, 1, p.close_gaps, p.gate, p.drift_ppm,
                                 p.birth_cost, p.gap_cost * 0.5, p.intensity_sigma, p.intensity_cap)
    match = solve_links(ca, cb, cg, n_raw, n_raw)
    # union chains end -> start (acyclic: a start always follows its end in time)
    root = np.arange(n_raw)
    has = match >= 0
    root[match[has]] = np.flatnonzero(has)
    while True:  # pointer jumping to chain heads
        nxt = root[root]
        if np.array_equal(nxt, root):
            break
        root = nxt
    return np.where(label >= 0, root[np.maximum(label, 0)], -1)


def detect_traces(run: MSRun, params: TrackParams | None = None, **overrides) -> Traces:
    """Detect m/z traces in a centroided run.

    Steps: (1) track with an initial precision model; (2) re-fit the precision
    model on the first-pass traces (``calibrate``); (3) track again; (4) close
    gaps; (5) drop traces shorter than ``min_points`` and relabel by m/z.
    """
    import time

    p = params or TrackParams()
    for k, v in overrides.items():
        setattr(p, k, v)
    timings = {}
    t0 = time.perf_counter()
    pm = p.precision or PrecisionModel.initial(run)
    q = p.q_ppm
    label = _track_once(run, p, pm, q)
    timings["track_pass1"] = time.perf_counter() - t0
    if p.calibrate and (p.precision is None or (p.motion == "velocity" and p.q_ppm is None)):
        t1 = time.perf_counter()
        counts = np.bincount(label)
        good = np.where(counts[label] >= 8, label, -1)
        if p.precision is None:
            pm = PrecisionModel.fit_from_traces(run, good, start=pm)
        if p.motion == "velocity" and p.q_ppm is None:
            q, _ = fit_motion(run, good, pm, p.v0_ppm)
        label = _track_once(run, p, pm, q)
        timings["calibrate_track_pass2"] = time.perf_counter() - t1
    if p.motion == "velocity":
        p = replace(p, q_ppm=q if q is not None else p.q_init)
    if p.close_gaps > 0 and label.size:
        t2 = time.perf_counter()
        label = close_gaps(run, label, p, pm)
        timings["close_gaps"] = time.perf_counter() - t2
    t3 = time.perf_counter()
    out = _finalise(run, label, pm, p, timings=timings)
    timings["finalise"] = time.perf_counter() - t3
    timings["total"] = time.perf_counter() - t0
    return out
