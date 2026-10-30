"""Co-elution graph of traces: elution-profile correlation between traces that apex together.

In DIA, all fragments of one precursor share its elution profile, so grouping fragment traces by profile
correlation (as DIA-Umpire's pseudo-spectra do) is the baseline any mass-relation method must beat.  The edges can
also be combined with the scan-level relation graph.
"""

from __future__ import annotations

import numba as nb
import numpy as np
import pandas as pd

from .._kernels import csr_from_labels


@nb.njit(cache=True)
def _profiles_corr(order_apex, apex, t_offs, t_scans, t_int, s0, s1, apex_tol, min_overlap):
    n = order_apex.shape[0]
    oi = []
    oj = []
    orr = []
    on = []
    buf = np.zeros(4096)
    for a in range(n):
        i = order_apex[a]
        b = a + 1
        while b < n and apex[order_apex[b]] - apex[i] <= apex_tol:
            j = order_apex[b]
            b += 1
            lo = max(s0[i], s0[j])
            hi = min(s1[i], s1[j])
            if hi - lo + 1 < min_overlap or hi - lo + 1 > 4096:
                continue
            # correlation of sqrt intensities on the overlap (absent = 0)
            L = hi - lo + 1
            for k in range(L):
                buf[k] = 0.0
            for q in range(t_offs[i], t_offs[i + 1]):
                s = t_scans[q]
                if lo <= s <= hi:
                    buf[s - lo] = np.sqrt(t_int[q])
            sx = 0.0
            sy = 0.0
            sxx = 0.0
            syy = 0.0
            sxy = 0.0
            ny = np.zeros(L)
            for q in range(t_offs[j], t_offs[j + 1]):
                s = t_scans[q]
                if lo <= s <= hi:
                    ny[s - lo] = np.sqrt(t_int[q])
            for k in range(L):
                x = buf[k]
                y = ny[k]
                sx += x
                sy += y
                sxx += x * x
                syy += y * y
                sxy += x * y
            vx = sxx - sx * sx / L
            vy = syy - sy * sy / L
            if vx <= 1e-12 or vy <= 1e-12:
                continue
            oi.append(i)
            oj.append(j)
            orr.append((sxy - sx * sy / L) / np.sqrt(vx * vy))
            on.append(L)
    return np.array(oi, np.int64), np.array(oj, np.int64), np.array(orr), np.array(on, np.int64)


def coelution_edges(traces, min_points: int = 5, apex_tol: int = 1, min_overlap: int = 5) -> pd.DataFrame:
    """Profile correlation between every pair of traces (>= ``min_points``) whose apexes are <= ``apex_tol`` scans apart."""
    run = traces.run
    order, offs = csr_from_labels(traces.label, len(traces))
    sc = run.scan_of_point[order].astype(np.int64)
    it = run.intensity[order].astype(np.float64)
    n = len(traces)
    apex = np.full(n, -10**9, np.int64)
    s0 = np.zeros(n, np.int64)
    s1 = np.zeros(n, np.int64)
    counts = np.diff(offs)
    ok = counts >= min_points
    for t in np.flatnonzero(ok):  # apex and span per trace
        a, b = offs[t], offs[t + 1]
        apex[t] = sc[a + int(np.argmax(it[a:b]))]
        s0[t], s1[t] = sc[a], sc[b - 1]
    cand = np.flatnonzero(ok)
    order_apex = cand[np.argsort(apex[cand], kind="stable")]
    i, j, r, nn = _profiles_corr(order_apex, apex, offs, sc, it, s0, s1, apex_tol, min_overlap)
    return pd.DataFrame(dict(i=i, j=j, r=r, n=nn))


def coelution_weights(edges: pd.DataFrame, r0: float = 0.9, clip: float = 6.0) -> np.ndarray:
    """Log-odds-like weights: Fisher z of r against r0, scaled by the overlap length's standard error."""
    r = np.clip(edges["r"].to_numpy(), -0.999, 0.999)
    z = (np.arctanh(r) - np.arctanh(r0)) * np.sqrt(np.maximum(edges["n"].to_numpy() - 3, 1))
    return np.clip(z, -clip, clip)


@nb.njit(cache=True)
def _seed_corr(cand, apex, t_offs, t_scans, t_int, s0, s1, cov_offs, cov_seed, seed_s0, seed_offs, seed_y,
               min_overlap, r_min):
    """Correlation of each candidate trace with every seed profile covering its apex; pairs with r >= r_min."""
    ot = []
    ok_ = []
    orr = []
    for a in range(cand.shape[0]):
        t = cand[a]
        sc = apex[t]
        lo_t, hi_t = s0[t], s1[t]
        for q in range(cov_offs[sc], cov_offs[sc + 1]):
            k = cov_seed[q]
            k0 = seed_s0[k]
            k1 = k0 + (seed_offs[k + 1] - seed_offs[k]) - 1
            lo = max(lo_t, k0)
            hi = min(hi_t, k1)
            L = hi - lo + 1
            if L < min_overlap:
                continue
            x = np.zeros(L)
            for v in range(t_offs[t], t_offs[t + 1]):
                s = t_scans[v]
                if lo <= s <= hi:
                    x[s - lo] = np.sqrt(t_int[v])
            sx = sy = sxx = syy = sxy = 0.0
            for u in range(L):
                y = np.sqrt(seed_y[seed_offs[k] + lo + u - k0])
                sx += x[u]
                sy += y
                sxx += x[u] * x[u]
                syy += y * y
                sxy += x[u] * y
            vx = sxx - sx * sx / L
            vy = syy - sy * sy / L
            if vx <= 1e-12 or vy <= 1e-12:
                continue
            r = (sxy - sx * sy / L) / np.sqrt(vx * vy)
            if r >= r_min:
                ot.append(t)
                ok_.append(k)
                orr.append(r)
    return np.array(ot, np.int64), np.array(ok_, np.int64), np.array(orr)


def assign_to_profiles(traces, seed_s0: np.ndarray, seed_y: list, r_min: float = 0.85, exclude=None,
                       shared: bool = False, min_points: int = 5, min_overlap: int = 5) -> pd.DataFrame:
    """Attach traces to seed elution profiles (``seed_y[k]`` = intensity on scans ``seed_s0[k]`` ..).

    Every trace with >= ``min_points`` (and not in ``exclude``) is correlated with each seed profile covering its
    apex.  Returns (trace, seed, r) rows: the best seed per trace, or with ``shared`` every seed reaching ``r_min``.
    """
    run = traces.run
    n_scans = len(run.rt)
    order, offs = csr_from_labels(traces.label, len(traces))
    sc = run.scan_of_point[order].astype(np.int64)
    it = run.intensity[order].astype(np.float64)
    lens = np.array([len(y) for y in seed_y], np.int64)
    seed_s0 = np.asarray(seed_s0, np.int64)
    seed_offs = np.concatenate([[0], np.cumsum(lens)]).astype(np.int64)
    flat = np.concatenate(seed_y) if len(seed_y) else np.zeros(0)
    cov_scan = np.repeat(seed_s0, lens) + np.arange(lens.sum()) - np.repeat(np.cumsum(lens) - lens, lens)
    cov_seed = np.repeat(np.arange(len(seed_y)), lens)
    keep = (cov_scan >= 0) & (cov_scan < n_scans)
    cov_scan, cov_seed = cov_scan[keep], cov_seed[keep]
    o = np.argsort(cov_scan, kind="stable")
    cov_offs = np.searchsorted(cov_scan[o], np.arange(n_scans + 1)).astype(np.int64)
    counts = np.diff(offs)
    ok = counts >= min_points
    if exclude is not None:
        ok &= ~exclude
    cand = np.flatnonzero(ok)
    apex = np.zeros(len(traces), np.int64)
    s0 = np.zeros(len(traces), np.int64)
    s1 = np.zeros(len(traces), np.int64)
    for t in cand:
        a, b = offs[t], offs[t + 1]
        apex[t] = sc[a + int(np.argmax(it[a:b]))]
        s0[t], s1[t] = sc[a], sc[b - 1]
    t, k, r = _seed_corr(cand.astype(np.int64), apex, offs, sc, it, s0, s1, cov_offs, cov_seed[o].astype(np.int64),
                         seed_s0, seed_offs, flat, min_overlap, r_min)
    df = pd.DataFrame(dict(trace=t, seed=k, r=r))
    if not shared and len(df):
        df = df.sort_values("r", ascending=False).drop_duplicates("trace")
    return df.reset_index(drop=True)


def merge_seeds(traces, members: pd.DataFrame, r0: float = 0.9, apex_tol: int = 1, min_overlap: int = 5):
    """Join seeds of one precursor (e.g. its b and y ladders) whose summed chromatograms co-elute and correlate.

    Seeds are treated as traces of their summed signal; pairs apexing within ``apex_tol`` scans get a Fisher-z
    weight against ``r0`` and are partitioned by GAEC.  Returns ``members`` with merged seed ids.
    """
    from .graph import gaec

    run = traces.run
    order, offs = csr_from_labels(traces.label, len(traces))
    sc = run.scan_of_point[order].astype(np.int64)
    it = run.intensity[order].astype(np.float64)
    seeds, k_of = np.unique(members["seed"].to_numpy(), return_inverse=True)
    # summed signal per (seed, scan)
    rows_k, rows_s, rows_i = [], [], []
    for k, r in zip(k_of, members.itertuples()):
        a, b = offs[r.trace], offs[r.trace + 1]
        m = (sc[a:b] >= r.s0) & (sc[a:b] <= r.s1)
        rows_k.append(np.full(int(m.sum()), k))
        rows_s.append(sc[a:b][m])
        rows_i.append(it[a:b][m])
    df = pd.DataFrame(dict(k=np.concatenate(rows_k), s=np.concatenate(rows_s), i=np.concatenate(rows_i)))
    agg = df.groupby(["k", "s"], sort=True)["i"].sum().reset_index()
    k_arr = agg["k"].to_numpy()
    s_offs = np.searchsorted(k_arr, np.arange(len(seeds) + 1))
    s_scans = agg["s"].to_numpy().astype(np.int64)
    s_int = agg["i"].to_numpy()
    apex = np.zeros(len(seeds), np.int64)
    s0 = np.zeros(len(seeds), np.int64)
    s1 = np.zeros(len(seeds), np.int64)
    for k in range(len(seeds)):
        a, b = s_offs[k], s_offs[k + 1]
        apex[k] = s_scans[a + int(np.argmax(s_int[a:b]))]
        s0[k], s1[k] = s_scans[a], s_scans[b - 1]
    order_apex = np.argsort(apex, kind="stable").astype(np.int64)
    i, j, r, nn = _profiles_corr(order_apex, apex, s_offs, s_scans, s_int, s0, s1, apex_tol, min_overlap)
    e = pd.DataFrame(dict(i=i, j=j, r=r, n=nn))
    lab = gaec(len(seeds), i, j, coelution_weights(e, r0)) if len(e) else np.arange(len(seeds))
    return members.assign(seed=lab[k_of])


def assign_to_seeds(traces, members: pd.DataFrame, r_min: float = 0.85, min_points: int = 5,
                    min_overlap: int = 5, shared: bool = False):
    """Grow seed groups (e.g. ladder-connected fragments) by elution profile.

    ``members`` lists seed member blocks (columns seed, trace, s0, s1).  Every other trace with >= ``min_points``
    joins the seed whose summed chromatogram covers its apex and correlates best with its profile, if that
    correlation reaches ``r_min``.  Returns (group of every trace, -1 if none; correlation of the assignment).
    With ``shared`` a trace joins every seed it follows (r >= ``r_min``), and the return value is instead a list
    of member-trace arrays, one per seed.
    """
    run = traces.run
    order, offs = csr_from_labels(traces.label, len(traces))
    sc = run.scan_of_point[order].astype(np.int64)
    it = run.intensity[order].astype(np.float64)
    group = np.full(len(traces), -1, np.int64)
    r_out = np.full(len(traces), np.nan)
    if not len(members):
        return group, r_out
    seeds, k_of = np.unique(members["seed"].to_numpy(), return_inverse=True)
    seed_s0 = np.full(len(seeds), len(run.rt), np.int64)
    seed_s1 = np.zeros(len(seeds), np.int64)
    np.minimum.at(seed_s0, k_of, members["s0"].to_numpy())
    np.maximum.at(seed_s1, k_of, members["s1"].to_numpy())
    seed_y = [np.zeros(int(b - a + 1)) for a, b in zip(seed_s0, seed_s1)]
    for k, r in zip(k_of, members.itertuples()):
        a, b = offs[r.trace], offs[r.trace + 1]
        m = (sc[a:b] >= r.s0) & (sc[a:b] <= r.s1)
        np.add.at(seed_y[k], sc[a:b][m] - seed_s0[k], it[a:b][m])
        group[r.trace] = k
    df = assign_to_profiles(traces, seed_s0, seed_y, r_min, exclude=group >= 0, min_points=min_points,
                            min_overlap=min_overlap, shared=shared)
    if shared:
        own = pd.DataFrame(dict(trace=members["trace"].to_numpy(), seed=k_of))
        allm = pd.concat([own, df[["trace", "seed"]]]).drop_duplicates()
        return [np.sort(d.trace.to_numpy()) for _, d in allm.groupby("seed")]
    group[df.trace.to_numpy()] = df.seed.to_numpy()
    r_out[df.trace.to_numpy()] = df.r.to_numpy()
    return group, r_out
