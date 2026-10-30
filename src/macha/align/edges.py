"""Alignment anchors from scan-graph relations: drift fitting without peak detection.

A single ion is ambiguous across runs wherever m/z is crowded, so feature anchors must be the only feature within
a few ppm over a wide RT window, and they need picked peaks.  A confident relation of the scan graph is keyed by a
rule and the m/z of both of its ions, which makes it far more specific, and it exists before any peak is picked:

1. **Relations.**  Confident trace-to-trace segments of one run are pooled per (anchor, partner, rule).  The
   relation says which two ions and roughly when; its centre comes from both ions' full traces in a fixed window
   of +/- 2 elution sd (the run's own, learned by the scan graph) around it, weighted by intensity squared and
   re-centred.  The scans where both ions were seen are not used directly: a weak partner is found only near the
   apex, so that window is truncated by an amount that depends on abundance, and so is its centre.  A graph block
   is not used either: its extent depends on what else the run happened to find.
2. **Precision.**  A relation's SNR is its height over the noise level in its two ions' own m/z windows
   (:mod:`macha.traces.noise`).  The centre's standard error is ``sqrt((k1 sd / snr)^2 + (k2 sd)^2)`` (``sd`` the
   elution spread); a first estimate of ``k1``, ``k2`` comes from relations of one compound that share no ion.
3. **Specificity.**  A relation anchors only if it elutes (most of its ions' signal near its centre, not a
   chemical-background plateau), has no competitor (no other ion within the matching tolerance reaches a fifth of
   its height in the matching window, counting the ion's own trace beyond a valley, so an isomer joined into one
   trace still counts), is unlikely by chance in its ions' noise-floor windows, and is unique in its run.  The
   matching window comes from the drift measured in a first pass (:func:`drift_window`).
4. **Matching.**  Relations are pooled over runs and chain-clustered by rule, then anchor m/z, partner m/z and RT;
   a cluster with at most one relation per run, present in ``min_frac`` of the runs, is an anchor.
5. **Calibration.**  Within a run, relations of different intensity are truncated differently, so their centre
   differences overstate the error that matters across runs.  After matching, the precision is refitted on double
   differences ``(t_a - t_b)_j - (t_a - t_b)_j'`` of two anchors of one compound with no shared ion: drift, the
   compound's RT error and each relation's own bias cancel.
6. **Duplicates.**  Relations that share an ion measure the same elution.  ``dedupe`` keeps one anchor per group of
   anchors linked through a shared trace; it is off by default, because with K and the smoothing chosen by
   cross-validation (:func:`macha.align.drift.cv_drift`) duplicates do not inflate the number of classes, while
   thinning leaves classes too sparse in places.  Averaging a group was worse still: one wrong link merges two
   compounds' centres.

:func:`relation_groups` groups *every* relation across runs in reference time (after the drift is known), with a
per-run abundance, which is cross-run grouping at the scan-graph stage before any peak detection.
"""

from __future__ import annotations

from dataclasses import dataclass

import numba as nb
import numpy as np
import pandas as pd

from ..traces.noise import floor_halfwidth_ppm, window_noise
from .anchors import Anchors

# apex_se in units of the elution sd: 1.5 fwhm / (2.355 snr) and 0.02 fwhm (fallback without enough pairs)
DEFAULT_PRECISION = (1.5, 0.02 * 2.355)


@nb.njit(cache=True)
def _apex_in(order, offs, scan_of, intensity, t, s0, s1):
    """Scan of the most intense point of trace ``t`` within scans s0..s1 (-1 if none)."""
    best = -1.0
    out = -1
    for k in range(offs[t], offs[t + 1]):
        sc = scan_of[order[k]]
        if s0 <= sc <= s1 and intensity[order[k]] > best:
            best = intensity[order[k]]
            out = sc
    return out


@nb.njit(cache=True)
def _window_centres(order, offs, scan_of, rt, intensity, ta, tb, c0, h, power, n_iter):
    """Two traces' points within +/- h of a centre (re-centred ``n_iter`` times with ``intensity^power`` weights).

    Columns: centre, sd, summed weight, max intensity of trace a, max intensity of trace b (in the final window).
    """
    n = ta.shape[0]
    out = np.full((n, 5), np.nan)
    for r in range(n):
        c = c0[r]
        for _ in range(n_iter):
            sw = 0.0
            swt = 0.0
            swtt = 0.0
            mx = np.zeros(2)
            for e in range(2):
                t = ta[r] if e == 0 else tb[r]
                lo, hi = offs[t], offs[t + 1]
                end = hi
                while lo < hi:  # first point with rt >= c - h (a trace's points are in scan order)
                    mid = (lo + hi) // 2
                    if rt[scan_of[order[mid]]] < c - h:
                        lo = mid + 1
                    else:
                        hi = mid
                k = lo
                while k < end and rt[scan_of[order[k]]] <= c + h:
                    x = rt[scan_of[order[k]]]
                    v = intensity[order[k]]
                    w = v ** power
                    sw += w
                    swt += w * x
                    swtt += w * x * x
                    if v > mx[e]:
                        mx[e] = v
                    k += 1
            if sw <= 0:
                break
            c = swt / sw
            out[r, 0] = c
            out[r, 1] = np.sqrt(max(swtt / sw - c * c, 0.0))
            out[r, 2] = sw
            out[r, 3] = mx[0]
            out[r, 4] = mx[1]
    return out


def relation_table(g, traces, p_min: float = 0.9, pad: int = 20, gate: float = 4.0, sys_ppm: float = 0.3,
                   precision: tuple | None = None, window: float = 2.0, power: float = 2.0) -> pd.DataFrame:
    """Confident trace-to-trace relations of one run with their own elution centre (one row per relation).

    ``g`` needs ``segments``, ``blocks`` (with ``compound``) and ``prior`` of a :class:`ScanGraph`.  Besides the
    centre: ``area`` (both ions' summed signal in the window), ``contain`` (share of the ions' signal within
    +/- ``window`` sd of all signal within three times that: ~1 for an eluting peak, ~1/3 for a plateau) and the
    heights of both ions.  ``precision`` = (k1, k2) fixes the centre standard error model; by default it is learned
    from the run (:func:`fit_centre_precision`, stored in ``attrs``).
    """
    seg = g.segments
    s = seg[(seg["posterior"].to_numpy() >= p_min) & (seg["partner"].to_numpy() >= 0)].reset_index(drop=True)
    cols = ["anchor", "partner", "target", "kind", "z", "steps"]
    if not len(s):
        return pd.DataFrame(columns=cols + ["mz_a", "mz_b", "rt", "rt_sd", "rt_se", "area", "contain", "height_a",
                                            "height_b", "height", "snr", "n", "s0", "s1", "posterior", "log_ratio",
                                            "compound"])
    run = traces.run
    rel = s.groupby(cols, sort=False).ngroup().to_numpy()
    n_rel = int(rel.max()) + 1
    first = pd.Series(np.arange(len(s))).groupby(rel).first().to_numpy()
    out = s.loc[first, cols].reset_index(drop=True)
    out["n"] = np.bincount(rel, s["n"].to_numpy(), n_rel).astype(np.int64)
    out["s0"] = s.groupby(rel)["s0"].min().to_numpy()
    out["s1"] = s.groupby(rel)["s1"].max().to_numpy()
    tab = traces.table
    out["mz_a"] = tab["mz"].to_numpy()[out["anchor"].to_numpy()]
    out["mz_b"] = tab["mz"].to_numpy()[out["partner"].to_numpy()]
    sigma = g.prior.sigma if g.prior is not None else 2.0 * float(run.cycle_time)
    h = window * sigma
    order, offs = traces.order, traces.offsets
    sc_of = run.scan_of_point.astype(np.int64)
    inten = run.intensity.astype(np.float64)
    ta, tb = out["anchor"].to_numpy(np.int64), out["partner"].to_numpy(np.int64)
    # seed: the anchor's most intense point within the relation's scans
    seed = np.array([_apex_in(order, offs, sc_of, inten, a, s0, s1)
                     for a, s0, s1 in zip(ta, out["s0"].to_numpy(), out["s1"].to_numpy())])
    c0 = np.where(seed >= 0, run.rt[np.maximum(seed, 0)], 0.5 * (run.rt[out["s0"].to_numpy()] + run.rt[out["s1"].to_numpy()]))
    W = _window_centres(order, offs, sc_of, run.rt, inten, ta, tb, c0, h, power, 3)
    out["rt"] = np.where(np.isfinite(W[:, 0]), W[:, 0], c0)
    out["rt_sd"] = np.maximum(np.nan_to_num(W[:, 1], nan=sigma), float(run.cycle_time) / np.sqrt(12))
    # abundance and containment from plain sums in the window and in three times the window
    Q = _window_centres(order, offs, sc_of, run.rt, inten, ta, tb, out["rt"].to_numpy(), h, 1.0, 1)
    Q3 = _window_centres(order, offs, sc_of, run.rt, inten, ta, tb, out["rt"].to_numpy(), 3 * h, 1.0, 1)
    out["area"] = np.nan_to_num(Q[:, 2])
    out["contain"] = np.where(Q3[:, 2] > 0, np.nan_to_num(Q[:, 2]) / Q3[:, 2], 0.0)
    out["height_a"] = np.nan_to_num(Q[:, 3])
    out["height_b"] = np.nan_to_num(Q[:, 4])
    height = out["height_a"].to_numpy() + out["height_b"].to_numpy()
    out["height"] = height
    # noise level of the summed signal from both ions' own m/z windows (run floor where a window is empty)
    hw = floor_halfwidth_ppm(traces, gate, sys_ppm)
    floor = float(np.log(np.quantile(run.intensity, 0.05)))
    lv = []
    for mzc, tr in (("mz_a", "anchor"), ("mz_b", "partner")):
        wn = window_noise(traces, out[mzc].to_numpy(), out["s0"].to_numpy() - pad, out["s1"].to_numpy() + pad, hw,
                          out[tr].to_numpy())
        lv.append(np.exp(wn["noise_log_mean"].fillna(floor).to_numpy()))
    out["snr"] = height / np.sqrt(lv[0] ** 2 + lv[1] ** 2)
    out["posterior"] = s.groupby(rel)["posterior"].max().to_numpy()
    out["log_ratio"] = (s["ratio_med"] * s["n"]).groupby(rel).sum().to_numpy() / out["n"].to_numpy()
    node = s["node_a"].to_numpy()
    comp = np.where(node >= 0, g.blocks["compound"].to_numpy()[np.maximum(node, 0)], -1)
    out["compound"] = pd.Series(comp).groupby(rel).first().to_numpy()
    k1, k2, n_pairs = precision + (0,) if precision is not None else fit_centre_precision(out)
    out["rt_se"] = centre_se(out, k1, k2)
    out.attrs["precision"] = dict(k1=k1, k2=k2, n_pairs=n_pairs)
    return out


def centre_se(rel: pd.DataFrame, k1: float, k2: float) -> np.ndarray:
    sd = rel["rt_sd"].to_numpy()
    return np.sqrt((k1 * sd / np.maximum(rel["snr"].to_numpy(), 1.0)) ** 2 + (k2 * sd) ** 2)


def _fit_precision(d, u2, v2, max_n: int = 20000, seed: int = 0):
    """ML (k1, k2) for differences ``d ~ (1 - eps) N(0, k1^2 u2 + k2^2 v2) + eps U`` on a grid (EM for eps)."""
    if d.size > max_n:
        sel = np.random.default_rng(seed).choice(d.size, max_n, replace=False)
        d, u2, v2 = d[sel], u2[sel], v2[sel]
    f0 = 1.0 / (2 * max(float(np.abs(d).max()), 10.0))
    best = (-np.inf, *DEFAULT_PRECISION)
    for k1 in np.logspace(-1.5, 1.5, 49):
        for k2 in np.concatenate([[0.0], np.logspace(-3.0, -0.7, 30)]):
            var = np.maximum(k1 * k1 * u2 + k2 * k2 * v2, 1e-8)
            f1 = np.exp(-0.5 * d * d / var) / np.sqrt(2 * np.pi * var)
            eps = 0.05
            for _ in range(15):
                rr = eps * f0 / (eps * f0 + (1 - eps) * f1)
                eps = float(np.clip(rr.mean(), 1e-3, 0.5))
            ll = float(np.log(eps * f0 + (1 - eps) * f1).sum())
            if ll > best[0]:
                best = (ll, float(k1), float(k2))
    return best[1], best[2]


def fit_centre_precision(rel: pd.DataFrame, min_snr: float = 3.0, per_compound: int = 12, seed: int = 0):
    """Learn (k1, k2) of the centre standard error from relations of one compound that share no ion.

    Their centre difference is ``N(0, se_1^2 + se_2^2)`` (same elution, independent noise) plus a uniform outlier
    component for impure compounds.  Relations of different intensity are truncated differently, so this
    within-run estimate includes relation-specific bias that cancels across runs; :func:`relation_anchors`
    recalibrates on double differences.  Returns (k1, k2, pairs).
    """
    r = rel[(rel["compound"].to_numpy() >= 0) & (rel["snr"].to_numpy() >= min_snr)]
    rng = np.random.default_rng(seed)
    pi, pj = [], []
    for _, d in r.groupby("compound"):
        idx = d.index.to_numpy()
        if idx.size > per_compound:
            idx = rng.choice(idx, per_compound, replace=False)
        a, b = rel.loc[idx, "anchor"].to_numpy(), rel.loc[idx, "partner"].to_numpy()
        for x in range(idx.size):
            for y in range(x + 1, idx.size):
                if len({a[x], b[x], a[y], b[y]}) == 4:
                    pi.append(idx[x])
                    pj.append(idx[y])
    if len(pi) < 30:
        return DEFAULT_PRECISION + (len(pi),)
    pi, pj = np.array(pi), np.array(pj)
    d = rel["rt"].to_numpy()[pi] - rel["rt"].to_numpy()[pj]
    sd, snr = rel["rt_sd"].to_numpy(), np.maximum(rel["snr"].to_numpy(), 1.0)
    u2 = (sd[pi] / snr[pi]) ** 2 + (sd[pj] / snr[pj]) ** 2
    v2 = sd[pi] ** 2 + sd[pj] ** 2
    return (*_fit_precision(d, u2, v2), int(pi.size))


def _double_differences(P: pd.DataFrame, per_compound: int = 12, seed: int = 0):
    """(t_a - t_b in run j) - (t_a - t_b in run j') for anchor clusters a, b in one scan-graph compound with no
    shared ion in both runs: drift, the compound's RT error and each relation's own bias cancel, leaving the
    measurement noise of four centres.  Returns (d, u2, v2)."""
    rng = np.random.default_rng(seed)
    sd = P["rt_sd"].to_numpy()
    u = (sd / np.maximum(P["snr"].to_numpy(), 1.0)) ** 2
    v = sd ** 2
    by_pair: dict = {}
    for (j, _), d in P[P["compound"] >= 0].groupby(["run", "compound"]):
        idx = d.index.to_numpy()
        if idx.size > per_compound:
            idx = rng.choice(idx, per_compound, replace=False)
        c = P.loc[idx, "cluster"].to_numpy()
        a, b, t = P.loc[idx, "anchor"].to_numpy(), P.loc[idx, "partner"].to_numpy(), P.loc[idx, "rt"].to_numpy()
        pos = P.index.get_indexer(idx)
        for x in range(idx.size):
            for y in range(x + 1, idx.size):
                if c[x] == c[y] or len({a[x], b[x], a[y], b[y]}) < 4:
                    continue
                key, sgn = ((c[x], c[y]), 1.0) if c[x] < c[y] else ((c[y], c[x]), -1.0)
                by_pair.setdefault(key, []).append((sgn * (t[x] - t[y]), u[pos[x]] + u[pos[y]], v[pos[x]] + v[pos[y]]))
    d, u2, v2 = [], [], []
    for obs in by_pair.values():
        for x in range(len(obs)):
            for y in range(x + 1, len(obs)):
                d.append(obs[x][0] - obs[y][0])
                u2.append(obs[x][1] + obs[y][1])
                v2.append(obs[x][2] + obs[y][2])
    return np.array(d), np.array(u2), np.array(v2)


NOISE_WINDOWS = (10.0, 30.0, 60.0)


def relation_noise(rel: pd.DataFrame, traces, windows=NOISE_WINDOWS, min_support: int = 3, gate: float = 4.0,
                   sys_ppm: float = 0.3, match_ppm: float = 6.0) -> pd.DataFrame:
    """Interference in both ions' m/z windows around the relation, for several RT half-windows (s): the chance of a
    spurious relation and the strongest competitor.

    With per-scan hit fractions ``f_a``, ``f_b`` in the two noise-floor windows, a chance relation needs both
    occupied in ``min_support`` consecutive scans, so over ``n`` scans ``p_chance ~ 1 - exp(-n (f_a f_b)^k)``.
    ``competitor`` is the most intense other centroid within ``match_ppm`` (what cross-run matching could confuse
    with the ion, e.g. an isomer) over the RT window, relative to the ion's own height; the larger of the two ions.
    Which window matters depends on the drift, so :func:`relation_anchors` picks the column for its matching
    window (``p_chance_<w>``, ``competitor_<w>``).
    """
    run = traces.run
    hw = floor_halfwidth_ppm(traces, gate, sys_ppm)
    rt = rel["rt"].to_numpy()
    order, offs = traces.order, traces.offsets
    sc_of = run.scan_of_point.astype(np.int64)
    inten = run.intensity.astype(np.float64)
    out = {}
    for wdw in windows:
        s0 = np.searchsorted(run.rt, rt - wdw, "left")
        s1 = np.searchsorted(run.rt, rt + wdw, "right") - 1
        a = window_noise(traces, rel["mz_a"].to_numpy(), s0, s1, hw, rel["anchor"].to_numpy())
        b = window_noise(traces, rel["mz_b"].to_numpy(), s0, s1, hw, rel["partner"].to_numpy())
        q = (a["hit_frac"].to_numpy() * b["hit_frac"].to_numpy()) ** min_support
        # other traces anywhere in the RT window; the ion's own trace only beyond a valley (an isomer joined into the
        # same trace still competes, the ion's own tails do not)
        comp = []
        for mzc, tr, hc in (("mz_a", "anchor", "height_a"), ("mz_b", "partner", "height_b")):
            other = window_noise(traces, rel[mzc].to_numpy(), s0, s1, match_ppm, rel[tr].to_numpy())["noise_max"]
            own = _beyond_valley(order, offs, sc_of, run.rt, inten, rel[tr].to_numpy(np.int64), rt, float(wdw))
            mx = np.fmax(np.nan_to_num(other.to_numpy()), own)
            comp.append(mx / np.maximum(rel[hc].to_numpy(), 1e-9))
        k = f"{wdw:g}"
        out[f"p_chance_{k}"] = -np.expm1(-a["n_scans"].to_numpy() * q)
        out[f"competitor_{k}"] = np.maximum(*comp)
    return pd.DataFrame(out, index=rel.index)


def _noise_column(noise: pd.DataFrame, name: str, rt_window: float) -> np.ndarray:
    """The stored column for the smallest half-window covering ``rt_window`` (the widest if none does)."""
    ws = sorted(float(c.split("_")[-1]) for c in noise.columns if c.startswith(name + "_"))
    w = next((x for x in ws if x >= rt_window), ws[-1])
    return noise[f"{name}_{w:g}"].to_numpy()


def drift_window(relations, noise, q: float = 0.95, margin: float = 2.0, floor: float = 10.0, **kw) -> float:
    """RT half-window for matching relations across runs, from a first pass at the widest window: the largest
    per-run ``q`` quantile of anchor deviations from their cross-run median, times ``margin`` (at least ``floor``
    seconds).  The worst run sets it: one injection can drift far more than the rest.  A fixed wide window is costly
    in dense data: on Astral proteomics MS1 a stronger ion within 6 ppm is almost always found somewhere in
    +/- 60 s, so the competitor filter would reject most relations."""
    A, _, _ = relation_anchors(relations, noise=noise, rt_window=max(NOISE_WINDOWS), recalibrate=False, **kw)
    dev = np.abs(A.rt - np.nanmedian(A.rt, 1, keepdims=True))
    if not np.isfinite(dev).any():
        return float(max(NOISE_WINDOWS))
    per_run = np.nanquantile(np.where(np.isfinite(dev), dev, np.nan), q, axis=0)
    return float(max(margin * np.nanmax(per_run), floor))


@nb.njit(cache=True)
def _beyond_valley(order, offs, scan_of, rt, intensity, traces, centre, R):
    """Most intense point of each trace within +/- R of ``centre`` but beyond the first valley on either side: after
    the signal has fallen below half of the peak, a rise to twice the running minimum starts another elution."""
    n = traces.shape[0]
    out = np.zeros(n)
    for r in range(n):
        t = traces[r]
        lo, hi = offs[t], offs[t + 1]
        if hi <= lo:
            continue
        # point nearest the centre
        a, b = lo, hi
        while a < b:
            mid = (a + b) // 2
            if rt[scan_of[order[mid]]] < centre[r]:
                a = mid + 1
            else:
                b = mid
        k0 = min(max(a, lo), hi - 1)
        # the local peak: climb from k0 to the nearby maximum
        peak = intensity[order[k0]]
        best = 0.0
        for step in (1, -1):
            vmin = peak
            top = peak
            past = False
            k = k0 + step
            while lo <= k < hi and abs(rt[scan_of[order[k]]] - centre[r]) <= R:
                v = intensity[order[k]]
                if not past:
                    if v > top and vmin >= 0.5 * top:
                        top = v  # still climbing the same peak
                        vmin = v
                    elif v < vmin:
                        vmin = v
                    if vmin < 0.5 * top and v >= 2.0 * vmin:
                        past = True
                if past and v > best:
                    best = v
                k += step
        out[r] = best
    return out


@nb.njit(cache=True)
def _ambiguous(key, mz_a, mz_b, rt, tol_ppm, rt_window):
    """True where another relation with the same key and both m/z within tol lies within rt_window (sorted by key,
    then mz_a)."""
    n = key.shape[0]
    out = np.zeros(n, np.bool_)
    for i in range(n):
        tol = mz_a[i] * tol_ppm * 1e-6
        j = i + 1
        while j < n and key[j] == key[i] and mz_a[j] - mz_a[i] <= tol:
            if abs(mz_b[j] - mz_b[i]) <= mz_b[i] * tol_ppm * 1e-6 and abs(rt[j] - rt[i]) <= rt_window:
                out[i] = True
                out[j] = True
            j += 1
    return out


def _chain(sorted_df, by, col, gap, rel: bool):
    """Cluster ids of a frame sorted by (by, col): a new cluster at every change of ``by`` or gap in ``col``."""
    v = sorted_df[col].to_numpy()
    b = sorted_df[by].to_numpy()
    new = np.ones(len(v), bool)
    d = np.diff(v)
    new[1:] = (b[1:] != b[:-1]) | (d > (gap * v[1:] * 1e-6 if rel else gap))
    return np.cumsum(new) - 1


def relation_anchors(relations: list[pd.DataFrame], mz_tol_ppm: float = 6.0, rt_window: float = 60.0,
                     min_frac: float = 0.6, min_snr: float = 5.0, p_chance_max: float = 0.01,
                     max_competitor: float = 0.2, min_contain: float = 0.6, dedupe: bool = False,
                     recalibrate: bool = True, noise: list[pd.DataFrame] | None = None):
    """Anchors (RT and m/z matrices, anchors x runs) from per-run relation tables.

    Returns ``(Anchors, table, members)``: ``table`` describes each anchor (rule, m/z of both ions, number of runs,
    relations linked to it), ``members`` lists the per-run relations behind it (``rel`` = row in its run's table).
    ``noise`` (from :func:`relation_noise`, aligned with ``relations``) drops relations whose windows are too busy
    for chance or hold a competitor above ``max_competitor`` of the ion's height (an isomer: when one of a pair is
    missed in some runs, the other would be matched in its place).  ``min_contain`` drops relations without a
    localised elution (chemical background), whose centre is arbitrary.  ``recalibrate`` refits the centre precision on cross-run double differences; ``dedupe`` keeps one anchor per
    group of anchors that share an ion.
    """
    J = len(relations)
    parts = []
    for j, r in enumerate(relations):
        r = r.reset_index(drop=True).assign(run=j, rel=np.arange(len(r)))
        keep = (r["snr"].to_numpy() >= min_snr) & (r["contain"].to_numpy() >= min_contain)
        if noise is not None:
            keep &= _noise_column(noise[j], "p_chance", rt_window) <= p_chance_max
            keep &= _noise_column(noise[j], "competitor", rt_window) <= max_competitor
        parts.append(r[keep].copy())
    keys = pd.concat([p[["target", "z", "steps"]] for p in parts]).drop_duplicates().reset_index(drop=True)
    kid = {tuple(v): i for i, v in enumerate(keys.itertuples(index=False))}
    for i, p in enumerate(parts):
        p["key"] = np.array([kid[(a, b, c)] for a, b, c in zip(p["target"], p["z"], p["steps"])], np.int64)
        p = p.sort_values(["key", "mz_a"], kind="stable")
        amb = _ambiguous(p["key"].to_numpy(np.int64), p["mz_a"].to_numpy(float), p["mz_b"].to_numpy(float),
                         p["rt"].to_numpy(float), mz_tol_ppm, rt_window)
        parts[i] = p[~amb]
    P = pd.concat(parts, ignore_index=True).sort_values(["key", "mz_a"], kind="stable").reset_index(drop=True)
    P["c1"] = _chain(P, "key", "mz_a", mz_tol_ppm, True)
    P = P.sort_values(["c1", "mz_b"], kind="stable").reset_index(drop=True)
    P["c2"] = _chain(P, "c1", "mz_b", mz_tol_ppm, True)
    P = P.sort_values(["c2", "rt"], kind="stable").reset_index(drop=True)
    P["cluster"] = _chain(P, "c2", "rt", rt_window, False)
    st = P.groupby("cluster").agg(n=("run", "size"), ns=("run", "nunique"), rtspan=("rt", lambda v: v.max() - v.min()))
    need = max(2, int(np.ceil(min_frac * J)))
    # two runs can deviate from the consensus in opposite directions, so the span may reach twice the window
    good = st[(st.n == st.ns) & (st.ns >= need) & (st.rtspan <= 2 * rt_window)].index
    P = P[P["cluster"].isin(good)].reset_index(drop=True)
    P["cluster"] = P["cluster"].map({c: i for i, c in enumerate(good)})
    prec = dict(relations[0].attrs.get("precision", {})) if J else {}
    if recalibrate and len(P):
        d, u2, v2 = _double_differences(P)
        if d.size >= 30:
            k1, k2 = _fit_precision(d, u2, v2)
            P["rt_se"] = centre_se(P, k1, k2)
            prec = dict(k1=k1, k2=k2, n_pairs=int(d.size), source="cross-run double differences")
    P["anchor_id"] = P["cluster"]
    P["linked"] = 1
    if dedupe and len(P):
        grp = _ion_groups(P)
        P["linked"] = pd.Series(P["cluster"]).groupby(grp).transform("nunique").to_numpy()
        score = P.groupby("cluster").agg(ns=("run", "nunique"), se=("rt_se", "median"))
        best = score.assign(g=pd.Series(grp, index=P.index).groupby(P["cluster"]).first()).sort_values(
            ["ns", "se"], ascending=[False, True]).groupby("g").head(1).index
        P = P[P["cluster"].isin(best)].reset_index(drop=True)
        P["anchor_id"] = P["cluster"].map({c: i for i, c in enumerate(sorted(best))})
    A, table, members = _matrices(P, J)
    table.attrs["precision"] = prec
    return A, table, members


@dataclass
class RelationGroups:
    """Relations grouped across runs (groups x runs matrices, NaN where a run has none)."""

    table: pd.DataFrame  # per group: rule, m/z of both ions, reference time, runs present
    rt: np.ndarray  # observed centre
    rt_ref: np.ndarray  # centre mapped to the reference (consensus) time
    area: np.ndarray  # both ions' summed signal in the fixed window
    snr: np.ndarray
    height: np.ndarray  # sum of both ions' maxima in the window
    members: pd.DataFrame  # group, run, rel (row in its run's relation table)


def to_reference(model, rt, j: int, n_iter: int = 4) -> np.ndarray:
    """Map run-``j`` times to reference time under the class-weighted mean warp.

    The class of a relation is not known before it is grouped, so the warp is averaged over classes with their
    weights; the residual is the spread between class trends at that time, which the grouping gap allows for.
    """
    rt = np.asarray(rt, float)
    mu = rt.copy()
    for _ in range(n_iter):
        mu = rt - np.einsum("k,kn->n", model.pi, model.warps(mu)[:, j])
    return mu


def relation_groups(relations: list[pd.DataFrame], model=None, mz_tol_ppm: float = 6.0, gap: float | None = None,
                    min_runs: int = 2, min_snr: float = 3.0) -> RelationGroups:
    """Group every relation across runs: cross-run grouping at the scan-graph stage, before peak detection.

    Relations are chain-clustered by rule, anchor m/z and partner m/z, then by reference time (observed centre
    mapped through ``model``, a drift mixture learned from :func:`relation_anchors`; raw time without one).  A new
    group starts at a time gap above ``gap`` (default: 4 x the predicted spread of one relation's reference time,
    at least 1 s), so isomers at one m/z separate once drift is removed.  A run with several relations in one group
    keeps the one closest to the group's median time.  Unlike anchors, nothing is required to be unique.
    """
    J = len(relations)
    parts = []
    for j, r in enumerate(relations):
        r = r.reset_index(drop=True).assign(run=j, rel=np.arange(len(r)))
        r = r[r["snr"].to_numpy() >= min_snr].copy()
        r["rt_ref"] = to_reference(model, r["rt"].to_numpy(), j) if model is not None and len(r) else r["rt"]
        parts.append(r)
    P = pd.concat(parts, ignore_index=True)
    if gap is None:
        se = float(np.nanmedian(P["rt_se"])) if len(P) else 0.0
        spread = np.sqrt((model.sigma if model is not None else 2.0) ** 2 + se * se)
        gap = max(4.0 * spread, 1.0)
    P["key"] = P.groupby(["target", "z", "steps"], sort=False).ngroup()
    P = P.sort_values(["key", "mz_a"], kind="stable").reset_index(drop=True)
    P["c1"] = _chain(P, "key", "mz_a", mz_tol_ppm, True)
    P = P.sort_values(["c1", "mz_b"], kind="stable").reset_index(drop=True)
    P["c2"] = _chain(P, "c1", "mz_b", mz_tol_ppm, True)
    P = P.sort_values(["c2", "rt_ref"], kind="stable").reset_index(drop=True)
    P["group"] = _chain(P, "c2", "rt_ref", gap, False)
    P["dev"] = np.abs(P["rt_ref"] - P.groupby("group")["rt_ref"].transform("median"))
    P = P.sort_values(["group", "run", "dev"], kind="stable").drop_duplicates(["group", "run"])
    ns = P.groupby("group")["run"].transform("size")
    P = P[ns.to_numpy() >= min_runs].copy()
    P["group"] = pd.factorize(P["group"], sort=True)[0]
    G = int(P["group"].max()) + 1 if len(P) else 0
    mats = {}
    for c in ("rt", "rt_ref", "area", "snr", "height"):
        M = np.full((G, J), np.nan)
        M[P["group"].to_numpy(), P["run"].to_numpy()] = P[c].to_numpy()
        mats[c] = M
    table = P.groupby("group").agg(target=("target", "first"), z=("z", "first"), steps=("steps", "first"),
                                   mz_a=("mz_a", "median"), mz_b=("mz_b", "median"), rt_ref=("rt_ref", "median"),
                                   n_runs=("run", "size")).reset_index(drop=True)
    table.attrs["gap"] = gap
    return RelationGroups(table, mats["rt"], mats["rt_ref"], mats["area"], mats["snr"], mats["height"],
                          P[["group", "run", "rel"]].reset_index(drop=True))


def _ion_groups(P: pd.DataFrame) -> np.ndarray:
    """Group anchor clusters whose relations share a trace, with agreeing centres, in at least two runs and in most
    of the runs where both are present (an ion with an isotope and an adduct partner anchors twice)."""
    ids = P["cluster"].to_numpy()
    n = int(ids.max()) + 1
    runs_of = [set() for _ in range(n)]
    for a, j in zip(ids.tolist(), P["run"].tolist()):
        runs_of[a].add(j)
    co: dict = {}
    Q = pd.concat([P[["run", "cluster", "rt", "rt_se"]].assign(trace=P["anchor"]),
                   P[["run", "cluster", "rt", "rt_se"]].assign(trace=P["partner"])])
    for j, d in Q.groupby("run"):
        linked = set()
        for _, e in d.groupby("trace"):
            if len(e) < 2:
                continue
            a, t, s = e["cluster"].to_numpy(), e["rt"].to_numpy(), e["rt_se"].to_numpy()
            for x in range(a.size):
                for y in range(x + 1, a.size):
                    if a[x] != a[y] and abs(t[x] - t[y]) <= 3 * np.hypot(s[x], s[y]) + 0.5:
                        linked.add((int(min(a[x], a[y])), int(max(a[x], a[y]))))
        for k in linked:
            co[k] = co.get(k, 0) + 1
    parent = np.arange(n)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for (a, b), c in co.items():
        if c >= 2 and c >= 0.5 * len(runs_of[a] & runs_of[b]):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
    return np.array([find(x) for x in range(n)])[ids]


def _matrices(P: pd.DataFrame, J: int):
    C = int(P["anchor_id"].max()) + 1 if len(P) else 0
    RT = np.full((C, J), np.nan)
    SE = np.full((C, J), np.nan)
    MZ = np.full((C, J), np.nan)
    if C:
        w = 1.0 / np.maximum(P["rt_se"].to_numpy(), 1e-3) ** 2
        g = P.assign(w=w, wt=w * P["rt"]).groupby(["anchor_id", "run"])
        agg = g.agg(w=("w", "sum"), wt=("wt", "sum"), se=("rt_se", "min"))
        a = agg.index.get_level_values(0).to_numpy()
        j = agg.index.get_level_values(1).to_numpy()
        RT[a, j] = (agg["wt"] / agg["w"]).to_numpy()
        SE[a, j] = agg["se"].to_numpy()
        # m/z for calibration: the member relation (cluster) present in most runs, its anchor ion
        rep = P.groupby(["anchor_id", "cluster"])["run"].size().reset_index().sort_values("run", ascending=False)
        rep = rep.drop_duplicates("anchor_id").set_index("anchor_id")["cluster"]
        R = P[P["cluster"].to_numpy() == rep.loc[P["anchor_id"]].to_numpy()]
        MZ[R["anchor_id"].to_numpy(), R["run"].to_numpy()] = R["mz_a"].to_numpy()
    table = P.groupby("anchor_id").agg(target=("target", "first"), mz_a=("mz_a", "median"), mz_b=("mz_b", "median"),
                                       rt=("rt", "median"), n_runs=("run", "nunique"),
                                       linked=("linked", "first")).reset_index(drop=True)
    members = P[["anchor_id", "run", "rel", "target", "z", "steps", "mz_a", "mz_b", "rt", "rt_se", "snr", "compound"]]
    return Anchors(RT, MZ, np.full((C, J), -1, np.int64), SE), table, members.reset_index(drop=True)
