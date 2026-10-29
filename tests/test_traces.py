import itertools

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy.optimize import linear_sum_assignment

from macha import MSRun
from macha.traces import (PrecisionModel, adjusted_rand, detect_traces, greedy_roi, rect_roi, slice_traces,
                          trace_metrics)
from macha.traces.tracking import hungarian, solve_links


@settings(max_examples=40, deadline=None)
@given(st.integers(1, 6), st.integers(0, 4), st.integers(0, 10_000))
def test_hungarian_matches_scipy(n, extra, seed):
    rng = np.random.default_rng(seed)
    C = rng.uniform(0, 10, (n, n + extra))
    a = hungarian(C)
    r, c = linear_sum_assignment(C)
    assert C[np.arange(n), a].sum() == pytest.approx(C[r, c].sum())
    assert len(set(a)) == n


@settings(max_examples=40, deadline=None)
@given(st.integers(1, 5), st.integers(1, 5), st.integers(0, 10_000))
def test_solve_links_is_optimal(na, nb, seed):
    rng = np.random.default_rng(seed)
    pairs = [(i, j) for i in range(na) for j in range(nb) if rng.random() < 0.6]
    if not pairs:
        return
    ca = np.array([p[0] for p in pairs], np.int64)
    cb = np.array([p[1] for p in pairs], np.int64)
    g = rng.uniform(0.1, 5, len(pairs))
    m = solve_links(ca, cb, g, na, nb)
    got = sum(g[k] for k, (i, j) in enumerate(pairs) if m[i] == j)
    # brute force over all matchings
    best = 0.0
    for r in range(1, min(na, nb, len(pairs)) + 1):
        for sub in itertools.combinations(range(len(pairs)), r):
            ai = [pairs[k][0] for k in sub]
            bj = [pairs[k][1] for k in sub]
            if len(set(ai)) == r and len(set(bj)) == r:
                best = max(best, sum(g[k] for k in sub))
    assert got == pytest.approx(best)
    used = m[m >= 0]
    assert len(set(used)) == used.size


def _two_close_traces(ppm=3.0, n=80, seed=0):
    rng = np.random.default_rng(seed)
    rt = np.arange(n) * 0.5
    scans, mzs, its, lab = [], [], [], []
    for k, (m0, h) in enumerate([(300.0, 2e5), (300.0 * (1 + ppm * 1e-6), 5e4)]):
        prof = h * np.exp(-0.5 * ((np.arange(n) - 40 - 6 * k) / 8.0) ** 2)
        keep = prof > 3e3
        s = np.flatnonzero(keep)
        sig = 0.4e-6 * m0
        scans.append(s)
        mzs.append(m0 + rng.normal(0, sig, s.size))
        its.append(prof[keep])
        lab.append(np.full(s.size, k))
    scan = np.concatenate(scans)
    run = MSRun.from_points(scan, np.concatenate(mzs), np.concatenate(its).astype(np.float32), rt)
    order = np.lexsort((np.concatenate(mzs), scan))
    return run, np.concatenate(lab)[order]


def test_tracker_separates_close_traces_that_a_slice_merges():
    run, truth = _two_close_traces(ppm=3.0)
    tr = detect_traces(run, calibrate=False, precision=PrecisionModel(0.4, 1.0, 0.5), min_points=3)
    m = trace_metrics(tr.label, truth, min_points=3)
    assert m["purity"] > 0.98 and m["completeness"] > 0.9
    sl = trace_metrics(slice_traces(run, ppm=10), truth, min_points=3)
    assert sl["purity"] < m["purity"]


def test_one_point_per_scan_and_partition(met_sample):
    tr = met_sample.traces
    sc = tr.run.scan_of_point
    lab = tr.label
    ok = lab >= 0
    key = lab[ok].astype(np.int64) * (tr.run.n_scans + 1) + sc[ok]
    assert np.unique(key).size == key.size  # no trace has two centroids in one scan
    assert tr.table["n_points"].sum() == ok.sum()
    assert np.all(np.diff(tr.table["mz"].to_numpy()) >= 0)
    for i in range(0, len(tr), max(len(tr) // 25, 1)):
        p = tr.points(i)
        assert np.all(np.diff(sc[p]) > 0)


def test_precision_model_is_calibrated(met_run, met_sample):
    pm = met_sample.traces.precision
    inst = met_run.scenario.instrument
    assert pm.a == pytest.approx(inst.ppm_floor, rel=0.35)
    assert pm.c == pytest.approx(1.0, abs=0.25)
    # predicted sigma at SNR 10 within 40 % of the simulator's
    assert pm.sigma_ppm(10 * inst.noise) == pytest.approx(float(inst.ppm_sigma(10 * inst.noise)), rel=0.4)


def test_tracker_beats_reference_methods(met_run, met_sample):
    m = trace_metrics(met_sample.traces.label, met_run.ion)
    assert m["purity"] > 0.98 and m["completeness"] > 0.93 and m["ari"] > 0.95
    for lab in (greedy_roi(met_run.run), slice_traces(met_run.run)):
        ref = trace_metrics(lab, met_run.ion)
        assert m["ari"] > ref["ari"]
    assert trace_metrics(rect_roi(met_run.run), met_run.ion)["purity"] < m["purity"]


def test_gap_closing_rejoins_dropouts():
    rng = np.random.default_rng(1)
    n = 60
    s = np.array([k for k in range(n) if not (25 <= k < 31)])
    mz = 400.0 + rng.normal(0, 1e-4, s.size)
    run = MSRun.from_points(s, mz, np.full(s.size, 1e5, np.float32), np.arange(n) * 0.5)
    pm = PrecisionModel(0.5, 1.0, 0.5)
    without = detect_traces(run, calibrate=False, precision=pm, close_gaps=0, min_points=3)
    with_ = detect_traces(run, calibrate=False, precision=pm, close_gaps=8, min_points=3)
    assert len(without) == 2 and len(with_) == 1


def test_trace_metrics_bounds():
    truth = np.repeat(np.arange(10), 20)
    assert trace_metrics(truth, truth)["ari"] == pytest.approx(1.0)
    perfect = trace_metrics(truth, truth)
    assert perfect["purity"] == 1.0 and perfect["completeness"] == 1.0
    rng = np.random.default_rng(0)
    assert adjusted_rand(rng.permutation(truth), truth) < 0.1


def test_eic_dense_view(met_sample):
    tr = met_sample.traces
    i = int(np.argmax(tr.table["n_points"].to_numpy()))
    scans, rt, inten, mz = tr.eic(i, pad=3)
    assert np.isfinite(inten).sum() == tr.table["n_points"].iat[i]
    assert np.isnan(inten[:3]).all()
