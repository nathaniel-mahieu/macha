import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy.signal import find_peaks, peak_prominences

from macha.baseline import als_baseline, censored_arpls, estimate_floor, whittaker_solve
from macha.peaks import cwt_peaks, detect_peaks_1d, fit_emg_mixture, persistence, sspp_peaks
from macha.synth import emg_height_normalised


@settings(max_examples=30, deadline=None)
@given(st.integers(5, 60), st.floats(0.1, 1e4), st.integers(0, 10_000))
def test_whittaker_matches_dense_solve(n, lam, seed):
    rng = np.random.default_rng(seed)
    y = rng.normal(size=n)
    w = rng.uniform(0.01, 1, n)
    D = np.diff(np.eye(n), 2, axis=0)
    dense = np.linalg.solve(np.diag(w) + lam * D.T @ D, w * y)
    assert np.allclose(whittaker_solve(y, w, lam), dense, atol=1e-8 * max(1, np.abs(dense).max()))


@settings(max_examples=40, deadline=None)
@given(st.integers(0, 10_000))
def test_persistence_equals_prominence(seed):
    rng = np.random.default_rng(seed)
    y = np.cumsum(rng.normal(size=80)) + rng.normal(size=80) * 0.3
    idx, pers = persistence(y)
    pk, _ = find_peaks(y)
    prom = peak_prominences(y, pk)[0]
    d = dict(zip(idx.tolist(), pers.tolist()))
    for p, pr in zip(pk, prom):
        # scipy treats the array ends as bases; persistence only dies at a saddle toward higher
        # ground.  The definitions coincide for peaks with higher ground on both sides.
        if (y[:p] > y[p]).any() and (y[p + 1:] > y[p]).any():
            assert d[int(p)] == pytest.approx(pr, rel=1e-9, abs=1e-9)


def _signal(rng, n=300, peaks=((80, 4, 1e5), (150, 3, 4e4), (160, 3, 3e4)), base=5e3, noise=1e3, thr=3e3):
    x = np.arange(n, dtype=float)
    true_base = base + 2e3 * np.sin(x / 60)
    y = true_base.copy()
    for c, s, h in peaks:
        y += h * np.exp(-0.5 * ((x - c) / s) ** 2)
    y += rng.normal(0, noise, n)
    y[y < thr] = np.nan
    return x, y, true_base


def test_censored_arpls_tracks_baseline_under_peaks():
    rng = np.random.default_rng(0)
    x, y, true_base = _signal(rng)
    obs = np.isfinite(y)
    z = censored_arpls(np.where(obs, y, 0), obs, np.full(y.size, 3e3), 1e4, 10)
    mask = np.ones(y.size, bool)
    for c in (80, 150, 160):
        mask[c - 15: c + 15] = False
    assert np.sqrt(np.mean((z[mask] - true_base[mask]) ** 2)) < 1.5e3
    # baseline never climbs into the peak apex region
    assert z[80] < 1e4
    # the locf/zero-filled classic ALS is worse here
    zl = als_baseline(y, 1e4, 0.01, fill="zero")
    assert np.sqrt(np.mean((zl[mask] - true_base[mask]) ** 2)) > np.sqrt(np.mean((z[mask] - true_base[mask]) ** 2))


def test_detect_peaks_1d_resolves_overlap():
    rng = np.random.default_rng(1)
    x, y, true_base = _signal(rng)
    pk = detect_peaks_1d(y, baseline=true_base, noise=1e3)
    apex = np.sort(pk["apex"].to_numpy())
    assert len(apex) == 3
    assert np.allclose(apex, [80, 150, 160], atol=2)
    assert pk.loc[pk["apex"].sub(80).abs().idxmin(), "height"] == pytest.approx(1e5, rel=0.1)


def test_reference_detectors_run():
    rng = np.random.default_rng(2)
    x, y, b = _signal(rng)
    assert np.any(np.abs(sspp_peaks(y, b, 1e3) - 80) <= 3)
    assert np.any(np.abs(cwt_peaks(y, b, 1e3) - 80) <= 3)


def test_emg_fit_recovers_two_components():
    t = np.linspace(0, 60, 240)
    y = 1e5 * emg_height_normalised(t, 25, 2.0, 1.5) + 4e4 * emg_height_normalised(t, 33, 2.0, 1.5)
    y = y + np.random.default_rng(3).normal(0, 300, t.size)
    fit = fit_emg_mixture(t, y, [(9e4, 24, 2.5, 1.0), (5e4, 34, 2.5, 1.0)])
    fit = fit.sort_values("mu")
    assert fit["mu"].to_numpy() == pytest.approx([25, 33], abs=0.5)
    assert fit["height"].to_numpy() == pytest.approx([1e5, 4e4], rel=0.12)
    assert fit.attrs["success"]


def test_floor_and_trace_baselines(met_run, met_sample):
    fl = estimate_floor(met_run.run)
    thr = met_run.scenario.instrument.threshold
    assert np.median(fl.grid) == pytest.approx(thr, rel=0.3)
    bl = met_sample.baselines
    assert bl.offsets[-1] == bl.intensity.size == bl.baseline.size
    assert np.all(bl.baseline >= 0) and np.all(bl.noise > 0)


def test_peak_recall_against_truth(met_sample, met_truth, met_run):
    f = met_sample.features
    assert (met_truth["ion"] >= 0).mean() > 0.9
    ok = met_truth["true_rt"].notna().to_numpy()
    err = np.abs(f["rt"].to_numpy()[ok] - met_truth["true_rt"].to_numpy()[ok])
    assert np.median(err) < 0.6
    mz_err = np.abs(f["mz"] - met_truth["true_mz"]) / f["mz"] * 1e6
    assert np.nanmedian(mz_err) < 1.0
    # standard errors are honest once the run's systematic calibration offset is removed: the
    # simulator adds a per-run ppm offset and slope, which no single-run estimate can see
    ppm = ((f["mz"] - met_truth["true_mz"]) / met_truth["true_mz"] * 1e6).to_numpy(float)
    se_ppm = (f["mz_se"] / f["mz"] * 1e6).to_numpy(float)
    rt = f["rt"].to_numpy()
    good = np.isfinite(ppm) & (met_truth["purity"].to_numpy() > 0.95)
    trend = np.polyval(np.polyfit(rt[good], ppm[good], 1), rt)
    z = (ppm - trend)[good] / se_ppm[good]
    assert 0.6 < np.median(np.abs(z)) / 0.674 < 1.8
