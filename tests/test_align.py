import numpy as np
import pytest

from macha.align import align, align_dtw, density_groups, dtw_path, fit_drift, select_drift
from macha.align.anchors import Anchors, fit_mz_calibration
from macha.align.metrics import correspondence_metrics, local_class_agreement, rt_dispersion


def _drift_matrix(n=300, J=6, K=2, amp=10.0, noise=0.2, seed=0):
    rng = np.random.default_rng(seed)
    mu = rng.uniform(20, 580, n)
    cls = rng.integers(0, K, n)
    order = np.linspace(-0.5, 0.5, J)
    prof = [lambda t: amp * np.sin(t / 90), lambda t: -amp * np.cos(t / 140), lambda t: amp * (t / 300 - 1)][:K]
    T = np.empty((n, J))
    for c in range(n):
        T[c] = mu[c] + order * prof[cls[c]](mu[c]) + rng.normal(0, noise, J)
    return T, cls, mu


def test_drift_mixture_selects_k_and_recovers_classes():
    T, cls, mu = _drift_matrix(K=2)
    best, fits = select_drift(T, 3, t_range=(0, 600))
    assert best.K == 2
    assert best.sigma == pytest.approx(0.2, rel=0.35)
    # classes are identifiable locally (labels may swap where the two trends cross)
    sep = np.abs(np.sin(mu / 90) + np.cos(mu / 140)) > 0.3
    assert local_class_agreement(best.resp.argmax(1)[sep], cls[sep], mu[sep], window=60.0) > 0.8
    one, _ = select_drift(_drift_matrix(K=1)[0], 3, t_range=(0, 600))
    assert one.K == 1


def test_drift_warps_are_centred_and_invertible():
    T, cls, mu = _drift_matrix(K=1, noise=0.05)
    m = fit_drift(T, 1, t_range=(0, 600))
    g, W = m.grid(200)
    assert np.allclose(W.sum(1), 0, atol=1e-8)
    t = np.array([100.0, 300.0])
    back = m.to_reference(t + m.warp(t, 0, 2), 0, 2)
    assert np.allclose(back, t, atol=1e-3)


def test_drift_is_robust_to_outliers():
    T, cls, mu = _drift_matrix(K=1, noise=0.1)
    T2 = T.copy()
    rng = np.random.default_rng(5)
    idx = rng.choice(T.size, T.size // 20, replace=False)
    T2.flat[idx] += rng.choice([-1, 1], idx.size) * 30.0  # 5 % gross mismatches
    m = fit_drift(T2, 1, t_range=(0, 600))
    assert m.sigma < 0.3 and m.eps > 0.02


def test_mz_calibration_recovers_offsets():
    rng = np.random.default_rng(2)
    C, J = 200, 4
    rt = rng.uniform(0, 600, (C, 1)) + rng.normal(0, 1, (C, J))
    true_mz = rng.uniform(100, 900, C)
    offs = np.array([-1.0, 0.0, 0.6, 1.2])
    ramp = np.array([0.0, 0.5, -0.5, 0.0])
    ppm = offs[None, :] + ramp[None, :] * rt / 600 - (offs.mean() + ramp.mean() * rt / 600)
    MZ = true_mz[:, None] * (1 + 1e-6 * (ppm + rng.normal(0, 0.2, (C, J))))
    cal = fit_mz_calibration(Anchors(rt, MZ, np.zeros((C, J), int)), (0, 600))
    for j in range(J):
        corr = cal.correct(MZ[:, j], rt[:, j], j)
        resid = (corr - true_mz) / true_mz * 1e6
        assert np.std(resid) < 0.35


def test_dtw_recovers_shift():
    t = np.arange(400.0)
    a = np.exp(-0.5 * ((t - 150) / 5) ** 2) + np.exp(-0.5 * ((t - 260) / 7) ** 2)
    b = np.exp(-0.5 * ((t - 160) / 5) ** 2) + np.exp(-0.5 * ((t - 270) / 7) ** 2)
    pi, pj, cost = dtw_path(a, b, 30)
    m = pi == 150
    assert np.median(pj[m]) == pytest.approx(160, abs=2)


def test_end_to_end_alignment_beats_references(batch, batch_processed):
    samples, truths, true_ion = batch_processed
    tables = [s.features for s in samples]
    al = align(tables)
    m = correspondence_metrics(al.feature_group, true_ion)
    assert m["precision"] > 0.97 and m["recall"] > 0.93
    rms = rt_dispersion(al.aligned_rt, true_ion)
    raw = rt_dispersion([t["rt"].to_numpy() for t in tables], true_ion)
    glob = align(tables, K=1)
    assert rms < 0.5 * raw
    assert rms < rt_dispersion(glob.aligned_rt, true_ion)
    dens = correspondence_metrics(density_groups(tables), true_ion)
    assert m["f1"] > dens["f1"]
    # structural guarantees: one feature per sample per group; probabilities in [0, 1]
    for j, g in enumerate(al.feature_group):
        assert np.unique(g).size == g.size
        assert np.all((al.membership[j] >= 0) & (al.membership[j] <= 1))
    assert al.intensity.shape == (len(al.groups), len(tables))


def test_membership_probabilities_are_informative(batch, batch_processed):
    samples, truths, true_ion = batch_processed
    al = align([s.features for s in samples])
    ion = np.concatenate(true_ion)
    grp = np.concatenate(al.feature_group)
    memb = np.concatenate(al.membership)
    ok = ion >= 0
    import pandas as pd

    df = pd.DataFrame(dict(g=grp[ok], ion=ion[ok], p=memb[ok]))
    maj = df.groupby("g")["ion"].agg(lambda x: x.value_counts().index[0])
    df["correct"] = df["ion"].to_numpy() == maj.loc[df["g"]].to_numpy()
    # confident memberships are more often correct than unconfident ones
    hi, lo = df[df.p > 0.9], df[df.p < 0.5]
    if len(lo) > 10:
        assert hi.correct.mean() > lo.correct.mean()
    assert hi.correct.mean() > 0.97


def test_dtw_reference_runs(batch, batch_processed):
    samples, _, true_ion = batch_processed
    chrom = [(sr.run.rt, sr.run.bpc()) for sr in batch.runs]
    al = align_dtw([s.features for s in samples], chrom)
    assert correspondence_metrics(al.feature_group, true_ion)["precision"] > 0.9


def test_consensus_annotation_uses_cross_sample_correlation(batch, batch_processed):
    from macha.align import consensus_features
    from macha.annotate import annotate

    samples, truths, true_ion = batch_processed
    tables = [s.features for s in samples]
    al = align(tables)
    cons = consensus_features(al, tables)
    assert len(cons) == len(al.groups) and np.isfinite(cons["mz_se"]).all()
    ann = annotate(cons, intensity_matrix=al.intensity, cycle=batch.runs[0].run.cycle_time)
    assert ann.model.use_x
    # map groups to true compounds (majority) and score accepted edges
    comp = np.concatenate([t["compound"].to_numpy() for t in truths])
    grp = np.concatenate(al.feature_group)
    import pandas as pd

    df = pd.DataFrame(dict(g=grp, c=comp))
    maj = df[df.c >= 0].groupby("g")["c"].agg(lambda x: x.value_counts().index[0])
    gc = np.full(len(cons), -1)
    gc[maj.index.to_numpy()] = maj.to_numpy()
    e = ann.edges[ann.edges.posterior >= 0.9]
    true = (gc[e.i] == gc[e.j]) & (gc[e.i] >= 0)
    assert len(e) > 100 and true.mean() > 0.95
    # cross-sample correlation is higher for true relations than for decoys
    assert ann.edges.xcorr[ann.edges.posterior > 0.9].median() > ann.decoys.xcorr.median()
