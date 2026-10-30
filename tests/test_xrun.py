from types import SimpleNamespace

import numpy as np
import pytest

from macha import MSRun
from macha.align import find_anchors, fit_drift
from macha.align.drift import cv_drift, support_range
from macha.align.edges import relation_anchors, relation_groups, relation_noise, relation_table
from macha.align.metrics import design_mixture, loo_prediction_error
from macha.annotate.scan_graph import scan_graph
from macha.pipeline import process_run
from macha.synth import make_batch, metabolomics_scenario
from macha.traces import detect_traces, window_noise

from .test_align import _drift_matrix


def test_window_noise_counts_interference_and_excludes_own_trace():
    n = 100
    scan = np.concatenate([np.arange(n), np.arange(0, n, 2), np.arange(n)])
    mz = np.concatenate([np.full(n, 300.0), np.full(n // 2, 300.0 + 300e-6 * 3), np.full(n, 310.0)])
    inten = np.concatenate([np.full(n, 1e5), np.full(n // 2, 2e3), np.full(n, 1e5)])
    o = np.lexsort((mz, scan))
    run = MSRun.from_points(scan[o], mz[o], inten[o], np.arange(n) * 0.5, 1)
    label = np.concatenate([np.zeros(n, int), np.full(n // 2, -1), np.ones(n, int)])[o]
    tr = SimpleNamespace(run=run, label=label)
    w = window_noise(tr, [300.0, 300.0, 310.0], 0, n - 1, 10.0, exclude=[0, -1, 1])
    assert w.hit_frac[0] == pytest.approx(0.5)  # own trace excluded, every other scan interfered
    assert w.hit_frac[1] == pytest.approx(1.0)  # nothing excluded: the trace itself occupies every scan
    assert w.hit_frac[2] == 0.0
    assert np.exp(w.noise_log_mean[0]) == pytest.approx(2e3)


def test_sigma_survives_heterogeneous_anchor_errors():
    # a few anchors with a large stated error used to drive the moment estimate negative and pin sigma to its floor
    T, cls, mu = _drift_matrix(K=2, noise=0.3)
    S = np.full(T.shape, 0.05)
    S[:30] = 5.0
    m = fit_drift(T, 2, t_range=(0, 600), S=S)
    assert m.sigma == pytest.approx(0.3, rel=0.3)


def test_cv_drift_chooses_classes_and_predicts_held_out_runs():
    T, cls, mu = _drift_matrix(K=2, noise=0.2)
    m, cv = cv_drift(T, None, K_max=3, t_range=(0, 600))
    assert m.K >= 2
    assert set(cv.columns) >= {"K", "lam", "rmse"}
    e = loo_prediction_error(m, T)
    assert np.sqrt(np.nanmean(e ** 2)) < 0.4  # noise 0.2 with one run held out


@pytest.fixture(scope="module")
def small_batch():
    sc = metabolomics_scenario(n_compounds=150, run_time=360.0, n_background=15, noise_points_per_scan=60.0)
    b = make_batch(sc, n_samples=4, n_drift_classes=2, class_drift=16.0, global_drift=6.0, seed=5)
    rels, noise, tables, comps = [], [], [], []
    for sr in b.runs:
        tr = detect_traces(sr.run)
        g = scan_graph(tr)
        rel = relation_table(g, tr)
        rels.append(rel)
        noise.append(relation_noise(rel, tr))
        tables.append(process_run(sr.run).features)
        comps.append((tr, sr.compound))
    return b, rels, noise, tables, comps


def test_relation_table_centres_and_precision(small_batch):
    b, rels, noise, tables, comps = small_batch
    r = rels[0]
    assert len(r) > 200
    assert {"rt", "rt_se", "snr", "mz_a", "mz_b", "compound"} <= set(r.columns)
    p = r.attrs["precision"]
    assert p["n_pairs"] >= 30 and 0 < p["k1"] < 10
    assert (r.rt_se > 0).all() and (r.snr > 0).all()
    assert (noise[0].p_chance_60.between(0, 1)).all() and (noise[0].competitor_10 <= noise[0].competitor_60).all()


def test_relation_anchors_match_runs_and_align_like_features(small_batch):
    b, rels, noise, tables, comps = small_batch
    A, tab, mem = relation_anchors(rels, noise=noise)
    assert A.rt.shape[1] == 4 and A.rt.shape[0] >= 50
    assert tab.attrs["precision"].get("source") == "cross-run double differences"
    # anchors are one compound: centres of a member relation agree with the compound's true RT up to its offset
    tr, comp = comps[0]
    m0 = mem[mem.run == 0]
    rel0 = rels[0].iloc[m0.rel.to_numpy()]
    truec = np.array([np.bincount(c[c >= 0]).argmax() if (c >= 0).any() else -1
                      for c in (comp[tr.points(int(a))] for a in rel0.anchor)])
    ok = truec >= 0
    assert ok.mean() > 0.9
    true = b.rt[b.drift_class >= 0]
    mr, _ = cv_drift(A.rt, A.se, K_max=3, n_basis=10)
    Af = find_anchors(tables)
    mf, _ = cv_drift(Af.rt, Af.se, K_max=3, n_basis=10)
    er = np.sqrt(np.nanmean(loo_prediction_error(mr, true) ** 2))
    ef = np.sqrt(np.nanmean(loo_prediction_error(mf, true) ** 2))
    raw = np.sqrt(np.nanmean((true - np.nanmean(true, 1, keepdims=True)) ** 2))
    assert er < 0.25 * raw
    assert er < 1.5 * ef


def test_relation_groups_group_one_compound_across_runs(small_batch):
    b, rels, noise, tables, comps = small_batch
    A, _, _ = relation_anchors(rels, noise=noise)
    m, _ = cv_drift(A.rt, A.se, K_max=2, n_basis=10)
    g = relation_groups(rels, m)
    assert len(g.table) > 200 and (g.table.n_runs >= 2).all()
    assert g.height.shape == (len(g.table), 4)
    # members of a group are the same relation in every run: same true compound
    truth = []
    for j, (tr, comp) in enumerate(comps):
        r = rels[j]
        truth.append(np.array([np.bincount(c[c >= 0]).argmax() if (c >= 0).any() else -1
                               for c in (comp[tr.points(int(a))] for a in r.anchor)]))
    mem = g.members.assign(tc=[truth[j][i] for j, i in zip(g.members.run, g.members.rel)])
    mem = mem[mem.tc >= 0]
    pure = mem.groupby("group").tc.nunique() == 1
    assert pure.mean() > 0.9


def test_design_mixture_recovers_species_modes():
    rng = np.random.default_rng(1)
    x = np.concatenate([rng.normal(0.1, 0.2, 3000), rng.normal(1.1, 0.3, 900), rng.normal(-1.9, 0.3, 300),
                        rng.uniform(-5, 4, 200)])
    m = design_mixture(x, [0, 1, -2], ["h", "y", "e"])
    assert m["shift"] == pytest.approx(0.1, abs=0.03)
    assert m["sd_h"] == pytest.approx(0.2, rel=0.15) and m["sd_y"] == pytest.approx(0.3, rel=0.2)
    assert m["weight_y"] == pytest.approx(900 / 4400, abs=0.03)
    f = design_mixture(x, [0, 1, -2], ["h", "y", "e"], free_means=True)
    assert f["mean_e"] == pytest.approx(-1.9, abs=0.1)
