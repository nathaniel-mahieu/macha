import numpy as np
import pytest

from macha import MSRun, chem
from macha.annotate import ProfileSource, annotate
from macha.annotate.rules import default_rules, make_decoys
from macha.annotate.graph import gaec, gaec_coherent
from macha.annotate.scan_graph import _instances, build_targets, repartition, scan_graph
from macha.annotate.scan_graph_eval import (background_mixed, block_ari, compound_recall, feature_graph_groups,
                                            group_majors, group_purity, scan_graph_groups, segment_truth)
from macha.synth import make_run, metabolomics_scenario
from macha.synth.truth import feature_truth
from macha.traces import PrecisionModel, detect_traces


def test_targets_merge_overlapping_isotope_steps_and_map_pairs():
    rules, pairs = default_rules(1)
    T = build_targets(rules, pairs)
    two = T[(T.z == 1) & (T.alpha == 1.0) & (T.beta > 1.99) & (T.beta < 2.01) & ~T.decoy.astype(bool)]
    assert len(two) == 1 and "13C" in two.target.iat[0] and "34S" in two.target.iat[0]
    assert (np.gcd(T.n.astype(int), T.z.astype(int)) == 1).all()
    # [M+H]+ -> [M+2H]2+ map
    pr = T[T.rule == "[M+H]/[M+2H]"].iloc[0]
    M = 500.0
    mh = chem.ADDUCTS["[M+H]+"].mz(M)
    m2h = chem.ADDUCTS["[M+2H]2+"].mz(M)
    assert pr.alpha * mh + pr.beta == pytest.approx(m2h, abs=1e-9)
    dr, dp = make_decoys(rules, pairs, 3)
    D = build_targets(dr, dp)
    assert D.decoy.astype(bool).all() and len(D) > len(T)


def _pair_run(missing_partner_scans=()):
    """An anchor trace and its 13C partner 1.0034 Da up; optionally drop partner centroids in some scans."""
    rng = np.random.default_rng(0)
    rt = np.arange(40) * 0.5
    sc, mz, it = [], [], []
    for s in range(40):
        h = 2e5 * np.exp(-0.5 * ((s - 20) / 5) ** 2)
        if h < 5e3:
            continue
        sc.append(s), mz.append(300.1 + rng.normal(0, 1e-4)), it.append(h)
        if s not in missing_partner_scans:
            sc.append(s), mz.append(300.1 + chem.NEUTRON_C13 + rng.normal(0, 1e-4)), it.append(0.15 * h)
    sc, mz, it = np.array(sc), np.array(mz), np.array(it, np.float32)
    return MSRun.from_points(sc, mz, it, rt)


def test_instances_bridge_missing_partner():
    run = _pair_run(missing_partner_scans={18, 22})
    pm = PrecisionModel(0.5, 10.0, 0.5)
    anchors = run.mz < 300.5
    args = (run.mz, run.intensity.astype(float), run.offsets, run.scan_of_point.astype(np.int64), anchors,
            np.array([1.0]), np.array([chem.NEUTRON_C13]), np.array([0.0]), pm.a, pm.b, pm.c, 0.3, 4.0)
    i0, j0, *_ = _instances(*args, 0, 1.0)
    i1, j1, *_ = _instances(*args, 1, 1.0)
    sc = run.scan_of_point
    assert set(sc[i0]) == set(sc[anchors]) - {18, 22}  # no bridge: anchors without a same-scan partner drop out
    assert set(sc[i1]) == set(sc[anchors])  # bridged: every anchor scan supported
    assert np.all(np.abs(sc[j1] - sc[i1]) <= 1)
    same = sc[i1] != 18
    assert np.all(sc[j1][same & (sc[i1] != 22)] == sc[i1][same & (sc[i1] != 22)])  # same scan preferred


def test_gaec_coherent_breaks_chains_across_time():
    # 0 -- 1 -- 2 strongly linked, but node 2 elutes far from 0 and 1
    i, j, w = np.array([0, 1]), np.array([1, 2]), np.array([3.0, 3.0])
    x = np.array([10.0, 10.5, 110.0])
    assert len(set(gaec(3, i, j, w))) == 1
    lab = gaec_coherent(3, i, j, w, x, np.ones(3))
    assert lab[0] == lab[1] != lab[2]
    # no precision = no opinion: identical to plain GAEC
    assert len(set(gaec_coherent(3, i, j, w, x, np.zeros(3)))) == 1


@pytest.fixture(scope="module")
def low_abundance():
    sr = make_run(metabolomics_scenario(n_compounds=250, run_time=500.0, n_background=15,
                                        log10_height=(4.3, 0.8)), seed=5)
    from macha.pipeline import process_run

    s = process_run(sr.run)
    g = scan_graph(s.traces)
    return sr, s, g


def test_scan_graph_finds_more_low_snr_compounds_than_feature_graph(low_abundance):
    sr, s, g = low_abundance
    ft = feature_truth(s.features, s.traces, sr.ion, sr.ions, sr.compounds)
    ann = annotate(s.features, ProfileSource(s.baselines), cycle=sr.run.cycle_time)
    noise = sr.scenario.instrument.noise
    gs, gf = scan_graph_groups(g, s.traces, sr.ion), feature_graph_groups(ann, ft)
    rs = compound_recall(gs, sr.ions, sr.compounds, noise, group_major=group_majors(gs, sr.ions)).set_index("snr_bin")
    rf = compound_recall(gf, sr.ions, sr.compounds, noise, group_major=group_majors(gf, sr.ions)).set_index("snr_bin")
    assert rs.loc["10-30", "recall"] > rf.loc["10-30", "recall"] + 0.15
    assert rs.loc["30-100", "recall"] >= rf.loc["30-100", "recall"] - 0.05


def test_scan_graph_relations_are_calibrated(low_abundance):
    sr, s, g = low_abundance
    seg = g.segments.assign(true=segment_truth(g, s.traces, sr.ion, sr.ions))
    hi = seg[seg.posterior >= 0.9]
    assert len(hi) > 200 and hi.true.mean() > 0.9
    # expected calibration error (count-weighted |predicted - observed| over deciles) is small
    b = np.clip((seg.posterior * 10).astype(int), 0, 9)
    c = seg.groupby(b).agg(p=("posterior", "mean"), f=("true", "mean"), n=("true", "size"))
    ece = float((np.abs(c.p - c.f) * c.n).sum() / c.n.sum())
    assert ece < 0.08
    # Bayesian and decoy FDR estimates agree at strict thresholds
    f = g.fdr[g.fdr.threshold > 0.85]
    assert np.all(np.abs(f.bayes_fdr - f.decoy_fdr) < 0.05)


def test_compound_chromatogram_and_blocks(low_abundance):
    sr, s, g = low_abundance
    assert len(g.compounds) > 50
    c = int(g.compounds.sort_values("n_edges").compound.iat[-1])
    scans, rt, y = g.compound_chromatogram(s.traces, c)
    assert y.max() > 0 and scans.size == rt.size == y.size
    b = g.blocks[g.blocks.trace >= 0]
    # blocks of one trace never overlap in time
    for t, d in b.groupby("trace"):
        d = d.sort_values("s0")
        assert np.all(d.s0.to_numpy()[1:] > d.s1.to_numpy()[:-1])


def test_time_prior_keeps_members_within_a_peak_width(low_abundance):
    sr, s, g = low_abundance
    g0 = repartition(g, s.traces.run, time_prior=False)
    # the learned peak spread matches the simulated elution (sd 1.5-3.5 s plus EMG tails)
    assert 1.0 < g.prior.sigma < 4.0 and g.prior.tau < g.prior.sigma
    c, c0 = (x.compounds[x.compounds.n_traces >= 2] for x in (g, g0))
    assert (c.center_span <= 3 * g.prior.fwhm).mean() > 0.99
    assert (c.center_span > 3 * g.prior.fwhm).mean() < (c0.center_span > 3 * g.prior.fwhm).mean()
    gs, gs0 = scan_graph_groups(g, s.traces, sr.ion), scan_graph_groups(g0, s.traces, sr.ion)
    assert group_purity(gs, sr.ions) > group_purity(gs0, sr.ions) + 0.03
    assert background_mixed(gs, sr.ions, sr.compounds) < 0.5 * background_mixed(gs0, sr.ions, sr.compounds)
    # grouping of eluting compounds is not traded away
    ari = block_ari(g, s.traces, sr.ion, sr.ions, sr.compounds, "eluting")
    assert ari >= block_ari(g0, s.traces, sr.ion, sr.ions, sr.compounds, "eluting") - 0.02


@pytest.mark.realdata
def test_real_scan_graph_fdr_agrees_with_decoys():
    from .conftest import real_files

    files = [f for f in real_files("metabolomics") if "multirange" not in f.name]
    if not files:
        pytest.skip("metabolomics data not downloaded")
    run = MSRun.load(files[0]).slice(rt_range=(200.0, 500.0))
    g = scan_graph(detect_traces(run))
    f = g.fdr[(g.fdr.threshold > 0.85) & (g.fdr.accepted > 100)]
    assert len(f) and np.all(f.decoy_fdr < 0.06)
    ratio = (f.bayes_fdr + 2e-3) / (f.decoy_fdr + 2e-3)
    assert np.all((ratio > 1 / 3) & (ratio < 3))
    # Astral peaks are a few seconds wide; compounds stay within a few peak widths
    assert 0.5 < g.prior.sigma < 6.0
    c = g.compounds[g.compounds.n_traces >= 2]
    assert (c.center_span > 3 * g.prior.fwhm).mean() < 0.06
