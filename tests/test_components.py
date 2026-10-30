import numpy as np
import pandas as pd
import pytest

from macha import MSRun
from macha.align.drift import cv_drift
from macha.align.edges import relation_anchors, relation_groups, relation_noise, relation_table
from macha.annotate.scan_graph import scan_graph
from macha.components import emg, fit_components, link_components, predict_apex, quantify
from macha.synth import make_batch, metabolomics_scenario
from macha.synth.lc import emg as synth_emg, emg_mode
from macha.traces import detect_traces


def test_emg_matches_the_simulator_and_integrates_to_one():
    t = np.linspace(-20, 60, 4001)
    for sig, tau in ((2.0, 0.05), (2.0, 1.5), (1.0, 6.0)):
        f = emg(t, 3.0, sig, tau)
        assert np.allclose(f, synth_emg(t, 3.0, sig, tau), rtol=1e-6, atol=1e-12)
        assert np.trapezoid(f, t) == pytest.approx(1.0, abs=1e-3)


def _merged_run(dt_sd, dm_ppm, partners, seed=0, n_scans=200, cycle=0.5, sig=2.0, tau=1.0):
    """Two sources whose main ions are ``dm_ppm`` apart (they coalesce into one centroid at the intensity-weighted
    m/z) and ``dt_sd`` peak sd apart in time.  ``partners``: list of (m/z per source, ratio per source) for the
    other ions; a partner of one source only has ratio 0 for the other."""
    rng = np.random.default_rng(seed)
    rt = np.arange(n_scans) * cycle
    mode1 = 45.0
    off = float(emg_mode(0.0, sig, tau))
    mus = np.array([mode1 - off, mode1 + dt_sd * sig - off])
    heights = np.array([1e6, 6e5])
    prof = np.array([synth_emg(rt, m, sig, tau) for m in mus])
    prof = prof / prof.max(1, keepdims=True) * heights[:, None]
    ions = [((400.0, 400.0 * (1 + dm_ppm * 1e-6)), (1.0, 1.0))] + list(partners)
    scans, mzs, ints = [], [], []
    for s in range(n_scans):
        for (m1, m2), (r1, r2) in ions:
            i1, i2 = prof[0, s] * r1, prof[1, s] * r2
            tot = i1 + i2
            if tot < 2e3:
                continue
            m = (m1 * i1 + m2 * i2) / tot
            sd = np.sqrt(0.3 ** 2 + (60.0 / np.sqrt(tot)) ** 2) * 1e-6 * m
            scans.append(s)
            mzs.append(m + rng.normal(0, sd))
            ints.append(tot * np.exp(rng.normal(0, 0.05)))
        k = rng.poisson(25)  # noise centroids
        scans += [s] * k
        mzs += list(rng.uniform(200, 800, k))
        ints += list(np.exp(rng.normal(np.log(3e3), 0.4, k)))
    scans, mzs, ints = map(np.asarray, (scans, mzs, ints))
    o = np.lexsort((mzs, scans))
    run = MSRun.from_points(scans[o], mzs[o], ints[o], rt, 1)
    tr = detect_traces(run)
    true_modes = mus + off
    return tr, ions, true_modes


def _relations(tr, ions):
    """Hand-made relation table: the main ion related to every partner, all one compound."""
    mz = tr.table["mz"].to_numpy()
    tid = [int(np.argmin(np.abs(mz - 0.5 * (m[0] + m[1])) + 1e3 * (tr.table["n_points"].to_numpy() < 10)))
           for m, _ in ions]
    rows = [dict(anchor=tid[0], partner=t, compound=0, s0=60, s1=120, height_a=1e6, height_b=1e5, height=1.1e6,
                 rt=46.0) for t in tid[1:]]
    rel = pd.DataFrame(rows)
    rel.attrs["sigma"] = 2.0
    return rel, tid


def test_split_with_partner_templates_recovers_both_apexes():
    # each source has a partner of its own: the partners pin the two profiles
    tr, ions, modes = _merged_run(1.0, 3.0, [((401.0034, 401.0034), (0.3, 0.0)), ((402.9871, 402.9871), (0.0, 0.4))])
    rel, tid = _relations(tr, ions)
    C = fit_components(rel, tr, cv=0.05, max_dm_ppm=10.0)
    assert C.fits["split"].iloc[0]
    rts = np.sort(C.table["rt"].to_numpy())
    assert np.abs(rts - modes).max() < 0.6
    # the merged ion is shared; each partner belongs to one side
    io = C.ions.set_index(["trace", "k"])["share"]
    assert io[(tid[0], 0)] > 0.3 and io[(tid[0], 1)] > 0.2
    assert max(io.get((tid[1], 0), 0), io.get((tid[1], 1), 0)) > 0.9


def test_mz_trajectory_separates_coalesced_sources():
    # both ions are coalesced pairs (main ion and its M+1), one peak sd apart: intensity alone barely changes, the
    # m/z trajectory decides, and the split recovers apexes, the amplitude share and the sources' m/z difference
    tr, ions, modes = _merged_run(1.0, 4.0, [((401.0034, 401.0034 * (1 + 4e-6)), (0.25, 0.25))], seed=3)
    rel, tid = _relations(tr, ions)
    with_mz = fit_components(rel, tr, cv=0.05, max_dm_ppm=10.0)
    f = with_mz.fits.iloc[0]
    assert f["split"] and f["gain_mz"] > 5 * f["penalty"] and f["gain_int"] < f["penalty"]
    t = with_mz.table.set_index("k")
    assert np.abs(np.sort(t["rt"].to_numpy()) - modes).max() < 0.4
    first, second = (0, 1) if t["rt"][0] < t["rt"][1] else (1, 0)
    assert t["area"][first] / t["area"].sum() == pytest.approx(1 / 1.6, abs=0.08)
    io = with_mz.ions[with_mz.ions.trace == tid[0]].set_index("k")
    dm = (io.loc[first, "mz_source"] - io.loc[second, "mz_source"]) / 400.0 * 1e6
    assert dm == pytest.approx(-4.0, abs=1.5)


def test_close_sources_are_flagged_merged_not_split():
    tr, ions, modes = _merged_run(0.4, 4.0, [((401.0034, 401.0034 * (1 + 4e-6)), (0.25, 0.25))], seed=1)
    rel, _ = _relations(tr, ions)
    f = fit_components(rel, tr, cv=0.05, max_dm_ppm=10.0).fits.iloc[0]
    assert f["merged"] and not f["split"]


def test_single_source_is_not_split():
    tr, ions, modes = _merged_run(0.0, 0.0, [((401.0034, 401.0034), (0.3, 0.3)), ((402.9871, 402.9871), (0.4, 0.4))],
                                  seed=5)
    rel, _ = _relations(tr, ions)
    C = fit_components(rel, tr, cv=0.05, max_dm_ppm=10.0)
    assert not C.fits["split"].iloc[0]
    assert C.table["rt"].iloc[0] == pytest.approx(modes[0], abs=0.4)
    assert C.table["sigma"].iloc[0] == pytest.approx(2.0, rel=0.15)


@pytest.fixture(scope="module")
def batch():
    sc = metabolomics_scenario(n_compounds=150, run_time=360.0, n_background=10, noise_points_per_scan=60.0)
    b = make_batch(sc, n_samples=4, n_drift_classes=2, class_drift=10.0, global_drift=5.0, seed=7)
    trs, rels, noise = [], [], []
    for sr in b.runs:
        tr = detect_traces(sr.run)
        rel = relation_table(scan_graph(tr), tr)
        trs.append(tr)
        rels.append(rel)
        noise.append(relation_noise(rel, tr))
    return b, trs, rels, noise


def test_components_link_across_runs_and_quantify_every_run(batch):
    b, trs, rels, noise = batch
    J = len(trs)
    comps = [fit_components(r, t) for r, t in zip(rels, trs)]
    assert all(0.01 < c.cv < 0.3 for c in comps)
    A, _, _ = relation_anchors(rels, noise=noise)
    model, _ = cv_drift(A.rt, A.se, K_max=2, n_basis=10)
    cr = link_components(comps, relation_groups(rels, model), model)
    assert len(cr.table) > 80 and (cr.table.n_runs >= 2).all()
    apex, sd, half = predict_apex(cr, model, J, comps)
    Q = [quantify(cr, j, trs[j].run, trs[j].precision, apex[:, j], sd[:, j], half[:, j], cv=comps[j].cv)
         for j in range(J)]
    area = np.stack([q.area for q in Q], 1)
    assert (area > 0).mean() > 0.97  # re-query reaches the runs the per-run fits missed
    # centred log areas follow the true per-run heights of the component's compound
    sr0, tr0, c0 = b.runs[0], trs[0], comps[0]
    m0 = cr.members[cr.members.run == 0]
    truth = []
    for x, comp in zip(m0.xc, m0.component):
        t = int(c0.ions[c0.ions.component == comp].sort_values("area").trace.iloc[-1])
        lab = sr0.compound[tr0.points(t)]
        lab = lab[lab >= 0]
        truth.append((x, np.bincount(lab).argmax() if lab.size else -1))
    errs = []
    for x, c in truth:
        if c < 0 or b.drift_class[c] < 0:
            continue
        la, lt = np.log(area[x]), np.log(b.height[c])
        if np.isfinite(la).all():
            errs.append(np.sqrt(np.mean(((la - la.mean()) - (lt - lt.mean())) ** 2)))
    assert len(errs) > 50
    assert np.median(errs) < 0.05
