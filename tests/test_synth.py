import numpy as np
import pytest

from macha import chem
from macha.synth import emg, emg_height_normalised, emg_mode, make_batch, make_run, metabolomics_scenario
from macha.synth.instrument import coalesce


def test_monoisotopic_and_adduct_masses():
    glc = chem.parse_formula("C6H12O6")
    assert chem.monoisotopic_mass(glc) == pytest.approx(180.063388, abs=1e-5)
    M = chem.monoisotopic_mass(glc)
    assert chem.ADDUCTS["[M+H]+"].mz(M) == pytest.approx(181.070665, abs=1e-5)
    assert chem.ADDUCTS["[M+Na]+"].mz(M) == pytest.approx(203.052609, abs=1e-5)
    assert chem.ADDUCTS["[M-H]-"].mz(M) == pytest.approx(179.056112, abs=1e-5)
    assert chem.ADDUCTS["[M+2H]2+"].mz(M) == pytest.approx((M + 2 * chem.PROTON) / 2, abs=1e-6)
    for name, a in chem.ADDUCTS.items():
        assert a.neutral_mass(a.mz(M)) == pytest.approx(M, abs=1e-8), name


def test_isotope_fine_structure_normalised_and_resolved():
    f = chem.parse_formula("C20H30N4O6S")
    m, p = chem.isotope_fine_structure(f)
    assert p.sum() == pytest.approx(1.0)
    assert m[0] == pytest.approx(chem.monoisotopic_mass(f), abs=1e-6)
    # M+1 fine structure: 13C and 15N peaks are distinct
    m1 = m[(m > m[0] + 0.9) & (m < m[0] + 1.1)]
    assert np.any(np.abs(m1 - (m[0] + 1.0033548)) < 1e-4)
    assert np.any(np.abs(m1 - (m[0] + 0.9970349)) < 1e-4)
    # at low resolving power they merge into one centroid near the averaged spacing
    cm, cp = chem.centroid_pattern(m, p, lambda x: 5_000.0)
    assert np.sum((cm > m[0] + 0.9) & (cm < m[0] + 1.1)) == 1


def test_emg_mode_and_normalisation():
    t = np.linspace(-20, 150, 17001)
    for mu, s, tau in [(0, 2, 0.1), (0, 2, 3), (5, 1, 8)]:
        y = emg(t, mu, s, tau)
        assert np.trapezoid(y, t) == pytest.approx(1.0, abs=2e-3)
        assert emg_mode(mu, s, tau) == pytest.approx(t[np.argmax(y)], abs=0.02)
        assert emg_height_normalised(t, mu, s, tau).max() == pytest.approx(1.0, abs=1e-3)


def test_coalescence_merges_unresolved_ions():
    scan = np.array([0, 0, 0, 1], np.int64)
    mz = np.array([100.0, 100.0005, 100.1, 100.0])
    it = np.array([10.0, 30.0, 5.0, 1.0])
    lab = np.array([1, 2, 3, 4], np.int64)
    s, m, i, l, merged = coalesce(scan, mz, it, lab, np.full(4, 1e-5))
    assert len(s) == 3
    assert m[0] == pytest.approx((100 * 10 + 100.0005 * 30) / 40)
    assert i[0] == 40 and l[0] == 2 and merged[0] and not merged[1]


def test_make_run_truth_consistency(met_run):
    sr = met_run
    run = sr.run
    assert run.n_points == sr.ion.size
    assert sr.ion.max() < len(sr.ions)
    sig = sr.ion >= 0
    # every reported analyte centroid is within a few ppm of its ion's true m/z
    err = (run.mz[sig] - sr.ions["mz"].to_numpy()[sr.ion[sig]]) / run.mz[sig] * 1e6
    assert np.median(np.abs(err[~sr.merged[sig]])) < 2.0
    assert run.intensity.min() >= sr.scenario.instrument.threshold * 0.999


def test_make_run_is_deterministic():
    sc = metabolomics_scenario(n_compounds=30, run_time=120.0, n_background=2)
    a, b = make_run(sc, seed=5), make_run(sc, seed=5)
    assert np.array_equal(a.run.mz, b.run.mz) and np.array_equal(a.ion, b.ion)


def test_make_batch_drift_structure():
    sc = metabolomics_scenario(n_compounds=60, run_time=200.0, n_background=3)
    b = make_batch(sc, n_samples=4, n_drift_classes=2, class_drift=10.0, global_drift=0.0, wobble=0.0,
                   jitter=0.0, seed=1)
    assert b.rt.shape == (len(b.compounds), 4)
    ok = b.drift_class >= 0
    dev = b.rt[ok] - b.compounds["rt"].to_numpy()[ok, None]
    # first and last injections drift in opposite directions (injection order -0.5 .. 0.5)
    assert np.corrcoef(dev[:, 0], dev[:, -1])[0, 1] < -0.9
    # each compound follows its class's trend exactly (no jitter / wobble / global drift here)
    t = b.compounds["rt"].to_numpy()[ok]
    cls = b.drift_class[ok]
    for k in (0, 1):
        expect = b.injection_order[-1] * b.warps["class"][k](t[cls == k])
        assert np.allclose(dev[cls == k, -1], expect, atol=1e-9)
    grid = np.linspace(0, 200, 50)
    assert np.max(np.abs(b.warps["class"][0](grid) - b.warps["class"][1](grid))) > 3.0
