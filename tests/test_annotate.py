import numpy as np
import pandas as pd
import pytest

from macha import chem
from macha.annotate import ProfileSource, annotate, default_rules, gaec, make_decoys
from macha.annotate.candidates import _modulo_residue, candidate_edges
from macha.annotate.graph import isotope_envelopes
from macha.annotate.rules import Rule
from macha.traces import adjusted_rand


def test_rule_masses():
    rules, pairs = default_rules(1)
    d = {r.name: r.delta for r in rules}
    assert d["13C"] == pytest.approx(1.0033548, abs=1e-7)
    assert d["Na-H"] == pytest.approx(21.9819442, abs=1e-6)
    assert d["H2O"] == pytest.approx(18.0105647, abs=1e-6)
    assert d["34S"] == pytest.approx(1.9957959, abs=1e-6)


def test_decoys_avoid_real_spacings():
    rules, pairs = default_rules(1, peptides=True)
    dr, dp = make_decoys(rules, pairs)
    assert dr and dp
    for d in dr:
        for r in rules:
            for z in range(1, 6):
                for n in range(1, 7):
                    assert abs(d.delta - n * r.delta / z) > 0.004


def _toy_features(rows):
    df = pd.DataFrame(rows, columns=["mz", "rt", "height"])
    df["mz_se"] = df["mz"] * 0.3e-6
    df["fwhm_s"] = 4.0
    df["snr"] = 100.0
    return df


def test_modulo_finds_chains_and_canonical_charge():
    m0 = 400.2
    c13 = chem.NEUTRON_C13
    rows = [(m0, 100.0, 1e6), (m0 + c13, 100.1, 2e5), (m0 + 2 * c13, 100.0, 3e4),  # z=1 envelope
            (600.3, 200.0, 1e6), (600.3 + c13 / 2, 200.1, 4e5), (600.3 + c13, 200.0, 1e5)]  # z=2 envelope
    f = _toy_features(rows)
    e = candidate_edges(f, [Rule("13C", c13, "isotope", 4, 3)], [], rt_window=3.0)
    got = set(zip(e.i, e.j, e.z, e.n))
    assert {(0, 1, 1, 1), (1, 2, 1, 1), (0, 2, 1, 2)} <= got
    assert {(3, 4, 2, 1), (4, 5, 2, 1)} <= got
    # (0, 1) at z=2, n=2 is the same spacing as z=1, n=1: only the canonical form is kept
    assert (0, 1, 2, 2) not in got
    # z=2 envelope also yields the z=1 n=1 edge 3->5, which is a legitimate (skip) relation
    assert (3, 5, 1, 1) in got


def test_residue_engine_finds_homologous_series():
    ch2 = 14.01565
    mz = np.array([200.1, 200.1 + ch2, 200.1 + 3 * ch2, 251.7, 200.1 + 5 * ch2])
    se = mz * 0.2e-6
    i, j, n, e, s = _modulo_residue(mz, se, 1.0, ch2, 12, 0.0, 0.3, 4.0)
    pairs = set(zip(i.tolist(), j.tolist(), n.tolist()))
    assert (0, 1, 1) in pairs and (1, 2, 2) in pairs and (0, 4, 5) in pairs
    assert not any(3 in (a, b) for a, b, _ in pairs)


def test_targeted_equals_sweep(prot_sample):
    rules, pairs = default_rules(1, peptides=True)
    dr, dp = make_decoys(rules, pairs)
    a = candidate_edges(prot_sample.features, rules + dr, pairs + dp, method="sweep")
    b = candidate_edges(prot_sample.features, rules + dr, pairs + dp)
    key = lambda d: set(zip(d.i, d.j, d.rule, d.z, d.n))  # noqa: E731
    assert len(a) > 100 and key(a) == key(b)


def test_gaec_respects_negative_edges():
    # two triangles joined by one weak positive and two strong negative edges
    i = np.array([0, 1, 0, 3, 4, 3, 2, 0, 1])
    j = np.array([1, 2, 2, 4, 5, 5, 3, 4, 5])
    w = np.array([5, 5, 5, 5, 5, 5, 1, -4, -4.0])
    lab = gaec(6, i, j, w)
    assert len(set(lab[:3])) == 1 and len(set(lab[3:])) == 1 and lab[0] != lab[3]
    # without the vetoes connected components would merge everything
    lab2 = gaec(6, i[:7], j[:7], w[:7])
    assert len(set(lab2)) == 1


def test_isotope_envelopes_pick_consistent_charge():
    edges = pd.DataFrame(dict(i=[0, 1, 0, 1], j=[1, 2, 2, 3], kind="isotope", target="13C", n=[1, 1, 1, 1],
                              z=[2, 2, 1, 1], posterior=[0.99, 0.98, 0.6, 0.55], log_odds=[6, 5, 0.4, 0.2],
                              nominal=[1, 1, 1, 1]))
    charge, root, nominal, has_lighter = isotope_envelopes(4, edges, np.zeros(4, int))
    assert list(charge[:3]) == [2, 2, 2] and list(nominal[:3]) == [0, 1, 2] and root[2] == 0


@pytest.mark.parametrize("which", ["met", "prot"])
def test_end_to_end_annotation_against_truth(which, met_run, met_sample, prot_run, prot_sample):
    from macha.synth.truth import feature_truth

    sr, s, pep = (met_run, met_sample, False) if which == "met" else (prot_run, prot_sample, True)
    truth = feature_truth(s.features, s.traces, sr.ion, sr.ions, sr.compounds)
    ann = annotate(s.features, ProfileSource(s.baselines), peptides=pep, cycle=sr.run.cycle_time)
    comp = truth["compound"].to_numpy()
    e = ann.edges
    true = (comp[e.i] == comp[e.j]) & (comp[e.i] >= 0)
    acc = e.posterior.to_numpy() >= 0.5
    assert true[acc].mean() > 0.97
    assert (acc & true).sum() / true.sum() > 0.85
    # the decoy-based and Bayesian FDR estimates agree at strict thresholds
    hi = ann.fdr[ann.fdr.threshold > 0.9]
    assert np.all(hi.bayes_fdr < 0.05) and np.all(hi.decoy_fdr < 0.05)
    a = comp >= 0
    assert adjusted_rand(ann.labels[a], comp[a]) > 0.9
    r = ann.roles
    ok = a & np.isfinite(r.neutral_mass.to_numpy(dtype=float))
    err = np.abs(r.neutral_mass[ok] - truth.neutral_mass[ok]) / truth.neutral_mass[ok] * 1e6
    assert (err < 10).mean() > 0.9
