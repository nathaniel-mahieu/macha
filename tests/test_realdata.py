"""Checks on public Orbitrap Astral data (run ``python -I scripts/fetch_data.py`` first).

There is no ground truth here, so these tests assert invariants and
internal-consistency properties: one centroid per scan per trace, a
physically plausible self-calibrated precision model, decoy-calibrated
annotation FDR, isotope-consistent alignment.
"""

import numpy as np
import pytest

from macha import MSRun
from macha.align import align
from macha.annotate import ProfileSource, annotate
from macha.pipeline import process_run

from .conftest import real_files

pytestmark = pytest.mark.realdata

PROT = real_files("proteomics_sc")
MET = real_files("metabolomics")


def _slice(path, rt=(500.0, 800.0)):
    return MSRun.load(path).slice(rt_range=rt)


@pytest.fixture(scope="module")
def prot_samples():
    if len(PROT) < 3:
        pytest.skip("proteomics_sc data not downloaded")
    return [process_run(_slice(p)) for p in PROT[:3]]


def test_real_traces_invariants(prot_samples):
    s = prot_samples[0]
    tr = s.traces
    sc = tr.run.scan_of_point
    ok = tr.label >= 0
    key = tr.label[ok].astype(np.int64) * (tr.run.n_scans + 1) + sc[ok]
    assert np.unique(key).size == key.size
    pm = tr.precision
    assert 0.05 < pm.a < 1.5  # Orbitrap precision floor, ppm
    assert 0.5 < pm.c < 1.5  # sigma ~ 1/SNR
    assert (tr.table["mz_sd_ppm"].median()) < 5.0


def test_real_annotation_is_calibrated(prot_samples):
    s = prot_samples[0]
    ann = annotate(s.features, ProfileSource(s.baselines), peptides=True, cycle=s.run.cycle_time)
    assert len(ann.edges) > 1000
    hi = ann.fdr[ann.fdr.threshold > 0.9]
    assert np.all(hi.decoy_fdr < 0.03)
    # Bayesian and decoy FDR within a factor of three at strict thresholds
    ratio = (hi.bayes_fdr + 1e-3) / (hi.decoy_fdr + 1e-3)
    assert np.all((ratio > 1 / 3) & (ratio < 3))
    # peptides: most isotope envelopes are multiply charged
    z = ann.roles["charge"].to_numpy()
    assert (z >= 2).sum() > 2 * (z == 1).sum()


def test_real_alignment_is_isotope_consistent(prot_samples):
    tables = [s.features for s in prot_samples]
    al = align(tables)
    assert al.drift is not None and al.drift.sigma < 3.0
    present = (al.feature_index >= 0).sum(1)
    assert (present == len(tables)).mean() > 0.2
    # isotopologues of one ion co-elute: their aligned group RTs must agree
    ann = annotate(tables[0], ProfileSource(prot_samples[0].baselines), peptides=True,
                   cycle=prot_samples[0].run.cycle_time)
    e = ann.edges[(ann.edges.target == "13C") & (ann.edges.posterior > 0.99)]
    g0 = al.feature_group[0]
    d = al.groups["rt"].to_numpy()[g0[e.i]] - al.groups["rt"].to_numpy()[g0[e.j]]
    assert np.median(np.abs(d)) < 0.5


@pytest.mark.skipif(len(MET) == 0, reason="metabolomics data not downloaded")
def test_real_metabolomics_pipeline():
    run = MSRun.load(MET[0]).slice(rt_range=(200.0, 400.0))
    s = process_run(run)
    assert len(s.features) > 100
    ann = annotate(s.features, ProfileSource(s.baselines), cycle=run.cycle_time)
    assert (ann.edges.posterior > 0.9).sum() > 50
