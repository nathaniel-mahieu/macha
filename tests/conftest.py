from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from macha.pipeline import process_run
from macha.synth import make_batch, make_run, metabolomics_scenario, proteomics_scenario
from macha.synth.truth import feature_truth

DATA = Path(__file__).resolve().parents[1] / "data" / "astral"


@pytest.fixture(scope="session")
def met_run():
    sc = metabolomics_scenario(n_compounds=200, run_time=400.0, n_background=20, noise_points_per_scan=60.0)
    return make_run(sc, seed=11)


@pytest.fixture(scope="session")
def met_sample(met_run):
    return process_run(met_run.run)


@pytest.fixture(scope="session")
def met_truth(met_run, met_sample):
    return feature_truth(met_sample.features, met_sample.traces, met_run.ion, met_run.ions, met_run.compounds)


@pytest.fixture(scope="session")
def prot_run():
    sc = proteomics_scenario(n_compounds=300, run_time=500.0, n_background=10, noise_points_per_scan=80.0)
    return make_run(sc, seed=12)


@pytest.fixture(scope="session")
def prot_sample(prot_run):
    return process_run(prot_run.run)


@pytest.fixture(scope="session")
def batch():
    sc = metabolomics_scenario(n_compounds=250, run_time=500.0, n_background=10, n_isomers=80,
                               noise_points_per_scan=50.0)
    return make_batch(sc, n_samples=6, n_drift_classes=2, class_drift=14.0, global_drift=5.0, seed=21)


@pytest.fixture(scope="session")
def batch_processed(batch):
    samples = [process_run(sr.run) for sr in batch.runs]
    truths = [feature_truth(s.features, s.traces, sr.ion, sr.ions, sr.compounds) for s, sr in zip(samples, batch.runs)]
    true_ion = [np.where(t["kind"].to_numpy() == "background", -2, t["ion"].to_numpy()) for t in truths]
    return samples, truths, true_ion


def real_files(dataset: str):
    d = DATA / dataset
    return sorted(d.glob("*.npz")) if d.exists() else []
