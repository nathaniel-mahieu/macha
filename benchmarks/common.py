"""Shared datasets and helpers for the module benchmarks.

Every benchmark writes ``reports/<module>.html`` (an explorer page) and
``reports/<module>.json`` (raw numbers).  Synthetic datasets carry ground
truth; real Orbitrap Astral datasets (``scripts/fetch_data.py``) are used for
truth-free diagnostics, runtime and internal-consistency checks.
"""

from __future__ import annotations

import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from macha import MSRun  # noqa: E402
from macha.synth import (hard_scenario, make_run, metabolomics_scenario,  # noqa: E402
                         proteomics_scenario)

REPORTS = ROOT / "reports"
CACHE = ROOT / "data" / "synthetic"
ASTRAL = ROOT / "data" / "astral"

NAV = [("Index", "index.html"), ("I/O & runs", "io.html"), ("Traces", "traces.html"),
       ("Baseline", "baseline.html"), ("Peaks", "peaks.html"), ("Annotation graph", "annotate.html"),
       ("Scan graph", "scan_graph.html"), ("DIA grouping", "dia.html"),
       ("Alignment", "align.html"), ("Cross-run graph", "xrun.html"),
       ("Components", "components.html")]

SYNTHETIC = {
    "syn-metabolomics": dict(make=lambda: make_run(metabolomics_scenario(), seed=1),
                             note="600 metabolites + 80 background ions, 120k, 0.5 s cycle"),
    "syn-hard": dict(make=lambda: make_run(hard_scenario(), seed=1),
                     note="900 metabolites, 150 near-isobaric (1-5 ppm) co-eluting pairs, 400 noise peaks/scan"),
    "syn-proteomics": dict(make=lambda: make_run(proteomics_scenario(), seed=1),
                           note="1200 tryptic peptides, z 1-5, 240k, 0.8 s cycle"),
    "syn-dense-proteomics": dict(make=lambda: make_run(proteomics_scenario(n_compounds=5000,
                                                                          noise_points_per_scan=600.0), seed=1),
                                 note="5000 peptides, 600 noise peaks/scan (Astral-like density)"),
}


def synthetic(name: str):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{name}.pkl"
    if path.exists():
        return pickle.loads(path.read_bytes())
    sr = SYNTHETIC[name]["make"]()
    path.write_bytes(pickle.dumps(sr))
    return sr


def real_runs(dataset: str) -> list[Path]:
    """Standard single-range runs of a dataset, samples first and blanks last."""
    d = ASTRAL / dataset
    if not d.exists():
        return []
    files = [f for f in d.glob("*.npz") if "multirange" not in f.name]
    return sorted(files, key=lambda f: (f.name.lower().startswith("blank"), f.name))


def real_run(dataset: str, i: int = 0, rt_range=None) -> MSRun | None:
    files = real_runs(dataset)
    if len(files) <= i:
        return None
    run = MSRun.load(files[i])
    return run.slice(rt_range=rt_range) if rt_range else run


REAL_LABEL = {"proteomics_sc": "Astral DIA 250 pg (PXD049412)", "proteomics": "Astral DIA 50 ng (ProteoBench M7)",
              "metabolomics": "Astral NIST fecal, RP+ (MSV000093526)"}


class Timer:
    def __enter__(self):
        self.t = time.perf_counter()
        return self

    def __exit__(self, *a):
        self.s = time.perf_counter() - self.t


def save_json(name: str, obj):
    REPORTS.mkdir(parents=True, exist_ok=True)

    def conv(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)

    (REPORTS / f"{name}.json").write_text(json.dumps(obj, indent=1, default=conv))
