"""Convenience wrappers chaining the modules.

Each step is also usable on its own; these functions only fix the order and
pass results along.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .baseline import BaselineResult, trace_baselines
from .peaks import PeakParams, detect_peaks
from .run import MSRun
from .traces import TrackParams, Traces, detect_traces


@dataclass
class SampleResult:
    run: MSRun
    traces: Traces
    baselines: BaselineResult
    features: pd.DataFrame
    timings: dict = field(default_factory=dict)


def process_run(run: MSRun, track: TrackParams | None = None, peaks: PeakParams | None = None,
                baseline_lambda: float = 1e4) -> SampleResult:
    """Traces -> baselines -> peaks for one run."""
    t0 = time.perf_counter()
    tr = detect_traces(run, track)
    t1 = time.perf_counter()
    bl = trace_baselines(tr, lam=baseline_lambda)
    t2 = time.perf_counter()
    feats = detect_peaks(tr, bl, peaks)
    t3 = time.perf_counter()
    return SampleResult(run, tr, bl, feats, dict(traces=t1 - t0, baseline=t2 - t1, peaks=t3 - t2, total=t3 - t0))


@dataclass
class BatchResult:
    samples: list
    alignment: object
    consensus: pd.DataFrame
    annotation: object
    timings: dict = field(default_factory=dict)


def process_batch(runs, peptides: bool = False, align_kw: dict | None = None, annotate_kw: dict | None = None,
                  **run_kw) -> BatchResult:
    """Process every run, align them, then annotate the consensus features.

    Consensus annotation adds cross-sample intensity correlation as evidence:
    isotopologues and adducts of one compound rise and fall together across
    samples, coincidental partners do not.
    """
    from .align import align, consensus_features
    from .annotate import annotate

    t0 = time.perf_counter()
    samples = [process_run(r, **run_kw) for r in runs]
    t1 = time.perf_counter()
    tables = [s.features for s in samples]
    al = align(tables, **(align_kw or {}))
    t2 = time.perf_counter()
    cons = consensus_features(al, tables)
    cycle = float(np.median([s.run.cycle_time for s in samples]))
    polarity = samples[0].run.polarity if samples else 1
    ann = annotate(cons, intensity_matrix=al.intensity, peptides=peptides, polarity=polarity, cycle=cycle,
                   **(annotate_kw or {}))
    t3 = time.perf_counter()
    return BatchResult(samples, al, cons, ann, dict(runs=t1 - t0, align=t2 - t1, annotate=t3 - t2, total=t3 - t0))
