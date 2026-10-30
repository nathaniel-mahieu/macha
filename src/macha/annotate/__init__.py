"""Isotope / adduct / neutral-loss graph with Bayesian edge evaluation.

>>> ann = annotate(sample.features, profiles=ProfileSource(sample.baselines), polarity=1)
>>> ann.edges        # target edges with evidence, log-odds and posterior
>>> ann.groups       # one row per inferred compound
>>> ann.roles        # per feature: group, charge, isotope index, adduct, neutral mass
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.special import expit

from .bayes import EdgeModel, fdr_curves
from .candidates import candidate_edges
from .evidence import ProfileSource, compute_evidence
from .graph import gaec, infer_compounds, pair_weights
from .rules import PairRule, Rule, default_rules, make_decoys


@dataclass
class Annotation:
    edges: pd.DataFrame
    decoys: pd.DataFrame
    series: pd.DataFrame
    model: EdgeModel
    labels: np.ndarray
    roles: pd.DataFrame
    groups: pd.DataFrame
    fdr: pd.DataFrame
    timings: dict = field(default_factory=dict)


def annotate(features: pd.DataFrame, profiles: ProfileSource | None = None, polarity: int = 1,
             peptides: bool = False, rules: list[Rule] | None = None, pair_rules: list[PairRule] | None = None,
             intensity_matrix=None, gate: float = 4.0, sys_ppm: float = 0.3, rt_window: float | None = None,
             cycle: float = 0.5, n_decoys: int = 8, series: bool = False, min_posterior: float = 0.5,
             compound_ppm: float = 5.0) -> Annotation:
    """Build, evaluate and partition the mass-difference graph of a feature table."""
    t0 = time.perf_counter()
    features = features.reset_index(drop=True)
    if rules is None or pair_rules is None:
        r0, p0 = default_rules(polarity, peptides, series)
        rules = rules if rules is not None else r0
        pair_rules = pair_rules if pair_rules is not None else p0
    drules, dprules = make_decoys(rules, pair_rules, n_decoys)
    cand = candidate_edges(features, rules + drules, pair_rules + dprules, gate, sys_ppm, rt_window)
    t1 = time.perf_counter()
    ser = cand[cand["kind"] == "series"].reset_index(drop=True)
    cand = cand[cand["kind"] != "series"].reset_index(drop=True)
    ev = compute_evidence(cand, features, profiles, intensity_matrix,
                          "peptide" if peptides else "metabolite", cycle)
    snr = features["snr"].to_numpy() if "snr" in features else np.full(len(features), 20.0)
    ev["snr_min"] = np.minimum(snr[ev["i"].to_numpy()], snr[ev["j"].to_numpy()]) if len(ev) else []
    t2 = time.perf_counter()
    targets = ev[~ev["decoy"].astype(bool)].reset_index(drop=True)
    decoys = ev[ev["decoy"].astype(bool)].reset_index(drop=True)
    model = EdgeModel(gate=gate).fit(targets, decoys) if len(decoys) and len(targets) else EdgeModel(gate=gate)
    for d in (targets, decoys):
        if len(d):
            d["log_odds"] = model.log_odds(d)
            d["posterior"] = expit(d["log_odds"].to_numpy())
        else:
            d["log_odds"] = []
            d["posterior"] = []
    t3 = time.perf_counter()
    if len(targets):
        pw = pair_weights(targets)
        labels = gaec(len(features), pw["a"].to_numpy(), pw["b"].to_numpy(), pw["w"].to_numpy())
    else:
        labels = np.arange(len(features))
    roles, groups = infer_compounds(features, targets, labels, polarity, compound_ppm, min_posterior, peptides)
    t4 = time.perf_counter()
    per_target = max(len(decoys) / max(len(targets), 1), 1e-9)
    k = (len(drules) + len(dprules)) / max(len(rules) + len(pair_rules), 1)  # decoy rules per target rule
    fdr = fdr_curves(targets, decoys, k) if len(targets) and len(decoys) else pd.DataFrame()
    return Annotation(edges=targets, decoys=decoys, series=ser, model=model, labels=labels, roles=roles,
                      groups=groups, fdr=fdr,
                      timings=dict(candidates=t1 - t0, evidence=t2 - t1, em=t3 - t2, graph=t4 - t3, total=t4 - t0,
                                   decoy_ratio=per_target))


__all__ = ["annotate", "Annotation", "ProfileSource", "EdgeModel", "candidate_edges", "compute_evidence", "gaec",
           "default_rules", "make_decoys", "Rule", "PairRule", "fdr_curves"]
