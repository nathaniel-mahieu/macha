"""Score correspondence against synthetic truth."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..traces.metrics import adjusted_rand


def _cross_pairs(df, keys):
    """Number of feature pairs from different samples sharing ``keys``."""
    per = df.groupby(keys + ["sample"]).size().rename("n").reset_index()
    tot = per.groupby(keys)["n"].agg(["sum", lambda x: (x**2).sum()])
    tot.columns = ["s", "s2"]
    return float(((tot["s"] ** 2 - tot["s2"]) / 2).sum())


def correspondence_metrics(groups: list, true_ion: list, min_samples: int = 2) -> dict:
    """Pairwise cross-sample precision / recall of grouping vs true ion identity.

    Features whose true ion is -1 (noise) count against precision when grouped
    with anything, and are ignored for recall.  Ions <= -2 are excluded.
    """
    df = pd.DataFrame(dict(sample=np.concatenate([np.full(len(g), j) for j, g in enumerate(groups)]),
                           group=np.concatenate(groups), ion=np.concatenate(true_ion)))
    df = df[df["ion"] >= -1]  # ion <= -2: features without a defined correspondence (e.g. background plateaus)
    pred = _cross_pairs(df, ["group"])
    sig = df[df["ion"] >= 0]
    truth = _cross_pairs(sig, ["ion"])
    tp = _cross_pairs(sig, ["group", "ion"])
    prec = tp / pred if pred else 1.0
    rec = tp / truth if truth else 1.0
    # ions seen in >= min_samples: fraction recovered as one group holding >= 90 % of their features
    g = sig.groupby("ion")
    ions = g.filter(lambda x: x["sample"].nunique() >= min_samples)
    comp = ions.groupby("ion").apply(lambda x: x["group"].value_counts().iloc[0] / len(x) >= 0.9, include_groups=False)
    return dict(precision=prec, recall=rec, f1=2 * prec * rec / (prec + rec) if prec + rec else 0.0,
                complete_ions=float(comp.mean()) if len(comp) else np.nan, n_groups=int(df["group"].nunique()))


def rt_dispersion(aligned_rt: list, true_ion: list) -> float:
    """RMS deviation of aligned RTs around each true ion's mean (lower = better alignment)."""
    rt = np.concatenate(aligned_rt)
    ion = np.concatenate(true_ion)
    ok = ion >= 0
    df = pd.DataFrame(dict(rt=rt[ok], ion=ion[ok]))
    dev = df["rt"] - df.groupby("ion")["rt"].transform("mean")
    n = df.groupby("ion")["rt"].transform("size")
    return float(np.sqrt((dev[n > 1] ** 2).mean()))


def class_recovery(groups_class: np.ndarray, group_compound_class: np.ndarray) -> float:
    ok = group_compound_class >= 0
    return adjusted_rand(groups_class[ok], group_compound_class[ok]) if ok.sum() > 1 else float("nan")


def local_class_agreement(pred_class: np.ndarray, true_class: np.ndarray, rt: np.ndarray, window: float = 60.0,
                          min_n: int = 10) -> float:
    """Count-weighted mean ARI of drift classes within RT windows.

    Drift classes are identifiable only locally: where two class trends cross,
    labels may swap between RT regions without changing any prediction, so a
    global ARI understates recovery.
    """
    ok = true_class >= 0
    pred_class, true_class, rt = pred_class[ok], true_class[ok], rt[ok]
    bins = np.floor((rt - rt.min()) / window).astype(int)
    num = den = 0.0
    for b in np.unique(bins):
        m = bins == b
        if m.sum() >= min_n and np.unique(true_class[m]).size > 1:
            num += adjusted_rand(pred_class[m], true_class[m]) * m.sum()
            den += m.sum()
    return num / den if den else float("nan")
