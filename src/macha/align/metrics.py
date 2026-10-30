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


def loo_prediction_error(model, T: np.ndarray) -> np.ndarray:
    """Leave-one-run-out RT prediction error of a drift model on an RT matrix (compounds x runs, NaN = missing).

    For every compound and run j, the class posterior and consensus time are estimated from the other runs and the
    RT in run j is predicted as the posterior mean.  On true compound RTs this scores the warps themselves,
    whatever anchors they were learned from; its floor is the compound-level RT jitter.
    """
    T = np.asarray(T, float)
    C, J = T.shape
    err = np.full((C, J), np.nan)
    logpi = np.log(np.maximum(model.pi, 1e-12))
    for j in range(J):
        obs = np.isfinite(T)
        obs[:, j] = False
        ok = np.isfinite(T[:, j]) & (obs.sum(1) >= 1)
        if not ok.any():
            continue
        Tm, om = T[ok], obs[ok]
        n = om.sum(1)
        pred = np.empty((model.K, Tm.shape[0]))
        ll = np.empty((model.K, Tm.shape[0]))
        for k in range(model.K):
            mu = np.where(om, Tm, 0.0).sum(1) / n
            for _ in range(5):
                W = model.warps(mu)[k]  # (J, n)
                mu = np.where(om, Tm - W.T, 0.0).sum(1) / n
            W = model.warps(mu)[k]
            r = np.where(om, Tm - mu[:, None] - W.T, 0.0)
            ll[k] = logpi[k] - 0.5 * (r * r).sum(1) / model.sigma ** 2
            pred[k] = mu + W[j]
        ll -= ll.max(0, keepdims=True)
        resp = np.exp(ll)
        resp /= resp.sum(0, keepdims=True)
        err[ok, j] = (resp * pred).sum(0) - Tm[:, j]
    return err


def design_mixture(x, offsets, names=None, free_means: bool = False, n_iter: int = 300) -> dict:
    """Fit log ratios to a mixture with known spacings: a common shift plus one Gaussian per offset, and a uniform
    background over the data range.

    For a mixed-species design (e.g. log2 A/B = 0, +1, -2 for human, yeast, E. coli) the components' widths
    measure quantification precision and the background weight the share of groups that fit no species, i.e.
    mis-grouped or noise.  ``free_means`` lets every component move on its own, which checks the design.
    """
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    off = np.asarray(offsets, float)
    names = list(names) if names is not None else [str(o) for o in off]
    S = off.size
    if x.size < 10 * S:
        return dict(n=int(x.size))
    lo, hi = float(x.min()), float(x.max())
    f0 = 1.0 / max(hi - lo, 1e-9)
    delta = float(np.median(x))
    mu = delta + off
    sd = np.full(S, 0.3)
    w = np.full(S + 1, 1.0 / (S + 1))
    for _ in range(n_iter):
        dens = np.vstack([w[s] * np.exp(-0.5 * ((x - mu[s]) / sd[s]) ** 2) / (np.sqrt(2 * np.pi) * sd[s])
                          for s in range(S)] + [np.full(x.size, w[S] * f0)])
        r = dens / np.maximum(dens.sum(0, keepdims=True), 1e-300)
        nk = r.sum(1)
        w = nk / x.size
        if free_means:
            mu = (r[:S] * x).sum(1) / np.maximum(nk[:S], 1e-12)
        else:
            prec = nk[:S] / sd ** 2
            delta = float(((r[:S] * (x[None, :] - off[:, None])).sum(1) / sd ** 2).sum() / max(prec.sum(), 1e-12))
            mu = delta + off
        sd = np.sqrt((r[:S] * (x[None, :] - mu[:, None]) ** 2).sum(1) / np.maximum(nk[:S], 1e-12))
        sd = np.maximum(sd, 0.02)
    out = dict(n=int(x.size), background=float(w[S]), loglik=float(np.log(dens.sum(0)).sum()))
    if not free_means:
        out["shift"] = delta
    for s, nm in enumerate(names):
        out[f"mean_{nm}"] = float(mu[s])
        out[f"sd_{nm}"] = float(sd[s])
        out[f"weight_{nm}"] = float(w[s])
    return out
