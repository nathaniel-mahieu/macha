"""Cross-sample alignment: drift-trend mixture + probabilistic correspondence.

>>> al = align([s.features for s in samples])
>>> al.groups          # consensus features: m/z, reference RT, drift class posterior
>>> al.intensity       # groups x samples matrix
>>> al.membership[j]   # per-feature membership probability in sample j
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .anchors import Anchors, MzCalibration, apex_se, find_anchors, fit_mz_calibration
from .correspond import mz_channels, solve_all
from .drift import DriftModel, fit_drift, select_drift, support_range
from .dtw import dtw_path, dtw_warp


@dataclass
class Alignment:
    groups: pd.DataFrame
    feature_group: list
    membership: list
    aligned_rt: list
    intensity: np.ndarray
    feature_index: np.ndarray
    drift: DriftModel | None
    drift_fits: list
    anchors: Anchors | None
    calibration: MzCalibration | None
    timings: dict = field(default_factory=dict)


def _zero_drift(J, t_range, sigma):
    from .drift import bspline_basis

    knots, degree = bspline_basis(t_range[0], t_range[1], 6)
    nb = len(knots) - degree - 1
    return DriftModel(1, knots, degree, np.zeros((1, J, nb)), np.ones(1), sigma, 0.0)


def align(tables: list[pd.DataFrame], K: int | None = None, K_max: int = 4, drift: str = "fit",
          mz_ppm: float = 5.0, anchor_ppm: float = 6.0, rt_window: float = 60.0, min_frac: float = 0.6,
          sigma_t: float | None = None, birth_cost: float = 10.0, gate: float = 4.0, sys_ppm: float = 0.5,
          n_refine: int = 2, refit: bool = True, n_basis: int = 12, lam: float = 1.0, calibrate_mz: bool = True,
          nu: float = 5.0, refit_max: int = 6000, value: str = "height", seed: int = 0,
          anchors: Anchors | None = None) -> Alignment:
    """Align feature tables from several runs.

    ``drift='fit'`` learns the drift-trend mixture (K fixed or chosen by BIC
    up to ``K_max``); ``drift='none'`` matches on raw RTs with tolerance
    ``sigma_t`` (used for reference comparisons and pre-corrected RTs).
    ``anchors`` replaces the feature anchors, e.g. relation anchors from the
    scan graph (:func:`macha.align.edges.relation_anchors`).
    """
    t0 = time.perf_counter()
    J = len(tables)
    tables = [t.reset_index(drop=True) for t in tables]
    allrt = np.concatenate([t["rt"].to_numpy() for t in tables])
    t_range = (float(allrt.min()) - 1.0, float(allrt.max()) + 1.0)
    anchors = anchors if anchors is not None else find_anchors(tables, anchor_ppm, rt_window, min_frac)
    cal = fit_mz_calibration(anchors, t_range) if (calibrate_mz and len(anchors.rt) >= 10) else None
    t1 = time.perf_counter()
    fits = []
    # the spline domain is where drift is observed (robust range); warps are held constant beyond it
    drift_range = support_range(anchors.rt) if len(anchors.rt) else t_range
    if drift == "fit" and len(anchors.rt) >= 10:
        if K is None:
            model, fits = select_drift(anchors.rt, K_max, t_range=drift_range, n_basis=n_basis, lam=lam, seed=seed,
                                       S=anchors.se)
        else:
            model = fit_drift(anchors.rt, K, t_range=drift_range, n_basis=n_basis, lam=lam, seed=seed, S=anchors.se)
            fits = [model]
    else:
        if sigma_t is None:
            dev = anchors.rt - np.nanmedian(anchors.rt, 1, keepdims=True) if len(anchors.rt) else np.array([2.0])
            sigma_t = float(1.4826 * np.nanmedian(np.abs(dev))) * 1.5 + 0.5 if np.isfinite(dev).any() else 3.0
        model = _zero_drift(J, t_range, sigma_t)
    t2 = time.perf_counter()

    rt = allrt
    smp = np.concatenate([np.full(len(t), j, np.int64) for j, t in enumerate(tables)])
    fidx = np.concatenate([np.arange(len(t)) for t in tables])
    mz = np.concatenate([t["mz"].to_numpy() for t in tables])
    se = np.concatenate([t["mz_se"].to_numpy() if "mz_se" in t else t["mz"].to_numpy() * 1e-6 for t in tables])
    if cal is not None:
        mz = np.concatenate([cal.correct(t["mz"].to_numpy(), t["rt"].to_numpy(), j) for j, t in enumerate(tables)])
    se = np.maximum(se, mz * 1e-8)
    rse = np.concatenate([_apex_se(t) for t in tables])

    def solve(m: DriftModel):
        G, W = m.grid()
        order, offs = mz_channels(mz, se, mz_ppm, gate)
        st = sigma_t if (drift != "fit" or sigma_t is not None) else m.sigma
        return solve_all(order, offs, rt, rse, mz, se, smp, J, G, np.ascontiguousarray(W), np.log(m.pi), float(st),
                         sys_ppm, birth_cost, gate, n_refine, nu)

    group, memb, gmu, gr, gmz, gn = solve(model)
    t_refit = 0.0
    if refit and drift == "fit" and len(anchors.rt) >= 10:
        # re-learn the drift mixture from confidently matched groups (many more than the anchors),
        # subsampled evenly over RT: a few thousand groups pin down smooth warps
        tr0 = time.perf_counter()
        T, S = _group_matrix(group, smp, rt, rse, gn, memb, J, min_frac)
        if T.shape[0] > refit_max:
            o = np.argsort(np.nanmedian(T, 1))
            keep = o[np.linspace(0, o.size - 1, refit_max).astype(int)]
            T, S = T[keep], S[keep]
        if T.shape[0] >= 20:
            model = fit_drift(T, model.K, t_range=support_range(T), n_basis=n_basis, lam=lam, seed=seed, S=S)
            t_refit = time.perf_counter() - tr0
            group, memb, gmu, gr, gmz, gn = solve(model)
    t3 = time.perf_counter()

    n_groups = gn.size
    FI = np.full((n_groups, J), -1, np.int64)
    FI[group, smp] = fidx
    vals = np.concatenate([t[value].to_numpy() if value in t else np.full(len(t), np.nan) for t in tables])
    X = np.full((n_groups, J), np.nan)
    X[group, smp] = vals
    MB = np.full((n_groups, J), np.nan)
    MB[group, smp] = memb
    cls = gr.argmax(1)
    mu = (gmu * gr).sum(1)
    groups = pd.DataFrame(dict(mz=gmz, rt=mu, n_samples=gn, drift_class=cls, class_prob=gr.max(1),
                               mean_membership=np.nanmean(MB, 1)))
    for k in range(gr.shape[1]):
        groups[f"p_class{k}"] = gr[:, k]
    # aligned (reference-time) RT of every feature under its group's most probable class
    G, W = model.grid()
    al_rt = rt.copy()
    for k in range(W.shape[0]):
        for j in range(J):
            m = (smp == j) & (cls[group] == k)
            if m.any():
                tt = rt[m]
                mm = tt.copy()
                for _ in range(4):
                    mm = tt - np.interp(mm, G, W[k, j])
                al_rt[m] = mm
    split = np.cumsum([len(t) for t in tables])[:-1]
    return Alignment(groups=groups, feature_group=np.split(group, split), membership=np.split(memb, split),
                     aligned_rt=np.split(al_rt, split), intensity=X, feature_index=FI,
                     drift=model if drift == "fit" else None, drift_fits=fits, anchors=anchors, calibration=cal,
                     timings=dict(anchors=t1 - t0, drift=t2 - t1, drift_refit=t_refit, correspond=t3 - t2 - t_refit,
                                  total=time.perf_counter() - t0))


def _apex_se(t: pd.DataFrame) -> np.ndarray:
    """Apex-location standard error from peak width and signal-to-noise."""
    if "fwhm_s" not in t or "snr" not in t:
        return np.zeros(len(t))
    return apex_se(t["fwhm_s"].to_numpy(), t["snr"].to_numpy())


def _group_matrix(group, smp, rt, rse, gn, memb, J, min_frac):
    need = max(2, int(np.ceil(min_frac * J)))
    T = np.full((gn.size, J), np.nan)
    S = np.full((gn.size, J), np.nan)
    ok = memb > 0.8
    T[group[ok], smp[ok]] = rt[ok]
    S[group[ok], smp[ok]] = rse[ok]
    keep = (gn >= need) & (np.isfinite(T).sum(1) >= need)
    return T[keep], S[keep]


def consensus_features(al: Alignment, tables: list[pd.DataFrame]) -> pd.DataFrame:
    """One row per aligned group with the columns :func:`macha.annotate.annotate` expects.

    m/z standard errors are combined across members, width / SNR / height are
    member means / maxima, and ``rt`` is the consensus time.
    """
    FI = al.feature_index
    G, J = FI.shape
    out = {c: np.full((G, J), np.nan) for c in ("mz_se", "fwhm_s", "snr", "height")}
    for j, t in enumerate(tables):
        ok = FI[:, j] >= 0
        rows = FI[ok, j]
        for c in out:
            if c in t:
                out[c][ok, j] = t[c].to_numpy()[rows]
    with np.errstate(invalid="ignore", divide="ignore"):
        se = 1.0 / np.sqrt(np.nansum(1.0 / out["mz_se"] ** 2, 1))
    return pd.DataFrame(dict(mz=al.groups["mz"].to_numpy(), mz_se=se, rt=al.groups["rt"].to_numpy(),
                             fwhm_s=np.nanmean(out["fwhm_s"], 1), snr=np.nanmax(out["snr"], 1),
                             height=np.nanmean(out["height"], 1), n_samples=al.groups["n_samples"].to_numpy()))


# ------------------------------------------------------------------ reference methods
def align_dtw(tables, chromatograms, reference: int | None = None, band_seconds: float = 90.0, **kw) -> Alignment:
    """Warp each run onto a reference by DTW of a chromatogram (e.g. BPC), then match without drift."""
    ref = reference if reference is not None else len(tables) // 2
    rt_ref, y_ref = chromatograms[ref]
    new = []
    for j, (t, (rtq, yq)) in enumerate(zip(tables, chromatograms)):
        t = t.copy()
        if j != ref:
            f = dtw_warp(rtq, yq, rt_ref, y_ref, band_seconds)
            t["rt"] = f(t["rt"].to_numpy())
        new.append(t)
    return align(new, drift="none", **kw)


def density_groups(tables, mz_ppm: float = 5.0, bw: float = 5.0):
    """xcms-like density grouping: m/z channels split at RT gaps wider than ``bw`` (no correction)."""
    rt = np.concatenate([t["rt"].to_numpy() for t in tables])
    mz = np.concatenate([t["mz"].to_numpy() for t in tables])
    se = np.concatenate([t["mz_se"].to_numpy() if "mz_se" in t else t["mz"].to_numpy() * 1e-6 for t in tables])
    order, offs = mz_channels(mz, se, mz_ppm)
    group = np.empty(rt.size, np.int64)
    base = 0
    for c in range(offs.size - 1):
        idx = order[offs[c]: offs[c + 1]]
        o = idx[np.argsort(rt[idx], kind="stable")]
        new = np.concatenate([[True], np.diff(rt[o]) > bw])
        lab = np.cumsum(new) - 1
        group[o] = base + lab
        base += int(lab.max()) + 1
    split = np.cumsum([len(t) for t in tables])[:-1]
    return np.split(group, split)


__all__ = ["align", "align_dtw", "density_groups", "consensus_features", "Alignment", "DriftModel", "fit_drift", "select_drift",
           "find_anchors", "fit_mz_calibration", "dtw_path", "dtw_warp"]
