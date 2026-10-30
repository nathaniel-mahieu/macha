"""Per-edge evidence: mass residual, co-elution, profile shape, isotope ratio, cross-sample correlation."""

from __future__ import annotations

import numba as nb
import numpy as np
import pandas as pd

from .. import chem


@nb.njit(cache=True)
def _value(t, s, offs, first_scan, y):
    p = s - first_scan[t]
    if p < 0 or p >= offs[t + 1] - offs[t]:
        return 0.0
    return y[offs[t] + p]


@nb.njit(cache=True)
def pair_profile_corr(ti, tj, s0, s1, offs, first_scan, y):
    """Pearson correlation of two traces' baseline-subtracted intensities over scans [s0, s1]."""
    m = ti.shape[0]
    out = np.full(m, np.nan)
    for k in range(m):
        n = s1[k] - s0[k] + 1
        if n < 4:
            continue
        sa = 0.0
        sb = 0.0
        saa = 0.0
        sbb = 0.0
        sab = 0.0
        for s in range(s0[k], s1[k] + 1):
            a = _value(ti[k], s, offs, first_scan, y)
            b = _value(tj[k], s, offs, first_scan, y)
            sa += a
            sb += b
            saa += a * a
            sbb += b * b
            sab += a * b
        va = saa - sa * sa / n
        vb = sbb - sb * sb / n
        if va > 0 and vb > 0:
            out[k] = (sab - sa * sb / n) / np.sqrt(va * vb)
    return out


class ProfileSource:
    """Dense baseline-subtracted traces used for shape correlations."""

    def __init__(self, baselines):
        self.offs = baselines.offsets
        self.first_scan = baselines.scans[baselines.offsets[:-1]]
        y = np.where(np.isfinite(baselines.intensity), baselines.intensity - baselines.baseline, 0.0)
        self.y = np.maximum(y, 0.0)

    def corr(self, ti, tj, s0, s1):
        return pair_profile_corr(np.asarray(ti, np.int64), np.asarray(tj, np.int64), np.asarray(s0, np.int64),
                                 np.asarray(s1, np.int64), self.offs, self.first_scan, self.y)


def apex_sigma(features: pd.DataFrame, cycle: float) -> np.ndarray:
    """Apex-location uncertainty (s): width / sqrt(SNR), floored at a quarter cycle."""
    fw = features["fwhm_s"].to_numpy() if "fwhm_s" in features else np.full(len(features), 4 * cycle)
    snr = features["snr"].to_numpy() if "snr" in features else np.full(len(features), 20.0)
    return np.sqrt((fw / (2.355 * np.sqrt(np.maximum(snr, 1.0)))) ** 2 + (0.25 * cycle) ** 2)


def compute_evidence(edges: pd.DataFrame, features: pd.DataFrame, profiles: ProfileSource | None = None,
                     intensity_matrix: np.ndarray | None = None, compound_class: str = "metabolite",
                     cycle: float = 0.5) -> pd.DataFrame:
    e = edges.copy()
    if e.empty:
        for c in ("u", "dt", "rho", "ratio", "xcorr"):
            e[c] = []
        return e
    i = e["i"].to_numpy()
    j = e["j"].to_numpy()
    e["u"] = e["resid"] / e["sigma"]
    rt = features["rt"].to_numpy()
    sa = apex_sigma(features, cycle)
    e["dt"] = (rt[j] - rt[i]) / np.sqrt(sa[i] ** 2 + sa[j] ** 2)
    # profile correlation over the union of both peak windows
    if profiles is not None and {"trace", "scan_start", "scan_end"} <= set(features.columns):
        tr = features["trace"].to_numpy()
        s0 = np.minimum(features["scan_start"].to_numpy()[i], features["scan_start"].to_numpy()[j])
        s1 = np.maximum(features["scan_end"].to_numpy()[i], features["scan_end"].to_numpy()[j])
        e["rho"] = profiles.corr(tr[i], tr[j], s0, s1)
    else:
        e["rho"] = np.nan
    # 13C isotope intensity ratio vs expectation from mass (Poisson carbon model)
    ratio = np.full(len(e), np.nan)
    iso = (e["target"] == "13C").to_numpy() & (e["n"].to_numpy() <= 2)
    if iso.any() and "height" in features:
        h = features["height"].to_numpy()
        snr = features["snr"].to_numpy() if "snr" in features else np.full(len(features), 20.0)
        M = features["mz"].to_numpy()[i[iso]] * e["z"].to_numpy()[iso]
        nn = e["n"].to_numpy()[iso]
        expct = np.array([chem.expected_isotope_ratio(m, int(k), compound_class) for m, k in zip(M, nn)])
        obs = np.log(h[j[iso]] / h[i[iso]])
        sd = np.sqrt(0.25**2 + 1 / snr[i[iso]] ** 2 + 1 / snr[j[iso]] ** 2)
        ratio[iso] = (obs - np.log(expct)) / sd
    e["ratio"] = ratio
    if intensity_matrix is not None:
        X = np.log(np.maximum(np.asarray(intensity_matrix, float), 1e-9))
        ok = np.isfinite(np.asarray(intensity_matrix, float)) & (np.asarray(intensity_matrix, float) > 0)
        e["xcorr"] = _row_corr(X, ok, i, j)
    else:
        e["xcorr"] = np.nan
    return e


def _row_corr(X, ok, i, j):
    out = np.full(i.size, np.nan)
    for k in range(i.size):
        m = ok[i[k]] & ok[j[k]]
        if m.sum() >= 4:
            a, b = X[i[k], m], X[j[k], m]
            if a.std() > 0 and b.std() > 0:
                out[k] = np.corrcoef(a, b)[0, 1]
    return out
