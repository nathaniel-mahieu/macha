"""Anchor selection and per-sample m/z calibration."""

from __future__ import annotations

from dataclasses import dataclass

import numba as nb
import numpy as np
import pandas as pd

from .drift import bspline_basis, design, diff_penalty


@nb.njit(cache=True)
def _unique_mask(mz, rt, tol_ppm, rt_window):
    """True where no other feature lies within tol_ppm and rt_window (mz sorted)."""
    n = mz.shape[0]
    out = np.ones(n, np.bool_)
    for i in range(n):
        tol = mz[i] * tol_ppm * 1e-6
        j = i + 1
        while j < n and mz[j] - mz[i] <= tol:
            if abs(rt[j] - rt[i]) <= rt_window:
                out[i] = False
                out[j] = False
            j += 1
    return out


def apex_se(fwhm_s, snr):
    """Apex-location standard error (s): ~ width / SNR, plus a small sampling floor."""
    fwhm_s = np.asarray(fwhm_s, float)
    snr = np.maximum(np.asarray(snr, float), 1.0)
    return np.sqrt((1.5 * fwhm_s / (2.355 * snr)) ** 2 + (0.02 * fwhm_s) ** 2)


@dataclass
class Anchors:
    rt: np.ndarray  # (C, J)
    mz: np.ndarray  # (C, J)
    feature: np.ndarray  # (C, J) feature index or -1
    se: np.ndarray | None = None  # (C, J) apex-location standard error


def find_anchors(tables: list[pd.DataFrame], mz_tol_ppm: float = 6.0, rt_window: float = 60.0,
                 min_frac: float = 0.6, min_snr: float = 10.0) -> Anchors:
    """Features that are unambiguous in m/z inside a wide RT window, matched across samples."""
    J = len(tables)
    parts = []
    for j, t in enumerate(tables):
        o = np.argsort(t["mz"].to_numpy(), kind="stable")
        mz = t["mz"].to_numpy()[o]
        rt = t["rt"].to_numpy()[o]
        uniq = _unique_mask(mz, rt, mz_tol_ppm, rt_window)
        snr = t["snr"].to_numpy()[o] if "snr" in t else np.full(len(t), np.inf)
        keep = uniq & (snr >= min_snr)
        if "fwhm_s" in t and "snr" in t:
            ase = apex_se(t["fwhm_s"].to_numpy()[o], snr)
        else:
            ase = np.zeros(len(t))
        parts.append(pd.DataFrame(dict(sample=j, feature=o[keep], mz=mz[keep], rt=rt[keep], se=ase[keep])))
    P = pd.concat(parts, ignore_index=True).sort_values("mz", kind="stable").reset_index(drop=True)
    mz = P["mz"].to_numpy()
    new = np.ones(len(P), bool)
    new[1:] = np.diff(mz) > mz[1:] * mz_tol_ppm * 1e-6
    cl = np.cumsum(new) - 1
    P["cluster"] = cl
    g = P.groupby("cluster")
    stats = g.agg(n=("sample", "size"), ns=("sample", "nunique"), rtspan=("rt", lambda x: x.max() - x.min()),
                  mzspan=("mz", lambda x: (x.max() - x.min()) / x.mean() * 1e6))
    good = stats[(stats.n == stats.ns) & (stats.ns >= max(2, int(np.ceil(min_frac * J)))) &
                 (stats.rtspan <= rt_window) & (stats.mzspan <= mz_tol_ppm)].index
    P = P[P["cluster"].isin(good)]
    ids = {c: i for i, c in enumerate(good)}
    C = len(good)
    RT = np.full((C, J), np.nan)
    MZ = np.full((C, J), np.nan)
    F = np.full((C, J), -1, np.int64)
    SE = np.full((C, J), np.nan)
    ci = P["cluster"].map(ids).to_numpy()
    sj = P["sample"].to_numpy()
    RT[ci, sj] = P["rt"].to_numpy()
    MZ[ci, sj] = P["mz"].to_numpy()
    F[ci, sj] = P["feature"].to_numpy()
    SE[ci, sj] = P["se"].to_numpy()
    return Anchors(RT, MZ, F, SE)


@dataclass
class MzCalibration:
    """Per-sample ppm offset as a smooth function of retention time."""

    knots: np.ndarray
    degree: int
    coef: np.ndarray  # (J, n_basis)

    def ppm(self, rt, j: int):
        return design(rt, self.knots, self.degree) @ self.coef[j]

    def correct(self, mz, rt, j: int):
        return np.asarray(mz) / (1.0 + 1e-6 * self.ppm(rt, j))


def fit_mz_calibration(anchors: Anchors, t_range, n_basis: int = 6, lam: float = 10.0, n_iter: int = 3) -> MzCalibration:
    knots, degree = bspline_basis(t_range[0], t_range[1], n_basis)
    MZ = anchors.mz
    ref = np.nanmedian(MZ, axis=1)
    D = (MZ - ref[:, None]) / ref[:, None] * 1e6
    J = MZ.shape[1]
    nb_ = len(knots) - degree - 1
    P = diff_penalty(nb_)
    coef = np.zeros((J, nb_))
    for j in range(J):
        m = np.isfinite(D[:, j])
        if m.sum() < 5:
            continue
        B = design(anchors.rt[m, j], knots, degree)
        y = D[m, j]
        w = np.ones(m.sum())
        for _ in range(n_iter):
            A = B.T @ (B * w[:, None]) + lam * P + 1e-6 * np.eye(nb_)
            c = np.linalg.solve(A, B.T @ (w * y))
            r = y - B @ c
            s = 1.4826 * np.median(np.abs(r)) + 1e-6
            w = (np.abs(r) < 3 * s).astype(float)
        coef[j] = c
    coef -= coef.mean(0, keepdims=True)
    return MzCalibration(knots, degree, coef)
