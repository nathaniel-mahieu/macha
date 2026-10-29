"""Two-groups empirical-Bayes evaluation of candidate edges.

Each candidate edge ``e`` (a co-eluting pair whose mass difference matches a
rule) is either a true relation (H1) or a coincidence (H0)::

    p(x_e) = pi_g f1(x_e) + (1 - pi_g) f0(x_e)

* ``f0`` is learned non-parametrically from **decoy edges**: the same search
  run with rule masses shifted to values no chemistry produces.  Decoys are
  coincidences by construction, so their evidence distribution is the null.
  The co-elution and shape evidence are modelled jointly under H0 because
  coincidentally co-eluting peaks are correlated more often the closer they
  are.
* ``f1`` is parametric: mass residual ~ N(0, s_u) truncated to the gate,
  |RT difference| ~ half-normal, Fisher-z of profile correlation ~ normal
  with an SNR-dependent mean, standardised 13C ratio residual ~ normal,
  Fisher-z of cross-sample correlation ~ normal.  Shape and ratio terms are
  robust mixtures ``(1 - eta) N + eta f0``: a true relation whose trace is
  contaminated by a co-eluting interferent looks like the null in that
  dimension, so it becomes uninformative rather than damning.
* ``pi_g`` is a prior per (rule, charge) with a Beta(1, 1) hyperprior.

EM estimates ``pi`` and the ``f1`` parameters with ``f0`` fixed.  The
posterior and log-odds of every edge follow.  Applying the fitted model to
the decoys gives an independent FDR estimate to check calibration against
the model's own (Bayesian) FDR.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.special import expit, ndtr

LOG_SQRT_2PI = 0.5 * np.log(2 * np.pi)


def _fz(r):
    return np.arctanh(np.clip(r, -0.995, 0.995))


class _Hist1D:
    """Smoothed histogram density shrunk toward a weak parametric prior (``prior_n`` pseudo-counts)."""

    def __init__(self, x, lo, hi, bins=30, smooth=1.0, prior=None, prior_n=40.0):
        x = x[np.isfinite(x)]
        self.edges = np.linspace(lo, hi, bins + 1)
        c, _ = np.histogram(np.clip(x, lo, hi), self.edges)
        c = gaussian_filter(c.astype(float), smooth, mode="nearest")
        mid = 0.5 * (self.edges[1:] + self.edges[:-1])
        w = np.diff(self.edges)
        p = np.ones(bins) if prior is None else prior(mid)
        c = c + prior_n * p * w / np.sum(p * w) + 1e-3
        self.logd = np.log(c / c.sum() / w)

    def __call__(self, x):
        k = np.clip(np.searchsorted(self.edges, x, "right") - 1, 0, self.logd.size - 1)
        return self.logd[k]


class _Hist2D:
    def __init__(self, x, y, xr, yr, bins=(14, 14), smooth=1.0, prior=None, prior_n=40.0):
        ok = np.isfinite(x) & np.isfinite(y)
        self.ex = np.linspace(*xr, bins[0] + 1)
        self.ey = np.linspace(*yr, bins[1] + 1)
        c, _, _ = np.histogram2d(np.clip(x[ok], *xr), np.clip(y[ok], *yr), [self.ex, self.ey])
        c = gaussian_filter(c, smooth, mode="nearest")
        area = np.outer(np.diff(self.ex), np.diff(self.ey))
        mx = 0.5 * (self.ex[1:] + self.ex[:-1])
        my = 0.5 * (self.ey[1:] + self.ey[:-1])
        p = np.ones(area.shape) if prior is None else prior(mx[:, None], my[None, :])
        c = c + prior_n * p * area / np.sum(p * area) + 1e-3
        self.logd = np.log(c / c.sum() / area)

    def __call__(self, x, y):
        kx = np.clip(np.searchsorted(self.ex, x, "right") - 1, 0, self.logd.shape[0] - 1)
        ky = np.clip(np.searchsorted(self.ey, y, "right") - 1, 0, self.logd.shape[1] - 1)
        return self.logd[kx, ky]


@dataclass
class EdgeModel:
    gate: float = 4.0
    s_u: float = 1.0
    s_t: float = 1.0
    rho_coef: np.ndarray = field(default_factory=lambda: np.array([1.2, 0.2]))
    s_r: float = 0.6
    m_q: float = 0.0
    s_q: float = 1.0
    m_x: float = 1.0
    s_x: float = 0.6
    eta_r: float = 0.1  # share of true edges whose profiles are contaminated (shape uninformative)
    eta_q: float = 0.1  # same for isotope ratios (interference on one isotopologue)
    pi: dict = field(default_factory=dict)
    history: list = field(default_factory=list)
    use_rho: bool = False
    use_x: bool = False

    # ------------------------------------------------------------ null densities
    def _fit_null(self, d: pd.DataFrame):
        adt = np.abs(d["dt"].to_numpy())
        self.dt_max = float(max(np.nanquantile(np.abs(d["dt"]), 0.995), 1.0)) if len(d) else 5.0
        # weak parametric nulls: coincidental co-elution ~ uniform within the window, unrelated
        # profiles ~ weakly positively correlated, isotope-ratio residuals and cross-sample
        # correlations broad.  Decoy counts dominate once there are more than ~prior_n of them.
        self.f0_dt = _Hist1D(adt, 0.0, self.dt_max, bins=24)
        null_rho = lambda z: np.exp(-0.5 * ((z - 0.3) / 0.7) ** 2)  # noqa: E731
        self.use_rho = bool(np.isfinite(d["rho"]).sum() >= 10)
        if self.use_rho:
            self.f0_dt_rho = _Hist2D(adt, _fz(d["rho"].to_numpy()), (0.0, self.dt_max), (-1.5, 3.0),
                                     prior=lambda x, z: np.ones_like(x) * null_rho(z))
            self.f0_rho = _Hist1D(_fz(d["rho"].to_numpy()), -1.5, 3.0, bins=24, prior=null_rho)
        q = d.loc[np.isfinite(d["ratio"]), "ratio"].to_numpy()
        self.use_q = q.size >= 10
        if self.use_q:
            self.f0_q = _Hist1D(q, -8.0, 8.0, bins=32, prior=lambda v: np.exp(-0.5 * (v / 3.0) ** 2))
        self.use_x = bool(np.isfinite(d["xcorr"]).sum() >= 10)
        if self.use_x:
            self.f0_x = _Hist1D(_fz(d["xcorr"].to_numpy()), -2.0, 3.0, bins=25,
                                prior=lambda v: np.exp(-0.5 * (v / 0.6) ** 2))

    # ------------------------------------------------------------ likelihood ratio
    def log_lr(self, e: pd.DataFrame) -> np.ndarray:
        u = e["u"].to_numpy()
        g = self.gate
        lr = (-0.5 * (u / self.s_u) ** 2 - LOG_SQRT_2PI - np.log(self.s_u)
              - np.log(ndtr(g / self.s_u) - ndtr(-g / self.s_u))) + np.log(2 * g)
        adt = np.abs(e["dt"].to_numpy())
        f1_dt = np.log(2.0) - 0.5 * (adt / self.s_t) ** 2 - LOG_SQRT_2PI - np.log(self.s_t)
        rho = e["rho"].to_numpy() if "rho" in e else np.full(len(e), np.nan)
        have_rho = np.isfinite(rho) & self.use_rho
        z = _fz(np.where(have_rho, rho, 0.0))
        snr = e["snr_min"].to_numpy() if "snr_min" in e else np.full(len(e), 20.0)
        m_r = self.rho_coef[0] + self.rho_coef[1] * np.log(np.maximum(snr, 1.0))
        f1_r = -0.5 * ((z - m_r) / self.s_r) ** 2 - LOG_SQRT_2PI - np.log(self.s_r)
        if self.use_rho:
            # robust H1: a contaminated true profile is distributed like the null (uninformative)
            f1_r = np.logaddexp(np.log1p(-self.eta_r) + f1_r, np.log(self.eta_r) + self.f0_rho(z))
            joint0 = self.f0_dt_rho(adt, z)
            lr += np.where(have_rho, f1_dt + f1_r - joint0, f1_dt - self.f0_dt(adt))
        else:
            lr += f1_dt - self.f0_dt(adt)
        if self.use_q:
            q = e["ratio"].to_numpy()
            hq = np.isfinite(q)
            qq = np.where(hq, q, 0.0)
            l1 = -0.5 * ((qq - self.m_q) / self.s_q) ** 2 - LOG_SQRT_2PI - np.log(self.s_q)
            l0 = self.f0_q(qq)
            lq = np.logaddexp(np.log1p(-self.eta_q) + l1, np.log(self.eta_q) + l0) - l0
            lr += np.where(hq, lq, 0.0)
        if self.use_x:
            x = e["xcorr"].to_numpy()
            hx = np.isfinite(x)
            zx = _fz(np.where(hx, x, 0.0))
            lx = -0.5 * ((zx - self.m_x) / self.s_x) ** 2 - LOG_SQRT_2PI - np.log(self.s_x) - self.f0_x(zx)
            lr += np.where(hx, lx, 0.0)
        return lr

    def _group(self, e):
        return (e["target"].astype(str) + "|z" + e["z"].astype(str)).to_numpy()

    def log_prior_odds(self, e: pd.DataFrame) -> np.ndarray:
        g = self._group(e)
        pi = np.array([self.pi.get(k, self.pi_default) for k in g])
        return np.log(pi) - np.log1p(-pi)

    def log_odds(self, e: pd.DataFrame) -> np.ndarray:
        return self.log_prior_odds(e) + self.log_lr(e)

    # ------------------------------------------------------------ EM
    def fit(self, targets: pd.DataFrame, decoys: pd.DataFrame, max_iter: int = 200, tol: float = 1e-5):
        self._fit_null(decoys)
        t = targets
        groups = self._group(t)
        ug, ginv = np.unique(groups, return_inverse=True)
        self.pi = {k: 0.3 for k in ug}
        self.pi_default = 0.05
        u = t["u"].to_numpy()
        adt = np.abs(t["dt"].to_numpy())
        rho = t["rho"].to_numpy()
        hr = np.isfinite(rho) & self.use_rho
        zr = _fz(np.where(hr, rho, 0.0))
        snr = t["snr_min"].to_numpy() if "snr_min" in t else np.full(len(t), 20.0)
        ls = np.log(np.maximum(snr, 1.0))
        q = t["ratio"].to_numpy()
        hq = np.isfinite(q)
        x = t["xcorr"].to_numpy()
        hx = np.isfinite(x) & self.use_x
        zx = _fz(np.where(hx, x, 0.0))
        gamma = np.full(len(t), 0.5)
        for it in range(max_iter):
            lo = self.log_odds(t)
            new = expit(lo)
            delta = float(np.max(np.abs(new - gamma))) if len(t) else 0.0
            gamma = new
            w = gamma
            sw = w.sum() + 1e-12
            self.s_u = float(np.clip(np.sqrt((w * u**2).sum() / sw), 0.3, self.gate))
            self.s_t = float(np.clip(np.sqrt((w * adt**2).sum() / sw), 0.15, 10.0))
            if hr.sum() > 10:
                A = np.column_stack([np.ones(hr.sum()), ls[hr]])
                lN = -0.5 * ((zr[hr] - A @ self.rho_coef) / self.s_r) ** 2 - np.log(self.s_r) - LOG_SQRT_2PI
                l0 = self.f0_rho(zr[hr])
                clean = expit(np.log1p(-self.eta_r) + lN - np.log(self.eta_r) - l0)
                W = w[hr] * clean
                coef = np.linalg.lstsq(A * np.sqrt(W)[:, None], zr[hr] * np.sqrt(W), rcond=None)[0]
                self.rho_coef = coef
                res = zr[hr] - A @ coef
                self.s_r = float(np.clip(np.sqrt((W * res**2).sum() / (W.sum() + 1e-12)), 0.15, 2.0))
                self.eta_r = float(np.clip(1 - W.sum() / (w[hr].sum() + 1e-12), 0.01, 0.5))
            if self.use_q and hq.sum() > 10:
                l1 = -0.5 * ((q[hq] - self.m_q) / self.s_q) ** 2 - np.log(self.s_q) - LOG_SQRT_2PI
                l0 = self.f0_q(q[hq])
                clean = expit(np.log1p(-self.eta_q) + l1 - np.log(self.eta_q) - l0)
                W = w[hq] * clean
                self.m_q = float((W * q[hq]).sum() / (W.sum() + 1e-12))
                self.s_q = float(np.clip(np.sqrt((W * (q[hq] - self.m_q) ** 2).sum() / (W.sum() + 1e-12)), 0.3, 5))
                self.eta_q = float(np.clip(1 - W.sum() / (w[hq].sum() + 1e-12), 0.01, 0.5))
            if hx.sum() > 10:
                W = w[hx]
                self.m_x = float((W * zx[hx]).sum() / (W.sum() + 1e-12))
                self.s_x = float(np.clip(np.sqrt((W * (zx[hx] - self.m_x) ** 2).sum() / (W.sum() + 1e-12)), 0.15, 2))
            sums = np.bincount(ginv, weights=w, minlength=ug.size)
            cnts = np.bincount(ginv, minlength=ug.size)
            self.pi = {k: float(np.clip((s + 1) / (c + 2), 1e-4, 1 - 1e-4)) for k, s, c in zip(ug, sums, cnts)}
            self.history.append(dict(iter=it, delta=delta, s_u=self.s_u, s_t=self.s_t, s_r=self.s_r,
                                     n_true=float(sw)))
            if delta < tol:
                break
        return self

    def summary(self) -> dict:
        return dict(s_u=self.s_u, s_t=self.s_t, rho_coef=list(map(float, self.rho_coef)), s_r=self.s_r,
                    eta_rho=self.eta_r, eta_ratio=self.eta_q,
                    m_q=self.m_q, s_q=self.s_q, m_x=self.m_x, s_x=self.s_x, use_rho=self.use_rho,
                    use_ratio=self.use_q, use_xcorr=self.use_x, pi=self.pi, iterations=len(self.history))


def fdr_curves(targets: pd.DataFrame, decoys: pd.DataFrame, n_decoy_per_target: float, thresholds=None):
    """Bayesian FDR (mean 1-posterior of accepted targets) vs decoy-estimated FDR."""
    thresholds = np.linspace(0.05, 0.99, 40) if thresholds is None else thresholds
    rows = []
    pt = targets["posterior"].to_numpy()
    pdc = decoys["posterior"].to_numpy()
    for th in thresholds:
        acc = pt >= th
        n = acc.sum()
        bay = float(np.mean(1 - pt[acc])) if n else np.nan
        dec = float((pdc >= th).sum() / n_decoy_per_target / n) if n else np.nan
        rows.append(dict(threshold=th, accepted=int(n), bayes_fdr=bay, decoy_fdr=dec))
    return pd.DataFrame(rows)
