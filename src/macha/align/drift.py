"""Drift-trend mixture model.

Each compound ``c`` belongs to a latent drift class ``k``.  Its retention
time in sample ``j`` is::

    t_cj = mu_c + w_kj(mu_c) + e,      e ~ (1 - eps) N(0, s^2) + eps U(outlier)

``w_kj`` is a P-spline (cubic B-splines, second-difference penalty) and the
warps of each class are centred across samples (``sum_j w_kj(t) = 0``) so
``mu_c`` is the cross-sample consensus time.  EM alternates class
responsibilities, penalised weighted least squares for every (class, sample)
warp, consensus times, the noise scale and class weights.  The number of
classes is chosen by BIC.

Rather than one warp per sample (global alignment) or one per m/z channel
(the original warpgroup/DTW approach), warps are shared by every compound
that drifts the same way, so a class with few anchors at some RT borrows
strength from the penalty and from its other compounds.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import BSpline
from scipy.special import logsumexp


def bspline_basis(t0: float, t1: float, n_basis: int, degree: int = 3):
    inner = np.linspace(t0, t1, n_basis - degree + 1)
    knots = np.concatenate([[t0] * degree, inner, [t1] * degree])
    return knots, degree


def design(x, knots, degree):
    x = np.clip(np.asarray(x, float), knots[0], knots[-1] - 1e-9)
    return BSpline.design_matrix(x, knots, degree).toarray()


def diff_penalty(n: int, order: int = 2):
    D = np.diff(np.eye(n), order, axis=0)
    return D.T @ D


@dataclass
class DriftModel:
    K: int
    knots: np.ndarray
    degree: int
    beta: np.ndarray  # (K, J, n_basis)
    pi: np.ndarray  # (K,)
    sigma: float
    eps: float
    loglik: float = np.nan
    bic: float = np.nan
    resp: np.ndarray | None = None  # anchors x K
    mu: np.ndarray | None = None  # anchor consensus times
    history: list = field(default_factory=list)

    @property
    def n_samples(self) -> int:
        return self.beta.shape[1]

    def warp(self, t, k: int, j: int):
        return design(t, self.knots, self.degree) @ self.beta[k, j]

    def warps(self, t):
        """All warps at times ``t``: array (K, J, len(t))."""
        B = design(t, self.knots, self.degree)
        return np.einsum("kjb,nb->kjn", self.beta, B)

    def to_reference(self, rt, k: int, j: int, n_iter: int = 4):
        """Invert t = mu + w_kj(mu) by fixed-point iteration."""
        mu = np.asarray(rt, float).copy()
        for _ in range(n_iter):
            mu = rt - self.warp(mu, k, j)
        return mu

    def grid(self, n: int = 2048):
        g = np.linspace(self.knots[0], self.knots[-1], n)
        return g, self.warps(g)


def _fit_em(T, K, knots, degree, lam, n_iter, init_resp, outlier_width, seed, tol=1e-6, S=None):
    C, J = T.shape
    obs = np.isfinite(T)
    S2 = np.zeros_like(T) if S is None else np.where(np.isfinite(S), S, 0.0) ** 2
    sigma_floor = 0.05  # s; keeps the scale identifiable when apex errors explain all scatter
    nb = len(knots) - degree - 1
    P = diff_penalty(nb)
    mu = np.nanmedian(T, axis=1)
    resp = init_resp
    beta = np.zeros((K, J, nb))
    sigma = max(float(np.nanstd(T - mu[:, None])), 0.5)
    eps = 0.02
    pi = resp.mean(0)
    lo_u = -np.log(outlier_width)
    out_p = np.zeros_like(T)  # outlier probability per anchor and sample
    history = []
    prev = -np.inf
    for it in range(n_iter):
        B = design(mu, knots, degree)
        # M-step: warps per class and sample.  A ridge (scale-aware, tiny where data exist) pulls each class warp
        # towards the run's class-weighted mean warp: where a class has no members it follows the other classes
        # (global drift is shared) instead of extrapolating a trend or falling back to no drift at all
        mean_beta = np.einsum("k,kjb->jb", pi, beta)
        for k in range(K):
            for j in range(J):
                m = obs[:, j]
                w = resp[m, k] * (1 - out_p[m, j]) / (sigma**2 + S2[m, j])
                Bm = B[m]
                ridge = 1e-2 * (np.median(w) if w.size else 1.0)
                A = Bm.T @ (Bm * w[:, None]) + lam * P + (ridge + 1e-8) * np.eye(nb)
                rhs = Bm.T @ (w * (T[m, j] - mu[m])) + ridge * mean_beta[j]
                beta[k, j] = np.linalg.solve(A, rhs)
            beta[k] -= beta[k].mean(0, keepdims=True)
        W = np.einsum("kjb,cb->ckj", beta, B)  # (C, K, J)
        # consensus times
        pred = np.einsum("ck,ckj->cj", resp, W)
        wmu = np.where(obs, (1 - out_p) / (sigma**2 + S2), 0.0)
        wmu = np.where(wmu.sum(1, keepdims=True) > 0, wmu, obs.astype(float))
        mu = np.where(obs, T - pred, 0.0).__mul__(wmu).sum(1) / wmu.sum(1)  # outlier-weighted consensus
        B = design(mu, knots, degree)
        W = np.einsum("kjb,cb->ckj", beta, B)
        R = np.where(obs[:, None, :], T[:, None, :] - mu[:, None, None] - W, 0.0)  # (C, K, J)
        # robust likelihood per observation: inlier gaussian vs uniform outlier
        V = sigma**2 + S2[:, None, :]
        lin = -0.5 * R**2 / V - 0.5 * np.log(V) - 0.5 * np.log(2 * np.pi) + np.log(1 - eps)
        lout = lo_u + np.log(eps)
        lobs = np.logaddexp(lin, lout)
        ll_ck = np.where(obs[:, None, :], lobs, 0.0).sum(2) + np.log(pi)[None, :]
        ll_c = logsumexp(ll_ck, axis=1)
        resp = np.exp(ll_ck - ll_c[:, None])
        # outlier probabilities (per anchor and sample, averaged over classes)
        pout = np.exp(lout - lobs)
        out_p = np.einsum("ck,ckj->cj", resp, pout) * obs
        win = np.einsum("ck,ckj->cj", resp, 1 - pout) * obs
        r2 = np.einsum("ck,ckj->cj", resp, (1 - pout) * R**2) * obs
        # noise scale by Fisher scoring on the inlier likelihood with variances sigma^2 + S^2: the moment estimate
        # sum(r^2 - S^2) / n goes negative (and pins sigma to its floor) when a few anchors carry a large S
        s2 = sigma**2
        for _ in range(30):
            V = s2 + S2
            step = float(((r2 - win * V) / V**2).sum() / max((win / V**2).sum(), 1e-12))
            s2 = max(s2 + step, sigma_floor**2)
            if abs(step) < 1e-6 * s2:
                break
        sigma = float(np.sqrt(s2))
        eps = float(np.clip(out_p.sum() / max(obs.sum(), 1), 1e-4, 0.3))
        pi = np.clip(resp.mean(0), 1e-6, None)
        pi /= pi.sum()
        ll = float(ll_c.sum())
        history.append(dict(iter=it, loglik=ll, sigma=sigma, eps=eps, pi=pi.tolist()))
        if abs(ll - prev) < tol * max(abs(ll), 1.0) and it > 3:
            break
        prev = ll
    return beta, pi, sigma, eps, ll, resp, mu, history


def _kmeans_init(D, K, seed, n_iter=50):
    rng = np.random.default_rng(seed)
    X = np.nan_to_num(D)
    C = X[rng.choice(len(X), K, replace=False)]
    for _ in range(n_iter):
        d = ((X[:, None, :] - C[None]) ** 2).sum(2)
        lab = d.argmin(1)
        for k in range(K):
            if (lab == k).any():
                C[k] = X[lab == k].mean(0)
    resp = np.full((len(X), K), 0.05 / max(K - 1, 1))
    resp[np.arange(len(X)), lab] = 0.95 if K > 1 else 1.0
    return resp / resp.sum(1, keepdims=True)


def _split_init(prev: DriftModel, T: np.ndarray, S=None):
    """Responsibilities for K+1 classes: split the worst-fitting class along its residuals' first PC."""
    K = prev.K
    B = design(prev.mu, prev.knots, prev.degree)
    lab = prev.resp.argmax(1)
    W = np.einsum("kjb,cb->ckj", prev.beta, B)
    R = T - prev.mu[:, None] - W[np.arange(len(T)), lab]
    R = np.where(np.isfinite(R), R, 0.0)
    sse = np.array([np.sum(R[lab == k] ** 2) for k in range(K)])
    k0 = int(np.argmax(sse))
    m = lab == k0
    resp = np.zeros((len(T), K + 1))
    resp[np.arange(len(T)), lab] = 1.0
    if m.sum() >= 4:
        Rk = R[m] - R[m].mean(0)
        v = np.linalg.svd(Rk, full_matrices=False)[2][0]
        side = Rk @ v > 0
        idx = np.flatnonzero(m)
        resp[idx[side], k0] = 0.0
        resp[idx[side], K] = 1.0
    resp = 0.9 * resp + 0.1 / (K + 1)
    return resp / resp.sum(1, keepdims=True)


def support_range(T: np.ndarray, q: float = 0.002, pad: float = 1.0):
    """Robust RT domain of a drift fit: quantiles of the consensus times (warps are constant beyond)."""
    mu = np.nanmedian(T, axis=1)
    mu = mu[np.isfinite(mu)]
    return float(np.quantile(mu, q)) - pad, float(np.quantile(mu, 1 - q)) + pad


def fit_drift(T: np.ndarray, K: int, t_range=None, n_basis: int = 10, lam: float = 1.0, n_iter: int = 200,
              n_init: int = 6, seed: int = 0, S: np.ndarray | None = None, init_resps=()) -> DriftModel:
    """Fit the K-class drift mixture to an anchor RT matrix ``T`` (anchors x samples, NaN = missing)."""
    T = np.asarray(T, float)
    C, J = T.shape
    t0, t1 = t_range if t_range is not None else support_range(T)
    knots, degree = bspline_basis(t0, t1, n_basis)
    span = float(np.nanmax(np.abs(T - np.nanmedian(T, 1, keepdims=True)))) * 4 + 10.0
    best = None
    # residual drift vectors after a one-class fit seed the class split
    base = _fit_em(T, 1, knots, degree, lam, n_iter, np.ones((C, 1)), span, seed, S=S)
    if K == 1:
        beta, pi, sigma, eps, ll, resp, mu, hist = base
        best = (ll, beta, pi, sigma, eps, resp, mu, hist)
    else:
        Bm = design(base[6], knots, degree)
        W1 = np.einsum("jb,cb->cj", base[0][0], Bm)
        D = T - base[6][:, None] - W1
        D = np.where(np.isfinite(D), D, 0.0)
        inits = [_kmeans_init(D, K, seed + s) for s in range(n_init)] + list(init_resps)
        for s, init in enumerate(inits):
            beta, pi, sigma, eps, ll, resp, mu, hist = _fit_em(T, K, knots, degree, lam, n_iter, init, span,
                                                               seed + s, S=S)
            if best is None or ll > best[0]:
                best = (ll, beta, pi, sigma, eps, resp, mu, hist)
    ll, beta, pi, sigma, eps, resp, mu, hist = best
    nb = len(knots) - degree - 1
    edf = _edf(nb, lam, T, mu, knots, degree)
    n_par = K * (J - 1) * edf + (K - 1) + 2 + C
    bic = -2 * ll + n_par * np.log(np.isfinite(T).sum())
    return DriftModel(K, knots, degree, beta, pi, sigma, eps, ll, bic, resp, mu, hist)


def _edf(nb, lam, T, mu, knots, degree):
    """Effective degrees of freedom of one warp (trace of the P-spline hat matrix)."""
    B = design(mu, knots, degree)
    P = diff_penalty(nb)
    BtB = B.T @ B
    A = BtB + lam * P + 1e-8 * np.eye(nb)
    return float(np.trace(np.linalg.solve(A, BtB)))


def select_drift(T: np.ndarray, K_max: int = 4, **kw) -> tuple[DriftModel, list[DriftModel]]:
    """Fit K = 1..K_max and return (best by BIC, all fits)."""
    fits = []
    for K in range(1, K_max + 1):
        if K > 1 and np.isfinite(T).all(axis=1).sum() < 10 * K:
            break
        extra = [_split_init(fits[-1], T)] if fits else []
        fits.append(fit_drift(T, K, init_resps=extra, **kw))
    best = min(fits, key=lambda m: m.bic)
    return best, fits


def cv_drift(T: np.ndarray, S: np.ndarray | None = None, K_max: int = 4, lams=(0.03, 0.3, 3.0), n_folds: int = 4,
             block: float = 10.0, seed: int = 0, max_anchors: int = 8000, **kw):
    """Choose the number of classes and the smoothing penalty by blocked cross-validation; refit on all anchors.

    The penalty is not scale-free: the data term grows with the number of anchors and their precision, so a fixed
    ``lam`` over-smooths a small anchor set and under-smooths a large one, and BIC over-counts anchors that are
    not independent (several ions of one compound share its RT error).  Anchors are assigned to RT blocks of
    ``block`` seconds (co-eluting siblings stay together) and blocks to folds; each held-out anchor's RT in every
    run is predicted from its other runs (:func:`macha.align.metrics.loo_prediction_error`).  Returns
    ``(model, cv)`` with ``cv`` the held-out RMSE per (K, lam).  The choice is made on at most ``max_anchors``
    anchors (a random subset); the final model uses all of them.
    """
    from .metrics import loo_prediction_error

    T_all = np.asarray(T, float)
    S_all = S
    t_range = kw.pop("t_range", None) or support_range(T_all)
    if T_all.shape[0] > max_anchors:
        sel = np.sort(np.random.default_rng(seed).choice(T_all.shape[0], max_anchors, replace=False))
        T, S = T_all[sel], (None if S_all is None else S_all[sel])
    else:
        T = T_all
    mu = np.nanmedian(T, axis=1)
    blocks = np.floor((mu - np.nanmin(mu)) / block).astype(np.int64)
    perm = np.random.default_rng(seed).permutation(blocks.max() + 1)
    fold = perm[blocks] % n_folds
    rows = []
    for K in range(1, K_max + 1):
        if K > 1 and np.isfinite(T).all(axis=1).sum() < 10 * K:
            break
        for lam in lams:
            err = []
            for f in range(n_folds):
                tr, te = fold != f, fold == f
                if te.sum() == 0 or tr.sum() < 10 * K:
                    continue
                m = fit_drift(T[tr], K, t_range=t_range, lam=lam, seed=seed, S=None if S is None else S[tr], **kw)
                err.append(loo_prediction_error(m, T[te]).ravel())
            e = np.concatenate(err) if err else np.array([np.nan])
            e = e[np.isfinite(e)]
            # a robust RMSE: gross mismatches would otherwise dominate the choice
            rows.append(dict(K=K, lam=lam, rmse=float(np.sqrt(np.mean(np.minimum(e * e, 25.0)))) if e.size else np.nan,
                             n=int(e.size)))
    import pandas as pd

    cv = pd.DataFrame(rows)
    best = cv.loc[cv["rmse"].idxmin()]
    model = fit_drift(T_all, int(best.K), t_range=t_range, lam=float(best.lam), seed=seed, S=S_all, **kw)
    return model, cv
