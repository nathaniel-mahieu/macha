"""Score trace assignments against synthetic ground truth."""

from __future__ import annotations

import numpy as np
from scipy import sparse


def _contingency(pred, truth):
    pu, pi = np.unique(pred, return_inverse=True)
    tu, ti = np.unique(truth, return_inverse=True)
    m = sparse.coo_matrix((np.ones(pred.size), (pi, ti)), shape=(pu.size, tu.size)).tocsr()
    return m, pu, tu


def adjusted_rand(pred, truth) -> float:
    m, _, _ = _contingency(pred, truth)
    n = pred.size

    def c2(x):
        x = np.asarray(x, dtype=np.float64)
        return float((x * (x - 1) / 2).sum())

    sij = c2(m.data)
    a = c2(np.asarray(m.sum(1)).ravel())
    b = c2(np.asarray(m.sum(0)).ravel())
    exp = a * b / (n * (n - 1) / 2)
    mx = 0.5 * (a + b)
    return float((sij - exp) / (mx - exp)) if mx != exp else 1.0


def trace_metrics(label: np.ndarray, ion: np.ndarray, min_points: int = 5) -> dict:
    """Purity / completeness / fragmentation of traces against true ion labels.

    ``label`` is the predicted trace per centroid (-1 = unassigned) and
    ``ion`` the true ion (-1 = noise).  Only traces and ions with at least
    ``min_points`` centroids are scored.
    """
    label = np.asarray(label)
    ion = np.asarray(ion)
    valid = label >= 0
    pc = np.bincount(label[valid]) if valid.any() else np.zeros(0, int)
    big_trace = np.zeros(label.size, bool)
    big_trace[valid] = pc[label[valid]] >= min_points

    # purity of predicted traces (noise counts as impurity)
    lt = label[big_trace]
    it = ion[big_trace]
    m, pu, tu = _contingency(lt, it)
    m = m.tocsr()
    noise_col = np.flatnonzero(tu == -1)
    row_max = np.zeros(m.shape[0])
    second = np.zeros(m.shape[0])
    for r in range(m.shape[0]):
        row = m.getrow(r)
        data = row.data.copy()
        if noise_col.size:
            cols = row.indices
            data = np.where(cols == noise_col[0], 0, data)
        s = np.sort(data)[::-1]
        row_max[r] = s[0] if s.size else 0
        second[r] = s[1] if s.size > 1 else 0
    sizes = np.asarray(m.sum(1)).ravel()
    purity = float(row_max.sum() / max(sizes.sum(), 1))
    merged = float(np.mean(second >= 0.1 * sizes)) if sizes.size else 0.0

    # completeness of true ions
    sig = ion >= 0
    ic = np.bincount(ion[sig]) if sig.any() else np.zeros(0, int)
    big_ion = np.zeros(ion.size, bool)
    big_ion[sig] = ic[ion[sig]] >= min_points
    li = np.where(big_trace, label, -1)[big_ion]
    ii = ion[big_ion]
    m2, pu2, tu2 = _contingency(ii, li)  # rows ions, cols traces (-1 = unassigned/short)
    m2 = m2.tocsr()
    unassigned = np.flatnonzero(tu2 == -1)
    best = np.zeros(m2.shape[0])
    frags = np.zeros(m2.shape[0])
    for r in range(m2.shape[0]):
        row = m2.getrow(r)
        data = row.data
        cols = row.indices
        if unassigned.size:
            data = np.where(cols == unassigned[0], 0, data)
        tot = row.data.sum()
        best[r] = data.max() / tot if tot else 0
        frags[r] = np.sum(data >= 0.1 * tot)
    isz = np.asarray(m2.sum(1)).ravel()
    completeness = float((best * isz).sum() / max(isz.sum(), 1))
    recall = float(np.mean(best >= 0.5)) if best.size else 0.0

    noise = ion == -1
    noise_in_traces = float(big_trace[noise].mean()) if noise.any() else 0.0
    keep = big_trace & sig
    ari = adjusted_rand(label[keep], ion[keep]) if keep.sum() > 1 else float("nan")
    return dict(
        n_traces=int(pu.size), purity=purity, completeness=completeness, ion_recall=recall,
        fragments_per_ion=float(frags.mean()) if frags.size else 0.0, merged_trace_frac=merged,
        noise_in_traces=noise_in_traces, ari=ari,
    )
