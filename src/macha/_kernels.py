"""Small numba kernels shared across modules."""

from __future__ import annotations

import numba as nb
import numpy as np


@nb.njit(cache=True)
def scanwise_searchsorted(mz, offsets, value, left):
    """For each scan, the absolute index where ``value`` would be inserted."""
    n = offsets.shape[0] - 1
    out = np.empty(n, dtype=np.int64)
    for s in range(n):
        lo, hi = offsets[s], offsets[s + 1]
        while lo < hi:
            mid = (lo + hi) // 2
            if (mz[mid] < value) if left else (mz[mid] <= value):
                lo = mid + 1
            else:
                hi = mid
        out[s] = lo
    return out


@nb.njit(cache=True)
def group_offsets(sorted_labels, n_groups):
    """CSR offsets for an array of labels that is already sorted."""
    off = np.zeros(n_groups + 1, dtype=np.int64)
    for v in sorted_labels:
        off[v + 1] += 1
    for g in range(n_groups):
        off[g + 1] += off[g]
    return off


def csr_from_labels(labels: np.ndarray, n_groups: int | None = None):
    """Return (order, offsets) so that ``order[offsets[g]:offsets[g+1]]`` lists members of g.

    Negative labels are dropped.
    """
    labels = np.asarray(labels)
    valid = np.flatnonzero(labels >= 0)
    if n_groups is None:
        n_groups = int(labels.max()) + 1 if valid.size else 0
    order = valid[np.argsort(labels[valid], kind="stable")]
    return order, group_offsets(labels[order].astype(np.int64), n_groups)
