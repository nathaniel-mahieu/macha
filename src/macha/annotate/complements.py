"""Precursor mass of a fragment group from complementary ion pairs.

A precursor of neutral mass M breaks into b_i and y_(L-i), and for singly charged fragments
``mz(b_i) + mz(y_(L-i)) = M + 2 * proton``.  Among the fragments of one group, the pairwise sums therefore pile up
at one value, which must place the precursor inside the isolation window for some charge.  The pile-up height is
compared with the same search in windows the precursor cannot come from (the null), and the inferred mass can be
checked against MS1, which plays no part in grouping.
"""

from __future__ import annotations

import numpy as np

from .. import chem


def complement_mass(mz: np.ndarray, window: tuple, charges=(2, 3, 4), tol_ppm: float = 10.0,
                    iso_slack: int = 2) -> tuple[float, int, int]:
    """Best precursor neutral mass from pairwise fragment sums: (mass, charge, supporting pairs).

    A candidate must put the monoisotopic precursor m/z inside ``window`` (allowing it ``iso_slack`` isotope
    spacings below the window, where an isolated M+1 / M+2 would leave it).  Returns (nan, 0, 0) if none.
    """
    mz = np.sort(np.asarray(mz, float))
    if mz.size < 2:
        return float("nan"), 0, 0
    i, j = np.triu_indices(mz.size, 1)
    S = mz[i] + mz[j] - 2 * chem.PROTON  # candidate neutral masses
    lo, hi = window
    best = (float("nan"), 0, 0)
    for z in charges:
        pm = (S + z * chem.PROTON) / z
        ok = (pm >= lo - iso_slack * chem.NEUTRON_C13 / z) & (pm < hi)
        if ok.sum() < 2:
            continue
        x = np.sort(S[ok])
        # densest cluster within tol: two-pointer sweep
        tol = tol_ppm * 1e-6 * x
        r = np.searchsorted(x, x + tol, "right")
        cnt = r - np.arange(x.size)
        k = int(np.argmax(cnt))
        if cnt[k] > best[2]:
            best = (float(np.median(x[k:r[k]])), z, int(cnt[k]))
    return best


def complement_null(mz: np.ndarray, window: tuple, shifts=(-60.0, -40.0, 40.0, 60.0), **kw) -> int:
    """Largest pile-up found when the precursor is required to sit in shifted (wrong) windows."""
    return max(complement_mass(mz, (window[0] + s, window[1] + s), **kw)[2] for s in shifts)
