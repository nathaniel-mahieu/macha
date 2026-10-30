"""Candidate edge generation.

Modulo engine
    For an additive rule with step ``delta`` and charge ``z``, features with
    charge-scaled masses ``M = mz * z`` related by ``n`` steps satisfy
    ``M_j mod delta == M_i mod delta`` and ``n = round((M_j - M_i) / delta)``.
    Series rules (homologues, which do not co-elute) sort by residue and sweep
    the circular residue axis, so a whole series falls into one residue bin.
    Co-eluting rules use the same residue test; ``method="sweep"`` scans each
    RT window pairwise (reference implementation) while the default
    ``"targeted"`` pre-sorts each RT block by m/z and binary-searches the
    ``n * delta / z`` targets, which is far faster on dense data and returns
    the identical edge set.

Neutral-mass engine
    For pairs of adduct hypotheses that are not a constant mass apart (charge
    states, multimers) both features are mapped to a neutral mass and matched.
"""

from __future__ import annotations

import numba as nb
import numpy as np
import pandas as pd

from .rules import PairRule, Rule, affine


@nb.njit(cache=True)
def _modulo_coelute(mz, se, rt, order, z, delta, n_max, sigma_extra, sys_ppm, gate, rt_window):
    n_feat = mz.shape[0]
    oi = []
    oj = []
    on = []
    oe = []
    os_ = []
    lo = 0
    for qi in range(n_feat):
        i = order[qi]
        while rt[order[lo]] < rt[i] - rt_window:
            lo += 1
        Mi = mz[i] * z
        for qj in range(lo, n_feat):
            j = order[qj]
            if rt[j] > rt[i] + rt_window:
                break
            Mj = mz[j] * z
            d = Mj - Mi
            if d <= 0.0:
                continue
            n = int(np.floor(d / delta + 0.5))
            if n < 1 or n > n_max:
                continue
            e = d - n * delta
            sig = np.sqrt(z * z * (se[i] ** 2 + se[j] ** 2) + (sys_ppm * 1e-6 * Mj) ** 2 + n * sigma_extra**2)
            if abs(e) <= gate * sig:
                oi.append(i)
                oj.append(j)
                on.append(n)
                oe.append(e)
                os_.append(sig)
    return (np.array(oi, np.int64), np.array(oj, np.int64), np.array(on, np.int64), np.array(oe), np.array(os_))


@nb.njit(cache=True)
def _modulo_residue(mz, se, z, delta, n_max, sigma_extra, sys_ppm, gate):
    """Series search: sweep features sorted by residue ``M mod delta`` (circular)."""
    n_feat = mz.shape[0]
    M = mz * z
    r = M - np.floor(M / delta) * delta
    order = np.argsort(r)
    tol = 0.0
    for i in range(n_feat):
        t = gate * np.sqrt(2 * z * z * se[i] ** 2 + (sys_ppm * 1e-6 * M[i]) ** 2 + n_max * sigma_extra**2)
        if t > tol:
            tol = t
    oi = []
    oj = []
    on = []
    oe = []
    os_ = []
    for qa in range(n_feat):
        a = order[qa]
        for step in range(1, n_feat):
            qb = qa + step
            wrap = 0.0
            if qb >= n_feat:
                qb -= n_feat
                wrap = delta
            if qb == qa:
                break
            b = order[qb]
            if r[b] + wrap - r[a] > tol:
                break
            i, j = (a, b) if M[a] < M[b] else (b, a)
            d = M[j] - M[i]
            n = int(np.floor(d / delta + 0.5))
            if n < 1 or n > n_max:
                continue
            e = d - n * delta
            sig = np.sqrt(z * z * (se[i] ** 2 + se[j] ** 2) + (sys_ppm * 1e-6 * M[j]) ** 2 + n * sigma_extra**2)
            if abs(e) <= gate * sig:
                oi.append(i)
                oj.append(j)
                on.append(n)
                oe.append(e)
                os_.append(sig)
    return (np.array(oi, np.int64), np.array(oj, np.int64), np.array(on, np.int64), np.array(oe), np.array(os_))


@nb.njit(cache=True)
def _neutral_pairs(mz, se, rt, order, Aa, Ba, na, Ab, Bb, nb_, shift, sys_ppm, gate, rt_window):
    n_feat = mz.shape[0]
    oi = []
    oj = []
    oe = []
    os_ = []
    lo = 0
    for qi in range(n_feat):
        i = order[qi]
        while rt[order[lo]] < rt[i] - rt_window:
            lo += 1
        Ma = (mz[i] * Aa + Ba) / na
        for qj in range(lo, n_feat):
            j = order[qj]
            if rt[j] > rt[i] + rt_window:
                break
            if j == i:
                continue
            Mb = (mz[j] * Ab + Bb) / nb_
            e = Mb - Ma - shift
            sig = np.sqrt((Aa * se[i] / na) ** 2 + (Ab * se[j] / nb_) ** 2 + (sys_ppm * 1e-6 * Ma) ** 2)
            if abs(e) <= gate * sig:
                oi.append(i)
                oj.append(j)
                oe.append(e)
                os_.append(sig)
    return np.array(oi, np.int64), np.array(oj, np.int64), np.array(oe), np.array(os_)


def block_index(mz, rt, width):
    """RT blocks of ``width``; each block lists the features of itself and both neighbours, sorted by m/z."""
    blk = np.floor((rt - rt.min()) / width).astype(np.int64) if rt.size else np.zeros(0, np.int64)
    n_blk = int(blk.max()) + 1 if rt.size else 0
    key = np.concatenate([blk - 1, blk, blk + 1])
    idx = np.tile(np.arange(rt.size), 3)
    ok = (key >= 0) & (key < n_blk)
    key, idx = key[ok], idx[ok]
    o = np.lexsort((mz[idx], key))
    key, idx = key[o], idx[o]
    off = np.searchsorted(key, np.arange(n_blk + 1))
    return blk, mz[idx].copy(), idx.astype(np.int64), off.astype(np.int64)


@nb.njit(cache=True)
def _targeted_additive(mz, se, rt, blk, seg_mz, seg_idx, seg_off, z, delta, n_max, sigma_extra, sys_ppm, gate,
                       rt_window, se_max):
    """Same edges as :func:`_modulo_coelute` (canonical n only), via binary search of target masses."""
    oi = []
    oj = []
    on = []
    oe = []
    os_ = []
    for i in range(mz.shape[0]):
        a0, a1 = seg_off[blk[i]], seg_off[blk[i] + 1]
        Mi = mz[i] * z
        for n in range(1, n_max + 1):
            a, b = n, int(z)
            while b:
                a, b = b, a % b
            if a != 1:
                continue
            target = mz[i] + n * delta / z
            tol = gate * np.sqrt(se[i] ** 2 + se_max**2 + (sys_ppm * 1e-6 * target) ** 2
                                 + n * sigma_extra**2 / (z * z))
            q = a0 + np.searchsorted(seg_mz[a0:a1], target - tol)
            while q < a1 and seg_mz[q] <= target + tol:
                j = seg_idx[q]
                q += 1
                if abs(rt[j] - rt[i]) > rt_window:
                    continue
                Mj = mz[j] * z
                e = Mj - Mi - n * delta
                sig = np.sqrt(z * z * (se[i] ** 2 + se[j] ** 2) + (sys_ppm * 1e-6 * Mj) ** 2 + n * sigma_extra**2)
                if abs(e) <= gate * sig:
                    oi.append(i)
                    oj.append(j)
                    on.append(n)
                    oe.append(e)
                    os_.append(sig)
    return (np.array(oi, np.int64), np.array(oj, np.int64), np.array(on, np.int64), np.array(oe), np.array(os_))


@nb.njit(cache=True)
def _targeted_pairs(mz, se, rt, blk, seg_mz, seg_idx, seg_off, Aa, Ba, na, Ab, Bb, nb_, shift, sys_ppm, gate,
                    rt_window, se_max):
    oi = []
    oj = []
    oe = []
    os_ = []
    for i in range(mz.shape[0]):
        a0, a1 = seg_off[blk[i]], seg_off[blk[i] + 1]
        Ma = (mz[i] * Aa + Ba) / na
        target = ((Ma + shift) * nb_ - Bb) / Ab
        sig_max = np.sqrt((Aa * se[i] / na) ** 2 + (Ab * se_max / nb_) ** 2 + (sys_ppm * 1e-6 * Ma) ** 2)
        tol = gate * sig_max * nb_ / Ab
        q = a0 + np.searchsorted(seg_mz[a0:a1], target - tol)
        while q < a1 and seg_mz[q] <= target + tol:
            j = seg_idx[q]
            q += 1
            if j == i or abs(rt[j] - rt[i]) > rt_window:
                continue
            Mb = (mz[j] * Ab + Bb) / nb_
            e = Mb - Ma - shift
            sig = np.sqrt((Aa * se[i] / na) ** 2 + (Ab * se[j] / nb_) ** 2 + (sys_ppm * 1e-6 * Ma) ** 2)
            if abs(e) <= gate * sig:
                oi.append(i)
                oj.append(j)
                oe.append(e)
                os_.append(sig)
    return np.array(oi, np.int64), np.array(oj, np.int64), np.array(oe), np.array(os_)


def candidate_edges(features: pd.DataFrame, rules: list[Rule], pair_rules: list[PairRule], gate: float = 4.0,
                    sys_ppm: float = 0.3, rt_window: float | None = None, method: str = "targeted") -> pd.DataFrame:
    """All rule-consistent co-eluting feature pairs (targets and decoys alike).

    ``features`` needs columns ``mz``, ``mz_se``, ``rt`` and ``fwhm_s``.  The RT
    window defaults to the median FWHM (co-elution is a necessary condition).
    """
    mz = features["mz"].to_numpy(np.float64)
    se = features["mz_se"].to_numpy(np.float64) if "mz_se" in features else mz * 1e-6
    rt = features["rt"].to_numpy(np.float64)
    if rt_window is None:
        rt_window = float(np.median(features["fwhm_s"])) if "fwhm_s" in features and len(features) else 5.0
    order = np.argsort(rt, kind="stable")
    se_max = float(se.max()) if se.size else 0.0
    if method == "targeted" and mz.size:
        blk, seg_mz, seg_idx, seg_off = block_index(mz, rt, rt_window)
    frames = []
    for r in rules:
        for z in range(1, r.z_max + 1):
            if r.coelute and method == "targeted":
                i, j, n, e, s = _targeted_additive(mz, se, rt, blk, seg_mz, seg_idx, seg_off, float(z), r.delta,
                                                   r.n_max, r.sigma_extra, sys_ppm, gate, rt_window, se_max)
            elif r.coelute:
                i, j, n, e, s = _modulo_coelute(mz, se, rt, order, float(z), r.delta, r.n_max, r.sigma_extra,
                                                sys_ppm, gate, rt_window)
            else:
                i, j, n, e, s = _modulo_residue(mz, se, float(z), r.delta, r.n_max, r.sigma_extra, sys_ppm, gate)
            # an edge n*delta/z with gcd(n, z) = g > 1 is the same m/z spacing as (n/g, z/g):
            # keep only the canonical form, so a higher charge needs its own half-spaced peaks
            canon = np.gcd(n, z) == 1
            i, j, n, e, s = i[canon], j[canon], n[canon], e[canon], s[canon]
            if i.size:
                frames.append(pd.DataFrame(dict(i=i, j=j, rule=r.name, kind=r.kind, z=z, n=n, resid=e, sigma=s,
                                                decoy=r.decoy, target=r.target or r.name,
                                                nominal=n * int(round(r.delta)))))
    for p in pair_rules:
        Aa, Ba, na = affine(p.a)
        Ab, Bb, nb_ = affine(p.b)
        if method == "targeted":
            i, j, e, s = _targeted_pairs(mz, se, rt, blk, seg_mz, seg_idx, seg_off, Aa, Ba, na, Ab, Bb, nb_,
                                         p.decoy_shift, sys_ppm, gate, rt_window, se_max)
        else:
            i, j, e, s = _neutral_pairs(mz, se, rt, order, Aa, Ba, na, Ab, Bb, nb_, p.decoy_shift, sys_ppm, gate,
                                        rt_window)
        if i.size:
            za = abs(int(round(Aa)))
            frames.append(pd.DataFrame(dict(i=i, j=j, rule=p.name, kind="pair", z=za, n=1, resid=e, sigma=s,
                                            decoy=p.decoy, target=p.target or p.name, nominal=0)))
    cols = ["i", "j", "rule", "kind", "z", "n", "resid", "sigma", "decoy", "target", "nominal"]
    if not frames:
        return pd.DataFrame(columns=cols)
    return pd.concat(frames, ignore_index=True)[cols]
