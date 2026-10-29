"""Probabilistic correspondence under the drift-trend mixture.

Features from all samples are pooled and cut into m/z channels (after m/z
calibration).  Inside a channel, groups (one feature per sample at most) are
built and refined by per-sample linear assignment.  The cost of putting
feature ``f`` (sample ``j``) into group ``g`` is the negative log posterior
predictive density of its RT and m/z:

* RT: a mixture over drift classes.  For class ``k`` the group's consensus
  time ``mu_gk`` predicts ``mu_gk + w_kj(mu_gk)``; the class posterior
  ``r_gk`` comes from how well the group's members fit each class.  The
  predictive variance adds the between-class spread, so a group whose class
  is still ambiguous (few members) gates loosely and a well-determined group
  gates tightly.
* m/z: Gaussian with the feature's standard error plus the group's.

A feature is left unmatched (starts its own group) when no assignment beats
``birth_cost``.  After the final pass every feature gets a leave-one-out
membership probability against all compatible alternatives in its channel,
including "belongs to no group".
"""

from __future__ import annotations

import numba as nb
import numpy as np

from ..traces.tracking import solve_links


@nb.njit(cache=True)
def _fit_group(idx, nmem, rt, smp, G, W, logpi, sigma_t, mu_prev, mu_k, r, nu):
    K = W.shape[0]
    ll = np.empty(K)
    for k in range(K):
        m = mu_prev[k]
        for _ in range(2):
            acc = 0.0
            for q in range(nmem):
                f = idx[q]
                acc += rt[f] - np.interp(m, G, W[k, smp[f]])
            m = acc / nmem
        mu_k[k] = m
        s = 0.0
        for q in range(nmem):
            f = idx[q]
            d = rt[f] - m - np.interp(m, G, W[k, smp[f]])
            s += 0.5 * (nu + 1.0) * np.log1p(d * d / (nu * sigma_t * sigma_t))
        ll[k] = logpi[k] - s
    mx = ll.max()
    tot = 0.0
    for k in range(K):
        r[k] = np.exp(ll[k] - mx)
        tot += r[k]
    for k in range(K):
        r[k] /= tot


@nb.njit(cache=True)
def _predict(mu_k, r, j, nmem, G, W, sigma_t):
    K = W.shape[0]
    tk = np.empty(K)
    mean = 0.0
    for k in range(K):
        tk[k] = mu_k[k] + np.interp(mu_k[k], G, W[k, j])
        mean += r[k] * tk[k]
    var = sigma_t * sigma_t * (1.0 + 1.0 / max(nmem, 1))
    for k in range(K):
        var += r[k] * (tk[k] - mean) ** 2
    return mean, var


@nb.njit(cache=True)
def _new_group_mu(t, j, G, W, out):
    for k in range(W.shape[0]):
        m = t
        for _ in range(4):
            m = t - np.interp(m, G, W[k, j])
        out[k] = m


@nb.njit(cache=True)
def _cost(rt_f, rse_f, mz_f, se_f, j, mu_k, r, nmem, gmz, gmzvar, G, W, sigma_t, sys_ppm, gate, nu):
    mean, var = _predict(mu_k, r, j, nmem, G, W, sigma_t)
    var += rse_f * rse_f
    dt2 = (rt_f - mean) ** 2 / var
    sm2 = se_f * se_f + gmzvar + (sys_ppm * 1e-6 * mz_f) ** 2
    dm2 = (mz_f - gmz) ** 2 / sm2
    if dm2 > gate * gate:
        return 1e18
    # Student-t in RT (robust to occasional warp misfit / peak-picking jitter), Gaussian in m/z
    return 0.5 * (nu + 1.0) * np.log1p(dt2 / nu) + 0.5 * np.log(var / (sigma_t * sigma_t)) + 0.5 * dm2


@nb.njit(cache=True)
def _far(rt_f, rse_f, mu_k, r, sigma_t, wbin, g0, gstep, dmax):
    """True if no class prediction of this group can come within ``dmax`` predictive sds of rt_f.

    Predictions lie within mu_k +/- the largest |warp| near mu_k (``wbin``: max
    over classes and runs per RT bin) and the predictive sd is at least
    sqrt(sigma_t^2 + rse_f^2), so this is a conservative bound: skipped pairs
    would cost more than ``birth_cost`` anyway.
    """
    lo = mu_k[0]
    hi = mu_k[0]
    for k in range(1, mu_k.shape[0]):
        lo = min(lo, mu_k[k])
        hi = max(hi, mu_k[k])
    nb_ = wbin.shape[0]
    b0 = min(max(int((lo - g0) / gstep) - 1, 0), nb_ - 1)
    b1 = min(max(int((hi - g0) / gstep) + 1, 0), nb_ - 1)
    wmax = 0.0
    for b in range(b0, b1 + 1):
        wmax = max(wmax, wbin[b])
    d = max(lo - wmax - rt_f, rt_f - hi - wmax, 0.0)
    return d > dmax * np.sqrt(sigma_t * sigma_t + rse_f * rse_f)


@nb.njit(cache=True)
def _refit(g, gmem, rt, mz, se, smp, G, W, logpi, sigma_t, gmu, gr, gmz, gmzvar, gcnt, idx_buf, nu):
    nmem = 0
    for j in range(gmem.shape[1]):  # a group holds at most one feature per run
        if gmem[g, j] >= 0:
            idx_buf[nmem] = gmem[g, j]
            nmem += 1
    gcnt[g] = nmem
    if nmem == 0:
        return
    _fit_group(idx_buf, nmem, rt, smp, G, W, logpi, sigma_t, gmu[g].copy(), gmu[g], gr[g], nu)
    sw = 0.0
    swx = 0.0
    for q in range(nmem):
        f = idx_buf[q]
        w = 1.0 / (se[f] * se[f])
        sw += w
        swx += w * mz[f]
    gmz[g] = swx / sw
    gmzvar[g] = 1.0 / sw


@nb.njit(cache=True)
def _solve_channel(rt, rse, mz, se, smp, J, G, W, logpi, sigma_t, sys_ppm, birth_cost, gate, n_refine, nu):
    # local warp bound per RT bin (max |w| over classes and runs) for the pruning test in _far
    n_bins = 64
    gstep = (G[-1] - G[0]) / n_bins + 1e-9
    wbin = np.zeros(n_bins)
    for q in range(G.shape[0]):
        b = min(int((G[q] - G[0]) / gstep), n_bins - 1)
        for k in range(W.shape[0]):
            for jj in range(W.shape[1]):
                wbin[b] = max(wbin[b], abs(W[k, jj, q]))
    g0 = G[0]
    # Student-t cost exceeds birth_cost beyond d^2 > nu (exp(2 birth / (nu + 1)) - 1); add margin
    dmax = 1.5 * np.sqrt(nu * (np.exp(2.0 * birth_cost / (nu + 1.0)) - 1.0)) + 1.0
    n = rt.shape[0]
    K = W.shape[0]
    pi = np.exp(logpi)
    asg = np.full(n, -1, np.int64)
    cap = n * (n_refine + 2) + 1  # groups born in every pass; dead slots are not reused
    gmu = np.zeros((cap, K))
    gr = np.zeros((cap, K))
    gmz = np.zeros(cap)
    gmzvar = np.zeros(cap)
    gcnt = np.zeros(cap, np.int64)
    gmem = np.full((cap, J), -1, np.int64)  # member feature of each group per run
    ng = 0
    idx_buf = np.empty(max(n, J), np.int64)
    counts = np.zeros(J, np.int64)
    for f in range(n):
        counts[smp[f]] += 1
    order = np.argsort(-counts, kind="mergesort")
    feats = np.empty(n, np.int64)
    for p in range(n_refine + 1):
        for oj in range(J):
            j = order[oj] if p == 0 else oj
            if counts[j] == 0:
                continue
            nf = 0
            for f in range(n):
                if smp[f] == j:
                    feats[nf] = f
                    nf += 1
            if p > 0:  # leave this sample out, refit touched groups
                for q in range(nf):
                    g = asg[feats[q]]
                    asg[feats[q]] = -1
                    if g >= 0:
                        gmem[g, j] = -1
                        _refit(g, gmem, rt, mz, se, smp, G, W, logpi, sigma_t, gmu, gr, gmz, gmzvar, gcnt, idx_buf, nu)
            # candidate links to live groups
            ca = []
            cb = []
            cg = []
            for q in range(nf):
                f = feats[q]
                for g in range(ng):
                    if gcnt[g] == 0 or _far(rt[f], rse[f], gmu[g], gr[g], sigma_t, wbin, g0, gstep, dmax):
                        continue
                    c = _cost(rt[f], rse[f], mz[f], se[f], j, gmu[g], gr[g], gcnt[g], gmz[g], gmzvar[g], G, W, sigma_t,
                              sys_ppm, gate, nu)
                    if c < birth_cost:
                        ca.append(q)
                        cb.append(g)
                        cg.append(birth_cost - c)
            match = solve_links(np.array(ca, np.int64), np.array(cb, np.int64), np.array(cg, np.float64), nf, max(ng, 1))
            touched = np.zeros(cap, np.bool_)
            for q in range(nf):
                f = feats[q]
                if match[q] >= 0:
                    asg[f] = match[q]
                    gmem[match[q], j] = f
                    touched[match[q]] = True
                else:
                    g = ng
                    ng += 1
                    asg[f] = g
                    gmem[g, j] = f
                    _new_group_mu(rt[f], j, G, W, gmu[g])
                    for k in range(K):
                        gr[g, k] = pi[k]
                    gmz[g] = mz[f]
                    gmzvar[g] = se[f] * se[f]
                    gcnt[g] = 1
            for g in range(ng):
                if touched[g]:
                    _refit(g, gmem, rt, mz, se, smp, G, W, logpi, sigma_t, gmu, gr, gmz, gmzvar, gcnt, idx_buf, nu)
    # leave-one-out membership probabilities
    memb = np.zeros(n)
    mu_tmp = np.zeros(K)
    r_tmp = np.zeros(K)
    for f in range(n):
        j = smp[f]
        g0 = asg[f]
        # LOO fit of own group
        nmem = 0
        for jj in range(J):
            h = gmem[g0, jj]
            if h >= 0 and h != f:
                idx_buf[nmem] = h
                nmem += 1
        num = 0.0
        if nmem > 0:
            for k in range(K):
                mu_tmp[k] = gmu[g0, k]
            _fit_group(idx_buf, nmem, rt, smp, G, W, logpi, sigma_t, gmu[g0].copy(), mu_tmp, r_tmp, nu)
            sw = 0.0
            swx = 0.0
            for q in range(nmem):
                h = idx_buf[q]
                w = 1.0 / (se[h] * se[h])
                sw += w
                swx += w * mz[h]
            c0 = _cost(rt[f], rse[f], mz[f], se[f], j, mu_tmp, r_tmp, nmem, swx / sw, 1.0 / sw, G, W, sigma_t, sys_ppm, gate, nu)
            num = np.exp(-c0) if c0 < 1e17 else 0.0
        den = num + np.exp(-birth_cost)
        for g in range(ng):
            if g == g0 or gcnt[g] == 0:
                continue
            if gmem[g, j] >= 0:  # groups that already hold a feature from this run are not alternatives
                continue
            if _far(rt[f], rse[f], gmu[g], gr[g], sigma_t, wbin, g0, gstep, dmax):
                continue
            c = _cost(rt[f], rse[f], mz[f], se[f], j, gmu[g], gr[g], gcnt[g], gmz[g], gmzvar[g], G, W, sigma_t, sys_ppm, gate, nu)
            if c < 1e17:
                den += np.exp(-c)
        memb[f] = (num if nmem > 0 else np.exp(-birth_cost)) / den
    return asg, memb, ng, gmu[:ng], gr[:ng], gmz[:ng], gcnt[:ng]


@nb.njit(cache=True)
def solve_all(order, offs, rt, rse, mz, se, smp, J, G, W, logpi, sigma_t, sys_ppm, birth_cost, gate, n_refine,
              nu):
    n_all = rt.shape[0]
    K = W.shape[0]
    group = np.full(n_all, -1, np.int64)
    memb = np.zeros(n_all)
    cap = n_all
    g_mu = np.zeros((cap, K))
    g_r = np.zeros((cap, K))
    g_mz = np.zeros(cap)
    g_n = np.zeros(cap, np.int64)
    base = 0
    for c in range(offs.shape[0] - 1):
        idx = order[offs[c]: offs[c + 1]]
        asg, mb, ng, gmu, gr, gmz, gcnt = _solve_channel(rt[idx], rse[idx], mz[idx], se[idx], smp[idx], J, G, W, logpi,
                                                         sigma_t, sys_ppm, birth_cost, gate, n_refine, nu)
        # compact (drop emptied groups)
        remap = np.full(ng, -1, np.int64)
        k = 0
        for g in range(ng):
            if gcnt[g] > 0:
                remap[g] = base + k
                g_mu[base + k] = gmu[g]
                g_r[base + k] = gr[g]
                g_mz[base + k] = gmz[g]
                g_n[base + k] = gcnt[g]
                k += 1
        for q in range(idx.shape[0]):
            group[idx[q]] = remap[asg[q]]
            memb[idx[q]] = mb[q]
        base += k
    return group, memb, g_mu[:base], g_r[:base], g_mz[:base], g_n[:base]


def mz_channels(mz, se, ppm: float = 5.0, gate: float = 4.0, max_size: int = 400):
    """Cut sorted pooled m/z values into channels at gaps beyond tolerance; split oversized channels."""
    order = np.argsort(mz, kind="stable")
    m = mz[order]
    s = se[order]
    gap = np.diff(m)
    tol = np.maximum(ppm * 1e-6 * m[1:], gate * np.sqrt(s[1:] ** 2 + s[:-1] ** 2))
    new = np.concatenate([[True], gap > tol])
    starts = list(np.flatnonzero(new)) + [m.size]
    bounds = []
    for a, b in zip(starts[:-1], starts[1:]):
        stack = [(a, b)]
        while stack:
            lo, hi = stack.pop()
            if hi - lo <= max_size:
                bounds.append((lo, hi))
                continue
            cut = lo + 1 + int(np.argmax(gap[lo: hi - 1]))
            stack += [(lo, cut), (cut, hi)]
    bounds.sort()
    offs = np.array([b[0] for b in bounds] + [m.size], np.int64)
    return order.astype(np.int64), offs
