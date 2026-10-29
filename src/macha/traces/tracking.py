"""Numba core of the multi-target trace tracker.

Every scan is one "frame".  Active traces carry an inverse-variance
weighted m/z estimate with exponential forgetting; centroids of the next scan
are linked to traces by solving a linear assignment problem whose costs are
negative log-likelihoods (m/z residual under an intensity-dependent
precision model, intensity continuity, gap length).  A link is accepted only
if it beats the cost of starting a new trace (``birth_cost``).  Conflicting
candidates are resolved exactly with the Hungarian algorithm on each
connected conflict component, which is almost always tiny.
"""

from __future__ import annotations

import numba as nb
import numpy as np

INF = 1e18


@nb.njit(cache=True)
def hungarian(cost):
    """Minimum-cost assignment of every row (n <= m columns). Returns row -> col."""
    n, m = cost.shape
    u = np.zeros(n + 1)
    v = np.zeros(m + 1)
    p = np.zeros(m + 1, np.int64)
    way = np.zeros(m + 1, np.int64)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = np.full(m + 1, INF)
        used = np.zeros(m + 1, np.bool_)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = INF
            j1 = 0
            for j in range(1, m + 1):
                if not used[j]:
                    cur = cost[i0 - 1, j - 1] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j] = cur
                        way[j] = j0
                    if minv[j] < delta:
                        delta = minv[j]
                        j1 = j
            for j in range(m + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break
    ans = np.full(n, -1, np.int64)
    for j in range(1, m + 1):
        if p[j] != 0:
            ans[p[j] - 1] = j - 1
    return ans


@nb.njit(cache=True)
def _find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


@nb.njit(cache=True)
def solve_links(cand_a, cand_b, gain, n_a, n_b):
    """Maximum-gain matching between a-nodes and b-nodes given candidate pairs.

    Each node may stay unmatched (gain 0).  Returns ``match_a`` (a -> b or -1).
    Connected components of the candidate graph are solved independently: a
    single pair is accepted directly, larger components go to Hungarian.
    """
    match_a = np.full(n_a, -1, np.int64)
    k = cand_a.shape[0]
    if k == 0:
        return match_a
    parent = np.arange(n_a + n_b)
    for e in range(k):
        ra = _find(parent, cand_a[e])
        rb = _find(parent, n_a + cand_b[e])
        if ra != rb:
            parent[ra] = rb
    root = np.empty(k, np.int64)
    for e in range(k):
        root[e] = _find(parent, cand_a[e])
    order = np.argsort(root)
    start = 0
    while start < k:
        stop = start
        r = root[order[start]]
        while stop < k and root[order[stop]] == r:
            stop += 1
        if stop - start == 1:
            e = order[start]
            match_a[cand_a[e]] = cand_b[e]
        else:
            # local index maps
            la = np.empty(stop - start, np.int64)
            lb = np.empty(stop - start, np.int64)
            na = 0
            nb_ = 0
            for q in range(start, stop):
                e = order[q]
                found = False
                for t in range(na):
                    if la[t] == cand_a[e]:
                        found = True
                        break
                if not found:
                    la[na] = cand_a[e]
                    na += 1
                found = False
                for t in range(nb_):
                    if lb[t] == cand_b[e]:
                        found = True
                        break
                if not found:
                    lb[nb_] = cand_b[e]
                    nb_ += 1
            # rows = a nodes, cols = b nodes + one "stay unmatched" column per row
            cost = np.full((na, nb_ + na), 1e9)
            for i in range(na):
                cost[i, nb_ + i] = 0.0
            for q in range(start, stop):
                e = order[q]
                ia = 0
                while la[ia] != cand_a[e]:
                    ia += 1
                ib = 0
                while lb[ib] != cand_b[e]:
                    ib += 1
                if -gain[e] < cost[ia, ib]:
                    cost[ia, ib] = -gain[e]
            ans = hungarian(cost)
            for i in range(na):
                if ans[i] < nb_:
                    match_a[la[i]] = lb[ans[i]]
        start = stop
    return match_a


@nb.njit(cache=True)
def track(mz, intensity, offsets, ppm_a, ppm_b, ppm_c, max_gap, gate, max_ppm, birth_cost,
          gap_cost, int_sigma, int_cap, forget, drift_ppm):
    """Assign every centroid to a trace id. Returns int64 labels (trace ids are dense)."""
    n = mz.shape[0]
    n_scans = offsets.shape[0] - 1
    label = np.full(n, -1, np.int64)
    # per-trace state (capacity n)
    S = np.zeros(n)
    P = np.zeros(n)
    last_scan = np.zeros(n, np.int64)
    last_int = np.zeros(n)
    n_tracks = 0
    active = np.empty(n, np.int64)
    n_active = 0
    tmp = np.empty(n, np.int64)
    means = np.empty(n)
    cand_p = np.empty(1024, np.int64)
    cand_t = np.empty(1024, np.int64)
    cand_g = np.empty(1024)
    for s in range(n_scans):
        a, b = offsets[s], offsets[s + 1]
        # retire traces whose gap would exceed max_gap
        w = 0
        for i in range(n_active):
            t = active[i]
            if s - last_scan[t] - 1 <= max_gap:
                active[w] = t
                w += 1
        n_active = w
        # keep active sorted by current mean (insertion sort, nearly sorted)
        for i in range(n_active):
            means[i] = S[active[i]] / P[active[i]]
        for i in range(1, n_active):
            key_t = active[i]
            key_m = means[i]
            j = i - 1
            while j >= 0 and means[j] > key_m:
                active[j + 1] = active[j]
                means[j + 1] = means[j]
                j -= 1
            active[j + 1] = key_t
            means[j + 1] = key_m
        npts = b - a
        if npts == 0:
            continue
        # candidate generation
        nc = 0
        lo = 0
        for ip in range(npts):
            x = mz[a + ip]
            it = intensity[a + ip]
            sig_ppm = np.sqrt(ppm_a * ppm_a + (ppm_b * it ** (-ppm_c)) ** 2)
            vp = (sig_ppm * 1e-6 * x) ** 2
            win = max_ppm * 1e-6 * x
            while lo < n_active and means[lo] < x - win:
                lo += 1
            j = lo
            while j < n_active and means[j] <= x + win:
                t = active[j]
                gap = s - last_scan[t] - 1
                vt = 1.0 / P[t] + (drift_ppm * 1e-6 * x) ** 2 * (gap + 1)
                d = x - means[j]
                d2 = d * d / (vp + vt)
                if d2 < gate * gate:
                    r = np.log(it / last_int[t])
                    ci = 0.5 * r * r / (int_sigma * int_sigma * (gap + 1))
                    if ci > int_cap:
                        ci = int_cap
                    cost = 0.5 * d2 + 0.5 * np.log((vp + vt) / vp) + ci + gap_cost * gap
                    if cost < birth_cost:
                        if nc >= cand_p.shape[0]:
                            newcap = cand_p.shape[0] * 2
                            cp2 = np.empty(newcap, np.int64)
                            ct2 = np.empty(newcap, np.int64)
                            cg2 = np.empty(newcap)
                            cp2[:nc] = cand_p[:nc]
                            ct2[:nc] = cand_t[:nc]
                            cg2[:nc] = cand_g[:nc]
                            cand_p, cand_t, cand_g = cp2, ct2, cg2
                        cand_p[nc] = ip
                        cand_t[nc] = j
                        cand_g[nc] = birth_cost - cost
                        nc += 1
                j += 1
        match = solve_links(cand_p[:nc], cand_t[:nc], cand_g[:nc], npts, n_active)
        # apply links; unmatched points open new traces
        n_new = 0
        for ip in range(npts):
            x = mz[a + ip]
            it = intensity[a + ip]
            sig_ppm = np.sqrt(ppm_a * ppm_a + (ppm_b * it ** (-ppm_c)) ** 2)
            prec = 1.0 / (sig_ppm * 1e-6 * x) ** 2
            if match[ip] >= 0:
                t = active[match[ip]]
                decay = forget ** (s - last_scan[t])
                S[t] = S[t] * decay + x * prec
                P[t] = P[t] * decay + prec
            else:
                t = n_tracks
                n_tracks += 1
                S[t] = x * prec
                P[t] = prec
                tmp[n_new] = t
                n_new += 1
            last_scan[t] = s
            last_int[t] = it
            label[a + ip] = t
        # merge new traces (already in m/z order) into the active list
        if n_new:
            for i in range(n_active):
                means[i] = S[active[i]] / P[active[i]]
            merged = np.empty(n_active + n_new, np.int64)
            i = 0
            j = 0
            k = 0
            while i < n_active or j < n_new:
                if j >= n_new or (i < n_active and means[i] <= S[tmp[j]] / P[tmp[j]]):
                    merged[k] = active[i]
                    i += 1
                else:
                    merged[k] = tmp[j]
                    j += 1
                k += 1
            active[: n_active + n_new] = merged
            n_active += n_new
    return label


@nb.njit(cache=True)
def trace_stats(label, order, offsets, scan, mz, intensity, ppm_a, ppm_b, ppm_c):
    """Per-trace summary: weighted mean m/z, its SE, scan span, apex, end/start states."""
    nt = offsets.shape[0] - 1
    mean = np.empty(nt)
    se = np.empty(nt)
    sd_ppm = np.empty(nt)
    s0 = np.empty(nt, np.int64)
    s1 = np.empty(nt, np.int64)
    npt = np.empty(nt, np.int64)
    apex_scan = np.empty(nt, np.int64)
    apex_int = np.empty(nt)
    total = np.empty(nt)
    first_int = np.empty(nt)
    last_int = np.empty(nt)
    for t in range(nt):
        a, b = offsets[t], offsets[t + 1]
        sw = 0.0
        swx = 0.0
        best = -1.0
        tot = 0.0
        for q in range(a, b):
            i = order[q]
            sig = np.sqrt(ppm_a * ppm_a + (ppm_b * intensity[i] ** (-ppm_c)) ** 2) * 1e-6 * mz[i]
            w = 1.0 / (sig * sig)
            sw += w
            swx += w * mz[i]
            tot += intensity[i]
            if intensity[i] > best:
                best = intensity[i]
                apex_scan[t] = scan[i]
        m = swx / sw
        mean[t] = m
        se[t] = np.sqrt(1.0 / sw)
        ss = 0.0
        for q in range(a, b):
            i = order[q]
            ss += (mz[i] - m) ** 2
        sd_ppm[t] = np.sqrt(ss / max(b - a - 1, 1)) / m * 1e6
        s0[t] = scan[order[a]]
        s1[t] = scan[order[b - 1]]
        npt[t] = b - a
        apex_int[t] = best
        total[t] = tot
        first_int[t] = intensity[order[a]]
        last_int[t] = intensity[order[b - 1]]
    return mean, se, sd_ppm, s0, s1, npt, apex_scan, apex_int, total, first_int, last_int


@nb.njit(cache=True)
def track_kf(mz, intensity, offsets, ppm_a, ppm_b, ppm_c, max_gap, gate, max_ppm, birth_cost,
             gap_cost, int_sigma, int_cap, q_ppm, v0_ppm):
    """Like :func:`track` but each trace is a constant-velocity Kalman filter.

    State is (m/z, drift per scan) with covariance P.  Prediction over ``ds``
    scans uses F = [[1, ds], [0, 1]] and white-acceleration process noise of
    spectral density ``(q_ppm * 1e-6 * mz)^2`` per scan^3; a new trace starts
    with zero drift and prior sd ``v0_ppm`` ppm/scan.  Gating and costs use
    the innovation variance, so traces follow steady drift without widening
    the m/z window that keeps neighbouring ions apart.
    """
    n = mz.shape[0]
    n_scans = offsets.shape[0] - 1
    label = np.full(n, -1, np.int64)
    m = np.zeros(n)
    v = np.zeros(n)
    p00 = np.zeros(n)
    p01 = np.zeros(n)
    p11 = np.zeros(n)
    last_scan = np.zeros(n, np.int64)
    last_int = np.zeros(n)
    n_tracks = 0
    active = np.empty(n, np.int64)
    n_active = 0
    tmp = np.empty(n, np.int64)
    pred = np.empty(n)
    cand_p = np.empty(1024, np.int64)
    cand_t = np.empty(1024, np.int64)
    cand_g = np.empty(1024)
    for s in range(n_scans):
        a, b = offsets[s], offsets[s + 1]
        w = 0
        for i in range(n_active):
            t = active[i]
            if s - last_scan[t] - 1 <= max_gap:
                active[w] = t
                w += 1
        n_active = w
        for i in range(n_active):
            t = active[i]
            pred[i] = m[t] + v[t] * (s - last_scan[t])
        for i in range(1, n_active):
            key_t = active[i]
            key_m = pred[i]
            j = i - 1
            while j >= 0 and pred[j] > key_m:
                active[j + 1] = active[j]
                pred[j + 1] = pred[j]
                j -= 1
            active[j + 1] = key_t
            pred[j + 1] = key_m
        npts = b - a
        if npts == 0:
            continue
        nc = 0
        lo = 0
        for ip in range(npts):
            x = mz[a + ip]
            it = intensity[a + ip]
            sig_ppm = np.sqrt(ppm_a * ppm_a + (ppm_b * it ** (-ppm_c)) ** 2)
            vp = (sig_ppm * 1e-6 * x) ** 2
            win = max_ppm * 1e-6 * x
            while lo < n_active and pred[lo] < x - win:
                lo += 1
            j = lo
            while j < n_active and pred[j] <= x + win:
                t = active[j]
                ds = s - last_scan[t]
                q = (q_ppm * 1e-6 * x) ** 2
                s00 = p00[t] + 2 * ds * p01[t] + ds * ds * p11[t] + q * ds ** 3 / 3.0
                var = s00 + vp
                d = x - pred[j]
                d2 = d * d / var
                if d2 < gate * gate:
                    gap = ds - 1
                    r = np.log(it / last_int[t])
                    ci = 0.5 * r * r / (int_sigma * int_sigma * ds)
                    if ci > int_cap:
                        ci = int_cap
                    cost = 0.5 * d2 + 0.5 * np.log(var / vp) + ci + gap_cost * gap
                    if cost < birth_cost:
                        if nc >= cand_p.shape[0]:
                            newcap = cand_p.shape[0] * 2
                            cp2 = np.empty(newcap, np.int64)
                            ct2 = np.empty(newcap, np.int64)
                            cg2 = np.empty(newcap)
                            cp2[:nc] = cand_p[:nc]
                            ct2[:nc] = cand_t[:nc]
                            cg2[:nc] = cand_g[:nc]
                            cand_p, cand_t, cand_g = cp2, ct2, cg2
                        cand_p[nc] = ip
                        cand_t[nc] = j
                        cand_g[nc] = birth_cost - cost
                        nc += 1
                j += 1
        match = solve_links(cand_p[:nc], cand_t[:nc], cand_g[:nc], npts, n_active)
        n_new = 0
        for ip in range(npts):
            x = mz[a + ip]
            it = intensity[a + ip]
            sig_ppm = np.sqrt(ppm_a * ppm_a + (ppm_b * it ** (-ppm_c)) ** 2)
            R = (sig_ppm * 1e-6 * x) ** 2
            if match[ip] >= 0:
                t = active[match[ip]]
                ds = s - last_scan[t]
                q = (q_ppm * 1e-6 * x) ** 2
                a00 = p00[t] + 2 * ds * p01[t] + ds * ds * p11[t] + q * ds ** 3 / 3.0
                a01 = p01[t] + ds * p11[t] + q * ds * ds / 2.0
                a11 = p11[t] + q * ds
                mp = m[t] + v[t] * ds
                S = a00 + R
                k0 = a00 / S
                k1 = a01 / S
                y = x - mp
                m[t] = mp + k0 * y
                v[t] = v[t] + k1 * y
                p00[t] = (1 - k0) * a00
                p01[t] = (1 - k0) * a01
                p11[t] = a11 - k1 * a01
            else:
                t = n_tracks
                n_tracks += 1
                m[t] = x
                v[t] = 0.0
                p00[t] = R
                p01[t] = 0.0
                p11[t] = (v0_ppm * 1e-6 * x) ** 2
                tmp[n_new] = t
                n_new += 1
            last_scan[t] = s
            last_int[t] = it
            label[a + ip] = t
        if n_new:
            for i in range(n_active):
                t = active[i]
                pred[i] = m[t]
            merged = np.empty(n_active + n_new, np.int64)
            i = 0
            j = 0
            k = 0
            while i < n_active or j < n_new:
                if j >= n_new or (i < n_active and pred[i] <= m[tmp[j]]):
                    merged[k] = active[i]
                    i += 1
                else:
                    merged[k] = tmp[j]
                    j += 1
                k += 1
            active[: n_active + n_new] = merged
            n_active += n_new
    return label


@nb.njit(cache=True)
def kf_loglik(order, offs, scan, mz, intensity, ppm_a, ppm_b, ppm_c, q_ppm, v0_ppm, cap):
    """Innovation log-likelihood of traces under a constant-velocity model (per-point terms capped)."""
    total = 0.0
    count = 0
    for t in range(offs.shape[0] - 1):
        a, b = offs[t], offs[t + 1]
        if b - a < 6:
            continue
        i0 = order[a]
        x0 = mz[i0]
        sig0 = np.sqrt(ppm_a * ppm_a + (ppm_b * intensity[i0] ** (-ppm_c)) ** 2) * 1e-6 * x0
        m_ = x0
        v_ = 0.0
        p00_ = sig0 * sig0
        p01_ = 0.0
        p11_ = (v0_ppm * 1e-6 * x0) ** 2
        last = scan[i0]
        for qq in range(a + 1, b):
            i = order[qq]
            x = mz[i]
            ds = scan[i] - last
            q = (q_ppm * 1e-6 * x) ** 2
            a00 = p00_ + 2 * ds * p01_ + ds * ds * p11_ + q * ds ** 3 / 3.0
            a01 = p01_ + ds * p11_ + q * ds * ds / 2.0
            a11 = p11_ + q * ds
            mp = m_ + v_ * ds
            R = (np.sqrt(ppm_a * ppm_a + (ppm_b * intensity[i] ** (-ppm_c)) ** 2) * 1e-6 * x) ** 2
            S = a00 + R
            y = x - mp
            term = 0.5 * (y * y / S + np.log(S / R))
            total -= min(term, cap)
            count += 1
            k0 = a00 / S
            k1 = a01 / S
            m_ = mp + k0 * y
            v_ = v_ + k1 * y
            p00_ = (1 - k0) * a00
            p01_ = (1 - k0) * a01
            p11_ = a11 - k1 * a01
            last = scan[i]
    return total, count
