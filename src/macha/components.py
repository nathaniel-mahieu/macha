"""Components: peak models of scan-graph relation clusters, per run and across runs.

A *component* is one compound's elution in one run: the ions its confident relations join, sharing one elution
profile.  Relations say which ions co-elute; a component says how: one exponentially modified Gaussian (EMG) shape
for all of its ions, with an amplitude per ion,

    y_i(t) = a_i * f(t; mu, sigma, tau) + b_i                         (one component)
    y_i(t) = a_i1 * f(t; mu_1, sigma, tau) + a_i2 * f(t; mu_2, sigma, tau) + b_i   (two)

fitted by variable projection (amplitudes and the ion's constant chemical-background level ``b_i`` by non-negative
least squares in closed form, the shape by Levenberg-Marquardt).  Every partner trace
is a separate noisy view of the same profile, so the shape is pinned by all ions together.  The intensity variance
is ``noise_i^2 + (cv * model)^2`` with ``noise_i`` measured in the ion's own m/z window
(:mod:`macha.traces.noise`); a scan inside the window where the ion has no centroid is censored at that level (it
costs only where the model rises above it), which is the dropout floor.

**Merged traces.**  Ions closer than about one resolution FWHM are not resolved: the instrument reports a single
centroid at their intensity-weighted mean m/z.  With two components in one ion, the centroid m/z follows

    m(t) = m_2 + (m_1 - m_2) * g(t) + beta * log I(t),    g = a_1 f_1 / (a_1 f_1 + a_2 f_2)

so it moves from one source m/z to the other as the dominant source changes, while an intensity-dependent shift of
one ion (the ``beta`` term) gives the same m/z at the same intensity on both flanks.  Centroid precision
(``traces.precision``, ~0.3 ppm at high intensity) is far below the resolution limit (several ppm), so two sources
a few ppm apart that never separate in time still leave a measurable trajectory.  The m/z residuals enter the same
objective, weighted by the precision model, so the split is decided by intensity shape and m/z together.

The split is kept when the two-component fit beats one component by more than a BIC-style penalty
(``split_penalty`` per extra parameter, times log of the number of observations).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numba as nb
import numpy as np
import pandas as pd

from .traces.noise import floor_halfwidth_ppm, window_noise

SQRT2 = math.sqrt(2.0)
SQRT_PI = math.sqrt(math.pi)
MZ_CLIP = 4.0  # m/z residuals beyond this many sd count as this many (heavy-tailed low-intensity centroids)


# ------------------------------------------------------------------------------------------------ EMG
@nb.njit(cache=True)
def _erfcx(u):
    if u < 25.0:
        return math.exp(u * u) * math.erfc(u)
    u2 = u * u
    return (1.0 - 0.5 / u2 + 0.75 / (u2 * u2)) / (u * SQRT_PI)


@nb.njit(cache=True)
def _emg(t, mu, sig, tau):
    """Area-normalised EMG density (Gaussian limit as tau -> 0), stable for any argument."""
    z = (t - mu) / sig
    u = (sig / tau - z) / SQRT2
    if u >= 0.0:
        return 0.5 / tau * math.exp(-0.5 * z * z) * _erfcx(u)
    return 0.5 / tau * math.exp(u * u - 0.5 * z * z) * math.erfc(u)


@nb.njit(cache=True)
def _emg_mode(mu, sig, tau):
    lo = mu - sig
    hi = mu + 0.5 * sig + min(tau, 4.0 * sig)
    n = 48
    best = lo
    fb = -1.0
    for k in range(n + 1):
        t = lo + (hi - lo) * k / n
        f = _emg(t, mu, sig, tau)
        if f > fb:
            fb = f
            best = t
    h = (hi - lo) / n
    a, b = best - h, best + h
    for _ in range(40):
        c1 = b - 0.618 * (b - a)
        c2 = a + 0.618 * (b - a)
        if _emg(c1, mu, sig, tau) > _emg(c2, mu, sig, tau):
            b = c2
        else:
            a = c1
    return 0.5 * (a + b)


@nb.njit(cache=True)
def _emg_bounds(mu, sig, tau, mode, frac):
    """Times on either side of the mode where the profile falls to ``frac`` of its apex."""
    top = _emg(mode, mu, sig, tau) * frac
    out = np.empty(2)
    for side in range(2):
        d = -1.0 if side == 0 else 1.0
        a = mode
        step = 0.25 * sig
        b = mode + d * step
        k = 0
        while _emg(b, mu, sig, tau) > top and k < 400:
            a = b
            b = b + d * step
            k += 1
        for _ in range(40):
            m = 0.5 * (a + b)
            if _emg(m, mu, sig, tau) > top:
                a = m
            else:
                b = m
        out[side] = 0.5 * (a + b)
    return out


def emg(t, mu, sig, tau):
    """Vectorised area-normalised EMG density."""
    t = np.asarray(t, float)
    out = np.empty(t.shape)
    flat = out.reshape(-1)
    for k, x in enumerate(t.reshape(-1)):
        flat[k] = _emg(float(x), float(mu), float(sig), float(tau))
    return out


# ------------------------------------------------------------------------------------------------ data layout
@nb.njit(cache=True)
def _gather(order, offs, scan_of, mz, intensity, ion_trace, ion_off, s_lo, L, dat_off, n_dat):
    """Dense ion x scan arrays (intensity, m/z, observed) of every compound's window from its traces' points."""
    Y = np.zeros(n_dat)
    MZ = np.full(n_dat, np.nan)
    OBS = np.zeros(n_dat, np.bool_)
    for q in range(s_lo.shape[0]):
        for ii in range(ion_off[q], ion_off[q + 1]):
            t = ion_trace[ii]
            base = dat_off[q] + (ii - ion_off[q]) * L[q]
            for k in range(offs[t], offs[t + 1]):
                p = order[k]
                s = scan_of[p] - s_lo[q]
                if 0 <= s < L[q]:
                    if intensity[p] > Y[base + s]:
                        Y[base + s] = intensity[p]
                        MZ[base + s] = mz[p]
                    OBS[base + s] = True
    return Y, MZ, OBS


# ------------------------------------------------------------------------------------------------ fitting kernel
@nb.njit(cache=True)
def _solve(A, b):
    """Small dense solve by Gaussian elimination with partial pivoting (returns zeros if singular)."""
    n = b.shape[0]
    M = A.copy()
    x = b.copy()
    for c in range(n):
        p = c
        for r in range(c + 1, n):
            if abs(M[r, c]) > abs(M[p, c]):
                p = r
        if abs(M[p, c]) < 1e-300:
            return np.zeros(n)
        if p != c:
            for k in range(n):
                tmp = M[c, k]
                M[c, k] = M[p, k]
                M[p, k] = tmp
            tmp = x[c]
            x[c] = x[p]
            x[p] = tmp
        for r in range(c + 1, n):
            f = M[r, c] / M[c, c]
            for k in range(c, n):
                M[r, k] -= f * M[c, k]
            x[r] -= f * x[c]
    for c in range(n - 1, -1, -1):
        s = x[c]
        for k in range(c + 1, n):
            s -= M[c, k] * x[k]
        x[c] = s / M[c, c]
    return x


@nb.njit(cache=True)
def _nnls_small(G, r):
    """min x'Gx - 2 r'x subject to x >= 0 for a handful of unknowns, by enumerating active sets."""
    n = r.shape[0]
    best = 0.0
    out = np.zeros(n)
    for mask in range(1, 1 << n):
        idx = np.empty(n, np.int64)
        m = 0
        for a in range(n):
            if mask & (1 << a):
                idx[m] = a
                m += 1
        Gs = np.empty((m, m))
        rs = np.empty(m)
        for a in range(m):
            rs[a] = r[idx[a]]
            for b in range(m):
                Gs[a, b] = G[idx[a], idx[b]]
            Gs[a, a] += 1e-12 * (Gs[a, a] + 1e-300)
        x = _solve(Gs, rs)
        ok = True
        obj = 0.0
        for a in range(m):
            if not (x[a] >= 0.0):
                ok = False
                break
            obj -= rs[a] * x[a]
        if ok and obj < best:
            best = obj
            out[:] = 0.0
            for a in range(m):
                out[idx[a]] = x[a]
    return out


@nb.njit(cache=True)
def _resid(theta, K, rt, L, n_ion, Y, MZ, OBS, W, U, LI, floor, use_mz, dm_sd, out, amp, mzp):
    """Residual vector (intensity part, then m/z part) of a K-component shared-shape model.

    ``theta`` = (mu_1[, mu_2], log sigma, log tau).  Amplitudes are profiled out per ion by weighted non-negative
    least squares; m/z parameters (source m/z, and with two components the m/z difference) by weighted least
    squares.  Returns (intensity chi2, m/z chi2).
    """
    sig = math.exp(theta[K])
    tau = math.exp(theta[K + 1])
    f = np.empty((K, L))
    for k in range(K):
        for s in range(L):
            f[k, s] = _emg(rt[s], theta[k], sig, tau)
    ci = 0.0
    cm = 0.0
    N = n_ion * L
    for i in range(n_ion):
        b0 = i * L
        # amplitudes and a constant baseline (the ion's chemical-background level), non-negative
        nc = K + 1
        G = np.zeros((nc, nc))
        rr = np.zeros(nc)
        col = np.empty(nc)
        for s in range(L):
            if OBS[b0 + s]:
                w = W[b0 + s]
                for k in range(K):
                    col[k] = f[k, s]
                col[K] = 1.0
                for a in range(nc):
                    rr[a] += w * Y[b0 + s] * col[a]
                    for b in range(nc):
                        G[a, b] += w * col[a] * col[b]
        x = _nnls_small(G, rr)
        amp[i, 0] = x[0]
        amp[i, 1] = x[1] if K == 2 else 0.0
        amp[i, 2] = x[K]
        # intensity residuals (censored where the ion has no centroid)
        for s in range(L):
            m = amp[i, 0] * f[0, s] + amp[i, 2]
            if K == 2:
                m += amp[i, 1] * f[1, s]
            idx = b0 + s
            if OBS[idx]:
                r = (Y[idx] - m) * math.sqrt(W[idx])
            else:
                r = max(m - floor[i], 0.0) * math.sqrt(W[idx])
            out[idx] = r
            ci += r * r
        # m/z residuals
        for s in range(L):
            out[N + b0 + s] = 0.0
        mzp[i, 0] = np.nan
        mzp[i, 1] = np.nan
        mzp[i, 2] = 0.0
        if not use_mz:
            continue
        two = K == 2 and amp[i, 0] > 0.0 and amp[i, 1] > 0.0
        nobs = 0
        sw = 0.0
        sg = 0.0
        sgg = 0.0
        for s in range(L):
            idx = b0 + s
            if OBS[idx] and MZ[idx] == MZ[idx]:
                nobs += 1
                if two:
                    g = amp[i, 0] * f[0, s] / max(amp[i, 0] * f[0, s] + amp[i, 1] * f[1, s], 1e-300)
                    sw += U[idx]
                    sg += U[idx] * g
                    sgg += U[idx] * g * g
        if nobs < 3:
            continue
        if two and sw > 0 and (sgg / sw - (sg / sw) ** 2) < 1e-4:
            two = False  # profiles too alike in this ion: the m/z difference is not identified
        use_b = nobs >= 8
        p = 1 + (1 if two else 0) + (1 if use_b else 0)
        A = np.zeros((p, p))
        bb = np.zeros(p)
        x = np.zeros(p)
        ref = 0.0
        for s in range(L):
            idx = b0 + s
            if OBS[idx] and MZ[idx] == MZ[idx]:
                ref = MZ[idx]
                break
        for s in range(L):
            idx = b0 + s
            if not (OBS[idx] and MZ[idx] == MZ[idx]):
                continue
            x[0] = 1.0
            c = 1
            if two:
                x[c] = amp[i, 0] * f[0, s] / max(amp[i, 0] * f[0, s] + amp[i, 1] * f[1, s], 1e-300)
                c += 1
            if use_b:
                x[c] = LI[idx]
            u = U[idx]
            for a in range(p):
                bb[a] += u * x[a] * (MZ[idx] - ref)
                for b in range(p):
                    A[a, b] += u * x[a] * x[b]
        for a in range(p):
            A[a, a] += 1e-12 * (A[a, a] + 1.0)
        # two sources only coalesce when closer than the resolution limit: Gaussian prior on their m/z difference
        prior = 0.0
        if two:
            prior = 1.0 / (dm_sd * 1e-6 * ref) ** 2
            A[1, 1] += prior
        beta = _solve(A, bb)
        if two:
            cm += beta[1] * beta[1] * prior
        for s in range(L):
            idx = b0 + s
            if not (OBS[idx] and MZ[idx] == MZ[idx]):
                continue
            pred = beta[0]
            c = 1
            if two:
                pred += beta[c] * amp[i, 0] * f[0, s] / max(amp[i, 0] * f[0, s] + amp[i, 1] * f[1, s], 1e-300)
                c += 1
            if use_b:
                pred += beta[c] * LI[idx]
            z = (MZ[idx] - ref - pred) * math.sqrt(U[idx])
            z = min(max(z, -MZ_CLIP), MZ_CLIP)
            out[N + idx] = z
            cm += z * z
        # source m/z at the ion's mean log intensity (LI is centred): component 2 at g = 0, component 1 at g = 1
        mzp[i, 1] = ref + beta[0]
        mzp[i, 0] = ref + beta[0] + (beta[1] if two else 0.0)
        mzp[i, 2] = 1.0 if two else 0.0
    return ci, cm


@nb.njit(cache=True)
def _lm(theta0, K, lo, hi, rt, L, n_ion, Y, MZ, OBS, W, U, LI, floor, use_mz, dm_sd, amp, mzp, n_iter):
    p = theta0.shape[0]
    N = 2 * n_ion * L
    theta = theta0.copy()
    for j in range(p):
        theta[j] = min(max(theta[j], lo[j]), hi[j])
    r = np.empty(N)
    r2 = np.empty(N)
    amp2 = np.empty_like(amp)
    mzp2 = np.empty_like(mzp)
    ci, cm = _resid(theta, K, rt, L, n_ion, Y, MZ, OBS, W, U, LI, floor, use_mz, dm_sd, r, amp, mzp)
    chi = ci + cm
    lam = 1e-2
    J = np.empty((N, p))
    steps = np.empty(p)
    sig = math.exp(theta[K])
    for j in range(p):
        steps[j] = 0.02 * sig if j < K else 0.01
    for it in range(n_iter):
        for j in range(p):
            th = theta.copy()
            h = steps[j]
            if th[j] + h > hi[j]:
                h = -h
            th[j] += h
            _resid(th, K, rt, L, n_ion, Y, MZ, OBS, W, U, LI, floor, use_mz, dm_sd, r2, amp2, mzp2)
            for n in range(N):
                J[n, j] = (r2[n] - r[n]) / h
        A = np.zeros((p, p))
        g = np.zeros(p)
        for n in range(N):
            if r[n] == 0.0:
                zero_row = True
                for j in range(p):
                    if J[n, j] != 0.0:
                        zero_row = False
                        break
                if zero_row:
                    continue
            for a in range(p):
                g[a] += J[n, a] * r[n]
                for b in range(a, p):
                    A[a, b] += J[n, a] * J[n, b]
        for a in range(p):
            for b in range(a):
                A[a, b] = A[b, a]
        improved = False
        rel = 0.0
        for _ in range(8):
            M = A.copy()
            for a in range(p):
                M[a, a] += lam * (A[a, a] + 1e-9)
            delta = _solve(M, -g)
            th = theta + delta
            for j in range(p):
                th[j] = min(max(th[j], lo[j]), hi[j])
            a_, b_ = _resid(th, K, rt, L, n_ion, Y, MZ, OBS, W, U, LI, floor, use_mz, dm_sd, r2, amp2, mzp2)
            c = a_ + b_
            if c < chi:
                rel = (chi - c) / max(chi, 1e-300)
                theta = th
                r[:] = r2
                chi = c
                lam = max(lam * 0.3, 1e-7)
                improved = True
                break
            lam *= 10.0
        if not improved or rel < 1e-7:
            break
    ci, cm = _resid(theta, K, rt, L, n_ion, Y, MZ, OBS, W, U, LI, floor, use_mz, dm_sd, r, amp, mzp)
    return theta, ci, cm


@nb.njit(cache=True)
def _weights(theta, K, amp, rt, L, n_ion, Y, OBS, noise, cv, W):
    """Intensity weights 1 / (noise^2 + (cv * model)^2) from a fitted model (data where the model is empty)."""
    sig = math.exp(theta[K])
    tau = math.exp(theta[K + 1])
    for i in range(n_ion):
        for s in range(L):
            m = amp[i, 2]
            for k in range(K):
                m += amp[i, k] * _emg(rt[s], theta[k], sig, tau)
            idx = i * L + s
            v = noise[i] * noise[i] + cv * cv * m * m
            W[idx] = 1.0 / v


@nb.njit(cache=True, parallel=True)
def _fit_all(rt_all, s_lo, L, ion_off, dat_off, Y, MZ, OBS, U, LI, noise, floor, mu0, sig0, cv, use_mz, dm_sd,
             try2, n_iter):
    Q = s_lo.shape[0]
    n_ions = noise.shape[0]
    P1 = np.full((Q, 3), np.nan)
    P2 = np.full((Q, 4), np.nan)
    C1 = np.zeros((Q, 2))
    C2 = np.full((Q, 2), np.nan)
    NOBS = np.zeros(Q, np.int64)
    A1 = np.zeros(n_ions)
    A2 = np.zeros((n_ions, 2))
    B = np.zeros((n_ions, 2))  # baseline per ion, one- and two-component fit
    M1 = np.full(n_ions, np.nan)
    M2 = np.full((n_ions, 3), np.nan)
    HI = np.zeros((Q, 2))  # sum of squared relative residuals and count at high intensity (for cv)
    for q in nb.prange(Q):
        n_ion = ion_off[q + 1] - ion_off[q]
        Lq = L[q]
        d0 = dat_off[q]
        d1 = d0 + n_ion * Lq
        rt = rt_all[s_lo[q]:s_lo[q] + Lq]
        Yq = Y[d0:d1]
        MZq = MZ[d0:d1]
        OBSq = OBS[d0:d1]
        Uq = U[d0:d1]
        LIq = LI[d0:d1]
        nq = noise[ion_off[q]:ion_off[q + 1]]
        fq = floor[ion_off[q]:ion_off[q + 1]]
        nobs = 0
        for k in range(OBSq.shape[0]):
            if OBSq[k]:
                nobs += 1
        NOBS[q] = nobs
        if nobs < 4:
            continue
        s0 = sig0[q]
        lo1 = np.array([rt[0], math.log(0.3 * s0), math.log(0.02 * s0)])
        hi1 = np.array([rt[Lq - 1], math.log(3.0 * s0), math.log(4.0 * s0)])
        th = np.array([mu0[q], math.log(s0), math.log(0.5 * s0)])
        W = np.empty(n_ion * Lq)
        for k in range(W.shape[0]):
            y = Yq[k] if OBSq[k] else 0.0
            W[k] = 1.0 / (nq[k // Lq] ** 2 + cv * cv * y * y)
        amp = np.zeros((n_ion, 3))
        mzp = np.zeros((n_ion, 3))
        # one component: fit, re-weight from the model, refit
        for _ in range(2):
            th, ci, cm = _lm(th, 1, lo1, hi1, rt, Lq, n_ion, Yq, MZq, OBSq, W, Uq, LIq, fq, use_mz, dm_sd, amp, mzp, n_iter)
            _weights(th, 1, amp, rt, Lq, n_ion, Yq, OBSq, nq, cv, W)
        th, ci, cm = _lm(th, 1, lo1, hi1, rt, Lq, n_ion, Yq, MZq, OBSq, W, Uq, LIq, fq, use_mz, dm_sd, amp, mzp, n_iter)
        P1[q, 0] = th[0]
        P1[q, 1] = math.exp(th[1])
        P1[q, 2] = math.exp(th[2])
        C1[q, 0] = ci
        C1[q, 1] = cm
        sig = math.exp(th[1])
        tau = math.exp(th[2])
        for i in range(n_ion):
            A1[ion_off[q] + i] = amp[i, 0]
            B[ion_off[q] + i, 0] = amp[i, 2]
            M1[ion_off[q] + i] = mzp[i, 1]
            for s in range(Lq):
                idx = i * Lq + s
                m = amp[i, 0] * _emg(rt[s], th[0], sig, tau) + amp[i, 2]
                if OBSq[idx] and m > 10.0 * nq[i]:
                    rr = (Yq[idx] - m) / m
                    HI[q, 0] += rr * rr
                    HI[q, 1] += 1.0
        if not try2 or nobs < 8:
            continue
        # two components, same weights (so the chi2 values compare), several starts around the one-component fit
        mu = th[0]
        lo2 = np.array([rt[0], rt[0], lo1[1], lo1[2]])
        hi2 = np.array([rt[Lq - 1], rt[Lq - 1], hi1[1], hi1[2]])
        best = 1e300
        amp2 = np.zeros((n_ion, 3))
        mzp2 = np.zeros((n_ion, 3))
        for start in range(3):
            if start == 0:
                a, b = mu - 0.8 * sig, mu + 0.8 * sig
            elif start == 1:
                a, b = mu - 1.6 * sig, mu + 0.4 * sig
            else:
                a, b = mu - 0.4 * sig, mu + 1.6 * sig
            t2 = np.array([a, b, th[1], th[2]])
            t2, ci2, cm2 = _lm(t2, 2, lo2, hi2, rt, Lq, n_ion, Yq, MZq, OBSq, W, Uq, LIq, fq, use_mz, dm_sd, amp2,
                               mzp2, n_iter)
            if ci2 + cm2 < best:
                best = ci2 + cm2
                P2[q, 0] = t2[0]
                P2[q, 1] = t2[1]
                P2[q, 2] = math.exp(t2[2])
                P2[q, 3] = math.exp(t2[3])
                C2[q, 0] = ci2
                C2[q, 1] = cm2
                for i in range(n_ion):
                    A2[ion_off[q] + i, 0] = amp2[i, 0]
                    A2[ion_off[q] + i, 1] = amp2[i, 1]
                    B[ion_off[q] + i, 1] = amp2[i, 2]
                    M2[ion_off[q] + i, 0] = mzp2[i, 0]
                    M2[ion_off[q] + i, 1] = mzp2[i, 1]
                    M2[ion_off[q] + i, 2] = mzp2[i, 2]
    return P1, P2, C1, C2, NOBS, A1, A2, M1, M2, HI, B


@nb.njit(cache=True, parallel=True)
def _shape_stats(P1, P2, split):
    """Mode, FWHM and 5 % bounds of each fitted component (two rows per compound; NaN for an absent second)."""
    Q = P1.shape[0]
    out = np.full((Q, 2, 4), np.nan)
    for q in nb.prange(Q):
        for k in range(2):
            if split[q]:
                mu, sig, tau = P2[q, k], P2[q, 2], P2[q, 3]
            elif k == 0:
                mu, sig, tau = P1[q, 0], P1[q, 1], P1[q, 2]
            else:
                continue
            if not (mu == mu):
                continue
            mode = _emg_mode(mu, sig, tau)
            half = _emg_bounds(mu, sig, tau, mode, 0.5)
            b5 = _emg_bounds(mu, sig, tau, mode, 0.05)
            out[q, k, 0] = mode
            out[q, k, 1] = half[1] - half[0]
            out[q, k, 2] = b5[0]
            out[q, k, 3] = b5[1]
    return out


@nb.njit(cache=True)
def _peak_density(P1, P2, split):
    """EMG density at the mode per component (height = amplitude x this)."""
    Q = P1.shape[0]
    out = np.full((Q, 2), np.nan)
    for q in range(Q):
        for k in range(2):
            if split[q]:
                mu, sig, tau = P2[q, k], P2[q, 2], P2[q, 3]
            elif k == 0:
                mu, sig, tau = P1[q, 0], P1[q, 1], P1[q, 2]
            else:
                continue
            if mu == mu:
                out[q, k] = _emg(_emg_mode(mu, sig, tau), mu, sig, tau)
    return out


# ------------------------------------------------------------------------------------------------ per run
@dataclass
class Components:
    """Components of one run.

    ``table``: one row per component (compound, k = 0 / 1 after a split) with apex, shape, area, height, SNR and
    fit statistics.  ``ions``: one row per (component, ion) with amplitude, height, share of the ion's signal and
    the ion's source m/z (for a split ion, the m/z of this component's source).  ``fits``: one row per compound
    with both fits' chi2 and the split evidence.  ``relation_component``: for every row of the relation table, the
    component it belongs to (-1 if none).
    """

    table: pd.DataFrame
    ions: pd.DataFrame
    fits: pd.DataFrame
    relation_component: np.ndarray
    cv: float


def _layout(rel, traces, sigma, pad_sd, max_ions, max_scans):
    run = traces.run
    r = rel[rel["compound"].to_numpy() >= 0]
    comp_ids = np.unique(r["compound"].to_numpy())
    cyc = float(run.cycle_time)
    pad = int(math.ceil(pad_sd * sigma / max(cyc, 1e-6)))
    s0 = r.groupby("compound")["s0"].min().reindex(comp_ids).to_numpy() - pad
    s1 = r.groupby("compound")["s1"].max().reindex(comp_ids).to_numpy() + pad
    s0 = np.maximum(s0, 0)
    s1 = np.minimum(s1, run.n_scans - 1)
    # very long windows (background) are cut around the compound's most intense relation
    long = (s1 - s0 + 1) > max_scans
    if long.any():
        top = r.sort_values("height", ascending=False).drop_duplicates("compound").set_index("compound")
        c = run.scan_of_rt(top.loc[comp_ids[long], "rt"].to_numpy()) if hasattr(run, "scan_of_rt") else \
            np.searchsorted(run.rt, top.loc[comp_ids[long], "rt"].to_numpy())
        s0[long] = np.maximum(c - max_scans // 2, 0)
        s1[long] = np.minimum(s0[long] + max_scans - 1, run.n_scans - 1)
    ions = pd.concat([r[["compound", "anchor", "height_a"]].set_axis(["compound", "trace", "h"], axis=1),
                      r[["compound", "partner", "height_b"]].set_axis(["compound", "trace", "h"], axis=1)])
    ions = ions.groupby(["compound", "trace"], sort=True)["h"].max().reset_index()
    ions = ions.sort_values(["compound", "h"], ascending=[True, False])
    ions = ions.groupby("compound").head(max_ions).sort_values(["compound", "trace"]).reset_index(drop=True)
    cidx = np.searchsorted(comp_ids, ions["compound"].to_numpy())
    ion_off = np.zeros(len(comp_ids) + 1, np.int64)
    np.add.at(ion_off, cidx + 1, 1)
    ion_off = np.cumsum(ion_off)
    L = (s1 - s0 + 1).astype(np.int64)
    n_per = np.diff(ion_off)
    dat_off = np.concatenate([[0], np.cumsum(n_per * L)]).astype(np.int64)
    return comp_ids, s0.astype(np.int64), L, ions, ion_off, dat_off


def coalescence_ppm(run, q: float = 0.002, max_scans: int = 300, seed: int = 0) -> float:
    """Closest spacing (ppm) at which two centroids still appear separately in one scan: below it, ions coalesce.

    A low quantile of within-scan neighbour spacings, scaled to m/z 400 assuming an Orbitrap-like resolution
    (FWHM in ppm growing with sqrt(m/z)); evaluated at the median m/z of the run's centroids.
    """
    rng = np.random.default_rng(seed)
    scans = rng.choice(run.n_scans, min(max_scans, run.n_scans), replace=False) if run.n_scans else []
    d, m = [], []
    for s in scans:
        a, b = run.offsets[s], run.offsets[s + 1]
        mz = run.mz[a:b]
        if mz.size > 1:
            d.append(np.diff(mz) / mz[:-1] * 1e6)
            m.append(mz[:-1])
    if not d:
        return 10.0
    d, m = np.concatenate(d), np.concatenate(m)
    k = float(np.quantile(d / np.sqrt(m / 400.0), q))
    return float(np.clip(k * np.sqrt(np.median(run.mz) / 400.0), 1.0, 30.0))


def fit_components(rel: pd.DataFrame, traces, sigma: float | None = None, pad_sd: float = 4.0, cv: float | None = None,
                   use_mz: bool = True, split_penalty: float = 1.0, min_split: float = 0.0, max_ions: int = 24,
                   max_scans: int = 600, gate: float = 4.0, sys_ppm: float = 0.3, n_iter: int = 40,
                   try_split: bool = True, max_dm_ppm: float | None = None, min_sep_sd: float = 0.75) -> Components:
    """Fit every relation cluster (scan-graph compound) of one run with one and two shared-shape EMG components.

    ``rel`` is the run's relation table (:func:`macha.align.edges.relation_table`); a compound's ions are the
    traces its confident relations join, its window the relations' scans +/- ``pad_sd`` elution sd (``sigma``,
    default the scan graph's learned value stored in ``rel.attrs`` or 2 cycles).  ``cv`` is the relative intensity
    error (estimated from high-intensity residuals of the one-component fits when ``None``).  A split is kept when
    the chi2 gain exceeds ``split_penalty`` x (extra parameters) x log(observations) + ``min_split``.
    ``max_dm_ppm`` is the closest two ions can be and still appear as separate centroids, the scale of the prior on
    the m/z difference of two coalesced sources (estimated from the run's own centroid spacing when ``None``).
    Two sources closer than ``min_sep_sd`` peak sd stay one component, flagged ``merged``: they are detected, but
    their amplitudes and source m/z are not identified separately.
    """
    run = traces.run
    if max_dm_ppm is None:
        max_dm_ppm = coalescence_ppm(run)
    if sigma is None:
        sigma = float(rel.attrs.get("sigma", 2.0 * run.cycle_time))
    comp_ids, s_lo, L, ions, ion_off, dat_off = _layout(rel, traces, sigma, pad_sd, max_ions, max_scans)
    Q = len(comp_ids)
    ion_trace = ions["trace"].to_numpy(np.int64)
    Y, MZ, OBS = _gather(traces.order.astype(np.int64), traces.offsets.astype(np.int64),
                         run.scan_of_point.astype(np.int64), run.mz, run.intensity.astype(np.float64), ion_trace,
                         ion_off, s_lo, L, dat_off, int(dat_off[-1]))
    # noise level and censoring floor of each ion: interference in its own m/z window over the compound's scans
    q_of_ion = np.repeat(np.arange(Q), np.diff(ion_off))
    hw = floor_halfwidth_ppm(traces, gate, sys_ppm)
    mz_ion = traces.table["mz"].to_numpy()[ion_trace]
    wn = window_noise(traces, mz_ion, s_lo[q_of_ion], s_lo[q_of_ion] + L[q_of_ion] - 1, hw, ion_trace)
    run_floor = float(np.quantile(run.intensity, 0.05)) if run.n_points else 1.0
    noise = np.exp(wn["noise_log_mean"].fillna(np.log(run_floor)).to_numpy())
    noise = np.maximum(noise, 0.25 * run_floor)
    floor = noise.copy()
    # m/z weights from the precision model at each centroid's intensity; centred log intensity per ion
    U = np.zeros_like(Y)
    LI = np.zeros_like(Y)
    pm = traces.precision
    ok = OBS & np.isfinite(MZ)
    sd = pm.sigma_ppm(np.maximum(Y[ok], 1.0)) * 1e-6 * MZ[ok]
    U[ok] = 1.0 / np.maximum(sd, 1e-9) ** 2
    ion_of_dat = np.repeat(np.arange(len(ions)), L[q_of_ion])
    li = np.log(np.maximum(Y, 1.0))
    sums = np.bincount(ion_of_dat[ok], li[ok], len(ions))
    cnts = np.maximum(np.bincount(ion_of_dat[ok], minlength=len(ions)), 1)
    LI[ok] = li[ok] - (sums / cnts)[ion_of_dat[ok]]
    # starting apex: signal-weighted centre of each window
    rt_dat = np.concatenate([run.rt[s:s + n] for s, n in zip(s_lo, L)]) if Q else np.zeros(0)
    rt_rep = np.repeat(np.arange(Q), L)  # window id per scan-slot
    w_sig = np.zeros(int(L.sum()))
    slot = np.concatenate([np.arange(n) for n in L]) if Q else np.zeros(0, int)
    slot_off = np.concatenate([[0], np.cumsum(L)])
    for q in range(Q):
        blk = Y[dat_off[q]:dat_off[q + 1]].reshape(-1, L[q]).sum(0)
        w_sig[slot_off[q]:slot_off[q + 1]] = blk ** 2
    mu0 = np.bincount(rt_rep, w_sig * rt_dat, Q) / np.maximum(np.bincount(rt_rep, w_sig, Q), 1e-300)
    sig0 = np.full(Q, float(sigma))
    rt_all = run.rt.astype(np.float64)

    def run_fit(cv_, try2):
        return _fit_all(rt_all, s_lo, L, ion_off, dat_off, Y, MZ, OBS, U, LI, noise, floor, mu0, sig0, cv_, use_mz,
                        0.5 * max_dm_ppm, try2, n_iter)

    if cv is None:
        res = run_fit(0.1, False)
        HI = res[-2]
        good = HI[:, 1] >= 5
        cv = float(np.sqrt(np.median(HI[good, 0] / HI[good, 1]))) if good.sum() >= 20 else 0.1
        cv = float(np.clip(cv, 0.03, 0.5))
    P1, P2, C1, C2, NOBS, A1, A2, M1, M2, HI, B = run_fit(cv, try_split)
    del slot

    # split decision: chi2 gain against extra parameters (second apex, a second amplitude and m/z per split ion)
    chi1 = C1.sum(1)
    chi2 = C2.sum(1)
    two_ion = (A2[:, 0] > 0) & (A2[:, 1] > 0)
    n_two = np.bincount(q_of_ion, two_ion.astype(float), Q)
    n_mz = np.bincount(q_of_ion, (np.nan_to_num(M2[:, 2]) > 0).astype(float), Q)
    extra = 1.0 + n_two + n_mz
    gain = chi1 - chi2
    pen = split_penalty * extra * np.log(np.maximum(NOBS, 2)) + min_split
    # a plateau (chemical background) is no peak: its one-component fit runs into the width bounds or out of the
    # window, and any two peaks fit it better, so it is never split
    one = _shape_stats(P1, P2, np.zeros(Q, np.bool_))
    edge = 2.0 * float(run.cycle_time)
    capped = L >= max_scans
    plateau = (capped | (P1[:, 1] >= 0.95 * 3.0 * sigma) | (P1[:, 2] >= 0.95 * 4.0 * sigma)
               | (one[:, 0, 2] <= run.rt[s_lo] + edge) | (one[:, 0, 3] >= run.rt[s_lo + L - 1] - edge))
    plateau &= np.isfinite(P1[:, 0])
    split = np.isfinite(gain) & (gain > pen) & ~plateau
    # each side must carry signal above noise in at least one ion
    dens = _peak_density(P1, P2, split)
    ion_noise_q = np.sqrt(np.bincount(q_of_ion, noise ** 2, Q))
    h2 = np.zeros((Q, 2))
    for k in range(2):
        h2[:, k] = np.bincount(q_of_ion, A2[:, k], Q) * np.nan_to_num(dens[:, k])
    split &= (h2[:, 0] >= 3 * ion_noise_q) & (h2[:, 1] >= 3 * ion_noise_q)
    # closer than ~0.75 sd, two sources are detected (mostly by their m/z trajectory) but their amplitude ratio and
    # source m/z trade off against each other: one component, flagged as holding two sources
    merged = split & (np.abs(P2[:, 0] - P2[:, 1]) < min_sep_sd * P2[:, 2])
    split &= ~merged
    dens = _peak_density(P1, P2, split)
    shape = _shape_stats(P1, P2, split)

    # tables
    fits = pd.DataFrame(dict(compound=comp_ids, n_ions=np.diff(ion_off), n_obs=NOBS, scans=L, rt_lo=run.rt[s_lo],
                             rt_hi=run.rt[s_lo + L - 1],
                             chi2_int_1=C1[:, 0], chi2_mz_1=C1[:, 1], chi2_int_2=C2[:, 0], chi2_mz_2=C2[:, 1],
                             gain=gain, gain_int=C1[:, 0] - C2[:, 0], gain_mz=C1[:, 1] - C2[:, 1], penalty=pen,
                             plateau=plateau, split=split, merged=merged))
    rows = []
    ion_rows = []
    ion_tab = traces.table
    for k in range(2):
        has = split if k == 1 else np.isfinite(P1[:, 0])
        q = np.flatnonzero(has)
        amp_ion = np.where(split[q_of_ion], A2[:, k], A1 if k == 0 else 0.0)
        area = np.bincount(q_of_ion, amp_ion, Q)
        height = area * np.nan_to_num(dens[:, k])
        mu = np.where(split, P2[:, k], P1[:, 0] if k == 0 else np.nan)
        sg = np.where(split, P2[:, 2], P1[:, 1])
        ta = np.where(split, P2[:, 3], P1[:, 2])
        rows.append(pd.DataFrame(dict(compound=comp_ids[q], k=k, rt=shape[q, k, 0], mu=mu[q], sigma=sg[q], tau=ta[q],
                                      fwhm=shape[q, k, 1], rt_start=shape[q, k, 2], rt_end=shape[q, k, 3],
                                      area=area[q], height=height[q], snr=height[q] / ion_noise_q[q],
                                      n_ions=np.bincount(q_of_ion, (amp_ion > 0).astype(float), Q)[q].astype(int),
                                      plateau=plateau[q], split=split[q], merged=merged[q])))
        src = np.where(split[q_of_ion], M2[:, 0] if k == 0 else M2[:, 1], M1 if k == 0 else np.nan)
        tot = np.where(split[q_of_ion], A2.sum(1), A1)
        sel = has[q_of_ion] & (amp_ion > 0)
        ion_rows.append(pd.DataFrame(dict(compound=comp_ids[q_of_ion[sel]], k=k, trace=ion_trace[sel],
                                          mz=mz_ion[sel], mz_source=src[sel], area=amp_ion[sel],
                                          height=amp_ion[sel] * np.nan_to_num(dens[q_of_ion[sel], k]),
                                          share=amp_ion[sel] / np.maximum(tot[sel], 1e-300), noise=noise[sel],
                                          baseline=np.where(split[q_of_ion], B[:, 1], B[:, 0])[sel],
                                          coalesced=(split[q_of_ion] & (np.nan_to_num(M2[:, 2]) > 0))[sel],
                                          apex_intensity=ion_tab["apex_intensity"].to_numpy()[ion_trace[sel]])))
    table = pd.concat(rows, ignore_index=True).sort_values(["compound", "k"]).reset_index(drop=True)
    table.insert(0, "component", np.arange(len(table)))
    ions_df = pd.concat(ion_rows, ignore_index=True)
    key = pd.Series(table["component"].to_numpy(), index=pd.MultiIndex.from_arrays([table["compound"], table["k"]]))
    ions_df.insert(0, "component", key.reindex(pd.MultiIndex.from_arrays([ions_df["compound"], ions_df["k"]])).to_numpy())
    ions_df = ions_df.sort_values(["component", "area"], ascending=[True, False]).reset_index(drop=True)
    relc = _relation_components(rel, ions_df, table)
    table.attrs.update(cv=cv, sigma=sigma)
    return Components(table, ions_df, fits, relc, cv)


def _relation_components(rel, ions, table):
    """Component of each relation: for a split compound, the side holding most of both ions' signal."""
    out = np.full(len(rel), -1, np.int64)
    comp = rel["compound"].to_numpy()
    first = table.drop_duplicates("compound").set_index("compound")["component"]
    has = np.isin(comp, first.index.to_numpy())
    out[has] = first.loc[comp[has]].to_numpy()
    sp = table[table["split"]]
    if not len(sp):
        return out
    share = ions.set_index(["compound", "trace", "k"])["share"]
    sel = np.flatnonzero(np.isin(comp, sp["compound"].to_numpy()))
    a = rel["anchor"].to_numpy()[sel]
    b = rel["partner"].to_numpy()[sel]
    c = comp[sel]
    s1 = np.zeros(sel.size)
    for t in (a, b):
        idx = pd.MultiIndex.from_arrays([c, t, np.ones(sel.size, int)])
        s1 += np.nan_to_num(share.reindex(idx).to_numpy())
    k = (s1 > 1.0).astype(int)
    key = pd.Series(table["component"].to_numpy(), index=pd.MultiIndex.from_arrays([table["compound"], table["k"]]))
    out[sel] = key.reindex(pd.MultiIndex.from_arrays([c, k])).to_numpy()
    return out


# ------------------------------------------------------------------------------------------------ across runs
@dataclass
class CrossRun:
    """Components linked across runs, and their targeted quantification in every run.

    ``table``: one row per cross-run component (reference apex, median shape, runs where it was fitted, the
    relation groups behind it, ``web_consistency``: the mean share of them present in a run that found it, and
    ``anchor_group``: the strongest of them, whose two ions quantify it).
    ``members``: (xc, run, component) for the runs where the per-run fit found it.  ``ions``: the ions used to
    quantify each cross-run component (m/z, relative amplitude ``ratio``, share of the ion's signal).  The
    matrices (components x runs) come from :func:`quantify`: summed area of the quantification ions, apex, the
    likelihood-ratio evidence for presence over the ions' own noise floor, and the noise-scaled height.
    ``area`` is the fitted amplitude (a matched filter: best at low SNR); ``area_int`` integrates the observed
    signal within the fitted peak (model where a scan is censored), which a shape misfit does not bias.
    """

    table: pd.DataFrame
    members: pd.DataFrame
    ions: pd.DataFrame
    siblings: pd.DataFrame | None = None  # xc pairs that were the two sides of a split in some run
    area: np.ndarray | None = None
    area_int: np.ndarray | None = None
    rt: np.ndarray | None = None
    llr: np.ndarray | None = None
    snr: np.ndarray | None = None


@nb.njit(cache=True)
def _close_pairs(mz, grp, t, tol_t, tol_ppm):
    """Pairs of entries (sorted by m/z) within ``tol_ppm`` of each other, of different ``grp``, whose times differ by
    at most the larger of their tolerances: a windowed join, so dense m/z ranges never form long chains."""
    n = mz.shape[0]
    cnt = 0
    for pass_ in range(2):
        if pass_ == 1:
            a = np.empty(cnt, np.int64)
            b = np.empty(cnt, np.int64)
            cnt = 0
        for i in range(n):
            lim = mz[i] * (1.0 + tol_ppm * 1e-6)
            j = i + 1
            while j < n and mz[j] <= lim:
                if grp[j] != grp[i] and abs(t[i] - t[j]) <= max(tol_t[i], tol_t[j]):
                    if pass_ == 1:
                        a[cnt] = min(grp[i], grp[j])
                        b[cnt] = max(grp[i], grp[j])
                    cnt += 1
                j += 1
    return a, b


def class_spread(model, t) -> np.ndarray:
    """Largest (over runs) sd of the class warps around their weighted mean at reference times ``t``: how far a
    compound can sit from the mean-warp prediction before its drift class is known."""
    if model is None:
        return np.zeros(np.size(t))
    W = model.warps(np.atleast_1d(np.asarray(t, float)))  # (K, J, n)
    mean = np.einsum("k,kjn->jn", model.pi, W)
    var = np.einsum("k,kjn->jn", model.pi, (W - mean[None]) ** 2)
    return np.sqrt(var.max(0))


def _merge_fragments(N, comps, model, mz_tol_ppm, min_share, n_sd: float = 4.0):
    """Join cross-run components found in disjoint sets of runs that share their ions and agree in reference time.

    Relation groups (and so the links between components) break where a compound's drift class departs from the
    mean warp; the pieces then cover complementary runs.  Two pieces are the same compound when no run holds both,
    at least two of their ions (the only one, for single-ion pieces) match within ``mz_tol_ppm``, and their
    reference apexes differ by less than ``n_sd`` x the drift model's residual and class spread.
    """
    from .align.edges import to_reference

    xc = N["xc"].to_numpy().copy()
    ref = np.full(len(N), np.nan)
    ion_rows = []
    for j, c in enumerate(comps):
        m = N["run"].to_numpy() == j
        comp = N["component"].to_numpy()[m]
        rt = c.table["rt"].to_numpy()[comp]
        ref[m] = to_reference(model, rt, j) if model is not None and len(rt) else rt
        io = c.ions[c.ions["component"].isin(comp)]
        io = io[io["share"].to_numpy() >= min_share]
        node = pd.Series(N["node"].to_numpy()[m], index=comp)
        ion_rows.append(pd.DataFrame(dict(node=node.loc[io["component"].to_numpy()].to_numpy(),
                                          mz=io["mz"].to_numpy())))
    I = pd.concat(ion_rows, ignore_index=True)
    X = int(xc.max()) + 1 if len(xc) else 0
    node_x = pd.Series(xc, index=N["node"].to_numpy())
    I["xc"] = node_x.loc[I["node"].to_numpy()].to_numpy()
    if not len(I):
        return xc
    # one entry per ion of a piece (its members in several runs list the same ions)
    I = I.sort_values(["xc", "mz"]).reset_index(drop=True)
    m_, x_ = I["mz"].to_numpy(), I["xc"].to_numpy()
    new = np.concatenate([[True], (x_[1:] != x_[:-1]) | (np.diff(m_) / m_[:-1] * 1e6 > mz_tol_ppm)])
    I = I.groupby(np.cumsum(new)).agg(xc=("xc", "first"), mz=("mz", "median")).sort_values("mz").reset_index(drop=True)
    xref = pd.Series(ref).groupby(xc).median().reindex(np.arange(X)).to_numpy()
    sig = float(model.sigma) if model is not None else 2.0
    tol_x = n_sd * np.sqrt(sig ** 2 + class_spread(model, np.nan_to_num(xref)) ** 2)
    xi = I["xc"].to_numpy(np.int64)
    a, b = _close_pairs(I["mz"].to_numpy(float), xi, np.nan_to_num(xref[xi], nan=-1e9), tol_x[xi], mz_tol_ppm)
    if not a.size:
        return xc
    n_ion = np.bincount(xi, minlength=X)
    P = pd.DataFrame(dict(xc_a=a, xc_b=b)).groupby(["xc_a", "xc_b"]).size().rename("shared").reset_index()
    need = np.minimum(2, np.minimum(n_ion[P["xc_a"]], n_ion[P["xc_b"]]))
    P = P[P["shared"].to_numpy() >= need]
    runs = pd.DataFrame(dict(xc=xc, run=N["run"].to_numpy()))
    mask = runs.groupby("xc")["run"].apply(lambda r: int(np.sum(2 ** np.unique(r).astype(np.int64)))).reindex(
        np.arange(X), fill_value=0).to_numpy()
    d = np.abs(xref[P["xc_a"]] - xref[P["xc_b"]])
    P = P.assign(d=d)[(mask[P["xc_a"]] & mask[P["xc_b"]]) == 0].sort_values("d")
    parent = np.arange(X)
    pm = mask.copy()

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in zip(P["xc_a"].to_numpy(), P["xc_b"].to_numpy()):
        ra, rb = find(a), find(b)
        if ra != rb and (pm[ra] & pm[rb]) == 0:
            lo, hi = min(ra, rb), max(ra, rb)
            parent[hi] = lo
            pm[lo] |= pm[hi]
    return np.array([find(x) for x in xc])


def _relation_ions(ions, M, N, groups, mz_tol_ppm):
    """Keep, per cross-run component, the ions matching both ions of its strongest relation group (in most of the
    component's runs, then highest median SNR); components without a match keep their ions."""
    node_xc = pd.Series(N["xc"].to_numpy(), index=N["node"].to_numpy())
    G = M.assign(xc=node_xc.loc[M["node"].to_numpy()].to_numpy())[["xc", "group", "run"]].drop_duplicates()
    snr = np.nanmedian(groups.snr, 1)
    G = G.groupby(["xc", "group"]).size().rename("runs").reset_index()
    G["snr"] = snr[G["group"].to_numpy()]
    best = G.sort_values(["xc", "runs", "snr"], ascending=[True, False, False]).drop_duplicates("xc")
    tab = groups.table
    best = best.assign(mz_a=tab["mz_a"].to_numpy()[best["group"].to_numpy()],
                       mz_b=tab["mz_b"].to_numpy()[best["group"].to_numpy()]).set_index("xc")
    has = ions["xc"].isin(best.index).to_numpy()
    ma = np.full(len(ions), np.nan)
    mb = np.full(len(ions), np.nan)
    ma[has] = best.loc[ions["xc"].to_numpy()[has], "mz_a"].to_numpy()
    mb[has] = best.loc[ions["xc"].to_numpy()[has], "mz_b"].to_numpy()
    mz = ions["mz"].to_numpy()
    hit = (np.abs(mz - ma) / mz * 1e6 <= mz_tol_ppm) | (np.abs(mz - mb) / mz * 1e6 <= mz_tol_ppm)
    n_hit = pd.Series(hit).groupby(ions["xc"].to_numpy()).transform("sum").to_numpy()
    out = ions[hit | (n_hit == 0)]
    return out.assign(exclusive=out["exclusive"] | (n_hit[hit | (n_hit == 0)] > 0)), best["group"]


def predict_rt(model, T: np.ndarray):
    """Posterior-mean RT in every run from the runs where each row was observed (NaN = missing), and its sd.

    The drift class of a row is unknown: each class's consensus time is fitted to the observed runs, classes are
    weighted by their fit, and the prediction is the class-weighted mean; its sd adds the spread between classes
    to the model's residual sd.
    """
    T = np.asarray(T, float)
    C, J = T.shape
    obs = np.isfinite(T)
    n = np.maximum(obs.sum(1), 1)
    logpi = np.log(np.maximum(model.pi, 1e-12))
    pred = np.empty((model.K, C, J))
    ll = np.empty((model.K, C))
    for k in range(model.K):
        mu = np.where(obs, T, 0.0).sum(1) / n
        for _ in range(5):
            W = model.warps(mu)[k]
            mu = np.where(obs, T - W.T, 0.0).sum(1) / n
        W = model.warps(mu)[k]
        r = np.where(obs, T - mu[:, None] - W.T, 0.0)
        ll[k] = logpi[k] - 0.5 * (r * r).sum(1) / model.sigma ** 2
        pred[k] = mu[:, None] + W.T
    ll -= ll.max(0, keepdims=True)
    resp = np.exp(ll)
    resp /= resp.sum(0, keepdims=True)
    mean = np.einsum("kc,kcj->cj", resp, pred)
    var = np.einsum("kc,kcj->cj", resp, (pred - mean[None]) ** 2)
    return mean, np.sqrt(model.sigma ** 2 + var)


def link_components(comps: list[Components], groups, model=None, min_runs: int = 2, max_ions: int = 4,
                    mz_tol_ppm: float = 6.0, min_share: float = 0.8, quant_ions: str = "relation") -> CrossRun:
    """Link the runs' components through the relation groups they share.

    Every relation group (:func:`macha.align.edges.relation_groups`: one relation matched across runs in reference
    time) votes for the pairs of components its members belong to.  Two components of different runs are linked
    when each is the other's best-voted partner; linked components form a cross-run component, keeping one
    component per run (the best-voted).  Quantification ions are the member ions present in most runs that the
    component holds alone (share >= ``min_share``, so not the coalesced ion of a split), ranked by amplitude.

    ``quant_ions="relation"`` (default) quantifies from the two ions of the component's strongest relation group
    (present in most of its runs, then highest SNR), which certainly belong together; a scan-graph compound in
    dense data can hold several co-eluting analytes, and its most abundant ions need not be one of them.
    ``"top"`` uses up to ``max_ions`` of the component's ions as described above.
    """
    from .align.edges import to_reference

    J = len(comps)
    M = groups.members.copy()
    M["comp"] = -1
    for j in range(J):
        m = M["run"].to_numpy() == j
        M.loc[m, "comp"] = comps[j].relation_component[M.loc[m, "rel"].to_numpy()]
    M = M[M["comp"].to_numpy() >= 0]
    off = np.concatenate([[0], np.cumsum([len(c.table) for c in comps])])
    M["node"] = off[M["run"].to_numpy()] + M["comp"].to_numpy()
    A = M[["group", "run", "node"]]
    V = A.merge(A, on="group", suffixes=("_a", "_b"))
    V = V[V["run_a"].to_numpy() < V["run_b"].to_numpy()]
    V = V.groupby(["run_a", "node_a", "run_b", "node_b"]).size().rename("votes").reset_index()
    best_b = V.sort_values("votes", ascending=False).drop_duplicates(["node_a", "run_b"])
    best_a = V.sort_values("votes", ascending=False).drop_duplicates(["node_b", "run_a"])
    L = best_b.merge(best_a[["node_a", "node_b"]], on=["node_a", "node_b"])
    n = int(off[-1])
    parent = np.arange(n)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in zip(L["node_a"].to_numpy(), L["node_b"].to_numpy()):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    root = np.array([find(x) for x in range(n)])
    run_of = np.repeat(np.arange(J), np.diff(off))
    # one component per run: keep the best-linked, the others become cross-run components of their own
    deg = np.bincount(np.concatenate([L["node_a"], L["node_b"]]).astype(np.int64),
                      np.concatenate([L["votes"], L["votes"]]).astype(float), n)
    N = pd.DataFrame(dict(node=np.arange(n), root=root, run=run_of, deg=deg))
    N = N.sort_values(["root", "run", "deg"], ascending=[True, True, False])
    dup = N.duplicated(["root", "run"]).to_numpy()
    N.loc[N.index[dup], "root"] = n + N.loc[N.index[dup], "node"]  # fresh ids: a cluster's root may be detached
    N = N.sort_values("node")
    N["xc"] = pd.factorize(N["root"], sort=True)[0]
    N["component"] = N["node"].to_numpy() - off[N["run"].to_numpy()]
    N["xc"] = _merge_fragments(N, comps, model, mz_tol_ppm, min_share)
    members = N[["xc", "run", "component"]].reset_index(drop=True)
    # consensus apex (reference time) and shape
    parts = []
    for j, c in enumerate(comps):
        t = c.table
        mj = members[members["run"] == j]
        r = t.set_index("component").loc[mj["component"].to_numpy()]
        rr = r["rt"].to_numpy()
        ref = to_reference(model, rr, j) if model is not None and len(rr) else rr
        parts.append(pd.DataFrame(dict(xc=mj["xc"].to_numpy(), run=j, rt=rr, rt_ref=ref, sigma=r["sigma"].to_numpy(),
                                       tau=r["tau"].to_numpy(), snr=r["snr"].to_numpy(), area=r["area"].to_numpy(),
                                       plateau=r["plateau"].to_numpy(), split=r["split"].to_numpy())))
    P = pd.concat(parts, ignore_index=True)
    table = P.groupby("xc").agg(n_runs=("run", "size"), rt_ref=("rt_ref", "median"), sigma=("sigma", "median"),
                                tau=("tau", "median"), snr=("snr", "median"), plateau=("plateau", "mean"),
                                split=("split", "max")).reset_index()
    table = table[table["n_runs"].to_numpy() >= min_runs].reset_index(drop=True)
    # the relation web: relation groups behind each component, and the share of them present in each run found
    node_xc = pd.Series(N["xc"].to_numpy(), index=N["node"].to_numpy())
    Wb = M.assign(xc=node_xc.loc[M["node"].to_numpy()].to_numpy())[["xc", "run", "group"]].drop_duplicates()
    G = Wb.groupby("xc")["group"].nunique()
    per = Wb.groupby(["xc", "run"])["group"].nunique()
    cons = (per / G.reindex(per.index.get_level_values(0)).to_numpy()).groupby(level=0).mean()
    table["relation_groups"] = G.reindex(table["xc"]).fillna(0).astype(int).to_numpy()
    table["web_consistency"] = cons.reindex(table["xc"]).to_numpy()
    keep = set(table["xc"].tolist())
    members = members[members["xc"].isin(keep)]
    # quantification ions: member ions clustered by m/z across runs
    ip = []
    for j, c in enumerate(comps):
        mj = members[members["run"] == j]
        io = c.ions[c.ions["component"].isin(mj["component"])]
        tot = io.groupby("component")["area"].transform("sum")
        io = io.assign(rel=io["area"] / tot.replace(0, np.nan), run=j,
                       xc=mj.set_index("component").loc[io["component"].to_numpy(), "xc"].to_numpy())
        ip.append(io[["xc", "run", "mz", "rel", "share", "area"]])
    I = pd.concat(ip, ignore_index=True).sort_values(["xc", "mz"]).reset_index(drop=True)
    mz = I["mz"].to_numpy()
    xc = I["xc"].to_numpy()
    new = np.ones(len(I), bool)
    new[1:] = (xc[1:] != xc[:-1]) | ((mz[1:] - mz[:-1]) / mz[:-1] * 1e6 > mz_tol_ppm)
    I["ion"] = np.cumsum(new) - 1
    ions = I.groupby("ion").agg(xc=("xc", "first"), mz=("mz", "median"), ratio=("rel", "median"),
                                share=("share", "median"), n_runs=("run", "nunique")).reset_index(drop=True)
    nr = table.set_index("xc")["n_runs"]
    ions["frac_runs"] = ions["n_runs"] / nr.loc[ions["xc"]].to_numpy()
    ions["exclusive"] = ions["share"] >= min_share
    anchor = None
    if quant_ions == "relation":
        ions, anchor = _relation_ions(ions, M, N, groups, mz_tol_ppm)
    ions = ions[ions["frac_runs"] >= 0.5]
    # a coalesced ion carries a neighbour's signal: quantify from the ions the component holds alone if it has any
    has_ex = ions.groupby("xc")["exclusive"].transform("any").to_numpy()
    ions = ions[ions["exclusive"].to_numpy() | ~has_ex]
    ions = ions.sort_values(["xc", "frac_runs", "ratio"], ascending=[True, False, False])
    ions = ions.groupby("xc").head(max_ions).reset_index(drop=True)
    ions["ratio"] = ions["ratio"] / ions.groupby("xc")["ratio"].transform("sum")
    table = table[table["xc"].isin(set(ions["xc"]))].reset_index(drop=True)
    # siblings: the two sides of a split in some run (distinct co-eluting sources sharing an ion)
    sib = []
    for j, c in enumerate(comps):
        sp = c.table[c.table["split"].to_numpy()]
        if not len(sp):
            continue
        node_xc = pd.Series(N["xc"].to_numpy(), index=N["node"].to_numpy())
        pr = sp.pivot_table(index="compound", columns="k", values="component", aggfunc="first").dropna()
        if len(pr):
            a = node_xc.loc[off[j] + pr[0].to_numpy().astype(np.int64)].to_numpy()
            b = node_xc.loc[off[j] + pr[1].to_numpy().astype(np.int64)].to_numpy()
            sib.append(pd.DataFrame(dict(xc_a=np.minimum(a, b), xc_b=np.maximum(a, b))))
    # dense ids
    ids = {x: k for k, x in enumerate(table["xc"].tolist())}
    table["anchor_group"] = anchor.reindex(table["xc"]).fillna(-1).astype(np.int64).to_numpy() if anchor is not None \
        else -1
    table["xc"] = np.arange(len(table))
    members = members[members["xc"].isin(ids)]
    members = members.assign(xc=members["xc"].map(ids).astype(np.int64)).reset_index(drop=True)
    ions = ions.assign(xc=ions["xc"].map(ids)).sort_values(["xc", "ratio"], ascending=[True, False]).reset_index(drop=True)
    S = pd.concat(sib, ignore_index=True) if sib else pd.DataFrame(dict(xc_a=[], xc_b=[]))
    S = S[S["xc_a"].isin(ids) & S["xc_b"].isin(ids) & (S["xc_a"] != S["xc_b"])]
    S = S.assign(xc_a=S["xc_a"].map(ids), xc_b=S["xc_b"].map(ids)).groupby(["xc_a", "xc_b"]).size()
    siblings = S.rename("runs").reset_index()
    return CrossRun(table, members, ions, siblings=siblings)


@nb.njit(cache=True, parallel=True)
def _xic(offsets, mz, intensity, q_mz, q_s0, q_len, q_off, pa, pb, pc, gate, sys_ppm, hw_max, closest):
    """Per query (ion window): one centroid per scan within the intensity-dependent m/z tolerance
    ``gate * sqrt(sigma_ppm(I)^2 + sys^2)`` (the tracker's acceptance) and its m/z: the closest in units of its
    tolerance (``closest``), else the most intense."""
    n = q_off[-1]
    Y = np.zeros(n)
    MZ = np.full(n, np.nan)
    n_scans = offsets.shape[0] - 1
    for q in nb.prange(q_mz.shape[0]):
        m0 = q_mz[q]
        lo_mz = m0 * (1.0 - hw_max * 1e-6)
        hi_mz = m0 * (1.0 + hw_max * 1e-6)
        for k in range(q_len[q]):
            s = q_s0[q] + k
            if s < 0 or s >= n_scans:
                continue
            a, b = offsets[s], offsets[s + 1]
            lo, hi = a, b
            while lo < hi:
                mid = (lo + hi) // 2
                if mz[mid] < lo_mz:
                    lo = mid + 1
                else:
                    hi = mid
            j = lo
            best = 0.0
            bm = np.nan
            bz = 1e300
            while j < b and mz[j] <= hi_mz:
                inten = intensity[j]
                sp = math.sqrt(pa * pa + (pb * max(inten, 1e-9) ** (-pc)) ** 2)
                tol = gate * math.sqrt(sp * sp + sys_ppm * sys_ppm)
                z = abs(mz[j] - m0) / m0 * 1e6 / tol
                if z <= 1.0 and ((closest and z < bz) or (not closest and inten > best)):
                    best = inten
                    bm = mz[j]
                    bz = z
                j += 1
            Y[q_off[q] + k] = best
            MZ[q_off[q] + k] = bm
    return Y, MZ


@nb.njit(cache=True, parallel=True)
def _fit_fixed(rt_all, x_s0, x_len, x_ion0, x_nion, ratio, Y, y_off, noise, base0, fixed, apex0, half, sd_pred,
               sig, tau, cv, n_grid):
    """Amplitude of a fixed-shape component (ion ratios fixed) over a constant baseline per ion, at the best apex
    within +/- ``half`` of the prediction (Gaussian prior with sd ``sd_pred``), censored at each ion's noise level.
    Returns amplitude, apex, chi2(baselines only) - chi2(amplitude), the number of observed points and the
    integrated area (observed signal above the baseline inside the fitted peak's 1 % bounds, the model where a scan
    is censored).  Where ``fixed`` is set for an ion, its baseline stays at ``base0`` (measured away from the peak)
    instead of being refitted with the amplitude."""
    X = x_s0.shape[0]
    out = np.full((X, 5), np.nan)
    for x in nb.prange(X):
        L = x_len[x]
        n_i = x_nion[x]
        if L <= 0 or not (sig[x] == sig[x]):
            continue
        rt = rt_all[x_s0[x]:x_s0[x] + L]
        off = _emg_mode(0.0, sig[x], tau[x])
        best = 1e300
        nobs = 0
        for i in range(n_i):
            b0 = y_off[x_ion0[x] + i]
            for s in range(L):
                if Y[b0 + s] > 0:
                    nobs += 1
        # null model: baselines only
        chi0 = 0.0
        for i in range(n_i):
            ii = x_ion0[x] + i
            b0 = y_off[ii]
            nl = noise[ii]
            sw = 0.0
            sy = 0.0
            for s in range(L):
                if Y[b0 + s] > 0:
                    sw += 1.0
                    sy += Y[b0 + s]
            bl = sy / sw if sw > 0 else 0.0
            if fixed[ii]:
                bl = base0[ii]
            for s in range(L):
                y = Y[b0 + s]
                w = 1.0 / (nl * nl + cv * cv * bl * bl)
                if y > 0:
                    chi0 += (y - bl) * (y - bl) * w
                elif bl > nl:
                    chi0 += (bl - nl) * (bl - nl) * w
        ng = n_grid if half[x] > 0 else 1
        f = np.empty(L)
        bl = np.empty(n_i)
        best_bl = np.zeros(n_i)
        for gi in range(ng):
            d = 0.0 if ng == 1 else -half[x] + 2.0 * half[x] * gi / (ng - 1)
            apex = apex0[x] + d
            mu = apex - off
            for s in range(L):
                f[s] = _emg(rt[s], mu, sig[x], tau[x])
            c = 0.0
            for i in range(n_i):
                bl[i] = base0[x_ion0[x] + i]
            for it in range(4):
                # amplitude given baselines, then baselines given amplitude (weights from the current model)
                num = 0.0
                den = 0.0
                for i in range(n_i):
                    ii = x_ion0[x] + i
                    b0 = y_off[ii]
                    nl = noise[ii]
                    for s in range(L):
                        y = Y[b0 + s]
                        if y > 0:
                            m = c * ratio[ii] * f[s] + bl[i]
                            w = 1.0 / (nl * nl + cv * cv * (m * m if it > 0 else y * y))
                            num += w * (y - bl[i]) * ratio[ii] * f[s]
                            den += w * ratio[ii] * ratio[ii] * f[s] * f[s]
                c = max(num / den, 0.0) if den > 0 else 0.0
                for i in range(n_i):
                    ii = x_ion0[x] + i
                    if fixed[ii]:
                        continue
                    b0 = y_off[ii]
                    nl = noise[ii]
                    sw = 0.0
                    sy = 0.0
                    for s in range(L):
                        y = Y[b0 + s]
                        if y > 0:
                            m = c * ratio[ii] * f[s] + bl[i]
                            w = 1.0 / (nl * nl + cv * cv * m * m)
                            sw += w
                            sy += w * (y - c * ratio[ii] * f[s])
                    bl[i] = max(sy / sw, 0.0) if sw > 0 else 0.0
            chi = 0.0
            for i in range(n_i):
                ii = x_ion0[x] + i
                b0 = y_off[ii]
                nl = noise[ii]
                for s in range(L):
                    y = Y[b0 + s]
                    m = c * ratio[ii] * f[s] + bl[i]
                    w = 1.0 / (nl * nl + cv * cv * m * m)
                    if y > 0:
                        chi += (y - m) * (y - m) * w
                    elif m > nl:
                        chi += (m - nl) * (m - nl) * w
            pen = (d / sd_pred[x]) ** 2 if sd_pred[x] > 0 else 0.0
            if chi + pen < best:
                best = chi + pen
                out[x, 0] = c
                out[x, 1] = apex
                out[x, 2] = chi0 - chi
                out[x, 3] = nobs
                best_bl[:] = bl
        if not (out[x, 1] == out[x, 1]):
            continue
        mu = out[x, 1] - off
        top = _emg(out[x, 1], mu, sig[x], tau[x])
        area = 0.0
        for s in range(L):
            fs = _emg(rt[s], mu, sig[x], tau[x])
            if fs < 0.01 * top:
                continue
            dt = (rt[min(s + 1, L - 1)] - rt[max(s - 1, 0)]) / (min(s + 1, L - 1) - max(s - 1, 0))
            for i in range(n_i):
                ii = x_ion0[x] + i
                y = Y[y_off[ii] + s]
                area += ((y - best_bl[i]) if y > 0 else out[x, 0] * ratio[ii] * fs) * dt
        out[x, 4] = area
    return out


@nb.njit(cache=True, parallel=True)
def _joint(rt_all, g_off, g_mem, x_s0, x_len, x_ion0, x_nion, ratio, chan, Y, y_off, noise, apex, sig, tau, cv,
           n_sweep):
    """Joint amplitudes of components that share an m/z channel and overlap in time (non-negative coordinate
    descent, weighted by each ion's noise and the relative error).  Every member's extracted chromatograms are
    modelled with the contributions of all members holding an ion in the same channel.  Returns per member the
    amplitude and the integrated area apportioned by the model."""
    out = np.full((g_mem.shape[0], 2), np.nan)
    for g in nb.prange(g_off.shape[0] - 1):
        a, b = g_off[g], g_off[g + 1]
        M = b - a
        c = np.zeros(M)
        off = np.empty(M)
        for u in range(M):
            x = g_mem[a + u]
            off[u] = _emg_mode(0.0, sig[x], tau[x])
        # observations: every member's ions' chromatograms; design columns: the members
        n_obs = 0
        for u in range(M):
            x = g_mem[a + u]
            n_obs += x_nion[x] * x_len[x]
        yv = np.zeros(n_obs)
        nv = np.zeros(n_obs)
        Dm = np.zeros((n_obs, M))
        dt = np.zeros(n_obs)
        own = np.full(n_obs, -1)
        k = 0
        for u in range(M):
            x = g_mem[a + u]
            L = x_len[x]
            for i in range(x_nion[x]):
                ii = x_ion0[x] + i
                for s in range(L):
                    t = rt_all[x_s0[x] + s]
                    yv[k] = Y[y_off[ii] + s]
                    nv[k] = noise[ii]
                    lo = rt_all[x_s0[x] + max(s - 1, 0)]
                    hi = rt_all[x_s0[x] + min(s + 1, L - 1)]
                    dt[k] = (hi - lo) / max(min(s + 1, L - 1) - max(s - 1, 0), 1)
                    own[k] = u
                    for v in range(M):
                        y = g_mem[a + v]
                        r = 0.0
                        for jj in range(x_ion0[y], x_ion0[y] + x_nion[y]):
                            if chan[jj] == chan[ii]:
                                r = ratio[jj]
                                break
                        if r > 0.0:
                            Dm[k, v] = r * _emg(t, apex[y] - off[v], sig[y], tau[y])
                    k += 1
        for u in range(M):
            x = g_mem[a + u]
            c[u] = 0.0
        w = np.empty(n_obs)
        for it in range(3):
            for kk in range(n_obs):
                m = 0.0
                for v in range(M):
                    m += c[v] * Dm[kk, v]
                ref = m if it > 0 else yv[kk]
                w[kk] = 1.0 / (nv[kk] * nv[kk] + cv * cv * ref * ref)
            for sweep in range(n_sweep):
                for v in range(M):
                    num = 0.0
                    den = 0.0
                    for kk in range(n_obs):
                        if yv[kk] <= 0:
                            continue
                        d = Dm[kk, v]
                        if d == 0.0:
                            continue
                        rest = 0.0
                        for q in range(M):
                            if q != v:
                                rest += c[q] * Dm[kk, q]
                        num += w[kk] * d * (yv[kk] - rest)
                        den += w[kk] * d * d
                    c[v] = max(num / den, 0.0) if den > 0 else 0.0
        for u in range(M):
            out[a + u, 0] = c[u]
        # integrated area: observed signal within each member's own peak, apportioned by the model
        for u in range(M):
            x = g_mem[a + u]
            area = 0.0
            top = _emg(apex[x], apex[x] - off[u], sig[x], tau[x])
            for kk in range(n_obs):
                if own[kk] != u:
                    continue
                fu = Dm[kk, u]
                if fu <= 0.0:
                    continue
                s_rt_f = fu  # ratio * density at this scan
                tot = 0.0
                for q in range(M):
                    tot += c[q] * Dm[kk, q]
                mine = c[u] * fu
                # skip scans outside the 1 % bounds of the member's profile (density relative to its apex)
                if s_rt_f < 0.01 * top * _max_ratio(ratio, x_ion0[x], x_nion[x]):
                    continue
                if yv[kk] > 0 and tot > 0:
                    area += yv[kk] * mine / tot * dt[kk]
                else:
                    area += mine * dt[kk]
            out[a + u, 1] = area
    return out


@nb.njit(cache=True)
def _max_ratio(ratio, i0, n):
    m = 0.0
    for i in range(i0, i0 + n):
        if ratio[i] > m:
            m = ratio[i]
    return m


def _joint_groups(I, apex, sig, tau, chan_ppm, siblings):
    """Sibling components (the two sides of a split in some run) that share an m/z channel (quantification ions
    within ``chan_ppm``, close enough to coalesce) and overlap in time in this run, as CSR (offsets, members).

    Only siblings: components with the same ions and apex are usually one compound the graph cut in two, and
    fitting those jointly would divide its signal arbitrarily."""
    mzs = I["mz"].to_numpy(float)
    order = np.argsort(mzs, kind="stable")
    xi = I["xc"].to_numpy(np.int64)[order]
    reach = 3.0 * sig + 2.0 * tau
    t = np.where(np.isfinite(apex), apex, -1e9)
    a, b = _close_pairs(mzs[order], xi, t[xi], np.nan_to_num(reach[xi]), chan_ppm)
    P = pd.DataFrame(dict(xc_a=a, xc_b=b)).drop_duplicates()
    if siblings is None or not len(siblings):
        P = P.iloc[:0]
    else:
        P = P.merge(siblings[["xc_a", "xc_b"]], on=["xc_a", "xc_b"])
    xa, xb = P["xc_a"].to_numpy(np.int64), P["xc_b"].to_numpy(np.int64)
    # m/z channel of every ion: ions within chan_ppm of a neighbour of another sibling component share it
    new = np.concatenate([[True], np.diff(mzs[order]) / mzs[order][:-1] * 1e6 > chan_ppm]) if len(mzs) else np.zeros(0, bool)
    chan = np.empty(len(mzs), np.int64)
    chan[order] = np.cumsum(new) - 1
    X = apex.shape[0]
    parent = np.arange(X)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in zip(xa.tolist(), xb.tolist()):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    root = np.array([find(x) for x in range(X)])
    size = np.bincount(root, minlength=X)
    sel = np.flatnonzero(size[root] >= 2)
    order = sel[np.argsort(root[sel], kind="stable")]
    r = root[order]
    starts = np.flatnonzero(np.concatenate([[True], r[1:] != r[:-1]])) if r.size else np.zeros(0, int)
    g_off = np.concatenate([starts, [r.size]]).astype(np.int64)
    return g_off, order.astype(np.int64), chan


def quantify(cr: CrossRun, j: int, run, precision, apex: np.ndarray, sd_pred: np.ndarray, half: np.ndarray,
             cv: float = 0.1, gate: float = 4.0, sys_ppm: float = 0.3, n_sd: float = 4.0, n_grid: int = 25,
             chan_ppm: float = 6.0, joint: bool = False, fix_baseline: bool = False,
             closest: bool = True) -> pd.DataFrame:
    """Quantify every cross-run component in run ``j`` with its fixed shape and ion ratios.

    ``apex`` is the expected apex in this run (the run's own fit where it found the component, else the
    reference apex mapped through the drift model), searched within +/- ``half`` under a Gaussian prior of sd
    ``sd_pred``.  Each ion's extracted chromatogram takes, per scan, the centroid closest to the ion's m/z (in units
    of the tracker's intensity-dependent tolerance; ``closest=False``: the most intense one inside it, which in
    dense data is often a stronger neighbour) and is modelled over a constant baseline; its baseline starts at the median of
    what the chromatogram holds away from the peak and its noise is the robust scatter there (the ion's own window,
    not the run), and scans without a centroid are censored at that noise level.  The baseline is refitted with the
    amplitude (``fix_baseline``: kept at the level measured away from the peak; worse on real data).

    ``llr`` (chi2 of baselines alone minus chi2 with the component) is the evidence that the component is present.
    Where it is absent, the amplitude measures the noise and interference in its ions' windows, which compresses
    ratios between runs towards 1: ratios should use runs with enough evidence (e.g. ``llr >= 25``, about 5 sd).  Every run is
    quantified the same way, found or not, so the areas of a component compare across runs.

    Components whose quantification ions fall in one m/z channel (within ``chan_ppm``: close enough to coalesce)
    and overlap in time are then fitted jointly at their apexes, so a near-isobaric neighbour's signal in a shared
    chromatogram is assigned to the neighbour (``joint``).
    """
    T, I = cr.table, cr.ions
    X = len(T)
    cyc = float(run.cycle_time)
    sig = T["sigma"].to_numpy(float)
    tau = T["tau"].to_numpy(float)
    reach = half + n_sd * sig + 3.0 * tau
    s0 = np.clip(np.searchsorted(run.rt, apex - reach), 0, run.n_scans - 1)
    s1 = np.clip(np.searchsorted(run.rt, apex + reach), 0, run.n_scans - 1)
    ok = np.isfinite(apex)
    s0 = np.where(ok, s0, 0).astype(np.int64)
    ln = np.where(ok, s1 - s0 + 1, 0).astype(np.int64)
    x_of_ion = I["xc"].to_numpy()
    x_ion0 = np.searchsorted(x_of_ion, np.arange(X)).astype(np.int64)
    x_nion = np.bincount(x_of_ion, minlength=X).astype(np.int64)
    q_len = ln[x_of_ion]
    q_off = np.concatenate([[0], np.cumsum(q_len)]).astype(np.int64)
    floor_i = float(np.quantile(run.intensity, 0.05)) if run.n_points else 1.0
    hw_max = gate * float(np.sqrt(precision.sigma_ppm(floor_i * 0.5) ** 2 + sys_ppm ** 2))
    Y, MZ = _xic(run.offsets.astype(np.int64), run.mz, run.intensity.astype(np.float64), I["mz"].to_numpy(float),
                 s0[x_of_ion], q_len, q_off, float(precision.a), float(precision.b), float(precision.c), gate, sys_ppm,
                 hw_max, closest)
    # noise per ion: median intensity of its chromatogram away from the peak core
    rt_slot = np.concatenate([run.rt[a:a + n] for a, n in zip(s0[x_of_ion], q_len)]) if len(I) else np.zeros(0)
    ion_of = np.repeat(np.arange(len(I)), q_len)
    core = np.abs(rt_slot - apex[x_of_ion][ion_of]) <= (half + 3.0 * sig + 2.0 * tau)[x_of_ion][ion_of]
    side = (Y > 0) & ~core
    # each ion's baseline (median away from the peak) and noise (robust scatter around it, at least the floor)
    noise = np.full(len(I), np.nan)
    base = np.zeros(len(I))
    if side.any():
        sdf = pd.DataFrame(dict(ion=ion_of[side], y=Y[side]))
        med = sdf.groupby("ion")["y"].median()
        mad = (sdf["y"] - med.reindex(sdf["ion"]).to_numpy()).abs().groupby(sdf["ion"]).median() * 1.4826
        n_side = sdf.groupby("ion").size()
        base[med.index.to_numpy()] = med.to_numpy()
        ok = n_side.to_numpy() >= 5
        noise[mad.index.to_numpy()[ok]] = mad.to_numpy()[ok]
        fixed = np.zeros(len(I), np.bool_)
        fixed[n_side.index.to_numpy()[ok]] = fix_baseline
    if not side.any():
        fixed = np.zeros(len(I), np.bool_)
    noise = np.where(np.isfinite(noise), noise, floor_i)
    noise = np.maximum(noise, 0.5 * floor_i)
    hit = np.bincount(ion_of[side], minlength=len(I)) / np.maximum(np.bincount(ion_of[~core], minlength=len(I)), 1)
    R = _fit_fixed(run.rt.astype(np.float64), s0, ln, x_ion0, x_nion, I["ratio"].to_numpy(float), Y,
                   q_off[:-1].copy(), noise, base, fixed, np.nan_to_num(apex), np.nan_to_num(half), sd_pred, sig, tau,
                   cv, n_grid)
    n_joint = 0
    if joint and X:
        ap = np.where(np.isfinite(R[:, 1]), R[:, 1], np.nan)
        g_off, g_mem, chan = _joint_groups(I, ap, sig, tau, chan_ppm, cr.siblings)
        if g_mem.size:
            JR = _joint(run.rt.astype(np.float64), g_off, g_mem, s0, ln, x_ion0, x_nion, I["ratio"].to_numpy(float),
                        chan, Y, q_off[:-1].copy(), noise, np.nan_to_num(ap), sig, tau, cv, 30)
            R[g_mem, 0] = JR[:, 0]
            R[g_mem, 4] = JR[:, 1]
            n_joint = int(g_mem.size)
    dens = np.array([_emg(_emg_mode(0.0, s, t), 0.0, s, t) if np.isfinite(s) else np.nan for s, t in zip(sig, tau)])
    nz = np.sqrt(np.bincount(x_of_ion, (I["ratio"].to_numpy() > 0) * noise ** 2, X))
    height = R[:, 0] * dens * np.bincount(x_of_ion, I["ratio"].to_numpy(), X)
    return pd.DataFrame(dict(area=R[:, 0], area_int=R[:, 4], rt=R[:, 1], llr=R[:, 2], n_obs=R[:, 3], height=height,
                             snr=height / np.maximum(nz, 1e-300),
                             hit_frac=np.bincount(x_of_ion, hit, X) / np.maximum(x_nion, 1))).assign(
        joint=np.isin(np.arange(X), g_mem) if joint and X and n_joint else False)


def predict_apex(cr: CrossRun, model, J: int, comps: list[Components], min_half: float = 1.0,
                 max_half: float = 10.0, found_half: float = 0.5):
    """Expected apex of every cross-run component in every run, with the search half-width and prior sd.

    Where the run's own fit found the component, its fitted apex is searched within +/- ``found_half`` peak sd (a
    relation cluster can hold co-eluting analytes, so its apex need not be that of the quantification ions).
    Elsewhere the apex is predicted from the runs that found it under the drift model, with the component's drift
    class inferred from them (:func:`predict_rt`).
    """
    X = len(cr.table)
    T = np.full((X, J), np.nan)
    for j, c in enumerate(comps):
        mj = cr.members[cr.members["run"] == j]
        T[mj["xc"].to_numpy(), j] = c.table["rt"].to_numpy()[mj["component"].to_numpy()]
    if model is not None:
        apex, sd = predict_rt(model, T)
    else:
        apex = np.repeat(np.nanmedian(T, 1)[:, None], J, 1)
        sd = np.full((X, J), 2.0)
    half = np.clip(3.0 * sd, min_half, max_half)
    found = np.isfinite(T)
    sig = np.repeat(cr.table["sigma"].to_numpy(float)[:, None], J, 1)
    apex[found] = T[found]
    half[found] = found_half * sig[found]
    sd[found] = found_half * sig[found]
    return apex, sd, half
