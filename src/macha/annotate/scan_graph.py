"""Scan-level modulo graph: compound detection before peak detection.

The feature-level graph (:func:`macha.annotate.annotate`) relates picked
peaks.  Here the same mass-difference rules are tested between *centroids*,
scan by scan, so relations are found wherever both ions are measured, with no
chromatographic model:

1. **Instances.**  Every centroid of a trace is an anchor.  For each rule
   target (``n * delta / z`` for additive rules; a neutral-mass map for charge
   states / multimers) the best partner centroid is searched in the same scan
   and, failing that, in the ``bridge`` neighbouring scans (a partner near the
   noise floor is often missing in one scan).  Partners need not belong to a
   trace, so centroids too sparse to form one can be rescued by a relation.
2. **Relation segments.**  Instances of one (anchor trace, target) are cut into
   runs of consecutive supporting scans.  A real isotope / adduct relation
   persists across the anchor's elution with a constant residual, a constant
   intensity ratio and one consistent partner; a coincidence flickers.  Each
   segment's evidence: support and coverage, the residual's inverse-variance
   mean in units of its standard error (shrinks as 1/sqrt(n)), residual
   dispersion, partner consistency, intensity-ratio constancy, anchor/partner
   correlation over the supporting scans, fraction bridged, and for 13C the
   ratio against expectation.
3. **Posteriors.**  The same search is run with decoy rules.  A weighted
   logistic model of target vs decoy segments (decoys weighted by the
   target/decoy rule ratio, so they represent the expected false targets)
   gives ``P = f_T / (f_T + f_false)`` and hence the posterior that a target is
   real, ``(2P - 1) / P``: target-decoy local FDR.
4. **Compounds in graph space.**  Each trace is cut into blocks covering the
   scan ranges of its relations (so isomers sharing one m/z trace at different
   times stay apart).  Blocks are nodes, segments are edges weighted by their
   log-odds, and correlation clustering (GAEC) gives compounds with a time
   extent and a summed compound chromatogram, on which peak detection can run
   if wanted.
5. **Elution-time prior.**  Relations are local, so on their own they let a
   long-lived ion (chemical background) act as a hub joining peaks minutes
   apart.  The partition therefore expects a compound's members to elute
   together, at a scale learned from the run: the bulk of each node's signal
   must lie within one peak width of the relation that links it (an edge is
   discounted by the log of the share it explains), and member signal centres
   scatter around a shared centre by a learned ``tau`` plus their own standard
   error (a Ward penalty inside GAEC).  Peaks near the noise floor need not
   apex together; their centres are simply less precise.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field, replace

import numba as nb
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit

from .. import chem
from .._kernels import csr_from_labels
from .graph import gaec, gaec_coherent
from .rules import PairRule, Rule, affine, default_rules, make_decoys


# ------------------------------------------------------------------------------------------------ targets
def merge_families(T: pd.DataFrame, tol: float = 0.012) -> pd.DataFrame:
    """Merge additive targets of one charge whose steps differ by < ``tol`` Da into one family target.

    At centroid precision the ~2 Da isotope steps (13C2, 18O, 34S, 37Cl, 81Br, 41K) overlap, so testing them
    separately multiplies the chances of a coincidence without adding information.  Each family keeps the
    mid step and widens its error term to cover the members.  Targets and decoys merge separately.
    """
    keep = []
    add = T[T.alpha != 1.0]
    for (z, decoy), d in T[T.alpha == 1.0].groupby(["z", "decoy"]):
        d = d.sort_values("beta")
        b = d.beta.to_numpy()
        cut = np.concatenate([[0], np.flatnonzero(np.diff(b) > tol / z) + 1, [len(d)]])
        for lo, hi in zip(cut[:-1], cut[1:]):
            g = d.iloc[lo:hi]
            if len(g) == 1:
                keep.append(g.iloc[0].to_dict())
                continue
            row = g.iloc[np.argmax(g.target == "13C") if (g.target == "13C").any() else 0].to_dict()
            names = "/".join(sorted(set(g.target)))
            row.update(beta=float(0.5 * (b[lo] + b[hi - 1])), extra=float(max(g.extra.max(), 0.5 * (b[hi - 1] - b[lo]))),
                       rule=("decoy:" if decoy else "") + names, target=names,
                       kind="isotope" if (g.kind == "isotope").all() else g.kind.iloc[0])
            keep.append(row)
    return pd.concat([pd.DataFrame(keep), add], ignore_index=True)


def build_targets(rules: list[Rule], pair_rules: list[PairRule]) -> pd.DataFrame:
    """Every rule as an affine map partner_mz = alpha * anchor_mz + beta (plus an absolute error term)."""
    rows = []
    for r in rules:
        if not r.coelute:
            continue
        for z in range(1, r.z_max + 1):
            for n in range(1, r.n_max + 1):
                if np.gcd(n, z) != 1:
                    continue
                rows.append(dict(rule=r.name, target=r.target or r.name, kind=r.kind, decoy=r.decoy, z=z, n=n,
                                 alpha=1.0, beta=n * r.delta / z, extra=np.sqrt(n) * r.sigma_extra / z,
                                 nominal=n * int(round(r.delta))))
    for p in pair_rules:
        Aa, Ba, na = affine(p.a)
        Ab, Bb, nb_ = affine(p.b)
        alpha = Aa * nb_ / (na * Ab)
        beta = ((Ba / na + p.decoy_shift) * nb_ - Bb) / Ab
        rows.append(dict(rule=p.name, target=p.target or p.name, kind="pair", decoy=p.decoy, z=int(round(Aa)), n=1,
                         alpha=alpha, beta=beta, extra=0.0, nominal=0))
    return merge_families(pd.DataFrame(rows))


# ------------------------------------------------------------------------------------------------ instances
@nb.njit(cache=True)
def _instances(mz, intensity, offsets, scan_of, anchor_ok, alpha, beta, extra, pa, pb, pc, sys_ppm, gate, bridge,
               sig_floor):
    """Best partner of every anchor centroid for every target (same scan first, then +/- bridge scans)."""
    n_t = alpha.shape[0]
    n_scans = offsets.shape[0] - 1
    oi = []
    oj = []
    ot = []
    oe = []
    os_ = []
    for i in range(mz.shape[0]):
        if not anchor_ok[i]:
            continue
        s = scan_of[i]
        si = np.sqrt(pa * pa + (pb * intensity[i] ** (-pc)) ** 2) * 1e-6 * mz[i]
        for t in range(n_t):
            tgt = alpha[t] * mz[i] + beta[t]
            if tgt <= 0:
                continue
            wide = gate * np.sqrt((alpha[t] * si) ** 2 + (sig_floor * 1e-6 * tgt) ** 2 + extra[t] ** 2
                                  + (sys_ppm * 1e-6 * tgt) ** 2)
            best_j = -1
            best_z = 1e18
            best_e = 0.0
            best_s = 0.0
            for k in range(0, 2 * bridge + 1):
                ds = (k + 1) // 2 * (1 if k % 2 == 1 else -1)  # 0, +1, -1, +2, -2, ...
                s2 = s + ds
                if s2 < 0 or s2 >= n_scans:
                    continue
                a, b = offsets[s2], offsets[s2 + 1]
                lo, hi = a, b
                while lo < hi:  # first centroid >= tgt - wide
                    mid = (lo + hi) // 2
                    if mz[mid] < tgt - wide:
                        lo = mid + 1
                    else:
                        hi = mid
                j = lo
                while j < b and mz[j] <= tgt + wide:
                    if j != i:
                        sj = np.sqrt(pa * pa + (pb * intensity[j] ** (-pc)) ** 2) * 1e-6 * mz[j]
                        sig = np.sqrt((alpha[t] * si) ** 2 + sj * sj + extra[t] ** 2 + (sys_ppm * 1e-6 * tgt) ** 2)
                        e = mz[j] - tgt
                        z = abs(e) / sig
                        if z <= gate and z < best_z:
                            best_z = z
                            best_j = j
                            best_e = e
                            best_s = sig
                    j += 1
                if best_j >= 0 and ds == 0:
                    break  # a same-scan partner beats any bridged one
            if best_j >= 0:
                oi.append(i)
                oj.append(best_j)
                ot.append(t)
                oe.append(best_e)
                os_.append(best_s)
    return (np.array(oi, np.int64), np.array(oj, np.int64), np.array(ot, np.int64), np.array(oe), np.array(os_))


# ------------------------------------------------------------------------------------------------ segments
@nb.njit(cache=True)
def _segments(key_order, keys, anchor_scan, partner_scan, anchor_label, partner_label, ai, aj, e, sig, intensity,
              trace_scans_offs, trace_scans, max_gap):
    """Cut instances (sorted by key then anchor scan) into runs and summarise each run.

    Columns: start, stop (instance index range), s0, s1, n, n_anchor, wmean, se, chi2, partner_major,
    partner_frac, ratio_med, ratio_mad, rho, bridged.
    """
    n = key_order.shape[0]
    out = np.zeros((n, 15))
    m = 0
    q = 0
    lab_buf = np.empty(n, np.int64)
    r_buf = np.empty(n)
    while q < n:
        k = keys[key_order[q]]
        r0 = q
        while q < n and keys[key_order[q]] == k:
            q += 1
        r1 = q
        # split this key's instances at anchor-scan gaps > max_gap
        a = r0
        while a < r1:
            b = a + 1
            while b < r1 and anchor_scan[key_order[b]] - anchor_scan[key_order[b - 1]] - 1 <= max_gap:
                b += 1
            sw = 0.0
            swe = 0.0
            nb_ = 0
            nbr = 0
            sa = 0.0
            sp = 0.0
            saa = 0.0
            spp = 0.0
            sap = 0.0
            for u in range(a, b):
                x = key_order[u]
                w = 1.0 / (sig[x] * sig[x])
                sw += w
                swe += w * e[x]
                lab_buf[nb_] = partner_label[x]
                r_buf[nb_] = np.log(intensity[aj[x]] / intensity[ai[x]])
                if partner_scan[x] != anchor_scan[x]:
                    nbr += 1
                ia = np.log(intensity[ai[x]])
                ip = np.log(intensity[aj[x]])
                sa += ia
                sp += ip
                saa += ia * ia
                spp += ip * ip
                sap += ia * ip
                nb_ += 1
            wmean = swe / sw
            chi2 = 0.0
            for u in range(a, b):
                x = key_order[u]
                chi2 += ((e[x] - wmean) / sig[x]) ** 2
            # partner consistency: share of the most common partner trace; "unassigned" (-1) is one class,
            # so a relation whose partner never formed a trace is still consistent
            labs = np.sort(lab_buf[:nb_])
            best = 1
            run = 1
            major = labs[0]
            for v in range(1, nb_):
                if labs[v] == labs[v - 1]:
                    run += 1
                else:
                    run = 1
                if run > best:
                    best = run
                    major = labs[v]
            rs = np.sort(r_buf[:nb_])
            med = rs[nb_ // 2]
            dev = np.abs(rs - med)
            mad = np.sort(dev)[nb_ // 2]
            va = saa - sa * sa / nb_
            vp = spp - sp * sp / nb_
            rho = (sap - sa * sp / nb_) / np.sqrt(va * vp) if va > 1e-12 and vp > 1e-12 else np.nan
            s0 = anchor_scan[key_order[a]]
            s1 = anchor_scan[key_order[b - 1]]
            A = anchor_label[key_order[a]]
            ta, tb = trace_scans_offs[A], trace_scans_offs[A + 1]
            n_anchor = 0
            for v in range(ta, tb):
                if s0 <= trace_scans[v] <= s1:
                    n_anchor += 1
            out[m, 0] = a
            out[m, 1] = b
            out[m, 2] = s0
            out[m, 3] = s1
            out[m, 4] = nb_
            out[m, 5] = n_anchor
            out[m, 6] = wmean
            out[m, 7] = np.sqrt(1.0 / sw)
            out[m, 8] = chi2 / max(nb_ - 1, 1)
            out[m, 9] = major
            out[m, 10] = best / nb_
            out[m, 11] = med
            out[m, 12] = mad
            out[m, 13] = rho
            out[m, 14] = nbr / nb_
            m += 1
            a = b
    return out[:m]


# ------------------------------------------------------------------------------------------------ model
FEATURES = ["log_n", "coverage", "abs_u", "log_chi2", "partner_frac", "ratio_mad", "rho_z", "bridged", "iso_dev",
            "has_iso", "is_isotope", "is_pair", "z2plus", "orphan", "rho_z*log_n", "ratio_mad*log_n", "log_tri",
            "log_tri*log_n"]


def _design(seg: pd.DataFrame) -> np.ndarray:
    log_n = np.log(seg["n"].to_numpy(float))
    rho_z = np.arctanh(np.clip(seg["rho"].fillna(0).to_numpy(), -0.99, 0.99))
    rmad = np.minimum(seg["ratio_mad"].to_numpy(), 3)
    tri = np.log1p(seg["tri"].to_numpy(float)) if "tri" in seg else np.zeros(len(seg))
    X = np.column_stack([
        log_n, seg["coverage"].clip(0, 1), np.minimum(np.abs(seg["u"]), 30),
        np.log(np.maximum(seg["chi2"], 1e-3)), seg["partner_frac"], rmad, rho_z, seg["bridged"],
        np.minimum(np.abs(seg["iso_dev"].fillna(0)), 5), seg["iso_dev"].notna().astype(float),
        (seg["kind"] == "isotope").astype(float), (seg["kind"] == "pair").astype(float),
        (seg["z"] >= 2).astype(float), (seg["partner"] < 0).astype(float), rho_z * log_n, rmad * log_n, tri,
        tri * log_n])
    return X


def _node_ids(B: pd.DataFrame, trace, s0, start: int, pad: int = 3):
    """Block index of (trace, scan); queries outside every block get fresh ids from ``start``."""
    trace = np.asarray(trace)
    s0 = np.asarray(s0)
    idx = _block_of(B, trace, s0)
    ok = (idx >= 0) & (trace >= 0)
    idx_c = np.clip(idx, 0, max(len(B) - 1, 0))
    if len(B):
        ok &= (B["trace"].to_numpy()[idx_c] == trace) & (s0 <= B["s1"].to_numpy()[idx_c] + pad)
    else:
        ok[:] = False
    out = np.where(ok, idx_c, -1)
    miss = np.flatnonzero(out < 0)
    out[miss] = start + np.arange(miss.size)
    return out


def triangle_support(targets: pd.DataFrame, decoys: pd.DataFrame, p_min: float = 0.5):
    """For every edge, the number of third nodes confidently related to both of its ends.

    Real relations sit in dense structures (M, M+1, M+2, adducts and their isotopes); a coincidental edge
    rarely closes a triangle.  The provisional graph uses targets with posterior >= ``p_min``.
    """
    B = _blocks(targets[targets.posterior >= 0.05], 0)
    nxt = len(B)
    ends = {}
    for name, d in (("t", targets), ("d", decoys)):
        a = _node_ids(B, d["anchor"], d["s0"], nxt)
        nxt += len(d)
        b = _node_ids(B, d["partner"], d["s0"], nxt)
        nxt += len(d)
        ends[name] = (a, b)
    ta, tb = ends["t"]
    conf = targets["posterior"].to_numpy() >= p_min
    nbrs: dict = {}
    for x, y in zip(ta[conf].tolist(), tb[conf].tolist()):
        if x != y:
            nbrs.setdefault(x, set()).add(y)
            nbrs.setdefault(y, set()).add(x)
    out = {}
    for name in ("t", "d"):
        a, b = ends[name]
        tri = np.zeros(a.size, np.int64)
        for k, (x, y) in enumerate(zip(a.tolist(), b.tolist())):
            nx, ny = nbrs.get(x), nbrs.get(y)
            if nx and ny:
                tri[k] = len(nx & ny)
        out[name] = tri
    return out["t"], out["d"]


@dataclass
class SegmentModel:
    coef: np.ndarray
    mean: np.ndarray
    scale: np.ndarray
    k: float  # decoy targets per target target

    def prob_target(self, seg):
        X = (_design(seg) - self.mean) / self.scale
        return expit(X @ self.coef[1:] + self.coef[0])

    def posterior(self, seg):
        P = np.clip(self.prob_target(seg), 1e-9, 1 - 1e-9)
        return np.clip((2 * P - 1) / P, 0.0, 1.0)


def fit_segment_model(targets: pd.DataFrame, decoys: pd.DataFrame, k: float, l2: float = 1.0) -> SegmentModel:
    X = np.vstack([_design(targets), _design(decoys)])
    y = np.concatenate([np.ones(len(targets)), np.zeros(len(decoys))])
    w = np.concatenate([np.ones(len(targets)), np.full(len(decoys), 1.0 / k)])
    mu, sd = X.mean(0), X.std(0) + 1e-9
    Xs = (X - mu) / sd

    def f(c):
        z = Xs @ c[1:] + c[0]
        p = expit(z)
        ll = w * (y * np.log(p + 1e-12) + (1 - y) * np.log(1 - p + 1e-12))
        g0 = w * (p - y)
        return -ll.sum() + 0.5 * l2 * (c[1:] ** 2).sum(), np.concatenate([[g0.sum()], Xs.T @ g0 + l2 * c[1:]])

    res = minimize(f, np.zeros(X.shape[1] + 1), jac=True, method="L-BFGS-B")
    return SegmentModel(res.x, mu, sd, k)


# ------------------------------------------------------------------------------------------------ graph
def _blocks(seg: pd.DataFrame, n_traces: int, pad: int = 3):
    """Cut each trace into blocks covering the union of its relation scan ranges."""
    ends = pd.concat([seg[["anchor", "s0", "s1"]].rename(columns={"anchor": "trace"}),
                      seg[seg.partner >= 0][["partner", "s0", "s1"]].rename(columns={"partner": "trace"})])
    ends = ends.sort_values(["trace", "s0"]).to_numpy()
    blocks = []  # trace, s0, s1
    cur = None
    for t, s0, s1 in ends:
        if cur is not None and cur[0] == t and s0 <= cur[2] + pad:
            cur[2] = max(cur[2], s1)
        else:
            if cur is not None:
                blocks.append(cur)
            cur = [t, s0, s1]
    if cur is not None:
        blocks.append(cur)
    B = pd.DataFrame(blocks, columns=["trace", "s0", "s1"])
    return B


def _block_of(B: pd.DataFrame, trace, s0):
    """Block index containing (trace, scan) for arrays of queries."""
    key = B["trace"].to_numpy() * 10_000_000 + B["s0"].to_numpy()
    q = np.asarray(trace) * 10_000_000 + np.asarray(s0)
    return np.searchsorted(key, q, side="right") - 1


# ------------------------------------------------------------------------------------------------ time prior
class _SignalIndex:
    """Prefix sums of every trace's signal in scan order: intensity moments of any (trace, s0..s1) range."""

    def __init__(self, traces):
        run = traces.run
        order, _ = csr_from_labels(traces.label, len(traces))
        sc = run.scan_of_point[order].astype(np.int64)
        self.key = traces.label[order].astype(np.int64) * 10_000_000 + sc
        w = run.intensity[order].astype(float)
        t = run.rt[sc]
        self.C = np.vstack([np.zeros((1, 4)), np.cumsum(np.column_stack([w, w * t, w * t * t, w * w]), 0)])

    def sums(self, trace, s0, s1):
        """Columns sum(w), sum(w t), sum(w t^2), sum(w^2) of trace points in scans [s0, s1] (0 if empty)."""
        trace = np.asarray(trace, np.int64) * 10_000_000
        lo = np.searchsorted(self.key, trace + np.asarray(s0), "left")
        hi = np.searchsorted(self.key, trace + np.asarray(s1), "right")
        return self.C[np.maximum(hi, lo)] - self.C[lo]


def _moments(S, floor):
    """Signal weight, intensity-weighted centre, spread and effective point count from summed moments."""
    W = S[:, 0]
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = S[:, 1] / W
        sd = np.sqrt(np.maximum(S[:, 2] / W - mu * mu, 0.0))
        ne = W * W / S[:, 3]
    return W, mu, np.maximum(sd, floor), ne


@dataclass
class TimePrior:
    """Elution-time prior of a run, learned from its confident relations.

    ``sigma`` is the typical intensity-weighted RT spread of one ion's elution and ``fwhm`` the matching peak
    width: the bulk of every member's signal must lie within ``margin * fwhm`` of the relation that links it.
    ``tau`` is the scatter of member signal centres around their compound's centre beyond each centre's own
    standard error (weak ions near the noise floor need not apex together, so it is learned, not fixed).
    """

    sigma: float
    fwhm: float
    tau: float
    margin: float
    n_edges: int
    inlier: float

    def as_dict(self):
        return dict(sigma=self.sigma, fwhm=self.fwhm, tau=self.tau, margin=self.margin, n_edges=self.n_edges,
                    inlier=self.inlier)


def _fit_tau(d, v0, sigma):
    """Robust scale of centre differences: d ~ pi N(0, 2 tau^2 + v0) + (1 - pi) U(-R, R), by grid + EM."""
    R = max(float(np.abs(d).max()), 10 * sigma) if d.size else 10 * sigma
    best = (-np.inf, 0.25 * sigma, 1.0)
    for tau in sigma * np.logspace(-2, 1, 61):
        var = 2 * tau * tau + v0
        f1 = np.exp(-0.5 * d * d / var) / np.sqrt(2 * np.pi * var)
        f0 = 1.0 / (2 * R)
        pi = 0.8
        for _ in range(50):
            r = pi * f1 / (pi * f1 + (1 - pi) * f0)
            pi = float(np.clip(r.mean(), 1e-3, 1 - 1e-3))
        ll = float(np.log(pi * f1 + (1 - pi) * f0).sum())
        if ll > best[0]:
            best = (ll, float(tau), pi)
    return best[1], best[2]


FWHM_PER_SD = 2 * np.sqrt(2 * np.log(2))


def peak_sigma(nodes, acc, conf, cycle: float = 0.5) -> float:
    """Typical elution spread: median intensity-weighted RT sd of the trace nodes of confident relations."""
    used = np.zeros(len(nodes), bool)
    used[acc["node_a"].to_numpy()[conf]] = True
    used[acc["node_b"].to_numpy()[conf]] = True
    ok = used & (nodes["trace"].to_numpy() >= 0) & (nodes["n_eff"].to_numpy() >= 3)
    return float(np.median(nodes["spread"].to_numpy()[ok])) if ok.sum() >= 20 else 2.0 * cycle


def learn_time_prior(nodes, acc, conf, sigma: float, margin: float = 1.0) -> TimePrior:
    """Centre scatter ``tau`` from confident relations whose relation window holds the bulk of both nodes."""
    tn = nodes["trace"].to_numpy() >= 0
    a, b = acc["node_a"].to_numpy()[conf], acc["node_b"].to_numpy()[conf]
    fa, fb = acc["bulk_a"].to_numpy()[conf], acc["bulk_b"].to_numpy()[conf]
    sel = tn[a] & tn[b] & (fa >= 0.8) & (fb >= 0.8)
    mu, se2 = nodes["center"].to_numpy(), nodes["se"].to_numpy() ** 2
    d = mu[a[sel]] - mu[b[sel]]
    v0 = se2[a[sel]] + se2[b[sel]]
    good = np.isfinite(d) & np.isfinite(v0)
    tau, pi = _fit_tau(d[good], v0[good], sigma) if good.sum() >= 20 else (0.25 * sigma, float("nan"))
    return TimePrior(sigma, float(FWHM_PER_SD * sigma), tau, margin, int(good.sum()), pi)


def _node_moments(nodes, sig, acc, orphan, inst, run, floor):
    """Signal centre / spread / effective points of every graph node (orphans from their partner centroids)."""
    n = len(nodes)
    W = np.zeros(n)
    mu = np.full(n, np.nan)
    sd = np.full(n, np.nan)
    ne = np.zeros(n)
    tn = nodes["trace"].to_numpy() >= 0
    if tn.any():
        S = sig.sums(nodes["trace"].to_numpy()[tn], nodes["s0"].to_numpy()[tn], nodes["s1"].to_numpy()[tn])
        W[tn], mu[tn], sd[tn], ne[tn] = _moments(S, floor)
    if orphan.any():
        o = acc[orphan]
        i0, i1 = o["i0"].to_numpy(), o["i1"].to_numpy()
        L = i1 - i0
        grp = np.repeat(np.arange(len(o)), L)
        pos = np.repeat(i0 - np.concatenate([[0], np.cumsum(L)[:-1]]), L) + np.arange(L.sum())
        pts = inst["partner"][inst["order"][pos]]
        w = run.intensity[pts].astype(float)
        t = run.rt[run.scan_of_point[pts]]
        S = np.column_stack([np.bincount(grp, w * u, len(o)) for u in (1.0, t, t * t, w)])
        idx = o["node_b"].to_numpy()
        W[idx], mu[idx], sd[idx], ne[idx] = _moments(S, floor)
    return W, mu, sd, ne


def _bulk_fraction(sig, nodes, node, s0, s1, rt, margin_s):
    """Share of each node's signal inside the relation window [rt(s0) - margin, rt(s1) + margin]."""
    tr = nodes["trace"].to_numpy()[node]
    out = np.ones(node.size)
    m = tr >= 0
    if not m.any():
        return out
    b0, b1 = nodes["s0"].to_numpy()[node[m]], nodes["s1"].to_numpy()[node[m]]
    w0 = np.searchsorted(rt, rt[s0[m]] - margin_s, "left")
    w1 = np.searchsorted(rt, rt[s1[m]] + margin_s, "right") - 1
    lo, hi = np.maximum(b0, w0), np.minimum(b1, w1)
    inside = np.where(lo <= hi, sig.sums(tr[m], lo, hi)[:, 0], 0.0)
    total = nodes["W"].to_numpy()[node[m]]
    with np.errstate(invalid="ignore", divide="ignore"):
        out[m] = np.where(total > 0, inside / total, 1.0)
    return np.clip(out, 0.0, 1.0)


@dataclass
class ScanGraph:
    segments: pd.DataFrame  # target relation segments with evidence and posterior
    decoys: pd.DataFrame
    blocks: pd.DataFrame  # graph nodes: trace time-blocks (+ orphan partner nodes)
    compounds: pd.DataFrame
    instances: dict
    model: SegmentModel
    fdr: pd.DataFrame
    prior: TimePrior | None = None
    timings: dict = field(default_factory=dict)

    def compound_chromatogram(self, traces, c: int):
        """Summed intensity of a compound's member blocks per scan (its graph-space chromatogram)."""
        run = traces.run
        mem = self.blocks[self.blocks.compound == c]
        s0, s1 = int(mem.s0.min()), int(mem.s1.max())
        y = np.zeros(s1 - s0 + 1)
        for r in mem.itertuples():
            if r.trace < 0:
                continue
            p = traces.points(int(r.trace))
            sc = run.scan_of_point[p]
            m = (sc >= r.s0) & (sc <= r.s1)
            np.add.at(y, sc[m] - s0, run.intensity[p][m])
        return np.arange(s0, s1 + 1), run.rt[s0:s1 + 1], y


def _partition(nodes, acc, run, prior, p_merge, min_posterior, time_prior):
    """Correlation clustering of the nodes into compounds, and the compound table."""
    # merging two blocks needs a confident relation: weights are log-odds relative to p_merge; with the time
    # prior a relation is discounted by how little of each node's signal it explains, and members must
    # co-elute within tau (plus each centre's own error)
    w = acc["log_odds"].to_numpy() - np.log(p_merge / (1 - p_merge))
    a, b = acc["node_a"].to_numpy(), acc["node_b"].to_numpy()
    if time_prior:
        w = w + np.log(np.maximum(acc["bulk_a"].to_numpy(), 1e-6)) + np.log(np.maximum(acc["bulk_b"].to_numpy(), 1e-6))
        prec = 1.0 / (prior.tau ** 2 + nodes["se"].to_numpy() ** 2)
        lab = gaec_coherent(len(nodes), a, b, np.clip(w, -12, 12), nodes["center"].to_numpy(),
                            np.where(nodes["W"].to_numpy() > 0, prec, 0.0))
    else:
        lab = gaec(len(nodes), a, b, np.clip(w, -12, 12))
    good = acc["posterior"].to_numpy() >= max(min_posterior, p_merge)
    la, lb = lab[a[good]], lab[b[good]]
    edges_per = pd.Series(la[la == lb]).value_counts()  # confident relations inside each compound
    nodes = nodes.assign(compound=lab)
    comp = nodes.groupby("compound").agg(n_nodes=("trace", "size"), n_traces=("trace", lambda x: (x >= 0).sum()),
                                         s0=("s0", "min"), s1=("s1", "max")).reset_index()
    span = nodes[nodes.trace >= 0].groupby("compound")["center"].agg(lambda x: float(x.max() - x.min()))
    cen = nodes.assign(wc=nodes.W * nodes.center).groupby("compound")[["wc", "W"]].sum()
    comp["n_edges"] = comp["compound"].map(edges_per).fillna(0).astype(int)
    comp["rt_center"] = comp["compound"].map(cen.wc / cen.W)
    comp["center_span"] = comp["compound"].map(span).fillna(0.0)  # spread of member signal centres (s)
    comp = comp[comp.n_edges >= 1].reset_index(drop=True)
    comp["rt_start"] = run.rt[comp.s0.to_numpy()]
    comp["rt_end"] = run.rt[comp.s1.to_numpy()]
    return lab, comp


def repartition(g: ScanGraph, run, time_prior: bool = True, p_merge: float = 0.9,
                min_posterior: float = 0.5) -> ScanGraph:
    """Re-cut the compounds of a built graph (another merge threshold, or with / without the time prior)."""
    acc = g.segments[g.segments["node_a"] >= 0]
    lab, comp = _partition(g.blocks, acc, run, g.prior, p_merge, min_posterior, time_prior)
    return replace(g, blocks=g.blocks.assign(compound=lab), compounds=comp)


def scan_graph(traces, polarity: int = 1, peptides: bool = False, rules=None, pair_rules=None, gate: float = 4.0,
               sys_ppm: float = 0.3, bridge: int = 1, max_gap: int = 2, min_support: int = 3,
               min_trace_points: int = 5, n_decoys: int = 3, min_posterior: float = 0.5,
               p_merge: float = 0.9, time_prior: bool = True, bulk_margin: float = 1.0) -> ScanGraph:
    """Build, score and partition the scan-level modulo graph of a :class:`~macha.traces.Traces`.

    With ``time_prior`` the partition expects a compound's members to elute together: the bulk of each
    node's signal must lie within ``bulk_margin`` peak widths (learned from the run) of the relation joining
    it, and member signal centres must agree within a learned scatter (see :class:`TimePrior`).
    """
    t0 = time.perf_counter()
    run = traces.run
    pm = traces.precision
    if rules is None or pair_rules is None:
        r0, p0 = default_rules(polarity, peptides)
        rules = rules if rules is not None else r0
        pair_rules = pair_rules if pair_rules is not None else p0
    drules, dprules = make_decoys(rules, pair_rules, n_decoys)
    T = build_targets(rules + drules, pair_rules + dprules)
    label = traces.label
    ok = label >= 0
    counts = np.bincount(label[ok], minlength=len(traces))
    anchor_ok = ok & (counts[np.maximum(label, 0)] >= min_trace_points)
    scan_of = run.scan_of_point.astype(np.int64)
    sig_floor = float(pm.sigma_ppm(np.quantile(run.intensity, 0.01)))
    it = run.intensity.astype(np.float64)
    ai, aj, at, ae, asg = _instances(run.mz, it, run.offsets, scan_of, anchor_ok, T.alpha.to_numpy(),
                                     T.beta.to_numpy(), T.extra.to_numpy(), pm.a, pm.b, pm.c, sys_ppm, gate, bridge,
                                     sig_floor)
    t1 = time.perf_counter()
    a_lab = label[ai]
    p_lab = label[aj]
    # an anchor's own trace is never its partner
    keep = a_lab != p_lab
    ai, aj, at, ae, asg, a_lab, p_lab = ai[keep], aj[keep], at[keep], ae[keep], asg[keep], a_lab[keep], p_lab[keep]
    key = a_lab.astype(np.int64) * T.shape[0] + at
    order = np.lexsort((scan_of[ai], key))
    t_order, t_offs = csr_from_labels(label, len(traces))
    t_scans = scan_of[t_order]
    S = _segments(order, key, scan_of[ai], scan_of[aj], a_lab, p_lab, ai, aj, ae, asg, it, t_offs, t_scans, max_gap)
    cols = ["i0", "i1", "s0", "s1", "n", "n_anchor", "wmean", "se", "chi2", "partner", "partner_frac", "ratio_med",
            "ratio_mad", "rho", "bridged"]
    seg = pd.DataFrame(S, columns=cols)
    for c in ("i0", "i1", "s0", "s1", "n", "n_anchor", "partner"):
        seg[c] = seg[c].astype(np.int64)
    first = order[seg["i0"].to_numpy()]
    seg["anchor"] = a_lab[first]
    tid = at[first]
    for c in ("rule", "target", "kind", "decoy", "z", "n", "nominal"):
        seg[c if c != "n" else "steps"] = T[c].to_numpy()[tid]
    seg["coverage"] = seg["n"] / np.maximum(seg["n_anchor"], 1)
    seg["u"] = seg["wmean"] / seg["se"]
    # 13C ratio against the Poisson carbon expectation (steps 1-2)
    iso = seg["target"].str.contains("13C").to_numpy() & (seg["steps"].to_numpy() <= 2)
    amz = traces.table["mz"].to_numpy()[seg["anchor"].to_numpy()]
    dev = np.full(len(seg), np.nan)
    kind = "peptide" if peptides else "metabolite"
    if iso.any():
        M = amz[iso] * seg["z"].to_numpy()[iso]
        expct = np.array([chem.expected_isotope_ratio(m, int(k), kind) for m, k in zip(M, seg["steps"].to_numpy()[iso])])
        dev[iso] = (seg["ratio_med"].to_numpy()[iso] - np.log(expct)) / 0.35
    seg["iso_dev"] = dev
    seg = seg[seg["n"] >= min_support].reset_index(drop=True)
    t2 = time.perf_counter()
    targets = seg[~seg["decoy"].astype(bool)].reset_index(drop=True)
    decoys = seg[seg["decoy"].astype(bool)].reset_index(drop=True)
    k = max(int(T.decoy.sum()), 1) / max(int((~T.decoy.astype(bool)).sum()), 1)
    # pass 1: local evidence only; pass 2: add triangle support from the provisional graph
    for d in (targets, decoys):
        d["tri"] = 0
    model = fit_segment_model(targets, decoys, k)
    for d in (targets, decoys):
        d["posterior"] = model.posterior(d) if len(d) else []
    if len(targets) and len(decoys):
        targets["tri"], decoys["tri"] = triangle_support(targets, decoys)
        model = fit_segment_model(targets, decoys, k)
    for d in (targets, decoys):
        d["posterior"] = model.posterior(d) if len(d) else []
        d["log_odds"] = np.log(np.clip(d["posterior"], 1e-6, 1 - 1e-6)) - np.log(np.clip(1 - d["posterior"], 1e-6, 1))
    t3 = time.perf_counter()
    # graph: trace blocks + orphan partner nodes
    acc = targets[targets.posterior >= 0.05].copy()
    B = _blocks(acc, len(traces))
    acc["node_a"] = _block_of(B, acc["anchor"], acc["s0"])
    orphan = acc["partner"].to_numpy() < 0
    node_b = np.full(len(acc), -1, np.int64)
    node_b[~orphan] = _block_of(B, acc["partner"].to_numpy()[~orphan], acc["s0"].to_numpy()[~orphan])
    n_blocks = len(B)
    node_b[orphan] = n_blocks + np.arange(orphan.sum())
    acc["node_b"] = node_b
    orph = acc[orphan][["s0", "s1"]].assign(trace=-1)
    nodes = pd.concat([B, orph[["trace", "s0", "s1"]]], ignore_index=True)
    # elution of every node: the bulk of its signal (intensity-weighted centre and spread)
    sig = _SignalIndex(traces)
    floor = float(run.cycle_time) / np.sqrt(12)
    W_, mu_, sd_, ne_ = _node_moments(nodes, sig, acc, orphan, dict(partner=aj, order=order), run, floor)
    nodes["W"], nodes["center"], nodes["spread"], nodes["n_eff"] = W_, mu_, sd_, ne_
    nodes["se"] = nodes["spread"] / np.sqrt(np.maximum(nodes["n_eff"], 1.0))
    # the data's peak width sets the window around a relation that must hold the bulk of both nodes' signal
    conf = acc["posterior"].to_numpy() >= p_merge
    sigma = peak_sigma(nodes, acc, conf, run.cycle_time)
    margin_s = bulk_margin * FWHM_PER_SD * sigma
    s0a, s1a = acc["s0"].to_numpy(), acc["s1"].to_numpy()
    acc["bulk_a"] = _bulk_fraction(sig, nodes, acc["node_a"].to_numpy(), s0a, s1a, run.rt, margin_s)
    acc["bulk_b"] = _bulk_fraction(sig, nodes, acc["node_b"].to_numpy(), s0a, s1a, run.rt, margin_s)
    prior = learn_time_prior(nodes, acc, conf, sigma, bulk_margin)
    lab, comp = _partition(nodes, acc, run, prior, p_merge, min_posterior, time_prior)
    nodes["compound"] = lab
    for c in ("node_a", "node_b"):
        targets[c] = -1
        targets.loc[acc.index, c] = acc[c].to_numpy()
    for c in ("bulk_a", "bulk_b"):
        targets[c] = np.nan
        targets.loc[acc.index, c] = acc[c].to_numpy()
    t4 = time.perf_counter()
    thr = np.linspace(0.05, 0.99, 40)
    pt, pd_ = targets.posterior.to_numpy(), decoys.posterior.to_numpy()
    fdr = pd.DataFrame([dict(threshold=th, accepted=int((pt >= th).sum()),
                             bayes_fdr=float(np.mean(1 - pt[pt >= th])) if (pt >= th).any() else np.nan,
                             decoy_fdr=float((pd_ >= th).sum() / k / max((pt >= th).sum(), 1))) for th in thr])
    inst = dict(anchor=ai, partner=aj, target=at, resid=ae, sigma=asg, order=order)  # segments index into order
    return ScanGraph(targets, decoys, nodes, comp, inst, model, fdr, prior,
                     dict(instances=t1 - t0, segments=t2 - t1, model=t3 - t2, graph=t4 - t3, total=t4 - t0,
                          n_instances=int(ai.size)))
