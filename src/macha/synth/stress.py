"""Stress tests for trace detection: controlled m/z drift and near-isobaric pairs.

Each *case* is one or two ions placed in its own m/z neighbourhood (cases are
3 Da apart) and RT slot, so every case is scored independently.  Case kinds:

* ``linear``    one ion whose m/z drifts at ``mag`` ppm/min
* ``sine``      one ion whose m/z oscillates with amplitude ``mag`` ppm over its peak
* ``step``      one ion whose m/z jumps by ``mag`` ppm at the apex (as at coalescence onset)
* ``charge``    one ion shifted by ``mag`` ppm x I / I_apex (space-charge-like)
* ``coelute``   two ions ``mag`` ppm apart, apex offset 0-2 sigma, intensity ratio 1-10
* ``handoff``   ion B (``mag`` ppm away) starts as ion A ends
* ``cross``     ion A drifts from -``mag`` to +``mag`` ppm across a static ion B while both elute

All centroids then pass through the same Orbitrap-like measurement model as
:func:`macha.synth.simulate_run` (coalescence, precision, threshold, noise).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..run import MSRun
from .instrument import Instrument, coalesce
from .lc import emg, emg_mode

DEFAULT_CASES = {
    "linear": [0, 5, 15, 30, 60],
    "sine": [1, 3, 6],
    "step": [2, 4, 8],
    "charge": [1, 3, 6],
    "coelute": [3, 5, 8, 12, 20],
    "handoff": [3, 6, 12],
    "cross": [5, 10, 20],
}


@dataclass
class StressRun:
    run: MSRun
    ion: np.ndarray  # ion id per centroid (-1 noise)
    merged: np.ndarray
    ions: pd.DataFrame  # ion_id, case, kind, mag, role, mz0, rt, sigma, height
    cases: pd.DataFrame


def make_stress_run(seed: int = 0, cases: dict | None = None, replicates: int = 6, cycle: float = 0.5,
                    instrument: Instrument | None = None, noise_per_scan: float = 80.0) -> StressRun:
    rng = np.random.default_rng(seed)
    inst = instrument or Instrument(resolution_200=240_000, noise=1_000.0, mz_range=(150.0, 1400.0))
    cases = cases or DEFAULT_CASES
    specs = []
    case_rows = []
    cid = 0
    for kind, mags in cases.items():
        for mag in mags:
            for rep in range(replicates):
                case_rows.append(dict(case=cid, kind=kind, mag=mag, rep=rep))
                cid += 1
    n_cases = cid
    slot = 70.0  # seconds of RT per case slot
    per_col = 14  # cases sharing one RT slot are spread over m/z
    run_time = slot * int(np.ceil(n_cases / per_col)) + 60
    rt = np.arange(0.0, run_time, cycle)
    rows = []
    for c in case_rows:
        k = c["case"]
        mz0 = 200.0 + 3.0 * (k % per_col) * 20 + 3.0 * (k // per_col % 20)  # unique neighbourhood per case
        rt0 = 40.0 + slot * (k // per_col) + rng.uniform(-5, 5)
        sigma = rng.uniform(2.0, 4.0)
        tau = sigma * rng.uniform(0.1, 0.6)
        h = 10 ** rng.uniform(4.3, 6.5)  # SNR ~20-3000
        kind, mag = c["kind"], float(c["mag"])
        base = dict(case=k, kind=kind, mag=mag, rt=rt0, sigma=sigma, tau=tau, height=h)
        if kind in ("linear", "sine", "step", "charge"):
            rows.append(dict(base, role="A", mz0=mz0))
        elif kind == "coelute":
            off = rng.uniform(0, 2) * sigma * rng.choice([-1, 1])
            ratio = 10 ** rng.uniform(0, 1)
            rows.append(dict(base, role="A", mz0=mz0))
            rows.append(dict(base, role="B", mz0=mz0 * (1 + mag * 1e-6), rt=rt0 + off, height=h / ratio))
        elif kind == "handoff":
            rows.append(dict(base, role="A", mz0=mz0, tau=0.1 * sigma))
            rows.append(dict(base, role="B", mz0=mz0 * (1 + mag * 1e-6), rt=rt0 + 5.5 * sigma, tau=0.1 * sigma))
        elif kind == "cross":
            rows.append(dict(base, role="A", mz0=mz0))
            rows.append(dict(base, role="B", mz0=mz0, rt=rt0 + rng.uniform(-0.5, 0.5) * sigma,
                             height=h * 10 ** rng.uniform(-0.5, 0.5)))
    ions = pd.DataFrame(rows)
    ions.insert(0, "ion_id", np.arange(len(ions)))
    # a pair is resolvable when its separation exceeds the instrument's coalescence distance
    sep_needed = inst.merge_fraction / inst.resolution(ions.mz0.to_numpy()) * 1e6
    ions["resolvable"] = (ions.kind.isin(["coelute", "handoff", "cross"]).to_numpy() & (ions.mag.to_numpy() > sep_needed)) \
        | ~ions.kind.isin(["coelute", "handoff", "cross"]).to_numpy()

    scans, mzs, its, labs = [], [], [], []
    for r in ions.itertuples(index=False):
        mode = float(emg_mode(r.rt, r.sigma, r.tau))
        lo, hi = np.searchsorted(rt, [r.rt - 5 * r.sigma, r.rt + 5 * r.sigma + 8 * r.tau])
        s = np.arange(lo, hi)
        t = rt[s]
        shape = emg(t, r.rt, r.sigma, r.tau) / emg(mode, r.rt, r.sigma, r.tau)
        if r.kind == "handoff" and r.role == "A":
            shape = np.where(t <= r.rt + 3 * r.sigma, shape, 0.0)  # A stops before B starts
        if r.kind == "handoff" and r.role == "B":
            shape = np.where(t >= r.rt - 2.5 * r.sigma, shape, 0.0)
        inten = r.height * shape
        ppm = np.zeros_like(t)
        if r.kind == "linear":
            ppm = r.mag * (t - mode) / 60.0
        elif r.kind == "sine":
            ppm = r.mag * np.sin(2 * np.pi * (t - mode) / (4 * r.sigma))
        elif r.kind == "step":
            ppm = np.where(t >= mode, r.mag, 0.0)
        elif r.kind == "charge":
            ppm = r.mag * shape
        elif r.kind == "cross" and r.role == "A":
            ppm = np.clip((t - mode) / (2.5 * r.sigma), -1, 1) * r.mag
        keep = inten > 0.2 * inst.threshold
        scans.append(s[keep])
        mzs.append(r.mz0 * (1 + 1e-6 * ppm[keep]))
        its.append(inten[keep])
        labs.append(np.full(keep.sum(), r.ion_id))
    scan = np.concatenate(scans).astype(np.int64)
    mz = np.concatenate(mzs)
    inten = np.concatenate(its)
    lab = np.concatenate(labs).astype(np.int64)
    o = np.lexsort((mz, scan))
    scan, mz, inten, lab = scan[o], mz[o], inten[o], lab[o]
    tol = inst.merge_fraction / inst.resolution(mz)
    scan, mz, inten, lab, merged = coalesce(scan, mz, inten, lab, tol)
    obs_mz, obs_i, ok = inst.measure(rng, scan, mz, inten, rt)
    scan, obs_mz, obs_i, lab, merged = scan[ok], obs_mz[ok], obs_i[ok], lab[ok], merged[ok]
    nn = rng.poisson(noise_per_scan, rt.size)
    ns = np.repeat(np.arange(rt.size), nn)
    nmz = rng.uniform(*inst.mz_range, ns.size)
    ni = inst.threshold * (1 + rng.exponential(0.5, ns.size))
    S = np.concatenate([scan, ns])
    M = np.concatenate([obs_mz, nmz])
    I = np.concatenate([obs_i, ni])
    L = np.concatenate([lab, np.full(ns.size, -1)])
    G = np.concatenate([merged, np.zeros(ns.size, bool)])
    o = np.lexsort((M, S))
    run = MSRun.from_points(S[o], M[o], I[o], rt, 1, meta={"synthetic": True, "scenario": "trace-stress"})
    return StressRun(run, L[o], G[o], ions, pd.DataFrame(case_rows))


def score_stress(sr: StressRun, label: np.ndarray, min_points: int = 3, share: float = 0.1) -> pd.DataFrame:
    """Per-case outcome: split (ion spread over >1 trace), merge (one trace holds two ions)."""
    lab = np.asarray(label)
    ok = lab >= 0
    counts = np.bincount(lab[ok]) if ok.any() else np.zeros(1, int)
    lab = np.where(ok & (counts[np.maximum(lab, 0)] >= min_points), lab, -1)
    df = pd.DataFrame(dict(ion=sr.ion, trace=lab, merged=sr.merged))
    df = df[df.ion >= 0]
    ions = sr.ions.set_index("ion_id")
    out = []
    for case, g in df.groupby(ions.loc[df.ion, "case"].to_numpy()):
        ion_ids = np.unique(g.ion)
        split = False
        completeness = []
        for i in ion_ids:
            gi = g[g.ion == i]
            tc = gi[gi.trace >= 0].trace.value_counts()
            n = len(gi)
            big = tc[tc >= share * n]
            split |= len(big) > 1
            completeness.append(tc.iloc[0] / n if len(tc) else 0.0)
        merge = False
        if len(ion_ids) == 2:
            # a merge: one trace holds >= 20 % of each ion's resolved (non-coalesced) centroids
            res = g[~g.merged & (g.trace >= 0)]
            for t, gt in res.groupby("trace"):
                fr = [np.sum(gt.ion == i) / max(np.sum((g.ion == i) & ~g.merged), 1) for i in ion_ids]
                if min(fr) >= 0.2:
                    merge = True
        row = ions[ions.case == case].iloc[0]
        out.append(dict(case=case, kind=row.kind, mag=row.mag, n_ions=len(ion_ids), split=split, merge=merge,
                        resolvable=bool(ions[ions.case == case].resolvable.all()),
                        completeness=float(np.mean(completeness))))
    return pd.DataFrame(out)
