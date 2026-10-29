"""Graph partitioning and compound inference.

Compound groups are the maximum a posteriori partition of the feature graph
under independent edge log-odds, i.e. correlation clustering.  It is solved
with greedy additive edge contraction (GAEC; Keuper et al. 2015): repeatedly
contract the edge with the largest positive summed log-odds, summing parallel
edges, until none is positive.  Unlike connected components on thresholded
edges, a negative edge between two clusters can veto a merge that a single
weak positive edge would otherwise cause.
"""

from __future__ import annotations

import heapq

import numpy as np
import pandas as pd

from .._kernels import csr_from_labels
from ..chem import NEUTRON_C13
from .rules import adduct_prior, affine


def gaec(n: int, i: np.ndarray, j: np.ndarray, w: np.ndarray) -> np.ndarray:
    adj: list[dict] = [dict() for _ in range(n)]
    for a, b, ww in zip(i.tolist(), j.tolist(), w.tolist()):
        if a == b:
            continue
        v = adj[a].get(b, 0.0) + ww
        adj[a][b] = v
        adj[b][a] = v
    heap = [(-v, a, b) for a in range(n) for b, v in adj[a].items() if a < b and v > 0]
    heapq.heapify(heap)
    parent = np.arange(n)
    alive = np.ones(n, bool)
    while heap:
        negv, a, b = heapq.heappop(heap)
        if not (alive[a] and alive[b]):
            continue
        cur = adj[a].get(b)
        if cur is None or abs(cur + negv) > 1e-9:
            continue
        if cur <= 0:
            break
        if len(adj[a]) < len(adj[b]):
            a, b = b, a
        del adj[a][b]
        del adj[b][a]
        for c, wc in adj[b].items():
            v = adj[a].get(c, 0.0) + wc
            adj[a][c] = v
            adj[c][a] = v
            del adj[c][b]
            if v > 0:
                heapq.heappush(heap, (-v, min(a, c), max(a, c)))
        adj[b] = {}
        alive[b] = False
        parent[b] = a
    while True:
        nxt = parent[parent]
        if np.array_equal(nxt, parent):
            break
        parent = nxt
    _, lab = np.unique(parent, return_inverse=True)
    return lab


def pair_weights(edges: pd.DataFrame, clip: float = 12.0) -> pd.DataFrame:
    """Collapse multiple rule edges between the same two features to their best explanation."""
    a = np.minimum(edges["i"], edges["j"])
    b = np.maximum(edges["i"], edges["j"])
    d = edges.assign(a=a, b=b, w=np.clip(edges["log_odds"], -clip, clip))
    idx = d.groupby(["a", "b"])["w"].idxmax()
    return d.loc[idx].reset_index(drop=True)


def _mono_offset(heights, steps, noise_root, lam, sd=0.35):
    """Choose how many isotopologues (0-2) lie below the lightest observed one.

    Compares observed intensities (relative to the lightest) with a Poisson
    isotope pattern of mean ``lam``; hypotheses that predict a missing lighter
    peak well above the noise are penalised.
    """
    from math import lgamma, log

    def logp(k):
        return k * log(lam) - lgamma(k + 1)

    h0 = heights[steps == 0].max()
    best_o, best_s = 0, np.inf
    for o in (0, 1, 2):
        s = 0.0
        for h, k in zip(heights, steps):
            if k == 0:
                continue
            s += ((np.log(h / h0) - (logp(o + k) - logp(o))) / sd) ** 2
        for k in range(1, o + 1):
            pred = h0 * np.exp(logp(o - k) - logp(o))
            if pred > 5 * noise_root:
                s += (np.log(pred / (5 * noise_root)) / sd) ** 2
        if s < best_s:
            best_o, best_s = o, s
    return best_o


def isotope_envelopes(n: int, edges: pd.DataFrame, labels: np.ndarray, min_posterior: float = 0.5):
    """Greedy selection of non-overlapping isotope envelopes.

    For each charge, 13C steps of 1/z (accepted, inside one compound group) are
    linked when they are each other's best successor / predecessor, giving
    chains.  Chains are accepted in order of total log-odds; a feature joins at
    most one envelope (a chain overlapping an accepted one keeps its longest
    free run).  Other isotope rules then attach fine-structure / heavy-atom
    peaks to envelope members of matching charge.

    Returns (charge, root, nominal index, has_lighter) per feature; charge 0 =
    no isotope evidence.
    """
    acc = edges[(edges["posterior"] >= min_posterior) & (edges["kind"] == "isotope")]
    if len(acc):
        acc = acc[labels[acc["i"].to_numpy()] == labels[acc["j"].to_numpy()]]
    c13 = acc[(acc["target"] == "13C") & (acc["n"] == 1)]
    envs = []
    for z, ez in c13.groupby("z"):
        ez = ez.sort_values("log_odds", ascending=False)
        succ_best = ez.drop_duplicates("i")
        pred_best = ez.drop_duplicates("j")
        mutual = succ_best.merge(pred_best[["i", "j"]], on=["i", "j"])
        succ = dict(zip(mutual["i"].tolist(), mutual["j"].tolist()))
        lo = dict(zip(zip(mutual["i"].tolist(), mutual["j"].tolist()), mutual["log_odds"].tolist()))
        for s0 in set(succ) - set(succ.values()):
            chain, score = [s0], 0.0
            while chain[-1] in succ:
                nxt = succ[chain[-1]]
                score += lo[(chain[-1], nxt)]
                chain.append(nxt)
            envs.append((score, int(z), chain))
    envs.sort(key=lambda t: -t[0])
    taken = np.zeros(n, bool)
    charge = np.zeros(n, np.int64)
    root = np.arange(n)
    nominal = np.zeros(n, np.int64)
    for _, z, chain in envs:
        runs, cur = [], []
        for f in chain:
            if taken[f]:
                if cur:
                    runs.append(cur)
                cur = []
            else:
                cur.append(f)
        if cur:
            runs.append(cur)
        run = max(runs, key=len) if runs else []
        if len(run) < 2:
            continue
        for k, f in enumerate(run):
            taken[f] = True
            charge[f] = z
            root[f] = run[0]
            nominal[f] = k
    other = acc[~((acc["target"] == "13C") & (acc["n"] == 1))].sort_values("log_odds", ascending=False)
    for r in other.itertuples(index=False):
        if taken[r.i] and not taken[r.j] and charge[r.i] == r.z:
            taken[r.j] = True
            charge[r.j] = r.z
            root[r.j] = root[r.i]
            nominal[r.j] = nominal[r.i] + r.nominal
        elif not taken[r.i] and not taken[r.j]:
            taken[r.i] = taken[r.j] = True
            charge[r.i] = charge[r.j] = r.z
            root[r.j] = r.i
            nominal[r.j] = r.nominal
    return charge, root, nominal, root != np.arange(n)


def infer_compounds(features: pd.DataFrame, edges: pd.DataFrame, labels: np.ndarray, polarity: int = 1,
                    ppm: float = 5.0, min_posterior: float = 0.5, peptides: bool = False,
                    default_charge: int | None = None):
    """Charge, isotope roles, adducts and neutral mass for each group.

    Returns (feature_roles DataFrame, groups DataFrame).
    """
    if default_charge is None:
        default_charge = 2 if peptides else 1
    n = len(features)
    mz = features["mz"].to_numpy()
    height = features["height"].to_numpy() if "height" in features else np.ones(n)
    charge, root, total_n, has_lighter = isotope_envelopes(n, edges, labels, min_posterior)
    # peptides: the monoisotopic peak may be below detection (averagine pattern check)
    offset = np.zeros(n, np.int64)
    if peptides:
        from .. import chem

        noise = features["noise"].to_numpy() if "noise" in features else np.full(n, np.median(height) / 100)
        chain_roots = np.unique(root[has_lighter])
        order_r, offs_r = csr_from_labels(np.where(np.isin(root, chain_roots), root, -1), n)
        for r0 in chain_roots:
            mem = order_r[offs_r[r0]: offs_r[r0 + 1]]
            zc = max(int(charge[r0]), 1)
            lam = chem.expected_isotope_ratio(mz[r0] * zc, 1, "peptide")
            offset[mem] = _mono_offset(height[mem], total_n[mem], noise[r0], lam)
    total_n = total_n + offset
    prior = adduct_prior(polarity)
    adducts_by_z = {}
    for name, p in prior.items():
        A, B, nm = affine(name)
        adducts_by_z.setdefault(int(round(A)), []).append((name, A, B, nm, p))
    role_adduct = np.array([""] * n, dtype=object)
    neutral = np.full(n, np.nan)
    rows = []
    order_g, offs_g = csr_from_labels(labels)
    for g in range(offs_g.size - 1):
        members = order_g[offs_g[g]: offs_g[g + 1]]
        if members.size == 0:
            continue
        monos = [m for m in members if not has_lighter[m]]
        zs = {m: (charge[m] if charge[m] > 0 else default_charge) for m in monos}
        hyps = []
        for m in monos:
            mono_mz = mz[m] - offset[m] * NEUTRON_C13 / zs[m]
            for name, A, B, nm, p in adducts_by_z.get(zs[m], []):
                hyps.append((m, name, (mono_mz * A + B) / nm, p))
        bestM, bestscore, bestassign = np.nan, -1.0, {}
        for _, _, Mc, _ in hyps:
            tol = Mc * ppm * 1e-6
            assign = {}
            score = 0.0
            for m, name, M, p in hyps:
                if abs(M - Mc) <= tol and p * np.log1p(height[m]) > assign.get(m, ("", 0.0, 0))[1]:
                    assign[m] = (name, p * np.log1p(height[m]), M)
            score = sum(v[1] for v in assign.values())
            if score > bestscore:
                bestscore, bestassign = score, assign
                bestM = np.average([v[2] for v in assign.values()], weights=[height[k] for k in assign])
        for m in members:
            r = root[m]
            if r in bestassign:
                role_adduct[m] = bestassign[r][0]
                neutral[m] = bestM
        rows.append(dict(group=int(g), n_features=members.size, neutral_mass=bestM,
                         charge=int(max(charge[members].max(), 1)), n_adducts=len(bestassign),
                         adducts=";".join(sorted(v[0] for v in bestassign.values())),
                         rt=float(np.average(features["rt"].to_numpy()[members], weights=height[members])),
                         height=float(height[members].max())))
    roles = pd.DataFrame(dict(group=labels, charge=charge, isotope=total_n, mono=root, adduct=role_adduct,
                              neutral_mass=neutral), index=features.index)
    return roles, pd.DataFrame(rows)
