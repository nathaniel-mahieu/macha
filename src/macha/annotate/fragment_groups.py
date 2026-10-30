"""Library-free grouping of DIA fragments into precursors: ladders seed, co-elution grows.

Per isolation window:

1. fragment traces (:func:`macha.traces.detect_traces` on the window's MS2 time series);
2. the scan-level graph with residue-ladder, fragment-isotope, loss and 1+/2+ rules and defect-matched residue
   decoys (:mod:`macha.annotate.ladders`).  Its compounds are near-pure fragment groups, but a precursor's b and
   y ladders are separate components and fragments without an observed neighbour stay out;
3. seeds of one precursor are merged when their summed chromatograms co-elute and correlate
   (:func:`~macha.annotate.coelution.merge_seeds`), and every remaining fragment trace joins the seed whose
   chromatogram it follows best (:func:`~macha.annotate.coelution.assign_to_seeds`);
4. each group gets its ladder evidence (confident residue steps, longest step chain = sequence tag) and a
   precursor mass from complementary b / y pairs (:mod:`macha.annotate.complements`), with a wrong-window null.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .. import chem
from .coelution import assign_to_seeds, merge_seeds
from .complements import complement_mass, complement_null
from .ladders import fragment_decoys, fragment_rules
from .scan_graph import ScanGraph, scan_graph


@dataclass
class FragmentGroups:
    traces: object
    graph: ScanGraph
    group_of_trace: np.ndarray  # -1 = not in a group (best group per trace)
    seed_member: np.ndarray  # trace was in a ladder seed (vs assigned by profile)
    groups: pd.DataFrame
    members: list | None = None  # shared membership: every trace following each group's profile


def _longest_chain(rel: pd.DataFrame, mz: np.ndarray):
    """Longest path of residue steps (anchor -> heavier partner) and its residue string."""
    if not len(rel):
        return 0, ""
    nodes = np.unique(np.concatenate([rel.anchor.to_numpy(), rel.partner.to_numpy()]))
    best = {int(n): (0, "") for n in nodes}
    for r in rel.assign(m=mz[rel.anchor.to_numpy()]).sort_values("m").itertuples():
        a, b = int(r.anchor), int(r.partner)
        cand = (best[a][0] + 1, best[a][1] + (r.target if len(r.target) == 1 else "L"))
        if cand[0] > best[b][0]:
            best[b] = cand
    return max(best.values())


def fragment_groups(traces, isolation: tuple, merge_r0: float = 0.9, assign_r: float = 0.8,
                    p_conf: float = 0.9, graph=None, shared: bool = True, **graph_kw) -> FragmentGroups:
    """Group the fragment traces of one isolation window into precursors (see module docstring)."""
    if graph is None:
        rules, pairs = fragment_rules()
        graph = scan_graph(traces, peptides=True, rules=rules, pair_rules=pairs,
                           decoys=fragment_decoys(rules, pairs), **graph_kw)
    g = graph
    B = g.blocks[g.blocks.trace >= 0]
    multi = B.groupby("compound").size()
    members = B[B.compound.isin(multi[multi >= 2].index)].rename(columns={"compound": "seed"})[
        ["seed", "trace", "s0", "s1"]]
    seed_member = np.zeros(len(traces), bool)
    seed_member[members.trace.to_numpy()] = True
    if len(members):
        members = merge_seeds(traces, members, r0=merge_r0)
    grp, _ = assign_to_seeds(traces, members, r_min=assign_r)
    shared_members = assign_to_seeds(traces, members, r_min=assign_r, shared=True) if shared and len(members) else None
    table = traces.table
    mz = table["mz"].to_numpy()
    seg = g.segments
    conf = seg[(seg.posterior >= p_conf) & (seg.partner >= 0)]
    iso_partner = np.zeros(len(traces), bool)
    iso_partner[conf.partner[conf.kind == "isotope"].to_numpy()] = True
    res = conf[conf.kind == "residue"]
    res = res[(grp[res.anchor.to_numpy()] >= 0) & (grp[res.anchor.to_numpy()] == grp[res.partner.to_numpy()])]
    res_by = {k: d for k, d in res.groupby(grp[res.anchor.to_numpy()])}
    rows = []
    for k in np.unique(grp[grp >= 0]):
        mem = np.flatnonzero(grp == k)
        frag = mem[~iso_partner[mem]]
        mass, z, n_pairs = complement_mass(mz[frag], isolation)
        rel = res_by.get(k, res.iloc[:0])
        chain, tag = _longest_chain(rel, mz)
        apex_t = table["apex_rt"].to_numpy()[mem]
        w = table["apex_intensity"].to_numpy()[mem]
        rows.append(dict(group=int(k), n_traces=mem.size, n_seed=int(seed_member[mem].sum()), n_fragments=frag.size,
                         rt=float(np.average(apex_t, weights=w)), residue_steps=len(rel), longest_chain=chain,
                         tag=tag, mass=mass, charge=z, complement_pairs=n_pairs,
                         complement_null=complement_null(mz[frag], isolation)))
    return FragmentGroups(traces, g, grp, seed_member, pd.DataFrame(rows), shared_members)


def ms1_precursors(features: pd.DataFrame, isolation: tuple, charges=(2, 3, 4), ppm: float = 10.0) -> pd.DataFrame:
    """MS1 features that can be the monoisotopic precursor of a window, with charge from isotope spacing.

    A feature whose +1/z isotope co-elutes gets that charge; one that is itself such an isotope of a stronger
    co-eluting feature is dropped; with no isotope seen, charges 2 and 3 are both kept as hypotheses.
    """
    lo, hi = isolation
    f = features.sort_values("mz").reset_index(drop=True)
    mz, rt = f["mz"].to_numpy(), f["rt"].to_numpy()
    h = f["height"].to_numpy()
    tol_rt = 0.5 * np.maximum(f["fwhm_s"].to_numpy(), 1.0)

    def partner(i, dz):
        t = mz[i] + dz
        j = np.flatnonzero(np.abs(mz - t) <= ppm * 1e-6 * t)
        j = j[np.abs(rt[j] - rt[i]) <= tol_rt[i]]
        return j

    rows = []
    cand = np.flatnonzero((mz >= lo - 2 * chem.NEUTRON_C13 / 2) & (mz < hi))
    for i in cand:
        zs = [z for z in charges if partner(i, chem.NEUTRON_C13 / z).size]
        if any((h[partner(i, -chem.NEUTRON_C13 / z)] > h[i]).any() for z in charges):
            continue  # an isotope of a stronger precursor
        for z in (zs[:1] if zs else (2, 3)):
            if mz[i] >= lo - 2 * chem.NEUTRON_C13 / z:
                rows.append(dict(feature=int(i), mz=float(mz[i]), z=int(z), mass=float((mz[i] - chem.PROTON) * z),
                                 rt=float(rt[i]), charge_seen=bool(zs), **{c: f[c].iat[i] for c in
                                                                           ("trace", "scan_start", "scan_end")}))
    return pd.DataFrame(rows)


def precursor_groups(traces, ms1_traces, precursors: pd.DataFrame, r_min: float = 0.8) -> tuple[list, pd.DataFrame]:
    """Classic precursor-centric grouping: fragments whose profiles follow each MS1 precursor's (shared allowed).

    The precursor's MS1 elution profile is interpolated onto the window's MS2 scan times.  Returns the member
    traces of every precursor and the (trace, seed, r) assignment table.
    """
    run2, run1 = traces.run, ms1_traces.run
    seed_s0, seed_y = [], []
    for p in precursors.itertuples():
        pts = ms1_traces.points(int(p.trace))
        sc = run1.scan_of_point[pts]
        m = (sc >= p.scan_start) & (sc <= p.scan_end)
        t1, y1 = run1.rt[sc[m]], run1.intensity[pts][m].astype(float)
        a = int(np.searchsorted(run2.rt, t1[0])) if t1.size else 0
        b = int(np.searchsorted(run2.rt, t1[-1], "right")) if t1.size else 0
        seed_s0.append(a)
        seed_y.append(np.interp(run2.rt[a:b], t1, y1) if b > a else np.zeros(0))
    df = assign_to_profiles_shared(traces, np.array(seed_s0, np.int64), seed_y, r_min)
    members = [df.trace[df.seed == k].to_numpy() for k in range(len(precursors))]
    return members, df


def assign_to_profiles_shared(traces, seed_s0, seed_y, r_min):
    from .coelution import assign_to_profiles

    keep = [len(y) >= 3 for y in seed_y]
    idx = np.flatnonzero(keep)
    df = assign_to_profiles(traces, seed_s0[idx], [seed_y[i] for i in idx], r_min, shared=True)
    return df.assign(seed=idx[df.seed.to_numpy()]) if len(df) else df
