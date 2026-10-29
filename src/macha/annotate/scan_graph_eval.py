"""Ground-truth scoring of the scan-level graph against the feature-level graph (synthetic data)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..traces.metrics import adjusted_rand


def _trace_compound(traces, ion_of_point, ions):
    lab = traces.label
    ok = (lab >= 0) & (ion_of_point >= 0)
    comp_of_ion = ions.set_index("ion_id")["compound_id"]
    df = pd.DataFrame(dict(t=lab[ok], c=comp_of_ion.loc[ion_of_point[ok]].to_numpy()))
    maj = df.groupby("t")["c"].agg(lambda x: x.value_counts().index[0])
    out = np.full(len(traces), -1, np.int64)
    out[maj.index.to_numpy()] = maj.to_numpy()
    return out


def segment_truth(g, traces, ion_of_point, ions) -> np.ndarray:
    """True when anchor and partner centroids of a segment share their majority compound (locally, in time)."""
    comp_of_ion = ions.set_index("ion_id")["compound_id"].to_numpy()
    seg = g.segments
    anchor, partner, order = g.instances["anchor"], g.instances["partner"], g.instances["order"]

    def major(pts):
        i = ion_of_point[pts]
        c = np.where(i >= 0, comp_of_ion[np.maximum(i, 0)], -1)
        vals, cnt = np.unique(c, return_counts=True)
        return vals[np.argmax(cnt)]

    out = np.zeros(len(seg), bool)
    for k, r in enumerate(seg.itertuples()):
        idx = order[r.i0: r.i1]
        ca, cp = major(anchor[idx]), major(partner[idx])
        out[k] = ca >= 0 and ca == cp
    return out


def compound_recall(groups_of_ion: dict, ions: pd.DataFrame, compounds: pd.DataFrame, noise: float,
                    bins=((3, 10), (10, 30), (30, 100), (100, np.inf)), group_major: dict | None = None) -> pd.DataFrame:
    """Share of true compounds found as a group joining >= 2 of their ions, by compound SNR.

    ``groups_of_ion`` maps ion id -> set of predicted group ids that contain that ion.  With
    ``group_major`` (group -> its plurality true compound) a compound only counts when it is also the
    plurality compound of that group, so giant merged groups do not inflate recall.
    """
    snr = compounds.set_index("compound_id")["height"] / noise
    multi = ions.groupby("compound_id")["ion_id"].apply(list)
    rows = []
    for c, ion_list in multi.items():
        if len(ion_list) < 2 or compounds.set_index("compound_id").loc[c, "kind"] == "background":
            continue
        seen = {}
        for i in ion_list:
            for g in groups_of_ion.get(i, ()):
                seen[g] = seen.get(g, 0) + 1
        ok = [g for g, v in seen.items() if v >= 2 and (group_major is None or group_major.get(g) == c)]
        rows.append(dict(compound=c, snr=float(snr.loc[c]), found=bool(ok)))
    df = pd.DataFrame(rows)
    out = []
    for lo, hi in bins:
        d = df[(df.snr >= lo) & (df.snr < hi)]
        out.append(dict(snr_bin=f"{lo:g}-{hi:g}" if np.isfinite(hi) else f">{lo:g}", n=len(d),
                        recall=d.found.mean() if len(d) else np.nan))
    return pd.DataFrame(out)


def scan_graph_groups(g, traces, ion_of_point):
    """ion id -> compounds (scan graph) through the member blocks' centroids."""
    run = traces.run
    out = {}
    blocks = g.blocks
    multi = blocks.groupby("compound").size()
    keep = set(multi[multi >= 2].index)
    for r in blocks[blocks.trace >= 0].itertuples():
        if r.compound not in keep:
            continue
        p = traces.points(int(r.trace))
        sc = run.scan_of_point[p]
        ions = ion_of_point[p[(sc >= r.s0) & (sc <= r.s1)]]
        for i in np.unique(ions[ions >= 0]):
            out.setdefault(int(i), set()).add(int(r.compound))
    return out


def feature_graph_groups(ann, truth):
    out = {}
    lab = ann.labels
    sizes = np.bincount(lab)
    for f, i in enumerate(truth["ion"].to_numpy()):
        if i >= 0 and sizes[lab[f]] >= 2:
            out.setdefault(int(i), set()).add(int(lab[f]))
    return out


def block_ari(g, traces, ion_of_point, ions, compounds=None, kind=None):
    """ARI of node compounds against each trace's true compound; ``kind`` "eluting" / "background" restricts it."""
    tc = _trace_compound(traces, ion_of_point, ions)
    b = g.blocks[g.blocks.trace >= 0]
    c = tc[b.trace.to_numpy()]
    ok = c >= 0
    if kind is not None:
        bg = compounds.set_index("compound_id")["kind"].reindex(c).to_numpy() == "background"
        ok &= bg if kind == "background" else ~bg
    return adjusted_rand(b.compound.to_numpy()[ok], c[ok])


def background_mixed(groups_of_ion: dict, ions: pd.DataFrame, compounds: pd.DataFrame) -> float:
    """Share of multi-ion predicted groups joining a long-lived background ion with an eluting compound's ion.

    The typical time-incoherent group: a background ion spans minutes, a coincidental relation during one
    peak's elution attaches that peak to it.
    """
    comp = ions.set_index("ion_id")["compound_id"]
    bg = compounds.set_index("compound_id")["kind"] == "background"
    members = {}
    for i, gs in groups_of_ion.items():
        for g in gs:
            members.setdefault(g, []).append(bool(bg.loc[comp.loc[i]]))
    out = [any(m) and not all(m) for m in members.values() if len(m) >= 2]
    return float(np.mean(out)) if out else float("nan")


def group_majors(groups_of_ion: dict, ions: pd.DataFrame) -> dict:
    """Plurality true compound of every predicted group (by ion count)."""
    comp = ions.set_index("ion_id")["compound_id"]
    tally = {}
    for i, gs in groups_of_ion.items():
        for g in gs:
            tally.setdefault(g, {}).setdefault(int(comp.loc[i]), 0)
            tally[g][int(comp.loc[i])] += 1
    return {g: max(d, key=d.get) for g, d in tally.items()}


def group_purity(groups_of_ion: dict, ions: pd.DataFrame) -> float:
    """Ion-weighted share of each multi-ion group's ions that belong to its plurality compound."""
    comp = ions.set_index("ion_id")["compound_id"]
    tally = {}
    for i, gs in groups_of_ion.items():
        for g in gs:
            tally.setdefault(g, []).append(int(comp.loc[i]))
    tot = top = 0
    for g, cs in tally.items():
        if len(cs) < 2:
            continue
        vals, cnt = np.unique(cs, return_counts=True)
        tot += len(cs)
        top += cnt.max()
    return top / tot if tot else float("nan")
