"""Ground-truth scoring of fragment ladders and fragment grouping on synthetic DIA windows."""

from __future__ import annotations

import numpy as np
import pandas as pd


def _major(vals):
    v, c = np.unique(vals, return_counts=True)
    return v[np.argmax(c)]


def fragment_segment_truth(g, ion_of_point, ions: pd.DataFrame, compounds: pd.DataFrame) -> pd.DataFrame:
    """Per relation segment: same peptide; a ladder relation; and a single-residue step.

    ``ladder`` = both ions in one series (or a/b) at one charge: their difference is then a sum of residues,
    so a matched residue rule is a real sequence relation, possibly an isobaric one (GG = N, GA = Q, S + CO = D).
    ``step`` = consecutive ions with exactly the matched residue.
    """
    seg = g.segments
    anchor, partner, order = g.instances["anchor"], g.instances["partner"], g.instances["order"]
    ion = ions.set_index("ion_id")
    seq = compounds.set_index("compound_id")["sequence"]
    same = np.zeros(len(seg), bool)
    ladder = np.zeros(len(seg), bool)
    step = np.zeros(len(seg), bool)
    for k, r in enumerate(seg.itertuples()):
        idx = order[r.i0: r.i1]
        ia, ip = _major(ion_of_point[anchor[idx]]), _major(ion_of_point[partner[idx]])
        if ia < 0 or ip < 0:
            continue
        a, b = ion.loc[ia], ion.loc[ip]
        if a.compound_id != b.compound_id:
            continue
        same[k] = True
        if r.kind != "residue" or a.frag_z != b.frag_z:
            continue
        ladder[k] = a.series == b.series or {a.series, b.series} == {"a", "b"}
        if a.series != b.series or a.series == "a":
            continue
        if abs(int(a["index"]) - int(b["index"])) != 1:
            continue
        s = seq.loc[a.compound_id]
        n = max(int(a["index"]), int(b["index"]))
        aa = s[n - 1] if a.series == "b" else s[len(s) - n]
        step[k] = (aa in ("L", "I") and r.target == "L/I") or aa == r.target
    return pd.DataFrame(dict(same_peptide=same, ladder=ladder, ladder_step=step), index=seg.index)


def trace_peptide(traces, ion_of_point, ions: pd.DataFrame) -> np.ndarray:
    """Majority peptide of each trace (-1 for noise traces)."""
    lab = traces.label
    ok = (lab >= 0) & (ion_of_point >= 0)
    comp = ions["compound_id"].to_numpy()
    df = pd.DataFrame(dict(t=lab[ok], c=comp[ion_of_point[ok]]))
    maj = df.groupby("t")["c"].agg(lambda x: x.value_counts().index[0])
    out = np.full(len(traces), -1, np.int64)
    out[maj.index.to_numpy()] = maj.to_numpy()
    return out


def trace_fragment(traces, ion_of_point, ions: pd.DataFrame) -> np.ndarray:
    """Majority fragment key (peptide, series, index, charge) of each trace as an integer code (-1 for noise)."""
    key = (ions["compound_id"].astype(np.int64) * 10_000 + ions["series"].map({"b": 0, "y": 1, "a": 2}).astype(np.int64)
           * 1000 + ions["index"].astype(np.int64) * 10 + ions["frag_z"].astype(np.int64)).to_numpy()
    lab = traces.label
    ok = (lab >= 0) & (ion_of_point >= 0)
    df = pd.DataFrame(dict(t=lab[ok], k=key[ion_of_point[ok]]))
    maj = df.groupby("t")["k"].agg(lambda x: x.value_counts().index[0])
    out = np.full(len(traces), -1, np.int64)
    out[maj.index.to_numpy()] = maj.to_numpy()
    return out


def grouping_scores(group_of_trace: np.ndarray, pep: np.ndarray, frag: np.ndarray, compounds: pd.DataFrame,
                    noise: float, min_fragments: int = 3, bins=((0, 10), (10, 30), (30, 100), (100, np.inf))):
    """Peptide recall (a group holding >= min_fragments distinct fragments of it, as its plurality) and purity.

    ``group_of_trace`` is -1 for traces left out of every group.
    """
    ok = (group_of_trace >= 0) & (pep >= 0)
    df = pd.DataFrame(dict(g=group_of_trace[ok], p=pep[ok], f=frag[ok]))
    sizes = df.groupby("g").size()
    multi = df[df.g.isin(sizes[sizes >= 2].index)]
    plural = multi.groupby("g")["p"].agg(lambda x: x.value_counts().index[0])
    top = multi.groupby("g")["p"].agg(lambda x: x.value_counts().iloc[0])
    purity = float(top.sum() / sizes[sizes >= 2].sum()) if len(top) else float("nan")
    m = multi.assign(plural=multi.g.map(plural))
    nf = m[m.p == m.plural].groupby(["g", "p"])["f"].nunique()
    found = set(nf[nf >= min_fragments].index.get_level_values(1))
    snr = compounds.set_index("compound_id")["height"] / noise
    pres = np.unique(pep[pep >= 0])
    rows = []
    for lo, hi in bins:
        sel = [p for p in pres if lo <= snr.loc[p] < hi]
        rows.append(dict(snr_bin=f"{lo:g}-{hi:g}" if np.isfinite(hi) else f">{lo:g}", n=len(sel),
                         recall=float(np.mean([p in found for p in sel])) if sel else np.nan))
    return pd.DataFrame(rows), purity, len(plural)
