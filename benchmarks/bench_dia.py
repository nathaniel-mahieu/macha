"""Benchmark + explorer: which way of grouping DIA fragments identifies the most peptides.

Every grouping feeds the same database search (macha.search) and the same peptide-level target-decoy FDR, so the
only difference is how fragments are assembled into pseudo-spectra:

* co-elution only    - traces grouped by elution-profile correlation (DIA-Umpire-style pseudo-spectra)
* ladder graph       - components of the scan-level graph with residue-ladder rules (fragment isotopes, b/y steps)
* ladder-seeded      - ladder components merged when their chromatograms match, then grown by profile correlation
* MS1-anchored       - fragments collected around each MS1 precursor feature by profile correlation (classic)
* MS1-anchored + ladder rescue - MS1-anchored groups plus the ladder-seeded groups whose fragments no MS1
                       precursor claimed (precursors MS1 missed), searched together under one FDR

Window-only groupings get a precursor mass from complementary b / y pairs when the pile-up is supported; the
MS1-anchored grouping uses the MS1 precursor mass.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from common import ASTRAL, NAV, REPORTS, ROOT, Timer, save_json
from macha.annotate.coelution import coelution_edges, coelution_weights
from macha.annotate.complements import complement_mass, complement_null
from macha.annotate.fragment_groups import fragment_groups, ms1_precursors, precursor_groups
from macha.annotate.graph import gaec
from macha.io import DIARun
from macha.pipeline import process_run
from macha.report import SERIES, Page, style
from macha.search import WindowIndex, digest, peptide_fdr, read_fasta, search_groups, with_decoys
from macha.synth import proteomics_scenario
from macha.synth.dia import DIAScenario, make_dia
from macha.traces import detect_traces

APPROACHES = ["co-elution only", "ladder graph", "ladder-seeded", "MS1-anchored", "MS1-anchored + ladder rescue"]
POOLED = {"MS1-anchored + ladder rescue": ("MS1-anchored", "ladder rescue")}
BASE = APPROACHES[:4]
COLORS = dict(zip(APPROACHES, SERIES))
FASTA = [ROOT / "data" / "fasta" / "human.fasta.gz", ROOT / "data" / "fasta" / "yeast.fasta.gz"]
REAL = ASTRAL / "dia" / "sc250pg_A_r1"


def _mass(mz, isolation):
    m, _, n = complement_mass(mz, isolation)
    return m if n >= 3 and n > complement_null(mz, isolation) else np.nan


def window_groups(traces, isolation, ms1=None, min_size: int = 3):
    """Fragment groups (member trace arrays) and precursor masses for every approach in one window."""
    mz = traces.table["mz"].to_numpy()
    out = {}
    with Timer() as t:
        fg = fragment_groups(traces, isolation)
    timing = {"ladder graph": t.s, "ladder-seeded": t.s}
    g = fg.graph
    B = g.blocks[g.blocks.trace >= 0]
    multi = B.groupby("compound")["trace"].nunique()
    lad = [np.unique(d.trace.to_numpy()) for c, d in B[B.compound.isin(multi[multi >= 2].index)].groupby("compound")]
    out["ladder graph"] = lad
    out["ladder-seeded"] = fg.members if fg.members is not None else [
        np.flatnonzero(fg.group_of_trace == k) for k in np.unique(fg.group_of_trace[fg.group_of_trace >= 0])]
    with Timer() as t:
        ce = coelution_edges(traces)
        lab = gaec(len(traces), ce.i.to_numpy(), ce.j.to_numpy(), coelution_weights(ce)) if len(ce) else np.arange(len(traces))
    timing["co-elution only"] = t.s
    order = np.argsort(lab, kind="stable")
    cuts = np.flatnonzero(np.diff(lab[order])) + 1
    out["co-elution only"] = [m for m in np.split(order, cuts)]
    groups, masses = {}, {}
    for k, v in out.items():
        v = [m for m in v if m.size >= min_size]
        groups[k] = [mz[m] for m in v]
        masses[k] = [_mass(mz[m], isolation) for m in v]
    if ms1 is not None:
        ms1_traces, feats = ms1
        with Timer() as t:
            prec = ms1_precursors(feats, isolation)
            mem, _ = precursor_groups(traces, ms1_traces, prec) if len(prec) else ([], None)
        timing["MS1-anchored"] = t.s
        keep = [i for i, m in enumerate(mem) if m.size >= min_size]
        groups["MS1-anchored"] = [mz[mem[i]] for i in keep]
        masses["MS1-anchored"] = [prec.mass.iat[i] for i in keep]
        out["MS1-anchored"] = [mem[i] for i in keep]
        # ladder rescue: ladder-seeded groups mostly made of fragments no MS1 precursor claimed
        claimed = np.zeros(len(traces), bool)
        for m in out["MS1-anchored"]:
            claimed[m] = True
        lad = [(m, x) for m, x in zip([m for m in out["ladder-seeded"] if m.size >= min_size], masses["ladder-seeded"])]
        res = [(m, x) for m, x in lad if claimed[m].mean() < 0.5]
        groups["ladder rescue"] = [mz[m] for m, _ in res]
        masses["ladder rescue"] = [x for _, x in res]
        out["ladder rescue"] = [m for m, _ in res]
        timing["ladder rescue"] = timing["ladder graph"]
    return groups, masses, out, timing


def run_dataset(name, windows, traces_of, isolation_of, ms1, db):
    """Search every approach's groups in every window; return per-approach hits and timing."""
    hits = {a: [] for a in APPROACHES + ["ladder rescue"]}
    sizes = {a: [] for a in APPROACHES + ["ladder rescue"]}
    tim = {a: 0.0 for a in APPROACHES + ["ladder rescue"]}
    for w in windows:
        t0 = time.time()
        tr = traces_of(w)
        iso = isolation_of(w)
        groups, masses, members, timing = window_groups(tr, iso, ms1)
        ix = WindowIndex.build(db, iso)
        for a in groups:
            tim[a] += timing.get(a, 0.0)
            h = search_groups(ix, db, groups[a], masses[a])
            h["window"] = w
            hits[a].append(h)
            sizes[a].extend(len(x) for x in groups[a])
        for p, parts in POOLED.items():
            hits[p].extend(hits[q][-1] for q in parts if q in groups)
            sizes[p].extend(len(x) for q in parts for x in groups.get(q, []))
            tim[p] += timing.get("MS1-anchored", 0.0) + timing.get("ladder graph", 0.0)
        print(f"  {name} window {w}: {len(tr)} traces, groups " +
              ", ".join(f"{a} {len(groups.get(a, []))}" for a in APPROACHES) + f" ({time.time() - t0:.0f}s)", flush=True)
    return {a: pd.concat(v, ignore_index=True) if v else pd.DataFrame() for a, v in hits.items()}, sizes, tim


def summarise(name, hits, sizes, tim, sample_peptides=None):
    rows = []
    for a in APPROACHES:
        h = hits[a]
        if not len(h):
            continue
        ids = peptide_fdr(h, 0.01)
        row = dict(dataset=name, approach=a, groups=len(sizes[a]), median_fragments=float(np.median(sizes[a])),
                   peptides_1pct=len(ids), seconds=tim[a])
        if sample_peptides is not None:
            key = ids.peptide.str.replace("I", "L")
            true = key.isin(sample_peptides)
            row.update(true_ids=int(true.sum()), actual_fdr=float(1 - true.mean()) if len(ids) else np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def synthetic(width, npep=40000, noise=1000.0, seed=3, db_base=None, max_windows=None):
    s = make_dia(DIAScenario(precursors=proteomics_scenario(n_compounds=npep, n_background=0), window_width=width,
                             ms2_noise_points=noise), seed=seed)
    sample = set(s.compounds.sequence.str.replace("I", "L"))
    db = with_decoys(sorted(set(db_base) | set(s.compounds.sequence)))
    ms1 = process_run(s.ms1.run)
    seqs = s.compounds.set_index("compound_id")["sequence"].str.replace("I", "L")
    wins = range(len(s.windows) if max_windows is None else min(max_windows, len(s.windows)))
    hits, sizes, tim = run_dataset(f"synthetic {width:g} Th", wins, lambda w: detect_traces(s.ms2[w].run),
                                   lambda w: (float(s.windows.lo[w]), float(s.windows.hi[w])),
                                   (ms1.traces, ms1.features), db)
    # recall by precursor SNR: sample peptides with fragments in a searched window, found at 1% FDR
    present = set()
    for w in wins:
        present |= set(s.ms2[w].ions.compound_id.unique())
    snr = s.compounds.set_index("compound_id")["height"] / s.scenario.precursors.instrument.noise
    rec = []
    for a in APPROACHES:
        ids = set(peptide_fdr(hits[a], 0.01).peptide.str.replace("I", "L")) if len(hits[a]) else set()
        for lo, hi in ((0, 10), (10, 30), (30, 100), (100, np.inf)):
            sel = [c for c in present if lo <= snr.loc[c] < hi]
            rec.append(dict(dataset=f"synthetic {width:g} Th", approach=a, snr_bin=f"{lo:g}-{hi:g}" if np.isfinite(hi) else f">{lo:g}",
                            n=len(sel), recall=float(np.mean([seqs.loc[c] in ids for c in sel])) if sel else np.nan))
    return summarise(f"synthetic {width:g} Th", hits, sizes, tim, sample), pd.DataFrame(rec), hits


def real(windows, db):
    d = DIARun.load(REAL, windows=list(windows))
    ms1 = process_run(d.ms1)
    hits, sizes, tim = run_dataset("Astral 250 pg", windows, lambda w: detect_traces(d.ms2[w]),
                                   lambda w: (float(d.windows.lo[w]), float(d.windows.hi[w])),
                                   (ms1.traces, ms1.features), db)
    return summarise("Astral 250 pg (HeLa + yeast)", hits, sizes, tim), hits


def ids_figure(S):
    fig = go.Figure()
    for a in APPROACHES:
        d = S[S.approach == a]
        fig.add_trace(go.Bar(x=d.dataset, y=d.peptides_1pct, name=a, marker=dict(color=COLORS[a]),
                             text=d.peptides_1pct, textposition="outside", cliponaxis=False))
    fig.update_layout(barmode="group", bargap=0.3, bargroupgap=0.06)
    fig.update_yaxes(title="peptides identified at 1 % FDR", rangemode="tozero")
    return style(fig, 380, title="Peptides identified, same search for every grouping")


def recall_figure(R, dataset):
    fig = go.Figure()
    for a in APPROACHES:
        d = R[(R.approach == a) & (R.dataset == dataset)]
        fig.add_trace(go.Scatter(x=d.snr_bin, y=d.recall, mode="lines+markers", name=a,
                                 line=dict(width=2, color=COLORS[a]), marker=dict(size=8)))
    fig.update_xaxes(title="true precursor apex SNR (MS1)")
    fig.update_yaxes(title="simulated peptides identified (1 % FDR)", range=[0, 1.02])
    return style(fig, 340, title=f"Identification rate by abundance, {dataset}")


def main(real_windows=(3, 7, 11, 15)):
    t0 = time.time()
    base = digest(read_fasta(*FASTA))
    S, R, H = [], [], {}
    for width in (2.0, 20.0):
        s_, r_, h_ = synthetic(width, db_base=base)
        S.append(s_)
        R.append(r_)
        H[f"synthetic {width:g} Th"] = h_
        print(s_.round(3).to_string(), flush=True)
    s_real, h_real = real(real_windows, with_decoys(base))
    S.append(s_real)
    print(s_real.round(3).to_string(), flush=True)
    S, R = pd.concat(S, ignore_index=True), pd.concat(R, ignore_index=True)
    # overlap of real identifications between approaches
    sets = {a: set(peptide_fdr(h_real[a], 0.01).peptide.str.replace("I", "L")) for a in BASE if len(h_real[a])}
    union = set().union(*sets.values()) if sets else set()
    overlap = pd.DataFrame([dict(approach=a, identified=len(v), only_this=len(v - set().union(*[w for b, w in sets.items() if b != a])),
                                 share_of_union=len(v) / max(len(union), 1)) for a, v in sets.items()])
    rescue = set(peptide_fdr(h_real["MS1-anchored + ladder rescue"], 0.01).peptide.str.replace("I", "L"))
    overlap = pd.concat([overlap, pd.DataFrame([dict(approach="MS1-anchored + ladder rescue", identified=len(rescue),
                                                      only_this=len(rescue - sets.get("MS1-anchored", set())),
                                                      share_of_union=len(rescue) / max(len(union), 1))])],
                        ignore_index=True)
    save_json("dia", dict(summary=S.to_dict("records"), recall=R.to_dict("records"), overlap=overlap.to_dict("records"),
                          real_windows=list(real_windows), union=len(union)))
    write_page(S, R, overlap, len(union), real_windows)
    print(overlap.to_string())
    print("total", round(time.time() - t0), "s")


def write_page(S, R, overlap, n_union, real_windows):
    page = Page("DIA fragment grouping", "Which way of assembling Astral DIA fragments into pseudo-spectra identifies the "
                "most peptides, when every grouping goes through the same search and FDR.", NAV)
    get = lambda ds, a: int(S[(S.dataset.str.startswith(ds)) & (S.approach == a)].peptides_1pct.sum())  # noqa: E731
    page.tiles([
        ("Astral 250 pg: MS1-anchored + ladder rescue", f"{get('Astral', 'MS1-anchored + ladder rescue')}",
         f"MS1-anchored alone {get('Astral', 'MS1-anchored')}, peptides at 1 % FDR"),
        ("Astral 250 pg: ladder-seeded vs co-elution only", f"{get('Astral', 'ladder-seeded')} / "
         f"{get('Astral', 'co-elution only')}", "without any MS1 information"),
        ("Synthetic 20 Th: ladder-seeded vs co-elution only", f"{get('synthetic 20', 'ladder-seeded')} / "
         f"{get('synthetic 20', 'co-elution only')}", "dense 40,000-peptide sample"),
    ])
    page.section("What is compared", """
<p>A DIA isolation window visited once per cycle is an LC-MS run of fragment ions, so each window's MS2 time series
goes through trace detection unchanged. Five groupings then assemble fragment traces into pseudo-spectra:</p>
<ul class="prose">
<li><b>co-elution only</b>: correlation clustering of traces whose apexes coincide, on elution-profile correlation
(DIA-Umpire-style pseudo-spectra without MS1);</li>
<li><b>ladder graph</b>: the <a href="scan_graph.html">scan-level graph</a> with amino-acid residue steps as the modulo
rules (plus fragment 13C, H2O / NH3 / CO losses and 1+/2+ pairs). Its residue decoys sit on the peptide mass-defect
line, more than 0.02 Da from every real step, each keeping one residue's offset from the line, because unrelated
fragment pairs also differ by peptide-like masses;</li>
<li><b>ladder-seeded</b>: ladder components (near-pure, but one per ion series) merged when their summed
chromatograms co-elute and correlate, then grown by every trace that follows a seed's profile (shared);</li>
<li><b>MS1-anchored</b>: the classic precursor-centric route: MS1 features in the window, charge from isotope
spacing, and the fragments whose profiles follow each precursor's MS1 profile (shared between precursors);</li>
<li><b>MS1-anchored + ladder rescue</b>: the MS1-anchored groups plus those ladder-seeded groups whose fragments no
MS1 precursor claimed, i.e. precursors MS1 missed, searched together under one FDR.</li>
</ul>
<p>All groups go through one search (<code>macha.search</code>): tryptic human + yeast peptides (1 missed cleavage,
7-30 residues) and pseudo-reversed decoys whose precursor fits the window at charge 2-4, b / y fragment matching at
10 ppm, and a binomial score against the chance-match rate implied by the group's own fragment density. Groups
without MS1 get a precursor mass from complementary b / y pairs (b<sub>i</sub> + y<sub>L-i</sub> = M + 2H<sup>+</sup>)
when that pile-up beats the same search in wrong windows. Identifications are counted at 1 % peptide-level FDR
(target-decoy competition). On synthetic data the database is the human + yeast digest plus the simulated
peptides, so any accepted human or yeast peptide is a known false identification and the actual FDR is measured.
The search is deliberately simple (no spectral library, intensity or retention-time prediction), so absolute numbers
are far below DIA-NN-class tools; it is the same for every grouping, which is what the comparison needs.</p>""")
    page.figure(ids_figure(S))
    page.table(S, "Per dataset and grouping. Synthetic: 40,000-peptide sample, precursors 500-580 m/z, 2 or 20 Th "
               "windows; true_ids / actual_fdr from the known sample. Real: PXD049412 250 pg HeLa + yeast, Astral "
               "nDIA with 20 Th windows, windows " + ", ".join(str(w) for w in real_windows) + ". Seconds: grouping "
               "time (ladder-seeded includes building the ladder graph).",
               better=dict(peptides_1pct="max", true_ids="max", actual_fdr="min"), group_col="dataset")
    for ds in ("synthetic 2 Th", "synthetic 20 Th"):
        page.figure(recall_figure(R, ds))
    if len(overlap):
        page.table(overlap, "Astral 250 pg: identifications per grouping, those no other grouping found (for the "
                   "rescue arm: those MS1-anchored alone did not find), and the share of the union of the four "
                   "single groupings.", better=dict(identified="max", only_this="max"))
    a = lambda ds: S[S.dataset.str.startswith(ds)].set_index("approach").peptides_1pct  # noqa: E731
    syn20, syn2, rl = a("synthetic 20"), a("synthetic 2 "), a("Astral")
    page.section("What the ladders buy", f"""
<ul class="prose">
<li><b>MS1 is the strongest anchor when the precursor is visible.</b> A precursor mass and a smooth, high-signal
elution template make MS1-anchored groups the most identifiable: {syn2['MS1-anchored']} / {syn20['MS1-anchored']}
peptides on the synthetic 2 / 20 Th runs and {rl['MS1-anchored']} on the Astral windows.</li>
<li><b>Without MS1, ladders beat co-elution where windows are crowded.</b> Ladder-seeded groups identify
{syn20['ladder-seeded']} vs {syn20['co-elution only']} peptides at 20 Th and {rl['ladder-seeded']} vs
{rl['co-elution only']} on the Astral data, where co-isolated precursors make pure co-elution groups chimeric. In
narrow 2 Th windows co-elution is already clean ({syn2['co-elution only']} vs {syn2['ladder-seeded']}).</li>
<li><b>Ladders alone are too sparse.</b> A ladder component only holds fragments whose neighbouring ion is also
observed, and a precursor's b and y series are separate components, so the bare ladder graph identifies little; its
value is as a seed.</li>
<li><b>As a rescue for precursors MS1 misses, ladders add a little.</b> Adding only ladder-seeded groups that
no MS1 precursor claimed gives {rl['MS1-anchored + ladder rescue']} vs {rl['MS1-anchored']} on the Astral windows
({syn20['MS1-anchored + ladder rescue']} vs {syn20['MS1-anchored']} and {syn2['MS1-anchored + ladder rescue']} vs
{syn2['MS1-anchored']} on synthetic). Adding all ladder groups instead lowered identifications in a first run:
every extra group is another chance for a decoy.</li>
</ul>""")
    page.write(REPORTS / "dia.html")


if __name__ == "__main__":
    main(tuple(int(x) for x in sys.argv[1:]) or (3, 7, 11, 15))
