"""Benchmark + explorer: the scan-level modulo graph (compounds detected before peak detection)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from common import NAV, REAL_LABEL, REPORTS, Timer, real_run, save_json
from macha.annotate import ProfileSource, annotate
from macha.annotate.scan_graph import repartition, scan_graph
from macha.annotate.scan_graph_eval import (background_mixed, block_ari, compound_recall, feature_graph_groups,
                                            group_majors, group_purity, scan_graph_groups, segment_truth)
from macha.pipeline import process_run
from macha.report import LIGHT, SERIES, Page, method_color, style
from macha.synth import make_run, metabolomics_scenario, proteomics_scenario
from macha.synth.truth import feature_truth

SCENARIOS = {
    "metabolomics, low abundance": (lambda: metabolomics_scenario(log10_height=(4.3, 0.8)), False,
                                    "600 metabolites, median apex SNR ~20"),
    "metabolomics": (lambda: metabolomics_scenario(), False, "600 metabolites, median apex SNR ~100"),
    "proteomics": (lambda: proteomics_scenario(), True, "1200 peptides, z 1-5, 240k"),
}
LABELS = {"scan": "scan-level modulo graph", "scan_noprior": "scan-level graph, no elution-time prior",
          "feature": "feature-level graph (after peak detection)"}


def span_share(g, k=3.0):
    """Share of multi-trace compounds whose member signal centres span more than k learned peak widths."""
    c = g.compounds[g.compounds.n_traces >= 2]
    return float((c.center_span > k * g.prior.fwhm).mean()) if len(c) else float("nan")


def evaluate(name, seed=2):
    make, pep, _ = SCENARIOS[name]
    sr = make_run(make(), seed=seed)
    s = process_run(sr.run)
    with Timer() as t:
        g = scan_graph(s.traces, peptides=pep)
    g0 = repartition(g, s.traces.run, time_prior=False)
    seg = g.segments.assign(true=segment_truth(g, s.traces, sr.ion, sr.ions))
    ft = feature_truth(s.features, s.traces, sr.ion, sr.ions, sr.compounds)
    with Timer() as t2:
        ann = annotate(s.features, ProfileSource(s.baselines), peptides=pep, cycle=sr.run.cycle_time)
    noise = sr.scenario.instrument.noise
    gs, gf = scan_graph_groups(g, s.traces, sr.ion), feature_graph_groups(ann, ft)
    gs0 = scan_graph_groups(g0, s.traces, sr.ion)
    rec = []
    for m, groups in (("scan", gs), ("scan_noprior", gs0), ("feature", gf)):
        r = compound_recall(groups, sr.ions, sr.compounds, noise, group_major=group_majors(groups, sr.ions))
        rec.append(r.assign(method=m, scenario=name))
    fe = ann.edges
    comp = ft["compound"].to_numpy()
    ftrue = (comp[fe.i] == comp[fe.j]) & (comp[fe.i] >= 0)
    prec9 = seg.true[seg.posterior >= 0.9].mean()
    summ = [dict(scenario=name, method=m, relation_precision_p09=prec9, relations_p09=int((seg.posterior >= 0.9).sum()),
                 group_purity=group_purity(groups, sr.ions), background_mixed=background_mixed(groups, sr.ions, sr.compounds),
                 node_ari_eluting=block_ari(x, s.traces, sr.ion, sr.ions, sr.compounds, "eluting"),
                 node_ari_background=block_ari(x, s.traces, sr.ion, sr.ions, sr.compounds, "background"),
                 span_gt_3_widths=span_share(x), seconds=t.s)
            for m, x, groups in (("scan", g, gs), ("scan_noprior", g0, gs0))]
    summ.append(dict(scenario=name, method="feature", relation_precision_p09=ftrue[fe.posterior >= 0.9].mean(),
                     relations_p09=int((fe.posterior >= 0.9).sum()), group_purity=group_purity(gf, sr.ions),
                     background_mixed=background_mixed(gf, sr.ions, sr.compounds), node_ari_eluting=np.nan,
                     node_ari_background=np.nan, span_gt_3_widths=np.nan,
                     seconds=t2.s + s.timings["peaks"] + s.timings["baseline"]))
    prior = dict(dataset=f"synthetic: {name}", **g.prior.as_dict())
    return pd.concat(rec), pd.DataFrame(summ), seg, prior


def recall_figure(R):
    fig = go.Figure()
    for k, sc in enumerate(SCENARIOS):
        for m, dash_slot in (("scan", 0), ("feature", 1)):
            d = R[(R.scenario == sc) & (R.method == m)]
            fig.add_trace(go.Scatter(x=d.snr_bin, y=d.recall, mode="lines+markers",
                                     name=f"{sc}: {LABELS[m]}", line=dict(width=2 if m == "scan" else 1,
                                                                          color=SERIES[k]),
                                     marker=dict(size=9 if m == "scan" else 6,
                                                 symbol="circle" if m == "scan" else "square")))
    fig.update_xaxes(title="true compound apex SNR")
    fig.update_yaxes(title="compounds found (≥ 2 ions, plurality of their group)", range=[0, 1.03])
    return style(fig, 400, title="Compound recall: thick lines with circles = scan graph, thin with squares = "
                                 "feature graph")


def calibration_figure(segs):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", name="perfect calibration",
                             line=dict(width=1, color=LIGHT["axis"])))
    for k, (name, seg) in enumerate(segs.items()):
        b = np.clip((seg.posterior * 10).astype(int), 0, 9)
        c = seg.groupby(b).agg(p=("posterior", "mean"), f=("true", "mean"), n=("true", "size"))
        c = c[c.n >= 20]
        fig.add_trace(go.Scatter(x=c.p, y=c.f, mode="lines+markers", name=name, line=dict(width=2, color=SERIES[k]),
                                 marker=dict(size=8)))
    fig.update_xaxes(title="relation posterior (binned)", range=[0, 1])
    fig.update_yaxes(title="fraction true (synthetic truth)", range=[0, 1.02])
    return style(fig, 360, title="Calibration of relation-segment posteriors")


def span_figure(rows):
    """Share of compounds spanning > 3 peak widths, with and without the prior, per dataset."""
    d = pd.DataFrame(rows)
    fig = go.Figure()
    for col, lab, color in (("without", "no elution-time prior", SERIES[1]),
                            ("with", "elution-time prior", SERIES[0])):
        fig.add_trace(go.Bar(x=d.dataset, y=d[col], name=lab, marker=dict(color=color),
                             text=[f"{v:.1%}" for v in d[col]], textposition="outside", cliponaxis=False))
    fig.update_layout(barmode="group", bargap=0.35, bargroupgap=0.08)
    fig.update_yaxes(title="multi-ion compounds whose member centres<br>span > 3 peak widths", tickformat=".0%",
                     rangemode="tozero")
    return style(fig, 360, title="Time-incoherent compounds, by dataset")


def split_view(g0, g, traces, c0, title):
    """Members of one compound found without the prior, coloured by the compound each joins with the prior."""
    run = traces.run
    mem = g0.blocks[(g0.blocks.compound == c0) & (g0.blocks.trace >= 0)]
    new = g.blocks.compound.to_numpy()[mem.index.to_numpy()]
    slot = {c: k for k, c in enumerate(pd.unique(new[np.argsort(mem.center.to_numpy())]))}
    fig = go.Figure()
    shown = set()
    for r, c in zip(mem.sort_values("center").itertuples(), new[np.argsort(mem.center.to_numpy())]):
        p = traces.points(int(r.trace))
        sc = run.scan_of_point[p]
        m = (sc >= r.s0) & (sc <= r.s1)
        it = run.intensity[p][m]
        k = slot[c]
        fig.add_trace(go.Scatter(x=run.rt[sc[m]], y=it / it.max(), mode="lines+markers",
                                 name=f"compound {k + 1} with the prior", legendgroup=str(k),
                                 showlegend=k not in shown, line=dict(width=1.5, color=SERIES[k % 8]),
                                 marker=dict(size=4),
                                 hovertext=f"m/z {traces.table['mz'].iat[int(r.trace)]:.4f}"))
        shown.add(k)
    fig.update_xaxes(title="retention time (s)")
    fig.update_yaxes(title="intensity (each block normalised to its max)")
    return style(fig, 340, title=title)


def compound_view(g, traces, c, title):
    """Member blocks of one graph-space compound: normalised intensities over time, plus the summed chromatogram."""
    run = traces.run
    mem = g.blocks[(g.blocks.compound == c) & (g.blocks.trace >= 0)]
    fig = go.Figure()
    scans, rt, y = g.compound_chromatogram(traces, c)
    fig.add_trace(go.Scatter(x=rt, y=y / max(y.max(), 1), mode="lines", name="compound chromatogram (sum)",
                             line=dict(width=3, color=LIGHT["ink"])))
    for k, r in enumerate(mem.sort_values("trace").head(8).itertuples()):
        p = traces.points(int(r.trace))
        sc = run.scan_of_point[p]
        m = (sc >= r.s0) & (sc <= r.s1)
        it = run.intensity[p][m]
        fig.add_trace(go.Scatter(x=run.rt[sc[m]], y=it / it.max(), mode="lines+markers",
                                 name=f"m/z {traces.table['mz'].iat[int(r.trace)]:.4f}",
                                 line=dict(width=1.5, color=SERIES[k % 8]), marker=dict(size=5)))
    fig.update_xaxes(title="retention time (s)")
    fig.update_yaxes(title="intensity (each normalised to its max)")
    return style(fig, 320, title=title)


def main():
    R, S, segs, priors, spans = [], [], {}, [], []
    for name in SCENARIOS:
        r, s_, seg, pr_ = evaluate(name)
        R.append(r)
        S.append(s_)
        segs[name] = seg
        priors.append(pr_)
        q = s_.set_index("method").span_gt_3_widths
        spans.append({"dataset": f"synthetic: {name}", "without": q["scan_noprior"], "with": q["scan"]})
        print(name, "done")
    R, S = pd.concat(R), pd.concat(S)
    real = []
    examples = []
    split = None
    for ds, pep in (("metabolomics", False), ("proteomics_sc", True), ("proteomics", True)):
        run = real_run(ds, 0)
        if run is None:
            continue
        smp = process_run(run)
        with Timer() as t:
            g = scan_graph(smp.traces, peptides=pep)
        g0 = repartition(g, run, time_prior=False)
        priors.append(dict(dataset=REAL_LABEL[ds], **g.prior.as_dict()))
        spans.append({"dataset": REAL_LABEL[ds], "without": span_share(g0), "with": span_share(g)})
        f = g.fdr
        f9 = f.iloc[(f.threshold - 0.9).abs().argmin()]
        seg = g.segments
        acc = seg.posterior >= 0.9
        peaks = {}
        for tt, a, b in zip(smp.features.trace, smp.features.scan_start, smp.features.scan_end):
            peaks.setdefault(int(tt), []).append((a, b))
        nopeak = np.array([not any(a <= s1 and b >= s0 for a, b in peaks.get(int(tt), []))
                           for tt, s0, s1 in zip(seg.anchor[acc], seg.s0[acc], seg.s1[acc])])
        real.append(dict(dataset=REAL_LABEL[ds], centroids=run.n_points, relation_segments=len(seg),
                         decoy_segments=len(g.decoys), relations_p09=int(acc.sum()), bayes_fdr_09=f9.bayes_fdr,
                         decoy_fdr_09=f9.decoy_fdr, graph_compounds=len(g.compounds),
                         graph_compounds_no_prior=len(g0.compounds), span_gt_3_widths=span_share(g),
                         span_gt_3_widths_no_prior=span_share(g0),
                         relations_without_picked_peak=float(nopeak.mean()), seconds=t.s))
        if ds == "metabolomics":
            real_fdr = (ds, f)
            comp = g.compounds.sort_values("n_edges", ascending=False)
            # one well-supported compound and one with no picked peak on any member
            pick = [comp.compound.iat[min(5, len(comp) - 1)]]
            for c in comp.compound.head(400):
                mem = g.blocks[(g.blocks.compound == c) & (g.blocks.trace >= 0)]
                if len(mem) >= 3 and not any(any(a <= r.s1 and b >= r.s0 for a, b in peaks.get(int(r.trace), []))
                                             for r in mem.itertuples()):
                    pick.append(c)
                    break
            examples = [(g, smp.traces, c, lab) for c, lab in zip(pick, ["well-supported compound",
                                                                         "compound with no picked peak on any member"])]
            # the most time-incoherent small compound without the prior, and how the prior splits it
            traced = g0.blocks.trace >= 0

            def n_split(c):  # compounds its traced members fall into with the prior
                return g.blocks.compound[traced & (g0.blocks.compound == c)].nunique()

            c0 = g0.compounds[(g0.compounds.n_traces >= 3) & (g0.compounds.n_traces <= 8)]
            c0 = c0[c0.compound.map(n_split) >= 2]
            if len(c0):
                top = c0.sort_values("center_span").iloc[-1]
                split = (g0, g, smp.traces, int(top.compound),
                         f"{REAL_LABEL[ds]}: a compound whose {int(top.n_traces)} traced members span "
                         f"{top.center_span:.0f} s without the prior, split into {n_split(top.compound)} by it")
        print(ds, "done")
    real = pd.DataFrame(real)
    priors, spans = pd.DataFrame(priors), pd.DataFrame(spans)
    save_json("scan_graph", dict(recall=R.to_dict("records"), summary=S.to_dict("records"),
                                 real=real.to_dict("records"), priors=priors.to_dict("records"),
                                 spans=spans.to_dict("records")))
    page = Page("Scan-level modulo graph",
                "Compounds detected in graph space before any peak picking: mass-difference rules tested between "
                "centroids scan by scan, relations aggregated over time, decoy-calibrated, and partitioned into "
                "compounds.", NAV)
    sc = R[R.scenario == "metabolomics, low abundance"].set_index(["method", "snr_bin"]).recall
    pr = R[R.scenario == "proteomics"].set_index(["method", "snr_bin"]).recall
    page.tiles([
        ("compounds found, SNR 10-30 (metabolomics)", f"{sc[('scan', '10-30')]:.0%}",
         f"feature graph {sc[('feature', '10-30')]:.0%}"),
        ("compounds found, SNR 3-10 (proteomics)", f"{pr[('scan', '3-10')]:.0%}",
         f"feature graph {pr[('feature', '3-10')]:.0%}"),
        ("Astral: Bayesian vs decoy FDR @0.9", f"{100 * real.bayes_fdr_09.iat[0]:.1f}% / "
                                                f"{100 * real.decoy_fdr_09.iat[0]:.1f}%" if len(real) else "–",
         real.dataset.iat[0] if len(real) else ""),
        ("Astral compounds spanning > 3 peak widths", f"{real.span_gt_3_widths.iat[0]:.1%}" if len(real) else "–",
         f"without the elution-time prior {real.span_gt_3_widths_no_prior.iat[0]:.1%}" if len(real) else ""),
    ])
    page.section("Feature level vs scan level", """
<p>The feature-level graph (<a href="annotate.html">annotation</a>) relates picked peaks, so it inherits every miss of
peak detection and judges co-elution from apex times. Here the modulo test runs between <b>centroids in every
scan</b>: each centroid of a trace is an anchor, the best partner for every rule target (n·Δ/z, or a neutral-mass
map for charge states and multimers) is searched in the same scan and, when absent, in the neighbouring scans
(<b>±1 scan bridging</b> for partners near the noise floor). Partners need not belong to a trace, so isotopes
too weak to form one are rescued. Overlapping ~2 Da isotope targets are merged into one family so every M+2 pair is
tested once.</p>
<p><b>Relations over time.</b> A real relation persists across the anchor's elution: one consistent partner, a
residual whose inverse-variance mean shrinks as 1/√n, a constant intensity ratio, co-varying intensities. Instances
are cut into segments of consecutive supporting scans and summarised by those statistics. The same search is run
with decoy rules; a weighted logistic model of targets vs decoys gives each segment a posterior (target-decoy local
FDR). A second pass adds <b>triangle support</b>: the number of third nodes confidently related to both ends of an
edge, since real relations form dense structures (M, M+1, M+2, adducts and their isotopes) and coincidences
rarely close triangles.</p>
<p><b>Compounds in graph space.</b> Traces are cut into blocks covering the time ranges of their relations (isomers
sharing an m/z stay apart), blocks become nodes, segments become edges weighted by log-odds relative to a merge
threshold, and correlation clustering yields compounds with a time extent and a summed compound chromatogram.
Peak detection, if wanted, then runs on that chromatogram, where every member ion adds signal. An
<a href="#elution-time-prior">elution-time prior</a> keeps each compound's members within a peak width of each
other.</p>""")
    page.section("Accuracy against ground truth", "<p>A compound counts as found when one group joins at least two "
                 "of its ions and that compound is the group's plurality, so oversized merged groups do not count. "
                 "Relation precision is the fraction of relations at posterior ≥ 0.9 linking two ions of one "
                 "compound.</p>")
    page.figure(recall_figure(R))
    page.table(S.assign(method=S.method.map(LABELS)), "Per scenario: " + "; ".join(
        f"<b>{k}</b>: {v[2]}" for k, v in SCENARIOS.items()),
               better=dict(relation_precision_p09="max", group_purity="max", background_mixed="min",
                           node_ari_eluting="max", span_gt_3_widths="min"), group_col="scenario")
    page.figure(calibration_figure(segs))
    if len(real):
        page.section("Real Orbitrap Astral runs", "<p>No truth: the decoy FDR checks calibration. "
                     "<i>relations_without_picked_peak</i> is the share of confident relations whose anchor has no "
                     "picked peak overlapping the relation's scans, i.e. evidence that exists only in graph space "
                     "(some of it from chemical background, whose ions also carry real isotopes).</p>")
        page.table(real, better=dict(seconds="min"))
        ds, f = real_fdr
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=f.threshold, y=f.bayes_fdr, mode="lines", name="Bayesian FDR",
                                 line=dict(width=2, color=method_color("macha"))))
        fig.add_trace(go.Scatter(x=f.threshold, y=f.decoy_fdr, mode="lines", name="decoy-estimated FDR",
                                 line=dict(width=2, color=method_color("greedy"))))
        fig.update_xaxes(title="posterior threshold")
        fig.update_yaxes(title="false discovery rate", rangemode="tozero")
        page.figure(style(fig, 320, title=f"Relation-segment FDR, {REAL_LABEL[ds]}"))
        for g, traces, c, lab in examples:
            page.figure(compound_view(g, traces, c, f"{REAL_LABEL['metabolomics']}: {lab} (graph compound {c})"))
    page.section("Elution-time prior", """
<p>Relations are local in time: a segment says two ions are related <i>at these scans</i>. Partitioning on relations
alone therefore lets a long-lived ion act as a hub: chemical background spans minutes, and one coincidental
relation during a peak's elution attaches that peak to it, so a compound ends up holding members whose signal sits
minutes apart. Ions of one compound elute together, so the partition carries a prior at the run's own peak width:</p>
<ul class="prose">
<li><b>Peak width from the data.</b> Every node gets the intensity-weighted centre and spread of its signal. The
typical spread σ is the median over nodes of confident relations, and the peak width is its FWHM (2.355 σ, a
moment width, so tails count).</li>
<li><b>Bulk containment.</b> The bulk of each node's signal must lie within one peak width of the relation linking
it. An edge is discounted by log f<sub>a</sub> + log f<sub>b</sub>, where f is the share of the node's signal inside
that window. A 10-second relation explains a peak but not a plateau that lasts minutes.</li>
<li><b>Co-elution.</b> Member signal centres scatter around the compound's centre by a learned τ plus each centre's
own standard error (spread / √n<sub>eff</sub>), so peaks near the noise floor need not apex together. Inside GAEC
merging two clusters costs the Ward increment ½·P<sub>A</sub>P<sub>B</sub>/(P<sub>A</sub>+P<sub>B</sub>)·Δcentre²,
and chains of individually good links can no longer reach across the run. τ is fitted by a Gaussian + uniform
mixture to centre differences of confident, contained relations.</li></ul>
<p>The prior acts only on the partition. Relation posteriors and their FDR are unchanged. On synthetic data
long-lived background ions are deliberately no longer joined across time, which is why the background ARI drops
while eluting compounds group as well or better.</p>""")
    page.figure(span_figure(spans), "Member centre span measured against each run's learned peak width; "
                "multi-ion compounds only.")
    page.table(priors.assign(tau_over_sigma=priors.tau / priors.sigma)[
        ["dataset", "sigma", "fwhm", "tau", "tau_over_sigma", "inlier", "n_edges"]],
        "Learned per run: σ = typical intensity-weighted elution spread (s), fwhm = 2.355 σ, τ = centre scatter of "
        "co-eluting members beyond their standard errors (s), inlier = share of confident relations in the "
        "Gaussian component. Simulated ions co-elute exactly (τ ≈ 0); real Astral ions scatter by about σ/2.")
    if split is not None:
        page.figure(split_view(*split))
    page.section("Limitations", """<ul class="prose">
<li>With the elution-time prior, group purity is close to the feature graph's. The remaining impurity comes from
coincidental relations between ions that really do co-elute. Negative evidence between co-eluting blocks with poor
ratio constancy and an isotope-envelope constraint per block would address it.</li>
<li>Long-lived background ions are split into time-local pieces, one per stretch of relations. When background
matters (e.g. for subtraction) its pieces could be re-joined by trace identity in a separate step.</li>
<li>Relations are anchored on the lighter ion; a relation whose lighter partner is too weak to form a trace is
only found if the heavier ion's reverse search is added.</li>
<li>Posteriors are mildly overconfident mid-range (≈0.1) and well calibrated at the extremes.</li></ul>""")
    page.write(REPORTS / "scan_graph.html")
    print(R.round(3).to_string())
    print(S.round(4).to_string())
    print(real.round(4).to_string())


if __name__ == "__main__":
    main()
