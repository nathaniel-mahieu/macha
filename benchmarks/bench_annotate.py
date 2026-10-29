"""Benchmark + explorer: modulo mass-difference graph with Bayesian edge evaluation."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from common import NAV, REAL_LABEL, REPORTS, Timer, real_run, save_json, synthetic
from macha import chem
from macha.annotate import ProfileSource, annotate
from macha.pipeline import process_run
from macha.report import LIGHT, SERIES, Page, binned_hist, method_color, style
from macha.synth.truth import feature_truth
from macha.traces import adjusted_rand

DATASETS = {"syn-metabolomics": False, "syn-hard": False, "syn-proteomics": True, "syn-dense-proteomics": True}
LABELS = {"macha": "Bayesian posterior + graph partition (macha)", "threshold": "rule + RT window only",
          "threshold+corr": "rule + RT window + shape correlation ≥ 0.8"}


def components(n, i, j):
    if len(i) == 0:
        return np.arange(n)
    g = coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n))
    return connected_components(g, directed=False)[1]


def evaluate(name, pep):
    sr = synthetic(name)
    s = process_run(sr.run)
    truth = feature_truth(s.features, s.traces, sr.ion, sr.ions, sr.compounds)
    with Timer() as t:
        ann = annotate(s.features, ProfileSource(s.baselines), peptides=pep, cycle=sr.run.cycle_time)
    comp = truth.compound.to_numpy()
    e = ann.edges.copy()
    e["true"] = (comp[e.i] == comp[e.j]) & (comp[e.i] >= 0)
    n = len(s.features)
    a = comp >= 0
    rows = []
    rules = {
        "macha": e.posterior >= 0.5,
        "threshold": (e.u.abs() <= 3) & (e.dt.abs() <= 2),
        "threshold+corr": (e.u.abs() <= 3) & (e.dt.abs() <= 2) & (e.rho.fillna(0) >= 0.8),
    }
    o = np.argsort(-e.posterior.to_numpy())
    cum = np.cumsum(e.true.to_numpy()[o])
    sweep_prec = cum / np.arange(1, len(o) + 1)
    sweep_rec = cum / e.true.sum()
    for m, acc in rules.items():
        acc = acc.to_numpy()
        prec = e.true[acc].mean()
        rec = (acc & e.true).sum() / e.true.sum()
        lab = ann.labels if m == "macha" else components(n, e.i[acc].to_numpy(), e.j[acc].to_numpy())
        rows.append(dict(dataset=name, method=m, edge_precision=prec, edge_recall=rec, edge_f1=2 * prec * rec / (prec + rec),
                         macha_recall_at_same_precision=float(sweep_rec[sweep_prec >= prec].max()),
                         group_ari=adjusted_rand(lab[a], comp[a]), seconds=t.s if m == "macha" else np.nan))
    r = ann.roles
    ok = a & np.isfinite(r.neutral_mass.to_numpy(float))
    err = np.abs(r.neutral_mass[ok] - truth.neutral_mass[ok]) / truth.neutral_mass[ok] * 1e6
    dz = 2 if pep else 1
    extra = dict(neutral_mass_assigned=ok.sum() / a.sum(), neutral_mass_correct=(err < 10).mean(),
                 charge_correct=(r.charge[a].replace(0, dz) == truth.charge[a]).mean(),
                 isotope_index_correct=(r.isotope[a] == truth.isotope[a]).mean(),
                 n_features=n, n_targets=len(e), n_decoys=len(ann.decoys))
    rows[0].update(extra)
    return rows, e, ann, s, truth


def pr_figure(results):
    fig = go.Figure()
    for k, (name, (e, rows)) in enumerate(results.items()):
        o = np.argsort(-e.posterior.to_numpy())
        tp = np.cumsum(e.true.to_numpy()[o])
        prec = tp / np.arange(1, len(o) + 1)
        rec = tp / e.true.sum()
        step = max(len(o) // 400, 1)
        fig.add_trace(go.Scatter(x=rec[::step], y=prec[::step], mode="lines", name=f"{name}: posterior sweep",
                                 line=dict(width=2, color=SERIES[k]), legendgroup=name))
        for rr in rows:
            if rr["method"] == "macha":
                continue
            fig.add_trace(go.Scatter(x=[rr["edge_recall"]], y=[rr["edge_precision"]], mode="markers",
                                     marker=dict(size=11, color=SERIES[k], symbol="diamond" if rr["method"] ==
                                                 "threshold" else "square", line=dict(width=2, color=LIGHT["surface"])),
                                     name=f"{name}: {LABELS[rr['method']]}", legendgroup=name, showlegend=True))
    fig.update_xaxes(title="edge recall (true same-compound relations found)", range=[0.5, 1.005])
    fig.update_yaxes(title="edge precision")
    return style(fig, 460, title="Precision–recall of candidate edges (lines: posterior threshold sweep; symbols: "
                                 "rule-based baselines)")


def calibration_figure(results):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", line=dict(width=1, color=LIGHT["axis"]),
                             name="perfect calibration"))
    for k, (name, (e, _)) in enumerate(results.items()):
        bins = np.linspace(0, 1, 11)
        b = np.clip(np.digitize(e.posterior, bins) - 1, 0, 9)
        g = e.groupby(b).agg(p=("posterior", "mean"), f=("true", "mean"), n=("true", "size"))
        g = g[g.n >= 10]
        fig.add_trace(go.Scatter(x=g.p, y=g.f, mode="lines+markers", name=name, line=dict(width=2, color=SERIES[k]),
                                 marker=dict(size=8), customdata=g.n,
                                 hovertemplate="mean posterior %{x:.2f}<br>fraction true %{y:.2f}<br>n=%{customdata}"
                                               "<extra></extra>"))
    fig.update_xaxes(title="predicted posterior (binned)", range=[0, 1])
    fig.update_yaxes(title="observed fraction of true edges", range=[0, 1.02])
    return style(fig, 380, title="Calibration: does a posterior of 0.7 mean 70 % of such edges are real?")


def fdr_figure(ann, e=None, title=""):
    f = ann.fdr
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=f.threshold, y=f.bayes_fdr, mode="lines", name="Bayesian FDR (mean 1 − posterior)",
                             line=dict(width=2, color=method_color("macha"))))
    fig.add_trace(go.Scatter(x=f.threshold, y=f.decoy_fdr, mode="lines", name="decoy-estimated FDR",
                             line=dict(width=2, color=method_color("greedy"))))
    if e is not None:
        true_fdr = [1 - e.true[e.posterior >= t].mean() for t in f.threshold]
        fig.add_trace(go.Scatter(x=f.threshold, y=true_fdr, mode="lines", name="true FDR (synthetic truth)",
                                 line=dict(width=2, color=method_color("slice"))))
    fig.update_xaxes(title="posterior threshold")
    fig.update_yaxes(title="false discovery rate", rangemode="tozero")
    return style(fig, 340, title=title)


def evidence_figure(ann, title):
    t, d = ann.edges, ann.decoys
    fig = go.Figure()
    for nm, df, col in (("target edges", t, method_color("macha")), ("decoy edges (null)", d, method_color("greedy"))):
        fig.add_trace(go.Bar(**binned_hist(df.rho, 50, (-1, 1)), name=nm, marker=dict(color=col), opacity=0.7))
    fig.update_layout(barmode="overlay")
    fig.update_xaxes(title="profile shape correlation between the two features' traces")
    fig.update_yaxes(title="density")
    return style(fig, 320, title=title)


def residue_figure(features, ann, rt_window, title, z=1):
    lo, hi = rt_window
    sel = (features.rt >= lo) & (features.rt <= hi)
    f = features[sel]
    idx = np.flatnonzero(sel.to_numpy())
    d = chem.NEUTRON_C13
    res = np.mod(z * f.mz.to_numpy(), d)
    fig = go.Figure()
    e = ann.edges
    e = e[(e.target == "13C") & (e.z == z) & (e.posterior >= 0.9) & e.i.isin(idx) & e.j.isin(idx)]
    mz = features.mz.to_numpy()
    xs, ys = [], []
    for i, j in zip(e.i, e.j):
        xs += [mz[i], mz[j], None]
        ys += [np.mod(z * mz[i], d), np.mod(z * mz[j], d), None]
    fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(width=2, color=method_color("macha")),
                             name=f"accepted 13C edges (posterior ≥ 0.9, z = {z})", hoverinfo="skip"))
    fig.add_trace(go.Scattergl(x=f.mz, y=res, mode="markers", name="features",
                               marker=dict(size=5, color=LIGHT["muted"]),
                               hovertemplate="m/z %{x:.4f}<br>residue %{y:.4f}<extra></extra>"))
    fig.update_xaxes(title="m/z")
    fig.update_yaxes(title=("" if z == 1 else f"{z}·") + "m/z mod 1.003355 (13C residue)")
    return style(fig, 420, title=title)


def group_figure(features, ann, g, title):
    r = ann.roles
    mem = np.flatnonzero(ann.labels == g)
    f = features.iloc[mem]
    roles = r.iloc[mem]
    fig = go.Figure()
    kinds = np.where(roles.isotope.to_numpy() > 0, "isotopologue", "monoisotopic")
    for k, col in (("monoisotopic", method_color("macha")), ("isotopologue", method_color("slice"))):
        m = kinds == k
        xs, ys = [], []
        for x, y in zip(f.mz[m], f.height[m]):
            xs += [x, x, None]
            ys += [0, y, None]
        fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(width=2, color=col), name=k, hoverinfo="skip"))
        lbl = [f"{a or '?'} M+{i} z={z}" for a, i, z in zip(roles.adduct[m], roles.isotope[m], roles.charge[m])]
        fig.add_trace(go.Scatter(x=f.mz[m], y=f.height[m], mode="markers", marker=dict(size=8, color=col,
                      line=dict(width=2, color=LIGHT["surface"])), showlegend=False, text=lbl,
                      hovertemplate="m/z %{x:.4f}<br>height %{y:.3s}<br>%{text}<extra></extra>"))
    e = ann.edges[ann.edges.i.isin(mem) & ann.edges.j.isin(mem) & (ann.edges.posterior >= 0.5)]
    top = float(f.height.max()) * 1.08
    for k, row in enumerate(e.sort_values("posterior", ascending=False).head(14).itertuples()):
        x0, x1 = features.mz.iat[row.i], features.mz.iat[row.j]
        yy = top * (1 + 0.05 * (k % 6))
        fig.add_trace(go.Scatter(x=[x0, x0, x1, x1], y=[features.height.iat[row.i], yy, yy,
                                                         features.height.iat[row.j]],
                                 mode="lines", line=dict(width=1, color=LIGHT["axis"]), showlegend=False,
                                 hovertemplate=f"{row.rule} z={row.z} n={row.n}<br>posterior {row.posterior:.3f}"
                                               f"<br>ρ={row.rho:.2f}<extra></extra>"))
    fig.update_xaxes(title="m/z")
    fig.update_yaxes(title="apex height")
    return style(fig, 340, title=title)


def main():
    rows = []
    results = {}
    keep = {}
    for name, pep in DATASETS.items():
        rr, e, ann, s, truth = evaluate(name, pep)
        rows += rr
        results[name] = (e, rr)
        keep[name] = (ann, s, truth)
        print(name, "done")
    df = pd.DataFrame(rows)
    page = Page("Annotation graph", "Isotopes, adducts, neutral losses, charge states and multimers are discovered "
                "with modulo arithmetic on charge-scaled masses, every candidate edge is scored by a two-groups "
                "empirical-Bayes model whose null is learned from decoy rules, and compounds are the MAP partition "
                "of the resulting graph.", NAV)
    m = df[df.method == "macha"].set_index("dataset")
    page.tiles([
        ("edge precision at posterior ≥ 0.5", f"{m.edge_precision.min():.3f}–{m.edge_precision.max():.3f}",
         "across four synthetic datasets"),
        ("compound grouping ARI", f"{m.group_ari.min():.3f}–{m.group_ari.max():.3f}",
         "rule + RT only: " + f"{df[df.method == 'threshold'].group_ari.min():.2f}–"
                              f"{df[df.method == 'threshold'].group_ari.max():.2f}"),
        ("neutral mass correct (<10 ppm)", f"{m.neutral_mass_correct.min():.3f}–{m.neutral_mass_correct.max():.3f}",
         "of features given a neutral mass"),
    ])
    page.section("What it does", """
<p><b>Modulo candidates (mz.unity, generalised).</b> For a rule with mass step Δ (¹³C, ¹⁵N, ³⁴S, ³⁷Cl, ⁸¹Br, ¹⁸O,
⁴¹K, Na–H, K–H, NH₃, H₂O, CO₂, …) and charge z, the charge-scaled masses M = z·m/z of related features differ by
n·Δ, so they share the residue M mod Δ and n = round(ΔM/Δ) counts the steps; chains of any length (M, M+1, M+2, …) fall
out of one search. Only canonical edges with gcd(n, z) = 1 are kept, so a z = 2 claim needs its own half-spaced peaks.
Charge states and multimers, which are not constant mass differences, are matched by mapping both features to a
neutral mass. Homologous series (CH₂, C₂H₄O) use the same residue test without co-elution. Co-eluting candidates
are found by binary search of n·Δ/z targets in m/z-sorted RT blocks, which gives exactly the edges of a pairwise
residue sweep at a fraction of the cost.</p>
<p><b>Bayesian edge evaluation.</b> Every edge carries evidence: mass residual / its propagated σ, RT difference / apex
uncertainty, Pearson correlation of the two raw traces over the union of both peak windows, the ¹³C intensity ratio
against a Poisson carbon-count expectation, and (for consensus features) cross-sample intensity correlation. Under
H₀ (coincidence) the evidence density f₀ is learned from <b>decoy rules</b>, the same search with Δ shifted to values no
chemistry produces, with co-elution and shape modelled jointly because co-eluting coincidences are correlated more
often. Under H₁ the densities are parametric. EM estimates a prior π per (rule, charge) and the H₁ parameters; each
edge gets a posterior and a log-odds. Applying the fitted model to the decoys yields an independent FDR estimate,
which checks calibration on real data where there is no truth.</p>
<p><b>Graph inference.</b> Compound groups are the MAP partition under the edge log-odds (correlation clustering by
greedy additive edge contraction), so a strongly negative edge can veto a merge that a single weak positive edge
would cause in connected components. Within a group, isotope envelopes are chosen greedily by total log-odds
(non-overlapping, per charge), other isotopes attach as fine-structure peaks, peptide envelopes are checked against
an averagine pattern for a missing monoisotopic peak, and adducts / neutral mass come from the hypothesis explaining
most monoisotopic features.</p>""")
    page.section("Accuracy against ground truth", "<p>An edge is true when both features come from the same "
                 "compound. Baselines use the same candidates: \"rule + RT\" is the CAMERA-style deterministic rule, "
                 "and adding a correlation cut is CAMERA with correlation grouping; groups are their connected "
                 "components.</p>")
    page.figure(pr_figure(results))
    tab = df[["dataset", "method", "edge_precision", "edge_recall", "macha_recall_at_same_precision", "edge_f1",
              "group_ari"]].copy()
    tab["method"] = tab["method"].map(LABELS)
    page.table(tab, "<b>macha_recall_at_same_precision</b>: recall of the posterior ranking when its threshold is set to "
               "match that row's precision. It is the like-for-like comparison: the posterior ranking dominates every "
               "rule-based baseline. On syn-hard (4-15 ppm co-eluting interferents, many coalesced) the default "
               "posterior cut of 0.5 is conservative, because coalesced features degrade mass, RT and shape evidence "
               "together.", better=dict(edge_precision="max", edge_recall="max", edge_f1="max", group_ari="max"),
               group_col="dataset")
    page.table(m[["n_features", "n_targets", "n_decoys", "neutral_mass_assigned", "neutral_mass_correct",
                  "charge_correct", "isotope_index_correct", "seconds"]].reset_index(),
               "Compound inference (macha only)")
    page.figure(calibration_figure(results))
    name = "syn-dense-proteomics"
    page.figure(fdr_figure(keep[name][0], results[name][0], f"FDR estimates vs truth, {name}"))
    # explorer: groups on synthetic proteomics
    ann, s, truth = keep["syn-proteomics"]
    sizes = pd.Series(ann.labels).value_counts()
    page.section("Explorer: compound groups", "<p>Stick spectra of inferred groups (hover for roles and edges). "
                 "Brackets join features linked by an accepted edge; hover a bracket for its rule, posterior and "
                 "shape correlation.</p>")
    for g in sizes.index[[0, 3, 10]]:
        nm = ann.groups.set_index("group").loc[g]
        page.figure(group_figure(s.features, ann, g, f"syn-proteomics group {g}: neutral mass "
                                                     f"{nm.neutral_mass:.4f}, {nm.n_features} features, "
                                                     f"adducts {nm.adducts or '–'}"))
    real_rows = []
    for ds, pep in (("proteomics_sc", True), ("proteomics", True), ("metabolomics", False)):
        run = real_run(ds, 0)
        if run is None:
            continue
        s = process_run(run)
        with Timer() as t:
            ann = annotate(s.features, ProfileSource(s.baselines), peptides=pep, cycle=run.cycle_time)
        f9 = ann.fdr.iloc[(ann.fdr.threshold - 0.9).abs().argmin()]
        z = ann.roles.charge
        real_rows.append(dict(dataset=REAL_LABEL[ds], features=len(s.features), targets=len(ann.edges),
                              decoys=len(ann.decoys), accepted_090=int(f9.accepted), bayes_fdr_090=f9.bayes_fdr,
                              decoy_fdr_090=f9.decoy_fdr, multi_feature_groups=int((ann.groups.n_features > 1).sum()),
                              envelopes_z1=int((z == 1).sum()), envelopes_z2plus=int((z >= 2).sum()), seconds=t.s))
        if ds == "proteomics_sc":
            real_ann, real_s = ann, s
        if ds == "metabolomics":
            met_ann, met_s = ann, s
    if real_rows:
        page.section("Real Orbitrap Astral data", "<p>No truth: the decoy-estimated FDR is the check. Where both "
                     "estimates agree, the posterior is calibrated on real data.</p>")
        page.table(pd.DataFrame(real_rows), better=dict(seconds="min"))
        page.figure(fdr_figure(real_ann, None, f"Bayesian vs decoy FDR, {REAL_LABEL['proteomics_sc']}"))
        page.figure(evidence_figure(real_ann, "Shape correlation separates real relations from decoy coincidences"))
        page.figure(residue_figure(real_s.features, real_ann, (700, 720),
                                   f"The modulo view, {REAL_LABEL['proteomics_sc']}, RT 700-720 s: doubly charged "
                                   "isotope partners share the residue of 2·m/z, so accepted edges are horizontal", z=2))
        if "met_s" in locals():
            page.figure(residue_figure(met_s.features, met_ann, (300, 330),
                                       f"The modulo view, {REAL_LABEL['metabolomics']}, RT 300-330 s (z = 1)", z=1))
    save_json("annotate", dict(synthetic=df.to_dict("records"), real=real_rows))
    page.write(REPORTS / "annotate.html")
    print(df.round(4).to_string())
    print(pd.DataFrame(real_rows).round(4).to_string())


if __name__ == "__main__":
    main()
