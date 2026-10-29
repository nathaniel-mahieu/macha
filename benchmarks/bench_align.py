"""Benchmark + explorer: drift-trend mixture alignment with probabilistic correspondence."""

from __future__ import annotations

import pickle

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from common import CACHE, NAV, REAL_LABEL, REPORTS, Timer, real_runs, save_json
from macha import MSRun
from macha.align import align, align_dtw, density_groups
from macha.align.metrics import correspondence_metrics, local_class_agreement, rt_dispersion
from macha.annotate import ProfileSource, annotate
from macha.pipeline import process_run
from macha.report import LIGHT, SERIES, Page, method_color, style
from macha.synth import make_batch, metabolomics_scenario, proteomics_scenario, true_warps
from macha.synth.truth import feature_truth

BATCHES = {
    "one-trend": dict(sc=lambda: metabolomics_scenario(n_compounds=500, n_isomers=250), n_samples=8,
                      n_drift_classes=1, class_drift=12.0, global_drift=6.0,
                      note="8 runs, a single drift trend (classic global warp is the right model)"),
    "three-trends": dict(sc=lambda: metabolomics_scenario(n_compounds=500, n_isomers=250), n_samples=8,
                         n_drift_classes=3, class_drift=16.0, global_drift=6.0,
                         note="8 runs, 3 latent compound classes with different drift trends, many isomers"),
    "three-trends-strong": dict(sc=lambda: metabolomics_scenario(n_compounds=500, n_isomers=250), n_samples=8,
                                n_drift_classes=3, class_drift=30.0, global_drift=8.0,
                                note="as above with drift up to ±15 s between first and last injection"),
    "peptides-two-trends": dict(sc=lambda: proteomics_scenario(n_compounds=1500, n_background=20), n_samples=6,
                                n_drift_classes=2, class_drift=16.0, global_drift=6.0,
                                note="6 runs, 1500 peptides (z 1-5), 2 drift classes"),
}
LABELS = {"density": "density grouping, no correction (xcms-like)", "none": "probabilistic matching, no correction",
          "dtw": "BPC dynamic time warping + matching", "global": "single smooth warp per run + matching",
          "macha": "drift-trend mixture + Bayesian membership (macha)"}


def load_batch(name):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"batch-{name}.pkl"
    if path.exists():
        return pickle.loads(path.read_bytes())
    spec = BATCHES[name]
    b = make_batch(spec["sc"](), n_samples=spec["n_samples"], n_drift_classes=spec["n_drift_classes"],
                   class_drift=spec["class_drift"], global_drift=spec["global_drift"], seed=3)
    samples = [process_run(sr.run) for sr in b.runs]
    truths = [feature_truth(s.features, s.traces, sr.ion, sr.ions, sr.compounds) for s, sr in zip(samples, b.runs)]
    tables = [s.features for s in samples]
    chrom = [(sr.run.rt, sr.run.bpc()) for sr in b.runs]
    grid = np.linspace(0, b.scenario.run_time, 600)
    out = dict(tables=tables, truths=truths, chrom=chrom, compounds=b.compounds, drift_class=b.drift_class,
               injection_order=b.injection_order, true_grid=grid, true_W=true_warps(b, grid),
               run_time=b.scenario.run_time)
    path.write_bytes(pickle.dumps(out))
    return out


def evaluate(name):
    d = load_batch(name)
    tables, truths = d["tables"], d["truths"]
    ion = [np.where(t["kind"].to_numpy() == "background", -2, t["ion"].to_numpy()) for t in truths]
    comp = [t["compound"].to_numpy() for t in truths]
    rows = []
    fits = {}
    for m in LABELS:
        with Timer() as t:
            if m == "density":
                groups, aligned, al = density_groups(tables, bw=5.0), [x["rt"].to_numpy() for x in tables], None
            elif m == "dtw":
                al = align_dtw(tables, d["chrom"])
                groups, aligned = al.feature_group, al.aligned_rt
            else:
                al = align(tables, **({"drift": "none"} if m == "none" else {"K": 1} if m == "global" else {}))
                groups, aligned = al.feature_group, al.aligned_rt
        met = correspondence_metrics(groups, ion)
        row = dict(batch=name, method=m, **met, rt_rms_s=rt_dispersion(aligned, ion), seconds=t.s)
        if m in ("macha", "global") and al is not None and al.drift is not None:
            row["K"] = al.drift.K
            row["sigma_s"] = al.drift.sigma
            g = np.concatenate(al.feature_group)
            c = np.concatenate(comp)
            ok = c >= 0
            df = pd.DataFrame(dict(g=g[ok], c=c[ok]))
            maj = df.groupby("g")["c"].agg(lambda x: x.value_counts().index[0])
            gcls = np.full(len(al.groups), -1)
            gcls[maj.index.to_numpy()] = d["drift_class"][maj.to_numpy()]
            if d["drift_class"].max() > 0:
                row["local_class_agreement"] = local_class_agreement(al.groups["drift_class"].to_numpy(), gcls,
                                                                     al.groups["rt"].to_numpy())
        rows.append(row)
        fits[m] = al
    return rows, fits, d, ion


def drift_explorer(al, d, title):
    """Per-run anchor deviations and learned class warps vs truth; a menu picks the run."""
    model = al.drift
    T = model.mu
    A = al.anchors
    resp = model.resp.argmax(1)
    grid = np.linspace(model.knots[0], model.knots[-1], 300)
    W = model.warps(grid)
    TW = np.stack([[np.interp(grid, d["true_grid"], w) for w in Wk] for Wk in d["true_W"]])
    J = model.n_samples
    fig = go.Figure()
    groups = []
    for j in range(J):
        start = len(fig.data)
        # anchors: deviation of this run from consensus, coloured by inferred class
        dev = A.rt[:, j] - T if A.rt.shape[0] == T.shape[0] else None
        for k in range(model.K):
            m = (resp == k) & np.isfinite(dev) if dev is not None else np.zeros(0, bool)
            if dev is not None:
                fig.add_trace(go.Scattergl(x=T[m], y=dev[m], mode="markers", name=f"anchors, inferred class {k}",
                                           marker=dict(size=5, color=SERIES[k], opacity=0.55), visible=j == 0,
                                           legendgroup=f"c{k}"))
            fig.add_trace(go.Scatter(x=grid, y=W[k, j], mode="lines", name=f"learned warp, class {k}",
                                     line=dict(width=2, color=SERIES[k]), visible=j == 0, legendgroup=f"c{k}"))
        for k in range(TW.shape[0]):
            fig.add_trace(go.Scatter(x=grid, y=TW[k, j], mode="lines", name="true class warps" if k == 0 else None,
                                     showlegend=k == 0, line=dict(width=1, color=LIGHT["ink"]), visible=j == 0,
                                     hoverinfo="skip"))
        groups.append((start, len(fig.data)))
    buttons = []
    for j, (a, b) in enumerate(groups):
        vis = [a <= q < b for q in range(len(fig.data))]
        buttons.append(dict(label=f"run {j + 1} (injection order {d['injection_order'][j]:+.2f})", method="update",
                            args=[{"visible": vis}]))
    fig.update_layout(updatemenus=[dict(buttons=buttons, active=0, x=1, xanchor="right", y=1.2, yanchor="top",
                                        bgcolor=LIGHT["surface"], bordercolor=LIGHT["axis"],
                                        font=dict(color=LIGHT["ink"]))])
    fig.update_xaxes(title="consensus retention time (s)")
    fig.update_yaxes(title="RT in this run − consensus (s)")
    return style(fig, 440, title=title)


def membership_figure(results):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode="lines", line=dict(width=1, color=LIGHT["axis"]),
                             name="perfect calibration"))
    for k, (name, (al, ion)) in enumerate(results.items()):
        g = np.concatenate(al.feature_group)
        p = np.concatenate(al.membership)
        i = np.concatenate(ion)
        ok = i >= 0
        df = pd.DataFrame(dict(g=g[ok], i=i[ok], p=p[ok]))
        maj = df.groupby("g")["i"].agg(lambda x: x.value_counts().index[0])
        size = df.groupby("g").size()
        df["correct"] = (df["i"].to_numpy() == maj.loc[df["g"]].to_numpy()) & (size.loc[df["g"]].to_numpy() > 1)
        bins = np.linspace(0, 1, 11)
        b = np.clip(np.digitize(df.p, bins) - 1, 0, 9)
        s = df.groupby(b).agg(p=("p", "mean"), f=("correct", "mean"), n=("p", "size"))
        s = s[s.n >= 20]
        fig.add_trace(go.Scatter(x=s.p, y=s.f, mode="lines+markers", name=name, line=dict(width=2, color=SERIES[k]),
                                 marker=dict(size=8), customdata=s.n,
                                 hovertemplate="membership %{x:.2f}<br>correctly grouped %{y:.2f}<br>n=%{customdata}"
                                               "<extra></extra>"))
    fig.update_xaxes(title="membership probability (binned)", range=[0, 1])
    fig.update_yaxes(title="fraction placed with its true ion", range=[0, 1.02])
    return style(fig, 360, title="Membership calibration")


def real_alignment(dataset, max_runs=6):
    files = real_runs(dataset)[:max_runs]
    if len(files) < 2:
        return None
    samples = [process_run(MSRun.load(f)) for f in files]
    tables = [s.features for s in samples]
    out = {}
    for m, kw in (("none", dict(drift="none")), ("global", dict(K=1)), ("macha", {})):
        with Timer() as t:
            al = align(tables, **kw)
        present = (al.feature_index >= 0).sum(1)
        # isotope consistency: 13C partners in run 1 should land in groups with the same consensus RT
        ann = annotate(tables[0], ProfileSource(samples[0].baselines), peptides=dataset != "metabolomics",
                       cycle=samples[0].run.cycle_time) if m == "none" else out["none"]["ann"]
        e = ann.edges[(ann.edges.target == "13C") & (ann.edges.posterior > 0.99)]
        g0 = al.feature_group[0]
        grt = al.groups["rt"].to_numpy()
        iso = np.abs(grt[g0[e.i]] - grt[g0[e.j]])
        both = (present[g0[e.i]] >= len(files) - 1) & (present[g0[e.j]] >= len(files) - 1)
        out[m] = dict(al=al, ann=ann, row=dict(
            dataset=REAL_LABEL[dataset], runs=len(files), method=m, features=int(sum(len(t) for t in tables)),
            groups=len(al.groups), groups_in_all_runs=int((present == len(files)).sum()),
            K=al.drift.K if al.drift is not None else np.nan,
            sigma_s=al.drift.sigma if al.drift is not None else np.nan,
            anchors=al.anchors.rt.shape[0], isotope_rt_disagreement_s=float(np.median(iso[both])) if both.any() else np.nan,
            isotope_pairs_same_group_pattern=float(np.mean([
                np.array_equal(al.feature_index[g0[i]] >= 0, al.feature_index[g0[j]] >= 0)
                for i, j in zip(e.i[both][:3000], e.j[both][:3000])])) if both.any() else np.nan,
            seconds=t.s))
    return out


def main():
    rows = []
    keep = {}
    for name in BATCHES:
        r, fits, d, ion = evaluate(name)
        rows += r
        keep[name] = (fits, d, ion)
        print(name, "done")
    df = pd.DataFrame(rows)
    page = Page("Alignment", "Retention-time drift is modelled as a mixture of drift trends shared by subsets of "
                "compounds, learned by EM on unambiguous anchors; correspondence is then solved in each m/z channel "
                "with the posterior predictive RT of each group, giving every feature a membership probability.", NAV)
    m = df[df.method == "macha"].set_index("batch")
    g = df[df.method == "global"].set_index("batch")
    page.tiles([
        ("aligned RT spread, three-trends-strong", f"{m.loc['three-trends-strong', 'rt_rms_s']:.2f} s",
         f"global warp {g.loc['three-trends-strong', 'rt_rms_s']:.2f} s · raw "
         f"{df[(df.batch == 'three-trends-strong') & (df.method == 'none')].rt_rms_s.iloc[0]:.2f} s"),
        ("correspondence F1, three-trends-strong", f"{m.loc['three-trends-strong', 'f1']:.3f}",
         f"global {g.loc['three-trends-strong', 'f1']:.3f} · density "
         f"{df[(df.batch == 'three-trends-strong') & (df.method == 'density')].f1.iloc[0]:.3f}"),
        ("drift classes found (BIC)", " / ".join(str(int(m.loc[b, 'K'])) for b in BATCHES),
         "true: " + " / ".join(str(BATCHES[b]["n_drift_classes"]) for b in BATCHES)),
    ])
    page.section("What it does", """
<p><b>Drift trends, not m/z channels.</b> The original warpgroup/warpcombine aligned each m/z group separately with DTW
(noisy, expensive, blind to what other compounds do), and <code>grtgmzcor</code> fitted one global curve per run. Here
each compound belongs to a latent <i>drift class</i>; its RT in run j is μ + w<sub>kj</sub>(μ), where w<sub>kj</sub> is a
P-spline warp shared by every compound of class k, centred across runs so μ is the consensus time. EM on anchors
(features unambiguous in m/z within a wide RT window, matched across runs) alternates class responsibilities,
penalised weighted least squares per (class, run), consensus times, a heteroscedastic noise scale (each anchor's apex
uncertainty from its width and SNR) and an outlier component for gross mismatches. BIC picks the number of
classes. Per-run m/z calibration (smooth ppm offset vs RT) is fitted on the same anchors.</p>
<p><b>Bayesian group membership.</b> Pooled, calibrated features are cut into m/z channels. In each channel, groups are
grown run by run and then refined by leaving each run out and re-assigning it, always by linear assignment, with at
most one feature per run and group. The cost is the negative log posterior predictive density: m/z is Gaussian, and
RT follows a Student-t whose centre and spread come from the mixture over drift classes. A group's class posterior
comes from how well its members fit each class, so a young or ambiguous group predicts broadly and a well-determined
one predicts tightly. Every feature then receives a leave-one-out membership probability against all compatible
alternatives, including "belongs to no group". After a first pass, the drift model is re-learned from all confident
groups (many more than the anchors) and correspondence is solved again.</p>
<p><b>Identifiability.</b> Where two class trends cross, class labels can swap between RT regions without changing
any prediction. Class recovery is therefore scored locally (ARI within 60 s windows); the RT predictions themselves
are unaffected.</p>""")
    page.section("Accuracy against ground truth", "<p>Correspondence is scored on pairs of features from different "
                 "runs: a pair is correct when both come from the same true ion. RT spread is the RMS deviation of "
                 "each method's aligned RTs around their true ion's mean.</p>")
    for metric, ttl in (("f1", "Correspondence F1"), ("rt_rms_s", "Aligned RT spread (s, lower is better)")):
        fig = go.Figure()
        for meth in LABELS:
            dd = df[df.method == meth].set_index("batch").reindex(list(BATCHES))
            fig.add_trace(go.Scatter(x=dd[metric], y=list(BATCHES), mode="markers", name=LABELS[meth],
                                     marker=dict(size=12, color=method_color(meth),
                                                 line=dict(width=2, color=LIGHT["surface"]))))
        fig.update_yaxes(autorange="reversed")
        fig.update_xaxes(title=metric, type="log" if metric == "rt_rms_s" else "linear")
        page.figure(style(fig, 280, title=ttl))
    tab = df[["batch", "method", "precision", "recall", "f1", "complete_ions", "rt_rms_s", "K", "sigma_s",
              "local_class_agreement", "seconds"]].copy()
    tab["method"] = tab["method"].map(LABELS)
    page.table(tab, "Batches: " + "; ".join(f"<b>{k}</b>: {v['note']}" for k, v in BATCHES.items()),
               better=dict(precision="max", recall="max", f1="max", complete_ions="max", rt_rms_s="min"),
               group_col="batch")
    page.section("Explorer: learned drift trends", "<p>Anchors of each run plotted as their deviation from the "
                 "consensus time, coloured by inferred class, with the learned class warps (coloured lines) and the "
                 "simulator's true class warps (thin dark lines). Pick a run from the menu; early and late "
                 "injections drift in opposite directions.</p>")
    fits, d, ion = keep["three-trends-strong"]
    page.figure(drift_explorer(fits["macha"], d, "three-trends-strong: anchors and warps per run"))
    bic = []
    for name in BATCHES:
        for f in keep[name][0]["macha"].drift_fits:
            bic.append(dict(batch=name, K=f.K, bic=f.bic))
    bic = pd.DataFrame(bic)
    fig = go.Figure()
    for k, name in enumerate(BATCHES):
        dd = bic[bic.batch == name]
        rel = dd.bic - dd.bic.min()
        fig.add_trace(go.Scatter(x=dd.K, y=rel, mode="lines+markers", name=name, line=dict(width=2, color=SERIES[k]),
                                 marker=dict(size=8)))
    fig.update_xaxes(title="number of drift classes K", dtick=1)
    fig.update_yaxes(title="BIC − min BIC", type="log")
    page.figure(style(fig, 320, title="Model selection: BIC by number of drift classes (minimum = selected)"))
    page.figure(membership_figure({n: (keep[n][0]["macha"], keep[n][2]) for n in BATCHES}))
    real_rows = []
    real_fig = None
    for ds in ("proteomics_sc", "proteomics", "metabolomics"):
        res = real_alignment(ds)
        if res is None:
            continue
        real_rows += [res[k]["row"] for k in ("none", "global", "macha")]
        if real_fig is None:
            al = res["macha"]["al"]
            model = al.drift
            grid = np.linspace(model.knots[0], model.knots[-1], 300)
            W = model.warps(grid)
            fig = go.Figure()
            for k in range(model.K):
                for j in range(model.n_samples):
                    fig.add_trace(go.Scatter(x=grid, y=W[k, j], mode="lines", line=dict(width=2, color=SERIES[k]),
                                             name=f"class {k} (π={model.pi[k]:.2f})", legendgroup=f"c{k}",
                                             showlegend=j == 0, hovertemplate=f"class {k}, run {j + 1}<br>RT %{{x:.0f}} s"
                                                                              "<br>shift %{y:.2f} s<extra></extra>"))
            fig.update_xaxes(title="consensus retention time (s)")
            fig.update_yaxes(title="learned shift (s)")
            real_fig = style(fig, 360, title=f"Learned drift trends, {REAL_LABEL[ds]} ({model.n_samples} runs, "
                                             f"K = {model.K})")
    if real_rows:
        page.section("Real Orbitrap Astral data", "<p>No ground truth. Isotope consistency is the internal check: "
                     "¹³C partners (posterior > 0.99 in run 1) are the same molecule, so their groups should have the "
                     "same consensus RT and the same presence pattern across runs.</p>")
        rt = pd.DataFrame(real_rows)
        rt["method"] = rt["method"].map(LABELS)
        page.table(rt, better=dict(groups_in_all_runs="max", isotope_rt_disagreement_s="min",
                                   isotope_pairs_same_group_pattern="max", seconds="min"), group_col="dataset")
        page.figure(real_fig)
    save_json("align", dict(synthetic=df.to_dict("records"), real=real_rows, bic=bic.to_dict("records")))
    page.write(REPORTS / "align.html")
    print(df.round(3).to_string())
    print(pd.DataFrame(real_rows).round(3).to_string())


if __name__ == "__main__":
    main()
