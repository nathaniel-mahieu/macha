"""Benchmark + explorer: cross-run alignment from the scan-level graph, before peak detection.

First step of grouping across runs at the scan-graph stage: can the drift-trend mixture be learned without picking
a single peak?  Relation anchors (confident scan-graph relations matched across runs by rule and the m/z of both
ions) are compared with feature anchors (picked peaks unique in m/z) on the alignment batches, which carry known
drift classes and a per-compound RT error shared by all of a compound's ions, and on real Astral runs, where the
two anchor sets predict each other.  Noise is measured per ROI, in each ROI's own m/z acceptance window.
"""

from __future__ import annotations

import pickle
from types import SimpleNamespace

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from bench_align import BATCHES
from common import CACHE, NAV, REAL_LABEL, REPORTS, Timer, real_runs, save_json
from macha import MSRun
from macha.align import align, find_anchors, fit_drift, select_drift
from macha.align.anchors import Anchors
from macha.align.drift import cv_drift, support_range
from macha.align.edges import drift_window, relation_anchors, relation_groups, relation_noise, relation_table
from macha.align.metrics import correspondence_metrics, design_mixture, loo_prediction_error, rt_dispersion
from macha.annotate.scan_graph import scan_graph
from macha.baseline import trace_baselines
from macha.peaks import detect_peaks
from macha.report import Page, method_color, style
from macha.synth import make_batch
from macha.synth.truth import feature_truth
from macha.traces import Traces, detect_traces, roi_noise
from macha.traces.noise import floor_halfwidth_ppm

LABELS = {"features": "feature anchors (picked peaks)", "relations": "relation anchors (scan graph, no peaks)"}
FITS = {"bic": "BIC, λ = 1 (current default)", "cv": "K and λ by blocked cross-validation"}


def relation_truth(rel, traces, compound_of_point):
    """True compound of each relation: the majority compound of both ions' centroids in its scans (-2 if they
    disagree, -1 if noise)."""
    scan = traces.run.scan_of_point
    out = np.full(len(rel), -1, np.int64)
    for k, (a, p, s0, s1) in enumerate(zip(rel["anchor"], rel["partner"], rel["s0"], rel["s1"])):
        cs = []
        for t in (a, p):
            pts = traces.points(int(t))
            sc = scan[pts]
            c = compound_of_point[pts[(sc >= s0) & (sc <= s1)]]
            c = c[c >= 0]
            cs.append(int(np.bincount(c).argmax()) if c.size else -1)
        out[k] = cs[0] if cs[0] == cs[1] else -2
    return out


def process(run, peptides):
    """Both paths for one run: traces -> scan graph -> relations, and traces -> baselines -> peaks."""
    with Timer() as t_tr:
        tr = detect_traces(run)
    with Timer() as t_g:
        g = scan_graph(tr, peptides=peptides)
        rel = relation_table(g, tr)
        rn = relation_noise(rel, tr)
    with Timer() as t_p:
        feats = detect_peaks(tr, trace_baselines(tr))
    return tr, g, rel, rn, feats, dict(traces=t_tr.s, scan_graph=t_g.s, peaks=t_p.s)


def load_batch(name):
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"xrun-{name}.pkl"
    if path.exists():
        return pickle.loads(path.read_bytes())
    spec = BATCHES[name]
    b = make_batch(spec["sc"](), n_samples=spec["n_samples"], n_drift_classes=spec["n_drift_classes"],
                   class_drift=spec["class_drift"], global_drift=spec["global_drift"], seed=3)
    peptides = b.scenario.kind == "peptide"
    rels, noise, tables, truths, timings = [], [], [], [], []
    for sr in b.runs:
        tr, g, rel, rn, feats, tm = process(sr.run, peptides)
        rel["true_compound"] = relation_truth(rel, tr, sr.compound)
        rels.append(rel)
        noise.append(rn)
        tables.append(feats)
        truths.append(feature_truth(feats, tr, sr.ion, sr.ions, sr.compounds))
        timings.append(tm)
        print(name, sr.run.meta.get("name"), {k: round(v, 1) for k, v in tm.items()}, len(rel), "relations")
    out = dict(rels=rels, noise=noise, tables=tables, truths=truths, timings=timings, rt=b.rt,
               drift_class=b.drift_class, note=spec["note"])
    path.write_bytes(pickle.dumps(out))
    return out


def anchor_truth(A, members, d, kind):
    """Per anchor: majority true compound and purity (all members from one compound)."""
    C, J = A.rt.shape
    comp = np.full((C, J), -1)
    if kind == "features":
        for j in range(J):
            ok = A.feature[:, j] >= 0
            comp[ok, j] = d["truths"][j]["compound"].to_numpy()[A.feature[ok, j]]
    else:
        tc = np.array([d["rels"][j]["true_compound"].to_numpy()[i] for j, i in zip(members.run, members.rel)])
        comp[members.anchor_id.to_numpy(), members.run.to_numpy()] = tc
    maj = np.array([np.bincount(r[r >= 0]).argmax() if (r >= 0).any() else -1 for r in comp])
    pure = np.array([(r[r != -1] == m).all() and m >= 0 for r, m in zip(comp, maj)])
    return maj, pure


def anchor_quality(A, maj, pure, RT, dc):
    """Measurement error of anchor centres against true RTs (the compound's RT error cancels), and SE calibration."""
    elut = (maj >= 0) & (dc[np.maximum(maj, 0)] >= 0)
    res = A.rt[elut] - RT[maj[elut]]
    res = res - np.nanmean(res, 1, keepdims=True)
    n = np.isfinite(res).sum(1, keepdims=True)
    res = res * np.sqrt(n / np.maximum(n - 1, 1))
    z = res / A.se[elut] if A.se is not None else np.full_like(res, np.nan)
    r = res[np.isfinite(res)]
    zz = z[np.isfinite(z)]
    return dict(anchors=A.rt.shape[0], compounds=int(len(set(maj[elut]))), purity=float(pure.mean()),
                background=float(np.mean((maj >= 0) & (dc[np.maximum(maj, 0)] < 0))),
                centre_error_s=float(1.4826 * np.median(np.abs(r))), gross_gt_1s=float(np.mean(np.abs(r) > 1.0)),
                se_calibration=float(1.4826 * np.median(np.abs(zz))) if zz.size else np.nan)


def evaluate(name):
    d = load_batch(name)
    RT, dc = d["rt"], d["drift_class"]
    true = RT[dc >= 0]
    tables, truths = d["tables"], d["truths"]
    ion = [np.where(t["kind"].to_numpy() == "background", -2, t["ion"].to_numpy()) for t in truths]
    Af = find_anchors(tables)
    Ar, tab, mem = relation_anchors(d["rels"], noise=d["noise"], rt_window=drift_window(d["rels"], d["noise"]))
    sets = {"features": (Af, None), "relations": (Ar, mem)}
    rows, quality, models = [], [], {}
    for kind, (A, members) in sets.items():
        maj, pure = anchor_truth(A, members, d, kind)
        quality.append(dict(batch=name, anchors_from=kind, **anchor_quality(A, maj, pure, RT, dc)))
        for fit in FITS:
            with Timer() as t:
                if fit == "bic":
                    m, _ = select_drift(A.rt, 4, t_range=support_range(A.rt), n_basis=12, lam=1.0, S=A.se)
                    lam = 1.0
                else:
                    m, cv = cv_drift(A.rt, A.se, n_basis=12)
                    lam = float(cv.loc[cv.rmse.idxmin(), "lam"])
            e = loo_prediction_error(m, true)
            al = align(tables, anchors=A, K=m.K, lam=lam, refit=False)
            met = correspondence_metrics(al.feature_group, ion)
            rows.append(dict(batch=name, anchors_from=kind, fit=fit, K=m.K, lam=lam, sigma_s=m.sigma,
                             loo_rmse_s=float(np.sqrt(np.nanmean(e ** 2))), loo_median_abs_s=float(np.nanmedian(np.abs(e))),
                             f1=met["f1"], rt_rms_s=rt_dispersion(al.aligned_rt, ion), fit_seconds=t.s))
            models[(kind, fit)] = m
    ref = align(tables)  # current pipeline: feature anchors, BIC, refit on all groups
    met = correspondence_metrics(ref.feature_group, ion)
    e = loo_prediction_error(ref.drift, true)
    rows.append(dict(batch=name, anchors_from="features", fit="pipeline", K=ref.drift.K, lam=1.0,
                     sigma_s=ref.drift.sigma, loo_rmse_s=float(np.sqrt(np.nanmean(e ** 2))),
                     loo_median_abs_s=float(np.nanmedian(np.abs(e))), f1=met["f1"],
                     rt_rms_s=rt_dispersion(ref.aligned_rt, ion), fit_seconds=ref.timings["total"]))
    tm = pd.DataFrame(d["timings"]).sum()
    timing = dict(batch=name, runs=len(tables), traces_s=tm.traces, scan_graph_s=tm.scan_graph, peaks_s=tm.peaks,
                  precision_k1=tab.attrs["precision"].get("k1"), precision_k2=tab.attrs["precision"].get("k2"))
    return rows, quality, timing, dict(Ar=Ar, Af=Af, mem=mem, d=d, models=models)


# ------------------------------------------------------------------------------------------------ real data
def tile_prediction(traces, hw_ppm, mz_bin=5.0, rt_bin=60.0):
    """Hit fraction a wide tile (``mz_bin`` Da x ``rt_bin`` s) would predict for each ROI's window."""
    run = traces.run
    sc = run.scan_of_point
    mzb = np.floor(run.mz / mz_bin).astype(np.int64)
    rtb = np.floor(run.rt[sc] / rt_bin).astype(np.int64)
    key = mzb * 100_000 + rtb
    uk, cnt = np.unique(key, return_counts=True)
    scans_per = np.bincount(np.floor(run.rt / rt_bin).astype(np.int64))
    t = traces.table
    apex_rt = t["apex_rt"].to_numpy()
    q = np.floor(t["mz"].to_numpy() / mz_bin).astype(np.int64) * 100_000 + np.floor(apex_rt / rt_bin).astype(np.int64)
    i = np.clip(np.searchsorted(uk, q), 0, uk.size - 1)
    n = np.where(uk[i] == q, cnt[i], 0)
    dens = n / (mz_bin * np.maximum(scans_per[np.floor(apex_rt / rt_bin).astype(np.int64)], 1))  # per Da per scan
    return -np.expm1(-dens * 2 * hw_ppm * 1e-6 * t["mz"].to_numpy())


def real_files(dataset):
    files = [f for f in real_runs(dataset) if not f.name.lower().startswith("blank")]
    return files[:5] if dataset == "metabolomics" else files


SEG_COLS = ["anchor", "partner", "target", "kind", "z", "steps", "n", "s0", "s1", "posterior", "ratio_med", "node_a"]


def process_real(dataset, f, first=False):
    """Both paths for one real run.  Tracking, the scan graph and peaks are cached per run (trace arrays, the
    confident segments, block compounds and time prior); relations are recomputed from that cache, so the relation
    step can change without re-tracking millions of centroids."""
    path = CACHE / f"xrun-real-{dataset}-{f.stem}-stage1.pkl"
    run = MSRun.load(f)
    if path.exists():
        c = pickle.loads(path.read_bytes())
    else:
        with Timer() as t_tr:
            tr = detect_traces(run)
        with Timer() as t_g:
            g = scan_graph(tr, peptides=dataset != "metabolomics")
        with Timer() as t_p:
            feats = detect_peaks(tr, trace_baselines(tr))
        seg = g.segments
        seg = seg[(seg.posterior >= 0.9) & (seg.partner >= 0)][SEG_COLS].reset_index(drop=True)
        c = dict(label=tr.label.astype(np.int32), order=tr.order.astype(np.int32), offsets=tr.offsets,
                 table=tr.table, precision=tr.precision, params=tr.params, seg=seg,
                 compound=g.blocks["compound"].to_numpy(), prior=g.prior, feats=feats,
                 timings=dict(traces=t_tr.s, scan_graph=t_g.s, peaks=t_p.s), roi=None)
        if first:
            nz = roi_noise(tr)
            c["roi"] = nz.assign(tile=tile_prediction(tr, nz["halfwidth_ppm"].to_numpy()), mz=tr.table["mz"].to_numpy(),
                                 apex_intensity=tr.table["apex_intensity"].to_numpy(),
                                 n_points=tr.table["n_points"].to_numpy(), floor_hw_ppm=floor_halfwidth_ppm(tr))
        path.write_bytes(pickle.dumps(c))
        print(dataset, f.stem, {k: round(v, 1) for k, v in c["timings"].items()}, len(seg), "confident segments")
    tr = Traces(run, c["label"].astype(np.int64), c["order"].astype(np.int64), c["offsets"], c["table"],
                c["precision"], c["params"])
    g = SimpleNamespace(segments=c["seg"], blocks=pd.DataFrame(dict(compound=c["compound"])), prior=c["prior"])
    with Timer() as t_r:
        rel = relation_table(g, tr)
        rn = relation_noise(rel, tr)
    tm = dict(c["timings"], relations=t_r.s)
    return dict(rel=rel, noise=rn, feats=c["feats"], timings=tm, name=f.stem, n_points=run.n_points, roi=c["roi"])


def real_batch(dataset):
    files = real_files(dataset)
    if len(files) < 3:
        return None
    rs = [process_real(dataset, f, first=(k == 0)) for k, f in enumerate(files)]
    return dict(rels=[r["rel"] for r in rs], noise=[r["noise"] for r in rs], tables=[r["feats"] for r in rs],
                timings=[r["timings"] for r in rs], names=[r["name"] for r in rs], roi=rs[0]["roi"])


def paired_cv(sets: dict, fits: dict, block: float = 10.0, n_folds: int = 4, seed: int = 0) -> pd.DataFrame:
    """Every anchor set's drift model scored on the same held-out RT blocks of every set.

    Blocks are cut on one absolute RT grid for all sets and assigned to the same folds, so for a given held-out
    block both models have seen neither set's anchors there; each held-out anchor's RT in every run is predicted
    from its other runs.  Rows sharing ``predicted`` are measured on the same anchors and compare directly.
    """
    lo = min(float(np.nanmin(A.rt)) for A in sets.values())
    nb = int(max(float(np.nanmax(A.rt)) for A in sets.values()) - lo) // int(block) + 2
    perm = np.random.default_rng(seed).permutation(nb)
    fold = {k: perm[np.floor((np.nanmedian(A.rt, 1) - lo) / block).astype(int)] % n_folds for k, A in sets.items()}
    err = {(src, tgt): [] for src in sets for tgt in sets}
    for f in range(n_folds):
        models = {}
        for src, A in sets.items():
            K, lam = fits[src]
            tr = fold[src] != f
            models[src] = fit_drift(A.rt[tr], K, t_range=support_range(A.rt), n_basis=12, lam=lam, S=A.se[tr])
        for tgt, A in sets.items():
            te = fold[tgt] == f
            for src, m in models.items():
                e = loo_prediction_error(m, A.rt[te]).ravel()
                err[(src, tgt)].append(e[np.isfinite(e)])
    rows = []
    for (src, tgt), es in err.items():
        e = np.concatenate(es)
        rows.append(dict(predicted=tgt, model_from=src, K=fits[src][0], lam=fits[src][1], n=int(e.size),
                         rmse_s=float(np.sqrt(np.mean(np.minimum(e * e, 25.0)))), median_abs_s=float(np.median(np.abs(e))),
                         gt_2s=float(np.mean(np.abs(e) > 2))))
    return pd.DataFrame(rows)


def cross_prediction(dataset, r):
    """No truth: feature and relation anchors (the relation matching window from the measured drift), each with K
    and lambda by blocked CV, then both models scored on the same held-out blocks of both anchor sets."""
    Af = find_anchors(r["tables"])
    with Timer() as tw:
        win = drift_window(r["rels"], r["noise"])
        Ar, tab, _ = relation_anchors(r["rels"], noise=r["noise"], rt_window=win)
    fitted = {}
    for kind, A in (("features", Af), ("relations", Ar)):
        with Timer() as t:
            m, cv = cv_drift(A.rt, A.se, n_basis=12)
        fitted[kind] = (m, cv, t.s)
    fits = {k: (v[0].K, float(v[1].sort_values("rmse").lam.iat[0])) for k, v in fitted.items()}
    df = paired_cv({"features": Af, "relations": Ar}, fits).assign(dataset=REAL_LABEL[dataset])
    tm = pd.DataFrame(r["timings"]).sum()
    info = dict(runs=len(r["names"]), names=r["names"], traces_s=tm.traces, scan_graph_s=tm.scan_graph,
                peaks_s=tm.peaks, relation_step_s=tm.get("relations", np.nan),
                relations=int(sum(len(x) for x in r["rels"])), precision=tab.attrs["precision"], rt_window=win,
                relation_anchors=int(Ar.rt.shape[0]), feature_anchors=int(Af.rt.shape[0]))
    return df, info, dict(Af=Af, Ar=Ar, fitted=fitted, fits=fits)


SPECIES = dict(names=["human", "yeast", "E. coli"], offsets=[0.0, 1.0, -2.0])


def log2_ratio(X, cond, min_per=2):
    """Mean log2 in A minus mean log2 in B for rows with at least ``min_per`` values in each condition."""
    L = np.log2(np.where(X > 0, X, np.nan))
    a, b = L[:, cond == "A"], L[:, cond == "B"]
    ok = (np.isfinite(a).sum(1) >= min_per) & (np.isfinite(b).sum(1) >= min_per)
    with np.errstate(invalid="ignore"):
        return np.nanmean(a, 1)[ok] - np.nanmean(b, 1)[ok], ok


SNR_BINS = ((0, np.inf, "all"), (3, 10, "3-10"), (10, 30, "10-30"), (30, np.inf, "30+"))


def design_test(dataset, r, bundle):
    """Mixed-species design (ProteoBench M7): A/B ratios of cross-run groups should sit at the species' known
    log2 fold changes; the mixture's widths and background weight score grouping + quantification together,
    stratified by the groups' median SNR (relations: height over the noise in the ions' own m/z windows; features:
    peak SNR)."""
    cond = np.array([n.split("_")[0] for n in r["names"]])
    if not ((cond == "A").sum() >= 2 and (cond == "B").sum() >= 2):
        return None
    fitted, fits = bundle["fitted"], bundle["fits"]
    mr = fitted["relations"][0]
    rows, ratios = [], {}
    with Timer() as t:
        grp = relation_groups(r["rels"], mr)
    x, ok = log2_ratio(grp.height, cond)
    snr = np.nanmedian(grp.snr, 1)[ok]
    lab = "relation groups (scan graph, no peaks)"
    ratios[lab] = (x, snr)
    for lo, hi, b in SNR_BINS:
        sel = (snr >= lo) & (snr < hi)
        rows.append(dict(groups=lab, drift="relation anchors, CV", snr=b, total=len(grp.table), seconds=t.s,
                         **_mix(x[sel])))
    for drift, extra in (("feature anchors, BIC (pipeline)", {}),
                         ("relation anchors, CV", dict(anchors=bundle["Ar"], K=fits["relations"][0],
                                                       lam=fits["relations"][1], refit=False))):
        with Timer() as t:
            al = align(r["tables"], **extra)
        x, ok = log2_ratio(al.intensity, cond)
        FI = al.feature_index
        S = np.full(FI.shape, np.nan)
        for j, tb in enumerate(r["tables"]):
            m = FI[:, j] >= 0
            S[m, j] = tb["snr"].to_numpy()[FI[m, j]]
        snr = np.nanmedian(S, 1)[ok]
        lab = "feature groups (picked peaks)"
        if drift.startswith("feature"):
            ratios[lab] = (x, snr)
        for lo, hi, b in SNR_BINS:
            sel = (snr >= lo) & (snr < hi)
            rows.append(dict(groups=lab, drift=drift, snr=b, total=len(al.groups), seconds=t.s, **_mix(x[sel])))
    return pd.DataFrame(rows), ratios, grp


def _mix(x):
    m = design_mixture(x, SPECIES["offsets"], SPECIES["names"])
    f = design_mixture(x, SPECIES["offsets"], SPECIES["names"], free_means=True)
    out = dict(quantified=m.get("n", 0), background=m.get("background", np.nan), shift=m.get("shift", np.nan))
    for nm in SPECIES["names"]:
        out[f"sd_{nm}"] = m.get(f"sd_{nm}", np.nan)
        out[f"weight_{nm}"] = m.get(f"weight_{nm}", np.nan)
        out[f"free_mean_{nm}"] = f.get(f"mean_{nm}", np.nan) - m.get("shift", 0.0)
    return out


# ------------------------------------------------------------------------------------------------ figures
def loo_figure(R):
    fig = go.Figure()
    for (kind, fit), d in R[R.fit != "pipeline"].groupby(["anchors_from", "fit"]):
        fig.add_trace(go.Bar(x=d.batch, y=d.loo_rmse_s, name=f"{LABELS[kind].split(' (')[0]}, {FITS[fit]}",
                             marker=dict(color=method_color("macha" if kind == "relations" else "global"),
                                         opacity=1.0 if fit == "cv" else 0.45)))
    p = R[R.fit == "pipeline"]
    fig.add_trace(go.Scatter(x=p.batch, y=p.loo_rmse_s, mode="markers", name="current pipeline (features, refit)",
                             marker=dict(symbol="line-ew-open", size=28, line=dict(width=2, color="#888"))))
    fig.update_layout(barmode="group")
    fig.update_yaxes(title="held-out RT prediction error (s, RMSE)", rangemode="tozero")
    return style(fig, 380, title="Drift model quality on true compound RTs (leave one run out)")


def precision_figure(mem, rels, RT, dc):
    """Relation-anchor centre error against SNR, with the learned standard error."""
    tc = np.array([rels[j]["true_compound"].to_numpy()[i] for j, i in zip(mem.run, mem.rel)])
    m = mem.assign(tc=tc)
    m = m[(m.tc >= 0)]
    m = m[dc[m.tc.to_numpy()] >= 0]
    res = m.rt.to_numpy() - RT[m.tc.to_numpy(), m.run.to_numpy()]
    m = m.assign(res=res)
    m["res"] = m.res - m.groupby("anchor_id").res.transform("mean")
    bins = np.logspace(1, 4, 13)
    m["bin"] = np.digitize(m.snr, bins)
    g = m.groupby("bin").agg(snr=("snr", "median"), err=("res", lambda x: 1.4826 * np.median(np.abs(x))),
                             se=("rt_se", "median"), n=("res", "size"))
    g = g[g.n >= 20]
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=g.snr, y=g.err, mode="markers+lines", name="observed centre error (robust sd)",
                             line=dict(color=method_color("macha"), width=2)))
    fig.add_trace(go.Scatter(x=g.snr, y=g.se, mode="lines", name="learned standard error (median)",
                             line=dict(color=method_color("global"), width=2, dash="dash")))
    fig.update_xaxes(title="relation SNR (height over the noise in its ions' own m/z windows)", type="log")
    fig.update_yaxes(title="seconds", rangemode="tozero")
    return style(fig, 340, title="Relation centres: error vs SNR, three-trends-strong")


def ratio_figure(ratios, D, min_snr=10.0):
    """log2 A/B of cross-run groups (median SNR >= ``min_snr``) with the expected species positions (the relation
    groups' fitted common shift + the known offsets)."""
    fig = go.Figure()
    colors = [method_color("macha"), method_color("global"), method_color("greedy")]
    edges = np.arange(-5, 4.0001, 0.05)
    for (lab, (x, snr)), c in zip(ratios.items(), colors):
        x = x[np.isfinite(x) & (snr >= min_snr)]
        h, _ = np.histogram(x, edges, density=True)  # binned here: the page should not carry every group
        fig.add_trace(go.Bar(x=0.5 * (edges[1:] + edges[:-1]), y=h, width=0.05, name=f"{lab} ({x.size:,})",
                             opacity=0.55, marker=dict(color=c, line=dict(width=0))))
    shift = float(D[(D.snr == "all")]["shift"].iloc[0])
    for nm, o in zip(SPECIES["names"], SPECIES["offsets"]):
        fig.add_vline(x=shift + o, line=dict(color="#888", dash="dot"), annotation_text=nm,
                      annotation_position="top")
    fig.update_layout(barmode="overlay")
    fig.update_xaxes(title="log2 (mean A / mean B)", range=[-4.5, 3])
    fig.update_yaxes(title="density")
    return style(fig, 380, title=f"ProteoBench M7: A/B ratios of cross-run groups, median SNR ≥ {min_snr:g}")


def noise_figure(roi):
    x = np.log10(np.maximum(roi["tile"].to_numpy(), 1e-4))
    y = np.log10(np.maximum(roi["hit_frac"].to_numpy(), 1e-4))
    edges = np.linspace(-4.05, 0.05, 61)
    H, _, _ = np.histogram2d(x, y, [edges, edges])
    mid = 0.5 * (edges[1:] + edges[:-1])
    fig = go.Figure(go.Heatmap(x=mid, y=mid, z=np.log10(1 + H.T), colorscale="Blues", showscale=False,
                               hovertemplate="tile %{x:.2f}, ROI %{y:.2f}: log10(1+n) %{z:.2f}<extra></extra>"))
    fig.add_trace(go.Scatter(x=[-4, 0], y=[-4, 0], mode="lines", line=dict(color="#888", dash="dot"),
                             name="tile = ROI"))
    fig.update_xaxes(title="log10 hit fraction predicted by a 5 Da x 60 s tile")
    fig.update_yaxes(title="log10 hit fraction in the ROI's own window")
    return style(fig, 380, legend=False, title="Interference per ROI vs a wide tile (0 hits shown at 1e-4; colour log count)")


def _m7_tiles(real):
    v = real.get("proteomics")
    if v is None or v["design"] is None:
        return []
    D = v["design"][0]
    hi = D[D.snr == "30+"].set_index("groups")
    rg, fg = hi.loc["relation groups (scan graph, no peaks)"], hi[hi.drift.str.startswith("feature anchors")].loc[
        "feature groups (picked peaks)"]
    return [("M7 yeast mode width, SNR ≥ 30", f"{rg['sd_yeast']:.2f} vs {fg['sd_yeast']:.2f}",
             f"log2 A/B, relation groups (no peaks) vs feature groups; yeast at {rg['free_mean_yeast']:+.2f} "
             f"vs {fg['free_mean_yeast']:+.2f} (design +1)")]


def compute():
    """Every number and figure input of the page, cached so the page can be rebuilt without recomputing."""
    rows, quality, timing, prec = [], [], [], None
    for name in BATCHES:
        r, q, t, k = evaluate(name)
        rows += r
        quality += q
        timing.append(t)
        if name == "three-trends-strong":
            d = k["d"]
            prec = (k["mem"], [x[["true_compound"]] for x in d["rels"]], d["rt"], d["drift_class"])
        print(name, "done")
    R, Q, T = pd.DataFrame(rows), pd.DataFrame(quality), pd.DataFrame(timing)
    real = {}
    for ds in ("metabolomics", "proteomics"):
        rb = real_batch(ds)
        if rb is None:
            continue
        df, info, bundle = cross_prediction(ds, rb)
        design = design_test(ds, rb, bundle)
        if design is not None:
            D, ratios, grp = design
            design = (D, ratios, dict(groups=len(grp.table), gap=grp.table.attrs["gap"]))
        real[ds] = dict(df=df, info=info, design=design, roi=rb["roi"])
        print(ds, "done")
    res = dict(R=R, Q=Q, T=T, real=real, prec=prec)
    (CACHE / "xrun-results.pkl").write_bytes(pickle.dumps(res))
    return res


def main(rebuild_page_only: bool = False):
    path = CACHE / "xrun-results.pkl"
    res = pickle.loads(path.read_bytes()) if rebuild_page_only and path.exists() else compute()
    R, Q, T, real = res["R"], res["Q"], res["T"], res["real"]
    save_json("xrun", dict(drift=R.to_dict("records"), anchors=Q.to_dict("records"), timing=T.to_dict("records"),
                           real={ds: dict(paired_cv=v["df"].to_dict("records"), info=v["info"],
                                          design=None if v["design"] is None else v["design"][0].to_dict("records"))
                                 for ds, v in real.items()}))

    page = Page("Cross-run graph: alignment before peaks",
                "First step towards grouping across runs at the scan-graph stage: the drift-trend mixture is "
                "learned from relations of the scan-level graph, matched across runs by rule and the m/z of both "
                "ions, without picking a single peak.", NAV)
    cv = R[R.fit == "cv"].pivot(index="batch", columns="anchors_from", values="loo_rmse_s")
    qq = Q.pivot(index="batch", columns="anchors_from", values="centre_error_s")
    page.tiles([
        ("drift error, relation vs feature anchors", f"{cv['relations'].mean():.2f} / {cv['features'].mean():.2f} s",
         "mean over batches of held-out RT error, K and λ by cross-validation"),
        ("anchor centre error", f"{qq['relations'].median():.3f} vs {qq['features'].median():.3f} s",
         "relation centre vs feature apex, against truth (median over batches)"),
    ] + _m7_tiles(real))
    page.section("What it does", """
<p><b>Relations as anchors.</b> A feature anchor must be the only feature within 6 ppm over a 60 s window in every
run; that fails wherever m/z is crowded, and it needs picked peaks. A confident scan-graph relation is keyed by its
rule and the m/z of <i>both</i> ions, so it is far more specific, and it exists before peak detection. Relations
are matched across runs by chain clustering on rule, anchor m/z, partner m/z and RT, at most one per run, in at
least 60 % of runs. A relation is kept as an anchor only if:</p>
<ul class="prose">
<li><b>it elutes:</b> at least 60 % of its ions' signal within ±3σ lies within ±σ of its centre (σ = 2 elution sd).
A relation between chemical-background ions has no localised elution and an arbitrary centre.</li>
<li><b>it has no competitor:</b> no other ion within 6 ppm of either of its ions reaches 20 % of its height in the
matching window. That other ion can be another trace, or the same trace beyond a valley, as with an isomer joined
into one trace. Otherwise, when one of an isomer pair is missed in some runs, the other is matched in its place.</li>
<li><b>chance is negligible:</b> its two ions' noise-floor windows are quiet enough that a chance relation of three
consecutive scans is implausible.</li>
<li><b>it is unique</b> in its run (no other confident relation with the same key in the window) and has SNR ≥ 5.</li></ul>
<p>The matching window is not fixed. A first pass at ±60 s measures the drift (2.5 × the 95th percentile of anchor
deviations, at least 10 s), and the filters are evaluated for that window. In dense Astral proteomics MS1 a
stronger ion within 6 ppm is almost always present somewhere in ±60 s.</p>
<p><b>Centres.</b> The relation says which two ions and roughly when. Its centre comes from both ions' full traces in
a fixed window of ±2 elution sd (the run's own, from the scan graph), weighted by intensity² and re-centred. The
scans where both ions were seen are not used directly. A weak partner is seen only near the apex, so that window
is truncated by an amount that depends on abundance, and its centre has twice the error.</p>
<p><b>Precision, learned.</b> SNR is height over the noise level in the ions' own m/z windows. The standard error is
√((k₁·sd/SNR)² + (k₂·sd)²). k₁, k₂ are fitted on <i>double differences</i>: for two relations of one compound with
no shared ion, (t<sub>a</sub> − t<sub>b</sub>) in run j minus the same in run j′. Drift, the compound's own RT error
and each relation's bias cancel, so only measurement noise is left. A within-run estimate overstates the error
about twofold, because truncation differs between relations.</p>
<p><b>Duplicates.</b> Relations sharing an ion measure the same elution: an ion with an isotope, an adduct and a dimer
partner anchors three times. Keeping one anchor per such group, or averaging the group, was tried and is worse.
Averaging lets one wrong link merge two compounds' centres. Keeping one leaves some drift classes too sparse in
places, and their warps wander there. With K and λ chosen by cross-validation instead of BIC, duplicates no longer
inflate the number of classes, so all of them are kept.</p>
<p><b>Choosing K and λ.</b> The drift penalty is not scale-free. The data term grows with the number and precision of
anchors, so λ = 1 suits one anchor set and over-smooths another. BIC also over-counts anchors that share a
compound's RT error. Both are therefore chosen by blocked cross-validation: anchors are assigned to 10 s RT blocks,
so co-eluting siblings are held out together. Each held-out anchor's RT in every run is predicted from its other
runs.</p>""")
    page.section("Noise per ROI", """
<p>Whether a match is a coincidence depends on how busy the window it was accepted in is. Noise is therefore
measured for every ROI in its own m/z window: the noise-floor acceptance width of the relation search, or the ROI's
own scatter if wider, over its own scans. Centroids of the ROI itself are excluded. The map gives the per-scan
chance of a spurious partner (hit fraction) and the local noise level (log-intensity mean and sd of interfering
centroids). Relation SNR and the chance-relation filter use it. A wide tile averages a window next to a persistent
background ion with an empty one; on real data the two disagree by orders of magnitude in both directions.</p>""")
    for ds, v in real.items():
        roi = v["roi"]
        if roi is None:
            continue
        hf, tl = roi["hit_frac"].to_numpy(), roi["tile"].to_numpy()
        both = (hf > 0) & (tl > 0)
        ratio = hf[both] / tl[both]
        page.figure(noise_figure(roi), f"{REAL_LABEL[ds]}, first run, {len(roi):,} ROIs. Of ROIs with interference, "
                    f"{np.mean(ratio > 3):.0%} see more than 3× and {np.mean(ratio < 1 / 3):.0%} less than ⅓ of the "
                    f"tile's prediction; {np.mean(hf == 0):.0%} have none at all.")
    page.section("Against ground truth", "<p>Alignment batches as on the <a href=\"align.html\">alignment page</a>, "
                 "same seeds. <i>Held-out RT error</i>: for every eluting compound and run, its true RT is predicted "
                 "from its true RTs in the other runs with the fitted warps and class posterior. This scores the "
                 "warps themselves, independent of which anchors they were learned from; the floor is the "
                 "compound-level RT jitter (0.25 s · √(1 + 1/(J−1)) ≈ 0.27 s). F1 and aligned RT spread come from "
                 "correspondence of picked features with that drift model and no refit, so only the anchors "
                 "differ.</p>")
    page.figure(loo_figure(R))
    page.table(R.assign(anchors_from=R.anchors_from.map(LABELS), fit=R.fit.map({**FITS, "pipeline": "pipeline"})),
               "Drift fits per anchor source. <i>pipeline</i> is the current default (feature anchors, BIC, refit on "
               "all matched groups).", better=dict(loo_rmse_s="min", f1="max", rt_rms_s="min"), group_col="batch")
    page.table(Q.assign(anchors_from=Q.anchors_from.map(LABELS)),
               "Anchors against truth. <i>centre_error_s</i>: robust sd of anchor RT minus true compound RT after "
               "removing each anchor's mean, so the compound's own RT error cancels and only measurement error is "
               "left. <i>se_calibration</i>: robust sd of that error over the stated standard error (1 = calibrated). "
               "<i>compounds</i>: eluting compounds covered.",
               better=dict(centre_error_s="min", gross_gt_1s="min", purity="max", compounds="max"), group_col="batch")
    if res["prec"] is not None:
        page.figure(precision_figure(*res["prec"]))
    page.table(T, "Summed per-batch time (s). The scan-graph path replaces baselines and peaks; precision k₁, k₂ "
               "from cross-run double differences.")
    for ds, v in real.items():
        df, info = v["df"], v["info"]
        design = {"metabolomics": "Sa04 ×3 injections, Sa09, Sa17", "proteomics": "conditions A and B ×3"}[ds]
        page.section(f"Real runs: {REAL_LABEL[ds]}", f"""
<p>{info['runs']} runs ({design}), no truth. Both drift models (K and λ chosen by blocked CV on their own anchors)
are refitted on the same RT-block folds and scored on the held-out blocks of <i>both</i> anchor sets, so rows that
share <i>predicted</i> are measured on the same anchors and neither model has seen anchors in a held-out block.
The relation matching window, ±{info['rt_window']:.0f} s, is set from the drift measured in a first pass.
{info['relation_anchors']:,} relation anchors from {info['relations']:,} confident relations,
{info['feature_anchors']:,} feature anchors. Learned centre precision k₁ = {info['precision'].get('k1', float('nan')):.2f},
k₂ = {info['precision'].get('k2', float('nan')):.3f}. Time for all runs: traces {info['traces_s']:.0f} s, scan graph
{info['scan_graph_s']:.0f} s, relations {info['relation_step_s']:.0f} s, baselines + peaks {info['peaks_s']:.0f} s.</p>""")
        page.table(df.drop(columns=["dataset"]), better=dict(rmse_s="min", median_abs_s="min", gt_2s="min"),
                   group_col="predicted")
        if v["design"] is not None:
            D, ratios, grp = v["design"]
            page.section("Mixed-species design (ProteoBench M7)", f"""
<p>Conditions A and B mix human, yeast and E. coli digests at A/B = 1, 2 and 0.25 (log2 0, +1, −2; ProteoBench's
scoring configuration). With no identifications, the A/B ratio of every cross-run group should still fall on one of
three modes. A mixture with a common shift, one Gaussian per species at the known spacing and a uniform background
scores grouping and quantification together: narrow modes mean precise ratios, background weight means groups that
fit no species. <i>free_mean</i> lets each mode move on its own, which checks the design.</p>
<p><b>Relation groups</b> are cross-run groups formed at the scan-graph stage, with no peak detection.
{grp['groups']:,} relation groups are formed by rule and the m/z of both ions in reference time (drift from relation
anchors; grouping gap {grp['gap']:.1f} s). Each is quantified per run by the sum of its two ions' maxima
in the fixed window. <b>Feature groups</b> are picked peaks matched by the current pipeline, and again with the relation
drift. Groups need two values in each condition and are stratified by median SNR.</p>""")
            page.figure(ratio_figure(ratios, D))
            cols = ["groups", "drift", "snr", "quantified", "background", "shift", "sd_human", "weight_human",
                    "sd_yeast", "weight_yeast", "free_mean_yeast", "sd_E. coli", "weight_E. coli", "free_mean_E. coli"]
            page.table(D[cols], better={"sd_human": "min", "sd_yeast": "min", "background": "min"}, group_col="snr")
    page.section("Next steps", """<ul class="prose">
<li><b>Correspondence of graph nodes.</b> Here the drift model is used to match picked features. Next comes matching
scan-graph nodes (trace blocks with centres from the same fixed-window estimator) across runs under this model,
neighbourhood-aware, so isomers at one m/z do not swap.</li>
<li><b>Cross-run edges with dropout.</b> Lift relations to node-group pairs. Re-query runs where a relation is
missing at its predicted place in raw centroids, and score presence and absence against detectability, using the
per-ROI noise level. The chance rate is the RT-shifted decoy hit rate in the ions' own windows.</li>
<li><b>Sample design.</b> A runs × sources matrix (replicates, spikes, matrix, blank, dilution) for expected presence
and intensity ratios.</li>
<li><b>Vector-coherent partition.</b> GAEC with one centre per run (Ward increment summed over runs, missing runs at
zero precision), so ions that move together across runs group even without a mass rule.</li></ul>""")
    page.write(REPORTS / "xrun.html")
    print(R.round(3).to_string())
    print(Q.round(3).to_string())
    print(T.round(1).to_string())
    for ds, v in real.items():
        print(v["df"].round(3).to_string())
        print(v["info"])
        if v["design"] is not None:
            print(v["design"][0].round(3).to_string())


if __name__ == "__main__":
    import sys

    main(rebuild_page_only="--page" in sys.argv)
