"""Benchmark + explorer: chromatographic peak detection."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from common import NAV, REAL_LABEL, REPORTS, Timer, real_run, save_json, synthetic
from eic_sim import NOISE, make_eic, make_pair
from macha.pipeline import process_run
from macha.peaks import cwt_peaks, detect_peaks_1d, fit_emg_mixture, sspp_peaks
from macha.report import LIGHT, Page, binned_hist, method_color, style
from macha.synth.truth import feature_truth

LABELS = {"macha": "persistence (macha)", "sspp": "scale-space (original macha)", "cwt": "CWT (centWave-like)"}
SNR_BINS = [3, 5, 10, 20, 50, 100, 300]


def detect(method, e):
    y, b = e["y"], e["base"]
    if method == "macha":
        d = detect_peaks_1d(y, baseline=b, noise=NOISE, floor=np.full(y.size, 3 * NOISE))
        return d["apex_frac"].to_numpy()
    if method == "sspp":
        return sspp_peaks(y, b, NOISE).astype(float)
    return cwt_peaks(y, b, NOISE).astype(float)


def match(found, peaks, tol_scans):
    used = set()
    tp = []
    for p in peaks:
        d = np.abs(np.asarray(found) - p["apex"])
        ok = [k for k in np.argsort(d) if d[k] <= max(tol_scans, 2.355 * p["sigma"] * 0.5) and k not in used]
        if ok:
            used.add(ok[0])
            tp.append((p, float(found[ok[0]] - p["apex"])))
    return tp, len(found) - len(used)


def run_benchmark(n_rep=400, seed=1):
    rng = np.random.default_rng(seed)
    rows, apex = [], []
    times = {m: 0.0 for m in LABELS}
    for r in range(n_rep):
        e = make_eic(rng, kind=["none", "constant", "drift"][r % 3], snr_range=(3, 300))
        for m in LABELS:
            with Timer() as t:
                f = detect(m, e)
            times[m] += t.s
            tp, fp = match(f, e["peaks"], 2.0)
            hit = {id(p) for p, _ in tp}
            for p in e["peaks"]:
                rows.append(dict(method=m, snr=p["snr"], found=id(p) in hit))
            rows.append(dict(method=m, snr=np.nan, found=np.nan, fp=fp))
            apex += [dict(method=m, err=err) for _, err in tp]
    df = pd.DataFrame(rows)
    pk = df[df.snr.notna()].copy()
    pk["bin"] = pd.cut(pk.snr, SNR_BINS)
    rec = pk.groupby(["method", "bin"], observed=True).found.mean().reset_index()
    fp = df[df.snr.isna()].groupby("method").fp.sum()
    n_found = pk.groupby("method").found.sum()
    summ = pd.DataFrame(dict(recall=pk.groupby("method").found.mean(),
                             precision=n_found / (n_found + fp),
                             apex_mae_scans=pd.DataFrame(apex).groupby("method").err.apply(lambda x: np.mean(np.abs(x))),
                             ms_per_trace=pd.Series(times) / n_rep * 1000)).reset_index().rename(columns={"index": "method"})
    summ["f1"] = 2 * summ.precision * summ.recall / (summ.precision + summ.recall)
    # resolution of overlapping pairs
    res = []
    for sep in np.arange(1.0, 6.01, 0.5):
        for _ in range(60):
            e = make_pair(rng, sep)
            for m in LABELS:
                f = detect(m, dict(y=e["y"], base=e["base"]))
                ok = all(np.min(np.abs(np.asarray(f) - a)) <= 2.5 for a in e["apexes"]) if len(f) else False
                res.append(dict(method=m, sep=sep, resolved=ok))
    res = pd.DataFrame(res).groupby(["method", "sep"]).resolved.mean().reset_index()
    return summ, rec, res


def main():
    summ, rec, res = run_benchmark()
    page = Page("Peak detection", "Topological persistence on lightly smoothed traces: every local maximum's prominence "
                "is computed in one union-find sweep, peaks are kept by prominence and height in noise units, and "
                "boundaries are the valleys between kept peaks trimmed to the baseline.", NAV)
    s = summ.set_index("method")
    page.tiles([("F1, simulated traces", f"{s.loc['macha', 'f1']:.3f}",
                 f"scale-space {s.loc['sspp', 'f1']:.3f} · CWT {s.loc['cwt', 'f1']:.3f}"),
                ("apex error", f"{s.loc['macha', 'apex_mae_scans']:.2f} scans",
                 f"scale-space {s.loc['sspp', 'apex_mae_scans']:.2f} · CWT {s.loc['cwt', 'apex_mae_scans']:.2f}"),
                ("time per trace", f"{s.loc['macha', 'ms_per_trace']:.3f} ms",
                 f"CWT {s.loc['cwt', 'ms_per_trace']:.2f} ms")])
    page.section("What it does", """
<p>The original macha seeded peaks with scale-space peak picking (Hamming smoothing at several widths, tracking maxima
across scales) and then iteratively fitted and re-seeded skew-normal mixtures with numerical-gradient optimisers. Here
detection is <b>0-dimensional persistent homology</b> of the smoothed trace: a maximum is born at its height and dies
at the saddle where it merges into a taller peak, so persistence equals topographic prominence, and all maxima are
scored in O(n log n). Peaks need prominence ≥ 3σ and height ≥ 5σ above the baseline (σ = trace noise from the
baseline module) and ≥ 4 observed centroids. Each feature then gets an apex RT (parabolic refinement), height, area,
FWHM, SNR and an inverse-variance m/z with its standard error from the trace module's precision model, which the
annotation and alignment modules use as uncertainty. Overlapping peaks can be refined by an EMG mixture fitted with an
analytic Jacobian (<code>fit_emg_mixture</code>).</p>""")
    page.section("Detection on simulated traces", "<p>400 traces with 1-5 EMG peaks (SNR 3-300, tailing τ/σ 0.1-1) "
                 "on zero, constant or drifting backgrounds, censored below 3σ. A detection matches a true apex "
                 "within max(2 scans, half an FWHM).</p>")
    tab = summ.copy()
    tab["method"] = tab["method"].map(LABELS)
    page.table(tab, better=dict(recall="max", precision="max", f1="max", apex_mae_scans="min", ms_per_trace="min"))
    fig = go.Figure()
    for m in LABELS:
        d = rec[rec.method == m]
        fig.add_trace(go.Scatter(x=[f"{iv.left:g}-{iv.right:g}" for iv in d["bin"]], y=d.found, mode="lines+markers",
                                 name=LABELS[m], line=dict(width=2, color=method_color(m)), marker=dict(size=8)))
    fig.update_xaxes(title="true SNR (height / noise σ)")
    fig.update_yaxes(title="recall", range=[0, 1.02])
    page.figure(style(fig, 320, title="Recall by signal-to-noise"))
    fig = go.Figure()
    for m in LABELS:
        d = res[res.method == m]
        fig.add_trace(go.Scatter(x=d.sep, y=d.resolved, mode="lines+markers", name=LABELS[m],
                                 line=dict(width=2, color=method_color(m)), marker=dict(size=8)))
    fig.update_xaxes(title="apex separation (σ of each peak; chromatographic resolution ≈ separation / 4)")
    fig.update_yaxes(title="both peaks found", range=[0, 1.02])
    page.figure(style(fig, 320, title="Resolving overlapping peaks (second peak at half height)"))
    # EMG refinement example
    rng = np.random.default_rng(7)
    e = make_pair(rng, 2.6, snr=80)
    d = detect_peaks_1d(e["y"], noise=NOISE)
    init = [(h, a, 3.0, 1.0) for h, a in zip(d.height, d.apex_frac)]
    x = e["x"]
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x, y=e["y"], mode="markers", name="centroids", marker=dict(size=6, color=LIGHT["muted"])))
    if len(init) >= 2:
        fit = fit_emg_mixture(x, e["y"], init[:2])
        from macha.synth import emg
        for k, r in fit.iterrows():
            fig.add_trace(go.Scatter(x=x, y=r.amplitude * emg(x, r.mu, r.sigma, r.tau), mode="lines",
                                     name=f"EMG component {k + 1}", line=dict(width=2, color=method_color(
                                         ["macha", "greedy"][k % 2]))))
    fig.update_xaxes(title="scan")
    fig.update_yaxes(title="intensity")
    page.figure(style(fig, 300, title="EMG mixture refinement of an overlapping pair (analytic-Jacobian least squares)"))
    # pipeline-level recall on synthetic runs
    rows = []
    for name in ("syn-metabolomics", "syn-proteomics", "syn-hard"):
        sr = synthetic(name)
        with Timer() as t:
            smp = process_run(sr.run)
        truth = feature_truth(smp.features, smp.traces, sr.ion, sr.ions, sr.compounds)
        ions = sr.ions.merge(sr.compounds[["compound_id", "height", "kind"]], on="compound_id")
        ions = ions[ions.kind != "background"]
        snr = ions.height * ions.rel / sr.scenario.instrument.noise
        found = set(truth.ion[truth.ion >= 0])
        for lo, hi in [(5, 10), (10, 30), (30, 1e9)]:
            sel = ions[(snr >= lo) & (snr < hi)]
            rows.append(dict(dataset=name, true_snr=f"{lo}-{hi if hi < 1e9 else '∞'}",
                             ion_recall=sel.ion_id.isin(found).mean(), n_ions=len(sel)))
        ok = truth.true_rt.notna()
        rows.append(dict(dataset=name, true_snr="all", ion_recall=np.nan, n_ions=len(smp.features),
                         rt_mae_s=float(np.mean(np.abs(smp.features.rt[ok] - truth.true_rt[ok]))),
                         feature_purity=float(truth.purity.mean()), pipeline_s=t.s))
    pr = pd.DataFrame(rows)
    page.section("Whole runs", "<p>Traces → baselines → peaks on the synthetic runs: fraction of true ions found as "
                 "a feature by true SNR, and apex RT error against the EMG mode.</p>")
    page.table(pr, better=dict(ion_recall="max", rt_mae_s="min"))
    real = []
    for ds in ("proteomics_sc", "proteomics", "metabolomics"):
        run = real_run(ds, 0)
        if run is None:
            continue
        with Timer() as t:
            smp = process_run(run)
        f = smp.features
        real.append(dict(dataset=REAL_LABEL[ds], centroids=run.n_points, traces=len(smp.traces), features=len(f),
                         median_fwhm_s=float(f.fwhm_s.median()), median_snr=float(f.snr.median()),
                         median_mz_se_ppm=float((f.mz_se / f.mz * 1e6).median()), seconds=t.s))
        if ds == "proteomics_sc":
            fig = go.Figure()
            fig.add_trace(go.Bar(**binned_hist(np.log10(f.snr), 60, density=False), name="features",
                                 marker=dict(color=method_color("macha"))))
            fig.update_xaxes(title="log10 SNR")
            fig.update_yaxes(title="features")
            hist_fig = style(fig, 280, legend=False, title=f"Feature SNR, {REAL_LABEL[ds]}")
    if real:
        page.section("Real Orbitrap Astral runs", "<p>Full runs through traces, baselines and peaks.</p>")
        page.table(pd.DataFrame(real), better=dict(seconds="min"))
        page.figure(hist_fig)
    save_json("peaks", dict(summary=summ.to_dict("records"), recall=rec.astype(str).to_dict("records"),
                            pairs=res.to_dict("records"), runs=pr.to_dict("records"), real=real))
    page.write(REPORTS / "peaks.html")
    print(summ.round(3).to_string())
    print(pr.round(3).to_string())
    print(pd.DataFrame(real).round(3).to_string())


if __name__ == "__main__":
    main()
