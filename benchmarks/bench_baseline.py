"""Benchmark + explorer: censoring-aware baselines."""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from common import NAV, REAL_LABEL, REPORTS, Timer, real_run, save_json
from eic_sim import NOISE, THRESH, make_eic
from macha.baseline import als_baseline, censored_arpls, rolling_min_baseline, trace_baselines
from macha.report import LIGHT, Page, method_color, style
from macha.traces import detect_traces

KINDS = ["none", "constant", "drift", "step"]
LABELS = {"macha": "censored arPLS (macha)", "als_zero": "ALS, missing = 0", "als_locf": "ALS, last value "
          "carried forward (original macha)", "rolling_min": "rolling minimum + smoothing"}


def estimate(method, y):
    obs = np.isfinite(y)
    if method == "macha":
        return censored_arpls(np.where(obs, y, 0.0), obs, np.full(y.size, THRESH), 1e4, 10)
    if method == "als_zero":
        return als_baseline(y, 1e4, 0.01, fill="zero")
    if method == "als_locf":
        return als_baseline(y, 1e4, 0.01, fill="locf")
    return rolling_min_baseline(y, 41, 21)


def peak_mask(e, n):
    m = np.zeros(n, bool)
    for p in e["peaks"]:
        a = int(p["apex"] - 4 * p["sigma"])
        b = int(p["apex"] + 4 * p["sigma"] + 6 * p["tau"])
        m[max(a, 0): min(b, n)] = True
    return m


def run_benchmark(n_rep=150, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    examples = {}
    for kind in KINDS:
        for r in range(n_rep):
            e = make_eic(rng, kind=kind)
            pm = peak_mask(e, e["x"].size)
            for m in LABELS:
                z = estimate(m, e["y"])
                err = z - e["base"]
                h_err = [((e["signal"][int(round(p["apex"]))] - z[int(round(p["apex"]))]) - p["height"]) / p["height"]
                         for p in e["peaks"]]
                rows.append(dict(background=kind, method=m, rmse_noise=np.sqrt(np.mean(err**2)) / NOISE,
                                 rmse_under_peaks=np.sqrt(np.mean(err[pm] ** 2)) / NOISE if pm.any() else np.nan,
                                 height_bias_pct=100 * float(np.median(h_err))))
            if r == 0:
                examples[kind] = e
    df = pd.DataFrame(rows)
    summ = df.groupby(["background", "method"]).agg(rmse_noise=("rmse_noise", "median"),
                                                    rmse_under_peaks=("rmse_under_peaks", "median"),
                                                    height_bias_pct=("height_bias_pct", "median")).reset_index()
    # runtime: batch of traces
    ys = [make_eic(rng, kind="drift")["y"] for _ in range(1000)]
    times = {}
    for m in LABELS:
        with Timer() as t:
            for y in ys:
                estimate(m, y)
        times[m] = t.s
    summ["ms_per_1000_traces"] = summ["method"].map(lambda m: times[m] * 1000)
    return summ, examples


def example_figure(e, kind):
    fig = go.Figure()
    x = e["x"]
    fig.add_trace(go.Scatter(x=x, y=e["y"], mode="markers", name="observed centroids",
                             marker=dict(size=5, color=LIGHT["muted"])))
    fig.add_trace(go.Scatter(x=x, y=e["base"], mode="lines", name="true baseline",
                             line=dict(width=2, color=LIGHT["ink"])))
    for m in LABELS:
        fig.add_trace(go.Scatter(x=x, y=estimate(m, e["y"]), mode="lines", name=LABELS[m],
                                 line=dict(width=2, color=method_color(m))))
    fig.add_hline(y=THRESH, line=dict(width=1, color=LIGHT["axis"]))
    fig.update_xaxes(title="scan")
    fig.update_yaxes(title="intensity", range=[0, np.nanpercentile(e["y"], 92) * 1.4 + 5 * THRESH])
    return style(fig, 340, title=f"Background: {kind} (horizontal rule = detection threshold)")


def main():
    summ, examples = run_benchmark()
    save_json("baseline", summ.to_dict("records"))
    page = Page("Baseline & noise", "Orbitrap centroid data never contains zeros: anything below the noise threshold is "
                "simply absent. The baseline treats an absent point as an upper bound (≤ detection floor), not as a "
                "zero and not as the last observed value.", NAV)
    piv = summ.pivot(index="background", columns="method", values="rmse_noise")
    page.tiles([
        (f"baseline RMSE, {k} background", f"{piv.loc[k, 'macha']:.2f}σ",
         f"ALS-zero {piv.loc[k, 'als_zero']:.2f}σ · LOCF {piv.loc[k, 'als_locf']:.2f}σ") for k in ("drift", "step", "constant")])
    page.section("What it does", """
<p>Asymmetric Whittaker smoothing (arPLS weights) with a second-difference penalty, solved in O(n) by a banded
LDLᵀ factorisation compiled with numba and run over every trace in one call. Observed points get arPLS weights, so
peaks are down-weighted. <b>Censored points</b> (scans where the trace has no centroid) are handled as one-sided
information: they gently pull the baseline toward a low value while it sits below the detection floor, and become
hard constraints at the floor if the baseline tries to rise above it, since a baseline above the floor would have been
observed. The detection floor itself is estimated per scan block and m/z block as a low quantile of reported
intensities. The original macha filled gaps by carrying the last value forward (LOCF) before an IRLS baseline, which
props the baseline up after every peak.</p>""")
    page.section("Accuracy on simulated traces", "<p>150 traces per background type, each with 1-5 EMG peaks "
                 "(SNR 5-200), multiplicative and additive noise and censoring below 3σ. RMSE in units of the noise σ "
                 "(median over traces); height bias is the median relative error of baseline-subtracted apex "
                 "heights.</p>")
    tab = summ.copy()
    tab["method"] = tab["method"].map(LABELS)
    page.table(tab, better=dict(rmse_noise="min", rmse_under_peaks="min", ms_per_1000_traces="min"),
               group_col="background")
    fig = go.Figure()
    for m in LABELS:
        d = summ[summ.method == m].set_index("background").reindex(KINDS)
        fig.add_trace(go.Scatter(x=d.rmse_noise, y=KINDS, mode="markers", name=LABELS[m],
                                 marker=dict(size=12, color=method_color(m), line=dict(width=2, color=LIGHT["surface"]))))
    fig.update_xaxes(title="baseline RMSE (noise σ units, log)", type="log")
    fig.update_yaxes(autorange="reversed")
    page.figure(style(fig, 280, title="Baseline error by background type"))
    page.section("Examples", "<p>One simulated trace per background type. Gray points are the centroids an "
                 "Orbitrap would report; gaps are censored scans.</p>")
    for k in KINDS:
        page.figure(example_figure(examples[k], k))
    run = real_run("proteomics_sc", 0, rt_range=(300, 900))
    if run is not None:
        tr = detect_traces(run)
        with Timer() as t:
            bl = trace_baselines(tr)
        page.section("Real Orbitrap Astral traces", f"<p>{len(tr):,} traces of {REAL_LABEL['proteomics_sc']} "
                     f"(RT 300-900 s) baselined in {t.s:.2f} s. Shown: the longest traces, which are background "
                     "ions with peaks on top - the case the method is built for.</p>")
        longest = np.argsort(-tr.table["n_points"].to_numpy())[:3]
        for i in longest:
            d = bl.trace(int(i))
            rt = run.rt[d["scans"]]
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=rt, y=d["intensity"], mode="markers", name="centroids",
                                     marker=dict(size=4, color=LIGHT["muted"])))
            fig.add_trace(go.Scatter(x=rt, y=d["baseline"], mode="lines", name="baseline",
                                     line=dict(width=2, color=method_color("macha"))))
            fig.add_trace(go.Scatter(x=rt, y=d["floor"], mode="lines", name="detection floor",
                                     line=dict(width=1, color=LIGHT["axis"])))
            fig.update_xaxes(title="retention time (s)")
            fig.update_yaxes(title="intensity", type="log")
            page.figure(style(fig, 300, title=f"trace {i}: m/z {tr.table['mz'].iat[i]:.4f}, "
                                              f"{tr.table['n_points'].iat[i]} centroids, noise {d['noise']:.3g}"))
    page.write(REPORTS / "baseline.html")
    print(summ.round(3).to_string())


if __name__ == "__main__":
    main()
