"""Run every module benchmark and build ``reports/index.html``.

    python benchmarks/run_all.py            # all benchmarks, then the index
    python benchmarks/run_all.py --index    # rebuild the index from existing JSON
"""

from __future__ import annotations

import json
import runpy
import sys

import pandas as pd

from common import NAV, REPORTS
from macha.report import Page

MODULES = ["bench_io", "bench_traces", "bench_baseline", "bench_peaks", "bench_annotate", "bench_scan_graph",
           "bench_dia", "bench_align", "bench_xrun", "bench_components"]


def load(name):
    p = REPORTS / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def index():
    page = Page("macha — module explorers",
                "A Python re-implementation of macha as a library of independent LC-MS processing modules. Every "
                "page explains the method, benchmarks it against reference approaches on ground-truth synthetic "
                "data and on public Orbitrap Astral proteomics and metabolomics runs, and lets you explore the "
                "results.", NAV)
    tiles = []
    tr = load("traces")
    if tr:
        s = pd.DataFrame(tr["synthetic"])
        h = s[s.dataset == "syn-hard"].set_index("method")
        tiles.append(("Traces: ARI, hard synthetic", f"{h.loc['macha', 'ari']:.3f}",
                      f"greedy ROI {h.loc['greedy', 'ari']:.3f}, slice {h.loc['slice', 'ari']:.3f}"))
    sg = load("scan_graph")
    if sg:
        r = pd.DataFrame(sg["recall"])
        r = r[r.scenario == "metabolomics, low abundance"].set_index(["method", "snr_bin"]).recall
        tiles.append(("Scan graph: compounds found, SNR 10-30", f"{r[('scan', '10-30')]:.0%}",
                      f"feature-level graph {r[('feature', '10-30')]:.0%} (low-abundance metabolomics)"))
    dia = load("dia")
    if dia:
        d = pd.DataFrame(dia["summary"])
        d = d[d.dataset.str.startswith("Astral")].set_index("approach")
        if len(d):
            best = d.peptides_1pct.idxmax()
            tiles.append(("DIA 250 pg: peptides at 1 % FDR", f"{int(d.peptides_1pct.max())}",
                          f"{best}; co-elution only {int(d.loc['co-elution only', 'peptides_1pct'])}"))
    bl = load("baseline")
    if bl:
        b = pd.DataFrame(bl)
        d = b[b.background == "drift"].set_index("method")
        tiles.append(("Baseline RMSE, drifting background", f"{d.loc['macha', 'rmse_noise']:.2f}σ",
                      f"ALS with last-value fill {d.loc['als_locf', 'rmse_noise']:.2f}σ"))
    pk = load("peaks")
    if pk:
        p = pd.DataFrame(pk["summary"]).set_index("method")
        tiles.append(("Peaks: F1", f"{p.loc['macha', 'f1']:.3f}", f"scale-space {p.loc['sspp', 'f1']:.3f}, "
                                                                   f"CWT {p.loc['cwt', 'f1']:.3f}"))
    an = load("annotate")
    if an:
        a = pd.DataFrame(an["synthetic"])
        m = a[a.method == "macha"]
        tiles.append(("Annotation: edge precision", f"{m.edge_precision.min():.3f}+",
                      f"compound ARI {m.group_ari.min():.2f}–{m.group_ari.max():.2f}"))
        if an["real"]:
            r = an["real"][0]
            tiles.append(("Astral: Bayesian vs decoy FDR @0.9", f"{100 * r['bayes_fdr_090']:.1f}% / "
                                                                 f"{100 * r['decoy_fdr_090']:.1f}%",
                          "posterior calibrated on real data"))
    al = load("align")
    if al:
        s = pd.DataFrame(al["synthetic"])
        x = s[s.batch == "three-trends-strong"].set_index("method")
        tiles.append(("Alignment: RT spread, 3 drift trends", f"{x.loc['macha', 'rt_rms_s']:.2f} s",
                      f"global warp {x.loc['global', 'rt_rms_s']:.2f} s, DTW {x.loc['dtw', 'rt_rms_s']:.2f} s"))
    xr = load("xrun")
    if xr:
        d = pd.DataFrame(xr["drift"])
        c = d[d.fit == "cv"].groupby("anchors_from").loo_rmse_s.mean()
        tiles.append(("Drift without peaks: held-out RT error", f"{c['relations']:.2f} s",
                      f"relation anchors from the scan graph; feature anchors {c['features']:.2f} s"))
    cp = load("components")
    if cp:
        r = cp["deconvolution"]["resolved"]
        tiles.append(("Merged near-isobaric pairs resolved", f"{r['mz']:.0%}",
                      f"components with the m/z trajectory; intensity only {r['int']:.0%}, picked peaks {r['peaks']:.0%}"))
    page.tiles(tiles)
    page.section("Modules", """
<p>Each module is usable on its own; <code>macha.pipeline</code> only fixes the order.</p>
<ul class="prose">
<li><a href="io.html"><b>macha.io / MSRun</b></a>: streaming mzML (local, gzip, HTTP), MS-Numpress, CSR run
container, multi-range stitching, mzML export.</li>
<li><a href="traces.html"><b>macha.traces</b></a>: ROI detection as multi-target tracking with linear assignment
and a self-calibrated, intensity-dependent m/z precision model (replaces rectroi / webtrace / dengroup).</li>
<li><a href="baseline.html"><b>macha.baseline</b></a>: censoring-aware asymmetric Whittaker baselines with Tobit
imputation of unreported scans, detection-floor and noise estimates (replaces LOCF + IRLS / alpha-hull).</li>
<li><a href="peaks.html"><b>macha.peaks</b></a>: persistence (prominence) peak detection, feature table with m/z
standard errors, optional EMG mixture refinement (replaces scale-space seeding + skew-normal refits).</li>
<li><a href="annotate.html"><b>macha.annotate</b></a>: modulo mass-difference graph (isotopes, adducts, losses,
charge states, multimers, homologous series) with decoy-calibrated two-groups Bayesian edge posteriors and MAP graph
partitioning into compounds (extends mz.unity).</li>
<li><a href="scan_graph.html"><b>macha.annotate.scan_graph</b></a>: the same rules tested between centroids scan by
scan with ±1 scan bridging; relations aggregated over time, decoy-calibrated, triangle-supported, and partitioned
into compounds in graph space before any peak detection, under an elution-time prior learned from the run's peak
width.</li>
<li><a href="dia.html"><b>DIA fragment grouping</b></a>: isolation windows as fragment LC-MS runs; residue-ladder
scan graph, ladder-seeded vs co-elution vs MS1-anchored pseudo-spectra, scored by one database search
(<code>macha.search</code>) at 1 % FDR.</li>
<li><a href="align.html"><b>macha.align</b></a>: drift-trend mixture model of retention-time warps with
probabilistic, one-per-run correspondence and membership probabilities (replaces grtgmzcor / warpgroup /
warpcombine DTW).</li>
<li><a href="xrun.html"><b>cross-run graph</b></a> (<code>macha.align.edges</code>): drift learned from scan-graph
relations matched across runs, before peak detection, with centre precision calibrated on cross-run double
differences and noise measured per ROI (<code>macha.traces.noise</code>).</li>
<li><a href="components.html"><b>macha.components</b></a>: peak models of scan-graph relation clusters: one shared
EMG per compound, split where two sources share an ion (with the m/z trajectory of coalesced centroids), linked
across runs and quantified in every run with the learned shape, censored at each ion's own noise floor.</li>
<li><b>macha.synth</b>: the ground-truth simulator behind every accuracy number: fine-structure isotopes, adducts,
EMG elution, Orbitrap-like coalescence, precision and censoring, multi-run batches with class-specific drift.</li>
</ul>""")
    page.write(REPORTS / "index.html")


if __name__ == "__main__":
    if "--index" not in sys.argv:
        for m in MODULES:
            print("==", m)
            runpy.run_module(m, run_name="__main__")
    index()
