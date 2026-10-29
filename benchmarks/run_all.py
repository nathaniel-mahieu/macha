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

MODULES = ["bench_io", "bench_traces", "bench_baseline", "bench_peaks", "bench_annotate", "bench_align"]


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
<li><a href="align.html"><b>macha.align</b></a>: drift-trend mixture model of retention-time warps with
probabilistic, one-per-run correspondence and membership probabilities (replaces grtgmzcor / warpgroup /
warpcombine DTW).</li>
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
