"""Benchmark + explorer: reading data, the MSRun container and the datasets used throughout."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from common import ASTRAL, NAV, REAL_LABEL, REPORTS, SYNTHETIC, Timer, real_runs, save_json, synthetic
from macha import MSRun, stitch_scans
from macha.io import read_mzml, write_mzml
from macha.report import SERIES, Page, style


def main():
    page = Page("I/O & runs", "Streaming mzML reading (local, gzip or HTTP) that keeps only MS1 centroids, and the "
                "compact CSR run container every module operates on.", NAV)
    # local parse throughput on a synthetic run written to mzML
    sr = synthetic("syn-proteomics")
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for suffix in (".mzML", ".mzML.gz"):
            path = Path(tmp) / f"x{suffix}"
            write_mzml(sr.run, path)
            size = path.stat().st_size
            with Timer() as t:
                back = read_mzml(path)
            assert back.n_points == sr.run.n_points
            rows.append(dict(file=f"synthetic proteomics{suffix}", MB=size / 1e6, spectra=back.n_scans,
                             centroids=back.n_points, seconds=t.s, MB_per_s=size / 1e6 / t.s,
                             centroids_per_s=back.n_points / t.s))
        with Timer() as t:
            sr.run.save(Path(tmp) / "x.npz")
        size = (Path(tmp) / "x.npz").stat().st_size
        with Timer() as t2:
            MSRun.load(Path(tmp) / "x.npz")
        rows.append(dict(file="same run as MSRun .npz", MB=size / 1e6, spectra=sr.run.n_scans,
                         centroids=sr.run.n_points, seconds=t2.s, MB_per_s=size / 1e6 / t2.s,
                         centroids_per_s=sr.run.n_points / t2.s))
    io = pd.DataFrame(rows)
    # real datasets
    real = []
    for ds in ("proteomics_sc", "proteomics", "metabolomics"):
        man = ASTRAL / ds / "manifest.json"
        if not man.exists():
            continue
        m = json.loads(man.read_text())
        for name, f in m.get("files", {}).items():
            real.append(dict(dataset=REAL_LABEL[ds], run=name, scans=f["n_scans"], centroids=f["n_points"],
                             rt_end_s=f["rt_range"][1], streamed_MB=(f.get("bytes_streamed") or 0) / 1e6,
                             stream_seconds=f.get("seconds"), source=f["url"].split("/")[2]))
    page.tiles([
        ("local mzML parse", f"{io.centroids_per_s.iloc[0] / 1e6:.1f} M centroids/s", "uncompressed, zlib arrays"),
        ("MSRun .npz load", f"{io.centroids_per_s.iloc[2] / 1e6:.0f} M centroids/s", "12 bytes per centroid"),
        ("public Astral runs", f"{len(real)}", "MS1 extracted while streaming"),
    ])
    page.section("What it does", """
<p><code>read_mzml</code> parses mzML in a single forward pass with an event-driven XML parser, decoding binary
arrays (zlib, 32/64-bit, MS-Numpress linear/pic/slof) only for MS1 spectra and skipping MS2 without decoding, so DIA
and DDA files reduce to their MS1 content as they stream. <code>fetch_ms1</code> does this over HTTP with retry and
backoff, which is how the public Astral data below was obtained without storing raw files. Spectra are sorted by m/z
and packed into an <code>MSRun</code>: flat m/z / intensity arrays with per-scan offsets (CSR), retention times,
scan windows and metadata. <code>stitch_scans</code> merges multi-range acquisitions (several MS1 scans per cycle
covering different m/z windows, the original <code>rawdata(splits=…)</code>) into one scan per cycle.
<code>write_mzml</code> exports a run back to mzML.</p>""")
    page.table(io, "Read / write throughput on this machine", better=dict(centroids_per_s="max"))
    if real:
        page.section("Public Orbitrap Astral datasets", "<p>Fetched by <code>python -I scripts/fetch_data.py</code>. "
                     "PXD049412 / ProteoBench module 9: 250 pg HeLa+yeast single-cell-level nDIA, MS1 at 240k; "
                     "ProteoBench module 7: 50 ng LFQ DIA, 15 min gradient; MSV000093526 (Dorrestein lab): NIST fecal "
                     "reference material, RP positive, MS1 at 120k with Astral DDA, Sa04 injected repeatedly across "
                     "the batch.</p>")
        page.table(pd.DataFrame(real))
        for ds in ("proteomics_sc", "proteomics", "metabolomics"):
            files = [f for f in real_runs(ds) if "multirange" not in f.name]
            if not files:
                continue
            fig = go.Figure()
            for k, f in enumerate(files[:8]):
                run = MSRun.load(f)
                b = run.bpc()
                kern = np.ones(9) / 9
                fig.add_trace(go.Scatter(x=run.rt / 60, y=np.convolve(b, kern, "same"), mode="lines", name=f.stem,
                                         line=dict(width=2, color=SERIES[k])))
            fig.update_xaxes(title="retention time (min)")
            fig.update_yaxes(title="base peak intensity (smoothed)", type="log")
            page.figure(style(fig, 320, title=f"Base peak chromatograms, {REAL_LABEL[ds]}"))
        dens = []
        for ds in ("proteomics_sc", "proteomics", "metabolomics"):
            for f in real_runs(ds)[:1]:
                run = MSRun.load(f)
                dens.append((REAL_LABEL[ds], np.diff(run.offsets)))
        for name in SYNTHETIC:
            dens.append((name, np.diff(synthetic(name).run.offsets)))
        fig = go.Figure()
        for k, (name, d) in enumerate(dens):
            fig.add_trace(go.Box(x=d, name=name, marker=dict(color=SERIES[k % 8]), line=dict(width=1),
                                 boxpoints=False))
        fig.update_xaxes(title="centroids per MS1 scan", type="log")
        page.figure(style(fig, 80 + 40 * len(dens), legend=False,
                          title="Data density: real Astral runs vs the synthetic scenarios"))
        multi = sorted((ASTRAL / "metabolomics").glob("*multirange*.npz"))
        if multi:
            run = MSRun.load(multi[0])
            w = run.scan_window
            n_seg = int(np.unique(np.round(w[:20, 0])).size)
            st = stitch_scans(run, n_seg)
            page.section("Multi-range acquisition", f"<p>{multi[0].stem}: {n_seg} MS1 windows per cycle "
                         f"({', '.join(f'{a:.0f}-{b:.0f}' for a, b in np.unique(np.round(w[:n_seg]), axis=0))}); "
                         f"{run.n_scans} scans stitched into {st.n_scans}.</p>")
    save_json("io", dict(throughput=io.to_dict("records"), real=real))
    page.write(REPORTS / "io.html")
    print(io.round(2).to_string())


if __name__ == "__main__":
    main()
