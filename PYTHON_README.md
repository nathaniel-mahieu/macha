# macha (Python)

A re-implementation of macha as a Python library of independent LC-MS processing modules. It keeps the
ideas worth keeping from the R package: ROIs that are more than an m/z slice, modulo mass-difference
discovery (mz.unity), and cross-run warping. Each one is rebuilt on a probabilistic footing, made fast with
numba, and benchmarked against ground truth and public Orbitrap Astral data.

```
src/macha/
  run.py          MSRun: CSR container for centroided MS1 data (+ stitch_scans for multi-range cycles)
  io/             streaming mzML reader/writer (local, .gz, HTTP; MS-Numpress), fetch_ms1
  synth/          ground-truth simulator: fine-structure isotopes, adducts, EMG elution, Orbitrap-like
                  coalescence / precision / censoring, multi-run batches with class-specific drift
  traces/         ROI detection by multi-target tracking (LAP + self-calibrated precision model)
  baseline.py     censoring-aware asymmetric Whittaker baselines, detection floor, noise
  peaks.py        persistence peak detection, feature table with m/z standard errors, EMG mixture fits
  annotate/       modulo mass-difference graph + decoy-calibrated Bayesian edges + MAP graph partition;
                  scan_graph: the same graph built from centroids scan by scan, before peak detection
  align/          drift-trend mixture model + probabilistic correspondence (+ DTW, references)
  pipeline.py     process_run / process_batch convenience wrappers
  report.py       themed, self-contained HTML explorers
tests/            pytest suite (synthetic ground truth, property-based checks, real-data checks)
benchmarks/       one benchmark + explorer per module -> reports/*.html
scripts/          fetch_data.py: stream public Astral runs, keep MS1
```

## Install

```bash
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python -e ".[dev]"
```

Python ≥ 3.11 with numpy, scipy, numba and pandas. Plotly is only needed for the explorers (`[report]`).

## Use

```python
from macha import MSRun
from macha.io import read_mzml
from macha.pipeline import process_run
from macha.annotate import annotate, ProfileSource
from macha.align import align

run = read_mzml("sample.mzML.gz")                    # MS1 centroids only
s = process_run(run)                                 # traces -> baselines -> peaks
s.features                                           # mz, mz_se, rt, height, area, snr, fwhm_s, trace, ...

ann = annotate(s.features, ProfileSource(s.baselines), polarity=1, peptides=False, cycle=run.cycle_time)
ann.edges        # candidate relations with evidence, log-odds, posterior
ann.groups       # one row per inferred compound: neutral mass, charge, adducts
ann.roles        # per feature: group, charge, isotope index, adduct, neutral mass
ann.fdr          # Bayesian vs decoy-estimated FDR by posterior threshold

samples = [process_run(read_mzml(p)) for p in paths]
al = align([x.features for x in samples])
al.groups        # consensus features: m/z, consensus RT, drift class posterior
al.intensity     # groups x runs
al.membership    # per run: membership probability of every feature
al.drift         # the fitted drift-trend mixture (K classes chosen by BIC)

# or the whole chain, ending with annotation of consensus features, where cross-sample
# intensity correlation joins the evidence (isotopes/adducts co-vary across samples)
from macha.pipeline import process_batch
b = process_batch([read_mzml(p) for p in paths])
b.alignment, b.consensus, b.annotation
# compounds in graph space, before peak detection: rules tested between centroids scan by scan
from macha.annotate.scan_graph import repartition, scan_graph
g = scan_graph(s.traces, polarity=1, peptides=False, bridge=1)
g.segments       # relation segments over time with evidence and posterior
g.compounds      # graph-space compounds with time extent and member-centre span
g.prior          # elution-time prior learned from the run: peak spread, FWHM, member-centre scatter
g.compound_chromatogram(s.traces, c)   # summed chromatogram of compound c
repartition(g, s.traces.run, p_merge=0.95)   # re-cut compounds without recomputing relations
```

Every step is a plain function on plain data. Run `detect_traces`, `trace_baselines`, `detect_peaks`,
`annotate` or `align` on their own, with your own parameters, or on features from another tool:
`annotate` and `align` only need a feature table.

## Data

```bash
python -I scripts/fetch_data.py                  # all datasets
python -I scripts/fetch_data.py proteomics_sc    # one dataset
```

Raw mzML is streamed and parsed on the fly, and only compact MS1 `.npz` files are written to `data/astral/`
(git-ignored):

| key | source | content |
|---|---|---|
| `proteomics_sc` | PRIDE PXD049412 / ProteoBench module 9 (OpenMS archive mirror) | 250 pg HeLa + yeast, Orbitrap Astral nDIA, MS1 240k, 2 conditions × 3 replicates |
| `proteomics` | ProteoBench module 7 (OpenMS archive mirror) | 50 ng LFQ, 15 min gradient, Orbitrap Astral DIA |
| `metabolomics` | MassIVE MSV000093526 (Dorrestein lab) | NIST fecal reference material, RP positive, MS1 120k + Astral DDA; Sa04 injected repeatedly; multi-range MS1 files |

MassIVE rate-limits downloads (HTTP 429). The fetcher backs off and retries, and finished files are skipped on
re-runs.

## Tests

```bash
.venv/bin/python -m pytest              # ~1.5 min; real-data tests skip if data/astral is absent
.venv/bin/python -m pytest -m realdata  # only the public-data checks
```

The suite tests algorithms against brute force (Hungarian vs scipy, persistence vs `peak_prominences`,
banded Whittaker vs dense solve, targeted candidate search vs the modulo sweep, MS-Numpress against reference
encoders), checks invariants with hypothesis (one centroid per scan per trace, optimal matchings), scores every
module end to end against synthetic truth, and on real Astral data asserts properties that need no truth: a
physically plausible self-calibrated precision model, decoy-calibrated annotation FDR, and isotope-consistent
alignment.

## Benchmarks and explorers

```bash
cd benchmarks && ../.venv/bin/python run_all.py    # writes reports/*.html and reports/*.json
```

`reports/index.html` links one page per module. Each page explains the method, compares it with reference
methods (including the original macha algorithms) on synthetic datasets with ground truth and on the public
Astral runs, and has interactive views: method-switchable trace windows, baseline examples, peak-detection
curves, annotation precision–recall / calibration / FDR / modulo views and compound groups, and drift trends per
run with membership calibration.

## Results at a glance

Synthetic numbers are scored against ground truth. Real numbers come from the public Astral runs, where no truth
exists. Full tables are in `reports/*.html` and `reports/*.json`.

| module | headline | reference methods |
|---|---|---|
| traces | ARI 0.971 / completeness 0.951 on the hard set (near-isobaric co-eluting ions, dense noise) | greedy ROI 0.938 / 0.857, m/z slice 0.884 / 0.735, rectroi 0.949 / 0.957 (purity 0.948) |
| traces (real) | 2.6M centroids in 4.4 s; learned precision floor 0.28 ppm, σ ∝ I^-1.03 | — |
| baseline | RMSE 0.35σ (drift), 0.35σ (constant), 1.15σ (step), 0.30σ (none) | ALS with last-value fill 1.58σ / 2.04σ / 2.82σ / 3.29σ |
| peaks | F1 0.899, precision 1.00, apex error 0.44 scans | scale-space 0.682, CWT 0.880 |
| annotation | edge precision 0.987–0.993, compound ARI 0.95–0.97 (3 of 4 sets); at matched precision the posterior ranking out-recalls rule + RT on every set | rule + RT ARI 0.87–0.97 |
| scan graph | compounds found at SNR 10-30: 59-60 % vs 22-28 % for the feature graph (metabolomics), 86 % vs 67 % (proteomics); at SNR 3-10 in proteomics 42 % vs 5 % | group purity 0.99 (0.93 without the elution-time prior) vs 0.99-1.00 for the feature graph |
| scan graph (real) | 2.4M centroids in 22 s; Bayesian vs decoy FDR at posterior ≥ 0.9: 1.9 % vs 1.8 % (metabolomics), 4.6 % vs 4.6 % (250 pg); 25-59 % of confident relations have no picked peak; compounds whose members span > 3 learned peak widths: 3.4 % / 0.8 % / 3.4 % | 18.6 % / 9.7 % / 10.8 % without the elution-time prior |
| annotation (real) | Bayesian vs decoy FDR at posterior ≥ 0.9: 1.6 % vs 1.6 % (250 pg), 2.3 % vs 1.8 % (50 ng), 2.7 % vs 0.5 % (metabolomics, conservative) | — |
| alignment | three drift trends, strong: F1 0.973, aligned RT spread 1.18 s | global warp 0.951 / 4.91 s, uncorrected 0.948 / 5.70 s, BPC-DTW 0.941 / 6.80 s, density 0.886 |
| alignment (real) | 850k features from 6 Astral runs in ~90 s; tighter anchor fit than a global warp (σ 0.32 vs 0.47 s) | isotope consistency equal to a global warp |

## Design notes: what changed and why

**Traces (ROI).** `rectroi` + `webtrace` (iterative min-cuts on a k-NN graph) is replaced by scan-by-scan
multi-target tracking. Each active trace has an inverse-variance m/z estimate with exponential forgetting. A link
costs a negative log-likelihood under an intensity-dependent precision model `σ(I) = sqrt(a² + (b·I^-c)²)`, plus
intensity continuity and gap terms, and is accepted only if it beats the cost of a new trace. Conflicts are
solved exactly by the Hungarian algorithm per conflict component, and a second assignment closes short gaps. The
precision model is fitted by maximum likelihood on consecutive-centroid differences, which are drift-free. On
Astral data it recovers an Orbitrap-like `σ ∝ 1/SNR` law with a ~0.3 ppm floor.

**Baselines.** Absent centroids are upper bounds (≤ the detection floor), not zeros and not the last value.
Asymmetric Whittaker smoothing (arPLS weights, banded O(n) solver) is combined with Tobit-style imputation of
censored scans from the local fraction of baseline-like observations.

**Peaks.** Scale-space seeding plus iterative skew-normal refits are replaced by 0-dimensional persistence
(prominence) on lightly smoothed traces, with an optional EMG mixture fit that uses an analytic Jacobian. Features
carry m/z standard errors and apex uncertainty, which downstream modules use as measurement error.

**Annotation (mz.unity, generalised).** Candidates come from residues of charge-scaled masses modulo each rule's
step, which gives isotope chains, adduct exchanges, losses and homologous series in one formalism. Charge states
and multimers are matched through neutral masses. Each edge is evaluated by a two-groups empirical-Bayes model:
the null comes from decoy rules, the alternative is parametric with robust shape/ratio terms, and EM fits a prior
per rule and charge. Compounds are the MAP partition (correlation clustering) of the graph, followed by greedy
isotope-envelope selection, an averagine check for a missing monoisotopic peak, and adduct / neutral-mass
inference. The decoy FDR gives an independent calibration check on real data.

**Scan-level graph (compounds before peaks).** The feature-level graph inherits every miss of peak detection. The
scan graph runs the modulo test between *centroids*: every centroid of a trace is an anchor, the best partner for
each rule target is searched in the same scan and, if absent, ±1 scan (bridging partners near the noise floor), and
partners need not form a trace of their own. Overlapping ~2 Da isotope targets are merged into one family. Per-scan
instances are cut into relation segments over consecutive scans. A real relation persists: one partner, a residual
whose mean shrinks as 1/√n, a constant intensity ratio, co-varying intensities. Segments are scored by a weighted
logistic model of targets vs decoys (target-decoy local FDR), refitted with *triangle support* from the provisional
graph, because real relations form dense structures and coincidences rarely close triangles. Traces are cut into
time blocks covering their relations; blocks are nodes, segments are edges, and correlation clustering gives
compounds with a time extent and a summed compound chromatogram on which peaks can be picked. Relations are local
in time, so the partition carries an **elution-time prior** at the run's own peak width. Without it, long-lived
background ions act as hubs that join peaks minutes apart. The prior has two parts. First, the bulk of each
node's signal must lie within one learned peak width of the relation that links it: an edge is discounted by the
log share of each node's signal it explains. Second, member signal centres scatter around a shared centre by a
learned τ plus their own standard error, which adds a Ward penalty inside GAEC. Peaks near the noise floor need
not apex together; their centres are simply less precise.

**Alignment.** Per-m/z-channel DTW (warpgroup/warpcombine) and a single global curve (`grtgmzcor`) are replaced
by a **drift-trend mixture**: each compound belongs to a latent class whose P-spline warp is shared by all its
members. It is learned by EM on anchors with heteroscedastic noise, an outlier component, and BIC to choose K.
Correspondence is solved per m/z channel by repeated linear assignment under the posterior predictive (Student-t)
RT of the class mixture, with leave-one-out membership probabilities.

## Known limitations / next steps

- Drift classes are identifiable only locally. Where two trends cross, labels can swap between RT regions
  without changing predictions, so class recovery is reported per RT window.
- BIC (and ICL) choose one more drift class than the simulator used (3 → 4, 2 → 3; 1 → 1 is right). The extra
  class buys real likelihood from structure the model cannot otherwise represent. On real replicate Astral runs
  the mixture fits anchors more tightly than a global warp, but correspondence quality, measured by isotope
  consistency, is equal or slightly lower. Next step: choose K by held-out-run predictive likelihood, and
  constrain warps to a low-rank form (shared run loadings × class RT profiles), which also makes class identity
  global.
- Traces are deliberately permissive at low intensity, where centroids genuinely are imprecise. Peak detection is
  the filter.
- Under heavy coalescence (unresolved near-isobaric co-eluting ions), annotation posteriors at the default 0.5
  cut are conservative. The posterior ranking still dominates the rule-based baselines at equal precision.
- The scan graph's elution-time prior splits long-lived background ions into time-local pieces (background ARI
  0.28 → 0.02 on synthetic data, while eluting compounds improve 0.86-0.90 → 0.91-0.92). The remaining impurity
  comes from coincidental relations between ions that truly co-elute. Next: negative evidence between co-eluting
  blocks with inconsistent intensity ratios, a one-envelope-per-block constraint, and a reverse (heavy-anchored)
  search so relations whose lighter ion never forms a trace are found too.
- Profile-mode data must be centroided first. MS2 is not used.
