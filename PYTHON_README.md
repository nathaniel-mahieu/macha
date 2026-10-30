# macha (Python)

A re-implementation of macha as a Python library of independent LC-MS processing modules. It keeps the
ideas worth keeping from the R package: ROIs that are more than an m/z slice, modulo mass-difference
discovery (mz.unity), and cross-run warping. Each one is rebuilt on a probabilistic footing, made fast with
numba, and benchmarked against ground truth and public Orbitrap Astral data.

```
src/macha/
  run.py          MSRun: CSR container for centroided MS1 data (+ stitch_scans for multi-range cycles)
  io/             streaming mzML reader/writer (local, .gz, HTTP; MS-Numpress), fetch_ms1;
                  dia.py: DIA as MS1 + one MS2 time series per isolation window, fetch_dia
  synth/          ground-truth simulator: fine-structure isotopes, adducts, EMG elution, Orbitrap-like
                  coalescence / precision / censoring, multi-run batches with class-specific drift;
                  dia.py: DIA windows with b / y fragment ladders, co-isolation and fragment truth
  traces/         ROI detection by multi-target tracking (LAP + self-calibrated precision model);
                  noise.py: interference and noise level in each ROI's own m/z window
  baseline.py     censoring-aware asymmetric Whittaker baselines, detection floor, noise
  peaks.py        persistence peak detection, feature table with m/z standard errors, EMG mixture fits
  annotate/       modulo mass-difference graph + decoy-calibrated Bayesian edges + MAP graph partition;
                  scan_graph: the same graph built from centroids scan by scan, before peak detection;
                  ladders / coelution / complements / fragment_groups: DIA fragment grouping
  search.py       tryptic digest, per-window fragment index, binomial scoring, peptide-level target-decoy FDR
  align/          drift-trend mixture model + probabilistic correspondence (+ DTW, references);
                  edges.py: anchors from scan-graph relations, so drift is learned before peak detection
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

# cross-run drift from scan-graph relations, before peak detection
from macha.align.drift import cv_drift
from macha.align.edges import relation_anchors, relation_noise, relation_table
rels = [relation_table(g, tr) for g, tr in zip(graphs, traces)]      # one per run
A, table, members = relation_anchors(rels, noise=[relation_noise(r, tr) for r, tr in zip(rels, traces)])
model, cv = cv_drift(A.rt, A.se)    # K and smoothing chosen by blocked cross-validation
al = align(tables, anchors=A, K=model.K, lam=cv.loc[cv.rmse.idxmin(), "lam"])   # optional: match features
# DIA: every isolation window is an LC-MS run of fragments; group them into pseudo-spectra and search
from macha.io import read_dia
from macha.annotate.fragment_groups import fragment_groups, ms1_precursors, precursor_groups
from macha.search import digest, read_fasta, with_decoys, WindowIndex, search_groups, peptide_fdr
d = read_dia("run.mzML.gz")                     # d.ms1, d.windows, d.ms2[w]
tr = detect_traces(d.ms2[w])
fg = fragment_groups(tr, (d.windows.lo[w], d.windows.hi[w]))   # residue ladders seed, co-elution grows
fg.groups                                       # ladder steps, longest tag, complement precursor mass
db = with_decoys(digest(read_fasta("human.fasta.gz")))
hits = search_groups(WindowIndex.build(db, (d.windows.lo[w], d.windows.hi[w])), db,
                     [tr.table.mz.to_numpy()[m] for m in fg.members])
peptide_fdr(hits, 0.01)
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
| DIA grouping (synthetic) | peptides at 1 % FDR, same search for every grouping (2 / 20 Th windows): MS1-anchored + ladder rescue 4114 / 3462, MS1-anchored 3963 / 3473, ladder-seeded 2384 / 1717; measured FDR 0.7-1.4 % | co-elution only 2601 / 1444, ladder graph alone 228 / 5 |
| DIA grouping (real) | Astral 250 pg, 4 of 20 windows: MS1-anchored + ladder rescue 438, MS1-anchored 423, ladder-seeded 318 (51 found by no other grouping) | co-elution only 229 |
| annotation (real) | Bayesian vs decoy FDR at posterior ≥ 0.9: 1.6 % vs 1.6 % (250 pg), 2.3 % vs 1.8 % (50 ng), 2.7 % vs 0.5 % (metabolomics, conservative) | — |
| alignment | three drift trends, strong: F1 0.973, aligned RT spread 1.17 s | global warp 0.951 / 4.91 s, uncorrected 0.948 / 5.70 s, BPC-DTW 0.941 / 6.80 s, density 0.886 |
| alignment (real) | 850k features from 6 Astral runs in ~90 s; tighter anchor fit than a global warp (σ 0.32 vs 0.47 s) | isotope consistency equal to a global warp |
| cross-run graph | drift learned from scan-graph relations without peak detection: held-out RT error 0.28 / 0.36 / 0.39 / 0.29 s on the four alignment batches; relation centres 0.06 s vs feature apexes 0.30 s from truth | feature anchors 0.28 / 0.34 / 0.36 / 0.30 s (K, λ by CV); current pipeline 0.28 / 0.41 / 0.77 / 0.29 s |
| cross-run graph (real) | metabolomics, same held-out blocks: relation model predicts feature anchors 0.81 s and relation anchors 0.55 s; ProteoBench M7: 1.14 / 0.79 s; M7 relation groups, no peaks, SNR ≥ 30: yeast mode at log2 A/B +1.03 (design +1), sd 0.25 | feature model 0.87 / 0.59 s (metabolomics), 1.08 / 0.82 s (M7); feature groups yeast +1.01, sd 0.19; human mode sd 0.23 vs 0.40 |

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

**DIA fragments.** A DIA isolation window visited once per cycle is an LC-MS run of fragment ions, so trace
detection and the scan-level graph apply unchanged, with amino-acid residue steps as the modulo rules. Residue
decoys need care: fragments of unrelated peptides also differ by peptide-like masses, so decoys shifted off the
mass-defect line (the MS1 scheme) would be hit too rarely. They are placed on the line, more than 0.02 Da from
every residue, residue pair, isotope or loss step, keeping each residue's offset from the line. Persistence does
not separate ladder steps from coincidences here: every fragment of a precursor co-elutes. Ladder components are
therefore near-pure but sparse, with one component per ion series. They work best as seeds: merged when their summed chromatograms
match, then grown by the fragments that follow them. Complementary pairs (b<sub>i</sub> + y<sub>L-i</sub> = M + 2H⁺)
give a precursor mass without MS1. `benchmarks/bench_dia.py` scores every grouping with one database search
(`macha.search`) at 1 % peptide-level FDR.

**Alignment.** Per-m/z-channel DTW (warpgroup/warpcombine) and a single global curve (`grtgmzcor`) are replaced
by a **drift-trend mixture**: each compound belongs to a latent class whose P-spline warp is shared by all its
members. It is learned by EM on anchors with heteroscedastic noise, an outlier component, and BIC to choose K.
Correspondence is solved per m/z channel by repeated linear assignment under the posterior predictive (Student-t)
RT of the class mixture, with leave-one-out membership probabilities. The noise scale is updated by Fisher
scoring on the heteroscedastic likelihood. The earlier moment estimate, Σ(r² − S²)/n, was pinned to its floor
whenever a few anchors stated a large error.

**Cross-run graph (alignment before peaks).** This is the first step towards grouping across runs at the
scan-graph stage. Feature anchors must be unique in m/z over a wide RT window and need picked peaks. A confident
scan-graph relation is keyed by its rule and the m/z of both ions, so it is more specific and exists earlier.
Relations are matched across runs by chain clustering on rule, anchor m/z, partner m/z and RT.
- **Anchor filters.** A relation anchors only if it elutes (its ions' signal is localised: chemical-background
  plateaus have arbitrary centres) and has no competitor. A competitor is another ion within 6 ppm reaching 20 % of
  its height in the matching window; the ion's own trace counts beyond a valley, so isomers joined into one trace
  are caught. Without this, when one of an isomer pair is missed in some runs, the other is matched in its place. A
  relation must also be unlikely by chance in its ions' noise-floor windows and unique in its run.
- **Matching window from the drift.** A first pass at ±60 s measures the drift. The window is 2 × the worst run's
  95th-percentile deviation, at least 10 s, and the cross-run span may reach twice that. In dense Astral
  proteomics MS1, a stronger ion within 6 ppm sits somewhere in ±60 s for most relations.
- **Centres.** The relation fixes which ions and roughly when. The centre comes from both ions' full traces in a
  fixed ±2σ window (σ learned by the scan graph), weighted by intensity². The scans where both ions were seen are
  truncated by an amount that depends on abundance, and their centre had twice the error.
- **Precision.** The standard error is √((k₁·sd/SNR)² + (k₂·sd)²), with SNR measured against the noise in the
  ions' own m/z windows. k₁ and k₂ are fitted on cross-run double differences of two relations of one compound
  that share no ion. Drift, the compound's RT error and per-relation bias cancel; the within-run estimate
  overstates the error about twofold.
- **Duplicates kept.** Relations sharing an ion measure one elution. Keeping one per group leaves drift classes
  too sparse in places, and averaging lets one wrong link merge two compounds. With K and λ chosen by blocked
  cross-validation (`cv_drift`, the P-spline penalty is not scale-free), duplicates no longer inflate K.
- **Drift model.** Where a class has no anchors, its warp is shrunk towards the run's class-weighted mean warp
  instead of towards no drift.
- **Cross-run groups without peaks.** Once the drift is known, `relation_groups` groups every relation across runs
  in reference time and quantifies it per run. This is cross-run grouping at the scan-graph stage.
- **Noise per ROI.** Noise is measured in each ROI's own acceptance window, not in m/z × RT tiles, which miss
  local structure by orders of magnitude on real data.
- **Mixed-species test.** On ProteoBench M7 (A/B = 1, 2, 0.25 for human, yeast and E. coli), log2 A/B of cross-run
  groups is fitted as a three-mode mixture with known spacing. Mode widths measure grouping and quantification
  without identifications.
- **Tried and dropped.** Flank-background-corrected heights: no gain.

## Known limitations / next steps

- Drift classes are identifiable only locally. Where two trends cross, labels can swap between RT regions
  without changing predictions, so class recovery is reported per RT window.
- BIC (and ICL) choose one more drift class than the simulator used (3 → 4, 2 → 3; 1 → 1 is right). The extra
  class buys real likelihood from structure the model cannot otherwise represent. On real replicate Astral runs
  the mixture fits anchors more tightly than a global warp, but correspondence quality, measured by isotope
  consistency, is equal or slightly lower. `cv_drift` now chooses K and the smoothing penalty by blocked
  cross-validation (used for relation anchors; `align` still defaults to BIC with λ = 1). Next: constrain warps
  to a low-rank form (shared run loadings × class RT profiles), which also makes class identity global.
- The cross-run graph learns drift and groups relations across runs, but it has no presence/absence model yet.
  On M7, ratios of low-SNR relation groups are compressed towards 0, and E. coli (log2 A/B −2, so 4× lower in A)
  is not resolved by relation groups or feature groups. This is the noise floor and dropout the design matrix and
  detectability model should handle: score presence and absence against the per-ROI noise level, re-query missing
  runs in raw centroids, and use RT-shifted decoys in the ions' own windows. At high SNR, relation groups place
  the yeast mode correctly, but the human mode is wider than for picked peaks (sd 0.40 vs 0.23). Window sums
  without peak boundaries pick up co-eluting signal in the same trace. Still open: node-level correspondence and a
  partition with one centre per run.
- The Fisher-scoring noise scale lowers σ on real data (metabolomics 0.30 → 0.18 s), which tightens the RT gate.
  There, complete groups drop 5 % and isotope RT disagreement rises from 0.055 to 0.061 s; synthetic and 250 pg
  results are equal or better.
- Traces are deliberately permissive at low intensity, where centroids genuinely are imprecise. Peak detection is
  the filter.
- Under heavy coalescence (unresolved near-isobaric co-eluting ions), annotation posteriors at the default 0.5
  cut are conservative. The posterior ranking still dominates the rule-based baselines at equal precision.
- The scan graph's elution-time prior splits long-lived background ions into time-local pieces (background ARI
  0.28 → 0.02 on synthetic data, while eluting compounds improve 0.86-0.90 → 0.91-0.92). The remaining impurity
  comes from coincidental relations between ions that truly co-elute. Next: negative evidence between co-eluting
  blocks with inconsistent intensity ratios, a one-envelope-per-block constraint, and a reverse (heavy-anchored)
  search so relations whose lighter ion never forms a trace are found too.
- DIA: ladder relations stay clean in narrow windows but get muddled in wide windows of dense samples (20 Th: ~70 %
  of confident residue steps are real, and the decoys under-estimate the rest), and ladder components alone are too
  sparse to identify peptides. As seeds without MS1 they beat co-elution-only grouping where windows are crowded,
  and as a rescue next to MS1-anchored grouping they add ~3-4 %. The search (`macha.search`) is a deliberately
  simple binomial scorer with no library, intensity or retention-time prediction, so absolute identification
  numbers are far below DIA-NN-class tools.
- Profile-mode data must be centroided first. MS2 is used only for DIA fragment grouping.
