"""Ground-truth synthetic LC-MS runs and multi-sample batches.

Every emitted centroid carries the id of the ion that produced it, every ion
the id of its compound, and every compound its drift class, so each module
can be scored exactly.

Model summary
-------------
* compounds: random metabolite formulas or tryptic peptides; elution is an
  EMG (or a long plateau for chemical background); all ions of a compound
  share its elution profile.
* ions: adducts / charge states x isotopic fine structure merged at the
  instrument's resolving power.
* instrument: same-scan coalescence of unresolved ions, intensity-dependent
  m/z precision, multiplicative + additive intensity noise, a detection
  threshold (absent points rather than zeros, as in Orbitrap data), smooth
  calibration drift and random noise centroids.
* batches: per-sample retention-time warps made of a global component plus
  a drift trend specific to each latent compound class, growing with
  injection order; biological abundance variation; per-sample calibration.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd

from .. import chem
from ..run import MSRun
from .instrument import Instrument, coalesce
from .lc import emg, emg_mode, smooth_random_function

POS_ADDUCTS = [  # name, probability present, ratio range relative to primary
    ("[M+Na]+", 0.45, (0.03, 0.9)),
    ("[M+K]+", 0.15, (0.02, 0.2)),
    ("[M+NH4]+", 0.25, (0.03, 0.6)),
    ("[M+H-H2O]+", 0.25, (0.02, 0.4)),
    ("[2M+H]+", 0.10, (0.02, 0.2)),
    ("[M+2H]2+", 0.06, (0.05, 0.4)),
]
NEG_ADDUCTS = [
    ("[M+Cl]-", 0.2, (0.03, 0.4)),
    ("[M+HCOO]-", 0.35, (0.05, 0.8)),
    ("[M-H-H2O]-", 0.15, (0.02, 0.3)),
    ("[2M-H]-", 0.10, (0.02, 0.2)),
    ("[M-2H]2-", 0.04, (0.05, 0.3)),
]


@dataclass
class Scenario:
    name: str = "metabolomics"
    kind: str = "metabolite"  # or "peptide"
    polarity: int = 1
    n_compounds: int = 600
    run_time: float = 900.0
    cycle: float = 0.5
    rt_margin: float = 20.0
    sigma_range: tuple = (1.5, 3.5)
    tau_ratio_range: tuple = (0.1, 1.2)
    log10_height: tuple = (5.0, 0.8)  # mean, sd of apex height (noise=1e3 -> SNR ~100)
    n_background: int = 80
    noise_points_per_scan: float = 120.0
    n_isobaric_pairs: int = 0
    isobaric_ppm: tuple = (1.0, 5.0)
    n_isomers: int = 0  # same formula, eluting isomer_dt seconds apart
    isomer_dt: tuple = (3.0, 20.0)
    instrument: Instrument = field(default_factory=Instrument)


def metabolomics_scenario(**kw) -> Scenario:
    return replace(Scenario(), **kw)


def proteomics_scenario(**kw) -> Scenario:
    base = Scenario(
        name="proteomics", kind="peptide", n_compounds=1200, run_time=1200.0, cycle=0.8,
        sigma_range=(2.5, 6.0), tau_ratio_range=(0.1, 0.8), log10_height=(4.4, 0.7),
        n_background=40, noise_points_per_scan=200.0,
        instrument=Instrument(resolution_200=240_000, noise=300.0, mz_range=(400.0, 1000.0)),
    )
    return replace(base, **kw)


def hard_scenario(**kw) -> Scenario:
    base = metabolomics_scenario(
        name="metabolomics_hard", n_compounds=900, n_isobaric_pairs=150, isobaric_ppm=(4.0, 15.0),
        noise_points_per_scan=400.0,
        n_background=200, log10_height=(4.6, 0.9), sigma_range=(1.5, 3.0),
    )
    return replace(base, **kw)


# --------------------------------------------------------------------------- compounds & ions
def make_compounds(sc: Scenario, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    t0, t1 = sc.rt_margin, sc.run_time - sc.rt_margin
    if sc.kind == "peptide":
        seqs = [chem.random_peptide(rng) for _ in range(sc.n_compounds)]
        hyd = np.array([chem.peptide_hydrophobicity(s) / len(s) ** 0.3 for s in seqs])
        ranks = (np.argsort(np.argsort(hyd)) + rng.uniform(-0.5, 0.5, len(hyd)) * len(hyd) * 0.05) / len(hyd)
        rts = t0 + np.clip(ranks, 0, 1) * (t1 - t0)
        for s, rt in zip(seqs, rts):
            f = chem.peptide_formula(s)
            rows.append(dict(kind="peptide", sequence=s, formula=chem.format_formula(f),
                             neutral_mass=chem.monoisotopic_mass(f), rt=float(rt)))
    else:
        for _ in range(sc.n_compounds):
            f = chem.random_metabolite_formula(rng)
            rows.append(dict(kind="metabolite", sequence="", formula=chem.format_formula(f),
                             neutral_mass=chem.monoisotopic_mass(f), rt=float(rng.uniform(t0, t1))))
    df = pd.DataFrame(rows)
    n = len(df)
    df["sigma"] = rng.uniform(*sc.sigma_range, n)
    df["tau"] = df["sigma"] * rng.uniform(*sc.tau_ratio_range, n)
    df["height"] = 10 ** np.clip(rng.normal(*sc.log10_height, n), 3.3, 9.0) * (sc.instrument.noise / 1000.0)
    df["profile"] = "emg"
    df["mass_offset_ppm"] = 0.0

    extra = []
    for _ in range(sc.n_isobaric_pairs):  # near-isobaric, partly co-eluting interferents
        src = df.iloc[int(rng.integers(len(df)))]
        r = src.to_dict()
        r["mass_offset_ppm"] = float(rng.choice([-1, 1]) * rng.uniform(*sc.isobaric_ppm))
        r["rt"] = float(src["rt"] + rng.uniform(-2.5, 2.5) * src["sigma"])
        r["height"] = float(src["height"] * 10 ** rng.uniform(-1.0, 1.0))
        extra.append(r)
    for _ in range(sc.n_isomers):  # structural isomers: identical m/z, nearby retention
        src = df.iloc[int(rng.integers(len(df)))]
        r = src.to_dict()
        r["rt"] = float(np.clip(src["rt"] + rng.choice([-1, 1]) * rng.uniform(*sc.isomer_dt), t0, t1))
        r["sigma"] = float(rng.uniform(*sc.sigma_range))
        r["tau"] = float(r["sigma"] * rng.uniform(*sc.tau_ratio_range))
        r["height"] = float(10 ** np.clip(rng.normal(*sc.log10_height), 3.3, 9.0) * (sc.instrument.noise / 1000.0))
        extra.append(r)
    for _ in range(sc.n_background):  # chemical background: long plateaus
        f = chem.random_metabolite_formula(rng)
        start = rng.uniform(0, sc.run_time * 0.6)
        dur = rng.uniform(0.25, 1.0) * (sc.run_time - start)
        extra.append(dict(kind="background", sequence="", formula=chem.format_formula(f),
                          neutral_mass=chem.monoisotopic_mass(f), rt=start + dur / 2, sigma=dur / 2,
                          tau=0.0, height=10 ** rng.uniform(0.7, 2.3) * sc.instrument.noise,
                          profile="plateau", mass_offset_ppm=0.0))
    if extra:
        df = pd.concat([df, pd.DataFrame(extra)], ignore_index=True)
    df.insert(0, "compound_id", np.arange(len(df)))
    return df


def _ion_pattern(formula: chem.Formula, adduct: chem.Adduct, inst: Instrument):
    ion_f = adduct.ion_formula(formula)
    masses, probs = chem.isotope_fine_structure(ion_f, prune=1e-5)
    mz = (masses - adduct.charge * chem.ELECTRON) / abs(adduct.charge)
    cmz, cp = chem.centroid_pattern(mz, probs, inst.resolution, merge=inst.merge_fraction, min_rel=2e-3)
    mono = (chem.monoisotopic_mass(ion_f) - adduct.charge * chem.ELECTRON) / abs(adduct.charge)
    iso = chem.nominal_isotope_index(cmz, mono, adduct.charge)
    return cmz, cp, iso


def make_ions(compounds: pd.DataFrame, sc: Scenario, rng: np.random.Generator) -> pd.DataFrame:
    inst = sc.instrument
    rows = []
    for c in compounds.itertuples(index=False):
        f = chem.parse_formula(c.formula)
        if c.kind == "peptide":
            nb = 1 + sum(c.sequence.count(a) for a in "KRH")
            z0 = 2 if len(c.sequence) <= 14 else 3
            zmax = min(nb + 1, 5)
            choices = [(z0, 1.0)]
            if z0 + 1 <= zmax and rng.random() < 0.6:
                choices.append((z0 + 1, rng.uniform(0.1, 0.7)))
            if z0 - 1 >= 1 and rng.random() < 0.5:
                choices.append((z0 - 1, rng.uniform(0.05, 0.5)))
            adducts = [(chem.ADDUCTS[f"[M+{z}H]{z}+" if z > 1 else "[M+H]+"], r) for z, r in choices]
        elif c.kind == "background":
            adducts = [(chem.ADDUCTS["[M+H]+" if sc.polarity > 0 else "[M-H]-"], 1.0)]
        else:
            primary = "[M+H]+" if sc.polarity > 0 else "[M-H]-"
            adducts = [(chem.ADDUCTS[primary], 1.0)]
            for name, p, (lo, hi) in (POS_ADDUCTS if sc.polarity > 0 else NEG_ADDUCTS):
                ad = chem.ADDUCTS[name]
                if abs(ad.charge) > 1 and c.neutral_mass < 300:
                    continue
                if rng.random() < p:
                    adducts.append((ad, float(np.exp(rng.uniform(np.log(lo), np.log(hi))))))
        ref_max = None
        for ad, ratio in adducts:
            try:
                mz, p, iso = _ion_pattern(f, ad, inst)
            except ValueError:  # e.g. water loss from a compound without oxygen
                continue
            if ref_max is None:
                ref_max = p.max()
            mz = mz * (1 + 1e-6 * c.mass_offset_ppm)
            for k in range(len(mz)):
                if not (inst.mz_range[0] <= mz[k] <= inst.mz_range[1]):
                    continue
                rows.append((c.compound_id, ad.name, ad.charge, int(iso[k]), float(mz[k]),
                             float(ratio * p[k] / ref_max)))
    ions = pd.DataFrame(rows, columns=["compound_id", "adduct", "charge", "isotope", "mz", "rel"])
    ions.insert(0, "ion_id", np.arange(len(ions)))
    return ions


# --------------------------------------------------------------------------- simulation
@dataclass
class SyntheticRun:
    run: MSRun
    ion: np.ndarray  # ion id per centroid, -1 for noise
    merged: np.ndarray  # centroid formed by coalescence of >1 ion
    ions: pd.DataFrame
    compounds: pd.DataFrame  # with this sample's rt / height
    scenario: Scenario

    @property
    def compound(self) -> np.ndarray:
        cid = self.ions["compound_id"].to_numpy()
        out = np.full(self.ion.shape, -1, dtype=np.int64)
        ok = self.ion >= 0
        out[ok] = cid[self.ion[ok]]
        return out


def simulate_run(compounds: pd.DataFrame, ions: pd.DataFrame, sc: Scenario, rng: np.random.Generator,
                 calibration_ppm=None, name: str = "sample") -> SyntheticRun:
    inst = sc.instrument
    rt = np.arange(0.0, sc.run_time, sc.cycle) + rng.uniform(-0.02, 0.02) * sc.cycle
    rt.sort()
    comp = compounds.set_index("compound_id")
    ion_c = ions["compound_id"].to_numpy()
    mu = comp["rt"].to_numpy()[ion_c]
    sig = comp["sigma"].to_numpy()[ion_c]
    tau = comp["tau"].to_numpy()[ion_c]
    height = comp["height"].to_numpy()[ion_c] * ions["rel"].to_numpy()
    plateau = (comp["profile"].to_numpy() == "plateau")[ion_c]

    lo = np.where(plateau, mu - sig, mu - 5 * sig)
    hi = np.where(plateau, mu + sig, mu + 5 * sig + 8 * tau)
    s_lo = np.searchsorted(rt, lo)
    s_hi = np.searchsorted(rt, hi)
    counts = np.maximum(s_hi - s_lo, 0)
    ion_idx = np.repeat(np.arange(len(ions)), counts)
    starts = np.repeat(s_lo, counts)
    scan = (starts + np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts)).astype(np.int64)
    t = rt[scan]

    cmu, csig, ctau = comp["rt"].to_numpy(), comp["sigma"].to_numpy(), comp["tau"].to_numpy()
    emg_c = comp["profile"].to_numpy() == "emg"
    norm = np.ones(len(comp))
    if emg_c.any():
        mode = emg_mode(cmu[emg_c], csig[emg_c], ctau[emg_c])
        norm[emg_c] = emg(mode, cmu[emg_c], csig[emg_c], ctau[emg_c])
    shape = np.empty(t.shape)
    pl = plateau[ion_idx]
    e = ~pl
    ii = ion_idx[e]
    shape[e] = emg(t[e], mu[ii], sig[ii], tau[ii]) / norm[ion_c[ii]]
    if pl.any():
        ip = ion_idx[pl]
        # smooth edges + slow multiplicative wander shared by a background compound's ions
        edge = np.minimum(t[pl] - lo[ip], hi[ip] - t[pl]) / 10.0
        cid = ion_c[ip]
        phase = (cid * 2.399963) % (2 * np.pi)
        wander = np.exp(0.35 * np.sin(t[pl] / (60.0 + (cid % 7) * 25.0) + phase))
        shape[pl] = np.clip(edge, 0, 1) * wander
    inten = height[ion_idx] * shape
    keep = inten > 0.2 * inst.threshold
    scan, ion_idx, inten = scan[keep], ion_idx[keep], inten[keep]
    mz_true = ions["mz"].to_numpy()[ion_idx]

    order = np.lexsort((mz_true, scan))
    scan, mz_true, inten, ion_idx = scan[order], mz_true[order], inten[order], ion_idx[order]
    tol = inst.merge_fraction / inst.resolution(mz_true)
    scan, mz_true, inten, ion_idx, merged = coalesce(scan, mz_true, inten, ion_idx.astype(np.int64), tol)

    obs_mz, obs_i, ok = inst.measure(rng, scan, mz_true, inten, rt, calibration_ppm)
    scan, obs_mz, obs_i, ion_idx, merged = scan[ok], obs_mz[ok], obs_i[ok], ion_idx[ok], merged[ok]

    n_noise = rng.poisson(sc.noise_points_per_scan, len(rt))
    nscan = np.repeat(np.arange(len(rt)), n_noise)
    nmz = rng.uniform(*inst.mz_range, nscan.size)
    ni = inst.threshold * (1.0 + rng.exponential(0.5, nscan.size))

    all_scan = np.concatenate([scan, nscan])
    all_mz = np.concatenate([obs_mz, nmz])
    all_i = np.concatenate([obs_i, ni])
    all_ion = np.concatenate([ion_idx, np.full(nscan.size, -1)])
    all_merged = np.concatenate([merged, np.zeros(nscan.size, bool)])
    order = np.lexsort((all_mz, all_scan))
    run = MSRun.from_points(all_scan[order], all_mz[order], all_i[order], rt, sc.polarity,
                            meta={"synthetic": True, "scenario": sc.name, "name": name})
    return SyntheticRun(run=run, ion=all_ion[order], merged=all_merged[order], ions=ions,
                        compounds=compounds, scenario=sc)


def make_run(scenario: Scenario | None = None, seed: int = 0) -> SyntheticRun:
    sc = scenario or metabolomics_scenario()
    rng = np.random.default_rng(seed)
    compounds = make_compounds(sc, rng)
    ions = make_ions(compounds, sc, rng)
    cal = _calibration(rng, sc, 0.5)
    return simulate_run(compounds, ions, sc, rng, calibration_ppm=cal)


def _calibration(rng, sc, scale):
    c0, c1 = rng.normal(0, scale), rng.normal(0, scale / 2)
    T = sc.run_time

    def cal(t, mz):
        return c0 + c1 * (np.asarray(t) / T)

    return cal


# --------------------------------------------------------------------------- batches
def _distinct_profiles(rng, t0, t1, k, amplitude, n_knots=5):
    """k smooth drift profiles, Gram-Schmidt orthogonalised on a grid so classes differ over most of the run."""
    from scipy.interpolate import CubicSpline

    knots = np.linspace(t0, t1, n_knots)
    V = rng.uniform(-1, 1, (k, n_knots))
    for i in range(k):
        for q in range(i):
            V[i] -= (V[i] @ V[q]) / (V[q] @ V[q]) * V[q]
        V[i] *= amplitude / (np.abs(V[i]).max() + 1e-12)
    splines = [CubicSpline(knots, v, bc_type="natural") for v in V]
    return [(lambda cs: (lambda t: cs(np.clip(t, t0, t1))))(cs) for cs in splines]



@dataclass
class SyntheticBatch:
    runs: list
    compounds: pd.DataFrame  # reference table incl. drift_class
    ions: pd.DataFrame
    rt: np.ndarray  # (n_compounds, n_samples) true apex-location per sample
    height: np.ndarray  # (n_compounds, n_samples)
    drift_class: np.ndarray
    injection_order: np.ndarray
    scenario: Scenario
    warps: dict = field(default_factory=dict)


def make_batch(scenario: Scenario | None = None, n_samples: int = 8, n_drift_classes: int = 3,
               class_drift: float = 6.0, global_drift: float = 4.0, wobble: float = 0.8,
               jitter: float = 0.25, bio_cv: float = 0.3, calibration_ppm: float = 0.8,
               class_weights=None, seed: int = 0) -> SyntheticBatch:
    """Simulate ``n_samples`` injections of one compound library.

    The RT of compound ``c`` (drift class ``k``) in sample ``j`` is::

        rt_c + o_j * (G(rt_c) + D_k(rt_c)) + w_j(rt_c) + N(0, jitter)

    with ``o_j`` the scaled injection order, ``G`` a global drift profile,
    ``D_k`` the class drift trend and ``w_j`` a small sample-specific wobble.
    """
    sc = scenario or metabolomics_scenario()
    rng = np.random.default_rng(seed)
    compounds = make_compounds(sc, rng)
    ions = make_ions(compounds, sc, rng)
    n = len(compounds)
    w = np.ones(n_drift_classes) / n_drift_classes if class_weights is None else np.asarray(class_weights)
    drift_class = rng.choice(n_drift_classes, size=n, p=w / w.sum())
    drift_class[compounds["kind"].to_numpy() == "background"] = -1
    compounds["drift_class"] = drift_class

    T0, T1 = 0.0, sc.run_time
    G = smooth_random_function(rng, T0, T1, 5, global_drift)
    D = _distinct_profiles(rng, T0, T1, n_drift_classes, class_drift)
    order = np.linspace(-0.5, 0.5, n_samples) if n_samples > 1 else np.zeros(1)
    runs, rts, hs, wobbles = [], [], [], []
    for j in range(n_samples):
        wob = smooth_random_function(rng, T0, T1, 8, wobble)
        wobbles.append(wob)
        base = compounds["rt"].to_numpy()
        shift = order[j] * G(base) + wob(base)
        for k in range(n_drift_classes):
            m = drift_class == k
            shift[m] += order[j] * D[k](base[m])
        rt_j = base + shift + rng.normal(0, jitter, n)
        rt_j[drift_class < 0] = base[drift_class < 0]
        h_j = compounds["height"].to_numpy() * np.exp(rng.normal(0, bio_cv, n))
        h_j[drift_class < 0] = compounds["height"].to_numpy()[drift_class < 0]
        cj = compounds.copy()
        cj["rt"], cj["height"] = rt_j, h_j
        cal = _calibration(rng, sc, calibration_ppm)
        runs.append(simulate_run(cj, ions, sc, rng, calibration_ppm=cal, name=f"S{j:02d}"))
        rts.append(rt_j)
        hs.append(h_j)
    return SyntheticBatch(runs=runs, compounds=compounds, ions=ions, rt=np.stack(rts, 1),
                          height=np.stack(hs, 1), drift_class=drift_class, injection_order=order,
                          scenario=sc, warps={"global": G, "class": D, "wobble": wobbles})


def true_warps(batch: SyntheticBatch, t) -> np.ndarray:
    """(K, J, len(t)) true class warps, centred across samples like the fitted model's."""
    t = np.asarray(t, float)
    G, D, wob = batch.warps["global"], batch.warps["class"], batch.warps["wobble"]
    o = batch.injection_order
    W = np.stack([np.stack([o[j] * (G(t) + D[k](t)) + wob[j](t) for j in range(len(o))]) for k in range(len(D))])
    return W - W.mean(1, keepdims=True)
