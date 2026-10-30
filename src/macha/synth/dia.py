"""Synthetic data-independent acquisition with fragment-level ground truth.

Peptides elute as in :func:`macha.synth.generate.make_run`.  The precursor m/z range of interest is cut into
isolation windows visited once per cycle.  Every precursor charge state contributes to the window(s) holding its
isotopes, in proportion to the isolated share, and fragments into b / y ion series (1+, and 2+ for 3+ and higher
precursors) with their isotopes.  Fragments share the precursor's elution profile.  Each window is rendered by
:func:`~macha.synth.generate.simulate_run`, so the same instrument model (precision, threshold, noise centroids)
applies, with an Astral-like analyser for MS2.

Truth: in window ``w`` the centroid ``ms2[w].ion[i]`` is fragment ion ``ms2[w].ions.ion_id`` whose
``compound_id`` is the peptide; ``series`` / ``index`` / ``frag_z`` locate it in its ladder.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
import pandas as pd

from .. import chem
from ..io.dia import DIARun
from .generate import Scenario, SyntheticRun, make_compounds, make_ions, proteomics_scenario, simulate_run
from .instrument import Instrument

H2O = chem.monoisotopic_mass({"H": 2, "O": 1})
CO = chem.monoisotopic_mass({"C": 1, "O": 1})
RES = {aa: chem.monoisotopic_mass(chem.parse_formula(f)) for aa, f in chem.AA_FORMULA.items()}


def astral_ms2() -> Instrument:
    """Astral-like MS2 analyser: ~80k resolving power at m/z 524, ~3 ppm at moderate signal."""
    return Instrument(resolution_200=80_000 * np.sqrt(524 / 200), noise=100.0, ppm_floor=0.6, ppm_snr1=12.0,
                      rel_noise=0.15, mz_range=(150.0, 1600.0))


@dataclass
class DIAScenario:
    precursors: Scenario = field(default_factory=lambda: proteomics_scenario(n_compounds=6000, n_background=0))
    window_width: float = 20.0
    mz_span: tuple = (500.0, 580.0)  # precursor m/z covered by the simulated windows
    ms2: Instrument = field(default_factory=astral_ms2)
    gain: float = 1.0  # base-peak fragment SNR relative to precursor apex SNR
    ms2_noise_points: float = 150.0  # random noise centroids per MS2 spectrum


@dataclass
class SyntheticDIA:
    ms1: SyntheticRun
    windows: pd.DataFrame
    ms2: list[SyntheticRun]
    compounds: pd.DataFrame
    scenario: DIAScenario

    def as_dia_run(self) -> DIARun:
        return DIARun(self.ms1.run, self.windows, [s.run for s in self.ms2], {"synthetic": True})


def fragments(seq: str, prec_z: int, rng: np.random.Generator, mz_range=(150.0, 1600.0)) -> pd.DataFrame:
    """b / y fragments (with a2 and 13C isotopes) of a peptide and random but realistic relative intensities."""
    res = np.array([RES[a] for a in seq])
    L = len(seq)
    pre = np.cumsum(res)
    rows = []
    zs = (1, 2) if prec_z >= 3 else (1,)
    for i in range(1, L):  # cleavage after residue i (1-based): b_i and y_{L-i}
        pro = 4.0 if seq[i] == "P" else 1.0  # enhanced cleavage N-terminal to proline
        asp = 2.0 if seq[i - 1] == "D" else 1.0
        for series, n, mass in (("b", i, pre[i - 1]), ("y", L - i, pre[-1] - pre[i - 1] + H2O)):
            if series == "b" and n == 1:
                continue  # b1 ions are rarely observed
            base = (1.0 if series == "y" else (0.8 if n == 2 else 0.3)) * pro * asp
            base *= np.exp(rng.normal(0, 0.9))
            for fz in zs:
                if fz == 2 and n < 4:
                    continue
                rel = base * (1.0 if fz == 1 else 0.35)
                mz0 = (mass + fz * chem.PROTON) / fz
                for k in range(3):
                    r = rel * (1.0 if k == 0 else chem.expected_isotope_ratio(mass, k, "peptide"))
                    m = mz0 + k * chem.NEUTRON_C13 / fz
                    if mz_range[0] <= m <= mz_range[1] and r > 0.01 * rel:
                        rows.append((series, n, fz, k, m, r))
            if series == "b" and n == 2:  # a2 = b2 - CO
                m = (mass - CO + chem.PROTON)
                if mz_range[0] <= m <= mz_range[1]:
                    rows.append(("a", 2, 1, 0, m, 0.3 * base))
    df = pd.DataFrame(rows, columns=["series", "index", "frag_z", "isotope", "mz", "rel"])
    if len(df):
        df["rel"] /= df["rel"].max()
    return df


def make_dia(sc: DIAScenario | None = None, seed: int = 0) -> SyntheticDIA:
    sc = sc or DIAScenario()
    rng = np.random.default_rng(seed)
    psc = sc.precursors
    lo, hi = sc.mz_span
    # peptides from a full sample; keep those with a charge state near the simulated span (fast, same density)
    comp = make_compounds(psc, rng)
    M = comp["neutral_mass"].to_numpy()
    near = np.zeros(len(comp), bool)
    for z in range(1, 6):
        mz = (M + z * chem.PROTON) / z
        near |= (mz > lo - 3) & (mz < hi + 1)
    comp = comp[near].reset_index(drop=True)
    comp["compound_id"] = np.arange(len(comp))
    p_sc = replace(psc, instrument=replace(psc.instrument, mz_range=(lo - 5, hi + 5)))
    prec = make_ions(comp, p_sc, rng)
    ms1 = simulate_run(comp, prec, p_sc, rng, name="ms1")

    edges = np.arange(lo, hi + 1e-9, sc.window_width)
    windows = pd.DataFrame(dict(window=np.arange(len(edges) - 1), lo=edges[:-1], hi=edges[1:]))
    windows.insert(1, "target", 0.5 * (windows.lo + windows.hi))
    # isolated share of each precursor charge state per window
    prec = prec.assign(window=np.searchsorted(edges, prec["mz"].to_numpy(), side="right") - 1)
    tot = prec.groupby(["compound_id", "charge"])["rel"].transform("sum")
    prec = prec.assign(share=prec["rel"] / tot, total=tot)
    iso = prec[(prec.window >= 0) & (prec.window < len(windows))]
    got = iso.groupby(["window", "compound_id", "charge"]).agg(share=("share", "sum"), total=("total", "first"))
    seqs = comp.set_index("compound_id")["sequence"]
    ms2_sc = replace(psc, instrument=sc.ms2, noise_points_per_scan=sc.ms2_noise_points)
    frag_cache = {}
    ms2 = []
    for w in range(len(windows)):
        rows = []
        if w in got.index.get_level_values(0):
            for (c, z), r in got.loc[w].iterrows():
                key = (c, z)
                if key not in frag_cache:
                    frag_cache[key] = fragments(seqs.loc[c], int(z), rng, sc.ms2.mz_range)
                f = frag_cache[key]
                # fragment apex SNR = gain * precursor apex SNR * isolated abundance * relative fragment intensity
                scale = sc.gain * (sc.ms2.noise / psc.instrument.noise) * r.share * r.total
                rows.append(f.assign(compound_id=c, prec_z=int(z), rel=f["rel"] * scale))
        ions = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
            columns=["series", "index", "frag_z", "isotope", "mz", "rel", "compound_id", "prec_z"])
        ions.insert(0, "ion_id", np.arange(len(ions)))
        ions["adduct"] = ions["series"] + ions["index"].astype(str) + "+" * 1
        ions["charge"] = ions["frag_z"]
        # windows are visited in turn: shift elution so the window's scans sit at its time in the cycle
        offset = psc.cycle * (w + 1) / (len(windows) + 1)
        cw = comp.assign(rt=comp["rt"] - offset)
        s = simulate_run(cw, ions, ms2_sc, rng, name=f"window{w}")
        s.run.rt = s.run.rt + offset
        s.run.meta["isolation"] = [float(windows.lo[w]), float(windows.hi[w])]
        ms2.append(replace(s, compounds=comp))
    return SyntheticDIA(ms1, windows, ms2, comp, sc)
