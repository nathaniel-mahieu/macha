"""Elemental masses, formulas, isotopic fine structure and ion definitions.

Shared by the synthetic generator (to emit realistic ions) and by the
annotation module (to build mass-difference rules and expected isotope
ratios).
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

PROTON = 1.007276466812
ELECTRON = 0.00054857990946
NEUTRON_C13 = 1.0033548378  # 13C - 12C

# (mass, abundance) per isotope; IUPAC 2013 / AME2016
ISOTOPES: dict[str, list[tuple[float, float]]] = {
    "C": [(12.0, 0.9893), (13.0033548378, 0.0107)],
    "H": [(1.00782503207, 0.999885), (2.0141017778, 0.000115)],
    "N": [(14.0030740048, 0.99636), (15.0001088982, 0.00364)],
    "O": [(15.99491461956, 0.99757), (16.99913170, 0.00038), (17.9991610, 0.00205)],
    "P": [(30.97376163, 1.0)],
    "S": [(31.97207100, 0.9499), (32.97145876, 0.0075), (33.96786690, 0.0425), (35.96708076, 0.0001)],
    "Cl": [(34.96885268, 0.7576), (36.96590259, 0.2424)],
    "Br": [(78.9183371, 0.5069), (80.9162906, 0.4931)],
    "Na": [(22.9897692809, 1.0)],
    "K": [(38.96370668, 0.932581), (39.96399848, 0.000117), (40.96182576, 0.067302)],
    "F": [(18.99840322, 1.0)],
    "Si": [(27.9769265325, 0.92223), (28.976494700, 0.04685), (29.97377017, 0.03092)],
}
ELEMENTS = list(ISOTOPES)
MONO = {e: iso[int(np.argmax([a for _, a in iso]))][0] for e, iso in ISOTOPES.items()}

Formula = dict  # element -> count

_FORMULA_RE = re.compile(r"([A-Z][a-z]?)(\d*)")


def parse_formula(s: str) -> Formula:
    out: Formula = {}
    for el, n in _FORMULA_RE.findall(s):
        if el not in ISOTOPES:
            raise ValueError(f"unknown element {el}")
        out[el] = out.get(el, 0) + (int(n) if n else 1)
    return out


def format_formula(f: Formula) -> str:
    order = ["C", "H"] + sorted(e for e in f if e not in ("C", "H"))
    return "".join(f"{e}{f[e] if f[e] != 1 else ''}" for e in order if f.get(e, 0) > 0)


def add_formula(a: Formula, b: Formula, sign: int = 1) -> Formula:
    out = dict(a)
    for e, n in b.items():
        out[e] = out.get(e, 0) + sign * n
    if any(v < 0 for v in out.values()):
        raise ValueError(f"negative element count in {out}")
    return {e: n for e, n in out.items() if n}


def scale_formula(a: Formula, k: int) -> Formula:
    return {e: n * k for e, n in a.items()}


def monoisotopic_mass(f: Formula) -> float:
    return float(sum(MONO[e] * n for e, n in f.items()))


# --------------------------------------------------------------------------- isotopes
def _merge(masses, probs, decimals=6):
    key = np.round(masses, decimals)
    uniq, inv = np.unique(key, return_inverse=True)
    p = np.bincount(inv, weights=probs)
    m = np.bincount(inv, weights=masses * probs) / p
    return m, p


def _convolve(m1, p1, m2, p2, prune):
    m = (m1[:, None] + m2[None, :]).ravel()
    p = (p1[:, None] * p2[None, :]).ravel()
    keep = p >= prune * p.max()
    return _merge(m[keep], p[keep])


def _element_distribution(el: str, n: int, prune: float):
    iso = ISOTOPES[el]
    m1 = np.array([m for m, _ in iso])
    p1 = np.array([a for _, a in iso])
    rm, rp = np.array([0.0]), np.array([1.0])
    bm, bp = m1, p1
    while n:
        if n & 1:
            rm, rp = _convolve(rm, rp, bm, bp, prune)
        n >>= 1
        if n:
            bm, bp = _convolve(bm, bp, bm, bp, prune)
    return rm, rp


def isotope_fine_structure(f: Formula, prune: float = 1e-6):
    """Isotopologue masses and probabilities (fine structure, pruned)."""
    m, p = np.array([0.0]), np.array([1.0])
    for el, n in f.items():
        if n <= 0:
            continue
        em, ep = _element_distribution(el, n, prune)
        m, p = _convolve(m, p, em, ep, prune)
    order = np.argsort(m)
    return m[order], p[order] / p.sum()


def centroid_pattern(mz: np.ndarray, prob: np.ndarray, resolution_at_mz, merge: float = 1.0,
                     min_rel: float = 1e-4):
    """Merge fine-structure peaks an Orbitrap could not resolve.

    ``resolution_at_mz(mz)`` returns resolving power (m/dm FWHM).  Peaks closer
    than ``merge * FWHM`` are combined into one centroid at the
    intensity-weighted m/z.
    """
    keep = prob >= min_rel * prob.max()
    mz, prob = mz[keep], prob[keep]
    out_m, out_p = [], []
    cm, cp = mz[0], prob[0]
    for m, p in zip(mz[1:], prob[1:]):
        fwhm = cm / resolution_at_mz(cm)
        if m - cm < merge * fwhm:
            cm = (cm * cp + m * p) / (cp + p)
            cp += p
        else:
            out_m.append(cm)
            out_p.append(cp)
            cm, cp = m, p
    out_m.append(cm)
    out_p.append(cp)
    return np.array(out_m), np.array(out_p)


def nominal_isotope_index(mz: np.ndarray, mono_mz: float, z: int) -> np.ndarray:
    return np.rint((mz - mono_mz) * abs(z) / NEUTRON_C13).astype(int)


# --------------------------------------------------------------------------- ions
@dataclass(frozen=True)
class Adduct:
    name: str
    n_mol: int
    add: str
    remove: str
    charge: int

    def ion_formula(self, f: Formula) -> Formula:
        out = scale_formula(f, self.n_mol)
        if self.add:
            out = add_formula(out, parse_formula(self.add))
        if self.remove:
            out = add_formula(out, parse_formula(self.remove), -1)
        return out

    def mz(self, neutral_mass: float) -> float:
        add = monoisotopic_mass(parse_formula(self.add)) if self.add else 0.0
        rem = monoisotopic_mass(parse_formula(self.remove)) if self.remove else 0.0
        return (self.n_mol * neutral_mass + add - rem - self.charge * ELECTRON) / abs(self.charge)

    def neutral_mass(self, mz: float) -> float:
        add = monoisotopic_mass(parse_formula(self.add)) if self.add else 0.0
        rem = monoisotopic_mass(parse_formula(self.remove)) if self.remove else 0.0
        return (mz * abs(self.charge) + self.charge * ELECTRON - add + rem) / self.n_mol


ADDUCTS = {a.name: a for a in [
    Adduct("[M+H]+", 1, "H", "", 1),
    Adduct("[M+Na]+", 1, "Na", "", 1),
    Adduct("[M+K]+", 1, "K", "", 1),
    Adduct("[M+NH4]+", 1, "NH4", "", 1),
    Adduct("[M+H-H2O]+", 1, "H", "H2O", 1),
    Adduct("[M+H-NH3]+", 1, "H", "NH3", 1),
    Adduct("[M+2H]2+", 1, "H2", "", 2),
    Adduct("[M+H+Na]2+", 1, "HNa", "", 2),
    Adduct("[M+3H]3+", 1, "H3", "", 3),
    Adduct("[M+4H]4+", 1, "H4", "", 4),
    Adduct("[M+5H]5+", 1, "H5", "", 5),
    Adduct("[2M+H]+", 2, "H", "", 1),
    Adduct("[2M+Na]+", 2, "Na", "", 1),
    Adduct("[M-H]-", 1, "", "H", -1),
    Adduct("[M+Cl]-", 1, "Cl", "", -1),
    Adduct("[M+HCOO]-", 1, "HCO2", "", -1),
    Adduct("[M-H-H2O]-", 1, "", "H3O", -1),
    Adduct("[M-2H]2-", 1, "", "H2", -2),
    Adduct("[2M-H]-", 2, "", "H", -1),
]}


# --------------------------------------------------------------------------- generators
def random_metabolite_formula(rng: np.random.Generator) -> Formula:
    """A plausible small-molecule CHNOPS(+halogen) formula (nitrogen rule respected)."""
    while True:
        c = int(np.clip(rng.gamma(4.0, 4.0) + 2, 2, 60))
        h = int(round(c * rng.uniform(0.9, 2.1)))
        n = int(min(rng.poisson(c * 0.08 + 0.3), c // 2 + 1))
        o = int(min(rng.poisson(c * 0.32 + 0.5), c + 2))
        f = {"C": c, "H": h, "N": n, "O": o}
        if rng.random() < 0.10:
            f["S"] = int(rng.integers(1, 3))
        if rng.random() < 0.07:
            f["P"] = int(rng.integers(1, 3))
            f["O"] = f["O"] + 3 * f["P"]
        if rng.random() < 0.04:
            f["Cl"] = int(rng.integers(1, 3))
        if rng.random() < 0.01:
            f["Br"] = 1
        odd = f["H"] + f["N"] + f.get("P", 0) + f.get("Cl", 0) + f.get("Br", 0)
        if odd % 2:
            f["H"] += 1
        hal = f.get("Cl", 0) + f.get("Br", 0)
        rdbe = c - (f["H"] + hal) / 2 + (n + f.get("P", 0)) / 2 + 1
        if 0 <= rdbe <= c and f["H"] > 0:
            mass = monoisotopic_mass(f)
            if 60 < mass < 1200:
                return {k: v for k, v in f.items() if v}


AA_FORMULA = {
    "G": "C2H3NO", "A": "C3H5NO", "S": "C3H5NO2", "P": "C5H7NO", "V": "C5H9NO", "T": "C4H7NO2",
    "C": "C5H8N2O2S",  # carbamidomethyl
    "L": "C6H11NO", "I": "C6H11NO", "N": "C4H6N2O2", "D": "C4H5NO3", "Q": "C5H8N2O2",
    "K": "C6H12N2O", "E": "C5H7NO3", "M": "C5H9NOS", "H": "C6H7N3O", "F": "C9H9NO",
    "R": "C6H12N4O", "Y": "C9H9NO2", "W": "C11H10N2O",
}
AA_FREQ = {"A": 7.0, "R": 5.6, "N": 3.6, "D": 4.7, "C": 2.3, "Q": 4.7, "E": 7.1, "G": 6.6, "H": 2.6,
           "I": 4.3, "L": 10.0, "K": 5.7, "M": 2.1, "F": 3.6, "P": 6.3, "S": 8.3, "T": 5.4, "W": 1.2,
           "Y": 2.7, "V": 6.0}
# Guo et al. 1986 retention coefficients (RP, pH 2), used for a crude RT model
AA_HYDRO = {"W": 8.8, "F": 8.1, "L": 8.1, "I": 7.4, "M": 5.5, "V": 5.0, "Y": 4.5, "C": 2.6, "P": 2.0,
            "A": 2.0, "E": 1.1, "T": 0.6, "D": 0.2, "Q": 0.0, "S": -0.2, "G": -0.2, "R": -0.6,
            "N": -0.6, "H": -2.1, "K": -2.1}


def random_peptide(rng: np.random.Generator) -> str:
    aas = list(AA_FREQ)
    p = np.array([AA_FREQ[a] for a in aas])
    p /= p.sum()
    length = int(np.clip(rng.gamma(5.0, 2.6), 6, 30))
    body = "".join(rng.choice(aas, size=length - 1, p=p))
    body = body.replace("K", "Q").replace("R", "N")  # tryptic: basic residue only at C-term
    return body + ("K" if rng.random() < 0.55 else "R")


def peptide_formula(seq: str) -> Formula:
    f = {"H": 2, "O": 1}
    for aa in seq:
        f = add_formula(f, parse_formula(AA_FORMULA[aa]))
    return f


def peptide_hydrophobicity(seq: str) -> float:
    return float(sum(AA_HYDRO[a] for a in seq))


def expected_isotope_ratio(mass: float, n: int = 1, kind: str = "metabolite") -> float:
    """Expected M+n / M intensity ratio from an averaged composition model.

    Uses carbon count estimates: ~mass/14.5 for metabolites, averagine
    (4.9384 C per 111.1254 Da) for peptides; Poisson approximation.
    """
    if kind == "peptide":
        n_c = mass / 111.1254 * 4.9384
        lam = n_c * 0.0107 + mass / 111.1254 * (1.3577 * 0.00364 + 7.7583 * 0.000115 + 0.0417 * 0.0075)
    else:
        n_c = mass / 14.5 * 0.85
        lam = n_c * 0.0107 + 0.0015
    from math import factorial

    return lam**n / factorial(n)
