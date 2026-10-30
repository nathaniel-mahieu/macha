"""Fragment-ladder rules for MS2 (DIA) data: amino-acid residue steps as modulo relations.

Consecutive b (or y) ions of one precursor differ by one residue mass, so a fragment ion series is a chain of
residue relations, and the scan-level graph (:func:`macha.annotate.scan_graph.scan_graph`) can test them between
fragment centroids exactly as it tests isotopes and adducts between precursors.

Decoys need care here.  Every fragment is a peptide fragment, so the differences between *unrelated* fragments
also cluster near the peptide mass-defect line (mass ~ 1.000495 x nominal).  Shifting a residue mass by a fraction
of a dalton (the MS1 decoy scheme) moves it off that line, where random differences are rare, and would make the
false-discovery estimate optimistic.  Residue decoys are therefore placed on the line, more than 0.02 Da from
every residue, residue pair, isotope or neutral-loss step, each keeping a residue's own deviation from the line.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from .. import chem
from .rules import PairRule, Rule, isotope_rules

PEPTIDE_LINE = 1.000495  # average monoisotopic mass per nominal mass unit of peptide residues


def residue_masses(merge_il: bool = True) -> dict:
    """Monoisotopic residue masses (carbamidomethyl C), with I/L merged into one isobaric step."""
    out = {}
    for aa, f in chem.AA_FORMULA.items():
        if merge_il and aa == "I":
            continue
        out["L/I" if (merge_il and aa == "L") else aa] = chem.monoisotopic_mass(chem.parse_formula(f))
    return out


def residue_rules(z_max: int = 2) -> list[Rule]:
    return [Rule(name, m, "residue", 1, z_max) for name, m in sorted(residue_masses().items(), key=lambda x: x[1])]


def loss_rules(z_max: int = 2) -> list[Rule]:
    m = lambda f: chem.monoisotopic_mass(chem.parse_formula(f))  # noqa: E731
    return [Rule("H2O", m("H2O"), "loss", 1, z_max), Rule("NH3", m("NH3"), "loss", 1, z_max),
            Rule("CO", m("CO"), "loss", 1, z_max)]


def fragment_rules(z_max: int = 2) -> tuple[list[Rule], list[PairRule]]:
    """Residue steps, fragment isotopes (13C) and common losses; plus the 1+/2+ fragment charge-state pair."""
    rules = isotope_rules(z_max=z_max, n_max=2)[:1] + residue_rules(z_max) + loss_rules(z_max)
    return rules, [PairRule("[M+H]/[M+2H]", "[M+H]+", "[M+2H]2+")]


def real_steps(rules: list[Rule]) -> np.ndarray:
    """Mass differences real fragment pairs produce: residues, residue pairs, each +/- isotope and loss steps."""
    res = [r.delta for r in rules if r.kind == "residue"]
    steps = [0.0] + [r.delta for r in rules if r.kind == "loss"] + [chem.NEUTRON_C13, 2 * chem.NEUTRON_C13]
    base = np.array(res + [a + b for a in res for b in res])
    return np.unique(np.concatenate([base + s for s in steps] + [base - s for s in steps[1:]]))


def residue_decoys(rules: list[Rule], n_per_rule: int = 3, mass_range=(50, 200), tol: float = 0.02,
                   seed: int = 7) -> list[Rule]:
    """Decoy residue rules on the peptide mass-defect line, at exact masses no real step comes close to.

    Each decoy keeps the deviation of one residue from the line (so coincidences between unrelated fragments
    hit it as often as they hit that residue) and sits more than ``tol`` Da from every real step.
    """
    rng = np.random.default_rng(seed)
    real = real_steps(rules)
    res = [r for r in rules if r.kind == "residue"]
    taken: list[float] = []
    out = []
    for k in range(n_per_rule):
        for r in res:
            dev = r.delta - round(r.delta) * PEPTIDE_LINE
            for n in rng.permutation(np.arange(*mass_range)):
                d = float(n * PEPTIDE_LINE + dev)
                if np.abs(real - d).min() > tol and all(abs(d - t) > 0.005 for t in taken):
                    break
            else:
                raise ValueError(f"no free decoy mass for {r.name}")
            taken.append(d)
            out.append(replace(r, name=f"decoy{k}:{r.name}", delta=d, decoy=True, target=r.name))
    return out


def fragment_decoys(rules: list[Rule], pair_rules: list[PairRule], n_per_rule: int = 3):
    """Defect-matched residue decoys plus the usual shifted decoys for isotope / loss / charge-state rules."""
    from .rules import make_decoys

    other, pdecoys = make_decoys([r for r in rules if r.kind != "residue"], pair_rules, n_per_rule)
    return residue_decoys(rules, n_per_rule) + other, pdecoys
