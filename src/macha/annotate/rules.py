"""Mass-difference rules and decoys.

Additive rules relate features of equal charge ``z`` through their
charge-scaled masses ``M = mz * z``: ``M_j - M_i = n * delta``.  Pair rules
relate two adduct hypotheses of one neutral molecule (different charge or
multimer), where the relation is not a constant mass difference.
"""

from __future__ import annotations

from dataclasses import dataclass, replace


from .. import chem

H = chem.MONO["H"]
N = chem.MONO["N"]
O = chem.MONO["O"]
C = chem.MONO["C"]


@dataclass(frozen=True)
class Rule:
    name: str
    delta: float  # Da in charge-scaled mass units
    kind: str  # isotope | adduct | loss | series
    n_max: int = 1
    z_max: int = 1  # largest charge for which the rule is searched
    coelute: bool = True
    sigma_extra: float = 0.0  # Da per step, e.g. unresolved isotopic fine structure
    decoy: bool = False
    target: str = ""  # for decoys: the rule they mimic

    def with_z(self, z_max: int) -> "Rule":
        return replace(self, z_max=z_max)


@dataclass(frozen=True)
class PairRule:
    """Two adduct hypotheses of the same neutral molecule (charge states, multimers)."""

    name: str
    a: str  # adduct name of the lighter-hypothesis feature
    b: str
    decoy_shift: float = 0.0  # decoys offset neutral mass agreement
    decoy: bool = False
    target: str = ""


def isotope_rules(z_max: int = 4, n_max: int = 4, sigma_extra: float = 6e-4) -> list[Rule]:
    iso = chem.ISOTOPES
    return [
        Rule("13C", chem.NEUTRON_C13, "isotope", n_max, z_max, True, sigma_extra),
        Rule("15N", iso["N"][1][0] - iso["N"][0][0], "isotope", 1, z_max),
        Rule("34S", iso["S"][2][0] - iso["S"][0][0], "isotope", 1, z_max),
        Rule("37Cl", iso["Cl"][1][0] - iso["Cl"][0][0], "isotope", 2, min(z_max, 2)),
        Rule("81Br", iso["Br"][1][0] - iso["Br"][0][0], "isotope", 1, min(z_max, 2)),
        Rule("18O", iso["O"][2][0] - iso["O"][0][0], "isotope", 1, min(z_max, 2)),
        Rule("41K", iso["K"][2][0] - iso["K"][0][0], "isotope", 1, 1),
    ]


def adduct_rules(polarity: int = 1, z_max: int = 2) -> list[Rule]:
    na, k, cl = chem.MONO["Na"], chem.MONO["K"], chem.MONO["Cl"]
    common = [
        Rule("H2O", 2 * H + O, "loss", 1, z_max),
        Rule("NH3", N + 3 * H, "loss", 1, z_max),
        Rule("CO2", C + 2 * O, "loss", 1, 1),
        Rule("CO", C + O, "loss", 1, 1),
        Rule("HCOOH", C + 2 * H + 2 * O, "adduct", 1, 1),
    ]
    if polarity > 0:
        return common + [
            Rule("Na-H", na - H, "adduct", 1, z_max),
            Rule("K-H", k - H, "adduct", 1, z_max),
            Rule("K-Na", k - na, "adduct", 1, 1),
            Rule("Na-NH4", na - N - 4 * H, "adduct", 1, 1),
        ]
    return common + [Rule("Cl+H", cl + H, "adduct", 1, 1)]


def series_rules() -> list[Rule]:
    return [Rule("CH2", C + 2 * H, "series", 12, 1, False), Rule("C2H4O", 2 * C + 4 * H + O, "series", 12, 1, False)]


def pair_rules(polarity: int = 1, peptides: bool = False) -> list[PairRule]:
    if polarity > 0:
        out = [PairRule("[M+H]/[2M+H]", "[M+H]+", "[2M+H]+"), PairRule("[M+Na]/[2M+Na]", "[M+Na]+", "[2M+Na]+"),
               PairRule("[M+H]/[M+2H]", "[M+H]+", "[M+2H]2+")]
        if peptides:
            out += [PairRule("[M+2H]/[M+3H]", "[M+2H]2+", "[M+3H]3+"),
                    PairRule("[M+3H]/[M+4H]", "[M+3H]3+", "[M+4H]4+"),
                    PairRule("[M+2H]/[M+4H]", "[M+2H]2+", "[M+4H]4+")]
        return out
    return [PairRule("[M-H]/[2M-H]", "[M-H]-", "[2M-H]-"), PairRule("[M-H]/[M-2H]", "[M-H]-", "[M-2H]2-")]


def default_rules(polarity: int = 1, peptides: bool = False, series: bool = False):
    z_iso = 5 if peptides else 3
    rules = isotope_rules(z_max=z_iso, n_max=5 if peptides else 4) + adduct_rules(polarity, z_max=1 if peptides else 2)
    if series:
        rules += series_rules()
    return rules, pair_rules(polarity, peptides)


# ------------------------------------------------------------------ decoys
DECOY_SHIFTS = (0.0731, -0.0857, 0.1173, -0.1291, 0.1613, -0.1777, 0.2311, -0.2459)


def _collides(delta, rules, tol=0.004, z_max=5, n_max=6):
    for r in rules:
        for z in range(1, z_max + 1):
            for n in range(1, n_max + 1):
                if abs(delta - n * r.delta / z) < tol:
                    return True
    # sums of an isotope step and an adduct difference
    return False


def make_decoys(rules: list[Rule], pair_rules_: list[PairRule], n_per_rule: int = 6):
    """Decoy rules: deltas shifted by amounts that match no rule (or multiple / z of one)."""
    decoys = []
    for r in rules:
        k = 0
        for s in DECOY_SHIFTS:
            d = r.delta + s
            if d <= 0.3 or _collides(d, rules):
                continue
            decoys.append(replace(r, name=f"decoy{k}:{r.name}", delta=d, decoy=True, target=r.name))
            k += 1
            if k >= n_per_rule:
                break
    pdecoys = []
    for p in pair_rules_:
        for k, s in enumerate(DECOY_SHIFTS[:n_per_rule]):
            pdecoys.append(PairRule(f"decoy{k}:{p.name}", p.a, p.b, decoy_shift=s * 7, decoy=True, target=p.name))
    return decoys, pdecoys


def affine(adduct_name: str):
    """(A, B, n) such that the neutral mass is (mz * A + B) / n for this adduct."""
    a = chem.ADDUCTS[adduct_name]
    add = chem.monoisotopic_mass(chem.parse_formula(a.add)) if a.add else 0.0
    rem = chem.monoisotopic_mass(chem.parse_formula(a.remove)) if a.remove else 0.0
    return float(abs(a.charge)), float(a.charge * chem.ELECTRON - add + rem), float(a.n_mol)


DEFAULT_ADDUCT_PRIOR_POS = {"[M+H]+": 1.0, "[M+Na]+": 0.5, "[M+NH4]+": 0.35, "[M+K]+": 0.2, "[M+H-H2O]+": 0.25,
                            "[M+H-NH3]+": 0.1, "[2M+H]+": 0.1, "[2M+Na]+": 0.05, "[M+2H]2+": 0.3,
                            "[M+3H]3+": 0.25, "[M+4H]4+": 0.1, "[M+5H]5+": 0.05, "[M+H+Na]2+": 0.05}
DEFAULT_ADDUCT_PRIOR_NEG = {"[M-H]-": 1.0, "[M+HCOO]-": 0.4, "[M+Cl]-": 0.2, "[M-H-H2O]-": 0.15, "[2M-H]-": 0.1,
                            "[M-2H]2-": 0.05}


def adduct_prior(polarity: int) -> dict:
    return DEFAULT_ADDUCT_PRIOR_POS if polarity > 0 else DEFAULT_ADDUCT_PRIOR_NEG


__all__ = ["Rule", "PairRule", "default_rules", "isotope_rules", "adduct_rules", "series_rules", "pair_rules",
           "make_decoys", "affine", "adduct_prior"]
