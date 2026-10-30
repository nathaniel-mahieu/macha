"""Peptide identification of DIA fragment groups: one search for every grouping method.

Each group (a set of fragment m/z values, optionally with a precursor mass) is matched against the tryptic
peptides whose precursor fits the group's isolation window.  The score is the binomial tail of the number of
matched b / y ions given the chance-match rate implied by the group's own fragment density, so large groups do not
win by size.  Pseudo-reversed decoys compete with targets, and identifications are reported at a target-decoy
false-discovery rate on the peptide level.
"""

from __future__ import annotations

import gzip
import math
from dataclasses import dataclass
from pathlib import Path

import numba as nb
import numpy as np
import pandas as pd

from . import chem

AA = "ACDEFGHIKLMNPQRSTVWY"
RES = np.zeros(128)
for _a in AA:
    RES[ord(_a)] = chem.monoisotopic_mass(chem.parse_formula(chem.AA_FORMULA[_a]))
H2O = chem.monoisotopic_mass({"H": 2, "O": 1})


def read_fasta(*paths) -> list[tuple[str, str]]:
    out = []
    for p in paths:
        p = Path(p)
        op = gzip.open if p.suffix == ".gz" else open
        name, seq = None, []
        with op(p, "rt") as f:
            for line in f:
                if line.startswith(">"):
                    if name is not None:
                        out.append((name, "".join(seq)))
                    name, seq = line[1:].split()[0], []
                else:
                    seq.append(line.strip())
        if name is not None:
            out.append((name, "".join(seq)))
    return out


def digest(proteins, missed: int = 1, min_len: int = 7, max_len: int = 30) -> list[str]:
    """Unique tryptic peptides (cleave after K / R, not before P), standard residues only."""
    peps = set()
    allowed = set(AA)
    for _, s in proteins:
        cuts = [0] + [i + 1 for i in range(len(s) - 1) if s[i] in "KR" and s[i + 1] != "P"] + [len(s)]
        for a in range(len(cuts) - 1):
            for m in range(missed + 1):
                b = a + 1 + m
                if b >= len(cuts):
                    break
                p = s[cuts[a]:cuts[b]]
                if min_len <= len(p) <= max_len and set(p) <= allowed:
                    peps.add(p)
    return sorted(peps)


def with_decoys(peptides: list[str]) -> pd.DataFrame:
    """Targets plus pseudo-reversed decoys (body reversed, C-terminal residue kept); I and L are one residue."""
    t = pd.Series(peptides)
    d = t.map(lambda p: p[-2::-1] + p[-1])
    key = lambda s: s.str.replace("I", "L")  # noqa: E731
    d = d[~key(d).isin(set(key(t)))]
    df = pd.DataFrame(dict(peptide=pd.concat([t, d], ignore_index=True),
                           decoy=np.r_[np.zeros(len(t), bool), np.ones(len(d), bool)]))
    df["mass"] = [RES[np.frombuffer(p.encode(), np.uint8)].sum() + H2O for p in df.peptide]
    return df


@dataclass
class WindowIndex:
    """Fragment index of the candidates (peptide, charge) whose precursor fits one isolation window."""

    cand_pep: np.ndarray
    cand_z: np.ndarray
    cand_mass: np.ndarray
    cand_n: np.ndarray  # theoretical fragments per candidate
    frag_mz: np.ndarray  # sorted
    frag_cand: np.ndarray
    span: float

    @classmethod
    def build(cls, db: pd.DataFrame, window: tuple, charges=(2, 3, 4), iso_slack: int = 2,
              mz_range=(150.0, 1600.0)) -> "WindowIndex":
        lo, hi = window
        M = db["mass"].to_numpy()
        pep, zs = [], []
        for z in charges:
            mz = (M + z * chem.PROTON) / z
            ok = np.flatnonzero((mz >= lo - iso_slack * chem.NEUTRON_C13 / z) & (mz < hi))
            pep.append(ok)
            zs.append(np.full(ok.size, z))
        pep, zs = np.concatenate(pep), np.concatenate(zs)
        seqs = db["peptide"].to_numpy()[pep]
        L = np.array([len(s) for s in seqs])
        Lmax = int(L.max()) if L.size else 1
        R = np.zeros((len(seqs), Lmax))
        for k, s in enumerate(seqs):
            R[k, : len(s)] = RES[np.frombuffer(s.encode(), np.uint8)]
        pre = np.cumsum(R, 1)  # b_i neutral residue sums
        tot = pre[np.arange(len(seqs)), L - 1]
        f_mz, f_c = [], []
        for i in range(1, Lmax):
            ok = L > i  # cleavage after residue i
            b = pre[ok, i - 1]
            y = tot[ok] - b + H2O
            c = np.flatnonzero(ok)
            for mass, keep in ((b, i >= 2), (y, True)):
                if not keep:
                    continue
                f_mz.append((mass + chem.PROTON))
                f_c.append(c)
                two = zs[c] >= 3
                f_mz.append(((mass + 2 * chem.PROTON) / 2)[two])
                f_c.append(c[two])
        f_mz, f_c = np.concatenate(f_mz), np.concatenate(f_c)
        keep = (f_mz >= mz_range[0]) & (f_mz <= mz_range[1])
        f_mz, f_c = f_mz[keep], f_c[keep]
        o = np.argsort(f_mz, kind="stable")
        return cls(pep, zs, db["mass"].to_numpy()[pep], np.bincount(f_c, minlength=len(seqs)).astype(np.int64),
                   f_mz[o], f_c[o].astype(np.int64), float(mz_range[1] - mz_range[0]))


@nb.njit(cache=True)
def _log10_binom_tail(k, n, p):
    """log10 P(X >= k), X ~ Binomial(n, p)."""
    if k <= 0:
        return 0.0
    lp, lq = math.log(p), math.log1p(-p)
    m = -1e300
    terms = np.empty(n - k + 1)
    for i in range(k, n + 1):
        t = math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * lp + (n - i) * lq
        terms[i - k] = t
        if t > m:
            m = t
    s = 0.0
    for t in terms:
        s += math.exp(t - m)
    return (m + math.log(s)) / math.log(10.0)


@nb.njit(cache=True)
def _search(g_offs, g_mz, g_mass, frag_mz, frag_cand, cand_n, cand_mass, tol_ppm, prec_ppm, span):
    n_g = g_offs.shape[0] - 1
    n_c = cand_n.shape[0]
    counts = np.zeros(n_c, np.int64)
    touched = np.empty(n_c, np.int64)
    best = np.full(n_g, -1, np.int64)
    score = np.zeros(n_g)
    kk = np.zeros(n_g, np.int64)
    for g in range(n_g):
        a, b = g_offs[g], g_offs[g + 1]
        nt = 0
        msum = 0.0
        for q in range(a, b):
            m = g_mz[q]
            msum += m
            tol = m * tol_ppm * 1e-6
            lo = np.searchsorted(frag_mz, m - tol)
            j = lo
            while j < frag_mz.shape[0] and frag_mz[j] <= m + tol:
                c = frag_cand[j]
                if counts[c] == 0:
                    touched[nt] = c
                    nt += 1
                counts[c] += 1
                j += 1
        nf = b - a
        if nf == 0:
            continue
        p = min(0.5, nf * 2.0 * tol_ppm * 1e-6 * (msum / nf) / span)
        gm = g_mass[g]
        for t in range(nt):
            c = touched[t]
            k = counts[c]
            counts[c] = 0
            if gm == gm:  # precursor mass known: allow -1 / 0 / +1 isotope errors
                ok = False
                for e in (-1.0, 0.0, 1.0):
                    if abs(cand_mass[c] + e * 1.0033548378 - gm) <= prec_ppm * 1e-6 * gm:
                        ok = True
                if not ok:
                    continue
            n = cand_n[c]
            if k > n:
                k = n
            s = -_log10_binom_tail(k, n, p)
            if s > score[g]:
                score[g] = s
                best[g] = c
                kk[g] = k
    return best, score, kk


def search_groups(index: WindowIndex, db: pd.DataFrame, groups: list, masses=None, tol_ppm: float = 10.0,
                  prec_ppm: float = 10.0) -> pd.DataFrame:
    """Best candidate per group.  ``groups`` is a list of fragment m/z arrays; ``masses`` optional precursor masses."""
    lens = np.array([len(g) for g in groups], np.int64)
    offs = np.concatenate([[0], np.cumsum(lens)])
    mz = np.concatenate([np.asarray(g, float) for g in groups]) if len(groups) else np.zeros(0)
    gm = np.full(len(groups), np.nan) if masses is None else np.asarray(masses, float)
    best, score, k = _search(offs, mz, gm, index.frag_mz, index.frag_cand, index.cand_n, index.cand_mass,
                             tol_ppm, prec_ppm, index.span)
    ok = best >= 0
    pep = np.where(ok, index.cand_pep[np.maximum(best, 0)], -1)
    return pd.DataFrame(dict(group=np.arange(len(groups)), score=score, matched=k, n_fragments=lens,
                             peptide=np.where(ok, db["peptide"].to_numpy()[np.maximum(pep, 0)], None),
                             decoy=np.where(ok, db["decoy"].to_numpy()[np.maximum(pep, 0)], False),
                             charge=np.where(ok, index.cand_z[np.maximum(best, 0)], 0)))[ok].reset_index(drop=True)


def peptide_fdr(hits: pd.DataFrame, q: float = 0.01, conservative: bool = False) -> pd.DataFrame:
    """Peptide-level target-decoy q-values (best score per peptide); returns targets with q <= ``q``.

    FDR at a threshold is decoys / targets above it (target-decoy competition); ``conservative`` uses
    (decoys + 1) / targets, which cannot certify fewer than 1/q identifications.
    """
    if not len(hits):
        return hits.assign(q=[])
    h = hits.assign(key=hits.peptide.str.replace("I", "L"))
    best = h.sort_values("score", ascending=False).drop_duplicates("key")
    d = best.decoy.to_numpy()
    fdr = (np.cumsum(d) + int(conservative)) / np.maximum(np.cumsum(~d), 1)
    best = best.assign(q=np.minimum.accumulate(fdr[::-1])[::-1])
    return best[(~best.decoy) & (best.q <= q)]
