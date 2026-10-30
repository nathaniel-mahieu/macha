import numpy as np
import pandas as pd
import pytest

from macha import chem
from macha.annotate.complements import complement_mass, complement_null
from macha.annotate.fragment_groups import fragment_groups
from macha.annotate.ladder_eval import fragment_segment_truth, trace_peptide
from macha.annotate.ladders import PEPTIDE_LINE, fragment_decoys, fragment_rules, real_steps
from macha.io import DIARun
from macha.search import WindowIndex, digest, peptide_fdr, search_groups, with_decoys
from macha.synth import proteomics_scenario
from macha.synth.dia import RES, DIAScenario, make_dia
from macha.traces import detect_traces


def _frags(seq):
    r = np.array([RES[a] for a in seq])
    pre = np.cumsum(r)
    b = [pre[i - 1] + chem.PROTON for i in range(2, len(seq))]
    y = [pre[-1] - pre[i - 1] + chem.monoisotopic_mass({"H": 2, "O": 1}) + chem.PROTON for i in range(1, len(seq))]
    return np.array(b), np.array(y)


def test_residue_decoys_avoid_real_steps_and_match_mass_defects():
    rules, pairs = fragment_rules()
    dec, _ = fragment_decoys(rules, pairs)
    res = np.array([r.delta for r in rules if r.kind == "residue"])
    dres = np.array([r.delta for r in dec if r.kind == "residue"])
    assert len(dres) == 3 * len(res)
    real = real_steps(rules)
    assert min(np.abs(real - d).min() for d in dres) > 0.02
    dev = lambda m: m - np.round(m) * PEPTIDE_LINE  # noqa: E731
    assert abs(dev(res).mean() - dev(dres).mean()) < 1e-3 and abs(dev(res).std() - dev(dres).std()) < 1e-3


def test_complementary_pairs_give_precursor_mass():
    seq = "LVNELTEFAK"
    b, y = _frags(seq)
    M = sum(RES[a] for a in seq) + chem.monoisotopic_mass({"H": 2, "O": 1})
    window = ((M + 2 * chem.PROTON) / 2 - 5, (M + 2 * chem.PROTON) / 2 + 5)
    rng = np.random.default_rng(0)
    mz = np.concatenate([b[1::2], y[::2], rng.uniform(150, 1100, 20)])
    m, z, n = complement_mass(mz, window)
    assert z == 2 and abs(m - M) < 1e-3 and n >= 3
    assert complement_null(mz, window) < n


def test_digest_decoys_and_search_retrieve_a_planted_peptide():
    prot = [("p1", "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVK")]
    peps = digest(prot)
    assert "GVFRR" not in peps and "SEIAHRFK" in peps  # 1 missed cleavage, min length 7
    db = with_decoys(peps + ["LVNELTEFAK"])
    assert db.decoy.sum() > 0 and not set(db[db.decoy].peptide) & set(db[~db.decoy].peptide)
    target = "LVNELTEFAK"
    M = db.mass[db.peptide == target].iat[0]
    pm = (M + 2 * chem.PROTON) / 2
    ix = WindowIndex.build(db, (pm - 10, pm + 10))
    b, y = _frags(target)
    rng = np.random.default_rng(1)
    h = search_groups(ix, db, [np.concatenate([b[:3], y[:4], rng.uniform(150, 1500, 40)])])
    assert h.peptide.iat[0] == target and not h.decoy.iat[0]
    ids = peptide_fdr(pd.concat([h.assign(score=50.0)] + [h.assign(peptide=f"X{k}", decoy=False, score=40.0 - k)
                                                           for k in range(5)]), 0.01)
    assert target in set(ids.peptide)


@pytest.fixture(scope="module")
def small_dia():
    return make_dia(DIAScenario(precursors=proteomics_scenario(n_compounds=3000, n_background=0), window_width=20.0,
                                mz_span=(500.0, 540.0), ms2_noise_points=100), seed=4)


def test_dia_run_round_trip(small_dia, tmp_path):
    d = small_dia.as_dia_run()
    d.save(tmp_path / "dia")
    e = DIARun.load(tmp_path / "dia")
    assert len(e.ms2) == len(d.ms2) == 2
    assert e.ms2[1].n_points == d.ms2[1].n_points and np.allclose(e.ms2[1].rt, d.ms2[1].rt)
    assert e.window_of(510.0) == 0 and e.window_of(530.0) == 1 and e.window_of(600.0) == -1


def test_ladders_and_grouping_on_simulated_dia(small_dia):
    s = small_dia
    sr = s.ms2[0]
    tr = detect_traces(sr.run)
    fg = fragment_groups(tr, (float(s.windows.lo[0]), float(s.windows.hi[0])))
    seg = fg.graph.segments.join(fragment_segment_truth(fg.graph, sr.ion, sr.ions, s.compounds))
    res = seg[(seg.kind == "residue") & (seg.posterior >= 0.9)]
    assert len(res) > 100 and res.ladder.mean() > 0.9
    # ladder-seeded groups: mostly one peptide each
    pep = trace_peptide(tr, sr.ion, sr.ions)
    pur = [np.bincount(p[p >= 0]).max() / max((p >= 0).sum(), 1) for p in (pep[m] for m in fg.members) if (p >= 0).sum() >= 3]
    assert np.mean(pur) > 0.85
