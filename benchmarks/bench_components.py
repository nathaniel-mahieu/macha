"""Components: peak models of scan-graph relation clusters, per run and across runs (macha.components).

1. Deconvolution (synthetic, one run): near-isobaric pairs 0.5-6 ppm apart (closer than the resolution, so their
   centroids coalesce) eluting within 2.5 peak sd.  A pair is resolved when two different components (or peaks)
   on its main ion have apexes within 0.5 sd of both true apexes.  Picked peaks, components from intensity alone,
   and components with the m/z trajectory.  False splits on eluting compounds with one true source.
2. Cross-run (synthetic batch, true heights per run): linking purity, fragments per compound, coverage before and
   after targeted re-query, and the error of centred log abundances against truth, for components, relation
   groups and picked-peak feature groups on the same compounds.
3. Real data: metabolomics (Sa04 injected three times: replicate spread without truth) and ProteoBench M7
   (mixed-species design: log2 A/B at known offsets), per-run split statistics, and the m/z difference the
   deconvolution assigns to coalesced splits.
"""

from __future__ import annotations

import pickle
import sys
from types import SimpleNamespace

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from bench_xrun import SEG_COLS, SPECIES, _mix, log2_ratio, real_files
from common import CACHE, NAV, REAL_LABEL, REPORTS, Timer, save_json
from macha import MSRun
from macha.align import align
from macha.align.drift import cv_drift
from macha.align.edges import drift_window, relation_anchors, relation_groups, relation_noise, relation_table
from macha.annotate.scan_graph import scan_graph
from macha.baseline import trace_baselines
from macha.components import fit_components, link_components, predict_apex, quantify
from macha.peaks import detect_peaks
from macha.report import Page, method_color, style
from macha.synth import make_batch, make_run, metabolomics_scenario
from macha.synth.lc import emg_mode
from macha.synth.truth import feature_truth
from macha.traces import Traces, detect_traces

DECONV = dict(n_compounds=300, n_isobaric_pairs=120, isobaric_ppm=(0.5, 6.0), n_background=20)
BATCH = dict(sc=dict(n_compounds=300, n_isobaric_pairs=40, isobaric_ppm=(0.5, 6.0), n_background=20),
             n_samples=6, n_drift_classes=2, class_drift=10.0, global_drift=5.0, seed=4)
METHODS = {"peaks": "picked peaks", "int": "components, intensity only", "mz": "components, intensity + m/z"}


# ------------------------------------------------------------------------------------------------ synthetic truth
def _pairs(compounds, mode):
    cp = compounds.reset_index(drop=True)
    cop = cp[cp.mass_offset_ppm != 0]
    src = cp[(cp.mass_offset_ppm == 0) & (cp.kind != "background")].drop_duplicates("formula").set_index("formula")
    P = pd.DataFrame(dict(c=cop.compound_id.to_numpy(), s=src.loc[cop.formula, "compound_id"].to_numpy()))
    P["dm_ppm"] = np.abs(cop.mass_offset_ppm.to_numpy())
    P["sig"] = cp.sigma.to_numpy()[P.s]
    P["dt_sd"] = np.abs(mode[P.c] - mode[P.s]) / P.sig
    P["log10_ratio"] = np.log10(cp.height.to_numpy()[P.c] / cp.height.to_numpy()[P.s])
    return P


def _resolved(P, main_mz, mode, sigma, rt, mz, ident, tol_ppm=10.0):
    """Two different entries (components or peaks) on the pair's main ion, each within 0.5 sd of one true apex."""
    out = []
    for s, c in zip(P.s.to_numpy(), P.c.to_numpy()):
        m = main_mz[s]
        near = np.abs(mz - m) / m * 1e6 < tol_ppm
        hits = [set(ident[near & (np.abs(rt - mode[k]) <= max(0.5 * sigma[k], 1.0))].tolist()) for k in (s, c)]
        out.append(bool(hits[0]) and bool(hits[1]) and len(hits[0] | hits[1]) >= 2)
    return np.array(out)


def _purity(C, tr, label, merged, pair_set, kind):
    """Per compound: true sources of its ions' signal in its window (>= 10 % of the largest)."""
    run = tr.run
    sc = run.scan_of_point
    win = C.fits.set_index("compound")
    rows = []
    for comp, g in C.ions[C.ions.k == 0].groupby("compound"):
        lo, hi = win.loc[comp, "rt_lo"], win.loc[comp, "rt_hi"]
        pts = np.concatenate([tr.points(int(t)) for t in g.trace])
        r = run.rt[sc[pts]]
        pts = pts[(r >= lo) & (r <= hi)]
        s = pd.Series(run.intensity[pts]).groupby(label[pts]).sum()
        s = s[s.index >= 0]
        top = int(s.idxmax()) if len(s) else -1
        rows.append(dict(compound=comp, n_true=int((s >= 0.1 * s.max()).sum()) if len(s) else 0, top=top,
                         coalesced_frac=float(merged[pts].mean()) if pts.size else 0.0))
    P = pd.DataFrame(rows).merge(C.fits, on="compound")
    P["background"] = kind[np.maximum(P.top, 0)] == "background"
    P["clean"] = (P.n_true == 1) & ~P.top.isin(pair_set) & (P.coalesced_frac < 0.02) & ~P.background
    return P


def deconvolution():
    sr = make_run(metabolomics_scenario(**DECONV), seed=2)
    with Timer() as t_tr:
        tr = detect_traces(sr.run)
    g = scan_graph(tr)
    rel = relation_table(g, tr)
    cp = sr.compounds.reset_index(drop=True)
    mode = emg_mode(cp.rt.to_numpy(), cp.sigma.to_numpy(), np.maximum(cp.tau.to_numpy(), 1e-3))
    P = _pairs(cp, mode)
    main = sr.ions.sort_values("rel", ascending=False).drop_duplicates("compound_id").set_index("compound_id")["mz"]
    main_mz = main.reindex(np.arange(len(cp))).to_numpy()
    sigma = cp.sigma.to_numpy()
    feats = detect_peaks(tr, trace_baselines(tr))
    P["peaks"] = _resolved(P, main_mz, mode, sigma, feats.rt.to_numpy(), feats.mz.to_numpy(), np.arange(len(feats)))
    kind = cp.kind.to_numpy()
    pair_set = set(P.s) | set(P.c)
    fp, comps, timing = [], {}, {}
    for key, kw in (("int", dict(use_mz=False)), ("mz", {})):
        fit_components(rel, tr, **kw)  # compile / warm
        with Timer() as t:
            C = fit_components(rel, tr, **kw)
        timing[key] = t.s
        comps[key] = C
        ci = C.ions.merge(C.table[["component", "rt"]], on="component")
        ci = ci[ci.share >= 0.02]
        P[key] = _resolved(P, main_mz, mode, sigma, ci.rt.to_numpy(), ci.mz.to_numpy(), ci.component.to_numpy())
        Q = _purity(C, tr, sr.compound, sr.merged, pair_set, kind)
        fp.append(dict(method=METHODS[key], compounds=len(Q), eluting_clean=int(Q.clean.sum()),
                       false_split=float(Q.split[Q.clean].mean()), false_merged=float(Q.merged[Q.clean].mean()),
                       merged_with_2_sources=float(Q.merged[Q.n_true >= 2].mean()), background=int(Q.background.sum()),
                       background_flagged_plateau=float(Q.plateau[Q.background].mean()),
                       background_split=float(Q.split[Q.background].mean()),
                       eluting_flagged_plateau=float(Q.plateau[~Q.background].mean()),
                       split_with_2_sources=float(Q.split[Q.n_true >= 2].mean()), seconds=timing[key]))
    # m/z the deconvolution assigns to the two sources of a coalesced ion, against the true difference
    C = comps["mz"]
    dm = []
    for s, c, d_true in zip(P.s[P.dt_sd <= 2.5], P.c[P.dt_sd <= 2.5], P.dm_ppm[P.dt_sd <= 2.5]):
        m = main_mz[s]
        io = C.ions[(np.abs(C.ions.mz - m) / m * 1e6 < 10) & C.ions.coalesced]
        for comp, gi in io.groupby("compound"):
            if set(gi.k) == {0, 1}:
                a = gi[gi.k == 0].mz_source.iloc[0]
                b = gi[gi.k == 1].mz_source.iloc[0]
                rt = C.table.set_index(["compound", "k"]).rt
                # orient by apex: the source eluting first is the one closer to the earlier true apex
                first_true = s if mode[s] < mode[c] else c
                sign = 1.0 if (rt[(comp, 0)] < rt[(comp, 1)]) else -1.0
                est = (a - b) / m * 1e6 * sign  # m/z of the earlier source minus the later
                true = (main_mz[first_true] - main_mz[c if first_true == s else s]) / m * 1e6
                dm.append(dict(true_ppm=true, est_ppm=est, dt_sd=P.dt_sd[(P.s == s) & (P.c == c)].iloc[0]))
                break
    return dict(pairs=P, fp=pd.DataFrame(fp), dm=pd.DataFrame(dm), trace_s=t_tr.s)


def _component_truth(C, tr, label):
    run = tr.run
    sc = run.scan_of_point
    out = np.full(len(C.table), -1)
    T = C.table.set_index("component")
    for comp, g in C.ions.groupby("component"):
        lo, hi = T.loc[comp, "rt_start"], T.loc[comp, "rt_end"]
        ls, ws = [], []
        for t, sh in zip(g.trace, g.share):
            p = tr.points(int(t))
            r = run.rt[sc[p]]
            k = (r >= lo) & (r <= hi)
            ls.append(label[p[k]])
            ws.append(run.intensity[p[k]] * sh)
        s = pd.Series(np.concatenate(ws)).groupby(np.concatenate(ls)).sum()
        s = s[s.index >= 0]
        out[comp] = int(s.idxmax()) if len(s) else -1
    return out


def _relation_truth(rel, tr, label):
    sc = tr.run.scan_of_point
    out = np.full(len(rel), -1)
    for k, (a, s0, s1) in enumerate(zip(rel.anchor, rel.s0, rel.s1)):
        p = tr.points(int(a))
        c = label[p[(sc[p] >= s0) & (sc[p] <= s1)]]
        c = c[c >= 0]
        out[k] = np.bincount(c).argmax() if c.size else -1
    return out


def _centred_error(M, rows, truth, H, min_n=3):
    """Per true compound (best of the entries mapped to it): RMS of centred log abundance minus centred log truth."""
    best = {}
    for r, c in zip(rows, truth):
        x = M[r]
        f = np.isfinite(x) & (x > 0)
        if f.sum() < min_n:
            continue
        lx, lt = np.log(x[f]), np.log(H[c][f])
        e = float(np.sqrt(np.mean(((lx - lx.mean()) - (lt - lt.mean())) ** 2)))
        if c not in best or e < best[c]:
            best[c] = e
    return best


def cross_run():
    sc = metabolomics_scenario(**BATCH["sc"])
    b = make_batch(sc, **{k: v for k, v in BATCH.items() if k != "sc"})
    J = len(b.runs)
    trs, rels, noise, comps, feats = [], [], [], [], []
    for sr in b.runs:
        tr = detect_traces(sr.run)
        g = scan_graph(tr)
        rel = relation_table(g, tr)
        trs.append(tr)
        rels.append(rel)
        noise.append(relation_noise(rel, tr))
        feats.append(detect_peaks(tr, trace_baselines(tr)))
    with Timer() as t_fit:
        comps = [fit_components(r, t) for r, t in zip(rels, trs)]
    A, _, _ = relation_anchors(rels, noise=noise)
    model, _ = cv_drift(A.rt, A.se, K_max=3, n_basis=10)
    grp = relation_groups(rels, model)
    with Timer() as t_link:
        cr = link_components(comps, grp, model)
    apex, sd, half = predict_apex(cr, model, J, comps)
    with Timer() as t_q:
        Q = [quantify(cr, j, trs[j].run, trs[j].precision, apex[:, j], sd[:, j], half[:, j], cv=comps[j].cv)
             for j in range(J)]
    kind = b.compounds.kind.to_numpy()
    H = b.height
    # truth of every cross-run component: majority over its members
    mt = [_component_truth(c, trs[j], b.runs[j].compound) for j, c in enumerate(comps)]
    mem = cr.members.assign(true=[mt[j][c] for j, c in zip(cr.members.run, cr.members.component)])
    xt = mem.groupby("xc")["true"].agg(lambda v: v.value_counts().index[0])
    purity = float(mem.assign(maj=mem.xc.map(xt)).pipe(lambda x: (x.true == x.maj).groupby(x.xc).mean()).mean())
    el = xt[(xt >= 0)]
    el = el[kind[el.to_numpy()] != "background"]
    found = np.zeros((len(cr.table), J), bool)
    for j in range(J):
        found[cr.members[cr.members.run == j].xc.to_numpy(), j] = True
    area = np.stack([q.area for q in Q], 1)
    area_int = np.stack([q.area_int for q in Q], 1)
    area_found = np.where(found, area, np.nan)
    # relation groups and feature groups with truth
    rt_ = [_relation_truth(r, trs[j], b.runs[j].compound) for j, r in enumerate(rels)]
    gm = grp.members.assign(tc=[rt_[j][i] for j, i in zip(grp.members.run, grp.members.rel)])
    gt = gm.groupby("group").tc.agg(lambda v: v.value_counts().index[0])
    gt = gt[gt >= 0]
    gt = gt[kind[gt.to_numpy()] != "background"]
    al = align(feats)
    ft = [feature_truth(f, trs[j], b.runs[j].ion, b.runs[j].ions, b.runs[j].compounds)["compound"].to_numpy()
          for j, f in enumerate(feats)]
    FI = al.feature_index
    fc = np.full(FI.shape, -1)
    for j in range(J):
        m = FI[:, j] >= 0
        fc[m, j] = ft[j][FI[m, j]]
    fmaj = np.array([np.bincount(r[r >= 0]).argmax() if (r >= 0).any() else -1 for r in fc])
    fsel = np.flatnonzero((fmaj >= 0) & (kind[np.maximum(fmaj, 0)] != "background"))
    pair_ids = set(b.compounds[b.compounds.mass_offset_ppm != 0].compound_id)
    pair_ids |= set(b.compounds[b.compounds.formula.isin(b.compounds[b.compounds.mass_offset_ppm != 0].formula)]
                    .compound_id)
    sets = {
        "components, re-queried in every run": (area, el.index.to_numpy(), el.to_numpy()),
        "components, integrated": (area_int, el.index.to_numpy(), el.to_numpy()),
        "components, runs where the fit found them": (area_found, el.index.to_numpy(), el.to_numpy()),
        "relation groups, window height": (grp.height, gt.index.to_numpy(), gt.to_numpy()),
        "relation groups, window area": (grp.area, gt.index.to_numpy(), gt.to_numpy()),
        "feature groups (picked peaks)": (al.intensity, fsel, fmaj[fsel]),
    }
    errs = {k: _centred_error(M, r, t, H) for k, (M, r, t) in sets.items()}
    common = set.intersection(*[set(e) for e in errs.values()])
    rows = []
    for k, e in errs.items():
        for lab, cs in (("all compounds", common), ("with a near-isobaric twin", common & pair_ids),
                        ("without", common - pair_ids)):
            v = np.array([e[c] for c in cs])
            rows.append(dict(quantity=k, compounds=lab, n=len(v), median=float(np.median(v)),
                             q90=float(np.quantile(v, 0.9)), rms=float(np.sqrt(np.mean(v ** 2)))))
    nr = cr.table.set_index("xc").loc[el.index, "n_runs"]
    link = dict(cross_run_components=len(cr.table), eluting=len(el), true_compounds=int(el.nunique()),
                fragments_per_compound=len(el) / max(el.nunique(), 1), member_purity=purity,
                in_all_runs=float((nr == J).mean()), found_cells=float(found[el.index].mean()),
                quantified_cells=float((area[el.index] > 0).mean()), fit_s=t_fit.s, link_s=t_link.s,
                quantify_s=t_q.s, runs=J)
    return dict(errors=pd.DataFrame(rows), link=link)


# ------------------------------------------------------------------------------------------------ real data
def real_stage(dataset, f):
    """Relations and per-run components of one real run, from the cross-run benchmark's tracking + scan-graph
    cache (``bench_xrun.process_real`` must have run first)."""
    path = CACHE / f"comp-real-{dataset}-{f.stem}.pkl"
    if path.exists():
        return pickle.loads(path.read_bytes())
    c = pickle.loads((CACHE / f"xrun-real-{dataset}-{f.stem}-stage1.pkl").read_bytes())
    run = MSRun.load(f)
    tr = Traces(run, c["label"].astype(np.int64), c["order"].astype(np.int64), c["offsets"], c["table"],
                c["precision"], c["params"])
    g = SimpleNamespace(segments=c["seg"], blocks=pd.DataFrame(dict(compound=c["compound"])), prior=c["prior"])
    with Timer() as t_r:
        rel = relation_table(g, tr)
        rn = relation_noise(rel, tr)
    with Timer() as t_c:
        comp = fit_components(rel, tr)
    out = dict(rel=rel, noise=rn, comp=comp, feats=c["feats"], precision=tr.precision, name=f.stem,
               timings=dict(relations=t_r.s, components=t_c.s, n_points=run.n_points))
    path.write_bytes(pickle.dumps(out))
    print(dataset, f.stem, {k: round(v, 1) for k, v in out["timings"].items()}, len(comp.table), "components",
          flush=True)
    return out


def real_link(dataset, files, **link):
    """Relation groups, linking and apex prediction for a real dataset (cached per linking option)."""
    tag = "".join(f"-{k}{v}" for k, v in sorted(link.items()))
    path = CACHE / f"comp-real-{dataset}-link{tag}.pkl"
    if path.exists():
        return pickle.loads(path.read_bytes())
    gpath = CACHE / f"comp-real-{dataset}-groups.pkl"
    st = []
    for f in files:
        s = real_stage(dataset, f)
        if dataset != "metabolomics":
            s["feats"] = None  # the design test does not use picked peaks; M7 needs the memory
        st.append(s)
    comps = [s["comp"] for s in st]
    J = len(st)
    if gpath.exists():
        D = pickle.loads(gpath.read_bytes())
    else:
        rels = [s["rel"] for s in st]
        noise = [s["noise"] for s in st]
        with Timer() as t_d:
            win = drift_window(rels, noise)
            A, _, _ = relation_anchors(rels, noise=noise, rt_window=win)
            model, _ = cv_drift(A.rt, A.se, n_basis=12)
            grp = relation_groups(rels, model)
        del rels, noise, A
        D = dict(model=model, grp=grp, seconds=t_d.s)
        gpath.write_bytes(pickle.dumps(D))
    for s in st:
        s.pop("rel", None)
        s.pop("noise", None)
    model, grp = D["model"], D["grp"]
    with Timer() as t_l:
        cr = link_components(comps, grp, model, **link)
    apex, sd, half = predict_apex(cr, model, J, comps)
    found = np.zeros((len(cr.table), J), bool)
    for j in range(J):
        found[cr.members[cr.members.run == j].xc.to_numpy(), j] = True
    out = dict(names=[s["name"] for s in st], cr=cr, found=found, apex=apex, sd=sd, half=half,
               precision=[s["precision"] for s in st], cv=[c.cv for c in comps],
               fits=[s["comp"].fits.assign(run=j) for j, s in enumerate(st)],
               tables=[s["comp"].table.assign(run=j) for j, s in enumerate(st)],
               ions_coalesced=[s["comp"].ions[s["comp"].ions.coalesced].assign(run=j) for j, s in enumerate(st)],
               grp_height=grp.height, grp_snr=grp.snr, groups=len(grp.table), feats=[s["feats"] for s in st],
               timings=dict(components=[s["timings"]["components"] for s in st],
                            relations=[s["timings"]["relations"] for s in st], drift=D["seconds"], link=t_l.s))
    path.write_bytes(pickle.dumps(out))
    return out


def real_dataset(dataset, link=None, **quant):
    """Linking (cached) and targeted quantification of every cross-run component in every run of a real dataset."""
    files = real_files(dataset)
    if len(files) < 3:
        return None
    link = link or {}
    tag = "".join(f"-{k}{v}" for k, v in sorted({**link, **quant}.items()))
    path = CACHE / f"comp-real-{dataset}-cross{tag}.pkl"
    if path.exists():
        return pickle.loads(path.read_bytes())
    out = real_link(dataset, files, **link)
    cr, apex, sd, half = out["cr"], out["apex"], out["sd"], out["half"]
    Q = []
    with Timer() as t_q:
        for j, f in enumerate(files):
            run = MSRun.load(f)
            parts = []
            for a in range(0, len(cr.table), 50_000):  # chunks bound the chromatogram arrays
                b = a + 50_000
                parts.append(quantify(_subset(cr, a, b), j, run, out["precision"][j], apex[a:b, j], sd[a:b, j],
                                      half[a:b, j], cv=out["cv"][j], **quant))
            Q.append(pd.concat(parts, ignore_index=True))
            del run
    out = dict(out, area=np.stack([q.area for q in Q], 1), area_int=np.stack([q.area_int for q in Q], 1),
               llr=np.stack([q.llr for q in Q], 1), snr=np.stack([q.snr for q in Q], 1),
               rt=np.stack([q.rt for q in Q], 1))
    out["timings"] = dict(out["timings"], quantify=t_q.s)
    path.write_bytes(pickle.dumps(out))
    return out


def _subset(cr, a, b):
    t = cr.table.iloc[a:b].copy()
    t["xc"] = t["xc"] - a
    ions = cr.ions[(cr.ions.xc >= a) & (cr.ions.xc < b)].assign(xc=lambda d: d.xc - a).reset_index(drop=True)
    return type(cr)(t.reset_index(drop=True), cr.members, ions)


PRESENT = 25.0  # likelihood-ratio evidence (chi2 gain, ~5 sd) for using a run's value in ratios


def replicate_spread(R, idx):
    """Spread of log abundance over replicate injections (sd of natural log), per SNR bin, for components,
    relation groups and feature groups."""
    rows = []
    snr_c = np.nanmedian(np.where(R["found"], R["snr"], np.nan), 1)
    sets = {"components, present (LLR ≥ 25)": (np.where(R["llr"] >= PRESENT, R["area"], np.nan), snr_c),
            "components, every run": (R["area"], snr_c),
            "relation groups (height)": (R["grp_height"], np.nanmedian(R["grp_snr"], 1))}
    al = align(R["feats"])
    FI = al.feature_index
    S = np.full(FI.shape, np.nan)
    for j, tb in enumerate(R["feats"]):
        m = FI[:, j] >= 0
        S[m, j] = tb["snr"].to_numpy()[FI[m, j]]
    sets["feature groups (picked peaks)"] = (al.intensity, np.nanmedian(S, 1))
    for lab, (M, snr) in sets.items():
        X = M[:, idx]
        ok = np.isfinite(X).all(1) & (X > 0).all(1)
        sd = np.log(X[ok]).std(1, ddof=1)
        s = snr[ok]
        for lo, hi, bn in ((0, np.inf, "all"), (3, 10, "3-10"), (10, 30, "10-30"), (30, np.inf, "30+")):
            m = (s >= lo) & (s < hi)
            rows.append(dict(groups=lab, snr=bn, n=int(m.sum()), median_sd=float(np.median(sd[m])) if m.any() else np.nan,
                             q90_sd=float(np.quantile(sd[m], 0.9)) if m.any() else np.nan))
    return pd.DataFrame(rows)


def design(R):
    cond = np.array([n.split("_")[0] for n in R["names"]])
    rows, ratios = [], {}
    snr_c = np.nanmedian(np.where(R["found"], R["snr"], np.nan), 1)
    for lab, M, snr in (("components, present (LLR ≥ 25)", np.where(R["llr"] >= PRESENT, R["area"], np.nan), snr_c),
                        ("components, every run", R["area"], snr_c),
                        ("relation groups (height)", R["grp_height"], np.nanmedian(R["grp_snr"], 1))):
        x, ok = log2_ratio(M, cond)
        s = snr[ok]
        ratios[lab] = (x, s)
        for lo, hi, bn in ((0, np.inf, "all"), (3, 10, "3-10"), (10, 30, "10-30"), (30, np.inf, "30+")):
            m = (s >= lo) & (s < hi)
            rows.append(dict(groups=lab, snr=bn, **_mix(x[m])))
    return pd.DataFrame(rows), ratios


def split_stats(R):
    F = pd.concat(R["fits"], ignore_index=True)
    T = pd.concat(R["tables"], ignore_index=True)
    io = pd.concat(R["ions_coalesced"], ignore_index=True)
    piv = io.pivot_table(index=["run", "compound", "trace"], columns="k", values="mz_source", aggfunc="first").dropna()
    dm = ((piv[0] - piv[1]) / piv[1] * 1e6).to_numpy() if len(piv) else np.zeros(0)
    el = ~F.plateau
    return dict(compounds=len(F), plateau=float(F.plateau.mean()), split=float(F.split[el].mean()),
                split_mz_driven=float((F.split & (F.gain_mz > F.gain_int))[el].mean()),
                components=len(T), coalesced_ions=int(len(piv)), dm_ppm=dm,
                median_abs_dm_ppm=float(np.median(np.abs(dm))) if dm.size else np.nan)


# ------------------------------------------------------------------------------------------------ figures
def resolution_figure(P):
    fig = go.Figure()
    bins = [0, 0.5, 1.0, 1.5, 2.5]
    P = P[P.dt_sd <= 2.5].assign(bin=pd.cut(P.dt_sd, bins))
    x = [f"{iv.left:g}-{iv.right:g}" for iv in P.bin.cat.categories]
    for key, nm in METHODS.items():
        y = P.groupby("bin", observed=False)[key].mean().to_numpy()
        fig.add_trace(go.Bar(x=x, y=y, name=nm, marker=dict(color=method_color(["dtw", "global", "macha"][list(METHODS).index(key)]))))
    fig.update_xaxes(title="apex separation of the pair (peak sd)")
    fig.update_yaxes(title="pairs resolved", range=[0, 1], tickformat=".0%")
    return style(fig, 360, title="Near-isobaric pairs (0.5-6 ppm, coalesced) resolved")


def dm_figure(dm):
    fig = go.Figure(go.Scatter(x=dm.true_ppm, y=dm.est_ppm, mode="markers", name="split coalesced ion",
                               marker=dict(color=method_color("macha"), size=7, opacity=0.8)))
    lim = 10.0  # a handful of estimates beyond are drawn at the edge
    dm = dm.assign(est_ppm=dm.est_ppm.clip(-lim, lim))
    fig.data[0].update(x=dm.true_ppm, y=dm.est_ppm)
    fig.add_trace(go.Scatter(x=[-lim, lim], y=[-lim, lim], mode="lines", line=dict(color="#888", dash="dot"),
                             name="estimate = truth"))
    fig.update_xaxes(title="true m/z difference of the sources (ppm, earlier minus later)", range=[-lim, lim])
    fig.update_yaxes(title="estimated from the m/z trajectory (ppm)", range=[-lim, lim])
    return style(fig, 380, title="Source m/z of a coalesced centroid series, recovered by the split")


def ratio_figure(ratios, D, min_snr=10.0):
    fig = go.Figure()
    colors = [method_color("macha"), method_color("greedy"), method_color("global")]
    edges = np.arange(-5, 4.0001, 0.05)
    for (lab, (x, snr)), c in zip(ratios.items(), colors):
        x = x[np.isfinite(x) & (snr >= min_snr)]
        h, _ = np.histogram(x, edges, density=True)
        fig.add_trace(go.Bar(x=0.5 * (edges[1:] + edges[:-1]), y=h, width=0.05, name=f"{lab} ({x.size:,})",
                             opacity=0.55, marker=dict(color=c, line=dict(width=0))))
    shift = float(D[(D.snr == "all") & (D.groups == "components, present (LLR ≥ 25)")]["shift"].iloc[0])
    for nm, o in zip(SPECIES["names"], SPECIES["offsets"]):
        fig.add_vline(x=shift + o, line=dict(color="#888", dash="dot"), annotation_text=nm, annotation_position="top")
    fig.update_layout(barmode="overlay")
    fig.update_xaxes(title="log2 (mean A / mean B)", range=[-4.5, 3])
    fig.update_yaxes(title="density")
    return style(fig, 380, title=f"ProteoBench M7: A/B ratios, median SNR ≥ {min_snr:g}")


def real_dm_figure(dm, label):
    edges = np.arange(-12, 12.001, 0.25)
    h, _ = np.histogram(dm, edges)
    fig = go.Figure(go.Bar(x=0.5 * (edges[1:] + edges[:-1]), y=h, width=0.25, name="split coalesced ions",
                           marker=dict(color=method_color("macha"), line=dict(width=0))))
    fig.update_xaxes(title="m/z of component 0's source minus component 1's (ppm)")
    fig.update_yaxes(title="ions")
    return style(fig, 340, legend=False, title=f"{label}: source m/z difference assigned to coalesced ions")


# ------------------------------------------------------------------------------------------------ page
def compute():
    res = dict(deconv=deconvolution())
    print("deconvolution done", flush=True)
    res["cross"] = cross_run()
    print("cross-run done", flush=True)
    real = {}
    for ds in ("metabolomics", "proteomics"):
        R = real_dataset(ds)
        if R is None:
            continue
        v = dict(names=R["names"], split=split_stats(R), timings=R["timings"], groups=R["groups"],
                 components=len(R["cr"].table), cv=R["cv"],
                 n_runs=R["cr"].table.n_runs.value_counts().sort_index().to_dict(),
                 found=float(R["found"].mean()), quantified=float((R["area"] > 0).mean()))
        if ds == "metabolomics":
            idx = [k for k, n in enumerate(R["names"]) if n.startswith("Sa04")]
            v["replicates"] = replicate_spread(R, idx) if len(idx) >= 3 else None
        else:
            v["design"] = design(R)
        real[ds] = v
        print(ds, "done", flush=True)
    res["real"] = real
    (CACHE / "components-results.pkl").write_bytes(pickle.dumps(res))
    return res


def main(rebuild_page_only: bool = False):
    path = CACHE / "components-results.pkl"
    res = pickle.loads(path.read_bytes()) if rebuild_page_only and path.exists() else compute()
    dc, xr, real = res["deconv"], res["cross"], res["real"]
    P = dc["pairs"]
    co = P[P.dt_sd <= 2.5]
    save_json("components", dict(
        deconvolution=dict(resolved={k: float(co[k].mean()) for k in METHODS}, n_pairs=len(co),
                           false_split=dc["fp"].to_dict("records"),
                           dm_abs_error_ppm=float(np.median(np.abs(dc["dm"].est_ppm - dc["dm"].true_ppm)))
                           if len(dc["dm"]) else None),
        cross_run=dict(link=xr["link"], errors=xr["errors"].to_dict("records")),
        real={ds: dict(split={k: v for k, v in r["split"].items() if k != "dm_ppm"}, timings=r["timings"],
                       components=r["components"], n_runs=r["n_runs"], found=r["found"], quantified=r["quantified"],
                       replicates=None if r.get("replicates") is None else r["replicates"].to_dict("records"),
                       design=None if r.get("design") is None else r["design"][0].to_dict("records"))
              for ds, r in real.items()}))

    page = Page("Components: peaks of relation clusters",
                "A component is one compound's elution in one run: the ions its confident scan-graph relations "
                "join, fitted with one shared exponentially modified Gaussian. Components split where two sources "
                "share an ion, are linked across runs through relation groups, and are quantified in every run "
                "with their learned shape.", NAV)
    E = xr["errors"]
    e_all = E[E.compounds == "all compounds"].set_index("quantity")
    tiles = [("merged pairs resolved", f"{co['mz'].mean():.0%}",
              f"picked peaks {co['peaks'].mean():.0%}; intensity only {co['int'].mean():.0%} (synthetic, ≤ 2.5 sd apart)"),
             ("false splits", f"{dc['fp'].set_index('method').loc[METHODS['mz'], 'false_split']:.1%}",
              "eluting compounds with one true source (synthetic)"),
             ("abundance error, median", f"{e_all.loc['components, re-queried in every run', 'median']:.3f}",
              f"natural log, centred over 6 runs; relation-group height "
              f"{e_all.loc['relation groups, window height', 'median']:.3f}, picked peaks "
              f"{e_all.loc['feature groups (picked peaks)', 'median']:.3f}")]
    if "proteomics" in real and real["proteomics"].get("design") is not None:
        D = real["proteomics"]["design"][0]
        d = D[(D.groups == "components, present (LLR ≥ 25)") & (D.snr == "30+")].iloc[0]
        tiles.append(("M7 yeast at SNR ≥ 30", f"{d['free_mean_yeast']:+.2f}",
                      f"log2 A/B (design +1), sd {d['sd_yeast']:.2f}; human sd {d['sd_human']:.2f}"))
    page.tiles(tiles)
    page.section("What a component is", """
<p><b>Shared shape.</b> The ions joined by a compound's confident relations are fitted together:
y<sub>i</sub>(t) = a<sub>i</sub>·f(t; μ, σ, τ) + b<sub>i</sub>, with f an exponentially modified Gaussian shared by all
ions, one amplitude per ion and a constant baseline per ion (real traces sit on a chemical-background level).
Amplitudes and baselines are solved in closed form by non-negative least squares (variable projection); the shape
is fitted by Levenberg-Marquardt. Every partner trace is a separate noisy view of the same elution, so all of them pin the shape.
The intensity variance is noise<sub>i</sub>² + (cv·model)². The noise level is measured in each ion's own m/z
window, and cv is learned from high-intensity residuals. A scan in the window where an ion has no centroid is
censored at that noise level: it costs only where the model rises above it. This is the dropout floor.</p>
<p><b>Splits.</b> A second component (own apex, shared width, its own amplitude per ion) is kept when it improves
the fit by more than a BIC-style penalty. Most ions then belong to one side: a partner trace that peaks with one source
becomes a template for it. A component whose one-component fit runs into the width bounds or out of its window has
no bounded peak (chemical background, or a weak cluster at the noise level): it is flagged as a plateau and never
split. Two sources closer than 0.75 sd are detected but not separable (their amplitude ratio and source m/z trade
off), so they stay one component, flagged <i>merged</i>.</p>
<p><b>Merged traces.</b> Ions closer than about one resolution FWHM (several ppm) are not resolved, and the
instrument reports their intensity-weighted mean m/z. With two sources in one ion, the centroid m/z follows
m(t) = m<sub>2</sub> + (m<sub>1</sub> − m<sub>2</sub>)·g(t) + β·log I(t), where g(t) is source 1's share of the
signal at time t. It moves from one source to the other as dominance changes. An intensity-dependent shift of a
single ion (the β term) gives the same m/z at the same intensity on both flanks; two sources trace a loop.
Centroid precision (~0.3 ppm at high intensity) is far below the resolution limit, so these m/z residuals enter the
same objective, weighted by the learned precision model. The m/z difference of the two sources has a Gaussian prior
scaled by the coalescence limit, estimated from the run's own closest within-scan centroid spacing.</p>
<p><b>Across runs.</b> Each relation group (one relation matched across runs in reference time) votes for the
components its members belong to. Mutual best-voted components of different runs are linked, keeping one per run.
Pieces found in disjoint sets of runs with the same ions and matching reference time are joined. Relation groups
break where a compound's drift class departs from the mean warp, and the pieces then cover complementary runs.
Each cross-run component is then quantified in <i>every</i> run with its median shape and ion ratios:
within ±0.5 sd of its fitted apex where the run's own fit found it, elsewhere at the apex predicted from the runs
that did (drift class inferred from them), searched under a Gaussian prior. Quantification uses the two ions of
the component's strongest relation group. A relation cluster in dense data can hold several co-eluting analytes,
and a confident relation's two ions certainly belong together. Each ion's chromatogram takes the closest centroid
per scan within the tracker's intensity-dependent tolerance, over a refitted baseline, censored at the ion's own
noise.</p>
<p><b>Presence.</b> Every value comes with a likelihood ratio: χ² of the baselines alone minus χ² with the
component. Where the component is absent, its amplitude measures noise and interference, which compresses ratios
between runs towards 1. Ratios therefore use the runs with LLR ≥ 25 (about 5 sd); "every run" shows what happens
without that.</p>""")
    page.section("Deconvolution on ground truth", f"""
<p>One synthetic run with {len(P)} near-isobaric pairs placed 0.5-6 ppm apart: closer than the resolution, so
co-present ions coalesce into one centroid. They elute within 2.5 peak sd of each other, with heights up to 10×
apart. A pair counts as resolved when two different components (or picked peaks) on its main ion each have an apex
within 0.5 peak sd of one of the true apexes.</p>""")
    page.figure(resolution_figure(P), f"{len(co)} pairs within 2.5 sd. Below ~0.5 sd no method separates them; the m/z "
                "trajectory doubles the resolved share at 0.5-1 sd.")
    fp = dc["fp"].copy()
    page.table(fp[["method", "eluting_clean", "false_split", "false_merged", "split_with_2_sources",
                   "merged_with_2_sources", "background",
                   "background_flagged_plateau", "background_split", "eluting_flagged_plateau", "seconds"]],
               "False splits and plateaus. *_with_2_sources: compounds whose window really holds two true sources "
               "(≥ 10 % each); merged: two sources detected but closer than 0.75 sd, kept as one component.",
               better=dict(false_split="min", split_with_2_sources="max"))
    if len(dc["dm"]):
        page.figure(dm_figure(dc["dm"]), f"{len(dc['dm'])} split coalesced ions. Median absolute error "
                    f"{np.median(np.abs(dc['dm'].est_ppm - dc['dm'].true_ppm)):.2f} ppm.")
    lk = xr["link"]
    page.section("Across runs, on ground truth", f"""
<p>A batch of {lk['runs']} runs with two drift classes, 40 near-isobaric pairs and per-run abundance variation
(log-normal, sd 0.3). {lk['cross_run_components']} cross-run components, {lk['eluting']} of them eluting compounds
({lk['true_compounds']} distinct: {lk['fragments_per_compound']:.2f} per compound). Member purity
{lk['member_purity']:.3f}. {lk['in_all_runs']:.0%} were found by the per-run fit in all runs. The per-run fits
found them in {lk['found_cells']:.0%} of (component, run) cells; targeted re-query quantifies
{lk['quantified_cells']:.0%}.</p>
<p>Error: per true compound, the RMS over runs of centred log abundance minus centred log true height. The best
entry per compound is kept, on the compounds every method quantified in at least 3 runs.</p>""")
    page.table(E, "Abundance error against truth (natural log, centred over runs).", better=dict(median="min", q90="min",
                                                                                                  rms="min"),
               group_col="compounds")
    for ds, r in real.items():
        s = r["split"]
        tm = r["timings"]
        page.section(REAL_LABEL[ds], f"""
<p>{len(r['names'])} runs. Per run: {s['compounds'] // len(r['names']):,} relation clusters fitted
({np.mean(tm['components']):.0f} s per run); {s['plateau']:.0%} plateaus; {s['split']:.1%} of eluting clusters
split, {s['split_mz_driven']:.1%} where the m/z trajectory gave more evidence than intensity. {r['components']:,}
cross-run components ({r['groups']:,} relation groups). The per-run fits found them in {r['found']:.0%} of cells, and
re-query quantified {r['quantified']:.0%}. Learned relative intensity error cv {np.mean(r['cv']):.3f}.
Coalesced ions given two sources: {s['coalesced_ions']:,}, median |Δm| {s['median_abs_dm_ppm']:.1f} ppm.</p>""")
        if s["dm_ppm"].size:
            page.figure(real_dm_figure(s["dm_ppm"], REAL_LABEL[ds]), "If splits of coalesced ions were spurious, the m/z assigned to the "
                        "two sources would differ by ~0. A real near-isobaric pair differs by up to the resolution limit.")
        if r.get("replicates") is not None:
            page.table(r["replicates"], "Sa04 injected three times: sd of natural-log abundance over the triplicate "
                       "(no truth needed; lower is better). Each method's own groups and SNR.",
                       better=dict(median_sd="min", q90_sd="min"), group_col="snr")
        if r.get("design") is not None:
            D, ratios = r["design"]
            page.figure(ratio_figure(ratios, D), "Dotted lines: species positions at the fitted common shift.")
            cols = ["groups", "snr", "quantified", "background", "shift", "sd_human", "sd_yeast", "sd_E. coli",
                    "free_mean_yeast", "free_mean_E. coli", "weight_E. coli"]
            page.table(D[cols], "Three-mode mixture at the design offsets (0, +1, −2); free_mean: modes fitted freely, "
                       "relative to the shift.", better={"sd_human": "min", "sd_yeast": "min"}, group_col="snr")
    page.write(REPORTS / "components.html")


if __name__ == "__main__":
    main(rebuild_page_only="--page" in sys.argv)
