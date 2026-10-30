"""Data-independent acquisition: MS1 plus one MS2 time series per isolation window.

A DIA cycle visits every isolation window once, so the MS2 spectra of one window form an ordinary LC-MS run of
fragment ions (one "scan" per cycle).  Each window is stored as an :class:`~macha.run.MSRun`, which lets trace
detection, baselines, peaks and the scan-level graph run on fragments unchanged.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from ..run import MSRun
from .mzml import iter_spectra, open_stream


@dataclass
class DIARun:
    ms1: MSRun
    windows: pd.DataFrame  # window, target, lo, hi (isolation bounds, m/z), n_spectra
    ms2: list[MSRun]  # ms2[w]: the MS2 spectra of window w over time
    meta: dict = field(default_factory=dict)

    def window_of(self, mz: float) -> int:
        """Index of the isolation window holding precursor m/z (-1 if none)."""
        w = self.windows
        hit = np.flatnonzero((w["lo"].to_numpy() <= mz) & (mz < w["hi"].to_numpy()))
        return int(hit[0]) if hit.size else -1

    def save(self, path: str | Path) -> None:
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        self.ms1.save(path / "ms1.npz")
        for k, r in enumerate(self.ms2):
            r.save(path / f"w{k:04d}.npz")
        self.windows.to_csv(path / "windows.csv", index=False)
        (path / "meta.json").write_text(json.dumps(self.meta, default=str))

    @classmethod
    def load(cls, path: str | Path, windows=None) -> "DIARun":
        """Load a saved run; ``windows`` (indices) loads only those MS2 windows (others are None)."""
        path = Path(path)
        win = pd.read_csv(path / "windows.csv")
        keep = set(range(len(win)) if windows is None else windows)
        ms2 = [MSRun.load(path / f"w{k:04d}.npz") if k in keep else None for k in range(len(win))]
        meta = json.loads((path / "meta.json").read_text()) if (path / "meta.json").exists() else {}
        return cls(MSRun.load(path / "ms1.npz"), win, ms2, meta)


def _to_run(spectra, polarity, meta) -> MSRun:
    rts = np.array([s[0] for s in spectra], dtype=np.float64)
    order = np.argsort(rts, kind="stable")
    mzs = [spectra[i][1] for i in order]
    its = [spectra[i][2] for i in order]
    counts = np.array([len(m) for m in mzs], dtype=np.int64)
    return MSRun(rt=rts[order], offsets=np.concatenate([[0], np.cumsum(counts)]),
                 mz=np.concatenate(mzs) if mzs else np.empty(0),
                 intensity=np.concatenate(its) if its else np.empty(0, np.float32), polarity=polarity,
                 scan_number=np.array([spectra[i][3] for i in order], dtype=np.int64), meta=meta)


def read_dia(source, rt_range=None, min_intensity: float = 0.0, on_progress=None) -> DIARun:
    """Read a DIA mzML (centroided) into MS1 and per-isolation-window MS2 runs.

    Windows are identified by their isolation bounds; ``rt_range`` (seconds) stops parsing past its end.
    """
    stream = open_stream(source)
    ms1, by_win = [], {}
    pols = []
    try:
        for sp in iter_spectra(stream, (1, 2), on_progress):
            if sp["centroid"] is False:
                raise ValueError("profile-mode spectrum found; centroid the data first (e.g. msconvert peakPicking)")
            rt = sp["rt"]
            if rt_range is not None:
                if rt < rt_range[0]:
                    continue
                if rt > rt_range[1]:
                    break
            mz, it = sp["mz"], sp["intensity"]
            keep = it > min_intensity
            order = np.argsort(mz[keep], kind="stable")
            row = (rt, mz[keep][order], it[keep][order].astype(np.float32), sp["index"])
            pols.append(sp["polarity"])
            if sp["level"] == 1:
                ms1.append(row)
            elif sp.get("isolation") is not None:
                lo, hi = sp["isolation"]
                by_win.setdefault((round(lo, 3), round(hi, 3)), []).append(row)
    finally:
        if not hasattr(source, "read"):
            stream.close()
    pol_values = [p for p in pols if p != 0]
    pol = int(np.sign(np.sum(pol_values))) if pol_values else 1
    keys = sorted(by_win)
    meta = {"source": str(source) if not hasattr(source, "read") else "stream"}
    windows = pd.DataFrame(dict(window=np.arange(len(keys)), lo=[k[0] for k in keys], hi=[k[1] for k in keys],
                                n_spectra=[len(by_win[k]) for k in keys]))
    windows.insert(1, "target", 0.5 * (windows.lo + windows.hi))
    ms2 = [_to_run(by_win[k], pol, dict(meta, isolation=list(k))) for k in keys]
    return DIARun(_to_run(ms1, pol, dict(meta, level=1)), windows, ms2, meta)
