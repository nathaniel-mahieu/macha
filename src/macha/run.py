"""Centroided MS1 run container.

An :class:`MSRun` stores every centroid of a run in flat arrays using a CSR
layout: the centroids of scan ``s`` live in ``[offsets[s], offsets[s + 1])``
and are sorted by m/z within the scan.  This keeps memory compact (12 bytes
per centroid) and lets every algorithm in the library operate on contiguous
arrays without Python-level loops over scans.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class MSRun:
    rt: np.ndarray
    offsets: np.ndarray
    mz: np.ndarray
    intensity: np.ndarray
    polarity: int = 1
    scan_number: np.ndarray | None = None
    scan_window: np.ndarray | None = None
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.rt = np.ascontiguousarray(self.rt, dtype=np.float64)
        self.offsets = np.ascontiguousarray(self.offsets, dtype=np.int64)
        self.mz = np.ascontiguousarray(self.mz, dtype=np.float64)
        self.intensity = np.ascontiguousarray(self.intensity, dtype=np.float32)
        if self.offsets.shape[0] != self.rt.shape[0] + 1:
            raise ValueError("offsets must have length n_scans + 1")
        if self.offsets[-1] != self.mz.shape[0] or self.mz.shape != self.intensity.shape:
            raise ValueError("offsets[-1] must equal the number of centroids")
        if self.scan_number is None:
            self.scan_number = np.arange(self.n_scans, dtype=np.int64)
        self._scan_of_point = None

    # ------------------------------------------------------------------ basics
    @property
    def n_scans(self) -> int:
        return int(self.rt.shape[0])

    @property
    def n_points(self) -> int:
        return int(self.mz.shape[0])

    @property
    def scan_of_point(self) -> np.ndarray:
        """Scan index (0-based, into ``rt``) of every centroid."""
        if self._scan_of_point is None:
            self._scan_of_point = np.repeat(
                np.arange(self.n_scans, dtype=np.int32), np.diff(self.offsets)
            )
        return self._scan_of_point

    def scan(self, s: int) -> tuple[np.ndarray, np.ndarray]:
        a, b = self.offsets[s], self.offsets[s + 1]
        return self.mz[a:b], self.intensity[a:b]

    @property
    def cycle_time(self) -> float:
        return float(np.median(np.diff(self.rt))) if self.n_scans > 1 else 1.0

    def tic(self) -> np.ndarray:
        return np.add.reduceat(self.intensity.astype(np.float64), self.offsets[:-1]) * (
            np.diff(self.offsets) > 0
        )

    def bpc(self) -> np.ndarray:
        out = np.zeros(self.n_scans)
        nonempty = np.diff(self.offsets) > 0
        out[nonempty] = np.maximum.reduceat(self.intensity, self.offsets[:-1][nonempty])
        return out

    def xic(self, mz: float, ppm: float = 5.0) -> np.ndarray:
        """Summed intensity per scan within +/- ppm of ``mz`` (vectorised)."""
        tol = mz * ppm * 1e-6
        lo = _scanwise_searchsorted(self.mz, self.offsets, mz - tol, "left")
        hi = _scanwise_searchsorted(self.mz, self.offsets, mz + tol, "right")
        cs = np.concatenate([[0.0], np.cumsum(self.intensity, dtype=np.float64)])
        return cs[hi] - cs[lo]

    # ----------------------------------------------------------------- slicing
    def slice(self, mz_range=None, rt_range=None) -> "MSRun":
        s0, s1 = 0, self.n_scans
        if rt_range is not None:
            s0 = int(np.searchsorted(self.rt, rt_range[0], "left"))
            s1 = int(np.searchsorted(self.rt, rt_range[1], "right"))
        a, b = self.offsets[s0], self.offsets[s1]
        mz, it = self.mz[a:b], self.intensity[a:b]
        offs = self.offsets[s0 : s1 + 1] - a
        if mz_range is not None:
            keep = (mz >= mz_range[0]) & (mz <= mz_range[1])
            counts = np.add.reduceat(keep.astype(np.int64), offs[:-1]) * (np.diff(offs) > 0)
            offs = np.concatenate([[0], np.cumsum(counts)])
            mz, it = mz[keep], it[keep]
        sw = None if self.scan_window is None else self.scan_window[s0:s1]
        meta = dict(self.meta)
        meta["slice"] = {"mz_range": mz_range, "rt_range": rt_range}
        return MSRun(
            rt=self.rt[s0:s1], offsets=offs, mz=mz, intensity=it, polarity=self.polarity,
            scan_number=self.scan_number[s0:s1], scan_window=sw, meta=meta,
        )

    def filter_intensity(self, min_intensity: float) -> "MSRun":
        keep = self.intensity >= min_intensity
        counts = np.add.reduceat(keep.astype(np.int64), self.offsets[:-1]) * (np.diff(self.offsets) > 0)
        offs = np.concatenate([[0], np.cumsum(counts)])
        return MSRun(self.rt, offs, self.mz[keep], self.intensity[keep], self.polarity,
                     self.scan_number, self.scan_window, dict(self.meta))

    # --------------------------------------------------------------------- I/O
    def save(self, path: str | Path, compress: bool = True) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        arrays = dict(rt=self.rt, offsets=self.offsets, mz=self.mz, intensity=self.intensity,
                      scan_number=self.scan_number,
                      polarity=np.array([self.polarity]),
                      meta=np.array([json.dumps(self.meta, default=str)]))
        if self.scan_window is not None:
            arrays["scan_window"] = self.scan_window
        (np.savez_compressed if compress else np.savez)(path, **arrays)

    @classmethod
    def load(cls, path: str | Path) -> "MSRun":
        with np.load(path, allow_pickle=False) as z:
            return cls(
                rt=z["rt"], offsets=z["offsets"], mz=z["mz"], intensity=z["intensity"],
                polarity=int(z["polarity"][0]), scan_number=z["scan_number"],
                scan_window=z["scan_window"] if "scan_window" in z.files else None,
                meta=json.loads(str(z["meta"][0])),
            )

    @classmethod
    def from_points(cls, scan: np.ndarray, mz: np.ndarray, intensity: np.ndarray, rt: np.ndarray,
                    polarity: int = 1, meta: dict | None = None) -> "MSRun":
        """Build a run from unordered (scan, mz, intensity) triples."""
        order = np.lexsort((mz, scan))
        scan, mz, intensity = scan[order], mz[order], intensity[order]
        counts = np.bincount(scan, minlength=len(rt))
        offsets = np.concatenate([[0], np.cumsum(counts)])
        return cls(rt=rt, offsets=offsets, mz=mz, intensity=intensity, polarity=polarity,
                   meta=meta or {})

    def __repr__(self) -> str:
        rt = f"{self.rt[0]:.1f}-{self.rt[-1]:.1f}s" if self.n_scans else "empty"
        return f"MSRun({self.n_scans} scans, {self.n_points:,} centroids, rt {rt}, polarity {self.polarity:+d})"


def stitch_scans(run: MSRun, n_segments: int, mz_ranges=None) -> MSRun:
    """Merge every ``n_segments`` consecutive scans into one scan.

    Instruments can acquire a cycle as several MS1 scans that each cover a
    different m/z range (the original ``rawdata(splits=...)``).  Each merged
    scan takes the mean retention time of its segments.  ``mz_ranges``
    optionally clips segment ``i`` to ``mz_ranges[i]`` before merging.  When it
    is omitted and the run carries scan windows, the windows are used.
    """
    n = (run.n_scans // n_segments) * n_segments
    if mz_ranges is None and run.scan_window is not None:
        win = run.scan_window[:n]
        if np.isfinite(win).all():
            # the acquisition must cycle through the same windows in a fixed order
            ref = win[:n_segments]
            if not np.allclose(win.reshape(-1, n_segments, 2), ref[None], atol=1e-3):
                raise ValueError("scan windows do not repeat every n_segments scans; pass mz_ranges explicitly "
                                 "or select the cycle-complete part of the run")
            mz_ranges = [list(r) for r in ref]
            # windows that overlap are cut at the midpoint so shared m/z is not counted twice
            order = np.argsort(ref[:, 0])
            for a, b in zip(order[:-1], order[1:]):
                if mz_ranges[a][1] > mz_ranges[b][0]:
                    mid = 0.5 * (mz_ranges[a][1] + mz_ranges[b][0])
                    mz_ranges[a][1] = mid
                    mz_ranges[b][0] = mid
    scan = run.scan_of_point
    keep = scan < n
    if mz_ranges is not None:
        seg = scan % n_segments
        lo = np.array([r[0] for r in mz_ranges])[seg]
        hi = np.array([r[1] for r in mz_ranges])[seg]
        keep &= (run.mz >= lo) & (run.mz < hi)
    new_scan = (scan[keep] // n_segments).astype(np.int64)
    rt = run.rt[:n].reshape(-1, n_segments).mean(axis=1)
    meta = dict(run.meta, stitched=n_segments)
    return MSRun.from_points(new_scan, run.mz[keep], run.intensity[keep], rt, run.polarity, meta)


def _scanwise_searchsorted(mz, offsets, value, side):
    from ._kernels import scanwise_searchsorted

    return scanwise_searchsorted(mz, offsets, float(value), side == "left")
