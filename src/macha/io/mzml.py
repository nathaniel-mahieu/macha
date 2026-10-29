"""Streaming mzML reader.

Reads MS1 centroids from local files, ``.gz`` files or any binary stream
(including an HTTP response) in a single forward pass.  MS2+ spectra are
skipped without decoding, so DIA/DDA files can be reduced to their MS1
content while they download.  Supports zlib / uncompressed 32- and 64-bit
arrays and MS-Numpress (linear, pic, slof).
"""

from __future__ import annotations

import base64
import gzip
import io
import zlib
from pathlib import Path
from typing import BinaryIO, Callable, Iterator
from xml.etree import ElementTree as ET

import numpy as np

from ..run import MSRun
from . import numpress

_NS = "{http://psi.hupo.org/ms/mzml}"

# PSI-MS accessions used below
MS_LEVEL = "MS:1000511"
SCAN_START = "MS:1000016"
POSITIVE, NEGATIVE = "MS:1000130", "MS:1000129"
CENTROID, PROFILE = "MS:1000127", "MS:1000128"
MZ_ARRAY, INT_ARRAY = "MS:1000514", "MS:1000515"
F32, F64 = "MS:1000521", "MS:1000523"
ZLIB, NOCOMP = "MS:1000574", "MS:1000576"
NUMPRESS = {
    "MS:1002312": ("linear", False), "MS:1002313": ("pic", False), "MS:1002314": ("slof", False),
    "MS:1002746": ("linear", True), "MS:1002747": ("pic", True), "MS:1002748": ("slof", True),
}
WIN_LO, WIN_HI = "MS:1000501", "MS:1000500"
FILTER = "MS:1000512"


def _tag(e) -> str:
    t = e.tag
    return t[len(_NS):] if t.startswith(_NS) else t


def _decode_array(bda) -> tuple[str | None, np.ndarray]:
    kind, dtype, comp, npress = None, np.float64, "none", None
    text = None
    for child in bda:
        t = _tag(child)
        if t == "cvParam":
            acc = child.get("accession")
            if acc == MZ_ARRAY:
                kind = "mz"
            elif acc == INT_ARRAY:
                kind = "intensity"
            elif acc == F32:
                dtype = np.float32
            elif acc == F64:
                dtype = np.float64
            elif acc == ZLIB:
                comp = "zlib"
            elif acc in NUMPRESS:
                npress = NUMPRESS[acc]
        elif t == "binary":
            text = child.text or ""
    if kind is None:
        return None, np.empty(0)
    raw = base64.b64decode(text) if text else b""
    if npress is not None:
        method, zipped = npress
        if zipped or comp == "zlib":
            raw = zlib.decompress(raw)
        return kind, numpress.decode(method, raw)
    if comp == "zlib" and raw:
        raw = zlib.decompress(raw)
    return kind, np.frombuffer(raw, dtype=dtype).astype(np.float64)


def iter_spectra(stream: BinaryIO, ms_levels=(1,), on_progress: Callable | None = None) -> Iterator[dict]:
    """Yield dicts for each spectrum whose MS level is in ``ms_levels``.

    Each dict has ``rt`` (seconds), ``mz``, ``intensity``, ``polarity``,
    ``centroid``, ``index``, ``id`` and optional ``window`` (lo, hi) and
    ``filter`` (Thermo filter string).
    """
    container = None
    n_seen = 0
    for event, elem in ET.iterparse(stream, events=("start", "end")):
        if event == "start":
            if _tag(elem) == "spectrumList":
                container = elem
            continue
        if _tag(elem) != "spectrum":
            continue
        n_seen += 1
        level, rt, pol, cent = None, None, 0, None
        lo = hi = None
        filt = None
        arrays = {}
        for node in elem.iter():
            t = _tag(node)
            if t == "cvParam":
                acc = node.get("accession")
                if acc == MS_LEVEL:
                    level = int(node.get("value"))
                    if level not in ms_levels:
                        break
                elif acc == SCAN_START:
                    rt = float(node.get("value"))
                    if node.get("unitName", "minute") == "minute" or node.get("unitAccession") == "UO:0000031":
                        rt *= 60.0
                elif acc == POSITIVE:
                    pol = 1
                elif acc == NEGATIVE:
                    pol = -1
                elif acc == CENTROID:
                    cent = True
                elif acc == PROFILE:
                    cent = False
                elif acc == WIN_LO:
                    lo = float(node.get("value"))
                elif acc == WIN_HI:
                    hi = float(node.get("value"))
                elif acc == FILTER:
                    filt = node.get("value")
            elif t == "binaryDataArray":
                k, arr = _decode_array(node)
                if k is not None:
                    arrays[k] = arr
        if level in ms_levels:
            mz = arrays.get("mz", np.empty(0))
            it = arrays.get("intensity", np.empty(0))
            yield dict(index=int(elem.get("index", n_seen - 1)), id=elem.get("id"), level=level,
                       rt=rt, polarity=pol, centroid=cent, mz=mz, intensity=it,
                       window=(lo, hi) if lo is not None else None, filter=filt)
        if on_progress is not None and n_seen % 2000 == 0:
            on_progress(n_seen)
        elem.clear()
        if container is not None and len(container):
            container.clear()


def open_stream(source: str | Path | BinaryIO) -> BinaryIO:
    """Open a path (``.gz`` aware) or pass a stream through, sniffing gzip magic."""
    if hasattr(source, "read"):
        stream = source
    else:
        stream = open(source, "rb")
    peekable = io.BufferedReader(stream) if not hasattr(stream, "peek") else stream
    head = peekable.peek(2)[:2]
    if head == b"\x1f\x8b":
        return gzip.GzipFile(fileobj=peekable)
    return peekable


def read_mzml(source, ms_level: int = 1, rt_range=None, mz_range=None, min_intensity: float = 0.0,
              polarity: int | None = None, on_progress=None) -> MSRun:
    """Read the MS1 (default) centroids of an mzML source into an :class:`MSRun`.

    ``rt_range`` (seconds) stops parsing once past the upper bound, which saves
    work when streaming.  Profile spectra are rejected, as the library works
    on centroids.
    """
    stream = open_stream(source)
    rts, mzs, its, nums, wins = [], [], [], [], []
    pols = []
    filters = []
    try:
        for sp in iter_spectra(stream, (ms_level,), on_progress):
            if sp["centroid"] is False:
                raise ValueError("profile-mode spectrum found; centroid the data first (e.g. msconvert peakPicking)")
            if polarity is not None and sp["polarity"] not in (0, polarity):
                continue
            rt = sp["rt"]
            if rt_range is not None:
                if rt < rt_range[0]:
                    continue
                if rt > rt_range[1]:
                    break
            mz, it = sp["mz"], sp["intensity"]
            keep = it > min_intensity
            if mz_range is not None:
                keep &= (mz >= mz_range[0]) & (mz <= mz_range[1])
            mz, it = mz[keep], it[keep]
            order = np.argsort(mz, kind="stable")
            rts.append(rt)
            mzs.append(mz[order])
            its.append(it[order].astype(np.float32))
            nums.append(sp["index"])
            pols.append(sp["polarity"])
            wins.append(sp["window"] or (np.nan, np.nan))
            filters.append(sp["filter"])
    finally:
        if not hasattr(source, "read"):
            stream.close()
    counts = np.array([len(m) for m in mzs], dtype=np.int64)
    pol_values = [p for p in pols if p != 0]
    pol = int(np.sign(np.sum(pol_values))) if pol_values else 1
    meta = {"source": str(source) if not hasattr(source, "read") else "stream",
            "n_distinct_filters": len(set(f for f in filters if f))}
    if filters and filters[0]:
        meta["first_filters"] = list(dict.fromkeys(f for f in filters[:20] if f))[:6]
    return MSRun(
        rt=np.array(rts, dtype=np.float64),
        offsets=np.concatenate([[0], np.cumsum(counts)]),
        mz=np.concatenate(mzs) if mzs else np.empty(0),
        intensity=np.concatenate(its) if its else np.empty(0, np.float32),
        polarity=pol or 1,
        scan_number=np.array(nums, dtype=np.int64),
        scan_window=np.array(wins, dtype=np.float64) if wins else None,
        meta=meta,
    )


def write_mzml(run: MSRun, path, compress: bool = True) -> None:
    """Write an :class:`MSRun` as a minimal MS1 centroid mzML (64-bit m/z, 32-bit intensity)."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open

    def enc(a, dt):
        raw = np.ascontiguousarray(a, dtype=dt).tobytes()
        return base64.b64encode(zlib.compress(raw) if compress else raw).decode()

    comp = ZLIB if compress else NOCOMP
    pol = POSITIVE if run.polarity >= 0 else NEGATIVE
    with opener(path, "wt") as fh:
        fh.write('<?xml version="1.0" encoding="utf-8"?>\n<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">'
                 f'<run id="macha"><spectrumList count="{run.n_scans}">\n')
        for s in range(run.n_scans):
            mz, it = run.scan(s)
            fh.write(f'<spectrum index="{s}" id="scan={int(run.scan_number[s]) + 1}" defaultArrayLength="{mz.size}">'
                     f'<cvParam cvRef="MS" accession="{MS_LEVEL}" value="1"/><cvParam cvRef="MS" accession="{pol}"/>'
                     f'<cvParam cvRef="MS" accession="{CENTROID}"/><scanList count="1"><scan>'
                     f'<cvParam cvRef="MS" accession="{SCAN_START}" value="{run.rt[s]:.6f}" unitName="second" '
                     f'unitAccession="UO:0000010"/></scan></scanList><binaryDataArrayList count="2">'
                     f'<binaryDataArray><cvParam cvRef="MS" accession="{F64}"/><cvParam cvRef="MS" accession="{comp}"/>'
                     f'<cvParam cvRef="MS" accession="{MZ_ARRAY}"/><binary>{enc(mz, np.float64)}</binary></binaryDataArray>'
                     f'<binaryDataArray><cvParam cvRef="MS" accession="{F32}"/><cvParam cvRef="MS" accession="{comp}"/>'
                     f'<cvParam cvRef="MS" accession="{INT_ARRAY}"/><binary>{enc(it, np.float32)}</binary>'
                     '</binaryDataArray></binaryDataArrayList></spectrum>\n')
        fh.write("</spectrumList></run></mzML>\n")
