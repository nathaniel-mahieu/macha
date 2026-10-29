import base64
import gzip
import io
import struct
import zlib

import numpy as np
import pytest

from macha import MSRun, stitch_scans
from macha.io import numpress, read_mzml


def _toy_run():
    rng = np.random.default_rng(0)
    scans, mzs, its = [], [], []
    for s in range(20):
        m = np.sort(rng.uniform(100, 500, 30))
        scans.append(np.full(30, s))
        mzs.append(m)
        its.append(rng.uniform(1e3, 1e5, 30))
    return MSRun.from_points(np.concatenate(scans), np.concatenate(mzs), np.concatenate(its).astype(np.float32),
                             np.arange(20) * 0.5)


def test_csr_invariants_and_accessors():
    run = _toy_run()
    assert run.n_scans == 20 and run.n_points == 600
    assert np.all(np.diff(run.offsets) == 30)
    for s in range(run.n_scans):
        mz, _ = run.scan(s)
        assert np.all(np.diff(mz) >= 0)
    assert np.array_equal(run.scan_of_point, np.repeat(np.arange(20), 30))
    tic = run.tic()
    assert np.allclose(tic[3], run.intensity[90:120].astype(float).sum(), rtol=1e-6)
    assert np.allclose(run.bpc()[3], run.intensity[90:120].max())


def test_xic_matches_bruteforce():
    run = _toy_run()
    target = run.mz[45]
    x = run.xic(target, ppm=20)
    tol = target * 20e-6
    brute = np.array([run.intensity[(run.scan_of_point == s) & (np.abs(run.mz - target) <= tol)].sum()
                      for s in range(run.n_scans)])
    assert np.allclose(x, brute, rtol=1e-5)


def test_slice_and_filter():
    run = _toy_run()
    sl = run.slice(mz_range=(200, 300), rt_range=(2.0, 6.0))
    assert sl.rt.min() >= 2.0 and sl.rt.max() <= 6.0
    assert sl.mz.min() >= 200 and sl.mz.max() <= 300
    assert sl.offsets[-1] == sl.n_points
    f = run.filter_intensity(5e4)
    assert f.intensity.min() >= 5e4


def test_save_load_roundtrip(tmp_path):
    run = _toy_run()
    run.meta["note"] = "x"
    run.save(tmp_path / "r.npz")
    back = MSRun.load(tmp_path / "r.npz")
    assert np.array_equal(back.mz, run.mz) and np.array_equal(back.intensity, run.intensity)
    assert back.meta["note"] == "x"


def test_stitch_scans_merges_segments():
    # alternating low / high mass segments, as in multi-range acquisitions
    scans, mzs = [], []
    for s in range(10):
        lo = s % 2 == 0
        m = np.linspace(100, 390, 5) if lo else np.linspace(410, 900, 5)
        scans.append(np.full(5, s))
        mzs.append(m)
    run = MSRun.from_points(np.concatenate(scans), np.concatenate(mzs), np.ones(50, np.float32), np.arange(10.0))
    st = stitch_scans(run, 2, mz_ranges=[(0, 400), (400, 1000)])
    assert st.n_scans == 5
    assert np.all(np.diff(st.offsets) == 10)
    assert np.allclose(st.rt, [0.5, 2.5, 4.5, 6.5, 8.5])


def _b64(arr, dtype, compress=True):
    raw = np.asarray(arr, dtype=dtype).tobytes()
    return base64.b64encode(zlib.compress(raw) if compress else raw).decode()


def _mzml(spectra):
    out = ['<?xml version="1.0"?><mzML xmlns="http://psi.hupo.org/ms/mzml"><run><spectrumList count="%d">' % len(spectra)]
    for i, sp in enumerate(spectra):
        out.append(f'<spectrum index="{i}" id="scan={i + 1}">')
        out.append(f'<cvParam accession="MS:1000511" value="{sp["level"]}"/>')
        out.append('<cvParam accession="MS:1000130" value=""/><cvParam accession="MS:1000127" value=""/>')
        out.append(f'<scanList><scan><cvParam accession="MS:1000016" value="{sp["rt_min"]}" unitName="minute"/>'
                   '</scan></scanList><binaryDataArrayList count="2">')
        out.append('<binaryDataArray><cvParam accession="MS:1000523"/><cvParam accession="MS:1000574"/>'
                   f'<cvParam accession="MS:1000514"/><binary>{_b64(sp["mz"], np.float64)}</binary></binaryDataArray>')
        out.append('<binaryDataArray><cvParam accession="MS:1000521"/><cvParam accession="MS:1000576"/>'
                   f'<cvParam accession="MS:1000515"/><binary>{_b64(sp["i"], np.float32, False)}</binary>'
                   '</binaryDataArray></binaryDataArrayList></spectrum>')
    out.append("</spectrumList></run></mzML>")
    return "".join(out).encode()


@pytest.mark.parametrize("gz", [False, True])
def test_read_mzml_skips_ms2_and_converts_minutes(tmp_path, gz):
    spectra = [
        dict(level=1, rt_min=0.5, mz=[300.0, 100.0, 200.0], i=[3, 1, 2]),
        dict(level=2, rt_min=0.51, mz=[50.0], i=[9]),
        dict(level=1, rt_min=1.0, mz=[150.0, 250.0], i=[5, 6]),
    ]
    data = _mzml(spectra)
    path = tmp_path / ("x.mzML.gz" if gz else "x.mzML")
    path.write_bytes(gzip.compress(data) if gz else data)
    run = read_mzml(path)
    assert run.n_scans == 2
    assert np.allclose(run.rt, [30.0, 60.0])
    assert np.allclose(run.mz, [100, 200, 300, 150, 250])  # sorted within scan
    assert np.allclose(run.intensity, [1, 2, 3, 5, 6])
    run2 = read_mzml(io.BytesIO(data), rt_range=(0, 45))
    assert run2.n_scans == 1


# ---------------------------------------------------------------- numpress (reference encoders)
def _encode_int(x):
    x &= 0xFFFFFFFF
    mask = 0xF0000000
    init = x & mask
    if init == 0:
        l = 8
        for i in range(8):
            if x & (mask >> (4 * i)):
                l = i
                break
        return [l] + [(x >> (4 * (i - l))) & 0xF for i in range(l, 8)]
    if init == mask:
        l = 7
        for i in range(8):
            m = mask >> (4 * i)
            if (x & m) != m:
                l = i
                break
        return [l + 8] + [(x >> (4 * (i - l))) & 0xF for i in range(l, 8)]
    return [0] + [(x >> (4 * i)) & 0xF for i in range(8)]


def _encode_linear(values, fixed):
    ints = [int(round(v * fixed)) for v in values]
    out = bytearray(struct.pack(">d", fixed))
    out += struct.pack("<i", ints[0]) + struct.pack("<i", ints[1])
    half = []
    for k in range(2, len(ints)):
        diff = ints[k] - (2 * ints[k - 1] - ints[k - 2])
        half += _encode_int(diff)
    for k in range(0, len(half) - 1, 2):
        out.append(((half[k] & 0xF) << 4) | (half[k + 1] & 0xF))
    if len(half) % 2:
        out.append((half[-1] & 0xF) << 4)
    return bytes(out)


def test_numpress_linear_roundtrip():
    rng = np.random.default_rng(3)
    vals = np.sort(rng.uniform(100, 1500, 200))
    fixed = 1e5
    dec = numpress.decode("linear", _encode_linear(vals, fixed))
    assert dec.size == vals.size
    assert np.max(np.abs(dec - vals)) < 1.0 / fixed


def test_numpress_slof_roundtrip():
    vals = np.array([0.0, 10.0, 1234.5, 1e5])
    fixed = 1000.0
    enc = struct.pack(">d", fixed) + b"".join(struct.pack("<H", int(round(np.log(v + 1) * fixed))) for v in vals)
    dec = numpress.decode("slof", enc)
    assert np.allclose(dec, vals, rtol=2e-3, atol=1e-3)


@pytest.mark.parametrize("suffix", [".mzML", ".mzML.gz"])
def test_write_read_roundtrip(tmp_path, suffix):
    from macha.io import write_mzml

    run = _toy_run()
    write_mzml(run, tmp_path / f"r{suffix}")
    back = read_mzml(tmp_path / f"r{suffix}")
    assert back.n_scans == run.n_scans
    assert np.allclose(back.rt, run.rt, atol=1e-5)
    assert np.array_equal(back.mz, run.mz) and np.array_equal(back.intensity, run.intensity)


def test_stitch_uses_scan_windows_and_splits_overlap():
    scans, mzs, wins = [], [], []
    for s in range(6):
        lo = s % 2 == 0
        m = np.array([150.0, 299.5, 299.7]) if lo else np.array([299.2, 299.6, 800.0])
        scans.append(np.full(3, s))
        mzs.append(m)
        wins.append((100.0, 300.0) if lo else (299.0, 1500.0))
    run = MSRun.from_points(np.concatenate(scans), np.concatenate(mzs), np.ones(18, np.float32), np.arange(6.0))
    run.scan_window = np.array(wins)
    st = stitch_scans(run, 2)
    # overlap 299-300 is split at 299.5: low window keeps < 299.5, high window keeps >= 299.5
    assert np.allclose(st.scan(0)[0], [150.0, 299.6, 800.0])
    run.scan_window = np.array(wins[1:] + wins[:1])
    run.scan_window[2] = (100.0, 300.0)  # break the cycle
    with pytest.raises(ValueError):
        stitch_scans(run, 2)
