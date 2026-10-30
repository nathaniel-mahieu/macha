"""MS-Numpress decoders (Teleman et al., MCP 2014), ported for numba.

Only decoding is implemented; it is needed to read files written with
msconvert's ``--numpress*`` options.  Follows the reference MSNumpress.cpp.
"""

from __future__ import annotations

import numba as nb
import numpy as np


def decode(method: str, data: bytes) -> np.ndarray:
    buf = np.frombuffer(data, dtype=np.uint8)
    if method == "linear":
        return _decode_linear(buf)
    if method == "pic":
        return _decode_pic(buf)
    if method == "slof":
        return _decode_slof(buf)
    raise ValueError(f"unknown numpress method {method}")


@nb.njit(cache=True)
def _nibble(buf, pos, half):
    b = buf[pos]
    return (b >> 4) & 0xF if half == 0 else b & 0xF


@nb.njit(cache=True)
def _decode_int(buf, pos, half):
    head = _nibble(buf, pos, half)
    if half == 1:
        pos += 1
    half = 1 - half
    res = np.int64(0)
    if head <= 8:
        n = head
    else:
        n = head - 8
        for i in range(n):
            res |= np.int64(0xF0000000) >> (4 * i)
    for i in range(n, 8):
        hb = _nibble(buf, pos, half)
        if half == 1:
            pos += 1
        half = 1 - half
        res |= np.int64(hb) << ((i - n) * 4)
    res = np.int64(np.int32(np.uint32(res & 0xFFFFFFFF)))
    return res, pos, half


@nb.njit(cache=True)
def _fixed_point(buf):
    bits = np.uint64(0)
    for i in range(8):
        bits = (bits << np.uint64(8)) | np.uint64(buf[i])
    return np.array([bits]).view(np.float64)[0]


@nb.njit(cache=True)
def _int32_le(buf, start):
    v = np.uint32(0)
    for i in range(4):
        v |= np.uint32(buf[start + i]) << np.uint32(8 * i)
    return np.int64(np.int32(v))


@nb.njit(cache=True)
def _decode_linear(buf):
    size = buf.shape[0]
    if size < 12:
        return np.empty(0)
    fixed = _fixed_point(buf)
    out = np.empty(max(2, (size - 16) * 2 + 2))
    i1 = _int32_le(buf, 8)
    out[0] = i1 / fixed
    if size < 16:
        return out[:1]
    i2 = _int32_le(buf, 12)
    out[1] = i2 / fixed
    n, pos, half = 2, 16, 0
    while pos < size:
        if pos == size - 1 and half == 1 and (buf[pos] & 0xF) == 0:
            break
        diff, pos, half = _decode_int(buf, pos, half)
        i0, i1 = i1, i2
        i2 = i1 + (i1 - i0) + diff
        out[n] = i2 / fixed
        n += 1
    return out[:n]


@nb.njit(cache=True)
def _decode_pic(buf):
    size = buf.shape[0]
    out = np.empty(size * 2)
    n, pos, half = 0, 0, 0
    while pos < size:
        if pos == size - 1 and half == 1 and (buf[pos] & 0xF) == 0:
            break
        v, pos, half = _decode_int(buf, pos, half)
        out[n] = v
        n += 1
    return out[:n]


@nb.njit(cache=True)
def _decode_slof(buf):
    fixed = _fixed_point(buf)
    n = (buf.shape[0] - 8) // 2
    out = np.empty(n)
    for i in range(n):
        x = np.int64(buf[8 + 2 * i]) | (np.int64(buf[9 + 2 * i]) << 8)
        out[i] = np.exp(x / fixed) - 1.0
    return out
