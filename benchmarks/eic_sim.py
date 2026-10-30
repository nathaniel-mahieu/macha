"""Ground-truth single-trace simulator for the baseline and peak benchmarks."""

from __future__ import annotations

import numpy as np

from macha.synth import emg_height_normalised, emg_mode
from macha.synth.lc import smooth_random_function

NOISE = 1_000.0
THRESH = 3_000.0


def background(rng, kind, n):
    x = np.arange(n, dtype=float)
    if kind == "none":
        return np.zeros(n)
    if kind == "constant":
        return np.full(n, rng.uniform(5, 30) * NOISE)
    if kind == "drift":
        f = smooth_random_function(rng, 0, n, 5, 1.0)
        return np.maximum(rng.uniform(2, 12) * NOISE * (1.2 + f(x)), 0)
    if kind == "step":
        lvl = np.where(x > rng.uniform(0.3, 0.7) * n, rng.uniform(6, 25), rng.uniform(0, 2)) * NOISE
        return np.convolve(lvl, np.ones(15) / 15, mode="same")
    raise ValueError(kind)


def make_eic(rng, n=400, kind="drift", n_peaks=None, snr_range=(5, 200), sigma_range=(2.0, 5.0), tau_ratio=(0.1, 1.0),
             min_sep=None):
    x = np.arange(n, dtype=float)
    base = background(rng, kind, n)
    k = rng.integers(1, 6) if n_peaks is None else n_peaks
    peaks = []
    y = base.copy()
    centers = np.sort(rng.uniform(30, n - 40, k))
    for c in centers:
        s = rng.uniform(*sigma_range)
        tau = s * rng.uniform(*tau_ratio)
        h = rng.uniform(np.log(snr_range[0]), np.log(snr_range[1]))
        h = np.exp(h) * NOISE
        y += h * emg_height_normalised(x, c, s, tau)
        peaks.append(dict(apex=float(emg_mode(c, s, tau)), height=h, sigma=s, tau=tau, snr=h / NOISE))
    sd = np.sqrt(NOISE**2 + (0.05 * y) ** 2)
    obs = y + rng.normal(0, 1, n) * sd
    obs[obs < THRESH] = np.nan
    return dict(x=x, y=obs, base=base, signal=y, peaks=peaks)


def make_pair(rng, sep_sigma, n=200, snr=60.0, ratio=0.5):
    """Two EMG peaks ``sep_sigma`` sigmas apart (resolution ~ sep_sigma / 4)."""
    x = np.arange(n, dtype=float)
    s = 3.0
    c1 = 80.0
    c2 = c1 + sep_sigma * s
    y = snr * NOISE * emg_height_normalised(x, c1, s, 0.3 * s) + ratio * snr * NOISE * emg_height_normalised(
        x, c2, s, 0.3 * s)
    obs = y + rng.normal(0, NOISE, n)
    obs[obs < THRESH] = np.nan
    return dict(x=x, y=obs, base=np.zeros(n), apexes=[float(emg_mode(c1, s, 0.3 * s)), float(emg_mode(c2, s, 0.3 * s))])
