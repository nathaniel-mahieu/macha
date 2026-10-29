"""Chromatographic peak shapes and smooth drift functions."""

from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import erfc, erfcx

SQRT2 = np.sqrt(2.0)


def emg(t, mu, sigma, tau):
    """Exponentially modified Gaussian density, numerically stable everywhere.

    Broadcasts over all arguments.  ``tau -> 0`` approaches a Gaussian.
    """
    t, mu, sigma, tau = np.broadcast_arrays(*(np.asarray(a, dtype=np.float64) for a in (t, mu, sigma, tau)))
    tau = np.maximum(tau, 1e-6)
    x = t - mu
    z = (sigma / tau - x / sigma) / SQRT2
    out = np.empty_like(x)
    pos = z >= 0
    out[pos] = np.exp(-0.5 * (x[pos] / sigma[pos]) ** 2) * erfcx(z[pos]) / (2 * tau[pos])
    neg = ~pos
    if neg.any():
        a = 0.5 * (sigma[neg] / tau[neg]) ** 2 - x[neg] / tau[neg]
        out[neg] = np.exp(a) * erfc(z[neg]) / (2 * tau[neg])
    return out


def emg_mode(mu, sigma, tau, n_iter: int = 30):
    """Location of the EMG maximum (Newton on the log-density derivative)."""
    mu, sigma, tau = (np.asarray(a, dtype=np.float64) for a in (mu, sigma, tau))
    t = mu + np.minimum(tau, sigma) * 0.5
    h = 1e-4 * sigma
    for _ in range(n_iter):
        f0 = np.log(emg(t, mu, sigma, tau) + 1e-300)
        fp = np.log(emg(t + h, mu, sigma, tau) + 1e-300)
        fm = np.log(emg(t - h, mu, sigma, tau) + 1e-300)
        d1 = (fp - fm) / (2 * h)
        d2 = (fp - 2 * f0 + fm) / h**2
        step = np.where(d2 < 0, -d1 / np.minimum(d2, -1e-12), d1 * sigma)
        t = t + np.clip(step, -sigma, sigma)
    return t


def emg_height_normalised(t, mu, sigma, tau):
    """EMG scaled so its apex equals 1."""
    mode = emg_mode(mu, sigma, tau)
    return emg(t, mu, sigma, tau) / emg(mode, mu, sigma, tau)


def smooth_random_function(rng: np.random.Generator, t0: float, t1: float, n_knots: int = 6,
                           amplitude: float = 1.0):
    """A random smooth function on [t0, t1] with values roughly in +/- amplitude."""
    knots = np.linspace(t0, t1, n_knots)
    vals = rng.uniform(-1, 1, n_knots) * amplitude
    cs = CubicSpline(knots, vals, bc_type="natural")

    def f(t):
        return cs(np.clip(t, t0, t1))

    return f
