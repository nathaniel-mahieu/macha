"""An Orbitrap-like measurement model for centroided MS1 data."""

from __future__ import annotations

from dataclasses import dataclass

import numba as nb
import numpy as np


@dataclass
class Instrument:
    resolution_200: float = 120_000.0  # FWHM resolving power at m/z 200
    noise: float = 1_000.0  # additive noise sigma (intensity units)
    threshold_snr: float = 3.0  # centroids below threshold_snr * noise are not reported
    ppm_floor: float = 0.25  # m/z precision at high signal-to-noise
    ppm_snr1: float = 12.0  # sigma_ppm = sqrt(floor^2 + (ppm_snr1 / snr)^2)
    rel_noise: float = 0.06  # multiplicative intensity noise
    shot_noise: float = 1.5  # variance term proportional to intensity
    merge_fraction: float = 0.9  # ions closer than this * FWHM coalesce into one centroid
    mz_range: tuple = (70.0, 1200.0)

    def resolution(self, mz):
        return self.resolution_200 * np.sqrt(200.0 / np.asarray(mz))

    def fwhm(self, mz):
        mz = np.asarray(mz)
        return mz / self.resolution(mz)

    @property
    def threshold(self) -> float:
        return self.threshold_snr * self.noise

    def ppm_sigma(self, intensity):
        snr = np.maximum(np.asarray(intensity) / self.noise, 1e-3)
        return np.sqrt(self.ppm_floor**2 + (self.ppm_snr1 / snr) ** 2)

    def measure(self, rng, scan, mz, intensity, rt, calibration_ppm=None):
        """Add intensity and m/z noise and apply the detection threshold.

        Returns observed (mz, intensity) and the boolean mask of reported points.
        """
        sd = np.sqrt(self.noise**2 + (self.rel_noise * intensity) ** 2 + self.shot_noise * intensity)
        obs_i = intensity + rng.normal(size=intensity.shape) * sd
        ppm = self.ppm_sigma(np.maximum(obs_i, self.noise)) * rng.normal(size=mz.shape)
        if calibration_ppm is not None:
            ppm = ppm + calibration_ppm(rt[scan], mz)
        obs_mz = mz * (1.0 + 1e-6 * ppm)
        keep = obs_i > self.threshold
        return obs_mz, obs_i, keep


@nb.njit(cache=True)
def coalesce(scan, mz, intensity, label, merge_tol_rel):
    """Merge same-scan ions closer than ``merge_tol_rel(mz) * mz``.

    Inputs must be sorted by (scan, mz).  ``merge_tol_rel`` is an array giving
    the relative tolerance at each point.  Returns arrays for the merged
    centroids plus a flag marking merged points.  The label of the most
    intense contributor is kept.
    """
    n = scan.shape[0]
    out_scan = np.empty(n, scan.dtype)
    out_mz = np.empty(n)
    out_i = np.empty(n)
    out_lab = np.empty(n, label.dtype)
    merged = np.zeros(n, np.bool_)
    k = -1
    best = 0.0
    for i in range(n):
        if k >= 0 and scan[i] == out_scan[k] and (mz[i] - out_mz[k]) < merge_tol_rel[i] * mz[i]:
            tot = out_i[k] + intensity[i]
            out_mz[k] = (out_mz[k] * out_i[k] + mz[i] * intensity[i]) / tot
            out_i[k] = tot
            if intensity[i] > best:
                best = intensity[i]
                out_lab[k] = label[i]
            merged[k] = True
        else:
            k += 1
            out_scan[k] = scan[i]
            out_mz[k] = mz[i]
            out_i[k] = intensity[i]
            out_lab[k] = label[i]
            best = intensity[i]
    k += 1
    return out_scan[:k], out_mz[:k], out_i[:k], out_lab[:k], merged[:k]
