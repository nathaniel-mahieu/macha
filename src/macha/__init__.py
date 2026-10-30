"""macha: modular LC-MS feature detection.

Modules
-------
macha.io        streaming mzML reading (local, gzip, remote)
macha.synth     ground-truth synthetic LC-MS generator (metabolomics & proteomics)
macha.traces    trace (ROI) detection by multi-target tracking
macha.baseline  censoring-aware baseline and noise estimation
macha.peaks     persistence peak detection and EMG fitting
macha.annotate  modulo mass-difference graph with Bayesian edge evaluation
macha.align     drift-trend mixture alignment with probabilistic correspondence
macha.pipeline  convenience end-to-end wrappers
"""

from .run import MSRun, stitch_scans

__version__ = "0.2.0"
__all__ = ["MSRun", "stitch_scans"]
