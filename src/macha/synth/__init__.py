"""Ground-truth synthetic LC-MS data."""

from .generate import (
    Scenario,
    SyntheticBatch,
    SyntheticRun,
    hard_scenario,
    make_batch,
    make_compounds,
    make_ions,
    make_run,
    metabolomics_scenario,
    proteomics_scenario,
    simulate_run,
    true_warps,
)
from .instrument import Instrument
from .lc import emg, emg_height_normalised, emg_mode

__all__ = [
    "Scenario", "SyntheticRun", "SyntheticBatch", "Instrument", "make_run", "make_batch",
    "make_compounds", "make_ions", "simulate_run", "metabolomics_scenario", "proteomics_scenario",
    "hard_scenario", "true_warps", "emg", "emg_mode", "emg_height_normalised",
]
