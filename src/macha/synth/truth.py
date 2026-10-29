"""Map detected features back to synthetic ground truth."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .lc import emg_mode


def feature_truth(features: pd.DataFrame, traces, ion_of_point: np.ndarray, ions: pd.DataFrame,
                  compounds: pd.DataFrame) -> pd.DataFrame:
    """Majority true ion of the centroids inside each feature's peak window."""
    scan = traces.run.scan_of_point
    out_ion = np.full(len(features), -1, np.int64)
    purity = np.zeros(len(features))
    for k, (t, s0, s1) in enumerate(zip(features["trace"].to_numpy(), features["scan_start"].to_numpy(),
                                        features["scan_end"].to_numpy())):
        pts = traces.points(int(t))
        sc = scan[pts]
        pts = pts[(sc >= s0) & (sc <= s1)]
        lab = ion_of_point[pts]
        lab = lab[lab >= 0]
        if lab.size:
            vals, cnt = np.unique(lab, return_counts=True)
            out_ion[k] = vals[np.argmax(cnt)]
            purity[k] = cnt.max() / pts.size
    ion_tab = ions.set_index("ion_id")
    comp = compounds.set_index("compound_id")
    ok = out_ion >= 0
    df = pd.DataFrame(index=features.index)
    df["ion"] = out_ion
    df["purity"] = purity
    df["compound"] = -1
    df.loc[ok, "compound"] = ion_tab.loc[out_ion[ok], "compound_id"].to_numpy()
    for col in ("adduct", "charge", "isotope", "mz"):
        df[col] = None if col == "adduct" else np.nan
        df.loc[ok, col] = ion_tab.loc[out_ion[ok], col].to_numpy()
    df = df.rename(columns={"mz": "true_mz"})
    c = df["compound"].to_numpy()
    okc = c >= 0
    df["kind"] = ""
    df.loc[okc, "kind"] = comp.loc[c[okc], "kind"].to_numpy()
    df["true_rt"] = np.nan
    emgc = okc & (df["kind"] != "background").to_numpy()
    cc = comp.loc[c[emgc]]
    df.loc[emgc, "true_rt"] = emg_mode(cc["rt"].to_numpy(), cc["sigma"].to_numpy(), cc["tau"].to_numpy())
    df["neutral_mass"] = np.nan
    df.loc[okc, "neutral_mass"] = comp.loc[c[okc], "neutral_mass"].to_numpy()
    return df
