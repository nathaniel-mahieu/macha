"""Download public Orbitrap Astral runs and keep only their MS1 centroids.

Raw mzML is streamed and parsed on the fly; only compact ``.npz`` MS1 files
are written to ``data/astral/<dataset>/``.

Datasets
--------
metabolomics  MassIVE MSV000093526 (Dorrestein lab): NIST fecal reference
              material, RP positive, Orbitrap Astral DDA (MS1 120k).  Sa04 is
              injected repeatedly across the batch, so its replicates carry
              real retention-time drift.
proteomics_sc ProteoBench module 9 / PRIDE PXD049412: 250 pg HeLa+yeast,
              Orbitrap Astral nDIA, MS1 240k, two conditions x 3 replicates.
proteomics    ProteoBench module 7: 50 ng HeLa/yeast/E. coli LFQ, 15 min
              gradient, Orbitrap Astral DIA.

Usage::

    python -I scripts/fetch_data.py                 # everything
    python -I scripts/fetch_data.py metabolomics    # one dataset
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "astral"

MASSIVE = "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000093526/peak/mzML/{}&forceDownload=true"
PB9 = "https://archive.openms.de/openms/benchmarks/pride-benchmarks/dia/OrbitrapAstral/ProteoBench_Module_9/20231123_DIA_240k_20Th_40ms_FAIMSCV-48_gas3p8_{}.mzML.gz"
PB7 = "https://archive.openms.de/openms/benchmarks/pride-benchmarks/dia/OrbitrapAstral/ProteoBench_Module_7/LFQ_Astral_DIA_15min_50ng_{}.mzML.gz"

DATASETS = {
    "metabolomics": {
        "description": "MSV000093526 NIST fecal RM, RP positive, Orbitrap Astral (MS1 120k)",
        "files": {
            "Sa04_inj005": MASSIVE.format("005_Sa04_Water_POS.mzML"),
            "Sa09_inj006": MASSIVE.format("006_Sa09_Water_POS.mzML"),
            "Sa17_inj007": MASSIVE.format("007_Sa17_Water_POS.mzML"),
            "Sa04_inj018": MASSIVE.format("018_Sa04_Water_POS.mzML"),
            "Sa04_inj023": MASSIVE.format("023_Sa04_Water_POS.mzML"),
            "Blank_inj001": MASSIVE.format("001_Blk_Water_POS.mzML"),
            "Sa04_inj005_multirange": MASSIVE.format("MultiMSranges/005_Sa04_Water_POS_Multi-MS1-Ranges.mzML"),
        },
    },
    "proteomics_sc": {
        "description": "PXD049412 / ProteoBench M9: 250 pg HeLa+yeast, Orbitrap Astral nDIA, MS1 240k",
        "files": {
            f"{cond}_r{r}": PB9.format(f"{pg}_H_Y_r{r}")
            for cond, pg in (("A", "240pg_10pg"), ("B", "200pg_50pg"))
            for r in (1, 2, 3)
        },
    },
    "proteomics": {
        "description": "ProteoBench M7: 50 ng LFQ, 15 min gradient, Orbitrap Astral DIA",
        "files": {
            "A_r1": PB7.format("Condition_A_REP1"),
            "A_r2": PB7.format("Condition_A_REP2"),
            "B_r1": PB7.format("Condition_B_REP1"),
        },
    },
}


def fetch(dataset: str) -> None:
    from macha.io import fetch_ms1

    spec = DATASETS[dataset]
    out = OUT / dataset
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = out / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest["description"] = spec["description"]
    files = manifest.setdefault("files", {})
    for name, url in spec["files"].items():
        dest = out / f"{name}.npz"
        if dest.exists():
            print(f"[{dataset}] {name}: present")
            continue
        print(f"[{dataset}] {name}: streaming {url}")
        t0 = time.time()
        massive = "massive.ucsd.edu" in url
        run = fetch_ms1(url, retries=14 if massive else 6, backoff=120.0 if massive else 30.0)
        run.meta.update(dataset=dataset, name=name)
        run.save(dest)
        files[name] = {
            "url": url, "n_scans": run.n_scans, "n_points": run.n_points,
            "rt_range": [float(run.rt[0]), float(run.rt[-1])],
            "bytes_streamed": run.meta.get("bytes_streamed"), "seconds": round(time.time() - t0, 1),
        }
        manifest_path.write_text(json.dumps(manifest, indent=2))
        print(f"[{dataset}] {name}: {run}")
        if "massive" in url:
            time.sleep(20)  # be gentle with MassIVE's rate limiter


if __name__ == "__main__":
    targets = sys.argv[1:] or list(DATASETS)
    for t in targets:
        fetch(t)
