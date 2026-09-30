"""
Add lrr_m_yr and lrr_r2 to run rate CSVs written before the LRR estimator existed.

    python scripts/input_prep/7-source-sink/1-prepare/backfill_run_lrr.py
    python scripts/input_prep/7-source-sink/1-prepare/backfill_run_lrr.py --check

Recomputes from each finished run's saved trajectory; --check reports
without writing. Details: scripts/input_prep/7-source-sink/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-22
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_REPO_ROOT = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO_ROOT / "scripts"))

from cascade_pipeline.run_layout import resolve  # noqa: E402
from cascade_pipeline.shoreline import compute_change_rate, compute_lrr  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW_RUNS = _REPO_ROOT / "output" / "raw_runs"

# The endpoint column must reproduce from the matrix to this tolerance for the pair to count as the same run
ENDPOINT_TOLERANCE_M_YR = 1e-9
# -----------------------------------------------------------------------------


# Every (rate_csv, shoreline_matrix) pair under a raw-runs tree
def find_pairs(root=RAW_RUNS):
    pairs = []
    for meta_path in sorted(root.rglob("*_run_metadata.json")):
        run_dir = meta_path.parent
        run_name = meta_path.name[: -len("_run_metadata.json")]
        csv_path = resolve(run_dir, "rate_csv", run_name)
        npy_path = resolve(run_dir, "matrix", run_name)
        if not csv_path.is_file():
            print(f"  SKIP  {run_dir.name}: no rate CSV")
        elif npy_path.is_file():
            pairs.append((csv_path, npy_path))
        else:
            print(f"  SKIP  {run_dir.name}: no shoreline matrix")
    return pairs


# Recomputes one run's LRR columns from its shoreline matrix
def backfill_one(csv_path, npy_path, geometry=HATTERAS_DOMAINS, check=False):
    frame = pd.read_csv(csv_path)
    matrix = np.load(npy_path)
    span = matrix.shape[0] - 1
    real = slice(geometry.start_real_index, geometry.end_real_index)

    # Guard: does the matrix reproduce the column the run actually wrote?
    endpoint = compute_change_rate(matrix, span_years=span)[real]
    if len(endpoint) != len(frame):
        return f"MISMATCH rows: csv {len(frame)}, matrix {len(endpoint)}"
    drift = float(np.nanmax(np.abs(endpoint - frame["change_rate_m_yr"].values)))
    if drift > ENDPOINT_TOLERANCE_M_YR:
        return (f"MISMATCH endpoint column differs from matrix by "
                f"{drift:.3e} m/yr -- CSV and matrix are from different runs")

    lrr, r2 = compute_lrr(matrix, span_years=span)
    if "lrr_m_yr" in frame.columns:
        existing = float(np.nanmax(np.abs(frame["lrr_m_yr"].values - lrr[real])))
        if existing <= ENDPOINT_TOLERANCE_M_YR:
            return "already"

    frame["lrr_m_yr"] = lrr[real]
    frame["lrr_r2"] = r2[real]
    if check:
        return "would write"
    frame.to_csv(csv_path, index=False)
    return "written"


# Run: every run missing the columns, or report with --check
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="report what would change and write nothing")
    parser.add_argument("--root", default=str(RAW_RUNS),
                        help="raw-runs tree to walk")
    args = parser.parse_args()

    pairs = find_pairs(Path(args.root))
    print(f"\n{len(pairs)} run(s) with both a rate CSV and a shoreline matrix\n")

    tally = {}
    for csv_path, npy_path in pairs:
        status = backfill_one(csv_path, npy_path, check=args.check)
        tally[status.split()[0]] = tally.get(status.split()[0], 0) + 1
        flag = "  " if status in ("written", "would write", "already") else "! "
        # The matrix stays at the run root in both layouts, so its parent is the run folder
        print(f"{flag}{npy_path.parent.name:58s} {status}")

    print("\n" + "  ".join(f"{k}={v}" for k, v in sorted(tally.items())))
    return 1 if any(k.startswith("MISMATCH") for k in tally) else 0


if __name__ == "__main__":
    raise SystemExit(main())
