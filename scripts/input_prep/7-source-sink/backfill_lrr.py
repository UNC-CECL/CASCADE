"""Adds lrr_m_yr / lrr_r2 to run rate CSVs written before the LRR estimator.

WHY THIS EXISTS RATHER THAN A RE-RUN. Every finished run already carries the
whole trajectory it was scored from: `*_shoreline_matrix.npy` is the
[state x padded domain] array that `compute_change_rate` reduced to a single
column. The LRR is a second reduction of that same array, so it can be
recovered exactly -- bit for bit identical to what the run would have written
had the estimator existed at the time -- without re-simulating anything.

It is a backfill and not a fix: `change_rate_m_yr` is left exactly as the run
wrote it, and this script refuses to touch a CSV whose endpoint column does
not reproduce from the matrix. If those two disagree, the CSV and the matrix
came from different runs, and guessing which one is authoritative is not this
script's job.

Usage:
    python scripts/input_prep/7-source-sink/backfill_lrr.py [--check]

    --check  report what would change and write nothing.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_REPO_ROOT = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO_ROOT / "scripts"))

from cascade_pipeline.shoreline import compute_change_rate, compute_lrr  # noqa: E402
from hatteras_site_config import HATTERAS_DOMAINS  # noqa: E402

RAW_RUNS = _REPO_ROOT / "output" / "raw_runs"

# The endpoint column must reproduce from the matrix to this tolerance for the
# pair to count as the same run. Both sides are float64 reductions of the same
# array, so the only expected difference is CSV round-tripping.
ENDPOINT_TOLERANCE_M_YR = 1e-9


def find_pairs(root=RAW_RUNS):
    """Every (rate_csv, shoreline_matrix) pair under a raw-runs tree.

    Args:
        root: Directory to walk. Period and preset nesting is not assumed --
            the two files are matched by living in the same directory.

    Returns:
        A list of (csv_path, npy_path) tuples, sorted by run directory.
    """
    pairs = []
    for csv_path in sorted(root.rglob("*_shoreline_change_rate.csv")):
        stem = csv_path.name[:-len("_shoreline_change_rate.csv")]
        npy_path = csv_path.parent / f"{stem}_shoreline_matrix.npy"
        if npy_path.exists():
            pairs.append((csv_path, npy_path))
        else:
            print(f"  SKIP  {csv_path.parent.name}: no shoreline matrix")
    return pairs


def backfill_one(csv_path, npy_path, geometry=HATTERAS_DOMAINS, check=False):
    """Recomputes one run's LRR columns from its shoreline matrix.

    Args:
        csv_path: The run's `*_shoreline_change_rate.csv`.
        npy_path: The run's `*_shoreline_matrix.npy`.
        geometry: DomainGeometry, for the real-domain slice.
        check: Report only; do not write.

    Returns:
        A status string: "written", "would write", "already", or a reason the
        run was left alone.
    """
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
        print(f"{flag}{csv_path.parent.name:58s} {status}")

    print("\n" + "  ".join(f"{k}={v}" for k, v in sorted(tally.items())))
    return 1 if any(k.startswith("MISMATCH") for k in tally) else 0


if __name__ == "__main__":
    raise SystemExit(main())
