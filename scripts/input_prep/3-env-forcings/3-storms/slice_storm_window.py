"""
Cut one window's storm series out of the full 1996-2025 series, re-indexed so model year 1 is the
window's first storm year.

    python scripts/input_prep/3-env-forcings/3-storms/slice_storm_window.py --start-year 1997 --end-year 2009
    python scripts/input_prep/3-env-forcings/3-storms/slice_storm_window.py --start-year 2010 --end-year 2025

Both years inclusive: the window holds storms from 1 Jan start to 31 Dec end. The 1996_2009 and
2009_2025 files were checked to be exactly such slices of 1996_2025 (2026-10-12), so slicing gives
the same events, Rhigh, Rlow, period and duration as building the window from the records again.
Writes <start>_<end>_storms_<variant>.npy, .csv and _summary.csv into the window's hindcast_storms
folder.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-12
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer import hat_env_forcings as _env  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
SOURCE = (1996, 2025)                           # the full series every window is cut from
# -----------------------------------------------------------------------------


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--start-year", type=int, required=True, help="first storm year (model year 1)")
    ap.add_argument("--end-year", type=int, required=True, help="last storm year, inclusive")
    ap.add_argument("--variant", default=_env.DEFAULT_STORM_VARIANT)
    a = ap.parse_args(argv)
    if not SOURCE[0] <= a.start_year <= a.end_year <= SOURCE[1]:
        sys.exit(f"window {a.start_year}-{a.end_year} is outside the source series {SOURCE[0]}-{SOURCE[1]}")

    src = _env.storm_series_file(*SOURCE, a.variant)
    storms = np.load(src)
    summary = pd.read_csv(_env.storm_summary_file(*SOURCE, a.variant))
    first, last = a.start_year - SOURCE[0] + 1, a.end_year - SOURCE[0] + 1   # model years in the source
    keep = (storms[:, 0] >= first) & (storms[:, 0] <= last)
    cut = storms[keep].copy()
    cut[:, 0] -= first - 1
    s = summary[(summary.calendar_year >= a.start_year) & (summary.calendar_year <= a.end_year)].copy()
    s["time"] = s["time"] - (first - 1)
    if len(s) != len(cut) or not np.array_equal(s["time"].to_numpy(), cut[:, 0]):
        sys.exit("the summary and the series disagree after the cut")

    out = _env.storm_series_file(a.start_year, a.end_year, a.variant)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.save(out, cut)
    pd.DataFrame(cut, columns=["time", "Rhigh", "Rlow", "period", "duration"]).astype(
        {"time": int, "duration": int}).to_csv(out.with_suffix(".csv"), index=False)
    s.to_csv(_env.storm_summary_file(a.start_year, a.end_year, a.variant), index=False)
    years = a.end_year - a.start_year + 1
    print(f"{out.name}: {len(cut)} storms in {years} model years "
          f"({int(cut[:, 0].min())}-{int(cut[:, 0].max())}), cut from {src.name}")


if __name__ == "__main__":
    main()
