# ==============================================================================
# HAT_compare_rerun.py
#
# What changed between a stored run and the same run made under today's code?
#
# Pairs each run in an ARM against the stored run of the same name, and reports
# the differences that matter: the skill metrics, whether a road drowned, how
# many relocations fired, and the largest per-domain rate change.
#
# WHY IT IS NOT ENOUGH TO COMPARE SKILL. An island-wide RMSE can sit still
# while a road drowns: drowning stops roadway management for that domain from
# that year on, which changes what happens to the interior without necessarily
# moving the shoreline much. So the road table is differenced too.
#
#     python HAT_compare_rerun.py --arm recode-20260914
#
# Author: Hannah A. Henry, UNC CECL
# ==============================================================================

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
REPO = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
RAW = REPO / "output" / "raw_runs"
INDEX = RAW / "run_index.csv"

METRICS = ["mean_bias_interior_m_yr", "rmse_interior_m_yr",
           "roads_drowned", "roads_reloc_blocked"]


def rate_table(run_dir):
    p = run_dir / "tables" / "shoreline_change_rate.csv"
    return pd.read_csv(p) if p.is_file() else None


def road_table(run_dir):
    p = run_dir / "tables" / "road_management.csv"
    return pd.read_csv(p) if p.is_file() else None


def find_dir(name, arm, period, preset):
    if arm == "calibration":
        return RAW / period / preset / name
    return RAW / "arms" / arm / period / preset / name


def main():
    ap = argparse.ArgumentParser(description="stored run against its re-run")
    ap.add_argument("--arm", required=True)
    ap.add_argument("--full", action="store_true",
                    help="list every domain that differs, not just the count")
    args = ap.parse_args()

    index = pd.read_csv(INDEX)
    rerun = index[index["arm"].astype(str) == args.arm]
    if rerun.empty:
        print(f"no runs in arm {args.arm!r}")
        return 1

    print("=" * 96)
    print(f"STORED vs RE-RUN   arm {args.arm}")
    print("=" * 96)
    print(f"{'run':<46} {'bias':>16} {'RMSE':>16} {'drowned':>10}")

    moved, identical = [], []
    for _, new in rerun.iterrows():
        name = new["run_name"]
        old_rows = index[(index["run_name"] == name)
                         & (index["arm"].astype(str) == "calibration")]
        if old_rows.empty:
            print(f"{name:<46}   no stored run to compare")
            continue
        old = old_rows.iloc[-1]
        period = f"{int(old['start_year'])}_{int(old['end_year'])}"
        preset = str(old["source_sink_preset"])

        bits = []
        for col in METRICS:
            a, b = old.get(col), new.get(col)
            try:
                same = abs(float(a) - float(b)) < 1e-9
            except (TypeError, ValueError):
                same = a == b
            if not same:
                bits.append(f"{col}: {a} -> {b}")

        ra = rate_table(find_dir(name, "calibration", period, preset))
        rb = rate_table(find_dir(name, args.arm, period, preset))
        max_rate = None
        if ra is not None and rb is not None:
            m = ra.merge(rb, on="gis_domain", suffixes=("_a", "_b"))
            max_rate = (m["lrr_m_yr_a"] - m["lrr_m_yr_b"]).abs().max()

        da = road_table(find_dir(name, "calibration", period, preset))
        db = road_table(find_dir(name, args.arm, period, preset))
        road_diff = []
        if da is not None and db is not None:
            m = da.merge(db, on="gis", suffixes=("_a", "_b"))
            for col in ("drowned", "relocations"):
                d = m[m[f"{col}_a"] != m[f"{col}_b"]]
                if len(d):
                    road_diff.append(f"{col} at GIS {sorted(d.gis.tolist())}")

        print(f"{name:<46} "
              f"{float(old['mean_bias_interior_m_yr']):>7.4f} ->"
              f"{float(new['mean_bias_interior_m_yr']):>8.4f} "
              f"{float(old['rmse_interior_m_yr']):>7.4f} ->"
              f"{float(new['rmse_interior_m_yr']):>8.4f} "
              f"{int(old['roads_drowned']):>4} ->{int(new['roads_drowned']):>4}")
        if max_rate is not None and max_rate > 1e-9:
            print(f"{'':<46}   max per-domain rate change {max_rate:.4f} m/yr")
        for r in road_diff:
            print(f"{'':<46}   {r}")
        (moved if (bits or road_diff or (max_rate or 0) > 1e-9)
         else identical).append(name)

    print("\n" + "=" * 96)
    print(f"{len(identical)} run(s) reproduce exactly, {len(moved)} moved")
    if moved:
        print("\nMOVED. The stored runs and today's code disagree; the stored "
              "run is a faithful record of the code it was made with, not a "
              "wrong answer. Decide deliberately which to publish.")
        for n in moved:
            print("   ", n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
