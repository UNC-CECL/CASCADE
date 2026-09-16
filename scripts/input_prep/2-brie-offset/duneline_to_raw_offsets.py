"""
Dune line -> raw per-transect offsets (the step ArcGIS used to do)
==================================================================

Every file in data/hatteras_init/2-brie-offset/raw_offsets/ was, until
2026-09-15, an ArcGIS export: buffer the digitised dune line by 1.5 m, intersect
the buffer with 1 m points generated along the 100 m transects, and export the
attribute table. The distance that matters is ORIG_LEN, the station of the
point along its transect measured from the transect's start on the offshore
datum line. island_offset_hybrid.py then keeps one point per transect, averages
the ~5 transects in each 500 m domain, and pads for CASCADE.

This script does the same intersection in shapely, so a re-digitised line can
be turned into a raw file from the repo alone, with no GIS session and no
external drive. One row per transect, holding the exact intersection station.

VALIDATED 2026-09-15 against the ArcGIS export of the v1 1997 line
(raw_offsets/1997_duneline_offset_raw.csv): all 450 transects in GIS 1-90
match, mean difference -1.01 m, sd 0.32 m, worst 1.9 m. The constant metre is
the GIS convention, not geometry: of the ~3 one-metre stations inside the
1.5 m buffer the export lists the LANDWARD-most first, and the downstream
scripts take the first row per transect. The exact intersection sits ~1 m
seaward of it. This cancels in island_offset_hybrid.py (each year is zeroed
on its own minimum) and cancels in an end-year difference of two files built
by THIS script; a difference between a GIS-built and a shapely-built file
carries the metre. Pass --validate-against to reproduce those numbers.

WHICH CROSSING when a transect meets the line more than once: the landward-most
(largest station), which is what the GIS first-row convention returned. The
count is written to n_crossings so those transects can be found.

USAGE
    python duneline_to_raw_offsets.py --duneline duneline_1997_v2.geojson \
        --out 1997_v2_duneline_offset_raw.csv \
        --validate-against 1997_duneline_offset_raw.csv

    --duneline   a file under 2-brie-offset/dunelines/, or a path
    --out        a file name under 2-brie-offset/raw_offsets/, or a path
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
BRIE_ROOT = INIT_ROOT / "2-brie-offset"
RAW_DIR = BRIE_ROOT / "raw_offsets"
DUNELINE_DIR = INIT_ROOT / "2-brie-offset" / "dunelines"

# The 100 m transects, 10 km long, each starting on the offshore datum line
# (x = 460198 in EPSG:3725) and running west across the island. domain_id is
# the ArcGIS spatial join onto the 500 m domain polygons; 450 of the 622
# transects fall in GIS 1-90, five per domain. Copied 2026-09-15 from
# hard-structures/groin/HAT-groin-gis-analysis/gis_data/, see transects/README.md.
TRANSECT_FILE = BRIE_ROOT / "transects" / "transects_100m.geojson"

FIRST_DOMAIN, LAST_DOMAIN = 1, 90

# Metadata columns copied from the dune-line feature when it carries them
# (the 1997 lines do; 1984 and 1967 carry none).
LINE_META = ("feature_type", "year", "source_type", "method", "editor",
             "edit_date", "notes")


def _resolve(name_or_path, default_dir):
    p = Path(name_or_path)
    if p.exists():
        return p.resolve()
    q = default_dir / name_or_path
    if q.exists():
        return q.resolve()
    raise FileNotFoundError(f"{name_or_path}: not a path, and not under {default_dir}")


def load_transects():
    t = gpd.read_file(TRANSECT_FILE)
    # The layer is an ArcGIS join export: every column is prefixed with the
    # table it came from ("Transects_100m.LineID"). Strip to the leaf name and
    # keep the first of any duplicates (OBJECTID and Shape_Length appear twice).
    t.columns = [c.split(".")[-1] for c in t.columns]
    t = t.loc[:, ~t.columns.duplicated()]
    t = t[t["domain_id"].notna()].copy()
    t["domain_id"] = t["domain_id"].astype(int)
    t["LineID"] = t["LineID"].astype(int)
    t = t[(t["domain_id"] >= FIRST_DOMAIN) & (t["domain_id"] <= LAST_DOMAIN)]
    return t.sort_values(["domain_id", "LineID"]).reset_index(drop=True)


def intersect(transects, line):
    """One row per transect: station of the landward-most crossing, or NaN."""
    rows = []
    for _, tr in transects.iterrows():
        geom = tr.geometry
        x = geom.intersection(line)
        if x.is_empty:
            rows.append((tr.domain_id, tr.LineID, np.nan, 0, np.nan, np.nan))
            continue
        pts = [x] if x.geom_type == "Point" else [g for g in x.geoms
                                                  if g.geom_type == "Point"]
        stations = [geom.project(p) for p in pts]
        k = int(np.argmax(stations))
        rows.append((tr.domain_id, tr.LineID, float(stations[k]), len(pts),
                     pts[k].x, pts[k].y))
    return pd.DataFrame(rows, columns=["domain_id", "LineID", "ORIG_LEN",
                                       "n_crossings", "x", "y"])


def validate(out_df, gis_path):
    """Compare per-transect stations with a GIS export (first row per transect)."""
    gis = pd.read_csv(gis_path)
    gis = gis.drop_duplicates(["domain_id", "LineID"])[["domain_id", "LineID", "ORIG_LEN"]]
    gis["LineID"] = gis["LineID"].astype(int)
    m = out_df.merge(gis, on=["domain_id", "LineID"], how="left",
                     suffixes=("", "_gis"))
    m["diff_m"] = m["ORIG_LEN"] - m["ORIG_LEN_gis"]
    d = m["diff_m"].dropna()
    dom = m.groupby("domain_id")[["ORIG_LEN", "ORIG_LEN_gis"]].mean()
    dom_diff = dom["ORIG_LEN"] - dom["ORIG_LEN_gis"]
    print(f"\nValidation against {gis_path.name}:")
    print(f"  transects matched      : {d.size} of {len(m)}")
    print(f"  per-transect diff (m)  : mean {d.mean():+.2f}  sd {d.std():.2f}  "
          f"min {d.min():+.2f}  max {d.max():+.2f}")
    print(f"  per-domain mean diff   : mean {dom_diff.mean():+.2f}  "
          f"sd {dom_diff.std():.2f}  range {dom_diff.min():+.2f} .. {dom_diff.max():+.2f}")
    return m


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--duneline", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--validate-against", default=None,
                    help="a GIS-exported raw CSV of the SAME line, to check the method")
    a = ap.parse_args(argv)

    line_path = _resolve(a.duneline, DUNELINE_DIR)
    out_path = Path(a.out) if Path(a.out).parent != Path(".") else RAW_DIR / a.out

    transects = load_transects()
    lines = gpd.read_file(line_path)
    if len(lines) != 1:
        sys.exit(f"{line_path.name}: expected one feature, found {len(lines)}")
    if lines.crs != transects.crs:
        print(f"  reprojecting {line_path.name} {lines.crs} -> {transects.crs}")
        lines = lines.to_crs(transects.crs)
    feat = lines.iloc[0]

    print(f"Dune line : {line_path}")
    print(f"Transects : {TRANSECT_FILE}  ({len(transects)} in GIS "
          f"{FIRST_DOMAIN}-{LAST_DOMAIN})")

    df = intersect(transects, feat.geometry)
    for col in LINE_META:
        df[col] = feat[col] if col in lines.columns else ""
    df["duneline_file"] = line_path.name
    df["transect_file"] = TRANSECT_FILE.name
    df["built"] = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    df["built_by"] = Path(__file__).name

    missing = df[df["ORIG_LEN"].isna()]
    multi = df[df["n_crossings"] > 1]
    print(f"  transects with no crossing : {len(missing)}"
          + (f"  -> LineID {missing.LineID.tolist()}" if len(missing) else ""))
    print(f"  transects crossed >1 times : {len(multi)}"
          + (f"  -> LineID {multi.LineID.tolist()}" if len(multi) else ""))
    per_dom = df.groupby("domain_id")["ORIG_LEN"].count()
    short = per_dom[per_dom < 5]
    if len(short):
        print(f"  domains with <5 transects  : {short.to_dict()}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"\nWrote {out_path}  ({len(df)} rows)")

    if a.validate_against:
        gis_path = _resolve(a.validate_against, RAW_DIR)
        m = validate(df, gis_path)
        vpath = out_path.with_name(out_path.stem + f"_validation_vs_{gis_path.stem}.csv")
        m.to_csv(vpath, index=False)
        print(f"  per-transect comparison -> {vpath.name}")


if __name__ == "__main__":
    main()
