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
    python duneline_to_raw_offsets.py --duneline duneline_1997.geojson \
        --out 1997_v2_duneline_offset_raw.csv \
        --validate-against 1997_duneline_offset_raw.csv

    --duneline   a file under 2-brie-offset/dunelines/, or a path
    --out        a file name under 2-brie-offset/raw_offsets/, or a path

EXTENSION MODE (2026-09-16, the Pea Island extension experiment)
    python duneline_to_raw_offsets.py --duneline duneline_1997.geojson --extension

    The 172 transects the surveyed polygon join left without a domain -- Pea
    Island north of GIS 90, and the last kilometre south of GIS 1 -- are given
    one by the SAME line-intersects-polygon join onto Hannah's whole-island
    polygons (hat_extension_domains.join_lines), and intersected the same
    way. A transect no polygon covers (the kilometre south of GIS 1, and the
    slivers between polygons) is dropped, as the surveyed join dropped it. Written
    to raw_offsets/ext/<vintage>_duneline_offset_raw_ext.csv, the same columns
    as the surveyed file, and the transect-to-domain table once to
    transects/transects_100m_ext.csv. The surveyed file is not touched.
"""

from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from site_layer.hat_extension_domains import EXTENDED_DOMAIN_POLYGONS, join_lines  # noqa: E402
# Resolved through hat_topo_version (2026-09-18).
from site_layer.hat_topo_version import (BRIE_ROOT, DUNELINE_DIR,  # noqa: E402
                              RAW_OFFSET_DIR as RAW_DIR, TRANSECT_EXT_TABLE,
                              TRANSECT_FILE_100M)
# Extension mode writes beside, never into, the surveyed raw files.
from site_layer.hat_topo_version import RAW_OFFSET_EXT_DIR as RAW_EXT_DIR  # noqa: E402

# The 100 m transects, 10 km long, each starting on the offshore datum line
# (x = 460198 in EPSG:3725) and running west across the island. domain_id is
# the ArcGIS spatial join onto the 500 m domain polygons; 450 of the 622
# transects fall in GIS 1-90, five per domain. Copied 2026-09-15 from
# hard-structures/groin/HAT-groin-gis-analysis/gis_data/, see transects/README.md.
TRANSECT_FILE = TRANSECT_FILE_100M

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


def load_transects(extension=False):
    """The 100 m transects with a domain number each.

    Surveyed mode: the 450 the ArcGIS polygon join placed in GIS 1-90.
    Extension mode: the 172 it left unplaced, numbered by their northing on
    their whole-island polygon (hat_extension_domains.join_lines); the table is
    also written to transects/transects_100m_ext.csv so the numbering is on
    disk beside the layer it extends.
    """
    t = gpd.read_file(TRANSECT_FILE)
    # The layer is an ArcGIS join export: every column is prefixed with the
    # table it came from ("Transects_100m.LineID"). Strip to the leaf name and
    # keep the first of any duplicates (OBJECTID and Shape_Length appear twice).
    t.columns = [c.split(".")[-1] for c in t.columns]
    t = t.loc[:, ~t.columns.duplicated()]
    t["LineID"] = t["LineID"].astype(int)
    if extension:
        t = t[t["domain_id"].isna()].copy()
        # A transect runs due west from the datum line, so either end's
        # northing is the transect's.
        t["northing_m"] = t.geometry.apply(lambda g: g.coords[0][1])
        # The same line-intersects-polygon join that placed the surveyed
        # transects, onto Hannah's whole-island polygons (2026-09-16). A
        # transect no polygon covers is dropped, as the surveyed join
        # dropped those in the slivers between polygons.
        t["domain_id"] = pd.Series(join_lines(t), index=t.index, dtype=float)
        dropped = t["domain_id"].isna()
        if dropped.any():
            print(f"  {int(dropped.sum())} transect(s) in no polygon skipped: LineID "
                  f"{t.loc[dropped, 'LineID'].tolist()}")
            t = t[~dropped].copy()
        t["domain_id"] = t["domain_id"].astype(int)
        t = t.sort_values(["domain_id", "LineID"]).reset_index(drop=True)
        print(f"  placed {len(t)} transects by polygon join onto "
              f"{EXTENDED_DOMAIN_POLYGONS.name}")
        TRANSECT_EXT_TABLE.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"LineID": t["LineID"], "domain_id": t["domain_id"],
                      "northing_m": t["northing_m"].round(3),
                      "crs": str(t.crs), "method": "polygon_join",
                      "polygons": EXTENDED_DOMAIN_POLYGONS.name}).to_csv(
            TRANSECT_EXT_TABLE, index=False)
        return t
    t = t[t["domain_id"].notna()].copy()
    t["domain_id"] = t["domain_id"].astype(int)
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
    ap.add_argument("--out", default=None,
                    help="required unless --extension, which has its own home")
    ap.add_argument("--validate-against", default=None,
                    help="a GIS-exported raw CSV of the SAME line, to check the method")
    ap.add_argument("--extension", action="store_true",
                    help="the transects beyond GIS 1-90, to raw_offsets/ext/")
    a = ap.parse_args(argv)

    line_path = _resolve(a.duneline, DUNELINE_DIR)
    if a.extension:
        if a.validate_against:
            ap.error("--validate-against has no GIS export to check in extension mode")
        vintage = re.search(r"duneline_(\d{4})", line_path.name)
        if a.out:
            out_path = Path(a.out) if Path(a.out).parent != Path(".") else RAW_EXT_DIR / a.out
        elif vintage:
            out_path = RAW_EXT_DIR / f"{vintage.group(1)}_duneline_offset_raw_ext.csv"
        else:
            ap.error(f"{line_path.name}: no duneline_<year> in the name; pass --out")
    elif a.out:
        out_path = Path(a.out) if Path(a.out).parent != Path(".") else RAW_DIR / a.out
    else:
        ap.error("--out is required (or --extension)")

    transects = load_transects(extension=a.extension)
    lines = gpd.read_file(line_path)
    if len(lines) != 1:
        sys.exit(f"{line_path.name}: expected one feature, found {len(lines)}")
    if lines.crs != transects.crs:
        print(f"  reprojecting {line_path.name} {lines.crs} -> {transects.crs}")
        lines = lines.to_crs(transects.crs)
    feat = lines.iloc[0]

    print(f"Dune line : {line_path}")
    print(f"Transects : {TRANSECT_FILE}  ({len(transects)} "
          + (f"beyond GIS {FIRST_DOMAIN}-{LAST_DOMAIN}, numbered "
             f"{transects.domain_id.min()}..{transects.domain_id.max()} by northing"
             if a.extension else f"in GIS {FIRST_DOMAIN}-{LAST_DOMAIN}") + ")")

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
