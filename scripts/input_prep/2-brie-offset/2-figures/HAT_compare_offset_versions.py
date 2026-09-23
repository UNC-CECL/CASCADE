"""
Two builds of one start year: where did the island offset move?
===============================================================

Compares two versions under data/hatteras_init/2-brie-offset/<year>/ (the
unpadded 90-domain files each version's island_offset_hybrid.py run wrote) and,
when both versions have a raw per-transect file, the ABSOLUTE distances behind
them. The unpadded files are each zeroed on their own minimum, so their
difference is the change the MODEL sees; the raw files share the offshore
datum, so their difference is where the dune line was actually moved.

Written 2026-09-15 for 1996 v1 (the ArcGIS intersection of duneline_1997) vs
v2 (the shapely intersection of duneline_1997_v2, local corrections only).
The raw comparison uses the v1 line re-intersected by the SAME shapely script,
so the 1 m station convention of the GIS export (see
duneline_to_raw_offsets.py) does not appear as a change.

Outputs, in <year>/<b>/:
    offset_<year>_<a>_vs_<b>.csv    per domain: both versions, both frames
    offset_<year>_<a>_vs_<b>.png/.pdf   two panels, caption in CAPTIONS.md

USAGE
    python HAT_compare_offset_versions.py --year 1996 --a v1 --b v2 \
        --raw-a <path to the v1 line's shapely raw> --raw-b 1997_v2_duneline_offset_raw.csv

    # a superseded build: --a is its folder, --label-a what the outputs call it
    python HAT_compare_offset_versions.py --year 1996 \
        --a superseded_20260919_pre-redigitized/v2 --label-a superseded_v2 --b v1 \
        --raw-a <its raw> --raw-b <v1's raw>

--label-a/--label-b (2026-09-23) name a build in the file stem, the columns,
the legend and the caption. They default to --a/--b; they exist because a
superseded build's folder is a path, and a slash cannot go in a file name.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from site_layer.hat_figure_style import (  # noqa: E402
    C_1984, C_1997, DOMAIN_AXIS_LABEL, INK_MUTED, _title, apply_style,
    caption, figsize, save, town_bands)

import sys as _tvsys
from pathlib import Path as _TVP
_tvsys.path.insert(0, str(next(_q for _q in _TVP(__file__).resolve().parents
                               if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_topo_version as _tv  # noqa: E402
BRIE_ROOT = _tv.BRIE_ROOT
RAW_DIR = _tv.RAW_OFFSET_DIR


def _unpadded(year, version, source=None):
    p = _tv.offset_file(year, "unpadded", version=version,
                        source=source or _tv.DEFAULT_OFFSET_SOURCE)
    df = pd.read_csv(p)
    return df.set_index("Domain_ID")[str(year)]


def _raw_domain_means(path):
    raw = pd.read_csv(path)
    per_transect = raw.drop_duplicates(["domain_id", "LineID"])
    return per_transect.groupby("domain_id")["ORIG_LEN"].mean()


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--source", default=_tv.DEFAULT_OFFSET_SOURCE,
                    choices=_tv.OFFSET_SOURCES,
                    help="which source's versions to compare "
                         f"(default {_tv.DEFAULT_OFFSET_SOURCE})")
    ap.add_argument("--a", required=True, help="the earlier version, e.g. v1")
    ap.add_argument("--b", required=True, help="the later version, e.g. v2")
    ap.add_argument("--label-a", default=None,
                    help="what the outputs call --a (default: --a itself)")
    ap.add_argument("--label-b", default=None,
                    help="what the outputs call --b (default: --b itself)")
    ap.add_argument("--raw-a", default=None, help="raw per-transect CSV behind --a")
    ap.add_argument("--raw-b", default=None, help="raw per-transect CSV behind --b")
    args = ap.parse_args(argv)
    la = args.label_a or args.a
    lb = args.label_b or args.b
    for lab in (la, lb):
        if "/" in lab or "\\" in lab:
            ap.error(f"{lab!r} would put a folder in the file name; pass --label-a/--label-b")

    ua = _unpadded(args.year, args.a, args.source)
    ub = _unpadded(args.year, args.b, args.source)
    out = pd.DataFrame({f"model_{la}_m": ua, f"model_{lb}_m": ub})
    out["model_diff_m"] = ub - ua

    have_raw = bool(args.raw_a and args.raw_b)
    if have_raw:
        ra = _raw_domain_means(Path(args.raw_a) if Path(args.raw_a).exists() else RAW_DIR / args.raw_a)
        rb = _raw_domain_means(Path(args.raw_b) if Path(args.raw_b).exists() else RAW_DIR / args.raw_b)
        out[f"abs_{la}_m"] = ra.reindex(out.index)
        out[f"abs_{lb}_m"] = rb.reindex(out.index)
        # + = the line moved LANDWARD (a larger station from the offshore datum)
        out["abs_diff_m"] = out[f"abs_{lb}_m"] - out[f"abs_{la}_m"]
    out.index.name = "gis_domain"

    # Under the SOURCE's folder since 2026-09-22 -- this joined <year>/ and the
    # version by hand, which after the split would have written the comparison
    # into a <year>/v2/ that no longer exists.
    out_dir = _tv.offset_build_dir(args.year, args.b, args.source)
    stem = f"offset_{args.year}_{la}_vs_{lb}"
    out.to_csv(out_dir / f"{stem}.csv", float_format="%.2f")

    d = out["abs_diff_m"] if have_raw else out["model_diff_m"]
    moved = d[d.abs() >= 0.5]
    print(f"{args.year} {la} -> {lb}")
    if have_raw:
        print(f"  absolute (datum frame): {len(moved)} of {len(d)} domains moved >= 0.5 m; "
              f"mean {d.mean():+.2f} m, min {d.min():+.2f} (GIS {d.idxmin()}), "
              f"max {d.max():+.2f} (GIS {d.idxmax()})")
        print(f"  moved domains: {moved.round(1).to_dict()}")
        print(f"  baseline shift ({lb} - {la} of the per-year minimum): "
              f"{(out['model_diff_m'] - out['abs_diff_m']).mean():+.2f} m")
    m = out["model_diff_m"]
    print(f"  model frame (each zeroed on its own minimum): mean {m.mean():+.2f} m, "
          f"range {m.min():+.2f} .. {m.max():+.2f}")

    # ---- figure ---------------------------------------------------------- #
    apply_style()
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", aspect=0.62),
                             sharex=True, constrained_layout=True)
    x = out.index.to_numpy()

    ax = axes[0]
    ax.plot(x, out[f"model_{la}_m"], color=C_1984, lw=1.2, label=la)
    ax.plot(x, out[f"model_{lb}_m"], color=C_1997, lw=1.2, label=lb)
    ax.set_ylabel("Island offset (m)")
    _title(ax, 0, f"Island offset read by the model, {args.year} start")
    ax.legend(loc="lower left")   # the profile is high on the left, and the village labels sit at the top

    ax = axes[1]
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    if have_raw:
        ax.bar(x, out["abs_diff_m"], color=C_1997, width=0.8,
               label=f"{lb} − {la}, fixed datum")
        ax.set_ylabel("Dune line moved (m, + landward)")
        _title(ax, 1, "Where the re-digitised line differs")
    else:
        ax.bar(x, m, color=C_1997, width=0.8)
        ax.set_ylabel(f"{lb} − {la} (m)")
        _title(ax, 1, "Difference in the model frame")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_xlim(x.min() - 0.5, x.max() + 0.5)
    for a_ in axes:
        town_bands(a_, where="top", strip=0.08, label=(a_ is axes[0]))

    cap = (f"Two builds of the {args.year} island offset. (a) The unpadded offset each "
           f"build hands the model, {la} in red and {lb} in blue, each zeroed on "
           f"its own most seaward domain. ")
    if have_raw:
        cap += (f"(b) The change in the dune line itself, {lb} minus {la}, measured "
                f"from the shared offshore datum along the 100 m transects and averaged per "
                f"500 m domain; positive is landward. {len(moved)} of {len(d)} domains "
                f"differ by 0.5 m or more (mean over all domains {d.mean():+.1f} m; the "
                f"largest, {d[d.abs().idxmax()]:+.1f} m, at GIS {d.abs().idxmax()}). Both raw "
                f"files were produced by the same shapely intersection, so the metre-scale "
                f"station convention of the earlier ArcGIS export is not part of the difference.")
    else:
        cap += f"(b) {lb} minus {la} in that model frame."
    caption(fig, cap)
    paths = save(fig, out_dir / stem, close=True)
    print(f"  wrote {out_dir / (stem + '.csv')}")
    for p in paths:
        print(f"  wrote {p}")


if __name__ == "__main__":
    main()
