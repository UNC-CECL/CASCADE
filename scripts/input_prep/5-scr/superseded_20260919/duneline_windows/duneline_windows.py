"""
RETIRED 2026-09-19: its figure duplicated 3-rates/duneline/endpoint/1996_2024
and the duneline_endpoint chain figure; its output is in
5-scr/archive/2026-09-19_4-comparisons_duplicates/duneline_windows/.

duneline_windows.py
==============================================================================
Net dune-line change per GIS domain over a long window and its two halves,
in METRES -- the endpoint mirror of the CoastSat halves figure
(coastsat_lrr_windows.py --overlay), built 2026-09-18 (Hannah, by interview).

WHAT IS MEASURED
    Three digitized dune lines, one per period year through
    hat_topo_version.DUNE_LINE_FOR_YEAR: 1996 -> the 1997 line, 2010 -> 2009,
    2024 -> 2023. Each line's per-transect stations are read from
    2-brie-offset/raw_offsets/<vintage>_duneline_offset_raw.csv exactly as the
    hindcast loader reads them: the first row per transect, ORIG_LEN (distance
    from the offshore datum, growing LANDWARD), averaged over the ~5
    transects of each 500 m domain. The change is end minus start with the
    sign flipped, so SEAWARD IS POSITIVE.

    NET POSITION CHANGE, NOT A RATE AND NOT A FIT (Hannah: "we are only
    looking at net position change"). No survey dates are needed, which
    matters because the 2023 imagery date is unknown. And because the two
    halves share the middle line, they ADD UP to the whole exactly, domain by
    domain; the script checks it. Unlike an LRR, the whole therefore always
    lies between (or at the sum of) its parts -- a storm jump cannot hide at
    the join.

WHAT IS DRAWN (the CoastSat figure's layout, on purpose, so the two sit side
by side)
    (a) the whole window, filled blue seaward and red landward;
    (b) the two halves, the earlier grey, the later black.
    Model-input fill footprints as bars above (a), the offshore shoals as
    faint hatched boxes, village bands and structure hairlines. Labels use
    the REAL line years (1997, 2009, 2023); the caption says what they stand
    in for. y axis: the tightest multiple of 10 m holding every line.

OUTPUT   data/hatteras_init/5-scr/4-comparisons/duneline_windows/<start>_<end>/
    duneline_change_<a>_<b>_halves.png     also published to
                                           output/figures/shoreline/
    supporting/
        duneline_change_<a>_<b>_halves.csv the three positions and three
                                           changes per domain, plus the
                                           additivity residual
        <figure>.pdf, CAPTIONS.md

USAGE
    python scripts/input_prep/5-scr/duneline_windows/duneline_windows.py              # 1996 2010 2024
    python scripts/input_prep/5-scr/duneline_windows/duneline_windows.py --years 1996 2010 2024
==============================================================================
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
# The CoastSat halves figure's drawing code, reused so the two figures match.
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "CoastSat"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

import coastsat_lrr_windows as cw  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, _title, apply_style, caption, figsize, figure_dir,
    save, support_dir,
)
from site_layer.hat_observed_rates import DUNELINE_WINDOWS  # noqa: E402
from site_layer.hat_topo_version import dune_line_for_year, dune_raw_file  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402

N_DOMAINS = 90
Y_LABEL = "Net dune-line change (m)"
Y_STEP_M = 10.0      # the bound is a multiple of this
Y_TICK_M = 20.0


def dune_position(vintage: int) -> pd.Series:
    """Mean ORIG_LEN per GIS domain, first row per transect -- the hindcast
    loader's reading (duneline_vs_coastsat.dune_position_by_domain does the
    same). Grows LANDWARD."""
    raw = pd.read_csv(dune_raw_file(vintage), encoding="utf-8-sig")
    per_transect = raw.drop_duplicates(subset=["domain_id", "LineID"])
    return (per_transect.groupby("domain_id")["ORIG_LEN"].mean()
            .reindex(range(1, N_DOMAINS + 1)))


def as_frame(change: pd.Series) -> pd.DataFrame:
    """The shape the CoastSat drawing code reads (mean_lrr is the value)."""
    return pd.DataFrame({"domain_number": change.index.to_numpy(),
                         "mean_lrr": change.to_numpy(float),
                         "std_lrr": np.nan})


def caption_text(v0, v1, v2, years, half, fills) -> str:
    y0, y1, y2 = years
    fill_txt = "; ".join(f"{y} at GIS {lo}–{hi}" for y, lo, hi in fills)
    shoal_txt = "; ".join(f"{n} GIS {lo}–{hi}" for n, (lo, hi)
                          in HATTERAS_ANNOTATIONS.shoal_zones.items())
    return ("Net change in dune-line position by GIS domain (1 at Cape Point, "
            f"90 at Pea Island). (a) {v0}–{v2}: the {v2} digitized dune line "
            f"minus the {v0} line, each measured along the 100 m transects "
            "from a fixed offshore datum and averaged over the transects of "
            "each 500 m domain; blue and filled where the dune moved seaward, "
            f"red where it moved landward. (b) The same over {v0}–{v1} (grey) "
            f"and {v1}–{v2} (black); the two add up to (a) exactly, since they "
            f"share the {v1} line. The lines stand in for the model years "
            f"{y0}, {y1} and {y2} (no island-wide imagery at those dates). This "
            "is net displacement, not a rate: no survey dates enter it. Black "
            "bars above (a) mark the beach fills placed inside the window at "
            f"the footprint the hindcast uses ({fill_txt}); hatched amber "
            f"boxes mark the offshore shoals ({shoal_txt}). Village spans are "
            "shaded; the solid hairline is the Buxton groin and the dotted "
            "hairlines are the Avon and Rodanthe piers. Both panels share a y "
            f"axis of ±{half:g} m, the smallest multiple of {Y_STEP_M:g} m "
            "that holds every value.")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 3)[2])
    ap.add_argument("--years", type=int, nargs=3, default=(1996, 2010, 2024),
                    metavar=("START", "MID", "END"),
                    help="period years; each reads its line through "
                         "DUNE_LINE_FOR_YEAR")
    a = ap.parse_args(argv)
    y0, y1, y2 = a.years
    v0, v1, v2 = (dune_line_for_year(y) for y in a.years)

    p0, p1, p2 = dune_position(v0), dune_position(v1), dune_position(v2)
    # ORIG_LEN grows landward, so seaward-positive change is start minus end
    whole, first, second = p0 - p2, p0 - p1, p1 - p2
    resid = (first + second - whole).abs().max()
    assert resid < 1e-6, f"halves do not add up to the whole: {resid}"
    missing = int(whole.isna().sum())

    out = DUNELINE_WINDOWS / f"{y0}_{y2}"
    stem = f"duneline_change_{v0}_{v2}_halves"
    table = pd.DataFrame({
        "domain_number": whole.index,
        f"position_{v0}_m": p0.round(2), f"position_{v1}_m": p1.round(2),
        f"position_{v2}_m": p2.round(2),
        f"change_{v0}_{v1}_m": first.round(2),
        f"change_{v1}_{v2}_m": second.round(2),
        f"change_{v0}_{v2}_m": whole.round(2),
    })
    table.to_csv(support_dir(out) / f"{stem}.csv", index=False)

    apply_style()
    frames = [as_frame(whole), as_frame(first), as_frame(second)]
    extreme = max(float(np.nanmax(np.abs(f["mean_lrr"]))) for f in frames)
    half = float(math.ceil(extreme / Y_STEP_M) * Y_STEP_M)
    fills = cw.fills_in(v0, v2)

    fig, (ax_a, ax_b) = plt.subplots(
        2, 1, sharex=True, sharey=True, constrained_layout=True,
        figsize=figsize("double", height=4.9))
    cw.draw_panel(ax_a, frames[0], half, std=False)
    cw.draw_fills(ax_a, fills, half)
    cw.draw_shoals(ax_a, label=True)
    _title(ax_a, 0, f"{v0}–{v2}")
    handles = cw.draw_halves(ax_b, [(v0, v1), (v1, v2)], frames[1:], half)
    cw.draw_shoals(ax_b, label=False)
    _title(ax_b, 1, "The two periods")
    for ax in (ax_a, ax_b):
        ax.yaxis.set_major_locator(MultipleLocator(Y_TICK_M))
    ax_b.set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel(Y_LABEL, fontsize=9)
    fig.legend(handles=handles, loc="outside lower center", ncol=2,
               frameon=False)
    caption(fig, caption_text(v0, v1, v2, (y0, y1, y2), half, fills))
    written = save(fig, out / stem)
    written += save(fig, figure_dir("shoreline") / stem)
    plt.close(fig)

    print(f"lines     {v0} -> {v1} -> {v2}  (for {y0}, {y1}, {y2})")
    print(f"domains   {N_DOMAINS - missing}/{N_DOMAINS} with all three lines; "
          f"halves add to whole within {resid:.1e} m")
    for label, s in ((f"{v0}-{v2}", whole), (f"{v0}-{v1}", first),
                     (f"{v1}-{v2}", second)):
        print(f"{label}  island mean {s.mean():+7.1f} m, landward in "
              f"{int((s < 0).sum())}/90, range {s.min():+.0f} .. {s.max():+.0f}")
    print(f"y bounds  +/-{half:g} m")
    for p in written:
        print("wrote    ", p.relative_to(_REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
