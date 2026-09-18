#!/usr/bin/env python3
"""
plot_dsas_calibration_periods.py
==============================================================================
Observed shoreline change rate per domain from the DSAS transect record --
the INDEPENDENT CHECK on the CoastSat figure, not the target.

CoastSat is what the model is graded against (see
coastsat_calibration_periods.png, and hat_observed_rates), so this one lives
under supporting/ and is drawn in exactly the same style, so the two can be
laid side by side and only the data differs (Hannah, 2026-09-17).

TWO THINGS THIS FIGURE CANNOT MATCH EXACTLY, and both are stated in the
caption rather than smoothed over:

  THE WINDOWS. DSAS has five digitised shorelines -- 1978, 1987, 1997, 2009,
  2019 -- so it cannot be cut to the run periods. The nearest pairs are used:
  1997-2009 stands in for 1996-2010 (one year in at each end) and 2009-2019
  for 2010-2024 (one year late, five years short). Until 2026-09-17 it drew
  1978-1997 and 1997-2019, which straddle the run periods rather than
  approximating them.

  THE ESTIMATOR. This is an END-POINT RATE: the first and last shoreline of
  the window, over the elapsed years. CoastSat is an OLS slope through every
  transect observation in the window, which is what the model is scored with
  (see hat_observed_rates and the LRR note in HAT_hindcast_methods). With two
  shorelines an OLS slope IS the end-point rate, so the two agree in form
  here; they would not if a third DSAS vintage fell inside a window.
==============================================================================
"""
import matplotlib
matplotlib.use("Agg")
import pandas as pd

# HOUSE STYLE: one typeface and one palette across every figure in this
# project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
# found by searching upward (ORGANIZATION.md rule 5), so this block is
# independent of whatever this script calls its own repository variable.
import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(next(_q for _q in _P(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer.hat_figure_style import (apply_style, C, C_1984, C_1997, C_1984_FILL,
                              C_1997_FILL, INK, INK_MUTED, GRID_C, figsize,
                              record_caption, town_bands, structures,
                              open_frame, DOMAIN_AXIS_LABEL, support_dir)
apply_style()
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from site_layer.hatteras_site_config import HATTERAS_PERIODS, HATTERAS_ANNOTATIONS

_REPO = next(_p for _p in _P(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
from site_layer.hat_observed_rates import DSAS_ROOT  # noqa: E402
INPUT_CSV = DSAS_ROOT / "All_Shoreline_Transect_Intersections.csv"
from site_layer import hat_figure_style as _hs  # noqa: E402
FIG_DIR = _hs.figure_dir("shoreline")

# The run periods this is checking, and the DSAS pair that stands in for each.
# Both are stated so a reader never has to infer the offset.
PERIOD_STARTS = (1996, 2010)
RUN_PERIODS = [(st, HATTERAS_PERIODS[st]["end_year"]) for st in PERIOD_STARTS]
DSAS_WINDOWS = ((1997, 2009), (2009, 2019))

TRANSECT_ID_COL = "Transects_100m_LineID"
DOMAIN_ID_COL = "Transects_100m_AddSpatialJoin_domain_id"
YEAR_COL, DISTANCE_COL = "Year", "NEAR_DIST"
DOMAIN_MIN, DOMAIN_MAX = 1, 90
PERIOD_COLOURS = ((C_1984, C_1984_FILL), (C_1997, C_1997_FILL))


def domain_rates():
    """Domain-mean end-point rate for each DSAS window, positive seaward."""
    frame = pd.read_csv(INPUT_CSV)
    frame = frame[[TRANSECT_ID_COL, DOMAIN_ID_COL, YEAR_COL, DISTANCE_COL]].dropna()
    domain_map = frame[[TRANSECT_ID_COL, DOMAIN_ID_COL]].drop_duplicates()
    # a shoreline can cross one transect twice; average those first
    frame = frame.groupby([TRANSECT_ID_COL, YEAR_COL], as_index=False)[DISTANCE_COL].mean()
    wide = frame.pivot(index=TRANSECT_ID_COL, columns=YEAR_COL, values=DISTANCE_COL)

    have = set(wide.columns)
    missing = [y for w in DSAS_WINDOWS for y in w if y not in have]
    if missing:
        raise SystemExit(
            f"\nDSAS shoreline year(s) {sorted(set(missing))} are not in "
            f"{INPUT_CSV.name}. Present: {sorted(have)}\n")

    rates = pd.DataFrame(index=wide.index)
    for i, (y1, y2) in enumerate(DSAS_WINDOWS, 1):
        # NEAR_DIST grows landward, so the sign flips to make + seaward
        rates[f"P{i}"] = -1.0 * (wide[y2] - wide[y1]) / (y2 - y1)

    joined = rates.merge(domain_map, left_index=True, right_on=TRANSECT_ID_COL)
    out = joined.groupby(DOMAIN_ID_COL)[[f"P{i}" for i in (1, 2)]].mean()
    out = out.reset_index().rename(columns={DOMAIN_ID_COL: "domain"})
    out = out[(out["domain"] >= DOMAIN_MIN) & (out["domain"] <= DOMAIN_MAX)]
    return out.sort_values("domain").reset_index(drop=True)


def main():
    rates = domain_rates()
    fig, ax = plt.subplots(figsize=figsize("double", height=3.6))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.235, top=0.96)

    # THE LIMITS GO FIRST. town_bands() skips any span outside the current
    # view and clamps a label to the visible part of its span, so calling it
    # before set_xlim silently dropped Buxton (GIS 7-8): the axes had
    # autoscaled to the shoal spans and 6.5-8.5 fell outside them. Band and
    # label both vanished, with no error (2026-09-17).
    ax.set_xlim(DOMAIN_MIN - 0.5, DOMAIN_MAX + 0.5)

    for name, (lo, hi) in HATTERAS_ANNOTATIONS.shoal_zones.items():
        ax.axvspan(lo - 0.5, hi + 0.5, facecolor=C["ADDED"], alpha=0.13,
                   lw=0, zorder=0.5)
        # a second row, clear of the village names at 0.985: Avon Shoals
        # spans Avon and Wimble Shoals spans Tri-Village, so the two sets
        # of labels overlap in x and must differ in y
        ax.text((lo + hi) / 2, 0.925, name, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=7, color="#8a620e", zorder=7)

    # NAMED, like the shoals and the structures: if a band is worth
    # drawing it is worth naming (Hannah, 2026-09-17). town_bands puts
    # these at the top of the panel, and structures() already knows to
    # tuck its own labels under them.
    town_bands(ax, shade="0.93")
    ax.axhline(0, color=INK_MUTED, lw=0.7, ls=(0, (4, 3)), zorder=2)

    coverage = {}
    x = rates["domain"].values
    for i, ((line_c, fill_c), (y1, y2)) in enumerate(
            zip(PERIOD_COLOURS, DSAS_WINDOWS), 1):
        col = rates[f"P{i}"]
        y = col.values
        # A domain is NaN where no transect in it carries BOTH of the window's
        # shorelines; matplotlib breaks the line and the fill there rather than
        # bridging a gap that has no data behind it.
        ax.fill_between(x, 0, y, color=fill_c, alpha=0.55, lw=0, zorder=3)
        ax.plot(x, y, color=line_c, lw=1.4, zorder=4)
        coverage[i] = (int(col.notna().sum()), len(col))
        print(f"  DSAS {y1}-{y2}: {coverage[i][0]}/{coverage[i][1]} domains, "
              f"island mean {col.mean():+.2f} m/yr")

    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("shoreline change rate\n(m/yr, + seaward)")
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.grid(axis="y", color=GRID_C, lw=0.5, zorder=1)
    open_frame(ax)
    structures(ax)

    # the key names the DSAS window AND the run period it stands in for
    handles = [Line2D([], [], color=c, lw=1.6,
                      label=f"DSAS {y1}–{y2}  (for {rs}–{re})")
               for (c, _), (y1, y2), (rs, re)
               in zip(PERIOD_COLOURS, DSAS_WINDOWS, RUN_PERIODS)]
    handles.append(Patch(facecolor=C["ADDED"], alpha=0.13, label="shoal influence"))
    handles.append(Patch(facecolor="0.93", label="village / community zone"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.005),
               ncol=4, frameon=False, fontsize=8, handlelength=1.8,
               columnspacing=1.6)

    # A CHECK, SO IT LIVES UNDER supporting/. The caption is keyed to the
    # figure's name in the folder's one CAPTIONS.md, beside the primary's.
    out_dir = support_dir(FIG_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "dsas_calibration_periods.png"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(out_dir / "dsas_calibration_periods.pdf", facecolor="white")
    plt.close(fig)

    record_caption(FIG_DIR / "dsas_calibration_periods.png",
        "THE INDEPENDENT CHECK on coastsat_calibration_periods.png, under "
        "supporting/ and drawn in the same style so the two can be compared "
        "directly. Observed shoreline change rate per model domain from the "
        "DSAS transect record, positive seaward. Two differences from the "
        "primary, neither of which can be removed. (1) THE WINDOWS: DSAS has "
        "five digitised shorelines (1978, 1987, 1997, 2009, 2019), so it "
        f"cannot be cut to the run periods; {DSAS_WINDOWS[0][0]}–"
        f"{DSAS_WINDOWS[0][1]} stands in for {RUN_PERIODS[0][0]}–"
        f"{RUN_PERIODS[0][1]} and {DSAS_WINDOWS[1][0]}–{DSAS_WINDOWS[1][1]} "
        f"for {RUN_PERIODS[1][0]}–{RUN_PERIODS[1][1]}, the latter five "
        "years short at the end. (2) THE ESTIMATOR: this is an end-point rate "
        "between the window's two shorelines, where CoastSat is an OLS slope "
        "through every observation in the window. With two shorelines the two "
        "forms coincide, so the comparison is fair as drawn. Bands and "
        "structures are as in the primary figure. COVERAGE is not complete: "
        f"the earlier window resolves {coverage[1][0]} of {coverage[1][1]} "
        f"domains and the later {coverage[2][0]} of {coverage[2][1]}, because "
        "no transect at the cape end carries both of the later window's "
        "shorelines; the curve is broken there rather than interpolated.")
    print(f"Saved: {png}")


if __name__ == "__main__":
    main()
