"""
Does each fill footprint show up in CoastSat? Shoreline position the year before placement against the year after.

    python scripts/input_prep/4-mgmt-forcings/nourishment_extent_coastsat.py

Per transect, the median CoastSat position in the 12 months before a project's
first placement date and the 12 months after its last; the change is set against
the transects outside every fill footprint. Dates come from
nourishment/1-sources/nourishment_placement_dates.csv, footprints from
HATTERAS_NOURISHMENT_PROJECTS. A check on the forcing; nothing here feeds a run.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-04
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "input_prep" / "5-scr" / "lib"))

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

import coastsat_lrr as cl  # noqa: E402
from site_layer import hat_observed_rates as obs  # noqa: E402
from site_layer.hat_figure_style import (C, DOMAIN_AXIS_LABEL, INK, INK_MUTED, _title, apply_style,  # noqa: E402
                                         figsize, open_frame, record_caption, save, structures, town_bands)
from site_layer.hat_topo_version import MGMT_ROOT  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_NOURISHMENT_PROJECTS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
DATES_CSV = MGMT_ROOT / "nourishment" / "1-sources" / "nourishment_placement_dates.csv"
OUT_DIR = MGMT_ROOT / "nourishment" / "4-extent-checks" / "coastsat"
WINDOW_DAYS = 365        # length of the before and after windows; a full year so the seasons cancel
MIN_OBS = 5              # passes a transect needs in each window
PAD_DOMAINS = 10         # domains drawn either side of a footprint
RUN_TRANSECTS = 5        # running-median width for the drawn line and the extent
N_MAD = 2.0              # width of the drawn background band, in robust SDs
HALF = 0.5               # the extent is where the running median stays above this share of its peak over background
# -----------------------------------------------------------------------------


# Every transect with its domain and a fractional alongshore position, south to north
def transect_frame():
    lk = pd.read_csv(obs.TRANSECT_DOMAINS / "transect_domain_lookup.csv").dropna(subset=["domain_number"])
    lk["domain_number"] = lk.domain_number.astype(int)
    lk["k"] = lk.groupby("domain_number").cumcount()
    n = lk.groupby("domain_number").transect_id.transform("size")
    lk["x"] = lk.domain_number - 0.5 + (lk.k + 0.5) / n
    return lk.set_index("transect_id")


def timeseries_path(tid):
    site = tid.rsplit("_", 1)[0]
    return obs.COASTSAT_TIMESERIES / f"{site}_timeseries" / f"{tid}.csv"


# Median position in a window, NaN with too few passes
def window_median(df, lo, hi):
    w = df[(df.date >= lo) & (df.date < hi)]
    return float(w.chainage_m.median()) if len(w) >= MIN_OBS else np.nan


def project_change(p, dates, frame):
    start = pd.Timestamp(dates.start_date, tz="UTC")
    end = pd.Timestamp(dates.end_date, tz="UTC") + pd.Timedelta(days=1)
    span = pd.Timedelta(days=WINDOW_DAYS)
    lo, hi = min(p.gis_domains) - PAD_DOMAINS, max(p.gis_domains) + PAD_DOMAINS
    sub = frame[frame.domain_number.between(lo, hi)].copy()
    rows = []
    for tid in sub.index:
        f = timeseries_path(tid)
        if not f.exists():
            continue
        df = cl.load_timeseries(str(f))
        b, a = window_median(df, start - span, start), window_median(df, end, end + span)
        rows.append(dict(transect_id=tid, before_m=b, after_m=a))
    out = sub.join(pd.DataFrame(rows).set_index("transect_id"), how="inner")
    out["change_m"] = out.after_m - out.before_m
    return out.sort_values("x"), (start - span, start, end, end + span)


# Contiguous run around the peak inside the model footprint where the running median stays above
# background + HALF x (peak - background): the half-maximum width of the fill
def observed_extent(t, background, first, last):
    sm = t.change_m.rolling(RUN_TRANSECTS, center=True, min_periods=3).median()
    in_fp = t.domain_number.between(first, last).to_numpy()
    peak_pos = np.where(in_fp, sm.to_numpy(), np.nan)
    if np.all(np.isnan(peak_pos)):
        return None, sm, np.nan
    i = int(np.nanargmax(peak_pos))
    level = background + HALF * (sm.iloc[i] - background)
    above = (sm > level).to_numpy()
    if sm.iloc[i] <= background:
        return None, sm, level
    j0 = i
    while j0 > 0 and above[j0 - 1]:
        j0 -= 1
    j1 = i
    while j1 < len(above) - 1 and above[j1 + 1]:
        j1 += 1
    return (t.iloc[j0], t.iloc[j1]), sm, level


def main():
    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    frame = transect_frame()
    dates = pd.read_csv(DATES_CSV).set_index("project")
    filled_domains = {g for q in HATTERAS_NOURISHMENT_PROJECTS for g in q.gis_domains}
    projects = sorted(HATTERAS_NOURISHMENT_PROJECTS, key=lambda q: (q.year, min(q.gis_domains)))

    fig, axes = plt.subplots(len(projects), 1, figsize=figsize("double", height=2.3 * len(projects)))
    fig.subplots_adjust(left=0.08, right=0.99, top=0.96, bottom=0.11, hspace=0.45)
    summary, per_transect = [], []
    for i, (ax, p) in enumerate(zip(axes, projects)):
        d = dates.loc[p.name]
        t, (b0, b1, a0, a1) = project_change(p, d, frame)
        first, last = min(p.gis_domains), max(p.gis_domains)
        ctrl = t[~t.domain_number.isin(filled_domains)].change_m.dropna()
        med = float(ctrl.median())
        sd = float(1.4826 * (ctrl - med).abs().median())
        thr = med + N_MAD * sd
        ext, sm, level = observed_extent(t, med, first, last)
        inside = t[t.domain_number.between(first, last)].change_m.dropna()
        per_transect.append(t.assign(project=p.name, model_year=p.year, smoothed_change_m=sm)
                            [["project", "model_year", "domain_number", "x", "before_m", "after_m",
                              "change_m", "smoothed_change_m"]])

        ax.axvspan(first - 0.5, last + 0.5, color=C["ACCENT_FILL"], alpha=0.35, lw=0, zorder=0)
        ax.axhspan(med - N_MAD * sd, thr, color=C["BASE_FILL"], alpha=0.6, lw=0, zorder=1)
        ax.axhline(med, color=C["BASE"], lw=0.8, ls=(0, (3, 2)), zorder=2)
        ax.axhline(0, color=INK_MUTED, lw=0.5, zorder=2)
        ax.scatter(t.x, t.change_m, s=4, color=C["BASE"], lw=0, zorder=3)
        ax.plot(t.x, sm, color=INK, lw=1.4, zorder=4)
        if ext is not None:
            ax.hlines(level, ext[0].x, ext[1].x, color=C["ACCENT"], lw=0.8, ls=(0, (1, 1.5)), zorder=5)
            ax.annotate("", xy=(ext[0].x, 0.93), xytext=(ext[1].x, 0.93),
                        xycoords=("data", "axes fraction"),
                        arrowprops=dict(arrowstyle="|-|", color=C["ACCENT"], lw=1.4, mutation_scale=3))
            ax.text((ext[0].x + ext[1].x) / 2, 0.96, f"CoastSat GIS {ext[0].domain_number}–{ext[1].domain_number}",
                    transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=8, color=C["ACCENT"])
        ax.set_xlim(t.x.min() - 0.5, t.x.max() + 0.5)
        ax.set_ylabel("After − before (m)")
        ax.grid(axis="y")
        open_frame(ax)
        _title(ax, i, f"{p.name} {p.year}: model GIS {first}–{last}, placed "
                      f"{d.start_date} to {d.end_date}")
        structures(ax, label=i == 0)
        summary.append(dict(
            project=p.name, model_year=p.year, model_first_gis=first, model_last_gis=last,
            placement_start=d.start_date, placement_end=d.end_date,
            before_window=f"{b0.date()} to {(b1 - pd.Timedelta(days=1)).date()}",
            after_window=f"{a0.date()} to {(a1 - pd.Timedelta(days=1)).date()}",
            observed_first_gis=None if ext is None else int(ext[0].domain_number),
            observed_last_gis=None if ext is None else int(ext[1].domain_number),
            median_change_in_footprint_m=round(float(inside.median()), 1),
            background_median_m=round(med, 1), background_robust_sd_m=round(sd, 1),
            half_peak_level_m=round(float(level), 1), background_band_top_m=round(thr, 1), n_transects=int(t.change_m.notna().sum()),
            n_background=int(len(ctrl))))
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    handles = [Line2D([], [], color=C["BASE"], marker="o", ms=3, ls="none", label="Transect"),
               Line2D([], [], color=INK, lw=1.4, label=f"Running median, {RUN_TRANSECTS} transects"),
               Line2D([], [], color=C["BASE"], lw=0.8, ls=(0, (3, 2)), label="Background median"),
               Patch(color=C["BASE_FILL"], label=f"Background ±{N_MAD:g} robust SD"),
               Patch(color=C["ACCENT_FILL"], alpha=0.6, label="Model footprint"),
               Line2D([], [], color=C["ACCENT"], lw=1.4, label="CoastSat extent (half-peak width)")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    out = save(fig, OUT_DIR / "nourishment_extent_coastsat", close=True)[0]

    s = pd.DataFrame(summary)
    s.to_csv(OUT_DIR / "nourishment_extent_coastsat_summary.csv", index=False)
    pd.concat(per_transect).to_csv(OUT_DIR / "nourishment_extent_coastsat_transects.csv", index=False,
                                   float_format="%.2f")
    record_caption(out, (
        "Shoreline position change across each fill the hindcast applies, from CoastSat. For every transect, "
        f"the median position over the {WINDOW_DAYS} days before the first placement date minus the median "
        f"over the {WINDOW_DAYS} days after the last (at least {MIN_OBS} passes in each), positive seaward. "
        f"Placement dates are from nourishment_placement_dates.csv. The black line is a {RUN_TRANSECTS}-transect "
        "running median. Background is every transect in the panel outside all fill footprints: its median "
        f"(dashed) and ±{N_MAD:g} robust SD (1.4826 × MAD, grey band). The CoastSat extent is the half-peak width: "
        "the contiguous run around the running median's peak inside the model footprint where it stays above "
        "the background median plus half of (peak − background median), drawn as the dotted level. Shaded columns are the model footprints "
        "(hatteras_site_config.HATTERAS_NOURISHMENT_PROJECTS). Buxton and Avon 2022 were placed at the same "
        "time, so each sits in the other's panel window and both are left out of the background."))
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(s.to_string(index=False))
    print(out)


if __name__ == "__main__":
    main()
