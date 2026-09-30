"""
coastsat_obx_lrr.py
==============================================================================
Full-record CoastSat LRR for every transect from Cape Point to the Virginia
line: the table handed to the Murray lab (2026-09-27) for their diffusivity
work, plus a README and one figure.

WHY A SIBLING OF coastsat_domain_lrr.py
    That script is driven by transect_domain_lookup.csv, so it only sees the
    transects inside the 90 GIS domains. This hand-off covers 153 km of coast,
    ~90 km of it north of the model domain where no domain exists, so it is driven by the CoastSat
    transect layer instead. The FIT is the same one -- coastsat_lrr.compute_lrr
    on the calendar window, every point, no filter, no weights, >= 3 points --
    so a transect inside the domain gets the number the model target uses.

THE SPEC (decided 2026-09-27)
    window    1984-01-01 -> 2025-12-31 (whole record, last full calendar year)
    extent    usa_NC_0032_0021 (south end of the model domain, GIS 1) through
              usa_NC_0049_0230 (ends at the NC/VA line, 36.550 N)
    columns   the hand-off CSV is what was asked for: transect ID, alongshore
              km, origin lon/lat, LRR, its 95% CI, flag. Every other fit stat
              and CoastSat's beach slope go in supporting/..._full.csv
    flags     flagged, never withheld: fewer than 50 positions, under 20 yr
              between first and last position, within 2 km of Oregon Inlet or
              of the south end at Cape Point (only the two location flags
              fire on this record)

ALONGSHORE DISTANCE
    km from the origin of usa_NC_0032_0021, accumulated transect to transect
    along the coast. Each step between neighbouring origins is projected onto
    the local shore-parallel direction (perpendicular to the two transects'
    mean bearing), because the origins wander cross-shore by up to ~300 m on
    the Currituck Banks and a straight point-to-point sum would add that
    wander as length. Oregon Inlet is the one real gap (~1.1 km).

USAGE
    python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_obx_lrr.py

OUTPUT  data/hatteras_init/5-scr/3-rates/coastsat/lrr/1984_2025_obx/
    coastsat_lrr_obx_1984_2025.csv        the hand-off table (7 columns)
    README.md                             written for the recipients
    lrr_obx_1984_2025.png                 rate vs alongshore km, place names,
                                          1 km median line (a guide only)
    supporting/coastsat_lrr_obx_1984_2025_full.csv   all 19 columns
    supporting/ PDF and CAPTIONS.md
    The maps are drawn from the full table by coastsat_obx_lrr_maps.py; a
    one-file version of the fit for others is coastsat_lrr_standalone.py.
==============================================================================

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-27
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from scipy import stats

_REPO = next(p for p in Path(__file__).resolve().parents
             if (p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)
from coastsat_lrr import load_timeseries, filter_dates, compute_lrr, _empty_lrr  # noqa: E402
from site_layer import hat_observed_rates as obs  # noqa: E402

# ============================================================
# CONFIG
# ============================================================
START_YEAR, END_YEAR = 1984, 2025
START_DATE = "{0}-01-01".format(START_YEAR)
END_DATE = "{0}-12-31".format(END_YEAR)
MIN_OBS = 3                     # below this no fit, as in coastsat_domain_lrr

FIRST_ID = "usa_NC_0032_0021"   # GIS 1, the south end of the model domain
SITES = ["usa_NC_{0:04d}".format(n) for n in range(32, 50)]   # 0032-0049

FLAG_MIN_OBS = 50               # flag: fewer positions than this
FLAG_MIN_SPAN_YR = 20.0         # flag: first-to-last position shorter than this
FLAG_END_KM = 2.0               # flag: this close to Oregon Inlet or Cape Point
INLET_GAP_M = 800.0             # an origin-to-origin step this long is an inlet

UTM = 32618                     # UTM 18N, metres
TRANSECT_LAYER = obs.TRANSECT_DOMAINS / "CoastSat_transect_layer.geojson"
OUT_DIR = obs.COASTSAT_LRR_ROOT / "{0}_{1}_obx".format(START_YEAR, END_YEAR)
TABLE_NAME = "coastsat_lrr_obx_{0}_{1}.csv".format(START_YEAR, END_YEAR)
# The hand-off table is the short one (what was asked for: rate by transect
# ID); every fit and CoastSat field is kept in supporting/ (Hannah, 09-27).
FULL_NAME = "coastsat_lrr_obx_{0}_{1}_full.csv".format(START_YEAR, END_YEAR)
MAIN_COLS = ["transect_id", "alongshore_km", "origin_lon", "origin_lat",
             "lrr_m_yr", "unc_m_yr", "flag"]
FIG_NAME = "lrr_obx_{0}_{1}.png".format(START_YEAR, END_YEAR)
Y_HALF = 8.0                    # m/yr, the lrr tree's fixed axis
V_HALF = 3.0                    # m/yr, colour scale; the maps use the same one

# Places named along the coast, by latitude. The one list for both this
# script's profile and coastsat_obx_lrr_maps.py, which imports it.
PLACES = {
    "Cape Point": 35.231,          # the south-end transect, just north of the tip
    "Buxton": 35.267, "Avon": 35.352, "Salvo": 35.542, "Waves": 35.566, "Rodanthe": 35.594,
    "Oregon Inlet": 35.770, "Nags Head": 35.957, "Kill Devil Hills": 36.030,
    "Kitty Hawk": 36.064, "Southern Shores": 36.130, "Duck": 36.163,
    "Corolla": 36.377, "Carova": 36.520,
}


# ============================================================
# TRANSECTS
# ============================================================

def load_transects():
    """The CoastSat layer's transects for SITES, south to north, from FIRST_ID.

    Order is site then transect number; both run south to north here (checked
    2026-09-27: four steps go <= 31 m south, all local wiggles)."""
    g = gpd.read_file(TRANSECT_LAYER)
    g = g[g["site_id"].isin(SITES)].copy()
    g["transect_id"] = g["id"].str.replace("-", "_")
    # CoastSat's own per-transect beach slope, the one its tidal correction
    # used, carried through as given (added 2026-09-27 at the recipients' use).
    g = g.rename(columns={"cil": "beach_slope_ci_lower", "ciu": "beach_slope_ci_upper"})
    g["num"] = g["transect_id"].str[-4:].astype(int)
    g = g.sort_values(["site_id", "num"]).reset_index(drop=True)
    start = g.index[g["transect_id"] == FIRST_ID][0]
    g = g.loc[start:].reset_index(drop=True)

    ll = np.array([ls.coords[0] + ls.coords[-1] for ls in g.geometry])
    g["origin_lon"], g["origin_lat"] = ll[:, 0], ll[:, 1]
    g["seaward_lon"], g["seaward_lat"] = ll[:, 2], ll[:, 3]

    u = g.to_crs(UTM).geometry
    o = np.array([ls.coords[0] for ls in u])
    e = np.array([ls.coords[-1] for ls in u])
    shore_normal = (e - o) / np.hypot(*(e - o).T)[:, None]
    along = np.column_stack([-shore_normal[:, 1], shore_normal[:, 0]])  # normal rotated +90: northward on an east-facing coast
    step = np.diff(o, axis=0)
    unit = along[:-1] + along[1:]
    unit /= np.hypot(*unit.T)[:, None]
    ds = np.einsum("ij,ij->i", step, unit)
    g["alongshore_km"] = np.round(np.concatenate([[0.0], np.cumsum(ds)]) / 1000.0, 3)
    g["_gap_m"] = np.concatenate([[0.0], np.hypot(*step.T)])
    return g


# ============================================================
# FIT AND FLAGS
# ============================================================

def fit_all(g):
    csv_root = obs.COASTSAT_TIMESERIES
    rows = []
    for tid, site in zip(g["transect_id"], g["site_id"]):
        path = csv_root / "{0}_timeseries".format(site) / "{0}.csv".format(tid)
        if not path.is_file():
            rows.append(_empty_lrr(0))
            continue
        df = filter_dates(load_timeseries(str(path)), START_DATE, END_DATE)
        if len(df) < MIN_OBS:
            rows.append(_empty_lrr(len(df)))
            continue
        r = compute_lrr(df)
        # compute_lrr rounds p to 6 dp, which printed 2,436 of them as 0.0.
        # Same x and y as its fit (years since the first position), unrounded.
        yrs = (df["date"] - df["date"].min()).dt.total_seconds() / 86400.0 / 365.25
        r["p_value"] = float(stats.linregress(yrs.values, df["chainage_m"].values).pvalue)
        rows.append(r)
    fit = pd.DataFrame(rows)
    fit["span_yr"] = ((pd.to_datetime(fit["end_date"]) - pd.to_datetime(fit["start_date"]))
                      .dt.days / 365.25).round(2)
    return fit


def flags(t):
    """Semicolon-joined reasons, empty when none. Nothing is withheld."""
    gaps = t.index[t["_gap_m"] > INLET_GAP_M]
    if len(gaps) != 1:
        raise RuntimeError("expected one inlet gap (Oregon Inlet), found {0}".format(len(gaps)))
    i = gaps[0]
    inlet_km = 0.5 * (t.at[i - 1, "alongshore_km"] + t.at[i, "alongshore_km"])
    out = []
    for _, r in t.iterrows():
        f = []
        if not np.isfinite(r["lrr_m_yr"]):
            f.append("no_fit")
        if r["n_obs"] < FLAG_MIN_OBS:
            f.append("fewer_than_{0}_positions".format(FLAG_MIN_OBS))
        if not (r["span_yr"] >= FLAG_MIN_SPAN_YR):
            f.append("span_under_{0:g}_yr".format(FLAG_MIN_SPAN_YR))
        if abs(r["alongshore_km"] - inlet_km) <= FLAG_END_KM:
            f.append("within_{0:g}_km_of_oregon_inlet".format(FLAG_END_KM))
        if r["alongshore_km"] <= FLAG_END_KM:
            f.append("within_{0:g}_km_of_cape_point".format(FLAG_END_KM))
        out.append(";".join(f))
    return out, inlet_km


# ============================================================
# FIGURE
# ============================================================

def draw(t, inlet_km, path):
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.lines import Line2D
    from site_layer.hat_figure_style import (
        apply_style, figsize, open_frame, save, caption, mark_offaxis,
        INK, INK_MUTED, GRID_C)
    apply_style()

    # Flagged transects are drawn like the rest (Hannah, 2026-09-27: the hollow
    # markers came off the figures); the flag lives in the table.
    x = t["alongshore_km"].to_numpy()
    y = t["lrr_m_yr"].to_numpy()
    # Running median over +/-0.5 km of coast, each side of the inlet on its
    # own. A count window (21 transects) was tried first: it spans well over
    # 1 km where transects are sparse and reached across the inlet.
    side = x > inlet_km
    med = np.full_like(y, np.nan)
    for i in range(len(x)):
        w = (np.abs(x - x[i]) <= 0.5) & (side == side[i]) & np.isfinite(y)
        if w.sum() >= 5:
            med[i] = np.median(y[w])

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.52), constrained_layout=True)
    ax.axhline(0, color=INK, lw=0.6, zorder=2)
    # Place names along the top edge, at the transect nearest each latitude;
    # Oregon Inlet at the gap itself, and the two ends of the reach.
    marks = {"Cape Point": 0.0}
    for name, lat in PLACES.items():
        marks[name] = (inlet_km if name == "Oregon Inlet" else
                       float(t["alongshore_km"].iloc[(t["origin_lat"] - lat).abs().argmin()]))
    marks["NC/VA line"] = float(x.max())
    # every named place gets the same faint dashed line (Oregon Inlet was the
    # only one until 09-27); the ends of the reach are the axis limits
    for name, km in marks.items():
        if name not in ("Cape Point", "NC/VA line"):
            ax.axvline(km, color=INK_MUTED, lw=0.5, ls=(0, (3, 2)), alpha=0.6, zorder=1)
    top = ax.secondary_xaxis("top")
    top.set_xticks(list(marks.values()), labels=list(marks.keys()))
    # vertical: at 40 degrees the pairs 4 km apart (Cape Point/Buxton,
    # Carova/NC-VA line) ran into each other
    top.tick_params(axis="x", labelsize=7, length=3, width=0.5, colors=INK_MUTED,
                    labelcolor=INK, labelrotation=90)
    top.spines["top"].set_visible(False)
    # Coloured by rate on the maps' scale (Hannah, 2026-09-27: erosion red,
    # accretion blue, deeper with severity); saturates at +/-V_HALF.
    cmap = plt.get_cmap("RdBu")
    norm = mcolors.Normalize(-V_HALF, V_HALF)
    # thin grey edge so rates near zero (near-white) stay visible (09-27)
    ax.scatter(x, y, s=6, c=np.nan_to_num(y), cmap=cmap, norm=norm,
               edgecolors="0.72", linewidths=0.12, zorder=3)   # lighter than the maps': 2,984 overlap here
    for s_ in (~side, side):
        ax.plot(x[s_], med[s_], color=INK, lw=0.9, zorder=4)
    ax.set_xlim(0, np.ceil(x.max() / 10) * 10)
    ax.set_ylim(-Y_HALF, Y_HALF)
    ax.set_yticks(np.arange(-Y_HALF, Y_HALF + 1, 2))
    off = mark_offaxis(ax, x, y, Y_HALF, color=cmap(0.0))
    ax.set_xlabel("Alongshore distance from Cape Point (km, south → north)")
    ax.set_ylabel("Shoreline change rate, LRR (m/yr)")
    ax.grid(True, axis="y", color=GRID_C, lw=0.4)
    open_frame(ax)
    ax.legend(handles=[Line2D([], [], color=INK, lw=0.9, label="1 km median")],
              loc="lower right", frameon=False)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                      orientation="horizontal", location="bottom", extend="both",
                      shrink=0.45, aspect=35, pad=0.02)
    cb.set_label("Transect rate, LRR (m/yr)")
    cb.outline.set_linewidth(0.5)

    clause = ""
    if off:
        ks = [k for k, _ in off]
        clause = (" {n} transects just north of Oregon Inlet ({a:.1f}-{b:.1f} km) fall "
                  "beyond -{h:g} m/yr, down to {m:+.1f} m/yr; they are marked with a "
                  "triangle at the axis edge.").format(
                      n=len(off), a=min(ks), b=max(ks), h=Y_HALF, m=min(v for _, v in off))
    caption(fig, (
        "CoastSat shoreline change rate, Cape Point to the North Carolina / Virginia "
        "line, {s}-{e}. Each point is one CoastSat transect ({n} transects, sites "
        "usa_NC_0032-0049, from {first}): the ordinary least-squares slope of every "
        "tidally corrected shoreline position between 1 January {s} and 31 December "
        "{e}, with no outlier filter and no weighting; positive is seaward "
        "(accretion). Distance is measured along the coast from the origin of "
        "{first} at Cape Point. Points are coloured by rate, red for erosion and blue "
        "for accretion, on the same scale as the maps, saturating at ±{v:g} m/yr. "
        "Every transect is drawn, including the ones the table flags. The black line is the median of the transect rates within 0.5 km either "
        "side (a 1 km window, about 21 transects), computed separately north and south "
        "of Oregon Inlet; it is left blank where fewer than 5 transects fall in the "
        "window, which happens only at about 144-148 km, where CoastSat's transects "
        "are 300-570 m apart. Place names along the top edge, each with a dashed "
        "line, mark the transect nearest the town's latitude; Oregon Inlet is marked "
        "at the gap in the transects.{clause}").format(
            s=START_YEAR, e=END_YEAR, n=len(t), first=FIRST_ID, v=V_HALF, clause=clause))
    save(fig, path, close=True, dpi=300)
    return off


# ============================================================
# README
# ============================================================

def write_readme(t, inlet_km, off, path):
    v = t["lrr_m_yr"]
    n_fit = int(v.notna().sum())
    counts = pd.Series([f for s in t["flag"] for f in s.split(";") if f]).value_counts()
    for name in ("fewer_than_{0}_positions".format(FLAG_MIN_OBS),
                 "span_under_{0:g}_yr".format(FLAG_MIN_SPAN_YR),
                 "within_{0:g}_km_of_oregon_inlet".format(FLAG_END_KM),
                 "within_{0:g}_km_of_cape_point".format(FLAG_END_KM), "no_fit"):
        counts[name] = int(counts.get(name, 0))
    n_flag = int((t["flag"] != "").sum())
    north = t.iloc[-1]
    south = t[t["alongshore_km"] < inlet_km]
    nth = t[t["alongshore_km"] > inlet_km]
    last_obs = t["end_date"].dropna().max()
    first_obs = t["start_date"].dropna().min()
    flag_lines = "\n".join("| `{0}` | {1} |".format(k, int(n)) for k, n in counts.items())

    text = f"""# CoastSat shoreline change rates, Cape Point to the Virginia line, {START_YEAR}-{END_YEAR}

One rate per CoastSat transect over the whole satellite record, for the
{len(t):,} transects from the south end of the Hatteras model domain
(`{FIRST_ID}`, Cape Point) north to the North Carolina / Virginia line
(`{north['transect_id']}`, {north['origin_lat']:.3f} N). Prepared 2026-09-27
for the Murray lab.

**File:** `{TABLE_NAME}`, one row per transect, south to north: the
transect ID, where it is, its rate and the rate's uncertainty. Every other
fit statistic and CoastSat field is in `supporting/{FULL_NAME}`.

## How the rate is measured

- **Data:** CoastSat shoreline time series (Vos et al.), downloaded per
  transect from the CoastSat portal (coastsat.space), sites `usa_NC_0032` to
  `usa_NC_0049`: one cross-shore position (chainage, metres along the
  transect) per usable satellite image. CoastSat tidally corrects these
  positions itself, with the FES2022 tide model and a satellite-derived beach
  slope for each transect. The positions used run {first_obs} to
  {last_obs}; the downloaded files continue into January 2026, so they come
  from the continuously updated portal and extend past the archived US East
  Coast release (Zenodo v1.0, 9 June 2025,
  doi:10.5281/zenodo.15626280, CC-BY-4.0, whose description gives the
  processing above). Please cite CoastSat if you use these numbers.
- **Window:** 1 January {START_YEAR} to 31 December {END_YEAR}, the whole
  record through the last full calendar year. The first images are from 1984;
  the few January 2026 positions are left out.
- **Rate:** linear regression rate (LRR), the ordinary least-squares slope of
  shoreline position against time, in m/yr. Every position in the window is
  used: **we apply no outlier filter and no weighting** on top of CoastSat's
  own processing. A transect needs at least
  {MIN_OBS} positions to get a rate.
- **Sign:** CoastSat chainage increases seaward, so **positive = seaward
  movement (accretion), negative = landward (erosion).**
- This is the same fit, on the same software, that produces the rate the
  Hatteras CASCADE model is scored against, so inside the model domain these
  numbers are directly comparable to that work (over a different window).

Things worth knowing about the record:

- Sampling is uneven in time. Averaged over these transects: about 8
  positions a year in 1984-1998, 18 a year in 1999-2020, and 29 a year in
  2021-2025. An unweighted OLS fit therefore leans toward recent years.
- `unc_m_yr` and `p_value` assume independent residuals. Shoreline position
  has seasonal and storm-driven memory, so consecutive positions are not
  independent and the true uncertainty is wider than stated; treat
  `unc_m_yr` as a lower bound and `p_value` as optimistic.
- Along Hatteras Island there is a real, abrupt seaward shift of about 17 m
  in 2021 (confirmed against the dune line, not a nourishment artefact). Over
  a {END_YEAR - START_YEAR + 1}-year fit it moves the rate much less than it moves
  a 10-15 year window; if you cut the record into sub-periods, it matters.

## Columns

`{TABLE_NAME}`:

| column | meaning |
|---|---|
| `transect_id` | CoastSat transect ID, `usa_NC_<site>_<transect>` |
| `alongshore_km` | distance north along the coast from the origin of `{FIRST_ID}` (km) |
| `origin_lon`, `origin_lat` | landward end of the transect (WGS84, decimal degrees) |
| `lrr_m_yr` | shoreline change rate (m/yr; positive = seaward/accretion, negative = landward/erosion) |
| `unc_m_yr` | 95 % confidence half-width on the rate (m/yr); see the caveat above |
| `flag` | reasons to look twice, `;`-separated; empty when none (see Flags) |

`supporting/{FULL_NAME}` has those columns plus:

| column | meaning |
|---|---|
| `site_id` | CoastSat site |
| `seaward_lon`, `seaward_lat` | seaward end of the transect (WGS84) |
| `r_squared` | R² of the fit (dimensionless) |
| `p_value` | p-value of the slope (two-sided, null: zero trend), 3 significant figures |
| `n_obs` | positions used in the fit |
| `start_date`, `end_date` | first and last position used (`YYYY-MM-DD`, UTC) |
| `span_yr` | years between them |
| `beach_slope` | CoastSat's satellite-derived beach-face slope (tan β, dimensionless), the value CoastSat used to tidally correct this transect |
| `beach_slope_ci_lower`, `beach_slope_ci_upper` | the confidence interval CoastSat reports on that slope (its `cil` / `ciu`) |

The three beach-slope columns are copied unchanged from the CoastSat
transect layer; nothing here estimates them. CoastSat reports the slope in
discrete steps (median 0.07, range 0.015-0.18 over these transects), and
neighbouring transects often share a value.

## Flags

Flagged transects are **kept, with their rate**; the flag says why a reader
might treat them with care. {n_flag} of {len(t):,} transects carry one or more.

| flag | transects |
|---|---|
{flag_lines}

- `fewer_than_{FLAG_MIN_OBS}_positions`: too few positions for a stable slope.
- `span_under_{FLAG_MIN_SPAN_YR:g}_yr`: the positions cover less than
  {FLAG_MIN_SPAN_YR:g} years of the window, so the rate is not the full-record rate.
- `within_{FLAG_END_KM:g}_km_of_oregon_inlet`: inlet-flank shoreline (Bodie Island
  spit and the north end of Pea Island), where inlet migration and the
  terminal groin, not open-coast processes, set the rate.
- `within_{FLAG_END_KM:g}_km_of_cape_point`: the first {FLAG_END_KM:g} km north of the
  tip of Cape Hatteras, where the shoreline turns and the cape shoals shelter it.
- `no_fit`: fewer than {MIN_OBS} positions; no rate.

## Alongshore distance

Each step between neighbouring transect origins is projected onto the local
shore-parallel direction (perpendicular to the transects' bearing) and
summed. That keeps the landward origins' cross-shore wander (up to ~300 m on
the Currituck Banks) from being counted as length. Oregon Inlet sits at
{inlet_km:.1f} km, the one real gap in the transects (~1.1 km). The transect
spacing is nominally 50 m; around `usa_NC_0049_0113`-`0116`, in the
northernmost site, CoastSat's transects are 300-570 m apart.

## Summary

| | south of Oregon Inlet | north of Oregon Inlet | all |
|---|---|---|---|
| transects | {len(south)} | {len(nth)} | {len(t)} |
| km | 0-{south['alongshore_km'].max():.1f} | {nth['alongshore_km'].min():.1f}-{nth['alongshore_km'].max():.1f} | {t['alongshore_km'].max():.1f} |
| median LRR (m/yr) | {south['lrr_m_yr'].median():+.2f} | {nth['lrr_m_yr'].median():+.2f} | {v.median():+.2f} |
| eroding (%) | {100 * (south['lrr_m_yr'] < 0).mean():.0f} | {100 * (nth['lrr_m_yr'] < 0).mean():.0f} | {100 * (v < 0).mean():.0f} |

{n_fit:,} of {len(t):,} transects have a rate.

## Figures

- `{FIG_NAME}`: every transect's rate against alongshore distance,
  coloured red (erosion) to blue (accretion), with place names along the top.
  The black line ("1 km median") is drawn over the points as a guide; it
  is not applied to the data. At each transect it is the median rate of all
  transects within 0.5 km either side (about 21), computed separately on each
  side of Oregon Inlet. It is left blank where fewer than 5 transects fall in
  that window, which happens only at about 144-148 km (see Alongshore
  distance). Twelve transects just north of Oregon Inlet are faster than
  -8 m/yr and sit off the axis, marked with a triangle; their values are in the
  CSV.
- `lrr_obx_{START_YEAR}_{END_YEAR}_map_overview.png`: map of the whole reach,
  transects coloured by rate, beside the rates plotted against northing;
  boxes A-D mark the regional maps.
- `..._map_A_cape_point_to_salvo.png`, `..._map_B_rodanthe_to_south_nags_head.png`,
  `..._map_C_nags_head_to_duck.png`, `..._map_D_corolla_to_virginia.png`:
  the same two panels zoomed to ~35-40 km each (1 km overlap). The colour
  scale is fixed at +/-3 m/yr in every map, so colours compare across them.
- The figures draw every transect the same way; the flags are in the table
  only.

Captions in `supporting/CAPTIONS.md`, PDFs in `supporting/`. The maps are
drawn by `coastsat_obx_lrr_maps.py` (basemap: Esri World Shaded Relief,
which carries no labels, shown in greyscale; every name on a map is one
placed here). Latitude is marked on each map's left edge; it is exact at
the coast and the (b) panels share that axis.

## Reproduce

```
python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_obx_lrr.py        # table, README, profile figure
python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_obx_lrr_maps.py   # the five maps (needs internet for the basemap)
```

**To run the same analysis on your own CoastSat data**, use the one-file
version, `coastsat_lrr_standalone.py` (in the repository at
`scripts/input_prep/5-scr/3-rates/coastsat/lrr/`, and sent alongside this
folder; needs only pandas and scipy):

```
python coastsat_lrr_standalone.py <folder of CoastSat CSVs> 1984 2025 rates.csv
```

It writes `transect_id, lrr_m_yr, unc_m_yr, n_obs` for every transect it
finds, and reproduces those columns of this table exactly.

in the CASCADE repository (branch `hannahaline/hatteras-cascade`). The fit
is `coastsat_lrr.compute_lrr` (`scripts/input_prep/5-scr/lib/`); the
transect geometry is the CoastSat global transect layer.
"""
    path.write_text(text, encoding="utf-8")


# ============================================================
# MAIN
# ============================================================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    g = load_transects()
    print("{0} transects, {1} -> {2}".format(len(g), g["transect_id"].iloc[0],
                                           g["transect_id"].iloc[-1]))
    fit = fit_all(g)
    t = pd.concat([g.drop(columns="geometry").reset_index(drop=True), fit], axis=1)
    t["flag"], inlet_km = flags(t)

    cols = ["transect_id", "site_id", "alongshore_km",
            "origin_lon", "origin_lat", "seaward_lon", "seaward_lat",
            "lrr_m_yr", "unc_m_yr", "r_squared", "p_value", "n_obs",
            "start_date", "end_date", "span_yr",
            "beach_slope", "beach_slope_ci_lower", "beach_slope_ci_upper", "flag"]
    for c in ["origin_lon", "origin_lat", "seaward_lon", "seaward_lat"]:
        t[c] = t[c].round(6)
    out = t[cols].copy()
    out["p_value"] = out["p_value"].map(lambda v: "" if pd.isna(v) else "{0:.3g}".format(v))
    (OUT_DIR / "supporting").mkdir(exist_ok=True)
    out.to_csv(OUT_DIR / "supporting" / FULL_NAME, index=False)
    out[MAIN_COLS].to_csv(OUT_DIR / TABLE_NAME, index=False)
    print("Saved", OUT_DIR / TABLE_NAME, "and supporting/" + FULL_NAME)

    off = draw(t, inlet_km, OUT_DIR / FIG_NAME)
    write_readme(t, inlet_km, off, OUT_DIR / "README.md")
    print("Saved", OUT_DIR / FIG_NAME, "and README.md")
    print("fitted {0}/{1}; flagged {2}; Oregon Inlet at {3:.2f} km; off-axis {4}".format(
        int(t["lrr_m_yr"].notna().sum()), len(t), int((t["flag"] != "").sum()),
        inlet_km, off))
    return t


if __name__ == "__main__":
    main()
