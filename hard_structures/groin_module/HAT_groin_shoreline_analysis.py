"""
HAT_groin_shoreline_analysis_v2.py
================================
Buxton Groin Field effectiveness analysis using multi-source shoreline
observations, all measured against CoastSat's own transect network.

ARCHITECTURE NOTE: wet-dry lines, NC Coastal Management shorelines, and
CoastSat shoreline positions are all extracted onto CoastSat's own
transect network (COASTSAT_TRANSECT_GEOM) -- CoastSat chainage is read
directly from its per-transect CSVs (load_coastsat_chainage()), and
wet-dry/NC state chainage comes from geometric intersection with the
same transects (extract_chainage_by_intersection()). Alongshore
position is built from a trusted along-coast ID order (see
_compute_alongshore_positions()), and CASCADE domain assignment comes
from the authoritative HAT_domains.json reference (see
load_domain_reference() / assign_domain_from_northing()), not a
formula.

Analyses (each per-transect, aggregated alongshore):

    1. PRE-INSTALLATION LRR (all shorelines dated before
       PRE_INSTALLATION_YEAR_CUTOFF, currently 1970)
       Establishes background shoreline change rate before groin
       installation, using every historical/aerial shoreline available
       pre-1970 rather than a fixed hand-picked list of years, pooled
       across every source that has pre-1970 data (a per-source
       breakdown is printed). The 1966 Buxton nourishment does bias the
       1966/1967 shorelines at Buxton-area transects -- those
       observations are KEPT in the regression (not dropped) and only
       FLAGGED (nourishment_zone=True), so a reader can see exactly
       which points are nourishment-influenced rather than having them
       silently removed. See NOURISHMENT_EXCLUSIONS.

    2. FULL-PERIOD REGRESSION with SIMPLE PIECEWISE CHANGE-POINT DETECTION
       Fits a single-slope regression over EVERY observation a
       transect has (no year cutoff), and separately searches
       BREAKPOINT_SEARCH_START..END for the best two-segment
       (piecewise linear) breakpoint year. This is a data-driven check
       independent of the assumed era boundaries below: does the
       record itself suggest a regime change, and if so, roughly when.
       (The RSS-reduction fraction reported is descriptive, not a
       significance test -- an exhaustive breakpoint search will
       always find *some* "best" breakpoint even under a true
       single-slope null, so a large reduction is suggestive, not
       confirmatory, without a permutation/Davies-type test.)

    3. PER-ERA LRR
       One rate per transect per structural era. Three eras (see
       ERAS below and HAT_groin_zone_investigation.py for the full
       cited maintenance history): Pre-install (no structure),
       Functional groin (1970-1995: built, damaged/repaired during
       construction, 1975 steel-pile repair, 1980/82 anti-flanking
       extensions, 1994 Gordon damage, 1995 south-groin repair -- the
       LAST documented repair), and Deteriorated (1996-2024: no
       further repairs after 1995, including through Hurricane Isabel
       in 2003, which is annotated on the decadal plots as a notable
       event but is not an era boundary since no maintenance bracket
       it). Shown as per-transect scatter only -- no binning, no
       smoothing, no connecting line, for any era.

    4. DECADAL LRR + ANOMALY (LRR minus pre-install baseline)
       At each transect, fits LRR within fixed, NON-OVERLAPPING
       DECADE_LENGTH_YEARS-year bins starting at DECADE_START_YEAR
       (1960s, 1970s, ...) -- replaces an earlier 10-year window that
       slid one year at a time (adjacent "center years" shared ~90% of
       their underlying observations, which read as smoother/more
       continuous than the survey frequency actually supports). Each
       decade here is a distinct, non-overlapping slice of time.
       The regional baseline subtracted to get the anomaly connects the
       sorted per-transect pre-install points directly (linear
       interpolation, no binning/smoothing) -- see
       interp_preinstall_baseline(). The decadal/anomaly/distance-band/
       signal-extent CSVs are all still written even though their
       companion plots were dropped as unhelpful (see Outputs below);
       the numbers may still be useful even without a chart.

    5. SIGNAL EXTENT OVER TIME (computed, not plotted -- see Outputs)
       The CONTIGUOUS alongshore distance, measured outward from the
       groin, over which the anomaly exceeds SIGNAL_ANOMALY_THRESHOLD_M_YR
       -- i.e. "how far" the groin's influence reaches on each side,
       tracked by decade. The threshold is checked against an empirical
       noise floor (far-field pre-install rate standard deviation,
       printed by compute_preinstall_lrr) rather than assumed outright
       -- see the config comment above SIGNAL_ANOMALY_THRESHOLD_M_YR.

    6. SHORELINE EVOLUTION GIF (groin area)
       Animated view of RAW shoreline position through time, zoomed to
       domain 1 to PLOT_UPDRIFT_MAX_KM updrift of the groin, using the
       same multi-source
       chainage table as everything else (so it covers the full
       historical record, not just the CoastSat era). No geometric
       correction of any kind -- CoastSat's transects are shore-normal
       and run parallel to the coast, so chainage is used exactly as
       provided. GIF_FRAME_MODE picks "date" (default -- one frame per
       unique observation date, using every observation) or "year"
       (pools all sources/dates within a calendar year -- fuller
       per-frame coverage, coarser resolution). Adapted from the
       whole-island version in HAT_shoreline_chainage_alldata_evolution.py.

Data sources (identical to shoreline inventory):
    Aerial wet-dry lines   - user-digitized (GeoJSON)
    NC Coastal Management  - historical shoreline shapefile (GeoJSON)
    CoastSat               - satellite-derived, one CSV per transect
    (CoastSat imagery doesn't reach back before ~1984, so the
    pre-installation baseline is necessarily wet-dry + NC state only --
    that's expected, not a bug; see the per-source breakdown printed by
    compute_preinstall_lrr.)

Chainage extraction:
    CoastSat  - read directly from its per-transect CSVs
    Wet-dry   - geometric intersection with transect line (shapely)
    NC state  - same intersection approach

All shoreline sources are treated as equal-weight observations in
regression. Wet-dry (~1-2 m precision) is more accurate than CoastSat
(~10 m), so weighted regression is a reasonable extension for later.

Distance from groin is measured ALONG THE SHORE (curvilinear, the same
cumulative alongshore distance used everywhere else), not as a
straight-line northing difference -- see assign_distance_from_groin().
The origin (0 km) is the NORTHERNMOST individual groin feature, since
updrift = north; the whole groin field (first to last groin) is shaded
on the profile plots. CASCADE domain number is the PRIMARY x-axis on
the profile plots and the GIF, with distance-from-groin (km) on a
secondary top axis -- see set_domain_primary_axis().

Curves vs. smoothing: no binning or averaging anywhere -- every series
(pre-install baseline, full-period, post-install, each era, era-to-era
differences) is real per-transect data, connected in distance order so
the trend reads as a line rather than a loose cloud of dots. Binning
was deliberately removed: even a non-smoothing fixed-width bin pools
multiple transects together and can shift or mute a real, spatially
localized groin signal.

The per-era profile plot (plot_era_lrr_profile) additionally draws a
LOWESS-smoothed overlay (SMOOTHED_OVERLAY_LOWESS_FRAC bandwidth, a
lighter shade of each line's own color) as a supplementary visual aid
-- it is drawn ALONGSIDE the raw connected points, never in place of
them, so the actual per-transect data is always visible underneath.

The regional pre-install baseline used for the DECADAL ANOMALY
calculation (a different thing from the plots above) connects the
sorted per-transect pre-install points directly via linear
interpolation -- see interp_preinstall_baseline().

Outputs
-------
  groin_analysis_chainage_all.csv                 - unified (transect, date, source, chainage) records
  groin_analysis_preinstall_lrr.csv               - per-transect pre-install LRR (nourishment-flagged, not dropped)
  groin_analysis_full_period_lrr.csv              - per-transect full-period LRR + piecewise breakpoint
  groin_analysis_post_install_lrr.csv             - per-transect post-install-only (1970-present) LRR
  groin_analysis_decadal_lrr.csv                  - long-format decadal (fixed, non-overlapping) LRR [no plot]
  groin_analysis_era_lrrs.csv                     - per-transect per-era LRR
  groin_analysis_decadal_anomaly.csv              - decadal LRR minus interpolated pre-install baseline [no plot]
  groin_analysis_signal_extent.csv                - contiguous signal extent by decade [no plot]
  groin_analysis_lrr_distance_band_series.csv     - decadal LRR aggregated into distance bands x decade [no plot]
  groin_analysis_anomaly_distance_band_series.csv - decadal anomaly aggregated into distance bands x decade [no plot]
  groin_analysis_alongshore_profile.png           - LRR vs. CASCADE domain / distance from groin (pre-install vs. full-period)
  groin_analysis_era_lrr_profile.png              - LRR vs. CASCADE domain / distance from groin, one curve per era
                                                     (+ pre-2022-nourishment subset, + LOWESS-smoothed overlay per line)
  groin_analysis_era_difference_profile.png       - era-to-era Δ LRR, both comparisons together (Functional-Pre-install, Deteriorated-Functional)
  groin_analysis_diff_pre_to_functional.png       - Δ LRR, Functional groin minus Pre-install baseline (focused single-comparison version)
  groin_analysis_diff_functional_to_deteriorated.png - Δ LRR, Deteriorated minus Functional groin (focused single-comparison version)
  groin_analysis_diff_deteriorated_to_preinstall.png - Δ LRR, Deteriorated minus Pre-install baseline (whole historical arc, skipping the middle era)
  groin_analysis_decade_increment_lrrs_{N}yr.csv  - per-transect LRR for each N-year window since installation (one file per DECADE_PLOT_INCREMENTS_YEARS entry)
  groin_analysis_decade_increment_profile_{N}yr.png - LRR profile per post-install window (see above), color = chronological order
  groin_analysis_zone_panels.png                  - one stacked panel per era, downdrift/updrift analysis zones shaded with mean-rate labels
  groin_analysis_shoreline_evolution.gif          - animated shoreline position (chainage) through time, groin area --
                                                     water/island shaded around the curve, y-range fixed to the full record
  groin_analysis_shoreline_evolution_zoomed.gif    - same, but zoomed to domain 1-GIF_ZOOM_WINDOW_DOMAIN_MAX for a closer look around the groin
  groin_analysis_metadata.txt                     - run metadata

(The distance-band, anomaly-profile, signal-extent, and transect-
diagnostic PLOTS were removed entirely -- not just unused -- since they
weren't helpful. Their underlying compute_* functions and CSV exports
are still here; only the plotting functions are gone.)

Usage
-----
Edit CONFIG below, then run:
    python HAT_groin_shoreline_analysis_v2.py
"""

# ============================================================
# CONFIG
# ============================================================

# ─── Paths ────────────────────────────────────────────────────────────
BASE = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\groin_module"

STUDY_AREA_FILTER_PATH = r"/scripts/groin_module_noBE/gis_analysis/gis_data/cascade_area.geojson"
STUDY_AREA_BUFFER_M    = 500

WET_DRY_PATH     = r"/scripts/groin_module_noBE/gis_analysis/gis_data/wet_dry_groin.geojson"
WET_DRY_DATE_COL = "date"

NC_STATE_PATH     = r"/scripts/groin_module_noBE/gis_analysis/gis_data/nc_shorelines.geojson"
NC_STATE_DATE_COL = "DATE_"

COASTSAT_TRANSECT_GEOM   = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\CoastSat_transect_layer.geojson"
COASTSAT_TRANSECT_ID_COL = "id"
COASTSAT_ROOT_DIR        = r"/scripts/input_prep/CoastSat/coastsat_timeseries"

# ─── CASCADE domain reference (authoritative) ──────────────────────────
# A GeoJSON of the 90 real CASCADE domain boxes (D1-D90), each a 500 m
# polygon with a "domain_id" property. This REPLACES the original formula-
# based domain assignment (CAPE_POINT_NORTHING + fixed 500 m spacing),
# which was measurably wrong -- domain 1's real midpoint northing is
# 3,899,029 per this file, not the ~3,897,750 the original formula assumed
# (over 2 domains off). See load_domain_reference() /
# assign_domain_from_northing().
DOMAINS_JSON_PATH = r"/scripts/groin_module_noBE/gis_analysis/gis_data/HAT_domains.json"

OUTPUT_DIR = r"/scripts/groin_module_noBE/gis_analysis/shoreline_output_grid100m"

# CRS for all spatial operations (UTM 18N covers NC Outer Banks)
PROJECTED_CRS = "EPSG:32618"

# ─── Plot x-axis limit ─────────────────────────────────────────────────
# Distance-from-groin plots (alongshore profile, era profile) are
# clipped to start exactly at CASCADE domain 1 (the southernmost real
# domain -- see domain_dist_km()) on the downdrift side, and to
# PLOT_UPDRIFT_MAX_KM on the updrift side.
PLOT_UPDRIFT_MAX_KM = 30
PLOT_X_MAX_KM = 60   # generic safety cap; kept as the GIF/other fallback default

# ─── Groin location ────────────────────────────────────────────────────
# The Buxton Groin Field is provided as a shapefile/GeoJSON with
# multiple LineString features (the individual groins). The script
# loads the geometry and uses:
#   - the NORTHERNMOST groin's northing as the origin for "distance
#     from groin" (0 km), since updrift = north here
#   - the full north-south extent (southernmost to northernmost groin)
#     as a shaded band on plots, showing the whole groin field's footprint
GROIN_GEOJSON_PATH = r"/scripts/groin_module_noBE/gis_analysis/gis_data/groins_hatteras.geojson"

# Fallback northing if the geojson can't be loaded. Used only if the
# geojson is missing. Set to None to error out instead.
GROIN_NORTHING_FALLBACK = 3_901_580.0    # from groins_hatteras.geojson center

# Direction of "updrift" — the side where sediment approaches from.
# On Hatteras Island, net alongshore transport is southward, so updrift
# is NORTH of the groin. Signed distance from the groin is then
# positive = north = updrift, negative = south = downdrift.
UPDRIFT_DIRECTION = "north"    # "north" or "south"

# ─── Groin installation year (used as fixed breakpoint reference,
#     and as the cutoff for the pre-installation baseline below) ──────
GROIN_INSTALLATION_YEAR = 1970

# ─── Pre-installation baseline ────────────────────────────────────────
# ALL shoreline observations dated strictly before GROIN_INSTALLATION_YEAR
# are used as the pre-groin baseline -- not a hand-picked list of years,
# and pooled across EVERY source that has pre-1970 data (aerial wet-dry
# lines and/or NC Coastal Management shorelines; CoastSat doesn't reach
# back this far). Nothing is dropped for being pre-1970 -- the code
# prints a per-source breakdown so it's visible exactly what went in.
PRE_INSTALLATION_YEAR_CUTOFF = GROIN_INSTALLATION_YEAR

# Nourishment zones by (name, year, y_min, y_max): the 1966 Buxton
# nourishment placed sand directly in front of the groin field, so the
# 1966 and 1967 shorelines there reflect a fill event, not the natural
# pre-install trend. Per-committee direction: keep this data in the
# regression (don't silently drop real observations) -- these
# transects are only FLAGGED (nourishment_zone=True in the comparison
# tables, shown as a distinct marker on the profile plot) so a reader
# can see exactly which points are nourishment-influenced and judge for
# themselves, rather than having them removed upstream.
NOURISHMENT_EXCLUSIONS = [
    ("Buxton 1966",   1966, 3_899_000, 3_907_000),
    ("Buxton 1966 → 1967 shoreline", 1967, 3_899_000, 3_907_000),  # 1966 nourishment biases 1967
]

# ─── Analysis periods (for era-averaged LRR profiles) ────────────────
# Three eras derived from groin structural/maintenance condition
# (CSE 2013 documented history -- see HAT_groin_zone_investigation.py
# for full citations):
#   Pre-install     - 1849 to 1969 (baseline, no structure)
#   Functional groin - 1970 to 1995 (structure actively maintained:
#                      built 1970, storm/construction damage same
#                      year, 1975 steel-pile repair on groin #1,
#                      1980/82 anti-flanking extensions, 1994 Gordon
#                      damage, 1995 south-groin repair)
#   Deteriorated    - 1996 to 2024 (LAST documented repair was the
#                      1995 south-groin steel sheet piling -- nothing
#                      after that. Hurricane Isabel (2003) is a real,
#                      notable event WITHIN this era, not the era
#                      boundary itself, since no maintenance happened
#                      either right before or after it; it's annotated
#                      on the decadal plots instead of used as a split)
#
# NOTE: this merges what were previously two separate eras
# ("Post-Gordon" 1995-2003 and "Deteriorated" 2004-2024) into one,
# since the maintenance history doesn't actually support treating them
# as functionally distinct -- both are "no more repairs happened."
# If the piecewise-breakpoint search (compute_fullperiod_lrr) turns up
# a real, data-supported rate change clustering near 2003 across many
# transects, that would be evidence FOR re-splitting at Isabel --
# worth checking rather than assuming either way.
ERAS = [
    ("Pre-install",      1849, PRE_INSTALLATION_YEAR_CUTOFF - 1),
    ("Functional groin",  1970, 1995),
    ("Deteriorated",      1996, 2024),
]

# Notable storm/repair events -- shown as annotated vertical lines on
# the decadal plots even though they no longer define era boundaries.
NOTABLE_EVENTS = [
    (1970, "Groins built"),
    (1994, "Gordon damage"),
    (1995, "Last repair (S. groin)"),
    (2003, "Isabel"),
    (2022, "Jug Handle Bridge → mgmt ends"),
]

# ─── Decadal window parameters (replaces the original sliding rolling window) ──
# Fixed, NON-OVERLAPPING 10-year bins are much easier to read and
# report than a continuously-slid 1-year-step window: in the original
# scheme, adjacent "center years" shared ~90% of their underlying
# observations (a 10-year window stepped by 1 year overlaps 9 of its
# 10 years with its neighbor), which made the resulting time series
# look smoother and more continuous than the underlying survey
# frequency actually supports, and made it hard to say what was
# genuinely new information between e.g. center-year 1994 and 1995.
# Decades start at DECADE_START_YEAR and run every DECADE_LENGTH_YEARS
# until the data runs out -- each number on a plot is now a real,
# distinct slice of time with no double-counted observations.
DECADE_START_YEAR      = 1960
DECADE_LENGTH_YEARS    = 10
DECADE_MIN_OBSERVATIONS = 5   # minimum observations in a decade to fit LRR

# ─── Signal-extent detection ──────────────────────────────────────────
# The alongshore extent of the groin signal is defined as the largest
# CONTIGUOUS run of domain-bins, starting immediately at the groin and
# working outward, whose |anomaly| (decadal LRR - pre-install baseline)
# exceeds SIGNAL_ANOMALY_THRESHOLD. Requiring contiguity from the groin
# (rather than "furthest point anywhere that crosses the threshold")
# stops a single noisy far-field transect from being reported as part
# of the groin signal. Signal extent is computed separately updrift and
# downdrift. See compute_signal_extent_over_time().
#
# Is 1.0 m/yr justified? It's a round, literature-typical number, and
# it IS checked against this site's own noise floor: compute_preinstall_lrr()
# computes the standard deviation of pre-install LRR at transects
# beyond SIGNAL_MAX_SEARCH_DISTANCE_M (far enough that the groin can't
# plausibly reach) and prints how many multiples of that std dev this
# threshold represents. A defensible threshold is roughly 1-2x that
# natural background standard deviation, so the detector isn't flagging
# ordinary shoreline variability as "groin signal." READ THAT PRINTED
# NUMBER on your first real run -- if the ratio isn't in a reasonable
# 1-2x range, adjust this constant to match your site's actual noise
# floor rather than leaving it at the round default.
SIGNAL_ANOMALY_THRESHOLD_M_YR = 1.0     # m/yr above baseline required
SIGNAL_MAX_SEARCH_DISTANCE_M  = 15_000  # cap on search distance from groin
SIGNAL_EXTENT_BIN_WIDTH_M     = 500     # bin width for contiguity check
                                          # (matches CASCADE domain length)
SIGNAL_EXTENT_MAX_GAP_BINS    = 1       # allow this many below-threshold
                                          # bins in a row before stopping
                                          # (absorbs single-bin noise)
# If a transect lacks a pre-install LRR (no historical shorelines on
# it), fall back to a regional baseline built from the transects that
# do have pre-install LRR -- see interp_preinstall_baseline(). This
# connects the sorted per-transect pre-install values directly (linear
# interpolation), not a binned step function or a smoothed curve. No
# LOWESS/kernel smoother is used anywhere in this script.

# Distance bands used for the "signal through time" curve plots that
# replace the transect x year heatmaps. Each band gets its own curve
# (decadal LRR, or anomaly, vs. decade) applied
# symmetrically updrift and downdrift.
DISTANCE_BAND_EDGES_M = [0, 1000, 3000, 6000, 10000, 15000]

# ─── Piecewise breakpoint search ──────────────────────────────────────
# The piecewise-linear regression per transect uses ALL observations
# (1849–present). This range is only the SEARCH space for the single
# breakpoint year. Extending it later lets the analysis catch regime
# changes near the recent record.
BREAKPOINT_SEARCH_START = 1970
BREAKPOINT_SEARCH_END   = 2020
BREAKPOINT_MIN_POINTS_PER_SEGMENT = 5

# ─── Transect chainage extraction ─────────────────────────────────────
# Length of transect line (from origin, seaward) used for shoreline
# intersection. CoastSat transects are typically ~500 m; a longer value
# handles cases where the actual shoreline sits further seaward than
# the CoastSat mean line.
TRANSECT_INTERSECTION_LENGTH_M = 800

# Maximum absolute chainage retained. Anything beyond this is treated
# as a spurious intersection (e.g., two shorelines happened to overlap
# far from the intended cross-shore location).
CHAINAGE_MAX_ABS_M = 700

# Minimum observations per transect to include it in analysis
MIN_OBSERVATIONS_PER_TRANSECT = 8

# ─── Shoreline evolution GIF (groin area) ──────────────────────────────
# Animated view of RAW shoreline position (chainage) through time,
# zoomed to the groin area -- no geometric correction of any kind:
# CoastSat's transects are shore-normal and run parallel to the coast,
# so chainage is used exactly as provided. Adapted from the whole-
# island version in HAT_shoreline_chainage_alldata_evolution.py, but
# built on THIS script's multi-source unified chainage table (wet-dry +
# NC state + CoastSat) so it can show the full historical record, not
# just the CoastSat era (1984+). See create_groin_evolution_gif().
#
# GIF_FRAME_MODE controls what one "frame" is:
#   "date" -- one frame per unique observation DATE, using every
#             observation the data supports (not collapsed into
#             years). Maximum temporal resolution, but a single
#             wet-dry or NC-state date, or an early/partial CoastSat
#             pass, often only covers PART of the window, so
#             individual frames can look spatially sparse -- that's
#             real data sparsity on that date, not a bug.
#   "year"  -- one frame per calendar year, pooling every source/date
#             within that year (median). Fuller per-frame spatial
#             coverage, coarser temporal resolution.
GIF_FRAME_MODE            = "date"    # "date" or "year"
# A date/year needs to clear BOTH of these bars to get a frame --
# lowered the fraction from the original 0.15, which was excluding most
# of the historical record: a single wet-dry/NC-state survey, or an
# early CoastSat pass, usually covers a small fraction of the window on
# its own, so a 15%-of-window bar was throwing out the majority of real
# observations. The absolute floor just guarantees a frame has at
# least a few real points to draw a meaningful curve from, even when
# the fraction bar alone would let through 1-2 lone transects.
GIF_MIN_TRANSECT_FRACTION = 0.75   # min fraction of in-window transects valid to draw a frame
GIF_MIN_TRANSECT_ABS      = 3      # AND at least this many actual transects
# Before CoastSat existed, every historical shoreline is precious and
# inherently sparse -- so ANY frame with at least GIF_MIN_TRANSECT_ABS
# transects gets kept regardless of coverage fraction, for any date
# before this year. GIF_MIN_TRANSECT_FRACTION only applies from this
# year onward, where CoastSat's dense modern coverage makes a
# meaningful fraction bar reasonable to enforce.
GIF_PRE_COASTSAT_CUTOFF_YEAR = 1984
# Second GIF, zoomed to domain 1 through this domain (instead of the
# full domain-1-to-PLOT_UPDRIFT_MAX_KM window) for a closer look at
# change right around the groin.
GIF_ZOOM_WINDOW_DOMAIN_MAX = 25
GIF_FRAME_DURATION_S      = 0.15   # seconds per frame (many frames now -- keep this brisk)
GIF_TRAIL_FRAMES          = 5      # preceding frames shown as a fading trail
GIF_DPI                   = 130
# Last pre-install shoreline in the record -- drawn as a fixed,
# always-visible reference line on every frame dated AFTER this year
# (not on frames before/during it, since it's meant as a baseline for
# what came after, not a comparison for the earlier record). Built the
# same way as any other frame -- median chainage per transect for this
# specific year, gaps left where no observation exists that year -- so
# it's exactly as sparse/honest as the real 1967 coverage actually is.
GIF_REFERENCE_YEAR        = 1967
# Background shading so each frame reads as a real coastal cross-
# section rather than an abstract line chart: everything seaward of
# the shoreline (larger chainage) shaded like water, everything
# landward (smaller chainage, back toward the transect origin) shaded
# like the island itself.
GIF_WATER_COLOR = "#8FC1E3"
GIF_LAND_COLOR  = "#E3D5A8"
# Requires Pillow: pip install pillow

# ─── GIF transect source ────────────────────────────────────────────────
# Which transects define alongshore position and/or measure shoreline
# position, for the GIF SPECIFICALLY -- every other analysis in this
# script always uses CoastSat's own shore-normal transects, since those
# measure true cross-shore distance regardless of which way they
# happen to face, which is what makes them correct for computing rates.
#
#   "coastsat" (default) -- CoastSat's own transects for BOTH the
#       chainage measurement (correct, shore-normal) and the alongshore
#       x-axis position (via its trusted ID order).
#   "hybrid" -- keeps CoastSat's own (correctly-measured) chainage
#       values exactly as they are, but repositions each one along the
#       x-axis using the nearest transects_100m.geojson transect's own
#       alongshore position instead of CoastSat's ID-order-based one --
#       perfectly even, precisely-known 100 m spacing on the x-axis,
#       without changing how the shoreline itself was measured.
#   "grid100m" -- fully re-measures the shoreline against the 100 m
#       grid's OWN transects (parallel, not shore-normal): each
#       chainage observation gets reconstructed into a real (x, y)
#       point using its original CoastSat transect, then re-projected
#       onto the nearest 100 m-grid transect. WARNING: the 100 m grid's
#       transects all point the same fixed direction rather than
#       following the coast's curve, which reintroduces a known
#       geometric artifact (a spurious drift in raw position wherever
#       the true coast diverges from that fixed direction -- see the
#       "IMPORTANT" note this caused several turns of debugging
#       earlier in this project). This option exists so you can see
#       that for yourself, not because it's recommended.
GIF_TRANSECT_SOURCE = "grid100m"   # "coastsat", "hybrid", or "grid100m"
TRANSECTS_100M_PATH = r"/scripts/groin_module_noBE/gis_analysis/gis_data/transects_100m.geojson"

# ─── CASCADE domain assignment ────────────────────────────────────────
# Domains are ~500 m long, numbered from south (D1) to north (D90).
# Real boundaries come from DOMAINS_JSON_PATH (see load_domain_reference()
# / assign_domain_from_northing()) -- NOT a formula anymore. The original
# formula (D1 midpoint at a fixed "CAPE_POINT_NORTHING") was measurably
# wrong: real D1 midpoint is at northing 3,899,029 per the domains file,
# over 2 domains away from what the formula assumed.
NUM_REAL_DOMAINS   = 90
DOMAIN_LENGTH_M    = 500

# ─── Plot styling ─────────────────────────────────────────────────────
SOURCE_COLORS = {
    "wet_dry":  "#9C4D51",   # muted brick red
    "nc_state": "#3E5F7E",   # steel blue
    "coastsat": "#C27B49",   # warm terracotta
}
SOURCE_LABELS = {
    "wet_dry":  "Aerial wet-dry lines",
    "nc_state": "NC Coastal Mgmt",
    "coastsat": "CoastSat (satellite)",
}

ERA_COLORS = {
    "Pre-install":       "#5A5A5A",
    "Functional groin":  "#2E7D32",
    "Deteriorated":      "#C62828",
}

# Extra display-only sub-period: Deteriorated era before the 2022
# nourishments (Buxton 2022, Avon 2022) -- shown alongside the full
# Deteriorated line on the era profile plot, NOT part of the official
# ERAS list (so it doesn't touch the anomaly/signal-extent pipeline,
# which is built on the 3 official eras only).
PRE_NOURISHMENT_PERIOD = ("1996-2021 (pre-2022 nourishment)", 1996, 2021)
PRE_NOURISHMENT_COLOR = "#F57C00"   # orange -- distinguishable from Deteriorated's red

# ─── Decade-increment LRR evolution plot ───────────────────────────────
# A finer-grained complement to the 3-era profile: shows the LRR
# profile across successive fixed-length windows since installation
# (e.g. 1970-1980, 1980-1990, ...), to see whether/how the groin's
# effect has changed continuously over its lifespan, not just as a
# single before/after-deterioration split.
DECADE_PLOT_START_YEAR      = GROIN_INSTALLATION_YEAR   # 1970
DECADE_PLOT_INCREMENTS_YEARS = [5, 10]   # generates one plot per increment in this list
DECADE_PLOT_END_YEAR        = 2024
DECADE_PLOT_COLORMAP        = "Greens"   # monochromatic ramp -- color = chronological order

# ─── Zone-panel profile plot ────────────────────────────────────────────
# Adapted from HAT_groin_zone_investigation.py's alongshore_profile_panels
# figure -- one stacked panel per era (using THIS script's own era
# definitions/data, not that script's separate rate pipeline), with
# shaded analysis zones and a mean-rate label per zone per panel.
ZONE_PANEL_DOWNDRIFT_DOMAINS = (1, 4)     # downdrift analysis zone
ZONE_PANEL_UPDRIFT_DOMAINS   = (7, 20)    # updrift analysis zone
ZONE_PANEL_X_MAX_DOMAIN      = 22         # a bit past the updrift zone, for padding
ZONE_PANEL_COLOR_DOWNDRIFT   = "#E67E22"   # orange
ZONE_PANEL_COLOR_UPDRIFT     = "#1565C0"   # blue
# Extra zone-panel rows: the Functional groin era restricted to just
# its first N years (starting at GROIN_INSTALLATION_YEAR), to see
# whether the groin looked most effective early on rather than across
# its whole 1970-1995 span. Empty list to skip.
ZONE_PANEL_FUNCTIONAL_EARLY_WINDOWS_YEARS = [5, 10]

# Used for any DIFFERENCE/delta line (era-to-era Δ LRR) -- deliberately
# NOT any color already used for a raw era curve (grey/green/red/
# orange), so a delta is never visually confusable with one of the two
# raw curves it was computed from.
DELTA_COLOR = "#6A3D9A"   # purple

# LOWESS bandwidth (fraction of points) for the smoothed-overlay lines
# on most profile plots -- a VISUAL AID only, drawn as a lighter
# shade of each line's own color, in addition to (not instead of) the
# raw per-transect points + connecting line. Lowered from 0.08 so the
# smoothed line tracks the raw per-transect data more closely (less
# aggressive smoothing) while still reading as a clean line.
SMOOTHED_OVERLAY_LOWESS_FRAC = 0.04

# Dedicated (lower) bandwidth for the structural-era profile
# specifically -- the sharp rate jumps right around the groin get
# smeared out at the shared 0.04 bandwidth, and this plot in
# particular is the one where those sharp local jumps are the point.
ERA_PROFILE_LOWESS_FRAC = 0.02

# ─── Text sizing (applies across all plots) ────────────────────────────
FONT_TITLE      = 15
FONT_AXIS_LABEL = 13
FONT_TICK       = 11
FONT_LEGEND     = 10.5
FONT_ANNOTATION = 10
FONT_LEGEND_GIF = 8.5   # the GIF's legend is denser (data + Ocean/Island + all annotations)
FONT_LEGEND_GIF_ZOOMED = 7   # zoomed GIF is narrower, so the same legend needs to be smaller still

# =============================================================================
# SECTION 6: GEOGRAPHIC ANNOTATION STYLING (annotated publication figure)
# =============================================================================
# Community/feature positions given as CASCADE domain numbers (or
# domain ranges), converted to distance-from-groin (km) via
# build_domain_to_dist_km() so they land in the right place regardless
# of which plot or x-window they're drawn on. See
# add_geographic_annotations().

ANN_TOWN_SPANS = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    # Tri-Village removed -- it and Rodanthe sit right at/beyond the
    # edge of the 30 km window and were rendering intersected with the
    # plot boundary rather than cleanly inside it.
}
ANN_VILLAGE_LINES = {}   # Salvo, Waves, Rodanthe all removed -- see notes above

# General cutoff: any annotation (village line, town span, shoal zone)
# past this domain is skipped or clipped -- e.g. "Waves" (domain 74)
# sits past this and is dropped entirely; a span like Wimble Shoals
# that starts before the cutoff but extends past it gets clipped to
# end AT the cutoff rather than being dropped outright.
ANN_MAX_DOMAIN = 70

ANN_PIER_LABEL_Y  = 0.76   # default rotated label y for any pier (0=bottom, 1=top axes fraction)
ANN_GROIN_LABEL_Y = 0.68   # rotated label y for groin lines
ANN_PIERS = {
    "Avon Pier":     (26, ANN_PIER_LABEL_Y),   # (domain, label_y) - adjust per pier
    # Rodanthe Pier removed -- same edge-of-window collision as above.
}
# NOTE: the Buxton Groin itself is NOT re-added here via ANN_GROINS --
# it's already the central, uniquely-marked feature of every plot in
# this script (the "Groin" line + "Groin field" shaded band), so a
# second label for the same structure would just be a duplicate.
ANN_WIMBLE_SHOALS = (60, 74)
ANN_AVON_SHOALS   = (24, 39)   # Avon Shoals influence zone

# Accretion / Erosion side labels - set to None for auto-computed midpoint,
# or a 0-1 axes fraction to pin to a fixed position.
LABEL_ACCRETION_Y = None
LABEL_EROSION_Y   = None

ANN_C_TOWN_SPAN    = "#90AFC5"
ANN_C_WIMBLE       = "#E0A800"   # amber - both shoal zones share this color
ANN_C_AVON_SHOALS  = "#E0A800"   # same amber as Wimble Shoals (same feature type)
ANN_C_VILLAGE_LINE = "0.40"
ANN_C_PIER         = "#1565C0"
ANN_C_GROIN        = "#B71C1C"
# Darker outline colors for the town-span / shoal shaded regions --
# their fill colors are close to the GIF's ocean/island background, so
# a plain fill alone can be hard to see; a darker border makes the
# extent of each region clear regardless of what's underneath it.
ANN_C_TOWN_SPAN_EDGE = "#2F5A73"
ANN_C_SHOAL_EDGE     = "#8A6800"

ANN_MODEL_COLOR = "#FF8C00"   # warm orange - modeled shoreline change rate


# ============================================================
# IMPORTS
# ============================================================
import os
import glob
import json
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import LineString, Point, MultiPoint, MultiLineString
from shapely.strtree import STRtree

import matplotlib.pyplot as plt
import matplotlib.collections
import matplotlib.dates as mdates
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import Normalize, to_rgb
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory

warnings.filterwarnings("ignore")


# ============================================================
# UTILITIES
# ============================================================

def _require_file(path: str, description: str):
    """Check a configured file path exists before handing it to
    geopandas/pyogrio, which otherwise raises a deep, unhelpful GDAL
    traceback for a plain "you haven't saved this file here yet" typo.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"\n\n  Could not find {description} at:\n    {path}\n"
            f"  This path is set in the CONFIG section near the top of "
            f"this script. Either save the file to that exact location, "
            f"or update the config variable to point at wherever you "
            f"actually saved it.\n")


_CRS_LOG = []   # populated by fix_crs(); printed via print_crs_summary()


def fix_crs(gdf: gpd.GeoDataFrame, source_name: str) -> gpd.GeoDataFrame:
    """Reproject to PROJECTED_CRS, handling missing or invalid CRS.
    Logs native CRS + reprojected bounds to _CRS_LOG for the
    consolidated spatial-reference confirmation (see print_crs_summary())."""
    native_crs = gdf.crs
    if native_crs is None:
        print(f"  ! {source_name}: no CRS defined, assuming EPSG:4326")
        gdf = gdf.set_crs("EPSG:4326")
        native_crs = gdf.crs
    try:
        reprojected = gdf.to_crs(PROJECTED_CRS)
        bounds = reprojected.total_bounds if len(reprojected) else None
        _CRS_LOG.append({
            "source": source_name, "native_crs": str(native_crs),
            "target_crs": PROJECTED_CRS, "ok": True, "bounds": bounds,
        })
        return reprojected
    except Exception as e:
        print(f"  ! {source_name}: CRS reprojection failed ({e}), "
              f"falling back to EPSG:32618")
        fallback = gdf.set_crs("EPSG:32618", allow_override=True)
        _CRS_LOG.append({
            "source": source_name, "native_crs": str(native_crs),
            "target_crs": "EPSG:32618 (fallback -- reprojection failed)",
            "ok": False, "bounds": fallback.total_bounds if len(fallback) else None,
        })
        return fallback


def print_crs_summary():
    """Consolidated confirmation of the spatial reference system for
    every shoreline/geometry file loaded this run: native CRS as read
    from the file, confirmation it was reprojected to PROJECTED_CRS,
    and the reprojected bounding box so a wrong/garbled reprojection
    (wildly out-of-range coordinates) is obvious at a glance rather
    than silently corrupting downstream distances.
    """
    print(f"\n{'='*72}\nSpatial reference system check")
    print(f"  Working CRS for all analysis: {PROJECTED_CRS}")
    if not _CRS_LOG:
        print("  ! no files logged yet")
        return
    for entry in _CRS_LOG:
        status = "OK" if entry["ok"] else "FAILED -- used fallback, verify manually"
        b = entry["bounds"]
        bounds_str = (f"E {b[0]:,.0f}-{b[2]:,.0f}  N {b[1]:,.0f}-{b[3]:,.0f}"
                      if b is not None else "n/a")
        print(f"  {entry['source']:<30} native={entry['native_crs']:<18} "
              f"reprojected bounds: {bounds_str}  [{status}]")
    print(f"  (All bounds should fall in a similar range if every layer "
          f"is genuinely covering the same stretch of coast -- a layer "
          f"with wildly different Easting/Northing here means its "
          f"reprojection is suspect, check that file's native CRS.)")


def parse_date_column(series: pd.Series, column_name: str) -> pd.Series:
    """Auto-detect date encoding and parse to pandas datetime.
    Handles: string dates, Unix milliseconds since epoch, year-only ints.
    """
    if pd.api.types.is_numeric_dtype(series):
        max_val = series.abs().max()
        if max_val > 1e9:
            # Unix milliseconds since epoch
            print(f"  '{column_name}': Unix milliseconds → parsing with unit='ms'")
            return pd.to_datetime(series, unit="ms", errors="coerce")
        elif max_val < 3000:
            # Bare year integer
            print(f"  '{column_name}': bare year integer → constructing YYYY-01-01")
            return pd.to_datetime(series.astype("Int64").astype(str) + "-01-01",
                                  errors="coerce")
    # String dates (assume ISO or common US formats)
    return pd.to_datetime(series, errors="coerce")


def decimal_year(dt: pd.Series) -> pd.Series:
    """Convert datetime to decimal year (e.g., 1985-07-02 → ~1985.50)."""
    dt = pd.to_datetime(dt)
    years = dt.dt.year
    year_starts = pd.to_datetime(years.astype(str) + "-01-01")
    year_ends = pd.to_datetime((years + 1).astype(str) + "-01-01")
    frac = (dt - year_starts) / (year_ends - year_starts)
    return years + frac


def load_domain_reference() -> pd.DataFrame:
    """Load the authoritative CASCADE domain boundaries (DOMAINS_JSON_PATH)
    -- 90 real domain polygons (D1-D90), each ~500 m alongshore, with a
    domain_id property. Replaces the original CAPE_POINT_NORTHING-based
    formula, which was measurably wrong (real D1 midpoint is at northing
    3,899,029 per this file, over 2 domains away from what the formula
    assumed).

    Returns DataFrame: domain_id, y_min, y_max, y_mid -- sorted by
    domain_id.
    """
    print(f"\nLoading CASCADE domain reference: "
          f"{os.path.basename(DOMAINS_JSON_PATH)}")
    _require_file(DOMAINS_JSON_PATH,
                  "the CASCADE domain reference (DOMAINS_JSON_PATH)")
    gdf = gpd.read_file(DOMAINS_JSON_PATH)
    print(f"  {len(gdf)} domain polygons, native CRS: {gdf.crs}")
    gdf = fix_crs(gdf, "CASCADE domains")

    id_col = "domain_id" if "domain_id" in gdf.columns else gdf.columns[0]
    records = []
    for _, row in gdf.iterrows():
        minx, miny, maxx, maxy = row.geometry.bounds
        records.append({"domain_id": int(row[id_col]), "y_min": miny, "y_max": maxy})
    domain_table = pd.DataFrame(records).sort_values("domain_id").reset_index(drop=True)
    domain_table["y_mid"] = (domain_table["y_min"] + domain_table["y_max"]) / 2.0

    print(f"  Domain range: D{domain_table['domain_id'].min()} to "
          f"D{domain_table['domain_id'].max()}")
    print(f"  Northing range: {domain_table['y_min'].min():,.0f} to "
          f"{domain_table['y_max'].max():,.0f}")
    return domain_table


def assign_domain_from_northing(y, domain_table: pd.DataFrame):
    """Assign northing(s) y (scalar or array-like) to a CASCADE domain
    number using the authoritative domain midpoints from
    load_domain_reference() -- nearest-midpoint assignment, so the
    small (~few meter) gaps/overlaps between the source polygons don't
    leave any northing unassigned. Values beyond the outermost real
    domain (buffer-zone transects) clip to D1 or D90 rather than
    extrapolating.
    """
    y_arr = np.atleast_1d(np.asarray(y, dtype=float))
    mids = domain_table["y_mid"].values
    ids = domain_table["domain_id"].values
    order = np.argsort(mids)
    mids_sorted = mids[order]
    ids_sorted = ids[order]
    idx = np.searchsorted(mids_sorted, y_arr)
    idx = np.clip(idx, 1, len(mids_sorted) - 1)
    left = idx - 1
    choose_left = np.abs(y_arr - mids_sorted[left]) <= np.abs(y_arr - mids_sorted[idx])
    nearest = np.where(choose_left, left, idx)
    result = ids_sorted[nearest].astype(int)
    return int(result[0]) if np.isscalar(y) or np.ndim(y) == 0 else result


def domain_dist_km(transects: gpd.GeoDataFrame, domain_number: int) -> float:
    """Mean dist_from_groin_m (km) of transects assigned to the given
    CASCADE domain -- used to anchor a plot's x-limit to a specific
    domain (e.g. "start the x-axis at domain 1") instead of wherever
    the data happens to start or stop. Returns None if that domain has
    no transects in this dataset.

    Prints the transect count and computed position so it's directly
    checkable in the run log -- if only one or two transects happen to
    be assigned to that domain (e.g. near the Cape Point bend, where
    domain assignment can be less certain), the resulting position may
    look further out than a naive "domain_number x 500 m" estimate
    would suggest; that's a real reflection of sparse/awkward coverage
    there, not necessarily a bug.
    """
    sub = transects[transects["domain"] == domain_number]
    if len(sub) == 0:
        print(f"  ! domain_dist_km: no transects assigned to domain "
              f"D{domain_number} -- falling back to auto-detected x-limit")
        return None
    pos_km = float(sub["dist_from_groin_m"].mean()) / 1000.0
    print(f"  domain_dist_km: D{domain_number} -> {pos_km:+.2f} km "
          f"({len(sub)} transect(s): "
          f"{sorted(sub['transect_id'].astype(str).tolist())[:5]}"
          f"{'...' if len(sub) > 5 else ''})")
    return pos_km


def build_domain_to_dist_km(transects: gpd.GeoDataFrame):
    """Build a callable domain_number -> dist_from_groin_m (km) mapping
    (linear interpolation across whichever domains have transects), and
    its inverse dist_km -> domain_number. Shared by
    set_domain_primary_axis() and add_geographic_annotations() so every
    domain-based reference on every plot uses the exact same mapping.

    Returns (domain_to_dist, dist_to_domain), both vectorized.
    """
    agg = (transects.groupby("domain")["dist_from_groin_m"]
                     .mean().sort_index())
    domains = agg.index.values.astype(float)
    dist_km = agg.values / 1000.0
    order = np.argsort(dist_km)
    dist_sorted = dist_km[order]
    dom_sorted = domains[order]

    def dist_to_domain(x):
        return np.interp(x, dist_sorted, dom_sorted)

    def domain_to_dist(d):
        return np.interp(d, domains, dist_km)

    return domain_to_dist, dist_to_domain


def set_domain_primary_axis(ax, transects: gpd.GeoDataFrame):
    """Make CASCADE domain number the PRIMARY (bottom) x-axis label,
    moving distance-from-groin (km) to a secondary (top) axis instead.
    The underlying plotted x-COORDINATE is unchanged -- still distance-
    from-groin in km, which is what all the binning/curve positions use
    -- this only changes which unit is the main axis label.

    Uses the domain assignment already computed per transect
    (transects['domain'], from the authoritative HAT_domains.json --
    see assign_domain_from_northing()), so it stays perfectly
    consistent with every other domain reference in this analysis.

    Call this AFTER the caller has already set the final x-limits
    (e.g. via apply_distance_xlim_and_autoscale_y), so only domains
    actually visible in that range get a tick.
    """
    agg = (transects.groupby("domain")["dist_from_groin_m"]
                     .mean().sort_index())
    if len(agg) < 2:
        return None
    dist_km = agg.values / 1000.0
    domain_to_dist, dist_to_domain = build_domain_to_dist_km(transects)

    x_lo, x_hi = ax.get_xlim()
    d_lo, d_hi = float(dist_to_domain(x_lo)), float(dist_to_domain(x_hi))
    d_lo, d_hi = min(d_lo, d_hi), max(d_lo, d_hi)
    step = 5
    domain_ticks = np.arange(np.floor(d_lo / step) * step, d_hi + step, step)
    domain_ticks = domain_ticks[domain_ticks >= 1]
    tick_positions_km = domain_to_dist(domain_ticks)
    # Exclude ticks within a small margin of the plot's edges -- a
    # label sitting right at the boundary (e.g. D70 at the very edge
    # of a 30 km window) tends to collide with geographic annotation
    # text placed near the edge, rather than reading as a normal tick.
    edge_margin_km = 0.6
    visible = ((tick_positions_km >= x_lo + edge_margin_km) &
                (tick_positions_km <= x_hi - edge_margin_km))

    # D65 is explicitly always shown if it's anywhere in the true
    # visible range (not shrunk by edge_margin_km) -- the margin
    # exclusion above is meant for ticks landing right at the very
    # edge (like D70), not for D65, which sits well within the window
    # for the standard 30 km-updrift view. IMPORTANT: domain_to_dist()
    # clips (doesn't extrapolate) queries outside the domains actually
    # present in `transects` -- on a narrower window (e.g. the zoomed
    # GIF, which only has domains up to ~25), domain_to_dist(65) would
    # silently return the position of whatever the real max domain is,
    # mislabeling it "D65". Guard against that by checking 65 is
    # actually within this transects table's own domain range first.
    domain_min, domain_max = agg.index.min(), agg.index.max()
    if 65 not in domain_ticks[visible] and domain_min <= 65 <= domain_max:
        km_65 = float(domain_to_dist(65))
        if x_lo <= km_65 <= x_hi:
            domain_ticks = np.append(domain_ticks, 65.0)
            tick_positions_km = domain_to_dist(domain_ticks)
            visible = ((tick_positions_km >= x_lo + edge_margin_km) &
                        (tick_positions_km <= x_hi - edge_margin_km))
            visible = visible | (domain_ticks == 65.0)

    ax.set_xticks(tick_positions_km[visible])
    ax.set_xticklabels([f"D{int(t)}" for t in domain_ticks[visible]], fontsize=FONT_TICK)
    ax.set_xlabel("CASCADE domain", fontsize=FONT_AXIS_LABEL)
    ax.tick_params(axis="y", labelsize=FONT_TICK)

    # Thin vertical line at EVERY individual domain's position (not
    # just the labeled ticks above, which only show every 5th/10th),
    # so each ~500 m domain boundary is visible on the plot itself.
    all_positions_km = dist_km  # one position per domain (order doesn't matter here)
    all_visible = (all_positions_km >= x_lo) & (all_positions_km <= x_hi)
    for pos in all_positions_km[all_visible]:
        ax.axvline(pos, color="#bbbbbb", lw=0.4, alpha=0.6, zorder=0)

    secax = ax.secondary_xaxis("top", functions=(lambda x: x, lambda x: x))
    secax.set_xlabel("Distance from groin (km, alongshore)  [+ = updrift]",
                       fontsize=FONT_AXIS_LABEL)
    secax.tick_params(labelsize=FONT_TICK)
    return secax


def lighten_color(hex_color: str, amount: float = 0.45):
    """Return a lighter shade of hex_color, blended toward white by
    `amount` (0 = unchanged, 1 = white). Used to give a smoothed-curve
    overlay a distinguishable shade of the same base color as its raw/
    connected-dots series, without hand-picking a second hex code for
    every line.
    """
    rgb = np.array(to_rgb(hex_color))
    return tuple(rgb + (np.array([1.0, 1.0, 1.0]) - rgb) * amount)


def add_geographic_annotations(ax, transects: gpd.GeoDataFrame,
                                gif_mode: bool = False):
    """Add community/feature context to a distance-from-groin plot:
    shaded town spans (Buxton, Avon), pier markers, and shoal-zone
    shading (Wimble Shoals, Avon Shoals) -- so a reader can see what's
    actually around a given stretch of coast, not just an abstract
    domain number. Positions come from ANN_* config (domain numbers/
    ranges), mapped to km via build_domain_to_dist_km() so they land
    correctly on whatever x-window this particular plot is showing.

    Styling matches the established convention from
    HAT_groin_hindcast_1984_2024.py: labels use a blended transform (x in
    data coordinates, y in AXES-FRACTION coordinates), which keeps
    label height fixed relative to the panel regardless of the y-data
    range or whether the axis is inverted (e.g. the GIF), and a white
    bbox behind every label so it stays legible even where a groin/pier
    line passes directly behind it.

    gif_mode=True switches to the bolder styling requested specifically
    for the GIF (higher alpha, and a higher zorder than the water/
    island fill_between so the shading renders ON TOP of it instead of
    being hidden underneath) -- the static plots use the subtler,
    hindcast-matching style by default.

    The groin itself is NOT re-labeled here -- it's already the
    central, uniquely-marked feature on every plot in this script.

    Call this AFTER the x-limits are set (same requirement as
    set_domain_primary_axis), and safe to call even if none of the
    annotated features fall within the current visible window (those
    are just skipped).
    """
    domain_to_dist, _ = build_domain_to_dist_km(transects)
    x_lo, x_hi = ax.get_xlim()
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    def visible(km):
        return x_lo <= km <= x_hi

    span_alpha = 0.30 if gif_mode else 0.14
    shoal_alpha = 0.30 if gif_mode else 0.10
    base_zorder = 2 if gif_mode else 0   # above the GIF's water/island fill (zorder=1)
    bbox = dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.82)

    # ─── Shoal zones (hatched fill + label) -- drawn first/underneath ──
    for name, (d_lo, d_hi), color in [
        ("Wimble Shoals", ANN_WIMBLE_SHOALS, ANN_C_WIMBLE),
        ("Avon Shoals", ANN_AVON_SHOALS, ANN_C_AVON_SHOALS),
    ]:
        if d_lo > ANN_MAX_DOMAIN:
            continue
        d_hi = min(d_hi, ANN_MAX_DOMAIN)
        km_lo, km_hi = float(domain_to_dist(d_lo)), float(domain_to_dist(d_hi))
        km_lo, km_hi = min(km_lo, km_hi), max(km_lo, km_hi)
        if km_hi < x_lo or km_lo > x_hi:
            continue
        clipped_lo, clipped_hi = max(km_lo, x_lo), min(km_hi, x_hi)
        ax.axvspan(clipped_lo, clipped_hi, color=color, alpha=shoal_alpha,
                    zorder=base_zorder, hatch="///", edgecolor=color, linewidth=0)
        mid_km = (clipped_lo + clipped_hi) / 2.0
        ax.text(mid_km, 0.04, name, transform=trans,
                 ha="center", va="bottom", fontsize=FONT_ANNOTATION,
                 color="#7A5800", style="italic", bbox=bbox, zorder=base_zorder + 10)

    # ─── Town spans (shaded band + label) -- centered by default, but ───
    # nudged off-center for any span where a pier sits near enough to
    # the midpoint to collide with a centered label (e.g. Avon Pier at
    # domain 26, the middle of Avon's 21-31 span).
    for name, (d_lo, d_hi) in ANN_TOWN_SPANS.items():
        if d_lo > ANN_MAX_DOMAIN:
            continue
        d_hi = min(d_hi, ANN_MAX_DOMAIN)
        km_lo, km_hi = float(domain_to_dist(d_lo)), float(domain_to_dist(d_hi))
        km_lo, km_hi = min(km_lo, km_hi), max(km_lo, km_hi)
        if km_hi < x_lo or km_lo > x_hi:
            continue
        clipped_lo, clipped_hi = max(km_lo, x_lo), min(km_hi, x_hi)
        ax.axvspan(clipped_lo, clipped_hi, color=ANN_C_TOWN_SPAN,
                    alpha=span_alpha, zorder=base_zorder)
        mid_domain = (d_lo + d_hi) / 2.0
        pier_conflict = any(abs(pd_ - mid_domain) < 0.15 * (d_hi - d_lo)
                             for pd_, _ in ANN_PIERS.values())
        frac = 0.25 if pier_conflict else 0.5
        label_km = clipped_lo + frac * (clipped_hi - clipped_lo)
        ax.text(label_km, 0.90, name, transform=trans,
                 ha="center", va="top", fontsize=FONT_ANNOTATION, color="0.25",
                 fontweight="bold", bbox=bbox, zorder=base_zorder + 10)

    # ─── Village center lines ───────────────────────────────────────────
    for name, d in ANN_VILLAGE_LINES.items():
        if d > ANN_MAX_DOMAIN:
            continue
        km = float(domain_to_dist(d))
        if not visible(km):
            continue
        ax.axvline(km, color=ANN_C_VILLAGE_LINE, lw=0.9, linestyle="--",
                    alpha=0.65, zorder=base_zorder + 1)
        ax.text(km, 0.84, name, transform=trans,
                 ha="center", va="top", fontsize=FONT_ANNOTATION, color="0.30",
                 bbox=bbox, zorder=base_zorder + 10)

    # ─── Piers ──────────────────────────────────────────────────────────
    for name, (d, label_y_frac) in ANN_PIERS.items():
        if d > ANN_MAX_DOMAIN:
            continue
        km = float(domain_to_dist(d))
        if not visible(km):
            continue
        ax.axvline(km, color=ANN_C_PIER, lw=1.0, linestyle="-.",
                    alpha=0.80, zorder=base_zorder + 2)
        ax.text(km, label_y_frac, name, transform=trans, rotation=90,
                 ha="center", va="top", fontsize=FONT_ANNOTATION,
                 color=ANN_C_PIER, bbox=bbox, zorder=base_zorder + 10)


def annotation_legend_handles():
    """Consolidated proxy legend artists for add_geographic_annotations()
    -- one entry per LAYER TYPE (not one per named feature), matching
    the established convention from HAT_groin_hindcast_1984_2024.py. Append
    these to a plot's own legend handles rather than relying on labels
    attached to the individual axvspan/axvline calls, which would
    otherwise add one legend entry per town/shoal/pier and make an
    already-busy legend worse.
    """
    handles = [Patch(fc=ANN_C_TOWN_SPAN, alpha=0.30, label="Community")]
    if ANN_WIMBLE_SHOALS[0] <= ANN_MAX_DOMAIN or ANN_AVON_SHOALS[0] <= ANN_MAX_DOMAIN:
        handles.append(Patch(fc=ANN_C_WIMBLE, alpha=0.25, hatch="///",
                               edgecolor=ANN_C_WIMBLE, linewidth=0,
                               label="Shoals position (Avon / Wimble)"))
    if ANN_VILLAGE_LINES:
        handles.append(Line2D([0], [0], color=ANN_C_VILLAGE_LINE, lw=0.9,
                                ls="--", label="Village center"))
    if ANN_PIERS:
        handles.append(Line2D([0], [0], color=ANN_C_PIER, lw=1.0, ls="-.",
                                label="Pier"))
    return handles


def _iter_axes_real_xy(ax):
    """Yield (xdata, ydata) arrays for every REAL data artist on ax --
    line plots and scatter points via their offsets, fill_between bands
    via their polygon vertices -- while skipping axhline/axvline
    reference lines (2-point constant lines using an axes-fraction
    blended transform on one axis, not real data)."""
    for line in ax.get_lines():
        xd, yd = np.asarray(line.get_xdata()), np.asarray(line.get_ydata())
        if len(xd) == 0:
            continue
        if len(xd) <= 2 and (np.ptp(xd) == 0 or np.ptp(yd) == 0):
            continue
        yield xd, yd
    for coll in ax.collections:   # scatter points + fill_between IQR bands
        if isinstance(coll, matplotlib.collections.PathCollection):
            # scatter: real data lives in offsets, not the marker path
            offsets = np.asarray(coll.get_offsets())
            if offsets.ndim == 2 and offsets.shape[0] > 0:
                yield offsets[:, 0], offsets[:, 1]
        else:
            # fill_between / PolyCollection: real data is the polygon
            # vertices themselves
            verts = np.concatenate([p.vertices for p in coll.get_paths()], axis=0) \
                    if coll.get_paths() else np.empty((0, 2))
            if len(verts) > 0:
                yield verts[:, 0], verts[:, 1]


def apply_distance_xlim_and_autoscale_y(ax, x_max_km: float = None,
                                         x_min_km: float = None,
                                         y_margin_frac: float = 0.10,
                                         x_data_margin_km: float = 1.0,
                                         y_percentile: tuple = (1, 99)):
    """Clip a distance-from-groin plot's x-axis and set the y-axis from
    the data actually VISIBLE within that window -- rather than an
    arbitrary fixed clip like the original +/-8 m/yr, or a symmetric +/-
    x_max_km that wastes space on the downdrift (short-coverage) side.

    x_max_km defaults to PLOT_X_MAX_KM (updrift coverage runs close to
    the full island). x_min_km defaults to None, which means: auto-
    detect the most negative (downdrift) x-value actually present in
    the plotted data, minus a small margin -- i.e. start the axis right
    where the real data starts instead of extending all the way to
    -x_max_km when downdrift coverage is much shorter than updrift.
    Pass an explicit x_min_km (e.g. -10) to override the auto-detect.

    y_percentile: the y-range is set from these percentiles of the
    visible raw data (not the strict min/max) -- a single extreme
    outlier point (which raw, unsmoothed per-transect data can have)
    would otherwise stretch the axis far past where nearly everything
    else sits, wasting most of the figure on empty space. Pass
    (0, 100) to restore strict min/max if you want every point
    guaranteed in view.

    Call this AFTER all plot/fill_between calls on the axes.
    """
    if x_max_km is None:
        x_max_km = PLOT_X_MAX_KM

    if x_min_km is None:
        data_x_min = np.inf
        for xd, _ in _iter_axes_real_xy(ax):
            finite = np.isfinite(xd)
            if finite.any():
                data_x_min = min(data_x_min, float(xd[finite].min()))
        x_min_km = (max(-x_max_km, data_x_min - x_data_margin_km)
                    if np.isfinite(data_x_min) else -x_max_km)

    ax.set_xlim(x_min_km, x_max_km)

    all_y_vis = []
    for xd, yd in _iter_axes_real_xy(ax):
        mask = (xd >= x_min_km) & (xd <= x_max_km)
        yd_vis = yd[mask]
        finite = np.isfinite(yd_vis)
        if finite.any():
            all_y_vis.append(yd_vis[finite])

    if all_y_vis:
        all_y_vis = np.concatenate(all_y_vis)
        y_min, y_max = np.percentile(all_y_vis, y_percentile)
        if y_max > y_min:
            pad = (y_max - y_min) * y_margin_frac
            ax.set_ylim(y_min - pad, y_max + pad)
    # else leave matplotlib's own autoscale in place


# ============================================================
# GROIN GEOMETRY
# ============================================================

def load_groin_geometry(domain_table: pd.DataFrame = None) -> dict:
    """Load the groin GeoJSON/shapefile and derive:
        reference_x, reference_y - centroid of the NORTHERNMOST
                       individual groin feature -- used as the origin
                       for distance-from-groin (0 km), since updrift =
                       north here
        y_min, y_max - alongshore extent across ALL groin features
                       (southernmost to northernmost), used to shade
                       the whole groin field's footprint on plots
        geometry    - the unioned geometry (for map overlays)
        n_features  - number of individual groin features
    Falls back to GROIN_NORTHING_FALLBACK if the file is missing.

    If domain_table is provided (see load_domain_reference()), prints
    which CASCADE domain the groin actually landed in -- a quick sanity
    check against wherever you expect the groin to be (e.g. "the groin
    is around domain 6/7"). A mismatch here usually means
    GROIN_GEOJSON_PATH points at the wrong file.
    """
    print(f"\nLoading groin geometry: {os.path.basename(GROIN_GEOJSON_PATH)}")
    if not os.path.exists(GROIN_GEOJSON_PATH):
        if GROIN_NORTHING_FALLBACK is None:
            raise FileNotFoundError(
                f"Groin file not found: {GROIN_GEOJSON_PATH}")
        print(f"  ! Not found — using GROIN_NORTHING_FALLBACK = "
              f"{GROIN_NORTHING_FALLBACK:,.0f}")
        return dict(reference_x=None, reference_y=GROIN_NORTHING_FALLBACK,
                    y_min=GROIN_NORTHING_FALLBACK,
                    y_max=GROIN_NORTHING_FALLBACK,
                    southernmost_x=None, southernmost_y=GROIN_NORTHING_FALLBACK,
                    middle_x=None, middle_y=GROIN_NORTHING_FALLBACK,
                    geometry=None,
                    n_features=0)

    gdf = gpd.read_file(GROIN_GEOJSON_PATH)
    print(f"  {len(gdf)} feature(s), native CRS: {gdf.crs}")
    gdf = fix_crs(gdf, "groin")
    geom_types = sorted(set(gdf.geometry.geom_type))
    print(f"  Geometry types: {geom_types}")
    if not any("Line" in gt for gt in geom_types):
        print(f"  ! WARNING: expected LineString feature(s) for a groin "
              f"structure, got {geom_types} instead. This usually means "
              f"GROIN_GEOJSON_PATH points at the wrong file -- double "
              f"check it's the groin shapefile, not e.g. the domain "
              f"reference or study area filter.")
    if len(gdf) > 20:
        print(f"  ! WARNING: {len(gdf)} features is a lot for a single "
              f"groin field ({GROIN_GEOJSON_PATH}) -- verify this is "
              f"really the groin file and not something else entirely.")

    merged = gdf.unary_union
    minx, miny, maxx, maxy = merged.bounds

    # Per-feature centroid, to find the NORTHERNMOST individual groin
    # feature specifically -- this becomes the distance-from-groin
    # origin (0 km), since updrift = north here. The overall span
    # (miny to maxy, across every feature) is kept separately for
    # shading the whole groin field's footprint.
    centroids = gdf.geometry.centroid
    northernmost_idx = centroids.y.idxmax()
    reference_x = float(centroids.x.loc[northernmost_idx])
    reference_y = float(centroids.y.loc[northernmost_idx])

    print(f"  Groin field extent (northing): {miny:,.0f} → {maxy:,.0f} "
          f"({(maxy - miny):.0f} m span across {len(gdf)} feature(s))")
    print(f"  Northernmost groin feature centroid (distance-from-groin "
          f"origin): x={reference_x:,.0f}, y={reference_y:,.0f}")

    if domain_table is not None:
        ref_domain = assign_domain_from_northing(reference_y, domain_table)
        span_domains = sorted(set(
            int(d) for d in assign_domain_from_northing(
                np.array([miny, maxy]), domain_table)))
        print(f"  → Northernmost groin feature is in domain D{ref_domain}; "
              f"the whole groin field spans domain(s) "
              f"D{span_domains[0]}-D{span_domains[-1]}"
              if len(span_domains) > 1 else
              f"  → Northernmost groin feature is in domain D{ref_domain}; "
              f"the whole groin field is within domain D{span_domains[0]}")
        print(f"  ^ if this doesn't match where you expect the groin to "
              f"be, GROIN_GEOJSON_PATH most likely points at the wrong "
              f"file.")

    # Southernmost and "middle" individual groin features too -- used
    # for lighter reference lines showing where each groin sits within
    # the field, not just the northernmost origin and the overall span.
    order = np.argsort(centroids.y.values)
    sorted_x = centroids.x.values[order]
    sorted_y = centroids.y.values[order]
    southernmost_x, southernmost_y = float(sorted_x[0]), float(sorted_y[0])
    mid_idx = len(sorted_y) // 2
    middle_x, middle_y = float(sorted_x[mid_idx]), float(sorted_y[mid_idx])

    return dict(reference_x=reference_x, reference_y=reference_y,
                y_min=miny, y_max=maxy,
                southernmost_x=southernmost_x, southernmost_y=southernmost_y,
                middle_x=middle_x, middle_y=middle_y,
                geometry=merged, n_features=len(gdf))


def assign_distance_from_groin(transects: gpd.GeoDataFrame,
                                groin: dict) -> gpd.GeoDataFrame:
    """Add a signed distance-from-groin column (in meters), measured
    ALONG THE SHORE (curvilinear alongshore_m), not as a straight-line
    northing difference. Hatteras Island bends sharply through Cape
    Point/Buxton; a straight north-south (northing) distance is not
    monotonic with true alongshore position there, which corrupts the
    x-axis of every distance-from-groin plot near the cape (visible as
    a fold-back/box artifact in profile plots). Using alongshore_m --
    the same cumulative along-path distance used everywhere else in
    this script -- keeps "distance from groin" consistent and
    monotonic with alongshore position.

    Positive = updrift, negative = downdrift, per UPDRIFT_DIRECTION.
    """
    sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
    tx = transects.copy()

    # Locate the groin's own alongshore position: the transect whose
    # NEAR-SHORE reference point (shore_x/shore_y) is closest to the
    # NORTHERNMOST groin feature's centroid (groin["reference_x"/"_y"]),
    # not the centroid of the whole groin field.
    if "shore_x" in tx.columns:
        ref_x, ref_y = tx["shore_x"], tx["shore_y"]
    else:
        ref_x, ref_y = tx["origin_x"], tx["origin_y"]

    gx = groin.get("reference_x")
    gy = groin["reference_y"]
    if gx is not None:
        d2 = (ref_x - gx) ** 2 + (ref_y - gy) ** 2
    else:
        # No geometry available -- fall back to nearest-by-northing
        d2 = (ref_y - gy) ** 2
    groin_idx = d2.idxmin()
    groin_alongshore_m = float(tx.loc[groin_idx, "alongshore_m"])

    tx["dist_from_groin_m"] = sign * (tx["alongshore_m"] - groin_alongshore_m)
    print(f"  Groin alongshore position: {groin_alongshore_m/1000:.2f} km "
          f"(nearest transect id={tx.loc[groin_idx, 'transect_id']}, "
          f"{d2.loc[groin_idx]**0.5:.0f} m away)")
    return tx


# ============================================================
# STUDY AREA FILTER
# ============================================================

def load_study_area_filter() -> gpd.GeoDataFrame:
    """Load and buffer the study area filter polygon."""
    print(f"\nLoading study area filter: {os.path.basename(STUDY_AREA_FILTER_PATH)}")
    _require_file(STUDY_AREA_FILTER_PATH, "the study area filter (STUDY_AREA_FILTER_PATH)")
    gdf = gpd.read_file(STUDY_AREA_FILTER_PATH)
    gdf = fix_crs(gdf, "study area filter")

    geom_types = set(gdf.geometry.geom_type)
    if geom_types & {"Polygon", "MultiPolygon"}:
        # Merge all polygons, then buffer by STUDY_AREA_BUFFER_M
        merged = gdf.unary_union.buffer(STUDY_AREA_BUFFER_M)
    else:
        # Lines: buffer to create corridor
        merged = gdf.unary_union.buffer(STUDY_AREA_BUFFER_M)

    filter_gdf = gpd.GeoDataFrame(geometry=[merged], crs=PROJECTED_CRS)
    bounds = filter_gdf.total_bounds
    print(f"  Filter polygon area: {filter_gdf.geometry.area.sum() / 1e6:.2f} km²")
    print(f"  Extent: {(bounds[2]-bounds[0])/1000:.1f} × "
          f"{(bounds[3]-bounds[1])/1000:.1f} km")
    return filter_gdf


# ============================================================
# LOADER: Shapefile / GeoJSON shorelines (wet-dry, NC state)
# ============================================================

def load_shoreline_lines(path: str, date_col: str, source_key: str,
                         filter_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Load a shapefile/GeoJSON of shoreline lines, project, spatially
    filter, and parse dates. Returns one row per feature with columns:
        geometry, date, year, source
    Multiple features on the same date are kept (needed for intersection).
    """
    print(f"\n[{SOURCE_LABELS[source_key]}] loading: {os.path.basename(path)}")
    if not os.path.exists(path):
        print(f"  ! File not found: {path}")
        return gpd.GeoDataFrame(columns=["geometry", "date", "year", "source"],
                                 geometry="geometry", crs=PROJECTED_CRS)

    gdf = gpd.read_file(path)
    print(f"  {len(gdf)} raw features, CRS: {gdf.crs}")
    gdf = fix_crs(gdf, source_key)

    # Spatial filter
    filter_geom = filter_gdf.geometry.iloc[0]
    keep_mask = gdf.geometry.intersects(filter_geom)
    gdf = gdf.loc[keep_mask].copy()
    print(f"  After spatial filter: {len(gdf)} features")

    if len(gdf) == 0:
        return gpd.GeoDataFrame(columns=["geometry", "date", "year", "source"],
                                 geometry="geometry", crs=PROJECTED_CRS)

    # Parse dates
    if date_col not in gdf.columns:
        available = ", ".join(gdf.columns[:15])
        raise ValueError(f"Date column '{date_col}' not found. "
                         f"Available: {available}")

    gdf["_parsed_date"] = parse_date_column(gdf[date_col], date_col)
    gdf = gdf.dropna(subset=["_parsed_date"])
    gdf["date"] = gdf["_parsed_date"]
    gdf["year"] = gdf["_parsed_date"].dt.year
    gdf["source"] = source_key

    print(f"  → {len(gdf)} valid features, "
          f"{gdf['year'].nunique()} unique years, "
          f"{gdf['date'].nunique()} unique dates")

    return gdf[["geometry", "date", "year", "source"]].reset_index(drop=True)


# ============================================================
# LOADER: CoastSat transect network (the analysis backbone)
# ============================================================

def load_coastsat_transects(filter_gdf: gpd.GeoDataFrame,
                             domain_table: pd.DataFrame) -> gpd.GeoDataFrame:
    """Load CoastSat's own transect layer and use it as the analysis
    backbone for EVERY source (wet-dry, NC state, and CoastSat itself)
    -- per your direction to use CoastSat's transects for all shoreline
    change analyses, replacing the earlier shared-100m-grid approach.

    Returns GeoDataFrame with columns:
        transect_id, geometry, origin_x, origin_y, dir_x, dir_y,
        shore_x, shore_y, alongshore_m, domain
    """
    print(f"\n[CoastSat transects] loading: "
          f"{os.path.basename(COASTSAT_TRANSECT_GEOM)}")
    _require_file(COASTSAT_TRANSECT_GEOM, "the CoastSat transect layer (COASTSAT_TRANSECT_GEOM)")
    gdf = gpd.read_file(COASTSAT_TRANSECT_GEOM)
    print(f"  {len(gdf)} global transects, native CRS: {gdf.crs}")
    gdf = fix_crs(gdf, "CoastSat transects")

    # Spatial filter using transect origin (start point)
    def get_origin(geom):
        coords = list(geom.coords)
        return Point(coords[0])

    gdf["origin"] = gdf.geometry.apply(get_origin)
    filter_geom = filter_gdf.geometry.iloc[0]
    keep_mask = gpd.GeoSeries(gdf["origin"], crs=PROJECTED_CRS).within(filter_geom)
    gdf = gdf.loc[keep_mask].copy()
    print(f"  In study area: {len(gdf)} transects")

    if len(gdf) == 0:
        raise RuntimeError("No CoastSat transects fell within study area — "
                           "check STUDY_AREA_FILTER_PATH and transect CRS.")

    # Compute origin and unit direction from geometry
    def unpack_geom(geom):
        coords = list(geom.coords)
        x0, y0 = coords[0]
        x1, y1 = coords[-1]
        dx, dy = x1 - x0, y1 - y0
        length = np.hypot(dx, dy)
        if length == 0:
            return x0, y0, 1.0, 0.0
        return x0, y0, dx / length, dy / length

    unpacked = gdf.geometry.apply(unpack_geom)
    gdf["origin_x"] = unpacked.apply(lambda t: t[0])
    gdf["origin_y"] = unpacked.apply(lambda t: t[1])
    gdf["dir_x"]    = unpacked.apply(lambda t: t[2])
    gdf["dir_y"]    = unpacked.apply(lambda t: t[3])
    gdf["transect_id"] = gdf[COASTSAT_TRANSECT_ID_COL].astype(str)
    # Near-shore reference point for groin-proximity matching (see
    # assign_distance_from_groin) -- for this short, beach-hugging
    # network, origin itself is already near shore.
    gdf["shore_x"] = gdf["origin_x"]
    gdf["shore_y"] = gdf["origin_y"]

    # Alongshore position: trusted ID order if possible, else NN fallback
    gdf = _compute_alongshore_positions(gdf, id_col=COASTSAT_TRANSECT_ID_COL)

    gdf["domain"] = assign_domain_from_northing(gdf["origin_y"].values, domain_table)

    print(f"  → alongshore extent: 0 → "
          f"{gdf['alongshore_m'].max()/1000:.1f} km")
    print(f"  → domain range: D{gdf['domain'].min()} to D{gdf['domain'].max()}")

    return gdf[[
        "transect_id", "geometry",
        "origin_x", "origin_y", "dir_x", "dir_y",
        "shore_x", "shore_y", "alongshore_m", "domain",
    ]].reset_index(drop=True)


def load_100m_grid_transects(domain_table: pd.DataFrame) -> gpd.GeoDataFrame:
    """Load the alternate 100 m transect grid (transects_100m.geojson)
    -- evenly-spaced, PARALLEL transects along a straight offshore
    reference line -- used ONLY by the GIF's "hybrid"/"grid100m"
    transect-source options (see GIF_TRANSECT_SOURCE). Not used
    anywhere else in this script; every other analysis always uses
    CoastSat's own shore-normal transects.

    Returns GeoDataFrame: transect_id, geometry, shore_x, shore_y
    (seaward end, on the straight reference line), origin_x, origin_y
    (landward end), dir_x, dir_y, alongshore_m (this grid's own
    northing-based alongshore position -- since it's evenly spaced
    along a STRAIGHT line, alongshore position is just northing raw_offset
    from the southernmost transect, no curvilinear ordering needed),
    domain (via the same authoritative domain_table used everywhere
    else in this script, so domain labels/annotations stay consistent
    regardless of which transect source the GIF is using).
    """
    print(f"\n[100m grid transects] loading: "
          f"{os.path.basename(TRANSECTS_100M_PATH)}")
    _require_file(TRANSECTS_100M_PATH,
                  "the 100 m transect grid (TRANSECTS_100M_PATH) -- "
                  "only needed for GIF_TRANSECT_SOURCE != 'coastsat'")
    gdf = gpd.read_file(TRANSECTS_100M_PATH)
    print(f"  {len(gdf)} transects, native CRS: {gdf.crs}")
    gdf = fix_crs(gdf, "100m grid transects")

    id_col = ("Transects_100m.LineID" if "Transects_100m.LineID" in gdf.columns
               else gdf.columns[0])
    gdf["transect_id"] = ("grid100m_" +
        pd.to_numeric(gdf[id_col], errors="coerce").astype("Int64").astype(str))

    def unpack(geom):
        coords = list(geom.coords)
        sx, sy = coords[0]    # seaward end (on the straight reference line)
        lx, ly = coords[-1]   # landward end
        # Direction points from the landward origin TOWARD the sea --
        # matching CoastSat's own convention (origin=landward,
        # direction=seaward), so "+chainage = accretion" means the same
        # thing for both transect sources. (seaward - landward), NOT
        # the reverse -- that sign was backwards before, which is why
        # grid100m chainage came out inverted.
        dx, dy = sx - lx, sy - ly
        length = np.hypot(dx, dy)
        if length == 0:
            return sx, sy, lx, ly, 1.0, 0.0
        return sx, sy, lx, ly, dx / length, dy / length

    unpacked = gdf.geometry.apply(unpack)
    gdf["shore_x"]  = unpacked.apply(lambda t: t[0])
    gdf["shore_y"]  = unpacked.apply(lambda t: t[1])
    gdf["origin_x"] = unpacked.apply(lambda t: t[2])
    gdf["origin_y"] = unpacked.apply(lambda t: t[3])
    gdf["dir_x"]    = unpacked.apply(lambda t: t[4])
    gdf["dir_y"]    = unpacked.apply(lambda t: t[5])

    gdf = gdf.sort_values("shore_y").reset_index(drop=True)
    gdf["alongshore_m"] = gdf["shore_y"] - gdf["shore_y"].min()
    gdf["domain"] = assign_domain_from_northing(gdf["shore_y"].values, domain_table)

    print(f"  → alongshore extent: 0 → "
          f"{gdf['alongshore_m'].max() / 1000:.1f} km")
    return gdf[["transect_id", "geometry", "shore_x", "shore_y",
                "origin_x", "origin_y", "dir_x", "dir_y",
                "alongshore_m", "domain"]].reset_index(drop=True)


def build_hybrid_position_map(transects: gpd.GeoDataFrame,
                               grid100m: gpd.GeoDataFrame,
                               groin: dict) -> pd.Series:
    """"hybrid" GIF transect-source: for each CoastSat transect, find
    the alongshore position it would have under the 100 m grid's
    coordinate system -- via nearest-neighbor matching (by real
    distance between each CoastSat transect's shore-side point and the
    100 m grid's shore-side points), NOT by reprojecting or re-
    measuring anything. This lets CoastSat's own (correctly-measured)
    chainage values be repositioned onto an evenly-spaced 100 m
    alongshore axis for display, without touching how the shoreline
    itself was measured.

    Returns a Series indexed by transect_id (str): dist_from_groin_m
    in the 100 m grid's coordinate system.
    """
    gx, gy = groin.get("reference_x"), groin["reference_y"]
    if gx is not None:
        d2 = (grid100m["shore_x"] - gx) ** 2 + (grid100m["shore_y"] - gy) ** 2
    else:
        d2 = (grid100m["shore_y"] - gy) ** 2
    groin_alongshore_100m = float(grid100m.loc[d2.idxmin(), "alongshore_m"])
    sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1

    grid_xy = grid100m[["shore_x", "shore_y"]].values
    grid_along = grid100m["alongshore_m"].values
    cs_xy = transects[["shore_x", "shore_y"]].values
    cs_ids = transects["transect_id"].astype(str).values

    results = {}
    for i, tid in enumerate(cs_ids):
        dists = np.hypot(grid_xy[:, 0] - cs_xy[i, 0], grid_xy[:, 1] - cs_xy[i, 1])
        nearest = int(np.argmin(dists))
        results[tid] = sign * (grid_along[nearest] - groin_alongshore_100m)
    return pd.Series(results, name="dist_from_groin_m_100mgrid")


def reconstruct_and_remeasure_on_grid100m(chainage_sub: pd.DataFrame,
                                           coastsat_transects: gpd.GeoDataFrame,
                                           grid100m: gpd.GeoDataFrame
                                           ) -> pd.DataFrame:
    """"grid100m" GIF transect-source: reconstruct each chainage
    observation into a real-world (x, y) point using its ORIGINAL
    (CoastSat) transect's origin/direction, then re-measure chainage
    against the NEAREST 100 m-grid transect's own origin/direction --
    restoring the shared-100m-grid approach from earlier in this
    project. WARNING: reintroduces the known parallel-transect
    geometric artifact (see GIF_TRANSECT_SOURCE config comment).

    Returns a new DataFrame with "transect_id" replaced by the matched
    100 m-grid transect ID and "chainage_m" replaced by the re-measured
    value.
    """
    tx = coastsat_transects.set_index("transect_id")[
        ["origin_x", "origin_y", "dir_x", "dir_y"]].copy()
    tx.index = tx.index.astype(str)

    sub = chainage_sub.copy()
    sub["transect_id"] = sub["transect_id"].astype(str)
    sub = sub[sub["transect_id"].isin(tx.index)].reset_index(drop=True)
    if len(sub) == 0:
        return sub

    row_info = tx.loc[sub["transect_id"]]
    ox, oy = row_info["origin_x"].values, row_info["origin_y"].values
    dx, dy = row_info["dir_x"].values, row_info["dir_y"].values
    ch = sub["chainage_m"].values
    px = ox + ch * dx
    py = oy + ch * dy

    grid_lines = grid100m.geometry.tolist()
    tree = STRtree(grid_lines)
    grid_origin = grid100m[["origin_x", "origin_y"]].values
    grid_dir = grid100m[["dir_x", "dir_y"]].values
    grid_ids = grid100m["transect_id"].values

    new_chainage = np.full(len(sub), np.nan)
    new_tid = np.full(len(sub), None, dtype=object)
    match_dist = np.full(len(sub), np.nan)
    for i in range(len(sub)):
        pt = Point(px[i], py[i])
        # True nearest-neighbor query (no distance cap) -- a buffer-
        # limited search (the previous approach) can silently drop
        # points that are legitimately far from any 100m-grid line,
        # which is exactly what happens near Cape Point/Buxton where
        # the parallel grid transects diverge most from CoastSat's
        # true shore-normal ones. This always finds SOME match.
        best_idx = int(tree.nearest(pt))
        match_dist[i] = grid_lines[best_idx].distance(pt)
        gox, goy = grid_origin[best_idx]
        gdx, gdy = grid_dir[best_idx]
        new_chainage[i] = (px[i] - gox) * gdx + (py[i] - goy) * gdy
        new_tid[i] = grid_ids[best_idx]

    if len(match_dist) > 0:
        print(f"  grid100m match distance (m): median {np.median(match_dist):.0f}, "
              f"95th pct {np.percentile(match_dist, 95):.0f}, "
              f"max {np.max(match_dist):.0f} "
              f"-- large values mean the reconstructed point sat far "
              f"from any 100m-grid transect (most likely near sharp "
              f"coastline curvature), not that data is missing")

    out = sub.copy()
    out["transect_id"] = new_tid
    out["chainage_m"] = new_chainage
    return out.dropna(subset=["chainage_m", "transect_id"]).reset_index(drop=True)


def _compute_alongshore_positions(gdf: gpd.GeoDataFrame,
                                   id_col: str = None) -> gpd.GeoDataFrame:
    """Cumulative alongshore distance for each transect.

    PRIMARY METHOD: if id_col parses cleanly as a numeric sequence, sort
    transects by that ID and accumulate REAL geometric distance between
    ID-CONSECUTIVE origins. CoastSat generates transect IDs in order
    along the digitized reference shoreline, so ID order IS true
    alongshore order -- trusting it sidesteps the failure mode below
    entirely.

    FALLBACK: greedy nearest-neighbor path from the southernmost
    origin. This is a heuristic over a raw point cloud and CAN go
    wrong: if the path runs out of nearby unvisited points it jumps
    across the island instead of along it, silently inflating
    alongshore_m for whatever transects get swept up afterward (this
    produced a real artifact previously — a handful of transects
    reporting distances-from-groin of 100+ km on an island that's
    actually ~45 km long). Any large jump is now flagged below
    regardless of which method was used, so it's visible in the run
    log rather than silently corrupting downstream plots.
    """
    gdf = gdf.copy()
    order = None
    if id_col is not None and id_col in gdf.columns:
        numeric_id = pd.to_numeric(gdf[id_col], errors="coerce")
        if numeric_id.notna().all() and numeric_id.nunique() == len(gdf):
            order = np.argsort(numeric_id.values)
            print(f"  Alongshore order: sorted by '{id_col}' "
                  f"(trusted sequential ID order)")
        else:
            # CoastSat's global transect IDs are often compound,
            # "{site}-{index}" (e.g. "usa_NC_0031-0160"), since a long
            # coastline gets split into multiple site chunks -- the
            # whole string doesn't parse as a plain number, but sorting
            # by (site prefix, numeric suffix) still recovers true
            # alongshore order, since CoastSat's own site numbering
            # runs along the coast too. This is trusted the same way
            # the plain-numeric case above is.
            id_str = gdf[id_col].astype(str).reset_index(drop=True)
            suffix_num = pd.to_numeric(
                id_str.str.extract(r"(\d+)$")[0], errors="coerce")
            prefix = id_str.str.replace(r"[-_]?\d+$", "", regex=True)
            if suffix_num.notna().all():
                sort_df = pd.DataFrame({"prefix": prefix, "suffix_num": suffix_num})
                order = sort_df.sort_values(["prefix", "suffix_num"]).index.values
                print(f"  Alongshore order: sorted by '{id_col}' site prefix + "
                      f"numeric suffix (trusted compound ID order, "
                      f"{prefix.nunique()} site chunk(s): "
                      f"{sorted(prefix.unique().tolist())})")

    coords = np.column_stack([gdf["origin_x"].values, gdf["origin_y"].values])
    n = len(coords)

    if order is None:
        print(f"  ! '{id_col}' didn't parse as a clean, unique numeric "
              f"sequence -- falling back to greedy nearest-neighbor path")
        start_idx = int(np.argmin(coords[:, 1]))
        order = [start_idx]
        visited = np.zeros(n, dtype=bool)
        visited[start_idx] = True
        for _ in range(n - 1):
            current = coords[order[-1]]
            dists = np.hypot(coords[:, 0] - current[0], coords[:, 1] - current[1])
            dists[visited] = np.inf
            nxt = int(np.argmin(dists))
            order.append(nxt)
            visited[nxt] = True
        order = np.array(order)

    ordered_coords = coords[order]
    segment_lengths = np.hypot(np.diff(ordered_coords[:, 0]),
                                np.diff(ordered_coords[:, 1]))

    # Flag suspiciously large jumps so they're visible in the run log
    # rather than silently inflating alongshore_m for later transects.
    if len(segment_lengths) > 0:
        typical = float(np.median(segment_lengths))
        bad = np.where(segment_lengths > max(20 * typical, 500))[0]
        if len(bad) > 0:
            print(f"  ! {len(bad)} suspiciously large alongshore jump(s) "
                  f"(>20x median spacing of {typical:.0f} m) -- inspect "
                  f"these transects before trusting far-field distances: "
                  f"{[str(gdf.iloc[order[b+1]].get('transect_id', order[b+1])) for b in bad[:10]]}"
                  f"{'...' if len(bad) > 10 else ''}")

    cumdist = np.concatenate([[0], np.cumsum(segment_lengths)])
    alongshore_map = dict(zip(order, cumdist))
    gdf["alongshore_m"] = [alongshore_map[i] for i in range(n)]
    return gdf


# ============================================================
# LOADER: CoastSat time-series CSVs
# ============================================================

def collect_csv_map(root_dir: str) -> dict:
    """Walk one level of subfolders under root_dir and return
    {csv_stem: full_filepath} for every CSV found."""
    csv_map = {}
    if not os.path.isdir(root_dir):
        print(f"  ! CoastSat root not found: {root_dir}")
        return csv_map
    for subfolder in sorted(os.listdir(root_dir)):
        sub_path = os.path.join(root_dir, subfolder)
        if not os.path.isdir(sub_path):
            continue
        for csv_file in sorted(glob.glob(os.path.join(sub_path, "*.csv"))):
            stem = Path(csv_file).stem
            csv_map[stem] = csv_file
    return csv_map


def _detect_coastsat_columns(df: pd.DataFrame) -> tuple:
    """Return (date_col, chainage_col), auto-detecting CoastSat's
    (slightly version-dependent) column names, falling back to
    "first column is date, second is chainage"."""
    date_col = next((c for c in ("dates", "date", "time", "datetime")
                       if c in df.columns), None)
    chainage_col = next((c for c in ("chainage", "chainage_m",
                                       "shoreline_position", "distance",
                                       "distance_m") if c in df.columns), None)
    if date_col is None:
        date_col = df.columns[0]
    if chainage_col is None:
        chainage_col = df.columns[1] if len(df.columns) > 1 else None
    return date_col, chainage_col


def load_coastsat_chainage(transects: gpd.GeoDataFrame,
                            root_dir: str) -> pd.DataFrame:
    """Read every CoastSat per-transect CSV directly, keyed by
    CoastSat's own transect ID -- straightforward again now that
    CoastSat's transects ARE the analysis backbone, so no point-
    reconstruction/re-snapping step is needed.

    Returns long-format DataFrame: transect_id, date, chainage_m, source
    """
    print(f"\n[CoastSat chainage] reading CSVs: {root_dir}")
    csv_map = collect_csv_map(root_dir)
    print(f"  Found {len(csv_map)} CSVs in the data folder")
    tids = transects["transect_id"].astype(str).tolist()

    records = []
    n_found = 0
    n_missing = 0
    for i, tid in enumerate(tids):
        if (i + 1) % 100 == 0 or i == len(tids) - 1:
            print(f"    Reading CSVs: {i+1}/{len(tids)}", end="\r")
        csv_path = csv_map.get(tid)
        if csv_path is None:
            n_missing += 1
            continue
        try:
            df = pd.read_csv(csv_path)
        except Exception:
            n_missing += 1
            continue
        n_found += 1

        date_col, chainage_col = _detect_coastsat_columns(df)
        if chainage_col is None:
            continue
        df["_dt"] = pd.to_datetime(df[date_col], errors="coerce",
                                     utc=True).dt.tz_localize(None)
        df = df.dropna(subset=["_dt", chainage_col])
        for _, r in df.iterrows():
            records.append({
                "transect_id": tid,
                "date": r["_dt"],
                "chainage_m": float(r[chainage_col]),
                "source": "coastsat",
            })
    print()
    print(f"  Read {n_found} CSVs, {n_missing} missing")
    out = pd.DataFrame(records) if records else pd.DataFrame(
        columns=["transect_id", "date", "chainage_m", "source"])
    print(f"  → {len(out)} chainage observations")
    if len(out) > 0:
        print(f"  → date range: {out['date'].min().date()} → "
              f"{out['date'].max().date()}")
    return out


# ============================================================
# CHAINAGE EXTRACTION — geometric intersection
# ============================================================

def extract_chainage_by_intersection(shorelines: gpd.GeoDataFrame,
                                      transects: gpd.GeoDataFrame,
                                      source_key: str,
                                      transect_length_m: float = None,
                                      valid_range_m: tuple = None) -> pd.DataFrame:
    """For each (shoreline feature, transect) pair whose bounding boxes
    overlap, compute the intersection point and the chainage (signed
    distance from transect origin along transect direction).

    transect_length_m / valid_range_m default to
    TRANSECT_INTERSECTION_LENGTH_M / ±CHAINAGE_MAX_ABS_M (CoastSat's
    short, beach-hugging transects).

    Returns long-format DataFrame:
        transect_id, date, chainage_m, source
    """
    if transect_length_m is None:
        transect_length_m = TRANSECT_INTERSECTION_LENGTH_M
    if valid_range_m is None:
        valid_range_m = (-CHAINAGE_MAX_ABS_M, CHAINAGE_MAX_ABS_M)

    print(f"\n[{SOURCE_LABELS[source_key]} chainage] extracting via intersection...")
    print(f"  {len(shorelines)} shoreline features × {len(transects)} transects, "
          f"transect length={transect_length_m:.0f} m, "
          f"valid chainage range={valid_range_m}")

    if len(shorelines) == 0 or len(transects) == 0:
        return pd.DataFrame(columns=["transect_id", "date", "chainage_m", "source"])

    # Build transect LineStrings + STRtree for fast bbox pre-filter
    transect_lines = []
    transect_meta = []
    for _, row in transects.iterrows():
        ox, oy = row["origin_x"], row["origin_y"]
        dx, dy = row["dir_x"], row["dir_y"]
        end_x = ox + transect_length_m * dx
        end_y = oy + transect_length_m * dy
        line = LineString([(ox, oy), (end_x, end_y)])
        transect_lines.append(line)
        transect_meta.append({
            "transect_id": row["transect_id"],
            "origin_x":     ox,
            "origin_y":     oy,
            "dir_x":        dx,
            "dir_y":        dy,
        })

    tree = STRtree(transect_lines)
    # Map from id(line) → index into transect_meta
    line_to_idx = {id(ln): i for i, ln in enumerate(transect_lines)}

    records = []
    n_dropped_range = 0
    n_shorelines = len(shorelines)
    for i, (_, sh_row) in enumerate(shorelines.iterrows()):
        if (i + 1) % 20 == 0 or i == n_shorelines - 1:
            print(f"    Processing shoreline {i+1}/{n_shorelines}", end="\r")
        shore_geom = sh_row.geometry
        shore_date = sh_row["date"]

        # STRtree.query returns candidate line geometries
        candidates = tree.query(shore_geom)
        # Handle both API styles (indices in newer shapely, geoms in older)
        for cand in candidates:
            if isinstance(cand, (int, np.integer)):
                idx = int(cand)
                t_line = transect_lines[idx]
            else:
                idx = line_to_idx.get(id(cand))
                t_line = cand
                if idx is None:
                    continue

            try:
                intersect = shore_geom.intersection(t_line)
            except Exception:
                continue
            if intersect.is_empty:
                continue

            meta = transect_meta[idx]
            pt = _choose_intersection_point(intersect, meta)
            if pt is None:
                continue
            # Chainage = signed dot product along direction
            chainage = ((pt.x - meta["origin_x"]) * meta["dir_x"]
                        + (pt.y - meta["origin_y"]) * meta["dir_y"])
            if not (valid_range_m[0] <= chainage <= valid_range_m[1]):
                n_dropped_range += 1
                continue
            records.append({
                "transect_id": meta["transect_id"],
                "date":         shore_date,
                "chainage_m":   float(chainage),
                "source":       source_key,
            })
    print()
    print(f"  {n_dropped_range} intersections dropped as outside valid range")

    result = pd.DataFrame(records)
    if len(result) > 0:
        pct = np.percentile(result["chainage_m"].values, [1, 25, 50, 75, 99])
        print(f"  Chainage percentiles [1,25,50,75,99]: "
              f"{np.round(pct, 1).tolist()}")
        # If multiple features from same date intersect same transect
        # (e.g., NC state segmented shorelines), collapse to per-date mean
        result = (result.groupby(["transect_id", "date", "source"])
                          .agg(chainage_m=("chainage_m", "mean"))
                          .reset_index())
    print(f"  → {len(result)} per-(transect, date) chainage observations")
    return result


def _choose_intersection_point(intersect, meta):
    """A single shapely intersection can be a Point, MultiPoint, LineString,
    or MultiLineString. Pick the single 2D coordinate best interpreted as
    the transect–shoreline crossing:
      - Point: use directly
      - MultiPoint: pick the one nearest to transect origin
      - LineString: pick midpoint (shoreline coincident with transect)
    """
    origin = Point(meta["origin_x"], meta["origin_y"])
    if isinstance(intersect, Point):
        return intersect
    if isinstance(intersect, MultiPoint):
        pts = list(intersect.geoms)
        pts.sort(key=lambda p: origin.distance(p))
        return pts[0]
    if isinstance(intersect, LineString):
        return intersect.interpolate(0.5, normalized=True)
    if isinstance(intersect, MultiLineString):
        best = None
        best_dist = np.inf
        for ln in intersect.geoms:
            mid = ln.interpolate(0.5, normalized=True)
            d = origin.distance(mid)
            if d < best_dist:
                best_dist = d
                best = mid
        return best
    # GeometryCollection or other — try first Point-like element
    if hasattr(intersect, "geoms"):
        for g in intersect.geoms:
            if isinstance(g, Point):
                return g
        for g in intersect.geoms:
            if isinstance(g, (LineString, MultiLineString)):
                return _choose_intersection_point(g, meta)
    return None


# ============================================================
# UNIFIED CHAINAGE TABLE
# ============================================================

def build_unified_chainage_table(
    coastsat_chainage: pd.DataFrame,
    wet_dry_chainage:  pd.DataFrame,
    nc_state_chainage: pd.DataFrame,
    transects: gpd.GeoDataFrame,
) -> pd.DataFrame:
    """Concatenate all three sources and enrich each row with alongshore
    position, domain, and decimal year. Also drops transects with too few
    observations to be useful.
    """
    all_chainage = pd.concat(
        [coastsat_chainage, wet_dry_chainage, nc_state_chainage],
        ignore_index=True,
    )
    if len(all_chainage) == 0:
        return all_chainage

    # Ensure transect_id is string on both sides
    all_chainage["transect_id"] = all_chainage["transect_id"].astype(str)
    tx = transects[["transect_id", "alongshore_m", "domain",
                    "origin_x", "origin_y"]].copy()
    tx["transect_id"] = tx["transect_id"].astype(str)

    merged = all_chainage.merge(tx, on="transect_id", how="inner")
    merged["decimal_year"] = decimal_year(merged["date"])
    merged["year"] = merged["date"].dt.year

    # Drop transects with too few observations
    counts = merged.groupby("transect_id").size()
    good = counts[counts >= MIN_OBSERVATIONS_PER_TRANSECT].index
    n_before = merged["transect_id"].nunique()
    merged = merged[merged["transect_id"].isin(good)].copy()
    n_after = merged["transect_id"].nunique()
    print(f"\n[merged] {len(merged)} total obs, "
          f"{n_after} transects kept "
          f"(dropped {n_before - n_after} with <{MIN_OBSERVATIONS_PER_TRANSECT} obs)")
    return merged


# ============================================================
# ANALYSIS — linear regression building block
# ============================================================

def regress_lrr(dec_years: np.ndarray, chainages: np.ndarray) -> dict:
    """Simple unweighted OLS linear regression.
    Returns dict with slope (m/yr), intercept, r_squared, n, se_slope."""
    x = np.asarray(dec_years, dtype=float)
    y = np.asarray(chainages, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    n = len(x)
    if n < 2:
        return dict(slope=np.nan, intercept=np.nan, r_squared=np.nan,
                    n=n, se_slope=np.nan)
    mx = x.mean()
    my = y.mean()
    dx = x - mx
    dy = y - my
    denom = (dx * dx).sum()
    if denom == 0:
        return dict(slope=np.nan, intercept=np.nan, r_squared=np.nan,
                    n=n, se_slope=np.nan)
    slope = (dx * dy).sum() / denom
    intercept = my - slope * mx
    y_pred = slope * x + intercept
    ss_res = ((y - y_pred) ** 2).sum()
    ss_tot = ((y - my) ** 2).sum()
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    if n > 2 and ss_res > 0:
        residual_var = ss_res / (n - 2)
        se_slope = np.sqrt(residual_var / denom)
    else:
        se_slope = np.nan
    return dict(slope=slope, intercept=intercept, r_squared=r_squared,
                n=n, se_slope=se_slope)


# ============================================================
# ANALYSIS — pre-installation LRR
# ============================================================

def compute_preinstall_lrr(chainage: pd.DataFrame,
                            transects: gpd.GeoDataFrame) -> pd.DataFrame:
    """For each transect, fit LRR using ALL observations dated before
    PRE_INSTALLATION_YEAR_CUTOFF (not a fixed list of years), pooled
    across every source that has pre-1970 data. Nourishment-affected
    observations (1966 Buxton fill, and the 1967 shoreline it biases)
    are NOT dropped -- they're included in the regression and simply
    FLAGGED (nourishment_zone=True) so a reader can see exactly which
    points are nourishment-influenced rather than having them silently
    removed.

    Returns DataFrame:
        transect_id, alongshore_m, domain, slope, intercept, r_squared,
        n, se_slope, nourishment_zone
    """
    print(f"\n{'='*72}\nPre-installation LRR")
    print(f"  Cutoff: all shorelines dated before {PRE_INSTALLATION_YEAR_CUTOFF}")

    # Filter to everything before the cutoff -- this picks up whatever
    # pre-1970 years actually exist in the data, rather than a
    # hand-picked subset.
    sub = chainage[chainage["year"] < PRE_INSTALLATION_YEAR_CUTOFF].copy()
    preinstall_years = sorted(sub["year"].unique().tolist())
    print(f"  Years present in data: {preinstall_years}")
    print(f"  → {len(sub)} observations across {sub['transect_id'].nunique()} "
          f"transects")
    print(f"  Per-source breakdown of pre-install observations:")
    for src, cnt in sub["source"].value_counts().items():
        print(f"    {SOURCE_LABELS.get(src, src):24s} {cnt}")

    # Flag (do NOT drop) observations in a documented nourishment zone,
    # so downstream plots can mark them distinctly instead of hiding them.
    def in_any_nourishment_zone(y):
        for name, year, y_min, y_max in NOURISHMENT_EXCLUSIONS:
            if year < PRE_INSTALLATION_YEAR_CUTOFF and y_min <= y <= y_max:
                return True
        return False

    for name, year, y_min, y_max in NOURISHMENT_EXCLUSIONS:
        if year >= PRE_INSTALLATION_YEAR_CUTOFF:
            continue
        n_flagged = int(((sub["year"] == year) &
                          (sub["origin_y"] >= y_min) &
                          (sub["origin_y"] <= y_max)).sum())
        print(f"  '{name}': {n_flagged} observations FLAGGED "
              f"(kept in the regression, marked nourishment_zone=True)")

    # Regress per transect (all data retained)
    tx_lookup = transects.set_index("transect_id")[["alongshore_m", "domain",
                                                      "origin_y",
                                                      "dist_from_groin_m"]].copy()
    tx_lookup.index = tx_lookup.index.astype(str)

    records = []
    for tid, grp in sub.groupby("transect_id"):
        fit = regress_lrr(grp["decimal_year"].values,
                          grp["chainage_m"].values)
        try:
            tx_row = tx_lookup.loc[str(tid)]
        except KeyError:
            continue
        records.append({
            "transect_id": tid,
            "alongshore_m": tx_row["alongshore_m"],
            "domain":       tx_row["domain"],
            "origin_y":     tx_row["origin_y"],
            "dist_from_groin_m": tx_row["dist_from_groin_m"],
            "n_obs":        fit["n"],
            "slope_m_yr":   fit["slope"],
            "intercept":    fit["intercept"],
            "r_squared":    fit["r_squared"],
            "se_slope":     fit["se_slope"],
            "nourishment_zone": in_any_nourishment_zone(tx_row["origin_y"]),
        })
    out = pd.DataFrame(records).sort_values("alongshore_m").reset_index(drop=True)
    n_good = int((out["n_obs"] >= 3).sum())
    print(f"  → {n_good}/{len(out)} transects have ≥3 pre-install observations "
          f"({int(out['nourishment_zone'].sum())} flagged nourishment-zone "
          f"transects included among them)")

    # ─── Noise-floor diagnostic, to help justify SIGNAL_ANOMALY_THRESHOLD_M_YR ───
    # Far-field (beyond SIGNAL_MAX_SEARCH_DISTANCE_M, i.e. outside where
    # the groin could plausibly reach) and not nourishment-flagged: this
    # is "ordinary" pre-install shoreline variability, unrelated to the
    # groin. Compare SIGNAL_ANOMALY_THRESHOLD_M_YR against this.
    control = out[(out["n_obs"] >= 3) & (~out["nourishment_zone"]) &
                  (out["dist_from_groin_m"].abs() > SIGNAL_MAX_SEARCH_DISTANCE_M)]
    if len(control) >= 5:
        noise_sd = float(control["slope_m_yr"].std())
        print(f"\n  [noise floor] far-field pre-install slope std dev: "
              f"{noise_sd:.2f} m/yr  (n={len(control)} transects, "
              f"beyond {SIGNAL_MAX_SEARCH_DISTANCE_M/1000:.0f} km)")
        print(f"  → SIGNAL_ANOMALY_THRESHOLD_M_YR is currently "
              f"{SIGNAL_ANOMALY_THRESHOLD_M_YR:.2f} m/yr = "
              f"{SIGNAL_ANOMALY_THRESHOLD_M_YR / noise_sd:.2f}x that noise floor "
              f"(a common rule of thumb: use 1-2x the far-field std dev, "
              f"so ordinary variability isn't flagged as groin signal)")
    else:
        print(f"\n  [noise floor] too few far-field control transects "
              f"({len(control)}) to estimate -- widen the check or use a "
              f"different control definition")

    out.attrs["preinstall_years"] = preinstall_years
    return out


def compute_post_install_lrr(chainage: pd.DataFrame,
                              transects: gpd.GeoDataFrame) -> pd.DataFrame:
    """For each transect, fit LRR using ALL observations dated on/after
    GROIN_INSTALLATION_YEAR (1970-present) -- the direct complement to
    compute_preinstall_lrr(). This isolates the post-groin rate
    cleanly, which compute_fullperiod_lrr() does NOT do: that function
    blends pre- and post-install years into a single slope per
    transect, which dilutes the apparent post-groin signal whenever a
    transect has pre-install observations mixed in with its post-
    install ones. Comparing THIS to the pre-install baseline is a
    cleaner before/after contrast than comparing the full-period line
    to the pre-install baseline.

    Returns DataFrame:
        transect_id, alongshore_m, domain, dist_from_groin_m, n_obs,
        slope_m_yr, intercept, r_squared, se_slope
    """
    print(f"\n{'='*72}\nPost-installation LRR")
    print(f"  Cutoff: all shorelines dated {GROIN_INSTALLATION_YEAR} or later")

    sub = chainage[chainage["year"] >= GROIN_INSTALLATION_YEAR].copy()
    print(f"  → {len(sub)} observations across {sub['transect_id'].nunique()} "
          f"transects")

    tx_lookup = transects.set_index("transect_id")[["alongshore_m", "domain",
                                                      "dist_from_groin_m"]].copy()
    tx_lookup.index = tx_lookup.index.astype(str)

    records = []
    for tid, grp in sub.groupby("transect_id"):
        fit = regress_lrr(grp["decimal_year"].values, grp["chainage_m"].values)
        try:
            tx_row = tx_lookup.loc[str(tid)]
        except KeyError:
            continue
        records.append({
            "transect_id": tid,
            "alongshore_m": tx_row["alongshore_m"],
            "domain":       tx_row["domain"],
            "dist_from_groin_m": tx_row["dist_from_groin_m"],
            "n_obs":        fit["n"],
            "slope_m_yr":   fit["slope"],
            "intercept":    fit["intercept"],
            "r_squared":    fit["r_squared"],
            "se_slope":     fit["se_slope"],
        })
    out = pd.DataFrame(records).sort_values("alongshore_m").reset_index(drop=True)
    n_good = int((out["n_obs"] >= 3).sum())
    print(f"  → {n_good}/{len(out)} transects have ≥3 post-install observations")
    return out


# ============================================================
# ANALYSIS — full-period LRR + piecewise breakpoint
# ============================================================

def compute_fullperiod_lrr(chainage: pd.DataFrame,
                            transects: gpd.GeoDataFrame) -> pd.DataFrame:
    """Per transect: single-slope OLS over ALL observations, and a
    simple piecewise-linear (two-segment) breakpoint search across
    BREAKPOINT_SEARCH_START..END. Returns both fits.
    """
    print(f"\n{'='*72}\nFull-period LRR + piecewise breakpoint")

    tx_lookup = transects.set_index("transect_id")[["alongshore_m", "domain",
                                                      "origin_y",
                                                      "dist_from_groin_m"]].copy()
    tx_lookup.index = tx_lookup.index.astype(str)

    records = []
    all_tids = sorted(chainage["transect_id"].unique())
    for i, tid in enumerate(all_tids):
        if (i + 1) % 100 == 0 or i == len(all_tids) - 1:
            print(f"    Regressing transect {i+1}/{len(all_tids)}", end="\r")
        grp = chainage[chainage["transect_id"] == tid]
        if len(grp) < 4:
            continue
        # Single-slope fit
        full = regress_lrr(grp["decimal_year"].values,
                           grp["chainage_m"].values)
        # Piecewise fit
        pw = _piecewise_breakpoint(grp["decimal_year"].values,
                                    grp["chainage_m"].values)
        try:
            tx_row = tx_lookup.loc[str(tid)]
        except KeyError:
            continue
        records.append({
            "transect_id": tid,
            "alongshore_m": tx_row["alongshore_m"],
            "domain":       tx_row["domain"],
            "origin_y":     tx_row["origin_y"],
            "dist_from_groin_m": tx_row["dist_from_groin_m"],
            "n_obs":        full["n"],
            "slope_full":   full["slope"],
            "r2_full":      full["r_squared"],
            "breakpoint_year":    pw["breakpoint_year"],
            "slope_before":       pw["slope_before"],
            "slope_after":        pw["slope_after"],
            "rss_full":     pw["rss_full"],
            "rss_piecewise": pw["rss_piecewise"],
            "rss_reduction_frac": pw["rss_reduction_frac"],
        })
    print()
    return pd.DataFrame(records).sort_values("alongshore_m").reset_index(drop=True)


def _piecewise_breakpoint(dec_years: np.ndarray, chainages: np.ndarray) -> dict:
    """Search for the two-segment breakpoint year (integer) that
    minimizes total residual sum of squares (RSS). Both segments must
    have BREAKPOINT_MIN_POINTS_PER_SEGMENT observations to be valid.
    """
    x = np.asarray(dec_years, dtype=float)
    y = np.asarray(chainages, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) < 2 * BREAKPOINT_MIN_POINTS_PER_SEGMENT:
        return dict(breakpoint_year=np.nan, slope_before=np.nan,
                    slope_after=np.nan, rss_full=np.nan,
                    rss_piecewise=np.nan, rss_reduction_frac=np.nan)
    # Single-slope RSS as baseline
    full_fit = regress_lrr(x, y)
    y_full_pred = full_fit["slope"] * x + full_fit["intercept"]
    rss_full = float(((y - y_full_pred) ** 2).sum())

    best_rss = np.inf
    best = dict(breakpoint_year=np.nan, slope_before=np.nan,
                slope_after=np.nan)
    for bp in range(BREAKPOINT_SEARCH_START, BREAKPOINT_SEARCH_END + 1):
        m_before = x < bp
        m_after  = x >= bp
        if m_before.sum() < BREAKPOINT_MIN_POINTS_PER_SEGMENT:
            continue
        if m_after.sum() < BREAKPOINT_MIN_POINTS_PER_SEGMENT:
            continue
        fit_b = regress_lrr(x[m_before], y[m_before])
        fit_a = regress_lrr(x[m_after],  y[m_after])
        if not (np.isfinite(fit_b["slope"]) and np.isfinite(fit_a["slope"])):
            continue
        pred_b = fit_b["slope"] * x[m_before] + fit_b["intercept"]
        pred_a = fit_a["slope"] * x[m_after]  + fit_a["intercept"]
        rss = float(((y[m_before] - pred_b) ** 2).sum()
                    + ((y[m_after]  - pred_a) ** 2).sum())
        if rss < best_rss:
            best_rss = rss
            best = dict(breakpoint_year=bp,
                        slope_before=fit_b["slope"],
                        slope_after=fit_a["slope"])
    if not np.isfinite(best_rss):
        return dict(breakpoint_year=np.nan, slope_before=np.nan,
                    slope_after=np.nan, rss_full=rss_full,
                    rss_piecewise=np.nan, rss_reduction_frac=np.nan)
    frac_reduction = (rss_full - best_rss) / rss_full if rss_full > 0 else np.nan
    return dict(
        breakpoint_year=best["breakpoint_year"],
        slope_before=best["slope_before"],
        slope_after=best["slope_after"],
        rss_full=rss_full,
        rss_piecewise=best_rss,
        rss_reduction_frac=frac_reduction,
    )


# ============================================================
# ANALYSIS — decadal LRR (fixed, non-overlapping 10-yr bins)
# ============================================================

def compute_decadal_lrr(chainage: pd.DataFrame) -> pd.DataFrame:
    """For each transect, fit LRR within each fixed, NON-OVERLAPPING
    DECADE_LENGTH_YEARS-year bin starting at DECADE_START_YEAR (default:
    1960s, 1970s, 1980s, ...) -- replaces a continuously-slid 1-year-
    step rolling window. The original version reported one number per year
    from 1975-2020, but adjacent years shared ~90% of their underlying
    observations (a 10-yr window stepped by 1 yr overlaps 9 of its 10
    years with its neighbor), so the resulting series looked smoother
    and more continuous than the actual survey frequency supports, and
    it was hard to say what was genuinely new information between e.g.
    center-year 1994 and 1995. Fixed decades are simpler to read and
    report: each number is a distinct, non-overlapping slice of time.

    Returns long-format DataFrame:
        transect_id, decade_start, decade_label, n_obs, slope_m_yr, r_squared
    """
    print(f"\n{'='*72}\nDecadal LRR ({DECADE_LENGTH_YEARS}-yr non-overlapping bins)")
    last_year = int(np.ceil(chainage["year"].max())) if len(chainage) else DECADE_START_YEAR
    decade_starts = list(range(DECADE_START_YEAR, last_year + 1, DECADE_LENGTH_YEARS))
    print(f"  Decades: {decade_starts[0]}s … {decade_starts[-1]}s "
          f"({DECADE_LENGTH_YEARS} yr each, non-overlapping), "
          f"min obs {DECADE_MIN_OBSERVATIONS}")

    records = []
    all_tids = sorted(chainage["transect_id"].unique())
    for i, tid in enumerate(all_tids):
        if (i + 1) % 100 == 0 or i == len(all_tids) - 1:
            print(f"    Transect {i+1}/{len(all_tids)}", end="\r")
        grp = chainage[chainage["transect_id"] == tid]
        dy = grp["decimal_year"].values
        yr = grp["year"].values
        ch = grp["chainage_m"].values
        for d0 in decade_starts:
            d1 = d0 + DECADE_LENGTH_YEARS
            mask = (yr >= d0) & (yr < d1)
            if mask.sum() < DECADE_MIN_OBSERVATIONS:
                continue
            fit = regress_lrr(dy[mask], ch[mask])
            records.append({
                "transect_id":  tid,
                "decade_start":  d0,          # decade START year
                "decade_label": f"{d0}s",
                "n_obs":        fit["n"],
                "slope_m_yr":   fit["slope"],
                "r_squared":    fit["r_squared"],
            })
    print()
    result = pd.DataFrame(records)
    print(f"  → {len(result)} (transect, decade) fits produced across "
          f"{result['decade_start'].nunique() if len(result) else 0} decades")
    return result


# ============================================================
# ANALYSIS — per-era LRR (one rate per transect per era)
# ============================================================

def compute_era_lrrs(chainage: pd.DataFrame,
                      transects: gpd.GeoDataFrame,
                      extra_periods: list = None) -> pd.DataFrame:
    """For each (transect × era in ERAS), fit an OLS linear regression
    of chainage vs. decimal year using ONLY observations within that
    era's year range. This is the structured complement to piecewise
    detection: instead of asking "when did the rate change", it asks
    "what was the rate during each pre-defined era?"

    extra_periods: optional list of additional (name, y_lo, y_hi)
    tuples computed the SAME way but NOT part of the official ERAS
    list -- e.g. a narrower display-only sub-period (see
    PRE_NOURISHMENT_PERIOD) shown for visual comparison in one plot,
    without disturbing the anomaly/signal-extent pipeline built on the
    3 official eras.

    Returns long-format DataFrame:
        transect_id, alongshore_m, dist_from_groin_m, domain, era,
        y_lo, y_hi, n_obs, slope_m_yr, intercept, r_squared, se_slope
    Rows are omitted where n_obs < 3 (regression not meaningful).
    """
    print(f"\n{'='*72}\nPer-era LRR (one rate per transect per era)")
    tx = transects.set_index("transect_id")[["alongshore_m", "domain",
                                                "origin_y",
                                                "dist_from_groin_m"]]
    tx.index = tx.index.astype(str)

    periods_to_compute = list(ERAS) + list(extra_periods or [])

    records = []
    all_tids = sorted(chainage["transect_id"].unique())
    for i, tid in enumerate(all_tids):
        if (i + 1) % 200 == 0 or i == len(all_tids) - 1:
            print(f"    Transect {i+1}/{len(all_tids)}", end="\r")
        grp = chainage[chainage["transect_id"] == tid]
        try:
            tx_row = tx.loc[str(tid)]
        except KeyError:
            continue
        for era_name, y_lo, y_hi in periods_to_compute:
            sub = grp[(grp["year"] >= y_lo) & (grp["year"] <= y_hi)]
            if len(sub) < 3:
                continue
            fit = regress_lrr(sub["decimal_year"].values,
                               sub["chainage_m"].values)
            records.append({
                "transect_id":       tid,
                "alongshore_m":      tx_row["alongshore_m"],
                "dist_from_groin_m": tx_row["dist_from_groin_m"],
                "domain":            tx_row["domain"],
                "origin_y":          tx_row["origin_y"],
                "era":               era_name,
                "y_lo":              y_lo,
                "y_hi":              y_hi,
                "n_obs":             fit["n"],
                "slope_m_yr":        fit["slope"],
                "intercept":         fit["intercept"],
                "r_squared":         fit["r_squared"],
                "se_slope":          fit["se_slope"],
            })
    print()
    out = pd.DataFrame(records)
    if len(out) > 0:
        summary = out.groupby("era")["transect_id"].nunique()
        for era_name, n_tx in summary.items():
            n_total = len(out[out["era"] == era_name])
            print(f"  {era_name:>15s}: {n_tx} transects, {n_total} fits")
    return out




def interp_preinstall_baseline(preinstall_lrr: pd.DataFrame) -> callable:
    """Build the regional pre-install baseline by linearly interpolating
    ("connecting the dots") between the sorted per-transect pre-install
    LRR values themselves -- no binning, no smoothing. Query distances
    beyond the data's range are clamped to the nearest end value.

    Used for BOTH the anomaly calculation and the reference curve drawn
    on the profile plots, so there's exactly one baseline definition,
    not two. No LOWESS/kernel smoother is used anywhere in this script.
    (An earlier version used a 500 m binned-median step function
    instead of interpolating; the sharp box-like jumps at each bin
    edge read poorly on the profile plots, so this connects the actual
    per-transect points directly instead.)

    Returns a function baseline(dist_from_groin_m) -> float (vectorized).
    """
    good = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                           preinstall_lrr["slope_m_yr"].notna() &
                           (~preinstall_lrr["nourishment_zone"])].copy()
    if len(good) < 2:
        print(f"  ! Only {len(good)} transects with usable pre-install LRR; "
              f"baseline will be a flat constant (median)")
        median = good["slope_m_yr"].median() if len(good) else 0.0
        return lambda x: np.full_like(np.asarray(x, dtype=float),
                                        float(median))

    good = good.sort_values("dist_from_groin_m")
    xs = good["dist_from_groin_m"].values
    ys = good["slope_m_yr"].values

    def baseline(query_dist_from_groin_m):
        q = np.asarray(query_dist_from_groin_m, dtype=float)
        return np.interp(q, xs, ys, left=ys[0], right=ys[-1])

    print(f"  Built connect-the-dots baseline (linear interpolation, "
          f"NOT binned/smoothed) over {len(good)} pre-install transects, "
          f"median = {np.median(ys):+.2f} m/yr")
    return baseline


# ============================================================
# ANALYSIS — decadal anomaly relative to baseline
# ============================================================

def compute_decadal_anomaly(decadal_lrr: pd.DataFrame,
                             transects: gpd.GeoDataFrame,
                             baseline_fn) -> pd.DataFrame:
    """Enrich the decadal LRR table with distance-from-groin and the
    anomaly relative to the pre-install regional baseline at each
    transect. Anomaly = LRR - baseline. Positive anomaly means
    the transect is accreting (or eroding less) faster than baseline;
    negative means the opposite.
    """
    print(f"\n{'='*72}\nDecadal anomaly relative to baseline")
    tx = transects[["transect_id", "alongshore_m", "dist_from_groin_m",
                     "origin_y", "domain"]].copy()
    tx["transect_id"] = tx["transect_id"].astype(str)
    merged = decadal_lrr.merge(tx, on="transect_id", how="inner")
    merged["baseline_lrr"] = baseline_fn(merged["dist_from_groin_m"].values)
    merged["anomaly_m_yr"] = merged["slope_m_yr"] - merged["baseline_lrr"]
    print(f"  {len(merged)} decadal anomaly rows across "
          f"{merged['transect_id'].nunique()} transects and "
          f"{merged['decade_start'].nunique()} decades")
    return merged


# ============================================================
# ANALYSIS — signal extent over time
# ============================================================

def _contiguous_extent_one_side(dist_abs: np.ndarray, anomaly: np.ndarray,
                                 thresh: float, bin_w: float, max_d: float,
                                 max_gap_bins: int, direction: int) -> tuple:
    """Walk outward from the groin (distance 0) in bin_w-wide bins,
    taking the median anomaly per bin. The signal is considered to
    extend through the bin as long as direction*median > thresh,
    tolerating up to max_gap_bins consecutive non-exceeding bins before
    stopping (absorbs single-bin noise). direction=+1 looks for
    anomaly > +thresh (updrift/accretion signal); direction=-1 looks
    for anomaly < -thresh (downdrift/erosion signal).
    Returns (extent_m, peak_anomaly_within_extent).
    """
    n_bins = int(np.ceil(max_d / bin_w))
    if len(dist_abs) == 0:
        return 0.0, np.nan
    bin_idx = np.clip((dist_abs // bin_w).astype(int), 0, n_bins - 1)
    bin_medians = np.full(n_bins, np.nan)
    for b in range(n_bins):
        vals = anomaly[bin_idx == b]
        if len(vals) > 0:
            bin_medians[b] = np.median(vals)

    extent_bins = 0
    gap = 0
    for b in range(n_bins):
        v = bin_medians[b]
        exceeds = np.isfinite(v) and (direction * v > thresh)
        if exceeds:
            extent_bins = b + 1
            gap = 0
        else:
            gap += 1
            if gap > max_gap_bins:
                break
    extent_m = extent_bins * bin_w

    if extent_bins > 0:
        within = bin_medians[:extent_bins]
        within = within[np.isfinite(within)]
        peak = float(within.max() if direction > 0 else within.min()) \
               if len(within) else np.nan
    else:
        peak = np.nan
    return extent_m, peak


def compute_signal_extent_over_time(decadal_anomaly: pd.DataFrame) -> pd.DataFrame:
    """For each decade, find how far the groin
    signal extends updrift (positive anomaly = accretion beyond
    baseline) and downdrift (negative anomaly = erosion beyond
    baseline).

    Method: CONTIGUOUS run from the groin outward. Distances are
    binned (SIGNAL_EXTENT_BIN_WIDTH_M) and the run of bins whose
    median |anomaly| exceeds SIGNAL_ANOMALY_THRESHOLD_M_YR, starting
    at the groin, defines the extent (small gaps up to
    SIGNAL_EXTENT_MAX_GAP_BINS are tolerated). This is a deliberate
    change from "the single furthest transect anywhere that happens to
    cross the threshold" -- that definition let one noisy far-field
    transect get reported as the edge of the groin's influence, which
    produced a signal-extent time series with implausible sudden km-
    scale jumps unrelated to any physical process near the groin.

    Returns DataFrame:
        decade_start,
        updrift_extent_m,       (contiguous extent where anomaly > +threshold)
        updrift_peak_anomaly,   (peak anomaly WITHIN that extent)
        downdrift_extent_m,     (contiguous extent where anomaly < -threshold)
        downdrift_peak_anomaly,
        n_updrift_transects, n_downdrift_transects
    """
    print(f"\n{'='*72}\nSignal extent over time")
    print(f"  Threshold: |anomaly| ≥ {SIGNAL_ANOMALY_THRESHOLD_M_YR:.2f} m/yr, "
          f"contiguous from groin, bin={SIGNAL_EXTENT_BIN_WIDTH_M:.0f} m, "
          f"max gap={SIGNAL_EXTENT_MAX_GAP_BINS} bin(s)")
    print(f"  Max search distance: {SIGNAL_MAX_SEARCH_DISTANCE_M/1000:.1f} km "
          f"either side of groin")

    records = []
    thresh = SIGNAL_ANOMALY_THRESHOLD_M_YR
    max_d = SIGNAL_MAX_SEARCH_DISTANCE_M
    bin_w = SIGNAL_EXTENT_BIN_WIDTH_M
    max_gap = SIGNAL_EXTENT_MAX_GAP_BINS
    for c, grp in decadal_anomaly.groupby("decade_start"):
        grp = grp[grp["dist_from_groin_m"].abs() <= max_d]
        up = grp[grp["dist_from_groin_m"] > 0]
        dn = grp[grp["dist_from_groin_m"] < 0]

        updrift_extent, updrift_peak = _contiguous_extent_one_side(
            up["dist_from_groin_m"].values, up["anomaly_m_yr"].values,
            thresh, bin_w, max_d, max_gap, direction=+1)
        downdrift_extent, downdrift_peak = _contiguous_extent_one_side(
            dn["dist_from_groin_m"].abs().values, dn["anomaly_m_yr"].values,
            thresh, bin_w, max_d, max_gap, direction=-1)

        records.append({
            "decade_start":            int(c),
            "updrift_extent_m":       updrift_extent,
            "updrift_peak_anomaly":   updrift_peak,
            "downdrift_extent_m":     downdrift_extent,
            "downdrift_peak_anomaly": downdrift_peak,
            "n_updrift_transects":   int(len(up)),
            "n_downdrift_transects": int(len(dn)),
        })
    out = pd.DataFrame(records).sort_values("decade_start").reset_index(drop=True)
    return out


# ============================================================
# PLOTS
# ============================================================

def plot_alongshore_lrr_profile(preinstall_lrr: pd.DataFrame,
                                 fullperiod_lrr: pd.DataFrame,
                                 postinstall_lrr: pd.DataFrame,
                                 transects: gpd.GeoDataFrame,
                                 groin: dict,
                                 output_path: str):
    """LRR (m/yr) vs distance from groin. Every series is per-transect
    data connected in distance order -- no binning, no smoothing, no
    averaging -- just the real points joined so the trend reads as a
    line instead of a loose cloud of dots.
        - Grey: pre-install baseline, per-transect LRR
          (nourishment-flagged transects are included here too, not
          singled out with a different marker, though still flagged
          in the data)
        - Dark grey/black: full-period (1849-2024, every observation
          blended into one slope per transect)
        - Blue: post-install only (1970-2024) -- a cleaner
          before/after contrast against the pre-install baseline than
          the full-period points, since it isn't diluted by any
          pre-install years mixed into the same transect's regression
        - Vertical shaded band: groin footprint
        - Top axis: CASCADE domain number
    """
    print(f"\n[plot] alongshore LRR profile → {output_path}")

    fig, ax = plt.subplots(figsize=(14, 7))

    # ─── Full-period per-transect points + connecting line ────────────────
    # (still every real point, in order -- nothing averaged/binned)
    fs = fullperiod_lrr[fullperiod_lrr["slope_full"].notna()].copy()
    fs = fs.sort_values("dist_from_groin_m")
    if len(fs) > 0:
        ax.plot(fs["dist_from_groin_m"] / 1000, fs["slope_full"],
                 color="#333", lw=0.7, alpha=0.45, marker=".", ms=4,
                 markeredgewidth=0, zorder=2,
                 label=f"Full-period per-transect LRR  "
                        f"(n={len(fs)} transects)")

    # ─── Post-install (1970-2024) per-transect points + connecting line ───
    ps = postinstall_lrr[postinstall_lrr["slope_m_yr"].notna()].copy()
    ps = ps.sort_values("dist_from_groin_m")
    if len(ps) > 0:
        ax.plot(ps["dist_from_groin_m"] / 1000, ps["slope_m_yr"],
                 color="#1565C0", lw=0.7, alpha=0.45, marker=".", ms=4,
                 markeredgewidth=0, zorder=2,
                 label=f"Post-install ({GROIN_INSTALLATION_YEAR}-present) "
                        f"per-transect LRR  (n={len(ps)} transects)")

    # ─── Pre-install per-transect points + connecting line (nourishment- ───
    # flagged transects included, same marker as everything else -- still
    # kept in the regression and in preinstall_lrr's nourishment_zone
    # column for anyone who wants to filter on it, just not called out
    # with a different marker here)
    ok_pre = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                             preinstall_lrr["slope_m_yr"].notna()].sort_values(
        "dist_from_groin_m")
    preinstall_years = preinstall_lrr.attrs.get("preinstall_years", [])
    yr_range_str = (f"{min(preinstall_years)}–{max(preinstall_years)}"
                    if preinstall_years else "no years found")
    if len(ok_pre) > 0:
        ax.plot(ok_pre["dist_from_groin_m"] / 1000, ok_pre["slope_m_yr"],
                 color=ERA_COLORS["Pre-install"], lw=1.0, alpha=0.6,
                 marker=".", ms=6, markeredgewidth=0, zorder=5,
                 label=f"Pre-install baseline  (per-transect LRR, "
                        f"n={len(ok_pre)}, {yr_range_str}, "
                        f"years={preinstall_years})")

    # ─── Reference lines ──────────────────────────────────────────────────
    ax.axhline(0, color="#999", lw=0.5, linestyle=":")
    ax.axvline(0, color="black", lw=1.0, alpha=0.7,
                label=f"Groin  (updrift = "
                       f"{'north' if UPDRIFT_DIRECTION == 'north' else 'south'})")

    # Groin footprint band -- spans the real southernmost-to-northernmost
    # extent of the groin field, relative to the northernmost-groin
    # origin (0 km) -- generally NOT symmetric around 0, since 0 sits at
    # the north end of the field, not its middle.
    if groin.get("geometry") is not None:
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        g_a = sign * (groin["y_min"] - groin["reference_y"]) / 1000
        g_b = sign * (groin["y_max"] - groin["reference_y"]) / 1000
        ax.axvspan(min(g_a, g_b), max(g_a, g_b),
                    color="black", alpha=0.12, zorder=0,
                    label="Groin field  (first to last groin)")

    # Clip x-axis to start exactly at domain 1 (downdrift) and run to
    # PLOT_UPDRIFT_MAX_KM updrift, rather than wherever the data
    # happens to start/stop
    d1_km = domain_dist_km(transects, 1)
    apply_distance_xlim_and_autoscale_y(ax, x_max_km=PLOT_UPDRIFT_MAX_KM,
                                          x_min_km=d1_km)
    x_lo, x_hi = ax.get_xlim()
    set_domain_primary_axis(ax, transects)
    add_geographic_annotations(ax, transects)

    ax.set_ylabel("LRR  (m/yr, + = accretion)", fontsize=FONT_AXIS_LABEL)
    ax.set_title(
        "Shoreline Change Rate vs. Distance from Groin",
        fontsize=FONT_TITLE, fontweight="bold", pad=10,
    )
    ax.grid(True, alpha=0.25, linewidth=0.4)
    handles, labels = ax.get_legend_handles_labels()
    handles += annotation_legend_handles()
    ax.legend(handles=handles, loc="lower right", fontsize=FONT_LEGEND, framealpha=0.92)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


_LOWESS_WARNED = False   # print the "statsmodels missing" warning only once


def _lowess_overlay(x: np.ndarray, y: np.ndarray, frac: float):
    """LOWESS-smoothed (x, y), sorted by x -- a VISUAL AID overlay
    only, drawn as a lighter shade of a line's own color, alongside
    (never instead of) the raw per-transect points + connecting line.
    Returns (x_smooth, y_smooth), or (None, None) if statsmodels isn't
    available or there's too little data to smooth meaningfully.
    """
    global _LOWESS_WARNED
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess
    except ImportError:
        if not _LOWESS_WARNED:
            print("  ! statsmodels not installed (pip install statsmodels) "
                  "-- skipping smoothed overlay lines")
            _LOWESS_WARNED = True
        return None, None
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 10:
        return None, None
    order = np.argsort(x)
    smoothed = lowess(y[order], x[order], frac=frac, return_sorted=True)
    return smoothed[:, 0], smoothed[:, 1]


def plot_era_lrr_profile(era_lrrs: pd.DataFrame,
                          preinstall_lrr: pd.DataFrame,
                          transects: gpd.GeoDataFrame,
                          groin: dict,
                          output_path: str):
    """Alongshore LRR profile with one color per era. Every series is
    per-transect data connected in distance order -- no binning, no
    averaging -- just the real points joined so the trend reads as a
    line. A LOWESS-smoothed overlay (lighter shade of the same color,
    ERA_PROFILE_LOWESS_FRAC bandwidth, deliberately tighter than the
    shared SMOOTHED_OVERLAY_LOWESS_FRAC so it captures the sharp jumps
    right around the groin) is drawn alongside each line as a visual
    aid, never replacing the raw connected points.

    Includes PRE_NOURISHMENT_PERIOD (1996-2021, pre-2022 nourishment)
    as an extra display-only line alongside the full Deteriorated era,
    to check whether the 2022 nourishments are skewing that era's
    overall rate.
    """
    print(f"\n[plot] per-era LRR profile → {output_path}")
    if len(era_lrrs) == 0:
        print("  ! no era LRR data")
        return

    fig, ax = plt.subplots(figsize=(14, 7.5))

    # Pre-install baseline as reference
    ok_pre = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                             preinstall_lrr["slope_m_yr"].notna()].sort_values(
        "dist_from_groin_m")
    preinstall_years = preinstall_lrr.attrs.get("preinstall_years", [])
    yr_range_str = (f"{min(preinstall_years)}–{max(preinstall_years)}"
                    if preinstall_years else "no years found")
    pre_color = ERA_COLORS["Pre-install"]
    if len(ok_pre) > 0:
        ax.plot(ok_pre["dist_from_groin_m"] / 1000, ok_pre["slope_m_yr"],
                 color=pre_color, lw=0.6, alpha=0.35, marker=".", ms=4,
                 markeredgewidth=0, zorder=3,
                 label=f"Pre-install baseline  (per-transect LRR, "
                        f"n={len(ok_pre)}, {yr_range_str})")
        xs, ys = _lowess_overlay(ok_pre["dist_from_groin_m"].values / 1000,
                                   ok_pre["slope_m_yr"].values,
                                   ERA_PROFILE_LOWESS_FRAC)
        if xs is not None:
            ax.plot(xs, ys, color=pre_color, lw=2.8, zorder=7)

    # Post-install eras (official) + the pre-nourishment sub-period:
    # per-transect points + connecting line, both in one faint call
    # (light background texture) + smoothed overlay (bold, the primary
    # signal to read).
    periods_to_plot = [(n, lo, hi) for n, lo, hi in ERAS if n != "Pre-install"]
    periods_to_plot.append(PRE_NOURISHMENT_PERIOD)
    for era_name, y_lo, y_hi in periods_to_plot:
        sub = era_lrrs[era_lrrs["era"] == era_name].sort_values("dist_from_groin_m")
        if len(sub) < 5:
            continue
        is_pre_nourishment = (era_name == PRE_NOURISHMENT_PERIOD[0])
        color = PRE_NOURISHMENT_COLOR if is_pre_nourishment else ERA_COLORS.get(era_name, "#333")

        ax.plot(sub["dist_from_groin_m"] / 1000, sub["slope_m_yr"],
                 color=color, lw=0.6, alpha=0.35, marker=".", ms=3.5,
                 markeredgewidth=0, zorder=3,
                 label=f"{era_name}  ({y_lo}–{y_hi}, n={len(sub)})")
        xs, ys = _lowess_overlay(sub["dist_from_groin_m"].values / 1000,
                                   sub["slope_m_yr"].values,
                                   ERA_PROFILE_LOWESS_FRAC)
        if xs is not None:
            # Pre-nourishment overlaps Deteriorated's range and color
            # family -- dash it so the two smoothed lines don't blur
            # together where they nearly coincide.
            ax.plot(xs, ys, color=color, lw=2.8, zorder=6,
                     linestyle="--" if is_pre_nourishment else "-")

    # Reference lines
    ax.axhline(0, color="#999", lw=0.5, linestyle=":")
    ax.axvline(0, color="black", lw=1.0, alpha=0.7,
                label=f"Groin  (updrift = "
                       f"{'north' if UPDRIFT_DIRECTION == 'north' else 'south'})")
    if groin.get("geometry") is not None:
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        g_a = sign * (groin["y_min"] - groin["reference_y"]) / 1000
        g_b = sign * (groin["y_max"] - groin["reference_y"]) / 1000
        ax.axvspan(min(g_a, g_b), max(g_a, g_b),
                    color="black", alpha=0.12, zorder=0,
                    label="Groin field  (first to last groin)")

    d1_km = domain_dist_km(transects, 1)
    apply_distance_xlim_and_autoscale_y(ax, x_max_km=PLOT_UPDRIFT_MAX_KM,
                                          x_min_km=d1_km)
    x_lo, x_hi = ax.get_xlim()
    set_domain_primary_axis(ax, transects)
    add_geographic_annotations(ax, transects)
    ax.set_ylabel("LRR  (m/yr, + = accretion)", fontsize=FONT_AXIS_LABEL)
    ax.set_title(
        "Shoreline Change Rate by Structural Era",
        fontsize=FONT_TITLE, fontweight="bold", pad=10,
    )
    ax.grid(True, alpha=0.25, linewidth=0.4)
    handles, labels = ax.get_legend_handles_labels()
    handles += annotation_legend_handles()
    ax.legend(handles=handles, loc="best", fontsize=FONT_LEGEND, framealpha=0.92, ncol=1)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _build_decade_periods(start_year: int, increment_years: int,
                           end_year: int) -> list:
    """Build a list of (name, y_lo, y_hi) non-overlapping period tuples
    from start_year to end_year in increment_years steps -- e.g.
    1970-1980, 1980-1990, ... The LAST period may be shorter than
    increment_years if end_year isn't an exact multiple of steps away
    from start_year -- that's expected (the data just runs out there),
    not an error.
    """
    periods = []
    y = start_year
    while y < end_year:
        y_hi = min(y + increment_years - 1, end_year)
        periods.append((f"{y}-{y_hi}", y, y_hi))
        y += increment_years
    return periods


def plot_decade_lrr_profile(decade_period_lrrs: pd.DataFrame,
                             transects: gpd.GeoDataFrame,
                             groin: dict,
                             periods: list,
                             increment_years: int,
                             output_path: str):
    """How has the groin's effect changed continuously over its
    lifespan? Shows the LRR profile for each fixed-length post-install
    window in `periods` (see DECADE_PLOT_START_YEAR /
    DECADE_PLOT_INCREMENTS_YEARS) as its own line, colored along a
    sequential colormap (DECADE_PLOT_COLORMAP) so chronological order
    reads directly from color -- a finer-grained complement to the
    3-era profile (which only splits pre/functional/deteriorated),
    for the specific question of whether effectiveness has been
    trending rather than just stepping between two states.
    """
    print(f"\n[plot] decade-increment LRR profile ({increment_years}-yr) → {output_path}")
    fig, ax = plt.subplots(figsize=(14, 7.5))

    n = len(periods)
    cmap = plt.get_cmap(DECADE_PLOT_COLORMAP)
    # Sample from 0.45-1.0: both 0.15 and 0.28 were still reported too
    # light for the earliest period's line to read clearly against a
    # white background. 0.45 gives a solid, clearly-visible medium
    # green as the lightest color, still distinguishable from the
    # darkest (1.0) end.
    colors = [cmap(0.45 + 0.55 * (i / max(n - 1, 1))) for i in range(n)]

    any_plotted = False
    for i, (period_name, y_lo, y_hi) in enumerate(periods):
        sub = decade_period_lrrs[decade_period_lrrs["era"] == period_name] \
            .sort_values("dist_from_groin_m")
        if len(sub) < 5:
            print(f"  ! {period_name}: only {len(sub)} transects -- skipping line")
            continue
        color = colors[i]
        xs, ys = _lowess_overlay(sub["dist_from_groin_m"].values / 1000,
                                   sub["slope_m_yr"].values,
                                   SMOOTHED_OVERLAY_LOWESS_FRAC)
        if xs is not None:
            ax.plot(xs, ys, color=color, lw=1.5, zorder=6,
                     label=f"{period_name}  (n={len(sub)})")
        any_plotted = True

    if not any_plotted:
        print("  ! no period had enough data -- skipping plot")
        plt.close(fig)
        return

    ax.axhline(0, color="#999", lw=0.5, linestyle=":")
    ax.axvline(0, color="black", lw=1.0, alpha=0.7,
                label=f"Groin  (updrift = "
                       f"{'north' if UPDRIFT_DIRECTION == 'north' else 'south'})")
    if groin.get("geometry") is not None:
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        g_a = sign * (groin["y_min"] - groin["reference_y"]) / 1000
        g_b = sign * (groin["y_max"] - groin["reference_y"]) / 1000
        ax.axvspan(min(g_a, g_b), max(g_a, g_b),
                    color="black", alpha=0.12, zorder=0,
                    label="Groin field  (first to last groin)")

    d1_km = domain_dist_km(transects, 1)
    apply_distance_xlim_and_autoscale_y(ax, x_max_km=PLOT_UPDRIFT_MAX_KM,
                                          x_min_km=d1_km)
    x_lo, x_hi = ax.get_xlim()
    set_domain_primary_axis(ax, transects)
    add_geographic_annotations(ax, transects)
    ax.set_ylabel("LRR  (m/yr, + = accretion)", fontsize=FONT_AXIS_LABEL)
    ax.set_title(
        f"Shoreline Change Rate — {increment_years}-Year Windows\n"
        f"Color: light → dark = {periods[0][0].split('-')[0]} → present",
        fontsize=FONT_TITLE, fontweight="bold", pad=10,
    )
    ax.grid(True, alpha=0.25, linewidth=0.4)
    handles, labels = ax.get_legend_handles_labels()
    handles += annotation_legend_handles()
    ax.legend(handles=handles, loc="best", fontsize=FONT_LEGEND, framealpha=0.92, ncol=2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _zone_mean_and_label(ax, df: pd.DataFrame, zone_domains: tuple,
                          color: str, domain_to_dist):
    """Draw a heavy horizontal segment at a zone's mean rate, plus a
    numeric label above/below it -- adapted from
    HAT_groin_zone_investigation.py's draw_zone_mean_rates(), using
    domain-number zone boundaries (converted to km via domain_to_dist)
    instead of that script's alongshore-meters convention.
    """
    d_lo, d_hi = zone_domains
    km_lo, km_hi = float(domain_to_dist(d_lo)), float(domain_to_dist(d_hi))
    km_lo, km_hi = min(km_lo, km_hi), max(km_lo, km_hi)
    mask = ((df["dist_from_groin_m"] / 1000 >= km_lo) &
            (df["dist_from_groin_m"] / 1000 <= km_hi))
    if mask.sum() == 0:
        return
    zone_mean = df.loc[mask, "slope_m_yr"].mean()
    ax.hlines(zone_mean, km_lo, km_hi, colors=color, lw=2.5, alpha=0.85, zorder=6)
    ax.annotate(f"{zone_mean:+.1f}", xy=((km_lo + km_hi) / 2, zone_mean),
                 xytext=(0, 6 if zone_mean >= 0 else -12),
                 textcoords="raw_offset points", ha="center",
                 fontsize=8, color=color, fontweight="bold", zorder=7)


def plot_zone_panels(preinstall_lrr: pd.DataFrame, era_lrrs: pd.DataFrame,
                      transects: gpd.GeoDataFrame, groin: dict,
                      output_path: str,
                      extra_functional_panels: list = None):
    """Multi-panel alongshore profile, one stacked panel per era --
    adapted from HAT_groin_zone_investigation.py's
    alongshore_profile_panels.png, but using THIS script's own era
    definitions and data (Pre-install, Functional groin, Deteriorated
    -- already computed elsewhere, not recomputed here) and a
    simplified two-zone structure (ZONE_PANEL_DOWNDRIFT_DOMAINS,
    ZONE_PANEL_UPDRIFT_DOMAINS) instead of that script's four-zone
    system, since these are the specific zones asked for.

    extra_functional_panels: optional list of (name, df) tuples --
    e.g. "Functional groin (first 5yr)" -- inserted right after the
    full "Functional groin" panel, for a closer look at whether the
    groin's early years looked different from its full functional era.

    Each panel shows the per-transect points connected in order, a
    LOWESS-smoothed overlay (matching this script's own convention
    elsewhere, not the reference script's rolling mean), a faint
    pre-install reference line (on every panel except its own), and a
    mean-rate horizontal segment + label for each zone.
    """
    print(f"\n[plot] zone panels → {output_path}")
    domain_to_dist, _ = build_domain_to_dist_km(transects)

    ok_pre = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                             preinstall_lrr["slope_m_yr"].notna()]
    panels = [
        ("Pre-install", ok_pre),
        ("Functional groin", era_lrrs[era_lrrs["era"] == "Functional groin"]),
    ]
    panels.extend(extra_functional_panels or [])
    panels.append(("Deteriorated", era_lrrs[era_lrrs["era"] == "Deteriorated"]))
    n_panels = len(panels)

    d1_km = domain_dist_km(transects, 1)
    x_min_km = d1_km if d1_km is not None else -3.0
    x_max_km = float(domain_to_dist(ZONE_PANEL_X_MAX_DOMAIN))

    fig, axes = plt.subplots(n_panels, 1, figsize=(12, 3.4 * n_panels),
                               sharex=True)
    if n_panels == 1:
        axes = [axes]

    preinstall_df = panels[0][1].sort_values("dist_from_groin_m")

    zones = [
        ("Downdrift zone", ZONE_PANEL_DOWNDRIFT_DOMAINS, ZONE_PANEL_COLOR_DOWNDRIFT),
        ("Updrift zone", ZONE_PANEL_UPDRIFT_DOMAINS, ZONE_PANEL_COLOR_UPDRIFT),
    ]

    for i, (era_name, df) in enumerate(panels):
        ax = axes[i]
        if era_name in ERA_COLORS:
            color = ERA_COLORS[era_name]
        elif "Functional" in era_name:
            color = lighten_color(ERA_COLORS["Functional groin"], amount=0.35)
        else:
            color = "#333"
        if len(df) == 0:
            ax.text(0.5, 0.5, f"{era_name}: no data", transform=ax.transAxes,
                     ha="center", va="center")
            continue
        df_sorted = df.sort_values("dist_from_groin_m")

        for zname, zdomains, zcolor in zones:
            zd_lo, zd_hi = zdomains
            zkm_lo, zkm_hi = float(domain_to_dist(zd_lo)), float(domain_to_dist(zd_hi))
            zkm_lo, zkm_hi = min(zkm_lo, zkm_hi), max(zkm_lo, zkm_hi)
            ax.axvspan(zkm_lo, zkm_hi, color=zcolor, alpha=0.12, zorder=0,
                        label=zname if i == 0 else None)

        ax.axhline(0, color="k", lw=0.5, zorder=2)
        ax.axvline(0, color="black", lw=1.0, alpha=0.7, zorder=2,
                    label="Groin" if i == 0 else None)

        # Pre-install reference line (faint), on every panel except its own
        if era_name != "Pre-install" and len(preinstall_df) > 0:
            xs_pc, ys_pc = _lowess_overlay(
                preinstall_df["dist_from_groin_m"].values / 1000,
                preinstall_df["slope_m_yr"].values,
                ERA_PROFILE_LOWESS_FRAC)
            if xs_pc is not None:
                ax.plot(xs_pc, ys_pc, color=ERA_COLORS["Pre-install"], lw=1.2,
                         ls="--", alpha=0.7, zorder=3,
                         label="Pre-install (smoothed)" if i == 1 else None)

        ax.plot(df_sorted["dist_from_groin_m"] / 1000, df_sorted["slope_m_yr"],
                 marker="o", ms=2.5, lw=0.8, color=color, alpha=0.6, zorder=4)

        xs, ys = _lowess_overlay(df_sorted["dist_from_groin_m"].values / 1000,
                                   df_sorted["slope_m_yr"].values,
                                   ERA_PROFILE_LOWESS_FRAC)
        if xs is not None:
            ax.plot(xs, ys, color=color, lw=2.2, alpha=0.95, zorder=5,
                     label="Smoothed" if i == 0 else None)

        ax.set_xlim(x_min_km, x_max_km)
        ax.set_ylabel(f"{era_name}\n(m/yr)", fontsize=FONT_ANNOTATION)
        ax.grid(True, alpha=0.3)

        for _, zdomains, zcolor in zones:
            _zone_mean_and_label(ax, df_sorted, zdomains, zcolor, domain_to_dist)

        if i == 0:
            domain_ticks = np.arange(1, ZONE_PANEL_X_MAX_DOMAIN + 1, 2)
            tick_km = domain_to_dist(domain_ticks)
            secax = ax.secondary_xaxis("top")
            secax.set_xticks(tick_km)
            secax.set_xticklabels([f"D{int(d)}" for d in domain_ticks],
                                   fontsize=FONT_TICK - 2)

    axes[-1].set_xlabel("Distance from groin (km, alongshore)  [+ = updrift]",
                          fontsize=FONT_AXIS_LABEL)
    fig.suptitle("Alongshore Shoreline Change Rate by Era, with Analysis Zones",
                  fontsize=FONT_TITLE, fontweight="bold", y=1.0)

    handles, labels = axes[0].get_legend_handles_labels()
    if len(axes) > 1:
        h2, l2 = axes[1].get_legend_handles_labels()
        for h, l in zip(h2, l2):
            if l not in labels:
                handles.append(h)
                labels.append(l)
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.955),
                ncol=len(handles), fontsize=FONT_LEGEND - 1, frameon=True)

    plt.tight_layout(rect=[0, 0, 1, 0.90])
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    print(f"  → saved {output_path}")
    plt.close(fig)


def plot_era_difference_profile(era_lrrs: pd.DataFrame,
                                 preinstall_lrr: pd.DataFrame,
                                 transects: gpd.GeoDataFrame,
                                 groin: dict,
                                 output_path: str):
    """Difference profile: the CHANGE in shoreline change rate between
    consecutive eras, computed PER TRANSECT (matched by transect_id) --
    no binning of any kind. A transect only appears in a given
    difference if it has a valid rate in BOTH eras being compared.
    Both comparisons use the same "later minus earlier" convention:
        - Functional groin − Pre-install baseline
          (did the groin change the rate relative to background?)
        - Deteriorated − Functional groin
          (positive = the rate was higher after deterioration than
          while the groin was still functional; negative = the rate
          was higher during the functional era, i.e. effectiveness
          declined)

    Points are connected in distance order, same as the other profile
    plots. Worth keeping in mind here specifically: a difference of two
    independently-noisy per-transect estimates can swing more sharply
    transect-to-transect than either estimate alone, so this line may
    look jumpier than the era profile's -- that's the data, not a
    rendering choice.
    """
    print(f"\n[plot] era difference profile → {output_path}")

    pre = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                          preinstall_lrr["slope_m_yr"].notna()][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    func = era_lrrs[era_lrrs["era"] == "Functional groin"][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    det = era_lrrs[era_lrrs["era"] == "Deteriorated"][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()

    if len(pre) == 0 or len(func) == 0 or len(det) == 0:
        print("  ! missing one of pre-install/functional/deteriorated -- skipping")
        return

    pre["transect_id"] = pre["transect_id"].astype(str)
    func["transect_id"] = func["transect_id"].astype(str)
    det["transect_id"] = det["transect_id"].astype(str)

    print("  [diagnostic] Pre-install vs Functional groin:")
    _diagnose_era_overlap("Pre-install", pre["transect_id"],
                           "Functional groin", func["transect_id"])
    print("  [diagnostic] Functional groin vs Deteriorated:")
    _diagnose_era_overlap("Functional groin", func["transect_id"],
                           "Deteriorated", det["transect_id"])

    # Functional − Pre-install: only transects present in BOTH.
    # func's own dist_from_groin_m is dropped before merging -- it's
    # the same transect's same position as pre's copy, so keeping both
    # would just get them suffixed apart (dist_from_groin_m_pre/_func)
    # instead of leaving one plain dist_from_groin_m column to plot.
    m1 = pre.merge(func.drop(columns=["dist_from_groin_m"]),
                    on="transect_id", suffixes=("_pre", "_func"))
    m1["diff"] = m1["slope_m_yr_func"] - m1["slope_m_yr_pre"]

    # Deteriorated − Functional: only transects present in BOTH
    m2 = func.merge(det.drop(columns=["dist_from_groin_m"]),
                     on="transect_id", suffixes=("_func", "_det"))
    m2["diff"] = m2["slope_m_yr_det"] - m2["slope_m_yr_func"]

    print(f"  Functional − Pre-install: {len(m1)} transects with a valid "
          f"rate in both eras")
    print(f"  Deteriorated − Functional: {len(m2)} transects with a valid "
          f"rate in both eras")
    _diagnose_alongshore_gaps(m1, "dist_from_groin_m",
                               "Pre-install", pre, "Functional groin", func)
    _diagnose_alongshore_gaps(m2, "dist_from_groin_m",
                               "Functional groin", func, "Deteriorated", det)

    fig, ax = plt.subplots(figsize=(14, 7))

    m1 = m1.sort_values("dist_from_groin_m")
    m2 = m2.sort_values("dist_from_groin_m")

    ax.plot(m1["dist_from_groin_m"] / 1000, m1["diff"],
             color=ERA_COLORS["Functional groin"], lw=1.1, alpha=0.55, zorder=4)
    ax.scatter(m1["dist_from_groin_m"] / 1000, m1["diff"],
                s=14, color=ERA_COLORS["Functional groin"], alpha=0.45,
                edgecolor="none",
                label=f"Functional groin − Pre-install baseline  "
                       f"(n={len(m1)}; did the groin change the rate?)")
    ax.plot(m2["dist_from_groin_m"] / 1000, m2["diff"],
             color=ERA_COLORS["Deteriorated"], lw=1.1, alpha=0.55, zorder=4)
    ax.scatter(m2["dist_from_groin_m"] / 1000, m2["diff"],
                s=14, color=ERA_COLORS["Deteriorated"], alpha=0.45,
                edgecolor="none",
                label=f"Deteriorated − Functional groin  "
                       f"(n={len(m2)}; positive = rate higher after "
                       f"deterioration)")

    ax.axhline(0, color="#999", lw=0.7, linestyle="-")
    ax.axvline(0, color="black", lw=1.0, alpha=0.7,
                label=f"Groin  (updrift = "
                       f"{'north' if UPDRIFT_DIRECTION == 'north' else 'south'})")
    if groin.get("geometry") is not None:
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        g_a = sign * (groin["y_min"] - groin["reference_y"]) / 1000
        g_b = sign * (groin["y_max"] - groin["reference_y"]) / 1000
        ax.axvspan(min(g_a, g_b), max(g_a, g_b),
                    color="black", alpha=0.12, zorder=0,
                    label="Groin field  (first to last groin)")

    d1_km = domain_dist_km(transects, 1)
    apply_distance_xlim_and_autoscale_y(ax, x_max_km=PLOT_UPDRIFT_MAX_KM,
                                          x_min_km=d1_km)
    x_lo, x_hi = ax.get_xlim()
    set_domain_primary_axis(ax, transects)
    add_geographic_annotations(ax, transects)

    ax.set_ylabel("Δ LRR  (m/yr)", fontsize=FONT_AXIS_LABEL)
    ax.set_title(
        "Era-to-Era Difference in Shoreline Change Rate",
        fontsize=FONT_TITLE, fontweight="bold", pad=10,
    )
    ax.grid(True, alpha=0.25, linewidth=0.4)
    handles, labels = ax.get_legend_handles_labels()
    handles += annotation_legend_handles()
    ax.legend(handles=handles, loc="best", fontsize=FONT_LEGEND, framealpha=0.92, ncol=2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _diagnose_era_overlap(name_a: str, ids_a, name_b: str, ids_b):
    """Print how many transects overlap between two era-filtered sets
    before merging them for a difference calculation. ids_a/ids_b
    should largely be the SAME underlying transect set (just filtered
    to whichever era each has enough observations in) -- if one era's
    valid-transect count is close to the full transect total, the
    overlap with ANY other era's valid set is mathematically guaranteed
    to be at least that other set's own size. If the printed overlap is
    well below that, something's off in how the two sets are being
    matched, not just "this is how the data happens to look" -- this
    print exists to make that visible rather than silently producing a
    plot that only covers part of the coast.
    """
    set_a, set_b = set(ids_a), set(ids_b)
    overlap = set_a & set_b
    print(f"  {name_a}: {len(set_a)} unique transects with a valid rate")
    print(f"  {name_b}: {len(set_b)} unique transects with a valid rate")
    print(f"  In both: {len(overlap)} transects")
    smaller = min(len(set_a), len(set_b))
    if smaller > 0 and len(overlap) < smaller * 0.9:
        only_a = sorted(set_a - set_b)[:5]
        only_b = sorted(set_b - set_a)[:5]
        print(f"  ! overlap ({len(overlap)}) is well below the smaller "
              f"individual count ({smaller}) -- if {name_a} and {name_b} "
              f"are supposed to be nearly the same transects at "
              f"different times, this is worth investigating rather "
              f"than assumed correct. Examples only in {name_a}: "
              f"{only_a}; only in {name_b}: {only_b}")


def _diagnose_alongshore_gaps(merged_df: pd.DataFrame, dist_col: str,
                               name_a: str, full_a: pd.DataFrame,
                               name_b: str, full_b: pd.DataFrame,
                               gap_threshold_km: float = 1.5):
    """After merging two era-filtered tables on transect_id, check for
    ALONGSHORE GAPS in the result -- stretches with no matched data at
    all. On a connected-line plot, a gap like this shows up as a
    misleadingly smooth, dead-straight "bridge" jumping from the last
    real point before it to the first real point after it, which can
    look like a plotting bug even though it's actually just a gap in
    the merged data being connected across.

    For each gap found, checks whether name_a's or name_b's FULL
    (pre-merge, not-yet-intersected) table has any data in that same
    stretch -- so instead of leaving "why is there a gap" a mystery,
    this identifies whether ONE era genuinely lacks coverage there
    (e.g. pre-install has no historical shorelines pre-dating CoastSat
    in that stretch) or whether both sides have data that simply isn't
    matching by transect_id (which WOULD point to a real bug worth
    investigating further).
    """
    if len(merged_df) < 2:
        return
    d = np.sort(merged_df[dist_col].values) / 1000.0
    gaps = np.diff(d)
    gap_idx = np.where(gaps > gap_threshold_km)[0]
    if len(gap_idx) == 0:
        return
    print(f"  [diagnostic] {len(gap_idx)} alongshore gap(s) > "
          f"{gap_threshold_km} km in the {name_a}/{name_b} match "
          f"(shows up as a straight-line bridge on the plot):")
    for i in gap_idx:
        lo, hi = float(d[i]), float(d[i + 1])
        # STRICTLY inside the gap (exclusive bounds) -- lo and hi are
        # themselves merged/intersected points, so they're guaranteed
        # to exist in BOTH full_a and full_b; inclusive bounds would
        # always count at least those 2 boundary points and falsely
        # suggest "both sides have data" regardless of the interior.
        a_has = int(((full_a[dist_col] / 1000 > lo) &
                      (full_a[dist_col] / 1000 < hi)).sum())
        b_has = int(((full_b[dist_col] / 1000 > lo) &
                      (full_b[dist_col] / 1000 < hi)).sum())
        if a_has == 0 and b_has > 0:
            cause = f"{name_a} has no coverage there -- likely a real data gap in that era"
        elif b_has == 0 and a_has > 0:
            cause = f"{name_b} has no coverage there -- likely a real data gap in that era"
        elif a_has == 0 and b_has == 0:
            cause = "neither era has coverage there -- real gap in both"
        else:
            cause = (f"BOTH sides have data there ({name_a}={a_has}, "
                      f"{name_b}={b_has}) but it isn't matching by "
                      f"transect_id -- worth checking transect_id "
                      f"consistency between these two tables")
        print(f"    {lo:+.2f} to {hi:+.2f} km: {name_a} alone has "
              f"{a_has} transect(s) there, {name_b} alone has {b_has} "
              f"-- {cause}")



def _plot_single_difference(dist_km, diff_values, n: int, color: str,
                             series_label: str, question_text: str,
                             title_prefix: str, transects: gpd.GeoDataFrame,
                             groin: dict, output_path: str,
                             curve_a_values=None, curve_a_label: str = None,
                             curve_a_color: str = None,
                             curve_b_values=None, curve_b_label: str = None,
                             curve_b_color: str = None):
    """Shared rendering for a single (one-line) era-difference plot --
    used by plot_pre_to_functional_difference() and
    plot_functional_to_deteriorated_difference(), which each isolate
    ONE of the two comparisons already shown together on
    plot_era_difference_profile(), as their own focused figure.

    `color` is the dedicated DELTA color (see DELTA_COLOR) -- used for
    the fill between the two curves, the Δ line, and its smoothed
    overlay. It's deliberately NOT the same as either curve_a_color or
    curve_b_color, so the difference is never visually confusable with
    one of the two raw curves it was computed from.

    When curve_a/curve_b are provided, the two ORIGINAL curves being
    subtracted are drawn too (thin reference lines in their own era
    colors), with the gap between them shaded -- so the difference is
    visible both as an explicit Δ line AND as "how far apart these two
    curves actually sit", which is what makes a delta concrete.
    """
    fig, ax = plt.subplots(figsize=(14, 7))

    # The two original curves + shaded gap between them (drawn first,
    # low zorder, so the Δ line and its smoothed overlay stay on top)
    if curve_a_values is not None and curve_b_values is not None:
        ax.fill_between(dist_km, curve_a_values, curve_b_values,
                          color=color, alpha=0.12, zorder=1)
        ax.plot(dist_km, curve_a_values, color=curve_a_color, lw=1.3,
                 alpha=0.75, zorder=2, label=curve_a_label)
        ax.plot(dist_km, curve_b_values, color=curve_b_color, lw=1.3,
                 alpha=0.75, zorder=2, label=curve_b_label)

    ax.plot(dist_km, diff_values, color=color, lw=0.7, alpha=0.35, zorder=4)
    ax.scatter(dist_km, diff_values, s=8, color=color, alpha=0.35,
                edgecolor="none", label=f"{series_label}  (n={n})")
    xs, ys = _lowess_overlay(np.asarray(dist_km), np.asarray(diff_values),
                               SMOOTHED_OVERLAY_LOWESS_FRAC)
    if xs is not None:
        ax.plot(xs, ys, color=color, lw=2.8, zorder=6,
                 label=f"{series_label}  (smoothed)")

    ax.axhline(0, color="#999", lw=0.7, linestyle="-")
    ax.axvline(0, color="black", lw=1.0, alpha=0.7,
                label=f"Groin  (updrift = "
                       f"{'north' if UPDRIFT_DIRECTION == 'north' else 'south'})")
    if groin.get("geometry") is not None:
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        g_a = sign * (groin["y_min"] - groin["reference_y"]) / 1000
        g_b = sign * (groin["y_max"] - groin["reference_y"]) / 1000
        ax.axvspan(min(g_a, g_b), max(g_a, g_b), color="black", alpha=0.12,
                    zorder=0, label="Groin field  (first to last groin)")

    d1_km = domain_dist_km(transects, 1)
    apply_distance_xlim_and_autoscale_y(ax, x_max_km=PLOT_UPDRIFT_MAX_KM,
                                          x_min_km=d1_km)
    x_lo, x_hi = ax.get_xlim()
    set_domain_primary_axis(ax, transects)
    add_geographic_annotations(ax, transects)

    ax.set_ylabel("LRR (m/yr)  /  Δ LRR (m/yr)", fontsize=FONT_AXIS_LABEL)
    ax.set_title(
        title_prefix,
        fontsize=FONT_TITLE, fontweight="bold", pad=10,
    )
    ax.grid(True, alpha=0.25, linewidth=0.4)
    handles, labels = ax.get_legend_handles_labels()
    handles += annotation_legend_handles()
    ax.legend(handles=handles, loc="best", fontsize=FONT_LEGEND, framealpha=0.92, ncol=2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_pre_to_functional_difference(era_lrrs: pd.DataFrame,
                                       preinstall_lrr: pd.DataFrame,
                                       transects: gpd.GeoDataFrame,
                                       groin: dict, output_path: str):
    """Single-comparison plot: Functional groin era LRR minus pre-
    install baseline LRR, matched per transect (both required). Isolates
    just this ONE comparison from plot_era_difference_profile() as its
    own focused figure -- did the groin change the rate relative to
    background?
    """
    print(f"\n[plot] pre-to-functional difference profile → {output_path}")
    pre = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                          preinstall_lrr["slope_m_yr"].notna()][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    func = era_lrrs[era_lrrs["era"] == "Functional groin"][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    if len(pre) == 0 or len(func) == 0:
        print("  ! missing pre-install or functional era data -- skipping")
        return
    pre["transect_id"] = pre["transect_id"].astype(str)
    func["transect_id"] = func["transect_id"].astype(str)
    _diagnose_era_overlap("Pre-install", pre["transect_id"],
                           "Functional groin", func["transect_id"])
    m = pre.merge(func.drop(columns=["dist_from_groin_m"]),
                   on="transect_id", suffixes=("_pre", "_func"))
    m["diff"] = m["slope_m_yr_func"] - m["slope_m_yr_pre"]
    m = m.sort_values("dist_from_groin_m")
    print(f"  Functional − Pre-install: {len(m)} transects with a valid "
          f"rate in both eras")
    _diagnose_alongshore_gaps(m, "dist_from_groin_m",
                               "Pre-install", pre, "Functional groin", func)

    _plot_single_difference(
        m["dist_from_groin_m"] / 1000, m["diff"], len(m),
        DELTA_COLOR,
        "Functional groin − Pre-install baseline",
        "did the groin change the rate relative to background?",
        "Did the groin change the shoreline change rate?",
        transects, groin, output_path,
        curve_a_values=m["slope_m_yr_pre"],
        curve_a_label=f"Pre-install baseline  (n={len(m)})",
        curve_a_color=ERA_COLORS["Pre-install"],
        curve_b_values=m["slope_m_yr_func"],
        curve_b_label=f"Functional groin  (n={len(m)})",
        curve_b_color=ERA_COLORS["Functional groin"])


def plot_deteriorated_to_preinstall_difference(era_lrrs: pd.DataFrame,
                                                preinstall_lrr: pd.DataFrame,
                                                transects: gpd.GeoDataFrame,
                                                groin: dict, output_path: str):
    """Single-comparison plot: Deteriorated era LRR minus pre-install
    baseline LRR (later minus earlier, same convention as the other two
    focused difference plots), matched per transect (both required).
    Spans the WHOLE historical arc -- before the groin existed at all,
    to its current deteriorated state -- skipping over the middle
    "Functional groin" era entirely. Positive means the rate is higher
    now than it was before the groin was ever built; negative means
    the opposite.
    """
    print(f"\n[plot] deteriorated-to-preinstall difference profile → {output_path}")
    pre = preinstall_lrr[(preinstall_lrr["n_obs"] >= 3) &
                          preinstall_lrr["slope_m_yr"].notna()][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    det = era_lrrs[era_lrrs["era"] == "Deteriorated"][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    if len(pre) == 0 or len(det) == 0:
        print("  ! missing pre-install or deteriorated era data -- skipping")
        return
    pre["transect_id"] = pre["transect_id"].astype(str)
    det["transect_id"] = det["transect_id"].astype(str)
    _diagnose_era_overlap("Pre-install", pre["transect_id"],
                           "Deteriorated", det["transect_id"])
    m = pre.merge(det.drop(columns=["dist_from_groin_m"]),
                   on="transect_id", suffixes=("_pre", "_det"))
    m["diff"] = m["slope_m_yr_det"] - m["slope_m_yr_pre"]
    m = m.sort_values("dist_from_groin_m")
    print(f"  Deteriorated − Pre-install: {len(m)} transects with a valid "
          f"rate in both eras")
    _diagnose_alongshore_gaps(m, "dist_from_groin_m",
                               "Pre-install", pre, "Deteriorated", det)

    _plot_single_difference(
        m["dist_from_groin_m"] / 1000, m["diff"], len(m),
        DELTA_COLOR,
        "Deteriorated − Pre-install baseline",
        "positive = the rate today is higher than it was before the "
        "groin was ever built",
        "How has the rate changed across the whole historical record?",
        transects, groin, output_path,
        curve_a_values=m["slope_m_yr_pre"],
        curve_a_label=f"Pre-install baseline  (n={len(m)})",
        curve_a_color=ERA_COLORS["Pre-install"],
        curve_b_values=m["slope_m_yr_det"],
        curve_b_label=f"Deteriorated  (n={len(m)})",
        curve_b_color=ERA_COLORS["Deteriorated"])


def plot_functional_to_deteriorated_difference(era_lrrs: pd.DataFrame,
                                                transects: gpd.GeoDataFrame,
                                                groin: dict, output_path: str):
    """Single-comparison plot: Deteriorated era LRR minus Functional
    groin era LRR (later minus earlier, matching the same convention as
    Functional − Pre-install), matched per transect (both required).
    Isolates just this ONE comparison from
    plot_era_difference_profile() as its own focused figure -- positive
    means the rate was higher after deterioration than while the groin
    was still functional; negative means the rate was higher during the
    functional era (effectiveness declined).
    """
    print(f"\n[plot] functional-to-deteriorated difference profile → {output_path}")
    func = era_lrrs[era_lrrs["era"] == "Functional groin"][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    det = era_lrrs[era_lrrs["era"] == "Deteriorated"][
        ["transect_id", "dist_from_groin_m", "slope_m_yr"]].copy()
    if len(func) == 0 or len(det) == 0:
        print("  ! missing functional or deteriorated era data -- skipping")
        return
    func["transect_id"] = func["transect_id"].astype(str)
    det["transect_id"] = det["transect_id"].astype(str)
    _diagnose_era_overlap("Functional groin", func["transect_id"],
                           "Deteriorated", det["transect_id"])
    m = func.merge(det.drop(columns=["dist_from_groin_m"]),
                    on="transect_id", suffixes=("_func", "_det"))
    m["diff"] = m["slope_m_yr_det"] - m["slope_m_yr_func"]
    m = m.sort_values("dist_from_groin_m")
    print(f"  Deteriorated − Functional: {len(m)} transects with a valid "
          f"rate in both eras")
    _diagnose_alongshore_gaps(m, "dist_from_groin_m",
                               "Functional groin", func, "Deteriorated", det)

    _plot_single_difference(
        m["dist_from_groin_m"] / 1000, m["diff"], len(m),
        DELTA_COLOR,
        "Deteriorated − Functional groin",
        "positive = rate was higher after deterioration than during "
        "the functional era",
        "Did the rate change again as the groin deteriorated?",
        transects, groin, output_path,
        curve_a_values=m["slope_m_yr_func"],
        curve_a_label=f"Functional groin  (n={len(m)})",
        curve_a_color=ERA_COLORS["Functional groin"],
        curve_b_values=m["slope_m_yr_det"],
        curve_b_label=f"Deteriorated  (n={len(m)})",
        curve_b_color=ERA_COLORS["Deteriorated"])


def compute_distance_band_series(decadal_df: pd.DataFrame,
                                  value_col: str) -> pd.DataFrame:
    """Groups transects into named alongshore distance bands
    (DISTANCE_BAND_EDGES_M), separately updrift and downdrift, and
    computes the median value_col per band per decade. CSV-only now
    (the companion plot was removed as unhelpful) -- kept because the
    numbers may still be useful even without a chart.

    Returns long-format DataFrame:
        decade_start, side ('updrift'/'downdrift'), band_label,
        band_lo_m, band_hi_m, n, median
    """
    edges = DISTANCE_BAND_EDGES_M
    df = decadal_df.copy()
    df["side"] = np.where(df["dist_from_groin_m"] >= 0, "updrift", "downdrift")
    df["abs_dist_m"] = df["dist_from_groin_m"].abs()

    records = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        label = f"{lo/1000:.0f}\u2013{hi/1000:.0f} km"
        band = df[(df["abs_dist_m"] >= lo) & (df["abs_dist_m"] < hi)]
        for (side, cy), sgrp in band.groupby(["side", "decade_start"]):
            vals = sgrp[value_col].dropna()
            if len(vals) == 0:
                continue
            records.append({
                "decade_start": cy, "side": side, "band_label": label,
                "band_lo_m": lo, "band_hi_m": hi,
                "n": len(vals), "median": float(vals.median()),
            })
    return pd.DataFrame(records)


# ============================================================
# SHORELINE EVOLUTION GIF (groin area)
# ============================================================

def create_groin_evolution_gif(chainage: pd.DataFrame,
                                transects: gpd.GeoDataFrame,
                                groin: dict,
                                output_dir: str,
                                domain_table: pd.DataFrame = None,
                                grid100m: gpd.GeoDataFrame = None,
                                window_domain_max: int = None,
                                output_suffix: str = ""):
    """Animated GIF of RAW shoreline position (chainage, m -- no
    correction/centering of any kind) through time, zoomed to
    domain 1 to PLOT_UPDRIFT_MAX_KM updrift by default (same window as
    the static profile plots). Pass window_domain_max (e.g. 25) for a
    tighter, domain-bounded window instead -- domain 1 to whatever
    km position that domain maps to -- e.g. for a closer look at
    change right around the groin rather than the whole island stretch.
    output_suffix (e.g. "_zoomed") distinguishes the comparison filename/
    frame folder when generating more than one GIF in the same run.
    Adapted from the whole-island
    version in HAT_shoreline_chainage_alldata_evolution.py, but built
    on this script's multi-source unified chainage table (wet-dry + NC
    state + CoastSat, already merged in `chainage`) instead of CoastSat
    alone, so the animation can show the full historical record instead
    of starting at 1984.

    GIF_TRANSECT_SOURCE controls which transects define alongshore
    position and/or measure shoreline position (see its config comment
    for the full explanation of "coastsat" / "hybrid" / "grid100m").
    domain_table and grid100m are only needed for "hybrid"/"grid100m"
    -- pass None (default) when using "coastsat".

    CoastSat's transects are shore-normal and run parallel to the
    coast, so the chainage values are used exactly as provided -- no
    geometric correction of any kind.

    GIF_FRAME_MODE picks "date" (default -- one frame per unique
    observation date, using every observation the data supports) or
    "year" (pools all sources/dates within a calendar year into one
    frame -- fuller per-frame spatial coverage, but coarser temporal
    resolution). In "date" mode, a single wet-dry or NC-state date, or
    an early/partial CoastSat pass, often only covers PART of the
    window, so individual frames can look spatially sparse -- that's
    real data sparsity on that date, not a bug.

    x-axis is CASCADE domain number (primary), with distance from groin
    (km) on a secondary top axis. y-axis is fixed to the TRUE min/max
    chainage across the entire animation (not a per-frame or percentile
    range), so the real shoreline position is always fully in view and
    the axis never rescales frame to frame -- and INVERTED, so the
    island (smaller chainage, landward) renders at the top and the
    ocean (larger chainage, seaward) at the bottom. Everything seaward
    of the current shoreline curve is shaded like water
    (GIF_WATER_COLOR), everything landward like the island itself
    (GIF_LAND_COLOR) -- the curve is the actual land/water boundary at
    that alongshore position and moment in time, not just an abstract
    line, and the shading updates every frame as that boundary moves.

    The current frame's line is colored by structural era (matches
    ERA_COLORS elsewhere); the preceding GIF_TRAIL_FRAMES frames are
    drawn as a fading grey trail (unshaded, so they don't visually
    compete with the current frame's water/island fill). A fixed,
    always-visible dashed reference line shows the GIF_REFERENCE_YEAR
    (1967) shoreline -- the last shoreline before the groin was built --
    on every frame dated after it.

    Requires Pillow (pip install pillow). Writes nothing and returns
    None if Pillow is unavailable or there's too little data in the
    window to animate.
    """
    try:
        from PIL import Image
    except ImportError:
        print("\n[gif] Pillow not installed (pip install pillow) -- "
              "skipping shoreline evolution GIF")
        return None

    if GIF_TRANSECT_SOURCE not in ("coastsat", "hybrid", "grid100m"):
        raise ValueError(f"GIF_TRANSECT_SOURCE must be 'coastsat', "
                          f"'hybrid', or 'grid100m', got "
                          f"{GIF_TRANSECT_SOURCE!r}")

    # x_max_km: either the standard PLOT_UPDRIFT_MAX_KM (full-island
    # window, same as the static plots) or, if window_domain_max is
    # given, wherever that specific domain maps to -- e.g. a tighter
    # domain 1-25 window zoomed in on the groin area.
    if window_domain_max is not None:
        x_max_km = domain_dist_km(transects, window_domain_max)
        if x_max_km is None:
            print(f"  ! domain {window_domain_max} has no transects -- "
                  f"falling back to PLOT_UPDRIFT_MAX_KM for the window")
            x_max_km = PLOT_UPDRIFT_MAX_KM
        window_desc = f"domain 1 to domain {window_domain_max} (~{x_max_km:+.1f} km)"
    else:
        x_max_km = PLOT_UPDRIFT_MAX_KM
        window_desc = f"domain 1 to +{PLOT_UPDRIFT_MAX_KM} km updrift"

    print(f"\n{'='*72}\nShoreline evolution GIF (groin area, {window_desc})")
    print(f"  Transect source: {GIF_TRANSECT_SOURCE}")

    frame_dir = os.path.join(output_dir, f"gif_frames_groin_area{output_suffix}")
    os.makedirs(frame_dir, exist_ok=True)

    effective_source = GIF_TRANSECT_SOURCE
    if GIF_TRANSECT_SOURCE in ("hybrid", "grid100m") and grid100m is None:
        print(f"  ! GIF_TRANSECT_SOURCE={GIF_TRANSECT_SOURCE!r} but no 100m "
              f"grid was loaded -- falling back to 'coastsat'")
        effective_source = "coastsat"

    # ─── Restrict to the window computed above: domain 1 (downdrift) ──
    # to x_max_km (updrift) of the groin, in CoastSat's own coordinate
    # system (transects already has this applied from main()'s Cape
    # Point exclusion).
    d1_km = domain_dist_km(transects, 1)
    x_min_km = d1_km if d1_km is not None else -PLOT_UPDRIFT_MAX_KM
    tx = transects[(transects["dist_from_groin_m"] / 1000 >= x_min_km) &
                   (transects["dist_from_groin_m"] / 1000 <= x_max_km)].copy()
    tx = tx.sort_values("dist_from_groin_m").reset_index(drop=True)
    print(f"  {len(tx)} CoastSat transects from domain 1 ({x_min_km:+.2f} km) "
          f"to {x_max_km:+.1f} km updrift of the groin")
    if len(tx) < 5:
        print("  ! too few transects in window -- skipping GIF")
        return None

    sub = chainage[chainage["transect_id"].astype(str).isin(
        tx["transect_id"].astype(str))].copy()
    sub["transect_id"] = sub["transect_id"].astype(str)
    sub["date"] = pd.to_datetime(sub["date"])
    if len(sub) == 0:
        print("  ! no chainage data in the groin-area window -- skipping GIF")
        return None

    if effective_source == "coastsat":
        tids_ordered = tx["transect_id"].astype(str).tolist()
        dist_km = tx["dist_from_groin_m"].values / 1000.0

    elif effective_source == "hybrid":
        # CoastSat's own (correctly-measured) chainage values, `sub`,
        # are untouched -- only WHERE each transect sits on the x-axis
        # changes, via nearest-neighbor matching to the 100 m grid.
        print("  Repositioning CoastSat transects onto the 100 m grid's "
              "alongshore coordinates (measurements unchanged)...")
        tx["transect_id"] = tx["transect_id"].astype(str)
        pos_map = build_hybrid_position_map(tx, grid100m, groin)
        tx["dist_from_groin_m"] = tx["transect_id"].map(pos_map)
        tx = tx.dropna(subset=["dist_from_groin_m"])
        tx = tx.sort_values("dist_from_groin_m").reset_index(drop=True)
        tids_ordered = tx["transect_id"].tolist()
        dist_km = tx["dist_from_groin_m"].values / 1000.0
        print(f"  → repositioned {len(tx)} transects "
              f"({dist_km.min():+.2f} to {dist_km.max():+.2f} km)")

    else:   # "grid100m"
        print("  ! reconstructing shoreline points and re-measuring "
              "against the 100 m grid's own (parallel) transects -- "
              "this reintroduces a known geometric artifact, see "
              "GIF_TRANSECT_SOURCE config comment")
        sub = reconstruct_and_remeasure_on_grid100m(sub, tx, grid100m)
        if len(sub) == 0:
            print("  ! no chainage survived re-measurement -- skipping GIF")
            return None

        # Window/order the 100 m-grid transects the same way: domain 1
        # to PLOT_UPDRIFT_MAX_KM, converted directly via the domain
        # table's authoritative boundary -- exact here (no curvilinear
        # ambiguity), since this grid's alongshore position IS northing.
        gx, gy = groin.get("reference_x"), groin["reference_y"]
        if gx is not None:
            d2 = (grid100m["shore_x"] - gx) ** 2 + (grid100m["shore_y"] - gy) ** 2
        else:
            d2 = (grid100m["shore_y"] - gy) ** 2
        groin_along_100m = float(grid100m.loc[d2.idxmin(), "alongshore_m"])
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        grid100m = grid100m.copy()
        grid100m["dist_from_groin_m"] = sign * (grid100m["alongshore_m"] - groin_along_100m)

        x_min_km_100m = -PLOT_UPDRIFT_MAX_KM
        if domain_table is not None:
            d1_row = domain_table[domain_table["domain_id"] == 1]
            if len(d1_row) > 0:
                d1_y_min = float(d1_row["y_min"].iloc[0])
                grid_y_min = float(grid100m["shore_y"].min())
                d1_along_100m = d1_y_min - grid_y_min
                x_min_km_100m = sign * (d1_along_100m - groin_along_100m) / 1000

        grid_window = grid100m[
            (grid100m["dist_from_groin_m"] / 1000 >= x_min_km_100m) &
            (grid100m["dist_from_groin_m"] / 1000 <= x_max_km)
        ].sort_values("dist_from_groin_m").reset_index(drop=True)
        tids_ordered = grid_window["transect_id"].astype(str).tolist()
        dist_km = grid_window["dist_from_groin_m"].values / 1000.0
        tx = grid_window   # used downstream by set_domain_primary_axis / add_geographic_annotations
        print(f"  → {len(tids_ordered)} 100m-grid transects in window "
              f"({x_min_km_100m:+.2f} to {x_max_km:.0f} km)")

    # ─── One row of raw chainage per frame (date or year) ─────────────
    if GIF_FRAME_MODE == "year":
        frame_keys_all = sorted(sub["year"].unique())
        print(f"  {len(frame_keys_all)} unique years in the window "
              f"({frame_keys_all[0]} → {frame_keys_all[-1]}), "
              f"pooling all sources/dates within each year")
    elif GIF_FRAME_MODE == "date":
        frame_keys_all = sorted(sub["date"].unique())
        print(f"  {len(frame_keys_all)} unique observation dates in the "
              f"window ({frame_keys_all[0].date()} → "
              f"{frame_keys_all[-1].date()})")
    else:
        raise ValueError(f"GIF_FRAME_MODE must be 'year' or 'date', "
                          f"got {GIF_FRAME_MODE!r}")

    frames_data = {}
    n_transects_per_key = {}
    for k in frame_keys_all:
        day_data = sub[sub["year"] == k] if GIF_FRAME_MODE == "year" \
                   else sub[sub["date"] == k]
        agg = day_data.groupby("transect_id")["chainage_m"].median()
        row = np.full(len(tids_ordered), np.nan)
        for i, tid in enumerate(tids_ordered):
            if tid in agg.index:
                row[i] = agg.loc[tid]
        n_present = int(np.sum(~np.isnan(row)))
        n_transects_per_key[k] = n_present
        frame_year = k if GIF_FRAME_MODE == "year" else k.year
        required_fraction = (0.0 if frame_year < GIF_PRE_COASTSAT_CUTOFF_YEAR
                             else GIF_MIN_TRANSECT_FRACTION)
        if (n_present >= required_fraction * len(tids_ordered)
                and n_present >= GIF_MIN_TRANSECT_ABS):
            frames_data[k] = row

    valid_keys = sorted(frames_data.keys())
    unit = "years" if GIF_FRAME_MODE == "year" else "dates"
    n_pre_coastsat_kept = sum(1 for k in valid_keys
                                if (k if GIF_FRAME_MODE == "year" else k.year)
                                < GIF_PRE_COASTSAT_CUTOFF_YEAR)
    print(f"  {len(valid_keys)}/{len(frame_keys_all)} {unit} have enough "
          f"transect coverage to render a frame: before "
          f"{GIF_PRE_COASTSAT_CUTOFF_YEAR} every {unit[:-1]} with >= "
          f"{GIF_MIN_TRANSECT_ABS} transects is kept regardless of "
          f"coverage fraction ({n_pre_coastsat_kept} such {unit} kept); "
          f"{GIF_PRE_COASTSAT_CUTOFF_YEAR} onward requires "
          f">= {GIF_MIN_TRANSECT_FRACTION:.0%} of "
          f"{len(tids_ordered)} transects")
    if len(frame_keys_all) > 0:
        counts = np.array(list(n_transects_per_key.values()))
        print(f"  Coverage per {unit[:-1]}: median {int(np.median(counts))}, "
              f"a few dates would additionally qualify at lower "
              f"thresholds -- e.g. >=1: {int((counts >= 1).sum())}, "
              f">=5: {int((counts >= 5).sum())}, "
              f">=20: {int((counts >= 20).sum())} "
              f"(adjust GIF_MIN_TRANSECT_FRACTION / GIF_MIN_TRANSECT_ABS "
              f"if you want more or fewer frames)")
    if len(valid_keys) < 3:
        print(f"  ! too few valid {unit} -- skipping GIF "
              "(try lowering GIF_MIN_TRANSECT_FRACTION / GIF_MIN_TRANSECT_ABS)")
        return None

    # ─── Persistent reference baseline: the GIF_REFERENCE_YEAR shoreline ──
    # Built exactly like any other frame -- median chainage per transect
    # for that specific year, NaN (gap) where no observation exists that
    # year -- so it's exactly as sparse/honest as the real coverage.
    # Drawn on every frame dated after GIF_REFERENCE_YEAR as a fixed
    # comparison line (not part of the fading trail).
    ref_data = sub[sub["date"].dt.year == GIF_REFERENCE_YEAR]
    ref_row = np.full(len(tids_ordered), np.nan)
    if len(ref_data) > 0:
        ref_agg = ref_data.groupby("transect_id")["chainage_m"].median()
        for i, tid in enumerate(tids_ordered):
            if tid in ref_agg.index:
                ref_row[i] = ref_agg.loc[tid]
    n_ref = int(np.sum(~np.isnan(ref_row)))
    print(f"  {GIF_REFERENCE_YEAR} reference baseline: {n_ref}/{len(tx)} "
          f"transects in window have a {GIF_REFERENCE_YEAR} observation"
          + ("" if n_ref > 0 else f"  ! no {GIF_REFERENCE_YEAR} data in "
                                    f"this window -- reference line will "
                                    f"be blank"))

    # Diagnostic: per-transect total observation count across the
    # WHOLE historical record in this window. A gap that recurs across
    # many frames at the same alongshore position is usually because
    # those specific transects simply have much sparser coverage than
    # their neighbors (a real data characteristic worth knowing about),
    # not a processing bug -- printed here so that's checkable
    # directly rather than left a mystery.
    coverage_counts = sub.groupby("transect_id").size()
    coverage_df = pd.DataFrame({
        "transect_id": tids_ordered,
        "dist_km": dist_km,
        "n_obs": [coverage_counts.get(tid, 0) for tid in tids_ordered],
    })
    sparsest = coverage_df.nsmallest(10, "n_obs")
    print(f"  Per-transect total observation count in this window: "
          f"median {int(coverage_df['n_obs'].median())}; "
          f"sparsest 10 transects (total obs across the whole record):")
    for _, row in sparsest.iterrows():
        print(f"    {row['transect_id']}  ({row['dist_km']:+.2f} km)  "
              f"n_obs={row['n_obs']}")

    # Shared y-axis: TRUE min/max across every frame in the whole time
    # period (not a percentile-trimmed range), so the shoreline's real
    # movement is always fully in view, in every single frame -- never
    # rescaling between frames, and never clipping a real extreme.
    # Includes the reference row too.
    all_vals = np.concatenate(
        [v[~np.isnan(v)] for v in frames_data.values()] + [ref_row[~np.isnan(ref_row)]])
    y_lo, y_hi = float(np.nanmin(all_vals)), float(np.nanmax(all_vals))
    pad = (y_hi - y_lo) * 0.08 if y_hi > y_lo else 1.0
    ylim = (y_lo - pad, y_hi + pad)

    def era_for_year(yr):
        for name, lo, hi in ERAS:
            if lo <= yr <= hi:
                return name
        return None

    x_lo_plot, x_hi_plot = x_min_km, x_max_km

    # Groin field footprint (southernmost to northernmost groin),
    # relative to the northernmost-groin origin (0 km) -- generally NOT
    # symmetric around 0, since 0 sits at the north end of the field.
    groin_span = None
    if groin.get("geometry") is not None:
        sign = 1 if UPDRIFT_DIRECTION.lower().startswith("n") else -1
        g_a = sign * (groin["y_min"] - groin["reference_y"]) / 1000
        g_b = sign * (groin["y_max"] - groin["reference_y"]) / 1000
        groin_span = (min(g_a, g_b), max(g_a, g_b))

    print(f"  Rendering {len(valid_keys)} frames...")
    frame_paths = []
    for fi, k in enumerate(valid_keys):
        if (fi + 1) % 50 == 0 or fi == len(valid_keys) - 1:
            print(f"    {fi+1}/{len(valid_keys)}", end="\r")

        if GIF_FRAME_MODE == "year":
            frame_year, label_str = int(k), str(int(k))
        else:
            frame_year, label_str = k.year, str(k.date())

        fig, ax = plt.subplots(figsize=(12, 5.5))

        # Persistent reference baseline (drawn every frame after the
        # reference year, not part of the fading trail)
        if n_ref > 0 and frame_year > GIF_REFERENCE_YEAR:
            ax.plot(dist_km, ref_row, color="#222", lw=0.9, linestyle="--",
                     alpha=0.6, zorder=3,
                     label=f"{GIF_REFERENCE_YEAR} reference baseline")

        # Fading trail of preceding frames
        for tr in range(1, GIF_TRAIL_FRAMES + 1):
            j = fi - tr
            if j < 0:
                continue
            alpha = 0.30 * (1 - tr / (GIF_TRAIL_FRAMES + 1))
            ax.plot(dist_km, frames_data[valid_keys[j]], color="#999",
                     lw=1.0, alpha=alpha, zorder=2)

        row = frames_data[k]
        era_name = era_for_year(frame_year)
        color = ERA_COLORS.get(era_name, "#1565C0")

        # Water (seaward of the shoreline, larger chainage) / island
        # (landward of it, smaller chainage) shading -- the shoreline
        # curve itself is the boundary between them. Uses an
        # INTERPOLATED version of row for the shading only, so missing-
        # data gaps don't flash white between frames -- the actual line
        # below still uses the raw row and keeps its honest gaps.
        row_series = pd.Series(row)
        if row_series.notna().sum() >= 2:
            row_filled = row_series.interpolate(
                limit_direction="both").to_numpy()
        else:
            row_filled = row
        ax.fill_between(dist_km, row_filled, ylim[1], color=GIF_WATER_COLOR,
                          alpha=0.55, zorder=1, linewidth=0)
        ax.fill_between(dist_km, ylim[0], row_filled, color=GIF_LAND_COLOR,
                          alpha=0.65, zorder=1, linewidth=0)

        ax.plot(dist_km, row, color=color, lw=1.0, marker="o", ms=1.3,
                 zorder=5,
                 label=f"{label_str}" + (f"  ({era_name})" if era_name else ""))

        ax.axvline(0, color="black", lw=1.0, alpha=0.7,
                    label=f"Groin  (updrift = "
                           f"{'north' if UPDRIFT_DIRECTION == 'north' else 'south'})")
        if groin_span is not None:
            ax.axvspan(groin_span[0], groin_span[1], color="black",
                        alpha=0.12, zorder=0,
                        label="Groin field  (first to last groin)")

        ax.set_xlim(x_lo_plot, x_hi_plot)
        ax.set_ylim(ylim[1], ylim[0])   # inverted: island (small chainage) on top, ocean (large chainage) on bottom
        set_domain_primary_axis(ax, tx)
        add_geographic_annotations(ax, tx, gif_mode=True)

        ax.set_ylabel("Shoreline position\n(chainage, m)", fontsize=FONT_AXIS_LABEL)
        ax.set_title(f"Groin-area shoreline position  —  {label_str}",
                       fontsize=FONT_TITLE, fontweight="bold")
        handles, labels = ax.get_legend_handles_labels()
        handles += [Patch(facecolor=GIF_WATER_COLOR, alpha=0.55, label="Ocean"),
                     Patch(facecolor=GIF_LAND_COLOR, alpha=0.65, label="Island")]
        handles += annotation_legend_handles()
        if window_domain_max is not None:
            ax.legend(handles=handles, loc="upper center",
                       fontsize=FONT_LEGEND_GIF_ZOOMED, ncol=2)
        else:
            ax.legend(handles=handles, loc="upper right",
                       fontsize=FONT_LEGEND_GIF, ncol=2)
        ax.grid(True, alpha=0.25, linewidth=0.4)

        plt.tight_layout()
        frame_path = os.path.join(frame_dir, f"frame_{fi:04d}_{label_str}.png")
        fig.savefig(frame_path, dpi=GIF_DPI, bbox_inches="tight")
        plt.close(fig)
        frame_paths.append(frame_path)
    print()

    gif_path = os.path.join(
        output_dir, f"groin_analysis_shoreline_evolution{output_suffix}.gif")
    print(f"  Assembling {len(frame_paths)}-frame GIF...")
    frames = [Image.open(fp).convert("RGBA") for fp in frame_paths]
    frames[0].save(
        gif_path, save_all=True, append_images=frames[1:],
        duration=int(GIF_FRAME_DURATION_S * 500),  # *500 corrects a Pillow duration bug
        loop=0,
    )
    print(f"  → saved {gif_path}")
    return gif_path


# ============================================================
# METADATA
# ============================================================

def write_metadata(output_path: str, stats: dict):
    """Write a plain-text metadata file summarizing the run."""
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("HAT_groin_shoreline_analysis_v2.py — run metadata\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Run at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("--- CONFIG ---\n")
        f.write(f"Pre-install year cutoff:  < {PRE_INSTALLATION_YEAR_CUTOFF} "
                f"(all shorelines before this year used)\n")
        f.write(f"Groin geojson: {GROIN_GEOJSON_PATH}\n")
        f.write(f"Updrift direction: {UPDRIFT_DIRECTION}\n")
        f.write(f"Signal anomaly threshold: "
                f"{SIGNAL_ANOMALY_THRESHOLD_M_YR} m/yr\n")
        f.write(f"Signal max search distance: "
                f"{SIGNAL_MAX_SEARCH_DISTANCE_M/1000:.1f} km\n")
        f.write(f"Signal extent bin width: "
                f"{SIGNAL_EXTENT_BIN_WIDTH_M:.0f} m, "
                f"max gap {SIGNAL_EXTENT_MAX_GAP_BINS} bin(s), contiguous from groin\n")
        f.write(f"Baseline method: connect-the-dots linear interpolation "
                f"across per-transect pre-install rates (no binning/smoothing)\n")
        f.write(f"Decadal bins: {DECADE_LENGTH_YEARS} yr, "
                f"start {DECADE_START_YEAR}, "
                f"min obs {DECADE_MIN_OBSERVATIONS}\n")
        f.write(f"Breakpoint search: {BREAKPOINT_SEARCH_START}–"
                f"{BREAKPOINT_SEARCH_END}, "
                f"min pts/segment {BREAKPOINT_MIN_POINTS_PER_SEGMENT}\n")
        f.write(f"Shoreline evolution GIF window: domain 1 to "
                f"+{PLOT_UPDRIFT_MAX_KM} km updrift, "
                f"{GIF_FRAME_DURATION_S}s/frame, mode={GIF_FRAME_MODE}\n")
        f.write(f"CoastSat transect layer: {COASTSAT_TRANSECT_GEOM}\n")
        f.write(f"Transect length: {TRANSECT_INTERSECTION_LENGTH_M} m, "
                f"chainage max abs {CHAINAGE_MAX_ABS_M} m\n")
        f.write(f"CASCADE domain reference: {DOMAINS_JSON_PATH}\n")
        f.write(f"Plot x-axis clip: ±{PLOT_X_MAX_KM} km\n")
        f.write(f"Min obs per transect: {MIN_OBSERVATIONS_PER_TRANSECT}\n\n")
        f.write("--- STATS ---\n")
        for k, v in stats.items():
            f.write(f"{k}: {v}\n")


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 72)
    print("HAT GROIN SHORELINE ANALYSIS — Hatteras Island")
    print("=" * 72)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 1. Study area filter
    filter_gdf = load_study_area_filter()

    # 2. Load the authoritative CASCADE domain reference (replaces the
    # original CAPE_POINT_NORTHING formula) -- loaded before the groin so its
    # position can be sanity-checked against a known domain.
    domain_table = load_domain_reference()

    # 3. Load groin geometry (sets the origin for distance-from-groin)
    groin = load_groin_geometry(domain_table)

    # 4. Load shoreline sources
    wet_dry_gdf = load_shoreline_lines(WET_DRY_PATH, WET_DRY_DATE_COL,
                                        "wet_dry", filter_gdf)
    nc_state_gdf = load_shoreline_lines(NC_STATE_PATH, NC_STATE_DATE_COL,
                                         "nc_state", filter_gdf)

    # 5. Load CoastSat's own transect network -- the analysis backbone
    # for EVERY source below (wet-dry, NC state, and CoastSat itself),
    # per direction to use CoastSat's transects for all shoreline
    # change analyses.
    transects = load_coastsat_transects(filter_gdf, domain_table)
    transects = assign_distance_from_groin(transects, groin)
    print(f"  → distance-from-groin range: "
          f"{transects['dist_from_groin_m'].min()/1000:+.2f} to "
          f"{transects['dist_from_groin_m'].max()/1000:+.2f} km")
    n_updrift = int((transects["dist_from_groin_m"] > 0).sum())
    n_downdrift = int((transects["dist_from_groin_m"] < 0).sum())
    print(f"  → {n_updrift} updrift transects, {n_downdrift} downdrift "
          f"transects")

    # 5b. Exclude transects south of domain 1's TRUE boundary (per the
    # authoritative domains.json, not just "nearest domain by midpoint"
    # -- assign_domain_from_northing() has no distance cutoff, so a
    # transect far south of domain 1 still gets labeled domain=1 simply
    # because domain 1 is the closest available bucket). This is a
    # REAL data exclusion, not a display clip: it happens before any
    # chainage extraction or regression below, so it affects every
    # curve, every smoothing pass, and the GIF -- not just where the
    # plots' x-axis starts. Per direction: this stretch is close enough
    # to Cape Point's sharp curvature that alongshore ordering and
    # chainage baselines are much less reliable there, and it isn't
    # useful for this analysis anyway.
    #
    # IMPORTANT: the cutoff is applied to dist_from_groin_m (alongshore
    # position, built from the trusted CoastSat ID-order) -- the SAME
    # coordinate every plot's x-axis actually uses -- not to northing
    # directly. Northing and alongshore position are exactly the two
    # things that disagree near Cape Point's sharp curve (that's the
    # reason this exclusion exists at all), so filtering on northing
    # while plotting against alongshore position let some transects
    # through the filter that still showed up south of domain 1 on the
    # actual plots. The domain-1 northing boundary is converted to an
    # alongshore-position threshold via a linear fit across ALL
    # transects (robust to the local Cape Point anomalies, since the
    # fit is dominated by the well-behaved majority of the island),
    # rather than by averaging just the handful of transects nearest-
    # midpoint-labeled "domain 1" (which is itself an unstable quantity
    # right in the area being excluded).
    domain1_row = domain_table[domain_table["domain_id"] == 1]
    if len(domain1_row) > 0:
        domain1_y_min = float(domain1_row["y_min"].iloc[0])
        slope, intercept = np.polyfit(transects["origin_y"].values,
                                        transects["dist_from_groin_m"].values, 1)
        dist_cutoff_m = slope * domain1_y_min + intercept
        n_before = len(transects)
        transects = transects[transects["dist_from_groin_m"] >= dist_cutoff_m].reset_index(drop=True)
        n_dropped = n_before - len(transects)
        print(f"  Domain 1's true northing boundary ({domain1_y_min:,.0f}) "
              f"maps to {dist_cutoff_m/1000:+.2f} km alongshore (via a "
              f"linear fit across all transects, robust to Cape Point's "
              f"local northing-vs-alongshore disagreement)")
        print(f"  Excluded {n_dropped}/{n_before} transect(s) south of "
              f"that alongshore position -- too close to Cape Point to "
              f"be reliable for this analysis")
        if len(transects) > 0:
            print(f"  → distance-from-groin range after exclusion: "
                  f"{transects['dist_from_groin_m'].min()/1000:+.2f} to "
                  f"{transects['dist_from_groin_m'].max()/1000:+.2f} km")
    else:
        print("  ! domain 1 not found in domain reference -- skipping "
              "Cape Point exclusion")

    # 6. CoastSat chainage -- read directly, keyed by CoastSat's own
    # transect ID (no reconstruction/snapping needed now that CoastSat's
    # transects ARE the backbone).
    coastsat_chainage = load_coastsat_chainage(transects, COASTSAT_ROOT_DIR)

    # 6. Wet-dry & NC state chainage via intersection with the SAME
    # CoastSat transects (short, ~800 m; default TRANSECT_INTERSECTION_
    # LENGTH_M / CHAINAGE_MAX_ABS_M).
    wet_dry_chainage = extract_chainage_by_intersection(
        wet_dry_gdf, transects, "wet_dry")
    nc_state_chainage = extract_chainage_by_intersection(
        nc_state_gdf, transects, "nc_state")

    # Sanity check: how much wet-dry coverage exists on each side of
    # the groin?
    wd_check = wet_dry_chainage.merge(
        transects[["transect_id", "dist_from_groin_m"]], on="transect_id", how="left")
    print(f"\n  [sanity check] wet-dry observations: "
          f"{int((wd_check['dist_from_groin_m'] > 0).sum())} updrift, "
          f"{int((wd_check['dist_from_groin_m'] < 0).sum())} downdrift "
          f"of the groin")

    # Consolidated spatial-reference confirmation for every file loaded
    # so far (study area, groin, domains, CoastSat transects, wet-dry,
    # NC state)
    print_crs_summary()

    # 7. Merge into unified per-(transect, date) table
    chainage = build_unified_chainage_table(
        coastsat_chainage, wet_dry_chainage, nc_state_chainage, transects)
    chainage.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_chainage_all.csv"),
        index=False)
    print(f"→ saved chainage_all.csv ({len(chainage)} rows)")

    # 8. Pre-installation LRR (per transect where computable)
    preinstall_lrr = compute_preinstall_lrr(chainage, transects)
    preinstall_lrr.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_preinstall_lrr.csv"),
        index=False)

    # 9. Regional baseline: connects the sorted per-transect pre-install
    # points directly (linear interpolation), used for the anomaly
    # calculation specifically -- see interp_preinstall_baseline().
    print(f"\n{'='*72}\nRegional baseline (connect-the-dots, not binned/smoothed)")
    baseline_fn = interp_preinstall_baseline(preinstall_lrr)

    # 10. Full-period LRR + piecewise breakpoint
    fullperiod_lrr = compute_fullperiod_lrr(chainage, transects)
    fullperiod_lrr.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_full_period_lrr.csv"),
        index=False)

    # 10b. Post-install-only LRR (1970-present) -- cleaner post-groin
    # rate than the full-period line, since it isn't diluted by any
    # pre-install years mixed into the same transect's regression.
    postinstall_lrr = compute_post_install_lrr(chainage, transects)
    postinstall_lrr.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_post_install_lrr.csv"),
        index=False)

    # 11. Decadal LRR (fixed 10-yr non-overlapping bins, starting 1960)
    decadal_lrr = compute_decadal_lrr(chainage)
    decadal_lrr.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_decadal_lrr.csv"),
        index=False)

    # 11b. Per-era LRR (one rate per transect per era)
    era_lrrs = compute_era_lrrs(chainage, transects,
                                  extra_periods=[PRE_NOURISHMENT_PERIOD])
    era_lrrs.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_era_lrrs.csv"),
        index=False)

    # 12. Decadal anomaly (LRR minus the interpolated baseline)
    # with distance-from-groin
    decadal_anomaly = compute_decadal_anomaly(decadal_lrr, transects,
                                                baseline_fn)
    decadal_anomaly.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_decadal_anomaly.csv"),
        index=False)

    # 13. Signal extent over time
    signal_extent = compute_signal_extent_over_time(decadal_anomaly)
    signal_extent.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_signal_extent.csv"),
        index=False)

    # 13b. Distance-band time series (curve replacement for the original
    # transect x year heatmaps)
    lrr_band_series = compute_distance_band_series(decadal_lrr.merge(
        transects[["transect_id", "dist_from_groin_m"]]
            .assign(transect_id=lambda d: d["transect_id"].astype(str)),
        on="transect_id", how="inner"), "slope_m_yr")
    lrr_band_series.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_lrr_distance_band_series.csv"),
        index=False)
    anomaly_band_series = compute_distance_band_series(
        decadal_anomaly, "anomaly_m_yr")
    anomaly_band_series.to_csv(
        os.path.join(OUTPUT_DIR, "groin_analysis_anomaly_distance_band_series.csv"),
        index=False)

    # 14. Plots
    # Per your feedback: only the two profile plots below were actually
    # useful. The distance-band, anomaly-profile, signal-extent, and
    # transect-diagnostic PLOTS have been dropped (their underlying CSVs
    # from steps 12-13b are still written, in case the numbers are
    # useful later even without a chart).
    print("\n" + "=" * 72)
    print("GENERATING PLOTS")
    print("=" * 72)
    plot_alongshore_lrr_profile(
        preinstall_lrr, fullperiod_lrr, postinstall_lrr, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_alongshore_profile.png"))
    plot_era_lrr_profile(
        era_lrrs, preinstall_lrr, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_era_lrr_profile.png"))
    plot_era_difference_profile(
        era_lrrs, preinstall_lrr, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_era_difference_profile.png"))
    plot_pre_to_functional_difference(
        era_lrrs, preinstall_lrr, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_diff_pre_to_functional.png"))
    plot_functional_to_deteriorated_difference(
        era_lrrs, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_diff_functional_to_deteriorated.png"))
    plot_deteriorated_to_preinstall_difference(
        era_lrrs, preinstall_lrr, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_diff_deteriorated_to_preinstall.png"))

    # Decade-increment LRR evolution: how has the groin's effect
    # changed continuously since installation, not just stepped
    # between two states? One plot per increment in
    # DECADE_PLOT_INCREMENTS_YEARS (e.g. both 5-yr and 10-yr).
    # Reuses compute_era_lrrs() with the decade windows as
    # extra_periods (the 3 official eras get computed alongside but
    # aren't used by this plot -- harmless, avoids a second near-
    # duplicate per-transect loop function).
    for increment_years in DECADE_PLOT_INCREMENTS_YEARS:
        decade_plot_periods = _build_decade_periods(
            DECADE_PLOT_START_YEAR, increment_years, DECADE_PLOT_END_YEAR)
        decade_period_lrrs = compute_era_lrrs(chainage, transects,
                                               extra_periods=decade_plot_periods)
        decade_period_lrrs.to_csv(
            os.path.join(OUTPUT_DIR,
                          f"groin_analysis_decade_increment_lrrs_{increment_years}yr.csv"),
            index=False)
        plot_decade_lrr_profile(
            decade_period_lrrs, transects, groin, decade_plot_periods,
            increment_years,
            os.path.join(OUTPUT_DIR,
                          f"groin_analysis_decade_increment_profile_{increment_years}yr.png"))

    # Zone-panel profile: one stacked panel per era (Pre-install,
    # Functional groin, Deteriorated), with the downdrift/updrift
    # analysis zones shaded and a mean-rate label per zone per era --
    # adapted from HAT_groin_zone_investigation.py's
    # alongshore_profile_panels.png. Also inserts a panel for each
    # early-functional-era window in
    # ZONE_PANEL_FUNCTIONAL_EARLY_WINDOWS_YEARS (e.g. "first 5yr",
    # "first 10yr" of the groin's functional life), to see whether it
    # looked most effective early on.
    early_functional_periods = [
        (f"Functional groin (first {n}yr)", GROIN_INSTALLATION_YEAR,
         GROIN_INSTALLATION_YEAR + n - 1)
        for n in ZONE_PANEL_FUNCTIONAL_EARLY_WINDOWS_YEARS
    ]
    early_functional_lrrs = compute_era_lrrs(chainage, transects,
                                              extra_periods=early_functional_periods)
    extra_functional_panels = [
        (name, early_functional_lrrs[early_functional_lrrs["era"] == name])
        for name, _, _ in early_functional_periods
    ]
    plot_zone_panels(
        preinstall_lrr, era_lrrs, transects, groin,
        os.path.join(OUTPUT_DIR, "groin_analysis_zone_panels.png"),
        extra_functional_panels=extra_functional_panels)

    # 14b. Shoreline evolution GIF -- re-enabled: now shaded as water/
    # island around the actual shoreline curve, with a y-range fixed to
    # the full historical extent so nothing goes off-frame at any point
    # in the animation.
    grid100m = (load_100m_grid_transects(domain_table)
                if GIF_TRANSECT_SOURCE in ("hybrid", "grid100m") else None)
    create_groin_evolution_gif(chainage, transects, groin, OUTPUT_DIR,
                                domain_table=domain_table, grid100m=grid100m)
    create_groin_evolution_gif(chainage, transects, groin, OUTPUT_DIR,
                                domain_table=domain_table, grid100m=grid100m,
                                window_domain_max=GIF_ZOOM_WINDOW_DOMAIN_MAX,
                                output_suffix="_zoomed")

    # 15. Metadata
    stats = {
        "groin_reference_northing (northernmost groin)": groin["reference_y"],
        "groin_field_extent_min":  groin["y_min"],
        "groin_field_extent_max":  groin["y_max"],
        "n_groin_features":        groin["n_features"],
        "n_coastsat_transects":    len(transects),
        "n_transects_analyzed":    chainage["transect_id"].nunique(),
        "n_total_observations":    len(chainage),
        "n_wetdry_obs_updrift":    int((wd_check["dist_from_groin_m"] > 0).sum()),
        "n_wetdry_obs_downdrift":  int((wd_check["dist_from_groin_m"] < 0).sum()),
        "n_preinstall_transects":  int((preinstall_lrr["n_obs"] >= 3).sum()),
        "preinstall_years_used":   preinstall_lrr.attrs.get("preinstall_years", []),
        "n_piecewise_bp_found":    int(fullperiod_lrr["breakpoint_year"].notna().sum()),
        "median_pre_install_LRR_m_yr":
            float(preinstall_lrr["slope_m_yr"].median()),
        "median_full_period_LRR_m_yr":
            float(fullperiod_lrr["slope_full"].median()),
        "median_updrift_extent_km":
            float(signal_extent["updrift_extent_m"].median() / 1000),
        "median_downdrift_extent_km":
            float(signal_extent["downdrift_extent_m"].median() / 1000),
        "max_updrift_peak_anomaly_m_yr":
            float(signal_extent["updrift_peak_anomaly"].max()),
        "max_downdrift_peak_anomaly_m_yr":
            float(signal_extent["downdrift_peak_anomaly"].min()),
    }
    write_metadata(
        os.path.join(OUTPUT_DIR, "groin_analysis_metadata.txt"), stats)

    print("\n" + "=" * 72)
    print("DONE.")
    print("=" * 72)


if __name__ == "__main__":
    main()
