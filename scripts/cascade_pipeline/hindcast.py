#!/usr/bin/env python3
"""Shared machinery for the Hatteras hindcast: loaders, diagnostics, runner.

WHY THIS MODULE EXISTS
    `HAT_hindcast_1984_2024.ipynb`, its headless mirror
    `HAT_hindcast_1984_2024.py`, and `HAT_groin_sweep_worker.py` each carried
    their own copy of these functions -- 667 lines duplicated character for
    character between the notebook and the .py, and a third copy of
    `build_cascade` in the sweep worker that the worker's own comment warned
    "can drift". It had already drifted in the docstrings.

    They live here now, defined once. The notebook and the .py keep everything
    that describes THIS study -- the paths, the switches, the reports, the
    figures -- and import the machinery that is the same either way.

WHAT IS NOT HERE
    Nothing that decides what is simulated. Every run-selecting value still
    comes from `HAT_hindcast_config` (which reads `hat_run.yaml`, overridden by
    the environment) and every path still comes from section 2 of the calling
    file, so a run remains fully described by the settings file plus the file
    you can read top to bottom. The functions below take those
    values as arguments rather than reaching for module globals, which is the
    only change made to any body while moving it.

THE SYNC RULE
    Editing a function here changes the notebook, the .py, and the sweep at
    once. That is the point. Verify with a full hindcast run, not by reading.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# True: sandbox copy with the pre-AST groin hook. False: real package, hook
# folded in. Resolved here so the notebook, the .py and the sweep worker cannot
# disagree about which Cascade they built.
#
# Read from the environment rather than from `HAT_hindcast_config`, which is
# the deliberate direction of the dependency: this package is shared, and the
# Hatteras settings file is not something it should know about. The runner
# reads `use_sandbox_cascade` from hat_run.yaml and exports it as
# CASCADE_USE_SANDBOX in its section 1, BEFORE this module is imported -- the
# choice selects an import, so it cannot be made after the fact.
USE_SANDBOX_CASCADE = os.environ.get(
    "CASCADE_USE_SANDBOX", "1").strip().lower() not in (
        "0", "false", "no", "off")
if USE_SANDBOX_CASCADE:
    from cascade.cascade_groin import Cascade
else:
    from cascade.cascade import Cascade

# scripts/ is the parent of this package, so the resolver is importable
# without a path hack. It owns the dune-topo directory AND the array names.
from hat_topo_version import domain_arrays

from cascade_pipeline import roadway as roadway_module
from cascade_pipeline.coastsat_loess import compute_domain_means

__all__ = [
    "USE_SANDBOX_CASCADE", "Cascade",
    "DAM_TO_M", "WATER_CLAMP_DAM", "MAX_PLAUSIBLE_DAM",
    "STORM_COLUMNS", "RAW_OFFSET_COLUMNS",
    "build_domain_file_paths", "load_barrier3d_contract", "check_domain_units",
    "load_island_offset_dam", "build_island_offset",
    "island_offset_tilts", "ISLAND_OFFSET_MODES", "load_storm_series", "build_background_erosion",
    "groin_trapping_schedule", "implied_interception_m3_yr",
    "build_target_table", "groin_differential", "scenario_run_name",
    "load_absolute_dune_distance", "build_shoreline_target",
    "build_cascade", "run_cascade_simulation", "brie_r_ipl",
    "measure_groin_extent",
]

# --- unit and file-format constants ------------------------------------------
# Fixed by Barrier3D's input contract and by the extractor's RUN_MANIFEST, not
# by the scenario, so they are the same for every run and every period.
DAM_TO_M = 10.0          # Barrier3D works in decameters
WATER_CLAMP_DAM = -0.3   # SENTINEL_WATER_M / WATER_CLAMP_M from RUN_MANIFEST.txt
MAX_PLAUSIBLE_DAM = 2.0  # 20 m; a barrier island in metres would blow past this

STORM_COLUMNS = ("time", "Rhigh", "Rlow", "period", "duration")
RAW_OFFSET_COLUMNS = dict(domain="domain_id", distance="ORIG_LEN",
                          transect="LineID")


# =============================================================================
# INPUT FILES AND THE UNITS CONTRACT
# =============================================================================

def build_domain_file_paths(geometry, product=None, override=None):
    """Builds one elevation and one dune file path per padded domain.

    Delegates to hat_topo_version.domain_arrays(), which owns BOTH the
    directory and the filename. This function used to spell the names itself -
    f"domain_{gis_id}_topography_{init_year}.npy" - and that body was copied
    verbatim into the groin sweep worker and the notebook, each with its own
    `init_year`. On 2026-08-26 the year was dropped from the names entirely and
    the three copies collapsed into one call; see the note at the top of
    hat_topo_version.py for why a per-period tag was tried and reverted.

    Args:
        geometry: DomainGeometry describing the padded domain array.
        product: topography product, e.g. "1984-start". None resolves the
            resolver's DEFAULT_PRODUCT.
        override: pin a dune-topo version instead of resolving it.

    Returns:
        An (elevation_paths, dune_paths) tuple of string lists. Each list is
        geometry.total_domains long and index-aligned with the padded array.
    """
    return domain_arrays(
        product=product,
        override=override,
        first_gis=geometry.first_gis_id,
        last_gis=geometry.last_gis_id,
        n_buffer=geometry.num_buffer_domains,
    )


def load_barrier3d_contract(parameter_path):
    """Reads the unit-relevant Barrier3D parameters from the CASCADE YAML.

    Mirrors the conversions in barrier3d/load_input.py so the values can be
    compared against the raw arrays on disk.

    Args:
        parameter_path: Path to the CASCADE parameter YAML.

    Returns:
        A dict with barrier_length_cells, mhw_dam and berm_el_dam.
    """
    with open(parameter_path, encoding="utf-8") as handle:
        params = yaml.safe_load(handle)

    mhw_dam = params["MHW"] / 10.0
    return {
        "barrier_length_cells": int(params["BarrierLength"] / 10.0),
        "mhw_dam": mhw_dam,
        "berm_el_dam": params["BermEl"] / 10.0 - mhw_dam,
    }


def check_domain_units(elevation_dam, dune_dam, contract):
    """Checks one domain's raw arrays against the Barrier3D input contract.

    Args:
        elevation_dam: Raw elevation array as stored on disk.
        dune_dam: Raw dune array as stored on disk.
        contract: Mapping from load_barrier3d_contract.

    Returns:
        A dict of check name -> bool, True when the array satisfies the check.
    """
    alongshore_cells = (elevation_dam.shape[1] if elevation_dam.ndim == 2
                        else None)
    return {
        "elevation is 2-D (cross_shore, alongshore)": elevation_dam.ndim == 2,
        "dune is 1-D (one crest per alongshore cell)": dune_dam.ndim == 1,
        "alongshore cells match dune length":
            alongshore_cells == dune_dam.size,
        f"alongshore >= BarrierLength ({contract['barrier_length_cells']})":
            alongshore_cells is not None
            and alongshore_cells >= contract["barrier_length_cells"],
        "dtype is floating point":
            np.issubdtype(elevation_dam.dtype, np.floating)
            and np.issubdtype(dune_dam.dtype, np.floating),
        "magnitudes are decameters, not metres":
            np.abs(elevation_dam).max() <= MAX_PLAUSIBLE_DAM,
        f"seaward floor at water clamp ({WATER_CLAMP_DAM} dam)":
            np.isclose(elevation_dam.min(), WATER_CLAMP_DAM),
        "dune heights positive (above berm)": dune_dam.min() > 0,
    }


# =============================================================================
# PERIOD FORCINGS
# =============================================================================

def load_island_offset_dam(offset_path, geometry):
    """Loads a padded BRIE island-offset file and converts it to decameters.

    Args:
        offset_path: Path to a single-column padded offset CSV, in meters.
        geometry: DomainGeometry the file must match in length.

    Returns:
        A 1-D array of offsets in decameters, one per padded domain, ready to
        pass to Cascade() as shoreline_offset.

    Raises:
        ValueError: If the file is not one value per padded domain.
    """
    offset_m = np.loadtxt(offset_path, skiprows=1, delimiter=",")
    if offset_m.ndim != 1 or offset_m.size != geometry.total_domains:
        raise ValueError(f"{offset_path.name}: expected "
                         f"{geometry.total_domains} values, got shape "
                         f"{offset_m.shape}")
    return offset_m / DAM_TO_M


ISLAND_OFFSET_MODES = ("asrun", "metres", "detrended")


def _close_offset_ring(real_m, buffers):
    """Buffer values carrying the shoreline from real_m[-1] back to real_m[0].

    BRIE's alongshore domain is PERIODIC -- `x_s[np.r_[1:ny, 0]]` makes the
    last padded domain a neighbour of the first -- so the padded offset has to
    come back to where it started. Hatteras does not: it runs 5.4 km across
    the shore from Cape Point to Oregon Inlet and keeps going. The buffers are
    the invented coast that closes the loop.

    A cubic Hermite is used, with end tangents matched to the island's own
    local slopes, so there is no kink where buffer meets real domain. The path
    runs real[-1] -> right buffer -> (wrap) -> left buffer -> real[0], which is
    2 * buffers + 1 steps.

    Why it matters that this be smooth: BRIE's diffusivity term
    `(1.2 sin^2 - cos^2)` changes SIGN near 42 degrees. Below that a shoreline
    smooths itself; above it, bumps GROW (high-angle instability). The original
    linear-bridge padding leaves a 1154 m step at the seam, which is 66.6
    degrees once the offset is in its correct unit -- unstable. It went
    unnoticed because the decameter bug divided it to a harmless 13 degrees.

    Args:
        real_m: The real-domain offsets, in meters.
        buffers: Buffer domains per side.

    Returns:
        A (right_buffer, left_buffer) tuple of arrays, each `buffers` long.
    """
    n = 2 * buffers + 1
    p0, p1 = real_m[-1], real_m[0]
    m0 = (real_m[-1] - real_m[-2]) * n
    m1 = (real_m[1] - real_m[0]) * n
    t = np.arange(1, n) / n
    path = ((2 * t ** 3 - 3 * t ** 2 + 1) * p0
            + (t ** 3 - 2 * t ** 2 + t) * m0
            + (-2 * t ** 3 + 3 * t ** 2) * p1
            + (t ** 3 - t ** 2) * m1)
    return path[:buffers], path[buffers:]


def build_island_offset(offset_path, geometry, mode="asrun"):
    """Builds the `shoreline_offset` array Cascade is handed.

    UNITS. `Cascade(shoreline_offset=...)` must be in METRES.
    `brie_coupler.offset_shoreline` adds the values straight onto
    `brie.x_t` / `brie.x_s` with no conversion, and those are metres -- see
    `brie_coupler.py:390` ("convert from dam to meters") and line 344
    (`barrier3d.x_s = brie.x_s / 10`). The measurement file is already metres,
    so it needs no conversion at all.

    Args:
        offset_path: Padded offset CSV, one value per padded domain, meters.
        geometry: DomainGeometry the file must match in length.
        mode: Which variant to build.
            "asrun"     - meters / 10, reproducing the historical unit error
                          exactly. Kept so earlier runs stay reproducible.
            "metres"    - the measurement as-is, with the ring re-closed.
                          Carries the island's full planform including its
                          8 degree lean.
            "detrended" - the measurement with its linear trend removed, ring
                          re-closed. Keeps all 2052 m of real curvature and
                          drops the lean, which is what lets the periodic
                          domain close without a steep buffer.

    Returns:
        A 1-D array, one value per padded domain, ready to pass to Cascade.

    Raises:
        ValueError: If the file length does not match the geometry, or the
            mode is unknown.
    """
    if mode not in ISLAND_OFFSET_MODES:
        raise ValueError(f"offset mode {mode!r} must be one of "
                         f"{list(ISLAND_OFFSET_MODES)}")
    offset_m = np.loadtxt(offset_path, skiprows=1, delimiter=",")
    if offset_m.ndim != 1 or offset_m.size != geometry.total_domains:
        raise ValueError(f"{offset_path.name}: expected "
                         f"{geometry.total_domains} values, got shape "
                         f"{offset_m.shape}")
    if mode == "asrun":
        return offset_m / DAM_TO_M

    lo, hi = geometry.start_real_index, geometry.end_real_index
    real_m = offset_m[lo:hi]
    if mode == "detrended":
        x = np.arange(real_m.size) * geometry.domain_spacing_m
        real_m = real_m - np.polyval(np.polyfit(x, real_m, 1), x)

    right, left = _close_offset_ring(real_m, geometry.num_buffer_domains)
    return np.concatenate([left, real_m, right])


def island_offset_tilts(offset, geometry):
    """Local shoreline tilts BRIE will read from an offset array, in degrees.

    Reproduces `brie.py:820`, `theta = atan2(diff(x_s), dy)`, including the
    periodic wrap, so a run can report the angles it is about to impose rather
    than leaving them to be discovered in the output.

    Args:
        offset: The padded array passed to Cascade, meters.
        geometry: DomainGeometry.

    Returns:
        A dict with island_max_deg, buffer_max_deg, seam_deg and unstable.
    """
    theta = np.degrees(np.arctan2(np.diff(np.r_[offset, offset[0]]),
                                  geometry.domain_spacing_m))
    lo, hi = geometry.start_real_index, geometry.end_real_index
    island = float(np.abs(theta[lo:hi - 1]).max())
    buffer_max = float(np.abs(np.r_[theta[:lo], theta[hi - 1:]]).max())
    return {"island_max_deg": island, "buffer_max_deg": buffer_max,
            "seam_deg": float(theta[-1]),
            # BRIE's (1.2 sin^2 - cos^2) changes sign here; above it the
            # shoreline is anti-diffusive.
            "unstable": bool(max(island, buffer_max) > 42.0)}


def load_storm_series(storm_path):
    """Loads a Barrier3D storm series into a DataFrame.

    Args:
        storm_path: Path to the storm .npy file.

    Returns:
        A DataFrame with the raw columns plus Rhigh_m and Rlow_m in meters.

    Raises:
        ValueError: If the file does not have the expected column count.
    """
    storms = np.load(storm_path)
    if storms.ndim != 2 or storms.shape[1] != len(STORM_COLUMNS):
        raise ValueError(f"{storm_path.name}: expected "
                         f"(n, {len(STORM_COLUMNS)}), got {storms.shape}")

    series = pd.DataFrame(storms, columns=STORM_COLUMNS)
    series["Rhigh_m"] = series["Rhigh"] * DAM_TO_M
    series["Rlow_m"] = series["Rlow"] * DAM_TO_M
    return series


def build_background_erosion(be_rates, geometry):
    """Expands sparse per-GIS background-erosion rates onto the padded array.

    Args:
        be_rates: Mapping of GIS domain ID to rate in m/yr. Domains absent
            from the mapping get 0.0.
        geometry: DomainGeometry describing the padded array.

    Returns:
        A list of geometry.total_domains rates, ready to pass to Cascade().

    Raises:
        ValueError: If a GIS ID falls outside the padded array.
    """
    rates = [0.0] * geometry.total_domains
    for gis_id, rate in be_rates.items():
        pad_index = geometry.gis_to_pad(gis_id)
        if not 0 <= pad_index < geometry.total_domains:
            raise ValueError(f"GIS {gis_id} -> pad index {pad_index}, outside "
                             f"0-{geometry.total_domains - 1}")
        rates[pad_index] = float(rate)
    return rates


# =============================================================================
# GROIN DIAGNOSTICS
# =============================================================================

def groin_trapping_schedule(callback, start_year, end_year):
    """Effective trapping rate for every model year of a period.

    Reads the callback's deterioration curve without running it.
    `_effective_trapping_rate` is a pure function of the year -- it does not
    touch `_call_count` -- so querying it here leaves the callback's year
    counter at zero for the actual run.

    Args:
        callback: A GroinCallback.
        start_year: First model year.
        end_year: Last model year, exclusive.

    Returns:
        A (years, M_eff) tuple of 1-D arrays.
    """
    years = np.arange(start_year, end_year)
    m_eff = np.array([callback._effective_trapping_rate(int(y)) for y in years])
    return years, m_eff


def implied_interception_m3_yr(m_per_yr, profile_height_m, geometry):
    """Sediment volume the dipole transfers across the groin each year.

    A shoreline displacement of `m_per_yr` over one domain implies a volume of
    `m_per_yr * dy * (d_sf + h_b)`. The dipole moves that from the downdrift
    side to the updrift side; the transfer is volume-neutral overall.

    Args:
        m_per_yr: Trapping rate M, in meters per year.
        profile_height_m: Active profile height (d_sf + h_b), in meters.
        geometry: DomainGeometry supplying the alongshore domain width.

    Returns:
        The implied transfer in cubic meters per year.
    """
    return m_per_yr * geometry.domain_spacing_m * profile_height_m


def measure_groin_extent(shoreline_m, baseline_m, geometry, updrift_gis,
                         downdrift_gis, threshold_frac):
    """Alongshore extent of the groin's effect, from a paired baseline run.

    The effect is this run's final shoreline minus the no-groin baseline's,
    sign-flipped so positive means seaward of the baseline. Extent is the
    contiguous run of domains, outward from the groin, where the effect holds
    at or above `threshold_frac` of its own peak magnitude.

    Args:
        shoreline_m: [time x domain] matrix from this run.
        baseline_m: [time x domain] matrix from the no-groin run.
        geometry: DomainGeometry describing the padded array.
        updrift_gis, downdrift_gis: The groin's flanking domains.
        threshold_frac: Fraction of the peak effect defining the edge.

    Returns:
        A dict with effect (padded array, m), peak_m, threshold_m, and the
        updrift/downdrift extents in domains and meters.
    """
    effect = -(np.asarray(shoreline_m)[-1] - np.asarray(baseline_m)[-1])
    peak = float(np.nanmax(np.abs(effect))) if effect.size else 0.0
    threshold = threshold_frac * peak

    def span(start_gis, step):
        # A zero peak means the two runs are identical, so the threshold is
        # also zero and every "abs(value) < 0" test is False -- the walk would
        # run to the end of the grid and report the whole island as fillet.
        if not np.isfinite(peak) or peak <= 0.0:
            return 0
        count, gis = 0, start_gis
        while geometry.first_gis_id <= gis <= geometry.last_gis_id:
            value = effect[geometry.gis_to_pad(gis)]
            if not np.isfinite(value) or abs(value) < threshold:
                break
            count += 1
            gis += step
        return count

    up = span(updrift_gis, +1)
    down = span(downdrift_gis, -1)
    return dict(
        effect=effect, peak_m=peak, threshold_m=threshold,
        updrift_domains=up, updrift_m=up * geometry.domain_spacing_m,
        downdrift_domains=down, downdrift_m=down * geometry.domain_spacing_m,
    )


# =============================================================================
# COASTSAT TARGETS AND THE SURVEYED END POSITION
# =============================================================================

def build_target_table(cs, loess_config, geometry, window):
    """Per-domain target rate, labelled with where each value came from.

    The target is not one curve: GIS 1..skip_southern_domains are raw
    per-domain means (LOESS is suppressed there), and the rest is the
    LOESS-smoothed reference window. Every row records which.

    Args:
        cs: One entry from build_coastsat_series.
        loess_config: LoessConfig used to build it.
        geometry: DomainGeometry describing the real-domain span.
        window: Reference window width, in domains.

    Returns:
        A DataFrame with gis_domain, target_lrr_m_yr, source, n_transects.

    Raises:
        ValueError: If `window` was not computed for this dataset.
    """
    match = [w for w in cs["windows"] if w["window"] == window]
    if not match:
        raise ValueError(
            f"window {window} not in {[w['window'] for w in cs['windows']]}")
    smoothed = dict(zip(match[0]["gis_x"], match[0]["smoothed"]))

    skip = loess_config.skip_southern_domains
    raw_x, raw_y = compute_domain_means(
        cs["transect_domains"], cs["transect_rates"],
        geometry.first_gis_id, skip)
    raw = dict(zip(raw_x, raw_y))
    counts = pd.Series(cs["transect_domains"]).value_counts()

    rows = []
    for gis in range(geometry.first_gis_id, geometry.last_gis_id + 1):
        if gis <= skip:
            value, source = raw.get(gis, np.nan), f"raw mean (D1-{skip})"
        else:
            value, source = smoothed.get(gis, np.nan), f"LOESS {window}-dom"
        rows.append(dict(gis_domain=gis, target_lrr_m_yr=value, source=source,
                         n_transects=int(counts.get(gis, 0))))
    return pd.DataFrame(rows)


def groin_differential(cs, updrift_gis, downdrift_gis):
    """Observed updrift-minus-downdrift mean LRR for one CoastSat period.

    The observational check on the dipole section 7 imposes: a positive
    differential means the updrift domain is retreating more slowly (or
    advancing faster) than the downdrift one, which is the signature a
    functioning groin should leave.

    Args:
        cs: One entry from build_coastsat_series.
        updrift_gis: GIS domain on the updrift side.
        downdrift_gis: GIS domain on the downdrift side.

    Returns:
        A dict of label, updrift, downdrift, differential, and per-domain
        transect counts. Rates are NaN where a domain has no transects.
    """
    lo, hi = min(updrift_gis, downdrift_gis), max(updrift_gis, downdrift_gis)
    gis_x, means = compute_domain_means(
        cs["transect_domains"], cs["transect_rates"], lo, hi)
    by_domain = dict(zip(gis_x, means))
    counts = pd.Series(cs["transect_domains"]).value_counts()
    up = by_domain.get(updrift_gis, np.nan)
    down = by_domain.get(downdrift_gis, np.nan)
    return dict(
        label=cs["label"], period_start=cs["period_start"],
        active=cs["active"], updrift=up, downdrift=down,
        differential=up - down,
        n_updrift=int(counts.get(updrift_gis, 0)),
        n_downdrift=int(counts.get(downdrift_gis, 0)),
    )


def load_absolute_dune_distance(year, geometry, raw_dir,
                                columns=RAW_OFFSET_COLUMNS):
    """Mean distance from the offshore datum to the dune line, per GIS domain.

    The absolute quantity the padded offset files are built from, BEFORE
    island_offset_hybrid.py subtracts that year's own minimum. Absolute
    distances share a fixed datum across years, so differencing two of them is
    a real shoreline change; differencing the padded files is not.

    One row per transect (`LineID`) is kept before averaging, mirroring
    island_offset_hybrid.py:108-113, so a transect sampled at many points does
    not outweigh one sampled at few.

    Args:
        year: Survey year; reads <raw_dir>/<year>_duneline_offset_raw.csv.
        geometry: DomainGeometry supplying the real-domain GIS range.
        raw_dir: Directory holding the raw per-year CSVs.
        columns: Mapping with "domain", "distance" and "transect" keys.

    Returns:
        1-D array of length geometry.num_real_domains, in meters, increasing
        LANDWARD, indexed by GIS domain. NaN where a domain has no transects.

    Raises:
        FileNotFoundError: If that year has no raw file.
        KeyError: If an expected column is missing.
    """
    path = raw_dir / f"{year}_duneline_offset_raw.csv"
    raw = pd.read_csv(path)
    missing = [c for c in columns.values() if c not in raw.columns]
    if missing:
        raise KeyError(f"{path.name}: missing column(s) {missing}; "
                       f"have {list(raw.columns)[:12]}...")

    per_transect = raw.drop_duplicates(
        subset=[columns["domain"], columns["transect"]])
    means = per_transect.groupby(columns["domain"])[columns["distance"]].mean()

    gis_ids = np.arange(geometry.first_gis_id, geometry.last_gis_id + 1)
    return means.reindex(gis_ids).to_numpy(dtype=float)


def build_shoreline_target(model_year0_m, start_year, end_year, geometry,
                           raw_dir):
    """Surveyed end-year shoreline position, in the model's own x_s frame.

    Takes the model's year-0 position and adds the OBSERVED start->end change,
    so the two share a base and the gap between them is the model's misfit.
    Returns positions rather than a change because that is what
    make_shoreline_gif's `target_m` expects.

    Args:
        model_year0_m: shoreline_m[0], the run's year-0 position (raw x_s_TS
            convention, meters, increasing landward), one value per PADDED
            domain.
        start_year: Run start year; must have a raw offset file.
        end_year: Run end year; returns None if it has no raw offset file.
        geometry: DomainGeometry.

    Returns:
        (target_m, observed_change_m):
            target_m: padded array of target positions, NaN in the buffers.
            observed_change_m: real-domain change, + = landward.
        Both None if the end year was never surveyed.
    """
    if not (raw_dir / f"{end_year}_duneline_offset_raw.csv").exists():
        return None, None

    start_dist = load_absolute_dune_distance(start_year, geometry, raw_dir)
    end_dist = load_absolute_dune_distance(end_year, geometry, raw_dir)
    observed_change_m = end_dist - start_dist   # + = landward, as x_s_TS is

    target_m = np.full(geometry.total_domains, np.nan)
    real = slice(geometry.start_real_index, geometry.end_real_index)
    target_m[real] = np.asarray(model_year0_m)[real] + observed_change_m
    return target_m, observed_change_m


# =============================================================================
# RUN IDENTITY
# =============================================================================

def scenario_run_name(switches, stem, **overrides):
    """Run name for a variant of the current scenario.

    Rebuilds the name from the same switch tokens section 7.5 used, with named
    switches overridden. Built from the token list rather than by editing the
    name string, so flipping "groin" cannot accidentally match inside "nogroin".

    Args:
        switches: SCENARIO_SWITCHES, as (label, value, token) triples.
        stem: RUN_NAME_STEM, the period prefix.
        **overrides: label -> replacement token. A token of None or "" drops
            that switch from the name, matching how 7.5 omits default states.

    Returns:
        The run name for that variant.

    Raises:
        KeyError: If an override names a switch that does not exist -- a typo
            would otherwise silently produce the unmodified name.
    """
    labels = {label for label, _, _ in switches}
    unknown = set(overrides) - labels
    if unknown:
        raise KeyError(f"no such switch: {sorted(unknown)}; have {sorted(labels)}")
    tokens = [overrides.get(label, token) for label, _, token in switches]
    return f"{stem}_{'_'.join(t for t in tokens if t)}"


# =============================================================================
# BUILD AND RUN
# =============================================================================
# The split is what lets the caller hold a built-but-unstepped Cascade. BRIE's
# diffusivity and the groin's fillet prediction are only meaningful as initial
# conditions, and a prediction printed after the run is not one.
#
# `run_years` is TRANSITIONS, not states. Barrier3D seeds _x_s_TS = [x_s] at
# init and appends one entry per update, so N updates produce N+1 annual states
# spanning N years. The original signature took `nt` and looped `range(nt - 1)`
# with `nt = END_YEAR - START_YEAR`, which ran 19 updates for a 20-year period
# while dividing by 20 -- every rate came out low by 19/20, and storm years 20
# and 21 were never applied. Here time_step_count = run_years + 1 and the loop
# runs exactly run_years updates.
#
# Nourishment goes to cascade.nourishment_volume, which is where CASCADE reads
# it. Writing it onto the manager instead hits the attribute CASCADE overwrites
# one line before the manager reads it, so the fill spends the init default.

def build_cascade(
    run_years, name, storm_file, alongshore_section_count, num_cores,
    rmin, rmax, elevation_file, dune_file,
    dune_design_elevation, dune_minimum_elevation,
    road_ele, road_width, road_setback,
    overwash_filter, overwash_to_dune,
    nourishment_volume, background_erosion,
    roadway_management_on, beach_dune_manager_on,
    sea_level_rise_rate, sea_level_constant,
    sandbag_management_on, sandbag_elevation,
    enable_shoreline_offset, shoreline_offset,
    wave_height, wave_period, wave_asymmetry, wave_angle_high_fraction,
    berm_elevation, MHW, data_base, parameter_file, groin_callback=None,
):
    """Constructs a Cascade and attaches the groin, without stepping it.

    Separate from run_cascade_simulation so section 11 can inspect the model
    before any time step runs: BRIE's diffusivity and the shoreface depth are
    only meaningful as initial conditions, and the fillet prediction they feed
    is only a prediction if it is printed before the run.

    Args:
        run_years: Annual transitions the run will simulate. The model is built
            with time_step_count = run_years + 1, giving run_years + 1 annual
            states. See the section 10 comment on why this is not run_years.
        name: Run name, passed to Cascade.
        storm_file, elevation_file, dune_file: Barrier3D input paths.
        alongshore_section_count: Padded domain count.
        num_cores: Cores for the parallel Barrier3D step.
        rmin, rmax: Dune growth rate bounds, per domain.
        dune_design_elevation: Rebuild target, m MHW. roadway_manager raises
            this to berm + 1.0 m on the first step if it is lower.
        dune_minimum_elevation: Rebuild trigger, m MHW. roadway_manager raises
            this to berm + 0.3 m on the first step if it is lower.
        road_ele, road_width, road_setback: Padded roadway forcing.
        overwash_filter, overwash_to_dune: beach_dune_manager forcing.
        nourishment_volume: Per-domain init volume. Every scheduled year
            overwrites it via the schedule; see section 6.
        background_erosion: Padded source/sink rates.
        roadway_management_on, beach_dune_manager_on: Per-domain module flags.
        sea_level_rise_rate, sea_level_constant: RSLR forcing.
        sandbag_management_on, sandbag_elevation: Sandbag forcing.
        enable_shoreline_offset, shoreline_offset: Island orientation, in dam.
        wave_height, wave_period, wave_asymmetry, wave_angle_high_fraction:
            Wave climate. wave_height also sets BRIE's shoreface depth,
            d_sf = 8.9 * Hs (brie.py:270).
        berm_elevation, MHW: Barrier3D datums.
        data_base: Datadir handed to Cascade -- the directory the
            parameter file and the padded input files are resolved
            against.
        parameter_file: Barrier3D parameter yaml name, resolved by
            CASCADE inside data_base.
        groin_callback: A GroinCallback to attach, or None.

    Returns:
        The constructed Cascade, before any update().
    """
    cascade = Cascade(
        data_base,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file=parameter_file,

        berm_elevation=berm_elevation,
        MHW=MHW,

        wave_height=wave_height,
        wave_period=wave_period,
        wave_asymmetry=wave_asymmetry,
        wave_angle_high_fraction=wave_angle_high_fraction,

        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_rise_constant=sea_level_constant,

        background_erosion=background_erosion,
        alongshore_section_count=alongshore_section_count,

        # run_years transitions need run_years + 1 states. TMAX is set from
        # this (brie_coupler.py:117), and the loop runs exactly run_years
        # updates, so the last write lands on the final valid index.
        time_step_count=run_years + 1,

        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=num_cores,

        roadway_management_module=roadway_management_on,
        beach_nourishment_module=beach_dune_manager_on,
        sandbag_management_on=sandbag_management_on,
        alongshore_transport_module=True,
        community_economics_module=False,

        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,

        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
        sandbag_elevation=sandbag_elevation,

        overwash_filter=overwash_filter,
        overwash_to_dune=overwash_to_dune,

        enable_shoreline_offset=enable_shoreline_offset,
        shoreline_offset=shoreline_offset,

        nourishment_volume=nourishment_volume,
        nourishment_interval=None,
    )

    if groin_callback is not None:
        cascade._groin_callback = groin_callback

    return cascade


def run_cascade_simulation(
    cascade, run_years, name, run_dir, start_year, geometry,
    alongshore_section_count,
    historical_road_events=(), relocations_enabled=True, setback_check=None,
    nourishment_schedule=None, groin_callback=None, progress=None,
    save_model_state=True,
):
    """Steps a built Cascade through its period and writes the run artifacts.

    Args:
        cascade: A Cascade from build_cascade, not yet stepped.
        run_years: Annual transitions to simulate; the loop runs exactly this
            many updates.
        name: Run name, used for output filenames.
        run_dir: Directory for the saved model and logs.
        start_year: Calendar year of the run's first state.
        geometry: DomainGeometry, for GIS <-> pad translation.
        alongshore_section_count: Padded domain count.
        historical_road_events: RelocationEvent / BridgeEvent sequence.
        relocations_enabled: Global toggle for relocation events.
        setback_check: {gis: measured_setback_m} reported beside relocations.
        nourishment_schedule: A NourishmentSchedule, or None for no fills.
        groin_callback: The attached GroinCallback, for diagnostics output.
        progress: Optional tqdm-like object with update() and write().
            When given, the year counter advances the bar and every event
            message is written through it, so the messages scroll above a
            live bar instead of shredding it. None falls back to a plain
            carriage-return counter, so callers without tqdm still work.
        save_model_state: Whether to write the ~160 MB model pickle.
            The caller passes RUN_CONFIG.save_model_state; everything
            else is written either way.

    Returns:
        The same Cascade, after the run.
    """
    nourishment_log = []
    # One emitter for everything the loop says, so the display mechanism is
    # the caller's choice and this function does not care which it is.
    emit = progress.write if progress is not None else print

    for time_step in range(run_years):
        current_year = start_year + time_step

        # --- historical beach nourishment -----------------------------------
        # apply_to_cascade rewrites BOTH nourish_now and nourishment_volume in
        # full every year, so nothing carries over. It writes the volume to
        # cascade.nourishment_volume, which stock Cascade.update() copies into
        # each BeachDuneManager before calling it -- see section 6.
        if nourishment_schedule is not None:
            applied = nourishment_schedule.apply_to_cascade(
                cascade, current_year)
            if applied:
                print(f"\n  -> nourishment {current_year}:")
                for _row in applied:
                    print(f"       GIS {_row['gis']:>3} (pad {_row['pad']:>3})  "
                          f"{_row['volume_m3_per_m']:.1f} m^3/m")
                    nourishment_log.append(dict(
                        run_name=name, time_step=time_step, **_row))
        else:
            cascade.nourish_now = np.zeros(alongshore_section_count)

        # --- historical roadway events --------------------------------------
        for _event in historical_road_events or ():
            if current_year != _event.year:
                continue
            rows = roadway_module.apply_historical_event(
                cascade, _event, geometry,
                relocations_enabled=relocations_enabled,
                setback_check=setback_check)
            if not rows:
                print(f"\n  -> {current_year} event skipped: {_event.note}")
                continue
            print(f"\n  -> {current_year}: {_event.note}")
            for _row in rows:
                if _row["kind"] == "bridge":
                    emit(f"       GIS {_row['gis']:>3} (pad {_row['pad']:>3})  "
                         f"roadway management OFF")
                    continue
                _chk = (f" | 2004 measured {_row['check_m']:.0f} m"
                        if _row["check_m"] is not None else "")
                emit(f"       GIS {_row['gis']:>3} (pad {_row['pad']:>3})  "
                     f"setback {_row['old_setback_m']:.1f} + "
                     f"{_row['displacement_m']:.0f} -> "
                     f"{_row['new_setback_m']:.1f} m{_chk}")
                for _warn in _row["warnings"]:
                    emit(f"                   [warn] {_warn}")
                if _row["warnings"]:
                    emit(f"                   -> applied anyway "
                         f"(prescribed historical event)")

        cascade.update()

        if progress is not None:
            progress.update(1)
        else:
            print(f"\rYear {time_step + 1}/{run_years}", end="", flush=True)

        if getattr(cascade, "b3d_break", False):
            print(f"\nModel stopped at year {time_step + 1} (b3d_break)")
            break

    # --- did the run do what it was configured to do? ------------------------
    _states = len(cascade.barrier3d[0].x_s_TS)
    if _states != run_years + 1:
        print(f"\n  NOTE: {_states} annual states for {run_years} run_years "
              f"(expected {run_years + 1}) -- the run ended early")
    if groin_callback is not None and len(groin_callback.year_TS) == 0:
        print("\n" + "!" * 74)
        print("WARNING: the groin callback was never called. The pre-AST hook")
        print("in cascade_groin.py is missing, so this run is identical to a")
        print("no-groin run despite GROIN_ENABLED being True.")
        print("!" * 74)

    # --- artifacts -----------------------------------------------------------
    os.makedirs(run_dir, exist_ok=True)
    # The .npz model pickle is ~160 MB and is what lets a figure be re-derived
    # without re-running, so it is written by default. Everything downstream
    # (rate CSV, shoreline matrix, metadata) is written regardless, so a run
    # skipped here is still a complete row in run_index.csv -- just one that
    # cannot be re-plotted from state.
    if save_model_state:
        cascade.save(run_dir)
        print(f"\n  saved: {run_dir}")
    else:
        print(f"\n  model state NOT saved (HAT_SAVE_MODEL_STATE=false)")

    if nourishment_log:
        _bn_csv = os.path.join(run_dir, f"{name}_nourishment_log.csv")
        pd.DataFrame(nourishment_log).to_csv(_bn_csv, index=False)
        print(f"  nourishment log ({len(nourishment_log)} events): {_bn_csv}")

    if groin_callback is not None and groin_callback.year_TS:
        _groin_csv = os.path.join(run_dir, f"{name}_groin_diagnostics.csv")
        pd.DataFrame(groin_callback.diagnostics_frame()).to_csv(
            _groin_csv, index=False)
        print(f"  groin diagnostics ({len(groin_callback.year_TS)} yrs): "
              f"{_groin_csv}")

    return cascade


def brie_r_ipl(cascade, theta_deg=0.0):
    """BRIE's diffusion number at the freshly built model's initial state.

    Reproduces `brie.py:1294`, which computes r_ipl as a local variable inside
    update() and never stores it:

        r_ipl = coast_diff[clip(round(90 - theta))] * dt / 2 / dy**2

    `_coast_diff` is the wave-climate-averaged shoreline diffusivity, built once
    in BRIE's __init__ and deleted by its finalize(). Nothing in CASCADE calls
    finalize, but the angle-dependent index changes as the shoreline evolves, so
    this is only the initial-condition value.

    Args:
        cascade: A built Cascade.
        theta_deg: Shoreline angle to evaluate at. 0 is shore-normal, the
            reference used for the fillet scaling.

    Returns:
        The dimensionless diffusion number.
    """
    brie = cascade._brie_coupler._brie
    index = int(np.clip(round(90 - theta_deg), 1, brie._wave_climl))
    return float(brie._coast_diff[index] * brie._dt / 2.0 / brie._dy ** 2)
