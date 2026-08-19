"""CoastSat transect loading and LOESS smoothing for the rate-comparison figures.

Transect-level LRR (linear regression rate) values are loaded from
transect_lrr_full.csv, LOESS-smoothed at transect resolution (physical
along-coast distance as x), then aggregated to GIS-domain resolution for
comparison against the CASCADE model. This module only computes -- nothing
here touches matplotlib; see cascade_pipeline.plotting.rate_comparison for the
figures that consume its output.
"""

import dataclasses
import os

import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

from cascade_pipeline.domains import DEFAULT_DOMAINS


@dataclasses.dataclass(frozen=True)
class CoastSatDataset:
    """One CoastSat transect CSV to load and smooth.

    Attributes:
        label: Legend label, e.g. "CoastSat LRR (1984-2004)".
        period_start: Calendar year this dataset's period begins; compared
            against the active run's start year to decide solid vs faded
            styling downstream (build_coastsat_series' active flag).
        csv_path: Path to transect_lrr_full.csv.
        domain_col: Column holding the GIS domain ID.
        rate_col: Column holding the LRR value (m/yr).
        transect_id_col: Column holding a per-transect ID, used to keep
            each domain's transects in a stable order.
    """

    label: str
    period_start: int
    csv_path: str
    domain_col: str = "domain_number"
    rate_col: str = "lrr_m_yr"
    transect_id_col: str = "transect_id"


@dataclasses.dataclass(frozen=True)
class LoessConfig:
    """LOESS smoothing settings shared by every CoastSat dataset.

    Attributes:
        window_domains: One or two window widths, in domain units
            (1 domain ~= 500 m). The largest is treated as the primary
            reference window wherever a single window is needed.
        skip_southern_domains: Domains 1..N shown as raw per-domain means
            instead of LOESS-smoothed -- boundary effects near Oregon Inlet
            dominate this zone and smoothing can obscure the sharp
            gradient there. 0 disables the splice (LOESS used everywhere).
    """

    window_domains: tuple = (7, 10)
    skip_southern_domains: int = 10


DEFAULT_LOESS = LoessConfig()


def estimate_transect_spacing(along_coast_m):
    """Median spacing between consecutive transects, in meters.

    Args:
        along_coast_m: Along-coast distance for each transect.

    Returns:
        Median of the positive consecutive differences, or 50.0 if there
        are fewer than two distinct positions.
    """
    arr = np.sort(along_coast_m)
    diffs = np.diff(arr)
    pos = diffs[diffs > 0]
    return float(np.median(pos)) if len(pos) else 50.0


def load_transect_data(dataset, domains=DEFAULT_DOMAINS):
    """Load individual transect LRR values for one CoastSatDataset.

    Along-coast distance is derived by spreading each domain's transects
    evenly across its domain_spacing_m band.

    Args:
        dataset: A CoastSatDataset.
        domains: DomainGeometry; first_gis_id/last_gis_id/domain_spacing_m
            are used to filter and space transects.

    Returns:
        (domain_ids, lrr_values, along_coast_m): int/float arrays, or
        (None, None, None) if the CSV doesn't exist.
    """
    if not os.path.exists(dataset.csv_path):
        print(f"  WARNING: Transect CSV not found: {dataset.csv_path}")
        return None, None, None

    df = pd.read_csv(dataset.csv_path)
    df.columns = [c.split(".csv")[-1] if ".csv" in c else c for c in df.columns]

    domain_col, rate_col, id_col = dataset.domain_col, dataset.rate_col, dataset.transect_id_col
    for col in (domain_col, rate_col):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna(subset=[domain_col, rate_col])
    df[domain_col] = df[domain_col].astype(int)
    df = df[(df[domain_col] >= domains.first_gis_id) & (df[domain_col] <= domains.last_gis_id)]

    sort_cols = [domain_col, id_col] if id_col in df.columns else [domain_col]
    df = df.sort_values(sort_cols).reset_index(drop=True)

    # Spread each domain's transects evenly across its domain_spacing_m band.
    df["_rank"] = df.groupby(domain_col).cumcount()
    df["_n"] = df.groupby(domain_col)[domain_col].transform("count")
    df["along_coast_m"] = (
        (df[domain_col] - 1) * domains.domain_spacing_m
        + (df["_rank"] + 0.5) * (domains.domain_spacing_m / df["_n"])
    )
    df = df.drop(columns=["_rank", "_n"])

    domain_ids = df[domain_col].values.astype(int)
    lrr_values = df[rate_col].values.astype(float)
    along_coast_m = df["along_coast_m"].values.astype(float)

    spacing = estimate_transect_spacing(along_coast_m)
    print(f"  {dataset.label}: {len(df)} transects  "
          f"est. spacing {spacing:.0f} m  "
          f"LRR range {np.nanmin(lrr_values):+.2f}-{np.nanmax(lrr_values):+.2f} m/yr")
    return domain_ids, lrr_values, along_coast_m


def loess_smooth_transect_to_domains(along_coast_m, lrr, domain_ids, window_domains,
                                      domains=DEFAULT_DOMAINS):
    """Apply LOESS at transect resolution, then aggregate to domain resolution.

    Args:
        along_coast_m, lrr, domain_ids: Output of load_transect_data.
        window_domains: LOESS window width, in domain units.
        domains: DomainGeometry; only domain_spacing_m is used.

    Returns:
        (gis_x, smoothed, frac): GIS domain IDs with at least one transect,
        the domain-averaged smoothed LRR (m/yr), and the LOESS frac used
        (for logging). (None, None, frac) if fewer than 5 valid transects.
    """
    window_km = window_domains * domains.domain_spacing_m / 1000.0
    spacing_m = estimate_transect_spacing(along_coast_m)
    n = len(along_coast_m)
    frac = float(np.clip((window_km * 1000.0 / spacing_m) / n, 0.02, 1.0))

    valid = np.isfinite(lrr)
    if valid.sum() < 5:
        print(f"  WARNING: Too few valid transects ({valid.sum()}) - skipping LOESS")
        return None, None, frac

    result = lowess(lrr[valid], along_coast_m[valid], frac=frac, return_sorted=True)
    smoothed_t = np.full(n, np.nan)
    smoothed_t[valid] = np.interp(along_coast_m[valid], result[:, 0], result[:, 1])

    dom_agg = (pd.DataFrame({"domain": domain_ids, "smoothed": smoothed_t})
               .groupby("domain")["smoothed"].mean()
               .dropna())

    return dom_agg.index.values.astype(int), dom_agg.values, frac


def compute_domain_means(domain_ids, lrr_values, gis_min, gis_max):
    """Mean LRR per GIS domain within [gis_min, gis_max].

    Used to substitute raw per-domain averages for LOESS smoothing in the
    southernmost domains, where boundary effects dominate.

    Args:
        domain_ids, lrr_values: Per-transect arrays (e.g. from
            load_transect_data).
        gis_min, gis_max: Inclusive GIS domain ID range.

    Returns:
        (gis_x, means): GIS domain IDs with at least one transect, mean
        LRR per domain.
    """
    df = pd.DataFrame({"domain": domain_ids, "lrr": lrr_values})
    sub = df[(df["domain"] >= gis_min) & (df["domain"] <= gis_max) & df["lrr"].notna()]
    if len(sub) == 0:
        return np.array([], dtype=int), np.array([], dtype=float)
    agg = sub.groupby("domain")["lrr"].mean()
    return agg.index.values.astype(int), agg.values


def splice_loess_with_raw_south(win_gis_x, win_smoothed,
                                 transect_domain_ids, transect_lrr_values,
                                 skip_n=DEFAULT_LOESS.skip_southern_domains,
                                 is_widest_window=False):
    """Return (plot_x, plot_y) for one LOESS window, LOESS line starting at skip_n+1.

    Domains 1..skip_n are omitted from the returned line -- they're shown as
    raw scatter only, no smoothed line, since boundary effects near Oregon
    Inlet dominate that zone.

    Args:
        win_gis_x, win_smoothed: Output of loess_smooth_transect_to_domains.
        transect_domain_ids, transect_lrr_values: Accepted for interface
            stability with callers that pass per-window transect data; not
            read by the current splice.
        skip_n: Domains 1..skip_n excluded from the LOESS line.
        is_widest_window: Accepted for interface stability; not read by the
            current splice.

    Returns:
        (plot_x, plot_y): arrays to plot, LOESS domains > skip_n only.
    """
    if skip_n == 0:
        return win_gis_x, win_smoothed
    mask = win_gis_x > skip_n
    return win_gis_x[mask], win_smoothed[mask]


def build_coastsat_series(datasets, active_period_start, loess_config=DEFAULT_LOESS,
                           domains=DEFAULT_DOMAINS):
    """Load every CoastSat dataset and LOESS-smooth it at each configured window.

    Replaces the per-dataset loop previously inlined in main(): loads
    transects, applies loess_smooth_transect_to_domains at each window in
    loess_config.window_domains, and tags each dataset active/reference by
    comparing its period_start to active_period_start (the run's START_YEAR).

    Args:
        datasets: Sequence of CoastSatDataset.
        active_period_start: The run's START_YEAR; datasets with a matching
            period_start are drawn solid/full-opacity downstream.
        loess_config: LoessConfig.
        domains: DomainGeometry.

    Returns:
        List of dicts, one per dataset that loaded successfully:
            label, period_start, active (bool), transect_domains,
            transect_rates, transect_along_coast,
            windows: list of dicts (window, gis_x, smoothed, frac).
    """
    series = []
    for ds in datasets:
        domain_ids, lrr_values, along_coast_m = load_transect_data(ds, domains=domains)
        if domain_ids is None:
            continue
        windows = []
        for w in loess_config.window_domains:
            gis_x, smoothed, frac = loess_smooth_transect_to_domains(
                along_coast_m, lrr_values, domain_ids, w, domains=domains
            )
            if gis_x is None:
                continue
            print(f"  LOESS applied: window={w} domains "
                  f"({w * domains.domain_spacing_m / 1000.0:.1f} km)  "
                  f"frac={frac:.3f}  ({ds.label})")
            windows.append(dict(window=w, gis_x=gis_x, smoothed=smoothed, frac=frac))
        series.append(dict(
            label=ds.label,
            period_start=ds.period_start,
            active=(ds.period_start == active_period_start),
            transect_domains=domain_ids,
            transect_rates=lrr_values,
            transect_along_coast=along_coast_m,
            windows=windows,
        ))
    return series
