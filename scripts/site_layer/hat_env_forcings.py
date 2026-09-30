"""
Where are the sea-level and storm forcings, and the records they come from?

    from site_layer.hat_env_forcings import storm_series_file, DUCK_GAUGE_FILE

Resolves data/hatteras_init/3-env-forcings/ (1-records, 2-rslr, 3-storms) once;
a window is <start>_<end>. Details: scripts/site_layer/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

from pathlib import Path

if __package__ in (None, ""):
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from site_layer.hat_topo_version import INIT_ROOT, init_relpath  # noqa: E402,F401

ENV_ROOT = INIT_ROOT / "3-env-forcings"

RECORDS = ENV_ROOT / "1-records"
WATER_LEVEL_DIR = RECORDS / "water_level"
WATER_LEVEL_CACHE = WATER_LEVEL_DIR / "noaa_cache_8651370"
DUCK_GAUGE_FILE = WATER_LEVEL_DIR / "8651370_DUCK_19840101_20241231_NAVD.csv"
WIS_DIR = RECORDS / "WIS_raw_data"
WIS_FILE = WIS_DIR / "ST63228-Generic_Export-20260427T11T11_48.csv"
STORM_RECORD_DIR = RECORDS / "storm_record"
WAVE_CLIMATE_DIR = RECORDS / "wave_climate_duke"

RSLR_ROOT = ENV_ROOT / "2-rslr"
RSLR_RECORD_DIR = RSLR_ROOT / "record"
RSLR_FITS_DIR = RSLR_ROOT / "fits"
RSLR_FIGURES_DIR = RSLR_ROOT / "figures"
RSLR_RECORD_FILE = RSLR_RECORD_DIR / "duck_8651370_meantrend.csv"
RSLR_RATES_CSV = RSLR_FITS_DIR / "duck_rslr_rates.csv"

STORMS_ROOT = ENV_ROOT / "3-storms"
HINDCAST_STORMS = STORMS_ROOT / "hindcast_storms"
STORM_VALIDATION = STORMS_ROOT / "validation"
STORM_FIGURES = STORMS_ROOT / "figures"
# The 1984-2024 series spliced for the full-span sweep. NOT a period.
SPLICED_1984_2024 = HINDCAST_STORMS / "1984_2024" / "1984_2024_storms_spliced.npy"

ARCHIVE = ENV_ROOT / "archive"
# benton_storms/ and testing_storms/, retired; the storm_check validators still read base_storms/
SUPERSEDED_STORMS = ARCHIVE / "storms_superseded_20260914"


def window_tag(start_year: int, end_year: int) -> str:
    return f"{int(start_year)}_{int(end_year)}"


def storm_window_dir(start_year: int, end_year: int) -> Path:
    """One window's model-facing storm series."""
    return HINDCAST_STORMS / window_tag(start_year, end_year)


# The hindcast's storm series: 24 h grouping, split at >= 12 h below the berm, each event cut to 24 h
DEFAULT_STORM_VARIANT = "v3_split12_trim24"


def storm_series_file(start_year: int, end_year: int,
                      variant: str = DEFAULT_STORM_VARIANT) -> Path:
    """The .npy CASCADE reads for a window: <window>_storms_<variant>.npy.
    v3_split12_trim24 (the default since 2026-09-29): grouped events split at
    >= 12 h below the berm, each trimmed to 24 h around its peak. v3_trim24
    (2026-09-28): trimmed, not split. v3_72: events over 72 h dropped."""
    tag = window_tag(start_year, end_year)
    return storm_window_dir(start_year, end_year) / f"{tag}_storms_{variant}.npy"


def storm_summary_file(start_year: int, end_year: int,
                       variant: str = DEFAULT_STORM_VARIANT) -> Path:
    """The per-storm table beside the .npy: dates, Rhigh, Rlow, duration."""
    f = storm_series_file(start_year, end_year, variant)
    return f.with_name(f.stem + "_summary.csv")


def storm_validation_dir(start_year: int, end_year: int) -> Path:
    """Where both storm validators write for one window."""
    return STORM_VALIDATION / window_tag(start_year, end_year)


if __name__ == "__main__":
    for name in ("RECORDS", "WATER_LEVEL_DIR", "DUCK_GAUGE_FILE", "WIS_FILE",
                 "STORM_RECORD_DIR", "WAVE_CLIMATE_DIR", "RSLR_ROOT",
                 "RSLR_RECORD_FILE", "RSLR_RATES_CSV", "HINDCAST_STORMS",
                 "STORM_VALIDATION", "STORM_FIGURES", "SPLICED_1984_2024",
                 "SUPERSEDED_STORMS"):
        path = globals()[name]
        print(f"{'ok' if path.exists() else 'MISSING':8} {name:18} "
              f"{path.relative_to(INIT_ROOT).as_posix()}")
    for w in ((1984, 2004), (1996, 2010), (2004, 2024), (2010, 2024)):
        p = storm_series_file(*w)
        print(f"{'ok' if p.exists() else 'MISSING':8} {'storms ' + window_tag(*w):18} "
              f"{p.relative_to(INIT_ROOT).as_posix()}")
