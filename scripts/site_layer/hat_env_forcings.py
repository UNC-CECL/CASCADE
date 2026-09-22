# ==============================================================================
# hat_env_forcings.py
#
# WHERE ARE THE SEA LEVEL AND STORM FORCINGS, AND THE RECORDS THEY COME FROM?
#
# WHY THIS EXISTS
#   A dozen scripts typed data/hatteras_init/3-env-forcings/... themselves,
#   and hatteras_site_config spelled every period's storm file out in full.
#   Same answer as hat_observed_rates.py for 5-scr: resolved ONCE (2026-09-18).
#
# THE LAYOUT, grouped by job (2026-09-18). Before that the records sat beside
# the forcings built from them, the WIS export inside storms/, and the storm
# validation in two places (storms/storm_check/ -- retired 2026-09-22 -- and a validation/ folder
# inside a model-input window).
#
#     data/hatteras_init/3-env-forcings/
#         1-records/            what the forcings are built from, as downloaded
#             water_level/      Duck gauge 8651370, hourly, m NAVD88 (+ cache)
#             WIS_raw_data/     WIS station 63228 wave export (git-ignored)
#             storm_record/     the hurricane history this is checked against
#             wave_climate_duke/  e_phi_0_OBX_yearly.nc; no script reads it
#         2-rslr/               sea level: record/ fits/ figures/ (09-15 layout)
#         3-storms/
#             hindcast_storms/<window>/   the MODEL INPUTS, one per window,
#                                         plus 1984_2024/ (spliced, not a period)
#             validation/<window>/        both validators' output
#             figures/
#         archive/              retired storm series, and Roya's figure
#
# The RSLR record sits in 2-rslr/record/, not 1-records/: rslr/ was laid out as
# record -> fits -> figures on 2026-09-15 and is kept whole.
#
# A WINDOW IS <start>_<end> (the end is a boundary; the model spends
# start..end-1), the same naming rule as the rest of the init tree.
# ==============================================================================

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
# benton_storms/ and testing_storms/, retired 2026-09-14; the storm_check
# validators still read testing_storms/base_storms/.
SUPERSEDED_STORMS = ARCHIVE / "storms_superseded_20260914"


def window_tag(start_year: int, end_year: int) -> str:
    return f"{int(start_year)}_{int(end_year)}"


def storm_window_dir(start_year: int, end_year: int) -> Path:
    """One window's model-facing storm series."""
    return HINDCAST_STORMS / window_tag(start_year, end_year)


def storm_series_file(start_year: int, end_year: int,
                      variant: str = "v3_72") -> Path:
    """The .npy CASCADE reads for a window: <window>_storms_<variant>.npy.
    v3_72 is the v3 method with a 72 h maximum storm duration."""
    tag = window_tag(start_year, end_year)
    return storm_window_dir(start_year, end_year) / f"{tag}_storms_{variant}.npy"


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
