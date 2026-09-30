"""
Where does the source/sink (background erosion) calibration live?

    from site_layer.hat_source_sink import CALIBRATE_ROOT, PASS0_BACKUP

The fit's products, figures and export under data/hatteras_init/7-source-sink/,
one folder per calibration pair; the field the model runs on is in hatteras_site_config. Details: scripts/site_layer/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-22
"""

from __future__ import annotations

from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents
                    if (_p / "pyproject.toml").exists())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

BE_ROOT = INIT_ROOT / "7-source-sink"
CALIBRATE_ROOT = BE_ROOT / "2-calibrate"
FIGURES_ROOT = BE_ROOT / "3-figures"
EXPORT_DIR = BE_ROOT / "4-export"
ARCHIVE = BE_ROOT / "archive"

# Config backups written before each apply pass, shared by every pair
PREBE_DIR = CALIBRATE_ROOT / "prebe"

# The pass-0 field, which plot_be_zones.py and export_be_calibration.py must both read
PASS0_BACKUP = PREBE_DIR / "hatteras_site_config_prebe_20260914_180700.py"

DEFAULT_PERIOD_STARTS = (1984, 2004)
DEFAULT_PAIR_TAG = "1984_2004__2004_2024"


def pair_tag(p1_start, p1_end, p2_start, p2_end):
    """The folder name for a jointly fitted pair of periods."""
    return f"{p1_start}_{p1_end}__{p2_start}_{p2_end}"


def calibrate_dir(tag=DEFAULT_PAIR_TAG):
    """Where one pair's fit writes its tables and convergence history."""
    return CALIBRATE_ROOT / tag


def figures_dir(tag=DEFAULT_PAIR_TAG):
    """Where one pair's figures go (1-field/, 2-method/, 3-limits/ inside)."""
    return FIGURES_ROOT / tag


if __name__ == "__main__":
    for name in ("BE_ROOT", "CALIBRATE_ROOT", "FIGURES_ROOT", "EXPORT_DIR",
                 "ARCHIVE", "PREBE_DIR", "PASS0_BACKUP"):
        path = globals()[name]
        print(f"{'ok' if path.exists() else 'MISSING':8} {name:15} "
              f"{path.relative_to(PROJECT_ROOT).as_posix()}")
    pairs = sorted(p.name for p in CALIBRATE_ROOT.glob("*__*") if p.is_dir())
    print(f"{'':8} {'pairs':15} {', '.join(pairs) or 'none'}")
