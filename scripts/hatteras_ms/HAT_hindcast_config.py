#!/usr/bin/env python3
"""
Run-selecting settings for the Hatteras hindcast, in one place: hat_run.yaml, overridable by HAT_* variables.

    python scripts/hatteras_ms/HAT_hindcast_config.py   # prints the resolved settings

Imported by the runner, the notebook and every driver, so a run is chosen
without editing tracked source. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

__all__ = [
    "RUN_CONFIG", "RunConfig", "load_run_config", "describe", "preflight",
    "field_default", "ENV_PREFIX", "IGNORE_ENV", "SETTINGS_PATH",
]

# --- CONFIG ------------------------------------------------------------------
ENV_PREFIX = "HAT_"
IGNORE_ENV = "HAT_IGNORE_SETTINGS"

SETTINGS_PATH = Path(__file__).resolve().parent / "hat_run.yaml"

# For the runtime estimate in preflight(). scripts/hatteras_ms -> repo root.
_PROJECT_BASE_DIR = next(_p for _p in Path(__file__).resolve().parents
                         if (_p / "pyproject.toml").exists())
RUN_INDEX_PATH = _PROJECT_BASE_DIR / "output" / "raw_runs" / "run_index.csv"

# The model .npz, for the preflight cost line
_MODEL_STATE_MB = 160
# -----------------------------------------------------------------------------


# Casts

# Each raises ValueError on malformed input

# Casts a boolean, strictly
def _as_bool(raw) -> bool:
    if isinstance(raw, bool):
        return raw
    lowered = str(raw).strip().lower()
    if lowered in ("1", "true", "yes", "on"):
        return True
    if lowered in ("0", "false", "no", "off"):
        return False
    raise ValueError(f"expected a boolean spelling, got {raw!r}")


# Casts a boolean that may also be explicitly unset
def _as_opt_bool(raw) -> Optional[bool]:
    if raw is None:
        return None
    if isinstance(raw, str) and raw.strip().lower() in ("", "none", "null"):
        return None
    return _as_bool(raw)


# Casts a float that may also be explicitly unset
def _as_opt_float(raw) -> Optional[float]:
    if raw is None:
        return None
    if isinstance(raw, str) and raw.strip().lower() in (
            "", "none", "null", "measured"):
        return None
    return _as_float(raw)


# An int from a yaml or environment value
def _as_int(raw) -> int:
    return int(str(raw).strip())


# A float from a yaml or environment value
def _as_float(raw) -> float:
    return float(str(raw).strip())


# A stripped string from a yaml or environment value
def _as_str(raw) -> str:
    return str(raw).strip()


# The fields

# (attribute, yaml path, cast, default); env name = HAT_ + attribute, never from the yaml path

_FIELDS: Tuple[Tuple[str, Tuple[str, ...], object, object], ...] = (
    # The project ran first; any other start is reachable by naming it
    ("start_year",                   ("start_year",),        _as_int,      1996),
    ("source_sink_preset",           ("source_sink",),       _as_str,      "zeroBE"),
    ("scenario",                     ("scenario",),          _as_str,      "full_management"),
    ("relocations",                  ("relocations",),       _as_opt_bool, None),
    # METRES SINCE 2026-09-24 (Hannah)
    ("offset_mode",                  ("offset_mode",),       _as_str,      "metres"),
    # The reach (2026-09-16): a name from hat_extension_domains.GEOMETRIES
    ("geometry",                     ("geometry",),          _as_str,      "base"),

    # WHERE A RUN IS FILED (2026-09-16)
    ("run_kind",                     ("run_kind",),          _as_str,      "matrix"),
    ("run_tag",                      ("run_tag",),           _as_str,      ""),

    ("groin_enabled",                ("groin", "enabled"),         _as_bool,  False),
    # M: the dipole's amplitude, unused by the pinned blocking groin (2026-08-30 value, stale under option A)
    ("groin_trapping_rate_m_yr",     ("groin", "trapping_M"),      _as_float, 60.0),
    ("groin_deterioration_fraction", ("groin", "deterioration_f"), _as_float, 0.6),
    # PINNED 2026-10-05: blocking b 0.6, f 0.6, fitted on the 1996-2009 calibration period
    # (hard-structures/groin/groin-module-test/1-dem-to-dem/2026-10-05-blocking-fit-calibration)
    ("groin_kind",                   ("groin", "kind"),            _as_str,   "blocking"),
    ("groin_blocking_fraction",      ("groin", "blocking_b"),      _as_float, 0.6),

    # Wave climate: option A since 2026-09-27; moving one earns a name token (README)
    ("hs",                           ("physics", "wave_height_Hs"), _as_float, 2.0),
    ("wave_period_s",                ("physics", "wave_period_s"),  _as_float, 7.5),
    ("wave_asymmetry",               ("physics", "wave_asymmetry"), _as_float, 0.6),
    ("wave_angle_high_fraction",     ("physics", "wave_angle_high_fraction"),
                                                                    _as_float, 0.5),

    # WHERE A RELOCATED ROAD GOES, in metres behind the dune line
    ("relocation_setback_m",         ("relocation_setback_m",),   _as_opt_float, 20.0),

    # Management, not physics: a defence someone decides to build
    ("sandbags",                     ("sandbags",),                 _as_bool,  False),

    ("show_figures",                 ("output", "show_figures"),    _as_opt_bool, None),
    ("make_gifs",                    ("output", "make_gifs"),       _as_bool,     True),
    ("save_model_state",             ("output", "save_model_state"), _as_bool,    True),
    ("overwrite",                    ("output", "overwrite"),       _as_bool,     False),

    # Effectively fixed, and deliberately absent from hat_run.yaml
    ("use_sandbox_cascade",          ("use_sandbox_cascade",), _as_bool, True),
)

# Aliases kept so an environment variable that predates the yaml still works
_ENV_ALIASES: Dict[str, Tuple[str, ...]] = {
    "source_sink_preset": ("SOURCE_SINK_PRESET",),
    "hs": ("HS",),
    "sandbags": ("SANDBAGS", "ENABLE_SANDBAG_PLACEMENT"),
}


# The code default for one field, ignoring the yaml and the environment
def field_default(name: str):
    for field, _, _, default in _FIELDS:
        if field == name:
            return default
    raise KeyError(f"no such setting: {name!r}; "
                   f"have {sorted(f for f, _, _, _ in _FIELDS)}")


# Reading the settings file

# Flattens a nested yaml mapping to {path tuple
def _flatten(mapping, prefix=()) -> Dict[Tuple[str, ...], object]:
    flat: Dict[Tuple[str, ...], object] = {}
    for key, value in (mapping or {}).items():
        path = prefix + (str(key),)
        if isinstance(value, dict):
            flat.update(_flatten(value, path))
        else:
            flat[path] = value
    return flat


# Reads hat_run.yaml, or returns nothing if it is absent or suppressed
def _load_settings_file(path: Path) -> Tuple[Dict[Tuple[str, ...], object], Optional[Path]]:
    if _as_bool(os.environ.get(IGNORE_ENV, "0")):
        return {}, None
    if not path.exists():
        return {}, None

    try:
        import yaml
    except ImportError as exc:                          # pragma: no cover
        raise RuntimeError(
            f"{path.name} exists but pyyaml is not installed, so its settings "
            f"would be silently ignored. `pip install pyyaml`, or delete the "
            f"file to run on defaults.") from exc

    with open(path, "r", encoding="utf-8") as handle:
        raw = yaml.safe_load(handle)
    if raw is None:                                     # an empty file
        return {}, path
    if not isinstance(raw, dict):
        raise ValueError(f"{path} must hold a mapping, got {type(raw).__name__}")

    flat = _flatten(raw)
    known = {yaml_path for _, yaml_path, _, _ in _FIELDS}
    unknown = sorted(".".join(p) for p in flat if p not in known)
    if unknown:
        raise ValueError(
            f"{path.name} has {len(unknown)} key(s) this runner does not "
            f"know: {', '.join(unknown)}\n"
            f"  known keys: {', '.join(sorted('.'.join(p) for p in known))}\n"
            f"A misspelled key is not ignored here on purpose -- it would "
            f"leave the run using the default while the file claims "
            f"otherwise.")
    return flat, path


# The configuration

# The values that select which run the hindcast performs
class RunConfig:

    def __init__(self, settings_path: Optional[Path] = None) -> None:
        path = SETTINGS_PATH if settings_path is None else Path(settings_path)
        file_values, read_from = _load_settings_file(path)
        self.settings_path: Optional[Path] = read_from
        self.origins: Dict[str, str] = {}

        for name, yaml_path, cast, default in _FIELDS:
            value, origin = self._resolve(
                name, yaml_path, cast, default, file_values)
            setattr(self, name, value)
            self.origins[name] = origin

    # Applies the precedence: environment, then the file, then default
    def _resolve(self, name, yaml_path, cast, default, file_values):
        for env_name in (name.upper(),) + _ENV_ALIASES.get(name, ()):
            raw = os.environ.get(ENV_PREFIX + env_name, "")
            if raw != "":
                return self._cast(cast, raw,
                                  f"{ENV_PREFIX + env_name}"), "env"

        if yaml_path in file_values:
            raw = file_values[yaml_path]
            # A yaml `key:` with nothing after it parses to None
            if raw is None and cast is not _as_opt_bool:
                raise ValueError(
                    f"{'.'.join(yaml_path)} in {self.settings_path} is empty. "
                    f"Give it a value, or delete the line to use the default "
                    f"({default!r}).")
            return self._cast(cast, raw,
                              f"{'.'.join(yaml_path)} in "
                              f"{getattr(self.settings_path, 'name', '?')}"), "file"

        return default, "default"

    @staticmethod
    def _cast(cast, raw, source):
        try:
            return cast(raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{source} = {raw!r} could not be read as "
                f"{cast.__name__.lstrip('_').replace('as_', '')}: {exc}"
            ) from exc

    # Returns the settings as a plain dict, for run metadata
    def as_dict(self) -> Dict[str, object]:
        return {name: getattr(self, name) for name, _, _, _ in _FIELDS}


# Reads the settings afresh
def load_run_config(settings_path: Optional[Path] = None) -> RunConfig:
    return RunConfig(settings_path)


RUN_CONFIG = load_run_config()


# Reporting

# Renders the settings and their provenance as a printable block
def describe(config: Optional[RunConfig] = None) -> str:
    config = RUN_CONFIG if config is None else config

    if config.settings_path is not None:
        header = f"run settings   ({config.settings_path.name})"
    elif _as_bool(os.environ.get(IGNORE_ENV, "0")):
        header = f"run settings   ({IGNORE_ENV}=1 -- the file was not read)"
    else:
        header = f"run settings   (no {SETTINGS_PATH.name}; env and defaults)"

    lines = [header]
    for name, yaml_path, _, _ in _FIELDS:
        origin = config.origins[name]
        if origin == "env":
            where = "env " + ENV_PREFIX + name.upper()
        elif origin == "file":
            where = "file " + ".".join(yaml_path)
        else:
            where = "default"
        lines.append(f"  {name:<30} {getattr(config, name)!r:<18} ({where})")

    counts: List[str] = []
    for label in ("file", "env", "default"):
        n = sum(1 for o in config.origins.values() if o == label)
        if n:
            counts.append(f"{n} {label}")
    lines.append("  " + ", ".join(counts))
    return "\n".join(lines)


# Median wall-clock of prior runs of this period, from run_index.csv
def _runtime_estimate(start_year: int, index_path: Optional[Path] = None):
    path = RUN_INDEX_PATH if index_path is None else Path(index_path)
    if not path.exists():
        return None, 0
    try:
        import pandas as pd
        frame = pd.read_csv(path)
        if not {"start_year", "runtime_min"} <= set(frame.columns):
            return None, 0
        same = frame.loc[frame["start_year"] == start_year, "runtime_min"]
        same = same.dropna()
        if same.empty:
            return None, 0
        return float(same.median()), int(same.size)
    except Exception:
        # An unreadable index must not stop a run: this is an advisory line.
        return None, 0


# Renders what this run will produce, before it produces it
def preflight(run_name: str, run_dir, config: Optional[RunConfig] = None,
              index_path: Optional[Path] = None) -> str:
    config = RUN_CONFIG if config is None else config
    run_dir = Path(run_dir)

    lines = ["preflight",
             f"  run name              {run_name}",
             f"  output directory      {run_dir}"]

    existing = sorted(p for p in run_dir.glob("*") if p.is_file()) \
        if run_dir.exists() else []
    if not existing:
        lines.append("  directory             new")
    elif config.overwrite:
        lines.append(f"  directory             EXISTS, {len(existing)} file(s) "
                     f"-- overwrite: true EMPTIES it in section 11")
    else:
        lines.append(f"  directory             EXISTS, {len(existing)} file(s) "
                     f"-- overwrite: false STOPS the run in section 11")

    minutes, n = _runtime_estimate(config.start_year, index_path)
    if minutes is None:
        lines.append("  estimated runtime     unknown (no prior run of this "
                     "period in run_index.csv)")
    else:
        lines.append(f"  estimated runtime     ~{minutes:.1f} min "
                     f"(median of {n} prior {config.start_year} run"
                     f"{'s' if n != 1 else ''})")

    disk = [f"~{_MODEL_STATE_MB} MB model .npz"] if config.save_model_state \
        else ["no model .npz (save_model_state: false)"]
    if config.make_gifs:
        disk.append("shoreline GIFs")
    lines.append("  writes                " + ", ".join(disk))

    return "\n".join(lines)


if __name__ == "__main__":
    print(describe())
