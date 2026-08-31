#!/usr/bin/env python3
"""Run-selecting settings for the Hatteras hindcast, in one place.

WHERE A SETTING IS TYPED
    `hat_run.yaml`, beside this file. Edit it, save it, run the hindcast.
    That file is the interface; this module is the machinery that reads it,
    and its own values are only the fallbacks.

WHY THIS MODULE EXISTS
    `HAT_hindcast_1984_2024.ipynb` and its headless mirror
    `HAT_hindcast_1984_2024.py` used to carry these values as literals in
    sections 1, 3, 7, 9 and 11. Running the scenario matrix therefore meant
    hand-editing a source file between every run, which is how the two
    published edgeBE runs in the retired `run_index.csv` ended up disagreeing
    with each other on the background-erosion values they were fit under.

    Both files now read the values from here, so a run can be selected by
    editing one untracked-by-the-model settings file, or driven through the
    environment without editing anything at all.

PRECEDENCE
    environment variable  >  hat_run.yaml  >  the default in this file

    `HAT_IGNORE_SETTINGS=1` drops the middle term. `HAT_run_all.py` sets it on
    every run it launches, so a batch run is described entirely by the driver
    plus this file's defaults, and a half-finished experiment left in
    hat_run.yaml can never reach the comparison matrix or the sweep.

    `describe()` reports which of the three each value came from, so a run log
    states how it was driven rather than leaving it to be inferred.

RE-READING, FOR THE NOTEBOOK
    `load_run_config()` re-reads hat_run.yaml on every call. Section 3 of both
    files calls it rather than using the module-level RUN_CONFIG, so editing
    the yaml and re-running the cell picks the change up without restarting
    the kernel -- a module-level constant would be cached by the import system
    and the edit would silently not apply.

WHAT IS *NOT* VALIDATED HERE
    Value legality: scenario names, periods and offset modes each have exactly
    one home in the code (`SCENARIOS` in section 3, `HATTERAS_PERIODS`,
    `ISLAND_OFFSET_MODES`), and a second copy here is a second thing to drift.
    A bad value raises in section 3 with the valid list attached. What IS
    validated here is what this module owns: that every key in the yaml is a
    key it knows, and that each value casts to the right type.

THE SYNC RULE STILL APPLIES
    This module is imported by BOTH the notebook and the .py. Adding a field
    here is only part of the change -- the yaml needs a documented key and the
    matching section of both files has to read it, or the two drift again.
    See the module docstring of the .py.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

__all__ = [
    "RUN_CONFIG", "RunConfig", "load_run_config", "describe", "preflight",
    "ENV_PREFIX", "IGNORE_ENV", "SETTINGS_PATH",
]

ENV_PREFIX = "HAT_"
IGNORE_ENV = "HAT_IGNORE_SETTINGS"

SETTINGS_PATH = Path(__file__).resolve().parent / "hat_run.yaml"

# For the runtime estimate in preflight(). scripts/hatteras_ms -> repo root.
_PROJECT_BASE_DIR = Path(__file__).resolve().parents[2]
RUN_INDEX_PATH = _PROJECT_BASE_DIR / "output" / "raw_runs" / "run_index.csv"

# The model .npz, for the preflight cost line. Measured, not guessed: the
# saved states dominate and the file is ~160 MB for either period.
_MODEL_STATE_MB = 160


# =============================================================================
# CASTS
# =============================================================================
# Each raises ValueError on malformed input. Deliberately fatal rather than
# falling back to the default: a driver or a yaml that misspells a value would
# otherwise run the default configuration under the name of the one it asked
# for.

def _as_bool(raw) -> bool:
    """Casts a boolean, strictly.

    Accepts real booleans (which is what the yaml parser produces) and the
    unambiguous string spellings a shell driver emits. Anything else raises
    rather than being silently truthy, because `bool("False")` is True and
    that failure mode is invisible in a run log.
    """
    if isinstance(raw, bool):
        return raw
    lowered = str(raw).strip().lower()
    if lowered in ("1", "true", "yes", "on"):
        return True
    if lowered in ("0", "false", "no", "off"):
        return False
    raise ValueError(f"expected a boolean spelling, got {raw!r}")


def _as_opt_bool(raw) -> Optional[bool]:
    """Casts a boolean that may also be explicitly unset.

    `null` in the yaml and "" or "none" in the environment mean "leave the
    decision to whoever reads this" -- for `relocations` that is the named
    scenario, for `show_figures` it is the notebook/.py split.
    """
    if raw is None:
        return None
    if isinstance(raw, str) and raw.strip().lower() in ("", "none", "null"):
        return None
    return _as_bool(raw)


def _as_int(raw) -> int:
    return int(str(raw).strip())


def _as_float(raw) -> float:
    return float(str(raw).strip())


def _as_str(raw) -> str:
    return str(raw).strip()


# =============================================================================
# THE FIELDS
# =============================================================================
# (attribute, yaml path, cast, default)
#
# The yaml path is a tuple so a nested block reads naturally in the file
# (`groin.enabled`) while staying a flat attribute in code (`groin_enabled`).
# The environment name is always ENV_PREFIX + the attribute upper-cased, and
# is NOT derived from the yaml path -- the ten names that existed before this
# file did are load-bearing in `HAT_run_all.py` and in shell history, and
# renaming them to match a yaml layout would break both.

_FIELDS: Tuple[Tuple[str, Tuple[str, ...], object, object], ...] = (
    ("start_year",                   ("start_year",),        _as_int,      1984),
    ("source_sink_preset",           ("source_sink",),       _as_str,      "zeroBE"),
    ("scenario",                     ("scenario",),          _as_str,      "full_management"),
    ("relocations",                  ("relocations",),       _as_opt_bool, None),
    ("offset_mode",                  ("offset_mode",),       _as_str,      "asrun"),

    ("groin_enabled",                ("groin", "enabled"),         _as_bool,  False),
    # The decided pair, 2026-08-30. M from period 1 (D4-D8 demeaned), f from
    # the 1967-2018 rig -- NOT jointly fitted. These defaults matter only when
    # hat_run.yaml is absent or omits the block; they used to read 50.0 / 0.9,
    # which were placeholders and f = 0.9 was never fitted at all.
    ("groin_trapping_rate_m_yr",     ("groin", "trapping_M"),      _as_float, 60.0),
    ("groin_deterioration_fraction", ("groin", "deterioration_f"), _as_float, 0.6),

    ("hs",                           ("physics", "wave_height_Hs"), _as_float, 2.5),

    # Management, not physics: a defence someone decides to build. Top-level
    # in the yaml for that reason, and not in the `scenario` table because no
    # historical sandbag campaign is reconstructed for either period.
    ("sandbags",                     ("sandbags",),                 _as_bool,  False),

    ("show_figures",                 ("output", "show_figures"),    _as_opt_bool, None),
    ("make_gifs",                    ("output", "make_gifs"),       _as_bool,     True),
    ("save_model_state",             ("output", "save_model_state"), _as_bool,    True),
    ("overwrite",                    ("output", "overwrite"),       _as_bool,     False),

    # Effectively fixed, and deliberately absent from hat_run.yaml -- see the
    # "not settable here" block at the foot of that file. It stays a field so
    # HAT_USE_SANDBOX_CASCADE can still force the installed package for a
    # one-off A/B, and so describe() records which model a run actually built.
    # It must NOT be derived from groin_enabled: section 12.3's paired
    # baseline has to be the same model as the groin run in every respect.
    ("use_sandbox_cascade",          ("use_sandbox_cascade",), _as_bool, True),
)

# Aliases kept so an environment variable that predates the yaml still works.
# HAT_SOURCE_SINK_PRESET is the name HAT_run_all.py sets; the attribute is the
# same, so this is only about the env spelling being longer than the yaml key.
_ENV_ALIASES: Dict[str, Tuple[str, ...]] = {
    "source_sink_preset": ("SOURCE_SINK_PRESET",),
    "hs": ("HS",),
    "sandbags": ("SANDBAGS", "ENABLE_SANDBAG_PLACEMENT"),
}


# =============================================================================
# READING THE SETTINGS FILE
# =============================================================================

def _flatten(mapping, prefix=()) -> Dict[Tuple[str, ...], object]:
    """Flattens a nested yaml mapping to {path tuple: value}.

    A block that is itself a known field's parent (e.g. `groin`) flattens into
    its leaves; a block that is not is reported as an unknown key by the
    caller, path and all, rather than being silently skipped.
    """
    flat: Dict[Tuple[str, ...], object] = {}
    for key, value in (mapping or {}).items():
        path = prefix + (str(key),)
        if isinstance(value, dict):
            flat.update(_flatten(value, path))
        else:
            flat[path] = value
    return flat


def _load_settings_file(path: Path) -> Tuple[Dict[Tuple[str, ...], object], Optional[Path]]:
    """Reads hat_run.yaml, or returns nothing if it is absent or suppressed.

    Returns:
        (flat mapping of yaml path -> value, the path actually read or None).

    Raises:
        RuntimeError: If the file exists but pyyaml is not installed -- the
            settings would otherwise be silently ignored and the run would use
            defaults under the name of whatever the file asked for.
        ValueError: If the file holds a key this module does not know. A
            misspelled key is the failure this catches: ignoring it produces a
            run that used the default while its settings file says otherwise.
    """
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


# =============================================================================
# THE CONFIGURATION
# =============================================================================

class RunConfig:
    """The values that select which run the hindcast performs.

    Attributes:
        start_year: 1984 or 2004; selects a period from HATTERAS_PERIODS.
        source_sink_preset: "zeroBE", "edgeBE" or "calibBE".
        scenario: A key of the SCENARIOS table in section 3.
        relocations: Overrides the scenario's historical-relocation switch.
            None leaves the scenario preset in charge.
        offset_mode: Which shoreline_offset variant to build ("asrun",
            "metres", "detrended").
        groin_enabled: Whether the groin callback is attached.
        groin_trapping_rate_m_yr: M, the groin amplitude knob.
        groin_deterioration_fraction: f, the post-deterioration floor as a
            fraction of M.
        hs: Significant wave height, m.
        sandbags: Whether sandbag placement is enabled.
        show_figures: True renders figures inline. None means the reader
            decides -- the .py uses False, the notebook True.
        make_gifs: Whether section 9's shoreline animations are built.
        save_model_state: Whether the ~160 MB model .npz is written.
        overwrite: True empties an existing run directory and reuses it.
        use_sandbox_cascade: True imports cascade.cascade_groin, which carries
            the hook the groin callback needs.
        origins: attribute -> "default" | "file" | "env", for describe().
        settings_path: The yaml actually read, or None.
    """

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

    def _resolve(self, name, yaml_path, cast, default, file_values):
        """Applies the precedence: environment, then the file, then default."""
        for env_name in (name.upper(),) + _ENV_ALIASES.get(name, ()):
            raw = os.environ.get(ENV_PREFIX + env_name, "")
            if raw != "":
                return self._cast(cast, raw,
                                  f"{ENV_PREFIX + env_name}"), "env"

        if yaml_path in file_values:
            raw = file_values[yaml_path]
            # A yaml `key:` with nothing after it parses to None. For a
            # nullable field that is a deliberate "unset"; for any other it is
            # an unfinished edit, and taking the default silently would run
            # something the file does not say.
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

    def as_dict(self) -> Dict[str, object]:
        """Returns the settings as a plain dict, for run metadata."""
        return {name: getattr(self, name) for name, _, _, _ in _FIELDS}


def load_run_config(settings_path: Optional[Path] = None) -> RunConfig:
    """Reads the settings afresh.

    Call this rather than using the module-level RUN_CONFIG when the file may
    have changed since import -- which in a notebook is every time, since the
    kernel caches the module and an edit to the yaml would otherwise not
    apply until a restart.
    """
    return RunConfig(settings_path)


RUN_CONFIG = load_run_config()


# =============================================================================
# REPORTING
# =============================================================================

def describe(config: Optional[RunConfig] = None) -> str:
    """Renders the settings and their provenance as a printable block.

    Section 3 of both the notebook and the .py prints this, so every run log
    records not just what was run but whether each value was typed in the
    settings file, driven from the environment, or left at this module's
    default -- the distinction that matters when a matrix run and a
    hand-iterated run land in the same index.

    Returns:
        A multi-line string, no trailing newline.
    """
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


def _runtime_estimate(start_year: int, index_path: Optional[Path] = None):
    """Median wall-clock of prior runs of this period, from run_index.csv.

    Measured rather than assumed: the index records `runtime_min` for every
    run that has completed, so the estimate is this machine's own history for
    this period and no constant has to be maintained here.

    Returns:
        (median minutes, number of runs it was taken over), or (None, 0) when
        the index is absent or holds no run of this period.
    """
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


def preflight(run_name: str, run_dir, config: Optional[RunConfig] = None,
              index_path: Optional[Path] = None) -> str:
    """Renders what this run will produce, before it produces it.

    Answers the three questions worth asking before a long run starts: what
    will it be called, where will it land, and does something already live
    there. The name is the section 3 preview, not a name built here -- section
    7.5 derives the authoritative one from what the modules actually built and
    raises if the two disagree, so this can never quietly become the thing
    that names the directory.

    The collision line WARNS rather than raises. `guard_run_dir` in section 11
    is the authority on that, and a second gate here would be a second place
    for the rule to live.

    Args:
        run_name: RUN_NAME_PREVIEW from section 3.
        run_dir: Directory the run will write to.
        config: Settings to report. Defaults to the module-level RUN_CONFIG.
        index_path: run_index.csv, for the runtime estimate. Defaults to the
            repo's own.

    Returns:
        A multi-line string, no trailing newline.
    """
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
