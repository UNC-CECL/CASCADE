#!/usr/bin/env python3
"""Run-selecting constants for the Hatteras hindcast, in one place.

WHY THIS MODULE EXISTS
    `HAT_hindcast_1984_2024.ipynb` and its headless mirror
    `HAT_hindcast_1984_2024.py` used to carry these seven values as literals
    in section 3 (plus 7 and 11). Running the scenario matrix therefore meant
    hand-editing a source file between every run, which is how the two
    published edgeBE runs in the retired `run_index.csv` ended up disagreeing
    with each other on the background-erosion values they were fit under.

    Both files now read the values from here, so a driver can select a run
    through the environment without editing tracked source, and an interactive
    notebook session can still override any of them by typing over the
    assignment in section 3.

PRECEDENCE
    environment variable  >  the default in this file

    Nothing else. There is no config file, no CLI parsing, and no per-period
    table: a run is fully described by the seven values below, and `describe()`
    reports which of them came from the environment so a run log states how it
    was driven rather than leaving it to be inferred.

THE SYNC RULE STILL APPLIES
    This module is imported by BOTH the notebook and the .py. Adding a value
    here is only half the change -- section 3 of both files has to read it, or
    the two drift again. See the module docstring of the .py.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import os
from typing import Dict, List

__all__ = ["RUN_CONFIG", "RunConfig", "describe", "ENV_PREFIX"]

ENV_PREFIX = "HAT_"

# Records which names were taken from the environment, for describe(). Filled
# by the readers below as a side effect, which is the only way to distinguish
# "the default happens to equal the environment value" from "not set at all".
_FROM_ENV: List[str] = []


def _read(name: str, default, cast):
    """Reads ENV_PREFIX + name from the environment, or returns the default.

    Args:
        name: Bare constant name, e.g. "START_YEAR". The environment variable
            read is ENV_PREFIX + name.
        default: Value used when the variable is unset or empty.
        cast: Callable converting the raw string to the wanted type. It is
            expected to raise ValueError on malformed input.

    Returns:
        The cast environment value, or `default`.

    Raises:
        ValueError: If the variable is set but cannot be cast. Deliberately
            fatal rather than falling back to the default: a driver that
            misspells a value would otherwise run the default configuration
            under the name of the one it asked for.
    """
    raw = os.environ.get(ENV_PREFIX + name, "")
    if raw == "":
        return default
    try:
        value = cast(raw)
    except ValueError as exc:
        raise ValueError(
            f"{ENV_PREFIX + name}={raw!r} could not be read as "
            f"{cast.__name__}: {exc}") from exc
    _FROM_ENV.append(name)
    return value


def _as_bool(raw: str) -> bool:
    """Casts a boolean environment string, strictly.

    Accepts only unambiguous spellings. "0"/"1" and "true"/"false" are the
    two conventions a shell driver is likely to emit; anything else raises
    rather than being silently truthy, because `bool("False")` is True and
    that failure mode is invisible in a run log.
    """
    lowered = raw.strip().lower()
    if lowered in ("1", "true", "yes", "on"):
        return True
    if lowered in ("0", "false", "no", "off"):
        return False
    raise ValueError(f"expected a boolean spelling, got {raw!r}")


class RunConfig:
    """The seven values that select which run the hindcast performs.

    Attributes:
        start_year: 1984 or 2004; selects a period from HATTERAS_PERIODS.
        source_sink_preset: "zeroBE", "edgeBE" or "calibBE".
        scenario: A key of the SCENARIOS table in section 3.
        groin_enabled: Whether the groin callback is attached.
        groin_trapping_rate_m_yr: M, the groin amplitude knob.
        groin_deterioration_fraction: f, the post-deterioration floor as a
            fraction of M.
        overwrite: True empties an existing run directory and reuses it.
    """

    def __init__(self) -> None:
        self.start_year: int = _read("START_YEAR", 1984, int)
        self.source_sink_preset: str = _read(
            "SOURCE_SINK_PRESET", "zeroBE", str)
        self.scenario: str = _read("SCENARIO", "full_management", str)
        self.groin_enabled: bool = _read("GROIN_ENABLED", False, _as_bool)
        self.groin_trapping_rate_m_yr: float = _read(
            "GROIN_TRAPPING_RATE_M_YR", 50.0, float)
        self.groin_deterioration_fraction: float = _read(
            "GROIN_DETERIORATION_FRACTION", 0.9, float)
        self.overwrite: bool = _read("OVERWRITE", False, _as_bool)
        # The .npz model pickle is ~160 MB per run and is the only large
        # artifact a run produces. It is what lets a figure be re-derived
        # without re-running, so it defaults to on; a long unattended matrix
        # that only needs the rate curves can turn it off and save ~6 GB.
        self.save_model_state: bool = _read(
            "SAVE_MODEL_STATE", True, _as_bool)

    def as_dict(self) -> Dict[str, object]:
        """Returns the configuration as a plain dict, for run metadata."""
        return dict(vars(self))


RUN_CONFIG = RunConfig()


def describe() -> str:
    """Renders the configuration and its provenance as a printable block.

    Section 3 of both the notebook and the .py prints this, so every run log
    records not just what was run but whether the values were driven from the
    environment or are this file's defaults -- the distinction that matters
    when a matrix run and a hand-iterated run land in the same index.

    Returns:
        A multi-line string, no trailing newline.
    """
    lines = ["run configuration"]
    for name, value in RUN_CONFIG.as_dict().items():
        env_name = name.upper()
        origin = ("env " + ENV_PREFIX + env_name if env_name in _FROM_ENV
                  else "default")
        lines.append(f"  {name:<32} {value!r:<20} ({origin})")
    if not _FROM_ENV:
        lines.append("  -- all defaults; nothing was set in the environment")
    return "\n".join(lines)


if __name__ == "__main__":
    print(describe())
