"""
Run names: built from a run's switches, and read back in either spelling.

    from cascade_pipeline.run_names import build_run_name, parse_run_name, other_spellings

Since 2026-10-29 a name is one cell of the experimental design, one field per slot:

    HAT_<start>-<end>_<scenario>_<groin>_<preset>[_<perturbation>...]
    HAT_1996-2009_managed_groin_edgeBE
    HAT_2009-2025_managed-nofill_nogroin_edgeBE_waveHs1p2

Runs made before then keep their old names on disk
(HAT_1996_2009_edgeBE_offsetmetres_road_bdm_groin); parse_run_name reads both, and
other_spellings is how find_run_dir finds a run under the spelling it was not asked for.
Details: scripts/cascade_pipeline/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-29
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

# --- CONFIG ------------------------------------------------------------------
SITE = "HAT"
# (roadway, beach_dune) -> the scenario word; the manuscript's scenario names
SCENARIO_TOKENS = {
    (False, False): "natural",
    (True, False): "road",
    (False, True): "dune",
    (True, True): "managed",
}
# Departures from the scenario, appended to its word with "-" so the field count stays fixed
NOFILL = "nofill"
RELOC = "reloc"
GROIN, NOGROIN = "groin", "nogroin"
# -----------------------------------------------------------------------------

_SCENARIO_FLAGS = {word: flags for flags, word in SCENARIO_TOKENS.items()}
_NEW_PERIOD = re.compile(r"(\d{4})-(\d{4})")
_YEAR = re.compile(r"\d{4}")


@dataclass
class RunName:
    """One run name, read into its fields.

    fills is True when the name records fills applied, False when it records
    them withheld, and None when the name does not say -- the new spelling never
    says, because whether a period has fills follows from the period.
    """
    start_year: int
    end_year: int
    roadway: bool
    beach_dune: bool
    groin: bool
    preset: str
    relocations: bool = False
    fills_withheld: bool = False
    fills: Optional[bool] = None
    perturbations: tuple = field(default_factory=tuple)
    legacy: bool = False

    @property
    def scenario(self):
        """The scenario word with its departures, e.g. "managed-nofill"."""
        return scenario_field(self.roadway, self.beach_dune,
                              self.fills_withheld, self.relocations)

    @property
    def name(self):
        """The name in the current spelling."""
        return build_run_name(self.start_year, self.end_year, self.roadway,
                              self.beach_dune, self.groin, self.preset,
                              self.fills_withheld, self.relocations,
                              self.perturbations)


def scenario_field(roadway, beach_dune, fills_withheld=False, relocations=False):
    """The scenario slot: the scenario word, then any departure from it.

    Args:
        roadway: Roadway management on.
        beach_dune: Beach/dune management on.
        fills_withheld: The period has fills and a manager to place them, but
            they were withheld (the full_no_fill scenario).
        relocations: The historical 1989/1999 road relocations fire.

    Returns:
        e.g. "managed", "managed-nofill", "road-reloc".
    """
    parts = [SCENARIO_TOKENS[(bool(roadway), bool(beach_dune))]]
    if fills_withheld:
        parts.append(NOFILL)
    if relocations:
        parts.append(RELOC)
    return "-".join(parts)


def build_run_name(start_year, end_year, roadway, beach_dune, groin, preset,
                   fills_withheld=False, relocations=False, perturbations=()):
    """The run name, which is also the run's directory name.

    Args:
        start_year, end_year: The hindcast window.
        roadway, beach_dune, fills_withheld, relocations: See scenario_field.
        groin: The groin is on.
        preset: The canonical source/sink preset, e.g. "edgeBE".
        perturbations: Sensitivity tokens ("rset40", "waveHs1p2"); None and ""
            are skipped. The wave token goes last: sweep_family reads the axis
            off the last field.

    Returns:
        e.g. "HAT_1996-2009_managed_groin_edgeBE".
    """
    fields = [SITE, f"{int(start_year)}-{int(end_year)}",
              scenario_field(roadway, beach_dune, fills_withheld, relocations),
              GROIN if groin else NOGROIN, preset]
    fields += [p for p in perturbations if p]
    return "_".join(fields)


def _parse_new(tokens):
    period = _NEW_PERIOD.fullmatch(tokens[1]) if len(tokens) >= 5 else None
    if not period:
        return None
    scenario, *departures = tokens[2].split("-")
    if (scenario not in _SCENARIO_FLAGS or tokens[3] not in (GROIN, NOGROIN)
            or any(d not in (NOFILL, RELOC) for d in departures)):
        return None
    roadway, beach_dune = _SCENARIO_FLAGS[scenario]
    return RunName(int(period[1]), int(period[2]), roadway, beach_dune,
                   tokens[3] == GROIN, tokens[4], relocations=RELOC in departures,
                   fills_withheld=NOFILL in departures,
                   fills=False if NOFILL in departures else None,
                   perturbations=tuple(tokens[5:]))


# HAT_<s>_<e>_<preset>_offsetmetres_(no)road_[reloc_](no)bdm_[(no)nourish_](no)groin[_rset..][_wave..]
def _parse_legacy(tokens):
    if (len(tokens) < 8 or not _YEAR.fullmatch(tokens[1])
            or not _YEAR.fullmatch(tokens[2]) or tokens[4] != "offsetmetres"):
        return None
    rest = tokens[5:]

    def take(*choices):
        if rest and rest[0] in choices:
            return rest.pop(0)
        return None

    road = take("road", "noroad")
    reloc = take(RELOC)
    bdm = take("bdm", "nobdm")
    nourish = take("nourish", "nonourish")
    groin = take(GROIN, NOGROIN)
    if road is None or bdm is None or groin is None:
        return None
    return RunName(int(tokens[1]), int(tokens[2]), road == "road", bdm == "bdm",
                   groin == GROIN, tokens[3], relocations=bool(reloc),
                   fills_withheld=nourish == "nonourish",
                   fills={"nourish": True, "nonourish": False}.get(nourish),
                   perturbations=tuple(rest), legacy=True)


def parse_run_name(name):
    """Reads a run name, in either spelling, into its fields.

    The old spelling is read only with its `offsetmetres` token: a name without
    it is a run from before the offsets were in metres, a different model, and
    is not translated.

    Args:
        name: A run name or directory name.

    Returns:
        A RunName, or None if the name is neither spelling (experiment
        controls, simple_* runs, pre-metres runs).
    """
    tokens = str(name).split("_")
    if not tokens or tokens[0] != SITE:
        return None
    return _parse_new(tokens) or _parse_legacy(tokens)


def legacy_run_names(parsed):
    """The old spellings a run could have been made under, most likely first.

    The old spelling wrote `nourish` only when fills were applied in the
    period, which the new spelling leaves to the period; both are offered.
    """
    stem = [SITE, str(parsed.start_year), str(parsed.end_year), parsed.preset,
            "offsetmetres", "road" if parsed.roadway else "noroad"]
    if parsed.relocations:
        stem.append(RELOC)
    stem.append("bdm" if parsed.beach_dune else "nobdm")
    tail = [GROIN if parsed.groin else NOGROIN, *parsed.perturbations]
    if parsed.fills_withheld:
        fills = [["nonourish"]]
    elif parsed.beach_dune:
        fills = [["nourish"], []] if parsed.fills is not False else [[]]
    else:
        fills = [[]]
    return ["_".join(stem + f + tail) for f in fills]


def other_spellings(name):
    """Every other name the same run could be filed under.

    Args:
        name: A run name in either spelling.

    Returns:
        A list of names, empty when the name cannot be read.
    """
    parsed = parse_run_name(name)
    if parsed is None:
        return []
    candidates = [parsed.name] + legacy_run_names(parsed)
    return [c for c in dict.fromkeys(candidates) if c != name]
