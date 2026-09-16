"""Run provenance, output-directory guarding, and the cross-run index.

The hindcast is run as a matrix -- period x source/sink preset x groin -- and
the three failure modes that costs are all bookkeeping ones:

  * a re-run silently overwriting the outputs of the run it was meant to be
    compared against (guard_run_dir),
  * a finished run whose metadata cannot say which code produced it, because
    the sandbox flag or the extractor version was flipped days earlier
    (git_provenance, and the [identity] section callers pass),
  * twelve runs whose results can only be compared by opening twelve
    hand-formatted text files (rebuild_run_index, which derives run_index.csv
    from every run's metadata; append_run_index is the pre-2026-09-16 writer
    and is kept for the groin sweep).

Metadata is written twice from ONE structure: a .txt to read and a .json to
parse. Rendering both from the same `sections` mapping is what keeps them from
disagreeing -- the previous inline version built only the prose, so anything
downstream had to re-parse it.

Used by both HAT_hindcast_1984_2024.ipynb and its headless mirror
HAT_hindcast_1984_2024.py, so the two cannot drift apart.
"""

import datetime
import hashlib
import json
import re
import subprocess
import shutil
from pathlib import Path

from cascade_pipeline.run_layout import SUBFOLDERS

import numpy as np
import pandas as pd

# Files whose presence means a run directory holds real output. A directory
# containing only these is treated as empty, so a stray .gitkeep or an
# editor's .DS_Store does not block a run.
_IGNORABLE_NAMES = {".gitkeep", ".gitignore", ".DS_Store", "Thumbs.db"}

RUN_INDEX_FILENAME = "run_index.csv"


# Paths a run WROTE BACK into the repository, so their being modified says
# nothing about whether the run's code was committed. See git_provenance.
# Since 2026-09-16 the runner copies the parameters yaml into the run
# directory and CASCADE rewrites only the copy, so this file should stay
# clean; it is still excluded so a tree dirtied by an older run reports right.
EXCLUDED_FROM_DIRTY = (
    "data/hatteras_init/Hatteras-CASCADE-parameters.yaml",
)


def git_provenance(repo_root):
    """Records which commit produced a run, and whether the tree was clean.

    A dirty tree is not an error -- most runs happen mid-edit -- but it does
    mean the commit hash alone will not reproduce the run, so the flag is
    recorded beside it rather than inferred later.

    ONE PATH IS EXCLUDED, AND WITHOUT IT THE FLAG IS USELESS.
    data/hatteras_init/Hatteras-CASCADE-parameters.yaml is TRACKED and is
    REWRITTEN BY EVERY CASCADE CONSTRUCTION -- it is the shared file behind
    the "never run two sweep orchestrators at once" rule. So the moment any
    run starts the tree is dirty and stays dirty, and `dirty` was True on
    every report this pipeline had ever written, including runs whose code
    was fully committed. A flag that cannot be False carries no information.

    Excluding it makes DIRTY TREE mean what it is read as meaning:
    uncommitted CODE or INPUTS, not the run's own scratch output. Found
    2026-08-31, when five relocation comparisons were stamped DIRTY TREE and
    the only run-relevant dirty paths turned out to be this yaml and one
    benign helper addition.

    Anything else volatile that a run writes back into the repository belongs
    in EXCLUDED_FROM_DIRTY too -- otherwise it re-breaks the flag silently.

    Args:
        repo_root: Path to the repository root.

    Returns:
        A dict with commit, branch and dirty. Values are the string "unknown"
        (and dirty None) if git is unavailable or this is not a repository, so
        provenance capture can never be the thing that fails a run.
    """
    def _git(*args):
        return subprocess.run(
            ["git", *args], cwd=str(repo_root), capture_output=True,
            text=True, timeout=15, check=True).stdout.strip()

    try:
        return {
            "commit": _git("rev-parse", "HEAD"),
            "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
            "dirty": bool(_git("status", "--porcelain", "--", ".",
                               *(f":!{path}" for path in EXCLUDED_FROM_DIRTY))),
        }
    except (subprocess.SubprocessError, OSError):
        return {"commit": "unknown", "branch": "unknown", "dirty": None}


def values_digest(mapping, length=12):
    """Fingerprints a {domain: rate} mapping.

    A preset's NAME and its VALUES are separate facts. The source/sink table
    is edited between runs -- re-solving an end domain, testing a different
    edge value -- while the preset keeps the name "edgeBE", so a run that
    records only the name cannot be told apart from the trial before it. Two
    scalar columns cover the end domains; this covers everything else,
    including the interior domains of a 90-domain calibrated fit that no
    column carries.

    Args:
        mapping: Mapping of domain id to rate.
        length: Characters of hex digest to keep. 12 is ~10^-14 collision
            odds over any realistic number of runs.

    Returns:
        A hex string, or "empty" when the mapping is empty -- a preset that
        imposes nothing has nothing to fingerprint, and saying so directly
        beats a hash of the empty string that looks like a real value.
    """
    if not mapping:
        return "empty"
    # Sorted and formatted rather than hashed off repr(): dict order and float
    # repr are not things to make a run's identity depend on.
    payload = ";".join(f"{key}:{float(value):.6g}"
                       for key, value in sorted(mapping.items()))
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()[:length]


# =============================================================================
# WHERE A RUN LIVES
# =============================================================================
# One run's outputs are at
#
#     <raw_runs>/[<arm>/]<start>_<end>/<preset>/<run_name>/
#
# and the ARM COMPONENT IS ABSENT for the calibration arm. That asymmetry is
# deliberate -- every run made before forcing arms existed is at the short
# path, and emitting the component unconditionally would rename all of them --
# but it does mean the path cannot be built by joining a fixed number of parts.
#
# THIS IS THE ONLY PLACE THAT SPELLING BELONGS. It was previously rebuilt by
# hand in six scripts, five of which predate the arm component and so join
# <period>/<preset>/<name> with no slot for it: an arm-scoped run is simply
# invisible to them. HAT_plot_sensitivity.py skipped its target-window check
# silently whenever the path did not resolve, which is the quiet-wrong-path
# failure hat_topo_version.py exists to end for the domain arrays, one tree
# over. A name that is not on disk is an error here, and the error names the
# arms the run IS under.

# =============================================================================
# THE TREE (2026-09-16): filed by PURPOSE
# =============================================================================
# A run is filed by what it is FOR, then by period and preset:
#
#   matrix/<period>/<preset>/<run_name>/                  the production runs
#   sensitivity/<axis>/<period>/<preset>/<run_name>_<token>/   sweep cells
#   experiments/<tag>/<period>/<preset>/<run_name>/       one question each
#   versions/<tag>/<period>/<preset>/<run_name>/          input-version pairs
#   archive/<tag>/<period>/<preset>/<run_name>/           superseded, intact
#
# The run NAME still describes the scenario and is derived from the switches
# the runner built. What used to be an "arm" -- one string that meant a
# forcing value, an input version or an ad hoc experiment label without saying
# which -- is now a KIND (one of KINDS) and a TAG. A matrix run has no tag. A
# sensitivity cell's tag is its axis folder, derived from the trailing token on
# its name, so the value it was swept to is in the NAME and the axis is in the
# PATH. An experiment's or version's tag is the folder it was filed under,
# "<set>/<member>" for a set of related runs.
#
# Why (Hannah, 2026-09-16): a wave sweep had fanned out into twelve top-level
# arms/waveHs<x>/ folders, each holding one run; nineteen of thirty arms were
# finished one-off experiments nothing marked as finished; and the layout
# could not tell those from the version comparisons that are kept on purpose.
#
# BOTH EARLIER LAYOUTS STILL READ. find_run_dir tries the purpose path, then
# the 2026-09-10 layout (arms/<arm>/..., <period>/<preset>/sweeps/<family>/),
# then the flat pre-09-10 one, so a tree that is half migrated resolves. The
# `arm=` keyword is accepted everywhere as a legacy spelling and translated
# through LEGACY_ARMS.

KINDS = ("matrix", "sensitivity", "experiment", "version", "archive")
MATRIX_KIND = "matrix"
KIND_DIR = {
    "matrix": "matrix",
    "sensitivity": "sensitivity",
    "experiment": "experiments",
    "version": "versions",
    "archive": "archive",
}
# The legacy name of the unscoped tree. Still what a pre-09-16 index row says
# in its `arm` column, and what old call sites pass.
CALIBRATION_ARM = "calibration"
ARMS_DIR = "arms"
SWEEPS_DIR = "sweeps"
# Trailing name tokens that mark a sensitivity cell, and the axis folder each
# files under. The runner derives the token (cascade_pipeline.hindcast's
# wave_climate_token / relocation_setback_token); a token combining two wave
# fields ("waveHs3Tp10") files under the FIRST family it starts with.
SWEEP_FAMILIES = ("waveHs", "waveTp", "waveahf", "waveasym", "rset")

# Where each pre-09-16 arm was filed on 2026-09-16, decided by Hannah in the
# same interview: the three version comparisons keep their names under
# versions/; everything else was a one-off experiment and is filed under
# experiments/ by the date it was run, with the arm's own name as the member.
# The twelve waveHs<x> arms are ABSENT on purpose: those were the 1996 wave
# cells filed by the 2026-09-01 rule, re-run as sensitivity cells and deleted.
LEGACY_ARMS = {
    "offsetv1": ("version", "offsetv1"),
    "version-check/v1": ("version", "version-check/v1"),
    "version-pair/v2": ("version", "version-pair/v2"),
    "version-pair/v3": ("version", "version-pair/v3"),
    "pea1989base": ("experiment", "2026-09-02-pea1989/base"),
    "pea1989basenoreloc": ("experiment", "2026-09-02-pea1989/basenoreloc"),
    "behindroad-copy": ("experiment", "2026-09-08-behindroad-copy"),
    "paramsplit-check": ("experiment", "2026-09-14-paramsplit/check"),
    "paramsplit-final": ("experiment", "2026-09-14-paramsplit/final"),
    "paramsplit-oldcode": ("experiment", "2026-09-14-paramsplit/oldcode"),
    "paramsplit-syncon": ("experiment", "2026-09-14-paramsplit/syncon"),
    "probe-natural": ("experiment", "2026-09-14-probe/natural"),
    "probe-paired": ("experiment", "2026-09-14-probe/paired"),
    "probe-roadonly": ("experiment", "2026-09-14-probe/roadonly"),
    "probe-unrounded": ("experiment", "2026-09-14-probe/unrounded"),
    "probe-v1setbacks": ("experiment", "2026-09-14-probe/v1setbacks"),
    "recode-20260914": ("experiment", "2026-09-14-recode"),
    "currency-20260914": ("experiment", "2026-09-14-currency"),
}

# A period directory is exactly <4 digits>_<4 digits>. Anything else directly
# under a kind folder is a tag, so the two levels can be told apart without a
# registry of tag names.
_PERIOD_DIR = re.compile(r"\d{4}_\d{4}")


def period_component(period):
    """The <start>_<end> path component, from either spelling of a period.

    Callers hold a period as a string in some places and as two integers in
    others; accepting both is what lets every call site pass what it already
    has rather than reformatting at each one.

    Args:
        period: Either "1984_2004" or a (start_year, end_year) pair.

    Returns:
        The directory component as a string.

    Raises:
        ValueError: If the string is not <4 digits>_<4 digits>, or the pair is
            not two values. A malformed period would otherwise build a path
            that cannot exist, and be reported as a missing run.
    """
    if isinstance(period, str):
        if not _PERIOD_DIR.fullmatch(period):
            raise ValueError(
                f"period {period!r} is not <start>_<end>, e.g. '1984_2004'")
        return period
    try:
        start, end = period
    except (TypeError, ValueError):
        raise ValueError(
            f"period must be '1984_2004' or (1984, 2004), got {period!r}"
        ) from None
    return f"{int(start)}_{int(end)}"


def sweep_family(run_name):
    """The sensitivity axis a run belongs to, or "" if it is a scenario run.

    Args:
        run_name: The run's derived name, which is also its directory name.

    Returns:
        One of SWEEP_FAMILIES, or "" when the trailing token names none of
        them -- which is every scenario run.
    """
    token = run_name.rpartition("_")[2]
    for family in SWEEP_FAMILIES:
        if token.startswith(family):
            return family
    return ""


def check_tag(tag):
    """Validates a tag: one path component, or two joined by '/'.

    Two levels is a set and its members -- the four probes of one experiment,
    the two versions of one pair -- and a third would be a taxonomy nobody
    has asked for. The value reaches here from an environment variable and is
    joined onto the output root, so it must not escape it.

    Args:
        tag: The tag string. None and "" mean no tag.

    Returns:
        The tag, stripped; "" for none.

    Raises:
        ValueError: If it is not one or two clean path components.
    """
    tag = (tag or "").strip()
    if not tag:
        return ""
    parts = tag.split("/")
    if ("\\" in tag or len(parts) > 2
            or any(not p or p.startswith(".") for p in parts)):
        raise ValueError(
            f"tag {tag!r} must be one path component, or two joined by '/' "
            f"for a set and its member -- it is joined onto the raw_runs "
            f"root and must not escape it.")
    return tag


def check_kind(kind):
    """Validates a kind, defaulting None and "" to matrix."""
    kind = (kind or MATRIX_KIND).strip()
    if kind not in KINDS:
        raise ValueError(f"kind {kind!r} is not one of {KINDS}")
    return kind


def legacy_arm_to_kind_tag(arm):
    """(kind, tag) for a pre-2026-09-16 arm name.

    Args:
        arm: The `arm` column of an old index row, or the value an old call
            site passes as arm=. None, "" and "calibration" are the matrix.

    Returns:
        (kind, tag). A wave arm (`waveHs1p2`) maps to a sensitivity cell
        under the family its name starts with; an arm in LEGACY_ARMS maps as
        that table says; anything else is an experiment tagged with the arm
        name itself, so an unknown arm still resolves somewhere sensible.
    """
    arm = (arm or "").strip()
    if not arm or arm == CALIBRATION_ARM:
        return MATRIX_KIND, ""
    if arm in LEGACY_ARMS:
        return LEGACY_ARMS[arm]
    for family in SWEEP_FAMILIES:
        if arm.startswith(family):
            return "sensitivity", family
    return "experiment", arm


def _resolve_identity(kind, tag, arm, run_name=None):
    """Normalises the three spellings a caller may use into (kind, tag).

    `arm=` wins when given, because a call site still passing it is one that
    has not been updated and means the OLD thing. For a sensitivity cell the
    tag is the axis folder, derived from the name when not given.
    """
    if kind not in (None, "") and kind not in KINDS:
        # A pre-09-16 caller passing the ARM positionally, where kind now
        # sits (find_run_dir(raw, name, period, preset, "version-pair/v2")).
        if arm not in (None, ""):
            raise ValueError(f"{kind!r} is not a kind, and arm= was also given")
        arm, kind = kind, None
    if arm is not None and arm != "":
        if kind not in (None, "", MATRIX_KIND) or tag:
            raise ValueError("pass either arm= (legacy) or kind=/tag=, not both")
        kind, tag = legacy_arm_to_kind_tag(arm)
    kind = check_kind(kind)
    tag = check_tag(tag)
    if kind == "sensitivity":
        family = sweep_family(run_name or "")
        if not tag:
            if not family:
                raise ValueError(
                    f"a sensitivity run needs an axis: {run_name!r} carries no "
                    f"sweep token and no tag was given")
            tag = family
    if kind == MATRIX_KIND and tag:
        raise ValueError(f"a matrix run carries no tag, got {tag!r}")
    return kind, tag


def preset_dir_for(raw_runs, period, preset, kind=MATRIX_KIND, tag="",
                   arm=None):
    """The directory holding every run of one period, preset, kind and tag.

    This is the runner's OUTPUT_BASE_DIR. It is the level anything that
    ENUMERATES runs works at -- the scenario grid, the relocation comparison,
    the source/sink calibration -- as against `run_dir_for`, which is for a run
    already named.

    Args:
        raw_runs: The output/raw_runs root.
        period: "1984_2004" or (1984, 2004).
        preset: Source/sink preset, e.g. "calibBE".
        kind: One of KINDS; the default is the matrix.
        tag: The experiment/version tag, or the axis for a sensitivity run.
        arm: LEGACY. A pre-09-16 arm name, translated through LEGACY_ARMS.

    Returns:
        The preset directory as a Path. Does not check it exists.
    """
    kind, tag = _resolve_identity(kind, tag, arm)
    base = Path(raw_runs) / KIND_DIR[kind]
    if tag:
        base = base / tag
    return base / period_component(period) / preset


def run_dir_for(raw_runs, run_name, period, preset, kind=MATRIX_KIND, tag="",
                arm=None):
    """The directory one run's output belongs in. Does not check it exists.

    The inverse of the runner's RUN_DIR, and the only place the layout is
    spelled. Use `find_run_dir` to READ a finished run; this builds the path a
    run would be WRITTEN to, which is what a writer and a collision guard need
    and what a reader should not be doing by hand.

    Args:
        raw_runs: The output/raw_runs root.
        run_name: The run's derived name, which is also its directory name.
        period: "1984_2004" or (1984, 2004).
        preset: Source/sink preset, e.g. "calibBE".
        kind: One of KINDS; the default is the matrix.
        tag: The experiment/version tag. For a sensitivity run it is derived
            from the name's token when not given.
        arm: LEGACY. A pre-09-16 arm name.

    Returns:
        The run directory as a Path.
    """
    kind, tag = _resolve_identity(kind, tag, arm, run_name)
    return preset_dir_for(raw_runs, period, preset, kind, tag) / run_name


def legacy_run_dirs_for(raw_runs, run_name, period, preset, arm=None):
    """Every place this run could sit under the two earlier layouts.

    2026-09-10 layout:  arms/<arm>/<period>/<preset>/<run>  for an arm,
                        <period>/<preset>/sweeps/<family>/<run>  for a cell,
                        <period>/<preset>/<run>  otherwise.
    Before that:        <arm>/<period>/<preset>/<run>  and  <period>/<preset>/<run>.

    Args:
        raw_runs: The output/raw_runs root.
        run_name: The run's directory name.
        period: "1984_2004" or (1984, 2004).
        preset: Source/sink preset.
        arm: The legacy arm the run was filed under; None/"" for calibration.

    Returns:
        Candidate Paths, most recent layout first. Not checked for existence.
    """
    root = Path(raw_runs)
    per = period_component(period)
    arm = (arm or "").strip()
    if arm == CALIBRATION_ARM:
        arm = ""
    family = sweep_family(run_name)
    out = []
    if arm:
        out.append(root / ARMS_DIR / arm / per / preset / run_name)
        out.append(root / arm / per / preset / run_name)
    else:
        if family:
            out.append(root / per / preset / SWEEPS_DIR / family / run_name)
        out.append(root / per / preset / run_name)
    return out


def _legacy_arm_for(kind, tag, run_name):
    """The arm name a (kind, tag) run would have carried before 09-16."""
    if kind == MATRIX_KIND:
        return ""
    if kind == "sensitivity":
        return ""            # cells sat in the calibration arm, token-named
    for arm, (k, t) in LEGACY_ARMS.items():
        if (k, t) == (kind, tag):
            return arm
    return tag


def kinds_holding(raw_runs, run_name, period, preset):
    """Every (kind, tag) under which this run exists on disk, either layout.

    A run name describes the SCENARIO, so one name can legitimately exist
    under several tags -- the matrix run and the version pair made from it.
    This is what makes that discoverable rather than a surprise: it is what
    `find_run_dir` reports when the place it was asked for holds nothing.

    Args:
        raw_runs: The output/raw_runs root.
        run_name: The run's directory name.
        period: "1984_2004" or (1984, 2004).
        preset: Source/sink preset.

    Returns:
        Sorted list of (kind, tag) pairs. Empty if the name is nowhere under
        this period and preset.
    """
    root = Path(raw_runs)
    if not root.is_dir():
        return []
    per = period_component(period)
    found = set()
    # The purpose layout: kind folders, each holding tags (one or two levels)
    # or, for the matrix, periods directly.
    for kind, folder in KIND_DIR.items():
        top = root / folder
        if not top.is_dir():
            continue
        if kind == MATRIX_KIND:
            if (top / per / preset / run_name).is_dir():
                found.add((kind, ""))
            continue
        for child in sorted(top.iterdir()):
            if not child.is_dir() or _PERIOD_DIR.fullmatch(child.name):
                continue
            candidates = [child.name]
            if not any(_PERIOD_DIR.fullmatch(g.name) for g in child.iterdir()
                       if g.is_dir()):
                candidates = [f"{child.name}/{g.name}"
                              for g in sorted(child.iterdir()) if g.is_dir()]
            for tag in candidates:
                if (top / tag / per / preset / run_name).is_dir():
                    found.add((kind, tag))
    # The earlier layouts, translated.
    arms = [""]
    for top in ([root / ARMS_DIR] if (root / ARMS_DIR).is_dir() else []) + [root]:
        for child in sorted(top.iterdir()):
            if (not child.is_dir() or _PERIOD_DIR.fullmatch(child.name)
                    or child.name in KIND_DIR.values() or child.name == ARMS_DIR):
                continue
            arms.append(child.name)
            if not any(_PERIOD_DIR.fullmatch(g.name) for g in child.iterdir()
                       if g.is_dir()):
                arms.extend(f"{child.name}/{g.name}"
                            for g in sorted(child.iterdir()) if g.is_dir())
    for arm in arms:
        if any(p.is_dir() for p in
               legacy_run_dirs_for(root, run_name, per, preset, arm)):
            found.add(legacy_arm_to_kind_tag(arm) if arm else
                      ("sensitivity", sweep_family(run_name)) if sweep_family(run_name)
                      else (MATRIX_KIND, ""))
    return sorted(found)


def find_run_dir(raw_runs, run_name, period, preset, kind=MATRIX_KIND, tag="",
                 arm=None):
    """Locates a finished run, raising with what IS on disk if it is absent.

    KIND DEFAULTS TO THE MATRIX RATHER THAN SEARCHING. A search would let a
    figure silently draw a run forced or versioned differently from the one it
    names -- exactly what the kind and tag exist to prevent. Naming nothing
    means the matrix, which is also what every call site did before arms
    existed, so routing an existing script through this cannot change which
    run it reads.

    Args:
        raw_runs: The output/raw_runs root.
        run_name: The run's directory name.
        period: "1984_2004" or (1984, 2004).
        preset: Source/sink preset.
        kind: One of KINDS. Defaults to the matrix.
        tag: The experiment/version tag, or the axis of a sensitivity cell
            (derived from the name when not given).
        arm: LEGACY. A pre-09-16 arm name; translated, and also tried at its
            old path.

    Returns:
        The run directory as a Path, which exists.

    Raises:
        FileNotFoundError: If that directory is absent. The message names the
            places the run DOES exist, so a run filed elsewhere reads as "it
            is over there" rather than as "it was never made".
    """
    directory = run_dir_for(raw_runs, run_name, period, preset, kind, tag, arm)
    if directory.is_dir():
        return directory
    kind, tag = _resolve_identity(kind, tag, arm, run_name)
    legacy_arm = arm if arm else _legacy_arm_for(kind, tag, run_name)
    for candidate in legacy_run_dirs_for(raw_runs, run_name, period, preset,
                                         legacy_arm):
        if candidate.is_dir():
            return candidate

    elsewhere = [f"{k}:{t}" if t else k
                 for k, t in kinds_holding(raw_runs, run_name, period, preset)
                 if (k, t) != (kind, tag)]
    hint = (f"\n  It exists under: {', '.join(elsewhere)} -- pass kind=/tag= "
            f"to read one of those."
            if elsewhere else
            "\n  It exists nowhere; the run has not been made.")
    where = f"{kind}:{tag}" if tag else kind
    raise FileNotFoundError(
        f"no run directory for {run_name!r} in {where}.\n  {directory}" + hint)


def run_dir_for_index_row(raw_runs, row):
    """The run directory named by one run_index.csv row.

    The index carries `kind`, `tag`, `start_year`, `end_year` and
    `source_sink_preset` for exactly this: a row and a directory can be
    matched without either side reconstructing the other's spelling. A row
    from before 2026-09-16 carries `arm` instead, which is translated.

    Args:
        raw_runs: The output/raw_runs root.
        row: A mapping or pandas Series with run_name, start_year, end_year,
            source_sink_preset, and kind/tag (or the legacy arm).

    Returns:
        The run directory as a Path, resolved through find_run_dir so an
        unmigrated tree still answers.
    """
    kind = row["kind"] if "kind" in row and str(row["kind"]).strip() else None
    tag = row["tag"] if "tag" in row and str(row["tag"]).strip() else ""
    arm = None
    if kind is None:
        arm = row["arm"] if "arm" in row else CALIBRATION_ARM
        if isinstance(arm, float):           # NaN from pandas
            arm = CALIBRATION_ARM
    if isinstance(tag, float):
        tag = ""
    return find_run_dir(
        raw_runs, row["run_name"],
        (int(row["start_year"]), int(row["end_year"])),
        row["source_sink_preset"], kind or MATRIX_KIND, tag, arm)


def arm_component(arm):
    """LEGACY. The path component a pre-09-16 arm contributed.

    Kept only so an old caller importing it still imports. New code files by
    kind and tag; see preset_dir_for.
    """
    arm = (arm or CALIBRATION_ARM).strip()
    return "" if arm == CALIBRATION_ARM else check_tag(arm)


def sweep_component(run_name):
    """LEGACY. The sweeps/<family> component of the 2026-09-10 layout."""
    family = sweep_family(run_name)
    return f"{SWEEPS_DIR}/{family}" if family else ""


def legacy_run_dir_for(raw_runs, run_name, period, preset, arm=CALIBRATION_ARM):
    """LEGACY. The flat pre-2026-09-10 location of a run."""
    return legacy_run_dirs_for(raw_runs, run_name, period, preset, arm)[-1]


def arms_holding(raw_runs, run_name, period, preset):
    """LEGACY. Old arm names this run exists under; see kinds_holding."""
    out = []
    for kind, tag in kinds_holding(raw_runs, run_name, period, preset):
        out.append(_legacy_arm_for(kind, tag, run_name) or CALIBRATION_ARM)
    return sorted(set(out))


# =============================================================================
# THE DERIVED INDEX
# =============================================================================
# run_index.csv is a RESTATEMENT of every run's metadata JSON in one table, so
# a question across runs is one read. Since 2026-09-16 no run appends to it:
# each run writes the row it would have appended INTO its metadata, under
# "index row", and the runner then calls rebuild_run_index, which regenerates
# the whole file from every metadata on disk and replaces it atomically. Two
# runs finishing at once both rebuild the same complete table, so the last
# writer wins nothing -- which is what lets runs be concurrent. Rows for runs
# that predate the "index row" section are carried over from the existing
# file by their old key, so nothing is lost in the changeover.

INDEX_KEY = ("run_name", "kind", "tag")
_LEGACY_INDEX_KEY = ("run_name", "Hs_m", "arm")
INDEX_SECTION = "index row"
_LEADING_COLUMNS = ("run_name", "kind", "tag", "status", "timestamp",
                    "start_year", "end_year", "source_sink_preset")


def _read_index_rows(index_path):
    """The index as a list of dicts of strings; [] if absent."""
    import csv
    index_path = Path(index_path)
    if not index_path.is_file():
        return []
    with open(index_path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _kind_tag_from_path(raw_runs, run_dir):
    """(kind, tag) read off a run directory's place in the tree.

    Either layout: a purpose folder says so directly; an old arms/<arm>/ path
    or a token-named cell is translated as legacy_arm_to_kind_tag would.
    """
    rel = Path(run_dir).relative_to(raw_runs).parts
    if not rel:
        return MATRIX_KIND, ""
    for kind, folder in KIND_DIR.items():
        if rel[0] == folder:
            if kind == MATRIX_KIND:
                return kind, ""
            # <folder>/<tag parts...>/<period>/<preset>/<run>
            i = next((k for k, part in enumerate(rel)
                      if _PERIOD_DIR.fullmatch(part)), None)
            tag = "/".join(rel[1:i]) if i else ""
            return kind, tag
    if rel[0] == ARMS_DIR:
        i = next((k for k, part in enumerate(rel)
                  if _PERIOD_DIR.fullmatch(part)), len(rel))
        return legacy_arm_to_kind_tag("/".join(rel[1:i]))
    if _PERIOD_DIR.fullmatch(rel[0]):
        family = sweep_family(rel[-1])
        return ("sensitivity", family) if family else (MATRIX_KIND, "")
    # Pre-09-10 loose arm: <arm parts>/<period>/...
    i = next((k for k, part in enumerate(rel)
              if _PERIOD_DIR.fullmatch(part)), len(rel))
    return legacy_arm_to_kind_tag("/".join(rel[:i]))


def _legacy_arm_from_path(raw_runs, run_dir):
    """The `arm` a pre-09-16 index row spelled for a run at this path.

    "calibration" for the unscoped tree (matrix runs AND token-named sweep
    cells, either old layout); the arm folder(s) for arms/<arm>/ or a loose
    pre-09-10 arm; and for the purpose layout, the arm LEGACY_ARMS maps the
    (kind, tag) back to. This is what finds an old row for a moved run.
    """
    rel = Path(run_dir).relative_to(raw_runs).parts
    i = next((k for k, part in enumerate(rel) if _PERIOD_DIR.fullmatch(part)),
             len(rel))
    if not rel:
        return CALIBRATION_ARM
    if rel[0] == ARMS_DIR:
        return "/".join(rel[1:i]) or CALIBRATION_ARM
    if rel[0] in KIND_DIR.values():
        kind, tag = _kind_tag_from_path(raw_runs, run_dir)
        return _legacy_arm_for(kind, tag, rel[-1]) or CALIBRATION_ARM
    return "/".join(rel[:i]) or CALIBRATION_ARM


def run_status(kind, product, version, current_versions):
    """current / superseded / archived, for the index's `status` column.

    A matrix or sensitivity run on a topography that is no longer the
    product's CURRENT is superseded: comparing it with a run made today would
    attribute the pick difference to whatever the figure is about. A version
    or experiment run is judged against nothing -- it names its own inputs
    deliberately -- and an archived run says so by where it sits.
    """
    if kind == "archive":
        return "archived"
    if kind in (MATRIX_KIND, "sensitivity"):
        want = current_versions.get(product)
        if want and version and str(version) != str(want):
            return "superseded"
    return "current"


def rebuild_run_index(raw_runs, index_path=None, current_versions=None):
    """Regenerates run_index.csv from every run's metadata; returns the rows.

    Args:
        raw_runs: The output/raw_runs root.
        index_path: Where to write; default raw_runs/run_index.csv.
        current_versions: {topo_product: CURRENT dune-topo version}, for the
            status column. None leaves status "current" for everything but
            the archive; pass hat_topo_version's answer to get supersession.

    Returns:
        The list of row dicts written, in file order.

    Notes:
        A run whose metadata has no "index row" section (made before
        2026-09-16) keeps the row the existing file holds for it, found by the
        old (run_name, Hs_m, arm) key and given kind/tag from its path. A run
        with neither is indexed by its identity columns alone, so it is not
        lost. Rows whose run is gone from disk are dropped here; recording
        them is HAT_index_runs.py's job (retired_runs.csv), which calls this.
    """
    import csv
    import os
    import tempfile
    raw_runs = Path(raw_runs)
    index_path = Path(index_path) if index_path else raw_runs / RUN_INDEX_FILENAME
    current_versions = current_versions or {}

    existing = _read_index_rows(index_path)
    by_new_key = {tuple(r.get(k, "") for k in INDEX_KEY): r for r in existing
                  if r.get("kind")}
    # The old key only means something in a file that still HAS the arm
    # column. On a rebuilt file every row's arm reads "", so three runs of
    # one name collapse onto one old key and the last wins -- which handed
    # the 1996 matrix run a version row's skill on 2026-09-16.
    by_old_key = {tuple(r.get(k, "") for k in _LEGACY_INDEX_KEY): r
                  for r in existing if "arm" in r}

    rows = []
    for meta in sorted(raw_runs.rglob("*_run_metadata.json")):
        try:
            data = json.loads(meta.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            continue
        run_dir = meta.parent
        kind, tag = _kind_tag_from_path(raw_runs, run_dir)
        ident = data.get("identity", {})
        row = None
        if INDEX_SECTION in data:
            row = {k: ("" if v is None else v)
                   for k, v in data[INDEX_SECTION].items()}
        else:
            name = str(ident.get("run_name") or run_dir.name)
            hs = str(data.get("wave climate", {}).get("wave_height_m", ""))
            row = by_new_key.get((name, kind, tag))
            row = dict(row) if row else None
            if row is None:
                old_arm = _legacy_arm_from_path(raw_runs, run_dir)
                for probe in ((name, hs, old_arm), (name, hs, ""),
                              (name, "", old_arm)):
                    if probe in by_old_key:
                        row = dict(by_old_key[probe])
                        break
            if row is None:
                row = {
                    "run_name": name, "timestamp": ident.get("timestamp", ""),
                    "start_year": data.get("period", {}).get("start_year", ""),
                    "end_year": data.get("period", {}).get("end_year", ""),
                    "source_sink_preset": data.get("source/sink", {}).get("preset", ""),
                    "Hs_m": hs,
                    "topo_product": ident.get("topo_product", ""),
                    "topo_dune_version": ident.get("topo_dune_version", "")}
        row.pop("arm", None)
        row["kind"], row["tag"] = kind, tag
        row["status"] = run_status(kind, row.get("topo_product", ""),
                                   row.get("topo_dune_version", ""),
                                   current_versions)
        rows.append(row)
        # BACKFILL: a legacy run's row, once found, is written into its own
        # metadata JSON so the next rebuild derives it from the run and not
        # from whatever the index file happens to hold. Without this, a
        # rebuild that read a file lacking the arm column matched three
        # runs of one name to one row (2026-09-16). The .txt is left alone.
        if INDEX_SECTION not in data and row.get("rmse_interior_m_yr", "") != "":
            data[INDEX_SECTION] = {k: v for k, v in row.items()
                                   if k not in ("kind", "tag", "status")}
            try:
                meta.write_text(json.dumps(data, indent=2, ensure_ascii=False)
                                + "\n", encoding="utf-8")
            except OSError:
                pass

    # Column order: identity first, then everything else in first-seen order,
    # so the file stays readable as columns accumulate.
    columns = [c for c in _LEADING_COLUMNS]
    for row in rows:
        for column in row:
            if column not in columns:
                columns.append(column)
    rows.sort(key=lambda r: tuple(str(r.get(k, "")) for k in
                                  ("start_year", "kind", "tag", "source_sink_preset",
                                   "run_name")))

    # csv, not pandas: pandas round-trips every float through repr, which once
    # rewrote unrelated rows. Atomic: written beside, then replaced, so a
    # reader never sees a half-written file and two writers cannot interleave.
    index_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".run_index_", suffix=".csv",
                               dir=str(index_path.parent))
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns,
                                    extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow({k: ("" if v is None else v) for k, v in row.items()})
        os.replace(tmp, index_path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    return rows


def load_run_index(index_path):
    """run_index.csv as a DataFrame with kind/tag/status guaranteed present.

    A file from before 2026-09-16 has `arm` instead; it is translated so a
    reader written against the new columns works on either.
    """
    frame = pd.read_csv(index_path, dtype=str, keep_default_na=False)
    if "kind" not in frame.columns:
        pairs = [legacy_arm_to_kind_tag(a) for a in frame.get("arm", "")]
        frame["kind"] = [k for k, _ in pairs]
        frame["tag"] = [t for _, t in pairs]
        # A token-named cell sat in the calibration arm; say what it is.
        fam = frame["run_name"].map(sweep_family)
        cells = (frame["kind"] == MATRIX_KIND) & (fam != "")
        frame.loc[cells, "kind"] = "sensitivity"
        frame.loc[cells, "tag"] = fam[cells]
    if "status" not in frame.columns:
        frame["status"] = "current"
    return frame


def run_dir_contents(run_dir):
    """Lists what a run directory holds that counts as output.

    The single definition of "this directory holds a result", so the guard
    below and anything reporting on it cannot disagree about whether a stray
    .gitkeep means the directory is occupied.

    Args:
        run_dir: Directory to inspect. Need not exist.

    Returns:
        Sorted list of Paths, empty if the directory is absent or holds
        nothing but ignorable files.
    """
    run_dir = Path(run_dir)
    if not run_dir.exists():
        return []
    return sorted(p for p in run_dir.glob("*")
                  if p.name not in _IGNORABLE_NAMES)


def guard_run_dir(run_dir, overwrite=False):
    """Refuses to write into a run directory that already holds output.

    Called before the model is stepped, not after, so a name collision costs
    nothing. Overwriting is what makes a paired comparison quietly wrong: the
    baseline the groin run measures against is resolved by directory name, so
    replacing one run's outputs silently redefines the other run's answer.
    That is why the default refuses rather than asks.

    With overwrite=True the directory is EMPTIED, not written over in place.
    Writing over in place leaves behind any file the previous run produced and
    the new one does not -- a difference GIF from a run that had a paired
    baseline, a road summary from a run whose roadway manager was on -- and a
    leftover file in a run directory is indistinguishable from a current one.
    Emptying first means the directory holds exactly one run's output.

    Args:
        run_dir: Directory this run will write to.
        overwrite: True to empty the directory and reuse it. Intended for
            iterating on one scenario -- tweak a value, re-run, read the
            figures, tweak again -- where only the current state matters. The
            previous run's output is deleted and is NOT recoverable.

    Returns:
        The run_dir as a Path, created if it did not exist.

    Raises:
        RuntimeError: If the directory holds output and overwrite is False,
            or if overwrite is True and the directory holds a subdirectory
            that is not one of the run's own (run_layout.SUBFOLDERS plus
            gif_frames). An unexpected subdirectory means this is not the
            directory it is taken to be, and deleting its contents is
            refused rather than guessed at.
    """
    run_dir = Path(run_dir)
    existing = run_dir_contents(run_dir)

    if existing and not overwrite:
        listed = ", ".join(sorted(p.name for p in existing)[:3])
        raise RuntimeError(
            f"{run_dir.name} already holds {len(existing)} file(s) "
            f"({listed}{', ...' if len(existing) > 3 else ''}).\n"
            f"  This scenario has been run before. Set OVERWRITE = True to "
            f"replace it, delete the directory, or change a switch so the "
            f"derived name differs.\n"
            f"  {run_dir}")

    if existing:
        # Refused rather than handled: every file a run writes is flat, so a
        # subdirectory here means run_dir is not pointing where it is thought
        # to be -- a half-built path, a period directory, the output root. The
        # cost of guessing wrong is a recursive delete of someone's runs.
        # A run directory holds a KNOWN set of subfolders (run_layout's
        # figures/ animations/ tables/, plus the gif_frames scratch) and
        # nothing else. Any OTHER subdirectory still means run_dir is not
        # pointing where it is thought to be -- a half-built path, a period
        # directory, the output root -- and the cost of guessing wrong is a
        # recursive delete of someone's runs, so that stays refused. Before
        # 2026-09-10 EVERY subdirectory was refused, which stopped OVERWRITE
        # working at all once the layout gained folders.
        known = set(SUBFOLDERS) | {"gif_frames"}
        unknown = sorted(p.name for p in existing
                         if p.is_dir() and p.name not in known)
        if unknown:
            raise RuntimeError(
                f"{run_dir.name} holds subdirector"
                f"{'ies' if len(unknown) > 1 else 'y'} "
                f"({', '.join(unknown[:3])}) that a run directory never "
                f"does.\n"
                f"  Refusing to empty it -- check that RUN_DIR points where "
                f"you think, then remove it by hand if it is really the "
                f"directory you meant.\n"
                f"  A run's own subfolders are: {', '.join(sorted(known))}.\n"
                f"  {run_dir}")
        for path in existing:
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()

    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def render_metadata_text(sections, header):
    """Renders the metadata sections as the human-readable .txt.

    Args:
        sections: Ordered mapping of section name to a mapping of key to
            either a value, or a (value, comment) tuple. The comment is shown
            in the .txt and dropped from the .json.
        header: Comment lines placed at the top of the file, without the
            leading "# ".

    Returns:
        The file contents as a string.
    """
    # One column width for the whole file, wide enough for the longest key, so
    # the "=" stays aligned across sections instead of jogging in and out.
    width = max([22] + [len(key) + 1
                        for entries in sections.values() for key in entries])

    lines = [f"# {line}" for line in header]
    for name, entries in sections.items():
        lines += ["", f"[{name}]"]
        for key, entry in entries.items():
            value, comment = entry if isinstance(entry, tuple) else (entry, None)
            rendered = f"{key:<{width}}= {value}"
            lines.append(f"{rendered}   # {comment}" if comment else rendered)
    return "\n".join(lines) + "\n"


def _json_safe(value):
    """Converts numpy scalars and Paths to types json.dump can write."""
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def write_run_metadata(run_dir, run_name, sections, header):
    """Writes the run metadata as both .txt and .json.

    Args:
        run_dir: Directory to write into.
        run_name: Run name, used for the filenames.
        sections: As accepted by render_metadata_text.
        header: Comment lines for the top of the .txt.

    Returns:
        A (txt_path, json_path) tuple of Paths.
    """
    run_dir = Path(run_dir)
    txt_path = run_dir / f"{run_name}_run_metadata.txt"
    json_path = run_dir / f"{run_name}_run_metadata.json"

    txt_path.write_text(render_metadata_text(sections, header), encoding="utf-8")

    payload = {
        name: {key: _json_safe(entry[0] if isinstance(entry, tuple) else entry)
               for key, entry in entries.items()}
        for name, entries in sections.items()
    }
    json_path.write_text(json.dumps(payload, indent=2, default=str) + "\n",
                         encoding="utf-8")
    return txt_path, json_path


def skill_vs_target(change_rate, target_table, geometry,
                    interior_margin=1):
    """Model-minus-observed skill over the real domains.

    Reported over two spans, deliberately. The end domains carry the locked
    source/sink values (tens of m/yr), so an island-wide RMSE for a calibBE or
    edgeBE run is dominated by two domains that were pinned rather than
    predicted -- and a zeroBE run has no such term. Comparing presets on the
    island-wide number alone would mostly compare the boundary treatment.
    The interior number excludes them, which is the same exclusion the
    source/sink QC plot makes for its zoom panel.

    Args:
        change_rate: Padded per-domain model rate array, m/yr, already
            sign-flipped so (+) is seaward.
        target_table: DataFrame with gis_domain and target_lrr_m_yr columns
            (COASTSAT_TARGET).
        geometry: DomainGeometry describing the padded array.
        interior_margin: Domains excluded from each end for the interior
            metrics. 1 drops the two locked end domains.

    Returns:
        A dict of mean bias and RMSE in m/yr, island-wide and interior, plus
        the domain count each was computed over. NaN where no domains remain.
    """
    model = np.asarray(change_rate)[
        geometry.start_real_index:geometry.end_real_index]
    gis_ids = np.arange(geometry.first_gis_id, geometry.last_gis_id + 1)

    target_by_gis = dict(zip(target_table["gis_domain"],
                             target_table["target_lrr_m_yr"]))
    target = np.array([target_by_gis.get(g, np.nan) for g in gis_ids],
                      dtype=float)

    def _metrics(mask):
        residual = model[mask] - target[mask]
        residual = residual[np.isfinite(residual)]
        if residual.size == 0:
            return np.nan, np.nan, 0
        return (float(residual.mean()),
                float(np.sqrt((residual ** 2).mean())),
                int(residual.size))

    everywhere = np.ones(gis_ids.size, dtype=bool)
    interior = everywhere.copy()
    if interior_margin > 0:
        interior[:interior_margin] = False
        interior[-interior_margin:] = False

    bias_all, rmse_all, n_all = _metrics(everywhere)
    bias_in, rmse_in, n_in = _metrics(interior)
    return {
        "mean_bias_m_yr": bias_all,
        "rmse_m_yr": rmse_all,
        "n_domains": n_all,
        "mean_bias_interior_m_yr": bias_in,
        "rmse_interior_m_yr": rmse_in,
        "n_domains_interior": n_in,
    }


def append_run_index(index_path, row, key="run_name"):
    """Adds one run to the cross-run index CSV, replacing any earlier row.

    Replaces rather than appends on a repeat run name so the index tracks what
    is currently on disk. A run directory holds exactly one result, so two
    index rows for one name could only ever mean one of them is stale.

    New columns are unioned in, so adding a field later does not invalidate an
    index written before it existed -- older rows get NaN for it.

    Args:
        index_path: Path to the index CSV. Created if absent.
        row: Mapping of column name to value for this run.
        key: Column, or sequence of columns, identifying a run uniquely.

    Returns:
        The full index as a DataFrame, as written.

    Note:
        The index is a DERIVED view of the runs on disk: every value in it is
        restated from a run's own metadata. This keeps it in step as runs are
        made; `HAT_index_runs.py` rebuilds it from the runs, checks it against
        them, and records anything whose run has been deleted.
    """
    index_path = Path(index_path)
    new = pd.DataFrame([row])
    # A COMPOSITE key is allowed because the run name no longer identifies a
    # run on its own: forcing that is not part of the scenario -- Hs -- scopes
    # the output DIRECTORY instead of adding a name token, so two runs can share
    # a name and differ in what they were forced with. Replacing on name alone
    # would silently drop one of them from the index.
    keys = (key,) if isinstance(key, str) else tuple(key)

    if index_path.exists():
        # AS TEXT. Parsing the existing rows into pandas and writing them
        # back rewrites every float through repr, which on 2026-09-10 silently
        # truncated the last digit of five columns in a row this call was not
        # touching. Only the row being added should change.
        existing = pd.read_csv(index_path, dtype=str, keep_default_na=False)
        usable = [k for k in keys if k in existing.columns]
        if usable:
            same = pd.Series(True, index=existing.index)
            for k in usable:
                same &= existing[k].astype(str) == str(row.get(k))
            existing = existing[~same]
        combined = pd.concat([existing, new], ignore_index=True)
    else:
        combined = new

    # Stable ordering: the identity columns first, then whatever else exists,
    # so the file stays readable as columns accumulate.
    leading = [c for c in keys + ("timestamp", "start_year", "end_year")
               if c in combined.columns]
    combined = combined[leading + [c for c in combined.columns
                                   if c not in leading]]
    sort_on = [k for k in keys if k in combined.columns]
    if sort_on:
        combined = combined.sort_values(sort_on,
                                        kind="stable").reset_index(drop=True)

    index_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(index_path, index=False)
    return combined


def timestamp():
    """Current local time, formatted for metadata and the index."""
    return f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S}"
