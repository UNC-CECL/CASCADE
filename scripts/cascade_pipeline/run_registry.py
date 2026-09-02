"""Run provenance, output-directory guarding, and the cross-run index.

The hindcast is run as a matrix -- period x source/sink preset x groin -- and
the three failure modes that costs are all bookkeeping ones:

  * a re-run silently overwriting the outputs of the run it was meant to be
    compared against (guard_run_dir),
  * a finished run whose metadata cannot say which code produced it, because
    the sandbox flag or the extractor version was flipped days earlier
    (git_provenance, and the [identity] section callers pass),
  * twelve runs whose results can only be compared by opening twelve
    hand-formatted text files (append_run_index).

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
from pathlib import Path

import numpy as np
import pandas as pd

# Files whose presence means a run directory holds real output. A directory
# containing only these is treated as empty, so a stray .gitkeep or an
# editor's .DS_Store does not block a run.
_IGNORABLE_NAMES = {".gitkeep", ".gitignore", ".DS_Store", "Thumbs.db"}

RUN_INDEX_FILENAME = "run_index.csv"


# Paths a run WRITES BACK into the repository, so their being modified says
# nothing about whether the run's code was committed. See git_provenance.
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

CALIBRATION_ARM = "calibration"

# A period directory is exactly <4 digits>_<4 digits>. Anything else directly
# under raw_runs/ is an arm, so the two levels can be told apart without a
# registry of arm names -- an arm is created by setting HAT_ARM_TAG, and no
# hardcoded list of them would stay current.
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


def arm_component(arm):
    """The path component an arm contributes, which is "" for calibration.

    Args:
        arm: Arm name, as the `arm` column of run_index.csv spells it. None
            and "" are read as the calibration arm.

    Returns:
        The arm name, or "" for the calibration arm.

    Raises:
        ValueError: If the arm is not a single path component. The value
            reaches here from an environment variable and is joined onto the
            output root.
    """
    arm = (arm or CALIBRATION_ARM).strip()
    if arm == CALIBRATION_ARM:
        return ""
    if "/" in arm or "\\" in arm or arm.startswith("."):
        raise ValueError(
            f"arm {arm!r} must be a single path component -- it is joined "
            f"onto the raw_runs root and must not escape it.")
    return arm


def preset_dir_for(raw_runs, period, preset, arm=CALIBRATION_ARM):
    """The directory holding every run of one period, preset and arm.

    This is the runner's OUTPUT_BASE_DIR. It is the level anything that
    ENUMERATES runs works at -- the scenario grid, the relocation comparison,
    the source/sink calibration -- as against `run_dir_for`, which is for a run
    already named.

    Args:
        raw_runs: The output/raw_runs root.
        period: "1984_2004" or (1984, 2004).
        preset: Source/sink preset, e.g. "calibBE".
        arm: Forcing arm. The default is the calibration arm, which
            contributes no path component.

    Returns:
        The preset directory as a Path. Does not check it exists.
    """
    root = Path(raw_runs)
    tag = arm_component(arm)
    base = root / tag if tag else root
    return base / period_component(period) / preset


def run_dir_for(raw_runs, run_name, period, preset, arm=CALIBRATION_ARM):
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
        arm: Forcing arm. The default is the calibration arm, which
            contributes no path component.

    Returns:
        The run directory as a Path.
    """
    return preset_dir_for(raw_runs, period, preset, arm) / run_name


def arms_holding(raw_runs, run_name, period, preset):
    """Every arm under which this run exists on disk.

    A run name describes the SCENARIO and an arm describes the FORCING, so one
    name can legitimately exist in several arms -- and three currently exist in
    four each. This is what makes that discoverable rather than a surprise: it
    is what `find_run_dir` reports when the arm it was asked for holds nothing.

    Args:
        raw_runs: The output/raw_runs root.
        run_name: The run's directory name.
        period: "1984_2004" or (1984, 2004).
        preset: Source/sink preset.

    Returns:
        Sorted list of arm names, using CALIBRATION_ARM for the unscoped tree.
        Empty if the name is nowhere under this period and preset.
    """
    root = Path(raw_runs)
    if not root.is_dir():
        return []
    candidates = [CALIBRATION_ARM] + sorted(
        child.name for child in root.iterdir()
        if child.is_dir() and not _PERIOD_DIR.fullmatch(child.name))
    return [arm for arm in candidates
            if run_dir_for(root, run_name, period, preset, arm).is_dir()]


def find_run_dir(raw_runs, run_name, period, preset, arm=CALIBRATION_ARM):
    """Locates a finished run, raising with what IS on disk if it is absent.

    ARM DEFAULTS TO CALIBRATION RATHER THAN SEARCHING. A search would let a
    figure silently draw a run forced at a wave climate other than the one it
    names -- exactly what the arm component exists to prevent -- and with three
    names currently present in four arms each it would have to guess between
    them. Naming no arm means the calibration arm, which is also what every
    call site did before arms existed, so routing an existing script through
    this cannot change which run it reads.

    Args:
        raw_runs: The output/raw_runs root.
        run_name: The run's directory name.
        period: "1984_2004" or (1984, 2004).
        preset: Source/sink preset.
        arm: Forcing arm to read from. Defaults to the calibration arm.

    Returns:
        The run directory as a Path, which exists.

    Raises:
        FileNotFoundError: If that directory is absent. The message names the
            arms the run DOES exist under, so an arm-scoped run reads as "it
            is over there" rather than as "it was never made".
    """
    directory = run_dir_for(raw_runs, run_name, period, preset, arm)
    if directory.is_dir():
        return directory

    elsewhere = [a for a in arms_holding(raw_runs, run_name, period, preset)
                 if a != arm]
    hint = (f"\n  It exists under arm(s): {', '.join(elsewhere)} -- pass "
            f"arm=... to read one of those."
            if elsewhere else
            "\n  It exists under no arm; the run has not been made.")
    raise FileNotFoundError(
        f"no run directory for {run_name!r} in arm {arm!r}.\n  {directory}"
        + hint)


def run_dir_for_index_row(raw_runs, row):
    """The run directory named by one run_index.csv row.

    The index carries `arm`, `start_year`, `end_year` and `source_sink_preset`
    for exactly this: a row and a directory can be matched without either side
    reconstructing the other's spelling.

    Args:
        raw_runs: The output/raw_runs root.
        row: A mapping or pandas Series with run_name, start_year, end_year,
            source_sink_preset, and optionally arm.

    Returns:
        The run directory as a Path. NOT checked for existence -- a row whose
        arm disagrees with the tree is a bookkeeping fault to be seen, not
        smoothed over here.
    """
    arm = row["arm"] if "arm" in row else CALIBRATION_ARM
    return run_dir_for(
        raw_runs, row["run_name"],
        (int(row["start_year"]), int(row["end_year"])),
        row["source_sink_preset"], arm)


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
            or if overwrite is True and the directory holds a subdirectory. A
            run directory is flat, so a subdirectory means this is not the
            directory it is taken to be, and deleting its contents is refused
            rather than guessed at.
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
        subdirs = [p.name for p in existing if p.is_dir()]
        if subdirs:
            raise RuntimeError(
                f"{run_dir.name} holds subdirector"
                f"{'ies' if len(subdirs) > 1 else 'y'} "
                f"({', '.join(sorted(subdirs)[:3])}), which a run directory "
                f"never does.\n"
                f"  Refusing to empty it -- check that RUN_DIR points where "
                f"you think, then remove it by hand if it is really the "
                f"directory you meant.\n"
                f"  {run_dir}")
        for path in existing:
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
        existing = pd.read_csv(index_path)
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
