r"""
build_island_offset.py -- a digitised dune line in, a model input out
==============================================================================
One command from a geojson under 2-brie-offset/dunelines/ to
the padded 120-domain file a hindcast start reads, following the same steps
that were run by hand until 2026-09-15 (Hannah: "one driver that chains the
two steps"):

    1. duneline_to_raw_offsets.py   line x 100 m transects -> per-transect
                                    stations from the offshore datum, written
                                    to raw_offsets/<vintage>_duneline_offset_raw.csv
    2. island_offset_hybrid.py      first row per transect, domain mean, zeroed
                                    on the minimum, padded to 120 -> <year>/v<n>/
    3. HAT_compare_offset_versions  the new build against the previous one,
                                    in both frames (model, and fixed datum)
    4. CURRENT <- v<n>              unless --no-current
    5. PROVENANCE.md                written from the run, not by hand

The two scripts are unchanged and still run on their own.

VINTAGE VS PERIOD YEAR
    The geojson is named for the IMAGERY vintage of the line (duneline_1997);
    --year is the hindcast start it serves (1996). The pairing is
    hat_topo_version.DUNE_LINE_FOR_YEAR, the only place it is spelled, and
    this driver refuses a pair the table does not hold rather than guess.

VERSIONS
    Every build is a version: <year>/v1/, v2/, ... and a CURRENT file naming
    the one the runner reads (hatteras_site_config._island_offset_file).
    A version is never overwritten; ask for the next one. Each version folder
    keeps a copy of the raw file it was built from, so it can be rebuilt with
    island_offset_hybrid.py --raw-file after the vintage's raw has moved on.

USAGE
    python build_island_offset.py --duneline duneline_1984.geojson --year 1984
    python build_island_offset.py --duneline duneline_1997_v2.geojson --year 1996 \
        --compare-with v1
    python build_island_offset.py --duneline duneline_2009.geojson --year 2010   # once in the table

    --version vN       name the version (default: the next free one)
    --compare-with X   a version folder under <year>/ to compare against
                       (default: the previous version, else a superseded_*
                       folder if one exists, else no comparison)
    --validate-against a GIS export of the same line, passed to step 1
    --no-current       build but leave CURRENT as it is
==============================================================================
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = next(_p for _p in Path(__file__).resolve().parents
                    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from hat_topo_version import (BRIE_ROOT, DUNE_LINE_FOR_YEAR, DUNELINE_DIR,  # noqa: E402
                              RAW_OFFSET_DIR, dune_line_for_year, dune_raw_file)

HERE = Path(__file__).resolve().parent
STEP_RAW = HERE / "duneline_to_raw_offsets.py"
STEP_PAD = HERE / "island_offset_hybrid.py"
STEP_CMP = HERE / "HAT_compare_offset_versions.py"

# geojson properties copied into the provenance when present; imagery_date is
# the one every new line should carry (the 1984 and 2004 dates had to be
# recovered from a metadata file and from memory).
LINE_PROPS = ("feature_type", "year", "imagery_date", "source_type", "method",
              "editor", "edit_date", "notes")


def _echo(text):
    """Print a step's output on whatever console this is (Windows cp1252
    included) without dying on a character it cannot show."""
    enc = sys.stdout.encoding or "utf-8"
    sys.stdout.write(text.encode(enc, errors="replace").decode(enc))
    sys.stdout.flush()


def _run(cmd, log):
    """Run a step, echo it, keep its stdout for the provenance."""
    print("\n$ " + " ".join(str(c) for c in cmd))
    r = subprocess.run([sys.executable, *map(str, cmd)], capture_output=True,
                       text=True, encoding="utf-8", errors="replace",
                       env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    _echo(r.stdout)
    if r.returncode != 0:
        _echo(r.stderr)
        sys.exit(f"step failed: {cmd[0]}")
    log.append((Path(cmd[0]).name, r.stdout))
    return r.stdout


def _resolve_geojson(name_or_path):
    p = Path(name_or_path)
    if p.exists():
        return p.resolve()
    q = DUNELINE_DIR / name_or_path
    if q.exists():
        return q.resolve()
    sys.exit(f"{name_or_path}: not a path and not under {DUNELINE_DIR}")


def _line_properties(path):
    import geopandas as gpd
    g = gpd.read_file(path)
    if len(g) != 1:
        sys.exit(f"{path.name}: expected one feature, found {len(g)}")
    props = {k: g.iloc[0][k] for k in LINE_PROPS if k in g.columns}
    # ArcGIS writes dates as epoch milliseconds
    for k in ("edit_date", "imagery_date"):
        v = props.get(k)
        if isinstance(v, (int, float)) and v > 1e11:
            props[k] = dt.datetime.fromtimestamp(v / 1000, dt.timezone.utc).strftime("%Y-%m-%d")
    return props, str(g.crs)


def _vintage_of(path, props):
    if props.get("year") not in (None, "") and not pd.isna(props.get("year")):
        return int(props["year"])
    m = re.search(r"duneline_(\d{4})", path.name)
    if not m:
        sys.exit(f"{path.name}: no 'year' property and no duneline_<year> in the name")
    return int(m.group(1))


def _versions(year_dir):
    return sorted((p.name for p in year_dir.iterdir()
                   if p.is_dir() and re.fullmatch(r"v\d+", p.name)),
                  key=lambda v: int(v[1:])) if year_dir.is_dir() else []


def main(argv=None):
    ap = argparse.ArgumentParser(description="a dune line in, a model input out")
    ap.add_argument("--duneline", required=True)
    ap.add_argument("--year", type=int, required=True, help="hindcast start year")
    ap.add_argument("--version", default=None)
    ap.add_argument("--compare-with", default=None)
    ap.add_argument("--validate-against", default=None)
    ap.add_argument("--no-current", action="store_true")
    a = ap.parse_args(argv)

    line = _resolve_geojson(a.duneline)
    props, crs = _line_properties(line)
    vintage = _vintage_of(line, props)
    expected = dune_line_for_year(a.year)
    if vintage != expected:
        sys.exit(f"\n{line.name} is a {vintage} line; DUNE_LINE_FOR_YEAR pairs the "
                 f"{a.year} start with the {expected} line. Change the table in "
                 f"hat_topo_version.py if the pairing is wrong; the driver does not guess.\n")

    year_dir = BRIE_ROOT / str(a.year)
    year_dir.mkdir(parents=True, exist_ok=True)
    have = _versions(year_dir)
    version = a.version or f"v{(int(have[-1][1:]) + 1) if have else 1}"
    if not re.fullmatch(r"v\d+", version):
        sys.exit(f"--version must look like v3, not {version!r}")
    if (year_dir / version).exists():
        sys.exit(f"{year_dir / version} exists; a version is never overwritten. "
                 f"Have {have}; the next free one is the default.")
    compare_with = a.compare_with
    if compare_with is None:
        sup = sorted(p.name for p in year_dir.iterdir()
                     if p.is_dir() and p.name.startswith("superseded_"))
        compare_with = have[-1] if have else (sup[-1] if sup else None)
    if compare_with and not (year_dir / compare_with).is_dir():
        sys.exit(f"--compare-with {compare_with}: no {year_dir / compare_with}")

    print(f"line     {line.name}  (vintage {vintage}, {crs})")
    print(f"start    {a.year}  ->  {year_dir.name}/{version}/"
          + (f"   compared with {compare_with}" if compare_with else ""))

    log = []
    raw_path = dune_raw_file(vintage)

    # 1. line -> raw stations
    cmd = [STEP_RAW, "--duneline", line, "--out", raw_path]
    if a.validate_against:
        cmd += ["--validate-against", a.validate_against]
    _run(cmd, log)
    raw = pd.read_csv(raw_path)
    n_missing = int(raw["ORIG_LEN"].isna().sum())
    n_multi = int((raw["n_crossings"] > 1).sum())
    per_dom = raw.groupby("domain_id")["ORIG_LEN"].count()

    # 2. raw -> padded model input
    out_pad = _run([STEP_PAD, "--year", a.year, "--version", version], log)
    m = re.search(r"Baseline distance = ([\d.]+) m", out_pad)
    baseline = float(m.group(1)) if m else float("nan")
    ver_dir = year_dir / version
    dom_means = raw.drop_duplicates(["domain_id", "LineID"]).groupby("domain_id")["ORIG_LEN"].mean()
    zero_domain = int(dom_means.idxmin())

    # each version keeps the raw it was built from
    shutil.copy2(raw_path, ver_dir / raw_path.name)

    # 3. against the previous build
    cmp_text = ""
    if compare_with:
        raw_a = year_dir / compare_with / raw_path.name
        if not raw_a.exists():
            alt = RAW_OFFSET_DIR / "superseded_20260915_gis_exports" / raw_path.name
            raw_a = alt if alt.exists() else None
        cmd = [STEP_CMP, "--year", a.year, "--a", compare_with, "--b", version]
        if raw_a is not None:
            cmd += ["--raw-a", raw_a, "--raw-b", ver_dir / raw_path.name]
        cmp_text = _run(cmd, log)

    # 4. CURRENT
    current = year_dir / "CURRENT"
    if not a.no_current:
        current.write_text(version + "\n", encoding="utf-8")

    # 5. provenance, from the run
    stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    short = per_dom[per_dom < 5]
    props_rows = "\n".join(f"| {k} | {v} |" for k, v in props.items()) or "| (none) | the file carries no properties |"
    lines = [
        f"# {a.year} island offsets, {version}",
        "",
        f"Built {stamp} by `scripts/input_prep/2-brie-offset/build_island_offset.py` "
        f"from `{line.name}`" + (f" ({vintage} imagery, standing in for the {a.year} start "
                                 f"through `DUNE_LINE_FOR_YEAR`)" if vintage != a.year else
                                 f" (a {a.year} survey, no stand-in)") + ".",
        "",
        "## The line",
        "",
        "| property | value |",
        "|---|---|",
        f"| file | `2-brie-offset/dunelines/{line.name}` |",
        f"| crs | {crs} |",
        props_rows,
        "",
    ]
    if "imagery_date" not in props:
        lines += ["**No `imagery_date` property.** Add one to the geojson so the "
                  "date does not have to be recovered later (the 1984 date came "
                  "from a USGS metadata file, the 2004 date from memory).", ""]
    lines += [
        "## The intersection (step 1)",
        "",
        f"`duneline_to_raw_offsets.py` against `transects/transects_100m.geojson`: "
        f"{len(raw)} transects in GIS 1-90, {n_missing} with no crossing, "
        f"{n_multi} crossed more than once"
        + (f"; domains with fewer than five transects: {short.to_dict()}" if len(short) else "")
        + f". Written to `raw_offsets/{raw_path.name}` (the vintage's current raw) and "
        f"copied here as `{raw_path.name}`.",
        "",
        "## The model input (step 2)",
        "",
        f"`island_offset_hybrid.py --year {a.year} --version {version}`: first row per "
        f"transect, mean of the transects in each domain, zeroed on the most seaward "
        f"domain (GIS {zero_domain}, {baseline:.2f} m from the datum), padded to 120 "
        f"with the slope-and-bridge buffer. Files: `Island_Dune_Offsets_{a.year}_"
        f"PADDED_120.csv` (read by the model), `_CASCADE_Input.csv`, "
        f"`_CASCADE_Input_unpadded.csv`, `_buffer_diagnostic.png`.",
        "",
    ]
    if compare_with:
        lines += [
            f"## Against {compare_with} (step 3)",
            "",
            "```",
            cmp_text.strip(),
            "```",
            "",
            f"`offset_{a.year}_{compare_with}_vs_{version}.csv` and `.png` beside this "
            f"file (PDF and caption under `supporting/`). In the fixed-datum frame a "
            f"positive difference is the line moved LANDWARD.",
            "",
        ]
    lines += [
        "## CURRENT",
        "",
        (f"`../CURRENT` = `{version}` since this build." if not a.no_current else
         f"Built with `--no-current`; `../CURRENT` still names "
         f"`{current.read_text().strip() if current.exists() else '(none)'}`."),
        " `hatteras_site_config._island_offset_file({0})` resolves it; env "
        f"`HAT_OFFSET_VERSION_{a.year}` outranks the file for one run.".format(a.year),
        "",
        "## Rebuild",
        "",
        "```",
        f"python scripts/input_prep/2-brie-offset/island_offset_hybrid.py --year {a.year} "
        f"--version {version} --raw-file data/hatteras_init/2-brie-offset/{a.year}/{version}/{raw_path.name}",
        "```",
        "",
        "## Step output",
        "",
    ]
    for name, text in log:
        lines += [f"### {name}", "", "```", text.strip(), "```", ""]
    (ver_dir / "PROVENANCE.md").write_text("\n".join(lines), encoding="utf-8")

    # the year's index: a build log the driver appends to
    idx = year_dir / "PROVENANCE.md"
    row = (f"| `{version}` | {stamp[:10]} | `{line.name}` | {vintage} | "
           f"GIS {zero_domain} | {compare_with or '—'} | "
           f"{'CURRENT' if not a.no_current else ''} |")
    header = ("\n## Builds\n\nWritten by `build_island_offset.py`, one row per build; "
              "each version's own `PROVENANCE.md` has the detail.\n\n"
              "| version | built | line | vintage | zero domain | compared with | |\n"
              "|---|---|---|---|---|---|---|\n")
    if idx.exists():
        text = idx.read_text(encoding="utf-8")
        if "\n## Builds\n" not in text:
            text = text.rstrip("\n") + "\n" + header
        text = text.rstrip("\n") + "\n" + row + "\n"
    else:
        text = (f"# {a.year} island offsets — version index\n\n"
                f"`CURRENT` names the build every reader takes; "
                f"`hatteras_site_config._island_offset_file({a.year})` resolves it "
                f"(env `HAT_OFFSET_VERSION_{a.year}` outranks the file).\n" + header + row + "\n")
    idx.write_text(text, encoding="utf-8")

    print(f"\n{year_dir.name}/{version}/  zero domain GIS {zero_domain}  "
          f"{'CURRENT' if not a.no_current else 'not current'}")
    print(f"provenance {ver_dir / 'PROVENANCE.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
