# How this repository is organised

Seven rules. They are not new: every one of them is a description of what the
project already does in the places that work, written down so the places that
do not can be brought into line, and so the next folder starts right.

`python scripts/repo_tools/hat_layout_check.py` reports every departure from them. It is
advisory and always exits zero — it produces a worklist, not an obstacle.

---

## 1. Code lives in `scripts/`, data in `data/`, products in `output/`

A script belongs beside other scripts; the files it reads and the files it
writes belong in the data or output trees. The test is simple: **if you
archived `data/hatteras_init/` and handed it to someone, could they see what
the model eats?** Until 2026-09-12 they could not, because the observed rate
targets — a model input, read on every run — lived in the scripts tree.

`README.md` is the one non-code file that belongs beside code.

**A `.py` extension is not a defence.** Data can wear one: the BE rate tables
are spelled as Python, and until 2026-09-18 four `hatteras_site_config_prebe_*.py`
snapshots of the solved BE field sat at the root of `scripts/`, written there by
the calibrate step before each pass. They looked like stray copies of the site
config, and the equivalent 2026-08-24 snapshot was discarded on exactly that
reading — taking the only record of the one-shot solve with it, and leaving
`plot_be_zones.py` undrawable for three weeks. They now file under
`data/hatteras_init/7-source-sink/2-calibrate/prebe/`.

The root of `scripts/` therefore holds **folders and `README.md`, nothing
else**, and the layout check reports any file that appears there. A module
belongs in `scripts/site_layer/`, a repo tool in `scripts/repo_tools/`, and
anything that is neither — a snapshot, an export, a copy left behind — belongs
in the data tree.

Only the root is judged. Deeper folders hold task scripts that are run by path
and imported by nobody, which is correct for them.

## 2. A survey is named for its year; a window for its span

```
4-mgmt-forcing/road_offset/dunestart_offset/derived/1996/      a survey: one moment
3-env-forcings/3-storms/hindcast_storms/1996_2010/     a window: an interval
```

A dune line, a road alignment or a lidar survey is a **moment**, so it carries
one year. A storm series, a rate fit or a model run spans an **interval**, so
it carries both ends. Two conventions, deliberately, because they mean
different things.

The end year is a **boundary, not a simulated year**: the model spends
`start..end-1`, so `1996_2010` drives a run of 1996 through 2009.

## 3. Where a file is named for a period, its vintage is recorded

Files are resolved by the **period start year**, so a survey from a nearby year
is copied under the period's name rather than the loader being taught about
vintages. That is a deliberate trade, and it costs nothing only if the real
vintage is written down beside it. Every such folder carries a `PROVENANCE.md`
saying what the survey actually is and what follows from the difference.

## 4. Retirement is a dated folder with a note

```
superseded_20260907/
    WHY.md        what this held, and why it stopped being used
```

Not `old_`, not `_ARCHIVE_`, not `_backup`, not a bare `old/`. A date says when
the decision was taken; the note says what the decision was. **Nothing is
deleted for being superseded** — a retired script is often the only record of
how something was done.

## 5. Find the root by searching upward, never by counting

```python
REPO = next(p for p in Path(__file__).resolve().parents
            if (p / "pyproject.toml").exists())
```

A counted depth (`parents[2]`) is correct only while the file stays where it
was written, and a wrong count does not raise — it resolves to a real
directory one level off and reads the wrong thing. Fourteen scripts were
converted before the `hatteras_ms` reorganisation on 2026-09-13 precisely so
the move could not break them silently.

Never write an absolute path into a home directory, and never a path rooted at
the drive.

**There is a residue of about 23 files the check will keep reporting**, and it
is not work left undone. Each names something that no longer exists and cannot
be mapped to anything that does:

* runs that were deleted (`HAT_1978_1997_natural`, `HAT_2004_2024_base_newbufferv3`);
* trees deleted before this repository's current shape (`data/hatteras_init/topography/2009_FIXED/`);
* another person's machine entirely, in the Ocracoke and Pea Island work.

Guessing a replacement would be worse than leaving the literal visible, because
a plausible wrong path fails silently and an obviously broken one does not. One
of the 23 is a comment describing this very bug, which the check cannot tell
from the bug itself.

**2026-09-22: the residue was triaged, and "cannot be mapped" was too strong.**
It is down to 19. Sorting the literals by what they actually name:

* **another machine** — the Ocracoke work points into
  `C:\Users\frank\OneDrive - University of North Carolina...`. Nothing to do.
* **a target that is genuinely gone** — the two deleted runs and
  `topography/2009_FIXED/`. The paragraph above holds for these: leave the
  literal visible.
* **a target that MOVED and kept its name** — one file was in this class.
  `figure_making/shoreline/dsas/dsas_from_gis.py` named
  `data/hatteras_init/shoreline_change/dsas_1997_2019_{rates,domain_means}.csv`;
  both files exist today under `5-scr/1-observations/dsas_1978_2019/` with the
  **same filenames**. That is a resolved relocation, not a guess, so it was
  fixed — through `hat_observed_rates.DSAS_ROOT`, per rule 6. The script had
  been dead, raising FileNotFoundError on import, and now runs.

The lesson for the next sweep: check whether the basename still exists
somewhere before filing a literal under "unmappable". An exact filename match
after a folder move is evidence, not a guess. A *similar* name is still a
guess, and the paragraph above still applies to it.

## 6. A location has one owner

Where something lives is decided once, in a resolver, and everything else asks:

| Resolver | Owns |
|---|---|
| `scripts/site_layer/hat_topo_version.py` | which Barrier3D domains a period reads |
| `scripts/site_layer/hat_observed_rates.py` | the CoastSat record and the rate fits |
| `scripts/site_layer/hat_elevation_products.py` | the DEM products |
| `scripts/site_layer/hatteras_site_config.py` | the period table and the forcing it names |

Twenty-two files built the rate path by hand, which is why moving it once cost
a day. A resolver makes the next move one edit, and turns a name that does not
exist into a loud error instead of a silently stale read.

## 7. Every folder a person opens has a `README.md`

One short file saying what is here, what produced it, and what not to trust.
The four data folders with the most complicated histories had one; the rest did
not, which is backwards — the gap was widest where a newcomer starts.

**Orientation is inherited.** A folder whose parent carries a README is
explained there, so the rule asks at the frontier: a folder holding files whose
parent explains nothing. Otherwise documenting a tree never finishes — one
sweep directory alone holds 488 machine-named cells, and nobody reads a README
for any of them. The corollary is that a parent README has to actually name its
children, because the check can see that one exists but not that it is any
good.

---

## Where things are

Two indexes sit beside this file and answer the questions this one does not:

- [`FIGURES.md`](FIGURES.md) — which figure answers which question, across
  `3-rates/`, `4-comparisons/` and `output/comparisons/`. Generated by
  `scripts/repo_tools/hat_write_figure_index.py`, which checks every path it
  prints, so it cannot rot silently.
- [`data/hatteras_init/5-scr/WINDOWS.md`](data/hatteras_init/5-scr/WINDOWS.md)
  — which window is which. Five window folders sit as peers across two chains
  and one of them is context only; the folder name says none of that.
  Generated from `hat_observed_rates.WINDOW_ROLE`.


```
cascade/          the model package
scripts/          all code
    hat_*.py          the resolvers and the site config
    input_prep/       builds model inputs, one folder per stage
    hatteras_ms/      the hindcast: the run, tools, experiments
    figure_making/    figures drawn from finished runs
    analyze_output/   cross-run analysis
    cascade_pipeline/ the shared library the runner imports
data/hatteras_init/   every model input, numbered by stage 0-9
output/               everything the model produced
tests/                the suite
hard-structures/      the groin study, which keeps its own conventions
```
