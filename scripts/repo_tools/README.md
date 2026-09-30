# `repo_tools/` — tools that act on the repository itself

Not site content and not part of any analysis: these read the tree and report
on it. They are **run**, never imported.

```
hat_layout_check.py          audit the tree against the seven rules in ORGANIZATION.md
hat_write_figure_index.py    write FIGURES.md: which shoreline-change figure to open for a question
style_guide_check.py         report where scripts depart from scripts/STYLE.md
style_equivalence_check.py   prove a restyled script does exactly what it did (parse-tree check)
style_run_compare.py         run scripts before and after a restyle in a scratch copy, compare outputs
style_readme_coverage.py     check that explanation removed from a script reached its README
style_restyle.py             the reassembly helper the restyle used
template/                    a general project layout and a script that builds it empty, to share
```

```
python scripts/repo_tools/hat_layout_check.py            every rule
python scripts/repo_tools/hat_layout_check.py --rule 5   just one
python scripts/repo_tools/hat_layout_check.py --full     every offender
```

**It is advisory and always exits zero.** Chosen deliberately (Hannah,
2026-09-13): a check that blocks you mid-experiment gets disabled, and a
disabled check reports nothing. This one is meant to be run when you want a
picture, and it produces a worklist rather than an obstacle.

It finds the repo root by searching upward for `pyproject.toml`, so it does not
care where under the tree it is invoked from or where it is moved to.

The `style_*` tools were written for the 2026-09-30 sweep that brought every
script in line with `scripts/STYLE.md`; each one's header says how to run it.

## The scripts in detail

For the two older tools, below is the original header and any notes that were
in the code, kept word for word when they were brought in line with
`scripts/STYLE.md` (2026-09-30).

### hat_layout_check.py

Check the repository against the seven layout rules in ORGANIZATION.md.

Notes that were in the code:

```text
hat_layout_check.py

Does the repository still follow the seven rules in ORGANIZATION.md?

ADVISORY, ALWAYS EXITS ZERO. It produces a worklist, not an obstacle. Every
tidy-up in this project so far has been someone noticing a case by hand,
which does not scale and decays between surveys. This is the same survey,
run on demand.

python scripts/repo_tools/hat_layout_check.py            every rule
python scripts/repo_tools/hat_layout_check.py --rule 5   just one
python scripts/repo_tools/hat_layout_check.py --full     every offender, not the
first few

WHY IT DOES NOT FAIL THE BUILD
Chosen deliberately (Hannah, 2026-09-13): a check that blocks you mid
experiment gets disabled, and a disabled check reports nothing. This one
is meant to be run when you want a picture.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-22
```

```text
Trees that hold the project's own work. Vendored code and caches are not
judged by these rules.
```

```text
.log added 2026-09-22. A 100 kB run log sat in
scripts/hatteras_ms/ for over a year and nothing reported
it: a log is a PRODUCT, and products live in output/
(output/README.md says where). Being gitignored hides it
from a diff, which is exactly why the check has to see it.
```

```text
A few data-shaped files belong beside code because they ARE code's input in
the sense of configuration, or documentation of it.
```

```text
Rule 1, second case: the root of scripts/ holds no files but README.md.

It held eleven on 2026-09-18 -- six site modules, this checker, and four
hatteras_site_config_prebe_<stamp>.py snapshots of the solved BE field that
the calibrate step had written beside the config it copied. The snapshots are
data wearing a .py extension, which is why the suffix test above never saw
them; they read as stray copies, which is how the equivalent 2026-08-24 one
came to be discarded, taking the only record of the one-shot solve with it
and leaving plot_be_zones.py undrawable for three weeks.

The modules then moved into site_layer/ and the checker into repo_tools/, so
the root is now a list of folders and one README. That is worth holding: a
reader who opens scripts/ should see the map, not the map plus whatever was
most recently left lying on it. Anything new at this level belongs in one of
the folders, or is a snapshot that belongs in the data tree.
```

```text
Rule 4: the one retirement idiom.
A retirement folder is `superseded_<date>`, and MAY carry a reason after the
date: superseded_20260919_pre-redigitized. Relaxed 2026-09-22 -- the bare form
was the rule until then, and it could not express the case the data tree
actually has, which is two retirements in one folder. A date alone cannot tell
`1996/superseded_20260915_flat` from the next one; the suffix is what
distinguishes them, so demanding a bare date asked for information to be
thrown away. The DATE STILL COMES FIRST, so the folders sort chronologically
and the rule's point -- a date says when the decision was taken -- survives.
```

```text
`output/archive/` files retired material as `YYYY-MM-DD_<what>/`, its own
documented idiom (output/README.md). What sits inside one of those is filed,
not stray, so it is not asked to be a superseded_ folder as well.
```

```text
Rule 5: a counted depth used to find a root. Only flagged when the name being
assigned looks like a root, so ordinary parents[] use is not noise.
```

```text
Retired code is not maintained -- every superseded folder says so -- and a
rule 5 complaint about a script nobody will run again is noise that makes the
real ones harder to see. Rules 4 and 7 still apply to these folders; rule 5
does not (2026-09-14).
```

```text
A file with no extension has no suffix to match, so a suffix check cannot
see it at all. Two sat in the trees untouched -- `HAT_hindcast_plan` (a
planning note the hatteras_ms README described as a FOLDER) and `Notes` (the
only record of the storm max-duration test). Both are now typed. These names
are the extensionless files that are meant to be here.
```

```text
Untyped: it renders nowhere, no tool can classify it, and it is
invisible to every check that works by suffix. Give it one.
```

```text
A notebook is JSON holding the same statements, one per string,
so the unanchored pattern reads both without a JSON parse.
```

<details><summary>Function notes (the original docstrings)</summary>

**`imported_module_names()`**

```text
Every bare name any code in the repo imports.

Read once, not once per candidate. Only the code trees are scanned: a .py
under data/ is itself data -- the BE rate tables are spelled that way --
and treating one as an importer would let a snapshot vouch for a snapshot.
```

**`rule_3_provenance_beside_derived()`**

```text
Period folders holding a file named for a year, with no PROVENANCE.md.

Only checked where the folder name IS a bare year, which is where the
period-start naming convention applies.
```

**`inside_compliant_retirement()`**

```text
Is this already filed under a dated superseded folder?

What sits INSIDE one keeps its original name on purpose -- that is what it
was called when it was in use, and renaming it would erase that. So the
convention applies to the retirement folder, not to its contents.
```

**`rule_7_readmes()`**

```text
Folders that need a README of their own.

ORIENTATION IS INHERITED (2026-09-14). A folder whose PARENT carries a
README is already explained there, so asking for one in every child turns
the rule into noise -- one sweep directory alone holds 488 machine-named
cells, none of which anyone reads a README for. The rule asks at the
FRONTIER: a folder holding files whose parent explains nothing.

That makes documenting a tree from the top genuinely finish, rather than
exposing a new row of demands each time.
```

</details>

### hat_write_figure_index.py

Write FIGURES.md at the repo root: which figure to open for each shoreline-change question.

From the script's original header:

```text
Write FIGURES.md at the repo root: one page that answers "which figure do I
open for X" across the three trees that hold shoreline-change work --
`data/hatteras_init/5-scr/3-rates/`, `.../4-comparisons/` and
`output/comparisons/`.

WHY THIS EXISTS (Hannah, 2026-09-21).  Each of those trees has a good README
of its own and none of them spans the others, so answering a question meant
already knowing which tree it lived in. Three things in particular were easy
to get wrong and are stated here once:

  * the estimator            rate (m/yr) vs distance (m), and which window
                             the rate was FITTED on -- see the vocabulary in
                             3-rates/README.md
  * the window               five window folders sit as peers across two
                             chains; see WINDOWS.md
  * the units                model_vs_observed/ is in m/yr; every other
                             comparison tree is in METRES. Two figures of the
                             same comparison in different units is the single
                             most confusable pair in the project.

HOW IT STAYS HONEST.  The map below is editorial -- a person decides which
question a figure answers -- but every path in it is CHECKED against the disk
when this runs, and a missing one is reported and marked in the output rather
than silently written. So the index cannot quietly rot the way a hand-kept
list does. Run it after adding or renaming a figure.

USAGE
    python scripts/repo_tools/hat_write_figure_index.py
    python scripts/repo_tools/hat_write_figure_index.py --check   # no write
```

Notes that were in the code:

```text
question -> [(what you get, path, note)]. A path ending in "/" is a folder.
<w> is substituted per window where a row covers several.
```

```text
Item 4 of the 2026-09-21 tidy: model_vs_observed slices by OBSERVATION while
every other tree slices by ESTIMATOR. Rather than rename it -- its axis
genuinely is a different question -- say which estimator each folder uses.
```

