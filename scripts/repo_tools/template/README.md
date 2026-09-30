# template — how to lay out a research project, to share

A folder layout for a modelling or data-analysis project, and a script that
builds it empty. Meant to be **copied out of this repository** and handed to a
colleague starting a project. It is the general form of this repository's
`ORGANIZATION.md`, without anything specific to Hatteras.

```
project_layout_template.py   builds the skeleton below, with a README in
                             every folder; never overwrites a file
```

```
python project_layout_template.py my_project --stages observations model-inputs runs
```

Standard library only.

## The skeleton

```
my_project/
    README.md              what the project is, where to start
    pyproject.toml         marks the root (see rule 6)
    scripts/               all code
        lib/paths.py       the one place locations are decided
        1-observations/    code for each stage, numbered in run order
        2-model-inputs/
        3-runs/
    data/                  everything the code reads
        1-observations/    the same numbers as scripts/
        2-model-inputs/
        3-runs/
    output/                everything the code produces
        figures/           PNGs; PDFs and CAPTIONS.md under supporting/
        logs/
        archive/           YYYY-MM-DD_<what>/
```

## The rules

**1. Code, data and products are three trees.** Scripts live in `scripts/`,
what they read in `data/`, what they make in `output/`. The test: could you
hand someone `data/` alone and have them see everything the analysis uses? A
log, a snapshot or an export is data or output even when it sits next to code,
and even when it has a `.py` extension.

**2. Stages are numbered, and the numbers match across trees.** `scripts/2-model-inputs/`
makes `data/2-model-inputs/`. The number is the run order, so nothing in stage
3 can depend on stage 4. Folders carry numbers; file names don't.

**3. A moment is named for its year, an interval for its span.** A survey is
`1996/`; a rate fit or a model run is `1996_2010/`. Write down whether the end
year is included.

**4. Versions restart with new source data.** Fixes to the same source continue
the count (`v1`, `v2`, ...). New source data (a re-digitized line, a new
survey) starts again at `v1`, and the old builds move to a retirement folder.
A `CURRENT` file in the folder names the version in use.

**5. Retire, with a date and a reason.** Old work moves to
`superseded_YYYYMMDD[_reason]/` with a `WHY.md` saying what it was and why it
stopped being used; output goes to `output/archive/YYYY-MM-DD_<what>/`. Not
`old/`, `backup/` or `_final2`. The date comes first so folders sort in order.
If you delete instead, say in the README how to get it back from git.

**6. Find the root by searching upward; never type a path.** Every script
finds the project root by walking up to `pyproject.toml`, and gets locations
from `scripts/lib/paths.py`:

```python
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
```

A counted depth (`parents[2]`) breaks silently when a file moves; a path typed
into a home directory breaks on every other machine.

**7. A location has one owner.** If three scripts build the same path by hand,
moving that folder breaks all three. Decide each location once, in
`paths.py`, and make an unknown name raise an error that lists the real ones.

**8. Every folder a person opens has a README.** Three short answers: what is
here, what made it, and what not to trust. A folder whose parent README
already explains it can skip its own.

**9. Products say where they came from.** A data or output folder that a
script writes carries a `PROVENANCE.md`: the script, the date, the inputs,
and anything assumed.

## Companion templates

In the same repository:

- `scripts/STYLE.md` — how a script is laid out, commented and signed.
- `scripts/input_prep/5-scr/template/` — the shoreline-rate scripts.
- `scripts/figure_making/template/` — a figure in the house style.

---

Hannah A. Henry, Coastal Environmental Change Lab, University of North Carolina at Chapel Hill  
hahenry@unc.edu · version 2026-09-30
