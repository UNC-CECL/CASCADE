# groin_init - the groin study's input products

The products of `../input_prep/` (island offset, storms, target) and one repair
tool. Map: `../README.md`.

```
fix_cascade_yaml.py    repair a CASCADE parameter YAML that numpy scalars were written into
island_offset/         1967 D2-D12 offsets: input geojsons, raw offsets, unpadded and padded CSVs
storms/                storm files: input_storms/, 1967_1997/, 1967_2017/, 1967_2024/
target/                HAT_target_1967_1997_*: the observed dune-line target
```

## The scripts in detail

Each script's header says what it does and how to run it; the reasoning and choices behind it are here. Moved out of the scripts on 2026-10-01, when they were brought in line with `scripts/STYLE.md` (code unchanged, proven with `style_equivalence_check.py`; path fixes listed per script).

### fix_cascade_yaml.py

A one-shot repair, not part of any pipeline. When a numpy value (`np.float64(...)`)
is serialised into the parameter YAML, CASCADE's `yaml.full_load` refuses to
rebuild it and every run crashes on load. The script backs the file up to
`.bak`, loads it with `yaml.UnsafeLoader` (which can rebuild the numpy
objects), converts every numpy scalar and array back to plain Python, rewrites
it with `yaml.safe_dump` (key order kept; CASCADE reads by key, so order and
comments do not matter to it) and proves `full_load` now succeeds.

`PARAM_FILE` points at `data/hatteras_init/Hatteras-CASCADE-parameters.yaml`,
anchored on the repo root since 2026-09-14.

From the script's original header:

```text
fix_cascade_yaml.py
===================
One-shot repair for a CASCADE parameter YAML that has been "poisoned" with
NumPy scalar objects, e.g.:

    yaml.constructor.ConstructorError: could not determine a constructor for the
    tag 'tag:yaml.org,2002:python/object/apply:numpy._core.multiarray.scalar'

That happens when a numpy value (np.float64(...), etc.) got serialized into the
file instead of a plain number. CASCADE reads the file with yaml.full_load,
which refuses to reconstruct those objects, so every run crashes on load.

This script loads the file with a loader that CAN rebuild the numpy objects,
converts every numpy scalar/array back to a plain Python float/int/list, and
rewrites the file with plain numbers only. A .bak copy is made first.

USAGE
-----
Just run it. Edit PARAM_FILE below if your path differs.
    python fix_cascade_yaml.py
```

<details><summary>Function notes (the original docstrings)</summary>

**`to_plain()`**

```text
Recursively convert numpy scalars/arrays to plain Python types.
```

</details>

<details><summary>Notes that were comments in the code</summary>

- Anchored 2026-09-14: this named a home directory, or a tree renamed since. Rule 5 of ORGANIZATION.md.
- 2. Load with UnsafeLoader -- unlike full_load, it can reconstruct the numpy scalar objects (numpy must be importable, which it is here).
- 4. Rewrite with safe_dump so only plain scalars are emitted. Key order is preserved; CASCADE reads by key, so order/comments don't matter to it.

</details>
