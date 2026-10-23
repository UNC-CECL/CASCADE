# CASCADE patch for the Hatteras hindcast

The hindcast runs on this repository's `cascade/` package (`obx-cascade` at
`e6a6642a`) with one fix that the lab has not merged yet. The clone holds
scripts only, so the fix is not edited into `cascade/` here. It reaches the
runs in two ways:

- **At run time:** [`../cascade_patch.py`](../cascade_patch.py) swaps the fixed
  function in when the runner starts, and every run records the list in its
  metadata (`cascade_patches`).
- **For the lab:** `resize_interior_domain_fix.patch` is the same change as a
  diff against `cascade/beach_dune_manager.py`, ready to offer as a bug fix.
  Not yet proposed (2026-10-23).

`resize_interior_domain_fix.patch` changes `cascade/beach_dune_manager.py` only:

| | where | status |
|---|---|---|
| 1 | `resize_interior_domain` | **patched at run time.** When the pre- and post-storm interior grids differ in cross-shore length, the obx version rebuilds the seaward shift from `dune_migration` and trims trailing all-bay rows. Neither step closes the gap reliably, so it can return mismatched grids and `filter_overwash` dies on a numpy broadcast error. The first Hatteras run on obx-cascade hit this (95 vs 97 rows, 2026-10-23). The fix anchors on the edge that did not move: seaward rows lost, keep the last rows of the pre-storm grid; bay rows gained, pad the pre-storm grid with bay cells |

An idea for later: cap only the sand the bulldozer adds to a dune, not the
whole dune cell.

The fix comes from the old branch `hannahaline/hatteras-cascade` (commit
`9b89dfee`, 2026-09-07), where the full record with evidence is
`HATTERAS_CASCADE_CHANGES.md` section 3. The Barrier3D side is
[`../barrier3d-patch/`](../barrier3d-patch/README.md).

## Offering it to the lab

From the repository root:

```
git apply --check scripts/real-world-cases/Hatteras/cascade-patch/resize_interior_domain_fix.patch
```

It was checked against `e6a6642a` on 2026-10-23. Once the lab merges it,
delete the function and its `PATCHES` entry from `../cascade_patch.py`.
