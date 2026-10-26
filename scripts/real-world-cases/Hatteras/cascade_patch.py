"""
Fixes to obx-cascade's cascade/ package that the Hatteras runs need, applied at run time.

    import cascade_patch
    cascade_patch.apply()

The clone holds scripts only, so a fix to the shared model code cannot go into
cascade/ here until the lab takes it. Until then it is patched in from this
file when the runner starts, the way barrier3d-patch/ carries the Barrier3D
fix. PATCHES names each one; the runner records the list in every run's
metadata, so a run says which code it actually used.

1. beach_dune_manager.resize_interior_domain -- the crash fix of 2026-09-07
   (old branch hannahaline/hatteras-cascade, commit 9b89dfee;
   HATTERAS_CASCADE_CHANGES.md section 3). obx-cascade still has the version
   that returns mismatched grids, and filter_overwash then dies on a numpy
   broadcast error (95 vs 97 rows, first obx-cascade run 2026-10-23). The
   function below is the old branch's, copied unchanged. The same change as a
   diff for the lab is cascade-patch/resize_interior_domain_fix.patch; delete
   it here once the lab merges it.

An idea for later: cap only the sand the bulldozer adds to a dune, not the
whole dune cell.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-23
"""
from __future__ import annotations

import numpy as np

import cascade.beach_dune_manager as _bdm
from cascade.beach_dune_manager import CascadeError

# What apply() replaces, for the run metadata
PATCHES = ("beach_dune_manager.resize_interior_domain (bay-edge anchor, 9b89dfee)",)


def resize_interior_domain(
    pre_storm_interior, post_storm_interior, bay_depth, dune_migration
):
    """
    Align the pre- and post-storm interior domains onto the same cross-shore rows
    so the two can be differenced.

    Barrier3D indexes the interior seaward-to-landward: row 0 sits behind the
    dunes, the last row is the bay. Between the pre-storm snapshot and the
    post-storm domain the grid can change length at either end, and WHICH end
    tells you how to line the two up:

      * rows are consumed at the SEAWARD end as the shoreline retreats into the
        island. The bay edge does not move, so it is the anchor: keep the LAST
        `n_post` rows of the pre-storm domain.
      * rows are appended at the BAY end as overwash rolls the island landward.
        The seaward edge is then the anchor, so the pre-storm domain is padded
        at the back with bay cells.

    This replaces a version that tried to reconstruct the seaward shift from
    `dune_migration` (`barrier3d.ShorelineChangeTS`) and then mop up any
    remainder by trimming trailing all-bay rows. That decomposition was wrong in
    both halves. The shoreline change is not the number of interior rows the
    front lost -- measured over one Hs = 1.2 m Hatteras run, `dune_migration`
    read -1.0 at three steps where the front had actually given up 2, 3 and 2
    rows -- and the trailing-row mop-up required the excess to be ENTIRELY at or
    below bay depth, which a partly-filled back-barrier row is not. When neither
    branch closed the gap the function returned mismatched arrays with nothing
    checking them, and `filter_overwash` died 140 lines later on a numpy
    broadcast error that named neither this function nor the cause.

    Anchoring on the bay edge is not a new rule, it is the rule the old code was
    approximating: over a run that completes under the old version, all 45
    `pre > post` calls produce a byte-identical array under both, so this cannot
    move a result that already exists. It only resolves the cases that used to
    raise.

    KNOWN LIMIT: a step that loses seaward rows AND gains bay rows at the same
    time has no single anchor, and this will mis-align it by the number of rows
    gained. Not seen in any Hatteras run measured so far (48 of 48 length
    changes were one-ended), and it cannot be detected from the two grids alone
    -- it needs Barrier3D to report the two counts separately.

    Parameters
    ----------
    pre_storm_interior: grid
        Pre-storm Barrier3D interior domain [dam * dam * dam]
    post_storm_interior: grid
        Post-storm Barrier3D interior domain [dam * dam * dam]
    bay_depth: float
        Bay depth [dam]
    dune_migration: float
        The number of grid cells [1 dam each] that the dunes migrated. Retained
        for call compatibility and no longer read: it is the quantity whose
        disagreement with the actual row change caused the failure above.

    Returns
    -------
    grid
        pre_storm_interior: resized to the post-storm row count
        post_storm_interior: unchanged
    """

    pre_storm_interior = np.asarray(pre_storm_interior)
    post_storm_interior = np.asarray(post_storm_interior)

    # The alongshore length is fixed for a Barrier3D domain, and nothing here
    # could align a change in it. Checked so that a violation is reported as
    # itself rather than as a broadcast error inside filter_overwash.
    if np.size(pre_storm_interior, 1) != np.size(post_storm_interior, 1):
        raise CascadeError(
            "pre- and post-storm interior domains disagree on alongshore "
            f"length ({np.size(pre_storm_interior, 1)} vs "
            f"{np.size(post_storm_interior, 1)} cells); "
            "resize_interior_domain only reconciles the cross-shore axis"
        )

    rows_pre = np.size(pre_storm_interior, 0)
    rows_post = np.size(post_storm_interior, 0)

    # seaward rows lost: anchor on the bay edge
    if rows_pre > rows_post:
        pre_storm_interior = pre_storm_interior[rows_pre - rows_post :]

    # bay rows gained: anchor on the seaward edge, pad the pre-storm bay
    elif rows_post > rows_pre:
        pre_storm_interior = np.concatenate(
            (
                pre_storm_interior,
                (
                    np.zeros([rows_post - rows_pre, np.size(post_storm_interior, 1)])
                    - bay_depth
                ),
            ),
            axis=0,
        )

    return pre_storm_interior, post_storm_interior




def apply():
    """Swaps the fixed functions into the cascade modules; returns PATCHES.

    BeachDuneManager.update looks resize_interior_domain up in its module at
    call time, so replacing the module attribute is enough. Safe to call twice.
    """
    _bdm.resize_interior_domain = resize_interior_domain
    return PATCHES
