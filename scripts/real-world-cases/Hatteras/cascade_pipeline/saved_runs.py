"""
Open a run's saved model state (the run folder's {run}.npz), including runs saved before the groin
became Cascade's GroinManager module.

    from cascade_pipeline.saved_runs import load_saved_cascade
    cascade = load_saved_cascade(run_dir / f"{run_dir.name}.npz")

Cascade.save pickles the whole model, and pickle stores each object by its module and class name. Runs
saved before 2026-10-23 with the groin on hold the groin under a name cascade.groin no longer has
(GroinCallback; CorrectedGroinCallback before 2026-10-20), attached as cascade._groin_callback, and
most were saved by the hatteras-cascade sandbox Cascade (cascade.cascade_groin), which this repo does
not have; a plain np.load fails on both. This module is the one place those old names are kept, only
to read old runs. While loading it reads cascade.cascade_groin as cascade.cascade and answers any name
cascade.groin lacks with a plain stand-in class; it then turns the old wave-driven groin into a
GroinManager and attaches it as cascade.groin, so its record, strength history, summary() and
diagnostics_frame() read as for a new run. The loaded groin is for reading: an old run is not stepped
on. Older groin kinds (the dipole and blocking groins) stay as their stand-in, their saved attributes
readable. Runs without a groin load unchanged.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-23
"""
from __future__ import annotations

import importlib.util
import sys

import numpy as np

import cascade.cascade as cascade_module
import cascade.groin as groin_module
from cascade.groin import GroinManager

SANDBOX_MODULE = "cascade.cascade_groin"   # the hatteras-cascade branch's Cascade, which saved most old runs

# The old groin's attribute names (the 2026-10-20 groin's, then the sandbox's), and the GroinManager's
_RENAMED = {
    "updrift_pad": "_updrift_pad",
    "downdrift_pad": "_downdrift_pad",
    "strength_history": "_strength_history",
    "install_year": "_install_year",
    "start_year": "_start_year",
    "n_domains": "_n_domains",
    "half_width": "_half_width",
    "spread_weights": "_spread_weights",
    "record": "_record",
    "_lo": "_lo_pad",
    "_hi": "_hi_pad",
    "_coast_diff0": "_coast_diff_at_start",
    "_q_net": "_drift_table_m3_yr",
}


class _SavedObject:
    """Stand-in for a class an old run pickled from cascade.groin; holds the saved attributes."""


def _stand_in(name):
    """A plain class under ``name``, for pickle to rebuild an old object into."""
    return type(name, (_SavedObject,), {"__module__": groin_module.__name__})


def _as_groin_manager(old):
    """The old wave-driven groin as a GroinManager (for reading), or None if it is an older kind."""
    state = dict(vars(old))
    if "strength_history" not in state or "record" not in state:
        return None
    for old_name, new_name in _RENAMED.items():
        if old_name in state:
            state[new_name] = state.pop(old_name)
    state.setdefault("_lo_pad", min(state["_updrift_pad"], state["_downdrift_pad"]))
    state.setdefault("_hi_pad", max(state["_updrift_pad"], state["_downdrift_pad"]))
    state.setdefault("_spread_weights", GroinManager.NO_SPREAD)
    away = 1 if state["_updrift_pad"] == state["_hi_pad"] else -1
    state.setdefault("_spread_pads", tuple(state["_updrift_pad"] + away * step for step in range(3)))
    state.setdefault("_footprint_m", None)
    # the spread came in 2026-10-21; years recorded before it had none
    state["_record"] = [{"spread_loss_m": 0.0, **row} for row in state["_record"]]
    manager = GroinManager.__new__(GroinManager)
    manager.__dict__.update(state)
    return manager


def _attach_groin(cascade):
    """Give an old Cascade the groin attributes the current Cascade has."""
    if "_groin" in vars(cascade):
        return cascade
    # the attribute the old runner attached the groin through, before Cascade built its GroinManager
    old = vars(cascade).pop("_groin_callback", None)
    groin = old
    if isinstance(old, _SavedObject):
        groin = _as_groin_manager(old) or old   # an older kind stays as its stand-in
    cascade._groin = groin
    cascade._groin_module = groin is not None
    return cascade


def load_saved_cascade(path):
    """Load the Cascade saved in a run's .npz, whichever groin version and Cascade saved it.

    Args:
        path: The run's {run}.npz (Cascade.save output).

    Returns:
        The saved Cascade, with cascade.groin set (None for a run without a groin).
    """
    had_getattr = "__getattr__" in vars(groin_module)

    def missing_name(name):
        if name.startswith("__"):
            raise AttributeError(name)
        stand_in = _stand_in(name)
        setattr(groin_module, name, stand_in)
        return stand_in

    map_sandbox = SANDBOX_MODULE not in sys.modules and importlib.util.find_spec(SANDBOX_MODULE) is None
    if map_sandbox:
        sys.modules[SANDBOX_MODULE] = cascade_module
    if not had_getattr:
        groin_module.__getattr__ = missing_name
    names_before = set(vars(groin_module))
    try:
        # the object array is unpickled when it is indexed, so index inside the try
        with np.load(path, allow_pickle=True) as saved:
            cascade = saved["cascade"].flat[0]
    finally:
        if not had_getattr:
            del groin_module.__getattr__
        for name in set(vars(groin_module)) - names_before:
            delattr(groin_module, name)
        if map_sandbox:
            del sys.modules[SANDBOX_MODULE]
    return _attach_groin(cascade)
