"""Roadway forcing for a CASCADE run: setbacks, elevations, events, drowning.

Everything CASCADE's `roadway_manager` needs, prepared before the run and
checked against the interior it will actually be spent on. Nothing here is
site-specific -- domain geometry arrives as a `DomainGeometry`, and the
Hatteras instances (community zones, historical events, file names) live in
`hatteras_site_config`.

Three Barrier3D conventions this module depends on:

- `road_setback` and `road_width` are METERS; `bulldoze` divides them by the
  cell size to get row indices, so a setback is only ever as precise as one
  cell.
- `road_ele` is METERS, MHW-RELATIVE. `bulldoze` writes it straight into the
  interior grid, which the extractor stores MHW-relative -- see the datum note
  on `load_road_elevations`.
- `bulldoze` tests the rows FLANKING the road, never the road's own cells, and
  drowns the road when either flank is more than `percent_water` water.
  `predict_drowning` reproduces that test rather than approximating it.
"""

import dataclasses

import numpy as np


@dataclasses.dataclass(frozen=True)
class RoadwayConfig:
    """Unit and threshold settings matching CASCADE's roadway_manager.

    Attributes:
        road_width_m: Roadway width in meters, as passed to Cascade.
        cell_size_m: Grid cell size in meters (Barrier3D's dx/dy).
        dam_to_m: Decameters-to-meters factor for the stored arrays.
        drown_threshold_m: Elevation at or below which a cell counts as water.
            roadway_manager passes 0 and calls it "m MSL", but compares it
            against MHW-relative elevations, so it is effectively 0 m MHW.
        percent_water: Fraction of flanking cells that may be water before the
            road drowns (roadway_manager's percent_water_cells_touching_road).
        sea_level_dam: Barrier3D's SL. Fixed at 0 in the Lagrangian frame --
            sea-level rise lowers the domain rather than raising SL.
        elevation_fallback_m: road_ele for domains with no measured value.
    """

    road_width_m: float = 20.0
    cell_size_m: float = 10.0
    dam_to_m: float = 10.0
    drown_threshold_m: float = 0.0
    percent_water: float = 0.2
    sea_level_dam: float = 0.0
    elevation_fallback_m: float = 1.45


DEFAULT_ROADWAY = RoadwayConfig()


@dataclasses.dataclass(frozen=True)
class RelocationEvent:
    """A historical roadway relocation, as a DISPLACEMENT.

    The road is moved landward by a measured distance; the resulting setback is
    that displacement added to whatever setback the model is carrying at the
    time. It is deliberately not an absolute setback: an absolute value
    referenced to an older dune line double-counts the shoreline retreat
    between that line and the topography's own dune line, which is what put
    NC-12 into the sound in earlier versions of this pipeline.

    Attributes:
        year: Calendar year the relocation happened.
        displacement_m: Landward displacement per GIS domain, in meters.
        note: Human-readable description for run metadata.
        enabled: Whether to apply this event.
    """

    year: int
    displacement_m: dict
    note: str = ""
    enabled: bool = True


@dataclasses.dataclass(frozen=True)
class BridgeEvent:
    """A bridge or alternate route replacing the road surface.

    Roadway management stops for the listed domains from this year on, with no
    setback change -- the road is gone, not moved.

    Attributes:
        year: Calendar year the bridge opened.
        gis_domains: GIS domains whose road is removed.
        note: Human-readable description for run metadata.
        enabled: Whether to apply this event.
    """

    year: int
    gis_domains: tuple
    note: str = ""
    enabled: bool = True


# =============================================================================
# Loading
# =============================================================================

def load_padded_series(path, geometry, first_gis, last_gis, fill=0.0):
    """Reads a 2-row CASCADE forcing file into a padded per-domain array.

    The 2-row format is (GIS IDs, values). CASCADE consumes forcings as arrays
    indexed by padded position, so the IDs are used to place values rather than
    being discarded -- reading positionally would silently shift every domain
    after a gap.

    Args:
        path: Path to the 2-row CSV.
        geometry: DomainGeometry describing the padded array.
        first_gis: First GIS domain the forcing covers.
        last_gis: Last GIS domain the forcing covers.
        fill: Value for domains the file does not cover.

    Returns:
        A (values, missing) tuple. values is a float array of length
        geometry.total_domains; missing lists GIS domains in the expected span
        that the file did not supply.

    Raises:
        ValueError: If the file is not 2 rows.
    """
    raw = np.loadtxt(path, delimiter=",")
    if raw.ndim != 2 or raw.shape[0] != 2:
        raise ValueError(
            f"{path}: expected a 2-row (IDs, values) file, got shape "
            f"{raw.shape}")

    by_gis = dict(zip(raw[0].astype(int).tolist(), raw[1].tolist()))
    values = np.full(geometry.total_domains, float(fill))
    missing = []
    for gis in range(first_gis, last_gis + 1):
        pad = geometry.gis_to_pad(gis)
        if gis in by_gis and 0 <= pad < values.size:
            values[pad] = by_gis[gis]
        else:
            missing.append(gis)
    return values, missing


def load_road_setbacks(path, geometry, first_gis, last_gis):
    """Loads road setbacks, in meters landward of the dune line.

    Args:
        path: Path to the 2-row RoadSetback_<year>.csv.
        geometry: DomainGeometry describing the padded array.
        first_gis: First GIS domain carrying road.
        last_gis: Last GIS domain carrying road.

    Returns:
        A (setbacks_m, missing) tuple; setbacks are 0.0 outside the road span.
    """
    return load_padded_series(path, geometry, first_gis, last_gis, fill=0.0)


def load_road_elevations(path, geometry, first_gis, last_gis,
                         config=DEFAULT_ROADWAY):
    """Loads per-domain road elevation, in meters MHW-relative.

    Datum: `bulldoze` writes road_ele into the interior grid after dividing by
    dz, and the extractor stores that grid MHW-relative. So road_ele must be
    MHW-relative meters -- NOT NAVD88. The file this reads is already in that
    frame; do not subtract MHW again.

    Not period-dependent. Road elevation is a property of the surveyed surface,
    and there is one topography for every period, so one elevation set serves
    all of them.

    Args:
        path: Path to the 2-row RoadElevation.csv. No year in the name --
            see the note above on why one set serves every period.
        geometry: DomainGeometry describing the padded array.
        first_gis: First GIS domain carrying road.
        last_gis: Last GIS domain carrying road.
        config: RoadwayConfig supplying the fallback elevation.

    Returns:
        A (elevations_m, missing) tuple. Domains the file omits, and every
        buffer domain, get config.elevation_fallback_m.
    """
    return load_padded_series(path, geometry, first_gis, last_gis,
                              fill=config.elevation_fallback_m)


# =============================================================================
# Management flags
# =============================================================================

def build_roadway_management_on(geometry, first_gis, last_gis,
                                community_zones=(), enabled=True):
    """Builds the per-domain roadway-management flag CASCADE consumes.

    Management runs on the road span minus the permanent community zones.
    Buffer domains are always off.

    Note CASCADE constructs a RoadwayManager for EVERY domain regardless of
    this flag ("always initialize just in case we want to add a road during the
    simulation"), so the existence of a manager says nothing about whether the
    road is managed -- only this flag does.

    Args:
        geometry: DomainGeometry describing the padded array.
        first_gis: First GIS domain carrying road.
        last_gis: Last GIS domain carrying road.
        community_zones: (first_gis, last_gis) pairs to exclude.
        enabled: Global switch; False turns management off everywhere.

    Returns:
        A list of bools, length geometry.total_domains.
    """
    flags = [False] * geometry.total_domains
    if not enabled:
        return flags

    for gis in range(first_gis, last_gis + 1):
        pad = geometry.gis_to_pad(gis)
        if 0 <= pad < len(flags):
            flags[pad] = True
    for zone_first, zone_last in community_zones:
        for gis in range(zone_first, zone_last + 1):
            pad = geometry.gis_to_pad(gis)
            if 0 <= pad < len(flags):
                flags[pad] = False
    return flags


# =============================================================================
# Drowning prediction
# =============================================================================

def interior_widths(interior_m, config=DEFAULT_ROADWAY):
    """Land run from row 0 to the first water cell, per alongshore profile.

    Transcribes barrier3d.FindWidths, including its `- 1` and its clamp at
    zero. This is the only width Barrier3D uses to decide where the island is;
    land behind a water cell is invisible to it.

    Args:
        interior_m: Interior elevations in meters MHW, (cross_shore, along).
        config: RoadwayConfig supplying the sea level.

    Returns:
        An int array of land-run lengths, one per alongshore profile.
    """
    n_rows, n_along = interior_m.shape
    sea_level_m = config.sea_level_dam * config.dam_to_m
    widths = np.empty(n_along, dtype=int)
    for along in range(n_along):
        column = interior_m[:, along]
        first_water = next(
            (row for row, value in enumerate(column) if value <= sea_level_m),
            n_rows)
        widths[along] = max(first_water - 1, 0)
    return widths


def predict_drowning(interior_m, setback_m, config=DEFAULT_ROADWAY):
    """Reproduces roadway_manager.bulldoze's drowning test at t=0.

    bulldoze checks the rows FLANKING the road -- `road_end + 1` on the bay
    side and `road_start - 1` on the sea side -- and never inspects the road's
    own cells. The row at `road_end` is skipped entirely: it supplies only the
    cell count. All three are reported so a caller can see the difference.

    Args:
        interior_m: Interior elevations in meters MHW, (cross_shore, along).
        setback_m: Road setback in meters.
        config: RoadwayConfig supplying widths and thresholds.

    Returns:
        A dict with road_start, road_end, border_row, n_rows, sea_water,
        bay_water, road_cells_water, drowns and wall. `wall` is None unless the
        setback would corrupt rather than drown the run: "NEGATIVE" (numpy
        wraps the index to the back of the barrier) or "OVERRUN" (the border
        row is past the array, which raises IndexError inside bulldoze).
    """
    n_rows, n_along = interior_m.shape
    road_start = int(setback_m / config.cell_size_m)
    road_end = road_start + int(config.road_width_m / config.cell_size_m)
    border = road_end + 1
    threshold = config.drown_threshold_m

    result = {
        "road_start": road_start, "road_end": road_end, "border_row": border,
        "n_rows": n_rows, "sea_water": np.nan, "bay_water": np.nan,
        "road_cells_water": np.nan, "drowns": False, "wall": None,
    }
    if road_start < 0:
        result["wall"] = "NEGATIVE"
        return result
    if border >= n_rows:
        result["wall"] = "OVERRUN"
        return result

    footprint = interior_m[road_start:road_end, :]
    result["road_cells_water"] = float((footprint <= threshold).mean())
    result["bay_water"] = float((interior_m[border, :] <= threshold).mean())
    result["sea_water"] = (
        float((interior_m[road_start - 1, :] <= threshold).mean())
        if road_start > 0 else 0.0)
    result["drowns"] = bool(result["sea_water"] > config.percent_water
                            or result["bay_water"] > config.percent_water)
    return result


def audit_setbacks(elevation_paths, setbacks_m, geometry, first_gis, last_gis,
                   management_on=None, config=DEFAULT_ROADWAY):
    """Predicts, before the run, which road_offset will not survive year one.

    A drowned road is not a warning: roadway_manager sets _drown_break and
    returns immediately on every later year, so the domain gets no overwash
    removal, no dune rebuilding and no relocation for the rest of the run. It
    becomes an unmanaged barrier wearing a road label, which is why these
    domains have to be named before any managed-vs-unmanaged comparison.

    Args:
        elevation_paths: Padded-order interior .npy paths, one per domain.
        setbacks_m: Padded per-domain setbacks, in meters.
        geometry: DomainGeometry describing the padded array.
        first_gis: First GIS domain carrying road.
        last_gis: Last GIS domain carrying road.
        management_on: Optional padded management flags; when given, the
            `managed` column records whether CASCADE will actually run the
            manager for that domain.
        config: RoadwayConfig supplying widths and thresholds.

    Returns:
        A list of per-domain dicts sorted by GIS id, each carrying gis, pad,
        setback_m, the predict_drowning fields, interior-width statistics and
        `managed`.
    """
    rows = []
    for gis in range(first_gis, last_gis + 1):
        pad = geometry.gis_to_pad(gis)
        if not 0 <= pad < len(elevation_paths):
            continue
        interior_m = np.load(elevation_paths[pad]) * config.dam_to_m
        prediction = predict_drowning(interior_m, float(setbacks_m[pad]),
                                      config)
        widths = interior_widths(interior_m, config)
        rows.append({
            "gis": gis, "pad": pad, "setback_m": float(setbacks_m[pad]),
            **prediction,
            "width_min": int(widths.min()),
            "width_median": int(np.median(widths)),
            "width_max": int(widths.max()),
            "managed": (True if management_on is None
                        else bool(management_on[pad])),
        })
    return rows


def summarise_audit(rows, config=DEFAULT_ROADWAY):
    """Reduces audit_setbacks output to the counts worth reporting.

    Args:
        rows: Output of audit_setbacks.
        config: RoadwayConfig supplying the water threshold, for the message.

    Returns:
        A dict with n_domains, n_managed, drowning (GIS ids that will drown AND
        are managed), drowning_unmanaged, and blocking (GIS ids whose setback
        would corrupt the run rather than drown it).
    """
    drowning = [r["gis"] for r in rows if r["drowns"] and r["managed"]]
    unmanaged = [r["gis"] for r in rows if r["drowns"] and not r["managed"]]
    blocking = [(r["gis"], r["wall"]) for r in rows if r["wall"]]
    return {
        "n_domains": len(rows),
        "n_managed": sum(1 for r in rows if r["managed"]),
        "drowning": drowning,
        "drowning_unmanaged": unmanaged,
        "blocking": blocking,
    }


# =============================================================================
# Historical events
# =============================================================================

def relocated_setbacks(setbacks_m, event, geometry):
    """Applies a relocation event to a padded setback array.

    Adds the measured displacement to the CURRENT setback. At t=0 the current
    setback is the initial one, so this is what the road's position would be if
    the event fired immediately -- the upper bound on the relocated setback,
    since CASCADE decrements the setback by dune migration between t=0 and the
    event year.

    Args:
        setbacks_m: Padded per-domain setbacks, in meters.
        event: A RelocationEvent.
        geometry: DomainGeometry describing the padded array.

    Returns:
        A new padded setback array; domains the event does not touch are
        unchanged.
    """
    out = np.asarray(setbacks_m, dtype=float).copy()
    for gis, displacement in event.displacement_m.items():
        pad = geometry.gis_to_pad(gis)
        if 0 <= pad < out.size:
            out[pad] += float(displacement)
    return out


def summarise_road_management(cascade, geometry, first_gis, last_gis):
    """Reports which domains kept their road, from CASCADE's post-run state.

    Reads state CASCADE already exposes rather than instrumenting the model:
    roadway_management_module says which domains were managed at all,
    drown_break and relocation_break say why management stopped, and
    _road_ele_TS is written every managed year so its last non-zero entry dates
    the last managed year.

    Args:
        cascade: A Cascade instance after its run.
        geometry: DomainGeometry describing the padded array.
        first_gis: First GIS domain carrying road.
        last_gis: Last GIS domain carrying road.

    Returns:
        A list of per-domain dicts for the managed domains only, each with gis,
        pad, drowned, relocation_blocked, last_managed_year, reason,
        overwash_removed_m3, dunes_rebuilt and relocations.
    """
    roadways = getattr(cascade, "roadways", None)
    if roadways is None:
        return []
    management = getattr(cascade, "roadway_management_module", None)

    rows = []
    for gis in range(first_gis, last_gis + 1):
        pad = geometry.gis_to_pad(gis)
        if not 0 <= pad < len(roadways) or roadways[pad] is None:
            continue
        if management is not None and not management[pad]:
            continue
        manager = roadways[pad]

        elevation_ts = np.asarray(getattr(manager, "_road_ele_TS", []),
                                  dtype=float)
        nonzero = np.nonzero(elevation_ts)[0]
        drowned = bool(getattr(manager, "drown_break", 0))
        blocked = bool(getattr(manager, "relocation_break", 0))
        rows.append({
            "gis": gis, "pad": pad,
            "drowned": drowned, "relocation_blocked": blocked,
            "last_managed_year": int(nonzero[-1]) + 1 if nonzero.size else 0,
            "reason": ("drowned" if drowned else
                       "no room to relocate" if blocked else
                       "managed throughout"),
            "overwash_removed_m3": float(np.sum(
                getattr(manager, "_road_overwash_volume", [0.0]))),
            "dunes_rebuilt": int(np.sum(
                getattr(manager, "_dunes_rebuilt_TS", [0]))),
            "relocations": int(np.sum(
                getattr(manager, "_road_relocated_TS", [0]))),
        })
    return rows


# =============================================================================
# Mid-run historical events
# =============================================================================

def apply_historical_event(cascade, event, geometry, relocations_enabled=True,
                           setback_check=None):
    """Applies one historical roadway event to a running CASCADE model.

    Called from inside the time loop when `event.year` comes up. Both event
    types mutate live model state, and each does so in a way that needs
    explaining:

    `RelocationEvent` adds its measured DISPLACEMENT to the setback the model
    is currently carrying. `roadway_manager` has been decrementing that setback
    by dune migration since t=0, so it already holds the modelled retreat;
    adding the displacement counts the retreat exactly once. An absolute
    setback referenced to an older dune line would count it twice.

    `road_ele` is deliberately left alone. A real relocation rebuilds the road
    at grade on new ground, so resetting it looks right -- but it would be an
    exact no-op. `road_ele` is initialised from the 2004 alignment, which IS
    the post-relocation road for events preceding 2004, and `roadway_manager`
    decrements it in the Lagrangian frame, so at year t it already holds
    `measured_2004 - sum(RSLR[0:t]) * 10`. Rebuilding at grade on that same
    alignment gives the identical number. This only breaks if a
    post-relocation alignment has a different measured elevation than the
    initial one, which is not the case here by construction.

    CASCADE's own relocation guards are evaluated but NOT obeyed. The
    model-driven relocation path refuses a move when the road would drown or
    the island is too narrow; this prescribed path does not, because these
    relocations actually happened and refusing them would be historically
    wrong. Each refusal is recorded in the returned row so the disagreement
    goes on the record rather than disappearing.

    `BridgeEvent` switches `cascade.roadway_management_module` off for its
    domains. It does NOT set `cascade.roadways[pad] = None`: `Cascade.update()`
    reads `self._roadways[iB3D].drown_break` on every domain with no None
    check, so nulling the object raises AttributeError on the next step.

    Args:
        cascade: A Cascade instance mid-run.
        event: A RelocationEvent or BridgeEvent whose year has arrived.
        geometry: DomainGeometry describing the padded array.
        relocations_enabled: Global toggle for relocation events; a False
            value skips RelocationEvents and reports them as skipped. Bridge
            events ignore it.
        setback_check: Optional {gis: measured_setback_m} for an independent
            cross-check, printed by the caller alongside the new setback.

    Returns:
        A list of per-domain dicts. Relocation rows carry gis, pad, kind,
        old_setback_m, displacement_m, new_setback_m, check_m and warnings;
        bridge rows carry gis, pad and kind. An event that is skipped or
        touches no managed domain returns an empty list.
    """
    import cascade.roadway_manager as _rm

    roadways = getattr(cascade, "roadways", None)
    if roadways is None:
        return []

    if isinstance(event, RelocationEvent):
        gis_range = sorted(event.displacement_m)
    else:
        gis_range = list(event.gis_domains)
    n_domains = len(cascade.barrier3d)
    pads = [(g, geometry.gis_to_pad(g)) for g in gis_range]
    pads = [(g, p) for g, p in pads if 0 <= p < n_domains]

    if isinstance(event, RelocationEvent):
        if not (relocations_enabled and event.enabled):
            return []
        return _apply_relocation(cascade, event, pads, roadways, _rm,
                                 setback_check or {})

    if not event.enabled:
        return []
    return _apply_bridge(cascade, pads, roadways)


def _apply_relocation(cascade, event, pads, roadways, _rm, setback_check):
    """Applies a RelocationEvent's displacement, probing CASCADE's guards."""
    rows = []
    b3d_all = getattr(cascade, "barrier3d", None)
    for gis, pad in pads:
        manager = roadways[pad]
        if manager is None or not hasattr(manager, "_road_setback"):
            continue

        old_setback = manager._road_setback
        displacement = event.displacement_m.get(gis)
        new_setback = (old_setback + float(displacement)
                       if displacement is not None
                       else manager._road_relocation_setback)

        warnings = []
        b3d = b3d_all[pad] if b3d_all is not None else None
        if b3d is not None:
            try:
                probe_ele, probe_drown = _rm.get_road_relocation_elevation(
                    time_index=manager._time_index,
                    xyz_interior_grid=b3d.InteriorDomain,
                    road_setback=new_setback,
                    road_width=manager._road_relocation_width,
                    dx=10, dy=10, dz=10,
                )
                if probe_drown:
                    warnings.append(
                        f"get_road_relocation_elevation REFUSES (road would "
                        f"be {probe_ele:.2f} m MSL)")
            except IndexError:
                warnings.append("get_road_relocation_elevation: setback is "
                                "past the end of the interior grid")
            island_width_m = float(b3d.InteriorWidth_AvgTS[-1]) * 10.0
            needed = new_setback + 2 * manager._road_relocation_width
            if needed > island_width_m:
                warnings.append(
                    f"road_relocation_checks REFUSES (needs {needed:.0f} m, "
                    f"island is {island_width_m:.0f} m)")

        manager._road_setback = new_setback
        manager._road_relocation_setback = new_setback   # keep the two in sync
        manager._road_setback_TS[manager._time_index - 1] = new_setback

        rows.append(dict(
            gis=gis, pad=pad, kind="relocation",
            old_setback_m=float(old_setback),
            displacement_m=float(displacement or 0.0),
            new_setback_m=float(new_setback),
            check_m=setback_check.get(gis),
            warnings=warnings,
        ))
    return rows


def _apply_bridge(cascade, pads, roadways):
    """Switches roadway management off for a BridgeEvent's domains.

    Flips `cascade.roadway_management_module` rather than tracking the pads in
    a set the model never reads -- which is what an earlier version did, so the
    bridged domains kept full road management for the rest of the run.
    """
    flags = list(cascade.roadway_management_module)
    rows = []
    for gis, pad in pads:
        if roadways[pad] is None:
            continue
        flags[pad] = False
        rows.append(dict(gis=gis, pad=pad, kind="bridge"))
    cascade.roadway_management_module = flags
    return rows
