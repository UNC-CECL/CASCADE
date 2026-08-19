"""Beach and dune forcing for a CASCADE run: nourishment schedule, overwash filter.

Everything CASCADE's `beach_dune_manager` needs, prepared before the run and
checked against what the module actually did afterwards. Nothing here is
site-specific -- domain geometry arrives as a `DomainGeometry`, and the
Hatteras instances (project extents, volumes, community zones) live in
`hatteras_site_config`.

Four things about `BeachDuneManager` this module exists to get right:

- **`overwash_filter` is a PERCENT, not a fraction.** `filter_overwash` divides
  it by 100, and the docstring cites 40-90 % from Rogers et al. (2015):
  residential to commercial. A value of 0.4 filters 0.4 % of overwash, which
  is indistinguishable from no filtering. `BeachDuneConfig` rejects the
  fraction scale rather than letting it pass silently.

- **The per-year volume must be written to the Cascade object, not the
  manager.** `cascade.update()` copies its own `_nourishment_volume[iB3D]`
  into each manager on every time step, immediately before calling
  `BeachDuneManager.update()` (cascade.py:772). A volume assigned straight
  onto `cascade.nourishments[i]` is therefore overwritten before it is ever
  used, and every nourishment silently applies the Cascade init default
  instead. `NourishmentSchedule.apply_to_cascade` writes
  `cascade.nourishment_volume`, which survives.

- **Enabling the module is not the same as scheduling a nourishment.** With
  `beach_nourishment_module=True`, a domain runs overwash filtering, the 50 m
  community-width drowning check, and fixed-dune-line dynamics EVERY year of
  the run -- `nourish_now` only controls whether sand is added. There is no
  way to get the fill without the rest.

- **The manager's time series are offset by one.** `update_dune_domain`
  increments `barrier3d.time_index` before the managers run, so a manager
  writing at `time_index - 1` lands on index `year - start_year + 1`.
  `NourishmentSchedule.time_index` is the single place that conversion lives.
"""

import dataclasses

import numpy as np

CUBIC_YARDS_TO_CUBIC_METERS = 0.764555


@dataclasses.dataclass(frozen=True)
class BeachDuneConfig:
    """Percent-scale settings for CASCADE's beach_dune_manager.

    Attributes:
        community_overwash_filter_pct: Percent of overwash deposition removed
            from the interior in a developed community and returned to the
            shoreface. Rogers et al. (2015) give 40-90 %, residential to
            commercial.
        default_overwash_filter_pct: Percent for every other domain --
            undeveloped ground filters nothing, so this is normally 0.
        overwash_to_dune_pct: Percent of the remaining overwash bulldozed onto
            the dunes. filter + to_dune must stay under 100.
        minimum_community_width_m: BeachDuneManager's own hard-coded threshold,
            repeated here so pre-flight checks can use it. Below this average
            interior width the community is abandoned and management stops.
    """

    community_overwash_filter_pct: float = 40.0
    default_overwash_filter_pct: float = 0.0
    overwash_to_dune_pct: float = 9.0
    minimum_community_width_m: float = 50.0

    def __post_init__(self):
        for name in ("community_overwash_filter_pct",
                     "default_overwash_filter_pct"):
            value = getattr(self, name)
            if 0.0 < value < 1.0:
                raise ValueError(
                    f"{name}={value} looks like a fraction. CASCADE's "
                    f"overwash_filter is a PERCENT (filter_overwash divides "
                    f"by 100), so {value} filters {value:.1f}% of overwash, "
                    f"not {value * 100:.0f}%. Pass {value * 100:.0f}.")
        total = self.community_overwash_filter_pct + self.overwash_to_dune_pct
        if total >= 100:
            raise ValueError(
                f"overwash_filter + overwash_to_dune must be under 100%, "
                f"got {total}")


DEFAULT_BEACH_DUNE = BeachDuneConfig()


@dataclasses.dataclass(frozen=True)
class NourishmentProject:
    """One historical beach nourishment project.

    Volume is the project total as reported by the permitting agency, in cubic
    yards, spread evenly across the project's domains. Even spreading is an
    assumption: real fill templates taper at the ends, but the reported totals
    are not resolved by domain.

    Attributes:
        name: Project name, used in printed summaries and run metadata.
        year: Calendar year the fill was placed.
        gis_domains: GIS domains the project covered.
        volume_cubic_yards: Reported project total.
        note: Human-readable description for run metadata.
        enabled: Whether to include this project in a schedule.
    """

    name: str
    year: int
    gis_domains: tuple
    volume_cubic_yards: float
    note: str = ""
    enabled: bool = True

    @property
    def volume_m3_total(self):
        """Project total in cubic meters."""
        return self.volume_cubic_yards * CUBIC_YARDS_TO_CUBIC_METERS

    def volume_m3_per_domain(self):
        """Cubic meters delivered to each of the project's domains."""
        return self.volume_m3_total / len(self.gis_domains)

    def volume_m3_per_m(self, domain_spacing_m):
        """Volume per unit alongshore length, the unit CASCADE wants.

        Args:
            domain_spacing_m: Alongshore width of one domain, in meters.

        Returns:
            Nourishment volume in m^3/m.
        """
        return self.volume_m3_per_domain() / domain_spacing_m


@dataclasses.dataclass(frozen=True)
class NourishmentSchedule:
    """Per-year nourishment flags and volumes on the padded array.

    Attributes:
        nourish_now: Mapping of calendar year to a 0/1 float array of length
            geometry.total_domains, matching CASCADE's nourish_now.
        volume_m3_per_m: Mapping of calendar year to a float array of length
            geometry.total_domains, in m^3/m.
        projects: Projects that fall inside the period.
        skipped: (project, reason) pairs for projects that do not.
        geometry: DomainGeometry the arrays are indexed by.
        start_year: First model year.
        end_year: Last model year.
    """

    nourish_now: dict
    volume_m3_per_m: dict
    projects: tuple
    skipped: tuple
    geometry: object
    start_year: int
    end_year: int

    @property
    def years(self):
        """Calendar years with at least one nourished domain, sorted."""
        return tuple(sorted(self.nourish_now))

    @property
    def nourished_gis(self):
        """Every GIS domain nourished at any point in the period, sorted."""
        gis_ids = set()
        for project in self.projects:
            gis_ids.update(project.gis_domains)
        return tuple(sorted(gis_ids))

    @property
    def is_empty(self):
        """True when no project falls inside this period."""
        return not self.projects

    def time_index(self, year):
        """Index into a BeachDuneManager time series for a calendar year.

        Barrier3D's `update_dune_domain` increments `time_index` before the
        management modules run, so a manager writing at `time_index - 1` for
        the first model year lands on index 1, not 0.

        Args:
            year: Calendar year.

        Returns:
            The index into `_nourishment_TS`, `_beach_width`, and the other
            per-time-step arrays on BeachDuneManager.
        """
        return year - self.start_year + 1

    def events(self):
        """Every scheduled nourishment, one row per domain-year.

        Returns:
            A list of dicts with year, time_index, pad, gis, and
            volume_m3_per_m, sorted by year then GIS domain.
        """
        rows = []
        for year in self.years:
            flags = self.nourish_now[year]
            volumes = self.volume_m3_per_m[year]
            for pad in np.flatnonzero(flags == 1):
                pad = int(pad)
                rows.append(dict(
                    year=year,
                    time_index=self.time_index(year),
                    pad=pad,
                    gis=int(self.geometry.pad_to_gis(pad)),
                    volume_m3_per_m=float(volumes[pad]),
                ))
        return sorted(rows, key=lambda row: (row["year"], row["gis"]))

    def apply_to_cascade(self, cascade, year, default_volume_m3_per_m=0.0):
        """Sets this year's nourishment flags and volumes on a live Cascade.

        Call this immediately BEFORE `cascade.update()` for the matching model
        year. Both arrays are rewritten in full every year, so nothing carries
        over from the previous one.

        The volume goes to `cascade.nourishment_volume` rather than to the
        BeachDuneManager instances: `cascade.update()` copies its own list into
        each manager before calling it, so a value written onto the manager is
        overwritten before use. See the module docstring.

        Args:
            cascade: A live Cascade instance.
            year: Calendar year about to be stepped.
            default_volume_m3_per_m: Volume carried by domains with no event.
                Unused by the model -- `nourish_now` is 0 there and no
                nourishment interval is set -- so 0 keeps it unambiguous.

        Returns:
            A list of dicts (year, pad, gis, volume_m3_per_m) for the domains
            nourished this year, in GIS order. Empty when nothing is scheduled.

        Raises:
            ValueError: If the schedule is wider than the Cascade grid.
        """
        n_domains = len(cascade.barrier3d)
        if n_domains != self.geometry.total_domains:
            raise ValueError(
                f"Cascade has {n_domains} domains, schedule was built for "
                f"{self.geometry.total_domains}")

        nourish_now = np.zeros(n_domains)
        volumes = [float(default_volume_m3_per_m)] * n_domains
        applied = []

        if year in self.nourish_now:
            flags = self.nourish_now[year]
            year_volumes = self.volume_m3_per_m[year]
            for pad in np.flatnonzero(flags == 1):
                pad = int(pad)
                nourish_now[pad] = 1
                volumes[pad] = float(year_volumes[pad])
                applied.append(dict(
                    year=year,
                    pad=pad,
                    gis=int(self.geometry.pad_to_gis(pad)),
                    volume_m3_per_m=volumes[pad],
                ))

        cascade.nourish_now = nourish_now
        cascade.nourishment_volume = volumes
        return sorted(applied, key=lambda row: row["gis"])


def build_schedule(projects, geometry, start_year, end_year):
    """Expands nourishment projects onto per-year padded arrays.

    Projects dated outside [start_year, end_year] are skipped rather than
    raising, so one project list serves every period -- a 1984-2004 run simply
    returns an empty schedule.

    Args:
        projects: Iterable of NourishmentProject.
        geometry: DomainGeometry describing the padded array.
        start_year: First model year.
        end_year: Last model year.

    Returns:
        A NourishmentSchedule.

    Raises:
        ValueError: If a project covers a GIS domain outside the real span, or
            names no domains, or carries a negative volume.
    """
    nourish_now = {}
    volume_m3_per_m = {}
    in_period = []
    skipped = []

    for project in projects:
        if not project.enabled:
            skipped.append((project, "disabled"))
            continue
        if not project.gis_domains:
            raise ValueError(f"{project.name} names no GIS domains")
        if project.volume_cubic_yards < 0:
            raise ValueError(
                f"{project.name} has a negative volume "
                f"({project.volume_cubic_yards} cy)")
        for gis in project.gis_domains:
            if not geometry.first_gis_id <= gis <= geometry.last_gis_id:
                raise ValueError(
                    f"{project.name} covers GIS {gis}, outside the real span "
                    f"{geometry.first_gis_id}-{geometry.last_gis_id}")
        if not start_year <= project.year <= end_year:
            skipped.append(
                (project, f"{project.year} outside {start_year}-{end_year}"))
            continue

        in_period.append(project)
        if project.year not in nourish_now:
            nourish_now[project.year] = np.zeros(geometry.total_domains)
            volume_m3_per_m[project.year] = np.zeros(geometry.total_domains)

        per_m = project.volume_m3_per_m(geometry.domain_spacing_m)
        for gis in project.gis_domains:
            pad = geometry.gis_to_pad(gis)
            # Two projects overlapping in one domain-year would silently
            # clobber each other; sum instead, which is what placing both
            # volumes on the same beach means.
            nourish_now[project.year][pad] = 1
            volume_m3_per_m[project.year][pad] += per_m

    return NourishmentSchedule(
        nourish_now=nourish_now,
        volume_m3_per_m=volume_m3_per_m,
        projects=tuple(in_period),
        skipped=tuple(skipped),
        geometry=geometry,
        start_year=start_year,
        end_year=end_year,
    )


def build_overwash_filter(geometry, community_zones,
                          config=DEFAULT_BEACH_DUNE):
    """Per-domain overwash filter percent, for Cascade's `overwash_filter`.

    Only developed ground filters overwash, so the community percent applies
    inside the community zones and the default (normally 0) everywhere else --
    including nourished domains outside a village, which are road corridor and
    refuge rather than development.

    Args:
        geometry: DomainGeometry describing the padded array.
        community_zones: Iterable of inclusive (first_gis, last_gis) spans.
        config: BeachDuneConfig supplying the two percentages.

    Returns:
        A list of geometry.total_domains percentages.
    """
    values = [config.default_overwash_filter_pct] * geometry.total_domains
    for first_gis, last_gis in community_zones:
        for gis in range(first_gis, last_gis + 1):
            values[geometry.gis_to_pad(gis)] = (
                config.community_overwash_filter_pct)
    return values


def build_beach_dune_management_on(geometry, community_zones, nourished_gis,
                                   enabled=True):
    """Which domains run CASCADE's beach_dune_manager.

    The union of two footprints drawn from different sources: the permanent
    community zones, which need overwash filtering every year, and the
    nourishment project extents, which need the module on at all to receive
    fill. The project extents are wider than the villages, so the union is
    larger than either.

    Everything in this list runs the module's always-on behaviour for the whole
    run -- fixed dune line, community-width drowning check -- not just in
    nourishment years. See the module docstring.

    Args:
        geometry: DomainGeometry describing the padded array.
        community_zones: Iterable of inclusive (first_gis, last_gis) spans.
        nourished_gis: Iterable of GIS domains that receive fill.
        enabled: Global switch; False turns the module off everywhere, which
            is what a natural-dynamics run wants. Mirrors the same argument on
            roadway.build_roadway_management_on.

    Returns:
        A list of geometry.total_domains booleans.
    """
    on = [False] * geometry.total_domains
    if not enabled:
        return on

    for first_gis, last_gis in community_zones:
        for gis in range(first_gis, last_gis + 1):
            on[geometry.gis_to_pad(gis)] = True
    for gis in nourished_gis:
        on[geometry.gis_to_pad(gis)] = True
    return on


def find_double_managed(beach_dune_on, roadway_on, geometry):
    """GIS domains where the roadway and beach-dune managers both run.

    CASCADE does not guard against this. Both loops run in `cascade.update()`,
    roadway first, and the combination has two consequences worth naming:
    overwash is removed twice over (bulldozed, then the survivors filtered),
    and BeachDuneManager holds `dune_migration_on` False, which leaves
    `ShorelineChangeTS` at 0 -- the value RoadwayManager reads to decrement the
    road setback. The road therefore stops retreating in these domains.

    Args:
        beach_dune_on: Per-domain booleans for beach_nourishment_module.
        roadway_on: Per-domain booleans for roadway_management_module.
        geometry: DomainGeometry describing the padded array.

    Returns:
        A sorted tuple of GIS domain IDs, empty if the footprints are disjoint.
    """
    gis_ids = []
    for pad in range(geometry.total_domains):
        if beach_dune_on[pad] and roadway_on[pad]:
            gis = int(geometry.pad_to_gis(pad))
            if geometry.first_gis_id <= gis <= geometry.last_gis_id:
                gis_ids.append(gis)
    return tuple(sorted(gis_ids))


def audit_schedule(schedule, beach_dune_on, config=DEFAULT_BEACH_DUNE):
    """Pre-flight check that the schedule can actually reach the model.

    Catches the failures that produce a plausible-looking run rather than an
    error: a domain scheduled for fill whose module is off (the nourishment is
    silently dropped), a volume outside any real project's range, or a schedule
    year outside the run.

    Args:
        schedule: The NourishmentSchedule about to be run.
        beach_dune_on: Per-domain booleans for beach_nourishment_module.
        config: BeachDuneConfig, for the reported drowning threshold.

    Returns:
        A dict with 'ok', 'module_off' (scheduled but unreachable, as GIS IDs),
        'implausible_volume' (event rows), 'out_of_period' (event rows), and
        'notes' (human-readable strings).
    """
    module_off = []
    implausible = []
    out_of_period = []
    notes = []

    for row in schedule.events():
        if not beach_dune_on[row["pad"]]:
            module_off.append(row["gis"])
        # Real projects run roughly 20-2000 m^3/m; outside that is a unit slip
        # (total m^3 not divided by domain length, or cy never converted).
        if not 1.0 <= row["volume_m3_per_m"] <= 5000.0:
            implausible.append(row)
        if not schedule.start_year <= row["year"] <= schedule.end_year:
            out_of_period.append(row)

    if schedule.is_empty:
        notes.append(
            f"No nourishment projects fall in {schedule.start_year}-"
            f"{schedule.end_year}; the module still runs its always-on "
            f"behaviour wherever beach_dune_management_on is True.")
    for project, reason in schedule.skipped:
        notes.append(f"skipped {project.name}: {reason}")

    return dict(
        ok=not (module_off or implausible or out_of_period),
        module_off=sorted(set(module_off)),
        implausible_volume=implausible,
        out_of_period=out_of_period,
        notes=notes,
    )


def verify_nourishment(cascade, schedule, beach_dune_on,
                       tolerance_m3_per_m=0.01):
    """Checks what BeachDuneManager actually did against the schedule.

    Reads each manager's own `_nourishment_TS` and `_nourishment_volume_TS`,
    which are written inside the nourishment branch of
    `BeachDuneManager.update` from the volume passed to `shoreface_nourishment`
    (beach_dune_manager.py:794-795). So this reports the volume the model
    used -- not the volume the run script intended, which is what a
    schedule-side log records.

    That distinction is the point of the function. Writing the volume onto the
    manager instead of onto the Cascade object produces a run where every
    nourishment applies the Cascade init default while the log looks correct;
    only the manager's own record disagrees.

    Args:
        cascade: The Cascade instance after the run.
        schedule: The NourishmentSchedule the run was driven by.
        beach_dune_on: Per-domain booleans for beach_nourishment_module.
        tolerance_m3_per_m: Allowed volume difference before it counts as
            wrong.

    Returns:
        A dict with 'ok', 'confirmed', 'missing', 'wrong_volume',
        'unexpected', 'abandoned', and 'notes'. Every list holds dicts with
        year, gis, and the expected/actual volumes where relevant.
    """
    confirmed = []
    missing = []
    wrong_volume = []
    unexpected = []
    abandoned = []
    notes = []

    expected = {(row["gis"], row["year"]): row for row in schedule.events()}

    for pad in range(schedule.geometry.total_domains):
        if not beach_dune_on[pad]:
            continue
        manager = cascade.nourishments[pad]
        gis = int(schedule.geometry.pad_to_gis(pad))

        if manager.narrow_break:
            abandoned.append(dict(gis=gis, pad=pad))

        nourished_TS = np.asarray(manager._nourishment_TS)
        volume_TS = np.asarray(manager._nourishment_volume_TS)

        for index in np.flatnonzero(nourished_TS == 1):
            index = int(index)
            year = schedule.start_year + index - 1
            key = (gis, year)
            actual = float(volume_TS[index])
            if key not in expected:
                unexpected.append(
                    dict(gis=gis, year=year, volume_m3_per_m=actual))
                continue
            wanted = expected[key]["volume_m3_per_m"]
            row = dict(gis=gis, year=year, expected_m3_per_m=wanted,
                       actual_m3_per_m=actual)
            if abs(actual - wanted) > tolerance_m3_per_m:
                wrong_volume.append(row)
            else:
                confirmed.append(row)

    seen = {(row["gis"], row["year"])
            for row in confirmed + wrong_volume}
    for key, row in expected.items():
        if key not in seen:
            missing.append(dict(gis=row["gis"], year=row["year"],
                                expected_m3_per_m=row["volume_m3_per_m"]))

    if wrong_volume:
        uniform = {round(row["actual_m3_per_m"], 6) for row in wrong_volume}
        if len(uniform) == 1:
            notes.append(
                f"Every wrong volume is the same value "
                f"({uniform.pop()} m^3/m), which is the signature of the "
                f"Cascade init default overwriting the schedule -- the "
                f"per-year volume was written onto the BeachDuneManager "
                f"instead of onto cascade.nourishment_volume.")
    if abandoned:
        notes.append(
            f"{len(abandoned)} domain(s) hit the 50 m community-width "
            f"drowning check and stopped being managed; any nourishment "
            f"scheduled after that year could not be delivered.")

    return dict(
        ok=not (missing or wrong_volume or unexpected),
        confirmed=sorted(confirmed, key=lambda r: (r["year"], r["gis"])),
        missing=sorted(missing, key=lambda r: (r["year"], r["gis"])),
        wrong_volume=sorted(wrong_volume, key=lambda r: (r["year"], r["gis"])),
        unexpected=sorted(unexpected, key=lambda r: (r["year"], r["gis"])),
        abandoned=abandoned,
        notes=notes,
    )


def verify_setbacks_frozen(cascade, double_managed_gis, geometry):
    """Confirms the road setback stopped moving in double-managed domains.

    The predicted consequence of running both managers in one domain, checked
    against the run rather than assumed: BeachDuneManager holds
    `dune_migration_on` False, Barrier3D then leaves `ShorelineChangeTS` at 0,
    and RoadwayManager only decrements the setback when that value is non-zero.

    Args:
        cascade: The Cascade instance after the run.
        double_managed_gis: GIS domains reported by `find_double_managed`.
        geometry: DomainGeometry describing the padded array.

    Returns:
        A list of dicts (gis, setback_start_m, setback_end_m, frozen), one per
        double-managed domain, sorted by GIS ID.
    """
    rows = []
    for gis in double_managed_gis:
        pad = geometry.gis_to_pad(gis)
        setback_TS = np.asarray(cascade.roadways[pad]._road_setback_TS)
        # Trailing zeros are unwritten years after the road stopped being
        # managed, not a setback of zero.
        written = setback_TS[: np.max(np.flatnonzero(setback_TS != 0)) + 1] \
            if np.any(setback_TS != 0) else setback_TS[:1]
        rows.append(dict(
            gis=gis,
            setback_start_m=float(written[0]),
            setback_end_m=float(written[-1]),
            frozen=bool(np.allclose(written, written[0])),
        ))
    return sorted(rows, key=lambda row: row["gis"])
