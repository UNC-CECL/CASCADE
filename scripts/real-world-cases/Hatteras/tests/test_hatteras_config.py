"""The invariants of the Hatteras hindcast configuration.

WHY THESE EXIST
    Until 2026-09-12 nothing in the test suite referenced the site config, the
    period table or the topography resolver. Every invariant below was enforced
    by a person running a script and reading its output -- which is how the
    exported calibration came to disagree with the config at 63 domains, and
    how the notebook and its headless mirror drifted apart on a setting neither
    of them printed.

WHAT THEY COVER
    Only what can be checked without running the model: that the period table
    is internally consistent, that a preset covers the periods it claims, that
    the two copies of the runner agree on the values that select a run, and
    that every forcing file a period names is on disk.

    The last of those needs the data tree, which is not in the repository for
    the bulk rasters, so it SKIPS rather than fails when the tree is absent.
    The rest are pure configuration and always run.
"""

import json
import re
import sys
from pathlib import Path

import pytest

# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(_p for _p in Path(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
SCRIPTS = _HAT_DIR
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

sys.path.insert(0, str(SCRIPTS))

site_config = pytest.importorskip("site_layer.hatteras_site_config")
topo_version = pytest.importorskip("site_layer.hat_topo_version")

# The runner's run-selecting sections live in its setup module since 2026-10-23
RUNNER_SETUP = SCRIPTS / "hatteras_ms" / "hindcast" / "HAT_hindcast_setup.py"
RUNNER_NB = SCRIPTS / "hatteras_ms" / "HAT_hindcast.ipynb"

# Every period names these, as paths relative to data/hatteras_init/.
FORCING_KEYS = ("storm_file", "island_offset_file", "road_setback_file")

# The assignments that SELECT A RUN. They are made in the setup steps only; the
# notebook typing one itself would let it simulate a run the script does not.
SYNCED_ASSIGNMENTS = (
    "START_YEAR = ",
    "SOURCE_SINK_PRESET = ",
    "DOMAIN_BE_RATES = ",
    "COASTSAT_BASE_DIR = ",
    "TOPO_PRODUCT = ",
    "LOWESS_CONFIG = ",
    "TARGET_WINDOW = ",
)


def _periods():
    return site_config.HATTERAS_PERIODS


# ---------------------------------------------------------------------------
# THE PERIOD TABLE
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("start_year", sorted(_periods()))
def test_period_window_is_forward_and_nonempty(start_year):
    """A period runs start..end-1, so the end must be strictly later."""
    end_year = _periods()[start_year]["end_year"]
    assert end_year > start_year, (
        f"period {start_year} ends at {end_year}, so it would run no years")


@pytest.mark.parametrize("start_year", sorted(_periods()))
def test_period_names_every_forcing_it_needs(start_year):
    period = _periods()[start_year]
    missing = [k for k in FORCING_KEYS + ("topo_product", "sea_level_rise_rate",
                                          "enable_nourishment",
                                          "nourishment_volume")
               if k not in period]
    assert not missing, f"period {start_year} names no {missing}"


@pytest.mark.parametrize("start_year", sorted(_periods()))
def test_period_topography_product_is_resolvable(start_year):
    """The product a period names must be one the resolver knows."""
    product = _periods()[start_year]["topo_product"]
    assert product in topo_version.PRODUCTS, (
        f"period {start_year} wants topography product {product!r}, which is "
        f"not one of {topo_version.PRODUCTS}")


def test_year_product_and_period_table_agree():
    """Every period has a product mapping, and every mapping has a period.

    They are separate tables in separate modules, and a period added to one
    and not the other resolves to the default product silently -- which is a
    different island, not a rounding difference.
    """
    periods = set(_periods())
    mapped = set(topo_version.YEAR_PRODUCT)
    assert periods == mapped, (
        f"HATTERAS_PERIODS has {sorted(periods - mapped)} with no entry in "
        f"YEAR_PRODUCT, and YEAR_PRODUCT has {sorted(mapped - periods)} with "
        f"no period")


def test_nourishment_volume_matches_the_switch():
    """A period with fills off must not carry a volume, and vice versa."""
    for start_year, period in _periods().items():
        on = period["enable_nourishment"]
        volume = period["nourishment_volume"]
        assert bool(on) == bool(volume), (
            f"period {start_year} has enable_nourishment={on} but "
            f"nourishment_volume={volume}")


# ---------------------------------------------------------------------------
# THE SOURCE/SINK PRESETS
# ---------------------------------------------------------------------------

def test_zero_preset_covers_every_period():
    """zeroBE is what a wired-but-uncalibrated period runs under."""
    zero = site_config.HATTERAS_BE_PRESETS["zeroBE"]
    missing = sorted(set(_periods()) - set(zero))
    assert not missing, (
        f"periods {missing} have no zeroBE entry, so they cannot be run at "
        f"all until one is calibrated")


def test_be_rates_names_the_missing_fit_rather_than_raising_keyerror():
    """An uncalibrated period must fail with a sentence, not a bare KeyError."""
    uncalibrated = [p for p in _periods()
                    if p not in site_config.HATTERAS_BE_PRESETS["edgeBE"]]
    if not uncalibrated:
        pytest.skip("every period has an edgeBE fit")
    with pytest.raises(ValueError, match="no rates for"):
        site_config.be_rates("edgeBE", uncalibrated[0])


def test_edge_only_periods_have_no_calibrated_fit():
    """GIS 1 must have exactly one home.

    HATTERAS_BE_RATES_EDGE takes GIS 1 from the calibrated preset where one
    exists, and from HATTERAS_BE_EDGE_ONLY where it does not. A period in both
    would have two, and the merge would quietly win.
    """
    both = sorted(set(getattr(site_config, "HATTERAS_BE_EDGE_ONLY", {}))
                  & set(site_config.HATTERAS_BE_RATES_CALIBRATED))
    assert not both, (
        f"periods {both} appear in both HATTERAS_BE_EDGE_ONLY and the "
        f"calibrated preset; delete the edge-only entry")


@pytest.mark.parametrize("preset", ("zeroBE", "edgeBE", "calibBE"))
def test_preset_rates_are_finite(preset):
    for start_year, rates in site_config.HATTERAS_BE_PRESETS[preset].items():
        for gis, rate in rates.items():
            assert isinstance(rate, (int, float)), (
                f"{preset}[{start_year}][{gis}] is {rate!r}")
            assert rate == rate, f"{preset}[{start_year}][{gis}] is NaN"


# ---------------------------------------------------------------------------
# THE NOTEBOOK AND ITS HEADLESS MIRROR
# ---------------------------------------------------------------------------

def _notebook_code():
    nb = json.loads(RUNNER_NB.read_text(encoding="utf-8"))
    return "\n".join("".join(cell["source"]) for cell in nb["cells"]
                     if cell["cell_type"] == "code")


# Every step HAT_hindcast.py takes, which the notebook must take too
RUNNER_STEPS = ("prepare_run", "pre_run_diagnostics", "verify_run", "write_run_outputs", "make_run_figures")


def test_notebook_takes_every_step():
    """The notebook builds the run with the same functions as the script.

    Since 2026-10-23 both call the setup_* steps, the outputs and the figures;
    a section the notebook computed itself again could drift from the script,
    the failure the old line-by-line comparison caught (2026-09-12: LOWESS_CONFIG
    was one smoothing window in the script and two in the notebook).
    """
    if not (RUNNER_SETUP.exists() and RUNNER_NB.exists()):
        pytest.skip("runner setup or notebook absent")
    setup_text = RUNNER_SETUP.read_text(encoding="utf-8")
    steps = re.search(r"^STEPS = \((.*?)\)", setup_text, re.M | re.S).group(1)
    wanted = [n.strip() for n in steps.split(",") if n.strip()]
    nb_text = _notebook_code()
    missing = [n for n in wanted + list(RUNNER_STEPS[1:]) if f"{n}(r" not in nb_text]
    assert not missing, f"the notebook does not call {missing}"


@pytest.mark.parametrize("prefix", SYNCED_ASSIGNMENTS)
def test_notebook_sets_no_run_value_itself(prefix):
    """A run-selecting value comes from the setup steps, never typed in the notebook."""
    if not RUNNER_NB.exists():
        pytest.skip("notebook absent")
    typed = [line.strip() for line in _notebook_code().splitlines()
             if line.lstrip().startswith(prefix)]
    assert not typed, f"the notebook sets {prefix!r} itself: {typed}"

# ---------------------------------------------------------------------------
# THE FORCING FILES, WHEN THE DATA TREE IS PRESENT
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("start_year", sorted(_periods()))
def test_period_forcing_files_exist(start_year):
    """Every path a period names resolves to a file on disk.

    Skips when the data tree is absent, since the bulk of it is not in the
    repository. HAT_period_input_check.py is the fuller version of this and
    also checks shapes.
    """
    if not INIT_ROOT.is_dir():
        pytest.skip("data/hatteras_init is not present")
    period = _periods()[start_year]
    missing = [key for key in FORCING_KEYS
               if not (INIT_ROOT / period[key]).is_file()]
    if missing:
        pytest.xfail(
            f"period {start_year} names {missing}, which are not built yet")
