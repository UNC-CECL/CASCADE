# Model units: what the hindcast supplies, and what the models expect

Reference for every quantity that crosses into or out of CASCADE / Barrier3D /
BRIE in the Hatteras hindcast. Audited 2026-09-27 against the option A matrix
(`scripts/hatteras_ms/HAT_hindcast_1984_2024.py`, the headless mirror of the
notebook) and the base-model source. **Result: consistent end to end.**

Re-runnable check for the elevation rows:
`python scripts/input_prep/HAT_units_datum_check.py`. It checks only the
elevation rows. The other rows were traced by hand; the line references below
are where to look if one changes.

Base models, as installed (editable):
- Barrier3D: `C:\Users\hanna\PycharmProjects\Barrier3D` (commit 49fd069)
- BRIE: `C:\Users\hanna\PycharmProjects\brie`
- CASCADE: `cascade/` in this repo

---

## The one rule to remember

Barrier3D computes in **decametres (dam) relative to MHW**. It converts
*scalars* from the parameter yaml for you (`load_input.py`: `/10`, and
`/10 - MHW` for elevations). It converts **nothing grid-shaped**: elevation and
dune arrays must arrive already in dam. BRIE works in **metres**; CASCADE's
coupler multiplies by 10 going to BRIE and divides by 10 coming back.

Datums: MHW = 0.36 m NAVD88 (Duck, NOAA 8651370). Berm = 1.7 m NAVD88
= 1.34 m MHW = 0.134 dam MHW.

---

## Inputs

### Elevations (datum-sensitive)

| Quantity | Hindcast supplies | Model expects | Converted where | Status |
|---|---|---|---|---|
| `elevation_file` (.npy per domain) | dam, MHW-relative | dam MHW | nowhere (`load_elevation` is `np.load`) | ✅ checked on arrays |
| `dune_file` (.npy per domain) | dam, height **above berm** | dam above berm | nowhere (`load_input.py:249`) | ✅ checked on arrays |
| `berm_elevation` = 1.7 | m NAVD88 | m NAVD88 | `load_input.py:241` `/10 - MHW` | ✅ |
| `MHW` = 0.36 | m NAVD88 | m NAVD88 | `load_input.py:227` `/10` | ✅ |
| dune ceiling (since 2026-09-28): `DuneCeilingFromStart: true`, `DuneCeilingFloor: 0.5` | each dune cell's own starting crest (from the dune `.npy`), floored 0.5 m above the berm | dam above the berm, per cell | `Barrier3d.__init__`: `max(DuneDomain[0].max(axis=1), floor / 10)` | ✅ requires Barrier3D `hatteras/adopted` |
| `Dmaxel` = 5.5 | m NAVD88 | m NAVD88 | `load_input.py:304` `/10 - MHW`; replaced by the median cell ceiling when the per-cell ceiling is on | ✅ explicit since 2026-09-28. Before that it was never set, so Barrier3D's 3.4 m NAVD88 default (3.04 m MHW, Virginia) capped every dune |
| `road_ele` (per domain) | m, MHW-relative | m MHW | `roadway_manager.bulldoze`: `/dz` | ✅ matches topo at the road, all 4 periods |
| `dune_design_elevation` = 3.0 | m MHW | m MHW | roadway_manager; floored at `BermEl*10 + 1.0` | ✅ |
| `dune_minimum_elevation` = 0.0 | m MHW | m MHW | floored at `BermEl*10 + 0.3` = 1.64 m | ✅ (floor governs) |
| BRIE `h_b_crit` | derived: berm − MHW = 1.34 m | m | `cascade.py` → `BrieCoupler` | ✅ |

### Storms

| Column | Hindcast supplies | Model expects | Status |
|---|---|---|---|
| `time` | model year, 1-based | model year | ✅ |
| `Rhigh`, `Rlow` | (TWL − 0.36 m) / 10 → dam MHW | dam MHW (compared with dune crest + berm, `barrier3d.py:1407`) | ✅ |
| `period` | s (WIS Tp at peak TWL) | s | ✅ |
| `duration` | hours with TWL above berm | hours | ✅ |

Built by `scripts/input_prep/3-env-forcings/3-storms/historical_storm_creation_v3_HAT.py`
from Duck water levels (m NAVD88, metric) + WIS Hs (m) / Tp (s). The storm threshold
(berm 1.7, MHW 0.36) matches the runner's datums.

### Horizontal geometry and rates

| Quantity | Hindcast supplies | Model expects | Converted where | Status |
|---|---|---|---|---|
| `shoreline_offset` | m (`offset_mode: metres`) | m, added to BRIE `x_t`, `x_s` | `brie_coupler.offset_shoreline`, no conversion | ✅ (`asrun` mode is the old ÷10 error) |
| `background_erosion` (BE) | m/yr, **+ = accretion/source** | m/yr, (−) erosion | `load_input.py:311` `Rat / -10` → `Qat`; `x_s_dt` is + landward, so + input = seaward | ✅ |
| `sea_level_rise_rate` | m/yr (0.004 / 0.006 / 0.007) | m/yr | `load_input.py:319` `/10` | ✅ |
| `road_setback`, `road_width` = 20, relocation setback | m | m | `bulldoze`: `/dy`, `/dx` | ✅ |
| Groin trapping rate M | m/yr | m per step on `x_s_dt` | `x_s_dt` is in m at the callback (`brie_coupler.batchB3D` ×10) | ✅ |
| Domain length | 500 m (BRIE `dy`, fixed) | m | `BarrierLength` `/10` | ✅ real spacing ~502.6 m; model uses 500 throughout |

### Waves (BRIE)

| Quantity | Value (option A) | Unit | Status |
|---|---|---|---|
| `wave_height` (Hs) | 2.0 | m, deep water | ✅ also sets shoreface depth `d_sf = 8.9 * Hs` |
| `wave_period` (Tp) | 7.5 | s | ✅ |
| `wave_asymmetry` | 0.6 | fraction from the left, looking onshore | ✅ |
| `wave_angle_high_fraction` | 0.5 | fraction | ✅ |

### Management

| Quantity | Hindcast supplies | Model expects | Converted where | Status |
|---|---|---|---|---|
| Nourishment volume | yd³ × 0.764555 ÷ n domains ÷ 500 m → **m³/m** | m³/m | `beach_dune_manager.py:829` `/100` → dam³/dam | ✅ |
| `overwash_filter` | **percent** | percent | `filter_overwash` `/100` | ✅ values in (0, 1) are rejected |
| `overwash_to_dune` | percent | percent | `/100` | ✅ |
| `rmin`, `rmax` | 0.55 / 0.95 | unitless growth rates | none | ✅ |
| bulldozer dune cap (`_artificial_maximum_dune_height`) | 4 (fixed in CASCADE, Nags Head value) | m above the berm | `filter_overwash`: `/10` → dam | ✅ since 2026-09-28 it limits only the overwash sand the manager adds to the dunes (`DUNE_CAP_APPLIES_TO`); before, it clipped the whole dune cell every year and cut natural village dunes |

---

## Outputs

| Quantity | Model native | Hindcast reports | Converted where | Status |
|---|---|---|---|---|
| Shoreline position `x_s_TS` | dam, **increases landward** | m | `cascade_pipeline/shoreline.py` `build_shoreline_matrix` ×10 | ✅ |
| `*_shoreline_matrix.npy` | — | m, raw sign (+ landward) | saved after ×10 | ✅ |
| Endpoint rate | — | m/yr, **+ = seaward** | `compute_change_rate`, sign flipped, ÷ run years | ✅ same sign as CoastSat |
| LRR | — | m/yr, **+ = seaward** | `compute_lrr`, OLS slope, sign flipped | ✅ same sign as CoastSat |
| Storm report `Rhigh_m`, `Rlow_m` | dam MHW | m MHW | `load_storm_series` ×10 | ✅ |
| Barrier height, shoreface depth diagnostics | dam | m | `DAM_TO_M = 10` | ✅ |

---

## Known quirks (not unit errors)

- **Storm run-up slope ≠ model `beta`.** Storms were built with beach slope
  0.06; `beta` is not passed to `Cascade()`, so the model uses the yaml
  default 0.04. In this setup `beta` only sets CASCADE's initial beach width
  (`int(BermEl/beta)*10`: 30 m at 0.04, 20 m at 0.06) and the shrub model,
  which is off. So it reaches the beach/dune manager only. Left as is:
  changing it changes the model inputs.
- **`drown_threshold = 0`** is commented "m MSL" in `roadway_manager.py`,
  but it is compared against MHW-relative elevations, so it is effectively
  0 m MHW (~0.26 m stricter than the comment).
- **The dune rebuild trigger** is always floored to 1.64 m MHW, so the value
  passed in does nothing.
- **`SANDBAG_ELEVATION = 0`** has no effect while sandbags are off.
- **`barrier3d.py` pops `MHW`** into `_MHW` and never uses it again. Nothing
  grid-shaped is ever shifted by MHW inside the model.
- **Barrier3D's `_SL` stays 0** (Lagrangian frame): sea-level rise lowers the
  domain instead of raising the water.
