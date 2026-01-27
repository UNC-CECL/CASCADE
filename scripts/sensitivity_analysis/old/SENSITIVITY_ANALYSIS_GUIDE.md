# CASCADE Sensitivity Analysis Guide

## Quick Start

This automated script runs CASCADE multiple times with different wave climate parameters to understand their influence on model output.

## Setup Instructions

### 1. File Placement
Place `sensitivity_analysis_wave_climate.py` in your CASCADE project directory:
```
C:\Users\hanna\PycharmProjects\CASCADE\
```

### 2. DSAS Data (Should Work Automatically!)
The script is pre-configured to use your DSAS file:
```
data/hatteras_init/shoreline_change/dsas_1978_1997_domain_means_SIMPLE.csv
```

**Expected columns:**
- `domain_id`: GIS domain numbers (1-90)
- `annual_rate_m_per_yr`: LRR shoreline change rate (m/yr)

✓ Your file has exactly this format with all 90 domains!

If this file exists, the script will automatically calculate R² and RMSE. If not found, the script still runs but only reports model statistics.

### 3. Data Extraction Method
The script now uses **your proven data extraction methods** from:
- `AllShoreline_with_DSAS_HAH_V2_LRR.py` for shoreline position extraction
- `overwash_HAH_V2.py` patterns for CASCADE object loading

This means it will handle:
- Both public and private attribute naming (`x_s_TS` vs `_x_s_TS`)
- Proper dam → meter conversion (×10)
- Correct buffer domain exclusion (only analyzes domains 15-104)
- NaN handling for domains without DSAS observations

## Running the Script

### Option 1: Run Full Analysis (All 4 Parameters)
```bash
python sensitivity_analysis_wave_climate.py
```

This will test:
- `wave_height`: 5 values (0.5, 0.75, 1.0, 1.25, 1.5 m)
- `wave_period`: 5 values (6, 7, 8, 9, 10 s) - includes your current value (7)
- `wave_asymmetry`: 7 values (0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8) - includes your current value (0.8)
- `wave_angle_high_fraction`: 5 values (0.1, 0.2, 0.3, 0.4, 0.5) - includes your current value (0.2)

**Total runs:** 22 CASCADE simulations (5 + 5 + 7 + 5)
**Estimated time:** ~44-66 minutes (if each run takes ~2-3 minutes)

**Baseline values** (what other parameters stay at when testing one):
- wave_height = 1.0 m
- wave_period = 7 s (your current value)
- wave_asymmetry = 0.8 (your current value)
- wave_angle_high_fraction = 0.2 (your current value)

### Option 2: Test Specific Parameters Only
Edit lines 74-77 to skip parameters:
```python
TEST_WAVE_HEIGHT = True              # Set to False to skip
TEST_WAVE_PERIOD = True              # Set to False to skip
TEST_WAVE_ASYMMETRY = False          # Example: skipping this one
TEST_WAVE_ANGLE_HIGH_FRACTION = True # Set to False to skip
```

### Option 3: Test Fewer Values
Modify the sensitivity ranges (lines 46-51):
```python
SENSITIVITY_RANGES = {
    'wave_height': [0.75, 1.0, 1.25],  # Only 3 values instead of 5
    'wave_period': [7, 8, 9],           # Only 3 values
    # ... etc
}
```

## What the Script Does

### 1. Automated Runs
For each parameter value:
- Creates a unique run name (e.g., `sens_wave_height_1.0_20260126_143022`)
- Runs CASCADE with that parameter while keeping all others at **your current baseline**
- Saves output to `output/sensitivity_analysis/`
- **Skips runs that already completed** (safe to restart if interrupted)

**Important:** The baseline values match your current CASCADE setup (wave_period=7, wave_asymmetry=0.8, wave_angle_high_fraction=0.2). This means:
- When testing wave_period from 6→10, all other params stay at (height=1.0, **asymmetry=0.8**, **angle_fraction=0.2**)
- You can see the isolated effect of changing each parameter from your actual starting point
- Your current configuration will be tested (e.g., when wave_period=7 is tested)

### 2. Metrics Calculated
For each run:
- **Mean change**: Average shoreline change rate
- **Median change**: Median shoreline change rate
- **Std dev**: Spatial variability
- **Min/max/range**: Extent of variability
- **R²**: Correlation with observations (if DSAS data available)
- **RMSE**: Root mean square error (if DSAS data available)
- **MAE**: Mean absolute error (if DSAS data available)
- **Bias**: Systematic over/under-prediction (if DSAS data available)

### 3. Outputs Generated

**Results CSV** (`sensitivity_results_TIMESTAMP.csv`):
- One row per run
- All metrics in spreadsheet format
- Easy to import into Excel or further analysis

**Plots** (in `plots_TIMESTAMP/` folder):
1. `sensitivity_mean_change.png` - How mean change varies with each parameter
2. `sensitivity_spatial_variability.png` - How std dev varies
3. `sensitivity_r_squared.png` - Model fit for each parameter (if DSAS data)
4. `tornado_diagram.png` - Relative sensitivity comparison
5. `summary_table.png` - Summary statistics table

## Interpreting Results

### Key Questions to Answer:

1. **Which parameter has the strongest influence?**
   - Look at tornado diagram
   - Parameter with largest bar = most influential

2. **Should you change from your current values?**
   - Your current setup: wave_period=7, wave_asymmetry=0.8, wave_angle_high_fraction=0.2
   - Look at R² plots: do suggested values (period=8-9, asymmetry=0.55-0.6, angle=0.3-0.4) improve fit?
   - If yes → update to better values
   - If no → keep current values (they're already optimal)

3. **Can wave parameters alone create enough variability?**
   - Look at range values
   - Compare to observed range (your data shows ~-5 to +5 m/yr, range of ~10 m/yr)
   - If modeled range << observed range → background erosion is essential

4. **What's the best parameter combination?**
   - Look at CSV file
   - Sort by R² (if available) or RMSE
   - Highest R² / lowest RMSE = best fit

### Expected Findings (Based on Your Advisor's Approach):

You'll likely find that:
- Wave parameters change the **magnitude** slightly
- Wave parameters create **some spatial variability** from transport gradients
- But **cannot reproduce the full 10 m/yr range** you see in observations
- Moving from current values (7, 0.8, 0.2) to suggested values (8-9, 0.55-0.6, 0.3-0.4) may improve fit slightly
- **But even with optimal parameters, background erosion is still essential**

This justifies calculating background erosion as a residual!

## Customization Options

### Change Baseline Parameters
Edit lines 40-45:
```python
BASELINE_PARAMS = {
    'wave_height': 1.0,           # Change these if you want
    'wave_period': 8,             # different baseline values
    'wave_asymmetry': 0.55,
    'wave_angle_high_fraction': 0.3,
}
```

### Add More Test Values
Edit lines 46-51 to add more granular testing:
```python
'wave_height': [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5],
```

### Change Number of CPU Cores
Edit line 68:
```python
NUM_CORES = 4  # Increase if you have more cores available
```

## Troubleshooting

### "DSAS data file not found"
- This is just a warning
- Script will run without R² and RMSE calculations
- To fix: provide the path to your DSAS CSV file

### "CASCADE output not found"
- Make sure CASCADE completed successfully
- Check the output directory for .npz files
- Look for error messages during CASCADE run

### Script Interrupted
- Just restart it!
- Already-completed runs are automatically skipped
- Progress will resume where it left off

### Out of Memory
- Reduce `NUM_CORES` 
- Test fewer parameters at once
- Close other applications

### Wrong Results Format
The script now uses your proven extraction methods from `AllShoreline_with_DSAS_HAH_V2_LRR.py`, so it should work correctly with your CASCADE setup. It handles:
- `x_s_TS` or `_x_s_TS` attributes (whichever your CASCADE version uses)
- Automatic dam → meter conversion
- Buffer domain exclusion

If you still get errors:
- Check that CASCADE completed successfully (look for CASCADE.npz file)
- Verify the .npz file contains the cascade object
- Make sure cascade.barrier3d list is populated

## Next Steps After Running

1. **Review the CSV file** - Open in Excel, sort by metrics
2. **Examine the plots** - Identify most influential parameters
3. **Select best parameter values** - Based on R² or physical reasoning
4. **Calculate background erosion residual**:
   ```python
   background_erosion = observed_rate - modeled_rate_with_best_params
   ```
5. **Re-run CASCADE** with best wave params + calculated background erosion
6. **Document findings** for your thesis/proposal

## For Your Thesis Methods Section

You can write something like:

> "A sensitivity analysis was conducted to identify optimal wave climate parameters for the Hatteras Island CASCADE implementation. Each of four wave parameters (wave height, wave period, wave asymmetry, and wave angle high fraction) was systematically varied across physically realistic ranges while holding all other parameters constant. Model performance was evaluated using spatial correlation (R²), root mean square error (RMSE), and spatial variability metrics against DSAS-derived shoreline change rates.
>
> Results showed that [parameter X] had the strongest influence on model output (Figure X), with [describe results]. However, wave climate parameters alone could not reproduce the full spatial variability observed in the data (modeled range: ±X m/yr vs. observed range: ±5 m/yr), indicating that background erosion processes not captured by storm events and wave-driven transport are essential for explaining Hatteras Island evolution. This finding justified the calculation of background erosion as a residual quantity: the difference between observed change rates and modeled change from wave and storm processes."

## Questions?

If you run into issues or need to modify the script, let me know! I can help adjust:
- Parameter ranges
- Metrics calculated
- Plot styles
- CASCADE configuration
- Anything else

---
**Script Version:** 1.0  
**Last Updated:** January 26, 2026  
**Author:** Hannah Henry (with Claude assist)
