# CASCADE Storm Viewer - Quick Guide

## What It Does

A clean, standardized viewer for CASCADE storm files. Does everything your old scripts did, but in one consistent place.

## Quick Start

1. **Set your paths** in the config section:
```python
STORM_PATH = Path(r"C:\path\to\storms_1978_1997.npy")
START_YEAR = 1978
```

2. **Run it**:
```bash
python storm_viewer_simple.py
```

That's it!

## Configuration

```python
# Required
STORM_PATH = Path(r"...")    # Your storm .npy file
START_YEAR = 1978            # First year of hindcast

# Optional
SHOW_PLOTS = True           # Display plots on screen
SAVE_PLOTS = False          # Save plots to disk  
OUTPUT_DIR = Path(r"...")   # Where to save (if SAVE_PLOTS=True)

# Year verification (for checking against historical records)
CHECK_YEAR = None           # Set to a year (e.g., 1985) to see all storms that year
```

## Verifying Storms Against Historical Records

Need to check if your storm file matches historical records for a specific year?

### Option 1: Set CHECK_YEAR in config
```python
CHECK_YEAR = 1985  # Will show all storms in 1985 when you run the script
```

### Option 2: Check multiple years interactively
```python
# After running the script, uncomment these lines at the bottom:
if __name__ == "__main__":
    df = main()
    
    # Check specific years
    check_year(df, 1985)  # Hurricane Gloria
    check_year(df, 1996)  # Hurricane Fran
    check_year(df, 1999)  # Hurricanes Dennis, Floyd, Irene
```

### Example Output:
```
================================================================================
STORMS IN 1985
================================================================================

Total storms: 3

Rhigh range: 1.234 to 2.456
Average Rhigh: 1.789
Total storm duration: 72.5 hours

--------------------------------------------------------------------------------
ALL 3 STORMS IN 1985:
--------------------------------------------------------------------------------
   Calendar_Year  Rhigh   Rlow  Wave Period  Duration
1           1985  2.456  1.234         12.3      36.0
2           1985  1.890  0.987         11.5      24.5
3           1985  1.234  0.678         10.2      12.0

Most intense storm: Rhigh=2.456, Duration=36.0 hrs
================================================================================
```

## Configuration

## What You Get

### Console Output:
- File info (shape, columns)
- Summary statistics (min, max, mean, percentiles)
- First 10 events
- Top 5 most intense storms
- Storms per year

### Plots:
1. **Four-panel distribution plot**:
   - Rhigh histogram
   - Wave period histogram
   - Duration histogram
   - Rlow vs Rhigh scatter

2. **Time series plot**:
   - Individual storm intensities (gray dots)
   - Yearly average intensity (purple line)

## Examples

### View Hatteras storms
```python
STORM_PATH = Path(r"C:\...\storms_1978_1997.npy")
START_YEAR = 1978
SHOW_PLOTS = True
SAVE_PLOTS = False
CHECK_YEAR = None
```

### Verify 1996 storms against historical records
```python
STORM_PATH = Path(r"C:\...\storms_1978_1997.npy")
START_YEAR = 1978
SHOW_PLOTS = False
SAVE_PLOTS = False
CHECK_YEAR = 1996  # Hurricane Fran year
```

### Create figures for proposal
```python
STORM_PATH = Path(r"C:\...\storms_1978_1997.npy")
START_YEAR = 1978
SHOW_PLOTS = False
SAVE_PLOTS = True
OUTPUT_DIR = Path(r"C:\...\output\figures")
CHECK_YEAR = None
```

### Compare to Ocracoke storms
```python
# Run once for Hatteras
STORM_PATH = Path(r"C:\...\storms_1978_1997.npy")
START_YEAR = 1978
SAVE_PLOTS = True
OUTPUT_DIR = Path(r"C:\...\figures\hatteras")

# Then run again for Ocracoke
STORM_PATH = Path(r"C:\...\OCR_storms.npy")
START_YEAR = 1974
SAVE_PLOTS = True
OUTPUT_DIR = Path(r"C:\...\figures\ocracoke")
```

## Notes

- **Column structure**: Assumes CASCADE standard format: `[Year_Index, Rhigh, Rlow, Wave Period, Duration]`
- **Calendar year**: Automatically converts `Year_Index` (0,1,2...) to `Calendar_Year` (1978,1979,1980...)
- **File format**: Works with `.npy` files (standard CASCADE storm format)

## What Changed From Your Old Scripts

This replaces:
- ✅ `check_storm_npy.py` (v3)
- ✅ `check_storms_v2.py`

Keep separate:
- ⚠️ `storm_diagnostic_HAH_V1.py` (different purpose - analyzes CASCADE model output)

## Relationship to Diagnostic Dashboard

**This script**: Inspects storm *input* files  
**Diagnostic dashboard**: Analyzes CASCADE model *output*

They're complementary tools for different stages of your workflow:
1. Use **storm viewer** to verify storm file before running CASCADE
2. Use **diagnostic dashboard** to analyze CASCADE results after running

## Tips

- Check the console output first - it tells you everything about the file
- For proposal figures, use `SAVE_PLOTS = True` and `SHOW_PLOTS = False`
- Keep `START_YEAR` correct or your time series will be wrong!
