# Domain Mapping Reference - Hatteras Island CASCADE Model

## Overview

Your CASCADE model uses **two different numbering systems** that need to be mapped correctly:

1. **GIS/DSAS domain numbers** (from ArcGIS shoreline analysis)
2. **CASCADE/Python domain indices** (model array indices)

## Domain Structure

### GIS/DSAS Domain Numbers
```
Domain 1-90:   YOUR Hatteras study area (South → North)
               ↓ These are the domains you're modeling
Domain 91-120: Your collaborator's study area
               ↓ NOT part of your CASCADE model
               ↓ Should be filtered out before processing
```

### CASCADE/Python Domain Indices
```
Index 0-14:     Left buffer domains (15 total)
                ↓ Fake domains to prevent boundary effects
                ↓ Always set to 0.0 background rate

Index 15-104:   Real island domains (90 total)
                ↓ Your actual Hatteras Island
                ↓ Maps to GIS domains 1-90

Index 105-119:  Right buffer domains (15 total)
                ↓ Fake domains to prevent boundary effects
                ↓ Always set to 0.0 background rate
```

## The Mapping

```
GIS Space          CASCADE Space         Geography
─────────────────────────────────────────────────────────
Domain 1      →    Index 15          South end (Cape Hatteras area)
Domain 2      →    Index 16
Domain 3      →    Index 17
...
Domain 45     →    Index 59          Mid-island
...
Domain 89     →    Index 103
Domain 90     →    Index 104         North end (toward Oregon Inlet)

Domain 91-120 →    [NOT USED]        Collaborator's area
```

## Formula

```python
# GIS to CASCADE:
cascade_index = gis_domain + 14

# CASCADE to GIS:
gis_domain = cascade_index - 14
```

**Example:**
- GIS Domain 1 → CASCADE Index 15 (because 1 + 14 = 15)
- GIS Domain 90 → CASCADE Index 104 (because 90 + 14 = 104)

## Geographic Direction

**South to North** (same in both systems):
- Lower numbers = South (Cape Hatteras)
- Higher numbers = North (Oregon Inlet direction)

```
        North ↑
              |
     Oregon Inlet area
              |
    GIS 90 / CASCADE 104
              |
         [your study]
              |
    GIS 45 / CASCADE 59
              |
         [your study]
              |
    GIS 1 / CASCADE 15
              |
    Cape Hatteras area
              |
        South ↓
```

## File Processing Rules

### When Cleaning DSAS Data (`dsas_clean_CORRECTED.py`)

**ALWAYS filter to domains 1-90:**

```python
# Remove collaborator's domains (91-120)
df = df[df['domain_id'] <= 90].copy()
```

**Result:** Output file has exactly 90 records (GIS domains 1-90)

### When Generating Background Rates (`generate_background_rates.py`)

**Configuration should be:**

```python
TOTAL_DOMAINS = 120      # 15 + 90 + 15
NUM_REAL_DOMAINS = 90    # Your actual island
NUM_LEFT_BUFFER = 15     # South buffer
NUM_RIGHT_BUFFER = 15    # North buffer
```

**Result:** Output array has 120 values:
- Indices 0-14: 0.0 (left buffer)
- Indices 15-104: Calculated rates (real island)
- Indices 105-119: 0.0 (right buffer)

### When Comparing Model to Observations (`AllShoreline_with_DSAS_HAH_V2_LRR.py`)

**Configuration should be:**

```python
LEFT_BUFFER_END = 14        # Last index of left buffer
RIGHT_BUFFER_START = 105    # First index of right buffer

# When loading DSAS data:
gis_domains = dsas['domain_id'].to_numpy()  # Values 1-90
obs_domains = gis_domains + 14              # Convert to CASCADE indices 15-104
```

**Result:** Observed data at CASCADE indices 15-104 for comparison

## Common Mistakes to Avoid

❌ **DON'T:**
- Include domains 91-120 in your cleaned DSAS file
- Set NUM_REAL_DOMAINS = 120
- Compare GIS domain numbers directly to CASCADE indices

✓ **DO:**
- Filter DSAS data to domains 1-90 before processing
- Use NUM_REAL_DOMAINS = 90
- Add 14 when converting GIS domains to CASCADE indices

## Verification Checklist

Before running your model, verify:

- [ ] DSAS cleaned file has exactly 90 records (domains 1-90)
- [ ] Background rates array has 120 values
- [ ] Indices 0-14 are all 0.0 (left buffer)
- [ ] Indices 15-104 have calculated rates (real island)
- [ ] Indices 105-119 are all 0.0 (right buffer)
- [ ] Domain 1 characteristics match your south end geography
- [ ] Domain 90 characteristics match your north end geography

## Quick Reference Table

| Description | GIS/DSAS | CASCADE | Count |
|-------------|----------|---------|-------|
| Study area start | Domain 1 | Index 15 | - |
| Study area end | Domain 90 | Index 104 | - |
| Total study area | Domains 1-90 | Indices 15-104 | 90 |
| Left buffer | N/A | Indices 0-14 | 15 |
| Right buffer | N/A | Indices 105-119 | 15 |
| Total CASCADE domains | N/A | Indices 0-119 | 120 |
| Collaborator area (unused) | Domains 91-120 | N/A | 30 |

## Example Output Verification

When you run `generate_background_rates.py` successfully, you should see:

```
LOADED DSAS DATA
Records: 90                    ← Should be 90, not 120!
Domain range: 1 to 90          ← Should stop at 90!

FINAL ARRAY STRUCTURE
Total domains: 120
  Left buffer (0.0): domains 0-14
  Real island: domains 15-104
  Right buffer (0.0): domains 105-119

SUMMARY STATISTICS
Real island domains (90 total):  ← Should be 90!
  Mean rate: 0.000 m/yr
  ...
```

If you see "Records: 120" or "Domain range: 1 to 120", something is wrong!
