# 4-split_windows — what the two halves look like at one transect

Asked for by colleagues (2026-10-02): the straight OLS lines at individual
transects, with year on x and shoreline position on y. Each panel shows the
1996–2024 rate and the two windows either side of a cutoff, so you can see
whether the halves the model is run on agree with the long-term rate.

```
split_windows_2010.png       eight transects, cut at 2010 (the model's two legs)
split_windows_cutoffs.png    the same eight, one row each (north at top), cut at 2005 / 2010 / 2015 / 2020,
                             with a locator map numbering each row 1 (Buxton) to 8 (Rodanthe)
split_windows_picks.csv      the eight: group, rank in the group, the three rates, the 2021 step
supporting/                  PDFs, CAPTIONS.md, split_windows_transects.csv (every transect scored)
interactive/                 split_windows_explorer.html + split_windows_data.json:
                             all 906 transects, a map picker and a cutoff slider
```

Interactive page (private until shared from its Share menu):
https://claude.ai/artifact/7Bhfzp2eeSU1v1Dx3s5vGt

## What every panel shows

- Grey dots: every CoastSat position, 1996–2024. Open circles: the annual
  median, a guide for the eye only. No line is fitted to it.
- Purple: OLS over 1996–2024. Red: the first window, 1996 to the cutoff.
  Blue: the second, the cutoff to 2024. The cutoff year belongs to **both**
  windows (calendar years, inclusive), as in every other window product here.
- Each line is fitted to the raw positions and drawn over the dates its window
  saw. The rates are the `1-rate_profiles/` fits, so they are not refitted.
  The page refits in the browser with the same arithmetic (time = days since
  1 Jan 1996 / 365.25), and it matches the stored fits to 0.0002 m/yr.
- Position is CoastSat chainage on an arbitrary origin per transect, so only
  slopes compare between panels.

## How the eight were picked (Hannah, 2026-10-02)

Not evenly spaced (`3-settling_window/.../a-eight_sites/` already does that).
Instead, two transects from each of four behaviour groups, so every kind of
behaviour appears on the page. The groups use the 2010 cutoff:

| group | in the group | ranked by |
|---|---|---|
| halves agree | both halves within ±0.5 m/yr of 1996–2024 | the closest agreement |
| halves disagree | not a sign flip | the largest gap between the halves |
| sign flip | halves of opposite sign, both at least 0.5 m/yr | the largest gap between the halves |
| 2021 step | domains not nourished (`step_2021_by_domain_prefill.csv`) | 2021 median minus 2019 median, largest |

The groups are filled in that order. Each takes the top-ranked transect that
is at least two domains from every earlier pick, and GIS 1 and 90 (the end
buffers) are never picked. Because of that rule, some picks are not rank 1;
the rank is in the CSV. "Disagree" excludes the sign flips, otherwise both
groups would show the same transects.

## What it shows

- **Agree** (GIS 56, 58): flat records, where all three lines lie on top of
  each other at any cutoff.
- **Disagree** (GIS 81, 66): a trend that changes partway through. At GIS 81,
  −6.9 m/yr to 2010 and then flat, so the long-term −2.5 describes neither
  half.
- **Sign flip** (GIS 36, 11): erosion to about 2010, then accretion. The
  long-term rate is near zero and matches neither window.
- **2021 step** (GIS 33, 18): +44 to +50 m between 2019 and 2021. This roughly
  doubles the second-window rate while the first window agrees with 1996–2024.
- In `split_windows_cutoffs.png`, a cut at 2020 leaves a five-year second
  window that the 2021 step dominates, and its rate swings by several m/yr.

## Producer

`scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/coastsat_split_windows.py`
reads `1-rate_profiles/*/window_profiles_transects.csv`, so run
`coastsat_window_profiles.py` first after any change to the fits. Running the
producer rewrites the data file. To update the page, republish
`interactive/split_windows_explorer.html` to the same URL with the data file
alongside it.
