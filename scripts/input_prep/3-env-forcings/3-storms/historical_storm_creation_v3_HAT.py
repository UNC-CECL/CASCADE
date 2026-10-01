"""
Build the historical storm series CASCADE reads for one hindcast window, from Duck water levels and WIS waves.

    python scripts/input_prep/3-env-forcings/3-storms/historical_storm_creation_v3_HAT.py --start-year 1996 --end-year 2010

Total water level = gauge + Stockdon R2%; events above the berm, split,
capped and trimmed; writes the storm .npy and its summary CSV. Its four
functions are also read by hatteras_ms/experiments/HAT_storm_max_duration.py. Details: scripts/input_prep/3-env-forcings/README.md.

Adapted from: from_lexi/historical_storm_creation_v3.ipynb, by Lexi
Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

# load necessary packages
import numpy as np
import pandas as pd
import os
from datetime import datetime, timedelta


# User inputs

# The two source paths are anchored on this file now (history in README)
import argparse as _argparse
from pathlib import Path as _Path

# Repo root, found by searching upward
_PATH_REPO = next(_p for _p in _Path(__file__).resolve().parents
                  if (_p / "pyproject.toml").exists())

_REPO = next(_p for _p in _Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
_ap = _argparse.ArgumentParser(
    description="CASCADE storm series for one hindcast window")
_ap.add_argument("--start-year", type=int, required=True,
                 help="period start, the year that becomes model step 1")
_ap.add_argument("--end-year", type=int, required=True,
                 help="period end. The model loop runs start..end-1, so the "
                      "last storm year the run can spend is end-1")
# THE LONG-EVENT RULE (2026-09-28, Hannah adopted "trim24")
_ap.add_argument("--max-duration", type=int, default=72,
                 help="hours; the drop limit, or the trim length with --long-events trim")
_ap.add_argument("--long-events", choices=("drop", "trim"), default="drop",
                 help="drop events longer than --max-duration (v3_<N>), or trim them (v3_trim<N>)")
# THE SPLIT RULE (2026-09-29, Hannah adopted "split12")
_ap.add_argument("--split-gap", type=int, default=None,
                 help="hours below the berm that split a grouped event (v3_split<H>_...)")
_ap.add_argument("--save-dir", default=None,
                 help="write here instead of the window's hindcast_storms folder (for checks)")
_args = _ap.parse_args()

# --- CONFIG ------------------------------------------------------------------
START_YEAR = _args.start_year
END_YEAR = _args.end_year
PERIOD_TAG = "{0}_{1}".format(START_YEAR, END_YEAR)

start_time = '{0}-01-01 00:00:00'.format(START_YEAR)  # date to start the storms
end_time = '{0}-12-31 23:00:00'.format(END_YEAR)      # date to end the storms
import sys as _envsys
from pathlib import Path as _EnvP
_envsys.path.insert(0, str(next(_q for _q in _EnvP(__file__).resolve().parents
                                if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_env_forcings as _env  # noqa: E402
water_levels_file = str(_env.DUCK_GAUGE_FILE)   # NOAA gauge
wis_file = str(_env.WIS_FILE)                     # WIS gauge
t_name_water = "t"       # name of the column that contains the datetimes in the water levels file
water_name = "v"         # name of the column that contains the water levels [m NAVD88] in the water levels file
t_name_wis = "time"  # name of the column that contains the datetimes in the WIS file
waveHs_name = "waveHs"   # name of the column that contains the significant wave heights [m] in the WIS file
waveTp_name = "waveTp"   # name of the column that contains the wave periods [s] in the WIS file
beach_slope = 0.06       # average beach slope for calculating run-up
berm_elevation = 1.7    # average berm elevation [m NAVD88]
weather_grouping = 24    # if storms occur within the specified limit, they are assumed part of same weather system and are grouped into one event [hrs]
MHW = 0.36              # conversion from NAVD88 to MHW [m]: 0 m NAVD88 = X m MHW
min_storm_dur = 8        # minimum duration that is considered a storm event [hrs]
max_storm_dur = _args.max_duration   # hours: the drop limit, or the trim length (see --long-events)
long_events = _args.long_events      # "drop" (v3_<N>) or "trim" (v3_trim<N>)
split_gap = _args.split_gap          # None, or hours below the berm that split an event
save_dfs = True         # determine whether to save the dataframes as csv and npy files
save_dir = _args.save_dir or str(_env.storm_window_dir(START_YEAR, END_YEAR))  # derived from the window
save_name = "{0}_storms_v3_{1}{2}{3}".format(  # the rules are in the name: v3_72, v3_trim24, v3_split12_trim24
    PERIOD_TAG, "split{0}_".format(split_gap) if split_gap else "",
    "trim" if long_events == "trim" else "", max_storm_dur)
# -----------------------------------------------------------------------------
_Path(save_dir).mkdir(parents=True, exist_ok=True)


# Do not edit code blocks below - all inputs are updated automatically


# Step 1: Load data and combine into a single dataframe

# this function finds all values that are NaN and groups them by consecutive time periods NOTE
def find_time_gaps(df, col_name):
    
    # check for NaN value and create a datetime list of all corresponding NaNs
    nan_index_numbers = np.where(df[col_name].isnull())[0]  # returns numeric index
    time_list_nans = df.index[nan_index_numbers]  # returns the datetime index

    # A step of more than 1 hour between NaN times ends one gap and starts the next
    datetime_nans = []
    start = time_list_nans[0]
    datetime_nans.append(start)  # add the first value to the list
    for t in range(len(time_list_nans)-1):
        # if more than 1 hr between this time (t) and the next time in the list, we have transition from one consecutive time period of nans to another
        if time_list_nans[t+1]-time_list_nans[t]!=timedelta(seconds=3600):  
            datetime_nans.append(time_list_nans[t])  # this is the last time step in a consecutive period
            datetime_nans.append(time_list_nans[t+1])  # this is the first time step in the next group of NaNs
    stop = time_list_nans[-1]
    datetime_nans.append(stop)  # add last value to list
    
    # Start and end datetimes: starts at odd indices, ends at even
    start_dt = []
    end_dt = []
    n_days = []
    for t in range(len(datetime_nans)):
        if t%2 == 0:
            start_dt.append(datetime_nans[t])  # just to make sure it always outputs in the Timestamp format
        else:
            end_dt.append(datetime_nans[t])
            
    # create comparison: missing days and hours to stay separate (e.g. 7 days and 12 hours):
    missing_days = [(end_dt[dt] - start_dt[dt]).days for dt in range(len(start_dt))]
    missing_hours = [(end_dt[dt] - start_dt[dt]).seconds/3600 + 1 for dt in range(len(start_dt))]
    missing_data_ranges = pd.DataFrame(columns=['start date', 'end date', "days", "hours"])
    missing_data_ranges["start date"] = start_dt
    missing_data_ranges["end date"] = end_dt
    missing_data_ranges["days"] = missing_days
    missing_data_ranges["hours"] = missing_hours
    
    print("{0} column is missing the following date ranges".format(col_name))
    print(missing_data_ranges)


# Load the records into one dataframe over start..stop, with rows holding NaNs dropped
def load_data(
    start_time,
    end_time,
    water_levels_file, 
    wis_file, 
    t_name_water="t", 
    water_name="v", 
    t_name_wis="time",
    waveHs_name="waveHs", 
    waveTp_name="waveTp"
):
    
    # initialize index with continuous datetimes
    dt_list = []
    start_dt = datetime.strptime(start_time, "%Y-%m-%d %H:%M:%S")
    end_dt = datetime.strptime(end_time, "%Y-%m-%d %H:%M:%S")
    # total hours
    hours = int(((end_dt - start_dt).seconds/3600))  # number of hours
    days = ((end_dt - start_dt).days)  # number of days
    total_hours = hours + days*24
    # create an hourly datetime list
    for n in range(total_hours + 1):
        dt_list.append(start_dt + timedelta(hours=n))

    # load water levels
    df_water = pd.read_csv(water_levels_file, index_col=t_name_water)
    df_water.index = pd.to_datetime(df_water.index)
    water_levels = df_water[water_name]
    
    # load WIS values
    df_wis = pd.read_csv(wis_file, index_col=t_name_wis)
    df_wis.index = pd.to_datetime(df_wis.index)
    waveHs = df_wis[waveHs_name]
    waveTp = df_wis[waveTp_name]

    # Create a merged dataframe   
    df_merged = pd.DataFrame(index=dt_list)  # the datatimes of both dfs should be the same
    df_merged["water_level"] = water_levels
    df_merged["Hs"] = waveHs
    df_merged["Tp"] = waveTp
    
    # Check for missing data water levels
    nan_index_numbers = np.where(df_merged["water_level"].isnull())[0]  # returns numeric index
    if len(nan_index_numbers) > 0:
        find_time_gaps(df_merged, "water_level")
    # wave height
    nan_index_numbers = np.where(df_merged["Hs"].isnull())[0]  # returns numeric index
    if len(nan_index_numbers) > 0:
        find_time_gaps(df_merged, "Hs")
    # wave period
    nan_index_numbers = np.where(df_merged["Tp"].isnull())[0]  # returns numeric index
    if len(nan_index_numbers) > 0:
        find_time_gaps(df_merged, "Tp")
    
    # Remove rows with missing data
    df_merged = df_merged.dropna()
    
    return df_merged


df_merged = load_data(
    start_time=start_time,
    end_time=end_time,
    water_levels_file=water_levels_file, 
    wis_file=wis_file,
    t_name_water=t_name_water, 
    water_name=water_name,
    t_name_wis=t_name_wis,
    waveHs_name=waveHs_name, 
    waveTp_name=waveTp_name
)

# review merged dataframe
print(df_merged)


# Step 2: Calculate the runup to use for the Total Water Level (twl)

# Calculate R2% (2% exceedance runup) using Stockdon et al
def calculate_r2_percent(Hs, Tp, slope):
    g = 9.81
    # Deep water wavelength
    L0 = (g * Tp**2) / (2 * np.pi)
    
    # Setup component
    setup = 0.35 * slope * np.sqrt(Hs * L0)
    
    # Swash components
    S_incident = 0.75 * slope * np.sqrt(Hs * L0)
    S_infragravity = 0.06 * np.sqrt(Hs * L0)
    S_total = np.sqrt(S_incident**2 + S_infragravity**2)
    
    # R2% runup
    R2 = 1.1 * (setup + S_total/2)
    
    return R2


# calculate R2 based on the beach slope
r2_values = calculate_r2_percent(df_merged['Hs'], df_merged['Tp'], beach_slope)

# add runup and TWL to the combined df
df_merged['R2'] = r2_values
df_merged['TWL'] = df_merged['water_level'] + df_merged['R2']

# review merged dataframe
print(df_merged)


# Step 3: Create the CASCADE storms based on berm elevation and twl

# df_merged
def create_storms(
    df_merged, 
    berm_elevation, 
    weather_grouping=24, 
    MHW=0.421,
    min_storm_dur=8,
    max_storm_dur=240,
    save_dfs=True,
    save_dir="",
    save_name="",
    window_start_year=None,
    long_events="drop",
    split_gap=None,
):
    
    # Step 1 identify storms

    # initialize new dataframe
    df = pd.DataFrame({
        "Time": pd.to_datetime(df_merged.index),
        "TWL": pd.to_numeric(df_merged["TWL"], errors="coerce"),
        "Tp": pd.to_numeric(df_merged["Tp"], errors="coerce")
    })
    df = df.sort_values("Time").reset_index(drop=True)

    # determine time steps that are storms (TWL > berm elevation)
    df["AboveBerm"] = df["TWL"] > berm_elevation

    # assign storm start based on consecutive times when TWL > berm elevation
    df["StormStart"] = df["AboveBerm"] & (~df["AboveBerm"].shift(1, fill_value=False)) 

    # adjust StormStart for continuous weather system
    temp_df = df[df["AboveBerm"]]
    index_vals = temp_df.index
    for i in range(1, len(index_vals)):  # skip the first row
        current_row = temp_df.loc[index_vals[i]]
        prev_row = temp_df.loc[index_vals[i-1]]
        if current_row.StormStart==True:  # if storm start is set to true, check if it is within limit hours of the previous False
            current_time = current_row.Time
            prev_row_time = prev_row.Time
            time_delta_hrs = (current_time - prev_row_time).days*24 + (current_time - prev_row_time).seconds/3600
            if time_delta_hrs < weather_grouping:  # if it is within limit, we want to make it part of the previous storm instead of a new storm
                df.loc[index_vals[i], "StormStart"] = False  # index values should be the same in the main df

    df["StormID"] = df["StormStart"].cumsum()
    df.loc[~df["AboveBerm"], "StormID"] = np.nan
    
    
    # Step 2 create storms

    storms = []
    storm_groups = df.dropna(subset=["StormID"]).groupby("StormID")

    # split_gap: cut events where the water stays below the berm that long; short pieces rejoin
    if split_gap:
        def _split(group):
            group = group.sort_values("Time")
            gaps = group["Time"].diff().dt.total_seconds().div(3600).fillna(0).values
            cuts = [i for i, g in enumerate(gaps) if g >= split_gap]
            bounds = [0] + cuts + [len(group)]
            pieces = [list(range(a, b)) for a, b in zip(bounds[:-1], bounds[1:])]
            i = 0
            while len(pieces) > 1 and i < len(pieces):
                if len(pieces[i]) < min_storm_dur:
                    j = i - 1 if i > 0 else i + 1
                    pieces[j] = sorted(pieces[j] + pieces[i])
                    del pieces[i]
                    i = 0
                    continue
                i += 1
            return [group.iloc[piece] for piece in pieces]

        storm_groups = [((sid, k), piece) for sid, group in storm_groups
                        for k, piece in enumerate(_split(group))]

    for sid, group in storm_groups:

        group = group.sort_values("Time").copy()
        start_time = group["Time"].iloc[0]
        end_time   = group["Time"].iloc[-1]
        event_year = start_time.year          # the year the EVENT starts, kept when trimmed

        # Duration (hours)
        duration = len(group)  # each row is 1 hour where an exceedance occured

        # 'trim': keep a long event, cut to max_storm_dur hours around its peak
        trimmed_from = 0
        if long_events == "trim" and duration > max_storm_dur:
            k = int(np.argmax(group["TWL"].values))
            lo = min(max(0, k - max_storm_dur // 2), duration - max_storm_dur)
            trimmed_from = duration
            group = group.iloc[lo:lo + max_storm_dur]
            start_time = group["Time"].iloc[0]
            end_time = group["Time"].iloc[-1]
            duration = len(group)

        # Rhigh and Rlow from TWL during the storm (in m NAVD88)
        rhigh = group["TWL"].max()
        rlow  = group["TWL"].min()
        # convert to decameters MHW
        rhigh = (rhigh - MHW) / 10
        rlow = (rlow - MHW) / 10

        # Period from Tp at peak TWL
        peak_idx = group["TWL"].idxmax()
        period = group.loc[peak_idx, "Tp"]
        

        # only add storms > specified duration (Magliocca et al., 2011) but less than maximum 
        if duration >= min_storm_dur and duration <= max_storm_dur:
            storms.append({
                "calendar_year": event_year,
                "StartTime": start_time,
                "EndTime": end_time,
                "Rhigh": rhigh,
                "Rlow": rlow,
                "period": period,
                "duration": duration,
                **({"trimmed_from": trimmed_from} if long_events == "trim" else {}),
            })

    storms_df = pd.DataFrame(storms)

    # Convert year → time (CASCADE)

    # Model step 1 is the window's first year, not the first stormy year
    if not storms_df.empty:
        first_storm_year = int(storms_df["calendar_year"].min())
        anchor = first_storm_year if window_start_year is None             else int(window_start_year)
        if anchor != first_storm_year:
            print("NOTE: no storm in {0}; anchoring model step 1 on the window "
                  "start anyway, so step 1 is {0} and the first storm falls at "
                  "step {1}.".format(anchor, first_storm_year - anchor + 1))
        storms_df["time"] = storms_df["calendar_year"] - anchor + 1
    else:
        storms_df["time"] = pd.Series(dtype=float)

    # Final format

    cascade_df = storms_df[["time", "Rhigh", "Rlow", "period", "duration"]].copy()
    cascade_df = cascade_df.reset_index(drop=True)
    casc_input_array = cascade_df.to_numpy()
    total_storms = len(cascade_df)
    avg_duration = round(np.mean(cascade_df.duration.values),0)

    # Save

    if save_dir == "":
        save_dir = os.getcwd()
    
    if save_name == "":    
        save_storm = os.path.join(save_dir, "storm_summary.csv")
        save_casc = os.path.join(save_dir, "cascade_storms.csv")
        save_npy = os.path.join(save_dir, "cascade_storms.npy")
    else:
        save_storm = os.path.join(save_dir, "{0}_summary.csv".format(save_name))
        save_casc = os.path.join(save_dir, "{0}.csv".format(save_name))
        save_npy = os.path.join(save_dir, "{0}.npy".format(save_name))
    
    if save_dfs:
        storms_df.to_csv(save_storm, index=False)
        cascade_df.to_csv(save_casc, index=False)
        np.save(save_npy, casc_input_array)
        print("dataframes saved to {0}".format(save_dir))
    else:
        print("Preview of storms below. Set save_dfs to True to save full versions.")
        print("Total storms: {0}".format(total_storms))
        print("Average duration: {0} hours".format(avg_duration))
        print("\nStorm summary:")
        print(storms_df.head())

        print("\nCASCADE CSV format (will not have index column):")
        print(cascade_df.head(5))

        print("\nCASCADE NPY format:")
        print(casc_input_array[0:5, :])


create_storms(
    df_merged=df_merged,
    berm_elevation=berm_elevation, 
    weather_grouping=weather_grouping, 
    MHW=MHW,
    min_storm_dur=min_storm_dur,
    max_storm_dur=max_storm_dur,
    save_dfs=save_dfs,
    save_dir=save_dir,
    save_name=save_name,
    window_start_year=START_YEAR,
    long_events=long_events,
    split_gap=split_gap,
)
