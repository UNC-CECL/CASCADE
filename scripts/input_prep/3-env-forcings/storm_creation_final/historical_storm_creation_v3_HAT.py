"""
Creating a historical storm series for CASCADE
================================================
CASCADE requires an input storm series consisting of the model year the storm
occurs, maximum run-up, minimum run-up, wave period, and storm duration.
Synthetic storms can be created using the multi-variate sea storm model of
Wahl et al. (2016), but this script converts recorded water level and wave
data into a CASCADE storm input.

Input files should be csv files with headers. Any missing data will be filled
with nans and dropped. Water levels can be pulled from a NOAA tide gauge
(https://tidesandcurrents.noaa.gov/) and wave information can be obtained
from a WIS Station (https://wisportal.erdc.dren.mil/).

Input files include:
- water levels file: should have datetimes and water levels in m NAVD88
- wave information file: should have datetime, significant wave height (m),
  and wave period (s) columns.

Input variables include:
- beach slope (to calculate runup)
- berm elevation in m NAVD88
- conversion from m NAVD88 to m MHW for your tidal gauge (check the NOAA
  datums for your gauge)
"""

# load necessary packages
import numpy as np
import pandas as pd
import os
from datetime import datetime, timedelta


# =============================================================================
# user inputs
# =============================================================================
start_time = '2004-01-01 00:00:00'  # date to start the storms
end_time = '2024-12-31 23:00:00'    # date to end the storms
water_levels_file = r"/scripts/input_prep/NOAA_water_level/8651370_DUCK_19840101_20241231_NAVD.csv"  # file that contains the water levels from NOAA gauge
wis_file = r"/data/hatteras_init/storms/WIS_raw_data/ST63228-Generic_Export-20260427T11T11_48.csv"  # file that contains the wave height and period from WIS gauge
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
max_storm_dur = 72      # maximum duration to include in storm events [hrs]
save_dfs = True         # determine whether to save the dataframes as csv and npy files
save_dir = r"/data/hatteras_init/storms/hindcast_storms/1984_2004"  # directory (no file name or extension) to save the storm files
save_name = "2004_2024_storms_v3_72"  # name of saved cascade storm files, if saved


# -----------------------------------------------------------------------------
# DO NOT EDIT CODE BLOCKS BELOW - all inputs are updated automatically. To
# change results, please change your inputs in the block above. However, you
# should review the comparison below to ensure it makes sense.
# -----------------------------------------------------------------------------


# =============================================================================
# Step 1: Load data and combine into a single dataframe
# =============================================================================
def find_time_gaps(df, col_name):
    """
    this function finds all values that are NaN and groups them by consecutive time periods
    NOTE: this function does not fill or modify missing data. it just informs you where the data is missing
    df: dataframe
        has datetime indeces and a column with data to evaluate
    col_name: string
        name of col to evaluate in df
    
    """
    
    # check for NaN value and create a datetime list of all corresponding NaNs
    nan_index_numbers = np.where(df[col_name].isnull())[0]  # returns numeric index
    time_list_nans = df.index[nan_index_numbers]  # returns the datetime index

    # loop through the NaN times and check if the difference between times is more than 1 hour
    # if so, that means this is the transition from one consecutive time period of nans to another
    # add the end value plus the next start value to our list of datetimes
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
    
    # separate out by start and end datetime
    # start dates will be odd indeces, end dates will be even indeces
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
    """
    This function loads the data, combines everything into a single dataframe with the specified start and stop values, 
    and drops rows with nans. It returns the combined dataframe.
    If data is missing from water level, wave height, or wave period, it will call the function above.
    
    start_time: string
        date to start the storms
    end_time: string
        date to end the storms
    water_levels_file: string
        file that contains the water levels from NOAA gauge
    wis_file: string
        file that contains the wave height and period from WIS gauge
    t_name_water: string, optional
        name of the column that contains the datetimes in the water levels file
    water_name: string, optional
        name of the column that contains the water levels in the water levels file [m NAVD88]
    t_name_wis: string, optional
        name of the column that contains the datetimes in the WIS file    
    waveHs_name: string, optional
        name of the column that contains the significant wave heights in the WIS file [m]
    waveTp_name: string, optional
        name of the column that contains the wave periods in the WIS file [s]
    """
    
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
    
    # check for missing data
    # water levels
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


# =============================================================================
# Step 2: Calculate the runup to use for the Total Water Level (TWL)
# =============================================================================
def calculate_r2_percent(Hs, Tp, slope):
    """Calculate R2% (2% exceedance runup) using Stockdon et al. (2006)"""
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


# =============================================================================
# Step 3: Create the CASCADE storms based on berm elevation and TWL
# =============================================================================
def create_storms(
    df_merged, 
    berm_elevation, 
    weather_grouping=24, 
    MHW=0.421,
    min_storm_dur=8,
    max_storm_dur=240,
    save_dfs=True,
    save_dir="",
    save_name=""
):
    
    """    
    df_merged: dataframe
        comparison dataframe from the previous function
    berm_elevation: float
        berm elevation of your location of interest [m NAVD88]
    weather_grouping: int, optional
        if storms occur within the specified limit, they are assumed part of same weather system and are grouped into one event [hrs]
    MHW: float, optional
        conversion from m NAVD88 to m MHW for your tidal gauge [m]
    min_storm_dur: int, optional
        minimum storm duration that is considered a storm event [hrs]
    max_storm_dur: int, optional
        maximum duration to include in storm events [hrs]
    save_dfs: bool, optional
        determine whether to save the dataframes as csv and npy files
    savedir: string, optional
        directory to save the storm files. if blank, use the current working directory
    save_name: string, optional
        name of this storm run
    
    """
    # -----------------------------------------------------------------------------
    # ---------------------------STEP 1 IDENTIFY STORMS ---------------------------
    # -----------------------------------------------------------------------------
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
    
    
    # -----------------------------------------------------------------------------
    # ---------------------------STEP 2 CREATE STORMS ---------------------------
    # -----------------------------------------------------------------------------
    storms = []
    storm_groups = df.dropna(subset=["StormID"]).groupby("StormID")

    for sid, group in storm_groups:

        group = group.sort_values("Time").copy()
        start_time = group["Time"].iloc[0]
        end_time   = group["Time"].iloc[-1]

        # Duration (hours)
        duration = len(group)  # each row is 1 hour where an exceedance occured

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
                "calendar_year": start_time.year,
                "StartTime": start_time,
                "EndTime": end_time,
                "Rhigh": rhigh,
                "Rlow": rlow,
                "period": period,
                "duration": duration
            })

    storms_df = pd.DataFrame(storms)

    # =========================
    # CONVERT YEAR → TIME (CASCADE)
    # =========================
    if not storms_df.empty:
        start_year = storms_df["calendar_year"].min()
        storms_df["time"] = storms_df["calendar_year"] - start_year + 1
    else:
        storms_df["time"] = pd.Series(dtype=float)

    # =========================
    # FINAL FORMAT
    # =========================
    cascade_df = storms_df[["time", "Rhigh", "Rlow", "period", "duration"]].copy()
    cascade_df = cascade_df.reset_index(drop=True)
    casc_input_array = cascade_df.to_numpy()
    total_storms = len(cascade_df)
    avg_duration = round(np.mean(cascade_df.duration.values),0)

    # =========================
    # SAVE
    # =========================
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
    save_name=save_name
)
