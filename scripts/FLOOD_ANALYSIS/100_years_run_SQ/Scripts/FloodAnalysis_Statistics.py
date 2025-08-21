### Author: Roya Sahraei
### Flood Analysis for Ocracoke, Summer 2025

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import argparse
import matplotlib.patches as patches
from matplotlib.patches import Patch


#Main Analysis Function

def perform_detailed_analysis(npz_file_path, water_level_file_path):
    """
    Processes one .npz file and returns a DataFrame with flood days,
    flood depth, and abandonment status for each year and domain.
    """
    print(f"Processing: {os.path.basename(npz_file_path)}")

    #Load Data
    try:
        data = np.load(npz_file_path, allow_pickle=True)
        cascade_obj = data['cascade'][0]
        MHW_ELEVATION_NAVD88 = 0.26
        water_NAVD_df = pd.read_csv(water_level_file_path)
        daily_tidal_levels_NAVD = water_NAVD_df['v'].values
        daily_tidal_MHW = daily_tidal_levels_NAVD - MHW_ELEVATION_NAVD88
    except Exception as e:
        print(f"  ERROR: Failed to load data. Skipping. Details: {e}")
        return None

    if not hasattr(cascade_obj, 'roadways') or not cascade_obj.roadways:
        print(f"  WARNING: 'roadways' data not found. Skipping.")
        return None

    abandonment_years_dict = {
        i: roadway._time_index - 1
        for i, roadway in enumerate(cascade_obj.roadways) if roadway._time_index is not None
    }

    #Calculate All Metrics
    results_list = []
    for year in range(100):
        for i, road_segment in enumerate(cascade_obj.roadways):
            if 15 <= i <= 55 and year < len(road_segment._road_ele_TS):
                road_elevation = road_segment._road_ele_TS[year]

                is_flooded = road_elevation < daily_tidal_MHW
                flooded_days = np.sum(is_flooded)

                max_flood_depth = 0
                if flooded_days > 0:
                    flood_depths = daily_tidal_MHW[is_flooded] - road_elevation
                    max_flood_depth = np.max(flood_depths)

                abandon_year = abandonment_years_dict.get(i)
                is_abandoned = 1 if (abandon_year is not None and year >= abandon_year) else 0

                results_list.append({
                    'Year': year,
                    'Domain': i,
                    'Flooded_Days': flooded_days,
                    'Max_Flood_Depth_m_MHW': max_flood_depth,
                    'Is_Abandoned': is_abandoned,
                    'Abandonment_Year': abandon_year if abandon_year is not None else np.nan
                })

    return pd.DataFrame(results_list)


#Plotting Function

def plot_flooddays_heatmap(df, value_col, title, cbar_label, output_filename, output_dir, abandonment_data=None, cmap='YlOrRd', fmt = None,):
    """
    Generates and saves a heatmap for any given aggregated metric.
    Optionally overlays abandonment data as gray patches.
    """
    print(f"Generating flood days heatmap for: {title}...")
    heatmap_data = df.pivot_table(index='Domain', columns='Year', values=value_col)

    # Create a boolean mask where the data is 0. This will be used to make these cells white.
    mask = heatmap_data == 0


    plt.figure(figsize=(47, 30), dpi=150)
    ax = sns.heatmap(heatmap_data,
                     annot=True,
                     annot_kws={'size': 9},
                     fmt = ".0f",
                     cmap=cmap,
                     linewidths=.5,
                     cbar_kws={'label': cbar_label, 'pad': 0.01},
                     mask=mask)
    if abandonment_data:
        y_labels = list(heatmap_data.index)
        x_labels = list(heatmap_data.columns)

        for domain, abandon_year in abandonment_data.items():
            if domain in y_labels and abandon_year in x_labels:
                try:
                    row_idx = y_labels.index(domain)
                    start_col_idx = x_labels.index(abandon_year)
                    width = len(x_labels) - start_col_idx
                    ax.add_patch(patches.Rectangle(
                        (start_col_idx, row_idx), width, 1,
                        facecolor='gray', alpha=0.6, lw=0
                    ))
                except (ValueError, IndexError):
                    continue

        legend_patch = Patch(facecolor='gray', alpha=0.6, label='Mean Year of Abandonment')
        plt.legend(handles=[legend_patch], loc='upper left', fontsize=20)

    ax.collections[0].colorbar.set_label(cbar_label, fontsize=20)
    ax.tick_params(axis='y', labelsize=15)
    ax.tick_params(axis='x', labelsize=15)

    plt.title(title, fontsize=25)
    plt.ylabel('Domain Number', fontsize=20)
    plt.xlabel('Simulation Year', fontsize=20)

    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  -> Saved to {output_filename}")

def plot_flood_depth_heatmap(df, value_col, title, cbar_label, output_filename, output_dir, abandonment_data=None, cmap='YlOrRd', fmt = None,):
    """
    Generates and saves a heatmap for any given aggregated metric.
    Optionally overlays abandonment data as gray patches.
    """
    print(f"Generating flood depth heatmap for: {title}...")
    heatmap_data = df.pivot_table(index='Domain', columns='Year', values=value_col)

    mask = heatmap_data == 0



    plt.figure(figsize=(47, 30), dpi=150)
    ax = sns.heatmap(heatmap_data,
                     annot=True,
                     annot_kws={'size': 9},
                     fmt = ".1f",
                     cmap=cmap,
                     linewidths=.5,
                     cbar_kws={'label': cbar_label, 'pad': 0.01},
                     mask=mask)

    if abandonment_data:
        y_labels = list(heatmap_data.index)
        x_labels = list(heatmap_data.columns)

        for domain, abandon_year in abandonment_data.items():
            if domain in y_labels and abandon_year in x_labels:
                try:
                    row_idx = y_labels.index(domain)
                    start_col_idx = x_labels.index(abandon_year)
                    width = len(x_labels) - start_col_idx
                    ax.add_patch(patches.Rectangle(
                        (start_col_idx, row_idx), width, 1,
                        facecolor='gray', alpha=0.6, lw=0
                    ))
                except (ValueError, IndexError):
                    continue

        legend_patch = Patch(facecolor='gray', alpha=0.6, label='Mean Year of Abandonment')
        plt.legend(handles=[legend_patch], loc='upper left', fontsize=20)

    ax.collections[0].colorbar.set_label(cbar_label, fontsize=20)
    ax.tick_params(axis='y', labelsize=15)
    ax.tick_params(axis='x', labelsize=15)

    plt.title(title, fontsize=25)
    plt.ylabel('Domain Number', fontsize=20)
    plt.xlabel('Simulation Year', fontsize=20)

    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  -> Saved to {output_filename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate flood analysis for CASCADE scenarios.")
    parser.add_argument("project_dir", help="The project directory containing 'Input_Data' and 'Daily_WL_NAVD.csv'.")
    args = parser.parse_args()

    # Setup paths
    input_data_dir = os.path.join(args.project_dir, "Input_Data")
    results_dir = os.path.join(args.project_dir, "Results")
    wl_path = os.path.join(args.project_dir, "Daily_WL_NAVD.csv")
    os.makedirs(results_dir, exist_ok=True)

    # Validate paths
    if not os.path.isdir(input_data_dir) or not os.path.exists(wl_path):
        print(f"FATAL ERROR: Required files/folders not found in '{args.project_dir}'")
        exit()

    print(f"Starting analysis in: {input_data_dir}")
    print("-" * 60)

    # Collect results from all .npz files
    all_results_list = []
    npz_files = glob.glob(os.path.join(input_data_dir, '**', '*.npz'), recursive=True)

    if not npz_files:
        print("FATAL ERROR: No .npz files found.")
        exit()

    for file_path in npz_files:
        single_run_results = perform_detailed_analysis(file_path, wl_path)
        if single_run_results is not None and not single_run_results.empty:
            all_results_list.append(single_run_results)

    # Aggregate all metrics (mean and median)
    if not all_results_list:
        print("Processing complete, but no valid data was generated.")
    else:
        print("\n" + "-" * 60)
        print(f"Aggregating results from {len(all_results_list)} files...")

        combined_df = pd.concat(all_results_list, ignore_index=True)

        aggregation_rules = {
            'Flooded_Days': ['mean', 'median'],
            'Max_Flood_Depth_m_MHW': ['mean', 'median'],
            'Is_Abandoned': 'mean'
        }

        aggregated_data = combined_df.groupby(['Year', 'Domain']).agg(aggregation_rules).reset_index()

        aggregated_data.columns = ['_'.join(col).strip() if isinstance(col, tuple) and col[1] else col[0] for col in
                                   aggregated_data.columns.values]
        aggregated_data = aggregated_data.rename(columns={'Is_Abandoned_mean': 'Abandonment_Probability'})

        output_csv_path = os.path.join(results_dir, 'Aggregated_Results_All_Metrics.csv')
        aggregated_data.to_csv(output_csv_path, index=False)
        print(f"Full aggregated results saved to: {os.path.basename(output_csv_path)}")

        # Calculate the Mean abandonment year for each domain
        mean_abandon_years = combined_df.dropna(subset=['Abandonment_Year']).groupby('Domain')[
            'Abandonment_Year'].mean().round().astype(int).to_dict()
        print("Calculated mean abandonment year for each domain.")

        # mean (averaged) heatmaps for each metric ---
        print("\n" + "-" * 60)
        plot_flooddays_heatmap(df=aggregated_data,
                                value_col='Flooded_Days_mean',
                                fmt=".0f",
                                title='Average Flood Days Throughout 100 Years',
                                cbar_label='Average Annual Flood Days',
                                output_filename='Average_Flood_Days.png',
                                output_dir=results_dir)
        #Median Flood Days
        plot_flooddays_heatmap(df=aggregated_data,
                                value_col='Flooded_Days_median',
                                fmt=".0f",
                                title='Median Flood Days Throughout 100 Years',
                                cbar_label='Median Annual Flood Days',
                                output_filename='Median_Flood_Days.png',
                                output_dir=results_dir)


        plot_flooddays_heatmap(df=aggregated_data,
                                value_col='Abandonment_Probability',
                                title='Average year of Road Abandonment Throughout 100 Years',
                                cbar_label='Average Abandonment Year',
                                output_filename='Average_Abandonment_Year.png',
                                output_dir=results_dir,
                                abandonment_data=mean_abandon_years)

        plot_flood_depth_heatmap(df=aggregated_data,
                                value_col='Max_Flood_Depth_m_MHW_mean',
                                title='Average Max Flood Depth Throughout 100 Years',
                                cbar_label='Average Max Flood Depth (m MHW)',
                                cmap='YlGnBu',
                                output_filename='Average_Max_Flood_Depth.png',
                                output_dir=results_dir)

    print("\n" + "-" * 60)
    print("All processing complete.")