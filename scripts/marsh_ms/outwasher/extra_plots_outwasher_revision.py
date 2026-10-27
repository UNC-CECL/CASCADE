# Lexi Van Blunk
# 11/22/2025
# plots for revisions to journal article:
# Barrier vulnerability following outwash: A balance of overwash and dune gap recovery
# which include shoreline position transects over time for one overwash storm and both dune growth rates
# and maybe I will also have plan-view plots

# --------------------------------- description of variables ----------------------------------------------------------
# there are eight main scenarios: overwash only, 100, 50, or 0% washout to shoreface and two dune growth rates each
# each scenario was run with 100 different storm sequences for a model duration of 100 years or until drowning
# a single storm sequence consists of an overwash storm series and potentially outwash storm series
# summary of each scenario:
# overwash only, r=0.25, 100 storm sequences
# overwash only, r=0.35, 100 storm sequences
# 100% WTS, r=0.25, 100 storm sequences
# 100% WTS, r=0.35, 100 storm sequences
# 50% WTS, r=0.25, 100 storm sequences
# 50% WTS, r=0.35, 100 storm sequences
# 0% WTS, r=0.25, 100 storm sequences
# 0% WTS, r=0.35, 100 storm sequences


import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colormaps

from cascade.tools.plotters import plot_ElevAnimation_CASCADE, plot_dune_domain, plot_ModelTransects
from tools.outwash_plotters import (
    plot_Elev_CASCADE_subplots,
    plot_Elev_CASCADE_subplots_all_scenarios,
    plot_Elev_CASCADE_subplots_all_scenarios_overwash,
)

plot_transects = False
plot_subplots = False
plot_subplots_all = False
plot_subplots_all_Qow = True

# ---------------------------------- set model parameters that change per run ------------------------------------------
# rname_array = ["r025", "r035"]
rname_array = ["r025"]
for rname in rname_array:

    if rname=="r025":
        rval = 0.25
    else:
        rval = 0.35

    # naming parameters
    storm_interval = 20
    config = 4  

    # location of the npz files
    datadir_b3d = r"C:\Users\agfig\Downloads\paper_results\paper_results/{}/".format(
        rname
    )
    datadir_100 = r"C:\Users\agfig\Downloads\paper_results\paper_results/{}/".format(
        rname
    )
    datadir_50 = r"C:\Users\agfig\Downloads\paper_results\paper_results/{}/".format(
        rname
    )
    datadir_0 = r"C:\Users\agfig\Downloads\paper_results\paper_results/{}/".format(
        rname
    )

    for storm_num in range(1, 2):  # only plotting the first storm

        # b3d variables ------------------------------------------------------------------------------------------------
        filename_b3d = (
            "config{}_b3d_startyr1_interval{}yrs_Slope0pt03_{}.npz".format(
                config, storm_interval, storm_num
            )
        )
        file_b3d = datadir_b3d + filename_b3d
        b3d = np.load(file_b3d, allow_pickle=True)
        b3d_obj = b3d["cascade"][0]
        tmax_b3d = b3d_obj.barrier3d[0].TMAX + 1

        # drowning variable - did this barrier drown?
        drowning_array_b3d = b3d_obj.barrier3d[0].drown_break

        # shoreline change (movement of the dunes/interior, NOT based on the beach)
        sc_TS_b3d = b3d_obj.barrier3d[0].ShorelineChangeTS 

        # shoreline position - if coming from barrier3d, they are in dam, if brie, they are in meters
        # m_xsTS_b3d = np.subtract(
        #     b3d_obj.barrier3d[0].x_s_TS, b3d_obj.barrier3d[0].x_s_TS[0]
        # )  
        m_xsTS_b3d = b3d_obj.barrier3d[0].x_s_TS
        m_xsTS_b3d = np.multiply(m_xsTS_b3d, 10)

        # toe position
        # m_xtTS_b3d = np.subtract(
        #     b3d_obj.barrier3d[0].x_t_TS, b3d_obj.barrier3d[0].x_t_TS[0]
        # ) 
        m_xtTS_b3d = b3d_obj.barrier3d[0].x_t_TS
        m_xtTS_b3d = np.multiply(m_xtTS_b3d, 10)

        # shoreface depth - dam
        d_shore_b3d =  b3d_obj.barrier3d[0].DShoreface * 10

        # overwash
        QowTS_b3d = b3d_obj.barrier3d[0].QowTS  # m^3/m, all model years for this storm sequence

        # barrier interior
        domain_b3d = b3d_obj.barrier3d[0].DomainTS  # dam


        # 100% variables  ---------------------------------------------------------------------------------------------
        filename_100 = (
            "config{}_outwash100_startyr1_interval{}yrs_Slope0pt03_{}.npz".format(
                config, storm_interval, storm_num
            )
        )
        file_100 = datadir_100 + filename_100
        outwash100 = np.load(file_100, allow_pickle=True)
        outwash100_obj = outwash100["cascade"][0]
        tmax_100 = outwash100_obj.barrier3d[0].TMAX + 1

        # drowning variable
        drowning_array_100 = outwash100_obj.barrier3d[0].drown_break

        # shoreline change (movement of the dunes/interior, NOT based on the beach)
        sc_TS_100 = outwash100_obj.barrier3d[0].ShorelineChangeTS 
        
        # shoreline position - only for barriers that do not drown
        # m_xsTS_100 = np.subtract(
        #     outwash100_obj.barrier3d[0].x_s_TS, outwash100_obj.barrier3d[0].x_s_TS[0]
        # )
        m_xsTS_100 = outwash100_obj.barrier3d[0].x_s_TS
        m_xsTS_100 = np.multiply(m_xsTS_100, 10)
               
        # toe position
        # m_xtTS_100 = np.subtract(
        #     outwash100_obj.barrier3d[0].x_t_TS, outwash100_obj.barrier3d[0].x_t_TS[0]
        # )
        m_xtTS_100 =outwash100_obj.barrier3d[0].x_t_TS
        m_xtTS_100 = np.multiply(m_xtTS_100, 10)

        # shoreface depth - dam
        d_shore_100 =  outwash100_obj.barrier3d[0].DShoreface * 10

        # overwash
        QowTS_100 = outwash100_obj.barrier3d[0].QowTS  # m^3/m, all model years for this storm sequence

        # outwash
        QoutTS_100 = outwash100_obj.outwash[0]._outwash_flux_TS  # m^3/m, all model years for this storm sequence
        QoutTS_100 = QoutTS_100[QoutTS_100 != 0]  # remove all nonzero values (they should have been initialized as nans)

        # barrier interior
        domain_100 = outwash100_obj.barrier3d[0].DomainTS


        # 50% variables ------------------------------------------------------------------------------------------------
        filename_50 = (
            "config{}_outwash50_startyr1_interval{}yrs_Slope0pt03_{}.npz".format(
                config, storm_interval, storm_num
            )
        )
        file_50 = datadir_50 + filename_50
        outwash50 = np.load(file_50, allow_pickle=True)
        outwash50_obj = outwash50["cascade"][0]
        tmax_50 = outwash50_obj.barrier3d[0].TMAX + 1

        # drowning variable
        drowning_array_50 = outwash50_obj.barrier3d[0].drown_break

        # shoreline change (movement of the dunes/interior, NOT based on the beach)
        sc_TS_50 = outwash50_obj.barrier3d[0].ShorelineChangeTS 
        
        # shoreline position - only for barriers that do not drown
        # m_xsTS_50 = np.subtract(
        #     outwash50_obj.barrier3d[0].x_s_TS, outwash50_obj.barrier3d[0].x_s_TS[0]
        # )
        m_xsTS_50 = outwash50_obj.barrier3d[0].x_s_TS
        m_xsTS_50 = np.multiply(m_xsTS_50, 10)

        # toe position
        # m_xtTS_50 = np.subtract(
        #     outwash50_obj.barrier3d[0].x_t_TS, outwash50_obj.barrier3d[0].x_t_TS[0]
        # )
        m_xtTS_50 = outwash50_obj.barrier3d[0].x_t_TS
        m_xtTS_50 = np.multiply(m_xtTS_50, 10)

        # shoreface depth - dam
        d_shore_50 =  outwash50_obj.barrier3d[0].DShoreface * 10

        # overwash
        QowTS_50 = outwash50_obj.barrier3d[0].QowTS  # m^3/m, all model years for this storm sequence

        # outwash
        QoutTS_50 = outwash50_obj.outwash[0]._outwash_flux_TS  # m^3/m, all model years for this storm sequence
        QoutTS_50 = QoutTS_50[QoutTS_50 != 0]

        # barrier interior
        domain_50 = outwash50_obj.barrier3d[0].DomainTS


        # 0% variables -------------------------------------------------------------------------------------------------
        filename_0 = (
            "config{}_outwash0_startyr1_interval{}yrs_Slope0pt03_{}.npz".format(
                config, storm_interval, storm_num
            )
        )
        file_0 = datadir_0 + filename_0
        outwash0 = np.load(file_0, allow_pickle=True)
        outwash0_obj = outwash0["cascade"][0]
        tmax_0 = outwash0_obj.barrier3d[0].TMAX + 1

        # drowning variable
        drowning_array_0 = outwash0_obj.barrier3d[0].drown_break

        # shoreline change (movement of the dunes/interior, NOT based on the beach)
        sc_TS_0 = outwash0_obj.barrier3d[0].ShorelineChangeTS 

        # shoreline position - only for barriers that do not drown
        # m_xsTS_0 = np.subtract(
        #     outwash0_obj.barrier3d[0].x_s_TS, outwash0_obj.barrier3d[0].x_s_TS[0]
        # )
        m_xsTS_0 = outwash0_obj.barrier3d[0].x_s_TS
        m_xsTS_0 = np.multiply(m_xsTS_0, 10)

        # toe position
        # m_xtTS_0 = np.subtract(
        #     outwash0_obj.barrier3d[0].x_t_TS, outwash0_obj.barrier3d[0].x_t_TS[0]
        # )
        m_xtTS_0 = outwash0_obj.barrier3d[0].x_t_TS
        m_xtTS_0 = np.multiply(m_xtTS_0, 10)

        # shoreface depth - dam
        d_shore_0 =  outwash0_obj.barrier3d[0].DShoreface * 10


        # overwash
        QowTS_0 = outwash0_obj.barrier3d[0].QowTS  # m^3/m, all model years for this storm sequence

        # outwash
        QoutTS_0 = outwash0_obj.outwash[0]._outwash_flux_TS  # m^3/m, all model years for this storm sequence
        QoutTS_0 = QoutTS_0[QoutTS_0 != 0]

        # barrier interior
        domain_0 = outwash0_obj.barrier3d[0].DomainTS


        # ------------------------------------------------------------------------------------------------------------------
        # now that we have a value for each scenario, make the plots
        plt.rcParams["font.size"] = 14
        length = 10
        height = 5  
        max_year = 21
        legend_cols = 3
        interval = 2  # years

        # plot all in single color for ease of passing time
        n = int(np.ceil(max_year / interval))
        blues_cmap = colormaps['Blues']
        blues = blues_cmap(np.linspace(0.2, 1, n))  # n evenly spaced colors

        # the cascade plotter for transects. plots the 10th transect
        if plot_transects:
            plot_ModelTransects(
                cascade=b3d_obj,
                time_step=np.arange(0,max_year, interval, dtype=int),
                iB3D=0)

        # subplots
        if plot_subplots:
            version = "altered_cmap"
            cols = 6
            years = [0, 1, 5, 10, 15, 20, 21, 41, 61, 81, 100]
            skip=False
            plot_overwash_text=False
            overwash_text_years=[5, 10, 15, 20]
            # barrier3d
            plot_Elev_CASCADE_subplots(
                cascade=b3d_obj,
                directory=r"C:\Users\agfig\Downloads\paper_results\paper_results\new_planview_subplots",
                TMAX_MGMT=0,
                name="{0}_storm{1}_b3d_{2}".format(rname, storm_num, version),
                TMAX_SIM=tmax_b3d,
                plot_timesteps=years,  # array with the timesteps you want to plot
                main_plot_title="baseline",  # title of the main plot (string)
                ny=1,
                beach_management_ny=None,  # list of bool the length of ny, or None for all False
                roadway_management_ny=None,
                y_lim=(150, 210),
                z_lim=6,
                z_bot=-3,
                fig_size=(15, 6),
                fig_eps=False,
                km_on=False,
                meters_on=True,
                outwash_scenario=False,
                n_cols=cols,
                skip_7=skip,
                overwash_text=plot_overwash_text,
                overwash_text_years=overwash_text_years
            )
            # 100%
            plot_Elev_CASCADE_subplots(
                cascade=outwash100_obj,
                directory=r"C:\Users\agfig\Downloads\paper_results\paper_results\new_planview_subplots",
                TMAX_MGMT=tmax_100,
                name="{0}_storm{1}_100_{2}".format(rname, storm_num, version),
                TMAX_SIM=tmax_100,
                plot_timesteps=years, # array with the timesteps you want to plot
                main_plot_title="100% washout to shoreface",  # title of the main plot (string)
                ny=1,
                beach_management_ny=None,  # list of bool the length of ny, or None for all False
                roadway_management_ny=None,
                y_lim=(150, 210),
                z_lim=6,
                z_bot=-3,
                fig_size=(15, 6),
                fig_eps=False,
                km_on=False,
                meters_on=True,
                outwash_scenario=True,
                n_cols=cols,
                skip_7=skip,
                overwash_text=plot_overwash_text,
                overwash_text_years=overwash_text_years
            )
            # 50%
            plot_Elev_CASCADE_subplots(
                cascade=outwash50_obj,
                directory=r"C:\Users\agfig\Downloads\paper_results\paper_results\new_planview_subplots",
                TMAX_MGMT=tmax_50,
                name="{0}_storm{1}_50_{2}".format(rname, storm_num, version),
                TMAX_SIM=tmax_50,
                plot_timesteps=years,  # array with the timesteps you want to plot
                main_plot_title="50% washout to shoreface",  # title of the main plot (string)
                ny=1,
                beach_management_ny=None,  # list of bool the length of ny, or None for all False
                roadway_management_ny=None,
                y_lim=(150, 210),
                z_lim=6,
                z_bot=-3,
                fig_size=(15, 6),
                fig_eps=False,
                km_on=False,
                meters_on=True,
                outwash_scenario=True,
                n_cols=cols,
                skip_7=skip,
                overwash_text=plot_overwash_text,
                overwash_text_years=overwash_text_years
            )
            # 0%
            plot_Elev_CASCADE_subplots(
                cascade=outwash0_obj,
                directory=r"C:\Users\agfig\Downloads\paper_results\paper_results\new_planview_subplots",
                TMAX_MGMT=tmax_0,
                name="{0}_storm{1}_0_{2}".format(rname, storm_num, version),
                TMAX_SIM=tmax_0,
                plot_timesteps=years,  # array with the timesteps you want to plot
                main_plot_title="0% washout to shoreface",  # title of the main plot (string)
                ny=1,
                beach_management_ny=None,  # list of bool the length of ny, or None for all False
                roadway_management_ny=None,
                y_lim=(150, 210),
                z_lim=6,
                z_bot=-3,
                fig_size=(15, 6),
                fig_eps=False,
                km_on=False,
                meters_on=True,
                outwash_scenario=True,
                n_cols=cols,
                skip_7=skip,
                overwash_text=plot_overwash_text,
                overwash_text_years=overwash_text_years
            )

        if plot_subplots_all:
            version = "all_scenarios_1row"
            cols = 11
            years = [0, 1, 5, 10, 15, 20, 21, 41, 61, 81, 100]
            plot_overwash_text=False
            overwash_text_years=[5, 10, 15, 20]
            plot_Elev_CASCADE_subplots_all_scenarios(
                cascade_list=[b3d_obj, outwash100_obj, outwash50_obj, outwash0_obj],
                directory=r"C:\Users\agfig\Downloads\paper_results\paper_results\new_planview_subplots",
                TMAX_MGMT_list=[tmax_b3d,tmax_100,tmax_50,tmax_0],
                name="{0}_storm{1}_{2}".format(rname, storm_num, version),
                TMAX_SIM_list=[tmax_b3d,tmax_100,tmax_50,tmax_0],
                plot_timesteps=years,  # array with the timesteps you want to plot
                main_plot_title=["baseline","100% washout to shoreface","50% washout to shoreface", "0% washout to shoreface"],  # title of the main plot (string)
                ny=1,
                beach_management_ny=None,  # list of bool the length of ny, or None for all False
                roadway_management_ny=None,
                y_lim=(150, 210),
                z_lim=6,
                z_bot=-3,
                fig_size=(22, 12),  # (13, 20) for 2 rows, 6 col
                fig_eps=False,
                km_on=False,
                meters_on=True,
                outwash_scenario=True,
                n_cols=cols,
                overwash_text=plot_overwash_text,
                overwash_text_years=overwash_text_years,
            )

        if plot_subplots_all_Qow:
            version = "all_scenarios_1row_qow"
            cols = 11
            years = [0, 1, 5, 10, 15, 20, 21, 41, 61, 81, 100]
            plot_overwash_text=False
            overwash_text_years=[5, 10, 15, 20]
            plot_Elev_CASCADE_subplots_all_scenarios_overwash(
                cascade_list=[b3d_obj, outwash100_obj, outwash50_obj, outwash0_obj],
                directory=r"C:\Users\agfig\Downloads\paper_results\paper_results\new_planview_subplots",
                TMAX_MGMT_list=[tmax_b3d,tmax_100,tmax_50,tmax_0],
                name="{0}_storm{1}_{2}".format(rname, storm_num, version),
                TMAX_SIM_list=[tmax_b3d,tmax_100,tmax_50,tmax_0],
                plot_timesteps=years,  # array with the timesteps you want to plot
                main_plot_title=["baseline","100% washout to shoreface","50% washout to shoreface", "0% washout to shoreface"],  # title of the main plot (string)
                ny=1,
                beach_management_ny=None,  # list of bool the length of ny, or None for all False
                roadway_management_ny=None,
                y_lim=(150, 210),
                z_lim=6,
                z_bot=-3,
                fig_size=(22, 15),  # (13, 20) for 2 rows, 6 col
                fig_eps=False,
                km_on=False,
                meters_on=True,
                outwash_scenario=True,
                n_cols=cols,
                overwash_text=plot_overwash_text,
                overwash_text_years=overwash_text_years,
                r_value=rval
            )

