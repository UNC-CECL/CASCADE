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
# from tools.outwash_plotters import plot_Elev_CASCADE_subplots


# ---------------------------------- set model parameters that change per run ------------------------------------------
# rname_array = ["r025", "r035"]
rname_array = ["r025"]
for rname in rname_array:

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
        plot_ModelTransects(
            cascade=b3d_obj,
            time_step=np.arange(0,max_year, interval, dtype=int),
            iB3D=0)


        # # SHOREFACE PLOT
        # fig1, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(12, 12), sharex="all")  # create 4 subplots, one per scenario
        # for y in np.arange(0, max_year, interval):
        #     # x values are the shoreline and toe positions
        #     # y values are the SL (0) and shoreface depth?
        #     ax1.plot([m_xsTS_b3d[y], m_xtTS_b3d[y]], [0, -d_shore_b3d], label="year {0}".format(y), marker="o", color=blues[y])
        #     ax2.plot([m_xsTS_100[y], m_xtTS_100[y]], [0, -d_shore_100], label="year {0}".format(y), marker="o", color=blues[y])
        #     ax3.plot([m_xsTS_50[y], m_xtTS_50[y]], [0, -d_shore_50], label="year {0}".format(y), marker="o", color=blues[y])
        #     ax4.plot([m_xsTS_0[y], m_xtTS_0[y]], [0, -d_shore_0], label="year {0}".format(y), marker="o", color=blues[y])
        #     # plot features
        #     ax1.legend(ncols=legend_cols)
        #     ax1.set_title("overwash only")
        #     ax2.set_title("100% washout to shoreface")
        #     ax3.set_title("50% washout to shoreface")
        #     ax4.set_title("0% washout to shoreface")
        #
        #     ax4.set_xlabel("shoreface toe ---> shoreline (m)")
        #     ax1.set_ylabel("Elevation (m MHW)")
        #     ax2.set_ylabel("Elevation (m MHW)")
        #     ax3.set_ylabel("Elevation (m MHW)")
        #     ax4.set_ylabel("Elevation (m MHW)")
        #
        #
        #
        # fig1.tight_layout()
        # plt.show()
        #
        # # TRANSECT FOR ALL 4 SCENARIOS
        # fig2, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(12, 12), sharex="all")  # create 4 subplots, one per scenario
        # for y in np.arange(0, max_year, interval):
        #     # current ts domains converted to averages along the rows
        #     ts_domain_b3d = np.mean(domain_b3d[y], axis=1) * 10
        #     ts_domain_100 = np.mean(domain_100[y], axis=1) * 10
        #     ts_domain_50 = np.mean(domain_50[y], axis=1) * 10
        #     ts_domain_0 = np.mean(domain_0[y], axis=1) * 10
        #
        #     # only use elevations > 0
        #     pos_ts_domain_b3d = ts_domain_b3d[ts_domain_b3d>0]
        #     pos_ts_domain_100 = ts_domain_100[ts_domain_100>0]
        #     pos_ts_domain_50 = ts_domain_50[ts_domain_50>0]
        #     pos_ts_domain_0 = ts_domain_0[ts_domain_0>0]
        #
        #
        #     # x-values - shift with movement of the dunes
        #     # b3d
        #     x_start_b3d = abs(sum(sc_TS_b3d[0:y+1]))
        #     x_end_b3d = len(ts_domain_b3d) + x_start_b3d
        #     x_values_b3d = np.arange(x_start_b3d, x_end_b3d, 1)
        #     # 100%
        #     x_start_100 = abs(sum(sc_TS_100[0:y+1]))
        #     x_end_100 = len(ts_domain_100) + x_start_100
        #     x_values_100 = np.arange(x_start_100, x_end_100, 1)
        #     # 50%
        #     x_start_50 = abs(sum(sc_TS_50[0:y+1]))
        #     x_end_50 = len(ts_domain_50) + x_start_50
        #     x_values_50 = np.arange(x_start_50, x_end_50, 1)
        #     # 0%
        #     x_start_0 = abs(sum(sc_TS_0[0:y+1]))
        #     x_end_0 = len(ts_domain_0) + x_start_0
        #     x_values_0 = np.arange(x_start_0, x_end_0, 1)
        #
        #     # plot
        #     if y == 0:
        #         ls = "solid"
        #     else:
        #         ls = "solid"
        #     ax1.plot(x_values_b3d, ts_domain_b3d, label="year {0}".format(y), ls=ls, color=blues[y])
        #     ax2.plot(x_values_100, ts_domain_100, label="year {0}".format(y), ls=ls, color=blues[y])
        #     ax3.plot(x_values_50, ts_domain_50, label="year {0}".format(y), ls=ls, color=blues[y])
        #     ax4.plot(x_values_0, ts_domain_0, label="year {0}".format(y), ls=ls, color=blues[y])
        #     # plot features
        #     ax1.legend(ncols=legend_cols)
        #     ax1.set_title("overwash only")
        #     ax2.set_title("100% washout to shoreface")
        #     ax3.set_title("50% washout to shoreface")
        #     ax4.set_title("0% washout to shoreface")
        #
        #     ax4.set_xlabel("cross-shore distance of the barrier interior (ocean ---> bay) (m)")
        #     ax1.set_ylabel("Elevation (m MHW)")
        #     ax2.set_ylabel("Elevation (m MHW)")
        #     ax3.set_ylabel("Elevation (m MHW)")
        #     ax4.set_ylabel("Elevation (m MHW)")
        #
        # fig2.tight_layout()
        # plt.show()
        #
        # # TRANSECT FOR ALL 4 SCENARIOS - ONLY POS ELEVATIONS
        # fig3, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, figsize=(12, 12), sharex="all")  # create 4 subplots, one per scenario
        # for y in np.arange(0, max_year, interval):
        #     # current ts domains converted to averages along the rows
        #     ts_domain_b3d = np.mean(domain_b3d[y], axis=1) * 10
        #     ts_domain_100 = np.mean(domain_100[y], axis=1) * 10
        #     ts_domain_50 = np.mean(domain_50[y], axis=1) * 10
        #     ts_domain_0 = np.mean(domain_0[y], axis=1) * 10
        #
        #     # only use elevations > 0
        #     pos_ts_domain_b3d = ts_domain_b3d[ts_domain_b3d>0]
        #     pos_ts_domain_100 = ts_domain_100[ts_domain_100>0]
        #     pos_ts_domain_50 = ts_domain_50[ts_domain_50>0]
        #     pos_ts_domain_0 = ts_domain_0[ts_domain_0>0]
        #
        #
        #     # x-values - shift with movement of the dunes
        #     # b3d
        #     x_start_b3d = abs(sum(sc_TS_b3d[0:y+1]))
        #     x_end_b3d = len(pos_ts_domain_b3d) + x_start_b3d
        #     x_values_b3d = np.arange(x_start_b3d, x_end_b3d, 1)
        #     # 100%
        #     x_start_100 = abs(sum(sc_TS_100[0:y+1]))
        #     x_end_100 = len(pos_ts_domain_100) + x_start_100
        #     x_values_100 = np.arange(x_start_100, x_end_100, 1)
        #     # 50%
        #     x_start_50 = abs(sum(sc_TS_50[0:y+1]))
        #     x_end_50 = len(pos_ts_domain_50) + x_start_50
        #     x_values_50 = np.arange(x_start_50, x_end_50, 1)
        #     # 0%
        #     x_start_0 = abs(sum(sc_TS_0[0:y+1]))
        #     x_end_0 = len(pos_ts_domain_0) + x_start_0
        #     x_values_0 = np.arange(x_start_0, x_end_0, 1)
        #
        #     # plot
        #     if y == 0:
        #         ls = "solid"
        #     else:
        #         ls = "solid"
        #     ax1.plot(x_values_b3d, pos_ts_domain_b3d, label="year {0}".format(y), ls=ls, color=blues[y])
        #     ax2.plot(x_values_100, pos_ts_domain_100, label="year {0}".format(y), ls=ls, color=blues[y])
        #     ax3.plot(x_values_50, pos_ts_domain_50, label="year {0}".format(y), ls=ls, color=blues[y])
        #     ax4.plot(x_values_0, pos_ts_domain_0, label="year {0}".format(y), ls=ls, color=blues[y])
        #     # plot features
        #     ax1.legend(ncols=legend_cols)
        #     ax1.set_title("overwash only")
        #     ax2.set_title("100% washout to shoreface")
        #     ax3.set_title("50% washout to shoreface")
        #     ax4.set_title("0% washout to shoreface")
        #
        #     ax4.set_xlabel("cross-shore distance of the barrier interior (ocean ---> bay) (m)")
        #     ax1.set_ylabel("Elevation (m MHW)")
        #     ax2.set_ylabel("Elevation (m MHW)")
        #     ax3.set_ylabel("Elevation (m MHW)")
        #     ax4.set_ylabel("Elevation (m MHW)")
        #
        # fig3.tight_layout()
        # plt.show()
