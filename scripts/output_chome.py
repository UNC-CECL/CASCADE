import numpy as np
import os
from scipy.io import loadmat
import matplotlib.pyplot as plt


def get_statistics_4_chome(iB3D, tmax_mgmt, name_prefix):

    os.chdir(
        "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output"
    )

    output = np.load(name_prefix + ".npz", allow_pickle=True)
    cascade = output["cascade"]
    cascade = cascade[0]
    barrier3d = cascade.barrier3d

    # Barrier3D in decameters --> convert to meters for CHOME
    barrier_height = []
    for t in range(0, tmax_mgmt):
        bh_array = np.array(barrier3d[iB3D].DomainTS[t]) * 10  # m MHW
        barrier_height.append(bh_array[bh_array > 0].mean())

    # shoreline change rate (m/yr) [negative numbers correspond to nourishments]
    scts = [
        (x - barrier3d[iB3D].x_s_TS[0]) * 10
        for x in barrier3d[iB3D].x_s_TS[0:tmax_mgmt]
    ]
    rate = [0]
    for k in range(1, len(scts)):
        rate.append(scts[k] - scts[k - 1])
    sc_rate = rate

    # maximum height of each row in DuneDomain, then average (m MHW)
    DuneDomainCrest = barrier3d[iB3D].DuneDomain[0:tmax_mgmt, :, :].max(axis=2)
    DuneRestart = barrier3d[iB3D].DuneRestart
    DuneDomainCrest[DuneDomainCrest < DuneRestart] = DuneRestart
    dune_crest_mean = (
        np.mean(DuneDomainCrest, axis=1) + barrier3d[iB3D].BermEl
    ) * 10  # m MHW

    # barrier width (meters)
    barrier_width = (
        np.array(barrier3d[iB3D].x_b_TS[0:tmax_mgmt])
        - np.array(barrier3d[iB3D].x_s_TS[0:tmax_mgmt])
    ) * 10  # includes beach width
    barrier_width[1:tmax_mgmt] = np.array(barrier_width[1:tmax_mgmt]) - np.array(
        cascade.nourishments[iB3D].beach_width[1:tmax_mgmt]
    )

    return barrier_width, dune_crest_mean, sc_rate, barrier_height


# name = "scenario_Split_Mgmt.mat"
name = "scenario_Rodanthe.mat"
name = "scenario_Rodanthe_acceratedSLR"
name = "scenario_subsidy_50p"


def plot_chome_output(name_prefix):

    os.chdir(
        "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output"
    )
    X_OF = loadmat(name)["X_OF"]
    X_NOF = loadmat(name)["X_NOF"]
    MMT = loadmat(name)["MMT"]
    ACOM = loadmat(name)["ACOM"]
    oceanfront_price = np.squeeze(X_OF["price"][0][0])
    nonoceanfront_price = np.squeeze(X_NOF["price"][0][0])
    oceanfront_marketshare = np.squeeze(X_OF["mkt"][0][0])
    nonoceanfront_marketshare = np.squeeze(X_NOF["mkt"][0][0])
    beach_width = np.squeeze(MMT["bw"][0][0])
    expected_beach_width = np.squeeze(ACOM["Ebw"][0][0])

    # each run has 100 years of warm up
    plt.plot(oceanfront_price[100:-1])
    plt.plot(nonoceanfront_price[100:-1])
    plt.legend(["oceanfront", "non-oceanfront"])
    plt.xlabel("time (yr))")
    plt.ylabel("property value (USD))")

    # if it’s 60% then 60% of residents are renters (and the other 40% is homeowners who live in the home)
    plt.plot(oceanfront_marketshare[100:-1])
    plt.plot(nonoceanfront_marketshare[100:-1])
    plt.legend(["oceanfront", "non-oceanfront"])
    # plt.legend(["oceanfront", "non-oceanfront", "OF - acc SLR", "NOF - acc SLR", "OF - 50% subsidy", "NOF - 50% subsidy"])
    plt.xlabel("time (yr))")
    plt.ylabel("investor market share (% renters)")

    # beach width through time is
    plt.plot(beach_width[100:-1])
    plt.plot(expected_beach_width[100:-1])
    plt.legend(["actual", "expected"])
    plt.xlabel("time (yr))")
    plt.ylabel("beach width (m)")

    # t1 = 2
    # t2 = 100 + length(er_rate)
    # for ti=t1:t2
    #     NOFtauO(ti) = mean(SV_NOF.tau_incO{ti});
    #     NOFtauR(ti) = mean(SV_NOF.tau_incR{ti});
    #     OFtauO(ti)  = mean(SV_OF.tau_incO{ti});
    #     OFtauR(ti)  = mean(SV_OF.tau_incR{ti});
    # end
    # plot(OFtauO) % oceanfront resident home owner average income
    # hold on
    # plot(OFtauR) % oceanfront renter average income
    # plot(NOFtauO) % nonoceanfront resident home owner average income
    # plot(NOFtauR) % nonoceanfront renter average income


# split management - 1 m background erosion and 0.004 m/yr SLR ----------------------------------------
name_prefix = "9-CASCADE_Rave_pt75_Nourish_2mDune_lowEle_comm_BE1m_RT1m_6AST_3roads"
iB3D = 1  # middle community
tmax_mgmt = 80  # stopped community management at 80 years

[barrier_width, dune_crest_mean, sc_rate, barrier_height] = get_statistics_4_chome(
    iB3D, tmax_mgmt, name_prefix
)

# save to text file for use in Matlab
np.savetxt(
    "Split_Mgmt.txt",
    (
        np.array(barrier_width),
        np.array(dune_crest_mean),
        np.array(sc_rate),
        np.array(barrier_height),
    ),
)

# Rodanthe scenario - 1 m background erosion and 0.004 m/yr SLR ----------------------------------------
name_prefix = name = "9-CASCADE_AST_3domains_BE1m"
iB3D = 1  # middle community
tmax_mgmt = 199  # never stopped community management

[barrier_width, dune_crest_mean, sc_rate, barrier_height] = get_statistics_4_chome(
    iB3D, tmax_mgmt, name_prefix
)

# save to text file for use in Matlab
np.savetxt(
    "Rodanthe.txt",
    (
        np.array(barrier_width),
        np.array(dune_crest_mean),
        np.array(sc_rate),
        np.array(barrier_height),
    ),
)

# Rodanthe scenario -  # 1 m background erosion and accelerated SLR ----------------------------------------
name_prefix = name = "9-CASCADE_AST_3domains_BE1m_AccSLR"
iB3D = 1  # middle community
tmax_mgmt = 88  # middle roadway drowned

[barrier_width, dune_crest_mean, sc_rate, barrier_height] = get_statistics_4_chome(
    iB3D, tmax_mgmt, name_prefix
)

# save to text file for use in Matlab
np.savetxt(
    "Rodanthe_acceratedSLR.txt",
    (
        np.array(barrier_width),
        np.array(dune_crest_mean),
        np.array(sc_rate),
        np.array(barrier_height),
    ),
)
