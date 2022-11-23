import csv
import numpy as np

def generate_outwash_hydrographs(
        directory="./data/outwash_data",
        input_file="sound_data.txt",
        output_file="outwash_storm_series_Dorian"
):
    # directory = "D:/NC State/Outwasher/chris stuff/sound_data.txt"
    fid = directory + "/" + input_file
    output_fid = directory + "/" + output_file

    with open(directory, newline='') as csvfile:
        sound_data = list(csv.reader(csvfile))[0]

    # [dam MHW] Chris' sound elevations were in m MSL, so converted to NAVD88 then MHW and dam
    sound_data = [float(s) / 10 - 0.054 for s in sound_data]
    sound_data = [s + 0.05 for s in sound_data]  # [dam MHW] just increasing the values

    # setting all negative values to 0
    sound_data = sound_data[20:44]
    sound_data[0] = 0

    # for now, make a synthetic storm series: year, hydrograph, duration
    storm_series = [1, sound_data, len(sound_data)]

    np.save(output_fid, storm_series)
