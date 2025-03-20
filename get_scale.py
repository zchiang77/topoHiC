import pandas as pd


def get_scale(find_thresh_output, scale_perc, scale_dist):

    perc_data = pd.read_csv(find_thresh_output, header=0, index_col=0)
    birth_least = perc_data[scale_perc][scale_dist]

    return birth_least

