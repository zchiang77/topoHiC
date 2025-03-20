from numba import njit
import numpy as np
import h5py
import math
import pandas as pd
import time
import logging


@njit
def circulation(b1, b2, val, w, n, r1, r2, lim):
    
    nan_count = 0
    est_bin_1 = []
    est_bin_2 = []
    est_bin_3 = []
    est_bin_4 = []
    est_bin_5 = []
    est_bin_6 = []
    est_bin_7 = []
    est_bin_8 = []
    est_bin_9 = []
    est_bin_10 = []
    est_bin_20 = []
    est_bin_50 = []
    est_bin_100 = []
    est_bin_200 = []
    for i in range(n):
        if ((b1[i] < r1) or (b2[i] > r2)) or ((b2[i] - b1[i] > lim) and (lim != -1)):
            continue
        else:
            if b2[i] - b1[i] == 1:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_1.append(est_dis)
            elif b2[i] - b1[i] == 2:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_2.append(est_dis)
            elif b2[i] - b1[i] == 3:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_3.append(est_dis)
            elif b2[i] - b1[i] == 4:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_4.append(est_dis)
            elif b2[i] - b1[i] == 5:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_5.append(est_dis)
            elif b2[i] - b1[i] == 6:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_6.append(est_dis)
            elif b2[i] - b1[i] == 7:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_7.append(est_dis)
            elif b2[i] - b1[i] == 8:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_8.append(est_dis)
            elif b2[i] - b1[i] == 9:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_9.append(est_dis)
            elif b2[i] - b1[i] == 10:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_10.append(est_dis)
            elif b2[i] - b1[i] == 20:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_20.append(est_dis)
            elif b2[i] - b1[i] == 50:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_50.append(est_dis)
            elif b2[i] - b1[i] == 100:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_100.append(est_dis)
            elif b2[i] - b1[i] == 200:
                bal_val = val[i] * w[b1[i]] * w[b2[i]]
                est_dis = 1 / bal_val
                if math.isnan(est_dis):
                    nan_count += 1
                est_bin_200.append(est_dis)
            else:
                continue
    return est_bin_1, est_bin_2, est_bin_3, est_bin_4, est_bin_5, \
        est_bin_6, est_bin_7, est_bin_8, est_bin_9, est_bin_10, \
        est_bin_20, est_bin_50, est_bin_100, est_bin_200, \
        nan_count


def find_thresh(cool_file, exp_name, weight_col, range_start, range_end, save_path, edge_lim):

    logging.basicConfig(filename=exp_name + '.log', level=logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s:%(message)s')
    logger = logging.getLogger(__name__)

    logger.info("Step1: Finding thresholds.")

    range1 = range_start
    range2 = range_end

    start_time = time.time()
    h5 = h5py.File(cool_file, 'r')

    n_elements = h5['pixels']['bin1_id'].shape[0]
    logger.info("Total {} edges in the .cool file.".format(n_elements))

    weights = h5['bins'][weight_col][:]
    bin1_all = h5['pixels']['bin1_id'][:]
    bin2_all = h5['pixels']['bin2_id'][:]
    pixels_all = h5['pixels']['count'][:]

    (est_bin1, est_bin2, est_bin3, est_bin4, est_bin5,
     est_bin6, est_bin7, est_bin8, est_bin9, est_bin10,
     est_bin20, est_bin50, est_bin100, est_bin200,
     n_count) = circulation(bin1_all, bin2_all, pixels_all, weights, n_elements, range1, range2, edge_lim)

    dis_perc = pd.DataFrame([np.nanpercentile(est_bin1, [5, 10, 25, 50, 75, 90, 95]),
                             np.nanpercentile(est_bin2, [5, 10, 25, 50, 75, 90, 95]),
                             np.nanpercentile(est_bin3, [5, 10, 25, 50, 75, 90, 95]),
                             np.nanpercentile(est_bin4, [5, 10, 25, 50, 75, 90, 95]),
                             np.nanpercentile(est_bin5, [5, 10, 25, 50, 75, 90, 95]),
                             np.nanpercentile(est_bin6, [5, 10, 25, 50, 75, 90, 95]),
                             np.nanpercentile(est_bin7, [5, 10, 25, 50, 75, 90, 95]),
                             np.nanpercentile(est_bin8, [5, 10, 25, 50, 75, 90, 95]),
                             np.nanpercentile(est_bin9, [5, 10, 25, 50, 75, 90, 95]),
                             np.nanpercentile(est_bin10, [5, 10, 25, 50, 75, 90, 95]),
                             np.nanpercentile(est_bin20, [5, 10, 25, 50, 75, 90, 95]),
                             np.nanpercentile(est_bin50, [5, 10, 25, 50, 75, 90, 95]),
                             np.nanpercentile(est_bin100, [5, 10, 25, 50, 75, 90, 95]),
                             np.nanpercentile(est_bin200, [5, 10, 25, 50, 75, 90, 95])],
                            columns=['5%', '10%', '25%', '50%', '75%', '90%', '95%'],
                            index=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                                   20, 50, 100, 200])

    mean_list = [np.nanmean(est_bin1), np.nanmean(est_bin2),
                 np.nanmean(est_bin3), np.nanmean(est_bin4),
                 np.nanmean(est_bin5), np.nanmean(est_bin6),
                 np.nanmean(est_bin7), np.nanmean(est_bin8),
                 np.nanmean(est_bin9), np.nanmean(est_bin10),
                 np.nanmean(est_bin20), np.nanmean(est_bin50),
                 np.nanmean(est_bin100), np.nanmean(est_bin200)]

    dis_perc['mean'] = mean_list

    logger.info("{}({:.2%}) estimated distances are NaN.".format(n_count, n_count / n_elements))
    logger.info("Calculation of percentiles and means do not consider these NaN values.")

    results_save_path = save_path + exp_name + "-perc.csv"
    logger.info("Calculation results save to {}.".format(results_save_path))

    dis_perc.to_csv(results_save_path, float_format='%.2f', header=True, index=True)

    end_time = time.time()
    logger.info("Step1 completed! Time taken(s): {:.2}.".format(end_time - start_time))

    return results_save_path
