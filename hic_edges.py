import numpy as np
import math
import pandas as pd
import h5py
from numba import njit
import time
import logging


@njit
def search_edges(n_elements, bin1, bin2, pixels, weights, thresh, edges, count, r1, r2, lim):
    nanc = 0

    for i in range(n_elements):
        val = pixels[i]
        bin1_e = bin1[i]
        bin2_e = bin2[i]

        if not (bin1_e < bin2_e):
            continue

        if (bin1_e < r1) or (bin2_e > r2):
            continue

        if (bin2_e - bin1_e > lim) and (lim != -1):
            continue

        if math.isnan(weights[bin1_e]) or math.isnan(weights[bin2_e]):
            nanc += 1
            continue

        bal_val = val * weights[bin1_e] * weights[bin2_e]
        bal_val = 1 / bal_val

        # edges[count, 0] = bin1_e
        # edges[count, 1] = bin2_e
        # edges[count, 2] = bal_val
        # count += 1

        if bal_val < thresh:
            edges[count, 0] = bin1_e
            edges[count, 1] = bin2_e
            edges[count, 2] = bal_val
            count += 1

    return edges, count, nanc


def save_edges(h5, exp_name, thresh, output_path, r1, r2, weight_col, lim):

    # print("\nFinding edges for", exp_name)
    edge_file = output_path + exp_name + '-edges.csv'

    n_elements = h5['pixels']['bin1_id'].shape[0]
    weights = h5['bins'][weight_col][:]

    chunk_size = 50000

    edges = np.ones((500000000, 3))

    i_start = 0
    i_end = chunk_size

    count = 0
    nantotal = 0
    while i_start < n_elements:
        if i_end > n_elements:
            i_end = n_elements

        # print(f'process: {i_start / n_elements:.2%}, count = {count}', end='\r')

        bin1 = h5['pixels']['bin1_id'][i_start:i_end]
        bin2 = h5['pixels']['bin2_id'][i_start:i_end]
        pixels = h5['pixels']['count'][i_start:i_end]

        edges, count, nanc = search_edges(i_end - i_start, bin1, bin2, pixels, weights, thresh, edges, count,
                                          r1, r2, lim)

        i_start += chunk_size
        i_end += chunk_size
        nantotal += nanc

    # print("\nSaving edges", count)
    # print(f"Nan percentage: {nantotal / count:.2%}")
    edges = edges[:count]
    np.savetxt(edge_file, edges, delimiter=',', fmt='%.3f')

    return edge_file, count, nantotal


def hic_edges(cool_file, exp_name, weight_col, find_thresh_output, range_start, range_end, output_path,
              perc, perc_dist, edge_lim):

    logging.basicConfig(filename=exp_name + '.log', level=logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s:%(message)s')
    logger = logging.getLogger(__name__)

    range1 = range_start
    range2 = range_end

    start_time = time.time()

    logger.info("Step2: Finding edges.")

    perc_data = pd.read_csv(find_thresh_output, header=0, index_col=0)
    thresh = perc_data[perc][perc_dist]

    logger.info("{}: {} quantile of the estimated distance of all bins of distance {}".format(
        thresh, perc, perc_dist
    ))
    logger.info("Edges of distance above {} will be abandoned".format(thresh))
    logger.info("Hi-C edges found ranging from {} to {}".format(range1, range2))

    h5 = h5py.File(cool_file, 'r')
    edge_file, count, nan_count = save_edges(h5, exp_name, thresh, output_path, range1, range2, weight_col,
                                             edge_lim)

    logger.info("Save {} edges.".format(count))
    logger.info("{} NaN not included.".format(nan_count))

    end_time = time.time()
    logger.info("Step2 completed! Time taken(s): {:.2}".format(end_time - start_time))

    return edge_file, thresh
