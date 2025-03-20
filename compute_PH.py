import pydory as dory
import check_connected_helper as cch
import helper_functions as hf
import networkx as nx
import pandas as pd
import time
import logging


def compute_ph(edge_file, exp_name, target, filetype, dim, begin_thresh, threads, end_thresh, scale):

    logging.basicConfig(filename=exp_name + '.log', level=logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s:%(message)s')
    logger = logging.getLogger(__name__)

    logger.info("Step3: (scale {}) Computing persistent homology.".format(scale))
    logger.info("Persistent homology computation ends at {}.".format(end_thresh))

    start_time = time.time()
    dory.compute_PH(edge_file, begin_thresh, end_thresh, filetype, threads, target, dim, 1, 1, scale, 1, 1000, 100)
    # dory.compute_PH(source, lower_thresh, thresh,
    #        filetype, threads, target, dim,
    #        compute_cycles, reduce_cyc_lengths, cyc_thresh,
    #        suppress_output, hom_batch_size, cohom_batch_size)

    logger.info("Start minimizing cycles using networkx(python) ...")

    # Return the basis of each cycle
    # that is the minimal collection of nodes that forms a cycle
    cyc_file = target + 'minimal_V_birth_H' + str(dim) + '.txt'
    cch.modify_minimal(cyc_file, target, dim)

    end_time = time.time()
    logger.info("Step3 (scale {}) completed! Time taken(s): {:.2}.".format(scale, end_time - start_time))

    return cyc_file
