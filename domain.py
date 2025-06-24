import get_domain_helper
import cooler
import os
import time
import pandas as pd
import numpy as np
from itertools import combinations

pair_merge = True
parent_dir = ""
exp_list = []
chrom_list = []
score_threshold = 0.5
utest_significance = 0
resolution = 5000
domain_limit = 50000
max_iteration = 0
scale_dist = [4]

if pair_merge:
    os.system("mkdir " + parent_dir + "topohic/merge/")
    exp_combination = list(combinations(exp_list, 2))
    print("==================  CURRENT PROCESS  ==================")
    for exp_comb in exp_combination:
        start_time = time.time()
        exp1 = exp_comb[0]
        exp2 = exp_comb[1]
        exp_example = cooler.Cooler(parent_dir + exp1 + "-5kb.cool")
        for chrom in chrom_list:
            cooler_ext = exp_example.extent(chrom)
            bin_range = list(cooler_ext)
            bin_range[1] -= 1
            for scd in scale_dist:
                directory1 = (parent_dir + "topohic/" + exp1 + "/topohic-" + chrom + "/" +
                              exp1 + "-sc" + str(scd) + "/score.csv")
                directory2 = (parent_dir + "topohic/" + exp2 + "/topohic-" + chrom + "/" +
                              exp2 + "-sc" + str(scd) + "/score.csv")
                score1 = pd.read_csv(directory1)
                score2 = pd.read_csv(directory2)
                y1 = np.array(score1['loop_enrichment'])
                y2 = np.array(score2['loop_enrichment'])
                score = pd.DataFrame({'bin_id': list(score1['bin_id']),
                                      'loop_enrichment': np.mean([y1, y2], axis=0).tolist()})
                directory_merge = parent_dir + "topohic/merge/" + exp1 + "_" + exp2 + "_" + chrom + "_sc" + str(scd)
                score.to_csv(directory_merge + "_score.csv", header=True, index=False)
                get_domain_helper.get_domain(directory_merge, score_threshold, bin_range, utest_significance,
                                             max_iteration, resolution, domain_limit, pair_merge)
        end_time = time.time()
        print("Finding TADs for {}&{} takes {} seconds.".format(exp1, exp2, end_time - start_time))
    print("=======================================================\n")
else:
    print("==================  CURRENT PROCESS  ==================")
    for exp_name in exp_list:
        start_time = time.time()
        exp_cooler = cooler.Cooler(parent_dir + exp_name + "-5kb.cool")
        for chrom in chrom_list:
            cooler_ext = exp_cooler.extent(chrom)
            bin_range = list(cooler_ext)
            bin_range[1] -= 1
            for scd in scale_dist:
                directory = (parent_dir + "topohic/" + exp_name + "/topohic-" + chrom + "/" +
                             exp_name + "-sc" + str(scd))
                get_domain_helper.get_domain(directory, score_threshold, bin_range, utest_significance,
                                             max_iteration, resolution, domain_limit, pair_merge)
        end_time = time.time()
        print("Finding TADs for {} takes {} seconds.".format(exp_name, end_time - start_time))
    print("=======================================================\n")





