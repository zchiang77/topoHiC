import find_thresh
import hic_edges
import compute_PH
import get_scale
import get_score
import cooler
import os
import time

parent_dir = ""
exp_list = []
weight_col = ""
chrom_list = []
resolution = 0
edge_limit_bp = 0

scale_dist = [4]
compute_end_perc = "25%"
compute_end_perc_dist = 50
scale_perc = "75%"
significance = 0

filetype = 2  # This indicates that the source file is in (vertex, vertex, edge) form
dim = 1  # PH computation up to dim
begin_thresh = 0
threads = 64

os.system("mkdir " + parent_dir + "topohic/")

for exp_name in exp_list:

    start_time = time.time()
    os.system("mkdir " + parent_dir + "topohic/" + exp_name)

    for chrom in chrom_list:

        cool_file = parent_dir + exp_name + "-5kb.cool"
        output_path = parent_dir + "topohic/" + exp_name + "/topohic-" + chrom + "/"

        cooler_obj = cooler.Cooler(cool_file)
        cooler_ext = cooler_obj.extent(chrom)
        bin_range = list(cooler_ext)
        bin_range[1] -= 1
        print("NOTE: USING BINS RANGING FROM {} TO {} !!".format(bin_range[0], bin_range[1]))

        os.system("mkdir " + output_path)

        if edge_limit_bp == -1:
            edge_limit = -1
        else:
            edge_limit = edge_limit_bp / resolution
        find_thresh_output = find_thresh.find_thresh(cool_file, exp_name, weight_col, bin_range[0], bin_range[1],
                                                     output_path, edge_limit)
        edge_output, end_thresh = hic_edges.hic_edges(cool_file, exp_name, weight_col, find_thresh_output,
                                                      bin_range[0], bin_range[1], output_path,
                                                      compute_end_perc, compute_end_perc_dist, edge_limit)

        for scd in scale_dist:
            sc = get_scale.get_scale(find_thresh_output, scale_perc, scd)
            directory_for_ph = output_path + exp_name + "-sc" + str(scd) + "/"
            os.system("mkdir " + directory_for_ph)
            cyc_file = compute_PH.compute_ph(edge_output, exp_name, directory_for_ph, filetype,
                                             dim, begin_thresh, threads, end_thresh, sc)
            get_score.get_score(exp_name, directory_for_ph, bin_range, resolution)

    end_time = time.time()

    print("==================  CURRENT PROCESS  ==================")
    print("topohic for {} takes {} seconds.".format(exp_name, end_time - start_time))
    print("=======================================================\n")




