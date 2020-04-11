import sys
import os
import numpy as np
from tabulate import tabulate



# parameters ########################################

# command line arguments
L_list = [int(L) for L in sys.argv[1][1:-1].split(",")] # list of considered sequence lengths



# directory and file structures ########################################

# input directory names
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name_1 = current_directory_name_oneback + "/7_fRNAdb_NNSE_NN_size_estimation/results/"
input_directory_name_2 = current_directory_name_oneback + "/8_fRNAdb_NNSE_NN_robustness_estimation/results/"

# results directory name
results_directory_name = "results"
results_directory_name += "/"

# create results directory (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass



# program ########################################

# results lists
results_L_list = []
results_number_not_converging_list = []

# go through all sequence lengths:
for L in L_list:

    ## 1. NN size estimates

    # results lists
    results_input_ID_list = []
    results_input_sequence_list = []
    results_input_pred_structure_list = []
    results_NNSE_NN_size_est_list = []

    # load results
    section_number = 1
    while os.path.exists(input_directory_name_1 + "L.%s/" % (L) + "L.%s_fRNAdb_NNSE_NN_size_estimates_section.%s.txt" % (L,section_number))==True:
        input_file_name_1 = input_directory_name_1 + "L.%s/" % (L) + "L.%s_fRNAdb_NNSE_NN_size_estimates_section.%s.txt" % (L,section_number)
        results_input_ID_list.extend(list(np.loadtxt(input_file_name_1, usecols=(0,), dtype=str, skiprows=2, unpack=True)))
        results_input_sequence_list.extend(list(np.loadtxt(input_file_name_1, usecols=(1,), dtype=str, skiprows=2, unpack=True)))
        results_input_pred_structure_list.extend(list(np.loadtxt(input_file_name_1, usecols=(2,), dtype=str, skiprows=2, unpack=True)))
        results_NNSE_NN_size_est_list.extend(list(np.loadtxt(input_file_name_1, usecols=(3,), dtype=str, skiprows=2, unpack=True)))
        section_number += 1

    # write results to .txt file
    table = [[results_input_ID_list[i],results_input_sequence_list[i],results_input_pred_structure_list[i],results_NNSE_NN_size_est_list[i]] for i,x in enumerate(results_input_ID_list)]
    file_name = results_directory_name + "L.%s_fRNAdb_NNSE_NN_size_estimates.txt" % (L)
    f = open(file_name,"w")
    f.write(tabulate(table,headers=["ID","sequence","pred. structure","NNSE: NN size est."]))
    f.close()

    # count number of sequences for which NNSE does not converge
    number_not_converging = results_NNSE_NN_size_est_list.count("n.c.")
    results_L_list.append(L)
    results_number_not_converging_list.append(number_not_converging)


    ## 2. NN robustness estimates

    # results lists
    results_input_ID_list = []
    results_input_sequence_list = []
    results_input_pred_structure_list = []
    results_NNSE_NN_rob_est_list = []

    # load results
    section_number = 1
    while os.path.exists(input_directory_name_2 + "L.%s/" % (L) + "L.%s_fRNAdb_NNSE_NN_rob_estimates_section.%s.txt" % (L,section_number))==True:
        input_file_name_2 = input_directory_name_2 + "L.%s/" % (L) + "L.%s_fRNAdb_NNSE_NN_rob_estimates_section.%s.txt" % (L,section_number)
        results_input_ID_list.extend(list(np.loadtxt(input_file_name_2, usecols=(0,), dtype=str, skiprows=2, unpack=True)))
        results_input_sequence_list.extend(list(np.loadtxt(input_file_name_2, usecols=(1,), dtype=str, skiprows=2, unpack=True)))
        results_input_pred_structure_list.extend(list(np.loadtxt(input_file_name_2, usecols=(2,), dtype=str, skiprows=2, unpack=True)))
        results_NNSE_NN_rob_est_list.extend(list(np.loadtxt(input_file_name_2, usecols=(3,), skiprows=2, unpack=True)))
        section_number += 1

    # write results to .txt file
    table = [[results_input_ID_list[i],results_input_sequence_list[i],results_input_pred_structure_list[i],results_NNSE_NN_rob_est_list[i]] for i,x in enumerate(results_input_ID_list)]
    file_name = results_directory_name + "L.%s_fRNAdb_NNSE_NN_rob_estimates.txt" % (L)
    f = open(file_name,"w")
    f.write(tabulate(table,headers=["ID","sequence","pred. structure","NNSE: NN rob. est."]))
    f.close()



# write overall results to .txt file ########################################

table = [[results_L_list[i],results_number_not_converging_list[i]] for i,x in enumerate(results_L_list)]
file_name = results_directory_name + "NNSE_NN_size_estimates_not_converging_summary.txt"
f = open(file_name,"w")
f.write(tabulate(table,headers=["L","# seq. not converging"]))
f.close()



print "done!"

quit()
