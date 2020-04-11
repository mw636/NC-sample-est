import sys
import os
import numpy as np
import subprocess
from tabulate import tabulate



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length
input_list_section_total_number = int(sys.argv[2]) # total number of sections the input list is split
input_list_section_number = int(sys.argv[3]) # section number considered in this run (starting from 1)

# fixed
number_measurements = 1 # number of measurements with NNSE per input



# directory and file structures ########################################

# input directory name
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name = current_directory_name_oneback + "/6_fRNAdb_input_filtering/results/"

# results directory name
results_directory_name = "results"
results_directory_name += "/L.%s" % (L)
results_directory_name += "/"

# create results directory and subdirectory (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass
try:
    os.makedirs(results_directory_name + "NNSE_output/")
except OSError:
    pass



# program ########################################

# load input
input_file_name = input_directory_name + "L.%s_filtered_fRNAdb_input.txt" % (L)
input_IDs = list(np.loadtxt(input_file_name, usecols=(0,), dtype=str, skiprows=2, unpack=True))
input_sequences = list(np.loadtxt(input_file_name, usecols=(1,), dtype=str, skiprows=2, unpack=True))
input_pred_structures = list(np.loadtxt(input_file_name, usecols=(2,), dtype=str, skiprows=2, unpack=True))

# restrict to considered input list section
input_list_section_length = len(input_IDs) / input_list_section_total_number
input_list_section_min_index = (input_list_section_number-1)*input_list_section_length
input_list_section_max_index = (input_list_section_number)*input_list_section_length
if input_list_section_number==input_list_section_total_number:
    input_list_section_max_index=len(input_IDs) # last section can be longer if len(input_IDs) / input_list_section_total_number not an integer

# considered input
input_IDs = input_IDs[input_list_section_min_index:input_list_section_max_index]
input_sequences = input_sequences[input_list_section_min_index:input_list_section_max_index]
input_pred_structures = input_pred_structures[input_list_section_min_index:input_list_section_max_index]

# results lists
results_input_ID_list = []
results_input_sequence_list = []
results_input_pred_structure_list = []
results_NNSE_NN_size_est_list = []

# go through input
for index,input_ID in enumerate(input_IDs):

    # execute program 'get_size_NN'
    NNSE_output_file_name = results_directory_name + "NNSE_output/" + "NNSE_output_%s.txt" % (input_IDs[index])
    subprocess.call("./get_size_NN -s '%s' -m %s > '%s'" % (input_pred_structures[index],number_measurements,NNSE_output_file_name),shell=True)
    
    # read output file (three lines below would need to be updated if number_measurements!=1)
    NNSE_NN_size_est = float(np.loadtxt(NNSE_output_file_name, usecols=(4,), skiprows=6+L+1, unpack=True))
    if NNSE_NN_size_est==0: # not converging NNSE
        NNSE_NN_size_est="n.c."

    results_input_ID_list.append(input_IDs[index])
    results_input_sequence_list.append(input_sequences[index])
    results_input_pred_structure_list.append(input_pred_structures[index])
    results_NNSE_NN_size_est_list.append(NNSE_NN_size_est)



# write results to .txt file ########################################

table = [[results_input_ID_list[i],results_input_sequence_list[i],results_input_pred_structure_list[i],results_NNSE_NN_size_est_list[i]] for i,x in enumerate(results_input_ID_list)]
file_name = results_directory_name + "L.%s_fRNAdb_NNSE_NN_size_estimates_section.%s.txt" % (L,input_list_section_number)
f = open(file_name,"w")
f.write(tabulate(table,headers=["ID","sequence","pred. structure","NNSE: NN size est."]))
f.close()



print "done!"

quit()
