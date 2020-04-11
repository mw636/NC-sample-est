import random
import os
import numpy as np
from tabulate import tabulate

import RNA



# parameters ########################################

# fixed
L_list = [20,40,45,50,55,60,65,70,75,80,85,90,95,100] # sequence lengths
alphabet = ["A","C","G","U"] # RNA alphabet
random_seed = 1 # seed for random number generator
random.seed(random_seed)
number_seq_NNSE_rob_est = 100



# directory and file structures ########################################

# input directory name
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name = current_directory_name_oneback + "/fRNAdb_input/"

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
results_number_unfiltered_list = []
results_number_filtered_putative_list = []
results_number_filtered_piRNA_list = []
results_number_filtered_incompatible_list = []
results_number_final_list = []

# go through all given sequence lengths
for L in L_list:



    ## 1. load and process fRNAdb input file ########################################

    # load fRNAdb input file
    input_file_name = input_directory_name + "L.%s_frnadb_summary.csv" % (L)

    # csv file can have column entries with a , within " ", such that loading columns according to , separation is not feasible
    # process fRNAdb input file: replace , within " " to a .
    f = open(input_file_name,"rb")
    f_processed = []
    for line in f:
        line_processed = ''
        quotation_check = False
        for character in line:
            if character=='"':
                if quotation_check==False:
                    quotation_check = True
                else:
                    quotation_check = False
            if character==',' and quotation_check==True:
                line_processed += '.'
            else:
                line_processed += character
        f_processed.append(line_processed)

    # load fRNAdb input
    input_IDs = np.loadtxt(f_processed, dtype=str, delimiter=',', usecols=(0,), skiprows=1, unpack=True)
    input_accessions = np.loadtxt(f_processed, dtype=str, delimiter=',', usecols=(1,), skiprows=1, unpack=True)
    input_descriptions = np.loadtxt(f_processed, dtype=str, delimiter=',', usecols=(2,), skiprows=1, unpack=True)
    input_SO_names = np.loadtxt(f_processed, dtype=str, delimiter=',', usecols=(3,), skiprows=1, unpack=True)
    input_organisms = np.loadtxt(f_processed, dtype=str, delimiter=',', usecols=(4,), skiprows=1, unpack=True)
    input_xrefs = np.loadtxt(f_processed, dtype=str, delimiter=',', usecols=(5,), skiprows=1, unpack=True)
    input_lengths = np.loadtxt(f_processed, dtype=str, delimiter=',', usecols=(6,), skiprows=1, unpack=True)
    input_sequences = np.loadtxt(f_processed, dtype=str, delimiter=',', usecols=(7,), skiprows=1, unpack=True)

    # process fRNAdb input: replace T in sequences by U
    for i,input_ID in enumerate(input_IDs):
        sequence_processed = ""
        for character in input_sequences[i]:
            if character=="T":
                sequence_processed += "U"
            else:
                sequence_processed += character
        input_sequences[i] = sequence_processed



    ## 2. filter fRNAdb input ########################################

    filtered_input_IDs = []
    filtered_input_sequences = []
    filtered_input_pred_structures = []

    number_unfiltered = len(input_IDs)

    number_filtered_putative = 0
    number_filtered_piRNA = 0
    number_filtered_incompatible = 0

    for i,input_ID in enumerate(input_IDs):
    
        # putative?
        if "Putative" in input_descriptions[i] or "Putative" in input_SO_names[i] or "putative" in input_descriptions[i] or "putative" in input_SO_names[i]:
            number_filtered_putative += 1
        else:

            # piRNA?
            if "PiRNA" in input_descriptions[i] or "PiRNA" in input_SO_names[i] or "piRNA" in input_descriptions[i] or "piRNA" in input_SO_names[i]:
                number_filtered_piRNA += 1
            else:
            
                # compatible?
                compatible_check = True
                
                # 1. sequence only composed of letters from RNA alphabet?
                for character in input_sequences[i]:
                    if character not in alphabet:
                        compatible_check = False
                # 2. structure defined?
                if compatible_check==True:
                    (structure, mfe) = RNA.fold(input_sequences[i])
                    if structure=="."*L:
                        compatible_check = False
                
                if compatible_check==False:
                    number_filtered_incompatible += 1
                else:
                    
                    filtered_input_IDs.append(input_IDs[i])
                    filtered_input_sequences.append(input_sequences[i])
                    filtered_input_pred_structures.append(structure)
    
    number_final = len(filtered_input_IDs)

    results_L_list.append(L)
    results_number_unfiltered_list.append(number_unfiltered)
    results_number_filtered_putative_list.append(number_filtered_putative)
    results_number_filtered_piRNA_list.append(number_filtered_piRNA)
    results_number_filtered_incompatible_list.append(number_filtered_incompatible)
    results_number_final_list.append(number_final)

    # write filtered fRNAdb input to .txt file
    table = [[filtered_input_IDs[i],filtered_input_sequences[i],filtered_input_pred_structures[i]] for i,filtered_input_ID in enumerate(filtered_input_IDs)]
    file_name = results_directory_name + "L.%s_filtered_fRNAdb_input.txt" % (L)
    f = open(file_name,"w")
    f.write(tabulate(table,headers=["ID","sequence","pred. structure"]))
    f.close()



    ## 3. randomly choose filtered fRNAdb input considered for NNSE robustness estimation  ########################################

    # choose random indices
    random_indices = random.sample(range(len(filtered_input_IDs)),number_seq_NNSE_rob_est)

    # write filtered fRNAdb input for NNSE robustness estimation to .txt file
    table = [[filtered_input_IDs[i],filtered_input_sequences[i],filtered_input_pred_structures[i]] for i in random_indices]
    file_name = results_directory_name + "L.%s_filtered_fRNAdb_input_for_NNSE_robustness_estimation.txt" % (L)
    f = open(file_name,"w")
    f.write(tabulate(table,headers=["ID","sequence","pred. structure"]))
    f.close()



# write overall results to .txt file ########################################

table = [[results_L_list[i],results_number_unfiltered_list[i],-results_number_filtered_putative_list[i],-results_number_filtered_piRNA_list[i],-results_number_filtered_incompatible_list[i],results_number_final_list[i]] for i,x in enumerate(results_L_list)]
file_name = results_directory_name + "fRNAdb_input_filtering_summary.txt"
f = open(file_name,"w")
f.write(tabulate(table,headers=["L","# seq.","# putative","# piRNA","# incomp.","# final"]))
f.close()



print "done!"

quit()
