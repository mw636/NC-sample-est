import sys
import os
import numpy as np
from tabulate import tabulate



# parameters ########################################

# command line arguments
L_list = [int(L) for L in sys.argv[1][1:-1].split(",")] # list of considered sequence lengths
sample_size_list = [int(sample_size) for sample_size in sys.argv[2][1:-1].split(",")] # list of considered sample sizes
random_subsample_size_list = [int(sample_size) for sample_size in sys.argv[3][1:-1].split(",")] # list of considered random subsample sizes

# fixed
alphabet = ["A","C","G","U"] # RNA alphabet



# directory and file structures ########################################

# input directory name
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name = current_directory_name_oneback + "/10_fRNAdb_NC_sampling/results/"

# results directory names
results_directory_name = "results"
results_directory_name_1 = results_directory_name + "/processed_fRNAdb_NC_sampling/"
results_directory_name_2 = results_directory_name + "/RNAfold_counting/"

# create results directories (if not exist)
try:
    os.makedirs(results_directory_name_1)
except OSError:
    pass
try:
    os.makedirs(results_directory_name_2)
except OSError:
    pass



# program ########################################

# (overall) results lists
overall_results_L_list = []
overall_results_sample_size_list = []
overall_results_random_subsample_size_list = []
overall_results_avg_RNAfold_counter_total_list = []
overall_results_avg_beta_s_list = []
overall_results_avg_gamma_list = []

# go through all sequence lengths:
for L in L_list:

    ## 1. NC estimates from samples
    
    # results lists
    results_input_ID_list = []
    results_input_sequence_list = []
    results_input_pred_structure_list = []
    results_sample_size_list = []
    results_random_subsample_size_list = []
    results_RNAfold_counter_sampling_list = []
    results_RNAfold_counter_nb_measurement_list = []
    results_random_subsample_avg_neutral_mut_per_site_list = []
    results_random_subsample_SD_neutral_mut_per_site_list = []

    # load results
    section_number = 1
    while os.path.exists(input_directory_name + "L.%s/" % (L) + "L.%s_fRNAdb_NC_sampling_sample_average_neutral_mutations_per_site_section.%s.txt" % (L,section_number))==True:
        input_file_name_1 = input_directory_name + "L.%s/" % (L) + "L.%s_fRNAdb_NC_sampling_sample_average_neutral_mutations_per_site_section.%s.txt" % (L,section_number)
        
        results_input_ID_list.extend(list(np.loadtxt(input_file_name_1, usecols=(0,), dtype=str, skiprows=2, unpack=True)))
        results_input_sequence_list.extend(list(np.loadtxt(input_file_name_1, usecols=(1,), dtype=str, skiprows=2, unpack=True)))
        results_input_pred_structure_list.extend(list(np.loadtxt(input_file_name_1, usecols=(2,), dtype=str, skiprows=2, unpack=True)))
        results_sample_size_list.extend(list(np.loadtxt(input_file_name_1, usecols=(3,), dtype=int, skiprows=2, unpack=True)))
        results_random_subsample_size_list.extend(list(np.loadtxt(input_file_name_1, usecols=(4,), dtype=int, skiprows=2, unpack=True)))
        results_RNAfold_counter_sampling_list.extend(list(np.loadtxt(input_file_name_1, usecols=(5,), dtype=int, skiprows=2, unpack=True)))
        results_RNAfold_counter_nb_measurement_list.extend(list(np.loadtxt(input_file_name_1, usecols=(6,), dtype=int, skiprows=2, unpack=True)))
        
        random_subsample_avg_neutral_mut_per_site_list = [[] for x in list(np.loadtxt(input_file_name_1, usecols=(7+0,), skiprows=2, unpack=True))]
        for l in range(L):
            load = list(np.loadtxt(input_file_name_1, usecols=(7+l,), dtype=str, skiprows=2, unpack=True))
            for n in range(len(load)):
                random_subsample_avg_neutral_mut_per_site_list[n].append(load[n])
        results_random_subsample_avg_neutral_mut_per_site_list.extend(random_subsample_avg_neutral_mut_per_site_list)
        
        section_number += 1

    section_number = 1
    while os.path.exists(input_directory_name + "L.%s/" % (L) + "L.%s_fRNAdb_NC_sampling_sample_SD_neutral_mutations_per_site_section.%s.txt" % (L,section_number))==True:
        input_file_name_2 = input_directory_name + "L.%s/" % (L) + "L.%s_fRNAdb_NC_sampling_sample_SD_neutral_mutations_per_site_section.%s.txt" % (L,section_number)
        
        random_subsample_SD_neutral_mut_per_site_list = [[] for x in list(np.loadtxt(input_file_name_2, usecols=(7+0,), skiprows=2, unpack=True))]
        for l in range(L):
            load = list(np.loadtxt(input_file_name_2, usecols=(7+l,), dtype=str, skiprows=2, unpack=True))
            for n in range(len(load)):
                random_subsample_SD_neutral_mut_per_site_list[n].append(load[n])
        results_random_subsample_SD_neutral_mut_per_site_list.extend(random_subsample_SD_neutral_mut_per_site_list)
        
        section_number += 1

    # write results to .txt file
    table = [[results_input_ID_list[i],
              results_input_sequence_list[i],
              results_input_pred_structure_list[i],
              results_sample_size_list[i],
              results_random_subsample_size_list[i]] for i,x in enumerate(results_input_ID_list)]
    for i,x in enumerate(results_input_ID_list):
        table[i].extend(results_random_subsample_avg_neutral_mut_per_site_list[i])
    file_name = results_directory_name_1 + "L.%s_fRNAdb_NC_sampling_sample_average_neutral_mutations_per_site.txt" % (L)
    f = open(file_name,"w")
    headers=["ID","sequence","pred. structure","sample size","random subsample size"] + [""]*L
    f.write(tabulate(table,headers=headers))
    f.close()

    table = [[results_input_ID_list[i],
              results_input_sequence_list[i],
              results_input_pred_structure_list[i],
              results_sample_size_list[i],
              results_random_subsample_size_list[i]] for i,x in enumerate(results_input_ID_list)]
    for i,x in enumerate(results_input_ID_list):
        table[i].extend(results_random_subsample_SD_neutral_mut_per_site_list[i])
    file_name = results_directory_name_1 + "L.%s_fRNAdb_NC_sampling_sample_SD_neutral_mutations_per_site.txt" % (L)
    f = open(file_name,"w")
    headers=["ID","sequence","pred. structure","sample size","random subsample size"] + [""]*L
    f.write(tabulate(table,headers=headers))
    f.close()
    
    ## 2. RNAfold counting analysis
    
    # results lists
    results_RNAfold_counter_total_list = []
    results_beta_s_list = []
    results_beta_r_list = []
    results_gamma_list = []
    
    # analysis
    for index,input_ID in enumerate(results_input_ID_list):
    
        RNAfold_counter_total = results_RNAfold_counter_sampling_list[index] + results_RNAfold_counter_nb_measurement_list[index]

        beta_s = float(results_RNAfold_counter_sampling_list[index])/results_sample_size_list[index]
        
        beta_r = float(results_RNAfold_counter_nb_measurement_list[index])/results_random_subsample_size_list[index]/(len(alphabet)-1)/L
        
        n_bp = results_input_pred_structure_list[index].count("(") # number of base pairs
        gamma = 3-float((1-beta_r)*(len(alphabet)-1)*L)/2/n_bp
        
        results_RNAfold_counter_total_list.append(RNAfold_counter_total)
        results_beta_s_list.append(beta_s)
        results_beta_r_list.append(beta_r)
        results_gamma_list.append(gamma)
        
    # write results to .txt file
    table = [[results_input_ID_list[i],
              results_input_sequence_list[i],
              results_input_pred_structure_list[i],
              results_sample_size_list[i],
              results_random_subsample_size_list[i],
              results_RNAfold_counter_sampling_list[i],
              results_RNAfold_counter_nb_measurement_list[i],
              results_RNAfold_counter_total_list[i],
              results_beta_s_list[i],
              results_beta_r_list[i],
              results_gamma_list[i]] for i,x in enumerate(results_input_ID_list)]
    file_name = results_directory_name_2 + "L.%s_fRNAdb_NC_sampling_RNAfold_counting_analysis.txt" % (L)
    f = open(file_name,"w")
    headers=["ID","sequence","pred. structure","sample size","random subsample size","# RNAfold: sampling","# RNAfold: nb meas.","# RNAfold: total","# RNAfold: beta_s","# RNAfold: beta_r","# RNAfold: gamma"]
    f.write(tabulate(table,headers=headers))
    f.close()

    ## 3. RNAfold counting analysis averaging

    for sample_size in sample_size_list:

        for random_subsample_size in random_subsample_size_list:
        
            if sample_size>=random_subsample_size:

                # filtering
                f_results_RNAfold_counter_total_list = []
                f_results_beta_s_list = []
                f_results_gamma_list = []
                for index,input_ID in enumerate(results_input_ID_list):
                    if results_sample_size_list[index]==sample_size and results_random_subsample_size_list[index]==random_subsample_size:
                        f_results_RNAfold_counter_total_list.append(results_RNAfold_counter_total_list[index])
                        f_results_beta_s_list.append(results_beta_s_list[index])
                        f_results_gamma_list.append(results_gamma_list[index])
                # averaging
                avg_RNAfold_counter_total = float(sum(f_results_RNAfold_counter_total_list))/len(f_results_RNAfold_counter_total_list)
                avg_beta_s = float(sum(f_results_beta_s_list))/len(f_results_beta_s_list)
                avg_gamma = float(sum(f_results_gamma_list))/len(f_results_gamma_list)

                overall_results_L_list.append(L)
                overall_results_sample_size_list.append(sample_size)
                overall_results_random_subsample_size_list.append(random_subsample_size)
                overall_results_avg_RNAfold_counter_total_list.append(avg_RNAfold_counter_total)
                overall_results_avg_beta_s_list.append(avg_beta_s)
                overall_results_avg_gamma_list.append(avg_gamma)



# write overall results to .txt file ########################################

table = [[overall_results_L_list[i],
          overall_results_sample_size_list[i],
          overall_results_random_subsample_size_list[i],
          overall_results_avg_RNAfold_counter_total_list[i],
          overall_results_avg_beta_s_list[i],
          overall_results_avg_gamma_list[i]] for i,x in enumerate(overall_results_L_list)]
file_name = results_directory_name_2 + "fRNAdb_NC_sampling_RNAfold_counting_overall_analysis.txt"
f = open(file_name,"w")
headers=["L","sample size","random subsample size","# RNAfold: avg. total","# RNAfold: avg. beta_s","# RNAfold: avg. gamma"]
f.write(tabulate(table,headers=headers))
f.close()



print "done!"

quit()
