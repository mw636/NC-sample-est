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

# input directory names
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name_1 = current_directory_name_oneback + "/13_opt_alpha_fits/results/"
input_directory_name_2 = current_directory_name_oneback + "/9_fRNAdb_NNSE_NN_estimates_processing/results/"
input_directory_name_3 = current_directory_name_oneback + "/11_fRNAdb_NC_sampling_processing/results/processed_fRNAdb_NC_sampling/"

# results directory name
results_directory_name = "results/"

# create results directory (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass



# help functions ########################################

# help function: alpha from fit function
input_file_name_1 = input_directory_name_1 + "averaged_optimised_alpha_fit_results.txt"
fit_A = float(np.loadtxt(input_file_name_1, usecols=(0,), skiprows=2, unpack=True))
fit_B = float(np.loadtxt(input_file_name_1, usecols=(1,), skiprows=2, unpack=True))
def alpha_L(fit_A,fit_B,L):
    return fit_A*(1-np.exp(-fit_B*L))

# help function: NC size estimation and extrapolated NN size estimation
def NC_size_and_NN_size_estimation(input_ID_list,input_pred_structure_list,random_subsample_avg_neutral_mut_per_site_list,random_subsample_SD_neutral_mut_per_site_list,L,alpha):
    NC_size_ests = []
    NN_size_ests = []
    for i,input_ID in enumerate(input_ID_list):
        NC_size_est = 1.0
        for l in range(L):
            if input_pred_structure_list[i][l]==".":
                if (1+random_subsample_avg_neutral_mut_per_site_list[i][l]+alpha*random_subsample_SD_neutral_mut_per_site_list[i][l])>4:
                    NC_size_est = NC_size_est*4.0
                else:
                    NC_size_est = NC_size_est*(1+random_subsample_avg_neutral_mut_per_site_list[i][l]+alpha*random_subsample_SD_neutral_mut_per_site_list[i][l])
            if input_pred_structure_list[i][l]!=".":
                if (1+random_subsample_avg_neutral_mut_per_site_list[i][l]+alpha*random_subsample_SD_neutral_mut_per_site_list[i][l])>2:
                    NC_size_est = NC_size_est*2.0
                else:
                    NC_size_est = NC_size_est*(1+random_subsample_avg_neutral_mut_per_site_list[i][l]+alpha*random_subsample_SD_neutral_mut_per_site_list[i][l])
        NC_size_ests.append(NC_size_est)
        NN_size_est = NC_size_est*2**(input_pred_structure_list[i].count("("))
        NN_size_ests.append(NN_size_est)
    return NC_size_ests, NN_size_ests

# help function: NC robustness estimation and extrapolated NN robustness estimation
def NC_rob_and_NN_rob_estimation(input_ID_list,input_pred_structure_list,random_subsample_avg_neutral_mut_per_site_list,alphabet,L):
    NC_rob_ests = []
    NN_rob_ests = []
    for i,input_ID in enumerate(input_ID_list):
        NC_rob_est = float(sum(random_subsample_avg_neutral_mut_per_site_list[i]))/(len(alphabet)-1)/L
        NC_rob_ests.append(NC_rob_est)
        NN_rob_est = NC_rob_est
        NN_rob_ests.append(NN_rob_est)
    return NC_rob_ests, NN_rob_ests

# help function: error calculation for extrapolated NN size estimates
def calc_error_NN_size_est(input_ID_list,NN_size_ests,NNSE_input_ID_list,NNSE_NN_size_ests):
    
    OE_average = 0.0
    OE_RMSD = 0.0
    OE_rel_RMSD = 0.0
    
    counter = 0
    
    for i,input_ID in enumerate(input_ID_list):
        # check if same input ID is considered
        if input_ID_list[i]!=NNSE_input_ID_list[i]:
            print "error: input IDs not matching"
            quit()
        else:
            if NNSE_NN_size_ests[i]!="n.c.":
                OE_average += np.log(NN_size_ests[i])/np.log(10)-np.log(float(NNSE_NN_size_ests[i]))/np.log(10)
                OE_RMSD += (np.log(NN_size_ests[i])/np.log(10)-np.log(float(NNSE_NN_size_ests[i]))/np.log(10))**2
                OE_rel_RMSD += (float(np.log(NN_size_ests[i])/np.log(10)-np.log(float(NNSE_NN_size_ests[i]))/np.log(10))/(np.log(float(NNSE_NN_size_ests[i]))/np.log(10)))**2
                
                counter += 1
    
    OE_average = float(OE_average)/counter
    OE_RMSD = np.sqrt(float(OE_RMSD)/counter)
    OE_rel_RMSD = np.sqrt(float(OE_rel_RMSD)/counter)

    OE_SD = 0.0

    for i,input_ID in enumerate(input_ID_list):
        if NNSE_NN_size_ests[i]!="n.c.":
            OE_SD += ((np.log(NN_size_ests[i])/np.log(10)-np.log(float(NNSE_NN_size_ests[i]))/np.log(10))-OE_average)**2

    OE_SD = np.sqrt(float(OE_SD)/counter)
    
    return OE_average, OE_SD, OE_RMSD, OE_rel_RMSD

# help function: error calculation for extrapolated NN robustness estimates
def calc_error_NN_rob_est(input_ID_list,NN_rob_ests,NNSE_input_ID_list,NNSE_NN_rob_ests):
    
    OE_average = 0.0
    OE_RMSD = 0.0
    
    counter = 0
    
    for i,NNSE_input_ID in enumerate(NNSE_input_ID_list): # NNSE estimates are restricting: only 100 NN robustness estimates per sequence length considered
    
        index = input_ID_list.index(NNSE_input_ID)
    
        OE_average += NN_rob_ests[index]-NNSE_NN_rob_ests[i]
        OE_RMSD += (NN_rob_ests[index]-NNSE_NN_rob_ests[i])**2
        
        counter += 1
    
    OE_average = float(OE_average)/counter
    OE_RMSD = np.sqrt(float(OE_RMSD)/counter)
    
    OE_SD = 0.0
    
    for i,NNSE_input_ID in enumerate(NNSE_input_ID_list):
    
        index = input_ID_list.index(NNSE_input_ID)
    
        OE_SD += ((NN_rob_ests[index]-NNSE_NN_rob_ests[i])-OE_average)**2

    OE_SD = np.sqrt(float(OE_SD)/counter)
    
    return OE_average, OE_SD, OE_RMSD



# program ########################################

# results lists
size_results_L_list = []
size_results_sample_size_list = []
size_results_random_subsample_size_list = []
size_results_OE_average_list = []
size_results_OE_SD_list = []
size_results_OE_RMSD_list = []
size_results_OE_rel_RMSD_list = []


rob_results_L_list = []
rob_results_sample_size_list = []
rob_results_random_subsample_size_list = []
rob_results_OE_average_list = []
rob_results_OE_SD_list = []
rob_results_OE_RMSD_list = []


# go through all sequence lengths:
for L in L_list:


    ## 1. NN and NN size estimates and error quantification
    
    # load input
    
    # NNSE NN size estimates
    input_file_name_2 = input_directory_name_2 + "L.%s_fRNAdb_NNSE_NN_size_estimates.txt" % (L)
    NNSE_input_ID_list = list(np.loadtxt(input_file_name_2, usecols=(0,), dtype=str, skiprows=2, unpack=True))
    NNSE_input_sequence_list = list(np.loadtxt(input_file_name_2, usecols=(1,), dtype=str, skiprows=2, unpack=True))
    NNSE_input_pred_structure_list = list(np.loadtxt(input_file_name_2, usecols=(2,), dtype=str, skiprows=2, unpack=True))
    NNSE_NN_size_est_list = list(np.loadtxt(input_file_name_2, usecols=(3,), dtype=str, skiprows=2, unpack=True))

    # NC estimates from samples
    input_file_name_3 = input_directory_name_3 + "L.%s_fRNAdb_NC_sampling_sample_average_neutral_mutations_per_site.txt" % (L)
    input_ID_list = list(np.loadtxt(input_file_name_3, usecols=(0,), dtype=str, skiprows=2, unpack=True))
    input_sequence_list = list(np.loadtxt(input_file_name_3, usecols=(1,), dtype=str, skiprows=2, unpack=True))
    input_pred_structure_list = list(np.loadtxt(input_file_name_3, usecols=(2,), dtype=str, skiprows=2, unpack=True))
    input_sample_size_list = list(np.loadtxt(input_file_name_3, usecols=(3,), dtype=int, skiprows=2, unpack=True))
    input_random_subsample_size_list = list(np.loadtxt(input_file_name_3, usecols=(4,), dtype=int, skiprows=2, unpack=True))

    random_subsample_avg_neutral_mut_per_site_list = [[] for x in list(np.loadtxt(input_file_name_3, usecols=(7+0,), skiprows=2, unpack=True))]
    for l in range(L):
        load = list(np.loadtxt(input_file_name_3, usecols=(5+l,), skiprows=2, unpack=True))
        for n in range(len(load)):
            random_subsample_avg_neutral_mut_per_site_list[n].append(load[n])

    input_file_name_4 = input_directory_name_3 + "L.%s_fRNAdb_NC_sampling_sample_SD_neutral_mutations_per_site.txt" % (L)

    random_subsample_SD_neutral_mut_per_site_list = [[] for x in list(np.loadtxt(input_file_name_4, usecols=(7+0,), skiprows=2, unpack=True))]
    for l in range(L):
        load = list(np.loadtxt(input_file_name_4, usecols=(5+l,), skiprows=2, unpack=True))
        for n in range(len(load)):
            random_subsample_SD_neutral_mut_per_site_list[n].append(load[n])

    # go through sample sizes
    for sample_size in sample_size_list:

        # go through random subsample sizes
        for random_subsample_size in random_subsample_size_list:
        
            if sample_size>=random_subsample_size:
 
                # filtering
                f_input_ID_list = []
                f_input_sequence_list = []
                f_input_pred_structure_list = []
                f_random_subsample_avg_neutral_mut_per_site_list = []
                f_random_subsample_SD_neutral_mut_per_site_list = []
                for i,input_ID in enumerate(input_ID_list):
                    if input_sample_size_list[i]==sample_size and input_random_subsample_size_list[i]==random_subsample_size:
                        f_input_ID_list.append(input_ID_list[i])
                        f_input_sequence_list.append(input_sequence_list[i])
                        f_input_pred_structure_list.append(input_pred_structure_list[i])
                        f_random_subsample_avg_neutral_mut_per_site_list.append(random_subsample_avg_neutral_mut_per_site_list[i])
                        f_random_subsample_SD_neutral_mut_per_site_list.append(random_subsample_SD_neutral_mut_per_site_list[i])

                # alpha from function
                alpha = alpha_L(fit_A,fit_B,L)

                # results
                NC_size_ests, NN_size_ests = NC_size_and_NN_size_estimation(f_input_ID_list,f_input_pred_structure_list,f_random_subsample_avg_neutral_mut_per_site_list,f_random_subsample_SD_neutral_mut_per_site_list,L,alpha)
                OE_average, OE_SD, OE_RMSD, OE_rel_RMSD = calc_error_NN_size_est(f_input_ID_list,NN_size_ests,NNSE_input_ID_list,NNSE_NN_size_est_list)

                size_results_L_list.append(L)
                size_results_sample_size_list.append(sample_size)
                size_results_random_subsample_size_list.append(random_subsample_size)
                size_results_OE_average_list.append(OE_average)
                size_results_OE_SD_list.append(OE_SD)
                size_results_OE_RMSD_list.append(OE_RMSD)
                size_results_OE_rel_RMSD_list.append(OE_rel_RMSD)

                # write NC size and extrapolated NN size estimates to .txt file
                table = [[f_input_ID_list[i],f_input_sequence_list[i],f_input_pred_structure_list[i],NC_size_ests[i],NN_size_ests[i]] for i,x in enumerate(f_input_ID_list)]
                file_name = results_directory_name + "L.%s_fRNAdb_NC_sampling_sample_size.%s_random_subsample_size.%s_neutral_component_and_network_size_estimates.txt" % (L,sample_size,random_subsample_size)
                f = open(file_name,"w")
                f.write(tabulate(table,headers=["ID","sequence","pred. structure","NC size est.","NN size est."]))
                f.close()


    ## 2. NN and NN robustness estimates and error quantification

    # load input

    # NNSE NN robustness estimates
    input_file_name_5 = input_directory_name_2 + "L.%s_fRNAdb_NNSE_NN_rob_estimates.txt" % (L)
    NNSE_input_ID_list = list(np.loadtxt(input_file_name_5, usecols=(0,), dtype=str, skiprows=2, unpack=True))
    NNSE_input_sequence_list = list(np.loadtxt(input_file_name_5, usecols=(1,), dtype=str, skiprows=2, unpack=True))
    NNSE_input_pred_structure_list = list(np.loadtxt(input_file_name_5, usecols=(2,), dtype=str, skiprows=2, unpack=True))
    NNSE_NN_rob_est_list = list(np.loadtxt(input_file_name_5, usecols=(3,), skiprows=2, unpack=True))

    # go through sample sizes
    for sample_size in sample_size_list:

        # go through random subsample sizes
        for random_subsample_size in random_subsample_size_list:
        
            if sample_size>=random_subsample_size:
 
                # filtering
                f_input_ID_list = []
                f_input_sequence_list = []
                f_input_pred_structure_list = []
                f_random_subsample_avg_neutral_mut_per_site_list = []
                f_random_subsample_SD_neutral_mut_per_site_list = []
                for i,input_ID in enumerate(input_ID_list):
                    if input_sample_size_list[i]==sample_size and input_random_subsample_size_list[i]==random_subsample_size:
                        f_input_ID_list.append(input_ID_list[i])
                        f_input_sequence_list.append(input_sequence_list[i])
                        f_input_pred_structure_list.append(input_pred_structure_list[i])
                        f_random_subsample_avg_neutral_mut_per_site_list.append(random_subsample_avg_neutral_mut_per_site_list[i])
                        f_random_subsample_SD_neutral_mut_per_site_list.append(random_subsample_SD_neutral_mut_per_site_list[i])

                # results
                NC_rob_ests, NN_rob_ests = NC_rob_and_NN_rob_estimation(f_input_ID_list,f_input_pred_structure_list,f_random_subsample_avg_neutral_mut_per_site_list,alphabet,L)
                OE_average, OE_SD, OE_RMSD = calc_error_NN_rob_est(f_input_ID_list,NN_rob_ests,NNSE_input_ID_list,NNSE_NN_rob_est_list)

                rob_results_L_list.append(L)
                rob_results_sample_size_list.append(sample_size)
                rob_results_random_subsample_size_list.append(random_subsample_size)
                rob_results_OE_average_list.append(OE_average)
                rob_results_OE_SD_list.append(OE_SD)
                rob_results_OE_RMSD_list.append(OE_RMSD)

                # write NC robustness and extrapolated NN robustness estimates to .txt file
                table = [[f_input_ID_list[i],f_input_sequence_list[i],f_input_pred_structure_list[i],NC_rob_ests[i],NN_rob_ests[i]] for i,x in enumerate(f_input_ID_list)]
                file_name = results_directory_name + "L.%s_fRNAdb_NC_sampling_sample_size.%s_random_subsample_size.%s_neutral_component_and_network_robustness_estimates.txt" % (L,sample_size,random_subsample_size)
                f = open(file_name,"w")
                f.write(tabulate(table,headers=["ID","sequence","pred. structure","NC rob. est.","NN rob. est."]))
                f.close()



# write overall results ########################################

# write fRNAdb NC sampling NN size estimates comparison results to .txt file
table = [[size_results_L_list[i],
          size_results_sample_size_list[i],
          size_results_random_subsample_size_list[i],
          size_results_OE_average_list[i],
          size_results_OE_SD_list[i],
          size_results_OE_RMSD_list[i],
          size_results_OE_rel_RMSD_list[i]] for i,x in enumerate(size_results_L_list)]
file_name = results_directory_name + "fRNAdb_NC_sampling_neutral_network_size_estimates_comparison_results.txt"
f = open(file_name,"w")
f.write(tabulate(table,headers=["L","sample size","random subsample size","OE average","OE SD","OE RMSD","OE rel. RMSD"]))
f.close()

# write fRNAdb NC sampling NN robustness estimates comparison results to .txt file
table = [[rob_results_L_list[i],
          rob_results_sample_size_list[i],
          rob_results_random_subsample_size_list[i],
          rob_results_OE_average_list[i],
          rob_results_OE_SD_list[i],
          rob_results_OE_RMSD_list[i]] for i,x in enumerate(rob_results_L_list)]
file_name = results_directory_name + "fRNAdb_NC_sampling_neutral_network_robustness_estimates_comparison_results.txt"
f = open(file_name,"w")
f.write(tabulate(table,headers=["L","sample size","random subsample size","OE average","OE SD","OE RMSD"]))
f.close()



print "done!"

quit()
