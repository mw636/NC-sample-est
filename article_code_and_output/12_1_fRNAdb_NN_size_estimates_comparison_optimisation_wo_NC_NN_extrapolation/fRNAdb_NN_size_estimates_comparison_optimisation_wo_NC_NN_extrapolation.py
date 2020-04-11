import sys
import os
import numpy as np
from scipy.optimize import minimize
from tabulate import tabulate



# parameters ########################################

# command line arguments
L_list = [int(L) for L in sys.argv[1][1:-1].split(",")] # list of considered sequence lengths
sample_size_list = [int(sample_size) for sample_size in sys.argv[2][1:-1].split(",")] # list of considered sample sizes
random_subsample_size_list = [int(sample_size) for sample_size in sys.argv[3][1:-1].split(",")] # list of considered random subsample sizes
alpha_list = [int(alpha) for alpha in sys.argv[4][1:-1].split(",")] # list of considered alpha values

# fixed
opt_method = 'Nelder-Mead'
opt_tol = 1e-8



# directory and file structures ########################################

# input directory names
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name_1 = current_directory_name_oneback + "/9_fRNAdb_NNSE_NN_estimates_processing/results/"
input_directory_name_2 = current_directory_name_oneback + "/11_fRNAdb_NC_sampling_processing/results/processed_fRNAdb_NC_sampling/"

# results directory name
results_directory_name = "results/"

# create results directory (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass



# help functions ########################################

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
        NN_size_est = NC_size_est # ! no NC-NN extrapolation
        NN_size_ests.append(NN_size_est)
    return NC_size_ests, NN_size_ests

# help function: error calculation
def calc_error(input_ID_list,NN_size_ests,NNSE_input_ID_list,NNSE_NN_size_ests):
    
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



# program ########################################

# results lists
results_L_list = []
results_sample_size_list = []
results_random_subsample_size_list = []
results_alpha_list = []
results_OE_average_list = []
results_OE_SD_list = []
results_OE_RMSD_list = []
results_OE_rel_RMSD_list = []
results_opt_check_list = []

# go through all sequence lengths:
for L in L_list:

    # load input
    
    # NNSE NN size estimates
    input_file_name_1 = input_directory_name_1 + "L.%s_fRNAdb_NNSE_NN_size_estimates.txt" % (L)
    NNSE_input_ID_list = list(np.loadtxt(input_file_name_1, usecols=(0,), dtype=str, skiprows=2, unpack=True))
    NNSE_input_sequence_list = list(np.loadtxt(input_file_name_1, usecols=(1,), dtype=str, skiprows=2, unpack=True))
    NNSE_input_pred_structure_list = list(np.loadtxt(input_file_name_1, usecols=(2,), dtype=str, skiprows=2, unpack=True))
    NNSE_NN_size_est_list = list(np.loadtxt(input_file_name_1, usecols=(3,), dtype=str, skiprows=2, unpack=True))

    # NC estimates from samples
    input_file_name_2 = input_directory_name_2 + "L.%s_fRNAdb_NC_sampling_sample_average_neutral_mutations_per_site.txt" % (L)
    input_ID_list = list(np.loadtxt(input_file_name_2, usecols=(0,), dtype=str, skiprows=2, unpack=True))
    input_sequence_list = list(np.loadtxt(input_file_name_2, usecols=(1,), dtype=str, skiprows=2, unpack=True))
    input_pred_structure_list = list(np.loadtxt(input_file_name_2, usecols=(2,), dtype=str, skiprows=2, unpack=True))
    input_sample_size_list = list(np.loadtxt(input_file_name_2, usecols=(3,), dtype=int, skiprows=2, unpack=True))
    input_random_subsample_size_list = list(np.loadtxt(input_file_name_2, usecols=(4,), dtype=int, skiprows=2, unpack=True))

    random_subsample_avg_neutral_mut_per_site_list = [[] for x in list(np.loadtxt(input_file_name_2, usecols=(7+0,), skiprows=2, unpack=True))]
    for l in range(L):
        load = list(np.loadtxt(input_file_name_2, usecols=(5+l,), skiprows=2, unpack=True))
        for n in range(len(load)):
            random_subsample_avg_neutral_mut_per_site_list[n].append(load[n])

    input_file_name_3 = input_directory_name_2 + "L.%s_fRNAdb_NC_sampling_sample_SD_neutral_mutations_per_site.txt" % (L)

    random_subsample_SD_neutral_mut_per_site_list = [[] for x in list(np.loadtxt(input_file_name_3, usecols=(7+0,), skiprows=2, unpack=True))]
    for l in range(L):
        load = list(np.loadtxt(input_file_name_3, usecols=(5+l,), skiprows=2, unpack=True))
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


                ## go through all alpha values
                for alpha in alpha_list:
                
                    NC_size_ests, NN_size_ests = NC_size_and_NN_size_estimation(f_input_ID_list,f_input_pred_structure_list,f_random_subsample_avg_neutral_mut_per_site_list,f_random_subsample_SD_neutral_mut_per_site_list,L,alpha)
                    OE_average, OE_SD, OE_RMSD, OE_rel_RMSD = calc_error(f_input_ID_list,NN_size_ests,NNSE_input_ID_list,NNSE_NN_size_est_list)

                    results_L_list.append(L)
                    results_sample_size_list.append(sample_size)
                    results_random_subsample_size_list.append(random_subsample_size)
                    results_alpha_list.append(alpha)
                    results_OE_average_list.append(OE_average)
                    results_OE_SD_list.append(OE_SD)
                    results_OE_RMSD_list.append(OE_RMSD)
                    results_OE_rel_RMSD_list.append(OE_rel_RMSD)
                    results_opt_check_list.append(False)
                

                ## alpha optimisation (via OE_RMSD)

                # function of which the result is optimised (minimised)
                def g_OE_RMSD(alpha,f_input_ID_list,f_input_pred_structure_list,f_random_subsample_avg_neutral_mut_per_site_list,f_random_subsample_SD_neutral_mut_per_site_list,NNSE_NN_size_est_list,L):

                    NC_size_ests, NN_size_ests = NC_size_and_NN_size_estimation(f_input_ID_list,f_input_pred_structure_list,f_random_subsample_avg_neutral_mut_per_site_list,f_random_subsample_SD_neutral_mut_per_site_list,L,alpha)
                    OE_average, OE_SD, OE_RMSD, OE_rel_RMSD = calc_error(f_input_ID_list,NN_size_ests,NNSE_input_ID_list,NNSE_NN_size_est_list)

                    return OE_RMSD

                # optimisation
                res = minimize(g_OE_RMSD, 1.0, args=(f_input_ID_list,f_input_pred_structure_list,f_random_subsample_avg_neutral_mut_per_site_list,f_random_subsample_SD_neutral_mut_per_site_list,NNSE_NN_size_est_list,L), method=opt_method, tol=opt_tol) # initialise optimisation with parameters = 1.0
                alpha_opt = res.x[0]

                # results for optimal alpha
                NC_size_ests_opt, NN_size_ests_opt = NC_size_and_NN_size_estimation(f_input_ID_list,f_input_pred_structure_list,f_random_subsample_avg_neutral_mut_per_site_list,f_random_subsample_SD_neutral_mut_per_site_list,L,alpha_opt)
                OE_average, OE_SD, OE_RMSD_opt, OE_rel_RMSD = calc_error(f_input_ID_list,NN_size_ests_opt,NNSE_input_ID_list,NNSE_NN_size_est_list)

                results_L_list.append(L)
                results_sample_size_list.append(sample_size)
                results_random_subsample_size_list.append(random_subsample_size)
                results_alpha_list.append(alpha_opt)
                results_OE_average_list.append(OE_average)
                results_OE_SD_list.append(OE_SD)
                results_OE_RMSD_list.append(OE_RMSD_opt)
                results_OE_rel_RMSD_list.append(OE_rel_RMSD)
                results_opt_check_list.append(True)

                # write optimised NC size and extrapolated NN size estimates to .txt file
                table = [[f_input_ID_list[i],f_input_sequence_list[i],f_input_pred_structure_list[i],NC_size_ests_opt[i],NN_size_ests_opt[i]] for i,x in enumerate(f_input_ID_list)]
                file_name = results_directory_name + "L.%s_fRNAdb_NC_sampling_sample_size.%s_random_subsample_size.%s_optimised_neutral_component_and_network_size_estimates_wo_NC_NN_extrapolation.txt" % (L,sample_size,random_subsample_size)
                f = open(file_name,"w")
                f.write(tabulate(table,headers=["ID","sequence","pred. structure","opt. NC size est.","opt. NN size est."]))
                f.close()



# write overall results ########################################

# write fRNAdb NC sampling optimised NN size estimates comparison results to .txt file
table = [[results_L_list[i],
         results_sample_size_list[i],
         results_random_subsample_size_list[i],
         results_alpha_list[i],
         results_OE_average_list[i],
         results_OE_SD_list[i],
         results_OE_RMSD_list[i],
         results_OE_rel_RMSD_list[i],
         results_opt_check_list[i]] for i,x in enumerate(results_L_list)]
file_name = results_directory_name + "fRNAdb_NC_sampling_optimised_neutral_network_size_estimates_comparison_results_wo_NC_NN_extrapolation.txt"
f = open(file_name,"w")
f.write(tabulate(table,headers=["L","sample size","random subsample size","alpha","OE average","OE SD","OE RMSD","OE rel. RMSD","optimised?"]))
f.close()



print "done!"

quit()
