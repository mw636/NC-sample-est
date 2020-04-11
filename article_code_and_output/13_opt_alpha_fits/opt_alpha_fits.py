import sys
import os
import numpy as np
from scipy.optimize import curve_fit
from tabulate import tabulate



# parameters ########################################

# command line arguments
exhaustive_L_list = [int(L) for L in sys.argv[1][1:-1].split(",")] # list of sequence lengths for which optimal alpha is considered from exhaustive GP map analysis
fRNAdb_L_list = [int(L) for L in sys.argv[2][1:-1].split(",")] # list of fRNAdb sequence lengths for which optimal alpha is considered
fRNAdb_sample_size_list = [int(sample_size) for sample_size in sys.argv[3][1:-1].split(",")] # list of considered sample sizes for fRNAdb alpha optimisation
fRNAdb_random_subsample_size_list = [int(sample_size) for sample_size in sys.argv[4][1:-1].split(",")] # list of considered random subsample sizes for fRNAdb alpha optimisation
average_fRNAdb_sample_size_list = [int(sample_size) for sample_size in sys.argv[5][1:-1].split(",")] # list of considered sample sizes from fRNAdb used for fit results averaging
average_fRNAdb_random_subsample_size_list = [int(sample_size) for sample_size in sys.argv[6][1:-1].split(",")] # list of considered random subsample sizes from fRNAdb used for fit results averaging

# fixed
maxfev = 1000000 # maximum number of function calls (for fit)



# directory and file structures ########################################

# input directory names
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name_1 = current_directory_name_oneback + "/3_NC_size_estimation_validation/results/"
input_directory_name_2 = current_directory_name_oneback + "/12_fRNAdb_NN_size_estimates_comparison_optimisation/results/"

# results directory name
results_directory_name = "results/"

# create results directory (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass



# load input data ########################################

# optimal alpha from exhaustive GP map analysis
exhaustive_alpha_opt_list = []
for L in exhaustive_L_list:
    input_file_name_1 = input_directory_name_1 + "L.%s/" % (L) + "L.%s_neutral_component_size_estimation_validation_results.txt" % (L)
    alpha_list = list(np.loadtxt(input_file_name_1, usecols=(0,), skiprows=2, unpack=True))
    opt_check_list = list(np.loadtxt(input_file_name_1, usecols=(4,), dtype=str, skiprows=2, unpack=True))
    if opt_check_list[-1]=="True":
        alpha_opt = alpha_list[-1]
    else:
        print "error: optimal alpha not found in input"
        quit()
    exhaustive_alpha_opt_list.append(alpha_opt)

# optimal alpha from fRNAdb sampling
input_file_name_2 = input_directory_name_2 + "fRNAdb_NC_sampling_optimised_neutral_network_size_estimates_comparison_results.txt"
input_fRNAdb_L_list = list(np.loadtxt(input_file_name_2, usecols=(0,), dtype=int, skiprows=2, unpack=True))
input_fRNAdb_sample_size_list = list(np.loadtxt(input_file_name_2, usecols=(1,), dtype=int, skiprows=2, unpack=True))
input_fRNAdb_random_subsample_size_list = list(np.loadtxt(input_file_name_2, usecols=(2,), dtype=int, skiprows=2, unpack=True))
input_fRNAdb_alpha_opt_list = list(np.loadtxt(input_file_name_2, usecols=(3,), skiprows=2, unpack=True))



# program ########################################

# results lists
results_sample_size_list = []
results_random_subsample_size_list = []
results_fit_A_list = []
results_fit_A_SD_list = []
results_fit_B_list = []
results_fit_B_SD_list = []

# fit function
def fit_function(x,A,B):
    return A*(1-np.exp(-B*x))

# go through fRNAdb sample sizes
for sample_size in fRNAdb_sample_size_list:

    # go through fRNAdb random subsample sizes
    for random_subsample_size in fRNAdb_random_subsample_size_list:
    
        if sample_size>=random_subsample_size:

            # filtering fRNAdb sampling input
            f_input_fRNAdb_L_list = []
            f_input_fRNAdb_alpha_opt_list = []
            for i,fRNAdb_L in enumerate(input_fRNAdb_L_list):
                if input_fRNAdb_L_list[i] in fRNAdb_L_list and input_fRNAdb_sample_size_list[i]==sample_size and input_fRNAdb_random_subsample_size_list[i]==random_subsample_size:
                    f_input_fRNAdb_L_list.append(input_fRNAdb_L_list[i])
                    f_input_fRNAdb_alpha_opt_list.append(input_fRNAdb_alpha_opt_list[i])
        
            # data considered for fit
            fit_L_list = exhaustive_L_list + f_input_fRNAdb_L_list
            fit_alpha_opt_list = exhaustive_alpha_opt_list + f_input_fRNAdb_alpha_opt_list
            
            # fit
            popt, pcov = curve_fit(fit_function,fit_L_list,fit_alpha_opt_list,p0=[1,1],maxfev=maxfev)
            fit_A = popt[0]
            fit_B = popt[1]
            
            # standard deviation of fit parameters
            perr = np.sqrt(np.diag(pcov))
            fit_A_SD = perr[0]
            fit_B_SD = perr[1]
            
            results_sample_size_list.append(sample_size)
            results_random_subsample_size_list.append(random_subsample_size)
            results_fit_A_list.append(fit_A)
            results_fit_A_SD_list.append(fit_A_SD)
            results_fit_B_list.append(fit_B)
            results_fit_B_SD_list.append(fit_B_SD)
            
# averaging of fit results
avg_fit_A = 0
avg_fit_B = 0
counter = 0
for sample_size in average_fRNAdb_sample_size_list:
    for random_subsample_size in average_fRNAdb_random_subsample_size_list:
        if sample_size>=random_subsample_size:
            for i,x in enumerate(results_sample_size_list):
                if results_sample_size_list[i]==sample_size and results_random_subsample_size_list[i]==random_subsample_size:
                    avg_fit_A += results_fit_A_list[i]
                    avg_fit_B += results_fit_B_list[i]
                    counter += 1
avg_fit_A = round(float(avg_fit_A)/counter,2)
avg_fit_B = round(float(avg_fit_B)/counter,3)



# write results ########################################

# write optimised alpha fit results to .txt file
table = [[results_sample_size_list[i],results_random_subsample_size_list[i],results_fit_A_list[i],results_fit_A_SD_list[i],results_fit_B_list[i],results_fit_B_SD_list[i]] for i,x in enumerate(results_sample_size_list)]
file_name = results_directory_name + "optimised_alpha_fit_results.txt"
f = open(file_name,"w")
f.write(tabulate(table,headers=["fRNAdb sample size","fRNAdb random subsample size","fit A","fit A SD","fit B","fit B SD"]))
f.close()

# write averaged alpha fit results to .txt file
table = [[avg_fit_A,avg_fit_B]]
file_name = results_directory_name + "averaged_optimised_alpha_fit_results.txt"
f = open(file_name,"w")
f.write(tabulate(table,headers=["avg. fit A","avg. fit B"]))
f.close()



print "done!"

quit()
