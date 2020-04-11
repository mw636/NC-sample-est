import sys
import os
import numpy as np
from scipy.optimize import minimize
from tabulate import tabulate



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length

# fixed
alpha_list = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
opt_method = 'Nelder-Mead'
opt_tol = 1e-8 



# directory and file structures ########################################

# input directory name
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name = current_directory_name_oneback + "/2_NC_analysis/step_3/results/L.%s/" % (L)

# results directory name
results_directory_name = "results"
results_directory_name += "/L.%s" % (L)
results_directory_name += "/"

# create results directory (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass



# load input data ########################################

# NC characteristics
input_file_name_1 = input_directory_name + "L.%s_neutral_component_characteristics.txt" % (L)

NC_ranks = np.loadtxt(input_file_name_1, usecols=(0,), dtype=int, skiprows=2, unpack=True)
NC_indices = np.loadtxt(input_file_name_1, usecols=(1,), dtype=int, skiprows=2, unpack=True)
NC_phen_indices = np.loadtxt(input_file_name_1, usecols=(2,), dtype=int, skiprows=2, unpack=True)
NC_phen_seqs = np.loadtxt(input_file_name_1, usecols=(3,), dtype=str, skiprows=2, unpack=True)
NC_sizes = np.loadtxt(input_file_name_1, usecols=(4,), dtype=int, skiprows=2, unpack=True)

# NC average neutral mutations per site
input_file_name_2 = input_directory_name + "L.%s_neutral_component_average_neutral_mutations_per_site.txt" % (L)

NC_average_neutral_mut_per_site_list = [[] for rank in NC_ranks]
for l in range(L):
    load = np.loadtxt(input_file_name_2, usecols=(4+l,), skiprows=2, unpack=True)
    for i,rank in enumerate(NC_ranks):
        NC_average_neutral_mut_per_site_list[i].append(load[i])

# NC SD neutral mutations per site
input_file_name_3 = input_directory_name + "L.%s_neutral_component_SD_neutral_mutations_per_site.txt" % (L)

NC_SD_neutral_mut_per_site_list = [[] for rank in NC_ranks]
for l in range(L):
    load = np.loadtxt(input_file_name_3, usecols=(4+l,), skiprows=2, unpack=True)
    for i,rank in enumerate(NC_ranks):
        NC_SD_neutral_mut_per_site_list[i].append(load[i])



# program ########################################

# help function: NC size estimation
def NC_size_estimation(NC_ranks,NC_phen_seqs,NC_average_neutral_mut_per_site_list,NC_SD_neutral_mut_per_site_list,L,alpha):
    NC_size_ests = []
    for i,NC_rank in enumerate(NC_ranks):
        NC_size_est = 1.0
        for l in range(L):
            if NC_phen_seqs[i][l]==".":
                if (1+NC_average_neutral_mut_per_site_list[i][l]+alpha*NC_SD_neutral_mut_per_site_list[i][l])>4:
                    NC_size_est = NC_size_est*4.0
                else:
                    NC_size_est = NC_size_est*(1+NC_average_neutral_mut_per_site_list[i][l]+alpha*NC_SD_neutral_mut_per_site_list[i][l])
            if NC_phen_seqs[i][l]!=".":
                if (1+NC_average_neutral_mut_per_site_list[i][l]+alpha*NC_SD_neutral_mut_per_site_list[i][l])>2:
                    NC_size_est = NC_size_est*2.0
                else:
                    NC_size_est = NC_size_est*(1+NC_average_neutral_mut_per_site_list[i][l]+alpha*NC_SD_neutral_mut_per_site_list[i][l])
        NC_size_ests.append(NC_size_est)
    return NC_size_ests

# help function: error calculation
def calc_error(NC_ranks,NC_sizes,NC_size_ests):
    
    OE_average = 0.0
    OE_RMSD = 0.0
    
    for i,NC_rank in enumerate(NC_ranks):
        OE_average += np.log(NC_size_ests[i])/np.log(10)-np.log(NC_sizes[i])/np.log(10)
        OE_RMSD += (np.log(NC_size_ests[i])/np.log(10)-np.log(NC_sizes[i])/np.log(10))**2
    
    OE_average = float(OE_average)/len(NC_ranks)
    OE_RMSD = np.sqrt(float(OE_RMSD)/len(NC_ranks))

    OE_SD = 0.0

    for i,NC_rank in enumerate(NC_ranks):
        OE_SD += ((np.log(NC_size_ests[i])/np.log(10)-np.log(NC_sizes[i])/np.log(10))-OE_average)**2

    OE_SD = np.sqrt(float(OE_SD)/len(NC_ranks))
    
    return OE_average, OE_SD, OE_RMSD

# results lists
results_alpha_list = []
results_OE_average_list = []
results_OE_SD_list = []
results_OE_RMSD_list = []
results_opt_check_list = []



## 1. go through all given alpha values ########################################

for alpha in alpha_list:

    NC_size_ests = NC_size_estimation(NC_ranks,NC_phen_seqs,NC_average_neutral_mut_per_site_list,NC_SD_neutral_mut_per_site_list,L,alpha)
    OE_average, OE_SD, OE_RMSD = calc_error(NC_ranks,NC_sizes,NC_size_ests)

    results_alpha_list.append(alpha)
    results_OE_average_list.append(OE_average)
    results_OE_SD_list.append(OE_SD)
    results_OE_RMSD_list.append(OE_RMSD)
    results_opt_check_list.append(False)



## 2. alpha optimisation (via OE_RMSD) ########################################

# function of which the result is optimised (minimised)
def g_OE_RMSD(alpha,NC_ranks,NC_phen_seqs,NC_average_neutral_mut_per_site_list,NC_SD_neutral_mut_per_site_list,NC_sizes,L):

    NC_size_ests = NC_size_estimation(NC_ranks,NC_phen_seqs,NC_average_neutral_mut_per_site_list,NC_SD_neutral_mut_per_site_list,L,alpha)
    OE_average, OE_SD, OE_RMSD = calc_error(NC_ranks,NC_sizes,NC_size_ests)

    return OE_RMSD

# optimisation
res = minimize(g_OE_RMSD, 1.0, args=(NC_ranks,NC_phen_seqs,NC_average_neutral_mut_per_site_list,NC_SD_neutral_mut_per_site_list,NC_sizes,L), method=opt_method, tol=opt_tol) # initialise optimisation with parameters = 1.0
alpha_opt = res.x[0]

# results for optimal alpha
NC_size_ests_opt = NC_size_estimation(NC_ranks,NC_phen_seqs,NC_average_neutral_mut_per_site_list,NC_SD_neutral_mut_per_site_list,L,alpha_opt)
OE_average_opt, OE_SD_opt, OE_RMSD_opt = calc_error(NC_ranks,NC_sizes,NC_size_ests_opt)

results_alpha_list.append(alpha_opt)
results_OE_average_list.append(OE_average_opt)
results_OE_SD_list.append(OE_SD_opt)
results_OE_RMSD_list.append(OE_RMSD_opt)
results_opt_check_list.append(True)



# write results ########################################

# write optimised NC size estimates to .txt file
table = [[NC_ranks[i],NC_indices[i],NC_phen_indices[i],NC_phen_seqs[i],NC_sizes[i],NC_size_ests_opt[i]] for i,x in enumerate(NC_ranks)]
file_name = results_directory_name + "L.%s_neutral_component_optimised_NC_size_estimates.txt" % (L)
f = open(file_name,"w")
f.write(tabulate(table,headers=["NC rank","NC index","NC phen. index","NC phen. seq.","NC size","opt. NC size est."]))
f.close()

# write optimised NC size estimation validation results to .txt file
table = [[results_alpha_list[i],results_OE_average_list[i],results_OE_SD_list[i],results_OE_RMSD_list[i],results_opt_check_list[i]] for i,x in enumerate(results_alpha_list)]
file_name = results_directory_name + "L.%s_neutral_component_size_estimation_validation_results.txt" % (L)
f = open(file_name,"w")
f.write(tabulate(table,headers=["alpha","OE average","OE SD","OE RMSD","optimised?"]))
f.close()



print "done!"

quit()
