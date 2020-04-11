import sys
import os
import numpy as np
from tabulate import tabulate



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length
max_NC_rank = int(sys.argv[2]) # maximum considered NC rank
sampling_method_list = [str(sampling_method) for sampling_method in sys.argv[3][1:-1].split(",")] # list of considered sampling methods
sample_size_list = [int(sample_size) for sample_size in sys.argv[4][1:-1].split(",")] # list of considered sample sizes
random_subsample_size_list = [int(sample_size) for sample_size in sys.argv[5][1:-1].split(",")] # list of considered random subsample sizes (not for random sampling)
number_samples = int(sys.argv[6]) # number of samples considered for each sampling method, sample size (combination), and NC, respectively

# fixed
alphabet = ["A","C","G","U"] # RNA alphabet



# directory and file structures ########################################

# input directory names
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name_1 = current_directory_name_oneback + "/2_NC_analysis/step_3/results/L.%s/" % (L)
input_directory_name_2 = current_directory_name_oneback + "/3_NC_size_estimation_validation/results/L.%s/" % (L)
input_directory_name_3 = current_directory_name_oneback + "/4_NC_sampling/results/L.%s/" % (L)

# results directory names
results_directory_name = "results"
results_directory_name += "/L.%s" % (L)
results_directory_name_1 = results_directory_name + "/NC_estimates_from_samples/"
results_directory_name_2 = results_directory_name + "/error_quantification/"

# create results directories (if not exist)
try:
    os.makedirs(results_directory_name_1)
except OSError:
    pass
try:
    os.makedirs(results_directory_name_2)
except OSError:
    pass



# load (general) input data ########################################

# NC characteristics
input_file_name_1 = input_directory_name_1 + "L.%s_neutral_component_characteristics.txt" % (L)

NC_ranks = list(np.loadtxt(input_file_name_1, usecols=(0,), dtype=int, skiprows=2, unpack=True))
NC_indices = list(np.loadtxt(input_file_name_1, usecols=(1,), dtype=int, skiprows=2, unpack=True))
NC_phen_indices = list(np.loadtxt(input_file_name_1, usecols=(2,), dtype=int, skiprows=2, unpack=True))
NC_phen_seqs = list(np.loadtxt(input_file_name_1, usecols=(3,), dtype=str, skiprows=2, unpack=True))
NC_sizes = list(np.loadtxt(input_file_name_1, usecols=(4,), dtype=int, skiprows=2, unpack=True))
NC_robs = list(np.loadtxt(input_file_name_1, usecols=(5,), skiprows=2, unpack=True))

# optimal alpha for NC size estimation
input_file_name_2 = input_directory_name_2 + "L.%s_neutral_component_size_estimation_validation_results.txt" % (L)

alpha_list = list(np.loadtxt(input_file_name_2, usecols=(0,), skiprows=2, unpack=True))
opt_check_list = list(np.loadtxt(input_file_name_2, usecols=(4,), dtype=str, skiprows=2, unpack=True))
if opt_check_list[-1]=="True":
    alpha_opt = alpha_list[-1]
else:
    print "error: optimal alpha not found in input"
    quit()



# help functions ########################################

# help function: NC size and robustness estimation from samples
def NC_size_and_rob_estimation_from_samples(NC_phen_seqs,NC_sample_avg_neutral_mut_per_site_storages,NC_sample_SD_neutral_mut_per_site_storages,number_samples,alpha_opt,alphabet,L):

    NC_size_est_storages = []
    NC_rob_est_storages = []
    
    for index,NC_sample_avg_neutral_mut_per_site_storage in enumerate(NC_sample_avg_neutral_mut_per_site_storages):
    
        NC_size_est_storage = []
        NC_rob_est_storage = []
        
        for n_sample in range(number_samples):
        
            NC_size_est = 1.0
            for l in range(L):
                if NC_phen_seqs[index][l]==".":
                    if (1+NC_sample_avg_neutral_mut_per_site_storages[index][n_sample][l]+alpha_opt*NC_sample_SD_neutral_mut_per_site_storages[index][n_sample][l])>4:
                        NC_size_est = NC_size_est*4.0
                    else:
                        NC_size_est = NC_size_est*(1+NC_sample_avg_neutral_mut_per_site_storages[index][n_sample][l]+alpha_opt*NC_sample_SD_neutral_mut_per_site_storages[index][n_sample][l])
                if NC_phen_seqs[index][l]!=".":
                    if (1+NC_sample_avg_neutral_mut_per_site_storages[index][n_sample][l]+alpha_opt*NC_sample_SD_neutral_mut_per_site_storages[index][n_sample][l])>2:
                        NC_size_est = NC_size_est*2.0
                    else:
                        NC_size_est = NC_size_est*(1+NC_sample_avg_neutral_mut_per_site_storages[index][n_sample][l]+alpha_opt*NC_sample_SD_neutral_mut_per_site_storages[index][n_sample][l])
            NC_size_est_storage.append(NC_size_est)
        
            NC_rob_est = float(sum(NC_sample_avg_neutral_mut_per_site_storages[index][n_sample]))/L/(len(alphabet)-1)
            NC_rob_est_storage.append(NC_rob_est)
            
        NC_size_est_storages.append(NC_size_est_storage)
        NC_rob_est_storages.append(NC_rob_est_storage)

    return NC_size_est_storages, NC_rob_est_storages
    
# help function: error quantification NC size estimation from samples
def calc_error_NC_size_estimation_from_samples(NC_ranks,NC_sizes,NC_size_est_storages,number_samples):

    OE_average = 0.0
    OE_RMSD = 0.0
    
    for index,NC_size_est_storage in enumerate(NC_size_est_storages):
        for n_sample,NC_size_est in enumerate(NC_size_est_storage):
            OE_average += np.log(NC_size_est_storages[index][n_sample])/np.log(10)-np.log(NC_sizes[index])/np.log(10)
            OE_RMSD += (np.log(NC_size_est_storages[index][n_sample])/np.log(10)-np.log(NC_sizes[index])/np.log(10))**2
            
    OE_average = float(OE_average)/len(NC_size_est_storages)/number_samples
    OE_RMSD = np.sqrt(float(OE_RMSD)/len(NC_size_est_storages)/number_samples)
    
    OE_SD = 0.0
    
    for index,NC_size_est_storage in enumerate(NC_size_est_storages):
        for n_sample,NC_size_est in enumerate(NC_size_est_storage):
            OE_SD += ((np.log(NC_size_est_storages[index][n_sample])/np.log(10)-np.log(NC_sizes[index])/np.log(10))-OE_average)**2

    OE_SD = np.sqrt(float(OE_SD)/len(NC_size_est_storages)/number_samples)
    
    return OE_average, OE_SD, OE_RMSD
            
# help function: error quantification NC robustness estimation from samples
def calc_error_NC_rob_estimation_from_samples(NC_ranks,NC_robs,NC_rob_est_storages,number_samples):

    OE_average = 0.0
    OE_RMSD = 0.0
    
    for index,NC_rob_est_storage in enumerate(NC_rob_est_storages):
        for n_sample,NC_rob_est in enumerate(NC_rob_est_storage):
            OE_average += NC_rob_est_storages[index][n_sample]-NC_robs[index]
            OE_RMSD += (NC_rob_est_storages[index][n_sample]-NC_robs[index])**2
            
    OE_average = float(OE_average)/len(NC_rob_est_storages)/number_samples
    OE_RMSD = np.sqrt(float(OE_RMSD)/len(NC_rob_est_storages)/number_samples)
    
    OE_SD = 0.0
    
    for index,NC_rob_est_storage in enumerate(NC_rob_est_storages):
        for n_sample,NC_rob_est in enumerate(NC_rob_est_storage):
            OE_SD += ((NC_rob_est_storages[index][n_sample]-NC_robs[index])-OE_average)**2

    OE_SD = np.sqrt(float(OE_SD)/len(NC_rob_est_storages)/number_samples)
    
    return OE_average, OE_SD, OE_RMSD



# program ########################################

# results lists
results_sampling_method_list = []
results_w_random_subsampling_list = []
results_sample_size_list = []
results_random_subsample_size_list = []
results_number_samples_list = []

results_NC_size_est_OE_average_list = []
results_NC_size_est_OE_SD_list = []
results_NC_size_est_OE_RMSD_list = []

results_NC_rob_est_OE_average_list = []
results_NC_rob_est_OE_SD_list = []
results_NC_rob_est_OE_RMSD_list = []

# go through all sampling methods
for sampling_method in sampling_method_list:


    # go through all sample sizes
    for sample_size in sample_size_list:
    
        # create results subdirectories (if not exist)
        try:
            os.makedirs(results_directory_name_1 + "%s/sample_size.%s/" % (sampling_method,sample_size))
        except OSError:
            pass
            
        if sampling_method!="random":
        
            try:
                os.makedirs(results_directory_name_1 + "%s_w_random_subsampling/sample_size.%s/" % (sampling_method,sample_size))
            except OSError:
                pass
        
        
        ## without random subsampling
        NC_sample_avg_neutral_mut_per_site_storages = []
        NC_sample_SD_neutral_mut_per_site_storages = []
        
        # go through all NCs (up to and including maximum NC rank and those larger or equal sample size)
        for index,NC_rank in enumerate(NC_ranks):
        
            if NC_ranks[index]<=max_NC_rank and NC_sizes[index]>=sample_size:
            
                # load sampling results
                input_file_name_3 = input_directory_name_3 + "%s/sample_size.%s/" % (sampling_method,sample_size)
                input_file_name_3 += "L.%s_%s_sampling_sample_size.%s_sample_average_neutral_mutations_per_site_NC_rank.%s.txt" % (L,sampling_method,sample_size,NC_rank)
                sample_avg_neutral_mut_per_site_storage = [[] for n_sample in range(number_samples)]
                for l in range(L):
                    load = list(np.loadtxt(input_file_name_3, usecols=(1+l,), skiprows=2, unpack=True))
                    for n_sample in range(number_samples):
                        sample_avg_neutral_mut_per_site_storage[n_sample].append(load[n_sample])
                        
                NC_sample_avg_neutral_mut_per_site_storages.append(sample_avg_neutral_mut_per_site_storage)
                
                input_file_name_4 = input_directory_name_3 + "%s/sample_size.%s/" % (sampling_method,sample_size)
                input_file_name_4 += "L.%s_%s_sampling_sample_size.%s_sample_SD_neutral_mutations_per_site_NC_rank.%s.txt" % (L,sampling_method,sample_size,NC_rank)
                sample_SD_neutral_mut_per_site_storage = [[] for n_sample in range(number_samples)]
                for l in range(L):
                    load = list(np.loadtxt(input_file_name_4, usecols=(1+l,), skiprows=2, unpack=True))
                    for n_sample in range(number_samples):
                        sample_SD_neutral_mut_per_site_storage[n_sample].append(load[n_sample])
                        
                NC_sample_SD_neutral_mut_per_site_storages.append(sample_SD_neutral_mut_per_site_storage)

        # NC size and robustness estimations
        NC_size_est_storages, NC_rob_est_storages = NC_size_and_rob_estimation_from_samples(NC_phen_seqs,NC_sample_avg_neutral_mut_per_site_storages,NC_sample_SD_neutral_mut_per_site_storages,number_samples,alpha_opt,alphabet,L)

        # write results to .txt file
        for index,NC_size_est_storage in enumerate(NC_size_est_storages):
            table = []
            for n_sample in range(number_samples):
                table.append([n_sample])
                table[n_sample].extend([NC_size_est_storages[index][n_sample]])
                table[n_sample].extend([NC_rob_est_storages[index][n_sample]])
            file_name = results_directory_name_1 + "%s/sample_size.%s/" % (sampling_method,sample_size)
            file_name += "L.%s_%s_sampling_sample_size.%s_sample_NC_estimates_NC_rank.%s.txt" % (L,sampling_method,sample_size,NC_ranks[index])
            f = open(file_name,"w")
            headers=["sample #","NC size est.","NC rob est."]
            f.write(tabulate(table,headers=headers))
            f.close()
        
        # error quantification
        NC_size_est_OE_average, NC_size_est_OE_SD, NC_size_est_OE_RMSD = calc_error_NC_size_estimation_from_samples(NC_ranks,NC_sizes,NC_size_est_storages,number_samples)
        NC_rob_est_OE_average, NC_rob_est_OE_SD, NC_rob_est_OE_RMSD = calc_error_NC_rob_estimation_from_samples(NC_ranks,NC_robs,NC_rob_est_storages,number_samples)
        
        results_sampling_method_list.append(sampling_method)
        results_w_random_subsampling_list.append(False)
        results_sample_size_list.append(sample_size)
        results_random_subsample_size_list.append("x")
        results_number_samples_list.append(number_samples)

        results_NC_size_est_OE_average_list.append(NC_size_est_OE_average)
        results_NC_size_est_OE_SD_list.append(NC_size_est_OE_SD)
        results_NC_size_est_OE_RMSD_list.append(NC_size_est_OE_RMSD)

        results_NC_rob_est_OE_average_list.append(NC_rob_est_OE_average)
        results_NC_rob_est_OE_SD_list.append(NC_rob_est_OE_SD)
        results_NC_rob_est_OE_RMSD_list.append(NC_rob_est_OE_RMSD)
                
                
        ## with random subsampling
        if sampling_method!="random":
        
            for j,random_subsample_size in enumerate(random_subsample_size_list):
            
                if random_subsample_size<=sample_size:
        
                    NC_random_subsample_avg_neutral_mut_per_site_storages = []
                    NC_random_subsample_SD_neutral_mut_per_site_storages = []
                    
                    # go through all NCs (up to and including maximum NC rank and those larger or equal sample size)
                    for index,NC_rank in enumerate(NC_ranks):
                    
                        if NC_ranks[index]<=max_NC_rank and NC_sizes[index]>=sample_size:
                        
                            # load sampling results
                            input_file_name_3 = input_directory_name_3 + "%s_w_random_subsampling/sample_size.%s/" % (sampling_method,sample_size)
                            input_file_name_3 += "L.%s_%s_w_random_subsampling_sample_size.%s_random_subsample_size.%s_sample_average_neutral_mutations_per_site_NC_rank.%s.txt" % (L,sampling_method,sample_size,random_subsample_size,NC_rank)
                            random_subsample_avg_neutral_mut_per_site_storage = [[] for n_sample in range(number_samples)]
                            for l in range(L):
                                load = list(np.loadtxt(input_file_name_3, usecols=(1+l,), skiprows=2, unpack=True))
                                for n_sample in range(number_samples):
                                    random_subsample_avg_neutral_mut_per_site_storage[n_sample].append(load[n_sample])
                               
                            NC_random_subsample_avg_neutral_mut_per_site_storages.append(random_subsample_avg_neutral_mut_per_site_storage)
                            
                            input_file_name_4 = input_directory_name_3 + "%s_w_random_subsampling/sample_size.%s/" % (sampling_method,sample_size)
                            input_file_name_4 += "L.%s_%s_w_random_subsampling_sample_size.%s_random_subsample_size.%s_sample_SD_neutral_mutations_per_site_NC_rank.%s.txt" % (L,sampling_method,sample_size,random_subsample_size,NC_rank)
                            random_subsample_SD_neutral_mut_per_site_storage = [[] for n_sample in range(number_samples)]
                            for l in range(L):
                                load = list(np.loadtxt(input_file_name_4, usecols=(1+l,), skiprows=2, unpack=True))
                                for n_sample in range(number_samples):
                                    random_subsample_SD_neutral_mut_per_site_storage[n_sample].append(load[n_sample])
                               
                            NC_random_subsample_SD_neutral_mut_per_site_storages.append(random_subsample_SD_neutral_mut_per_site_storage)

                    # NC size and robustness estimations
                    NC_size_est_storages, NC_rob_est_storages = NC_size_and_rob_estimation_from_samples(NC_phen_seqs,NC_random_subsample_avg_neutral_mut_per_site_storages,NC_random_subsample_SD_neutral_mut_per_site_storages,number_samples,alpha_opt,alphabet,L)
                    
                    # write results to .txt file
                    for index,NC_size_est_storage in enumerate(NC_size_est_storages):
                        table = []
                        for n_sample in range(number_samples):
                            table.append([n_sample])
                            table[n_sample].extend([NC_size_est_storages[index][n_sample]])
                            table[n_sample].extend([NC_rob_est_storages[index][n_sample]])
                        file_name = results_directory_name_1 + "%s_w_random_subsampling/sample_size.%s/" % (sampling_method,sample_size)
                        file_name += "L.%s_%s_w_random_subsampling_sample_size.%s_random_subsample_size.%s_sample_NC_estimates_NC_rank.%s.txt" % (L,sampling_method,sample_size,random_subsample_size,NC_ranks[index])
                        f = open(file_name,"w")
                        headers=["sample #","NC size est.","NC rob est."]
                        f.write(tabulate(table,headers=headers))
                        f.close()
                    
                    # error quantification
                    NC_size_est_OE_average, NC_size_est_OE_SD, NC_size_est_OE_RMSD = calc_error_NC_size_estimation_from_samples(NC_ranks,NC_sizes,NC_size_est_storages,number_samples)
                    NC_rob_est_OE_average, NC_rob_est_OE_SD, NC_rob_est_OE_RMSD = calc_error_NC_rob_estimation_from_samples(NC_ranks,NC_robs,NC_rob_est_storages,number_samples)
                    
                    results_sampling_method_list.append(sampling_method)
                    results_w_random_subsampling_list.append(True)
                    results_sample_size_list.append(sample_size)
                    results_random_subsample_size_list.append(random_subsample_size)
                    results_number_samples_list.append(number_samples)

                    results_NC_size_est_OE_average_list.append(NC_size_est_OE_average)
                    results_NC_size_est_OE_SD_list.append(NC_size_est_OE_SD)
                    results_NC_size_est_OE_RMSD_list.append(NC_size_est_OE_RMSD)

                    results_NC_rob_est_OE_average_list.append(NC_rob_est_OE_average)
                    results_NC_rob_est_OE_SD_list.append(NC_rob_est_OE_SD)
                    results_NC_rob_est_OE_RMSD_list.append(NC_rob_est_OE_RMSD)
                


# write (overall) results to .txt file ########################################

table = [[results_sampling_method_list[i],
          results_w_random_subsampling_list[i],
          results_sample_size_list[i],
          results_random_subsample_size_list[i],
          results_number_samples_list[i],
          results_NC_size_est_OE_average_list[i],
          results_NC_size_est_OE_SD_list[i],
          results_NC_size_est_OE_RMSD_list[i],
          results_NC_rob_est_OE_average_list[i],
          results_NC_rob_est_OE_SD_list[i],
          results_NC_rob_est_OE_RMSD_list[i]] for i,x in enumerate(results_sampling_method_list)]
file_name = results_directory_name_2 + "L.%s_NC_estimations_from_samples_error_quantification.txt" % (L)
f = open(file_name,"w")
headers=["sampling method",
         "random subsampling?",
         "sample size",
         "random subsample size",
         "# samples",
         "NC size est.: OE average",
         "NC size est.: OE SD",
         "NC size est.: OE RMSD",
         "NC rob est.: OE average",
         "NC rob est.: OE SD",
         "NC rob est.: OE RMSD"]
f.write(tabulate(table,headers=headers))
f.close()



print "done!"

quit()
