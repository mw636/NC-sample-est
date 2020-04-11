import random
import sys
import os
import numpy as np
import json
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
random_seed = 1 # seed for random number generator
random.seed(random_seed)



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

NC_ranks = list(np.loadtxt(input_file_name_1, usecols=(0,), dtype=int, skiprows=2, unpack=True))
NC_indices = list(np.loadtxt(input_file_name_1, usecols=(1,), dtype=int, skiprows=2, unpack=True))
NC_phen_indices = list(np.loadtxt(input_file_name_1, usecols=(2,), dtype=int, skiprows=2, unpack=True))
NC_phen_seqs = list(np.loadtxt(input_file_name_1, usecols=(3,), dtype=str, skiprows=2, unpack=True))
NC_sizes = list(np.loadtxt(input_file_name_1, usecols=(4,), dtype=int, skiprows=2, unpack=True))

# dictionary: genotype index -> NC index
input_file_name_2 = input_directory_name + "L.%s_dictionary_genotype_index_to_NC_index.json" % (L)

with open(input_file_name_2, 'r') as handle:
    gNC_map = json.load(handle)

# set up NCs = lists of genotype indices (up to and including maximum NC rank)
NCs = [[] for NC_index in NC_indices]
for gen_index,NC_index in gNC_map.items():
    i = NC_indices.index(NC_index)
    if NC_ranks[i]<=max_NC_rank:
        NCs[i].append(int(gen_index))



# help functions ########################################

# help function (base transfer): genotype sequence -> genotype index
def gen_seq_to_gen_index(gen_seq,alphabet,L):
    gen_index = 0
    for site_index in range(L):
       gen_index += alphabet.index(gen_seq[site_index])*len(alphabet)**(L-site_index-1)
    return gen_index

# help function (base transfer): genotype index -> genotype sequence
def gen_index_to_gen_seq(gen_index,alphabet,L):
    gen_seq = []
    for site_index in range(L):
        value = gen_index // len(alphabet)**(L-site_index-1)
        letter = alphabet[value]
        gen_seq.append(letter)
        gen_index -= value*len(alphabet)**(L-site_index-1)
    return ''.join(gen_seq)

# help function: 'RW' neutral mutation for a genotype index
def RW_neutral_mutation(gen_index,ref_NC_index,gNC_map,alphabet,L):
    gen_seq = gen_index_to_gen_seq(gen_index,alphabet,L)
    gen_seq_array = list(gen_seq)
    success = False
    while success==False:
        # random position
        mut_position = random.randint(0,L-1)
        # mutation alphabet
        mut_alphabet = [letter for letter in alphabet if letter!=gen_seq_array[mut_position]]
        # random letter
        mut_letter = random.choice(mut_alphabet)
        mut_gen_seq_array = [j for j in gen_seq_array] # to preserve "gen_seq_array"
        # mutate
        mut_gen_seq_array[mut_position] = mut_letter
        mut_gen_index = gen_seq_to_gen_index(''.join(mut_gen_seq_array),alphabet,L)
        # check if neutral
        try:
            mut_NC_index = gNC_map["%s" % mut_gen_index]
            if mut_NC_index==ref_NC_index:
                success=True
        except KeyError:
            pass
    return mut_gen_index

# help function: 'site scanning' neutral mutation for a genotype index
def site_scanning_neutral_mutation(gen_index,ref_mut_position,ref_NC_index,gNC_map,alphabet,L):
    gen_seq = gen_index_to_gen_seq(gen_index,alphabet,L)
    gen_seq_array = list(gen_seq)
    success = False
    while success==False:
        # position
        mut_position = ref_mut_position
        # mutation alphabet
        mut_alphabet = [letter for letter in alphabet if letter!=gen_seq_array[mut_position]]
        # randomly test mutation alphabet until success
        while len(mut_alphabet)>0:
            # random letter
            mut_letter = random.choice(mut_alphabet)
            mut_gen_seq_array = [j for j in gen_seq_array] # to preserve "gen_seq_array"
            # mutate
            mut_gen_seq_array[mut_position] = mut_letter
            mut_gen_index = gen_seq_to_gen_index(''.join(mut_gen_seq_array),alphabet,L)
            # check if neutral
            try:
                mut_NC_index = gNC_map["%s" % mut_gen_index]
                if mut_NC_index==ref_NC_index:
                    success=True
                    break
            except KeyError:
                pass
            # if no success, update mutation alphabet
            mut_alphabet.remove(mut_letter)
        # if no success, go to next position
        if len(mut_alphabet)==0 and success==False:
            ref_mut_position = (ref_mut_position+1) % L
    return mut_gen_index, mut_position

# help function: find neutral mutations per site for a genotype index
def find_neutral_mut_per_site(gen_index,ref_NC_index,gNC_map,alphabet,L):
    gen_seq = gen_index_to_gen_seq(gen_index,alphabet,L)
    gen_seq_array = list(gen_seq)
    gen_neutral_mut_per_site = [0]*L
    # go through all characters of the genotype sequence
    for i,character in enumerate(gen_seq_array):
        # go through all letters of the alphabet
        for letter in alphabet:
            # mutate character to "value" of letter if character != letter
            if character!=letter:
                mut_gen_seq_array = [j for j in gen_seq_array] # to preserve "gen_seq_array"
                # mutate
                mut_gen_seq_array[i] = letter
                mut_gen_index = gen_seq_to_gen_index(''.join(mut_gen_seq_array),alphabet,L)
                # check if neutral
                try:
                    mut_NC_index = gNC_map["%s" % mut_gen_index]
                    if mut_NC_index==ref_NC_index:
                        # if neutral mutation, add it to the specific site
                        gen_neutral_mut_per_site[i] += 1
                except KeyError:
                    pass
    return gen_neutral_mut_per_site

# help function: sample measurement (averaging)
def sample_measurement(sample_gen_neutral_mut_per_site_storage,L):

    sample_avg_neutral_mut_per_site = [0 for l in range(L)]
    sample_SD_neutral_mut_per_site = [0 for l in range(L)]

    for k in range(len(sample_gen_neutral_mut_per_site_storage)):
        for l in range(L):
            sample_avg_neutral_mut_per_site[l] += sample_gen_neutral_mut_per_site_storage[k][l]
    for l in range(L):
        sample_avg_neutral_mut_per_site[l] = float(sample_avg_neutral_mut_per_site[l])/len(sample_gen_neutral_mut_per_site_storage)

    for k in range(len(sample_gen_neutral_mut_per_site_storage)):
        for l in range(L):
            sample_SD_neutral_mut_per_site[l] += (sample_gen_neutral_mut_per_site_storage[k][l]-sample_avg_neutral_mut_per_site[l])**2
    for l in range(L):
        if len(sample_gen_neutral_mut_per_site_storage)>1:
            sample_SD_neutral_mut_per_site[l] = np.sqrt(float(sample_SD_neutral_mut_per_site[l])/(len(sample_gen_neutral_mut_per_site_storage)-1))
        if len(sample_gen_neutral_mut_per_site_storage)==1:
            sample_SD_neutral_mut_per_site[l] = 0

    return sample_avg_neutral_mut_per_site, sample_SD_neutral_mut_per_site



# program ########################################

# go through all sampling methods
for sampling_method in sampling_method_list:


    # go through all sample sizes
    for sample_size in sample_size_list:
    
        # create results subdirectories (if not exist)
        try:
            os.makedirs(results_directory_name + "%s/sample_size.%s/" % (sampling_method,sample_size))
        except OSError:
            pass
            
        if sampling_method!="random":
        
            try:
                os.makedirs(results_directory_name + "%s_w_random_subsampling/sample_size.%s/" % (sampling_method,sample_size))
            except OSError:
                pass
    
    
        # go through all NCs (up to and including maximum NC rank and those larger or equal sample size)
        for index,NC_rank in enumerate(NC_ranks):
        
            ref_NC_index = NC_indices[index]
            
            # NC = list of NC genotype indices
            NC = NCs[index]

            if NC_ranks[index]<=max_NC_rank and NC_sizes[index]>=sample_size:

                # storages of measurement results for all samples considered for this sampling method, sample size, and NC, respectively
                sample_avg_neutral_mut_per_site_storage = []
                sample_SD_neutral_mut_per_site_storage = []
                
                if sampling_method!="random":
                
                    random_subsample_avg_neutral_mut_per_site_storage = [[] for random_subsample_size in random_subsample_size_list]
                    random_subsample_SD_neutral_mut_per_site_storage = [[] for random_subsample_size in random_subsample_size_list]


                # go through all samples
                for n_sample in range(number_samples):


                    # sample creation + one-point mutational neighbourhood measurements of sample genotypes
                    sample = []
                    sample_gen_neutral_mut_per_site_storage = []


                    if sampling_method=="random":
                    
                        sample = random.sample(NC,sample_size)

                        for gen_index in sample:
                            gen_neutral_mut_per_site = find_neutral_mut_per_site(gen_index,ref_NC_index,gNC_map,alphabet,L)
                            sample_gen_neutral_mut_per_site_storage.append(gen_neutral_mut_per_site)


                    if sampling_method=="RW":
                    
                        # randomly choose first genotype index
                        ref_gen_index = random.sample(NC,1)[0]
                        sample = [ref_gen_index]
                    
                        # fill sample
                        while len(sample) < sample_size:
                            mut_gen_index = RW_neutral_mutation(ref_gen_index,ref_NC_index,gNC_map,alphabet,L)
                            sample.append(mut_gen_index)
                            ref_gen_index = mut_gen_index
                        
                        for gen_index in sample:
                            gen_neutral_mut_per_site = find_neutral_mut_per_site(gen_index,ref_NC_index,gNC_map,alphabet,L)
                            sample_gen_neutral_mut_per_site_storage.append(gen_neutral_mut_per_site)


                    if sampling_method=="site_scanning":
                    
                        # randomly choose first genotype index
                        ref_gen_index = random.sample(NC,1)[0]
                        sample = [ref_gen_index]
                        ref_mut_position = 0

                        # fill sample
                        while len(sample) < sample_size:
                            mut_gen_index, mut_position = site_scanning_neutral_mutation(ref_gen_index,ref_mut_position,ref_NC_index,gNC_map,alphabet,L)
                            sample.append(mut_gen_index)
                            ref_gen_index = mut_gen_index
                            ref_mut_position = (mut_position+1) % L
                        
                        for gen_index in sample:
                            gen_neutral_mut_per_site = find_neutral_mut_per_site(gen_index,ref_NC_index,gNC_map,alphabet,L)
                            sample_gen_neutral_mut_per_site_storage.append(gen_neutral_mut_per_site)
                            

                    # sample measurement (averaging)
                    sample_avg_neutral_mut_per_site, sample_SD_neutral_mut_per_site = sample_measurement(sample_gen_neutral_mut_per_site_storage,L)

                    sample_avg_neutral_mut_per_site_storage.append(sample_avg_neutral_mut_per_site)
                    sample_SD_neutral_mut_per_site_storage.append(sample_SD_neutral_mut_per_site)


                    # random subsampling
                    if sampling_method!="random":

                        for j,random_subsample_size in enumerate(random_subsample_size_list):
                        
                            if random_subsample_size<=sample_size:
                            
                                # subsample creation (indices indicating "position" within original sample)
                                random_indices = random.sample(range(sample_size),random_subsample_size)
                                
                                random_subsample_gen_neutral_mut_per_site_storage = []
                                
                                for random_index in random_indices:
                                    random_subsample_gen_neutral_mut_per_site_storage.append(sample_gen_neutral_mut_per_site_storage[random_index])
                                
                                # random subsample measurement (averaging)
                                random_subsample_avg_neutral_mut_per_site, random_subsample_SD_neutral_mut_per_site = sample_measurement(random_subsample_gen_neutral_mut_per_site_storage,L)

                                random_subsample_avg_neutral_mut_per_site_storage[j].append(random_subsample_avg_neutral_mut_per_site)
                                random_subsample_SD_neutral_mut_per_site_storage[j].append(random_subsample_SD_neutral_mut_per_site)


                # write results to .txt file
                
                table = []
                for n_sample in range(number_samples):
                    table.append([n_sample])
                    table[n_sample].extend(sample_avg_neutral_mut_per_site_storage[n_sample])
                file_name = results_directory_name + "%s/sample_size.%s/" % (sampling_method,sample_size)
                file_name += "L.%s_%s_sampling_sample_size.%s_sample_average_neutral_mutations_per_site_NC_rank.%s.txt" % (L,sampling_method,sample_size,NC_rank)
                f = open(file_name,"w")
                headers=["sample #"] + [""]*L
                f.write(tabulate(table,headers=headers))
                f.close()
                
                table = []
                for n_sample in range(number_samples):
                    table.append([n_sample])
                    table[n_sample].extend(sample_SD_neutral_mut_per_site_storage[n_sample])
                file_name = results_directory_name + "%s/sample_size.%s/" % (sampling_method,sample_size)
                file_name += "L.%s_%s_sampling_sample_size.%s_sample_SD_neutral_mutations_per_site_NC_rank.%s.txt" % (L,sampling_method,sample_size,NC_rank)
                f = open(file_name,"w")
                headers=["sample #"] + [""]*L
                f.write(tabulate(table,headers=headers))
                f.close()

                if sampling_method!="random":

                        for j,random_subsample_size in enumerate(random_subsample_size_list):
                        
                            if random_subsample_size<=sample_size:

                                table = []
                                for n_sample in range(number_samples):
                                    table.append([n_sample])
                                    table[n_sample].extend(random_subsample_avg_neutral_mut_per_site_storage[j][n_sample])
                                file_name = results_directory_name + "%s_w_random_subsampling/sample_size.%s/" % (sampling_method,sample_size)
                                file_name += "L.%s_%s_w_random_subsampling_sample_size.%s_random_subsample_size.%s_sample_average_neutral_mutations_per_site_NC_rank.%s.txt" % (L,sampling_method,sample_size,random_subsample_size,NC_rank)
                                f = open(file_name,"w")
                                headers=["sample #"] + [""]*L
                                f.write(tabulate(table,headers=headers))
                                f.close()

                                table = []
                                for n_sample in range(number_samples):
                                    table.append([n_sample])
                                    table[n_sample].extend(random_subsample_SD_neutral_mut_per_site_storage[j][n_sample])
                                file_name = results_directory_name + "%s_w_random_subsampling/sample_size.%s/" % (sampling_method,sample_size)
                                file_name += "L.%s_%s_w_random_subsampling_sample_size.%s_random_subsample_size.%s_sample_SD_neutral_mutations_per_site_NC_rank.%s.txt" % (L,sampling_method,sample_size,random_subsample_size,NC_rank)
                                f = open(file_name,"w")
                                headers=["sample #"] + [""]*L
                                f.write(tabulate(table,headers=headers))
                                f.close()



print "done!"

quit()
